"""
======================COPYRIGHT/LICENSE START==========================

PippFormat.py: Contains functions specific to Pipp conversions.

Copyright (C) 2005-2007 Wim Vranken (European Bioinformatics Institute)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)
- PDBe website (http://www.ebi.ac.uk/pdbe/)

- contact Wim Vranken (wim@ebi.ac.uk)
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================
"""

import copy

from memops.api import Implementation
from memops.universal.Util import returnInt

from ccp.api.nmr import Nmr
from ccp.general.Util import setUniqueAppData

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords
from ccpnmr.format.general.Constants import height_kw
from ccpnmr.format.general.Constants import defaultMolCode
from ccpnmr.format.general.Util import getResName

#
# Add some information to IOkeywords...
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readShifts']['getSeqFlag'] = (False,False,"Read sequence from the chemical shift file.")

class PippFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'pipp'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    
    self.getPeaks = self.getPeaksGeneric
    
    self.getMeasurements = self.getMeasurementsGeneric
    self.createMeasurementFile = self.createMeasurementFileGeneric

  #
  # Deviations from generic import stuff
  #

  def getMeasurementsSequenceRead(self):
    
    sequenceText = ''
    
    # If shifts, check whether sequence has to be created as well...
    if not self.project.molSystems and self.measurementType == 'shift' and self.getSeqFlag:
      self.setSequenceFileClass()
      self.sequenceFile = self.SequenceFileClass(self.fileName)
      self.sequenceFile.readFromShifts(self.measurementFile)
      self.file = 'temporary'
      self.readSequence(minimalPrompts = self.minimalPrompts)
      self.file = None
      
      sequenceText = 'and sequence '
    
    return sequenceText

  #
  # Functions different to default functions in DataFormat
  #

  def setPeakIntensity(self):

    # PeakIntensity attributes
    if self.rawPeak.intensity != 0:
      peakInt = Nmr.PeakIntensity(self.peak,value = self.rawPeak.intensity, method = self.methods[self.format]['Intensity'])
    
      # Assuming always intensity for pipp
      peakInt.intensityType = height_kw

  def setPeakExtras(self):
  
    if self.rawPeak.assign[0]:
      
      #
      # Make peakContrib...
      #
      
      self.peakContribs.append(Nmr.PeakContrib(self.peak))

  def setPeakDim(self):

    # TODO: rawPeak.dimSlice[]

    dataDimRef = self.dataDimRefs[self.rawPeakDimIndex]

    self.peakDim = self.peak.findFirstPeakDim(dim = dataDimRef.dataDim.dim)

    self.peakDim.dataDimRef = dataDimRef

    self.peakDim.value = self.rawPeak.ppm[self.rawPeakDimIndex]

  def getPeakResNames(self):
  
    self.resNames = []

    if self.rawPeak.assign[self.rawPeakDimIndex] not in [None]:
    
      try:
        self.getSinglePeakResName(self.rawPeak.assign[self.rawPeakDimIndex])
        
      except:
        print "  Problems handling %s assignment name %s... ignored." % (self.format,self.rawPeak.assign[self.rawPeakDimIndex])

  def getSinglePeakResName(self,assignName):
  
    (residue,atomName) = assignName.split(self.peakFile.assignSep)
    
    resLabel = residue[0]
    seqCode = returnInt(residue[1:])
    
    if atomName.count('|'):
      (atomName,otherAtomName) = atomName.split('|')
    
    self.resNames.append(getResName(defaultMolCode,seqCode,atomName))

  def createShift(self,resonance,chemShift):
    
    #
    # Also set ambiguity code if listed
    #
    
    if chemShift.ambCode != None:
      setUniqueAppData('AppDataInt',resonance,self.format,'ambiguityCode',chemShift.ambCode)
      
    #
    # Create shift linked to resonance
    #

    self.convertCount[self.mainCode][1][self.valueCode] += 1

    return Nmr.Shift(self.measurementList, value = chemShift.value, resonance = resonance)

  def setChemShiftFileValue(self):
  
    self.setMeasurementFileGenericValue()

  def setMeasurementFileGenericValue(self):
      
    if self.atomNamesDict.has_key(self.atomName):
    
      resonanceToAtom = self.atomNamesDict[self.atomName]
      measurement = self.origAtomMeasurements[resonanceToAtom]
      residueLabel = self.getResidueLabelThreeLetter(resonanceToAtom)

      otherAtomNames = []
      
      # Get a list of the other prochiral atoms if not stereospecifically assigned...
      if len(resonanceToAtom.resonance.resonanceSet.atomSets) > 1:

        if resonanceToAtom.chemAtom:
          chemAtomOrSet = resonanceToAtom.chemAtom
        else:
          chemAtomOrSet = resonanceToAtom.chemAtomSet
                
        # Also catch more unusual situations (single prochiral, side chain N GLN, ...)
        # TODO THIS IS NOT WORKING PROPERLY!
        if not chemAtomOrSet:
          chemAtomOrSet = resonanceToAtom.resonance.resonanceSet.findFirstAtomSet().findFirstAtom().chemAtom.chemAtomSet
        
        if chemAtomOrSet and chemAtomOrSet.chemAtomSet:
          
          chemComp = chemAtomOrSet.chemComp
          namingSystem = chemComp.findFirstNamingSystem(name = self.namingSystemName)
          
          chemAtomSet = chemAtomOrSet.chemAtomSet
          
          if resonanceToAtom.chemAtom:
            chemAtomOrSets = chemAtomSet.chemAtoms
          else:
            chemAtomOrSets = chemAtomSet.chemAtomSets
          
          for chemAtomOrSet in chemAtomOrSets:
            chemAtomSysName = self.getChemAtomOrSetSysName((chemAtomOrSet,namingSystem))
          
            if not chemAtomSysName:
              iupacNamingSystem = chemComp.findFirstNamingSystem(name = 'IUPAC')
              chemAtomSysName =  self.getChemAtomOrSetSysName((chemAtomOrSet,iupacNamingSystem))
            
            if chemAtomSysName and chemAtomSysName.sysName != self.atomName:
              otherAtomNames.append(chemAtomSysName.sysName)
      
      rawChemShift = self.rawMeasurementClass(None)
      rawChemShift.setValue('seqCode',self.seqCode)
      rawChemShift.setValue('resLabel',residueLabel)
      rawChemShift.setValue('atomName',self.atomName)
      rawChemShift.setValue('otherAtomNames',tuple(otherAtomNames))
      rawChemShift.setValue('value',measurement.value)
                                              
      self.measurementFileValues.append(rawChemShift)
      
      if self.atomMeasurements.has_key(resonanceToAtom):
        del(self.atomMeasurements[resonanceToAtom])

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)

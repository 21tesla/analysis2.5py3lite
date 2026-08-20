"""
======================COPYRIGHT/LICENSE START==========================

TalosFormat.py: Contains functions specific to Talos conversions.

Copyright (C) 2005-2008 Wim Vranken (European Bioinformatics Institute)

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


from ccpnmr.format.general.Util import setOriginalData

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccpnmr.format.general.Constants import volume_kw, height_kw

import ccp.api.nmr.Nmr as Nmr

#
# Additional IOkeywords definitions
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readDihedralConstraints']['validCount'] =        (None,False,'The constraint is only considered valid if\n its TALOS count is higher than this number.')
IOkeywords['readDihedralConstraints']['validClasses'] =      (None,False,'The constraint is only considered valid if\n its TALOS class is one of the ones selected.')
IOkeywords['readDihedralConstraints']['multiplyDeviation'] = (2.0,False,'The standard angle deviation reported by TALOS will be multiplied\n by this factor to get the dihedral angle constraint limits.')
IOkeywords['writeProject']['chain'] =               (None,True,'The chain to be exported to the TALOS input file.')
IOkeywords['writeProject']['shiftList'] =           (None,True,'The shift list to be exported to the TALOS input file.')
IOkeywords['writeProject']['useOriginalData'] =     (False,False,'Use original data for export (if available).')

class TalosFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'talos'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    self.createSequenceFile = self.createSequenceFileGeneric

    self.getMeasurements = self.getMeasurementsGeneric
    self.createMeasurementFile = self.createMeasurementFileGeneric

    self.getConstraints = self.getConstraintsGeneric

  #
  # Deviations from generic import stuff
  #

  def getFullProject(self,fileName,chemShiftKeyWdList = None, sequenceKeyWds = None):
  
    if self.verbose == 1:
      print "Reading %s project from file %s" % (self.formatLabel,fileName)
   
    self.file = self.projectIO.TalosProjectFile(fileName)
    self.file.read()
    
    #
    # Process elements - call functions 'as if' reading real file, but pass
    # on objects created during project import
    #
    
    self.fileName = fileName
    
    if self.file.sequenceFile:
      self.sequenceFile = self.file.sequenceFile
      
      if sequenceKeyWds:
        keywds = sequenceKeyWds
      else:
        keywds = {}
        
      self.readSequence(minimalPrompts = self.minimalPrompts, **keywds)
    
    for i in range(0,len(self.file.chemShiftFiles)):
      
      self.measurementFile = self.file.chemShiftFiles[i]

      if chemShiftKeyWdList:
        keywds = chemShiftKeyWdList[i]
      else:
        keywds = {}

      self.readShifts(minimalPrompts = self.minimalPrompts, **keywds)
    
    #
    # Avoid reading project again
    #
    
    self.file = None

    return self.project


  def createFullProject(self,fileName):
  
    if self.verbose == 1:
      print "Writing %s project from file %s" % (self.formatLabel,fileName)
   
    self.file = self.projectIO.TalosProjectFile(fileName)
   
    self.writeSequence(fileName, chains = [self.chain], noWrite = True, resetIOkeywords = False)
    self.file.sequenceFile = self.sequenceFile

    self.writeShifts(fileName, measurementList = self.shiftList, chains = [self.chain], noWrite = True,  resetMapping = False, resetIOkeywords = False)
    self.file.chemShiftFiles = [self.measurementFile]
    
    self.file.write()
    
    return self.file

  #
  # Functions different to default functions in DataFormat
  #

  def createSequence(self):
  
    self.sequenceFile.sequences.append(self.sequenceIO.TalosSequence(molName = self.chain.molecule.name))
    self.sequence = self.sequenceFile.sequences[-1]

  def setSequenceFileElements(self):
    
    self.sequence.elements.append(self.sequenceIO.TalosSequenceElement(self.seqCode,self.getResidueCode1Letter(self.residue)))

  def rawConstraintValid(self):
  
    (validCount,validClass) = (True,True)
    
    if self.validCount != None:
      validCount = False
      if self.rawConstraint.count >= self.validCount:
        validCount = True
    
    if self.validClasses != None:
      validClass = False
      if self.rawConstraint.className in self.validClasses:
        validClass = True
    
    if validCount and validClass:
      constraintValid = True
    else:
      constraintValid = False
  
    return constraintValid

  def setDihedralConstraintInfo(self):
    
    #
    # Keep track of the variables that are not used in the data model
    #
    
    for origKeywd in ['count','className','dist']:
    
      setOriginalData(self.format,self.constraint,self.rawConstraint,origKeywd)
      
  def getDihedralAngleErrorRange(self):
    
    return self.rawConstraint.devAngle * self.multiplyDeviation

  def setChemShiftFileValue(self):
  
    if self.atomNamesDict.has_key(self.atomName):
    
      resonanceToAtom = self.atomNamesDict[self.atomName]
      chemShift = self.origAtomMeasurements[resonanceToAtom]

      residueLabel = self.getResidueLabelOneLetter(resonanceToAtom)

      rawMeasurement = self.rawMeasurementClass(self.measurementFile)

      rawMeasurement.seqCode = self.seqCode
      rawMeasurement.resLabel = residueLabel
      rawMeasurement.atomName = self.atomName
      rawMeasurement.value = chemShift.value

      self.measurementFileValues.append(rawMeasurement)

      if self.atomMeasurements.has_key(resonanceToAtom):
        del(self.atomMeasurements[resonanceToAtom])

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)

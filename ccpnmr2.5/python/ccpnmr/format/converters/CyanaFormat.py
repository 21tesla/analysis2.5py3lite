"""
======================COPYRIGHT/LICENSE START==========================

CyanaFormat.py: Contains functions specific to Cyana conversions.

Copyright (C) 2005 Wim Vranken (European Bioinformatics Institute)

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

import traceback, sys, copy

from ccpnmr.format.converters.DyanaFormat import DyanaFormat, IOkeywords

from ccpnmr.format.general.Constants import defaultMolCode
from ccpnmr.format.general.Util import getResName

#
# Add some information to IOkeywords...
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readPeaks']['cyanaFormat'] = (True,False,'Read as CYANA format (with multiple assignments).')
IOkeywords['readPeaks']['cyanaTextAssignments'] = (True,False,'Read assignments as text, not serials.')
IOkeywords['readPeaks']['forceHeight'] = (None,False,'Force import of XEASY volumes as heights in CCPN.')
IOkeywords['writePeaks']['cyanaFormat'] = (True,False,'Write as CYANA format (with multiple assignments).')
IOkeywords['writePeaks']['integrationMethod'] = ('volume',False,'Type of integration method to use in output file.')
IOkeywords['writePeaks']['cyanaTextAssignments'] = (True,False,'Write assignments as text, not serials.')

class CyanaFormat(DyanaFormat):

  def setFormat(self):
  
    self.format = 'cyana'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getCoordinates = self.getCoordinatesGeneric
    self.createCoordinateFile = self.createCoordinateFileGeneric 

    self.createPeakFile = self.createPeakFileGeneric

  #
  # Deviations from generic import stuff
  #

  def getCoordinatesSetSequenceFile(self):
    
    self.setSequenceFileClass()
    self.sequenceFile = self.SequenceFileClass(self.fileName, version =  self.version)
    self.sequenceFile.readFromCoordinates(self.coordinateFile)

  def getSequence(self):
    
    # Not generic - passing in version    
    self.sequenceFile = self.SequenceFileClass(self.fileName, version =  self.version)
    self.sequenceFile.read()

    if self.verbose == 1:
      print "Reading sequence from %s file %s" % (self.formatLabel,self.fileName)

  def createSequenceFile(self):
  
    # Not generic - passing in version
    if self.verbose == 1:
      print "Writing sequence to %s file %s" % (self.formatLabel,self.fileName)

    self.sequenceFile = self.SequenceFileClass(self.fileName, version =  self.version)

  def createPeakFileFormatSpecific(self):

    self.peakFile.setSpectrumInfo(self.fileName,self.numPeakDim)
    self.writeKeywds['cyanaFormat'] = self.cyanaFormat
    self.writeKeywds['cyanaTextAssignments'] = self.cyanaTextAssignments

  def getPeaks(self):
  
    if self.verbose == 1:
      print "Reading peak list from %s file %s" % (self.formatLabel,self.fileName)

    # TODO HERE: have to figure out what to do if project file read...
   
    # Read in whole file (even if only experiment wanted later on...)
    
    self.peakFile = self.PeakFileClass(self.fileName)
    self.peakFile.read(cyanaFormat=self.cyanaFormat,cyanaTextAssignments=self.cyanaTextAssignments)

  def getConstraints(self):
    
    #
    # Format slightly different for constraints....
    #
    
    try:
      
      addKeywords =self.getLowerLimitFile()
    
      self.constraintFile = self.ConstraintFileClass(self.fileName, version =  self.version, **addKeywords)
      self.constraintFile.read(self.fileName)

      if self.verbose == 1:
        print "Reading %s constraint list from %s file %s" % (self.constraintType,self.formatLabel,self.fileName)

    except:
      errorMessage = traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
      self.messageReporter.showWarning("Warning"," Cannot read %s constraints for %s...:\n%s" % (self.constraintApiCode,self.formatLabel,errorMessage),self.guiParent)
      self.constraintFile = None

      return traceback.format_exception(sys.exc_type,sys.exc_value,sys.exc_info()[2]) 

  def createConstraintFile(self):
  
    try:

      addKeywords = self.getLowerLimitFile(messageType = 'write')

      self.constraintFile = self.ConstraintFileClass(self.fileName, version =  self.version, **addKeywords)
    
      if self.verbose == 1:
        print "Writing %s constraints to %s file %s" % (self.constraintType,self.formatLabel,self.fileName)
        
      self.setSpecificConstraintFileWriteInfo()

    except:
    
      self.messageReporter.showWarning("Warning"," No write%sConstraints available for %s..." % (self.constraintApiCode,self.formatLabel),self.guiParent)
      self.quit()

  #
  # Functions different to default functions in DataFormat
  #
    
  def createSequence(self):
  
    self.sequenceFile.sequences.append(self.sequenceIO.CyanaSequence(molName = self.chain.molecule.name))
    self.sequence = self.sequenceFile.sequences[-1]
    
  # setSequenceFileElements is same as Dyana - bit hacky...

  def setPeakDim(self):
  
    dataDimRef = self.dataDimRefs[self.rawPeakDimIndex]

    self.peakDim = self.peak.findFirstPeakDim(dim = dataDimRef.dataDim.dim)

    self.peakDim.dataDimRef = dataDimRef
    
    self.peakDim.value = self.rawPeak.ppm[self.rawPeakDimIndex]

  def getPeakResNames(self):
  
    self.resNames = []

    if self.rawPeak.assign[self.rawPeakDimIndex] not in ['-']:
      
      for assignment in [self.rawPeak.assign] + self.rawPeak.ambiguousAssign:
      
        (atomName,seqCode) = assignment[self.rawPeakDimIndex].split('.')
        
        self.resNames.append(getResName(defaultMolCode,seqCode,atomName))
        
  def setPeakFilePeakExtras(self):
  
    # Taken directly from XEasyFormat
    
    self.colour = 1              # TODO: application data?
    self.userCode = "T"          # TODO: application data?
    self.assign = []             # Can only do this if atomSerial assigned...
    self.ambiguousAssign = []    # Can only do this if atomSerial assigned and cyanaFormat option selected.
    self.ppm = []
    
  def setPeakFilePeakIntensity(self):
  
    self.volume = 0.0
    self.volumeError = 0.0
    self.intCode = '-'

    intensity = self.peak.findFirstPeakIntensity(intensityType = self.integrationMethod)
    
    if intensity:
      
      self.volume = intensity.value
      
      if intensity.error != None:
        self.volumeError = intensity.error
      else:
        self.volumeError = 0.0
        
      intMethod = intensity.method.procedure
      
      #
      # XEasy specific codes...
      # Now assuming 'automatic' as default if a volume is available
      #
      
      self.intCode = self.peakFile.translateIntMethod('ccpn',intMethod)

  def setPeakFilePeakDimInfo(self):
    
    self.ppm.append(self.peakDim.getValue())

    peakDimIndex = list(self.peak.sortedPeakDims()).index(self.peakDim)
    
    assignLen = len(self.peakAssignmentList)
    resName = '-'
    
    if assignLen:

      resToAtom = self.peakAssignmentList[0][peakDimIndex]
      resName = self.getCyanaResName(resToAtom)

      for i in range(1,assignLen):

        resToAtom = self.peakAssignmentList[i][peakDimIndex]

        if len(self.ambiguousAssign) < assignLen - 1:
          self.ambiguousAssign.append([])

        self.ambiguousAssign[i-1].append(self.getCyanaResName(resToAtom))
        
    self.assign.append(resName)      

  def getCyanaResName(self,resToAtom):
    
    resName = '-'
    
    if resToAtom:
      chain = resToAtom.chain
      seqId = resToAtom.seqId

      seqCode = self.getExportSeqCode(self.chainDict[chain][1],chain.findFirstResidue(seqId = seqId))

      resName = "%s.%s" % (resToAtom.atomName,str(seqCode))
    
    return resName

  def createPeakFilePeak(self):
    
    self.peakFile.peaks.append(self.peaksIO.CyanaPeak(self.peakNum,self.ppm,self.colour,self.userCode,self.volume,self.volumeError,self.intCode,self.assign,textAssign=True)) 
    self.peakFile.peaks[-1].ambiguousAssign = self.ambiguousAssign[:]

  def setPeakFileInfo(self):    
    
    #
    # Get isotopecode info
    #

    for dimIndex in range(0,len(self.dataDimRefs)):
      expDimRef = self.dataDimRefs[dimIndex].expDimRef
      
      elementType = expDimRef.isotopeCodes[0][-1]
      
      expTransfers = list(expDimRef.expTransfers)
 
      if len(expTransfers) == 1 and expTransfers[0].transferType == 'NOESY':
        elementType = elementType.lower()
      
      self.peakFile.dimCodes[dimIndex] = elementType
      
    #
    # Fix in homonuclear 2D case...
    #
    
    if self.peakFile.dimCodes.count(self.peakFile.dimCodes[0]) == len(self.peakFile.dimCodes):
      if self.peakFile.dimCodes[0] == self.peakFile.dimCodes[0].lower():
        self.peakFile.dimCodes[0] = self.peakFile.dimCodes[0].upper()

    #
    # Also do a check for which intensity type is available, reset if not good
    #
    
    intensityTypes = []
    
    for peak in self.peakList.peaks: 
      for intensity in peak.findAllPeakIntensities():
        if not intensity.intensityType in intensityTypes:
          intensityTypes.append(intensity.intensityType)
    
    if self.integrationMethod not in intensityTypes and len(intensityTypes) == 1:
      print "  Warning: resetting peak integration method to %s - no %s values available." % (intensityTypes[0],self.integrationMethod)
      self.integrationMethod = intensityTypes[0]
    

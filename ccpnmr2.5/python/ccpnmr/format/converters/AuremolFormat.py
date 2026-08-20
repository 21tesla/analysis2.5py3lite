"""
======================COPYRIGHT/LICENSE START==========================

AuremolFormat.py: Contains functions specific to Auremol conversions.

Copyright (C) 2007 Wim Vranken (European Bioinformatics Institute)

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

import string, copy, os, traceback, sys

from memops.universal.Io import joinPath, splitPath

from memops.api import Implementation
from ccp.api.nmr import Nmr

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords
from ccpnmr.format.general.Constants import ccpNmr_kw, volume_kw, height_kw

from ccp.general.Util import setUniqueAppData

from ccp.format.general.Constants import defaultSeqInsertCode

class AuremolFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'auremol'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    
    self.getPeaks = self.getPeaksGeneric

    self.getCoordinates = self.getCoordinatesGeneric
    self.createCoordinateFile = self.createCoordinateFileGeneric 

  #
  # Deviations from generic import stuff
  #

  def getCoordinatesSetSequenceFile(self):
    
    self.setSequenceFileClass()
    self.sequenceFile = self.SequenceFileClass(self.fileName)
    self.sequenceFile.readFromCoordinates(self.coordinateFile)

  def getMeasurements(self):
  
    #
    # Get a measurement list
    #
    
    try:

      if self.verbose == 1:
        print "Reading %ss from %s file %s" % (self.measurementType,self.formatLabel,self.fileName)

      # TODO HERE: have to figure out what to do if project file read...
  
      # Read initial metafile, make sure all comp files exist!
      
      metaFile = self.generalIO.AuremolMetaFile(self.fileName)
      metaFile.read()
      
      for compound in metaFile.compounds:
        while not os.path.exists(compound.compFile):
          (dirName,fileName) = os.path.split(compound.compFile)
          print dirName, fileName
          interaction = self.multiDialog.FileName(self.guiParent,file = 'auremol.comp', component = 'chemComps', selectionText = 'Select file for compound %s' % compound.label, format = self.format)
          compound.compFile = interaction.file

      self.measurementFile = self.MeasurementFileClass(self.fileName)
      self.measurementFile.read(metaFile = metaFile)

    except StandardError:

      errorMessage = traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
    
      self.messageReporter.showWarning("Warning"," Cannot read %ss for %s...:\n%s." % (self.measurementType,self.formatLabel,errorMessage),self.guiParent)
      self.measurementFile = None
      
      return traceback.format_exception(sys.exc_type,sys.exc_value,sys.exc_info()[2]) 

  #
  # Functions different to default functions in DataFormat
  #

  def createShift(self,resonance,chemShift):
    
    keywds = {}
    
    #
    # Also set ambiguity code if listed
    #
    
    if chemShift.ambCode != None:
      setUniqueAppData('AppDataInt',resonance,self.format,'ambiguityCode',chemShift.ambCode)
    
    #
    # Shift error could be missing...
    #
    
    if hasattr(chemShift,'error'):

      keywds['error'] = chemShift.error
      
    #
    # Create shift linked to resonance
    #

    self.convertCount[self.mainCode][1][self.valueCode] += 1

    return Nmr.Shift(self.measurementList, value = chemShift.value, resonance = resonance, **keywds)


  def setPeakIntensity(self):
  
    # PeakIntensity attributes
    # Intensity and volume are always set for an nmrView peak

    if self.rawPeak.volume != 0:
      if self.rawPeak.volumeError:
        error = self.rawPeak.volumeError
      else:
        error = None
    
      peakInt = Nmr.PeakIntensity(self.peak,value = self.rawPeak.volume, error = error, method = self.methods[self.format]['Volume'])
      peakInt.intensityType = volume_kw

    if self.rawPeak.intensity != 0:
      peakInt = Nmr.PeakIntensity(self.peak,value = self.rawPeak.intensity, method = self.methods[self.format]['Intensity'])
      peakInt.intensityType = height_kw

  def setSpectrumInfoFromPeakList(self):
    
    for i in range(0,len(self.dataDimRefs)):
      
      sw = self.peakFile.spectrumInfo['swsHz'][i]
      sf = self.peakFile.spectrumInfo['sfs'][i]
      numPoints = self.peakFile.spectrumInfo['numPoints'][i]
      
      # TODO USE OFFSET!!!!! 
      
      dataDimRef = self.dataDimRefs[i]
      expDimRef = dataDimRef.expDimRef
      freqDataDim = dataDimRef.dataDim
      
      expDimRefDefault = expDimRef.findFirstApplicationData(application = ccpNmr_kw, keyword = 'default')
      freqDataDimDefault = freqDataDim.findFirstApplicationData(application = ccpNmr_kw, keyword = 'default')
      
      #
      # Only reset if clear that experiment was created as default - redundant as checking for new datasource but...
      #
      
      if expDimRefDefault and expDimRefDefault.value and freqDataDimDefault and freqDataDimDefault.value:
        
        expDimRef.sf = sf
        expDimRef.baseFrequency = sf
        freqDataDim.numPointsOrig = numPoints
        freqDataDim.numPoints = numPoints
        freqDataDim.valuePerPoint = (2.0 * sw) / ((freqDataDim.isComplex + 1) * freqDataDim.numPointsOrig)
        expDimRefDefault.delete()
        freqDataDimDefault.delete()

  def setPeakExtras(self):
  
    self.peak.annotation = self.rawPeak.label
    self.peak.details = self.rawPeak.comment
    
    if self.rawPeak.probability and self.rawPeak.probability != -1.0:
      self.peak.figOfMerit = self.rawPeak.probability
  
    if self.rawPeak.assign:
      
      #
      # Make peakContribs for this one...
      #
      
      pass
      
      # TODO get link to atom info to set correct info here!
      
      #numContribs = len(self.rawPeak.assign[0].split(self.assignTagSep))
      
      #for i in range(numContribs):
      #  self.peakContribs.append(Nmr.PeakContrib(self.peak))
        
    if hasattr(self.rawPeak,'subPeaks'):
      # TODO HANDLE THIS!!
      pass

  def setPeakDim(self):
  
    dataDimRef = self.dataDimRefs[self.rawPeakDimIndex]

    self.peakDim = self.peak.findFirstPeakDim(dim = dataDimRef.dataDim.dim)

    self.peakDim.dataDimRef = dataDimRef

    self.peakDim.value = self.rawPeak.ppm[self.rawPeakDimIndex]
    
    if self.rawPeak.width[self.rawPeakDimIndex] != -1.0:
      #rawPeak.width[i] is in points?
      # TODO: make sure this is used correctly!!!
      self.peakDim.decayRate = self.rawPeak.width[self.rawPeakDimIndex] / dataDimRef.valuePerPoint

  def thisPeakValid(self):
    
    # TODO add some info here?
    if not self.rawPeak.ppm:
      return 0
    
    return 1

  def getPeakResNames(self):
  
    self.resNames = []
   
    # TODO need to handle this at one stage!!
    pass
    
    
  def getExportChainCode(self,exportChainCode,chain):
    
    # Returns original chain code - this format uses segId so OK if longer than 1 char
    
    return exportChainCode
  
  def getPresetChainMapping(self,chainList):
  
    return self.getMultiChainFormatPresetChainMapping(chainList)

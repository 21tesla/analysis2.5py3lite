

"""
======================COPYRIGHT/LICENSE START==========================

NmrPipeFormat.py: Contains functions specific to NmrPipe conversions.

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
from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccpnmr.format.converters.VarianFormat import hasVarianSwitchedPhaseOrder
from ccpnmr.format.converters.BrukerFormat import getDMXtypeSpectrometer
from ccpnmr.format.converters.BrukerFormat import getBrukerDQDdata

from ccpnmr.format.general.Util import getNominalRef

class NmrPipeFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'nmrPipe'

  def createProcParsFile(self):
  
    if self.verbose == 1:
      print "Writing conversion parameters to %s file %s" % (self.formatLabel,self.fileName)
      
    if not self.outputFile:
    
      # Need selection!
      self.messageReporter.showError("Error","No output file given for %s conversion parameters" % (self.format),parent = self.guiParent)
     
    #
    # Set conversion command
    #
    
    if self.expOrigin == 'varian':
      convCommand = 'var2pipe'

    elif self.expOrigin == 'bruker':
      convCommand = 'bruk2pipe'        

    else:
      self.messageReporter.showError("Error","%s conversion parameters unknown for %s" % (self.format,self.expOrigin),parent = self.guiParent)
      return False

    self.procParsFile = self.conversionParsIO.NmrPipeConversionFile(self.fileName,self.inputFile,self.outputFile,self.numDim,convCommand)

    return None

  def createProcScriptFile(self):
  
    if self.verbose == 1:
      print "Writing processing script to %s file %s" % (self.formatLabel,self.fileName)
      
    if not self.outputFile:
    
      # Need selection!
      self.messageReporter.showError("Error","No output file given for %s processing script" % (self.format),parent = self.guiParent)
    
    self.procScriptFile = self.processingScriptIO.NmrPipeProcessingScriptFile(self.fileName,self.inputFile,self.outputFile,self.numDim)

  #
  # Functions different to default functions in DataFormat
  #

  def setProcParsSpecific(self): 
    
    #
    # Set extra commands at beginning and end of conversion script
    #
    
    dqdData = None
    
    if hasVarianSwitchedPhaseOrder(self.rawData):
    
      self.procParsFile.addEndCommands.append("-aqORD%15s" % '1')

    if self.expOrigin == 'bruker':
    
      specModel = getDMXtypeSpectrometer(self.experiment.spectrometer)
    
      if specModel:
      
        self.procParsFile.addBeginCommands.append("-%s" % specModel)
    
        dqdData = getBrukerDQDdata(self.rawData)
    
        if dqdData:

          self.procParsFile.addBeginCommands.append("-dspfvs " + dqdData['DSPFVS'])
          self.procParsFile.addBeginCommands.append("-decim " + dqdData['DECIM'])

    #
    # Get values
    #
    
    for dataDimRef in self.dataDimRefs:
    
      expDimRef = dataDimRef.expDimRef
      expDim = expDimRef.expDim
      expDimIndex = expDim.dim - 1

      fidDataDim = self.rawData.findFirstDataDim(expDim = expDim)
      #freqDataDim = dataDimRef.dataDim

      # Get label
      applData = dataDimRef.findFirstApplicationData(application = self.format, keyword = 'label')

      if applData:
        label = applData.value
      else:
        # Only use first isotope in list for label: has to be 'normal' experiment
        label = expDimRef.isotopeCodes[0]

      # Get complexFactor (multiplies number of points)
      complexFactor = self.getComplexFactor(fidDataDim)

      # Get nominal ref value (ppm on expDimRef.sf) from datadimref
      nominalRefValue = getNominalRef(dataDimRef)

      # Default settings - this list has to be in order with self.procParsFile.dataRefs
      values = [fidDataDim.numPointsValid * complexFactor,fidDataDim.numPointsValid,'Complex',fidDataDim.spectralWidth,expDimRef.sf,nominalRefValue,label]

      # Set values
      for i in range(0,len(values)):
        value = values[i]
        self.procParsFile.dataValues[expDimIndex].append(self.conversionParsIO.NmrPipeConversion(self.procParsFile,i,value))

      # Change values if necessary
      if fidDataDim.isComplex == False:
        # Non-complex data
        self.procParsFile.setDataValue(expDimIndex,'AcquisitionMode','Real')

    # Change values if necessary
    if not self.rawData.isNormalStorage:
      
      method = self.rawData.storageDetails
      
      #Sensitivity enhanced data

      for dataDimRef in self.dataDimRefs:
    
        expDimRef = dataDimRef.expDimRef
        expDim = expDimRef.expDim
        expDimIndex = expDim.dim - 1
        
        fidDataDim = self.rawData.findFirstDataDim(expDim = expDim)
        
        if fidDataDim.findFirstApplicationData(application = 'allProcessing', keyword = method, value = 'on'):
          
          if method == "echoAntiEcho":
            self.procParsFile.setDataValue(expDimIndex,'AcquisitionMode','Echo-AntiEcho')

          elif method == "ranceKay":
            self.procParsFile.setDataValue(expDimIndex,'AcquisitionMode','Rance-Kay')
            
          method = None

      if method:
        print "  Error: Unknown storage method %s" % method

    if dqdData:
      # Bruker DQD: acquisition dim (assumed to be first) set to DQD
      self.procParsFile.setDataValue(0,'AcquisitionMode','DQD')

 
  def setProcScriptSpecific(self):
      
    for dataDimRef in self.dataDimRefs:
    
      expDimRef = dataDimRef.expDimRef
      expDim = expDimRef.expDim
      expDimIndex = expDim.dim - 1

      fidDataDim = self.rawData.findFirstDataDim(expDim = expDim)
      freqDataDim = dataDimRef.dataDim

      # Default settings
      values = [None,
                None,
                '-off 0.3333 -pow 2.0',
                "-size %d" % (freqDataDim.numPointsOrig),
                '',
                "-p0 %.1f -p1 %.1f" % (freqDataDim.phase0,freqDataDim.phase1)]

      # Set values
      for i in range(0,len(values)):
        value = values[i]
        self.procScriptFile.dataValues[expDimIndex].append(self.processingScriptIO.NmrPipeProcessingScript(self,i,value))

      # Change values if necessary
      if fidDataDim.isComplex == False:
        # Non-complex data
        self.procScriptFile.setDataValue(expDimIndex,'FourierTransform','-real')

      if fidDataDim.alternateSign == True:
        # Negate alternate points (spectrum shifted by sw/2)
        self.procScriptFile.setDataValue(expDimIndex,'NegateAlternate','-alt')

      if fidDataDim.negateImaginary == True:
        # Conjugate of data (flip axis direction)
        self.procScriptFile.setDataValue(expDimIndex,'Conjugate','-i')

      #
      # Reset window function if available
      #

      keyword = 'windowFunction'
      applData = freqDataDim.findFirstApplicationData(application = self.format, keyword = keyword)

      if applData:
        windowFunc = applData.value
        argString = " "
        findWord = keyword + windowFunc

        for applData in freqDataDim.findAllApplicationData(application = self.format):
          if applData.keyword.find(findWord) == 0:
            argString += "-" + applData.keyword[len(findWord):] + " " + str(applData.value) + " "

        argString = argString[:-1]

        self.procScriptFile.setDataValue(expDimIndex,'WindowFunction',argString, dataRefCode = windowFunc)

      #
      # Set deleteImaginary (is AppDataInt)
      #

      keyword = 'deleteImaginary'
      applData = freqDataDim.findFirstApplicationData(keyword = keyword)

      if applData:
        self.procScriptFile.deleteImaginary[expDimIndex] = applData.value

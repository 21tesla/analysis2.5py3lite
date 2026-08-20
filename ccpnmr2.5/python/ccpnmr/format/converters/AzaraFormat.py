
"""
======================COPYRIGHT/LICENSE START==========================

AzaraFormat.py: Contains functions specific to Azara conversions.

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
import copy

from memops.universal.Util import returnInt #, isBigEndian
from memops.universal.Util import getUpperPowerTwo, getPowerTwoExp

from ccpnmr.format.general.Constants import newline
from ccp.format.general.Constants import indent

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords
from ccpnmr.format.converters.BrukerFormat import getBrukerDQDdata
from ccpnmr.format.converters.VarianFormat import hasVarianSwitchedPhaseOrder

#
# Format specific keywords...
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['writeProcScript']['parFile'] = (None,False,'The file the processing script will use for writing spectrum parameters to.')

class AzaraFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'azara'
    self.IOkeywords = IOkeywords

  def getProcPars(self):
  
    # TODO HERE: have to figure out what to do if project file read...
  
    self.procParsFile = self.processingParsIO.AzaraProcessingParsFile(self.fileName)
    self.procParsFile.read()

    if self.verbose == 1:
      print "Reading processing pars from %s file %s" % (self.formatLabel,self.fileName)

  def createProcParsFile(self):
  
    #
    # Check for outputfile, command, ask if necessary!!
    #
  
    if self.verbose == 1:
      print "Writing processing parameters to %s file %s" % (self.formatLabel,self.fileName)
  
    self.procParsFile = self.processingParsIO.AzaraProcessingParsFile(self.fileName)


  def createProcScriptFile(self):
  
    if self.verbose == 1:
      print "Writing processing script to %s file %s" % (self.formatLabel,self.fileName)
      
    if not self.outputFile:
    
      # Need selection!
      self.messageReporter.showError("Error","No output file given for %s processing script" % (self.format),parent = self.guiParent)
    
    self.procScriptFile = self.processingScriptIO.AzaraProcessingScriptFile(self.fileName,self.inputFile,self.outputFile,self.parFile,self.numDim)


  #
  # Functions different to default functions in DataFormat
  #

  def setProcParsSpecific(self): 
    
    self.procParsFile.fPars['inputFile'] = self.inputFile
    self.procParsFile.fPars['numDim'] = self.numDim
    self.procParsFile.initDims()
    
    # Take care of integer/float setting (WARNING - ALSO WHEN READING INFO!!)
    
    if self.rawData.numberType == 'int':
      self.procParsFile.fPars['integer'] = True

    # Set extra commands at beginning and end of processing pars
    if self.expOrigin == 'varian': 
      if self.numDim == 3:
        if hasVarianSwitchedPhaseOrder(self.rawData):
          value = "3 2"
        else:
          value = "2 3"

      elif self.numDim == 2:
        value = "2"
      else:
        self.messageReporter.showError("Error","only up to 3D experiment provided for azara Varian processing",parent = self.guiParent)
        return False

      self.procParsFile.fPars['varian'] = value

    #
    # Set whether data needs to be swapped
    #

    swapApplData = self.rawData.findFirstApplicationData(application = 'allProcessing', keyword = 'swap')

    if swapApplData:
      if returnInt(swapApplData.value):  #Swap is a string
        self.procParsFile.fPars['swap'] = True

    # Take care of endian setting (WARNING - ALSO WHEN READING INFO!!)
    
    if self.rawData.isBigEndian:
      self.procParsFile.fPars['bigEndian'] = True
      
      # TODO: set based on default like in Wayne's datasource creation? basically isBigEndian() for non-swapped, 'not isBigEndian()' for swapped
    
    # Get values
    self.procParsFile.fPars['dim'] = True
    
    for dataDimRef in self.dataDimRefs:
    
      expDimRef = dataDimRef.expDimRef
      expDim = expDimRef.expDim
      expDimIndex = expDim.dim - 1

      fidDataDim = self.rawData.findFirstDataDim(expDim = expDim)

      # Get complexFactor (multiplies number of points)
      complexFactor = self.getComplexFactor(fidDataDim)

      self.procParsFile.fPars['numPoints'][expDimIndex] = fidDataDim.numPointsValid * complexFactor
      self.procParsFile.fPars['spectralWidth'][expDimIndex] = fidDataDim.spectralWidth
      self.procParsFile.fPars['spectrometerFreq'][expDimIndex] = expDimRef.sf
      self.procParsFile.fPars['refPpm'][expDimIndex] = dataDimRef.refValue
      self.procParsFile.fPars['refPoint'][expDimIndex] = dataDimRef.refPoint * complexFactor
      self.procParsFile.fPars['nucleus'][expDimIndex] = expDimRef.isotopeCodes[0]

    return None
 
  def setProcScriptSpecific(self):
    
    #
    # Determine acquisition dimension
    #
    
    acqExpDim = self.experiment.findFirstExpDim(isAcquisition = True)
    
    if not acqExpDim: 
      acqExpDim = self.experiment.sortedExpDims()[0]
    
    acquisitionDimIndex = list(self.experiment.sortedExpDims()).index(acqExpDim)
    acquisitionDim = self.rawData.findFirstDataDim(expDim = acqExpDim)  

    #
    # Loop over all dims
    #
      
    for dataDimRef in self.dataDimRefs:
    
      expDimRef = dataDimRef.expDimRef
      expDim = expDimRef.expDim
      expDimIndex = expDim.dim - 1

      fidDataDim = self.rawData.findFirstDataDim(expDim = expDim)
      freqDataDim = dataDimRef.dataDim

      # Default settings
      values = ['','complex','','','','','sinebell2 60','zerofill 1','fft','',"phase %.1f %.1f" % (freqDataDim.phase0,freqDataDim.phase1),'']

      # Set values
      for i in range(0,len(values)):
        value = values[i]
        self.procScriptFile.dataValues[expDimIndex].append(self.processingScriptIO.AzaraProcessingScript(self.procScriptFile,i,value))

      # Change values if necessary
      if self.expOrigin == 'varian':
        # Varian: to keep axis in line with nmrPipe processing
        self.procScriptFile.setDataValue(expDimIndex,'Varian','conjugate')

      if fidDataDim.isComplex == False:
        # Non-complex data
        self.procScriptFile.setDataValue(expDimIndex,'DataType','real')
        self.procScriptFile.setDataValue(expDimIndex,'FourierTransform','rft')

      if fidDataDim.alternateSign == True:
        # Negate alternate points (spectrum shifted by sw/2)
        self.procScriptFile.setDataValue(expDimIndex,'NegateAlternate',"real%smask_ppmm" % (newline + indent))

      if fidDataDim.negateImaginary == True:
        # Conjugate of data (flip axis direction)
        self.procScriptFile.setDataValue(expDimIndex,'Conjugate','conjugate')

      #
      # Setting zerofill value
      #

      fidComplexFactor = self.getComplexFactor(fidDataDim)
      freqComplexFactor = self.getComplexFactor(freqDataDim)
      ratioPowerTwo = (freqDataDim.numPoints * freqComplexFactor)/getUpperPowerTwo(fidDataDim.numPointsValid * fidComplexFactor)

      if ratioPowerTwo > 1:

        #
        # In Azara number of real points counts. If one has 2048 real points for complex FID
        # data, then zerofill 1 will make this 4096 real points during FT. If the imaginary
        # data is deleted (isComplex = false), one then ends up with 2048 real points... 
        #

        zeroFillValue = getPowerTwoExp(ratioPowerTwo) + 1# In Azara number real points counts...
        self.procScriptFile.setDataValue(expDimIndex,'Zerofill','zerofill %d' % zeroFillValue)


      #
      # Reset window function if available (note that setup different from nmrPipe!)
      #

      keyword = 'windowFunction'
      applData = freqDataDim.findFirstApplicationData(application = self.format, keyword = keyword)

      if applData:
        windowFunc = applData.value
        argString = windowFunc
        findWord = keyword + windowFunc

        for applData in freqDataDim.findAllApplicationData(application = self.format):
          if applData.keyword.find(findWord) == 0:
            argString += " " + str(applData.value)

        self.procScriptFile.setDataValue(expDimIndex,'WindowFunction',argString)

      #
      # Set deleteImaginary (note that setup different from nmrPipe!)
      #

      keyword = 'deleteImaginary'  
      applData = freqDataDim.findFirstApplicationData(keyword = keyword)

      if applData:
        deleteImaginary = applData.value

        if deleteImaginary:
          self.procScriptFile.setDataValue(expDimIndex,'Reduce','reduce')

    #
    # Change values if necessary
    #
    
    if not self.rawData.isNormalStorage:
      
      #
      # Sensitivity enhanced data
      #

      method = self.rawData.storageDetails
      sensEnhanced = False

      for dataDimRef in self.dataDimRefs:
    
        expDimRef = dataDimRef.expDimRef
        expDim = expDimRef.expDim
        expDimIndex = expDim.dim - 1
        
        fidDataDim = self.rawData.findFirstDataDim(expDim = expDim)
        
        if fidDataDim.findFirstApplicationData(application = 'allProcessing', keyword = method, value = 'on'):
          
          sensEnhanced = True
          dimRef = self.procScriptFile.dimRefs[expDimIndex]
          fidDataDim = self.rawData.findFirstDataDim(className = "FidDataDim",expDim = expDim)


      if sensEnhanced:
      
        if method == "echoAntiEcho":
          mode = 'Echo-AntiEcho'

          # Set multipliers for Echo-antiecho dims
          self.procScriptFile.setDataValue(expDimIndex,'Rescale',"exchange                 ! Echo-antiecho type data")
          self.procScriptFile.setDataValue(acquisitionDimIndex,'Rescale',"scale2 1 %-d 0 1         ! Echo-antiecho type data" % acquisitionDim.numPointsValid)

        elif method == "ranceKay":
          mode = 'Rance-Kay'

          # Set multipliers for Echo-antiecho dims
          self.procScriptFile.setDataValue(expDimIndex,'Rescale',"scale2 1 %-d 0 1         ! Rance-Kay type data" % fidDataDim.numPointsValid)

          # Set multiplication for all other dims
          for tempDataDimRef in self.dataDimRefs:
    
            tempExpDimRef = tempDataDimRef.expDimRef
            tempExpDim = tempExpDimRef.expDim
            tempExpDimIndex = tempExpDim.dim - 1

            if tempExpDimIndex != expDimIndex and tempExpDimIndex != acquisitionDimIndex:
              tempFidDim = self.rawData.findFirstDataDim(className = "FidDataDim",expDim = tempExpDim)
              self.procScriptFile.setDataValue(tempExpDimIndex,'Rescale',"scale2 1 %-d -1 0         ! Rance-Kay type data" % tempFidDim.numPointsValid)

      else:
        print "  Error: Unknown sensitivity enhancement method %s" % method
        mode = ''

      self.procScriptFile.addBeginCommands.append('interlace ' + dimRef + '         ! %s type data' % mode)
      self.procScriptFile.setDataValue(acquisitionDimIndex,'Rescale',"scale2 1 %-d 0 1         ! %s type data" % (acquisitionDim.numPointsValid,mode))

    #
    # Bruker DQD data
    #
    
    dqdData = getBrukerDQDdata(self.rawData)

    if dqdData:
      self.procScriptFile.setDataValue(acquisitionDimIndex,'BrukerAvance','avance '+ dqdData['DSPFVS'] + " " + dqdData['DECIM'])
      self.procScriptFile.setDataValue(acquisitionDimIndex,'BrukerAvancePhase','avance_phase')

"""
======================COPYRIGHT/LICENSE START==========================

NmrDrawFormat.py: Contains functions specific to NmrDraw conversions.

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

import copy

from memops.api import Implementation

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords, filteredNameInfo

from ccpnmr.format.general.Constants import volume_kw, height_kw
from ccpnmr.format.general.Constants import tagSep

from ccpnmr.format.general.Util import getResName, getResNameText

import ccp.api.nmr.Nmr as Nmr

#
# Additional IOkeywords definitions
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readPeaks']['usePointsInfo'] = (False,False,'If selected, information about the width and peak box will be directly imported as points. Only do this if your spectrum in the data model is exactly the same as the original NMRPIPE one!')
IOkeywords['writePeaks']['useAliasing'] = (True,False,'If deselected, will export the peak chemical shifts at their real unaliased positions (i.e. outside the spectrum limits)')

class NmrDrawFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'nmrDraw'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
   
    self.getPeaks = self.getPeaksGeneric
    self.createPeakFile = self.createPeakFileGeneric

  #
  # Deviations from generic import stuff
  #


  #
  # Functions different to default functions in DataFormat
  #

  def thisPeakValid(self):
  
    #
    # Only use real peaks: no noise or artefacts
    #
  
    if self.rawPeak.type != 1:
      return 0

    return 1

  def setRawPeakFileFormatSpecific(self):
    
    #
    # Try to get the relationship between points/shifts from the peak positions
    # This is later used to calculate the peak width and box width in ppm
    #
    
    self.pointToPpm = []

    refPeak = self.peakFile.peaks[0]
    peakNumDims = len(refPeak.points)    
    
    highValues = []
    lowValues = []
    
    for rawPeakDimIndex in range(peakNumDims):
      refValue = (refPeak.points[rawPeakDimIndex],refPeak.ppm[rawPeakDimIndex])
      highValues.append(refValue)
      lowValues.append(refValue)
    
    for self.rawPeak in self.peakFile.peaks[1:]:
    
      if not self.thisPeakValid():
        continue
      
      for rawPeakDimIndex in range(peakNumDims):
        if self.rawPeak.points[rawPeakDimIndex] > highValues[rawPeakDimIndex][0]:
          highValues[rawPeakDimIndex] = (self.rawPeak.points[rawPeakDimIndex],self.rawPeak.ppm[rawPeakDimIndex])
        elif self.rawPeak.points[rawPeakDimIndex] < lowValues[rawPeakDimIndex][0]:
          lowValues[rawPeakDimIndex] = (self.rawPeak.points[rawPeakDimIndex],self.rawPeak.ppm[rawPeakDimIndex])
    
    for rawPeakDimIndex in range(peakNumDims):
      
      pointToPpm = None
      
      if not highValues[rawPeakDimIndex] == lowValues[rawPeakDimIndex]:
        pointDiff = highValues[rawPeakDimIndex][0] - lowValues[rawPeakDimIndex][0]
        ppmDiff = highValues[rawPeakDimIndex][1] - lowValues[rawPeakDimIndex][1]
        
        # Since values are relative
        pointToPpm = abs(ppmDiff / pointDiff)
      
      self.pointToPpm.append(pointToPpm)
    
  def setPeakIntensity(self):
    
    #
    # PeakIntensity attributes
    #
  
    if self.rawPeak.intensity != 0:
    
      peakInt = Nmr.PeakIntensity(self.peak,value = self.rawPeak.intensity, error = self.rawPeak.intensityError, method = self.methods[self.format]['Intensity'])
      peakInt.intensityType = height_kw
      
    if self.rawPeak.volume != 0:
    
      peakInt = Nmr.PeakIntensity(self.peak,value = self.rawPeak.volume, method = self.methods[self.format]['Volume'])
      peakInt.intensityType = volume_kw


  def setPeakExtras(self):
  
    #
    # chiProb is 1 for pure noise in nmrDraw
    #
  
    self.peak.figOfMerit = 1 - self.rawPeak.chiProb
    
    clustIdAppData = Implementation.AppDataInt(application=self.format,keyword='clustId',value=self.rawPeak.clustId)
    self.peak.addApplicationData(clustIdAppData)

    numPeaksClustAppData = Implementation.AppDataInt(application=self.format,keyword='numPeaksClust',value=self.rawPeak.numPeaksClust)
    self.peak.addApplicationData(numPeaksClustAppData)

  def setPeakDim(self):

    dataDimRef = self.dataDimRefs[self.rawPeakDimIndex]

    self.peakDim = self.peak.findFirstPeakDim(dim = dataDimRef.dataDim.dim)

    self.peakDim.dataDimRef = dataDimRef


    self.peakDim.value = self.rawPeak.ppm[self.rawPeakDimIndex]
    
    #
    # Use the points/ppm conversion from the peak positions, unless the 
    # usePointsInfo flag was set - only useful if spectrum info correct
    #
    
    boxWidth = self.rawPeak.limits[self.rawPeakDimIndex][1] - self.rawPeak.limits[self.rawPeakDimIndex][0]

    if self.usePointsInfo:
    
      self.peakDim.positionError = self.rawPeak.pointsDev[self.rawPeakDimIndex]      
      self.peakDim.boxWidth = boxWidth
      self.peakDim.decayRate = self.rawPeak.fwhm[self.rawPeakDimIndex]

    elif self.pointToPpm[self.rawPeakDimIndex]:
      
      self.peakDim.valueError = self.pointToPpm[self.rawPeakDimIndex] * self.rawPeak.pointsDev[self.rawPeakDimIndex]

      boxWidthPpm = self.pointToPpm[self.rawPeakDimIndex] * boxWidth
      self.peakDim.boxWidth = boxWidthPpm / self.dataDimRefs[self.rawPeakDimIndex].valuePerPoint
      
      decayRatePpm = self.pointToPpm[self.rawPeakDimIndex] * self.rawPeak.fwhm[self.rawPeakDimIndex]
      self.peakDim.decayRate = decayRatePpm / self.dataDimRefs[self.rawPeakDimIndex].valuePerPoint

  def setPeakFileInfo(self):
      
    self.peakFile.specNames = [self.dataSource.experiment.name]
    self.peakFile.setDimensions(self.dataSource.experiment.numDim)
    
    self.peakFile.refDimCodes = 'XYZA'[:self.dataSource.experiment.numDim]

    self.peakFile.refDimInfo['Label'] = []
    self.peakFile.refDimInfo['highPoint'] = []
    self.peakFile.refDimInfo['lowPoint'] = []
    self.peakFile.refDimInfo['highPpm'] = []
    self.peakFile.refDimInfo['lowPpm'] = []
    
    for dimIndex in range(0,len(self.dataDimRefs)):
      
      dataDimRef = self.dataDimRefs[dimIndex]
      expDimRef = dataDimRef.expDimRef
      dataDim = dataDimRef.dataDim
      
      self.peakFile.refDimInfo['Label'].append(expDimRef.isotopeCodes[0])
      self.peakFile.refDimInfo['highPoint'].append(1)
      self.peakFile.refDimInfo['lowPoint'].append(dataDim.numPoints)
      
      highPpm = dataDimRef.pointToValue(1)
      lowPpm = dataDimRef.pointToValue(dataDim.numPoints)

      self.peakFile.refDimInfo['highPpm'].append(highPpm)
      self.peakFile.refDimInfo['lowPpm'].append(lowPpm)

  def setPeakFilePeakExtras(self):
  
    self.points = []
    self.pointsDev = []
    self.ppm = []
    self.hz = []
    self.fwhm = []
    self.fwhm_Hz = []
    self.limits = []
    
    self.chiProb = 0.0    
    self.type = 1
    self.assignment = ''
    self.numPeaksClust = 0
    self.clustId = 0

    clustIdAppData = self.peak.findFirstApplicationData(application=self.format,keyword='clustId')
    if clustIdAppData:
      self.clustId = clustIdAppData.value

    numPeaksClustAppData = self.peak.findFirstApplicationData(application=self.format,keyword='numPeaksClust')
    if numPeaksClustAppData:
      self.numPeaksClust = numPeaksClustAppData.value
    
  def setPeakFilePeakIntensity(self):
    
    self.intensity = 0.0
    self.intensityError = 0.0
    self.volume = 0.0

    if self.peak.peakIntensities:
      
      # Only use obvious intensity/volume. Ignore otherwise (can be redone
      # in nmrView after all...
      
      for peakIntensity in self.peak.peakIntensities:
        if peakIntensity.intensityType == volume_kw:
          self.volume = peakIntensity.value
        elif peakIntensity.intensityType == height_kw:
          self.intensity = peakIntensity.value
          
          self.intensityError = peakIntensity.error
       
  def setPeakFilePeakDimInfo(self):

    point = self.peakDim.position

    self.points.append(point)
    
    if self.peakDim.positionError:
      self.pointsDev.append(self.peakDim.positionError)
    else:
      self.pointsDev.append(2)
  
    if self.useAliasing:
      ppm = self.peakDim.dataDimRef.pointToValue(point)
    else:
      ppm = self.peakDim.value
          
    self.ppm.append(ppm)
    
    self.hz.append((point * self.dataDimRef.dataDim.spectralWidth) / self.dataDimRef.dataDim.numPoints)
        
    if self.peakDim.boxWidth:
      boxWidth = self.peakDim.boxWidth
    else:
      boxWidth = 1
      
    self.limits.append([int(point-boxWidth),int(point+boxWidth)])
      
    if self.peakDim.decayRate:
      fwhm_Hz = self.peakDim.decayRate
      fwhm = fwhm_Hz * self.dataDimRef.dataDim.numPoints / self.dataDimRef.dataDim.spectralWidth
    elif hasattr(self.peakDim,'lineWidth') and self.peakDim.lineWidth:
      fwhm_Hz = self.peakDim.lineWidth
      fwhm = fwhm_Hz * self.dataDimRef.dataDim.numPoints / self.dataDimRef.dataDim.spectralWidth      
    else:
      fwhm = 0
      fwhm_Hz = 0
    
    self.fwhm.append(fwhm)
    self.fwhm_Hz.append(fwhm_Hz)
    
    #
    # Transfer assignment - TODO: currently lose information on specific assignment (peakContrib) items here!
    #

    peakDimContribs = self.peakDim.peakDimContribs

    resNameText = ""
    
    if self.writeAssignments and peakDimContribs:
          
      resNames = []
      resonanceToAtoms = []
      
      for peakDimContrib in peakDimContribs:
      
        if isinstance(peakDimContrib,Nmr.PeakDimContribN):
          print "   Warning: cannot handle multiple resonances for one peakDim contribution. Ignored."
          continue
     
        # get resonanceToAtoms
        ll = self.resonanceToAtoms.get(peakDimContrib.resonance)
        if ll is None:
          print ("  Can't handle resonance %s for peak %d" 
                 % (getResNameText(peakDimContrib.resonance), 
                 self.peak.serial))
        else:
          resonanceToAtoms.extend(ll)
     
      #  sort, filter, and transform
      if self.presetResonanceMapping:
        nameInfo = filteredNameInfo(resonanceToAtoms, mode='preset')
      elif self.compressResonances:
        nameInfo = filteredNameInfo(resonanceToAtoms, mode='compress')
      else:
        nameInfo = filteredNameInfo(resonanceToAtoms)
      
      # make names

      #
      # Don't use chainCode for nmrDraw output
      # Using chain and seqCode mapping for output
      # Note that code same whether original res names
      # are being used or not...
      #
     
      for resonanceToAtom, atomName in nameInfo:
              
            chain = resonanceToAtom.chain
            seqId = resonanceToAtom.seqId
              
            seqCode = self.getExportSeqCode(self.chainDict[chain][1],chain.findFirstResidue(seqId = seqId))
            
            resName = str(seqCode) + tagSep + atomName
          
            if resName not in resNames:
              resNames.append(resName)
              
      # Make resNameText
      resNameText = ','.join(resNames)
      
      
      ##
      ## Use atomSet if same for all (and more than one different contrib)
      ##
    
      #if (len(resNames) > 1 and len(resAtomSetNames) == 1 
      #    and None not in resAtomSetNames and self.compressResonances):

      #  resNameText = string.join(resAtomSetNames,',')
      
      #else:
    
      #  resNameText = string.join(resNames,',')
          
    self.assignment += resNameText + ';'
    
  
  def createPeakFilePeak(self):
    
    self.peakFile.peaks.append(self.peaksIO.NmrDrawPeak(self.peakFile))
    
    peak = self.peakFile.peaks[-1]
    
    peak.num = self.peakNum
    
    if self.assignment == '':
      self.assignment = 'None'
    else:
      self.assignment = self.assignment[:-1]
    
    for parName in ['points','pointsDev','ppm','hz','fwhm','fwhm_Hz',
                    'limits','chiProb','type','assignment','clustId',
                    'numPeaksClust','intensity','intensityError','volume']:
      
      value = getattr(self,parName)
      if value != None or (type(value) == type([]) and value.count(None) == len(value)):
        setattr(peak,parName,value)


  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)

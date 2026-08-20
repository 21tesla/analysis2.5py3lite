"""
======================COPYRIGHT/LICENSE START==========================

NmrViewFormat.py: Contains functions specific to NmrView conversions.

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

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords
from ccpnmr.format.converters.DataFormat import filteredNameInfo

from ccpnmr.format.general.Constants import tagSep
from ccpnmr.format.general.Constants import defaultMolCode
from ccpnmr.format.general.Constants import volume_kw, height_kw, ccpNmr_kw
from ccpnmr.format.general.Util import getResName, getResNameText

from ccp.lib.nmr.Nmr.DataSource import getXEasyDimCodes

import ccp.api.nmr.Nmr as Nmr

#
# Add some information to IOkeywords...
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readPeaks']['assignTagSep'] = (' ',False,"Separator for ambiguous assignments (e.g. the ' ' in '3.HA 4.HA').")
IOkeywords['readPeaks']['processingParsFile'] = (None,False,"NmrView .par file for setting parameters during experiment creation.")
IOkeywords['writePeaks']['writeMissingPeaks'] = (True,False,'Write out peaks with status -1 in between\n peaks with non-sequential peak numbers.')
IOkeywords['writePeaks']['assignTagSep'] = (' ',False,"Separator for ambiguous assignments (e.g. the ' ' in '3.HA 4.HA').")

class NmrViewFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'nmrView'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    self.createSequenceFile = self.createSequenceFileGeneric
        
    #self.getPeaks = self.getPeaksGeneric
    self.createPeakFile = self.createPeakFileGeneric
    
    self.getMeasurements = self.getMeasurementsGeneric
    self.createMeasurementFile = self.createMeasurementFileGeneric

  #
  # Deviations from generic import stuff
  #

  def createPeakFileFormatSpecific(self):

    self.peakFile.setSpectrumInfo(self.peakLists[0].name,self.numPeakDim)

  def getProcPars(self):
  
    # TODO HERE: have to figure out what to do if project file read...
  
    self.procParsFile = self.processingParsIO.NmrViewProcessingParsFile(self.fileName)
    self.procParsFile.read()

    if self.verbose == 1:
      print "Reading processing pars from %s file %s" % (self.formatLabel,self.fileName)

  def getPeaks(self):
  
    if self.processingParsFile and not self.dataSource:
    
      #
      # Read in spectrum parameters file if available
      #
      
      self.peaksFileName = self.fileName
      self.fileName = self.processingParsFile
      
      self.getProcPars()

      self.file = 'temporary'
      
      #
      # Set flags
      #
      
      readKeywds = {}
      
      for readKeywd in ('minimalPrompts',):
        
        value = getattr(self,readKeywd)
        
        if value:
          readKeywds[readKeywd] = value
      
      #
      # Create Experiment and DataSource in the data model
      #
      
      self.dataSource = self.readProcPars(**readKeywds)
      self.file = None
      self.fileName = self.peaksFileName

    if self.verbose == 1:
      print "Reading peak list from %s file %s" % (self.formatLabel,self.fileName)

    # TODO HERE: have to figure out what to do if project file read...
   
    # Assuming only one peaklist per file  

    self.peakFile = self.PeakFileClass(self.fileName, assignTagSep = self.assignTagSep)
    self.peakFile.read()
    
  def getFullProject(self,fileName,peakKeyWdList = None, sequenceKeyWds = None):
  
    if self.verbose == 1:
      print "Reading %s star project from file %s" % (self.formatLabel,fileName)
   
    self.file = self.projectIO.NmrViewStarFile(fileName)
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
    
    for i in range(0,len(self.file.peakFiles)):
      
      self.peakFile = self.file.peakFiles[i]

      if peakKeyWdList:
        keywds = peakKeyWdList[i]
      else:
        keywds = {}

      self.readPeaks(minimalPrompts = self.minimalPrompts, **keywds)
    
    #
    # Avoid reading project again
    #
    
    self.file = None
    
    return self.project
    

  #
  # Functions different to default functions in DataFormat
  #
    
  def createSequence(self):
  
    self.sequenceFile.sequences.append(self.sequenceIO.NmrViewSequence(molName = self.chain.molecule.name))
    self.sequence = self.sequenceFile.sequences[-1]
    
  def setSequenceFileElements(self):
  
    self.sequence.elements.append(self.sequenceIO.NmrViewSequenceElement(self.seqCode,self.residue.molResidue.chemComp.ccpCode))

  def setPeakIntensity(self):
  
    # PeakIntensity attributes
    # Intensity and volume are always set for an nmrView peak

    if self.rawPeak.volume != 0:
      peakInt = Nmr.PeakIntensity(self.peak,value = self.rawPeak.volume, method = self.methods[self.format]['Volume'])
      peakInt.intensityType = volume_kw

    if self.rawPeak.intensity != 0:
      peakInt = Nmr.PeakIntensity(self.peak,value = self.rawPeak.intensity, method = self.methods[self.format]['Intensity'])
      peakInt.intensityType = height_kw

  def setSpectrumInfoFromPeakList(self):
    
    for i in range(0,len(self.dataDimRefs)):
      
      sw = self.peakFile.sws[i]
      sf = self.peakFile.sfs[i]
      label = self.peakFile.dimCodes[i]
      
      dataDimRef = self.dataDimRefs[i]
      expDimRef = dataDimRef.expDimRef
      freqDataDim = dataDimRef.dataDim
      
      expDimRefDefault = expDimRef.findFirstApplicationData(application = ccpNmr_kw, keyword = 'default')
      freqDataDimDefault = freqDataDim.findFirstApplicationData(application = ccpNmr_kw, keyword = 'default')
      
      #
      # Only reset if clear that experiment was created as default
      #
      
      if expDimRefDefault and expDimRefDefault.value and freqDataDimDefault and freqDataDimDefault.value:
        
        expDimRef.sf = sf
        expDimRef.baseFrequency = sf
        freqDataDim.valuePerPoint = (2.0 * sw) / ((freqDataDim.isComplex + 1) * freqDataDim.numPointsOrig)
        
        expDimRef.removeApplicationData(expDimRefDefault)
        freqDataDim.removeApplicationData(freqDataDimDefault)

  def setPeakExtras(self):
  
    if self.peakFile.hasAssignItems:
      
      #
      # Find maximum number (could in principle be different?)
      #
      
      maxNumContribs = 0
      
      for assignIndex in range(len(self.rawPeak.assign)):
        if self.rawPeak.assign[assignIndex]:
          numContribs = len(self.rawPeak.assign[assignIndex].split(self.assignTagSep))
        else:
          numContribs = 0
        
        if maxNumContribs < numContribs:
          maxNumContribs = numContribs
      
      for i in range(maxNumContribs):
        self.peakContribs.append(Nmr.PeakContrib(self.peak))

  def setPeakDim(self):
  
    # TODO: rawPeak.shape[i]
    dataDimRef = self.dataDimRefs[self.rawPeakDimIndex]

    self.peakDim = self.peak.findFirstPeakDim(dim = dataDimRef.dataDim.dim)

    self.peakDim.dataDimRef = dataDimRef

    self.peakDim.value = self.rawPeak.ppm[self.rawPeakDimIndex]

    self.peakDim.boxWidth = self.rawPeak.box[self.rawPeakDimIndex] / dataDimRef.valuePerPoint
        
    #rawPeak.width[i] is related to decayrate, HERE GIVEN IN POINTS!!!!
    # TODO: make this the correct measurement at one point...
    if self.rawPeak.width[self.rawPeakDimIndex] and dataDimRef.valuePerPoint:
      self.peakDim.decayRate = self.rawPeak.width[self.rawPeakDimIndex] / dataDimRef.valuePerPoint

  def thisPeakValid(self):
  
    if self.rawPeak.status == -1:

      return 0
      
    return 1
    
  def getPeakResNames(self):
  
    self.resNames = []

    if self.rawPeak.assign[self.rawPeakDimIndex] not in [None,'?','null','']:
      
      try:
        self.getSinglePeakResName(self.rawPeak.assign[self.rawPeakDimIndex])
        
      except:
        
        try:
          assignNames = self.rawPeak.assign[self.rawPeakDimIndex].split(self.assignTagSep)
          
          for assignName in assignNames:
            self.getSinglePeakResName(assignName)
        
        except:
          print "  Problems handling %s assignment name %s... ignored." % (self.format,self.rawPeak.assign[self.rawPeakDimIndex])

  def getSinglePeakResName(self,assignName):
  
    (seqCode,atomName) = assignName.split(tagSep)
    self.resNames.append(getResName(defaultMolCode,seqCode,atomName))

  def setPeakFileInfo(self):
      
    self.peakFile.specName = self.dataSource.experiment.name

    #
    # Set the NmrView dimension labels...
    #
    
    peakDimCodes = []
    peakDimInfo = []
    
    for dataDimRef in self.dataDimRefs:
      expDimRef = dataDimRef.expDimRef
      sw = dataDimRef.dataDim.spectralWidth
      isotopeCode = expDimRef.isotopeCodes[0]
      peakDimCodes.append(isotopeCode)
      peakDimInfo.append("Label for dimension %d (%s, %.1f)" % (expDimRef.expDim.dim,isotopeCode,sw))
      
    if self.useXeasyDimCodes:
      # NBNB Rasmus 5 Aug 2012. Replaced with XEASY style, so that CYANA can read it
      # This would work globally, but the if statement avoids breaking normal use.
      peakDimCodes = getXEasyDimCodes(self.dataSource)
    
    
    if not self.minimalPrompts:
    
      interaction = self.multiDialog.EntryList(

                                 self.guiParent,
                                 peakDimCodes,
                                 peakDimInfo,
                                 title = 'Select labels for %s' % self.formatLabel,
                                 topText = 'Select dimension labels for %s peak export' % self.formatLabel)
      
      if interaction.values:
        peakDimCodes = interaction.values
        
      
    #
    # Now set the rest...
    #
    
    for dimIndex in range(0,len(self.dataDimRefs)):
      
      expDimRef = self.dataDimRefs[dimIndex].expDimRef
      self.peakFile.dimCodes[dimIndex] = peakDimCodes[dimIndex]
      self.peakFile.sfs.append(expDimRef.sf)      
 
      dataDim = self.dataDimRefs[dimIndex].dataDim
      self.peakFile.sws.append(dataDim.spectralWidth)
      
    #
    # This is used to put peaks with -1 in between peaks with non-sequential numbering
    #
    
    self.previousPeakNum = None

  def setPeakFilePeakExtras(self):
  
    self.status = 1
    
    if self.peak.isDeleted:
      # This should never happen... 
      self.status = -1
         
    self.ppm = []
    self.assign = []
    self.width = []
    self.box = []
    self.shape = []

  def setPeakFilePeakIntensity(self):
    
    self.volume = 0.0
    self.intensity = 0.0
         
    # Only use obvious intensity/volume. Ignore otherwise (can be redone
    # in nmrView after all...
    
    for peakIntensity in self.peak.peakIntensities:
      if peakIntensity.intensityType == volume_kw:
        self.volume = peakIntensity.value
      elif peakIntensity.intensityType == height_kw:
        self.intensity = peakIntensity.value
       
  def setPeakFilePeakDimInfo(self):

    ppm = self.peakDim.getValue()
        
    if self.peakDim.boxWidth:
      boxWidth = self.peakDim.boxWidth * self.dataDimRef.valuePerPoint
    else:
      boxWidth = 0.05
      
    if self.peakDim.decayRate:
      decayRate = self.peakDim.decayRate  * self.dataDimRef.valuePerPoint  # TODO: this is NOT correct - shouldn't be stored in points in data model
    else:
      decayRate = 0.05
    
    self.ppm.append(ppm)
    self.box.append(boxWidth)
    self.width.append(decayRate)
    self.shape.append('?')             # TODO: if appldata for nmrView, get from there?

    #
    # Transfer assignment. Code by Gary Thompson, should handle NmrView correctly!
    #
    # Refactored by rasmus fogh, July 2012
    #
    
    resNameTexts = []
    peakContribs = self.peakDim.peak.sortedPeakContribs()
    
    if (self.writeAssignments 
        and not (self.skipMultiAssignments and len(peakContribs) > 1)):

      peakDim = self.peakDim.dim
      peakDimContribSets = []

      for peakContrib in peakContribs:

        if len(peakContrib.peakDimContribs) == 0:
          continue

        peakDimContribSet = []
        peakDimContribSets.append(peakDimContribSet)

        for peakDimContrib in peakContrib.peakDimContribs:
          if peakDimContrib.dim ==  peakDim:
            peakDimContribSet.append(peakDimContrib)
      
      
      # get and filter resonanceToAtoms
      for peakDimContribs in peakDimContribSets:
        # loop over contribs
        
        if not peakDimContribs:
          continue
        
        resonanceToAtoms = []
        for peakDimContrib in peakDimContribs:
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
 
        # now get names
        resNames = []
        for resonanceToAtom, atomName in nameInfo:
          chain = resonanceToAtom.chain
          seqId = resonanceToAtom.seqId
          seqCode = self.getExportSeqCode(self.chainDict[chain][1],chain.findFirstResidue(seqId = seqId))
          resName = str(seqCode) + tagSep + atomName
          resNames.append(resName)
        resNames.sort()
        
        if not (len(resNames) > 1 and self.skipMultiAssignments):
        
          # Make resNameText
          resNameText = self.assignTagSep.join(resNames)

          if len(resNames) > 1:
          #if lenResNames > 1 and len(peakDimContribSets) > 1: # NBNB Rasmus try-out
            resNameText = "{%s}" % resNameText

          resNameTexts.append(resNameText)

    result  = " ".join(resNameTexts)

    if not result.strip():
      result = "?"
      
    self.assign.append(result)

  def createPeakFilePeak(self):
      
    #
    # To ensure that peaks are read in with the correct number, have to create
    # fake peaks...
    #
    
    if hasattr(self,'writeMissingPeaks') and self.writeMissingPeaks:
      
      if self.previousPeakNum != None and self.previousPeakNum != self.peakNum - 1:
      
        numDim = len(self.dataDimRefs)
        
        for tempPeakNum in range(self.previousPeakNum + 1, self.peakNum):
          
          self.peakFile.peaks.append(self.peaksIO.NmrViewPeak(tempPeakNum,
                                                         [''] * numDim,
                                                         [0.0] * numDim,
                                                         [0.0] * numDim,
                                                         [0.0] * numDim,
                                                         ['?'] * numDim,
                                                         0.0,
                                                         0.0,
                                                         -1,
                                                         isReading = False))
      
      self.previousPeakNum = self.peakNum    
    
    #
    # Now write the real peak...
    #

    self.peakFile.peaks.append(self.peaksIO.NmrViewPeak(self.peakNum,
                                                   self.assign,
                                                   self.ppm,
                                                   self.width,
                                                   self.box,
                                                   self.shape,
                                                   self.volume,
                                                   self.intensity,
                                                   self.status,
                                                   isReading = False))
                                                   
    
  def setChemShiftFileValue(self):
  
    self.setMeasurementFileGenericValue()

  def setMeasurementFileGenericValue(self):
  
    #
    # TODO: keep track of original code in nmrView (as applData to shift?)
    # 
    
    #
    # For NMRView all atomSet names should be written out in full...
    #
    
    if self.atomNamesDict.has_key(self.atomName):
    
      resonanceToAtom = self.atomNamesDict[self.atomName]
      measurement = self.origAtomMeasurements[resonanceToAtom]
      
      atomNames = []
      
      if resonanceToAtom.chemAtomSet and resonanceToAtom.chemAtomSet.isEquivalent:
      
        self.chemComp = resonanceToAtom.chemAtomSet.chemComp
        atomNames = self.getChemAtomSysNamesFromSet(resonanceToAtom.chemAtomSet,resonanceToAtom.getResidue().chemCompVar)
            
      if not atomNames:
      
        atomNames = [self.atomName]
        
      for atomName in atomNames:
    
        self.measurementFileValues.append(
      
               self.rawMeasurementClass(measurement.value,
               		                       atomName,
               		                       self.seqCode,
			                                     defaultMolCode))
      
      if self.atomMeasurements.has_key(resonanceToAtom):
        del(self.atomMeasurements[resonanceToAtom])
              		         
    #else:
      
    #  self.measurementFileValues.append(self.rawMeasurementClass(
    #		                            999.000,
    #          		                  self.atomName,
    #          		                  self.seqCode,defaultMolCode))

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)

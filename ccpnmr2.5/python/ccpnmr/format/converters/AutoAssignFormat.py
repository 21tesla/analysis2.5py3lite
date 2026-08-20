#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

AutoAssignFormat.py: Contains functions specific to AutoAssign conversions.

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

from ccpnmr.format.general.Constants import volume_kw, height_kw

from memops.api import Implementation
from ccp.api.nmr import Nmr

import copy, os

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readShifts']['getSeqFlag'] = (False,False,"Read sequence from the chemical shift file.")
IOkeywords['writeProject']['chain'] = (None,True,'The chain to be analyzed in AutoAssign.')

class AutoAssignFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'autoAssign'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    self.createSequenceFile = self.createSequenceFileGeneric
    
    self.getPeaks = self.getPeaksGeneric
    self.createPeakFile = self.createPeakFileGeneric
    
    self.getMeasurements = self.getMeasurementsGeneric

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

  def createPeakFileFormatSpecific(self):

    self.peakFile.setSpectrumInfo(self.peakLists[0].name,self.numPeakDim)

  def getFullProject(self,fileName, peakKeyWdList = None, sequenceKeyWds = None):
  
    if self.verbose == 1:
      print "Reading %s project from file %s" % (self.formatLabel,fileName)
   
    self.file = self.projectIO.AutoAssignProjectFile(fileName)
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

  def createFullProject(self,fileName, peakListInfo = None, tolerances = None, setTag = None):
  
    if self.verbose == 1:
      print "Writing %s project from file %s" % (self.formatLabel,fileName)
      
    self.file = self.projectIO.AutoAssignProjectFile(fileName)

    #
    # Need a GUI to set this up... but allow variables to be passed in as well, using:
    #
   
    if self.guiParent:
      from ccpnmr.format.gui.AutoAssignProjectInfoPopup import AutoAssignProjectInfoPopup
      interaction = AutoAssignProjectInfoPopup(self.guiParent,self.project,self.file)
      
      peakListInfo = interaction.peakListInfo
      tolerances = interaction.tolerances
    
    #
    # Quick check...
    #
    
    if not tolerances or not peakListInfo:
      self.messageReporter.showError("Error","Cannot write %s project in non-GUI enviroment without passing in peakListInfo list and tolerances dictionary." % self.format,parent = self.guiParent)
      return
    
    #
    # Set tag for export/import cycle
    #
    
    self.setDateTimeStamp(self.project,setTag,isUnique = False)
    self.setDateTimeStamp(self.chain,setTag,isUnique = True) # Only one export per chain!
    
    #    
    # Have to create an autoassign seq file with seq info, then all the peak lists...
    #
    
    self.writeSequence(fileName, chains = [self.chain], noWrite = True, resetIOkeywords = False)
    self.file.sequenceFile = self.sequenceFile
    
    self.file.tolerances = tolerances
    
    (baseName,projectFileName) = os.path.split(fileName)
    
    for (specType,peakList,phase) in peakListInfo:
      
      self.setDateTimeStamp(peakList,setTag,isUnique = False) # Spectra can be used for multiple chains...
      
      specFileName = os.path.join(baseName,"%s.pks" % specType)
      
      # For now assuming that only one dataDimRef per peak dimension... should be OK for sequential-type experiments
      # Otherwise have to do selection in export window...
      if not peakList.peaks:
        self.messageReporter.showError("Error","Cannot export peak list of type %s - no peaks available!" % specType,parent = self.guiParent)
        continue
    
      dataDimRefs = []
      for peakDim in peakList.sortedPeaks()[0].peakDims:
        dataDimRefs.append(peakDim.dataDimRef)
      
      
      writeKeywds = {}
      if specType == self.file.defaultSpectra[0][0]:
        writeKeywds['tagIoCycle'] = setTag

      self.writePeaks(specFileName, peakLists = [peakList], dataDimRefs = dataDimRefs, resetMapping = False, noWrite = True, resetIOkeywords = False,**writeKeywds)
      
      # Need to add some format-specific information
      dimInfo = []
      correction = 0
      
      for dataDimRef in dataDimRefs:
        
        freqDataDim = dataDimRef.dataDim
        expDimRef = dataDimRef.expDimRef
        
        atomInfo = []

        if expDimRef.refExpDimRef:
          # Doing some educated (hopefully) guessing...
          expMeasurement = expDimRef.refExpDimRef.expMeasurement
          for atomSite in expMeasurement.atomSites:
            if atomSite.name in tolerances.keys():
              atomInfo.append(atomSite.name)
            elif atomSite.name[0] == 'N':
              atomInfo.append('N15')
            elif atomSite.name[0] == 'H':
              transfers = atomSite.findAllExpTransfers(transferType = 'onebond')
              for transfer in transfers:
                tAtomSites = list(transfer.atomSites)
                otherAtomSite = tAtomSites[not tAtomSites.index(atomSite)]
                if otherAtomSite.name[0] == 'N':
                  atomInfo.append('HN')
                elif otherAtomSite.name == 'CA':
                  atomInfo.append('HA')
                elif otherAtomSite.name[0] == 'C':
                  atomInfo.append('HA')
                  atomInfo.append('HB') # TODO remove this one except if clear that CB involved?
                else:
                  continue
                break
            elif atomSite.name[0] == 'C':
              atomInfo.append('CA')
              atomInfo.append('CB')
            else:
              print "  Warning: no %s atom info set for atomSite '%s'." % (self.format,atomSite.name)
          
          if not atomInfo:
            print "  Warning: no atomSite info for expDimRef dim %d of reference experiment %s - guessing information from expDimRef." % (expDimRef.expDim.dim,expMeasurement.nmrExpPrototype.name)

        if not atomInfo:
          #
          # Guessing
          #
          if '1H' in expDimRef.isotopeCodes:
            if specType.count('HA') and dataDimRefs.index(dataDimRef) > 0:
              atomInfo.append('HA')
            else:
              atomInfo.append('HN')
          elif '13C' in expDimRef.isotopeCodes:
            for atomName in ('CA','CB','CO'): 
              if specType.count(atomName) > 0:
                atomInfo.append(atomName)
            if not atomInfo:
              print "  Warning: no specific 13C atom type found in experiment type '%s' - set to 'C13'." % specType
              atomInfo.append('C13')
          elif '15N' in expDimRef.isotopeCodes:
            atomInfo.append('N15')
          else:
            print "  Warning: no %s atom info set for expDimRef with isotope codes '%s'." % (self.format,str(expDimRef.isotopeCodes))
            atomInfo.append('X')
        
        if expDimRef.isFolded:
          foldStatus = 'folded'
        else:
          foldStatus = 'unfolded'
        
        if expDimRef.maxAliasedFreq:
          highPpm = expDimRef.maxAliasedFreq
        else: # TODO is this correct?
          #print dataDimRef.refPoint, dataDimRef.refValue
          highPpm = (dataDimRef.refPoint - 1) * dataDimRef.valuePerPoint + dataDimRef.refValue
          
        if expDimRef.minAliasedFreq:
          lowPpm = expDimRef.minAliasedFreq
        else: # TODO Is this correct?
            
          #print freqDataDim.numPoints, dataDimRef.valuePerPoint
          lowPpm = dataDimRef.refValue + (dataDimRef.refPoint - 1 - freqDataDim.numPoints) * dataDimRef.valuePerPoint
        
        #
        # Get the tolerances
        #
        
        tolerance = 0.001
        for atomName in atomInfo:
          if tolerances.has_key(atomName):
            if tolerances[atomName] > tolerance:
              tolerance = tolerances[atomName]
                   
        #
        # Track info per dimension
        #
        
        dimInfo.append((atomInfo,lowPpm,highPpm,correction,tolerance,foldStatus))
      
      for i in range(len(self.file.defaultSpectra)):
        if specType == self.file.defaultSpectra[i][0]:
          (specType,intra,sequential,throughSpace) = self.file.defaultSpectra[i]
          if i == 0:
            rootSpecName = 'ROOT'
            self.setIoCycleTag(peakList,str((setTag,'isRoot')),isUnique = False)
          else:
            rootSpecName = self.file.defaultSpectra[0][0]
      
      self.peakFile.setProjectSpectrumInfo(specType,rootSpecName,intra,sequential,throughSpace,phase,dimInfo)
 
      self.peakFile.write(**self.writeKeywds)
      self.file.peakFiles.append(self.peakFile)
      del self.peakFile # TODO Why this?

    
    self.file.write()
    
    return self.file
    
  #
  # Functions different to default functions in DataFormat
  #

  def createSequence(self):
  
    self.sequenceFile.sequences.append(self.sequenceIO.AutoAssignSequence(molName = self.chain.molecule.name))
    self.sequence = self.sequenceFile.sequences[-1]

  def setSequenceFileElements(self):
    
    self.sequence.elements.append(self.sequenceIO.AutoAssignSequenceElement(self.seqCode,self.getResidueCode1Letter(self.residue)))

  def setFitMethodSerial(self):
  
    self.peakList.fitMethodSerial = self.methods['unknown'][self.peakPickMethod]

  def setPeakIntensity(self):

    # PeakIntensity attributes
    if self.rawPeak.intensity != 0:
    
      # Software method not known: could come from anywhere
      peakInt = Nmr.PeakIntensity(self.peak,value = self.rawPeak.intensity, method = self.methods['unknown']['Intensity'])
    
      # Assuming height for autoassign
      peakInt.intensityType = height_kw
    
    # TODO: peak status!!!

  def setPeakDim(self):

    dataDimRef = self.dataDimRefs[self.rawPeakDimIndex]

    self.peakDim = self.peak.findFirstPeakDim(dim = dataDimRef.dataDim.dim)

    self.peakDim.dataDimRef = dataDimRef

    self.peakDim.value = self.rawPeak.ppm[self.rawPeakDimIndex]

  def setPeakFileInfo(self):
      
    self.peakFile.specName = self.dataSource.experiment.name

  def setPeakFilePeakExtras(self):
  
    self.status = self.peak.isDeleted
         
    self.ppm = []

  def setPeakFilePeakIntensity(self):
    
    self.intensity = 0.0
   
    if self.peak.peakIntensities:
      
      # Only use obvious intensity/volume... TODO: height has preference (??)
      
      for peakIntensity in self.peak.peakIntensities:
        if peakIntensity.intensityType == height_kw:
          self.intensity = peakIntensity.value
          
      if not self.intensity:
        for peakIntensity in self.peak.peakIntensities:
          if peakIntensity.intensityType == volume_kw:
            self.intensity = peakIntensity.value
       
  def setPeakFilePeakDimInfo(self):

    ppm = self.peakDim.getValue()
         
    self.ppm.append(ppm)
  
  def createPeakFilePeak(self):
      
    self.peakFile.peaks.append(self.peaksIO.AutoAssignPeak(self.peakNum,
                                                   self.ppm,
                                                   self.intensity))

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)

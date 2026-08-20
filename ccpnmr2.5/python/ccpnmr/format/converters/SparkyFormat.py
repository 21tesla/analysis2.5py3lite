"""
======================COPYRIGHT/LICENSE START==========================

SparkyFormat.py: Contains functions specific to Sparky conversions.

Copyright (C) 2007-2008 Wim Vranken (European Bioinformatics Institute)

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

import os, copy

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords, filteredNameInfo
from ccp.format.general.Util import getRegularExpressions
from ccp.general.Io import getDataSourceFileName

from ccpnmr.format.general.Util import getResName, getResNameText
from ccpnmr.format.general.Util import getNmrIsotope
from ccpnmr.format.general.Constants import volume_kw, height_kw, peakNum_kw
from ccpnmr.format.general.Constants import defaultMolCode

from memops.universal.Io import joinPath

import ccp.api.nmr.Nmr as Nmr
import memops.api.Implementation as Implementation

#
# Additional IOkeywords definitions
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readProject']['ignoreShifts'] = (None,False,'Ignore reading of chemical shift values.')
IOkeywords['writeProject']['entry'] = (None,False,'The entry to be exported to the Sparky file.')
IOkeywords['writeProject']['molSystem'] = (None,False,'The molecular system to be exported to the Sparky file.')
IOkeywords['writeProject']['useOriginalResNames'] = (False,False,'If set to True (on), will use the original\n resonance names for export (where possible).')
IOkeywords['writeProject']['useOriginalNumbers'] = (True,False,'If set to False (off), will not use the original\n peak numbers for export (if available).')
IOkeywords['writePeaks']['volumeIntensityType'] = ([height_kw,volume_kw],False,"Values to output to Sparky peak list. Give as list with values 'volume' or 'intensity', or both.")

class SparkyFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'sparky'
    self.IOkeywords = IOkeywords
    
    self.patt = getRegularExpressions(format = self.format)

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    self.createSequenceFile = self.createSequenceFileGeneric

    self.getPeaks = self.getPeaksGeneric
    self.createPeakFile = self.createPeakFileGeneric
    # NOTE: writing peak lists separately is not possible (?) - have to use
    #       UCSF proprietary format. Also use this for reading projects, ... !!!
    
    self.getMeasurements = self.getMeasurementsGeneric
    self.createMeasurementFile = self.createMeasurementFileGeneric

  #
  # Deviations from generic import stuff
  #

  def getFullProject(self,fileName,peakKeyWdList = None, sequenceKeyWds = None, shiftList = None):
  
    if self.verbose:
      print "Reading %s project from file %s" % (self.formatLabel,fileName)
   
    fileType = 'projectFile'
    self.file = self.projectIO.SparkyProjectFile(fileName)
    self.file.read(verbose = self.verbose)
    
    #
    # Detect if no information present - is probably saveFile so try reading it like this.
    #
    
    if not self.file.saveFiles and not self.file.sequenceFile and not self.file.chemShiftFiles:
      self.file.read(fileType = 'save', saveFileKey = fileName, verbose = self.verbose)
      if self.file.peaksLocation.has_key(fileName):
        fileType = 'saveFile'
        self.file.saveFiles.append(fileName)
        print "  Warning: file is of 'save' type - reading only single spectrum information, no sequences or chemical shifts."
    
    #
    # Process elements - call functions 'as if' reading real file, but pass
    # on objects created during project import
    #
    
    self.fileName = fileName
    
    self.sequenceFile = self.file.sequenceFile
 
    if self.sequenceFile:

      if sequenceKeyWds:
        keywds = sequenceKeyWds
      else:
        keywds = {}

      self.readSequence(minimalPrompts = self.minimalPrompts, **keywds)
    
    #
    # Select the peakFiles to load...
    #
    
    peakFiles = []
    peakFileDict = {}
    spectrumNames = []
    
    for saveFile in self.file.saveFiles:
    
      if self.file.peaksLocation.has_key(saveFile):
      
        for i in range(0,len(self.file.peaksLocation[saveFile])):
          peakFile = self.file.peaksLocation[saveFile][i]
          
          spectrumName = str(peakFile.specNames)[1:-1]
          peakFileDict[spectrumName] = peakFile
          peakFiles.append(peakFile)
          
          spectrumNames.append(spectrumName)
          
      else:      
        print "  Error: save file %s could not be read (not present in peaksLocation dictionary)" % saveFile
    
    #
    # Select the relevant spectra (assuming that names are unique!)
    #
    
    if not self.minimalPrompts:
    
      interaction = self.multiDialog.MultiSelectionList(

                         self.guiParent,
                         spectrumNames,
                         selectionDict = peakFileDict,
                         title = "Project '%s': " % self.project.name + 'Select Sparky peak lists to read in.',
                         text = "Spectrum names:",
                         dismissText = 'Select all'
                        )

      if interaction.isSelectedList != None:
        peakFiles = interaction.isSelectedList[:]
    
    #
    # Read selected peak files
    #

    for i in range(len(peakFiles)):
    
      self.peakFile = peakFiles[i]
      
      # TODO Probably better to do this by keyword...
      if peakKeyWdList:
        keywds = peakKeyWdList[i]
      else:
        keywds = {}

      self.readPeaks(minimalPrompts = self.minimalPrompts, **keywds)
    
    #
    # Read chemical shift files
    #
    
    (pathName,baseName) = os.path.split(self.fileName)

    if not self.ignoreShifts:
      for listSerial in range(len(self.file.chemShiftFiles)):
        chemShiftFile = self.file.chemShiftFiles[listSerial]
        self.measurementFile = chemShiftFile
        self.readShifts(minimalPrompts = self.minimalPrompts, measurementList = shiftList, measurementListName = "Sparky %s project list %d" % (baseName,listSerial))
    
    #
    # Avoid reading project again
    #
    
    self.file = None
    
    return self.project

  def createFullProject(self,fileName,peakKeyWdList = None):
    
    if not self.entry and not self.molSystem:
      print "  ERROR: need valid entry or molSystem to proceed with Sparky project file writing."
      return False
    
    if self.verbose == 1:
      print "Writing %s project to file %s" % (self.formatLabel,fileName)
    
    self.file = self.projectIO.SparkyProjectFile(fileName)
        
    if self.entry and self.entry.molSystem:
      # Entry takes precedence...
      self.molSystem = self.entry.molSystem 
      
    if not self.molSystem:
      print "  ERROR: need a valid molSystem from entry to proceed." 
      return False
    
    #
    # Default for now is to write it all for this project...
    # Am using the molSystem as Sparky 'molecule'
    #
    
    # First sort the chemical shift lists
    csListDict = {}
    for csList in self.project.currentNmrProject.findAllMeasurementLists(className = 'ShiftList'):
      csListDict[csList.serial] = csList

    serials = csListDict.keys()
    serials.sort()
    
    csLists = []
    for serial in serials:
      csLists.append(csListDict[serial])
      
    # Now loop over the shift lists in order... 
    
    csListNum = 0

    molName = self.molSystem.code
    self.file.molNames.append(molName)
    self.file.moleculeConditions[molName] = []
    self.file.chemShiftLocation[molName] = {}

    for csListNum in range(0,len(csLists)):

      conditionName = 'condition%d' % csListNum
      self.file.moleculeConditions[molName].append(conditionName)

      self.writeShifts('',measurementList = csLists[csListNum],chains = list(self.molSystem.sortedChains()), useOriginalResNames = self.useOriginalResNames, minimalPrompts = self.minimalPrompts, noWrite = True)
      self.file.chemShiftLocation[molName][conditionName] = self.measurementFile
    
    #
    # Can only get conditions if peak list connected to experiment!!
    # Otherwise I have to create conditions for everything?
    #
    
    expConditions = {}

    for csListNum in range(0,len(csLists)):

      conditionName = 'condition%d' % csListNum
      
      for exp in csLists[csListNum].experiments:
      
        if not expConditions.has_key(exp):
          expConditions[exp] = []
          
        expConditions[exp].append(conditionName)

    #
    # Now go over all experiments, see if there are any peaklists...
    # only write experiment if peak list available.
    #
    # There should be a good selection process for this... now writes it all
    # if no entry given...
    #
    # Could in principle also use Sample but that is not necessarily there either.
    #
    
    if self.entry:
    
      experimentList = self.entry.experiments
    
    else:
    
      experimentList = self.project.currentNmrProject.sortedExperiments()
        
    #
    # Loop over experiments
    #
    
    fileCount = 0
    
    for experiment in experimentList:
    
      if expConditions.has_key(experiment):
        
        #
        # Use first one in list...
        #
        
        conditionName = expConditions[experiment][0] 
        
      else:
      
        csListNum += 1
        conditionName = 'condition%d' % csListNum
        expConditions[experiment] = [conditionName]
        
      #
      # Loop over datasources
      #
    
      for ds in experiment.sortedDataSources():
      
        if not ds.peakLists:
        
          continue
          
        if ds.activePeakList:
        
          peakLists = [ds.activePeakList]
          
        else:
        
          peakLists = ds.sortedPeakLists()
        
        #
        # Loop over peaklists
        # 
        
        for peakList in peakLists:
           
          if peakKeyWdList:
            keywds = peakKeyWdList[fileCount]
          else:
            keywds = {}

          fileCount += 1

          self.writePeaks('',peakLists = [peakList], useOriginalResNames = self.useOriginalResNames, minimalPrompts = self.minimalPrompts, useOriginalNumbers = self.useOriginalNumbers, noWrite = True, **keywds)
          
          saveName = fileName + '.' + str(fileCount)
          if not self.peakFile.peaks:
            print "  Warning aborting export of peak list %s for Sparky - has no peaks." % saveName
            continue

          self.file.saveFiles.append(saveName)
          self.file.peaksLocation[saveName] = [self.peakFile]
          self.file.peaksCondition[self.peakFile] = [molName,conditionName]
           
    #
    # Write it all out...
    #
    
    self.file.write()
    
    return True

  #
  # Functions different to default functions in DataFormat
  #

    
  def createSequence(self):
  
    self.sequenceFile.sequences.append(self.sequenceIO.SparkySequence(molName = self.chain.molecule.name))
    self.sequence = self.sequenceFile.sequences[-1]
    
  def setSequenceFileElements(self):
        
    self.sequence.elements.append(self.sequenceIO.SparkySequenceElement(self.seqCode,code3Letter = self.residue.molResidue.chemComp.ccpCode))

  def setPeak(self):
  
    """
    
    Cannot overwrite for Sparky: have to delete all old peaks and create
    all new ones.
    
    """
    
    if self.rawPeak.num != None:
      details = "Original number %d" % self.rawPeak.num
    else:
      details = "No original number"
  
    self.peak = self.peakList.newPeak(details = details)
    self.convertCount[self.mainCode][1]['peaks'] += 1
    
    #
    # Also set original number in applicationdata so can be reused for writing
    #

    if self.rawPeak.num != None:

      self.peak.addApplicationData(Implementation.AppDataInt(application = self.format, keyword = peakNum_kw, value = self.rawPeak.num))

    # Sparky specific - has to be reset for every peak
    self.resLabelCodes = []

  def setPeakIntensity(self):
  
    # PeakIntensity attributes (only set if they were reported in the file)
    if self.rawPeak.intensity != None and (self.rawPeak.intensity or not self.rawPeak.volume):
      peakInt = Nmr.PeakIntensity(self.peak,value = self.rawPeak.intensity, method = self.methods[self.format]['Intensity (data)'])
      peakInt.intensityType = height_kw

    if self.rawPeak.intensityFit != None:
      peakInt = Nmr.PeakIntensity(self.peak,value = self.rawPeak.intensityFit, method = self.methods[self.format]['Intensity (fitted)'])
      peakInt.intensityType = height_kw

    if self.rawPeak.volume != None and (self.rawPeak.volume or not self.rawPeak.intensity):
      methodKey = "Volume (%s)" % self.rawPeak.volumeMethod
      peakInt = Nmr.PeakIntensity(self.peak,value = self.rawPeak.volume, method = self.methods[self.format][methodKey])
      peakInt.intensityType = volume_kw

  def setPeakExtras(self):
  
    # Use ppm if available, Hz if not
    if self.rawPeak.ppm:
      self.unit = 'ppm'
    elif self.rawPeak.hz:
      self.unit = 'hz'
    else:
      self.unit = None
      self.messageReporter.showError("Error","No peak position available from sparky peak file %s" % file,self.guiParent)
   

  def setPeakDim(self):
  
    if not self.unit:
      return

    dataDimRef = self.dataDimRefs[self.rawPeakDimIndex]

    self.peakDim = self.peak.findFirstPeakDim(dim = dataDimRef.dataDim.dim)

    self.peakDim.dataDimRef = dataDimRef
          
    if self.unit == 'ppm':
      self.peakDim.value = self.rawPeak.ppm[self.rawPeakDimIndex]

    elif self.unit == 'hz':
      # TODO: how can I define this?
      #(position,numAliasing) = hzToPoint(dataDimRef,rawPeak.hz[i])
      #peak.sortedPeakDims()[peakDimIndex].position = position
      #peak.sortedPeakDims()[peakDimIndex].numAliasing = numAliasing 
      return
      

  def getPeakResNames(self):
  
    self.resNames = []
    
    if self.rawPeak.assign:
    
      # Have to do some interpretation of the assigns...
      assign = self.rawPeak.assign[self.rawPeakDimIndex]
      
      # Special handling for Montelione group
      if '/' in assign:
        addAssignStrings = assign.split('/')
        assign = addAssignStrings.pop(0)
      else:
        addAssignStrings = []

      if assign != '?':

        # Ignore non-assigned peak dims
        assignString = self.patt[self.format + 'LabelCodeName'].search(assign)

        #
        # In Sparky with K18,19HA-R20HN style assignments multiple resNames
        # are possible here... making list of resNames to accomodate
        #

        if assignString:
          self.resLabelCodes = assignString.group(2).split(',')
          atomName = assignString.group(4)
          
        else:
          atomName = assign

        for resLabelCode in self.resLabelCodes: 
          self.resNames.append(getResName(defaultMolCode,resLabelCode,atomName))      
          
        # Add any other info (Montelione group again)
        for addAssignString in addAssignStrings:
          # This should always work in this case
          assignString = self.patt[self.format + 'LabelCodeName'].search(addAssignString)
          resName = getResName(defaultMolCode,assignString.group(2),assignString.group(4))
          self.resNames.append(resName)

  def setResNameExportFormat(self):
  
    #
    # Sparky is free format
    #
  
    if self.resNameFormat == None and self.useOriginalResNames == 1:
      self.resNameFormat = self.format
      
      
  def setPeakFileInfo(self):
  
    dataFileName = getDataSourceFileName(self.dataSource)
    self.peakFile.specFiles = [dataFileName]
      
    self.peakFile.specNames = [self.dataSource.experiment.name]
    self.peakFile.numDims = [len(self.dataSource.dataDims)]
    self.peakFile.specNumPoints = [[]]
      
    for dimIndex in range(0,len(self.dataDimRefs)):
      
      dataDim = self.dataDimRefs[dimIndex].dataDim
      self.peakFile.specNumPoints[0].append(dataDim.numPoints)

    #
    # Also do a check for which intensity type is available, reset if not good
    #
    
    intensityTypes = {height_kw: 0, volume_kw: 0}
    self.peakIntDict = {}
    
    for peak in self.peakList.peaks:
      
      self.peakIntDict[peak] = {'volume': None, 'volumeMethod': 'Unknown method', 'intensity': 0.0}#, 'intensityFit': 0.0}

      bestVolume = None
      bestVolumeMethod = None
      
      volume = None
     
      if peak.peakIntensities:
        
        # Only use obvious intensity/volume. Ignore otherwise (can be redone
        # after all...)
        
        for peakIntensity in peak.peakIntensities:
        
          if peakIntensity.intensityType == volume_kw:
            
            volume = peakIntensity.value
            method = peakIntensity.method
  
            if method.software.name == 'sparky':
  
              methodStringSearch = self.patt[self.format + 'BracketBetween'].search(method.procedure)
  
              if methodStringSearch:
                
                bestVolume = volume
                bestVolumeMethod = methodStringSearch.group(1)
  
            elif not bestVolumeMethod:
              
              volumeMethodParam = method.findFirstParameter(name = 'method')
              
              if volumeMethodParam:
                if volumeMethodParam.value == 'box sum':
                  bestVolumeMethod = 'Sum over box'
                # TODO This might not be correct!
                elif volumeMethodParam.value == 'truncated box sum':
                  bestVolumeMethod = 'Sum over ellipse'
                elif volumeMethodParam.value == 'parabolic fit':
                  bestVolumeMethod = 'Gaussian'
                  
                if bestVolumeMethod:
                  bestVolume = volume
              
          elif peakIntensity.intensityType == height_kw:
          
            method = peakIntensity.method
  
            if method.software.name == 'sparky':
  
              methodStringSearch = self.patt[self.format + 'BracketBetween'].search(method.procedure)
  
              if methodStringSearch.group(1) == 'data':
  
                self.peakIntDict[peak]['intensity'] = peakIntensity.value
                intensityTypes[height_kw] += 1
                
              elif  methodStringSearch.group(1) == 'fitted':
  
                self.peakIntDict[peak]['intensity'] = peakIntensity.value  # Used to be intensityFit
                intensityTypes[height_kw] += 1
            
            elif not self.peakIntDict[peak]['intensity']:
              self.peakIntDict[peak]['intensity'] = peakIntensity.value
              intensityTypes[height_kw] += 1
      
      #
      # Set volume method only if good info found!
      #
        
      if bestVolumeMethod:
        self.peakIntDict[peak]['volume'] = bestVolume
        self.peakIntDict[peak]['volumeMethod'] = bestVolumeMethod
        intensityTypes[volume_kw] += 1
      elif volume != None:
        self.peakIntDict[peak]['volume'] = volume
        intensityTypes[volume_kw] += 1

    #
    # First remove columns if no values
    #
       
    volumeIntensityType = self.volumeIntensityType[:]
    
    for intensityType in intensityTypes.keys():
      if intensityTypes[intensityType] == 0:
        if intensityType in volumeIntensityType:
          if self.verbose:
            print "  Warning: removing peak integration method %s - no values available." % (intensityType)
          volumeIntensityType.pop(volumeIntensityType.index(intensityType))
    
    #
    # Now set what has to be written out...
    #
    
    self.peakFile.volumeIntensityType = volumeIntensityType
          
  def setPeakFilePeakExtras(self):
         
    self.ppm = []
    self.assign = []
    self.widths = []

  def setPeakFilePeakIntensity(self):
    
    for attrName in self.peakIntDict[self.peak]:
      setattr(self,attrName,self.peakIntDict[self.peak][attrName])
      
  def setPeakFilePeakDimInfo(self):

    ppm = self.peakDim.getValue()
        
    if self.peakDim.decayRate:
      decayRate = self.peakDim.decayRate  * self.dataDimRef.valuePerPoint  # TODO: this is NOT correct - shouldn't be stored in points in data model
    else:
      decayRate = None
      
    self.ppm.append(ppm)
    self.widths.append(decayRate)

    #
    # Transfer assignment
    #

    peakDimContribs = self.peakDim.peakDimContribs

    resNameText = ""

    #if self.writeAssignments and peakDimContribs != ():
    if self.writeAssignments and peakDimContribs:
          
      resNames = []
      resonanceToAtoms = []
      
      # get and filter resonanceToAtoms
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
      
      for resonanceToAtom, atomName in nameInfo:

        #
        # Don't use chainCode for sparky output
        # Using chain and seqCode mapping for output
        # Note that code same whether original res names
        # are being used or not...
        #
          
        chain = resonanceToAtom.chain
        seqId = resonanceToAtom.seqId
          
        seqCode = self.getExportSeqCode(self.chainDict[chain][1],chain.findFirstResidue(seqId = seqId))
        residueLabel = self.getResidueLabelOneLetter(resonanceToAtom)
        
        resName = residueLabel + str(seqCode) + atomName
        
        if resName not in resNames:
          resNames.append(resName)

      if resNames: 
        resNameText = resNames[0]
      
        if len(resNames) > 1:   
          print "  Warning: Multiple assignments currently not supported for %s. Only using first one." % self.format
                
    self.assign.append(resNameText)
      
  def createPeakFilePeak(self):

    peak = self.peaksIO.SparkyPeak(self.peakFile,{})
    self.peakFile.peaks.append(peak)
    
    peak.num = self.peakNum
    peak.intensityFit = 0.0 #self.intensityFit
    peak.intensity = self.intensity
    peak.ppm = self.ppm
    peak.assign = self.assign[:]
    peak.volume = self.volume
    peak.volumeMethod = self.volumeMethod
    
    if None not in self.widths:
    
      self.peakFile.widths = self.widths[:]
 
  def setChemShiftFileValue(self):
          
    if self.atomNamesDict.has_key(self.atomName):
    
      resonanceToAtom = self.atomNamesDict[self.atomName]
      chemShift = self.origAtomMeasurements[resonanceToAtom]
         
      residueLabel = self.getResidueLabelOneLetter(resonanceToAtom)
      
      atom = resonanceToAtom.resonance.resonanceSet.findFirstAtomSet().findFirstAtom()
      atomType = getNmrIsotope(self.project,atom.chemAtom.elementSymbol)
      
      values = {'residue':   [residueLabel + str(self.seqCode)],
                'atom':      [self.atomName],
                'atomType':  [atomType],
                'chemShift': [str(chemShift.value)]
               }
    
      self.measurementFileValues.append(
          self.rawMeasurementClass(self.measurementFile,values)
                                       )
      
      if self.atomMeasurements.has_key(resonanceToAtom):
        del(self.atomMeasurements[resonanceToAtom])

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)

"""
======================COPYRIGHT/LICENSE START==========================

Conversion.py: Set of classes to simplify conversion between format files

Copyright (C) 2008-2011 Wim Vranken (Vrije Universiteit Brussel)

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

from memops.api import Implementation

from memops.universal.Util import drawBox
from memops.universal.Io import joinPath, convertStringToFileName

from memops.general.Util import returnMemopsWord

from ccp.general.Util import setCurrentStore

from ccpnmr.format.general.Constants import IoSetupList, projectExportRequirements
from ccpnmr.format.general.Constants import fileTypeDict
from ccpnmr.format.general.Util import createExperiment, createPpmFreqDataSource
from ccpnmr.format.general.Util import createSelection, getValidExportFormats

from ccpnmr.format.process.linkResonances import linkResonances

from ccpnmr.format.process.sequenceCompare import SequenceCompare

import os, time, sys, traceback

fileAttrNames = {

'peaks':  'peakFile',
'shifts': 'chemShiftFile', # Should this be measurementfile?
'sequence': 'sequenceFile',
'coordinates': 'coordinateFile',

'distanceConstraints':  'constraintFile',
'hBondConstraints':     'constraintFile',
'dihedralConstraints':  'constraintFile',
'rdcConstraints':       'constraintFile',
'jCouplingConstraints': 'constraintFile',

}

def getFullErrorString(errorMessage):
  
  returnStr = "\n" + drawBox("Error details") + "\n"
  returnStr += "".join(errorMessage)
  returnStr += "\n"
  
  return returnStr

class FormatConversionError(StandardError):

  def __init__(self, value):
    
    self.value = value
    
    if sys.exc_type:
      self.details = traceback.format_exception(sys.exc_type,sys.exc_value,sys.exc_info()[2]) 
    else:
      self.details = None
    
  def __str__(self):
    
    if self.value.count("\n"):
      returnStr = self.value
    else:
      returnStr = repr(self.value)
    
    if self.details:
      returnStr += getFullErrorString(self.details)
    
    return returnStr
 
class FormatConversion(object):

  """
  
  This class provides defaults for import/export of data, but can be customised.
  
  It works from one single project, and uses an NmrEntry as the link to all top objects.
  
  TODO: if, for example, required keyword for 'peaks' import missing, say that required, also list
        all available dataSources (if any) with the right 'key' to find them to link to the peaklist!
  
  
  Special (required) keywords:
  
  - 'peaks':  dataSourceSerial = <integer>    # Use if want to link directly to existing dataSource
  
              OR
              
              refExperimentIdString = <string>    # Use if want to create mock dataSource from scratch. Use 'old' experiment names,
                                                  # should be changed to BMRB ones!
                                                  
              OR

              refExperiment = <CCPN ReferenceExperiment>
  
  """

  def __init__(self,ccpnProject = None,
                    ccpnProjectName = None,
                    identifier = None,
                    useGui = False,
                    guiRoot = None,
                    resonancesLinked = False,
                    silentRead = True):
  
    self.formats = {}
    self.preparsedFiles = {}
    self.ccpnObjectDict = {}
    
    # This variable tracks if linkResonances was already run or not
    self.resonancesLinked = resonancesLinked
    
    if not identifier:
      identifier = time.strftime("%Yo%mo%do%Ho%Mo%S")
    
    if not ccpnProjectName:
      ccpnProjectName = identifier
    
    if ccpnProject:
      self.ccpnProject = ccpnProject
    else:
      validCcpnProjectName = returnMemopsWord(ccpnProjectName)
      
      if not validCcpnProjectName == ccpnProjectName:
        print "  Warning: Project name changed from '%s' to '%s' for CCPN validity" % (ccpnProjectName,validCcpnProjectName)
      
      self.ccpnProject = Implementation.MemopsRoot(name = validCcpnProjectName)
    
    #
    # Set the reference NmrProject
    #
    
    setCurrentStore(self.ccpnProject,'NmrProject')
    self.nmrProject = self.ccpnProject.currentNmrProject
    
    #
    # Set 'silent reading' - will not print out read messages when identifying a format based on dataType
    #
    
    self.silentRead = silentRead
    
    # Keep stdout setting to avoid breaking incoming stdout redirection
    self.stdstdout = sys.stdout
    self.stdstderr = sys.stderr
    
    #
    # Set reference NMR experiment information
    #
    
    self.setNmrRefExperimentInfo()
    
    #
    # Set the reference NmrEntry
    #
    
    setCurrentStore(self.ccpnProject,'NmrEntryStore')
    
    self.nmrEntry = self.ccpnProject.currentNmrEntryStore.findFirstEntry(name = identifier)
    
    if not self.nmrEntry:
      self.nmrEntry = self.ccpnProject.currentNmrEntryStore.newEntry(name = identifier)
    
    #
    # Set graphical interface, if required
    #
    
    if guiRoot:
      self.guiRoot = guiRoot
    elif useGui:
      import Tkinter 
      self.guiRoot = Tkinter.Tk()
    else:
      self.guiRoot = None      
      
    if self.guiRoot:
      from memops.gui.ProgressBar import ProgressBar
      self.ProgressBar = ProgressBar
    
    #
    # Set information types dictionary
    #
    
    self.informationTypeDict = {}
    
    for (informationType,readFunc,writeFunc,fileDict,importSetup,exportSetup,exportFunc,importFunc) in IoSetupList:
    
      self.informationTypeDict[informationType] = {}
      
      funcs = (readFunc,writeFunc)
      
      for formatName in fileDict.keys():
      
        self.informationTypeDict[informationType][formatName] = [None,None]
        
        # Check whether import/export possible at all
        for i in range(1,3):
          if fileDict[formatName][i]:
            self.informationTypeDict[informationType][formatName][i-1] = funcs[i-1]
  
  def getFormatClass(self,formatName):
  
    if not self.formats.has_key(formatName):
      formatClass = self.createFormatClass(formatName)
      self.formats[formatName] = formatClass
    else:
      formatClass = self.formats[formatName]
    
    return formatClass
    
  def createFormatClass(self,formatName):

    from ccpnmr.format.general.Constants import allFormatsDict
    
    if allFormatsDict.has_key(formatName):

      formatLabel = allFormatsDict[formatName]
      formatModule = __import__('ccpnmr.format.converters.%sFormat' % formatLabel,{},{},['%sFormat' % formatLabel]) # TODO: Also need IOKeywords?
    
      FormatClass = getattr(formatModule,'%sFormat' % formatLabel)
      formatClass = FormatClass(self.ccpnProject,self.guiRoot,verbose = 0)
      
    else:
      raise FormatConversionError("No '%s' format available!" % formatName)
 
    return formatClass
    
  def saveCcpnProject(self):
  
    self.ccpnProject.saveModified()
    
  def packCcpnProject(self):
    
    # TODO put in code to do a tar cvfz on the project
    # Get this from extend-NMR pipeline code!
    pass
    
  #
  # Simple import functions
  #
  
  def importFile(self, informationType, formatName, fileName, addKeywords=None, addObjectStrings=None, isMulti=True):
    
    if addKeywords == None:
      self.addKeywords = {}
    else:
      self.addKeywords = addKeywords
      
    self.convertObjectStrings(addObjectStrings)
       
    fileName = self.getImportFileName(informationType,fileName,isMulti=isMulti)

    self.executeIO('import',fileName,informationType,formatName)
    
    return self.conversionReturn
    
  def importFiles(self, informationType, formatName, fileNames, addKeywords=None, addObjectStrings=None):

    fileNames = self.getImportFileName(informationType,fileNames,isMulti=True)
    
    returnList = []
    
    for fileName in fileNames:
      # Make sure don't do multi again on files when importing!
      self.importFile(informationType,formatName,fileName,addKeywords=addKeywords,addObjectStrings=addObjectStrings,isMulti=False)
      returnList.append(self.conversionReturn)

    return returnList
        
  def getImportFileName(self,informationType,fileName,isMulti=False):

    if informationType == 'coordinates' and isMulti:
      # Note that if fileName is a list, will do import on this list at once when running self.importFile!
      fileName = (fileName,)

    return fileName

  #
  # Simple export functions
  #
      
  def exportFile(self,informationType,formatName,fileName,addKeywords = None, addObjectStrings = None):
    
    if addKeywords == None:
      self.addKeywords = {}
    else:
      self.addKeywords = addKeywords
     
    self.convertObjectStrings(addObjectStrings)
          
    self.executeIO('export',fileName,informationType,formatName)
    
  def exportFiles(self,informationType,formatName,fileNames,addKeywords = None, addObjectStrings = None):

    for fileName in fileNames:
      self.exportFile(fileName,informationType,formatName,addKeywords=addKeywords,addObjectStrings=addObjectStrings)
      
  def exportAllObjects(self,informationType,formatName,exportDir,addKeywords = None, addObjectStrings = None):
    
    if addKeywords == None:
      addKeywords = {}
  
    # TODO: Could in principle pass this in to avoid having to do all over again every time
    (exportFormats, ccpnObjects) = getValidExportFormats(self.ccpnProject) 
    
    if not exportFormats.has_key(formatName):
      print "   Warning: no information for %s. Aborting." % formatName
    elif informationType not in exportFormats[formatName]:
      print "   Warning: no information to export %s for %s. Aborting." % (informationType,formatName)
    else:
      
      ccpnObjectsLoop = []
      
      if informationType == 'peaks':
        ccpnObjectsLoop = [[peakList] for peakList in ccpnObjects['PeakList']]
        objectAddKeywordName = 'peakLists'
        
      elif informationType == 'sequence': # Or chain by chain?
        ccpnObjectsLoop = [molSystem.sortedChains() for molSystem in ccpnObjects['MolSystem']]
        objectAddKeywordName = 'chains'
        
      elif informationType == 'shifts':
        ccpnObjectsLoop = ccpnObjects['ShiftList']
        objectAddKeywordName = 'measurementList'
        
      elif informationType == 'coordinates':
        ccpnObjectsLoop = [strucEns.sortedModels() for strucEns in ccpnObjects['StructureEnsemble']]
        objectAddKeywordName = 'structures'
        
      elif informationType[-11:] == 'Constraints':
        constraintListClassName = "%s%sList" % (informationType[0].upper(),informationType[1:-1])
        ccpnObjectsLoop = ccpnObjects[constraintListClassName]
        objectAddKeywordName = 'constraintList'
        
      elif informationType == 'project':
        print "  Error: cannot export project data in exportAllObjects."
      else:
        print "  Error: unknown information type %s for exportAllObjects." % informationType
      
      #
      # Now do the export.
      #
      
      exportFileNames = []
      
      for i in range(len(ccpnObjectsLoop)):
        ccpnObjectOrList = ccpnObjectsLoop[i]
        addKeywords[objectAddKeywordName] = ccpnObjectOrList
        
        #
        # Custom set suffix and change middle part if required. Not handling anything more complicated at the moment.
        #
        
        fileType = ".%s" % informationType
        suffix = 'export'
        
        if fileTypeDict.has_key(informationType):
          if fileTypeDict[informationType].has_key(formatName):
            filePattern = fileTypeDict[informationType][formatName][0]
            if filePattern.count('*') == 1 and filePattern[:2] == '*.':
              suffix = filePattern[2:]
              
              if suffix == informationType[:len(suffix)]:
                fileType = ""
              elif suffix == 'pks' and informationType == 'peaks':
                fileType = ""
              elif suffix == 'ppm' and informationType == 'shifts':
                fileType = ""
        
        #
        # Custom set filename
        #
        # Note: using serial to distinguish files in some cases, so if there are empty peakLists (or constraint lists),
        # the number in the fileName will be higher than 1, even if there is only one valid list!
        #
        
        identifier = str(i+1)
        if informationType == 'peaks':
          if len(ccpnObjectOrList) == 1:
            peakList = ccpnObjectOrList[0]
            experiment = peakList.dataSource.experiment
            identifier = experiment.name       

            dataSources = experiment.findAllDataSources(dataType = 'processed')
            if len(dataSources) == 1:
              if len(dataSources.pop().peakLists) != 1:
                identifier += ".%d" % peakList.serial
            else:
              identifier += ".%d.%d" % (dataSource.serial,peakList.serial)
        
        elif informationType == 'sequence':
          identifier = ccpnObjectOrList[0].molSystem.code
        
        elif informationType == 'shifts' or informationType[-11:] == 'Constraints':
          identifier = str(ccpnObjectOrList.serial)
            
        identifier = '.' + identifier
        
        #
        # Make sure filename is OK, also make sure it is unique!
        #
        
        fileName = convertStringToFileName("%s%s%s.%s" % (formatName,fileType,identifier,suffix))
        
        addId = 1
        while fileName in exportFileNames:
          fileName = convertStringToFileName("%s%s%s.%d.%s" % (formatName,fileType,identifier,addId,suffix))
          addId += 1
        exportFileNames.append(fileName)
        
        filePath = joinPath(exportDir,fileName)        
 
        self.exportFile(informationType,formatName,filePath,addKeywords=addKeywords,addObjectStrings=addObjectStrings)
      
  #
  # Pre-parsing functions - try to get relevant information from the external file
  #
  # TODO: might want to have this in DataFormat.py for automatic peakDim matching?!
  #
  
  def preparseFile(self,informationType,formatName,fileName,addKeywords = None):
    
    formatClass = self.getFormatClass(formatName)
    
    fileName = self.getImportFileName(informationType,fileName,isMulti=True)
    formatClass.fileName = fileName
    
    defaultIODict = {'minimalPrompts': True}
     
    try:
      
      if informationType == 'peaks':

        formatClass.setIOkeywords('readPeaks', defaultIODict, verbose = False)
        
        formatClass.setPeakFileClass()
        formatClass.getPeaks()
        formatClassFile = formatClass.peakFile
           
        #
        # Get information for the peak dimensions
        #
        
        if not formatClassFile.peaks:
          returnInfo = (False,"Input not valid - no peak file or peaks found.")
        else:
          returnInfo = (True,self.getPeakDimInfo(formatClass))
  
      elif informationType == 'sequence':
        
        formatClass.setIOkeywords('readSequence', defaultIODict, verbose = False)

        formatClass.setSequenceFileClass()
        formatClass.getSequence()
        formatClassFile = formatClass.sequenceFile
        
        # Remove sequences with no information.
        for sequence in formatClassFile.sequences[:]:
          if not sequence.elements:
            formatClassFile.sequences.pop(formatClassFile.sequences.index(sequence))
            
        if not formatClassFile.sequences:
          returnInfo = (False,"No sequences found in input file.")
        else:
          returnInfo = (True,None)
        
      elif informationType == 'coordinates':
        
        formatClass.setIOkeywords('readCoordinates', defaultIODict, verbose = False)

        formatClass.setCoordinateFileClass()
        
        # Here have to hack to be able to handle multiple files at once...
        hasCoordinates = False
        formatClassFile = []
        for singleFileName in fileName:
          formatClass.fileName = singleFileName    
          formatClass.getCoordinates(noSequenceRead=True)
          
          if formatClass.coordinateFile and formatClass.coordinateFile.modelCoordinates:
            for modelNum in formatClass.coordinateFile.modelCoordinates.keys():
              if formatClass.coordinateFile.modelCoordinates[modelNum]:
                hasCoordinates = True
                formatClassFile.append(formatClass.coordinateFile)
                break
        
        if not hasCoordinates:
          returnInfo = (False,"Input not valid - no coordinate file or coordinates found.")
        else:
          returnInfo = (True,None)
        
      elif informationType in ('shifts',):
        
        if informationType == 'shifts':
          measurementType = 'shift'
          formatClassListName = 'chemShifts'
        
        IODict = {'measurementType': measurementType}
        IODict.update(defaultIODict)
        
        formatClass.setIOkeywords('readMeasurements', IODict, verbose = False)
        formatClass.setMeasurementCodes('read', {}, verbose = False)

        formatClass.setMeasurementFileClass()
        errorMessage = formatClass.getMeasurements()
        
        if errorMessage:
          returnInfo = (False,getFullErrorString(errorMessage))
          
        else:
          formatClassFile = formatClass.measurementFile     
          
          if not getattr(formatClassFile,formatClassListName):
            returnInfo = (False,"Input not valid - no measurement file or measurements found.")
          else:
            returnInfo = (True,None)
        
      elif informationType in ('distanceConstraints','hBondConstraints','dihedralConstraints','rdcConstraints','jCouplingConstraints'):
        
        constraintType = informationType.replace('Constraints','')
        
        IODict = {'constraintType': constraintType}
        IODict.update(defaultIODict)
        
        formatClass.setIOkeywords('readConstraints', IODict, verbose = False)
        formatClass.setConstraintCodes('read', {}, verbose = False)
        formatClass.setConstraintFileClass()
        
        errorMessage = formatClass.getConstraints()
        
        if errorMessage:
          returnInfo = (False,getFullErrorString(errorMessage))
          
        else:
          formatClassFile = formatClass.constraintFile     
          
          if not formatClassFile.constraints:
            returnInfo = (False,"Input not valid - no constraints file or constraints found.")
          else:
            returnInfo = (True,None)

      elif informationType == 'project':
        
        #
        # Note pre-parsing is not really possible on the project level, as it often depends on many
        # files, and the handling is definitely format-specific.
        #
        
        formatClass.setIOkeywords('readProject', defaultIODict, verbose = False)
        fileInformation = None
        
        
        if formatName == 'sparky':

          formatClassFile = formatClass.projectIO.SparkyProjectFile(fileName)
          formatClassFile.read(readSaveFiles = False)
          
          if hasattr(formatClass,'sequenceFile') and formatClass.sequenceFile and formatClassFile.sequenceFile.sequences and formatClassFile.saveFiles:
          
            fileInformation = {'sequenceFile': formatClassFile.sequenceFile,
                               'saveFilesList':formatClassFile.saveFiles}
               
            returnInfo = (True,fileInformation)
          
          else:
            
            returnInfo = (False,'Not a valid Sparky project file.')
        
        else:
          returnInfo = (False,'Unable to pre-parse %s project files.' % formatName)

      #
      # Keep track of the file so can re-use it later on (no need to reparse if OK)
      #
      
      if returnInfo[0]:

        if not self.preparsedFiles.has_key(informationType):
          self.preparsedFiles[informationType] = {}
        self.preparsedFiles[informationType][fileName] = formatClassFile

    except:

      returnInfo = (False,FormatConversionError("Problem pre-parsing %s file '%s' for format %s." % (informationType,fileName,formatName)).__str__())

    return returnInfo
 
  def getPeakDimInfo(self,formatClass):

    peakFile = formatClass.peakFile
    
    #
    # Check whether multiple spectra in one peaklist... TODO also need to make
    # use of keyword different from refExperimentName to set ref experiments if multiple spectra!
    #
    
    checkPeakSpectrum = False
    if len(formatClass.peakFile.numDims) > 1:
      checkPeakSpectrum = True
    
    #
    # Loop over the spectra in peaklist (usually only one)
    #
        
    spectraDict = {}
    
    for specIndex in range(len(formatClass.peakFile.numDims)):

      numDim = formatClass.peakFile.numDims[specIndex]
      if not numDim:
        continue
      
      specName = formatClass.peakFile.specNames[specIndex]
      
      peakPpmRanges = []
      
      #
      # Get information about the ppm range for each dimension... 
      #

      for dim in range(numDim):
        peakPpmRanges.append([999,-999])

      for peak in formatClass.peakFile.peaks:
      
        if checkPeakSpectrum:
          if peak.specName != specName:
            continue
      
        for dim in range(numDim):
          if peak.ppm[dim] < peakPpmRanges[dim][0]:
            peakPpmRanges[dim][0] = peak.ppm[dim]
          elif peak.ppm[dim] > peakPpmRanges[dim][1]:
            peakPpmRanges[dim][1] = peak.ppm[dim]

      if hasattr(formatClass.peakFile,'dimCodes'):
        dimCodes = formatClass.peakFile.dimCodes
      else:
        dimCodes = ['unknown'] * len(peakPpmRanges)
 
      spectraDict[specName] = {'peakPpmRanges': peakPpmRanges, 'dimensionCodes': dimCodes}
   
    return spectraDict
   
  def getRefExperimentInfo(self,refExperiment,ppmRangeDict=None):

    #
    # In order to match up, need to know dimensions... have to access reference experiment
    #
    # TODO: can refine further by adding in atomSite.name (if available!)
    #
    
    if not ppmRangeDict:
      ppmRangeDict = {}
      
    for (nucleus,shiftRanges) in (
                                  
      ('1H',  {'all': (-10,15), '15N': (5,13), '13C': (-10,7.0)}),
      ('15N', {'all': (100,140)}),
      ('13C', {'all': (5,200), '1H': (0,140), 'CA': (50,70)})
      
      ):
    
      if nucleus not in ppmRangeDict.keys():
        ppmRangeDict[nucleus] = shiftRanges
        
    refExpDimPpmRange = {}
    refExpDimInfo = {}  # Atomtype, acqu dimension?
    
    numDims = len(refExperiment.refExpDims)
    
    for refExpDim in refExperiment.sortedRefExpDims():
      refExpDimPpmRange[refExpDim.dim] = [999,-999]
      
      # TODO: Could in principle also get min/max shift range if available, but leave out for now
      atomIds = []
      
      for refExpDimRef in refExpDim.refExpDimRefs:
        for atomSite in refExpDimRef.expMeasurement.atomSites:
        
          if ppmRangeDict.has_key(atomSite.isotopeCode):
            ppmRangeInfo = ppmRangeDict[atomSite.isotopeCode]
            
            # Also set information per dimension for display
            if atomSite.name:
              atomId = atomSite.name
            else:
              atomId = atomSite.isotopeCode
            
            atomIds.append(atomId)
          
          else:
            continue
            
          #print "    ", atomSite.isotopeCode, atomSite.name, atomSite.minShift, atomSite.maxShift
          
          #
          # Set minShift and maxShift, if given - continue if so (is prime reference)
          #
          
          if atomSite.minShift:
            refExpDimPpmRange[refExpDim.dim][0] = atomSite.minShift
          
          if atomSite.maxShift:
            refExpDimPpmRange[refExpDim.dim][1] = atomSite.maxShift
            
          if atomSite.minShift and atomSite.maxShift:
            continue
          
          #
          # Set by name of atom, if possible - continue if so
          #
          
          if atomSite.name and ppmRangeInfo.has_key(atomSite.name):
            self.changePpmRange(refExpDimPpmRange[refExpDim.dim],ppmRangeInfo[atomSite.name])
            continue
          
          else:
            
            #
            # Otherwise make an educated guess based on dictionary.
            #
            # TODO: in principle these values could be part of the reference experiments,
            #       if we do an extensive chemical shift analysis, use reasonable cutoffs... 
            #       should not be too much hassle
            #
            
            setByTransfer = False
          
            for expTransfer in atomSite.expTransfers:

              if expTransfer.transferType == 'onebond':

                #
                # Use expTransfer to determine whether linked to N, C, ... . Might narrow ppm window.
                #

                expTransferAtomSites = expTransfer.sortedAtomSites()
                otherAtomSite = expTransferAtomSites[not expTransferAtomSites.index(atomSite)]

                #print "       linked to: %s %s via %s" % (otherAtomSite.isotopeCode,otherAtomSite.name,expTransfer.transferType)
                
                if ppmRangeInfo.has_key(otherAtomSite.isotopeCode):
                  setByTransfer = True
                  self.changePpmRange(refExpDimPpmRange[refExpDim.dim],ppmRangeInfo[otherAtomSite.isotopeCode])

            if not setByTransfer:            
              self.changePpmRange(refExpDimPpmRange[refExpDim.dim],ppmRangeInfo['all'])
      
      #
      # Finally set information for each dimension.
      #
      
      if not refExperiment.isReversed and refExpDim.dim == 1:
        acqString = " (acquisition dimension)"
      elif refExperiment.isReversed and refExpDim.dim == numDims:
        acqString = " (acquisition dimension)"
      else:
        acqString = ""
        
      atomIds.sort()      
      dimensionInfoString = "%s%s" % (",".join(atomIds),acqString)

      refExpDimInfo[refExpDim.dim] = dimensionInfoString
    
    return {'ppmRange': refExpDimPpmRange, 'dimInfo': refExpDimInfo}
    
  def changePpmRange(self,ppmRange,checkPpmRange):
    
    if ppmRange[0] > checkPpmRange[0]:
      ppmRange[0] = checkPpmRange[0]
    
    if ppmRange[1] < checkPpmRange[1]:
      ppmRange[1] = checkPpmRange[1]
      
  def matchFormatAndCcpnPeakDims(self,peakPpmRanges,dimCodes,refExpDimPpmRange):
  
    """
    
    Returns tuple of two dictionaries. First one is 'extended' match between the expDims and all
    matching peakDims (by index), the second one is only populated in case of a unique match and
    gives the direct match between expDim and peakDim (by index).

    """
    
    #
    # TODO: also use dimCodes, if possible. Could give better info, could also give worse info, though...
    #
    
    #
    # Determine matches
    #
    
    matches = {}
    matchesCopy = {}
    
    for expDim in refExpDimPpmRange.keys():
      matches[expDim] = []
      matchesCopy[expDim] = []
      expDimPpmRange = refExpDimPpmRange[expDim]
      
      for peakDim in range(len(peakPpmRanges)):
        peakDimPpmRange = peakPpmRanges[peakDim]
        
        if peakDimPpmRange[0] > expDimPpmRange[0] and peakDimPpmRange[1] < expDimPpmRange[1]:
          matches[expDim].append(peakDim)
          matchesCopy[expDim].append(peakDim)
    
    #
    # Now see if can get unique match...
    #
    
    uniqueMatches = {}
    expDims = matches.keys()
    
    expDimIndex = 0
    dataChanged = False
    
    while(1):
      
      expDim = expDims[expDimIndex]
      
      if len(matchesCopy[expDim]) == 1:
        peakDim = matchesCopy[expDim][0]
        if not peakDim in matchesCopy.values():
          uniqueMatches[expDim] = peakDim
          expDims.pop(expDimIndex)
          del(matchesCopy[expDim])
          for tmpExpDim in matchesCopy.keys():
            if peakDim in matchesCopy[tmpExpDim]:
              matchesCopy[tmpExpDim].pop(matchesCopy[tmpExpDim].index(peakDim))
            
          dataChanged = True
        else:
          # TODO what should happen here? 
          print "ERROR - two unique matches!"
          uniqueMatches = {}
          break
      
      else:
        expDimIndex += 1
      
      # Start loop again, only exit if all done (no more changes or all assigned)
      if not expDims:
        break
      elif expDimIndex >= len(expDims):
        #print expDimIndex, len(expDims), dataChanged
        if dataChanged:
          dataChanged = False
          expDimIndex = 0
        else:
          break
    
    #
    # If unique matches 1 less than total dims, set the final one!
    #
    if len(uniqueMatches) == len(refExpDimPpmRange) - 1:
      missingDim = None
      missingExpDim = None
      for expDimIndex in range(len(refExpDimPpmRange)):
        print expDimIndex,
        if not expDimIndex in uniqueMatches.values():
          missingExpDim = expDimIndex
        if not expDimIndex + 1 in uniqueMatches.keys():
          missingDim = expDimIndex + 1

      if missingExpDim != None and missingDim != None:
        uniqueMatches[missingDim] = missingExpDim
      
    #print matches
    #print uniqueMatches
    
    return (matches,uniqueMatches)

  #
  # Lower level functions to do import/export
  #
  
  def executeIO(self,funcType,fileName,informationType,formatName):
        
    # Exceptions...
    if formatName == 'cyana' and funcType == 'export' and informationType in ('shifts',):
      if informationType == 'peaks':
        self.addKeywords['cyanaFormat'] = True
      elif informationType == 'shifts':
        if funcType == 'export':
          self.addKeywords['cyana21Naming'] = True
      self.exportFile(informationType, 'xeasy', fileName, addKeywords = self.addKeywords)
  
    else:

      self.funcType = funcType
      self.fileName = fileName
      self.informationType = informationType
      self.formatName = formatName
      
      formatFunc = self.getFormatFunction()

      keywords = self.getKeywords()
    
      self.conversionReturn = formatFunc(**keywords)
      
      (self.conversionInfo,self.conversionSuccess) = self.getConversionOutcome()
      
  def getConversionOutcome(self):
          
    formatObject = self.getFormatClass(self.formatName)

    conversionSuccess = True
    if formatObject.convertCount.keys():
      conversionSuccess = False
      for convertMainObj in formatObject.convertCount.keys():
        conversionInfo = formatObject.convertCount[convertMainObj]
        if conversionInfo[1].keys():
          for convertSubObj in conversionInfo[1].keys():
            if conversionInfo[1][convertSubObj]:
              conversionSuccess = True
              break
        elif conversionInfo[0]:
          conversionSuccess = True
          
        if conversionSuccess:
          break
      
    convertInfoString = formatObject.getConvertCountString(startIndent = "")
    
    return (convertInfoString,conversionSuccess)

  def getFormatFunction(self):
    
    formatClass = self.getFormatClass(self.formatName)
        
    # NMR-STAR export exception!
    
    if self.formatName == 'nmrStar' and self.funcType == 'export':
      
      formatFunc = self.exportNmrStar
      
    # This is for all FormatConverter based stuff...    
    elif self.informationTypeDict.has_key(self.informationType):
    
      if self.informationTypeDict[self.informationType].has_key(self.formatName):
        
        if self.funcType == 'import':
          index = 0              
          
        elif self.funcType == 'export':
          index = 1
        
        if self.informationTypeDict[self.informationType][self.formatName][index]:
          funcName = self.informationTypeDict[self.informationType][self.formatName][index]
          formatFunc = getattr(formatClass,funcName)
         
        else:
          raise FormatConversionError("No '%s' function available for '%s' information type, format '%s'!" % (self.funcType,self.informationType,self.formatName))
      
      else:
        raise FormatConversionError("No '%s' information type available for format '%s'!" % (self.informationType,self.formatName))
    
    else:
      raise FormatConversionError("No '%s' information type available!" % self.informationType)
    
    return formatFunc    

  def determineFormatNamesForFile(self,informationType,filePath,formatNameList = None):
  
    """
    This is a more generic function that determines valid formats for a given file and informationType
    """
  
    formatNameSuggestions = []
    progress = None
    testMode = False
    
    # Fix for handling coordinates
    #if informationType == 'coordinates' and type(filePath) == type(''):
    #  filePath = (filePath,)

    # Silent mode.. hide any prints or errors that are raised by the import testing
    if self.silentRead and not testMode:
      self.goSilent()
    
    try:
    
      if not formatNameList:
        formatNameList = self.informationTypeDict[informationType]
        
      if self.guiRoot:
        progress = self.ProgressBar(self.guiRoot,text="Determining %s format..." % informationType,total=len(formatNameList))

      for formatName in formatNameList:

        if testMode:
          self.goSilent()
          
        (fileRead,fileInformation) = self.preparseFile(informationType,formatName,filePath)
        
        if testMode:
          self.endSilent()   
          print informationType,formatName, fileRead, fileInformation
        
        if fileRead:
          formatNameSuggestions.append(formatName)      
          
        if self.guiRoot:
          progress.increment()
          
    except:
      if testMode:
        raise

    finally:
      # Return to normal print and error reporting
      self.endSilent()
      
      if self.guiRoot and progress:
        progress.destroy()
  
    return formatNameSuggestions
    
  def goSilent(self):
    
    # for stdout, stderr management
    self.stdstdout = sys.stdout
    self.stdstderr = sys.stderr
  
    sys.stderr = open('/dev/null', 'w')
    sys.stdout = open('/dev/null', 'w')

  def endSilent(self):

    #sys.stderr = sys.__stderr__
    #sys.stdout = sys.__stdout__

    sys.stderr = self.stdstderr
    sys.stdout = self.stdstdout
    
  def getKeywords(self):
  
    keywords = {}
    
    #
    # Set preparsed file, if possible
    #    
    
    formatClass = self.getFormatClass(self.formatName)

    if self.preparsedFiles.has_key(self.informationType):
      if self.preparsedFiles[self.informationType].has_key(self.fileName):
        fileAttrName = fileAttrNames[self.informationType]
        setattr(formatClass,fileAttrName,self.preparsedFiles[self.informationType][self.fileName])
        formatClass.file = getattr(formatClass,fileAttrName)
        
        # Don't reset fileName when reading list of filenames for coordinates!
        if self.informationType != 'coordinates' or type(self.fileName) == type(''):
          self.fileName = None  
        
        print "  Warning: using preparsed file to read %s information..." % self.informationType
    
    #
    # Set single filename, or multiple ones if not string (assuming tuple or list)
    #
    
    if self.fileName:

      if type(self.fileName) == type(''):
        keywords['fileName'] = self.fileName
      else:
        keywords['fileNames'] = self.fileName
    
    #
    # Check for required keywords
    #
    
    self.checkRequiredKeywords(keywords)
    
    #
    # Set defaults depending on import/export, format and information type
    #  
    
    self.setDefaultKeywords(keywords)
      
    #
    # Set added keywords
    #
    
    for keyword in self.addKeywords.keys():
      if keywords.has_key(keyword):
        print "Warning: overwriting default keyword '%s' by set value!" % keyword
      keywords[keyword] = self.addKeywords[keyword]
      
    #
    # Hack to set whether molecular system created
    #

    if hasattr(self,'molSystemCreated') and self.molSystemCreated:
      formatClass.molSystemCreated = True
      self.molSystemCreated = False
      
    return keywords
    
  def checkRequiredKeywords(self,keywords):
  
    if self.funcType == 'import':
      if self.informationType == 'peaks':
        
        #
        # TODO: parts of this code have to be moved elsewhere for 'pre-parsing' information!!!
        #
        
        referenceInfoFound = False
        
        if 'dataSource' in self.addKeywords:
          
          keywords['dataSource'] = self.addKeywords['dataSource']
          referenceInfoFound = True
        
        elif 'dataSourceId' in self.addKeywords:
          
          dataSourceId = self.addKeywords.pop('dataSourceId')
          
          (selectionList,selectionDict) = self.getDataSourceInfo()
          
          if selectionDict.has_key(dataSourceId):
            keywords['dataSource'] = selectionDict[dataSourceId]
            referenceInfoFound = True
          else:
            raise FormatConversionError("Data source Id '%s' not recognised. Pick one from this list:\n%s" % (dataSourceId,'\n'.join(selectionList)))
            
        elif 'refExperiment' in self.addKeywords or 'refExperimentIdString' in self.addKeywords:          
          
          refExperiment = None
          if 'refExperimentIdString' in self.addKeywords:
            refExpIdString = self.addKeywords.pop('refExperimentIdString')

            for numDim in range(1,5):
              self.getNmrRefExperimentNames(numDim)
              if self.allRefExpIdStrings.has_key(refExpIdString):
                refExperiment = self.allRefExpIdStrings[refExpIdString]
                break
              
          else:
            refExperiment = self.addKeywords.pop('refExperiment')
          # TODO Might want to use other name than refExperiment.name
          
          if refExperiment:
            experiment = createExperiment(self.ccpnProject,refExperiment.name,refExperiment)
    
            numDim = len(experiment.expDims)
    
            keywords['dataSource'] = createPpmFreqDataSource(experiment,refExperiment.name,'processed',numDim)    
            
            referenceInfoFound = True

        if not referenceInfoFound:
          
          (selectionList,selectionDict) = self.getDataSourceInfo()
          
          expNames = []
          for numDim in range(1,5):
            self.getNmrRefExperimentNames(numDim)
            for expName in self.allRefExpIdStrings.keys():
              if not expName in expNames:
                expNames.append(expName)
          expNames.sort()
          
          raise FormatConversionError("\n\n*** Need either 'dataSourceId' or 'refExperiment' or 'refExperimentIdString' keyword ***\nFor 'dataSourceId', pick one from this list:\n%s\n\nFor 'refExperimentIdString', pick one from this list:\n%s" % ('\n'.join(selectionList),'\n'.join(expNames)))

  def setDefaultKeywords(self,keywords):
  
    # TODO: set some info in here! Also: how to handle object selection for export? Set latest import objects, see if those available
    # If not, then have to select somehow - again try to be smart about what the user wants here.

    keywords['minimalPrompts'] = True

    if self.funcType == 'import':
    
      if self.informationType == 'coordinates':
        
        keywords['strucGen'] = self.getStructureGeneration()
        keywords['linkAtoms'] = False
        
      elif self.informationType == 'sequence':
      
        keywords['molSystem'] = self.getMolSystem()
        
      elif self.informationType == 'shifts':
        
        pass
        
      elif self.informationType == 'peaks':
        
        #
        # DataSource already set higher up - TODO should I do this here instead?
        #
            
        keywords['dataDimRefs'] = self.getDataDimRefs(keywords)
        
      elif self.informationType in ('distanceConstraints','hBondConstraints','dihedralConstraints','rdcConstraints','jCouplingConstraints'):
        
        pass
    
    elif self.funcType == 'export':
    
      if self.informationType == 'coordinates':
      
        keywords['structures'] = self.getModels()
        keywords['forceDefaultChainMapping'] = True
        
      elif self.informationType == 'sequence':
        
        keywords['chains'] = self.getChains()
        keywords['forceDefaultChainMapping'] = True
        
      elif self.informationType == 'shifts':
        
        if not self.resonancesLinked:
          self.linkResonances()
        
        keywords['measurementList'] = self.getMeasurementList(self.informationType)
        keywords['forceDefaultChainMapping'] = True
        
      elif self.informationType == 'peaks':
      
        if not self.resonancesLinked:
          self.linkResonances()
        
        # Function will pop info from self.addKeywords!
        keywords['peakLists'] = self.getPeakLists()
        keywords['dataDimRefs'] = self.getDataDimRefsForExport(keywords['peakLists'])
        keywords['forceDefaultChainMapping'] = True
        
      elif self.informationType in ('distanceConstraints','hBondConstraints','dihedralConstraints','rdcConstraints','jCouplingConstraints'):
      
        if not self.resonancesLinked:
          self.linkResonances()
          
        keywords['constraintList'] = self.getConstraintList(self.informationType)
        keywords['forceDefaultChainMapping'] = True
        
      elif self.informationType in ('project',):
        
        if not self.resonancesLinked:
          self.linkResonances()
        
        # The keywords are specific depending on the format - handled in this function
        self.setProjectExportKeywords(keywords)

  
  def getStructureGeneration(self,entry=None):
    
    if not entry:
      entry = self.nmrEntry
    
    # TODO this has to be more sophisticated, by working with iterations and selecting the 'right' strucGen?
    
    structureGeneration = entry.findFirstStructureGeneration()
    if not structureGeneration:
      if len(self.nmrProject.structureGenerations) == 1:
        # TODO different behaviour when reading or writing here?
        structureGeneration = self.nmrProject.findFirstStructureGeneration()
      else:
        structureGeneration = self.nmrProject.newStructureGeneration()
      entry.addStructureGeneration(structureGeneration)
      
    return structureGeneration

  def getMolSystem(self,entry=None):
    
    if not entry:
      entry = self.nmrEntry
  
    molSystem = entry.molSystem
    
    if not molSystem:
      
      if len(self.ccpnProject.molSystems) == 1:
        molSystem = self.ccpnProject.findFirstMolSystem()
        entry.molSystem = molSystem
        
      if not molSystem:
       
        for i in range(10):
          molSystemCode = 'Molecular_system_%d' % i
          molSystem = self.ccpnProject.findFirstMolSystem(code = molSystemCode)
          if not molSystem:
            molSystem = self.ccpnProject.newMolSystem(code = molSystemCode)
            entry.molSystem = molSystem
            self.molSystemCreated = True
            break

    return molSystem
    
    
  def getChains(self):

    if 'chains' in self.addKeywords.keys():
      chains = self.addKeywords.pop('chains')
    
    else:
    
      molSystem = self.getMolSystem()
      
      chains = molSystem.sortedChains()
    
    return chains
    
  def getModels(self):
  
    if 'structures' in self.addKeywords.keys():
      models = self.addKeywords.pop('structures')
  
    else:
      strucGen = self.getStructureGeneration()
      structureEnsemble = strucGen.structureEnsemble
      
      if not structureEnsemble:
        structureEnsemble = self.ccpnProject.findFirstStructureEnsemble()
      
      if structureEnsemble:
        models = structureEnsemble.sortedModels()
      else:
        raise FormatConversionError("No structure ensemble available!")
      
    return models
    
  def getMeasurementList(self,informationType):
    
    if 'measurementList' in self.addKeywords and self.addKeywords['measurementList']:
      # TODO this might be redundant...
      measurementList = self.addKeywords.pop('measurementList')
            
    else:
      
      if informationType == 'shifts':
        listType = 'ShiftList'
      else:
        # TODO add other ones!
        pass
        
      
      listKeyword = listType[0].lower() + listType[1:]
      if listKeyword in self.addKeywords and self.addKeywords[listKeyword]:
        measurementList = self.addKeywords.pop(listKeyword)
        
      else:

        measurementLists = self.nmrEntry.findAllMeasurementLists(className = listType)

        if not measurementLists:       
          measurementLists = self.nmrProject.findAllMeasurementLists(className = listType)

        if not measurementLists:
          raise FormatConversionError("No %s present - cannot export." % listType)
        elif len(measurementLists) > 1:  
          # TODO can here have option to just select first one?
          # TODO Or, alternatively, export all relevant lists?
          raise FormatConversionError("Multiple %ss present - have to select one and pass it as 'measurementList' using addKeywords" % listType)
        else:
          measurementList = measurementLists.pop()

    return measurementList

  def getConstraintList(self,informationType):
    
    if 'constraintList' in self.addKeywords:
      # TODO this might be redundant...
      constraintList = self.addKeywords.pop('constraintList')
      
    else:
      
      # TODO: should probably get this from some central dictionary! Also useful in DataFormat.py?
      if informationType == 'distanceConstraints':
        listType = 'DistanceConstraintList'
      elif informationType == 'hBondConstraints':
        listType = 'HBondConstraintList'
      elif informationType == 'dihedralConstraints':
        listType = 'DihedralConstraintList'
      elif informationType == 'rdcConstraints':
        listType = 'RdcConstraintList'
      elif informationType == 'jCouplingConstraints':
        listType = 'JCouplingConstraintList'
      else:
        # TODO add other ones!
        pass
      
      constraintLists = set()
      
      for nmrConstraintStore in self.ccpnProject.sortedNmrConstraintStores():
        constraintLists.update(nmrConstraintStore.findAllConstraintLists(className = listType))
      
      if not constraintLists:
        raise FormatConversionError("No %s present - cannot export." % listType)
      elif len(constraintLists) > 1:  
        # TODO can here have option to just select first one?
        # TODO Or, alternatively, export all relevant lists?
        print constraintLists
        raise FormatConversionError("Multiple %ss present - have to select one and pass it as 'constraintList' using addKeywords" % listType)
      else:
        constraintList = constraintLists.pop()

    return constraintList
                        
  def getPeakLists(self):
    
    if 'peakLists' in self.addKeywords and self.addKeywords['peakLists']:
      # TODO this might be redundant...
      peakLists = self.addKeywords.pop('peakLists')
      
    else:
      
      peakLists = []
      
      for nmrExp in self.nmrProject.experiments:
        for dataSource in nmrExp.dataSources:
          for peakList in dataSource.peakLists:
            peakLists.append(peakList)
      
      if not peakLists:
        raise FormatConversionError("No peakLists present - cannot export.")
      elif len(peakLists) > 1:  
        # TODO can here have option to just select first one? Or one based on dataSource/NMR experiment?
        # TODO Or, alternatively, export all peakLists?
        raise FormatConversionError("Multiple peakLists present - have to select one and pass it as 'peakList' using addKeywords")

    return peakLists
  
  def getDataDimRefsForExport(self,peakLists):
    
    dataDimRefs = []
    
    if 'dataDimRefs' in self.addKeywords.keys():
      dataDimRefs = self.addKeywords['dataDimRefs']
    
    else:     
      refPeak = peakLists[0].findFirstPeak()
      
      if not refPeak:
        raise FormatConversionError("PeakList has no peaks - cannot export.")
  
      else:
        # TODO this should be customisable!
        for peakDim in refPeak.sortedPeakDims():
          # Ignore SampledDataDims!
          if peakDim.dataDim.className == 'SampledDataDim':
            continue
          dataDimRefs.append(peakDim.dataDimRef)
        
        if None in dataDimRefs:
          raise FormatConversionError("Not all peakList dimensions linked to experiment ones - cannot export.")
    
    return dataDimRefs
        
  def setProjectExportKeywords(self,keywords):
    
    #
    # Chain keyword
    #
    
    if self.formatName in projectExportRequirements['chain']:

      if 'chain' in self.addKeywords:
        chain = self.addKeywords.pop('chain')
      else:
        chains = self.getChains()
        
        if len(chains):
          if len(chains) > 1:
            print "  Warning: multiple chains available - picking first one only."
          chain = chains.pop()
        else:
          raise FormatConversionError("No chain available - cannot export project.")
      
      # Can only do one chain at a time - again best to select this...
      keywords['chain'] = chain
    
    #
    # MolSystem
    #
    
    if self.formatName in projectExportRequirements['molSystem']:
      
      molSystem = self.getMolSystem()
      
      if molSystem:
        keywords['molSystem'] = molSystem

    #
    # Shiftlist
    #

    if self.formatName in projectExportRequirements['shiftList']:
      
      keywords['shiftList'] = self.getMeasurementList('shifts')
      
    #
    # Entry
    #

    if self.formatName in projectExportRequirements['entry']:

      if 'entry' in self.addKeywords and self.addKeywords['entry']:
        entry = self.addKeywords.pop('entry')
      elif self.formatName == 'sparky' and 'molSystem' in keywords.keys():
        # Not necessary for sparky if MolSystem selected.
        entry = None
      else:
        # This is always set.
        entry = self.nmrEntry

      """
      else:
        entries = set()
        
        for nmrEntryStore in self.ccpnProject.nmrEntryStores:
          entries.update(nmrEntryStore.entries)
        
        if len(entries):
          if len(entries) > 1:
            print "  Warning: multiple entries available - picking first one only."
          entry = entries.pop()
        else:
          raise FormatConversionError("No entry available - cannot export project.")
      
      # Create an entry if none exists
      if not entry:
        nmrEntryStore = self.ccpnProject.currentNmrEntryStore
        if not nmrEntryStore:
          nmrEntryStore = self.ccpnProject.findFirstNmrEntryStore()
          if not nmrEntryStore:
            nmrEntryStore = self.ccpnProject.newNmrEntryStore(name='Format conversion %s' % self.ccpnProject.name)

        entry = nmrEntryStore.newEntry(name='Format conversion %s' % self.ccpnProject.name)
      """
      
      # Now make sure most important connections are set!
      # TODO probably have to make this generic code to create a default entry!! Chris P. has something?
      if not entry.molSystem:
        # Will set link automatically
        self.getMolSystem(entry=entry)
      
      if not entry.structureGenerations:
        # Will set link automatically
        self.getStructureGeneration(entry=entry)
      
      if not entry.authors and self.ccpnProject.currentAffiliationStore:
        for person in self.ccpnProject.currentAffiliationStore.sortedPersons():
          entry.addAuther(person)
      
      if not entry.experiments:
        for nmrExperiment in self.nmrProject.experiments:
          entry.addExperiment(nmrExperiment)
      
      if not entry.measurementLists:
        for measurementList in self.nmrProject.measurementLists:
          entry.addMeasurementList(measurementList)
      
      if not entry.peakLists:
        for nmrExperiment in self.nmrProject.experiments:
          for dataSource in nmrExperiment.dataSources:
            for peakList in dataSource.peakLists:
              entry.addPeakList(peakList)
      
      # Set to keywords
      
      keywords['entry'] = entry
      
    #
    # Directory to write program input files to (so to export into)
    #
    
    if self.formatName in ('mars',):
      
      if 'inputDirectory' in self.addKeywords:
        inputDirectory = self.addKeywords.pop('inputDirectory')
      else:
        inputDirectory = "./%sInput" % self.formatName
        
      keywords['inputDirectory'] = inputDirectory

  def exportNmrStar(self,**keywords):
    
    try:
      from pdbe.nmrStar.IO.NmrStarExport import NmrStarExport
     
      nmrStarExport = NmrStarExport(keywords['entry'], nmrStarVersion = '3.1')  
      nmrStarExport.createFile(self.fileName)
      nmrStarExport.writeFile()
      
    except:
      raise FormatConversionError("NMR-STAR export not available!")
  
  def checkUnlinkedResonancesPresent(self,fractionCutoff = 0.05):

    #
    # Check unlinked resonance percentages
    #
    
    unlinkedResonancesPresent = False
    
    #
    # First do NMR resonances
    #
    
    fractionUnlinkedResonances = self.getFractionUnlinkedResonances(self.nmrProject.resonances)

    if fractionUnlinkedResonances > fractionCutoff:
      unlinkedResonancesPresent = True
        
    #
    # Then NMR constraint stores, if necessary
    #  
    
    if not unlinkedResonancesPresent:
      
      for nmrConstraintStore in self.ccpnProject.sortedNmrConstraintStores():
        
        fractionUnlinkedResonances = self.getFractionUnlinkedResonances(nmrConstraintStore.fixedResonances)
        
        if fractionUnlinkedResonances > fractionCutoff:
          unlinkedResonancesPresent = True
          break
      
    return unlinkedResonancesPresent
      
  def getFractionUnlinkedResonances(self,resonances):

    numUnlinkedResonances = 0
    
    if resonances:
      
      numResonances = len(resonances)
    
      for resonance in resonances:
  
        resonanceSet = resonance.resonanceSet
  
        if not resonanceSet:
          
          numUnlinkedResonances += 1
    
    else:
      numResonances = 1
    
    return (numUnlinkedResonances * 1.0 / numResonances)
    
  def linkResonances(self,fractionCutoff = 0.05, **keywds):
  
    #
    # Determine whether NMR resonances or NMRConstraint ones
    #
    # TODO: could pass this in, in principle, otherwise do everything
    #
    
    if not 'minimalPrompts' in keywds.keys():
      keywds['minimalPrompts'] = True
    
    if not 'allowPopups' in keywds.keys():
      if not self.guiRoot:
        keywds['allowPopups'] = False
      else:
        keywds['allowPopups'] = True
        
    if self.guiRoot:
      keywds['guiParent'] = self.guiRoot

    if not 'useLinkResonancePopup' in keywds.keys():
      keywds['useLinkResonancePopup'] = False
      
    if not 'forceChainMappings' in keywds.keys():
      keywds['forceDefaultChainMapping'] = True

    #
    # Now run linkResonances, but only if there are a reasonable amount of unlinked resonance present
    #
    
    self.numResonancesLinked = {'linked': [], 'unlinked': [], 'info': [], 'origUnlinked': [], 'origLinked': []}

    for nmrConstraintStore in self.ccpnProject.sortedNmrConstraintStores():
      
      numResonances = len(nmrConstraintStore.fixedResonances)
      fractionUnlinkedRes = self.getFractionUnlinkedResonances(nmrConstraintStore.fixedResonances)
      
      if fractionUnlinkedRes > fractionCutoff:    
        localKeywds = self.setForceChainMappings(keywds,nmrConstraintStore.fixedResonances)
        lr = linkResonances(self.ccpnProject, nmrConstraintStore = nmrConstraintStore, **localKeywds)
        self.numResonancesLinked['linked'].append(lr.numResonancesLinked)
        self.numResonancesLinked['unlinked'].append(self.getFractionUnlinkedResonances(nmrConstraintStore.fixedResonances) * numResonances)
        self.numResonancesLinked['origUnlinked'].append(fractionUnlinkedRes * numResonances)
        self.numResonancesLinked['origLinked'].append(numResonances - (fractionUnlinkedRes * numResonances))
        self.numResonancesLinked['info'].append("NMR constraint store %d" % nmrConstraintStore.serial)
    
    numResonances = len(self.nmrProject.resonances)
    fractionUnlinkedRes = self.getFractionUnlinkedResonances(self.nmrProject.resonances)
    
    if fractionUnlinkedRes > fractionCutoff:    
      lr =linkResonances(self.ccpnProject, nmrProject = self.nmrProject, **keywds)
      self.numResonancesLinked['linked'].append(lr.numResonancesLinked)
      self.numResonancesLinked['unlinked'].append(self.getFractionUnlinkedResonances(self.nmrProject.resonances) * numResonances)
      self.numResonancesLinked['origUnlinked'].append(fractionUnlinkedRes * numResonances)
      self.numResonancesLinked['origLinked'].append(numResonances - (fractionUnlinkedRes * numResonances))
      self.numResonancesLinked['info'].append("NMR project %s" % self.nmrProject.name)
    
    if self.numResonancesLinked['linked'].count(0) != len(self.numResonancesLinked['linked']):
      self.resonancesLinked = True
    
    
    self.allResonancesLinked = True

    for i in range(len(self.numResonancesLinked['linked'])):
      if self.numResonancesLinked['origUnlinked'][i]:
        print "Linked %d resonances for %s." % (self.numResonancesLinked['linked'][i],self.numResonancesLinked['info'][i]),
        if self.numResonancesLinked['unlinked'][i]:
          print " Unable to link %d resonances." % self.numResonancesLinked['unlinked'][i]
          self.allResonancesLinked = False
        else:
          print " All resonances now linked."

  def setForceChainMappings(self,keywds,resonances):

    if 'forceChainMappings' not in keywds.keys():
      sequenceComparison = SequenceCompare()
      sequenceComparison.createCcpnChainInformation(self.ccpnProject.currentMolSystem.sortedChains())
      sequenceComparison.getFormatFileInformation(resonances,self.formatName)
  
      sequenceComparison.createFormatFileChainInformation()

      forceChainMappings = sequenceComparison.compareFormatFileToCcpnInfo()        
      if forceChainMappings:
        keywds['forceChainMappings'] = forceChainMappings
        print forceChainMappings
        
    return keywds.copy()
             
  def getValidExportFormats(self):
    
    (self.validExportFormats,self.ccpnObjects) = getValidExportFormats(self.ccpnProject)
    
    #
    # Now make strings out of objects for selection lists...    
    #
    
    ccpnObjectSelectionStrings = {}
    
    for ccpnClass in self.ccpnObjects.keys():
      (selectionList,selectionDict) = createSelection(self.ccpnObjects[ccpnClass])
      
      selectionList.sort()

      self.updateCcpnObjectDict(ccpnClass,selectionDict)
      ccpnObjectSelectionStrings[ccpnClass] = selectionList
      
    return (self.validExportFormats,ccpnObjectSelectionStrings)
  
  def getCcpnObjectSelectionStrings(self,ccpnObjects):
    
    ccpnObjectSelectionStrings = []
    
    if ccpnObjects:
      ccpnClass = ccpnObjects[0].className
      
      (selectionList,selectionDict) = createSelection(ccpnObjects)
      
      # Update, not replace...
      self.updateCcpnObjectDict(ccpnClass,selectionDict)
      ccpnObjectSelectionStrings = selectionList
      
    return ccpnObjectSelectionStrings

  def updateCcpnObjectDict(self,ccpnClass,selectionDict):
    
    if not self.ccpnObjectDict.has_key(ccpnClass):
      self.ccpnObjectDict[ccpnClass] = {}
    
    # OK to overwrite existing values...
    self.ccpnObjectDict[ccpnClass].update(selectionDict)
  
  def getCcpnObjectFromSelectionString(self,ccpnClass,selectionString):
    
    ccpnObject = None
    
    if ccpnClass:
      if self.ccpnObjectDict.has_key(ccpnClass):
        if self.ccpnObjectDict[ccpnClass].has_key(selectionString):
          ccpnObject = self.ccpnObjectDict[ccpnClass][selectionString]
        else:
          raise FormatConversionError("No object for selection string '%s', CCPN class '%s'." % (selectionString,ccpnClass))
      else:
        raise FormatConversionError("No CCPN class '%s' objects." % (ccpnClass))
      
    return ccpnObject
  
  def convertObjectStrings(self,objectStrings):
    
    if objectStrings:
      for addKeyword in objectStrings.keys():
        objectInfo = objectStrings[addKeyword]      
        if type(objectInfo) == type([]):
          self.addKeywords[addKeyword] = []
          for (ccpnClass,selectionString) in objectInfo:
            object = self.getCcpnObjectFromSelectionString(ccpnClass,selectionString)
            if object:
              self.addKeywords[addKeyword].append(object)
        else:
          object = self.getCcpnObjectFromSelectionString(objectInfo[0],objectInfo[1])
          self.addKeywords[addKeyword] = object
    
  def getDataSourceInfo(self):

    dataSources = []
    for exp in self.nmrProject.experiments:
      for ds in exp.dataSources:
        dataSources.append(ds)

    return createSelection(dataSources)
    
  def getDataDimRefs(self,keywords):
  
    #
    # Get from input, if given, if not try to set automatically.
    #
    
    dataSource = keywords['dataSource']
    
    dataDimRefs = []
    
    if self.addKeywords.has_key('expDimToPeakDim') and self.addKeywords['expDimToPeakDim']:
      uniqueMatches = self.addKeywords['expDimToPeakDim']
      for i in range(0,max(uniqueMatches.values()) + 1):
        for expDim in uniqueMatches.keys():
          if uniqueMatches[expDim] == i:
            dataDim = dataSource.findFirstDataDim(dim = expDim)
            dataDimRefs.append(dataDim.findFirstDataDimRef())  # TODO: Guess need to qualify this by atom type...    
            #print i, expDim, dataDimRefs[-1]
      
      del(self.addKeywords['expDimToPeakDim'])
  
    return dataDimRefs

  #
  # NMR reference experiment handling...
  #

  def setNmrRefExperimentInfo(self):
  
    self.nmrRefExperimentDict = {}
    self.allRefExpIdStrings = {}  

    for nep in self.ccpnProject.sortedNmrExpPrototypes():

      nmrExpPrototypeName = nep.name

      for re in nep.sortedRefExperiments():
        
        #
        # Get atom information
        #
        
        atomDimList = []
        for refExpDim in re.refExpDims:
          atomsPerDimList = []
          for refExpDimRef in refExpDim.refExpDimRefs:
            expMeasurement = refExpDimRef.expMeasurement
            atoms = []
            addText = ""
            for atomSite in expMeasurement.atomSites:
              atoms.append(atomSite.isotopeCode[-1])
            if expMeasurement.measurementType != "Shift":
              addText = "(%s)" % expMeasurement.measurementType
            atomsPerDimList.append("%s%s" % (''.join(atoms),addText))
          atomDimList.append(','.join(atomsPerDimList))
        
        #
        # Get number of dimensions
        #
        
        numDim = len(re.refExpDims)
        
        #
        # Get a useable name...
        #
        
        systematicName = re.findFirstSystematicName(namingSystem = 'BMRB')
        
        if systematicName:
          bmrbName = systematicName.name
        else:
          bmrbName = None
          
        #
        # Set the reference information
        #
        
        self.nmrRefExperimentDict[re] = (atomDimList,numDim,bmrbName,nep.synonym,nep.name,re.name)

  def getNmrRefExperimentNames(self,numDim,onlyCommon = True):
    
    nmrRefExperiments = {}
    
    for refExp in self.nmrRefExperimentDict.keys():
    
      (atomDimList,expNumDim,bmrbName,synonym,protocolName,refExpName) = self.nmrRefExperimentDict[refExp]
      
      if numDim == expNumDim and (not onlyCommon or (onlyCommon and bmrbName)):
        
        if bmrbName:
          useName = bmrbName
        elif synonym:
          useName = synonym
        else:
          useName = refExpName
        
        idString = "%s (%s)" % (useName,'->'.join(atomDimList))
        nmrRefExperiments[idString] = refExp
        
        if not self.allRefExpIdStrings.has_key(idString):
          self.allRefExpIdStrings[idString] = refExp

    return nmrRefExperiments

  def runFromXml(self,xmlFileName,allowedFormats=None):
  
    from memops.universal.ElementTree import ElementTree as ET
   
    root = ET.ElementTree(file=xmlFileName)
    
    topElement = root.getroot()
    
    assert topElement.tag == 'conversionProject'
    assert topElement.get("status") == 'unconverted'
    
    # First handle input data...
    # Should probably have some input checks here
    for inputElement in topElement.findall("input"):
      format = inputElement.get("format")
      
      # This can be extended, of course
      assert not allowedFormats or format in allowedFormats
      
      for fileElement in inputElement.findall("file"):
        fileType = fileElement.get("type")
        fileName = fileElement.text
        
        self.importFile(fileType,format,fileName)
    
    
    #
    # Now decide what to do on CCPN level
    #
    
    ccpnElement = topElement.find("ccpnProject")
    
    if ccpnElement.get("linkResonances") == 'yes':
      self.linkResonances(setSinglePossEquiv = True, setSingleProchiral = True)
    
    # TODO could also add 'packCcpnProject' here! Might make more sense...
    if ccpnElement.get("saveData") == 'yes':
      self.saveCcpnProject()
      
    #
    # Now export
    #
  
    for outputElement in topElement.findall("output"):
      format = outputElement.get("format")
      
      # This can be extended, of course
      assert not allowedFormats or format in allowedFormats
      
      for fileElement in outputElement.findall("file"):
        fileType = fileElement.get("type")
        fileName = fileElement.text
        
        self.exportFile(fileType,format,fileName)

    #
    # Finished OK, reset status tag and re-save (?)
    #
    
    topElement.set("status",'finished')
    
    root.write(xmlFileName + ".test")

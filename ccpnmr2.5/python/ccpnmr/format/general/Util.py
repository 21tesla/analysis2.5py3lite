
"""
======================COPYRIGHT/LICENSE START==========================

Util.py: Useful functions for scripts in ccpnmr/format and subdirectories

Copyright (C) 2005-2010 Wim Vranken (European Bioinformatics Institute)

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

from memops.universal.Util import returnInt, returnFloat, returnLong

from memops.general.Util import returnMemopsLine

from memops.api import Implementation
from ccp.api.general import Method, Instrument, Citation
from ccp.api.general import Affiliation
from ccp.api.molecule import Molecule, MolSystem, ChemComp, MolStructure
from ccp.api.nmr import Nmr, NmrConstraint, NmrEntry


from ccpnmr.format.general.Constants import tagSep, oldToNewExpDict
from ccpnmr.format.general.Constants import commonNmrActiveIsotopes
from ccpnmr.format.general.Constants import softwareDefinitions, projectExportRequirements
from ccpnmr.format.general.Constants import INDIRECT_FREQ_CONV
from ccpnmr.format.general.Constants import allFormatsDict, IoSetupList
from ccpnmr.format.general.Constants import ccpNmr_kw, assign_kw, allResidueAtoms_kw

from ccp.format.general.Constants import defaultSeqInsertCode
from ccp.format.general.Util import getSeqAndInsertCode
from ccp.general.Util import setUniqueAppData, setCurrentStore

#
# Simple error class
#

class SimpleError(StandardError):

  def __init__(self, value):
    
    self.value = value
        
  def __str__(self):
    
    if self.value.count("\n"):
      returnStr = self.value
    else:
      returnStr = repr(self.value)
    
    return returnStr  

#
# Fc info dictionaries
#

def setNonDictKeywords(function,IOkeywordList,numMandatory):

  funcName = function.func_name
  
  if funcName.count('write'):
    actionType = 'write to'
  else:
    actionType = 'read in'

  for (IOkeyword,default,mandatory,infoText) in (('fileName',None,True,'Always first keyword, is name of file to %s. The file name includes the directory path.' % actionType),
                                                 ('fileNames',None,True,'Always first keyword, is list of names of file to %s. The file names include the directory path.' % actionType),
                                                 ('verbose',False,False,'Print out verbose information')):

    if IOkeyword in function.func_code.co_varnames:
      if mandatory:

        IOkeywordList.insert(numMandatory,(IOkeyword,default,mandatory,False,infoText))
        numMandatory += 1

      else:

        IOkeywordList.append((IOkeyword,default,mandatory,False,infoText))

def getWidgetInfo(functionWidgetInfo,widgetInfo):
  
  for (widgetType,fcKeyword,label,defaults,void) in widgetInfo:
  
    functionWidgetInfo[fcKeyword] = (widgetType,defaults,label)

def getFormatConverterInfoDict():

  automaticKeywords = {'Constraints': 'constraintType',
                       'Measurements': 'measurementType'}

  ccpnProject = Implementation.MemopsRoot(name = 'test')
  
  formatNames = allFormatsDict.keys()
  formatNames.sort()
  
  for removeFormatName in ('generic','marvin','regine'):
    formatNames.pop(formatNames.index(removeFormatName))
    
    
  formatInfoDict = {}
  subClassedFunctions = {}

  for formatName in formatNames:
  
    formatLabel = allFormatsDict[formatName]
    
    formatModule = __import__('ccpnmr.format.converters.%sFormat' % formatLabel,{},{},'%sFormat' % formatLabel)
    formatClass = getattr(formatModule,'%sFormat' % formatLabel)
 
    formatObject = formatClass(ccpnProject)    
    
    formatObject.IOkeywords = copy.deepcopy(formatObject.IOkeywords)
        
    functionNames = formatObject.IOkeywords.keys()
    functionNames.sort()    
    
    subClasses = formatObject.IOkeywords['subClasses']
    functionNames.pop(functionNames.index('subClasses'))
    
    for functionName in functionNames[:]:
      for superClass in subClasses.keys():
        if functionName in ('read%s' % superClass,'write%s' % superClass):

          function = getattr(formatObject,functionName)
          addIOkeywordList = []
          setNonDictKeywords(function,addIOkeywordList,0)
          
          automaticKeyword = automaticKeywords[superClass]
            
          for subClass in subClasses[superClass]:
            subClassFunctionName = functionName.replace(superClass,subClass)
            
            if not subClassedFunctions.has_key(subClassFunctionName):
              subClassedFunctions[subClassFunctionName] = functionName

            formatObject.IOkeywords[subClassFunctionName].update(formatObject.IOkeywords[functionName])
            
            if formatObject.IOkeywords[subClassFunctionName].has_key(automaticKeyword):
              (default,mandatory,infoText) = formatObject.IOkeywords[subClassFunctionName][automaticKeyword]
              formatObject.IOkeywords[subClassFunctionName][automaticKeyword] = (default,mandatory,True,infoText)
 
            for (IOkeyword,default,mandatory,isAutomatic,infoText) in addIOkeywordList:
              formatObject.IOkeywords[subClassFunctionName][IOkeyword] = (default,mandatory,False,infoText)
      
          functionNames.pop(functionNames.index(functionName))

        
    functionNameInfo = {}
    functionWidgetInfo = {}
    fileNameInfo = {}
    
    for functionName in functionNames:
    
      isValidFunction = False
      
      for (component,readName,writeName,formatsDict,importWidgets,exportWidgets,void,void) in IoSetupList:
        validFunctionNames = [readName,writeName]
        
        if functionName in validFunctionNames:          
          functionIndex = validFunctionNames.index(functionName)
          if formatsDict.has_key(formatName):
          
            # Track 'default' file name
            fileNameInfo[functionName] = formatsDict[formatName][0]
        
            if formatsDict[formatName][functionIndex + 1]:
              isValidFunction = True
              
              # Now get widget information
              if not functionWidgetInfo.has_key(functionName):
                functionWidgetInfo[functionName] = {}
              
              # Generic info
              if functionIndex == 0:
                widgetInfo = importWidgets
              else:
                widgetInfo = exportWidgets
                
              getWidgetInfo(functionWidgetInfo[functionName],widgetInfo)
              
              # Format specific info
              getWidgetInfo(functionWidgetInfo[functionName],formatsDict[formatName][functionIndex + 3])                
              
              break

              
      if not isValidFunction:
        continue
                
      
      IOkeywords = formatObject.IOkeywords[functionName].keys()
      IOkeywords.sort()
      
      numMandatory = 0
      
      IOkeywordList = []
      
      for IOkeyword in IOkeywords:
 
        if len(formatObject.IOkeywords[functionName][IOkeyword]) == 3:
          (default,mandatory,infoText) = formatObject.IOkeywords[functionName][IOkeyword]
          isAutomatic = False
        else:
          (default,mandatory,isAutomatic,infoText) = formatObject.IOkeywords[functionName][IOkeyword]

        if mandatory:
                  
          IOkeywordList.insert(numMandatory,(IOkeyword,default,mandatory,isAutomatic,infoText))
          numMandatory += 1
          
        else:
        
          IOkeywordList.append((IOkeyword,default,mandatory,isAutomatic,infoText))
          
      # Add info (if required)
      
      function = getattr(formatObject,functionName)
      
      setNonDictKeywords(function,IOkeywordList,numMandatory)
      
      if functionName == readName:
        actionType = 'import'
      elif functionName == writeName:
        actionType = 'export'
        
      if not functionNameInfo.has_key(actionType):
        functionNameInfo[actionType] = {}

      functionNameInfo[actionType][component] = (functionName,IOkeywordList)
    
    #
    # Set the info
    #
    
    formatInfoDict[formatName] = {}
    
    for actionType in functionNameInfo.keys():

      formatInfoDict[formatName][actionType] = {}

      for component in functionNameInfo[actionType].keys():

        functionName = functionNameInfo[actionType][component][0]

        if fileNameInfo.has_key(functionName):
          defaultFileName = fileNameInfo[functionName]
        else:
          defaultFileName = 'noDefaultFile.xxx'

        formatInfoDict[formatName][actionType][component] = (functionName,defaultFileName,[])

        for (IOkeyword,default,mandatory,isAutomatic,infoText) in functionNameInfo[actionType][component][1]:
          if functionWidgetInfo[functionName].has_key(IOkeyword):
            (widgetType,defaults,label) = functionWidgetInfo[functionName][IOkeyword]
          else:
            (widgetType,defaults,label) = (None,None,None)
            
          formatInfoDict[formatName][actionType][component][2].append((IOkeyword,default,mandatory,isAutomatic,infoText,widgetType,defaults,label))

  return (formatInfoDict,subClassedFunctions)

#
# Gives a list of valid format names for exporting, and a dictionary with information on CCPN objects
# that are relevant for the export (so a selection list can be made).
#

def getValidExportFormats(project):
  
  #
  # Sequence information
  #
  # this works for a single molSystem, with a single chain, single peak list and single chemical shift list
  # the setup of a normal format to format conversion without existing ccpn project
  
  ccpnObjects={'MolSystem':     [],
               'Chain':         [],
               'NmrProject':    [],
               'NmrExperiment': [],
               'DataSource':    [],
               'PeakList':      [],
               'StructureEnsemble': [],
               'Entry': [],
               # Measurement and Constraint list stuff is set on the fly
               }
  
  hasChains = False
    
  for molSystem in project.molSystems:
    ccpnObjects['MolSystem'].append(molSystem)
    if molSystem.chains:
      hasChains = True
      for chain in molSystem.chains:
        ccpnObjects['Chain'].append(chain)
  
  #
  # Coordinates
  #
  
  hasCoordinates = False
   
  for structureEnsemble in project.structureEnsembles:    
    for model in structureEnsemble.models:
      if model.coords:
        hasCoordinates = True
        break
    
    ccpnObjects['StructureEnsemble'].append(structureEnsemble)
  
  #
  # ChemComps
  #
  
  # TODO
  
  #
  # NMR information
  #
  
  hasPeakLists = False
  
  measurementListClassNames = set(('ShiftList',))
  hasMeasurements = {}
  
  for nmrProject in project.nmrProjects:
    #print 'dataDimRefs: ',nmrProject.dataDimRefs
    ccpnObjects['NmrProject'].append(nmrProject)
    #
    # Peak lists and peak assignments
    #
    nExp=0
    for nmrExperiment in nmrProject.experiments:
      ccpnObjects['NmrExperiment'].append(nmrExperiment)
      nExp+=1
      nDataSources=0
      #print 'experiment %i %s'%(nExp,str(nmrExperiment))
      #print 'dataDimRefs: ',nmrExperiment.dataDimRefs
      for dataSource in nmrExperiment.dataSources:
        #print 'dataDimRefs: ',dataSource.dataDimRefs
        ccpnObjects['DataSource'].append(dataSource)  
        nDataSources+=1
        #print 'dataSource %i: %s'%(nDataSources,str(dataSource))
        for peakList in dataSource.peakLists:
          #print 'dataDimRefs: ',peakList.dataDimRefs
          
          if peakList.peaks:
            ccpnObjects['PeakList'].append(peakList)
            hasPeakLists = True

    #
    # Measurements (only shifts, currently)
    #

    for className in measurementListClassNames:
      for measurementList in nmrProject.findAllMeasurementLists(className = className):
        if measurementList.measurements:
          if not ccpnObjects.has_key(className):
            ccpnObjects[className] = []
          ccpnObjects[className].append(measurementList)
          hasMeasurements[measurementList.findFirstMeasurement().className] = True

  #
  # Constraint lists, by type
  #
  
  constraintListClassNames = set(('DistanceConstraintList','HBondConstraintList','DihedralConstraintList','RdcConstraintList',))
  hasConstraints = {}
  
  for nmrConstraintStore in project.nmrConstraintStores:
  
    #
    # Constraints
    #

    for className in constraintListClassNames:
    
      if hasConstraints.has_key(className):
        continue
    
      for constraintList in  nmrConstraintStore.findAllConstraintLists(className = className):
        if constraintList.constraints:
          if not ccpnObjects.has_key(className):
            ccpnObjects[className] = []
          ccpnObjects[className].append(constraintList)
          hasConstraints[constraintList.findFirstConstraint().className] = True 
  
  #
  # Project writing
  #
  
  hasProjectInfo = False
  
  #
  # More sophisticated determination if objects available for project export
  #
  
  hasProjectInfoPerFormat = {}
  
  for objectType in projectExportRequirements.keys():
    
    hasInfo = True
    if objectType == 'chain' and not hasChains:
      hasInfo = False
    elif objectType == 'molSystem' and not hasChains:
      hasInfo = False
    elif objectType == 'shiftList' and not hasMeasurements.has_key('Shift'):
      hasInfo = False
    elif objectType == 'entry':
      hasInfo = False
      for nmrEntryStore in project.nmrEntryStores:
        for entry in nmrEntryStore.sortedEntries():
          ccpnObjects['Entry'].append(entry)
          hasInfo = True
        
    for formatName in projectExportRequirements[objectType]:
      if not hasProjectInfoPerFormat.has_key(formatName):
        hasProjectInfoPerFormat[formatName] = hasInfo
      
      # Reset if something is missing.
      if hasProjectInfoPerFormat[formatName] and not hasInfo:
        hasProjectInfoPerFormat[formatName] = hasInfo 
      
  #
  # Now set the information
  #
  
  (fcInfo,subClassedFunctions) = getFormatConverterInfoDict()  
  
  exportFormats = {}
  
  for formatName in fcInfo.keys():
    
    exportFormats[formatName] = []
    
    actionTypes = fcInfo[formatName].keys()
    
    if 'export' in actionTypes:
    
      components = fcInfo[formatName]['export'].keys()

      if hasChains and 'sequence' in components:
        exportFormats[formatName].append('sequence')

      if hasCoordinates and 'coordinates' in components:
        exportFormats[formatName].append('coordinates')

      if hasPeakLists and 'peaks' in components:
        exportFormats[formatName].append('peaks')

      if 'project' in components and hasProjectInfoPerFormat.has_key(formatName) and hasProjectInfoPerFormat[formatName]:
        exportFormats[formatName].append('project')

      for measurementClass in hasMeasurements.keys():
        funcName = "write%ss" % measurementClass
        for component in components:          
          if funcName == fcInfo[formatName]['export'][component][0]:
            exportFormats[formatName].append(component)
            break


      for constraintClass in hasConstraints.keys():
        funcName = "write%ss" % constraintClass
        for component in components:          
          if funcName == fcInfo[formatName]['export'][component][0]:
            exportFormats[formatName].append(component)
            break

    if not exportFormats[formatName]:
      del(exportFormats[formatName])


  return exportFormats, ccpnObjects

#
# Get the regular expressions
#

from ccp.format.general.Util import getRegularExpressions
ccpPatt = getRegularExpressions()

def getResName(chainCode,seqCode,atomName,seqInsertCode = '',spinSystemId = None, fixChainCode = True):
  """ NBNB TBD Rsmus July 2012:
  What happens if atomName is an empty string?
  Can come about e.g. in Ansig if only sequence is asssigned for a PeakDim
  """
  if seqCode == None and spinSystemId != None:
    seqCode = "SS%d" % spinSystemId
  
  if None in (chainCode,seqCode):

    print "  Error creating resName: chainCode or seqCode element is None ('%s')" % (str((chainCode,seqCode)))
    resName = None

  else:
  
    if fixChainCode:
      chainCode = '%s' % chainCode
      if chainCode.count(tagSep):
        oldChainCode = chainCode
        chainCode = chainCode.replace(tagSep,'_')
        print "  Warning: replacing illegal format chain code '%s' by '%s'." % (oldChainCode,chainCode)

    # This to handle interresidue constraints
    if atomName == None:
      atomName = allResidueAtoms_kw
    
    resName = chainCode + tagSep + str(seqCode) + seqInsertCode.strip() + tagSep + atomName.upper()
  
    # Check if it doesn't have any illegal chars
    if len(resName.split(tagSep)) > 3:
      print "  Warning: string containing separator '%s' in '%s'" % (tagSep,resName)

  return resName
  

def createNonStereoName(atomSetName,code):

  # See if it's a group
  asterisk = atomSetName.find('*')
  
  if asterisk == -1:
    nonStereoName = atomSetName[:-1] + code

  elif asterisk == len(atomSetName) - 1:
    nonStereoName = atomSetName[:-2] + code + atomSetName[-1]
  
  else:
    print "  Warning can't create non stereo name for %s" % atomSetName
    nonStereoName = atomSetName
  
  return nonStereoName

def createSelection(ccpnObjectList):

  """
  Create selection - this creates a selectionList and Dict based on the CCPN
  object type...
  """

  selectionList = []
  selectionDict = {}

  if ccpnObjectList:
  
    if not type(ccpnObjectList) == type([]):
      ccpnObjectList = list(ccpnObjectList)
  
    ccpnRefObject = ccpnObjectList[0]

    if type(ccpnRefObject) in [Nmr.Experiment,
                               Affiliation.Organisation]:
      for obj in ccpnObjectList:
        label = "%d:%s" % (obj.serial,obj.name)
        selectionList.append(label)
        selectionDict[label] = obj

    elif type(ccpnRefObject) in (Nmr.StructureGeneration,Nmr.StructureAnalysis):
      
      if ccpnRefObject.className == 'StructureGeneration':
        prefix = 'Structure generation'
      else:
        prefix = 'Structure analysis'
      
      for obj in ccpnObjectList:
        if obj.name:
          label = "%s: %s" % (prefix,obj.name)
        else:
          label = '%s cycle #%d' % (prefix,obj.serial)
        
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) == MolStructure.StructureEnsemble:
      for obj in ccpnObjectList:
        label = 'Structure ensemble #%d (%d models)' % (obj.ensembleId,len(obj.coordChains))
        
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) == NmrConstraint.NmrConstraintStore:
      for obj in ccpnObjectList:
        # Add info on connected structure generation/analysis
        strucGenTexts = []
        for strucGen in obj.structureGenerations:
          if strucGen.name:
            strucGenTexts.append("'%s'" % strucGen.name)
            
        strucAnalTexts = []
        for strucAnal in obj.structureAnalyses:
          if strucAnal.name:
            strucAnalTexts.append("'%s'" % strucAnal.name)

        addText = ""
        if strucGenTexts:
          addText+= "Structure generation(s): %s" % (", ".join(strucGenTexts))
        if strucAnalTexts:
          if addText:
            addText += ', '
          addText+= "Structure analysis: %s" % (", ".join(strucAnalTexts))
        if addText:
          addText = " (%s)" % addText  
        
        label = 'Constraint set #%d%s' % (obj.serial,addText)
        
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) == Nmr.DataSource:
      for obj in ccpnObjectList:
        label = '%s:%d:%s (%s)' % (obj.experiment.name,obj.serial,obj.name,obj.dataType)
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) == Nmr.PeakList:
      for obj in ccpnObjectList:
        ds = obj.dataSource
        label = '%s:%d:%s:%d:%s' % (ds.experiment.name,ds.serial,ds.name,obj.serial,obj.name)
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) in [NmrConstraint.RdcConstraintList,
                                 NmrConstraint.DistanceConstraintList,
                                 NmrConstraint.ChemShiftConstraintList,
                                 NmrConstraint.HBondConstraintList,
                                 NmrConstraint.DihedralConstraintList]:
      
      for obj in ccpnObjectList:
        label = '%d:%d:%s' % (obj.parent.serial,obj.serial,obj.name)
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) == Molecule.Molecule:
      for obj in ccpnObjectList:
        label = obj.name
        selectionList.append(label)
        selectionDict[label] = obj

    elif type(ccpnRefObject) == MolSystem.MolSystem:
      for obj in ccpnObjectList:
        label = obj.code
        selectionList.append(label)
        selectionDict[label] = obj
        
    elif type(ccpnRefObject) == MolSystem.Chain:
      for obj in ccpnObjectList:
        label = "%s:'%s'" % (obj.molSystem.code,obj.code)
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) == MolStructure.Model:
      for obj in ccpnObjectList:
        label = "%s:%s:model_%d" % (obj.parent.molSystem.code,obj.parent.ensembleId,obj.serial)
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) == ChemComp.ChemCompVar:
      for obj in ccpnObjectList:
        label = "%s,%s:%s,%s" % (obj.chemComp.molType,obj.chemComp.ccpCode,obj.linking,obj.descriptor)
        
        # Could get the same chemComp multiple times...    
        if not label in selectionList:
          selectionList.append(label)
          selectionDict[label] = obj
        
    elif type(ccpnRefObject) == Instrument.NmrSpectrometer:
      for obj in ccpnObjectList:
        label = "%d:%s,%s (%s)" % (obj.serial,obj.manufacturer,obj.model,obj.nominalFreq)
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) == Affiliation.Person:
      for obj in ccpnObjectList:
        label = "%s, %s" % (obj.familyName,obj.givenName)
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) in [NmrEntry.Entry,NmrEntry.Study]:
      for obj in ccpnObjectList:
        label = obj.name
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) == Method.Software:
      for obj in ccpnObjectList:
        label = "%s (version %s)" % (obj.name,obj.version)
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) == Method.Method:
      for obj in ccpnObjectList:
        label = "%d: %s (task: '%s', procedure: '%s')" % (obj.serial,obj.name,obj.task,obj.procedure)
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) == Affiliation.Group:
      for obj in ccpnObjectList:
        label = "%s (%s)" % (obj.name,obj.parent.name)
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) == Affiliation.PersonInGroup:
      for obj in ccpnObjectList:
        label = "%s, %s (%s)" % (obj.person.familyName,obj.person.givenName,obj.group.name)
        selectionList.append(label)
        selectionDict[label] = obj  

    elif type(ccpnRefObject) in [Citation.BookCitation,
                                 Citation.ConferenceCitation,
                                 Citation.JournalCitation,
                                 Citation.ThesisCitation]:
      for obj in ccpnObjectList:
        refText = ""
        if obj.firstPage:
          refText += " %s" % obj.firstPage
          if obj.lastPage:
            refText += "-%s" % obj.lastPage
        if obj.year:
          refText += " (%d)" % obj.year
          
        label = "%d:%s '%s'%s" % (obj.serial,obj.className[:-8],obj.title,refText)
        selectionList.append(label)
        selectionDict[label] = obj
                               
    elif hasattr(ccpnRefObject,'measurements'): # This is a hack to get at all measurement lists
      
      for obj in ccpnObjectList:
        label = '%d:%s' % (obj.serial,obj.name)
        selectionList.append(label)
        selectionDict[label] = obj  
        
  return (selectionList,selectionDict)

def setCcpNmrResonanceName(resonance,refAtom,atomName = None, forceChange = 0):
  
  #
  # Set the CCP names
  #
  # Full 'resonance' name on application data level
  # Just atom (group) name on resonance level!
  #
  
  residue = refAtom.residue
    
  seqCode = residue.seqCode
  chainCode = residue.chain.code

  if atomName == None:
    atomName = refAtom.name
  
  resName = chainCode + tagSep + str(seqCode) + tagSep + atomName
  
  # Set resonance name
  if resonance.name == None or forceChange == 1:
    resonance.name = atomName
    
    # Here only look for ONE match - don't want multiple matches for ccp names...

    assignData = resonance.findFirstApplicationData(application = ccpNmr_kw, keyword = assign_kw)
    
    if assignData:
      if assignData.value != resName:
        if forceChange == 0:
          print "  Warning: couldn't set ccpNmr resonance name %s (name already set to %s)" % (resName,assignData.value)
        else:
          assignData.delete()
          resonance.addApplicationData(Implementation.AppDataString( application = ccpNmr_kw, keyword = assign_kw, value = resName))
          
    else:
      resonance.addApplicationData(Implementation.AppDataString( application = ccpNmr_kw, keyword = assign_kw, value = resName))

  else:
    if resonance.name != atomName:
      print "  Warning: couldn't set ccpNmr name %s for resonance %s (name already set)" % (atomName,resonance.name)
      resName = resonance.name
      
  return resName

def getStereoCode(resonance,assignFormat):
  
  stereoCode = resonance.findFirstApplicationData(application = assignFormat, keyword = 'ambiguityCode')

  if stereoCode:
    if assignFormat == 'nmrStar':
      return stereoCode.value
  
  else:
    return None

def getIsotopes(atomType):
  if atomType.upper() == 'X':
    return ['13C','15N']
  elif atomType.upper() == 'H':
    return ['1H']
  elif atomType.upper() == 'D':
    return ['2H']
  elif atomType.upper() == 'C':
    return ['13C']
  elif atomType.upper() == 'N':
    return ['15N']
  else:
    print "  Unknown atomType %s!" % atomType
    return []

def getNmrIsotope(proj,elementSymbol):

  chemElement = None
  
  if not proj.currentChemElementStore:
    proj.currentChemElementStore = proj.findFirstChemElementStore()
    
  chemElement = proj.currentChemElementStore.findFirstChemElement(symbol = elementSymbol)
    
  isotopeCode = 'unknown'
  
  if chemElement:
    
    isotopes = chemElement.sortedIsotopes()
    
    if isotopes:
    
      #
      # For common ones pick the relevant NMR active one
      #
      
      if chemElement.symbol in commonNmrActiveIsotopes.keys():
      
        selectedIsotope = chemElement.findFirstIsotope(massNumber = commonNmrActiveIsotopes[chemElement.symbol])
    
      elif len(isotopes) > 1:

        #
        # Silently pick most abundant one
        #
        #
        
        maxAbundance = 0
        
        for isotope in isotopes:
        
          if isotope.abundance > maxAbundance:
            maxAbundance = isotope.abundance
            selectedIsotope = isotope
      
      else:
        
        selectedIsotope = isotopes[0]
      
      isotopeCode = str(selectedIsotope.massNumber) + elementSymbol
    
    else:
    
      print "  Element symbol %s has no nmr active isotopes" % elementSymbol
  
  else:
  
    print "  Unknown element symbol %s." % elementSymbol

  return isotopeCode

def getAtomCode(isotope):
  if isotope.upper() == '1H':
    return ['h']
  elif isotope.upper() == '2H':
    return ['d']
  elif isotope.upper() == '13C':
    return ['c']
  elif isotope.upper() == '15N':
    return ['n']
  else:
    print "  Unknown isotope %s!" % isotope
    return []
 

def createFixedResonance(strucGen,resonance):

  return NmrConstraint.FixedResonance(

              strucGen.nmrConstraintStore,
              resonanceSerial = resonance.serial,
              isotopeCode =     resonance.isotopeCode,
              name =            resonance.name,
              details =         resonance.details

                                    )

def getApplResNames(format,resonances, keyword = assign_kw):
  
  resonanceNames = {}
  
  for resonance in resonances:
    applData = resonance.findAllApplicationData(application = format, keyword = keyword)
    
    for applDatum in applData:
      resName = applDatum.value
      
      if not resonanceNames.has_key(resName):
        resonanceNames[resName] = []
       
      resonanceNames[resName].append(resonance)
    
  return resonanceNames
          
def getResNameText(resonance):

  resonanceNames = getApplResNames(ccpNmr_kw,[resonance])
 
  if resonanceNames:
    resonanceNameText = str(resonanceNames.keys()) + " (ccpNmr names)"

  else:
    applNames = resonance.findAllApplicationData(keyword = assign_kw)
    
    if applNames:
    
      resonanceNameText = ""
      
      for applName in applNames:
      
        resonanceNameText += "'" + applName.value + "' (%s names), " % applName.application
     
      resonanceNameText = resonanceNameText[:-2]
      
    else:
     
      if resonance.name:
        resonanceNameText = resonance.name + " (resonance name)"
      else:
        resonanceNameText = 'Unknown'
  
  return resonanceNameText

def getAtomNameText(atom):

  atomName = atom.name
  seqCode = atom.residue.seqCode
  insertCode = atom.residue.seqInsertCode
  chainCode = atom.residue.chain.code
  
  atomNameText = "%s.%d%s.%s" % (chainCode,seqCode,insertCode.strip(),atomName)

  return atomNameText

def getNameInfo(resName, verbose = 1):
  
  #
  # Assumes chainCode(molCode?).seqCode(seqInsertCode).atomName format
  # TODO: Should this be format-dependent?
  #

  try:
    (chainCode,seqCode,atomName) = resName.split(tagSep)
    spinSystemId = None
    
    if seqCode[:2] == 'SS':
      (spinSystemId,seqInsertCode) = getSeqAndInsertCode(seqCode[2:])
      seqCode = None
    else:
      (seqCode,seqInsertCode) = getSeqAndInsertCode(seqCode)
      spinSystemId = None
    
    if seqCode != None or spinSystemId != None:
      if atomName.strip() == "":
        atomName = None
      else:
        atomName = atomName.upper()
    else:
      if verbose:
        print "  Error: getNameInfo Can't decompose resonance name '%s'." % (resName) 
      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = (None,None,None,defaultSeqInsertCode,None)
    
  except:
    if verbose:
      print "  Error: getNameInfo Can't decompose resonance name '%s'." % (resName) 
    (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = (None,None,None,defaultSeqInsertCode,None)
    
  return (chainCode,seqCode,spinSystemId,seqInsertCode,atomName)

def getShiftRange(shift):

  shiftValue = shift.value
  shiftError = shift.error
  
  if shiftError == None:
    shiftError = 0.002
  
  return (shiftValue,shiftValue-shiftError,shiftValue+shiftError)
  
def getAtomSet(atomSetsAtomDict,atomName,atomsFound):
  
  #
  # Return atomSet if exists
  #
  # Have to resort every time to compare since atom order in data
  # model might be different from order in which atoms were found...
  #
  
  atomSet = None
  
  atomsFound.sort()
  
  for checkAtomSet in atomSetsAtomDict.keys():
  
    if atomsFound == atomSetsAtomDict[checkAtomSet]:
    
      atomSet = checkAtomSet
      
      if atomName != atomSet.name:
      
        print " Atomset %s with same atoms and different name (%s) already exists." % (atomName,atomSet.name)
      
      break
      
  return atomSet

def returnProchiralAtom(atom,chemCompVar = None):

  chemAtom = atom.chemAtom
  chemAtomSet = chemAtom.chemAtomSet
    
  # Only prochiral chemAtomSets!
  if chemAtomSet != None and chemAtomSet.isProchiral == True:
    
    # Loop over other atoms for this residue and see if there's a match...
    for otherAtom in atom.residue.sortedAtoms():
      if otherAtom != atom:
          
        if otherAtom.chemAtom in chemAtomSet.chemAtoms:
          return otherAtom
          
  return None

def returnProchiralAtomGroup(atoms,chemCompVar = None):

  chemAtom = atoms[0].chemAtom
  chemAtomSet = chemAtom.chemAtomSet
  
  if not chemAtomSet:
    return (None,None)
  
  deepChemAtomSet = chemAtomSet.chemAtomSet

  if deepChemAtomSet:
    if deepChemAtomSet.isProchiral == True:
      otherAtoms = []
      
      # Get other chematomset
      chemAtomSetList = deepChemAtomSet.sortedChemAtomSets()
      if chemAtomSetList[0] == chemAtomSet:
        otherChemAtomSet = chemAtomSetList[1]
      else:
        otherChemAtomSet = chemAtomSetList[0]
     
      # The name for the other 'atom' (is really an atomset here)
      # comes from the otherchematomset
      otherAtomName = otherChemAtomSet.name

      # Loop over other atoms for this residue and see which ones match...
      otherChemAtomSetAtoms = otherChemAtomSet.sortedChemAtoms()
     
      for otherAtom in atoms[0].residue.sortedAtoms():
        if not otherAtom in atoms:
  
          if otherAtom.chemAtom in otherChemAtomSetAtoms:
            otherAtoms.append(otherAtom)

      return (otherAtoms,otherAtomName)
  
  return (None,None)


def updateResonanceNamesDict(resonance,resonanceNames,assignFormat = None):
  
  resonanceNamesFound = []
  
  # Get application data name for resonance (if any)
  if assignFormat != None:
    assignData = resonance.findAllApplicationData(application = assignFormat, keyword = assign_kw)
  else:
    assignData = ()

  # Set ccpn format resonance name if available
  if not assignData:
    assignData = resonance.findAllApplicationData(application = ccpNmr_kw, keyword = assign_kw)

    if not assignData:
      # Default to actual resonance name if nothing else available??!
      resonanceNamesFound.append(resonance.name)
      
  if assignData:
    for assignDatum in assignData:
      resonanceNamesFound.append(assignDatum.value)
  

  for resonanceNameFound in resonanceNamesFound:
  
    if not resonanceNames.has_key(resonanceNameFound):
      resonanceNames[resonanceNameFound] = []
  
    resonanceNames[resonanceNameFound].append(resonance)

  return resonanceNamesFound

##########################################
# Functions to handle original data tags #
##########################################

def getOriginalData(format,parent,rawData,formatName,setFormatName = None):

  #
  # This could become generic - could be used for other formats as well
  #

  if not setFormatName:

    setFormatName = formatName

  keyword = 'orig' + formatName[0].upper() + formatName[1:]
  applData = parent.findFirstApplicationData(application = format, keyword = keyword)

  if applData:

    #
    # Bit dangerous but because star makes everything strings it shouldn't matter...
    # is mostly for setting None!
    #

    try:

      value = eval(applData.value)

      if type(value) == type([]) and hasattr(rawData,setFormatName) and type(getattr(rawData,setFormatName)) == type([]):

        getattr(rawData,setFormatName).extend(value)

      else:

        setattr(rawData,setFormatName,value)

    except:     

      setattr(rawData,setFormatName,applData.value)

    return True

  return False

def getOriginalDataValue(format,parent,formatName,setFormatName = None):

  #
  # This could become generic - could be used for other formats as well
  #
  
  value = None

  if not setFormatName:

    setFormatName = formatName

  keyword = 'orig' + formatName[0].upper() + formatName[1:]
  applData = parent.findFirstApplicationData(application = format, keyword = keyword)

  if applData:

    value = applData.value

  return value

def setOriginalData(format,parent,rawData,name,isUnique = True):

  #
  # Code to set original data in the data model as application data
  #

  if hasattr(rawData,name):

    keyword = 'orig' + name[0].upper() + name[1:]
    value = getattr(rawData,name)
    appDataClass = None
    
    if value != None:

      if type(returnInt(value,default=None,verbose=0)) == type(value):
        appDataClass = 'AppDataInt'
      elif type(returnFloat(value,default=None,verbose=0)) == type(value):
        appDataClass = 'AppDataFloat'
      elif type(returnLong(value,default=None,verbose=0)) == type(value):
        appDataClass = 'AppDataLong'
    
    if not appDataClass:
      appDataClass = 'AppDataString'
      value = str(value)
    
    if isUnique:
      setUniqueAppData(appDataClass,parent,format,keyword,value)
    else:
      appData = getattr(Implementation,appDataClass)(application = format,keyword = keyword,value = value)
      parent.addApplicationData(appData)

    return True

  return False


###########################################################
# Experiment setup - now based on reference data in model #
###########################################################

def getRefExpFromOldExpType(project,expType):

  refExperiment = None
  
  if oldToNewExpDict.has_key(expType):
    (expPrototypeName,refExpName) = oldToNewExpDict[expType]
    
    setCurrentStore(project,'ChemElementStore')
        
    expPrototype = project.findFirstNmrExpPrototype(name = expPrototypeName)
    
    if expPrototype:
      refExperiment = expPrototype.findFirstRefExperiment(name = refExpName)

  return refExperiment

def createExperiment(proj,expName,refExperiment,unit = 'ppm'):
  
  #
  # Create a default experiment based on a reference experiment... 
  #

  numDim = len(refExperiment.refExpDims)
  isReversed = refExperiment.isReversed
  
  # Make sure experiment name fits Ccp Line
  expName = returnMemopsLine(expName)

  # Create the experiment
  exp = Nmr.Experiment(proj.currentNmrProject,name = expName, numDim = numDim, refExperiment = refExperiment)
  
  # Tag to mark that this is a 'default' experiment
  exp.addApplicationData(Implementation.AppDataBoolean(application = ccpNmr_kw, keyword = 'default', value = True))

  #
  # Set expDimRefs
  #
  
  atomSiteList = []
  acquisitionDim = None
  
  for i in range(0,numDim):
  
    expDim = exp.findFirstExpDim(dim = i + 1)
   
    refExpDim = refExperiment.findFirstRefExpDim(dim = i + 1)
    expDim.refExpDim = refExpDim

    atomSiteList.append([])

    for refExpDimRef in refExpDim.refExpDimRefs:

      # TODO Here catch things like ConstantTimePeriod for editing!!!!?!?!
      
      expMeasurement = refExpDimRef.expMeasurement
      
      isotopeCodes = []
      for atomSite in expMeasurement.atomSites:
        isotopeCodes.append(atomSite.isotopeCode)

      if isotopeCodes.count('1H'):
        sf = baseFrequency = 500.0
      elif isotopeCodes.count('13C'):
        sf = baseFrequency = 125.7367377
      elif isotopeCodes.count('15N'):
        sf = baseFrequency = 50.670637
      else:
        sf = baseFrequency = 1.0
        
      #
      # Setting paramaters correctly goes via GUI interface... these are 'defaults'
      #

      expDimRef = Nmr.ExpDimRef(expDim, unit = unit,
                                        sf = sf,
                                        baseFrequency = baseFrequency,
                                        measurementType = expMeasurement.measurementType,
                                        isotopeCodes = isotopeCodes,
                                        refExpDimRef = refExpDimRef)

      # Tag to mark that this is a 'default' expDimRef
      expDimRef.addApplicationData(Implementation.AppDataBoolean(application = ccpNmr_kw, keyword = 'default', value = True))

      #
      # Keep track for expTransfer settings...
      #
      
      atomSiteList[-1].append((expDimRef,expMeasurement.atomSites))

  #
  # Set acquisition dimension. Is the expMeasurement connected to the last
  # expStep (or first expStep in case isReversed).
  #
  
  expDim = None
  
  for expGraph in refExperiment.nmrExpPrototype.expGraphs:

    expStepNums = []
    for expStep in expGraph.expSteps:
      expStepNums.append(expStep.stepNumber)
    expStepNums.sort()

    if not isReversed:
      expStepNums.reverse()
    
    for expStepNum in expStepNums:
      expStep = expGraph.findFirstExpStep(stepNumber = expStepNum)

      for refExpDimRef in expStep.expMeasurement.refExpDimRefs:
         
        refExpDim = refExpDimRef.refExpDim
        expDim = exp.findFirstExpDim(refExpDim = refExpDim)

        if expDim:
          break
      
      if expDim:
        break
    
    if expDim:
      break
  
  if expDim:  
    expDim.isAcquisition = True
    
  #
  # Set expTransfers...
  #
  
  for i in range(0,len(atomSiteList)-1):

    for (expDimRef,atomSites) in atomSiteList[i]:
      for otherAtomSiteList in atomSiteList[i+1:]:
        for (otherExpDimRef,otherAtomSites) in otherAtomSiteList:
          for expGraph in atomSites[0].nmrExpPrototype.expGraphs:

            #
            # Transfers are only set in simple case: direct transfer between single atomSites corresponding to expDimRefs
            #

            expDimRefSet = set([expDimRef,otherExpDimRef])
            
            sites = frozenset([atomSites[0],otherAtomSites[0]])
            refExpTransfer = expGraph.findFirstExpTransfer(atomSites = sites)

            if refExpTransfer:
              expTransfer = exp.findFirstExpTransfer(transferType = refExpTransfer.transferType, expDimRefs = expDimRefSet)
              if not expTransfer:
                expTransfer = exp.newExpTransfer(transferType = refExpTransfer.transferType, expDimRefs = expDimRefSet)        

  return exp
  
def createPpmFreqDataSource(exp,dsName,dsType,numDim):

  ds = Nmr.DataSource(exp,name = dsName, numDim = numDim, dataType = dsType)

  for dim in range(0,numDim):
    expDim = exp.findFirstExpDim(dim = (dim + 1))
    expDimRefs = expDim.sortedExpDimRefs()
    refIsotopeCode = expDimRefs[0].isotopeCodes[0]
    
    # User first expDimRef to set number of points
    # Set numPoints, numPointsOrig to 20 for H, 250 for C,N!
    if refIsotopeCode == '1H':
      numPoints = 1024
      valuePerPoint = 6.5
    else:
      numPoints = 256
      if refIsotopeCode == '13C':
        valuePerPoint = 100.0
      elif refIsotopeCode== '15N':
        valuePerPoint = 40.0
      else:
        valuePerPoint = 1.0

    # Setting up empty freqdatadim and associated datadimref...
    # Assuming that imaginary data will be deleted!
    freqDataDim = Nmr.FreqDataDim(ds, expDim = expDim, dim = (dim + 1), isComplex = True,valuePerPoint = valuePerPoint, numPoints = numPoints, numPointsOrig = numPoints)
      
    # Tag to mark that this is a 'default' freqDataDim
    freqDataDim.addApplicationData(Implementation.AppDataBoolean(application = ccpNmr_kw, keyword = 'default', value = True))
  
  # Set sensible dataDimRef info.
  setDataDimRefs(ds)

  return ds

def setDataDimRefs(dsFreq):
  
  # TODO: base this on temperature!?!
  
  protonFreqRefPpm = 4.773  
  otherFreqRefPpm = {'13C': 100.0, '15N': 120.0} # Arbitrary!

  # Get the default center proton frequency
  
  protonFreq = 0
  otherFreqs = {}
  
  protonFreqRefZero = 0.0
  
  for dataDim in dsFreq.dataDims:
    expDim = dataDim.expDim  
    
    for expDimRef in expDim.expDimRefs:
      isotope = expDimRef.isotopeCodes[0]
      
      # Set default for first proton dim found
      # TODO: Should really have a popup here!
      if isotope == '1H':
        protonFreq = expDimRef.sf
        break
      elif isotope not in otherFreqs.keys():
        otherFreqs[isotope] = expDimRef.sf
    
    if protonFreq != 0:
      break
  
  if protonFreq == 0:
  
    if otherFreqs:
      isotopeCodes = otherFreqs.keys()
      isotopeCodes.sort()
      
      for isotopeCode in isotopeCodes:
        
        if otherFreqRefPpm.has_key(isotopeCode) and INDIRECT_FREQ_CONV.has_key(isotopeCode):
          print "  No 1H expDimRef found, trying to use %s dimension for referencing." % isotopeCode
          heteroFreqRefPpm = otherFreqRefPpm[isotopeCode]
          heteroFreq = otherFreqs[isotopeCode]
          heteroFreqRefZero = heteroFreq / (heteroFreqRefPpm / 1000000 + 1)
          protonFreqRefZero = heteroFreqRefZero * 100.0 / INDIRECT_FREQ_CONV[isotopeCode]
          break
    
    if not protonFreqRefZero:
      print "  No relevant expDimRef info found for referencing - aborting."
      return

  else:
    protonFreqRefZero = protonFreq / (protonFreqRefPpm / 1000000 + 1)
    
  for dataDim in dsFreq.dataDims:
    expDim = dataDim.expDim
    
    for expDimRef in expDim.expDimRefs:
      isotope = expDimRef.isotopeCodes[0]
      
      freqRefZero = protonFreqRefZero * INDIRECT_FREQ_CONV[isotope] / 100
        
      freqRefPpm = ((expDimRef.sf - freqRefZero) / freqRefZero) * 1000000
      freqRefPoint = returnInt(dataDim.numPoints / 2) + 1
        
      dataDimRef = Nmr.DataDimRef(dataDim,refPoint = freqRefPoint, refValue = freqRefPpm, expDimRef = expDimRef)

def getAzaraBlockSizes(dataDims):

  blockSize = []
  numDim = len(dataDims)
  
  #
  # Set default blockSize
  #
  
  for i in range(0,numDim):
    blockSize.append(2)
    
  
  s = 2 ** numDim  # This is also assuming default blockSize is 2
    
  #
  # Loop while s < 4096
  #
  
  while (s < 4096):
  
    rMax = 0
    
    for i in range(0,numDim):
    
      r = dataDims[i].numPointsValid / blockSize[i]
      
      if r > rMax:
        
        iMax = i
        rMax = r
      
    blockSize[iMax] *= 2

    s *= 2
    
  return blockSize


def getNominalRef(dataDimRef):

  freqDataDim = dataDimRef.dataDim
  
  # This assumes that:
  #   - the first point is 1 and is at highest frequency (sw)
  #   - the the last point is (npoints + 1) and is 0 Hz 
  #   - Works correctly with nmrPipe and Azara
  
  firstPoint = 1

  centerPoint = (freqDataDim.numPointsOrig) / 2  + firstPoint

  centerPpmValue = ((dataDimRef.refPoint - centerPoint) * dataDimRef.valuePerPoint) + dataDimRef.refValue
  
  return centerPpmValue

def getSoftwareDefinitions(project):

  for softwareName in softwareDefinitions.keys():
  
    softwareDefs = softwareDefinitions[softwareName]
    
    if not project.currentMethodStore:
      methodStore = project.newMethodStore(name = project.name)
    
    software = project.currentMethodStore.findFirstSoftware(name = softwareName)
    
    if not software:
    
      software = Method.Software(project.currentMethodStore, name = softwareName, version = softwareDefs['version'])
      software.tasks = softwareDefs['tasks']
    
    for taskName in softwareDefs['methods'].keys():
    
      for procedureName in softwareDefs['methods'][taskName]:
    
        method = project.currentMethodStore.findFirstMethod(software = software, task = taskName, procedure = procedureName)
      
        if not method:
      
          method = Method.Method(project.currentMethodStore, software = software, task = taskName, procedure = procedureName)

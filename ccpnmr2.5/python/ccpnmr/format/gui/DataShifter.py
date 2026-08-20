"""
======================COPYRIGHT/LICENSE START==========================

DataShifter.py: Allows copying of data between projects.

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

import Tkinter

import sys, traceback

from memops.universal.Io import joinPath

from memops.general.Io import loadProject

from memops.api import Implementation

from memops.gui.Label import Label
from memops.gui.BasePopup import BasePopup
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo
from memops.gui.Separator import Separator

from memops.editor.OpenProjectPopup import OpenProjectPopup

from ccp.general.Util import setCurrentStore

from ccpnmr.format.gui.ImportExportFormatPopup import setIoInfo, getIoInfo, getCcpnObjects
from ccpnmr.format.gui.FormatFilePopup import FormatFilePopup
from ccpnmr.format.general.userInteraction import setupMultiDialog, setupMessageReporter, setupDataEntry

from ccpnmr.format.general.Constants import assign_kw, newResonances_kw
from ccpnmr.format.general.Util import createSelection, getResName
from ccpnmr.format.general.Io import getHelpUrlDir

###################
# Copy info dicts #
###################

# 1:  full copy of linked objects
# 0:  only copy over attributes and relevant links to existing objects
# -1: check if object exists, if so then do '0', if not, then do '1'

copyOrderDict = {

  'Citation':     [('authors',0),('editors',0),('chemComps',0),('molecules',0),('nmrSpectrometers',0),('molSystems',0),('methods',0),('software',0),('molSystems',0),('experiments',0)],
  'DataSource':   [('dataDims',1),('dataDims.dataDimRefs',1),('compressMethod',-1),('processMethod',-1),('processedFrom',0),('processedTo',0),('dataStore',1)],
  'DataUrl':      [],
  'Entry':        [('measurementLists',0),('experiments',0),('molSystem',0),('relatedEntries',0),('software',-1),('structureGenerations',0),('structureGroups',0),('study',-1)],
  'Experiment':   [('expDims',0),('expDims.expDimRefs',1),('expTransfers',1)],
  'Molecule':     [('moleculeSysNames',1),('molSeqFragments',1),('molResidues',1),('molResidues.molResLinkEnds',0),('molResLinks',1)],
  'MolSystem':    [('molSystemSysNames',1),('chains',1),('chainInteractions',1),('chains.residues',0),('chains.chainFragments',-1),('chains.residues.atoms',0),('chains.residues.molSystemLinkEnds',1),('molSystemLinks',1)],
  'Method':       [('parameters',1)],
  'AbstractMeasurementList': [('measurements',1),('experiments',0)],
  'StructureGeneration':     [('structureEnsemble',1),('nmrConstraintStore',-1),('entries',0)],
  'NmrConstraintStore':      [('structureGenerations',0),('violationLists',1),('constraintLists',1),('constraintLists.constraints',1)], #('fixedResonances',1),('constraintLists.constraints.resonances',1)
  'NumericMatrix':   [('dataUrl',-1)],
  'Organisation': [('groups',1)],
  'PeakList':     [('peaks',1),('peaks.peakIntensities',1),('peaks.peakDims',0),('peaks.peakContribs',1),('peaks.peakDims.peakDimContribs',1),('peaks.peakDims.peakDimComponents',1)],
  'Person':       [('personInGroups',1),('authorCitations',-1),('editorCitations',-1)],
  'Software':     [('methods',-1)],

  }
  
classLinkNames = {

  'Method': 'methods',
  'Software': 'software',
  'Molecule': 'molecules',
  'BlockedBinaryMatrix': 'dataStore',
  
  }

confirmMappingClasses = {
  
  # The from->to mapping for these classes needs to be confirmed
  
  'Molecule': ['name','name',None],
  'MolSystem': ['code','code',None],
  'Entry': ['name','name',None],
  'ShiftList': ['name','serial','MeasurementList']                   
                         
}

linkDataTypes = (type(frozenset()),type(set()),type(tuple()))

# Classes autocreated by the API
autoCreatedClasses = ('ExpDim',)

# Ignore these (className,attribute) combinations when creating a list of key names to compare objects
ignoreForObjectKeys = (('PeakDim','numAliasing'),)

###################
# Copy class defs #
###################

def addClassToCreatedObjectDict(createdObjectDict,className):

  createdObjectDict[className] = {'fromObjects': [], 'toObjects': []}

def addToCreatedObjectDict(createdObjectDictForClass,fromObject,toObject):
  
  createdObjectDictForClass['fromObjects'].append(fromObject)
  createdObjectDictForClass['toObjects'].append(toObject)

class CopyObjectData:

  def __init__(self,ccpnObject,toProject,dialogInfo,toLinkedObjectData,createdObjectDict = None,ignorePresets = False,ignoreAttrs = None,extendedCopy = True, verbose = False, postponedObjects = None):
   
    if ccpnObject:
      
      #
      # Only copy if something is passed in... can be None for empty links.
      #
      
      self.ccpnObject = ccpnObject
      
      self.objectClassName = getClassOrSuperClassName(ccpnObject)
      
      self.toProject = toProject
      
      # Make sure chemElements are set...
      setCurrentStore(self.toProject, 'ChemElementStore')

      if ccpnObject.parent.className == toProject.className:
        self.parentToObject = toProject
      else:
        self.parentToObject = None

      self.parentClassName = ccpnObject.parent.className

      self.dialogInfo = dialogInfo
      
      self.toLinkedObjectData = toLinkedObjectData

      self.createdObjectDict = createdObjectDict
      self.ignorePresets = ignorePresets
      self.extendedCopy = extendedCopy
      
      self.postponedObjects = []
      
      self.ignoreAttrs = ignoreAttrs
      
      self.verbose = verbose
      
      self.toProjectResInfo = {}

      if copyOrderDict.has_key(self.objectClassName):
        self.copyObject()
      else:
        print "No copy information available for %s... aborting." % self.objectClassName
      
  def copyObject(self):
         
    copyOrder = copyOrderDict[self.objectClassName]

    if self.createdObjectDict == None:
      self.createdObjectDict = {}
  
    self.initializeParent()
    
    self.confirmCopy()
    
    self.initializeSpecificInfo()

    self.toObject = createInProject(self.ccpnObject,self.toProject,copyOrder,self.createdObjectDict, self.toProjectResInfo, self.dialogInfo, self.toLinkedObjectData, postponedObjects = self.postponedObjects, ignorePresets = self.ignorePresets, verbose = self.verbose, ignoreAttrs = self.ignoreAttrs)
    
    self.setSpecificInfo()
    
    # Create postponed objects  
    for fromObject in self.postponedObjects:
      #print fromObject
      toObject = createInProject(fromObject, self.toProject,[],self.createdObjectDict,self.toProjectResInfo, self.dialogInfo, self.toLinkedObjectData, ignorePresets = True, doMapping = False,verbose = self.verbose)
      #print toObject.resonances

  def setSpecificInfo(self):
    
    pass
    
  def initializeSpecificInfo(self):
  
    pass
  
  def confirmCopy(self):
    # TODO NOT USED YET - need to sort this out as well?
    pass
    
  def initializeParent(self):
    
    self.setStoreInfo()

    self.setCurrentStores()
        
    if not self.createdObjectDict.has_key(self.parentClassName):
      self.setParent()
    
  def setParent(self):
    
    addClassToCreatedObjectDict(self.createdObjectDict,self.parentClassName)
    addToCreatedObjectDict(self.createdObjectDict[self.parentClassName],self.ccpnObject.parent,self.parentToObject)

  def setStoreInfo(self):
    
    self.storeClassName = None
    self.storeLinkName = None

  def setCurrentStores(self):
    
    if not self.storeClassName:
      return
    
    self.fromStore = setCurrentStore(self.ccpnObject.root, self.storeClassName, linkName = self.storeLinkName)
    self.toStore = setCurrentStore(self.toProject, self.storeClassName, linkName = self.storeLinkName)

    addClassToCreatedObjectDict(self.createdObjectDict,self.storeClassName)
    addToCreatedObjectDict(self.createdObjectDict[self.storeClassName],self.fromStore,self.toStore)
    
#
# These are all 'default' copies...
#

CopyMolSystem = CopyObjectData

class CopyNmrProjectData(CopyObjectData):

  def setStoreInfo(self):
    
    self.storeClassName = 'NmrProject'
    self.storeLinkName = None

CopyNmrConstraintStore = CopyNmrProjectData
CopyAbstractMeasurementList = CopyNmrProjectData
CopyStructureGeneration = CopyNmrProjectData

class CopyAffiliationStoreData(CopyObjectData):

  def setStoreInfo(self):
    
    self.storeClassName = 'AffiliationStore'
    self.storeLinkName = None

CopyOrganisation = CopyAffiliationStoreData

class CopyMethodStoreData(CopyObjectData):

  def setStoreInfo(self):
    
    self.storeClassName = 'MethodStore'
    self.storeLinkName = None

CopyMethod = CopySoftware = CopyMethodStoreData


class CopyCitation(CopyObjectData):

  def setStoreInfo(self):
    
    self.storeClassName = 'CitationStore'
    self.storeLinkName = None

class CopyMolecule(CopyObjectData):
  
  def initializeSpecificInfo(self):
  
    self.ignoreAttrs = ['isFinalised']
 
  def setSpecificInfo(self):
    
    # Set isFinalised at end of copy process only
    self.toObject.isFinalised = self.ccpnObject.isFinalised
      
class CopyExperiment(CopyNmrProjectData):
  
  def setSpecificInfo(self):

    if self.extendedCopy:

      for dataSource in self.ccpnObject.dataSources:

        CopyDataSource(dataSource,self.toObject.root, self.dialogInfo, self.toLinkedObjectData, createdObjectDict = self.createdObjectDict, verbose = self.verbose)

      # TODO: HAVE TO ASK if need to copy associated objects... 
      # 
      # expChainStates?
      # dipolarRelaxList, ...
      # expSeriesConditions?
      # sampleConditionSet?
      #
      # probe:        ccp.api.nmr.general.Instrument.NmrProbe
      # sample:       ccp.api.lims.Sample.Sample
      # spectrometer: ccp.api.nmr.general.Instrument.NmrSpectrometer

class CopyEntry(CopyObjectData):

  def setStoreInfo(self):
    
    self.storeClassName = 'NmrEntryStore'
    self.storeLinkName = None

  def setSpecificInfo(self):

    if self.extendedCopy:
    
      for organisation in self.ccpnObject.laboratories:
        
        CopyOrganisation(organisation,self.toProject,self.dialogInfo, self.toLinkedObjectData, createdObjectDict = self.createdObjectDict, verbose = self.verbose)

      for person in list(self.ccpnObject.authors) + self.ccpnObject.sortedContactPersons():

        CopyPerson(person,self.toProject, self.dialogInfo,self.toLinkedObjectData, createdObjectDict = self.createdObjectDict, verbose = self.verbose)

class CopyPerson(CopyAffiliationStoreData):

  def initializeSpecificInfo(self):
  
    self.ignoreAttrs = ['currentPersonInGroup']
    
    organisations = []
        
    for personInGroup in self.ccpnObject.personInGroups:
      organisation = personInGroup.group.organisation
      toOrganisation = getToProjectObject(organisation,self.createdObjectDict)
      if not toOrganisation:
        organisations.append(organisation)

    for organisation in organisations:
      CopyOrganisation(organisation,self.toProject,self.dialogInfo, self.toLinkedObjectData, createdObjectDict = self.createdObjectDict, verbose = self.verbose)

  def setSpecificInfo(self):
    
    if self.ccpnObject.currentPersonInGroup:
      toCurrentPersonInGroup = getToProjectObject( self.ccpnObject.currentPersonInGroup,self.createdObjectDict)
      if toCurrentPersonInGroup:
        self.toObject.currentPersonInGroup = toCurrentPersonInGroup
      else:
        print "  Warning: could not copy current group for Person."

class CopyDataSource(CopyNmrProjectData):

  def setParent(self):
  
    expCopyOrder = copyOrderDict[self.parentClassName]

    addClassToCreatedObjectDict(self.createdObjectDict,self.parentClassName)
    addToCreatedObjectDict(self.createdObjectDict[self.parentClassName],self.ccpnObject.parent,self.parentToObject)
    
    createInProject(self.ccpnObject.parent,self.toProject,expCopyOrder,self.createdObjectDict,self.toProjectResInfo, self.dialogInfo,self.toLinkedObjectData, setMappingOnly = True, verbose = self.verbose)
  
  def setSpecificInfo(self):
 
    if self.extendedCopy:

      (selectionList,selectionDict) = createSelection(self.ccpnObject.sortedPeakLists())

      #
      # Copy over the peakList(s) if required...
      #

      if selectionList and showYesNo("Copy peakLists","Do you also want to copy the peakList(s) associated with dataSource '%d:%s'?" % (self.ccpnObject.serial,self.ccpnObject.name)):

        interaction = self.dialogInfo[1].MultiSelectionList(

                           self.dialogInfo[0],
                           selectionList,
                           selectionDict = selectionDict,
                           title = "Project '%s': " % self.toProject.name + 'Select peakLists to copy',
                           text = "Peaklists belonging to dataSource:",
                           dismissText = 'Cancel copy'

                           )

        if interaction.isSelectedList:
          for peakList in interaction.isSelectedList:
            CopyPeakList(peakList,self.toObject.root, self.dialogInfo, self.toLinkedObjectData, createdObjectDict = self.createdObjectDict, verbose = self.verbose)

class CopyBlockedBinaryMatrix(CopyObjectData):

  def setStoreInfo(self):
    
    self.storeClassName = 'DataLocationStore'
    self.storeLinkName = None

  def initializeSpecificInfo(self):

    CopyDataUrl(self.ccpnObject.dataUrl,self.toProject, self.dialogInfo, self.toLinkedObjectData, createdObjectDict = self.createdObjectDict, verbose = self.verbose)

class CopyDataUrl(CopyObjectData):

  def setStoreInfo(self):
    
    self.storeClassName = 'DataLocationStore'
    self.storeLinkName = None

  def initializeSpecificInfo(self):

    self.ignoreAttrs = ['dataStores']
 
class CopyPeakList(CopyNmrProjectData):

  def initializeParent(self):

    self.setStoreInfo()

    self.setCurrentStores()

    #
    # Setup...
    #

    experiment = getToProjectObject(self.ccpnObject.dataSource.experiment,self.createdObjectDict)
    isCreated = False
          
    if not experiment:
  
      (selectionList,selectionDict) = createSelection(self.toProject.currentNmrProject.sortedExperiments())

      #
      # Get or create the relevant NMR experiment
      #

      if selectionList and not showYesNo("New experiment","Do you also want to create a new NMR experiment from scratch? Press No if you want to link to an existing one."):

        interaction = self.dialogInfo[1].SelectionList(

                           self.dialogInfo[0],
                           selectionList,
                           selectionDict = selectionDict,
                           title = "Project '%s': " % self.toProject.name + 'Select experiment to link peak list to',
                           text = "Existing experiments:",
                           dismissText = 'Create new'

                           )

        if interaction.isSelected:
          experiment = interaction.selection

          if experiment:
            addClassToCreatedObjectDict(self.createdObjectDict,'Experiment')
            addToCreatedObjectDict(self.createdObjectDict['Experiment'],self.ccpnObject.dataSource.experiment,experiment)

            expCopyOrder = copyOrderDict['Experiment']
            createInProject(self.ccpnObject.dataSource.experiment,self.toProject,expCopyOrder,self.createdObjectDict,self.toProjectResInfo, self.dialogInfo, self.toLinkedObjectData, setMappingOnly = True, verbose = self.verbose)

      if not experiment:
        # Copy experiment and relevant datasource only!
        CopyExperiment(self.ccpnObject.dataSource.experiment,self.toProject, self.dialogInfo, self.toLinkedObjectData, createdObjectDict = self.createdObjectDict, extendedCopy = False, verbose = self.verbose)
        toExperiment = getToProjectObject(self.ccpnObject.dataSource.experiment,self.createdObjectDict)
        CopyDataSource(self.ccpnObject.dataSource,toExperiment.root, self.dialogInfo, self.toLinkedObjectData, createdObjectDict = self.createdObjectDict, extendedCopy = False, verbose = self.verbose)
        isCreated = True

    #
    # Get or create the relevant dataSource - only if existing experiment selected
    #

    if not isCreated or (experiment and not experiment.dataSources):

      dataSource = getToProjectObject(self.ccpnObject.dataSource,self.createdObjectDict)
      
      if not dataSource:
      
        (selectionList,selectionDict) = createSelection(experiment.sortedDataSources())

        if selectionList and not showYesNo("New dataSource","Do you also want to create a new dataSource from scratch? Press No if you want to link to an existing one."):

          interaction = self.dialogInfo[1].SelectionList(

                             self.dialogInfo[0],
                             selectionList,
                             selectionDict = selectionDict,
                             title = "Experiment '%s': " % experiment.name + 'Select dataSource to link peak list to',
                             text = "Existing dataSources:",
                             dismissText = 'Create new'

                             )

          if interaction.isSelected:
            dataSource = interaction.selection

            if dataSource:
              dsCopyOrder = copyOrderDict['DataSource']
              addClassToCreatedObjectDict(self.createdObjectDict,'DataSource')
              addToCreatedObjectDict(self.createdObjectDict['DataSource'],self.ccpnObject.dataSource,dataSource)

              createInProject(self.ccpnObject.dataSource,self.toProject,dsCopyOrder,self.createdObjectDict,self.toProjectResInfo,self.dialogInfo,self.toLinkedObjectData, setMappingOnly = True, verbose = self.verbose)

        if not dataSource:
          CopyDataSource(self.ccpnObject.dataSource,experiment.root,self.dialogInfo, self.toLinkedObjectData, createdObjectDict = self.createdObjectDict, extendedCopy = False, verbose = self.verbose)


##########################
# Generic copy functions #
##########################

def createInProject(ccpnObject,toProject,copyOrder,createdObjectDict,toProjectResInfo,dialogInfo,toLinkedObjectData,postponedObjects = None,setMappingOnly = False, ignorePresets = False, verbose = False, ignoreAttrs = None, doMapping = True):
  
  #
  # This dictionary keeps track of links that have to be set at end of script
  #
  
  postponedLinkDict = {}
  
  #
  # Tracks postponed objects (with links back to object to be copied, e.g. ResonanceGroup)
  #
  
  if postponedObjects == None:
    postponedObjects = []

  #
  # Copy over the main object
  #
  
  metaInfo = getMetaInfo(ccpnObject)
  
  if not setMappingOnly:
    linkName = metaInfo[0].parentRole.otherRole.name
    newObject = None
    
    if doMapping:
      # Sometimes have to force creation of new object - for resonances for example
      # TODO: make sure that this does not create spurious resonances!
      newObject = mapObject(ccpnObject,toProject,createdObjectDict,linkName,metaInfo,dialogInfo, toLinkedObjectData, verbose = verbose)
      
    if not newObject:
      newObject = copyObject(ccpnObject,toProject,createdObjectDict,postponedLinkDict,metaInfo,toProjectResInfo,postponedObjects,dialogInfo,toLinkedObjectData,ignorePresets = ignorePresets, verbose = verbose, ignoreAttrs = ignoreAttrs)
      if verbose:
        print "  Copied during createInProject:", newObject
    elif verbose:
      print "  Mapped during createInProject:", newObject
    
  else:
    newObject = None
  
  #
  # Now do all the listed links...
  #
  
  for (linkName,doCopy) in copyOrder:

    ccpnObjectStrings = linkName.split('.')      
    ccpnObjects = [ccpnObject]

    for i in range(len(ccpnObjectStrings)):
      
      # Quick check if sorted available - use if so. Pretty hacky this!
      # TODO Ask Rasmus if can get this from metamodel!
      sortedLinkName = "sorted%s%s" % (ccpnObjectStrings[i:i+1][0][0].upper(),ccpnObjectStrings[i:i+1][0][1:])
      if not hasattr(ccpnObjects[0],sortedLinkName):
        sortedLinkName = None
    
     
      newCcpnObjects = []

      for tempCcpnObject in ccpnObjects:
        try:
          if sortedLinkName:
            newCcpnObjects.extend(getattr(tempCcpnObject,sortedLinkName)())
          else:
            # To many link (e.g. project.currentNmrProject.sortedExperiments())
            newCcpnObjects.extend(list(reduce(getattr,ccpnObjectStrings[i:i+1],tempCcpnObject)))
        except:
          try:
            # Single attribute (e.g. strucGen.nmrConstraintStore)
            newCcpnObjects.append(getattr(tempCcpnObject,ccpnObjectStrings[i]))
          except:
            print "Warning: object %s does not have attribute/link %s... ignored." % (tempCcpnObject.className,ccpnObjectStrings[i])

      ccpnObjects = newCcpnObjects[:]

    #
    # Only go ahead if anything exists!
    #

    if ccpnObjects and ccpnObjects != [None]:

      #
      # Have to keep track of className in case there's subclasses... need to reset metaInfo then
      #
      
      refClassName = None

      #
      # Now loop over the objects in the fromProject
      #

      for linkedCcpnObject in ccpnObjects:
        
        #
        # Check if same class
        #
        
        if linkedCcpnObject.className != refClassName:
          metaInfo = getMetaInfo(linkedCcpnObject)
          refClassName = linkedCcpnObject.className
        
        #
        # Copy for new objects, map for objects that were automatically created (e.g. ExpDims, ...)
        #
        
        findObject = None
        
        if doCopy < 1 or setMappingOnly:
          
          #
          # Also have situation where some stuff is autocreated, but can also be manually created... what a mess!
          #
          
          if not setMappingOnly:
            setInfo = True
          else:
            setInfo = False
          
          findObject = mapObject(linkedCcpnObject,toProject,createdObjectDict,ccpnObjectStrings[-1],metaInfo,dialogInfo,toLinkedObjectData, setInfo = setInfo,verbose = verbose)
        
        if not setMappingOnly and (doCopy == 1 or (doCopy == -1 and not findObject)):
        
          copyObject(linkedCcpnObject,toProject,createdObjectDict,postponedLinkDict,metaInfo,toProjectResInfo,postponedObjects,dialogInfo,toLinkedObjectData, ignorePresets = ignorePresets, verbose = verbose)
  
  #
  # Set all the postponed links...
  #
  
  if not setMappingOnly:
  
    for fromCcpnObject in postponedLinkDict.keys():

      toCcpnObject = getToProjectObject(fromCcpnObject,createdObjectDict)

      if not toCcpnObject:

        #
        # Try to get it the 'slow way'
        #

        linkedMetaInfo = getMetaInfo(fromCcpnObject)
        linkedLinkName = linkedMetaInfo[0].parentRole.otherRole.name
        if verbose:
          print "  Setting postponed %s object in createInProject." % linkedMetaInfo[0].name
        toCcpnObject = mapObject(fromCcpnObject,toProject,createdObjectDict,linkedLinkName,linkedMetaInfo, dialogInfo,toLinkedObjectData, verbose = verbose)
        
        if not toCcpnObject:
          continue

      linkNames = postponedLinkDict[fromCcpnObject].keys()
      linkNames.sort()

      for linkName in linkNames:

        (hicard,fromLinkValue) = postponedLinkDict[fromCcpnObject][linkName]

        currentToLinkValue = getattr(toCcpnObject,linkName)

        if hicard == 1:

          toLinkValue = getToProjectObject(fromLinkValue,createdObjectDict)

          if toLinkValue:
            if currentToLinkValue != toLinkValue:
              setattr(toCcpnObject,linkName,toLinkValue)
              if verbose:
                print "  Link set in createInProject:", linkName, toLinkValue
            del(postponedLinkDict[fromCcpnObject][linkName])

        else:

          nonMatchingLinkValues = []
          toLinkValue = []
          currentToLinkValue = list(currentToLinkValue)

          for fromLinkItem in fromLinkValue:
            toLinkItem = getToProjectObject(fromLinkItem,createdObjectDict)
            if toLinkItem and toLinkItem not in toLinkValue:
              toLinkValue.append(toLinkItem)
              if toLinkItem not in currentToLinkValue:
                nonMatchingLinkValues.append(toLinkValue)

          if toLinkValue:
            if nonMatchingLinkValues:
              #print linkName
              #print currentToLinkValue
              #print toLinkValue
              #print tempToLinkValue
              # TODO: add items from nonMatchingLinkValues if there is a match!?!
              setattr(toCcpnObject,linkName,toLinkValue)
              if verbose:
                print linkName, toLinkValue
            del(postponedLinkDict[fromCcpnObject][linkName])

      if not postponedLinkDict[fromCcpnObject]:
        del(postponedLinkDict[fromCcpnObject])

  if verbose:
    print "  Postponed during createInProject: %s" % str(postponedLinkDict)

  return newObject
            
def getMetaInfo(ccpnObject):

  #
  # Get meta information...
  #

  metaClass = ccpnObject._metaclass
  metaInfo = (metaClass,[],[])

  # getAllAttributes() includes those from superClasses
  #print metaClass.name
  for metaAttr in metaClass.getAllAttributes():

    if not metaAttr.isDerived and not metaAttr.isAutomatic and not metaAttr.isImplementation:
      metaInfo[1].append(metaAttr)
      #print metaAttr.name, metaAttr.isImplementation

  for metaRole in metaClass.getAllRoles():

    if not metaRole.isDerived and not metaRole.hierarchy == 'child' and not metaRole.hierarchy == 'parent' and not metaRole.isImplementation:
      metaInfo[2].append(metaRole)
  
  return metaInfo
      
def setToProjectResInfo(toProjectResonances,toProjectResInfo,createdObjectDict):

  toResonanceAtoms = {}
  resNamesList = []

  for toResonance in toProjectResonances:

    toResonanceMapped = False

    if toResonance in createdObjectDict[toResonance.className]['toObjects']:
      toResonanceMapped = True
      value = getResonanceResName(toResonance)
      if value:
        resNamesList.append(value)

      break

    if not toResonanceMapped:
      
      setResonanceAtoms(toResonance,toResonanceAtoms)
        
  toProjectResInfo['toResonanceAtoms'] = toResonanceAtoms
  toProjectResInfo['resNamesList'] = resNamesList

def setResonanceAtoms(resonance,resonanceAtoms):

  resonanceSet = resonance.resonanceSet
  if resonanceSet:
    resonanceAtoms[resonance] = []
    for atomSet in resonanceSet.atomSets:
      atoms = list(atomSet.atoms)
      atoms.sort()
      resonanceAtoms[resonance].append(tuple(atoms))
    resonanceAtoms[resonance].sort()

def getResonanceResName(resonance):
  
  value = None
  appData = resonance.findFirstApplicationData(application = 'ccpn', keyword = assign_kw)
  if appData:
    value = appData.value
  return value

def setPresetLinks(metaRole,otherObjects,toProject,keywds,createdObjectDict,toProjectResInfo,postponedObjects,dialogInfo,toLinkedObjectData,ignorePresets = False, verbose = False):

  newObjects = []
  
  if otherObjects and not (len(otherObjects) == 1 and otherObjects[0] == None):

    refObject = otherObjects[0]
    
    if refObject.className[-8:] == 'ChemComp':

      newObjects.append(toProject.findFirstChemComp(molType = refObject.molType,ccpCode = refObject.ccpCode))
        
    #elif refObject.className == 'DataLocation':
    
    #  url = getUrl(toProject, name='dataLocation', path=refObject.url.path)
    #  newObjects.append(Implementation.DataLocation(toProject, path=refObject.path, url=url))
    
      """  
    elif refObject.className == 'ContentStorage':
    
      fromStoragePath = refObject.url.path
      fromProjectPath = refObject.root.url.path
      toProjectPath = toProject.url.path
      
      newPath = fromStoragePath.replace(fromProjectPath,toProjectPath)
      
      if fromProjectPath == toProjectPath or newPath != fromStoragePath:      
        url = getUrl(toProject, refObject.url.name, newPath)
        newObjects.append(Implementation.ContentStorage(toProject, package = refObject.package, path=refObject.path, url=url))
      """  
  
    elif refObject.className in classLinkNames:
    
      for refObject in otherObjects:

        if refObject.className == 'Method':

          #
          # For Method also have to set software... if not set already!!
          #
          
          setPresetLinks(refObject._metaclass.getElement('software'),[refObject.software],toProject,{},createdObjectDict,toProjectResInfo,postponedObjects, dialogInfo, toLinkedObjectData, ignorePresets = True, verbose = verbose)
          forceAttrOrRoles = ['software']

        else:
          forceAttrOrRoles = []

        parentLinkName = classLinkNames[refObject.className]

        metaInfo = getMetaInfo(refObject)

        if createdObjectDict == None:
          createdObjectDict = {}
        if not createdObjectDict.has_key('Project'):
          addClassToCreatedObjectDict(createdObjectDict,'Project')
          addToCreatedObjectDict(createdObjectDict['Project'],refObject.root,toProject)

        newObject = mapObject(refObject,toProject,createdObjectDict,parentLinkName,metaInfo,dialogInfo,toLinkedObjectData, forceAttrOrRoles = forceAttrOrRoles,verbose = verbose)

        if newObject:
          newObjects.append(newObject)

        else:

          CopyClass = globals()["Copy%s" % refObject.className]
          # Using dialogInfo here, but probably will not give right results when asking questions...
          CopyClass(refObject,toProject,dialogInfo,toLinkedObjectData,createdObjectDict = createdObjectDict,ignorePresets = ignorePresets, verbose = verbose)

          if refObject in createdObjectDict[refObject.className]['fromObjects']:
            refIndex = createdObjectDict[refObject.className]['fromObjects'].index(refObject)
            newObjects.append(createdObjectDict[refObject.className]['toObjects'][refIndex])
            break

        #if refObject.className == 'Method':
        #  print refObject.task, newObjects[0].task
        #  print refObject.procedure, newObjects[0].procedure
        #  print refObject.software, newObjects[0].software
        
    elif refObject.className in ('Resonance','FixedResonance'):
      
      resClassName = refObject.className
      
      #
      # Setup
      #
      
      if createdObjectDict == None:
        createdObjectDict = {}
      if not createdObjectDict.has_key('Project'):
        addClassToCreatedObjectDict(createdObjectDict,'Project')
        addToCreatedObjectDict(createdObjectDict['Project'],refObject.root,toProject)
      if not createdObjectDict.has_key(refObject.parent.className):
        # This is for NmrConstraintStore in case of FixedResonances...
        addClassToCreatedObjectDict(createdObjectDict,refObject.parent.className)
        #addToCreatedObjectDict(createdObjectDict[refObject.parent.className],refObject.parent,toProject)
      if not createdObjectDict.has_key(resClassName):
        addClassToCreatedObjectDict(createdObjectDict,resClassName)
      
      #
      # Check if mapping exists already...
      #
      
      checkResonances = []
        
      for refObject in otherObjects:
        
        resonanceFound = False
        
        if refObject in createdObjectDict[resClassName]['fromObjects']:
          refIndex = createdObjectDict[resClassName]['fromObjects'].index(refObject)
          newObjects.append(createdObjectDict[resClassName]['toObjects'][refIndex])
          if verbose:
            print "  Direct %s map: %s -> %s" % (resClassName,ccpnObject,newObject)
          resonanceFound = True
          break
        
        if not resonanceFound:
          checkResonances.append(refObject)
          newObjects.append(None) # Have to stick to original order...
      
      #
      # If not mapped, try to set the mapping...
      #
      
      if checkResonances:
      
        #
        # Keep track of molecules, molSystems and chains...
        #

        linkedObjects = [[],[],[],[]]
        parentLinkNames = ('molecules','molSystems','chains','residues')
        # If change above, also set linkedInfo below correctly!!

        fromObjectDict = {'Molecule': {}, 'MolSystem': {}, 'Chain': {}, 'Residue': {}}

        for resonance in checkResonances:

          resonanceSet = resonance.resonanceSet
          if resonanceSet:
            for atomSet in resonanceSet.atomSets:
              refAtom = atomSet.findFirstAtom()
              residue = refAtom.residue
              chain = residue.chain
              linkedInfo = (chain.molecule,chain.molSystem,chain,residue)
              for i in range(4):
                linkedObject = linkedInfo[i]
                
                # Check if not already there!
                objectFound = False
               
                if linkedObject and createdObjectDict.has_key(linkedObject.className):
                  if linkedObject in createdObjectDict[linkedObject.className]['fromObjects']:
                    refIndex = createdObjectDict[linkedObject.className]['fromObjects'].index(linkedObject)
                    toObject = createdObjectDict[linkedObject.className]['toObjects'][refIndex]
                    fromObjectDict[linkedObject.className][linkedObject] = toObject
                    objectFound = True
                
                if not objectFound and linkedObject not in linkedObjects[i]:
                  linkedObjects[i].append(linkedObject)
        
        #
        # Make sure unmapped molecules, molSystems and chains are mapped...
        #
        
        for i in range(4):
          if linkedObjects[i]:
            linkName = parentLinkNames[i]
            metaInfo = getMetaInfo(linkedObjects[i][0])
            for linkedObject in linkedObjects[i]:
              if linkedObject:
                toObject = mapObject(linkedObject,toProject,createdObjectDict,linkName,metaInfo,dialogInfo,toLinkedObjectData,verbose = verbose)  
                if toObject:
                  fromObjectDict[linkedObject.className][linkedObject] = toObject      

        #
        # Make sure a list of resonances and mapped info is available (only done once)
        #
        # AtomSet ordering is also taken into account
        #
        
        if not toProjectResInfo:
          if resClassName == 'Resonance':
            fullReslist = []
            if toProject.currentNmrProject:
              fullReslist = toProject.currentNmrProject.sortedResonances()
          else:
            toNmrConstraintStore = getToProjectObject(refObject.parent,createdObjectDict)
            fullReslist = toNmrConstraintStore.sortedFixedResonances()
            
          setToProjectResInfo(fullReslist,toProjectResInfo,createdObjectDict)
        
        toResonanceAtoms = toProjectResInfo['toResonanceAtoms']
        resNamesList = toProjectResInfo['resNamesList']
  
        #
        # Now check whether the resonance can be mapped unambiguously
        #
        
        addClassToCreatedObjectDict(createdObjectDict,'Atom')
                
        for resonance in checkResonances:
          
          resonanceSet = resonance.resonanceSet
          toResonance = None
          toAtomList = []
          
          if resonanceSet:
            for atomSet in resonanceSet.atomSets:
              tempAtomList = []
              for atom in atomSet.atoms:
                toAtom = None
                fromResidue = atom.residue
                if fromObjectDict['Residue'].has_key(fromResidue):
                  for tempToAtom in fromObjectDict['Residue'][fromResidue].atoms:
                    if tempToAtom.name == atom.name:
                      toAtom = tempToAtom
                      addToCreatedObjectDict(createdObjectDict['Atom'],atom,toAtom)
                      break
                if toAtom:
                  tempAtomList.append(toAtom)
              tempAtomList.sort()
              if tempAtomList:
                toAtomList.append(tuple(tempAtomList))
            toAtomList.sort()
           
            for tempToResonance in toResonanceAtoms.keys():
              if toAtomList == toResonanceAtoms[tempToResonance]:
                toResonance = tempToResonance
                # Only map if not mapped already!
                if toResonance not in createdObjectDict[resonance.className]['toObjects']:
                  if verbose:
                    print "  Atom based %s map: %s -> %s" % (resonance,resClassName,toResonance)
                    print "    ",toAtomList
                  addToCreatedObjectDict(createdObjectDict[resonance.className],resonance,toResonance)
                  break
          
          #
          # Create a new resonance if no direct match found
          #
          
          if not toResonance:
            
            if resClassName == 'Resonance' and resonance.resonanceGroup:
              resGroup = resonance.resonanceGroup
              toResGroup = getToProjectObject(resGroup,createdObjectDict)

              if not toResGroup and resGroup not in postponedObjects:
                postponedObjects.append(resGroup)
                
            toResonance = createInProject(resonance,toProject,[],createdObjectDict,toProjectResInfo, dialogInfo, toLinkedObjectData, verbose = verbose, doMapping = False)            

            if toAtomList:

              #
              # No match, but have list of atoms - link to the corresponding ones
              # This is the case where there already is a matching molecular system present!!!
              #
              
              for atomSet in resonanceSet.atomSets:
                
                toAtomSet = createInProject(atomSet,toProject,[],createdObjectDict,toProjectResInfo, dialogInfo, toLinkedObjectData, verbose = verbose)
              
              toResonanceSet = createInProject(resonanceSet,toProject,[],createdObjectDict,toProjectResInfo, dialogInfo, toLinkedObjectData, verbose = verbose)             
                            
            else:

              #
              # Set mapping info for linkResonances... now also works for ambiguous stuff...
              # This is the case where there is no molecular system (or no matching one)
              #

              if resonanceSet:
                
                resNames = []

                for atomSet in resonanceSet.atomSets:
                  atomName = atomSet.name
                  residue = atomSet.findFirstAtom().residue
                  chainCode = residue.chain.code

                  resNames.append(getResName(chainCode,residue.seqCode,atomName,seqInsertCode = residue.seqInsertCode))
                  
                # Specific Analysis case...
                # TODO need to put these strings in ccpnmr.general.Constants or something! Shared...
                analysisAppData = resonance.findFirstApplicationData(application ='Analysis',keyword='ambigProchiralLabel')
                
                if analysisAppData and len(resNames) == 2:
                  # TODO REMOVE PRINT!
                  #print "CHECK ON DATASHIFTER",resNames,analysisAppData.value
                  if analysisAppData.value == 'a':
                    resNames.pop(1)
                  else:
                    resNames.pop(0)
                    
                for resName in resNames:

                  if resName in resNamesList:
                    continue

                  # Try to avoid copying info to same resonance in case of ambiguity
                  if not toResonance.findFirstApplicationData(application = 'ccpn', keyword = assign_kw):

                    resNamesList.append(resName)
                    toResonance.addApplicationData(Implementation.AppDataString(application = 'ccpn', keyword = assign_kw, value = resName))
                    break

                if not toResonance.parent.findFirstApplicationData(keyword = newResonances_kw):
                  toResonance.parent.addApplicationData(Implementation.AppDataBoolean(application = 'ccpn', keyword = newResonances_kw, value = True))
          
          newObjects[newObjects.index(None)] = toResonance
          
          #print "TPRESINFO", toProjectResInfo
        
        #print checkResonances
        #print newObjects
          
    #
    # Make sure the link gets set
    #
     
    if newObjects:
      if metaRole.hicard == 1:
        keywds[metaRole.name] = newObjects[0]
      else:
        keywds[metaRole.name] = newObjects
  
  return newObjects

def mapObject(ccpnObject,toProject,createdObjectDict,linkName,metaInfo,dialogInfo,toLinkedObjectData,setInfo = False,forceAttrOrRoles = None, verbose = False, parentSearch = True):
  
  debug = False
  #if ccpnObject.className in ('ResonanceSet','AtomSet'):
  #  debug = True
  
  if debug:
    print " MO",setInfo, forceAttrOrRoles
  
  objectMapped = None
  
  #
  # If objects automatically created in new project, have to map them
  #
  
  metaClass = metaInfo[0]
  keyNames = metaClass.keyNames
  
  if debug:
    print metaClass.name
  
  #
  # Check if keyNames fit with metaInfo that was passed in (root out serials and such)
  #
  
  validNames = []
  
  for metaAttrOrRole in metaInfo[1] + metaInfo[2]:
    validNames.append(metaAttrOrRole.name)
    # Also include obligatory stuff as part of key! Except if class was set automatically - then have to transfer value later on. 
    if metaClass.name not in autoCreatedClasses and metaAttrOrRole.locard >= 1 and metaAttrOrRole.name not in keyNames:
      #print metaClass.name, metaAttrOrRole.name, ignoreForObjectKeys
      if (metaClass.name,metaAttrOrRole.name) not in ignoreForObjectKeys:
        keyNames.append(metaAttrOrRole.name)
      
  for keyName in keyNames[:]:
    if keyName not in validNames:
      keyNames.pop(keyNames.index(keyName))
      
  if not keyNames:
    for metaAttr in metaInfo[1]:
      keyNames.append(metaAttr.name)

  if forceAttrOrRoles:
    keyNames.extend(forceAttrOrRoles)
  
  #
  # Get info from fromProject object
  #
  
  fromData = []

  for keyName in keyNames:
  
    linkedData = getattr(ccpnObject,keyName)
    
    if type(linkedData) in linkDataTypes:
      linkedData = list(linkedData)
    
    if type(linkedData) == type([]) and linkedData and hasattr(linkedData[0],'parent'):
      linkedData.sort()
    
    if not (keyName,linkedData) in fromData:
      fromData.append((keyName,linkedData))
  
  if debug:
    print keyNames
    print fromData
  
  #
  # Get the relevant parent object in toProject
  #
  
  fromParent = ccpnObject.parent
  toParent = getToProjectObject(fromParent,createdObjectDict)
  
  if debug:
    print "FROM PARENT", fromParent
    print "FOUND TO PARENT", toParent
  
  if not toParent and parentSearch:
  
    #
    # Try to get it the 'slow way'
    #
    
    parentMetaInfo = getMetaInfo(fromParent)
    
    if parentMetaInfo[0].parentRole:
      # TODO TODO: is this supposed to be allowed to fail?!? The .parentRole?
      parentLinkName = parentMetaInfo[0].parentRole.otherRole.name
      if verbose or debug:
        print "  mapObjects trying to find %s" % parentMetaInfo[0].name
      toParent = mapObject(ccpnObject.parent,toProject,createdObjectDict,parentLinkName,parentMetaInfo,dialogInfo, toLinkedObjectData,verbose = False, parentSearch = False)

  if toParent:
    
    if debug:
      print " TO PARENT", toParent
  
    #
    # Now compare...
    #
    
    if not hasattr(toParent,linkName):
      linkName = metaInfo[0].parentRole.otherRole.name
      
    if debug:
      print "  LINKNAME", linkName

    for toLinkedObject in getattr(toParent,linkName):
    
      doBreak = False
      
      if debug:
        print "   TO LINKED OBJECT", toLinkedObject
      
      # TODO: this is slightly dangerous, but shouldn't be changing info of existing objects, so is fine?
      if toLinkedObjectData.has_key((toParent,linkName,toLinkedObject)):
        toData = toLinkedObjectData[(toParent,linkName,toLinkedObject)]
        if debug:
          print "     DIRECT TO DATA", toData
      
      else:
        toData = []
        if debug:
          print "     INDIRECT TO DATA"

        for keyName in keyNames:

          linkedData = getattr(toLinkedObject,keyName)
          
          if debug:
             print "          LINKED DATA",keyName,linkedData
          
          if type(linkedData) in linkDataTypes:
            linkedData = list(linkedData)

          if hasattr(linkedData,'parent'):

            linkedData = getFromProjectObject(linkedData,createdObjectDict)

          elif type(linkedData) == type([]) and linkedData and hasattr(linkedData[0],'parent'):

            for i in range(len(linkedData)):
              linkedData[i] = getFromProjectObject(linkedData[i],createdObjectDict)

            linkedData.sort()
            
          if not (keyName,linkedData) in toData:
            toData.append((keyName,linkedData))
        
        toLinkedObjectData[(toParent,linkName,toLinkedObject)] = toData
        #print "SET", (toParent,linkName,toLinkedObject), toData
        #if doBreak:
        #  sys.exit()

      if toData == fromData:
        
        #
        # For specific high-level classes, make sure that mapping is OK with user...
        #
        
        doMapping = True
        
        if ccpnObject.className in confirmMappingClasses.keys():
          className = ccpnObject.className
          (attrName,keyName,superClass) = confirmMappingClasses[className]
          fromAttrValue = getattr(ccpnObject,attrName)
          toAttrValue = getattr(toLinkedObject,attrName)
          if not dialogInfo[2].showYesNo("Mapping to existing object","%s object with %s '%s' is being mapped to '%s' in the other project. Is this mapping correct?" % (className,attrName,str(fromAttrValue),str(toAttrValue))):
            doMapping = False
        
        if debug:
          print "  DOING MAPPING: %s" % doMapping
        
        if doMapping:
          if not createdObjectDict.has_key(metaClass.name):
            addClassToCreatedObjectDict(createdObjectDict,metaClass.name)
  
          addToCreatedObjectDict(createdObjectDict[metaClass.name],ccpnObject,toLinkedObject)
  
          objectMapped = toLinkedObject
  
          break

    #
    # Now also have to transfer information between objects if required... it's a mess!!
    #
    
    if debug:
        print "  OBJECT MAPPED %s and SETINFO %s" % (str(objectMapped),setInfo)

    if objectMapped and setInfo:

      for metaAttrOrRole in metaInfo[1] + metaInfo[2]:

        if metaAttrOrRole.changeability == 'changeable':

          attrOrRoleName = metaAttrOrRole.name

          if attrOrRoleName not in keyNames:
            value = getattr(ccpnObject,attrOrRoleName)
            
            # Convert (frozen)sets to lists
            if type(value) in linkDataTypes:
              value = list(value)
              
            #print ccpnObject, attrOrRoleName, value, type(value), linkDataTypes 

            if value or value == 0.0 or value == False:
              existingValue = getattr(objectMapped,attrOrRoleName)

              if not existingValue:

                #
                # Get relevant objects from toProject if required...
                #

                if hasattr(value,'parent'):

                  value = getToProjectObject(value,createdObjectDict)

                elif type(value) == type([]) and value and hasattr(value[0],'parent'):

                  for i in range(len(value)):
                    value[i] = getToProjectObject(value[i],createdObjectDict)

                  if value.count(None) == len(value):
                    value = None
                
                elif type(value) not in (type(0.0),type(""),type(1),type(False)):
                  
                  print "  Error: unknown original data type %s - ignored." % type(value)
                  value = None

                #
                # Set the information
                #

                if value != None:
                  
                  if debug:
                    print "TRYING TO SET LINK",objectMapped,attrOrRoleName,value

                  try:
                    setattr(objectMapped,attrOrRoleName,value)
                    if verbose:
                      print "  mapObjects set %s.%s to '%s'" % (objectMapped.className,attrOrRoleName,str(value))
                  except:
                    print "  Error: mapObjects cannot set %s for %s... API error '%s'" % (attrOrRoleName,objectMapped.className,traceback.format_exception_only(sys.exc_type,sys.exc_value))

              else:
                print "  Warning: mapObjects NOT setting %s.%s, is '%s'" % (objectMapped.className,attrOrRoleName,str(existingValue))
    
  #
  # Return mapped object...
  #
  
  if verbose:
    if objectMapped:
      print "  mapObjects mapped:", objectMapped
    else:
      print "  WARNING: mapObjects could not map %s from original object %s." % (metaClass.name,ccpnObject)
      #print keyNames
  
  return objectMapped
      
def getToProjectObject(fromObject,createdObjectDict):

  metaClassName = fromObject._metaclass.name
  toObject = None
  
  if createdObjectDict.has_key(metaClassName):
    if fromObject in createdObjectDict[metaClassName]['fromObjects']:
      refIndex = createdObjectDict[metaClassName]['fromObjects'].index(fromObject)
      toObject = createdObjectDict[metaClassName]['toObjects'][refIndex]
      
  return toObject

def getFromProjectObject(toObject,createdObjectDict):

  metaClassName = toObject._metaclass.name
  fromObject = None
  
  if createdObjectDict.has_key(metaClassName):
    if toObject in createdObjectDict[metaClassName]['toObjects']:
      refIndex = createdObjectDict[metaClassName]['toObjects'].index(toObject)
      fromObject = createdObjectDict[metaClassName]['fromObjects'][refIndex]
      
  return fromObject

def getClassOrSuperClassName(ccpnObject):
        
  ccpnObjectMetaClass = ccpnObject._metaclass
  superType = ccpnObjectMetaClass.getSupertype()

  if superType and superType.name not in ('TopObject','DataObject'):
    objectClassName = superType.name
  else:
    objectClassName = ccpnObject.className
    
  return objectClassName

def copyObject(ccpnObject,toProject,createdObjectDict,postponedLinkDict,metaInfo,toProjectResInfo,postponedObjects,dialogInfo,toLinkedObjectData,ignorePresets = False, verbose = False, ignoreAttrs = None):

  metaClass = metaInfo[0]
  
  keywds = {}
  
  if not ignoreAttrs:
    ignoreAttrs = []
  
  for metaAttr in metaInfo[1]:
    
    attrName = metaAttr.name
    
    if attrName not in ignoreAttrs:
      keywds[attrName] = getattr(ccpnObject,attrName)
  
  for metaRole in metaInfo[2]:
    
    roleName = metaRole.name
    otherMetaClassName = None

    if roleName not in ignoreAttrs:
      
      otherMetaClassName = metaRole.valueType.name

      if metaRole.hicard > 1 or metaRole.hicard == -1:
        otherObjects = getattr(ccpnObject,roleName)        
      else:
        otherObjects = [getattr(ccpnObject,roleName)]

    #
    # See if linked object already exists in toProject
    #
    
    validToProjectObjects = [] 

    if metaClass.name != otherMetaClassName: # Have to postpone links to same class.
      for otherObject in otherObjects:
        if createdObjectDict.has_key(otherMetaClassName):
          if otherObject in createdObjectDict[otherMetaClassName]['fromObjects']:
            refIndex = createdObjectDict[otherMetaClassName]['fromObjects'].index(otherObject)
            toProjectObject = createdObjectDict[otherMetaClassName]['toObjects'][refIndex]
            validToProjectObjects.append(toProjectObject)
          else:
            validToProjectObjects.append(None)
        else:
          validToProjectObjects.append(None)

    #
    # If nothing, use hardcoded links to reference data
    #

    if not ignorePresets and None in validToProjectObjects:
      otherObjects = list(otherObjects)

      for otherIndex in range(len(otherObjects)):
        if not validToProjectObjects[otherIndex]:
        
          otherObject = otherObjects[otherIndex]
          newObjects = setPresetLinks(metaRole,[otherObject],toProject,keywds,createdObjectDict,toProjectResInfo,postponedObjects, dialogInfo, toLinkedObjectData, verbose = verbose)

          if newObjects:
            if verbose:
              print "  copyObject; original:'%s' to new:'%s'." % (str(otherObject),str(newObjects))
            if not createdObjectDict.has_key(otherMetaClassName):
              addClassToCreatedObjectDict(createdObjectDict,otherMetaClassName)

            if not otherObject in createdObjectDict[otherMetaClassName]['fromObjects']:
              addToCreatedObjectDict(createdObjectDict[otherMetaClassName],otherObject,newObjects[0])
              
            validToProjectObjects[otherIndex] = newObjects[0]

            continue
     
    #
    # Set the info for the link...
    #
    
    for projObj in validToProjectObjects[:]:
      if projObj == None:
        validToProjectObjects.pop(validToProjectObjects.index(projObj))

    if validToProjectObjects:
      if metaRole.hicard == 1:
        keywds[roleName] = validToProjectObjects[0]
      else:
        keywds[roleName] = validToProjectObjects

    #
    # If not, postpone setting this link till later...
    #

    if not keywds.has_key(roleName):

      linkValue = getattr(ccpnObject,roleName)
      
      if linkValue:
      
        if not postponedLinkDict.has_key(ccpnObject):
          postponedLinkDict[ccpnObject] = {}

        if not postponedLinkDict[ccpnObject].has_key(roleName):
          postponedLinkDict[ccpnObject][roleName] = None
          
        postponedLinkDict[ccpnObject][roleName] = (metaRole.hicard,linkValue)
                  
  #
  # First check if object wasn't created during above process... is possible!
  #
  
  newObject = getToProjectObject(ccpnObject,createdObjectDict)
          
  #
  # If it doesn't exist yet, create it...
  #
  
  if not newObject:
  
    createClass = ccpnObject.__class__  

    objectParent = ccpnObject.parent
    newParent = getToProjectObject(objectParent,createdObjectDict)
    
    if newParent:
      try:
        if ccpnObject.className in confirmMappingClasses.keys():
          className = ccpnObject.className
          (attrName,keyName,superClassName) = confirmMappingClasses[className]
          attrValue = getattr(ccpnObject,attrName)
          keyValue = getattr(ccpnObject,keyName)
          
          searchKwd = {keyName: keyValue}
          
          searchClass = className
          if superClassName:
            searchClass = superClassName
          
          while (getattr(newParent,'findFirst%s' % searchClass)(**searchKwd)):
            
            attrValue = dialogInfo[3].askString("Enter %s %s" % (className,attrName),"Enter a new %s for the copy of the %s in the other project:" % (attrName,className),attrValue,dialogInfo[0])
            searchKwd = {attrName: attrValue}
            
          keywds[attrName] = attrValue
      
        newObject = createClass(newParent,**keywds)

        if verbose and newObject:
          print "  copyObjects created copy:", newObject

      except:
        # Was it missed in the mapping? Try this...
        linkName = metaInfo[0].parentRole.otherRole.name
        newObject = mapObject(ccpnObject,toProject,createdObjectDict,linkName,metaInfo, dialogInfo, toLinkedObjectData, verbose = verbose)
        
        if not newObject:
          print "  ERROR: copyObjects could not create or map new %s..." % ccpnObject.className
          print newParent
          print keywds
          print ccpnObject
          raise

      if not createdObjectDict.has_key(metaClass.name):
        addClassToCreatedObjectDict(createdObjectDict,metaClass.name)

      addToCreatedObjectDict(createdObjectDict[metaClass.name],ccpnObject,newObject)

    else:
      print "  ERROR: copyObjects could not create new %s... parent missing in other project." % ccpnObject.className
      print objectParent
    
  return newObject

############
# GUI code #
############

class DataShifter(BasePopup):

  """
  
  GUI to copy 'data entities' between different CCPN projects.
  
  """
  
  help_url = joinPath(getHelpUrlDir(),'DataShifter.html')
  applicationName = 'DataShifter'
  
  def __init__(self, parent, project1 = None, project2 = None, testMode = False, verbose = False):
  
    self.guiParent = parent
    self.project1 = project1
    self.project2 = project2
    self.verbose = verbose
    self.testMode = testMode
    
    self.multiDialog = setupMultiDialog(self.guiParent)
    self.messageReporter = setupMessageReporter(self.guiParent)
    self.dataEntry = setupDataEntry(self.guiParent)

    #
    # Below are defs that will autocreate the window
    #

    self.setCurrentStores = ['NmrProject','MethodStore','AffiliationStore','CitationStore','NmrEntryStore','DataLocationStore']
    self.widgetDefs = [
    
      #['FileButton','fileName','File for project:','Select project file',None],
      ['SelectionListWithNone','molSystem','Select a molecular system to copy:','project.sortedMolSystems()',None],
      ['SelectionListWithNone','molecule','Select a molecule to copy:','project.sortedMolecules()',None],
      ['Separator',None,None,None,None],
      ['SelectionListWithNone','experiment','Select NMR experiment to copy:',"project.currentNmrProject.sortedExperiments()",None],
      ['SelectionListWithNone','peakList','Select peak list to copy:',"project.currentNmrProject.sortedExperiments().sortedDataSources().sortedPeakLists()",None],
      ['SelectionListWithNone','measurementList','Select measurement list to copy:',"project.currentNmrProject.sortedMeasurementLists()",None],
      ['Separator',None,None,None,None],
      ['SelectionListWithNone','Structure generation','Select Structure Generation to copy:','project.currentNmrProject.sortedStructureGenerations()',None],      
      ['SelectionListWithNone','Constraint set','Select NmrConstraintStore to copy:','project.currentNmrProject.sortedNmrConstraintStores()',None],
      ['Separator',None,None,None,None],
      #['SelectionListWithNone','method','Select a method to copy:','project.currentMethodStore.methods',None],
      ['SelectionListWithNone','software','Select software to copy:','project.currentMethodStore.sortedSoftware()',None],
      ['SelectionListWithNone','organisations','Select organisation to copy:','project.currentAffiliationStore.sortedOrganisations()',None],
      ['SelectionListWithNone','persons','Select person to copy:','project.currentAffiliationStore.sortedPersons()',None],
      ['SelectionListWithNone','citations','Select citation to copy:','project.currentCitationStore.sortedCitations()',None],
      ['Separator',None,None,None,None],
      ['SelectionListWithNone','entries','Select entries to copy:','project.currentNmrEntryStore.sortedEntries()',None],
     
      ]

    BasePopup.__init__(self,parent = parent, title = "Shift data between projects", modal = False, transient=False)
      
  def body(self, master):
    
    #
    # Set project stuff
    #

    self.widgets = {}
    self.widgetInfo = {}
    
    widgetRowList = []
    
    for (projectNum,colNum) in ((1,0),(2,3)):
     
      projectId = "project%d" % projectNum
      project = getattr(self,projectId)
      
      if project:
        for currentStoreName in self.setCurrentStores:
          setCurrentStore(project,currentStoreName,createNew = False)
        
      self.widgets[projectId] = {}
      self.widgetInfo[projectId] = {}
          
      row = 0
      
      label = Label(master, text = "Details for project %d" % projectNum)
      label.grid(row=row, column=colNum, sticky=Tkinter.EW)

      row += 1
      
      separator = Separator(master,height = 3)
      separator.setColor('black', bgColor = 'black')
      separator.grid(row=row, column=colNum, sticky=Tkinter.EW)

      row += 1
      self.widgetStartRow = row
      
      #
      # Project load
      #
      
      buttonKeyword = 'directoryName'
      label = Label(master, text = "CCPN project")
      widget = Tkinter.Button(master, text = "Project directory", command = lambda projectId = projectId: self.selectProjectDir(projectId=projectId))
      
      self.widgets[projectId][buttonKeyword + "Label"] = label
      self.widgets[projectId][buttonKeyword + "Label"].grid(row=row, column=colNum, sticky=Tkinter.EW)

      row += 1
        
      widgetRowList.append(row)
      self.widgets[projectId][buttonKeyword] = widget
      self.widgetInfo[projectId][buttonKeyword] = {}

      self.widgets[projectId][buttonKeyword].grid(row=row, column=colNum, sticky=Tkinter.EW)

      row += 1
      
      #
      # Autoset the widgets on both sides
      #
    
      for widgetIndex in range(len(self.widgetDefs)):

        widgetInfo = self.widgetDefs[widgetIndex]
        
        if widgetInfo[0] == 'Separator':
          separator = Separator(master,height = 1)
          separator.setColor('black', bgColor = 'black')
          separator.grid(row=row, column=colNum, sticky=Tkinter.EW)
          row += 1
          continue
        
        buttonKeyword = widgetInfo[1]

        (label,widget,selectionDict) = setIoInfo(master,widgetInfo,project,(projectId,'Project',self.selectFile,None), forceCreation = True, appendNone = True)

        self.widgets[projectId][buttonKeyword + "Label"] = label
        self.widgets[projectId][buttonKeyword + "Label"].grid(row=row, column=colNum, sticky=Tkinter.EW)

        row += 1
        
        widgetRowList.append(row)
        self.widgets[projectId][buttonKeyword] = widget
        self.widgetInfo[projectId][buttonKeyword] = selectionDict

        self.widgets[projectId][buttonKeyword].grid(row=row, column=colNum, sticky=Tkinter.EW)

        row += 1
          
      widgetEndRow = row
      
      button = Tkinter.Button(master, text = "Save %s" % projectId, command = lambda projectId = projectId:self.saveProject(projectId))
      button.grid(row=row, column=colNum, sticky=Tkinter.EW)
      
      if self.testMode:
        row += 1
        button = Tkinter.Button(master, text = "Test", command = lambda projectId = projectId:self.doTest(projectId))
        button.grid(row=row, column=colNum, sticky=Tkinter.EW)

      if project:
        self.setNewProjectName(projectId,'directoryName')
    
    #
    # Set Help button
    #
    
    row += 1
    buttons = createDismissHelpButtonList(master, dismiss_text = 'Close', help_url=self.help_url, expands = True)
    buttons.grid(row=row, column=0, columnspan=4)
    
    #
    # Set shifter arrows
    #
    
    label = Label(master, text = "Copy")
    label.grid(row=self.widgetStartRow + 1, column=1, columnspan = 2, sticky=Tkinter.EW)

    self.shiftButtons = {}
    
    if self.project1 and self.project2:
      buttonState = "normal"
    else:
      buttonState = "disabled"

    widgetRowIndex = 0
    
    for widgetIndex in range(len(self.widgetDefs)):
    
      widgetInfo = self.widgetDefs[widgetIndex]

      if widgetInfo[0] == 'Separator':
        continue
      
      widgetRowIndex += 1
      buttonKeyword = widgetInfo[1]      
      self.shiftButtons[buttonKeyword] = []
      
      row = widgetRowList[widgetRowIndex]

      for (shiftCol,shiftLabelText) in ((1,'<-'),(2,'->')):
        
        button = Tkinter.Button(master, text = shiftLabelText, command = lambda buttonKeyword = buttonKeyword, toProjectNum = shiftCol:self.shiftData(buttonKeyword,toProjectNum), state = buttonState)
        button.grid(row=row, column=shiftCol, sticky=Tkinter.EW)
        
        self.shiftButtons[buttonKeyword].append(button)

  
  def saveProject(self,projectId):
  
    project = getattr(self,projectId)
    
    if project:

      try:
          
        project.checkAllValid(complete = True)
        
        try:
          project.saveModified()
        except:
          showError("Project could not be saved","The project could not be saved. Sorry - try again from scratch.")
          print traceback.format_exception_only(sys.exc_type,sys.exc_value)
      except:
        showError("Project not valid","The project is not valid and will not be saved. Sorry - try again from scratch.")
        print traceback.format_exception_only(sys.exc_type,sys.exc_value)
   
  def doTest(self,projectId):
  
    project = getattr(self,projectId)
    
    if project:
      for repository in project.repositories:
        print repository.name, repository.url
        for packageLocator in repository.sortedStored():
          print "     ", packageLocator.targetName
          
  def shiftData(self,buttonKeyword,toProjectNum):
  
    if toProjectNum == 1:
      fromProjectId = 'project2'
    else:    
      fromProjectId = 'project1'
      
    toProjectId = "project%d" % toProjectNum
    
    setupInfo = None
    
    for widgetIndex in range(len(self.widgetDefs)):
    
      widgetInfo = self.widgetDefs[widgetIndex]
      
      if widgetInfo[0] == 'Separator':
        continue

      if buttonKeyword == widgetInfo[1]:
        setupInfo = widgetInfo
        break
      
    ccpnObject = getIoInfo(self.widgets[fromProjectId][buttonKeyword],setupInfo,getattr(self,fromProjectId),self.widgetInfo[fromProjectId][buttonKeyword])

    if ccpnObject:
    
      try:
        toProject = getattr(self,toProjectId)
        
        objectClassName = getClassOrSuperClassName(ccpnObject)
                
        try:
          CopyClass = globals()["Copy%s" % objectClassName]
        except:
          showError("Not implemented","Cannot copy %s objects yet..." % objectClassName)
        
        dialogInfo = (self.guiParent,self.multiDialog,self.messageReporter,self.dataEntry)
        copyClass = CopyClass(ccpnObject,toProject,dialogInfo,{},verbose = self.verbose)  
        
        for resonanceParent in toProject.findAllNmrProjects().union(toProject.findAllNmrConstraintStores()):

          newResonancesFlag = resonanceParent.findFirstApplicationData(keyword = newResonances_kw)
        
          if newResonancesFlag and toProject.molSystems:
            if showYesNo("Run linkResonances","Unlinked resonances are present - do you want to link them now?"):
              if resonanceParent.className == 'NmrConstraintStore':
                nmrConstraintStore = resonanceParent
              else:
                nmrConstraintStore = None
              self.multiDialog.LinkResonancesSetup(self,toProject,nmrConstraintStore)
              resonanceParent.removeApplicationData(newResonancesFlag)
        
        print "Checking new object validity - this might take a while."
        copyClass.toObject.checkAllValid(complete = True)
        print "Checking done."
        
      except:
        raise
        print traceback.format_exception_only(sys.exc_type,sys.exc_value)
        showError("Could not copy","Copying resulted in a '%s'. See shell for details." % sys.exc_info()[0])
   
      self.updateProjectWidgets(toProjectId)

  def updateProjectWidgets(self,projectId):
  
    project = getattr(self,projectId)
    
    for widgetIndex in range(len(self.widgetDefs)):
    
      (widgetType,buttonKeyword,labelText,IOdefault,extraDefs) = self.widgetDefs[widgetIndex]

      if widgetType == 'Separator':
        continue

      widget = self.widgets[projectId][buttonKeyword]
      
      if widgetType[:13] == 'SelectionList' or widgetType[:18] == 'MultiSelectionList':

        ccpnObjects = list(getCcpnObjects(project,IOdefault))
        
        (selectionList,selectionDict) = createSelection(ccpnObjects)

        widget.setup(selectionList + ['None'],None,None)
        selectionDict['None'] = None

        self.widgetInfo[projectId][buttonKeyword] = selectionDict
        
    self.updateShiftButtons()
    
  def updateShiftButtons(self):
  
    buttonState = "disabled"
    
    if self.project1 and self.project2:
      buttonState = "normal"
    
    for buttonKey in self.shiftButtons.keys():
      for shiftButton in self.shiftButtons[buttonKey]:
        shiftButton.config(state = buttonState)

  def selectFile(self,projectId,void2,buttonKeyword,void3,void4,multiSelect = False):
  
    popup = FormatFilePopup(self, component = 'project', format = 'ccpn', title = 'Select CCPN project directory')
    
    if popup.fileSelected:
      try:
        project = loadProject(popup.file)
        
        setattr(self,projectId,project)
        
        self.setNewProjectName(projectId,buttonKeyword)
        self.updateProjectWidgets(projectId)
                
      except:
        raise
        showError('Invalid file','Could not load project from file %s' % popup.file)

  def selectProjectDir(self,projectId):
  
    popup = OpenProjectPopup(self, callback = lambda project, projectId = projectId: self.initProject(project,projectId=projectId))
    
  def initProject(self,project, projectId = None):
    
    if project:
      setattr(self,projectId,project)
      self.setNewProjectName(projectId,'directoryName')
      self.updateProjectWidgets(projectId)
    
  def setNewProjectName(self,projectId,buttonKeyword):

    self.widgets[projectId][buttonKeyword].config(text = getattr(self,projectId).name)
    self.widgets[projectId][buttonKeyword + "Label"].set("Project name:")
  
  def apply(self):
  
    return True

if __name__ == '__main__':
  import os

  projDirs = sys.argv[1:]
  
  addKeywds = {}
  if projDirs:
    for i in range(len(projDirs)):
      projDir = projDirs[i]
      try:
        project = loadProject(projDir)
        addKeywds["project%d" % (i+1)] = project
      except:
        print "No project found in %s - not loading" % projDir

  guiRoot = Tkinter.Tk()
  
  ds = DataShifter(guiRoot, **addKeywds)
  
  guiRoot.mainloop()
  os._exit(0)
"""
======================COPYRIGHT/LICENSE START==========================

Converters1.py: Data compatibility handling

Copyright (C) 2007 Rasmus Fogh (CCPN project)
 
=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license.
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

======================COPYRIGHT/LICENSE END============================

To obtain more information about this code:

- CCPN website (http://www.ccpn.ac.uk)

- contact Rasmus Fogh (ccpn@bioc.cam.ac.uk)

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following reference:

===========================REFERENCE START=============================
Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and 
automated software development. Bioinformatics 21, 1678-1684.
===========================REFERENCE END===============================
 
"""

import gc
import time
# imports
import os
from memops.general.Implementation import ApiError
from memops.universal import Io as uniIo
from memops.metamodel import ImpConstants

from memops.general.Constants import currentModelVersion
currentVersionStr = str(currentModelVersion)
   
def majorUpgradeToCurrent(oldRoot, oldVersionStr, newName=None, newDir=None,
                          oldTags=None, doSave=True):
  """
  Do major upgrade of oldRoot tree to current version.
  NB oldRoot tree must be all loaded in memory
  """
  
  if oldVersionStr == currentVersionStr:
    return oldRoot
  
  from memops.api import Implementation
  
  # make new project and move old directory if necessary
  oldName = oldRoot.name
  if newName is None:
    newName = oldName
  
  newRoot = Implementation.MemopsRoot(name=newName)
  
  if newDir is None:
    newPath = getDataDir(oldRoot, oldVersionStr)
  else:
    newPath = uniIo.normalisePath((os.path.join(newDir, newName)), 
                                   makeAbsolute=True)
  newRoot.packageLocator.findFirstRepository().url = (
   Implementation.Url(path=newPath)
  )
  newRoot.findFirstRepository(name='backup').url = (
   Implementation.Url(path=newPath+'_backup')
  )
  
  upgrader = DataUpgrader(oldVersionStr, newRoot, oldRoot, oldTags=oldTags)
  
  # make upgrade
  return upgrader.majorUpgrade(doSave=doSave)


class DataUpgrader(object):
  """class to handle backwards and forwards compatibility for data
  
  NB currently handles only compatibility to whichever version is current
  as determined by the Python import path
  """
  
  class DummyObject:
    """Superclass for dummy objects made during major upgrades
    """
  
    def getFullKey(self):
      """ needed for mapping of reference objects
      """
      raise ApiError("getFullKey must be implemented in subclass")
    
    def getLocalKey(self):
      raise ApiError("getLocalKey must be implemented in subclass")
  
  
  def __init__(self, oldVersionStr, newRoot=None, oldRoot=None,
                newVersionStr=None, oldTags=None):
    """old and new root, old and new versionStr are pre and post upgrade
    memopsRoot and version string, respectively.
    globalMapping is the top level I/O map to use 
    - normally not passed but set here.
    oldTags are the names of the links from oldRoot to all package topObjects,
    in import order.
    """
    
    self.oldRoot = oldRoot
    self.newRoot = newRoot
    self.oldVersionStr = oldVersionStr
    if newVersionStr is None:
      newVersionStr = currentVersionStr
    elif newVersionStr != currentVersionStr:
      raise ApiError(" DataUpgrader to non-current version %s not implemented"
                     % newVersionStr)
    self.newVersionStr = newVersionStr
    self.oldTags = oldTags
    
    # Version-specific compatibility code module.
    self.compatibilityModule = getCompatibilityModule(oldVersionStr)
    

  def majorUpgrade(self, doSave=False):
    """
    """
    
    try:
    
      xtime = time.time()
      wasEnabled = gc.isenabled()
      if wasEnabled:
        gc.disable()
 
      newRoot = self.newRoot
 
      # old to new data object map
      self.oldToNew = {self.oldRoot:newRoot}
 
      # Modify old object hierarchy in place
      self.massageData()
 
      # load all new data before modifying IO map
      from memops.general.Util import loadAllData
      loadAllData(newRoot)
      
      tt = xtime
      xtime = time.time()
      print('+++ Done load new data - time %s' % (xtime - tt))
    
 
      # check for existing files, move old file as needed, and save
      oldPath = getDataDir(self.oldRoot, self.oldVersionStr)
      newPath = getDataDir(newRoot, self.newVersionStr)
 
      # save destination tests. NB done regardless to avoid accidental overwrites
      if oldPath == newPath:
        # using same name - change old version to '_bak'
        dest = oldPath + '_bak'
        if os.path.exists(dest):
          raise ApiError("Cannot move old data, %s already exists" % dest)
 
        # move old data - must be done here to avoid accidental reading instead
        # files on oldPath location
        if os.path.exists(oldPath):
          os.rename(oldPath, dest)
 
      # check that new locations are free:
      if os.path.exists(newPath):
        raise ApiError("Cannot save new data, %s already exists" % newPath)
      
      # Make new data and do minor upgrade
      self.transferData()
    
    finally:
      if wasEnabled:
        gc.enable()
      tt = xtime
      xtime = time.time()
      print('+++ Done transfer data - time %s ' % (xtime - tt))
    
    # check validity
    # NB complete=True checking will be slow - may be relaxed later
    newRoot.checkAllValid(complete=True)
    tt = xtime
    xtime = time.time()
    print('+++ Done complete checkValid time %s' % (xtime - tt))
 
    if doSave:
      # save new data
      newRoot.saveModified()
   
      tt = xtime
      xtime = time.time()
      print('+++ Done save result - time %s' % (xtime - tt))
    #
    return newRoot


  def massageData(self):
    """ Change old data in place and prepare new data so that they match
    """
    
    self.compatibilityModule.Major.modifyData(self)
      
  
  def mapAllReferenceData(self, topObjList, globalMapping):
    """ map reference data and remove mappable TopObjects from topObjList
    NB it is assumed that topObjList is in revese inheritance order
    """
    newRoot = self.newRoot
    newTopObjByGuid = newRoot.__dict__['topObjects']
    
    mapsByGuid = globalMapping['mapsByGuid']
    
    fullKeyConverters = dd = {}
    if self.compatibilityModule is not None:
      dd.update(self.compatibilityModule.General.fullKeyConverters)

    
    self.adaptIoMap(globalMapping)
    
    
    for ii in range(len(topObjList)-1, -1, -1):
      # loop over topObjList from end, so we can remove processed objs.
      oldObj = topObjList[ii]
      obj = None
      
      rootName = oldObj.root.__class__.__name__
      
      # look for corresponding new object
      # version 2+ - map based on guid
      obj = newTopObjByGuid.get(oldObj.guid)
      
      if obj is not None:
        # new object already there - map to existing data
        self.mapReferenceData(oldObj, obj, fullKeyConverters, globalMapping)
        # remove from oldTopObjects as no minor compatibility needed
        topObjList.remove(oldObj)
        
        
  def mapReferenceData(self, oldTopObj, newTopObj, fullKeyConverters, 
                       globalMapping):
    """ Map reference data for a single TopObject.
    Part 1 - connect old objects with their fullKeys
    NB it is assumed that all imported packages are already OK.
    """
    packageShortName = oldTopObj.packageShortName
    dataTypeMap = globalMapping[packageShortName]['abstractTypes']
    delayObjects = []
    fullKeys = []
    

    # stack of child objects to map - these are old objects
    oldChildStack = [[oldTopObj]]
    
    while oldChildStack:
      ll = oldChildStack[-1]
      if ll:
        oldObj = ll.pop()
        curMap = dataTypeMap[oldObj.__class__.__name__]
        
        if curMap.get('proc') == 'skip':
          oldObj.__dict__.clear()
        
        else:
          childList = []
          oldChildStack.append(childList)
          contDict = curMap['content']
          tags = curMap.get('children')
          if tags:
            for tag in tags:
              val = getattr(oldObj, tag)
              if contDict[tag]['hicard'] == 1:
                if val is not None:
                  childList.append(val)
              else:
                childList.extend(val)
 
          if curMap.get('_transf'):
            # Object or one of its descendants needed in oldToNew map
            # keep for processing
            delayObjects.append(oldObj)
            fullKeys.append(oldObj.getFullKey())
          else:
            # object or children will not be needed. Liberate memory
            oldObj.__dict__.clear()
 
      else:
        # no children left - go up a step
        oldChildStack.pop()
    
    # clear delayObjects (NB object itself still needed for oldToNew dict)
    # henceforth no per-object information is available
    for obj in delayObjects:
      obj.__dict__.clear()
    
    self.connectReferenceData(delayObjects, fullKeys, packageShortName, 
                              fullKeyConverters, str(oldTopObj), globalMapping)
    
    
  def connectReferenceData(self, delayObjects, fullKeys, packageShortName,
                           fullKeyConverters, topObjString, globalMapping):
    """ Map reference data for a single TopObject.
    Part 2 - set up new-to-old and old-to-new maps
    NB it is assumed that all imported packages are already OK,
    and that delayObjects are child-before-parent
    """
    
    oldToNew = self.oldToNew
    
    # get OldMemopsObject class
    if delayObjects:
      oldObj = delayObjects[-1]
      ll = [x for x in oldObj.__class__.__mro__ 
            if x.__name__ in ('MemopsObject','MemopsBaseClass') ]
      if len(ll) == 1:
        # MemopsObject for old object tree 
        # - done this way to avoid import trouble
        OldMemopsObject = ll[0]
      else:
        raise ApiError("superType 'MemopsObject' not found from object %s" 
                       % oldObj)
    else:
      return
 
    # dereference fullKeys that include links, as long as we can
    found = True
    while found:
      # keep looping over list while we find new things to map
      found = False
      
      for ii in range(len(delayObjects) - 1, -1, -1):
        # NB loop over list from end to be able to shorten during loop
        
        oldObj = delayObjects[ii]
        
        fullKey = fullKeys[ii]
        for jj, xx in enumerate(fullKey):
          if isinstance(xx, OldMemopsObject):
            # is there an old CCPN object in the key?
            # NB new objects will not be found as their MemopsObject class
            # is a different object
            xx2 = oldToNew.get(xx)
            if xx2 is None:
              # object not mapped yet - try again later
              break
            else:
              # replace with image
              fullKey[jj] = xx2
              
        else:
          newObj = self.newObjFromOldFullKey(oldObj, fullKey, packageShortName,
                                             fullKeyConverters, globalMapping)
          if newObj is not None:
            oldToNew[oldObj] = newObj
            del delayObjects[ii]
            del fullKeys[ii]
            found = True
        
        
    # print diagnostics
    nObjs = len(delayObjects)
    if nObjs == 1:
      print ("""\nFinished mapping reference objects
  %s:
  %s had unmappable keys""" % (topObjString, delayObjects[0]))
    elif nObjs:
      ll = list(set(x.qualifiedName for x in delayObjects))
      print ("""\nFinished mapping reference objects
  %s:
  %s objects of classes %s had unmappable keys""" % (topObjString, nObjs, ll))
  
  
  def transferData(self):
    """ transfer data from old to new object tree, 
    and map between existing and new reference data.
    """
    #from memops.xml import Implementation as xmlImplementation
    # temporary fix pending final move over
    from memops.format.xml import XmlIO as xmlImplementation
 
    # set up
    oldToNew = self.oldToNew
    oldRoot = self.oldRoot
    # not needed after this. Clean up.
    del self.oldRoot
    newRoot = self.newRoot
    
    # new I/O map modified for reding old data
    from memops.xml.Implementation import getGlobalMap
    globalMapping = getGlobalMap(oldVersionStr=self.oldVersionStr)
    
    # get topObjects to start with:
    transferTopObjects = []
    curMap = globalMapping['mapsByGuid'][self.newRoot.metaclass.guid]
    for tag in self.oldTags:
      val = getattr(oldRoot, tag)
      if val:
        if isinstance(val, (set, frozenset, list, tuple)):
          transferTopObjects.extend(x for x in val 
                               if x.packageName != 'memops.Implementation')
        elif val.packageShortName != 'IMPL':
          transferTopObjects.append(val)
    
    # transfer reference data (those that exist already)
    # and remove mappable TopObjects from transferTopObjects
    transferTopObjects.reverse()
    self.mapAllReferenceData(transferTopObjects, globalMapping)
    transferTopObjects.reverse()
    
    # full garbage collect before final object creation
    print('+++ Done full garbage collection - removed %s objects' 
          % gc.collect())
    
    # create new data
    from memops.general.Util import transferData
    for topObj in transferTopObjects:
      transferData(newRoot, topObj, oldToNew=oldToNew, 
                  oldVersionStr=self.oldVersionStr)
      # NBNB TBD check this OK.
      topObj.__dict__['isModified'] = True

  def newObjFromOldFullKey(self, oldObj, fullKey, packageShortName,
                           fullKeyConverters, globalMapping):
    """get object using fullKey
    NB by the time we get here all objects in the fullKey are mapped
    to new objects already.
    """
    oldMeta = oldObj.metaclass
    guid = oldMeta.guid
    # convert fullKey, if needed
    func = fullKeyConverters.get(guid)
    if func is not None:
      fullKey = func(fullKey)
 
    # get new object
    newMap = globalMapping['mapsByGuid'].get(guid)

    #
    if newMap.get('proc') == 'skip':
      return None
    else:
      result = newMap['class'].getByKey(self.newRoot, fullKey)
      
    return result
  
  
  def adaptIoMap(self, globalMapping):
    """  tweak I/O map etc. for use in mapAllReferenceData and
    functions below it. Done only once for each package
    """
    
    transferClasses = set()
    mapsByGuid = globalMapping['mapsByGuid']
    
    # loop over all MetaRoles
    leafPackages = []
    packages = [self.newRoot.metaclass.topPackage()]
    for pp in packages:
      childPackages = pp.containedPackages
      if childPackages:
        packages.extend(childPackages)
      else:
        if pp.name == 'Implementation':
          implPackage = pp
        leafPackages.append(pp)
    
    for pp in leafPackages:
      importedPackages = pp.importedPackages
      for clazz in pp.classes:
        keyNames = clazz.keyNames
        for role in clazz.roles:
          if not role.isDerived:
            
            # add targets of out-of-package links and key links
            vt = role.valueType
            pp2 =  vt.container
            if (pp2 in importedPackages and pp2 is not implPackage or
                role.name in keyNames):
              
              if (pp2 in importedPackages and pp2 is not implPackage):
                # target of interpackage link needed in all contexts
                for xx in vt.getAllSubtypes():
                  # NB subtypes includes vt itself
                  # necessary for mapReferenceData but not for self.oldToNew
                  transferClasses.add(xx)
                  mapsByGuid[xx.guid]['_transf'] = 1
                
                # get parent to avoid looping over all this again below
                pr = vt.parentRole
                if pr is None:
                  continue
                else:
                  vt = pr.valueType
              
              # add valueType, its parents and subtypes
              while vt.container is pp2:
                for xx in vt.getAllSubtypes():
                  if xx not in transferClasses:
                    # necessary for mapReferenceData but not for self.oldToNew
                    transferClasses.add(xx)
                    mapsByGuid[xx.guid]['_transf'] = 2
                pr = vt.parentRole
                if pr is not  None:
                  vt = pr.valueType


def modifyIoMap(oldVersionStr, globalMapping):
  """ Adapt globalMapping to read oldVersion data
  """
  
  compatibilityModule = getCompatibilityModule(oldVersionStr)
  
  # Model information by guid
  anApiClass = globalMapping['IMPL']['abstractTypes']['MemopsObject']['class']
  newElementsByGuid = getElementsByGuid(anApiClass._metaclass.topPackage())
  
  # new elements retrofitted by hand and therefore treated as old 
  # used for e.g. changes of parent tree.
  elemsTreatedAsOld = set(compatibilityModule.General.elemsTreatedAsOld)
  
  loadMaps = globalMapping['loadMaps']
  mapsByGuid = globalMapping['mapsByGuid']
  diffMap = compatibilityModule.MapInfo

  # correct maps for newly introduced elements
  for tag in ('newElements', 'newMandatories', 'neutraliseElements'):
    ll = getattr(diffMap, tag)
    for (prefix, typeName, elemName, newGuid) in ll:
      if prefix in globalMapping and newGuid not in elemsTreatedAsOld:
        removeElementName(newGuid, globalMapping, newElementsByGuid)

  # correct maps for skipped elements
  for (prefix, typeName, elemName, newGuid, elemType) in diffMap.skipElements:
    if prefix in globalMapping and newGuid not in elemsTreatedAsOld:
      fixElementMap(newGuid, 'skip', prefix, globalMapping, newElementsByGuid,
                    typeName, elemName)
      if elemType == 'MetaClass':
        fixExoLinkMap(newGuid, 'skip', prefix, globalMapping, typeName)

  # correct maps for changed valueTypes
  for tt in diffMap.typeChanges:
    (prefix, typeName, elemName, action, newGuid, elemMap, valueTypeGuid) = tt
    if prefix in globalMapping and newGuid not in elemsTreatedAsOld:
      fixElementMap(newGuid, action, prefix, globalMapping, newElementsByGuid,
                         typeName, elemName, valueTypeGuid=valueTypeGuid,
                         elemMap=elemMap, overrideGuidMap=True)

  # correct maps for renamed elements
  # NB must be done after valueTypes corrections, 
  # as elements may be both renamed and 'delay'ed or 'skip'ped
  # by changed valueType 
  for (prefix, typeName, elemName, newName, newGuid) in diffMap.renames:
    if prefix in globalMapping and newGuid not in elemsTreatedAsOld:
      fixElementMap(newGuid, 'rename', prefix, globalMapping, newElementsByGuid,
                    typeName, elemName)
      if elemName is None:
        fixExoLinkMap(newGuid, 'rename', prefix, globalMapping, typeName)
      

  # correct maps for delayed elements
  # NB, must be done after renaming, as we use type names to find 
  # appropriate maps
  for tt in diffMap.delayElements:
    (prefix, typeName, elemName, newGuid, elemMap, valueTypeGuid) = tt
    if prefix in globalMapping and newGuid not in elemsTreatedAsOld:
      fixElementMap(newGuid, 'delay', prefix, globalMapping, newElementsByGuid,
                    typeName, elemName, valueTypeGuid=valueTypeGuid,
                    elemMap=elemMap)
  
  # make extra changes
  compatibilityModule.General.extraMapChanges(globalMapping)
  
def removeElementName(guid, globalMapping, newElementsByGuid):
  """ Remove element from name lists used for selecting input
  and remove element from loadMaps (in case another element with the
  same signature needs the slot).
  """
  
  # set up
  useObj = newElementsByGuid.get(guid)
  if useObj is None:
    # guid does not exist. Probably was not included in model generation. 
    # Anyway, nothing to do.
    return
  newName = useObj.name
  mapsByGuid = globalMapping['mapsByGuid']
  
  # remove from lists in container and its subtypes
  if useObj.__class__.__name__ in ('MetaAttribute', 'MetaRole'):
    
    if (useObj.isAbstract 
        or (useObj.isDerived and useObj.changeability == ImpConstants.frozen)):
      # These elements have no maps, so nothing to do
      return
      
    subtypes = useObj.container.getAllSubtypes()
    for metaclass in subtypes:
      containerMap =  mapsByGuid.get(metaclass.guid)
      if containerMap is not None:
        # NB there might not be a map if this is (e.g.) an abstract superclass 
        # that does not inherit from DataObject 
        for tag in ('simpleAttrs', 'headerAttrs', 'optLinks', 'cplxAttrs',
                    'children'):
          ll = containerMap.get(tag)
          if ll and newName in ll:
            ll.remove(newName)
  
    # remove from loadMaps
    oldMap = mapsByGuid.get(guid)
    if oldMap is None:
      # check if there is need for concern
      if (useObj.__class__.__name__ == 'MetaRole' and 
          useObj.container.container in useObj.valueType.container.importedPackages):
        # interpackage link in backwards direction
        pass
      elif (useObj.__class__.__name__ == 'MetaRole' and 
            useObj.hierarchy == 'parent'):
        pass
      elif useObj.container.guid not in mapsByGuid:
        # The class is not mapped. 
        # Could be (is?) abstract multi-inheritance class
        pass
      else:
        print 'WARNING, no previous map found for %s (%s)' % (useObj, guid)
    else:
      del globalMapping['loadMaps'][oldMap['tag']]
  
  
def fixElementMap(newGuid, action, prefix, globalMapping, newElementsByGuid, 
                  typeName=None, elemName=None, overrideGuidMap=False, 
                  valueTypeGuid=None, elemMap=None):
  """ insert/modify I/O map in mapping 
  - newGuid is the guid to use in the new map 
  (new guid if element exists, old guid otherwise) 
  - prefix is the package shortname
  - typeName is the name of the class/DataObjType (if any) 
  - elemName is the name of the ClassElement (if any) 
  - overrideGuidMap controls if an existing map should be overridden
  - valueTypeGuid is the guid for a new valueType (for attr)
  - elemMap is a passed-in ({typ, proc, name, eType}) I/O map for the element
  
  """
  
  from memops.format.xml import XmlGen
  
  if action == 'ignore':
    return
  
  elif action not in ('rename', 'skip', 'delay'):
    raise ApiError("action: %s not supported" % action)
  
  loadMaps = globalMapping['loadMaps']
  mapsByGuid = globalMapping['mapsByGuid']
  mapping = globalMapping[prefix]
  
  if typeName is None:
    return
  
  xmlTag = XmlGen.xmlTag(prefix, typeName, elemName)
  if elemName is None:
    typ = 'dataType'
  else:
    typ = 'elem'  
  
  if action == 'skip':
    # we are skipping the element
    if newGuid in mapsByGuid:
      raise ApiError("%s: guid already has map" % xmlTag)
    else:
      mapsByGuid[newGuid] = newMap = {'guid':newGuid, 'proc':'skip'}
    
    if typ == 'dataType':
      # put map in abstractTypes
      abstractTypes = mapping['abstractTypes']
      if typeName in abstractTypes:
        raise ApiError("%s: name %s already in use" % (xmlTag, typeName))
      else:
        abstractTypes[typeName] = newMap
  
  elif action == 'rename':
    # this is a renaming
    
    # get previous I/O map (must exist)
    newMap = mapsByGuid.get(newGuid)
    if newMap is None:
      # guid does not exist. Probably was not included in model generation. 
      # Anyway, nothing to do.
      return
    
    # clean up loadMaps
    del loadMaps[newMap['tag']]
    
    if typ == 'elem':
      # ClassElement
      newObj = newElementsByGuid[newGuid]
      newName = newObj.name
      
      if newName != elemName:
        # they may be equal if the container class only is being renamed
        
        # fix attribute lists and contents for container and its subtypes
        subtypes = newObj.container.getAllSubtypes()
        for metaclass in subtypes:
          fixMap = mapsByGuid[metaclass.guid]
 
          # fix content map
          contDict = fixMap['content']
          if elemName in contDict:
            raise ApiError("%s: name %s already in use" % (xmlTag, elemName))
          contDict[elemName] = contDict[newName]
          del contDict[newName]
 
          # fix attribute lists
          for tag in ('simpleAttrs', 'headerAttrs', 'optLinks', 'cplxAttrs',
                      'children'):
            ll = fixMap.get(tag)
            if ll and newName in ll:
              ll[ll.index(newName)] = elemName
    else:
      # classifier - stil needs renaming
      # as we use the type name to find the right map
      # put map in abstractTypes
      abstractTypes = mapping['abstractTypes']
      if typeName in abstractTypes:
        raise ApiError("%s: name %s already in use" % (xmlTag, typeName))
      else:
        abstractTypes[typeName] = newMap
    
  else:
    # action == 'delay'
    # we are delaying storage for use in compatibility
    if typ != 'elem':
      raise ApiError("%s: only ClassElements may be 'delay'ed" % xmlTag)
    
    if valueTypeGuid is not None:
      elemMap['data'] = mapsByGuid[valueTypeGuid]
    
    if overrideGuidMap:
      newMap = mapsByGuid[newGuid]
    
    elif newGuid in mapsByGuid:
      raise ApiError("%s: guid %s already in use" % (xmlTag, newGuid))
    
    elif xmlTag in loadMaps:
      raise ApiError("%s: xml tag already in use" % (xmlTag,))
    
    else:
      mapsByGuid[newGuid] = newMap = {}
      containerGuid = mapping['abstractTypes'][typeName]['guid']
      newContainer = newElementsByGuid[containerGuid]
      subtypes = newContainer.getAllSubtypes()
      for metaclass in subtypes:
        # put new map in content dictionary
        fixMap = mapsByGuid[metaclass.guid]
        contDict = fixMap['content']
        if elemName in contDict:
          # we are overriding the previous definition.The tag is already there
          pass
        else:
          # put tag on cplxAttrs (it does not matter which list it is on)
          contDict[elemName] = newMap
          ll = fixMap.setdefault('cplxAttrs', [])
          ll.append(elemName)
          
    # update map contents
    newMap.update(elemMap)
    newMap['proc'] = 'delay'
  
  # fix loadMaps
  loadMaps[xmlTag] = newMap
  
  
def fixExoLinkMap(newGuid, action, prefix, globalMapping, typeName):
  """ Fix exoLink maps
  if action == 'rename' enter newGuid map under (prefix, typeName) name
  if action == 'skip', enter skip record for guid under (prefix, typeName) name
  """
  
  from memops.format.xml import XmlGen
  
  exoTag = XmlGen.xmlTag(prefix, typeName, var='exo')
  
  # set up
  exolinks = globalMapping[prefix]['exolinks']
  loadMaps = globalMapping['loadMaps']
  
  # check if slots are occupied
  if typeName in exolinks:
    raise ApiError("%s exolink: name %s already in use" % (exoTag, typeName))
  if exoTag in loadMaps:
    raise ApiError("%s exolink: name already in use" % exoTag)
  
  if action == 'skip':
    # we are not replacing anything, but skippping
    newMap = {
     'type':'exo',
     'guid':newGuid,
     'eType':'cplx', 
     'proc':'skip'
    }
  
  else:
    # replacing an existing map
    newClassMap = globalMapping['mapsByGuid'].get(newGuid)
    newExoTag =  XmlGen.xmlTag(prefix, newClassMap['class'].__name__, var='exo')
    #newMap = globalMapping['mapsByGuid'].get(newGuid)
    newMap = loadMaps[newExoTag]
    del exolinks[newMap['class'].__name__]
    #del loadMaps[newMap['tag']]
  
  # set new map
  exolinks[typeName] = newMap
  loadMaps[exoTag] = newMap



def getDataDir(rootObj, versionStr):
  """Get directory that contains user data
  """
  
  locator = rootObj.findFirstPackageLocator(targetName='any')
  result = locator.findFirstRepository().url.path
  #
  return result
  
  
def getElementsByGuid(rootPackage):
  """ get guid:obj dictionary for MetaModelElements
  in MetaModel corresponding to rootPackage
  
  Used in mapping of reference data.
  NBNB may need version-specific expansion later to handle 
  difficult reference data modifications.
  """
  from memops.metamodel.MetaModel import makeObjDict
  #
  return makeObjDict(rootPackage, ignoreImplicit=True, crucialOnly=True)

def getCompatibilityModule(fromVersion):
  """ import relevant compatibility module
  """
  compModules = ['Major', 'Minor', 'General', 'MapInfo']
  
  from memops.general.Version import cmpVersionStrings
  
  newVersion = currentVersionStr
  
  compar = cmpVersionStrings(fromVersion, newVersion)
  if compar < 0:
    # fromVersion is older 
    srcDir = 'upgrade'
  
  elif compar > 0:
    # fromVersion is newer 
    srcDir = 'downgrade'
  
  else:
    # versions identical
    return
    
  ss = fromVersion.replace('.', '_')
  moduleDir = "memops.format.compatibility.%s.v_%s" % (srcDir, ss)
  try:
    result = __import__(moduleDir, {}, {}, compModules)
  except ImportError:
    if compar < 0:
      raise ApiError("""Could not import %s
Backwards compatibility from %s to %s missing.""" 
       % (moduleDir, fromVersion, newVersion))
    else:
      raise ApiError("""Could not import %s
Forwards compatibility from %s to %s missing.
Trying to load new-version data with old-version code?
""" 
       % (moduleDir, fromVersion, newVersion))
      
  
  for ss in compModules:
    if not hasattr(result, ss):
      raise ApiError("Compatibility error: No module %s found in %s" 
                     % (ss, moduleDir))
  #
  return result

def minorPostProcess(fromVersionStr, topObj, delayDataDict, toNewObjDict,
                     mapping=None, topObjElem=None):
  """ postProcess - update newRoot object tree
  May be used either to postprocess a file load (minor upgrade)
  or as part of an in-memory data transfer (major upgrade)
  
  topObj is a TopObject in the new tree
  delayDataDict is a {newobj:{tag:value}} dictionary of delayed set elements
  and children
  toNewObjDict is _ID:newObj for minor and oldObj:newObj for major upgrades
  """
  
  mod = getCompatibilityModule(fromVersionStr)
  if mod is not None:
    if topObjElem:
      # Done like this to avoid rewriting correctData functions that do
      # not need topObjElem
      mod.Minor.correctData(topObj, delayDataDict, toNewObjDict, mapping, 
                            topObjElem)
    else:
      mod.Minor.correctData(topObj, delayDataDict, toNewObjDict, mapping)



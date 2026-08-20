"""
======================COPYRIGHT/LICENSE START==========================

Util.py: code for CCPN data model and code generation framework

Copyright (C) 2005  (CCPN Project)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
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

- email: ccpn@bioc.cam.ac.uk

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

"""

######################################################################
# hack for Python 2.1 compatibility  NBNB                                                      #
######################################################################
try:
  junk = True
  junk = False
except:
  True = not 0
  False = not 1
  
# miscellaneous useful functions

import os
import time

from memops.metamodel import ImpConstants
from memops.general import Constants as genConstants
# NB cannot import Impelmentation as it is not compatible with Python 2.1 
# (for ObjectDomain)
#from memops.general.Implementation import ApiError

from memops.universal import Io as uniIo

def returnMemopsWord(value):
 
  """
  Converts a string to a memops dataType 'Word'
  max length is 32, no whitespace allowed (is taken out!)
  """
 
  if value:
    
    # NB changed - did not fulfillits promise
    #wordString = value.replace(' ','')
    wordString = ''.join(value.split())
    wordString = wordString[:32]
 
    return wordString
 
  else:
 
    return value

def returnMemopsLine(value):
 
  """
  Converts a string to a memops dataType 'Line'
  max length is 80
  """
 
  if value:
 
    wordString = value.replace(os.linesep,' ')
    wordString = wordString[:80]
 
    return wordString
 
  else:
 
    return value

def returnMemopsText(value):
 
  """
  Converts a string to a memops dataType 'Text'
  max length is 254
  """
 
  if value:
 
    wordString = value[:254]
 
    return wordString
 
  else:
 
    return value

def returnMemopsTexts(value):
 
  """
  Converts a string or list/tuple of strings to a memops
  dataType 'Text'
  """
 
  if type(value) == type(()):
    value = list(value)
  elif value and type(value) != type([]):
    value = [value]
  
  if value:
    for i in range(0,len(value)):
      value[i] = returnMemopsText(value[i])
 
  return value

def isProjectModified(project):
  """
  Checks whether any of project has been modified.
  """

  topObjects = (project,) + tuple(project.topObjects)

  for topObject in topObjects:
    if topObject.isModified:
      return True

  return False

def getRepositoryRevision(metaObj):
  """
  Get repository revision string for MetaObj
  """
  
  unknown = '?'
  
  idString = metaObj.taggedValues.get('repositoryId')
  if idString:
    return getRepositoryInfo(idString).get('revision', unknown)
  else:
    return unknown

def getRepositoryInfo(repositoryId):
  """
  Takes Id string as created (and maintained) by CVS and returns a dictionary
  with keys "date", "file", "revision", "state" and values the corresponding
  ones determined from the Id string.
  """
  fields = repositoryId.split()
  if (len(fields) >= 7):
    date = fields[3] + " " + fields[4]
    file = fields[1][:-2]  # strips off irritating ",v"
    revision = fields[2]
    source = fields[5]
    state = fields[6]
    info = { "date": date, "file":file, "revision":revision, "source":source,
     "state":state 
    }
  else:
    info = {}

  return info

def getRepositoryTag(repositoryTag):
  """
  Takes tag (Name) string as created (and maintained) by CVS and returns
  the same with the garbage stripped out and '_' replaced with '.'.
  Assumes input string of form '$-N-a-m-e: something $' (without the '-')
  """

  from memops.general.Constants import repositoryTagPrefix

  repository_prefix = '$Name:'
  repository_suffix = '$'
  tag = repositoryTag[len(repository_prefix):-len(repository_suffix)].strip()

  tag = tag[len(repositoryTagPrefix):]

  tag = tag.replace('_', '.')

  return tag


def configParameter(keyWord):
  """Get configuration parameter corresponding to 'keyWord'.
  
  NBNB TBD this should go to a configuration file to get its data
  The format of this file is still up for grabs.
  """
  
  rootDir = uniIo.normalisePath(uniIo.os.getcwd())
  
  
  if keyWord == 'repositories':
    # list of (name,urlPath) pairs
    return []

  
  elif keyWord == 'packageLocators':
    # list of (packageName, repositoryName) Triplets
    # If the packageName appears more than once, the named repositories
    # are added in the order given
    return [
('ccp.molecule.ChemComp', 'userData'),
('ccp.molecule.ChemComp', 'refData'),
('ccp.molecule.ChemCompCharge','userData'),
('ccp.molecule.ChemCompCharge','refData'),
('ccp.molecule.ChemCompCoord', 'userData'),
('ccp.molecule.ChemCompCoord', 'refData'),
('ccp.molecule.ChemCompLabel', 'userData'),
('ccp.molecule.ChemCompLabel', 'refData'),
('ccp.molecule.ChemElement', 'refData'),
('ccp.molecule.StereoChemistry', 'refData'),
('ccp.nmr.NmrExpPrototype', 'userData'),
('ccp.nmr.NmrExpPrototype', 'refData'),
('ccp.nmr.NmrReference', 'refData'),
('ccpnmr.AnalysisProfile', 'generalData'),
('ccpnmr.AnalysisProfile', 'refData'),
]
  
  else:
    return None


def copySubTree(sourceObj, newParent, maySkipCrosslinks=False,
 topObjectParameters = None, objectMap=None
):
  """ Descrn: Copies subtree rooted in sourceObj to a subtree rooted in 
              a new targetObj that is a child of newParent
              Recursively copies all children.
     Inputs: - top object to be copied, 
             - parent of new copy, 
             - Boolean: if a link cannot be validly copied, can the function
               omit it (or must it throw an error),
             - dictionary of parameters to be passed to copy of top object
             - dictionary of oldObject:newObject mappings. 
             This is for cases where links from the old tree to certain objects
             should be replaced with links from the new tree not to the same
             objects but to a a different set of (pre-existing) objects. E.g.
             when copying peak lists from one experiment to another you have to
             give the mapping from old DataDimRefs to new DataDimRefs.

     Output: new copy of top object.
  
  (Parts of) crosslinks to objects within the subtree are copied to
  links to the new object copies;
  (Parts of) crosslinks to objects not within the subtree are copied to
  links to the old objects, provided this can be done without cutting 
  pre-existing links.
  If the above will not work *and* maySkipCrosslinks is True, the routine
  tries to set the crosslink using only the objects within the subtree.
  If none of the above works, an error is thrown.
  
  The key,val pairs in the topObjectParameters dictionary are passed to 
  the top object constructor, and the pre-existing values in the sourceObj 
  are ignored. This can be used to set new values for the keys of sourceObj.
  
  If the top object has 'serial' as the key and no valid serial is passed in
  topObjectParameters, the routine will set the serial to the next available
  value.
  
  Note that the function first builds all objects, then connects crosslinks,
  then connects parent-to-child links. Finally all notifiers are called but in
  random order. If there is an error the routine tries to delete all created
  objects before re-raising the original error. A a failed function call may
  consume serial numbers if the key of the sourceObj is 'serial'. Also, there is
  a relatively high bug risk, as is always the case with functions that have to
  clean up after an error. 
  """

  from memops.metamodel.MetaModel import MemopsError
  if sourceObj.root is sourceObj:
    raise MemopsError("copySubTree cannot be used to copy entire projects")
      
  result = transferData(newParent, sourceObj, oldToNew=objectMap, 
                        targetObjParams=topObjectParameters,
                        ignoreMissing=maySkipCrosslinks, useOptLinks=True)
  #
  return result
  
 
def newGuid(prefix = ''):
 
  if prefix:
    prefix = '%s_' % prefix

  # slightly modelled on memops.api.Implementation.MemopsRoot.newGuid()

  timeStamp = time.strftime("%Y_%m_%d_%H_%M_%S")
  user = os.environ.get('USER', 'unknown')

  guid = '%s%s_%s' % (prefix, user, timeStamp)

  return guid


def transferData(newParent, sourceObj, oldToNew=None, 
                 oldVersionStr=None, targetObjParams=None,
                 ignoreMissing=True, useOptLinks=False):
  """ Copy sourceObj and recursively all its children,
  to a new tree where the new targetObj is a child of newRoot
  - If oldVersionStr is set, do as  backwards compatibility, 
  including minor post-processing, otherwise do as subtree copying
  - targetObjParams: parameters to be passed to the copy of sourceObj.
    Only meaningful for subtree copy, and ignored for
    backwards compatibility.
  - oldToNew is an old-to-new-object dictionary, serves for either
  - useOptLinks controls if optional links (basically the -to-one
    direction of one-to-many links) should be followed. For compatibility
    this is awaste of time (but harmless), but for copySubTree it is
    necessary
  """
  
  from memops.metamodel import ImpConstants
  
  from memops.xml import Implementation as xmlImplementation
  # temporary fix pending final move over:
  #from memops.format.xml import XmlIO as xmlImplementation
  
  serialTag = ImpConstants.serial_attribute
  # NB serialDictTag hardwired to avoid using the varNames dictionary
  serialDictTag = '_' + ImpConstants.serialdict_attribute
    
  globalMapping = xmlImplementation.getGlobalMap(oldVersionStr)
  
  mapsByGuid = globalMapping['mapsByGuid']
  
  if targetObjParams is None:
    targetObjParams = {}
  
  if oldToNew is None:
    oldToNew = {}
  
  # decide which links to follow
  if useOptLinks:
    followTags = ('headerAttrs', 'simpleAttrs', 'cplxAttrs', 'optLinks')
  else:
    followTags = ('headerAttrs', 'simpleAttrs', 'cplxAttrs')
  
  localOldToNew = {}
  emptyList = []
  oneElemList = [None]
  emptySet = set()
  crossLinkData = []
  appCrossLinkData = crossLinkData.append
  delayDataDict = {}
  # stack of child objects to map - these are old objects
  oldChildStack = [[sourceObj]]
  # stack of parent objects to attach to - these are new objects
  newParentStack = [newParent]
  # objects to notify on - for correct ordering of notifiers when copying subtree
  # should not put newParent in notifyObjects because already exists
  #notifyObjects = [newParent]
  notifyObjects = []
  
  targetObj = None
  try:
 
    nextDd = {}
    while oldChildStack:
      ll = oldChildStack[-1]
      if ll:
        oldObj = ll.pop()
 
        # current object map
        ss = oldObj.packageShortName
        curMap = globalMapping[ss]['abstractTypes'][oldObj.__class__.__name__]

        if curMap.get('proc') == 'skip':
          # skip this one
          continue
 
        # create or get new object
        parent = newParentStack[-1]
 
        if parent is newParent:
          # this is the target object - special case
             
          # fix serial key for sourceObj if copying subtree
          if oldVersionStr is None:
            # we are copying a subtree
            if serialTag in sourceObj.metaclass.keyNames:
              # serial key for top
              if serialTag not in targetObjParams:
                # not being passed in explicitly - we must fix it (special case)
                # NB _serialDict is hardwired to avoid using the varNames di
                serialDict = newParent.__dict__.setdefault(serialDictTag,{})
                oldSerial = serialDict.get(curMap['fromParent'],0)
                targetObjParams[serialTag] = oldSerial + 1
          
          # new object not already there - make it and transfer from old
          if parent.root is newParent:
            # targetObj is a TopObject
            obj = targetObj = curMap['class'](parent, isReading=True,
                                              **targetObjParams)
          else:
            # targetObj is not a TopObject (and we are doing subtree copying)
            parent.topObject.__dict__['isReading'] = True
            obj = targetObj = curMap['class'](parent, **targetObjParams)
          
          notifyObjects.append(obj)
 
          objId = obj
          delayDataDict[objId] = nextDd
          if curMap.get('_transf') == 1:
            oldToNew[oldObj] = obj
          else:
            localOldToNew[oldObj] = obj
          for tag in curMap.get('children', emptyList):
            nextDd[tag] = []
 
          content = curMap['content']
          for tag in curMap.get('cplxAttrs', emptyList):
            if content[tag]['type'] == 'dobj':
              nextDd[tag] = []
 
        else:
          # normal object
 
          if curMap['type'] == 'cplx':
            # complex data type
            obj = curMap['class'](override=True)
 
            objId = id(obj)
            delayDataDict[objId] = nextDd
            localOldToNew[id(oldObj)] = obj
 
          else:
            # type class
            obj = curMap['class'](parent)
 
            objId = obj
            delayDataDict[objId] = nextDd
            if curMap.get('_transf') == 1:
              oldToNew[oldObj] = obj
            else:
              localOldToNew[oldObj] = obj
            delayDataDict[parent][curMap['fromParent']].append(obj)
 
            for tag in curMap.get('children', emptyList):
              nextDd[tag] = []
            
            notifyObjects.append(obj)
 
          # add list for complex data type attrs
          content = curMap['content']
          for tag in curMap.get('cplxAttrs', emptyList):
            if content[tag]['type'] == 'dobj':
              nextDd[tag] = []
 
        # put objects on stack
        childList = []
        oldChildStack.append(childList)
        newParentStack.append(obj)
        contDict = curMap['content']
 
        # transfer object contents
        for ss in followTags:
          tags = curMap.get(ss, emptyList)
          for tag in tags:
            if obj is targetObj and tag in targetObjParams:
              # special case: parameters passed in directly to targetObj
              # needed for tree copying only
              continue
 
            tmpMap = contDict[tag]
            name = tmpMap['name']
            val = getattr(oldObj, tag)
 
            if val is None:
              # no value - skip empties
              continue
 
            elif isinstance(val, (tuple, frozenset)):
              # convert to list for future processing
              if val:
                valIsList = True
              else:
                # no values - skip empties
                continue
 
            else:
              valIsList = False
 
            typ = tmpMap['type']
            
            if typ == 'attr':
              # simple type attribute
 
              proc = tmpMap.get('proc')
 
              if proc == 'delay':
                # pass to compatibility processing - making sure it is a list
                if valIsList:
                  delayDataDict[objId][name] = val
                else:
                  delayDataDict[objId][name] = [val]
 
              else:
                # fix list/nonlist and set
                if tmpMap['hicard'] == 1:
                  if valIsList:
                    for vv in val:
                      break
                    val = vv
                  if proc == 'direct':
                    # direct setting if simple non-constrained attribute
                    obj.__dict__[name] = val
                  else:
                    setattr(obj, name, val)
                
                else:
                  if not valIsList:
                    # optimisation - avoid creating temporary lists
                    oneElemList[0] = val
                    val = oneElemList
                  setattr(obj, name, val)
 
            elif typ == 'child':
              # normal child
              if valIsList:
                childList.extend(val)
              else:
                childList.append(val)
 
            elif typ == 'dobj':
              # normal DataTypeObject
 
              if not valIsList:
                val = [val]
 
              # put on stack for further processing
              childList.extend(val)
 
              if tmpMap.get('proc') == 'delay':
                # delayed - put in delayDataDict
                delayDataDict[objId][name] = val
 
              else:
                # convert to ID and put in crossLinkData to resolve links later
                appCrossLinkData(obj)
                appCrossLinkData([id(xx) for xx in val])
                appCrossLinkData(tmpMap)
 
            else:
              # typ in ('link', 'exolink', 'exotop')
 
              if not valIsList:
                val = [val]
 
              if tmpMap.get('proc') == 'delay':
                delayDataDict[objId][name] = val
 
              else:
                # put in crossLinkData to resolve links later
                appCrossLinkData(obj)
                appCrossLinkData(val)
                appCrossLinkData(tmpMap)
 
        if oldVersionStr is None:
          # this is tree copying
          nextDd = {}
          
        else:
          # this is backwards compatibility
          # clear old object to keep memory use down
          nextDd = oldObj.__dict__
          nextDd.clear()
 
 
      else:
        # no children left - go up a step
        oldChildStack.pop()
        lastParent = newParentStack.pop()
        if lastParent.metaclass.__class__.__name__ == 'MetaDataObjType':
          # parent is complex data type
          lastParent.endOverride()
 
        if not newParentStack:
          # back at root - put root back in
          newParentStack.append(lastParent)
      
    # update local oldToNew map
    localOldToNew.update(oldToNew)

    # postprocess objects - set links now all objects are done
    if oldVersionStr is None:
      # copy subtree
      delayedLoadLinksCopy(localOldToNew, crossLinkData, 
                       ignoreMissing=ignoreMissing)
 
    else:
      # backwards compatibility.

      # first dereference links
      delayedLoadLinksComp(localOldToNew, crossLinkData)

      # minor post-processing
      from memops.format.compatibility.part1.Converters1 import minorPostProcess
      minorPostProcess(oldVersionStr, targetObj, delayDataDict, localOldToNew)
 
      # set TopObjects into TopObjects dictionary.
      newTopObjByGuid = newParent.__dict__['topObjects']
      guid = targetObj.__dict__['guid']
      if guid not in newTopObjByGuid:
        newTopObjByGuid[guid] = targetObj
      else:
        raise Exception("CCPN API error: %s: guid %s already in use"
                       % (targetObj, targetObj.__dict__['guid']))
 
    newTopObj = targetObj.topObject
    
    # set parent-to-child links
    mapping = globalMapping[newTopObj.metaclass.container.shortName]
    if not targetObj.isDeleted:
      # filter out deleted targetObj - could happen in minor postprocessing
      xmlImplementation.linkChildData(delayDataDict, targetObj, mapping,
                                      linkTopToParent=True)
  
  
  except:
    # try cleaning up
    import sys
    exc_info = sys.exc_info()
    # NB '[]' only put in for Python 2.1
    objsToBeDeleted = set([x for x in delayDataDict if not isinstance(x, int,)])
    deleteFailed = False
    for obj in objsToBeDeleted:
      try:
        obj._singleDelete(objsToBeDeleted)
      except:
        deleteFailed = True
        print ("WARNING Error in deleting object of class %s, id %s"
               % (xx.__class__, id(xx)))
  
    if targetObj is not None:
      try:
        topObj = targetObj.topObject
        topObj.__dict__['isReading'] = False
      except:
        deleteFailed = True
    
    if deleteFailed:
      print '''WARNING Error in clean-up of incorrectly copied data tree. 
      Data may be left in an illegal state'''
    else:
      print "NOTE created objects deleted without error"
    
    # re-raise original exception
    raise exc_info[0], exc_info[1], exc_info[2]
  
  # unset isReading and set to modified
  newTopObj.__dict__['isReading'] = False
  newTopObj.__dict__['isModified'] = True
    
 
  if oldVersionStr is None:
    # we are copying a subtree
    
    if targetObj is newTopObj:
      # root is a TopObject - we need to set isLoaded
      root = newTopObj.root
      newTopObj.__dict__['isLoaded'] = True
      guid = root.newGuid()
      newTopObj.__dict__['guid'] = guid
      root.__dict__['topObjects'][guid] = newTopObj
      
    #check validity
    targetObj.checkAllValid()
 
    # notify - list gicves you parent-before-child order
    for xx in notifyObjects:
      for notify in xx.__class__._notifies.get('__init__', ()):
        notify(xx)
    del notifyObjects
    
    
  else:
    # we are doing backwards compatibility.
    # Set to loaded and leave teh rest to others
    newTopObj.__dict__['isLoaded'] = True
 
  # clean up
  delayDataDict.clear()
  
  #
  return targetObj

def delayedLoadLinksComp(objectDict, linkData):
  """ Set links (of whatever kind) derefencing as you go using objectDict.
  Skips objects not found in the map.
  For backwards compatibility rather than compatibility 
  """
  
  popLinkData = linkData.pop
  getObj = objectDict.get

  try:
    while linkData:
      # setup
      curMap = popLinkData()
      val = popLinkData()
      obj = popLinkData()

      # map values
      valueList = list()
      for vv in val:
        oo = getObj(vv)
        if (oo is not None):
          valueList.append(oo)

      if valueList:
        name = curMap.get('name')
        hicard = curMap.get('hicard')
        # set element
        if (hicard == 1):
          ov = valueList[0]
        elif (hicard > 1):
          ov = valueList[:hicard]
        else:
          ov = valueList

        setattr(obj, name, ov)

  except:
    print 'Error during link dereferencing. Object was: ', obj
    print 'values were: ', val
    print 'tag name was: ', name
    raise


def delayedLoadLinksCopy(objectDict, linkData, ignoreMissing=False):
  """ Set links (of whatever kind) derefencing as you go using objectDict.
  For copySubTree rather than compatibility 
  """
  
  popLinkData = linkData.pop
  getObj = objectDict.get
   
  while linkData:
    
    # set up 
    curMap = popLinkData()
    val = popLinkData()
    obj = popLinkData()
    name = curMap['name']
    hicard = curMap['hicard']
    
    locard = curMap['locard']
    copyOverride = curMap.get('copyOverride')
    
    # NB copyOverride determines whether you are allowed to set links
    # that modifies object outside the subtree
    
    if hicard == 1:
    
      # get new value
      vv = val[0]
      newVal = getObj(vv)
      
      if newVal:
        # linked-to object replaced by new object
        setattr(obj, name, newVal)
      
      elif copyOverride:
        # try setting link to old object
        try:
          setattr(obj, name, vv)
          
        except:
          if ignoreMissing:
            # we can skip the link
            pass
          else:
            # link could not be handled
            raise Exception(
             "%s: Out-of-subtree -to-one link %s cannot be copied"
             % (obj, name)
            )
    
    else:
      # -to-many link
      
      # get new values
      others = []
      newies = []
      foundAll = True
      for vv in val:
        other = getObj(vv)
        if other:
          newies.append(other)
          others.append(other)
        else:
          foundAll = False
          others.append(vv)
      
      done = False
          
      if foundAll or copyOverride:
        # try making link to full set of objects
        try:
          setattr(obj, name, others)
          done = True
        except:
          pass
      
      if not done and not foundAll and ignoreMissing:
        # try making link to the subset of objects found in subtree
        try:
          setattr(obj, name, newies)
          done = True
        except:
          pass
          
      if not done:
        # link could not be handled
        raise Exception(
         "%s: Out-of-subtree link %s cannot be copied"
         % (obj, name)
        )
        

def downlinkTagsByImport(root):
  """ gives you the role names of links from MemopsRoot to TopObjects 
  in import order, so that imported packages come before importing packages
  """

  from memops.metamodel import Util as metaUtil

  leafPackages = []
  packages = [root.metaclass.container.topPackage()]
  for pp in packages:
    childPackages = pp.containedPackages
    if childPackages:
      packages.extend(childPackages)
    else:
      leafPackages.append(pp)

  # sort leafPackages by import (imported before importing)
  leafPackages = metaUtil.topologicalSortSubgraph(leafPackages,
                                                  'accessedPackages')
  tags = []
  for pp in leafPackages:
    cc = pp.topObjectClass
    if cc is not None:
      pr = cc.parentRole
      if pr is not None:
        tags.append(pr.otherRole.name)
  #
  return tags
  
  
def loadAllData(root):
  """ Load all data for a given root (version >= 2.0)
  """
  # load all new data before modifying IO map
  for tag in downlinkTagsByImport(root):
    for topObj in getattr(root, tag):
      if not topObj.isLoaded:
        topObj.load()


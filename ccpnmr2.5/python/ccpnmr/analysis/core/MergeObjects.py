LICENSE = """
======================COPYRIGHT/LICENSE START==========================

MergeObjects.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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
# A generic method to merge data model objects
# Transfer simple and link attributes from the source object to target object
# Does not transfer derived, automatic or immutable attributes
# Links will be transferred where possible
# Where necessary the Api is bypassed

# Logical analysis and design by R.H. Fogh

# Coding and testing by T.J. Stevens

#Definitions:
#Objects targetObj and sourceObj of class O
#link O.a (a) to class A, with backlink A.o ( o)

# A note on checks:
# Where the API is bypassed, the function does validity checks at each step,
# and rolls back the last step if the checks fail.
# The checks are done on souceObj, targetObj, objects on the other end of
# links, and the parents of the latter. The check on parents is done because
# this includes a check on the keys of the children - the merge cannot change
# the keys of either source or target, but can change the key of linked-to objects.

from memops.general import Constants as genConstants
from memops.general.Implementation import ApiError
from memops.universal import Util as uniUtil
from memops.metamodel import ImpConstants

def mergeObjects(sourceObj,targetObj):
  """Merges sourceObj into targetObj, deleting sourceObj.
  Attributes and links from sourceObj are added to targetObj
  provided 1) that they are not there already, and  
  2) that there is room.
  
  WARNING this function bypasses the API. 
  WARNING This function just might leave the data in an illegal state
  The function performs a number of checks for each individual change.
  If a check fails, the latest change is undone before the error exit,
  in an attempt to leave the data in a state that is legal. Note that only
  the latest change is undone - in case of error the data state will not be
  brought back to the state from before the execution of the command.
  Note that sourceObj is likely to be in an illegal state during execution,
  so that an error may well leave sourceObj in an illegal state. If this happens,
  deleting sourceObj may bring the data back to a legal state, and is unlikely 
  to cause further problems.
  In spite of the checks, some objects (not limited to sourceObj and targetObj)
  may be left in an illegal state, even if no error is raised.
  It is recommended to use this function with caution, 
  and to run checkAllValid after it has been used.  """

  #exclude parent # done! Mandatory
  
  #same class check
  if sourceObj.qualifiedName != targetObj.qualifiedName:
    return

  #ATTRIBUTES:
  #Objects targetObj, sourceObj, with attribute a
  
  objClass = targetObj.metaclass

  for a in objClass.getAllAttributes():
    attrName = a.name
 
    if a.isDerived or a.isAutomatic or a.changeability == ImpConstants.frozen:
      continue
      
    elif a.hicard == a.locard:
      continue
 
    elif a.hicard == 1:
      if targetObj.__dict__[attrName] is None:
        setattr(targetObj,attrName,sourceObj.__dict__[attrName])

    else:
      # find add operation
      addfunc = getattr(targetObj, 'add' + uniUtil.upperFirst(a.baseName))
      
      # Thisd is OK, as we iterate over list2 but modify list1
      attrList1 = targetObj.__dict__[attrName]
      attrList2 = sourceObj.__dict__[attrName]
      
      if a.isUnique:
        # no duplicates = might be list or set
        if a.hicard > 1:
          nSpaces = max(0, a.hicard - len(attrList1))
        else:
          nSpaces = -1
          
        for aVal in attrList2:
          if nSpaces == 0:
            break
          else:
            nSpaces -= 1
            addfunc(aVal)
            
      else:
        # might have duplicates (and must be an internal list)
        for aVal in attrList2:
          if len(attrList1) >= a.hicard and a.hicard != genConstants.infinity:
            break
          # keep adding while there is room
          if attrList1.count(aVal) < attrList2.count(aVal):
            addfunc(aVal)
 
  #LINKS:
  niceLinks = []
  nastyLinks = []
  childLinks = []
  for a in objClass.getAllRoles():
    linkName = a.name
    
    # select links and how to treat them
    
    if a.hicard == a.locard or a.isDerived or a.isAutomatic:
      continue
 
    if a.changeability == ImpConstants.frozen:
      continue
      # This is probably right; it could be changed if we bypassed the API.
      
    o = a.otherRole
    
    if a.hierarchy == ImpConstants.child_hierarchy:
      childLinks.append(a)
    
    elif o is None or o.changeability != ImpConstants.frozen:
      # links that can be handled without bypassing API
      niceLinks.append(a)
    else:
      # links that require bypassing API
      nastyLinks.append(a)
  
  for a in niceLinks:
    # links that can be handled without bypassing API
 
    linkName = a.name
    o = a.otherRole
    if o is not None:
      backName = o.name
      #print linkName, a.locard, a.hicard, o.locard, o.hicard  
      
    if o is None or o.hicard != o.locard:
      #
      #print "C3", linkName
      #
      # NB this does NOT break API
      #
      # We are setting/adding/removing from the .a side.
      # if o is None there will be no problems.
      # If o.hicard == 1, attrObj.o can be overwritten
      # Otherwise, as o.hicard != o.locard it will always be possible either to
      # remove sourceObj from attrObj or to add targetObj to attrObj 
      # regardless of the exact cardinalities and of len(attrObj.o)

      if a.hicard == 1:
        #assert a.locard == 0
        #
        # there will be no problems on the source/target side as we only
        # make changes when the link is unset in the target
        
        # Whatever the number of objects, you will eitehr be able to
        # remove sourceObj from attrObj or to add targetObj.
        #
        if getattr(targetObj,linkName) is None:
          attrObj = getattr(sourceObj,linkName)
          try:
            setattr(sourceObj,linkName,None)
          except:
            pass
          setattr(targetObj,linkName,attrObj)
      
      else:
        # There will be no problems on the attrObj side (see above).
        # On the source/target side we will get the desired result as
        # the API simply passes if you try to add an existing object,
        # It also deletes the old link to sourceObj where appropriate
        
        # find add operation
        ss = uniUtil.upperFirst(a.baseName)
        addfunc = getattr(targetObj, 'add' + ss)
        removefunc = getattr(sourceObj, 'remove' + ss)
        
        # NB we cannot use teh raw list as we modify it during the loop
        for attrObj in getattr(sourceObj,linkName):
        
          try:
            removefunc(attrObj)
          except ApiError:
            pass
            #print 'Failed to remove %s for %s' % (linkName,sourceObj.className)
          
          
          try:
            # Adds objects to targetObj.a as long as there is room
            # (a.hicard could be e.g. 2)
            addfunc(attrObj)
          except ApiError:
            pass
            #print 'Failed to add %s for %s' % (linkName,targetObj.className)
            break

    elif o.hicard == 1 and o.locard == 1:
      
      if a.hicard == 1:
        #assert a.locard == 0
        oldVal = getattr(targetObj,linkName)
        if oldVal is None:
          newVal = getattr(sourceObj,linkName)
          setattr(newVal, backName, targetObj)
      
      else:
        # assert a.hicard != a.locard
        # asser a.hicard != 1
        for attrObj in getattr(sourceObj,linkName):
          try:
            setattr(attrObj, backName, targetObj)
          except ApiError:
            pass
      
    else:
      #
      #print "C4", linkNam
      #
      # NB this does NOT break API
      #
      # we know that o.hicard == o.locard > 1
      # The trick is that since o.hicard == o.locard > 1 and a.changeability != frozen,
      # it must be possible to set attrObj.o to an appropriate tuple without
      # getting into trouble.
      if a.hicard == 1:
        #assert a.locard == 0
        attrObj = getattr(sourceObj,linkName)
        linkList = list(getattr(attrObj,backName))
        linkList[linkList.index(sourceObj)] = targetObj
        setattr(attrObj,backName,linkList)
        
      else:
        # assert a.hicard != a.locard
        # asser a.hicard != 1
        for attrObj in getattr(sourceObj,linkName):
          linkList = list(getattr(attrObj,backName))
          linkList[linkList.index(sourceObj)] = targetObj
          setattr(attrObj, backName, linkList)
   
  # make sure we are valid before going into the tough part
  targetObj.checkValid()
  
  
  if nastyLinks:
    root = targetObj.root
    try:
      root.override = True 
     
      for a in nastyLinks:
        # links that can *NOT* be handled without bypassing API
 
        linkName = a.name
        o = a.otherRole
        backName = o.name
        keyNames = o.container.keyNames
        attrObjClass = a.valueType
        downlink = attrObjClass.parentRole.otherRole.name
        #print linkName, a.locard, a.hicard, o.locard, o.hicard
      
        if a.hicard == 1:
          #print "C1", linkName

          if  getattr(targetObj,linkName) is None:

            #do
            attrObj = getattr(sourceObj,linkName)
            if backName in keyNames:
              oldKey = attrObj.getLocalKey()
            setattr(sourceObj, linkName, None)
            setattr(targetObj, linkName, attrObj)
            if backName in keyNames:
              newKey = attrObj.getLocalKey()
              # this changes key for attrObj - fix it.
              childDict = attrObj.parent.__dict__[downlink]
              if newKey in childDict:
                # key already taken - undo
                setattr(targetObj, linkName, None)
                setattr(sourceObj, linkName, attrObj)
                raise ApiError("Merge failure: %s key %s already in use"
                               % (attrObj.qualifiedName(), newKey))
              else:
                del childDict[oldKey]
                childDict[newKey] = attrObj

            # test
            try:
              attrObj.checkValid()
              targetObj.checkValid()

            # undo
            except:
              setattr(targetObj, linkName, None)
              setattr(sourceObj, linkName, attrObj)
              if backName in keyNames:
                del childDict[newKey]
                childDict[oldKey] = attrObj
              print ("Merge failure: %s, %s result is not valid"
                     % (targetObj, attrObj))
              raise

        else:
          #
          # assert a.hicard != 1
          #
          # NB if a.locard > 0 the code below could create an illegal
          # sourceObj. Which would not be a problem if all went well,
          # but would render the final state illegal if the merge ran into an 
          # error somewhere else later
          # We ignore this as links that are locard>0 in one direction and
          # frozen in the other direction would make both objects impossible
          # to create except under override conditions. The problem is *very*
          # unlikely ever to arise.
          
          #
          #print "C2", linkName
          # set up
          keepList = list(getattr(targetObj, linkName))
          ll = list(getattr(sourceObj, linkName))
          
          if a.hicard == genConstants.infinity:
            moveList = ll
            ignoreList = []
          else:
            nSpaces = a.hicard - len(keepList)
            if nSpaces > 0:
              moveList = ll[:nSpaces]
              ignoreList = ll[nSpaces:]
            else:
              continue
            
          # do
          if backName in keyNames:
            oldKeys = [x.getLocalKey() for x in moveList]
          setattr(sourceObj, linkName, ignoreList)
          setattr(targetObj, linkName, keepList + moveList)
          if backName in keyNames:
            newKeys = []
            for ii, attrObj in enumerate(moveList):
              childDict = attrObj.parent.__dict__[downlink]
              newKey = attrObj.getLocalKey()
              if newKey in childDict:
                # key already taken - undo
                setattr(targetObj, linkName, None)
                setattr(sourceObj, linkName, attrObj)
                for jj, nk in enumerate(newKeys):
                  ao = moveList[jj]
                  cd = ao.parent.__dict__[downlink]
                  cd[oldKeys[jj]] = ao
                  del cd[nk]
                raise ApiError("Merge failure: %s key %s already in use"
                               % (attrObj.qualifiedName(), newKey))
              else:
                newKeys.append(newKey)
                # del childDict[oldKey] 
                # TJS edit: to be checked
                del childDict[oldKeys[ii]]
                childDict[newKey] = attrObj

          # test
          try:
            targetObj.checkValid()
            for attrObj in moveList:
              attrObj.checkValid()

          # undo
          except:
            setattr(targetObj, linkName, keepList)
            setattr(sourceObj, linkName, moveList + ignoreList)
            if backName in keyNames:
              for jj, nk in enumerate(newKeys):
                ao = moveList[jj]
                cd = ao.parent.__dict__[downlink]
                cd[oldKeys[jj]] = ao
                del cd[nk]
            raise
          
 
    finally:
      root.override = False 
  
  
  if childLinks:
    # now move children. This is a full bypass, no overrides
    for a in childLinks:
      parentName = a.otherRole.name
      sourceDd = sourceObj.__dict__[a.name]
      targetDd = targetObj.__dict__[a.name]
      topObj = targetObj.topObject
      
      if a.hicard == 1:
        # single kid (rare case)
        targetObj.__dict__[a.name] = oo = sourceObj.__dict__[a.name]
        sourceObj.__dict__[a.name] = None
        oo.__dict__[parentName] = targetObj
        oo.__dict__['topObject'] = topObj
 
      elif a.valueType.keyNames == ['serial']:
        # multiple kid with serial key
        nextSerial = targetObj.__dict__['_serialDict'][a.name] + 1
        for junk, oo in sorted(sourceDd.items()):
          targetDd[nextSerial] = oo
          del sourceDd[nextSerial]
          oo.__dict__[parentName] = targetObj
          oo.__dict__['topObject'] = topObj
        targetObj.__dict__['_serialDict'][a.name] = nextSerial
 
 
      else:
        # multiple kid with normal key
        for localKey,oo in sorted(sourceDd.items()):
          if localKey in targetDd:
            # key is taken - skip object
            continue
          else:
            targetDd[localKey] = oo
            del sourceDd[localKey]
            oo.__dict__[parentName] = targetObj
            oo.__dict__['topObject'] = topObj
      
  targetObj.checkValid()
  #print "S1", sourceObj
  sourceObj.delete()
  #print "S2", sourceObj
  #print "T1", targetObj
  return targetObj
   

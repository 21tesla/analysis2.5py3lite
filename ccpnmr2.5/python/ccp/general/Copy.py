"""
======================COPYRIGHT/LICENSE START==========================

Copy.py: Copy functions for CCPN objects

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

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
- PDBe website (http://www.ebi.ac.uk/pdbe/)

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

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================
"""
import traceback, sys


import memops.api.Implementation as Implementation

def isValidMetaAttr(metaAttr):
  
  isValid = True
  
  # Don't use this one... TODO CHECK ALL BELOW WITH RASMUS!!!
  if metaAttr.name == 'inConstructor':
    isValid = False

  # Ignore derived attributes and automatic ones
  if metaAttr.isDerived == True or metaAttr.isAutomatic == True or metaAttr.isImplementation or metaAttr.isImplicit:
    isValid = False
    
  return isValid

def copyAttributeInfo(fromInstance,toInstance,overWrite = 1,verbose = True):

  #
  # Copy attribute information
  #
  
  instanceMetaClass = fromInstance.metaclass

  for metaAttr in instanceMetaClass.getAllAttributes():
   
    attrName = metaAttr.name
    
    if not isValidMetaAttr(metaAttr):
      continue
    
    # Get value...
    fromInstanceAttrValue = getattr(fromInstance,attrName)
   
    if not fromInstanceAttrValue:
      # Ignore empty values
      continue

    else:
      # Get other value
      toInstanceAttrValue = getattr(toInstance,attrName)
      
      # Ignore if not overwrite and value set
      if not overWrite and toInstanceAttrValue:
        continue
      
      # Ignore if values are the same
      if fromInstanceAttrValue == toInstanceAttrValue:
        continue

      # If fromValue not empty and the values are not the same, try to
      # set the value in toInstance
      
      # Is a single value
      if metaAttr.hicard == 1:
      
        # Value already set and can't be added in toInstance... print warning
        if toInstanceAttrValue != None:
          if verbose:
            print "    Warning: can't relink attribute/role %s value %s: already set" % (attrName,fromInstanceAttrValue)
        
        else:
          setattr(toInstance,attrName,fromInstanceAttrValue)
          if verbose:
            print "    Copied attribute %s (now %s, was %s)." % (attrName,fromInstanceAttrValue,toInstanceAttrValue)
          
      # Is a list - just add the value (if not already present!)
      else:
        for fromInstanceAttrValueElement in fromInstanceAttrValue:
          if not fromInstanceAttrValueElement in toInstanceAttrValue:
            toInstanceAttrValue += (fromInstanceAttrValueElement,)
            setattr(toInstance,attrName,toInstanceAttrValue)
            if verbose:
              print "    Added %s to attribute %s." % (fromInstanceAttrValueElement,attrName)

def getAttributeParameterDict(instance):

  parameters = {}

  instanceMetaClass = instance.metaclass

  for metaAttr in instanceMetaClass.getAllAttributes():
   
    attrName = metaAttr.name

    # Ignore derived attributes and automatic ones
    if not isValidMetaAttr(metaAttr):
      continue
    
    # Get value...
    instanceAttrValue = getattr(instance,attrName)
    parameters[attrName] = instanceAttrValue

  return parameters
  
def nmrShiftSelect(shift,otherShift,toResName):
  
  #
  # Dummy function - should be set on application level with
  # correct user interface
  #
  
  pass
  
def copyResonanceInfo(fromRes,toRes,nmrShiftSelect = nmrShiftSelect,deleteOldInfo = False,toResName = None,copyTrack = None, verbose = True):

  if fromRes.className != toRes.className or fromRes.className not in ['Resonance','FixedResonance']:
    if verbose:
      print "  Error: CopyResonanceInfo only for resonances. Aborting."
    return
  
  copyAttributeInfo(fromRes,toRes,verbose = verbose)
  
  #
  # Copy link info... (WARNING: are returned as [], not as () as for attributes!)
  #
  # Ignoring links to resonanceSet, fixedResonanceSet:
  #   should be set in existing resonance (assuming that this one not assigned yet)
  #
  # Parent of both classes assumed to be same project
  #
  # Other classes should all be 0...* links
  #
  # Treating cases differently according to link back to resonance from linked class
  #
  
  ignoreLinkNames = ['resonanceSet','root']
  
  resMetaClass = fromRes.metaclass

  for metaLink in resMetaClass.getAllRoles():
   
    linkName = metaLink.name

    # Ignore derived attributes and automatic ones
    if metaLink.isDerived == True or metaLink.isAutomatic == True or linkName in ignoreLinkNames:
      continue
    
    # Get value...
    linkedClass = getattr(fromRes,linkName)
    
    if not linkedClass:
      # Ignore empty values
      continue
    
    else:
    
      if type(linkedClass) == type(set()) or type(linkedClass) == type(frozenset()):
        linkedClassList = list(linkedClass)
        refLinkedClass = linkedClassList[0]
      
      else:
        refLinkedClass = linkedClass
        linkedClassList = [linkedClass]
      
      # Ignore parent link
      if refLinkedClass == fromRes.parent:
        continue

      linkedMetaClass = refLinkedClass.metaclass
      linkedClassName = refLinkedClass.metaclass.name   # Could also use className
      linkedClassParentRole = refLinkedClass.metaclass.parentRole.name

      #
      # All others are linked to resonance 'normally' 
      #
      
      isChild = 0
      
      resonanceRole = linkedMetaClass.getRole(name = 'resonance')
      resonancesRole = linkedMetaClass.getRole(name = 'resonances')
      
      if resonanceRole:
        linkedMetaRole_hicard = resonanceRole.hicard
        linkedMetaRole_locard = resonanceRole.locard
        
        if resonanceRole.changeability == 'changeable':
          linkedMetaRole_isChangeable = True
        else:
          linkedMetaRole_isChangeable = False
          
      elif resonancesRole:
        linkedMetaRole_hicard = resonancesRole.hicard
        linkedMetaRole_locard = resonancesRole.locard
        
        if resonancesRole.changeability == 'changeable':
          linkedMetaRole_isChangeable = True
        else:
          linkedMetaRole_isChangeable = False
          
      else:
        # fromRes is parent
        linkedMetaRole_hicard = 1
        linkedMetaRole_locard = 1
        linkedMetaRole_isChangeable = False
        isChild = 1
      
      if verbose:
        print "    Copying %s..." % linkName

      #
      # Loop over all linked classes for this linkName
      #

      for linkedClass in linkedClassList:
        
        #
        # First check if class is not already deleted
        # (this occurs in cases where there's a distConstrItem with
        # the same resonance for both resonances for example)
        #
        
        if linkedClass.isDeleted == True:
          continue
          
          
        #
        # Get the parent class for this linkedClass
        #
        
        linkedClassParent = getattr(linkedClass,linkedClassParentRole)

        #
        # If linkedClass has toRes as parent... special case
        # (only really valid for ApplicationData in this case!)
        #

        if isChild == 1:
          if linkedClassName[0:7] != 'AppData':
            if verbose:
              print "      Can currently only handle applicationData children, not %s!" % linkedClassName
            continue
            
          parameters = getAttributeParameterDict(linkedClass)
        
          applData = toRes.findFirstApplicationData(**parameters)
          
          if not applData and parameters != {}:
          
            try:
              newLinkedClass = linkedClass.__class__(toRes, **parameters)
              
            except:
              if verbose:
                print "      Problems copying %s:" % linkedClassName
                print "        %s" % traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
              
          else:
            # TODO: should really only check if key attributes are the same...
            pass

        #
        # Link back to resonance if has hicard 1
        #
 
        elif linkedMetaRole_hicard == 1:
 
          #
          # If changeable, then reset resonance from same object
          #
 
          if linkedMetaRole_isChangeable == True:
            linkedClass.resonance = toRes

          #
          # If not changeable, then make new class with link to fromRes,
          # copy attributes from the old class to new one
          #
          # Assuming that only link from other class is to resonance and its parent
          #
 
          else:
            parameters = getAttributeParameterDict(linkedClass)
            
            if linkedClassName[:14] == 'PeakDimContrib':
              if linkedClass.peakContribs:
                parameters['peakContribs'] = linkedClass.peakContribs
            
            try:
              newLinkedClass = linkedClass.__class__(linkedClassParent, resonance = toRes, **parameters)
              # Delete old class if it's a 'move' (not 'copy')
              if deleteOldInfo == True:
                linkedClass.delete()
 
              newLinkedClass.checkValid()
              
            except:
              #
              # TODO: user interaction to pick 'correct' shift?
              #
              if verbose:
                print "      Problems creating new class %s:" % linkedClassName
                print "        %s" % traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
              
              #
              # Exception for shifts: take average value of shifts and report!
              #
              
              if linkedClass.className == 'Shift':
                
                otherShift = linkedClassParent.findFirstMeasurement(resonance = toRes)
                
                nmrShiftSelect(linkedClass,otherShift,toResName)
                                                
        #
        # Link back to multiple resonances
        #
 
        else:

          #
          # Setup a list of resonances with the toRes instead of fromRes
          #
          
          doCopy = 1
          
          resonanceList = list(linkedClass.resonances)
            
          resonanceIndex = resonanceList.index(fromRes)
          if toRes not in resonanceList:
            resonanceList[resonanceIndex] = toRes
          else:
            if linkedMetaRole_hicard == linkedMetaRole_locard:
              doCopy = 0
            else:
              resonanceList.pop(resonanceIndex)
              
          #
          # Same resonance occurs twice in list with locard = hicard - this is not allowed
          #
          
          if not doCopy:
          
            #
            # Delete old class if it's a 'move'... cannot be used if same resonance occurs twice
            #

            if deleteOldInfo == True:
              linkedClass.delete()

 
          #
          # If changeable, then reset resonance from same object
          # Take care not to set same resonance twice!
          #
 
          elif linkedMetaRole_isChangeable == True:
            linkedClass.resonances = resonanceList

          #
          # If not changeable, then make new class with link to same resonances
          # except toRes replaced by fromRes,
          # copy attributes from the old class to new one
          #
          # Assuming that only link from other class is to resonance and its parent
          #
          # TODO: will crash if it needs some value beside the parent!
          #
 
          else:
          
            parameters = getAttributeParameterDict(linkedClass)
              
            try:
              newLinkedClass = linkedClass.__class__(linkedClassParent,resonances = resonanceList, **parameters)
            except:
              if verbose:
                print "      Problems creating %s:" % linkedClass.className
                print "        %s" % traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
              continue

            #
            # Exception for dihedralConstraints: have items linked to them as well!
            # Have to copy these separately...
            #

            if linkedClassName == 'DihedralConstraint':

              if verbose:
                print  "      Copying dihedralConstraintItems..."

              for dihedralConstraintItem in linkedClass.items:

                parameters = getAttributeParameterDict(dihedralConstraintItem)
                newDihedralConstraintItem = dihedralConstraintItem.__class__(newLinkedClass,**parameters)

            else:
              
              #
              # Make sure to retain order for PairwiseConstraintItems
              #
              
              if linkedClassName in ('DistanceConstraintItem', 'HBondConstraintItem', 'JCouplingConstraintItem', 'RdcConstraintItem'):
                firstResonance = linkedClass.firstResonance
                if firstResonance == fromRes:
                  newLinkedClass.firstResonance = toRes
                else:
                  newLinkedClass.firstResonance = firstResonance

              #
              # Tag copied items if necessary (not for dihedrals)
              #

              if copyTrack:

                resonanceList = list(newLinkedClass.resonances)
                resonanceIndex = resonanceList.index(toRes)

                newLinkedClass.addApplicationData(Implementation.AppDataInt(application = 'ccp',keyword = 'copied',value = resonanceIndex))

            #
            # Also need to copy application data...
            #

            for applData in linkedClass.applicationData:

              parameters = getAttributeParameterDict(applData)

              try:
                newApplData = newLinkedClass.findFirstApplicationData(**parameters)
                if not newApplData:
                  newApplData = applData.__class__(**parameters)
                  newLinkedClass.addApplicationData(newApplData)

              except:
                if verbose:
                  print "      Problems copying %s:" % applData.className
                  print "        %s" % traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]


            #
            # Delete old class if it's a 'move' (not 'copy')
            #

            if deleteOldInfo == True:
              linkedClass.delete()

            newLinkedClass.checkValid()

            #except:

            #  print "      Problems creating new class %s:" % linkedClassName
            #  print "        %s" % traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
             
  # Check if valid!
  toRes.checkValid()
 
#
# Code to make a copy of a constraint list, based on the original file.
# Directly changes information by bypassing the API - use with care
#

import re, os, shutil

def createNewNmrConstraintStoreFromCopy(project,
                                        nmrProject,
                                        ncs,
                                        structureEnsemble,
                                        entry = None,
                                        newStrucAnalName = 'New structure analysis',
                                        newStrucGenName = None # If this is set, will make structure generation instead of analysis
                                        ):
  
  """
  
  Makes a copy of a constraint list based on the original file. Use with care, as this
  code bypasses the API.
  
  Input:
  
   args:
    project:      A CCPN project
    nmrProject:   A CCPN Nmr.NmrProject
    ncs:          A CCPN Nmr.NmrConstraintStore on which the new one will be based
    structureEnsemble:   A MolStructure.StructureEnsemble
   
   kywds:
    entry:        A CCPN NmrEntry.Entry (the new structure generation will be linked to this)                  
    newStrucAnalName:   Name of the new structure analysis associated with this constraint store (default)
    newStrucGenName:    Name of the new structure generation associated with this constraint store.
  
  Output:
  
    A new Nmr.NmrConstraintStore object
  
  """
    
  #
  # Make a new constraint list set based on the old one - this one will be recalibrated
  #
  
  repository = ncs.findFirstActiveRepository()
  package = 'ccp.nmr.NmrConstraint'
  subDir = repository.getFileLocation(package)
  
  fileName = "%d+%s.xml" % (ncs.serial,ncs.guid)
  
  origConstraintsFile = os.path.join(subDir,fileName)
  
  newNmrConstraintStore = project.newNmrConstraintStore(nmrProject = nmrProject)
  
  newConstraintsFileName = '%d+%s.xml' % (newNmrConstraintStore.serial,newNmrConstraintStore.guid)
  newConstraintsFile = os.path.join(subDir,newConstraintsFileName)
  
  shutil.copy(origConstraintsFile,newConstraintsFile)
  
  resetNmrConstraintStoreFile(newNmrConstraintStore,newConstraintsFile)
  
  #newNmrConstraintStore.structureAnalyses = []
  newNmrConstraintStore.structureGeneration = None
  
  for strucAnal in newNmrConstraintStore.structureAnalyses:
    if strucAnal.nmrConstraintStore == newNmrConstraintStore:
      strucAnal.nmrConstraintStore = ncs
  
  creationKeywds = {'nmrConstraintStore': newNmrConstraintStore}
  
  if newStrucGenName:
    creationKeywds['name'] = newStrucGenName
    creationKeywds['structureEnsemble'] = structureEnsemble
    creationClass = nmrProject.newStructureGeneration
  else:
    creationKeywds['name'] = newStrucAnalName
    creationKeywds['structureEnsembles'] = [structureEnsemble]
    creationClass = nmrProject.newStructureAnalysis
  
  strucGenOrAnal = creationClass(**creationKeywds)

  if entry:
    strucGenOrAnal.addEntry(entry)
  
  # Make sure all is OK
  try:
    newNmrConstraintStore.constraintLists
  except:
    os.remove(newConstraintsFile)
    newNmrConstraintStore = None
  
  return newNmrConstraintStore

def resetNmrConstraintStoreFile(newNmrConstraintStore,newConstraintsFile):
  
  #
  # NOTE: these are 'expert' Data Model hacks - DO NOT USE AS EXAMPLE!!
  #
  
  newGuid = newNmrConstraintStore.guid
  newCreatedBy = newNmrConstraintStore.createdBy
  newSerial = newNmrConstraintStore.serial
  
  newNmrConstraintStore.__dict__['isLoaded'] = False
  newNmrConstraintStore.__dict__['isModified'] = False
  newNmrConstraintStore.__dict__['isStored'] = True
  
  cfile = open(newConstraintsFile)
  
  cfileLines = cfile.readlines()
   
  guidPatt = re.compile('guid\="([^"]+)"')
  userPatt = re.compile('createdBy\="([^"]+)"')
  serialPatt = re.compile('(serial="\d+")')
  
  strucGenPatt = re.compile("structureGenerations")
  strucAnalPatt = re.compile("structureAnalyses")
    
  pattFound = 0
  lineNum = 0
  
  while lineNum < len(cfileLines):
    line = cfileLines[lineNum]
    
    guidSearch = guidPatt.search(line)
    if guidSearch:
      line = line.replace(guidSearch.group(1),newGuid)
      cfileLines[lineNum] = line
      pattFound += 1
    
    userSearch = userPatt.search(line)
    if userSearch:
      line = line.replace(userSearch.group(1),newCreatedBy)
      cfileLines[lineNum] = line
      pattFound += 1

    serialSearch = serialPatt.search(line)
    if serialSearch:
      line = line.replace(serialSearch.group(1),'serial="%d"' % newSerial)
      cfileLines[lineNum] = line
      pattFound += 1
    
    """
    for strucPatt in (strucGenPatt,strucAnalPatt):
      strucSearch = strucPatt.search(line)
      if strucSearch:
        strucSearch = None
        cfileLines.pop(lineNum)
        print line,
        line = cfileLines[lineNum]
        while not strucSearch:
          line = cfileLines[lineNum]
          strucSearch = strucPatt.search(line)
          print line,
          cfileLines.pop(lineNum)
    """
    
    lineNum += 1
    
    if pattFound == 3:
      break
  
  cfile.close()
       
  cfile = open(newConstraintsFile,'w')
  for line in cfileLines:
    cfile.write(line)
  cfile.close()

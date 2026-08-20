"""
======================COPYRIGHT/LICENSE START==========================

ChemCompEditor.py: Allows editing of chemComp information.

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in license/LGPL.license
 
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
import os, string, time

import Tkinter


from memops.universal.Io import joinPath

from memops.universal.Util import returnInt, returnFloat

from memops.gui.BasePopup import BasePopup
from memops.gui.ButtonList import ButtonList
from memops.gui.Frame import Frame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Entry import Entry
from memops.gui.Label import Label
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.Util import createDismissButton
from memops.gui.ScrolledListbox import ScrolledListbox

from memops.gui.MessageReporter import showError, showYesNo, showInfo

from ccp.general.Util import getDescriptorDict, mergeDescriptorDicts, getDescriptorFromDict

from ccp.api.molecule import ChemComp
from memops.api import Implementation

from memops.metamodel import ImpConstants

from memops.gui.MultiSelectionList import MultiSelectionListPopup
  
# CHECK IF something similar ALREADY EXISTS!

"""
def getCcpnValue(dataType,value):

  if value == 'None':
    value = None
    
  # Can also use dataType.name: gives more specific info (e.g. Word, ...)
  if dataType.typeCode == 'String':
    value = value
    
  elif dataType.typeCode == 'Boolean':
    if value == '0':
      value = 0
    if value:
      value = 1
      
  else:
    raise "ERROR: %s not handled" % dataType.typeCode

  return value
"""

#
# TODO: need a text popup window that lists info on the current chemcomp (e.g.
# list a chemAtom, all the bonds, angles, torsions that it's in...).
#

def setCcpnObjectValue(guiParent,ccpnObject,metaAttrOrRole,value,widget,updateAfterSet,checkParent = True, checkValid = True, attrIndex = None, attrName = None, ignoreList = False):

  setValue = 1

  if not attrName:
    attrName = metaAttrOrRole.name
    
  origValue = getattr(ccpnObject,attrName)
  
  if type(origValue) == type(('',)):
    origValue = list(origValue)
  
  # Handle atom_1, atom_2, .. type stuff
  if attrIndex != None:
    tvalue = origValue[:]
    if attrIndex < len(tvalue):
      tvalue.pop(attrIndex)
      tvalue.insert(attrIndex,value)
    else:
      tvalue.append(value)
    value = tvalue

  try:

    if not ignoreList and metaAttrOrRole.hicard > 1 or metaAttrOrRole.hicard == -1:
    
      try:
        if not type(value) in [type([]),type(('',))]:
          value = list(eval(value))        

        valueIndex = 0
        while valueIndex < len(value):
          if value[valueIndex] != None:
            valueIndex += 1
          else:
            value.pop(valueIndex)

      except:
        showError('Invalid list form', "Value '%s' is not a python list/tuple!" % value, parent = guiParent)
        widget.set(str(origValue))
        setValue = 0

    if setValue and value != origValue:
    
      if metaAttrOrRole.changeability != ImpConstants.changeable:

        makeChange = True

        if checkValid:
          makeChange = showYesNo("Confirm change","This a non-changeable attribute. Really change?",parent = guiParent)

        if makeChange:
          # THIS IS DANGEROUS! Immediately run checkValid...
          
          ccpnObject.__dict__[attrName] = value
          
          try:
            if checkValid:
              if checkParent:
                ccpnObject.parent.checkAllValid()
              else:
                ccpnObject.checkAllValid()
                
          except Implementation.ApiError, e:
            ccpnObject.__dict__[attrName] = origValue
            showError('Cannot set non-changeable value', e.error_msg, parent = guiParent)

      else:
        
        try:
        
          setattr(ccpnObject,attrName,value)
          
          if checkValid:
            ccpnObject.checkAllValid()
            
        except Implementation.ApiError, e:
          setattr(ccpnObject,attrName,origValue)
          showError('Problems setting value', e.error_msg, parent = guiParent)
    
    if value != origValue:
      updateAfterSet()

  except Implementation.ApiError, e:
    showError('Setting value', e.error_msg, parent = guiParent)


def makeObjectKey(object,currentKeys):

  objectKey = "%s:" % object.className
  
  #
  # Hacks to give specific classes a good key...
  #
  
  if object.className == 'ChemAtomSetVar':
    currentKeys = list(currentKeys)
    currentKeys.append(('chemAtoms',object.chemAtoms))
  
  #
  # Now generic setup
  #
  
  for (attrName,value) in currentKeys:
    curAttr = getattr(object,attrName)

    if type(curAttr) in [type([]),type(('',))]:
      attrString = ""
      for curAttrItem in curAttr:
        if isinstance(curAttrItem,ChemComp.ChemAtom) or isinstance(curAttrItem,ChemComp.LinkAtom):
          curAttrString = curAttrItem.name
        else:
          curAttrString = str(curAttrItem)

        attrString += "%s, " % curAttrString

      attrString = attrString[:-2]

    else:
      if isinstance(curAttr,ChemComp.ChemAtom) or isinstance(curAttr,ChemComp.LinkAtom):
        attrString = curAttr.name
      else:
        attrString = str(curAttr)

    objectKey += "[%s], " % attrString

  objectKey = objectKey[:-2]

  return objectKey

class MultiRowScrolledMatrix(ScrolledMatrix):

  def getEditWidget(self,colNum):
    
    if hasattr(self.editWidgets[colNum],'widget'):
      editWidget = self.editWidgets[colNum].getWidget()
    else:
      editWidget = self.editWidgets[colNum]
    
    return editWidget

# relies on data being stored in molType/ccpCode.xml

class MultiWidget:

  def __init__(self,window,setFunc,parDict):
  
    self.window = window
    self.widget = None
    self.setFunc = setFunc
    self.chemComp = parDict['chemComp']
  
  def createPulldownMenu(self, widgetList, widgetDict = None):
  
    self.widgetList = widgetList
    self.widgetDict = widgetDict
    self.widget = PulldownMenu(self.window, entries = self.widgetList, callback = self.setFunc, width=25, do_initial_callback = False, force_callback = True)
    self.widgetType = 'PulldownMenu'    

  def createScrolledList(self, widgetList, widgetDict):
  
    self.widgetList = widgetList
    self.widgetDict = widgetDict
    self.widget = ScrolledListbox(self.window,width = 50,height = 5,
                                              selectmode = Tkinter.MULTIPLE,
                                              initial_list = self.widgetList)
    self.widgetType = 'ScrolledList'    

  def createValueList(self, widgetList):
  
    self.widgetList = widgetList
    self.widgetDict = {}
    self.widget = ScrolledListbox(self.window,width = 50,height = 5,selectmode = Tkinter.MULTIPLE,initial_list = widgetList,addDeleteButtons = True)
    self.widgetType = 'ValueList'    
  
  def createEntry(self):
  
    self.widgetList = []
    self.widgetDict = {}
    self.widget = Entry(self.window, returnCallback=self.setFunc, width=25)
    self.widgetType = 'Entry'
  
  def getWidget(self):
    
    #print "WidgetCreate"
    
    # DON'T FORGET: also put info for setter, getter!!
    
    self.widgetList = None
    self.widget = None
  
    attrName = self.window.getCurrentAttributeName()
    
    if self.window.metaRoleOrder[attrName]:
      (attrName,attrIndex) = self.window.metaRoleOrder[attrName]
      hasOrder = True
    else:
      hasOrder = False
      
    metaAttrOrRole = self.window.metaAttrOrRoleDict[attrName]
    
    dataType = None
    if hasattr(metaAttrOrRole,'dataType'):
      dataType = metaAttrOrRole.dataType
    
    # MAKE SURE LISTS IN PULLDOWNS ARE LOCALLY TRACKED!
    if dataType:
    
      if metaAttrOrRole.locard == 0:
        addKeyList = [None]
      else:
        addKeyList = []
        
      #
      # Different types...
      #

      if dataType.typeCode == 'Boolean':
        self.createPulldownMenu(['1','0','None']) 
    
      # TODO Can also use dataType.name: gives more specific info (e.g. Word, ...)
      elif dataType.typeCode == 'String':
        
        if metaAttrOrRole.hicard > 1 or metaAttrOrRole.hicard == -1:

          if dataType.enumeration:
            self.createScrolledList(list(dataType.enumeration) + addKeyList, {})
          
          else:
            attrList = list(getattr(self.window.ccpnObject,attrName))
            self.createValueList(attrList + addKeyList)
      
        else:

          if dataType.enumeration:
            self.createPulldownMenu(list(dataType.enumeration) + addKeyList) 

          else:
            self.createEntry()
         
      else:
      
        if dataType.enumeration:
          self.createPulldownMenu(list(dataType.enumeration) + addKeyList) 

        else:
          self.createEntry()
    
    else:
    
      #
      # Links... have to set up a selection list
      #
      
      otherClassName = metaAttrOrRole.otherClass.name
      otherClass = getattr(ChemComp,otherClassName)
      (keyNum,upLinks,downLinks,allKeys) = otherClass.keyPath
      currentKeys = allKeys[-1]
      
      #
      # First get a list of all the relevant objects...
      #
      
      chemCompIndex = list(upLinks).index('chemComp')
      for i in range(chemCompIndex,-1,-1):
        downLinkName = downLinks[len(downLinks)-i-1]
        
        if i == chemCompIndex:
          objects = list(getattr(self.chemComp,downLinkName))
        else:
          parentObjects = objects
          objects = []
          for parentObject in parentObjects:
            objects.extend(getattr(parentObject,downLinkName))
      
      #
      # Then make sure a good key is created for them
      #
      
      objectKeyList = []
      objectDict = {}
      
      for object in objects:
        objectKey = makeObjectKey(object,currentKeys)
       
        objectKeyList.append(objectKey)
        objectDict[objectKey] = object
      
      objectKeyList.sort()
      
      #
      # Only add None if can be 0!
      #
      
      if metaAttrOrRole.locard == 0:
        objectKeyList.append('None')
        objectDict['None'] = None
      
      #
      # Create the relevant widget...
      #
            
      if not hasOrder and (metaAttrOrRole.hicard > 1 or metaAttrOrRole.hicard == -1):
        # Multilist
        self.createScrolledList(objectKeyList, objectDict)
      
      else:
        # Single list
        self.createPulldownMenu(objectKeyList, widgetDict = objectDict) 

    return self.widget
  
  def set(self,value):
    #print "WSet"
    #print value
       
    attrName = self.window.getCurrentAttributeName()

    if self.window.metaRoleOrder[attrName]:
      (attrName,attrIndex) = self.window.metaRoleOrder[attrName]
      hasOrder = True
    else:
      hasOrder = False

    metaAttrOrRole = self.window.metaAttrOrRoleDict[attrName]
    
    dataType = None
    if hasattr(metaAttrOrRole,'dataType'):
      dataType = metaAttrOrRole.dataType
      
    if value == None:
      value = 'None'

    if dataType:
      if dataType.typeCode == 'Boolean':
        value = str(value)
    
      # TODO Can also use dataType.name: gives more specific info (e.g. Word, ...)
      elif dataType.typeCode == 'String':
        pass
        
      elif type(value) in [type([]),type('',)]:
        pass
      
      else:
        value = str(value)
    
    else:
      pass
    
    
    if self.widgetType in ['Entry','PulldownMenu']:
      if self.widgetDict:
        for widgetKey in self.widgetDict.keys():
          if self.widgetDict[widgetKey] == value:
            value = widgetKey

      self.widget.set(value)

    elif self.widgetType in ['ScrolledList']:
      selectedItems = []
      if value:
        if self.widgetDict:
          for widgetKey in self.widgetDict.keys():
            if self.widgetDict[widgetKey] in value:
              selectedItems.append(widgetKey)
        else:
          for tval in value:
            selectedItems.append(tval)
     
      if self.widgetType == 'ScrolledList':
        self.widget.setSelectedItems(selectedItems)
        
    elif self.widgetType in ['ValueList']:
      
      pass
      
    else:
      raise "Unknown widget type %s" % self.widgetType
    
    return 1
  
  def get(self):

    #print "WGet"

    if self.widgetType in ['Entry','PulldownMenu']:
      value = self.widget.get()
    elif self.widgetType == 'ScrolledList':
      value = self.widget.getSelectedItems()
    elif self.widgetType == 'ValueList':
      value = self.widget.getAllItems()
    else:
      raise "Unknown widget type %s" % self.widgetType
    
    attrName = self.window.getCurrentAttributeName()

    if self.window.metaRoleOrder[attrName]:
      (attrName,attrIndex) = self.window.metaRoleOrder[attrName]
      hasOrder = True
    else:
      hasOrder = False

    metaAttrOrRole = self.window.metaAttrOrRoleDict[attrName]
    
    dataType = None
    if hasattr(metaAttrOrRole,'dataType'):
      dataType = metaAttrOrRole.dataType
    
    # TODO MESSY: clean up!
    if value in ['None','']:
      value = None
    
    elif dataType:
      if dataType.typeCode in ['Boolean','Int','Long']:
        value = returnInt(value)

      elif dataType.typeCode in ['Float','Double']:
        value = returnFloat(value)

      elif dataType.typeCode == 'String':
        if type(value) == type([]):
          returnValue = []
          for tval in value:
            if tval in ['None','']:
              continue
            returnValue.append(tval)
          value = returnValue
        
      else:
        pass
    
    else:
      if type(value) == type([]):
      
        returnValue = []
        
        for tval in value:
        
          if tval in ['None','']:
            continue
            
          if self.widgetDict:
            tval = self.widgetDict[tval]
          returnValue.append(tval)
          
        value = returnValue

      else:
        if self.widgetDict:
          value = self.widgetDict[value]
    
    return value
  
  def configure(self, *args, **keywds):
  
    self.widget.configure(*args,**keywds)

class CcpnAttributeEditPopup(BasePopup):
    
  def setupMetaInfo(self,metaClass = None, keyNames = None, getRoles = False, ignoreNames = [], orderRoles = [],setDummy = False):
  
    """
    metaclass:
    ['__doc__', '__init__', '__initialiser__', '__module__', 'addAttribute', 'addConstraint', 'addKeyName', 'addOperation', 'addRole', 'addTaggedValue', 'allowedTags', 'attributeAllNames', 'attributeNames', 'booleanMetaParameters', 'checkValid', 'constraints', 'container', 'documentation', 'elementDict', 'fixedMetaParameters', 'getAllSupertypes', 'getChildClasses', 'getNonabstractSubclasses', 'getParentClass', 'getSubclasses', 'hasSpecialConstructor', 'hasSpecialDestructor', 'isAbstract', 'isLeaf', 'isRoot', 'isSingleton', 'keyNames', 'metaObjFromQualName', 'metaParameterDefaults', 'metaParameterNames', 'name', 'operationAllNames', 'operationNames', 'package', 'parentRole', 'qualifiedName', 'roleAllNames', 'roleNames', 'setSuperclass', 'sortElementNames', 'subclasses', 'superclass', 'taggedValues', 'visibility']
    
    metaattr:
    ['__doc__', '__init__', '__initialiser__', '__module__', 'addConstraint', 'addTaggedValue', 'allowedTags', 'baseName', 'booleanMetaParameters', 'checkValid', 'clazz', 'constraints', 'container', 'dataType', 'defaultValue', 'documentation', 'fixedMetaParameters', 'hicard', 'isAutomatic', 'changeability', 'isDerived', 'isOrdered', 'isUnique', 'locard', 'metaObjFromQualName', 'metaParameterDefaults', 'metaParameterNames', 'name', 'operationName', 'qualifiedName', 'scope', 'taggedValues', 'visibility']    
    
    metadatatype:
    ['__doc__', '__init__', '__initialiser__', '__module__', 'addConstraint', 'addTaggedValue', 'allowedTags', 'booleanMetaParameters', 'checkValid', 'constraints', 'container', 'documentation', 'elementDict', 'enumeration', 'fixedMetaParameters', 'getAllSupertypes', 'isAbstract', 'isLeaf', 'isOpen', 'isRoot', 'isValid', 'javaSimpleTypeCode', 'javaTypeCode', 'jdbcTypeCode', 'length', 'metaObjFromQualName', 'metaParameterDefaults', 'metaParameterNames', 'name', 'package', 'qualifiedName', 'setSupertype', 'subtypes', 'supertype', 'taggedValues', 'typeCode', 'visibility', 'xmlTypeCode']
    """
    
    if metaClass:
      self.metaClass = metaClass
    else:
      self.metaClass = self.ccpnObject._metaclass
      
    self.metaAttrOrRoleDict = {}
    self.metaRoleOrder = {}
    
    attrNames = self.metaClass.attributeAllNames

    if getRoles:
      for roleName in self.metaClass.roleAllNames:
        role = self.metaClass.elementDict[roleName]
        
        if not role.hierarchy in ['child','parent'] and role.otherClass.package.name in ['ChemComp','RefStereochemistry']:
          if roleName not in attrNames:
            attrNames.append(roleName)
    
    self.validAttrNames = []

    if keyNames == None:
      keyNames = self.metaClass.keyNames
    
    for attrName in attrNames:
      
      if attrName in ignoreNames:
        continue

      if attrName in keyNames: #TODO: BUT SET AS INDICATOR! ALSO ON CREATION NEW!
        continue
      
      metaAttr = self.metaClass.elementDict[attrName]
      
      if metaAttr.isDerived:
        continue
      
      self.metaAttrOrRoleDict[attrName] = metaAttr
      self.metaRoleOrder[attrName] = []
      
      if attrName in orderRoles:
        baseName = metaAttr.baseName
        for i in range(0,metaAttr.hicard):
          orderAttrName = "%s_%d" % (baseName,i+1)

          self.validAttrNames.append(orderAttrName)
          self.metaRoleOrder[orderAttrName] = (attrName,i)
          
          if setDummy:
            self.setDummyValue(1,orderAttrName)
      
      else:
        self.validAttrNames.append(attrName)
        if setDummy:
          self.setDummyValue(metaAttr.hicard,attrName)
   
  def setDummyValue(self, hicard, attrName):
  
    if hicard > 1 or hicard == -1:
      defaultValue = []
    else:
      defaultValue = ''

    setattr(self.ccpnObject,attrName,defaultValue)

  def getValue(self, attrName):

    #print "Get %s" % attrName
    
    if not hasattr(self.ccpnObject,'dummy') and self.metaRoleOrder[attrName]:
      (attrName,attrIndex) = self.metaRoleOrder[attrName]
      valueList = getattr(self.ccpnObject,attrName)
      
      if attrIndex < len(valueList):
        value = valueList[attrIndex]
      else:
        value = None

    else:
      value = getattr(self.ccpnObject,attrName)
    
    self.objectValueEdit['Value'].set(value)
    
  def getInfo(self, attrName):
    
    showInfo('Attribute information',self.metaAttrOrRoleDict[attrName].documentation,parent = self)    
    
  def setValue(self, *extra):

    #print "Set"
    
    value = self.objectValueEdit['Value'].get()
    attrName = self.getCurrentAttributeName()
    realAttrName = attrName
    ignoreList = False
    attrIndex = None
    
    if self.metaRoleOrder[attrName]:
      (realAttrName,attrIndex) = self.metaRoleOrder[attrName]
      
    metaAttrOrRole = self.metaAttrOrRoleDict[realAttrName]
    
    if hasattr(self.ccpnObject,'dummy'):
      attrIndex = None
      if self.metaRoleOrder[attrName]:
        ignoreList = True
      
    else:
      attrName = None
    
    # WARNING: THIS DOES UPDATE OF THE TABLE IF VALUE CHANGED!!
    if value == '':
      value = None
      
    setCcpnObjectValue(self,self.ccpnObject,metaAttrOrRole,value,self.objectValueEdit['Value'],self.updateTables,checkParent = False, checkValid = self.attributeTable.checkValid, attrIndex = attrIndex, attrName = attrName, ignoreList = ignoreList)
      
  def getCurrentAttributeName(self):

    if (self.attributeTable):
      attrName = self.attributeTable.currentObject
    else:
      attrName = None

    return attrName

  def setupTable(self,master, parDict, headings):

    numCols = len(headings)

    editWidgets = []
    editGetCallbacks = []
    editSetCallbacks = []

    x = y = None
    
    for attrName in ['objectValueEdit','objectValueEditType']:
      if not hasattr(self,attrName):
        setattr(self,attrName,{})

    for heading in headings:
    
      #
      # These are for attribute type tables
      #
    
      if heading == 'Attribute':
      
        setFunc = None
        getFunc = None
        self.objectValueEdit[heading] = None
      
      elif heading == 'Value':
      
        setFunc = self.setValue
        getFunc = self.getValue
        self.objectValueEdit[heading] = MultiWidget(self,setFunc,parDict)

      elif heading == 'Documentation':
      
        setFunc = None
        getFunc = self.getInfo
        self.objectValueEdit[heading] = None

      #
      # These are for sysName type tables
      #

      elif heading == 'namingSystem':
        setFunc = lambda selInd = x, selVal = y, attrName = heading: self.setObjectValue(selInd,selVal,attrName = attrName)
        getFunc = lambda object = x, attrName = heading: self.getObjectValue(object,attrName)
        self.objectValueEdit[heading] = PulldownMenu(self, entries = parDict['sysNameList'], callback= setFunc, width=25, do_initial_callback = False)
        self.objectValueEditType[heading] = 'PulldownMenu'

      elif heading == 'specificChemCompVars': 

        setFunc = None
        getFunc = lambda object = x, attrName = heading: self.getMultiListObjectValue(object,attrName)
        self.objectValueEdit[heading] = None

      elif heading == 'altSysNames': 

        setFunc = lambda object = x, attrName = heading: self.setObjectValue(object,attrName = attrName)
        getFunc = lambda object = x, attrName = heading: self.getObjectValue(object,attrName)
        self.objectValueEdit[heading] = ScrolledListbox(self,width = 50,height = 5,selectmode = Tkinter.MULTIPLE,addDeleteButtons = True)
        self.objectValueEditType[heading] = 'ValueList'

      else:
        setFunc = lambda object = x, attrName = heading: self.setObjectValue(object,attrName = attrName)
        getFunc = lambda object = x, attrName = heading: self.getObjectValue(object,attrName)
        self.objectValueEdit[heading] = Entry(self, returnCallback = setFunc, width=25)
        self.objectValueEditType[heading] = 'Entry'

      editWidgets.append(self.objectValueEdit[heading])

      editGetCallbacks.append( getFunc )
      editSetCallbacks.append( setFunc )
        
    table = MultiRowScrolledMatrix(master, headingList=headings,
                                   editWidgets=editWidgets,
                                   editGetCallbacks=editGetCallbacks,
                                   editSetCallbacks=editSetCallbacks)
                                  
    if parDict.has_key('checkValid'):
      table.checkValid = parDict['checkValid']
    else:
      table.checkValid = True

    return table

  def setupAttributeTable(self,master,row,column,chemComp = None, checkValid = True):
    
    parDict = {'checkValid': checkValid, 'chemComp': chemComp}
    headings = ('Attribute', 'Value', 'Documentation')
    self.attributeTable = self.setupTable(master,parDict,headings)
    self.attributeTable.grid(row=row, column=0, sticky=Tkinter.NSEW)

  def updateAttributeTable(self):

    textMatrix = []
    
    if hasattr(self,'ccpnObject'):
      for attrName in self.validAttrNames:

        # Get real value
        if not hasattr(self.ccpnObject,'dummy'):
          # If split up...
          if self.metaRoleOrder[attrName]:
            (attrName,attrIndex) = self.metaRoleOrder[attrName]
          else:
            attrIndex = None

          value = getattr(self.ccpnObject,attrName)

          # If split up...
          if attrIndex != None:
            if attrIndex < len(value):
              value = value[attrIndex]

              if value and type(value) != type(''):
                (keyNum,upLinks,downLinks,allKeys) = value.keyPath
                currentKeys = allKeys[-1]
                value = makeObjectKey(value,currentKeys)

            else:
              value = None
              
          elif hasattr(value,'className'):
            (keyNum,upLinks,downLinks,allKeys,currentStore) = value.keyPath
            currentKeys = allKeys[-1]
            value = makeObjectKey(value,currentKeys)
            

        # Dummy object, just get value straight away
        else:
          value = getattr(self.ccpnObject,attrName)

          if value and type(value) not in [type(1),type(1.1),type([]),type(('',)),type('')]:
            (keyNum,upLinks,downLinks,allKeys) = value.keyPath
            currentKeys = allKeys[-1]
            value = makeObjectKey(value,currentKeys)

        if type(value) in [type([]),type(('',))]:
          value = 'Click (%d selected)' % len(value)
        else:
          value = str(value)

        text = [attrName,value,'Info']
        textMatrix.append(text)

    self.attributeTable.update(objectList=self.validAttrNames, textMatrix=textMatrix)

  def setupSysNameTable(self,master,row,column,headings,sysNameList,chemComp = None,checkValid = True):
  
    parDict= {'sysNameList': sysNameList, 'chemComp': chemComp,'checkValid': checkValid}
    
    self.sysNameHeadings = headings
     
    self.sysNameTable = self.setupTable(master,parDict, headings)
    self.sysNameTable.grid(row=row, column=0, sticky=Tkinter.NSEW)
    
    row += 1
    
    texts = [ 'Add sysName', 'Delete sysName' ]
    commands = [ self.addSysName, self.deleteSysName ]
    buttons = ButtonList(master, texts=texts, commands=commands)
    buttons.grid(row=row, column=0, sticky=Tkinter.EW)

  def resetSysNameWidgetList(self,sysNameList,selectedName):

    selectedIndex = sysNameList.index(selectedName)
    self.objectValueEdit['namingSystem'].replace(sysNameList,selected_index = selectedIndex)
    self.updateSysNameTable()

  def deleteSysName(self):
  
    rowObject = self.getCurrentObject()
    
    if rowObject:
      rowObject.delete()
      self.updateSysNameTable()
      
  def addSysName(self):
    
    # HAS TO BE DEFINED ON LOCAL LEVEL!
    pass 

  def getObjectValue(self, rowObject, attrName):
    
    if self.objectValueEditType[attrName] == 'ValueList':
      items = getattr(rowObject,attrName)
      self.objectValueEdit[attrName].setItems(items)
      
    else:    
      self.objectValueEdit[attrName].set(str(getattr(rowObject,attrName)))
  
  def getMultiListObjectValue(self, rowObject, attrName):
    
    value = []

    if rowObject:
      if attrName == 'specificChemCompVars':
        parentAttrName = 'chemCompVars'
      else:
        parentAttrName = attrName
      value = self.multiSelectPopup(getattr(rowObject,attrName),getattr(self.ccpnObject,parentAttrName))
      setattr(rowObject,attrName,value)
      self.updateTables()
        
  def setObjectValue(self, *extra, **kywds):
    
    attrName = kywds['attrName']
    
    if self.objectValueEditType[attrName] == 'ValueList':
      value = self.objectValueEdit[attrName].getAllItems()
    else:
      value = self.objectValueEdit[attrName].get()

    rowObject = self.getCurrentObject()
    
    if rowObject:
      metaAttr = rowObject._metaclass.elementDict[attrName]
    
      setCcpnObjectValue(self,rowObject,metaAttr,value,self.objectValueEdit[attrName],self.updateTables, checkValid = self.sysNameTable.checkValid)
        
  def getCurrentObject(self):

    if (self.sysNameTable):
      object = self.sysNameTable.currentObject
    else:
      object = None

    return object


  def multiSelectPopup(self,ccpnSelectedObjects,ccpnObjects):
  
    selectionList = []
    selectionDict = {}
    selectedItems = []
    
    if ccpnObjects:

      metaClass = ccpnObjects[0]._metaclass
      keyNames = metaClass.keyNames

      for ccpnObject in ccpnObjects:
        objectKeys = []
        for keyName in keyNames:
          objectKeys.append(getattr(ccpnObject,keyName))
        objectKey = string.join(objectKeys,', ')
        selectionDict[objectKey] = ccpnObject
        selectionList.append(objectKey)
        
        if ccpnObject in ccpnSelectedObjects:
          selectedItems.append(objectKey)

      tempPopup = MultiSelectionListPopup(self,selectionList,selectionDict,selectedItems = selectedItems, text = 'Select chemCompVars', modal = True)

      if hasattr(tempPopup,'isSelectedList') and tempPopup.isSelectedList != None:
        ccpnSelectedObjects = []
        if tempPopup.isSelectedList:
          ccpnSelectedObjects = tempPopup.isSelectedList
    
    return ccpnSelectedObjects

class DummyClass:

  def __init__(self):
  
    self.dummy = True
    
class CreateCcpnObject(CcpnAttributeEditPopup):

  def __init__(self, parent, ccpnParent,classType, chemComp = None, fromTable = None,ignoreNames = [], orderRoles = []):
  
    self.parent = parent
    self.ccpnParent = ccpnParent
    self.classType = classType
    self.fromTable = fromTable
    self.ignoreNames = ignoreNames
    self.orderRoles = orderRoles
    
    if chemComp:
      self.chemComp = chemComp
    else:
      self.chemComp = ccpnParent
    
    BasePopup.__init__(self, parent=parent, title='Edit chemComp (%s, %s) window' % (self.chemComp.molType,self.chemComp.ccpCode), modal = True)

  def body(self, master):

    self.master_frame = master

    self.ccpnObject = DummyClass()
    self.setupMetaInfo(metaClass = self.classType._metaclass, keyNames = [], getRoles = True, setDummy = True,ignoreNames = self.ignoreNames, orderRoles = self.orderRoles)
    
    master.grid_columnconfigure(0, weight=1)

    row = 0
    
    #
    # Attributes frame
    #
    
    master.grid_rowconfigure(row, weight=0)
    frame1 = LabelFrame(master, text="Create new %s" % self.classType._metaclass.name) 
    frame1.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_columnconfigure(0, weight=1)
    
    self.setupAttributeTable(frame1,row,0, chemComp = self.chemComp, checkValid = False)

    row = row + 1
    texts = [ 'Create', 'Dismiss' ]
    commands = [ self.createObject, self.close ]
    buttons = ButtonList(self, texts=texts, commands=commands)
    buttons.grid(row=row, column=0, sticky=Tkinter.EW)

    self.updateTables()
  
  def updateTables(self):
  
    self.updateAttributeTable()
    
  def createObject(self):
  
    message = 'Could not create new %s.' % self.classType._metaclass.name
    created = 0
  
    try:
      createDict = {}

      for attrName in self.validAttrNames:
        if attrName[:2] != '__':
          value = getattr(self.ccpnObject,attrName)
          if value in ['None','']:
            value = None
            
          if self.metaRoleOrder[attrName]:
            if value:
            
              (attrName,index) = self.metaRoleOrder[attrName]
              if not createDict.has_key(attrName):
                createDict[attrName] = []
              createDict[attrName].append(value)
            
          else:
            createDict[attrName] = value
      
      #print createDict
      message = 'Creation error for %s.' % self.classType._metaclass.name
      ccpnObject = self.classType(self.ccpnParent,**createDict)
      
      try:
        ccpnObject.checkAllValid()
        created = 1
      except:
        message = 'Validation error for %s.' % self.classType._metaclass.name
        ccpnObject.delete()
        raise
    
    except Implementation.ApiError, e:
      showError(message, e.error_msg, parent = self)
    
    
    if created:
       
      self.parent.updateTables()

      if self.fromTable:
        self.fromTable.selectObject(ccpnObject)
    
      self.close()
    
  def cancel(self):
    # TODO : WHAT IS THE DEFAULT AGAIN?
    pass

class ChemCompEditPopup(CcpnAttributeEditPopup):

  def __init__(self, parent, chemComp, *args, **kw):

    if not chemComp:
      return
    
    self.ccpnObject = chemComp
    self.chemComp = chemComp
    self.chemComp.isModifiable = True
  
    self.sysNameList = []
    self.setSysNameList()
    
    self.testMode = False
    if kw.has_key('test'):
      self.testMode = kw['test']
      del(kw['test'])
    
    if self.testMode:
      # TODO CHANGE THIS
      #refUrl = getUrl(self.ccpnObject.root,'localTest',os.path.join(os.getcwd(),'local'))
      #self.ccpnObject.chemCompHead.contentStorage.url = refUrl
      print
      print " WARNING NOT WORKING!"
      print
    
    self.popups = {}
    
    BasePopup.__init__(self, parent=parent, title='Edit chemComp (%s, %s) window' % (chemComp.molType,chemComp.ccpCode), **kw)

  def body(self, master):

    self.master_frame = master

    self.setupMetaInfo(keyNames = self.chemComp._metaclass.keyNames, getRoles = True)
    
    master.grid_columnconfigure(0, weight=1)
    row = 0
    
    #
    # ChemComp attributes frame
    #
    
    frame1 = LabelFrame(master, text="Edit chemComp attributes")
    frame1.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_columnconfigure(0, weight=1)
    master.grid_rowconfigure(row, weight=1)
    
    self.setupAttributeTable(frame1,row,0,chemComp = self.chemComp)

    row += 1
    
    #
    # ChemCompSysName frame
    #
    
    frame2 = LabelFrame(master, text="Edit chemCompSysNames")
    frame2.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame2.grid_rowconfigure(0, weight=1)
    frame2.grid_columnconfigure(0, weight=1)
    master.grid_rowconfigure(row, weight=1)

    headings = ('namingSystem', 'sysName', 'specificChemCompVars')
    self.setupSysNameTable(frame2,0,0,headings,self.sysNameList,chemComp = self.chemComp)

    row = row + 1
    texts = [ 'Edit ChemAtoms','Edit ChemBonds', 'Edit ChemAngles','Edit ChemTorsions' ]
    commands = [ self.editChemAtoms,  self.editChemBonds, self.editChemAngles, self.editChemTorsions ]
    buttons = ButtonList(self, texts=texts, commands=commands, spacers = 15)
    buttons.grid(row=row, column=0)

    row = row + 1
    texts = [ 'Edit ChemAtomSets', 'Edit ChemCompVars','Edit LinkEnds', 'Edit Stereochemistries' ]
    commands = [ self.editChemAtomSets, self.editChemCompVar , self.editLinkEnd, self.editStereochem ]
    buttons = ButtonList(self, texts=texts, commands=commands, spacers = 15)
    buttons.grid(row=row, column=0)

    row = row + 1
    texts = [ 'Save', 'Check validity' ]
    commands = [ self.save, self.checkValid ]
    buttons = ButtonList(self, texts=texts, commands=commands, spacers = 15)
    buttons.grid(row=row, column=0)

    self.updateTables()    
  
  def setSysNameList(self):
    
    newNames = []
    
    ctsns = []
    for ct in self.chemComp.chemTorsions:
      for ctsn in ct.sysNames:
        ctsns.append(ctsn)
        
    for sn in self.chemComp.chemAtomSysNames + self.chemComp.sysNames + tuple(ctsns):
      if sn.namingSystem not in self.sysNameList:
        newNames.append(sn.namingSystem)
        self.sysNameList.append(sn.namingSystem)
        
    if newNames:
      self.sysNameList.sort()
      
    return newNames

  def updateSysNameTable(self):
  
    self.sysNames = self.ccpnObject.sysNames

    textMatrix = []
    
    for sysName in self.sysNames:
    
      text = []
      
      for attrName in self.sysNameHeadings:
        value = getattr(sysName,attrName)
        
        if attrName == 'specificChemCompVars':
          numObjects = len(getattr(sysName,'chemCompVars'))
          stringValue = "%s (%d objects)" % (attrName,numObjects) # Special case...
          text.append(stringValue)
        else:
          text.append(str(value))
                
      textMatrix.append(text)

    self.sysNameTable.update(objectList=self.sysNames, textMatrix=textMatrix)

  def addSysName(self):
    
    CreateCcpnObject(self,self.chemComp,ChemComp.ChemCompSysName,fromTable = self.sysNameTable)
    
    newNames = self.setSysNameList()
    if newNames:
      self.resetSysNameWidgetList(self.sysNameList,newNames[0])

  def checkValid(self):
    
    isValid = False
    
    try:
      self.ccpnObject.checkAllValid(complete = True)
      isValid = True
      showInfo('All valid','All valid',parent = self)
    except Implementation.ApiError, e:
      showError('Error in validity check', e.error_msg, parent = self)
    
    return isValid

  def editLinkEnd(self):

    if self.popups.has_key('linkEndEdit'):
      self.popups['linkEndEdit'].open()
    
    else:
      self.popups['linkEndEdit'] = LinkEndEditPopup(self)

  def editChemAtoms(self):
    
    if self.popups.has_key('chemAtomEdit'):
      self.popups['chemAtomEdit'].open()
    
    else:
      self.popups['chemAtomEdit'] = ChemAtomEditPopup(self)

  def editChemAtomSets(self):
    
    if self.popups.has_key('chemAtomSetEdit'):
      self.popups['chemAtomSetEdit'].open()
    
    else:
      self.popups['chemAtomSetEdit'] = ChemAtomSetEditPopup(self)
      
  def editChemBonds(self):
    
    if self.popups.has_key('chemBondEdit'):
      self.popups['chemBondEdit'].open()
    
    else:
      self.popups['chemBondEdit'] = ChemBondEditPopup(self)

  def editChemAngles(self):
    
    if self.popups.has_key('chemAngleEdit'):
      self.popups['chemAngleEdit'].open()
    
    else:
      self.popups['chemAngleEdit'] = ChemAngleEditPopup(self)
 
  def editChemTorsions(self):
    
    if self.popups.has_key('chemTorsionEdit'):
      self.popups['chemTorsionEdit'].open()
    
    else:
      self.popups['chemTorsionEdit'] = ChemTorsionEditPopup(self)

  def editStereochem(self):
    
    if self.popups.has_key('stereochemEdit'):
      self.popups['stereochemEdit'].open()
    
    else:
      self.popups['stereochemEdit'] = StereochemistryEditPopup(self)

  def editChemCompVar(self):
    
    if self.popups.has_key('chemCompVarEdit'):
      self.popups['chemCompVarEdit'].open()
    
    else:
      self.popups['chemCompVarEdit'] = ChemCompVarEditPopup(self)

  def updateTables(self, *extra):

    self.updateAttributeTable()
    self.updateSysNameTable()

  def updateChemCompTable(self, *extra):

    pass
    
  def save(self):
  
    if self.checkValid():
    
      doWrite = True
      
      if self.testMode:
        # TODO THIS NEEDS TO CHANGE - different GUID
        origPath = self.chemComp.contentStorage.path
        self.chemComp.contentStorage.path += time.strftime(".%y%m%d-%H:%M")
      else:
       doWrite = showYesNo("Confirm overwrite","Do you really want to overwrite the existing file?")
    
      if doWrite:

        self.chemComp.save()
        
        if self.testMode:
          # TODO THIS NEEDS TO CHANGE - different GUID
          self.chemComp.path = origPath
      

class ChemAtomEditPopup(CcpnAttributeEditPopup):

  def __init__(self, parent, *args, **kw):
    
    self.popups = {}
    
    self.parent = parent
    self.chemComp = self.parent.ccpnObject
    
    self.ignoreNames = ['linkEnds','linkEnd','chemTorsions','chemBonds','chemAngles']
    self.orderRoles = []
    
    self.chemAtomType = None
    self.setChemAtomList()
        
    BasePopup.__init__(self, parent=parent, title='Edit chemAtoms window', **kw)

  def body(self, master):

    self.master_frame = master
    
    master.grid_columnconfigure(0, weight=1)
    
    #
    # TODO: NEED LINKS TO CHEMBONDS HERE? Or leave like is and only have them listed?
    #
    
    row = 0
    
    #
    # Top Frame
    #
    
    frame = Frame(master)
    frame.grid(row=row, column=0, sticky=Tkinter.W)
    
    frameRow = 0
    label = Label(frame, text='ChemAtom type:')
    label.grid(row=frameRow, column=0, sticky=Tkinter.W)

    self.chemAtomTypeWidget = PulldownMenu(frame, entries=self.chemAtomTypes, callback=self.updateChemAtomType)
    self.chemAtomTypeWidget.grid(row=frameRow, column=1, sticky=Tkinter.W)

    frameRow += 1
    label = Label(frame, text='Currently editing chemAtom:')
    label.grid(row=frameRow, column=0, sticky=Tkinter.W)

    self.chemAtomListWidget = PulldownMenu(frame, entries=self.chemAtomList[self.chemAtomType], callback=self.updateChemAtomList)
    self.chemAtomListWidget.grid(row=frameRow, column=1, sticky=Tkinter.W)
    
    frameRow += 1
    
    label = Label(frame, text='ChemAtomSet:')
    label.grid(row=frameRow, column=0, sticky=Tkinter.W)

    self.chemAtomSetButton = Tkinter.Button(frame, text = '', command = self.selectChemAtomSet)
    self.chemAtomSetButton.grid(row=frameRow, column=1, sticky=Tkinter.W)   

    master.grid_rowconfigure(row, weight=1)

    row += 1
    
    #
    # ChemAtom attributes frame
    #
    
    master.grid_rowconfigure(row, weight=0)
    frame1 = LabelFrame(master, text="Edit chemAtom attributes")
    frame1.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_columnconfigure(0, weight=1)

    self.setupAttributeTable(frame1,0,0,chemComp = self.chemComp)
    master.grid_rowconfigure(row, weight=1)

    row += 1
    
    #
    # ChemAtomSysName frame
    #
    
    master.grid_rowconfigure(row, weight=0)
    frame2 = LabelFrame(master, text="Edit chemAtomSysNames")
    frame2.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame2.grid_rowconfigure(0, weight=1)
    frame2.grid_columnconfigure(0, weight=1)

    headings = ('namingSystem', 'sysName', 'altSysNames', 'specificChemCompVars')
    self.setupSysNameTable(frame2,0,0,headings,self.parent.sysNameList, chemComp = self.chemComp)
    master.grid_rowconfigure(row, weight=1)
    

    row = row + 1
    texts = [ 'Add atom', 'Delete current' ]
    commands = [ self.addChemAtom, self.deleteChemAtom ]
    buttons = ButtonList(self, texts=texts, commands=commands)
    buttons.grid(row=row, column=0, sticky=Tkinter.EW)

    self.updateChemAtomSet()
    self.updateTables()
  
  def updateTables(self, *extra):

    self.updateAttributeTable()
    self.updateSysNameTable()
    
  def setChemAtomList(self):

    self.chemAtomList = {}
    self.chemAtomDict = {}
    self.chemAtomTypes = ['ChemAtom','LinkAtom']
    
    for chemAtomType in self.chemAtomTypes[:]:
      
      self.chemAtomList[chemAtomType] = []
      self.chemAtomDict[chemAtomType] = {}
      
      for chemAtom in self.parent.ccpnObject.findAllChemAtoms(className = chemAtomType):
        chemAtomKey = '%s' % (chemAtom.name)
        self.chemAtomList[chemAtomType].append(chemAtomKey)
        self.chemAtomDict[chemAtomType][chemAtomKey] = chemAtom

      self.chemAtomList[chemAtomType].sort()

  def updateChemAtomType(self,index,chemAtomType):
    
    if not self.chemAtomList[chemAtomType]:
      chemAtomKey = None
      
      for ccChild in self.parent.metaClass.getChildClasses():
        if ccChild.name == 'AbstractChemAtom':
          for acaSub in ccChild.getSubClasses():
            if acaSub.name == chemAtomType:
              metaClass = acaSub
      
    else:
    
      chemAtomKey = self.chemAtomList[chemAtomType][0]
      self.ccpnObject = self.chemAtomDict[chemAtomType][chemAtomKey]

      metaClass = None

    self.setupMetaInfo(metaClass = metaClass,getRoles = True,ignoreNames = self.ignoreNames, orderRoles = self.orderRoles) 
    
    oldChemAtomType = self.chemAtomType
    self.chemAtomType = chemAtomType

    if hasattr(self,'chemAtomListWidget'):
    
      if oldChemAtomType != chemAtomType:
        self.resetChemAtomList()
    
      self.updateChemAtomList(0,chemAtomKey)
      
  def updateChemAtomSet(self):
  
    text = 'None'
    
    if self.ccpnObject and hasattr(self.ccpnObject,'chemAtomSet') and self.ccpnObject.chemAtomSet:
      chemAtomNames = []
      for chemAtom in self.ccpnObject.chemAtomSet.chemAtoms:
        chemAtomNames.append(chemAtom.name)
      text = "%s (atoms '%s')" % (self.ccpnObject.chemAtomSet.name,string.join(chemAtomNames,','))
        
    self.chemAtomSetButton.config(text = text)

  def selectChemAtomSet(self):
  
    if self.popups.has_key('chemAtomSetEdit'):
      self.popups['chemAtomSetEdit'].open()
    
    else:
      # TODO: make specific!!
      self.popups['chemAtomSetEdit'] = ChemAtomSetEditPopup(self.parent)

  def updateSysNameTable(self):
  
    self.sysNames = self.parent.ccpnObject.findAllChemAtomSysNames(atomName = self.ccpnObject.name)

    textMatrix = []
    
    for chemAtomSysName in self.sysNames:
    
      text = []
      
      for attrName in self.sysNameHeadings:
        value = getattr(chemAtomSysName,attrName)
        
        if attrName == 'specificChemCompVars':
          numObjects = len(getattr(chemAtomSysName,'chemCompVars'))
          stringValue = "%s (%d objects)" % (attrName,numObjects) # Special case...
          text.append(stringValue)
        else:
          text.append(str(value))
      
      textMatrix.append(text)

    self.sysNameTable.update(objectList=self.sysNames, textMatrix=textMatrix)

  def addSysName(self):
    
    CreateCcpnObject(self,self.chemComp,ChemComp.ChemAtomSysName,fromTable = self.sysNameTable)

    newNames = self.parent.setSysNameList()
    if newNames:
      self.resetSysNameWidgetList(self.parent.sysNameList,newNames[0])
      self.parent.resetSysNameWidgetList(self.parent.sysNameList,newNames[0])
    
  def resetChemAtomList(self):
  
    self.chemAtomListWidget.replace(self.chemAtomList[self.chemAtomType])
    
  def updateChemAtomList(self,index,chemAtomKey):

    if chemAtomKey:

      chemAtom = self.chemAtomDict[self.chemAtomType][chemAtomKey]
      self.ccpnObject = chemAtom
      
    else:
      
      self.ccpnObject = None
    
    if hasattr(self,'chemAtomSetButton'):
      self.updateChemAtomSet()
    
    if hasattr(self,'attributeTable'):
      self.updateTables()

  def addChemAtom(self):
    
    currentClass = getattr(ChemComp,self.chemAtomType)
    CreateCcpnObject(self,self.chemComp,currentClass,ignoreNames = self.ignoreNames, orderRoles = self.orderRoles)
    self.setChemAtomList()
    self.resetChemAtomList()
      
  def deleteChemAtom(self):
    
    chemAtomKey = self.chemAtomListWidget.getSelected()
    chemAtom = self.chemAtomDict[self.chemAtomType][chemAtomKey]

    doDelete = showYesNo("Confirm atom delete","Do you really want to delete atom %s?" % chemAtomKey)
    
    if doDelete:
      chemAtom.delete()
      self.setChemAtomList()
      self.resetChemAtomList()

class ChemAtomSetEditPopup(CcpnAttributeEditPopup):

  def __init__(self, parent, *args, **kw):
    
    self.popups = {}
    
    self.parent = parent
    self.chemComp = self.parent.ccpnObject
    
    self.ignoreNames = []
    self.orderRoles = []
    
    self.setChemAtomSetList()
        
    BasePopup.__init__(self, parent=parent, title='Edit chemAtomSets window', **kw)

  def body(self, master):

    self.master_frame = master
    
    master.grid_columnconfigure(0, weight=1)
    
    row = 0
    
    #
    # Top Frame
    #
    
    frame = Frame(master)
    frame.grid(row=row, column=0, sticky=Tkinter.W)
    
    frameRow = 0
    label = Label(frame, text='Currently editing chemAtomSet:')
    label.grid(row=frameRow, column=0, sticky=Tkinter.W)

    self.chemAtomSetListWidget = PulldownMenu(frame, entries=self.chemAtomSetList, callback=self.updateChemAtomSetList)
    self.chemAtomSetListWidget.grid(row=frameRow, column=1, sticky=Tkinter.W)
    
    master.grid_rowconfigure(row, weight=1)

    self.setupMetaInfo(metaClass = self.ccpnObject._metaclass,getRoles = True,ignoreNames = self.ignoreNames, orderRoles = self.orderRoles) 

    row += 1
    
    #
    # ChemAtomSet attributes frame
    #
    
    master.grid_rowconfigure(row, weight=0)
    frame1 = LabelFrame(master, text="Edit chemAtomSet attributes")
    frame1.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_columnconfigure(0, weight=1)

    self.setupAttributeTable(frame1,0,0,chemComp = self.chemComp)
    master.grid_rowconfigure(row, weight=1)

    row += 1
    
    #
    # ChemAtomSetSysName frame
    #
    
    master.grid_rowconfigure(row, weight=0)
    frame2 = LabelFrame(master, text="Edit chemAtomSetSysNames")
    frame2.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame2.grid_rowconfigure(0, weight=1)
    frame2.grid_columnconfigure(0, weight=1)

    headings = ('namingSystem', 'sysName', 'altSysNames', 'specificChemCompVars')
    self.setupSysNameTable(frame2,0,0,headings,self.parent.sysNameList, chemComp = self.chemComp)
    master.grid_rowconfigure(row, weight=1)
    

    row = row + 1
    texts = [ 'Add atomSet', 'Delete current', 'Edit chemAtomSetVars' ]
    commands = [ self.addChemAtomSet, self.deleteChemAtomSet, self.editChemAtomSetVars ]
    buttons = ButtonList(self, texts=texts, commands=commands)
    buttons.grid(row=row, column=0, sticky=Tkinter.EW)

    self.updateTables()
  
  def updateTables(self, *extra):

    self.updateAttributeTable()
    self.updateSysNameTable()
    
  def setChemAtomSetList(self):

    self.chemAtomSetList = []
    self.chemAtomSetDict = {}

    for chemAtomSet in self.parent.ccpnObject.chemAtomSets:
      chemAtomSetKey = '%s' % (chemAtomSet.name)
      self.chemAtomSetList.append(chemAtomSetKey)
      self.chemAtomSetDict[chemAtomSetKey] = chemAtomSet

    self.chemAtomSetList.sort()

  def updateSysNameTable(self):
  
    self.sysNames = self.parent.ccpnObject.findAllChemAtomSysNames(atomName = self.ccpnObject.name)

    textMatrix = []
    
    for chemAtomSetSysName in self.sysNames:
    
      text = []
      
      for attrName in self.sysNameHeadings:
        value = getattr(chemAtomSetSysName,attrName)
        
        if attrName == 'specificChemCompVars':
          numObjects = len(getattr(chemAtomSetSysName,'chemCompVars'))
          stringValue = "%s (%d objects)" % (attrName,numObjects) # Special case...
          text.append(stringValue)
        else:
          text.append(str(value))
      
      textMatrix.append(text)

    self.sysNameTable.update(objectList=self.sysNames, textMatrix=textMatrix)

  def addSysName(self):
    
    CreateCcpnObject(self,self.chemComp,ChemComp.ChemAtomSetSysName,fromTable = self.sysNameTable)

    newNames = self.parent.setSysNameList()
    if newNames:
      self.resetSysNameWidgetList(self.parent.sysNameList,newNames[0])
      self.parent.resetSysNameWidgetList(self.parent.sysNameList,newNames[0])
    
  def resetChemAtomSetList(self):
  
    self.chemAtomSetListWidget.replace(self.chemAtomSetList)
    
  def updateChemAtomSetList(self,index,chemAtomSetKey):

    if chemAtomSetKey:

      chemAtomSet = self.chemAtomSetDict[chemAtomSetKey]
      self.ccpnObject = chemAtomSet
      
    else:
      
      self.ccpnObject = None
    
    if hasattr(self,'attributeTable'):
      self.updateTables()

  def addChemAtomSet(self):
    
    CreateCcpnObject(self,self.chemComp,ChemComp.ChemAtomSet,ignoreNames = self.ignoreNames, orderRoles = self.orderRoles)
    self.setChemAtomSetList()
    self.resetChemAtomSetList()
      
  def deleteChemAtomSet(self):
    
    chemAtomSetKey = self.chemAtomSetListWidget.getSelected()
    chemAtomSet = self.chemAtomDict[chemAtomSetKey]

    doDelete = showYesNo("Confirm atom set delete","Do you really want to delete atomSet %s?" % chemAtomSetKey)
    
    if doDelete:
      chemAtomSet.delete()
      self.setChemAtomSetList()
      self.resetChemAtomSetList()

  def editChemAtomSetVars(self):
  
    if self.popups.has_key('chemAtomSetVarEdit') and self.popups['chemAtomSetVarEdit'].has_key(self.ccpnObject):
      self.popups['chemAtomSetVarEdit'][self.ccpnObject].open()
    
    else:
      if not self.popups.has_key('chemAtomSetVarEdit'):
        self.popups['chemAtomSetVarEdit'] = {}
      self.popups['chemAtomSetVarEdit'][self.ccpnObject] = ChemAtomSetVarEditPopup(self,chemComp = self.parent.ccpnObject)

class GenericEditPopup(CcpnAttributeEditPopup):

  def __init__(self, parent, *args, **kw):
    
    self.popups = {}
    
    self.parent = parent
    
    if kw.has_key('chemComp'):
      self.chemComp = kw['chemComp']
      del(kw['chemComp'])
    else:
      self.chemComp = self.parent.ccpnObject
    
    self.setSpecificInit()
    
    self.setGenericList()
      
    for ccChild in self.parent.metaClass.getChildClasses():
      if ccChild.name == self.className:
        metaClass = ccChild
        
    self.setupMetaInfo(metaClass = metaClass,getRoles = True,ignoreNames = self.ignoreNames, orderRoles = self.orderRoles)

    BasePopup.__init__(self, parent=parent, title=self.titleText, **kw)
    
  def body(self, master):

    self.master_frame = master
    
    master.grid_columnconfigure(0, weight=1)
        
    row = 0
    
    #
    # Top Frame
    #
    
    frame = Frame(master)
    frame.grid(row=row, column=0, sticky=Tkinter.W)
    
    frameRow = 0

    label = Label(frame, text='Currently editing %s:' % self.className)
    label.grid(row=frameRow, column=0, sticky=Tkinter.W)

    self.listWidget = PulldownMenu(frame, entries=self.genericList, callback=self.updateGenericList)
    self.listWidget.grid(row=frameRow, column=1, sticky=Tkinter.W)
    master.grid_rowconfigure(row, weight=1)
    
    row += 1
    
    #
    # Attributes frame
    #
    
    master.grid_rowconfigure(row, weight=0)
    frame1 = LabelFrame(master, text="Edit %s attributes" % self.className)
    frame1.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_columnconfigure(0, weight=1)

    self.setupAttributeTable(frame1,0,0,chemComp = self.chemComp)
    master.grid_rowconfigure(row, weight=1)
    
    #
    # Possible sysName frame...
    #
    
    row = self.setupSysNameFrame(row,master)

    row += 1
        
    texts = [ 'Add %s' % self.className, 'Delete current' ]
    commands = [ self.addGeneric, self.deleteGeneric ]
    buttons = ButtonList(self, texts=texts, commands=commands)
    buttons.grid(row=row, column=0, sticky=Tkinter.EW)

    self.updateTables() 

  def setupSysNameFrame(self,row,master):
  
    return row
  
  def setSpecificInit(self):
    pass

  def updateTables(self, *extra):

    self.updateAttributeTable()
    self.updateSysNameTable()
  
  def updateSysNameTable(self):
  
    pass
    
  def setGenericList(self):

    self.genericList = []
    self.genericDict = {}
    
    objects = self.getSpecificObjects()
    for object in objects:
      genericKey = self.getSpecificKey(object)
      self.genericList.append(genericKey)
      self.genericDict[genericKey] = object

    self.genericList.sort()
  
  def updateGenericList(self,index,value):

    listKey = self.genericList[index]
    self.ccpnObject = self.genericDict[listKey]
    
    if hasattr(self,'attributeTable'):
      self.updateTables()

  def resetGenericList(self):
  
    self.listWidget.replace(self.genericList)

  def addGeneric(self):
    
    CreateCcpnObject(self,self.chemComp,getattr(ChemComp,self.className),ignoreNames = self.ignoreNames, orderRoles = self.orderRoles)
    self.setGenericList()
    self.resetGenericList()
      
  def deleteGeneric(self):
    
    listKey = self.listWidget.getSelected()
    object = self.genericDict[listKey]

    doDelete = showYesNo("Confirm delete","Do you really want to delete %s %s?" % (self.className,listKey))
    
    if doDelete:
      object.delete()
      self.setGenericList()
      self.resetGenericList()

#
# ChemCompVar Editor
#

class ChemCompVarEditPopup(GenericEditPopup):

  #
  # ChemCompVarEdit - allow (de)selection based on *separate* descriptor list (e.g. prot:HD2, deprot:H3, ...!!)
  #
  # Also allow ADDING of new descriptor to ONLY selected vars...
  #

  def setSpecificInit(self):
    
    self.ignoreNames = []
    self.orderRoles = []
    self.className = 'ChemCompVar'
    self.titleText = 'Edit chemCompVars window'
  
  def getSpecificObjects(self):
  
    return self.parent.ccpnObject.chemCompVars
  
  def getSpecificKey(self,ccv):
  
    return '%s:%s' % (ccv.linking,ccv.descriptor)

  def addGeneric(self):
    
    addMode = showYesNo("Creation type","Do you want to add a descriptor to existing chemCompVars?")
    
    if addMode:
      ChemCompVarAddPopup(self,self.chemComp)
    else:
      CreateCcpnObject(self,self.chemComp,getattr(ChemComp,self.className),ignoreNames = self.ignoreNames, orderRoles = self.orderRoles)

    self.setGenericList()
    self.resetGenericList()

#
# Add ChemCompVar window
#

class ChemCompVarAddPopup(ChemCompVarEditPopup):

  def __init__(self, parent, chemComp, *args, **kw):
    
    self.popups = {}
    
    self.parent = parent
    self.chemComp = chemComp
    
    self.setSpecificInit()
    
    self.setupChemAtomList()
    self.setGenericList()
      
    BasePopup.__init__(self, parent=parent, title='Add a chemCompVar', **kw)
    
  def body(self, master):

    self.master_frame = master
    
    master.grid_columnconfigure(0, weight=1)
        
    row = 0
    
    #
    # Top Frame
    #
    
    frame = Frame(master)
    frame.grid(row=row, column=0, sticky=Tkinter.W)
    
    frameRow = 0
    
    #
    # Define the new descriptor
    #

    label = Label(frame, text='New descriptor:')
    label.grid(row=frameRow, column=0, sticky=Tkinter.W)
    
    self.descTypeWidget = PulldownMenu(frame, entries=['prot','deprot','link'], callback=self.updateDescType)
    self.descTypeWidget.grid(row=frameRow, column=1, sticky=Tkinter.W)

    self.descAtomWidget = PulldownMenu(frame, entries=self.chemAtomList, callback=self.updateDescAtom)
    self.descAtomWidget.grid(row=frameRow, column=2, sticky=Tkinter.W)
    
    frameRow += 1
    
    #
    # Multiselection for chemCompVars this descriptor should be copied to...
    #
    
    self.chemCompVarsWidget = ScrolledListbox(frame,width = 50,height = 5,
                                              selectmode = Tkinter.MULTIPLE,
                                              initial_list = self.genericList)

    self.chemCompVarsWidget.grid(row=frameRow, column=0, columnspan = 3, sticky=Tkinter.W)
    
    #
    # Selector by descriptor...
    #
    
    row += 1
    
    texts = [ 'Replace selection by linking/descriptor']
              
    commands = [ self.replaceChemCompVars]
    buttons = ButtonList(self, texts=texts, commands=commands)
    buttons.grid(row=row, column=0, columnspan = 3, sticky=Tkinter.EW)

    row += 1
    
    texts = [ 'Add to selection by linking/descriptor']
              
    commands = [self.addSelectedChemCompVars]
    buttons = ButtonList(self, texts=texts, commands=commands)
    buttons.grid(row=row, column=0, columnspan = 3, sticky=Tkinter.EW)

    row += 1
    
    texts = [ 'Delete from selection by linking/descriptor']
              
    commands = [ self.deselectChemCompVars]
    buttons = ButtonList(self, texts=texts, commands=commands)
    buttons.grid(row=row, column=0, columnspan = 3, sticky=Tkinter.EW)
    
    #
    # TODO: THEN HAVE CREATION WHERE ATOMS, ... CAN BE SELECTED - ALL THESE WILL BE 
    # COPIED USING THE EXISTING DESCRIPTORS!! NEW VARS HAVE TO BE ADDED FOR EACH!!!! Not trivial.
    #
    
    #
    # Based on the existing ccv, have to find *common* sets, create using those.
    # Probably easiest to FIRST create the separate var, then use that as comparison to create
    # the new ones... 
    #
    
    row += 1
        
    texts = [ 'Create ChemCompVar' ]
    commands = [ self.createChemCompVar]
    buttons = ButtonList(self, texts=texts, commands=commands)
    buttons.grid(row=row, column=0, sticky=Tkinter.EW)

  def setupChemAtomList(self,elementSymbol = None,excludeElementSymbol = None):
  
    self.chemAtomList = []
    self.chemAtomDict = {}
    
    keywds = {}
    if elementSymbol:
      keywds['elementSymbol'] = elementSymbol
    
    objects = list(self.chemComp.findAllChemAtoms(className = 'ChemAtom',**keywds))
    
    if excludeElementSymbol:
      for object in objects[:]:
        if object.elementSymbol == excludeElementSymbol:
          objects.pop(objects.index(object))
    
    for object in objects:
      genericKey = "%s" % object.name
      self.chemAtomList.append(genericKey)
      self.chemAtomDict[genericKey] = object

    self.chemAtomList.sort()

  def updateDescType(self,index,value):
  
    if value in ['deprot','prot']:
      self.setupChemAtomList(elementSymbol = 'H')
    else:
      self.setupChemAtomList(excludeElementSymbol = 'H')
      
    self.resetDescAtom()
    
  def updateDescAtom(self,index,value):
  
    pass
    
  def resetDescAtom(self):
  
    if hasattr(self,'descAtomWidget'):
      self.descAtomWidget.replace(self.chemAtomList)
 
  def getSpecificObjects(self):
  
    return self.parent.parent.ccpnObject.chemCompVars
  
  def replaceChemCompVars(self):
  
    selectedCcvs = self.selectChemCompVars()
    self.chemCompVarsWidget.setSelectedItems(selectedCcvs)

  def addSelectedChemCompVars(self):
  
    selectedCcvs = self.selectChemCompVars()
    entries = list(self.chemCompVarsWidget.getItems())
    
    for selectedCcv in selectedCcvs:
      ind = entries.index(selectedCcv)
      self.chemCompVarsWidget.select_set(ind)

  def deselectChemCompVars(self):

    selectedCcvs = self.selectChemCompVars()
    entries = list(self.chemCompVarsWidget.getItems())
    
    for selectedCcv in selectedCcvs:
      ind = entries.index(selectedCcv)
      self.chemCompVarsWidget.select_clear(ind)
  
  def selectChemCompVars(self):
    
    linkings = []
    descriptors = []
    
    for ccv in self.getSpecificObjects():
      if ccv.linking not in linkings:
        linkings.append(ccv.linking)
      
      descriptorDict = getDescriptorDict(ccv.descriptor)
      
      for descKey in descriptorDict.keys():
        for atomName in descriptorDict[descKey]:
          descriptor = "%s:%s" % (descKey,atomName)
          if descriptor not in descriptors:
            descriptors.append(descriptor)
    
    linkings.sort()
    descriptors.sort()
    
    selectionList = []
    selectionDict = {}
    for link in linkings:
      selectionKey = "Linking: %s" % link
      selectionList.append(selectionKey)
      selectionDict[selectionKey] = ('linking',link)
    for descriptor in descriptors:
      selectionKey ="Descriptor: %s" % descriptor
      selectionList.append(selectionKey)
      selectionDict[selectionKey] = ('descriptor',descriptor)
      
    tempPopup = MultiSelectionListPopup(self,selectionList,selectionDict, text = 'Select linking and descriptors', modal = True)

    if hasattr(tempPopup,'isSelectedList') and tempPopup.isSelectedList != None:
      selectedCodes = []
      
      if tempPopup.isSelectedList:

        linkings = []
        descriptors = []
        for (codeType,codeName) in tempPopup.isSelectedList:
          if codeType == 'linking':
            linkings.append(codeName)
          else:
            descDict = getDescriptorDict(codeName)
            descKey = descDict.keys()[0]
            descriptors.append((descKey,descDict[descKey][0]))
        
        #
        # Now select the chemCompVars...
        #
        
        selectedCcvs = []

        for ccv in self.getSpecificObjects():
          if not linkings or ccv.linking in linkings:
            addCcv = None
            if not descriptors:
              addCcv = self.getSpecificKey(ccv)
            else:
              ccvDescriptorDict = getDescriptorDict(ccv.descriptor)
              ccvMatch = True
              for (descKey,descAtom) in descriptors:
                if ccvDescriptorDict.has_key(descKey):
                  if descAtom not in ccvDescriptorDict[descKey]:
                    ccvMatch = False
                else:
                  ccvMatch = False
              
              if ccvMatch:     
                addCcv = self.getSpecificKey(ccv)
  
            if addCcv and addCcv not in selectedCcvs:
              selectedCcvs.append(addCcv)
    
    return selectedCcvs
    
  def createChemCompVar(self):

    #
    # TODO: CREATION IS CURRENTLY TOO SIMPLE - HAVE TO USE INFORMATION FROM WINDOW!!!
    #
    
    #
    # So... do the following:
    #
    # 1. Based on the new linking and/or descriptor, make a list of new valid linking/descriptor
    # combinations (based on the selection). This should check for atom names and the likes.
    # 
    # 2. Pop up a window where the attributes valid for ALL the new chemCompVars can be set (??).
    #    This includes a list of chemAtoms? But has to be modified depending on the other
    #    parts of the descriptor/linking... not trivial
    #
    # Might be best to have a 'skeleton' chemComp attribute window FIRST, then say which ones
    # you want to 'add' this skeleton to... can then exclude possibilities from the list that
    # don't make sense based on the skeleton...
    #
    #

    CreateCcpnObject(self,self.chemComp,getattr(ChemComp,self.className),ignoreNames = self.ignoreNames, orderRoles = self.orderRoles)
    self.setGenericList()
    self.resetGenericList()

class LinkEndEditPopup(GenericEditPopup):

  def setSpecificInit(self):
    
    self.ignoreNames = []
    self.orderRoles = []
    self.className = 'LinkEnd'
    self.titleText = 'Edit linkEnds window'
  
  def getSpecificObjects(self):
  
    return self.parent.ccpnObject.linkEnds
  
  def getSpecificKey(self,linkEnd):
  
    return '%s' % (linkEnd.linkCode)

class ChemBondEditPopup(GenericEditPopup):

  def setSpecificInit(self):
    
    self.ignoreNames = []
    self.orderRoles = ['chemAtoms']
    self.className = 'ChemBond'
    self.titleText = 'Edit chemBonds window'
  
  def getSpecificObjects(self):
  
    return self.parent.ccpnObject.chemBonds
  
  def getSpecificKey(self,chemBond):
  
    return '%s-%s' % (chemBond.chemAtoms[0].name,chemBond.chemAtoms[1].name)

class ChemAtomSetVarEditPopup(GenericEditPopup):

  def setSpecificInit(self):
    
    self.ignoreNames = []
    self.orderRoles = ['chemAtoms']
    self.className = 'ChemAtomSetVar'
    self.titleText = 'Edit chemAtomSetVars window'
    
    self.parentChemAtomSet = self.parent.ccpnObject
      
  def getSpecificObjects(self):
  
    return self.parent.ccpnObject.chemAtomSetVars
  
  def getSpecificKey(self,chemAtomSetVar):
  
    chemAtomKey = []
    
    for ca in chemAtomSetVar.chemAtoms:
      chemAtomKey.append(ca.name)
    
    if chemAtomSetVar.distCorr == None:
      distCorrText = "(None)"
    else:
      distCorrText = " (%7.1f)" % chemAtomSetVar.distCorr

    return '%d:%s%s' % (chemAtomSetVar.serial,str(chemAtomKey),distCorrText)

  def addGeneric(self):
    
    CreateCcpnObject(self,self.parentChemAtomSet,getattr(ChemComp,self.className),ignoreNames = self.ignoreNames, orderRoles = self.orderRoles)
    self.setGenericList()
    self.resetGenericList()

class ChemAngleEditPopup(GenericEditPopup):

  def setSpecificInit(self):
    
    self.ignoreNames = []
    self.orderRoles = ['chemAtoms']
    self.className = 'ChemAngle'
    self.titleText = 'Edit chemAngles window'
  
  def getSpecificObjects(self):
  
    return self.parent.ccpnObject.chemAngles
  
  def getSpecificKey(self,chemAngle):
  
    return '%s-%s-%s' % (chemAngle.chemAtoms[0].name,chemAngle.chemAtoms[1].name,chemAngle.chemAtoms[2].name)


class ChemTorsionEditPopup(GenericEditPopup):

  def setSpecificInit(self):
    
    self.ignoreNames = []
    self.orderRoles = ['chemAtoms']
    self.className = 'ChemTorsion'
    self.titleText = 'Edit chemTorsions window'

    self.sysNameList = []
    self.parent.setSysNameList()
  
  def getSpecificObjects(self):
  
    return self.parent.ccpnObject.chemTorsions
  
  def getSpecificKey(self,chemTorsion):
  
    return '%s:%s-%s-%s-%s' % (chemTorsion.name,chemTorsion.chemAtoms[0].name,
                                                chemTorsion.chemAtoms[1].name,
                                                chemTorsion.chemAtoms[2].name,
                                                chemTorsion.chemAtoms[3].name)

  def setupSysNameFrame(self,row,master):
    
    row +=1
    
    #
    # ChemTorsionSysName frame
    #
    
    frame2 = LabelFrame(master, text="Edit chemTorsionSysNames")
    frame2.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame2.grid_rowconfigure(0, weight=1)
    frame2.grid_columnconfigure(0, weight=1)
    master.grid_rowconfigure(row, weight=1)

    headings = ('namingSystem', 'sysName')
    self.setupSysNameTable(frame2,0,0,headings,self.sysNameList,chemComp = self.chemComp)
    
    return row
    
  def addSysName(self):
    
    CreateCcpnObject(self,self.ccpnObject,ChemComp.ChemTorsionSysName, chemComp = self.chemComp,fromTable = self.sysNameTable)
    
    newNames = self.parent.setSysNameList()
    if newNames:
      self.resetSysNameWidgetList(self.sysNameList,newNames[0])

  def updateSysNameTable(self):
  
    self.sysNames = self.ccpnObject.sysNames

    textMatrix = []
    
    for chemTorsionSysName in self.sysNames:
    
      text = []
      
      for attrName in self.sysNameHeadings:
        text.append(str(getattr(chemTorsionSysName,attrName)))
      
      textMatrix.append(text)

    self.sysNameTable.update(objectList=self.sysNames, textMatrix=textMatrix)

class StereochemistryEditPopup(GenericEditPopup):
  
  #
  # TODO: still problem with setting refStereochemistry... fix when there's data for it.
  #

  def setSpecificInit(self):
    
    self.ignoreNames = []
    self.orderRoles = []
    self.className = 'Stereochemistry'
    self.titleText = 'Edit stereochemistries window'
  
  def getSpecificObjects(self):
  
    return self.parent.ccpnObject.stereochemistries
  
  def getSpecificKey(self,stereochem):
  
    return '%d:%s,%s' % (stereochem.serial,stereochem.value)

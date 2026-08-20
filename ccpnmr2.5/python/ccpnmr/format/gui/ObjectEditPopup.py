
"""
======================COPYRIGHT/LICENSE START==========================

ObjectEditPopup.py: generic GUI setup to modify CCPN objects

Copyright (C) 2005-2008 Wim Vranken (European Bioinformatics Institute)

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
import string

import Tkinter

from memops.universal.Io import joinPath
from memops.universal.Util import returnInt, returnFloat

from ccpnmr.format.general.Io import getHelpUrlDir

from memops.gui.Entry import Entry
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError, showInfo, showYesNo
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.ToggleLabel import ToggleLabel
from memops.gui.ScrolledListbox import ScrolledListbox
from memops.gui.Separator import Separator

from ccp.general.Util import setCurrentStore

from ccpnmr.format.general.Util import createSelection

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup
from ccpnmr.format.gui.SelectionListPopup import SelectionListPopup

# TODO TODO THIS SHOULD PROBABLY BE SPLIT OUT! Also could use widget creation stuff from here?
from ccpnmr.format.gui.ImportExportFormatPopup import getCcpnObjects

from memops.api import Implementation
from ccp.api.nmr import NmrEntry
from ccp.api.general import Citation
from ccp.api.general import Affiliation

class ObjectEditPopup(TemporaryBasePopup):
  
  # Need to set these!!
  help_url = None
  objectName = ''
  storeClassNameList = []
  
  def __init__(self, parent, project, selectedObject = None, createCallback = None, ccpnParent = None, ccpnParentSelect = None):

    self.project = project
    self.parent = parent
    self.selectedObject = selectedObject
    if selectedObject:
      self.metaClass = selectedObject._metaclass
    
    self.createCallback = createCallback
    
    self.setCurrentStores()

    if ccpnParentSelect:
      
      ccpnParent = None
      
      (objectName,linkName) = ccpnParentSelect
      parentObjects = getCcpnObjects(self.project,linkName)
      (parentList,parentDict) = createSelection(parentObjects)
 
      if parentList:
        interaction = SelectionListPopup(self.parent, parentList, title = 'Select associated %s' % objectName, text = 'Select associated %s' % objectName, dismissText = 'Create new', selectionDict = parentDict, dismissButton = True, modal = True)
          
        if interaction.isSelected:
          ccpnParent = interaction.selection
      
      if not ccpnParent:
        # Make a new one, assuming that all within one package (if dependencies), so use the affiliation object
        storeLinkName = '.'.join(linkName.split('.')[:2])
        storeObject = getCcpnObjects(self.project,storeLinkName)
        if storeObject:
          storeObject = storeObject[0]
          interaction = globals()['%sCreatePopup' % objectName](self.parent, self.project, storeObject, objectName)
          
          if hasattr(interaction,'object') and interaction.object:
            ccpnParent = interaction.object
        
    if ccpnParent:
      self.ccpnParent = ccpnParent
    else:
      self.ccpnParent = self.project
    
    self.popups = {}

    #
    # Generically handle toggle stuff...
    #    
    
    self.toggleWidgetInfo = {}
    self.toggledWidgets = []      # Need to track these because can get regrid sometimes - should not happen if toggled
    
    #
    # Make an object if none available!
    #
    
    (objectIds,self.objectDict) = self.getObjectIdList()

    self.initialCreate = True
    
    if not objectIds:
      self.createObject()
      (objectIds,self.objectDict) = self.getObjectIdList()
      
      if not objectIds:
        showError('No valid object', 'You have to create a valid object before this window can be opened.')
        self.initialCreate = False

    if self.initialCreate:

      self.initialCreate = False
      TemporaryBasePopup.__init__(self, parent=parent, title="Project '%s': Edit %s" % (project.name,self.objectName), modal=False, transient=True)
    
    else:
      self.close()

  def body(self, master):
 
    self.selfMaster = master

    self.widgets = {}
    self.subClassWidgets = {}
    self.otherWidgets = {}

    self.varLists = {}
    self.varDicts = {}
    
    self.names = {'attrs': [], 'attrsList': [], 'attrsEnum': [], 'pulldowns': [], 'selectionLists': [], 'children': []}
    self.nameFuncs = {'Int': [], 'Boolean': [], 'Float': [], 'Double': [], 'Line': [], 'Text': []}

    self.columnspan = 3
    self.row = 0
 
    (objectIds,self.objectDict) = self.getObjectIdList()
    
    
    label = Label(master, text='Pick %s: ' % (self.objectName.lower()))
    label.grid(row=self.row, column=0, sticky=Tkinter.E)
    master.grid_columnconfigure(0,weight = 1)
    
    self.objectWidget = PulldownMenu(master, entries = objectIds, callback = self.updateCallback, do_initial_callback = False)
    self.objectWidget.grid(row=self.row, column=1, columnspan = (self.columnspan - 1), sticky=Tkinter.W)
    master.grid_columnconfigure(1,weight = 1)
    
    master.grid_rowconfigure(self.row,weight = 1)
    
    self.setObjectInformation(initialize = 1)

    self.row += 1
    
    texts = [ 'Create new' , 'Delete current' , 'Update current']
    commands = [ self.createObject, self.deleteObject, self.updateObject ]
    
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, dismiss_text = 'Cancel', help_url=self.help_url)
    buttons.grid(row=self.row, column=0, columnspan=self.columnspan)
    
  def setCurrentStores(self):
    
    for storeClassName in self.storeClassNameList:
     
      setCurrentStore(self.project,storeClassName)
    
  def setObjectInformation(self,initialize = 0, objectId = None):
  
    self.changedObject = None
    
    #
    # If an object was passed in at creation time, then select it straight away!
    #
    
    if initialize and self.selectedObject:
 
      object = self.selectedObject
      
      for objectId in self.objectDict.keys():
        if self.objectDict[objectId] == object:
          break
      
      self.setObjectWidget(objectId)
    

    #
    # Otherwise try to use object Id...
    #
      
    else:

      if not objectId:
        objectId = self.objectWidget.getSelected()
      else:
        self.setObjectWidget(objectId)

      object = self.objectDict[objectId]
      self.selectedObject = object
      self.metaClass = object._metaclass
    
    #
    # Subclass specific...
    #
    
    self.setObjectSpecificWidgets(object,initialize)

  def setAttrSelectionList(self,object,attrName,attrLabel,linkName,initialize,listParentSelect = None, listParent = None, forceDisplay = False, createNew = None, isChild = False, parentObjectLink = None): 
  
    # Only works if unique link...
    if parentObjectLink:
      object = getattr(object,parentObjectLink)
  
    if not listParent and not isChild:
      listParent = self.project
    
    if listParent:
      linkedObjects = getCcpnObjects(listParent,linkName)
      (self.varLists[attrName],self.varDicts[attrName]) = createSelection(linkedObjects)
      
    else:
      self.varLists[attrName] = None
      self.varDicts[attrName] = None

    # Don't display if nothing available...
    if not self.varLists[attrName] and not forceDisplay:
      return

    # Here have to decide whether to use 'sorted' method or direct (ordered) link
    selectedList = []  
    if object and attrName:
      metaAttr = object.metaclass.getAttribute(attrName)
      if metaAttr.isOrdered:
        objectList = list(getattr(object,attrName))
      else:
        attrNameUpper = attrName[0].upper() + attrName[1:]
        objectList = getattr(object,'sorted%s' % attrNameUpper)()
        
      if objectList:
        (selectedList,selectedDict) = createSelection(objectList)
      
            
    if initialize:
      
      if isChild:
        self.names['children'].append((attrName,attrName,parentObjectLink))
      else:
        self.names['selectionLists'].append((attrName,attrName,parentObjectLink))
      self.row += 1
      
      if object != self.project:
        labelText = 'Select relevant %s:' % attrLabel
      else:
        labelText = '%s list:' % (string.upper(attrLabel[0]) + attrLabel[1:])
        
      label = Label(self.selfMaster, text = labelText)
      label.grid(row=self.row, column=0, sticky=Tkinter.E)

      self.widgets[attrName] = ScrolledListbox(self.selfMaster,
                                         width = 50,
	                                        height = 5,
	                                        selectmode = Tkinter.MULTIPLE,
                                         initial_list = self.varLists[attrName],
                                         select_callback = self.changesMade)                                         
      
      columnspan = (self.columnspan-1)
      self.widgets[attrName].grid(row=self.row, column=1, columnspan = columnspan, sticky=Tkinter.W)
      
      #Implementation.registerNotify(self.resetWidget, "%s.%s" % (object._packageName,object.className), attrName)

      if createNew:
        self.row += 1
        self.makeCreationButton(attrName, createNew, 1, columnspan, listParent, listParentSelect)
   
    else:
      
      if attrName:
        self.widgets[attrName].clear()
      
        for varListItem in self.varLists[attrName]:
          self.widgets[attrName].append(varListItem)
      
      # Have to update the parent...
      if createNew:
        creationWidgetName = "Create%s" % createNew
        self.widgets[creationWidgetName].config(command = self.getCreationCommand(createNew,attrName,listParent,listParentSelect))
      
    if selectedList:
      self.widgets[attrName].setSelectedItems(selectedList)

  def setAttrPulldown(self,object,attrName,attrLabel,linkName,initialize, forceDisplay = False, createNew = None, listParentSelect = None, listParent = None, parentObjectLink = None):
  
    # Only works if unique link...
    if parentObjectLink:
      object = getattr(object,parentObjectLink)
 
    if not listParent:
      listParent = self.project

    linkedObjects = getCcpnObjects(listParent,linkName)

    (self.varLists[attrName],self.varDicts[attrName]) = createSelection(linkedObjects)
      
    # Don't display if nothing available...
    if not self.varLists[attrName] and not forceDisplay:
      return

    self.addNone(attrName)
    
    if object and getattr(object,attrName):
      (selectedList,selectedDict) = createSelection([getattr(object,attrName)])
      if selectedList:
        selectedIndex = self.varLists[attrName].index(selectedList[0])
    else:
      selectedIndex = 0
      
    if initialize:
      self.names['pulldowns'].append((attrName,attrName,parentObjectLink))
      self.row += 1
      
      if createNew:
        columnspan = 1
      else:
        columnspan = (self.columnspan - 1)
      
      if object != self.project:
        labelText = 'Select relevant %s:' % attrLabel
      else:
        labelText = '%s list:' % (string.upper(attrLabel[0]) + attrLabel[1:])

      self.widgets[attrName] = self.setSelectionList(labelText,columnspan,self.varLists[attrName],selectedIndex)
      
      if createNew:
        self.makeCreationButton(attrName,createNew,2,1,listParent,listParentSelect)
    
    else:
      # Avoid callback...
      callback = self.widgets[attrName].callback
      self.widgets[attrName].callback = None
      self.widgets[attrName].replace(self.varLists[attrName])
      self.widgets[attrName].setSelectedIndex(selectedIndex)
      self.widgets[attrName].callback = callback

      # Have to update the parent...
      if createNew:
        creationWidgetName = "Create%s" % createNew
        self.widgets[creationWidgetName].config(command = self.getCreationCommand(createNew,attrName,listParent,listParentSelect))

  def setNameFuncs(self,metaAttr,attrName):
    
    #
    # Keeps track of dataType...
    #
  
    dataTypeName = metaAttr.valueType.name
    
    if not self.nameFuncs.has_key(dataTypeName):
      self.nameFuncs[dataTypeName] = []
      
    if not attrName in self.nameFuncs[dataTypeName]:
      self.nameFuncs[dataTypeName].append(attrName)
    
  def setAttribute(self,object,attrName,attrLabel,initialize,subClassName = None):
    
    #
    # Get the correct metaAttribute
    #
    
    isCorrectSubClass = True
    
    if not subClassName or object.metaclass.name == subClassName:
      metaAttr = self.metaClass.getAttribute(name = attrName)
    else:
      isCorrectSubClass = False
      for subClass in self.metaClass.getSupertype().getAllSubtypes():
        if subClass.name == subClassName:
          metaAttr = subClass.getAttribute(name = attrName)
    
    metaDataType = metaAttr.valueType
    
    #
    # Set a widget name...
    #
    
    if subClassName:
      widgetName = "%s.%s" % (subClassName,attrName)
    else:
      widgetName = attrName
    
    if not self.otherWidgets.has_key(widgetName): 
      self.otherWidgets[widgetName] = []
    
    if not metaDataType.enumeration:
      if metaAttr.hicard == 1:

        #
        # Normal attribute
        #

        if object and isCorrectSubClass and getattr(object,attrName):
          value = getattr(object,attrName)
        else:
          value = ""

        if initialize:
        
          self.setNameFuncs(metaAttr,attrName)
            
          self.names['attrs'].append((widgetName,attrName,None))
          self.row += 1

          label = Label(self.selfMaster, text='%s: ' % (attrLabel))
          label.grid(row=self.row, column=0, sticky=Tkinter.E)
          self.otherWidgets[widgetName].append(label)

          self.widgets[widgetName] = Entry(self.selfMaster, text = value, width = 50,returnCallback = self.changesMade)
          self.widgets[widgetName].grid(row=self.row, column=1, columnspan = (self.columnspan - 1), sticky=Tkinter.W)

        else:

          self.widgets[widgetName].set(value)   
      
      else:
      
        #
        # List of attributes
        #
        
        if object and isCorrectSubClass:
          varList = getattr(object,attrName)
        else:
          varList = []

        if initialize:

          self.setNameFuncs(metaAttr,attrName)
          self.names['attrsList'].append((widgetName,attrName,None))
          self.row += 1

          label = Label(self.selfMaster, text = '%s:' % attrLabel)
          label.grid(row=self.row, column=0, sticky=Tkinter.E)
          self.otherWidgets[widgetName].append(label)

          self.widgets[widgetName] = ScrolledListbox(self.selfMaster,
                                             width = 50,
	                                            height = 5,
	                                            selectmode = Tkinter.MULTIPLE,
                                             initial_list = varList,
                                             select_callback = self.changesMade)                                         

          columnspan = (self.columnspan-1)
          self.widgets[widgetName].grid(row=self.row, column=1, columnspan = columnspan, sticky=Tkinter.W)

          self.row += 1

          button = Tkinter.Button(self.selfMaster, text = "Add value", command = lambda x = widgetName: self.addAttrListValue(x))
          button.grid(row=self.row, column=1, columnspan = 1,sticky=Tkinter.EW)
          self.otherWidgets[widgetName].append(button)

          button = Tkinter.Button(self.selfMaster, text = "Delete selected", command = lambda x = widgetName: self.deleteAttrListSelected(x))
          button.grid(row=self.row, column=2, columnspan = 1,sticky=Tkinter.EW)
          self.otherWidgets[widgetName].append(button)

        else:
          self.widgets[widgetName].clear()

          for varListItem in varList:
            self.widgets[widgetName].append(varListItem)
    
    else:

      #
      # Enumerated attribute
      #
      
      valueList = ['None'] + list(metaDataType.enumeration)

      if object and isCorrectSubClass and getattr(object,attrName):
        value = getattr(object,attrName)
      else:
        value = "None"
        
      selectedIndex = valueList.index(value)

      if initialize:
        self.setNameFuncs(metaAttr,attrName)
        self.names['attrsEnum'].append((widgetName,attrName,None))
        self.row += 1

        columnspan = (self.columnspan - 1)

        self.widgets[widgetName] = self.setSelectionList('%s (select from list):' % attrLabel,columnspan,valueList,selectedIndex)

        # TODO ALLOW ADDITION OF NEW VALUES?

      else:
        # Avoid callback...
        callback = self.widgets[widgetName].callback
        self.widgets[widgetName].callback = None
        self.widgets[widgetName].replace(valueList)
        self.widgets[widgetName].setSelectedIndex(selectedIndex)
        self.widgets[widgetName].callback = callback
    
    #
    # Hide subclasses that are not relevant...
    #
    
    if subClassName and object._metaclass.name != subClassName:
      self.widgets[widgetName].grid_remove()
      for widget in self.otherWidgets[widgetName]:
        widget.grid_remove()

    elif not initialize:
      widget = self.widgets[widgetName]
      
      if not widget in self.toggledWidgets:
        self.widgets[widgetName].grid()
        for widget in self.otherWidgets[widgetName]:
          widget.grid()
    
  def addAttrListValue(self,widgetName):
    
    self.widgets[widgetName].addItem()
  
  def deleteAttrListSelected(self,widgetName):
   
    self.widgets[widgetName].deleteItems()

  def updateObject(self, object = None):
    
    #
    # Updates the CCPN objects based on the info in the widgets
    #

    if not object:
      objectId = self.objectWidget.getSelected()
      object = self.objectDict[objectId]
  
    if object:
    
      self.metaClass = object._metaclass
      
      for attrType in ('attrs','pulldowns','attrsEnum','selectionLists','attrsList'):
      
        for (widgetName,attrName,parentObjectLink) in self.names[attrType]:
        
          setValue = True
          
          if parentObjectLink:
            parent = getattr(object,parentObjectLink)
          else:
            parent = object
          
          if attrType == 'attrs':
            value = self.widgets[widgetName].get()
        
            if attrName in self.nameFuncs['Int'] + self.nameFuncs['Boolean']:
              value = returnInt(value)
            elif attrName in self.nameFuncs['Float'] + self.nameFuncs['Double']:
              value = returnFloat(value)
            elif attrName in self.nameFuncs['Line'] + self.nameFuncs['Text'] and value == '':
              # Can't have empty strings - should be None
              value = None
              
          elif attrType == 'pulldowns':
            value = self.varDicts[attrName][self.widgets[widgetName].getSelected()]
              
          elif attrType == 'attrsEnum':
            value = self.widgets[widgetName].getSelected()
            if value == 'None':
              setValue = False
              
          elif attrType in 'selectionLists':
            value = []
            for attrLabel in self.widgets[widgetName].getSelectedItems():
              attrItem = self.varDicts[attrName][attrLabel]
              value.append(attrItem)
            
          elif attrType in 'attrsList':
            value = []
            for attrLabel in self.widgets[widgetName].getItems():
              value.append(attrLabel)
         
          if parent and setValue:
            setattr(parent,attrName,value)

    self.changedObject = None
    
  def addNone(self,attrName):
  
    self.varLists[attrName].insert(0,'None')
    self.varDicts[attrName]['None'] = None
    
  def setSelectionList(self,text,columnspan,selectionList,selectedIndex = 0):
    
    label = Label(self.selfMaster, text = text)
    label.grid(row=self.row, column=0, sticky=Tkinter.E)

    widget = PulldownMenu(self.selfMaster, entries = selectionList, selected_index = selectedIndex,callback = self.changesMade,do_initial_callback = False)
    widget.grid(row=self.row, column=1, columnspan = columnspan, sticky=Tkinter.W)
    
    return widget
    
  def setSeparator(self,initialize):
    
    if initialize:
      separator = Separator(self.selfMaster,height = 3)

      self.row += 1
      separator.grid(row = self.row, columnspan = self.columnspan, sticky = Tkinter.EW)
      
  def changesMade(self,*args,**kywds):

    self.changedObject = self.objectDict[self.objectWidget.getSelected()]

  def updateCallback(self,a,b):

    if self.changedObject:
      for objectId in self.objectDict.keys():
        if self.objectDict[objectId] == self.changedObject:
          break
      if (showYesNo('Update %s' % self.objectName.lower(), "Update changes to %s '%s' first?" % (self.objectName.lower(),objectId))):
        self.updateObject(object = self.changedObject)

    self.setObjectInformation(0)
  
  def updateWidgets(self, object = None):
     
    (objectIds,self.objectDict) = self.getObjectIdList()

    self.objectWidget.replace(objectIds)
    
    objectId = None
    
    if object:
      for tObjectId in objectIds:
        if self.objectDict[tObjectId] == object:
          objectId = tObjectId
          break
    
    self.setObjectInformation(0, objectId = objectId)
  
  def deleteObject(self):
  
    objectId = self.objectWidget.getSelected()
    
    if objectId:
    
      object = self.objectDict[objectId]

      if object:
        if (showYesNo('Delete %s' % self.objectName.lower(), "Really delete %s '%s'?" % (self.objectName.lower(),objectId))):
          object.delete()

      else:
        showInfo("No such %s" % self.objectName.lower(),"%s '%s' does not exist." % (self.objectName,objectId))
            
    self.updateWidgets()


  def setObjectWidget(self,objectId):
    
    # Avoid callback...
    
    callback = self.objectWidget.callback
    self.objectWidget.callback = None
    self.objectWidget.set(objectId)
    self.objectWidget.callback = callback

  def apply(self):
    
    # TODO anything here??
    
    return True
  
  def makeCreationButton(self,attrName,objectText,column,columnspan,ccpnParent,ccpnParentSelect):
  
    # TODO: need to have better system of tracking possible parent objects... or just select when creating if any choice?
    
    # TODO: Probably also have to track the link this is set for - automatically set the new object?!! Or NOT!??!
    
    command = self.getCreationCommand(objectText,attrName,ccpnParent,ccpnParentSelect)
    
    widgetName = "Create%s" % objectText
    self.widgets[widgetName] = Tkinter.Button(self.selfMaster, text = "Edit/new %s" % objectText, command = command)
    self.widgets[widgetName].grid(row=self.row, column=column, columnspan = columnspan,sticky=Tkinter.EW)
    
  def getCreationCommand(self,objectText,attrName,ccpnParent,ccpnParentSelect):

    return (lambda x = objectText, y = attrName, z = ccpnParent, a = ccpnParentSelect: self.createNewObject(x,y,z,a))
  
  def createNewObject(self,objectName,attrName,ccpnParent,ccpnParentSelect):
    
    self.updateObject()
    
    popup = None

    if self.popups.has_key(objectName):

      try:
        self.popups[objectName].open(doWait = False)
        popup = self.popups[objectName]
      except:
        pass
    
    if not popup:

      try:
        self.popups[objectName] = globals()['%sEditPopup' % objectName](self,self.project, ccpnParent = ccpnParent, ccpnParentSelect = ccpnParentSelect, createCallback = lambda x = objectName, y = attrName, createdObject = None: self.setOtherObject(x,y,createdObject))
      except:
        raise

    elif popup.ccpnParent != ccpnParent:
      popup.ccpnParent = ccpnParent
      
  def setOtherObject(self,objectName,attrName,createdObject):

    if createdObject:
      newOtherObject = createdObject
    elif self.popups.has_key(objectName):
      newOtherObject = self.popups[objectName].selectedObject
    else:
      newOtherObject = None
    
    if newOtherObject:
      
      objectId = self.objectWidget.getSelected()
      object = self.objectDict[objectId]
  
      for (widgetName,tAttrName,parentObjectLink) in self.names['pulldowns']:
        if attrName == tAttrName:
          setattr(object,attrName,newOtherObject)
          break
  
      for (widgetName,tAttrName,parentObjectLink) in self.names['selectionLists']:
        if attrName == tAttrName:
          setattr(object,'add%s' % objectName,newOtherObject)
          break
      
      for (widgetName,tAttrName,parentObjectLink) in self.names['children']:
        if attrName == tAttrName:
          setattr(object,'new%s' % objectName,newOtherObject)
          break
      
      if not self.initialCreate:
        self.updateWidgets()
    
  def getObjectIdList(self):

    (objectIds,objectDict) = self.setSpecificObjectIds()

    objectIds.sort()
    
    return (objectIds,objectDict)

  def createObject(self):
  
    popup = self.createSpecificObject()
    
    if hasattr(popup,'object'):
      object = popup.object
    else:
      object = None

    if not self.initialCreate:
      self.updateWidgets(object = object)

    if self.createCallback:
      self.createCallback(createdObject = object)
  
  #
  # All below is user-defined for specific object!!! See EntryEditPopup as example...
  #
          
  def setObjectSpecificWidgets(self,object,initialize):
  
    pass

  def setSpecificObjectIds(self):
    
    pass
    
  def setMetaClass(self):
  
    pass
  
  def createSpecificObject(self):
  
    pass

  def startToggle(self,toggleName):
    #print "Setting up toggle"
    if not self.toggleWidgetInfo.has_key(toggleName):
      self.row += 1      

      self.toggleWidgetInfo[toggleName] = {}
      self.toggleWidgetInfo[toggleName]['toggleWidgets'] = {}
      self.toggleWidgetInfo[toggleName]['startRow'] = self.row + 1
    
      self.toggleWidgetInfo[toggleName]['widget'] = ToggleLabel(self.selfMaster,  text='Edit %s' % toggleName, callback= lambda hidden, name = toggleName: self.toggleInfo(hidden,name))
      self.toggleWidgetInfo[toggleName]['widget'].grid(row=self.row, column=0, columnspan = 1, sticky = Tkinter.W)
    
      self.otherWidgets[toggleName] = self.toggleWidgetInfo[toggleName]['widget']

  def endToggle(self,toggleName):
    
    if not self.toggleWidgetInfo[toggleName].has_key('endRow'):
      self.toggleWidgetInfo[toggleName]['endRow'] = self.row
      #print "CB from endToggle"
      self.toggleWidgetInfo[toggleName]['widget'].callback(1)

  def toggleInfo(self,hidden,toggleName):
    #print "*** INITIALISING %s, %s***" % (toggleName, str(hidden))
    if hidden:
      gridWidgets = self.selfMaster.grid_slaves()
      toggleStartRow = self.toggleWidgetInfo[toggleName]['startRow']
      toggleEndRow = self.toggleWidgetInfo[toggleName]['endRow']
      
      #print "HIDDEN: %d, %d, %d" % (len(gridWidgets),toggleStartRow,toggleEndRow)
      
      for widget in gridWidgets:
        widgetInfo = widget.grid_info()
        widgetRow = int(widgetInfo['row'])
        
        if toggleEndRow >= widgetRow >= toggleStartRow:
          if not self.toggleWidgetInfo[toggleName]['toggleWidgets'].has_key(widget):
            self.toggleWidgetInfo[toggleName]['toggleWidgets'][widget] = [widgetRow,int(widgetInfo['column']),int(widgetInfo['rowspan']),int(widgetInfo['columnspan']),widgetInfo['sticky']]
            self.toggledWidgets.append(widget)
          #print "  Hiding: %s" % self.toggleWidgetInfo[toggleName]['toggleWidgets'][widget]
          widget.grid_forget()

    else:
      
      for widget in self.toggleWidgetInfo[toggleName]['toggleWidgets']:        
        (row,column,rowspan,columnspan,sticky) = self.toggleWidgetInfo[toggleName]['toggleWidgets'][widget]
        #print "  Showing: %s" % self.toggleWidgetInfo[toggleName]['toggleWidgets'][widget]
        widget.grid(row = row, column = column, rowspan = rowspan, columnspan = columnspan, sticky = sticky)

class ObjectCreatePopup(TemporaryBasePopup):

  def __init__(self, parent, project, ccpnParent, objectName, metaClass = None):

    self.name = None
    self.project = project
    self.parent = parent
    self.ccpnParent = ccpnParent
    self.metaClass = metaClass
    
    self.widgets = {}
    
    self.setLabels()
    
    self.setStoreInfo()
    
    setCurrentStore(self.project,self.storeClassName, linkName = self.storeLinkName)
    #print "INIT", objectName, ccpnParent
    TemporaryBasePopup.__init__(self, parent=parent, title="Project '%s': " % project.name + 'Create %s' % objectName, modal=True, transient=True)

  def body(self, master):
    
    columnspan = 2

    row = 0
    
    for (labelText,widgetName) in self.labels:
    
      label = Label(master, text= labelText)
      label.grid(row=row, column=0, sticky=Tkinter.W)
      self.widgets[widgetName] = Entry(master, text = '')
      self.widgets[widgetName].grid(row=row, column=1, sticky=Tkinter.W)
      row = row + 1

    row = self.setSpecificWidgets(row,columnspan,master)

    texts = [ 'Create' ]
    commands = [ self.ok ]
    
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, dismiss_text = 'Skip', help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan=columnspan)
  
  #
  # Object specific
  #
  
  def setStoreInfo(self):
  
    self.storeClassName = None
    self.storeLinkName = None
  
  def setLabels(self):
  
    pass
  
  def setSpecificWidgets(self,row,columnspan,master):
  
    pass
  
  def apply(self):
  
    pass

class AffiliationObjectCreatePopup(ObjectCreatePopup):
  
  def setStoreInfo(self):
    
    self.storeClassName = 'AffiliationStore'
    self.storeLinkName = None

class CitationObjectCreatePopup(ObjectCreatePopup):

  def setStoreInfo(self):
    
    self.storeClassName = 'CitationStore'
    self.storeLinkName = None

class NmrEntryObjectCreatePopup(ObjectCreatePopup):

  def setStoreInfo(self):
    
    self.storeClassName = 'NmrEntryStore'
    self.storeLinkName = None

#####################
#                   #
#  NmrEntry.Study  #
#                   #
#####################

class StudyEditPopup(ObjectEditPopup):

  help_url = joinPath(getHelpUrlDir(),'StudyEdit.html')
  objectName = 'Study'
  storeClassNameList = ['NmrProject','NmrEntryStore']
    
  def setObjectSpecificWidgets(self,study,initialize):
     
    #
    # Normal attributes
    #
    
    for (attrLabel,attrName) in (
                     ('Type of study','studyType'),
                     ('Study details','details'),
                     ('Study keywords','keywords')):

      self.setAttribute(study,attrName,attrLabel,initialize)
  
  def setSpecificObjectIds(self):
    
    studies = []
    
    if self.project.currentNmrEntryStore:
      studies = self.project.currentNmrEntryStore.sortedStudies()
    
    return createSelection(studies)
    
  def createSpecificObject(self):
  
    return StudyCreatePopup(self.parent,self.project,self.ccpnParent,self.objectName)

class StudyCreatePopup(NmrEntryObjectCreatePopup):

  help_url = joinPath(getHelpUrlDir(),'StudyCreate.html')
  
  def setLabels(self):
  
    self.labels = [('Name for new study','studyName')]
  
  def apply(self):

    self.name = string.strip(self.widgets['studyName'].get())

    if (not self.name):
      showError('No study name', 'Need to enter name')
      return False

    names = []
    
    if self.project.currentNmrEntryStore:
      names = [study.name for study in self.project.currentNmrEntryStore.sortedStudies()]

    if (self.name in names):
      showError('Repeated study name', 'Name already used')
      return False

    #
    # Create the new study...
    #
    
    self.object = NmrEntry.Study(self.ccpnParent, name = self.name)
    
    return True

#####################
#                   #
#  NmrEntry.Entry  #
#                   #
#####################

class EntryEditPopup(ObjectEditPopup):

  help_url = joinPath(getHelpUrlDir(),'EntryEdit.html')
  objectName = 'Entry'
  storeClassNameList = ['NmrProject','NmrEntryStore','AffiliationStore','CitationStore']
  
  def setObjectSpecificWidgets(self,entry,initialize):

    self.startToggle('Entry details')

    #
    # Normal attributes
    #
    
    for (attrLabel,attrName) in (
                  
                  ('Title for this entry','title'),
                  ('Entry details',       'details')):
    
      self.setAttribute(entry,attrName,attrLabel,initialize)
        
    #
    # Study
    #
     
    self.setAttrPulldown(entry,'study','study for this entry','currentNmrEntryStore.sortedStudies()',initialize, listParent = self.project.currentNmrEntryStore, forceDisplay = True, createNew = 'Study')
         
    self.endToggle('Entry details')

    self.setSeparator(initialize)
    
    self.startToggle('Affiliations')

    #
    # Laboratories
    #
     
    self.setAttrSelectionList(self.project.currentAffiliationStore,'organisations','organisations in entry','currentAffiliationStore.sortedOrganisations()',initialize, listParent = self.project.currentAffiliationStore, forceDisplay = True,  createNew = 'Organisation')
    self.setAttrSelectionList(entry,'laboratories','groups generating this entry','currentAffiliationStore.sortedOrganisations().sortedGroups()',initialize, listParent = self.project.currentAffiliationStore, listParentSelect = ('Organisation','project.currentAffiliationStore.sortedOrganisations()'), forceDisplay = True, createNew = 'Group')
    
    #
    # Contact persons
    #
     
    self.setAttrSelectionList(entry,'contactPersons','contact persons for entry','currentAffiliationStore.sortedPersons()',initialize, listParent = self.project.currentAffiliationStore, forceDisplay = True, createNew = 'Person')

    #
    # Authors
    #
     
    self.setAttrSelectionList(entry,'authors','entry authors','currentAffiliationStore.sortedPersons()',initialize, listParent = self.project.currentAffiliationStore, forceDisplay = True, createNew = 'Person')
 
    self.endToggle('Affiliations')
    
    self.setSeparator(initialize)

    self.startToggle('MolSystem and citations')

    #
    # MolSystem
    #
     
    self.setAttrPulldown(entry,'molSystem','molecular system','project.sortedMolSystems()',initialize)

    #
    # Citations
    #
     
    self.setAttrPulldown(entry,'primaryCitation','primary citations relevant for this entry','currentCitationStore.sortedCitations()',initialize, listParent = self.project.currentCitationStore, forceDisplay = True, createNew = 'Citation')
    self.setAttrSelectionList(entry,'otherCitations','other citations relevant for this entry','currentCitationStore.sortedCitations()',initialize, listParent = self.project.currentCitationStore, forceDisplay = True, createNew = 'Citation')
    
    self.endToggle('MolSystem and citations')
    
    self.setSeparator(initialize)

    self.startToggle('NMR and structure')

    #
    # Experiments
    #
     
    self.setAttrSelectionList(entry,'experiments','NMR experiments','currentNmrProject.sortedExperiments()',initialize, listParent = self.project.currentNmrProject)

    #
    # StructureGenerations
    #
     
    self.setAttrSelectionList(entry,'structureGenerations','structure generations','currentNmrProject.sortedStructureGenerations()',initialize, listParent = self.project.currentNmrProject)

    #
    # Measurement lists
    #
     
    self.setAttrSelectionList(entry,'measurementLists','measurement lists (e.g. shifts)','currentNmrProject.sortedMeasurementLists()',initialize, listParent = self.project.currentNmrProject)

    #
    # Derived data lists
    #
     
    self.setAttrSelectionList(entry,'derivedDataLists',"derived data lists (e.g. pKa's)",'currentNmrProject.sortedDerivedDataLists()',initialize, listParent = self.project.currentNmrProject)

    self.endToggle('NMR and structure')
    
    self.setSeparator(initialize)

    #
    # StructureGroups (part of MolSystem! TODO update when molSYstem changed!!)
    #
     
    #self.setAttrSelectionList(entry,'structureGroups','structureGroups',initialize,listParent = entry.molSystem)
    
    # TODO: relatedEntries? 
    
    # TODO: naturalSource - BUT NEED POPUPS TO SET THIS INFORMATION!!!
    
  def setSpecificObjectIds(self):

    return createSelection(self.project.currentNmrEntryStore.sortedEntries())
    
  def createSpecificObject(self):
  
    return EntryCreatePopup(self.parent,self.project,self.ccpnParent,self.objectName)

class EntryCreatePopup(NmrEntryObjectCreatePopup):

  help_url = joinPath(getHelpUrlDir(),'EntryCreate.html')
  
  def setLabels(self):
  
    self.labels = [('Name for new entry','entryName')]

  def apply(self):

    self.name = string.strip(self.widgets['entryName'].get())

    if (not self.name):
      showError('No entry name', 'Need to enter name')
      return False
    
    names = []
    
    if self.project.currentNmrEntryStore:
      names = [entry.name for entry in self.project.currentNmrEntryStore.sortedEntries()]

    if (self.name in names):
      showError('Repeated entry name', 'Name already used')
      return False

    #
    # Create the new entry...
    #
    
    self.object = NmrEntry.Entry(self.ccpnParent.currentNmrEntryStore, name = self.name)
    
    return True

########################
#                      #
#  Affiliation.Person  #
#                      #
########################

class PersonEditPopup(ObjectEditPopup):

  help_url = joinPath(getHelpUrlDir(),'PersonEdit.html')
  objectName = 'Person'
  storeClassNameList = ['AffiliationStore']

  def setObjectSpecificWidgets(self,person,initialize):
     
    #
    # Attribute lists
    #
    
    for (attrLabel,attrName) in (
                                 
          ('Initials','middleInitials'),
          ('Family title','familyTitle'),
          ('Title','title')):

      self.setAttribute(person,attrName,attrLabel,initialize)

    #
    # Person in group... TODO is this the best way to set children?!?!?
    #
    
    self.setAttrSelectionList(person,'personInGroups','groups where this person works','person.sortedPersonInGroups()',initialize, listParent = person, forceDisplay = True, createNew = 'PersonInGroup', isChild = True)
  
    #
    # Current person in group
    #
     
    self.setAttrPulldown(person,'currentPersonInGroup','current group that employs this person','person.sortedPersonInGroups()',initialize, listParent = person, forceDisplay = True)

    """
    TODO???
    authorCitations 	ccp.api.Citation.Citation 	0..* 	Citations where Person is author  
    editorCitations 	ccp.api.Citation.Citation 	0..* 	Citations where Person is editor  
    targets 	ccp.api.Target.Target 	0..* 	Targets created by Person.  
    expBlueprints 	ccp.api.ExpBlueprint.ExpBlueprint 	0..* 	ExpBlueprints owned by a Person.  
    createdProtocols 	ccp.api.Protocol.Protocol 	0..* 	Protocols created by Person.  
    editedProtocols 	ccp.api.Protocol.Protocol 	0..* 	Protocols last modified by Person.  
    createdExps 	ccp.api.Experiment.Experiment 	0..* 	Experiments created by Person  
    editedExps 	ccp.api.Experiment.Experiment 	0..* 	Experiments last modified by Person  
    """
   
  def setSpecificObjectIds(self):
  
    return createSelection(self.project.currentAffiliationStore.sortedPersons())
    
  def createSpecificObject(self):
  
    return PersonCreatePopup(self.parent,self.project,self.ccpnParent,self.objectName)

class PersonCreatePopup(AffiliationObjectCreatePopup):

  help_url = joinPath(getHelpUrlDir(),'PersonCreate.html')

  def setLabels(self):
  
    self.labels = [('Family name :','familyName'),
                   ('Given name :','givenName')]

  def apply(self):

    self.familyName = string.strip(self.widgets['familyName'].get())
    self.givenName = string.strip(self.widgets['givenName'].get())

    if (not self.familyName):
      showError('No person family name', 'Need to enter family name')
      return False

    for person in self.project.currentAffiliationStore.persons:

      if person.familyName == self.familyName and person.givenName == self.givenName:
        showError('Repeated person name', 'Name already used')
        return False

    #
    # Create the new person...
    #
    
    self.object = Affiliation.Person(self.ccpnParent, familyName = self.familyName, givenName = self.givenName)
    
    return True

##############################
#                            #
#  Affiliation.Organisation  #
#                            #
##############################

class OrganisationEditPopup(ObjectEditPopup):

  help_url = joinPath(getHelpUrlDir(),'OrganisationEdit.html')
  objectName = 'Organisation'
  storeClassNameList = ['AffiliationStore']
  
  def setObjectSpecificWidgets(self,organisation,initialize):
     
    #
    # Attribute lists
    #
    
    for (attrLabel,attrName) in (
                     
            ('Type of organisation (e.g. academic)','organisationType'),
            ('List of addresses',                   'addresses'),
            ('City',                                'city'),
            ('Postal code',                         'postalCode'),
            ('Country',                             'country'),
            ('Central phone number',                'phoneNumber'),
            ('Central fax number',                  'faxNumber'),
            ('List of email addresses',             'emailAddress'),
            ('URL for organisation',                'url')):

      self.setAttribute(organisation,attrName,attrLabel,initialize)
      
    #
    # Groups in organisation... TODO is this the best way to set children?!?!?
    #
    
    self.setAttrSelectionList(organisation,'groups','groups in this organisation','organisation.sortedGroups()',initialize, listParent = organisation, forceDisplay = True, createNew = 'Group', isChild = True)
   
  def setSpecificObjectIds(self):
  
    return createSelection(self.project.currentAffiliationStore.sortedOrganisations())
    
  def createSpecificObject(self):
 
    return OrganisationCreatePopup(self.parent,self.project,self.ccpnParent,self.objectName)

class OrganisationCreatePopup(AffiliationObjectCreatePopup):

  help_url = joinPath(getHelpUrlDir(),'OrganisationCreate.html')

  def setLabels(self):
  
    self.labels = [('Name of organisation :','name')]

  def apply(self):

    self.name = string.strip(self.widgets['name'].get())

    if (not self.name):
      showError('No organisation name', 'Need to enter name for organisation')
      return False

    for organisation in self.ccpnParent.organisations:

      if organisation.name == self.name:
        showError('Repeated organisation name', 'Name already used')
        return False

    #
    # Create the new organisation...
    #
    
    self.object = Affiliation.Organisation(self.ccpnParent, name = self.name)
    
    return True


#######################
#                     #
#  Affiliation.Group  #
#                     #
#######################

class GroupEditPopup(ObjectEditPopup):

  help_url = joinPath(getHelpUrlDir(),'GroupEdit.html')
  objectName = 'Group'
  storeClassNameList = ['AffiliationStore']

  def setObjectSpecificWidgets(self,group,initialize):
     
    #
    # Attribute lists
    #
    for (attrLabel,attrName) in (('Group URL','url'),):

      self.setAttribute(group,attrName,attrLabel,initialize)
   
  def setSpecificObjectIds(self):
  
    return createSelection(self.ccpnParent.sortedGroups())
    
  def createSpecificObject(self):
  
    return GroupCreatePopup(self.parent,self.project,self.ccpnParent,self.objectName)

class GroupCreatePopup(AffiliationObjectCreatePopup):

  help_url = joinPath(getHelpUrlDir(),'GroupCreate.html')

  def setLabels(self):
  
    self.labels = [('Name of group :','name')]

  def apply(self):

    self.name = string.strip(self.widgets['name'].get())

    if (not self.name):
      showError('No group name', 'Need to enter name for group')
      return False

    for group in self.ccpnParent.groups:

      if group.name == self.name:
        showError('Repeated group name', 'Name already used')
        return False

    #
    # Create the new group...
    #
    
    self.object = Affiliation.Group(self.ccpnParent, name = self.name)
    
    return True


###############################
#                             #
#  Affiliation.PersonInGroup  #
#                             #
###############################

class PersonInGroupEditPopup(ObjectEditPopup):

  help_url = joinPath(getHelpUrlDir(),'PersonInGroupEdit.html')
  objectName = 'PersonInGroup'
  storeClassNameList = ['AffiliationStore']

  def setObjectSpecificWidgets(self,personInGroup,initialize):
     
    #
    # Attribute lists
    #
    
    for (attrLabel,attrName) in (
                     
          ('Position (e.g. postdoc)',  'position'),
          ('Mailing address',          'mailingAddress'),
          ('Delivery address',         'deliveryAddress'),
          ('Email address',            'emailAddress'),
          ('List of phone numbers',    'phoneNumbers'),
          ('Fax number',               'faxNumber'),
          #('Probable end of contract', 'endDate')
          
          ):

      self.setAttribute(personInGroup,attrName,attrLabel,initialize)
    
  def setSpecificObjectIds(self):
  
    return createSelection(self.ccpnParent.sortedPersonInGroups())
    
  def createSpecificObject(self):
  
    return PersonInGroupCreatePopup(self.parent,self.project,self.ccpnParent,self.objectName)

class PersonInGroupCreatePopup(AffiliationObjectCreatePopup):

  help_url = joinPath(getHelpUrlDir(),'PersonInGroupCreate.html')

  def setLabels(self):
  
    self.labels = []
    
  def setSpecificWidgets(self,row,columnspan,master):
  
    organisations = self.project.currentAffiliationStore.organisations

    if organisations:
      
      groups = []
      
      for organisation in organisations:
        groups.extend(list(organisation.sortedGroups()))

      if groups:
        
        (groupIds,self.groupDict) = createSelection(groups)
        label = Label(master, text= 'Select group: ')
        label.grid(row=row, column=0, sticky=Tkinter.W)
        self.widgets['group'] =  PulldownMenu(master, entries = groupIds)
        self.widgets['group'].grid(row=row, column=1, sticky=Tkinter.W)

      else:

        label = Label(master, text= "No groups present for any organisations - press Skip and create one first.")
        label.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)
      
    else:
    
      label = Label(master, text= 'No organisations present - press Skip and create one first.')
      label.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)
 
    row = row + 1
   
    return row

  def apply(self):

    self.group = self.groupDict[self.widgets['group'].getSelected()]

    if (not self.group):
      showError('No group selected', 'Need to select a valid group to create a personInGroup')
      return False

    for personInGroup in self.ccpnParent.personInGroups:

      if personInGroup.group == self.group:
        showError('Repeated personInGroup', 'A personInGroup for this person and group already exists.')
        return False

    #
    # Create the new personInGroup...
    #
    
    self.object = Affiliation.PersonInGroup(self.ccpnParent, group = self.group)
    
    return True


#######################
#                     #
#  Citation.Citation  #
#                     #
#######################

class CitationEditPopup(ObjectEditPopup):

  help_url = joinPath(getHelpUrlDir(),'CitationEdit.html')
  objectName = 'Citation'
  storeClassNameList = ['CitationStore']

  def setObjectSpecificWidgets(self,citation,initialize):
     
    #
    # Attribute lists
    #
    
    for (attrLabel,attrName) in (
                                 
            ('Publication status',  'status'),
            ('First page number',   'firstPage'),
            ('Last page number',    'lastPage'),
            ('Year of publication', 'year')
            
            ):

      self.setAttribute(citation,attrName,attrLabel,initialize)

    #
    # Subclass stuff
    #

    for (attrLabel,attrName) in (
                                 
            ('Title of book',      'bookTitle'),
            ('Volume number',      'volume'),
            ('Publisher',          'publisher'),
            ('City of publisher',  'publisherCity'),
            ('Name of book series','bookSeries'),
            ('ISBN number',        'isbn')):

      self.setAttribute(citation,attrName,attrLabel,initialize,subClassName = 'BookCitation')

    for (attrLabel,attrName) in (
                                 
            ('Number of abstract',              'abstractNumber'),
            ('Title of conference',             'conferenceTitle'),
            ('Site where conference was held',  'conferenceSite'),
            ('City where conference was held',  'city'),
            ('State or province for city',      'stateProvince'),
            ('Country name',                    'country'),
            ('Start date of conference',        'startDate'),
            ('End date of conference',          'endDate'),
            
            ):

      self.setAttribute(citation,attrName,attrLabel,initialize,subClassName = 'ConferenceCitation')

    for (attrLabel,attrName) in (
                                 
            ('Abbreviated name for journal',    'journalAbbreviation'),
            ('Full name for journal',           'journalFullName'),
            ('Journal volume',                  'volume'),
            ('Journal issue',                   'issue'),
            ('ASTM number',                     'astm'),
            ('ISSN number',                     'issn'),
            ('CSD number',                      'csd')):

      self.setAttribute(citation,attrName,attrLabel,initialize,subClassName = 'JournalCitation')

    for (attrLabel,attrName) in (
                                 
            ('Institution where thesis was completed',  'institution'),
            ('City of institution',                     'city'),
            ('State or province',                      'stateProvince'),
            ('Country',                                'country')):

      self.setAttribute(citation,attrName,attrLabel,initialize,subClassName = 'ThesisCitation')
    
    #
    # Back to superclass..
    #
    
    for (attrLabel,attrName) in (
                                 
            ('General comments',      'details'),
            ('List of keywords',      'keywords'),
            ('PUBMED ID code',        'pubMedId'),
            ('CAS abstract code',     'casAbstractCode'),
            ('MEDLINE UI code',       'medlineUiCode')):

      self.setAttribute(citation,attrName,attrLabel,initialize)

    #
    # Authors
    #
     
    self.setAttrSelectionList(citation,'authors','authors of this citation','currentAffiliationStore.persons',initialize, listParent = self.project.currentAffiliationStore, forceDisplay = True, createNew = 'Person')    

    #
    # Editors
    #
     
    self.setAttrSelectionList(citation,'editors','editors of this citation (if any)','currentAffiliationStore.persons',initialize, listParent = self.project.currentAffiliationStore, forceDisplay = True, createNew = 'Person')

  def setSpecificObjectIds(self):
  
    return createSelection(self.ccpnParent.sortedCitations())
    
  def createSpecificObject(self):
  
    return CitationCreatePopup(self.parent,self.project,self.ccpnParent,self.objectName, metaClass = Citation.Citation._metaclass)

class CitationCreatePopup(CitationObjectCreatePopup):

  help_url = joinPath(getHelpUrlDir(),'CitationCreate.html')

  def setLabels(self):
  
    self.labels = [('Title of citation :','title')]

  def setSpecificWidgets(self,row,columnspan,master):
    
    subClassNames = []
    
    for subClass in self.metaClass.getAllSubtypes():
      scName = subClass.name
      subClassNames.append(scName)
      
    label = Label(master, text= 'Select citation type: ')
    label.grid(row=row, column=0, sticky=Tkinter.W)
    self.widgets['subType'] =  PulldownMenu(master, entries = subClassNames)
    self.widgets['subType'].grid(row=row, column=1, sticky=Tkinter.W)
 
    row = row + 1
   
    return row

  def apply(self):
  
    self.subType = self.widgets['subType'].getSelected()
    self.title = string.strip(self.widgets['title'].get())

    if (not self.title):
      showError('No citation title', 'Need to enter title for citation')
      return False

    for citation in self.ccpnParent.sortedCitations():

      if citation.title == self.title:
        showError('Repeated citation name', 'Name already used')
        return False

    #
    # Create the new citation...
    #
    
    self.object = getattr(Citation,self.subType)(self.ccpnParent, title = self.title)
    
    return True


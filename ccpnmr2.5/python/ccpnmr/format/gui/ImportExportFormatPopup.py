"""
======================COPYRIGHT/LICENSE START==========================

ImportExportFormatPopup.py: GUI window for selecting which files to import/export (replaces old Import/Export windows)

Copyright (C) 2005-2006 Wim Vranken (European Bioinformatics Institute)

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

import string, Tkinter, threading, time, os

from memops.universal.Util import returnInt, returnFloat # Do memopsWord type stuff as well?
from memops.universal.Io import joinPath, splitPath

from memops.gui.Label import Label
from memops.gui.CheckButton import CheckButton
from memops.gui.BasePopup import BasePopup
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo
from memops.gui.Separator import Separator
from memops.gui.ToggleLabel import ToggleLabel
from memops.gui.Entry import Entry
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledListbox import ScrolledListbox
from memops.gui.DataEntry import askDir
from memops.gui.ScrolledFrame import ScrolledFrame

from ccpnmr.format.gui.AcqProcParsEditPopup import AcqProcParsEditPopup
from ccpnmr.format.gui.InfoPopup import InfoPopup
from ccpnmr.format.gui.FormatFilePopup import FormatFilePopup

from ccpnmr.format.general.Io import getHelpUrlDir
from ccpnmr.format.general.Constants import IoSetupList
from ccpnmr.format.general.Constants import allFormatsDict
from ccpnmr.format.general.Util import createSelection

#
# Thread subclass to get return value
#

class Thread(threading.Thread):

  def run(self):
    self.returnValue = None
    
    if self.__target:
      self.returnValue = self.__target(*self.__args, **self.__kwargs)

#
# Handlers for the IoDefinitions...
#

def getCcpnObjects(parent,IOdefault):

  ccpnObjectStrings = string.split(IOdefault,'.')      
  ccpnObjects = [parent]

  for i in range(1,len(ccpnObjectStrings)):
    newCcpnObjects = []
    for ccpnObject in ccpnObjects:
      try:
        # To many link (e.g. project.currentNmrProject.sortedExperiments())
        newCcpnObjects.extend(list(reduce(getattr,ccpnObjectStrings[i:i+1],ccpnObject)))
      except:
        try:
          # Single attribute (e.g. strucGen.nmrConstraintStore)
          newCcpnObjects.append(getattr(ccpnObject,ccpnObjectStrings[i]))
        except:
          try:
            # Function (e.g. nmrConstraintStore.findAllConstraintLists(className = 'DistanceConstraintList'))
            newCcpnObjects.extend(list(eval("ccpnObject.%s" % ccpnObjectStrings[i])))
          except:
            pass

    ccpnObjects = newCcpnObjects[:]

  return ccpnObjects
  
def setIoInfo(master,setupInfo,project,fileDefs,forceCreation = False, appendNone = False):
  
  (widgetType,IOkeywd,labelText,IOdefault,extraDefs) = setupInfo
  
  label = Label(master, text = labelText)
  widget = None
  selectionDict = None
  
  if widgetType[-5:] == 'Entry':
    
    addKeywds = {}
    if extraDefs and extraDefs.has_key('width'):
      addKeywds['width'] = extraDefs['width']
      
    widget = Entry(master, text='%s' % str(IOdefault),**addKeywds)
    
  elif widgetType == 'CheckButton':
  
    widget = CheckButton(master)
    widget.setSelected(IOdefault)
    
  elif widgetType[:13] == 'SelectionList' or widgetType[:18] == 'MultiSelectionList':

    if IOdefault[:7] == 'project':
       
      ccpnObjects = getCcpnObjects(project,IOdefault)

      (selectionList,selectionDict) = createSelection(ccpnObjects)
      
      if (selectionList or forceCreation) and widgetType[-8:] == 'WithNone':
        if appendNone:
          selectionList.append('None')
        else:
          selectionList.insert(0,'None')
        selectionDict['None'] = None
      
    else:
      selectionList = IOdefault # Should be list!
    
    if selectionList:
      if widgetType[:13] == 'SelectionList':
        widget = PulldownList(master, texts = selectionList)
      elif widgetType[:18] == 'MultiSelectionList':
        widget = ScrolledListbox(master,height = 5,selectmode = Tkinter.MULTIPLE,initial_list = selectionList)

  elif widgetType in ('FileButton','FileMultiButton','DirButton'):
 
    (component,format,selectFile,selectDir) = fileDefs
    fileComponent = None
    
    if extraDefs:
      if extraDefs.has_key('component'):
        fileComponent = extraDefs['component'] # Can't reset component 'cos necessary for handling widget!
      elif extraDefs.has_key('format'):
        format = extraDefs['format']
    
    if widgetType[:4] == 'File':
      if widgetType == 'FileMultiButton':
        multiSelect= True
      else:
        multiSelect = False
      
      widget = Tkinter.Button(master, text = IOdefault, command = lambda comp = component, format = format, buttonKeyword = IOkeywd, defaultText = IOdefault, fileComponent = fileComponent, multiSelect = multiSelect: selectFile(comp,format,buttonKeyword,defaultText,fileComponent,multiSelect = multiSelect))

    elif widgetType == 'DirButton':
      widget = Tkinter.Button(master, text = IOdefault, command = lambda comp = component, format = format, buttonKeyword = IOkeywd: selectDir(comp,format,buttonKeyword))

  return (label,widget,selectionDict)

def getIoInfo(widget,setupInfo,project,selectionDict):

  (widgetType,IOkeywd,labelText,IOdefault,extraDefs) = setupInfo
  
  value = None
  
  if widgetType == 'IntEntry':
    
    value = returnInt(widget.get())
    
    if value == None:
      value = IOdefault
    
  elif widgetType == 'FloatEntry':
    
    value = returnFloat(widget.get())
    
    if value == None:
      value = IOdefault

  elif widgetType == 'StringEntry':
    
    # Get value directly so spaces are not stripped
    value = Tkinter.Entry.get(widget)

  elif widgetType == 'CheckButton':

    value = widget.isSelected()

  elif widgetType[:13] == 'SelectionList':
  
    value = widget.getText()
    
    if value and selectionDict:
      value = selectionDict[value]

  elif widgetType[:18] == 'MultiSelectionList':
  
    selectionLabels = widget.getSelectedItems()
    
    if selectionLabels:
      value = []
      
      if selectionDict:
        for selectionLabel in selectionLabels:
          value.append(selectionDict[selectionLabel])
      else:
        value = selectionLabels

  elif widgetType in ('FileButton','FileMultiButton','DirButton'):
  
    if widget.__getitem__('text'):
      value = widget.__getitem__('text')
  
  return value

class GenericFormatPopup(BasePopup):

  def __init__(self, parent, format, component = None, threading = True):
  
    self.project = parent.project
    self.guiParent = parent
    self.format = format 
    self.formatObjectDict = parent.formatObjectDict
    self.component = component
    self.threading = threading
    
    self.setImportExportFlag()
    
    self.files = {}
    self.function = {}
    self.components = []
    self.widgetSetup = {}
    self.toggleInfoIndex = {}

    #
    # Do the import of the relevant format class and the IOkeywords dictionary
    #
    
    self.formatLabel = allFormatsDict[self.format]
    self.formatModule = __import__('ccpnmr.format.converters.%sFormat' % self.formatLabel,{},{},['%sFormat' % self.formatLabel,'IOkeywords'])

    #
    # Definitions for import/export windows, files and functions (come from IoDefinitions)
    # 
    
    self.exportFunc = {}
    self.importFunc = {}
    self.IOkeywords = {}
    
    for (component,readFunc,writeFunc,fileDict,importSetup,exportSetup,exportFunc,importFunc) in IoSetupList:
      
      #
      # If component passed in, only do that one.
      #
      
      if self.component and self.component != component:
        continue
      
      #
      # Set import or export info
      #
      
      self.toggleInfoIndex[component] = 0

      if self.importExportFlag == 'import':
        self.widgetSetup[component] = importSetup[:]
        self.function[component] = readFunc
        componentAction = 'read'
        
        #
        # For multiple file import, don't bother with project files
        #
        
        if not self.component and component == 'project':
          continue
        
      else:
        self.widgetSetup[component] = exportSetup[:]
        self.function[component] = writeFunc
        componentAction = 'write'
      
      #
      # Set some variables and functions...
      #
      
      componentKeyWord = writeFunc[5:]        # This is also OK for import...
      
      self.exportFunc[component] = exportFunc
      self.importFunc[component] = importFunc
      
      #
      # Now set format specific information
      #
      
      if fileDict.has_key(format):
      
        (fileName,importFlag,exportFlag,addImport,addExport) = fileDict[format]
        
        if self.importExportFlag == 'import':
          validFlag = importFlag
          addInfo = addImport
        else:
          validFlag = exportFlag
          addInfo = addExport
        
        if validFlag:
        
          self.files[component] = fileName
          self.components.append(component)          
          self.IOkeywords[component] = {}
          
          #
          # Normal import/export flag
          #
          
          if validFlag == 1:

            #
            # Set fileName(s) info if available (is not in IOkeywords from Format files!!)
            #

            for widgetInfo in self.widgetSetup[component]:
              if widgetInfo[1] == 'fileName':
                self.IOkeywords[component]['fileName'] = (None,True,'The file name to be %sed.' % self.importExportFlag)
                break
              elif widgetInfo[1] == 'fileNames':
                self.IOkeywords[component]['fileNames'] = (None,True,'The file names to be %sed.' % self.importExportFlag)
                break

            #
            # Get IOkeywords info from Format files...
            #

            if self.formatModule.IOkeywords.has_key(self.function[component]):
              self.IOkeywords[component].update(self.formatModule.IOkeywords[self.function[component]])

              for superClass in self.formatModule.IOkeywords['subClasses']:
                if componentKeyWord in self.formatModule.IOkeywords['subClasses'][superClass]:
                  superClassFunc = componentAction + superClass
                  self.IOkeywords[component].update(self.formatModule.IOkeywords[superClassFunc])

            else:
              print "  Warning: no IOkeywords information for %s!" % self.function[component]
              self.IOkeywords[component] = {}

            #
            # Also add to self.toggleInfoIndex[component] if required keywords...
            #

            for IoKeyword in self.IOkeywords[component].keys():
              if self.IOkeywords[component][IoKeyword][1]:
                for widgetInfo in self.widgetSetup[component]:
                  if widgetInfo[1] == IoKeyword:
                    self.toggleInfoIndex[component] += 1

          #
          # Display text for import that refers to other mechanism...
          #
          
          else:
          
            self.IOkeywords[component]['showText'] = "Please use %s to %s %s." % (validFlag,self.importExportFlag,component)
            addInfo = None

        if addInfo:
                
          #
          # Now add format specific widget info and make sure it's all valid...
          #

          for addItem in addInfo:
            inputKeyWord = addItem[1]
            
            if inputKeyWord not in self.IOkeywords[component]:
              print "  Error: invalid %s keyword '%s' for %s!" % (self.importExportFlag,inputKeyWord,component)
              continue

            # Check whether keyword already exists! If so, replace...
            for widgetIndex in range(0,len(self.widgetSetup[component])):
              widgetInfo = self.widgetSetup[component][widgetIndex]
              if inputKeyWord == widgetInfo[1]:
                self.widgetSetup[component][widgetIndex] = addItem
                inputKeyWord = None

            # If the keyword wasn't found, then append to the full list.
            if inputKeyWord:
              self.widgetSetup[component].insert(self.toggleInfoIndex[component],addItem)
              self.toggleInfoIndex[component] += 1
                 
    # TODO: use labels here instead of name    
    BasePopup.__init__(self,parent = parent, title = self.titleText, modal = False, transient=False)
 
  def body(self, master):

    # Ensure that the first row and column in popup expand
    master.grid_rowconfigure(0, weight=1)
    master.grid_columnconfigure(0, weight=1, minsize=200)

    # add a scrolled frame for the widgets
    scrolledFrame = ScrolledFrame(master, xscroll=False)
    scrolledFrame.grid(row=0, column=0, sticky='nsew')
    _frame = scrolledFrame.frame

    if not self.components:
      showError("Error","No valid %s implemented (yet)." % self.importExportFlag,_frame)
      self.close()
      
    # 
    # Set specific import/export info
    #

    if self.importExportFlag == 'import' and not self.component:
      self.checkButton = {}
    else:
      self.exportButton = {}

    #
    # Generic info
    #
    
    self.widgets = {}
    self.widgetInfo = {}

    row = -1
    
    for component in self.components:

      # NOTE:ED component may be {} if buggy importer
      if not self.IOkeywords[component]:
        continue

      self.fileDefs = (component,self.format,self.selectFile,self.selectDir)
      normalWidgets = []
      toggleWidgets = []
      allWidgets =[]
      
      localWidgetCount = 0
      
      #
      # This is for handling rereferring text
      #
      
      if self.IOkeywords[component].has_key('showText'):
        
        row += 1
        
        label = Label(_frame, text = self.IOkeywords[component]['showText'])
        label.grid(row=row, column=0, columnspan = 3, sticky=Tkinter.EW)
        
        self.widgets[component] = None
        
        row = row + 1
        separator = Separator(_frame,height = 3)
        separator.setColor('black', bgColor = 'black')
        separator.grid(row=row, column=0, columnspan = 3, sticky=Tkinter.EW)

        continue

      #
      # Import ONLY!
      #
      
      if self.importExportFlag == 'import' and not self.component:
        row += 1
        
        label = Label(_frame, text = 'Import %s:' % component)
        label.grid(row=row, column=0, sticky=Tkinter.E)
      
        self.checkButton[component] = CheckButton(_frame) 
        self.checkButton[component].grid(row=row, column=1, sticky=Tkinter.W)
            
      toggleInfoLabelRow = None
      
      if self.importExportFlag == 'export':
        row += 1
        componentToggleRow = row
      else:
        componentToggleRow = None
        
      #row += 1

      for widgetIndex in range(0,len(self.widgetSetup[component])):

        widgetInfo = self.widgetSetup[component][widgetIndex]
        buttonKeyword = widgetInfo[1]

        (label,widget,selectionDict) = setIoInfo(_frame,widgetInfo,self.project,self.fileDefs)

        #
        # Make sure that widget info does not appear if no valid objects
        # available for mandatory information!!
        #

        if self.IOkeywords[component].has_key(widgetInfo[1]) and self.IOkeywords[component][widgetInfo[1]][1]:
          if not widget:
            if self.widgets.has_key(component):
              showError("Error","Could not set up export menu for %s: no '%s' information available." % (component,widgetInfo[3]),self)
              for (label,widget,infoButton) in normalWidgets + toggleWidgets:
                label.grid_forget()
                widget.grid_forget()
                infoButton.grid_forget()
                row -=1
                
              del(self.widgets[component])
            break

        #
        # This can fail (e.g. nothing available to make a list...)
        #

        if widget:
          if not self.widgets.has_key(component):
            self.widgets[component] = {}
            self.widgetInfo[component] = {}
          
          if localWidgetCount == self.toggleInfoIndex[component]:
            row += 1
            toggleInfoLabelRow = row
            allWidgets.append('TOGGLEWIDGET')

          self.widgets[component][buttonKeyword] = widget
          self.widgetInfo[component][buttonKeyword] = (widgetInfo,selectionDict)

          try:
            # infoButton = Tkinter.Button(_frame,text = 'i', font = ('Courier','10','bold'), width = 0, height = 0, command = lambda text = self.IOkeywords[component][buttonKeyword][2]: self.doInfoPopup(text))
            labelText = self.IOkeywords[component][buttonKeyword][2]
          except Exception as es:
            labelText = 'Label Error'

          infoButton = Tkinter.Button(_frame,text = 'i', font = ('Courier','10','bold'), width = 0, height = 0, command = lambda text = labelText: self.doInfoPopup(text))

          row = row + 1
          label.grid(row=row, column=0, sticky=Tkinter.E)
          self.widgets[component][buttonKeyword].grid(row=row, column=1, sticky=Tkinter.W)
          infoButton.grid(row=row, column=2, sticky=Tkinter.E)
        
          tempWidgetInfo = (label,self.widgets[component][buttonKeyword],infoButton)
          allWidgets.append(tempWidgetInfo)
          if self.widgetSetup[component].index(widgetInfo) >= self.toggleInfoIndex[component]:
            toggleWidgets.append(tempWidgetInfo)
          else:
            normalWidgets.append(tempWidgetInfo)

          localWidgetCount += 1

      if self.widgets.has_key(component) and self.widgets[component]:
        if self.importExportFlag == 'export':

          row += 1
          #print 'export button row %d' % row

          self.exportButton[component] = Tkinter.Button(_frame, text = "Export %s file." % component, command = lambda comp = component: self.importExportFile(comp))
          self.exportButton[component].grid(row=row, column=0, columnspan = 3, sticky=Tkinter.EW)
          
          allWidgets.append((None,self.exportButton[component],None))

        if toggleInfoLabelRow != None:
          #print toggleInfoLabelRow
          toggleInfoLabel = ToggleLabel(_frame,  text='Additional options', callback= lambda hidden, curRow = toggleInfoLabelRow, toggleWidgets = toggleWidgets: self.toggleInfo(hidden,curRow,toggleWidgets))
          toggleInfoLabel.grid(row=toggleInfoLabelRow, column=0, sticky = Tkinter.E)
 
          if componentToggleRow != None:
            # Keep track of widget for hiding when encompassing widget is expanded
            allWidgetIndex = allWidgets.index('TOGGLEWIDGET')
            allWidgets[allWidgetIndex] = (toggleInfoLabel,'TOGGLE',None)
          else:
            # If not part of encompassing toggle, just minimize now...
            toggleInfoLabel.callback(1)
        if componentToggleRow != None:
          toggleInfoLabel = ToggleLabel(_frame,  text='Click here for the %s export menu.' % component, callback= lambda hidden, curRow = componentToggleRow, allWidgets = allWidgets: self.toggleInfo(hidden,curRow,allWidgets))
          toggleInfoLabel.grid(row=componentToggleRow, column=0, columnspan = 3, sticky = Tkinter.EW)
          toggleInfoLabel.callback(1)

        row = row + 1
        separator = Separator(_frame,height = 3)
        separator.setColor('black', bgColor = 'black')
        separator.grid(row=row, column=0, columnspan = 3, sticky=Tkinter.EW)
      
    if self.widgets:
      row = row + 1
      
      if self.importExportFlag == 'import':
        texts = [ 'IMPORT' ]
        commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
      else:
        texts = []
        commands = []
      
      buttons = createDismissHelpButtonList(_frame, texts=texts, commands=commands, help_url=self.help_url)
      buttons.grid(row=row, column=0, columnspan = 3)

    else:
      showError("Error","No valid objects to %s" % self.importExportFlag,_frame)
      self.close()
  
  def doInfoPopup(self,text):
  
    popup = InfoPopup(self,'Information popup',text)
  
  def toggleInfo(self,hidden,curRow,toggleWidgets):
    #print "Current row:", curRow
    if hidden:
      for twidgets in toggleWidgets:
        for twidget in twidgets:
          if twidget != None and type(twidget) != type(''):
            twidget.grid_forget()

        #print "Hiding 1 row..."
        
    else:
      
      doToggle = None

      for i in range(0,len(toggleWidgets)):
        row = i + curRow + 1
        (label,widget,infoButton) = toggleWidgets[i]
        
        if label == None and infoButton == None:
          widget.grid(row = row, column = 0, columnspan = 3, sticky = Tkinter.EW)
          #print "LONG"
        else:
          label.grid(row=row, column=0, sticky=Tkinter.E)
          #print "LABEL",
          if widget:
            if widget == 'TOGGLE':
              doToggle = label
            else:
              #print "WIDGET",
              widget.grid(row=row, column=1, sticky=Tkinter.W)
              infoButton.grid(row=row, column=2, sticky=Tkinter.E)
          #print

        #print "Show %d - %s,%s,%s (%d)" % (row,label,widget,infoButton,curRow)

      if doToggle:
        doToggle.callback(1)

  def importExportFile(self,checkComponent):
    
    returnFlag = True
    importFuncs = []
    
    for component in self.components:
    
      #
      # Pass if no info...
      #
     
      if not self.IOkeywords[component]:
        continue
      
      #
      # If component passed in, only do that one
      #
      
      if checkComponent and checkComponent != component:
        continue
      
      else:

        #
        # Otherwise, if import, then only do so if checkButton is pressed
        #
        
        if self.importExportFlag == 'import' and not self.component:
          if not self.checkButton.has_key(component):
            continue
          elif not self.checkButton[component].isSelected():
            continue
    
      #
      # Get or set up format classes
      #

      if self.formatObjectDict[self.format]:

        formatObject = self.formatObjectDict[self.format]

      else:

        formatClass = getattr(self.formatModule,'%sFormat' % self.formatLabel)

        #
        # ...and create the format class instance
        #

        formatObject = formatClass(self.project,self.guiParent,verbose = 0)
        self.formatObjectDict[self.format] = formatObject

      #
      # Add input keyword information...
      #

      addArgs = []
      addKeywds = {}
      inputValid = True
      
      if self.widgets[component]:

        for inputKeyWord in self.widgets[component].keys():

          (widgetInfo,selectionDict) = self.widgetInfo[component][inputKeyWord]
          value = getIoInfo(self.widgets[component][inputKeyWord],widgetInfo,self.project,selectionDict)

          #
          # File name is a special case...
          #

          if inputKeyWord in ('fileName','fileNames'):

            self.fileName = self.fileNames = False

            if value == widgetInfo[3]:
              setattr(self,inputKeyWord,None)
            else:
              if inputKeyWord == 'fileName':
                self.fileName = value
              else:
                self.fileNames = self.decomposeMultipleFileText(value)

            addArgs.append(getattr(self,inputKeyWord))

          elif value != None:
            # Do not set as keyword for FileButtons that haven't been set!
            if widgetInfo[0] not in ('FileButton','FileMultiButton','DirButton') or widgetInfo[3] != value:
              if widgetInfo[0] == 'FileMultiButton':
                addKeywds[inputKeyWord] = self.decomposeMultipleFileText(value)
              else:
                addKeywds[inputKeyWord] = value

          else:
            if self.IOkeywords[component][inputKeyWord][1]:
              showWarning("Error","Nothing selected for %s writing, mandatory input %s!" % (component,inputKeyWord),self)
              inputValid = False
              break
      
      if not inputValid:
        break
            
      #
      # exportFunctions
      #

      if self.importExportFlag == 'export' and self.exportFunc[component]:

        for exportFuncLine in self.exportFunc[component]:

          if exportFuncLine == 'return':
            return

          if type(exportFuncLine) == type(''):
            eval(exportFuncLine)
          else:
            exportFuncLine()

      #
      # First check if a file name was selected?
      # Here add other checks for mandatory input to function?
      #

      if hasattr(self,'fileNames') and self.fileNames == False and not self.fileName:

       showWarning("Error","Please select a valid file name.",self)

      elif hasattr(self,'fileName') and self.fileName == False and not self.fileNames:

       showWarning("Error","Please select valid file names.",self)

      else:

        #
        # Do the read/write
        #
        
        rwFunc = getattr(formatObject,self.function[component])
        
        returnValue = None
        
        if self.threading:
          curThread = Thread(target = rwFunc,args = addArgs,kwargs = addKeywds)
          curThread.start()
          curThread.join()
          
          if curThread and curThread.returnValue:
            returnValue = curThread.returnValue
        
        else:
          returnValue = rwFunc(*addArgs,**addKeywds)
          
        #
        # Do some prep here...
        #
        
        fileNameText = ""
        
        if hasattr(self,'fileNames'):
          if self.fileNames == False:
            (tPathName,tFileName) = os.path.split(self.fileName)
            fileNameText = " file with name: %s." % tFileName
          else:
            fileNameText = " multiple files."
        
        #
        # Now show popup and set info depending on return value
        #

        if returnValue:
          
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
            
          if conversionSuccess:
            if fileNameText:
              showInfo("Success","Successfully %sed %s\n\n%s" % (self.importExportFlag,fileNameText,convertInfoString),self)
            
            returnFlag = True
            
            if self.importExportFlag == 'import':
              # TODO: this is a hack... needs to be more sophisticated than this (format selection,
              # even format dependent on particular resonance (!?!?!)
              if component != 'sequence':
                self.guiParent.resonanceFormat = self.format
  
              if self.importFunc[component]:
                for importFunc in self.importFunc[component]:
                  if importFunc not in importFuncs:
                    importFuncs.append(importFunc)

          else:
            if fileNameText:
              showInfo("Not finished","The %s did not finish correctly for %s\n\n%s" % (self.importExportFlag,fileNameText,convertInfoString),self)

        elif self.winfo_exists():
          if fileNameText:
            showError("Error","Error while %sing %s" % (self.importExportFlag,fileNameText),self)
            
    #
    # TODO: this might not be the best way to handle this!!
    #
    
    if returnFlag:
    
      if self.importExportFlag == 'import':
        if importFuncs:

          for importFuncLine in importFuncs:

            if type(importFuncLine) == type(''):
              eval(importFuncLine)
            else:
              importFuncLine()

        if hasattr(self.guiParent,'setMenuState'):
          self.guiParent.setMenuState()

    return returnFlag
   
  def selectFile(self,component,format,buttonKeyword,defaultText,fileComponent,multiSelect = False):
    
    fileName = self.widgets[component][buttonKeyword].__getitem__('text')
    
    #
    # Initialise file name, if required
    #
    
    if fileName == defaultText:

      if not fileComponent:
        fileName = self.files[component]
        fileComponent = component
      else:
        fileName = ''
        
    #
    # Now set up keywords
    #
    
    popupKeywds = {'component': fileComponent,'format': format, 'multiSelect': multiSelect}

    if multiSelect:
      files = self.decomposeMultipleFileText(fileName)
      popupKeywds['files'] = files
    
    else:
      popupKeywds['file'] = fileName
  
    if self.importExportFlag == 'import':
      popup = FormatFilePopup(self, **popupKeywds)
    else:
      popup = FormatFilePopup(self, selectionText = 'Select')
    
    # TODO: some testing on file name? Is it correct (ie has % code in for fid/ft bits?)
    
    if popup.fileSelected:
      
      if multiSelect:
        text = self.composeMultipleFileText(popup.files)
      else:
        text = popup.file
        
      self.widgets[component][buttonKeyword].config(text = text)
      
      if self.importExportFlag == 'import' and not self.component:
        self.checkButton[component].set(1)
    
    popup.destroy()
    
  def decomposeMultipleFileText(self,text):

    files = []
    
    if text:
      textItems = text.split("\n")
      if len(textItems) == 2:
        (path,fileNameText) = textItems
        for fileName in fileNameText.split(', '):
          files.append(joinPath(path,fileName))  
    
    return files

  def composeMultipleFileText(self,files):

    fileNames = []
    for file in files:
      (path,fileName) = splitPath(file)
      fileNames.append(fileName)
      
    text = "%s\n%s" % (path,', '.join(fileNames))
    
    return text

  def selectDir(self,component,format,buttonKeyword):
    
    dirName = self.widgets[component][buttonKeyword].__getitem__('text')
    
    title = 'Select directory for %s' % format
    prompt = 'Select %s directory for %s' % (component,format)

    newDirName = askDir(title,prompt,initial_value = dirName,parent = self)
    
    if newDirName:
      self.widgets[component][buttonKeyword].config(text = newDirName)
      
      if self.importExportFlag == 'import' and not self.component:
        self.checkButton[component].set(1)
  
  def ok(self):
  
    #
    # Redefined because windows can be destroyed while in action
    #
    
    if (not self.apply()):
      self.initial_focus.focus_set() # put focus back
      return
      
    if self.winfo_exists():
      self.close()
   
  def apply(self):
    
    returnValue = True
    
    if self.importExportFlag == 'import':
      
      returnValue = self.importExportFile(self.component) # If none, will do all components

      if returnValue and self.project.molSystems:
        
        if self.project.currentNmrProject.resonances or self.project.currentNmrProject.nmrConstraintStores:
          
          runLinkResonances = showYesNo('Run linkResonances','You have a defined molecular system and NMR atom assignment information in your project. To connect this information to each other you have to run linkResonances. Do you want to run this program now?')
          
          if runLinkResonances:
            self.parent.linkResonances()

    return returnValue

class ExportFormatPopup(GenericFormatPopup):

  help_url = joinPath(getHelpUrlDir(),'ExportFormat.html')
  
  def setImportExportFlag(self):
  
    self.importExportFlag = 'export'
    self.titleText = "Project '%s': " % self.project.name + '%s %s' % (self.format,self.importExportFlag)
    
class ImportFormatPopup(GenericFormatPopup):
 
  help_url = joinPath(getHelpUrlDir(),'ImportFormat.html')

  def setImportExportFlag(self):
  
    self.importExportFlag = 'import'
    self.titleText = "Project '%s': " % self.project.name + '%s multiple %s' % (self.format,self.importExportFlag)

class ImportFormatComponentPopup(GenericFormatPopup):
 
  help_url = joinPath(getHelpUrlDir(),'ImportFormatComponent.html')

  def setImportExportFlag(self):
  
    self.importExportFlag = 'import'
    self.titleText = "Project '%s': " % self.project.name + '%s %s %s' % (self.format,self.component,self.importExportFlag)

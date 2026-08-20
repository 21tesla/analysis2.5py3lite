
"""
======================COPYRIGHT/LICENSE START==========================

FormatConverter.py: Main GUI to format conversion functions

Copyright (C) 2005 Wim Vranken (European Bioinformatics Institute)

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
import string, sys, os, time

import Tkinter

#
# Import ccpn stuff
#

from memops.universal.Io import getPythonDirectory
from memops.universal.Io import joinPath
from memops.universal.Util import returnInt

from memops.general.Implementation import ApiError

from memops.general.Io import loadProject

from memops.editor.OpenProjectPopup import OpenProjectPopup
from memops.editor.SaveProjectPopup import SaveProjectPopup

import memops.api.Implementation as Implementation
from ccp.api.nmr import Nmr

#
# Get GUI stuff
#

from ccpnmr.format.gui.MultiDialog import multiDialog

from ccpnmr.format.general.Constants import IoSetupList
from ccpnmr.format.general.Constants import formatComponents
from ccpnmr.format.general.Constants import allFormatsDict
from ccpnmr.format.general.Constants import tagSep
from ccpnmr.format.general.Io import getHelpUrlDir
from ccpnmr.format.general.Io import TextPipe
from ccpnmr.format.general.Util import createSelection
from ccpnmr.format.general.Version import version

from ccpnmr.format.process.cleanProjectAppData import cleanAppData
from ccpnmr.format.process.cleanResonanceAppData import cleanResonanceAppData

from ccpnmr.format.gui.ImportExportFormatPopup import ImportFormatPopup
from ccpnmr.format.gui.ImportExportFormatPopup import ImportFormatComponentPopup
from ccpnmr.format.gui.ImportExportFormatPopup import ExportFormatPopup

from ccpnmr.format.gui.CreateShiftListPopup import CreateShiftListPopup
from ccpnmr.format.gui.WriteMappingPopup import WriteMappingPopup
from ccpnmr.format.gui.WriteBmrbChemShiftDepPopup import WriteBmrbChemShiftDepPopup
from ccpnmr.format.gui.AcqProcParsEditPopup import AcqProcParsEditPopup
from ccpnmr.format.gui.TextOutputPopup import TextOutputPopup

from memops.gui.BasePopup import BasePopup
from memops.gui.Menu import Menu
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.DataEntry import askString, askDir, askFile
from memops.gui.MessageReporter import showError, showYesNo, showInfo, showWarning
from memops.gui.WebBrowser import WebBrowser

ProjectMenu = 'Project'
SetupMenu = 'Setup'
ProcessMenu = 'Process'
ImportMenu = 'Import'
ExportMenu = 'Export'
OtherMenu = 'Other'
HelpMenu = 'Help'

class FormatConverter(BasePopup):

  applicationName = 'FormatConverter'

  def __init__(self, root, project = None, mem_cache = None, threading = True, *args, **keywds):

    #
    # Setup depends on whether project object passed in or not
    #   
    
    self.project = project
    
    if self.project:
      self.externalProject = True
    else:
      self.externalProject = False
    
    self.threading = threading
    
    self.multiDialog = multiDialog
    
    if self.project:
    
      initLocation = '310x10+100+100'
      
    else:
    
      initLocation = '380x10+100+100'
          
    self.mem_cache = mem_cache

    self.components = formatComponents
        	 
    #
    # Check which import/export formats available
    #
    
    self.importFormats = []
    self.exportFormats = []
    self.importExportFormats = []
    
    for infoList in IoSetupList:
    
      IoDefs = infoList[3]
    
      for formatName in IoDefs.keys():
      
        if IoDefs[formatName][1]:
          
          if formatName not in self.importFormats:
          
            self.importFormats.append(formatName)
          
        if IoDefs[formatName][2]:
        
          if formatName not in self.exportFormats:
          
            self.exportFormats.append(formatName)
            
        if IoDefs[formatName][1] and IoDefs[formatName][2]:
          
          if formatName not in self.importExportFormats:

            self.importExportFormats.append(formatName)

    self.importFormats.sort()
    self.exportFormats.sort()
    self.importExportFormats.sort()
	 	         
    BasePopup.__init__(self, parent=root, title='FormatConverter', location=initLocation, *args, **keywds)

  def body(self, master):
  
    self.webBrowser = WebBrowser(self.top, popup=self)

    #
    # Set up the text output window for standard output...
    #
    
    if not self.externalProject:

      self.textOutputPopup = TextOutputPopup(self, title = 'FormatConverter text output')
      self.textOutputPopup.withdraw()
    
      self.textPipe = TextPipe(self.textOutputPopup.text.text_area, guiWindow = self.textOutputPopup)
      sys.stdout = self.textPipe
    
    #
    # Set up rest of stuff...
    #

    self.menus = {}
    self.menu_items = {}
    self.popups = {}

    self.protocol('WM_DELETE_WINDOW', self.quit)

    self.menubar = Menu(self)

    self.setProjectMenu()
    
    self.initProject(project = self.project)
    
    self.setSetupMenu()

    self.setProcessMenu()

    self.setImportMenu()
    
    self.setExportMenu()
    
    self.setOtherMenu()
    
    self.setHelpMenu()

    self.setMenuState() # need to do it again because of OtherMenu state

    self.config(menu=self.menubar)


  """
  MENU SETUP
  """

  def setProjectMenu(self):
  
    menu = Menu(self.menubar, tearoff=0)

    if not self.externalProject:

      menu.add_command(label='New', shortcut='N', command=self.newProject)
      menu.add_command(label='Open', shortcut='O', command=self.openProject)
      menu.add_command(label='Close', shortcut='C', command=self.closeProject)
      menu.add_command(label='Save', shortcut='S', command=self.saveProject)
      menu.add_command(label='Save As', shortcut='A', command=self.saveAsProject)

    menu.add_command(label='Quit', shortcut='Q', command=self.quit)
    self.menubar.add_cascade(label=ProjectMenu, shortcut='P', menu=menu)
 
    self.menus[ProjectMenu] = menu
    
    menuItemList = []
    if not self.externalProject:
      menuItemList = ['New', 'Open', 'Close', 'Save', 'Save As']
    menuItemList.append('Quit')
    
    self.menu_items[ProjectMenu] = menuItemList

  def void(self):
  
    pass

  def setImportMenu(self):
  
    menu_items = []

    for componentType in self.components:
    
      menu_items.append({'kind': 'cascade','label': '%s files' % string.capitalize(componentType),'submenu': []})
      menu_item = menu_items[-1]

      for (compLabel,compName) in self.components[componentType]:
    
        #
        # First check if any import available at all!
        #
      
        setSubMenu = False

        for formatName in self.importFormats:
          for infoList in IoSetupList:
            if compName == infoList[0]:
              fileDict = infoList[3]

              if formatName in fileDict.keys():
                if fileDict[formatName][1]:
                  setSubMenu = True
                  break

        if setSubMenu:
          menu_item['submenu'].append({'kind': 'cascade','label': compLabel,'submenu': []})
          sub_menu_item = menu_item['submenu'][-1]
          
        else:
          continue

        #
        # Now set per format...
        #

        for formatName in self.importFormats:

          formatLabel = allFormatsDict[formatName]

          for infoList in IoSetupList:
            if compName == infoList[0]:
              fileDict = infoList[3]

              if formatName in fileDict.keys():
                if fileDict[formatName][1]:
                  sub_menu_item['submenu'].append({'kind': 'command',
                                                   'label': formatLabel,
                                                   'command': lambda format = formatName, comp = compName: self.importExportFormat('import',format,comp)})
  

    menu_items.append({'kind': 'cascade','label': 'Multiple files','submenu': []})

    for formatName in self.importFormats:
      
      formatLabel = allFormatsDict[formatName]

      menu_items[-1]['submenu'].append({'kind': 'command',
                                        'label': formatLabel,
                                        'command': lambda format = formatName: self.importExportFormat('import',format,"")})

    menu = Menu(self.menubar, menu_items = menu_items, tearoff=0)
    self.menubar.add_cascade(label=ImportMenu, shortcut='I', menu=menu)
 
    self.menus[ImportMenu] = menu
    self.menu_items[ImportMenu] = [menu_item['label'] for menu_item in menu_items]

  def setExportMenu(self):
  
    menu_items = []

    for formatName in self.exportFormats:
      
      formatLabel = allFormatsDict[formatName]

      menu_items.append({'kind': 'command',
                         'label': formatLabel,
                         'command': lambda format = formatName: self.importExportFormat('export',format,"")})

    menu = Menu(self.menubar, menu_items = menu_items, tearoff=0)
    self.menubar.add_cascade(label=ExportMenu, shortcut='E', menu=menu)
 
    self.menus[ExportMenu] = menu
    self.menu_items[ExportMenu] = [menu_item['label'] for menu_item in menu_items]
   
  def setSetupMenu(self):
  
    menu = Menu(self.menubar, tearoff=0)
    menu.add_command(label='Edit default experiment types', shortcut='E', command=self.editExperimentTypes)
    menu.add_command(label='Create 2D experiment', shortcut='2', command= lambda: self.createExperiment(2))
    menu.add_command(label='Create 3D experiment', shortcut='3', command= lambda: self.createExperiment(3))
    menu.add_command(label='Create 4D experiment', shortcut='4', command= lambda: self.createExperiment(4))
    self.menubar.add_cascade(label=SetupMenu, shortcut='S', menu=menu)
 
    self.menus[SetupMenu] = menu
    self.menu_items[SetupMenu] = ['Edit default experiment types',
                                  'Create 2D experiment',
                                  'Create 3D experiment',
                                  'Create 4D experiment']
  def setOtherMenu(self):
  
    menu = Menu(self.menubar, tearoff=0)
    # wb104: commented out 14 May 2012 because General Editor is broken
    ###menu.add_command(label='General Editor', command=self.runGeneralEditor)
    menu.add_command(label='Create/edit Entry', command=self.editEntry)
    menu.add_command(label='Edit acq/proc pars', command=self.editAcqProcPars)
    self.menubar.add_cascade(label=OtherMenu, shortcut='O', menu=menu)
 
    self.menus[OtherMenu] = menu
    self.menu_items[OtherMenu] = ['General Editor',
                                  'Create/edit Entry',
                                  'Edit acq/proc pars']

  def setProcessMenu(self):
  
    menu = Menu(self.menubar, tearoff=0)
    menu.add_command(label='Run linkResonances', shortcut='R', command=self.linkResonances)
    menu.add_command(label='Check stereo status prochirals', command=self.checkProchiralStatus)
    menu.add_command(label='Create chemical shifts from peaklist(s)', command= self.createShiftList)
    menu.add_command(label='Group chemical shifts for deposition', command= self.groupChemShifts)
    menu.add_command(label='Write resonance->atom mapping file', command= self.writeMapping)
    menu.add_command(label='Write BMRB chem shift deposition', command=self.writeBmrbCsDep)
    menu.add_command(label='Clean up application data in project', command=self.cleanAppData)
    menu.add_command(label='Clean up resonance application data', command=self.cleanResonanceAppData)
    self.menubar.add_cascade(label=ProcessMenu, shortcut='P', menu=menu)
 
    self.menus[ProcessMenu] = menu
    self.menu_items[ProcessMenu] = ['Run linkResonances',
                                    'Check stereo status prochirals',
                                    'Create chemical shifts from peaklist(s)',
                                    'Group chemical shifts for deposition',
                                    'Write resonance->atom mapping file',
	                                   'Write BMRB chem shift deposition',
                                    'Clean up application data in project',
                                    'Clean up resonance application data',
                                    ]

  def setHelpMenu(self):
  
    menu = Menu(self.menubar, tearoff=0)
    menuItems = ['Glossary','Menus','Quick Start','Tutorial','About','Version']
    
    for menuItem in menuItems:
    
      menu.add_command(label=menuItem, command= lambda helpItem = menuItem: self.help(helpItem))

    self.menubar.add_cascade(label=HelpMenu, shortcut='H', menu=menu)

    self.menus[HelpMenu] = menu
    self.menu_items[HelpMenu] = menuItems
    
  
  def setMenuState(self):

    # This one blanks certain options out if no project exists yet!!!
    
    if (self.project):
      state = Tkinter.NORMAL
    else:
      state = Tkinter.DISABLED
      
    #
    # Possible that this item not valid (if project object passed in)
    #

    if self.menus.has_key(ProjectMenu):
    
      menu = self.menus[ProjectMenu]
      menu_items = self.menu_items[ProjectMenu]
      for item in ['Close', 'Save', 'Save As']:
        try:
          n = menu_items.index(item)
          menu.entryconfig(n, state=state)
        except:
          pass

    menu = self.menus[SetupMenu]
    menu_items = self.menu_items[SetupMenu]
    for item in menu_items:
      try:
        n = menu_items.index(item)
        menu.entryconfig(n, state=state)
      except:
        pass

    menu = self.menus[ProcessMenu]
    menu_items = self.menu_items[ProcessMenu]
    for item in menu_items:
      try:
        n = menu_items.index(item)
        menu.entryconfig(n, state=state)
      except:
        pass

    menu = self.menus[ImportMenu]
    menu_items = self.menu_items[ImportMenu]

    for item in menu_items:
      try:
        n = menu_items.index(item)
        menu.entryconfig(n, state=state)
      except:
        pass

    menu = self.menus[ExportMenu]
    menu_items = self.menu_items[ExportMenu]
    for item in menu_items:
      try:
        n = menu_items.index(item)
        menu.entryconfig(n, state=state)
      except:
        pass

    menu = self.menus[OtherMenu]
    menu_items = self.menu_items[OtherMenu]
    for item in menu_items:
      try:
        n = menu_items.index(item)
        menu.entryconfig(n, state=state)
      except:
        pass
        
    # Exceptions for Other menu
        
    state = Tkinter.DISABLED
    
    if self.project:
    
      for exp in self.project.currentNmrProject.sortedExperiments():
        if exp.rawData:
          state = Tkinter.NORMAL
          break
        
    for item in ['Edit acq/proc pars']:
      try:
        n = menu_items.index(item)
        menu.entryconfig(n, state=state)
      except:
        pass

  """
  Other functions called from above
  """
  
  def help(self,helpItem):
  
    if helpItem == 'Version':
    
      showInfo('FormatConverter version','This is version %s' % version)
    
    else:

      urlFileDict = {'Glossary': 'Glossary',
                     'About': 'FormatConverterGeneral',
                     'Menus': 'MenuOverview',
                     'Quick Start': 'QuickStart',
                     'Tutorial': 'Tutorial'}

      url = joinPath(getHelpUrlDir(),urlFileDict[helpItem] + '.html')

      self.webBrowser.open(url)

  def setupFormatObjectDict(self):
  
    self.formatObjectDict = {}
      
    for formatName in self.importFormats + self.exportFormats:
    
      self.formatObjectDict[formatName] = None

  def importExportFormat(self, action, format, component):
  
    action = action

    popup = self.popups.get(action + format + component)
    
    # Remake every time...
    if (popup):
      popup.close()
      popup.destroy()
      
    if action == 'import':

      if component:
        self.popups[action + format + component] = ImportFormatComponentPopup(self, format, threading = self.threading, component = component)

      else:
        self.popups[action + format] = ImportFormatPopup(self, format, threading = self.threading)

    elif action == 'export':

      self.popups[action + format] = ExportFormatPopup(self, format, threading = self.threading)

  def editEntry(self):

    popupName = 'edit_entry'

    popup = self.popups.get(popupName)

    if (popup):
      try:
        popup.open()
      except:
        popup = None
        
    if not popup:
      try:
        self.popups[popupName] = self.multiDialog.EntryEdit(self, self.project)
      except:
        raise

  def editAcqProcPars(self):
  
    self.selectExperiment(hasAttrs = ['rawData'])
    
    if self.experiment:

      popupName = 'edit_pars_%s' % self.experiment.name
    
      popup = self.popups.get(popupName)
    
      if (popup):
        popup.open()
      else:
        self.popups[popupName] = AcqProcParsEditPopup(self, newFidDataSources = ("nmrPipe","azara"))
        
    else:
     
      showWarning(self,'No experiments with raw data available... try reading acquisition parameters.')       

  def selectExperiment(self, hasAttrs = None):
  
    self.experiment = None
 
    experimentList = list(self.project.currentNmrProject.sortedExperiments())

    if experimentList:
      
      for experiment in experimentList[:]:
        for attr in hasAttrs:
          if not hasattr(experiment,attr) or not getattr(experiment,attr):
            experimentList.pop(experimentList.index(experiment))
    
      if experimentList:
      
        if len(experimentList) == 1:
        
          self.experiment = experimentList[0]
        
        else:

          (selectionList,selectionDict) = createSelection(experimentList)

          interaction = self.multiDialog.SelectionList(

                               self,
                               selectionList,
                               selectionDict = selectionDict,
                               title = "Project '%s': " % self.project.name + 'Select experiment for acq/proc par editing',
                               text = "Existing experiment:",
                               dismissText = 'Dismiss',
                               urlFile = 'SelectExperimentForAcq')

          if interaction.isSelected:
            self.experiment = interaction.selection

  def initProject(self, project = None):
    
    self.project = project
    
    # Set up a dict of what component can be read from which format
    # Reset every time (necessary for formatObjectDict)
    self.setupFormatObjectDict()
    
    try:
      self.setMenuState()
    except: # fails first time around because OtherMenu not set up yet
      pass
      
    self.update()

  def newProject(self):

    if (self.project):
      if (not self.closeProject()):
        return

    name = askString(title='Project name', prompt='Enter project name:',
                     parent=self)

    if (name):
      project = Implementation.MemopsRoot(name=name)
      nmrProject = Nmr.NmrProject(project, name = project.name)
      self.initProject(project)
    
  def openProject(self):

    if (self.project):
      if (not self.closeProject()):
        return

    popup = self.popups.get('open_project')
    if (popup):
      popup.open()
    else:
      self.popups['open_project'] = OpenProjectPopup(self, callback = self.initProject)

  def checkSaveStatus(self):

    saveStatus = False
    
    if self.project:
      for topObject in self.project.topObjects:
        if topObject.isModified:
          saveStatus = True
          break
    
    return saveStatus

  def closeProject(self):

    if (not showYesNo('Close project', 'Close current project?')):
      return False
    if (self.checkSaveStatus()):
      if (showYesNo('Save project', 'Save changes to project?')):
        if (not self.saveFile()):
          return False

    self.project = None
    
    # Destroy live popups
    self.destroyPopups()
    
    # Destroy all other children (expect for text output window)
    for childWindow in self.children.values():
      if childWindow != self.textOutputPopup:
        
        # Set stopProcess if 'close()' def
        if childWindow.winfo_exists() and hasattr(childWindow,'close'):
          childWindow.stopProcess = True
    
    # Clean up any running imports
    formatClassFormats = self.formatObjectDict.keys()
    for format in formatClassFormats:
      del(self.formatObjectDict[format])
      
    return True

  def saveProject(self):
    
    self.project.saveModified()
    #projectDir = os.path.join(self.project.findFirstRepository(name = 'userData').url.path,'memops','Implementation')

    #if (os.path.exists(projectDir)):
    #else:
      # TODO THIS IS NOT CORRECT!
    #  self.askSaveFile()

  def saveAsProject(self):

    self.askSaveFile()

  def close(self, *event):

    if not self.externalProject:
      if (self.checkSaveStatus()):
        if (showYesNo('Save project', 'Save changes to project?')):
          if (not self.saveFile()):
            return False

    if (self.state() == 'normal'):
      self.withdraw()

    self.grab_release()

    self.var.set('')
    self.update_idletasks()
    
    if (self.parent):
      self.parent.focus_set()
      sys.stdout = sys.__stdout__

  def closeNoSave(self, *event):

    if (self.state() == 'normal'):
      self.withdraw()

    self.grab_release()

    self.var.set('')
    self.update_idletasks()
    
    if (self.parent):
      self.parent.focus_set()
      sys.stdout = sys.__stdout__

  def quit(self):

    if not self.externalProject:
      if (not showYesNo('Quit FormatConverter', 'Quit FormatConverter?')):
        return

      if self.checkSaveStatus():
        if (showYesNo('Save project', 'Save changes to project?')):
          if (not self.saveFile()):
            return False

      self.destroy()
      os._exit(0)
    else:
      self.closeNoSave()

  def saveFile(self):

    if (not self.project):
      return

    try:
      self.project.saveModified()
      return True
    except IOError, e:
      showError('Saving file', str(e))
      return False
    except ApiError, e:
      showError('Saving file', str(e))
      return False


  def askSaveFile(self):

    popup = self.popups.get('save_project')
    if (popup):
      popup.open()
    else:
      self.popups['save_project'] = SaveProjectPopup(self,project = self.project)

  def getProject(self):

    return self.project
  
  def destroyPopups(self):

    # TBD: what about other popups (i.e. created elsewhere)
    for key in self.popups.keys():
      popup = self.popups[key]
    
      # Stop process from children of popups (necessary?)
      for childWindow in popup.children.values():
        childWindow.stopProcess = True

      popup.stopProcess = True
      popup.close()
      popup.destroy()
      
    #self.textOutputPopup.text.text_area.clear()

    self.popups = {}

    if (self.project):
      self.project.delete()
      self.initProject(None)
    
    # Also delete dictionary for dataformat

  def linkResonances(self):
    
    #
    # Check what can be linked (Nmr resonances or constraint list info)
    #
    
    optionList = []
    nmrResonanceText = 'NMR project resonances'
    
    for resonance in self.project.currentNmrProject.resonances:
      if not resonance.resonanceSet:
        optionList.append(nmrResonanceText)
        break
        
    if self.project.currentNmrProject.nmrConstraintStores:
      
      unlinkedNcs = []
      for ncs in self.project.currentNmrProject.nmrConstraintStores:
        for resonance in ncs.fixedResonances:
          if not resonance.resonanceSet:
            unlinkedNcs.append(ncs)
            break
      
      (selectionList,selectionDict) = createSelection(unlinkedNcs)
      optionList.extend(selectionList)
    
    #
    # Ask user what to do...
    #
    
    option = None
    
    if optionList:
    
    
      if len(optionList) > 1:

        interaction = self.multiDialog.SelectionList(

                         self,
                         optionList,
                         title = "Project '%s': " % self.project.name + 'Select resonances to link',
                         text = "Select whether to link NMR project resonances or 'fixed' resonances from NMR constraint lists",
                         urlFile = ''

                         )
      
        if interaction.isSelected:
          option = interaction.selection
      
      else:
      
        option = optionList[0]
      
      #
      # Select nmrConstraintStore if chosen
      #
      
      if option:
      
        if option == nmrResonanceText:
          nmrConstraintStore = None
          
        else:
          nmrConstraintStore = selectionDict[option]
    
        #
        # Popup to select format if not chosen correctly, can also set all stereospecifically assigned
        #

        popup = self.multiDialog.LinkResonancesSetup(self,self.project,nmrConstraintStore)
    
    #
    # Say something at end...
    #

    if not option:
    
      showWarning(self,'No valid resonances or nothing selected - aborting linkResonances') 
  
  def checkProchiralStatus(self):
    
    popup = self.multiDialog.ProchiralStatus(self,self.project)

  def groupChemShifts(self):
    
    popup = self.multiDialog.GroupChemShifts(self,self.project)

  def createShiftList(self):
  
    popup = CreateShiftListPopup(self,self.project)  
  
  def writeMapping(self):
  
    popup = WriteMappingPopup(self,self.project)

  def writeBmrbCsDep(self):
  
    popup = WriteBmrbChemShiftDepPopup(self,self.project)

  def cleanAppData(self):
  
    if showYesNo('Clean application data', 'Remove all FormatConverter related application data? This will reduce the size of your project, but all information on previous import/export to external formats will be lost.'):
      
      cleanAppData(self.project)

  def cleanResonanceAppData(self):
    
    interaction = self.multiDialog.SelectionList(

                     self,
                     self.importExportFormats,
                     title = "Project '%s': " % self.project.name + 'Select format name for removing resonance tags',
                     text = 'Select the format name that you want to remove the Resonance application data tags for:',
                     urlFile = ''
                     )
  
    if interaction.isSelected:
      formatName = interaction.selection
  
  
      if showYesNo('Clean resonance application data', 'Remove all FormatConverter related Resonance linked application data for format %s?' % formatName):
       
        cleanResonanceAppData(self.project,formatName)

  def editExperimentTypes(self):
    
    popup = self.popups.get('edit_expTypes')

    if (popup):
      popup.open()
    else:
      self.popups['edit_expTypes'] = self.multiDialog.ExperimentTypesList(self,self.project)

  def createExperiment(self,numDim):

    popup = self.multiDialog.ExperimentCreate(self,self.project,numDim = numDim)

if __name__ == '__main__':

  inputArgs = sys.argv[1:]

  root = Tkinter.Tk()
  root.withdraw()
  
  #p = Implementation.Project(name = 'test')
  #top = FormatConverter(root, project = p)
  
  threading = False
  if inputArgs:
    # Redundant, but keep in anyway for backward compat
    for keywd in ('-nothread','-noThread'):
      if keywd in inputArgs:
        threading = False
        inputArgs.pop(inputArgs.index(keywd))

    for keywd in ('-thread','-useThread'):
      if keywd in inputArgs:
        threading = True
        inputArgs.pop(inputArgs.index(keywd))

  top = FormatConverter(root,threading = threading)
  top.open()

  # If file given, read in!
  if inputArgs:
    fileDir = inputArgs[0]
    
    askdir = lambda title, prompt, initial_value, default_dir: askDir(title, prompt,
        initial_value, parent=root, extra_dismiss_text='Skip', default_dir = default_dir)
    askfile = lambda title, prompt, initial_value: askFile(title, prompt,
        initial_value, parent=root, extra_dismiss_text='Skip')
    
    # TODO this probably has to be better again...
    top.initProject(loadProject(fileDir))
    
  else:
    top.initProject()

  root.mainloop()
 


"""
======================COPYRIGHT/LICENSE START==========================

ApplicationTemplate.py: <write function here>

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

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

===========================REFERENCE END===============================
"""
import os
import sys
import Tkinter

from memops.api                    import Implementation
from ccp.api.nmr import Nmr
from memops.universal.Io           import normalisePath
from memops.general.Implementation import ApiError
from memops.editor.BasePopup       import BasePopup

from memops.general.Application    import Application
from memops.general.Io             import loadProject
from memops.general.Util           import isProjectModified

from memops.gui.DataEntry       import askString, askDir, askFile
from memops.gui.FontMenu        import FontMenu
from memops.gui.Frame           import Frame
from memops.gui.Menu            import Menu
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo
from memops.gui.Util            import createDismissHelpButtonList

from memops.editor.OpenProjectPopup     import OpenProjectPopup
from memops.editor.SaveProjectPopup     import SaveProjectPopup
from memops.editor.PopupTemplate        import PopupTemplate
from memops.editor.WebBrowser            import ProjectWebBrowser

class ApplicationPopup(BasePopup):

  """
  **Application Template**

  This popup still needs its documentation.
  """

  help_url = 'http://www.ccpn.ac.uk'

  def __init__(self, root, programName='ApplicationTemplate'):

    self.programName = programName
    # Application object needed to store application-specific data with project
    self.application = Application(name=programName) 
    self.versionInfo = 'Version 0.00a'
    
    # Application popup is a superclass of memops.editor.BasePoup
    BasePopup.__init__(self, parent=root, title=programName,
                       location='+100+100', class_=self.application.name)
                       
  def body(self, guiParent):

    # Ensure that the first row and column in popup expand
    guiParent.grid_rowconfigure(0, weight=1)
    guiParent.grid_columnconfigure(0, weight=1, minsize=200)
    frame = Frame(guiParent) # Body widgets can be put in this frame
    frame.grid()

    buttons = createDismissHelpButtonList(frame, dismiss_text='Quit',
                dismiss_cmd=self.quit, help_url=self.help_url, expands=1)
    buttons.grid(row=0, column=0, sticky=Tkinter.NSEW)

    self.helpButton = buttons.buttons[-1]
    self.helpButtonCommand = self.helpButton.buttonCommand  # the command with no project

    # Dictionary to store popups opened by this application - e.g. need to close upon quit
    self.popups = {}
    
    # Default font
    self.font   = 'Helvetica 10'

    # Closing the window from the desktop environment calls the proper quit
    self.protocol('WM_DELETE_WINDOW', self.quit)

    self.mainMenu    = Menu(self)
    self.projectMenu = self.makeProjectMenu()
    self.otherMenu   = self.makeOtherMenu()
    
    # Put the main menu 
    self.config(menu=self.mainMenu)

    self.initProject()

  def makeProjectMenu(self):

    # Submenu of the min menu
    menu = Menu(self.mainMenu, tearoff=0)
    # Add various options to the menu and state the functions they call
    menu.add_command(label='New',     shortcut='N', command=self.newProject,
                     tipText='Make a new, blank CCPN project. Closes any existing project')
    menu.add_command(label='Open',    shortcut='O', command=self.openProject,
                     tipText='Open a new CCPN project from disk, closing any existing project')
    menu.add_command(label='Close',   shortcut='C', command=self.closeProject,
                     tipText='Close the current CCPN project')
    menu.add_command(label='Save',    shortcut='S', command=self.saveProject,
                     tipText='Save the current CCPN prjoject at its pre-specified location')
    menu.add_command(label='Save As', shortcut='A', command=self.saveAsProject,
                     tipText='Save the current CCPN project at a new location')
    menu.add_command(label='Quit',    shortcut='Q', command=self.quit,
                     tipText='Exit from the DANGLE program')
    menu.add_command(label='Version', shortcut='V', command=self.showVersion)
    self.mainMenu.add_cascade(label='Project', shortcut='P', menu=menu)

    menu.options = ['New','Open','Close','Save','Save As','Quit','Version']  
    return menu


  def makeOtherMenu(self):
    
    # The fonts menu is a pre-created widget
    fontsMenu = FontMenu(self.mainMenu, self.setFont, sizes=(8,10,12),
                         doItalic=0, doBoldItalic=0, tearoff=0)

    # Submenu of the min menu
    menu = Menu(self.mainMenu, tearoff=0)
    # Only Fonts option so far
    menu.add_command(label='Popup Template', shortcut='P', command=self.popupTemplate)
    menu.add_cascade(label='Fonts',          shortcut='F', menu=fontsMenu)
    self.mainMenu.add_cascade(label='Other', shortcut='O', menu=menu)

    menu.options = ['Popup Template','Fonts']
    return menu

  def popupTemplate(self):
  
    popup = self.openPopup('popup_template', PopupTemplate)

  def newProject(self):

    if self.project:
      # Project already present
      if not self.closeProject():
        # If we don't close the current project do nothing
        return

    name = askString(title='Project name', prompt='Enter project name:',parent=self)

    if name:
      # Make the API Project object
      project = Implementation.MemopsRoot(name=name)
      nmrProject = Nmr.NmrProject(project, name = project.name)
      self.initProject(project)

  def openProject(self):

    if self.project:
      if not self.closeProject():
        return

    self.openPopup('open_project', OpenProjectPopup, callback=self.initProject)


  def openPopup(self, popup_name, clazz, *args, **kw):

    popup = self.popups.get(popup_name)
  
    if popup:
      popup.open()
      
    else:
      transient = getAreDialogsTransient(self.project)
      # Below automatically opens popup
      popup = self.popups[popup_name] = clazz(self, project=self.project, transient=transient, *args, **kw)

    return popup


  def closeProject(self, queryClose = True, querySave = True):

    if queryClose:
      if (not showYesNo('Close project', 'Close current project?')):
        return False

    if querySave:
      if (not self.checkSaving()):
        return False

    self.helpButton.config(command=self.helpButtonCommand)
    self.destroyPopups()
    self.initProject()
 
    return True

  def saveProject(self):

    self.saveFile()
    
    #if (self.project.isStored):
    #  self.saveFile()
    #else:
    #  self.askSaveFile()

  def saveAsProject(self):

    self.askSaveFile()

  def showVersion(self):

    showInfo('Version', self.versionInfo, parent=self)


  def setFont(self, fontName=None, popup=None):

    if fontName is None:
      ###if self.project:
      ###  # Get font specification saved as project application data
      ###  self.font = getProjectFont(self.project)
      pass
    
    else:
      self.font = fontName
      ###if self.project:
      ###  # Set font specification as project application data
      ###  setProjectFont(self.project, self.font)
    
    if not popup:
      popup = self
      
    childList = popup.children.values()
    
    classes = [Tkinter.Button,
               Tkinter.Label,
               Tkinter.Menu,
               Tkinter.Entry,
               Tkinter.Checkbutton,
               Tkinter.Radiobutton]
    
    for child in childList:
      for clazz in classes:
        if isinstance(child, clazz):
          if hasattr(child, 'font'):
            if not child.font:
              child.config( font=self.font)
          else:
            child.config( font=self.font )
          break
    
      if hasattr(child, 'children'):
        childList.extend( child.children.values() )

  def initProject(self, project=None):

    self.project = project
    self.setTitle(self.programName)

    if project:
      self.project.application = self.application
      ###self.setFont()
      webBrowser = ProjectWebBrowser(self.top, popup=self, project=project)
      self.helpButton.config(command = lambda url=self.help_url: webBrowser.open(url) )
    
    self.setMenuState()
    
  def setMenuState(self):

    if self.project:
      state = Tkinter.NORMAL
    else:
      state = Tkinter.DISABLED

    # Disable bits of the project menu if there's no project
    for option in ('Save','Save As','Quit'):
      i = self.projectMenu.options.index(option)
      self.projectMenu.entryconfig(i, state=state)

    # Disable other manus if there's no project
    for menu in [self.otherMenu,]: # Include more menus in this list
      for i in range(len(menu.options)):
        menu.entryconfig(i, state=state)
      
      
  def askSaveFile(self):

    popup = self.openPopup('save_project', SaveProjectPopup)
    #popup.refresh()

  def modalAskSaveFile(self):

    popup = SaveProjectPopup(self, project=self.project, dismiss_text='Cancel Quit',
                             modal=True)
    did_save = popup.did_save
    popup.destroy()

    return did_save

  def quitSaveProject(self):

    if (self.project.activeRepositories):
      return self.saveFile()
    else:
      return self.modalAskSaveFile()

  def destroyPopups(self):

    for key in self.popups.keys():
      popup = self.popups[key]
      popup.destroy()

    self.popups = {}

  def destroyPopup(self, popup_name):

    if popup_name in self.popups.keys():
      popup =  self.popups.pop(popup_name)
      popup.destroy()

  def copyModifiedStorages(self):

    modifiedStorages = [storage for storage in self.project.storages if storage.isModified]
    if (self.project.isModified):
      modifiedStorages.append(self.project)

    for storage in modifiedStorages:
      self.copyStorage(storage)

  def saveFile(self):

    if (not self.project):
      return False

    try:
      self.copyModifiedStorages()
      self.project.saveModified()
      print 'Successfully saved project'
      return True
      
    except IOError, e:
      showError('Saving file', str(e))
      return False

  def checkSaving(self):

    if (not self.project):
      return True

    if isProjectModified(self.project):
      if showYesNo('Save project', 'Save changes to project?'):
        if not self.quitSaveProject():
          return False

    return True

  def quit(self):

    if not showYesNo('Quit ' + self.programName, 'Quit ' + self.programName + '?'):
      return

    if self.project:
      if not self.checkSaving():
        return

    self.destroy()

    # Need the below code to fix the problem with Bash
    # where the commend line was getting screwed up on exit.
    if os.name == 'posix':
      os.system('stty sane')
      
    sys.exit(0)

  def destroy(self):
  
    BasePopup.destroy(self)

def setAreDialogsTransient(project, areTransient):

  if (project):
    application = project.application
    application.setValue(project, keyword='areDialogsTransient',
                         value=areTransient, isBoolean=True)

def getAreDialogsTransient(project):

  if (project):
    """ Below does not work in v2
    application = project.application
    return application.getValue(project, keyword='areDialogsTransient',
                                defaultValue=True, isBoolean=True)
"""
    return True
  else:
    return True

""" Below does not work in v2
def getProjectFont(project):

  application = project.application
  return application.getValue(project, keyword='projectFont', defaultValue='Helvetica 10')

def setProjectFont(project, font):

  application = project.application
  application.setValue(project, keyword='projectFont', value=font)
"""

def launchApplication(filename=None):

  global top

  root = Tkinter.Tk()
  root.withdraw()
  top  = ApplicationPopup(root)
 
  project = None
  if filename:
    file    = normalisePath(filename)
    askdir  = lambda title, prompt, initial_value, default_dir: askDir(title, prompt,
              initial_value, parent=top, extra_dismiss_text='Skip', default_dir=default_dir)
    askfile = lambda title, prompt, initial_value: askFile(title, prompt,
              initial_value, parent=top, extra_dismiss_text='Skip')
    try:
      project = loadProject(path=file, showWarning=showWarning, askDir=askdir,
                            askFile=askfile)
    except ApiError, e:
      showError('Reading project', e.error_msg)
 
  top.update_idletasks()
  
  if project:
    top.initProject(project)

if (__name__ == '__main__'):

  argv = sys.argv[:]
  n    = len(argv)

  if n > 1:
    filename = argv[1]
  else:
    filename = None

  launchApplication(filename)

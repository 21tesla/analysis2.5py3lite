"""
======================COPYRIGHT/LICENSE START==========================

dataFileImportGui.py: Main GUI code to handle import of files for PDBe deposition

Copyright (C) 2010 Wim Vranken (PDBe, EBI)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.

A copy of this license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wim@ebi.ac.uk
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

import os, sys
import Tkinter

from memops.general.Application    import Application
from memops.general.Io             import loadProject, saveProject
from memops.general.Util           import isProjectModified

from memops.gui.DataEntry       import askString #, askDir, askFile
from memops.gui.Menu            import Menu
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo
from memops.gui.ScrolledMatrix  import ScrolledMatrix

from memops.editor.BasePopup          import BasePopup
from memops.editor.OpenProjectPopup   import OpenProjectPopup
from memops.editor.SaveProjectPopup   import SaveProjectPopup
from memops.editor.WebBrowser         import ProjectWebBrowser

from pdbe.deposition.dataFileImport.dataFileImportFrame import DataFileImportFrame
from pdbe.deposition.dataFileImport.formatConverterWrapper import FormatConverterWrapper

DEFAULT_FONT = 'Helvetica 10'

PROGRAM_NAME = 'File import for PDBe deposition'

VERSION = '0.9'

class DataFileImportGui(BasePopup):

  # TODO CHANGE
  help_url = 'http://www.ebi.ac.uk/pdbe'

  def __init__(self, root):

    self.root = root

    self.font = DEFAULT_FONT
    # Application object needed to store application-specific data with project
    self.application = Application(name=PROGRAM_NAME)
    self.versionInfo = VERSION
    
    ccpnProjectName = askString("CCPN project name",'Please give a name to this import session:')
    
    self.fcWrapper = FormatConverterWrapper(ccpnProjectName=ccpnProjectName,guiRoot=self)
    self.project = self.fcWrapper.formatConversion.ccpnProject
    
    # Application popup is a superclass of memops.editor.BasePoup
    BasePopup.__init__(self, parent=root, title=PROGRAM_NAME,
                       location='+100+100', class_=self.application.name)

    self.setTitle(PROGRAM_NAME)
                     
  def body(self, guiParent):

    self.geometry('600x350')

    # Ensure that the first row and column in popup expand
    guiParent.expandGrid(0,0)
    
    self.importFrame = DataFileImportFrame(guiParent, basePopup=self, grid=(0,0))
     
    # Dictionary to store popups opened by this application - e.g. need to close upon quit
    self.popups = {}
    
    # Default font
    self.font = DEFAULT_FONT

    # Closing the window from the desktop environment calls the proper quit
    self.protocol('WM_DELETE_WINDOW', self.quit)

    self.mainMenu    = Menu(self)
    self.projectMenu = self.makeProjectMenu()
    
    self.mainMenu.add_command(label='Formats',shortcut='F',command=self.importFrame.showFormats)
    self.mainMenu.add_command(label='Help',shortcut='H',command=self.importFrame.showMainInstructions)
   
    # Put the main menu 
    self.config(menu=self.mainMenu)

    self.setMenuState()

  def makeProjectMenu(self):

    # Submenu of the min menu
    menu = Menu(self.mainMenu, tearoff=0)
    menu.add_command(label='Save',    shortcut='S', command=self.saveProject)
    menu.add_command(label='Save As', shortcut='A', command=self.saveAsProject)
    menu.add_command(label='Quit',    shortcut='Q', command=self.quit)
    menu.add_command(label='Version', shortcut='V', command=self.showVersion)
    self.mainMenu.add_cascade(label='Project', shortcut='P', menu=menu)

    menu.options = ['Save','Save As','Quit','Version']
    return menu

  def saveProject(self):

    self.saveFile()

  def saveAsProject(self):

    self.askSaveFile()

  def showVersion(self):

    showInfo('Version', 'Version ' + self.versionInfo, parent=self)

  
  def setMenuState(self):

    state = Tkinter.NORMAL

    for option in ('Save','Save As','Quit','Version'):
      i = self.projectMenu.options.index(option)
      self.projectMenu.entryconfig(i, state=state)      
      
  def askSaveFile(self):

    popup = self.openPopup('save_project', SaveProjectPopup)
    popup.refresh()

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

  """
  def destroyPopups(self):

    for key in self.popups.keys():
      popup = self.popups[key]
      popup.destroy()

    self.popups = {}
  """

  def saveFile(self):

    if not self.project:
      return False

    try:
      saveProject(self.project, createFallback=True)

      print 'successfully saved project'
      return True
    except IOError, e:
      showError('Saving file', str(e), parent=self)
      return False

  def checkSaving(self):

    if not self.project:
      return True

    if isProjectModified(self.project):
      if showYesNo('Save project', 'Save changes to project?', parent=self):
        if not self.quitSaveProject():
          return False

    return True

  def quit(self):

    msg = 'Quit ' + PROGRAM_NAME

    if not showYesNo('Confirm', msg + '?', parent=self):
      return

    if self.project:
      if not self.checkSaving():
        return

    self.destroy()

    # Need the below code to fix the problem with Bash
    # where the commend line was getting screwed up on exit.
    if os.name == 'posix':
      os.system('stty sane')
      
    os._exit(0)

  def destroy(self):
  
    BasePopup.destroy(self)


def launchDataFileImport():

  global top

  root = Tkinter.Tk()
  root.withdraw()
  top  = DataFileImportGui(root)
 
  top.update_idletasks()
  
  root.mainloop()  

if __name__ == '__main__':

  launchDataFileImport()

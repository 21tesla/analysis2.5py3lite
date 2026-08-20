
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

from memops.universal.Io           import normalisePath

from memops.general.Implementation import ApiError
from memops.general.Application    import Application
from memops.general.Io             import loadProject, saveProject
from memops.general.Util           import isProjectModified

from memops.gui.DataEntry       import askString, askDir, askFile
from memops.gui.Button          import Button
from memops.gui.ButtonList      import ButtonList
from memops.gui.FontMenu        import FontMenu
from memops.gui.Frame           import Frame
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FileSelect      import FileType
from memops.gui.Label           import Label
from memops.gui.Menu            import Menu
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo
from memops.gui.TabbedFrame     import TabbedFrame
from memops.gui.Util            import createDismissHelpButtonList
from memops.gui.ScrolledMatrix  import ScrolledMatrix

from memops.editor.BasePopup          import BasePopup
from memops.editor.OpenProjectPopup   import OpenProjectPopup
from memops.editor.SaveProjectPopup   import SaveProjectPopup
from memops.editor.WebBrowser         import ProjectWebBrowser

from memops.universal.Io import getTopDirectory


# WMS

from cambridge.wms.RepositoryFrame import RepositoryFrame
from cambridge.wms.ExtendNmrFrame import ExtendNmrFrame
from cambridge.wms.ProjectFrame import ProjectFrame
from cambridge.wms.TasksFrame import TasksFrame
from cambridge.wms.DataExchFrame import DataExchFrame
from cambridge.wms.WorkflowFrame import WorkflowFrame
from cambridge.wms.ProtocolFrame import ProtocolFrame
from cambridge.wms.UserFrame import UserFrame

from cambridge.wms.LoginFrame import LoginPopup
from cambridge.wms.RepositoryPropertiesFrame import RepositoryPropertiesPopup

from cambridge.wms.RepositoryList import RepositoryList

from cambridge.wms.TaskManager import TaskManager

DEFAULT_FONT = 'Helvetica 10'

PROGRAM_NAME = 'Ccpn WMS'

VERSION = '0.0.1'

class ApplicationPopup(BasePopup):

  help_url = 'http://www.ccpn.ac.uk'

  def __init__(self, root, programName='Ccpn Workflow Management System'):

    self.programName = programName
    self.font = DEFAULT_FONT
    # Application object needed to store application-specific data with project
    self.application = Application(name=PROGRAM_NAME)
    self.versionInfo = 'Version ' + VERSION
    self.ariaProjectFile = None

    self.username = None

    # what do these do??
    self.updateFuncs = []
    self.projButtons = []

    # may need these in the end
    self.ariaPaths = []
    self.isdPaths = []

    # FIXME JMCI

    # hard code this for now. this should be initialised from file
    # read during GUI login and then extended using the menu options

    # should set a placeholder. In fact, this placeholder will also have
    # a default user (guest) so that there is always a structure to
    # extend should the user wish to start from scratch and not log in
    
    self.repList = RepositoryList()
    self.currentTask = None

    # going to stick in shortcuts to frames for quick access
    self.frameShortcuts = {}

    # Application popup is a superclass of memops.editor.BasePoup
    BasePopup.__init__(self, parent=root, title=PROGRAM_NAME,
                       location='+100+100', class_=self.application.name)

    self.setTitle(PROGRAM_NAME)

                       
  def body(self, guiParent):

    self.geometry('700x600')

    # Ensure that the first row and column in popup expand
    guiParent.grid_rowconfigure(0, weight=1)
    guiParent.grid_columnconfigure(0, weight=1, minsize=200)
    frame = Frame(guiParent) # Body widgets can be put in this frame
    frame.grid()

    # this is going to set up the options on the tabs. Really we need
    # stronger bookeeping between the options and the frames.  it
    # might make sense to define a master dictionary with option_name
    # => index. We can use this dictionary to build up the software
    # options and can also use it to define any refreshes etc (although
    # these would be better done through notifiers)

    
    softwareOpts = ['Repository','Project', 'Process','Client Tasks', 'Server Tasks', 'Data Exchange','Workflow','Protocols','User']

    # FIXME JMCI

    # we should have a tidier tabbed frame. ned to discuss rules with
    # tjs but it would seem to make sense to have the tabs stretch
    # across the entire top of the frame. Check the existing options
    # first, though

    self.tabbedFrame = TabbedFrame(guiParent, options=softwareOpts,
                                   toggleOff=False, selected=0)
    self.tabbedFrame.grid(row=0, column=0, sticky='nsew')


    frames = self.tabbedFrame.frames 

    # Logos
    ccpnDir = getTopDirectory()    
    imageDir = os.path.join(ccpnDir,'python','extendNmr','images')    
    imageFile = os.path.join(imageDir,'CcpnLogo.gif')
    self.ccpnLogo = Tkinter.PhotoImage(file=imageFile)

    # Dictionary to store popups opened by this application - e.g. need to close upon quit
    self.popups = {}
    
    # Default font
    self.font   = 'Helvetica 10'

    # Closing the window from the desktop environment calls the proper quit
    self.protocol('WM_DELETE_WINDOW', self.quit)

    # jmci we are going to have to change these drastically
    self.mainMenu    = Menu(self)
    self.userMenu = self.makeUserMenu()
    self.repositoryMenu = self.makeRepositoryMenu()
    self.otherMenu   = self.makeOtherMenu()
    
    # Put the top level menu
    self.config(menu=self.mainMenu)
    self.initRepository(frames[0])
    self.initProjectFrame(frames[1])
    self.initWms(frames[2])
    self.initTasks(frames[3])
    self.initDataExch(frames[5])
    self.initWorkflow(frames[6])
    self.initProtocol(frames[7])
    self.initUserFrame(frames[8])


    self.initProject()

    if not self.project:
      for button in self.projButtons:
        button.disable()

    # the main client side daemon. This is going to run in the background
    # and control the task processing
    self.taskDaemon = TaskManager(self)
    self.taskDaemon.active = False
    self.taskDaemon.start()


  ## THE TABBED FRAMES ##

  def initRepository(self, frame):

    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
     
    repositoryFrame = RepositoryFrame(frame, basePopup=self)
    repositoryFrame.grid(row=1, column=0, sticky='nsew')

    self.updateFuncs.append(repositoryFrame.updateAll)
    

  def initProjectFrame(self, frame):

    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
    
    self.projectFrame = ProjectFrame(frame, basePopup=self)
    self.projectFrame.grid(row=1, column=0, sticky='nsew')

    self.updateFuncs.append(self.projectFrame.updateAll)
    
  def initWms(self, frame):

    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
     
    extendNmrFrame = ExtendNmrFrame(frame, basePopup=self)
    extendNmrFrame.grid(row=1, column=0, sticky='nsew')

    self.updateFuncs.append(extendNmrFrame.updateAll)
    

  def initTasks(self, frame):

    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
     
    tasksFrame = TasksFrame(frame, self)
    tasksFrame.grid(row=1, column=0, sticky='nsew')

    self.updateFuncs.append(tasksFrame.updateAll)


  def initDataExch(self, frame):

    row = 0
    frame.config(bg='#E6E6E6')

    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)

    dataExchFrame = DataExchFrame(frame, basePopup=self)
    dataExchFrame.grid(row=1, column=0, sticky='nsew')

    self.updateFuncs.append(dataExchFrame.updateAll)
    
  def initWorkflow(self, frame):

    row = 0
    frame.config(bg='#FFFFFF')

    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)

    workflowFrame = WorkflowFrame(frame, basePopup=self)
    workflowFrame.grid(row=1, column=0, sticky='nsew')

    self.updateFuncs.append(workflowFrame.updateAll)

  def initProtocol(self, frame):

    row = 0
    frame.config(bg='#FFFFFF')

    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)

    protocolFrame = ProtocolFrame(frame, basePopup=self)
    protocolFrame.grid(row=1, column=0, sticky='nsew')

    self.updateFuncs.append(protocolFrame.updateAll)

  def initUserFrame(self, frame):

    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)

    row = 0
    frame.config(bg='#FFFFFF')

    userFrame = UserFrame(frame, basePopup=self)
    userFrame.grid(row=1, column=0, sticky='nsew')

    self.updateFuncs.append(userFrame.updateAll)
    

  ## The Menu ##

  def makeUserMenu(self):
    menu = Menu(self.mainMenu, tearoff=0)
    menu.add_command(label='Log In',    shortcut='I', command=self.openLoginPopup)
    menu.add_command(label='Log Out',    shortcut='O', command=self.logout)
    menu.add_command(label='Save Profile',    shortcut='S', command=self.tmpCall)
    menu.add_command(label='View Profile',    shortcut='V', command=self.tmpCall)
    menu.add_command(label='Create Profile',    shortcut='C', command=self.tmpCall)
    menu.add_command(label='Edit Profile',    shortcut='E', command=self.tmpCall)
    self.mainMenu.add_cascade(label='User', shortcut='U', menu=menu)
    menu.options = ['Log In','Log Out','Save Profile','View Profile','Create Profile', 'Edit Profile']
    return menu 
    
  def makeRepositoryMenu(self):
    self.repMenu = Menu(self.mainMenu, tearoff=0)
    # Add various options to the menu and state the functions they call
    self.repMenu.add_command(label='Manager',    shortcut='M', command=self.tmpCall)
    self.repMenu.add_command(label='Connect',    shortcut='C', command=self.openRepository)
    
    # doesn't really work to have a disconnect; we do not have a list
    # of repositories and so would it refer to? Better to do this from
    # within the RepositoryFrame

    self.repMenu.add_command(label='Disconnect',    shortcut='D', command=self.tmpCall)

    # again, this is functionality that would make sense but which it is hard to put
    # in the top level menu because it is not obvious which repository it refers
    # to

    self.repMenu.add_command(label='Import Project', shortcut='I', command=self.tmpCall)
    self.repMenu.add_command(label='Export Project', shortcut='E', command=self.tmpCall)


    # FIXME JMCI
    
    # the quit option needs to be tied to a dialogue box that asks whether
    # the current connection state should be saved. This will have to be
    # secured and so means that the GUI should have a login itself (beyond
    # the logins for the various resources).
    
    self.repMenu.add_command(label='Quit',       shortcut='Q', command=self.quit)

    self.mainMenu.add_cascade(label='Repository', shortcut='R', menu=self.repMenu)
    self.repMenu.options = ['Connect','Disconnect','Import Project','Export Project','Quit']
    return self.repMenu

    
  def makeProjectMenu(self):

    # Submenu of the min menu
    menu = Menu(self.mainMenu, tearoff=0)
    # Add various options to the menu and state the functions they call
    menu.add_command(label='New',     shortcut='N', command=self.newProject)
    menu.add_command(label='Open',    shortcut='O', command=self.openProject)
    menu.add_command(label='Close',   shortcut='C', command=self.closeProject)
    menu.add_command(label='Save',    shortcut='S', command=self.saveProject)
    menu.add_command(label='Save As', shortcut='A', command=self.saveAsProject)
    menu.add_command(label='Quit',    shortcut='Q', command=self.quit)
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
    menu.add_cascade(label='Fonts',          shortcut='F', menu=fontsMenu)
    self.mainMenu.add_cascade(label='Options', shortcut='O', menu=menu)

    menu.options = ['Fonts',]
    return menu

  def makeCcpNmrMenu(self):

    try:
      from ccpnmr.analysis.AnalysisPopup import AnalysisPopup
      haveAnalysis = True
    except ImportError:
      haveAnalysis = False

    try:
      from ccpnmr.format.gui.FormatConverter import FormatConverter
      from ccpnmr.format.gui.DataShifter import DataShifter
      haveFormatConv = True
    except ImportError:
      haveFormatConv = False
  
    menu = Menu(self.mainMenu, tearoff=0)
    
    if haveAnalysis:
      menu.add_command(label='Analysis',     shortcut='A', command=self.ccpNmrAnalysis)
    else:
      menu.add_command(label='Analysis ** NOT INSTALLED **',     shortcut='A', command=None)
    
    if haveFormatConv:
      menu.add_command(label='Format Converter',  shortcut='F', command=self.ccpNmrFormatConverter)
      menu.add_command(label='Data Shifter',    shortcut='D', command=self.ccpNmrDataShifter)
    else:
      menu.add_command(label='Format Converter ** NOT INSTALLED **',  shortcut='F', command=None)
      menu.add_command(label='Data Shifter ** NOT INSTALLED **',    shortcut='D', command=None)
 
    self.mainMenu.add_cascade(label='CcpNmr', shortcut='C', menu=menu)

    menu.options = ['Analysis','Format Converter','Data Shifter']  
    return menu
  
  def ccpNmrAnalysis(self, cacheSize=64):
  
    from ccpnmr.analysis.AnalysisPopup import AnalysisPopup
    analysis = AnalysisPopup(self, cache_size=cacheSize, glDirect=None)
    analysis.update_idletasks() # fast
    if self.project:
       analysis.initProject(self.project)
       
  def ccpNmrFormatConverter(self):
  
    from ccpnmr.format.gui.FormatConverter import FormatConverter
    fc = FormatConverter(self,threading=False,project=self.project)
    fc.open()
    if self.project:
      fc.initProject(self.project)
      
  def ccpNmrDataShifter(self):
    
    from ccpnmr.format.gui.DataShifter import DataShifter
    ds = DataShifter(self)
    ds.open()

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
      nmrProject = project.newNmrProject(name=name)
      self.initProject(project)

  # placeholder method until we write code

  def tmpCall(self):

    return


  def openProject(self):

    if self.project:
      if not self.closeProject():
        return

    self.openPopup('open_project', OpenProjectPopup,
                   callback=self.initProject)

  def openLoginPopup(self):

    # at the moment we do this through methods in the popup. It
    # would be much more sensible to do this through some
    # sort of callback. 

    login_popup = LoginPopup(self)

    # really we want to imoblise the rest. how to do?


  def setRepository(self, rep):

    self.repList.currentRepository = rep

    if self.frameShortcuts.has_key('Repository'):
      self.frameShortcuts['Repository'].drawFrame()

    if self.frameShortcuts.has_key('Tasks'):
      self.frameShortcuts['Tasks'].drawFrame()
    #for ff in self.tabbedFrame.frames[0].children.values():
    #  ff.draw()    

    #for ff in self.tabbedFrame.frames[3].children.values():
    #  ff.draw()

    for ff in self.tabbedFrame.frames[8].children.values():
      ff.draw()

  def openRepository(self):

    rep_popup = RepositoryPropertiesPopup(self)


  def closeRepository(self):

    pass
  
  def logout(self):

    self.setUsername(None)
    self.repList = RepositoryList()

    print 'in logout ', self.repList
    print '#########################'

    print self.tabbedFrame.frames[0].__dict__
    print '#########################'

    if self.frameShortcuts.has_key('Repository'):
      self.frameShortcuts['Repository'].drawFrame()

    #for ff in self.tabbedFrame.frames[0].children.values():
    #  ff.draw()
    for ff in self.tabbedFrame.frames[8].children.values():
      ff.draw()


  def openPopup(self, popup_name, clazz, oldStyle=False, *args, **kw):

    popup = self.popups.get(popup_name)
  
    if (popup):
      popup.open()
    else:
      if self.project:
        analysisProfile = self.analysisProfile
      else:
        analysisProfile = None
      
      if analysisProfile:
        transient = analysisProfile.transientDialogs
      else:
        transient = True
      
      name = popup_name
      
      if (oldStyle):
        popup = self.popups[popup_name] = clazz(self, transient=transient, *args, **kw)
      else:
        popup = self.popups[popup_name] = clazz(self, project=self.project, popup_name=name,
                                                transient=transient, *args, **kw)
      # above automatically opens popup

    return popup


  def closeProject(self, queryClose = True, querySave = True):

    if queryClose:
      if not showYesNo('Close project', 'Close current project?', parent=self):
        return False

    if querySave:
      if not self.checkSaving():
        return False

    self.destroyPopups()
    self.initProject()
 
    return True

  def saveProject(self):

    self.saveFile()
    
    # NBNB Rasmus 22/8/08 changed to fit AnalysisPopup (and removed 'isStored' bug)
    #if (self.project.isStored):
    #  self.saveFile()
    #else:
    #  self.askSaveFile()

  def saveAsProject(self):

    self.askSaveFile()

  def showVersion(self):

    showInfo('Version', self.versionInfo, parent=self)


  def changeFont(self, analysisProfile):
  
    if analysisProfile is self.analysisProfile:
      if self.font != analysisProfile.font:
        self.font = analysisProfile.font or DEFAULT_FONT 
        self.setFont()
      
  def selectFont(self, font): # Only called from main menu option
    
    if self.project:
      self.analysisProfile.font = font    

  def setFont(self, font=None, popup=None):

    if font is None:
      font = self.font or DEFAULT_FONT 
          
    else:
      self.font = font

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
              child.config( font=font)
          else:
            child.config( font=font )
          break
    
      if isinstance(child, Tkinter.Frame):
        child.font = font
      
      elif isinstance(child, ScrolledMatrix):
        child.setFont(font)  
    
      if hasattr(child, 'children'):
        childList.extend( child.children.values() )
        
  def initProject(self, project=None):

    self.project = project
    if project:
      self.initTopObjects(self.project)
      self.project.application = self.application
      self.setFont()
      webBrowser = ProjectWebBrowser(self.top, popup=self, project=project)
      #self.helpButton.config(command = lambda url=self.help_url: webBrowser.open(url) )
      
      for i, func in enumerate(self.updateFuncs):
        func(project)
        
      for button in self.projButtons:
        button.enable()
      
      self.setupSoftware()
        
    else:
      for button in self.projButtons:
        button.disable()
    
    self.setMenuState()


  def setUsername(self, name):

    self.username = name

    print 'set username ', self.username

  # jmci FIXME

  # need to look into this carefully and analyse the possible
  # states
    
  def setMenuState(self):

    if self.project:
      state = Tkinter.NORMAL
    else:
      state = Tkinter.DISABLED

    # Disable bits of the project menu if there's no project
    #for option in ('Save','Save As','Close'):
    #  i = self.projectMenu.options.index(option)
    #  self.projectMenu.entryconfig(i, state=state)

    # Disable other manus if there's no project
    #for menu in [self.otherMenu,]: # Include more menus in this list
    #  for i in range(len(menu.options)):
    #    menu.entryconfig(i, state=state)
      
  
  def setupSoftware(self):
    
     project = self.project
    
     methodStore = project.currentMethodStore or \
                   project.findFirstMethodStore() or \
                   project.newMethodStore(name=project.name)
        
     software = methodStore.findFirstSoftware(name=PROGRAM_NAME, version=VERSION)
     if not software:
       software = methodStore.newSoftware(name=PROGRAM_NAME, version=VERSION)

     #software.vendorName = 'PDBe & CCPN'
     #software.vendorAddress = ''
     #software.vendorWebAddress = 'http:'
     #software.details = ''
      
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

  def destroyPopups(self):

    for key in self.popups.keys():
      popup = self.popups[key]
      popup.destroy()

    self.popups = {}


  ## WMS methods


  # should maybe merge this with setUsername

  # this will read a config file that stores user setups. How this
  # is done TBD so hack for now and just return
  def readConfig(self, userName):

    # this should return the full RepositoryList structure. In reality this
    # will parse from file. set it by hand from here for now for speed
    repList = RepositoryList(userName)

    rep1 = repList.newRepository('remote','localhost:8080/HalxService','ccpn','remote')
    repList.currentRepository = rep1
    # should these be directory objects?
    repList.current_export_dir = '/home/jionides/work/CCPN/test_JMCI_data'
    repList.current_import_dir = '/home/jionides/work/CCPN/test_EXP_data'
   

    self.repList = repList

    print 'in readconfig '

    self.setRepository(rep1)

    # really should be done through notifiers
    if self.frameShortcuts.has_key('Repository'):
      self.frameShortcuts['Repository'].drawFrame()

    if self.frameShortcuts.has_key('Tasks'):
      self.frameShortcuts['Tasks'].drawFrame()

    if self.frameShortcuts.has_key('DataExch'):
      self.frameShortcuts['DataExch'].drawFrame()

    for ff in self.tabbedFrame.frames[8].children.values():
      ff.draw()

    # want to update the top menu. The idea should be to add a bar item
    # and then to list all the repositories underneath. Clicking on a given
    # repository will set the home repository as current. That seems to
    # be a good sensible way of simplifying the interface in a expected
    # way. We could, alternatively, use tabs.

    opts = {}
    self.repMenu.add('separator',**opts)
    
    for rep in repList.repositories:
      item = {'command': lambda rep=rep: self.setRepository(rep),
              'label': rep.name}
      self.repMenu.add('comman', **item)

    self.repMenu.config()
              

  ## end WMS methods


  def copyModifiedStorages(self):

    modifiedStorages = [storage for storage in self.project.storages if storage.isModified]
    if (self.project.isModified):
      modifiedStorages.append(self.project)

    for storage in modifiedStorages:
      self.copyStorage(storage)

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

  def initTopObjects(self, project):

    if not project.currentInstrumentStore:
      instrumentStores = project.sortedInstrumentStores()
      if instrumentStores:
        project.currentInstrumentStore = instrumentStores[0]
      else:
        # this automatically sets project.currentInstrumentStore
        project.newInstrumentStore(name=project.name)

    if not project.currentNmrProject:
      nmrProjects = project.sortedNmrProjects()
      if nmrProjects:
        project.currentNmrProject = nmrProjects[0]
      else:
        # this automatically sets project.currentNmrProject
        project.newNmrProject(name=project.name)

    if not project.currentAnalysisProject:
      analysisProjects = project.sortedAnalysisProjects()
      if analysisProjects:
        project.currentAnalysisProject = analysisProjects[0]
      else:
        # this automatically sets project.currentAnalysisProject
        project.newAnalysisProject(name=project.name, nmrProject=project.currentNmrProject)

    if not project.currentAnalysisProfile:
      analysisProfiles = project.sortedAnalysisProfiles()
      if analysisProfiles:
        project.currentAnalysisProfile = analysisProfiles[0]
      else:
        # this automatically sets project.currentAnalysisProfile
        project.newAnalysisProfile(name=project.name)

    self.nmrProject = project.currentNmrProject
    self.analysisProject = project.currentAnalysisProject
    self.analysisProfile = project.currentAnalysisProfile
    self.analysisProject.isThreadingAllowed = False


  def quit(self):

    msg = 'Quit ' + self.programName

    if not showYesNo('Confirm', msg + '?', parent=self):
      return

    # turn off the daemon!
    self.taskDaemon.alive=False

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


def getProjectFont(project):

  analysisProject = project.currentAnalysisProject

  if analysisProject:
    font = analysisProject.font
  else:
    font = 'Helvetica 10'
  
  return font

def setProjectFont(project, font):

  
  analysisProject = project.currentAnalysisProject
  nmrProject = project.currentNmrProject

  if not analysisProject:
    analysisProject = project.newAnalysisProject(name=project.name,
                                                 nmrProject=nmrProject)
    
  analysisProject.font = font

def launchApplication(projectDir=None):

  global top

  root = Tkinter.Tk()
  root.withdraw()  # get rid of the root window Tkinter would otherwise put up
  top  = ApplicationPopup(root) # creates our controlled window (this class in fact)
 
  project = None
  if projectDir:
    projectDir = normalisePath(projectDir)
    askdir  = lambda title, prompt, initial_value, default_dir: askDir(title, prompt,
              initial_value, parent=top, extra_dismiss_text='Skip')
    askfile = lambda title, prompt, initial_value: askFile(title, prompt,
              initial_value, parent=top, extra_dismiss_text='Skip')
    try:
      project = loadProject(path=projectDir, showWarning=showWarning, askDir=askdir, askFile=askfile)
    except ApiError, e:
      showError('Reading project', e.error_msg, parent=top)
 
  top.update_idletasks()
  
  if project:
    top.initProject(project)
 
  # waiting for mouse and keyboard commands
  root.mainloop()  

if __name__ == '__main__':

  argv = sys.argv[:]
  n    = len(argv)

  ## we are going to have to adjust this because again, there is not really
  ## a single project root for WMS.

  ## we are almost certainly going to need a set of local config files that
  ## store various bits and pieces. It might make sense for this to be
  ## modelled and to use a Memops-generated API
  
  if n > 1:
    projectDir = argv[1]
  else:
    projectDir = None

  launchApplication(projectDir)

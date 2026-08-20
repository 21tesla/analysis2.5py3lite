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
from memops.gui.ButtonList      import UtilityButtonList
from memops.gui.FontMenu        import FontMenu
from memops.gui.Frame           import Frame
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FileSelect      import FileType
from memops.gui.Label           import Label
from memops.gui.Menu            import Menu
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo
from memops.gui.ScrolledMatrix  import ScrolledMatrix

from memops.editor.BasePopup          import BasePopup
from memops.editor.OpenProjectPopup   import OpenProjectPopup
from memops.editor.SaveProjectPopup   import SaveProjectPopup
from memops.editor.WebBrowser         import ProjectWebBrowser

from ccpnmr.format.converters.NmrStarFormat import NmrStarFormat
from ccpnmr.format.converters.FastaFormat import FastaFormat
from ccpnmr.format.converters.CnsFormat import CnsFormat
from ccpnmr.format.converters.PseudoPdbFormat import PseudoPdbFormat
from ccpnmr.eci.ReadPdb import ReadPdb

DEFAULT_FONT = 'Helvetica 10'

PROGRAM_NAME = 'CcpNmr Entry Completion Interface'

VERSION = '2.1'

from ccpnmr.eci.EntryCompletionFrame import EntryCompletionFrame

class EntryCompletionGui(BasePopup):

  '''
The standalone version of CcpNmr Entry Completion Interface (ECI).  With it, you can complete all the necessary information required for PDB and BMRB data depositions by adding an "Entry" object to your CCPN project. An "Entry" object contains all the information that you wish to deposit with your submission. You can also select chemical shift lists, peak lists, structural restraints, ensembles, etc. and all at the click of a button. In addition, you can add all the meta data that is required for submissions to the PDB and BMRB. This can be done securely on your desktop computer over the duration of your NMR project.
  '''

  help_url = 'http://www.ccpn.ac.uk'

  def __init__(self, root):

    self.root = root

    self.font = DEFAULT_FONT
    # Application object needed to store application-specific data with project
    self.application = Application(name=PROGRAM_NAME)
    self.versionInfo = VERSION
    
    # Application popup is a superclass of memops.editor.BasePoup
    BasePopup.__init__(self, parent=root, title=PROGRAM_NAME,
                       location='+100+100', class_=self.application.name)

    #self.setTitle(PROGRAM_NAME)

                       
  def body(self, guiParent):

    self.geometry('900x750')

    # Ensure that the first row and column in popup expand
    guiParent.expandGrid(0,0)
    
    self.eciFrame = EntryCompletionFrame(guiParent, basePopup=self, grid=(0,0))
 
    # Add help, close buttons
    #buttonList = UtilityButtonList(self.eciFrame.tabbedFrame.sideFrame,
    #                               grid=(0,0), sticky='e')
    
    # Dictionary to store popups opened by this application - e.g. need to close upon quit
    self.popups = {}
    
    # Default font
    self.font   = 'Helvetica 10'

    # Closing the window from the desktop environment calls the proper quit
    self.protocol('WM_DELETE_WINDOW', self.quit)

    self.mainMenu    = Menu(self)
    self.projectMenu = self.makeProjectMenu()
    self.ccpNmrMenu  = self.makeCcpNmrMenu()
    self.otherMenu   = self.makeOtherMenu()
    
    # Put the main menu 
    self.config(menu=self.mainMenu)

    self.initProject()

  def makeProjectMenu(self):

    importsMenu = Menu(self.mainMenu, tearoff=False)
    importsMenu.add_command(label='Fasta Sequence',command=self.importFastaSeq)
    importsMenu.add_command(label='NMR-STAR 2.1.1', command=self.importNmrStar211)
    importsMenu.add_command(label='NMR-STAR 3.1', command=self.importNmrStar31)
    importsMenu.add_command(label='PDB 3.20',command=self.importPdb)
    importsMenu.add_command(label='PDB Coords',command=self.importPdbCoords)
    #importsMenu.add_command(label='CNS Distance Restraints',command=self.importCnsDistance)

    exportsMenu = Menu(self.mainMenu, tearoff=False)
    exportsMenu.add_command(label='NMR-STAR 3.1', command=self.eciFrame.exportNmrStar31)
    exportsMenu.add_command(label='PDB Coords',command=self.eciFrame.exportPdb)
    #exportsMenu.add_command(label='CNS Distance Restraints',command=self.eciFrame.exportCnsDistance)

    # Submenu of the min menu
    menu = Menu(self.mainMenu, tearoff=0)
    # Add various options to the menu and state the functions they call
    menu.add_command(label='New',     shortcut='N', command=self.newProject)
    menu.add_command(label='Open',    shortcut='O', command=self.openProject)
    menu.add_command(label='Close',   shortcut='C', command=self.closeProject)
    menu.add_command(label='Save',    shortcut='S', command=self.saveProject)
    menu.add_command(label='Save As', shortcut='A', command=self.saveAsProject)
    menu.add_cascade(label='Import',  shortcut='I', menu=importsMenu)
    menu.add_cascade(label='Export',  shortcut='E', menu=exportsMenu)
    menu.add_command(label='Quit',    shortcut='Q', command=self.quit)
    menu.add_command(label='Version', shortcut='V', command=self.showVersion)
    self.mainMenu.add_cascade(label='Project', shortcut='P', menu=menu)

    menu.options = ['New','Open','Close','Save','Save As','Import','Export','Quit','Version']  
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

  def openProject(self):

    if self.project:
      if not self.closeProject():
        return

    self.openPopup('open_project', OpenProjectPopup,
                   callback=self.initProject)

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

  def saveAsProject(self):

    self.askSaveFile()

  def importNmrStar211(self):

    if not self.project:
      showWarning('Failure','Please create a new CCPN project.', parent=self)
      return

    fileTypes = [ FileType('STAR', ['*.str']),
                  FileType('All', ['*'])]

    fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                        title='Import NMR-STAR 2.1.1 file', dismiss_text='Cancel',
                        selected_file_must_exist=True, multiSelect=False,)

    fileName = fileSelectPopup.getFile()

    if not os.path.exists(fileName):
      showWarning('Failure','No such file.', parent=self)
      return

    if fileName:
      nmrStarObj = NmrStarFormat(self.project, self, verbose=True)
      nmrStarObj.readProject(fileName, minimalPrompts=True, version='2.1.1')

      self.eciFrame.updateAll()

      return nmrStarObj

  def importNmrStar31(self):

    if not self.project:
      showWarning('Failure','Please create a new CCPN project.', parent=self)
      return

    fileTypes = [ FileType('STAR', ['*.str']),
                  FileType('All', ['*'])]

    fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                        title='Import NMR-STAR 3.1 file', dismiss_text='Cancel',
                        selected_file_must_exist=True, multiSelect=False,)

    fileName = fileSelectPopup.getFile()

    if not os.path.exists(fileName):
      showWarning('Failure','No such file.', parent=self)
      return

    if fileName:
      nmrStarObj = NmrStarFormat(self.project, self, verbose=True)
      nmrStarObj.readProject(fileName, minimalPrompts=True, version='3.1')

      self.eciFrame.updateAll()

      return nmrStarObj

  def importPdb(self):

    if not self.project:
      showWarning('Failure','Please create a new CCPN project.', parent=self)
      return

    fileTypes = [ FileType('PDB', ['*.pdb']),
                  FileType('PDB Entry', ['*.ent']),
                  FileType('All', ['*']) ]

    fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                        title='Import PDB 3.20 file', dismiss_text='Cancel',
                        selected_file_must_exist=True, multiSelect=False,)

    fileName = fileSelectPopup.getFile()

    if not os.path.exists(fileName):
      showWarning('Failure','No such file.', parent=self)
      return

    if fileName:
      name = os.path.split(fileName)[1]
      if '.' in name:
        name = name.split('.')[0]

      pdbObj = ReadPdb(fileName, self.project, 'PDB_'+ name, self)

      self.eciFrame.updateAll()

      return pdbObj

  def importCnsDistance(self):

    if not self.project:
      showWarning('Failure','Please create a new CCPN project.', parent=self)
      return

    fileTypes = [ FileType('CNS', ['*.tbl']),
                  FileType('All', ['*']) ]

    fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                        title='Import CNS distance restraints file', dismiss_text='Cancel',
                        selected_file_must_exist=True, multiSelect=False,)

    fileName = fileSelectPopup.getFile()

    if not os.path.exists(fileName):
      showWarning('Failure','No such file.', parent=self)
      return

    if fileName:
      name = os.path.split(fileName)[1]
      if '.' in name:
        name = name.split('.')[0]

      cnsObj = CnsFormat(self.project, self, verbose=True)
      cnsObj.readDistanceConstraints(fileName, minimalPrompts=True)

      self.eciFrame.updateAll()

      return cnsObj

  def importFastaSeq(self):

    if not self.project:
      showWarning('Failure','Please create a new CCPN project.', parent=self)
      return

    fileTypes = [ FileType('Fasta', ['*.fsa']),
                  FileType('All', ['*']) ]

    fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                        title='Import Fasta sequence file', dismiss_text='Cancel',
                        selected_file_must_exist=True, multiSelect=False,)

    fileName = fileSelectPopup.getFile()

    if not os.path.exists(fileName):
      showWarning('Failure','No such file.', parent=self)
      return

    if fileName:
      name = os.path.split(fileName)[1]
      if '.' in name:
        name = name.split('.')[0]

      fastaObj = FastaFormat(self.project, self, verbose=True)
      fastaObj.readSequence(fileName, minimalPrompts=True)

      self.eciFrame.updateAll()

      return fastaObj

  def importPdbCoords(self):

    if not self.project:
      showWarning('Failure','Please create a new CCPN project.', parent=self)
      return

    fileTypes = [ FileType('PDB', ['*.pdb']),
                  FileType('PDB Entry', ['*.ent']),
                  FileType('All', ['*']) ]

    fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                        title='Import PDB coordinates file', dismiss_text='Cancel',
                        selected_file_must_exist=True, multiSelect=False,)

    fileName = fileSelectPopup.getFile()

    if not os.path.exists(fileName):
      showWarning('Failure','No such file.', parent=self)
      return

    if fileName:
      name = os.path.split(fileName)[1]
      if '.' in name:
        name = name.split('.')[0]

      pdbObj = PseudoPdbFormat(self.project, self, verbose=True)
      pdbObj.readCoordinates(fileName, minimalPrompts=True)

      self.eciFrame.updateAll()

      return pdbObj

  def showVersion(self):

    showInfo('Version', 'Version ' + self.versionInfo, parent=self)

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
      
      self.setupSoftware()
      self.eciFrame.updateAll(project)
      
    
    self.setMenuState()
  
  def setupSoftware(self):
    
     project = self.project
    
     methodStore = project.currentMethodStore or \
                   project.findFirstMethodStore() or \
                   project.newMethodStore(name=project.name)
        
     software = methodStore.findFirstSoftware(name=PROGRAM_NAME, version=VERSION)
     if not software:
       software = methodStore.newSoftware(name=PROGRAM_NAME, version=VERSION)

     software.vendorName = 'PDBe & CCPN'
     #software.vendorAddress = ''
     software.vendorWebAddress = 'http://www.ebi.ac.uk/pdbe/'
     software.details = 'The PDBe & CCPN Entry Completion Interface'
     software.tasks = ['data deposition']

  def setMenuState(self):

    if self.project:
      state = Tkinter.NORMAL
    else:
      state = Tkinter.DISABLED

    # Disable bits of the project menu if there's no project
    for option in ('Close','Save','Save As',):
      i = self.projectMenu.options.index(option)
      self.projectMenu.entryconfig(i, state=state)

    # Disable other menus if there's no project
    for menu in [self.otherMenu,]: # Include more menus in this list
      for i in range(len(menu.options)):
        menu.entryconfig(i, state=state)
      
      
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


def launchApplication(projectDir=None):

  global top

  root = Tkinter.Tk()
  root.withdraw()
  top  = EntryCompletionGui(root)
 
  project = None
  if projectDir:
    projectDir = normalisePath(projectDir)
    askdir  = lambda title, prompt, initial_value, default_dir: askDir(title, prompt,
              initial_value, parent=top, extra_dismiss_text='Skip')
    askfile = lambda title, prompt, initial_value: askFile(title, prompt,
              initial_value, parent=top, extra_dismiss_text='Skip')
    try:
      project = loadProject(path=projectDir, showWarning=showWarning, askDir=askdir, askFile=askfile)

      userRepos = project.findFirstRepository(name='userData')
      chemCompPackLoc = project.findFirstPackageLocator(targetName='ccp.molecule.ChemComp')

      if userRepos not in chemCompPackLoc.repositories:
        chemCompPackLoc.addRepository(userRepos)

    except ApiError, e:
      showError('Reading project', e.error_msg, parent=top)
 
  top.update_idletasks()
  
  if project:
    top.initProject(project)
    
  root.mainloop()  

if __name__ == '__main__':

  argv = sys.argv[:]
  n    = len(argv)

  if n > 1:
    projectDir = argv[1]
  else:
    projectDir = None

  launchApplication(projectDir)

"""
============================COPYRIGHT START=============================

DangleGui.py: Part of the DANGLE package (release v1.1)

DANGLE: Dihedral ANgles from Global Likelihood Estimates 
Copyright (C) 2008 Nicole Cheung, Tim Stevens, Bill Broadhurst (University of Cambridge)

========================================================================


If you make use of the software or any documents in this package, 
please give credit by citing this package, its authors and references
in the literature with the same authors.

We would appreciate hearing of any problems you may encounter, 
but the programs, the documents and any files created by the programs 
are provided WITHOUT ANY WARRANTY and without even the implied warranty of
CORRECTNESS, MERCHANTABILITY or FITNESS FOR A PARTICULAR OR GENERAL USE.

THE RESPONSIBILITY FOR ANY ADVERSE CONSEQUENCES FROM THE USE OF PROGRAMS OR
DOCUMENTS OR ANY FILE OR FILES CREATED BY USE OF THE PROGRAMS OR DOCUMENTS
LIES SOLELY WITH THE USERS OF THE PROGRAMS OR DOCUMENTS OR FILE OR FILES AND
NOT WITH AUTHORS OF THE PROGRAMS OR DOCUMENTS.

You are not permitted to use any pieces of DANGLE in other programs, make
modifications to DANGLE, or make what a lawyer would call a "derived work" 
in any other way without the consent from either author.
 

============================COPYRIGHT END===============================

for further information, please contact the authors:

- Nicole Cheung   : msc51@cam.ac.uk

- Tim Stevens     : tjs23@cam.ac.uk

- Bill Broadhurst : r.w.broadhurst@bioc.cam.ac.uk

========================================================================

If you are using this software for academic purposes, we suggest
quoting the following reference:

===========================REFERENCE START==============================

DANGLE: To be completed.

CCPN: 
Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END================================
"""

import os
import sys
import Tkinter

from memops.api                    import Implementation
from memops.universal.Io           import normalisePath
from memops.general.Implementation import ApiError

from memops.general.Application    import Application
from memops.general.Io             import saveProject, loadProject
from memops.general.Util           import isProjectModified

from memops.gui.DataEntry       import askString, askDir, askFile
from memops.gui.FontMenu        import FontMenu
from memops.gui.Frame           import Frame
from memops.gui.Menu            import Menu
from memops.gui.ScrolledGraph   import ScrolledGraph
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo

#from memops.editor.BackupProjectPopup   import BackupProjectPopup
from memops.editor.OpenProjectPopup     import OpenProjectPopup
from memops.editor.SaveProjectPopup     import SaveProjectPopup

from ccpnmr.analysis.popups.BasePopup import BasePopup
from cambridge.dangle.DangleFrame import DangleFrame

PROGRAM_NAME = 'DANGLE'

VERSION = '1.1'

class DangleGui(BasePopup):
  """
  **A Bayesian inferential prediction method for protein backbone dihedral angles**
  
  DANGLE (Dihedral ANgles from Global Likelihood Estimates) predicts protein
  backbone Phi and Psi angles and secondary structure assignments solely from
  amino acid sequence information, experimental chemical shifts and a database
  of known protein structures and their associated shifts. This approach
  uses Bayesian inferential logic to analyse the likelihood of conformations
  throughout Ramachandran space, paying explicit attention to the population
  distributions expected for different amino acid residue types.

  Simple filtering procedures can identify the most "predictable" residues,
  yielding 92% of all Phi and Psi predictions accurate to within +/-30 degrees.
  In contrast to previous approaches, more than 80% of Phi or Psi predictions
  for glycine and pre-proline are reliable. Furthermore, DANGLE provides
  meaningful upper and lower bounds for the predictions which are shown to
  represent the precision of the prediction. Over 90% of the experimental
  dihedral angles in the set of test proteins are within the boundary ranges
  suggested by DANGLE. At a lower resolution level, the program correctly
  assigns each residue to one of three secondary structure states (H, E or C) in
  85% of cases.

  DANGLE also provides an indication of the degeneracy in the relationship
  between shift measurements and conformation at each site.

  For more documentation see the `DANGLE web site`_ at SourceForge.
  
  .. _`DANGLE web site`: http://dangle.sourceforge.net/

  **Reference**
  
  *Ming-Sin Cheung, Mahon L. Maguire, Tim J. Stevens, and R. William Broadhurst.
  "DANGLE: A Bayesian inferential method for predicting protein backbone
  dihedral angles and secondary structure." Journal of Magnetic Resonance Volume
  202, Issue 2, February 2010, Pages 223-233*

  """

  def __init__(self, root, programName=PROGRAM_NAME, **kw):

    self.parent = root
    self.programName = programName
    # Application object needed to store application-specific data with project
    self.application = Application(name=PROGRAM_NAME) 
    self.versionInfo = 'Version ' + VERSION
    self.graphPopup  = None
    self.analysisProfile = None
    
    # Application popup is a superclass of memops.editor.BasePoup
    BasePopup.__init__(self, parent=root, title=PROGRAM_NAME,
                       location='+100+100', class_=self.application.name, **kw)
                       
    
  def body(self, guiParent):

    # Ensure that the first row and column in popup expand
    guiParent.grid_rowconfigure(1, weight=1)
    guiParent.grid_columnconfigure(0, weight=1, minsize=200)
    frame = Frame(guiParent) # Body widgets can be put in this frame
    frame.grid(row=1, column=0, sticky='nsew')
    
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1, minsize=200)

    self.dangleFrame = DangleFrame(frame, self, project=self.project)
    self.dangleFrame.grid(row=0, column=0, sticky='nsew')

    # Dictionary to store popups opened by this application - e.g. need to close upon quit
    self.popups = {}
    
    # Default font
    self.font   = 'Helvetica 10'

    # Closing the window from the desktop environment calls the proper quit
    self.protocol('WM_DELETE_WINDOW', self.quit)

    self.mainMenu    = Menu(self)
    self.projectMenu = self.makeProjectMenu()
    self.viewMenu    = self.viewFileMenu()
    self.otherMenu   = self.makeAppearanceMenu()
    
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
    menu.add_command(label='Backup',  shortcut='B', command=self.backupProject,
                     tipText='Manage options for automatic CCPN project backup')
    menu.add_command(label='Quit',    shortcut='Q', command=self.quit,
                     tipText='Exit from the DANGLE program')
    menu.add_command(label='Version', shortcut='V', command=self.showVersion)
    self.mainMenu.add_cascade(label='Project', shortcut='P', menu=menu)

    menu.options = ['New','Open','Close','Save','Save As','Backup','Quit','Version']  
    return menu


  def viewFileMenu(self):

    # Submenu of the min menu
    menu = Menu(self.mainMenu, tearoff=0)
    # Add various options to the menu and state the functions they call
    menu.add_command(label='Prediction Graph', shortcut='G', command=self.viewGraph,
                     tipText='Show a graphical view of the poredicted backbone dihedral angles')
    self.mainMenu.add_cascade(label='View', shortcut='V', menu=menu)

    menu.options = ['Prediction Graph']  
    return menu


  def makeAppearanceMenu(self):
    
    # The fonts menu is a pre-created widget
    fontsMenu = FontMenu(self.mainMenu, self.setFont, sizes=(8,10,12),
                         doItalic=0, doBoldItalic=0, tearoff=0)

    # Plot color schemes
    colorMenu = Menu(self.mainMenu, tearoff=0)
    colorMenu.add_radiobutton(label='Red', command=self.colorSchemeRed,)
    colorMenu.add_radiobutton(label='Rainbow', command=self.colorSchemeRainbow,)
    #colorMenu.add_radiobutton(label='Black',   shortcut='B', command=self.colorSchemeBlack)
    #colorMenu.add_radiobutton(label='White',   shortcut='W', command=self.colorSchemeWhite)

    # Submenu of the main menu
    menu = Menu(self.mainMenu, tearoff=0)
    # Only Fonts option so far
    menu.add_cascade(label='Plot Colour Scheme',shortcut='C', menu=colorMenu)
    menu.add_cascade(label='Fonts',             shortcut='F', menu=fontsMenu)
    self.mainMenu.add_cascade(label='Appearance', shortcut='A', menu=menu)

    menu.options = ['Plot Colour Scheme','Fonts']
    return menu

  def colorSchemeRed(self):
  
    self.dangleFrame.colorScheme = 'red'
    
  def colorSchemeRainbow(self):
  
    self.dangleFrame.colorScheme = 'rainbow'

  def colorSchemeBlack(self):
  
    self.dangleFrame.colorScheme = 'black'

  def colorSchemeWhite(self):
  
    self.dangleFrame.colorScheme = 'white'
  
  def viewGraph(self):
    
    phiData, psiData = self.dangleFrame.getPhiPsiPredictions()
    
    if self.graphPopup:
      self.graphPopup.open()
    else:
      self.graphPopup = DangleGraphPopup(self)
    
    dangleChain = self.dangleFrame.dangleChain
    
    good = []
    med = []
    bad = []
    if dangleChain:
      
      resNum = 0 # In case no residues
      
      for dResidue in dangleChain.dangleResidues:
        resNum = dResidue.residue.seqCode
        nIslands = dResidue.numIslands
        
        if nIslands < 2:
          good.append((resNum, nIslands))
        elif nIslands < 4:
          med.append((resNum, nIslands))
        else:
          bad.append((resNum, nIslands))
      
      good.append((resNum, 0)) # To bias the graph
      
    islandData = [good, med, bad]
    
    self.graphPopup.update(phiData,psiData,islandData)
    

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
      nmrProject = project.newNmrProject(name = project.name)
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
      if self.analysisProfile:
        transient = self.analysisProfile.transientWindows
      else:
        transient = True
        
      # Below automatically opens popup
      popup = self.popups[popup_name] = clazz(self, project=self.project, transient=transient, *args, **kw)

    return popup


  def closeProject(self, queryClose=True, querySave=True):

    if queryClose:
      if (not showYesNo('Close project', 'Close current project?')):
        return False

    if querySave:
      if (not self.checkSaving()):
        return False

    self.destroyPopups()
    self.initProject()
 
    return True

  def saveProject(self):

    #if (self.project.isStored):
    self.saveFile()
    #else:
    #  self.askSaveFile()

  def saveAsProject(self):

    self.askSaveFile()

  def backupProject(self):

    pass
    #self.openPopup('backup_project', BackupProjectPopup)

  def showVersion(self):

    showInfo('Version', self.versionInfo, parent=self)


  def setFont(self, fontName=None, popup=None):

    if fontName is None:
      if self.analysisProfile:
        # Get font specification saved as project application data
        self.font = self.analysisProfile.font
    
    else:
      self.font = fontName
      if self.analysisProfile:
        # Set font specification as project application data
        self.analysisProfile.font = fontName
    
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

    if project:
      self.project = project
    else:
      project = self.project
    
    self.setTitle(self.programName)

    if project:
      self.project.application = self.application
      self.setFont()
      self.dangleFrame.updateProject(project)
      self.setupSoftware()

      if not project.currentAnalysisProfile:
        analysisProfiles = project.sortedAnalysisProfiles()
        if analysisProfiles:
          project.currentAnalysisProfile = analysisProfiles[0]
        else:
          project.newAnalysisProfile(name=project.name)
 
      self.analysisProfile = project.currentAnalysisProfile
    
    self.setMenuState()
  
  def setupSoftware(self):
    
     project = self.project
    
     methodStore = project.currentMethodStore or \
                   project.findFirstMethodStore() or \
                   project.newMethodStore(name=project.name)
        
     software = methodStore.findFirstSoftware(name=PROGRAM_NAME, version=VERSION)
     if not software:
       software = methodStore.newSoftware(name=PROGRAM_NAME, version=VERSION)

     #software.vendorName = 'Nicole Cheung'
     #software.vendorAddress = ''
     #software.vendorWebAddress = 'http:'
     #software.details = ''
      
  def setMenuState(self):

    if self.project:
      state = Tkinter.NORMAL
    else:
      state = Tkinter.DISABLED

    # Disable bits of the project menu if there's no project
    for option in ('Save','Save As','Backup','Close'):
      i = self.projectMenu.options.index(option)
      self.projectMenu.entryconfig(i, state=state)

    # Disable bits of the project menu if not stand-alone
    if hasattr(self.parent, 'project'): 
      for option in ('New','Open','Close'):
        i = self.projectMenu.options.index(option)
        self.projectMenu.entryconfig(i, state=Tkinter.DISABLED)

    # Disable other manus and view menu if there's no project
    for menu in [self.otherMenu,self.viewMenu]: # Include more menus in this list
      for i in range(len(menu.options)):
        menu.entryconfig(i, state=state)
      
      
  def askSaveFile(self):

    popup = self.openPopup('save_project', SaveProjectPopup)

  def modalAskSaveFile(self):

    popup = SaveProjectPopup(self, project=self.project, dismiss_text='Cancel Quit',
                             modal=True)
    did_save = popup.did_save
    popup.destroy()

    return did_save

  def quitSaveProject(self):

    if self.project.activeRepositories:
      return self.saveFile()
    else:
      return self.modalAskSaveFile()

  def destroyPopups(self):

    for key in self.popups.keys():
      popup = self.popups[key]
      popup.destroy()

    self.popups = {}

  def saveFile(self):

    if (not self.project):
      return False

    try:
      self.project.saveModified()
      print 'Successfully saved project'
      return True
      
    except IOError, e:
      showError('Saving file', str(e))
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

    if not hasattr(self.parent, 'project'): 
      if not showYesNo('Quit ' + self.programName, 'Quit ' + self.programName + '?'):
        return

    #if self.project:
    #  if not self.checkSaving():
    #    return

    self.destroy()

    # Only exit system of not inside another app 
    if not hasattr(self.parent, 'project'): 
      # Need the below code to fix the problem with Bash
      # where the command line was getting screwed up on exit.
      if os.name == 'posix':
        os.system('stty sane')
 
      # sys.exit(0)
      os._exit(0)

  def destroy(self):
  
    BasePopup.destroy(self)


def launchDangle(filename=None):

  global top

  root = Tkinter.Tk()
  root.withdraw()
  top  = DangleGui(root)
 
  project = None
  if filename:
    file    = normalisePath(filename)
    askdir  = lambda title, prompt, initial_value, default_dir: askDir(title, prompt,
              initial_value, parent=top, extra_dismiss_text='Skip', default_dir=default_dir)
    askfile = lambda title, prompt, initial_value: askFile(title, prompt,
              initial_value, parent=top, extra_dismiss_text='Skip')
    try:
      applicationName = top.application.name
      project = loadProject(file, showWarning=showWarning, askDir=askdir,
                            askFile=askfile)
    except ApiError, e:
      showError('Reading project', e.error_msg)
 
  top.update_idletasks()
  
  if project:
    top.initProject(project)
    
  root.mainloop()  
  
  
  
class DangleGraphPopup(BasePopup):
  
  def __init__(self, parent):

    BasePopup.__init__(self, parent=parent, title=" Dangle Prediction Graph",
                       transient=False, borderwidth=6)

    parent.protocol("WM_DELETE_WINDOW", self.close)

    self.geometry('700x700+50+50')

  def body(self, guiFrame):
 
    guiFrame.expandGrid(0,0)
    guiFrame.expandGrid(1,0)
 
    self.anglesGraph = ScrolledGraph(guiFrame, xLabel='Residue',yLabel='Angle',
                                     relief='flat', symbolSize=1, grid=(0,0),
                                     zoom=1.0,width=500,height=200,
                                     dataColors=['#D00000','#0000D0'],
                                     dataNames=['Phi','Psi'])
    self.anglesGraph.draw()

    self.islandsGraph = ScrolledGraph(guiFrame, xLabel='Residue',yLabel='No. Islands',
                                     relief='flat', graphType='histogram',
                                     zoom=1.0,width=500,height=200 ,grid=(1,0),
                                     dataColors=['#8080FF','#FFFF00','#E00000'],
                                     dataNames=['1','2-3','4+'],)
    self.islandsGraph.draw()

  def update(self, phiData, psiData, islandData):
    
    self.anglesGraph.update(dataSets=[phiData,psiData])   
    self.islandsGraph.update(dataSets=islandData)
  
  
  

if (__name__ == '__main__'):

  argv = sys.argv[:]
  n    = len(argv)

  if n > 1:
    filename = argv[1]
  else:
    filename = None

  launchDangle(filename)

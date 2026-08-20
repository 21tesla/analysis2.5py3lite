
"""
======================COPYRIGHT/LICENSE START==========================

ExtendNmrGui.py: <write function here>

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

===========================REFERENCE END===============================
"""
import os
import sys
import Tkinter

from memops.api                    import Implementation

from memops.universal.Io           import normalisePath, getTopDirectory

from memops.general.Implementation import ApiError
from memops.general.Application    import Application
from memops.general.Io             import saveProject
from memops.general.Util           import isProjectModified

from memops.gui.Button          import Button
from memops.gui.ButtonList      import ButtonList
from memops.gui.Frame           import Frame
from memops.gui.Label           import Label
from memops.gui.Menu            import Menu
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo
from memops.gui.TabbedFrame     import TabbedFrame
from memops.gui.FontMenu        import FontMenu

from memops.editor.BasePopup          import BasePopup

from ccp.gui.Io import loadProject

# ARIA
from paris.aria.AriaExtendNmrFrame import AriaFrame

# AUREMOL
from regensburg.auremol.AuremolFrame import AuremolFrame

# CING
from nijmegen.cing.CingFrame import CingFrame

# HADDOCK
try:
  # Extend-NMR version
  from haddock.HaddockFrame import HaddockFrame
except ImportError:
  from utrecht.haddock.HaddockFrame import HaddockFrame

# ISD
from cambridge.isd.IsdFrame import IsdFrame

# MDD



# PRODEOMP
from gothenburg.prodecomp.ProdecompFrame import ProdecompFrame

# ECI
from ccpnmr.eci.EntryCompletionFrame import EntryCompletionFrame


DEFAULT_FONT = 'Helvetica 10'

PROGRAM_NAME = 'Extend-NMR GUI'

VERSION = '0.6'

from ccpnmr.analysis.AnalysisPopup import AnalysisPopup, window_popup_prefix
from ccpnmr.analysis.Analysis import Analysis

class ApplicationPopup(AnalysisPopup):

  help_url = 'http://www.ccpn.ac.uk'

  allowedCallbackFuncnames = ('init', 'save', 'close')

  def __init__(self, root, programName='Extend-NMR Interface'):

    self.font = DEFAULT_FONT
    # Application object needed to store application-specific data with project
    self.application = Application(name=PROGRAM_NAME)
    self.versionInfo = 'Version' + VERSION
    self.ariaProjectFile = None
    self.doneAnalysisInfo = False
    
    self.updateFuncs = []
    self.projButtons = []
    
    self.ariaPaths = []
    self.isdPaths = []

    AnalysisPopup.__init__(self, root)

    self.program_name = PROGRAM_NAME
    self.setTitle(PROGRAM_NAME)

  def printAnalysisCommandLineInfo(self, event):
    
    if not self.doneAnalysisInfo:
      Analysis.printCommandLineInfo(self)
      self.doneAnalysisInfo = True

  def printCommandLineInfo(self):
    
    print """
 For program documentation see:
 http://www.extend-nmr.eu   
    """

  def body(self, guiParent):
    
    self.menus = {}
    self.menu_items = {}
    
    self.fixedActiveMenus = {}

    self.popups = {}

    self.callbacksDict = {}

    self.selected_objects = []

    self.menubar = menubar = Menu(guiParent)

    self.font = DEFAULT_FONT
    
    self.setProjectMenu()
    # 

    menu = Menu(self.menubar, tearoff=0)
    menu.bind('<Button>', self.printAnalysisCommandLineInfo)
    self.menubar.add_cascade(label='CcpNmr Analysis', shortcut='C', menu=menu)
    self.menubar.add_command(label='FormatConverter', shortcut='F',
                             command=self.runFormatConverter)
                     
    self.menubar = menu
    self.initProject()
    self.setPeaksMenu()
    self.setMoleculeMenu()
    self.setAssignMenu()
    self.setResonanceMenu()
    self.setDataMenu()
    self.setStructureMenu()
    self.setChartMenu()
    self.setMacroMenu()
    self.setOtherMenu()
    self.setMenuState() # need to do it again because of OtherMenu state

    
    # Help Submenu
    
    helpMenu = Menu(self.menubar, tearoff=0)
    helpMenu.add_command(label='Version', shortcut='V',
                         command=self.showVersion)
    helpMenu.add_command(label='About',   shortcut='A',
                         command=self.showAbout)
    helpMenu.add_command(label='Help',    shortcut='H',
                         command=self.showHelp)

    menu.add_separator()
    menu.add_command(label='CCPN Updates', shortcut='U',    
                     image=self.iconRefresh, compound='left',
                     command=self.updateAnalysis,
                     tipText='Get any new patches and updates to CcpNmr')
    menu.add_cascade(label='CCPN Help', shortcut='H',   
                     image=self.iconHelp, compound='left',
                     menu=helpMenu)

    self.config(menu=menubar)

    # Ensure that the first row and column in popup expand
    guiParent.grid_rowconfigure(0, weight=1)
    guiParent.grid_columnconfigure(0, weight=1, minsize=200)
    frame = Frame(guiParent) # Body widgets can be put in this frame
    frame.grid()

    softwareOpts = ['Extend-NMR','ARIA 2','Auremol',
                    'CING',' ECI ','HADDOCK',' ISD ',
                    'PRODECOMP']

    self.tabbedFrame = TabbedFrame(guiParent, options=softwareOpts,
                                   toggleOff=False, selected=0,
                                   callback=self.toggleTab)
    self.tabbedFrame.grid(row=0, column=0, sticky='nsew')

    frames = self.tabbedFrame.frames 

    # Logos
    ccpnDir = getTopDirectory()
    
    imageDir = os.path.join(ccpnDir,'python','extendNmr','images')
    
    imageFile = os.path.join(imageDir,'Fp6Logo.gif')
    self.fp6Logo = Tkinter.PhotoImage(file=imageFile)
    imageFile = os.path.join(imageDir,'CingLogo.gif')
    self.cingLogo = Tkinter.PhotoImage(file=imageFile)
    imageFile = os.path.join(imageDir,'AriaLogo.gif')
    self.ariaLogo = Tkinter.PhotoImage(file=imageFile)
    imageFile = os.path.join(imageDir,'IsdLogo.gif')
    self.isdLogo = Tkinter.PhotoImage(file=imageFile)
    imageFile = os.path.join(imageDir,'HaddockLogo.gif')
    self.haddockLogo = Tkinter.PhotoImage(file=imageFile)
    imageFile = os.path.join(imageDir,'AuremolLogo.gif')
    self.auremolLogo = Tkinter.PhotoImage(file=imageFile)
    imageFile = os.path.join(imageDir,'CcpnLogo.gif')
    self.ccpnLogo = Tkinter.PhotoImage(file=imageFile)
    imageFile = os.path.join(imageDir,'ProdecompLogo.gif')
    self.prodecompLogo = Tkinter.PhotoImage(file=imageFile)
    imageFile = os.path.join(imageDir,'MddLogo.gif')
    self.mddLogo = Tkinter.PhotoImage(file=imageFile)
    imageFile = os.path.join(imageDir,'BrukerLogo.gif')
    self.brukerLogo = Tkinter.PhotoImage(file=imageFile)
    imageFile = os.path.join(imageDir,'MsdLogo.gif')
    self.msdLogo = Tkinter.PhotoImage(file=imageFile)

    self.initExtendNmr(frames[0])

    self.initAria(frames[1])

    self.initAuremol(frames[2])

    self.initCing(frames[3])

    self.initEci(frames[4])

    self.initHaddock(frames[5])
    
    self.initIsd(frames[6])
    
    self.initProdecomp(frames[7])

    self.initProject(self.project)

    if not self.project:
      for button in self.projButtons:
        button.disable()

    self.geometry('680x670')
    
  def setProjectMenu(self):
    ProjectMenu = 'Project'

    # Imports submenu    
    
    importsMenu = Menu(self.menubar, tearoff=False)
    importsMenu.add_command(label='Via Format Converter', shortcut='F', command=self.runFormatConverter)
    importsMenu.add_command(label='NMR-STAR 2.1.1', command=self.importNmrStar211 )
    importsMenu.add_command(label='NMR-STAR 3.1', shortcut='N', command=self.importNmrStar31 )
    importsMenu.add_command(label='PDB 3.20', shortcut='P',command=self.importPdb )
    importsMenu.add_command(label='Coordinates (PDB-style)', shortcut='C',command=self.importCoordinates )


    # Preferences submenu

    fontsMenu = FontMenu(self.menubar, self.selectFont, sizes=(8,10,12),
                         doItalic=False, doBoldItalic=False, tearoff=0)
    
    prefsMenu = Menu(self.menubar, tearoff=False)
    prefsMenu.add_cascade(label='Fonts', shortcut='F', 
                          image=self.iconFont, compound='left',
                          menu=fontsMenu,
                          tipText='Select font to use in the graphical interface')
    prefsMenu.add_command(label='Colour Schemes', 
                          image=self.iconTable, compound='left',
                          shortcut='C',
                          tipText='Edit and create colour schemes',
                          command=self.editColorSchemes)
    prefsMenu.add_command(label='Residue Codes',  
                          image=self.iconTable, compound='left',
                          shortcut='R', tipText='User-specified codes that override the standard residue names',
                          command=self.editResidueCodes)
    prefsMenu.add_command(label='User Options',  
                          image=self.iconTable, compound='left',
                          shortcut='U', tipText='General options for Analysis behaviour',
                          command=self.editProfiles)
    
   
    # 

    menu = Menu(self.menubar, tearoff=0)
    menu.add_command(label='New', shortcut='N',   
                     image=self.iconNewWindow, compound='left',
                     command=self.newProject,
                     tipText='Create a new, blank CCPN project (closes any open project)')
    menu.add_command(label='Open Project', shortcut='O',   
                     image=self.iconOpen, compound='left',
                     command=self.openProject,
                     tipText='Open a new CCPN project by selecting a project directory on disk')
    menu.add_command(label='Open Spectra', shortcut='p',   
                     image=self.iconOpenFile, compound='left',
                     command=self.openSpectrum,
                     tipText='Open spectrum data fom disk, creating a default CCPN project if needed')
    menu.add_command(label='Save', shortcut='S',   
                     image=self.iconSave, compound='left',
                     command=self.saveProject,
                     tipText='Save the current CCPN project on disk')
    menu.add_command(label='Save As', shortcut='A',   
                     image=self.iconSaveAs, compound='left',
                     command=self.saveAsProject,
                     tipText='Save the current CCPN project under a different name (project directory)')
    menu.add_cascade(label='Import', shortcut='I',   
                     image=self.iconImport, compound='left',
                     menu=importsMenu)
    menu.add_command(label='Close', shortcut='C',   
                     image=self.iconClose, compound='left',
                     command=self.closeProject,
                     tipText='Close the current CCPN project')
    menu.add_command(label='Quit', shortcut='Q',   
                     image=self.iconQuit, compound='left',
                     command=self.quit,
                     tipText='Quit Extend-NMR, closing any open CCPN project')
    menu.add_separator()
    menu.add_cascade(label='Preferences', shortcut='P',  
                     image=self.iconPrefs, compound='left',
                     menu=prefsMenu)
    menu.add_command(label='Validate',shortcut='V',    
                     image=self.iconTool, compound='left',
                     command=self.validateProject,
                     tipText='Check the current CCPN project for data model consistency errors')
    menu.add_command(label='Backup', shortcut='B',   
                     image=self.iconTool, compound='left',
                     command=self.backupProject,
                     tipText='Setup options for automated backup of CCPN project data')
    menu.add_command(label='Archive', shortcut='r',   
                     image=self.iconTool, compound='left',
                     command=self.archiveProject,
                     tipText='Save the current CCPN project in an archived form, e.g. tar gzipped')
    
    self.menubar.add_cascade(label=ProjectMenu, shortcut='j', menu=menu)
    self.menus[ProjectMenu] = menu
    self.menu_items[ProjectMenu] = ['New', 'Open Project', 'Open Spectra',
                                    'Save', 'Save As', 'Import', 'Close', 
				    'Quit', 'Preferences',
                                    'Validate',  'Backup', 'Archive',]
    
    # Menus that area ctive in absence of a project
    #for ii in (0,1,2,7,13,15):
    for ii in (0,1,2,7,):
      self.fixedActiveMenus[(ProjectMenu,ii)] = True

  def setPopupGeometries(self):

    for key in self.popups.keys():
      popup = self.popups[key]
      if not key.startswith(window_popup_prefix): # bit of a hack
        not self.setPopupGeometry(popup, key)

  def toggleTab(self, index):
    
    if (index > 0) and not self.project:
      showWarning('Warning', 'No active project', parent=self)
      self.tabbedFrame.select(0)
      return
     
    frame = self.tabbedFrame.frames[index]
    if hasattr(frame,'printOutDocString'):
      print frame.printOutDocString
      # only print it once
      del frame.printOutDocString

  def initCing(self, frame):

    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
    
    canvas = Tkinter.Canvas(frame, width=702, height=77)
    canvas.grid(row=0,column=0,sticky='ew')
    canvas.create_image(0,0, anchor='nw', image=self.cingLogo)
    
    cingFrame = CingFrame(frame, self)
    cingFrame.grid(row=1, column=0, sticky='nsew')

    #self.projButtons.extend(cingFrame.buttonBar.buttons)
    self.updateFuncs.append(cingFrame.updateAll)

  def initProdecomp(self, frame):
    
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
    frame.configure(bg='#FFFFFF')

    canvas = Tkinter.Canvas(frame, width=1024, height=90, bg='#FFFFFF')
    canvas.grid(row=0,column=0,sticky='ew')
    canvas.create_image(0,0, anchor='nw', image=self.prodecompLogo)
    
    prodecompFrame = ProdecompFrame(frame, basePopup=self, ccpnProject=self.project)
    prodecompFrame.grid(row=1, column=0, sticky='nsew')
    
    # set printOutDocString:
    frame.printOutDocString = prodecompFrame.printOutDocString

    #self.projButtons.extend(isdFrame.buttons)
    self.updateFuncs.append(prodecompFrame.updateAll)

  def initAuremol(self, frame):
  
    frame.expandGrid(1,0)
    
    refText = """Gronwald W, Brunner K, Kirchhofer R, Nasser A, Trenner J, Ganslmeier B,
Riepl H, Ried A, Scheiber J, Elsner R, Neidig K-P, Kalbitzer HR

AUREMOL, a New Program for the Automated Structure Elucidation of
Biological Macromolecules. Bruker Reports 2004; 154/155: 11-14
    """
    
    canvas = Tkinter.Canvas(frame, width=700, height=160, bg='#FFFFFF') 
    canvas.grid(row=0,column=0,sticky='ew')
    canvas.create_image(12,12, anchor='nw', image=self.auremolLogo)
    canvas.create_text(200, 10, anchor='nw', text=refText)
     
    auremolFrame = AuremolFrame(frame, self.project, grid=(1,0))

    #self.projButtons.extend(isdFrame.buttons)  
    self.updateFuncs.append(auremolFrame.updateAll)

  def initEci(self, frame):
    
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
    frame.parent = self # For notifiers
    frame.configure(bg='#FFFFFF')

    canvas = Tkinter.Canvas(frame, width=800, height=73, bd=3, bg='#FFFFFF') 
    canvas.grid(row=0,column=0,sticky='ew')
    canvas.create_image(10,10, anchor='nw', image=self.msdLogo)
    
    eciFrame = EntryCompletionFrame(frame, basePopup=self)
    eciFrame.grid(row=1, column=0, sticky='nsew')
    
    #self.projButtons.extend(isdFrame.buttons)
    self.updateFuncs.append(eciFrame.updateAll)

  def initHaddock(self, frame):
    
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
    frame.parent = self # For notifiers

    canvas = Tkinter.Canvas(frame, width=753, height=92) 
    canvas.grid(row=0,column=0,sticky='ew')
    canvas.create_image(0,0, anchor='nw', image=self.haddockLogo)
    
    haddockFrame = HaddockFrame(frame, self.project)
    haddockFrame.grid(row=1, column=0, sticky='nsew')

    #self.projButtons.extend(isdFrame.buttons)
    self.updateFuncs.append(haddockFrame.updateAll)

  def initIsd(self, frame):

    global isd
    
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)

    frame.configure(bg='#FFFFFF')
    
    canvas = Tkinter.Canvas(frame, width=800, height=91, bd=3, bg='#FFFFFF')
    canvas.grid(row=0,column=0,sticky='ew', padx=0)
    canvas.create_image(0,0, anchor='nw', image=self.isdLogo)
    
    isdFrame = IsdFrame(frame, self.project)
    isdFrame.grid(row=1, column=0, sticky='nsew')

    self.isd = isdFrame
    isd = self.isd

    #self.projButtons.extend(isdFrame.buttons)
    self.updateFuncs.append(isdFrame.updateAll)

  def initExtendNmr(self, frame):
    
    row = 0
    frame.config(bd=5)
    canvas = Tkinter.Canvas(frame, width=640, height=600,
                            bg='#FFFFFF', bd=5)
    canvas.grid(row=row,column=0,sticky='nsew')
    canvas.create_image(15,15,    anchor='nw', image=self.fp6Logo)
    canvas.create_image(330, 15, anchor='nw', image=self.prodecompLogo)
    canvas.create_image(330,115, anchor='nw', image=self.mddLogo)
    canvas.create_image(20, 200, anchor='nw', image=self.ccpnLogo)
    canvas.create_image(230,220, anchor='nw', image=self.msdLogo)
    canvas.create_image(450,200, anchor='nw', image=self.brukerLogo)
    canvas.create_image(5,  310, anchor='nw', image=self.isdLogo)
    canvas.create_image(230,300, anchor='nw', image=self.ariaLogo)
    canvas.create_image(450,290, anchor='nw', image=self.auremolLogo)
    canvas.create_image(5,  420, anchor='nw', image=self.haddockLogo)
    canvas.create_image(5,  520, anchor='nw', image=self.cingLogo)
    
    #row += 1
    #l1 = Label(frame, text='Welcome to the Extend NMR suite')
    #l1.grid(row=row, column=0, sticky='ew')

  def initAria(self, frame):

    welcomeMessage = \
    """ARIA Version 2.3. Authors: Benjamin Bardiaux, Michael Habeck,
Jens Linge, Therese Malliavin, Sean O'Donoghue, Wolfgang Rieping,
and Michael Nilges.

Rieping W., Habeck M., Bardiaux B., Bernard A., Malliavin T.E.,
Nilges M.(2007) ARIA2: automated NOE assignment and data
integration in NMR structure calculation. Bioinformatics 23:381-382"""
    
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
    
    canvas = Tkinter.Canvas(frame, width=700, height=114, bg='#FFFFFF')
    canvas.grid(row=0,column=0,sticky='ew')
    canvas.create_image(0,0, anchor='nw', image=self.ariaLogo)
    canvas.create_text(250, 10, anchor='nw', text=welcomeMessage)
    
    ariaFrame = AriaFrame(frame, self)
    ariaFrame.grid(row=1, column=0, sticky='nsew')

    self.projButtons.extend(ariaFrame.buttons)
    self.updateFuncs.append(ariaFrame.updateAll)

 
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
        if popup_name.startswith(window_popup_prefix): # bit of a hack
          transient = analysisProfile.transientWindows
        else:
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

  def initProject(self, project=None):

    AnalysisPopup.initProject(self, project, False)

    self.project = project
    if project:
      for i, func in enumerate(self.updateFuncs):
        func(project)

      for button in self.projButtons:
        button.enable()

    else:
      for button in self.projButtons:
        button.disable()

  def setupSoftware(self):
     
     return
    
     project = self.project
    
     methodStore = project.currentMethodStore or \
                   project.findFirstMethodStore() or \
                   project.newMethodStore(name=project.name)
        
     software = methodStore.findFirstSoftware(name=PROGRAM_NAME, version=VERSION)
     if not software:
       software = methodStore.newSoftware(name=PROGRAM_NAME,
                                          version=VERSION)

     software.task = 'data pipeline'
     software.vendorName = 'The Extend-NMR Project'
     #software.vendorAddress = ''
     software.vendorWebAddress = 'http://www.extend-nmr.eu'
     #software.details = ''
 


def launchApplication(projectDir=None):

  global top

  root = Tkinter.Tk()
  root.withdraw()
  top = ApplicationPopup(root)

  project = None
  if projectDir:
    projectDir = normalisePath(projectDir)
    try:
      project = loadProject(top, path=projectDir)
    except ApiError, e:
      showError('Reading project', e.error_msg, parent=top)

  top.update_idletasks()

  if project:
    top.initProject(project)

  #root.mainloop()

if __name__ == '__main__':

  argv = sys.argv[:]
  n  = len(argv)

  if n > 1:
    projectDir = argv[1]
  else:
    projectDir = None

  launchApplication(projectDir)

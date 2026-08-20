
# Edited version for protocol testing. We will use two dummy tasks
# Test1Task and Test2Task to test the basic infrastructure for the
# setting up an loggin of tasks. Then switch this back to the
# ExtendNMRFrame.


import os

# added jmci; do we need this??
import Tkinter
from memops.universal.Io import getTopDirectory

#from ccpnmr.analysis.core.ExperimentBasic import getOnebondExpDimRefs
from ccpnmr.analysis.popups.BasePopup import BasePopup

from memops.editor.Util import createDismissHelpButtonList

from memops.gui.Frame import Frame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.ButtonList import ButtonList
from memops.gui.Label import Label
from memops.gui.MessageReporter import showWarning, showInfo
from memops.gui.Entry import Entry
from memops.gui.IntEntry import IntEntry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FileSelect import FileType
from memops.gui.PulldownList import PulldownList
from memops.gui.LabelDivider import LabelDivider
from memops.gui.RadioButtons import RadioButtons

# FIXME JMCI

# THE METHODS. Each method is going to have its own frame.  For now,
# keep the test method to show the overall functionality while we work
# on the remainder

# We also need to decide how to merge advances in the ExtendNMR GUI
# code with changes we make here. Note that we are going to have to
# extend these interfaces a little becuase they will not contain all
# the workflow management widgets that we require

# for testing
from cambridge.wms.Test1Frame import Test1Frame

# ARIA
from paris.aria.AriaExtendNmrFrame import AriaFrame

# AUREMOL
from regensburg.auremol.AuremolFrame import AuremolFrame

# CING
from nijmegen.cing.CingFrame import CingFrame

# HADDOCK
#from utrecht.haddock.HaddockFrame import HaddockFrame

# ISD
from cambridge.isd.IsdFrame import IsdFrame

# PRODEOMP
from gothenburg.prodecomp.ProdecompFrame import ProdecompFrame

# ECI
#from ccpnmr.eci.EntryCompletionFrame import EntryCompletionFrame

# # # # #   B I G   T O - D O  P O I N T S   # # # # # 
#

# Need to decide how and when to open the project file itself. Each
# tab implicitly operates on a single ProjectVersion and the best
# approach seems to be to keep minimal information in the WMS layer,
# just storing a pointer in to the CCP layer.

# # # # #   L A T E R   T O - D O  P O I N T S    # # # # #
#

 
class ExtendNmrPopup(BasePopup):

  def __init__(self, parent, ccpnProject=None):

    self.parent      = parent
    self.ccpnProject = ccpnProject
    
    BasePopup.__init__(self, parent=parent, title='PRODECOMP')

                       
  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    frame = ExtendNmrFrame(guiFrame, basePopup=self)
    frame.grid(row=0, column=0, sticky='nsew')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    
    self.geometry('650x600')
    
    self.update_idletasks()
 

class ExtendNmrFrame(Frame):

  def __init__(self, guiParent, basePopup, project=None):
    # Base popup required to handle notification of data model changes
    # e.g. new peak lists, so that the GUI can update to the latest
    # state
    self.basePopup = basePopup
    self.guiParent = guiParent

    # FIXME JMCI
    # need to work through this far more carefully
    self.project = project

    # shortcuts to this layer
    self.basePopup.frameShortcuts['Task'] = self

    # shortcuts to layer below. We should set these in the
    # initMethod methods
    self.frameShortcuts = {}

    # TODO understand what these are doing
    self.registerNotify=basePopup.registerNotify
    self.unregisterNotify=basePopup.unregisterNotify

    Frame.__init__(self, guiParent)
  
    self.grid_rowconfigure(1, weight=1)
    self.grid_columnconfigure(0, weight=1)

    # Logos. Code nicked from ExtendNMR GUI
    ccpnDir = getTopDirectory()
    
    imageDir = os.path.join(ccpnDir,'python','extendNmr','images')
    

    # basically nicked from Extend-NMR
    softwareOpts = ['Test1Task','ARIA 2','Auremol',
                    'CING',' ECI ','HADDOCK',' ISD ',
                    'PRODECOMP']


    self.tabbedFrame = TabbedFrame(guiParent, options=softwareOpts,
                                   toggleOff=False, selected=0)
    self.tabbedFrame.grid(row=1,column=0,sticky='nsew')
    
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


    # don't know what this does
    self.updateFuncs = []
    self.projButtons = []
  
    # basically nicked from Extend-NMR
    self.initTest1(frames[0])
    self.initAria(frames[1])
    self.initAuremol(frames[2])
    self.initCing(frames[3])
    #self.initEci(frames[4])
    #self.initHaddock(frames[5])
    self.initIsd(frames[6])
    self.initProdecomp(frames[7])
    

  def initTest1(self, frame):
  
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
        
    # FIXME JMCI

    # Provide the basePopup for navigation purposes right now However,
    # there is a strong case for not making the method-specific Frames
    # dependent on the master. It might therefore be better to pass
    # the critical information directly, rather than passing a
    # reference to the basePopup which can then be used to extract it
    # indirectly

    test1Frame = Test1Frame(frame, self.basePopup)
    test1Frame.grid(row=1, column=0, sticky='nsew')

    self.frameShortcuts['Test1'] = test1Frame
    #self.projButtons.extend(cingFrame.buttonBar.buttons)  
    self.updateFuncs.append(test1Frame.updateAll)

  def initAria(self, frame):
  
    welcomeMessage = \
    """ARIA Version 2.2. Authors: Benjamin Bardiaux, Michael Habeck,
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

    self.frameShortcuts['ARIA'] = ariaFrame
    #self.projButtons.extend(ariaFrame.buttons)  
    self.updateFuncs.append(ariaFrame.updateAll)

  def initAuremol(self, frame):
  
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)
    
    auremolFrame = AuremolFrame(frame, self.project)
    auremolFrame.grid(row=0, column=0, sticky='nsew')

    self.frameShortcuts['AUREMOL'] = auremolFrame
    #self.projButtons.extend(isdFrame.buttons)  
    self.updateFuncs.append(auremolFrame.updateAll)


  def initCing(self, frame):
  
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
    
    canvas = Tkinter.Canvas(frame, width=702, height=77) 
    canvas.grid(row=0,column=0)
    canvas.create_image(0,0, anchor='nw', image=self.cingLogo)
    
    cingFrame = CingFrame(frame, self.project)
    cingFrame.grid(row=1, column=0, sticky='nsew')

    self.frameShortcuts['CING'] = cingFrame
    #self.projButtons.extend(cingFrame.buttonBar.buttons)  
    self.updateFuncs.append(cingFrame.updateAll)

    

  def initEci(self, frame):
    
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
    frame.parent = self # For notifiers

    canvas = Tkinter.Canvas(frame, width=753, height=73, bg='#FFFFFF') 
    canvas.grid(row=0,column=0)
    canvas.create_image(10,10, anchor='nw', image=self.msdLogo)
    
    eciFrame = EntryCompletionFrame(frame, basePopup=self)
    eciFrame.grid(row=1, column=0, sticky='nsew')
    
    self.frameShortcuts['ECI'] = eciFrame
    #self.projButtons.extend(isdFrame.buttons)  
    self.updateFuncs.append(eciFrame.updateAll)


  def initHaddock(self, frame):
    
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)
    frame.parent = self # For notifiers

    canvas = Tkinter.Canvas(frame, width=753, height=92) 
    canvas.grid(row=0,column=0)
    canvas.create_image(0,0, anchor='nw', image=self.haddockLogo)
    
    haddockFrame = HaddockFrame(frame, self.project)
    haddockFrame.grid(row=1, column=0, sticky='nsew')

    self.frameShortcuts['HADDOCK'] = haddockFrame
    #self.projButtons.extend(isdFrame.buttons)  
    self.updateFuncs.append(haddockFrame.updateAll)


  def initIsd(self, frame):
      
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)

    canvas = Tkinter.Canvas(frame, width=700, height=91, bg='#FFFFFF') 
    canvas.grid(row=0,column=0,sticky='ew')
    canvas.create_image(0,0, anchor='nw', image=self.isdLogo)
   
    isdFrame = IsdFrame(frame, self.project)
    isdFrame.grid(row=1, column=0, sticky='nsew')

    self.frameShortcuts['ISD'] = isdFrame
    #self.projButtons.extend(isdFrame.buttons)  
    self.updateFuncs.append(isdFrame.updateAll)


  def initProdecomp(self, frame):
  
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(1, weight=1)

    canvas = Tkinter.Canvas(frame, width=700, height=90, bg='#FFFFFF') 
    canvas.grid(row=0,column=0)
    canvas.create_image(0,0, anchor='nw', image=self.prodecompLogo)
     
    prodecompFrame = ProdecompFrame(frame, basePopup=self, ccpnProject=self.project)
    prodecompFrame.grid(row=1, column=0, sticky='nsew')

    self.frameShortcuts['PRODECOMP'] = prodecompFrame
    #self.projButtons.extend(isdFrame.buttons)  
    self.updateFuncs.append(prodecompFrame.updateAll)


  # need to understand this method a little bit better
  def administerNotifiers(self, notifyFunc):

      for func in ('__init__','delete','setName'):
        notifyFunc(self.updateAllAfter, 'ccp.nmr.Nmr.Experiment', func)
        notifyFunc(self.updateAllAfter, 'ccp.nmr.Nmr.DataSource', func)


  def updateAllAfter(self, obj):

    self.after_idle(self.updateAll)


  def updateAll(self, project=None):

    if project:
      self.project = project
      self.nmrProject = project.currentNmrProject

      if not self.nmrProject:
        self.nmrProject = project.newNmrProject(name=project.name)

    if not self.project:
      return


  def quit(self):
  
    self.guiParent.parent.destroy()
    
  def destroy(self):
  
    self.administerNotifiers(self.basePopup.unregisterNotify)
    Frame.destroy(self)


if __name__ == "__main__":

  import sys
  import Tkinter

  root = Tkinter.Tk()
  root.withdraw()
   
  if len(sys.argv) == 2:
    path = sys.argv[1]
    from ccp.gui.Io import loadProject
    ccpnProject = loadProject(root, path=path)
  else:
    ccpnProject = None

  popup = ExtendNmrPopup(root, ccpnProject=ccpnProject)
  
  root.mainloop()

"""
======================COPYRIGHT/LICENSE START==========================

AutoBackbonePopup.py: Part of the CcpNmr Nexus program

Copyright (C) 2003-2010 Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the author: tjs23@cam.ac.uk
=======================================================================
===========================REFERENCE START=============================
===========================REFERENCE END===============================

"""
import Tkinter, os
from os.path import isdir

try:
  import numpy
  HAVE_NUMPY = True
  from ccpnmr.nexus.AutoBackboneNexus import autoBackboneNexus
  
except ImportError:
  HAVE_NUMPY = False


from ccpnmr.analysis.core.AssignmentBasic import findConnectedSpinSystem, assignSpinSystemResidue
from ccpnmr.analysis.core.AssignmentBasic import mergeSpinSystems, mergeResonances
from ccpnmr.analysis.core.AssignmentBasic import assignTentativeSpinSystemResidues
from ccpnmr.analysis.core.AssignmentBasic import getSyntheticShiftLists
from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes
from ccpnmr.analysis.core.ExperimentBasic import getSeqAssignRefExperiments
from ccpnmr.analysis.core.MoleculeBasic import getLinkedResidue, getResidueCode
from ccpnmr.analysis.core.WindowBasic import getActiveWindows, isSpectrumInWindowPane
from ccpnmr.analysis.core.WindowBasic import displayStrips, displaySpinSystemStrips, getWindowPaneName
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.wrappers.Mars import runMars

from ccpnmr.format.converters.MarsFormat import MarsFormat

from ccpnmr.nexus.NexusBasic import guessAtomType, getSpinSystemInterIntraResonances
from ccpnmr.nexus.NexusBasic import getAmideSpinSystems, isResonanceAmide
from ccpnmr.nexus.NexusBasic import linkSpinSystemInterIntraResonances

from memops.editor.Util import createDismissHelpButtonList

from memops.general import Implementation

from memops.gui.Button import Button
from memops.gui.ButtonList import UtilityButtonList, ButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showOkCancel, showError
from memops.gui.ProgressBar import ProgressBar
from memops.gui.PulldownList import PulldownList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.RadioButtons import RadioButtons
from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.ScrolledFrame import ScrolledFrame
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

    
ALLOWED_ATOM_TYPES = ['CA','CB','C','HA','HA2','HA2','HB','HB2','HB3']

MARS_TOLERANCES= {'CO':0.25, 'CA':0.2, 'CB':0.5, 'HA':0.25, 'N':0.5, 'HN':0.3}

NEXUS = 'Nexus'
MARS = 'MARS'

TOLERANCES = {'13C':0.14,'1H':0.08,}

# # # # # # # # # # TBD # # # # # # # # # 
# Consider proton experiments - switch atom types
# Merge CA, CB, CO resonances

PROGRAM_ATTRIBUTES = [(NEXUS, 'nexusPrediction'),
                      (MARS, 'marsPrediction')]

PROGRAM_THRESHOLDS = {NEXUS: (0.75, 0.5),
                      MARS: (0.9, 0.4),}    # proposed by Zweckstetter
                      #MARS: (0.8, 0.4),}

SCORE_COLORS = ['#8080D0','#E0E000','#E00000']

COMPARISONS = ['Current Assignment',
               NEXUS + ' Prediction',
               MARS + ' Prediction']

INSTRUCTIONS = """Before you start its is assumed that:

a) You have amide-linked spectra from the allowed types: HNCA, HNCO, HNCACB,
   HNCB, HNC, HN(ca)CB, HN(ca)CO, HN(co)CA, HN(co)CACB, HN(coca)CB, CBCAcoNH,
   CAcoNH, CcoNH, HNHA, HNHB, HNHAHB, HN(co)HA, HN(co)HB, HN(co)HAHB, HNH-TOCSY,
   C(cco)NH-TOCSY, H(cco)NH-TOCSY.

b) You have picked peaks in your backbone spectra (although in general, genuinely
   missing  peaks will be tolerated).

c) The picked peaks have all been assigned to (anonymous) backbone nitrogen
   resonances, in their amide nitrogen dimensions. 

d) Multiple peaks assigned to the same amide nitrogen relate to the same residue.

e) The assigned amide nitrogen resonances are in spin systems.

Steps b)-e) may be easily achieved by running Assignment::Initialse Root
Resonancess on an HSQC then Assignment::Pick & Asisgn From Roots to transfer
amides and spin systems to backbone experiments.

To run the assignment:

1) Check initial setting: Make sure you are working with the correct molecule 
   chain, chemical shift list and window (to display sequentially assigned 
   strips in).

2) Set the ppm tolerances. Within these tolerances peaks may be attributed to 
   the same resonance, although if two match the closest will be used.

3) Toggle the "Use" column values to specify which spectra will be analysed. 

4) In the backbone Spin Systems tab click [Find Resonance From Peaks]

5) Correct any problems in the spin system table - indicated by red rows.
   Typical issues include: 

   * Having too many peaks e.g. 3 peaks for one amide in HNCA or having 2 CB 
     (according to ppm and sign) peaks in HNcoCA/CB
   * Peaks don't match within tolerances: e.g. an HNcoCA peak doesn't match any
     of the HNCA peaks for the same amide.
   * An assignment has already already been made but it doesn't match the logic
     of matching inter-residue specific peaks, e.g. an HNCA peak is assigned to
     the same residue as the amide, but this peak is closer to the HNcoCA peak
     than the other HNCA peak.

6) In the Automation tab adjust the parameters as required and click [Run Nexus].
   MARS functionality will be incorporated in due course.
   
   For Nexus you shouldn't need to change the number of iterations in normal
   circumstances, but longer runs are allowed just in case you want to eek out a
   bit more info. However, if it doesn't give a reasonable result on the minimum
   setting, going for more iterations probably won't help. All of the convergence
   measurements (which are used as confidence measures together with the peak
   matches) happen before the minimum number of iterations is reached so a user's
   choice cannot affect this. The only thing that can be affected is the 'best
   guess' assignment.

7) When the assignment algorithm is finished you will be presented with a graph of
   assignment reliability. 

8) In the Predictions tab you must confirm all assignments before you commit them.
   Either set a threshold and commit en masse by clicking on  
   [Confirm Above Score Threshold], or confirm your manual selections via the button
   [Confirm Selected]. Use the [Strip Selected Using:] option with a few sequential,
   confirmed assignments to check the tentative assignments in the spectra. Note that
   the adjasent pulldown menu states which spin systems should be used to form the
   strips.

9) Click [Commit Assignments] to do the real assignment of your peaks and spin 
   systems to residues of your molecular chain. To take things more cautiously you can 
   use [Assign Tentatively] which will mark the spin systems with the potential
   residue assignment without making full resonance assignments.
   
10) After having committed the assignments for regions of the protein sequence that
    appear most certain, the Nexus procedure can be repeated with the "Use existing
    assignments" option set on. This may allow convergence in regions that were
    uncertain in previous runs and hence give more complete results.
"""

# To-do
# Update buttons
#

STRIP_ACTIONS = ('Any Available', 'Confirmed Only',
                 'Best Match', '1st Alterative (if available)')

class AutoBackbonePopup(BasePopup):

  def __init__(self, parent, returnCallback=None, title = None, *args, **kw):

    if not title:
      title = "Automated Protein Sequence Assignment"

    self.guiParent = parent
    self.project = parent.getProject()
    self.returnCallback = returnCallback
    self.waiting = False
    self.spinSystems = []
    self.errorCount = 0
    self.spinSystem = None
    self.chain = None
    self.shiftList = None
    self.syntheticShiftList = None
    self.residue = None
    self.windowPane = None
    self.waitingPredictions = False
    self.waitingSpinSystems = False
    self.predictMethod = 'Nexus' #  'MARS'
    self.assignments = {}
    if self.project:
      self.nmrProject = self.project.currentNmrProject
    else:
      self.nmrProject = None
      
    self.allowedRefExps, coRefExpts = getSeqAssignRefExperiments(self.project)
    self.activeLists = {}

    BasePopup.__init__(self, parent=parent, title=title, **kw)
    


  def body(self, guiFrame):

    self.geometry('700x600')

    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)

    options = ['Settings','Spin Systems','Automation',
               'Score Graph','Predictions','Instructions',
               'Comparison']
    self.tabbedFrame = TabbedFrame(guiFrame, options=options, grid=(0,0))
    frameA, frameB, frameC, frameD, frameE, frameF, frameG = self.tabbedFrame.frames
    frameA.expandGrid(2,0)
    frameB.expandGrid(0,0) 
    frameC.expandGrid(2,0) 
    frameD.expandGrid(0,0) 
    frameD.expandGrid(1,0) 
    frameE.expandGrid(1,0) 
    frameF.expandGrid(0,0)    

    #
    # Options
    #

    # General 
    
    frame = LabelFrame(frameA, text='General', grid=(0,0))
    frame.expandGrid(None,6)
    
    label = Label(frame, text='Chain')
    label.grid(row=0, column=0, sticky='w')
    self.chainPulldown = PulldownList(frame, self.changeChain)
    self.chainPulldown.grid(row=0, column=1, sticky='w')
            
    label = Label(frame, text='Shift List:')
    label.grid(row=0, column=2, sticky='w')
    self.shiftListPulldown = PulldownList(frame, self.changeShiftList)
    self.shiftListPulldown.grid(row=0, column=3, sticky='w')

    label = Label(frame, text='Window:')
    label.grid(row=0, column=4, sticky='w')
    self.windowPulldown = PulldownList(frame, self.changeWindow)
    self.windowPulldown.grid(row=0, column=6, sticky='w')


    # Peak Matching Tolerances

    frame = LabelFrame(frameA, text='Peak-Peak Match Tolerances')
    frame.grid_columnconfigure(4, weight=1)
    frame.grid_rowconfigure(0, weight=1)
    frame.grid(row=1, column=0, sticky='nsew')

    label = Label(frame, text='13C (ppm)')
    label.grid(row=0, column=0, sticky='w')

    self.toleranceEntryC = FloatEntry(frame, width=8)
    self.toleranceEntryC.set(TOLERANCES['13C'])     # default 13C tolerance = 0.1 ppm
    self.toleranceEntryC.grid(row=0, column=1, sticky='w')


    label = Label(frame, text='1H (ppm)')
    label.grid(row=0, column=2, sticky='w')

    self.toleranceEntryH = FloatEntry(frame, width=8)
    self.toleranceEntryH.set(TOLERANCES['1H'])     
    self.toleranceEntryH.grid(row=0, column=3, sticky='w')

    # Peak Lists

    frame = LabelFrame(frameA, text='Peak Lists', grid=(2,0))
    frame.expandGrid(0,0)

    headingList      = ['Spectrum','Experiment Type','Use?','Num Peaks']
    editWidgets      = [None,None,None,None]
    editGetCallbacks = [None,None,self.toggleUseBox,None]
    editSetCallbacks = [None,None,None,None]
    self.spectraMatrix = ScrolledMatrix(frame,
                                        headingList=headingList,
                                        editWidgets=editWidgets,
                                        editGetCallbacks=editGetCallbacks,
                                        editSetCallbacks=editSetCallbacks)
    self.spectraMatrix.grid(row=0, column=0, sticky='nsew')

    #self.getBackboneButton = Button(frame, command=self.getAmideSpinSystems, borderwidth=1)
    #self.getBackboneButton.setText('Get Backbone Spin Systems')
    #self.getBackboneButton.grid(row=1, column=0, sticky='ew')

    #
    # Spin Systems
    #


    headingList = ['#','Residue','Intra residue\nresonances',
                   'Inter residue\nresonances','Comments']
    self.spinSystemsMatrix = ScrolledMatrix(frameB, headingList=headingList,
                                            callback=self.selectSpinSystem)
    self.spinSystemsMatrix.grid(row=0, column=0, sticky='nsew')
    
    
    texts = ['Find Resonances\nFrom Peaks','Show Peaks','Go To Root Location']
    commands = [self.setupResonances,self.showPeaks,self.findSpinSystemRoot]
    self.buttonList = ButtonList(frameB, texts=texts, commands=commands)
    self.buttonList.grid(row=1, column=0, sticky='ew')
        
    #
    # Assignment routines
    #

    # Nexus

    frame = LabelFrame(frameC, text='CcpNmr Nexus', grid=(0,0))
    frame.expandGrid(None, 4)

    if not HAVE_NUMPY:
      msg  = 'WARNING: Python NumPy module not installed or accessible.\n'
      msg += 'NumPy is required for Nexus, the CcpNmr automatic backbone assignment routine.'
      label = Label(frame, text=msg, grid=(0,0))
    
    else:

      label = Label(frame, text='Use existing assignments?', grid=(0,2))
      self.nexusKeepExistingSelect = CheckButton(frame, callback=None, grid=(0,3))
      self.nexusKeepExistingSelect.set(False)
      
      label = Label(frame, text='Number of iterations', grid=(1,2))
      
      objects = [10,20,30,40]
      texts = [str(x) for x in objects]
      self.nexusIterationsPulldown = PulldownList(frame, objects=objects, texts=texts,
                                                  index=1, grid=(1,3))
      

      label = Label(frame, text='Gaussian 13C variance (ppm)', grid=(0,0))
      self.nexusVarianceEntryC = FloatEntry(frame, text=0.1, width=8, grid=(0,1))

      label = Label(frame, text='Gaussian 1H variance (ppm)', grid=(1,0))
      self.nexusVarianceEntryH = FloatEntry(frame, text=0.04, width=8, 
                                            grid=(1,1))

      buttonList = ButtonList(frame, commands=[self.runNexus],
                              texts=['Run Nexus',], grid=(0,4), gridSpan=(2,1), 
                              sticky='nsew')

    # Mars
    
    frame = LabelFrame(frameC, text='MARS')
    frame.expandGrid(None, 4)
    self.marsFrame = frame
    
    marsDir = os.environ.get('MARSHOME', '')
    marsExec = os.path.join(marsDir, 'runmars')
    
    if not marsDir:
      msg  = 'WARNING: The MARSHOME environment variable is not set.\n'
      msg += 'Set this to the full path of the directory containing '
      msg += 'your "runmars" file'
      label = Label(frame, text=msg, grid=(0,0))

      if not os.environ.get('PSIPRED_DIR'):
        msg  = 'WARNING: The PSIPRED_DIR environment variable is not set.\n'
        msg += 'Set this to the full path of the directory which '
        msg += 'contains the "bin" and "data" directories for PSIPRED.'
        label = Label(frame, text=msg, grid=(1,0))
    
    elif not os.path.exists(marsExec):
      msg  = 'WARNING: Cannot find MARS program.\n'
      msg += 'Looking at file location: %s\n' % marsExec
      msg += 'Check the MARSHOME environment variable is set to'
      msg += ' the directory containing the "runmars" file.'
      label = Label(frame, text=msg, grid=(0,0))

      if not os.environ.get('PSIPRED_DIR'):
        msg  = 'WARNING: The PSIPRED_DIR environment variable is not set.\n'
        msg += 'Set this to the full path of the directory which '
        msg += 'contains the "bin" and "data" directories for PSIPRED.'
        label = Label(frame, text=msg, grid=(1,0))
     
    elif not os.environ.get('PSIPRED_DIR'):
      msg  = 'WARNING: The PSIPRED_DIR environment variable is not set.\n'
      msg += 'Set this to the full path of the directory which '
      msg += 'contains the "bin" and "data" directories for PSIPRED.'
      label = Label(frame, text=msg, grid=(0,0))
    
    else:

      row = 0
      label = Label(frame, text='Fragment size (residues):', grid=(row,0))
      self.marsfragSizeEntry = IntEntry(frame, text=5, width=8,  grid=(row,1))
      
      label = Label(frame, text='Is protein perdeuterated?', grid=(row,2))
      self.marsDeuteratedSelect = CheckButton(frame, callback=None, grid=(row,3), 
                                            selected=False, sticky='w')

      row += 1
      
      label = Label(frame, text='Keep existing assignments?', grid=(row,0))
      self.marsKeepAssignSelect = CheckButton(frame, callback=None,grid=(row,1),
                                              selected=False, sticky='w')
      
      label = Label(frame, text='Correct shifts for disulfides?', grid=(row,2))
      self.marsDisulfideSelect = CheckButton(frame, callback=None, grid=(row,3), 
                                            selected=False, sticky='w')

      row += 1
                                            
      label = Label(frame, text='Keep existing connectivity?', grid=(row,0))
      self.marsKeepConnSelect = CheckButton(frame, callback=None, grid=(row,1), 
                                            selected=False, sticky='w')
      
      #label = Label(frame, text='Is protein unfolded?', grid=(row,2))
      #self.marsUnfoldedSelect = CheckButton(frame, callback=None, grid=(row,3), 
      #                                      selected=False, sticky='w')
                                            
      #row += 1

      #label = Label(frame, text='Use theoretical shift list?', grid=(row,0))
      #self.syntheticShiftListPulldown = PulldownList(frame, 
      #                                               self.changeSyntheticShiftList, 
      #                                               grid=(row,1))
   
      #row += 1
      #label = Label(frame, text='Solution pH (if unfolded):', grid=(row,2))
      #self.marspHEntry = FloatEntry(frame, text=7.0, width=8,  grid=(row,3))
      
      row += 1
      subframe = Frame(frame, grid=(row,1), gridSpan=(1,3), sticky='w')
      label = Label(subframe, text='CO', grid=(0,0))
      label = Label(subframe, text='CA', grid=(0,1))
      label = Label(subframe, text='CB', grid=(0,2))
      label = Label(subframe, text='HA', grid=(0,3))
      #label = Label(subframe, text=' N', grid=(0,5))
      #label = Label(subframe, text='HN', grid=(0,6))
      label = Label(frame, text='Cutoff (ppm): ', grid=(row,0))
      self.marsCoCutoffEntry = FloatEntry(subframe, grid=(1,0), width=8,  
                                          text=MARS_TOLERANCES.get('CO'))
      self.marsCaCutoffEntry = FloatEntry(subframe, grid=(1,1), width=8,  
                                          text=MARS_TOLERANCES.get('CA'))
      self.marsCbCutoffEntry = FloatEntry(subframe, grid=(1,2), width=8, 
                                          text=MARS_TOLERANCES.get('CB'))
      self.marsHaCutoffEntry = FloatEntry(subframe,  grid=(1,3), width=8, 
                                          text=MARS_TOLERANCES.get('HA'))
      #self.marsNCutoffEntry = FloatEntry(subframe, grid=(1,5), width=8, 
      #                                    text=MARS_TOLERANCES.get('N'))
      #self.marsHnCutoffEntry = FloatEntry(subframe, grid=(1,6), width=8, 
      #                                    text=MARS_TOLERANCES.get('HN'))
      row += 1
      label = Label(frame, text=' ', grid=(row,0))
      row += 1
      label = Label(frame, text='Native MARS has more options, including RDC and coupling input.', grid=(row,0),
                    gridSpan=(1,4))
                          
      row += 1
      #label = Label(frame, text=' ', grid=(row,0))
      #row += 1
      buttonList = ButtonList(frame, commands=[self.runMars],
                              texts=['Run MARS',], grid=(0,4),
                              gridSpan=(6,1), sticky='nsew')
   


    #
    # Graph
    #

    self.scoreGraph = ScrolledGraph(frameD, graphType='histogram',
                                    width=700, height=200)
    self.scoreGraph.grid(row=0, column=0, sticky='nsew')

    self.scoreGraph2 = ScrolledGraph(frameD, graphType='histogram',
                                     width=700, height=200)
    self.scoreGraph2.grid(row=1, column=0, sticky='nsew')

    #
    # Predictions
    #
    
    frame = Frame(frameE, grid=(0,0))
    frame.expandGrid(0,2)
    
    label = Label(frame, text='Show only problem regions', grid=(0,0))
    self.probRegionsCheck = CheckButton(frame,  grid=(0,1),
                                        callback=self.updatePredictionsAfter)

    label = Label(frame, text='Prediction Method:', grid=(0,3))
    self.methodPulldown = PulldownList(frame,  grid=(0,4), texts=(NEXUS, MARS),
                                       index=0, callback=self.selectPredictionMethod)
                                        
    headingList = ['Residue','Assigned\nSpin System',
                   'Confirmed\nSpin System',
                   'Prediction\nScore',
                   'Best Matching\nSpin System',
                   '1st Alternative\nSpin System',
                   '2nd Alternative\nSpin System',]
                   
    self.assignPulldown = PulldownList(self, texts=[])
    editWidgets         = [None, None, self.assignPulldown,
                           None, None, None, None]
    editGetCallbacks    = [None, None, self.getResidueSpinSystems,
                           None, None, None, None]
    editSetCallbacks    = [None, None, self.setResidueSpinSystem,
                           None, None, None, None]
    self.assignmentMatrix = ScrolledMatrix(frameE, 
                                           headingList=headingList,
                                           editWidgets=editWidgets,
                                           editSetCallbacks=editSetCallbacks,
                                           editGetCallbacks=editGetCallbacks,
                                           callback=self.selectResidue,
                                           multiSelect=True, grid=(1,0))

    # Options

    frame = Frame(frameE, grid=(2,0))
    frame.expandGrid(0,0)
    
    texts = ['Clear Confirmation', 'Confirm Selected',
             'Confirm Above Score Threshold:']

    commands = [self.clearSelected,
                self.setSelectedResidues,
                self.setAboveThreshold]
     
    self.confirmButtons = ButtonList(frame, texts=texts,
                                     commands=commands, grid=(0,0))
    
    threshold = PROGRAM_THRESHOLDS[self.predictMethod][0]
    self.assignThresholdEntry = FloatEntry(frame, width=8, grid=(0,1),
                                           text=threshold)
    
    # Buttons

    frame = Frame(frameE, grid=(3,0))
    frame.expandGrid(0,0)

    texts = ['Commit Assignments','Assign Tentatively',
             'Strip Selected Using:']

    commands = [self.commitAssignments, self.tentativeAssignments,
                self.gotoPredictionStrips]
                
    self.assignButtons = ButtonList(frame, texts=texts,
                                    commands=commands, grid=(0,0))

    self.stripPulldown = PulldownList(frame, texts=STRIP_ACTIONS, grid=(0,1),
                                      objects=range(len(STRIP_ACTIONS)))

    #
    # Instructions
    #

    
    frame = LabelFrame(frameF, text='Basic Instructions', grid=(0,0))
    frame.expandGrid(0,0)
    
    sFrame = ScrolledFrame(frame, grid=(0,0))
    
    label = Label(sFrame.frame, text=INSTRUCTIONS, justify='left', grid=(0,0))
      
    # 
    # Comparisons
    #

    frameG.expandGrid(1,4)
    
    label = Label(frameG, text='X Axis:',  grid=(0,0))
    self.compXPulldown = PulldownList(frameG, texts=COMPARISONS, grid=(0,1),
                                      objects=[None, NEXUS, MARS],
                                      callback=self.updateComparison)
    
    label = Label(frameG, text='Y Axis:',  grid=(0,2))
    self.compYPulldown = PulldownList(frameG, texts=COMPARISONS, grid=(0,3),
                                      objects=[None, NEXUS, MARS],
                                      callback=self.updateComparison)
    
    self.compGraph = ScrolledGraph(frameG, graphType='scatter',
                                   width=500, height=500)
                                   
    self.compGraph.grid(row=1, column=0, columnspan=5, sticky='nsew')
    
    # Main

    buttonList = UtilityButtonList(self.tabbedFrame.sideFrame)
    buttonList.grid(row=0, column=0, sticky='e')
     
    self.setDefaultPeakLists()

    self.updateSpectraMatrixAfter()

    self.updateWindowPulldown()

    self.updateChainPulldown()
    
    self.updateShiftLists()

    self.updateSpinSystemsAfter()
    
    if HAVE_NUMPY:
      self.updateIterations()

    self.notifiers(self.registerNotify)  
  
    self.activateMars()  

  def changeShiftList(self, shiftList):
  
    self.shiftList = shiftList

  def changeSyntheticShiftList(self, shiftList):
  
    self.syntheticShiftList = shiftList
  
  def selectPredictionMethod(self, method):
    
    if method != self.predictMethod:
      self.predictMethod = method
      self.updatePredictions()
  
  def activateMars(self):
  
    self.marsFrame.grid(row=1, column=0, sticky='nsew')
    
  
  def open(self):
  
    self.setDefaultPeakLists()

    self.updateSpectraMatrixAfter()

    self.updateWindowPulldown()

    self.updateChainPulldown()
    
    self.updateShiftLists()

    self.updateSpinSystemsAfter()
  
    self.updateIterations()
    
    BasePopup.open(self)
  
  def notifiers(self, notifyFunc):
    
    for func in ('__init__', 'delete', 'setName'):
      notifyFunc(self.updateWindowPulldown, 'ccpnmr.Analysis.SpectrumWindow', func)
      
    for func in ('__init__', 'delete', 'setName'):
      notifyFunc(self.updateWindowPulldown, 'ccpnmr.Analysis.SpectrumWindowPane', func)
    
    for func in ('__init__', 'delete', 'setValue'):
      notifyFunc(self.updateSpinSystemsAfter, 'ccp.nmr.Nmr.PeakDim', func)
        
    for func in ('__init__', 'delete', 'removeResonance', 'addResonance', 'setResonances'):
      notifyFunc(self.updateSpinSystemsAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)

    for func in ('__init__', 'delete',):
      notifyFunc(self.updateSpinSystemsAfter, 'ccp.nmr.Nmr.ResonanceGroupProb', func)
        
    for func in ('delete', 'setAssignName', 'addAssignName'):
      notifyFunc(self.updateSpinSystemsAfter, 'ccp.nmr.Nmr.Resonance', func)
      
    for func in ('__init__', 'delete', 'setName', 'setActivePeakList'):
      notifyFunc(self.updateSpectraMatrixAfter, 'ccp.nmr.Nmr.DataSource', func)
    
    notifyFunc(self.updateSpectraMatrixAfter, 'ccp.nmr.Nmr.Experiment', 'setName')
    notifyFunc(self.updateSpectraMatrixAfter, 'ccp.nmr.Nmr.PeakList', 'delete')
    
    notifyFunc(self.updateShiftLists, 'ccp.nmr.Nmr.Experiment', 'setShiftList')
  
    for func in ('__init__', 'delete'):
      notifyFunc(self.updateChainPulldown, 'ccp.molecule.MolSystem.Chain', func)
    
    for func in ('__init__', 'delete','setResidue','setCcpCode'):
      notifyFunc(self.updatePredictionsAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)
  
  def updateIterations(self):
  
    prev = self.nexusIterationsPulldown.getObject()
  
    if self.chain:
      analysisSize = int(10.0 * (len(self.chain.residues)/50.0)**2) + 1

      objects = []
      v = 10
      while len(objects) < 4:
        if v > analysisSize:
          objects.append(v)
        v += 10   
      
    else:
      objects = [10,20,30,40]
    
    if prev in objects:
      index = objects.index(prev)
    else:
      index = 0
    
    texts = [str(x) for x in objects]                                           
    self.nexusIterationsPulldown.setup(texts, objects, index)
  
  def setSelectedResidues(self):

    attrName = self.getPredictionAttr()
    
    for residue in self.assignmentMatrix.currentObjects:
      spinSystem, score = self.getBestPrediction(residue, attrName)
        
      if spinSystem:
        self.setChosenAssignment(residue, spinSystem)

    self.updatePredictionsAfter()

  def setChosenAssignment(self, residue, spinSystem):
        
    spinSystem2 = self.assignments.get(residue)
    
    if spinSystem2:
      self.assignments[spinSystem2] = None
    
    if spinSystem:
      residue2 = self.assignments.get(spinSystem) 
   
      if residue2 and residue2 is not residue:
        self.assignments[residue2] = None
   
      self.assignments[residue] = spinSystem
      self.assignments[spinSystem] = residue
   
    else:
      self.assignments[residue] = None
				
  def clearSelected(self):
        
    for residue in self.assignmentMatrix.currentObjects:
      self.setChosenAssignment(residue, None)
      
    self.updatePredictionsAfter()  
   
  def setAboveThreshold(self):
    
    try:
      threshold = float(self.assignThresholdEntry.get())
    except:
      showError('Threshold Error', 'Threshold value invalid.', parent=self)
      return
    
    attrName = self.getPredictionAttr()
    
    for residue in self.chain.sortedResidues():
      spinSystem, score = self.getBestPrediction(residue, attrName)
        
      if spinSystem and score > threshold:
        self.setChosenAssignment(residue, spinSystem)
      
    self.updatePredictionsAfter()  
           
  def getPredictionAttr(self, method=None):
  
    if not method:
      method = self.predictMethod
  
    attrName = PROGRAM_ATTRIBUTES[0][1]
    for program, attrNameB in PROGRAM_ATTRIBUTES:
      if program == method:
        attrName = attrNameB
        break

    return attrName

  
  def getBestPrediction(self, residue, attrName):

    if hasattr(residue, attrName):
      score, spinSystems = getattr(residue, attrName)

      return spinSystems[0], score

    return None, None


  def setResidueSpinSystem(self, spinSystem):

    spinSystem = self.assignPulldown.getObject()
    self.setChosenAssignment(self.residue, spinSystem)
    self.updatePredictionsAfter()  
 
 
  def getResidueSpinSystems(self, residue):

    if residue:
      attrName = self.getPredictionAttr()
      spinSystems = [None, ]
 
      if hasattr(residue, attrName):
        score, spinSystemsB = getattr(residue, attrName)
        for spinSystemB in spinSystemsB:
          if spinSystemB:
            spinSystems.append(spinSystemB)
       
      spinSystem = None
      for spinSystemB in residue.resonanceGroups:
        if spinSystemB in spinSystems:
          spinSystem = spinSystemB
	  break

      if not spinSystem:
        spinSystem = self.assignments.get(residue)

      if not spinSystem:
        spinSystem, score = self.getBestPrediction(residue, attrName)
        
      names = ['<None>',] +  [getSpinSystemName(ss) for ss in spinSystems[1:]]
    
      index = 0
      if spinSystem in spinSystems:
        index = spinSystems.index(spinSystem)

      self.assignPulldown.setup(names, spinSystems, index)
         
      
  def updatePredictionsAfter(self, obj=None):
    
    if self.waitingPredictions:
      return
    else:
      self.waitingPredictions = True
      self.after_idle(self.updatePredictions)
  
  def updatePredictions(self):
    
    self.updateScoreGraph()
    self.updateButtons()
    self.updateComparison()
    
    onlyProbRegions = self.probRegionsCheck.get()
    
    textMatrix  = []
    objectList  = []
    colorMatrix = []
    
    if self.chain:
      residues = self.chain.sortedResidues()
    else:
      residues = []  
    
    attrName = self.getPredictionAttr()
    tGood, tBad = PROGRAM_THRESHOLDS[self.predictMethod]
  
    scores = []
    for residue in residues:
      if hasattr(residue, attrName):
        score, spinSystems = getattr(residue, attrName)
        while len(spinSystems) < 3:
          spinSystems.append(None)
      else:
        spinSystems = [None, None, None]
        score = None
        
      scores.append((score, spinSystems))  
      
    for i, residue in enumerate(residues):
      score, spinSystems = scores[i]
      
      colors = [None] * 7

      if score is not None:
        if score > tGood:
          if onlyProbRegions:
          
            prev = None
            if i > 1:
              prev = scores[i-1][0]
            
            next = None
            if i+1 < len(scores):
              next = scores[i+1][0]
          
            if not (((prev is not None) and (prev < tGood)) \
                    or ((next is not None) and (next < tGood))):
              continue  
          index = 0
        elif score < tBad:
          index = 2
        else:
          index = 1
       
        colors[3] = SCORE_COLORS[index]
        score = round(score, 3)
        
      elif onlyProbRegions:
        continue
      
      text = []
      text.append('%d%s' % (residue.seqCode, getResidueCode(residue)))
      
      assSpinSystems = []
      for spinSystemB in residue.resonanceGroups:
        for resonanceB in spinSystemB.resonances:
          if isResonanceAmide(resonanceB):
            assSpinSystems.append(spinSystemB)
            break
      
      if assSpinSystems:
        names = '\n'.join(['{%d}' % ss.serial for ss in assSpinSystems])
        text.append(names)
      else:
        text.append(None)
      
      if not self.assignments.has_key(residue):
        spinSystem = residue.findFirstResonanceGroup()
        self.setChosenAssignment(residue, spinSystem)
      
      spinSystem = self.assignments.get(residue)
     
      if spinSystem:
        text.append(getSpinSystemName(spinSystem))
      else:
        text.append(None)
        
      text.append(score)
      
      for spinSystemB in spinSystems:
        if spinSystemB:
          text.append(getSpinSystemName(spinSystemB))
        else:
          text.append(None)
        
      textMatrix.append(text)
      objectList.append(residue)
      
      if assSpinSystems:
        if spinSystem not in assSpinSystems:
          colors[1] = SCORE_COLORS[2]
 

      colorMatrix.append(colors)
    
    self.assignmentMatrix.update(objectList=objectList,
                                 textMatrix=textMatrix, 
                                 colorMatrix=colorMatrix)
    
    self.waitingPredictions = False    


  def tentativeAssignments(self):

    msg = 'Are you sure you wish to set the chosen tentative assignments?'
    if not showOkCancel('Query', msg, parent=self):
      return

    if not self.chain:
      return

    self.notifiers(self.unregisterNotify)  
    for residue in self.chain.residues:
      spinSystem = self.assignments.get(residue)
      
      if not spinSystem:
        continue
    
      assignTentativeSpinSystemResidues(spinSystem, [residue,])

      # check i-1 and i+1 to merge spin systems if necessary
      pResidue = getLinkedResidue(residue, linkCode='prev')
      nResidue = getLinkedResidue(residue, linkCode='next')

      pSpinSystem = findConnectedSpinSystem(spinSystem)

      if pSpinSystem and pResidue and (pSpinSystem.residue is not pResidue):
        assignTentativeSpinSystemResidues(pSpinSystem, [pResidue,])

    self.notifiers(self.registerNotify)  
    self.updatePredictionsAfter()

  def commitAssignments(self):

    msg = 'Are you sure you wish to commit the chosen assignments?'
    if not showOkCancel('Query', msg, parent=self):
      return

    if not self.chain:
      return

    self.notifiers(self.unregisterNotify)  
    for residue in self.chain.residues:
      spinSystem = self.assignments.get(residue)
      
      if not spinSystem:
        continue
    
      assignSpinSystemResidue(spinSystem, residue)

      # check i-1 and i+1 to merge spin systems if necessary
      pResidue = getLinkedResidue(residue, linkCode='prev')
      nResidue = getLinkedResidue(residue, linkCode='next')

      pSpinSystem = findConnectedSpinSystem(spinSystem)

      if pSpinSystem and pResidue and (pSpinSystem.residue is not pResidue):
        assignSpinSystemResidue(pSpinSystem, pResidue)

     	for ss in pResidue.resonanceGroups:
	  if ss is not pSpinSystem:
	    self.mergeResonances(pSpinSystem, ss)
     	    mergeSpinSystems(pSpinSystem, ss)

      pSpinSystem = None
      if nResidue:
     	for ss in nResidue.resonanceGroups:
     	  pSpinSystem = findConnectedSpinSystem(ss)
     	  break

      if pSpinSystem:
     	self.mergeResonances(pSpinSystem, spinSystem)
     	mergeSpinSystems(pSpinSystem, spinSystem)

    self.notifiers(self.registerNotify)  
    self.updatePredictionsAfter()

  def mergeResonances(self, ss1, ss2):

    typeDict = {}
    
    for resonance in ss1.resonances:
      assignNames = resonance.assignNames
      
      if len(assignNames) == 1: # No ambiguity
        atomType = assignNames[0]
	typeDict[atomType] = typeDict.get(atomType, []) + [resonance,]

    for resonance in ss2.resonances:
      assignNames = resonance.assignNames
      
      if len(assignNames) == 1: # No ambiguity
        atomType = assignNames[0]
	typeDict[atomType] = typeDict.get(atomType, []) + [resonance,]

    for atomType in typeDict.keys():
      resonances = typeDict[atomType]
      
      while len(resonances) > 1:
        resonance1 = resonances[0]
	resonance2 = resonances.pop()
        mergeResonances(resonance1, resonance2)

  def updateComparison(self, obj=None):

    if self.chain:
      residues = self.chain.sortedResidues()
    else:
      residues = []
  
    methodX = self.compXPulldown.getObject()
    methodY = self.compYPulldown.getObject()

    xMap = {}
    yMap = {}

    if methodX is None:
      index = 0
      for spinSystem in self.nmrProject.resonanceGroups:
        residue = spinSystem.residue
        if residue and (residue.chain is self.chain):
          xMap[residue] = spinSystem, index
          xMap[spinSystem] = residue
       
    else:
      attrName = self.getPredictionAttr(methodX)
      tGood, tBad = PROGRAM_THRESHOLDS[methodX]
      for residue in residues:
        if hasattr(residue, attrName):
          score, spinSystems = getattr(residue, attrName)

          if score > tGood:
            index = 0
          elif score < tBad:
            index = 2
          else:
            index = 1
           
          if spinSystems:
            spinSystem = spinSystems[0]
            xMap[residue] = spinSystem, index
            xMap[spinSystem] = residue
    
    if methodY is None:
      index = 0
      for spinSystem in self.nmrProject.resonanceGroups:
        residue = spinSystem.residue
        if residue and (residue.chain is self.chain):
          yMap[residue] = spinSystem, index
          yMap[spinSystem] = residue
    
    else:
      attrName = self.getPredictionAttr(methodY)
      tGood, tBad = PROGRAM_THRESHOLDS[methodY]
      for residue in residues:
        if hasattr(residue, attrName):
          score, spinSystems = getattr(residue, attrName)

          if score > tGood:
            index = 0
          elif score < tBad:
            index = 2
          else:
            index = 1
         
          if spinSystems:
            spinSystem = spinSystems[0]
            yMap[residue] = spinSystem, index
            yMap[spinSystem] = residue
        
    dataSets = [[],[],[]]
    
    for residue in residues:
      spinSystem, i = xMap.get(residue, (None, None))
      
      if spinSystem:
        residue2 = yMap.get(spinSystem)
        
        if residue2:
          dataSets[i].append((residue.seqCode, residue2.seqCode))

    for residue in residues:
      spinSystem, i = yMap.get(residue, (None, None))
      
      if spinSystem:
        residue2 = xMap.get(spinSystem)
        
        if residue2:
          dataSets[i].append((residue.seqCode, residue2.seqCode))
    
    self.compGraph.update(dataSets=dataSets, dataColors=SCORE_COLORS,
                          dataNames=['Good','Mediocre','Bad'], 
                          xLabel=COMPARISONS[self.compXPulldown.index],
                          yLabel=COMPARISONS[self.compYPulldown.index],
                          title='Residue Assignment Comparison')

  def updateScoreGraph(self):
       
    if self.chain:
      residues = self.chain.sortedResidues()
    else:
      residues = []
    
    # Nexus
    
    dataSets = [[],[],[]]
    attrName = self.getPredictionAttr(NEXUS)
    tGood, tBad = PROGRAM_THRESHOLDS[NEXUS]
    
    for residue in residues:
      if not hasattr(residue, attrName):
        score = 0.0
        index = 2
    
      else:
        score, spinSystems = getattr(residue, attrName)
 
        if score > tGood:
          index = 0
        elif score < tBad:
          index = 2
        else:
          index = 1
 
      dataSets[index].append((residue.seqCode, score))
    
    dataNames = ['Good','Mediocre','Bad']
    
    if dataSets[0] or dataSets[1] or dataSets[2]:
      dataSets.append([(0,0),])
    self.scoreGraph.update(dataSets=dataSets,dataNames=dataNames,
                           dataColors=SCORE_COLORS,
                           xLabel='Residue Number', yLabel='Score',
                           title='Nexus Assignment Reliability')
    self.scoreGraph.draw()
    
    # Mars
    
    dataSets = [[],[],[]]
    attrName = self.getPredictionAttr(MARS)
    tGood, tBad = PROGRAM_THRESHOLDS[MARS]
    
    for residue in residues:
      if not hasattr(residue, attrName):
        score = 0.0
        index = 2
      else:
 
        score, spinSystems = getattr(residue, attrName)
 
        if score > tGood:
          index = 0
        elif score < tBad:
          index = 2
        else:
          index = 1
    
      dataSets[index].append((residue.seqCode, score))
    
    dataNames = ['High','Medium','Low']
    
    if dataSets[0] or dataSets[1] or dataSets[2]:
      dataSets.append([(0,0),])
    self.scoreGraph2.update(dataSets=dataSets,dataNames=dataNames,
                           dataColors=SCORE_COLORS,
                           xLabel='Residue Number', yLabel='Score',
                           title='MARS Assignment Reliability')
    self.scoreGraph2.draw()
    
  def selectResidue(self, obj, col, row):
    
    self.residue = obj
    self.updateButtons()

  def runNexus(self):
    
    if not (self.chain and self.shiftList and self.spinSystems):
      return
    
    if self.errorCount:
      msg = 'Spin system setup has produced warnings '
      msg += 'Continue with automated assignment?'
      
      if not showOkCancel('Warning', msg, parent=self):
        return
    
    progressBar = ProgressBar(self, text="Assigning Backbone",
                              total=len(self.spinSystems))
      
    residues = self.chain.sortedResidues()
    peakLists = [pl for pl in self.getPeakLists() if self.activeLists.get(pl)]
    iterations = self.nexusIterationsPulldown.getObject() or 10
    varianceC = self.nexusVarianceEntryC.get() or 0.1
    varianceH = self.nexusVarianceEntryH.get() or 0.04
    keepExisting = self.nexusKeepExistingSelect.get()

    assign = autoBackboneNexus(self.chain, self.spinSystems, self.shiftList,
                               peakLists, iterations, varianceC, varianceH,
                               keepExisting, progressBar)
    
    for residue in assign.keys():
      score, spinSystems = assign[residue]
      residue.nexusPrediction = score, spinSystems
    
    progressBar.destroy()
    
    self.selectPredictionMethod(NEXUS)
    self.methodPulldown.setIndex(0)
    self.updatePredictions()

    self.tabbedFrame.select(3)

  def runMars(self):
    
    if not self.chain and self.shiftList and self.spinSystems:
      return
    
    if self.errorCount:
      msg = 'Spin system setup has produced warnings '
      msg += 'Continue with automated assignment?'
      
      if not showOkCancel('Warning', msg, parent=self):
        return
    
    dataRepository = self.shiftList.root.findFirstRepository(name='userData')
    projPath = dataRepository.url.dataLocation
    
    if not os.access(projPath, os.W_OK):
      msg = 'Current CCPN project directory must be writable '
      msg += 'to store MARS results. Cannot continue'
      
      if showError('Failure', msg, parent=self):
        return
    

    useAssignment = self.marsKeepAssignSelect.get()
    useConnections = self.marsKeepConnSelect.get()
    fragSize = self.marsfragSizeEntry.get() or 5
    isDeuterated = self.marsDeuteratedSelect.get()
    cutoffCO = self.marsCoCutoffEntry.get() or MARS_TOLERANCES.get('CO')
    cutoffCA = self.marsCaCutoffEntry.get() or MARS_TOLERANCES.get('CA')
    cutoffCB = self.marsCbCutoffEntry.get() or MARS_TOLERANCES.get('CB')
    cutoffHA = self.marsHaCutoffEntry.get() or MARS_TOLERANCES.get('HA')
    corrDisulfides = self.marsDisulfideSelect.get()
    # parameters not yet implemented
    #cutoffN  = self.marsNCutoffEntry.get()  or MARS_TOLERANCES.get('N')
    cutoffN = MARS_TOLERANCES.get('N')
    #cutoffHN = self.marsHnCutoffEntry.get() or MARS_TOLERANCES.get('HN')
    cutoffHN = MARS_TOLERANCES.get('HN')
    #isUnfolded = self.marsUnfoldedSelect.get()
    isUnfolded = False
    #pHvalue = marspHEntry.get() or 7.0
    pHvalue = 7.0
    #syntheticShiftList = self.syntheticShiftList
    syntheticShiftList = None
        
    result = runMars(self.shiftList, self.chain, fragSize,
                     cutoffCO, cutoffCA, cutoffCB, cutoffHA, cutoffN, cutoffHN, 
                     useConnections, useAssignment, isDeuterated,
                     isUnfolded, pHvalue, corrDisulfides, syntheticShiftList)
      
    for residue in result.keys():
      score, spinSystems = result[residue] or (0.0, [])
      residue.marsPrediction = score, spinSystems
    
    self.selectPredictionMethod(MARS)
    self.methodPulldown.setIndex(1)
    self.updatePredictions()

    self.tabbedFrame.select(3)
  
  def getChains(self):

    chains = []
    for molSystem in self.project.sortedMolSystems():
      for chain in molSystem.sortedChains():
        for residue in chain.residues:
          if residue.molType == 'protein':
            chains.append(chain)
            break
    
    return chains 
     
  def updateChainPulldown(self, chain=None):
  
    index = 0
    names = []
    chain = self.chain
    chains = self.getChains()
    
    if chains:

      if chain not in chains:
        chain = chains[0]
    
      index = chains.index(chain)
      names = ['%s:%s' % (ch.molSystem.code,ch.code) for ch in chains]
      
    if chain is not self.chain:
      self.chain = chain
      self.assigments = {}
      self.updatePredictionsAfter()  
              
    self.chainPulldown.setup(names,chains,index)
  
  def updateShiftLists(self, experiment=None):
    
    # main shift lists
    peakLists = self.getPeakLists()
    
    shiftLists = []
    listCounts = {}
    for peakList in peakLists:
      shiftList = peakList.dataSource.experiment.shiftList
      listCounts[shiftList] = listCounts.get(shiftList, 0) + 1
      
      if shiftList not in shiftLists:
        shiftLists.append(shiftList)
    
    sortList = [(listCounts[sl], sl) for sl in shiftLists]
    sortList.sort()
    sortList.reverse()
    
    shiftLists = [x[1] for x in sortList] 
    
    index = 0
    names = []
    shiftList = self.shiftList
    
    if shiftLists:
      
      if shiftList not in shiftLists:
        shiftList = shiftLists[0]
        
      index = shiftLists.index(shiftList)
      names = ['%d' % sl.serial for sl in shiftLists]
    
    else:
      shiftList = None
    
    if shiftList is not self.shiftList:
      self.shiftList = shiftList
    
    self.shiftListPulldown.setup(names,shiftLists,index)
    
    # synthetic shift lists (for MARS)
    #shiftLists = getSyntheticShiftLists(self.nmrProject)
    #shiftList = self.syntheticShiftList
    # 
    #if shiftList in shiftLists:
    #  index = shiftLists.index(shiftList)
    #
    #else:
    #  index = 0
    #  self.syntheticShiftList = None
    #  
    #names = ['<None>'] + [x.name or '%d' % x.serial for x in shiftLists]
    #shiftLists = [None] + shiftLists
    #
    #self.shiftListPulldown.setup(names,shiftLists,index)
    
  def changeChain(self, chain):
  
    if chain is not self.chain:
      self.assignments = {}
      self.chain = chain
      self.updatePredictionsAfter()  
      self.updateIterations()
   
  def setDefaultPeakLists(self):
  
    nmrProject = self.project.currentNmrProject
        
    for experiment in nmrProject.experiments:
      refExperiment = experiment.refExperiment
      
      if refExperiment and (refExperiment in self.allowedRefExps):
        for spectrum in experiment.dataSources:
          for peakList in spectrum.peakLists:
            if peakList is spectrum.activePeakList:
              if 'C' in refExperiment.name:
                self.activeLists[peakList] = True
              else:
                self.activeLists[peakList] = False
            else:
              self.activeLists[peakList] = False

   
  def getPeakLists(self):
  
    peakLists = []
  
    nmrProject = self.project.currentNmrProject
        
    for experiment in nmrProject.experiments:
      refExperiment = experiment.refExperiment
      
      if refExperiment and (refExperiment in self.allowedRefExps):
        for spectrum in experiment.dataSources:
          for peakList in spectrum.peakLists:
            peakLists.append(peakList)
    
    return peakLists
    
  def updateSpectraMatrix(self):
    
    objectList       = self.getPeakLists()
    textMatrix       = []
    colorMatrix      = []
    
    for peakList in objectList:
      
      if self.activeLists.get(peakList):
        isUsed = 'Yes'
        colors = (None,None,'#B0FFB0')
      else:
        isUsed = 'No'  
        colors = (None,None,'#FFB0B0')
      
      spectrum = peakList.dataSource
      experiment = spectrum.experiment
      
      listName   = '%s:%s:%d' % (experiment.name, spectrum.name, peakList.serial)
      
      textMatrix.append((listName, experiment.refExperiment.name,
                         isUsed, len(peakList.peaks)))
      colorMatrix.append(colors)
        
    self.spectraMatrix.update(objectList=objectList,
                              textMatrix=textMatrix,
                              colorMatrix=colorMatrix)  
    
    self.updateSpinSystemsAfter()
    self.updateWindowPulldown()
    
  def updateSpinSystemsMatrix(self):
    
    self.errorCount = 0
    
    textMatrix  = []
    objectList  = []
    colorMatrix = []
  
    peakLists = [pl for pl in self.activeLists.keys() if self.activeLists[pl]]
    
    spinSystems = getAmideSpinSystems(peakLists)
        
    self.spinSystems = []
    appendSpinSystem = self.spinSystems.append
  
    for spinSystem in spinSystems:

      comments = []
       
      data = getSpinSystemInterIntraResonances(spinSystem, peakLists)
      pRes, iRes, pPpm, iPpm, pTypes, iTypes = data 
      
      for atomType in ALLOWED_ATOM_TYPES:
        if iTypes.count(atomType) > 1:
          comments.append('Multiple %s' % atomType)
 
        if pTypes.count(atomType) > 1:
          comments.append('Multiple prev %s' % atomType)
 
      intra = []
      if iRes:
        for i, ppm in enumerate(iPpm):
          if iPpm is None:
            continue
          
          if isResonanceAmide(iRes[i]):
            continue
          
          atomType = iTypes[i]
          #if atomType and (atomType not in ALLOWED_ATOM_TYPES):
          #  continue
          
          if not atomType:
            atomType = iRes[i].isotopeCode
            
          intra.append('%s:%.3f' % (atomType, ppm))
          
      inter = []
      if pRes:
        for i, ppm in enumerate(pPpm):
          if pPpm is None:
            continue
          
          if isResonanceAmide(pRes[i]):
            continue
            
          atomType = pTypes[i]
          #if atomType and (atomType not in ALLOWED_ATOM_TYPES):
          #  continue
          
          if not atomType:
            atomType = pRes[i].isotopeCode
            
          inter.append('%s:%.3f' % (atomType, ppm))
      
      intra.sort()
      inter.sort()
      
      #if inter or intra:
      appendSpinSystem(spinSystem)
      
      residue = spinSystem.residue
      if residue:
        resCode = '%d%s' % (residue.seqCode, getResidueCode(residue))
      else:
        resCode = getResidueCode(spinSystem) or None
       
      datum = [getSpinSystemName(spinSystem),
               resCode,
               ' '.join(intra),
               ' '.join(inter),
               ','.join(comments)]
              
      textMatrix.append(datum)
      
      if not comments:
        colorMatrix.append([None] * 5)
      else:
        self.errorCount += 1 
        colorMatrix.append(['#FFB0B0']*5)
      
      objectList.append(spinSystem)
    
    if not spinSystems:
      self.errorCount = 1
    
    self.spinSystemsMatrix.update(textMatrix=textMatrix,
                                  objectList=objectList,
                                  colorMatrix=colorMatrix)
      
  def setupResonances(self):
  
    TOLERANCES['1H'] = self.toleranceEntryH.get() or 0.04
    TOLERANCES['13C'] = self.toleranceEntryC.get() or 0.1
  
    msg  = 'Continue and add any missing intra-'
    msg += ' and inter-residue resonances'
    msg += ' using selected peak lists?'
    msg += ' (Saving your project first may be wise)'
    
    if not showOkCancel('Query', msg, parent=self):
      return 

    self.notifiers(self.unregisterNotify)  
    
    spinSystems = self.spinSystemsMatrix.objectList
    
    progressBar = ProgressBar(self, text="Curating intra and inter residue resonances",
                              total=len(spinSystems))
    
    peakLists = set()
    for peakList in self.activeLists:
      if self.activeLists[peakList]:
        peakLists.add(peakList)
                              
    for spinSystem in spinSystems:
      linkSpinSystemInterIntraResonances(spinSystem, peakLists, TOLERANCES)
      progressBar.increment()

    progressBar.destroy()
    
    self.notifiers(self.registerNotify)  
    self.updateSpinSystemsAfter()
             
             
  def toggleUseBox(self, peakList):
    
    bool = self.activeLists.get(peakList, False)
    self.activeLists[peakList] = not bool
              
    self.updateSpectraMatrixAfter()
    self.updateSpinSystemsAfter()
    
    
  def selectSpinSystem(self, obj, row, col):
    
    self.spinSystem = obj
    
    self.updateButtons()
         
  
  def changeWindow(self, windowPane):
  
    if windowPane is not self.windowPane:
      self.windowPane = windowPane
      # Maybe updates
         
  def showPeaks(self):

    spinSystem = self.spinSystem

    peaks = set()
    for resonance in spinSystem.resonances:
      #if resonance.isotopeCode != '15N':
      #  continue
    
      #if not isResonanceAmide(resonance):
      #  continue
    
      for contrib in resonance.peakDimContribs:
        peak = contrib.peakDim.peak
        
        if self.activeLists.get(peak.peakList):
          peaks.add(peak)

    if peaks:
      self.guiParent.viewPeaks(list(peaks))    

  def gotoPredictionStrips(self):

    index = self.stripPulldown.getObject()
    
    options = ((True, True, False),
               (True, False, False),
               (False, True, False),
               (False, True, True))
    
    useConfirmed, useMatch, useAlternative = options[index]
 
    self.gotoSpinSystemStrips(useConfirmed, useMatch, useAlternative)
 

  def gotoSpinSystemStrips(self, useConfirmed=True, useMatch=True,
                           useAlternative=False):
  
    windowPane = self.windowPane
    if not windowPane:
      return
    
    shiftList = self.shiftList  
    axisH = None
    axisN = None
    for axisPanel in windowPane.sortedAxisPanels():
      axisType = axisPanel.axisType
 
      if axisType:
        if '1H' in axisType.isotopeCodes:
          axisH = axisPanel.label
 
        elif '15N' in axisType.isotopeCodes:
          axisN = axisPanel.label
    
    if not (axisH and axisN):
      return
    
    attrName = self.getPredictionAttr()
    
    positions = []
    for residue in self.assignmentMatrix.currentObjects:
      if useConfirmed:
        spinSystem = self.assignments.get(residue)
      else:
        spinSystem = None
      
      if (not spinSystem) and hasattr(residue, attrName):
        score, spinSystems = getattr(residue, attrName)
        
        if useAlternative and (spinSystems[1] is not None):
          spinSystem = spinSystems[1]
        elif useMatch:
          spinSystem = spinSystems[0]
      
      if not spinSystem:
        continue
      
      shiftN = None
      shiftH = None
      for resonance in spinSystem.resonances:
        if resonance.isotopeCode != '15N':
          continue
          
        if not isResonanceAmide(resonance):
          continue  
      
        bound = resonance.findFirstCovalentlyBound(isotopeCode='1H')
        shiftN = resonance.findFirstShift(parentList=shiftList)
        shiftH = bound.findFirstShift(parentList=shiftList)
        
        if shiftH and shiftN:
          position = {axisH:shiftH.value, axisN:shiftN.value}
          positions.append(position) 
    
    if positions:
    
      displayStrips(self.parent, positions, windowPane=windowPane)
      self.update_idletasks()
      displayStrips(self.parent, positions, windowPane=windowPane)
        
        
  def findSpinSystemRoot(self):
    
    if self.spinSystem and self.windowPane:
      spinSystems = [self.spinSystem,]
      displaySpinSystemStrips(self.parent, spinSystems, self.windowPane)


  def updateWindowPulldown(self, win=None):
  
    windowPane = self.windowPane
    windowPanes = self.getWindows()
    names = []
    index = 0

    if windowPanes:
      if windowPane not in windowPanes:
        windowPane = windowPanes[0]
    
      names = [getWindowPaneName(w) for w in windowPanes]
      index = windowPanes.index(windowPane)
    
    else:
      windowPane = None
    
    self.changeWindow(windowPane)  
    
    self.windowPulldown.setup(names, windowPanes, index)

  def getWindows(self):
  
    windowPanes = []
    for peakList in self.activeLists:
      project = peakList.root
      spectrum = peakList.dataSource
      tryWindows = getActiveWindows(project)
      
      for window in tryWindows:
        for windowPane in window.sortedSpectrumWindowPanes():
          if isSpectrumInWindowPane(windowPane, spectrum):
            if windowPane not in windowPanes:
              windowPanes.append(windowPane)
    
    return windowPanes
      
  def updateSpectraMatrixAfter(self, obj=None):  

    if obj:
      if obj.className == 'DataSource':
        refExperiment = obj.experiment.refExperiment
	if not refExperiment:
	  return
	
	if refExperiment not in self.allowedRefExps:
          return
 
      elif obj.className == 'Experiment':
        if obj.refExperiment not in self.allowedRefExps:
          return
 
      elif obj.className == 'PeakList':
        expType = obj.dataSource.experiment.refExperiment
        if expType not in self.allowedRefExps:
          return
 
        if obj.isDeleted and (self.activeLists.has_key(obj)):
          del self.activeLists[obj]
          self.updateSpinSystemsAfter() # was being used
          return
     
    self.after_idle(self.updateSpectraMatrix)
    
    
  def updateButtons(self):
    
    buttons = self.buttonList.buttons
    
    # check spinSystemsmMatrix
    #if self.errorCount == 0:
    #  buttons[0].enable()

    #else:
    #  buttons[0].disable()

    if self.spinSystem:
      buttons[1].enable()
      
      if not self.windowPane:
        buttons[2].disable()
     
      else:
        buttons[2].enable()
      
    else:
      buttons[1].disable()
      buttons[2].disable()

  def updateSpinSystemsAfter(self, obj=None):

    if self.waitingSpinSystems:
      return
    
    if obj:
      if obj.className == 'PeakDim':
        if not self.activeLists.has_key(obj.peak.peakList):
          return
          
      elif obj.className == 'ResonanceGroup':
        if obj not in self.spinSystems:
          return
      
      elif obj.className == 'Resonance':
        if obj.resonanceGroup not in self.spinSystems:
          return           

      elif obj.className == 'ResonanceGroupProb':
        if not obj.isSelected:
          return           
        if obj.linkType != 'sequential':
          return           
        if abs(obj.sequenceOffset) != 1:
          return           

      else:
        return

    self.waiting = True
    self.after_idle(self.updateSpinSystems)
    
    
  def updateSpinSystems(self):
    
    self.updateSpinSystemsMatrix()
    self.updateButtons()
    
    self.waitingSpinSystems = False
    

  def destroy(self):
      
    self.notifiers(self.unregisterNotify)  
      
    BasePopup.destroy(self)

def getSpinSystemName(spinSystem):

  residue = spinSystem.residue
  
  if residue:
    assignStr = ':%d%s' % (residue.seqCode, getResidueCode(residue))
  else:
    assignStr = ''

  return '{%s}%s' % (spinSystem.serial, assignStr)

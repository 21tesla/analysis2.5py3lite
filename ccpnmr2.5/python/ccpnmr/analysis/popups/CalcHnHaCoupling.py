
"""
======================COPYRIGHT/LICENSE START==========================

CalcHnHaCoupling.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.LabelDivider import LabelDivider
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showWarning, showOkCancel, showError
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.core.AssignmentBasic import assignAtomsToRes, assignResToDim
from ccpnmr.analysis.core.ChemicalShiftBasic import getAtomProbability
from ccpnmr.analysis.core.ConstraintBasic import makeNmrConstraintStore, getFixedResonance
from ccpnmr.analysis.core.CouplingBasic import getResidueJCoupling, setResidueJCoupling
from ccpnmr.analysis.core.ExperimentBasic import getOnebondDataDims, getSpectrumNoise
from ccpnmr.analysis.core.MarkBasic import createPeakMark
from ccpnmr.analysis.core.MoleculeBasic import getRandomCoilShift, getLinkedResidue
from ccpnmr.analysis.core.WindowBasic import getWindowPaneName, isSpectrumInWindowPane, getActiveWindows, zoomToShowPeaks
from ccpnmr.analysis.popups.BasePopup import BasePopup


from math import sqrt, atan, acos, radians, cos, exp

ALLOWED_REF_EXPTS = ['H[{N,HA}]','H{[N]+[HA]}']

ATOM_NAME_DICT = {'H[{N,HA}]':('H','N','CA','HA'),
                  'H{[N]+[HA]}':('H','N','CA','HA')}

ISOTOPE_DICT = {'H':'1H','N':'15N','C':'13C'}

FALLBACK_REF_SHIFTS = {'H':(8.3, 0.6),'HA':(4.2, 0.5)}

twoPi = 6.2831853071795864

def testCalcHnHaCouplingPopup(argServer):

  popup = CalcHnHaCouplingPopup(argServer.parent)
  popup.open()

# Future: Maybe show Karplus curve and predicted angles upon it

  
class CalcHnHaCouplingPopup(BasePopup):
  """
  **Calculate 3J HN-HA Couplings Using an HNHA Experiment**
  
  This system is designed to efficiently extract 3J coupling information from
  the relative intensities of HA and amide peaks recorded in a quantitative HNHA
  experiment. The couplings may be used, according to the Karplus relationship,
  to estimate the Phi protein backbone dihedral angle.

  The layout is broken down into two tabs: one for the results in the "Spin
  System Table" and an other for setting the various options. Although the user
  can adjust the options at any time the system will automatically perform the
  calculations and fill the results table using the first HNHA peak list it
  comes across.
  
  For the calculation to be performed the experiment must be set to the specific
  type "H{[N]+[HA]}", which may be set in the Experiment Types tab of the main
  Experiments_ table. It is also assumed that a peak list has been picked in the
  HNHA spectrum, for both the HA cross peaks and the homonuclear diagonal amide
  peaks. For a peak within the selected list to be considered it must be
  assigned to resonances in its amide dimensions *and a spin system*, for
  example using the `Pick & Assign From Roots`_ tool. This spin system
  assignment is the means to achieve the pairing of a diagonal amide peak with a
  corresponding HA peak; if two peaks for the same residue are not in the same
  spin system then the analysis will not be made. The peaks do not have to be
  assigned to particular residues, i.e. they can have anonymous spin system
  numbers. It is also unnecessary to assign the non-amide dimension (e.g. HA
  resonance) of the peaks; the amide and HA peak are automatically distinguished
  by sign and ppm location.

  In the "options" tab the user can switch between different HNHA peak lists,
  control how output data will be made; like where dihedral restraints of
  coupling measurements go (which restraint set), and setup the parameters for
  the estimation of Phi dihedral angles. Changing the Karplus coefficients, to
  change the below curve, will affect the prediction of the Phi dihedral angle.
  The user may wish to adjust this but the defaults are values commonly used
  for protein structures.   

  In the main results table clicking on a particular row, assuming the follow
  option and spectrum is selected, will cause the display to show the location
  of peaks for the selected spin system. The 3J coupling values are
  automatically extracted from the relative intensities if the HN (diagonal) and
  HA peaks, according to the method referred to below, whenever the table is
  updated. The coupling measurements may be saved in the CCPN project by using
  [Make Coupling List] - this creates a kind of measurement list which may be
  inspected at any time, without having to repeat any calculations. Making
  "Coupling Restraints" is specifically for 3D structure calculations that know
  how to interpret such information. The dihedral angle estimates that are
  presented are made according to the Karplus curve displayed in the "Options",
  noting that because of its oscillatory nature there may be more than one range
  of possible angles for a particular 3J coupling. The "Only likely angles"
  option can reduce the ambiguity in Phi angle prediction by allowing only
  values that are in the common regions of the Ramachandran plot. Angle
  predictions may be converted into dihedral angle restrains, for use in
  structure calculations, with the [Make Dihedral Restraints] button. 

  **Caveats & Tips**
  
  Because this method relies upon peak intensities, the user should be cautious
  in regions of the spectrum where peak overlap is significant enough to affect
  the size and shape of peaks. Where there is overlap using the peak 'height' 
  may perform better than 'volume', but this will not overcome the problem
  entirely.

  Residues may be excluded from restraint and coupling list generation by
  double-clicking in the "Use?" column.

  If a residue or spin system appears to be missing, check that the correct peak
  list is selected (there could be many for one spectrum) and that both peaks
  for a residue are assigned to a spin system in their amide dimensions. 

  No special provision is made for glycine residues, but they may be included
  in the analysis using the "Show Glycines" option.

  **Reference**
  
  The method used by this system closely follows that which is described in the
  following reference:
  
  *G. W. Vuister and A. Bax (1993). "Quantitative J correlation: a new approach
  for measuring homonuclear three-bond J(HNHa) coupling constants in
  15N-enriched proteins". J. Am. Chem. Soc. 115 (17): 7772-7777*
  
  .. _Experiments: EditExperimentPopup.html
  .. _`Pick & Assign From Roots`: LinkPeakListsPopup.html
  
  """

  def __init__(self, parent, *args, **kw):
    
    self.waiting       = False
    self.peakList      = None
    self.jCouplingList = None
    self.constraintSet = None
    self.spinSystem    = None
    self.peakPair      = None
    self.guiParent     = parent
    self.couplingData  = []
    self.windowPane = None
    
    BasePopup.__init__(self, parent=parent, title=u'Data Analysis : 3J H-H\u03B1 Coupling', **kw)

  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)

  def close(self):

    self.couplingData = []
    BasePopup.close(self)

  def body(self, guiFrame):

    self.geometry('600x500')
   
    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)
    
    options = ['Spin System Table','Options']
    tabbedFrame = TabbedFrame(guiFrame, options=options)
    tabbedFrame.grid(row=0, column=0, sticky='nsew')
    self.tabbedFrame = tabbedFrame
    frameB, frameA = tabbedFrame.frames
    
    row = 0
    frame = Frame(frameA, grid=(row,0))
    frame.grid_columnconfigure(8, weight=1)

    tipText = 'Selects the HNHA peak list to calculate couplings from'
    label = Label(frame, text='Peak List:')
    label.grid(row=0, column=0, sticky='e')
    self.peakListPulldown = PulldownList(frame, callback=self.changePeakList, tipText=tipText)
    self.peakListPulldown.grid(row=0, column=1, sticky='w')

    tipText = 'Selects whether to compare peak heights or volume integrals'
    label = Label(frame, text=' Intensity Type:')
    label.grid(row=0, column=2, sticky='e')
    self.intensityPulldown = PulldownList(frame, texts=['height','volume'],
                                          callback=self.changeIntensity,
                                          tipText=tipText)
    self.intensityPulldown.grid(row=0, column=3, sticky='w')

    tipText = 'The transfer period used in the HNHA experiment'
    label = Label(frame, text=' Transfer Time:')
    label.grid(row=0, column=4, sticky='e')
    self.timeEntry = FloatEntry(frame, text=0.01305, width=8, tipText=tipText)
    self.timeEntry.grid(row=0, column=5)

    tipText = 'Correction factor to compensate for differential relaxation'
    label = Label(frame, text=' Relax Correction:')
    label.grid(row=0, column=6, sticky='e')
    self.relaxCorrEntry = FloatEntry(frame, text=1.11, width=6, tipText=tipText, grid=(0,7))

    row += 1
    frame = Frame(frameA, grid=(row,0))
    frame.grid_columnconfigure(6, weight=1)

    tipText = 'Selects a J-coupling list to store results in'
    label = Label(frame, text='JCoupling List:')
    label.grid(row=1, column=0, sticky='e')
    self.jCouplingListPulldown = PulldownList(frame, callback=self.changeJCouplingList, tipText=tipText)
    self.jCouplingListPulldown.grid(row=1, column=1, sticky='w')

    tipText = 'Selects a restraint set to make dihedral (phi) angle restraints in'
    label = Label(frame, text=' Restraint Set:')
    label.grid(row=1, column=2, sticky='e')
    self.constraintSetPulldown = PulldownList(frame, callback=self.changeConstraintSet, tipText=tipText)
    self.constraintSetPulldown.grid(row=1, column=3, sticky='w')

    tipText = 'Selects a spectrum window to display peak locations in'
    self.windowLabel = Label(frame, text=' Window:')
    self.windowLabel.grid(row=1, column=4, sticky='e')
    self.windowPulldown = PulldownList(frame, self.selectWindowPane, tipText=tipText)
    self.windowPulldown.grid(row=1, column=5, sticky='w')

    row += 1
    frame = Frame(frameA, grid=(row,0))
    frame.grid_columnconfigure(6, weight=1)
 
    tipText = 'Whether to show Glycine residues in spin system table: defaults off due to two HAs'
    label = Label(frame, text='Show Glycines:')
    label.grid(row=2, column=0, sticky='e')
    self.glycineSelect = CheckButton(frame, callback=self.updateAfter, tipText=tipText)
    self.glycineSelect.grid(row=2, column=1, sticky='w')
    self.glycineSelect.set(False)

    tipText = 'Whether to mark peak positions in windows with lines'
    label = Label(frame, text='Mark Peaks:')
    label.grid(row=2, column=2, sticky='e')
    self.markSelect = CheckButton(frame, callback=None, tipText=tipText)
    self.markSelect.grid(row=2, column=3, sticky='w')
    self.markSelect.set(True)

    tipText = 'Whethe rto follow peak positions in a spectrum window when clicking in the spin system table'
    label = Label(frame, text='Follow Peaks:')
    label.grid(row=2, column=4, sticky='e')
    self.followSelect = CheckButton(frame, callback=None, tipText=tipText)
    self.followSelect.grid(row=2, column=5, sticky='w')
    self.followSelect.set(True)

    row += 1
    div = LabelDivider(frameA, text=u'\u03A6 Dihedral Angle Prediction')
    div.grid(row=row, column=0, sticky='ew')
    
    row += 1
    frame = Frame(frameA)
    frame.grid(row=row, column=0, sticky='ew')
    frame.grid_columnconfigure(10, weight=1)
    frame.grid_rowconfigure(0, weight=1)

    tipText = 'The first coefficient in the Karplus equation; for the cosine^2(angle) term'
    label = Label(frame, text='Karplus Coefficients A:')
    label.grid(row=0, column=0, sticky='w')
    self.karplusAEntry = FloatEntry(frame, text=6.51, width=6, tipText=tipText,
                                    returnCallback=self.updateAfter)
    self.karplusAEntry.grid(row=0, column=1, sticky='ew')
    self.karplusAEntry.bind('<Leave>', self.updateKarplusCurve, '+')

    tipText = 'The second coefficient in the Karplus equation; for the cosine(angle) term'
    label = Label(frame, text='B:')
    label.grid(row=0, column=2, sticky='w')
    self.karplusBEntry = FloatEntry(frame, text=-1.76, width=6, tipText=tipText,
                                    returnCallback=self.updateAfter)
    self.karplusBEntry.grid(row=0, column=3, sticky='ew')
    self.karplusBEntry.bind('<Leave>', self.updateKarplusCurve, '+')

    tipText = 'The third coefficient in the Karplus equation; the scalar constant'
    label = Label(frame, text='C:')
    label.grid(row=0, column=4, sticky='w')
    self.karplusCEntry = FloatEntry(frame, text=1.60, width=6, tipText=tipText,
                                    returnCallback=self.updateAfter)
    self.karplusCEntry.grid(row=0, column=5, sticky='ew')
    self.karplusCEntry.bind('<Leave>', self.updateKarplusCurve, '+')

    tipText = 'Whether to show only predicted phi angles for common protein conformations'
    label = Label(frame, text='Only likely angles:')
    label.grid(row=0, column=6, sticky='w')
    self.angleFilterSelect = CheckButton(frame, callback=self.updateAfter, tipText=tipText)
    self.angleFilterSelect.grid(row=0, column=7, sticky='w')
    self.angleFilterSelect.set(True)

    tipText = 'Width to derive upper and lower bound for dihedral angle restraints'
    label = Label(frame, text='Angle error:')
    label.grid(row=0, column=8, sticky='w')
    self.angleErrorEntry = FloatEntry(frame, text=30.0, width=6,
                                      returnCallback=None, tipText=tipText)
    self.angleErrorEntry.grid(row=0, column=9, sticky='ew')
 
    row += 1
    div = LabelDivider(frameA, text='Karplus Curve', grid=(row,0))
    
    row += 1
    frameA.expandGrid(row, 0)
    self.karplusCurve = ScrolledGraph(frameA, title='Karplus Curve', width=400,
                                      height=200, xLabel='Angle (degrees)',
                                      yLabel='3J Coupling (Hz)', grid=(row,0),
                                      symbolSize=1, graphType='scatter',
                                      dataColors=['#FF8080','#000000'],)


    #
    # Spin system table
    #
    
    frameB.expandGrid(0, 0)

    tipTexts = ['Assigned name of the root (amide) spin system',
                'Whether or not to use a given spin system when recording restraints & couplings',
                'The three bond amide H to alpha H coupling in Hz',
                'The estimated error in the three-bond HN-HA coupling',
                'Difference between the observed and (sequence adjusted) random coil HA chemical shift',
                'Predicted phi backbone torsion angles, in degrees, according to the Karplus equation',
                'Ratio of peak intensities used in calculation: HA crosspeak over diagonal',
                'Intensity of amide H diagonal peak',
                'Intensity of alpha H crosspeak']
    headingList = ['Root Spin\nSystem','Use?',
                   u'3J[H,H\u03B1]',u'Error\n3J[H,H\u03B1]',
                   u'\u0394\u03B4H\u03B1',u'\u03A6\nAngles',
                   'Intensity\nRatio','Amide\nIntensity',u'\u03B1\nIntensity']
    editSetCallbacks = [None] * 8
    editGetCallbacks = [None] * 8
    editWidgets      = [None] * 8
    editGetCallbacks[1] = self.toggleConfirm
    self.spinSystemMatrix = ScrolledMatrix(frameB, headingList=headingList,
                                           editSetCallbacks=editSetCallbacks,
                                           editGetCallbacks=editGetCallbacks,
                                           editWidgets=editWidgets,
                                           callback=self.selectSpinSystem,
                                           multiSelect=True, tipTexts=tipTexts)
    self.spinSystemMatrix.grid(row=0, column=0, sticky='nsew')
                                        
    tipTexts = ['Assign amide H and alpha H resonances to the indirect 1H peak dimensions',
                'Make a J-coupling measurement list from results; stored in CCPN project',
                'Make a dihedral angle restraint list from results',
                'Make a J-coupling restraint list from results',
                'Manually update the table & redo calculations; useful after adjusting intensities etc.']
    texts     = ['Assign non-root\ndimensions','Make Coupling\nList',
                 'Make Dihedral\nRestraints','Make Coupling\nRestraints',
                 'Update\nTable']
    commands  = [self.assignNonRootDim,self.makeJCouplingList,
                 self.makeAngleConstraints,self.makeCouplingConstraints,self.updateAfter]
    self.tableButtons = ButtonList(frameB, texts=texts, commands=commands,
                                   tipTexts=tipTexts, grid=(1,0))

    bottomButtons = UtilityButtonList(tabbedFrame.sideFrame, expands=True, helpUrl=self.help_url)
    bottomButtons.grid(row=0, column=0, sticky='e')
          
    self.administerNotifiers(self.registerNotify)

    self.inBody = True
    self.updatePeakLists()
    self.updateJCouplingLists()
    self.updateConstraintSets()
    self.updateKarplusCurve()
    self.updateButtons()
    self.inBody = False

  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.PeakDimContrib',):
        notifyFunc(self.contribUpdateAfter, clazz, func)

    for func in ('__init__', 'delete','setName'):
      for clazz in ('ccp.nmr.Nmr.PeakList',
                    'ccp.nmr.Nmr.DataSource',
                    'ccp.nmr.Nmr.Experiment'):
        notifyFunc(self.updatePeakListsAfter, clazz, func)

    notifyFunc(self.updatePeakListsAfter, 'ccp.nmr.Nmr.Experiment','setRefExperiment')

    for func in ('__init__','delete','setName'):
      notifyFunc(self.updateWindowListAfter, 'ccpnmr.Analysis.SpectrumWindow', func)

    for func in ('__init__','delete','setName'):
      notifyFunc(self.updateWindowListAfter, 'ccpnmr.Analysis.SpectrumWindowPane', func)

    for func in ('__init__','delete',):
      notifyFunc(self.peakUpdateAfter, 'ccp.nmr.Nmr.Peak', func)
      
    for func in ('setAnnotation','setPosition','setNumAliasing'):
      notifyFunc(self.peakDimUpdateAfter, 'ccp.nmr.Nmr.PeakDim', func)

    for func in ('__init__','delete','setValue'):
      notifyFunc(self.peakDimUpdateAfter, 'ccp.nmr.Nmr.PeakIntensity', func)

    for func in ('__init__', 'delete' 'setResidue', 'setName',
                 'setResonances', 'addResonance',
                 'removeResonance','setCcpCode', ):
      notifyFunc(self.updateSpinSystemAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)

    for func in ('__init__','delete'):
      notifyFunc(self.updateConstraintSets, 'ccp.nmr.NmrConstraint.NmrConstraintStore', func)

    for func in ('__init__','delete','setName'):
      notifyFunc(self.updateJCouplingLists, 'ccp.nmr.Nmr.JCouplingList', func)
  
  def updateKarplusCurve(self, event=None):
  
    kA, kB, kC = self.getCoefficients()
    
    if not self.couplingData:
      self.couplingData = couplingData = self.extractHnHaCouplingData()
    else:
      couplingData = self.couplingData
    
    dataSet = []
    dataSetAppend = dataSet.append
    
    for angle in range(-180,181,2):
      c = cos(radians(angle-60))
      y = (kA*c*c) + (kB*c) + kC
      dataSetAppend((angle,y))

    dataSetB = []
    dataSetBAppend = dataSetB.append
    for ss, peaks, coupling, couplingError, secShift, phis, ratio, iH, iHa, in  couplingData:
      for phi in phis:
        dataSetBAppend((phi, coupling))
    
    self.karplusCurve.update([dataSet])
    
  def changeIntensity(self, intensityType):
  
    if self.waiting or not self.peakList:
      return
  
    self.updateAfter()
    
  def getCoefficients(self):
  
    kA = self.karplusAEntry.get() or 6.51
    kB = self.karplusBEntry.get() or -1.76
    kC = self.karplusCEntry.get() or 1.60
  
    return kA, kB, kC
  
  def toggleConfirm(self, obj):
  
    spinSystem, peakPair, c, c1, c2 = obj
  
    boolean = not spinSystem.calcHnHaConfirmed
    
    spinSystem.calcHnHaConfirmed = boolean

    self.update(recalc=False)  
  
  def contribUpdateAfter(self, contrib):
  
    if self.waiting:
      return
  
    if contrib.peakDim.peak.peakList is self.peakList:
      self.updateAfter() 
  
  def updateSpinSystemAfter(self, spinSystem):
  
    if self.waiting:
      return
    
    for resonance in spinSystem.resonances:
      for contrib in resonance.peakDimContribs:
        peakList = contrib.peakDim.peak.peakList
      
        if peakList is self.peakList:
          self.updateAfter()
          return
    
  def updatePeakListsAfter(self, obj):
  
    self.after_idle(self.updatePeakLists)

  def updateWindowListAfter(self, obj):

    self.after_idle(self.updateWindowList)
  
  def updateWindowList(self):
  
    index = 0
    names = []
    windowPane = self.windowPane
    windowPanes = self.getWindows()
    
    if windowPanes:
      names = [getWindowPaneName(wp) for wp in windowPanes]
      if windowPane not in windowPanes:
        windowPane = windowPanes[0]
      
      index = windowPanes.index(windowPane)  
    
    else:
      windowPane = None
    
    if windowPane is not self.windowPane:
      self.selectWindowPane(windowPane)
    
    self.windowPulldown.setup(names, windowPanes, index)
  
  def selectWindowPane(self, windowPane):
  
    if windowPane is not self.windowPane:
      self.windowPane = windowPane
       
     
  def peakUpdateAfter(self, peak):
 
    if self.waiting:
      return
  
    if peak.peakList is self.peakList:
      self.updateAfter() 
  
  def peakDimUpdateAfter(self, peakDim):

    if self.waiting:
      return
  
    if peakDim.peak.peakList is self.peakList:
      self.updateAfter()

  def selectSpinSystem(self, obj, row, col):

    self.spinSystem = obj[0]
    self.peakPair   = obj[1]

    if self.followSelect.get():
      self.followPeaks()
    
    if self.markSelect.get():
      self.markPeaks()
      
    self.updateButtons()

  def markPeaks(self):

    if self.peakPair:
      peakA, peakB = self.peakPair
      createPeakMark(peakA, lineWidth=1.0)
      createPeakMark(peakB, lineWidth=1.0, remove=False)
 
  def followPeaks(self):

    if self.peakPair and self.windowPane:
      self.windowPane.getWindowFrame()
      zoomToShowPeaks(self.peakPair, self.windowPane)
      

  def getWindows(self):
  
    windows    = []
    if self.peakList:
      project    = self.peakList.root
      spectrum   = self.peakList.dataSource
      tryWindows = getActiveWindows(project)
      for window in tryWindows:
        for windowPane in window.sortedSpectrumWindowPanes():
          if isSpectrumInWindowPane(windowPane, spectrum):
            windows.append( windowPane )
    
    return windows
    
  def getPeakLists(self):
  
    peakLists = []
  
    for experiment in self.nmrProject.experiments:
      refExperiment = experiment.refExperiment
      
      if refExperiment and (refExperiment.name in ALLOWED_REF_EXPTS):
        for spectrum in experiment.dataSources:
          if spectrum.dataType == 'processed':
            for peakList in spectrum.peakLists:
              peakLists.append((self.getPeakListName(peakList), peakList))
  
    peakLists.sort()
    
    return [x[1] for x in peakLists]

  def getPeakListName(self, peakList):
  
    spectrum = peakList.dataSource
  
    return '%s:%s:%d' % (spectrum.experiment.name, spectrum.name, peakList.serial)
    
  def updatePeakLists(self):
  
    index = -1
    names = []
    peakList = self.peakList
    peakLists = self.getPeakLists()
    
    if peakLists:
      names = [self.getPeakListName(pl) for pl in peakLists]
    
      if peakList not in peakLists:
        peakList = peakLists[0]
    
      index = peakLists.index(peakList)
    
    else:
      peakList = None
    
     
    self.changePeakList(peakList)
    self.peakListPulldown.setup(names, peakLists, index)
    
    
  def changePeakList(self, peakList):
      
    if peakList is not self.peakList:
      self.peakList = peakList
      self.updateAfter()
      self.updateButtons()
      self.updateWindowList()
  
  
  def updateJCouplingLists(self, obj=None):

    names = ['<New>',]
    index = 0
    jCouplingLists = self.getJCouplingLists()
    
    if self.jCouplingList not in jCouplingLists:
      self.jCouplingList = None
    
    names.extend(['%d:%s' % (jcl.serial, jcl.name) for jcl in jCouplingLists[1:]])
    index = jCouplingLists.index(self.jCouplingList)
  
    self.jCouplingListPulldown.setup(names, jCouplingLists, index)
  
  
  def getJCouplingLists(self):
    
    jCouplingLists = [None,] + list(self.nmrProject.findAllMeasurementLists(className='JCouplingList'))
  
    return jCouplingLists
  
  
  def changeJCouplingList(self, jCouplingList):
  
    self.jCouplingList = jCouplingList
  
  
  def updateConstraintSets(self, obj=None):

    names = ['<New>',]
    index = 0
    constraintSets = self.getConstraintSets()
    
    if self.constraintSet not in constraintSets:
      self.constraintSet = None
    
    names.extend(['%d' % cs.serial for cs in constraintSets[1:]])
    index = constraintSets.index(self.constraintSet)

    self.constraintSetPulldown.setup(names, constraintSets, index)


  def getConstraintSets(self):
  
    constraintSets = [None,] + self.nmrProject.sortedNmrConstraintStores()

    return constraintSets
  
  
  def changeConstraintSet(self, constraintSet):

    self.constraintSet = constraintSet


  def makeJCouplingList(self, forConstraints=False):
    
    if not self.peakList:
      showError('No peak list', 'No peak list selected.',parent=self)
      return

    experiment = self.peakList.dataSource.experiment
    atomNames  = ATOM_NAME_DICT[experiment.refExperiment.name]
    atomNames  = (atomNames[0], atomNames[3])
  
    if forConstraints or (self.jCouplingList is None):
      name  = '3J[%s]' % (','.join(atomNames))
      self.jCouplingList = self.nmrProject.newJCouplingList(name=name)
  
    jCouplingList = self.jCouplingList
    
    if experiment not in jCouplingList.experiments:
      jCouplingList.addExperiment(experiment)
      
    for spinSystem, peaks, couplingValue, couplingError, phis in self.spinSystemMatrix.objectList:
      residue = spinSystem.residue
      if residue and spinSystem.calcHnHaConfirmed:
        coupling = setResidueJCoupling(jCouplingList, residue, atomNames, couplingValue)
        coupling.error = couplingError
        coupling.peaks = peaks
   
    if not forConstraints:
      self.guiParent.editMeasurements(jCouplingList)
    
    return jCouplingList
   
  def makeCouplingConstraints(self):
  
    if not self.peakList:
      showError('No peak list', 'No peak list selected.',parent=self)
      return

    nmrProject = self.nmrProject
    spectrum   = self.peakList.dataSource
    experiment = spectrum.experiment
    eSerial    = experiment.serial
    
    if not self.jCouplingList:
      msg  = 'No current J Coupling measurement list.'
      msg += ' Continue and make a J Coupling List to derive restraints from?'
      if showOkCancel('Coupling Restraints Warning',msg, parent=self):
        self.makeJCouplingList(forConstraints=True)
      else:
        return
    
    elif experiment not in self.jCouplingList.experiments:
      msg  = 'The current experiment (%s) does not match the current J Coupling List' % experiment.name
      msg += ' Continue and make a new matching J Coupling List?'
      if showOkCancel('Coupling Restraints Warning',msg, parent=self):
        self.makeJCouplingList(forConstraints=True)
      else:
        return

    if not self.jCouplingList.measurements:
      msg  = 'J Coupling measurement list is empty.'
      showWarning('Failure',msg, parent=self)
      return
   
    if self.constraintSet is None:
      self.constraintSet = makeNmrConstraintStore(self.nmrProject)
   
    constraintsPopup = self.guiParent.popups.get('browse_restraints')
    if constraintsPopup:
      constraintsPopup.turnOffNotifiers()

    sSerial    = spectrum.serial
    pSerial    = self.peakList.serial
    mSerial    = self.jCouplingList.serial
    atomNames  = ATOM_NAME_DICT[experiment.refExperiment.name]
    name       = '3J[%s]' % (','.join((atomNames[0], atomNames[3])))
    
    constraintList = self.constraintSet.newJCouplingConstraintList(name=name,unit='Hz',
                                                                   experimentSerials=(eSerial,),
                                                                   measureListSerials=(mSerial,))
     
    for jCoupling in self.jCouplingList.measurements:
      targetValue = jCoupling.value
      error       = jCoupling.error
      upperLimit  = targetValue + error
      lowerLimit  = targetValue - error
      constraint  = constraintList.newJCouplingConstraint(targetValue=targetValue, upperLimit=upperLimit,
                                                          lowerLimit=lowerLimit, error=error)
      
      fixedResonances = [getFixedResonance(self.constraintSet,r) for r in jCoupling.resonances]
      
      item = constraint.newJCouplingConstraintItem(resonances=fixedResonances)
                                                         
      for peak in jCoupling.peaks:
        peakContrib = constraint.newConstraintPeakContrib(experimentSerial=eSerial,
                                                          dataSourceSerial=sSerial,
                                                          peakListSerial=pSerial,
                                                          peakSerial=peak.serial)      
    
    if not constraintsPopup:
      constraintsPopup = self.guiParent.browseConstraints(constraintList)
    
    else:  
      constraintsPopup.update(constraintList)
      constraintsPopup.turnOnNotifiers() 
           
  def makeAngleConstraints(self):
  
    if not self.peakList:
      showError('No peak list', 'No peak list selected.',parent=self)
      return

    if self.constraintSet is None:
      self.constraintSet = makeNmrConstraintStore(self.nmrProject)

    constraintsPopup = self.guiParent.popups.get('browse_restraints')
    if constraintsPopup:
      constraintsPopup.turnOffNotifiers()
     
    def middleAngle(a1, a2):

      delta = abs(a1- a2) % 360
 
      if delta > 180:
        delta = 360 - delta
        return a1 + 0.5*delta
          
      else:
        return a2 + 0.5*delta
   
    nmrProject = self.nmrProject
    spectrum   = self.peakList.dataSource
    experiment = spectrum.experiment
    eSerial    = experiment.serial
    sSerial    = spectrum.serial
    pSerial    = self.peakList.serial
    atomNames  = ATOM_NAME_DICT[experiment.refExperiment.name]
    name       = '3J[%s]' % (','.join((atomNames[0], atomNames[3])))
    error      = self.angleErrorEntry.get() or 0.0
    constraintList = self.constraintSet.newDihedralConstraintList(name=name,unit='degrees',
                                                                    experimentSerials=(eSerial,))
    if self.jCouplingList:
      constraintList.setMeasureListSerials((self.jCouplingList.serial,))
                                                                    
    for spinSystem, peaks, couplingValue, couplingError, phis in self.spinSystemMatrix.objectList:
      
      if not phis:
        continue
    
      residue = spinSystem.residue
      
      if not residue:
        continue
      
      prev = getLinkedResidue(residue, 'prev')
      if not prev:
        continue
      
      atomC0 = prev.findFirstAtom(name='C')
      atomN  = residue.findFirstAtom(name='N')
      atomCa = residue.findFirstAtom(name='CA')
      atomC  = residue.findFirstAtom(name='C')

      resonances = []
    
      for atom in (atomC0,atomN,atomCa,atomC):
        atomSet = atom.atomSet
        
        resonance = None
        if atomSet:
          resonanceSet = atomSet.findFirstResonanceSet()
          
          if resonanceSet:
            resonance = resonanceSet.findFirstResonance()
          else:
            # TBD: Should really use sample info to determine 12C or 13C...
            resonance = nmrProject.newResonance(isotopeCode=ISOTOPE_DICT[atom.name[0]])
            assignAtomsToRes((atomSet,),resonance)
        
        else:
          resonance = nmrProject.newResonance(isotopeCode=ISOTOPE_DICT[atom.name[0]])
     
        resonances.append(resonance)
                  
      fixedResonances = [getFixedResonance(self.constraintSet,r) for r in resonances]
      
      constraint = constraintList.newDihedralConstraint(resonances=fixedResonances,
                                                        origData=couplingValue)
      
      for peak in peaks:
        peakContrib = constraint.newConstraintPeakContrib(experimentSerial=eSerial,
                                                          dataSourceSerial=sSerial,
                                                          peakListSerial=pSerial,
                                                          peakSerial=peak.serial)      
                                                          
      phiRanges = []
      for phi in phis:
        phiMin = phi-error
        phiMax = phi+error
        n = len(phiRanges)
        for i in range(n):
          phiMin2, pphi2, phiMax2 = phiRanges[i]
          
          if ((phiMin >= phiMin2) and (phiMin <= phiMax2)) or \
             ((phiMax >= phiMin2) and (phiMax<=phiMax2)):
            phiMin2 = min(phiMin,phiMin2)
            phiMax2 = max(phiMax,phiMax2)
            phi2    = middleAngle(phiMin2,phiMax2)
            phiRanges[i] = phiMin2, pphi2, phiMax2
       
        else:
          phiRanges.append((phiMin, phi, phiMax))
        
      for phiMin, phi, phiMax in phiRanges:
        constraint.newDihedralConstraintItem(targetValue=phi,
                                             upperLimit=phiMax,
                                             lowerLimit=phiMin)
 
    if constraintsPopup:
      constraintsPopup.update(constraintList)
      constraintsPopup.turnOnNotifiers()  
    else:
      self.guiParent.browseConstraints(constraintList)
            
  def updateButtons(self):
    
    # ['Assign non-root\ndimensions','Make J-Coupling\nList','Make Dihedral\nConstraints']

    buttons = self.tableButtons.buttons
 
    if self.spinSystemMatrix.objectList:
      buttons[0].enable()
      buttons[1].enable()
      buttons[2].enable()
      buttons[3].enable()
    else:
      buttons[0].disable()
      buttons[1].disable()
      buttons[2].disable()
      buttons[3].disable()
 
  def updateAfter(self, obj=None):

    if not self.waiting:
      self.waiting = True
      self.after_idle(self.update)
 
  
  def update(self, recalc=True):
  
    red = '#FF8080'
  
    textMatrix  = []
    objectList  = []
    colorMatrix = []

    if recalc or not self.couplingData:
      self.couplingData = couplingData = self.extractHnHaCouplingData()
    else:
      couplingData = self.couplingData
    
    dataList = []
    for spinSystem, peaks, coupling, couplingError, secShift, phis, ratio, iH, iHa, in  couplingData:
      colors = [None] * 9
      
      if spinSystem.calcHnHaConfirmed:
        useText   = 'Yes'
      else:
        colors[0] = red
        colors[1] = red
        colors[2] = red
        colors[3] = red
        useText   = 'No'

      #if secShift > 0.1:
      #  colors[4] = '#C080F0'
      #elif secShift < -0.2:
      #  colors[4] = '#208020'
      
      residue = spinSystem.residue
      if residue:
        chain = residue.chain
        chainCode = ''
        
        if len(chain.molSystem.chains) > 1:
          chainCode = chain.code + ' '
      
        spinSystText = '%s%d%s' % (chainCode, residue.seqCode, residue.ccpCode)
        key = '%s%6.6d%s' % (chainCode, residue.seqCode, residue.ccpCode) 

      elif spinSystem.ccpCode:
        spinSystText = '%s{%d}' % (spinSystem.ccpCode,spinSystem.serial)
        key = '%s{%6.6d}' % (spinSystem.ccpCode,spinSystem.serial)

      elif spinSystem.name:
        spinSystText = '{%d:%s}' % (spinSystem.serial,spinSystem.name)
        key = '{%6.6d:%s}' % (spinSystem.ccpCode,spinSystem.serial)
        
      else:
        spinSystText = '{%d}' % spinSystem.serial
        key = '{%6.6d}' % spinSystem.serial
        
      datum = [spinSystText, useText,
               coupling, couplingError,
               secShift, ','.join(['%.0f' % phi for phi in phis]),
               ratio, iH, iHa]
    
      dataList.append((key, datum,(spinSystem, peaks, coupling, couplingError, phis) , colors))
  
    dataList.sort()
    for key, datum, objs, colors in dataList:
      objectList.append(objs)
      textMatrix.append(datum)
      colorMatrix.append(colors)
 
    self.spinSystemMatrix.update(textMatrix=textMatrix,
                                 objectList=objectList,
                                 colorMatrix=colorMatrix)
    self.updateButtons()
    self.waiting = False
  

  def assignNonRootDim(self):
  
    if not self.peakList:
      showError('No peak list', 'No peak list selected.',parent=self)
      return

    spectrum  = self.peakList.dataSource
    boundDims = getOnebondDataDims(spectrum)
    dataDims  = spectrum.sortedDataDims()
    atomNames = ATOM_NAME_DICT[spectrum.experiment.refExperiment.name]
    atomNames = (atomNames[0], atomNames[3])
    nmrProject = self.nmrProject
   
    for dimA, dimB in boundDims:
      dataDims.remove(dimA)
      dataDims.remove(dimB)
  
    dim = dataDims[0].dim
    
    msg = 'Continue and assign dimension %d %s resonances?'
    if not showOkCancel('Confirm', msg % (dim,','.join(atomNames)), parent=self):
      return
  
    for spinSystem, peaks, coupling, couplingError, phis in self.spinSystemMatrix.objectList:
      
      residue    = spinSystem.residue
      resonanceA = None
      resonanceB = None
      
      for resonance in spinSystem.resonances:
        resonanceSet = resonance.resonanceSet
        
        if resonanceSet:
          assignNames = [ass.name for ass in resonanceSet.atomSets]
        else:
          assignNames = resonance.assignNames
      
        if atomNames[0] in assignNames:
          resonanceA = resonance
        if atomNames[1] in assignNames:
          resonanceB = resonance
  
      if not resonanceA:
        resonanceA = nmrProject.newResonance(isotopeCode=ISOTOPE_DICT[atomNames[0][0]])

        if residue:
          for atom in residue.atoms:
            atomSet = atom.atomSet
            
            if atomSet and (atomSet.name == atomNames[0]):
              assignAtomsToRes((atomSet,),resonanceA)
              break

      if not resonanceA:
        resonanceB = nmrProject.newResonance(isotopeCode=ISOTOPE_DICT[atomNames[1][0]])
        
        # Could change here for through carbonyl couplings - different residue
 
        if residue:
          for atom in residue.atoms:
            atomSet = atom.atomSet
            
            if atomSet and (atomSet.name == atomNames[1]):
              assignAtomsToRes((atomSet,),resonanceB)
              break
      
      peakDimA = peaks[0].findFirstPeakDim(dim=dim)
      peakDimB = peaks[1].findFirstPeakDim(dim=dim)
      
      assignResToDim(peakDimA,resonanceA)
      assignResToDim(peakDimB,resonanceB)
                
  def extractHnHaCouplingData(self):
  
    if not self.peakList:
      return []
  
    kA, kB, kC = self.getCoefficients()
    tTime      = self.timeEntry.get()
    peakList   = self.peakList
    iType      = self.intensityPulldown.getText()
    spectrum   = self.peakList.dataSource

    # NOTE:ED - not sure whether to just exit or carry on with an empty list
    getBonds = getOnebondDataDims(spectrum)
    if getBonds and len(getBonds) > 0:
      dims = [dataDim.dim for dataDim in getBonds[0]]
    else:
      dims = []
      msg = 'Contains no one-bond dataDims - possibly wrong experimentType.'
      showWarning('HnHa Coupling error', msg, parent=self)

    project    = self.project
    atomNames  = ATOM_NAME_DICT[spectrum.experiment.refExperiment.name]
    atomNames  = (atomNames[0], atomNames[3])
    invalid    = 0
    noise      = getSpectrumNoise(spectrum)
    useGly     = self.glycineSelect.get()
    aFilter    = self.angleFilterSelect.get()
    relaxCorr  = self.relaxCorrEntry.get() or 1.11
    invalidReasons = []

    f = 4*twoPi*twoPi*tTime*tTime
    
    degrees = 360/twoPi
    
    for dataDim in spectrum.dataDims:
      if dataDim.dim not in dims:
        nonRootDim = dataDim.dim
        break
    
    ssPeaks = {}
    for peak in peakList.peaks:
      intensity = peak.findFirstPeakIntensity(intensityType=iType)
      
      height = None
      if iType != 'height':
        heightObj = peak.findFirstPeakIntensity(intensityType='height')
        if heightObj:
          height = heightObj.value
        
      else:
        height = intensity.value
            
      if (not intensity) or (not intensity.value) or not height:
        if not intensity:
          invalidReasons.append('HNHA Error: no intensity [peak: %s]' % str(peak))
        elif not intensity.value:
          invalidReasons.append('HNHA Error: no intensity.value [peak: %s]' % str(peak))
        elif not height:
          invalidReasons.append('HNHA Error: no height [peak: %s]' % str(peak))
        invalid += 1
        continue
        
        
      spinSystem = None
      for dim in dims:
        peakDim = peak.findFirstPeakDim(dim=dim)
        
        for contrib in peakDim.peakDimContribs:
          spinSystem0 = contrib.resonance.resonanceGroup
          
          if spinSystem0:
            spinSystem = spinSystem0
            break
        
        else:
          continue
          
        break      

      if spinSystem:
        ppm = peak.findFirstPeakDim(dim=nonRootDim).value
        
        if not ssPeaks.has_key(spinSystem):
          ssPeaks[spinSystem] = []
        
        ssPeaks[spinSystem].append((peak, intensity.value, height, ppm))
    
    data = []
    for spinSystem in ssPeaks:
      peaks = ssPeaks[spinSystem]
      
      residue = spinSystem.residue
      context = [None] * 5
      
      if residue:
        ccpCode = residue.ccpCode
        molType = residue.molResidue.molType
        prev1   = getLinkedResidue(residue,'prev')
        next1   = getLinkedResidue(residue,'next')
        
        if next1:
          next2 = getLinkedResidue(next1,'next')
        else:
          next2 = None
          
        if prev1:
          prev2 = getLinkedResidue(prev1,'prev')
        else:
          prev2 = None
        
        context = [prev2,prev1,residue,next1,next2]
        
      elif spinSystem.ccpCode:
        ccpCode = spinSystem.ccpCode
        molType = spinSystem.molType
      else:
        ccpCode = 'Gln' # TBD: Use generic residue type
        molType = 'protein'
        
      if ccpCode in ('Pro','Hyp'):
        continue

      elif (ccpCode == 'Gly') and not useGly:
        continue
      
      chemComp = project.findFirstChemComp(molType=molType, ccpCode=ccpCode)
      atomName = atomNames[1]
      if ccpCode == 'Gly' and atomName == 'HA':
        chemAtom  = chemComp.findFirstChemAtom(name='HA2') 
      else:
        chemAtom  = chemComp.findFirstChemAtom(name=atomNames[1]) 
      
      coilShift = getRandomCoilShift(chemAtom, context=context)
        
      if not hasattr(spinSystem, 'calcHnHaConfirmed'):
        if ccpCode == 'Gly':
          spinSystem.calcHnHaConfirmed = False
        else:
          spinSystem.calcHnHaConfirmed = True
        
      if len(peaks) > 2:
        peaks2 = set()
        
        for atomName in atomNames:
          best = None
        
          for peak, intensity, height, ppm in peaks:
            probability = self.failsafeGetAtomProbability(ccpCode, atomName, ppm)
            # wb104: 10 Aug 2014: added abs() because intensities can be negative
            score = abs(intensity) * probability
          
            if (best is None) or (score > best[0]):
              best = (score, peak, intensity, height,  ppm)
          
          peaks2.add(best[1:])
          
        peaks = list(peaks2)

      if len(peaks) != 2:
        continue
      
      peakA, intensityA, heightA, ppmA = peaks[0]
      peakB, intensityB, heightB, ppmB = peaks[1]
      
      score1 = self.failsafeGetAtomProbability(ccpCode, atomNames[0], ppmA)
      score2 = self.failsafeGetAtomProbability(ccpCode, atomNames[0], ppmB)
     
      if score2 > score1:
        peakA, intensityA, heightA, ppmA = peaks[1]
        peakB, intensityB, heightB, ppmB = peaks[0]
      
      ratio = intensityB/intensityA
      
      if ratio > 0:
        invalidReasons.append('HNHA Error: bad ratio [value: %f spin-system: %s]' %(ratio, str(spinSystem)))
        invalid += 1
        continue
      
      coupling = relaxCorr * atan(sqrt(-ratio )) / (twoPi*tTime)
      
      sigmaSq = relaxCorr * noise * noise * (1+(ratio*ratio)) / ( -1 * f * ratio * heightA * heightA * (1-ratio) * (1-ratio) )
      
      # coupling to phis
      
      a = kB*kB
      b = 4*kA*(kC-coupling)
      
      if b > a:
        invalidReasons.append('HNHA Error: b > a [a: %f b: %f spin-system: %s]' % (a,b, str(spinSystem)))
        invalid += 1
        continue
      
      x = sqrt(a-b)
      c1 = (kB + x)/(2*kA)
      c2 = (kB - x)/(2*kA)
      
      phis = []
      
      if abs(c1) <= 1.0:
        acosC1 = acos(c1)
        phis.append(acosC1)
        phis.append(-acosC1)

      if abs(c2) <= 1.0:
        acosC2 = acos(c2)
        phis.append(acosC2)
        phis.append(-acosC2)
      
      secondaryShift = ppmB-coilShift
      
      phis = [((phi*degrees+60)%360)-180 for phi in phis]
      
      if aFilter:
        for phi in phis[::-1]:
          if (phi > 0):# and (secondaryShift < -0.2 or secondaryShift > 0.1):
            phis.pop()
          elif (secondaryShift < -0.2) and (phi < -120):
            phis.pop()
      
      phis.sort()
      
      datum = (spinSystem, (peakA, peakB),
               coupling, sqrt(sigmaSq),
               secondaryShift, phis, 
               ratio, intensityA, intensityB)
      
      data.append(datum)   
     
    if invalid:

      if invalidReasons:
        msg = 'reasons for HNHA errors'
        msg2 = '-'*len(msg)
        invalidReasons.insert(0, msg2)
        invalidReasons.insert(0, msg)

      if self.inBody:
        msg1 = 'Warning: %d peaks have an unusable %s value' % (invalid,iType)
        invalidReasons.insert(0, msg1)
        error = '\n'.join(invalidReasons)
        print(error)

      else:
        msg1 = '%d peaks have an unusable %s value' % (invalid,iType)
        invalidReasons.insert(0, msg1)
        error = '\n'.join(invalidReasons)
        showWarning('Warning',
                    error,
                    parent=self)
       
    return data

  
  def failsafeGetAtomProbability(self, ccpCode, atomName, ppm):
  
    probability = getAtomProbability(ccpCode, atomName, ppm)
    if probability is None:
      if atomName in FALLBACK_REF_SHIFTS:
        mean, sd = FALLBACK_REF_SHIFTS[atomName]
        e = (ppm-mean)/sd
        probability = exp(-0.5*e*e)/(sd*2.506628274631)
      else:
        probability = 0.0
      
    return probability
  
  def destroy(self): 

    self.administerNotifiers(self.unregisterNotify)
   
    BasePopup.destroy(self)

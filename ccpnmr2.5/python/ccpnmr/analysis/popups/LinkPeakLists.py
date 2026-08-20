
"""
======================COPYRIGHT/LICENSE START==========================

LinkPeakLists.py: Part of the CcpNmr Analysis program

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

import cPickle

from ccpnmr.analysis.core.AssignmentAdvanced import pickAssignSpecFromRoot, assignSpecNonRootResonances

from ccpnmr.analysis.popups.BasePopup     import BasePopup
from ccpnmr.analysis.core.ExperimentBasic import getOnebondExpDimRefs, getPrimaryDataDimRef, getSeqAssignRefExperiments
from ccpnmr.analysis.core.MarkBasic       import createPeakMark
from ccpnmr.analysis.core.WindowBasic import getWindowPaneName, getSpectrumViews, toggleSpectrum
from ccpnmr.analysis.core.WindowBasic import findOrthogonalWindows, getPeakWindowPosition
from ccpnmr.analysis.core.WindowBasic import getPeakDimAxisMapping
from ccpnmr.analysis.core.Util import getAnalysisDataDim

from memops.gui.ButtonList      import ButtonList, UtilityButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.IntEntry        import IntEntry
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Frame           import Frame
from memops.gui.Label           import Label
from memops.gui.LabelDivider    import LabelDivider
from memops.gui.MessageReporter import showWarning, showOkCancel
from memops.gui.ProgressBar     import ProgressBar
from memops.gui.PulldownList    import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.TabbedFrame     import TabbedFrame


# Shrink non-essential bits of window
# Equalise window sizes
# Equalise window ranges
# 3D roots

defaultTolerances = {'1H':0.05,'13C':0.20,'15N':0.3}

def testPopup(argServer):
  

  popup = LinkPeakListsPopup(argServer.parent)
  popup.open()


class LinkPeakListsPopup(BasePopup):
  """
  **Pick and Assign Peaks Based on Locations in Root Spectra**
  
  This popup window is designed to help in the early stages of peak picking and
  resonance assignment. The main idea is that certain "root" spectra, typically
  an HNCO or 15N HSQC, with initialised peak lists are used as a basis for
  locating and assigning peaks in other spectra that overlap in the relevant
  dimensions. For example you may use the amide peak positions in a 15N HSQC
  spectrum to pick peaks in restricted regions of a 3D HNCACB spectrum and copy
  amide resonance assignments to the amide dimensions of the 3D spectrum where
  they overlap with the HSQC peaks.

  Often the root spectrum, that is the source of peak position and assignment
  information, will only be assigned in an anonymous way; the peaks will have 
  spin system and resonance numbers (e.g. "{7}[4][5]") and not link to any
  particular residues or atoms. For this tool it doesn't matter whether the root
  peaks carry full assignments or not, all that is important is that there are
  some peaks to specify locations and some form of assignment to copy. The
  easiest way to setup initial assignments on a root spectrum, after picking
  some peaks, is to use the `Initialise Root Resonances`_ option.

  In normal operation the user chooses a peak list to act as the the source of
  the root assignments and positions, and maybe a spectrum window to view those
  peaks in. Next the assignable spectra that are to be the target for peak
  picking and/or assignment are chosen. Firstly, the user adds and removes
  spectrum windows to the "Target Windows" list, the spectra to be operated on
  are then selected from those that are visible in these windows. By selecting
  different kinds of window different kinds of spectra may be operated on, so
  for example the user could both work with HCN spectra and HHN spectra with the
  same HSQC root. This tool is largely visual and it is important to be able to
  see the locations that are bing considered (e.g. amide), hence the user is
  only presented with assignable spectra that can potentially be seen. The 
  "Assignable Spectra" section is filled with the spectra for the target
  windows that may be assigned, and the user double-clicks in the "Active?"
  to set whether the individual spectra should be operated on or not (for
  peak picking and assignment).

  One windows and spectra are setup is is advisable to consider the "Tolerances"
  tab that controls how wide a search region is used to pick and assign peaks
  relative to the positions of the reference, root peaks. The user can also
  setup exclusions to avoid picking peaks near the water signal (e.g. in 15N
  HSQC-NOESY) or homonuclear diagonal.

  The last "Link Peaks" tab is the one that remains active while the user is
  actually running the the peak picking and assignment functions. Here, the main
  table lists all of the peaks in the root peak list that are used for position
  and assignment references. Clicking on a row in this table, assuming the
  relevant navigations are checked above, will cause the location of any root
  and target spectrum windows to move in order to show the root (e.g. amide)
  location for the selected row. Even if peaks are not peaked or assigned this 
  tool may be used for efficient coordination in window. Typically clicking  on
  an HSQC peak will present an amide position in 3D target windows, locating the
  X and Z axes, so that the user can see the column of peaks that may be
  picked/assigned.

  For a single selected root location, the user may transfer assignments or pick
  peaks *and* transfer assignments, by clicking the appropriate buttons in the
  "Pick & Assign Functions" section. This is one-by-one way of working is the
  safest because the user is presented will the spectra for each location and
  will be able to view the result of the peak picking and assignment. The
  equivalent "Process-all Functions" will work though all root locations in the
  table picking and assigning peaks en masse, according to the set tolerances, 
  in a quick but less controlled manner.
  
  The "Assign Non-root Resonances" options are present so that you can give a
  starting assignment to the assigned spectra in the dimension that does not
  match the root location. For example, you could add 13C resonances to an
  HNcoCA spectrum or 1H resonances to the indirect dimension 1H dimension of a
  15N HSQC-TOCSY. These "non-root" resonance numbers will all be new and unique,
  thus this operation should only be used for types of experiment where the is
  one peak for each non-root resonance. For example HNcoCA has one peak per CA
  resonance by HNCA usually has two, so the function is only recommended for the
  former.

  **Caveats & Tips**
  
  This tool can be operated in two slightly different ways, according to the
  preference of the user and the quality of the spectra. The "Process-all"
  functions can be used to start with, given fairly strict tolerances, but the
  user should then go through each root position checking for and tidying up
  mistakes (picking noise & artefacts for example). Alternatively the root
  locations could be picked and assigned one-by-one so the user can spot
  problems as they occur.

  Although the same tolerances, set via the "Tolerances" tab, are used for both
  peak picking and for resonance assignment some peaks that can be picked may
  not be assigned with the same settings. Whereas the peak picking is done
  relative to the root location the assignment, in common with the rest if
  Analysis, is relative to the current chemical shift average and thus may
  differ from the root location. Widening the tolarances a little, or increasing
  the chemical shift weighting of the root spectrum can allow assignments to be
  made if they were previously out of bounds.

  .. _`Initialise Root Resonances`: InitRootAssignmentsPopup.html
  
  """

  def __init__(self, parent, *args, **kw):

    self.waiting    = False
    self.windowPane = None
    self.windowPanes = []
    self.rootPeak   = None
    self.peakList   = None
    self.rootPane = None
    self.guiParent  = parent
    self.targetPeakLists = []
    self.dimMappings     = {}
    self.selDimMapping   = {}
    self.nonRootDict     = {}
    self.marks           = []
    self.project         = parent.project
    
    BasePopup.__init__(self, parent=parent, title='Assignment : Pick & Assign From Roots', **kw)

  def open(self):
  
    BasePopup.open(self)
    self.updateRootPeakList()


  def close(self):
  
    self.setAppDataOptions()
    BasePopup.close(self)


  def body(self, guiFrame):

    self.geometry('500x600')
    
    guiFrame.expandGrid(0,0)
    
    self.progressBar = ProgressBar(self,text = '', progress = 0, total = 100, title='Progress')
    self.progressBar.close()
    
    row = 0
    tipTexts = ['Selection of source or "root" for assignments, the target spectra and which spectrum windows to navigate within',
                'Settings that relate to assignment and peak picking tolerances when going from "root" positions to related spectra',
                'The main assignment & peak picking functions; using a table of the "root" peak that are the source of peak (e.g. amide) positions and assignments.']
    options = ['Windows & Spectra','Tolerances','Link Peaks']
    tabbedFrame = TabbedFrame(guiFrame, options=options,
                              grid=(row,0), tipTexts=tipTexts)
    frameA, frameB, frameC = tabbedFrame.frames
    frameA.expandGrid(5,1)
    frameB.expandGrid(3,0)
    frameC.expandGrid(2,0)

    # Windows & Spectra
    
    frame = Frame(frameA, grid=(0,0), gridSpan=(1,2), sticky='ew')
    frame.grid_columnconfigure(3, weight=1)

    label = Label(frame, text='Root Peak List:', grid=(0,0))
    tipText = 'Selects which peak list is considered as the "root"; the positions and assignments that will be used to pick/assign related spectra'
    self.peakListPulldown = PulldownList(frame, callback=self.changeRootPeakList,
                                         grid=(0,1), tipText=tipText)

    label = Label(frame, text='Root Window:', grid=(0,2))
    tipText = 'Selects which spectrum window is used to navigate to the location of root peak selected in the "Link Peaks" table'
    self.rootPanePulldown = PulldownList(frame, callback=self.changeRootWindow,
                                         grid=(0,3), tipText=tipText)
   
    div = LabelDivider(frameA, text='Target Windows',
                       grid=(1,0), gridSpan=(1,2))
    
    tipTexts = ['Remove the selected target window from consideration, so that it is no longer used for navigation and assignment',
                'Add the window selected in the adjacent pulldown list as a source of assignable spectra and a target for navigation']
    texts    = ['Remove Target Window','Add Target Window:']
    commands = [self.removeWindow,self.addWindow]
    self.windowButtons = ButtonList(frameA, texts=texts, tipTexts=tipTexts,
                                    commands=commands, grid=(2,0))

    tipText = 'Selects a spectrum window, from those not already selected, that may be used for assignment/navigation'
    self.windowPulldown = PulldownList(frameA, callback=None,
                                       grid=(2,1), tipText=tipText)

    tipTexts = ['The serial number of the spectrum window',
                'The name of the  spectrum window used for navigation and providing assignable spectra',
                'When there are multiple options, states which axes of the window (X, Y, Z...) correspond to the directly bound "root" dimensions']
    headingList = ['#','Name','Selected\nMapping']
  
    self.mappingPulldown = PulldownList(self, callback=self.setMapping)
    editWidgets      = [None, None, self.mappingPulldown]
    editGetCallbacks = [None, None, self.getMapping]
    editSetCallbacks = [None, None, self.setMapping]
    self.windowMatrix = ScrolledMatrix(frameA, headingList=headingList,
                                       editSetCallbacks=editSetCallbacks,
                                       editGetCallbacks=editGetCallbacks,
                                       editWidgets=editWidgets, 
                                       callback=self.selectWindow,
                                       grid=(3,0), gridSpan=(1,2),
                                       tipTexts=tipTexts)
   
    div = LabelDivider(frameA, text='Assignable Spectra', grid=(4,0), gridSpan=(1,2))

    tipTexts = ['The "experiment:spectrum" name for the spectrum which may be used for peak picking & assignment; must be present in the above windows',
                'Whether the spectrum is considered active for the processes of peak picking and/or assignment; inactive ones will not be affected',
                'Whether spectrum dimensions that do not map to the root (typically non-amide) may be assigned; common for HNCO, HNH-TOCSY etc.',
                'The full CCPN experiment type of the spectrum; should be something with an obvious relation to the root peak list',
                'The kinds of isotope (or otherwise) present on the dimensions of the spectrum']
    headingList = ['Spectrum','Active?','Assign\nNon-root dim?',
                   'Experiment\nType','Dimensions']
    editWidgets      = [None, None, None, None, None,]
    editGetCallbacks = [None, self.togglePeakList, self.toggleNonRootAssign, None, None,]
    editSetCallbacks = [None, None, None, None, None,]
    self.spectrumMatrix = ScrolledMatrix(frameA, headingList=headingList,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets, callback=None,
                                         grid=(5,0), gridSpan=(1,2),
                                         tipTexts=tipTexts)
    
    # Tolerances
    


    #div = LabelDivider(frameB, text='Tolerances')
    #div.grid(row=0, column=0, sticky='nsew')
    
    frame = Frame(frameB, grid=(1,0), sticky='nsew')
    frame.grid_columnconfigure(5, weight=1)
    
    self.tolLabel1 = Label(frame, text='Root Dim 1:', grid=(0,0))
    tipText = 'Sets the upper limit to the dimension 1 assignment tolerance, within which peaks may be picked and assignments made, relative to the root location'
    self.tolEntry1 = FloatEntry(frame, returnCallback=self.setDataDimTolerances,
                                width=8, grid=(0,1), tipText=tipText)
    self.tolEntry1.bind('<Leave>', self.setDataDimTolerances, '+')

    self.tolLabel2 = Label(frame, text='Root Dim 2:', grid=(1,0))
    tipText = 'Sets the upper limit to the dimension 2 assignment tolerance, within which peaks may be picked and assignments made, relative to the root location'
    self.tolEntry2 = FloatEntry(frame, returnCallback=self.setDataDimTolerances,
                                width=8, grid=(1,1), tipText=tipText)
    self.tolEntry2.bind('<Leave>', self.setDataDimTolerances, '+')
    
    label = Label(frame, text='Min Water:', grid=(2,0))
    tipText = 'Sets the lower bound to the 1H ppm exclusion zone, usually representing a water signal, within which no peaks will be picked or assigned'
    self.minWaterEntry = FloatEntry(frame, width=8, text=4.95, 
                                    grid=(2,1), tipText=tipText)

    label = Label(frame, text='Max Water:', grid=(3,0))
    tipText = 'Sets the upper bound to the 1H ppm exclusion zone, usually representing a water signal, within which no peaks will be picked or assigned'
    self.maxWaterEntry = FloatEntry(frame, width=8, text=4.95, 
                                    grid=(3,1), tipText=tipText)
     
    label = Label(frame, text='Max Num Marks:', grid=(4,0))
    tipText = 'When using multi-dimensional cross marks to indicate peak positions, sets how many marks persist from subsequent peak navigations'
    self.marksEntry = IntEntry(frame, width=8, text=1, 
                               grid=(4,1), tipText=tipText)
    self.marksEntry.grid(row=4, column=1, sticky='w')
     
    label = Label(frame, text='1H Diagonal:', grid=(5,0))
    tipText = 'Sets the width of the exclusion zone around the 1H-1H homonuclear diagonal, within which no peaks will be picked or assigned '
    self.diagEntry = FloatEntry(frame, width=8, text=0.25,
                                grid=(5,1), tipText=tipText)
 
    # Peaks

    div = LabelDivider(frameC, text='Root Peaks', grid=(0,0))

    frame = Frame(frameC, grid=(1,0), sticky='ew')
    frame.expandGrid(1,3)
    
    #label = Label(frame, text='Find Spin System:')
    #label.grid(row=1, column=0, sticky='nw')

    label = Label(frame, text='Navigate to root', grid=(0,0))
    tipText = 'Sets whether clicking in the root peak table will cause the selected root spectrum window to display the selected peak'
    self.followRootSelect = CheckButton(frame, callback=None, tipText=tipText,
                                        grid=(0,1), selected=True)

    label = Label(frame, text='Navigate to targets', grid=(0,2))
    tipText = 'Sets whether clicking in the root peak table will cause the selected navigation window views to move to the selected root position'
    self.followTargetSelect = CheckButton(frame, callback=None, tipText=tipText,
                                          grid=(0,3), selected=True)

    label = Label(frame, text='Set Non-root Atom Types', grid=(0,4))
    tipText = 'Sets whether for appropriate experiments like HNCO, HNCA, HNHA etc., whether non-root assignment set resonance atom type'
    self.assignTypeSelect = CheckButton(frame, callback=None, tipText=tipText,
                                        grid=(0,5), selected=False)

    tipTexts = ['The serial number of the peak in the root peak list',
                'The assignment of the root peak in the first root dimension, that maps to the (typically) higher dimensionality assignment spectra',
                'The assignment of the root peak in the second root dimension, that maps to the (typically) higher dimensionality assignment spectra',
                'The location of the root peak in the first root dimension; the basis for peak picking and assignment zones',
                'The location of the root peak in the second root dimension; the basis for peak picking and assignment zones']
    self.rootPeakTipTexts = tipTexts
    headingList = ['#','Assign F1','Assign F2','Shift F1','Shift F2']
    self.rootPeakMatrix = ScrolledMatrix(frameC, headingList=headingList,
                                         callback=self.selectRootPeak,
                                         grid=(2,0), tipTexts=tipTexts)

    tipTexts = ['If a peak is selected, select the next peak in the table and navigate to that peaks root location',
                'If a peak is selected, select the previous peak in the table and navigate to that peaks root location']
    texts    = ['Next Root','Previous Root']
    commands = [self.nextTarget,self.prevTarget]
    self.prevNextButtons = ButtonList(frameC, texts=texts, tipTexts=tipTexts,
                                      commands=commands, grid=(3,0))

    div = LabelDivider(frameC, text='Pick & Assign Functions', grid=(4,0))

    tipTexts = ['Using the selected root peak as the source of assignment, spread assignments to the active spectra, within stated tolerances',
                'Using the selected root peak as the source of assignment and centre of peak pick zone, pick and assign peaks in the active spectra',
                'For the spectrum positions that match the selected root peak, add new (separate) non-root resonances to peaks (where missing)']
    texts    = ['Assign\nRoot Resonances',
                'Pick & Assign\nRoot Resonances',
                'Assign Non-root\nResonances']
    commands = [self.assignTarget,
                self.pickAssignTarget,
                self.assignNonRootTarget]
    self.rootPeakButtons = ButtonList(frameC, texts=texts, tipTexts=tipTexts,
                                      commands=commands, grid=(5,0))

    div = LabelDivider(frameC, text='Process-all Functions', grid=(6,0))
    
    tipTexts = ['Using all peaks in the root peak list as an assignment source, spread assignments to the active spectra, within stated tolerances',
                'Using all peaks in the root peak list as assignment sources and peak pick locations, pick and assign peaks in the active spectra',
                'For the spectrum positions that match the all root peaks, add new (separate) non-root resonances to peaks (where missing)']
    texts = ['Assign All\nRoot Resonances',
             'Pick All & Assign\nRoot Resonances',
             'Assign All Non-root\n Resonances']
    commands = [self.assignAllTarget,
                self.pickAssignAllTarget,
                self.assignNonRootAllTarget]
    self.bottomButtons = ButtonList(frameC, commands=commands, texts=texts,
                                    grid=(7,0), sticky='ew', tipTexts=tipTexts)

    buttons = UtilityButtonList(tabbedFrame.sideFrame, helpUrl=self.help_url,
                                grid=(0,0), sticky='e')
    
      
    self.getAppDataOptions()  
      
    self.updateRootPeakList()
    self.updateRootWindow()
    self.updateWindowPulldown()

    self.administerNotifiers(self.registerNotify)

  def administerNotifiers(self, notifyFunc):

  
    for func in ('__init__', 'delete','setName'):
      for clazz in ('ccp.nmr.Nmr.DataSource', 'ccp.nmr.Nmr.Experiment'):
        notifyFunc(self.updateSpectra, clazz, func)
      
      notifyFunc(self.updateWindowPaneAfter, 'ccpnmr.Analysis.SpectrumWindowPane', func)
    notifyFunc(self.updateWindowsAfter, 'ccpnmr.Analysis.SpectrumWindow', 'setName')
      
    notifyFunc(self.updateSpectra, 'ccpnmr.Analysis.SpectrumWindowView', 'setIsPosVisible')
    notifyFunc(self.updateSpectra, 'ccpnmr.Analysis.SpectrumWindowView', 'setIsNegVisible')

    for func in ('__init__', 'delete',):
      notifyFunc(self.updateAfter,   'ccp.nmr.Nmr.Peak',     func)
      notifyFunc(self.updateSpectra, 'ccp.nmr.Nmr.PeakList', func)

    for func in ('setAnnotation','setPosition','setNumAliasing'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.PeakDim', func)


  def toggleNonRootAssign(self, peakList):
        
    boolean = self.nonRootDict.get(self.getPeakListId(peakList), False)
    self.nonRootDict[self.getPeakListId(peakList)] = not boolean
    
    self.updateSpectra()

  def setMapping(self, null):
    
    
    mapping = self.mappingPulldown.getObject()
    if self.selDimMapping.get(self.windowPane) != mapping:
      self.selDimMapping[self.windowPane] = mapping
      self.updateWindows()
  
  
  def getMapping(self, pane):

    index = -1
  
    mapping  = self.selDimMapping.get(pane)
    mappings = self.dimMappings.get(pane) or []

    if mapping and (mapping in mappings):
      index = mappings.index(mapping)
     
    self.mappingPulldown.setup(mappings, mappings, index)

  def updateButtons(self):

    if self.windowPane:
      self.windowButtons.buttons[0].enable()
      
    else:
      self.windowButtons.buttons[0].disable()

    if self.peakList and self.peakList.peaks and self.targetPeakLists:
      self.bottomButtons.buttons[0].enable()
      self.bottomButtons.buttons[1].enable()
      self.bottomButtons.buttons[2].enable()
      
    else:  
      self.bottomButtons.buttons[0].disable()
      self.bottomButtons.buttons[1].disable()
      self.bottomButtons.buttons[2].disable()
   
    if self.rootPeak:
      self.prevNextButtons.buttons[0].enable()
      self.prevNextButtons.buttons[1].enable()

    else:  
      self.prevNextButtons.buttons[0].disable()
      self.prevNextButtons.buttons[1].disable()
    
   
    if self.rootPeak and self.targetPeakLists:
      self.rootPeakButtons.buttons[0].enable()
      self.rootPeakButtons.buttons[1].enable()
      self.rootPeakButtons.buttons[2].enable()

    else:  
      self.rootPeakButtons.buttons[0].disable()
      self.rootPeakButtons.buttons[1].disable()
      self.rootPeakButtons.buttons[2].disable()


  def getRootPeakLists(self):
  
    project = self.project
    
    peakLists = []
    for experiment in self.nmrProject.experiments:
      isRoot = False
      for expTransfer in experiment.expTransfers:
        if expTransfer.transferType in ('onebond','CP'):
          expDimRefs = expTransfer.sortedExpDimRefs()
          if '1H' in expDimRefs[0].isotopeCodes:
            if '1H' not in expDimRefs[1].isotopeCodes:
              isRoot = True
              break
              
          else:
            isRoot = True
            break
      
      if isRoot:
        for spectrum in experiment.dataSources:
          if len(spectrum.dataDims) != 2:
            if experiment.refExperiment:
	      if experiment.refExperiment.name not in ('H[N[CO]]','H[N[co[CA]]]'):
	        continue
	    else:
              continue
        
          for peakList in spectrum.peakLists:
	    data = (experiment.name,spectrum.name,peakList.serial)
            peakLists.append(['%s:%s:%d' % data, peakList])
    
    peakLists.sort()
    return peakLists


  def getRootDataDims(self):
  
    dims = []
    if self.peakList:
      spectrum = self.peakList.dataSource
      experiment = spectrum.experiment
      expDimRefPairs = []
      
      for expDimRef0, expDimRef1 in getOnebondExpDimRefs(experiment):
        if '1H' in expDimRef0.isotopeCodes:
          expDimRefPairs.append( (expDimRef0, expDimRef1) ) 
        elif '1H' in expDimRef1.isotopeCodes:
          expDimRefPairs.append( (expDimRef0, expDimRef1) ) 
          
      if not expDimRefPairs:
        expDimRefPairs = getOnebondExpDimRefs(experiment)
        
      for expDimRef0, expDimRef1 in expDimRefPairs:
        dataDim0 = spectrum.findFirstDataDim(expDim=expDimRef0.expDim)
        dataDim1 = spectrum.findFirstDataDim(expDim=expDimRef1.expDim)

        if dataDim0 and dataDim1:
          dims = [(dataDim0.dim, dataDim0), (dataDim1.dim, dataDim1)]
          dims.sort()
          dims = [x[1] for x in dims]
          break # First onebond pair only

    return dims


  def getTargetPeakLists(self):
  
    peakLists = set()
    for windowPane in self.windowPanes:
      for view in windowPane.spectrumWindowViews:
        spectrum = view.analysisSpectrum.dataSource
        freqDims = [dd for dd in spectrum.dataDims if dd.className == 'FreqDataDim']
        
        if len(freqDims) >= 2:
          peakList = spectrum.activePeakList or spectrum.findFirstPeakList()
          
          if peakList:
            peakLists.add(peakList)
         
    return list(peakLists)


  def getRootWindows(self):
  
    panes = set()
    if self.peakList:
      views = getSpectrumViews(self.peakList.dataSource)
      for view in views:
        panes.add( view.spectrumWindowPane )
    
    winData = [[getWindowPaneName(p), p] for p in panes]
    winData.sort()
    
    return winData


  def getTargetWindows(self):
      
    self.dimMappings = {}
    
    if self.rootPane:
      windowZplanes = findOrthogonalWindows(self.rootPane, [], minDims=2)
      for (planeName0, pane, positions) in windowZplanes:
        planeName = planeName0.split(' in ')[0]
        if self.dimMappings.get(pane) is None:
          self.dimMappings[pane] = []
          
        
        if self.selDimMapping.get(pane) is None:
          self.selDimMapping[pane] = planeName

        self.dimMappings[pane].append(planeName)

      for pane in self.selDimMapping:
        avail = self.dimMappings.get(pane, [])
        if self.selDimMapping[pane] not in avail:
          if avail:
            self.selDimMapping[pane] = avail[0]
          else:
            self.selDimMapping[pane] = None

    winData = [[getWindowPaneName(p), p] for p in self.dimMappings.keys()]
    winData.sort()
    
    return winData

  def changeRootPeakList(self, peakList):

    if peakList is not self.peakList:
      self.setAppDataOptions()
      self.peakList = peakList
      self.rootPeak = None
      self.updateRootPeakList()
      self.updateRootWindow()
      self.updateAfter()

  def changeRootWindow(self, rootPane):

    if rootPane is not self.rootPane:
      self.rootPane = rootPane
      self.updateWindowPulldown()


  def navigateToRoot(self):
    
    if self.rootPeak:
      maxMarks = self.marksEntry.get()
      rootPane  = self.rootPanePulldown.getObject()
      if rootPane:
                
        while len(self.marks) >= maxMarks:
          if self.marks:
            oldMark = self.marks.pop(0)
            if not oldMark.isDeleted:
              oldMark.delete()

          else:
            break

        if maxMarks > 0:
          mark = createPeakMark(self.rootPeak, lineWidth=2.0, remove=False)
          self.marks.append(mark)

        if self.followRootSelect.get():
          windowFrame = rootPane.getWindowFrame()
          windowFrame.gotoPeak(self.rootPeak)
 
 
  def navigateToTarget(self):
      
    if self.rootPeak and self.rootPane and self.followTargetSelect.get():

      position = getPeakWindowPosition(self.rootPeak, self.rootPane,
                                       useDefault=False)

      # Blot out non-root positions from HNCO, HNcoCA
      if len(self.rootPeak.peakDims) > 2:
        dimMapping = getPeakDimAxisMapping(self.rootPeak, self.rootPane)
        dataDims = self.getRootDataDims()
        axisPanels = self.rootPane.sortedAxisPanels()
        
        for i, axisPanel in enumerate(axisPanels):
          peakDim = dimMapping.get(axisPanel.label)
          
          if peakDim.dataDim not in dataDims:
            position[i] = None

      windowZplanes = findOrthogonalWindows(self.rootPane, [position,], minDims=2)
      for (planeName0, windowPane, positions) in windowZplanes:
        planeName = planeName0.split(' in ')[0]
        if (windowPane in self.windowPanes) and \
           (planeName == self.selDimMapping.get(windowPane)):
           
          windowFrame = windowPane.getWindowFrame()
          windowFrame.gotoPosition(position=positions[0])
  
  def getWaterExclusionRegion(self):
    
    waterMinPpm=self.minWaterEntry.get() or 4.95
    waterMaxPpm=self.maxWaterEntry.get() or 4.95     
    waterExclusion = [waterMinPpm, waterMaxPpm]
    waterExclusion.sort()
    
    return waterExclusion
           
  def pickAssignTarget(self):

    if self.rootPeak and self.targetPeakLists:
      tolerances = [self.tolEntry1.get() or 0.1, self.tolEntry2.get() or 0.1]
      diagTolerance=self.diagEntry.get() or 0.0
      waterExclusion = self.getWaterExclusionRegion()
      
      for peakList in self.targetPeakLists:
        pickAssignSpecFromRoot([self.rootPeak,], peakList, 
                               tolerances=tolerances,
                               progressBar=None,
                               diagTolerance=diagTolerance,
                               waterExclusion=waterExclusion)


  def pickAssignAllTarget(self):
  
    if showOkCancel('Confirm','Pick and assign all targets?', parent=self):
      if self.peakList and self.targetPeakLists:
        tolerances = [self.tolEntry1.get() or 0.1, self.tolEntry2.get() or 0.1]
        diagTolerance=self.diagEntry.get() or 0.0
        waterExclusion = self.getWaterExclusionRegion()
 
        for peakList in self.targetPeakLists:
          spectrum   = peakList.dataSource
          experiment = spectrum.experiment
 
          self.progressBar.setText('Working on %s:%s:%s' % (experiment.name,spectrum.name,peakList.serial))
          pickAssignSpecFromRoot(self.peakList.peaks, peakList,
                                 tolerances=tolerances,
                                 progressBar=self.progressBar,
                                 diagTolerance=diagTolerance,
                                 waterExclusion=waterExclusion)

          
  def assignTarget(self):

    if self.rootPeak and self.targetPeakLists:
      tolerances = [self.tolEntry1.get() or 0.1, self.tolEntry2.get() or 0.1]
      diagTolerance=self.diagEntry.get() or 0.0
      waterExclusion = self.getWaterExclusionRegion()
      
      for peakList in self.targetPeakLists:
        pickAssignSpecFromRoot([self.rootPeak,], peakList, pickNew=False,
                               tolerances=tolerances, progressBar=None,
                               diagTolerance=diagTolerance,
                               waterExclusion=waterExclusion)


  def assignAllTarget(self):
  
    if showOkCancel('Confirm','Assign root resonances to all targets?', parent=self):
      if self.peakList and self.targetPeakLists:
        tolerances = [self.tolEntry1.get() or 0.1, self.tolEntry2.get() or 0.1]
        diagTolerance=self.diagEntry.get() or 0.0
        waterExclusion = self.getWaterExclusionRegion()
 
        for peakList in self.targetPeakLists:
          spectrum   = peakList.dataSource
          experiment = spectrum.experiment
 
          self.progressBar.setText('Working on %s:%s:%s' % (experiment.name,spectrum.name,peakList.serial))
          pickAssignSpecFromRoot(self.peakList.peaks, peakList, pickNew=False,
                                 tolerances=tolerances, progressBar=self.progressBar,
                                 diagTolerance=diagTolerance,
                                 waterExclusion=waterExclusion)

  
  def nextTarget(self):
  
    if self.rootPeak:
      peaks = self.rootPeakMatrix.objectList
      index = peaks.index(self.rootPeak)
      peak  = peaks[(index+1) % len(peaks)]
      self.rootPeakMatrix.selectObject(peak)
 
  
  def prevTarget(self):

    if self.rootPeak:
      peaks = self.rootPeakMatrix.objectList
      index = peaks.index(self.rootPeak)
      peak  = peaks[(index-1) % len(peaks)]
      self.rootPeakMatrix.selectObject(peak)

  
  def assignNonRootTarget(self):
  
    if self.rootPeak and self.targetPeakLists:
      found = False
      for peakList in self.targetPeakLists:
        if self.nonRootDict.get(self.getPeakListId(peakList)):
          found = True
          break
      
      if not found:
        showWarning('Warning','No spectra are selected for non-root dim assignment.', parent=self)
        return
      
      refExps, refExpsCO = getSeqAssignRefExperiments(self.project)
      assignType = self.assignTypeSelect.get()
      tolerances = [self.tolEntry1.get() or 0.1, self.tolEntry2.get() or 0.1]
      
      for peakList in self.targetPeakLists:
        if self.nonRootDict.get(self.getPeakListId(peakList)):
          assignSpecNonRootResonances([self.rootPeak,], peakList,
                                      diagTolerance=self.diagEntry.get() or 0.3,
                                      waterMinPpm=self.minWaterEntry.get() or 4.95,
                                      waterMaxPpm=self.maxWaterEntry.get() or 4.95,
                                      tolerances=tolerances, progressBar=None,
                                      assignType=assignType, refExpsCO=refExpsCO)
      
        

  def assignNonRootAllTarget(self):
  
    if showOkCancel('Confirm','Assign non-root resonances to all targets?', parent=self):
      if self.peakList and self.targetPeakLists:
        found = False
        for peakList in self.targetPeakLists:
          if self.nonRootDict.get(self.getPeakListId(peakList)):
            found = True
            break
 
        if not found:
          showWarning('Warning','No spectra are selected for non-root dim assignment.', parent=self)
          return

      
        refExps, refExpsCO = getSeqAssignRefExperiments(self.project)
        assignType = self.assignTypeSelect.get()
        tolerances = [self.tolEntry1.get() or 0.1, self.tolEntry2.get() or 0.1]
        
        for peakList in self.targetPeakLists:
          if self.nonRootDict.get(self.getPeakListId(peakList)):
            spectrum   = peakList.dataSource
            experiment = spectrum.experiment
 
            self.progressBar.setText('Working on %s:%s:%s' % (experiment.name,spectrum.name,peakList.serial))
            assignSpecNonRootResonances(self.peakList.peaks, peakList,
                                        diagTolerance=self.diagEntry.get() or 0.3,
                                        waterMinPpm=self.minWaterEntry.get() or 4.95,
                                        waterMaxPpm=self.maxWaterEntry.get() or 4.95,
                                        tolerances=tolerances, progressBar=self.progressBar,
                                        assignType=assignType, refExpsCO=refExpsCO)
  
  
  def removeWindow(self):

    windowPane = self.windowPane
    if windowPane in self.windowPanes:
      self.windowPanes.remove(windowPane)
      self.windowPane = None
      self.updateWindows(windowPane)
  
  
  def addWindow(self):

    if self.peakList:
      windowPane = self.windowPulldown.getObject()
 
      if windowPane and (windowPane not in self.windowPanes):
        self.windowPanes.append(windowPane)
        self.updateWindows(windowPane)
      

  def selectRootPeak(self, obj, row, col):

    if obj:
      self.rootPeak = obj
      self.updateButtons()
      self.navigateToRoot() # Always done for marker

      if self.followTargetSelect.get():
        self.navigateToTarget()


  def selectWindow(self, obj, row, col):
  
    self.windowPane = obj
    self.updateButtons()

  
  def updateWindowPulldown(self):
  
    windowData = self.getTargetWindows()

    index  = 0
    names  = []
    windowPanes = []
 
    if windowData:
      names   = [x[0] for x in windowData]
      windowPanes = [x[1] for x in windowData]
      
      newWindows = []
      for windowPane in self.windowPanes:
        if windowPane in windowPanes:
          newWindows.append(windowPane)
          
          # Remove already present windows from 'to-add' list
          ii = windowPanes.index(windowPane)
          del names[ii]
          del windowPanes[ii]
      
      if newWindows != self.windowPanes:
        self.windowPanes = newWindows
      
      self.updateWindows()
    

    self.windowPulldown.setup(names, windowPanes, index)
  
  
  def updateRootWindow(self):
  
    windowData = self.getRootWindows()

    index  = -1
    names  = []
    rootPane = self.rootPane
     
    if windowData:
      names = [x[0] for x in windowData]
      panes = [x[1] for x in windowData]
 
      if rootPane not in panes:
        rootPane = panes[0]
 
      index = panes.index(rootPane)
 
    else:
      rootPane = None
      panes = []

    if rootPane is not self.rootPane:
      self.rootPane = rootPane
      self.updateWindowPulldown()
 
    self.rootPanePulldown.setup(names, panes, index)
 
  def getAppDataOptions(self):
    
    project = self.project
    analysisProject = self.analysisProject
  
    app  = project.application
    data = analysisProject.linkPeakListsData
    if data:
      try:
        options = cPickle.loads(data)
      except Exception:
        options = {}  
    else:
      options = {}  

    if options.get('peakList'):
      s1,s2,s3 = options['peakList']
      experiment = self.nmrProject.findFirstExperiment(serial=s1)
      if experiment:
        spectrum = experiment.findFirstDataSource(serial=s2)
        if spectrum:
          self.peakList = spectrum.findFirstPeakList(serial=s3)
          self.rootPeak = None

    if options.get('rootWindow'):
      serials = options['rootWindow']
      self.rootPane = None
      
      if type(serials) == type(1):
        window = analysisProject.findFirstSpectrumWindow(serial=serials)
        
        if window:
          self.rootPane = window.findFirstSpectrumWindowPane()
        
      else:
        serial1, serial2 = serials
        window = analysisProject.findFirstSpectrumWindow(serial=serial1)
     
        if window:
          self.rootPane = window.findFirstSpectrumWindowPane(serial=serial2)
 
    if options.get('windowPanes'):
      for serial1, serial2 in options['windowPanes']:
        window = analysisProject.findFirstSpectrumWindow(serial=serial1)
        
        if window:
          windowPane = window.findFirstSpectrumWindowPane(serial=serial2)
        
          if windowPane:
            self.windowPanes.append(windowPane)
 
    elif options.get('windows'):
      self.windowPanes = []
         
      for serial in options['windows']:
        window = analysisProject.findFirstSpectrumWindow(serial=serial)
        if window:
          self.windowPanes.append(window.findFirstSpectrumWindowPane())

    self.diagEntry.set(options.get('diagTol',0.3))
    self.minWaterEntry.set(options.get('minWater',4.95))
    self.maxWaterEntry.set(options.get('maxWater',4.95))
    self.followRootSelect.set(options.get('followRoot',True))
    self.followTargetSelect.set(options.get('followTarget',True))
    self.assignTypeSelect.set(options.get('assignType',False))
    
    mappings = options.get('selDimMapping', [])
    self.selDimMapping = {}
    if type(mappings) is type({}):
      for pane in self.windowPanes:
        mapping = mappings.get(pane.spectrumWindow.serial)
        avail = self.dimMappings.get(pane, [])
        
        if mapping not in avail:
          if avail:
            mapping = avail[0]
          else: 
            mapping = None
        
        self.selDimMapping[pane] = mapping
        
    else:
      for i, pane in enumerate(self.windowPanes):
        self.selDimMapping[pane] = mappings[i]
    
    self.nonRootDict   = options.get('nonRootDict', {})

    
  def setAppDataOptions(self):
    
    project = self.project
 
    options = {}
    options['diagTol']      = self.diagEntry.get()
    options['minWater']     = self.minWaterEntry.get()
    options['maxWater']     = self.maxWaterEntry.get()
    options['followRoot']   = self.followRootSelect.get()
    options['followTarget'] = self.followTargetSelect.get()
    options['assignType']   = self.assignTypeSelect.get()
    options['windowPanes']  = [(p.spectrumWindow.serial, p.serial) for p in self.windowPanes]
    options['selDimMapping']= [self.selDimMapping[w] for w in self.windowPanes]
    options['nonRootDict']  = self.nonRootDict
     
    if self.rootPane:
      options['rootWindow'] = (self.rootPane.spectrumWindow.serial,
                               self.rootPane.serial) 
    if self.peakList:
      options['peakList'] = self.getPeakListId(self.peakList)
    
    try:
      data = cPickle.dumps(options)
      self.analysisProject.linkPeakListsData=data
    except:
      pass  

  def setDataDimTolerances(self, *opt):
  
    dataDim1, dataDim2 = self.getRootDataDims()
        
    analysisDataDim1 = getAnalysisDataDim(dataDim1)
    analysisDataDim2 = getAnalysisDataDim(dataDim2)
    
    tol1 = analysisDataDim1.assignTolerance
    tol2 = analysisDataDim2.assignTolerance
 
    analysisDataDim1.assignTolerance = self.tolEntry1.get() or tol1
    analysisDataDim2.assignTolerance = self.tolEntry2.get() or tol2
    

  def updateRootPeakList(self):

    peakListData = self.getRootPeakLists()

    index    = -1
    names    = []
    peakList = self.peakList
 
    if peakListData:
      names     = [x[0] for x in peakListData]
      peakLists = [x[1] for x in peakListData]
 
      if peakList not in peakLists:
    
        peakList = peakLists[0]
 
      index = peakLists.index(peakList)
 
    else:
      peakList = None
      peakLists = []

    if peakList is not self.peakList:
      self.peakList = peakList
      self.rootPeak = None
      self.updateRootWindow()

    if self.peakList:      
      spectrum = self.peakList.dataSource
      dataDim1, dataDim2 = self.getRootDataDims()
 
      analysisDataDim1 = getAnalysisDataDim(dataDim1)
      analysisDataDim2 = getAnalysisDataDim(dataDim2)

      dataDims = self.peakList.dataSource.sortedDataDims()
 
      isotopeStr1 = ','.join(getPrimaryDataDimRef(dataDim1).expDimRef.isotopeCodes)
      isotopeStr2 = ','.join(getPrimaryDataDimRef(dataDim2).expDimRef.isotopeCodes)
 
      self.tolEntry1.set( analysisDataDim1.assignTolerance )
      self.tolLabel1.set('Root %s Dim %d' % (isotopeStr1,dataDim1.dim) )
      self.tolEntry2.set( analysisDataDim2.assignTolerance )
      self.tolLabel2.set('Root %s Dim %d' % (isotopeStr2,dataDim2.dim) )
     
    self.updateAfter()
    self.peakListPulldown.setup(names, peakLists, index)

  def getPeakListId(self, peakList):
  
    spectrum = peakList.dataSource

    return (spectrum.experiment.serial,spectrum.serial,peakList.serial)
 
  def getIsActive(self, peakList):
    
    for windowPane in self.windowPanes:
      for view in windowPane.spectrumWindowViews:
        if view.analysisSpectrum.dataSource is peakList.dataSource:
          if view.isPosVisible or view.isNegVisible:
            return True
    
    return False
 
  def togglePeakList(self, peakList):
 
    if peakList:
      views = set([])
      for windowPane in self.windowPanes:
        for view in windowPane.spectrumWindowViews:
          if view.analysisSpectrum.dataSource is peakList.dataSource:
            views.add(view)
            
      onViews = []
      offViews = []
      for view in views:
        if view.isPosVisible or view.isNegVisible:
          onViews.append(view)
        else:
          offViews.append(view)  
     
      if onViews:
        for view in onViews:
          toggleSpectrum(view.spectrumWindowPane.spectrumWindow,
                         spectrum=view.analysisSpectrum.dataSource)
     
      elif offViews:
        for view in offViews:
          toggleSpectrum(view.spectrumWindowPane.spectrumWindow,
                         spectrum=view.analysisSpectrum.dataSource)
        
      
  def updateSpectra(self, object=None):
 
    if object:
      self.updateRootPeakList()
      
    self.targetPeakLists = []
    peakLists = self.getTargetPeakLists()
 
    activeColors = ['#B0FFB0','#B0FFB0',None,None,None]
    inactiveColors = ['#FFB0B0','#FFB0B0',None,None,None]
 
    textMatrix = []
    objectList = []
    colorMatrix = []
    for peakList in peakLists:
      spectrum   = peakList.dataSource
      experiment = spectrum.experiment
      
      refExpName = '-'
      if experiment.refExperiment:
        refExpName = experiment.refExperiment.name
        
      isotopes   = []
      for dataDim in spectrum.dataDims:
        isotopeDict = {}
        for expDimRef in dataDim.expDim.expDimRefs:
          for isotope in expDimRef.isotopeCodes:
            isotopeDict[isotope] = None
            
        isotopes.append( ','.join( isotopeDict.keys() )  )
    
      yesNo = 'No'
      if self.nonRootDict.get(self.getPeakListId(peakList)):
        yesNo = 'Yes'
      
      if self.getIsActive(peakList):
        isActive = 'Yes'
        colorMatrix.append(activeColors)
        self.targetPeakLists.append(peakList)
      else:
        isActive = 'No'
        colorMatrix.append(inactiveColors)
      
      datum = ['%s:%s' % (experiment.name,spectrum.name),
               isActive,
               yesNo,
               refExpName,
                ' '.join(isotopes)]
      
      textMatrix.append(datum)
      objectList.append(peakList)

    self.spectrumMatrix.update(textMatrix=textMatrix,
                               colorMatrix=colorMatrix,
                               objectList=objectList)
    self.setAppDataOptions()
   
  def updateWindowsAfter(self, window=None):
    
    self.after_idle(self.updateWindows) # Axis panel notifiers need a bit of time
  
  def updateWindowPaneAfter(self, windowPane=None):
    
    if windowPane:
      self.after_idle(lambda :self.updateWindows(windowPane))

  def updateWindows(self, windowPane=None):
    
    if windowPane:
      self.updateRootWindow()
      self.updateWindowPulldown()
      if windowPane.isDeleted and (windowPane in self.windowPanes):
        if self.windowPane is windowPane:
          self.windowPane = None
        self.windowPanes.remove(windowPane)
    
    # '#','Name','Axes','Selected mapping'
    textMatrix = []
    objectList = []
    for windowPane in self.windowPanes:
      datum = []
      datum.append(windowPane.spectrumWindow.serial)
      datum.append(getWindowPaneName(windowPane))
      #datum.append( ' '.join([axisPanel.axisType.name for axisPanel in windowPane.axisPanels]) )
      
      mapping = self.selDimMapping.get(windowPane)
      avail = self.dimMappings.get(windowPane)
      
      if avail and (mapping not in avail):
        self.selDimMapping[windowPane] = avail[0]
      
      datum.append(self.selDimMapping[windowPane])
      
      textMatrix.append(datum)
      objectList.append(windowPane)
  
    self.updateSpectra()
    self.windowMatrix.update(textMatrix=textMatrix,
                             objectList=objectList)
     
  def updateAfter(self, obj=None):
 
    if obj:
      if obj.className == 'Peak':
        peak = obj
      else: # must be peakDim
        peak = obj.peak
    
      if peak.peakList is not self.peakList:
        return
         
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)

  def update(self):
  
    textMatrix = []
    objectList = []
    headingList = ['#','Assign F1','Assign F2','Shift F1','Shift F2']

    if self.peakList:
      dataDim1, dataDim2 = self.getRootDataDims()
       
      i = 1
      for dataDim in (dataDim1, dataDim2):
        dim = dataDim.dim
        headingList[i]   = 'Assign F%d' % (dim)
        headingList[i+2] = 'Shift F%d' %  (dim)
        i += 1
      
      for peak in self.peakList.peaks:
        peakDim0 = peak.findFirstPeakDim(dataDim=dataDim1)
        peakDim1 = peak.findFirstPeakDim(dataDim=dataDim2)
        
        datum = []
        datum.append(peak.serial)
        datum.append(peakDim0.annotation or '-')
        datum.append(peakDim1.annotation or '-')
        datum.append(peakDim0.value)
        datum.append(peakDim1.value)

        textMatrix.append(datum)
        objectList.append(peak)
  
    self.rootPeakMatrix.update(textMatrix=textMatrix,
                               objectList=objectList,
                               headingList=headingList,
                               tipTexts=self.rootPeakTipTexts)
  
    self.updateButtons()
    self.setAppDataOptions()
    self.waiting = False

  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)
  
    self.setAppDataOptions()
  
    BasePopup.destroy(self)

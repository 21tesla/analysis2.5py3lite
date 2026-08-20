"""
======================COPYRIGHT/LICENSE START==========================

LinkNoeResonances.py: Part of the CcpNmr Analysis program

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

from math import sqrt

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

from ccp.util.NmrExpPrototype import longRangeTransfers

from ccpnmr.analysis.popups.BasePopup     import BasePopup
from ccpnmr.analysis.core.AssignmentBasic import makeResonanceGuiName, clearPeakDim, makePeakDimAnnotation
from ccpnmr.analysis.core.AssignmentBasic import findMatchingPeakDimShifts, assignResToDim
from ccpnmr.analysis.core.ExperimentBasic import getThroughSpacePeakLists, getOnebondDataDims
from ccpnmr.analysis.core.MarkBasic       import createPeakMark
from ccpnmr.analysis.core.MoleculeBasic   import areResonancesBound
from ccpnmr.analysis.core.PeakBasic       import getPeakHeight, pickPeak
from ccpnmr.analysis.core.PeakBasic       import findSymmetryPeaks, getPeakVolume
from ccpnmr.analysis.core.StructureBasic  import getAtomSetsDistance
from ccpnmr.analysis.core.Util            import getAnalysisPeakList, getAnalysisDataDim
from ccpnmr.analysis.core.WindowBasic     import getPeakDimAxisMapping
from ccpnmr.analysis.core.WindowBasic     import getWindowPaneName, displayPeakStrips

def testNoePopup(argServer):

  popup = LinkNoeResonancesPopup(argServer.parent)
  popup.open()


assignColor = '#80D080'

class LinkNoeResonancesPopup(BasePopup):
  """
  **Assign Contributions to NOE Peaks Using Structures and Chemical Shifts**
  
  This popup window is designed to assist in the latter stages of assigning NOE
  and other through-space spectra that are used in the generation of
  macromolecular structures. The system uses the matching of chemical shifts to
  peak positions, and a structure if available, to suggest assignments for peaks.
  Guiding structures are used to eliminate assignment possibilities where the
  resonances would be too far apart in space to give a contribution to the
  signal. Such structures may include comparative (homology) models or
  unfinalised structures from part way though the structure generation process,
  for example as produced by ARIA or CYANA programs.

  **User Operation**

  To operate this tool the user selects a NOESY or other though-space peak list
  that will be assigned. Various options may be chosen, including how assignment
  possibilities are displayed on a graphical structure view and what the
  tolerances are for considering an assignment possibility. The distance
  tolerance is the maximum allowable separation between potential peak dimension
  assignments (connected by a through-space transfer) if a guide structure is
  present. The F1, F2 etc. dimensional isotope tolerances are the maximum
  allowed chemical shift differences for matching a peak position to a resonance
  of known chemical shift. The "Navigation Windows" table allows the user to
  specify which spectrum windows will be used for navigation; to show a selected
  peak's location. This table displays all windows that the peak list may be
  viewed in, and the "Use?" column may be toggled to activate a window for
  navigation.
  
  **Reciprocal Return Peaks**
  
  The active windows can not only display the location of a selected peak, but
  also a series of strips (window sub-divisions) that relate to any reciprocal,
  symmetry related (return) peaks. This is useful for the resolution of
  assignment ambiguity. If a peak is thought to be due to an NOE from resonance
  A to resonance B a reciprocal NOE (from B to A) will usually be present to
  confirm the assignment possibility, assuming the return peak would be visible
  in the spectrum. If strips for return peaks are required the user should
  select the "Strip Return Peaks" option at the top of the settings tab. 

  **Peak Assignments**

  In the "Peak Assignments" tab the top table lists all of the peaks that are
  present in the selected peak list. Clicking on a peak row in this table will
  cause the lower table to be filled with the assignment possibilities for the
  peak. If navigation windows and a guide structure display is active selecting
  a peak in this way will also navigate to the appropriate spot in the spectrum
  and show structural connections.

  Peak assignment possibilities in the lower table are, by default, displayed in
  a ranked order according to how close the peak position is to the chemical
  shifts and how distant in space the relevant atoms are in any guide structure.
  The chemical shift differences are weighted by the per-dimension tolerances so
  that values from dissimilar isotopes may be combined. A peak assignment is
  made by selecting the specific row for the resonances you wish to assign to
  the selected peak and clicking [Assign Selected]. An individual dimension
  assignment may be removed by clicking on the appropriate cell in the "Possible
  Assignments" table and clicking [Delete Assignment], or all assignments for
  the peak may be cleared with [Delete All Assignments].

  The [Predict Peaks] function is designed to show you where the listed
  assignment possibilities would make peaks in the spectrum, given their
  recorded chemical shifts. This option makes synthetic peak crosses to show
  where on the spectrum the resonance intersections would be, and if possible
  what structural distance they represent. From these the user may be able to
  make a better judgement of which possibilities are likely to cause a peak.
  Such synthetic peaks are made in an entirely separate peak list to the main
  one, and their list may be deleted at any time without adverse effect.

  **Caveats and Tips**
  
  The "Strip Return Peaks" option will override any existing strips within the
  selected spectrum windows. If the windows should not be affected in this way
  the option should be switched off.
    
  The "Tol. Weighted Combined Delta" value in the lower table is a weighted
  chemical shift distance (structure distance is not included) that can be used
  an assignment criterion, as long as your shift tolerances are set
  appropriately. The actual calculation is to take the shift difference for each
  dimension and divide by the tolerance for that dimension (normalise to be <=
  1.0). All dimensions are then combined by taking the square root of the sum of
  squares, i.e. to give a "shift distance". The maximum value will be the square
  root of the number of dimensions.

  It should be noted that the [Assign Whole Peak List] will assign all of the 
  possibilities, as listed in the lower table, for all the peaks in the entire
  peak list. Hence, when using this function it is recommended that suitably
  narrow distance and chemical shift tolerances are set.

  The chemical shift list that is used for matching resonances to peak positions
  is the one set on the experiment record for the selected peak list.

  The experiment record of a peak list must be set as a through-space type (e.g.
  NOESY) for it to be considered.

  This system will work with non-NOESY though-space experiments such as spin
  diffusion types used for solid-state samples. However, at present assignment
  contribution fractions are only calculated using the Intensity^(1/6) summation
  appropriate for NOE experiments.

  This system takes no account of missing assignments; resonances that may cause
  a peak but do not have a known chemical shift.
  """

  def __init__(self, parent, *args, **kw):

    self.waiting    = False
    self.waitingAssign = False
    self.windowPanes = []
    self.peakList   = None
    self.peak       = None
    self.guiParent  = parent
    self.project    = parent.project
    self.assignment = None
    self.structure  = None
    self.isotopes   = []
    self.dim = None
    
    self.tolLabels  = []
    self.tolEntries = []
    self.heterDims  = []
    self.hydroDims  = []
        
    BasePopup.__init__(self, parent=parent, title='Assignment : NOE Contributions', **kw)

  def open(self):
  
    BasePopup.open(self)
    self.updatePeakLists()

  def close(self):
  
    self.setAppDataOptions()
    BasePopup.close(self)

  def body(self, guiFrame):

    self.geometry('500x700')

    guiFrame.expandGrid(0,0)
    
    tipTexts = ['Sets the parameters required for matching peaks to NOE assignments, including which peak list and structure to use',
                'A table of the peaks that NOE assignments are considered for, and their structure & shift filtered assignment possibilities']
    options = ['Peak List & Display Settings','Peak Assignments']
    tabbedFrame = TabbedFrame(guiFrame, options=options,
                              grid=(0,0), tipTexts=tipTexts)
    
    frameA, frameB = tabbedFrame.frames
    frameA.grid_columnconfigure(0, weight=1)
    frameB.grid_columnconfigure(0, weight=1)
    
    #
    # Options
    #
    #
    
    
    row = 0
    frame = Frame(frameA, grid=(row,0), sticky='ew')
    frame.grid_columnconfigure(5, weight=1)

    label = Label(frame, text='NOE Peak List:', grid=(0,0), sticky='e')
    tipText = 'Selects which NOE peak list will be considered for the matching of close resonances; in terms of both shift and distance)'
    self.peakListPulldown = PulldownList(frame, callback=self.changePeakList,
                                         grid=(0,1), tipText=tipText)

    subFrame = Frame(frame, grid=(1,0), gridSpan=(1,6), sticky='ew')
    subFrame.grid_columnconfigure(4, weight=1)

    label = Label(subFrame, text='Structure Display:', grid=(0,0))
    tipText = 'Selects which kind of assignment possibilities, for a peak, to display on a graphical structure display'
    self.strucDispPulldown = PulldownList(subFrame, texts=['All','Assigned','None'],
                                          index=2, grid=(0,1), tipText=tipText)

    label = Label(subFrame, text=' Structure:', grid=(0,2), sticky='e')
    tipText = 'Selects which structure ensemble will be used for calculating atomic distances within the graphical structure display'
    self.structurePulldown = PulldownList(subFrame, callback=self.changeStructure,
                                          grid=(0,3), tipText=tipText)

    subFrame = Frame(frame, grid=(0,0))
    subFrame.grid(row=2, column=0, columnspan=6, sticky='ew')
    subFrame.grid_columnconfigure(7, weight=1)

    label = Label(subFrame, text='Aliased Possible:', grid=(0,0))
    tipText = 'Sets whether the peaks could be caused by aliased resonances; not at their real ppm value, but a whole number of sweep widths away'
    self.aliasSelect = CheckButton(subFrame, callback=None, grid=(0,1),
                                   selected=True, tipText=tipText)

    label = Label(subFrame, text='Focus Structure:', grid=(0,2))
    tipText = 'Sets whether to rotate the graphical structure display to focus on the current assignment possibilities'
    self.focusSelect = CheckButton(subFrame, callback=None, grid=(0,3),
                                   selected=True, tipText=tipText)

    label = Label(subFrame, text='Mark Peak:', grid=(0,4))
    tipText = 'Sets whether to mark the selected peak position, as displayed in the selected navigation windows'
    self.markSelect = CheckButton(subFrame, callback=None, grid=(0,5),
                                  selected=False, tipText=tipText)

    label = Label(subFrame, text='Auto Find Reciprocal Peaks:', grid=(0,6))
    tipText = 'Whether to display strips for possible (symmetry related) NOE reciprocating return peaks, e.g. for an A-B peak look for possible B-A peaks.'
    self.returnSelect = CheckButton(subFrame, callback=None, grid=(0,7),
                                  selected=False, tipText=tipText)
     
    row +=1
    div = LabelDivider(frameA, text='Tolerances', grid=(row,0))
    
    row +=1
    self.tolFrame = Frame(frameA, grid=(row,0), sticky='ew')
    self.tolFrame.grid_columnconfigure(10, weight=1)

    label = Label(self.tolFrame, text='Distance:', grid=(0,0))
    tipText = 'Sets the upper limit for structural distance, below which pairs of close atom sites may be considered for assignment'
    self.distThresholdEntry = FloatEntry(self.tolFrame, text=8.0, width=8,
                                         grid=(0,1), tipText=tipText)
    
    row +=1
    div = LabelDivider(frameA, text='Navigation Windows', grid=(row,0))
    
    
    row +=1
    frameA.grid_rowconfigure(row, weight=1)
    tipTexts = ['The serial number of the spectrum window, which may be used to navigate to NOE peak locations',
                'The name of the spectrum window, which may be used to navigate to NOE peak locations',
                'Sets whether to use the spectrum window for navigation, if so its view will be moved to display the selected peak location']
    headingList = ['#','Name','Use?']
    editWidgets      = [None, None, None]
    editGetCallbacks = [None, None, self.toggleUseWindow]
    editSetCallbacks = [None, None, None]
    self.windowMatrix = ScrolledMatrix(frameA, headingList=headingList,
                                       editSetCallbacks=editSetCallbacks,
                                       editGetCallbacks=editGetCallbacks,
                                       editWidgets=editWidgets,
                                       highlightType=None,
                                       callback=None,
                                       grid=(row,0), tipTexts=tipTexts)

    #
    # Peak Assignments
    #

    row = 0
    div = LabelDivider(frameB, text='NOE Peaks', grid=(row,0))
    
    row +=1
    #frameB.grid_rowconfigure(row, weight=1)
    tipTexts = ['The serial number of the peak in its peak list',
                'The position of the peak in all of its dimensions',
                'The current resonance assignment annotation of the peak',
                'A textual description of the peak; can be set from the main peak tables, e.g. to mark violated peaks',
                'The height intensity of the peak; usually from the spectrum data matrix,but possibly set manually',
                'The volume integral of the peak; calculated using the method stated in the peak list table']
    headingList = ['#','Position','Assignment','Details','Height','Volume']
    editWidgets      = [None, None, None, None, None, None,]
    editGetCallbacks = [None, None, None, None, None, None,]
    editSetCallbacks = [None, None, None, None, None, None,]
    self.peakMatrix = ScrolledMatrix(frameB, headingList=headingList, 
                                     editSetCallbacks=editSetCallbacks,
                                     editGetCallbacks=editGetCallbacks,
                                     editWidgets=editWidgets, 
                                     callback=self.selectPeak,
                                     grid=(row,0), tipTexts=tipTexts)

    row +=1
    tipTexts = ['Selects the previous peak in the table and display NOE assignment possibilities; causes navigation within spectrum windows if this is set',
                'Selects a peak in the table based upon the peak currently selected in spectrum windows; uses the first if many are selected',
                'Selects the next peak in the table and display NOE assignment possibilities; causes navigation within spectrum windows if this is set',
                'Display strips for possible (symmetry related) NOE reciprocating return peaks, e.g. for an A-B peak look for possible B-A peaks.']
    texts    = ['Previous Peak','Selected Peak','Next Peak','Find Reciprocal Peaks']
    commands = [self.prevPeak,self.selectedPeak,
                self.nextPeak,self.stripSymmetryRelatedPeaks]
    self.peakButtons = ButtonList(frameB, texts=texts, commands=commands,
                                  grid=(row,0), tipTexts=tipTexts)

    row +=1
    frame = Frame(frameB, grid=(row,0))
    frame.expandGrid(0,0)
    
    div = LabelDivider(frame, text='Possible Assignments', grid=(0,0))
    tipText = 'Indicates the equivalent NOE distance of the peak, using the combined assignment contributions and assuming an intensity^(-1/6) summation'
    label = Label(frame, text=u'\u03A3NOE dist: ',
                  grid=(0,1), tipText=tipText)
    self.noeSumLabel = Label(frame, text='', grid=(0,2), tipText=tipText)
    
    row +=1
    frameB.grid_rowconfigure(row, weight=1)
    headingList = ['#',
                   'Reson\nF1','Reson\nF2','Reson\nF3','Reson\nF4',
                   'Dist.',u'Tol. Weighted\nCombined \u0394',
                   u'\u0394F1',u'\u0394F2',u'\u0394F3',u'\u0394F4','Contrib\nFrac.'
                   'Shift\nF1','Shift\nF2','Shift\nF3','Shift\nF4',]
    self.assignMatrix = ScrolledMatrix(frameB, headingList=headingList,
                                       multiSelect=True, grid=(row,0),
                                       highlightType='cell',
                                       callback=self.selectAssignment)
 
    row +=1
    tipTexts = ['Assign the peak to the individual resonance+dimension possibility selected in the above table',
                'Remove the assignment of the individual resonance+dimension possibility selected in the above table',
                'Delete all resonance assignment for all dimensions of the selected peak',
                'Make simulated peaks (in a separate, temporary peak list) at the locations of the displayed resonance intersections',
                'Assign all the peaks in the selected peak list to all of their NOE assignment possibilities; only useful if strict tolerances are set']
    texts = ['Assign\nSelected', 'Delete\nAssignment',
             'Delete All\nAssignments',
             'Predict\nPeaks','Assign Whole\nPeak List']
    commands = [self.assignSelected, self.deleteAssignment,
                self.clearAssignments,
                self.predictPeaks,self.processAll]
    self.assignButtons = ButtonList(frameB, texts=texts, commands=commands,
                                    grid=(row,0), tipTexts=tipTexts)

    buttons = UtilityButtonList(tabbedFrame.sideFrame, helpUrl=self.help_url,
                                grid=(0,0), sticky='e')

    self.administerNotifiers(self.registerNotify)
      
    self.getAppDataOptions()  
      
    self.updatePeakLists()

  def administerNotifiers(self, notifyFunc):

    for func in ('__init__', 'delete','setName'):
      for clazz in ('ccp.nmr.Nmr.DataSource', 'ccp.nmr.Nmr.Experiment'):
        notifyFunc(self.updatePeakLists, clazz, func)
      
      notifyFunc(self.updateWindowsAfter, 'ccpnmr.Analysis.SpectrumWindow', func)
      

    for func in ('__init__', 'delete',):
      notifyFunc(self.updatePeaksAfter, 'ccp.nmr.Nmr.Peak',     func)
      notifyFunc(self.updatePeaksAfter, 'ccp.nmr.Nmr.PeakList', func)
      notifyFunc(self.updateAssignmentsAfter, 'ccp.nmr.Nmr.PeakDimContrib', func)
      notifyFunc(self.updateWindowsAfter, 'ccpnmr.Analysis.SpectrumWindowView', func)
  
    for func in ('setPosition','setNumAliasing'):
      notifyFunc(self.updateAssignmentsAfter, 'ccp.nmr.Nmr.PeakDim', func)

    for func in ('setAnnotation',):
      notifyFunc(self.updatePeaksAfter, 'ccp.nmr.Nmr.PeakDim', func)

    for func in ('setDetails',):
      notifyFunc(self.updatePeaksAfter, 'ccp.nmr.Nmr.Peak', func)

    for func in ('__init__', 'delete', 'setValue'):
      notifyFunc(self.updateIntensityAfter, 'ccp.nmr.Nmr.PeakIntensity', func)

    for func in ('delete','__init__'):
      notifyFunc(self.updateStructures, 'ccp.molecule.MolStructure.StructureEnsemble', func)

  def stripSymmetryRelatedPeaks(self):

    if self.peak:
      peaks = findSymmetryPeaks(self.peak, peakLists=[self.peakList,])
      peaks = [self.peak,] + peaks
      if len(peaks) > 9:
        peaks = peaks[:9]
 
      for windowPane in self.windowPanes:
        displayPeakStrips(self.guiParent, peaks, windowPane)
        #self.update_idletasks()
        #displayPeakStrips(self.guiParent, peaks, windowPane)
        windowFrame = windowPane.getWindowFrame()
        if windowPane.spectrumWindow.stripAxis == 'y':
          axis = 'x'
          setDividersTo = windowFrame.setColsTo
        else:
          axis = 'y'
          setDividersTo = windowFrame.setRowsTo
       
        dimMappingA = getPeakDimAxisMapping(peaks[0], windowPane)
        peakDimA = dimMappingA.get(axis)
        positions = [peakDimA.realValue,]
        
        if len(peaks) > 1:
          dimMappingB = getPeakDimAxisMapping(peaks[1], windowPane)
          peakDimB = dimMappingB.get(axis)
          positions.append(peakDimB.realValue)
        
        positions.sort()
        panel = windowPane.findFirstAxisPanel(label=axis)
        setDividersTo(len(positions))

        regions = panel.sortedAxisRegions()
        for i, ppm in enumerate(positions):
          begin, end = regions[i].region
          width = (end - begin)/2.0
          region = (ppm-width, ppm+width)
          regions[i].region = region
                
  def updateIntensityAfter(self, peakIntensity):
    
    if self.waiting:
      return
    
    if peakIntensity:
      self.updatePeaksAfter( peakIntensity.peak)

  def getPeakLists(self):
  
    return getThroughSpacePeakLists(self.project)
   
  def getStructures(self):

    structures = []
    if self.peakList:
      experiment = self.peakList.dataSource.experiment
      
      if not experiment.molSystems:
        resonances = {}
      
        for peak in self.peakList.peaks:
          for peakDim in peak.peakDims:
            for peakDimContrib in peakDim.peakDimContribs:
              resonance = peakDimContrib.resonance
              
              if resonance.resonanceSet:
                resonances[resonance] = True
     
        molSystems = {}
        for resonance in resonances.keys():
          molSystem = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue.chain.molSystem
          molSystems[molSystem] = True
    
        if molSystems:
          experiment.setMolSystems(molSystems.keys())
    
      molSystems = experiment.molSystems or self.project.molSystems
    
      for molSystem in molSystems:
        for structure in molSystem.sortedStructureEnsembles():
          structures.append(structure)
  
    return structures

  def toggleUseWindow(self, windowPane):
  
    if windowPane in self.windowPanes:
      self.windowPanes.remove(windowPane)
      
    else:
      self.windowPanes.append(windowPane)
      
    self.updateWindows()   

  def changePeakList(self, peakList):

    if peakList is not self.peakList:
      self.peakList = peakList
      self.peak     = None
      self.updatePeaksAfter()
      self.updateStructures()
      self.updateAssignmentsAfter()
      self.updateTolerances()
      self.windowPanes = [] # self.getWindows()
      self.updateWindows()
  
  def getPeakListId(self, peakList):
  
    spectrum = peakList.dataSource

    return '%s:%s:%d' % (spectrum.experiment.name,spectrum.name,peakList.serial)
  
  def updatePeakLists(self, obj = None):

    peakLists = self.getPeakLists()

    index    = -1
    names    = []
    peakList = self.peakList
 
    if peakLists:
      names     = [self.getPeakListId(pl) for pl in peakLists]
 
      if peakList not in peakLists:
        peakList = peakLists[0]
 
      index = peakLists.index(peakList)
 
    else:
      peakList = None

    if peakList is not self.peakList:
      self.peakList = peakList
      self.peak = None
      self.updatePeaksAfter()
      self.updateAssignmentsAfter()   
      self.updateTolerances()   
      self.updateStructures()
      self.windowPanes = [] # self.getWindows()
      self.updateWindows()

    self.peakListPulldown.setup(names, peakLists, index)
  
  def updateTolerances(self):
  
    if self.peakList:
      spectrum = self.peakList.dataSource

      boundDims = getOnebondDataDims(spectrum)
      
      self.heterDims = []
      self.hydroDims = []
      self.isotopes  = []
      
      for dataDim in spectrum.sortedDataDims():
        dim = dataDim.dim
        j   = dim-1
        
        if dataDim.className != 'SampledDataDim':
          isotopes = {}
          
          isDistDim = False
          for dataDimRef in dataDim.dataDimRefs:
            expDimRef = dataDimRef.expDimRef
          
            for isotopeCode in expDimRef.isotopeCodes:
              isotopes[isotopeCode] = True
              
            for expTransfer in expDimRef.expTransfers:
              if expTransfer.transferType in longRangeTransfers:
                isDistDim = True

          isotope = ','.join(isotopes.keys())
          self.isotopes.append(isotope)
          
          if isDistDim:
            self.hydroDims.append(dataDim.dim)
    
          text = 'F%d %s:' % (dim,isotope)
          tol  = 2.0 * getAnalysisDataDim(dataDim).assignTolerance
          tipText = 'The shift scaling factor for dimension %d (%s); used to give equivalency to ppm distances for different kinds of isotope' % (dim,isotope)
          
          if dim > len(self.tolLabels):
            label = Label(self.tolFrame, text=text)
            entry = FloatEntry(self.tolFrame, text=tol,
                               width=6, tipText=tipText)
 
            self.tolLabels.append(label)
            self.tolEntries.append(entry)
 
          else:
            self.tolLabels[j].set(text)
            self.tolEntries[j].set(tol)
            self.tolEntries[j].toolTip.text = tipText
 
          self.tolLabels[j].grid(row=0,column=(2*j)+2)
          self.tolEntries[j].grid(row=0,column=(2*j)+3)
          
        else:
          self.isotopes.append('Sampled')
          text = 'F%d Sampled' % (dim)

          if dim > len(self.tolLabels):
            label = Label(self.tolFrame, text=text)
            entry = FloatEntry(self.tolFrame, text=0,
                               width=5, tipText='**Tolerance not used for sampled dimensions**')
 
            self.tolLabels.append(label)
            self.tolEntries.append(entry)
 
          else:
            self.tolLabels[j].set(text)
 
          self.tolLabels[j].grid(row=0,column=2*j)
          self.tolEntries[j].grid_forget()          
        
      if len(self.tolLabels) > spectrum.numDim:
        for j in range(spectrum.numDim,len(self.tolLabels)):
          self.tolLabels[j].grid_forget()
          self.tolEntries[j].grid_forget()
    
      for dimH in self.hydroDims:
        dimX = None
      
        for dim1, dim2 in boundDims:
          if dimH == dim1.dim:
            dimX = dim2.dim
          elif dimH == dim2.dim:
            dimX = dim1.dim
         
        self.heterDims.append(dimX)
        
    else:
      self.heterDims = []
      self.hydroDims = []
      for label in self.tolLabels:
        label.grid_forget()
        
      for entry in self.tolEntries:
        entry.grid_forget()
  
  def changeStructure(self, structure):

    if structure is not self.structure:
      self.structure = structure
      self.updateAssignments()
  
  def updateStructures(self, object=None):
  
    names = []
    index = 0
    structures = self.getStructures()

    if structures:
      if self.structure in structures:
        structure = self.structure
      else:
        structure = structures[0]
 
      names = [str(x.ensembleId) for x in structures]
      index = structures.index(structure)
 
    else:
      structure = None

    self.structurePulldown.setup(names, structures, index)

    if structure is not self.structure:
      self.structure = structure
    
  
  def selectPeak(self, object, row, col):
  
    self.peak = object
    self.updateAssignmentsAfter()

    if self.markSelect.get():
      createPeakMark(self.peak, lineWidth=1.0)
    
    if self.returnSelect.get():
      self.stripSymmetryRelatedPeaks()
    
    else:
      for windowPane in self.windowPanes:
        windowFrame = windowPane.getWindowFrame()
        windowFrame.setColsTo(1)
        windowFrame.setRowsTo(1)
        windowFrame.gotoPeak(self.peak)
  
      
  def nextPeak(self):
  
    peaks = self.peakMatrix.objectList
    
    if self.peak in peaks:
      index = peaks.index(self.peak)
    else:
      index = -1
      
    index += 1  
    if index >= len(peaks):
      index = 0  
  
    peak = peaks[index]
    self.peakMatrix.selectObject(peak)
  
  def prevPeak(self):

    peaks = self.peakMatrix.objectList
    
    if self.peak in peaks:
      index = peaks.index(self.peak)
    else:
      index = 1
      
    index -= 1  
    if index < 0:
      index = len(peaks)-1  
  
    peak = peaks[index]
    self.peakMatrix.selectObject(peak)
  
  def selectedPeak(self):
 
    peaks = [p for p in self.guiParent.currentPeaks if p.peakList is self.peakList]
    
    if peaks:
      peak = peaks[-1]
      self.peakMatrix.selectObject(peak)
       
    elif self.guiParent.currentPeaks:
      showWarning('Warning', 'Selected peak not in current peak list', parent=self)  

    else:
      showWarning('Warning', 'No selected peak', parent=self)  
  
  def deleteAssignment(self):
  
    if self.peak and self.assignment:
      if self.dim is None:
        nDim = len(self.peak.peakDims)
        msg = 'Must select an assignment column F1-F%d' % nDim
        showWarning('Warning',msg , parent=self)  
        return
    
      score, resonances, dist = self.assignment
      
      peakDim = self.peak.findFirstPeakDim(dim=self.dim)
      
      if peakDim:
        resonance = resonances[self.dim-1]
        contrib = peakDim.findFirstPeakDimContrib(resonance=resonance)
        
        if contrib:
          clearPeakDim(peakDim, contrib)
      
      self.assignment = None
  
  def clearAssignments(self):

    if self.peak:
      for peakContrib in self.peak.peakContribs:
        peakContrib.delete()
    
      for peakDim in self.peak.peakDims:
        clearPeakDim(peakDim)

  def getWindows(self):
  
    windowPanes = []
    for window in self.analysisProject.spectrumWindows:
      for windowPane in window.spectrumWindowPanes:
        for view in windowPane.spectrumWindowViews:
          for winPeakList in view.windowPeakLists:
            if winPeakList.analysisPeakList.peakList is self.peakList:
              windowPanes.append(windowPane)
              break
          
          else:
            continue
          break
    
    return windowPanes

  def updateWindows(self, obj=None):
  
    windowPanes = self.getWindows()
    
    textMatrix = []
    objectList = []
    colorMatrix = []
    
    for windowPane in windowPanes:
      
      if windowPane in self.windowPanes:
        use = 'Yes'
        #colors = [assignColor] * 3
      else:
        use = 'No'
      
      colors = [None] * 3
      
      datum = [windowPane.spectrumWindow.serial,
               getWindowPaneName(windowPane),
               use,
              ]
              
      textMatrix.append(datum)
      objectList.append(windowPane)
      colorMatrix.append(colors)
    
    self.windowMatrix.update(textMatrix=textMatrix,
                             objectList=objectList,
                             colorMatrix=colorMatrix)
  
  def updatePeaksAfter(self, object=None):


    if self.waiting:
      return
      
    if object:
      className = object.className
    
      if className == 'Peak':
        if object.peakList is not self.peakList:
          return
     
      elif className == 'PeakList':
        if object is not self.peakList:
          return
   
    self.waiting = True
    self.after_idle(self.updatePeaks)

  def updatePeaks(self):
  
    textMatrix = []
    objectList = []
    #headingList = ['#','Position','Assignment','Height','Volume']
    
    if self.peakList:
      for peak in self.peakList.sortedPeaks():
        details = peak.details
        if details and (len(details)>32):
          details = details[:32]+'...'
      
        peakDims = peak.sortedPeakDims()
        datum = [peak.serial,
                 ' '.join(['%.3f' % pd.realValue for pd in peakDims]),
                 ' '.join([pd.annotation or '-' for pd in peakDims]),
                 details,
                 getPeakHeight(peak),
                 getPeakVolume(peak),]
 
        textMatrix.append(datum)
        objectList.append(peak)
        
    
    self.waiting = False
    self.peakMatrix.update(textMatrix=textMatrix,objectList=objectList)
  
  def selectAssignment(self, resonances, row, col):
  
    self.assignment = resonances
    self.dim = None 
    
    if self.peak:
      nDims = len(self.peak.peakDims)
      
      if 0 < col <= nDims:
        self.dim = col

  
  def assignSelected(self):
  
    if self.peak:
      for score, resonances, dist in self.assignMatrix.currentObjects:
        peakContribs = [self.peak.newPeakContrib(),]
        
        for peakDim in self.peak.sortedPeakDims():
          dim = peakDim.dim
          i = dim-1
          
          if dim not in self.hydroDims:
            if dim not in self.heterDims:
              continue
          
          resonance = resonances[i]
          contrib = peakDim.findFirstPeakDimContrib(resonance=resonance)
          
          tol = self.tolEntries[i].get() or getAnalysisDataDim(peakDim.dataDim).assignTolerance
          assignResToDim(peakDim, resonance, tolerance=tol,
                         contrib=contrib, peakContribs=peakContribs)
              
  def processAll(self):
  
    if not self.peakList and self.structure:
      return
    
    peaks = self.peakList.peaks
    msg = 'Really assign all %s peaks in list to their matching assignment ' % len(peaks)
    msg += 'possibilities using current distance and shift tolerances?'
    
    if not showOkCancel('Confirm', msg, parent=self):
      return
    
    structure = self.structure
    distThreshold = self.distThresholdEntry.get() or 8.0
    dimH1, dimH2 = self.hydroDims
    dimX1, dimX2 = self.heterDims
    
    progressBar = ProgressBar(self, text='Peak progress:', progress=0,
                              total=len(peaks), title='Assigning close resonances')

    for peak in peaks:
      tolerances, assignments = self.getPossibleAssignments(peak)
  
      for shiftData in assignments:

        resonance1, valueH1, deltaH1 = shiftData[dimH1-1]
        resonance2, valueH2, deltaH2 = shiftData[dimH2-1]
        
        resonanceSet1 = resonance1.resonanceSet
        resonanceSet2 = resonance2.resonanceSet
        
        dist = None
        if resonanceSet1 and resonanceSet2:
          dist = getAtomSetsDistance(resonanceSet1.atomSets,
                                     resonanceSet2.atomSets,
                                     structure)

        if (dist is not None) and (dist>distThreshold):
          continue

        for j, peakDim in enumerate(peak.sortedPeakDims()):
          dim = peakDim.dim
          if dim not in self.hydroDims:
            if dim not in self.heterDims:
              continue

          tol = self.tolEntries[j].get() or \
                getAnalysisDataDim(peakDim.dataDim).assignTolerance

          resonance, value, delta = shiftData[j]
          assignResToDim(peakDim, resonance, tolerance=tol)
      
      progressBar.increment()
    
    progressBar.destroy()
 
  def getPossibleAssignments(self, peak=None):
  
    if not peak:
      peak = self.peak

    assignments = []
    tolerances  = []

    if self.peakList and peak:
      spectrum    = self.peakList.dataSource
      numDim      = spectrum.numDim
      doAliasing  = self.aliasSelect.get()
      hResonances = [[],[]]
      
      for dataDim in spectrum.sortedDataDims():
        if dataDim.className != 'SampledDataDim':
          i = dataDim.dim-1
          tolerance =  self.tolEntries[i].get() or getAnalysisDataDim(dataDim).assignTolerance
          tolerances.append(tolerance)
          
        else:
          tolerances.append(None)  
          
      j = 0
      for dim in self.hydroDims:
        i = dim-1
        peakDim   = peak.findFirstPeakDim(dim=dim)
        value     = peakDim.realValue
        tolerance = tolerances[i]
        dimX      = self.heterDims[j]
        shiftsH   = findMatchingPeakDimShifts(peakDim, tolerance=tolerance,
                                              aliasing=doAliasing,
                                              findAssigned=True)
        
        if dimX is not None:
          peakDimX   = peak.findFirstPeakDim(dim=dimX)
          resonances = []
          shiftsX = findMatchingPeakDimShifts(peakDimX, tolerance=tolerances[dimX-1],
                                              aliasing=doAliasing,
                                              findAssigned=True)
          
          for shiftH in shiftsH:
            resonanceH = shiftH.resonance
          
            for shiftX in shiftsX:
              resonanceX = shiftX.resonance
              
              if areResonancesBound(resonanceH, resonanceX):
                valueH = shiftH.value
                valueX = shiftX.value
                shiftDataH = (resonanceH, valueH, abs(valueH-value))
                shiftDataX = (resonanceX, valueX, abs(valueX-peakDimX.realValue))
                hResonances[j].append((i, dimX-1, shiftDataH, shiftDataX))
                break
        
        else:
          for shiftH in shiftsH:
            valueH = shiftH.value
            shiftDataH = (shiftH.resonance, valueH, abs(valueH-value))
            hResonances[j].append((i,None,shiftDataH,None))
     
        j += 1
  
 
      for iH1, iX1, shiftDataH1, shiftDataX1 in hResonances[0]:
        resonanceSet1 = shiftDataH1[0].resonanceSet
        if resonanceSet1:
          molSystem1 = resonanceSet1.findFirstAtomSet().findFirstAtom().residue.chain.molSystem
        else:
          molSystem1 = None
      
        for iH2, iX2, shiftDataH2, shiftDataX2 in hResonances[1]:
          resonanceSet2 = shiftDataH2[0].resonanceSet
          if resonanceSet2:
            molSystem2 = resonanceSet2.findFirstAtomSet().findFirstAtom().residue.chain.molSystem
          else:
            molSystem2 = None
 
          if molSystem1 and (molSystem1 is not molSystem2):
            continue
 
          assignment = [None] * numDim
          assignment[iH1] = shiftDataH1
          assignment[iH2] = shiftDataH2
          
          if iX1 is not None:
            assignment[iX1] = shiftDataX1
            
          if iX2 is not None:
            assignment[iX2] = shiftDataX2   
 
          assignments.append(assignment)
 
  
    return tolerances, assignments

  def updateWindowsAfter(self, obj):
    
    self.after_idle(self.updateWindows)
  
  def updateAssignmentsAfter(self, obj=None):

    if self.waitingAssign:
      return

    if obj:
      if obj.className == 'PeakDim':
        if obj.peak is not self.peak:
          return
      
      elif obj.className == 'PeakDimContrib':
        if obj.peakDim.peak is not self.peak:
          return
    
    self.waitingAssign = True
    self.after_idle(self.updateAssignments)

  def getHeadingList(self, numDim):
  
    dimRange = range(numDim)
    headingList  = ['#']
    headingList += ['Assign\nF%d: %s' % (i+1,self.isotopes[i]) for i in dimRange]
    headingList += ['Dist',u'Tol. Weighted\nCombined \u0394',]
    headingList += [u'\u0394F%d' % (i+1) for i in dimRange]
    headingList += ['Contrib\nFrac',]
    headingList += ['Shift\nF%d' % (i+1) for i in dimRange]
  
    tipTexts = ['Rank of assignment possibility; best matching first',]
    tip = 'The %s resonance assignment possibility for dimension %d'
    tipTexts += [ tip % (self.isotopes[i],i+1) for i in dimRange]
    tipTexts += ['Distance, in the selected structure, between atom assignment possibilities, for the NOESY linked dimensions',
                 'The combined chemical shift difference, for all dimensions & weighted by isotope tolerance, from peak to possible resonance assignment',]
    tip = 'The chemical shift difference between peak position and resonance assignment for dimension %d'
    tipTexts += [tip % (i+1) for i in dimRange]
    tip = 'The fraction of the total NOE intensity that a committed assignment contributes; assuming NOE intensity proportional to distance to the power of -6'
    tipTexts += [tip,]
    tip = 'The chemical shift of the resonance assignment possibility for dimension %d'
    tipTexts += [tip % (i+1) for i in dimRange]
  
    return headingList,  tipTexts

  def updateAssignments(self):
  
    textMatrix    = []
    objectList    = []
    colorMatrix   = []
    distThreshold = self.distThresholdEntry.get() or 8.0
    noeSum = 0.0
    
    if self.structure:
      showConn = self.strucDispPulldown.getText()
    else:
      showConn = 'None'  
    
    if showConn != 'None':  
      self.guiParent.viewStructure(self.structure)
      popup = self.guiParent.popups['view_structure']
      popup.clearConnections()
        
    if self.peak:
      numDim   = len(self.peak.peakDims)
      dimRange = range(numDim)
      headingList, tipTexts = self.getHeadingList(numDim)
      
      dimResonances = []
      for peakDim in self.peak.sortedPeakDims():
        resonances  = []
        
        for contrib in peakDim.peakDimContribs:
          resonances.append(contrib.resonance)
        
        dimResonances.append(resonances)
            
      tolerances, assignments = self.getPossibleAssignments()
      
      sortedData = []
      atoms = {}
      dimH1, dimH2 = self.hydroDims
      dimX1, dimX2 = self.heterDims
      for shiftData in assignments:
      
        resonance1, valueH1, deltaH1 = shiftData[dimH1-1]
        resonance2, valueH2, deltaH2 = shiftData[dimH2-1]
        
        if self.structure:
          resonanceSet1 = resonance1.resonanceSet
          resonanceSet2 = resonance2.resonanceSet
 
          if resonanceSet1:
            for atomSet in resonanceSet1.atomSets:
              for atom in atomSet.atoms:
                atoms[atom] = atoms.get(atom, 0) + 1
 
          if resonanceSet2:
            for atomSet in resonanceSet2.atomSets:
              for atom in atomSet.atoms:
                atoms[atom] = atoms.get(atom, 0) + 1
 
          if resonanceSet1 and resonanceSet2:
            dist = getAtomSetsDistance(resonanceSet1.atomSets,
                                       resonanceSet2.atomSets,
                                       self.structure)
          else:
            dist = None
        else:
          dist = None
        
        if (dist is not None) and (dist>distThreshold):
          continue
        
        shiftDist = 0.0
        shiftDist += (deltaH1/tolerances[dimH1-1])**2 
        shiftDist += (deltaH2/tolerances[dimH2-1])**2 
        
        if dimX1 is not None:
          bound, valueX1, deltaX1 = shiftData[dimX1-1]
          shiftDist += (deltaX1/tolerances[dimX1-1])**2 
         
        if dimX2 is not None:
          bound, valueX2, deltaX2 = shiftData[dimX2-1]
          shiftDist += (deltaX2/tolerances[dimX2-1])**2
 
        if not dist:
          score = sqrt( shiftDist )
        else:
          score = sqrt( ((dist/distThreshold)**2) + shiftDist )  
        
          nAssign = 0
          for j in dimRange:
            resonance, value, delta = shiftData[j]
            if resonance in dimResonances[j]:
              nAssign += 1
 
          if nAssign == numDim:
            noeSum += dist ** -6.0
  
        sortedData.append((score, shiftData, dist, sqrt(shiftDist)))
        
      coordAtom  = None
      if self.focusSelect.get():
        atoms = [(atoms[atom],atom) for atom in atoms.keys()]
        if atoms:
          atoms.sort()
          atom    = atoms[-1][1]
          residue = atom.residue
          chain   = residue.chain
 
          coordChain = self.structure.findFirstCoordChain(chain=chain)
 
          if coordChain:
            coordResidue = coordChain.findFirstResidue(residue=residue)
 
            if coordResidue:
              coordAtom = coordResidue.findFirstAtom(atom=atom)
      
      sortedData.sort()
      
      i = 1
      for score, shiftData, dist, shiftDist, in sortedData:
        # resonance1, valueH1, deltaH1 =  data1

        colors = [None] * (3*numDim+3)
        names  = [None] * numDim
        deltas = [None] * numDim
        values = [None] * numDim
        resonances = []
        
        for j in dimRange:
          resonance, value, delta = shiftData[j]
          resonances.append(resonance)
          deltas[j] = '%.3f' % delta
          names[j]  = makeResonanceGuiName(resonance)
          values[j] = value
          
          f = min(1.0, delta/tolerances[j])
          colors[j+numDim+3] = self.getScaledRgb(f)
          
          if resonance in dimResonances[j]:
            colors[j+1] = assignColor
         
        if dist:
          f = min(dist/distThreshold, 1.0)
          color = self.getScaledRgb(f)
          
        else:
          color = None
        
        colors[numDim+1] = color

        colors[numDim+2] = self.getScaledRgb(min(shiftDist, 1.0))        
        
        if dist is None:
          distStr = None
          contribFrac = None
        else:
          distStr = '%.3f' % dist
        
          if (None not in colors[1:numDim+1]) and noeSum and dist:
            contribFrac = '%.3f' % ((dist**-6.0)/noeSum)
          else:
            contribFrac = None
        
        datum = [i,] + names + [distStr,shiftDist] + deltas + [contribFrac] + values
 
        textMatrix.append(datum)
        objectList.append([score, resonances, dist])
        colorMatrix.append(colors)
        i += 1
    
        resonance1 = resonances[self.hydroDims[0]-1]
        resonance2 = resonances[self.hydroDims[1]-1]
   
        if showConn == 'All':
          if coordAtom:
            popup.structFrame.focusOnAtom(coordAtom)
            
          popup.showResonancesConnection(resonance1,resonance2)
                    
        elif showConn == 'Assigned':
          if (resonance1 in dimResonances[0]) and (resonance2 in dimResonances[1]):
            # Must be from same peakContrib really...
          
            if coordAtom:
              popup.structFrame.focusOnAtom(coordAtom)
              
            popup.showResonancesConnection(resonance1,resonance2)
    
    elif self.peakList:
      headingList, tipTexts = self.getHeadingList(self.peakList.dataSource.numDim)

    else:
      headingList, tipTexts = self.getHeadingList(3)
    
    self.assignMatrix.update(textMatrix=textMatrix,
                             headingList=headingList,
                             tipTexts=tipTexts,
                             objectList=objectList,
                             colorMatrix=colorMatrix)
    
    if noeSum:                      
      self.noeSumLabel.set('%.3f' % (noeSum**(-1.0/6.0)))
    else:
      self.noeSumLabel.set('')
    
    self.waitingAssign = False
  
  def getScaledRgb(self, factor):
  
    r = min(1.0,max(0.0,factor))*255
    g = 255-r
    b = 32
    
    return '#%02x%02x%02x' % (r, g, b)
  
  def predictPeaks(self):
  
    if self.peak and self.assignMatrix.objectList:
      
      analysisProject = self.analysisProject
      prevDoDetails = analysisProject.doDetailAnnotations
      analysisProject.doDetailAnnotations = True
      
      distThreshold = self.distThresholdEntry.get() or 8.0
    
      spectrum  = self.peak.peakList.dataSource
      shiftList = spectrum.experiment.shiftList
      peakList  = spectrum.findFirstPeakList(isSimulated=True, name='StructurePredict')
    
      if not peakList:
        peakList = spectrum.newPeakList(isSimulated=True, name='StructurePredict')
        analysisPeakList = getAnalysisPeakList(peakList)
        analysisPeakList.symbolStyle = '+'
        analysisPeakList.symbolColor = '#FF0000'
        analysisPeakList.textColor = '#BB0000'
      
      peakList.details = 'Structure Prediction'
       
      for score, resonances, dist in self.assignMatrix.objectList:
        if dist < distThreshold:
          numRes = len(resonances)
          shifts = [r.findFirstShift(parentList=shiftList) for r in resonances]
          values = [s.value for s in shifts]
 
          if len(values) == numRes:
            peaks = {}
            i = 0
            for resonance in resonances:
              i += 1
 
              for contrib in resonance.peakDimContribs:
                peakDim = contrib.peakDim
 
                if peakDim.dim == i:
                  peak0 = peakDim.peak
                  
                  if peak0.peakList is peakList:
                    peaks[peak0] = peaks.get(peak0, 0) + 1

            peak = pickPeak(peakList, values, unit=shiftList.unit)
            for peak0 in peaks.keys():
              if peaks[peak0] == numRes:
                peak0.delete()
 
            if dist:
              peak.setDetails(' :%.2f' % dist)
            else:
              peak.setDetails(None)
            
            for i, peakDim in enumerate(peak.sortedPeakDims()):
              assignResToDim(peakDim,resonances[i])
              
      analysisProject.doDetailAnnotations = prevDoDetails
    
  def setAppDataOptions(self):

    pass

  def getAppDataOptions(self):

    pass
  
  def destroy(self): 
  
    self.administerNotifiers(self.unregisterNotify)      
  
    BasePopup.destroy(self)


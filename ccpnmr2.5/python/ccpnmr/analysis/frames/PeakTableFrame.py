
"""
======================COPYRIGHT/LICENSE START==========================

PeakTableFrame.py: Part of the CcpNmr Analysis program

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

from memops.gui.Button          import Button
from memops.gui.ButtonList      import ButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.DataEntry       import askInteger, askString
from memops.gui.Entry           import Entry
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Frame           import Frame
from memops.gui.Label           import Label
from memops.gui.MessageReporter import showOkCancel, showWarning, showError, showMulti, showYesNo
from memops.gui.PulldownList    import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix

from ccp.util.NmrExpPrototype import longRangeTransfers
 
from ccpnmr.analysis.core import ExperimentBasic
from ccpnmr.analysis.core import PeakBasic
from ccpnmr.analysis.core import UnitConverter
from ccpnmr.analysis.core import WindowBasic 

from ccpnmr.analysis.core.AssignmentBasic import clearPeakDim, propagatePeakAssignments
from ccpnmr.analysis.core.MarkBasic       import createPeakMark
from ccpnmr.analysis.core.StructureBasic  import getAtomSetsDistance
from ccpnmr.analysis.core.Util            import getMethod

UNIT_POINTS = 'point'
UNIT_PPM = 'ppm'
UNIT_HZ = 'Hz'
UNITS = [UNIT_PPM, UNIT_HZ, UNIT_POINTS]

def _getPeakDimAnnotation(peakDim):
  # 23 May 2016: add molSystem to peak table assignment information if doMolSystemsAnnotations
  analysisProject = peakDim.root.currentAnalysisProject
  annotation = peakDim.annotation
  if analysisProject.doMolSysAnnotations:
    molSystem = None
    for contrib in peakDim.peakDimContribs:
      if contrib.resonance.resonanceSet:
        atom = contrib.resonance.resonanceSet.findFirstAtomSet().findFirstAtom()
        molSystem = atom.topObject
        break
    if molSystem:
      molSysCode = '%s:' % (molSystem.code)
      annotation = molSysCode + annotation

  return annotation

class PeakTableFrame(Frame):

  def __init__(self, parent, analysis, nmrProject=None, peakList=None,
               peaks=None, simplified=False, *args, **kw):
    
    self.nmrProject = nmrProject
    self.analysisApp = analysis
    self.analysisProject = analysis.analysisProject
    
    widget = parent
    while widget.__class__.__name__ != 'AnalysisPopup':
      widget = widget.parent
    self.topPopup = widget
    
    self.peak = None
    self.peaks = peaks or []
    self.marks = []
    self.peakList = peakList
    self.windowPane = None
    self.selectStatus = 'Any'
    self.unit = UNIT_PPM
    self.waiting = False
    self.simplified = simplified
    
    Frame.__init__(self, parent,  **kw)

    self.grid_rowconfigure(1, weight=1)
    self.grid_columnconfigure(0, weight=1)
   
    self.structure = None
    self.topFrame = topFrame = Frame(self, borderwidth=2, relief='groove')
    
    ssFrame = Frame(topFrame)
    ssFrame.expandGrid(None,4)
    
    statusLabel = Label(ssFrame, text=' Status:', grid=(0,0))
    texts = ['Any','Fully Assigned','Part Assigned',
             'Unassigned','Intra Spin System']
             #'Sequential', 'Within 4 Residues', '2 to 4 Residues', 'Long Range']
    self.statusPulldown = PulldownList(ssFrame, callback=self.changeStatus,
                                       texts=texts, index=0, grid=(0,1),
                                       tipText='Restrict peak table to only those assigned a certain way')
    structLabel = Label(ssFrame, text=' Structure:', grid=(0,2))
    self.structPulldown = PulldownList(ssFrame, grid=(0,3),
                                       callback=self.setStructure,
                                       tipText='Structure to calculate distances from')
    
    if self.nmrProject:
      frame = Frame(topFrame, grid=(0,0))
      frame.expandGrid(None,7)

      label = Label(frame, text='Peak List:', grid=(0,0))
      self.listPulldown = PulldownList(frame, callback=self.changeList, grid=(0,1),
                                       tipText='Select peak list to display peaks from')
      
      label = Label(frame, text=' Position Unit:', grid=(0,2))
      self.unitPulldown = PulldownList(frame, callback=self.changeUnit, 
                                       texts=UNITS, index=0, grid=(0,3),
                                       tipText='Choose whether to display peak positions in ppm, Hz or points')
 
      ssFrame.grid(row=1, column=0, sticky='ew')
    
    else:
      self.topFrame.expandGrid(None,10)
      if not self.simplified:
        label = Label(topFrame, text=' Position Unit:', grid=(1,9))
        self.unitPulldown = PulldownList(topFrame, callback=self.changeUnit,
                                         texts=UNITS, index=0, grid=(1,10),
                                         tipText='Choose whether to display peak positions in ppm, Hz or points')

    self.markButton  = Button(topFrame, text=' Mark Selected',
                              borderwidth=1, padx=2, pady=1,
                              command=self.makePeakMarks,
                              tipText='Put multidimensional cross-marks through selected peaks')
    self.stripButton = Button(topFrame, text=' Strip Selected ', borderwidth=1,
                              padx=2, pady=1,command=self.makePeakStrips,
                              tipText='Use the positions of the selected peaks to specify strip locations in the selected window')
    """
    self.findButton  = Button(topFrame, text=' Find Peak ',
                              borderwidth=1, padx=2, pady=1,command=self.findPeak,
                              tipText='Locate the currently selected peak in the specified window')
"""
    self.findButton  = CheckButton(topFrame, text=' Find Peak ',
                              borderwidth=1, padx=2, pady=1, callback=lambda isSelected: self.findPeak(),
                              tipText='Locate the currently selected peak in the specified window')
    self.windowLabel = Label(topFrame, text=' Window:', fg = 'grey')
    self.windowPulldown = PulldownList(topFrame, callback=self.selectWindowPane,
                                       tipText='Choose the spectrum window for locating peaks or strips')
    self.markFoundLabel = Label(topFrame, text='Mark Found')
    self.markFoundCheck = CheckButton(topFrame,
                                      tipText='Whether to put a cross-mark though peaks found in a given window')
    self.markFoundCheck.set(True)

      
    self.orthStripButton = Button(topFrame, borderwidth=1, padx=2, pady=1,
                                  text = 'Strip Locations' , command=self.gotoOrthogStrips,
                                  tipText='Use the selected peak positions to specify strip positions in a higher dimensionality window.')
    """
    self.orthNavButton   = Button(topFrame, borderwidth=1, padx=2, pady=1,
                                  text = 'Go To Position' , command=self.gotoOrthogPlane,
                                  tipText='Use the selected peak position to navigate to a  location in a higher dimensionality window')
"""
    self.orthNavButton   = CheckButton(topFrame, borderwidth=1, padx=2, pady=1,
                                  text = 'Go To Position', callback=lambda isSelected: self.gotoOrthogPlane(),
                                  tipText='Use the selected peak position to navigate to a location in a higher dimensionality window')
    self.orthWinPulldown = PulldownList(topFrame, callback=None,
                                        tipText='Choose the higher dimensionality spectrum window to use as a navigation & strip target')
 
    if self.simplified:
      self.stripButton.grid(row=0, column=2, sticky='ew')
      self.findButton.grid(row=0, column=3, sticky='ew')
      self.windowLabel.grid(row=0, column=4, sticky='w')
      self.windowPulldown.grid(row=0, column=5, sticky='w')
      self.markFoundLabel.grid(row=0, column=6, sticky='w')
      self.markFoundCheck.grid(row=0, column=7,sticky='w')

      self.markButton.grid(row=0, column=8, sticky='ew', columnspan=2)
    
    else:
 
      self.stripButton.grid(row=0, column=2, sticky='ew')
      self.findButton.grid(row=0, column=3, sticky='ew')
      self.windowLabel.grid(row=0, column=4, sticky='w')
      self.windowPulldown.grid(row=0, column=5, sticky='w')
      self.markFoundLabel.grid(row=0, column=6, sticky='w')
      self.markFoundCheck.grid(row=0, column=7,sticky='w')

      self.orthStripButton.grid(row=1, column=2, sticky='ew')
      self.orthNavButton.grid(row=1, column=3, sticky='ew')
      self.orthWinPulldown.grid(row=1, column=4, sticky='w', columnspan=2)
      self.markButton.grid(row=1, column=6, sticky='ew', columnspan=2)
   
    
    row = 0
    topFrame.grid(row=row, column=0, sticky='ew')
    
    if self.nmrProject:
      topFrame.grid_columnconfigure(1, weight=1)
      
    else:
      topFrame.grid_columnconfigure(10, weight=1)
      
    row += 1
    self.shiftEntries = []
    self.lineWidthEntries = []
    
    headingList = ['#','Position\nF1','Position\nF2',
                   'Assign\nF1', 'Assign\nF2',
                   'Height','Volume',
                   'Line Width\nF1 (Hz)','Line Width\nF2 (Hz)',
                   'Merit','Details','Fit Method','Vol. Method',]
                   
    self.scrolledMatrix = ScrolledMatrix(self,
                                         headingList=headingList,
                                         multiSelect=True,
                                         callback=self.selectCell,
                                         deleteFunc=self.peakRemove)
    self.scrolledMatrix.grid(row=row, column=0, sticky='nsew')
    
    # These are deliberately defined after the scrolled matrix; for some reason
    # defining these first (and not passing them in like usual)
    # results in a widget that is not displayable for some bizarre reason
    self.heightEntry  = FloatEntry(self,text='', returnCallback=self.setHeight, width=10, formatPlaces=9)
    self.volumeEntry  = FloatEntry(self,text='', returnCallback=self.setVolume, width=10, formatPlaces=9)
    self.meritEntry   = FloatEntry(self,text='', returnCallback=self.setMerit,  width=6)
    self.detailsEntry = Entry(self,text='', returnCallback=self.setDetails, width=16)
    
    tipTexts = ['Add a new peak, specifying its position',
                'Edit the position of the currently selected peak',
                'Move the ppm position of a peak a number of sweep withs to its correct aliased/folded position',
                'Delete the currently selected peaks',
                'Assign the dimensions of the currently selected peak',
                'Remove all assignments from the currently selected peaks',
                'Set the details field of the currently selected peaks',
                'Set the peaks selected in the table as the ones selected in the spectrum windows',
                'Show a table of resonances assigned to the selected peaks']    
    
    texts = ['Add','Edit','Unalias','Delete',
             'Assign','Deassign','Set Details',
             'Set As Current', 'Resonances']
                
    commands = [self.peakAdd,self.peakEdit,self.aliasPeak,
                self.peakRemove,self.peakAssign,self.peakDeassign,
                self.peaksDetails,self.setAsCurrent,self.showResonances]
    
    row += 1            
    self.bottomButtons1 = ButtonList(self, texts=texts, commands=commands,
                                     tipTexts=tipTexts)
    self.bottomButtons1.grid(row=row, column=0, sticky='ew')
    
    tipTexts = ['Deassign a specified dimension of the selected peaks',
                'Recalculate the center, height and line width of the selected peaks',
                'Recalculate the volume of the selected peaks',
                'Show the assignment connections of the selected peaks on the selected structure',
                'Spread the resonance assignments of the peak last selected to all selected peaks',
                'Copy the merit value of the last selected peak to all selected peaks',
                'Copy the details of the last selected peak to all selected peaks']
    
    texts = ['Deassign Dim','Recalc Fit','Recalc Volume','Show On Structure',
             'Propagate Assign','Propagate Merit','Propagate Details']
                
    commands = [self.dimDeassign, self.recalcPeakFit, self.recalcPeakVolume,
                self.showStructConnections,self.propagatePeakAssign,self.propagateMerit,
                self.propagateDetails]
                
    row += 1            
    self.bottomButtons2 = ButtonList(self, texts=texts, expands=True,
                                    commands=commands, tipTexts=tipTexts)
    self.bottomButtons2.grid(row=row, column=0, sticky='ew')
    
    if self.nmrProject:
      self.updatePeakLists()
    
    self.updateWindowLists()
    self.updateStructures()
    self.updateButtons()
    self.updatePeaksAfter()

  def propagatePeakAssign(self):
  
    peaks = list(self.scrolledMatrix.currentObjects)
    if len(peaks) > 1:
      propagatePeakAssignments(peaks, warnUnalias=True)
    else:
      msg = 'More than one peak must be selected in the table for assignment propagation.'
      showWarning('Warning', msg, parent=self)

  def propagateMerit(self):
  
    peaks = list(self.scrolledMatrix.currentObjects)
    if len(peaks) > 1:
      merit = self.peak.figOfMerit
      for peak in peaks:
        peak.figOfMerit = merit

    else:
      msg = 'More than one peak must be selected in the table for merit propagation.'
      showWarning('Warning', msg, parent=self)

  def propagateDetails(self):
  
    peaks = list(self.scrolledMatrix.currentObjects)
    if len(peaks) > 1:
      details = self.peak.details
      for peak in peaks:
        peak.details = details

    else:
      msg = 'More than one peak must be selected in the table for details propagation.'
      showWarning('Warning', msg, parent=self)

  def setupOrthogWinPositions(self, *event):
    
    if not self.peak:
      self.orthWinPulldown.clear()
      return
    
    selected = self.orthWinPulldown.getText()
    
    windowPane0 = None
    spectrum = self.peak.peakList.dataSource
    nDims = spectrum.numDim

    for window in self.analysisProject.spectrumWindows:
      for windowPane in window.spectrumWindowPanes:
        if len(windowPane.axisPanels) == nDims:
          if WindowBasic.isSpectrumInWindowPane(windowPane, spectrum):
            windowPane0 = windowPane
            break
 
    
    if not windowPane0:
      for window in self.analysisProject.spectrumWindows:
        for windowPane in window.spectrumWindowPanes:
          if WindowBasic.isSpectrumInWindowPane(windowPane, spectrum):
            windowPane0 = windowPane
            break
        
    orthogonalDict = {}
    if windowPane0:
      xyz = []
      for peak in self.scrolledMatrix.currentObjects:
        xyz.append( WindowBasic.getPeakWindowPosition(peak, windowPane0) )

      windowZPlanes = WindowBasic.findOrthogonalWindows(windowPane0, xyz, excludeQuery=False)
      for windowPosition in windowZPlanes:
        if windowPosition:
          key, windowPane, xyz = windowPosition
          orthogonalDict[key] = (windowPane,xyz)
 
    locations = []
    index = 0
    xyPlaneKeys = orthogonalDict.keys()
    xyPlaneKeys.sort()
    
    if xyPlaneKeys:
      if selected in xyPlaneKeys:
        index = xyPlaneKeys.index(selected)
      locations = [orthogonalDict[key] for key in xyPlaneKeys]  
    
    self.orthWinPulldown.setup(xyPlaneKeys, locations, index)
  
  def gotoOrthogStrips(self):
  
    location = self.orthWinPulldown.getObject()
    if self.peak and location:
      (windowPane, positions) = location
      window = windowPane.spectrumWindow
      
      orthoAxis = 'y'
      if window.stripAxis == 'y':
        orthoAxis = 'x'

      N = len(positions)
      if N > 10:
        msg = '%d positions selected. Really make %d strips in window %s?' % (N,N,window.name)
        if not showOkCancel('Warning',msg, parent=self):
          return

      M = 0.0
      orthoPos = 0.0
      for xyz in positions:
        if xyz.get(orthoAxis):
          orthoPos += xyz.get(orthoAxis)
          M += 1.0

      if M > 0.0:
        orthoPos /= M
        for xyz in positions:
          xyz[orthoAxis] = orthoPos
 
      if positions:
        WindowBasic.displayStrips(self.analysisApp, positions,
                                  orthoPositions=None, spectrum=None,
                                  windowPane=windowPane)
        self.update_idletasks()
        WindowBasic.displayStrips(self.analysisApp, positions,
                                  orthoPositions=None, spectrum=None,
                                  windowPane=windowPane)
        

  def gotoOrthogPlane(self):
  
    location = self.orthWinPulldown.getObject()
    if self.peak and location:
      try:
        index = self.scrolledMatrix.currentObjects.index(self.peak)
      except ValueError, e:
        # this happens if self.peak is not on list
        # that should not be the case but sometimes is
        return
      (windowPane, positions) = location
      
      windowFrame = windowPane.getWindowFrame()
      windowFrame.gotoPosition(position=positions[index])

  def showStructConnections(self):
  
    peaks = list(self.scrolledMatrix.currentObjects)
    if peaks:
      self.analysisApp.viewStructure(self.structure)
      popup = self.analysisApp.popups['view_structure']
      popup.clearConnections()
      
      for peak in peaks:
        popup.showPeakConnection(peak)

  def showAllStructConnections(self):
  
    if self.peakList or self.peaks:
      if self.peakList and self.nmrProject:
        peaks = self.peakList.peaks
      else:
        peaks = self.peaks  
      
      self.analysisApp.viewStructure(self.structure)
      popup = self.analysisApp.popups['view_structure']
      popup.clearConnections()
      for peak in peaks:
        popup.showPeakConnection(peak)

  def updateStructures(self, *obj):
    
    index = 0
    names = ['<None>',]
    structures = self.getStructures()
    
    for structure in structures[1:]:
      name = '%s:%d' % (structure.molSystem.code, structure.ensembleId)
      names.append(name)
  
    self.structPulldown.setup(names, structures, index)


  def getStructures(self):

    if self.peakList and self.nmrProject:
      peakLists = [self.peakList,]
    else:
      peakLists = set()
      for peak in self.peaks:
        peakLists.add(peak.peakList)
    
      peakLists = list(peakLists)
    
    molSystems = set()
    for peakList in peakLists:  
      experiment = peakList.dataSource.experiment
      
      for molSystem in experiment.sortedMolSystems():
        molSystems.add(molSystem)


    if not molSystems:
      molSystems = self.analysisApp.project.sortedMolSystems()  
    
    structures = []  
    for molSystem in molSystems:
      code = molSystem.code
      
      for structure in molSystem.structureEnsembles:
        structures.append((code, structure.ensembleId, structure))

    structures.sort()
    structures = [x[2] for x in structures]

    return [None,] + structures


  def setStructure(self, structure):
    
    if structure is not self.structure:
      self.structure = structure
      self.updateButtons()
      self.updatePeaksAfter()


  def recalcPeakFit(self):

    texts = list(PeakBasic.PEAK_FIT_METHODS)
    texts.append('cancel')
    objects = texts
    fitMethod = showMulti('Peak fit method', 'Select peak fit method:', texts, objects, self)
    if not fitMethod or fitMethod == 'cancel':
      return

    doFit = showYesNo('Fit position?', 'Should position be fit (otherwise only height and linewidth are)?', parent=self)

    recalcPeakFit = PeakBasic.setupPeakFit
    
    for peak in self.scrolledMatrix.currentObjects:
      recalcPeakFit(peak, fitMethod, doFit)

  def recalcPeakVolume(self):

    texts = list(PeakBasic.PEAK_VOLUME_METHODS)
    texts.append('cancel')
    objects = texts
    volumeMethod = showMulti('Peak volume method', 'Select peak volume method:', texts, objects, self)
    if not volumeMethod or volumeMethod == 'cancel':
      return

    recalcPeakVolume = PeakBasic.setupPeakVolume
    
    for peak in self.scrolledMatrix.currentObjects:
      recalcPeakVolume(peak, volumeMethod)

  def aliasPeak(self):
  
    if self.peak:
      peaks = self.scrolledMatrix.currentObjects
      if len(peaks) > 1:
        message = 'Propagate unaliasing from primary peak %d to all selected peaks?\n' % self.peak.serial
 
        for peakDim in self.peak.sortedPeakDims():
          txt = 'Dimension %d: Num aliasing %d (Ref peak position %f)\n'
          message += txt % (peakDim.dim,peakDim.numAliasing, peakDim.value)
 
        if showOkCancel('Confirm',message, parent=self):
          PeakBasic.propagatePeakUnaliasing(self.peak, peaks)
 
      else:
        self.analysisApp.editPeakAliasing(peak=self.peak)

  def changeStatus(self, key):
  
    self.selectStatus = key
    
    self.updatePeaksAfter()

  def makePeakMarks(self):
  
    if self.peak:
      peaks = self.scrolledMatrix.currentObjects
      
      for mark in self.marks:
        if not mark.isDeleted:
          mark.delete()
      
      self.marks = [createPeakMark(peak, lineWidth=2.0, remove=False) for peak in peaks]

  def makePeakStrips(self):
  
    if self.peak and self.windowPane:
      # window is None means new window
      peaks = self.scrolledMatrix.currentObjects
      WindowBasic.displayPeakStrips(self.analysisApp, peaks, self.windowPane)
      self.update_idletasks()
      WindowBasic.displayPeakStrips(self.analysisApp, peaks, self.windowPane)
      

  def findPeak(self):

    if self.peak and self.windowPane:
      if self.markFoundCheck.get():
        createPeakMark(self.peak, lineWidth=2.0)
      
      windowFrame = self.windowPane.getWindowFrame()
      windowFrame.gotoPeak(self.peak)

  def selectWindowPane(self, windowPane):
  
    if windowPane is not self.windowPane:
      self.windowPane = windowPane

  def updateWindows(self):
  
    index = 0
    windowPane = self.windowPane
    windowPanes = []
    names = []
    
    peakList = None
    if self.peakList:
      peakList = self.peakList
    elif self.peak:
      peakList = self.peak.peakList
      
    if peakList:
      spectrum   = peakList.dataSource
      tryWindows = WindowBasic.getActiveWindows(peakList.root)
      windowData = []
      getName = WindowBasic.getWindowPaneName
      
      for window in tryWindows:
        for windowPane0 in window.spectrumWindowPanes:
          if WindowBasic.isSpectrumInWindowPane(windowPane0, spectrum):
            windowData.append( (getName(windowPane0), windowPane0) )
      
      windowData.sort()
      names = [x[0] for x in windowData]
      windowPanes = [x[1] for x in windowData]
    
    if windowPanes:
      if windowPane not in windowPanes:
        windowPane = windowPanes[0]
        
      index = windowPanes.index(windowPane)
        
    else:
      windowPane = None
    
    self.selectWindowPane(windowPane)
    
    self.windowPulldown.setup(names, windowPanes, index)

  def updateWindowListsAfter(self, *event):

    self.after_idle(self.updateWindowLists)

  def updateWindowLists(self, *event):

    self.updateWindows()
    self.setupOrthogWinPositions()

  def selectCell(self, peak, row, col):
 
    self.peak = peak
    if peak.peakList is not self.peakList:
      self.peakList = peak.peakList
      self.updateWindowLists()
    
    else:
      self.setupOrthogWinPositions()
      
    self.after(500, self.updateButtons)
  
    if self.findButton.get():
      self.findPeak()

    if self.orthNavButton.get():
      self.gotoOrthogPlane()

  def setShift(self, index, event):

    value = self.shiftEntries[index].get()
    if value:
      value = float(value)
    else:
      return
    
    peakDim = self.peak.findFirstPeakDim(dim=index+1)
    dataDimRef = peakDim.dataDimRef
    pnt2ppm = UnitConverter.pnt2ppm
    unit = self.unit
    
    if not dataDimRef:
      peakDim.position = value
    elif unit == UNIT_POINTS:
      peakDim.value = pnt2ppm(value, dataDimRef)
    elif unit == UNIT_PPM:
      peakDim.value = value
    else:  # unit == UNIT_HZ
      peakDim.value = value / dataDimRef.expDimRef.sf
    
    # self.updatePeaksAfter()

  def getShift(self, index, peak):

    if not peak:
      return
    
    peakDims = peak.sortedPeakDims()
    if index >= len(peakDims):
      self.scrolledMatrix.keyPressEscape()
      return
    
    peakDim = peakDims[index]
    dataDimRef = peakDim.dataDimRef
    ppm2pnt = UnitConverter.ppm2pnt
    unit = self.unit

    if not dataDimRef:
      position = peakDim.position
    elif unit == UNIT_POINTS:
      position = ppm2pnt(peakDim.value, dataDimRef)
    elif unit == UNIT_PPM:
      position = peakDim.value
    else:  # unit == UNIT_HZ
      position = peakDim.value * dataDimRef.expDimRef.sf

    self.shiftEntries[index].set(position)      

  def setLineWidth(self, index, event):

    value = self.lineWidthEntries[index].get()
    if value is not None:
      value = float(value)
    else:
      return
    
    peakDim = self.peak.findFirstPeakDim(dim=index+1)
    if peakDim:
      peakDim.lineWidth = value
    
    # self.updatePeaksAfter()

  def getLineWidth(self, index, peak):

    if not peak:
      return
    
    peakDims = peak.sortedPeakDims()
    if index >= len(peakDims):
      self.scrolledMatrix.keyPressEscape()
      return
    
    peakDim = peakDims[index]
    self.lineWidthEntries[index].set(peakDim.lineWidth)

  def setDetails(self, event):

    text = self.detailsEntry.get()
    if self.peak:
      self.peak.details = text or None
      self.updatePeaksAfter()
 
  def getDetails(self, peak):
    
    if peak and peak.details:
      self.detailsEntry.set(peak.details)
  
  def setMerit(self, event):
  
    value = self.meritEntry.get()
    if value is not None:
      v = max(min(1,value),0)
      self.peak.figOfMerit = v
      self.updatePeaksAfter()
  
  def getMerit(self, peak):

    if peak:
      self.meritEntry.set(peak.figOfMerit)

  def setHeight(self, event):
  
    value = self.heightEntry.get()
    if value:
      PeakBasic.setManualPeakIntensity(self.peak, value, intensityType='height')
 
  def getHeight(self, peak):

    if peak and peak.peakIntensities:
      intensity = peak.findFirstPeakIntensity(intensityType='height')
      value = intensity and intensity.value
      self.heightEntry.set(value)
  
  def setVolume(self, event):

    value = self.volumeEntry.get()
    if value:
      PeakBasic.setManualPeakIntensity(self.peak, value, intensityType='volume')
  
  def getVolume(self, peak):

    if peak and peak.peakIntensities:
      intensity = peak.findFirstPeakIntensity(intensityType='volume')
      value = intensity and intensity.value
      self.volumeEntry.set(value)
  
  def setFitMethod(self, event):

    method = self.fitMethodPulldown.getText()
    PeakBasic.setupPeakFit(self.peak, method)

  def getFitMethod(self, peak):

    if peak:
      fitMethod = peak.fitMethod
      if fitMethod:
        parameter = fitMethod.findFirstParameter(name='method')
        if parameter:
          method = parameter.value
          if method in PeakBasic.PEAK_FIT_METHODS:
            self.fitMethodPulldown.set(method)

  def setVolMethod(self, event):

    method = self.volMethodPulldown.getText()
    PeakBasic.setupPeakVolume(self.peak, method)

  def getVolMethod(self, peak):

    if peak:
      volumeIntensity = peak.findFirstPeakIntensity(intensityType='volume')
      if volumeIntensity:
        parameter = volumeIntensity.method.findFirstParameter(name = 'method')
        if parameter:
          method = parameter.value
          if method in PeakBasic.PEAK_VOLUME_METHODS:
            self.volMethodPulldown.set(method)

  def updatePeakLists(self):
  
    index = 0
    names = []    
    peakLists = self.getPeakLists()
    peakList = self.peakList
    
    if peakLists:
      if peakList not in peakLists:
        peakList = peakLists[0]
      
      for peakList0 in peakLists:
        spectrum = peakList0.dataSource
        experiment = spectrum.experiment
        names.append( '%s:%s:%d' % (experiment.name, spectrum.name, peakList0.serial) ) 
       
      index = peakLists.index(peakList)
      
    else:
      peakList = None  
  
    self.changeList(peakList)
  
    self.listPulldown.setup(names, peakLists, index)
    

  def updatePeakListsAfter(self, obj):
  
    if self.nmrProject:
      self.after_idle(self.updatePeakLists)
  
  def changeUnit(self, unit):
  
    if unit is not self.unit:
      self.unit = unit     
      self.updatePeaksAfter()
     
  def changeList(self, peakList):
    
    if peakList is not self.peakList:  
      self.peakList = peakList
      self.peak = None
      self.updatePeaksAfter()
      self.updateWindowListsAfter()

  def getDefaultPeakList(self):
  
    if self.peak:
      return self.peak.peakList
  
    peakList = None
    if self.nmrProject.experiments:
      expts = self.nmrProject.sortedExperiments()
    else:
      showError('PeakList Editor Error', 'Project has no experiments', parent=self)
      return None
 
    specs = None
    for expt in expts:
      if expt.dataSources:
        specs = expt.sortedDataSources()
        break
 
    if specs is None:
      showError('PeakList Editor Error', 'Experiments have no spectra', parent=self)
      return None

    spec = specs[0]
 
    peakList = spec.activePeakList
    
    if not peakList:
      peakList = spec.findFirstPeakList()
      
      if not peakList:
        peakList = spec.newPeakList(details = 'Default list')
    
    return peakList

  def updateButtons(self):

    buttons = list(self.bottomButtons1.buttons)
    buttons += self.bottomButtons2.buttons
    n = len(buttons)
    
    if self.peak:
      for i in range(1,n):
        buttons[i].enable()
        
      if not self.structure:
        buttons[12].disable()
        
      self.windowLabel.config(fg = 'black')
      #self.windowPulldown.activate()

    else:
      for i in range(1,n):
        buttons[i].disable()
      
  def update(self, peakList=None, peaks=None):

    if peaks:
      self.peaks = peaks

    if self.nmrProject:
      if peakList:
        self.peakList = peakList

      if not self.peakList:
        self.peakList = self.getDefaultPeakList()
 
      self.updatePeakLists()
      self.updateStructures()
    
    self.updateButtons()
    self.updateWindowLists()
    self.updatePeaksAfter()
    
  def getPeakLists(self):
    
    peakLists = []
    for expt in self.nmrProject.sortedExperiments():
      expName = expt.name
      
      for spec in expt.sortedDataSources():
        if spec.dataType == 'processed':
          specName = spec.name
          
          if not spec.peakLists:
            # wb104: 21 May 2012: if opening a new spectrum then this code
            # is called on that new spectrum but it is not ready yet for
            # a new peak list, and Analysis.py should always create a peak list
            pass
            #spec.newPeakList(details='Default List')
        
          for peakList in spec.sortedPeakLists():
            name = '%s:%s:%d' % (expName, specName, peakList.serial)
            peakLists.append((name, peakList))
  
    peakLists.sort()
  
    return [x[1] for x in peakLists]
  

  def intensityUpdateAfter(self, peakIntensity):
    
    if self.waiting:
      return
    
    if peakIntensity:
      peak = peakIntensity.peak
      self.updatePeaksAfter(peak)

  def contribUpdateAfter(self, contrib):

    if self.waiting:
      return
  
    if contrib:
      peak = contrib.peakDim.peak
      self.updatePeaksAfter(peak)

  def peakDimUpdateAfter(self, peakDim):
 
    if self.waiting:
      return
 
    if peakDim:
      peak = peakDim.peak
      self.updatePeaksAfter(peak)
  
  def updatePeaksAfter(self, peak=None):
  
    if self.waiting:
      return

    if self.nmrProject:
      if (peak is None) or (peak.peakList is self.peakList):
        self.waiting = True
        self.after_idle(self.updatePeaksFromList)
         
    
    else:
      if (peak is None) or (peak in self.peaks):
        self.waiting = True
        self.after_idle(self.updatePeaksFromSelection)
  
  def getDistanceDimensions(self, spectrum):
    
    distDims  = []
    for expTransfer in spectrum.experiment.expTransfers:
      if expTransfer.transferType in longRangeTransfers:
        for expDimRef in expTransfer.expDimRefs:
          dataDim = expDimRef.expDim.findFirstDataDim(dataSource=spectrum)
 
          if dataDim:
            distDims.append(dataDim.dim-1)

        break

    if len(distDims) < 2:
      distDims = []
   
    
    return distDims
    
  def updatePeaksFromList(self):
  
    if not self.peakList:
      self.scrolledMatrix.update(objectList=[],
                                 textMatrix=[])
      return                          
    
    unit = self.unit
    ppm2pnt = UnitConverter.ppm2pnt
    structure    = self.structure
    selectStatus = self.selectStatus
    spectrum     = self.peakList.dataSource
    nDim         = spectrum.numDim
    dataDims     = spectrum.sortedDataDims()
    notSimplified = not self.simplified
 
    distDims  = []
    if structure:
      distDims = self.getDistanceDimensions(spectrum)
    
    headingList, tipTexts = self.getHeadings(nDim)
        
    sampledDims  = {}
    for i in range(nDim):
      dataDim = dataDims[i]

      if dataDim.className == 'SampledDataDim':
        headingList[i+1] = 'Sampled\n%s' % dataDim.conditionVaried
        sampledDims[i] = dataDim

    if distDims:
      headingList.insert((2*nDim)+1, 'Dist.')
 
 
    nCols = len(headingList)
    textMatrix = []
    objects = []
    
    peaks = self.peakList.sortedPeaks()
    
    for peak in peaks:
      peakDims = peak.sortedPeakDims()
      assignedDims = 0
            
      for peakDim in peakDims:
        if len(peakDim.peakDimContribs) > 0:
          assignedDims +=1
                
      if selectStatus == 'Fully Assigned':
        if assignedDims != len(peakDims):
          continue
 
      elif selectStatus == 'Part Assigned':
        if assignedDims == len(peakDims) or assignedDims == 0:
          continue
  
      elif selectStatus == 'Unassigned':
        if assignedDims != 0:
          continue
  
      elif selectStatus == 'Intra Spin System':
        if assignedDims == nDim:
          n = 0
          spinSystems = {}
          for peakDim in peakDims:
            for contrib in peakDim.peakDimContribs:
              spinSystem = contrib.resonance.resonanceGroup
              if spinSystem:
                if spinSystems.get(spinSystem) is None:
                  spinSystems[spinSystem] = {}
                
                spinSystems[spinSystem][peakDim] = None

          found = False  
          for spinSystem in spinSystems.keys():
            if len(spinSystems[spinSystem].keys()) == nDim:
              found = True
              break
          
          if not found:
            continue    

      """
      elif selectStatus == 'Sequential':
        pass
      elif selectStatus == 'Within 4 Residues':
        pass
      elif selectStatus == '2 to 4 Residues':
        pass
      elif selectStatus == 'Long Range':
        pass
"""
 
      datum = [None] * (2*nDim+1)
      datum[0] = peak.serial

      lineWidths = []
      for i, peakDim in enumerate(peakDims):
        
        sampledDim = sampledDims.get(i)
        if peakDim.position is None:
          value = "NOT SET"

        elif sampledDim:
          value = sampledDim.pointValues[int(peakDim.position)-1]
        
        elif unit == UNIT_PPM:
          value = peakDim.value
        
        elif unit == UNIT_HZ:
          value = peakDim.value*peakDim.dataDimRef.expDimRef.sf
        
        else:
          value = ppm2pnt(peakDim.value, peakDim.dataDimRef)
 
        #annotation = peakDim.annotation
        annotation = _getPeakDimAnnotation(peakDim)
        if annotation:
          annotation = ' ' + annotation

        datum[i+1] = value
        datum[i+1+nDim] = annotation
        lineWidths.append(peakDim.lineWidth)

      if distDims:
        atomSets = []
        peakDimA = peakDims[distDims[0]]
        peakDimB = peakDims[distDims[1]]
        
        
        for peakDim in (peakDimA, peakDimB):
          assignedSets = []
          
          for contrib in peakDim.peakDimContribs:
            resonanceSet = contrib.resonance.resonanceSet
            
            if resonanceSet:
              assignedSets.extend( list(resonanceSet.atomSets) )
              
          atomSets.append( assignedSets )
 
        if atomSets and atomSets[0] and atomSets[1]:
          dist = getAtomSetsDistance(atomSets[0], atomSets[1], structure, method='min')
        else:
          dist = None
 
        datum.append(dist)

      heightIntensity = peak.findFirstPeakIntensity(intensityType='height')
      volumeIntensity = peak.findFirstPeakIntensity(intensityType='volume')

      height = None
      if heightIntensity:
        height = heightIntensity.value

      volume = None
      volMethod = None
      if volumeIntensity:
        volume = volumeIntensity.value
        parameter = volumeIntensity.method.findFirstParameter(name = 'method')
        if parameter:
          volMethod = parameter.value
 
      param = None
      fitMethod = peak.fitMethod
      if fitMethod:
        parameter = fitMethod.findFirstParameter(name = 'method')
        if parameter:
          param = parameter.value

      datum.extend([height, volume,])
      if notSimplified:
        datum.extend(lineWidths)
        datum.extend([peak.figOfMerit,
                      peak.details,
                      param, volMethod,])

      else:      
        datum.extend([peak.figOfMerit,
                      peak.details])
 
      textMatrix.append( datum )
      objects.append(peak)

    if self.peak:
      if self.peak not in objects:
        self.peak = None
   
        
    (editWidgets,editGetCallbacks,editSetCallbacks) = self.getEditData(nDim)
    self.scrolledMatrix.update(editWidgets=editWidgets,
                               editGetCallbacks=editGetCallbacks,
                               editSetCallbacks=editSetCallbacks,
                               headingList=headingList,
                               objectList=objects,
                               tipTexts=tipTexts,
                               textMatrix=textMatrix)

    self.waiting = False

  def updatePeaksFromSelection(self):

    self.peaks = [pk for pk in self.peaks if not pk.isDeleted]
    
    textMatrix = []
    peaks = self.peaks
    
    n = 0
    for peak in peaks:
      m = peak.peakList.dataSource.numDim
      if m > n:
        n = m
 
    simplified = self.simplified
    getPeakData = self.getPeakData
    textMatrix = [getPeakData(peak,n, simplified) for peak in peaks]

      
    (editWidgets,editGetCallbacks,editSetCallbacks) = self.getEditData(n)

    headingList, tipTexts = self.getHeadings(n)
    self.scrolledMatrix.update(objectList=peaks,
                               textMatrix=textMatrix,
                               editWidgets=editWidgets,
                               editGetCallbacks=editGetCallbacks,
                               editSetCallbacks=editSetCallbacks,
                               tipTexts=tipTexts,
                               headingList=headingList)
                               
    self.scrolledMatrix.refreshScrollbars()

    if len(peaks) == 1:
      self.scrolledMatrix.selectObject(peaks[0])

    self.waiting = False
  
  def getEditData(self, n):
  
    for widget in self.shiftEntries:
      widget.destroy()

    self.shiftEntries = []
    
    if self.nmrProject:
      editWidgets      = [None,]
      editGetCallbacks = [None,]
      editSetCallbacks = [None,]
    
    else:
      editWidgets      = [None,None,None,]
      editGetCallbacks = [None,None,None,]
      editSetCallbacks = [None,None,None,]

    for i in range(n):
      widget = FloatEntry(self,text='', returnCallback=lambda p=self, d=i: self.setShift(d,p), width=8 )
      self.shiftEntries.append(widget)
      
      editWidgets.append(widget)   
      editGetCallbacks.append(lambda p=self, d=i: self.getShift(d,p))  
      editSetCallbacks.append(lambda p=self, d=i: self.setShift(d,p))  
   
    for i in range(n):
      editWidgets.append(None)
      editGetCallbacks.append(self.peakAssign)
      editSetCallbacks.append(None)
    
    editWidgets += [self.heightEntry,self.volumeEntry]
    editGetCallbacks += [self.getHeight,self.getVolume]
    editSetCallbacks += [self.setHeight,self.setVolume]

    if not self.simplified:
      for widget in self.lineWidthEntries:
        widget.destroy()

      self.lineWidthEntries = []
      for i in range(n):
        widget = FloatEntry(self,text='', returnCallback=lambda p=self, d=i: self.setLineWidth(d,p), width=8 )
        self.lineWidthEntries.append(widget)

        editWidgets.append(widget)   
        editGetCallbacks.append(lambda p=self, d=i: self.getLineWidth(d,p))  
        editSetCallbacks.append(lambda p=self, d=i: self.setLineWidth(d,p))  

    editWidgets += [self.meritEntry,self.detailsEntry]
    editGetCallbacks += [self.getMerit,self.getDetails]
    editSetCallbacks += [self.setMerit,self.setDetails]
  
    if not self.simplified:
      self.fitMethodPulldown = PulldownList(self, callback=self.setFitMethod, texts=PeakBasic.PEAK_FIT_METHODS)
      editWidgets.append(self.fitMethodPulldown)
      editGetCallbacks.append(self.getFitMethod)
      editSetCallbacks.append(self.setFitMethod)

      self.volMethodPulldown = PulldownList(self, callback=self.setVolMethod, texts=PeakBasic.PEAK_VOLUME_METHODS)
      editWidgets.append(self.volMethodPulldown)
      editGetCallbacks.append(self.getVolMethod)
      editSetCallbacks.append(self.setVolMethod)

    return (editWidgets,editGetCallbacks,editSetCallbacks)

  def getPeakData(self, peak, n, simplified=False):
  
    unit = self.unit
    ppm2pnt = UnitConverter.ppm2pnt
    dims = peak.sortedPeakDims()
    peakList = peak.peakList
    spectrum = peakList.dataSource
    experiment = spectrum.experiment
    
    data = ['%s:%s' % (experiment.name,spectrum.name),peakList.serial,peak.serial ]
 
    linewidths = []
    for i in range(n):
      if i >= len(dims):
        data.append('N/A')
        linewidths.append(None)
      else:
        peakDim = dims[i]
        dataDimRef = peakDim.dataDimRef
        linewidths.append(peakDim.lineWidth)

        if not dataDimRef:
          dataDim = spectrum.findFirstDataDim(dim=peakDim.dim)
          
          if dataDim.className == 'SampledDataDim':
            value = dataDim.pointValues[int(peakDim.position)-1]
          else:  
            value = peakDim.position
            
        elif unit == UNIT_PPM:
          value = peakDim.value
        
        elif unit == UNIT_HZ:
          value = peakDim.value*peakDim.dataDimRef.expDimRef.sf
        
        else:
          value = ppm2pnt(peakDim.value, dataDimRef)
        
        data.append(value)
 
    for i in range(n):
      if i >= len(dims):
        data.append(None)
      else:
        #data.append( dims[i].annotation )
        annotation = _getPeakDimAnnotation(dims[i])
        data.append(annotation)

    heightIntensity = peak.findFirstPeakIntensity(intensityType='height')
    volumeIntensity = peak.findFirstPeakIntensity(intensityType='volume')

    if heightIntensity:
      hVal = heightIntensity.value 
    else:
      hVal = None
      
    vVal = None
    volMethod = None
    if volumeIntensity:
      vVal = volumeIntensity.value
      parameter = volumeIntensity.method.findFirstParameter(name = 'method')
      if parameter:
        volMethod = parameter.value

    param = None
    fitMethod = peak.fitMethod
    if fitMethod:
      parameter = fitMethod.findFirstParameter(name = 'method')
      if parameter:
        param = parameter.value

    details = peak.details
    
    if details and len(details) > 16:
      details = '%0.16s' % details
    
    data += [hVal, vVal]
    if simplified:
      data += [peak.figOfMerit, details]
    
    else:
      data += linewidths
      data += [peak.figOfMerit, details, param, volMethod]
    
    return data
  
  def getHeadings(self, n):
    
    if self.nmrProject:
      dataList = ['#',]
      tipTexts = ['Row number']
    else:
      dataList = ['Spectrum','List','Peak']
      tipTexts = ['Experiment:spectrum of peak',
                  'Peak list number',
                  'Peak serial number']
    
    
    for i in range(n):
      dataList.append('Position\nF%d' % (i+1))
      tipTexts.append('Peak location in F%d dimension' % (i+1))
   
    for i in range(n):
      dataList.append('Assign\nF%d' % (i+1))
      tipTexts.append('Resonance assignments of peak in F%d dimension' % (i+1))
      
    dataList += ['Height', 'Volume',]
    tipTexts += ['Magnitude of spectrum intensity at peak center (interpolated), unless user edited',
                 'Integral of spectrum intensity around peak location, according to chosen volume method']
       
    if not self.simplified:
      for i in range(n):
        dataList.append('Line Width\nF%d (Hz)' % (i+1))
        tipTexts.append('The linewidth in Hz of peak in F%d dimension' % (i+1))
        
    dataList += ['Merit', 'Details']
    tipTexts += ['Figure of merit value for peak; zero: "bad" one: "good"',
                 'User editable textual comment for peak']
    
    if not self.simplified:
      dataList += ['Fit Method', 'Vol. Method',]
      tipTexts += ['Method used to fit the peak location',
                   'Method used to calculate the peak volume integral']
    
    return dataList, tipTexts
         
  def peakAssign(self, *event):
  
    if self.peak:
      self.analysisApp.assignmentPanel()
      popup = self.analysisApp.popups['edit_assignment']
      popup.update(self.peak)
         
         
  def dimDeassign(self):
  
    peaks = list(self.scrolledMatrix.currentObjects)
    if peaks:
      dims = [i+1 for i in range(len(self.peak.peakDims))]
      opts = ','.join([str(d) for d in dims])
      msg = 'Clearing assignments on %d selected peaks\n' % len(peaks)
      msg += 'Enter dimension number to deassign [%s]' % opts
      dim = askInteger('Dimension Number', msg, 1, parent=self)

      if dim is None:
        return

      if dim not in dims:
        msg = '%d is not a valid dimension number '
        showWarning('Dim deassign failure', msg % dim, parent=self)
        return
 
      for peak in peaks:
        peakDim = peak.findFirstPeakDim(dim=dim)
        
        if peakDim:
          clearPeakDim(peakDim)
          
  def peaksDetails(self):
  
    peaks = list(self.scrolledMatrix.currentObjects)
    if peaks:
      text = askString('Enter text:', '', parent=self)

      if text:
        text = text.strip()
      else:
        text = None
      
      for peak in peaks:
        peak.details = text
        
  def showResonances(self):
  
    resonances = set()
    peaks = list(self.scrolledMatrix.currentObjects)
    
    for peak in peaks:
      for peakDim in peak.peakDims:
        for contrib in peakDim.peakDimContribs:
          if contrib.peakDimComponent:
            continue
          
          resonances.add(contrib.resonance) 
            
    if resonances:
      self.analysisApp.viewSelectedResonances(resonances)

  def setAsCurrent(self):

    peaks = list(self.scrolledMatrix.currentObjects)
    self.analysisApp.setCurrentPeaks(peaks)

  def peakDeassign(self):
  
    peaks = list(self.scrolledMatrix.currentObjects)
    if peaks:
      msg = 'Clear assignments for %d selected peaks?' 
      if showOkCancel('Warning', msg % len(peaks), parent=self):
        PeakBasic.removePeaksAssignment(peaks)
           
         
  def peakEdit(self):
  
    if self.peak:
      peak = self.peak
      self.analysisApp.editPeak(peak=peak)


  def peakAdd(self):
  
    if self.peakList and self.nmrProject:
      peakList = self.peakList
    elif self.peak:  
      peakList = self.peak.peakList
    else:
      peakList = None
  
    if peakList:
      self.analysisApp.editPeak(peak=None,peakList=peakList)
    else:
      msg = 'No peak or peak list selected'
      showWarning('Failure', msg, parent=self)
      
  def peakRemove(self, *event):

    if self.scrolledMatrix.currentObjects:
      # 3 Sep 2015: changed the below because the PeakBasic call
      # does not update the current objects in Analysis.py
      # and this crashes the programme
      ###PeakBasic.deletePeak(self.scrolledMatrix.currentObjects)
      self.topPopup.queryDeletePeaks(self.scrolledMatrix.currentObjects, self)
   
    """
   
      done =
      if done:
        self.updatePeaksAfter()
        for peak in self.scrolledMatrix.currentObjects:
          if peak in self.analysisApp.currentPeaks:
            self.analysisApp.removeSelected(peak)
 
      return True
 
    else:
      return False
    """

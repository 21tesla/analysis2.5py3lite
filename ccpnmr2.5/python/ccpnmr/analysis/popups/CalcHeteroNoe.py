
"""
======================COPYRIGHT/LICENSE START==========================

CalcHeteroNoe.py: Part of the CcpNmr Analysis program

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

import Tkinter

from math import sqrt

from memops.general import Implementation
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.LabelDivider import LabelDivider
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showWarning
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.core.AssignmentBasic import assignResToDim, propagatePeakAssignments, makeResonanceGuiName
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.DataAnalysisBasic import matchHnoePeaks
from ccpnmr.analysis.core.ExperimentBasic import findSpectrumDimsByIsotope, getSpectrumNoise
from ccpnmr.analysis.core.PeakBasic import getPeakIntensity, getPeakAnnotation, getPeakSeqCodes
from ccpnmr.analysis.frames.PeakTableFrame import PeakTableFrame

from ccp.api.nmr import Nmr
class CalcHeteroNoePopup(BasePopup):
  """
  **Calculate Heteronuclear NOE Values From Peak Intensities**
  
  The purpose of this popup window is to calculate the heteronuclear NOE for
  amide resonances based upon a comparison of the peak intensities in spectra
  that derive from an NOE saturated experiment and an unsaturated (reference)
  experiment. The basic idea of this tool is that three peak lists are chosen,
  two of which are for heteronuclear NOE experiments (H,N axes); unsaturated
  reference and saturated, and one which is the source of assignments and peak
  locations. This last "Assignment" peak list may be the same as one of the NOE
  peak lists, but may also be entirely separate.

  The "Assignment" peak list is used to specify which peak assignments and
  locations should be used for the calculation of the heteronuclear NOE values,
  and thus can be used to specify only a subset of the resonances for
  measurement. For example, it is common to copy an HQSC peak list for use as
  the "Assignment" peak list but remove overlapped and NH2 peaks so that the NOE
  values are only calculated for separated backbone amides. The calculation
  process involves taking each of these assigned peaks and finding peaks with
  the same assignment in the NOE peak lists, or if that fails finding peaks with
  a similar position (within the stated tolerances); new peaks may be picked if
  the "Pick new peaks?" option is set.

  The first "Peak Lists & Settings" tab allows the user to choose the peak lists
  and various options that will be used in the peak-finding and NOE calculation.
  The "Peaks" table allows the peaks from each of the three peak list selections
  to be displayed when one of the "Show" buttons is clicked. The [Separate Peak
  Table] function allows these peaks to be displayed in the main, separate `Peak
  Lists`_ table, which has many more peak manipulation options. The options
  below the table may be used to locate selected peaks within the spectrum
  window displays.

  The second "Peak Intensity Comparison" tab is where the heteronuclear NOE
  values are actually calculated. Assuming that two NOE experiment peak lists
  have been chosen and that some of their peaks match the assigned peak
  positions then the peak intensities are extracted and NOE values automatically
  calculated when the tab is opened. Although, a refresh of the table can be
  forced with [Find Matching Peaks] at the bottom

  If pairs of NOE saturated and reference peaks are found then the actual
  heteronuclear NOE value is displayed as the "Intensity Ratio" in the last,
  rightmost, column of the table. To store these values as a NOE measurement
  list; so that the data can be saved in the CCPN project without need for
  recalculation, the [Create Hetero NOE List] function can be used. The results
  are then available to view at any time via the `Measurement Lists`_ table.

  **Caveats & Tips**

  Erroneous peak intensity comparisons may be removed with the [Remove Pairs]
  function, but its is common to curate the "Assign" peak list first and
  avoid tidying afterwards.

  The "Closeness score" can be used to find peak positions where the compared
  NOE peaks are unexpectedly far from one another.

  .. _`Peak Lists`: EditPeakListsPopup.html
  .. _`Measurement Lists`: EditMeasurementListsPopup.html
  
  """

  def __init__(self, parent, *args, **kw):

    self.guiParent      = parent
    self.peakPairs      = []
    self.intensityType  = 'height'
    self.selectedPair   = None
    self.assignPeakList = None
    self.refPeakList    = None
    self.satPeakList    = None
    self.displayPeakList = None
    self.waiting        = 0
  
    BasePopup.__init__(self, parent, title="Data Analysis : Heteronuclear NOE", **kw)

  def body(self, guiFrame):

    self.geometry('700x700')
   
    guiFrame.expandGrid(0,0)
    
    options = ['Peak Lists & Settings','Peak Intensity Comparison']
    tabbedFrame = TabbedFrame(guiFrame, options=options, callback=self.changeTab)
    tabbedFrame.grid(row=0, column=0, sticky='nsew')
    self.tabbedFrame = tabbedFrame
    frameA, frameB = tabbedFrame.frames

    row = 0
    frameA.grid_columnconfigure(1, weight=1)
    frameA.grid_columnconfigure(3, weight=1)
    frameA.grid_columnconfigure(5, weight=1)
    frameA.grid_rowconfigure(5, weight=1)

    tipText = 'Number of reference peaks (no saturation)'
    self.peaksALabel = Label(frameA, text='Number of Ref Peaks: ', tipText=tipText)
    self.peaksALabel.grid(row=1,column=0,columnspan=2,sticky='w')

    tipText = 'Number of NOE saturation peaks'
    self.peaksBLabel = Label(frameA, text='Number of Sat Peaks: ', tipText=tipText)
    self.peaksBLabel.grid(row=1,column=2,columnspan=2,sticky='w')

    tipText = 'Number of peaks in assigned list'
    self.peaksCLabel = Label(frameA, text='Number of Assign Peaks: ', tipText=tipText)
    self.peaksCLabel.grid(row=1,column=4,columnspan=2,sticky='w')
    
    tipText = 'Selects which peak list is considered the NOE intensity reference (no saturation)'
    specALabel = Label(frameA, text='Ref Peak List: ')
    specALabel.grid(row=0,column=0,sticky='w')
    self.specAPulldown = PulldownList(frameA, callback=self.setRefPeakList, tipText=tipText)
    self.specAPulldown.grid(row=0,column=1,sticky='w')

    tipText = 'Selects which peak list is considered as NOE saturated.'
    specBLabel = Label(frameA, text='Sat Peak List: ')
    specBLabel.grid(row=0,column=2,sticky='w')
    self.specBPulldown = PulldownList(frameA, callback=self.setSatPeakList, tipText=tipText)
    self.specBPulldown.grid(row=0,column=3,sticky='w')

    tipText = 'Selects a peak list with assignments to use as a positional reference'
    specCLabel = Label(frameA, text='Assignment Peak List: ')
    specCLabel.grid(row=0,column=4,sticky='w')
    self.specCPulldown = PulldownList(frameA, callback=self.setAssignPeakList, tipText=tipText)
    self.specCPulldown.grid(row=0,column=5,sticky='w')

    frame0a = Frame(frameA)
    frame0a.grid(row=2,column=0,columnspan=6,sticky='nsew')
    frame0a.grid_columnconfigure(9, weight=1)
    
    tipText = '1H ppm tolerance for matching assigned peaks to reference & NOE saturation peaks'
    tolHLabel   = Label(frame0a, text='Tolerances: 1H')
    tolHLabel.grid(row=0,column=0,sticky='w')
    self.tolHEntry = FloatEntry(frame0a,text='0.02', width=6, tipText=tipText)
    self.tolHEntry .grid(row=0,column=1,sticky='w')  

    tipText = '15N ppm tolerance for matching assigned peaks to reference & NOE saturation peaks'
    tolNLabel   = Label(frame0a, text=' 15N')
    tolNLabel .grid(row=0,column=2,sticky='w')   
    self.tolNEntry = FloatEntry(frame0a,text='0.1', width=6, tipText=tipText)
    self.tolNEntry .grid(row=0,column=3,sticky='w')   

    tipText = 'Whether to peak new peaks in reference & NOE saturated lists (at assignment locations)'
    label = Label(frame0a, text=' Pick new peaks?', grid=(0,4)) 
    self.pickPeaksSelect = CheckButton(frame0a, tipText=tipText,
                                       grid=(0,5), selected=True)

    tipText = 'Whether to assign peaks in the peaks in the reference & NOE saturation lists, if not already assigned'
    label = Label(frame0a, text=' Assign peaks?')
    label.grid(row=0,column=6,sticky='w')   
    self.assignSelect = CheckButton(frame0a, tipText=tipText)
    self.assignSelect.set(1)
    self.assignSelect.grid(row=0,column=7,sticky='w')    

    tipText = 'Whether to consider peak height or volume in the heteronuclear NOE calculation'
    intensLabel = Label(frame0a, text=' Intensity Type:')
    intensLabel .grid(row=0,column=8,sticky='w')   
    self.intensPulldown = PulldownList(frame0a, texts=['height','volume'],
                                       callback=self.setIntensityType,
                                       tipText=tipText)
    self.intensPulldown.grid(row=0,column=9,sticky='w')    

    divider = LabelDivider(frameA, text='Peaks', grid=(3,0),
                           gridSpan=(1,6))

    tipTexts = ['Show the selected intensity reference peaks in the below table',
                'Show the selected NOE saturation peaks in the below table',
                'Show the selected assigned peak list in the below table',
                'Show the displayed peaks in a separate peak table, where assignments etc. may be adjusted']
    texts    = ['Show Ref Peaks','Show Sat Peaks',
                'Show Assign Peaks', 'Separate Peak Table']
    commands = [self.viewRefPeakList, self.viewSatPeakList,
                self.viewAssignPeakList, self.viewSeparatePeakTable]
    self.viewPeaksButtons = ButtonList(frameA, expands=True, tipTexts=tipTexts,
                                       texts=texts, commands=commands)
    self.viewPeaksButtons.grid(row=4,column=0,columnspan=6,sticky='nsew')

    self.peakTable = PeakTableFrame(frameA, self.guiParent, grid=(5,0),
                                    gridSpan=(1,6))
    self.peakTable.bottomButtons1.grid_forget()
    self.peakTable.bottomButtons2.grid_forget()
    #self.peakTable.topFrame.grid_forget()
    self.peakTable.topFrame.grid(row=2, column=0, sticky='ew')
    # Next tab

    frameB.expandGrid(0,0)
    
    tipTexts = ['Row number',
                'Assignment annotation for NOE saturation peak',
                'Assignment annotation for reference peak (no saturation)',
                '1H chemical shift of NOE saturation peak',
                '1H chemical shift of reference peak',
                '15N chemical shift of NOE saturation peak',
                '15N chemical shift of reference peak',
                'The separation between compared peaks: square root of the sum of ppm differences squared',
                'The intensity if the NOE saturation peak',
                'The intensity of the reference peak (no saturation)',
                'Ratio of peak intensities: saturated over reference',
                'Residue(s) for reference peak']
    colHeadings      = ['#','Sat Peak','Ref Peak','1H shift A',
                        '1H shift B','15N shift A','15N shift B',
                        'Closeness\nScore','Intensity A','Intensity B',
                        'Intensity\nRatio','Residue']
    self.scrolledMatrix = ScrolledMatrix(frameB, multiSelect=True, 
                                         headingList=colHeadings,
                                         callback=self.selectCell,
                                         tipTexts=tipTexts,
                                         grid=(0,0),
                                         deleteFunc=self.removePair)

    tipTexts = ['Force a manual update of the table; pair-up NOE saturation and reference peaks according to assigned peak positions',
                'Remove the selected rows of peak pairs',
                'Show peaks corresponding to the selected row in a table',
                'Save the Heteronuclear NOE values in the CCPN project as a data list']
    texts    = ['Refresh Table','Remove Pairs',
                'Show Peak Pair','Create Hetero NOE List']
    commands = [self.matchPeaks,self.removePair,
                self.showPeakPair,self.makeNoeList]
    self.pairButtons = ButtonList(frameB, tipTexts=tipTexts, grid=(1,0),
                                  texts=texts, commands=commands)


    bottomButtons = UtilityButtonList(tabbedFrame.sideFrame, helpUrl=self.help_url)
    bottomButtons.grid(row=0, column=0, sticky='e')
    
    self.updatePulldowns()
    self.updateAfter()

    self.administerNotifiers(self.registerNotify)

  def administerNotifiers(self, notifyFunc):

    for func in ('__init__', 'delete','setName'):
      for clazz in ('ccp.nmr.Nmr.DataSource', 'ccp.nmr.Nmr.Experiment',):
        notifyFunc(self.updatePulldowns, clazz, func)
        
    for func in ('__init__', 'delete'):
      notifyFunc(self.updatePulldowns,'ccp.nmr.Nmr.PeakList', func)

    for func in ('__init__', 'delete','setAnnotation','setFigOfMerit'):
      notifyFunc(self.updatePeaks, 'ccp.nmr.Nmr.Peak', func)
    for func in ('setAnnotation','setPosition','setNumAliasing'):
      notifyFunc(self.updatePeakChild, 'ccp.nmr.Nmr.PeakDim', func)
    for func in ('__init__', 'delete', 'setValue'):
      notifyFunc(self.updatePeakChild, 'ccp.nmr.Nmr.PeakIntensity', func)

  def changeTab(self, index):
  
    if index == 1:
      self.matchPeaks()

  def open(self):
  
    self.updatePulldowns()
    self.updateAfter()
    BasePopup.open(self)
    
    
  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)
 
  def updatePulldowns(self, *obj):

    index0 = 0
    index1 = 0
    index2 = 0
    names, peakLists = self.getPeakLists()
    
    if names:

      if self.refPeakList not in peakLists:
        self.refPeakList = peakLists[0]

      if  self.satPeakList not in peakLists:
        self.satPeakList = peakLists[0]

      if self.assignPeakList not in peakLists:
        self.assignPeakList = peakLists[0]

      index0 = peakLists.index(self.refPeakList)
      index1 = peakLists.index(self.satPeakList)
      index2 = peakLists.index(self.assignPeakList)
    
    self.specAPulldown.setup(names, peakLists, index0)
    self.specBPulldown.setup(names, peakLists, index1)
    self.specCPulldown.setup(names, peakLists, index2)

  def updatePeakChild(self,peakChild):
  
    if self.waiting:
      return
    
    self.updatePeaks(peakChild.peak)
  
  def updatePeaks(self, peak):
  
    if self.waiting:
      return
    
    if peak.peakList in (self.refPeakList,self.satPeakList,self.assignPeakList):
      if peak.isDeleted and (peak.peakList in (self.refPeakList,self.satPeakList) ):
        for peaks in self.peakPairs:
          if peak in peaks:
            self.peakPairs.remove(peaks)
            if self.selectedPair is peaks:
              self.selectedPair = None
            self.updateAfter()
            return
            
      self.updateAfter()
 
  def setIntensityType(self, intensityType):
    
    self.intensityType = intensityType
    self.updateAfter()
  
  def viewRefPeakList(self):
  
    if self.refPeakList:
      self.updatePeakTable(self.refPeakList)
   
  def viewSatPeakList(self):
  
    if self.satPeakList:
      self.updatePeakTable(self.satPeakList)
  
  def viewAssignPeakList(self):
  
    if self.assignPeakList:
      self.updatePeakTable(self.assignPeakList)
  
  def viewSeparatePeakTable(self):
  
    if self.displayPeakList:
      self.guiParent.editPeakList(peakList=self.displayPeakList)
  
  def setRefPeakList(self, refPeakList):
  
    if self.displayPeakList is self.refPeakList:
      self.updatePeakTable(refPeakList)

    self.refPeakList = refPeakList
    self.updateViewButtons()
    self.updateAfter()
  
  def setSatPeakList(self, satPeakList):
  
    if self.displayPeakList is self.satPeakList:
      self.updatePeakTable(satPeakList)
      
    self.satPeakList = satPeakList
    self.updateViewButtons()
    self.updateAfter()
  
  def setAssignPeakList(self, assignPeakList):
  
    if self.displayPeakList is self.assignPeakList:
      self.updatePeakTable(assignPeakList)
      
    self.assignPeakList = assignPeakList
    self.updateViewButtons()
    self.updateAfter()
  
  def getPeakListName(self, peakList):
  
    if peakList:
      spectrum = peakList.dataSource
      experiment = spectrum.experiment
      name = '%s:%s:%d' % (experiment.name, spectrum.name, peakList.serial)
    else:
      name = '<None>'
  
    return name
  
  def getPeakLists(self):
  
    names = []
    peakLists = []
    
    for experiment in self.nmrProject.sortedExperiments():
      for dataSource in experiment.sortedDataSources():
        if dataSource.numDim == 2:
          dimsN = findSpectrumDimsByIsotope(dataSource,'15N')
          dimsH = findSpectrumDimsByIsotope(dataSource,'1H')
          if len(dimsN) == 1 and len(dimsH) == 1:
            for peakList in dataSource.sortedPeakLists():
              name = self.getPeakListName(peakList)
              names.append( name )
              peakLists.append(peakList)
    
    return names, peakLists
  
  def showPeakPair(self):
  
    if self.selectedPair:
      self.guiParent.viewPeaks(self.selectedPair)
          
  def selectCell(self, object, row, col):
  
    self.selectedPair = object

    if self.selectedPair:
      self.pairButtons.buttons[1].enable()
      self.pairButtons.buttons[2].enable()
    else:
      self.pairButtons.buttons[1].disable()
      self.pairButtons.buttons[2].disable()
  
  def removePair(self, *event):
  
    pairs = self.scrolledMatrix.currentObjects
    
    if pairs:
      for pair in pairs:
        self.peakPairs.remove(pair)
        
      self.selectedPair = None
      self.updateAfter()
  
  def matchPeaks(self):
  
    # assign relative to reference
    
    if self.assignPeakList and self.assignPeakList.peaks and self.refPeakList and self.satPeakList:
 
      tolH = float( self.tolHEntry.get() )
      tolN = float( self.tolNEntry.get() )
      pickNewPeaks = self.pickPeaksSelect.get()
      doAssign     = self.assignSelect.get()
 
      dimH  = findSpectrumDimsByIsotope(self.assignPeakList.dataSource,'1H' )[0]
      dimHA = findSpectrumDimsByIsotope(self.refPeakList.dataSource,'1H' )[0]
      dimHB = findSpectrumDimsByIsotope(self.satPeakList.dataSource,'1H' )[0]
      dimN  = 1-dimH
      dimNA = 1-dimHA
      dimNB = 1-dimHB
 
      tolerancesA = [0,0]
      tolerancesA[dimHA] = tolH
      tolerancesA[dimNA] = tolN

      tolerancesB = [0,0]
      tolerancesB[dimHB] = tolH
      tolerancesB[dimNB] = tolN
 
      self.peakPairs = matchHnoePeaks(self.assignPeakList,self.refPeakList,
                                      self.satPeakList,tolerancesA,tolerancesB,
                                      pickNewPeaks,doAssign)
 
      self.updateAfter()
 
  def makeNoeList(self):

    if self.refPeakList is self.satPeakList:
      showWarning('Same Peak List',
                  'Ref Peak List and Sat Peak List cannot be the same',
                  parent=self)
      return

    if self.peakPairs:
      s1 = self.refPeakList.dataSource
      s2 = self.satPeakList.dataSource
      noiseRef = getSpectrumNoise(s1)
      noiseSat = getSpectrumNoise(s2)
      
      es = '%s-%s' % (s1.experiment.name,s2.experiment.name)
      if len(es) > 50:
        es = 'Expt(%d)-Expt(%d)' % (s1.experiment.serial,s2.experiment.serial)

      noeList = self.nmrProject.newNoeList(unit='None',name='Hetero NOE list for %s' % es)
      noeList.setExperiments([s1.experiment,])
      if s1.experiment is not s2.experiment:
        noeList.addExperiment( s2.experiment )
      # TBD: sf, noeValueType, refValue, refDescription
    
      resonancePairsSeen = set()
      for (peakA,peakB) in self.peakPairs: # peakA is sat
        intensA    = getPeakIntensity(peakA,self.intensityType)
        intensB    = getPeakIntensity(peakB,self.intensityType)
        value      = float(intensA)/intensB
        error      = abs(value) * sqrt((noiseSat/intensA)**2 + (noiseRef/intensB)**2)
        resonances = tuple(self.getPeakResonances(peakA))
        frozenResonances = frozenset(resonances)
        if len(resonances) < 2:
          pl = peakA.peakList
          sp = pl.dataSource
          msg = 'Skipping %s:%s:%d peak %d it has too few resonances assigned'
          data = (sp.experiment.name, sp.name, pl.serial, peakA.serial)
          showWarning('Warning',msg % data, parent=self)
        
        elif len(resonances) > 2:
          pl = peakA.peakList
          sp = pl.dataSource
          resonanceText = ' '.join([makeResonanceGuiName(r) for r in resonances])
          msg = 'Skipping %s:%s:%d peak %d it has too many resonances assigned (%s)'
          data = (sp.experiment.name, sp.name, pl.serial, peakA.serial, resonanceText)
          showWarning('Warning', msg % data, parent=self)
        
        elif frozenResonances not in resonancePairsSeen:
          resonancePairsSeen.add(frozenResonances)
          noeList.newNoe(value=value,resonances=resonances,peaks=[peakA,peakB],error=error)
          
        else:
          resonanceText = ' '.join([makeResonanceGuiName(r) for r in resonances])
          msg = 'Skipping duplicate entry for resonances %s' % resonanceText
          showWarning('Warning', msg, parent=self)

      self.parent.editMeasurements(measurementList=noeList)

  def getPeakResonances(self,peak):
  
    resonances = []
    for peakDim in peak.sortedPeakDims():
      for contrib in peakDim.sortedPeakDimContribs():
        resonances.append(contrib.resonance)
    
    return resonances

  def updateAfter(self, *opt):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
 
  def updateViewButtons(self):
  
    if self.refPeakList:
      self.viewPeaksButtons.buttons[0].enable()
    else:
      self.viewPeaksButtons.buttons[0].disable()
  
    if self.satPeakList:
      self.viewPeaksButtons.buttons[1].enable()
    else:
      self.viewPeaksButtons.buttons[1].disable()
  
    if self.assignPeakList:
      self.viewPeaksButtons.buttons[2].enable()
    else:
      self.viewPeaksButtons.buttons[2].disable()
  
  def updatePeakTable(self, peakList):
    
    if peakList is not self.displayPeakList:
      self.displayPeakList = peakList
      self.peakTable.update(peaks=peakList.sortedPeaks())

  
  def update(self):

    if self.refPeakList:
      self.peaksALabel.set( 'Number of Ref Peaks: %d' % len(self.refPeakList.peaks) )
    else:
      self.peaksALabel.set( 'Number of Ref Peaks: %d' % 0 )
    if self.satPeakList:
      self.peaksBLabel.set( 'Number of Sat Peaks: %d' % len(self.satPeakList.peaks) )
    else:
      self.peaksBLabel.set( 'Number of Sat Peaks: %d' % 0 )
    if self.assignPeakList:
      self.peaksCLabel.set( 'Number of Assign Peaks: %d' % len(self.assignPeakList.peaks) )
    else:
      self.peaksCLabel.set( 'Number of Assign Peaks: %d' % 0 )

    if self.refPeakList and self.satPeakList and self.assignPeakList:
      if self.refPeakList is self.satPeakList:
        self.pairButtons.buttons[0].disable()
      else:
        self.pairButtons.buttons[0].enable()
    else:
      self.pairButtons.buttons[0].disable()
 
    if self.selectedPair:
      self.pairButtons.buttons[1].enable()
      self.pairButtons.buttons[2].enable()
    else:
      self.pairButtons.buttons[1].disable()
      self.pairButtons.buttons[2].disable()

    if self.peakPairs:
      self.pairButtons.buttons[3].enable()
      dsA = self.peakPairs[0][0].peakList.dataSource
      dsB = self.peakPairs[0][1].peakList.dataSource
      dimHA = findSpectrumDimsByIsotope(dsA,'1H')[0]
      dimHB = findSpectrumDimsByIsotope(dsB,'1H')[0]
      dimNA = findSpectrumDimsByIsotope(dsA,'15N')[0]
      dimNB = findSpectrumDimsByIsotope(dsB,'15N')[0]
    else:
      self.pairButtons.buttons[3].disable()
    
    objectList  = []
    textMatrix  = []

    i = 0
    for (peakA,peakB) in self.peakPairs:
      i += 1

      peakDimsA = peakA.sortedPeakDims()
      peakDimsB = peakB.sortedPeakDims()

      ppm0 = peakDimsA[dimHA].value
      ppm1 = peakDimsB[dimHB].value
      ppm2 = peakDimsA[dimNA].value
      ppm3 = peakDimsB[dimNB].value
      d0 = abs(ppm0-ppm1)
      d1 = abs(ppm2-ppm3)
      intensA = getPeakIntensity(peakA,self.intensityType)
      intensB = getPeakIntensity(peakB,self.intensityType)
      datum = []
      datum.append( i )
      datum.append( getPeakAnnotation(peakA, doPeakDims=False) )
      datum.append( getPeakAnnotation(peakB, doPeakDims=False) )
      datum.append( ppm0 )
      datum.append( ppm1 )
      datum.append( ppm2 )
      datum.append( ppm3 )
      datum.append( sqrt((d0*d0)+(d1*d1)) )
      datum.append( intensA )
      datum.append( intensB )
      if intensB:
        datum.append( float(intensA)/intensB )
      else:
        datum.append( None )
      seqCodes = ','.join(['%s' % seqCode for seqCode in getPeakSeqCodes(peakB)])
      datum.append(seqCodes)
      
      objectList.append( (peakA,peakB) )
      textMatrix.append( datum )
      
    if not objectList:
      textMatrix.append([])

    self.scrolledMatrix.update(objectList=objectList, textMatrix=textMatrix)

    self.waiting = False

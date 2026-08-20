
"""
======================COPYRIGHT/LICENSE START==========================

InitRootAssignments.py: Part of the CcpNmr Analysis program

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

from memops.gui.ButtonList      import ButtonList, UtilityButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Frame           import Frame
from memops.gui.Label           import Label
from memops.gui.LabelDivider    import LabelDivider
from memops.gui.MessageReporter import showWarning, showError
from memops.gui.PulldownList    import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix

from ccpnmr.analysis.core.AssignmentBasic import assignResToDim, findSpinSystem, addSpinSystemResonance
from ccpnmr.analysis.popups.BasePopup     import BasePopup
from ccpnmr.analysis.core.ExperimentBasic import getOnebondDataDims, getSpectrumIsotopes
from ccpnmr.analysis.core.MarkBasic       import createPeakMark
from ccpnmr.analysis.core.WindowBasic     import isSpectrumInWindowPane, getActiveWindows
from ccpnmr.analysis.core.WindowBasic     import getWindowPaneName, zoomToShowPeaks

ALLOWED_REF_EXPTS = ['H[N]','H[N[CO]]','H[N[co[CA]]]','H[C]']

def testInitialiseRootPeakListPopup(argServer):

  popup = InitRootAssignmentsPopup(argServer.parent)
  popup.open()
  
class InitRootAssignmentsPopup(BasePopup):
  """
  **Add Starting Resonances and Spin Systems to Root Spectra**
  
  The purpose of this popup window is to setup certain kinds of spectra at the very
  start of the assignment process. This initialisation involves adding new, anonymous
  resonances and spin system (residue) assignments to picked, but unassigned peaks. The
  kinds of spectra that have their peak lists processed in this way are normally those
  that, on the whole, have one peak for each residue in the molecular system. This
  typically involves 15N HSQC and HNCO spectra where you get one peak for each NH
  group. In these instances two resonance are added for each peak; one for the 15N
  dimension and one for the 1H dimension, both resonances for the peak are then added
  to the same spin system (a resonance grouping that stands-in for a residue).

  The initial resonances and spin system groups that have been added to an initialised
  "root" spectrum then serve as the basis for assigning related spectra, often with
  higher dimensionality. The general principle is that the positions of the peaks, and
  hence the chemical shifts of their assigned resonances, serve as guides to pick,
  locate and assign peaks in other spectra that share the same resonances in some of
  their dimensions. For example peaks in a 15N HSQC spectrum can be used to link
  equivalent peaks in a 3D HNCA spectrum or 3D 15N HSQC-NOSEY spectrum of the same
  sample because peaks in the 3D spectra derive from the same amide groups. Once
  related spectra have been linked via assignments to a common set of resonances and
  sequential assignment has been performed the resonances and spin systems will no
  longer be anonymous; the resonances will be assigned to specific atoms, and the spin
  systems will be assigned to residues.

  One complication of the initialisation process is that not all peaks in the root 
  spectra relate to unique residues or spin systems. The most common such examples are
  peaks that derive from NH2 (amide) side chain groups of Glutamine and Asparagine
  residues. In an 15N HSQC spectrum an NH2 group will give rise to two peaks one for
  each hydrogen, but these peaks share the same nitrogen and thus appear at the same
  15N frequency. When a root spectrum is initialised such peaks must be taken care of
  in the appropriate way; the two peaks of an NH2 group share the same nitrogen
  resonance *and* the same spin system number.

  The "Amide Side Chain Peaks" table in the popup window lists all of the pairs of
  peaks that have similar 15N resonance positions, and which lie in the chemical shift
  regions for NH2 groups. The purpose of the table is for the user to confirm or deny
  that the pairs of peak really are NH2 groups, and not coincidental matches. The user
  can locate a potential NH2 peak pair by selecting an appropriate spectrum window and,
  assuming "Follow Amides" is set, clicking on the row for the peak pair. If the peaks
  look like a genuine NH2 pair (similar intensity, deuterium sub-peak etc) then the
  peaks are may be confirmed by double clicking in the "Confirmed?" column.
  
  With any NH2 groups confirmed, the peak list is initialised via [Initialise Peak
  List!] at the top; all peaks will have new resonances and spin systems added, and
  peaks from NH2 groups will be linked appropriately.

  **Caveats & Tips**
  
  This initialisation tool is currently limited to the following types of experiment:
  15N HSQC, HNCO, HNcoCA, 13C HSQC.
  
  Only potential NH2 peak pairs with 15N positions that are within the stated tolerance
  are shown. This tolerance can be reduced to make the amide search more specific,
  although this may be at the expense of detecting NH2 pairs that are distorted due to
  peak overlap.
  """

  def __init__(self, parent, *args, **kw):
    
    self.waiting    = False
    self.peakList   = None
    self.amidePair  = None
    self.amidePairs = []
    self.isAmide    = {}
    self.guiParent  = parent
    self.windowPane = None
    self.marks = []
    
    BasePopup.__init__(self, parent=parent, title='Assignment : Initialise Root Resonances', **kw)

  def open(self):
  
    self.update()
    BasePopup.open(self)

  def close(self):

    self.clearMarks()
    BasePopup.close(self)

  def body(self, guiFrame):

    self.geometry('550x600')
    guiFrame.grid_columnconfigure(0, weight=1)
    
    row = 0
    
    frame = Frame(guiFrame, grid=(row,0))
    frame.expandGrid(0,4)

    label = Label(frame, text='Peak List:', grid=(0,0))
    tipText = 'For spectra with (mostly) one peak per residue, selects which peak list to initialise; by adding anonymous resonance and spin system numbers'
    self.peakListPulldown = PulldownList(frame, self.changePeakList,
                                         grid=(0,1), tipText=tipText)

    label = Label(frame, text='15N tolerance (ppm):', grid=(0,2))
    tipText = 'The upper limit in the difference between 15N ppm positions for two peaks to be considered as a potential amide'
    self.tolEntry = FloatEntry(frame, text=0.02, width=8, grid=(0,3),
                               sticky='ew', tipText=tipText)

    frame2 = Frame(frame, grid=(1,0), gridSpan=(1,4), sticky='ew')
    frame2.grid_columnconfigure(5, weight=1)

    label = Label(frame2, text='Follow Amides:', grid=(0,0))
    tipText = 'Sets whether to follow the H-N locations of potential amide peaks when clicking in rows of the table'
    self.followSelect = CheckButton(frame2, callback=None, grid=(0,1),
                                    selected=True, tipText=tipText)

    self.windowLabel = Label(frame2, text='Window:', grid=(0,2))

    tipText = 'Selects the spectrum window in which to show positions of potential amide peaks'
    self.windowPulldown = PulldownList(frame2, self.selectWindowPane,
                                       grid=(0,3), tipText=tipText)

    label = Label(frame2, text='Mark Amides:', grid=(0,4))
    tipText = 'Whether to put a multi-dimensional cross mark though the H-N positions of the selected peak pair'
    self.markSelect = CheckButton(frame2, callback=None, grid=(0,5),
                                  selected=True, tipText=tipText)

                                         
    utilButtons = UtilityButtonList(guiFrame, grid=(row,1), sticky='ne',
                                     helpUrl=self.help_url)

    row += 1
    tipTexts = ['For the stated peak list, considering confirmed amide side chain peaks, add spin system and resonance numbers to all peaks',]
    commands = [self.initialisePeakList,]
    texts    = ['Initialise Peak List!',]
    self.initButtons = ButtonList(guiFrame, commands=commands,
                                  grid=(row, 0), texts=texts,
                                  gridSpan=(1,2), tipTexts=tipTexts)
    self.initButtons.buttons[0].config(bg='#B0FFB0')
    
    row += 1
    div = LabelDivider(guiFrame, text='Amide Side Chain Peaks',
                       gridSpan=(1,2), grid=(row,0))

    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)

    frame = Frame(guiFrame, gridSpan=(1,2), grid=(row,0))
    frame.expandGrid(0,0)
    
    editSetCallbacks = [None] * 8
    editGetCallbacks = [None] * 8
    editWidgets = [None] * 8
    editGetCallbacks[1] = self.toggleConfirm
    #self.userCodeEntry = Entry(self,text='', returnCallback=self.setResidueCode, width=6)
    tipTexts = ['Number of the peak pair for the table',
                'Whether the peak pair is confirmed as being from an amide side chain; a common nitrogen by different hydrogens',
                'The difference in 15N ppm position between the two peaks of the pair',
                'The assignment annotation of the first peak in the pair',
                'The assignment annotation of the second peak in the pair',
                'The average 15N position of the peaks in the pair',
                'The 1H chemical shift of the first peak in the pair',
                'The 1H chemical shift of the second peak in the pair']
    headingList      = ['#','Confirmed?',u'15N\n\u0394ppm',
                        'Peak 1','Peak 2',
                        'Shift N','Shift H1','Shift H2']
    self.amidePairMatrix = ScrolledMatrix(frame, headingList=headingList,
                                          editSetCallbacks=editSetCallbacks,
                                          editGetCallbacks=editGetCallbacks,
                                          editWidgets=editWidgets,
                                          callback=self.selectAmidePair,
                                          multiSelect=True,
                                          grid=(0,0), tipTexts=tipTexts)
                                        
    tipTexts = ['Confirm the selected peak pairs as being amide side chain peaks; a common nitrogen by different hydrogens',
                'remove any conformation for the selected peak pairs being amide side chain peaks',
                'Manually force an update the table of potential pairs of amide side chain peaks']
    commands = [self.confirmAmidePairs,
                self.unconfirmAmidePairs,
                self.predictAmidePairs]
    texts    = ['Confirm\nSelected',
                'Unconfirm\nSelected',
                'Update\nTable']
    self.amideButtons = ButtonList(frame, texts=texts, commands=commands, grid=(1,0), tipTexts=tipTexts)
 
    self.updatePeakLists()
    self.administerNotifiers(self.registerNotify)
  
  def toggleConfirm(self, peakPair):
  
    peakA, peakB = peakPair
  
    boolean = not peakA.amideConfirmed
    
    peakA.amideConfirmed = boolean
    peakB.amideConfirmed = boolean
    
    self.updateAfter()  
  
  def confirmAmidePairs(self):  
  
    for peakA, peakB in self.amidePairMatrix.currentObjects:
      peakA.amideConfirmed = True
      peakB.amideConfirmed = True
    
    self.updateAfter()  
  
  def unconfirmAmidePairs(self):  
  
    for peakA, peakB in self.amidePairMatrix.currentObjects:
      peakA.amideConfirmed = False
      peakB.amideConfirmed = False

    self.updateAfter()  
    
  def updatePeakListsAfter(self, obj):
  
    self.after_idle(self.updatePeakLists)

  def updateWindowListsAfter(self, obj):

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

  def contribUpdateAfter(self, contrib):

    if self.waiting:
      return
  
    if contrib.peakDim.peak.peakList is self.peakList:
      self.updateAfter()

  def markPeaks(self):

    if self.amidePair:
      peakA, peakB = self.amidePair
      markA = createPeakMark(peakA, lineWidth=2.0, remove=False)
      markB = createPeakMark(peakB, lineWidth=2.0, remove=False)
      self.clearMarks()
      self.marks = [markA, markB]
 
  def clearMarks(self):
  
    for mark in self.marks:
      if not mark.isDeleted:
        mark.delete()
    
    self.marks = []    
 
  def followPeaks(self):

    if self.amidePair and self.windowPane:
      self.windowPane.getWindowFrame()
      zoomToShowPeaks(self.amidePair, self.windowPane)

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
  
    index = 0
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
    
     
    if peakList is not self.peakList:
      self.changePeakList(peakList)
    
    self.peakListPulldown.setup(names, peakLists, index)
    
  def changePeakList(self, peakList):
    
    if peakList is not self.peakList:
      self.peakList = peakList
      self.predictAmidePairs()
      self.updateAfter()
      self.updateButtons()
      self.updateWindowList()
    

  def selectAmidePair(self, obj, row, col):

    self.amidePair = obj
    
    if self.followSelect.get():
      self.followPeaks()
    
    if self.markSelect.get():
      self.markPeaks()
      
    self.updateButtons()

  def updateAfter(self):

    if not self.waiting:
      self.waiting = True
      self.after_idle(self.update)

  def update(self):

    textMatrix  = []
    objectList  = []
    colorMatrix = []
    #['#','Delta\n15N ppm','Peak 1','Peak 2','Shift N','Shift H1','Shift H2']
    
    i = 0
    for amidePair in self.amidePairs:
      peakA, peakB, delta, ppmN, ppmH1, ppmH2 = amidePair
      i += 1
      
      colors = [None] * 8
      
      if peakA.amideConfirmed and peakB.amideConfirmed:
        colors[1] = '#B0FFB0'
        confText  = 'Yes'
      else:
        colors[1] = '#FFB0B0'
        confText  = 'No'
      
      datum = [i,confText,delta,
               '\n'.join([pd.annotation or '-' for pd in peakA.sortedPeakDims()]),
               '\n'.join([pd.annotation or '-' for pd in peakB.sortedPeakDims()]),
               ppmN, ppmH1, ppmH2]
               
      textMatrix.append(datum)
      objectList.append((peakA,peakB))
      colorMatrix.append(colors)
    
    self.amidePairMatrix.update(textMatrix=textMatrix,
                                objectList=objectList,
                                colorMatrix=colorMatrix)
    self.waiting = False

  def updateButtons(self):
  
    buttons  = self.amideButtons.buttons
    buttons2 = self.initButtons.buttons
  
    if self.peakList:
      buttons2[0].enable()
      buttons[2].enable()
     
      if self.amidePair:
        buttons[0].enable()
        buttons[1].enable()
      
      else:
        buttons[0].disable()
        buttons[1].disable()
     
    else:
      buttons[0].disable()
      buttons[1].disable()
      buttons[2].disable()
       
  def predictAmidePairs(self):

    self.amidePair = None
    isAmide    = self.isAmide = {}
    amidePairs = self.amidePairs = []
    peakList   = self.peakList

    if peakList:
      tol      = self.tolEntry.get() or 0.001
      spectrum = self.peakList.dataSource
      dimPairs = getOnebondDataDims(spectrum)

      if not dimPairs:
        return
      
      isotopes = getSpectrumIsotopes(spectrum)
      
      hDim = None
      for dimA, dimB in dimPairs:
        i = dimA.dim-1
        j = dimB.dim-1
        
        if ('1H' in isotopes[i]) and ('1H' not in isotopes[j]):
          hDim = i
          nDim = j
          
        elif ('1H' in isotopes[j]) and ('1H' not in isotopes[i]):
          hDim = j
          nDim = i
      
        if hDim is not None:
          break
      
      isAssigned = {}
      matchPeaks = []
      for peak in peakList.peaks:
        if isAmide.get(peak):
          continue
      
        peakDims = peak.sortedPeakDims()
        ppmN = peakDims[nDim].realValue
        ppmH = peakDims[hDim].realValue
        
        for contrib in peakDims[nDim].peakDimContribs:
          resonance = contrib.resonance
          
          for contrib2 in resonance.peakDimContribs:
            peak2 = contrib2.peakDim.peak
            if (peak2.peakList is peakList) and (peak2 is not peak):
              isAmide[peak]  = peak2
              isAmide[peak2] = peak
              
              if hasattr(peak, 'amideConfirmed'):
                 peak2.amideConfirmed = peak.amideConfirmed
              else:
                 peak.amideConfirmed = True
                 peak2.amideConfirmed = True
                 
              ppmH2  = peak2.findFirstPeakDim(dim=hDim+1).value
              ppmN2  = peak2.findFirstPeakDim(dim=nDim+1).value
              ppmNm  = 0.5*(ppmN+ppmN2)
              deltaN = abs(ppmN-ppmN2)
              amidePairs.append((peak,peak2,deltaN,ppmNm,ppmH, ppmH2))
              break
              
          else:
            continue
          break  
              
        else:
          if ppmN > 122.0:
            continue
          elif ppmN < 103:
            continue
 
          if ppmH > 9.0:
            continue
          elif ppmH < 5.4:
            continue
 
          if not hasattr(peak, 'amideConfirmed'):
            peak.amideConfirmed = False
          matchPeaks.append((ppmN, ppmH, peak))
 
      N = len(matchPeaks)
      unassignedPairs = []
      for i in range(N-1):
        ppmN, ppmH, peak = matchPeaks[i]
 
        for j in range(i+1,N):
          ppmN2, ppmH2, peak2 = matchPeaks[j]
 
          deltaN = abs(ppmN2-ppmN)
          if deltaN > tol:
            continue
 
          if abs(ppmH-ppmH2) > 1.50:
            continue
 
          ppmNm = 0.5*(ppmN+ppmN2)
 
          isAmide[peak]  = peak2
          isAmide[peak2] = peak
          unassignedPairs.append((deltaN,peak,peak2,ppmNm,ppmH, ppmH2))
      
      unassignedPairs.sort()
      done = {}
      for deltaN,peak,peak2,ppmNm,ppmH,ppmH2 in unassignedPairs:
        if done.get(peak):
          continue
        if done.get(peak2):
          continue
          
        done[peak]  = True
        done[peak2] = True
        
        amidePairs.append((peak,peak2,deltaN,ppmNm,ppmH, ppmH2))
      
      peaks = isAmide.keys()
      for peak in peaks:
        if done.get(peak) is None:
          del isAmide[peak]
       
      self.updateAfter()

  def initialisePeakList(self):

      
    isAmide  = self.isAmide
    peakList = self.peakList

    if not peakList:
      showWarning('Warning','No peak list available or selected', parent=self)
      return
 
    experimentType = peakList.dataSource.experiment.refExperiment.name

    for peak in peakList.peaks:

      peakDims = peak.sortedPeakDims()

      spinSystem = None
      resonances = []
 
      for peakDim in peakDims:
        dataDimRef = peakDim.dataDimRef
 
        if not dataDimRef:
          continue
 
        isotopes = dataDimRef.expDimRef.isotopeCodes
 
 
        if not peakDim.peakDimContribs:
          if '15N' in isotopes:
            if hasattr(peak, 'amideConfirmed') and  peak.amideConfirmed and isAmide.get(peak):
              peakDim2 = isAmide[peak].findFirstPeakDim(dim=peakDim.dim)
              contrib2 = peakDim2.findFirstPeakDimContrib()
 
              if contrib2:
                contrib = assignResToDim(peakDim, contrib2.resonance)
              else:
                contrib = assignResToDim(peakDim)
                assignResToDim(peakDim2, contrib.resonance)
 
            else:
              contrib = assignResToDim(peakDim)
          else:
            contrib = assignResToDim(peakDim)
        else:
          contrib = peakDim.findFirstPeakDimContrib()
        
        if not contrib:
          continue
        
        resonance = contrib.resonance

        if '13C' in isotopes:
          if experimentType == 'H[N[CO]]':
            resonance.setAssignNames(['C',])
          elif experimentType == 'H[N[co[CA]]]':
            resonance.setAssignNames(['CA',])
 
          continue
 
        resonances.append( resonance )
 
        if isAmide.get(peak) and hasattr(peak, 'amideConfirmed') and  peak.amideConfirmed:
          continue
 
        if ('1H' in isotopes) and (experimentType != 'H[C]'):
          resonance.setAssignNames(['H',])
        elif '15N' in isotopes:
          resonance.setAssignNames(['N',])
 
      for resonance in resonances:
        if resonance.resonanceGroup:
          spinSystem = resonance.resonanceGroup
 
      for resonance in resonances:
        if spinSystem is None:
          spinSystem = findSpinSystem(resonance)
        addSpinSystemResonance(spinSystem, resonance)
 
  def administerNotifiers(self, notifyFunc):
    
    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.PeakDimContrib',):
        notifyFunc(self.contribUpdateAfter, clazz, func)

    for func in ('__init__', 'delete','setName'):
      for clazz in ('ccp.nmr.Nmr.PeakList', 'ccp.nmr.Nmr.DataSource', 'ccp.nmr.Nmr.Experiment'):
        notifyFunc(self.updatePeakListsAfter, clazz, func)

    notifyFunc(self.updatePeakListsAfter, 'ccp.nmr.Nmr.Experiment','setRefExperiment')

    for func in ('__init__', 'delete','setName'):
      notifyFunc(self.updateWindowListsAfter, 'ccpnmr.Analysis.SpectrumWindow', func)

    for func in ('__init__', 'delete','setName'):
      notifyFunc(self.updateWindowListsAfter, 'ccpnmr.Analysis.SpectrumWindowPane', func)

    for func in ('__init__', 'delete','setAnnotation'):
      notifyFunc(self.peakUpdateAfter, 'ccp.nmr.Nmr.Peak', func)
      
    for func in ('setAnnotation','setPosition','setNumAliasing'):
      notifyFunc(self.peakDimUpdateAfter, 'ccp.nmr.Nmr.PeakDim', func)
      
  def destroy(self): 
  
    self.administerNotifiers(self.unregisterNotify)
    self.clearMarks()
    BasePopup.destroy(self)

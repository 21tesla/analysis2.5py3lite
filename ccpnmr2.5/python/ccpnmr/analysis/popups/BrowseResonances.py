
"""
======================COPYRIGHT/LICENSE START==========================

BrowseResonances.py: Part of the CcpNmr Analysis program

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
import re

from memops.general import Implementation


from memops.gui.Button          import Button
from memops.gui.ButtonList      import ButtonList, UtilityButtonList
from memops.gui.Entry           import Entry
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Frame           import Frame
from memops.gui.Label           import Label
from memops.gui.LabelDivider    import LabelDivider
from memops.gui.LabelFrame      import LabelFrame
from memops.gui.PulldownMenu    import PulldownMenu
from memops.gui.PulldownList    import PulldownList
from memops.gui.MessageReporter import showInfo, showOkCancel, showWarning, showYesNo
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.Text            import Text


from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.AssignmentBasic import assignResonanceType, getResonanceResidue, findResonanceSet, \
                                            deassignResonance, getResonanceAtomTuple, removeSpinSystemResonance, \
                                            makeResonanceGuiName, assignAtomsToRes, mergeResonances, \
                                            addSpinSystemResonance, removeSpinSystemResonance, assignResonanceResidue, \
                                            swapProchiralResonance, splitResonance, clearResonancePeakDimContribs, \
                                            getShiftLists, getResidueResonances, getOnebondResonances, \
                                            updateResonanceAnnotation

from ccpnmr.analysis.core.WindowBasic import getWindowPaneName, displayStrips, gotoResonancePosition, \
                                             windowPaneHasValueAxis
from ccpnmr.analysis.core.MarkBasic import createRuler
from ccpnmr.analysis.core.MoleculeBasic import getResidueCode
from ccpnmr.analysis.core.CouplingBasic import getPeakDimComponentSplitting

ISOTOPE_NAMES  = ['Any','1H','13C','15N','31P','2H','19F','unknown']
NORM_CONFIG = { 'fg': 'black', 'relief': 'raised' , 'bg': 'grey83'}
CANCEL_CONFIG = { 'fg': 'white', 'relief': 'raised' , 'bg': '#%02x%02x%02x' % (128, 128, 200)}
ATOM_COLOUR = {'H':'#a0d0a0','N':'#a0a0d0',
               'C':'#d0d0a0','O':'#d0a0a0',
               'P':'#d0a0d0','S':'#d0c090',
               'Se':'#e0b0a0','F':'#d0a0d0'}

# TBD: Allow shiftList pulldown with 
#      resonance selection mode?

class BrowseResonancesPopup(BasePopup):
  """
  **A Table of Resonances Organised by Chemical Shift List**
  
  This popup window is basically a table of all the resonances that have
  chemical shifts contained within a selected shift list. This data, in-keeping
  with the CCPN Resonance entities it represents, serves as the central record
  that ties atom assignments to peak assignments. An individual resonance many
  be anonymous and not (yet) carry an atomic assignment but still exists as an
  entity that carries a chemical shift record and connections to peaks in
  spectra. Each resonance row in the table may be considered as a "phenomenon
  number" that may be used to connect various things together when assignments
  are made. Within a given shift list each resonance will one have only one
  chemical shift value, but a resonance may have different chemical shift values
  in separate shift lists, for example corresponding to different conditions;
  where the resonance phenomenon is the same (from the same atom) but its
  position in spectra changes.

  The resonance table serves both as a source of information in its own right
  and a means of finding information related to a given set of assignments. For
  example the user may wish to sort the resonances according to residue sequence
  (click on the "Residue" heading), scroll to and select the desired residue and
  then [Show All Peaks] that are associated with the residue assignment. In
  normal operation the user selects the appropriate shift list (e.g. for the set
  of conditions) and any isotope, status or molecular subset to fill the table
  with a series of resonance records. Each resonance is shown with its chemical
  shift value, any atomic or residue assignment, isotope type and various other
  useful information. From the table, operations may be performed on a selected
  subset of resonances by selecting the relevant row; left click +/-
  <Ctrl>/<Shift>, and then clicking on the buttons.

  **Navigation**
  
  The upper "Navigation & Marks" section allows the user to display resonance
  positions in the spectrum windows. The user can build strips (sub-divisions)
  with the selected resonance positions (e.g. selecting H & N amide resonances),
  center the contour display and create marker lines at the locations. Although
  similar functionality is also available from the peak table, the user cam get
  more flexible control by selecting resonances, for example to navigate to a
  location in a spectrum where a peak is *expected*, but has not been picked. It
  should be noted that the [Strip Selected] and [Go To position] functions
  operate on the spectrum window selected in the "Window" pulldown menu. 

  **Atom Assignment**
   
  The penultimate row of buttons represent functions that change the atomic
  assignment if the selected resonance (or resonances, depending on the
  operation). Accordingly, a resonance may be disconnected from atoms, linked to
  atoms, given an atom type (e.g. a general "HA") and added to or removed from a
  spin system; a group of resonances belonging to a residue. The last two options
  relate only to resonances that are deemed to be prochiral, for example like
  Ser HB2/HB3, which often cannot be distinguished until the latter stages of 3D
  structure calculation.

  **General & Information**

  The lower row of buttons relate to functions that do not involve atom
  assignment. For example resonances may be deleted; either those that are
  selected or those that are considered to be "orphans" - without any current
  assignments to spectrum peaks. The two "Show Peaks" buttons display tables  of
  the spectrum peaks that are currently linked to the selected resonances,
  optionally considering only spectra that use the selected shift list. The last
  [Info] button opens the `Resonance Info`_ table for the last selected
  resonance to display the positions of all the peak dimensions to which it is
  assigned and all its chemical shift values, from multiple shift lists.

  **Caveats & Tips**

  The [Delete Orphans] function is commonly used to tidy up projects where there
  are large numbers of resonances that are no longer assigned to any peaks. This
  may occur for example resonances have been automatically added to newly picked
  peaks but once assignment is complete many of the resonances turned out to  be
  duplicate, e.g. separate resonances were initially added to HNCA and HNcoCA
  peaks even though the resonances should be shared.

  If a resonance has a large standard deviation in its chemical shift use the
  [Info] button to show the peak dimension positions that formed the chemical
  shift average.

  If the list of resonances in the table is too large to easily find a
  particular resonance, consider using the filtering selections at the top or
  click the "?" in the table headings to filter directly on table content.

  .. _`Resonance Info`: ResonanceInfoPopup.html
  
  """


  def __init__(self, parent, resonances=None, shiftList=None, *args, **kw):

    self.guiParent = parent
    self.isotope = 'Any'
    self.status = 'Any'
    self.shiftList = shiftList
    self.chain = None
    self.ccpCode = None
    self.rulers = []
    self.residueRulers = []
    self.windowPane = None
    self.resonance = None
    self.shift = None
    self.resonances = resonances
    
    self.waitingForAtoms = False
    self.wainingForType  = False
    self.waitingForSpinSystem = False
    self.waitingForTentative = False
    
    if resonances:
      title = "Resonance : Selected Resonances"
    else:
      title = "Resonance : Resonances"
    BasePopup.__init__(self, parent=parent, title=title, **kw)
   
  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)
    
  def body(self, guiFrame):
    
    self.geometry("+150+150")

    statusNames    = ['Any','Assigned','Unassigned','Orphaned']

    self.shiftValueEntry = FloatEntry(self, returnCallback=self.setShiftValue, width=10)
    self.nameEntry       = Entry(self, returnCallback=self.setName, width=12)

    if not self.shiftList:
      shiftLists = getShiftLists(self.nmrProject)
      if shiftLists:
        self.shiftList = shiftLists[0] 

    guiFrame.grid_columnconfigure(0, weight=1)
    
    row = 0
    
    if not self.resonances:
      frame = Frame(guiFrame, grid=(row,0))
      frame.expandGrid(0,0)
      div = LabelDivider(frame, text='Selection', grid=(0,0))
      utilButtons = UtilityButtonList(frame, helpUrl=self.help_url,
                                      grid=(0,1))

      row += 1
      frame = Frame(guiFrame, grid=(row,0))
      frame.expandGrid(0,9)

      tipText = 'Selects which chemical shift list to display resonances for, and hence which\n' + \
                'shift value to use for a resonance. Includes an option to show resonances without any shifts'
      self.shiftListLabel    = Label(frame, text ='Shift List:', grid=(0,0))
      self.shiftListPulldown = PulldownList(frame, self.changeShiftList, grid=(0,1), tipText=tipText)

      tipText = 'Selects subsets of resonances to list, according to their isotope type'
      self.isotopeLabel      = Label(frame, text =' Isotope:', grid=(0,2))
      self.isotopePulldown   = PulldownList(frame, self.changeIsotope,
                                            ISOTOPE_NAMES, grid=(0,3), tipText=tipText)

      tipText = 'Selects subsets of resonances to list, according to their assignment status'
      self.statusLabel       = Label(frame, text =' Status:', grid=(0,4))
      self.statusPulldown    = PulldownList(frame, self.changeStatus,
                                            statusNames, grid=(0,5), tipText=tipText)

      tipText = 'Allows restriction of the resonance list to only those assigned to a given molecular chain'
      self.chainLabel        = Label(frame, text =' Chain:', grid=(0,6))
      self.chainPulldown     = PulldownList(frame, self.changeChain, grid=(0,7), tipText=tipText)

      tipText = 'Allows restriction of the resonance list to only those assigned toa particular kind of residue'
      self.ccpCodeLabel      = Label(frame, text =' Residue Code:', grid=(0,8))
      self.ccpCodePulldown   = PulldownList(frame, self.changeCcpCode, grid=(0,9), tipText=tipText)
 
      row += 1
    
    frame = Frame(guiFrame, grid=(row,0))
    frame.expandGrid(0,0)
    tipText = 'Options to locate resonance positions in spectrum windows'
    div = LabelDivider(frame, text='Navigation & Marks',
                       tipText=tipText, grid=(0,0))
    
    if self.resonances:
      utilButtons = UtilityButtonList(frame, helpUrl=self.help_url,
                                      grid=(0,1))

    row += 1
    frame = Frame(guiFrame, grid=(row,0))
    frame.expandGrid(0,5)

    tipText = 'Make strips in the specified window using the positions of the selected resonances'    
    self.stripButton = Button(frame, text='Strip Selected', tipText=tipText,
                              command=self.stripSelected, bd=1)
    self.stripButton.grid(row=0,column=0,sticky='ew')

    tipText = 'Navigate to the selected resonance positions in the specified window'    
    self.gotoButton = Button(frame, text='Go To Position', tipText=tipText,
                             command=self.gotoPosition, bd=1)
    self.gotoButton.grid(row=0,column=1,columnspan=2,sticky='ew')
    
    tipText = 'Selects the spectrum window to use for resonance based navigation'    
    label = Label(frame, text='Window:', grid=(0,3))
    self.windowPulldown = PulldownList(frame, callback=self.changeWindow,
                                       grid=(0,4), tipText=tipText)
    
    tipText = 'Put multi-dimensional cross marks at the selected resonance positions'    
    self.markButton = Button(frame, text='Mark Selected', grid=(0,6),
                              command=self.markSelected, bd=1, sticky='ew', tipText=tipText)

    tipText = 'Remove the resonance-located cross marks'    
    self.clearMarksButton = Button(frame, text='Clear Marks', grid=(0,7,),
                                   command=self.clearMarks, bd=1, sticky='ew', tipText=tipText)

    row += 1
    frame = Frame(guiFrame, grid=(row,0))
    frame.expandGrid(0,0)
    tipText = 'Options to create strips for residue resonances in three spectrum windows'
    div = LabelDivider(frame, text='Strip Residue Resonances',
                       tipText=tipText, grid=(0,0))
    
    row += 1
    frame = Frame(guiFrame, grid=(row,0))
    frame.expandGrid(0,5)
    
    tipText = 'Selects the first spectrum window to use for residue resonances'    
    label = Label(frame, text='>2D Window:', grid=(0,0))
    self.windowRes1Pulldown = PulldownList(frame, grid=(0,1), tipText=tipText)
    
    tipText = 'Selects the second spectrum window to use for residue resonances, which should have y and z axes swapped relative to first window'    
    label = Label(frame, text='>2D Window:', grid=(0,2))
    self.windowRes2Pulldown = PulldownList(frame, grid=(0,3), tipText=tipText)
    
    tipText = 'Selects a 2D spectrum window to use for residue resonances'    
    label = Label(frame, text='2D Window:', grid=(0,4))
    self.windowRes3Pulldown = PulldownList(frame, grid=(0,5), tipText=tipText)
    
    tipText = 'Create strips at the residue resonance positions'    
    self.stripResButton = Button(frame, text='Make Residue Strips', grid=(1,0), gridSpan=(1,2),
                                command=self.stripResidueResonances, bd=1, sticky='ew', tipText=tipText)
    
    tipText = 'Put multi-dimensional cross marks at the residue resonance positions'    
    self.markResButton = Button(frame, text='Mark Residue Resonances', grid=(1,2), gridSpan=(1,2),
                                command=self.markResidueResonances, bd=1, sticky='ew', tipText=tipText)
    
    tipText = 'Remove the marks at the residue resonance positions'    
    self.clearMarksButton = Button(frame, text='Clear Residue Marks', grid=(1,4,), gridSpan=(1,2),
                                   command=self.clearResidueMarks, bd=1, sticky='w', tipText=tipText)
    
    row += 1
    tipText = 'The main table displaying the resonance records in the selected shift list, and considering filtering options'
    div = LabelDivider(guiFrame, text='Resonance Table',
                       tipText=tipText, grid=(row,0))

    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    
    justifyList = ['center','right','left','left','center',
                   'center','center','center','center','center']

    tipTexts = ['The serial number of the resonance in the NMR project',
                'The chemical shift value of the resonance in the selected shift list',
                'The standard deviation in the chemical shift value for the resonance',
                'The name bestowed upon the resonance by virtue of its atomic assignments',
                'The residue to which the resonance may be assigned',
                'The type of isotope responsible for the resonance',
                'A user-editable name for the resonance, often displayed if atom assignments are absent',
                'The total number of peaks assigned (via their dimensions) to a resonance',
                'The number of peaks assigned to a resonance, considering only the selected shift list',
                'The serial number of the spin system in which the resonance resides',
                'The resonance details']
                   
    colHeadings = ['#',"Shift","SD",'Assign\nName','Residue',
                   'Isotope','Other\nName','All\nPeaks',
                   'Shift List\nPeaks','Spin\nSystem #','Details']
                   
    self.detailsEntry = Entry(self, width=10, returnCallback=self.setDetails)
    
    editWidgets      = [None,self.shiftValueEntry,None,None,
                        None,None,self.nameEntry,None,None,None,self.detailsEntry]
    editGetCallbacks = [None,self.getShiftValue,  None,None,
                        None,None,self.getName  ,self.showPeaks,
                        self.showShiftListPeaks,None,self.getDetails]
    editSetCallbacks = [None,self.setShiftValue,  None,None,
                        None,None,self.setName  ,None,None,None,self.setDetails]

    self.scrolledMatrix = ScrolledMatrix(guiFrame, justifyList=justifyList,
                                         initialRows=16, headingList=colHeadings,
                                         editWidgets=editWidgets,
                                         editGetCallbacks=editGetCallbacks,
                                         editSetCallbacks=editSetCallbacks,
                                         callback=self.selectCell,
                                         multiSelect=True, grid=(row,0),
                                         gridSpan=(1,10), tipTexts=tipTexts)
    self.scrolledMatrix.doEditMarkExtraRules = self.doEditMarkExtraRules
    
    row += 1
    tipText = 'Functions to change how the selected resonances are assigned to atoms'
    div = LabelDivider(guiFrame, text='Atom Assignment',
                       grid=(row,0), tipText=tipText)

    row += 1
    tipTexts = ['Remove any atomic assignments from the selected resonances',
                'Assign or re-assign the last selected resonance to atoms (opens the atom browser)',
                'Set the assignment type of a resonance by selecting a representative atom',
                'Add the last selected resonance to a spin system grouping, if not already in one',
                'Remove the selected resonances form any spin system grouping',
                'If the last selected resonance is a prochiral, swap its atom assignment for the other stereochemical option',
                'For the last selected resonance, make any unambiguous prochiral assignment ambiguous, i.e HB2 -> HBa/HBb']
    texts    = ['Deassign','Assign',
                'Set Atom\nType',
                'Add to\nSpin System',
                'Remove from\nSpin System',
                'Swap\nProchirals',
                'Ambiguate\nProchirals']
    commands = [self.disconnectAtoms,
                self.waitForAtoms,
                self.waitForType,
                self.addToSpinSys,
                self.removeFromSpinSys,
                self.swapProchiral,
                self.ambiguateProchiral]
    buttonList = ButtonList(guiFrame, texts=texts, grid=(row,0), tipTexts=tipTexts,
                            commands=commands, expands=True, gridSpan=(1,10))

    self.assignButton         = buttonList.buttons[1]
    self.typeButton           = buttonList.buttons[2]
    self.assSpinSysButton     = buttonList.buttons[3]
    self.swapProchiralButton  = buttonList.buttons[5]
    self.ambigProchiralButton = buttonList.buttons[6]

    row += 1
    tipText = 'General functions to administer resonances and get connected peak information'
    div = LabelDivider(guiFrame, text='General & Information',
                       grid=(row,0), tipText=tipText)
    
    row += 1
    tipTexts = ['Delete the selected resonances form the CCPN project',
                'Delete all resonances that have no assignments to peaks',
                'Merge multiple resonances & any peak assignments into one. For a single resonance split into two duplicate resonances',
                'Show a table of all the peaks to which the selected resonances are assigned',
                'Considering only the selected shift list, show a table of peaks to which selected resonances are assigned',
                'Show a table of information for the last selected resonance, including a list of all peak dimension positions']
    texts    = ['Delete','Delete\nOrphans',
                'Merge or\nSplit','Show All\nPeaks',
                'Show Shift\nList Peaks','Info']
    commands = [self.deleteResonance,self.deleteOrphans,
                self.mergeOrSplitResonances,self.showPeaks,
                self.showShiftListPeaks,self.showInfo]
    self.bottomButtons = ButtonList(guiFrame, texts=texts, commands=commands,
                                    grid=(row,0), gridSpan=(1,10), tipTexts=tipTexts)
    self.mergeButton = self.bottomButtons.buttons[2]
    
    self.curateNotifiers(self.registerNotify)
    
    
    self.refresh = False 
    self.updateAfter()
  
  def getWindows(self):

    windowPanes = []
    for window in self.analysisProject.sortedSpectrumWindows():
      for windowPane in window.spectrumWindowPanes:
        if self.resonance and (self.isotope != 'Any'):
          axisPanel = windowPane.findFirstAxisPanel(label='x')
          
          if axisPanel: # should always be ok
            axisType = axisPanel.axisType
            if axisType and (self.isotope in axisType.isotopeCodes):
              windowPanes.append([getWindowPaneName(windowPane), windowPane])
        else:
          windowPanes.append([getWindowPaneName(windowPane), windowPane])
 
    windowPanes.sort()
    
    return [x[1] for x in windowPanes]

  def updateWindowsAfter(self, window=None):
    
    self.after_idle(lambda :self.updateWindows(window)) # Axis panel notifiers need a bit of time
 
 
  def updateWindows(self, *window):
      
    windowPanes = self.getWindows()
    index = 0
    names = []
    windowPane  = self.windowPane
 
    if windowPanes:
      names  = [getWindowPaneName(p) for p in windowPanes]
     
      if windowPane not in windowPanes:
        windowPane = windowPanes[0]
        
      index = windowPanes.index(windowPane) 
 
    if windowPane is not self.windowPane:
      self.windowPane = windowPane
 
    self.windowPulldown.setup(names, windowPanes, index)
  
  
  def getResWindows(self):

    windowPanes2D = []
    windowPanesND = []
    for window in self.analysisProject.spectrumWindows:
      for windowPane in window.spectrumWindowPanes:
        if windowPaneHasValueAxis(windowPane):
          continue
        n = len(windowPane.axisPanels)
        if n == 2:
          windowPanes2D.append((getWindowPaneName(windowPane), windowPane))
        elif n > 2:
          windowPanesND.append((getWindowPaneName(windowPane), windowPane))

    windowPanes2D.sort()
    windowPanesND.sort()

    return [x[1] for x in windowPanes2D], [x[1] for x in windowPanesND]

  def updateResWindowsAfter(self, *windowPane):

    self.after_idle(lambda :self.updateResWindows()) # Axis panel notifiers need a bit of time

  def updateResWindows(self, *windowPane):

    windowPanes2D, windowPanesND = self.getResWindows()

    for windowPanes, windowPulldown in ((windowPanesND, self.windowRes1Pulldown),
                                        (windowPanesND, self.windowRes2Pulldown),
                                        (windowPanes2D, self.windowRes3Pulldown)):
      names  = ['<None>'] + [getWindowPaneName(p) for p in windowPanes]
      windowPanes = [None] + windowPanes

      windowPane  = windowPulldown.getObject()
      if windowPane not in windowPanes:
        windowPane = None

      index = windowPanes.index(windowPane)
      windowPulldown.setup(names, windowPanes, index)

  def changeWindow(self, windowPane):

    if windowPane is not self.windowPane:
      self.windowPane = windowPane
  
   
  def stripSelected(self):
  
    if self.resonance and self.windowPane and self.shift:
      shiftList = self.shift.parentList
      stripAxis  = self.windowPane.spectrumWindow.stripAxis 
      shifts = [x[1] for x in self.scrolledMatrix.currentObjects if x[1]]
      
      stripAxisPanel = self.windowPane.findFirstAxisPanel(label=stripAxis)
      isotopes = stripAxisPanel.axisType.isotopeCodes
      stripShifts = [s for s in shifts if s.resonance.isotopeCode in isotopes]
      
      positions  = [{stripAxis:s.value} for s in stripShifts] 
      
      if not positions:
        msg = 'None of the selected resonances match the strip axis (%s)' % stripAxis
        showWarning('Abort', msg, parent=self)
        return 

      
      displayStrips(self.parent, positions,
                    windowPane=self.windowPane)
      
      stripResonances = [s.resonance for s in stripShifts]
      otherResonances = set([s.resonance for s in shifts if s not in stripShifts])
      
      prevActive = stripAxisPanel.findFirstAxisRegion(isActive=True)
      regions = stripAxisPanel.sortedAxisRegions()
      for region in regions:
        region.isActive = False
      
      for i, resonance in enumerate(stripResonances):
        bound = list(resonance.covalentlyBound or [])
        bound.sort(key=makeResonanceGuiName)
        resonances2 = set(bound) & otherResonances
        
        resonance2 = None
        
        if not resonances2:
          resonances2 = bound
        
        regions[i].isActive = True
        
        for resonance2 in resonances2:
        
          gotoResonancePosition(self.guiParent, resonance2,
                                self.windowPane, shiftList)
          
        gotoResonancePosition(self.guiParent, resonance,
                              self.windowPane, shiftList)
                              
        regions[i].isActive = False
      
      if prevActive:
        prevActive.isActive = True
          
 
  def gotoPosition(self):
  
    if self.resonance and self.windowPane and self.shift:
      shiftList = self.shift.parentList
      gotoResonancePosition(self.guiParent, self.resonance,
                            self.windowPane, shiftList)
      
      bound = list(self.resonance.covalentlyBound) or []
      if len(bound) == 1:
        gotoResonancePosition(self.guiParent, bound[0],
                              self.windowPane, shiftList)

  def gotoResonancesPosition(self, windowPane, resonances, **rowColDict):

    if self.shift:
      shiftList = self.shift.parentList
      for resonance in resonances:
        gotoResonancePosition(self.guiParent, resonance,
                              windowPane, shiftList, **rowColDict)

  def markSelected(self):

    if self.resonance and self.shift:
    
      shifts = [x[1] for x in self.scrolledMatrix.currentObjects if x[1]]
      shifts = [(s,s.resonance.isotopeCode) for s in shifts]
      
      panelTypeDict = {}
      
      if self.windowPane:
        panelTypes = [ap.panelType for ap in self.windowPane.axisPanels if ap.panelType]
      else:
        panelTypes = self.analysisProject.panelTypes
        
      for panelType in panelTypes:
        for isotope in panelType.axisType.isotopeCodes:
          if panelTypeDict.get(isotope) is None:
             panelTypeDict[isotope] = []
          
          panelTypeDict[isotope].append(panelType)   

      for shift, isotope in shifts:
        panelTypes = panelTypeDict.get(isotope, [])
        
        for panelType in panelTypes:
          ruler = createRuler(shift.value, panelType, remove=False)
          self.rulers.append(ruler)
 
  def clearMarks(self):
  
    for ruler in self.rulers:
      if not ruler.isDeleted:
        ruler.delete()
        
  def stripResidueResonances(self):
    
    if not self.resonance:
      return

    residue = getResonanceResidue(self.resonance)
    if not residue:
      return

    resonancePairs = []
    carbonResonances = getResidueResonances(residue, atomType='C')
    for carbonResonance in carbonResonances:
      (molSystem, chain, residue, name) = getResonanceAtomTuple(carbonResonance)
      protonResonances = getOnebondResonances(carbonResonance, isotopeCode='1H')
      if len(protonResonances) > 0:
        for protonResonance in protonResonances:
        # wb104: not sure why we need to do the below but Igor's macro did it
        # (molSystem1, chain1, residue1, name1) = getResonanceAtomTuple(protonResonance)
        # if name=='Cga' or name=='Cgb' or name=='Cda' or name=='Cdb':
        #   if name[1:3]!=name1[1:3]:
        #     continue
          resonancePairs.append([carbonResonance,protonResonance])
      else:
        resonancePairs.append([carbonResonance])
    for windowPulldown in (self.windowRes1Pulldown, self.windowRes2Pulldown, self.windowRes3Pulldown):
      windowPane = windowPulldown.getObject()
      if windowPane:
        windowFrame = windowPane.getWindowFrame()
        numStrips = len(resonancePairs)
        if windowPane.spectrumWindow.stripAxis == 'x':
          windowFrame.setColsTo(numStrips)
          attr = 'col'
        else:
          windowFrame.setRowsTo(numStrips)
          attr = 'row'
        for strip, resonancePair in enumerate(resonancePairs):
          rowColDict = {attr: strip}
          self.gotoResonancesPosition(windowPane, resonancePair, **rowColDict)
  
  def markResidueResonances(self):

    if not self.resonance:
      return

    residue = getResonanceResidue(self.resonance)
    if not residue:
      return

    carbonResonances = getResidueResonances(residue, atomType='C')
    for carbonResonance in carbonResonances:
      self.markResonance(carbonResonance)
      protonResonances = getOnebondResonances(carbonResonance, isotopeCode='1H')
      for protonResonance in protonResonances:
        self.markResonance(protonResonance)
  
  def clearResidueMarks(self):

    for ruler in self.residueRulers:
      if not ruler.isDeleted:
        ruler.delete()
  
  def markResonance(self, resonance):

    if not self.shiftList:
      return

    shift= resonance.findFirstShift(parentList=self.shiftList)
    if not shift:
      return

    for windowPulldown in (self.windowRes1Pulldown, self.windowRes2Pulldown, self.windowRes3Pulldown):
      windowPane = windowPulldown.getObject()
      if not windowPane:
        continue
      panelTypes = [ap.panelType for ap in self.windowPane.axisPanels if ap.panelType]
      panelTypeDict = {}

      for panelType in panelTypes:
        for isotope in panelType.axisType.isotopeCodes:
          if panelTypeDict.get(isotope) is None:
            panelTypeDict[isotope] = []
          panelTypeDict[isotope].append(panelType)

      panelTypes = panelTypeDict.get(shift.resonance.isotopeCode, [])

      for panelType in panelTypes:
        ruler = createRuler(shift.value, panelType, remove=False)
        self.residueRulers.append(ruler)

  def updateButtons(self):
  
    self.stripButton.disable()
    self.markButton.disable()
    self.swapProchiralButton.disable() 
    self.ambigProchiralButton.disable() 
    if self.resonance:
      resonanceSet = self.resonance.resonanceSet
      if resonanceSet and (len(resonanceSet.atomSets) == 1):
        chemAtom = resonanceSet.findFirstAtomSet().findFirstAtom().chemAtom
        if chemAtom.chemAtomSet:
          self.swapProchiralButton.enable()
          self.ambigProchiralButton.enable()
      
      if self.shift:
        self.stripButton.enable()
        self.markButton.enable()
      
      resonances = self.scrolledMatrix.currentObjects
      N = len(resonances)
      if N == 0:
        self.mergeButton.config(text='Merge or\nSplit')    
      elif N == 1:
        self.mergeButton.config(text='Split')    
      elif N > 1:
        self.mergeButton.config(text='Merge')    
  
  def swapProchiral(self):
  
    if self.resonance:
      resonances = [x[0] for x in self.scrolledMatrix.currentObjects]
      if len(resonances) > 1:
        msg = 'Really swap all prochirals within selection?'
        if not showOkCancel('Confirm', msg, parent=self):
          return
      
      for resonance in resonances:
        swapProchiralResonance(resonance)
  
  def ambiguateProchiral(self):
  
    if self.resonance:
      resonances = [x[0] for x in self.scrolledMatrix.currentObjects]
      msg = 'Really make all prochirals within selection ambiguous?'
      if not showOkCancel('Confirm', msg, parent=self):
        return
    
      for resonance in resonances:
        swapProchiralResonance(resonance, makeAmbiguous=True)
  
  def doEditMarkExtraRules(self, obj, row, col):
  
    resonance, shift = obj
    if (col == 1) and resonance and shift:
      shiftList = shift.parentList
      for contrib in resonance.peakDimContribs:
        peak = contrib.peakDim.peak
        if shiftList is peak.peakList.dataSource.experiment.shiftList:
          return False
                  
    return True
  
  def getShiftValue(self, obj):
    
    resonance, shift = obj
    if shift:
      self.shiftValueEntry.set(shift.value)
  
  def setShiftValue(self, *whatever):
  
    if self.shift:
      value = self.shiftValueEntry.get()
      self.shift.setValue(value)
  
  def getName(self, obj):

    resonance, shift = obj
    if resonance:
      self.nameEntry.set(resonance.name)
  
  def setName(self, *whatever):
  
    if self.resonance:
      value = self.nameEntry.get().strip() or None
      self.resonance.setName(value)
        
  def getDetails(self, resonanceShift=None):
    
    resonance, shift = resonanceShift
    if resonance:
      self.detailsEntry.set(resonance.details or '')
  
  def setDetails(self, *args):
    
    if self.resonance:
      details = self.detailsEntry.get().strip() or None
      if details != self.resonance.details:
        self.resonance.details = details
        updateResonanceAnnotation(self.resonance)
      
  def updateChains(self):
  
    if self.resonances is not None:
      return
      
    chains = self.getChains()
    names = [self.getChainName(x) for x in chains]
    chain = self.chain
    index = 0
    
    if chains:
      if chain not in chains:
        chain = chains[0]

      index = chains.index(chain)
    
    else:
      chain = None
    
    if chain is not self.chain:
      self.changeChain(chain)
    
    self.chainPulldown.setup(names, chains, index)

  def updateCcpCodes(self):
  
    ccpCodes = [None,] + self.getCcpCodes()
    ccpCode = self.ccpCode

    if ccpCode not in ccpCodes:
      ccpCode = None
    
    index = ccpCodes.index(ccpCode)
    
    if self.ccpCode is not ccpCode:
      self.changeCcpCode(ccpCode)      
    
    self.ccpCodePulldown.setup(['<Any>',] + ccpCodes[1:], ccpCodes, index)
    
  def changeChain(self, chain):
  
    if chain is not self.chain:
      self.chain = chain
      self.updateAfter()
  
  def changeCcpCode(self, ccpCode):
  
    if ccpCode != self.ccpCode:
      self.ccpCode = ccpCode
      self.updateAfter()
  
  def getChainName(self, chain):
  
    if chain is None:
      return '<Any>'
    
    data = (chain.molSystem.code, chain.code, chain.molecule.molType)
    return '%s:%s-%s' % data
    
  def getChains(self):
  
    chains = [None,]
    if self.project:
      for molSystem in self.project.sortedMolSystems():
        for chain in molSystem.sortedChains():
          chains.append(chain)
          
    return chains
  
  def getCcpCodes(self):
  
    ccpCodes = set()
    if self.chain:
      for residue in self.chain.residues:
        ccpCodes.add(getResidueCode(residue))
    
    else:
      for chain in self.getChains():
        if chain:
          for residue in chain.residues:
            ccpCodes.add(getResidueCode(residue))
    
    codes = list(ccpCodes)
    codes.sort()
    
    return codes
  
  def cancelSpinSystemWait(self):
  
    self.waitingForSpinSystem = None
    self.assSpinSysButton.config(text='Add to\nSpin System')
    self.assSpinSysButton.config(**NORM_CONFIG)  
  
  def addToSpinSys(self):
  
    self.cancelAtomsWait()
    self.cancelTypeWait()

    if self.waitingForSpinSystem:
      self.cancelSpinSystemWait()
    else:
      self.waitingForSpinSystem = True
      self.guiParent.editSpinSystems()
      self.assSpinSysButton.config(**CANCEL_CONFIG)
      self.assSpinSysButton.config(text='Cancel Add')
      
  def chooseSpinSystem(self, spinSystem):
  
    if self.resonance and self.waitingForSpinSystem:
      for resonance, shift in self.scrolledMatrix.currentObjects:
        addSpinSystemResonance(spinSystem,resonance)
  
  def removeFromSpinSys(self):
  
    objs = self.scrolledMatrix.currentObjects
    if objs:
      for resonance, shift in objs:
        spinSystem = resonance.resonanceGroup
        if spinSystem:
          deassignResonance(resonance, clearAssignNames=False)
          removeSpinSystemResonance(spinSystem, resonance)
    
  def changeStatus(self, status):
  
    if self.status != status:
      self.status = status
      self.updateAfter()
 
  def changeIsotope(self, isotope):
  
    if isotope != self.isotope:
      self.isotope = isotope
      self.updateAfter()
 
  def selectCell(self, obj, row, col):
 
    resonance, shift = obj
    self.resonance = resonance
    self.shift = shift
    self.updateButtons()
    self.updateResWindows()
            
  def getShiftListNames(self, shiftLists):
    
    shiftListNames = []
    for shiftList in shiftLists:
      name = shiftList.name or '<No name>'
      name = '%s [%d]' % (name, shiftList.serial)
      shiftListNames.append(name)

    return shiftListNames
 
  def changeShiftList(self, shiftList):
    
    if self.shiftList is not shiftList:
      self.shiftList = shiftList
      self.updateAfter()
  
  def deleteOrphans(self):
  
    self.cancelAllWaits()
    
    if self.shiftList:
      msg = 'Use only the selected shift list (#%d)? ' % self.shiftList.serial
      msg += 'Otherwise all lists will be considered' 
      response = showYesNo('Delete Orphans', msg, parent=self)
    else:
      response = False
    
    if response:
      shiftLists = [self.shiftList,]
    else:
      shiftLists = getShiftLists(self.nmrProject)
    
    orphans    = []
    for shiftList in shiftLists:
      for shift in shiftList.measurements:
        resonance = shift.resonance
        
        for contrib in resonance.peakDimContribs:
          if contrib.peakDim.peak.peakList.dataSource.experiment.shiftList is shiftList:
            break
        else:
          orphans.append(shift)
 
    if len(orphans) > 0:
      msg = 'Are you sure you want to\ndelete %d orphaned shifts?'
      if showOkCancel('Delete Orphans', msg % len(orphans) , parent=self):
        doomed = set()
        
        for shift in orphans:
          resonance = shift.resonance
          shift.delete()
        
          if not resonance.shifts:
            doomed.add(resonance)
        
        msg = 'Also delete %d resonances which now have no shifts?' % len(doomed)
        if doomed and showOkCancel('Delete Resonances', msg, parent=self):
          for resonance in doomed:
            resonanceSet = resonance.resonanceSet
            if resonanceSet:
              if len(resonanceSet.resonances) == 1:
                resonanceSet.delete()
              else:
                resonanceSet.removeResonance(resonance)

            resonance.delete()
        
        self.updateAfter()

  def mergeOrSplitResonances(self):
    
    resonances = [x[0] for x in self.scrolledMatrix.currentObjects]
    if len(resonances) == 1:
      msg = 'Are you sure you want to split the selected resonance?'
      if showOkCancel('Split resonance', msg, parent=self):
        resonanceB = splitResonance(resonances[0])
        n = len(resonanceB.peakDimContribs)
        msg = 'Split spawned new resonance [%d]; assigned to %d peak dimensions'
        showInfo('Info', msg % (resonanceB.serial, n), parent=self)
  
    elif len(resonances) > 1:
      msg =  'Are you sure you want to\nmerge %d resonances?'
      if showOkCancel('Merge resonances', msg % len(resonances) , parent=self):
        resonance1 = resonances[0]
        for resonance2 in resonances[1:]:
            mergeResonances(resonance2, resonance1)
      
  def cancelAtomsWait(self):
  
    self.waitingForAtoms = False
    self.assignButton.config(**NORM_CONFIG)
    self.assignButton.config(text = 'Assign')

  def cancelTypeWait(self):
  
    self.waitingForType = False
    self.typeButton.config(**NORM_CONFIG)
    self.typeButton.config(text = 'Set Atom\nType')
      
  def waitForType(self):

    self.cancelSpinSystemWait()
    self.cancelAtomsWait()
    if self.waitingForType:
      self.cancelTypeWait()
    else:
      self.typeButton.config(**CANCEL_CONFIG)
      self.typeButton.config(text = 'Cancel\nSet Type')
      self.guiParent.browseAtoms(requestor=self)
      self.waitingForType = True

  def waitForAtoms(self):

    self.cancelSpinSystemWait()
    self.cancelTypeWait()
    if self.waitingForAtoms:
      self.cancelAtomsWait()
    else:
      self.assignButton.config(**CANCEL_CONFIG)
      self.assignButton.config(text = 'Cancel\nAssign')
      self.guiParent.browseAtoms(requestor=self)
      self.waitingForAtoms = True
  
  def chooseAtoms(self, atomSetMapping):
    
    if self.waitingForAtoms:
      self.assignResonance(atomSetMapping)
    
    elif self.waitingForType:
      self.assignAtomType(atomSetMapping)
    
    self.cancelAllWaits()
    
  
  def assignAtomType(self, atomSetMapping):
  
    if self.waitingForType:
      resonances = [x[0] for x in self.scrolledMatrix.currentObjects]
      if resonances:
        if atomSetMapping.mappingType == 'ambiguous':
          msg = 'Cannot set atom type with an ambiguous atom selection'
          showWarning('Failure', msg, parent=self)
          self.waitingForType = False
          self.typeButton.config(**NORM_CONFIG)
          return

        failed   = 0
        atomSets = atomSetMapping.atomSets
        for resonance in resonances:
          element = resonance.isotope.chemElement.symbol
          if atomSetMapping.elementSymbol != element:
            failed += 1
            continue
 
          assignResonanceType(resonance,atomSets)
        
        if failed:
          msg = 'Some resonances (%d) were ignored '
          msg += 'because they had a conflicting isotope type'
          showWarning('Warning', msg % failed, parent=self)
            
    
  def assignResonance(self, atomSetMapping):
  
    if self.waitingForAtoms:
      if self.resonance:
        element = self.resonance.isotope.chemElement.symbol
        if atomSetMapping.elementSymbol == element:
          if atomSetMapping.mappingType == 'ambiguous':
            msg = 'Cannot assign resonance to\nambiguous atom selection'
            showWarning('Assign failed', msg, parent=self)
            self.waitingForAtoms = False
            self.assignButton.config(**NORM_CONFIG)
            return
 
          serials = list(atomSetMapping.resonanceSerials)
          serial  = self.resonance.serial
 
          atomSets = atomSetMapping.atomSets
          resonanceSet = findResonanceSet(self.resonance,atomSets)
          for atomSet in atomSets:
            if not self.checkSpinSystem(self.resonance, atomSet):
              return
              
          resonanceSet = assignAtomsToRes(atomSets,self.resonance,resonanceSet)


  def assignResidue(self,residue):
  
    if self.waitingForAtoms:
      if self.resonance:
        assignResonanceResidue(self.resonance,residue)

      self.waitingForAtoms = False
      self.assignButton.config(**NORM_CONFIG)

    
  def disconnectAtoms(self):
  
    self.cancelAllWaits()
    resonances = [x[0] for x in self.scrolledMatrix.currentObjects]
    if len(resonances) == 1:
      if self.resonance.resonanceSet:
        if showOkCancel('Confirm', 'Really deassign resonance?', parent=self):
          deassignResonance(self.resonance)

      elif self.resonance.assignNames:
        msg = 'Really remove atom type information?'
        if showOkCancel('Confirm', msg, parent=self):
          assignResonanceType(self.resonance, None)
    
    else:
      msg = 'Really deassign %d resonances?'
      msg2 = 'Remove atom type information instead?'
      if showOkCancel('Confirm', msg % len(resonances), parent=self):
        for resonance in resonances:
          if resonance.resonanceSet:
            deassignResonance(resonance)
 
      elif showOkCancel('Confirm', msg2, parent=self):
        for resonance in resonances:
          assignResonanceType(resonance, None)
    
    for resonance in resonances:
      if resonance.peakDimContribs:
        msg = 'Remove resonance assignments to peaks?'
        if showYesNo('Confirm', msg, parent=self):
          for resonance1 in resonances:
            clearResonancePeakDimContribs(resonance1)
        break
    
  def deleteResonance(self):
  
    self.cancelAllWaits()
    resonances = list(set([x[0] for x in self.scrolledMatrix.currentObjects]))
    if len(resonances) == 1:
      resonance = resonances[0]
      n = len(resonance.peakDimContribs)
      if n > 0:
        msg = 'Resonance still assigned to %d peaks'
        showWarning('Delete failed', msg % n, parent=self)
        return
      else:
        msg = 'Delete resonance %d?' % resonance.serial
        if showOkCancel('Delete Resonance', msg, parent=self):
          resonanceSet = resonance.resonanceSet
          if resonanceSet:
            if len(resonanceSet.resonances) == 1:
              resonanceSet.delete()
            else:
              resonanceSet.removeResonance(resonance)
          
          if not resonance.isDeleted:
            resonance.delete()
            
          self.updateAfter()
 
    elif  len(resonances) > 1:
      msg = 'Delete %d resonances' % len(resonances)
      if showOkCancel('Delete Resonance', msg, parent=self):
        toDelete = []
        stillAssigned = 0
        
        for r in resonances:
          if len(r.peakDimContribs) > 0:
            stillAssigned += 1
          elif not r.isDeleted:
            toDelete.append(r)
            
        for resonance in toDelete:
          resonanceSet = resonance.resonanceSet
          if resonanceSet:
            if len(resonanceSet.resonances) == 1:
              resonanceSet.delete()
            else:
              resonanceSet.removeResonance(resonance)
              
          resonance.delete()
        
        if stillAssigned > 0:
          msg = '%d resonances not deleted: still assigned to peaks' 
          showWarning('Delete failed', msg % stillAssigned, parent=self)
        
        self.updateAfter()   

  def showPeaks(self, resonanceShift=None):
  
    self.cancelAllWaits()
    
    if resonanceShift:
      resonances = [resonanceShift[0],]
    else:
      resonances = [x[0] for x in self.scrolledMatrix.currentObjects]
    
    if resonances:
      peaks = set()
      for resonance in resonances:
        for contrib in resonance.peakDimContribs:
          peaks.add(contrib.peakDim.peak)
      
      if len(peaks) > 0:
        self.guiParent.viewPeaks(list(peaks))

  def showShiftListPeaks(self, resonanceShift=None):
  
    self.cancelAllWaits()
    
    if resonanceShift:
      resonances = [resonanceShift[0],]
    else:
      resonances = [x[0] for x in self.scrolledMatrix.currentObjects]
    
    if resonances:
      peaks = set()
      
      for resonance in resonances:
        for contrib in resonance.peakDimContribs:
          peak = contrib.peakDim.peak
          
          if self.shiftList is None:
            peaks.add(peak)
          
          else:
            experiment = peak.peakList.dataSource.experiment
            if experiment.shiftList is self.shiftList:
               peaks.add(peak)
          
      if len(peaks) > 0:
        self.guiParent.viewPeaks(list(peaks))

  def showInfo(self):
    
    if self.resonance:
      popup = ResonanceInfoPopup(self.guiParent, self.resonance)
      popup.open()
      self.guiParent.popups['resonance_info'] = popup

  def cancelAllWaits(self):

    self.cancelSpinSystemWait()
    self.cancelAtomsWait()
    self.cancelTypeWait()
 
  def updateShiftLists(self, obj=None):
  
    if self.resonances:
      if obj is self.shiftList:
        self.updateAfter()
      return
  
    index = 0
    shiftLists = getShiftLists(self.nmrProject) + [False, None]
    shiftListNames = self.getShiftListNames(shiftLists[:-2]) + ['None','<Any>']
    shiftList = self.shiftList
    
    if shiftListNames:
      if shiftList not in shiftLists:
        shiftList = shiftLists[0]
      
      index = shiftLists.index(shiftList)
      
    else:
      shiftList = None
    
    if self.shiftList is not shiftList:
      self.changeShiftList(shiftList)
    
    self.shiftListPulldown.setup(shiftListNames,shiftLists,index)

  def close(self):
  
    self.cancelAllWaits()
    BasePopup.close(self)
 
  def curateNotifiers(self, notifyFunc):
  
    notifyFunc(self.updateShiftLists, 'ccp.nmr.Nmr.ShiftList', 'setName')
    notifyFunc(self.updateShiftLists, 'ccp.nmr.Nmr.ShiftList', '__init__')
    notifyFunc(self.updateShiftLists, 'ccp.nmr.Nmr.ShiftList', 'delete')
 
    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.PeakDimContrib',
                    'ccp.nmr.Nmr.Shift','ccp.nmr.Nmr.ResonanceSet',
                    'ccp.nmr.Nmr.Resonance'):
        notifyFunc(self.updateAfter, clazz, func)
        
    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.ResonanceSet', 'addResonance')
    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.Shift', 'setValue')
    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.Shift', 'setError')

    for func in ('delete','__init__','setWeight','setPossibility'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.ResidueProb', func)
    
    for func in ('setCcpCode', 'addResonance', 'setName',
                 'removeResonance', 'setResonances',
		 'addAtomSet', 'delete'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)

    for func in ('setDetails', 'setIsotopeCode', 'setName', 'setResonanceGroup',
                 'setAssignNames', 'addAssignName', 'removeAssignName'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.Resonance', func)

    for func in ('__init__', 'delete','setName'):
      notifyFunc(self.updateWindowsAfter, 'ccpnmr.Analysis.SpectrumWindow', func)
      notifyFunc(self.updateResWindowsAfter, 'ccpnmr.Analysis.SpectrumWindow', func)

    #notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.Shift', 'getResonance')
 
  def destroy(self):
  
    self.curateNotifiers(self.unregisterNotify)

    BasePopup.destroy(self) 
  
  def updateResonanceAfter(self, resonance):
    
    if (self.resonances is not None) and (resonance in self.resonances):
      self.refresh = True
      self.after_idle(self.update)
    
    elif resonance.findFirstShift(parentList=self.shiftList):
      self.refresh = True
      self.after_idle(self.update)
      
    elif self.shiftList is None:
      self.refresh = True
      self.after_idle(self.update)

  
  def updateResonanceSetAfter(self, resonanceSet):
     
    for resonance in resonanceSet.resonances:
      if self.resonances and (resonance in self.resonances):
        break
      elif resonance.findFirstShift(parentList=self.shiftList):
        break
    
    else:
      # No reason to update
      return
    
    self.refresh = True
    self.after_idle(self.update)  
  
  def updateShiftAfter(self, shift):
  
    if (self.resonances is not None) and (shift.resonance in self.resonances):
      self.refresh = True
      self.after_idle(self.update)

    elif shift.parentList is self.shiftList:
      self.refresh = True
      self.after_idle(self.update)      

    elif self.shiftList is None:
      self.refresh = True
      self.after_idle(self.update)      
        
  def updatePeakDimContribAfter(self, contrib):

    resonance = contrib.resonance
    if self.resonances and (resonance in self.resonances):
      self.refresh = True
      self.after_idle(self.update)

    elif self.shiftList is None:
      self.refresh = True
      self.after_idle(self.update)
    
    elif resonance.findFirstShift(parentList=self.shiftList):
      self.refresh = True
      self.after_idle(self.update)
           
  def updateResonanceGroupAfter(self, spinSystem):

    for resonance in spinSystem.resonances:
      if self.resonances and (resonance in self.resonances):
        break
      elif resonance.findFirstShift(parentList=self.shiftList):
        break
      elif self.shiftList is None:
        break
    
    else:
      # No reason to update
      return
    
    self.refresh = True
    self.after_idle(self.update)  
  
  def updateResidueProbAfter(self, residueProb):

    for resonance in residueProb.resonanceGroup.resonances:
      if self.resonances and (resonance in self.resonances):
        break
      elif resonance.findFirstShift(parentList=self.shiftList):
        break
      elif self.shiftList is None:
        break
    
    else:
      # No reason to update
      return
    
    self.refresh = True
    self.after_idle(self.update)  

  def updateAfter(self, obj=None):
    
    if self.refresh:
      return

    funcs = {'Resonance':self.updateResonanceAfter,
             'ResonanceSet':self.updateResonanceSetAfter,
             'Shift':self.updateShiftAfter,
             'PeakDimContrib':self.updatePeakDimContribAfter,
             'ResonanceGroup':self.updateResonanceGroupAfter,
             'ResidueProb':self.updateResidueProbAfter}

    if obj:
      func = funcs[obj.className]
      func(obj)

    else:        
      self.refresh = True
      self.after_idle(self.update)
    
  def update(self, resonances=None, shiftList=None):

    self.cancelAllWaits()
    self.updateWindows()
    self.updateResWindows()

    if resonances and (self.resonances is not None):
      self.resonances = resonances
      if shiftList:
        self.shiftList = shiftList
    
    if self.resonances is None:
      self.updateChains()
      self.updateCcpCodes()
      self.updateShiftLists()
    
    analysisProject = self.analysisProject
    doChain   = analysisProject.doChainAnnotations
    doMolSys  = analysisProject.doMolSysAnnotations
    project   = self.project
    status    = self.status
    chain     = self.chain
    isotope   = self.isotope
    ccpCode   = self.ccpCode
    shiftList = self.shiftList
    
    resonances = []
    resonancesAppend = resonances.append
    shifts  = []
    shiftsAppend = shifts.append
    textMatrix = []
    textMatrixAppend = textMatrix.append
    colorMatrix = []
    colorMatrixAppend = colorMatrix.append
    tryResonances = []
    tryResonancesAppend = tryResonances.append    
    tryShifts  = []
    tryShiftsAppend = tryShifts.append
    
    if self.resonances:
      tryResonances = [r for r in self.resonances if not r.isDeleted]
      self.resonances = tryResonances
      if shiftList:
        tryShifts = [[r.findFirstShift(parentList=shiftList)] for r in tryResonances]
      elif shiftList is None:
        for resonance in tryResonances:
          tryShiftsAppend(resonance.sortedShifts())
          
      else:
        tryShifts = [[r.findFirstShift()] for r in tryResonances]
    
    elif shiftList:
      for shift in shiftList.sortedMeasurements():
        tryResonancesAppend(shift.resonance)
        tryShiftsAppend([shift,])

    elif shiftList is None:
      tryResonances = self.nmrProject.sortedResonances()
      for resonance in tryResonances:
        tryShiftsAppend(resonance.sortedShifts())

    else:
      for resonance in self.nmrProject.sortedResonances():
        if not resonance.shifts:
          tryResonancesAppend(resonance)
          tryShiftsAppend([None,])
    
    for i, resonance in enumerate(tryResonances):
      shiftsB = tryShifts[i]
      
      if chain:
        residue = getResonanceResidue(resonance)
        if not residue or (residue.chain is not chain):
          continue
        
      if ccpCode:
        residue = getResonanceResidue(resonance)
        if residue:
          if getResidueCode(residue) != ccpCode:
            continue
          
        else:
          spinSystem = resonance.resonanceGroup
          if (not spinSystem) or (getResidueCode(spinSystem) != ccpCode):
            continue

      isotopeCode = resonance.isotopeCode
      if isotopeCode == 'unknown':
        peakDimContrib = resonance.findFirstPeakDimContrib()
        if peakDimContrib:
          expDimRef = peakDimContrib.peakDim.dataDimRef.expDimRef
          resonance.isotopeCode = expDimRef.isotopeCodes[0]
      
      #print self.isotope, resonance.isotope
      if (isotope == 'Any') or (isotope == isotopeCode):
        if status == 'Any':
          for shift in shiftsB:
            resonancesAppend(resonance)
            shiftsAppend((shift, isotopeCode))
        elif resonance.resonanceSet and (status == 'Assigned'):
          for shift in shiftsB:
            resonancesAppend(resonance)
            shiftsAppend((shift, isotopeCode))
        elif (not resonance.resonanceSet) and (status == 'Unassigned'):
          for shift in shiftsB:
            resonancesAppend(resonance)
            shiftsAppend((shift, isotopeCode))
        elif (len(resonance.peakDimContribs) == 0) and (status == 'Orphaned'):
          for shift in shiftsB:
            resonancesAppend(resonance)
            shiftsAppend((shift, isotopeCode))
    
    colorCache = {}
    colorCacheGet = colorCache.get
    objectList = []
    objectListAppend = objectList.append
    
    for i, resonance in enumerate(resonances):
      shift, isotopeCode = shifts[i]
      objectListAppend( (resonance, shift) )
      
      if shift is None:
        shiftValue = None
        shiftError = None
      else:
        shiftValue = shift.value
        shiftError = shift.error
      
      color = colorCacheGet(isotopeCode)
      
      if not color:
        isotope = resonance.isotope
        if isotope:
          color = ATOM_COLOUR.get(isotope.chemElement.symbol)
          isotopeText = isotopeCode
          colorCache[isotopeCode] = color
        else:
          color = None
          isotopeText = ''
        
      else:
        isotopeText = isotopeCode
      
      data = getResonanceAtomTuple(resonance)
      (molSystem, chain, residueText, resonanceName) = data

      if not doChain:
        chain = ''
      if not doMolSys:
        molSystem = ''
      
      if molSystem:
        residueText = '%s:%s %s' % (molSystem, chain, residueText)
      elif chain:
        residueText = '%s %s' % (chain, residueText)
      
      total = set()
      sList = set()
      for contrib in resonance.peakDimContribs:
        peak = contrib.peakDim.peak
        total.add(peak)
        if shiftList is peak.peakList.dataSource.experiment.shiftList:
          sList.add(peak)
      
      spinSystemSerial = None
      spinSystem = resonance.resonanceGroup
      if spinSystem:
        spinSystemSerial = spinSystem.serial
      
      data = [resonance.serial,
              shiftValue,
              shiftError,
              resonanceName,
              residueText,
              isotopeText,
              resonance.name,
              len(total),
              len(sList),
              spinSystemSerial,
              resonance.details]
 
      textMatrixAppend( data )
      colorMatrixAppend( [color, None, None, None, None,
                           None, None, None, None] )


    self.waitingForAtoms = False
    self.assignButton.config(**NORM_CONFIG)
    self.updateButtons()
    #self.resonance = None

    self.scrolledMatrix.update(objectList=objectList,
                               textMatrix=textMatrix,
                               colorMatrix=colorMatrix)
    self.refresh = False

  def checkSpinSystem(self, resonance, atomSet):
  
    spinSystem = resonance.resonanceGroup
    if spinSystem:
      if not spinSystem.residue:
        return True
      if len(spinSystem.resonances) < 2:
        return True
 
      residue = atomSet.findFirstAtom().residue
      if spinSystem.residue is not residue:
        msg = 'Move resonance %s to different residue %d%s?'
        data = (makeResonanceGuiName(resonance),residue.seqCode,
                getResidueCode(residue))
        if not showOkCancel('Confirm', msg % data, parent=self):
          return False
        
        msg = 'Re-assign all resonances in the same residue?'
        if not showYesNo('Question', msg, self):
          removeSpinSystemResonance(spinSystem, resonance)
    
    return True
          
# TBD:
#   PeakDim positions graph
#   Bound resonances
#   Shifts onto BMRB chart

class ResonanceInfoPopup(BasePopup):
  """
  **Detailed Information About a Resonance's Peak Assignments**
  
  This popup window shows the locations and identities of all of the peak 
  dimensions to which a particular resonance is assigned, together with
  the averaged chemical shift values for the resonances derived from
  these assigned locations

  The upper table lists the positions in the spectra that the resonance is
  assigned to. Thus, the user may see the complement of peaks that the resonance
  contributes to and what the chemical shift positions of these peaks are, in the
  relevant assigned dimension. This display is handy for identifying chemical
  shift outliers that may be indicative of a mistaken assignment.
 
  The lower table shows the averaged chemical shift values that derive from the
  resonance's assignment to the various peak dimensions. It should be noted that
  each peak dimension does not necessarily contribute equally to the shift
  average. The average is weighted by spectrum (and hence peak) dimension.  By
  default this is even, but imprecise or ambiguous peak dimensions are often
  down-weighted and precise ones are up-weighted. The spectrum dimension
  weightings are set via the "Tolerances" tab, "Shift Weighting" column of main
  `Spectra`_ option. 

  .. _`Spectra`: EditSpectrumPopup.html
    
  """

  def __init__(self, parent, resonance=None, *args, **kw):

    self.guiParent = parent
    self.resonance = resonance
    self.shift = None
    self.peakDim = None
    self.waiting = False
    
    title = "Resonance : Resonance Info"
    BasePopup.__init__(self, parent=parent, title=title, **kw)
    
    self.geometry("400x400+150+150")

  def body(self, guiFrame):

    guiFrame.expandGrid(3,0)
    
    tipText='The number and atomic assignment of the resonance record that information is displayed for'
    self.mainLabel = Label(guiFrame, text='', grid=(0,0),
                           tipText=tipText)
    
    tipTexts = ['For the peak assignments selected in the upper table show a table of the full peak records',]
    buttonList = UtilityButtonList(guiFrame, grid=(0,1),
                                   commands=[self.showPeaks,],
                                   texts=['Show Peaks',],
                                   helpUrl=self.help_url,
                                   tipTexts=tipTexts)
    buttonList.buttons[0].config(bg='#B0FFB0')
    
    detailsFrame = Frame(guiFrame, grid=(1,0), gridSpan=(1,2), sticky='ew')
    label = Label(detailsFrame, grid=(0,0), text='Details:')
    text = self.resonance.details or '' if self.resonance else ''
    self.detailsEntry = Entry(detailsFrame, text=text, returnCallback=self.updateResonanceDetails, grid=(0,1), sticky='ew')
    
    tipText='A table of the peak dimensions, and hence positions, to which the current resonance is assigned'
    div = LabelDivider(guiFrame, text='Peak Assignments',
                       gridSpan=(1,2), grid=(2,0), tipText=tipText)

    tipTexts = ['The peak list of the peak to which the current resonance is assigned (experiment:spectrum:list number)',
                'The serial number of the peak, within its peak list, to which the current resonance is assigned ',
                'The dimension number of the peak that the resonance is assigned to',
                'The position of the peak dimension, typically in ppm units',
                'The assignment annotation of the whole peak',
                'The shift list name and number that the peak uses (via its experiment record)']
    
    headingList = ['Peak List','Peak','Dim',
                   'Position','Assignment',
                   'Shift List',]

    self.peakDimTable = ScrolledMatrix(guiFrame, grid=(3,0),
                                       headingList=headingList,
                                       multiSelect=True, gridSpan=(1,2), 
                                       callback=self.selectPeakDim,
                                       tipTexts=tipTexts)

    tipText='A table of the chemical shift values associated with the current resonance'
    div = LabelDivider(guiFrame, text='Chemical Shifts',
                       gridSpan=(1,2), grid=(4,0), tipText=tipText)

    headingList = ['Shift List','Shift','SD']
    tipTexts = ['The shift list that records a chemical shift value for the resonance',
                'The (averaged) value of the chemical shift in this specific shift list',
                'The standard deviation in the averaged chemical shift value']
    self.shiftTable = ScrolledMatrix(guiFrame, grid=(5,0),
                                     headingList=headingList,
                                     multiSelect=False, gridSpan=(1,2), 
                                     callback=self.selectShift,
                                     tipTexts=tipTexts)
   
    self.curateNotifiers(self.registerNotify)
    
    self.updateAfter()
             
  def destroy(self):
  
    self.curateNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)
    
  def selectShift(self, obj, row, col):
  
    self.shift = obj
                                       
  def selectPeakDim(self, obj, row, col):
    
    self.peakDim = obj

  def curateNotifiers(self, notifyFunc):

    notifyFunc(self.updateDetailsEntry, 'ccp.nmr.Nmr.Resonance', 'setDetails')

  def updateContrib(self, peakDimContrib):

    if peakDimContrib.resonance is self.resonance:
      self.updateAfter()

  def updateShift(self, shift):

    if shift.resonance is self.resonance:
      self.updateAfter()
  
  def updateResonanceSet(self, resonanceSet):
    
    if self.resonance in resonanceSet.resonances:
      self.updateAfter()
  
  def updateResonance(self, resonance):
  
    if self.resonance is resonance:
      self.updateAfter()
  
  def updateResidueProb(self, residueProb):
  
    if self.resonance in residueProb.resonanceGroup.resonances:
      self.updateAfter() 
 
  def updateSpinSystem(self, spinSystem):

    if self.resonance in spinSystem.resonances:
      self.updateAfter()

  def updateAfter(self):
  
    if self.waiting:
      return
      
    else:
      self.waiting = True
      self.after_idle(self.update)
      
  def updateDetailsEntry(self, obj):
    
    if self.resonance is obj:
      self.detailsEntry.set(self.resonance.details or '')
      
  def updateResonanceDetails(self, *args, **kw):
    
    if self.resonance:
      details = self.detailsEntry.get().strip() or None  # empty string not allowed by API, hence None
      if details != self.resonance.details:
        self.resonance.details = details
        updateResonanceAnnotation(self.resonance)
      
  def update(self, resonance=None):
    
    self.updateResonanceDetails()
      
    if resonance:
      self.resonance = resonance
      self.updateDetailsEntry(resonance)

    if self.resonance:
      serial = self.resonance.serial
      name = makeResonanceGuiName(self.resonance)
      self.mainLabel.set('Resonance: %d %s' % (serial, name))
    else:
      self.mainLabel.set('Resonance:')

    textMatrix = []
    peakDimContribs = []
    
    if self.resonance:
      peakDimContribs = self.resonance.sortedPeakDimContribs()

      for contrib in peakDimContribs:
        peakDim = contrib.peakDim
        peak = peakDim.peak
        peakList = peak.peakList
        spectrum = peakList.dataSource
        experiment = spectrum.experiment
        shiftList = experiment.shiftList
        peakListData = (experiment.name, spectrum.name, peakList.serial)
        assigned = [c.resonance for c in peakDim.peakDimContribs]
        assigned = [makeResonanceGuiName(r) for r in assigned]
        assigned.sort()
        
        if contrib.peakDimComponent:
          ppm = getPeakDimComponentSplitting(contrib.peakDimComponent)
        else:
          ppm = peakDim.realValue
        
        datum = ['%s:%s:%d' % peakListData,
                 peak.serial,
                 peakDim.dim,
                 ppm,
                 ','.join(assigned),
                 '%s:%d' % (shiftList.name, shiftList.serial)]
        
        textMatrix.append(datum)

    self.peakDimTable.update(textMatrix=textMatrix,
                             objectList=peakDimContribs)

    textMatrix = []
    shifts = []
    
    if self.resonance:
      shifts = self.resonance.sortedShifts()
      
      for shift in shifts:
        shiftList = shift.parentList
        datum = ['%s:%d' % (shiftList.name, shiftList.serial),
                 shift.value,
                 shift.error]
        
        textMatrix.append(datum)

    self.shiftTable.update(textMatrix=textMatrix,
                           objectList=shifts)
  
    self.waiting = False          
  
  def showPeaks(self):
  
    contribs = self.peakDimTable.currentObjects
    
    peaks = set()
    
    for contrib in contribs:
      peaks.add(contrib.peakDim.peak)
    
    if len(peaks) > 0:
      self.guiParent.viewPeaks(list(peaks))

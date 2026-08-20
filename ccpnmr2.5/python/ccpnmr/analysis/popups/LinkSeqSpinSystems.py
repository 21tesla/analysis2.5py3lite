
"""
======================COPYRIGHT/LICENSE START==========================

LinkSeqSpinSystems.py: Part of the CcpNmr Analysis program

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
from math import log
import cPickle

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.ExperimentBasic import getOnebondExpDimRefs, getSpectrumIsotopes
from ccpnmr.analysis.core.ExperimentBasic import getPrimaryDataDimRef, getOnebondDataDims
from ccpnmr.analysis.core.ExperimentBasic import getSeqAssignRefExperiments
from ccpnmr.analysis.core.WindowBasic import getWindowPaneName, displayStrips
from ccpnmr.analysis.core.WindowBasic import getSpinSystemWindowShifts, getPeakDimAxisMapping
from ccpnmr.analysis.core.MarkBasic import createRuler
from ccpnmr.analysis.core.MoleculeBasic import getResidueCode, getLinkedResidue
from ccpnmr.analysis.core.AssignmentBasic import assignSpinSystemResidue, findConnectedSpinSystem, clearSeqSpinSystemLinks, deassignSpinSystem
from ccpnmr.analysis.core.Util import getSpectrumPosContourColor, getAnalysisDataDim
from ccpnmr.analysis.core.PeakBasic import findPositionMatches
from ccpnmr.analysis.core.ChemicalShiftBasic import getShiftsChainProbabilities

from ccpnmr.nexus.NexusBasic import linkSpinSystemInterIntraResonances

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Color import black, hexRepr
from memops.gui.FloatEntry import FloatEntry 
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.LabelDivider import LabelDivider
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showWarning, showOkCancel, showYesNo
from memops.gui.PartitionedSelector import PartitionedSelector
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame


SeqNoesyExperiments = ['H[N]_H.through-space', 'H_H[N].through-space', 
                       'H[{c|n}]_H[N].through-space', 'h[C]_H[N].through-space'
                      ]


DEFAULT_FOCUS_WIDTH = 5.0
DEFAULT_FOCUS_WIDTH_H = 0.5
MATCH_ORDERS = ['Total peak score', 'Average peak score']
PROLYL = set(('Hyp','Pro'))

X_Y = ['x','y']

def LinkSeqSpinSystemsTestMacro(argServer):
  
  popup = LinkSeqSpinSystemsPopup(argServer.parent)
  popup.open()


class LinkSeqSpinSystemsPopup(BasePopup):
  """
  **Find and Link Sequentially Related Protein Backbone Assignments**
  
  This popup window is designed to expedite the process of sequentially
  assigning a polypeptide backbone. The tool is semi-automated, in that
  possibilities for sequentially related peaks and spin systems are
  automatically found and presented, but it is still up to the user to determine
  which are really sequentially linked. Likewise, the tool will show where a
  sequential stretch of assignments best fits with a molecular sequence, but its
  is the users' decision to commit the final choice. A fully automated tool is
  available using the `Automated Seq. Assignment`_ option, although this is best
  suited to relatively complete and clean data sets.

  This system works with pairs of spectra where one spectrum has peaks that are
  through-carbonyl specific; to show the previous (i-1) residue relative to the
  amide (like with HNcoCA/CB), the other spectrum is non-specific in that it
  shows both peaks for the amides current residue (i) and the previous residue
  in the chain (i-1) (link with HNCA/CB). The general principle is that the
  combination of the carbonyl-specific and non-specific experiment enables
  identification of the resonances for the current residue and the residue one
  previous in the sequence. The sequential assignment is made by finding amide
  locations where the "current residue" resonances of one (hopefully
  unambiguously) match the "previous residue" resonances of another. Fitting of
  chains of amides linked in this way to the molecular sequence, by considering 
  chemical shift values, then gives the final residue assignment. 
  
  For any of the spectra to be used by this system they must be peak picked and 
  corresponding positions must be linked to the relevant spin system and
  resonance assignments, albeit in an anonymous way. This setup is readily
  achieved by the `Pick & Assign From Roots`_ option. The peaks of different
  spectra that relate to the  same amide are only known to be from the same
  amide because of assignment. Several different kinds of spectra may be used,
  including those that detect 13C resonances (like HNCO, HNcaCO, HNcoCA, HNCA
  etc) and those that detect 1H resonances (like HNcocaHA, HNHA, etc). This tool
  can match multiple types of spectra at the same time to improve the assignment
  process, for example matching HA, CO, CA and CB positions to minimise
  ambiguity. However, the spectra that are involved with matching 1H positions
  are naturally displayed in a different spectrum window to those using 13C
  positions.

  **User Operation**

  To perform the sequence assignment the user first selects spectra, and 
  spectrum windows to display them in, categorising them as "query"; selecting
  in the upper panel, or "match"; selecting in the lower panel (toggling the
  "Use?" column). The "query" spectra are those that remain at a set position
  according to a selected amide. The "match" spectra are those that the query
  peaks are compared against and represent the, potentially multiple, positions
  of candidates for amides that might be sequence neighbours to the query. Query
  spectra can be of through-carbonyl types (e.g, HNcoCA), in which case matches
  represent the preceding residue in the sequence (an i to i-1 step), or the
  query spectra can be non-specific, where the matches will be to the following
  residue in the sequence (i to i+1). When making sequential links the tool
  automatically knows which sequence direction to make links with.

  The "Options" tab may be viewed to control how the matching of peaks works, 
  but the default values will be acceptable in most instances. Most options
  control how the spectrum window strip regions are made and displayed.
  Occasionally it may be useful to consider the "Filter By Inter/Intra Type"
  options. By default matches between "previous residue" peaks and other
  "previous residue" peaks from different amides are discarded. However, if a
  "current residue" peak overlaps with a "previous residue" peak for the same
  amide (e.g. very similar 13C positions) then the default rule may miss a
  legitimate match. Toggling the "Filter By Inter/Intra Type" option off will
  allow positions with such overlap to be successfully matched, although this may
  increase the number of spurious matches and add ambiguity.

  The actual sequential assignment is performed in the "Spin System Table" tab.
  The general idea is that the user selects a spin system row, whereupon the
  query window navigates to the amide position of that row and the  match window
  displays all of the other amide positions, in as many strips as required, 
  that have the relevant matches to at least some of the query peaks. The
  identity of the potential sequential matches is also listed in the middle 
  "Matched Peak Positions" table. In this table the user may discard certain
  amide positions/strips from consideration, but the objective is to select one
  match position that represents a real sequential neighbour. If such a choice
  can be made [Set Seq Link] sets the appropriate connection between the query
  and matched amide spin systems. Not that this is merely a statement of
  sequence relationship, no residue assignment will be made unless one of the
  amides is assigned, whereupon the assignment will spread along the chain in
  the appropriate direction.

  Once the user has connected sufficient spin systems into regions of sequential
  connectivity  (visible in the "i+1" and "i-1" columns) the next task is to
  associate the linked but unassigned amides with particular residues. Full
  residue assignment is made by selecting a region of sequence in the lower left
  "Sequence Locations" table and clicking [Assign Selected]. This table and the
  "Residue Types" table to the right are designed to show you how the amides may
  fit into the residue sequence. The right hand table shows the probable residue
  types for the query amide, based upon the chemical shifts in the spin system
  (use [Predict Type] to see details of the prediction for the current spin
  system), and any sequentially linked neighbours. The left hand table shows
  where in the sequence the stretch of sequentially linked amides best fits,
  according to residue type predictions. Those residues that have already been
  assigned are coloured blue.

  **Caveats & Tips**

  When working with multiple polypeptide chains, make sure the correct one is
  selected in the "options" tab to give the right sequence.

  The peak-peak match tolerance is set on a per-spectrum basis (for the relevant
  dimension) in the "Tolerance" columns of the "Windows and Spectra" tab. Note
  that this is actually the same attribute as the regular peak assignment
  tolerance.

  It should be noted that mixing of carbonyl-specific and non-specific in the
  query or match spectra is not permitted; otherwise it is not clear in which
  direction the sequential linking is made.

  The ability to use spectra that match in a 15N dimension, e.g. HNcaN, HNcoCAN
  will be added in due course.

  .. _`Automated Seq. Assignment`: AutoBackbonePopup.html
  .. _`Pick & Assign From Roots`: LinkPeakListsPopup.html
  
  """

  def __init__(self, parent, *args, **kw):

    self.waiting = False
    self.refreshSeqSeg = False
    self.guiParent = parent
    
    self.rootPane = None
    self.queryPaneC = None
    self.matchPaneC = None
    self.queryPaneH = None
    self.matchPaneH = None
    self.spinSystem = None
    self.positions = None
    self.shifts = None
    self.chain = None
    self.residues = None
    self.rulers = []
    self.rootRulers = []
    
    self.queryToggleDict = {}
    self.matchToggleDict = {}
    
    self.peakList = None
    
    self.queryPeakLists = set()
    self.matchPeakLists = set()
    self.matches = []
    self.seqSegments = []
    self.match = None
    self.matchDim = 1
    self.project = parent.project
    self.matchOrder = MATCH_ORDERS[0]
    
    
    BasePopup.__init__(self, parent=parent, title='Assignment : Protein Sequence Assignment', **kw)

    self.refNames, self.refNamesCo = self.getRefExpNames()
    self.color = '#808080'

    self.updateChains()
    self.getAppDataOptions()
    self.updateQueryPanes()
    self.updateMatchPanes()
    self.updateQuerySpectra()
    self.updateMatchSpectra()
    self.updateSeqSegmentsAfter()

  def open(self):
  
    BasePopup.open(self)
    self.updateChains()
    self.getAppDataOptions()
    self.updateQueryPanes()
    self.updateSeqSegmentsAfter()

  def close(self):
  
    self.setAppDataOptions()
    BasePopup.close(self)

 
  def body(self, guiFrame):
  
    self.geometry('600x600')

    guiFrame.expandGrid(0,0)
   
    tipTexts = ['Selects which spectra are used as the source for sequential spin system queries and which are used to find match possibilities',
                'Lists the spin systems from the query spectra, how each matches to potential sequence neighbours and fits in a protein sequence',
                'Specifies the various parameters that are used when matching peaks to fins potential sequentially linked spin systems']
    options = ['Windows & Spectra','Spin System Table','Options']
      
    tabbedFrame = TabbedFrame(guiFrame, options=options, tipTexts=tipTexts,
                              callback=self.selectTab, grid=(0,0))
    self.tabbedFrame = tabbedFrame
    frameA, frameB, frameC = tabbedFrame.frames
    
    #
    # Settings
    #
    
    
    row = 0
    tipText = 'The spectra that are the fixed source of peak positions for a spin system, and the windows to display them in'
    div = LabelDivider(frameA, text='Query Windows and Spectra',
                       grid=(row,0), tipText=tipText)
    
    row +=1
    frameA.expandGrid(row,0)
    frame = Frame(frameA, grid=(row,0))
    frame.expandGrid(1,4)
    
    label = Label(frame, text='13C Window:', grid=(0,0))
    tipText = 'Selects which spectrum window is used to display peaks that provide the source 13C positions; to match against potential sequence neighbours'
    self.queryPanePulldownC = PulldownList(frame, self.changeQueryPaneC,
                                           grid=(0,1), tipText=tipText)

    label = Label(frame, text='1H Window:', grid=(0,2))
    tipText = 'Selects which spectrum window is used to display peaks that provide the source 1H positions; to match against potential sequence neighbours'
    self.queryPanePulldownH = PulldownList(frame, self.changeQueryPaneH,
                                           grid=(0,3), tipText=tipText)

    label = Label(frame, text='Root Window:', grid=(0,5))
    tipText = 'Selects which spectrum window is used to follow the root (amide) position of the query spin system'
    self.rootPanePulldown = PulldownList(frame, callback=self.changeRootPane,
                                         grid=(0,6), tipText=tipText)
    
    self.toleranceEntry = FloatEntry(self, returnCallback=self.setTolerance,
                                     text=0.1, width=6)
    
    tipTexts = ['The "experiment:spectrum" name of the spectrum for a potential query peak list; used as the source of peak positions for a spin system',
                'The serial number of a potential query peak list within its spectrum; used as the source of peak positions for a spin system',
                'Sets whether the peak list is currently used as a source of spin systems for peak-peak matching to find sequence neighbours',
                'Sets maximum ppm difference for peaks to be considered as potentially coming from the same resonance; in the spectrum dimension being matched',
                'The full CCPN experiment type for the peak list; should correspond to a partner experiment type in the match spectra']
    headingList = ['Spectrum', 'Peak List', 'Use?', 'Tolerance','Expt Type']
    editWidgets         = [None, None, None, self.toleranceEntry, None]
    editGetCallbacks    = [self.toggleQueryPeakList, None,
                           self.toggleQueryPeakList, self.getTolerance, None]
    editSetCallbacks    = [None, None, None, self.setTolerance, None]
    
    self.querySpecTable = ScrolledMatrix(frame, gridSpan=(1,7),
                                         highlightType=2, tipTexts=tipTexts,
                                         headingList=headingList,
                                         editWidgets=editWidgets,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         callback=self.selectPeakList,
                                         multiSelect=False, grid=(1,0))
    
    row +=1
    tipText = 'The spectra that matches to the query peaks are found within, and the windows to display them in'
    div = LabelDivider(frameA, text='Match Windows and Spectra',
                       grid=(row,0), tipText=tipText)
    
    row +=1
    frameA.expandGrid(row,0)
    frame = Frame(frameA, grid=(row,0))
    frame.expandGrid(1,4)
    
    label = Label(frame, text='13C Window:', grid=(0,0))
    tipText = 'Selects which spectrum window is used to display potential (sequence neighbour) spin system matches to 13C query locations'
    self.matchPanePulldownC = PulldownList(frame, callback=self.changeMatchPaneC,
                                           grid=(0,1), tipText=tipText)

    label = Label(frame, text='1H Window:', grid=(0,2))
    tipText = 'Selects which spectrum window is used to display potential (sequence neighbour) spin system matches to 1H query locations'
    self.matchPanePulldownH = PulldownList(frame, callback=self.changeMatchPaneH,
                                           grid=(0,3), tipText=tipText)
                                                    
    tipTexts = ['The "experiment:spectrum" name of the spectrum that may be searched for potential matches to the peak positions from the query spin systems',
                'The serial number of potential match peak list within its spectrum; may be searched for matches to query spin system peak positions',
                'Sets whether the peak list is currently used as a target for peak-peak matching to find sequence neighbours',
                'Sets maximum ppm difference for peaks to be considered as potentially coming from the same resonance; in the spectrum dimension being matched. The value also weights peak closeness scores',
                'The full CCPN experiment type for the peak list; should correspond to a partner experiment type in the query spectra']
    headingList = ['Spectrum', 'Peak List', 'Use?', 'Tolerance','Expt Type']
    editWidgets         = [None, None, None, self.toleranceEntry, None]
    editGetCallbacks    = [self.toggleMatchPeakList, None,
                           self.toggleMatchPeakList, self.getTolerance, None]
    editSetCallbacks    = [None, None, None, self.setTolerance, None]
    self.matchSpecTable = ScrolledMatrix(frame, gridSpan=(1,5),
                                         highlightType=2, tipTexts=tipTexts,
                                         headingList=headingList,
                                         editWidgets=editWidgets,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         callback=self.selectPeakList,
                                         multiSelect=False, grid=(1,0))
    

    #
    # Spin Systems
    #
    
    frameB.expandGrid(None, 0)
  
    row = 0
    
    #div = LabelDivider(frameB, text='Spin Systems', grid=(row,0))
    #row +=1

  
    frame = Frame(frameB, grid=(row,0), expandGrid=(0,0))
    
    row +=1
    frameB.grid_rowconfigure(row, weight=1)
    tipTexts = ['The serial number of the query spin system; used as a source of peak positions to match against',
                'The assignment annotation name for the query spin system',
                'The spin system that is currently linked as being one position later in the sequence, relative to the spin system for the row',
                'The spin system that is currently linked as being one position earlier in the sequence, relative to the spin system for the row',
                'The identity of any extra spin system that is considered as being at the same sequential position; could be an NH2 side chain group etc.',
                'For sequentially connected runs of spin systems, the number of the run and the position within the run at which the spin system is found',
                'The "root" location of the spin system (e.g. amide); the matching to find sequential matches used non-root peak positions']
    headingList = ['#','Name','i-1','i+1','i+0',
                   'Seq, Segment','Position']
    self.spinSystemMatrix = ScrolledMatrix(frameB, headingList=headingList,
                                           multiSelect=True, grid=(row,0),
                                           callback=self.selectSpinSystem,
                                           tipTexts=tipTexts)

    row +=1
    tipTexts = ['Show peak matches locations for the next spin system in the table, assuming a spin system is currently selected',
                'Show peak matches locations for the previous spin system in the table, assuming a spin system is currently selected',
                'Show peak matches locations for the spin system connected as one earlier in the sequence, relative to the currently selected one',
                'Show peak matches locations for the spin system connected as one later in the sequence, relative to the currently selected one']
    texts    = ['Next Row','Previous Row',
                'Goto i-1','Goto i+1']
    commands = [self.nextSpinSystem,self.prevSpinSystem,
                self.prevLinkedSpinSystem,self.nextLinkedSpinSystem]
    self.spinSystemButtons = ButtonList(frameB, texts=texts, tipTexts=tipTexts,
                                         commands=commands, grid=(row,0))
    
    row +=1
    tipTexts = ['For the selected spin system, clear any connections to spin systems set as being one earlier in the sequence',
                'For the selected spin system, clear any connections to spin systems set as being one later in the sequence',
                'For the selected spin system, clear all connections to any sequentially linked spin systems']
    texts    = ['Clear i-1 Links','Clear i+1 Links',
                'Clear All Links']
    commands = [self.clearPrevLinks, self.clearNextLinks,
                self.clearAllLinks,]
    self.spinSystemButtons2 = ButtonList(frameB, texts=texts, tipTexts=tipTexts,
                                        commands=commands, grid=(row,0))
    row +=1
    tipTexts = ['Considering how peaks in different spectra overlap, add non-root (HA/CO/CA/CB) resonance assignments to peaks in the selected spin system',
                'Show a table of peaks currently assigned to the selected spin system',
                'Show a table of the resonances currently contained within the selected spin system']
    texts    = ['Add HA/CO/CA/CB','Show Peaks',
                'Predict Type']
    commands = [self.addResonances, self.showPeaks,
                self.showResonances]
    self.spinSystemButtons3 = ButtonList(frameB, texts=texts, tipTexts=tipTexts,
                                         commands=commands, grid=(row,0))

    row +=1
    tipText = 'The locations of the amide spin systems that are potential sequence neighbours to the query spin system'
    div = LabelDivider(frameB, text='Matched Peak Positions',
                       grid=(row,0), tipText=tipText)

    row +=1
    frameB.grid_rowconfigure(row, weight=1)
    tipTexts = ['The rank of the spin system with regards to how well its peaks match the peaks in the query',
                'The average match score for the peaks matched in the spin system',
                'The total of the match score for all peaks matched in the spin system',
                'The identity (assignment annotation) for the spin system that was matched by considering peaks in the selected spectra',
                'The identity of any spin system connected to the matched spin system one earlier in the sequence',
                'The identity of any spin system connected to the matched spin system one later in the sequence',
                'The number of peaks involved in the position comparison between the query and the matched spin system']
    headingList = ['Rank','Mean Peak\nScore','Total Peak\nScore',
                   'Spin\nSystem','i-1','i+1','Num\nPeaks']
    self.matchMatrix = ScrolledMatrix(frameB, headingList=headingList,
                                      multiSelect=True, grid=(row,0),
                                      callback=self.selectMatch,
                                      tipTexts=tipTexts)

    
    row +=1
    tipTexts = ['Manually force a refresh of the above table',
                'Remove the currently selected match spin system from the table; only removes from display does not change anything else',
                'Set the selected spin system as being sequentially related to the query; sets up the appropriate sequence link and any inherited assignment',
                'Clear all ruler lines, that are use to mark peak positions, from the spectrum windows']
    texts = ['Refresh Matches', 'Discard Selected',
             'Set Seq Link', 'Clear Rulers']
    commands = [self.findMatches, self.discardMatches,
                self.setSeqLink ,self.clearRulers]
    self.matchButtons = ButtonList(frameB, texts=texts, tipTexts=tipTexts,
                                   commands=commands, grid=(row,0))

    row +=1
    frame = Frame(frameB, grid=(row,0))
    frame.expandGrid(0,0)
    frame.expandGrid(0,1)
    
    tipText = 'A ranked list of the sequence locations that best match the residue type possibilities of the query spin system and any sequential links'
    divL = LabelDivider(frame, text='Sequence Locations',
                        grid=(0,0), tipText=tipText)
    
    tipTexts = ['The ranking how well the query and selected match spin system fits a region of protein sequence',
                'The score for how well the query and selected match spin system, considering any set sequence links, fits the protein sequence',
                'For a potential sequence location, the residue two positions before the residue paired with the query spin system',
                'For a potential sequence location, the residue one position before the residue paired with the query spin system',
                'For a potential sequence location, the residue paired with the query spin system',
                'For a potential sequence location, the residue two positions after the residue paired with the query spin system',
                'For a potential sequence location, the residue two positions after the residue paired with the query spin system']
    headingList = ['Rank','Score','i-2','i-1','i','i+1','i+2']
    self.seqMatrix = ScrolledMatrix(frame, headingList=headingList,
                                    multiSelect=True, grid=(1,0),
                                    initialRows=4, tipTexts=tipTexts,
                                    callback=self.selectResidues)
    
    tipTexts = ['Assign the query spin system and any sequentially linked spin systems to the residue sequence at the selected location',
                'Remove any spin system assignments to the selected residue (position i), including any sequentially linked neighbours',
                'Display the spin systems in the selected sequence region as strips (within the match window)']
    texts = ['Assign','Deassign','Strip']
    commands = [self.assignSeq, self.deassignSeq, self.stripSeq]
    self.assignButtons = ButtonList(frame, texts=texts, commands=commands,
                                    grid=(2,0), tipTexts=tipTexts)
    
    tipText = 'The predicted residue types for the selected spin system and any sequence neighbours'
    divR = LabelDivider(frame, text='Residue Types',
                        grid=(0,1), tipText=tipText)
    
    tipTexts = ['The score for predicting the residue type of the i-2 spin system; connected as two sequence positions earlier than the query spin system',
                'The score for predicting the residue type of the i-1 spin system; connected as one sequence position earlier than the query spin system',
                'The score for predicting the residue type of the query spin system',
                'The score for predicting the residue type of the i+1 spin system; connected as one sequence position later than the query spin system',
                'The score for predicting the residue type of the i+2 spin system; connected as two sequence positions later than the query spin system']
    headingList = ['i-2','i-1','i','i+1','i+2']
    self.typeMatrix = ScrolledMatrix(frame, headingList=headingList,
                                     multiSelect=True, grid=(1,1),
                                     initialRows=5, tipTexts=tipTexts,
                                     gridSpan=(2,1), callback=None)
    
    #
    # Options
    #
    
    frameC.expandGrid(4,0)
    
    row = 0
    div = LabelDivider(frameC, text='Sequence', grid=(row,0))
    
    row +=1
    frame = Frame(frameC, grid=(row,0), expandGrid=(2,2))
    
    label = Label(frame, text='Chain:', grid=(0,0))
    tipText = 'Selects which molecular chain the assignment is done for; sets which sequence spin systems match and are assigned to'
    self.chainPulldown = PulldownList(frame, callback=self.changeChain,
                                      grid=(0,1), tipText=tipText)

    row += 1
    div = LabelDivider(frameC, text='Peak Matching', grid=(row,0))
    
    row +=1
    frame = Frame(frameC, grid=(row,0))
    frame.expandGrid(None, 4)
    
    label = Label(frame, text='Max Strips:', grid=(0,0))
    tipText = 'Sets the maximum number of matching spin systems to display in strips of the selected spectrum window'
    self.maxMatchesEntry = IntEntry(frame, width=2, text=7,
                                    grid=(0,1), tipText=tipText)

    label = Label(frame, text='13C Focus width (ppm):', grid=(1,0), tipText=tipText)
    tipText = 'When the "Focus 13C matches?" option is on, specifies how closely to zoom in on the 13C peak locations; the width of the window region'
    self.focusEntryC = FloatEntry(frame, width=8, text=DEFAULT_FOCUS_WIDTH,
                                  grid=(1,1), tipText=tipText)

    label = Label(frame, text='1H Focus width (ppm):', grid=(2,0))
    tipText = 'When the "Focus 1HC matches?" option is on, specifies how closely to zoom in on the 1H peak locations; the width of the window region'
    self.focusEntryH = FloatEntry(frame, width=8, text=DEFAULT_FOCUS_WIDTH_H,
                                  grid=(2,1), tipText=tipText)

    label = Label(frame, text='Auto Match', grid=(3,0))
    tipText = 'Whether to automatically get sequential spin system match candidates when selecting in the list of query spin systems'
    self.matchSelect = CheckButton(frame, callback=None, grid=(3,1), selected=True, tipText=tipText)

    label = Label(frame, text='Filter 13C By Inter/Intra Type', grid=(4,0))
    tipText = 'Whether to exclude peak matches for 13C positions where the i-1 resonance matches another i-1 resonance; gives more spurious matches but useful if i & i-1 peaks have the same shift'
    self.filterSelectC = CheckButton(frame, callback=self.findMatches,
                                     grid=(4,1), selected=True, tipText=tipText)

    label = Label(frame, text='Filter 1H By Inter/Intra Type', grid=(5,0))
    tipText = 'Whether to exclude peak matches for 1H positions where the i-1 resonance matches another i-1 resonance; gives more spurious matches but useful if i & i-1 peaks have the same shift'
    self.filterSelectH = CheckButton(frame, callback=self.findMatches,
                                     grid=(5,1), selected=True, tipText=tipText)
                                    
    label = Label(frame, text='Focus 13C matches?', grid=(6,0))
    tipText = 'Whether to zoom in on the 13C peak locations using the above stated focus width for the window region'
    self.focusSelectC = CheckButton(frame, callback=self.findMatches,
                                    grid=(6,1), selected=True, tipText=tipText)

    label = Label(frame, text='Focus 1H matches?', grid=(7,0))
    tipText = 'Whether to zoom in on the 1H peak locations using the above stated focus width for the window region'
    self.focusSelectH = CheckButton(frame, callback=self.findMatches,
                                    grid=(7,1), selected=True, tipText=tipText)

    label = Label(frame, text='Match Order:', grid=(8,0))
    index = MATCH_ORDERS.index(self.matchOrder)
    tipText = 'Sets whether potential sequential spin system matches should be ranked by the average peak position difference or total difference'
    self.matchOrderPulldown = PulldownList(frame, texts=MATCH_ORDERS, grid=(8,1), tipText=tipText,
                                           index=index, callback=self.changeMatchOrder)

    
    #
    # Main
    #
    
    self.bottomButtons = UtilityButtonList(tabbedFrame.sideFrame,
                                           helpUrl=self.help_url,
                                           grid=(0,0), sticky='e')

    self.administerNotifiers(self.registerNotify)
  
  def getRefExpNames(self):
    """ Added SeqNoesyExperiments to allow e.g. HSQC-TOCSY v. HSQC-NOESY
    """
  
    refExps, refExpsCO = getSeqAssignRefExperiments(self.project)
    #expNames = set([e.name for e in refExps])
    expNames = set([e.name for e in refExps] + SeqNoesyExperiments)
    expNamesCO = set([e.name for e in refExpsCO])
    
    return expNames, expNamesCO
    
  def checkExpTypeSanity(self):

    queryTypes = [pl.dataSource.experiment.refExperiment.name for pl in self.queryPeakLists]
    matchTypes = [pl.dataSource.experiment.refExperiment.name for pl in self.matchPeakLists]

    queryInter = [expType in self.refNamesCo for expType in queryTypes]
    matchInter = [expType in self.refNamesCo for expType in matchTypes]
    
    warn = ''
    
    if (True in queryInter) and (False in queryInter):
      warn += 'Query peak lists represent a mixture of i-1 specific'
      warn += ' and non-specific experiments\n'

    if (True in matchInter) and (False in matchInter):
      warn += 'Match peak lists represent a mixture of i-1 specific'
      warn += ' and non-specific experiments\n'
  
    if (not warn) and (True in matchInter) and (True in queryInter):
      warn += 'Both query and match spectra are i-1 specific;'
      warn += ' cannot set sequential links\n'
    
    if (not warn) and (False in matchInter) and (False in queryInter):
      warn += 'Neither query nor match spectra are i-1 specific;'
      warn += ' cannot set sequential links\n'
      
    if warn:
      showWarning('Experiment Type Error', warn, parent=self)
      return False
    else:
      return True  
      
  
  def changeMatchOrder(self, mode):
  
    if mode != self.matchOrder:
      self.matchOrder = mode
      self.findMatches()
      
  def getPeakListAnalysisDataDim(self, peakList):

    spectrum = peakList.dataSource
    
    for dataDim in spectrum.dataDims:
      transfers = []
      for expDimRef in dataDim.expDim.expDimRefs:
        transfers.extend(expDimRef.expTransfers)
      
      for expTransfer in transfers:
        if expTransfer.transferType in ('onebond','CP'):
          expDimRef1, expDimRef2 = expTransfer.expDimRefs
          
          if ('1H' in expDimRef1.isotopeCodes) and ('15N' in expDimRef2.isotopeCodes):
            # Its amide
            break
          
          elif ('1H' in expDimRef2.isotopeCodes) and ('15N' in expDimRef1.isotopeCodes):
            # Its amide
            break
    
      else:
        return getAnalysisDataDim(dataDim)

  def getTolerance(self, peakList):

    analysisDataDim = self.getPeakListAnalysisDataDim(peakList)

    self.toleranceEntry.set(analysisDataDim.assignTolerance)

  def setTolerance(self, *extra):
    
    analysisDataDim = self.getPeakListAnalysisDataDim(self.peakList)

    value = self.toleranceEntry.get()
    
    if value is None:
      value = analysisDataDim.assignTolerance
      
    analysisDataDim.assignTolerance = value 
 
  def selectTab(self, index):
  
    if index > 0:
      if not self.checkExpTypeSanity():
        self.tabbedFrame.select(0)
    
    if index == 1:
      self.updateAfter()

 
  def selectPeakList(self, obj, row, col):
  
    self.peakList = obj
  
  def addResonances(self):
    
    peakLists = self.queryPeakLists | self.matchPeakLists

    if self.spinSystem and peakLists:
      spinSystem = self.spinSystem;

      linkSpinSystemInterIntraResonances(spinSystem, peakLists) 
  
  def selectResidues(self, residues, row, col):
  
    self.residues = residues
  
  def assignSeq(self):
  
    if self.residues and self.spinSystem:
      spinSystem = self.spinSystem
      residue = self.residues[2]
      data = (spinSystem.serial, residue.seqCode, residue.ccpCode)
      msg = 'Really assign spin system %d to residue %d %s?' % data 
    
      if showOkCancel('Confirm', msg, parent=self):
        
        i = -1
        residueN = getLinkedResidue(residue, linkCode='prev')
        connSpinSystemN = findConnectedSpinSystem(spinSystem, delta=-1)
        while connSpinSystemN:
          
          if not residueN:
            msg = 'Linked spin systems would extend beyond the start of '
            msg += 'the sequence. Really continue with assignment?'
            
            if showOkCancel('Warning', msg, parent=self):
              break
            else:
              return  

          if residueN.ccpCode in PROLYL:
            msg = 'Residue at position i%d is a proline. ' % i
            msg += 'Really continue with assignment?'
            
            if not showOkCancel('Warning', msg, parent=self):
              return  
            
          i -= 1
          residueN = getLinkedResidue(residueN, linkCode='prev')
          connSpinSystemN = findConnectedSpinSystem(connSpinSystemN, delta=-1)
      
        residueC = getLinkedResidue(residue, linkCode='next')
        connSpinSystemC = findConnectedSpinSystem(spinSystem, delta=1)
        i = 1
        while connSpinSystemC:
          
          if not residueC:
            msg = 'Linked spin systems would extend beyond the end of '
            msg += 'the sequence. Really continue with assignment?'
            
            if showOkCancel('Warning', msg, parent=self):
              break
            else:
              return  

          if residueC.ccpCode in PROLYL:
            msg = 'Residue at position i+%d is a proline. ' % i
            msg += 'Really continue with assignment?'
            
            if not showOkCancel('Warning', msg, parent=self):
              return  

          i += 1
          residueC = getLinkedResidue(residueC, linkCode='next')
          connSpinSystemC = findConnectedSpinSystem(connSpinSystemC, delta=1)
      
      
        assignSpinSystemResidue(spinSystem, residue, warnMerge=True)
      
  def deassignSeq(self):
  
    if self.residues and self.spinSystem:
     
      residue = self.residues[2]
      data = (residue.seqCode, residue.ccpCode)
      msg = 'Really deassign spin system from residue %d %s and ' % data
      msg += 'deassign any sequentially connected connected spin systems?' 
      if showOkCancel('Confirm', msg, parent=self):
        for spinSystem in residue.resonanceGroups:
          deassignSpinSystem(spinSystem)
 
          connSpinSystemN = findConnectedSpinSystem(spinSystem, delta=-1)
          while connSpinSystemN:
            deassignSpinSystem(connSpinSystemN)
            connSpinSystemN = findConnectedSpinSystem(connSpinSystemN, delta=-1)

          connSpinSystemC = findConnectedSpinSystem(spinSystem, delta=1)
          while connSpinSystemC:
            deassignSpinSystem(connSpinSystemC)
            connSpinSystemC = findConnectedSpinSystem(connSpinSystemC, delta=1)
 
  def stripSeq(self):
  
    if self.residues and self.spinSystem:
      objectList = self.spinSystemMatrix.objectList
      querySpinSystems = [s[0] for s in objectList]
      querySpinSystemsSet = set(querySpinSystems)
      nullH = {}
      nullC = {}
 
      isotopesC = []
      if self.matchPaneC:
        nullC = {self.matchPaneC.spectrumWindow.stripAxis: 15, 'z1':99}
        stripAxisC = self.matchPaneC.spectrumWindow.stripAxis
        axisLabelsC = [stripAxisC,'z1','z2','z3']
        for label in axisLabelsC:
          axisPanel = self.matchPaneC.findFirstAxisPanel(label=label)
 
          if axisPanel:
            isotopesC.append(set(axisPanel.axisType.isotopeCodes))
          else:
            break
      
      isotopesH = []
      if self.matchPaneH:
        nullH = {self.matchPaneH.spectrumWindow.stripAxis: 15, 'z1':99}
        stripAxisH = self.matchPaneH.spectrumWindow.stripAxis
        axisLabelsH = [stripAxisH,'z1','z2','z3']
        for label in axisLabelsH:
          axisPanel = self.matchPaneH.findFirstAxisPanel(label=label)
 
          if axisPanel:
            isotopesH.append(set(axisPanel.axisType.isotopeCodes))
          else:
            break
            
       
      positionsH = []
      positionsC = []
      for residue in self.residues:
        if not residue:
          positionsC.append(nullC)
          positionsH.append(nullH)
          continue
        
        
        for spinSystem in residue.resonanceGroups:
          
          if spinSystem in querySpinSystemsSet:
            index = querySpinSystems.index(spinSystem)
            shifts = objectList[index][1]
            ppmIso = [(s.value, s.resonance.isotopeCode) for s in shifts]
              
            posC = {}
            posH = {}
            for i, isotopeC in enumerate(isotopesC):
              for ppm, isotope in ppmIso:
                if isotope in isotopeC:
                  posC[axisLabelsC[i]] = ppm
                  break

            for i, isotopeH in enumerate(isotopesH):
              for ppm, isotope in ppmIso:
                if isotope in isotopeH:
                  posC[axisLabelsH[i]] = ppm
                  break
 
            positionsC.append(posC)
            positionsH.append(posH)
             
          else: # E.g. Pro
            positionsC.append(nullC)
            positionsH.append(nullH)              
      
      if self.matchPaneC and positionsC:
        displayStrips(self.guiParent, positionsC, windowPane=self.matchPaneC)
      
      if self.matchPaneH and positionsH:
        displayStrips(self.guiParent, positionsH, windowPane=self.matchPaneH)
 
      
        
  def getSpinSystemScore(self, spinSystem, shifts, chain):
  
    scores = getShiftsChainProbabilities(shifts, chain)
    total = sum(scores.values()) 
    
    if total:
      for ccpCode in scores:
        scores[ccpCode] *= 100.0/total
    
    else:
      return scores
    
    return scores
    
  def updateResTypesAfter(self, resonance=None):

    if self.waiting:
      return 
    
    self.after_idle(self.updateResTypes)

  def updateResTypes(self):
  
    if not (self.spinSystem and self.queryPeakLists and self.chain):
      self.updateSeqLocations()
      return
    
    shiftList = list(self.queryPeakLists)[0].dataSource.experiment.shiftList
    spinSystem = self.spinSystem
    chain = self.chain
    getSpinSystemScore = self.getSpinSystemScore
    
    prev = None 
    next = None
    prev2 = None
    next2 = None

    if self.match:
      spinSystemM = self.match[2]
      experiment = self.match[4][0].peakList.dataSource.experiment
      refExperiment = experiment.refExperiment
      if refExperiment.name in self.refNamesCo:
        # Q i -> M i+1
        next = spinSystemM
      else:
        # Q i -> M i-1
        prev = spinSystemM
        
    if not prev:
      prev = findConnectedSpinSystem(spinSystem, delta=-1)
    if not next:
      next = findConnectedSpinSystem(spinSystem, delta=1)
    
    if prev:
      prev2 = findConnectedSpinSystem(prev, delta=-1)
    if next:
      next2 = findConnectedSpinSystem(next, delta=1)
    
    spinSystems = [prev2, prev, spinSystem, next, next2]
    scoreMatrix = [] 

    ccpCodes  = self.getCcpCodes(chain)
    N = len(ccpCodes)
    baseLevel = 100.0/N
    
    for spinSystem0 in spinSystems:
      scoreList = [None] * N
      
      if spinSystem0:
        shifts = []
        for resonance in spinSystem0.resonances:
          shift = resonance.findFirstShift(parentList=shiftList)
          if shift:
            shifts.append(shift)
 
        scores = getSpinSystemScore(spinSystem0, shifts, chain)
        for i, ccpCode in enumerate(ccpCodes):
          scoreList[i] = (scores[ccpCode], ccpCode)
 
        scoreList.sort()
        scoreList.reverse()

      scoreMatrix.append(scoreList)
  
    textMatrix = []
    colorMatrix = []
    objectList = []
    for j in range(N):
    
      datum = [None] * 5
      colors = [None] * 5
      for i, spinSystem0 in enumerate(spinSystems):
        scores = scoreMatrix[i][j]
        
        if scores:
          score, ccpCode = scores
          
          if score > 1.0:
            datum[i] = '%s %.1f' % (ccpCode, score)
 
            if score >= min(100.0,5*baseLevel):
              color = '#80ff80'
            elif score > 2*baseLevel:
              color = '#ffff80'
            elif score > baseLevel:
              color = '#ffc080'
            else:
              color = '#ff8080'
 
            colors[i] = color
 
      if scores == [None]*5:
        break

      objectList.append(j)
      textMatrix.append(datum)
      colorMatrix.append(colors)    
    
    self.typeMatrix.update(textMatrix=textMatrix,
                           colorMatrix=colorMatrix,
                           objectList=objectList)
                           
    self.updateSeqLocations(scoreMatrix, spinSystems)
    
  def updateSeqLocations(self, scoreMatrix=None, spinSystems=None):
  
    if self.waiting:
      return
  
    window = 5
  
    self.residues = None
  
    textMatrix = []
    colorMatrix = []
    objectList = []
    
    if self.chain and scoreMatrix:
      matches = []
    
      assignDict = {}
      for spinSystem in self.nmrProject.resonanceGroups:
        residue = spinSystem.residue
        if residue:
          assignDict[residue] = spinSystem
      
      residues = self.chain.sortedResidues()
      seq = [r.ccpCode for r in residues]
      
      seq = [None, None] + seq + [None, None]
      residues = [None, None] + residues + [None, None]
      nRes = len(seq)
      
      if nRes >= window:
        scoreDicts = []
        ccpCodes  = self.getCcpCodes(self.chain)
 
        for scores in scoreMatrix:
          scoreDict = {}
          for ccpCode in ccpCodes:
            scoreDict[ccpCode] = None
 
          for data in scores:
            if data: 
              score, ccpCode = data
              scoreDict[ccpCode] = score
 
          scoreDicts.append(scoreDict)
        
        sumScore = 0.0
        for i in range(nRes-window):
          if seq[i+2] and (seq[i+2] in PROLYL):
            continue
        
          score = 1.0
 
          for j in range(window):
            ccpCode = seq[i+j]
            score0 = scoreDicts[j].get(ccpCode)
            
            if (ccpCode is None) and (spinSystems[j]):
              break
            elif score0:
              score *= score0
            elif score0 == 0.0:
              break
           
          else:
            #score /= float(window)
            matches.append( (score, residues[i:i+window]) )
            sumScore += score
 
        matches.sort()
        matches.reverse()
 
        for i, data in enumerate(matches[:10]):
          score, residues = data
          score /= sumScore
          datum = [i+1, 100.0*score]
          color = hexRepr(1-score, score, 0)
           
          colors = [color, color] 
          for residue in residues:
            if residue:
              datum.append('%d%s' % (residue.seqCode, residue.ccpCode))
              if assignDict.get(residue):
                colors.append('#8080FF')
              else:
                colors.append(None)
            
            else:
              datum.append(None)
              colors.append(None)
 
          textMatrix.append(datum)
          objectList.append(residues)
          colorMatrix.append(colors)
    
    self.seqMatrix.update(textMatrix=textMatrix,
                          colorMatrix=colorMatrix,
                          objectList=objectList)
  
  
  def changeChain(self, chain):
  
    if chain is not self.chain:
      self.chain = chain
      
      if self.spinSystem:
        self.updateResTypes()
  
  def updateChains(self, obj=None):

    names = []
    chains = []
    index = 0
    chain = self.chain
    
    for molSystem in self.project.sortedMolSystems():
      msCode = molSystem.code
      
      for chain2 in molSystem.sortedChains():
        for residue in chain2.residues:
          if residue.molResidue.molType == 'protein':
            chains.append(chain2)
            names.append('%s:%s' % (msCode, chain2.code))
            break
 

    if chains:
      if chain not in chains:
        chain = chains[0]
        index = 0
      else:
        index = chains.index(chain)  
    
    else:
      chain = None

    self.changeChain(chain)

    self.chainPulldown.setup(names, chains, index)

  def getCcpCodes(self, chain):
  
    codeDict = {}
    for residue in chain.residues:
      codeDict[residue.ccpCode] = True

    ccpCodes = codeDict.keys()
    ccpCodes.sort()
    
    return ccpCodes
  
  def administerNotifiers(self, notifyFunc):

    for func in ('__init__', 'delete','setName'):
      for clazz in ('ccp.nmr.Nmr.DataSource', 'ccp.nmr.Nmr.Experiment'):
        notifyFunc(self.updateSpectra, clazz, func)
      notifyFunc(self.updateWindowsAfter, 'ccpnmr.Analysis.SpectrumWindow', func)
    notifyFunc(self.updateWindowsAfter, 'ccpnmr.Analysis.SpectrumWindow', 'setStripAxis')
    notifyFunc(self.updateSpectra, 'ccp.nmr.Nmr.Experiment', 'setRefExperiment')
   
    notifyFunc(self.updateSpectra, 'ccpnmr.Analysis.AnalysisDataDim', 'setAssignTolerance')
      
    notifyFunc(self.updateWindowsAfter, 'ccpnmr.Analysis.SpectrumWindowView', 'setIsPosVisible')
    notifyFunc(self.updateWindowsAfter, 'ccpnmr.Analysis.SpectrumWindowView', 'setIsNegVisible')

    for func in ('__init__', 'delete', 'setCcpCode', 'setResidue', 'setName',
                 'addResonance', 'removeResonance', 'setResonances'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)

    for func in ('__init__','delete'):
      # peakList needs at least one peak to appear in tables
      # so need to do notifier on Peak, not PeakList
      notifyFunc(self.updateSpectra, 'ccp.nmr.Nmr.Peak', func)

    for func in ('__init__', 'delete'):
      notifyFunc(self.updatePeakDimContribAfter, 'ccp.nmr.Nmr.PeakDimContrib', func)

    for func in ('setPosition',):
      notifyFunc(self.updatePeakDimAfter, 'ccp.nmr.Nmr.PeakDim', func)

    for func in ('setValue',):
      notifyFunc(self.updatePeakIntensityAfter, 'ccp.nmr.Nmr.PeakIntensity', func)

    for func in ('__init__', 'delete', 'setIsSelected', 'setSequenceOffset'):
      notifyFunc(self.updateLinks, 'ccp.nmr.Nmr.ResonanceGroupProb', func)

    for func in ('__init__', 'delete',):
      notifyFunc(self.updateChains, 'ccp.molecule.MolSystem.Chain', func)

    for func in ('setAssignNames','addAssignName','removeAssignName'):
      notifyFunc(self.updateResTypesAfter, 'ccp.nmr.Nmr.Resonance', func)

  def getAppDataOptions(self):
    
    project         = self.project
    analysisProject = self.analysisProject
    nmrProject      = self.nmrProject
  
    data = self.analysisProject.linkSeqSpinSystemsData
    if data:
      try:
        options = cPickle.loads(data)
      except:
        options = {} # Mal formed opt string - prob XML truncated
    
    else:
      options = {}

    for attrib in ('query%sC','match%sC','query%sH','match%sH'):
      serial1, serial2 = options.get(attrib % 'Pane', (None, None))
      serial = options.get(attrib[:-1] % 'Window', None)
            
      if serial1:
        window = analysisProject.findFirstSpectrumWindow(serial=serial1)
     
        if window:
          windowPane = window.findFirstSpectrumWindowPane(serial=serial2)
          setattr(self, attrib % 'Pane', windowPane)
          
      elif serial:
        window = analysisProject.findFirstSpectrumWindow(serial=serial)
        
        if window:
          windowPane = window.findFirstSpectrumWindowPane()
          setattr(self, attrib % 'Pane', windowPane)
          
      else:
        windowPane = getattr(self, attrib % 'Pane')
    
    for attrib in ('query','match',):
      peakLists = set()
      toggleDict = getattr(self, attrib + 'ToggleDict')

      for serials in options.get(attrib +'PeakLists', []):
      
        s1,s2,s3 = serials
        experiment = nmrProject.findFirstExperiment(serial=s1)
        if not experiment:
          continue

        spectrum = experiment.findFirstDataSource(serial=s2)
        if not spectrum:
          continue
 
        analysisSpectrum = spectrum.analysisSpectrum
        if not analysisSpectrum:
          continue
 
        getView = analysisSpectrum.findFirstSpectrumWindowView
        for elem in ('H','C'):
          windowPane = getattr(self, attrib + 'Pane' + elem)
         
          if getView(spectrumWindowPane=windowPane):
            peakList = spectrum.findFirstPeakList(serial=s3)
            if peakList:
              peakLists.add(peakList)
              toggleDict[peakList] = True
            
            break
 
      setattr(self, attrib + 'PeakLists', peakLists)

    self.matchSelect.set( options.get('matchSelect',True) )
    self.focusSelectH.set( options.get('focusSelectH',True) )
    self.focusSelectC.set( options.get('focusSelectC',True) )
    self.focusEntryC.set( options.get('focusEntry',DEFAULT_FOCUS_WIDTH) )
    self.focusEntryH.set( options.get('focusEntryH',DEFAULT_FOCUS_WIDTH_H) )
    self.maxMatchesEntry.set( options.get('maxMatches',7) )
    
  def getPeakListId(self, peakList):
  
    spectrum = peakList.dataSource

    return (spectrum.experiment.serial,spectrum.serial,peakList.serial)
    
  def setAppDataOptions(self):
    
    project = self.project
 
    options = {}
    
    if self.queryPaneC:  
      options['queryPaneC'] = (self.queryPaneC.spectrumWindow.serial,
                               self.queryPaneC.serial)
    if self.matchPaneC:  
      options['matchPaneC'] = (self.matchPaneC.spectrumWindow.serial,
                               self.matchPaneC.serial)
    if self.queryPaneH:  
      options['queryPaneH'] = (self.queryPaneH.spectrumWindow.serial,
                               self.queryPaneH.serial)
    if self.matchPaneH:  
      options['matchPaneH'] = (self.matchPaneH.spectrumWindow.serial,
                               self.matchPaneH.serial)

    ids = [self.getPeakListId(pl) for pl in self.queryPeakLists]
    ids.sort()
    options['queryPeakLists'] = ids 
    
    ids = [self.getPeakListId(pl) for pl in self.matchPeakLists]
    ids.sort()
    options['matchPeakLists'] = ids
    
    options['matchSelect'] = self.matchSelect.get()
    options['focusSelectH'] = self.focusSelectH.get()
    options['focusSelectC'] = self.focusSelectC.get()
    options['focusEntry'] = self.focusEntryC.get() or DEFAULT_FOCUS_WIDTH
    options['focusEntryH'] = self.focusEntryH.get() or DEFAULT_FOCUS_WIDTH_H
    options['maxMatches'] = self.maxMatchesEntry.get() or 7
    
    data = cPickle.dumps(options)
    self.analysisProject.linkSeqSpinSystemsData = data

  def updateLinks(self, resonanceGroupProb):
    
    if resonanceGroupProb.linkType == 'sequential':
      self.updateSeqSegmentsAfter()
    
    self.updateResTypesAfter()

  def clearRulers(self):
  
    for ruler in self.rulers:
      if not ruler.isDeleted:
        ruler.delete()
    
    self.rulers = []    

  def getRootWindows(self):
  
    windows = []
    project = self.project
    analysisProject = self.analysisProject
    
    if project:

      for window in analysisProject.spectrumWindows:
        windowPanes = window.sortedSpectrumWindowPanes()
      
        for windowPane in windowPanes:
          if len(windowPane.axisPanels) == 2:
            axisPanel = windowPane.findFirstAxisPanel(label='x')
            if axisPanel and axisPanel.axisType:
              xIsotopes = set(axisPanel.axisType.isotopeCodes)
            else:
              continue
            axisPanel = windowPane.findFirstAxisPanel(label='y')
            if axisPanel and axisPanel.axisType:
              yIsotopes = set(axisPanel.axisType.isotopeCodes)
            else:
              continue
 
            name = getWindowPaneName(windowPane)
            
            if ('15N' in xIsotopes) and ('1H' in yIsotopes):
              windows.append([name, windowPane])
              
            elif ('15N' in yIsotopes) and ('1H' in xIsotopes):
              windows.append([name, windowPane])
     
    windows.sort()
    
    return windows

  def getWindows(self, isMatch=False):
  
    windowsC = []
    windowsH = []
    project = self.project
    analysisProject = self.analysisProject
    
    if project:
      queryIsotopesC = []
      queryIsotopesH = []
      
      if isMatch and self.queryPaneC:
        for axisPanel in self.queryPaneC.sortedAxisPanels():
          queryIsotopesC.append(axisPanel.axisType.isotopeCodes)

      if isMatch and self.queryPaneH:
        for axisPanel in self.queryPaneH.sortedAxisPanels():
          queryIsotopesH.append(axisPanel.axisType.isotopeCodes)

      for window in analysisProject.spectrumWindows:
        windowPanes = window.sortedSpectrumWindowPanes()
        
        orthoAxis = X_Y[1-X_Y.index(window.stripAxis)]
      
        for windowPane in windowPanes:
          if len(windowPane.axisPanels) > 2:
            matchIsotopes = []
            stripIsotopes = []
            isotopes = []
 
            for axisPanel in windowPane.sortedAxisPanels():
              isotopeCodes = axisPanel.axisType.isotopeCodes
              isotopes.append(isotopeCodes)
 
              for isotope in isotopeCodes:
                if axisPanel.label == orthoAxis:
                  matchIsotopes.append(isotope)
                else:
                  stripIsotopes.append(isotope)
 
            name = getWindowPaneName(windowPane)
            
            if not queryIsotopesH or (isotopes == queryIsotopesH):
              if ('15N' in stripIsotopes) and ('1H' in matchIsotopes):
                windowsH.append([name, windowPane])
                        
            if not queryIsotopesC or (isotopes == queryIsotopesC):
              if ('15N' in stripIsotopes) and ('13C' in matchIsotopes):
                windowsC.append([name, windowPane])
     
    windowsC.sort()
    windowsH.sort()
    
    return windowsC, windowsH

  def changeQueryPaneC(self, windowPane):
  
    if windowPane is not self.queryPaneC:
      self.queryPaneC = windowPane
      self.updateQuerySpectra()
      self.updateMatchPanes()
      self.updateAfter()

  def changeQueryPaneH(self, windowPane):
  
    if windowPane is not self.queryPaneH:
      self.queryPaneH = windowPane
      self.updateQuerySpectra()
      self.updateMatchPanes()
      self.updateAfter()
  
  def updateQueryPanes(self):
  
    windowDataC, windowDataH = self.getWindows()

    index  = 0
    names  = []
    window = self.queryPaneC
    windows = []
 
    if windowDataC:
      names   = ['<None>',] + [x[0] for x in windowDataC]
      windows = [None,] + [x[1] for x in windowDataC]
      
      if window not in windows:
        window = windows[0]
        index = 0
 
      else:
        index = windows.index(window)
 
    else:
      window = None
      
    if window is not self.queryPaneC:
      self.queryPaneC = window
      if not self.matchPaneC:
        self.updateMatchPanes()
      self.updateQuerySpectra()
      self.updateAfter()
 
 
    self.queryPanePulldownC.setup(names, windows, index)
      
    index  = 0
    names  = []
    window = self.queryPaneH
    windows = []
 
    if windowDataH:
      names   = ['<None>',] + [x[0] for x in windowDataH]
      windows = [None,] + [x[1] for x in windowDataH]
      
      if window not in windows:
        window = windows[0]
        index  = 0
 
      else:
        index = windows.index(window)
 
    else:
      window = None
      
    if window is not self.queryPaneH:
      self.queryPaneH = window
      if not self.matchPaneH:
        self.updateMatchPanes()
      self.updateQuerySpectra()
      self.updateAfter()
 
    self.queryPanePulldownH.setup(names, windows, index)

    index  = 0
    names  = []
    window = self.rootPane
    windows = []
    
    windowData = self.getRootWindows()
    
    if windowData:
      names   = ['<None>',] + [x[0] for x in windowData]
      windows = [None,] + [x[1] for x in windowData]
      
      if window not in windows:
        window = windows[0]
        index  = 0
 
      else:
        index = windows.index(window)
 
    else:
      window = None
      
    if window is not self.rootPane:
      self.rootPane = window

    self.rootPanePulldown.setup(names, windows, index)

  def changeRootPane(self, window):
  
    if window is not self.rootPane:
      self.rootPane = window

  def changeMatchPaneH(self, window):
  
    if window is not self.matchPaneH:
      self.matchPaneH = window
      self.updateMatchSpectra()

  def changeMatchPaneC(self, window):
  
    if window is not self.matchPaneC:
      self.matchPaneC = window
      self.updateMatchSpectra()
  
  def updateMatchPanes(self):
    
    windowDataC, windowDataH = self.getWindows(isMatch=True)
   
    index = 0
    names = []
    windows = []
    window = self.matchPaneC
    
    if self.queryPaneC:
      windows = [x[1] for x in windowDataC]

    if windows:
      names = [x[0] for x in windowDataC]
      
      if window not in windows:
        window = windows[-1]
      
      index = windows.index(window)
 
    else:
      window = None
      
    if window is not self.matchPaneC:
      self.matchPaneC = window
      self.updateMatchSpectra()
 
    self.matchPanePulldownC.setup(names, windows, index)
    
    index = 0
    names = []
    windows = []
    window = self.matchPaneH
    
    if self.queryPaneH:
      windows = [x[1] for x in windowDataH]

    if windows:
      names = [x[0] for x in windowDataH]
      
      if window not in windows:
        window = windows[-1]
      
      index = windows.index(window)
 
    else:
      window = None
      
    if window is not self.matchPaneH:
      self.matchPaneH = window
      self.updateMatchSpectra()
 
    self.matchPanePulldownH.setup(names, windows, index)

  def toggleQueryPeakList(self, peakList):
      
    if peakList in self.queryPeakLists:
      self.queryPeakLists.remove(peakList)
      self.queryToggleDict[peakList] = False
    else:
      self.queryPeakLists.add(peakList)
      self.queryToggleDict[peakList] = True

    self.updateQuerySpectra()

  def toggleMatchPeakList(self, peakList):
   
    if peakList in self.matchPeakLists:
      self.matchPeakLists.remove(peakList)
      self.matchToggleDict[peakList] = False
    else:
      self.matchPeakLists.add(peakList)
      self.matchToggleDict[peakList] = True

    self.updateMatchSpectra()

  def getSpectra(self, windowPanes):
  
    spectra = set()
    
    for windowPane in windowPanes:
      if not windowPane:
        continue
        
      for view in windowPane.spectrumWindowViews:
        #if not (view.isPosVisible or view.isNegVisible):
        #  continue
        
        spectrum = view.analysisSpectrum.dataSource
 
        peakLists = spectrum.peakLists
        if peakLists:
          refExperiment = spectrum.experiment.refExperiment
          if refExperiment and refExperiment.name in self.refNames:
            spectra.add(spectrum)
          
          # remember prev selections 
          for peakList in peakLists:
            if not peakList.peaks:
              continue
          
            if self.queryToggleDict.get(peakList):
              if peakList not in self.queryPeakLists:
                self.queryPeakLists.add(peakList)
            if self.matchToggleDict.get(peakList):
              if peakList not in self.matchPeakLists:
                self.matchPeakLists.add(peakList)
    
    sortList = [('%s:%s' % (s.experiment.name,s.name), s) for s in spectra]
    sortList.sort()
          
    return [x[1] for x in sortList]

  def updateSpectra(self, obj=None):
    
    self.after_idle(self.updateQuerySpectra)
    self.after_idle(self.updateMatchSpectra)

  def updateQuerySpectra(self):
  
    spectra = self.getSpectra((self.queryPaneH, self.queryPaneC))
    textMatrix = []
    colorMatrix = []
    objectList = []
    
    for spectrum in spectra:
      hexColors = spectrum.analysisSpectrum.posColors
      hexColor = hexColors[int(0.7*len(hexColors))]
      experiment = spectrum.experiment
      refExperiment = experiment.refExperiment
      
      if refExperiment:
        expType = refExperiment.name
      else:
        expType = None
      
      specName = '%s:%s' % (experiment.name, spectrum.name)    
            
      for peakList in spectrum.sortedPeakLists():
        if not peakList.findFirstPeak(): # empty peakList
          continue
        
        if peakList in self.queryPeakLists:
          use = 'Yes'
          colors = [hexColor, None, None, None]
        else:
          use = 'No'
          colors = [None] * 4
        
        analysisDataDim = self.getPeakListAnalysisDataDim(peakList)
        
        datum = [specName, peakList.serial, use, analysisDataDim.assignTolerance, expType]
        
        colorMatrix.append(colors)
        textMatrix.append(datum)
        objectList.append(peakList)
    
    self.querySpecTable.update(textMatrix=textMatrix,
                               colorMatrix=colorMatrix,
                               objectList=objectList)

    for peakList in list(self.queryPeakLists):
      if peakList.dataSource not in spectra:
        self.queryPeakLists.remove(peakList)

    self.setAppDataOptions()

  def updateMatchSpectra(self):
  
    spectra = self.getSpectra((self.matchPaneH, self.matchPaneC))
    textMatrix = []
    colorMatrix = []
    objectList = []
    
    for spectrum in spectra:
      hexColors = spectrum.analysisSpectrum.posColors
      hexColor = hexColors[int(0.7*len(hexColors))]
      experiment = spectrum.experiment
      refExperiment = experiment.refExperiment
      
      if refExperiment:
        expType = refExperiment.name
      else:
        expType = None
      
      specName = '%s:%s' % (experiment.name, spectrum.name)    
          
      for peakList in spectrum.sortedPeakLists():
        if not peakList.findFirstPeak(): # empty peakList
          continue
        
        if peakList in self.matchPeakLists:
          use = 'Yes'
          colors = [hexColor, None, None, None]
        else:
          use = 'No'
          colors = [None] * 4
        
        analysisDataDim = self.getPeakListAnalysisDataDim(peakList)
        datum = [specName, peakList.serial, use, analysisDataDim.assignTolerance, expType]
        
        colorMatrix.append(colors)
        textMatrix.append(datum)
        objectList.append(peakList)
    
    self.matchSpecTable.update(textMatrix=textMatrix,
                               colorMatrix=colorMatrix,
                               objectList=objectList)

    for peakList in list(self.matchPeakLists):
      if peakList.dataSource not in spectra:
        self.matchPeakLists.remove(peakList)
 
    self.setAppDataOptions()


  def nextSpinSystem(self):
  
    if self.spinSystemMatrix.currentCell:
      row, col = self.spinSystemMatrix.currentCell
      index = row+1 % len(self.spinSystemMatrix.objectList)
      self.spinSystemMatrix.selectNthObject(index)
 
    
  def prevSpinSystem(self):
  
    if self.spinSystem:
      row, col = self.spinSystemMatrix.currentCell
      index = row-1 % len(self.spinSystemMatrix.objectList)
      self.spinSystemMatrix.selectNthObject(index)

  def selectMatch(self, obj, row, col):
      
    if obj:
      if obj != self.match:
        self.match = obj
        self.updateButtons()
        self.updateResTypesAfter()

  def getSpinSystemSeqSegment(self, spinSystem):
      
    for i, segment in enumerate(self.seqSegments):
      if spinSystem in segment:
        j = segment.index(spinSystem)
        return i, j
 
    self.seqSegments.append([spinSystem])
 
    return i, 0

  def updateSeqSegmentsAfter(self, *opt):
      
    if self.refreshSeqSeg:
      return
    else:
      self.refreshSeqSeg = True
      self.after_idle(self.updateSeqSegments)

  def updateSeqSegments(self):
        # called by notifier when links are made
    
    self.seqSegments = []
    if self.guiParent.project:
    
      seqSegments = []
      done = {}
      for spinSystem in self.nmrProject.resonanceGroups:   
        if not done.get(spinSystem):
          segment = [spinSystem,]
          done[spinSystem] = 1
          
          spinSystemN = findConnectedSpinSystem(spinSystem, delta=-1)
          while spinSystemN and not done.get(spinSystemN):
            segment.insert(0,spinSystemN)
            done[spinSystemN] = 1
            spinSystemN = findConnectedSpinSystem(spinSystemN, delta=-1)
                   
          spinSystemC = findConnectedSpinSystem(spinSystem, delta=1)
          while spinSystemC and not done.get(spinSystemC):
            segment.append(spinSystemC)
            done[spinSystemC] = 1
            spinSystemC = findConnectedSpinSystem(spinSystemC, delta=1)

          residue = segment[0].residue
          if residue:
            chain = residue.chain
            index = 'a%s%s%5.5d' % (chain.molSystem.code, chain.code, residue.seqCode)
          elif len(segment) > 1:
            index = 'b%d%d' % (9999-len(segment),segment[0].serial)
          else:
            index = 'c%d' % segment[0].serial
          
          seqSegments.append((index,segment))
  
      seqSegments.sort()
      self.seqSegments = [x[1] for x in seqSegments]
      self.updateAfter()
    

  def updateMatches(self):
      
    objectList = []
    textMatrix = []
    maxStrips = self.maxMatchesEntry.get() or 7
    
    windowPanes = ((self.matchPaneC, self.queryPaneC),
                   (self.matchPaneH, self.queryPaneH))
 
    getSpinSystemName = self.getSpinSystemName
    
    if not self.matches:
      self.setAppDataOptions()
      self.matchMatrix.update(textMatrix=textMatrix,objectList=objectList)
      return
    
    axesAndPanes = []
    for matchPane, queryPane in windowPanes:
      if not matchPane:
        continue

      stripAxis = matchPane.spectrumWindow.stripAxis
      axisLabels = [stripAxis,'z1','z2','z3']
      
      orthoAxis = 'y'
      if stripAxis == 'y':
        orthoAxis = 'x'
      axesAndPanes.append((matchPane, queryPane, orthoAxis, axisLabels, []))
      
    for j, match in enumerate(self.matches):
      score, scoreP, spinSystem, position, peaks = match

      if not spinSystem:
        continue

      spinSystemP = findConnectedSpinSystem(spinSystem, delta=-1)
      spinSystemN = findConnectedSpinSystem(spinSystem, delta=1)
 
      ssNameP = '?'
      ssNameN = '?'
      ssName = getSpinSystemName(spinSystem)
 
      if spinSystemP:
        ssNameP = getSpinSystemName(spinSystemP)
      if spinSystemN:
        ssNameN = getSpinSystemName(spinSystemN)
 
      datum = ['%d' % (j+1),
               '%.3f' % score,
               '%.3f' % scoreP,
               ssName  or None,
               ssNameP or None,
               ssNameN or None,
               '%d' % len(peaks)]
 
      textMatrix.append(datum)
      objectList.append(match)
 
      for matchPane, queryPane, orthoAxis, axisLabels, positions in axesAndPanes:
        if matchPane and j+1 <= maxStrips:
 
          pos = {}
          for i, ppm in enumerate(position):
            pos[axisLabels[i]] = ppm
 
          positions.append(pos)
    
    for matchPane, queryPane, orthoAxis, axisLabels, positions in axesAndPanes:
    
      if positions and not self.refreshSeqSeg:
        displayStrips(self.guiParent, positions,
                      spectrum=None, windowPane=matchPane)
                      
        # often fails first time...
        self.update_idletasks()
        displayStrips(self.guiParent, positions,
                      spectrum=None, windowPane=matchPane)
 
      if queryPane is self.queryPaneC:
        focus = self.focusSelectC.get()
        matchPane = self.matchPaneC
      else:
        focus = self.focusSelectH.get()
        matchPane = self.matchPaneH
 
      if focus:
        orthoPpms = self.getQueryOrthoPositions(queryPane, matchPane)
 
        if orthoPpms:
          self.defineWindowRegions(matchPane, orthoPpms, orthoAxis)
 

    self.setAppDataOptions()
    self.matchMatrix.update(textMatrix=textMatrix,objectList=objectList)
  
  def defineWindowRegions(self, windowPane, positions, axis, width=None):
  
  
    panel = windowPane.findFirstAxisPanel(label=axis)
    
    if width is None:
      if '1H' in panel.axisType.isotopeCodes:
        width = self.focusEntryH.get() or DEFAULT_FOCUS_WIDTH_H
      else:
        width = self.focusEntryC.get() or DEFAULT_FOCUS_WIDTH
      width /= 2.0
    
    windowFrame = windowPane.getWindowFrame()
    
    if axis == 'y':
      windowFrame.setRowsTo(len(positions))
    else:
      windowFrame.setColsTo(len(positions))
 
    regions = panel.sortedAxisRegions()
    for i, ppm in enumerate(positions):
      region = (ppm-width, ppm+width)
      regions[i].region = region
  
  
  def getQueryOrthoPositions(self, queryPane, matchPane):
     
    if self.spinSystem:
      getWindowView = queryPane.findFirstSpectrumWindowView
      getAnalysisSpectrum = self.analysisProject.findFirstAnalysisSpectrum
      
      orthoAxis = 'y'
      if queryPane.spectrumWindow.stripAxis == 'y':
        orthoAxis = 'x'
    
      peaks, dimDict, dimMappingDict = self.getMatchPeaks(self.shifts, queryPane, matchPane,
                                       self.queryPeakLists, self.matchPeakLists)
      orthoPos  = []
      tolerances = []
      for peak in peaks:
        analysisSpectrum = getAnalysisSpectrum(dataSource=peak.peakList.dataSource)
        
        if not analysisSpectrum:
          continue
          
        if not getWindowView(analysisSpectrum=analysisSpectrum):
          continue
      
        dimMapping = getPeakDimAxisMapping(peak, queryPane)
        peakDim = dimMapping.get(orthoAxis)
        analysisDataDim = peakDim.dataDim.analysisDataDim
        tolerances.append(analysisDataDim.assignTolerance)
        orthoPos.append(peakDim.value)
     
      if len(orthoPos) > 1:
        orthoPos.sort()
 
        while True:
           for i in range(len(orthoPos)-1):
             tolerance = tolerances[i]
             ppm1 = orthoPos[i]
             ppm2 = orthoPos[i+1]
             
             if abs(ppm1-ppm2) < tolerance:
               orthoPos[i] = (ppm1+ppm2)*0.5
               orthoPos.pop(i+1)
               break
 
           else:
             break
   
      return orthoPos
   
  def getShiftPeaks(self, shifts, peakLists):
      
    peakLists = set(peakLists)
    peaks = set()
    
    for shift in shifts:
      if not shift:
        continue
    
      resonance = shift.resonance
      for contrib in resonance.peakDimContribs:
        peak = contrib.peakDim.peak
        if peak.peakList in peakLists:
          if peak.findFirstPeakIntensity(intensityType='height'):
            peaks.add(peak)
 
    return list(peaks)

  def categorisePeakLists(self, peakLists):
  
    # Add HNN later
    
    peakListsH = []
    peakListsC = []
  
    for peakList in peakLists:
      isotopes = getSpectrumIsotopes(peakList.dataSource)
      
      if isotopes.count('1H') > 1:
        peakListsH.append(peakList)
      else:
        peakListsC.append(peakList)
  
    return peakListsH, peakListsC

  def getMatchPeaks(self, shifts, queryPane, matchPane, queryPeakLists, matchPeakLists):
  
    peaks = self.getShiftPeaks(shifts, queryPeakLists)
    isPeakPositive = self.isPeakPositive
    
    # All peaks at the mo - refine later?
    
    if not peaks:
      return [], {}, {}

    if matchPane.spectrumWindow.stripAxis == 'y':
      orthoAxisM = 'x'
    else:
      orthoAxisM = 'y'

    if queryPane.spectrumWindow.stripAxis == 'y':
      orthoAxisQ = 'x'
    else:
      orthoAxisQ = 'y'

    dimDict = {}
    dimMappingDict = {}
    for peakList in matchPeakLists:
      mPeak = peakList.findFirstPeak()
      dimMapping = getPeakDimAxisMapping(mPeak, matchPane)
      dim = dimMapping.get(orthoAxisM).dim
      dimDict[peakList] = dim
      tt = []
      nn = 0
      for label in sorted(dimMapping):
        if label != orthoAxisM:
          tt.append((dimMapping[label].dim, nn))
          nn += 1
      tt.sort()
      tt = [xx[1] for xx in tt]
      dimMappingDict[peakList] = tt
    
    for peakList in queryPeakLists:
      qPeak = peakList.findFirstPeak()
      dimMapping = getPeakDimAxisMapping(qPeak, queryPane)
      dim = dimMapping.get(orthoAxisQ).dim
      dimDict[peakList] = dim
     
    interResidueQuery = self.haveInterResiduePeakList(queryPeakLists)
    
    if not interResidueQuery: # E.g. query with HNCA
      peaks2 = self.getShiftPeaks(shifts, matchPeakLists)
      if len(peaks) > len(peaks2):
        for peak2 in peaks2:
          dimM = dimDict[peak2.peakList]
          ppm2 = peak2.findFirstPeakDim(dim=dimM).value
          posSign2 = isPeakPositive(peak2)
        
          overlap = []
          nonOverlap = 0
        
          for peak in peaks:
            dimQ = dimDict[peak.peakList]
            peakDim = peak.findFirstPeakDim(dim=dimQ)
            tolerance = getAnalysisDataDim(peakDim.dataDim).assignTolerance
            ppm = peakDim.value
            posSign = isPeakPositive(peak)
          
            if posSign == posSign2:
              if abs(ppm2-ppm) < tolerance:
                overlap.append(peak)
              else:
                nonOverlap += 1
         
          if nonOverlap and (len(overlap) == 1):
            peaks.remove(overlap[0])  
    
    # Setup data for ruler marks            

    panelPositions  = []
    for peak in peaks:
      dimMapping = getPeakDimAxisMapping(peak, queryPane)
      peakDim = dimMapping.get(orthoAxisQ)
    
      for axisPanel in queryPane.sortedAxisPanels():
        if dimMapping[axisPanel.label] is peakDim:
          panelPositions.append((peakDim.value, axisPanel.panelType))
          break

    # Make display lines

    for position, panelType in panelPositions:
      ruler = createRuler(position, panelType, dashLength=3,
                          gapLength=1, color=self.color, remove=False)
      self.rulers.append(ruler)
    
    return peaks, dimDict, dimMappingDict

  def haveInterResiduePeakList(self, peakLists):
  
    interResidueQuery = False
    for peakList in peakLists:
      if peakList.dataSource.experiment.refExperiment.name in self.refNamesCo:
        interResidueQuery = True
        break
    
    return  interResidueQuery

  def findMatches(self, event=None):
          
    self.match = None
    self.matches = []

    if not self.spinSystem:
      self.updateMatches()
      self.waiting = False
      return
    
    querySpinSystem = self.spinSystem
    queryPeakLists = self.queryPeakLists
    matchPeakLists = self.matchPeakLists
    
    if not queryPeakLists:
      showWarning('Warning','No query spectra selected', parent=self)
      self.waiting = False
      return

    if not matchPeakLists:
      #showWarning('Warning','No match spectra selected', parent=self)
      self.waiting = False
      return
    
    isPeakPositive = self.isPeakPositive
    queryPaneC = self.queryPaneC
    queryPaneH = self.queryPaneH
    matchPaneC = self.matchPaneC
    matchPaneH = self.matchPaneH
    
    if not (queryPaneC or queryPaneH):
      showWarning('Warning','No query windows selected', parent=self)      
      self.waiting = False
      return 
    
    if not self.haveInterResiduePeakList(list(queryPeakLists) + list(matchPeakLists)):
      msg = 'Neither query nor matches spectra represent '
      msg += 'any i->i-1 specific experiments'
      showWarning('Warning', msg , parent=self)
      return 
    
    popups = []
    for windowPane in (matchPaneC, matchPaneH):
      if windowPane:
        popup = windowPane.getWindowFrame().windowPopup
        popup.open()
        popups.append( (popup, popup.turnDrawRequestsOff()) )

    filterByTypeC = self.filterSelectC.get()
    filterByTypeH = self.filterSelectH.get()
    ssShifts = self.shifts
    
    # Split peak matching by H/C/N
    
    queryPeakListsH, queryPeakListsC = self.categorisePeakLists(queryPeakLists)
    matchPeakListsH, matchPeakListsC = self.categorisePeakLists(matchPeakLists)
    
    # Find Carbon matches
    
    
    data = ((queryPaneC, matchPaneC, queryPeakListsC, matchPeakListsC, filterByTypeC),
            (queryPaneH, matchPaneH, queryPeakListsH, matchPeakListsH, filterByTypeH))

    allMatches = []
    for queryPane, matchPane, queryPeakLists, matchPeakLists,  filterByType in data:
      if not (queryPane and queryPeakLists):
        continue
    
      peaks, dimDict, dimMappingDict = self.getMatchPeaks(ssShifts, queryPane, matchPane,
                                         queryPeakLists, matchPeakLists)
      
      if not peaks:
        name = self.getSpinSystemName(querySpinSystem)
        msg = 'No query peaks assigned to this spin system (%s)'
        showWarning('Warning',msg % name, parent=self)
        continue
       
      matches = findPositionMatches(peaks, dimDict, matchPeakLists,
                                    matchTolerance=None,
                                    separateSpinSystems=True)
      matches = self.filterMatches(matches, dimDict, querySpinSystem,
                                   queryPeakLists, matchPeakLists,
                                   filterByType)

      for score, position, peaks in matches:
        if peaks:
          # TBD: clunky this: assumes peakLists have same dim order; and since the match is a tuple it does
          # a big song and dance to get the position in the correct order
          peakList = peaks[0].peakList
          if peakList in dimMappingDict:
            dimMapping = dimMappingDict[peakList]
            pos = []
            # in 4D possibly the below is backwards, so possibly pos[dimMapping[n]] = position[n]
            # below is pos[n] = position[dimMapping[n]]
            for n, p in enumerate(position):
              pos.append(position[dimMapping[n]])
            for n in range(len(position)):
              position[n] = pos[n]

      allMatches.append((matches, dimDict))

    # Merge matches (also add spinSystem info)
    
    if allMatches:
      self.matches = self.mergeMatches(allMatches)
    else:
      self.matches = []
 
    self.updateMatches()
    self.setMatchPaneRegion()
    self.waiting = False

    # Reinstate drawing

    for popup, turnedOff in popups:
      if turnedOff:
        popup.turnDrawRequestsOn()
  
  def filterMatches(self, matches, dimDict, querySpinSystem, queryPeakLists,
                    matchPeakLists, filterByType):

    interResidueQuery = self.haveInterResiduePeakList(queryPeakLists)
    isPeakPositive = self.isPeakPositive
    getShiftPeaks = self.getShiftPeaks
    
    filteredMatches = []
    
    for match in matches:
      remove = False
      peaks  = match[2]
    
      for peak in peaks:
        peakList = peak.peakList
        experiment = peakList.dataSource.experiment
        shiftList = experiment.shiftList
        dim = dimDict[peakList]
        posSign = isPeakPositive(peak)
        peakDims = peak.sortedPeakDims()
        tolerance = getAnalysisDataDim(peakDims[dim-1].dataDim).assignTolerance
        resonances = set()
           
        for peakDim in peakDims:
          if peakDim.dim == dim:
            continue
        
          for contrib in peakDim.peakDimContribs:
            resonances.add(contrib.resonance)
        
        for resonance in resonances: 
          resonance  = contrib.resonance
          spinSystem = resonance.resonanceGroup
        
          if spinSystem and (spinSystem is querySpinSystem):
            remove = True
            break
    
        else:
          if 'N[coca]' in experiment.refExperiment.name:
            # Do not filter intra specific matches
            pass
        
          elif filterByType and interResidueQuery: # E.g. HNCOCA query
            shifts = [r.findFirstShift(parentList=shiftList) for r in resonances]
            qPeaks = getShiftPeaks(shifts, queryPeakLists)
 
            overlap = 0
            nonOverlap = 0
 
            for qPeak in qPeaks:
              if qPeak is peak:
                continue
            
              dimQ = dimDict[qPeak.peakList]
              peakDimQ = qPeak.findFirstPeakDim(dim=dimQ)
              ppmQ = peakDimQ.value
              posSignQ = isPeakPositive(qPeak)
 
              if (posSignQ == posSign) and (abs(peakDims[dimQ-1].value-ppmQ) < tolerance):
                overlap += 1
              
            if overlap:
              nonOverlap = 0
              mPeaks = getShiftPeaks(shifts, matchPeakLists)
 
              for mPeak in mPeaks:
                if (mPeak not in peaks) and (isPeakPositive(mPeak) == posSign):
                  dimM = dimDict[mPeak.peakList]
                  ppmM = mPeak.findFirstPeakDim(dim=dimM).value
                  if abs(peakDims[dim-1].value-ppmM) > tolerance:
                    nonOverlap += 1

              if nonOverlap:
                remove = True
                break

          if remove:
            break
       
        if remove:
          break
    
      if not remove:
        filteredMatches.append(match)

    return filteredMatches

  def mergeMatches(self, matchLists):
    # Merge by spin system
  
    spinSystemDict = {}
  
    for matches, dimDict in matchLists:
      for score, position, peaks in matches:
      
        spinSystems = set()
        
        peak = peaks[0]
        for peakDim in peak.peakDims:
          dim = dimDict[peak.peakList]
          if peakDim.dim != dim:
            for contrib in peakDim.peakDimContribs:
              spinSystems.add(contrib.resonance.resonanceGroup)
      
        spinSystems = frozenset(spinSystems)
        
        if spinSystemDict.get(spinSystems) is None:
          spinSystemDict[spinSystems] = []
        
        spinSystemDict[spinSystems].append( (score, position, peaks) )

    m = len(matchLists)
    mergedMatches = []
    for spinSystems in spinSystemDict:
      if not spinSystems:
        continue
        
      data = spinSystemDict[spinSystems]
      #if len(data) < m:
      #  continue

      scoreM = 0
      scoreP = 0
      positionM = [0.0] * len(data[0][1])
      peaksM = []
       
      for score, position, peaks in data:
        if not peaks:
          break
        
        scoreP += score
        scoreM += score/len(peaks)
        peaksM += peaks
        
        for i, pos in enumerate(position):
          positionM[i] += pos
      
      else:
        n = float(len(data))
 
        #scoreM /= n
        positionM = [p/n for p in positionM]
 
        spinSystem = list(spinSystems).pop()
        mergedMatches.append( (scoreM, scoreP, spinSystem, positionM, peaksM) )
 
    if self.matchOrder == MATCH_ORDERS[0]:
      def compFunc(a,b):
        return cmp(a[1], b[1])
      
    else:
      def compFunc(a,b):
        return cmp(a[0], b[0])
    
    mergedMatches.sort(compFunc)
    mergedMatches.reverse()
    
    return mergedMatches

  def isPeakPositive(self, peak):
  
    intensity  = peak.findFirstPeakIntensity(intensityType='height')
    if intensity and (intensity.value > 0):
      isPositive = True
    else:
      isPositive = False

    return isPositive

  def discardMatches(self):
    
    matches = self.matchMatrix.currentObjects
    
    for match in matches:
      self.matches.remove(match)

    if matches:
      self.updateMatches()

  def getPeakListsCompDims(self, peakLists):
  
    dimDict = {}
    
    for peakList in peakLists:
      spectrum = peakList.dataSource
      
      amide = set()
      for dimA, dimB in getOnebondDataDims(spectrum):
        dataDimRefA = getPrimaryDataDimRef(dimA)
        dataDimRefB = getPrimaryDataDimRef(dimB)
        isotopesA = dataDimRefA.expDimRef.isotopeCodes
        isotopesB = dataDimRefB.expDimRef.isotopeCodes
        
        if ('15N' in isotopesA) and ('1H' in isotopesB):
          amide.update((dimA, dimB))
        elif ('15N' in isotopesB) and ('1H' in isotopesA):
          amide.update((dimA, dimB))
      
      for dataDim in spectrum.dataDims:
        if dataDim not in amide:
          dimDict[peakList] = dataDim.dim

    return dimDict

  def setSeqLink(self):
      
    from ccpnmr.analysis.core.AssignmentBasic import addSpinSystemResonance, assignResToDim
    from ccpnmr.analysis.core.AssignmentBasic import addPeakResonancesToSeqSpinSystems, clearSeqSpinSystemLinks
    from ccpnmr.analysis.core.AssignmentBasic import mergeSpinSystems, mergeResonances
    self.administerNotifiers(self.unregisterNotify)
  
    if self.spinSystem and self.match:
      
      mScore, mScoreP, mSpinSystem, mPosition, matchPeaks = self.match
      if not mSpinSystem:
        msg = 'Matched peaks are not assigned to a spin system'
        showWarning('Failure', msg, parent=self)
        return
        
      qSpinSystem = self.spinSystem
      qPeaks = self.getShiftPeaks(self.shifts, self.queryPeakLists)
      nmrProject = qSpinSystem.topObject
      
      compDimDict = self.getPeakListsCompDims(list(self.queryPeakLists)+list(self.matchPeakLists))
      
      for peak in qPeaks:
        dim = compDimDict[peak.peakList]
        peak.cDim = peakDim = peak.findFirstPeakDim(dim=dim)
        peak.isotope = ','.join(peakDim.dataDimRef.expDimRef.isotopeCodes)
      
      queryPeakListsH, queryPeakListsC = self.categorisePeakLists(self.queryPeakLists)
              
      matchPairs  = []
      ambiguous   = []
      warnRefExpt = False
      spinSystems = []
      interMatches = 0
      
      for peak in matchPeaks:
        dim = compDimDict.get(peak.peakList)
        
        if dim is None:
          # Don't understand why this would happen, but it did in a course
          continue
        
        peakDim = peak.findFirstPeakDim(dim=dim)
        position = peakDim.value
        isotope = ','.join(peakDim.dataDimRef.expDimRef.isotopeCodes)
        peak.cDim = peakDim
        peak.isotope = isotope
        tolerance = getAnalysisDataDim(peakDim.dataDim).assignTolerance
        
        close = []
        for qPeak in qPeaks:
          if qPeak.isotope != isotope:
            continue
        
          delta = abs(position - qPeak.cDim.value)
          if delta < tolerance:
            close.append( [delta, qPeak] )
        
        if len(close) > 2:
          for i in range(len(close)-1):
            delta1, close1 = close[i]
            if delta1 is None:
              continue
          
            for j in range(i+1,len(close)):
              delta2, close2 = close[j]
              if delta2 is None:
                continue
          
              if close1.peakList is close2.peakList:
                if delta1 < delta2:
                  close[j][0] = None
                else:
                  close[i][0] = None  
       
        experiment = peak.peakList.dataSource.experiment
        expDimRefs = []
        for expDimRef1, expDimRef2 in getOnebondExpDimRefs(experiment):
          expDimRefs.append(expDimRef1)
          expDimRefs.append(expDimRef2)
        
        if not experiment.refExperiment:
          close = []
          warnRefExpt = True
          
        elif experiment.refExperiment.name in self.refNamesCo:
          interMatches += 1
        
        peak.close = []
        for delta, closePeak in close:
          if delta is not None:
            peak.close.append(closePeak)

      if mSpinSystem is qSpinSystem:
        msg = 'Query and match spin systems are the same'
        showWarning('Failure', msg, parent=self)
        return
      
      if ambiguous:
        msg = ''
        for peak in ambiguous:
          msg += ','.join([pd.annotation or '-' for pd in peak.sortedPeakDims])
          msg += ' '
 
        msg = '%d peak(s): %s matches multiple peaks' % (len(ambiguous),msg)
        msg += ' in the same peak list. These will not be linked'
        showWarning('Warning', msg, parent=self)
      
      if warnRefExpt:
        msg = 'Some of the matched peaks do not have reference'
        msg += ' experiment types set. These peaks cannot be linked'
        showWarning('Warning', msg, parent=self)
      
      for peak in matchPeaks:
        #print peak, peak.close
        if peak.close:
          experiment = peak.peakList.dataSource.experiment
          refExperiment = experiment.refExperiment
 
          if refExperiment:
            peakDims = peak.sortedPeakDims()
            tolerance = getAnalysisDataDim(peak.cDim.dataDim).assignTolerance
            
            cResonances = set()
            mResonances = set()
            qResonances = set()
            for contrib in peak.cDim.peakDimContribs:
              resonance = contrib.resonance
              cResonances.add(resonance)
              mResonances.add(resonance)

            for qPeak in peak.close:
              for contrib in qPeak.cDim.peakDimContribs:
                resonance = contrib.resonance
                cResonances.add(resonance)
                qResonances.add(resonance)

            resonances =  list(cResonances)
            mResonances = list(mResonances)
            qResonances = list(qResonances)
 
            isotopeCode = peak.isotope
            N = len(resonances)
            if N == 0:
              resonance = nmrProject.newResonance(isotopeCode=isotopeCode)

            elif N > 1:
              data = (N,isotopeCode,peak.cDim.value)
              msg  = 'There are %d resonances corresponding to %s position %.3f, ' % data
              msg += 'Merge resonances together?'
              if showYesNo('Warning', msg, parent=self):
                resonance = resonances[0]
                for resonanceB in resonances[1:]:
                  resonance = mergeResonances(resonance, resonanceB)
 
              else:
                if refExperiment.name in self.refNamesCo:
                  # Q i -> M i+1
                  
                  if qResonances:
                    resonance = qResonances[0]
                  else:
                    resonance = nmrProject.newResonance(isotopeCode=isotopeCode)
               
                else:
                  # Q i -> M i-1
                  
                  if mResonances:
                    resonance = mResonances[0]
                  else:
                    resonance = nmrProject.newResonance(isotopeCode=isotopeCode)

            else:
              resonance = resonances[0]
            
            resonance.resonanceGroup = None
            
            for qPeak in peak.close:
              if not qPeak.cDim.peakDimContribs:
                assignResToDim(qPeak.cDim,resonance,tolerance=tolerance)
             
            if not peak.cDim.peakDimContribs:
              assignResToDim(peak.cDim,resonance,tolerance=tolerance)
            
            if refExperiment.name in self.refNamesCo:
              # Q i -> M i+1
 
              for qPeak in peak.close:
                qRefExperiment = qPeak.peakList.dataSource.experiment.refExperiment
                if qRefExperiment and (qRefExperiment.name not in self.refNamesCo):
                  addSpinSystemResonance(qSpinSystem, resonance)
              
              clearSeqSpinSystemLinks(qSpinSystem, delta=+1)
              clearSeqSpinSystemLinks(mSpinSystem, delta=-1)
              
              seqOffsets = [None for pd in peak.sortedPeakDims()]
              seqOffsets[peak.cDim.dim-1] = -1
              
              addPeakResonancesToSeqSpinSystems(peak, seqOffsets)

            elif refExperiment.name not in self.refNamesCo:
              # Q i -> M i-1

              addSpinSystemResonance(mSpinSystem, resonance)
              
              clearSeqSpinSystemLinks(qSpinSystem, delta=-1)
              clearSeqSpinSystemLinks(mSpinSystem, delta=+1)

              for qPeak in peak.close:
                seqOffsets = [None for pd in qPeak.sortedPeakDims()]
                seqOffsets[qPeak.cDim.dim-1] = -1
                
                qRefExperiment = qPeak.peakList.dataSource.experiment.refExperiment
                if qRefExperiment and (qRefExperiment.name in self.refNamesCo):
                  addPeakResonancesToSeqSpinSystems(qPeak, seqOffsets)
                  
                elif interMatches == 0:
                  addPeakResonancesToSeqSpinSystems(qPeak, seqOffsets)

      for peak in qPeaks:
        del peak.cDim
        del peak.isotope
        
      for peak in matchPeaks:
        if hasattr(peak, 'cDim'):
          del peak.cDim
          del peak.close  
          del peak.isotope  

    self.updateSeqSegmentsAfter()
    self.administerNotifiers(self.registerNotify)
    self.update()

  def clearPrevLinks(self):
    
    self.clearLinks(-1)
    
  def clearNextLinks(self):
      
    self.clearLinks(1)
    
  def clearAllLinks(self):
  
    self.clearLinks(None)

  def clearLinks(self, delta):

    spinSystems = [s[0] for s in self.spinSystemMatrix.currentObjects]
    if spinSystems:
      if showOkCancel('Confirm','Are you sure?', parent=self):
        for spinSystem in spinSystems:
          clearSeqSpinSystemLinks(spinSystem, delta=delta)

  def nextLinkedSpinSystem(self):
      
    if self.spinSystem:
      nextSS  = findConnectedSpinSystem(self.spinSystem, delta=1)
      
      i = 0
      for spinSystem, shifts, pos in self.spinSystemMatrix.objectList:
        if nextSS is spinSystem:
          self.spinSystemMatrix.selectNthObject(i)
          break
        i += 1
    
  def prevLinkedSpinSystem(self):
      
    if self.spinSystem:
      prevSS  = findConnectedSpinSystem(self.spinSystem, delta=-1)
      i = 0
      for spinSystem, shifts, pos in self.spinSystemMatrix.objectList:
        if prevSS is spinSystem:
          self.spinSystemMatrix.selectNthObject(i)
          break
        i += 1
    
  def showPeaks(self):
      
    if self.spinSystem:
      peaksDict = {}
      for resonance in self.spinSystem.resonances:
        for contrib in resonance.peakDimContribs:
          peaksDict[contrib.peakDim.peak] = None

      peaks = peaksDict.keys() 
      if len(peaks) > 0:
        self.guiParent.viewPeaks(peaks)
    
  def showResonances(self):
      
    if self.spinSystem:
      self.parent.typeSpinSystem(spinSystem=self.spinSystem)    
          
  def selectSpinSystem(self, obj, row, col):
        
    if obj:
      self.spinSystem, self.shifts, self.positions = obj
      self.updateButtons()
      self.navigateToQuery()
      if self.spinSystem and self.matchSelect.get():
        self.findMatches()

      self.updateResTypesAfter()

  def makeRootRuler(self):
    
    if self.positions:
      for windowPane, posDict in self.positions:
        axisPanels = windowPane.axisPanels
 
        for label in posDict.keys():
          if label in ('x','y'):
            location = posDict[label].value
            axisPanel = windowPane.findFirstAxisPanel(label=label)
 
            if axisPanel:
               ruler = createRuler(location, axisPanel.panelType,
                                   color=self.color, lineWidth=1,
                                   dashLength=4, gapLength = 2, remove=False)
               self.rulers.append(ruler)

            break

  def navigateToRoot(self):
  
    for ruler in self.rootRulers:
      if not ruler.isDeleted:
        ruler.delete()
    
    self.rootRulers = []
    if self.rootPane and self.spinSystem and self.shifts:
      windowPane = self.rootPane
      position = {}
      
      xAxisPanel = windowPane.findFirstAxisPanel(label='x')
      yAxisPanel = windowPane.findFirstAxisPanel(label='y')
     
      for shift in self.shifts:
        isotope = shift.resonance.isotopeCode
        
        if isotope in xAxisPanel.axisType.isotopeCodes:
          position['x'] = shift.value
          ruler = createRuler(shift.value, xAxisPanel.panelType,
                              color=self.color, lineWidth=1,
                              dashLength=4, gapLength = 2, remove=False)
          self.rootRulers.append(ruler)

        elif isotope in yAxisPanel.axisType.isotopeCodes:
          position['y'] = shift.value
          ruler = createRuler(shift.value, yAxisPanel.panelType,
                              color=self.color, lineWidth=1,
                              dashLength=4, gapLength = 2, remove=False)
          self.rootRulers.append(ruler)
     
      windowFrame = windowPane.getWindowFrame()
      windowFrame.gotoPosition(position=position)

  def navigateToQuery(self):
      
    self.clearRulers()
  
    if self.spinSystem and self.queryPeakLists:
      queryPanes = [p for p in (self.queryPaneC, self.queryPaneH) if p]
      
      self.makeRootRuler()
      self.navigateToRoot()
      
      for windowPane in queryPanes:
        
        if windowPane is self.queryPaneC:
          focus = self.focusSelectC.get()
          matchPane = self.matchPaneC
        else:
          focus = self.focusSelectH.get()
          matchPane = self.matchPaneH
   
        orthoAxis = 'y'
        if windowPane.spectrumWindow.stripAxis == 'y':
          orthoAxis = 'x'
  
        windowFrame = windowPane.getWindowFrame()
 
        for windowPane0, position0 in self.positions:
          position = position0.copy()
          
          if windowPane0 is windowPane:
            for key in position:
              position[key] = position0[key].value
            
            windowFrame.gotoPosition(position=position)
            break
 
        if focus:
          orthoPpms = self.getQueryOrthoPositions(windowPane, matchPane)
          if orthoPpms:
            self.defineWindowRegions(windowPane, orthoPpms, orthoAxis)
 
        
  def setMatchPaneRegion(self):
     
    panes = ((self.queryPaneC, self.matchPaneC), (self.queryPaneH, self.matchPaneH))
    for queryPane, matchPane in panes:
 
      if queryPane and matchPane:
      
        orthoAxisQ = 'y'
        if queryPane.spectrumWindow.stripAxis == 'y':
          orthoAxisQ = 'x'

        orthoAxisM = 'y'
        if matchPane.spectrumWindow.stripAxis == 'y':
          orthoAxisM = 'x'
 
        axisPanelQ = queryPane.findFirstAxisPanel(label=orthoAxisQ)
        axisPanelM = matchPane.findFirstAxisPanel(label=orthoAxisM)
        if axisPanelQ and axisPanelM:
          regionsQ = [ar.region for ar in axisPanelQ.sortedAxisRegions()]
          regionsM = [ar.region for ar in axisPanelM.sortedAxisRegions()]
 
          if len(regionsQ) == len(regionsM):
            axisPanelOrtho = matchPane.findFirstAxisPanel(label=orthoAxisM)
 
            if axisPanelOrtho:
              for i, region in enumerate(regionsQ):
                if i >= len(axisPanelOrtho.axisRegions):
                  break
                axisPanelOrtho.sortedAxisRegions()[i].region = regionsQ[i]
                
 
  def updateWindowsAfter(self, window=None):
        
    self.after_idle(lambda :self.updateWindows(window)) # Axis panel notifiers need a bit of time


  def updateButtons(self):
    
    if self.spinSystem:
      for button in self.spinSystemButtons.buttons:
        button.enable()
      for button in self.spinSystemButtons2.buttons:
        button.enable()
      for button in self.spinSystemButtons3.buttons:
        button.enable()
    
    else:
      for button in self.spinSystemButtons.buttons:
        button.disable()
      for button in self.spinSystemButtons2.buttons:
        button.disable()
      for button in self.spinSystemButtons3.buttons:
        button.disable()

    if self.matches:
      self.matchButtons.buttons[3].enable()
      if self.match:
        self.matchButtons.buttons[1].enable()
        self.matchButtons.buttons[2].enable()
      else:
        self.matchButtons.buttons[1].disable()
        self.matchButtons.buttons[2].disable()
    
    else:
      self.matchButtons.buttons[3].disable()

    if self.spinSystem and (self.matchPaneH or self.matchPaneC):
      self.matchButtons.buttons[0].enable()
    else:
      self.matchButtons.buttons[0].disable()

  def updateWindows(self, window=None):
    
    self.updateQueryPanes()
    self.updateMatchPanes()

  def updatePeakDimContribAfter(self, contrib):
 
    if self.tabbedFrame.selected != 1:
      return

    if self.waiting:
      return
      
    self.updatePeakAfter(contrib.peakDim.peak)

  def updatePeakIntensityAfter(self, peakIntensity):
  
    if self.tabbedFrame.selected != 1:
      return

    if self.waiting:
      return
      
    self.updatePeakAfter(peakIntensity.peak)

  def updatePeakDimAfter(self, peakDim):
   
    if self.tabbedFrame.selected != 1:
      return

    if self.waiting:
      return
      
    self.updatePeakAfter(peakDim.peak)

  def updatePeakAfter(self, peak):
  
    if self.tabbedFrame.selected != 1:
      return

    if self.waiting:
      return
      
    peakList = peak.peakList
    
    if peakList in self.queryPeakLists:
      self.waiting = True
      self.navigateToQuery()
      self.after_idle(self.update)
      self.after_idle(self.findMatches)
    
    elif peakList in self.matchPeakLists:
      self.waiting = True
      self.after_idle(self.findMatches)


  def updateAfter(self, obj=None):
    
    if self.tabbedFrame.selected != 1:
      return
     
    if obj:
      if obj.className == 'Peak':
        if obj.peakList is not self.rootPeakList:
          return
          
      elif obj.className == 'PeakDim':
        if obj.peak.peakList is not self.rootPeakList:
          return

      else: # It's a resonance group
        self.updateResTypesAfter()
         
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)

  def getSpinSystemName(self, spinSystem):
        
    residue = spinSystem.residue
    if residue:
      tlc = getResidueCode(residue)
      
      chainCode = ''
      if len(residue.chain.molSystem.chains) > 1:
        chainCode = residue.chain.code + ' '
      
      name = '%s%d%s' % (chainCode,residue.seqCode,tlc)

    elif spinSystem.ccpCode:
      name = '{%d}%s' % (spinSystem.serial,spinSystem.ccpCode)
    
    elif spinSystem.name:
      name = '{%d:%s}' % (spinSystem.serial,spinSystem.name)
    
    else:
      name = '{%d}' % spinSystem.serial

    return name

  def update(self):
      
    textMatrix = []
    objectList = []

    getName = self.getSpinSystemName
    getSegment = self.getSpinSystemSeqSegment
    queryPeakLists = self.queryPeakLists
    getShiftPeaks = self.getShiftPeaks
    queryPanes = [p for p in (self.queryPaneC, self.queryPaneH) if p]
    
    if queryPeakLists:
      shiftList = list(queryPeakLists)[0].dataSource.experiment.shiftList
    
    else:
      shiftList = None
    
    if queryPanes:
      spinSystems = self.nmrProject.sortedResonanceGroups()
      
      for spinSystem in spinSystems:        
        positionDict = {}
        
        for queryPane in queryPanes:
          shiftDicts, yShifts = getSpinSystemWindowShifts([spinSystem,],
                                         queryPane, shiftList=shiftList)
          
          for shiftDict in shiftDicts:
            shifts = frozenset(shiftDict.values())
            peaks = getShiftPeaks(shifts, queryPeakLists)
            
            if not peaks:
              continue
          
            if positionDict.get(shifts) is None:
              positionDict[shifts] = []
              
            positionDict[shifts].append( (queryPane, shiftDict) )
        
        for shifts in positionDict:
          objectList.append((spinSystem, shifts, positionDict[shifts]))
 
    
      for spinSystem, shifts, positionDict in objectList:
        # '#','Name','i-1','i+1','i+0','Position'
        prevSS  = findConnectedSpinSystem(spinSystem, delta=-1)
        nextSS  = findConnectedSpinSystem(spinSystem, delta=1)
        otherSS = findConnectedSpinSystem(spinSystem, delta=0)
    
        prevName  = ''
        nextName  = ''
        otherName = ''
        
        if prevSS:
          prevName = getName(prevSS)
        if nextSS:
          nextName = getName(nextSS)
        if otherSS:
          otherName = getName(otherSS)    
    
        shiftValues = [s.value for s in shifts]
        shiftValues.sort()
        shiftValues = ['%.3f' % v for v in shiftValues]
    
        datum = [spinSystem.serial,
                 getName(spinSystem),
                 prevName,
                 nextName,
                 otherName,
                 '%d:%3.3d' % getSegment(spinSystem),
                 ','.join(shiftValues)]
       
        textMatrix.append(datum)
    

    self.spinSystemMatrix.update(textMatrix=textMatrix,objectList=objectList)
  
    self.updateMatches()
    self.updateButtons()
    self.updateResTypes()
    
    self.refreshSeqSeg = False
    self.waiting = False
    
  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)
    self.setAppDataOptions()
    BasePopup.destroy(self)

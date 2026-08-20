
"""
======================COPYRIGHT/LICENSE START==========================

EditPeakLists.py: Part of the CcpNmr Analysis program

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

from memops.general import Implementation

from memops.gui.Button          import Button
from memops.gui.ButtonList      import ButtonList, UtilityButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.Entry           import Entry
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Frame           import Frame
from memops.gui.Label           import Label
from memops.gui.LabelDivider    import LabelDivider
from memops.gui.MessageReporter import showOkCancel, showWarning, showError
from memops.gui.MultiWidget     import MultiWidget
from memops.gui.ProgressBar     import ProgressBar
from memops.gui.PulldownList    import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.TabbedFrame     import TabbedFrame
from memops.gui.Text            import Text

from memops.gui import Color

from ccpnmr.analysis.popups.BasePopup import BasePopup

from ccpnmr.analysis.core import ExperimentBasic 
from ccpnmr.analysis.core import PeakBasic
from ccpnmr.analysis.core import UnitConverter
#from ccpnmr.analysis.core import Util
from ccpnmr.analysis.core import WindowBasic 

from ccpnmr.analysis.core.AssignmentBasic import getShiftLists
from ccpnmr.analysis.core.ConstraintBasic import makePeaksFromConstraints
from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes

from ccpnmr.analysis.frames.PeakTableFrame import PeakTableFrame


def testEditPeakListsPopup(argServer):

  popup = EditPeakListsPopup(argServer.parent)
  popup.open()

def testSelectedPeaksPopup(argServer):

  popup = SelectedPeaksPopup(argServer.parent)
  popup.open()
   
THROUGH_SPACE_RESIDUE_LIMITS = (('Intra only',0), (u'\u00B11',1), (u'\u00B12',2))

class EditPeakListsPopup(BasePopup):
  """
  **Display, Edit and Create Peak Lists**
  
  This popup window collects together most of the data that relates to peaks
  picked within spectra and how they are organised within peak lists.
  Higher-level information that relates to the spectra or the experiments that
  generated the data are found elsewhere, in the Spectra_ and Experiments_
  tables. Each peak list relates to only one spectrum, but a spectrum can have
  multiple peak lists, although most spectra will normally have only one. Having
  more than one peak list for a spectrum allows different kinds of peaks to be
  separated, for example to distinguish between different peak picking
  techniques or real and synthetic peaks.

  This popup window is divided into three tab sections. The first shows all of  
  the peak lists within the current project; the second lists all of the peaks
  within a selected list and the third allows the user to make artificial or
  synthetic peak lists based upon things like chemical shift values and
  coordinate structures.

  **Peak Lists**
  
  This table shows all of the peaks lists within the project, the parameters
  that relate directly to them and how they fit with the spectrum/experiment
  hierarchy. Each peak list is identified by the experiment and spectrum to
  which it belongs, as well as a serial number for the list. Some of the values
  in the table may be edited to control a peak list's display and behavior. The
  most important of these is set in the "Active?" column. Each spectrum may have
  only one active peak list. The active peak list for a spectrum is the only
  peak list that will normally be assigned, have its peak selected in the
  spectrum display and will have newly picked peaks added to it. Inactive peak
  lists are displayed but do not usually participate in the manipulations that
  are done via the spectrum windows.

  Each peak list may be displayed in a different style, thus this table allows
  the user to change the symbol and colour of the peak list when it is displayed
  within spectra. The symbol refers to the cross mark that is placed to indicate
  the picked position of the peak; usually an extremum in spectrum intensity.
  The colour of the peak list relates to both the colour of the marker symbols
  and any textual annotation, e.g. to indicate resonance assignment.

  Below the peak list table a few functions that operate on whole peak lists are
  available. The [Add Sister List] is notable because it allows the user to make
  a new blank peak list, to which new peaks may be copied or added. Here, the
  new peak list is added to the same spectrum as the peak list that is currently
  selected in the table, hence the use of the word "Sister". The [Copy Peaks] is
  used to duplicate peaks, from one peak list to another, assuming the spectrum
  dimensions are compatible in terms of isotope types. The [Shift Whole Peak
  List] is useful if all of the peaks of a list need to be moved by the same
  amount relative to their spectrum. This is useful when peaks are imported with
  an offset, e.g. for matching TROSY and non-TROSY peaks, or if spectrum
  referencing has changed and peaks have not been properly located.

  **Peak Table**
  
  The idea behind this table is to give a textual listing of all the peaks in
  the selected peak list. The user can compare and edit peak parameters, e.g in
  terms of assignment or intensity, and can follow links to other kinds of
  information. For example the user can go from the selected peaks to the
  resonances that are assigned to those peaks. Double-clicking on the editable
  columns of the table allows many of the parameters to be directly edited,
  although discretion is recommended for certain operations like changing
  positions and intensity values. The pulldown menu above the peak table allows
  the user to specify which peak list will be displayed. The "Status" just below
  may be used to show only a subset of peaks, depending on how they are
  assigned. 

  Most of the options above the peak table provide functions that allow the user
  to locate and mark peaks within the spectrum display windows. The [Strip
  Selected] and [Find Peak] options are used to find particular peaks that have
  been selected in the table (with left-click +/- <Ctrl>/<Shift>). The former
  sub-divides the locations for separate peaks using strips (window dividers)
  and the latter locates only the last selected peak. The [Strip Locations] and
  [Go To Position] work in a similar way, except that they use only the
  *position information* from the selected peaks and don't necessarily display
  the actual selected peaks. The idea here is that peak positions can be used to
  move the display of any spectrum window that share at least some of the same
  kind of isotopes on their axes. For example a 2D 15N HSQC peak can be used to
  find an amide position (H & N) in 3D HNcoCA spectrum, thus locating two out of
  three axes. This functionality is very useful when some peak lists are used as
  "roots" to coordinate others.

  The buttons below the main peak table allow the user to edit and view many
  types of data associated with the peaks, although many parameters may be
  changed by double-clicking in the peak table. Some notable functions include: 
  [Unalias] which is used to move 'ghost' peaks to their real underlying ppm
  values by adding or subtracting a whole number of spectrum widths to their
  position. This is used when a peak lies outside the normal recorded bounds of
  the spectrum but nonetheless still appears within it because as an aliased
  signal;  [Assign] opens the `Assignment Panel`_ to control which resonances
  have been linked to the dimensions of the peak; to indicate what caused the
  peak. Such assignments may be to the resonances of specific atoms or
  resonances that are currently anonymous;  [Deassign] clears all of the
  resonance assignments to the peak dimensions. This does not affect how the
  resonances may be connected to atoms;  [Set Details] allows the user to set
  the "Details" column for all of the peaks selected in the table with a single
  operation; [Propagate Assign] spreads resonance assignments from one peak (the
  last selected) to the others selected in the table, which is useful even if
  not all of the peak dimensions are the same type, for example when spreading
  amide H & N resonances from an HSQC peak to 3D triple-resonance peaks.

  **Synthetic Lists**

  The final tab is used to make new peak lists where peak entities are predicted
  by artificial or synthetic means, i.e. not by direct inspection of spectrum
  intensities. If possible, peaks that are created by one of the synthetic
  methods and which are used as a source of evidence for NMR derived information
  should ultimately be related to the spectrum data; by re-centering the peaks
  and re-calculating their intensities. Four methods are currently available for
  synthesising new peak lists:

  The "From Shift Intersections" section allows the user to make a peak list
  based upon the intersection between chemical shifts that occur within the
  bounds of a given spectrum (the same spectrum that the peak list will belong
  to). A shift list is selected to specify which chemical shift values are used,
  and by default this shift list is the one that the spectrum uses during
  assignment (set at the experiment level). The synthesised peaks will correspond
  to all chemical shift values that match the spectrum data dimensions, in terms
  of isotope, and fit together to give a complete peak assignment. This shift
  matching process also considers whether the spectrum data dimensions represent
  only atom sites that are directly bound; thus removing many spurious peak
  locations. For example a 15N HSQC peak list may be made from 1H and 15N
  resonances that lie along the two  spectrum axes of the resonances and are
  covalently bound to one another. The through-space options control how far
  predictions are made for transfers like NOESY, i.e. in the absence of a
  structure you can limit the number of bonds and/or number of sequential
  residues to consider. It should be noted that this system cannot predict peaks
  for all kinds of spectra. It is currently limited to those kinds of experiment
  where all magnetisation transfers go through atom sites that are either
  recorded directly as a spectrum dimension (e.g. in HNCA, HH NOESY) or are
  intermediate between two recorded dimensions (e.g. HN(co)CA) and linked via
  J-coupling or one-bond transfer. More kinds of peak list may be predicted in
  the future. It is notable that for NOESY peak lists only peaks that represent
  connections between atoms that are either in the same residue or within five
  covalent bonds of one another are used.

  The "For Shifts and Structure" section makes synthetic peak lists using
  chemical shift intersections, as described above, but only allowing spectra
  that have a through-space connectivity (e.g. NOESY) between two of their
  dimensions. Also when comparing through-space connections, the possible peaks
  are filtered according to how close assigned atoms are within a structure
  ensemble. Any peaks that would be assigned to pairs of atoms that are further
  apart than the "Max Dist" value will not be made. The user can also reject
  certain peaks by insisting on a minimum spectrum intensity value at its
  position. 

  The "From Transposition" section is used to make peaks based upon reflection
  about a homonuclear diagonal. For example if peaks have been picked on only
  one side of a 2D NOESY experiment's diagonal, corresponding duplicates may be
  made on the other side of the diagonal where the assignments and positions of
  the two dimensions are swapped. Naturally, this functionality is limited to
  only spectra that have two data dimensions with of the same isotope type.

  The "From Distance Restraints" section is used to recreate NOE and other
  through-space peak lists based upon distance restraints, which may have been
  imported from outside CCPN. The individual distance restrains provide a list
  of atom pairs that are known to be close. The peaks are then synthesised from
  those by obtaining the chemical shift of the atom's resonances. Where
  relevant, the chemical shift of any covalently bound resonances are also
  obtained. The chemical shifts are then matched within the bounds of the
  selected spectrum, considering whether any of the data dimensions represent a
  'onebond' relationship that much be preserved in the peak assignments. It
  should be noted that a single distance restraint may give rise to more than
  one peak, when it is not possible to determine which restrained atom goes on
  which spectrum dimension. Also, ambiguous distance restraints will give rise
  to ambiguous peak assignments.

  **Caveats & Tips**
  
  If you need the peak positions to be displayed in Hz units or as data point
  positions, the "Position Unit" pulldown at the top may be changed.

  The [Set Details] function can be handy for marking peaks that cause violations
  in a structure calculation and thus need further attention.

  Peak lists may be merged by copying the peaks from one into another, although
  this takes no account of duplication.
  
  When making synthetic NOESY (or other though-space) peaks, having a residue
  limit of 1 or 2 will naturally not have an effect if the bond limit is too
  short, i.e. where the wouldn't be enough bonds to get to the next residues.

  .. _Spectra: EditSpectrumPopup.html
  .. _Experiments: EditExperimentPopup.html
  .. _`Edit Peak`: EditPeakPopup.html
  .. _`Assignment Panel`: EditAssignmentPopup.html

  """
  def __init__(self, parent, *args, **kw):

    self.popups = {}
    self.guiParent = parent
    self.shiftsSpectrum = None
    self.structSpectrum = None
    self.transposeSource = None
    self.structure = None
    self.constraintList = None
    self.constraintSet = None
    self.constraintSpectrum = None
    self.labellingScheme = True
    self.shiftList = None
    self.strucShiftList = None
    self.molSystem = None
 
    self.peakList = None
    self.cloneSourceList = None
    self.reproduceSourceList = None
    self.subtractPeakList = None
    self.peakListTable = None
    self.getShiftDeltasWidget = None
    
    BasePopup.__init__(self, parent=parent, title="Peak : Peak Lists", **kw)

  def open(self):
  
    self.toggleTab(self.tabbedFrame.selected)
    
    BasePopup.open(self)

  def body(self, guiFrame):

    self.geometry('770x500')
    
    row = 0
    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)
    
    tipTexts = ['A table of all the peak lists within the current CCPN project',
                'A table of the individual peaks within one peak list',
                'Tools to generate peak lists by artificial means; from shifts, structures, transposition & restraints']
    options = (' Peak Lists ',' Peak Table ',' Synthetic Lists ')
    tabbedFrame = TabbedFrame(guiFrame, options=options,
                              callback=self.toggleTab,
                              grid=(0,0), tipTexts=tipTexts)
    self.tabbedFrame = tabbedFrame
    frameA, frameB, frameC = tabbedFrame.frames
    
    
    self.editDetailsEntry = Entry(self, text='', returnCallback=self.setDetails, width=12)
    colors = Color.standardColors
    
    self.peakColorWidget = PulldownList(self, callback=self.setColor,
                                        texts=[c.name for c in colors],
                                        objects=[c.hex for c in colors],
                                        colors=[c.hex for c in colors])
                                        
    symbols = ['+','x','o','*']                                        
    self.peakSymbolWidget = PulldownList(self, callback=self.setSymbol,
                                         texts=symbols)
    
    frameA.grid_columnconfigure(0, weight=1)
    frameA.grid_rowconfigure(0, weight=1)
    
    tipTexts = ['The name of the experiment that gave rise to the spectrum & hence contains the peak list',
                'The name of the spectrum to which the peak list pertains',
                'The serial number of the peak list, within its spectrum',
                'Sets whether the peak list is the active one for its spectrum; always "yes" if there is only one list for a spectrum',
                'The colour to use for the peak cross/symbol and assignment annotation; each list may have a different colour',
                'The kind of symbol (a small mark on the screen) used to indicate picked peak locations for the list',
                'The number of peaks in the peak list',
                'The percentage of peak dimensions in the whole peak list that carry assignments',
                'Whether the peak list is synthetic/simulated; if so, peak assignments in the list have very little influence on chemical shift averages',
                'A textual comment for the peak, often user-supplied']
                
    colHeadings = ['Experiment','Spectrum','List','Active?','Color',
                   'Symbol','No. Peaks','% Assigned','Synthetic?','Details']
                   
    editWidgets      = [None, None,
                        None, None,
                        self.peakColorWidget,
                        self.peakSymbolWidget, None,
                        None, None,
                        self.editDetailsEntry]
    editGetCallbacks = [None, None, 
                        None, self.setActive,
                        self.getColor,
                        self.getSymbol, None,
                        None, None,
                        self.getDetails]
    editSetCallbacks = [None, None,
                        None, None,
                        self.setColor,
                        self.setSymbol, None,
                        None, None,
                        self.setDetails]
                        
    self.peakListTable = ScrolledMatrix(frameA, tipTexts=tipTexts,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        initialRows=10, grid=(0,0),
                                        headingList=colHeadings,
                                        callback=self.selectCell,
                                        deleteFunc=self.deleteList)
                                        
    self.peakListTable.doEditMarkExtraRules = self.isActivateable
    
    tipTexts = ['Show a table of the individual peaks within the selected peak list',
                'Delete the selected peak list and all its peaks',
                'Add a new, blank peak list in the same spectrum as the selected peak list',
                'Copy the peaks (positions, intensities and assignments) from the selected peak list to another with matching dimensions',
                'Create a new peak list which contains peaks in the first selected peak list which are not in the second selected peak list, based on position (assignments not copied)',
                'Move all of the peaks in the selected peak list by specified offsets for each dimension; offsets are prompted']
    texts    = ['Edit Peaks','Delete','Add Sister List',
                'Copy Peaks','Subtract Peaks','Shift Whole Peak List']
    commands = [self.editPeaks,self.deleteList,
                self.addList,self.reproduceList,
                self.subtractList,self.getShiftDeltas]
    self.peakListButtons  = ButtonList(frameA, texts=texts, tipTexts=tipTexts,
                                       commands=commands, grid=(1,0))
    self.getShiftDeltasButton = self.peakListButtons.buttons[4]                               
    
    # Peaks Table
    
    frameB.expandGrid(0,0)
    
    self.peakTableFrame = PeakTableFrame(frameB, self.guiParent,
                                         nmrProject=self.nmrProject,
                                         grid=(0,0))
    
    # Synthetic lists
    
    frameC.grid_columnconfigure(0, weight=1, minsize=150)
    frameC.grid_columnconfigure(1, weight=1, minsize=200)
    frameC.grid_rowconfigure(1, weight=1)
    frameC.grid_rowconfigure(3, weight=1)
    frameC.grid_rowconfigure(5, weight=1)
    frameC.grid_rowconfigure(7, weight=1)
    
    div = LabelDivider(frameC, text='From Shift Intersections', grid=(0,0), gridSpan=(1,2))
    
    tipText = 'Make a synthetic peak list using the intersection of chemical shift values, according to the specified settings'
    self.predFromShiftsButton = Button(frameC, borderwidth=1, text='Predict from shifts',
                                       command=self.synthesiseFromShifts,
                                       grid=(1,0), sticky='nsew', tipText=tipText)
    self.predFromShiftsButton.config(bg='#B0FFB0') 
    
    subFrame = Frame(frameC, grid=(1,1))
    subFrame.expandGrid(5,6)
    
    # From shifts

    label = Label(subFrame,text='Spectrum: ', grid=(0,0))
    tipText = 'Selects which spectrum to make synthetic peaks for; a new peak list is generated'
    self.shiftsSpectrumPulldown = PulldownList(subFrame, self.changeShiftsSpectrum,
                                               grid=(0,1), tipText=tipText)

    label = Label(subFrame,text='Shift List: ', grid=(1,0))
    tipText = 'Sets which list to use as the source of chemical shift values; peaks will be made at relevant shift intersections'
    self.shiftListPulldown = PulldownList(subFrame, self.changeShiftList,
                                          grid=(1,1), tipText=tipText)

    label = Label(subFrame,text='Mol System: ', grid=(2,0))
    tipText = 'Allows the considered shifts to be limited to only those assigned to a given molecular system (group of chains)'
    self.molSystemPulldown = PulldownList(subFrame, self.changeMolSystem,
                                          grid=(2,1), tipText=tipText)
    
    label = Label(subFrame,text='Use unassigned?', grid=(3,0))
    tipText = 'Whether to make peaks for chemical shift intersections involving one or more unassigned resonances'
    self.useUnassignedSelect = CheckButton(subFrame, grid=(3,1),
                                           selected=False, tipText=tipText)

    label = Label(subFrame, text='Through-space\nbond limit:', grid=(4,0))
    tipText = 'When dealing with through-space transfers like, NOESY or DARR, which atom pairs to include, based on the maximum number if intervening bonds'
    values = [None] + list(range(1,10))
    texts = [str(x) for x in values]
    self.throughSpaceBondsPulldown = PulldownList(subFrame, None, texts, values,
                                                  index=0, grid=(4,1), tipText=tipText)

    label= Label(subFrame, text = 'Isotope\nLabelling:', grid=(0,2), gridSpan=(2,1))
    tipText = 'If required, selects the isotope labelling specification to filter possible chemical shift intersections & hence peak creation'
    self.labellingSchemePulldownA = PulldownList(subFrame, self.changeLabellingScheme,
                                                grid=(0,3), gridSpan=(2,1), tipText=tipText)

    label= Label(subFrame, text='Min Isotope\nFraction:', grid=(2,2), gridSpan=(2,1))
    tipText = 'The minimum spin active isotope incorporation for an atom site and correlation to be considered, according to the selected labelling'
    self.minFractionEntryA = FloatEntry(subFrame, text=0.25, width=8,
                                       grid=(2,3), gridSpan=(2,1), tipText=tipText) 

    label = Label(subFrame, text='Through-space\nresidue limit:', grid=(4,2))
    texts = [x[0] for x in THROUGH_SPACE_RESIDUE_LIMITS]
    objs = [x[1] for x in THROUGH_SPACE_RESIDUE_LIMITS]
    tipText = 'When dealing with through-space transfers like, NOESY or DARR, which atom pairs to include based on residue separation'
    self.throughSpaceResiduePulldown = PulldownList(subFrame, None, texts, objs,
                                                    index=1, grid=(4,3), tipText=tipText)
                                       
    # from structure                                   
    
    div = LabelDivider(frameC, text='From Shifts and Structure', grid=(2,0), gridSpan=(1,2))
                                          
    tipText = 'Make a synthetic peak list using structural distance to filter possible chemical shift intersections'
    self.predFromStructButton = Button(frameC, borderwidth=1, text='Predict from structure',
                                       command=self.synthesiseFromStruct,
                                       grid=(3,0), sticky='nsew', tipText=tipText)
    self.predFromStructButton.config(bg='#B0B0FF')     
    
    subFrame = Frame(frameC, grid=(3,1))
    subFrame.expandGrid(5,6)

    label = Label(subFrame,text='Spectrum: ', grid=(0,0))
    tipText = 'Selects which spectrum to make distance filtered synthetic peaks for; a new peak list is generated'
    self.structSpectrumPulldown = PulldownList(subFrame, self.changeStructSpectrum,
                                               grid=(0,1), tipText=tipText)

    label = Label(subFrame,text='Structure: ', grid=(1,0))
    tipText = 'Selects which structure ensemble to use for calculating atomic distances'
    self.structurePulldown = PulldownList(subFrame, self.changeStructure,
                                          grid=(1,1), tipText=tipText)

    label = Label(subFrame,text='Max Dist: ', grid=(2,0))
    tipText = 'Sets the Angstrom threshold below which atoms may be considered for making peaks; if their shift intersection lies in the spectrum bounds'
    self.distThresholdEntry = FloatEntry(subFrame, text=5.0, width=8,
                                         grid=(2,1), tipText=tipText)

    label = Label(subFrame,text='Shift List: ', grid=(3,0))
    tipText = 'Sets which chemical shift list to use for generating potential peak locations'
    self.strucShiftListPulldown = PulldownList(subFrame, self.changeStrucShiftList,
                                               grid=(3,1), tipText=tipText)

    label= Label(subFrame, text='Min Spectrum\nValue (Height):', grid=(4,0))
    tipText = 'Selects a threshold for the spectrum intensity at a potential peak location, if the intensity is below the value no peak is made'
    self.minHeightEntry = FloatEntry(subFrame, text='', width=16, grid=(4,1), tipText=tipText)

    label= Label(subFrame, text = 'Isotope\nLabelling:', grid=(0,2), gridSpan=(2,1))
    
    tipText = 'If required, selects a isotope labelling scheme to apply further filtering for possible chemical shift intersections'
    self.labellingSchemePulldownB = PulldownList(subFrame, self.changeLabellingScheme,
                                                grid=(0,3), gridSpan=(2,1), tipText=tipText)

    label= Label(subFrame, text='Min Isotope\nFraction:', grid=(2,2), gridSpan=(2,1))
    tipText = 'The minimum spin active isotope incorporation for an atom site and correlation to be considered, according to the selected scheme'
    self.minFractionEntryB = FloatEntry(subFrame, text=0.25, width=8,
                                       grid=(2,3), gridSpan=(2,1), tipText=tipText)


    # Transpose
   
    div = LabelDivider(frameC, text='From Transposition', grid=(4,0), gridSpan=(1,2))

    tipText = 'Make a peak list, in the same spectrum as the selected list, where peak positions in homonuclear dimensions are swapped relative to the source'
    self.transposeButton = Button(frameC, borderwidth=1, text='Make transpose list',
                                  command=self.makeTranspose, 
                                  grid=(5,0), sticky='nsew', tipText=tipText)
    self.transposeButton.config(bg='#B0FFB0')   
                  
    subFrame = Frame(frameC, grid=(5,1))
    subFrame.expandGrid(3,1)

    label = Label(subFrame,text='Source Peak List: ', grid=(0,0))
    tipText = 'Selects which peak list to use as the source of the peak position information; the new peak list will share the same spectrum'
    self.transposeSourcePulldown = PulldownList(subFrame, self.changeTransposeSource,
                                                grid=(0,1), tipText=tipText)
    
    
    # From constraints
   
    div = LabelDivider(frameC, text='From Distance Restraints', grid=(6,0), gridSpan=(1,2))
    
    tipText = 'Make a synthetic peak list using information from a distance restraint list'
    self.constraintsButton = Button(frameC, borderwidth=1,
                                    text='Make from restraints',
                                    command=self.makeFromConstraints,
                                    grid=(7,0), sticky='nsew', tipText=tipText)
    self.constraintsButton.config(bg='#B0B0FF')  
    
    subFrame = Frame(frameC, grid=(7,1))
    subFrame.expandGrid(3,3)

    label = Label(subFrame,text='Spectrum: ', grid=(0,0))
    tipText = 'Selects which spectrum to make a new peak list in; the experiment for this spectrum states which shift list provides peak locations'
    self.constraintSpectrumPulldown = PulldownList(subFrame,  self.changeConstraintSpectrum,
                                                   grid=(0,1), tipText=tipText)

    label = Label(subFrame,text='Restraint Set: ', grid=(1,0))
    tipText = 'Selects which set of restraints the input distance restraint list is contained within'
    self.constraintSetPulldown = PulldownList(subFrame, self.changeConstraintSet,
                                              grid=(1,1), tipText=tipText)

    label = Label(subFrame,text='Distance List: ', grid=(2,0))
    tipText = 'Selects which distance restraint list will be used to provide information about which correlations will be made in the new peak list'
    self.constraintListPulldown = PulldownList(subFrame, self.changeConstraintList,
                                               grid=(2,1), tipText=tipText)
   
    # Main Table

    bottomButtons = UtilityButtonList(tabbedFrame.sideFrame,
                                      helpUrl=self.help_url)
    bottomButtons.grid(row=0, column=0, sticky='e')
 
    #self.cloneButton = self.peakListButtons.buttons[3]
    self.reproduceButton = self.peakListButtons.buttons[3]
    self.subtractButton = self.peakListButtons.buttons[4]
 
    self.waiting = False
    self.update()
    
    self.administerNotifiers(self.registerNotify)
  
  def toggleTab(self, i):
  
    if i == 0:
      self.update()
    
    elif i == 1:
      self.peakTableFrame.update()
    
    elif i == 2:
      self.updateShiftsSpectrum()
      self.updateShiftLists()
      self.updateMolSystems()
      self.updateStructSpectrum()
      self.updateConstraintSpectra()
      self.updateTransposePeakLists()
      self.updateStructure()
      self.updateStrucShiftLists()
      self.updateLabellingSchemes()
      self.updateConstraintSets()
      self.updateConstraintLists()
 
    
  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete',):
      for clazz in ('ccp.molecule.ChemCompLabel.LabelingScheme',):
        notifyFunc(self.updateLabellingSchemes, clazz, func)

    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.PeakDimContrib','ccp.nmr.Nmr.Peak',
                    'ccp.nmr.Nmr.PeakList', 'ccp.nmr.Nmr.DataSource',
                    'ccp.nmr.Nmr.Experiment'):
        notifyFunc(self.updateAfter, clazz, func)

    for func in ('setTextColor','setSymbolColor','setSymbolStyle',):
      notifyFunc(self.updateAfter, 'ccpnmr.Analysis.AnalysisPeakList')

    for func in ('setName'):
      for clazz in ('ccp.nmr.Nmr.PeakList',
                    'ccp.nmr.Nmr.DataSource',
                    'ccp.nmr.Nmr.Experiment'):
        notifyFunc(self.updateAfter, clazz, func)
    
    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.DataSource', 'setActivePeakList')
        
    for func in ('__init__', 'delete','setName'):
      notifyFunc(self.updateStructSpectrum, 'ccp.nmr.Nmr.DataSource', func)
      notifyFunc(self.updateShiftsSpectrum, 'ccp.nmr.Nmr.DataSource', func)
      notifyFunc(self.updateShiftLists, 'ccp.nmr.Nmr.ShiftList', func)
      notifyFunc(self.updateStrucShiftLists, 'ccp.nmr.Nmr.ShiftList', func)
      notifyFunc(self.updateConstraintSpectra, 'ccp.nmr.Nmr.DataSource', func)
      notifyFunc(self.updateConstraintLists, 'ccp.nmr.NmrConstraint.DistanceConstraintList', func)

    for func in ('__init__', 'delete'):
      notifyFunc(self.updateMolSystems, 'ccp.molecule.MolSystem.MolSystem', func)
      notifyFunc(self.updateStructure, 'ccp.molecule.MolStructure.StructureEnsemble', func)
      notifyFunc(self.updateTransposePeakLists, 'ccp.nmr.Nmr.PeakList', func)
      notifyFunc(self.updateConstraintSets, 'ccp.nmr.NmrConstraint.NmrConstraintStore', func)

    # Peak Table
    
    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.PeakDimContrib',):
        notifyFunc(self.peakTableFrame.contribUpdateAfter, clazz, func)

    for func in ('__init__', 'delete', 'setName'):
      for clazz in ('ccp.nmr.Nmr.DataSource', 'ccp.nmr.Nmr.Experiment'):
        notifyFunc(self.peakTableFrame.updatePeakListsAfter, clazz, func)
        
    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.PeakList',):
        notifyFunc(self.peakTableFrame.updatePeakListsAfter, clazz, func)

    for func in ('__init__', 'delete','setName'):
      notifyFunc(self.peakTableFrame.updateWindowListsAfter, 'ccpnmr.Analysis.SpectrumWindow', func)

    for func in ('__init__', 'delete','setAnnotation','setDetails','setFigOfMerit','setFitMethod'):
      notifyFunc(self.peakTableFrame.updatePeaksAfter, 'ccp.nmr.Nmr.Peak', func)
    for func in ('setAnnotation','setPosition','setNumAliasing','setLineWidth'):
      notifyFunc(self.peakTableFrame.peakDimUpdateAfter, 'ccp.nmr.Nmr.PeakDim', func)
    for func in ('__init__', 'delete', 'setValue'):
      notifyFunc(self.peakTableFrame.intensityUpdateAfter, 'ccp.nmr.Nmr.PeakIntensity', func)
    for func in ('__init__', 'delete'):
      notifyFunc(self.peakTableFrame.updateStructures, 'ccp.molecule.MolStructure.StructureEnsemble', func)

     
  
  def updateAnalysisPeakList(self, analysisPeakList):
  
    self.updateAfter()
      
  def changeLabellingScheme(self, scheme):

    if scheme is not self.labellingScheme:
      self.labellingScheme = scheme
      self.updateLabellingSchemes()
  
  def updateLabellingSchemes(self, notifyScheme=None):
  
    names = []
    index = 0
    
    scheme = self.labellingScheme
    
    names = ['Automatic from sample', '<None>']
    schemes = [True, None]
    
    schemes0 = self.getLabellingSchemes()
    schemes  += schemes0
    names += ['Scheme: %s' % s.name for s in schemes0]
    
    if scheme not in schemes:
      scheme = schemes[0]
    
    index = schemes.index(scheme)
    
    if scheme is not self.labellingScheme:
      self.labellingScheme = scheme
      self.updateAfter()  
    
    self.labellingSchemePulldownA.setup(names, schemes, index)     
    self.labellingSchemePulldownB.setup(names, schemes, index)     
  
  def getLabellingSchemes(self):
  
    return self.project.sortedLabelingSchemes()
    
  def makeFromConstraints(self):

    if self.constraintList and self.constraintSpectrum:
      peaks = makePeaksFromConstraints(self.constraintList.sortedConstraints(), self.constraintSpectrum)
  
      if peaks:
        self.editPeaks(peaks[0].peakList)
  
  def changeConstraintSpectrum(self, spectrum):
    
    if spectrum:
      self.constraintSpectrum = spectrum

  def changeConstraintSet(self, constraintSet):
    
    if constraintSet is not self.constraintSet:
      self.constraintSet = constraintSet
      self.updateConstraintLists()
  
  def changeConstraintList(self, constraintList):
     
    self.constraintList = constraintList
  
  def updateConstraintSets(self, obj=None):
    
    index = -1
    constraintSets = self.getConstraintSets()
    names = ['%d' % (cs.serial) for cs in constraintSets]
    
    if constraintSets:
    
      if self.constraintSet not in constraintSets:
        self.constraintSet = constraintSets[0]
    
      index = constraintSets.index(self.constraintSet)
    
    self.constraintSetPulldown.setup(names, constraintSets, index)
    self.updateConstraintLists()
    
  def updateConstraintLists(self, obj=None):
    
    index = -1
    constraintLists = self.getDistanceConstraintLists()
    names = ['%d:%s' % (cl.serial, cl.name or '') for cl in constraintLists]
    
    if constraintLists:
    
      if self.constraintList not in constraintLists:
        self.constraintList = constraintLists[0]
    
      index = constraintLists.index(self.constraintList)
    
    self.constraintListPulldown.setup(names, constraintLists, index)
  
  def getDistanceConstraintLists(self):

    constraintLists = []
    
    if self.constraintSet:
      for constraintList in self.constraintSet.constraintLists:
        if constraintList.className == 'DistanceConstraintList':
          constraintLists.append(constraintList)
    
    return constraintLists

  def getConstraintSets(self):
  
    return self.nmrProject.sortedNmrConstraintStores()
  
  def makeTranspose(self):
  
    if self.transposeSource:
      peakList = PeakBasic.makeTransposePeakList(self.transposeSource)

      if peakList:
        self.editPeaks(peakList)


  def updateTransposePeakLists(self, obj=None):
      
    index     = -1    
    peakLists = self.getTransposePeakLists()
    names     = ['%s:%s%d' % (pl.dataSource.experiment.name, pl.dataSource.name, pl.serial) for pl in peakLists]
    
    if peakLists:
    
      if self.transposeSource not in peakLists:
        self.transposeSource = peakLists[0]
    
      index = peakLists.index(self.transposeSource)
    
    else:
      self.transposeSource = None
      
    self.transposeSourcePulldown.setup(names, peakLists, index)


  def updateShiftsSpectrum(self, obj=None):
    
    index   = -1    
    spectra = self.getShiftsSpectra()
    names   = ['%s:%s' % (s.experiment.name, s.name) for s in spectra]
    
    if spectra:
    
      if self.shiftsSpectrum not in spectra:
        self.shiftsSpectrum = spectra[0]
    
      index = spectra.index(self.shiftsSpectrum)
    
    else:
      self.shiftsSpectrum = None
    
    self.shiftsSpectrumPulldown.setup(names, spectra, index)

  def updateShiftLists(self, obj=None):
    
    index = 0   
    shiftLists = getShiftLists(self.nmrProject)
    names = [sl.name or 'ShiftList %d' % sl.serial for sl in shiftLists]
    
    if shiftLists:
    
      if self.shiftList not in shiftLists:
        self.shiftList = shiftLists[0]
    
      index = shiftLists.index(self.shiftList)
    else:
      self.shiftList = None  
    
    self.shiftListPulldown.setup(names, shiftLists, index)

  def updateMolSystems(self, obj=None):
    
    index = 0   
    molSystems = self.project.sortedMolSystems()
    names = ['<Any>',] + [ms.code for ms in molSystems]
    molSystems.insert(0, None)
    
    if self.molSystem not in molSystems:
      self.molSystem = molSystems[0]
    
    index = molSystems.index(self.molSystem)
    
    self.molSystemPulldown.setup(names, molSystems, index)

  def updateStructSpectrum(self, obj=None):
    
    index   = -1
    spectra = self.getStructSpectra()
    names   = ['%s:%s' % (s.experiment.name, s.name) for s in spectra]
    
    if spectra:
    
      if self.structSpectrum not in spectra:
        self.structSpectrum = spectra[0]
    
      index = spectra.index(self.structSpectrum)
    
    else:
      self.structSpectrum = None
    
    self.structSpectrumPulldown.setup(names, spectra, index)

  def updateConstraintSpectra(self, obj=None):
    
    index   = -1
    spectra = self.getConstraintSpectra()
    names   = ['%s:%s' % (s.experiment.name, s.name) for s in spectra]
    
    if spectra:
    
      if self.constraintSpectrum not in spectra:
        self.constraintSpectrum = spectra[0]
    
      index = spectra.index(self.constraintSpectrum)
    
    else:
      self.constraintSpectrum = None  
    
    self.constraintSpectrumPulldown.setup(names, spectra, index)


  def updateStructure(self, obj=None):
    
    index      = -1
    structures = self.getStructures()
    names      = ['%s:%d' % (s.molSystem.code,s.ensembleId) for s in structures]
    
    if structures:
    
      if self.structure not in structures:
        self.structure = structures[0]
        
      index = structures.index(self.structure)  

    else:
      self.structure = None
    
    self.structurePulldown.setup(names, structures, index)

  def updateStrucShiftLists(self, obj=None):
    
    index = 0   
    shiftLists = getShiftLists(self.nmrProject)
    names = [sl.name or 'ShiftList %d' % sl.serial for sl in shiftLists]
    
    if shiftLists:
    
      if self.strucShiftList not in shiftLists:
        self.strucShiftList = shiftLists[0]
    
      index = shiftLists.index(self.strucShiftList)
    
    else:
      self.strucShiftList = None
    
    self.strucShiftListPulldown.setup(names, shiftLists, index)

  def changeTransposeSource(self, peakList):
  
    self.transposeSource = peakList
  
  def changeShiftsSpectrum(self, spectrum):
  
    self.shiftsSpectrum = spectrum
    self.shiftList = spectrum.experiment.shiftList 
    self.updateShiftLists()

  def changeShiftList(self, shiftList):
  
    self.shiftList = shiftList

  def changeMolSystem(self, molSystem):
  
    self.molSystem = molSystem

  def changeStructSpectrum(self, spectrum):
  
    self.structSpectrum = spectrum
    self.strucShiftList = spectrum.experiment.shiftList 
    self.updateStrucShiftLists()

  def changeStructure(self, structure):
  
    self.structure = structure

  def changeStrucShiftList(self, shiftList):
  
    self.strucShiftList = shiftList

  def getTransposePeakLists(self):
  
    peakLists = []
    for experiment in self.nmrProject.sortedExperiments():
      for spectrum in experiment.sortedDataSources():
        if ExperimentBasic.getEquivalentDataDims(spectrum):
          for peakList in spectrum.peakLists:
            peakLists.append(peakList)
            
    return peakLists

  def getShiftsSpectra(self):
  
    throughCoExceptions = set(['H[N[co[CA]]]', 'H[N[co[{CA|ca[C]}]]]', 'H[N[ca[CO]]]'])
    
    spectra = []
    for experiment in self.nmrProject.sortedExperiments():
      for spectrum in experiment.sortedDataSources():
        for dataDim in spectrum.dataDims:
          if dataDim.className != 'FreqDataDim':
            break
        
          dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)

          if not dataDimRef.expDimRef.expTransfers:
            if not experiment.refExperiment:
              break
            
            if experiment.refExperiment.name not in throughCoExceptions:
              break
        
        else:
          spectra.append(spectrum)
  
    return spectra

  def getStructSpectra(self):
  
    spectraDict = {}
    peakLists = ExperimentBasic.getThroughSpacePeakLists(self.nmrProject.root,
                                                       excludeSimulated=False)
    for peakList in peakLists:
      spectraDict[peakList.dataSource] = True
    
    spectra = spectraDict.keys()
      
    return spectra

  def getConstraintSpectra(self):
  
    spectraDict = {}
    peakLists = ExperimentBasic.getThroughSpacePeakLists(self.nmrProject.root,
                                                       excludeSimulated=False)
    for peakList in peakLists:
      spectraDict[peakList.dataSource] = True
    
    spectra = spectraDict.keys()
      
    return spectra

  def getStructures(self):
  
    structures = []
    if self.structSpectrum and self.structSpectrum.experiment.molSystems:
      for molSystem in self.structSpectrum.experiment.sortedMolSystems():
        for structure in molSystem.sortedStructureEnsembles():
          structures.append(structure)
      
    else:
      for molSystem in self.project.sortedMolSystems():
        for structure in molSystem.sortedStructureEnsembles():
          structures.append(structure)
  
    return structures

  def synthesiseFromShifts(self):

    if not self.shiftList.measurements:
      showError('Failure','Shift list is empty', parent=self)
      return

    if self.shiftsSpectrum:
      if not self.shiftsSpectrum.experiment.refExperiment:
        showError('Failure','Experiment does not have experiment type set', parent=self)
        return

      if (self.labellingScheme is True) and not self.shiftsSpectrum.experiment.labeledMixtures:
        labelling = None
      else:
        labelling = self.labellingScheme
    
      self.waiting = True
      progressBar   = ProgressBar(self, text='Synthesize from shifts')
      useUnassigned = self.useUnassignedSelect.get()
      lThreshold = self.minFractionEntryA.get() or 0.0
      
      bondLimit = self.throughSpaceBondsPulldown.getObject()
      residueLimit = self.throughSpaceResiduePulldown.getObject()
      
      peakList = PeakBasic.makePeakListFromShifts(self.shiftsSpectrum,
                                                  useUnassigned=useUnassigned,
                                                  progressBar=progressBar,
                                                  shiftList=self.shiftList,
                                                  molSystem=self.molSystem,
                                                  bondLimit=bondLimit,
                                                  residueLimit=residueLimit,
                                                  labelling=labelling,
                                                  labellingThreshold=lThreshold)
      progressBar.destroy()
      if peakList:
        self.editPeaks(peakList)
      
      
  
  def synthesiseFromStruct(self):
  
    if self.structSpectrum and self.structure:
      if (self.labellingScheme is True) and not self.structSpectrum.experiment.labeledMixtures:
        labelling = None
      else:
        labelling = self.labellingScheme
        
      self.waiting = True
      progressBar = ProgressBar(self, text='Synthesize from structure')
      threshold = self.distThresholdEntry.get() or 5.0
      lThreshold = self.minFractionEntryB.get() or 0.0
      minHeight =  self.minHeightEntry.get() or None
      peakList = PeakBasic.structurePredictNoePeakList(self.structure, self.structSpectrum, 
                                                       distThreshold=threshold,
                                                       progressBar=progressBar,
                                                       labelling=labelling,
                                                       labellingThreshold=lThreshold,
                                                       minHeight=minHeight,
                                                       shiftList=self.strucShiftList)
      progressBar.destroy()
      if peakList:
        self.editPeaks(peakList)

  def isActivateable(self, peakList, row, col):

    if col == 3:
      activePeakList = peakList.dataSource.activePeakList
      if peakList is activePeakList:
        return False
      else:
        return True
    else:
      return True

  def setActive(self, peakList):
  
    peakList.dataSource.activePeakList = peakList

  def getColor(self, peakList):
 
    color = peakList.analysisPeakList.symbolColor
    self.peakColorWidget.set(color)
 
  def setColor(self, null):

    color = self.peakColorWidget.getObject()
    peakList = self.peakList
    if peakList:
      peakList.analysisPeakList.symbolColor = color
      peakList.analysisPeakList.textColor = color
  
  def getSymbol(self, peakList):
 
    symbol = peakList.analysisPeakList.symbolStyle
    self.peakSymbolWidget.set(symbol)
 
  def setSymbol(self, null):

    symbol = self.peakSymbolWidget.getObject()

    peakList = self.peakList
    if peakList:
      peakList.analysisPeakList.symbolStyle = symbol

  def setDetails(self, event):
  
    text = self.editDetailsEntry.get()
    if text:
      self.peakList.details = text
      self.update()
  
  def getDetails(self,peakList):
  
    if peakList:
      self.editDetailsEntry.set(peakList.details)
  
  def updateAfter(self, obj=None):
    
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
  
  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)

  def cancelCopy(self):
  
    if self.cloneSourceList or self.reproduceSourceList:
      self.cloneSourceList = None
      self.reproduceSourceList = None
      #self.cloneButton.config(text='Clone Peaks', bg='grey90')
      self.reproduceButton.config(text='Copy Peaks', bg='grey90')

  def cancelSubtract(self):
  
    if self.subtractPeakList:
      self.subtractPeakList = None
      self.subtractButton.config(text='Subtract Peaks', bg='grey90')

  def selectCell(self, peakList, row, col):
  
    self.peakList = peakList
    if self.cloneSourceList:
      pointsA = [dd.numPoints for dd in self.peakList.dataSource.sortedDataDims()]
      pointsB = [dd.numPoints for dd in self.cloneSourceList.dataSource.sortedDataDims()]
      if pointsA != pointsB:
        showWarning('Stop','Peak clone not possible. Spectra must have the same number of points for each dimension', parent=self)
        return  
    
      N = len(self.cloneSourceList.peaks)
      if N > 2000:
        msg = 'Copy will duplicate %d peaks. This may take some time.'
        if not showOkCancel('Warning', msg % N, parent=self):
          self.cancelCopy()
          return 
          
      editPeaksPopup = self.parent.popups.get('edit_peak_list')
      if editPeaksPopup:
        editPeaksPopup.turnOffNotifiers()
      self.turnOffNotifiers()
      self.cloneButton.config(text='Clone Peaks', bg='grey90')

      progressBar = ProgressBar(self, text="Copying peaks ",total = len(self.cloneSourceList.peaks))
      PeakBasic.copyPeakList(self.cloneSourceList,self.peakList,progressBar=progressBar,rePick=False)
      progressBar.destroy()

      if editPeaksPopup:
        editPeaksPopup.update()
        editPeaksPopup.turnOnNotifiers()
      self.turnOnNotifiers()
      self.cancelCopy()
      self.update()
      
    elif self.reproduceSourceList:
      N = len(self.peakList.peaks)
      if N > 500:
        if not showOkCancel('Warning','Copy will duplicate %d peaks. This may take some time.' % N, parent=self):
          self.cancelCopy()
          return 

      editPeaksPopup = self.parent.popups.get('edit_peak_list')
      if editPeaksPopup:
        editPeaksPopup.turnOffNotifiers()
      self.turnOffNotifiers()
      self.reproduceButton.config(text='Copy Peaks', bg='grey90')

      rePick = True
      if self.reproduceSourceList.dataSource is self.peakList.dataSource:
        rePick = False

      progressBar = ProgressBar(self, text="Copying peaks ",total = len(self.reproduceSourceList.peaks))
      PeakBasic.copyPeakList(self.reproduceSourceList,self.peakList,progressBar=progressBar,rePick=rePick)
      progressBar.destroy()
      
      if editPeaksPopup:
        editPeaksPopup.update()
        editPeaksPopup.turnOnNotifiers()
        
      self.turnOnNotifiers()
      self.cancelCopy()
      self.update()
      
    elif self.subtractPeakList:
      if self.peakList.dataSource is not self.subtractPeakList.dataSource:
        showError
        showError('Different spectra','For now can only subtract peaks from peak list in same spectrum', parent=self)

      N = len(self.subtractPeakList.peaks)
      if N > 500:
        if not showOkCancel('Warning','Subtraction is checking %d peaks. This may take some time.' % N, parent=self):
          self.cancelSubtract()
          return 

      editPeaksPopup = self.parent.popups.get('edit_peak_list')
      if editPeaksPopup:
        editPeaksPopup.turnOffNotifiers()
      self.turnOffNotifiers()
      self.subtractButton.config(text='Subtract Peaks', bg='grey90')

      progressBar = ProgressBar(self, text="Subtracting peaks ",total = len(self.subtractPeakList.peaks))
      peakListNew = PeakBasic.subtractPeakLists(self.subtractPeakList,self.peakList, progressBar=progressBar)
      progressBar.destroy()
      
      if editPeaksPopup:
        editPeaksPopup.update()
        editPeaksPopup.turnOnNotifiers()
        
      self.turnOnNotifiers()
      self.cancelSubtract()
      self.update()
      
    for button in self.peakListButtons.buttons[0:5]:
      button.enable()

  def cloneList(self):
  
    if self.peakList and not self.cloneSourceList:
      if showWarning('Info','Now select a destination list', parent=self):
        self.cloneSourceList = self.peakList
        self.cloneButton.config(text='Cancel Cloning', bg='lightBlue')
        return
    
    self.cancelCopy()

  def reproduceList(self):
  
    if self.peakList and not self.reproduceSourceList:
      if showWarning('Info','Now select a destination list', parent=self):
        self.reproduceSourceList = self.peakList
        self.reproduceButton.config(text='Cancel Copy', bg='lightBlue')
        return
    
    self.cancelCopy()
  
  def subtractList(self):
  
    if self.peakList and not self.subtractPeakList:
      if showWarning('Info','Now select a list (from same spectrum) to subtract', parent=self):
        self.subtractPeakList = self.peakList
        self.subtractButton.config(text='Cancel Subtract', bg='lightBlue')
        return
    
    self.cancelSubtract()
  
  def shiftPeakList(self, deltas):
  
    if self.getShiftDeltasWidget and self.peakList:
      self.getShiftDeltasWidget.place_forget()
      self.getShiftDeltasWidget = None
      
      if deltas is not None:
        for peak in self.peakList.peaks:
          for i, peakDim in enumerate(peak.sortedPeakDims()):
            if peakDim.value is not None:
              peakDim.value += deltas[i]
      
  
  def getShiftDeltas(self):
    
    self.cancelCopy()
    
    if self.peakList:
      frame = self.tabbedFrame.frames[0]
      frame2 = self.peakListButtons
      button = self.getShiftDeltasButton
      x = frame.winfo_x() + frame2.winfo_x() + button.winfo_x() + button.winfo_width()
      y = frame.winfo_y() + frame2.winfo_y() + button.winfo_y() + button.winfo_height()
    
      spectrum = self.peakList.dataSource
      nDim = spectrum.numDim
      isotopes = getSpectrumIsotopes(spectrum)
      options = [u'\u0394 ppm F%d (%s):' % (i+1, iso) for i, iso in enumerate(isotopes)]
      
      widget = MultiWidget(self, FloatEntry, useImages=False,
                           callback=self.shiftPeakList,
                           minRows=nDim, maxRows=nDim,
                           relief='raised', borderwidth=2,
                           options=options, values=[0.0] * nDim,
                           defaultValue=0.0)
 
      widget.place(x=x, y=y, anchor='se')
      
      self.getShiftDeltasWidget = widget
    
  def addList(self):
                 
    self.cancelCopy()
    if (self.peakList):
      spectrum = self.peakList.dataSource
      spectrum.newPeakList(details = 'New list')
          
  def deleteList(self, *event):

    self.cancelCopy()
    
    if self.peakList:
      spectrum   = self.peakList.dataSource
      experiment = spectrum.experiment
      if len(spectrum.peakLists) == 1:
        data = (experiment.name, spectrum.name)
        showWarning('Delete failed',
                    'Cannot delete the only peak list in spectrum %s:%s' % data, parent=self)
        return
        
      data = (experiment.name, spectrum.name, self.peakList.serial)
      if not showOkCancel('Confirm','Really delete %s:%s peak list %d?' % data, parent=self):
        return
        
      self.peakList.delete()
      self.peakList = None
    
  def update(self):
  
    getColorScheme = self.analysisProfile.findFirstColorScheme
    peakLists = []
    textMatrix = []
    colorMatrix = []
    for expt in self.nmrProject.sortedExperiments():
      for spec in expt.sortedDataSources():
        if spec.dataType != 'processed':
          continue
      
  	for peakList in spec.sortedPeakLists():
          peakLists.append(peakList)
          
    for i, peakList in enumerate(peakLists):
      peakList = peakLists[i]
      spectrum = peakList.dataSource
      peaks    = peakList.peaks
      length   = len(peaks)
      totalPossAssign = spectrum.numDim * length
      actualAssigned = 0
      
      details = peakList.details
      if not details:
        details = ' '
      
      for peak in peaks:
        for dim in peak.peakDims:
          if dim.peakDimContribs:
            actualAssigned += 1
            
      if totalPossAssign > 0:
        percent = '%3.1f' % (100*actualAssigned/float(totalPossAssign) )
      else:
        percent = '0.0'
        
      colors   = [None] * 10
      if peakList.dataSource.activePeakList is peakList:
        isActive = 'Yes'
        colors[3] = '#a0d0a0'
      else:  
        isActive = 'No'
      
      if peakList.isSimulated:
        isSynthetic = 'Yes'
        colors[8] = '#ff8080'
      else:
        isSynthetic = 'No' 
      
      analysisPeakList = peakList.analysisPeakList
      
      scheme = getColorScheme(colors=(analysisPeakList.symbolColor,))
      if scheme:
        plColor = scheme.name
      else:
        plColor = analysisPeakList.symbolColor      
      
      symbol  = analysisPeakList.symbolStyle
      colors[4] = analysisPeakList.symbolColor
       
      colorMatrix.append(colors) 
      textMatrix.append([spectrum.experiment.name,
                         spectrum.name,
                         peakList.serial,
                         isActive,
                         plColor,
                         symbol,
                         length,
                         percent,
                         isSynthetic,
                         details])
    
    if not self.peakList: 
      for button in self.peakListButtons.buttons[0:5]:
        button.disable()
        
    self.peakListTable.update(objectList=peakLists,
                               textMatrix=textMatrix,
                               colorMatrix=colorMatrix)

    self.waiting = False

 
  def editPeaks(self, peakList=None):
  
    self.cancelCopy()
    peakList = peakList or self.peakList
    
    if peakList:
      self.peakTableFrame.update(peakList)
      self.tabbedFrame.select(1)


    
class SelectedPeaksPopup(BasePopup):
  """
  **Show a Subset of Peaks in Tabular Form**
  
  The idea behind this table is to give a textual listing of a particular subset
  of peaks. This table may be shown in two separate contexts; to represent all
  the peak entities which may be selected in one or more spectrum window
  displays inside Analysis, or it may be used by the various [Show Peaks]
  functions inside Analysis. When used for the display of in-spectrum
  selections, this table may be opened from the main Analysis menu or via the
  "Peak::Selection Table" option of the right-click window (shortcut is "s").
  Also, the contents of this table changes as the selected peaks change, thus
  representing the current selection status. 

  The listed peaks may be from several different spectra of any dimensionality
  or type. Putting the peak selection into a table of this kind allows the user
  to make easy comparisons and editing of their parameters, e.g in terms of
  assignment or intensity, and also provides links to other kinds of
  information. For example the user can go from the selected peaks to the
  resonances that are assigned to those peaks. Double-clicking on the editable
  columns of the table allows many of the parameters to be directly edited,
  although discretion is recommended for certain operations like changing
  positions and intensity values.

  **Spectrum Window Navigation**

  The options above the peak table provide functions that allow the user to 
  locate and mark peaks within the spectrum display windows. The [Strip
  Selected] and [Find Peak] options are used to find particular peaks that have
  been selected in the table (with left-click +/- <Ctrl>/<Shift>). The former
  sub-divides the locations for separate peaks using strips (window dividers)
  and the latter locates only the last selected peak. The [Strip Locations] and
  [Go To Position] work in a similar way, except that they use only the
  *position information* from the selected peaks and don't necessarily display
  the actual selected peaks. The idea here is that peak positions can be used to
  move the display of any spectrum window that share at least some of the same
  kind of isotopes on their axes. For example a 2D 15N HSQC peak can be used to
  find an amide position (H & N) in 3D HNcoCA spectrum, thus locating two out of
  three axes. This functionality is very useful when some peak lists are used as
  "roots" to coordinate others.

  **Button Functions**
  
  The buttons below the main peak table allow the user to edit and view many
  types of data associated with the peaks, although many parameters may be
  changed by double-clicking in the peak table. Some notable functions include: 
  [Unalias] which is used to move 'ghost' peaks to their real underlying ppm
  values by adding or subtracting a whole number of spectrum widths to their
  position. This is used when a peak lies outside the normal recorded bounds of
  the spectrum but nonetheless still appears within it because as an aliased
  signal;  [Assign] opens the `Assignment Panel`_ to control which resonances
  have been linked to the dimensions of the peak; to indicate what caused the
  peak. Such assignments may be to the resonances of specific atoms or resonances
  that are currently anonymous;  [Deassign] clears all of the resonance
  assignments to the peak dimensions. This does not affect how the resonances
  may be connected to atoms;  [Set Details] allows the user to set the "Details"
  column for all of the peaks selected in the table with a single operation;
  [Propagate Assign] spreads resonance assignments from one peak (the last
  selected) to the others selected in the table, which is useful even if not all
  of the peak dimensions are the same type, for example when spreading amide H &
  N resonances from an HSQC peak to 3D triple-resonance peaks.

  **Caveats & Tips**
  
  If you need the peak positions to be displayed in Hz units or as data point
  positions, the "Position Unit" pulldown at the top may be changed.

  The [Set Details] function can be handy for marking peaks that cause violations
  in a structure calculation and thus need further attention.

  .. _`Edit Peak`: EditPeakPopup.html
  .. _`Assignment Panel`: EditAssignmentPopup.html

  """

  def __init__(self, parent, peaks=None, *args, **kw):

    self.guiParent  = parent
    self.peaks = []
         
    BasePopup.__init__(self, parent=parent, title="Peak : Selected Peaks", **kw)

    if peaks is not None:
      self.update(peaks)

  def open(self):
  
    self.peakTableFrame.updatePeaksAfter()
    BasePopup.open(self)

  def update(self, peaks=None):
  
    if peaks is not None:
      peaks = list(peaks)
    
      self.peaks = peaks
      self.peakTableFrame.peaks = peaks
    
    self.peakTableFrame.selectStatus = 'Any'
    self.peakTableFrame.updatePeaksAfter()
    

  def body(self, guiFrame):
    
    self.geometry('700x500')
    
    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)
    
    self.peakTableFrame = PeakTableFrame(guiFrame, self.guiParent, peaks=self.peaks)
    self.peakTableFrame.grid(row=0, column=0,sticky='nsew')
    
    frame = Frame(self.peakTableFrame.topFrame, grid=(0,9), gridSpan=(1,2))
    frame.expandGrid(None,0)
    
    texts    = []
    commands = []
    bottomButtons = UtilityButtonList(frame, commands=commands, grid=(0,1),
                                      texts=texts, helpUrl=self.help_url)
    
    self.administerNotifiers(self.registerNotify)
    
  def administerNotifiers(self, notifyFunc):
  
    # Peak Table
    
    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.PeakDimContrib',):
        notifyFunc(self.peakTableFrame.contribUpdateAfter, clazz, func)

    for func in ('__init__', 'delete','setName'):
      for clazz in ('ccp.nmr.Nmr.PeakList', 'ccp.nmr.Nmr.DataSource', 'ccp.nmr.Nmr.Experiment'):
        notifyFunc(self.peakTableFrame.updatePeakListsAfter, clazz, func)

    for func in ('__init__', 'delete','setName'):
      notifyFunc(self.peakTableFrame.updateWindowListsAfter, 'ccpnmr.Analysis.SpectrumWindow', func)

    for func in ('__init__', 'delete','setAnnotation','setDetails','setFigOfMerit'):
      notifyFunc(self.peakTableFrame.updatePeaksAfter, 'ccp.nmr.Nmr.Peak', func)
    for func in ('setAnnotation','setPosition','setNumAliasing','setLineWidth'):
      notifyFunc(self.peakTableFrame.peakDimUpdateAfter, 'ccp.nmr.Nmr.PeakDim', func)
    for func in ('__init__', 'delete', 'setValue'):
      notifyFunc(self.peakTableFrame.intensityUpdateAfter, 'ccp.nmr.Nmr.PeakIntensity', func)
    for func in ('__init__', 'delete'):
      notifyFunc(self.peakTableFrame.updateStructures, 'ccp.molecule.MolStructure.StructureEnsemble', func)

  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)
     

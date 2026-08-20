"""
======================COPYRIGHT/LICENSE START==========================

BrowseConstraints.py: Part of the CcpNmr Analysis program

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
import Tkinter, re

from memops.general import Implementation

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.Entry import Entry
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Frame import Frame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.LabelDivider import LabelDivider
from memops.gui.Label import Label
from memops.gui.MessageReporter import showOkCancel, showWarning, showInfo, showYesNo
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.popups.BasePopup       import BasePopup
from ccpnmr.analysis.core.ConstraintBasic import makeNmrConstraintStore, \
                                     splitConstraintListAmbigUnambig, mergeConstraintLists, \
                                     makePeaksFromConstraints, splitConstraintListViol

from ccpnmr.analysis.core.AssignmentBasic import makeResonanceGuiName
from ccpnmr.analysis.core.ConstraintBasic import updateDistConstraintFromPeakAssign, mergeDuplicateConstraints
from ccpnmr.analysis.core.ConstraintBasic import updatePeaksFromConstraints, getStructureViolations
from ccpnmr.analysis.core.ConstraintBasic import exportAriaTbl, getConstraintStoreResonances
from ccpnmr.analysis.core.StructureBasic  import getAtomSetsDistance, getAtomSetsDihedral

# TBD add double click zooming
# Group and colour code restraint buttons

VIOL_COLOR_A = '#C0FF80'
VIOL_COLOR_B = '#FFFF80'
VIOL_COLOR_C = '#FFC080'
VIOL_COLOR_D = '#FF8080'

class BrowseConstraintsPopup(BasePopup):
  """
  **Display Structure Restraint and Restraint Violation Data**

  This popup window is used to manage and display structural restraints. Mostly
  such restraints are derived from NMR data and are applied during a structure
  calculation to restrict molecular conformations to those which are consistent
  with the observed data. Many kinds of NMR data can be used to generate
  restraints and the restraints themselves are of different types. This
  system will display tables for any kind or restraint described in the CCPN
  data model, but the most common kind of restraints are for distances and
  dihedral (torsion) angles.

  This system is used for viewing and managing restraint information on a large
  scale, it is not used directly for the initial generation of the restraints
  themselves. Structural restraints may be imported into CCPN via the 
  FormatConverter_ or they may be created with various various dedicated
  components available in CCPN software. For example distance restraints may be
  made via the `Make Distance Restraints`_ (using through space experiments like
  NOESY) and `Make H Bond Restraints`_ options. Dihedral angle restraints can be
  made from chemical shift information using DANGLE_ and from scalar couplings
  using the `3J H-Ha Coupling`_ option.

  **Restraint Sets**  

  The first "Restraint Sets" tabs lists all of the sets into which the restraint
  lists of the project are grouped into. The table indicates how many restraint
  and violation lists are in the set as well as an indication of the number of
  "fixed" resonance assignments that have been frozen in the set. Restraint sets
  may be deleted here, but are created elsewhere; when restraint lists are
  created they can be put into a new or existing set.

  Within CCPN, restraints are stored in restraint lists which group together
  restraints of the same type (i.e. angle and distance restraints cannot be
  mixed in a list). Restraint lists are further grouped together into restraint
  sets. A restraint set often represents all of the restraint lists that will be
  used together at the same time during a structure calculation. However, this
  grouping is slightly more subtle because a restraint set has important
  consequences with regards to NMR resonance assignment. In essence, a restraint
  set takes a fixed snapshot of the resonance assignment status at the time the
  restraints were made. Consequently, if atomic assignments change after
  restraints are generated the restraints will remain linked to the original
  atoms. This is helpful to the user because it is always possible to know what
  was actually restrained in a structure calculation. The user should however be
  careful when putting new restraints into old restraint sets; if atom
  assignments have changed then even the new restraints (actually between NMR
  resonances) will still restrain the old atom assignments. Thus it is good
  practice to always make a new restraint set if assignments have changed.

  **Restraint Lists**
  
  The second "Restraint Lists" tab displays all of the restraint lists that are
  available in the selected restraint set. The restraint set may be selected in
  the pulldown menu above the tabs. Generally this table is used to give an
  overview of the restraint lists but several high-level functions are
  accessible via the buttons below the table. As described above, the actual
  restraint generation takes place elsewhere in dedicated systems.

  Restraint lists my be deleted and merged/combined together into a single list
  (if they are of the same type). It is notable that restraint lists may  also
  be split. Currently there are two means of doing this, the first is to [Split
  Ambig/Unambig] which separates a distance restraint list into two, where one
  list contains only restraints that are equivalent to a single assignment (only
  two linked atoms or prochiral sets) and the other list has restraints that
  represent multiple, ambiguous assignments. The [Split Violated] function is
  helpful during violation analyses, which aim to determine which restraints are
  incompatible (e.g. mistaken assignment, minor conformation or artifact) with a
  structural model. This will separate out a list of those restraints that are
  violated in a structure (beyond the restraint bounds) for further inspection.

  **Restraints**
  
  The third "Restraints" tab allows the user to view all of the individual
  restraints within a restraint list, selected via the left hand "Restraint
  List" pulldown. If any structure violation analyses have been performed on the
  list the user can select from the "Violation List" pulldown. Any restraints
  that  were violated in the structural analysis will be coloured; red, orange
  or yellow depending on severity. Also where violations are recorded the
  columns for "Mean Viol"; average violation amount over a structure ensemble,
  and "Viol Fraction"; the proportion of an ensemble's models that were violated,
  are filled.

  Selecting from the "Structure" pulldown menu allows geometric information to
  be extracted from a structure ensemble, so that it may be compared with the
  values that are being restrained. Accordingly, if a structure violation
  analysis is performed using [Calculate Violations] then it is the selected
  structure that is compared with the restraints to find inconsistencies. Also,
  selecting as structure allows the "Struc Value" column to be filled with the
  value calculated from the resonance locations in the structure. Here   the kind
  of value naturally differs according to which type of restraint  list is
  displayed.

  The "Value Method" pulldown menu is notable for distance restraints because it
  allows the user to switch between "NOE sum" and "Minimum" options. These are
  relevant because a restraint, which restrains pairs of resonances, may link
  groups that contain multiple atoms (e.g. a methyl group). With multiple atoms
  at the end of a distance restraint there is flexibility about how the distance
  is defined. The "Minimum" method to measuring distances simply records the
  shortest distance between any pair of atoms from either side of the restraint
  (e.g. the closes atoms in two methyls). The "NOE sum" method uses r^-6
  distance summation to give a value that represents what the equivalent NOE
  estimated distance would be if a signal were recorded between the two atom
  groups. For example, if a restraint is to a methyl group all three methyl
  atoms will contribute to increase the intensity of an NOE signal, so the
  distance will appear to *shorter* than the single atom equivalent. Hence it is
  this shorter distance that will come from the NMR data and appear in the
  restraint, and thus the shorter distance that should be used in violation
  analyses. The get the "NOE sum" distance the pairs of atom distances
  are first converted into a kind of 'NOE intensity' (r^-6), these are added
  together to give the total 'NOE intensity', which is then converted back
  into a distance.

  In the main restraints table the restrains are listed as coloured rows, but
  there will often be grey rows present. A grey row represents an alternative
  set of restrained resonances (e.g. a different pair of atoms) but relates to
  the same restrained value. The main coloured row for a restraint and the grey
  alternatives are both referred to as the "items" for the restraint. Hence at
  the bottom there is the [Delete Items] button to remove specific restraint
  alternatives. In this regard there is no difference between the first green
  item and the other possibilities; it is fairly arbitrary which item comes
  first. Naturally deleting all of the items (possibilities) in a restraint
  removes the whole restraint. Although the first item is often coloured green
  it may have other colours indicating a structural violation. A restraint with
  more than one item may reflect either a genuine signal overlap, e.g. an NOE
  peak is caused by two or more close atom pairs that have similar chemical
  shifts and cannot be separated, or ambiguous trial assignments that come from
  speculatively matching chemical shifts to spectrum peaks.

  The various functions below the table allow restraints to be managed but the
  restrained resonances, and the parameter values they relate, to come from
  specialised restraint generation systems. In general the restraint information
  is derived from NMR experimentation and not adjusted by hand, Several functions
  allow the user to link between the restraints and any spectrum peaks that they
  are derived from. This can be done after selecting restraints with [Show
  Peaks] or after selecting peaks with [Show Restraints For Selected Peaks].
  There are also functions to coordinate assignments between the restraints and
  peaks. How this is done in practice tends to vary with personal preference,
  but it is possible to both refresh a restraint according to the latest peak
  assignment with [Update Assignment From Peak] or refresh the peaks according
  to the restraints (e.g. after an ARIA or CYANA run that removes restraint
  ambiguity) with [Update Peak Assignments].
      
  **Violation Lists**
  
  The forth "Violation Lists" tab gives details of all of the violation analyses
  that have been recorded in the project. Strictly speaking, each violation list
  can contain results from several analyses but this is usually not the case in
  normal Analysis operation. The individual violation records are presented in
  the last "Violations" tab, and also  affect the violation data presented in
  the "Restraints" table. Violation lists may be calculated, given a structure
  and a restraint list, inside Analysis using [Calculate Violations] in the
  "Restraints" tab, or they may be imported from external programs like ARIA.

  A violation list is a grouping of violation entities, and these record an
  inconsistency between the restrained values of a restraint and what is
  actually observed in a structure (or ensemble). Essentially this means a
  structural value is outside the bonds of the restraint. Clicking on a
  violation list in the upper table will fill the lower table with an overview
  of the structure models that were used in the restraint violation analysis
  (usually just a single ensemble).

  **Violations**
  
  The last tab lists all of the violation items that were recorded in a 
  selected violation list. This table will present some of the same information
  as the "Restraints" tab, if a violation list is selected there. The data is
  merely presented in a different way. Selecting a specific row in the table
  allows the user to [Show Restraint], to view the restraint that was violated
  in the "Restraints" tab. Other functions allow the user to jump to the peak,
  that gave rise to the restraint, which may be the underlying source of the
  problem. 

  **Caveats & Tips**

  Having large restraint sets in a CCPN project can make loading the CCPN data
  into Analysis slow, and takes up more memory. Thus it is advisable to delete
  any restraint sets that are no longer in use,

  If the restraint table is slow to update, consider switching off any structure
  selection. Extracting values from large structure ensembles can take a
  noticeable amount of time.

  The [Calculate Violations] function will not work until a Structure is
  selected.

  Although the CCPN violation analysis tries to do a good job it can never be
  exactly the same as analyses performed during a structure calculation. Most of
  the reason for this is due to the stereo-specific resolution of prochiral
  atoms. For example, Analysis can only determine which of HBa or HBb atoms
  actually goes with the HB2 or HB3 sites by taking a poll (based on minimising
  violation) after the structure calculation is complete, and this prochiral
  resolution may differ from what actually occurred in the calculation
  (minimising energy).

  .. _FormatConverter: FormatConverter.html
  .. _`Make Distance Restraints`: CalcDistConstraintsPopup.html
  .. _`Make H Bond Restraints`: MakeHbondRestraintsPopup.html
  .. _DANGLE: DangleGui.html
  .. _`3J H-Ha Coupling`: CalcHnHaCouplingPopup.html
    
  """

  def __init__(self, parent, *args, **kw):


    self.guiParent      = parent
    self.constraintSet  = None
    self.violationListA = None
    self.violationListB = None
    self.violationListC = None
    self.violation      = None
    self.constraintListA = None
    self.constraintListB = None
    self.structure      = None
    self.item           = None
    self.distMethod     = 'noe'
    self.quickStrucDict = {}
        
    self.waiting            = False
    self.waitingViolations  = False
    self.waitingRestraints  = False
    self.waitingViolationLists = False
    self.waitingConstraintLists = False

    BasePopup.__init__(self, parent=parent, title="Structure : Restraints and Violations", **kw)

  def body(self, guiFrame):
    
    self.geometry("700x600")
    
    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(1, weight=1)

    # Main - top

    frame = Frame(guiFrame)
    frame.grid(row=0, column=0, sticky='ew')
    frame.grid_columnconfigure(1, weight=1)

    tipText = 'Which set of restraints is currently active in the popup;\n' + \
              'a grouping for restraint and violation lists'
    label = Label(frame, text ='Restraint Set:', grid=(0,0))

    self.constraintSetPulldown = PulldownList(frame, tipText=tipText, grid=(0,1),
                                              callback=self.changeConstraintSet)

    utilButtons = UtilityButtonList(frame, grid=(0,2), sticky='e', 
                                    helpUrl=self.help_url )

    options = ['Restraint Sets','Restraint Lists',
               'Restraints','Violations Lists',
               'Violations']
    self.tabbedFrame = TabbedFrame(guiFrame, options=options, grid=(1,0),
                                   callback=self.updateTabbedFrame)
    frameA, frameB, frameC, frameD, frameE = self.tabbedFrame.frames
    
    
    # Restraint Sets
    
    frameA.grid_columnconfigure(0, weight=1)
    frameA.grid_rowconfigure(0, weight=1)
    
    tipTexts = ['Serial number of restraint set',
                'Number of restraint lists, of all types, in the set',
                'Number of violation lists in the set',
                'Number of resonances (fixed in assignment) used to describe restraints in the set',
                'Number of atom groups (fix in assignment) used in the set\'s restraints',
                'Number of fixed resonance to atom set assignments',
                'Number of alternate chain conformations explicitly covered']
    
    colHeadings = ['#','Restraint\nLists','Violation\nLists','Fixed\nResonances',
                   'Fixed\nAtom Sets','Fixed\nResonance Sets','Chain\nStates']

    self.constraintSetMatrix = ScrolledMatrix(frameA, multiSelect=True, grid=(0,0),
                                              headingList=colHeadings, tipTexts=tipTexts,
                                              callback=self.selectConstraintSetCell,
                                              deleteFunc=self.deleteRestraintSet)
                                              

    tipTexts = ['**Disabled option**, import restraints from file. Use FormatConverter',
                'Delete the selected restraint set with all its restraints and violations']
    texts    = ['Import','Delete']
    commands = [None,self.deleteRestraintSet]
    self.constraintSetButtons = ButtonList(frameA,texts=texts, tipTexts=tipTexts,
                                           expands=True, commands=commands, grid=(1,0))
    
    
    # Restraint Lists
    frameB.grid_columnconfigure(0, weight=1)
    frameB.grid_rowconfigure(0, weight=1)
    
    self.nameListEntry    = Entry(self,text='', returnCallback=self.setListName, width=10)
    self.detailsListEntry = Entry(self,text='', returnCallback=self.setListDetails, width=12)

    tipTexts = ['The restraint list serial number, within its containing restraint set',
                'The restraint list type; distance, dihedral, h-bond etc.',
                'The human name to identify the restraint list',
                'The number of individual restraints in the list',
                'The NMR experiments used to derive the restraints',
                'User editable textual comment for the list',
                'The units of the restrained quantity']
    colHeadings = ['#','Type','Name','Restraints','Experiments','Details','Unit']
    justifyList = ['center', 'center', 'left', 'center', 'left', 'left', 'center']
    editWidgets      = [None, None, self.nameListEntry, None, None, self.detailsListEntry, None]
    editGetCallbacks = [None, None, self.getListName,   None, None, self.getListDetails,   None]
    editSetCallbacks = [None, None, self.setListName,   None, None, self.setListDetails,   None]
    self.constraintListMatrix = ScrolledMatrix(frameB,  multiSelect=True,
                                           editSetCallbacks=editSetCallbacks,
                                           editGetCallbacks=editGetCallbacks,
                                           editWidgets=editWidgets,initialRows=10,
                                           headingList=colHeadings, grid=(0,0),
                                           justifyList=justifyList, tipTexts=tipTexts,
                                           callback=self.selectConstraintListCell,
                                           deleteFunc=self.deleteRestraintList)
    
    tipTexts = ['Show a table of individual restraints for the selected restraint list',
                'Export the selected restraint list (if distance type) as a CNS .tbl file for ARIA v1',
                'Make peaks equivalent to restraints. Requires that the restraints were initially made in CCPN',
                'Delete the selected restraint lists']                                       
    texts    = ['Show\nRestraints','Export\nARIA list',
                'Create\nEquiv. Peaks','Delete\nSelected Lists']
    commands = [self.viewConstraints,self.exportAriaList,
                self.makeRestraintPeaks,self.deleteRestraintList]
    self.restraintListButtons = ButtonList(frameB, texts=texts, grid=(1,0),
                                           tipTexts=tipTexts, commands=commands)
    #self.restraintListButtons.buttons[1].disable()

    tipTexts = ['Merge the selected restraint lists into a single list (if the same type)',
                'Split the last selected restraint list into two new lists; one which has ambiguous (excluding prochiral ambiguity) atom restraints and one which does not',
                'Split the last selected restraint list into two new lists; one with violated restraints, one with non-violated']
    texts    = ['Merge Lists','Split Ambig/Unambig','Split Violated']
    commands = [self.mergeRestraintLists,
                self.splitRestraintListAmbig,
                self.splitRestraintListViol]
    self.restraintListButtons2 = ButtonList(frameB, texts=texts, tipTexts=tipTexts, 
                                            commands=commands, grid=(2,0))

    # Restraints
    
    frameC.grid_columnconfigure(7, weight=1)

    tipText = 'Select the specific list to display restraints for; number:type:name'
    label = Label(frameC, text='Restraint List:', grid=(0,0), sticky='e')

    self.constraintListPulldown = PulldownList(frameC, grid=(0,1), tipText=tipText,
                                               callback=self.changeRestraintList)

    tipText = 'Select which violation list, if any, to display restraint violations from'
    label = Label(frameC, text=' Violation List:', grid=(0,2))
    self.violationListPulldownA = PulldownList(frameC, tipText=tipText, grid=(0,3),
                                               callback=self.changeViolationListA)
    
    tipText = 'Select which loaded structure to use for restraint value comparison'
    label = Label(frameC, text=' Structure:', grid=(0,4))
    self.structurePulldown = PulldownList(frameC, callback=self.changeStructure,
                                          grid=(0,5), tipText=tipText)
    
    tipText = 'Select how values are generated from a structure: "NOE sum" is used where\n' + \
    'ambiguous distance restraint options all contribute NOE peak intensity'
    label = Label(frameC, text=' Value Method:', grid=(0,6))
    self.valueMethodPulldown = PulldownList(frameC, callback=self.setValueMethod,
                                            grid=(0,7), tipText=tipText)
    
    
    frameC.grid_rowconfigure(1, weight=1)

    self.detailsRestraintEntry = Entry(self, text='', width=12,
                                       returnCallback=self.setRestraintDetails)
    
    tipTexts = ['Serial number of restraint & any ambiguous item (restraint:item)',
                'The name of the resonances restrained (fixed within the restraint set)',
                'The target or best value for the restraint, e.g. optimum distance',
                'The lower bound value for the restraint',
                'The upper bound value for the restraint',
                'An error value for the restraint, often just the with between bounds',
                'The number of CCPN peaks used to derive the restraint information',
                'The average value of a violation over a structure ensemble',
                'The fraction of a structure ensemble\'s conformations (models) violated by the restraint',
                'The value of the data being restrained in the selected structure ensemble',
                'The data value used to derive the restraint e.g. relative peak intensity',
                'The relative strength weighing of the restraint (not used by all structure calculation protocols) ',
                'User-editable textual comment for the restraint',
                'The number of ambiguous or alternative possibilities on the restraint']                                       
    colHeadings = ['#','Resonances','Value','Upper\nLimit','Lower\nLimit',
                   'Error','Peaks','Mean\nViol','Viol\nFraction','Struc\nValue','Orig.\nData',
                   'Weight','Details','Items']
    editWidgets      = [None] * len(colHeadings)
    editGetCallbacks = [None] * len(colHeadings)
    editSetCallbacks = [None] * len(colHeadings)
    
    editWidgets[12] = self.detailsRestraintEntry
    editGetCallbacks[12] = self.getRestraintDetails
    editSetCallbacks[12] = self.setRestraintDetails
    
    self.constraintsMatrix = ScrolledMatrix(frameC, multiSelect=True,
                                            editSetCallbacks=editSetCallbacks,
                                            editGetCallbacks=editGetCallbacks,
                                            editWidgets=editWidgets,
                                            headingList=colHeadings,
                                            tipTexts=tipTexts,
                                            callback=self.selectRestraintCell,
                                            deleteFunc=self.deleteConstraint)
    self.constraintsMatrix.grid(row=1, column=0, columnspan=8, sticky='nsew')
    self.constraintsMatrix.doEditMarkExtraRules = self.doRestraintEditMarkExtraRules
    
    tipTexts = ['Delete the selected restraint items; assignment possibility not necessarily whole restraint',
                'Delete all (whole) restraints for selected items, even where only one item of an ambiguous restraint is selected',
                'If present in the table, highlight the rows of the restraints derived from peaks selected in spectrum windows',
                'Update the items (assignment possibilities) for a restraint from any linked CCPN peak',
                'Show a table of peaks from which the selected restraints were derived',
                'Show the connections of the selected restraints on a graphical structure display',
                'Assign the peak from which a restraint is derived']                                       
    texts    = ['Delete\nItems','Delete\nRestraints',
                'Show Restraints for\nSelected Peaks',
                'Update Assignment\nFrom Peak',
                'Show\nPeaks','Show Selected\nOn Structure',
                'Assign\nPeak']
    commands = [self.deleteItem, self.deleteConstraint,
                self.viewForPeakSelection,
                self.updateAssignFromPeak,
                self.viewPeaks, self.showStructConnections,
                self.assignRestraintPeak]
    self.restraintButtons = ButtonList(frameC, tipTexts=tipTexts,
                                       texts=texts, commands=commands,
                                       grid=(2,0), gridSpan=(1,8))
    
    tipTexts = ['Show a table of individual violation records for selected list',
                'Calculate restraint violation records for restraints using selected structure',
                'Export the restraints as a CNS .tbl file for ARIA v1',
                'Delete the entire restraint list currently on display',
                'For imported restraints add any possible CCPN resonance links',
                'Merge any restraints that represent the same atoms pairings together, e.g. from reciprocating return NOE peaks',
                'Update the assignments of any linked CCPN peaks according to the current restraint items']                                       
    texts    = ['Show\nViolations','Calculate\nViolations',
                'Export\nARIA List','Delete\nList',
                'Setup\nResonances','Merge\nDuplicates',
                'Update Peak\nAssignments',]
    commands = [self.viewViolationList, self.calcViolations,
                self.exportAriaList,self.deleteRestraintListSingle,
                self.setupResonances, self.mergeDuplicates,
                self.updatePeakAssign]
                
    self.restraintButtons2 = ButtonList(frameC,commands=commands, grid=(3,0),
                                        texts=texts, tipTexts=tipTexts, gridSpan=(1,8))
    

    # Violation Lists
    frameD.grid_columnconfigure(0, weight=1)
    frameD.grid_rowconfigure(0, weight=1)
    
    
    tipTexts = ['The serial number of the violation list',
                'The number of violation records in the list',
                'The number of structural models used in violation calculation',
                'User-editable textual comment for violation list']                                       
    colHeadings = ['#','Violations','Structures','Details',]
    self.violationListMatrix = ScrolledMatrix(frameD, multiSelect=True,
                                         headingList=colHeadings, grid=(0,0),
                                         callback=self.selectViolationListCell,
                                         deleteFunc=self.deleteViolationList,
                                         tipTexts=tipTexts)
                                         
    tipTexts = ['Show a table of individual violation entries for the selected list',
                'Delete the selected violation lists']                                       
    texts    = ['Show Violations','Delete Lists']
    commands = [self.viewViolations,self.deleteViolationList]
    self.violationListButtons = ButtonList(frameD, texts=texts, tipTexts=tipTexts,
                                           grid=(1,0), commands=commands)

    div = LabelDivider(frameD, text='Structures Analysed', grid=(2,0), gridSpan=(1,4))
     
    tipTexts = ['Number of structure ensemble used in violation analysis',
                'The number of the conformational models used in the violation analysis',
                'The molecular chains present in the structure',
                'Any structure generation group to which the structure belongs'] 
    colHeadings = ['#','Models','Chains','Generation']
    self.structuresMatrix = ScrolledMatrix(frameD, headingList=colHeadings, tipTexts=tipTexts)
    self.structuresMatrix.grid(row=3, column=0, sticky='nsew')

    # Violations
    
    frameE.grid_columnconfigure(0, weight=1)
    frameE.grid_rowconfigure(1, weight=1)
    frameE1 = Frame(frameE)
    frameE1.grid(row=0, column=0, sticky='ew')
    frameE1.grid_columnconfigure(1, weight=1)
    
    tipText = 'Select which violation list to show violation items for'
    label = Label(frameE1, text = 'Violation List:', grid=(0,0), sticky='e')
    self.violationListPulldownB  = PulldownList(frameE1, grid=(0,1), tipText=tipText, 
                                                callback=self.changeViolationListB)
    
    tipTexts = ['The restraint list that the violations were calculated for',
                'The restraint that was analysed in terms of value & bounds',
                'The amount violated; difference from the nearest restraint bound',
                'The % of conformational models violated in structure',
                'The type of the restraint analysed; distance, dihedral etc.',
                'The number of peaks used to derive the constraint',
                'The structure-derived value use to calculate the violation',
                'The error in the violation amount']                                       
    colHeadings = ['Restraints\nList','Restraint','Amount',
                   '% Violated','Type','Peaks',
                   'Calc\nValue','Value\nError']
    self.violationsMatrix = ScrolledMatrix(frameE, tipTexts=tipTexts,
                                           headingList=colHeadings,
                                           callback=self.selectViolationCell)
    self.violationsMatrix.grid(row=1, column=0, sticky='nsew')

    tipTexts = ['Assign any peak from which the violated restraint may be derived',                                       
                'Show a table of peaks thare were used to create the restraint that was violated',
                'Open a table of restraints and highlight the one corresponding to the selected violation']
    texts    = ['Assign Peak', 'Show Peaks', 'Show Restraint']
    commands = [self.assignViolationPeak,
                self.showViolationPeaks,
                self.showViolationRestraint]
    self.violationButtons = ButtonList(frameE, expands=True,
                                       texts=texts, commands=commands)
    self.violationButtons.grid(row=2, column=0, sticky='ew')
     
    # Main
    
    self.updateConstraintSetsAfter()

    self.updateButtons()

    self.administerNotifiers(self.registerNotify)
  
  
  def doRestraintEditMarkExtraRules(self, item, row, col):
    
    if str(self.constraintsMatrix.textMatrix[row][0])[-2:] == ':0':
      #if item is item.constraint.sortedItems()[0]:
      return True
    
    return False  
  
  def administerNotifiers(self, notifyFunc):

    for func in ('__init__', 'delete','setChainStates','addChainState','removeChainState'):
        notifyFunc(self.updateConstraintSetsAfter, 'ccp.nmr.NmrConstraint.NmrConstraintStore', func)

    for func in ('__init__', 'delete','setName','setDetails','setUnit',
                 'setExperiments','addExperiment','removeExperiment'):
      for clazz in ('ccp.nmr.NmrConstraint.ChemShiftConstraintList',
                    'ccp.nmr.NmrConstraint.DihedralConstraintList',
                    'ccp.nmr.NmrConstraint.DistanceConstraintList',
                    'ccp.nmr.NmrConstraint.HBondConstraintList',
                    'ccp.nmr.NmrConstraint.JCouplingConstraintList',
                    'ccp.nmr.NmrConstraint.RdcConstraintList'):
        notifyFunc(self.updateConstraintListsAfter, clazz, func)

    for func in ('__init__', 'delete','setDetails'):   
      for clazz in ('ccp.nmr.NmrConstraint.ViolationList',):
        notifyFunc(self.updateViolationListsAfter, clazz, func)
             
    for func in ('__init__', 'delete','setError',
                 'setDetails','setOrigData','setWeight',
                 'setLowerLimit','setTargetValue','setUpperLimit'):
      for clazz in ('ccp.nmr.NmrConstraint.ChemShiftConstraint',
                    'ccp.nmr.NmrConstraint.DistanceConstraint',
                    'ccp.nmr.NmrConstraint.HBondConstraint',
                    'ccp.nmr.NmrConstraint.JCouplingConstraint',
                    'ccp.nmr.NmrConstraint.RdcConstraint'):
        notifyFunc(self.updateRestraintsAfter, clazz, func)

    for func in ('__init__', 'delete','setError',
                 'setLowerLimit','setTargetValue','setUpperLimit'):
       notifyFunc(self.updateItemsAfter, 'ccp.nmr.NmrConstraint.DihedralConstraintItem', func)
 
    for func in ('__init__', 'delete','setError','setValue'):
      for clazz in ('ccp.nmr.NmrConstraint.Violation',):
        notifyFunc(self.updateViolationsAfter, clazz, func)
        
    for func in ('__init__', 'delete'):
      notifyFunc(self.updateResonanceSetAfter,'ccp.nmr.NmrConstraint.FixedResonanceSet', '')

    for func in ('__init__', 'delete','setResonances'):
      notifyFunc(self.updateRestraintsAfter, 'ccp.nmr.NmrConstraint.DihedralConstraint', func)
        
    for func in ('__init__', 'delete','setResonances'):
      for clazz in ('ccp.nmr.NmrConstraint.DistanceConstraintItem',
                    'ccp.nmr.NmrConstraint.HBondConstraintItem',
                    'ccp.nmr.NmrConstraint.JCouplingConstraintItem',
                    'ccp.nmr.NmrConstraint.RdcConstraintItem'):
        notifyFunc(self.updateItemsAfter, clazz, func)

    notifyFunc(self.updateRestraintsAfter,'ccp.nmr.NmrConstraint.ChemShiftConstraint', 'setResonance')
    notifyFunc(self.updateResiduesAfter,'ccp.molecule.MolSystem.Residue', 'setSeqCode')
    notifyFunc(self.updateResiduesAfter,'ccp.molecule.MolSystem.Residue', 'setSeqCode')
    #notifyFunc(self.updateAfter,'ccp.nmr.NmrConstraint.ChemShiftConstraintList', 'setIsotopeCode')
    
    for func in ('delete','__init__'):
      notifyFunc(self.updateStructures, 'ccp.molecule.MolStructure.StructureEnsemble', func)
      notifyFunc(self.updateStructureModels, 'ccp.molecule.MolStructure.Model', func)

  def open(self):
  
    self.updateConstraintSetsAfter()
    self.updateConstraintListsAfter()
    self.updateViolationListsAfter()
    self.updateStructures()
    self.updateButtons()

    BasePopup.open(self)
  
      
  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)
  
  def updateTabbedFrame(self, index):
  
    # Updates when viewed
    funcDict = {0:self.updateConstraintSetsAfter,
                1:self.updateConstraintListsAfter,
                2:self.updateRestraintsAfter,
                3:self.updateViolationListsAfter,
                4:self.updateViolationsAfter}

    funcDict[index]()
    
    self.updateButtons()

  def update(self, obj):
  
    self.constraintSet = obj.nmrConstraintStore
    self.updateConstraintSetsAfter()
    
    if obj.className == 'ViolationList':
      self.tabbedFrame.select(3)
    else:
      self.tabbedFrame.select(1)
      
 
  # Update functions - pulldowns and tables

  def updateResonanceSetAfter(self, fixedResonanceSet):
  
    if self.tabbedFrame.selected != 2:
      return
      
    constraintListB = self.constraintListB
    
    if constraintListB and fixedResonanceSet.nmrConstraintStore is self.constraintSet:
      constraintType = self.constraintListB.className[:-14].lower()
      
      if constraintType == 'dihedral':
        for resonance in fixedResonanceSet.resonances:
          for constraint in resonance.dihedralConstraints:
            if constraint.constraintList is constraintListB:
              self.updateRestraintsAfter()
              return

      else:
        attribute = constraintType + 'ConstrItems'
      
        for resonance in fixedResonanceSet.resonances:
          for item in resonance.__dict__[attribute]:
            if item.constraint.constraintList is constraintListB:
              self.updateRestraintsAfter()
              return      
  
  def updateConstraintSetsAfter(self, obj=None):
  
    if self.waiting:
      return
      
    self.waiting=True

    self.after_idle(self.updateConstraintSets)  

  def updateConstraintSets(self):

    # Table
  
    if self.tabbedFrame.selected == 0:
      objectList  = []
      textMatrix  = []
 
      for constraintSet in self.nmrProject.sortedNmrConstraintStores():
        datum = [constraintSet.serial,
                 len(constraintSet.constraintLists),
                 len(constraintSet.violationLists),
                 len(constraintSet.fixedResonances),
                 len(constraintSet.fixedAtomSets),
                 len(constraintSet.fixedResonanceSets),
                 len(constraintSet.chainStates)]
 
        textMatrix.append(datum)
        objectList.append(constraintSet)
 
      self.constraintSetMatrix.update(objectList=objectList,
                                      textMatrix=textMatrix)
 
  
    # Pulldowns
  
    index = 0
    names = []
    constraintSet  = self.constraintSet
  
    constraintSets = self.nmrProject.sortedNmrConstraintStores()
    
    if constraintSets:
      if constraintSet not in constraintSets:
        constraintSet = constraintSets[0]       
          
      index =  constraintSets.index(constraintSet)
      names = ['%d' % cs.serial for cs in constraintSets]

    else:
      constraintSet = None
    
    if self.constraintSet is not constraintSet:
      self.constraintSet = constraintSet
      self.updateConstraintListsAfter()
      self.updateViolationListsAfter()
      self.updateStructures()
    
    self.constraintSetPulldown.setup(names, constraintSets, index)
    
    self.waiting = False  
  
  def updateConstraintListsAfter(self, restraintList=None):
     
    if self.tabbedFrame.selected in (0,1,2):
      if restraintList:
        constraintSet = restraintList.parent
        if constraintSet is self.constraintSet:
          self.updateConstraintLists()
      else:
        self.updateConstraintLists()
    
    if self.tabbedFrame.selected == 1:
      if not self.waitingConstraintLists:
        self.waitingConstraintLists = True
        self.after_idle(self.updateConstraintListTable)
 
  def getConstraintLists(self):
  
    constraintLists = []
    if self.constraintSet:
      for constraintList in self.constraintSet.sortedConstraintLists():
        constraintLists.append(constraintList)
    
    return constraintLists  
    
  def updateConstraintLists(self):
    
    index = 0
    names = []
    constraintList = self.constraintListB
    constraintLists = self.getConstraintLists()
    getName = self.getConstraintListName
    
    if constraintLists:
      if constraintList not in constraintLists:
        constraintList = constraintLists[0]
      
      names = [getName(cl) for cl in constraintLists]  
      index = constraintLists.index(constraintList)
      #label = constraintList.className[:-14]
      
    else:
      constraintList = None

    if constraintList is not self.constraintListB:
      self.constraintListB = constraintList 
      self.updateRestraintsAfter()
      self.updateValueMethods()

    self.constraintListPulldown.setup(names, constraintLists, index)

    self.waitingConstraintLists = False
 
  def updateConstraintListTable(self, *opt):  
    
    constraintSet = self.constraintSet
    
    objectList  = []
    textMatrix  = []
    
    if constraintSet:
      for constrList in constraintSet.sortedConstraintLists():
        datum = [constrList.serial,
                 constrList.className[0:-14],
                 constrList.name or ' ',
                 len(constrList.constraints),
                 ' '.join([e.name for e in constrList.experiments]),
                 constrList.details,
                 constrList.unit]
 
        textMatrix.append(datum)
        objectList.append(constrList)
     
    self.constraintListMatrix.update(objectList=objectList,
                                     textMatrix=textMatrix)
 
  def updateViolationListsAfter(self, violationList=None):
    
    if self.tabbedFrame.selected == 3:
      if not self.waitingViolationLists:
        self.waitingViolationLists = True
        self.after_idle(self.updateViolationListTable)
    
    if self.tabbedFrame.selected in (2,3,4):
      if violationList:
        if violationList.parent is self.constraintSet:
          self.updateViolationLists()
      else:
        self.updateViolationLists()
       

  def updateViolationLists(self):
  
    indexA = 0
    indexC = 0
    names  = []
    violationListsA = []
    violationListsC = []
    
    if self.constraintSet:
      violationListsC = self.constraintSet.sortedViolationLists()
      violationListsA = violationListsC + [None,]
      violationListA = self.violationListA 
      violationListC = self.violationListC 
 
      if violationListsC:
 
        if violationListA not in violationListsA:
          violationListA = violationListsA[-1]
          
        if violationListC not in violationListsC:
          violationListC = violationListsC[0]
 
        names = ['%d' % vl.serial for vl in violationListsC]
        indexA = violationListsA.index(violationListA)
        indexC = violationListsC.index(violationListC)
        
    else:
      violationListA = None
      violationListC = None
    
    if violationListA is not self.violationListA:
      self.violationListA = violationListA
      self.updateRestraintsAfter()
    
    if violationListC is not self.violationListC:
      self.violationListC = violationListC
      self.updateViolationsAfter()
    
    self.violationListPulldownB.setup(names, violationListsC, indexC)
    
    names.append('<New>')
    self.violationListPulldownA.setup(names, violationListsA, indexA)


  def updateViolationListTable(self):  
    
    constraintSet = self.constraintSet
    
    objectList  = []
    textMatrix  = []
    
    if constraintSet:
      for violationList in constraintSet.sortedViolationLists():
        strucText = None
        models = violationList.molStructures
        
        if models:
          structDict = {}
          
          for model in models:
            structure = model.structureEnsemble
            structDict[structure] = True
          
          data = ['%s:%s' % (s.molSystem.code,s.ensembleId) for s in structDict]
          strucText = ','.join(data)
         
        datum = [violationList.serial,
                 len(violationList.violations),
                 strucText,
                 violationList.details]
                 
        objectList.append(violationList)
        textMatrix.append(datum)
 
    self.violationListMatrix.update(objectList=objectList,
                               textMatrix=textMatrix)
    
    self.waitingViolationLists = False

  def updateItemsAfter(self, item):

   self.updateRestraintsAfter(item.constraint)

  def updateResiduesAfter(self, residue):
    
    selected = self.tabbedFrame.selected
        
    if selected == 2:
      if not self.waitingRestraints:
         self.waitingRestraints = True
         self.after_idle(self.updateRestraints)
    
    elif selected == 4:
      if not self.waitingViolations:
        self.waitingViolations = True
        self.after_idle(self.updateViolations)

  
  def updateRestraintsAfter(self, restraint=None):
     
    if self.tabbedFrame.selected == 1:
      if restraint and restraint.parent.parent is self.constraintSet:
        if not self.waitingConstraintLists:
          self.waitingConstraintLists = True
          self.after_idle(self.updateConstraintListTable)
        
    if self.tabbedFrame.selected == 2:
      if restraint:
        if restraint.parent is self.constraintListB:
          if not self.waitingRestraints:
            self.waitingRestraints = True
            self.after_idle(self.updateRestraints)
      else:
        self.waitingRestraints = True
        self.after_idle(self.updateRestraints)
      
      self.updateValueMethods()
      self.updateConstraintListsAfter()
      self.updateViolationListsAfter()
      self.updateStructures()

  def updateRestraints(self):  
    
    objectList  = []
    textMatrix  = []
    colorMatrix = []
    constraints = []
    constrType = None
    
    # For speed
    objectListAppend = objectList.append
    textMatrixAppend = textMatrix.append
    colorMatrixAppend = colorMatrix.append
    
    if self.constraintListB:
      constraints = self.constraintListB.sortedConstraints()
      constrType  = self.constraintListB.className[:-14]
    
    nCols = len(self.constraintsMatrix.headingList)
    colors0 = ['#A0D0A0'] * nCols
    colors0[0] = '#80A080'
    colorsA = [VIOL_COLOR_A] * nCols
    colorsB = [VIOL_COLOR_B] * nCols
    colorsC = [VIOL_COLOR_C] * nCols
    colorsD = [VIOL_COLOR_D] * nCols
    blank = [None] * nCols
    getResonanceName = self.makeResonanceGuiName
    violListA = self.violationListA
    struct = self.structure
    distMethod = self.distMethod
    dashJoin = '-'.join
    
    if struct:
      quickStrucDict = self.quickStrucDict.get(struct, {})
      self.quickStrucDict[struct] = quickStrucDict
      quickStrucDictGet = quickStrucDict.get
      
    # These conditions are unrolled for speed
    
    if constrType == 'Dihedral':
      
      for constraint in constraints:
        cPeaks = [cpc.peak for cpc in constraint.peakContribs]
        nPeaks = len(cPeaks)
        
        if None in cPeaks:
          if nPeaks == cPeaks.count(None):
            nPeaks = '%d (deleted)' % nPeaks
         
          else:
            nPeaks = '%d (%d deleted)' % (nPeaks, cPeaks.count(None))
        
        items  = constraint.sortedItems()
        nItems = len(items)
        cSerial = constraint.serial
        origdat = constraint.origData
        weight  = constraint.weight
 
        value      = None
        upperLimit = None
        lowerLimit = None
        error      = None

        resonances = constraint.resonances
        resonanceNames = [getResonanceName(r) for r in resonances]
        resonanceText = dashJoin(resonanceNames)
  
        strucVal = None
        if struct:
          atomSets = []
          for resonance in resonances:
            resonanceSet = resonance.resonanceSet

            if resonanceSet:
              atomSets.append(resonanceSet.atomSets)

          if len(atomSets) == 4:
            atomSets = tuple(atomSets)
            strucVal = quickStrucDictGet(atomSets) or \
                       getAtomSetsDihedral(atomSets, struct)
            quickStrucDict[atomSets] = strucVal
            
        viol   = None
        violFrac = None
        for violation in constraint.violations:
          if violation.violationList is violListA:
            viol = violation.violation
            violFrac = violation.fractionViolated
            break
 
        for j, item in enumerate(items):
          objectListAppend(item)
 
          targetValue = item.targetValue
          if targetValue is None:
            value = None
          else:
            value = '%5.3f' % targetValue
 
          upperLimit = '%5.3f' % item.upperLimit
          lowerLimit = '%5.3f' % item.lowerLimit
 
          error = item.error
          if error is not None:
            error = '%5.3f' % error
 
          if j == 0:
            serial  = '%d:0' % cSerial
            details = constraint.details
 
            if viol is not None:
              if violFrac <= 0.1:
                colorMatrixAppend(colorsA)
              elif violFrac <= 0.3:
                colorMatrixAppend(colorsB)
              elif violFrac <= 0.6:
                colorMatrixAppend(colorsC)
              else:
                colorMatrixAppend(colorsD)

            else:
              colorMatrixAppend(colors0)
 
          else:
            serial  = '%d:%2.2d' % (cSerial, j)
            details = ''
            colorMatrixAppend(blank)
 
          datum = [cSerial,
                   resonanceText,
                   value,
                   upperLimit,
                   lowerLimit,
                   error,
                   nPeaks,
                   viol,
                   violFrac,
                   strucVal,
                   origdat,
                   weight,
                   details,
                   nItems]
 
          textMatrixAppend(datum)    
    
    elif constrType == 'ChemShift':
    
      for constraint in constraints:
        objectListAppend(constraint)
        cPeaks = [cpc.peak for cpc in constraint.peakContribs]
        nPeaks = len(cPeaks)
        
        if None in cPeaks:
          if nPeaks == cPeaks.count(None):
            nPeaks = '%d (deleted)' % nPeaks
         
          else:
            nPeaks = '%d (%d deleted)' % (nPeaks, cPeaks.count(None))
        
        strucVal = None
        details = constraint.details
        nItems = None
        cSerial = constraint.serial
        origdat = constraint.origData
        weight  = constraint.weight
 
        value      = constraint.targetValue
        upperLimit = constraint.upperLimit
        lowerLimit = constraint.lowerLimit
        error      = constraint.error
        
        resonance = constraint.resonance
        resonanceText = getResonanceName(resonance)
        
        viol = None
        violFrac = None
        for violation in constraint.violations:
          if violation.violationList is violListA:
            viol = violation.violation
            violFrac = violation.fractionViolated
            break

        if viol is not None:
          if violFrac <= 0.1:
            colorMatrixAppend(colorsA)
          elif violFrac <= 0.3:
            colorMatrixAppend(colorsB)
          elif violFrac <= 0.6:
            colorMatrixAppend(colorsC)
          else:
            colorMatrixAppend(colorsD)

        else:
          colorMatrixAppend(colors0)
        
        datum = [cSerial,
                 resonanceText,
                 value,
                 upperLimit,
                 lowerLimit,
                 error,
                 nPeaks,
                 viol,
                 violFrac,
                 strucVal,
                 origdat,
                 weight,
                 details,
                 nItems]

        textMatrixAppend(datum)
     
    else:
 
      constraintItems = []
      for constraint in constraints:
        viol   = None
        violFrac = None
        
        items = constraint.items
        cSerial = constraint.serial
        origdat = constraint.origData
        weight  = constraint.weight
        value = constraint.targetValue
        upperLimit = constraint.upperLimit
        lowerLimit = constraint.lowerLimit
        error = constraint.error
        
        cPeaks = [cpc.peak for cpc in constraint.peakContribs]
        nPeaks = len(cPeaks)
        
        if None in cPeaks:
          if nPeaks == cPeaks.count(None):
            nPeaks = '%d (deleted)' % nPeaks
         
          else:
            nPeaks = '%d (%d deleted)' % (nPeaks, cPeaks.count(None))
          
        nItems = len(items)
 
        for violation in constraint.violations:
          if violation.violationList is violListA:
            viol = violation.violation
            violFrac = violation.fractionViolated
            break
 
        ambigNum = {}
        ambigNumGet = ambigNum.get
        itemList = []
        itemListAppend = itemList.append
 
        for item in items:

          resonances = item.resonances
          for resonance in resonances:
            ambigNum[resonance] = ambigNumGet(resonance, 0) + 1
 
          strucVal = None
          if struct:
            atomSets = [None, None]
            for i, resonance in enumerate(resonances):
              assignedSets = []
              resonanceSet = resonance.resonanceSet
 
              if resonanceSet:
                atomSets[i] = resonanceSet.atomSets
              else:
                break  
 
            else:
              atomSets = tuple(atomSets)
              strucVal = quickStrucDictGet(atomSets) or \
                         getAtomSetsDistance(atomSets[0], atomSets[1],
                                             struct, method=distMethod)
              quickStrucDict[atomSets] = strucVal
              
          datum = [resonances,
                   None,
                   value,
                   upperLimit,
                   lowerLimit,
                   error,
                   nPeaks,
                   viol,
                   violFrac,
                   strucVal,
                   origdat,
                   weight,
                   None,
                   nItems,
                   item]
 
          itemListAppend(datum)
 
        for datum in itemList:
          resonances = [(ambigNumGet(r, 0), r) for r in datum[0]]
          resonances.sort()
          resonances.reverse()
          resonanceNames = [getResonanceName(r[1]) for r in resonances]
          resonanceText = dashJoin(resonanceNames)
          datum[0] = resonanceText
          datum[1] = resonanceText
 
        itemList.sort()
        for j, datum in enumerate(itemList):
          if j == 0:
            serial = '%d:0' % cSerial
            details = constraint.details
 
            if viol is not None:
              if violFrac <= 0.1:
                colorMatrixAppend(colorsA)
              elif violFrac <= 0.3:
                colorMatrixAppend(colorsB)
              elif violFrac <= 0.6:
                colorMatrixAppend(colorsC)
              else:
                colorMatrixAppend(colorsD)

            else:
              colorMatrixAppend(colors0)
 
          else:
            serial = '%d:%2.2d' % (cSerial, j)
            details = ''
            colorMatrixAppend(blank)
 
          datum[0]  = serial
          datum[12] = details
          textMatrixAppend(datum[:-1])
          objectListAppend(datum[-1])
        
    self.constraintsMatrix.update(objectList=objectList,
                                  textMatrix=textMatrix,
                                  colorMatrix=colorMatrix)
    
    self.waitingRestraints = False
     
  
  def updateItemAfter(self, item=None):

    if self.tabbedFrame.selected != 2:
      return
    
    if item and item.constraint.parent is not self.constraintListB:
      return
    
    if not self.waitingRestraints:
      self.waitingRestraints = True
      self.after_idle(self.updateRestraints)

  def updateViolationsAfter(self, violation=None):
  
    if violation:
      if violation.violationList is self.violationListA:
        self.updateRestraintsAfter()
          
        if self.tabbedFrame.selected == 4:
          if violation.violationList is self.violationListC:
            if not self.waitingViolations:
              self.waitingViolations = True
              self.after_idle(self.updateViolations)

    elif self.tabbedFrame.selected == 4:
      if not self.waitingViolations:
        self.waitingViolations = True
        self.after_idle(self.updateViolations)

  def updateViolations(self):  
    
    objectList  = []
    textMatrix  = []

    if self.violationListC:
      for violation in self.violationListC.sortedViolations():
        constraint = violation.constraint
        constrList = constraint.parentList
        constrType = constrList.className[:-14]
 
        datum   = [constrList.serial,
                   constraint.serial,
                   violation.violation,
                   violation.fractionViolated * 100.0,
                   constrType,
                   len(constraint.peaks),
                   violation.calcValue,
                   violation.calcValueError]
 
        textMatrix.append(datum)
        objectList.append(violation)

        
    self.violationsMatrix.update(objectList=objectList,
                                 textMatrix=textMatrix)
    
    self.updateStructures()

    self.waitingViolations = False        
    
    
  def updateStructureModels(self, model):
  
    # Coords have changed, clear cache
    structure = model.structureEnsemble
    if self.quickStrucDict.get(structure):
      del self.quickStrucDict[structure]
  
    self.updateStructures()

  def updateStructures(self, structure=None):
  
    if structure and structure.isDeleted:
      if self.quickStrucDict.get(structure):
        del self.quickStrucDict[structure]
  
    # Pulldown
    
    if self.tabbedFrame.selected == 2:
      index = 0
      structures = [None,] + self.getStructures()

      if self.structure in structures:
        structure = self.structure
      else:
        structure = None
 
      names = ['<None>',] + [str(x.ensembleId) for x in structures[1:]]
      index = structures.index(structure)

      self.structurePulldown.setup(names, structures, index)

      if structure is not self.structure:
        self.structure = structure
    
    # Table
    if self.tabbedFrame.selected == 3:
      textMatrix = []
      objectList = []
 
      if self.violationListB:
        structDict = {}
        for model in self.violationListB.sortedMolStructures():
          structure = model.structureEnsemble
          if structDict.get(structure) is None:
            structDict[structure] = []
 
          structDict[structure].append(model)
 
        structures = [(s.ensembleId,s) for s in structDict.keys()]
 
        for eId, structure in structures:
          models = len(structDict[structure])
          chains = ' '.join([chain.code for chain in structure.coordChains])
          
          sg = structure.structureGeneration
                    
          if sg:
            group = '%d:%s' % (sg.serial,sg.name)
            
          else:
            group = None
            
          datum = [eId,models,chains,group]
 
          textMatrix.append(datum)
          objectList.append(structure)

      self.structuresMatrix.update(objectList=objectList, textMatrix=textMatrix)


  def updateValueMethods(self):

    if self.constraintListB and self.constraintListB.className == 'DistanceConstraintList':
      names = ['NOE sum','Minimum']
      if self.distMethod == 'noe':
        index = 0
      else:
        index = 1  
    else:
      names = ['Default',]
      index = 0
    
    self.valueMethodPulldown.setup(names, names, index)
    
  
  # Pulldown menu selection
  
      
  def changeConstraintSet(self, constraintSet):    
    
    if constraintSet is not self.constraintSet: 
      self.constraintSet = constraintSet
      self.violationListC = None
      self.violationListA = None
      self.constraintSetMatrix.hilightObject(constraintSet)
      self.updateConstraintListsAfter()
      self.updateViolationListsAfter()

  def changeViolationListA(self, violationList): 
     
    if violationList is not self.violationListA:
      self.violationListA = violationList
      
      if self.structure and violationList:
        if self.structure.models != violationList.molStructures:
          self.structure = None
      
      self.updateViolationListsAfter()
      self.updateRestraintsAfter()
    
  def changeViolationListB(self, violationList): 
    
    if violationList is not self.violationListC:
      self.violationListC = violationList
      self.updateViolationListsAfter()
      self.updateViolationsAfter()

  def changeRestraintList(self, constrList): 
    
    if constrList is not self.constraintListB:
      self.constraintListB = constrList
      self.updateValueMethods()
      self.updateRestraintsButtons()
      self.updateRestraintsAfter()

  def changeStructure(self, structure):

    if structure is not self.structure:
      self.structure = structure
      if structure and self.violationListA:
        if structure.models != self.violationListA.molStructures:
          self.violationListA = None
      
      else:
        self.violationListA = None
       
      self.violationListC = None
      self.updateRestraintsButtons()
      self.updateRestraintsAfter()

  def setValueMethod(self, name):

    if self.constraintListB.className == 'DistanceConstraintList':
      if name == 'NOE sum':
        distMethod = 'noe'
      else:
        distMethod = 'min'
     
      if distMethod != self.distMethod:
        if self.quickStrucDict.get(self.structure):
          del self.quickStrucDict[self.structure]

        self.distMethod = distMethod
        self.updateRestraintsAfter()

  # Table cell selection
  
  def selectViolationCell(self, violation, row, col):
  
    self.violation = violation
    self.updateViolationsButtons()

  def selectConstraintSetCell(self, constraintSet, row, col):
  
    self.constraintSet = constraintSet
    self.updateConstraintSets()
    self.updateConstraintSetButtons()

  def selectConstraintListCell(self, constrList, row, col):
  
    self.constraintListA = constrList
    self.constraintListB = constrList
    self.updateRestraintListButtons()
 
  def selectViolationListCell(self, violationList, row, col):
  
    self.violationListB = violationList
    self.updateStructures()
    self.updateViolationListButtons()
  
  def selectRestraintCell(self, item, row, col):
  
    self.item = item
    self.updateRestraintsButtons()
  
  def updateConstraintSetButtons(self):
  
    #    ['Import','Delete']
  
    buttons = self.constraintSetButtons.buttons
  
    buttons[0].disable()
    
    if self.constraintSet:
      buttons[1].enable()
    else:
      buttons[1].disable()
  
  def updateRestraintListButtons(self):
  
    #['View\nConstraints','Export\nARIA list',
    # 'Create\nEquiv. Peaks','Delete\nLists']
    # ['Merge Lists','Split Ambig/Unambig']
  
    buttons = self.restraintListButtons.buttons
    buttons2 = self.restraintListButtons2.buttons
  
    if self.constraintListA:
      buttons[0].enable()
      buttons[3].enable()
      buttons2[0].enable()
      buttons2[2].enable()
      if self.constraintListA.className == 'DistanceConstraintList':
        buttons[1].enable()
        buttons[2].enable()
        buttons2[1].enable()
    
      else:
        buttons[1].disable()
        buttons[2].disable()
        buttons2[1].disable()
    
    else:
      for button in buttons:
        button.disable()
      for button in buttons2:
        button.disable()
       
  
  def updateRestraintsButtons(self):
  
    # ['Delete\nItems','Delete\nConstraints',
    #  'Show Constraints for\nselected peaks',
    #  'Update Assignment\nFrom Peak',
    #  'View\nPeaks','View Selected\nOn Structure',
    #  'Assign\nPeak','Merge\nDuplicates']
    #  ['View\nViolations','Calculate\nViolations',
    #   'Export\nARIA List','Delete\nList','Setup\nResonances']
  
    buttons = self.restraintButtons.buttons
    buttons2 = self.restraintButtons2.buttons
  
    if self.constraintListB:
      buttons[2].enable()
      buttons2[3].enable()
      buttons2[4].enable()
      buttons2[5].enable()
    
      if self.constraintListB.className == 'DistanceConstraintList':
        buttons2[2].enable()
      else:  
        buttons2[2].disable()
        
      if self.violationListA:
        buttons2[0].enable()
      else:
        buttons2[0].disable()
      
      if self.structure:
        buttons2[1].enable()
      else:
        buttons2[1].disable()
      
    else:
      buttons2[0].disable()
      buttons2[1].disable()
      buttons2[2].disable()
      buttons2[3].disable()
      buttons2[4].disable()
      buttons2[5].disable()
      buttons[2].disable()
      
  
    if self.item:
      buttons[0].enable()
      buttons[1].enable()
      buttons2[6].enable()
      if self.constraintListB and \
         self.constraintListB.className == 'DistanceConstraintList':
        buttons[3].enable()
      else:  
        buttons[3].disable()

      if self.item.className == "ChemShiftConstraint":
        constraint = self.item
      else:
        constraint = self.item.constraint
      
      if constraint.peaks:
        buttons[4].enable()
        buttons[6].enable()
      else:  
        buttons[4].disable()
        buttons[6].disable()

      if self.structure:
        buttons[5].enable()
      else:
        buttons[5].disable()
      
    else:
      buttons[0].disable()
      buttons[1].disable()
      buttons[3].disable()
      buttons[4].disable()
      buttons[5].disable()
      buttons[6].disable()
      buttons2[6].disable()
        
  
  def updateViolationListButtons(self):
  
    buttons = self.violationListButtons.buttons
  
    if self.violationListB:
      buttons[0].enable()
      buttons[1].enable()
    else:
      buttons[0].disable()
      buttons[1].disable()
    
  
  def updateViolationsButtons(self):
  
    buttons = self.violationButtons.buttons
    
    if self.violation:
      buttons[0].enable()
      buttons[1].enable()
      buttons[2].enable()
    else:
      buttons[0].disable()
      buttons[1].disable()
      buttons[2].disable()
  
  def updateButtons(self):
  
    self.updateConstraintSetButtons()
    self.updateRestraintListButtons()
    self.updateRestraintsButtons()
    self.updateViolationListButtons() 
    self.updateViolationsButtons()
    
  # Table edit functions
    
  def getListName(self, constraintList):

    if constraintList :
      self.nameListEntry.set(constraintList.name)

  def setListName(self, event):

    text = self.nameListEntry.get()
    if self.constraintListA and text and text != ' ':
      self.constraintListA.setName( text )

  def getListDetails(self, constraintList):

    if constraintList :
      self.detailsListEntry.set(constraintList.details)

  def setListDetails(self, event):

    text = self.detailsListEntry.get()
    if self.constraintListA and text and text != ' ':
      self.constraintListA.setDetails(text)

  def getRestraintDetails(self, item):

    if item :
      self.detailsRestraintEntry.set(item.constraint.details)

  def setRestraintDetails(self, event):

    text = self.detailsRestraintEntry.get()
    
    if text:
      text = text.strip()
    else:
      text = None
    
    if self.item:
      self.item.constraint.setDetails(text)
      
  # Button functions
  
  def assignRestraintPeak(self):
  
    if self.item:
      peaks = list(self.item.constraint.peaks)
      if len(peaks) == 1:
        self.guiParent.assignmentPanel()
        self.parent.popups['edit_assignment'].update(peaks[0])

  def assignViolationPeak(self):
  
    if self.violation:
      peaks = list(self.violation.constraint.peaks)
      if len(peaks) == 1:
        self.guiParent.assignmentPanel()
        self.parent.popups['edit_assignment'].update(peaks[0])

  def showViolationPeaks(self):
  
    if self.violation:
      peaks = list(self.violation.constraint.peaks)

      if peaks:
        self.guiParent.viewPeaks(peaks)
  
  def showViolationRestraint(self):
  
    if self.violation:
      restraint = self.violation.constraint
      self.changeRestraintList(restraint.parent)
      items = list(restraint.items)
      self.constraintsMatrix.hilightObject(items[0])
      self.constraintsMatrix.selectObjects(items)
      self.tabbedFrame.select(2)
  
  def makeRestraintPeaks(self):
   
    numFail  = 0
    numPeaks = 0
    constraintLists = self.constraintListMatrix.currentObjects
    for constraintList in constraintLists:
      if constraintList.className == 'DistanceConstraintList':
        peaks = makePeaksFromConstraints(constraintList.sortedConstraints())
        numPeaks += len(peaks)
        
      else:
        numFail += 1
  
    if numFail:
      msg = 'Could not make peaks for (%d) lists: Not distance restraints' % numFail
      showWarning('Warning', msg, parent=self)
      
    showWarning('Complete','Made %d peaks' % numPeaks)
  
  def splitRestraintListAmbig(self):
  
     constraintLists = self.constraintListMatrix.currentObjects
     
     msg = 'Split %d restraint lists by ambiguity?' % len(constraintLists)
     if showOkCancel('Confirm', msg, parent=self):
       for constraintList in constraintLists:
        splitConstraintListAmbigUnambig(constraintList)
  
  def splitRestraintListViol(self):
   
     #if not self.violationListA:
     #  msg = 'No violation list selected'
     #  showWarning('Failure', msg, parent=self)
     #  return 
  
     constraintLists = self.constraintListMatrix.currentObjects
     
     msg = 'Split %d restraint lists by violation?' % len(constraintLists)
     if showOkCancel('Confirm', msg, parent=self):
       for constraintList in constraintLists:
        splitConstraintListViol(constraintList, None)


  def mergeRestraintLists(self):
  
    constraintLists = self.constraintListMatrix.currentObjects
    cName1 = self.constraintListA.className

    if len(constraintLists) > 1:
      msg = 'Really merge %d restraint lists?' % len(constraintLists)
      if showOkCancel('Confirm', msg, parent=self):
        for constraintList in constraintLists:
          if constraintList is self.constraintListA:
            continue
            
          cName2 = constraintList.className
          if cName2 != cName1:
            msg = 'Attempt to merge %s with %s. - %s ignored'
            showWarning('Warning', msg % (cName2, cName1, cName2), parent=self)
            continue
 
          mergeConstraintLists(self.constraintListA, constraintList)
  
  def cleanHaddockProjectConstraints(self, constraintSet):
  
    for haddockProject in self.project.haddockProjects:
      for run in haddockProject.runs:
        if run.nmrConstraintStore is constraintSet:
          for energyTerm in run.haddockEnergyTerms:
            energyTerm.constraintList = None
          run.__dict__['nmrConstraintStore'] = None
          run.checkValid(complete=True)
      
  def deleteRestraintSet(self, *event):

    constraintSets = self.constraintSetMatrix.currentObjects
    
    if len(constraintSets) == 1:
      if showOkCancel('Confirm','Really delete restraint set?', parent=self):
        self.cleanHaddockProjectConstraints(self.constraintSet)
        self.constraintSet.delete()
      else:
        return

    elif len(constraintSets) > 1:
      if showOkCancel('Confirm','Really delete %d restraint sets?' % len(constraintSets), parent=self):
        for constraintSet in constraintSets:
          self.cleanHaddockProjectConstraints(constraintSet)
          constraintSet.delete()
      else:
        return
    
    else:
      return
          
    self.constraintSet  = None
    self.violist    = None
    self.constraintSetButtons.buttons[0].disable()
    self.constraintSetButtons.buttons[1].disable()

  def deleteViolationList(self, *event):

    violationLists = self.violationListMatrix.currentObjects
    
    if len(violationLists) == 1:
      if showOkCancel('Confirm','Really delete violation list?', parent=self):
        self.violationListB.delete()
      else:
        return
   
    elif len(violationLists) > 1:
      if showOkCancel('Confirm','Really delete %d violation lists?' % len(violationLists), parent=self):
        for violationList in violationLists:
           violationList.delete()
      else:
        return
    else:
      return
           
    self.violationListB = None
    
    if self.violationListA and self.violationListA.isDeleted:
      self.violationListA = None
      
    if self.violationListC and self.violationListC.isDeleted:
      self.violationListC = None
    

  def deleteRestraintList(self, *event):

    constrLists = self.constraintListMatrix.currentObjects

    if len(constrLists) == 1:
      if showOkCancel('Confirm','Really delete restraint list?', parent=self):
        self.constraintListA.delete()
      else:
        return
  
    elif len(constrLists) > 1:
      if showOkCancel('Confirm','Really delete %d restraint lists?' % len(constrLists), parent=self):
        for constrList in constrLists:
          constrList.delete()
      else:
        return
    else:
      return
           
    self.constraintList = None
    
  
  def deleteRestraintListSingle(self):

    msg = 'Really delete whole restraint list?'
    if self.constraintListB and showOkCancel('Confirm',msg, parent=self):
      self.constraintListB.delete()
      self.constraintListB = None

  def deleteItem(self):
  
    items = self.constraintsMatrix.currentObjects
    
    if items:
      self.project.__dict__['override'] = True
 
      constraints = []
      for item in items:
        
        constraint = item.constraint
        #print constraint.serial, [makeResonanceGuiName(r) for r in item.resonances]
        item.delete()
        if not constraint.items:
          constraint.delete()
          constraints.append(constraint)
 
      self.constraintSet.checkAllValid()
      self.project.__dict__['override'] = False
      self.constraintSet.__dict__['isModified'] = True
 
 
      for item in items:
        for func in item._notifies.get('delete', []):
          func(item)
 
      for constraint in constraints:
        for func in constraint._notifies.get('delete', []):
          func(constraint)
 
      self.item = None
    
  def deleteConstraint(self, *event):
  
    items = self.constraintsMatrix.currentObjects
    
    if items:
      self.project.__dict__['override'] = True
      
      constraints = {}
      for item in items:
        constraints[item.constraint] = True

      for constraint in constraints.keys():
        constraint.delete()

      self.constraintSet.checkAllValid()
      self.project.__dict__['override'] = False
      self.constraintSet.__dict__['isModified'] = True
      
      for constraint in constraints:
        for func in constraint._notifies.get('delete', []):
          func(constraint)
          
      self.item = None
    
  def viewPeaks(self):
  
    items = self.constraintsMatrix.currentObjects
    
    if items:
      peaks = []
      constraints = {}
      for item in items:
        constraints[item.constraint] = True

      for constraint in constraints.keys():
        peaks.extend(constraint.peaks)
        
      if peaks:
        self.guiParent.viewPeaks(peaks)
  
  def viewViolationList(self):

    self.changeViolationListB(self.violationListA)
    self.tabbedFrame.select(4)

  def viewViolations(self):
    
    self.changeViolationListB(self.violationListB)
    self.tabbedFrame.select(4)
    
  def viewConstraints(self):
  
    self.constraintListB = self.constraintListA
    self.tabbedFrame.select(2)


  def viewForPeakSelection(self):
  
    if self.constraintListB:
      constraints = []
 
      peaks = self.guiParent.currentPeaks
 
      if peaks:
        for constraint in self.constraintListB.constraints:
          for peak in peaks:
            if peak in constraint.peaks:
              if constraint not in constraints:
                constraints.append(constraint)
                break
 
      else:
        showWarning('Warning','No peaks selected in any spectrum windows.')
        return
 
 
      if not constraints:
        showWarning('Warning','No constraints for selected peaks in this list.')

      else:
        items = []
        for constraint in constraints:
          for item in constraint.items:
            items.append(item)

        if items:
          self.constraintsMatrix.hilightObject(items[0], focus=True)
          for item in items[0:]:
            self.constraintsMatrix.hilightObject(item, focus=False)
  
  def mergeDuplicates(self):
  
    if self.constraintListB:
      merged = mergeDuplicateConstraints(self.constraintListB)
      msg = 'Performed %d restraint merges' % len(merged)
      showWarning('Completion', msg, parent=self)

  def updatePeakAssign(self):
    
    if self.constraintListB.className != 'DistanceConstraintList':
      msg = 'Can only update peak assignments for distance restraints'
      showWarning('Stop', msg, parent=self)
     
    
    constraints = set()
    items = self.constraintsMatrix.currentObjects
    
    for item in items:
      constraints.add(item.constraint)
    
    msg = 'Really re-assign peaks based on constraint items?'  
    if showOkCancel('Confirm', msg, parent=self):
      msg2 = 'Keep any existing assignments?'
      keepExisting = showYesNo('Question', msg2, parent=self)
      updatePeaksFromConstraints(constraints, replace=not keepExisting)

  def calcViolations(self):
  
    if self.constraintListB and self.structure:

      self.administerNotifiers(self.unregisterNotify)
      nViol = 0
      if self.violationListA:
        nViol = len(self.violationListA.violations)
    
      self.violationListA = getStructureViolations(self.constraintListB,
                                                   self.structure,
                                                   self.violationListA)
      
      if self.violationListA:
        nViol = len(self.violationListA.violations) - nViol
      else:
        nViol = 0
        
      self.administerNotifiers(self.registerNotify)

      self.updateRestraintsAfter()
      self.updateViolationsAfter()
      
      showInfo('Violation Analysis','Found %d new restraint violations' % (nViol,), parent=self)

  def updateAssignFromPeak(self):

    if self.constraintListB.className == 'DistanceConstraintList':
      constrDict = {}
      
      for item in self.constraintsMatrix.currentObjects:
        constraint = item.constraint
      
        if constraint.peakContribs:
          constrDict[constraint] = None
     
      for constraint in constrDict:
        updateDistConstraintFromPeakAssign(constraint)
     
      self.updateRestraintsAfter()

  def setupResonances(self):
  
    msg = 'Resonances will be created/matched for unlinked restraints'
    if self.constraintListB and showOkCancel('Warning', msg, parent=self):
      nmrConstraintStore = self.constraintListB.parent
      n = len(self.nmrProject.resonances)
      resonances = getConstraintStoreResonances(nmrConstraintStore)
      delta = len(self.nmrProject.resonances)

      showInfo('Done','Made %d new resonances from restraints' % (delta,), parent=self)

  def showStructConnections(self):
  
    items = list(self.constraintsMatrix.currentObjects)
    
    if self.structure and items:
      self.guiParent.viewStructure(self.structure)
      popup = self.guiParent.popups['view_structure']
      popup.clearConnections()
      
      if self.constraintListB.className == 'DihedralConstraintList':
        constrDict = {}
        for item in items:
          constrDict[item.constraint] = None
          
        for constraint in constrDict.keys():
          popup.showResonancesDihedral(constraint.resonances)
      
      else:
        for item in items:
          resonances = list(item.resonances)
          popup.showResonancesConnection(resonances[0], resonances[1])

  # # # Format IO

  """def importAria(self):
    
    fileTypes = [ FileType('All', ['*'])]
    fileSelectPopup = FileSelectPopup(self, file_types = None, show_file=False, 
               title = 'Locate ARIA2 iteration directory', dismiss_text = 'Cancel',
               selected_file_must_exist = False)

    dirName = fileSelectPopup.getDirectory()
    if dirName: 
      if self.constraintSet:
        constraintSet = self.constraintSet
      else:
        constraintSet = makeNmrConstraintStore(self.nmrProject)
      importAria2RunData(dirName, constraintSet)"""
         
    
  def exportAriaList(self):
  
    #TBD: Dihedrals, hbonds
  
    if not self.constraintListA:
      return
 
  
    fileTypes = [  FileType('Table', ['*.tbl']), FileType('All', ['*'])]
    fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
               title = 'Export ARIA file', dismiss_text = 'Cancel',
               selected_file_must_exist = False)

    fileName = fileSelectPopup.getFile() 
    if fileName:
      constraints = self.constraintListA.sortedConstraints()
      exportAriaTbl(constraints, fileName)
 
  def exportCnsList(self):
  
    if not self.constraintListB:
      return
  
    fileTypes = [  FileType('.tbl', ['*.tbl']), FileType('All', ['*'])]
    fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
               title = 'Export CNS file', dismiss_text = 'Cancel',
               selected_file_must_exist = False)

    fileName = fileSelectPopup.getFile() 
    if fileName :
      pass       
      
     
  # General utility
    

  def getStructures(self):

    structures = []
    if self.constraintListB:
      if self.constraintListB.experiments:
        
        molSystems = set()
        for experiment in self.constraintListB.experiments:
          for molSystem in experiment.molSystems:
            molSystems.add(molSystem)
            
      else:
        molSystems = self.project.molSystems
        
      for molSystem in molSystems:
        structures.extend(list(molSystem.structureEnsembles))
  
      structures = [(s.molSystem.code, s.ensembleId, s) for s in structures]
      structures.sort()
  
    return [s[02] for s in structures]
    
  def makeResonanceGuiName(self, resonance):
  
    # For speed - avoids func calls
    if not hasattr(resonance, 'guiName'):
      resonance.guiName = makeResonanceGuiName(resonance, doAtoms=True)
   
    return resonance.guiName
  
  def getConstraintListName(self, constraintList):
    
    if constraintList.name:
      listName = ':%s' % (constraintList.name) 
    else:
      listName = ''
      
    name = '%d:%s%s' % (constraintList.serial,constraintList.className[:-14],listName)
    return name

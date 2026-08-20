
"""
======================COPYRIGHT/LICENSE START==========================

MakeHbondRestraints.py: Part of the CcpNmr Analysis program

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
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.LabelDivider import LabelDivider
from memops.gui.Label import Label
from memops.gui.MessageReporter import showOkCancel, showWarning
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix 
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.AssignmentBasic import makeResonanceGuiName
from ccpnmr.analysis.core.ConstraintBasic import getFixedResonance
from ccpnmr.analysis.core.MoleculeBasic import getResidueCode, getLinkedResidue
from ccpnmr.analysis.core.MoleculeBasic import getBoundAtoms

# TBD
# Nucleic tab - Auto strands Rev Complement
# Fill base pair table automatically, allow some bonds to be switched off

def testHbonds(argServer):

  popup = MakeHbondRestraintsPopup(argServer.parent)
  popup.open()

HBOND_TYPES = (('N-H..O=C','N','O',True),('N-H..N=C','N','N',True),
               ('O-H..O=C','O','O',True),('O-H..N=C','O','N',True),
               ('C-H..O=C','C','O',True),('C-H..N=C','C','N',True),
               ('N-H..O','N','O',False),('N-H..N','N','N',False),
               ('O-H..O','O','O',False),('O-H..N','O','N',False),
               ('C-H..O','C','O',False),('C-H..N','C','N',False),
               ('<Any>',None,None,False))
"""
assign (resid 201 and name N1) (resid 232 and name N3) 2.95 0.2 0.2
assign (resid 201 and name H1) (resid 232 and name N3) 1.95 0.2 0.2
assign (resid 201 and name O6) (resid 232 and name N4) 2.91 0.2 0.2
assign (resid 201 and name O6) (resid 232 and name H42) 1.91 0.2 0.2
assign (resid 201 and name N2) (resid 232 and name O2) 2.86 0.2 0.2
assign (resid 201 and name H22) (resid 232 and name O2) 1.86 0.2 0.2

assign (resid 204 and name N1) (resid 229 and name N3) 2.82 0.2 0.2
assign (resid 204 and name N1) (resid 229 and name H3) 1.82 0.2 0.2
assign (resid 204 and name N6) (resid 229 and name O4) 2.95 0.2 0.2
assign (resid 204 and name H62) (resid 229 and name O4) 1.95 0.2 0.2

assign (resid 11 and name O)(resid 15 and name N ) 2.8 0.4 0.7
assign (resid 11 and name O)(resid 15 and name HN) 1.8 0.3 0.7
assign (resid 12 and name O)(resid 16 and name N ) 2.8 0.4 0.7
assign (resid 12 and name O)(resid 16 and name HN) 1.8 0.3 0.7
assign (resid 13 and name O)(resid 17 and name N ) 2.8 0.4 0.7

opt at 2.79
"""

# Secondary H-bonds

PROTEIN_MOLTYPE = 'protein'

NUCLEIC_MOLTYPES = set(['DNA','RNA','DNA/RNA'])

NUCLEIC_HBONDS = {frozenset(['G','C']):({'G':'O6', 'C':'H42'},
                                        {'G':'H1', 'C':'N3' },
                                        {'G':'H22','C':'O2' }),
                  frozenset(['A','T']):({'A':'N1', 'T':'H3' },
                                        {'A':'H62','T':'O4' }),
                  frozenset(['A','U']):({'A':'N1', 'U':'H3' },
                                        {'A':'H62','U':'O4' }),
                  frozenset(['I','A']):({'I':'06', 'A':'H62'},
                                        {'I':'H1', 'A':'N1' }),
                  frozenset(['I','C']):({'I':'06', 'C':'H42'},
                                        {'I':'H1', 'C':'N3' }),
                  frozenset(['I','U']):({'I':'06', 'U':'H3' },
                                        {'I':'H1', 'U':'O4' }),
                  frozenset(['I','T']):({'I':'06', 'T':'H3' },
                                        {'I':'H1', 'T':'O4' }),
                 }

ALPHA = u'\u03B1'
PI = u'\u03C0'
#THREE_TEN = u'3\u2080\u2081'
THREE_TEN = '3_10'

HELIX_TYPES = [ALPHA, PI, THREE_TEN]
SHEET_TYPES = ['anti-parallel','parallel']

HELIX_DELTAS = {3:THREE_TEN, 4:ALPHA, 5:PI}

DEFAULT_UPPER = 2.70
DEFAULT_TARGET = 2.20
DEFAULT_LOWER = 1.73  
DEFAULT_BOND = 1.0  

PROLYL = set(['Pro','Hyp'])

class MakeHbondRestraintsPopup(BasePopup):
  """
  **Create Hydrogen Bond Restraints**
  
  This popup allows the user to make simple distance-based hydrogen bond
  restraints for structure determination. These restraints are essentially the
  same as other distance restraints, although a single hydrogen bond is usually
  restrained with a complementary pair of restraints: one restraint is between
  the hydrogen and the acceptor atom site (usually a backbone carbonyl oxygen in
  proteins) and the other is between donor and acceptor heavy atoms (e.g.
  between protein backbone amide nitrogen and carbonyl oxygen). This combination
  of restraints work together to give an energy minimum for the restraint when
  the donor, hydrogen and acceptor atoms are co-linear; a known optimum for
  hydrogen bonds.

  To make hydrogen bond restraints the user first selects various options from
  the upper pulldown menus: a restraint set that groups the h-bond restraints
  with other restraint lists, an H-bond type restraint list that the individual
  restraints are place in, and a molecular system to specify which set of atoms
  the restraints are between. 

  Next the "New H Bond Type" should be checked, but the default "N-H..O=C" is
  already set for polypeptide backbone H-bonds. The various default distances can
  be adjusted to specify what the restrained values will initially be, although
  the user can make manual adjustments (including by setting and propagating new
  defaults) afterwards. The "Lower Limit", "Target" and "Upper Limit" distances
  refer to the separation between the hydrogen and its acceptor. The distance
  between the donor heavy atom and the acceptor atom are made automatically
  based on the "X-H" bond length setting. Specifically, the target distance for
  the heavy atoms' restraint will be the hydrogen to acceptor distance plus
  the bond length. Also, the upper and lower limits for these atoms will follow
  the same proportions as the hydrogen to acceptor values; the values will be
  the same fractions of the target distances.

  Hydrogen bond restraints can be added using [Add New] and then setting the
  Atom A and Atom B values by double clicking in the table and selecting the
  required atom sites. Note that generally the user only needs to set the
  hydrogen to acceptor restraints. The heavy atoms' restraint is made afterwards
  by pressing the [Make Co-linear Restraints] button. It notable that this
  function will actually add any missing restraint of a co-linear pair.

  **Evidence**
  
  This system relies upon the user having a reasonable and unbiased judgement
  about the location of hydrogen bond in a molecular structure. Whenever H-bond
  restraints are made the user should always be aware of the evidence for
  particular hydrogen bond. It is possible to get nice looking but otherwise
  fictitious structures by using hydrogen bond restraints too liberally.
  
  Often H-bond evidence often comes from knowledge of protein secondary
  structure (particularly alpha-helix for i to i+4 links ), which in turn can be
  deduced using chemical shift data, for example using DANGLE_, or confident
  comparative modelling of protein homologues. Also, there are some NMR
  experiments that will exploit spin couplings over hydrogen bonds themselves.
  Whatever the source of confidence about hydrogen bonds, using them as
  restraints generally requires that consistency with any other structural data
  is checked, for example by violation analysis. With NOE data for example it is
  prudent to generate initial protein structures without any H-bond information
  and only add-in H-bond restraints if they would be consistent with the initial
  NOE-derived conformations. In this way 'low evidence' hydrogen bonds can be
  added to improve the later stages of structure refinement without having an
  undue influence of the general fold of a structure.

  .. _DANGLE: DangleGui.html

  """

  def __init__(self, parent, *args, **kw):

    self.restraintSet = None
    self.restraintList = None
    self.restraint = None
    self.guiParent = parent
    self.waiting = False
    self.waitingNc = False
    self.molSystem = None
    self.tempResonanceH = None
    self.tempResonanceX = None
    self.hBondType = HBOND_TYPES[0]
    self.helix = None
    self.helices = []
    self.strand = None
    self.strands = []
  
    BasePopup.__init__(self, parent=parent, 
                       title='Structure : Make H Bond Restraints', **kw)

  def body(self, guiFrame):

    self.geometry('800x500')

    guiFrame.expandGrid(0,0)
    
    tipTexts = ['Allows the setup of individual hydrogen bond restraints',]
    options = ['H Bonds Table'] #,'Peptide','DNA/RNA']
    self.tabbedFrame = TabbedFrame(guiFrame, options=options, grid=(0,0), tipTexts=tipTexts)
    tFrame = self.tabbedFrame.frames[0]
    
    # Do not display until ready
    pFrame = Frame(guiFrame)
    nFrame = Frame(guiFrame)
    
    # Main H bonds table
    
    tFrame.expandGrid(2,0)
        
    frame = Frame(tFrame, grid=(0,0))
    frame.expandGrid(None, 8)
    
    div = LabelDivider(frame, text='Default Distances',grid=(0,0), gridSpan=(1,6))
    label = Label(frame, text='Lower Limit:', grid=(1,0))
    tipText = 'Set the lower distance limit, in Angstrom, to apply to new hydrogen bond restraints; each limit may be subsequently adjusted'
    self.lowerEntry = FloatEntry(frame, text=DEFAULT_LOWER, grid=(1,1),
                                 width=8, tipText=tipText)
    label = Label(frame, text='Target:', grid=(1,2))
    tipText = 'Set best-fit/target distance, in Angstrom, to apply to new hydrogen bond restraints; each limit may be subsequently adjusted'
    self.targetEntry = FloatEntry(frame, text=DEFAULT_TARGET, grid=(1,3),
                                  width=8, tipText=tipText)
    label = Label(frame, text='Upper Limit:', grid=(1,4))
    tipText = 'Set the upper distance limit, in Angstrom, to apply to new hydrogen bond restraints; each limit may be subsequently adjusted'
    self.upperEntry = FloatEntry(frame, text=DEFAULT_UPPER, grid=(1,5),
                                 width=8, tipText=tipText)

    label = Label(frame, text='X-H Bond Length:', grid=(1,6))
    tipText = 'The distance between the hydrogen and its covalently bound partner atom; user to make co-linearity restraints'
    self.bondEntry = FloatEntry(frame, text=DEFAULT_BOND, grid=(1,7),
                                width=8, tipText=tipText)

    frame = Frame(tFrame, grid=(0,1))
    frame.expandGrid(0,0)
    label = Label(frame, text='New H Bond Type:', grid=(1,1))
    
    texts = [x[0] for x in HBOND_TYPES]
    tipText = 'Selects which kind of hydrogen bond is being made, with respect to the kinds of atoms involved'
    self.hbondTypePulldown = PulldownList(frame, grid=(1,2), texts=texts,
                                          objects=HBOND_TYPES, tipText=tipText,
                                          callback=self.changeHbondType)
    
    div = LabelDivider(tFrame, text='Hydrogen Bond Restraints',
                       gridSpan=(1,2), grid=(1,0))

    self.distLowerEntry = FloatEntry(self, text=0.0, width=8,
                                returnCallback=self.setLower)
    
    self.distTargetEntry = FloatEntry(self, text=0.0, width=8,
                                returnCallback=self.setTarget)
    
    self.distUpperEntry = FloatEntry(self, text=0.0, width=8,
                                returnCallback=self.setUpper)
   
    self.donorPulldown = PulldownList(self, callback=self.setAtomA)
    
    self.acceptorPulldown = PulldownList(self, callback=self.setAtomB)
    
    tipTexts = ['The donor hydrogen atom involved in the hydrogen bond restraint',
                'The acceptor atom involved in the hydrogen bond restraint',
                'The lower distance bound value for the restraint, in Angstrom',
                'The best-fit/target distance for the hydrogen bond restraint, in Angstrom',
                'The upper distance bound value for the restraint, in Angstrom']
    headingList = ['Atom A', 'Atom B', 'Lower Limit',
                   'Target Distance', 'Upper Limit',]
    editWidgets = [self.donorPulldown, self.acceptorPulldown,
                   self.distLowerEntry, self.distTargetEntry,
                   self.distUpperEntry,]
    editGetCallbacks = [self.getAtomA, self.getAtomB,
                        self.getLower, self.getTarget, self.getUpper,]
    editSetCallbacks = [self.setAtomA, self.setAtomB,
                        self.setLower, self.setTarget, self.setUpper,]
    self.hbondTable = ScrolledMatrix(tFrame, headingList=headingList,
                                     callback=self.selectRestraint,
                                     editWidgets=editWidgets,
                                     editSetCallbacks=editSetCallbacks,
                                     editGetCallbacks=editGetCallbacks,
                                     deleteFunc=self.deleteRestraints,
                                     grid=(2,0), gridSpan=(1,2),
                                     multiSelect=True,
                                     tipTexts=tipTexts)
    
    tipTexts = ['Add a new hydrogen bond restraint to the selected list; the acceptor and donor atoms etc. are set subsequently',
                'Delete selected hydrogen bond restraint from the list',
                'Set the selected hydrogen bond restraints to have the default upper, lower and target distances',
                'Delete the whole hydrogen bond restraint list',
                'Make any complementary restraints which restrain co-linearity of X-H..Y bonds']
    texts = ['Add New','Delete Selected',
             'Set Default\nDistances',
             'Delete List','Make Co-linear\nRestraints']
    commands = [self.newRestraint, self.deleteRestraints,
                self.setDefaultDistances, self.deleteRestraintList,
                self.makeCoLin]
    self.restraintButtons = ButtonList(tFrame, texts=texts, gridSpan=(1,2),
                                       commands=commands, grid=(3,0),
                                       tipTexts=tipTexts)
    
    # # # # # # #  P E P T I D E S  # # # # # # #
    
    # Helices table

    self.helixStartPulldown = PulldownList(self, callback=self.setHelixStart)
    self.helixEndPulldown = PulldownList(self, callback=self.setHelixEnd)
    self.helixTypePulldown = PulldownList(self, callback=self.setHelixType,
                                          objects=[4,5,3], texts=HELIX_TYPES )
    
    pFrame.expandGrid(1,0)
    div = LabelDivider(pFrame, text='Helices', grid=(0,0))

    headingList = ['#', 'Start', 'End', 'Type',]
    editWidgets = [None,  self.helixStartPulldown,
                   self.helixEndPulldown,
                   self.helixTypePulldown]
    editGetCallbacks = [None, self.getHelixStart, 
                        self.getHelixEnd, self.getHelixType]
    editSetCallbacks = [None, self.setHelixStart, 
                        self.setHelixEnd, self.setHelixType]
    self.helixTable = ScrolledMatrix(pFrame, headingList=headingList,
                                     callback=self.selectHelix,
                                     editWidgets=editWidgets,
                                     editSetCallbacks=editSetCallbacks,
                                     editGetCallbacks=editGetCallbacks,
                                     grid=(1,0), 
                                     multiSelect=True)
    
    texts = ['Add Helix','Delete Helix']
    commands = [self.addHelix, self.removeHelix]
    self.helixButtons = ButtonList(pFrame, texts=texts,
                                   commands=commands) #, grid=(2,0))
    # Sheets table
    
    self.strandStart1Pulldown = PulldownList(self, callback=self.setStrand1Start)
    self.strandEnd1Pulldown = PulldownList(self, callback=self.setStrand1End)
    self.strandStart2Pulldown = PulldownList(self, callback=self.setStrand2Start)
    self.strandEnd2Pulldown = PulldownList(self, callback=self.setStrand2End)
    self.strandTypePulldown = PulldownList(self, callback=self.setStrandType,
                                           objects=[-1, 1], texts=SHEET_TYPES )
    pFrame.expandGrid(4,0)
    div = LabelDivider(pFrame, text='Paired Strands', grid=(3,0))

    headingList = ['#', 'Start A', 'End A', 'Start B', 'End B', 'Type',]
    editWidgets = [None,
                   self.strandStart1Pulldown, self.strandEnd1Pulldown,
                   self.strandStart2Pulldown, self.strandEnd2Pulldown,
                   self.strandTypePulldown]
    editGetCallbacks = [None, self.setStrand1Start, self.setStrand1End,
                        self.setStrand2Start, self.setStrand2End,
                        self.setStrandType]
    editSetCallbacks = [None, self.getStrand1Start, self.getStrand1End, 
                        self.getStrand2Start, self.getStrand2End,
                        self.getStrandType]
    self.sheetTable = ScrolledMatrix(pFrame, headingList=headingList,
                                     callback=self.selectStrands,
                                     editWidgets=editWidgets,
                                     editSetCallbacks=editSetCallbacks,
                                     editGetCallbacks=editGetCallbacks,
                                     grid=(4,0), 
                                     multiSelect=True)
    texts = ['Add Strands','Delete Strands']
    commands = [self.addStrands, self.removeStrands]
    self.sheetButtons = ButtonList(pFrame, texts=texts,
                                   commands=commands) #, grid=(5,0))
   
    div = LabelDivider(pFrame, text='Generate Restraints') #, grid=(6,0))
   
    texts = ['Helices', 'Strands', 'Helices and Strands',]
    commands = [self.makeRestraintsFromHelices,
                self.makeRestraintsFromSheets,
                self.makeRestraintsFromNonCovalent]
                
    self.pepRestraintButtons = ButtonList(pFrame, texts=texts, 
                                          commands=commands) #, grid=(7,0))
                                       
    # # # # # # #  N U C L E I C   A C I D S  # # # # # # #
    # Duplex regions user specified
    # Perfect RC filled in
    #
    # Base pairs default from alignment
    #
    #
    
    nFrame.expandGrid(1,0)
    div = LabelDivider(nFrame, text='Duplex Regions', grid=(0,0))

    headingList = ['#', 'Start A', 'End A', 'Start B', 'End B']
    editWidgets = [None, None, None, None]
    editGetCallbacks = [None, None, None, None]
    editSetCallbacks = [None, None, None, None]
    self.duplexTable = ScrolledMatrix(nFrame, headingList=headingList,
                                     callback=None,
                                     editWidgets=editWidgets,
                                     editSetCallbacks=editSetCallbacks,
                                     editGetCallbacks=editGetCallbacks,
                                     grid=(1,0), 
                                     multiSelect=True)

    div = LabelDivider(nFrame, text='Base Pairs', grid=(2,0))

    nFrame.expandGrid(3,0)
    headingList = ['Residue A','Residue B','H Bonds']
    editWidgets = [None, None, None]
    editGetCallbacks = [None, None, None]
    editSetCallbacks = [None, None, None]
    self.basePairTable = ScrolledMatrix(nFrame, headingList=headingList,
                                        callback=None,
                                        editWidgets=editWidgets,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        grid=(3,0),
                                        multiSelect=True)


    # main
    
    frame = self.tabbedFrame.sideFrame
    frame.grid_columnconfigure(0, weight=0)
    frame.grid_columnconfigure(6, weight=1)
    
    label = Label(frame, text='Restraint Set:', grid=(0,0))
    tipText = 'Selects which restraint set the selected, or new, hydrogen bond restraint list resides in'
    self.restraintSetPulldown = PulldownList(frame, grid=(0,1), tipText=tipText,
                                             callback=self.selectRestraintSet)

    label = Label(frame, text='H Bond List:', grid=(0,2))
    tipText = 'Selects which hydrogen bond restraint list to display; may be set to "<New>" to put restraints in a new list'
    self.restraintListPulldown = PulldownList(frame, grid=(0,3), tipText=tipText,
                                              callback=self.selectRestraintList)
                                              
    label = Label(frame, text='MolSystem:', grid=(0,4))
    tipText = 'Selects which molecular system (group of chains) to make restraints for; sets which atoms are available for selection'
    self.molSystemPulldown = PulldownList(frame, grid=(0,5), tipText=tipText,
                                      callback=self.selectMolSystem)
    
    buttons = UtilityButtonList(frame, expands=True)
    buttons.grid(row=0, column=7, sticky='e')
    
    self.updateAfter()
    
    self.administerNotifiers(self.registerNotify)

  def makeCoLin(self):

    bondLen = self.bondEntry.get() or 1.0
    restraintList = self.restraintList
    getItemResonances = self.getItemResonances
    findResonance = self.findResonance
    triples = {}
        
    if restraintList:
      for restraint in restraintList.constraints:
        item = restraint.findFirstItem()
        if not item:
          continue
        
        r1, r2 = getItemResonances(item)
        rs1 = r1.resonanceSet
        rs2 = r2.resonanceSet
        
        if not (rs1 and rs2):
          continue
        
        
        atom1 = rs1.findFirstAtomSet().findFirstAtom()
        atom2 = rs2.findFirstAtomSet().findFirstAtom()
        
        if r1.isotopeCode == '1H':
          atomD = getBoundAtoms(atom1)[0]
          atomH = atom1
          atomA = atom2
          index = 0
        
        elif r2.isotopeCode == '1H':
          atomD = getBoundAtoms(atom2)[0]
          atomH = atom2 
          atomA = atom1
          index = 0
         
        else:
          boundH1 = [a for a in getBoundAtoms(atom1) \
                     if a.chemAtom.elementSymbol == 'H']
          boundH2 = [a for a in getBoundAtoms(atom2) \
                     if a.chemAtom.elementSymbol == 'H']
          
          if len(boundH1 + boundH2) == 1:
            index = 1
            if boundH1:
              atomH = boundH1[0]
              atomD = atom1
              atomA = atom2
            else:
              atomH = boundH2[0]
              atomD = atom2
              atomA = atom1
          
          else:
            continue

        key = (atomD, atomH, atomA)
        
        if key not in triples:
          triples[key] = [None, None]
        
        triples[key][index] = restraint

  
    for key in triples:
      atomD, atomH, atomA = key
      restraintH, restraintX = triples[key]
    
      if restraintH is None:
        targetX = restraintX.targetValue \
                  or (0.5*(restraintX.upperLimit + restraintX.lowerLimit)) 
        
        target = targetX - bondLen
        upper = target * restraintX.upperLimit/targetX
        lower = target * restraintX.lowerLimit/targetX 
        restraintH = restraintList.newHBondConstraint(targetValue=target,
                                                      upperLimit=upper,
                                                      lowerLimit=lower)
        r1 = findResonance(atomH)
        r2 = findResonance(atomA)
        restraintH.newHBondConstraintItem(resonances=(r1,r2))

      elif restraintX is None:
        targetH = restraintH.targetValue \
                  or (0.5*(restraintH.upperLimit + restraintH.lowerLimit)) 
        
        target = targetH + bondLen
        upper = target * restraintH.upperLimit/targetH
        lower = target * restraintH.lowerLimit/targetH 
        restraintH = restraintList.newHBondConstraint(targetValue=target,
                                                      upperLimit=upper,
                                                      lowerLimit=lower)
        r1 = findResonance(atomD)
        r2 = findResonance(atomA)
        restraintH.newHBondConstraintItem(resonances=(r1,r2))
    


  def updateMolSystems(self, obj=None):
  
    index = 0
    molSystems = self.project.sortedMolSystems()
    molSystem = self.molSystem
    names = [ms.code for ms in molSystems]
    
    if molSystems:
      if molSystem not in molSystems:
        molSystem = molSystems[0]
      index = molSystems.index(molSystem)

    if self.molSystem is not molSystem:
      self.selectMolSystem(molSystem) 

    self.molSystemPulldown.setup(names, molSystems, index)

  def selectMolSystem(self, molSystem):
  
    if molSystem is not self.molSystem:
      self.helices = getHelicalHbondSegments(molSystem)
      self.strands = getPairedHbondSheets(molSystem)
      self.molSystem = molSystem
      self.updateAfter()

  def updateRestraintSets(self, obj=None):
  
    index = 0
    names = []
    restraintSets = self.nmrProject.sortedNmrConstraintStores()
    restraintSet = self.restraintSet
    
    if restraintSets:
      if restraintSet not in restraintSets:
        restraintSet = restraintSets[-1]
 
      index = restraintSets.index(restraintSet)
      names = ['%d' % rs.serial for rs in restraintSets]
    
    if restraintSet is not self.restraintSet:
      self.selectRestraintSet(restraintSet)
    
    self.restraintSetPulldown.setup(names, restraintSets, index)

  def updateRestraintLists(self, obj=None):
  
    index = 0
    restraintLists = []
    restraintList = self.restraintList
    
    if self.restraintSet:
      getLists = self.restraintSet.findAllConstraintLists
      hbondLists = getLists(className='HBondConstraintList')
      sortList = [('%d' % l.serial, l) for l in hbondLists]
      sortList.sort()
      
      restraintLists = [None,] + [x[1] for x in sortList]
      names = ['<New>',] + [x[0] for x in sortList]
    
      if restraintList not in restraintLists:
        restraintList = restraintLists[-1]
    
      index = restraintLists.index(restraintList)
    
    else:
      restraintList = None
       
    if restraintList is not self.restraintList:
      self.selectRestraintList(restraintList)
    
    self.restraintListPulldown.setup(names, restraintLists, index)

  def selectRestraintList(self, restraintList):
  
    if restraintList is not self.restraintList:
      self.restraintList = restraintList
      self.updateAfter()

  def selectRestraintSet(self, restraintSet):
  
    if restraintSet is not self.restraintSet:
      self.restraintSet = restraintSet
      getFirst = restraintSet.findFirstConstraintList
      self.restraintList = getFirst(className='HBondConstraintList')
      self.updateRestraintLists()
      self.updateAfter()

  # # # # # # # # Main table # # # # # # # # 
  
  def findResonance(self, atom):
  
    atomSet = atom.atomSet
    
    if not atomSet:
      atomSet = self.nmrProject.newAtomSet(name=atom.name, atoms=(atom,))
  
    resonance = None
    for resonanceSet in atomSet.resonanceSets:
      resonances = list(resonanceSet.resonances)
  
      if len(resonances) != 1:
        continue
        
      atomSets = resonanceSet.atomSets
      if len(atomSets) != 1:
        continue
      
      resonance = resonances.pop()  
  
    if not resonance:
      nmrProj = self.nmrProject
      chemElement = atom.chemAtom.chemElement

      isotopes = [(iso.abundance, iso.massNumber) for iso in chemElement.isotopes]
      isotopes.sort()
      
      isotope = '%d%s' % (isotopes[-1][1], chemElement.symbol)

      resonance  = nmrProj.newResonance(isotopeCode=isotope)
      resonanceSet = nmrProj.newResonanceSet(atomSets=(atomSet,),
                                             resonances=(resonance,))
    
    fixedResonance = getFixedResonance(self.restraintSet, resonance)
   
    return fixedResonance
  
  def setupAtomOptions(self, pulldownList, resonance, isDonor=True):
  
    index = 0 # No change unless selected
    names = []
    objects = []
    cats = []
    
    namesAppend = names.append
    objectsAppend = objects.append
    catsAppend = cats.append
    
    if resonance:
      namesAppend(self.makeResonanceGuiName(resonance))
    else:
      namesAppend('<None>')  
    
    objectsAppend(None) # First item (default) => no change
    cats.append(None)
    
    if self.molSystem:
      text, atomTypeA, atomTypeB, isDoubleBound = self.hBondType
    
      if isDonor:
        atomType = atomTypeA
      else:
        atomType = atomTypeB   
    
      atomList = []
      atomListAppend = atomList.append
      for chain in self.molSystem.sortedChains():
        chCode = chain.code
      
        for residue in chain.sortedResidues():
          seqCode = residue.seqCode
          ccpCode = getResidueCode(residue)
 
          for atom in residue.sortedAtoms():
            chemAtom = atom.chemAtom
            element = chemAtom.elementSymbol
 
            if isDonor:
              if element == 'H':
                if atomType:
                  bond = chemAtom.findFirstChemBond()
                  chemAtoms = list(bond.chemAtoms)
                  chemAtoms.remove(chemAtom)
 
                  if chemAtoms[0].elementSymbol == atomType:
                    atomListAppend((chCode, seqCode, ccpCode, atom))
 
                else:
                    atomListAppend((chCode, seqCode, ccpCode, atom))
 
            elif not atomType:
              if element != 'H':
                atomListAppend((chCode, seqCode, ccpCode, atom))
 
            elif element == atomType:
 
              if isDoubleBound:
                for bond in chemAtom.chemBonds:
                  if bond.bondType in ('double', 'aromatic'):
                    chemAtoms = list(bond.chemAtoms)
                    chemAtoms.remove(chemAtom)
 
                    if chemAtoms[0].elementSymbol == 'C':
                      atomListAppend((chCode, seqCode, ccpCode, atom))
                      break
 
              else:
                atomListAppend((chCode, seqCode, ccpCode, atom))
      
      #atomList.sort()
      for chCode, seqCode, ccpCode, atom in atomList:
        namesAppend('%s %d%s%s' % (chCode, seqCode, ccpCode, atom.name))
        objectsAppend(atom)
        x = 10*(seqCode//10)
        cats.append('%s %d-%d' % (chCode, x,x+9))
  
    pulldownList.setup(names, objects, index, categories=cats)
  
  def changeResonances(self, restraint, resonanceA, resonanceB): 
  
    resonances = set([resonanceA, resonanceB])
    
    if restraint.findFirstItem(resonances=resonances):
      return
    
    for item in restraint.items:
      item.delete()
    
    restraint.newHBondConstraintItem(resonances=resonances)
  
  def getAtomA(self, restraint):
  
    item = restraint.findFirstItem()
    
    if item:
      r1, r2 = self.getItemResonances(item)
    else:
      r1 = None
      
    self.setupAtomOptions(self.donorPulldown, r1, True) 
  
  def setAtomA(self, *opt):
    
    atom = self.donorPulldown.getObject()
    if self.restraint and atom:
      item = self.restraint.findFirstItem()
      
      if item:
        r1, r2 = self.getItemResonances(item)
      else:
        r2 = self.getTempResonanceX(self.restraintSet)
    
      r1 = self.findResonance(atom)
   
      self.changeResonances(self.restraint, r1, r2)
    
  def getAtomB(self, restraint):
  
    item = restraint.findFirstItem()
    
    if item:
      r1, r2 = self.getItemResonances(item)
    else:
      r2 = None
      
    self.setupAtomOptions(self.acceptorPulldown, r2, False) 
  
  def setAtomB(self, *opt):
    
    atom = self.acceptorPulldown.getObject()
    if self.restraint and atom:
      item = self.restraint.findFirstItem()
      
      if item:
        r1, r2 = self.getItemResonances(item)
      else:
        r1 = self.getTempResonanceH(self.restraintSet)
    
      r2 = self.findResonance(atom)
   
      self.changeResonances(self.restraint, r1, r2)
  
  def changeHbondType(self, selection):
  
    self.hBondType = selection    
  
  def getLower(self, restraint):

    if self.restraint:
      value = restraint.lowerLimit
      self.distLowerEntry.set(value)
  
  def setLower(self, *whatever):
  
    if self.restraint:
      value = self.distLowerEntry.get() or 0.0
      
      if (value < 0.0) or (value > 4.0):
        return
      
      if value > self.restraint.targetValue:
        self.restraint.targetValue = value

      if value > self.restraint.upperLimit:
        self.restraint.upperLimit = value
      
      self.restraint.lowerLimit = value
  
  def getTarget(self, restraint):

    if self.restraint:
      value = restraint.targetValue
      self.distTargetEntry.set(value)
  
  def setTarget(self, *whatever):
  
    if self.restraint:
      value = self.distTargetEntry.get() or 0.0
      
      if (value < 0.0) or (value > 4.0):
        return

      if value > self.restraint.upperLimit:
        self.restraint.upperLimit = value

      if value < self.restraint.lowerLimit:
        self.restraint.lowerLimit = value

      self.restraint.targetValue = value
  
  def getUpper(self, restraint):

    if self.restraint:
      value = restraint.upperLimit
      self.distUpperEntry.set(value)
  
  def setUpper(self, *whatever):
  
    if self.restraint:
      value = self.distUpperEntry.get() or 0.0
      
      if (value < 0.0) or (value > 4.0):
        return
      
      if value < self.restraint.targetValue:
        self.restraint.targetValue = value

      if value < self.restraint.lowerLimit:
        self.restraint.lowerLimit = value
        
      self.restraint.upperLimit = value

  def selectRestraint(self, obj, row, col):
  
    self.restraint = obj
    

  def newRestraint(self):
  
    if self.restraintSet:
      if not self.restraintList:
        self.restraintList = self.restraintSet.newHBondConstraintList()
        self.updateRestraintLists()
        
      lower, target, upper = self.getDefaultDistances()
      restraint = self.restraintList.newHBondConstraint(targetValue=target,
                                                        upperLimit=upper,
                                                        lowerLimit=lower)
                                                        
      self.update_idletasks()
      self.hbondTable.selectObject(restraint)
    else:
      showWarning('Failure', 'No restraint set selected', parent=self)

  def deleteRestraints(self, event=None):
  
    if self.restraint:
      restraints = self.hbondTable.currentObjects
      msg = 'Really remove %d hydrogen bond restraints?' % len(restraints)
    
      if showOkCancel('Confirm', msg, parent=self):
        for restraint in restraints:
          restraint.delete()
    
    else:
      showWarning('Failure', 'No selected restraints', parent=self)

  def deleteRestraintList(self, event=None):
  
    restraintList = self.restraintList

    if restraintList:
      msg = 'Really delete the current hydrogen bond restraint list?'
    
      if showOkCancel('Confirm', msg, parent=self):
        restraintList.delete()
        
 
  def getDefaultDistances(self):
  
    lower = self.lowerEntry.get() or DEFAULT_LOWER
    target = self.targetEntry.get() or DEFAULT_TARGET
    upper = self.upperEntry.get() or DEFAULT_UPPER
    
    return lower, target, upper

  def setDefaultDistances(self):

    if self.restraint:
      restraints = self.hbondTable.currentObjects
      msg = 'Really set %d hydrogen bond distances?' % len(restraints)
      
      lower, target, upper = self.getDefaultDistances()
    
      if showOkCancel('Confirm', msg, parent=self):
        for restraint in restraints:
          restraint.lowerLimit = lower
          restraint.targetValue = target
          restraint.upperLimit = upper

    else:
      showWarning('Failure', 'No selected restraints', parent=self)
  
  def getItemResonances(self, item):
  
    r1, r2 = item.resonances
    isotope1 = r1.isotopeCode
    isotope2 = r2.isotopeCode
    
    
    if isotope1 == '1H':
      rD, rA = r1, r2

    elif isotope2 == '1H':
      rA, rD = r1, r2

    else:
      sortList = [(isotope1, r1), (isotope2, r2)]
      sortList.sort()
      
      rD, rA = [x[1] for x in sortList] 
    
    return rD, rA
      
  def updateMain(self):
  
    textMatrix = []
    objectList = []
    textMatrixAppend = textMatrix.append
    objectListAppend = objectList.append
    getResonanceName = self.makeResonanceGuiName
    getItemResonances = self.getItemResonances
    
    if self.restraintList:
      
      for restraint in self.restraintList.sortedConstraints():
        items = restraint.items
        
        dRes = []
        aRes = []
      
        for item in items:
          rD, rA = getItemResonances(item)
     
          dRes.append(getResonanceName(rD))
          aRes.append(getResonanceName(rA))
      
        donor = '/'.join(dRes)
        acceptor = '/'.join(aRes)
      
        datum = [donor, acceptor,
                 restraint.lowerLimit,
                 restraint.targetValue,
                 restraint.upperLimit]
      
        textMatrixAppend(datum)
        objectListAppend(restraint)
        
        
        
    
    self.hbondTable.update(textMatrix=textMatrix,
                           objectList=objectList)
    
  def makeResonanceGuiName(self, resonance):
  
    # For speed - avoids func calls
    if not hasattr(resonance, 'guiName'):
      resonance.guiName = makeResonanceGuiName(resonance, fullName=True)
   
    return resonance.guiName
  
  def getTempResonanceH(self, restraintSet):
  
    if not self.tempResonanceH:
      fResonance = restraintSet.findFirstFixedResonance(name='HBondTemp',
                                                        isotopeCode='1H')
    
      if not fResonance:
        fResonance = restraintSet.newFixedResonance(name='HBondTemp',
                                                    isotopeCode='1H')
      fResonance.guiName = 'H?' 
      self.tempResonanceH = fResonance
    
    return self.tempResonanceH
  
  def getTempResonanceX(self, restraintSet):
  
    if not self.tempResonanceX:
      fResonance = restraintSet.findFirstFixedResonance(name='HBondTemp',
                                                        isotopeCode='unknown')
    
      if not fResonance:
        fResonance = restraintSet.newFixedResonance(name='HBondTemp',
                                                    isotopeCode='unknown')
      fResonance.guiName = '?' 
      self.tempResonanceX = fResonance
    
    return self.tempResonanceX

  # # # # # # # General # # # # # # # 

  def updateAfter(self):
  
    if self.waiting:
      return
    
    self.waiting = True  
    self.after_idle(self.updateAll)  

  def updateRestraintsAfter(self, restraint):
     
    if not self.waiting:
      if restraint.parent is self.restraintList:
        self.updateAfter()
        
  def updateItemsAfter(self, item):

    self.updateRestraintsAfter(item.constraint)

  def updateResonanceSetAfter(self, fixedResonanceSet):
  
    if self.waiting:
      return
  
    restraintList = self.restraintList
  
    if fixedResonanceSet.nmrConstraintStore is self.restraintSet:
      for resonance in fixedResonanceSet.resonances:
        for item in resonance.pairwiseConstraintItems:
          if item.constraint.constraintList is restraintList:
            self.updateAfter()
            return      

  def updateAll(self, obj=None):
    
    self.updateRestraintSets()
    
    self.updateMolSystems()
    
    self.updateMain()
        
    self.updateNonCovalent()
  
    self.waiting = False

  def administerNotifiers(self, notifyFunc):

    for func in ('__init__', 'delete',):
      notifyFunc(self.updateMolSystems, 'ccp.molecule.MolSystem.Chain', func)
      notifyFunc(self.updateRestraintSets, 'ccp.nmr.NmrConstraint.NmrConstraintStore', func)
      notifyFunc(self.updateRestraintLists, 'ccp.nmr.NmrConstraint.HBondConstraintList', func)
      notifyFunc(self.updateNonCovalentAfter, 'ccp.molecule.MolSystem.NonCovalentBond', func)
    
    notifyFunc(self.updateResonanceSetAfter, 'ccp.nmr.NmrConstraint.FixedResonanceSet', '')
              
    for func in ('__init__', 'delete','setLowerLimit',
                 'setTargetValue','setUpperLimit'):
      notifyFunc(self.updateRestraintsAfter, 'ccp.nmr.NmrConstraint.HBondConstraint', func)
        
    for func in ('__init__', 'delete','setResonances'):
      notifyFunc(self.updateItemsAfter, 'ccp.nmr.NmrConstraint.HBondConstraintItem', func)

    notifyFunc(self.updateAll,'ccp.molecule.MolSystem.Residue', 'setSeqCode')
    
  def open(self):
   
   BasePopup.open(self)
       
  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
    self.helices = []
    self.strands = []

    BasePopup.destroy(self)
    
  # # # # # # # # Peptides # # # # # # # #
  
  def makeRestraintsFromHelices(self):
  
    self.makeRestraintsFromNonCovalent(doHelices=True, doSheets=False)

  def makeRestraintsFromSheets(self):
  
    self.makeRestraintsFromNonCovalent(doHelices=False, doSheets=True)
  
  def makeRestraintsFromNonCovalent(self, doHelices=True, doSheets=True):
  
    if self.restraintSet and self.molSystem:
      if not self.restraintList:
        self.restraintList = self.restraintSet.newHBondConstraintList()
        self.updateRestraintLists()
      
      existing = set()
      for restraint in self.restraintList.constraints:
        item = restraint.findFirstItem()
        
        if item:
          existing.add(item.resonances) # Add frozenset
      
      findResonance = self.findResonance  
      lower, target, upper = self.getDefaultDistances()
      makeRestraint = self.restraintList.newHBondConstraint
      atomPairs = getProteinBackboneHbonds(self.molSystem)
      
      elements = []
      if doHelices:
        elements += self.helices
      if doSheets:
        elements += self.strands  
       
      for element in elements:
        for ncBond in element.nonCovalentBonds:
          atomA, atomB = ncBond.atoms
      
          resonances = frozenset([findResonance(atomA), findResonance(atomB)])
 
          if resonances in existing: # set lookup quick
            continue
 
          restraint = makeRestraint(targetValue=target,
                                    upperLimit=upper,
                                    lowerLimit=lower)
 
          restraint.newHBondConstraintItem(resonances=resonances)
      
      self.tabbedFrame.select(0)                                                  
      self.update_idletasks()
 
  
  def selectHelix(self, obj, row, col):
  
    self.helix = obj
  
  def selectStrands(self, obj, row, col):
  
    self.strand = obj
  
  def updateNonCovalentAfter(self, obj):
  
    if self.waitingNc:
      return
      
    self.waitingNc = True
    self.after_idle(self.updateNonCovalent)  
  
  def updateNonCovalent(self):

    textMatrixH = []
    textMatrixS = []
    
    sheetObjs = []
    if self.molSystem:
      
      for i, helix in enumerate(self.helices):
        atomO = helix.getStartAtom()
        atomH = helix.getEndAtom()
        
        residueA = atomO.residue
        residueB = atomH.residue
        
        datum = [i+1,
                 '%d%s%s' % (residueA.seqCode, getResidueCode(residueA), atomO.name),
                 '%d%s%s' % (residueB.seqCode, getResidueCode(residueB), atomH.name),
                 helix.getTypeString(),
                 ]
        
        textMatrixH.append(datum)
     
      for i, sheet in enumerate(self.strands):
        
        atomO1 = sheet.getStartAtomO()
        atomO2 = sheet.getEndtAtomO()
        atomH1 = sheet.getStartAtomH()
        atomH2 = sheet.getEndtAtomH()
        
        residueA1 = atomO1.residue
        residueB1 = atomH1.residue
        residueA2 = atomO2.residue
        residueB2 = atomH2.residue
                 
        datum = [i+1,
                 '%d%s%s' % (residueA1.seqCode, getResidueCode(residueA1), atomO1.name),
                 '%d%s%s' % (residueA2.seqCode, getResidueCode(residueA2), atomO2.name),
                 '%d%s%s' % (residueB1.seqCode, getResidueCode(residueB1), atomH1.name),
                 '%d%s%s' % (residueB2.seqCode, getResidueCode(residueB2), atomH2.name),
                 sheet.getTypeString(),
                ]
        
        delta = 'TBD'
        sheetObjs.append((residueA1, residueB1, residueA2, residueB2, delta))
        
        textMatrixS.append(datum)
        
     
    self.helixTable.update(textMatrix=textMatrixH, objectList=self.helices)
    self.sheetTable.update(textMatrix=textMatrixS, objectList=self.strands)
       
    self.waitingNc = False
  
  def addHelix(self):
  
    if self.molSystem:
      self.helices.append(Helix(self.molSystem))
      self.helices.sort()
      self.updateNonCovalent()
  
  def removeHelix(self):
  
    if self.helix and self.molSystem:
      self.helices.remove(self.helix)
      self.helix.delete()
      self.helix = None
      self.updateNonCovalent()
  
  def addStrands(self):
  
    if self.molSystem:
      self.strands.append(Sheet(self.molSystem))
      self.strands.sort()
      self.updateNonCovalent()
  
  def removeStrands(self):
  
    if self.strand and self.molSystem:
      self.strands.remove(self.strand)
      self.strand.delete()
      self.strand = None
      self.updateNonCovalent()
  
  def setHelixStart(self, residue):
  
    if self.helix and residue:
      self.helix.setStart(residue)   
      
  def setHelixEnd(self, residue):
      
    if self.helix and residue:
      self.helix.setEnd(residue)

  def setHelixType(self, delta):
      
    if self.helix:
      self.helix.setRegister(delta)
  
  def setStrand1Start(self, residue):

    if self.strand and residue:
      self.strand.setStartO(residue)        
  
  def setStrand1End(self, residue):
  
    if self.strand and residue:
      self.strand.setEndO(residue)        

  def setStrand2Start(self, residue):
  
    if self.strand and residue:
      self.strand.setStartH(residue)
  
  def setStrand2End(self, residue):
  
    if self.strand and residue:
      self.strand.setEndH(residue)

  def setStrandType(self, delta):
  
    if self.strand:
      self.strand.setRegister(delta)
  
  def _isValidPosition(self, residue, exclude=None):
  
    for helix in self.helices:
      if exclude is helix:
        # Don't know whether setting start or end
        # Dealt with in Helix obj
        
        continue
      
      if residue in helix.residues:
        return False
    
    for sheet in self.strands:
      
      # Could fall off the ends
      # Could be too small
      # Dealt with in Sheet obj
      
      if exclude is sheet:
        continue
        
      if residue in sheet.residuesO:
        return False
      
      if residue in sheet.residuesH:
        return False    
    
    return True    
        
  def getHelixStart(self, helix):
  
    self.getElementPosition(helix, helix.startResidue, self.helixStartPulldown)
   
  def getHelixEnd(self, helix):
      
    self.getElementPosition(helix, helix.endResidue, self.helixEndPulldown)
  
  def getElementPosition(self, element, residue0, pulldown):

    molSystem = self.molSystem
    names = ['%s%s' % (residue0.seqCode, getResidueCode(residue0)),]
    residues = [residue0,]
    cats = [None,]
    isValid = self._isValidPosition
    
    for chain in molSystem.sortedChains():
      chCode = chain.code
    
      for residue in chain.sortedResidues():
         
        if not isValid(residue, exclude=element):
          continue

        seqCode = residue.seqCode
        x = 10*(seqCode//10)
        
        cats.append('%s %d-%d' % (chCode, x,x+9))
        names.append('%s%s' % (seqCode, getResidueCode(residue)))
        residues.append(residue)

    pulldown.setup(names, residues, 0, categories=cats)
    
    
  def getHelixType(self, helix):
  
    objs = [4,5,3]
    residueA, residueB, delta = helix
    index = objs.index(delta)
  
    self.helixTypePulldown.setup(HELIX_TYPES, objs, index)
    
  def getStrand1Start(self, strand):
   
    self.getElementPosition(strand, strand.startResidueA, self.strandStart1Pulldown)
   
  def getStrand1End(self, strand):
   
    self.getElementPosition(strand, strand.endResidueA, self.strandEnd1Pulldown)
    
  def getStrand2Start(self, strand):
  
    self.getElementPosition(strand, strand.startResidueB, self.strandStart2Pulldown)
   
  def getStrand2End(self, strand):
  
    self.getElementPosition(strand, strand.endResidueB, self.strandEnd2Pulldown)
    
  def getStrandType(self, strand):
   
    if strand[-1] < 0:
      index = 0
    else:
      index = 1
 
    self.strandTypePulldown.setup(SHEET_TYPES, [-1,1], index)
    
  # # # # # # # # Nucleic acids # # # # # # # #





# # # # # SS Classes # # # # # 

class Duplex:

  err = 'MakeHbondREstraints Duplex Class Error: '

  def __init__(self, molSystem, startResidueA=None, endResidueA=None,
               startResidueB=None, endResidueB=None):
  
    self.molSystem = molSystem
    
    self._fragmentA = startResidueA.chainFragment
    self._fragmentB = startResidueB.chainFragment
    
    self.residues1 = None
    self.residues2 = None
    
    self.nonCovalentBonds = []
    
    if startResidueA and endResidueA and startResidueB and endResidueB:
      self._setRange(startResidueA, endResidueA, startResidueB, endResidueB)

  def __cmp__(self, other):
  
    if isinstance(Sheet, other):
      minA = min(self.startResidueA.ccpCode, self.startResidueB.ccpCode)
      minB = min(other.startResidueA.ccpCode, other.startResidueB.ccpCode)
    
      return cmp(minA, minB)
    else:
      return cmp(self, other)

  def delete(self):

    deleteBasePairHbonds(self.molSystem, self.residues1+self.residues2)
    del self
  
  def _setRange(self, startResidueA, endResidueA, startResidueB, endResidueB, adjust=1):
    
    self._fragmentA = startResidueA.chainFragment
    self._fragmentB = startResidueB.chainFragment
    
    if endResidueA.chainFragment is not self._fragmentA:
      raise Exception(self.err + 'Mismatching chain fragments on first strand')

    if endResidueB.chainFragment is not self._fragmentB:
      raise Exception(self.err + 'Mismatching chain fragments on second strand')
    
    if startResidueA.molType not in NUCLEIC_MOLTYPES:
      raise Exception(self.err + 'Start residue  on first strand must be in a nucleic acid')
    
    if endResidueA.molType  not in NUCLEIC_MOLTYPES:
      raise Exception(self.err + 'End residue on first strand must be in a nucleic acid')
    
    if startResidueB.molType  not in NUCLEIC_MOLTYPES:
      raise Exception(self.err + 'Start residue on second strand must be in a nucleic acid')
    
    if endResidueB.molType  not in NUCLEIC_MOLTYPES:
      raise Exception(self.err + 'End residue on second strand must be in a nucleic acid')
    
    fragresidues1 = list(self._fragmentA.residues)
    indexA1 = fragresidues1.index(self.startResidueA)
    indexB1 = fragresidues1.index(self.endResidueA)
  
    fragresidues2 = list(self._fragmentB.residues)
    indexA2 = fragresidues2.index(self.startResidueB)
    indexB2 = fragresidues2.index(self.endResidueB)
  
    if indexA1 > indexB1:
      indexA1, indexB1 = indexB1, indexA1
      self.startResidueA, self.endResidueA = self.endResidueA, self.startResidueA     

    if indexA2 < indexB2:
      indexA2, indexB2 = indexB2, indexA2
      self.startResidueB, self.endResidueB = self.endResidueB, self.startResidueB
    
    delta1 = indexB1 - indexA1
    delta2 = indexA2 - indexB2
    
    if delta1 != delta2:
      if adjust == 1: # Set 1 start 
        indexA2 = min(indexB2+delta1, len(fragresidues2)-1)
        startResidueB = fragresidues2[indexA2]
        
      elif adjust == 2: # Set 1 end
        indexB2 = max(0, indexA2-delta1)
        endResidueB = fragresidues2[indexB2]
      
      elif adjust == 3: # Set 2 start
        indexA1 = max(0, indexB1-delta2)
        startResidueA = fragresidues1[indexA1]
      
      else: # Set 2 end
        indexB1 = min(indexA1+delta2, len(fragresidues1)-1)
        endResidueA = fragresidues1[indexB1]
      
    if abs(indexB1-indexA1) < 1:
      showWarning(self.err, 'Duplex too short')
      return
      
    self.residues1 = [fragresidues1[indexA1:indexB1+1]]
    self.residues2 = [fragresidues1[indexB2:indexA2+1]]
    self.residues2.reverse()
    
    self.startResidueA = startResidueA
    self.endResidueA = endResidueA
    self.startResidueB = startResidueB
    self.endResidueB = endResidueB
    
    self._makeHbonds()
  
  def setStart1(self, residue):
  
    if self.endResidueA and self.startResidueB and self.endResidueB:
      self._setRange(residue, self.endResidueA, self.startResidueB, self.endResidueB, 1)
  
  def setEnd1(self, residue):
  
    if self.startResidueA and self.startResidueB and self.endResidueB:
      self._setRange(self.startResidueA, residue, self.startResidueB, self.endResidueB, 2)

  
  def setStart2(self, residue):
  
    if self.startResidueA and self.endResidueA and self.endResidueB:
      self._setRange(self.startResidueA, self.endResidueA, residue, self.endResidueB, 3)

  
  def setEnd2(self, residue):
  
    if self.startResidueA and self.endResidueA and self.startResidueB:
      self._setRange(self.startResidueA, self.endResidueA, self.startResidueB, residue, 4)
    
  def _makeHbonds(self):
  
    molSystem = self.molSystem
    deleteBackboneHbonds(molSystem, self.residues1+self.residues2)
    
    atomPairs = getProteinBackboneHbonds(molSystem)
    fragmentA = list(self._fragmentA.residues)
    fragmentB = list(self._fragmentB.residues)
 
    residueA = self.startResidueA
    residueB = self.startResidueB
    lastResidue = self.endResidueA
    
    self.nonCovalentBonds = []
    ncBondsAppend = self.nonCovalentBonds.append
    
    while residueA:
      indexA = fragmentA.index(residueA)
      indexB = fragmentB.index(residueB)
   
      ccpCodeA = residueA.ccpCode
      ccpCodeB = residueB.ccpCode
      ccpCodes = frozenset([ccpCodeA ,ccpCodeB])
   
      hbonds = NUCLEIC_HBONDS.get(ccpCodes, [])
      for hbondDict in hbonds:
        
        atomA = residueA.findFirstAtom(name=hbondDict[ccpCodeA])
        atomB = residueB.findFirstAtom(name=hbondDict[ccpCodeB])
        
        # TBD: wb104: 11 Feb 2011: added below; not sure if this is the intended order or whether O is before H (both seem to occur in NUCLEIC_HBONDS)
        atomH = atomA 
        atomO = atomB 

        if atomPairs.get(atomO) is not atomH:
          ncBond = molSystem.newNonCovalentBond(atoms=(atomO,atomH))
          ncBondsAppend(ncBond)
          atomPairs[atomO] = atomH
          atomPairs[atomH] = atomO
 
      if indexA+1 >= len(fragmentA):
        break
 
      elif 0 <= indexB-1:
        residueA = fragmentA[indexA+1]
        residueB = fragmentB[indexB-1]
 
      else:
        break
      
      if residueA is lastResidue:
        residueA = None  

class Sheet:
  
  err = 'MakeHbondREstraints Sheet Class Error: '

  def __init__(self, molSystem, startResidueA=None, endResidueA=None,
               startResidueB=None, endResidueB=None, register=-1):

    self.molSystem = molSystem
    self.register = register
    
    self._fragmentA = startResidueA.chainFragment
    self._fragmentB = startResidueB.chainFragment
    
    self.residuesO = None
    self.residuesH = None
    
    self.nonCovalentBonds = []
    
    if startResidueA and endResidueA and startResidueB and endResidueB:
      self._setRange(startResidueA, endResidueA, startResidueB, endResidueB)

  def __cmp__(self, other):
  
    if isinstance(Sheet, other):
      minA = min(self.startResidueA.ccpCode, self.startResidueB.ccpCode)
      minB = min(other.startResidueA.ccpCode, other.startResidueB.ccpCode)
    
      return cmp(minA, minB)
    else:
      return cmp(self, other)

  def delete(self):

    deleteBackboneHbonds(self.molSystem, self.residuesO+self.residuesH)
    del self

  def getTypeString(self):
  
    if self.register < 0:
      return SHEET_TYPES[0]
    else:
      return SHEET_TYPES[1]
 
   
  def _setRange(self, startResidueA, endResidueA, startResidueB, endResidueB, adjust=1):
    
    self._fragmentA = startResidueA.chainFragment
    self._fragmentB = startResidueB.chainFragment
    
    if endResidueA.chainFragment is not self._fragmentA:
      raise Exception(self.err + 'Mismatching chain fragments on first strand')

    if endResidueB.chainFragment is not self._fragmentB:
      raise Exception(self.err + 'Mismatching chain fragments on second strand')
    
    if startResidueA.molType != PROTEIN_MOLTYPE:
      raise Exception(self.err + 'Start residue  on first strand must be in a polypeptide')
    
    if endResidueA.molType != PROTEIN_MOLTYPE:
      raise Exception(self.err + 'End residue on first strand must be in a polypeptide')
    
    if startResidueB.molType != PROTEIN_MOLTYPE:
      raise Exception(self.err + 'Start residue on second strand must be in a polypeptide')
    
    if endResidueB.molType != PROTEIN_MOLTYPE:
      raise Exception(self.err + 'End residue on second strand must be in a polypeptide')
    
    if not startResidueA.findFirstAtom(name='O'):
      showWarning(self.err, 'Start residue on first strand has no carbonyl O')
      return

    if not endResidueA.findFirstAtom(name='O'):
      showWarning(self.err, 'End residue on first strand has no carbonyl O')
      return

    if not startResidueB.findFirstAtom(name='H'):
      showWarning(self.err, 'Start residue on second strand has no amide H')
      return

    if not endResidueB.findFirstAtom(name='H'):
      showWarning(self.err, 'End residue on second strand has no amide H')
      return

    if startResidueB.ccpCode in PROLYL:
      showWarning(self.err, 'Start residue on second strand cannot be prolyl')
      return

    if endResidueB.ccpCode in PROLYL:
      showWarning(self.err, 'End residue on second strand cannot be prolyl')
      return
    
    fragResiduesO = list(self._fragmentA.residues)
    indexA1 = fragResiduesO.index(self.startResidueA)
    indexB1 = fragResiduesO.index(self.endResidueA)
  
    fragResiduesH = list(self._fragmentB.residues)
    indexA2 = fragResiduesH.index(self.startResidueB)
    indexB2 = fragResiduesH.index(self.endResidueB)
  
    if indexA1 > indexB1:
      indexA1, indexB1 = indexB1, indexA1
      self.startResidueA, self.endResidueA = self.endResidueA, self.startResidueA     

    if indexA2 > indexB2:
      indexA2, indexB2 = indexB2, indexA2
      self.startResidueB, self.endResidueB = self.endResidueB, self.startResidueB
    
    deltaO = indexB1 - indexA1
    deltaH = abs(indexB2-indexA2)
    
    if deltaO != deltaH:
      if adjust == 1: # Set O start 
        if self.register < 0:
          indexA2 = min(indexB2+deltaO, len(fragResiduesH)-1)
        else:
          indexA2 = max(0, indexB2-deltaO)
        
        startResidueB = fragResiduesH[indexA2]
        
      elif adjust == 2: # Set O end
        if self.register < 0:
          indexB2 = max(0, indexA2-deltaO)
        else:
          indexB2 = min(indexA2+deltaO, len(fragResiduesH)-1)
        
        endResidueB = fragResiduesH[indexB2]
      
      elif adjust == 3: # Set H start
        indexA1 = max(0, indexB1-deltaH)
        startResidueA = fragResiduesO[indexA1]
      
      else: # Set H end
        indexB1 = min(indexA1+deltaH, len(fragResiduesO)-1)
        endResidueA = fragResiduesO[indexB1]
      
    if abs(indexB1-indexA1) < 2:
      showWarning(self.err, 'Strand too short')
      return
      
    self.residuesO = [fragResiduesO[indexA1:indexB1+1]]
     
    self.residuesH = []
    if indexA2 < indexB2: # Parallel
      self.residuesH =[fragResiduesO[indexA2:indexB2+1]]
    
    else: # Anti-parallel
      self.residuesH =[fragResiduesO[indexB2:indexA2+1]]
      self.residuesH.reverse()
     
    
    self.startResidueA = startResidueA
    self.endResidueA = endResidueA
    self.startResidueB = startResidueB
    self.endResidueB = endResidueB
    
    self._makeHbonds()
  
  def setStartO(self, residue):
  
    if self.endResidueA and self.startResidueB and self.endResidueB:
      self._setRange(residue, self.endResidueA, self.startResidueB, self.endResidueB, 1)
  
  def setEndO(self, residue):
  
    if self.startResidueA and self.startResidueB and self.endResidueB:
      self._setRange(self.startResidueA, residue, self.startResidueB, self.endResidueB, 2)

  
  def setStartH(self, residue):
  
    if self.startResidueA and self.endResidueA and self.endResidueB:
      self._setRange(self.startResidueA, self.endResidueA, residue, self.endResidueB, 3)

  
  def setEndH(self, residue):
  
    if self.startResidueA and self.endResidueA and self.startResidueB:
      self._setRange(self.startResidueA, self.endResidueA, self.startResidueB, residue, 4)

  
  def setRegister(self, register):
  
    if register != self.register:
      self.register = register
      
      fragmentH = list(self._fragmentB.residues)
      indexA = fragmentH.index(self.startResidueB)
      indexB = fragmentH.index(self.endResidueB)
      
      if (register < 0) and (indexA < indexB):
        self.startResidueB, self.endResidueB = self.endResidueB, self.startResidueB
      
      elif indexA > indexB:
        self.startResidueB, self.endResidueB = self.endResidueB, self.startResidueB
       
      self._makeHbonds()   
    
  def getStartAtomO(self, residue):
  
    return self.startResidueA.findFirstAtom(name='O')

  def getEndAtomO(self, residue):
  
    if len(self.residuesO) % 2 == 0:
      return self.endResidueA.findFirstAtom(name='H')
    
    else:
      return self.endResidueA.findFirstAtom(name='O')
    
  def getStartAtomH(self, residue):
  
    return self.startResidueB.findFirstAtom(name='O')

  def getEndAtomH(self, residue):
  
    if len(self.residuesO) % 2 == 0:
       return self.endResidueB.findFirstAtom(name='O')
   
    else:
      return self.endResidueB.findFirstAtom(name='H')
    
  def _makeHbonds(self):
  
    molSystem = self.molSystem
    deleteBackboneHbonds(molSystem, self.residuesO+self.residuesH)
    
    atomPairs = getProteinBackboneHbonds(molSystem)
    fragmentA = list(self._fragmentA.residues)
    fragmentB = list(self._fragmentB.residues)
 
    jump = self.register
 
    residueA = self.startResidueA
    residueB = self.startResidueB
    lastResidue = self.endResidueA
    
    self.nonCovalentBonds = []
    ncBondsAppend = self.nonCovalentBonds.append
    
    i = 0
    while residueA:
      i += 1
      if i % 2 == 0:
        atomO = residueB.findFirstAtom(name='O')
        atomH = residueA.findFirstAtom(name='H')
      
      else: # Setup to always start on the O side
        atomO = residueA.findFirstAtom(name='O')
        atomH = residueB.findFirstAtom(name='H')
 
      if not atomO:
        showWarning(self.err,'Residue missing carbonyl oxygen')
        return
 
      indexA = fragmentA.index(residueA)
      indexB = fragmentB.index(residueB)
   
      if atomH: # Not PRO, HYP
        if atomPairs.get(atomO) is not atomH:
          ncBond = molSystem.newNonCovalentBond(atoms=(atomO,atomH))
          ncBondsAppend(ncBond)
          atomPairs[atomO] = atomH
          atomPairs[atomH] = atomO
 
      if indexA+1 >= len(fragmentA):
        break
 
      elif 0 <= indexB+jump < len(fragmentB):
        residueA = fragmentA[indexA+1]
        residueB = fragmentB[indexB+jump]
 
      else:
        break
      
      if residueA is lastResidue:
        residueA = None  
        
class Helix:
  
  err = 'MakeHbondREstraints Helix Class Error: '

  def __init__(self, molSystem, startResidue=None, endResidue=None, register=4):

    self.molSystem = molSystem
    self.register = register
    self._fragment = startResidue.chainFragment
    self.residues = None
    self.nonCovalentBonds = []
    
    if startResidue and endResidue:
      self._setRange(startResidue, endResidue)

  def __cmp__(self, other):
  
    if isinstance(Helix, other):
      return cmp(self.startResidue.ccpCode, other.startResidue.ccpCode)
    else:
      return cmp(self, other)

  def delete(self):

    deleteBackboneHbonds(self.molSystem, self.residues)
    del self

  def getTypeString(self):
  
    return HELIX_DELTAS[self.register]
  
  def _setRange(self, startResidue, endResidue):
    
    self._fragment = startResidue.chainFragment
    
    if endResidue.chainFragment is not self._fragment:
      raise Exception(self.err + 'Mismatching chain fragments')
    
    if startResidue.molType != PROTEIN_MOLTYPE:
      raise Exception(self.err + 'Start residue must be in a polypeptide')
    
    if endResidue.molType != PROTEIN_MOLTYPE:
      raise Exception(self.err + 'End residue must be in a polypeptide')
    
    if not startResidue.findFirstAtom(name='O'):
      showWarning(self.err, 'Start residue has no carbonyl O')
      return

    if not endResidue.findFirstAtom(name='H'):
      showWarning(self.err, 'End residue has no amide H')
      return

    if startResidue.ccpCode in PROLYL:
      showWarning(self.err, 'Start residue cannot be prolyl')
      return

    if endResidue.ccpCode in PROLYL:
      showWarning(self.err, 'End residue cannot be prolyl')
      return
    
    fragResidues = list(self._fragment.residues)
    indexA = fragResidues.index(startResidue)
    indexB = fragResidues.index(endResidue)
    
    if abs(indexA-indexB) < (self.register+1):
      showWarning(self.err, 'Too short a span for a helix')
      return
     
    if indexA > indexB:
      indexA, indexB = indexB, indexA
      startResidue, endResidue = endResidue, startResidue
    
    self.residues = self._fragment[indexA:indexB+1-self.register]
    
    self.startResidue = startResidue
    self.endResidue = endResidue
    self._makeHbonds()
  
  def setStart(self, residue):
  
    if self.endResidue:
      self._setRange(residue, self.endResidue)
  
  def setEnd(self, residue):
  
    if self.startResidue:
      self._setRange(self.startResidue, residue)
  
  def setRegister(self, register):
  
    if self.register != register:
      self.register = register
      self._makeHbonds()   
  
  def getStartAtom(self, residue):
  
    return self.startResidue.findFirstAtom(name='O')

  def getEndAtom(self, residue):
  
    return self.endResidue.findFirstAtom(name='H')
  
  def _makeHbonds(self):

    molSystem = self.molSystem
    register = self.register
    fragment = list(self._fragment.residues)    
    residue = self.startResidue 
    endResidue = self.endResidue
    partnerResidue = None
    self.nonCovalentBonds = []
    ncBondsAppend = self.nonCovalentBonds.append

    deleteBackboneHbonds(molSystem, self.residues)
    atomPairs = getProteinBackboneHbonds(molSystem)
    
    while residue:
      atomO = residue.findFirstAtom(name='O')

      if not atomO:
        showWarning(self.err, 'Residue missing carbonyl oxygen')
        return
 
      index = fragment.index(residue)
      partnerResidue = fragment[index+register]
      atomH = partnerResidue.findFirstAtom(name='H')
 
      if atomH: # Not PRO, HYP
        if atomPairs.get(atomO) is not atomH:
          ncBond = molSystem.newNonCovalentBond(atoms=(atomO,atomH))
          ncBondsAppend(ncBond)
          atomPairs[atomO] = atomH
          atomPairs[atomH] = atomO
 
      if index+1 < len(fragment):
        residue = fragment[index+1]
      else:
        break
        
      if partnerResidue is endResidue:
        residue = None
        
# # # # # # # Utility # # # # # #

def deleteBasePairHbonds(molSystem, residues):

  residues = set(residues)

  for residueA in residues:
    atomO = residueA.findFirstAtom(name='O')
    
    if atomO:
      ncBonds = list(atomO.nonCovalentBonds)
      for ncBond in ncBonds:
        atomA, atomB = ncBond.atoms

        if atomA is atomO:
          atomH = atomB
        else:
          atomH = atomA
          
        if atomH.name != 'H':
          continue
        
        residueB = atomH.residue  
        if residueB.molType != PROTEIN_MOLTYPE:
          continue
        
        if residueB in residues:
          ncBond.delete()

def getBasePairHbonds(molSystem):

  pass

def deleteBackboneHbonds(molSystem, residues):

  residues = set(residues)

  for residueA in residues:
    atomO = residueA.findFirstAtom(name='O')
    
    if atomO:
      ncBonds = list(atomO.nonCovalentBonds)
      for ncBond in ncBonds:
        atomA, atomB = ncBond.atoms

        if atomA is atomO:
          atomH = atomB
        else:
          atomH = atomA
          
        if atomH.name != 'H':
          continue
        
        residueB = atomH.residue  
        if residueB.molType != PROTEIN_MOLTYPE:
          continue
        
        if residueB in residues:
          ncBond.delete()
                  

def getProteinBackboneHbonds(molSystem):

  atomPairs = {}
  
  for ncBond in molSystem.nonCovalentBonds: 
    
    atom1, atom2 = ncBond.atoms
    
    atomNames = [atom1.name, atom2.name]
    atomNames.sort()
    
    if atomNames != ['H','O']:
      continue
          
    molType1 = atom1.residue.molType
    if molType1 != PROTEIN_MOLTYPE:
      continue
    
    molType2 = atom2.residue.molType
    if molType2 != PROTEIN_MOLTYPE:
      continue
     
    atomPairs[atom1] = atom2
    atomPairs[atom2] = atom1
      
  return atomPairs
      
def getHelicalHbondSegments(molSystem, minHbonds=2):


  atomPairs = getProteinBackboneHbonds(molSystem)
  segmentsLong = []

  for chain in molSystem.sortedChains():
    for fragment in chain.sortedChainFragments():
      if fragment.molType != PROTEIN_MOLTYPE:
        continue
      
      if not fragment.isLinearPolymer:
        continue
    
      residues = list(fragment.residues)
      segments = []
      helixBonds = []
      for i, residue in enumerate(residues):
        atomO = residue.findFirstAtom(name='O')
        atomH = atomPairs.get(atomO)

        if not (atomO and atomH):
          continue
          
        residueB = atomH.residue
        if residueB in residues:
          posO = residues.index(residue)
          posH = residues.index(residueB)
          delta = posH - posO
          
          if not (2 < delta < 6):
            continue
        
	  helixBonds.append([delta, posO, posH, atomO, atomH])

      if not helixBonds:
        continue

      helixBonds.sort()
      for helixBond in helixBonds:
        if segments:
	  lastSeg = segments[-1]
	  delta1, posO1, posH1, atomO1, atomH1 = lastSeg
          delta2, posO2, posH2, atomO2, atomH2 = helixBond
  
          if delta1 != delta2:
	    segments.append(helixBond)
	  
	  elif posO2 != posO1+1:
	    segments.append(helixBond)
	  
	  else:
	    lastSeg[2] = posH2
	    lastSeg[4] = atomH2
        
	else:
	  segments.append(helixBond)
      
	      
  
      for delta, start, end, atomO, atomH in segments:
        if end-start >= delta+minHbonds-1:
          segmentsLong.append(Helix(molSystem, atomO.residue, atomH.residue, delta))

  return segmentsLong

def getPairedHbondSheets(molSystem):

  sheets = []

  atomPairs = getProteinBackboneHbonds(molSystem)

  clusters = []
  for chain in molSystem.sortedChains():
    for fragment in chain.sortedChainFragments():
      if fragment.molType != PROTEIN_MOLTYPE:
        continue
      
      if not fragment.isLinearPolymer:
        continue
    
      residues = list(fragment.residues)
      for i, residue in enumerate(residues):
        atomO = residue.findFirstAtom(name='O')
	
	atomH = atomPairs.get(atomO)
        if not (atomO and atomH):
          continue

        residueB = atomH.residue
        fragmentB = residueB.chainFragment
        if fragmentB.molType != PROTEIN_MOLTYPE:
          continue
	
	residuesB = list(fragmentB.residues)
	posO = residues.index(residue)
        posH = residuesB.index(residueB)
  
        cluster = (fragment, fragmentB, posO, posH, atomO, atomH, None, None)
        clusters.append([cluster,])
	
  changes = 0
  while changes:
    changes = 0
    
    n = len(clusters)
    for i in range(n):
      clusterI = clusters[i]
      
      if not clusterI:
        continue
    
      fragmentO1, fragmentH1, posO1, posH1, atomO1, atomH1, delta1 = clusterI[-1]
    
      for j in range(n):
        if i == j:
          continue
          
        clusterJ = clusters[j]
        if not clusterJ:
          continue
          
        fragmentO2, fragmentH2, posO2, posH2, atomO2, atomH2, delta2 = clusterJ[0]
        
        if (delta1 and delta1) and (delta1 != delta2):
          continue
        
        if fragmentO1 is not fragmentH2:
          continue
        
        if fragmentH1 is not fragmentO2:
          continue
        
        deltaA = posH2-posO1
        if abs(deltaA) != 1:
          continue
        
        deltaB = posH1-posO2
        if abs(deltaB) != 1:
          continue
          
        if deltaA == deltaB : # Anti
          delta = -1
           
        else: # Para
          delta = 1
        
        if delta1 and (delta1 != delta):
          continue
          
        if delta2 and (delta2 != delta):
          continue
        
        # Link
        if not delta1:
          clusterI[0][-1] = delta

        if not delta2:
          clusterJ[0][-1] = delta
           
        changes += 1
        clusters[i] = clusterI + clusterJ
        clusters[j] = []
  
    clusters = [c for c in clusters if c]
  
  for cluster in clusters:
    nRes = len(cluster)
  
    if nRes < 2:
      continue
  
    delta = cluster[0][-1]
    startA = cluster[0][4] # Start at O
    startB = cluster[0][5] # Start at H
    
    if nRes % 2 == 0:
      endA = cluster[-1][5] # Even -> Ends in H
      endB = cluster[-1][4] # Even -> Ends in O
    else:
      endA = cluster[-1][4] # Odd -> ends in O
      endB = cluster[-1][5] # Odd -> ends in H

    sheet = Sheet(molSystem, startA.residue, endA.residue,
                  startB.residue, endB.residue, delta)
                  
    sheets.append(sheet)
  
  return sheets



"""
======================COPYRIGHT/LICENSE START==========================

BrowseAtoms.py: Part of the CcpNmr Analysis program

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

from memops.general import Implementation

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.RadioButtons import RadioButtons
from memops.gui.CheckButton import CheckButton
from memops.gui.Color import scaleColor
from memops.gui.PartitionedSelector import PartitionedSelector
from memops.gui.ProgressBar import ProgressBar
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showOkCancel
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.AssignmentBasic import deassignResonance, getShiftLists
from ccpnmr.analysis.core.MoleculeBasic import greekSortAtomNames, getResidueMapping, getUnicodeAtomName
from ccpnmr.analysis.core.MoleculeBasic import getResidueCode, getResidueObservableAtoms
from ccpnmr.analysis.core.MoleculeBasic import makeResidueLocTag, makeResidueAtomSetsEquivalent
from ccpnmr.analysis.core.MoleculeBasic import PROTEIN_RESIDUE_CLASS_DICT, makeResidueAtomSetsNonEquivalent

FLIP_AROMATICS = ('Phe','Ptr','Tyr')

ATOM_COLOR_DICT = {'H':'#a0d0a0','N':'#a0a0d0',
                   'C':'#d0d0a0','O':'#d0a0a0',
                   'P':'#d0a0d0','S':'#d0c090',
                   'Se':'#e0b0a0',}
                   
MAP_TYPE_OPTS = {'nonstereo':1 ,'stereo':0,
                 'ambiguous':0, 'simple':1}

class BrowseAtomsPopup(BasePopup):
  """
  **Select Atoms and Residues for Assignment**

  This popup window allows the user to view and select specific atom, atom set
  and residue entities from any of the molecular systems within the CCPN
  project. Most importantly this allow the user to perform resonance assignment,
  of the signal phenomena observed in NMR spectra, to particular atomic groups,
  thus labelling the origin of a resonance. In a similar manner, spin systems
  (groups of resonances from the same residue) may be linked to residues. In
  Analysis such atom and residue assignment may be performed in a specific way,
  i.e. choosing a particular atom in a particular residue within the sequence.
  However, this system also allows the user to select atoms and  residues in
  order to specify the type of a resonance or spin system; to say what kind of
  atom or residue a resonance or spin system relates to without actually
  specifying which position in the molecular sequence. This type-only
  information is useful to reflect the current state of knowledge, particularly
  during early stages of assignment, without having to make notes outside of the
  CCPN system.

  This system has a second important function, to show the user the current
  assignment state for the atoms of the selected chemical element types,
  molecular chain (and bearing in mind any shift list) selection. The cells
  within the table will show a darker colour if they carry a resonance
  assignment (irrespective of whether the resonance is associated with peaks)
  and thus the user can quickly get a visual appreciation of the parts of the
  molecular chain that have been assigned and those that have yet to be assigned

  **Atom Table**

  In order to use this system the user must ensure that the various options have
  been chosen appropriate to what they wish to view. Generally in the first
  "Atom Table" tab the user selects which molecular chain in the upper left
  pulldown men and click the "Elements" buttons to specify which kinds of atom
  will be visible at this time. Naturally the selection of chemical elements
  will reflect the kind of operation or assignment that is being performed at
  the  time. Thus for example if the user is assigning a 15N HSQC NOESY
  experiment which shows only 15N and 1H resonances then only the [H] and [N]
  chemical element buttons need to be toggled on. If the user wants to reduce
  clutter and only look at a particular type of residue (amino acid, nucleotide,
  sugar) then the "Residue" pulldown can be used to display only a subset of
  the sequence. The adjacent "Nucleotides" option is only relevant for nucleic
  acids like DNA & RNA, and allow the user the option of showing each residue on
  a single line (long row), or have the sugar phosphate and aromatic bas
  components displayed independently on separate lines; minimising clutter and
  the need to scroll.

  There are a few functions at the bottom of the first tab that control
  assignment and link to other kinds of data entity in the CCPN project. The
  [Show ...] buttons allow the user to get a separate table of either the
  spectrum peaks or resonance phenomena that have been associated with the atoms
  or resonances that have been selected in the table above. Specifically the
  `Selected Peaks`_ and `Selected Resonances`_ popups are opened. This is useful
  to navigate from the molecule based information display into the spectrum
  display and to understand the links throughout the project. The [Remove Atom
  Assignments] button does as the name suggests and disconnects any resonances
  assigned to the selected atoms from those atoms. This merely removes the
  specific atom information for those resonances, it does not affect the linking
  of the resonances to any peaks. Lastly the [Set Ring Flip Equivalency] button
  is used on the selected residue for the rare occasions that aromatic rings (in
  residues like Phe & Tyr) display two distinct resonances for atoms that in
  most circumstances are indistinguishable in NMR spectra, because they are in
  very rapid exchange due to rotation of the aromatic ring around an axis of
  symmetry. In essence this makes Phe & Thr HD*, CD, HE* and CE separate into
  HD1, HD2, CD1, CD2, HE1, HE2, CE1 & CE2 selections. This operation is
  reversible but will naturally loose fidelity of peak assignments if a peak is
  assigned to only one side of the ring.

  **Options**
  
  The second tab gives the user further control over how atom information is
  displayed in the main table. The upper "Prochirals" section allows the choice
  of how to display stereochemical information for atom sets that are often  not
  distinguishable in stereospecific way in NMR spectra. For example a Ser residue
  has both HB2 and HB3 atoms in its side chain, but these can prochiral atoms
  cannot generally be distinguished from one another in NMR spectra. In such
  circumstances it is commonplace to use non-stereospecific assignment labels
  like "HBa" and "HBb" so that two different resonances, usually with different
  chemical shifts, can be identified in NMR spectra, without needing to specify
  any stereochemical information. In this example "HBa" and "HBb" resonance
  assignments will both be linked to both HB2 and HB3, but this is done in a way
  that is mutually exclusive; Analysis knows that one really only goes to HB2
  and the other to HB3, whichever way round that may be. If the  user has
  structural  information that allows a proper determination of which resonance
  goes with which stereochemical location, then the "Stereospecific" option can be
  set. This would make Ser HB2 and HB3 available as separate options to "HBa"
  and "HBb" which come from the "Non-stereospecific" option. The "Ambiguous"
  option is something of a remnant and is generally only used for situations
  when the atoms in the stereo group have exactly the same chemical shift. For
  Ser HB2 & HB3 this would provide the "HB*" option in the main table, however
  this is really just a pseudonym; assigning to this actually assigns to both
  the stereospecific HB2 and HB3 at the same time (not HBa, HBb). 

  The "Assignment Status" selection allows the number of atom options displayed
  in the main table to be reduced by only showing specific cells based upon what
  kind of assignment the underlying atoms have. For example the "unassigned"
  option may be selected to view more easily what has yet to be assigned. The
  "Shift List" pulldown controls which atom options are darkened, indicating
  resonance assignments, based upon the experiments that those assignments are
  made to and which set of conditions (and thus shift list) those assignments
  relate. By default this option is set to "<Ignore>", meaning that all
  resonance assignments are considered. Selecting a specific shift list will
  reduce the indication of assignment to only those peaks that use that shift
  list (via the experiment's connection).

  The lower "Isotope Labelling" section is used if the user has any experiments
  that were performed with selective spin-active isotope labelling, which would
  reduce the complement of atom sites that were deemed to be visible in spectra
  and thus available for assignment. Here, the user selects which type of isotope
  labelling scheme was used, and what level of isotopic incorporation (between
  0,0 and 1.0) must be present for an atom site to be visible.  Atom sites that
  are not deemed to be visible will not be displayed in the  main sequence + atom
  table.

  **Tips and Caveats**
  
  At some stage in the future it is planned that the Atom Table will become
  (optionally) context sensitive to the kind of atoms that *ought* to be visible
  in the specific spectrum being assigned at the time. For example if an HNCA
  experiment is assigned then only the H, N & CA atom options cane be displayed.

  .. _`Selected Peaks`: SelectedPeaksPopup.html
  .. _`Selected Resonances`: BrowseResonancesPopup.html

  """

  def __init__(self, parent, *args, **kw):

    self.nucleotideOpt   = 0
    self.chemElementDict = {}
    self.chemElements    = []
    self.residueASMLists = []
    self.chain           = None
    self.residue         = None
    self.atomSetMapping  = None
    self.ccpCode         = None
    self.shiftList       = None
    self.labellingScheme = None
    self.refExperiment   = None
    self.refresh         = False
    
    BasePopup.__init__(self, parent=parent, title="Molecule : Atom Browser", **kw)

  def body(self, guiParent):

    self.geometry('500x500')

    guiParent.expandGrid(0,0)
    
    options = ['Atom Table','Options',]
    self.tabbedFrame = TabbedFrame(guiParent, options=options,
                                   callback=self.toggleTab, grid=(0,0))
    frameA, frameB = self.tabbedFrame.frames
    frameA.expandGrid(1,0)
    frameB.expandGrid(4,0)

    frameT = self.tabbedFrame.sideFrame
    self.bottomButtons = UtilityButtonList(frameT, doClone=False, grid=(0,0),
                                           helpUrl=self.help_url, sticky='e')
    # Atom Table

    frame = Frame(frameA, grid=(0,0), sticky='ew')
    frame.expandGrid(None,5)

    tipText = 'Select the molecular system & chain to show atoms for'
    label = Label(frame, text='Chain:', grid=(0,0), sticky='e')
    self.chainPulldown = PulldownList(frame, self.changeChain,
                                      grid=(0,1), tipText=tipText)

    tipText = 'Select the type of residue to display atoms for, defaults to all types'
    label = Label(frame, text = 'Residue:', grid=(0,2))
    self.resPulldown = PulldownList(frame, self.changeResCode,
                                    grid=(0,3), tipText=tipText)

    tipText = 'Select how to display rows of base and sugar-phosphate atoms for DNA & RNA nucleotide residues'
    label= Label(frame, text = 'Nucleotides:', grid=(0,4))
    self.nucleotidePulldown = PulldownList(frame,
                                           callback=self.changeNucleotideOpt,
                                           texts=['Interlaced','Base','Sugar-P','Long Row'],
                                           objects=range(4),
                                           grid=(0,5), tipText=tipText,
                                           index=self.nucleotideOpt)

    tipText = 'The selection of which chemical elements to display atoms for'
    self.chemElementLabel = Label(frame, text=' Elements:', grid=(1,0))
    labels = ['H','N','C','P']
    self.elementSelect = PartitionedSelector(frame, labels=labels,
                           objects=labels, callback=self.toggleChemElement,
                           selected=[], grid=(1,1), gridSpan=(None,5),
                           tipText=tipText)

    tipTexts = ['Sequence number code',
                'Residue type code']
    self.scrolledMatrix = ScrolledMatrix(frameA, specialBg='black',
                                         initialCols=10, tipTexts=tipTexts,
                                         headingList=['#','Residue'],
                                         highlightType='cell', grid=(1,0),
                                         callback=self.selectCell)
    
    tipTexts = ['Toggle whether the symmetric aromatic ring of selected residue\nflips quickly on the NMR timescale',
                'Remove all resonance assignments from the selected atoms',
                'Show a table of all peaks assigned to the selected atoms',
                'Show a table of all resonances assigned to the selected atoms']
    texts    = ['Set Ring Flip\nEquivalency',
                'Remove Atom\nAssignments',
                'Show\nPeaks','Show\nResonances']
    #commands = [self.minimAtomMatrix,self.changeStereo,self.getInfo]
    commands = [self.setAromaticEquivalency,
                self.removeAssignments,
                self.showPeaks, self.showResonances]
    self.bottomButtons = ButtonList(frameA, texts=texts, tipTexts=tipTexts,
                                    commands=commands, grid=(2,0))

    # Options

    row = 0

    #label= Label(frameB, text = 'Experiment Type:', grid=(0,0))
    tipText = 'Restrict atom options to only those visible in a selected kind of NMR experiment'
    self.experimentTypePulldown = PulldownList(frameB, self.changeExperimentType,
                                               tipText=tipText)
    #self.experimentTypePulldown.grid(row=0,column=1,sticky='w')

    numbers = range(1,7)
    #label= Label(frameB, text = 'Num J coupling bonds:', grid=(0,2))
    tipText = 'Select the number of covalent bonds to consider J coupling when predicting observed atoms'
    self.jBondsSelect = PartitionedSelector(frameB, labels=numbers, objects=numbers,
                                            tipText=tipText, callback=self.toggleJbonds,
                                            selected=range(1,3))
    #self.jBondsSelect.grid(row=0,column=3,sticky='w')
    
    frame = LabelFrame(frameB, text='Prochirals', grid=(0,0), sticky='ew')
    frame.expandGrid(0,6)
    
    tipText = 'Whether to show distinct but stereoscopically ambiguous assignment options.\n' + \
              'For example you may not know if a Serine Hb is "Hb2" or "Hb3"' 
    labelNonstereo = Label(frame, text='Non-Stereospecific ', grid=(0,0))
    self.checkNonstereo = CheckButton(frame, callback=self.setNonstereo,
                                      grid=(0,1), tipText=tipText)

    tipText = 'Whether to show stereospecific atom options (IUPAC system) like Serine Hb2 & Hb3'
    labelStereo = Label(frame, text='Stereospecific ', grid=(0,2))
    self.checkStereo = CheckButton(frame, callback=self.setStereo,
                                   grid=(0,3), tipText=tipText)

    tipText = 'Whether to include multi-atom ambiguous options like Ser Hb*.\n' + \
              'Methyl groups are not affected by this option unless the whole group is prochiral.' 
    labelAmbiguous = Label(frame, text='Ambiguous', grid=(0,4))
    self.checkAmbiguous = CheckButton(frame, callback=self.setAmbiguous,
                                      grid=(0,5), tipText=tipText)
    

    frame = LabelFrame(frameB, text='Assignment Status', grid=(1,0), sticky='ew')
    frame.expandGrid(0,2)
    entries = ['Any ','Assigned ','Unassigned', 'Tentative']
    tipTexts = ['Show all atom options irrespective of assignments',
                'Show only atom options which are assigned to NMR resonances',
                'Show only atom options which have no resonance assignments',
                'Show only atom options that carry tentative resonance assignments']
    self.assignSelect = RadioButtons(frame, entries=entries, grid=(0,0),
                                     select_callback=self.selectAssign,
                                     gridSpan=(1,3), tipTexts=tipTexts)
    
    tipText = 'Select a specific shift list to determine an atoms assignment status or "<ignore>" to use any'
    label= Label(frame, text='Shift List:', grid=(1,0))
    self.shiftListPulldown = PulldownList(frame, self.changeShiftList,
                                          grid=(1,1), sticky='w', tipText=tipText)

    frame = LabelFrame(frameB, text='Isotope Filtering', grid=(2,0), sticky='ew')
    frame.expandGrid(0,4)
    
    tipText = 'Select which per-residue isotope labelling scheme/pattern to use'
    label= Label(frame, text = 'Isotope Labelling:', grid=(0,0))
    self.labellingSchemePulldown = PulldownList(frame, self.changeLabellingScheme,
                                                grid=(0,1), tipText=tipText)

    tipText = 'Specify what proportion of atoms must be isotopically labelled (spin active) to be considered as assignable'
    label= Label(frame, text = ' Min Isotope Fraction:', grid=(0,2))
    self.minFractionEntry = FloatEntry(frame, text=0.25, width=8,
                                       grid=(0,3), tipText=tipText)

    #

    self.curateNotifiers(self.registerNotify)
    
    self.updateShiftLists()
    self.updateChains()
    self.updateExperimentTypes()
    self.updateLabellingSchemes()
    self.checkNonstereo.set(True)
  
  def curateNotifiers(self, notifyFunc):

    notifyFunc(self.updateChains, 'ccp.molecule.MolSystem.Chain', 'delete')
    notifyFunc(self.updateChains, 'ccp.molecule.MolSystem.Chain', '__init__')
    notifyFunc(self.updateAfter, 'ccp.molecule.MolSystem.Residue', 'setSeqCode') 
    notifyFunc(self.updateAfter, 'ccpnmr.Analysis.AtomSetMapping', 'delete') 
    notifyFunc(self.updateAfter, 'ccpnmr.Analysis.AtomSetMapping', 'setResonanceSerials') 
    notifyFunc(self.updateExperimentTypes,  'ccp.nmr.Nmr.Experiment', 'setRefExperiment') 
    notifyFunc(self.updateLabellingSchemes, 'ccp.molecule.ChemCompLabel.LabelingScheme', '__init__') 
    notifyFunc(self.updateLabellingSchemes, 'ccp.molecule.ChemCompLabel.LabelingScheme', 'delete') 

    for func in ('__init__', 'delete') :
      notifyFunc(self.updateLabellingSchemes, 'ccp.molecule.ChemCompLabel.LabelingScheme', func)
      notifyFunc(self.updateLabellingSchemes, 'ccp.molecule.LabeledMolecule.LabeledMixture', func) 
    
    for func in ('__init__', 'delete', 'setWeight'): 
      notifyFunc(self.updateMoLabelFraction, 'ccp.molecule.LabeledMolecule.MolLabelFraction', func) 

    for func in ('__init__','delete','setPossibility'):
      notifyFunc(self.updateResidueProb, 'ccp.nmr.Nmr.ResidueProb', func) 
    
    for func in ('setAssignNames','addAssignName','removeAssignName'):
      notifyFunc(self.updateResonance, 'ccp.nmr.Nmr.Resonance', func) 

  def open(self):
  
    self.updateShiftLists()
    self.updateChains()
    self.updateExperimentTypes()
    self.updateLabellingSchemes()
    BasePopup.open(self)
  
  def toggleTab(self, index):
  
    if index == 0:
      self.updateAfter()
  
  def updateResonance(self, resonance):
  
    spinSystem = resonance.resonanceGroup
    if spinSystem and spinSystem.residueProbs:
      for residueProb in spinSystem.residueProbs:
        self.updateResidueProb(residueProb)
  
  def updateResidueProb(self, residueProb):
  
    if residueProb.possibility.chain is self.chain:
      self.updateAfter()
    
  def selectAssign(self, name):
  
    self.updateAfter()
    
  def changeExperimentType(self, refExperiment):
    
    if refExperiment is not self.refExperiment:
      self.refExperiment = refExperiment
      self.updateAfter()
  
  def updateExperimentTypes(self, notifyExperiment=None):
  
    names = []
    index = -1
    
    refExperiment  = self.refExperiment
    refExperiments = self.getRefExperiments()
    
    if refExperiments:
      names = ['<None>'] + [refExp.name for refExp in refExperiments[1:]]
      
      if refExperiment not in refExperiments:
        refExperiment = refExperiments[0]
      
      index = refExperiments.index(refExperiment)

    else:
      refExperiment = None  
    
    if refExperiment is not self.refExperiment:
      self.refExperiment = refExperiment
      self.updateAfter()
      
    self.experimentTypePulldown.setup(names, refExperiments, index)     
  
  def getRefExperiments(self):
  
    refExperiments = set()
    
    for experiment in self.nmrProject.experiments:
      refExperiment = experiment.refExperiment
      
      if refExperiment:
        refExperiments.add((refExperiment.name, refExperiment))

    refExperiments = list(refExperiments)
    refExperiments.sort()
    
    return [None,] + [x[1] for x in refExperiments]

  def changeLabellingScheme(self, scheme):
    
    if scheme is not self.labellingScheme:
      self.labellingScheme = scheme
      self.updateAfter()
      
  def updateMoLabelFraction(self, molLabelFraction=None):

    if molLabelFraction.labeledMixture is self.labellingScheme:
      self.updateLabellingSchemes()

  def updateLabellingSchemes(self, notifyScheme=None):
 
    names = []
    labellings = []
    index = 0
    
    labelling = self.labellingScheme
    
    for mixture in self.getLabelledMixtures():
      name = 'Sample: %s v%d' % (mixture.labeledMolecule.molecule.name, mixture.serial)
      
      labellings.append(mixture)
      names.append(name)
    
    for scheme in self.getLabellingSchemes():
      name = 'Scheme: %s' % (scheme.name)
       
      labellings.append(scheme)
      names.append(name)
      
    names.append('<None>')
    labellings.append(None)

    if labelling not in labellings:
      labelling = None
      index = len(labellings)-1
    else:  
      index = labellings.index(labelling)
      
    if labelling is not self.labellingScheme:
      self.labellingScheme = labelling  
      self.updateAfter()  
   
    self.labellingSchemePulldown.setup(names, labellings, index)  
  
  def getLabellingSchemes(self):
  
    return self.project.sortedLabelingSchemes()

  def getLabelledMixtures(self):
  
    mixtures = []
    if self.chain:
      molecule = self.chain.molecule
      labelledMolecules = self.project.findAllLabeledMolecules(molecule=molecule)
      
      for labelledMolecule in labelledMolecules:
        mixtures += labelledMolecule.sortedLabeledMixtures()
  
    return mixtures
 
  def showPeaks(self):
  
    if self.atomSetMapping and self.atomSetMapping.resonanceSerials:
    
      serials = self.atomSetMapping.resonanceSerials
      residue = self.atomSetMapping.residueMapping.residue
      resonances = []
      peakDict = {}
      for atom in residue.atoms:
        if atom.atomSet:
          for resonanceSet in atom.atomSet.resonanceSets:
            for resonance in resonanceSet.resonances:
              if resonance.serial in serials:
                resonances.append(resonance)
      
      for resonance in resonances:
        for contrib in resonance.peakDimContribs:
          peakDict[contrib.peakDim.peak] = 1
  
      self.parent.viewPeaks(peakDict.keys())

  def showResonances(self):
  
    resonances = []
    if self.atomSetMapping and self.atomSetMapping.resonanceSerials:
      serials = self.atomSetMapping.resonanceSerials
      residue = self.atomSetMapping.residueMapping.residue

      for atom in residue.atoms:
        if atom.atomSet:
          for resonanceSet in atom.atomSet.resonanceSets:
            for resonance in resonanceSet.resonances:
              if resonance.serial in serials:
                resonances.append(resonance)
  
    elif self.residue:
      for atom in self.residue.atoms:
        if atom.atomSet:
          for resonanceSet in atom.atomSet.resonanceSets:
            for resonance in resonanceSet.resonances:
              resonances.append(resonance)
              
    if resonances:
      self.parent.viewSelectedResonances(resonances, self.shiftList)

  def removeAssignments(self):
  
    if self.atomSetMapping and self.atomSetMapping.resonanceSerials:
      residue = self.atomSetMapping.residueMapping.residue
      ccpCode = getResidueCode(residue)
      text = '%d%s%s' % (residue.seqCode,ccpCode,self.getUnicodeAtomName(self.atomSetMapping, residue.molResidue.molType))
      if showOkCancel('Confirm','Are you sure? All assignments to %s will be removed.' % text, parent=self):
        serials    = self.atomSetMapping.resonanceSerials
        resonances = []
        for atom in residue.atoms:
          if atom.atomSet:
            for resonanceSet in atom.atomSet.resonanceSets:
              for resonance in resonanceSet.resonances:
                if resonance.serial in serials:
                  resonances.append(resonance)
 
        for resonance in resonances:
          deassignResonance(resonance)

  def setAromaticEquivalency(self):
  
    if self.residue and (self.residue.ccpCode in FLIP_AROMATICS ):
      atom = self.residue.findFirstAtom(name='CD1')
      if atom and atom.atomSet:
        if self.isAromaticDEAssigned(self.residue):
          if len(atom.atomSet.atoms) > 1:
            if showOkCancel('Confirm','Are you sure? Any current assignments will become non-stereospecific.', parent=self):
              makeResidueAtomSetsNonEquivalent(self.residue)
          else:
            if showOkCancel('Confirm','Are you sure? Any current assignments may be duplicated', parent=self):
              makeResidueAtomSetsEquivalent(self.residue)

        else:
          if len(atom.atomSet.atoms) > 1:
            makeResidueAtomSetsNonEquivalent(self.residue)
          else:
            makeResidueAtomSetsEquivalent(self.residue)
        
  def isAromaticDEAssigned(self, residue):
  
    for atom in residue.atoms:
      if atom.name in ('CD1','CE1','CD2','CE2','HD1','HE1','HD2','HE2',):
        if atom.atomSet and atom.atomSet.resonanceSets:
          return True
    
    return False

  def setNonstereo(self, boolean):
  
    MAP_TYPE_OPTS['nonstereo'] = boolean
    self.updateAfter()
  
  def setStereo(self,  boolean):

    MAP_TYPE_OPTS['stereo'] = boolean
    self.updateAfter()
  
  def setAmbiguous(self, boolean ):

    MAP_TYPE_OPTS['ambiguous'] = boolean
    self.updateAfter()

  def getShiftListNames(self, shiftLists):
    
    shiftListNames = []
    for shiftList in shiftLists:
      if not shiftList.name:
        shiftList.name = "ShiftList "+ str(shiftList.serial)
      shiftListNames.append(shiftList.name)

    return shiftListNames

  def updateShiftLists(self, *obj):
  
    nmrProject = self.nmrProject
    shiftLists0 = getShiftLists(nmrProject)
    shiftLists = [None,] + shiftLists0
    shiftListNames = ['<Ignore>',] + self.getShiftListNames(shiftLists0)
    
    if self.shiftList not in shiftLists:
      self.shiftList = shiftLists[0]
    
    index = shiftLists.index(self.shiftList)
       
    self.shiftListPulldown.setup(shiftListNames, shiftLists, index)

  def changeShiftList(self, shiftList):
    
    self.shiftList = shiftList
    self.updateAfter()

  def changeNucleotideOpt(self, i):
  
     self.nucleotideOpt = i
     self.updateAfter()

  def changeResCode(self, ccpCode):

    self.ccpCode = ccpCode
    self.updateAfter()

  def setChain(self, chain):
  
    if chain is not self.chain:
      self.chain = chain
      self.updateChains()
      self.updateChemElements()
      self.updateResidues()
      self.updateLabellingSchemes()
      self.updateAfter()
    
  def changeChain(self, chain):

    self.chain = chain
    self.updateChemElements()
    self.updateResidues()
    self.updateAfter()

  
  def toggleChemElement(self, name):

    self.chemElementDict = {}
    self.chemElements = self.elementSelect.getSelected()

    for chemElement in self.chemElements:
      self.chemElementDict[chemElement] = 1
    
    self.updateAfter()
  
  def toggleJbonds(self, name):

    if self.refExperiment:
      self.updateAfter()

  def updateAfter(self, *opt):
  
    if self.refresh:
      return
    else:
      self.refresh = True
      self.after_idle(self.update)

  def getChainChemElements(self, chain):
  
    if hasattr(chain, 'chemElements'):
      chemElements = chain.chemElements
    else:
      self.updateChainAttrs(self.chain)
      chemElements = chain.chemElements
    
    chemElements.sort()
      
    return chemElements

  def getChainName(self, chain):
   
    return '%s:%s(%s)' % (chain.molSystem.name,chain.code,chain.molecule.molType)

  def getSystemChains(self):
    
    chains = []
    for molSystem in self.project.sortedMolSystems():
      for chain in molSystem.sortedChains():
        chains.append(chain)
    
    return chains

  def getChainCcpCodes(self, chain):

    if not chain:
      return []
    
    codeDict  = {}
    for residue in chain.residues:
      codeDict[getResidueCode(residue)] = True
    
    codes = codeDict.keys()
    codes.sort()
  
    return codes
  
  def updateResidues(self):

    codes = self.getChainCcpCodes(self.chain)
    self.resPulldown.setup(['<All>'] + codes, [None,] + codes, 0)
    self.ccpCode = None
 
  def updateChains(self, *opt):
  
    chains = self.getSystemChains()
    
    names = []
    index = 0
    chain = self.chain
    
    if chains:
      if chain not in chains:
        chain = chains[0]
        
      names = [self.getChainName(x) for x in chains]
      index = chains.index(chain)
    else:
      chain = None
  
    self.chainPulldown.setup(names, chains, index)
 
    if chain is not self.chain:
      self.chain = chain
      if chain:
        self.updateChainAttrs(chain)
      self.updateChemElements()
      self.updateResidues()
      
    self.updateAfter()

  def updateChemElements(self):
      
    self.chemElements = []
    colors = []
    if self.chain:
      
      chemElements = self.getChainChemElements(self.chain)
        
      for chemElement in chemElements:
        colors.append(ATOM_COLOR_DICT.get(chemElement) or '#a0d0d0')
        if self.chemElementDict.get(chemElement):
          self.chemElements.append(chemElement)
       
      self.elementSelect.update(objects=chemElements,labels=chemElements,colors=colors,selected=self.chemElements)

  def updateButtons(self):
  
    button = self.bottomButtons.buttons[0]

    if self.residue and ( self.residue.ccpCode in FLIP_AROMATICS ):
      # TBD other funnies
      atom = self.residue.findFirstAtom(name='CD1')
      if (not atom) or (not atom.atomSet):
        return
      
      ccpCode = getResidueCode(self.residue)
      if len(atom.atomSet.atoms) > 1:
        text = 'Make %d %s\nring non-equivalent' % (self.residue.seqCode,ccpCode)
      else:
        text = 'Make %d %s\nring equivalent' % (self.residue.seqCode,ccpCode)
      button.enable()
      
    else:
      text = 'Set Ring Flip\nEquivalency'
      button.disable()

    button.config(text=text)
  
    buttons = self.bottomButtons.buttons
    if self.atomSetMapping and self.atomSetMapping.resonanceSerials:
      buttons[1].enable()
      buttons[2].enable()
      buttons[3].enable()
    else:
      buttons[1].disable()
      buttons[2].disable()
      
      if self.residue:
        buttons[3].enable()
      else:
        buttons[3].disable()
  
  def getUnicodeAtomName(self, mapping, molType):
  
    name = mapping.name
    
    if molType == 'protein':
      # Back compatibility
      if name == 'Hn':
        name = 'H'
      else:
        name = mapping.name
    
      return getUnicodeAtomName(name, mapping.elementSymbol)
        
    else:
      return name
      

  def update(self):
    
    self.residue = None
    self.atomSetMapping = None
    self.updateButtons()
    
    chemElements = self.chemElements
    minFraction  = self.minFractionEntry.get() or 0.0
    headingList  = ['#','Residue']
    objectList   = []
    textMatrix   = []
    colorMatrix  = []
    objectListAppend = objectList.append
    textMatrixAppend = textMatrix.append
    colorMatrixAppend = colorMatrix.append

    
    getUnicodeName  = self.getUnicodeAtomName
    jCouplingBonds  = self.jBondsSelect.getSelected()
    residueASMLists = []
    resonanceDict   = {}
    
    nmrProject = self.nmrProject
    nucleotideOpt = self.nucleotideOpt
    refExperiment = self.refExperiment
    labelling = self.labellingScheme
    shiftList = self.shiftList
    
    maxAtoms = 10
    if self.chain:
      self.updateChainAttrs(self.chain)
      if shiftList:
        for resonance in nmrProject.resonances:
          resonanceDict[resonance.serial] = resonance
        
      for residue in self.chain.sortedResidues():
        
        molType = residue.molResidue.molType
        tlc = getResidueCode(residue)
        if self.ccpCode:
          # To be resurrected at some stage (with tabs)
          # Should really be derived or modelled attrs
          """
          if PROTEIN_RESIDUE_CLASS_DICT.get(self.ccpCode):
            if residue.ccpCode not in PROTEIN_RESIDUE_CLASS_DICT[self.ccpCode]:
              continue
          elif tlc != self.ccpCode:
            continue
          """
          if tlc != self.ccpCode:
            continue
 
        obsAtomSets = {}
        if refExperiment or labelling:
          obsAtoms = getResidueObservableAtoms(residue, refExperiment,
                                               labelling,
                                               minFraction, jCouplingBonds,
                                               usePermissiveShifts=False,
                                               chemElements=chemElements)
          
          for obsAtom in obsAtoms:
            atomSet = obsAtom.atomSet
          
            if atomSet:
              obsAtomSets[atomSet.serial] = True
        
        residueMapping  = getResidueMapping(residue)
        atomSetMappings = []
        for atomSetMapping in residueMapping.atomSetMappings:
          if refExperiment or labelling:
            for atomSetSerial in atomSetMapping.atomSetSerials:
              if not obsAtomSets.get(atomSetSerial):
                break
            
            else:
              if MAP_TYPE_OPTS[atomSetMapping.mappingType]:
                unicodeName = getUnicodeName(atomSetMapping, molType)
                atomSetMappings.append( (unicodeName,atomSetMapping) )
          
          else:
            if atomSetMapping.elementSymbol in chemElements:
              if MAP_TYPE_OPTS[atomSetMapping.mappingType]:
                unicodeName = getUnicodeName(atomSetMapping, molType)
                atomSetMappings.append( (unicodeName,atomSetMapping) )
        
        objectListAppend( [residue,residue] )
        
        if (molType in ('RNA','DNA')) and (nucleotideOpt < 3):
          isSugarP = {}
          for name, atomSetMapping in atomSetMappings:
            if ("'" in name) or ("p" in name) or ("P" in name):
              isSugarP[atomSetMapping] = True
                       
          if nucleotideOpt == 0:
            objectListAppend( [residue,residue] )
            textMatrixAppend( ['%d%s' % (residue.seqCode,'b'),tlc])
            textMatrixAppend( ['%d%s' % (residue.seqCode,'s'),tlc])
            colorMatrixAppend(['#b0b0b0','#b0b0b0'])
            colorMatrixAppend(['#b090b0','#b090b0'])
            
            atomSetMappings2 = [nameAsm for nameAsm in atomSetMappings if not isSugarP.get(nameAsm[1])]
            maxAtoms = max(maxAtoms, len(atomSetMappings2)+1)
            atomSetMappings2 = greekSortAtomNames(atomSetMappings2, molType=molType)
            residueASMLists.append( atomSetMappings2 )
 
            atomSetMappings = [nameAsm for nameAsm in atomSetMappings if isSugarP.get(nameAsm[1])]
            
          elif nucleotideOpt == 1:
            textMatrixAppend( ['%d%s' % (residue.seqCode,'b'),tlc])
            colorMatrixAppend(['#b0b0b0','#b0b0b0'])
            atomSetMappings = [nameAsm for nameAsm in atomSetMappings if not isSugarP.get(nameAsm[1])]
            
          else:
            textMatrixAppend( ['%d%s' % (residue.seqCode,'s'),tlc])
            colorMatrixAppend(['#b090b0','#b090b0'])
            atomSetMappings = [nameAsm for nameAsm in atomSetMappings if isSugarP.get(nameAsm[1])]                      
          
        else:
          textMatrixAppend( [residue.seqCode,tlc])
          colorMatrixAppend(['#b0b0b0','#b0b0b0'])

        maxAtoms = max(maxAtoms, len(atomSetMappings)+1)
        atomSetMappings = greekSortAtomNames(atomSetMappings, molType=molType)
        residueASMLists.append( atomSetMappings )
  

    assignShow = self.assignSelect.getIndex()
    for i in range(len(residueASMLists)):
      j = 0
      
      tentAtomNames = set()
      if assignShow == 3:
        residue = objectList[i][0]
        for residueProb in residue.residueProbs:
          if not residueProb.weight:
            continue
          
          spinSystem = residueProb.resonanceGroup
          for resonance in spinSystem.resonances:
            tentAtomNames.update(resonance.assignNames)
      
      for (name,atomSetMapping) in residueASMLists[i]:
        color = ATOM_COLOR_DICT.get(atomSetMapping.elementSymbol) or '#a0d0d0'
        
        resonanceSerials = atomSetMapping.resonanceSerials
        if resonanceSerials and (assignShow < 2):
          if shiftList:
            for serial in resonanceSerials:
              resonance = resonanceDict.get(serial)
              if resonance and resonance.findFirstShift(parentList=shiftList):
                color = scaleColor(self,color,0.7)
                break
          
          else:
            color = scaleColor(self,color,0.7)
          
          textMatrix[i].append( '%-4s' % name )
          objectList[i].append( atomSetMapping )
          colorMatrix[i].append( color )
          j +=1
        
        elif assignShow == 0:
          textMatrix[i].append( '%-4s' % name )
          objectList[i].append( atomSetMapping )
          colorMatrix[i].append( color )
          j +=1
        
        elif (assignShow == 2) and not resonanceSerials:
          textMatrix[i].append( '%-4s' % name )
          objectList[i].append( atomSetMapping )
          colorMatrix[i].append( color )
          j +=1
        
        elif assignShow == 3: # Tentative
          atomSets = atomSetMapping.atomSets

          for atomSet in atomSets:
            if atomSet.name not in tentAtomNames:
              break
          
          else:  
            textMatrix[i].append( '%-4s' % name )
            objectList[i].append( atomSetMapping )
            colorMatrix[i].append( color )
            j +=1    
           
      for k in range(j, maxAtoms):
        textMatrix[i].append('    ')
        colorMatrix[i].append(None)
        objectList[i].append(None)

    for k in range(0, maxAtoms):
      headingList.append(' ')
  
    if not textMatrix:
      textMatrix = [[]]

    tipTexts = self.scrolledMatrix.tipTexts
    self.residueASMLists = residueASMLists
    self.scrolledMatrix.update(headingList=headingList,
                               objectList=objectList,
                               textMatrix=textMatrix,
                               colorMatrix=colorMatrix,
                               tipTexts=tipTexts)
   
    self.refresh = False
 
  def updateChainAttrs(self, chain):
         
    if hasattr(chain, 'chemElements'):
      return
    else:
      chain.chemElements = []
        
    for residue in chain.sortedResidues():
      residueMapping  = getResidueMapping(residue)
      atomSetMappings = residueMapping.atomSetMappings
  	  
      for atomSetMapping in atomSetMappings:
        if not atomSetMapping.elementSymbol:
          atom = atomSetMapping.findFirstAtomSet().findFirstAtom()
          atomSetMapping.elementSymbol = atom.chemAtom.elementSymbol

        if atomSetMapping.mappingType == 'ambiguous':
          if not hasattr(atomSetMapping, 'subSets'):
            atomSetMapping.subSets = []
            for otherMapping in atomSetMappings:
	      if len(otherMapping.atomSetSerials) == 1:
                if list(otherMapping.atomSetSerials)[0] in atomSetMapping.atomSetSerials:
                  atomSetMapping.subSets.append(otherMapping)
		  
      chemAtoms = residue.molResidue.chemCompVar.chemAtoms
      for chemAtom in chemAtoms:
        if chemAtom.className == 'ChemAtom':
          element = chemAtom.elementSymbol
          if element not in chain.chemElements:
            chain.chemElements.append(element)

  
  def selectCell(self, objects, row, col):
  
    if col < len(objects):
      object = objects[col]
    else:
      object = None
    
    if object:
      if object.className == 'Residue':
        self.residue = object
        self.atomSetMapping = None
        spinSystemPopup = self.parent.popups.get('edit_spin_system')
        if spinSystemPopup:
          molResidue = object.molResidue
          spinSystemPopup.chooseResidue( object )
          spinSystemPopup.chooseTentativeResidue( object )
          spinSystemPopup.chooseType( molResidue.ccpCode, molResidue.molType )

        resonanceBrowser = self.parent.popups.get('browse_resonances')
        if resonanceBrowser:
          resonanceBrowser.assignResidue( object )

        resonanceSelection = self.parent.popups.get('selected_resonances')
        if resonanceSelection:
          resonanceSelection.assignResidue( object )
        
        assignmentPanel = self.parent.popups.get('edit_assignment')
        if assignmentPanel:
          assignmentPanel.chooseResidue( object )
 
      elif object.className == 'AtomSetMapping':
        self.residue = None
        self.atomSetMapping = object
        assignmentPanel = self.parent.popups.get('edit_assignment')
        if assignmentPanel:
          assignmentPanel.chooseAtoms( object )

        resonanceBrowser = self.parent.popups.get('browse_resonances')
        if resonanceBrowser:
          resonanceBrowser.chooseAtoms( object )

        resonanceSelection = self.parent.popups.get('selected_resonances')
        if resonanceSelection:
          resonanceSelection.chooseAtoms( object )

    else:
      self.residue = None   
      self.atomSetMapping = None

    self.updateButtons()
 
  def setRequestor(self, requestor=None):
 
    spinSystemPopup = self.parent.popups.get('edit_spin_system')
    if spinSystemPopup and (spinSystemPopup is not requestor):
      spinSystemPopup.cancelAllWaits()

    resonanceBrowser = self.parent.popups.get('browse_resonances')
    if resonanceBrowser and (resonanceBrowser is not requestor):
      resonanceBrowser.cancelAllWaits()
    
    resonanceSelection = self.parent.popups.get('selected_resonances')
    if resonanceSelection:
      resonanceSelection.cancelAllWaits()
 
    assignmentPanel = self.parent.popups.get('edit_assignment')
    if assignmentPanel and (assignmentPanel is not requestor):
      assignmentPanel.cancelAllWaits()
  
  def close(self):
  
    self.setRequestor(None)
    BasePopup.close(self)
 
  def destroy(self):
  
    self.curateNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)
 

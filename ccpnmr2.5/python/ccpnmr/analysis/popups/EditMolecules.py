
"""
======================COPYRIGHT/LICENSE START==========================

EditMolecules.py: Part of the CcpNmr Analysis program

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

import string, re

from ccp.general.ChemCompOverview import chemCompStdDict
from ccp.general.Constants import ccpCodeToCode1LetterDict 
from ccp.general.Constants import code1LetterToCcpCodeDict

from ccp.general.Io import getChemComp, getChemCompCoord

from ccp.util.Molecule import addMolResidues, setMolResidueChemCompVar, setMolResidueCcpCode
from ccp.util.Molecule import makeChainCopy, makeChain, nextChainCode

from ccp.gui.ViewChemCompVarFrame import ViewChemCompVarFrame

from ccpnmr.analysis.core.AssignmentBasic import isChainAssigned
from ccpnmr.analysis.core.MoleculeBasic import getChemCompOverview, getMolTypeCcpCodes, copyMolecule, getResidueCode
from ccpnmr.analysis.popups.BasePopup import BasePopup

from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.DataEntry import askString
from memops.gui.Entry import Entry
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.LabelDivider import LabelDivider
from memops.gui.Label import Label
from memops.gui.MessageReporter import showWarning, showOkCancel, showYesNo
from memops.gui.PulldownList import PulldownList
from memops.gui.RadioButtons import RadioButtons
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.ScrolledText import ScrolledText
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.ToggleLabel import ToggleLabel
from memops.gui.Util import getPopup

from memops.universal.Io import joinPath

WHITESPACE = '\t\n\x0b\x0c\r '

RES_CODE_MODES = ['1-Letter','3-Letter/Ccp']

POLYMER_MOLTYPES = ['protein','DNA','RNA','carbohydrate','DNA/RNA']

ALPHABET = set(string.uppercase)

# TBD
#
# Tab swap "table edits"
# Small mol notional import funcs
# Util code to basic


# Simplification
#
# No Mol systems tab.
# No template molecules tab.
# No links tab

class EditMoleculesPopup(BasePopup):
  """
  **Define Molecular Information in a CCPN Project**
  
  This popup window is used to setup all of the molecular information that will
  be used in a CCPN project. Naturally, the objective is to define the
  molecular entities that were present in the sample when an NMR experiment
  was run so that the resonances that appear in the resulting spectra can
  be ascribed to particular atoms, residues and polymer chains.
  
  **Chemical Elements & Compounds**
  
  In general this system is concerned with several different concepts which
  combine to give the final molecular description. The first concept is that of
  the reference information that relates to chemical elements and chemical
  compounds. As far as the user is concerned this kind of information is mostly
  fixed; there should be only one invariant definition of of the element carbon
  or the amino acid histidine for example. A chemical compound in this respect
  contains a description of a group of atoms and how they are bound together, and
  also lists the different forms that the compound will take, like protonation
  states. For a given type of molecule (protein, DNA, RNA etc) a chemical
  compound will be identified by a residue code name, e.g. histidine is "His". In
  CCPN this name is usually referred to as the "Ccp code", and is often, but not
  always three letters. The common biological residues will also have a
  one-letter code, e.g. histidine is "H", that may be used in some
  circumstances. Defining the sequence of a biopolymer with residue code names,
  or letters, is really selecting a list of connected chemical compound
  descriptions. Most of the chemical compounds that are in common use will be
  pre-defined and available to any CCPN project; either directly or via
  automated download. It is only bespoke or unusual compounds that will need to
  be defined by the user and imported. For example the user can import PDB and
  mol2 format compound descriptions and import them into CCPN via the
  Format Converter.

  **Molecule Templates**

  The second major molecular concept in CCPN is the molecule template, which
  records a sequence of chemical compounds and how they are linked together. A
  template will also record which version of a chemical compound is used, for
  example whether a Cys residue is disulphide linked and how acidic and basic
  residues are protonated. A molecule template does not directly contain any
  atom information, but this is available indirectly by virtue of the template
  residues linking to a given chemical compound. As the name suggests, a
  molecule template is not the entity that is directly used during the NMR
  assignment process; it is just a sequence template that is used to make one or
  more molecular systems. If a single chemical compound is going to be used  in
  NMR assignment, e.g. a ligand or cofactor, then that compound will still be
  used within a molecule template, it is just that the template will only
  contain one sequence item/residue. 

  **Molecular Systems**
  
  Molecular systems are groups of atoms, residues and chains (protein, DNA
  etc..) that are involved in the assignment process. They are constructed
  using sequence information from molecule templates and atom information from
  the chemical compound descriptions. Molecular systems can be thought of as
  representing all of the atoms that are present within an NMR sample or
  biolmolecular assembly. If a molecular system is a homo-oligomeric protein
  complex, where there are several chains of the same sequence, then the sample
  molecule template will be used to generate each of the sub-units, but the
  result will be three separate polypeptide chains, and each can be assigned
  separately to NMR data. For example the bacteriorhodposin complex is a
  homo-trimer, and so to make this complex available in CCPN the user will enter
  the sequence of the protein chain once (to make the a template) and then use
  the sequence three times to make three chains A, B and C.

  **Chains Tab**
  
  The first tab in the popup window lists all on the assignable molecular chains
  that are available in the project. For example if the CCPN project was setup
  for a dimer then you might have chain A and chain B listed in a given
  molecular system. Here the user can delete chains, if they are not assigned,
  and create new chains. There are three basic ways of making chains, the first
  is to copy an existing chain, although naturally the user is asked to select a
  different identifying chain letter/code to distinguish it from the original.
  The second way is to make a chain from an existing molecule template ([Make
  Chain From Template]) and the third way is to make a new chain and a new
  molecular template at the same time (via [Add Sequence]).

  Most commonly, for simple monomeric protein systems, the user will click [Add
  Sequence]. After entering a sequence of residue codes, and making sure that a
  new molecular system and a new template molecule are to be made, a new
  molecular chain is created and ready for assignment. If the molecular system
  is more complicated and involves several different chains, then the user may
  use [Add Chain] several times for each of the chains within the complex,
  making sure that the same destination molecular system is chosen each time.
  Alternatively, the user can add all of the template molecules first, and then
  make the chains from the template molecules afterwards. This method gives most
  flexibility and is recommended if any of the molecules have a non-standard
  connectivity or pronation state; this cannot be changed in a template once a
  chain has bee created. For example to add insulin, which has two discrete
  regions of polypeptide that are disulphide linked, the user would make a
  molecule first by adding two polypeptide sections. change the Cys residues to
  the disulphide state and make the two disulphide links. Only after the molecule
  template definition was complete would the user [Make Chain From Template].

  The lower table lists the "Chain Fragments" that are present in the  chain
  entity selected in the upper table. Each fragment represents a discrete
  segment of biopolymer sequence, and although fragments are bound together in
  the same chain they are not connected by the usual biopolymer links. For
  example insulin will consist of two polypeptide chain fragments that are
  connected together by two disulphide bonds. Most of the time however a
  biopolymer will consist of only one chain fragment with a continuous backbone.
  The starting sequence number of a chain fragment can be edited by the user, as
  long as fragments do not overlap numerically. In this way the user can control
  the numbering that is present in a chain and thus used in assignment. This
  sequence numbering is reversible and can be changed at any time with immediate
  effect (although Analysis can be slow to update if there are many
  assignments).

  **Mol Systems Tab**
  
  This table list all of the molecular systems within the CCPN projects. These
  are used to contain the chains that go together to form a multi-component
  complex or NMR sample. The user can add new, empty molecular systems with a
  view to adding assignable chains to them at a later time. The user can also
  delete a molecular system, and thus all of the chains that it contains, but
  only if these chains do not have any assignments.

  **Sequences Tab**
  
  This table is used to list the sequences that are defined within the CCPN
  project. The user selects the molecular template in the upper pulldown menu
  and the sequence is listed in the table. If a sequence has not been used to
  make any chains (in a molecular system) then the sequence can be changed; the
  user can adjust the kind of residue at a given position, the compound version
  or (descriptor) and how the residue links to others. For example the user
  could setup protonation states, disulphide links and join/break sequence
  fragments. For a new molecule the user can [Add Polymer] and [Add Compound]
  multiple times, thus creating discrete bits of sequence that can be combined
  together into a larger molecule. In this way the user could add a polymer
  sequence and then a compound representing a cofactor. Alternatively, to
  construct insulin [Add Polymer] twice and then [Edit Links] after setting Cys
  residue "descriptors" to "link:SG".
  
  The [Copy Molecule] function is used to make a new molecule based on the
  sequence of an existing one, but this is intended to help add similar
  sequences with small differences in residues, protonation states and links.
  It is not intended to make identical copies of a molecule; identical copies
  are never needed since a single chain can be used multiple times in as many
  chains and molecular systems as is required. 

  The [Lock Molecule] function can be used to prevent any further alteration to
  a sequence. Any molecule that has been used to make any chains (for
  assignment) will be automatically be locked and cannot be unlocked until its
  chains have been deleted. This locking mechanism means that chain any
  assignments then may carry cannot be ruined by editing sequences. If the user
  genuinely made a mistake in entering a sequence after assignments have already
  been made to a chain then the best course of action is to make a new molecule
  template, with the correct sequence, a new molecular system, i.e. building
  chains using the new template, and copy the assignments from the old molecular
  system to the new one via `Copy Assignments`_; selecting "Between Molecule
  Chains".

  **Template Molecules Tab** 
  
  This section lists all of the available sequence templates in the upper table
  and a list of the selected molecule's attributes in the lower table. The user
  can copy, add, extend and delete molecules. Any fine scale editing of the
  molecule template should be done via the adjacent "Sequences" tab. In keeping
  with the above description [Delete Molecule], [Add Polymer Residues] and [Add
  Compound] may only be used if a molecule is unlocked and has not been used to
  make any chain entities, i.e. in a molecular system. Otherwise the [Add ...]
  functions operate in the same way as the equivalent in the "Sequences" tab.

  **Links**
  
  This tab is used to specify how the sequence elements of a molecule template
  are connected together. It should be noted that these links can only be edited
  if the template in question is not used in any chains and is unlocked (via the
  Sequences tab). The general idea is that in the sequences table you select a
  template residue and click [Edit Links] or you select the required template
  molecule and residue in the pulldown menus in the "Intra Molecule Links"
  section.

  With the required template residue selected the upper table is filled with
  rows that represent the covalent links that are, or could be, made from the
  selected residue to other "Destination" template residues. If the molecule is
  unlocked then the user can change or set the residue that is being linked to
  by double clicking in the "Destination Residue" column and selecting the 
  required residue option. Note that only template residues that have an
  unsatisfied link will be shown as possible targets. If the desired destination
  is not free you can go to that residue and remove an existing link (setting
  destination to "<None>") before making the required a new one. Most biopolymer
  residues will show the "prev" and "next" connections that represent the
  regular backbone connectivity; peptide for proteins, phosphodiester for
  nucleic acids. These connections are rarely changed, but can be adjusted to
  beak or join a biopolymer. More typically potential side chain and branching
  links are adjusted. For example if a Cys amino acid has had its descriptor
  set to "Link:SG" then an extra type of link will appear in the table. The user
  can then connect the SG of the Cys to any other free Cys SG, thus specifying
  a disulphide bridge.

  The "Inter Molecule Links" section is not currently used, but will be filled
  in due course.

  **Add Sequence**  

  This section is used to add new biopolymer sequences to a CCPN project. The
  sequence can be used to make a new molecule template, extend an existing
  molecule template or make a molecule template and an assignable chain at the
  same time. Initially the user should set the required polymer type to  say
  what kind of molecule is being made; "GCAT" in DNA is different to "GCAT" in
  protein. Then the user sets whether the sequence will be entered as one-letter
  codes or as Ccp codes (often three-letter residue codes). Note that
  three-letter codes could be interpreted as one-letter codes without error,
  i.e. "ALA" could be interpreted as "Ala, Leu, Ala" if the one-letter option
  was on. Once the sequence is entered correctly the user can switch between the
  one-letter and Ccp codes without a problem. The start number sets what the
  number of the first residue in the sequence will be, although this can be
  changed at any time after. The "Cyclic" option means that a biopolymer link
  will be made between the first and last residue, e.g. for a protein the amide
  which would otherwise be the N-terminus is peptide linked to the last carbonyl
  group. 

  The "Destination Molecule" and "Destination Mol System" settings are important
  because they govern how a sequence will be used in relation to existing 
  molecular descriptions. If the destination molecule is set to anything
  other than "<New>" then the system will attempt to add the entered sequence on
  to the end of the selected molecule template. This may fail if the selected
  template has chains or is otherwise locked. If a new section of sequence is
  added to an existing molecule template then the user can adjust how those
  template residues are connected to any others via the "Sequences" and "Links"
  tabs. The setting for the destination molecular system dictates whether the
  molecule template that the sequence will make (or will extend) is used
  immediately to make an assignable chain entity. For adding simple, single
  protein chains this is usually the case and the "Destination Mol System" is
  set to "<New>", although by selecting an existing molecular system the user can
  add a new subunit to the specification of a complex. Setting the molecular
  system to "<None>" allows the user to add a sequence template without making
  any chains at all. This is important if the molecule template needs further
  adjustment, for example to define non-backbone links or set protonation
  states. If the setting was not "<None>" and a chain was erroneously made, then
  the user should delete the chain in the "Chains" tab and [Unlock Molecule] in
  the "Sequences" tab before making alterations.

  The required biopolymer sequence is entered in the "Sequence Input" text
  panel. Usually this is done by using [Read File] to get the sequence data from
  disk, but it is also possible to place sequences into the text area on Linux
  systems via the mouse middle button, select-and-paste mechanism. Once sequence
  appears in the text area it may be further edited by typing. At any stage the
  user may click the [Tidy] function which will number the sequence and arrange
  it in neat rows. This also has the helpful effect of checking whether the
  sequence is valid and will be interpreted by the CCPN system without error.
  Once the sequence is correct and the other settings, including the destination
  molecular system, have been checked the clicking the [Add Sequence!] button 
  commits the sequence by making the specified molecule template and any new
  chains.

  **Small Compounds**

  Here the user can add small compounds to molecule templates. These may be any
  ligands, co-factors, solvents, metabolites etc. that are assignable. Once a
  compound had been added to a new or existing molecule template, the compound
  can be made accessible for assignment by using its template to create a chain.

  Generally the user selects which kind of biopolymer the compound relates to,
  or "other" of is not associated with any polymer type. In the table on the
  right the user scrolls and selects the required chemical compound row. If
  idealised atomic coordinates are known for the selected compound these are
  displayed in the 3D coordinate display on the left. An absence of such
  coordinates does not matter with regard to selecting a compound for
  assignment. If a compound row is coloured blue, then the compound is available
  locally and may be used immediately. If the row is grey, then the compound will
  be downloaded from a remote database before use. Any downloading of this kind
  is automatic, but naturally requires an active Internet connection. 

  Once the desired small molecule is selected the user next selects the required
  destination molecule template in the pulldown menu below the table. The
  compound can either be put into a new template, i.e. on its own, or added
  to an existing molecule, as would normally be done for a covalently linked
  cofactor. Finally [Add Compound] actually added the selected
  compound to the molecule template. To actually assign to this molecule
  the molecule template must be used to make a chain (i.e. in a given molecular
  system), which is readily achieved by using the [Make Chain From Template]
  function from the "Chains" tab.

  **Caveats & Tips**

  In order to delete a chain that carries assignments, the chain must be cleared
  of assignments first. The assignments could be moved to a different chain via
  the `Copy Assignments`_ system or removed completely. Completely removing a
  chain's assignments can be done via the `Spin Systems`_ or Resonances_ tables
  by selecting the rows that relate to the relevant chain (sorting rows by
  residue may help) and then deassigning the spin systems or resonances.

  The sequence numbering of a chain can be altered at any stage by editing the
  "Start Seq Number" in the "Chain Fragments" table of the "Chains" tab.
  Although the start number can be adjusted, currently the residue numbering
  will always be sequential within a given fragment. The residue template
  numbering in the "Sequences" tab cannot be changed, although this may change
  in the future. However, the template numbering system is separate from the
  chain numbering system.

  .. _`Copy Assignments`: CopyAssignmentsPopup.html
  .. _`Spin Systems`: EditSpinSystemPopup.html
  .. _Resonances: BrowseResonancesPopup.html
  
  """ 

  def __init__(self, parent, *args, **kw):

    self.guiParent  = parent
    self.waiting = False
    self.waitingMolRes = False
    self.waitingChain = False
    
    self.molSystem = None
    self.chain = None
    self.molecule = False # Not None
    self.moleculeB = None
    self.moleculeC = None
    self.molSystemC = None
    self.molResidue = None
    self.fragment = None

    self.molResLinkEnd     = None
    self.linkMolecule      = None
    self.waitingMolResLink = False
    self.linkMolResidue    = None
    self.linkCcpCode       = None

    self.molTypes     = ('other',)
    self.compoundKey  = None
    self.chemCompVar  = None
     
    self.inputSeq    = ''
    self.seqCodeMode = RES_CODE_MODES[0]
    self.polymerMolType = POLYMER_MOLTYPES[0]
    self.seqStartNum = 1
    
    self.moleculeS = None
    self.molSystemS = True
    self.moleculeComp = None

    # Safety code - could put elsewhere
                    
    BasePopup.__init__(self, parent=parent, title="Molecule : Molecules", **kw)
  
    for molecule in self.project.molecules:
      if molecule.chains:
        molecule.isFinalised = True    

  def body(self, guiFrame):

    self.geometry('720x500')
      
    guiFrame.grid_columnconfigure(0,weight=1)
    guiFrame.grid_rowconfigure(0,weight=1)


    tipTexts = ['A table of all the molecular chains used in assignment, possibly from various samples and complexes',
                'A table of the molecular system groups: collections of polymer chains and small molecules that are in the same sample or complex',
                'A table of the residue sequences for the molecules defined in the project, including bio-polymer links and protonation state',
                'A table of all of the molecules in the project; used as sequence templates to build assignable "chain" entities',
                'A table to view and setup links and bonds between residues that are not part of the normal bio-polymer connectivity',
                'Load or paste bio-polymer (protein, DNA, RNA, carbohydrate) sequences to build molecule templates',
                'A system to choose smaller chemical components for an assignable system; includes ligands, cofactors and unusual residues']
    options = ['Chains','Mol Systems','Sequences',
               'Template Molecules','Links',
               'Add Sequence','Small Compounds']
    tabbedFrame = TabbedFrame(guiFrame, options=options, tipTexts=tipTexts,
                              grid=(0,0), callback=self.updateTab)
    self.tabbedFrame = tabbedFrame
    
    frameA, frameB, frameC, frameD, frameE, frameF, frameG = tabbedFrame.frames
    
    #
    # Chains Frame
    #
    
    frameA.grid_columnconfigure(0,weight=1)
    frameA.grid_rowconfigure(0,weight=1)

    tipTexts = ['The short name of the molecular system to which the chain belongs; a grouping of chains, e.g. in a complex',
                'The short textual identifier for the chain, for graphical displays and assignment annotations (where required)',
                'The name of the molecular sequence template on which the chain is based',
                'The number of residue monomers (or other separate chemical compounds) that make up the chain',
                'The number of discontinuous polymer sequence fragments into which the chain is broken',
                'The kind of bio-polymer the chain represents; usually just a single kind, but hybrids like DNA/RNA are possible',
                'A user-editable textual comment for the molecular chain, which defaults to the initial snippet of one-letter sequence']

    headingList = ['Mol\nSystem','Chain\nCode','Molecule\nTemplate','Residues',
                   'Chain\nFragments','Molecular\nTypes','Details']
                   
    self.editChainDetlEntry = Entry(self, width=12, text='',
                                    returnCallback=self.setChainDetl)

    
    editWidgets      = [None, None, None, None, None, None, self.editChainDetlEntry]
    editGetCallbacks = [None, None, None, None, None, None, self.getChainDetl]
    editSetCallbacks = [None, None, None, None, None, None, self.setChainDetl]
    
    self.chainsMatrix = ScrolledMatrix(frameA, grid=(0,0),
                                       tipTexts=tipTexts, 
                                       editSetCallbacks=editSetCallbacks,
                                       editGetCallbacks=editGetCallbacks,
                                       editWidgets=editWidgets,
                                       headingList=headingList,
                                       callback=self.selectChain,
                                       deleteFunc=self.deleteChain)
    
    tipTexts = ['Show a table containing the residue sequence for the chain, together with associated linking and protonation states',
                'Delete the selected chain specification, if it is not assigned to NMR resonances (via its atoms)',
                'Make a duplicate of the selected chain within its molecular system, using the next available chain code letter; convenient for making homo-oligomers',
                'Add a new bio-polymer sequence, to make a molecule template, from which a new chain may be built',
                'Show a table of all the molecule templates available in the project, from which chains may be built']

    texts = ['Show\nSeq','Delete\nChain', 'Copy\nChain', 'Add\nSequence',
             'Edit Molecule\nTemplates']
    commands = [self.showChainSequence, self.deleteChain, self.copyChain,
                self.chainAddSequence, self.editMolecule]
    
    self.chainButtons = ButtonList(frameA, texts=texts, grid=(1,0),
                                   commands=commands, tipTexts=tipTexts)
     
    frame = Frame(frameA)
    frame.grid(row=2, column=0, sticky='ew')
    frame.grid_columnconfigure(4,weight=1)
 
    label = Label(frame, text='Mol System for new chain:', grid=(0,0))

    tipText = 'Selects which molecular system group to place new chains in (when constructed with a sequence template)'
    self.chainMolSystemPulldown = PulldownList(frame, grid=(0,1), tipText=tipText,
                                               callback=self.selectChainMolSystem)

    label = Label(frame, text='Template for new chain:', grid=(0,2))
    
    tipText = 'Selects which molecular sequence template is used to make a new chain (which can then be assigned)'
    self.chainMolPulldown = PulldownList(frame, grid=(0,3), tipText=tipText,
                                         callback=self.setChainMolecule)

    tipTexts = ['Using the selected molecule template make a new chain withing the stated molecular system group; the chain made can be used for NMR assignments',]
    texts = ['Make Chain From Template',]
    commands = [self.makeChain,]
    
    self.chainButtons2 = ButtonList(frame, texts=texts, grid=(0,4),
                                   commands=commands, tipTexts=tipTexts)
    self.chainButtons2.buttons[0].config(bg='#C0E0C0')
    
    # Fragments
    
    

    self.fragmentLabel = LabelDivider(frameA, text='Chain Fragments')
    self.fragmentLabel.grid(row=3, column=0, sticky='ew')
    
    self.editStartNumlEntry = IntEntry(self, text=0, width=6,
                                       returnCallback=self.setStartNum)
                                       
    frameA.grid_rowconfigure(4,weight=1)
    tipTexts = ['The serial number of the fragment within a chain',
                'The type of bio-polymer represented in the fragment',
                'The number of residues connected within the fragment of the chain',
                'The number of the first residue in the fragment, for graphical display; change this to set the residue numbering used in assignment',
                'Whether the fragment represents a standard linear bio-polymer',
                'A short snippet of the fragments sequence, for human appreciation']
    headingList = ['#','Mol Type','Residues','Start Seq\nNumber',
                   'Linear\nPolymer?','Sequence']
                   
    editWidgets      = [None, None, None, self.editStartNumlEntry, None, None]
    editGetCallbacks = [None, None, None, self.getStartNum,        None, None]
    editSetCallbacks = [None, None, None, self.setStartNum,        None, None]
    
    self.fragmentMatrix = ScrolledMatrix(frameA, grid=(4,0), tipTexts=tipTexts,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets,
                                         headingList=headingList,
                                         callback=self.selectFragment)
                                         
    
    #
    # MolSystem Frame
    #
    
    frameB.grid_columnconfigure(0,weight=1)
    frameB.grid_rowconfigure(0,weight=1)

    tipTexts = ['A sort textual code for the molecular system grouping; set when the molecular system is first created',
                'An editable, descriptive name for the molecular system, e.g. the name of a multimeric complex',
                'The identifying codes of the molecular chains that are present in the molecular system group',
                'The names of the molecule templates which have been used to make the chains of the molecular system',
                'The bio-polymer or other molecule types present in the molecular system',
                'The number of structure ensembles (coordinate sets) associated with the molecular system',
                'A user-edited verbose description or comment about the molecular system']
    colHeadings = ['Code', 'Name', 'Chains', 'Molecules','Mol Types',
                   'Structures','Details']
                   
    self.editMolSysNameEntry = Entry(self, text='', width=8,
                                     returnCallback=self.setMolSysName )
    self.editMolSysDetlEntry = Entry(self, text='', width=12,
                                     returnCallback=self.setMolSysDetl)
    
    editWidgets   = [None, self.editMolSysNameEntry, None, None, 
                     None, None, self.editMolSysDetlEntry]
                     
    editGetCallbacks = [None, self.getMolSysName, None, None, 
                        None, None, self.getMolSysDetl]
                        
    editSetCallbacks = [None, self.setMolSysName, None ,None, 
                        None, None, self.setMolSysDetl]
    
    self.molSysMatrix = ScrolledMatrix(frameB, tipTexts=tipTexts,
                                       editSetCallbacks=editSetCallbacks,
                                       editGetCallbacks=editGetCallbacks,
                                       editWidgets=editWidgets,
                                       headingList=colHeadings,
                                       callback=self.selectMolSystem,
                                       deleteFunc=self.delMolSys,
                                       grid=(0,0))
 
    tipTexts = ['Add a new, blank molecular system to the CCPN project; which will then go on to contain chains of residues and their atoms',
                'Delete the selected molecular system definition, but only if it carries no NMR assignments']
    texts = ['Add Mol System','Delete Mol System']
    commands = [self.addMolSys,self.delMolSys]
    self.molSysButtons = ButtonList(frameB, texts=texts, tipTexts=tipTexts,
                                    commands=commands, grid=(1,0))
     
    
    #
    # Mol Residue Frame
    #
    
    frameC.grid_columnconfigure(2,weight=1)
    frameC.grid_rowconfigure(1,weight=1)
    
    label = Label(frameC, text='Molecule: ')
    label.grid(row=0, column=0, sticky='w')

    tipText = 'The name of the molecular template to display the sequence for'
    self.moleculePulldown = PulldownList(frameC, tipText=tipText, grid=(0,1),
                                         callback=self.selectResMolecule)
    
    tipText = 'An indication of which chains the template has been used to construct (inside the molecular systems)'
    self.molResChainInfoLabel = Label(frameC, tipText=tipText, grid=(0,2),
                                      text='Chains with this sequence:')
                                    
    self.descriptorPulldown = PulldownList(self, callback=self.setDescriptor)
    self.linkingPulldown = PulldownList(self, callback=self.setLinking)
    self.ccpCodePulldown = PulldownList(self, callback=self.setCcpCode)
    
    tipTexts = ['The number of the residue in the sequence; separate to the chain (fragment) numbering which can be edited by the user',
                'The identifying code for the kind of residue; usually three letters, but not always so',
                'A description of the residue\'s protonation and stereochemical state; can only be set for unlocked templates, which have not yet been used to make chains',
                'The bio-polymer or other type of the residue',
                'The position of the residue relative to the bio-polymer backbone; a description of which kinds of links are made',
                'The identities of residues that are covalently linked, including both regular bio-polymer links and other non-standard kinds']

    colHeadings      = ['Seq\nCode','Ccp\nCode','Descriptor &\nStereochemistry',
                        'Mol Type','Polymer\nLinking','Linked\nResidues']
                        
    #editWidgets      = [None, self.ccpCodePulldown, self.descriptorPulldown, 
    #                    None, self.linkingPulldown, None]
    #editGetCallbacks = [None, self.getCcpCode, self.getDescriptor,      
    #                    None, self.getLinking, self.showResLinks]
    #editSetCallbacks = [None, self.setCcpCode, self.setDescriptor,      
    #                    None, self.setLinking, None]
    
    # TBD look at CcpCode setting
    
    editWidgets      = [None, None, self.descriptorPulldown, 
                        None, self.linkingPulldown, None]
    editGetCallbacks = [None, None, self.getDescriptor,      
                        None, self.getLinking, self.showResLinks]
    editSetCallbacks = [None, None, self.setDescriptor,      
                        None, self.setLinking, None]

    self.molResidueMatrix = ScrolledMatrix(frameC, grid=(1,0),
                                           gridSpan=(1,3), tipTexts=tipTexts,
                                           editSetCallbacks=editSetCallbacks,
                                           editGetCallbacks=editGetCallbacks,
                                           editWidgets=editWidgets,
                                           headingList=colHeadings,
                                           callback=self.selectMolResidue)
 
    tipTexts = ['Add a bio-polymer sequence to construct a new molecular template, from which chains may be constructed',
                'Copy the currently visible residue sequence into the "Add Sequence" tool; as a basis for making a similar sequence',
                'Add an isolated compound (small molecule) to the current template; which may include cofactors or sugar residues',
                'Display and edit the covalent links for the selected residue; only possible for unlocked molecules (without assignments)',
                'Unlock the molecular sequence template so that its links and protonation states etc make be modified. Unlocking is only possible if the template has not been used to make chains']

    texts    = ['Add Polymer','Copy Sequence','Add Compound',
                'Edit Links','Unlock Molecule']
                
    commands = [self.addResPolymer,self.copyResSequence,self.addResCompound,
                self.showResLinks,self.toggleMoleculeLock]
                
    self.molResidueButtons = ButtonList(frameC,texts=texts, tipTexts=tipTexts,
                                        commands=commands, grid=(2,0),
                                        gridSpan=(1,3))
    
    
    #
    # Molecules Frame
    #
    
    frameD.grid_columnconfigure(0,weight=1)
    frameD.grid_rowconfigure(0,weight=1)
    frameD.grid_rowconfigure(3,weight=1)
    
    self.longNameEntry = Entry(self, text='', width=16,
                               returnCallback=self.setLongName)
    self.molDetailsEntry = Entry(self, text='', width=16,
                                 returnCallback=self.setMolDetails)
    self.seqMolDetailsEntry = Entry(self, text='', width=16,
                                    returnCallback=self.setSeqDetails)

    tipTexts = ['The short name for the molecule template, for use in graphical displays; set when the template is first created',
                'A verbose name describing the molecular sequence; typically the protein or gene name',
                'The types of bio-polymer or other molecule represented by the sequence template',
                'Which chains the sequence template has been used to construct (inside molecular systems)',
                'User-supplied textual comment for the molecule template',
                'User-supplied comments specifically regarding the sequence of the molecule']
  
    colHeadings  = ['Name','Long Name','Mol Type',
                    'Chains','Details','Seq\nDetails']
    #                'Common\nNames','Keywords','Functions',
                    
    editWidgets      = [None, self.longNameEntry, None, None, 
                        self.molDetailsEntry , self.seqMolDetailsEntry]
    editGetCallbacks = [None, self.getLongName,   None, None, 
                        self.getMolDetails, self.getSeqDetails  ]
    editSetCallbacks = [None, self.setLongName,   None, None, 
                        self.setMolDetails, self.setSeqDetails  ]

    self.moleculeMatrix = ScrolledMatrix(frameD, grid=(0,0),
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets,
                                         initialRows=4, tipTexts=tipTexts,
                                         headingList=colHeadings,
                                         callback=self.selectMolecule,
                                         deleteFunc=self.deleteMolecule)

    tipTexts = ['Show a table containing the residue sequence of the selected molecule template',
                'Add a new, blank molecular template, to which bio-polymer sequence or small molecules may be added',
                'Make another molecular template with the same sequence as the selected; allows identical sequences but different protonation states',
                'Delete the selected molecular sequence template; only allowed if it is not used by any molecular system chains',
                'Add a sequence of bio-polymer residues to the selected template; only allowed if the template is unlocked',
                'Add a non-bio-polymer or small compound to the selected template; only allowed if the template is unlocked']
    texts    = ['Show\nSeq', 'Add New\nMolecule','Copy\nMolecule','Delete\nMolecule', 
                'Add Polymer\nResidues','Add Other\nCompound',]
    commands = [self.showMoleculeSequence, self.addMolecule, self.copyMolecule,
                self.deleteMolecule, self.addMolPolymer, self.addMoleculeCompound,]
    
    self.moleculeButtons = ButtonList(frameD,texts=texts, grid=(1,0),
                                      tipTexts=tipTexts, commands=commands)
    
    #
    
    div = LabelDivider(frameD, text='Molecule Data')
    div.grid(row=2, column=0, sticky='ew')
    
    tipTexts = ['Various characteristics and attributes of the selected molecule',
                'The value of the parameter described on each row']
    colHeadings = ['Attribute','Value']

    self.molDataMatrix = ScrolledMatrix(frameD, grid=(3,0),
                                        headingList=colHeadings,
                                        tipTexts=tipTexts)
                                        
    
    #
    # Links
    #

    frameE.grid_columnconfigure(0,weight=1)
    frameE.grid_rowconfigure(2,weight=1)
    frameE.grid_rowconfigure(4,weight=1)
    
    div = LabelDivider(frameE, text='Intra Molecule Links')
    div.grid(row=0, column=0, sticky='ew')
        
    self.nonSeqLinkPulldown = PulldownList(self, callback=self.setMolResLink)

    subFrame = Frame(frameE)
    subFrame.grid(row=1, column=0, sticky='ew')
    subFrame.grid_columnconfigure(5, weight=1)
    
    label = Label(subFrame, text='Molecule:', grid=(0,0))
    tipText = 'Selects which molecule template to display or edit covalent links for'
    self.linkMoleculePulldown=PulldownList(subFrame, tipText=tipText, grid=(0,1),
                                           callback=self.changeLinkMolecule)
    
    label = Label(subFrame, text='Residue Type:', grid=(0,2))
    tipText = 'Selects which kind of residue to show covalent links for; restricts the contents of the "Source Residue" pulldown menu'
    self.linkCcpCodesPulldown = PulldownList(subFrame, tipText=tipText, grid=(0,3),
                                             callback=self.setLinkCcpCode)

    label = Label(subFrame, text='Source Residue:', grid=(0,4))
    tipText = 'Selects which specific residue, within the selected molecule template, to show covalent links for'
    self.linkResiduePulldown = PulldownList(subFrame, tipText=tipText, grid=(0,5),
                                             callback=self.setLinkMolResidue)

    colHeadings  = ['Source\'s\nLink','Destination\n Residue','Destination\'s\nLink']
    
    tipTexts = ['A name for the kind of link relative to the selected "source" residue; the one which all links emanate from',
                'The sequence number and name of the residue to which the selected "source" residue is linked',
                'The name for the kind of link relative to the residue being linked to from the "source"']
    self.molResLinkMatrix = ScrolledMatrix(frameE, tipTexts=tipTexts,
                                           headingList=colHeadings,
                                           callback=self.selectMolResLinkEnd,
                                           grid=(2,0), gridSpan=(1,2))
                                         
    
    #
    
    tipText = '*This section is not currently used* - To be filled in the future'
    div = LabelDivider(frameE, text='Inter Molecule Links',
                       tipText=tipText, grid=(3,0))
    
    #
    # Polymer Seq Frame
    #
    
    frameF.expandGrid(3,0)
    
    frame = Frame(frameF, grid=(0,0), sticky='ew')
    frame.grid_columnconfigure(8,weight=1)

    label = Label(frame, text='Polymer Type:', grid=(0,0))
    
    texts  = ['Protein','DNA','RNA','Carbohydrate','DNA/RNA mix']
    tipText = 'Selects which kind of bio-polymer sequence is being added'
    self.polymerMolTypePulldown = PulldownList(frame, texts=texts, grid=(0,1),
                                               objects=POLYMER_MOLTYPES,
                                               callback=self.setPolymerSeqMolType,
                                               tipText=tipText)
    
    label = Label(frame, text='Start Number:', grid=(0,2))
    tipText = 'Sets the number of the first residue in the sequence; currently not subsequently editable, but can be adjusted when making chains'
    self.polymerStartEntry = IntEntry(frame, text=self.seqStartNum, grid=(0,3),
                                      width=5, tipText=tipText,
                                      returnCallback=self.tidyPolymerSeq)

    label = Label(frame, text='Cyclic:', grid=(0,4))
    tipText = 'Whether the polymer is circular, with the end of the sated sequence connecting to the beginning in the normal bio-polymer manner'
    self.cyclicSelect = CheckButton(frame, grid=(0,5),
                                    tipText=tipText, selected=False)
    
    label = Label(frame, text='Input Type:', grid=(0,6))
    tipTexts = ['The sequence will be specified in one-lette code form like "GCAT" or "TEAANDCAKE"',
                 'The sequence will be specified in three-letter (or full CCPN) code form like "ALA GLY MET SER"']
    self.seqCodeModeButtons = RadioButtons(frame, RES_CODE_MODES,
                                           self.changeCodeMode,
                                           selected_index=0, borderwidth=1,
                                           relief='flat', grid=(0,7),
                                           sticky = 'w', tipTexts=tipTexts)
    
    frame = Frame(frameF, grid=(1,0), sticky='ew')
    frame.grid_columnconfigure(3,weight=1)
   
    label = Label(frame, text='Destination Molecule:', grid=(0,0))
    tipText = 'Selects which molecule template the sequence will be added to; when making a new molecule you will be prompted for the name'
    self.seqMoleculePulldown = PulldownList(frame, grid=(0,1), tipText=tipText,
                                            callback=self.changeSeqMolecule)

    label = Label(frame, text='Destination Mol System:', grid=(0,2))
    tipText = 'If set, specifies that the sequence should be used to immediately make an assignable chain; this disallows further editing (e.g. of protonation state)'
    self.seqMolSystemPulldown = PulldownList(frame, grid=(0,3), tipText=tipText,
                                             callback=self.changeSeqMolSystem)
   
    label = Label(frame, text='Ccp Codes:', grid=(0,4), sticky='e')
    
    ccpCodes = self.getPolymerCcpCodes()
    tipText = 'Selects a residue type from a list of all of the CCPN residue codes available for the selected kind of bio-polymer molecule'
    self.ccpCodeSeqPulldown = PulldownList(frame, texts=ccpCodes, grid=(0,5),
                                        callback=self.setPolymerSeqCcpCode,
                                        sticky='e', tipText=tipText)
    
    tipText = 'Add a residue of the kind selected in the adjacent pulldown menu to the end of the sequence'
    self.appendButton = Button(frame, text='Append', grid=(0,6),
                               command=self.appendPolymerResidue,
                               borderwidth=1, pady=1, sticky='e',
                               tipText=tipText) 
    
    div = LabelDivider(frameF,text='Sequence Input', grid=(2,0))
    
    tipText = 'Cut, paste and edit bio-polymer residue sequences in this text box'
    self.seqTextBox = ScrolledText(frameF, width=60, height=10, grid=(3,0),
                                   text=self.inputSeq, xscroll=False,
                                   tipText=tipText)
  
    tipTexts = ['Add the states sequence to the selected, or new, molecule making a new chain if a molecular system is specified',
                'Tidy the text window containing the sequence; lines things up and checks if the residue codes will be interpreted properly',
                'Read residue sequences from a file; places the result in the above text window',
                'Clear the text window containing the sequence',
                'For DNA and RNA make a sequence that is the reverse complement of the visible sequence']

    texts    = ['Add Sequence!','Tidy','Read File','Clear','Reverse Complement']
    
    commands = [self.commitPolymerSeq,self.tidyPolymerSeq,
                self.readFile, self.clearPolymerSeq, self.reverseComplement]
                
    self.polymerSeqButtons = ButtonList(frameF, texts=texts, tipTexts=tipTexts,
                                       commands=commands, grid=(4,0))
                                       
    self.polymerSeqButtons.buttons[0].config(bg='#B0FFB0')
    self.reverseCompButton = self.polymerSeqButtons.buttons[-1]
    self.reverseCompButton.disable()
    
    #
    # Small Compound Frame
    #
    
    frameG.grid_columnconfigure(0,weight=0,minsize=400)
    frameG.expandGrid(1,2)
    
    self.chemCompView=ViewChemCompVarFrame(frameG, grid=(0,0), gridSpan=(4,1))
    
    
    label = Label(frameG, text='Mol Type: ', grid=(0,1))
    
    texts = ['Protein','DNA','RNA','Nucleic Acid','Carbohydrate','Other','All']
    molTypeSelections = [('protein',),
                         ('DNA',),
                         ('RNA',),
                         ('DNA','RNA'),
                         ('carbohydrate',),
                         ('other',),
                         ('protein','DNA','RNA','carbohydrate','other')
                        ]
            
    tipText = 'Selects which bio-polymer type to show compounds for; defaults to "other", representing things that do not fit into other categories'
    self.molTypePulldown = PulldownList(frameG, texts=texts, tipText=tipText, 
                                        objects=molTypeSelections,
                                        callback=self.selectMolTypes,
                                        index=5, grid=(0,2))

    tipTexts = ['Import a non-standard compound definition from a CCPN ChemComp XML file, e.g. something made by CcpNmr ChemBuild',]
    texts    = ['Import XML File',]
    commands = [self.importXml,]
    buttons  = ButtonList(frameG, texts=texts, commands=commands,
                          grid=(0,3), tipTexts=tipTexts)
    
    tipTexts = ['The CCPN code of the compound',
                'The PDB three-letter code of the compound',
                'The long chemical name of the compound. *Note this is not standardised*']
    headingList = ['Ccp\nCode', '3-Letter\nCode', 'Long Name']
    justifyList = ['center', 'center', 'left']
    
    self.compoundMatrix=ScrolledMatrix(frameG, tipTexts=tipTexts,
                                       justifyList=justifyList,
                                       headingList=headingList,
                                       callback=self.chooseChemComp,
                                       gridSpan=(1,3), grid=(1,1))
                                       
    
    frame = Frame(frameG, grid=(2,1), gridSpan=(1,3), sticky='ew')
    frame.grid_columnconfigure(1,weight=1)
   
    label = Label(frame, text='Destination Molecule:', grid=(0,0))
    tipText = 'Selects which molecul template the selected chemical compound may be added'
    self.compMoleculePulldown = PulldownList(frame, tipText=tipText,
                                             grid=(0,1),
                                             callback=self.changeCompMolecule)
        
    tipTexts = ['Add the selected chemical compound to the above selected, or new, molecule template; when done the template can be used to make an assignable chain',]
    texts    = [ 'Add Compound', ]
    commands = [ self.commitCompound,]
    buttons  = ButtonList(frameG, texts=texts, commands=commands,
                          gridSpan=(1,3), grid=(3,1), tipTexts=tipTexts)
    buttons.buttons[0].config(bg='#B0FFB0')
         
    
    #
    # Main Frame
    #

    self.bottomButtons = UtilityButtonList(tabbedFrame.sideFrame,
                                           helpUrl=self.help_url,
                                           grid=(0,0), sticky = 'e')
    
    self.updateTopLevel()
    self.updateChains()
    self.updateCompoundTable()
    self.compoundMatrix.selectObject(('Atp','other')) # Just a pretty default
    self.updateButtons()
    
    self.administerNotifiers(self.registerNotify)

  def open(self):
  
    self.updateTab(self.tabbedFrame.selected)    
    
    BasePopup.open(self)
  
  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)


  def administerNotifiers(self, notifyFunc):
  
    
    for func in ('setCommonNames','addCommonName','removeCommonName','setDetails',
                 'setFunctions','addFunction','removeFunction','setIsFinalised',
                 'setKeywords','addKeyword','removeKeyword', 'setCalcIsoelectricPoint',
                 'setLongName','setSeqDetails','setSmiles','setSmilesType'):    
      notifyFunc(self.updateMolDataAfter, 'ccp.molecule.Molecule.Molecule', func)
          
    for func in ('__init__', 'delete'):
      notifyFunc(self.updateTopLevelAfter, 'ccp.molecule.Molecule.Molecule', func)
      notifyFunc(self.updateTopLevelAfter, 'ccp.molecule.MolSystem.MolSystem', func)
         
    for func in ('setName','setDetails','setStructureEnsembles'):
      notifyFunc(self.updateMolSystemsAfter,'ccp.molecule.MolSystem.MolSystem', func)
    
    for func in ('__init__', 'delete'):
      notifyFunc(self.updateMolSystemsAfter, 'ccp.molecule.MolSystem.MolSystemLink', func)

    for func in ('setLongName','setDetails', 'setSeqDetails'):
      notifyFunc(self.updateMoleculesAfter,'ccp.molecule.Molecule.Molecule', func)
    
    for func in ('__init__', 'delete','setSeqInsertCode','setSeqCode'):
      notifyFunc(self.updateMolResiduesAfter, 'ccp.molecule.Molecule.MolResidue', func)

    for func in ('__init__', 'delete'):
      notifyFunc(self.updateMolResiduesAfter, 'ccp.molecule.Molecule.MolResLink', func)
      notifyFunc(self.updateMolResLinksAfter, 'ccp.molecule.Molecule.MolResLink', func)
    
    for func in ('__init__', 'delete'):
      notifyFunc(self.updateChainsAfter, 'ccp.molecule.MolSystem.Chain', func)
      notifyFunc(self.updateMolSystemsAfter, 'ccp.molecule.MolSystem.Chain', func)
      
    for func in ('setDetails', ):
      notifyFunc(self.updateChainsAfter, 'ccp.molecule.MolSystem.Chain', func)

    for func in ('setPhysicalState','setPdbOneLetterCode'):
      notifyFunc(self.updateChainsAfter,'ccp.molecule.MolSystem.Chain', func)
    notifyFunc(self.updateResidueAfter,'ccp.molecule.MolSystem.Residue', 'setSeqCode')    
    
     
    #for func in ('__init__', 'delete'):
    #  notifyFunc('ccp.nmr.Nmr.ResonanceSet', clazz, func)

  #
  # Table Selection 
  #
  
  def selectChain(self, chain, row, col):
  
    if chain is not self.chain:
      self.chain = chain
      self.updateFragments()
      self.updateButtons()  
  
  def selectFragment(self, fragment, row, col):
  
    if fragment is not self.fragment:
      self.fragment = fragment
      self.updateButtons()  
      
  #
  
  def selectMolSystem(self, molSystem, row, col):

    if molSystem is not self.molSystem:
      self.molSystem = molSystem
      self.updateButtons()  
  
  #
  
  def selectMolResidue(self, molResidue, row, col):
  
    if molResidue is not self.molResidue:
      self.molResidue = molResidue
      self.updateButtons()
  
  #
  
  def selectMolecule(self, molecule, row, col):
  
    if molecule is not self.moleculeB:
      self.moleculeB = molecule
      self.updateMolData()
      self.updateButtons()
 
  #
    
  def selectMolResLinkEnd(self, molResLinkEnd, row, col):
  
    if molResLinkEnd is not self.molResLinkEnd:
      self.molResLinkEnd = molResLinkEnd
  

  def _chooseChemCompCoord(self, obj, row, col):
    """New routine to try to load ChemCompCoord files from Url
    Not use which are required so check all in the imported source_names list
    """
    from ccp.gui.ViewChemCompVarFrame import SOURCE_NAMES

    self.compoundKey = obj
    (ccpCode, molType) = self.compoundKey

    chemCompVarCoord = []
    for sourceName in SOURCE_NAMES:
      chemCompCoord = self.project.findFirstChemCompCoord(sourceName=sourceName,
                                                     molType=molType,
                                                     ccpCode=ccpCode)

      if not chemCompCoord:
        chemCompVarCoord.append(getChemCompCoord(self.project, sourceName, molType, ccpCode, download=True))

  def chooseChemComp(self, obj, row, col):
  
    self.compoundKey = obj
    (ccpCode,molType) = self.compoundKey
    
    chemComp = self.project.findFirstChemComp(molType=molType, ccpCode=ccpCode)
    
    if not chemComp:
      msg = 'Chemical component not available locally. OK to download to project directory?'
      if showOkCancel('Query', msg, parent=self):
        chemComp = getChemComp(self.project, molType, ccpCode, download=True)

        # 20190321:ED read the chemCompCoord files before view update
        # assumes that existence of chemComp implies existence of Coord files, otherwise will test on every click
        self._chooseChemCompCoord(obj, row, col)

    if chemComp:
      chemCompVar = chemComp.findFirstChemCompVar(linking='middle', descriptor='neutral') \
                    or chemComp.findFirstChemCompVar(linking='none') \
                    or chemComp.findFirstChemCompVar()
                      
      self.chemCompVar = chemCompVar
      self.chemCompView.update(chemCompVar)
    else:
      showWarning('Warning','Chemical compound definition not available.', parent=self)
      self.chemCompVar = None
      self.chemCompView.update(None)

  #
  # Pulldowns
  #
  
  
  def selectChainMolSystem(self, molSystem):
  
    if molSystem is not self.molSystemC:
      self.molSystemC = molSystem
    
  def updateChainMolSysPulldown(self):

    index = 0
    molSystems = self.project.sortedMolSystems()
    molSystem = self.molSystemC
    names = [ms.code for ms in molSystems]
    
    molSystems.append(None)
    names.append('<New>')
    
    if molSystems:
      if molSystem not in molSystems:
        molSystem = molSystems[0]
      
      index = molSystems.index(molSystem)
        
    else:
      molSystem = None
      
    if molSystem is not self.molSystemC:
      self.molSystemC = molSystem
        
    self.chainMolSystemPulldown.setup(names, molSystems, index)
       
  def setChainMolecule(self, molecule):
  
    self.moleculeC = molecule
    self.updateButtons()

  def updateChainMolsPulldown(self):
  
    index = 0
    molecules = self.project.sortedMolecules()
    names = [m.name for m in molecules]
    
    if molecules:
      if self.moleculeC not in molecules:
        self.moleculeC = molecules[0]
      
      index = molecules.index(self.moleculeC)
        
    else:
      self.moleculeC = None
    
    self.chainMolPulldown.setup(names, molecules, index)


  def selectResMolecule(self, molecule):
  
    if molecule is not self.molecule:
      self.molecule = molecule
      self.molResidue = None
      self.updateButtons()
      self.updateMolResiduesAfter()

  def updateMoleculePulldown(self):

    index = 0
    molecules = self.project.sortedMolecules()
    molecule = self.molecule
    names = [m.name for m in molecules]
    
    molecules.append(None)
    names.append('<New>')
    
    if molecules:
      if molecule not in molecules:
        molecule = molecules[0]
      
      index = molecules.index(molecule)
        
    else:
      molecule = None
      
    if molecule is not self.molecule:
      self.molecule = molecule  
      self.molResidue = None
      self.updateButtons()
      self.updateMolResiduesAfter()
    
    
    self.moleculePulldown.setup(names, molecules, index)
    
  #
 
  def updateLinkMoleculePulldown(self, *opt):

    molResidue = self.linkMolResidue

    if molResidue:
      self.linkMolecule = molResidue.molecule
      self.linkCcpCode  = getResidueCode(molResidue)
  
    names = []
    index = 0
    molecules = self.project.sortedMolecules()

    for molecule in molecules:
      names.append(molecule.name)
    
    if molecules:
      if self.linkMolecule not in molecules:
        self.linkMolecule = molecules[0]
      index = molecules.index(self.linkMolecule)
    else:
      self.linkMolecule = None
    
    self.linkMoleculePulldown.setup(names, molecules, index)
    self.updateLinkCcpCodes()
  
  def updateLinkCcpCodes(self):
  
    names = []
    ccpCodes = []
    index = 0
    if self.linkMolecule and self.linkMolecule.molResidues:
      found = set([])
      
      for molResidue in self.linkMolecule.molResidues:
        found.add(getResidueCode(molResidue))
        
      ccpCodes = list(found)
      ccpCodes.sort()
      names = ['<All>',] + ccpCodes
      ccpCodes = [None,] + ccpCodes
      
      if self.linkCcpCode not in ccpCodes:
        self.linkCcpCode = None
      
      index = ccpCodes.index(self.linkCcpCode)
        
    else:
      self.linkCcpCode = None
      
    self.linkCcpCodesPulldown.setup(names, ccpCodes, index)
    self.updateLinkResiduesAfter()
    self.updateLinkMolResidues()
           

  def updateLinkMolResidues(self):

    names = []
    index = 0
    molResidues = []
    cats = []
    
    if self.linkMolecule and self.linkCcpCode:
      cats = None
      
      for mr in self.linkMolecule.sortedMolResidues():
        ccpCode = getResidueCode(mr)
        
        if ccpCode == self.linkCcpCode:
          molResidues.append(mr)
          names.append('%d%s' % (mr.seqCode,ccpCode) )
      
      if molResidues:
        if self.linkMolResidue not in molResidues:
          self.linkMolResidue = molResidues[0]

        index = molResidues.index(self.linkMolResidue)

    elif self.linkMolecule:
      for mr in self.linkMolecule.sortedMolResidues():
        ccpCode = getResidueCode(mr)
        molResidues.append(mr)
        names.append('%d%s' % (mr.seqCode,ccpCode) )
        seq = 10*(mr.seqCode//10)
        cats.append('%d-%d' % (seq,seq+9))
      
      if molResidues:
        if self.linkMolResidue not in molResidues:
          self.linkMolResidue = molResidues[0]

        index = molResidues.index(self.linkMolResidue)

    else:
      self.linkMolResidue = None
    
    self.linkResiduePulldown.setup(names, molResidues, index, categories=cats)
    self.updateMolResLinksAfter()


  def changeLinkMolecule(self, null):

    molecule = self.linkMoleculePulldown.getObject()
    
    if molecule is not self.linkMolecule:
      self.linkMolecule = molecule
      self.updateLinkCcpCodes()

  
  def setLinkCcpCode(self, null):
  
    ccpCode = self.linkCcpCodesPulldown.getObject()
    
    if ccpCode is not self.linkCcpCode:
      self.linkCcpCode = ccpCode
      self.updateLinkResiduesAfter()

  
  def setLinkMolResidue(self, molResidue):
    
    if self.linkMolecule:
      if molResidue is not self.linkMolResidue:
        self.linkMolResidue = molResidue
        self.updateMolResLinksAfter()
  
  def selectMolTypes(self, molTypes):
  
    if molTypes != self.molTypes:
      self.molTypes = molTypes
      self.updateCompoundTable()
  
  #
  

  def setPolymerSeqCcpCode(self, ccpCode):
  
    self.appendButton.config(text='Append "%s"' % ccpCode)

  def getPolymerCcpCodes(self):
  
  
    ccpCodes = getMolTypeCcpCodes(self.polymerMolType, self.project)
    
    return ccpCodes
  
  def changeCodeMode(self, mode):
    
    one, ccp = RES_CODE_MODES
    protein, dna, rna, carbo, dnaRna = POLYMER_MOLTYPES
    
    self.getInputPolymerSeq()
    
    if mode == one:
      if self.polymerMolType == carbo:
        self.seqCodeMode = ccp
        self.seqCodeModeButtons.setIndex(1)
        msg = 'Cannot use 1-letter codes for carbohydrates'
        showWarning('Warning', msg, parent=self)
    
      else:
        self.seqCodeMode = one

    else:
      self.seqCodeMode = ccp

    self.displayPolymerSeq()


  def setPolymerSeqMolType(self, molType):
  
    protein, dna, rna, carbo, dnaRna = POLYMER_MOLTYPES
      
    if molType in (dna, rna, dnaRna):
      if molType != dnaRna:
        self.reverseCompButton.enable()
      else:
        self.reverseCompButton.disable()
      
      if self.polymerMolType in (protein, carbo):
        self.clearPolymerSeq()
      
      else:
        self.getInputPolymerSeq()
      
        if (molType == dna) and (self.polymerMolType in (rna, dnaRna)):
          for i in range(len(self.inputSeq)):
            if self.inputSeq[i][0] == 'U':
              chars = self.inputSeq[i][1:]
              self.inputSeq[i] = 'T' + chars
 
        elif (molType == rna) and (self.polymerMolType in (dna, dnaRna)):
          for i in range(len(self.inputSeq)):
            if self.inputSeq[i][0] == 'T':
              chars = self.inputSeq[i][1:]
              self.inputSeq[i] = 'U' + chars
    
    
    
    else:
      self.reverseCompButton.disable()
      
      if self.polymerMolType != molType:
        self.clearPolymerSeq()
      
      else:
        self.getInputPolymerSeq()

      if molType == carbo:
        self.seqCodeModeButtons.setIndex(1)
       

    self.polymerMolType = molType
    self.displayPolymerSeq()

    ccpCodes = self.getPolymerCcpCodes()
    self.ccpCodeSeqPulldown.setup(ccpCodes, ccpCodes, 0)
  
  def changeSeqMolecule(self, molecule):
  
    if molecule is not self.moleculeS:
      self.moleculeS = molecule
      
  def changeSeqMolSystem(self, molSystem):
  
    if molSystem is not self.molSystemS:
      self.molSystemS = molSystem
  
  def updateSeqDestination(self):
  
    molSystems = self.project.sortedMolSystems()
    molecules = self.project.sortedMolecules()
  
    namesMs = [ms.code for ms in molSystems]
    namesMol = ['%s%s' % (m.name, m.isFinalised and ' (locked)' or '') for m in molecules]
    
    molSystems.insert(0, None)
    molSystems.append(True)
    molecules.append(None)
  
    namesMs.insert(0, '<None>')
    namesMs.append('<New>')
    namesMol.append('<New>')
    
    if self.molSystemS not in molSystems:
      self.molSystemS = molSystems[0]
    
    if self.moleculeS not in molecules:
      self.moleculeS = None

    if self.moleculeComp not in molecules:
      self.moleculeComp = None
  
    indexComp = molecules.index(self.moleculeComp)
    indexMol = molecules.index(self.moleculeS)
    indexMs = molSystems.index(self.molSystemS)
    
    self.seqMolSystemPulldown.setup(namesMs, molSystems, indexMs)
    self.seqMoleculePulldown.setup(namesMol, molecules, indexMol)
    self.compMoleculePulldown.setup(namesMol, molecules, indexComp)
  
  #
  
  def changeCompMolecule(self, molecule):
  
    if molecule is not self.moleculeComp:
      self.moleculeComp = molecule
  
  #
  # Update Functions
  #
  
  def updateButtons(self):
  
    buttons = self.chainButtons.buttons + self.chainButtons2.buttons
    #buttons[3].enable()
    #buttons[4].enable()
    
    if self.chain:
      buttons[0].enable()
      buttons[1].enable()
      buttons[2].enable()
    
    else:
      buttons[0].disable()  
      buttons[1].disable()  
      buttons[2].disable()  
    
    if self.moleculeC:
      buttons[5].enable()
    else:
      buttons[5].disable()
    
    buttons = self.molSysButtons.buttons
    # buttons[0].enable()
    
    if self.molSystem:
      buttons[1].enable()
    else:
      buttons[1].disable()  
    
    
    buttons = self.molResidueButtons.buttons
    # buttons[0].enable()
    # buttons[1].enable()
    
    if self.molResidue:
      if self.molResidue.molecule.isFinalised:
        buttons[3].disable()
      else:
        buttons[3].enable()
        
    else:
      buttons[3].disable()
     
    if self.molecule:
      if self.molecule.isFinalised:
        buttons[0].disable()
        buttons[2].disable()
        buttons[3].disable()
        buttons[4].config(text='Unlock Molecule')
      else:
        buttons[0].enable()
        buttons[2].enable()
        buttons[3].enable()
        buttons[4].config(text='Lock Molecule')
      
      if self.molecule.chains:
        buttons[4].disable()
      else:
        buttons[4].enable()
        
    else:
      buttons[0].enable()
      buttons[2].enable()
      buttons[3].disable()
      buttons[4].disable()
      buttons[4].config(text='Lock Molecule')
    
    buttons = self.moleculeButtons.buttons
    # buttons[0].enable()
  
    if self.moleculeB:
      buttons[0].enable()
      buttons[2].enable()
      buttons[3].enable()
      buttons[4].enable()
      buttons[5].enable()
    
    else:
      buttons[0].disable()  
      buttons[2].disable()  
      buttons[3].disable()  
      buttons[4].disable()  
      buttons[5].disable()  
         

  def updateTab(self, i):
  
    if i in (1,3,5,6):
      self.updateTopLevelAfter()
    elif i == 0:
      self.updateChainsAfter()
    elif i == 2:
      self.updateMolResiduesAfter()
    elif i == 4:
      self.updateLinkMoleculePulldown()
    
    self.updateButtons()
        
  def updateResidueAfter(self, residue):
  
    if self.waitingChain:
      return      
    
    if residue.chain is self.chain:
      self.updateChainsAfter()
 

  def updateMoleculesAfter(self, molecule=None):
  
    if self.waiting:
      return

    self.after_idle(self.updateMolecules)
    

  def updateMolSystemsAfter(self, molSystem=None):
 
    if self.waiting:
      return
    
    self.after_idle(self.updateMolSystems)


  def updateMolDataAfter(self, molecule=None):
  
    if self.waiting:
      return
      
    if molecule is self.moleculeB:
      self.after_idle(self.updateMolData)


  def updateMolResiduesAfter(self, obj=None):
  
    if self.waitingMolRes:
      return

    if (not obj) or (obj.molecule is self.molecule): 
      self.after_idle(self.updateMolResidues)

    if (not obj) or (obj.molecule is self.linkMolecule): 
      self.after_idle(self.updateLinkMolResidues)
    
  
  def updateChainsAfter(self, chain=None):
  
    if chain and chain.isDeleted and (chain is self.chain):
       self.chain = None
  
    if self.waitingChain:
      return
      
    self.waitingChain = True
    
    self.after_idle(self.updateChainMolSysPulldown)
    self.after_idle(self.updateChainMolsPulldown)
    self.after_idle(self.updateChains)
     
  def updateTopLevelAfter(self, obj=None):
   
    if obj and obj.isDeleted:
      if obj is self.molSystem:
        self.molSystem = None
    
      elif obj is self.moleculeB:
        self.moleculeB = None
  
    if self.waiting:
      return
    
    self.waiting = True
    self.after_idle(self.updateTopLevel)  
    
  
  def updateTopLevel(self):
  
    self.updateMolSystems()
    self.updateMolecules()
    self.updateMoleculePulldown()
    self.updateChainMolsPulldown()
    self.updateChainMolSysPulldown()
    self.updateLinkMoleculePulldown()
    self.updateSeqDestination()
    self.waiting = False


  def updateMolSystems(self):
      
    textMatrix = []
    objectList = []
    for molSystem in self.project.sortedMolSystems():
      chains = molSystem.sortedChains()
    
      dict = {}
      for chain in chains:
        dict[chain.molecule.molType or 'mixed'] = True
    
      molTypeText = ' '.join(dict.keys())
        
      datum = [molSystem.code,
               molSystem.name,
               ','.join([c.code for c in chains]),
               '\n'.join([c.molecule.name for c in chains]),
               molTypeText,
               len(molSystem.structureEnsembles),
               molSystem.details or ' ']
 
      textMatrix.append(datum)
      objectList.append(molSystem)
             
    self.molSysMatrix.update(objectList=objectList, textMatrix=textMatrix)
 
  def updateChains(self):
  
    objectList = []
    textMatrix = []
    chains = []
    
    for molSystem in self.project.sortedMolSystems():
      msCode = molSystem.code
      chains = molSystem.sortedChains()
      
      for chain in chains :
        details = chain.details or ' '
        if len(details) > 10:
          details = chain.details[:10] + '...'
      
        datum = [msCode, 
                 chain.code,
                 chain.molecule.name,
                 len(chain.residues),
                 len(chain.chainFragments),
                 chain.molecule.molType,
                 details or ' ']
 
        textMatrix.append(datum)
        objectList.append(chain)
     
    self.chainsMatrix.update(objectList=objectList, textMatrix=textMatrix)
    
    if chains and not self.chain:
      self.chainsMatrix.selectObject(chains[0])
    
    self.updateFragments()
    self.waitingChain = False  
    
  def updateFragments(self): 
    
    objectList = []
    textMatrix = []
     
    if self.chain:
      self.fragmentLabel.setText('Chain "%s" Fragments' % self.chain.code)
      
      fragments = []
      for fragment in self.chain.sortedChainFragments():
        fragments.append( (fragment.findFirstResidue().seqCode,fragment) )
      
      fragments.sort()
      objectList = [x[1] for x in fragments]
        
      for fragment in objectList:
        linPolText = 'No'
        if fragment.isLinearPolymer:
          linPolText = 'Yes'
        residues = fragment.residues
        datum = [fragment.serial,
                 fragment.molType,
                 len(residues),
                 residues[0].seqCode,
                 linPolText,
                 self.makeSequenceSnippet(fragment)]
 
        textMatrix.append(datum)
    else:
      self.fragmentLabel.setText('Chain Fragments')
      
     
    self.fragmentMatrix.update(objectList=objectList, textMatrix=textMatrix)
  
  
  def updateMolResidues(self):
  
    editWidgets      = [None] * 6 
    editGetCallbacks = [None] * 6 
    editSetCallbacks = [None] * 6 
        
    objectList = []
    textMatrix = []
    
    chainText = 'Chains with this sequence: '
    
    if self.molecule:
      chainCodes = ['%s:%s' % (c.molSystem.code, c.code) for c in self.molecule.chains]
      chainCodes.sort()
      chainText += ', '.join(chainCodes)
    
      for molResidue in self.molecule.sortedMolResidues():
      
        linkMolResidues = []
        for le in molResidue.molResLinkEnds:
          if le.molResLink:
            for le2 in le.molResLink.molResLinkEnds:
              if le2 is not le:
                linkMolResidue = le2.molResidue
                linkMolResidues.append((linkMolResidue.seqCode, linkMolResidue))
        
        linkMolResidues.sort()
        linkText = ' '.join(['%d%s' % (seqCode,getResidueCode(mr)) for (seqCode, mr) in linkMolResidues])
        
        datum = [molResidue.seqCode,
                 getResidueCode(molResidue),
                 self.getDescriptorStereoName(molResidue.chemCompVar),
                 molResidue.molType,
                 molResidue.linking,
                 linkText]
 
        objectList.append(molResidue)
        textMatrix.append(datum)
      
      if not self.molecule.isFinalised:
      
        editWidgets = [None, self.ccpCodePulldown,
                       self.descriptorPulldown, None,
                       self.linkingPulldown, None]
        editGetCallbacks = [None, self.getCcpCode,
                            self.getDescriptor, None,
                            self.getLinking, self.showResLinks]
        editSetCallbacks = [None, self.setCcpCode,
                            self.setDescriptor, None,
                            self.setLinking, None]
    
    self.molResChainInfoLabel.set(chainText)
 
    self.molResidueMatrix.update(objectList=objectList, textMatrix=textMatrix,
                                 editSetCallbacks=editSetCallbacks,
                                 editGetCallbacks=editGetCallbacks,
                                 editWidgets=editWidgets)
  
    self.waitingMolRes = False
  

  def updateMolecules(self):
  
    objectList = []
    textMatrix = []
    
    for molecule in self.project.sortedMolecules():
      chains = molecule.chains 
      chainCodes = ['%s:%s' % (c.molSystem.code, c.code) for c in chains ]
      chainCodes.sort()
           
      datum = [molecule.name,
               molecule.longName or ' ',
               molecule.molType,
               ', '.join(chainCodes),
               molecule.details or ' ',
               molecule.seqDetails or ' ']
      #         ' '.join([x for x in molecule.commonNames]),
      #         ' '.join([x for x in molecule.keywords]),
      #         ' '.join([x for x in molecule.functions]),
      
      objectList.append(molecule)
      textMatrix.append(datum)
    
    self.moleculeMatrix.update(objectList=objectList,
                               textMatrix=textMatrix)
  
    self.updateMolData()

  def updateMolData(self):

    objectList = []
    textMatrix = []
    
    molecule = self.moleculeB
    
    if molecule:    
            
      isAromatic = 'No'
      if molecule.isAromatic:
        isAromatic = 'Yes'
      
      isParamagnetic = 'No'
      if molecule.isParamagnetic:
        isParamagnetic = 'Yes'
      
      hasNonStdChemComp = 'No'
      if molecule.hasNonStdChemComp:
        hasNonStdChemComp = 'Yes'

      hasNonStdChirality = 'No'
      if molecule.hasNonStdChirality:
        hasNonStdChirality = 'Yes'

      isStdLinear = 'No'
      if molecule.isStdLinear:
        isStdLinear = 'Yes'

      isStdCyclic = 'No'
      if molecule.isStdCyclic:
        isStdCyclic = 'Yes'

      datum = molecule.name
      objectList.append(datum)
      textMatrix.append(['Name',datum])

      datum = molecule.seqLength
      objectList.append(datum)
      textMatrix.append(['Sequence Length',datum])

      datum = molecule.calcIsoelectricPoint
      objectList.append(datum)
      textMatrix.append(['Isoelectric Point',datum])

      datum = molecule.empiricalFormula
      objectList.append(datum)
      textMatrix.append(['Empirical Formula',datum])

      datum = molecule.molecularMass
      objectList.append(datum)
      textMatrix.append(['Molecular Mass',datum])

      datum = molecule.formalCharge
      objectList.append(datum)
      textMatrix.append(['Formal Charge',datum])

      datum = isAromatic
      objectList.append(datum)
      textMatrix.append(['Has Aromatics',datum])

      datum =  isParamagnetic
      objectList.append(datum)
      textMatrix.append(['Is Paramagnetic',datum])

      datum = molecule.seqString or '' # NBNB TBD fix this properly
      objectList.append(datum)
      if len(datum) > 10:
        seq = datum[:10] + '...'
      else:
        seq = datum
      objectList.append(datum)
      textMatrix.append(['Sequence',seq])

      datum = molecule.stdSeqString or '' # NBNB TBD fix this properly
      if len(datum) > 10:
        seq = datum[:10] + '...'
      else:
        seq = datum
      objectList.append(datum)
      textMatrix.append(['Std Sequence',seq])

      datum = hasNonStdChemComp
      objectList.append(datum)
      textMatrix.append(['Has Non-Standard Residues',datum])

      datum = hasNonStdChirality
      objectList.append(datum)
      textMatrix.append(['Has Non-Standard Chirality',datum])
  
      datum = isStdLinear
      objectList.append(datum)
      textMatrix.append(['Is Standard Linear',datum])

      datum = isStdCyclic
      objectList.append(datum)
      textMatrix.append(['Is Standard Cyclic',datum])
    
    self.molDataMatrix.update(objectList=objectList, textMatrix=textMatrix)

  #

  def updateLinkResiduesAfter(self, object=None):

    if object and object.molecule is not self.linkMolecule:
      return
      
    self.after_idle(self.updateLinkMolResidues) 
          
  def updateMolResLinksAfter(self, object=None):

    if self.waitingMolResLink:
      return
    
    if object:
      if object.molecule is not self.linkMolecule:
        return

      molResidues = [le.molResidue for le in object.molResLinkEnds]
      if self.linkMolResidue not in molResidues:
        return
 
    self.waitingMolResLink = True
    self.after_idle(self.updateMolResLinks)

  def getLinkEditData(self):

    if self.linkMolResidue.molecule.isFinalised:
      editWidgets      = [ None, None, None ]
      editGetCallbacks = [ None, None, None ]
      editSetCallbacks = [ None, None, None ]

    else:
      editWidgets      = [ None, self.nonSeqLinkPulldown, None ]
      editGetCallbacks = [ None, self.getMolResLink,      None ]
      editSetCallbacks = [ None, self.setMolResLink,      None ]

    return (editWidgets, editGetCallbacks, editSetCallbacks)
    
  def updateMolResLinks(self):
    
    objectList = []
    textMatrix = []
    editWidgets = []
    editGetCallbacks = []
    editSetCallbacks = []
    
    if self.linkMolResidue:
    
      for molResLinkEnd in self.linkMolResidue.molResLinkEnds:
        molResLink = molResLinkEnd.molResLink
        otherResidue = '<None>'
        otherCode    = ''
        if molResLink:
          if molResLinkEnd.linkCode: #  not in ('prev','next'):
            objectList.append(molResLinkEnd)
          else:
            objectList.append(None)
            
          for le in molResLink.molResLinkEnds:
            if le is not molResLinkEnd:
              mr = le.molResidue
              otherResidue = '%d%s' % (mr.seqCode,getResidueCode(mr))
              otherCode = le.linkCode
        else:
          objectList.append(molResLinkEnd)
          
        datum = []
        datum.append(molResLinkEnd.linkCode)
        datum.append(otherResidue)
        datum.append(otherCode)
        
        textMatrix.append(datum)
        
      (editWidgets, editGetCallbacks, editSetCallbacks) = self.getLinkEditData()
    
    self.molResLinkMatrix.update(textMatrix=textMatrix,
                                 objectList=objectList,
                                 editSetCallbacks=editSetCallbacks,
                                 editGetCallbacks=editGetCallbacks,
                                 editWidgets=editWidgets)
                                 
    self.waitingMolResLink = False
      
  #
  
 
  def updateCompoundTable(self):

    project = self.project
    findChemComp = project.findFirstChemComp
    textMatrix = []
    objects = []
    colorMatrix = []
    for molType in self.molTypes:     
      dictCC = getChemCompOverview(molType, project)
        
      ccpCodes = dictCC.keys()
      ccpCodes.sort()
      for ccpCode in ccpCodes:
        (code1Letter, code3Letter, name, mf) = dictCC.get(ccpCode, ('?','Error',None,None))
        if not name:
          continue
          
        if findChemComp(molType=molType, ccpCode=ccpCode):
          colors = ['#B0B0FF'] * 3
        else:
          colors = [None] * 3
          
        name = re.sub('(.{50}\w+\W)',r'\1\n',name)
        textMatrix.append( [ccpCode, code3Letter, name])
        objects.append((ccpCode, molType))
        colorMatrix.append(colors)
    
    self.compoundMatrix.update(objectList=objects,
                               textMatrix=textMatrix,
                               colorMatrix=colorMatrix)
   

  #
  # Actions
  #

  def addSequence(self):
  
    self.moleculeS = None
    self.molSystemS = True
    self.updateSeqDestination()
    
    self.tabbedFrame.select(5)
    

  def showChainSequence(self):
  
    if self.chain:
      self.molecule = self.chain.molecule
      self.updateMoleculePulldown()
      self.tabbedFrame.select(2)

  def chainAddSequence(self):

    if self.chain:
      molSystem = self.chain.molSystem
    else:
      molSystem = True
  
    self.moleculeS = None
    self.molSystemS = molSystem
    self.updateSeqDestination()
    
    self.tabbedFrame.select(5)

  def deleteChain(self, *event):
  
    chain = self.chain
  
    if chain:
      if isChainAssigned(chain):
        msg = 'Cannot delete: chain %s carries assignments'
        showWarning('Warning',msg % (chain.code), parent=self)
        return

      msg = 'Really delete mol system %s chain %s?'
      if showOkCancel('Confirm',msg % (chain.molSystem.code,chain.code), parent=self):
        self.chain = None
        
        molecule = chain.molecule
        chain.delete()
        
        if not molecule.chains:
          msg = 'Molecule %s now has no chains. Delete this too?'
          if showYesNo('Query', msg % molecule.name, parent=self):
            if molecule is self.molecule:
              self.molecule = None
            
            if molecule is self.moleculeB:
              self.moleculeB = None
            
            if molecule is self.moleculeC:
              self.moleculeC = None
            
            if molecule is self.moleculeS:
              self.moleculeS = None
            
            for molComponent in molecule.molComponents:
              molComponent.delete()
              
            molecule.delete()
      
  def copyChain(self):
  
    chain = self.chain
    if self.chain:
      makeChainCopy(self.chain)

        
  def makeChain(self, molSystem=None, molecule=None):

    if not molSystem:
      if not self.molSystemC:
        self.molSystemC = self.addMolSys()
        
      molSystem = self.molSystemC
    
    if not molecule:
      if not self.moleculeC:
        self.moleculeC = self.addMolecule()
        
      molecule = self.moleculeC


    if molSystem and molecule:        
      code = nextChainCode(molSystem)
      code = askString('Request','Chain Code:', code, parent=self)

      if not code:
        return
    
      if len(code.splitlines()) > 1:
        showWarning('Failure', 'Code cannot contain a line break' ,parent=self)
        return
    
      if molSystem.findFirstChain(code=code):
        showWarning('Failure', 'Code already used', parent=self)
        return
         	 	
      return makeChain(molSystem, molecule, code=code)

  #

  def addMolSys(self):
  
    project = self.project
    N = 1
    while project.findFirstMolSystem(code='MS%d' % (N)):
      N += 1
    
    name = 'MS%d' % N
    name = askString('Request','Mol System Code:',name,parent=self)

    if not name:
      return
    
    for character in WHITESPACE:
      if character in name:
        showWarning('Failure','Code cannot contain whitespace',parent=self)
        return

    if project.findFirstMolSystem(code=name):
      showWarning('Failure','Code already used',parent=self)
      return
    
    molSystem = project.newMolSystem(code=name)
    molSystem.name = molSystem.code
    
    return molSystem

  def delMolSys(self, *event):
  
    molSystem = self.molSystem
  
    if not molSystem:
      return
    
    if len(self.molSystem.chains) > 0:
      msg = 'Molecular System has chains. Delete these first.'
      showWarning('Failed', msg, parent=self)
      return
         
    self.molSystem = None  
    molSystem.delete()

  def editMolecule(self, molecule=None):
  
    if (not molecule) and self.chain:
      molecule = self.chain.molecule
  
    if molecule:
      self.moleculeMatrix.selectObject(molecule)
      
    self.tabbedFrame.select(3)

  #

  def showMoleculeSequence(self):
  
    if self.moleculeB:
      self.molecule = self.moleculeB
      self.updateMoleculePulldown()
      self.tabbedFrame.select(2)
  
  def toggleMoleculeLock(self):
  
    if self.molecule:
      boolean = self.molecule.isFinalised
      
      if boolean and self.molecule.chains:
        showWarning('Failure','Cannot unlock a molecule with chains.',parent=self)
        return
        
      self.molecule.isFinalised = not boolean
      self.updateButtons()
      self.updateMolResiduesAfter()

  def showResLinks(self, *opt):
 
    if self.molResidue:
      self.linkMolResidue = self.molResidue
      self.updateLinkMoleculePulldown()
      self.tabbedFrame.select(4)
      
  def addMolPolymer(self):
   
    molecule = self.moleculeB
    if not molecule:
      molecule = self.addMolecule()
    
    self.moleculeS = molecule
    self.molSystemS = None
    
    self.updateSeqDestination()
    
    self.tabbedFrame.select(5)
   
  def addResPolymer(self):
    
    molecule = self.molecule
    if not molecule:
      molecule = self.addMolecule()
    
    self.moleculeS = molecule
    self.molSystemS = None
    
    self.updateSeqDestination()
    
    self.tabbedFrame.select(5)
  
  
  def copyResSequence(self):
  
    if self.molecule:
      molResidues = self.molecule.sortedMolResidues()
      seq = [mr.ccpCode for mr in molResidues]
      one, ccp = RES_CODE_MODES

      molType = molResidues[0].molType
      
      self.clearPolymerSeq()
      self.polymerMolType = molType
      self.polymerMolTypePulldown.set(molType)
      self.seqCodeMode = ccp
      self.seqCodeModeButtons.setIndex(1)
      self.seqTextBox.setText(' '.join(seq))
      
      self.moleculeS = None
      self.molSystemS = None
      self.updateSeqDestination()
      self.tabbedFrame.select(5)
  
  def addMoleculeCompound(self):
   
    molecule = self.moleculeB
    if not molecule:
      molecule = self.addMolecule()

    self.moleculeComp = molecule
    self.updateSeqDestination()
    self.tabbedFrame.select(6)

  def addResCompound(self):
    
    molecule = self.molecule
    if not molecule:
      molecule = self.addMolecule()
    
    self.moleculeComp = molecule
    self.updateSeqDestination()
    self.tabbedFrame.select(6)
  
  def addMolecule(self, molName=None):
  
    i = len(self.project.molecules) + 1
    molName = molName or 'Molecule %d' % (i)
    while self.project.findFirstMolecule(name=molName):
      i += 1
      molName = 'Molecule %d' % (i)

    name = askString('Input text', 'New Molecule Name:', molName, parent=self)
    if name:
      if self.project.findFirstMolecule(name=name):
        showWarning('Failure','Molecule name already in use.',parent=self)
        return
        
      else:
        if len(name) > 79:
          name = name[:76] + '...'
        molecule = self.project.newMolecule(name=name)
        
    else:
      return
  
    return molecule
  
  def copyMolecule(self, molName=None):
  
    if self.moleculeB:
      i = len(self.project.molecules) + 1
      molName = molName or 'Molecule %d' % (i)
      while self.project.findFirstMolecule(name=molName):
        i += 1
        molName = 'Molecule %d' % (i)

      name = askString('Input text','New Molecule Name:',molName,parent=self)
 
      if name:
        if self.project.findFirstMolecule(name=name):
          showWarning('Failure','Molecule name already in use.',parent=self)
        else:
          molecule = copyMolecule(self.moleculeB, newName=name)
          molecule.isFinalised = False
          return molecule


  def deleteMolecule(self, *event):
  
    molecule = self.moleculeB
  
    if molecule:
      if molecule.isFinalised:
        showWarning('Warning', 'Cannot delete. Molecule is locked.', parent=self)
        return
    
      if molecule.chains:
        for chain in molecule.chains:
          if isChainAssigned(chain):
            msg = 'Cannot delete. Molecule has assigned chains'
            showWarning('Warning', msg, parent=self)
            return
        
        chainText = ' '.join([c.code for c in molecule.chains])
        
        msg = 'Are you sure you want to delete molecule: %s\n and its chains [%s]'
        if showOkCancel('Confirm', msg % (molecule.name,chainText), parent=self):
          for chain in molecule.chains:
            chain.delete()

          self.moleculeB = None
          for molComponent in molecule.molComponents:
            molComponent.delete()
          
          molecule.delete()
        
        return
      
      msg =  'Are you sure you want to delete molecule: %s'
      if showOkCancel('Confirm', msg % molecule.name, parent=self):
        self.moleculeB = None
        for molComponent in molecule.molComponents:
          molComponent.delete()
        
        molecule.delete()
  
  def importXml(self):
    
    if self.project:
      from os import path
      from memops.format.xml import XmlIO
      from memops.general.Implementation import ApiError
      
      dataRepository = self.project.findFirstRepository(name='userData')
      projDir = dataRepository.url.path
      projParentDir = path.dirname(projDir)
 
      fileTypes = [ FileType('CCPN XML', ['*.xml']), ]
      fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                                        title='Import ChemComp XML File',
                                        dismiss_text='Cancel',
                                        selected_file_must_exist=True,
                                        multiSelect=True,
                                        directory=projParentDir)

      filePaths = fileSelectPopup.file_select.getFiles()
      
      molTypes = set(POLYMER_MOLTYPES + ['other'])
      chemComp = None

      for filePath in filePaths:
        dirName, fileName = path.split(filePath)
        xmlNameComponents = fileName.split('+')
        
        if len(xmlNameComponents) != 3:
          msg = 'Non-standard CCPN XML file name "%s".\n' % fileName
          msg += 'ChemComp XML file names should be of the form "molType+ccpCode+GUID".\n'
          msg += '(If the file was renamed try the original name.)\n'
          showWarning('Import Failure', msg, parent=self)
          return
        
        molType, ccpCode, guid = xmlNameComponents
        
        if molType not in molTypes:
          msg = 'Molecule type "%s" not used by CCPN.\n' % (molType)
          msg += 'Valid types are:\n'
          msg += ', '.join(sorted(list(molTypes)))
          showWarning('Import Failure', msg, parent=self)
          return
        
        if self.project.findFirstChemComp(molType=molType, ccpCode=ccpCode):
          msg = 'Code "%s" for molecule type "%s" is already in use.\n' % (ccpCode, molType)
          msg += 'Compound codes must be unique for a given molecule type.\n'
          msg += 'This compound may already be available without an import.\n'
          showWarning('Import Failure', msg, parent=self)
          return
      
        try:
          chemComp = self.project.importData(filePath)      
          XmlIO.save(projDir, chemComp)
          
        except ApiError, e:
          msg = 'XML file could not be loaded, it may be corrupt.'
          msg += '\n\nUnderlying error:\n\n%s' % e
          showWarning('Import Failure', msg, parent=self)
          return
      
      if chemComp: # Last loaded
        self.molTypes = (chemComp.molType,)
        self.molTypePulldown.set(chemComp.molType)
        self.updateCompoundTable()
        key = (chemComp.ccpCode,chemComp.molType)
        self.compoundMatrix.selectObject(key) 


  def commitCompound(self):
  
    if self.chemCompVar:
      chemComp = self.chemCompVar.chemComp
    
      ccpCode = chemComp.ccpCode
      molType = chemComp.molType
      molecule = self.moleculeComp
      
      if not molecule:
        molecule = self.addMolecule(chemComp.name)
     
      if molecule: # May have cancelled above
      
        if molecule.isFinalised:
          msg = 'Cannot add residue to a locked molecule.'
          showWarning('Failure', msg, parent=self)
          return
 
        addMolResidues(molecule,molType,[ccpCode,],startNum=1)
 
        self.molecule = molecule
        self.updateMoleculePulldown()
        self.tabbedFrame.select(2)

  def appendPolymerResidue(self):
  
    ccpCode = self.ccpCodeSeqPulldown.getText()
     
    self.getInputPolymerSeq()
    self.inputSeq.append(ccpCode)
    self.displayPolymerSeq()

  def reverseComplement(self):

    if self.polymerMolType not in POLYMER_MOLTYPES[1:3]:
      return # safety, should not be here unless above true in any case

    if self.polymerMolType == POLYMER_MOLTYPES[1]:
      resMap = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G' }
    else:
      resMap = { 'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G' }

    self.getInputPolymerSeq()
    firstLetters = [ s[0] for s in self.inputSeq ]
    for s in firstLetters:
      if not resMap.has_key(s):
        msg =  'Cannot use reverse complement button when sequence contains %s'
        showWarning('Warning', msg % repr(s), parent=self)
        return

    outSeq = []
    firstLetters.reverse()
    for n in range(len(self.inputSeq)):
      s = resMap[firstLetters[n]] + self.inputSeq[n][1:]
      outSeq.append(s)

    self.inputSeq = outSeq
    self.displayPolymerSeq()

  #
      
  def readFile(self):
  
    from memops.gui.FileSelect import FileType
    from memops.gui.FileSelectPopup import FileSelectPopup

    fileSelectPopup = FileSelectPopup(self,
                                      title='Read sequence file',
                                      dismiss_text='Cancel',
                                      selected_file_must_exist=True)

    fileName = fileSelectPopup.getFile()
        
    if fileName:
      file = open(fileName, 'r')
      text = ''
      line = file.readline()
      while line:
        text += line
        line = file.readline()
      
      # Could maybe guess codeMode
      self.seqTextBox.setText( text )

  #
  # Table edit functions
  #

  def getChainDetl(self, chain):
  
    if chain and chain.details:
      self.editChainDetlEntry.set(chain.details)
  
  def setChainDetl(self, event):
  
    value = self.editChainDetlEntry.get()
    if value:
      self.chain.details = value
  
  # 
  
  def getStartNum(self, fragment):
  
    if fragment:
      # NB fragment.residues is an ordered list
      residues = fragment.residues
      if residues:
        start = residues[0].seqCode
        self.editStartNumlEntry.set(start)
    
  def setStartNum(self, event):

    self.configure(cursor="watch")
    self.update_idletasks()
  
    start = self.editStartNumlEntry.get()
    if (start is not None) and self.fragment:
      residues = self.fragment.residues
      initial  = residues[0].seqCode
      if start != initial:
        delta   = start-initial
        end     = residues[-1].seqCode+delta
        residueS = self.chain.findFirstResidue(seqCode=start)
        residueE = self.chain.findFirstResidue(seqCode=end)

        for i in range(start, end+1):
          testResidue = self.chain.findFirstResidue(seqCode=i)
          if testResidue and (testResidue not in residues):
            showWarning('Failure', 'This numbering would overlap with residues from other fragments.', parent=self)
            self.editStartNumlEntry.set(initial)
            return
          
        for residue in residues:
          i = residue.seqCode
          residue.seqCode = i + delta
    
    self.after_idle(lambda: self.configure(cursor="arrow"))
  
  #

  def getMolSysName(self, molSystem):
  
    if molSystem and molSystem.name:
      self.editMolSysNameEntry.set(molSystem.name)
  
  def setMolSysName(self, event):
  
    value = self.editMolSysNameEntry.get()
    if value:
      self.molSystem.name = value
  
  def getMolSysDetl(self, molSystem):
  
    if molSystem and molSystem.details:
      self.editMolSysDetlEntry.set(molSystem.details)
  
  def setMolSysDetl(self, event):
  
    value = self.editMolSysDetlEntry.get()
    if value:
      self.molSystem.details = value

  #
  

 
  def getCcpCode(self, molResidue):
  
    index = 0
    ccpCodes = getMolTypeCcpCodes(molResidue.molType, self.project)

    if ccpCodes:
      ccpCode = getResidueCode(molResidue)
      
      if ccpCode in ccpCodes:
        index = ccpCodes.index(ccpCode)
        
    self.ccpCodePulldown.setup(ccpCodes, ccpCodes, index)

  def setCcpCode(self, null):
  
    ccpCode = self.ccpCodePulldown.getText()
      
    if self.molResidue:
      if ccpCode != getResidueCode(self.molResidue):
        self.molResidue = setMolResidueCcpCode(self.molResidue,ccpCode)
        self.updateMolDataAfter()
  
  def getDescriptor(self, molResidue):

    index = 0
    data = self.getDescriptors(molResidue)
    names = [x[0] for x in data]
    descriptors = [x[1] for x in data]
    
    if names:
      index = descriptors.index(molResidue.chemCompVar.descriptor)
    
    self.descriptorPulldown.setup(names, descriptors, index)


  def setDescriptor(self, null):
  
    descriptor = self.descriptorPulldown.getObject()
    molResidue = self.molResidue
       
    if molResidue: 
      if not self.checkMoleculeMutable(molResidue.molecule):
        return
    
      chemComp = molResidue.chemComp
      
      if descriptor != molResidue.chemCompVar.descriptor:
        chemCompVar = chemComp.findFirstChemCompVar(descriptor=descriptor,
                                                    linking=molResidue.chemCompVar.linking)
        if not chemCompVar:
          chemCompVar = chemComp.findFirstChemCompVar(descriptor=descriptor)
        
        if chemCompVar:
          self.molResidue = setMolResidueChemCompVar(molResidue,chemCompVar)
          self.updateMolDataAfter()
          self.updateMolResiduesAfter()
 

  def getLinking(self, molResidue):

    index = 0
    names = self.getLinkings(molResidue)
    
    if names:
      index = names.index(molResidue.chemCompVar.linking)
    
    self.linkingPulldown.setup(names, names, index)

  def checkMoleculeMutable(self, molecule):
    
    chains = molecule.sortedChains()
    
    if chains:
      assigned = []
      for chain in chains:
        if isChainAssigned(chain):
          assigned.append(chain)
      
      if assigned:
        chainText = ['%s:%s' % (c.molSystem.code, c.code) for c in assigned]
        
        msg = 'Cannot change molecule residue state. '
        msg += 'Molecule is associated with assigned chains: %s ' % chainText
        msg += 'You can make a new molecule and transfer assignments.'
        showWarning('Failure', msg, parent=self)
        return False
      
      else:
        chainText = ['%s:%s' % (c.molSystem.code, c.code) for c in chains]
        
        msg = 'Cannot change molecule residue state. '
        msg += 'Molecule is associated with  chains: %s ' % chainText
        msg += 'delete these first.'
        showWarning('Failure', msg, parent=self)
        return False

    return True

  def setLinking(self, null):
  
    linking = self.linkingPulldown.getObject()
    molResidue = self.molResidue
    
    if molResidue:
      if not self.checkMoleculeMutable(molResidue.molecule):
        return  
        
      chemCompVar = molResidue.chemCompVar
      
      if linking != chemCompVar.linking:
        descriptor = chemCompVar.descriptor
        chemCompVar = molResidue.chemComp.findFirstChemCompVar(linking=linking)
        
        if chemCompVar:
          self.molResidue = setMolResidueChemCompVar(molResidue,chemCompVar)
          self.updateMolDataAfter()
          self.updateMolResiduesAfter()
  #
  
  def setSeqDetails(self, event):

    text = self.seqMolDetailsEntry.get()
    if text:
      self.moleculeB.seqDetails = text
 
  def getSeqDetails(self, molecule):

    if molecule:
      self.seqMolDetailsEntry.set(molecule.seqDetails)

  def setMolDetails(self, event):

    text = self.molDetailsEntry.get()
    if text:
      self.moleculeB.details = text
 
  def getMolDetails(self, molecule):

    if molecule:
      self.molDetailsEntry.set(molecule.details)

  def setLongName(self, event):

    text = self.longNameEntry.get()
    if text:
      self.moleculeB.longName = text
 
  def getLongName(self, molecule):

    if molecule:
      self.longNameEntry.set(molecule.longName)
  
  #
  
  def getMolResLink(self, molResLinkEnd):
  
    names = ['<None>',]
    links = [None,]
    index = 0
    
    linkEnds = []
    linkCode = molResLinkEnd.linkCode
    linkResidue = None
    if molResLinkEnd.molResLink:
      for molResLinkEnd2 in molResLinkEnd.molResLink.molResLinkEnds:
        if molResLinkEnd2 is not molResLinkEnd:
          linkResidue = molResLinkEnd2.molResidue
          linkEnds.append(molResLinkEnd2)
          index = 1

    linkEnds += self.getFreeMolResidueLinks(self.linkMolecule)
    for linkEnd in linkEnds:
      if (linkEnd.linkCode == linkCode) and (linkCode in ('prev','next')):
        continue
    
      molResidue = linkEnd.molResidue
      names.append('%d%s %s' % (molResidue.seqCode, getResidueCode(molResidue), linkEnd.linkCode))
      links.append(linkEnd)
      
    
    self.nonSeqLinkPulldown.setup(names, links, index)
  
  def setMolResLink(self, null):
    
    
    linkEnd = self.nonSeqLinkPulldown.getObject()

    if self.linkMolResidue and self.linkMolecule:
      selLinkEnd = self.molResLinkEnd
      
      if (linkEnd is None) and selLinkEnd and selLinkEnd.molResLink:
        selLinkEnd.molResLink.delete()
    
      elif linkEnd:
        molResidue = linkEnd.molResidue
        linkCode = linkEnd.linkCode
        
        if molResidue is self.linkMolResidue: # Shouldn't come to this anyway...
          return
        
        if selLinkEnd and selLinkEnd.molResLink:
          if selLinkEnd.molResLink.findFirstMolResLinkEnd(molResidue=molResidue,
                                                          linkCode=linkCode):
            return

        linkEndB = molResidue.findFirstMolResLinkEnd(linkCode=linkCode)
        if selLinkEnd and linkEndB:
          if selLinkEnd.molResLink:
            selLinkEnd.molResLink.delete()
          
          if linkEndB.molResLink:
            linkEndB.molResLink.delete()
          
          molResLink = self.linkMolecule.newMolResLink(molResLinkEnds=(selLinkEnd,linkEndB))

  
  #
  # Utility
  #

  def makeSequenceString(self, chain):
  
    protein, dna, rna, carbo, dnaRna = POLYMER_MOLTYPES
  
    string = ''
    residues = chain.sortedResidues()
    i = residues[0].seqId
    for j in range( (i-1) % 10 ):
      string = string + '    '
    
    molType = chain.molecule.molType
    for res in chain.molecule.sortedMolResidues():
      if molType in (dna, rna, dnaRna):
        code = res.chemCompVar.chemComp.code1Letter
      else:
        code = res.ccpCode
      string = string + '%-3.3s ' % code
      if ( i%10 == 0):
        string = string + ' %4.1d\n' % i
      i += 1
    
    return string


  
  def makeSequenceSnippet(self,fragment):
  
    residues = fragment.residues
    seq = ' '.join([r.ccpCode for r in residues[:10]])
    
    return seq + '...'
 
  def getLinkings(self, molResidue):
  
    #descriptor = molResidue.chemCompVar.descriptor
    linkings   = set()
    for chemCompVar in molResidue.chemComp.chemCompVars:
      linkings.add( chemCompVar.linking )
    
    linkings = list(linkings)
    linkings.sort()
    
    return linkings

  def getDescriptors(self, molResidue):
  
    linking = molResidue.chemCompVar.linking
    descriptors = set()
    
    for chemCompVar in molResidue.chemComp.findAllChemCompVars(linking=linking):
      name = self.getDescriptorStereoName(chemCompVar)
      descriptors.add( (name, chemCompVar.descriptor) )
    
    descriptors = list(descriptors)
    descriptors.sort()
    
    return descriptors

  def getDescriptorStereoName(self, chemCompVar):

    chemComp = chemCompVar.chemComp
    stereochemistries = chemComp.stereochemistries
    if stereochemistries:
      chemAtoms = chemCompVar.chemAtoms
      stereos   = []
      
      for stereochemistry in stereochemistries:
        
        use = True
        for chemAtom in stereochemistry.chemAtoms:
          if chemAtom not in chemAtoms:
            use = False
            break

        if use:
          stereoClass = stereochemistry
        
          if chemComp.molType == 'carbohydrate':
            if stereoClass[:8] == 'stereo_1':
              stereoClass = 'alpha'
            elif stereoClass[:8] == 'stereo_2':
              stereoClass = 'beta'
        
          stereos.append(stereoClass)

      name = ' '.join(stereos)
      name += chemCompVar.descriptor
    else:
      name = chemCompVar.descriptor
            
      if chemComp.molType == 'carbohydrate':
        if name[:8] == 'stereo_1':
          name = 'alpha' + name[8:]
        elif name[:8] == 'stereo_2':
          name = 'beta' + name[8:]     

    return name


  def getFreeMolResidueLinks(self, molecule):
  
    links = []
    for molResidue in molecule.sortedMolResidues():
      if molResidue is self.linkMolResidue:
        continue

      for linkEnd in molResidue.molResLinkEnds:
        if not linkEnd.molResLink:
          links.append( linkEnd )

    return links 
                  
  def checkSeq(self, seq, molType, checkChemComps=False):

    if molType not in POLYMER_MOLTYPES:
      showWarning('Warning', 'Invalid chain type', parent=self)
      return False, []
    
    project = self.project
    findChemComp = project.findFirstChemComp
    
    isValid = True
    invalid = [] 
    
    unknowns = set()
    knownCcpCodes = getMolTypeCcpCodes(molType=molType)
    
    n = len(seq)
    for i in range(n):
      if checkChemComps or (seq[i] not in knownCcpCodes):
        resCode = seq[i]
      
        chemComp = findChemComp(molType=molType, ccpCode=resCode)
        
        if not chemComp:
          chemComp = findChemComp(molType=molType, code3Letter=resCode)

        if not chemComp:
          chemComp = findChemComp(molType=molType, code3Letter=resCode.upper())
       
        if not chemComp:
          chemComp = findChemComp(molType='other', ccpCode=resCode)
        
        
        if not chemComp:
          unknowns.add(seq[i])
          invalid.append(i)
          isValid = False
        
        elif (i < n-1) and not chemComp.findFirstLinkEnd(linkCode='next'):
          invalid.append(i)
          isValid = False
        
        elif (i > 1) and not chemComp.findFirstLinkEnd(linkCode='prev'):
          invalid.append(i)
          isValid = False
            
        else:
          seq[i] = chemComp.ccpCode
    
    if not isValid:
      unknowns = list(unknowns)
      unknowns.sort()
      unknown = ','.join(unknowns)

      msg = 'Chemical components %s not available locally. OK to download to project directory?' % unknown
      if showOkCancel('Query', msg, parent=self):        
        for ccpCode in unknowns[:]:
          chemComp = getChemComp(project, molType, ccpCode, download=True)
 
          if chemComp:
            unknowns.remove(ccpCode)
    
        if not unknowns:
          isValid = True
    
        if not isValid:
          unknown = ','.join(unknowns)
          msg = 'Residue code(s) %s cannot be found for molecule type %s.'
          #if molType == 'protein':
          #  msg += ' Substituted code XXX'
 
          showWarning('Warning', msg % (unknown,molType) , parent=self)
 
    invalid.reverse()
    #for i in invalid:
    #  del seq[i]
    
    return isValid, seq

  #
  # Polymer Seq Entry Functions
  #
  def setPolymerStart(self):
  
    i = self.polymerStartEntry.get()
    if i is None:
      i = 0
      
    self.seqStartNum = i

  def getInputPolymerSeq(self, checkChemComps=False):
  
    # DNA/RNA mixed polymer types TBD
    protein, dna, rna, carbo, dnaRna = POLYMER_MOLTYPES
    one, ccp = RES_CODE_MODES
    
    #defaultResidue = {protein:'XXX'}
    
    text = self.seqTextBox.getText()
    seq = []
    
    if self.seqCodeMode == one:
      molType = self.polymerMolType
 
      dict = code1LetterToCcpCodeDict.get(molType, {})
      #if molType != carbo:
      #  dict = code1LetterToCcpCodeDict[molType]
      #else:
      #  dict = {}
      
      i = 0;
      while i < len(text):
        letter = text[i].upper()
        if letter == '[':
          code = ''
          i +=1
          while (text[i] != ']') and (i < len(text)):
            code += text[i]
            i +=1
          seq.append(code)
          i +=1
      
        elif letter in ALPHABET:
          code = dict.get(letter, letter)
          seq.append(code)
          i += 1

        else:
          i += 1
 
    elif self.seqCodeMode == ccp:
      if self.polymerMolType != carbo:
        text = re.sub('!.+\n', ' ', text)
        text = re.sub('\W\d+\W', ' ', text)
        text = re.sub('[,.]', ' ', text)
        text = re.sub('residue', ' ', text)
      text = re.sub('\s+', ' ', text)
      text = re.sub('^\s+', '', text)
      text = text.rstrip()
      if text == '':
        seq = []
      else:
        seq = text.split()  
        
    isValid, self.inputSeq = self.checkSeq(seq, self.polymerMolType, checkChemComps=checkChemComps)
    
    return isValid


  def clearPolymerSeq(self):
  
    self.inputSeq = []
    self.seqTextBox.setText('')
  
  def tidyPolymerSeq(self, *opt):
  
    self.setPolymerStart()
    self.getInputPolymerSeq()
    self.displayPolymerSeq()
  
  def displayPolymerSeq(self):
  
    one, ccp = RES_CODE_MODES
    protein, dna, rna, carbo, dnaRna = POLYMER_MOLTYPES
  
    offset = self.seqStartNum-1
    seq = ''
    
    if self.seqCodeMode == one:
      i = self.seqStartNum
      typ =  self.polymerMolType
      for j in range( (i-1) % 50 ):
        seq = seq + ' '
        if j%10 == 0:
          seq = seq + ' '
      
      if typ == carbo:  
        inverseDict = {}
      else:
        inverseDict = ccpCodeToCode1LetterDict[typ]
        # TBD: is it possible for typ == 'DNA/RNA' here?
      
      for j in range( len(self.inputSeq)):
        seq = seq + inverseDict.get(self.inputSeq[j], '[%s]' % self.inputSeq[j])
        
        if i%10 == 0:
          seq = seq + ' '
        if i%50 == 0:
          seq = seq + ' %4.1d\n' % i
          
        i += 1
     
    elif self.seqCodeMode == ccp:
      i = self.seqStartNum
      
      for j in range( (i-1) % 10 ):
        seq = seq + '    '
        
      for j in range( len(self.inputSeq)):
        seq = seq + '%-3.32s' % self.inputSeq[j] + ' '
        if i%10 == 0:
          seq = seq + ' %4.1d\n' % i
          
        i += 1

    self.seqTextBox.setText(seq)
  
  def commitPolymerSeq(self):
  
    self.setPolymerStart()
    
    isValid  = self.getInputPolymerSeq(checkChemComps=True)
    isCyclic = self.cyclicSelect.get()
    molType  = self.polymerMolType
    sequence = self.inputSeq
    startNum = self.seqStartNum

    if not sequence:
      return
    
    if not isValid:
      showWarning('Failure','Sequence not interpretable',parent=self)
      return
      
    molSystem = self.molSystemS
    molecule = self.moleculeS
    
    if not molecule:
      molecule = self.addMolecule()
    
    if not molecule: # Cancelled
      return

    if molecule.isFinalised:
      msg = 'Cannot add residues to a locked molecule template.'
      showWarning('Failure', msg, parent=self)
      return
    
    addMolResidues(molecule, molType, sequence,
                   startNum=startNum, isCyclic=isCyclic)

    if molSystem:
      if molSystem is True:
        molSystem = self.addMolSys()
    
      chain = self.makeChain(molSystem=molSystem, molecule=molecule)
      molecule.isFinalised = True
      self.updateChains()
      self.chainsMatrix.selectObject(chain)
      self.tabbedFrame.select(0)
                      
    else:
      self.molecule = molecule
      self.updateMoleculePulldown()
      self.tabbedFrame.select(2)

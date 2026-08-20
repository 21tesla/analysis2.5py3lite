
"""
======================COPYRIGHT/LICENSE START==========================

EditMolLabelling.py: Part of the CcpNmr Analysis program

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

# TTD
#
# Is there any way to specify correlactions?
# - atomLabels belong to a resLabel fraction?
#
# Adjust Analysis to use molLabelling
#  - Browse atoms
#  - CalcDistConstr
#  - edit assign
#  - edit peak list
# AssignmentBasic.py
# ConstraintBasic.py
# MoleculeBasic.py
# NetworkAnchoring.py
# PeakBasic.py

# Check we have equiv of scheme incorporation getters
#  inclusing for resonances

# Pass in labelledMisture to core function; work out
#  based of peak/spec before hand
#
# or
#
# Boolean useMolLabelling
#  and derive mixture from peakList etc.

# getPeakLabelledMixture(peak, molecule)

# Cache lookups for speed? At atom level?
# Notifier on mixtures to reset

# Future
#
# Coloured molecule representation
#
# Colour structure atoms

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.LabelDivider import LabelDivider
from memops.gui.Label import Label
from memops.gui.MessageReporter import showYesNo, showOkCancel, showWarning
from memops.gui.MultiWidget import MultiWidget
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.core.MoleculeBasic import greekSortAtomNames
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.popups.IsotopeSchemeEditor import getSortedIsotopes

UNICODE_SUPERSCRIPT = {'0':u'\u2070','1':u'\u00B9','2':u'\u00B2','3':u'\u00B3','4':u'\u2074',
                       '5':u'\u2075','6':u'\u2076','7':u'\u2077','8':u'\u2078','9':u'\u2079'}


ORD_0 = ord('A')-1

UNLABELLED_RGB = (0.82, 0.82, 0.82)
LABELLED_RGB = (0.5, 0.0, 0.0)

class EditMolLabellingPopup(BasePopup):
  """
  **Setup Molecule Isotope Labelling Patterns**
  
  This popup is used to specify the isotopic incorporation of molecules in
  different samples. This is typically done in the earlier stages of a project, 
  before the labelling information is required for assignment etc.
  
  Although the user can employ the standard isotope labelling schemes when 
  performing operations in Analysis, the setup of molecular isotope patterns
  allows a greater level of control and allows the isotopic labelling
  to be automatically set according to context. The standard labelling schemes,
  which are also used as references for molecule specific patterns, may be
  created and edited via the `Reference Isotope Schemes`_ popup.
  
  A given molecule can be represented by several different samples, each with
  distinct isotope patterns. Each sample is associated with the NMR experiments
  (and hence spectra) that were performed on that sample. When assigning a peak,
  if the peak's experiment has been linked to a sample with labelling
  information then the isotopic proportions hat should be considered are
  automatically known.

  **Labelled Molecules**

  The first tab allows the user to setup the samples (labelled versions) of 
  each molecule that are available. Each molecule version will appear as a row
  in the top table. The user should set the "Experiments" of each sample to
  specify which NMR experiments used that labelled version on the molecule.

  Each sample of a molecule may refer to only a single labelled version or a
  mixture of different labelled forms (of the same molecule), e.g. both fully
  labelled and unlabelled. These different forms are referred to as "Patterns"
  and can be used in any proportion, to create a sample with a mixed
  composition. 

  The Sample Composition table allows the user to add patterns (labelled
  versions of the molecule) to the currently selected mixture. Often  the user
  will add a single pattern based on a standard labelling scheme; this is the
  simplest setup. However, any number of patterns may be added and their
  proportions in a mixed sample may be adjusted by changing the "weight" of
  each. It is possible for a given pattern to be re-used in a number of
  different samples. For example a natural abundance version of the molecule may
  appear on its own in a sample, but as a mixture together with an enriched
  pattern in another sample.

  **Labelling Patterns**
  
  The second tab allows the user to specify any deviations of a labelling
  pattern (version of a molecule) from a standard labelling scheme. An initial
  pattern will contain residues that are associated with an initial scheme; the
  residues isotopomers are set to be the same as in this scheme. However, the
  user may adjust any of the residues in the pattern to use a different scheme.
  Here the term "isotopomer" is used to refer to a single isotope pattern of an
  individual residue, and there may be several such patterns present for a given
  type of residue, so that a polymer molecule will generally be a random mixture
  of the different individual forms.

  A residue within a pattern may deviate from a standard scheme by adjusting the
  proportions of the different isotopomers; this is set via the "Weighting" value
  after selecting the required isotopomer in the pulldown menu. Adjusting these
  weightings preserves the internal correlations between atoms within
  isotopomers; only the relative abundance of separate isotopomer forms is
  affected.

  The individual atomic incorporations for a residue may be adjusted if
  the "Isotopomer" pulldown is set to "Override Scheme", but this looses any
  atom pair information that was specified via an isotopomer; only average
  incorporations are set. For example there may be no single isotopomer with
  both CA and CB labelled at the same time, even though the average shows
  incorporation for both sites.

  **Examples**
  
  These examples illustrate how to setup various kinds of isotope labelling
  "sample", starting from a project with molecule specifications but without
  existing labelling information. After creation each sample should be
  associated with the required NMR experiments (which were run on that sample),
  by double clicking in the "Experiments" column  of the Labelled Samples table
  and checking the required options.

  *Simple Single Scheme Labelling*
  
  As an example of specifying a molecule with a single pattern from the standard
  "1,3 13C Glycerol" scheme, the user would first select the required molecule
  (to the right of "For molecule") and press [New Sample].  In the lower right
  pulldown menu the user then selects "13Glycerol" and then presses [New Pattern
  From Scheme]. The molecule in this sample will have the labelling specified by
  one pattern made according to the selected scheme.

  *Mixture of Labelled and Unlabelled*
  
  As an example of a mixture, say that a sample of a molecule was composed of a
  1:3 ratio of fully labelled and unlabelled molecules (two different patterns
  in a mixture). The user would make one sample for the molecule using [New
  Sample], then setup the Sample Composition to have two patterns, one based on a
  uniform labelling scheme (e.g. "uni_15N13C2H") and one for the "NatAbun"
  scheme; select each scheme in the lower right pulldown menu and for each make
  [New Pattern From Scheme]. Labelling patterns "A" and "B" will be created. For
  the pattern that corresponds to the "NatAbun" basis scheme set the  weight to
  3 (the weight for the other pattern will remain at 1.0). This labelled sample
  specification is nor ready to used because the individual residues do not need
  any further adjustment.

  *Labelling a Specific Residue*

  As an example of a molecule that is unlabelled except at one specific sequence
  position, where the residue has 15N, the user would first make a New Sample
  for the required molecule and add a New Pattern from the "NatAbun" scheme;
  this will make Pattern "A". 

  The user should then go to the second "Labelling Patterns" tab (at the top) and
  ensure that the correct molecule and Pattern "A" are selected in the upper
  pulldown menus. The sequence of the molecule will be displayed in the upper
  table. To set the labelling for the only labelled residue, the user scrolls
  down until the residue is visible, then selects the row for the residue in the
  table (highlighting the row). Next, in the pulldown menu to the bottom right
  of the residue table, the user selects the scheme "uni_15N" and then presses
  the adjacent [Set Residue Composition From Scheme:]. The residue will now be listed in
  the table as using the "uni_15N", and not "NatAbun" like all the others.
  
  .. _`Reference Isotope Schemes`: IsotopeSchemeEditor.html

  """

  def __init__(self, parent, *args, **kw):
    
    self.molecule = None
    self.mixture = None
    self.molFraction = None
    self.molLabel = None # pattern
    self.resLabel = None
    self.resLabelFraction = False # deliberate; None meaningful
    self.atomLabelTuple = None
    self.chemElement = None
    self.isotopomer = None
       
    self.waitingMolLabel = False
    self.waitingAtoms = False
    
    BasePopup.__init__(self, parent, title="Molecule : Isotope Labelling", **kw)

  def body(self, guiFrame):

    self.geometry('750x500')

    guiFrame.expandGrid(0,0)
    
    tipTexts = ['Setup isotope labelling at the molecule level; the versions of a present molecule in different samples',
                'Setup the specific residue and atom labelling patterns in a given molecule version']
                
    options = ['Labelled Molecules',
               'Labelling Patterns']
      
    tabbedFrame = TabbedFrame(guiFrame, options=options,
                              tipTexts=tipTexts, grid=(0, 0),
                              callback=self.toggleTab)

    self.bottomButtons = UtilityButtonList(tabbedFrame.sideFrame,
                                           helpUrl=self.help_url,
                                           grid=(0,0), sticky='e')
    
    self.tabbedFrame = tabbedFrame
    frameA, frameB = tabbedFrame.frames

    #
    # Labelled Molecules
    #

    frameA.expandGrid(1,0)
    div = LabelDivider(frameA, text='Labelled Samples', grid=(0,0))
    
    self.mixtureNameEntry = Entry(self, text='', width=10,
                                  returnCallback=self.setMixtureName)

    self.experimentSelect = MultiWidget(self, CheckButton, minRows=0,
                                       callback=self.setExperiments,
                                       useImages=False)

    tipTexts = ['The name of the molecule (sequence template) that the isotope labelled sample refers to',
                'The serial number of the sample (mixture of labelling patterns) in the CCPN project',
                'A user editable name for the isotope labelled sample',
                'The NMR experiments that are asscociated with the isotope labelling sample; can be set at any time by the user',
                'One or more labelling patterns (versions of a molecyle), in diifferent proportions, constiture the sample',
                'The molecular system chains (assignable entities in a complex) that are associated with the labelled molecule',
                'Part of the residue sequence of the labelled molecule',
                'The estimated average molecular mass of the molecule (albeit a mixture of different patterns)']
                
    headingList = ['Molecule', 'Sample', 'Sample Name',
                   'Experiments', 'Labelling Patterns',
                   'Chains', 'Sequence', 'Mol Mass']

    editWidgets      = [None, None, self.mixtureNameEntry,
                        self.experimentSelect,
                        None, None, None, None]
    editGetCallbacks = [None, None, self.getMixtureName,
                        self.getExperiments,
                        None, None, None, None]
    editSetCallbacks = [None, None, self.setMixtureName,
                        self.setExperiments,
                        None, None, None, None]
    
    self.mixtureMatrix = ScrolledMatrix(frameA, tipTexts=tipTexts, 
                                        multiSelect=False, grid=(1, 0), 
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks, 
                                        editWidgets=editWidgets,
                                        headingList=headingList, 
                                        callback=self.selectMixture,
                                        deleteFunc=self.deleteMixture)
    
    frame = Frame(frameA, grid=(2,0))
    frame.expandGrid(0,0)
    
    tipTexts = ['Delete the selected isotope labelling molecule specification',
                'Make a new, intially bank, isotope labelled sample, to which molecule labelling patterns may be added']
    texts = ['Delete Sample','New Sample']
    commands = [self.deleteMixture,self.addMixture]
    buttonList = ButtonList(frame, tipTexts=tipTexts, texts=texts,
                            commands=commands, grid=(0,0))
    self.sampleButtons = buttonList.buttons

    tipText = 'Selects which molecular sequence template to use when making a new specifcation of an isotope labelled sample'
    label = Label(frame, text='For molecule:', grid=(0,1))
    self.mixMoleculePulldown = PulldownList(frame, grid=(0,2), tipText=tipText)

    div = LabelDivider(frameA, text='Sample Composition', grid=(3,0))

    self.molFracWeightEntry = FloatEntry(self, text=1.0, width=8,
                                         returnCallback=self.setMolFracWeight)

    tipTexts = ['A code letter for a specific labelling pattern version of the molecule sample',
                'The standard isotope labelling scheme on which the labelling pattern is primarily based (for mosed residue)',
                'The relative weighting of this labelling patterns compated to other patterns in the sample sample; for a mixture of molecule versions',
                'The relative composition of this pattern compated to others in the sample',
                'The estimated molecular masss of the version of the isotopic molecule represented by the pattern']
                
    headingList = ['Labelling Pattern', 'Basis Scheme',
                   'Weight', '% Composition', 'Mol Mass']

    editWidgets      = [None, None, self.molFracWeightEntry, None, None]
    editGetCallbacks = [None, None, self.getMolFracWeight, None, None]
    editSetCallbacks = [None, None, self.setMolFracWeight, None, None]
    
    self.molFractionMatrix = ScrolledMatrix(frameA, tipTexts=tipTexts, 
                                            multiSelect=False, grid=(4, 0),
                                            editSetCallbacks=editSetCallbacks,
                                            editGetCallbacks=editGetCallbacks,
                                            editWidgets=editWidgets,
                                            headingList=headingList,
                                            callback=self.selectMolFraction,
                                            deleteFunc=self.deleteMolFraction)
    
    frame = Frame(frameA, grid=(5,0))
    frame.expandGrid(0,0)

    tipTexts = ['Remove the selected molecule labelling pattern from the sample/mixture; the pattern is deleted if it is no longer used in any sample',
                'Edit the selected molecule labelling pattern (isotopic version), in terms of residue and atom incorporations',
                'Add a molecule labelling pattern to the current sample specification, using a pre-constructed pattern']
    texts = ['Remove Pattern','Edit Pattern','Add Existing Pattern:']
    commands = [self.deleteMolFraction,
                self.editMolLabel,
                self.addMolFraction]
    buttonList = ButtonList(frame, tipTexts=tipTexts, texts=texts,
                            commands=commands, grid=(0,0))
    self.patternButtons = buttonList.buttons

    tipText = 'Selects which extsing molecule labelling pattern to add to the selecte sample'
    self.addMolFracPulldown = PulldownList(frame, grid=(0,1), tipText=tipText,
                                           callback=self.updateButtons)

    tipTexts = ['make a complete new isotope labelling pattern for the molecule using the selected isitope labelling scheme',]
    texts = ['New Pattern From Scheme:']
    commands = [self.newMolFraction,]
    buttonList = ButtonList(frame, tipTexts=tipTexts, texts=texts,
                            commands=commands, grid=(0,2))
    self.patternButtons += buttonList.buttons
    
    tipText = 'Selects which standard isotope labelling scheme to use as a template to make a new molecule specific labelling pattern'
    self.newMolFracSchemePulldown = PulldownList(frame, grid=(0,3), tipText=tipText)

    #
    # Labelling Patterns
    #
    
    frameB.expandGrid(2,0)
    frameB.expandGrid(6,0)

    frame = Frame(frameB, grid=(0,0))
    frame.expandGrid(0,4)
        
    tipText = 'Selects which molecule to view the residue labelling patterns for'
    label = Label(frame, text='Molecule:', grid=(0,0))
    self.moleculePulldown = PulldownList(frame, grid=(0,1), tipText=tipText,
                                         callback=self.changePatternMolecule)
    
    tipText = 'Selects which residue labelling to view from those available for the selected molecule'
    label = Label(frame, text='Pattern:', grid=(0,2))
    self.patternPulldown = PulldownList(frame, grid=(0,3), tipText=tipText,
                                        callback=self.changeMolLabel)
    
    tipTexts = ['make a new labelling pattern (isotope version of the molecule) using the selected standard labelling scheme',]
    texts = ['New Pattern From Scheme:']
    commands = [self.newPattern,]
    buttonList = ButtonList(frame, tipTexts=tipTexts, texts=texts,
                            commands=commands, grid=(0,5))

    tipText = 'Selects which standardised isotope labelling scheme to use as a template for creating a to make a new molecule specific labelling pattern'
    self.newMolPatternSchemePulldown = PulldownList(frame, grid=(0,6), tipText=tipText)

    div = LabelDivider(frameB, text='Residue Pattern', grid=(1,0))
    
    tipTexts = ['The sequence ID number of the labelled residue in the molecule; a separate numbering system to the numbering in chains',
                'The code of the residue identifying which type it is',
                'The names of the standard isotope labelling schemes on which the residues labelling is based within this specific pattern',
                'Tne numbers of the isotopomer versions of this residue that specify labelling; from the standard scheme']
                
    headingList = ['Seq Id','Residue Code','Schemes','Isotopomers']

    editWidgets      = [None, None, None, None]
    editGetCallbacks = [None, None, None, None]
    editSetCallbacks = [None, None, None, None]
    
    self.residueMatrix = ScrolledMatrix(frameB, tipTexts=tipTexts, 
                                        multiSelect=True, grid=(2, 0),
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        headingList=headingList,
                                        callback=self.selectResLabel)

    frame = Frame(frameB, grid=(3,0))
    frame.expandGrid(0,0)
    
    tipTexts = ['Set the selected residue labelling to have natural abundance isotopic incorporation labels.',
                'Set the selected residue to have labelling based upon isotopomers from the selected standard scheme']
    texts = ['Set Natural Abundance','Set Residue Composition From Scheme:']
    commands = [self.setResidueNatAbun,
                self.setResidueIsotopomers]
    buttonList = ButtonList(frame, tipTexts=tipTexts, texts=texts,
                            commands=commands, grid=(0,0))
    self.residueButtons = buttonList.buttons

    tipText = 'Selects which standardised isotope labelling scheme to base the residues isotopomer labelling upon'
    self.residueSchemePulldown = PulldownList(frame, grid=(0,1), tipText=tipText)
    
    self.residueDiv = LabelDivider(frameB, text='Residue Composition', grid=(4,0))
    
    frame = Frame(frameB, grid=(5,0))
    frame.expandGrid(0,4)
    
    tipText = 'Selects which chemical element to display atomic isotope incorporation levels for'
    label = Label(frame, text='Chem Element:', grid=(0,0))
    self.chemElementPulldown = PulldownList(frame, grid=(0,1), tipText=tipText,
                                            callback=self.changeChemElement,
                                            texts=['C','N','H'], index=0)
    
    tipText = 'Sets whether to display water exchangeable atoms in the atom isotope table'
    label = Label(frame, text='Show Exchangeable', grid=(0,2))
    self.waterExchangeSelect = CheckButton(frame, tipText=tipText,
                                           selected=False, grid=(0,3),
                                           callback=self.updateAtomMatrix)
    
    tipText = 'Selects which isotopomer version, from those set for the residue to show atom istopes for; can be set to override the schem istopomers'
    label = Label(frame, text='  Isotopomer:', grid=(0,5))
    self.isotopomerPulldown = PulldownList(frame, grid=(0,6), tipText=tipText,
                                           callback=self.changeIsotopomer)
    
    tipText = 'Specifies what relative weighing the selected residue isotopomer has compared to others for the same residue'
    self.isotopomerFractionLabel = Label(frame, text=' Weighting:', grid=(0,7))
    self.isotopomerFractionEntry = FloatEntry(frame, text=1.0, width=8, tipText=tipText,
                                              returnCallback=self.setIsotopomerWeight,
                                              grid=(0,8))
    self.isotopomerFractionEntry.bind('<Leave>', self.setIsotopomerWeight, '+')                                   
    
    
    self.atomLabelWeightEntry = FloatEntry(self, text=1.0, width=8,
                                           returnCallback=self.setAtomLabelWeight)

    headingList, tipTexts, editWidgets, editGetCallbacks, editSetCallbacks = self.getIsotopeHeadingLists()
    
    self.atomMatrix = ScrolledMatrix(frameB, tipTexts=tipTexts, 
                                     multiSelect=True, grid=(6, 0),
                                     highlightType='cell',
                                     editSetCallbacks=editSetCallbacks,
                                     editGetCallbacks=editGetCallbacks,
                                     editWidgets=editWidgets,
                                     headingList=headingList,
                                     callback=self.selectAtomLabel)

    tipTexts = ['Resets the selected residue to use atom isotope incorporation levels from the isotopomers of the set scheme; removes any atom level edits',
                'Spread the istotope incorporation levels from the last selected atom (in the table) to all selected atomes; overrides isotopomer levels',
                'Set the selected atoms to have natural abundance isotope incorporation levels']
    texts = ['Reset Residue To Orig Scheme', 'Propagate Abundances', 'Set Natural Abundance']
    commands = [self.resetResidue,
                self.propagateAbundances,
                self.setAtomsNatAbun]
    buttonList = ButtonList(frameB, tipTexts=tipTexts, texts=texts,
                            commands=commands, grid=(7,0))
    self.atomButtons = buttonList.buttons
    
    #
    
    self.updateMixtureMatrix()
    self.updateMoleculePulldown()
    self.updateMixMoleculePulldown()
    self.updateAddMolFracPulldown()
    self.updateNewMolFracSchemePulldown()
    self.updatePatternPulldown()
    self.updateResidueSchemePulldown()
    self.updateNewMolPatternSchemePulldown()
    self.updateChemElementPulldown()
    self.updateIsotopomerPulldown()
    self.updateButtons()
    self.administerNotifiers(self.registerNotify)
  
  def toggleTab(self, index):
  
    if index == 0:
      self.updateMolFractionMatrix()
      self.updateMixtureMatrix()
  
  def updateSchemes(self, labellingScheme):
  
    self.updateNewMolFracSchemePulldown()
    self.updateResidueSchemePulldown()
    self.updateNewMolPatternSchemePulldown()

  def selectMixture(self, mixture, row, col):
  
    if mixture is not self.mixture:
      self.mixture = mixture
      self.updateMolFractionMatrix()
      self.updateAddMolFracPulldown()
      self.updateButtons()

  def setMixtureName(self, event):
  
    name = self.mixtureNameEntry.get().strip() or None
    
    if name != self.mixture.name:
      self.mixture.name = name

  def getMixtureName(self, labelledMixture):
  
    self.mixtureNameEntry.set(labelledMixture.name)
 
  def getPossMixExperiments(self, mixture):
    
    experiments = self.nmrProject.sortedExperiments()
    molecule = mixture.labeledMolecule.molecule
    
    exclude = set()
    for labelledMolecule in self.project.labeledMolecules:
      if labelledMolecule.molecule is molecule:
        for labelledMixture in labelledMolecule.labeledMixtures:
          if labelledMixture is mixture:
            continue
        
          exclude.update(labelledMixture.experiments)
  
    return [e for e in experiments if e not in exclude]
  
  def getExperiments(self, mixture):
  
    experiments = self.getPossMixExperiments(mixture)
    
    names  = []
    values = []
    for experiment in experiments:
      names.append(experiment.name)
      
      if experiment in mixture.experiments:
        values.append(True)
      else:
        values.append(False)  
  
    self.experimentSelect.set(values=values, options=names)
  
  def setExperiments(self, obj):

    if self.mixture:
      if obj is None:
        self.mixtureMatrix.keyPressEscape()
      else:
        experiments = self.getPossMixExperiments(self.mixture)
        values = self.experimentSelect.get()
        selected = [experiments[i] for i in range(len(values)) if values[i]]
        self.mixture.setExperiments(selected)
        self.mixtureMatrix.keyPressEscape()

  def deleteMixture(self, *event):
  
    msg = 'Really delete selected sample?'
    if self.mixture and showWarning('Confirm', msg, parent=self):
      self.mixture.delete()
      self.mixture = None

  def getLabelledMolecule(self, molecule):
  
    if molecule:
      labelledMolecule = self.project.findFirstLabeledMolecule(molecule=molecule)
      
      if labelledMolecule is None:
        root = molecule.root
        labelledMolecule = root.newLabeledMolecule(name=molecule.name)

      return labelledMolecule

  def addMixture(self):
    
    molecule = self.mixMoleculePulldown.getObject()
    if molecule:
      labelledMolecule = self.getLabelledMolecule(molecule)
      labelledMolecule.newLabeledMixture()

  def setMolFracWeight(self, event):
  
    if self.molFraction:
      weight = abs(self.molFracWeightEntry.get() or 0.0)
      self.molFraction.weight = weight

  def getMolFracWeight(self, molFraction):
  
    self.molFracWeightEntry.set(molFraction.weight)

  def deleteMolFraction(self, *event):
  
    if self.molFraction:
      molLabel = self.molFraction.molLabel
    
      self.molFraction.delete()
      self.molFraction = None
      
      if not molLabel.molLabelFractions:
        if molLabel is self.molLabel:
          self.changeMolLabel(None)
        molLabel.delete()

  def selectMolFraction(self, molFraction, row, col):
  
    if molFraction is not self.molFraction:
      self.molFraction = molFraction
      self.updateButtons()
      

  def editMolLabel(self):
  
    if self.molFraction:
      self.changeMolLabel(self.molFraction.molLabel)
      self.updatePatternPulldown()
      self.tabbedFrame.select(1)

  def addMolFraction(self):
  
    molLabel = self.addMolFracPulldown.getObject()
    if self.mixture and molLabel:
      fraction = self.mixture.newMolLabelFraction(molLabel=molLabel, weight=1.0)

  def newMolFraction(self):
  
    scheme = self.newMolFracSchemePulldown.getObject()
    
    if not self.mixture:
      showWarning('Warning', 'No labelled sample selected', parent=self)
      return    
    
    if self.mixture and scheme:
      molLabel = getMolLabelFromScheme(self.mixture.labeledMolecule, scheme)
      fraction = self.mixture.newMolLabelFraction(molLabel=molLabel, weight=1.0)

  def changePatternMolecule(self, molecule):
  
    if molecule is not self.molecule:
      self.molecule = molecule
      labelledMolecule = self.getLabelledMolecule(molecule)
      molLabels = labelledMolecule.sortedMolLabels()
      
      if molLabels:
        self.molLabel = molLabels[0]
      else:
        self.molLabel = None
      
      self.resLabel = None
      self.updateResidueMatrix()
      self.updatePatternPulldown()
      self.updateIsotopomerPulldown()

  def changeMolLabel(self, molLabel):
    
    if molLabel is not self.molLabel:
      self.molLabel = molLabel
      self.updateResidueMatrix()
      self.resLabel = None
      self.updateIsotopomerPulldown()

  def newPattern(self):
  
    scheme = self.newMolPatternSchemePulldown.getObject()
     
    if self.molecule and scheme:
      labelledMolecule = self.getLabelledMolecule(self.molecule)
      molLabel = getMolLabelFromScheme(labelledMolecule, scheme)
      self.changeMolLabel(molLabel)

  def selectResLabel(self, resLabel, row, col):
  
    if self.resLabel is not resLabel:
      self.resLabel = resLabel
      self.updateIsotopomerPulldown()
      self.updateButtons()

  def setResidueNatAbun(self):
  
    if self.resLabel:
      resLabels = self.residueMatrix.currentObjects[:]
      
      for resLabel in resLabels:
        setResLabelNaturalAbundance(resLabel)

  def setResidueIsotopomers(self):
  
    if self.resLabel:
      resLabels = self.residueMatrix.currentObjects[:]
      scheme = self.residueSchemePulldown.getObject()

      for resLabel in resLabels:
        updateResLabelFractions(resLabel, scheme)

  def changeIsotopomer(self, resLabelFraction):
  
    if resLabelFraction is not self.resLabelFraction:
      self.resLabelFraction = resLabelFraction
      if resLabelFraction and self.resLabelFraction.isotopomerSerial:
        self.isotopomer = list(self.resLabelFraction.isotopomers)[0]
      else:
        self.isotopomer = None
        
      self.updateAtomMatrix()

  def changeChemElement(self, chemElement):
  
    if chemElement != self.chemElement:
      self.chemElement = chemElement
      self.updateAtomMatrix()

  def setIsotopomerWeight(self, event=None):
  
    if self.resLabelFraction:
      weight = self.isotopomerFractionEntry.get() or 0.0
      
      if weight != self.resLabelFraction.weight:
        self.resLabelFraction.weight = weight
        self.updateAtomMatrix()

    else:
      self.isotopomerFractionEntry.set(None)

  def getIsotopeHeadingLists(self):
  
    element  = self.chemElementPulldown.getText()
    
    headingList = ['Atom',]
    tipTexts = ['The name of the atom within its residue, for which to set isotopic abundances, within the selected labelling pattern',]
    
    isotopeCodes = []
    for code, isotope in  getSortedIsotopes(self.project, element):
      massNumber = isotope.massNumber
      superscript = u''
      
      for digit in str(massNumber):
        superscript += UNICODE_SUPERSCRIPT[digit]
    
      isotopeCodes.append(superscript+element)
    
    headingList.extend(['Weighting %s' % x for x in  isotopeCodes])
    tip = 'The amount of %s isotope incorporation; can be a ratio, percentage or fraction (the value used is relative to the sum of all weights)'
    tipTexts.extend([tip % x for x in isotopeCodes])
    
    tip = 'The percentage %s isotope incorporation for this atom, calculated using stated weights'
    headingList.extend(['%% %s' % x for x in  isotopeCodes])
    tipTexts.extend([tip % x for x in isotopeCodes])
    
    editWidgets      = [None,]
    editGetCallbacks = [None,]
    editSetCallbacks = [None,]
    editWidgets.extend([self.atomLabelWeightEntry for x in isotopeCodes])
    editGetCallbacks.extend([self.getAtomLabelWeight for x in isotopeCodes])
    editSetCallbacks.extend([self.setAtomLabelWeight for x in isotopeCodes])
    editWidgets.extend([None for x in isotopeCodes])
    editGetCallbacks.extend([None for x in isotopeCodes])
    editSetCallbacks.extend([None for x in isotopeCodes])

    return headingList, tipTexts, editWidgets, editGetCallbacks, editSetCallbacks

  def getAtomLabelWeight(self, obj):
    
    atomName, atomLabels, col = self.atomLabelTuple
  
    if atomLabels and (col > 0):
      atomLabel, isotope, weight = atomLabels[col-1]
      
      if atomLabel is None:
        weight = 0.0
            
      self.atomLabelWeightEntry.set(weight)
  
  def setAtomLabelWeight(self, event):

    newWeight = self.atomLabelWeightEntry.get()
   
    if newWeight is not None:
      atomName, atomLabels, col = self.atomLabelTuple
      
      for i, (atomLabel, isotope, weight) in enumerate(atomLabels):
        
        if i == col-1:
          weight = newWeight
      
        if (not atomLabel) or (atomLabel.className != 'SingleAtomLabel'):
          atomLabel = self.resLabel.newSingleAtomLabel(atomName=atomName,
                                                 massNumber=isotope.massNumber,
                                                 weight=weight)
        else:
          atomLabel.setWeight(weight)


  def selectAtomLabel(self, obj, row, col):
  
    atomName, atomLabels = obj
    self.atomLabelTuple = (atomName, atomLabels, col)
    self.updateButtons()

  def resetResidue(self):
  
    if self.resLabel:
      resLabelFraction = self.resLabel.findFirstResLabelFraction()
      scheme = resLabelFraction.labelingScheme
      updateResLabelFractions(self.resLabel, scheme)

  def propagateAbundances(self):
  
    resLabel = self.resLabel
    atomName, atomLabels, col = self.atomLabelTuple
    atomLabelTuples = self.atomMatrix.currentObjects[:]
  
    weightDict = {}
    for atomLabel, isotope, weight in atomLabels:
       
      if atomLabel: 
        weightDict[isotope.massNumber] = weight
    
    for atomName0, atomLabels0 in atomLabelTuples:
      if atomName0 != atomName:
      
        doneLabels = set()
        for massNumber in weightDict:
          weight = weightDict[massNumber]
          atomLabel = resLabel.findFirstAtomLabel(atomName=atomName0,
                                                  massNumber=massNumber)
          if atomLabel is None:
            atomLabel = resLabel.newSingleAtomLabel(atomName=atomName0, 
                                              massNumber=massNumber,
                                              weight=weight)
          else:
            atomLabel.weight = weight
            
          doneLabels.add(atomLabel)  
 
  def setAtomsNatAbun(self):
  
    atomLabelTuples = self.atomMatrix.currentObjects
     
    for atomName, atomLabels in atomLabelTuples:
      for atomLabel, isotope, weight in atomLabels:
        
        if (not atomLabel) or (atomLabel.className != 'SingleAtomLabel'):
          atomLabel = self.resLabel.newSingleAtomLabel(atomName=atomName,
                                                 massNumber=isotope.massNumber,
                                                 weight=isotope.abundance)
        else:
          atomLabel.setWeight(isotope.abundance)

  def updateMixMoleculePulldown(self):

    index = 0
    molecules = self.project.sortedMolecules()
    names = [m.name for m in molecules]
    
    if self.molecule in molecules:
      index = molecules.index(self.molecule)
    
    self.mixMoleculePulldown.setup(names, molecules, index)
 
  def updateAddMolFracPulldown(self):

    index = 0
    names = []
    molLabels = []

    if self.mixture:
      molLabels = self.mixture.labeledMolecule.sortedMolLabels()
      
      for molLabelFraction in self.mixture.molLabelFractions:
        molLabels.remove(molLabelFraction.molLabel)

      names = [chr(ml.serial+ORD_0) for ml in molLabels]

    self.addMolFracPulldown.setup(names, molLabels, index)
 
  def updateNewMolFracSchemePulldown(self):

    index = 0
    schemes = self.project.sortedLabelingSchemes()
    names = [s.name for s in schemes]

    # Could filter on what fractions are already present

    self.newMolFracSchemePulldown.setup(names, schemes, index)

  def updateMoleculePulldown(self):

    index = 0
    molecules = self.project.sortedMolecules()
    names = [m.name for m in molecules]
    molecule = self.molecule
    
    if molecules:
      if molecule in molecules:
        index = molecules.index(molecule)
      else:
        molecule = molecules[0]
    
    else:
      molecule = None

    self.changePatternMolecule(molecule)    
    self.moleculePulldown.setup(names, molecules, index)

  def updatePatternPulldown(self):

    index = 0
    molLabel = self.molLabel
    
    if self.molecule:
      labeledMolecule = self.getLabelledMolecule(self.molecule)
      molLabels = labeledMolecule.sortedMolLabels()
    else:
      molLabels = []
    
    if molLabels:
      names = [chr(ml.serial+ORD_0) for ml in molLabels]
    
      if molLabel in molLabels:
        index = molLabels.index(molLabel)
      else:
        molLabel = molLabels[0]
    
    else:
      names = []
      molLabel = None

    self.changeMolLabel(molLabel)
    self.patternPulldown.setup(names, molLabels, index)
    self.waitingMolLabel = False
    
  def updateNewMolPatternSchemePulldown(self):

    index = 0
    schemes = self.project.sortedLabelingSchemes()
    names = [s.name for s in schemes]

    # Could filter on what molLabels are already present

    self.newMolPatternSchemePulldown.setup(names, schemes, index)

  def updateResidueSchemePulldown(self):

    index = 0
    schemes = self.project.sortedLabelingSchemes()
    names = [s.name for s in schemes]

    # Could filter on what molLabels are already present

    self.residueSchemePulldown.setup(names, schemes, index)

  def updateIsotopomerPulldown(self):

    index = 0
    names = []
    resLabelFractions = []
    resLabelFraction = self.resLabelFraction
    
    if self.resLabel:
      resLabelFractions = self.resLabel.sortedResLabelFractions()
      resLabelFractions = [rlf for rlf in resLabelFractions if rlf.isotopomerSerial]
      
      for rlf in resLabelFractions:
        isotopomer = list(rlf.isotopomers)[0]
        chemCompLabel = isotopomer.chemCompLabel
        scheme = chemCompLabel.labelingScheme
        
        name = '%s:%s:%d' % (scheme.name, chemCompLabel.ccpCode, isotopomer.serial)
        names.append(name)
      
      resLabelFractions.append(None)
      names.append("*Override Scheme*")
      
      if resLabelFraction not in resLabelFractions:
        resLabelFraction = resLabelFractions[0]
    
      index = resLabelFractions.index(resLabelFraction)
      
    else:
      resLabelFraction = False
    
    
    self.changeIsotopomer(resLabelFraction)
    
    self.isotopomerPulldown.setup(names, resLabelFractions, index)

  def updateChemElementPulldown(self):

    chemElement = self.chemElement
    
    if self.isotopomer:    
      chemCompLabel = self.isotopomer.chemCompLabel
      elementDict = {}
      
      for chemAtom in chemCompLabel.chemComp.chemAtoms:
        symbol = chemAtom.elementSymbol
        if symbol:
          elementDict[symbol] = True
      
      names = elementDict.keys()
      names.sort()
      
    else:
      names = ['C','N','H']
      
    if chemElement not in names:
      chemElement  = 'C'
    
    index = names.index(chemElement)

    self.changeChemElement(chemElement)
    self.chemElementPulldown.setup(names, names, index)
  
  def updateAtomLabelAfter(self, atomLabel):
    
    if self.waitingAtoms:
      return
    
    resLabel = atomLabel.resLabel
    if resLabel is self.resLabel:
      self.waitingAtoms = True
      self.after_idle(self.updateAtomMatrix)

  def updateAtomMatrix(self, event=None):
        
    element  = self.chemElementPulldown.getText()
    doExchangeable = self.waterExchangeSelect.get()
    isotopes = [x[1] for x in getSortedIsotopes(self.project, element)]
    massNumbers = set([i.massNumber for i in  isotopes])
    isotopeDict = {}
    for isotope in isotopes:
      isotopeDict[isotope.massNumber] = isotope
    
    isotopeCodes = set(['%d%s' % (i.massNumber, element) for i in  isotopes])
    resLabel = self.resLabel

    textMatrix = []
    objectList = []
    colorMatrix = []
    headingList, tipTexts, editWidgets, editGetCallbacks, editSetCallbacks = self.getIsotopeHeadingLists()
    residueText = 'Residue Composition'
    
    if resLabel:
      
      # Get chem atoms
      
      molResidue = resLabel.molResidue
      chemAtoms = molResidue.chemCompVar.chemAtoms
      chemAtoms = [a for a in chemAtoms if a.className != 'LinkAtom']
      chemAtoms = [a for a in chemAtoms if a.elementSymbol == element]
      if not doExchangeable:
        chemAtoms = [a for a in chemAtoms if not a.waterExchangeable]
      
      data = {}
      atomNames = []
      residueText += ': %s %s' % (molResidue.seqCode, molResidue.ccpCode)
      
      if self.resLabelFraction and self.isotopomer:
        # # # # Specific isotopomer # # # #
        self.isotopomerFractionEntry.set(self.resLabelFraction.weight)
        self.isotopomerFractionLabel.grid(row=0, column=7)
        self.isotopomerFractionEntry.grid(row=0, column=8)
        
        isotopomer = self.isotopomer
        
        for chemAtom in chemAtoms:
          atomName = chemAtom.name
          atomNames.append(atomName)
          isotopeData = {}
         
          if not (resLabel.findFirstAtomLabel(atomName=atomName) or \
             resLabel.findFirstAtomLabel(elementName=element)):
            # Not Overrriden
         
            atomLabels = isotopomer.findAllAtomLabels(name=atomName)
            atomLabels = [al for al in atomLabels if al.isotopeCode in isotopeCodes]
            sumFrac = sum([al.weight for al in atomLabels])
 
            for al in atomLabels:
              massNumber = al.isotope.massNumber
              isotopeData[massNumber] = (al.weight/sumFrac, al)
          
          else:
            isotopeData = None
          
          data[atomName] = isotopeData, 1.0
          
      else:
        # # # # Overall/average # # # #
       
        self.isotopomerFractionEntry.set(None)
        self.isotopomerFractionLabel.grid_forget()
        self.isotopomerFractionEntry.grid_forget()
     
        # Get scheme proportions
 
        isotopomerFracs = []
        for resLabelFraction in resLabel.resLabelFractions:
          allIsotopomers = resLabelFraction.isotopomerSerial == 0
          frac = resLabelFraction.weight
 
          if allIsotopomers:
            isotopomerFracsA = [(i.weight, i) for i in resLabel.isotopomers]
            sumFracA = sum([x[0] for x in isotopomerFracsA])
            isotopomerFracs.extend([(frac*i.weight/sumFracA, i) for i in isotopomerFracsA])
 
          else:
            for isotopomer in resLabelFraction.isotopomers:
              isotopomerFracs.append((frac, isotopomer))
 
        sumFrac = sum([x[0] for x in isotopomerFracs])
 
        # Get atoms & weights
 
        for chemAtom in chemAtoms:
          atomName = chemAtom.name
          atomNames.append(atomName)
          isotopeData = {}
 
          # Try specific label
 
          atomLabels = resLabel.findAllAtomLabels(atomName=atomName)
          atomLabels = [al for al in atomLabels if al.massNumber in massNumbers]
 
          # Try uniform label
 
          if not atomLabels:
            atomLabels = resLabel.findAllAtomLabels(elementName=element)
            atomLabels = [al for al in atomLabels if al.massNumber in massNumbers]
 
          totalWeight = 0.0
          for al in atomLabels:
            isotopeData[al.massNumber] = (al.weight, al)
            totalWeight += al.weight
 
          # Try isotopomers from resLabelFractions
 
          if not isotopeData:
 
            for frac, isotopomer in isotopomerFracs:
              atomLabelsB = isotopomer.findAllAtomLabels(name=atomName)
              atomLabelsB = [al for al in atomLabelsB if al.isotopeCode in isotopeCodes]
              fracB = frac/sumFrac
              sumFracB = sum([al.weight for al in atomLabelsB])
 
              for al in atomLabelsB:
                massNumber = al.isotope.massNumber
 
                if massNumber in isotopeData:
                  weightB = isotopeData[massNumber][0]
                  isotopeData[massNumber] = (fracB*al.weight/sumFracB + weightB, al)
                else:
                  isotopeData[massNumber] = (fracB*al.weight/sumFracB, al)
 
            totalWeight = 1.0
 
          # Try natural abundance

          if not isotopeData:
            chemElement = chemAtom.chemElement
 
            for isotope in chemElement.isotopes:
 
             if isotope.abundance < 0.0003:
               if isotope.spin == '0': # Look into this again!
                 continue
 
             isotopeData[isotope.massNumber] = (isotope.abundance, None)
 
             totalWeight = 100.0
 
          #
 
          data[atomName] = isotopeData, totalWeight


      atomNames = greekSortAtomNames(data.keys())
      massNumbers = list(massNumbers)
      massNumbers.sort()
      
      for name in atomNames:
        isotopeData, totalWeight = data[name]
        row = [name,]
        atomLabels = []
       
        if isotopeData is not None:
          spinActiveWeight = 0.0
 
          for massNumber in massNumbers:
            weight, atomLabel = isotopeData.get(massNumber, (0.0, None))
            isotope = isotopeDict[massNumber]
            row.append(weight)
            atomLabels.append((atomLabel, isotope, weight))
 
            if isotope.spin == '1/2':
              spinActiveWeight += weight

          for massNumber in massNumbers:
            weight, atomLabel = isotopeData.get(massNumber, (0.0, None))
            row.append(100.0*weight/totalWeight)

          fraction = spinActiveWeight/totalWeight
          r1, g1, b1 = UNLABELLED_RGB
          r2, g2, b2 = LABELLED_RGB
 
          fr = 255 * fraction
          un = 255 - fr
          rgb = ((r1*un)+(fr*r2),
                 (g1*un)+(fr*g2),
                 (b1*un)+(fr*b2))
 
          colors = [None] * len(row)
          colors[0] = '#%02x%02x%02x' % rgb
        
        else:
          row += [None] * len(massNumbers)
          row += ['*Edited*'] * len(massNumbers)
          colors = ['#808080'] * len(row)
            

        textMatrix.append(row)
        objectList.append((name, atomLabels))
        colorMatrix.append(colors)

    self.residueDiv.setText(residueText)

    self.atomMatrix.update(textMatrix=textMatrix,
                           objectList=objectList,
                           colorMatrix=colorMatrix,
                           headingList=headingList,
                           tipTexts=tipTexts,
                           editWidgets=editWidgets,
                           editGetCallbacks=editGetCallbacks,
                           editSetCallbacks=editSetCallbacks)
    self.atomMatrix.doEditMarkExtraRules = self.doEditMarkExtraRules
    
    self.waitingAtoms = False
  
  def doEditMarkExtraRules(self, obj, row, col):
  
    if self.resLabelFraction:
      return False
      
    else:
      return True  

  def updateMixtureMatrix(self, mixture=None):

    textMatrix = []
    objectList = []
 
    for labelledMolecule in self.project.sortedLabeledMolecules():
      molecule = labelledMolecule.molecule
      chainText = ' '.join([c.code for c in molecule.chains])
      sequence = molecule.seqString
      if len(sequence) > 10:
        sequence = sequence[:10] + '...'
    
      for mixture in labelledMolecule.sortedLabeledMixtures():
        
        expNames = []
        for i, experiment in enumerate(mixture.experiments):
          expNames.append(experiment.name)
          if i and (i % 5 == 0):
            expNames.append('\n')
        
        molMass = 0.0
        molLabelTexts = []
        totalWeight = sum([mlf.weight for mlf in mixture.molLabelFractions]) or 1.0
        for molLabelFraction in mixture.molLabelFractions:
          molLabel  = molLabelFraction.molLabel
          weight = molLabelFraction.weight
          pc = 100.0 * weight / totalWeight
          letter = chr(ORD_0+molLabel.serial)
          molLabelTexts.append('%s(%.1f%%)' % (letter, pc))
          molMass += weight * getMolLabelMass(molLabel)
        
        molMass /= totalWeight
        molLabelTexts.sort()
         
        datum = [labelledMolecule.name,
                 mixture.serial,
                 mixture.name,
                 ' '.join(expNames).rstrip(),
                 ' '.join(molLabelTexts),
                 chainText,
                 sequence,
                 molMass]
        
        objectList.append(mixture)
        textMatrix.append(datum)
 
 
    self.mixtureMatrix.update(textMatrix=textMatrix,
                              objectList=objectList)
 
  def updateMolFractionMatrix(self):

    self.updateAddMolFracPulldown()
    
    textMatrix = []
    objectList = []

    mixture = self.mixture
    if mixture:
      totalweight = sum([mlf.weight for mlf in mixture.molLabelFractions]) or 1.0
    
      sortList = []
      for molLabelFraction in mixture.sortedMolLabelFractions():
        molLabel = molLabelFraction.molLabel
        pattern = chr(ORD_0+molLabel.serial)
        
        sortList.append((pattern, molLabelFraction.weight, molLabel, molLabelFraction))
      
      sortList.sort()
      
      for  pattern, serial, molLabel, molLabelFraction in sortList:
        schemes = set([])
        for resLabel in molLabel.resLabels:
          for resLabelFraction in resLabel.resLabelFractions:
            schemes.add(resLabelFraction.schemeName)
    
        schemes = list(schemes)
        schemes.sort()
        
        weight = molLabelFraction.weight
        comp = 100.0 * weight / totalweight
        
        molMass = getMolLabelMass(molLabel)
        
        datum = [pattern,
                 ','.join(schemes),
                 weight,
                 comp,
                 molMass]
        
        objectList.append(molLabelFraction)
        textMatrix.append(datum)
    
    self.molFractionMatrix.update(textMatrix=textMatrix,
                                  objectList=objectList)
  
  
  def updateResLabelFractionAfter(self, resLabelFraction):
    
     
    resLabel = resLabelFraction.resLabel
    if resLabel.molLabel is self.molLabel:
      self.after_idle(self.updateResidueMatrix)
          
    if resLabel is self.resLabel:
      self.after_idle(self.updateIsotopomerPulldown)
  

  def updateResidueMatrix(self):

    textMatrix = []
    objectList = []

    if self.molLabel:
    
      for resLabel in self.molLabel.sortedResLabels():
        molResidue = resLabel.molResidue
        
        schemes = set([])
        isotopomers = {}
        for resLabelFraction in resLabel.resLabelFractions:
          scheme = resLabelFraction.schemeName
          schemes.add(scheme)
          
          if scheme not in isotopomers:
            isotopomers[scheme] = []
            
          isotopomers[scheme].append(resLabelFraction.isotopomerSerial)
    
        schemes = list(schemes)
        schemes.sort()
        
        if len(schemes) == 1:
          isotopomers = isotopomers[schemes[0]]
          isotopomers = [str(x) or '<All>' for x in isotopomers]
          isotopomers.sort()
          isotopomerText = ','.join(isotopomers)
 
        else:
          isotopomerText = ''
          for scheme in schemes:
            isotopomers = isotopomers[scheme]
            isotopomers = [str(x) or '<All>' for x in isotopomers]
            isotopomers.sort()
 
            isotopomerText += '%s(%s) ' % (scheme, ','.join(isotopomers))
          
          isotopomerText = isotopomerText.rstrip()
          
        datum = [molResidue.seqCode,
                 molResidue.ccpCode,
                 ' '.join(schemes),
                 isotopomerText]
        
        objectList.append(resLabel)
        textMatrix.append(datum)
        
    
    self.residueMatrix.update(textMatrix=textMatrix,
                              objectList=objectList)
  
  def updateResLabelAfter(self, resLabel):
    
    if self.waitingMolLabel:
      return

    self.updateMolLabelAfter(resLabel.molLabel)
  
  def updateMolLabelAfter(self, molLabel):
    
    if self.waitingMolLabel:
      return
    
    else:
      self.waitingMolLabel = True
      self.after_idle(self.updatePatternPulldown)
      self.after_idle(self.updateMolFractionMatrix)
  
  def updateMixtureAfter(self, mixture):
  
    self.after_idle(self.updateMixtureMatrix)

    if mixture is self.mixture:
      self.after_idle(self.updateMolFractionMatrix)

  def updateMolLabelFractionAfter(self, molLabelFraction):
  
    self.after_idle(self.updateMixtureMatrix)
    
    if molLabelFraction.labeledMixture is self.mixture:
      self.after_idle(self.updateMolFractionMatrix)  
  
  def updateButtons(self, obj=None):

    molLabelB = self.addMolFracPulldown.getObject()
  
    if self.mixture:
      self.sampleButtons[0].enable()
      self.patternButtons[3].enable()
      
      if molLabelB:
        self.patternButtons[2].enable()
      else:
        self.patternButtons[2].disable()  
 
      if self.molFraction:
        self.patternButtons[0].enable()
        self.patternButtons[1].enable()
      else:
        self.patternButtons[0].disable()
        self.patternButtons[1].disable()

    else:
      self.sampleButtons[0].disable()
      
      for button in self.patternButtons:
        button.disable()
    
      
    if self.resLabel:
      self.residueButtons[0].enable()
      self.residueButtons[1].enable()
      self.atomButtons[0].enable()
    else:
      self.residueButtons[0].disable()
      self.residueButtons[1].disable()
      self.atomButtons[0].disable()
      
    if self.atomLabelTuple:
      self.atomButtons[0].enable()
      self.atomButtons[1].enable()
    else:
      self.atomButtons[0].disable()
      self.atomButtons[1].disable()
  
  
  def administerNotifiers(self, notifyFunc):
    
    clazz = 'ccp.molecule.LabeledMolecule.ResLabelFraction'
    for func in ('__init__', 'delete','setWeight'):
      notifyFunc(self.updateResLabelFractionAfter, clazz, func)

    clazz = 'ccp.molecule.LabeledMolecule.SingleAtomLabel'
    for func in ('__init__', 'delete','setWeight'):
      notifyFunc(self.updateAtomLabelAfter, clazz, func)

    clazz = 'ccp.molecule.LabeledMolecule.UniformAtomLabel'
    for func in ('__init__', 'delete','setWeight'):
      notifyFunc(self.updateAtomLabelAfter, clazz, func)
    
    clazz = 'ccp.molecule.LabeledMolecule.ResLabel'
    for func in ('__init__', 'delete'):
      notifyFunc(self.updateResLabelAfter, clazz, func)

    clazz = 'ccp.molecule.LabeledMolecule.MolLabel'
    for func in ('__init__', 'delete'):
      notifyFunc(self.updateMolLabelAfter, clazz, func)

    clazz = 'ccp.molecule.LabeledMolecule.MolLabelFraction'
    for func in ('__init__', 'delete', 'setWeight'):
      notifyFunc(self.updateMolLabelFractionAfter, clazz, func)

    clazz = 'ccp.molecule.LabeledMolecule.LabeledMixture'
    for func in ('__init__', 'delete', 'addExperiment',
                 'removeExperiment','setExperiments','setName'):
      notifyFunc(self.updateMixtureAfter, clazz, func)

    for func in ('__init__', 'delete', 'setLongName', 'setDetails'):
      for clazz in ('ccp.molecule.ChemCompLabel.LabelingScheme',):
        notifyFunc(self.updateSchemes, clazz, func)

  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
        
    BasePopup.destroy(self)
  
def getMolLabelFromScheme(labelledMolecule, scheme): # Move to core
  """
  Make a new labelled molecule pattern using an isotopomer labelling scheme.

  .. describe:: Input
  
  molecule.LabeledMolecule.LabeledMolecule
  molecule.ChemCompLabel.LabelingScheme  
  
  .. describe:: Output
  
  molecule.LabeledMolecule.MolLabel    
  """  

  molLabel = labelledMolecule.newMolLabel()
  molecule = labelledMolecule.molecule
  
  for molResidue in molecule.molResidues:
    resLabel = molLabel.newResLabel(resId=molResidue.serial)
    updateResLabelFractions(resLabel, scheme)

  return molLabel

def updateResLabelFractions(resLabel, scheme=None):
  """
  Set a residue label, within a labelled molecule pattern, to have
  the same isotopic incorporation as the corresponding isotopomers
  from a labelling scheme. If no scheme is specifed any existing
  linked schemes are used, refereshing isotope incorporations according
  to the weights of any resLabelFractions. Any existing AtomLabel
  specifications are removed.
  
  .. describe:: Input
  
  molecule.ChemCompLabel.LabelingScheme  
  molecule.LabeledMolecule.ResLabel    
  
  .. describe:: Output
  
  None
  """  

  # Delete any atomLabels that override the
  # resLabelFractions?

  molResidue = resLabel.molResidue
  getFraction = resLabel.findFirstResLabelFraction
  
  # Cleanup AtomLabels
  
  for atomLabel in resLabel.atomLabels:
    atomLabel.delete()

  if scheme:
    # Potentially new weights
    
    chemComp = molResidue.chemComp
    chemCompLabel = scheme.findFirstChemCompLabel(chemComp=chemComp)
  
    for resLabelFraction in resLabel.resLabelFractions:
      if resLabelFraction.labelingScheme is scheme:
        if resLabelFraction.isotopomerSerial is None:
          resLabelFraction.delete()
          
      else:
        resLabelFraction.delete()
  
    if chemCompLabel: # Do something else otherwise?
      for isotopomer in chemCompLabel.isotopomers:
        serial = isotopomer.serial
        resLabelFraction = getFraction(isotopomerSerial=serial,
                                       schemeName=scheme.name)
        weight = isotopomer.weight
 
        if resLabelFraction:
          resLabelFraction.weight = weight
        else:
          resLabel.newResLabelFraction(isotopomerSerial=serial,
                                       schemeName=scheme.name,
                                       weight=weight)

  else:
    for resLabelFraction in resLabel.resLabelFractions:
       resLabelFraction.delete()
  


def setResLabelNaturalAbundance(resLabel):
  """
  Set a residue label, within a labelled molecule pattern, to have
  natural abundance isotope incorporation.
  
  .. describe:: Input
  
  molecule.LabeledMolecule.ResLabel    
  
  .. describe:: Output
  
  None
  """  

  molResidue = resLabel.molResidue
  chemCompVar = molResidue.chemCompVar

  # Delete any resLabelFractions, really?
  
  for resLabelFraction in resLabel.resLabelFractions:
     resLabelFraction.delete()

  for atomLabel in resLabel.atomLabels:
    if atomLabel.className == 'SingleAtomLabel':
      atomLabel.delete()

  elements = set([])
  for chemAtom in chemCompVar.chemAtoms:
    if chemAtom.className == 'LinkAtom':
      continue
  
    elements.add(chemAtom.chemElement)
  
  for element in elements:
    symbol = element.symbol
    
    for isotope in element.isotopes:
      abundance = isotope.abundance
      
      if not abundance:
        continue
      
      if abundance < 0.0003:
        continue
    
      massNumber = isotope.massNumber
    
      atomLabel = resLabel.findFirstAtomLabel(elementName=symbol,
                                              massNumber=massNumber)

      if not atomLabel:
        atomLabel = resLabel.newUniformAtomLabel(elementName=symbol,
                                                 massNumber=massNumber)

      atomLabel.weight = 100.0*abundance

def getMolLabelMass(molLabel):
  """
  Get the molecular mass of a molecule with a specific labelling pattern

  .. describe:: Input
  
  molecule.LabeledMolecule.MolLabel    
  
  .. describe:: Output
  
  Float
  """  

  chemElementStore = molLabel.root.currentChemElementStore
  getElement = chemElementStore.findFirstChemElement
  molMass = 0.0

  for resLabel in molLabel.resLabels:
    chemCompVar = resLabel.molResidue.chemCompVar
    uniform = {}
    single = {}
    
    for atomLabel in resLabel.atomLabels:
      if atomLabel.className == 'UniformAtomLabel':
        element = atomLabel.elementName
        if element not in uniform:
          uniform[element] = []
          
        chemElement = getElement(symbol=atomLabel.elementName)
        isotope = chemElement.findFirstIsotope(massNumber=atomLabel.massNumber)
        uniform[element].append((isotope.mass, atomLabel.weight))
        
      else:
        name = atomLabel.atomName
        
        if name not in single:
          single[name] = []
          
        chemElement = atomLabel.chemAtom.chemElement
        isotope = chemElement.findFirstIsotope(massNumber=atomLabel.massNumber)
        single[name].append((isotope.mass, atomLabel.weight))  

    uniformMass = {}
    for element in uniform:
      atomMass = 0.0
      atomFrac = 0.0
 
      for mass, frac in uniform[element]:
        atomFrac += frac
        atomMass += frac*mass
 
      if atomFrac:
        atomMass /= atomFrac
        
      uniformMass[element] = atomMass

    singleMass = {}
    for name in single:
      atomMass = 0.0
      atomFrac = 0.0
 
      for mass, frac in single[name]:
        atomFrac += frac
        atomMass += frac*mass
 
      if atomFrac:
        atomMass /= atomFrac
        
      singleMass[name] = atomMass

    remainingAtoms = set()
    
    for chemAtom in chemCompVar.chemAtoms:
      name = chemAtom.name
      
      if chemAtom.className == 'LinkAtom':
        continue
        
      element = chemAtom.elementSymbol
      
      if name in singleMass:
        molMass += singleMass[name]
      
      elif element in uniformMass:
        molMass += uniformMass[element]
      
      else:
        remainingAtoms.add(name)
    
    if remainingAtoms:
      isotopomerAtoms = set()
      sumFrac = 0.0
    
      for resLabelFraction in resLabel.resLabelFractions:
        allIsotopomers = resLabelFraction.isotopomerSerial == 0
        frac = resLabelFraction.weight
        
        for isotopomer in resLabelFraction.isotopomers:
          if allIsotopomers:
            frac2 = isotopomer.weight
          else:
            frac2 = frac 
          
          sumFrac += frac2
      
      if sumFrac:
        for resLabelFraction in resLabel.resLabelFractions:
          allIsotopomers = resLabelFraction.isotopomerSerial == 0
          frac = resLabelFraction.weight
 
          for isotopomer in resLabelFraction.isotopomers:
            if allIsotopomers:
              frac2 = isotopomer.weight
            else:
              frac2 = frac
            
            frac2 /= sumFrac
            atomData = {}
            for atomLabel in isotopomer.atomLabels:
              name = atomLabel.name
              
              if name not in remainingAtoms:
                continue
 
              if name not in atomData:
                atomData[name] = []
 
              frac3 = atomLabel.weight
              mass = atomLabel.isotope.mass
              atomData[name].append((mass, frac3))
 
            for name in atomData:
              isotopomerAtoms.add(name)
              atomMass = 0.0
              atomFrac = 0.0
 
              for mass, frac in atomData[name]:
                atomFrac += frac
                atomMass += frac*mass
 
              if atomFrac:
                atomMass /= atomFrac
 
                molMass += frac2*atomMass
            
      unaccounted = remainingAtoms - isotopomerAtoms
      for name in unaccounted:
        chemAtom = chemCompVar.findFirstChemAtom(name=name)
        molMass += chemAtom.chemElement.mass # Natural abundance


  return molMass

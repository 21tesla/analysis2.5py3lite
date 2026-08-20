
"""
======================COPYRIGHT/LICENSE START==========================

CalcDistConstraints.py: Part of the CcpNmr Analysis program

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

from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.LabelDivider import LabelDivider
from memops.gui.MessageReporter import showOkCancel, showWarning, showMulti
from memops.gui.ProgressBar import ProgressBar
from memops.gui.PulldownList import PulldownList
from memops.gui.RadioButtons import RadioButtons
from memops.gui.ScrolledMatrix import ScrolledMatrix 
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.core.AssignmentBasic import getShiftLists, getPeakDimResonances
from ccpnmr.analysis.frames.NoeDistParamsFrame import NoeDistParamsFrame
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.ConstraintBasic import makeAmbigDistConstraints, makeDistConstraints
from ccpnmr.analysis.core.ConstraintBasic import findPeakConstraintPairs, getDistMinMax, getNoeDistance
from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes, getThroughSpacePeakLists, getNoesyPeakLists, getOnebondDataDims, getCompatibleHsqcPeakLists, getDiagonalBondedDataDims
from ccpnmr.analysis.core.PeakBasic import getPeakHeight, getPeakVolume, getSimplePeakAnnotation
from ccpnmr.analysis.core.Util import setSpectrumNoeShiftRanges, getSpectrumNoesShiftRanges
from ccpnmr.analysis.core.Util import setSpectrumNoeResidueRanges, getSpectrumNoesResidueRanges
from ccpnmr.analysis.core.Util import getAnalysisDataDim
from ccpnmr.analysis.core.NetworkAnchoring import networkAnchorAssign
                                 
from ccpnmr.analysis.popups.ViewPeakGroups import ViewPeakGroupsPopup

# TBD: UpdateButtons ADCs

hTol = 0.04 
nTol = 0.4  
cTol = 0.4  

TOLERANCE_DICT = {'1H':hTol,'15N':nTol,'13C':cTol,'13C,15N':cTol}

tip = 'Shows a table of '
GROUP_TIP_TEXTS = {'Matched':tip + 'the peaks that were successfully matched to chemical shifts, and which would result in a restraint',
                   'Assigned':tip + 'the peaks that were already assigned to specific resonances and so are not included for shift match restraint generation',
                   'Diagonal':tip + 'the peaks that were excluded because they were too close to the homonuclear diagonal',
                   'Out of range':tip + 'the peaks that were excluded because they only matched assignments were not in the stated residue range ',
                   'Poor Merit':tip + 'the peaks that were excluded because they have a figure-of-merit value below the stated threshold',
                   'Too Distal':tip + 'the peaks that were excluded because they only match assignments that are too far apart given the stated guide structure',
                   'Unmatchable':tip + 'the peaks that were excluded because they do not match any chemical shifts in the stated ranges'}

USE_REF_INTENSITY = 'Use ref intensity'
USE_TABLE_INTENSITY = 'Use table intensities'

class ResidueRowObject:
  def __init__(self, dataDims,chain,startRes,endRes):
    self.dataDims = dataDims
    self.chain = chain
    self.startRes = startRes
    self.endRes = endRes

class ToleranceRowObject:
  def __init__(self,dataDimInfo,minTol,multiplier,maxTol):
    self.dataDimInfo = dataDimInfo
    self.minTol = minTol
    self.multiplier = multiplier
    self.maxTol = maxTol

class ShiftRowObject:
  def __init__(self,dataDim,isotopes,start,end):
    self.dataDim = dataDim
    self.isotopes = isotopes
    self.start = start
    self.end = end

class CalcDistConstraintsPopup(BasePopup):
  """
  **Create Structural Distance Restraints From NMR Peaks**
  
  This system is used to make resonance to resonance distance restraints based
  on peak intensity information from through-space NMR experiments, such as
  NOESY. Three mechanisms are available with regard to the assignment of
  resonance pairs to individual restraints, but for all three methods the
  distance information is extracted from peak intensities in the same way.  The
  use of peak intensities and the three restraint generation methods are
  described below.

  **Restraint Distance Estimation**

  The options available in the "Restraint Distance Params section" if the
  "Settings" tab allow the user to specify how through-space peak intensities
  relate to the distance bounds of any any generated distance restraints. The
  default method is to calculate a target distance as peak height raised to the
  power of -1/6 multiplied by some scaling factor. However the user may also use
  peak volumes by changing the "Intensity Type" and other power relations 
  via the "Distance function" pulldown menu.
  
  The peak intensity scaling factor is defined such that the reference intensity
  (defaults to the peak list's average volume) exactly corresponds to the
  reference distance (default 3.2 Angstroms). The user may change either of
  these values to calibrate in a different way. For example by measuring the
  peak intensities of known atomic distances such as those in an aromatic ring.
  The default values are appropriate for most protein structure calculations
  based on NOESY experiments, and in any case structure calculation programs
  like ARIA and CYANA will normally re-calibrate peak intensity scaling during
  the structure generation process. Hence, the initial calibration is not too
  critical; covalent geometry will provide better global calibration.

  The upper and lower bounds of the distance restraints are calculated as
  fractional differences from the calculated target distance. For example if a
  target distance is 3.2 Angstroms an upper fractional error of 0.20 gives an
  upper limit of 3.84 (20% more) and 0.2 lower fractional error gives a lower
  limit of 2.56 (20% less). The upper and lower bounds are also confined the
  upper and lower distance limits. These are absolute minimum and maximum values
  for the bounds (1.72 & 8.00 Angstroms respectively by default).

  If none of the fractional power distance functions are appropriate for the NMR
  data then the user may select the "Distance Bins" option. Here the user
  defines peak intensity classes, each of which relates to a target distance and
  distance bounds. Each class has a "Min. Peak Intensity" value, which is a
  multiplying factor relative to the reference intensity (by default the peak
  list mean). When Analysis determines which distance class to use for a peak it
  first divides its intensity value by the reference intensity to get a relative
  value, then the distance classes are searched for the class with the largest
  "Min. Peak Intensity" that is below the peak's relative intensity. For
  example with two classes defined at 0.0 and 1.0 minimum intensities, a peak
  with a relative intensity value of 0.5 would fall in the first class, and a
  peak with 1.2 would fall in the second.

  **Restraints From Fully Assigned Peaks**
  
  The simplest way of generating distance restraints from through-space NMR
  experiments is to use the spectrum peaks that have been fully assigned. If a
  peak is caused by correlation between resonances and the atomic identities of
  those causal resonances is known (or at least suspected) a distance restraint
  can be created from the peak information. In CCPN terms the restraint will be
  between the resonances that are assigned to the through-space correlated peak
  dimensions, but only if those resonances are assigned to particular atoms (in
  a molecular system). Restraints are made for assigned peaks by selecting an
  appropriate peak list in the pulldown menu, at the very top of the popup,
  above the tabs, and pressing the [Make Assigned Restraints] at the bottom.
  Naturally, the user should make sure that parameters are specified for
  distance estimation and residue ranges and that the correct restraint set is
  selected in the "Settings" tab, to say where the generated restraint list
  should be placed.

  Restraints will be generated for peaks that carry multiple, ambiguous
  assignments. This happens where the observed peak in a spectrum covers
  multiple resonance intersections that are unresolvable due to the relative
  scale of linewidths compared to the similarity of chemical shifts. In such
  circumstances all of the resonances assigned to the peaks through-space
  correlated dimensions will appear in the restraint. By default, if a peak is
  assigned to resonances A & B in one dimension and C & D on the other then an
  ambiguous restraint will be generated with items for all four resonance pairs;
  A-C, A-D, B-C and B-D, and a restraint would be satisfied by any of these
  pairs being within distance bounds. However, if the user knows that only
  certain pairings are possible on an ambiguously assigned peak then the
  assignments may be grouped via the `Assignment Panel`_. Referring to the
  previous example if resonance A & C are in group 1 and B & D are in group 2
  then the generated restraint will only have two items; A-C and B-D.

  **Restraints From Matching Chemical Shifts** 
  
  There is a second common way to generate distance restraints, which is to
  match the chemical shifts of resonances to unassigned peak positions, thus
  generating potentially highly ambiguous distance restraints. Such restraints
  would typically be filtered to select only the correct contributing resonance
  pairs, by iterative structure generation and violation analysis in a program
  like ARIA or CYANA. Often the user can leave the matching of chemical shifts
  to the structure generation program, by passing the program peak lists rather
  than restraint lists. However, it is also possible to make such restraints in
  CCPN.  The [Make Shift Match Restraints] can be used to generate restraints
  for the unassigned peaks in the selected peak list, however the user should be
  aware of all of the parameters that will be applied, including the "Shift
  Match Tolerances", "Residue Ranges" and "Chem Shift Ranges".

  The "Shift Match Tolerances" tab is important because it specifies how the
  matching is done between (unassigned) peak positions and resonances with known
  chemical shifts. In essence, every unassigned peak dimension is compared to
  the chemical shifts contained within the shift list that is associated  with
  the NMR experiment from which the peak derive. The shift list that relates to
  an experiment, and hence to a peak list, may be adjusted in the main
  Experiments_ table. Resonances whose shifts are close to the peak position may
  be the cause of that peak and so can be used to make a restraint. Often
  several resonances will match in each peak dimension, resulting in a restraint
  with many alternative, ambiguous resonance pairings (restraint "items"). The
  incorrect possibilities may be removed by automated structure generation
  programs, which generate guide structures based on the unambiguous restraints.
  The "... PPM tolerance" values are used to control the maximum allowed
  difference from a peak location to a chemical shift for a resonance to be
  considered as an assignment possibility; the threshold width of the shift
  matching. In the absence of linewidth information only the "Min/default PPM
  tolerance" is used. If a peak has a recorded linewidth then its value may be
  used to scale the shift-match threshold, i.e. wider peaks match a wider PPM
  range. In such cases the shift-match threshold is the peak's linewidth
  multiplied by the "line width multiplier" value, but bounded by the minimum
  and maximum PPM tolerance values.

  The [Test Shift Match] function can be used to trail the peak to resonance
  matching procedure without generating any distance restraints. This is useful
  to test how many peaks would fail to be matched, and the resonance why.

  When shift-match restraints are made the resonance matches can be filtered by
  using a guide structure (e.g. homologue or rough structure) with a
  corresponding distance cutoff, so that resonance pairs that are far apart  in
  the guide will not be places in the restraints. Also, shift matching obeys any
  isotope labelling scheme setting; any resonance pair that has a spin-active
  isotope incorporation below the "Minimum Isotope Fraction" is deemed not to be
  visible and will not be present in a restraint.
  
  **Network Anchoring**
  
  In the case of the shift-matching method potentially ambiguous distance
  restraints are generated by simply matching peak positions to close chemical
  shifts. In the case of network anchoring method, chemical shifts are also
  matched to peaks, but the ambiguous possibilities are refined by selecting
  only NOE assignments from amongst the possibilities that are supported by
  other known close resonances or covalent structure. Say, for example, that a
  peak could have arisen from a number of resonance pairs with matching shifts.
  Two resonances A & B are more likely to be a correct assignment for the peak
  if we know that they are close to (or bound to) the same intermediary
  resonance, C and therefore must be close to each other. 

  The network anchoring procedure is directed using the "Network Anchoring" tab
  and   operates over several peak lists at the same time. These peak lists are
  selected in the "Peak List Selection" table by double clicking in the "Use?"
  column to indicate "Yes". Any assigned peaks in these peak lists, together
  with the covalent connectivity of the atoms, will form the initial
  through-space network. The restraint generation occurs after the iterative
  network refinement process converges (usually in fewer than five cycles). The
  parameters for the refinement are specified above the main table. It should be
  noted that only a shift list and set of shift match tolerances is currently
  used.

  The "Min net. score" is used during refinement to dictate how much network
  support an assignment possibility must have for it to be included. Setting
  this to a higher value will mean fewer mistakes are made but at the expense of
  making fewer restraints. The "Strictness" setting relates to which residues
  may be involved in the local network to help resolve restraint possibilities.
  Considering only connections within a single residue makes few mistakes but
  doesn't resolve much of the ambiguity. Using any residue makes most mistakes
  but most restraints. The default "Normal: Sequential residue support"
  is a useful compromise between these two extremes.

  **Residue and Chemical Shift Ranges**
  
  The "Residue Ranges" section is used to restrict the restraint generation to
  only specific regions of a molecular system.  For example, this is useful if a
  protein is known to have an unstructured tail and the chemical shifts from
  this should not be considered because they would add unnecessary ambiguity. A
  residue range applies to one molecular chain and the peak dimensions on one
  side of a through-space correlation. This means that different dimensions can
  be restricted to different bits of a molecular system; useful with X-filtered
  NMR experiments. Residue ranges are applied to all forms of restraint
  generation, i.e. for assigned peaks, shift matching and network anchoring.

  The "Chem Shift Ranges" tab is used to allow or disallow certain resonance
  signals from being used in the shift matching and network anchoring restraint
  generation. Different allowed ranges can be specified for each spectrum
  dimension (of the selected peak list) and a single dimension can have multiple
  ranges. In this way the user can control whether certain kind of resonances
  with distinct shifts are used, e.g. methyls. Often a 'water notch' is used,
  where two ranges with a small gap between are constructed to exclude peaks
  near the chemical shift of water; such peaks do not related to close macromolecule
  distances.
  
  **Caveats & Tips**
  
  Restraints are always made in a restraint list that is placed in the selected
  restraint set. By default a new restraint set is made, thus to  make restraint
  lists that go together in an existing set the right  one must be selected in
  the "Settings Tab".

  A new restraint set should always be used if any resonance to atom assignments
  have been changed, otherwise new restraints will still obey the original
  "fixed" assignments that were present in the old restraint set, not the new
  ones.

  .. _`Assignment Panel`: EditAssignmentPopup.html
  .. _Experiments: EditExperimentPopup.html
  
  """

  def __init__(self, parent, *args, **kw):

    self.guiParent = parent
  
    self.aliasing = True
    self.ignoreDiagonals = True
    self.spectrum = None
    self.peakList = None
    self.molSystem = None
    self.minPeakMerit = None
    self.residueRanges = []
    self.residueRange = None
    self.tolerances = []
    self.chemShiftRanges = []
    self.chemShiftRange = None
    self.constraintSet = None
    self.netPeakList = None
    self.netShiftList = None
    self.labellingScheme = True
    self.structure = None
    self.diagPeakMap = {}
    self.hsqcPeakMap = {}
    self.whichDistanceMap = {}

    BasePopup.__init__(self, parent=parent, title='Structure : Make Distance Restraints', **kw)

  def body(self, guiFrame):

    self.geometry('800x500')
  
    guiFrame.grid_columnconfigure(1, weight=1)
    guiFrame.grid_rowconfigure(1, weight=1)
    
    # Table widgets

    self.resDimsPulldown  = PulldownList(self, callback=self.setResDim)
    self.resChainPulldown = PulldownList(self, callback=self.setResChain)
    self.resStartEntry    = IntEntry(self, returnCallback=self.setResStart, width=5)
    self.resEndEntry      = IntEntry(self, returnCallback=self.setResEnd, width=5)
    
    self.tolMinPpmEntry     = FloatEntry(self, returnCallback=self.setTolMinPpm, width=5)
    self.tolMultiplierEntry = FloatEntry(self, returnCallback=self.setTolMultiplier, width=5)
    self.tolMaxPpmEntry     = FloatEntry(self, returnCallback=self.setTolMaxPpm, width=5)

    self.shiftDataDimPulldown = PulldownList(self, callback=self.setShiftDataDim)
    self.shiftStartEntry  = FloatEntry(self, returnCallback=self.setShiftStart, width=8)
    self.shiftEndEntry    = FloatEntry(self, returnCallback=self.setShiftEnd, width=8)
    

    frame = Frame(guiFrame, grid=(0,0))
    tipText = 'Selects which peak list is used in the generation of distance restraints'
    label = Label(frame, text='Peak List: ', grid=(0,0))
    self.peakListPulldown = PulldownList(frame, grid=(0,1), tipText=tipText,
                                         callback=self.changePeakList)

    frame = Frame(guiFrame, grid=(0,1))
    tipText = 'Whether simulated peak lists are included in peak lists'
    label = Label(frame, text='  Include simulated: ', grid=(0,0))
    self.includeSimulatedSelect = CheckButton(frame, grid=(0,1), tipText=tipText,
                                              callback=self.updatePeakLists)

    buttons = UtilityButtonList(guiFrame, helpUrl=self.help_url, grid=(0,2))
    
    #tipTexts = ['General settings that specify how restraints are made',
    #            'Allows restriction of restraints to specific sections of molecular sequence(s)',
    #            'Setup thresholds for generating restraints by matching peak positions to chemical shifts',
    #            'Allows restraints made with shift-peak matching to use restricted chemical shift values',
    #            'Use NOE network information to make shift-peak match generation of restraints less ambiguous']
    options = ['Settings','Residue Ranges','Shift Match Tolerances',
               'Chem Shift Ranges','Network Anchoring',]
    self.tabbedFrame = TabbedFrame(guiFrame, options=options)#, tipTexts=tipTexts)
    self.tabbedFrame.grid(row=1, column=0, columnspan=3, sticky='nsew')
    optFrame, resFrame, tolFrame, shiftFrame, netFrame = self.tabbedFrame.frames
    
    optFrame.grid_columnconfigure(3, weight=1)
    optFrame.grid_rowconfigure(3, weight=1)
    resFrame.grid_columnconfigure(1, weight=1)
    resFrame.grid_rowconfigure(1, weight=1)
    shiftFrame.grid_columnconfigure(0, weight=1)
    shiftFrame.grid_rowconfigure(0, weight=1)
    netFrame.grid_columnconfigure(0, weight=1)
    
    # Opt frame

    tipText = 'Selects the restraint set in which to make new restraint lists'
    label = Label(optFrame,text='Restraint Set: ', grid=(0,0))
    self.constraintSetPulldown = PulldownList(optFrame, grid=(0,1), tipText=tipText,
                                              callback=self.changeConstraintSet)

    tipText = 'Selects the isotope labelling scheme, if any, to filter restraint possibilities'
    label = Label(optFrame,text='Isotope Labelling: ', grid=(0,2))
    self.labellingSchemePulldown = PulldownList(optFrame, grid=(0,3), tipText=tipText,
                                                callback=self.changeLabellingScheme)
  
    tipText = 'Threshold for peak figure-of-merit value, below which peaks will not be considered for making restraints'
    label = Label(optFrame, text='Minimum\nPeak Merit: ', grid=(1,0))
    self.minMeritEntry = FloatEntry(optFrame, text='1.0', grid=(1,1), tipText=tipText,
                                    returnCallback=self.setMinPeakMerit, width=8)
   
    tipText = 'If using an isotope labelling scheme, the minimum proportion of spin active isotope to consider a resonance for restraints'
    label = Label(optFrame, text='Minimum\nIsotope Fraction: ', grid=(1,2))
    self.minIsoFractionEntry = FloatEntry(optFrame, text='0.1', width=8, grid=(1,3), tipText=tipText)
    
    
    tipText = 'Selects a structure that can be used to remove distal resonance pairs from distance restraints'
    label = Label(optFrame,text='Shift Match\nGuide Structure: ', grid=(0,4))
    self.structurePulldown = PulldownList(optFrame, grid=(0,5), tipText=tipText,
                                          callback=self.changeStructure)
  
    tipText = 'The maximum allowed distance if filtering restraint possibilities using a structure'
    label = Label(optFrame, text=u'Max Structure\nDistance (\u00C5): ', grid=(1,4))
    self.maxDistEntry = FloatEntry(optFrame, text='20.0', grid=(1,5), width=8, tipText=tipText)
    

    self.noeFrame = NoeDistParamsFrame(optFrame, self.nmrProject,
                                       intensityTypeCallback=self.updatePeakFrameMatrix,
                                       distanceFunctionCallback=self.updatePeakFrameMatrix)
    self.noeFrame.grid(row=2, column=0, columnspan=6, sticky='ew')
     
    peakFrame = self.peakFrame = LabelFrame(optFrame, text='Peak Normalisation')
    peakFrame.grid(row=3, column=0, columnspan=6, sticky='nsew')
    peakFrame.grid_rowconfigure(3, weight=1)
    peakFrame.grid_columnconfigure(1, weight=1)

    label = Label(peakFrame, text='Table below allows you to specify alternatives to ref intensity (cells highlighted in green are used)', grid=(0, 0), gridSpan=(1,2))
    
    tipTexts = ('Always use the ref intensity no matter what the table below says', 'Use the setting for each peak as specified in the table below')
    self.intensityChoiceButton = RadioButtons(peakFrame, entries=(USE_REF_INTENSITY, USE_TABLE_INTENSITY),
                                              tipTexts=tipTexts, grid=(1, 0))
    
    tipTexts = ['Update the table (but most changes should be recognised automatically)',
                'Set distance category to be used for selected rows']
    texts = ['Refresh\nTable', 'Set Distance\nTo Use...']
    commands = [self.updatePeakFrameMatrix, self.setWhichDistance]
    buttons = ButtonList(peakFrame, commands=commands, texts=texts, tipTexts=tipTexts)
    buttons.grid(row=2, column=0, sticky='w')

    frame2 = Frame(peakFrame)
    frame2.grid(row=2, column=1, sticky='e')
    frame2.grid_rowconfigure(0, weight=1)
    frame2.grid_columnconfigure(1, weight=1)
    
    label = Label(frame2, text='HSQC Peak List: ')
    label.grid(row=0, column=0, sticky='e')

    tipText = 'Selects which HSQC peak list is used for analysis'
    self.hsqcPeakListPulldown = PulldownList(frame2, grid=(0,1), tipText=tipText,
                                             callback=self.changeHsqcPeakList)

    self.whichUsePulldown = PulldownList(self, callback=self.setWhichUse)
    tipTexts = ['The peak',
                'The peak annotation',
                'The peak intensity (height or volume)',
                'The diagonal peak',
                'The diagonal peak annotation',
                'The diagonal peak intensity (height or volume)',
                'The HSQC peak',
                'The HSQC peak annotation',
                'The HSQC peak intensity (height or volume)',
                'The distance derived using reference intensity',
                'The distance derived using diagonal intensity',
                'The distance derived using HSQC intensity',
                'Which intensity to use for normalisation']
    headingList = ['Peak',
                   'Peak\nAnnotation',
                   'Distance\nTo Use',
                   'Ref Derived\nDistance',
                   'Diag Derived\nDistance',
                   'HSQC Derived\nDistance',
                   'Peak\nIntensity',
                   'Diagonal\nPeak',
                   'Diagonal Peak\nAnnotation',
                   'Diagonal Peak\nIntensity',
                   'HSQC Peak',
                   'HSQC Peak\nAnnotation',
                   'HSQC Peak\nIntensity']
    n = len(headingList) - 1
    n0 = 2
    n1 = n - n0
    editWidgets = n0*[None] + [self.whichUsePulldown] + n1*[None]
    editGetCallbacks = n0*[None] + [self.getWhichUse] + n1*[None]
    editSetCallbacks = n0*[None] + [self.setWhichUse] + n1*[None]
    self.peaksMatrix = ScrolledMatrix(peakFrame, multiSelect=True,
                                      headingList=headingList,
                                      editWidgets=editWidgets, tipTexts=tipTexts,
                                      editSetCallbacks=editSetCallbacks,
                                      editGetCallbacks=editGetCallbacks)
    self.peaksMatrix.grid(row=3, column=0, columnspan=2, sticky='nsew')

    # Res frame
 
    tipText = 'Selects which molecular chain to set sequence/residue ranges for'
    label = Label(resFrame, text='Molecular system: ', grid=(0,0))
    self.molSysPulldown = PulldownList(resFrame, callback=self.changeMolSystem,
                                       grid=(0,1), tipText=tipText)
    
    tipTexts = ['Which spectrum dimensions the residue range filter applies to',
                'The molecular chain of the allowed residue range',
                'The first allowable residue position to make restraints for',
                'The last allowable residue position to make restraints for']
    headingList = ['Dimensions','Chain','Start','End']
    editWidgets = [self.resDimsPulldown, self.resChainPulldown, self.resStartEntry, self.resEndEntry]
    editGetCallbacks = [self.getResDim, self.getResChain, self.getResStart, self.getResEnd]
    editSetCallbacks = [self.setResDim, self.setResChain, self.setResStart, self.setResEnd]
    self.rangeMatrix = ScrolledMatrix(resFrame, headingList=headingList,
                                      callback=self.selectResidueRange, initialRows=2,
                                      editWidgets=editWidgets, tipTexts=tipTexts,
                                      editSetCallbacks=editSetCallbacks,
                                      editGetCallbacks=editGetCallbacks,
                                      deleteFunc=self.deleteResidueRange)
    self.rangeMatrix.grid(row=1, column=0, columnspan=2, sticky='nsew')
  
    tipTexts = ['Add a new allowable residue range to the table',
                'Remove the currently selected residue range']
    texts = ['Add Residue Range','Delete Residue Range']
    commands = [self.addResidueRange,self.deleteResidueRange]
    self.resButtons = ButtonList(resFrame, commands=commands, texts=texts, tipTexts=tipTexts)
    self.resButtons.grid(row=2, column=0, columnspan=2, sticky='ew')

    
    # Tol frame
    tolFrame.expandGrid(1,4)
  
    tipText = 'Whether peaks that are matched to shifts might be aliased one or more sweep withs; i.e. not at their real PPM'
    label = Label(tolFrame, text='Consider aliased shifts?', grid=(0,0))
    self.aliasSelect = CheckButton(tolFrame, callback=self.setAlias, grid=(0,1), tipText=tipText)
    self.aliasSelect.set(1)
  
    tipText = 'Whether to avoid matching shifts to diagonal peaks'
    label = Label(tolFrame, text='Ignore diagonal peaks?', grid=(0,2))
    self.diagSelect = CheckButton(tolFrame, callback=self.setDiag, grid=(0,3), tipText=tipText)
    self.diagSelect.set(1)
    
    tipText = 'Copy shift match tolerances settings from a different spectrum (overwrites old values)'
    self.copySpecButton = Button(tolFrame, text='Copy From:', tipText=tipText,
                                 command=self.copySpecTols, bd=1, grid=(0,5))
    self.copySpecPulldown = PulldownList(tolFrame, callback=None, tipText=tipText, grid=(0,6))
    
    tipTexts = ['The spectrum dimension that the tolerance applies to',
                'If a peak has no recorded linewidth, how far in ppm shift values can be from peak positions to make a restraint. With a linewidth, the absolute lower limit to this ppm difference ',
                'A factor that governs how the peak-shift match tolerance scales with peak linewidth',
                'If a peak has a recorded linewidth, the absolute upper limit to the peak-shift ppm difference']
    headingList = ['Dimension','Min/default PPM\ntolerance',
                   'line width\nmultiplier','Max PPM\ntolerance']
    editWidgets = [None, self.tolMinPpmEntry, self.tolMultiplierEntry, self.tolMaxPpmEntry]
    editGetCallbacks = [None, self.getTolMinPpm, self.getTolMultiplier, self.getTolMaxPpm]
    editSetCallbacks = [None, self.setTolMinPpm, self.setTolMultiplier, self.setTolMaxPpm]
    self.toleranceMatrix = ScrolledMatrix(tolFrame, headingList=headingList, 
                                          editWidgets=editWidgets, tipTexts=tipTexts,
                                          editSetCallbacks=editSetCallbacks,
                                          editGetCallbacks=editGetCallbacks,
                                          grid=(1,0), gridSpan=(None,7))
    # Shift frame

    tipTexts = ['The spectrum dimension to which an allowed chemical shift range applies',
                'The isotope(s) relevant to the spectrum dimension',
                'The lower bound of the allowed chemical shift range for matching shifts to peaks',
                'The upper bound of the allowed chemical shift range for matching shifts to peaks']
    headingList = ['Dimension','Isotope','Start','End']
    editWidgets = [self.shiftDataDimPulldown, None, self.shiftStartEntry,self.shiftEndEntry]
    editGetCallbacks = [self.getShiftDataDim, None, self.getShiftStart, self.getShiftEnd]
    editSetCallbacks = [self.setShiftDataDim, None,self.setShiftStart, self.setShiftEnd]
    self.shiftsMatrix = ScrolledMatrix(shiftFrame, headingList=headingList,
                                       callback=self.selectChemShiftRange,
                                       editWidgets=editWidgets, tipTexts=tipTexts,
                                       editSetCallbacks=editSetCallbacks,
                                       editGetCallbacks=editGetCallbacks,
                                       deleteFunc=self.deleteChemShiftRange)
    self.shiftsMatrix.grid(row=0, column=0, sticky='nsew')
  
    tipTexts = ['Add a new chemical shift range for matching peak positions to chemical shifts',
                'Remove the selected chemical shift range from the table']
    texts = ['Add Chem Shift Range','Delete Chem Shift Range']
    commands = [self.addChemShiftRange,self.deleteChemShiftRange]
    self.shiftsButtons = ButtonList(shiftFrame, commands=commands, texts=texts, tipTexts=tipTexts)
    self.shiftsButtons.grid(row=1, column=0, sticky='ew')
  
    # Net anchor frame
    
    frame = Frame(netFrame)
    frame.grid(row=0, column=0, sticky='ew')
    frame.grid_columnconfigure(5, weight=1)
    
    tipText = 'The threshold NOE network sore, below which an assignment possibility will not be considered'
    label = Label(frame, text='Min net. score:', grid=(0,0))
    self.thresholdEntry = FloatEntry(frame, text=1.0, width=6, tipText=tipText, grid=(0,1))

    tipText = 'Whether to assign the peaks of the input peak lists according to the filtered assignment possibilities'
    label = Label(frame, text=' Assign Peaks:', grid=(0,2), sticky='e')
    self.peakListCheck = CheckButton(frame, tipText=tipText, selected=False, grid=(0,3))
    
    tipText = 'Selects which kinds of NOE network information can filter assignment options;\n' + \
              'higher strictness means fewer errors but fewer restraints'
    label = Label(frame, text='Strictness:')
    label.grid(row=0, column=4, stick='w')
    entries = ['Low: Any residue support',
               'Normal: Sequential residue support',
               'High: Same residue support']
    self.strictnessPulldown = PulldownList(frame, texts=entries, grid=(0,5),
                                           callback=None, index=1, tipText=tipText)

    label = Label(frame, text='Tolerances:', grid=(1,0))
    subFrame = Frame(frame)
    subFrame.grid(row=1, column=1, columnspan=5, sticky='w')
    
    tipText = 'When generating restraint possibilities, the maximum difference between peak position and 1H chemical shift'
    label = Label(subFrame, text='1H', grid=(0,0))
    self.tolHydrogenEntry = FloatEntry(subFrame, text=TOLERANCE_DICT['1H'],
                                       grid=(0,1), tipText=tipText, width=7)

    tipText = 'When generating restraint possibilities, the maximum difference between peak position and 15N chemical shift'
    label = Label(subFrame, text='15N', grid=(0,2))
    self.tolNitrogenEntry = FloatEntry(subFrame, text=TOLERANCE_DICT['15N'], tipText=tipText,
                                       grid=(0,3), width=7)

    tipText = 'When generating restraint possibilities, the maximum difference between peak position and 13C chemical shift'
    label = Label(subFrame, text='13C', grid=(0,4))
    self.tolCarbonEntry = FloatEntry(subFrame, text=TOLERANCE_DICT['13C'],
                                     grid=(0,5), tipText=tipText, width=7)
    
    tipText = 'Selects which list of chemical shifts to compare to peak positions'
    label = Label(subFrame, text='Shift List', grid=(0,6))
    self.shiftListPulldown = PulldownList(subFrame, tipText=tipText, grid=(0,7),
                                          callback=self.changeNetShiftList)
    
    div = LabelDivider(netFrame, text='Peak List Selection', grid=(1,0))
    
    netFrame.grid_rowconfigure(2, weight=1)
    
    tipTexts = ['The experiment: spectrum of the input peak list',
                'Peak list serial number',
                'Whether to use this peak list in the NOE network and restraint generation',
                'NUmber of peaks in the peak list']
    headingList = ['Spectrum','PeakList','Use?','Num Peaks']
    editWidgets = [None,None,None,None]
    editSetCallbacks = [None,None,None,None]
    editGetCallbacks = [None,None,self.toggleUseNetPeakList,None]
    self.networkPeakListsMatrix = ScrolledMatrix(netFrame, headingList=headingList,
                                         callback=self.selectNetPeakList,
                                         editWidgets=editWidgets, tipTexts=tipTexts,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks)
    self.networkPeakListsMatrix.grid(row=2, column=0, sticky='nsew') 
    
  
    # Outside buttons
  
    tipTexts = ['Make restraints using the assigned peaks in the selected peak list, according to the settings',
                'Make restraints by matching peaks in the selected list to chemical shifts',
                'Make an analysis of peak to chemical shift matches, without committing any restraints',
                'Run the network anchoring calculation to generate restraints']
    texts = ['Make Assigned\nRestraints',
             'Make Shift Match\nRestraints',
             'Test\nShift Match',
             'Network Anchor\nRestraints']
    commands = [self.calcNormRestraints,
                self.calcShiftMatchRestraints,
                self.testShiftMatch,
                self.runAnchoring]
    self.bottomButtons = ButtonList(guiFrame, commands=commands, texts=texts,
                                    grid=(2,0), gridSpan=(1,3), tipTexts=tipTexts)

    # Main

    self.updatePeakLists()
    self.updateMolSystems()
    self.updateConstraintSets()
    self.updateStructures()
    self.updateLabellingSchemes()
    
    self.updateResidueRanges()
    self.updateTolerances()
    self.updateChemShiftRanges()
    self.updateNetShiftListPulldown()
    self.updateNetPeakLists()

    self.updatePeakFrame()

    self.administerNotifiers(self.registerNotify)

  def administerNotifiers(self, notifyFunc):

    for func in ('__init__','delete'):
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.PeakList', func)
      notifyFunc(self.updateNetPeakLists, 'ccp.nmr.Nmr.PeakList', func)
      notifyFunc(self.updateConstraintSets, 'ccp.nmr.NmrConstraint.NmrConstraintStore', func)
      notifyFunc(self.updateStructures, 'ccp.molecule.MolStructure.StructureEnsemble', func)
      notifyFunc(self.updateMolSystems, 'ccp.molecule.MolSystem.MolSystem', func)
      notifyFunc(self.updateChain, 'ccp.molecule.MolSystem.Chain', func)
    notifyFunc(self.updateResidue,'ccp.molecule.MolSystem.Residue', 'setSeqCode')
 
    for func in ('__init__','delete','setName'):
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.Experiment', func)
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.DataSource', func)
      notifyFunc(self.noeFrame.binsFrame.updateSpectra, 'ccp.nmr.Nmr.Experiment', func)
      notifyFunc(self.noeFrame.binsFrame.updateSpectra, 'ccp.nmr.Nmr.DataSource', func)

    for func in ('setMolSystems','addMolSystem','removeMolSystem'):
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.Experiment', func)
    
    for func in ('__init__','delete'):
      notifyFunc(self.updateNetShiftListPulldown, 'ccp.nmr.Nmr.ShiftList', func)
    notifyFunc(self.updateNetShiftListPulldown, 'ccp.nmr.Nmr.Experiment', 'setShiftList')
    
    for func in ('delete','__init__'):
      notifyFunc(self.updateLabellingSchemes, 'ccp.molecule.ChemCompLabel.LabelingScheme', func)

  def open(self):
    
    self.updatePeakLists()
    self.updateMolSystems()
    self.updateConstraintSets()
    self.updateStructures()
    self.updateLabellingSchemes()

    self.updateResidueRanges()
    self.updateTolerances()
    self.updateChemShiftRanges()
    self.updateNetShiftListPulldown()
    self.updateNetPeakLists()
    
    self.updatePeakFrame()

    BasePopup.open(self)

  def selectNetPeakList(self, peakList, row, col):

    if peakList:
      self.netPeakList = peakList
      
  def toggleUseNetPeakList(self, *opt):

   if self.netPeakList:
     if not hasattr(self.netPeakList,'networkAnchorUse'):
       self.netPeakList.networkAnchorUse = False
     
     self.netPeakList.networkAnchorUse = not self.netPeakList.networkAnchorUse
     self.updateNetPeakLists()

  def changeNetShiftList(self, shiftList):
    
    if shiftList is not self.netShiftList:
      self.netShiftList = shiftList
      self.updateNetPeakLists()

  def updateNetShiftListPulldown(self, obj=None):

    index = 0
    names = []
    shiftLists = getShiftLists(self.nmrProject)
    shiftList = self.netShiftList
    
    
    if shiftLists:
      if shiftList not in shiftLists:
        shiftList = shiftLists[0]
    
      index = shiftLists.index(shiftList)
      names = ['%d' % sl.serial for sl in shiftLists]
    
    else:
      shiftList = None
     
    if self.netShiftList is not shiftList:  
      self.changeNetShiftList(shiftList)
    elif obj:
      self.updateNetPeakLists()
    
    self.shiftListPulldown.setup(names, shiftLists, index)  

  def updateNetPeakLists(self, obj=None):

    textMatrix = []
    objectList = []
    colorMatrix = []

    excludeSimulated = not self.includeSimulatedSelect.get()
    peakLists = getThroughSpacePeakLists(self.project, excludeSimulated)
    for peakList in peakLists:
      if not hasattr(peakList,'networkAnchorUse'):
        peakList.networkAnchorUse = False
        
      spectrum = peakList.dataSource
      experiment = spectrum.experiment
      
      if experiment.shiftList is not self.netShiftList:
        continue
      
      if experiment.findFirstExpTransfer(isDirect=False):
        continue
      
      name = '%s:%s' % (experiment.name,spectrum.name)
        
      boolean = peakList.networkAnchorUse
      inUse = 'No'
      colors = [None, None, '#FFB0B0']
      if boolean:
        inUse = 'Yes' 
        colors = [None, None, '#B0FFB0']       

      datum = [name, peakList.serial, inUse, len(peakList.peaks)]
  
      textMatrix.append(datum)
      objectList.append(peakList)
      colorMatrix.append(colors)

    self.networkPeakListsMatrix.update(colorMatrix=colorMatrix,
                                       objectList=objectList,
                                       textMatrix=textMatrix)
     
  def calcAdcs(self):
    
    if self.spectrum and self.constraintSet:
    
      spectrum = self.spectrum
      from ccpnmr.clouds.ResonanceIdentification import makeNoeAdcs 
    
      resonances = spectrum.topObject.sortedResonances()
      
      print "making ADCs"
      adcList    = makeNoeAdcs(resonances, spectrum,
                               self.constraintSet
                               ,allowedAtomTypes=['H',])
                               
      self.parent.browseConstraints(adcList)

    
  def setAlias(self, boolean):
  
    self.aliasing = boolean

  def setDiag(self, boolean):
  
    self.ignoreDiagonals = boolean

  def getResidueDimSelection(self):
  
    dims = []
    if self.spectrum:
      bonded   = {}
      dataDims = getOnebondDataDims(self.spectrum)
 
      for dataDim1, dataDim2 in dataDims:
        if (not dataDim1.expDim.expDimRefs) or (not dataDim2.expDim.expDimRefs):
          continue
        if ('1H' not in self.getDataDimIsotopes(dataDim1)) and ('1H' not in self.getDataDimIsotopes(dataDim2)):
          continue

        nameI = self.getDataDimName(dataDim1)
        nameJ = self.getDataDimName(dataDim2)
 
        dims.append( ['%s - %s' % (nameI,nameJ), [dataDim1,dataDim2]  ] )
        bonded[dataDim1] = True
        bonded[dataDim2] = True

      for dataDim in self.spectrum.sortedDataDims():
        if bonded.get(dataDim) is None:
          if dataDim.expDim.expDimRefs:
            dims.append( [ self.getDataDimName(dataDim),[dataDim,] ] )

    return dims

  def setResDim(self, dimSelection):

    i = self.resDimsPulldown.getSelectedIndex()
    
    if self.residueRange:
      self.residueRange.dataDims = tuple(self.getResidueDimSelection()[i])
    
    self.updateResidueRanges()
  
  def getResDim(self, resRange):
    
    if resRange:
      dimList = self.getResidueDimSelection()
      names = [x[0] for x in dimList]
      self.resDimsPulldown.setup(names, dimList, dimList.index(list(resRange.dataDims)))

  def getChains(self):
  
    chains = []
    if self.molSystem:
      for chain in self.molSystem.sortedChains():
        if chain.residues:
          chains.append(chain)
	
    return chains

  def setResChain(self, null):
    
    chain = self.resChainPulldown.getObject()
     
    if self.residueRange:
      if chain is not self.residueRange.chain:
        residues = chain.sortedResidues()
        self.residueRange.chain = chain
        self.residueRange.startRes = residues[0]
        self.residueRange.endRes = residues[-1]
    
    self.updateResidueRanges()
    
  def getResChain(self, resRange):
    
    chains = self.getChains()
    if resRange and chains:
      names = [c.code for c in chains]
      self.resChainPulldown.setup(names, chains, names.index(resRange.chain.code) )

  def setResStart(self, event):
    
    if self.residueRange:
      val = self.resStartEntry.get()
      chain   = self.residueRange.chain
      residue = chain.findFirstResidue(seqCode=val)
      if residue:
        self.residueRange.startRes = residue
      
    self.updateResidueRanges()
  
  def getResStart(self, resRange):
    
    if resRange:
      self.resStartEntry.set(resRange.startRes.seqCode)

  def setResEnd(self, event):
    
    if self.residueRange:
      val     = self.resEndEntry.get()
      chain   = self.residueRange.chain
      residue = chain.findFirstResidue(seqCode=val)
      if residue:
        self.residueRange.endRes = residue
      
    self.updateResidueRanges()
      
  def getResEnd(self, resRange):
    
    if resRange:
      self.resEndEntry.set(resRange.endRes.seqCode)
    
  def setTolMinPpm(self, event):
    
    obj = self.toleranceMatrix.currentObject
    if obj:
      obj.minTol = self.tolMinPpmEntry.get()

    self.updateTolerances()
  
  def getTolMinPpm(self, datum):
    
    if datum:
      self.tolMinPpmEntry.set(datum.minTol)

  def setTolMultiplier(self, event):
    
    obj = self.toleranceMatrix.currentObject
    if obj:
      obj.multiplier = self.tolMultiplierEntry.get()    
  
    self.updateTolerances()

  def getTolMultiplier(self, datum):
    
    if datum:
      self.tolMultiplierEntry.set(datum.multiplier)

  def setTolMaxPpm(self, event):
      
    obj = self.toleranceMatrix.currentObject
    if obj:
      obj.maxTol = self.tolMaxPpmEntry.get()

    self.updateTolerances()
      
  def getTolMaxPpm(self, datum):
    
    if datum:
      self.tolMaxPpmEntry.set(datum.maxTol)

  def getIsotopes(self):
  
    isotopes = []
    if self.spectrum:
      for isotope in getSpectrumIsotopes(self.spectrum):
        if isotope not in isotopes:
	  isotopes.append(isotope)
    
    return isotopes

  def getDataDimNames(self):
  
    names = []
    if self.spectrum:
      names = ['F%d' % (i+1) for i in range(self.spectrum.numDim)]
    
    return names

  def setShiftDataDim(self, dataDim):
    
    if self.chemShiftRange and (self.chemShiftRange.dataDim is not dataDim):
      isotopes = self.getDataDimIsotopes(dataDim)
      self.chemShiftRange.dataDim = dataDim
      self.chemShiftRange.isotopes = isotopes
      
      isotope   = isotopes[0]
      rangeDict = self.getShiftRangeDict()
      if rangeDict.get(isotope):
        start, end = rangeDict[isotope]
        self.chemShiftRange.start = start
        self.chemShiftRange.end = end

    self.updateChemShiftRanges()  
  
  def getShiftDataDim(self, shiftRange):
    
    if shiftRange:
      texts = self.getDataDimNames()
      dataDim = shiftRange.dataDim
      dataDims = list(dataDim.dataSource.sortedDataDims())
      index = dataDims.index(dataDim)
      
      self.shiftDataDimPulldown.setup(texts,dataDims,index)

  def setShiftStart(self, event):
    
    if self.chemShiftRange:
      self.chemShiftRange.start = self.shiftStartEntry.get()
    
    self.updateChemShiftRanges()  
  
  def getShiftStart(self, shiftRange):
    
    if shiftRange:
      start = shiftRange.start
      self.shiftStartEntry.set(start)

  def setShiftEnd(self, event):
    
    if self.chemShiftRange:
      self.chemShiftRange.end = self.shiftEndEntry.get()
    
    self.updateChemShiftRanges()  
  
  def getShiftEnd(self, shiftRange):
    
    if shiftRange:
      end = shiftRange.end
      self.shiftEndEntry.set(end)

  def getResidueRanges(self):
  
    # sets the defaults and checks for update
    if not self.residueRanges:

      residueRanges = self.getSpectrumResidueRanges()
      
      if residueRanges:
        # stored
        self.residueRanges = residueRanges
      else:
        # default 
        residueRanges = []
        if self.molSystem and self.molSystem.chains:
          for chain in self.molSystem.chains:
            residues = chain.sortedResidues()

            for dim in self.getResidueDimSelection():
              residueRanges.append( [ dim, chain, residues[0], residues[-1]] )
        
        self.residueRanges = [ResidueRowObject(dataDims,chain,startRes,endRes) for (dataDims,chain,startRes,endRes) in residueRanges]
        self.setSpectrumResidueRanges()
	  
    else:
      
      # check existing
      residueRanges = self.residueRanges[:]
      for obj in residueRanges:
        chain = obj.chain
        residues = chain.sortedResidues()
        chain2 = self.molSystem.findFirstChain(code=chain.code)
        if chain2:
          start = obj.startRes
          end = obj.endRes
          obj.chain = chain2
          
          if start.seqCode > residues[-1].seqCode:
            obj.startRes = residues[0]
            
          elif start.seqCode > residues[0].seqCode:
            obj.startRes = start
            
          else:
            obj.startRes = residues[0]
 
          if end.seqCode < obj.startRes.seqCode:
            obj.endRes = residues[-1]
            
          elif end.seqCode < residues[-1].seqCode:
            obj.endRes = end
            
          else:
            obj.endRes = residues[-1]
            
        else:
          self.residueRanges.remove(obj)
      
      
      # store changes  
      self.setSpectrumResidueRanges()
	
    
    return self.residueRanges


  def getSpectrumResidueRanges(self):
  
    specResidueRanges = None
    if self.spectrum:
      specResidueRanges = getSpectrumNoesResidueRanges(self.spectrum)
    
    residueRanges = []
    if specResidueRanges:
      for residueRange in specResidueRanges:
        (dims,chain,start,end) = residueRange
       
        if chain.molSystem is self.molSystem:
          dimsName = ' - '.join([self.getDataDimName(x) for x in dims])
          residueRange[0] = (dimsName, dims)
          residueRanges.append(ResidueRowObject(*tuple(residueRange)))

    return residueRanges

  def setSpectrumResidueRanges(self):
  
    if self.spectrum and self.residueRanges:
      
      specResidueRanges = []
      for obj in self.residueRanges:
        specResidueRange = []
        specResidueRange.append(obj.dataDims[1])
        specResidueRange.extend([obj.chain,obj.startRes,obj.endRes])
        specResidueRanges.append(specResidueRange)
      
      setSpectrumNoeResidueRanges(self.spectrum, specResidueRanges)
     
  def getTolerances(self):
  
    if not self.tolerances:
      # set default
      tolerances = []
    
      if self.spectrum:
        if hasattr(self.spectrum,'shiftMatchTolerances'):
          tolerances = self.spectrum.shiftMatchTolerances
        else:
          for dataDim in self.spectrum.sortedDataDims():
            expDimRefs = dataDim.expDim.expDimRefs
            if expDimRefs:
	      name = self.getDataDimName(dataDim)
	      isotopes =self.getDataDimIsotopes(dataDim)
	      tol = getAnalysisDataDim(dataDim).noeTolerance
              tolerances.append( [ [name,dataDim],tol,1,tol] )
                
	  self.spectrum.shiftMatchTolerances = tolerances
        self.tolerances = [ToleranceRowObject(dataDimInfo,minTol,multiplier,maxTol) for (dataDimInfo,minTol,multiplier,maxTol) in tolerances]

    else:
      for obj in self.tolerances:
        name,dataDim = obj.dataDimInfo
        getAnalysisDataDim(dataDim).noeTolerance = obj.minTol
    
    return self.tolerances
  
  def getShiftRangeDict(self):
  
    rangeDict = {}
    for axisType in self.parent.getAxisTypes():
      if len(axisType.isotopeCodes) == 1:
        rangeDict[axisType.isotopeCodes[0]] = axisType.region
    return rangeDict
  
  def getDataDimIsotopes(self, dataDim):
  
    isotopes = set()
    for expDimRef in dataDim.expDim.expDimRefs:
      isotopes.update(expDimRef.isotopeCodes)
    
    return list(isotopes)
  
  def getShiftRanges(self):

    if not self.chemShiftRanges:
      #default
      chemShiftRanges = []
      if self.spectrum:
        if hasattr(self.spectrum,'shiftMatchShiftRanges'):
          chemShiftRanges = self.spectrum.shiftMatchShiftRanges
        else:
          chemShiftRanges = getSpectrumNoesShiftRanges(self.spectrum)
          
          if not self.chemShiftRanges:
            rangeDict = self.getShiftRangeDict()

            for dataDim in self.spectrum.sortedDataDims():
              isotopes = self.getDataDimIsotopes(dataDim)
              if rangeDict.get(isotopes[0]):
                (start,end) = rangeDict[isotopes[0]]
                if dataDim.expDim.isAcquisition and (isotopes[0] == '1H'):
                  chemShiftRanges.append( [dataDim,isotopes,start,4.85] )
                  chemShiftRanges.append( [dataDim,isotopes, 4.85,end ] )
                else:
                  chemShiftRanges.append( [dataDim,isotopes,start,end] )
 
            self.spectrum.shiftMatchShiftRanges = chemShiftRanges
            setSpectrumNoeShiftRanges(self.spectrum, chemShiftRanges)
        self.chemShiftRanges = [ShiftRowObject(dataDim,isotopes,start,end) for (dataDim,isotopes,start,end) in chemShiftRanges]

    else:
      # update check
      for obj in self.chemShiftRanges:
        if obj.start > obj.end:
          (obj.start,obj.end) = (obj.end,obj.start)
          
      chemShiftRanges = [(obj.dataDim,obj.isotopes,obj.start,obj.end) for obj in self.chemShiftRanges]
      setSpectrumNoeShiftRanges(self.spectrum, chemShiftRanges)
      self.chemShiftRanges = [ShiftRowObject(dataDim,isotopes,start,end) for (dataDim,isotopes,start,end) in chemShiftRanges]
    
    return self.chemShiftRanges

  def updateResidueRanges(self):

    if self.residueRange and ( len(self.residueRanges)>1 ):
      self.resButtons.buttons[1].enable()
    else:
      self.resButtons.buttons[1].disable()

    if self.molSystem and self.spectrum:
      self.resButtons.buttons[0].enable()
    else:
      self.resButtons.buttons[0].disable()
  
    textMatrix = []
    objectList = self.getResidueRanges()
    for obj in objectList:
      datum = []
      (dimsName,dimsObjects) = obj.dataDims
      datum.append(dimsName)
      datum.append(obj.chain.code)
      datum.append(obj.startRes.seqCode)
      datum.append(obj.endRes.seqCode)
      textMatrix.append(datum)
    
    self.rangeMatrix.update(objectList=objectList,textMatrix=textMatrix)
    
  def updateTolerances(self):
  
    textMatrix = []
    objectList = self.getTolerances()
    for obj in objectList:
      (dataDimName,dataDim) = obj.dataDimInfo
      datum = []
      datum.append(dataDimName)
      datum.append(obj.minTol)
      datum.append(obj.multiplier)
      datum.append(obj.maxTol)
      textMatrix.append(datum)
    
    self.toleranceMatrix.update(objectList=objectList,textMatrix=textMatrix)
     
  def updateChemShiftRanges(self):
  
    if self.chemShiftRange:
      # fore the one per isotope check
      self.selectChemShiftRange(self.chemShiftRange, 0, 0)
    else:
      self.shiftsButtons.buttons[1].disable()
    
    textMatrix = []
    objectList = self.getShiftRanges()
    for obj in objectList:
      datum = []
      
      dataDimStr = 'F%d' % obj.dataDim.dim
      isotopeStr = ' '.join(obj.isotopes)
      datum.append(dataDimStr)
      datum.append(isotopeStr)
      datum.append(obj.start)
      datum.append(obj.end)
      textMatrix.append(datum)
    
    self.shiftsMatrix.update(objectList=objectList,textMatrix=textMatrix)
 
  def setMinPeakMerit(self, *opt):
  
    self.minPeakMerit = self.minMeritEntry.get()
 

  def changeLabellingScheme(self, scheme):
  
    self.labellingScheme = scheme
    

  def changeConstraintSet(self, constraintSet):
 
    self.constraintSet = constraintSet
  
  def changeStructure(self, structure):
  
    self.structure = structure
  
  def changeMolSystem(self, molSystem):
  
    self.molSystem = molSystem
    self.residueRanges = []
    self.updateResidueRanges()
  
  def changePeakList(self, peakList):
  
    if peakList is not self.peakList:
      self.peakList = peakList
      self.noeFrame.update(self.peakList)
      self.updatePeakFrame()
      
      spectrum = peakList.dataSource
      if spectrum is not self.spectrum:
        self.spectrum = spectrum
        self.residueRanges = []
        self.updateMolSystems()
        self.updateResidueRanges()
        self.updatePeakLists()
        self.tolerances = []
        self.updateTolerances()
        self.chemShiftRanges = []
        self.updateChemShiftRanges()
    
  def getMolSystems(self):
  
    molSystems = []
    
    if self.spectrum:
      for molSystem in self.spectrum.experiment.molSystems:
        if molSystem.chains:
          molSystems.append(molSystem)
        
    if not molSystems:
      for molSystem in self.project.molSystems:
        if molSystem.chains:
          molSystems.append(molSystem)
           
    return molSystems
   
  def addResidueRange(self):
  
    if self.molSystem and self.molSystem.chains:
      chain = self.molSystem.sortedChains()[0]
      residues = chain.sortedResidues()
      startRes = residues[0]
      endRes   = residues[-1]
    
      dims = self.getResidueDimSelection()
    
      residueRange = ResidueRowObject(dims[0],chain,startRes,endRes)
 
      self.residueRanges.append(residueRange)
 
      self.updateResidueRanges()

  def deleteResidueRange(self, *event):
  
    if self.residueRange:
      self.residueRanges.remove(self.residueRange)
      self.residueRange = []
      self.updateResidueRanges()
  
  def addChemShiftRange(self):
  
    if self.spectrum:
      rangeDict = self.getShiftRangeDict()
      dataDim   = self.spectrum.findFirstDataDim(dim=1)
      isotopes  = self.getDataDimIsotopes(dataDim)
      isotope   = isotopes[0]
      if not rangeDict.get(isotope):
        isotope = rangeDict.keys()[0]
 
      start,end = rangeDict[isotope]
      chemShiftRange= ShiftRowObject(dataDim,isotopes,start,end)
 
      self.chemShiftRanges.append(chemShiftRange)
 
      self.updateChemShiftRanges()
  
  def deleteChemShiftRange(self, *event):
  
    if self.chemShiftRange:
      dataDim = self.chemShiftRange.dataDim
      dataDims = [x.dataDim for x in self.chemShiftRanges]
      if dataDims.count(dataDim) < 2:
        msg = 'Cannot delete range. There must be at least one for each dimension'
        showWarning('Warning',msg, parent=self)
	return
      
      self.chemShiftRanges.remove(self.chemShiftRange)
      self.chemShiftRange = []
      self.updateChemShiftRanges()

  def getDataDimName(self,dataDim):
  
    isotopes = ''
    if dataDim.expDim.expDimRefs:
      for isotope in self.getDataDimIsotopes(dataDim):
        isotopes += isotope
  
    name = 'F%d %s' % (dataDim.dim, isotopes)
    return name
  
  def setDistanceClasses(self):
  
    self.guiParent.editNoeClasses()
  
  def testShiftMatch(self):
 
    if self.spectrum and self.molSystem:
      if not self.peakList.peaks:
        showWarning('Failure','Selected peak list has no peaks', parent=self)
        return
        
      intensityType = self.noeFrame.intensityType
      if not self.checkPeakListIntensities(self.peakList, intensityType):
        return
    
      spectrum = self.spectrum
      residueRanges = []
      for residueRange in self.residueRanges:
        datum = []
        datum.append( residueRange.dataDims[1] ) #dataDims
        datum.append( residueRange.chain ) # chain 
        datum.append( residueRange.startRes ) # start residue
        datum.append( residueRange.endRes ) # end residue
        
        residueRanges.append( datum )

      tolerances = []
      for obj in self.tolerances:
        (name, dataDim) = obj.dataDimInfo
        tolerances.append( [dataDim, obj.minTol, obj.maxTol, obj.multiplier] )
	
      self.setMinPeakMerit()
      minMerit        = self.minPeakMerit
      scale           = self.noeFrame.getScale()
      params          = self.noeFrame.getParams()
      minLabelFrac    = self.minIsoFractionEntry.get() or 0.0
      chemShiftRanges = self.getShiftRanges()
      maxDist = self.maxDistEntry.get() or None
      peakCategories  = {}
      progressBar     = ProgressBar(self, text="Matching resonances to peaks",total = len(self.peakList.peaks))
      
      chemShiftRanges = [(obj.dataDim, obj.isotopes, obj.start, obj.end) for obj in chemShiftRanges]
      constraintList  = makeAmbigDistConstraints(self.peakList, tolerances, chemShiftRanges, testOnly=True,
                                                 labelling=self.labellingScheme, minLabelFraction=minLabelFrac,
                                                 constraintSet=self.constraintSet,residueRanges=residueRanges,
                                                 minMerit=minMerit, doAliasing=self.aliasing, progressBar=progressBar, 
                                                 ignoreDiagonals=self.ignoreDiagonals, intensityType=intensityType,
                                                 structure=self.structure, maxDist=maxDist,
                                                 scale=scale, params=params, peakCategories=peakCategories)
 
      progressBar.destroy()
      
      total = 0
      msg = ''
      groupNames = peakCategories.keys()
      groupNames.sort()
      groups = []
      tipTexts = [GROUP_TIP_TEXTS['Matched'],]
      for key in groupNames :
        peaks = peakCategories[key]
        groups.append(peaks)
        tipTexts.append(GROUP_TIP_TEXTS.get(key, ''))
        n = len(peaks)
        msg += '%d peaks exclusions categorised as %s\n' % (n, key)
        total += n
      
      msg = 'Total of %d peaks matched\n' % (len(self.peakList.peaks)-total) + msg
      msg = 'Total of %d peaks excluded\n' % (total) + msg
      title  = 'Shift Match Report'
      
      peakDict = {}
      for group in groups:
        for peak in group:
          peakDict[peak] = True
        
      matchedPeaks = []
      for peak in self.peakList.peaks:
        if not peakDict.get(peak):
          matchedPeaks.append(peak)
      
      groups = [matchedPeaks,] + groups
      groupNames = ['Matched',] + groupNames
      
      popup = ViewPeakGroupsPopup(self.guiParent, title, msg, groups, groupNames, tipTexts)

  def checkPeakListIntensities(self, peakList, intensityType):
  
    peaks = peakList.peaks
    N = len(peaks)
    
    c = 0
    for peak in peaks:
      if not peak.findFirstPeakIntensity(intensityType=intensityType):
        c +=1
        
    if c == 0:
      return 1
      
    return showOkCancel('Warning','%d peaks from %d in the peak list are missing a %s value. Continue?' % (c,N,intensityType), parent=self)

  def calcShiftMatchRestraints(self):
  
    self.calculateConstraints(doShiftMatch=True)
  
  def calcNormRestraints(self):
  
    self.calculateConstraints(doShiftMatch=False)

          
  def calculateConstraints(self, doShiftMatch=False):
  
    if self.spectrum and self.molSystem and self.peakList:
      if not self.peakList.peaks:
        showWarning('Failure','Selected peak list has no peaks', parent=self)
        return

      intensityType = self.noeFrame.intensityType
      if not self.checkPeakListIntensities(self.peakList, intensityType):
        return
        
      constraintsPopup = self.parent.popups.get('browse_restraints')
      if constraintsPopup:
        constraintsPopup.turnOffNotifiers()

      spectrum = self.spectrum
      residueRanges = []
      for residueRange in self.residueRanges:
        datum = []
        datum.append( residueRange.dataDims[1] ) #dataDims
        datum.append( residueRange.chain ) # chain
        datum.append( residueRange.startRes ) # start residue
        datum.append( residueRange.endRes ) # end residue
        
        residueRanges.append( datum )
	
      self.setMinPeakMerit()
      minMerit = self.minPeakMerit
      scale    = self.noeFrame.getScale()
      params   = self.noeFrame.getParams()
      maxDist = self.maxDistEntry.get() or None

      if doShiftMatch:

        minLabelFrac = self.minIsoFractionEntry.get() or 0.0
        tolerances = []
        for obj in self.tolerances:
          (name, dataDim) = obj.dataDimInfo
          tolerances.append( [dataDim,obj.minTol, obj.maxTol, obj.multiplier] )
 
        progressBar = ProgressBar(self, text="Matching resonances to peaks",
                                  total=len(self.peakList.peaks))
                                  
        chemShiftRanges = [(obj.dataDim, obj.isotopes, obj.start, obj.end) for obj in self.chemShiftRanges]
        constraintList = makeAmbigDistConstraints(self.peakList, tolerances, chemShiftRanges,
                                                  constraintSet=self.constraintSet,
                                                  labelling=self.labellingScheme,
                                                  minLabelFraction=minLabelFrac,                                                  
                                                  progressBar=progressBar,
                                                  residueRanges=residueRanges,
                                                  minMerit=minMerit, doAliasing=self.aliasing,
                                                  scale=scale, params=params,
                                                  ignoreDiagonals=self.ignoreDiagonals,
                                                  intensityType=intensityType,
                                                  structure=self.structure,
                                                  maxDist=maxDist)
 
        progressBar.destroy()
      
      else:
        if intensityType == 'height': 
          getIntensity = getPeakHeight
        else:
          getIntensity = getPeakVolume
        scaleDict = {}
        for peak in self.peakList.peaks:
          whichDistance = self.whichDistanceMap.get(peak)
          if whichDistance == 'Diagonal':
            diagPeak = self.diagPeakMap.get(peak)
            diagIntensity = getIntensity(diagPeak)
            scaleDict[peak] = diagIntensity
          elif whichDistance == 'HSQC':
            hsqcPeak = self.hsqcPeakMap.get(peak)
            hsqcIntensity = getIntensity(hsqcPeak)
            scaleDict[peak] = hsqcIntensity
          elif whichDistance == 'None':
            scaleDict[peak] = None
 
        constraintList = makeDistConstraints(self.peakList,
                                             constraintSet=self.constraintSet,
                                             labelling=self.labellingScheme,
                                             intensityType=intensityType,
                                             residueRanges=residueRanges,
                                             minMerit=minMerit,
                                             params=params,
                                             scale=scale, scaleDict=scaleDict)
      
      self.updateConstraintSets()
      
      if constraintList: # Did not barf
        self.constraintSet = constraintList.nmrConstraintStore

        if constraintsPopup:
          constraintsPopup.update(constraintList)
          constraintsPopup.turnOnNotifiers()
          constraintsPopup.open()
        else:
          self.guiParent.browseConstraints(constraintList)
  
  def selectResidueRange(self, residueRange, row, col):
  
    if residueRange:
      self.residueRange = residueRange
      if len(self.residueRanges) < 2:
        self.resButtons.buttons[1].disable()
      else:
        self.resButtons.buttons[1].enable()

  def selectChemShiftRange(self, chemShiftRange, row, col):
  
    if chemShiftRange:
      self.chemShiftRange = chemShiftRange
      dataDim  = chemShiftRange.dataDim
      dataDims = [x.dataDim for x in self.chemShiftRanges]
      if dataDims.count(dataDim) < 2:
        self.shiftsButtons.buttons[1].disable()
      else:
        self.shiftsButtons.buttons[1].enable()

  def updateLabellingSchemes(self, obj=None):
  
    index = 0
    schemes = [True, None,] + self.project.sortedLabelingSchemes()
    scheme = self.labellingScheme
    names = ['Automatic from sample', '<None>',] + [sc.name for sc in schemes[2:]]
  
    if scheme not in schemes:
      scheme = schemes[0]
    
    index = schemes.index(scheme)
      
    if scheme is not self.labellingScheme:
      self.labellingScheme = scheme
      self.updateAfter()
  
    self.labellingSchemePulldown.setup(names, schemes, index)    

  def updateConstraintSets(self, *opt):
  
    index = 0
    names = ['<New>',]
    constraintSets = [None,]
    
    constraintSets += self.nmrProject.sortedNmrConstraintStores()
    if constraintSets:
      if self.constraintSet not in constraintSets[1:]:
        self.constraintSet = constraintSets[0]
        
      index = constraintSets.index(self.constraintSet)
      names += ['%d' % cs.serial for cs in constraintSets[1:]]
  
    self.constraintSetPulldown.setup(names, constraintSets, index)

  def updateStructures(self, *opt):
  
    index = 0
    names = ['<None>',]
    ensembles = [None,]
    
    ensembles += self.project.sortedStructureEnsembles()
    if self.structure not in ensembles[1:]:
      self.structure = ensembles[0]
       
    index = ensembles.index(self.structure)
    names += ['%s:%d' % (s.molSystem.code, s.ensembleId) for s in ensembles[1:]]
  
    self.structurePulldown.setup(names, ensembles, index)
  
  def updateMolSystems(self, *opt):
  
    names = []
    index = -1
    molSystems = self.getMolSystems()

    if molSystems:
      if self.molSystem not in molSystems:
        self.molSystem = molSystems[0]
    
      names = [ms.code for ms in molSystems]
      index = molSystems.index(self.molSystem)
    
    self.molSysPulldown.setup(names, molSystems, index)
  
  def updateResidue(self, residue):
  
    self.updateResidueRanges()
    
  def updateChain(self, chain):
  
    if self.molSystem and (chain.molSystem is self.molSystem):
      residueRanges = self.residueRanges[:]
      for residueRange in residueRanges:
        if residueRange.chain.isDeleted:
          self.residueRanges.remove(residueRange)
    else:
      self.updateMolSystems()
    
  def getSpectrumName(self, spectrum):
  
    name = '%s:%s' % (spectrum.experiment.name,spectrum.name)
    return name
  
  def updatePeakLists(self, *opt):
  
    index = 0
    names = []
    excludeSimulated = not self.includeSimulatedSelect.get()
    peakLists = getThroughSpacePeakLists(self.project, excludeSimulated)
    peakList = self.peakList
    getName = self.getSpectrumName
      
    if peakLists:
      names = ['%s:%s' % (getName(pl.dataSource), str(pl.serial)) for pl in peakLists]
 
      if peakList not in peakLists:
        peakList = peakLists[0]
 
      index = peakLists.index(peakList)
      spectrum = peakList.dataSource
         
    else:
      peakList = None
      spectrum = None
    
    if peakList is not self.peakList:
      self.peakList = peakList
      self.noeFrame.update(peakList)
    
    if spectrum is not self.spectrum:
      self.spectrum = spectrum
      self.residueRanges = []
      self.updateMolSystems()
      self.updateResidueRanges()
      self.updatePeakLists()
      self.tolerances = []
      self.updateTolerances()
      self.chemShiftRanges = []
      self.updateChemShiftRanges()
     
    self.peakListPulldown.setup(names, peakLists, index)
    
    index = 0
    spectra = []
    names = []
    
    if spectrum:
      availSpectra = []
      isotopes = getSpectrumIsotopes(spectrum)
      spectra = set([pl.dataSource for pl in peakLists])
      for spec in spectra:
        if spec is spectrum:
          continue
        if getSpectrumIsotopes(spec) != isotopes:
          continue
          
        availSpectra.append( ('%s:%s' % (spec.experiment.name,spec.name), spec) )
 
      availSpectra.sort()
      names = [x[0] for x in availSpectra]
      spectra = [x[1] for x in availSpectra]
    
    self.copySpecPulldown.setup(names, spectra, index)

    if spectra:
      self.copySpecButton.enable()
    else:
      self.copySpecButton.disable()
  
  def copySpecTols(self):
  
    if self.peakList:
      spectrum = self.copySpecPulldown.getObject()
      
      if spectrum:
        if hasattr(spectrum,'shiftMatchTolerances'):
          tolerances = spectrum.shiftMatchTolerances
        else:
          tolerances = []
          for dataDim in spectrum.sortedDataDims():
            expDimRefs = dataDim.expDim.expDimRefs
            if expDimRefs:
	      name = self.getDataDimName(dataDim)
	      isotopes = self.getDataDimIsotopes(dataDim)
	      tol = getAnalysisDataDim(dataDim).noeTolerance
              tolerances.append( [ [name,dataDim],tol,1,tol] )
        
        tolCopy = []
        for dataDimInfo, minTol, multi, maxTol in tolerances:
          nameOld, dataDimOld = dataDimInfo
          dataDim = self.spectrum.findFirstDataDim(dim=dataDimOld.dim)
          name = self.getDataDimName(dataDim)
          tolCopy.append( ToleranceRowObject([name,dataDim], minTol, multi, maxTol))
               
        self.spectrum.shiftMatchTolerances = tolerances
        self.tolerances = tolCopy
        self.updateTolerances()
      
  def runAnchoring(self):  
    
    TOLERANCE_DICT['1H']  = self.tolHydrogenEntry.get() or 0.04
    TOLERANCE_DICT['15N'] = self.tolNitrogenEntry.get() or 0.4
    TOLERANCE_DICT['13C'] = self.tolCarbonEntry.get() or 0.4
    
    if not self.networkPeakListsMatrix.objectList:
      msg = 'No peaks lists selected for network anchoring'
      showWarning('Failure',msg,parent=self)
      self.tabbedFrame.select(4)
      return   
    
    peakLists = []
    numPeaks  = 0
    for peakList in self.networkPeakListsMatrix.objectList:
      if peakList.networkAnchorUse:
        peakLists.append(peakList)
        numPeaks += len(peakList.peaks)
        
    if peakLists:
      if not numPeaks:
        showWarning('Failure','No peaks in selected peak lists',parent=self)
        return
      
      constraintsPopup = None  
      if self.constraintSet:
        constraintsPopup = self.parent.popups.get('browse_restraints')
        if constraintsPopup:
          constraintsPopup.turnOffNotifiers()
    
      assignPeakList = self.peakListCheck.get()
      intensityType  = self.noeFrame.intensityType
      params         = self.noeFrame.getParams()
      strictness     = self.strictnessPulldown.index
      minLabelFrac   = self.minIsoFractionEntry.get() or 0.0
      progressBar    = ProgressBar(self, text="Network anchoring",total=numPeaks)
      threshold      = self.thresholdEntry.get() or 1.0
      constraintList = networkAnchorAssign(peakLists, intensityType=intensityType,
                                           strictness=strictness,progressBar=progressBar,
                                           constraintSet=self.constraintSet, distParams=params,
                                           labelling=self.labellingScheme,
                                           minLabelFraction=minLabelFrac, 
                                           threshold=threshold, isotopeTolerances=TOLERANCE_DICT, 
                                           assignPeakList=assignPeakList,
                                           structure=None, nexus=None)

      progressBar.destroy()

      if constraintsPopup:
        if constraintList:
          constraintsPopup.update(constraintList)
        constraintsPopup.turnOnNotifiers()

      if constraintList:
        self.parent.browseConstraints(constraintList)

    else:
        showWarning('Failure','No peak lists selected',parent=self)
    
    
  def toggleIntensityChoice(self, entry):

    self.updatePeakFrameMatrix() 
    
  def changeHsqcPeakList(self, peakList):

    self.updatePeakFrameMatrix() 
  
  def setWhichUse(self, junk):

    peak = self.peaksMatrix.currentObject
    self.whichDistanceMap[peak] = self.whichUsePulldown.getText()
    self.updatePeakFrameMatrix() 

  def getWhichUse(self, peak):

    texts = ['None']
    if self.findPeakDistance(peak, 1.0):
      ###intensityType = self.intensityButtons.get().lower()
      intensityType = self.noeFrame.intensityType
      if peak.findFirstPeakIntensity(intensityType=intensityType):
        texts.append('Reference')
      if self.diagPeakMap.get(peak):
        texts.append('Diagonal')
      if self.hsqcPeakMap.get(peak):
        texts.append('HSQC')

    objects = None
    ind = 0

    self.whichUsePulldown.setup(texts, objects, ind)

  def setWhichDistance(self):

    peaks = self.peaksMatrix.currentObjects
    if not peaks:
      msg = 'No peaks selected in table'
      showWarning('Warning',msg, parent=self)
      return

    objects = texts = ['Reference', 'Diagonal', 'HSQC', 'Cancel']
    whichDistance = showMulti('Distance', 'Select which distance to use', texts, objects, parent=self)

    if whichDistance == 'Cancel':
      return
    
    for peak in peaks:
      if whichDistance == 'Diagonal':
        diagPeak = self.diagPeakMap.get(peak)
        if not diagPeak:
          continue
      elif whichDistance == 'HSQC':
        hsqcPeak = self.hsqcPeakMap.get(peak)
        if not hsqcPeak:
          continue
      self.whichDistanceMap[peak] = whichDistance
 
    self.updatePeakFrameMatrix()
 
  def updatePeakFrame(self):

    self.updateHsqcPeakLists()
    self.updatePeakFrameMatrix()

  def updateHsqcPeakLists(self, *opt):

    ind = 0
    peakList = self.peakList
    if peakList:
      peakLists = getCompatibleHsqcPeakLists(peakList)
      getName = self.getSpectrumName
      names = ['%s:%s' % (getName(pl.dataSource), str(pl.serial)) for pl in peakLists]
    else:
      peakLists = []
      names = []

    self.hsqcPeakListPulldown.setup(names, peakLists, ind)

  def calcDiagPeakMap(self, peakList):

    # this calculates a map from a peak to a corresponding diagonal peak

    dataDimPairs = getDiagonalBondedDataDims(peakList.dataSource)
    peaks = peakList.peaks

    self.diagPeakMap = {}
    for dataDim1, dataDim2 in dataDimPairs:
      dim1 = dataDim1.dim
      dim2 = dataDim2.dim
      dd = {}
      for peak in peaks:
        peakDim1 = peak.findFirstPeakDim(dim=dim1)
        peakDim2 = peak.findFirstPeakDim(dim=dim2)
        resonances1 = getPeakDimResonances(peakDim1)
        resonances2 = getPeakDimResonances(peakDim2)
        resonances = resonances1 & resonances2
        for resonance in resonances:
          dd[resonance] = peak
      for peak in peaks:
        peakDim1 = peak.findFirstPeakDim(dim=dim1)
        resonances1 = getPeakDimResonances(peakDim1)
        for resonance in resonances1:
          if resonance in dd:
            # just take first peak that works
            # TBD: if more than one peak that works
            # then probably should do something smarter
            # could happen in 3D in weird situation
            # could happen in 4D quite easily
            self.diagPeakMap[peak] = dd[resonance]

  def calcHsqcPeakMap(self, hsqcPeakList, peakList):

    # this calculates a map from a peakList peak to an hsqcPeakList peak
    # with the same resonances

    # this might not be the best/fastest way to get at the mapping
    # also ignoring which dims are onebond because it ought to work
    # ok but not guaranteed I guess
    numDim = hsqcPeakList.dataSource.numDim # should be 2 presumably
    for peak in hsqcPeakList.peaks:
      peakCountDict = {}
      for peakDim in peak.peakDims:
        peakSet = set()
        for peakDimContrib in peakDim.peakDimContribs:
          resonance = peakDimContrib.resonance
          for peakDimContrib2 in resonance.peakDimContribs:
            peak2 = peakDimContrib2.peakDim.peak
            if peak2.peakList is peakList:
              peakSet.add(peak2)
        for peak2 in peakSet:
          peakCountDict[peak2] = peakCountDict.get(peak2, 0) + 1
      for peak2 in peakCountDict:
        if peakCountDict[peak2] == numDim:
          self.hsqcPeakMap[peak2] = peak

  def findPeakDistance(self, peak, intensityScale):

    if intensityScale is None:
      return None

    ###intensityType = self.intensityButtons.get().lower()
    intensityType = self.noeFrame.intensityType
    intensity = peak.findFirstPeakIntensity(intensityType=intensityType)
    if not intensity:
      return None

    intensityValue = abs(intensity.value)
    intensityScale = abs(intensityScale)

    params   = self.noeFrame.getParams()
    distanceFunction = lambda val:getNoeDistance(val, params)
    resonances0, resonances1 = findPeakConstraintPairs(peak)
    if not resonances0 or not resonances1:
      return None

    (dist,minDist,maxDist) = getDistMinMax(intensityValue, intensityScale, resonances0, resonances1, distanceFunction)

    return dist

  def updatePeakFrameMatrix(self, *extra):

    self.diagPeakMap = {}
    self.hsqcPeakMap = {}
    peakList = self.peakList
    if peakList:
      peaks = peakList.sortedPeaks()
      self.calcDiagPeakMap(peakList)
      hsqcPeakList = self.hsqcPeakListPulldown.getObject()
      if hsqcPeakList:
        self.calcHsqcPeakMap(hsqcPeakList, peakList)
    else:
      peaks = []
    
    intensityType = self.noeFrame.intensityType
    if intensityType == 'height': 
      getIntensity = getPeakHeight
    else:
      getIntensity = getPeakVolume

    textMatrix = []
    colorMatrix = []
    distanceColIndex = {'Reference': 3, 'Diagonal': 4, 'HSQC': 5}
    for peak in peaks:
      texts = []

      peakIntensity = getIntensity(peak)
      intensityScale = self.noeFrame.getScale() or 1.0
      peakDist = self.findPeakDistance(peak, intensityScale)

      diagPeak = self.diagPeakMap.get(peak)
      if diagPeak:
        diagText = getSimplePeakAnnotation(diagPeak)
        diagSerial = diagPeak.serial
        diagIntensity = getIntensity(diagPeak)
        diagDist = self.findPeakDistance(peak, diagIntensity)
      else:
        diagText = None
        diagSerial = None
        diagIntensity = None
        diagDist = None

      hsqcPeak = self.hsqcPeakMap.get(peak)
      if hsqcPeak:
        hsqcText = getSimplePeakAnnotation(hsqcPeak)
        hsqcSerial = hsqcPeak.serial
        hsqcIntensity = getIntensity(hsqcPeak)
        hsqcDist = self.findPeakDistance(peak, hsqcIntensity)
      else:
        hsqcText = None
        hsqcSerial = None
        hsqcIntensity = None
        hsqcDist = None

      whichDistance = self.whichDistanceMap.get(peak)
      if whichDistance == 'Reference':
        if peakDist is None:
          whichDistance = None
      elif whichDistance == 'Diagonal':
        if diagDist is None:
          whichDistance = None
      elif whichDistance == 'HSQC':
        if hsqcDist is None:
          whichDistance = None

      if not whichDistance:
        if peakDist is not None:
          whichDistance = 'Reference'
        elif diagDist is not None:
          whichDistance = 'Diagonal'
        elif hsqcDist is not None:
          whichDistance = 'HSQC'
        else:
          whichDistance = 'None'
      texts = [peak.serial, getSimplePeakAnnotation(peak),
               whichDistance,
               peakDist, diagDist, hsqcDist,
               peakIntensity,
               diagSerial, diagText, diagIntensity,
               hsqcSerial, hsqcText, hsqcIntensity]
      textMatrix.append(texts)

      colors = len(texts) * [None]
      if whichDistance != 'None':
        colors[distanceColIndex[whichDistance]] = '#C0FFC0'
      colorMatrix.append(colors)

    self.peaksMatrix.update(textMatrix=textMatrix,
                            objectList=peaks,
                            colorMatrix=colorMatrix)
   
  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)


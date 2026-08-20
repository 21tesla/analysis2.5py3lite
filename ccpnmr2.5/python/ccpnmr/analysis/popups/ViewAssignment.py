
"""
======================COPYRIGHT/LICENSE START==========================

ViewAssignment.py: Part of the CcpNmr Analysis program

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

from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.Canvas import Canvas
from memops.gui.CheckButton import CheckButton
from memops.gui.Color import black, hsbToRgb, inverseGrey, standardColors
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.PartitionedSelector import PartitionedSelector
from memops.gui.PulldownList import PulldownList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.RadioButtons import RadioButtons
from memops.gui.ScrolledCanvas import ScrolledCanvas
from memops.gui.Spacer import Spacer
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.ToggleLabel import ToggleLabel

from ccpnmr.analysis.core.AssignmentBasic import isAtomAssigned
from ccpnmr.analysis.core.ExperimentBasic  import getOnebondDataDims
from ccpnmr.analysis.core.MoleculeBasic import STANDARD_ISOTOPES, unicodeGreek
from ccpnmr.analysis.core.MoleculeBasic import areAtomsBound
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.MoleculeBasic import getResidueCode
from ccpnmr.analysis.core.Util import stringFromExperimentSpectrum, getSpectrumPosContourColor

from ccp.util.LabeledMolecule import getIsotopomerSingleAtomFractions, singleAtomFractions

HILIGHT_COLOR = '#%02x%02x%02x' % (128, 128, 200)
TITLE_COLOR   = '#%02x%02x%02x' % (0, 0, 128)

H_PHOBE= {'Phe':-3.7,'Met':-3.4,'Ile':-3.1,'Leu':-2.8,
          'Val':-2.6,'Cys':-2.0,'Trp':-1.9,'Ala':-1.6,
          'Thr':-1.2,'Gly':-1.0,'Ser':-0.6,'Pro':0.2 ,
          'Tyr':0.7 ,'His':3.0 ,'Gln':4.1 ,'Asn':4.8 ,
          'Glu':8.2 ,'Lys':8.8 ,'Asp':9.2 ,'Arg':12.3}
          
RES_CHARGE= {'Phe':0, 'Met':0,'Ile':0, 'Leu':0,
             'Val':0, 'Cys':0,'Trp':0, 'Ala':0,
             'The':0, 'Gly':0,'Ser':0, 'Pro':0,
             'Tyr':0, 'His':1,'Gln':0, 'Asn':0,
             'Glu':-1,'Lys':1,'Asp':-1,'Arg':1}
             
class ViewAssignmentPopup(BasePopup):
  """
  **Graphical Display of Assignments and Connections**
  
  This popup window is designed to give a graphical display of the atom sites
  available in a given molecular chain. In this regard an atom site is a
  position along the side chain or scaffold of a residue, like the "beta"
  position in  amino acids, or the "three prime" position in nucleic acids.
  Accordingly, each atom site represents a carbon or nitrogen atom *and* any
  bound hydrogens. Representing atom sites rather than individual atoms provides
  an alternative to the `Atom Browser`_ with less clutter but still with
  sufficient information to give a clear indication of assignment status, and 
  in this particular case, of assignment connectivity between atom sites

  The layout involves two tabs, one on the left contains the actual graphical
  display of the atom sites with assignments, the right "Options" tab  specifies
  the details of how to display things; spectrum connections, atom sites and
  residues to include.

  In normal operation the user starts by slecting the appropriate molecular
  chain in the first tab; labelled symbols will appear to display the available
  atom sites. The user then switches to the "Options" tab and chooses which
  spectra, atom sites and residue types to display. Moving back to the first tab
  the user will see that the symbols of assigned atom sites are circled, with a
  colour appropriate to the spectrum they are assigned within. Also, if the
  selected spectra contain assignments that link different atom sites the
  graphical display will include lines that connect the atom site symbols. The
  user may use the "Connection Type" pulldown menu to display different kinds of
  connectivity lines. By default the displyed connectivity is "Within 4
  Residues". Having more than one set of spectrum massignments and connections
  displayed at the same time is possible, but this may obscure information in
  the chart.

  The background colour of the atom set symbol (a labelled circle) is normally
  determined according to whether there are *any resonance assignments* to that
  site, considering all spectra; black indicates assigned and grey unassigned. 
  If the user chooses an isotope labelling scheme/pattern in the "Options" tab,
  then this black/grey indicator is replaced by colours which represent the
  amount of spin-active isotope incorporation at the atom site, according the
  scheme or molecule patterns. The various options in the "Isotope Labelling"
  section control which colours to use, whether to use only two colours or a
  gradient, which isotope is considered (e.g. usually 13C or 1H) and what the
  minimum incorporation is for something to be considered as labelled.

  **Caveats & Tips**

  Clicking with the left mouse button on any of the atom site symbols for a
  residue will toggle the activation for the whole residue; flipping between the
  greyed-out incative state, that shows no assignments or connections, and the
  fully active state.

  Right click with the mouse over the main chart to save the graphical display
  as a PostScript image, which maye be printed or manipulated in an external
  program.

  The colours of the spectrum toggle buttons in the "Options" tab is determined
  by the positive contour colour for the relevant spectrum. If the spectrum uses
  a multicolour scheme the colour used is one from about two thirds of the way
  through the scheme.

  .. _`Atom Browser`: BrowseAtomsPopup.html

  """
  
  def __init__(self, parent, *args, **kw):
  
    self.font = 'Helvetica 10'
    self.sFont = 'Helvetica %d'
    self.project   = parent.project
    self.guiParent = parent
    self.waiting   = False
    
    BasePopup.__init__(self, parent, title= "Chart : Assignment Graph", **kw)

  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)

  def body(self, guiFrame):

    self.geometry('600x600')

    self.chain         = None
    self.spectraDict   = {}
    self.spectra       = []
    self.spectraNames  = []
    self.canvasDict    = {}
    self.guiNameDict   = {}
    self.atomNamesDict = {}
    self.atomNamesList = []
    self.residuesDict  = {}
    self.bindDict = {}
    self.residueColorMode = 'None'
    self.maxConnectionDist = 4
    self.minConnectionDist = 0
    self.atomSelectionColumns = 10
    self.specSelectionColumns = 5
    self.labellingScheme = None
    
    isotopes = self.getIsotopes()
    self.isotope = isotopes[1]

    guiFrame.expandGrid(0,0)

    tipTexts = ['Shows a chart illustrating selected atom groups for a given sequence, their assignment status and connections between them',
                'Allows specification of which atoms, residues and spectrum information to use']
    tabbedFrame = TabbedFrame(guiFrame,options=['Assignment & Connections Chart','Options'],
                              callback=self.toggleTab, grid=(0,0), tipTexts=tipTexts)
    self.tabbedFrame = tabbedFrame

    graphFrame, optFrame = tabbedFrame.frames
    graphFrame.expandGrid(1,0)
    optFrame.expandGrid(5,0)

    # Graph

    topFrame = Frame(graphFrame, grid=(0,0))
    topFrame.grid_columnconfigure(3, weight=1)
    
    label = Label(topFrame, text='Chain:', grid=(0,0))
    tipText = 'Selects which molecular chain to display the residue and atom groups for'
    self.molPulldown = PulldownList(topFrame, callback=self.changeMolecule,
                                    grid=(0,1), tipText=tipText)
    

    label = Label(topFrame, text=' Connection Type:', grid=(0,2))
    texts= ['Intra-residue','Sequential','Within 4 Residues',
            '2 to 4 Residues','Long Range','All Connections' ]
    tipText = 'Restricts the assignment connection lines to atom groups with particular sequential relationships'
    self.modePulldown = PulldownList(topFrame, callback=self.showSetConnections,
                                     texts=texts, objects=range(len(texts)), index=2,
                                     grid=(0,3), tipText=tipText)
 
    self.canvasFrame = ScrolledCanvas(graphFrame, relief='groove', borderwidth=2,
                                      resizeCallback=None, height=200, grid=(1,0))
    
    self.canvas = self.canvasFrame.canvas
    self.canvas.bind('<Button-1>', self.toggleResidue)
    
    
    # Options   
    
    row = 0
    self.specOptFrame = LabelFrame(optFrame, text='Spectra', grid=(row,0))
    self.specOptFrame.grid_columnconfigure(1, weight=1)
    
    label= Label(self.specOptFrame, text = 'Include predicted peak lists:', grid=(0,0))
    tipText = 'Whether to include predicted peak lists in graph'
    self.includePredictedCheck = CheckButton(self.specOptFrame, selected=False, grid=(0,1), tipText=tipText)

    tipText = 'Toggles individual spectra on and off in tha chart; when set on a spectrum will be used to highlight atom groups assigned and connected in that spectrum'
    self.specSelector=PartitionedSelector(self.specOptFrame, self.toggleSpectrum,
                                          maxRowObjects=5, grid=(1,0), tipText=tipText, gridSpan=(1,2))
    
    row += 1
    self.specFrame = Frame(optFrame, relief='flat', borderwidth=0, grid=(row,0))
    self.specFrame.grid_columnconfigure(0, weight=0)
    
    
    row += 1
    self.atomsOptFrame = LabelFrame(optFrame, text='Atoms Shown', grid=(row,0))
    self.atomsOptFrame.grid_columnconfigure(0, weight=1)
    self.atomsOptFrame.grid_columnconfigure(1, weight=1)
    self.atomsOptFrame.grid_columnconfigure(2, weight=1)

    tipText = 'Toggles particular kinds of atom group on and off in the chart (for all residues with that kind of group) '
    self.atomSelector=PartitionedSelector(self.atomsOptFrame, self.toggleAtom,
                                          maxRowObjects=10, grid=(0,0),
                                          gridSpan=(1,3), sticky='nsew',
                                          tipText=tipText)
    
    tipTexts = ['Switches on assignment display and connections for atoms in "set 1"',
                'Switches on assignment display and connections for atoms in "set 2"',
                'Switches on assignment display and connections for all atom types ']
    texts = ['Set1 Atoms','Set2 Atoms','All atoms']
    commands = [self.showSet1Atoms,self.showSet2Atoms,self.showAllAtoms]
    buttonList = ButtonList(self.atomsOptFrame, texts=texts,
                            commands=commands, grid=(1,0),
                            tipTexts=tipTexts)
    self.set1AtomsButton, self.set2AtomsButton, self.allAtomsButton = buttonList.buttons
    
    
    row += 1
    self.atomsFrame = Frame(optFrame, relief='flat', bd=0, grid=(row,0))
    self.atomsFrame.grid_columnconfigure(0, weight=0)
    
    
    row += 1
    self.residueFrame = LabelFrame(optFrame, text='Residues shown', grid=(row, 0))
    self.residueFrame.grid_columnconfigure(0, weight=1)
    self.residueFrame.grid_columnconfigure(1, weight=1)

    tipText = 'Toggles particular kinds of residue on and off in the chart, still observing the above atom site selections; "off" residues are merely greyed out not removes'
    self.residueSelector=PartitionedSelector(self.residueFrame, self.toggleResidueType,
                                             maxRowObjects=10, gridSpan=(1,3),
                                             grid=(0,0), tipText=tipText)
    
    tipTexts = ['Set all residue types, and any of the selected atom sites they contain, to be visible in the main chart',
                'Make an inverse residue type selection, so that types that were not visible in the main chart now are, and those that were are now greyed out']
    texts = ['All residues','Invert Selection']
    commands = [self.showAllResidues,self.showInverseResidues]
    buttonList = ButtonList(self.residueFrame, texts=texts, commands=commands,
                            grid=(1,0), tipTexts=tipTexts)
    self.residuesAllButton, self.residuesInvertButton = buttonList.buttons
    
    row += 1
    frame = LabelFrame(optFrame, text='Isotope Labelling', grid=(row, 0))
    frame.expandGrid(3,5)
    
    label= Label(frame, text = 'Pattern:', grid=(0,0))
    tipText = 'If required, selects an isotope labelling scheme to use for colouring the atom site nodes'
    self.labellingSchemePulldown = PulldownList(frame, self.changeLabellingScheme,
                                                grid=(0,1), tipText=tipText)
                                                
    label= Label(frame, text = 'Isotope:', grid=(1,0))
    tipText = 'Selects which kind of nuclear isotope atom site colours are derived from using the selected scheme; typically 13C for checking C-C connections'
    self.isotopePulldown = PulldownList(frame, texts=STANDARD_ISOTOPES, index=1,
                                        objects=isotopes,  grid=(1,1),
                                        callback=self.changeIsotope, tipText=tipText)
                                                
    label= Label(frame, text = ' Min Isotope Fraction:', grid=(0,2))
    tipText = 'Specifies the minimum fractional spin-active isotopic incorporation for an atom to be considered "labelled"'
    self.minFractionEntry = FloatEntry(frame, text=0.25, width=8,
                                       grid=(0,3), tipText=tipText)
    
    
    label= Label(frame, text = ' Colour Gradient:', grid=(1,2))
    tipText = 'Whether, for labelled atom sites, to use a colour gradient with various hues to show isotopic incorporation, or otherwise only two colours for labelled & unlabelled'
    self.colorGradCheck = CheckButton(frame, selected=True, grid=(1,3), tipText=tipText)
    
    l, u = 5, 0
    self.labelColor = standardColors[l]
    self.unlabelColor = standardColors[u]
    
    label= Label(frame, text = 'Labelled Colour:', grid=(0,4))
    tipText = 'Selects which colour is use to display the atom site nodes in the main chart if they are unlabelled, according to the stated settings'
    self.labelColorPulldown = PulldownList(frame, callback=self.setLabelColor,
                                           texts=[c.name for c in standardColors],
                                           objects=standardColors,
                                           colors=[c.hex for c in standardColors],
                                           grid=(0,5), index=l, tipText=tipText)
    
   
    label= Label(frame, text = 'Unlabelled Colour:', grid=(1,4))
    tipText = 'Selects which colour is use to display the atom site nodes in the main chart if they are isotopically labelled, according to the stated settings'
    self.unlabelColorPulldown = PulldownList(frame, callback=self.setUnlabelColor,
                                           texts=[c.name for c in standardColors],
                                           objects=standardColors,
                                           colors=[c.hex for c in standardColors],
                                           grid=(1,5), index=u, tipText=tipText)

    
    row += 1
    self.optionFrame = Frame(optFrame, relief='flat',  bd=0, grid=(row,0))
    self.optionFrame.grid_columnconfigure(0, weight=1)
    self.optionFrame.grid_columnconfigure(1, weight=1)
    self.optionFrame.grid_columnconfigure(2, weight=1)
 
    self.resColorPulldown = PulldownMenu(self.optionFrame, self.changeResColor,
                                         ['None','Rainbow','Charge','Hydrophobicity'],
                                         selected_index=0, label_color=TITLE_COLOR)
    #self.resColorPulldown.grid(row=0, column=0, sticky='ne' )



    tipTexts = ['Manually trigger a redraw of the main chart, in light of new assignments etc',]
    texts = ['Update',]
    commands = [self.updateAllConnections,]
    self.bottomButtons = UtilityButtonList(tabbedFrame.sideFrame, texts=texts,
                                           commands=commands, helpUrl=self.help_url,
                                           grid=(0,0), sticky = 'e', tipTexts=tipTexts)

    self.administerNotifiers(self.registerNotify)
    self.updateChains()

  def setLabelColor(self, color):
  
    self.labelColor = color

  def setUnlabelColor(self, color):
  
    self.unlabelColor = color

  def getIsotopes(self):
  
    store = self.project.currentChemElementStore
    isotopes = []
    
    for code in STANDARD_ISOTOPES:
      
      massNumber = 13
      symbol = 'C'
      for i in range(len(code)):
        if code[i].isalpha():
          massNumber = int(code[:i])
          symbol = code[i:]
          break
    
      element = store.findFirstChemElement(symbol=symbol)
      isotope = element.findFirstIsotope(massNumber=massNumber)
      isotopes.append(isotope)
      
    return isotopes

  def changeIsotope(self, isotope):
  
    self.isotope = isotope

  def toggleTab(self, index):
  
    if index == 0:
      self.updateMolecule()
      self.updateSpecSelection()
      self.updateAllConnections()
      self.canvasFrame.refresh()
      self.canvas.xview('moveto', 0)
      self.canvas.yview('moveto', 0)
    
    else:
      self.updateLabellingSchemes()
      self.updateSpecSelection()
      self.updateAtomSelection()
      self.updateResidueTypeSelection()
      self.updateButtons()

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

  def changeLabellingScheme(self, scheme):
    
    if scheme is not self.labellingScheme:
      self.labellingScheme = scheme
      
  def updateMoLabelFraction(self, molLabelFraction=None):

    if molLabelFraction.labeledMixture is self.labellingScheme:
      self.updateLabellingSchemes()

  def updateLabellingSchemes(self, notifyScheme=None):
  
    if notifyScheme and (notifyScheme is not self.labellingScheme):
      return
  
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
    
    self.labellingSchemePulldown.setup(names, labellings, index)  
    
    
  def administerNotifiers(self, notifyFunc):

    notifyFunc(self.updateAssign, 'ccp.nmr.Nmr.ResonanceSet', 'delete')
    notifyFunc(self.updateAssign, 'ccp.nmr.Nmr.ResonanceSet', '__init__')

    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.Experiment', 'setName')
    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.DataSource', 'setName')
    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.DataSource', 'delete')
    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.DataSource', '__init__')
    
    notifyFunc(self.updateChains, 'ccp.molecule.MolSystem.Chain', 'delete')
    notifyFunc(self.updateChains, 'ccp.molecule.MolSystem.Chain', '__init__')
    notifyFunc(self.updateResidue, 'ccp.molecule.MolSystem.Residue', 'setSeqCode')

    
    for func in ('__init__', 'delete') :
      notifyFunc(self.updateLabellingSchemes, 'ccp.molecule.ChemCompLabel.LabelingScheme', func)
      notifyFunc(self.updateLabellingSchemes, 'ccp.molecule.LabeledMolecule.LabeledMixture', func) 
    
    for func in ('__init__', 'delete', 'setWeight'): 
      notifyFunc(self.updateMoLabelFraction, 'ccp.molecule.LabeledMolecule.MolLabelFraction', func) 
    
     
    notifyFunc(self.updateSpectrumColors, 'ccpnmr.Analysis.AnalysisSpectrum', 'setPosLevels')
    notifyFunc(self.updateSpectrumColors, 'ccpnmr.Analysis.AnalysisSpectrum', 'setNegLevels')

  def updateResidue(self, residue):
  
    if residue.chain is self.chain:
      self.updateAfter()

  def updateSpectrumColors(self, object=None):
 
    self.updateSpecSelection()
    self.updateAllConnections()

  def toggleResidueType(self, ccpCode):
  
    if self.residuesDict.get(ccpCode) is None:
      self.residuesDict[ccpCode] = False
    else:
      self.residuesDict[ccpCode] = not self.residuesDict[ccpCode]

    for residue in self.chain.residues:
      if getResidueCode(residue) == ccpCode:
        self.canvasDict[residue]['inactive'] = not self.residuesDict[ccpCode]    
        self.updateResidueAtoms(residue)

  def getCcpCodes(self):
  
    ccpCodeList = []
    
    if self.chain:
      residues = self.chain.residues
      for residue in residues:
        ccpCode = getResidueCode(residue)
 
        if ccpCode not in ccpCodeList:
          ccpCodeList.append(ccpCode)
 
      ccpCodeList.sort()
      
    return ccpCodeList

  def showAllResidues(self):

    residues = self.chain.residues
    for residue in residues:
      self.residuesDict[getResidueCode(residue)] = True
      self.canvasDict[residue]['inactive'] = False
      self.updateResidueAtoms(residue)

    self.updateResidueSelector()
    self.updateAllConnections()
  
  def showInverseResidues(self):
  
    residues = self.chain.residues
    for residue in residues:
      boolean = self.canvasDict[residue].get('inactive')
      if boolean is not None:
        self.canvasDict[residue]['inactive'] = not boolean
      else:
        self.canvasDict[residue]['inactive'] = True
      self.updateResidueAtoms(residue)
	
    self.updateResidueSelector()
    self.updateAllConnections()

  def updateButtons(self):
 
    toolTipA = self.set1AtomsButton.toolTip
    toolTipB = self.set2AtomsButton.toolTip
    
    tip = 'Switches on assignment display and connections for atoms '
 
    if self.chain:
      self.set1AtomsButton.enable()
      self.set2AtomsButton.enable()
      self.allAtomsButton.enable()
      self.residuesAllButton.enable()
      self.residuesInvertButton.enable()
      
      if self.chain.molecule.molType == 'protein':
        self.set1AtomsButton.configure(fg='black', state='normal',
                                       text='backbone, alpha & beta')
        self.set2AtomsButton.configure(fg='black', state='normal',
                                       text='backbone, alpha - delta')
        toolTipA.text = tip + 'in the protein backbone and at alpha and beta positions'
        toolTipB.text = tip + 'in the protein backbone and at side chain positions from alpha to delta, inclusive'

      elif self.chain.molecule.molType == 'DNA':
        self.set1AtomsButton.configure(fg='black', state='normal', text='Base')
        self.set2AtomsButton.configure(fg='black', state='normal', text='Deoxyribose')
        toolTipA.text = tip + 'in the DNA bases'
        toolTipB.text = tip + 'in the deoxyribose sugar'


      elif self.chain.molecule.molType == 'RNA':
        self.set1AtomsButton.configure(fg='black', state='normal', text='Base')
        self.set2AtomsButton.configure(fg='black', state='normal', text='Ribose')
        toolTipA.text = tip + 'in the RNA bases'
        toolTipB.text = tip + 'in the ribose sugar'
        
      elif self.chain.molecule.molType == 'carbohydrate':
        self.set1AtomsButton.configure(fg='grey', state='normal', text='Pos 1,4,6')
        self.set2AtomsButton.configure(fg='grey', state='normal', text='Other Pos')
        toolTipA.text = tip + 'at 1, 4 and 6 positions in the sugar ring'
        toolTipB.text = tip + 'not at 1, 4 and 6 positions in the sugar ring'

      else:
        self.set1AtomsButton.configure(fg='grey', state='disabled', text='Set1 Atoms')
        self.set2AtomsButton.configure(fg='grey', state='disabled', text='Set2 Atoms')
        toolTipA.text = tip + 'within the first set of atoms'
        toolTipB.text = tip + 'within the second set of atoms'
        
    else:
      self.set1AtomsButton.disable()
      self.set2AtomsButton.disable()
      self.allAtomsButton.disable()
      self.residuesAllButton.disable()
      self.residuesInvertButton.disable()
      self.set1AtomsButton.configure(fg='grey',state='disabled',text='Set1 Atoms')
      self.set2AtomsButton.configure(fg='grey',state='disabled',text='Set2 Atoms')
      toolTipA.text = tip + 'within the first set of atoms'
      toolTipB.text = tip + 'within the second set of atoms'

  def refreshGridding(self, dx, dy):
  
     self.gridSpecSelection(int(dx/75)-1)

  def showAllAtoms(self):

    self.setAtomSelection([x[0] for x in self.atomNamesList])
    self.updateAllResidueAtoms()
  
  def showSet1Atoms(self):
  
    atomSelection = []
    molTypes = set([r.molType for r in self.chain.residues])

    if 'protein' in molTypes:
      atomSelection += ['CO','NH','a','b']
    
    if 'DNA' in molTypes:
      atomSelection += ['1','2','3','4','5','6','7','8','9']
    elif 'RNA' in molTypes:
      atomSelection += ['1','2','3','4','5','6','7','8','9']
    elif 'carbohydrate' in molTypes:
      atomSelection = ('1','4','6')

    self.setAtomSelection(atomSelection)
    self.updateAllResidueAtoms()
  
  def showSet2Atoms(self):
  
    atomSelection = []
    molTypes = set([r.molType for r in self.chain.residues])

    if 'protein' in molTypes:
      atomSelection += ['CO','NH','a','b','g','g1','g2','d','d1','d2']
    
    if 'DNA' in molTypes:
      atomSelection += ["1'","2'","3'","4'","5'"]
    elif 'RNA' in molTypes:
      atomSelection += ["1'","2'","3'","4'","5'"]
    
    if 'carbohydrate' in molTypes:
      atomSelection = ('2','3','5','7','8','9')

    self.setAtomSelection(atomSelection)
    self.updateAllResidueAtoms()

  def setAtomSelection(self, atomNames):
  
    i = 0
    for atomName, molType in self.atomNamesList:
      if atomName in atomNames:
        self.atomNamesDict[atomName] = True
        self.atomSelector.setButtonState(i, True)
      else:
        self.atomNamesDict[atomName] = None
        self.atomSelector.setButtonState(i, False)
      i += 1

  def toggleAtom(self, atomName):
  
    if self.atomNamesDict.get(atomName, True):
      self.atomNamesDict[atomName] = False
    else:
      self.atomNamesDict[atomName] = True

    for residue in self.chain.residues:
      if atomName in self.canvasDict[residue]['atm'].keys():
        self.clearResidueAtoms(residue)
        self.updateResidueAtoms(residue)

  def toggleAtomSelection(self, i):

    atomName, molType = self.atomNamesList[i]

    if self.atomNamesDict.get(atomName, True):
      self.atomNamesDict[atomName] = False
    else:
      self.atomNamesDict[atomName] = True

    for residue in self.chain.residues:
      if atomName in self.canvasDict[residue]['atm'].keys():
        self.clearResidueAtoms(residue)
        self.updateResidueAtoms(residue)
    
  def updateAllResidueAtoms(self):
    
    
    for residue in self.chain.residues:
      self.clearResidueAtoms(residue)
      self.updateResidueAtoms(residue)
      
    self.updateAllConnections()
    self.canvasFrame.refresh()
 
  def clearResidueAtoms(self,residue):
  
     for atomName2 in self.canvasDict[residue]['atm'].keys():
      self.canvas.delete(self.canvasDict[residue]['atm'][atomName2])
      self.canvas.delete(self.canvasDict[residue]['atmTxt'][atomName2])
      self.canvasDict[residue]['atm'][atomName2] = None
      self.canvasDict[residue]['atmTxt'][atomName2] = None
  
  def showSetConnections(self, option):
 
    boundsDict = {0:(0,0),1:(0,1),
                  2:(0,4),3:(2,4),
                  4:(5,None),}
 
    minDist, maxDist = boundsDict.get(option, (0, None))
   
    self.maxConnectionDist = maxDist
    self.minConnectionDist = minDist
    self.updateAllConnections()

  def changeResColor(self, i, option):
  
    if option != self.residueColorMode:
      self.residueColorMode = option
        
      for residue in self.chain.residues:
        tlc = getResidueCode(residue)
        num = residue.seqCode
        color =  self.getResidueColor(tlc, num)
        self.canvas.itemconfigure(self.canvasDict[residue]['tlc'], fill=color)
        self.canvas.itemconfigure(self.canvasDict[residue]['num'], fill=color)
   
  def changeMolecule(self, chain):
  
    if chain is not self.chain:
      self.chain = chain
      self.atomNamesList = []
      self.update()
    
  def toggleSpectrum(self, spectrum):
  
    if self.spectraDict[spectrum]:
      self.spectraDict[spectrum] = False
    else:
      self.spectraDict[spectrum] = True
       
  def greekSort(self,dataList):
  
    sortList = []
    for x, molType in dataList:
      sortName = x
      sortName = re.sub('(.+\')', 'zzz@\\1', sortName)
      sortName = re.sub('^(\d)','zz@\\1', sortName)
      sortName = re.sub('N(\S*)','\'\'@N\\1',sortName)
      sortName = re.sub('CO','\'\'@CO',sortName)
      sortName = re.sub('g', 'c@g', sortName)
      sortName = re.sub('z', 'f@z', sortName)
      sortList.append( (sortName,x,molType) )
    sortList.sort()
    dataList = [e[1:] for e in sortList]
    
    return dataList

  def updateResidueTypeSelection(self):
  
    ccpCodes = self.getCcpCodes()    
    
    colors = []
    for ccpCode in ccpCodes:
      colors.append( HILIGHT_COLOR )

    self.residueSelector.update(objects=ccpCodes,labels=ccpCodes,colors=colors)
   
  def updateAtomSelection(self):
    
    self.atomNamesList = self.greekSort(self.atomNamesList)
        
    nameList = []
    colors = []
    selected = []
    for atomName, molType in self.atomNamesList:
      colors.append( HILIGHT_COLOR )
      
      if (molType == 'protein') and (atomName not in ('NH','CO','N')):
        name = ''.join([unicodeGreek.get(x, x) for x in atomName])
        nameList.append(name)
      else:
        nameList.append(atomName)
      
      if self.atomNamesDict.get(atomName, True):
        selected.append(atomName)

    self.atomSelector.update(objects=[x[0] for x in self.atomNamesList],
                             labels=nameList,
                             colors=colors,
                             selected=selected)
        
  def updateSpecSelection(self):
  
    spectra      = []
    spectraNames = []
    for expt in self.nmrProject.sortedExperiments():
       eName = expt.name
       for spec in expt.sortedDataSources():
         if spec.dataType == 'processed':
           spectraNames.append('%s:%s' % (eName, spec.name))
           spectra.append( spec )
           if self.spectraDict.get( spec ) is None:
             self.spectraDict[spec] = False
        
    self.spectra = spectra
    analysisProject = self.analysisProject
    colors = []
    fonts  = []
    for i in range(len(spectra)):
      font = self.font
      
      hexColors = spectra[i].analysisSpectrum.posColors
      hexColor = hexColors[int(0.7*len(hexColors))]
        
      colors.append( hexColor )
      fonts.append(font)
    
    self.specSelector.update(objects=spectra,labels=spectraNames,
                             colors=colors,fonts=fonts)
    for i in range(len(spectra)):
      self.specSelector.setButtonState(i, self.spectraDict[spectra[i]])

  def updateConnections(self, spectrum):
  
    delete = self.canvas.delete
    canvasDict = self.canvasDict
    if canvasDict.get(spectrum):
      for item in canvasDict[spectrum]:
        delete(item)
      canvasDict[spectrum] = None

    if not self.spectraDict.get(spectrum):
       return
     
    hexColors = spectrum.analysisSpectrum.posColors
    hexColor = hexColors[int(0.7*len(hexColors))]

    chain = self.chain
    connections = []

    nDim = spectrum.numDim
    
    maxConnectionDist = self.maxConnectionDist
    minConnectionDist = self.minConnectionDist
    drawConnection = self.drawConnection
    
    boundDims = []
    for dataDim1, dataDim2 in getOnebondDataDims(spectrum):
      boundDims.append(set([dataDim1.dim-1, dataDim2.dim-1]))
    
    includePredicted = self.includePredictedCheck.get()
    for peakList in spectrum.peakLists:
      if (not includePredicted) and peakList.isSimulated:
        continue

      for peak in peakList.peaks:
        dimAtoms = {}
        for i in range(nDim):
          dimAtoms[i] = set()
          
        if peak.peakContribs:
          for peakContrib in peak.peakContribs:
             
            for contrib in peakContrib.peakDimContribs:
              dim = contrib.peakDim
              
              resonanceSet = contrib.resonance.resonanceSet
              if resonanceSet:
                for atomSet in resonanceSet.atomSets:
                  atom = atomSet.findFirstAtom()
                  residue = atom.residue
                  if residue.chain is chain:
                    dimAtoms[dim.dim-1].add((atom,residue))
          
        else:
          for dim in peak.peakDims:
            for contrib in dim.peakDimContribs:
              resonanceSet = contrib.resonance.resonanceSet
              if resonanceSet:
                for atomSet in resonanceSet.atomSets:
                  atom = atomSet.findFirstAtom()
                  residue = atom.residue
                  if residue.chain is chain:
                    dimAtoms[dim.dim-1].add((atom,residue))
    
        for i in range(nDim-1):
          for j in range(i+1, nDim):
            atomPairs = []
 
            if set([i,j]) in boundDims:
              for atomA, residueA in dimAtoms[i]:
                for atomB, residueB in dimAtoms[j]:
 
                   if areAtomsBound(atomA, atomB):
                     atomPairs.append((atomA, residueA, atomB, residueB))
 
            if not atomPairs:
              for atomA, residueA in dimAtoms[i]:
                for atomB, residueB in dimAtoms[j]:
                  atomPairs.append((atomA, residueA, atomB, residueB))
 
            for atomA, residueA, atomB, residueB in atomPairs:
 
              if canvasDict[residueA].get('inactive'):
                continue
 
              if canvasDict[residueB].get('inactive'):
                continue
 
              if atomA is atomB:
                continue
 
              diff = abs(residueA.seqId - residueB.seqId)
              if ((maxConnectionDist is None) or (diff <= maxConnectionDist)) and \
                  (diff >= minConnectionDist):
                drawConnection(hexColor, spectrum, atomA, atomB)
            
          

  def drawConnection(self, color, spectrum, atom1, atom2):
  
    residue1 = atom1.residue
    residue2 = atom2.residue
    
    canvasDict = self.canvasDict
    canvas     = self.canvas
    
    atomGroups = [e[1] for e in residue1.graphAtomList]
    atomNames  = [e[2] for e in residue1.graphAtomList]
    name1 = atomNames[0]
    # TBD possible failure below ?
    for j in range(len(atomNames)):
      if atom1 in atomGroups[j]:
        name1 = atomNames[j]
        break

    atomGroups = [e[1] for e in residue2.graphAtomList]
    atomNames  = [e[2] for e in residue2.graphAtomList]
    name2 = atomNames[0]
    for j in range(len(atomNames)):
      if atom2 in atomGroups[j]:
        name2 = atomNames[j]
        break
    
    if not canvasDict.get(spectrum):
      canvasDict[spectrum] = []

    if canvasDict[residue1]['atm'][name1]:
      coords = canvas.coords(canvasDict[residue1]['atm'][name1])
      coords1= [coords[0]+8, coords[1]+8]
      circ1 = self.createAtomCircle(canvas, coords1[0], coords1[1], 9, outline=color, width=3)
      canvasDict[spectrum].append( circ1 )

    if canvasDict[residue2]['atm'][name2]:
      coords = canvas.coords(canvasDict[residue2]['atm'][name2])
      coords2= [coords[0]+8, coords[1]+8]
      circ2 = self.createAtomCircle(canvas, coords2[0], coords2[1], 9, outline=color, width=3)
      canvasDict[spectrum].append( circ2 )

    if canvasDict[residue2]['atm'][name2] and canvasDict[residue1]['atm'][name1]:
      if canvasDict[residue2]['atm'][name2] != canvasDict[residue1]['atm'][name1]:
        x1 = coords1[0]
        y1 = coords1[1]
        x2 = coords2[0]
        y2 = coords2[1]
 
        dx =  x1 - x2
        dy =  y1 - y2
 
        hyp = pow((dx*dx) + (dy*dy), 0.5)

        if dx != 0:
          cosine= dx/hyp
          x1 -= 8 * cosine
          x2 += 8 * cosine
 
        if dy != 0:
          sine = dy/hyp
          y1 -= 8 * sine
          y2 += 8 * sine
 
        line1 = canvas.create_line(x1,y1,x2,y2,fill= color ,width=4)
        line2 = canvas.create_line(x1,y1,x2,y2,fill='black',width=1)
        canvas.lower(line2)
        canvas.lower(line1)
        canvasDict[spectrum].append( line1 )
        canvasDict[spectrum].append( line2 )
  
  def updateAssign(self, resonanceSet):
    
    residue = resonanceSet.findFirstAtomSet().findFirstAtom().residue
    resonances = resonanceSet.resonances
      
    if self.canvasDict.get(residue) is not None:
      spectra = {}
      for resonance in resonances:
        for contrib in resonance.peakDimContribs:
          spec = contrib.peakDim.peak.peakList.dataSource
          spectra[spec] = True
        
      self.updateResidueAtoms(residue)
      for spectrum in spectra.keys():
        self.updateConnections(spectrum)

  def getGroupName(self, atomName, molType=None):
  
    if molType:
      if  molType == 'protein':
        if atomName in ('NH','N','CO'):
          return atomName
        
        letters = [unicodeGreek.get(x, x) for x in atomName]
        return ''.join(letters)
          
      else:
        return atomName
    else:
      return atomName
 

  def getIsotopeColors(self, isoScheme, residue, atomGroup, minIsoFrac=0.1, colorGrad=True):

    chemElement = self.isotope.chemElement
    atomNames = [a.name for a in atomGroup if a.chemAtom.chemElement is chemElement]

    if not atomNames:
      fraction = 0.0
      
    else:
      isotope = '%d%s' % (self.isotope.massNumber,chemElement.symbol)
      
      if isoScheme.className == 'LabelingScheme':
        ccpCode = residue.ccpCode
        molType = residue.molType
 
        chemCompLabel = isoScheme.findFirstChemCompLabel(ccpCode=ccpCode,
                                                         molType=molType)
 
        if not chemCompLabel:
          natAbun = residue.root.findFirstLabelingScheme(name='NatAbun')
 
          if natAbun:
            chemCompLabel = natAbun.findFirstChemCompLabel(ccpCode=ccpCode,
                                                           molType=molType)
 
        if chemCompLabel:
          isotopomers = chemCompLabel.isotopomers
          fraction = 0.0
 
          for atomName in atomNames:
            fracDict = getIsotopomerSingleAtomFractions(isotopomers,
                                                        atomName, 1)
            fraction += fracDict.get(isotope, 1.0)
 
          fraction /= float(len(atomNames))
      
      else:
        fraction = 0.0
 
        for atomName in atomNames:
          fracDict = singleAtomFractions(isoScheme, residue.molResidue.serial, atomName)
          fraction += fracDict.get(isotope, 1.0)
 
        fraction /= float(len(atomNames))
         
        
    
    if fraction < minIsoFrac:
      color1 = self.unlabelColor.hex
      color2 = inverseGrey(color1)
      color3 = 'grey83'
      
    else:  
      if colorGrad:
        r1, g1, b1 = self.unlabelColor.rgb()
        r2, g2, b2 = self.labelColor.rgb()
        
        fr = 255 * fraction
        un = 255 - fr
        rgb = ((r1*un)+(fr*r2),
               (g1*un)+(fr*g2),
               (b1*un)+(fr*b2))
               
        color1 = '#%02x%02x%02x' % rgb
        color2 = inverseGrey(color1)
        color3 = 'grey83'
      
      else:
        color1 = self.labelColor.hex
        color2 = inverseGrey(color1)
        color3 = 'grey83'
          
    return color1, color2, color3


  def updateResidueAtoms(self, residue):

    sFont = self.sFont 
    canvas        = self.canvas
    canvasDict    = self.canvasDict
    bindDict      = self.bindDict
    circle        = self.createAtomCircle
    atomNamesDict = self.atomNamesDict
    atomNamesList = self.atomNamesList
    create_text   = canvas.create_text
    itemconfigure = canvas.itemconfigure
    cDelete       = canvas.delete
    rDict         = canvasDict[residue]
    getName       = self.getGroupName
    colorGrad     = self.colorGradCheck.get()
    minIsoFrac    = self.minFractionEntry.get() or 0.0
    isoScheme     = self.labellingScheme
    getIsotopeColors = self.getIsotopeColors
    
    if not hasattr(residue, 'graphAtomList'):
      residue.graphAtomList = self.makeGraphAtomList(residue)
      
    x = rDict['x']
    y = rDict['y']
      
    atomNames  = [e[2] for e in residue.graphAtomList]
    atomGroups = [e[1] for e in residue.graphAtomList]
    
    localY = y + 36
    for j, atomName in enumerate(atomNames):
      atomGroup = atomGroups[j]
      name = getName(atomName, residue.molType)
      if len(name) == 1:
        font = sFont % 8
      else:
        font = sFont % 7

      color1 = 'grey83'
      color2 = 'black'
      color3 = 'grey'
      if rDict.get('inactive'):
        color1 = 'grey90'
        color2 = 'grey'
        color3 = 'grey90'
        
      elif isoScheme:
        color1, color2, color3 = getIsotopeColors(isoScheme, residue, atomGroup,
                                                  minIsoFrac, colorGrad)
      
      else:
        for atom in atomGroup:
          if isAtomAssigned(atom, toPeaks=True):
            color1 = 'black'
            color2 = 'white'
            color3 = 'grey50'
            break

      aDict = rDict['atm']
      if aDict.get(atomName) is None:
        if (atomName, residue.molType) not in atomNamesList:
          atomNamesList.append( (atomName, residue.molType) )
        if atomNamesDict.get(atomName, True):
          circ = circle(canvas, x, localY, 8, fill=color1, outline=color3)
          text = create_text(x, localY+1, text=name, font=font, fill=color2)
          bindDict[circ] = residue
          bindDict[text] = residue
          aDict[atomName] = circ
          rDict['atmTxt'][atomName] = text
          localY += 23
      else:
        if atomNamesDict.get(atomName, True):
          circ = aDict[atomName]
          text = rDict['atmTxt'][atomName]
          bindDict[circ] = residue
          bindDict[text] = residue
          itemconfigure(circ,fill=color1, outline=color3)
          itemconfigure(text, fill=color2)
          localY += 23
        else:
          cDelete(aDict[atomName])
          cDelete(rDict['atmTxt'][atomName])
          aDict[atomName] = None
          rDict['atmTxt'][atomName] = None
 
  def toggleResidue(self, *opt):
  
    self.canvasFrame.removeMenu()
  
    if not self.chain:
      return

    event = opt[0]
    x = self.canvas.canvasx(event.x)
    y = self.canvas.canvasy(event.y)
    obj = self.canvas.find('closest',x,y)[0]
    residue = self.bindDict.get(obj)
    
    if not residue:
      return
    
    boolean = self.canvasDict[residue].get('inactive')
    
    if boolean is not None:
      self.canvasDict[residue]['inactive'] = not boolean
    else:
      self.canvasDict[residue]['inactive'] = True
    
    self.updateResidueSelector()
    self.updateResidueAtoms(residue)
    self.updateAllConnections()

  def updateResidueSelector(self):

    ccpCodes = self.getCcpCodes()
    ccpCodesDict = {}

    for residue in self.chain.residues:
      if not self.canvasDict[residue].get('inactive'):
        ccpCodesDict[getResidueCode(residue)] = True

    for i in range(len(ccpCodes)):
      if ccpCodesDict.get(ccpCodes[i]):
        self.residueSelector.setButtonState(i, True)
	self.residuesDict[ccpCodes[i]] = True
      else:
        self.residueSelector.setButtonState(i, False)
        self.residuesDict[ccpCodes[i]] = False

  def updateMolecule(self):
  
    if not self.chain:
      return
      
    canvas = self.canvas
    lower = canvas.lower
    itemconfigure = canvas.itemconfigure
    create_text = canvas.create_text
    getBbox = canvas.bbox
    setCoords = canvas.coords
    create_rectangle = canvas.create_rectangle
    cDict = self.canvasDict
    bDict = self.bindDict = {}
    updateResidueAtoms = self.updateResidueAtoms
    font = self.font
    
    x = 18
    i = 0
    y = 14
 
    residues = self.chain.sortedResidues()
    for residue in residues:
        
      tlc = getResidueCode(residue)
      num = residue.seqCode
      color = self.getResidueColor(tlc, num)
        
      if cDict.get(residue) is None:
        cDict[residue] = {}

      rDict = cDict[residue]
      if rDict.get('tlc') is None:
        textItemNum = create_text(x, y, text=num, font=font, fill=color, justify='center')
        textItemTlc = create_text(x, y+16, text=tlc, font=font, fill=color)
        
        bbox = getBbox(textItemTlc)
        dxA = bbox[2] - bbox[0]
 
        bbox = getBbox(textItemNum)
        dxB = bbox[2] - bbox[0]
 
        width = max((dxA, dxB, 24))+2
        
        box = create_rectangle(x,y-10,x+width,y+24, outline='grey', fill='grey90')
        setCoords(textItemNum, x+width/2, y)
        setCoords(textItemTlc, x+width/2, y+16)
        lower(box)
        
        bDict[box]         = residue
        bDict[textItemNum] = residue
        bDict[textItemTlc] = residue
        rDict['tlc'] = textItemTlc
        rDict['num'] = textItemNum
        rDict['box'] = box
        rDict['atm'] = {}
        rDict['atmTxt'] = {}
        rDict['y'] = y
        rDict['x'] = x+width/2
      else:
        itemconfigure(rDict['tlc'], fill=color)
        itemconfigure(rDict['num'], fill=color, text=num)
        width = 2*(rDict['x']-x)
        
      updateResidueAtoms(residue)
      x += width + 2
      i += 1
            
    # clear unwanted residue gfx
    delFunc = canvas.delete
    for ID in cDict.keys():
      if (ID not in residues) and (ID not in self.spectra):
        delFunc(cDict[ID]['box'])
        delFunc(cDict[ID]['tlc'])
        delFunc(cDict[ID]['num'])
        cDict[ID]['box'] = None
        cDict[ID]['tlc'] = None
        cDict[ID]['num'] = None
        for name in cDict[ID]['atm'].keys():
          delFunc(cDict[ID]['atm'][name])
          delFunc(cDict[ID]['atmTxt'][name])
          cDict[ID]['atm'][name] = None
          cDict[ID]['atmTxt'][name] = None
           
  def makeGraphAtomList(self,residue):

     atomList  = []
     heavies   = []
     hydrogens = []
     for atom in residue.atoms:
       chemAtom = atom.chemAtom
       #if chemAtom.elementSymbol in ('N','C','P'):
       if chemAtom.elementSymbol in ('N','C'):
         heavies.append((atom,chemAtom))
       elif chemAtom.elementSymbol == 'H':
         hydrogens.append((atom,chemAtom))
 
     hydrogenChemAtoms = [e[1] for e in hydrogens]
     hydrogenAtoms     = [e[0] for e in hydrogens]
     newNames = []
 
     for heavy in heavies:
       chemAtom  = heavy[1]
       atomGroup = [ heavy[0], ]
       name      = chemAtom.name
       symb      = chemAtom.elementSymbol

       if name == 'N':
         newName = '\'\'@N'
         if residue.ccpCode == 'PRO':
           name = 'N'
         else:
           name = 'NH'
       elif name == 'C':
         newName = '\'\'@C'
         name = 'CO'
       else:
         name    = re.sub(symb,'', name)
         newName = name
         newName = re.sub('(.+\')','zzz@\\1', newName)
         newName = re.sub('^(\d)','zz@\\1', newName)
         if residue.molType == 'protein':
           name = name.lower()
           newName = re.sub('G','C@G', newName)
 
       for bond in chemAtom.chemBonds:
         for chemAtom2 in bond.chemAtoms:
           if chemAtom2 in hydrogenChemAtoms:
             atom = hydrogenAtoms[hydrogenChemAtoms.index(chemAtom2)]
             atomGroup.append( atom )
             
       if newName in newNames:
         i = newNames.index(newName)
         newGroup = atomList[i][1]
         for atom in atomGroup:
           newGroup.append(atom)
         atomList[i] = (newName, newGroup,name)
       else:
         newNames.append(newName)
         atomList.append( (newName,atomGroup,name) )

     atomList.sort()
     return atomList

  def createAtomCircle(self, canvas, centerX, centerY, radius,
                       fill=None, outline='black', width=1):
  
    x1 = centerX - radius
    y1 = centerY - radius
    x2 = centerX + radius
    y2 = centerY + radius
    
    if not fill:
      circle = canvas.create_oval(x1,y1,x2,y2,outline=outline,width=width)
    else:
      circle = canvas.create_oval(x1,y1,x2,y2,fill=fill,outline=outline,width=width)
    
    return circle

  def getResidueColor(self, tlc, num):

    if self.residueColorMode == 'Rainbow':
       n = num % 10
       n = n / 10.
       rgb = hsbToRgb(n, 1, 0.8)
       color = '#%02x%02x%02x' % (rgb[0]*255, rgb[1]*255, rgb[2]*255)

    
    elif self.residueColorMode == 'Charge':
      charge = RES_CHARGE.get(tlc, 0)
      if charge < 0:
        color = 'red'
      elif charge > 0:
        color = 'blue'
      else:
        color = 'black'
        
    elif self.residueColorMode == 'Hydrophobicity':
       h = H_PHOBE.get(tlc, 0)
       if h > 1:
         color = 'blue'
       elif h > -1:
         color = 'darkgreen'
       else:
         color = 'black'
       
    else:
      color = 'black'

    return color
  
  def getChainName(self, chain):
  
    return '%s:%s (%s)' % (chain.molSystem.code,chain.code,chain.molecule.molType)
  
  
  def getChains(self):
  
    chains = []
    if self.project:
      for molSystem in self.project.sortedMolSystems():
        for chain in molSystem.sortedChains():
          if chain.residues:
            chains.append(chain)
   
    return chains
   
   
  def updateAfter(self, *opt):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
  
  
  def updateChains(self, *opt):
  
    index  = 0
    texts  = []
    chains = self.getChains()
    chain  = self.chain

    if chains:
      if chain not in chains:
        chain = chains[0]

      texts = [self.getChainName(c) for c in chains]
      index = chains.index(chain)
        
    else:
      chain = None  
      
    self.molPulldown.setup(texts,chains,index)
  
    if chain is not self.chain:
      self.chain = chain
      self.update()
  
  
  def update(self):
    
    self.toggleTab(self.tabbedFrame.selected) 
    self.waiting = False
      
  def updateAllConnections(self):
  
    for spectrum in self.spectra:
      self.updateConnections(spectrum)
  
  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
    
    BasePopup.destroy(self)

"""
======================COPYRIGHT/LICENSE START==========================

IsotopeSchemeEditor.py: Part of the CcpNmr Analysis program

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
from memops.api import Implementation


from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.DataEntry import askString
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.LabelDivider import LabelDivider
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showWarning, showOkCancel, showInfo
from memops.gui.MultiWidget import MultiWidget
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from memops.general.Util import copySubTree

from ccp.general.Constants import standardResidueCcpCodes

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.MoleculeBasic import greekSortAtomNames
from ccpnmr.analysis.frames.ViewIsotopomerFrame import ViewIsotopomerFrame

import os

# ------------------------------ To-do list --------------------------------
#
#  * Optionally include all isotopes (e.g. 14C - not spin-active and very, very low abundance)
#
#  * Rotatable residue viewer with colours for isotopes
#
#  * General repository (scheme) for all useful isotopomers
#
#  * Add unusual residue types, Hydroxyproline, phosphotyrosine for starters...
#
#  * Phosphorous

def testChemCompLabelEditorMacro(argServer):

  popup = IsotopeSchemeEditor(argServer.parent, argServer.getProject())

elementSymbols = ('H','N','C',)

class IsotopeSchemeEditor(BasePopup):
  """
  **Create and Edit Per-residue Reference Isotope Schemes**
  
  This system allows the user to create schemes that describe particular
  patterns of atomic isotope labelling in terms of combinations of isotopically
  labelled forms of residues. Once constructed, these schemes may then be
  applied to a molecule of known residue sequence to gauge the levels of
  spin-active isotope incorporation in an NMR experiment. This information is
  useful in several places withing Analysis, including giving more intelligent
  assignment options and in the generation of distance restraints by matching
  peak positions to chemical shifts. Although the schemes may be used directly
  they are typically used as reference information for configuring the `Isotope
  Labelling`_ system; where isotope labelling patterns are tied to particular
  molecules and experiments.

  Because all of the different isotope labelled versions (isotopomers) of each
  residue type are described independently, a scheme can be used to estimate the
  specific amounts of incorporation present at multiple atom sites at the same
  time. For example, although a residue type may have significant levels of 13C
  at the CA and CB positions on average, there may be no form of the residue
  where CA and CB are labelled at the same time, and thus CA-CB correlations
  would not be observed in NMR.

  This popup window is divided into three main tabs, the first describes the
  overall schemes that are available; that would be applied to a molecule in a
  given situation. The second tab details the residue isotopomer components
  within the selected scheme, i.e. which labelled residue forms are present. The
  last tab displays isotopomer labelling in a graphical, three-dimensional way.
  If any isotope labelling schemes have been created or edited the user may 
  immediately save these to disk via the [Save Schemes] button to the right of
  the tabs, although these will naturally be saved when the main CCPN project
  is.

  **Reference Schemes**

  This table lists all of the reference isotope schemes that are available to
  the project. A number of standard schemes are included by default, as part of
  the main CCPN installation. However, the user is free to create new schemes,
  either from a completely blank description or by copying and modifying one of
  the existing schemes. By selecting on a isttope scheme row in the table the
  scheme is selected to be active for the whole popup and the user can see the
  contents of the scheme via the other two tabs.

  It should be noted that the user cannot edit the standard schemes provided by
  CCPN, given that these are stored with the software. Any new or copied schemes
  that the user creates will be stored inside the current CCPN project. If a new
  scheme should be made available to multiple projects, its XML file can be
  copied into the main CCPN installation, if the user has appropriate write
  access.

  **Isotopomers**
  
  The middle tab allows the user to view, and where appropriate edit, the
  isotope labelling descriptions for the residues within the current scheme
  (selected in the pulldown menu at the top). An isotope scheme is constructed
  by specifying one or more isotopomers for each residue type. Each isotopomer
  represents a different way of incorporating spin-active atom labels into a
  given kind of residue. Often there will only be one labelled form of a
  residue, and hence one isotopomer. However, with some kinds of isotope
  enrichment, for example using glycerol 13C labelled at the C2 position,
  connected labelled and unlabelled atom sites can be incorporated in
  alternative ways, resulting in distinct forms of labelling patterns that are
  not the result of a pure random mix. Knowing which labels are present at the
  same time, in the same isotopomer form, can be very important for determining
  which NMR correlations are possible.

  In general use when looking through the default, immutable reference schemes
  that come with CCPN the user can scroll through the isotopomer versions of
  each residue in the upper table. By clicking on one of these rows the lower
  table is filled with details of the amount of each kind of isotope (on
  average) at each atom site. For the lower "Atom Labels" table only one kind of
  chemical element is shown at a time, but the user may switch to  a different
  one via the "Chemical Element" pulldown.

  **Editing Isotopomers**

  When dealing with copied or new isotope schemes the user is allowed to 
  edit all aspects of the scheme. With a completely new scheme there will be no
  isotopomer records to start with and it is common practice to fill in a
  standard set of isotopomers, one for each residue type, made with a base level
  of isotope incorporation. To set this base level the user can use [Set Default
  Abundances] to manually specify values, although the default is to use natural
  abundance levels, which is appropriate in most circumstances. With the base
  levels set the [Add Default Abundance Set] will automatically fill-in a
  starting set of isotopomers for the scheme. Extra isotopomers can be added for
  a specific type of residue via the [Add New:] function and adjacent pulldown
  menu or by copying existing ones; whichever is easier. Each isotopomer has an
  editable weight to enable the user to indicate the relative abundance within a
  given residue type.

  Once a new isotopomer specification is created clicking on its row allows the 
  user to specify the isotope labelling pattern in the lower "Atom Labels"
  table. Here the user selects which kind of chemical element to consider and
  then  double-clicks to edit the "Weighting" columns in the table. The
  weightings represent the relative abundance of a given nuclear isotope at a
  given atom site. The weightings could be set as ratios, fractions, or
  percentages; it is only the relative proportion that is important. For example
  if a carbon atom site was known to have 5% Carbon-12 and 95% Carbon-13
  isotopes then the respective weights could be entered as 1 & 19 or  0.05 &
  0.95; whatever is most convenient. For efficient setup of schemes the
  [Propagate Abundances] function can be used to spread the same levels of
  incorporation over several atom sites (from the last selected row).

  **Isotopomer Structure**
  
  The last tab is an alternative way of presenting the isotope patterns present
  within the residues of the current scheme (selected in either of the first two
  tabs). Here the user selects a residue type in the upper left pulldown menu
  and then a numbered isotopomer, or an average of all isotopomers, in the right
  hand pulldown menu. The structural display will show a moveable picture of the
  residue (in a standard conformation) where unlabelled atom sites are
  represented with grey spheres, labelled sites with yellow spheres and
  intermediate incorporation with shades in between.

  It should be noted that this kind of 3D display is only possible if there is
  an idealised structure available for a residue type. This data will be 
  present for all of the regular biopolymer residues, but may be missing for
  more unusual compounds; although a lack of coordinates does not impact upon
  the isotopomer setup.

  To move and rotate the three-dimensional residue display the following
  keyboard controls may be used:
  
  * Rotate: Arrow keys
  
  * Zoom: Page Up & Page Down keys

  * Translate: Arrow keys + Control key

  Or alternatively the following mouse controls:
  
  * Rotate: Middle button click & drag
  
  * Zoom: Mouse wheel or middle button click + Shift key & drag up/down

  * Translate: Middle button click & drag + Control key

  Also an options menu appears when the right mouse button is clicked.

  .. _`Isotope Labelling`: EditMolLabellingPopup.html
 
  """

  def __init__(self, parent, project=None, *args, **kw):
  
    if not project:
      self.project    = Implementation.MemopsRoot(name='defaultUtilityProject')
    else:
      self.project = project
    
    self.waiting     = False
    self.waitingAtom = False
    self.molType     = 'protein'
    self.scheme      = None
    self.isotopomer  = None
    self.isotopomerV = False # Not None
    self.ccpCodeV    = None
    self.element     = 'C'
    self.atomLabelTuple  = None
    self.isotopes    = [x[0] for x in getSortedIsotopes(self.project, 'C')]
    self.defaultAbun = {}
    
    BasePopup.__init__(self, parent=parent, title='Molecule : Reference Isotope Schemes', **kw)

  def body(self, guiFrame):
  
    self.geometry('700x600')
  
    guiFrame.expandGrid(0,0)
   
    tipTexts = ['A table of all of the reference isotope scheme definitions available to the project',
                'A list of the residue isotopomers that comprise the selected isotope labelling scheme',
                'A three-dimensional representation of residues and their isotopomer labelling']
                
    options = ['Reference Schemes',
               'Isotopomers',
               'Isotopomer Structure']
      
    tabbedFrame = TabbedFrame(guiFrame, options=options,
                              grid=(0,0), tipTexts=tipTexts)
                              
    self.tabbedFrame = tabbedFrame
    frameA, frameB, frameC = tabbedFrame.frames
    
    #
    # Schemes
    #
    
    frameA.expandGrid(0,0)

    tipTexts = ['A short textual code that identifies the reference isotope scheme in graphical displays',
                'The full name for the isotope scheme',
                'A detailed description of the isotope scheme including user comments',
                'The name of the CCPN data repository in which the isotope scheme is saved; "refData" is in the CCPn installation']
    headingList = ['Code','Name','Description','Save Location']
    self.schemeNameEntry = Entry(self, text='', returnCallback=self.setSchemeName, width=20)
    self.schemeDetailsEntry = Entry(self, text='', returnCallback=self.setSchemeDetails, width=20)
    editWidgets      = [None, self.schemeNameEntry,
                        self.schemeDetailsEntry, None]
    editGetCallbacks = [None, self.getSchemeName,   
                        self.getSchemeDetails,   None]
    editSetCallbacks = [None, self.setSchemeName,   
                        self.setSchemeDetails,   None]
    
    self.schemeMatrix = ScrolledMatrix(frameA,
                                       headingList=headingList,
                                       callback=self.selectScheme,
                                       editWidgets=editWidgets,
                                       editSetCallbacks=editSetCallbacks,
                                       editGetCallbacks=editGetCallbacks,
                                       multiSelect=False,
                                       grid=(0,0), tipTexts=tipTexts)
    self.schemeMatrix.doEditMarkExtraRules = self.schemeEditRules
    
    tipTexts = ['Make a new reference isotope scheme definition based on a copy of the scheme currently selected',
                'Delete the selected isotope scheme',
                'Make a new, blank isotope scheme']
    texts = ['Copy','Delete','New']
    commands = [self.copyScheme, self.removeScheme, self.makeNewScheme]
    self.schemeButtons = ButtonList(frameA, texts=texts, commands=commands,
                                    grid=(1,0), tipTexts=tipTexts)
    
    #
    # Isotopomers
    #
    
    frameB.expandGrid(3,0)

    row = 0
    frame = Frame(frameB, grid=(row,0))
    frame.expandGrid(0,2)
    
    tipText = 'Selects which of the available isotope schemes to view/edit'
    label = Label(frame, text='Reference Scheme:', grid=(0,0))
    self.schemePulldown = PulldownList(frame, callback=self.setLabellingScheme,
                                       grid=(0,1), tipText=tipText)
    row += 1
    div = LabelDivider(frameB, text='Isotopomers', grid=(row,0))

    row += 1
    frame = Frame(frameB, grid=(row,0))
    frame.expandGrid(1,2)

    self.isotopomerFrame = frame
    self.abundanceWidget = MultiWidget(self, FloatEntry, relief='raised',
                                       borderwidth=2,
                                       callback=self.setDefaultAbundances,
                                       useImages=False)

    tipText = 'Opens a panel that allows you to set the basis/default abundances for C, H & N isotopes; used as the starting point for new isotopomer definitions'
    self.abundanceButton = Button(frame, text='Set Default\nAbundances',
                                  borderwidth=1, command=self.enterDefaultAbundances,
                                  grid=(0,0), tipText=tipText)

    tipText = 'Sets the basis/default abundances for C, H & N isotopes to their natural abundance proportions'
    button = Button(frame, text='Set Natural\nAbundance Default',
                    borderwidth=1, command=self.resetDefaultAbundance,
                    grid=(0,1), sticky='ew', tipText=tipText)

    label = Label(frame, text='Molecule Type:', grid=(0,2), sticky='e')
    entries = standardResidueCcpCodes.keys()
    entries.sort()
    entries.reverse()
    tipText = 'Selects which type of bio-polymer to define residue isotopomer labelling for'
    self.molTypePulldown = PulldownList(frame, callback=self.setMolType,
                                        texts=entries, grid=(0,3), tipText=tipText)

    row += 1
    tipTexts = ['The CCPN code that identifies the kind of residue the isotopomer relates to',
                'The number of the particular isotopomer (isotope pattern) within its residue type',
                'The fraction of the total residues, of its kind, that the isotopomer make up']
    headingList = ['Ccp Code','Variant','Weight']
    self.isotopomerWeightEntry = FloatEntry(self, text='', returnCallback=self.setIsotopomerWeight, width=6)
    editWidgets      = [None, None,  self.isotopomerWeightEntry]
    editGetCallbacks = [None, None,  self.getIsotopomerWeight  ]
    editSetCallbacks = [None, None,  self.setIsotopomerWeight  ]
    
    self.isotopomerMatrix = ScrolledMatrix(frameB, tipTexts=tipTexts,
                                           headingList=headingList,
                                           callback=self.selectIsotopomer,
                                           editWidgets=editWidgets,
                                           editSetCallbacks=editSetCallbacks,
                                           editGetCallbacks=editGetCallbacks,
                                           multiSelect=True, grid=(row,0))
    self.isotopomerMatrix.doEditMarkExtraRules = self.isotopomerEditRules
    
    row += 1
    frame = Frame(frameB, grid=(row,0), sticky='ew')
    frame.expandGrid(0,0)
       
    tipTexts = ['Delete the selected residue isotopomers from the current isotope scheme',
                'Make a new residue isotopomer definition by copying the details of the last selected isotopomer',
                'Add a complete set of isotopomers to the isotope scheme, one for each residue type, based on the states default isotope abundances',
                'For all residue isotopomers in the scheme, set the labelling of one kind of atom (the user is prompted) to its default isotopic incorporation ',
                'Add a new residue isotopomer definition that uses the default isotopic incorporation']
                
    texts = ['Delete\nSelected','Copy\nSelected','Add Default\nAbundance Set',
             'Set Atom Type\nTo Default','Add\nNew:']
             
    commands = [self.removeIsotopomers,
                self.duplicateResidues,
                self.addDefaultIsotopomers,
                self.setAtomTypeDefault,
                self.addNewIsotopomer]
                                
    self.isotopomerButtons = ButtonList(frame, texts=texts, commands=commands,
                                        grid=(0,0), tipTexts=tipTexts)
    tipText = 'Selects which kind of residue isotopomer may be added to the current isotope scheme'
    self.ccpCodePulldown = PulldownList(frame, callback=None, grid=(0,1),
                                        sticky='e', tipText=tipText)
    
    
    row += 1
    div = LabelDivider(frameB, text='Atom Labels', grid=(row,0))
    
    row += 1
    frame = Frame(frameB, grid=(row,0))
    frame.expandGrid(1,3)

    label = Label(frame, text='Chemical Element:', grid=(0,0))
    tipText = 'Selects which kind of atoms to select from the selected residue isotopomer; to display isotopic incorporation in the below table'
    self.elementPulldown = PulldownList(frame, callback=self.changeChemElement,
                                        grid=(0,1), tipText=tipText)
    self.updateChemElements()                                    

    label = Label(frame, text='Water Exchangeable Atoms:', grid=(0,2))
    tipText = 'Sets whether to show atoms considered as being "water exchangeable"; their isotopic labelling will rapidly equilibrate with aqueous solvent'
    self.exchangeCheck = CheckButton(frame, callback=self.updateAtomLabelsAfter,
                                      grid=(0,3), selected=False, tipText=tipText)
    row += 1
    # Tip texts set on update
    headingList = ['Atom\nName','Weighting\n13C'
                   'Weighting\n12C','%12C','%13C']
    self.atomLabelTupleWeightEntry = FloatEntry(self, text='', width=6,
                                                returnCallback=self.setAtomLabelWeight)
    
    self.atomsMatrix = ScrolledMatrix(frameB,
                                      headingList=headingList,
                                      callback=self.selectAtomLabel,
                                      multiSelect=True, grid=(row,0))
    self.atomsMatrix.doEditMarkExtraRules = self.atomsEditRules

    row += 1
    tipTexts = ['For the selected atom sites, in the current isotopomer, set their isotopic incorporation to the default values',
                'Spread the isotopic incorporation values from the last selected atom site to all selected atoms sites']
    texts = ['Reset Selected to Default Abundance','Propagate Abundances']
    commands = [self.setAtomLabelsDefault,self.propagateAbundances]
    self.atomButtons = ButtonList(frameB, texts=texts, commands=commands, 
                                  grid=(row,0), tipTexts=tipTexts)
   
    #
    # View Frame
    #
    
    frameC.expandGrid(1,0)
    
    row = 0
    frame = Frame(frameC, grid=(row,0), sticky='ew')
    frame.grid_columnconfigure(3, weight=1)
    
    label = Label(frame, text='Residue Type:', grid=(0,0))
    tipText = 'Selects which kind of residue, within the current isotope scheme, to show isotopomer structures for'
    self.viewCcpCodePulldown = PulldownList(frame, callback=self.selectViewCcpcode, 
                                            grid=(0,1), tipText=tipText)
    
    label = Label(frame, text='Isotopomer:', grid=(0,2))
    tipText = 'Selects which kind of isotopomer (labelling pattern) to display, from the selected residue type.'
    self.viewIsotopomerPulldown = PulldownList(frame, callback=self.selectViewIsotopomer,
                                               grid=(0,3), tipText=tipText)
    
    row += 1
    self.viewIsotopomerFrame = ViewIsotopomerFrame(frameC, None, grid=(row,0))
    
    #
    # Main
    #
   
    tipTexts = ['Save all changes to the reference isotope scheme to disk; the saves ALL changes to the CCPN installation for all projects to use',]
    texts    = ['Save Schemes']
    commands = [self.saveSchemes]
    self.bottomButtons = UtilityButtonList(tabbedFrame.sideFrame, texts=texts, 
                                           commands=commands, helpUrl=self.help_url,
                                           grid=(0,0), sticky = 'e', tipTexts=tipTexts)

    self.updateChemElements()
    self.updateCcpCodes()
    self.updateSchemes()
    self.administerNotifiers(self.registerNotify)

  def atomsEditRules(self, atomLabel, row, col):

    if self.scheme:
      return isSchemeEditable(self.scheme)
    
    else:
      return False

  def isotopomerEditRules(self, isotopomer, row, col):
  
    if self.scheme:
      return isSchemeEditable(self.scheme)
    
    else:
      return False

  def schemeEditRules(self, scheme, row, col):
  
    return isSchemeEditable(scheme)
  
  def administerNotifiers(self, notifyFunc):
    
    for func in ('__init__', 'delete', 'setLongName', 'setDetails'):
      for clazz in ('ccp.molecule.ChemCompLabel.LabelingScheme',):
        notifyFunc(self.updateSchemes, clazz, func)

    for func in ('__init__', 'delete', 'setWeight'):
      notifyFunc(self.updateIsotopomersAfter, 'ccp.molecule.ChemCompLabel.Isotopomer', func)
      notifyFunc(self.updateAtomLabelsAfter, 'ccp.molecule.ChemCompLabel.AtomLabel', func)
 
  def getCcpCodeIsotopomers(self, ccpCode):
    
    chemCompLabel = self.scheme.findFirstChemCompLabel(molType=self.molType, 
                                                         ccpCode=ccpCode)
      
    if chemCompLabel:
      isotopomers = list(chemCompLabel.isotopomers)
      
    else:
      isotopomers = []
   
    return isotopomers
 
  def selectViewCcpcode(self, ccpCode):
  
    if ccpCode != self.ccpCodeV:
      self.ccpCodeV = ccpCode
      self.isotopomerV = False
      self.updateViewIsotopomerPulldown()
  
  def selectViewIsotopomer(self, isotopomer):
 
    self.isotopomerV = isotopomer
    
    if isotopomer is None:
      isotopomers = self.getCcpCodeIsotopomers(self.ccpCodeV)
    else:
      isotopomers = [isotopomer,]
 
    self.viewIsotopomerFrame.setIsotopomers(isotopomers)
 
  def updateViewCcpCodePulldown(self):
  
    if self.scheme:
      codes = self.getCcpCodes(self.molType)
 
      if self.ccpCodeV not in codes:
        self.ccpCodeV = codes[0]
        self.isotopomerV = False # Not None
        self.updateViewIsotopomerPulldown()

      index = codes.index(self.ccpCodeV)
    
    else:
      codes = []
      index = 0
    
    self.viewCcpCodePulldown.setup(codes, codes, index)

    
  def updateViewIsotopomerPulldown(self):

    index = 0
    isotopomers = []
    names = []
    
    if self.scheme:
      isotopomers = self.getCcpCodeIsotopomers(self.ccpCodeV)
      names = ['%d' % i.serial for i in isotopomers]

      isotopomers.insert(0, None)
      names.insert(0, '<All>')

      if self.isotopomerV not in isotopomers:
        self.isotopomerV = None
        isotopomers = self.getCcpCodeIsotopomers(self.ccpCodeV)
        self.viewIsotopomerFrame.setIsotopomers(isotopomers)
    
    self.viewIsotopomerPulldown.setup(names, isotopomers, index)
    
  def updateButtons(self):
  
    buttonsA = self.schemeButtons.buttons
    buttonsB = self.isotopomerButtons.buttons
    buttonsC = self.atomButtons.buttons
  
    if self.scheme:
      buttonsA[0].enable()
      isEditable = isSchemeEditable(self.scheme)
      
      if isEditable:
        buttonsA[1].enable()
      else:
        buttonsA[1].disable()

      buttonsB[2].enable()
      buttonsB[3].enable()
      self.bottomButtons.buttons[0].enable()

      if isEditable:
        if self.isotopomer:
          for button in buttonsB:
            button.enable()
          for button in buttonsC:
            button.enable()
        
        else:
          buttonsB[0].disable()
          buttonsB[1].disable()
          buttonsB[3].disable()
          buttonsC[0].disable()
          buttonsC[1].disable()

          buttonsB[2].enable()
          buttonsB[4].enable()

          
      else:
        for button in buttonsB:
          button.disable()
        for button in buttonsC:
          button.disable()
      
    else:  
      buttonsA[0].disable()
      buttonsA[1].disable()
      for button in buttonsB:
        button.disable()
      for button in buttonsC:
        button.disable()
      self.bottomButtons.buttons[0].disable()
 
      
 
  def resetDefaultAbundance(self):
  
    self.defaultAbun = {}
  
  def setDefaultAbundances(self, values):

    self.abundanceWidget.place_forget()
  
    if values is not None:
      i = 0
      for element in elementSymbols: # TBD getAllIsotopes
        for code, isotope in getSortedIsotopes(self.project, element):
          self.defaultAbun[isotope] = values[i]
          i += 1
  
  def enterDefaultAbundances(self):
  
    x  = self.isotopomerFrame.winfo_x()
    y  = self.isotopomerFrame.winfo_y()
    x0 = self.abundanceButton.winfo_x()
    y0 = self.abundanceButton.winfo_y()
        
    values  = []
    options = []
    
    for element in elementSymbols: # TBD getAllIsotopes
      for code, isotope in getSortedIsotopes(self.project, element):
        options.append(code+':')
        values.append(self.defaultAbun.get(isotope, 100.0*isotope.abundance))
        
    N = len(values)    
    self.abundanceWidget.maxRows = N
    self.abundanceWidget.minRows = N
    self.abundanceWidget.set(values=values, options=options)
    self.abundanceWidget.place(x=x+x0, y=y+y0)
    
  
  def selectAtomLabel(self, obj, row, col):
  
    self.atomLabelTuple = (obj, col)
 
  def setMolType(self, molType):
  
    if molType != self.molType:
      self.molType = molType
      self.isotopomer = None
      self.updateCcpCodes()
      self.updateIsotopomers()
    
  
  def getCcpCodes(self, molType):
  
    codes = []
    for code in standardResidueCcpCodes[molType]:
      codes.append( code )
    
    codes.sort()
  
    return codes
  
  def updateCcpCodes(self):
  
    codes = self.getCcpCodes(self.molType)
    
    if self.isotopomer:
      index = codes.index(self.isotopomer.ccpCode)
    else:
      index = 0
    
    self.ccpCodePulldown.setup(codes, codes, index)

  def setIsotopomerWeight(self, event):
  
    value = self.isotopomerWeightEntry.get()
    
    if value is not None:
      self.isotopomer.setWeight( abs(value) )
    
  def getIsotopomerWeight(self, isotopomer):
  
    if isotopomer:
      self.isotopomerWeightEntry.set(isotopomer.weight)
  
  def setSchemeName(self, event):

    text = self.schemeNameEntry.get()
     
    if text:
      text = text.strip() or None
    else:
      text = None
   
    self.scheme.setLongName(text)
  
  def getSchemeName(self, scheme):
  
    if scheme:
      self.schemeNameEntry.set(scheme.longName)
  
  def getSchemeDetails(self, scheme):
  
    if scheme:
      self.schemeDetailsEntry.set(scheme.details)
  
  def setSchemeDetails(self, event):

    text = self.schemeDetailsEntry.get()
    
    if text:
      text = text.strip() or None
    else:
      text = None
    
    self.scheme.setDetails(text)
  
  def updateSchemes(self, obj=None):
  
    textMatrix = []
    objectList = []
    
    for labelingScheme in self.project.sortedLabelingSchemes():
      
      repository = labelingScheme.findFirstActiveRepository()
      
      if repository:
        saveLocation = repository.name
      else:
        saveLocation = None
      
      line = [labelingScheme.name,
              labelingScheme.longName,
              labelingScheme.details,
              saveLocation]
    
      textMatrix.append(line)
      objectList.append(labelingScheme)
    
    self.schemeMatrix.update(textMatrix=textMatrix, objectList=objectList)
    
    self.updateSchemePulldown()  
    self.updateIsotopomers()

  def updateSchemePulldown(self):
  
    scheme = self.scheme
    schemes = self.project.sortedLabelingSchemes()
    names = [ls.longName or ls.name for ls in schemes]
    
    if names:
      if scheme not in schemes:
        scheme = schemes[0]
      
      index = schemes.index(scheme)  
    
    else:
      index = 0
      scheme = None
    
    self.setLabellingScheme(scheme)
    self.schemePulldown.setup(names, schemes, index)
   
  def copyScheme(self):
  
    if self.scheme:
      name = askString('Input text','New Scheme Code:', '', parent=self)
      scheme = self.project.findFirstLabelingScheme(name=name)
      if scheme:
        showWarning('Failure','Scheme name already in use')
        return
        
      if name:
        newScheme = copySubTree(self.scheme, self.project, 
                                topObjectParameters={'name':name})
         
      else:
        showWarning('Failure','No name specified')
      
    else:
      showWarning('Failure','No scheme selected to copy')

  def removeScheme(self):
  
    if self.scheme and isSchemeEditable(self.scheme):
      self.scheme.delete()
      self.scheme = None

  def makeNewScheme(self):
  
    name = askString('Input text','New Scheme Code:', '', parent=self)

    if name:
      scheme = self.project.findFirstLabelingScheme(name=name)
      if scheme:
        showWarning('Failure','Scheme name already in use')
      else:
        scheme = self.project.newLabelingScheme(name=name)
        self.scheme = scheme
      
    else:
      showWarning('Failure','No name specified')

  def setLabellingScheme(self, scheme):

    if scheme is not self.scheme:
      self.scheme = scheme
      self.isotopomerV = False
      self.isotopomer = None
      self.updateIsotopomers()

  def selectScheme(self, object, row, col):
  
    self.setLabellingScheme(object)
    self.updateSchemePulldown()

  def open(self):

    BasePopup.open(self)

    self.updateSchemes()
          
  def saveSchemes(self):
    
    schemes = [x for x in self.project.labelingSchemes if x.isModified]
    
    if schemes:
      for scheme in schemes:
        scheme.save()
      showInfo('Notice','Successfully saved %d schemes' % len(schemes))
      self.updateSchemes()
     
    else:
      showWarning('Notice','No modified schemes to save')

  def addNewIsotopomer(self):

    if self.scheme:
      ccpCode = self.ccpCodePulldown.getObject()
      chemCompLabel = self.getChemCompLabel(self.molType, ccpCode)
         
      self.makeIsotopomer(chemCompLabel)
  
  def makeIsotopomer(self, chemCompLabel, weight=1.0):
      
      isotopomer = chemCompLabel.newIsotopomer(weight=weight)
      
      chemComp = chemCompLabel.chemComp

      for chemAtom in chemComp.chemAtoms:
        if chemAtom.elementSymbol:
          chemElement = chemAtom.chemElement
          for isotope in chemElement.isotopes:
            code = '%d%s' % (isotope.massNumber,chemAtom.elementSymbol)
            weight = self.defaultAbun.get(isotope, 100.0*isotope.abundance)
            isotopomer.newAtomLabel(name=chemAtom.name, 
                                    subType=chemAtom.subType,
                                    isotopeCode=code,
                                    weight=weight)
      
      return isotopomer
      
  def getChemCompLabel(self, molType, ccpCode):
  
    chemCompLabel = None
    if self.scheme:
      chemCompLabel = self.scheme.findFirstChemCompLabel(molType=molType, 
                                                         ccpCode=ccpCode)
      if not chemCompLabel:
        chemCompLabel = self.scheme.newChemCompLabel(molType=molType, 
                                                     ccpCode=ccpCode)
    return chemCompLabel

  def selectIsotopomer(self, obj, row, col):
  
    self.isotopomer = obj
    self.updateChemElements()
    self.updateAtomLabels()
    
  def updateIsotopomersAfter(self, obj=None):

    if self.waiting:
      return
    else:  
      self.waiting = True
      self.after_idle(self.updateIsotopomers)
 
  def updateIsotopomers(self):

    self.updateViewCcpCodePulldown()
    self.updateViewIsotopomerPulldown()

    textMatrix = []
    objectList = []
    
    if self.scheme:
      chemCompLabels = [(label.ccpCode, label) 
                        for label in self.scheme.chemCompLabels]
      chemCompLabels.sort()
      
      for key, chemCompLabel in chemCompLabels:
        if chemCompLabel.molType == self.molType:
          for isotopomer in chemCompLabel.sortedIsotopomers():
 
            line = [chemCompLabel.ccpCode, isotopomer.serial, isotopomer.weight]
 
            textMatrix.append(line)
            objectList.append(isotopomer)
    
    self.isotopomerMatrix.update(textMatrix=textMatrix, objectList=objectList)
    self.updateAtomLabelsAfter()
    self.waiting = False

  def setAtomTypeDefault(self):
  
    if self.scheme:
      atomName = askString('Query','Specify atom name to set\ndefault abundance for','H', parent=self)
    
      if not atomName:
        return
    
      atomLabels = []
      for chemCompLabel in self.scheme.chemCompLabels:
        if chemCompLabel.molType == self.molType:
          for isotopomer in chemCompLabel.isotopomers:
            # Multiple because of isotopes and subTypes
            atomLabels += isotopomer.findAllAtomLabels(name=atomName)
      
      if atomLabels:
        for atomLabel in atomLabels:
          isotope = atomLabel.isotope
          weight  = self.defaultAbun.get(isotope, 100.0*isotope.abundance)
          atomLabel.weight=weight
      
      else:
        data = (atomName,self.scheme.name)
        msg  = 'Atom name %s does not match any atoms in %s scheme isotopomers' % data
        showWarning('Failure',msg)
          

  def addDefaultIsotopomers(self):

    if self.scheme:
      codes = self.getCcpCodes(self.molType)
      
      for ccpCode in codes:
        chemCompLabel = self.getChemCompLabel(self.molType, ccpCode)
        
        if not chemCompLabel.isotopomers:
          self.makeIsotopomer(chemCompLabel)
          
  def removeIsotopomers(self):

    isotopomers = self.isotopomerMatrix.currentObjects
    
    if isotopomers:
      for isotopomer in isotopomers:
        isotopomer.delete()
     
      self.isotopomer = None
      

  def duplicateResidues(self):

    isotopomers = self.isotopomerMatrix.currentObjects

    for isotopomer in isotopomers:
      chemCompLabel = isotopomer.chemCompLabel
      new = copySubTree(isotopomer, chemCompLabel)


  def updateChemElements(self):

    
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
      
    if self.element not in names:
      index = 0
    else:
      index = names.index(self.element)

    self.elementPulldown.setup(names, names, index)

  def updateAtomLabelsAfter(self, obj=None):

    if self.waitingAtom:
      return
    
    else:
      self.waitingAtom = True
      self.after_idle(self.updateAtomLabels)


  def updateAtomLabels(self):

    element  = self.elementPulldown.getText()
    
    textMatrix = []
    objectList = []
    headingList = ['Atom Name',]
    tipTexts = ['The name of the atom within its residue, for which to set isotopic abundances, within the selected isotopomer',]
    
    isotopeCodes = [x[0] for x in getSortedIsotopes(self.project, element)]
    headingList.extend(['Weighting\n%s' % x for x in  isotopeCodes])
    tip = 'The amount of %s isotope incorporation; can be a ratio, percentage or fraction (the value used is relative to the sum of all weights)'
    tipTexts.extend([tip % x for x in isotopeCodes])
    
    tip = 'The percentage %s isotope incorporation, calculated using stated weights'
    headingList.extend(['%%%s' % x for x in  isotopeCodes])
    tipTexts.extend([tip % x for x in isotopeCodes])
    
    editWidgets      = [None,]
    editGetCallbacks = [None,]
    editSetCallbacks = [None,]
    editWidgets.extend([self.atomLabelTupleWeightEntry for x in isotopeCodes])
    editGetCallbacks.extend([self.getAtomLabelWeight for x in isotopeCodes])
    editSetCallbacks.extend([self.setAtomLabelWeight for x in isotopeCodes])
    editWidgets.extend([None for x in isotopeCodes])
    editGetCallbacks.extend([None for x in isotopeCodes])
    editSetCallbacks.extend([None for x in isotopeCodes])
    
    doExchangeable = self.exchangeCheck.get()
    
    if self.isotopomer:
      
      atomDict = {}
      for atomLabel in self.isotopomer.sortedAtomLabels():
        if atomLabel.chemAtom.elementSymbol == element:
          if (not doExchangeable) and (atomLabel.chemAtom.waterExchangeable):
            continue
            
          name    = atomLabel.name
          subType = atomLabel.subType
          atomDict[(name,subType)] = True

      atomNames = atomDict.keys()
      atomNames = greekSortAtomNames(atomNames)
      
      for atomName, subType in atomNames:
      
        if subType == 1:
          name = atomName
        else:
          name = '%s:%d' % (atomName,subType)  
      
        line = [name,]
        atomLabels = []
        
        sumWeights = 0.0
        for isotope in isotopeCodes:
          atomLabel = self.isotopomer.findFirstAtomLabel(name=atomName, 
                                                         subType=subType, 
                                                         isotopeCode=isotope)
          atomLabels.append(atomLabel)
          
          if atomLabel:
            weight = atomLabel.weight
            sumWeights += weight
            line.append(weight)
            
          else:
            line.append(0.0)
        
        if sumWeights:
          for atomLabel in atomLabels:
            line.append(100.0*atomLabel.weight/sumWeights)
          
        else:
          for atomLabel in atomLabels:
            line.append(None)
            
        textMatrix.append(line)
        objectList.append(atomLabels)
    
    
    self.atomsMatrix.update(textMatrix=textMatrix,
                            objectList=objectList,
                            headingList=headingList,
                            tipTexts=tipTexts,
                            editWidgets=editWidgets,
                            editGetCallbacks=editGetCallbacks,
                            editSetCallbacks=editSetCallbacks)

    self.waitingAtom = False
    self.updateButtons()
  
  def setAtomLabelWeight(self, event):
  
    value = self.atomLabelTupleWeightEntry.get()

    if value is not None:
      atomLabels, col = self.atomLabelTuple
      
      chemAtom = None
      for label in atomLabels:
        if label:
          chemAtom = label.chemAtom
          break
      
      atomLabel = atomLabels[col-1]
      
      if chemAtom and (atomLabel is None):
        isotopeCode, isotope = getSortedIsotopes(self.project, chemAtom.elementSymbol)[col-1]
        atomLabel = self.isotopomer.newAtomLabel(name=chemAtom.name, 
                                                 subType=chemAtom.subType,
                                                 isotopeCode=isotopeCode,
                                                 weight=value)
      atomLabel.setWeight(value)
      
  def setAtomLabelsDefault(self):
  
    atomLabelTuples = self.atomsMatrix.currentObjects
    
    for atomLabels in atomLabelTuples:
      for atomLabel in atomLabels:
        isotope = atomLabel.isotope
        weight  = self.defaultAbun.get(isotope, 100.0*isotope.abundance)
        
        atomLabel.weight = weight
        
  def propagateAbundances(self): 
  
    atomLabels, col = self.atomLabelTuple
    atomLabelTuples = self.atomsMatrix.currentObjects
  
    weightDict = {}
    for atomLabel in atomLabels:
      weightDict[atomLabel.isotope] = atomLabel.weight
    
    for atomLabels0 in atomLabelTuples:
      if atomLabels0 != atomLabels:
        for atomLabel in atomLabels0:
          atomLabel.weight = weightDict[atomLabel.isotope]
  
  def getAtomLabelWeight(self, null):
  
    atomLabels, col = self.atomLabelTuple
  
    if atomLabels and (col > 0):
      atomLabel = atomLabels[col-1]
      
      if atomLabel is None:
        weight = 0.0
      else:
        weight = atomLabel.weight  
      
      self.atomLabelTupleWeightEntry.set(weight)
    
  def changeChemElement(self, name):
    
    self.element = name        
    self.isotopes = [x[0] for x in getSortedIsotopes(self.project, name)]
    self.updateAtomLabels()

  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)  
        
    BasePopup.destroy(self)

def getSortedIsotopes(project, element, minAbundance=0.0003):
  
  chemElem = project.currentChemElementStore.findFirstChemElement(symbol=element)
  
  codes = []
  for isotope in chemElem.sortedIsotopes():
    if isotope.abundance < minAbundance:
      if isotope.spin == '0': # Look into this again!
        continue
  
    code = '%d%s' % (isotope.massNumber,element)
    codes.append((code, isotope))
  
  codes.sort()
  #codes.reverse()
  
  return codes

def getPrimaryIsotope(project, element):

  chemElem = project.currentChemElementStore.findFirstChemElement(symbol=element)
  
  bestIsotope = None
  bestAbundance = None
  
  for isotope in chemElem.sortedIsotopes():
    abundance = isotope.abundance
    
    if (bestAbundance is None) or (abundance>bestAbundance):
      bestAbundance = abundance
      bestIsotope = isotope
  
  return '%d%s' % (bestIsotope.massNumber,element)

def isSchemeEditable(scheme):
  
  if scheme.findFirstActiveRepository(name='refData'):
    return False
  else:
    return True  

if __name__ == '__main__':
  
  from memops.format.xml import Util as xmlUtil
  from memops.format.xml import XmlIO
  from memops.universal import Io as uniIo
  from memops.api import Implementation
  import Tkinter
  
  
  # standard project name
  projectName = 'refIsotopeEditor'
  
  # get or make project
  repDir = uniIo.joinPath(uniIo.normalisePath(uniIo.os.getcwd()), projectName)
  projectFile = xmlUtil.getProjectFile(repDir, projectName=projectName)
  if os.path.isfile(projectFile):
    project = XmlIO.loadProjectFile(projectFile)
  else:
    project = Implementation.MemopsRoot(name=projectName)
    
    # reorder to make sure data are saved in refData
    labelingSchemeLocator = project.findFirstPackageLocator(
     targetName='ccp.molecule.ChemCompLabel'
    )
    repositories = list(labelingSchemeLocator.repositories)
    for rr in repositories:
      if rr.name=='refData':
        refRepository = rr
        break
    else:
      refRepository = None
    if refRepository is not None:
      repositories.remove(refRepository)
      repositories.insert(0, refRepository)
      labelingSchemeLocator.repositories = repositories
      
  root = Tkinter.Tk()
  root.withdraw()
  top = IsotopeSchemeEditor(root, project)

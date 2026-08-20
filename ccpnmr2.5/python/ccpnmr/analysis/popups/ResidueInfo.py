
"""
======================COPYRIGHT/LICENSE START==========================

ResidueInfo.py: Part of the CcpNmr Analysis program

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

from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.Spacer import Spacer

from ccpnmr.analysis.core.AssignmentBasic import isResidueAssigned, getShiftLists
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.frames.ViewResidueFrame import ViewResidueFrame

class ResidueInfoPopup(BasePopup):
  """
  **Residue Assignment Information**
  
  This popup is designed to give an alternative view of residue assignment
  status to the `Atom Browser`_. The idea is that the user selects a particular
  molecular chain from the upper-left pulldown menu and then a *type* of residue
  from that chain in the adjacent pulldown. The left hand table then shows all
  the residues of that type within the selected chain's sequence.

  The left hand table lists residues of the selected type in the middle column.
  The first and last columns show the previous and next residue in the sequence
  for each central residue, thus showing a little of the sequence context for
  the selected kind of residue. If a row of this table is clicked, the central
  "i" residue is selected and the assignments are displayed in a 3D structural
  view in the right hand panel.

  The right hand panel contains a structural display for the kind of residue
  selected, using *idealised* coordinates; not coordinates from a particular
  structure. With the "Show Assignment" option selected, this 3D view has the
  chemical shifts of the assigned atoms within the selected residue superimposed
  in the view. The chemical shifts are listed after the name of each atom, but
  naturally only if it has a resonance assignment in the selected shift list.

  It should be noted that where atoms are deemed to be equivalent, like the
  three hydrogens in a fast rotating methyl atom set, the same chemical shift
  values will be used for all of the atoms within the set, although strictly
  speaking the value only really applies to the set as a whole. Also, for
  assignments that are non-stereospecific, e.g. there may be an assignment to
  Ser HBa which doesn't commit to either HB2 or HB3 specifically, the display
  will show both possible chemical shift values (should they exist) for a given
  atom. For example, a Ser HB2 atom may be labeled as "HB" 3.72/3.85" because it
  potentially relates to either the "HBa" resonance at 3.72 ppm or the "HBb"
  resonance at 3.85 ppm.

  **3D View Controls**
  
  To move and rotate the three-dimensional residue display the following
  keyboard controls may be used:
  
  * Rotate: Arrow keys
  
  * Zoom: Page Up & Page Down keys

  * Translate: Arrow keys + Control key

  Or altenatively the following mouse controls:
  
  * Rotate: Middle button click & drag
  
  * Zoom: Mouse wheel or middle button click + Shift key & drag up/down

  * Translate: Middle button click & drag + Control key

  _`Atom Browser`: BrowseAtomsPopup.html
  
  """

  def __init__(self, parent, *args, **kw):

    self.residue = None
    self.chain   = None
    self.shiftList = None

    BasePopup.__init__(self, parent=parent, title="Molecule : Residue Information", **kw)

  def open(self):
  
    self.updateChains()
    BasePopup.open(self)
    
  def body(self, guiFrame):

    self.refresh = False
    self.hilightColor   = 'lightBlue'
    self.matrixCellGrey = 'grey82'
    self.chain = None
    self.ccpCode = None
     
    row = 0
    label = Label(guiFrame, text = 'Molecule Chain:', grid=(row,0))
    tipText = 'Selects which molecular chain to select and show residues from'
    self.chainPulldown = PulldownList(guiFrame, callback=self.changeChain,
                                      grid=(row,1), tipText=tipText)
                                      
    label = Label(guiFrame, text = 'Ccp Residue Code:', grid=(row,2))
    tipText = 'Selects which kind of residue to display sequence and assignment information for'
    self.ccpCodePulldown  = PulldownList(guiFrame, callback=self.changeCcpCode,
                                         grid=(row,3), tipText=tipText)

    label = Label(guiFrame, text = 'Shift List:', grid=(row,4))
    tipText = 'Selects which shiftlist to get chemical shift values from'
    self.shiftListPulldown  = PulldownList(guiFrame, callback=self.changeShiftList,
                                            grid=(row,5), tipText=tipText)
    
    row += 1
    guiFrame.expandGrid(row, 6)
    
    self.residueFrame = LabelFrame(guiFrame, text='%s Information' % (self.ccpCode),
                                   grid=(row,0), gridSpan=(1,3))
    self.residueFrame.expandGrid(1, 0)
                                   
    self.atomsFrame= LabelFrame(guiFrame, text='Atom Information & Idealised Structure',
                                grid=(row,3), gridSpan=(1,5))
    self.atomsFrame.expandGrid(0, 0)
    
    
    self.assignedText  = 'Assigned:'
    tipText = 'How many residues of the selected kind are assigned out of the total number available'
    self.assignedLabel = Label(self.residueFrame, text=self.assignedText,
                               grid=(0,0), sticky='ew', tipText=tipText)
    
    tipTexts = ['Identity of the previous residue in the sequence',
                'Locations of the selected kind of residue, considering the selected molecular chain',
                'Identity of the next residue in the sequence']
    headingList = ['i-1','i','i+1']
    self.neighbourMatrix=ScrolledMatrix(self.residueFrame, initialRows=8,
                             headingList=headingList,  minCellWidth=6, sorting=0,
                             highlightType=1,  callback=self.selectResidue,
                             grid=(1,0), gridSpan=(1,3), tipTexts=tipTexts)
 
    self.viewResidueFrame = ViewResidueFrame(self.atomsFrame,
                                             residue=self.residue,
                                             project=self.project,
                                             grid=(0,0), tipText=tipText)
    
    tipTexts = ['Print the three-dimensional coordinate display to a PostScript, EPS or PDF file',]
    texts = [ 'Print' ]
    commands = [ self.viewResidueFrame.printStructure ]
    self.utilButtons = UtilityButtonList(guiFrame, helpUrl=self.help_url,
                                         commands=commands, texts=texts,
                                         grid=(0,7), tipTexts=tipTexts)
 
    self.updateChains()
    self.updateShiftLists()
    self.updateAfter()
    self.administerNotifiers(self.registerNotify)

  def administerNotifiers(self, notifyFunc):

    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.ResonanceSet',):
        notifyFunc(self.updateAfter, clazz, func)
 
    for func in ('__init__', 'delete'):
      for clazz in ('ccp.molecule.MolSystem.Chain',):
        notifyFunc(self.updateChains, clazz, func)

    for func in ('__init__', 'delete', 'setName'):
      for clazz in ('ccp.nmr.Nmr.ShiftList',):
        notifyFunc(self.updateShiftLists, clazz, func)
 
    notifyFunc(self.updateAfter,'ccp.molecule.MolSystem.Residue','setSeqCode')

  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)

  def selectResidue(self, residue, row, col):
  
    self.residue = residue
    self.updateAfter()

  def getChains(self, project):

    chains = []
    chainCodes = []
    for molSystem in self.project.sortedMolSystems():
      for chain in molSystem.sortedChains():
        chains.append(chain)
        chainCodes.append(self.getChainName(chain))
    
    return [chains,chainCodes]
  
  def getCcpCodes(self, chain):

    if not chain:
      return []
      
    dict  = {}
    for residue in chain.residues:
      dict[residue.molResidue.ccpCode] = residue.molResidue.molType
    
    codes = []
    codesTemp = dict.keys()
    codesTemp.sort()   
    for code in codesTemp:
      codes.append( code )
  
    return codes

  def getChainName(self, chain):
   
    return '%s:%s(%s)' % (chain.molSystem.name,chain.code,chain.molecule.molType)
 
  def changeChain(self, chain):
      
    if self.chain is not chain:
      self.chain = chain
      self.updateCcpCodes()
      self.residue = None
  
  def changeCcpCode(self, ccpCode):
  
    if self.ccpCode != ccpCode:
      self.ccpCode = ccpCode
      self.residue = None
      self.updateAfter()

  def changeShiftList(self, shiftList):
  
    if shiftList is not self.shiftList:
      self.shiftList = shiftList
      self.updateAfter()
 
  def updateChains(self, *opt):
  
    index = []
    chain = self.chain
    [chains,chainCodes] = self.getChains(self.project)
    
  
    if chains:
      if chain not in chains:
        chain = chains[0]
      
      index = chains.index(chain)  
    
    else:
      chain = None
    
    if self.chain is not chain:
      self.chain = chain
      self.updateCcpCodes()
      
    self.chainPulldown.setup(chainCodes, chains, index)

  def updateCcpCodes(self, *opt):

    index = 0
    ccpCodes = []

    if self.chain:
      ccpCodes = self.getCcpCodes(self.chain)
      
      if self.ccpCode not in ccpCodes:
        self.ccpCode = ccpCodes[0]
        self.residue = None
      
      index = ccpCodes.index(self.ccpCode)
      
      if self.residue and (self.residue.chain is not self.chain):
        self.residue = None
      
    self.ccpCodePulldown.setup(ccpCodes, ccpCodes, index) 
    self.updateAfter()


  def updateShiftLists(self, *opt):

    index = 0
    names = []
    shiftLists = getShiftLists(self.nmrProject)
    shiftList = self.shiftList

    if shiftLists:
      if shiftList not in shiftLists:
        shiftList = shiftLists[0]
      
      names = ['%d:%s' % (sl.serial, sl.name) for sl in shiftLists]
      index = shiftLists.index(shiftList)
      
    else:
      shiftList = None
    
    if shiftList is not self.shiftList:
      self.changeShiftList(shiftList)
     
    self.shiftListPulldown.setup(names, shiftLists, index) 

  def updateAfter(self, *opt):
  
    if self.refresh:
      return
    else:
      self.refresh = True
      self.after_idle(self.update)
    
  def update(self):
  
    # problem is that selectCell updates on old residue list
    # selecting a cell doesn't normally update the matrix
    if not self.chain:
      self.refresh = False
      return
    
    chain   = self.chain
    ccpCode = self.ccpCode
    sameTypeResidues = []
    seqNeighbours    = []
    neighbourData    = []
    colorData        = []
    assigned = 0

    residues = chain.sortedResidues()
        
    for residue in residues:
      if residue.molResidue.ccpCode == ccpCode:
        if self.residue is None:
          self.residue = residue
        sameTypeResidues.append(residue)
        neighbourData.append( [None,None,None] )
        colorData.append( [self.matrixCellGrey,self.matrixCellGrey,self.matrixCellGrey] )

    if self.residue and (self.residue.molResidue.molType =='protein'):
      tlc = ccpCode[:1] + ccpCode[1:].lower()
    else:
      tlc = ccpCode

    i = 0
    for residue in sameTypeResidues:
        
      j = residues.index(residue)
      neighbourData[i][1] = '%d %s' % (residue.seqCode,tlc)
      if isResidueAssigned(residue):
        assigned += 1
        colorData[i][1] = self.hilightColor

      if j >= residues[0].seqId:
        prevRes = residues[j-1]
        tlc1 = prevRes.molResidue.ccpCode
        if prevRes.molResidue.molType == 'protein':
          tlc1 = tlc1[:1] + tlc1[1:].lower()
        neighbourData[i][0] = '%d %s' % (prevRes.seqCode,tlc1)
        if isResidueAssigned(prevRes):
          colorData[i][0] = self.hilightColor

      if j+1 < residues[-1].seqId:
        nextRes = residues[j+1]
        tlc2 = nextRes.molResidue.ccpCode
        if nextRes.molResidue.molType == 'protein':
          tlc2 = tlc2[:1] + tlc2[1:].lower()
        neighbourData[i][2] = '%d %s' % (nextRes.seqCode,tlc2)
        if isResidueAssigned(nextRes):
          colorData[i][2] = self.hilightColor
      i += 1
  
    self.residueFrame.setText('%s Information' % tlc)
    self.assignedLabel.set(self.assignedText + ' %d' % assigned + ' of' + ' %d' % len(sameTypeResidues))
    self.neighbourMatrix.update(objectList=sameTypeResidues, textMatrix=neighbourData, colorMatrix=colorData)
    if self.residue:
      #self.specificLabel.set('Chain %s Residue %s %s' % (chain.code,tlc,str(self.residue.seqCode)))
      self.neighbourMatrix.hilightObject(self.residue)     
 
    self.viewResidueFrame.update(self.residue, self.shiftList)
    self.refresh = False

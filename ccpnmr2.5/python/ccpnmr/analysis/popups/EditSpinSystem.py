
"""
======================COPYRIGHT/LICENSE START==========================

EditSpinSystem.py: Part of the CcpNmr Analysis program

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
from memops.gui.CheckButton import CheckButton
from memops.gui.Frame import Frame
from memops.gui.Entry import Entry
from memops.gui.Label import Label
from memops.gui.MessageReporter import showOkCancel, showWarning, showYesNo
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.Spacer import Spacer
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.Text import Text

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.WindowBasic import getWindowPaneName, getActiveWindows, displaySpinSystemStrips
from ccpnmr.analysis.core.AssignmentBasic import newSpinSystem, assignSpinSystemType, getShiftLists, \
                                            assignSpinSystemResidue, mergeSpinSystems, getResonanceName, \
                                            clearSeqSpinSystemLinks, makeSeqSpinSystemLink, \
                                            findConnectedSpinSystem, assignTentativeSpinSystemResidues
from ccpnmr.analysis.core.MoleculeBasic import getResidueCode
 
STATUS_ANY = 'Any'
STATUS_TYPED = 'Typed Only'
STATUS_ASSIGN = 'Assigned'
STATUS_UNASSIGN = 'Unassigned'
STATUS_TENTATIVE = 'Tentative'
STATUS_ORPHAN = 'Orphaned'

PROTEIN_BACKBONE = set(['CA', 'CB', 'C', 'H', 'N', 'HA', 'HA2', 'HA3'])

class EditSpinSystemPopup(BasePopup):
  """
  **Tabular Displays of Resonance Groups**
  
  This popup window is used to display tables of the spin systems within a
  project, from these the user can perform various residue-level assignment
  and navigation operations and follow links to several kinds of relevant
  information, like peak and resonance information. In CCPN parlance a "spin
  system" is a grouping of resonances that relate to a residue or part
  residue.  Typically in an assignment project you might start with a spin
  system that consists only of amide 1H and 15N resonances and then add side
  chain resonances like CA & HA. Such a grouping of resonances relates to the
  atoms within a given residue, but the identity of the residue need not be
  known before a spin system is generated. Accordingly, a spin system may be
  completely anonymous. As more information is gained the spin system could
  be allocated a residue type, and then finally a sequence position.

  The various tabs of this popup break the spin system information into
  several sections, with the aim of avoiding clutter, but the underlying spin
  system entities represented in all the tables is the same and all carry an
  indication of the current assignment status. The differences are only a
  matter of presentation. Selections made in any of the spin system tables
  dictate which are used for the upper "Display" buttons or the lower panel
  of buttons. The four tables sub-divide the spin system information as
  follows:

  * The "Assignments" tab lists the names of the individual resonances that are
    grouped by each spin system. Also the user is able to set a name for the
    spin system, shich will be displayed on spectra in the absense of any
    residue assignment information.

  * The "Seq. Links" tab is used to display how spin systems have been
    connected as sequential neighbours. Such links are independent of a
    full residue assignment, and usually derive from the peak matching
    performed by tools like the `Protein Sequence Assignment`_ option.
    Often spin systems are connected as sequential neighbours before
    they are fitted into a protein chain to give the final residue assignment.

  * The "Shifts" tab displays the chemical shifts of the resonances that are
    contained by each spin system.  The actual values displayed come from the
    from shift list selected at the top.

  * The "Details" tab allows the user to make verbose textual comments for
    each spin system. This is particularly useful for recording thoughts that
    relate to assignment possibilities.

  **Button Functions**

  The various functions available as buttons operate on the spin system  (or
  spin systems where appropriate) that have been selected in the current
  table; using left mouse click +/- <Ctrl>/<Shift>. The "Display" functions
  at the top are designed to locate the selected spectrum window so so that it
  displays the resonances positions from the selected spin system, although
  only appropriate resonances are used according to axis/isotope type and
  spectrum width. In this way the user can find amide strips in an HNH window
  or HC strips in a HCH window. The "Display Cells" option is the same as the
  "Display Strips" option in most regards but differs in that the strips are
  sub-divided into cells, if resonances from the spin system are visible on
  that axis. This is handy to show the points of resonance intersection,
  e.g. for amino acid side chain assignment.

  The various Assign/Deassign buttons naturally control how a spin system is
  linked to residue information. Assignments may be made (or broken) to
  specific residues, residue type only or in a tentative/speculative manner.
  Although these functions are often handy, spin systems will automatically
  be assigned to the relevant residue if one of the contained resonances
  is assigned elsewhere to a specific atom in a residue; typically
  residue assignments are make when resonances are assigned to atoms
  via the `Assignment Panel`_.

  The "Show" functions are used to show tables of the spectrum peaks and
  resonances that connect to selected spin systems, and from these the user
  may then navigate within spectrum windows and perform assignment from these
  other contexts.

  **Caveats & Tips**

  If a shift appears to be missing within a spin system I could be that a
  resonance is simply not assigned in an experiment that uses the selected
  shift list, in this instance there may be a shift in a different list.

  To expedite finding particular spin systems within a large table the user
  can filter the display to only particular spin system types using the
  "Status" pulldown menu at the top. Alternatively the user can click on a
  "?" in the headings of a column to filter the table rows according to a
  user-specified value.
  
  .. _`Protein Sequence Assignment`: LinkSeqSpinSystemsPopup.html
  .. _`Assignment Panel`: EditAssignmentPopup.html
  
  """

  def __init__(self, parent, *args, **kw):

    self.guiParent = parent
    self.shiftList = None
    self.windowPane = None
                       
    BasePopup.__init__(self, parent=parent, title="Resonance : Spin Systems", **kw)

  def body(self, guiFrame):
        
    self.geometry("700x600")

    guiFrame.expandGrid(1,0)

    self.detailsEntry = Entry(self, text='', returnCallback=self.setDetails, width=12)
    self.nameEntry = Entry(self, text='', returnCallback=self.setName, width=12)
    self.detailsText = Text(self, text='',width=64,height=2)
    self.seqSegPulldown = PulldownList(self, callback=self.setSeqSeg)
      
    self.assignmentTypes = [STATUS_ANY, STATUS_TYPED,
                            STATUS_ASSIGN, STATUS_TENTATIVE,
                            STATUS_UNASSIGN, STATUS_ORPHAN] 
    self.assignmentType  = STATUS_ANY
    
    self.untoggledDict = { 'fg': '#000000', 'bg': 'grey82'}
    self.toggledDict   = { 'fg': '#FFFFFF', 'bg': '#8080C8'}

    topFrame = Frame(guiFrame)
    topFrame.grid_columnconfigure(4, weight=1)
    label = Label(topFrame, text='Status:', grid=(0,0))
    
    tipText = 'Selects which kind of spin system to show, according to its assignment status'
    self.statusPulldown = PulldownList(topFrame, callback=self.changeStatus,
                                       texts=self.assignmentTypes, grid=(0,1),
                                       tipText=tipText)
    
    self.shiftListLabel = Label(topFrame, text='Shift List:', grid=(0,2))
    tipText = 'Selects which shift list the displayed chemical shift values are from'
    self.shiftListPulldown = PulldownList(topFrame, self.changeShiftList,
                                          grid=(0,3), tipText=tipText)
    
    #label = Label(topFrame, text='Filter by Shift List:', grid=(0,4))
    text = 'Filter table rows by Shift List'
    tipText = 'Whether table resonances and rows are filtered by Shift List for Assignments and Seq. Links tabs'
    self.shiftListCheckButton = CheckButton(topFrame, text=text,
                                            callback=self.changeFilter,
                                            grid=(0,4), tipText=tipText)
    
    utilButtons = UtilityButtonList(topFrame, helpUrl=self.help_url, grid=(0,5))
        
    topFrame.grid(row=0, column=0, sticky='nsew', padx=3, pady=2)
    
    
    frame = Frame(guiFrame, grid=(2, 0))
    frame.expandGrid(None, 1)
    frame.expandGrid(None, 3)
    
    label = Label(frame, text='Assign:', grid=(0,0))
    tipTexts = ['Assign the last selected spin system, and the resonances it contains, to a specific residue in a molecular chain',
                'Assign the last selected spin system, and the resonances it contains,  in a fuzz or putative manner to a specific residue',
                'Assign the residue type of the spin system; selecting any residue of the required kind']
    texts = ['Residue','Tentative','Type']
    commands = [self.assignResidue, self.assignTentative, self.assignType]
    buttonList = ButtonList(frame, commands=commands, tipTexts=tipTexts,
                            texts=texts, grid=(0,1))
    
    for button in buttonList.buttons:
      button.config(bg='#C0D0E0')
    

    self.assignResButton = buttonList.buttons[0]
    self.assignTentativeButton = buttonList.buttons[1]
    self.assignTypButton = buttonList.buttons[2]

    label = Label(frame, text='Deassign:', grid=(0,2))
    tipTexts = ['Removes the specific residue assignment form the selected spin systems, any contained resonances will be deassigned too, but will still carry atom type information',
                'Removes any fuzzy or putative residue assignment from the selected spin systems',
                'Removes any  residue type information from the selected spin systems']
    texts = ['Residue','Tentative','Type',]
    commands = [self.deassignSeq, self.deassignTentative,
                self.deassignType]
    buttonList = ButtonList(frame, commands=commands, tipTexts=tipTexts,
                            texts=texts, grid=(0,3))
    
    for button in buttonList.buttons:
      button.config(bg='#E0D0C0')
 
    self.deassignSeqButton = buttonList.buttons[0]
    self.deassignTentativeButton = buttonList.buttons[1]
    self.deassignTypeButton = buttonList.buttons[2]

    tipTexts = ['Merge the selected spin systems into one; all contained resonances will be put in the same group and the last selected spin system will have any overriding assignment',
                'Delete the selected spin systems; only allowed if the have no resonances, show and unlink resonances if cleanup is required',
                'Show a table of the peaks assigned to the resonances of the selected spin systems',
                'Show a table of the resonances contained in the selected spin systems',
                'Open a tool that predicts the residue type of the last selected spin system, based on the chemical shift values if its resonances',
                'Make a new, blank spin system, into which resonances may be placed']
    texts    = ['Merge','Delete','Show\nPeaks', 'Show\nResonances',
                'Predict\nType','New\nSpin System']
    commands = [self.merge, self.delete, self.showPeaks, self.showResonances,
                self.predictType, self.newSpinSystem]
    self.bottomButtons = ButtonList(guiFrame, commands=commands, grid=(3,0), 
                                    texts=texts, tipTexts=tipTexts)
    self.newButton = self.bottomButtons.buttons[5]
   

    self.spinSystems = []
    self.spinSystem  = None
    self.seqSegments = []
    self.waitingForSpinSystem = False
    self.waitingForResidue = False
    self.waitingForType = False
    self.waitingForTentative = False
    
    tipTexts = ['A table listing all the spin systems within the project, including the names of all contained resonances',
                'A table listing spin systems and any sequential connections they nay have',
                'A table listing all of the chemical shift values of the spin systems within the project',
                'A table listing allowing verbose textual comments for each spin system']
    options = ['Assignments','Seq. Links','Shifts','Details',]
    self.tabbedFrame = TabbedFrame(guiFrame, options=options, grid=(1,0),
                                   callback=self.update, tipTexts=tipTexts)
    
    genFrame, seqFrame, shiftFrame, detailFrame = self.tabbedFrame.frames
    
    tipTexts = ['Use the resonances in the selected spin systems to define positions for building strip panels in the selected spectrum window',
                'Use the resonances in the selected spin systems to define positions for building resonance intersection cells in the selected spectrum window']
    texts = ['Display Strips','Display Cells']
    commands = [self.displayStrips, self.displayStripCells]
    buttonList = ButtonList(self.tabbedFrame.sideFrame, grid=(0,0),
                            texts=texts, commands=commands,
                            sticky='e', tipTexts=tipTexts)
   
    self.stripsButton, self.cellsButton = buttonList.buttons
  
    self.windowLabel = Label(self.tabbedFrame.sideFrame, text='Window:',
                             fg='#808080', sticky='e', grid=(0,1))
    
    tipText = ''
    self.stripsPulldown = PulldownList(self.tabbedFrame.sideFrame,
                                       callback=self.selectWindow,
                                       sticky='e', grid=(0,2),
                                       tipText=tipText)
    #self.stripsPulldown.inactivate()

    # General

    genFrame.grid_rowconfigure(0, weight=1)
    genFrame.grid_columnconfigure(0, weight=1)
    
    tipTexts = ['The spin system serial number',
                'The code of the molecular chain to which the spin system is assigned',
                'The sequence number and name of the residue to which the spin system is assigned',
                'A user editable name for the spin system; displayed in spectra if there are no residue assignments',
                'The complement resonances that are present within the spin system',
                "The number of peaks is the project to which the spin system's resonances are assigned"]
    colHeadings = ['#','Chain','Residue','Other\nName',
                   'Resonances','Num\nPeaks']
    editWidgets      = [None, None, None, self.nameEntry, None, None]
    editGetCallbacks = [None, None, self.predictType, self.getName,
                        self.showResonances, None]
    editSetCallbacks = [None, None, None, self.setName, None, None]
    self.generalTable = ScrolledMatrix(genFrame, tipTexts=tipTexts,
                                       multiSelect=True, 
                                       headingList=colHeadings,
                                       callback=self.selectCell,
                                       deleteFunc=self.delete,
                                       editSetCallbacks=editSetCallbacks,
                                       editGetCallbacks=editGetCallbacks,
                                       editWidgets=editWidgets)
    self.generalTable.grid(row=0, column=0, sticky='nsew')

    # Seq Conn
 
    seqFrame.grid_rowconfigure(0, weight=1)
    seqFrame.grid_columnconfigure(0, weight=1)
   
    tipTexts = ['The spin system serial number',
                'The code of the molecular chain to which the spin system is assigned',
                'The residue linked as being one position before the spin system in chain sequence',
                'The sequence number and name of the residue to which the spin system is assigned',
                'The residue linked as being one position after the spin system in chain sequence',
                'The names of the assigned backbone resonances within the spin system',
                'For sequentially connected runs of spin systems, the number of the run and the position within the run at which the spin system is found']
    colHeadings = ['#','Chain','Previous','Current','Next',
                   'Backbone\nResonances','Seq.\nSegment']
    editWidgets      = [None, None, None, None,
                        None, None, self.seqSegPulldown]
    editGetCallbacks = [None, None, None, None, None,
                        self.showResonances, self.getSeqSeg]
    editSetCallbacks = [None, None, None, None,
                        None, None, self.setSeqSeg]
    self.seqTable = ScrolledMatrix(seqFrame, tipTexts=tipTexts,
                                   multiSelect=True, 
                                   headingList=colHeadings,
                                   callback=self.selectCell,
                                   deleteFunc=self.delete,
                                   editSetCallbacks=editSetCallbacks,
                                   editGetCallbacks=editGetCallbacks,
                                   editWidgets=editWidgets)
    self.seqTable.grid(row=0, column=0, sticky='nsew')
    
    # Shifts

    shiftFrame.grid_rowconfigure(0, weight=1)
    shiftFrame.grid_columnconfigure(0, weight=1)
    
    tipTexts = ['The spin system serial number',
                'The code of the molecular chain to which the spin system is assigned',
                'The sequence number and name of the residue to which the spin system is assigned']
    self.shiftTableTipTexts = tipTexts
    colHeadings = ['#','Chain','Residue']
    editWidgets      = [None, None, None]
    editGetCallbacks = [None, None, self.predictType]
    editSetCallbacks = [None, None, None]
    self.shiftTable = ScrolledMatrix(shiftFrame, tipTexts=tipTexts,
                                     multiSelect=True, 
                                     headingList=colHeadings,
                                     callback=self.selectCell,
                                     deleteFunc=self.delete,
                                     editSetCallbacks=editSetCallbacks,
                                     editGetCallbacks=editGetCallbacks,
                                     editWidgets=editWidgets)
    self.shiftTable.grid(row=0, column=0, sticky='nsew')
    
    # Details

    detailFrame.grid_rowconfigure(0, weight=1)
    detailFrame.grid_columnconfigure(0, weight=1)

    space = ' ' * 64
    tipTexts = ['The spin system serial number',
                'The code of the molecular chain to which the spin system is assigned',
                'The sequence number and name of the residue to which the spin system is assigned',
                'A verbose user-editable textual comment for the spin system']
    colHeadings = ['#','Chain','Residue',
                   'Details' + space]
    justifyList = ['center'] * 5
    justifyList.append('left')
    editWidgets      = [None, None, None,
                        self.detailsText]
    editGetCallbacks = [None, None, self.predictType, 
                        self.getDetailsText]
    editSetCallbacks = [None, None, None,      
                         self.setDetailsText ]
    self.detailTable = ScrolledMatrix(detailFrame, tipTexts=tipTexts,
                                      multiSelect=True, 
                                      headingList=colHeadings,
                                      callback=self.selectCell,
                                      deleteFunc=self.delete,
                                      justifyList=justifyList,
                                      editSetCallbacks=editSetCallbacks,
                                      editGetCallbacks=editGetCallbacks,
                                      editWidgets=editWidgets)
    self.detailTable.grid(row=0, column=0, sticky='nsew')

    # Main

    self.refresh = False 
    self.refreshSeqSeg = False
    self.updateShiftLists()
    self.updateSeqSegmentsAfter()
    self.administerNotifiers(self.registerNotify)

  def administerNotifiers(self, notify):

    for func in ('__init__', 'delete'):
      for clazz in ('ccp.molecule.MolSystem.Chain',
                    'ccp.nmr.Nmr.ResonanceSet',
                    'ccp.nmr.Nmr.Resonance'):
        notify(self.updateAfter, clazz, func)
      notify(self.updateShiftLists, 'ccp.nmr.Nmr.ShiftList', func)
    
    notify(self.updateShiftLists, 'ccp.nmr.Nmr.ShiftList', 'setName')
    
    for func in ('setName','setAssignNames','addAssignName',
                 'removeAssignName', 'setResonanceGroup'):
      notify(self.updateResonanceAfter, 'ccp.nmr.Nmr.Resonance', func)
        
    notify(self.updateAfter,'ccp.molecule.MolSystem.Residue','setSeqCode')

    for func in ('delete','__init__','setWeight','setPossibility'):
      notify(self.updateAfter, 'ccp.nmr.Nmr.ResidueProb', func)

    for func in ('__init__', 'delete', 'setNmrChains', 'addNmrChain', 'setName',
                 'removeNmrChain', 'setResidue', 'setResonances', 'addResonance',
                 'removeResonance', 'setCcpCode', 'setDetails', 'setMolType'):
      notify(self.updateAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)

    notify(self.updateSeqSegmentsAfter, 'ccp.nmr.Nmr.ResonanceGroup', 'setResidue')
    for func in ('__init__', 'delete', 'setIsSelected', 'setSequenceOffset'):
      notify(self.updateSeqSegmentsAfter, 'ccp.nmr.Nmr.ResonanceGroupProb', func)

  def updateResonanceAfter(self, resonance):
  
    if self.refresh:
      pass
      #return
      
    elif resonance.resonanceGroup:
      self.refresh = True
      self.after_idle(self.update)
  

  def updateShiftLists(self, obj=None):
  
    index = 0
    names = []
    shiftLists = getShiftLists(self.nmrProject)
    shiftList = self.shiftList
  
    if shiftLists:
      if shiftList not in shiftLists:
        shiftList = shiftLists[0]
    
      names = ['%s [%d]' % (sl.name or '<No name>', sl.serial) for sl in shiftLists]
      index = shiftLists.index(shiftList)
    
    else:
      shiftList = None
      
    if self.shiftList is not shiftList:
      self.changeShiftList(shiftList)  
  
    self.shiftListPulldown.setup(names, shiftLists, index)
  
  def changeShiftList(self, shiftList):
  
    if shiftList is not self.shiftList:
      self.shiftList = shiftList
      if self.tabbedFrame.selected in (0, 2):
        self.update()

  def changeFilter(self, *extra):

    self.updateAfter()
  
  def displayStrips(self):
  
    windowPane = self.windowPane

    if windowPane and self.spinSystem:
      spinSystems = self.getTableSelected()
      displaySpinSystemStrips(self.guiParent, spinSystems,
                              windowPane, self.shiftList)
      self.update_idletasks()
      displaySpinSystemStrips(self.guiParent, spinSystems,
                              windowPane, self.shiftList)
      

  def displayStripCells(self):
  
    windowPane = self.windowPane
    
    if windowPane and self.spinSystem:
      displaySpinSystemStrips(self.guiParent, [self.spinSystem,],
                              windowPane, splitIntoCells=True)
      self.update_idletasks()
      displaySpinSystemStrips(self.guiParent, [self.spinSystem,],
                              windowPane, splitIntoCells=True)
      

  def updateWindows(self):
  
    index = 0
    names = []
    windowPanes = []
    windowPane = self.windowPane
    
    if self.spinSystem:
      windows  = getActiveWindows(self.spinSystem.root)
    
      if windows:
        for window in windows:
          for windowPane0 in window.sortedSpectrumWindowPanes():
            windowPanes.append(windowPane0)
            names.append(getWindowPaneName(windowPane0))
    
        if windowPane not in windowPanes:
          windowPane = windowPanes[0]

        index = windowPanes.index(windowPane)
    
        if windowPane is not self.windowPane:
          self.selectWindow(windowPane)
    
    self.stripsPulldown.setup(names, windowPanes, index)
    
  def selectWindow(self, windowPane):
  
    if windowPane is not self.windowPane:
      self.windowPane = windowPane
      # Update something?

  def showResonances(self, *null):
  
    if self.spinSystem and self.shiftList:
      resonances = self.spinSystem.resonances
      self.guiParent.viewSelectedResonances(resonances, self.shiftList)
      
  def predictType(self, *event):
  
    if self.spinSystem:
      self.guiParent.typeSpinSystem(self.spinSystem)
      

  def getSeqSegs(self, spinSystem):
    
    segId = self.getSpinSystemSeqSegment(spinSystem)
    names = ['<New>',]
    i = 0
    for segment in self.seqSegments:
      j = 0
      for spinSystem2 in segment:
        if j > 0:
          if not spinSystem2.residue: # only unassigned central positions
            names.append('%d:%3.3d' % (i,j))
          elif segId == (i,j): # unless there is no move
            names.append('%d:%3.3d' % (i,j))
 
        else:
          names.append('%d:%3.3d' % (i,j))
        
        j += 1
      if (i,j-1) != segId:  
        names.append('%d:%3.3d' % (i,j))
      
      i += 1

    return names
  
  def getSeqSeg(self, spinSystem):

    names = self.getSeqSegs(spinSystem)
    index = 0
    if names:
      segId = self.getSpinSystemSeqSegment(spinSystem)
      if segId: 
        key = '%d:%3.3d' % segId
        if key in names:
          index = names.index(key)

    self.seqSegPulldown.setup(names, names, index)

  def setSeqSeg(self, obj=None):
    
    name = self.seqSegPulldown.getObject()
      
    if self.spinSystem:
      segId = self.getSpinSystemSeqSegment(self.spinSystem)
      if name == '<New>':
        i = len(self.seqSegments)
        j = 0
      else:
        i,j   = name.split(':')
        i     = int(i)
        j     = int(j)
      
      if segId and ( segId == (i,j) ):
        return
        
      self.setSpinSystemSequenceSegment(self.spinSystem, i,j)
      
  def getDetailsText(self, spinSystem):

    if spinSystem:
      self.detailsText.setText(spinSystem.details or '')
  
  def setDetailsText(self, event):

    if self.spinSystem:
      text = self.detailsText.getText().strip()
      self.spinSystem.details = text or None
      
  def getName(self, spinSystem):

    if spinSystem:
      self.nameEntry.set(spinSystem.name or '')
  
  def setName(self, event):

    if self.spinSystem:
      text = self.nameEntry.get() or ''
      text = text.strip()
      self.spinSystem.setName(text or None)
      
  def getDetails(self, spinSystem):

    if spinSystem:
      self.detailsEntry.set(spinSystem.details or '')
  
  def setDetails(self, event):

    if self.spinSystem:
      text = self.detailsEntry.get()
      if text and text != ' ':
        self.spinSystem.setDetails( text )
   
  def newSpinSystem(self):

    self.cancelAllWaits()
    self.spinSystem = newSpinSystem(self.project)
    self.updateAfter()
    
  def changeStatus(self, status):
  
    self.assignmentType = status
    self.updateAfter()
      
  def cancelResidueWait(self):
  
    self.assignResButton.config(**self.untoggledDict)
    self.assignResButton.config(text='Assign Residue')
    self.waitingForResidue = False

  def cancelTentativeWait(self):
  
    self.assignTentativeButton.config(**self.untoggledDict)
    self.assignTentativeButton.config(text='Assign Tentative')
    self.waitingForTentative = False

  def cancelTypeWait(self):
  
    self.assignTypButton.config(**self.untoggledDict)
    self.assignTypButton.config(text='Assign Type')
    self.waitingForType = False

  def cancelAllWaits(self):
  
    if self.waitingForType:
      self.cancelTypeWait()
    
    if self.waitingForResidue:
      self.cancelResidueWait()
      
    if self.waitingForTentative:
      self.cancelTentativeWait()
   
  def deassignSeq(self):
  
    self.cancelAllWaits()
    spinSystems = self.getTableSelected()
    
    for spinSystem in spinSystems:
      assignSpinSystemResidue(spinSystem, None)
            
  def deassignTentative(self):

    self.cancelAllWaits()
    spinSystems = self.getTableSelected()
    
    for spinSystem in spinSystems:
      if spinSystem.residue:
        continue
    
      residueProbs = list(spinSystem.residueProbs)

      for residueProb in residueProbs:
        residueProb.delete()


  def deassignType(self):

    self.cancelAllWaits()
    spinSystems = self.getTableSelected()
    
    for spinSystem in spinSystems:
      if spinSystem.residue:
        continue
      
      if spinSystem.residueProbs:
        continue
      
      assignSpinSystemType(spinSystem, None)
 
  def merge(self):
  
    self.cancelAllWaits()
    spinSystems = self.getTableSelected()
    
    if len(spinSystems) == 1:
      msg = 'More than one spin system\nmust be selected for merging'
      showWarning('Merge Failed', msg, parent=self)
    
    elif len(spinSystems) > 1:
      msg = 'Do you want to merge %d spin systems' % len(spinSystems)
    
      if showOkCancel('Merge Spin System', msg, parent=self):
        for spinSystem in spinSystems[1:]:
          mergeSpinSystems(spinSystems[0],spinSystem)
  
  def delete(self, *event):

    self.cancelAllWaits()
    spinSystems = self.getTableSelected()
    
    if len(spinSystems) >0:
      toDelete = []
      assigned = 0
      for spinSystem in spinSystems:
        if len(spinSystem.resonances) > 0:
          assigned += 1
        else:
          toDelete.append(spinSystem)
    
      if showOkCancel('Delete Spin Systems','Do you want to delete %d spin systems' % len(spinSystems), parent=self):
        for spinSystem in toDelete:
          spinSystem.delete()

      if assigned > 0:
        if assigned == 1:
          showWarning('Delete failed','Cannot delete %d Spin Systems: Spin Systems still contain resonances' % assigned, parent=self)
        else: 
          showWarning('Delete failed','Cannot delete Spin System: Spin System still contain resonances', parent=self)
      
  def showPeaks(self):
  
    self.cancelAllWaits()
    spinSystems = self.getTableSelected()
    
    if len(spinSystems) > 0:
      peaksDict = {}
      for spinSystem in spinSystems:
        for resonance in spinSystem.resonances:
          for contrib in resonance.peakDimContribs:
            peaksDict[contrib.peakDim.peak] = None

      peaks = peaksDict.keys() 
      if len(peaks) > 0:
        self.guiParent.viewPeaks(peaks)
  
  def assignResidue(self):
  
    self.cancelTypeWait()
    self.guiParent.browseAtoms(requestor=self)
  
    if self.waitingForResidue:
      self.cancelResidueWait()
    else:
      self.waitingForResidue = True
      self.assignResButton.config(**self.toggledDict)  
      self.assignResButton.config(text='Cancel Assign')

  def assignType(self):
  
    self.cancelResidueWait()
    self.guiParent.browseAtoms(requestor=self)

    if self.waitingForType:
      self.cancelTypeWait()
    else:
      self.waitingForType = True
      self.assignTypButton.config(**self.toggledDict)  
      self.assignTypButton.config(text='Cancel Assign')

  def assignTentative(self):
  
    self.cancelResidueWait()
    self.guiParent.browseAtoms(requestor=self)

    if self.waitingForTentative:
      self.cancelTentativeWait()
    else:
      self.waitingForTentative = True
      self.assignTentativeButton.config(**self.toggledDict)  
      self.assignTentativeButton.config(text='Cancel Assign')
  
  def chooseResidue(self, residue):
  
    if self.waitingForResidue and self.spinSystem:
      assignSpinSystemResidue(self.spinSystem, residue, warnMerge=True)
      self.cancelResidueWait()

  def chooseTentativeResidue(self, residue):
  
    if self.waitingForTentative and self.spinSystem:
      assignTentativeSpinSystemResidues(self.spinSystem, [residue,], doWarnings=True)
      self.cancelTentativeWait()

  def chooseType(self, ccpCode, molType):
  
    if self.waitingForType and self.spinSystem:
      assignSpinSystemType(self.spinSystem, ccpCode, molType)
      self.cancelTypeWait()
  
  def selectCell(self, spinSystem, row, col):
 
    self.cancelAllWaits()
    self.spinSystem = spinSystem
    self.updateWindows()
   
    resonanceBrowser = self.guiParent.popups.get('browse_resonances')
    if resonanceBrowser:
      resonanceBrowser.chooseSpinSystem( spinSystem )

    self.updateButtons()
  
  def updateButtons(self):
  
    if self.spinSystem:
      self.assignResButton.enable()
      self.assignTypButton.enable()
      self.assignTentativeButton.enable()
      
      if self.spinSystem.residue:
        self.deassignSeqButton.enable()
        self.deassignTypeButton.disable()
        self.deassignTentativeButton.disable()
      
      elif self.spinSystem.residueProbs:
        self.deassignSeqButton.disable()
        self.deassignTentativeButton.enable()
        self.deassignTypeButton.disable()
       
      elif self.spinSystem.ccpCode :
        self.deassignSeqButton.disable()
        self.deassignTentativeButton.disable()
        self.deassignTypeButton.enable()

      else:
        self.deassignSeqButton.disable()
        self.deassignTentativeButton.disable()
        self.deassignTypeButton.disable()
              
      self.stripsButton.enable()
      self.cellsButton.enable()
      self.windowLabel.config(fg = 'black')
      #self.stripsPulldown.activate()
      for button in self.bottomButtons.buttons[:5]:
        button.enable()
    
    else:
      self.assignResButton.disable()
      self.assignTypButton.disable()
      self.assignTentativeButton.disable()
      self.deassignSeqButton.disable()
      self.deassignTentativeButton.disable()
      self.deassignTypeButton.disable()
      self.stripsButton.disable()
      self.cellsButton.disable()
      self.windowLabel.config(fg = 'grey')
      #self.stripsPulldown.inactivate()
      for button in self.bottomButtons.buttons[:5]:
        button.disable()
 
                     
  def updateAfter(self, *opt):
  
    if self.refresh:
      pass
      #return
    else:
      self.refresh = True
      self.after_idle(self.update)
  
  def _getUbiquitousInfo(self, spinSystem, status):

    residueText = None
    moleculeText = None
    residue = spinSystem.residue
    if residue:
      if status in (STATUS_TYPED, STATUS_UNASSIGN,
                    STATUS_TENTATIVE, STATUS_ORPHAN):
        return 
      
      molResidue = residue.molResidue
      chain = residue.chain
      
      residueText = '%d%s' % (residue.seqCode,getResidueCode(molResidue))
        
      if len(spinSystem.root.molSystems) > 1:
        moleculeText = '%s:%s' % (chain.molSystem.code,chain.code)
      else:
        moleculeText = '%s' % (chain.code)
        
    else:
      if status == STATUS_ASSIGN:
        return
        
      chains = spinSystem.sortedChains()
      
      if spinSystem.residueProbs:
        if status == STATUS_UNASSIGN:
          return
        
        resTexts = []
        resSeqs = []
        resCodes = set()
 
        for residueProb in spinSystem.residueProbs:
          if not residueProb.weight:
            continue
            
          residue = residueProb.possibility
          seq = residue.seqCode
          resCode = getResidueCode(residue)  
          resText = '%d?%s' % (seq, resCode)

          resTexts.append(resText)
          resSeqs.append('%d?' % seq)
          resCodes.add(resCode)
 
        if len(resCodes) == 1:
          residueText = '/'.join(resSeqs) + resCodes.pop()
        else:
          residueText = '/'.join(resTexts)
      
      elif spinSystem.ccpCode:
        if status == STATUS_UNASSIGN:
          return
        if status == STATUS_TENTATIVE:
          return
        #if status == STATUS_ORPHAN:
        #  return
        
        residueText = getResidueCode(spinSystem)
        
      elif status == STATUS_TYPED:
        return
        
      elif status == STATUS_TENTATIVE:
        return
      
      if spinSystem.resonances and status == STATUS_ORPHAN:
        return
        
      if chains:
        moleculeText = ' '.join([chain.code for chain in chains])

      elif spinSystem.molType:
        moleculeText = spinSystem.molType
    
    return (spinSystem.serial, moleculeText, residueText)

  def _getSpinSystemName(self, spinSystem):
        
    residue = spinSystem.residue
    if residue:
      tlc = getResidueCode(residue)
      
      chainCode = ''
      if len(residue.chain.molSystem.chains) > 1:
        chainCode = residue.chain.code + ' '
      
      name = '%s%d%s' % (chainCode,residue.seqCode,tlc)

    elif spinSystem.ccpCode:
      name = '{%d}%s' % (spinSystem.serial,spinSystem.ccpCode)

    elif spinSystem.name:
      name = '{%d:%s}' % (spinSystem.serial,spinSystem.name)
    
    else:
      name = '{%d}' % spinSystem.serial

    return name
    
  def getTableSelected(self):

    tableDict = {
      0: self.generalTable,
      1: self.seqTable,
      2: self.shiftTable,
      3: self.detailTable
    }

    index = self.tabbedFrame.selected
    objs = tableDict[index].currentObjects

    objs.reverse()

    return objs

  def update(self, index=None):

    self.cancelAllWaits()
    self.spinSystems = self.nmrProject.sortedResonanceGroups()
    if self.spinSystem:
      if self.spinSystem.isDeleted or (self.spinSystem.topObject is not self.nmrProject):
        self.spinSystem = None
  
    if index is None:
      index = self.tabbedFrame.selected
    #else:
    #  self.spinSystem = None

    updateFuncs = {0:self._updateGeneral,
                   1:self._updateSeq,
                   2:self._updateShifts,
                   3:self._updateDetails}

    #if index == 2:
    #  self.shiftListLabel.grid(row=0, column=2, sticky = 'w')
    #  self.shiftListPulldown.grid(row=0, column=3, sticky='w')
    #else:
    #  self.shiftListPulldown.grid_forget()
    #  self.shiftListLabel.grid_forget() 

    update = updateFuncs[index]
    update()

    self.updateButtons()

    #if self.spinSystem in objectList:
    #  self.generalTable.selectObject(self.spinSystem)
      
    self.refresh = False

  def _updateDetails(self):

    getInfo = self._getUbiquitousInfo
    status = self.assignmentType

    textMatrix = []
    objectList = [] 
    for spinSystem in self.spinSystems:
      info = getInfo(spinSystem, status)
      
      if not info:
        continue
      
      serial, molText, resText = info
          
      data = [serial,
              molText,
              resText,
              spinSystem.details]
 
      objectList.append( spinSystem )
      textMatrix.append( data )

    self.detailTable.update(objectList=objectList,
                            textMatrix=textMatrix)
                            
    if self.spinSystem: 
      self.detailTable.selectObject(self.spinSystem)
  
  def getAtomKey(self, name):
  
    if name[-1] == "'":
      return name[:2] + "'"
    else:
      return name[:2]  

    
  def _updateShifts(self):

    headingList = ['#','Chain','Residue']
    tipTexts = self.shiftTableTipTexts[:]
    
    getInfo = self._getUbiquitousInfo
    status = self.assignmentType
    shiftList = self.shiftList
    getAtomKey = self.getAtomKey
    infoDict = {}
    
    spinSystems = []
    atomTypes = set()
    for spinSystem in self.spinSystems:
      info = getInfo(spinSystem, status)
      
      if not info:
        continue
        
      resonances = []
      for resonance in spinSystem.resonances:
        shift = resonance.findFirstShift(parentList=shiftList)
        
        if not shift:
          continue

        name = getResonanceName(resonance)
        resonanceSet = resonance.resonanceSet
        if resonanceSet:
          atomTypes.add(getAtomKey(name))
        
        resonances.append((resonance, shift, name, resonanceSet))

      infoDict[spinSystem] = (resonances, info)
      spinSystems.append(spinSystem)
    
    atomTypes = list(atomTypes)
    atomTypes.sort()
    
    for atomName in ('Cg','Cb','Ca','C','N','H'):
      if atomName in atomTypes:
        atomTypes.remove(atomName)
        atomTypes.insert(0, atomName)
        
    headingList.extend(atomTypes)
    tipFormat = 'The chemical shift of the resonance(s) assigned to the %s atom(s), according to the selected shift list'
    tipTexts.extend([tipFormat % an for an in atomTypes])
    atomTypeIndices = {}
    for i, atomType in enumerate(atomTypes):
      atomTypeIndices[atomType] = i
    
    nAtoms = len(atomTypes)
    
    textMatrix = []
    objectList = [] 

    for spinSystem in spinSystems:
      resonances, info = infoDict[spinSystem]
      data = list(info)
      unassigned = []
      atomShifts = [[] for x in xrange(nAtoms)]

      for resonance, shift, name, assigned in resonances:
        if assigned:
          index = atomTypeIndices[getAtomKey(name)]
          atomShifts[index].append(shift.value)
        else:
          unassigned.append('%s:%.2f' % (name, shift.value))

      for values in atomShifts:
        values.sort()
        shiftText = ','.join(['%.2f' % v for v in values])
        data.append(shiftText or None)
  
      unassigned.sort()
      data.append(' '.join(unassigned))
           
      objectList.append( spinSystem )
      textMatrix.append( data )

    headingList.append('Unassigned')
    tipTexts.append('The chemical shifts of the resonances in the spin system not assigned to atoms')

    self.shiftTable.update(objectList=objectList,
                           textMatrix=textMatrix,
                           headingList=headingList,
                           tipTexts=tipTexts)
    if self.spinSystem: 
      self.shiftTable.selectObject(self.spinSystem)

  def _updateSeq(self):  

    getInfo = self._getUbiquitousInfo
    status = self.assignmentType
    getSpinSystemSeqSegment = self.getSpinSystemSeqSegment
    getName = self._getSpinSystemName
    
    shiftList = self.shiftList
    isFiltered = self.shiftListCheckButton.get()
    textMatrix = []
    objectList = [] 
    for spinSystem in self.spinSystems:
      info = getInfo(spinSystem, status)
      
      if not info:
        continue
      
      resonanceNames =  []
      for resonance in spinSystem.resonances:
        if not resonance.shifts:
          continue
        if isFiltered and not resonance.findFirstShift(parentList=shiftList):
          continue
      
        name = getResonanceName(resonance)
        
        if name.upper() in PROTEIN_BACKBONE:
          resonanceNames.append(name)
 
      if isFiltered and not resonanceNames:
        continue

      resonanceNames.sort()
      resonanceNames = ' '.join(resonanceNames)
        
      serial, molText, resText = info
      prevText = None
      nextText = None
      prevSS  = findConnectedSpinSystem(spinSystem, delta=-1)
      nextSS  = findConnectedSpinSystem(spinSystem, delta=1)

      if prevSS:
        prevText = getName(prevSS)
        
      if nextSS:
        nextText = getName(nextSS)
      
      seqSegText = '%d.%3.3d' % getSpinSystemSeqSegment(spinSystem)
          
      data = [serial,
              molText,
              prevText,
              resText,
              nextText,
              resonanceNames,
              seqSegText]
 
      objectList.append( spinSystem )
      textMatrix.append( data )

    self.seqTable.update(objectList=objectList,
                         textMatrix=textMatrix)

    if self.spinSystem: 
      self.seqTable.selectObject(self.spinSystem)
      
  def _updateGeneral(self):
        
    getInfo = self._getUbiquitousInfo
    status = self.assignmentType
    
    shiftList = self.shiftList
    isFiltered = self.shiftListCheckButton.get()
    textMatrix = []
    objectList = [] 
    for spinSystem in self.spinSystems:
      info = getInfo(spinSystem, status)
      
      if not info:
        continue
        
      resonanceNames =  []
      peakDict = {}
      
      for resonance in spinSystem.resonances:
        if isFiltered and not resonance.findFirstShift(parentList=shiftList):
          continue

        name = getResonanceName(resonance)
        resonanceNames.append(name)
         
        for contrib in resonance.peakDimContribs:
          peakDict[contrib.peakDim.peak] = None
  
      if isFiltered and not resonanceNames:
        continue

      resonanceNames.sort()
      nNames = len(resonanceNames)
      if nNames > 15:
        resonanceNames.insert(int(nNames/2), '\n')
      
      resonanceNames = ' '.join(resonanceNames)
          
      numPeaks = len(peakDict.keys())
        
      serial, molText, resText = info
          
      data = [serial,
              molText,
              resText,
              spinSystem.name,
              resonanceNames,
              numPeaks]
 
      objectList.append( spinSystem )
      textMatrix.append( data )
      
    self.generalTable.update(objectList=objectList,
                             textMatrix=textMatrix)
      

    if self.spinSystem: 
      self.generalTable.hilightObject(self.spinSystem)

  def getSpinSystemSeqSegment(self, spinSystem):
    
    i = 0
    
    for i, segment in enumerate(self.seqSegments):
      if spinSystem in segment:
        j = segment.index(spinSystem)
        return i, j
 
    self.seqSegments.append([spinSystem])
 
    return i, 0

  def updateSeqSegmentsAfter(self, *opt):
  
    if not self.refreshSeqSeg:
      self.refreshSeqSeg = True
      self.after_idle(self.updateSeqSegments)

  def updateSeqSegments(self):
    # called by notifier when links are made
    
    self.seqSegments = []
    if self.nmrProject:
    
      seqSegments = []
      done = {}
      for spinSystem in self.nmrProject.resonanceGroups:   
        if not done.get(spinSystem):
          segment = [spinSystem,]
          done[spinSystem] = 1
          
          spinSystemN = findConnectedSpinSystem(spinSystem, delta=-1)
          while spinSystemN and not done.get(spinSystemN):
            segment.insert(0,spinSystemN)
            done[spinSystemN] = 1
            spinSystemN = findConnectedSpinSystem(spinSystemN, delta=-1)
                   
          spinSystemC = findConnectedSpinSystem(spinSystem, delta=1)
          while spinSystemC and not done.get(spinSystemC):
            segment.append(spinSystemC)
            done[spinSystemC] = 1
            spinSystemC = findConnectedSpinSystem(spinSystemC, delta=1)

          residue = segment[0].residue
          if residue:
            chain = residue.chain
            index = 'a%s%s%5.5d' % (chain.molSystem.code, chain.code, residue.seqCode)
          elif len(segment) > 1:
            index = 'b%d%d' % (9999-len(segment),segment[0].serial)
          else:
            index = 'c%d' % segment[0].serial
           
          seqSegments.append((index,segment))
          
      seqSegments.sort()
      self.seqSegments = [x[1] for x in seqSegments]
      self.update()   
    
    self.refreshSeqSeg = False
    
  def setSpinSystemSequenceSegment(self, spinSystem, segId, position=None):

    n = len(self.seqSegments)
    segment = self.getSpinSystemSeqSegment(spinSystem)
    if segment:
      i, j = segment
      
    i2 = segId
    j2 = position or 0
    
    if i2 >= n:
      i2 = n
      self.seqSegments.append([])
   
    m = len(self.seqSegments[i2])
    if j2 >= m:
      j2 = m
      self.seqSegments[i2]
    
    if j2 <  len(self.seqSegments[i2]) and j2 > 0:
      target = self.seqSegments[i2][j2]
      if target.residue:
        idT = '%d%s' % (target.residue.seqCode,getResidueCode(target.residue))
        msg = 'Sequence placement of spin system overlaps with %s.' % idT
        if not showWarning('Failure', msg, parent=self):
          return
        
    # isolate spin system
    clearSeqSpinSystemLinks(spinSystem, delta=1)
    clearSeqSpinSystemLinks(spinSystem, delta=-1)
         
    # join any source segment ends
    if segment:
      if (j-1 >= 0) and (j+1 < len(self.seqSegments[i])):
        spinSystemN = self.seqSegments[i][j-1]
        spinSystemC = self.seqSegments[i][j+1]
        newLink = None
        if not (spinSystemN.residue and spinSystemN.residue):
          newLink = makeSeqSpinSystemLink(spinSystemN, spinSystemC, delta=1)
        if not newLink:
          # must be broken
          self.seqSegments.append( list(self.seqSegments[i][j+1:]) )
          self.seqSegments[i] = self.seqSegments[i][:j+1]
        
      self.seqSegments[i].remove(spinSystem)
    
    if i == i2:
      if j < j2:
        j2 -= 1
    
    self.seqSegments[i2].insert(j2,spinSystem)
    
    # splice into destination
    if j2 > 0:
      spinSystemN = self.seqSegments[i2][j2-1]
      clearSeqSpinSystemLinks(spinSystemN , delta=1) # clean cut
      makeSeqSpinSystemLink(spinSystemN, spinSystem, delta=1) # splice

    if j2+1 < len(self.seqSegments[i2]):
      spinSystemC = self.seqSegments[i2][j2+1]
      clearSeqSpinSystemLinks(spinSystemC , delta=-1)
      makeSeqSpinSystemLink(spinSystemC, spinSystem, delta=-1)      

    if segment:
      if not self.seqSegments[i]:
        del self.seqSegments[i]
   
  def close(self):

    self.cancelAllWaits()
    BasePopup.close(self)


  def open(self):

    self.updateAfter()
    BasePopup.open(self)

              
  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)


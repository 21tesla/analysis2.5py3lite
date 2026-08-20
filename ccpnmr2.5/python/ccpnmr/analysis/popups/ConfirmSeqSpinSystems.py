
"""
======================COPYRIGHT/LICENSE START==========================

ConfirmSeqSpinSystems.py: Part of the CcpNmr Analysis program

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

from memops.gui.ButtonList          import UtilityButtonList
from memops.gui.Label               import Label
from memops.gui.Frame               import Frame
from memops.gui.LabelFrame          import LabelFrame
from memops.gui.MessageReporter     import showYesNo
from memops.gui.PartitionedSelector import PartitionedSelector
from memops.gui.PulldownMenu        import PulldownMenu
from memops.gui.ScrolledMatrix      import ScrolledMatrix


from ccpnmr.analysis.popups.BasePopup     import BasePopup
from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes
from ccpnmr.analysis.core.MoleculeBasic   import getResidueCode
from ccpnmr.analysis.core.AssignmentBasic import assignSpinSystemResidue, makeSeqSpinSystemLink, \
                                            findConnectedSpinSystem, mergeSpinSystems



# TBD merge resonances

class ConfirmSeqSpinSystemsPopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.guiParent   = parent
    self.spinSystems = []
    self.spectrum    = None
    self.spectra     = []
    self.waiting     = 0 
    self.shiftList   = None
    self.spinSystem  = None
    self.link        = '-1'

    BasePopup.__init__(self, parent=parent, title="Confirm Sequential Spin System", **kw)


  def body(self, guiFrame):

    guiFrame.grid_columnconfigure(0, weight=1)   
    
    row = 0

    frame = Frame(guiFrame)
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(3, weight=1)
     
    label = Label(frame, text='Shift List:')
    label.grid(row=0, column=0, sticky='w')
    
    self.shiftListPulldown = PulldownMenu(frame, callback=self.setShiftList)
    self.shiftListPulldown.grid(row=0, column=1, sticky='w')

    label = Label(frame, text='Sequential Link Type:')
    label.grid(row=0, column=2, sticky='w')
    
    entries = ['-1','-1,+1','+1']
    self.linkPulldown = PulldownMenu(frame, callback=self.setLink, entries=entries, do_initial_callback=False, selected_index=entries.index(self.link))
    self.linkPulldown.grid(row=0, column=3, sticky='w')
      
    row += 1
    frame = LabelFrame(guiFrame, text='Link Atoms:')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)

    labels   = ['C','CA','CB','CG','CD','H','HA','HB','HG','HD']
    selected = ['CA','CB']
    self.atomSelector = PartitionedSelector(frame, objects=labels, labels=labels, selected=selected, toggledBg='#808080',
                                            callback=self.changeAtoms, maxRowObjects=10)
    self.atomSelector.grid(row=0, column=0, sticky='ew')

    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)

    frame = LabelFrame(guiFrame, text='Predicted Residue Assignments')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)

    headingList = ['#','Predicted\nResidue','Prob.','Links','CA','CA -1','CB','CB -1']
    self.spinSystemMatrix = ScrolledMatrix(frame, headingList=headingList, callback=self.selectSpinSystem, multiSelect=1)
    self.spinSystemMatrix.grid(row=0, column=0, sticky='nsew')

    row += 1
    texts = ['Link Selected','Link All','Commit Assignment']
    commands = [self.linkSelectedSpinSystems,self.linkAllSpinSystems,self.commitAssignments]
    buttonList = UtilityButtonList(guiFrame, texts=texts,
                                   commands=commands, helpUrl=self.help_url)
    buttonList.grid(row=row, column=0, sticky='ew')
    
    self.buttons = buttonList.buttons

    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.ShiftList',):
        self.registerNotify(self.updateShiftLists, clazz, func)

    for func in ('__init__', 'delete', 'setNmrChains', 'setResidue',
                 'setResonances','addResonance','removeResonance'):
      self.registerNotify(self.updateSpinSystemsAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)
    
    self.updateShiftLists()
  
  def setLink(self, index, name):
  
    self.link = name
    self.updateSpinSystemsAfter()
  
  
  def updateButtons(self):
  
    if len(self.spinSystemMatrix.currentObjects) > 1:
      self.buttons[0].enable()
    else:
      self.buttons[0].disable()
    
    if self.spinSystemMatrix.objectList:
      self.buttons[1].enable()
    else:
      self.buttons[1].disable()
   
    if self.spinSystem:
      self.buttons[2].enable()
    else:
      self.buttons[2].disable()
    
    
  def changeAtoms(self, *opt):
  
    self.updateSpinSystemsAfter()


  def getShiftListNames(self,shiftLists):
    
    shiftListNames = []
    for shiftList in shiftLists:
      if not hasattr(shiftList, 'name'):
        shiftList.name = "ShiftList "+ str(shiftList.serial)
      elif not shiftList.name:
        shiftList.name = "ShiftList "+ str(shiftList.serial)
      shiftListNames.append(shiftList.name)

    return shiftListNames

  def updateShiftLists(self, *opt):
  
    shiftLists     = list(self.nmrProject.findAllMeasurementLists(className = 'ShiftList'))
    shiftListNames = self.getShiftListNames(shiftLists)
    shiftList      = None
    index          = -1
    
    if shiftListNames:
      
      if self.shiftList in shiftLists:
        shiftList = self.shiftList
      else:
        shiftList = shiftLists[0]  
      
      index = shiftLists.index(shiftList)
    
    self.shiftList = shiftList
    self.shiftListPulldown.setup(shiftListNames,index)


  def setShiftList(self, index, name=None):

    shiftLists = list(self.nmrProject.findAllMeasurementLists(className = 'ShiftList'))
    if shiftLists:
      self.shiftList = shiftLists[index]
    else:
      self.shiftList = None
    
    self.updateSpinSystemsAfter()
  
  def linkSelectedSpinSystems(self):

    spinSystems = self.spinSystemMatrix.currentObjects
    self.linkSpinSystems(spinSystems)

  def linkAllSpinSystems(self):

    spinSystems = self.spinSystemMatrix.objectList
    self.linkSpinSystems(spinSystems)

  def linkSpinSystems(self, spinSystems):
    
    data = []
    
    for spinSystem in spinSystems:
      residue, p = self.getProbableResidue(spinSystem)
      key = '%s%s%4.4d' % (residue.chain.molSystem.code,residue.chain.code,residue.seqCode)
      data.append([key,residue.seqCode, residue, spinSystem])
    
    data.sort()
    seqCodes    = [x[1] for x in data]
    residues    = [x[2] for x in data]
    spinSystems = [x[3] for x in data]
    
    N = len(data)
    for i in range(N):
      
      if i > 0:
        delta = seqCodes[i] - seqCodes[i-1]
        if delta == 1 and (residues[i].chain is  residues[i-1].chain): 
          ss = findConnectedSpinSystem( spinSystems[i], delta=-1)
          if ss:
            mergeSpinSystems(ss, spinSystems[i-1]) # copy resonances across & delete
        
          makeSeqSpinSystemLink(spinSystems[i-1], spinSystems[i], delta=1)
      
      if i < N-1:
        delta = seqCodes[i+1] - seqCodes[i]
        if delta == 1 and (residues[i].chain is  residues[i+1].chain):
          ss = findConnectedSpinSystem( spinSystems[i], delta=1)
          if ss:
            mergeSpinSystems(ss, spinSystems[i+1])# copy resonances across & delete
            
          makeSeqSpinSystemLink(spinSystems[i], spinSystems[i+1], delta=1)
             
    self.updateSpinSystemsAfter()

  def commitAssignments(self):

    merge = showYesNo('Query','Merge with any existing spin systems?', self)

    for spinSystem in self.spinSystemMatrix.currentObjects:
      if not spinSystem.residue:
        residue, p = self.getProbableResidue(spinSystem)
        assignSpinSystemResidue(spinSystem, residue, warnMerge=merge)

  def getProbableResidue(self, spinSystem):

    residue     = None
    probability = 0.0

    if spinSystem.residue:
      return spinSystem.residue, 1.0

    data = []
    for residueProb in spinSystem.residueProbs:
      data.append( (residueProb.weight,residueProb.possibility) )
    data.sort()

    if data:
      residue, probability = data[-1]

    return residue, probability

  def getTentativeSpinSystems(self):

   spinSystemsList = []
   
   if self.project:
     
    for spinSystem in self.nmrProject.sortedResonanceGroups():
      if spinSystem.residueProbs and not spinSystem.residue:
      #if spinSystem.residue:
        residue, p = self.getProbableResidue(spinSystem)
        key = '%s%s%4.4d' % (residue.chain.molSystem.code,residue.chain.code,residue.seqCode)
        spinSystemsList.append( (key,spinSystem) )
   
   spinSystemsList.sort()
   
   return [x[1] for x in spinSystemsList]

  def updateSpinSystemsAfter(self, spinSystem=None):

    if self.waiting:
      return
    else:
      if spinSystem:
        if not spinSystem.residueProbs:
          return
        if spinSystem.residue:
          return
    
      self.waiting = True
      self.after_idle(self.updateSpinSystems)
  
  def getHeadings(self):
  
    headingList = ['#','Predicted\nResidue','Prob.','Links']
  
    atoms = self.atomSelector.getSelected()
  
    for atom in atoms:
      
      headingList.append( atom )
      
      if self.link == '-1':
        headingList.append( '%s -1' % atom )
        headingList.append( 'D %s -1' % atom )
      
      elif self.link == '+1':
        headingList.append( '%s +1' % atom )
        headingList.append( 'D %s +1' % atom )
      
      else:
        headingList.append( '%s -1' % atom )
        headingList.append( 'D %s -1' % atom )
        headingList.append( '%s +1' % atom )
        headingList.append( 'D %s +1' % atom )
  
    return headingList    
  
  def addShiftData(self,spinSystem,datum):
  
    prevSS, nextSS = self.getConnectedSpinSystems(spinSystem)
    atoms = self.atomSelector.getSelected()
    dict  = {}
  
    residue,p = self.getProbableResidue(spinSystem)
    prevRes   = residue.chain.findFirstResidue(seqCode=residue.seqCode-1)
    nextRes   = residue.chain.findFirstResidue(seqCode=residue.seqCode+1)
    prevSS0   = None
    nextSS0   = None
    for ss in self.getTentativeSpinSystems():
      if self.getProbableResidue(ss)[0] is prevRes:
        prevSS0 = ss
      if self.getProbableResidue(ss)[0] is nextRes:
        nextSS0 = ss
      if nextSS0 and prevSS0:
        break  
  
    for atom in atoms:
    
      resonances = []
      resonancesPrev  = []
      resonancesNext  = []
      resonancesPrev0 = []
      resonancesNext0 = []
      for resonance0 in spinSystem.sortedResonances():
        for name in resonance0.assignNames:
          if name[:2] == atom:
            resonances.append(resonance0)
            break
    
      text = ''
      if resonances:
        text = self.getResonanceText(resonances)
      datum.append(text)

      if prevSS and ('-1' in self.link):
        for resonance0 in prevSS.sortedResonances():
          for name in resonance0.assignNames:
            if name[:2] == atom:
              resonancesPrev.append(resonance0)
              break

      deltasPrev = []
      if prevSS0 and resonancesPrev:
        for resonance1 in prevSS0.sortedResonances():
          for name1 in resonance1.assignNames:
            if name1[:2] == atom:
              shift1 = resonance1.findFirstShift(parentList = self.shiftList)
              deltas = []
              for resonance2 in resonancesPrev:
                shift2 = resonance2.findFirstShift(parentList = self.shiftList)
                if shift1 and shift2:
                  deltas.append(abs(shift1.value-shift2.value))
              
              if deltas:
                deltas.sort()
                deltasPrev.append('%.2f' % deltas[0])
              break

      if nextSS and ('+1' in self.link):
        for resonance0 in nextSS.sortedResonances():
          for name in resonance0.assignNames:
            if name[:2] == atom:
              resonancesNext.append(resonance0)
              break

      deltasNext = []
      if nextSS0 and resonancesNext:
        for resonance1 in nextSS0.sortedResonances():
          for name1 in resonance1.assignNames:
            if name1[:2] == atom:
              shift1 = resonance1.findFirstShift(parentList = self.shiftList)
              deltas = []
              for resonance2 in resonancesNext:
                shift2 = resonance2.findFirstShift(parentList = self.shiftList)
                if shift1 and shift2:
                  deltas.append(abs(shift1.value-shift2.value))
              
              if deltas:
                deltas.sort()
                deltasNext.append('%.2f' % deltas[0])
              break
           
      if self.link == '-1':
        ppms = ''
        diff = ''
        if resonancesPrev:
          ppms = self.getResonanceText(resonancesPrev)
          diff = ','.join(deltasPrev)
        datum.append( ppms )
        datum.append( diff )
      
      elif self.link == '+1':
        ppms = ''
        diff = ''
        if resonancesNext:
          ppms = self.getResonanceText(resonancesNext)
          diff = ','.join(deltasNext)
        datum.append( ppms )
        datum.append( diff )
      
      else:
        ppms = ''
        diff = ''
        if resonancesPrev:
          ppms = self.getResonanceText(resonancesPrev)
          diff = ','.join(deltasPrev)
        datum.append( ppms )
        datum.append( diff )
          
        ppms = ''
        diff = ''
        if resonancesNext:
          ppms = self.getResonanceText(resonancesNext)
          diff = ','.join(deltasNext)
        datum.append( ppms )
        datum.append( diff )
  
  def updateSpinSystems(self):

    textMatrix  = []
    objectList  = []
    colorMatrix = []
    headingList = self.getHeadings()

    for spinSystem in self.getTentativeSpinSystems():
      
      residueText = None
      residue, probability  = self.getProbableResidue(spinSystem)
      if residue:
        residueText = '%d%s' % (residue.seqCode, getResidueCode(residue)) 

      links = []
      color = '#D04040'
      
      if findConnectedSpinSystem(spinSystem, delta=-1):
        links.append('-1')
      
      if findConnectedSpinSystem(spinSystem, delta=1):
        links.append('+1')

      if len(links) == 2:
        color = '#40B040'
      elif len(links) == 1:      
        color = '#B0B040'

      datum = []
      datum.append(spinSystem.serial)
      datum.append(residueText)
      datum.append(probability)
      datum.append(' '.join(links))
      
      self.addShiftData(spinSystem, datum)
      
      colors    = [None] * len(headingList)
      colors[3] = color
      
      objectList.append(spinSystem)
      textMatrix.append(datum)
      colorMatrix.append(colors)

    if self.spinSystem not in objectList:
      self.spinSystem = None

    self.spinSystemMatrix.update(headingList=headingList, objectList=objectList, textMatrix=textMatrix, colorMatrix=colorMatrix)
    self.updateButtons()
    self.waiting = False


  def getResonanceText(self, resonances):

    shifts = []
    
    for resonance in resonances:
      shift = resonance.findFirstShift(parentList=self.shiftList)
      if shift:
        shifts.append('%.2f' % shift.value)

    return ','.join(shifts)

  def getConnectedSpinSystems(self, spinSystem):

    if self.link == '-1':
      prevSS = findConnectedSpinSystem(spinSystem, delta=-1)
      nextSS = None
    elif self.link == '+1':
      prevSS = None
      nextSS = findConnectedSpinSystem(spinSystem, delta=1)
    else:
      prevSS = findConnectedSpinSystem(spinSystem, delta=-1)
      nextSS = findConnectedSpinSystem(spinSystem, delta=1)
    
    return prevSS, nextSS

  def selectSpinSystem(self, object, row, col):
  
    self.spinSystem = object
    self.updateButtons()

  def destroy(self):
  
    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.ShiftList',):
        self.unregisterNotify(self.updateShiftLists, clazz, func)

    for func in ('__init__', 'delete', 'setNmrChains', 'setResidue',
                 'setResonances', 'addResonance','removeResonance'):
      self.unregisterNotify(self.updateSpinSystemsAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)

    BasePopup.destroy(self)



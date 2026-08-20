
"""
======================COPYRIGHT/LICENSE START==========================

CalcShiftDifference.py: Part of the CcpNmr Analysis program

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

from math import sqrt

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.AssignmentBasic import makeResonanceGuiName, getResonanceResidue
from ccpnmr.analysis.core.ExperimentBasic import getOnebondDataDims
from ccpnmr.analysis.core.MoleculeBasic import getResidueCode, getChainResidueMapping

from memops.gui.ButtonList import UtilityButtonList, ButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showWarning
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

# Allow user edit of alignment

NONE_SET = set([None,])

class CalcShiftDifferencePopup(BasePopup):
  """
  **Measure Differences Between Pairs of Chemical Shifts**
  
  This system is designed to calculate the differences between pairs of chemical
  shift values in two contexts: firstly directly comparing peak locations in
  related spectra and secondly by comparing averaged chemical shift values
  (averaged over potentially many peak locations) in two different shift lists.
  The former allows the comparison to be restricted to only certain peaks, i.e.
  those with common assignments in two peak lists, and gives an overall 'shift
  distance' for each peak, rather than just differences for the separate
  dimensions.

  **Peak List Comparison**
  
  Comparing peak positions is achieved with the functions presented in the first
  tab. In essence, the user selects two peak lists from the top pulldown menus
  and views the results in the table below. The "Atom Names" options allow the
  comparison to be restricted to only certain kinds of assignment. (The atom
  names have to be complete, no wild cards are allowed.)  The "Scale Factor"
  values are used in the calculation of the "Shift Sum" and "Shift Dist" values
  so that positional differences in dimensions with dissimilar isotopes may be
  compared; the difference for a dimension is multiplied by the scale factor
  for its isotope type.
  
  **Shift List Comparison**
  
  The second tab simply takes two chemical shift lists, specified in the
  top pulldown menus and compares the values of the contained chemical shifts.
  Each row corresponds to one resonance that has a shift measurement in the
  two shift lists. The user may optionally filter the list by adding specific 
  "Atom Names" so that only certain kinds of assigned resonance are compared.
  (The atom  names have to be complete, no wild cards are allowed.) It should
  be noted that the values in shift lists are typically based upon an average
  of the assigned peak positions for resonances that use that list.  Two
  different sets of experiments may be compared by allocating them among two
  separate shift lists; see the main Experiments_ option.

  **Sequence Alignments**

  When comparing peak lists or shift lists where assignments are made to two
  different but similar chains, the calculations can optionally use a sequence
  alignment. The sequence alignment that is used is presented in the third tab
  and represents a residue to residue mapping. Which chains are used in the
  alignment is automatically determined from those that are assigned in the
  selected peak lists or shift lists and at present only the highest scoring
  pair of sequences is considered. More fine-grain control may be added in the
  future.

  **Function Buttons**

  The lower function buttons relate to either the peak list or shift list
  comparison, depending on which is currently viewed. The [Make Shift
  Difference List] is handy to give a permanent record, stored in the CCPN
  project, of the calculated chemical shift differences. The shift  difference
  lists are a kind of measurement list and are viewable at any time via the
  `Measurement Lists`_ table. If a measurement list is not made the only way to
  get the same shift differences again is to repeat the calculation. The [Show
  On Structure] is designed to colour a graphical structure display so that the
  size and colour of atoms is determined by the shift difference of their
  resonances - atoms without a shift difference are left at their default
  rendering.

  **Caveats & Tips**

  For investigating how chemical shifts differ when considering more than two
  peaks, e.g. when following how peaks move during a titration the
  `Follow Shift Changes`_ tool should be used.
  
  To make a peak list that is based on averaged chemical shift values, which may
  be compared by this system, then for certain types of spectra the user may
  create a synthetic peak list. This is accessed via the `Peak Lists`_ popup
  option [Predict from Shifts].

  .. _`Follow Shift Changes`: FollowShiftChangesPopup.html
  .. _`Peak Lists`: EditPeakListsPopup.html
  .. _Experiments: EditExperimentPopup.html
  .. _`Measurement Lists`: EditMeasurementListsPopup.html
  
  """

  def __init__(self, parent, *args, **kw):

    self.guiParent = parent
    self.peakList1 = None
    self.peakList2 = None
    self.shiftList1 = None
    self.shiftList2 = None
    self.waiting = False
    self.alignment = {}
    self.residue = None
    self.currentObject = None
    
    BasePopup.__init__(self, parent=parent, title='Data Analysis : Shift Differences', **kw)

  def body(self, guiFrame):

    self.geometry('600x700')
    
    guiFrame.expandGrid(0,0)
    
    row = 0 
    
    tipTexts = ['Measure differences in chemical shifts by comparing peak positions',
                'Measure differences in the chemical shift values from two shift lists',
                'The sequence alignment used to compare shifts between two different molecules']
    options = ['Peak List Comparison',
               'Shift List Comparison',
               'Sequence Alignment']
    tabbedFrame = TabbedFrame(guiFrame, options=options, tipTexts=tipTexts,
                              grid=(row,0), callback=self.toggleTab)
    frameA, frameB, frameC = tabbedFrame.frames
    self.tabbedFrame = tabbedFrame
    
    utilButtons = UtilityButtonList(tabbedFrame.sideFrame, grid=(0,0),
                                    sticky='e', helpUrl=self.help_url)
    
    # # # # Peak Differences Frame # # # #
    
    frame = LabelFrame(frameA, text='Options', grid=(0,0))
    frame.grid_columnconfigure(4, weight=1)

    tipText = 'First peak list to use in comparison'
    label = Label(frame, text='Peak List A:', grid=(0,0))
    self.peakList1Pulldown = PulldownList(frame, self.changePeakList1, tipText=tipText)
    self.peakList1Pulldown.grid(row=0, column=1, sticky='nw') 

    tipText = 'Second peak list to use in comparison'
    label = Label(frame, text='Peak List B:', grid=(0,2))
    self.peakList2Pulldown = PulldownList(frame, self.changePeakList2, tipText=tipText)
    self.peakList2Pulldown.grid(row=0, column=3, sticky='nw') 

    tipText = 'Whether to compare positions for only fully assigned peaks'
    label = Label(frame, text='Fully assigned only:', grid=(0,5))
    self.assignCheck = CheckButton(frame, callback=self.updateAfter, tipText=tipText,
                                   grid=(0,6), sticky='e', selected=False)
    #self.assignCheck.set(False)

    tipText = 'Restricts the comparison to peaks assigned to particular kinds of atoms, in their first dimension'
    label = Label(frame, text='Atom Names 1:', grid=(1,0))
    self.atomEntry1 = Entry(frame, text='H,H1', width=10, grid=(1,1),
                            returnCallback=self.updateAfter, tipText=tipText)

    tipText = 'Restricts the comparison to peaks assigned to particular kinds of atoms, in their second dimension'
    label = Label(frame, text='Atom Names 2:', grid=(1,2))
    self.atomEntry2 = Entry(frame, text='N', width=10, grid=(1,3),
                            returnCallback=self.updateAfter, tipText=tipText)

    tipText = 'First dimension scaling factor use to compare shifts of dissimilar isotopes'
    label = Label(frame, text='Scale factor 1:', grid=(2,0))
    self.scaleEntry1 = FloatEntry(frame, text=1.0, width=10, grid=(2,1),
                                  returnCallback=self.updateAfter, tipText=tipText)

    tipText = 'Second dimension scaling factor use to compare shifts of dissimilar isotopes'
    label = Label(frame, text='Scale factor 2:', grid=(2,2))
    self.scaleEntry2 = FloatEntry(frame, text=0.15, width=10, grid=(2,3),
                                  returnCallback=self.updateAfter, tipText=tipText)

    tipText = 'Whether to use the residue sequence alignment in the shift comparison'
    label = Label(frame, text='Use Sequence\nAlignment:', grid=(1,5), gridSpan=(2,1))
    self.alignCheck = CheckButton(frame, callback=self.updateAlignment,
                                  grid=(1,6), gridSpan=(2,1), tipText=tipText,
                                  sticky='e', selected=False)


    frameA.expandGrid(1,0)
    tipTexts = ['The residue assignment(s) of the peaks',
                'Resonance assignment in the first dimension',
                'Chemical shift of first peak (A) in first dimension',
                'Chemical shift of second peak (B) in first dimension',
                'Chemical shift difference in ppm for first dimension',
                'Resonance assignment in the second dimension',
                'Chemical shift of first peak (A) in second dimension',
                'Chemical shift of second peak (B) in second dimension',
                'Chemical shift difference in ppm for second dimension',
                'Sum of the chemical shift differences in the two dimensions',
                'The square root of the sum of the isotope weighted shift differences squared',
                'Sequence number of any assigned residue (for easy graphing)']
    headingList = ['Residue(s)',
                   'Reson.\n1','Shift\n 1A','Shift\n1B',u'\u0394 1\n(ppm)',
                   'Reson.\n2','Shift\n 2A','Shift\n2B',u'\u0394 2\n(ppm)',
                   'Shift\nSum','Shift\nDist','Seq\nNum']
    self.peakCompTable = ScrolledMatrix(frameA, headingList=headingList, grid=(1,0),
                                         callback=self.selectRow,
                                         multiSelect=True, tipTexts=tipTexts)

    # # # # Shift List Differences # # # #
    
    frameB.expandGrid(1,0)
    
    frame = Frame(frameB, grid=(0,0))
    frame.expandGrid = (0,8)
    
    tipText = 'First chemical shift list used in comparison'
    label = Label(frame, text='Shift List A:', grid=(0,0))
    self.shiftList1Pulldown = PulldownList(frame, self.changeShiftList1,
                                           grid=(0,1), tipText=tipText)

    tipText = 'Second chemical shift list used in comparison'
    label = Label(frame, text='Shift List B:', grid=(0,2))
    self.shiftList2Pulldown = PulldownList(frame, self.changeShiftList2,
                                           grid=(0,3), tipText=tipText)

    tipText = 'Whether to use the residue sequence alignment to compare dissimilar molecules'
    label = Label(frame, text='Use Sequence\nAlignment:', grid=(0,4))
    self.alignCheck2 = CheckButton(frame, callback=self.updateAlignment, tipText=tipText,
                                   grid=(0,5), sticky='w', selected=False)

    tipText = 'Restricts the comparison to only specific atom types, e.g. "HA, HB"'
    label = Label(frame, text='Atom\nNames:', grid=(0, 6))
    self.atomEntryS = Entry(frame, text='', width=10, grid=(0,7),
                            returnCallback=self.updateShiftComp, tipText=tipText)

    tipTexts = ['Serial number of first resonance',
                'Serial number of second resonance',
                'Assignment of resonance(s)',
                'Residue assigned to first resonance',
                'Residue assigned to second resonance',
                'Chemical shift value of first resonance',
                'Chemical shift value of second resonance',
                'Chemical shift difference, in units of the shift list']
    headingList = ['#1','#2','Assign\nName',
                   'Residue\n1','Residue\n2',
                   'Shift 1','Shift 2',u'\u0394']
    self.shiftCompTable = ScrolledMatrix(frameB, headingList=headingList,
                                         grid=(1,0), callback=self.selectRow,
                                         multiSelect=True, tipTexts=tipTexts)

    # # # # Alignment Frame # # # #

    frameC.expandGrid(0,0)
    tipTexts = ['Residue in the first sequence',
                'Aligned residue in the second sequence',
                'Names of atom assigned resonances common to aligned residues']
    headingList = ['Residue', 'Aligned Residue', 'Common Assignments']
    self.alignmentTable = ScrolledMatrix(frameC, headingList=headingList, grid=(0,0),
                                         callback=self.selectResidue,
                                         multiSelect=True, tipTexts=tipTexts)

    # # # # Main # # # #
    
    row +=1
    tipTexts = ['Show a table of peaks corresponding/assigned to the selected rows',
                'Force a manual refresh of the shift difference calculation',
                'Save the results in the CCPN project as a data list',
                'Show the shift difference results on a graphical structure display']
    texts = ['Show Peaks','Update','Make Shift Difference List','Show On Structure']
    commands = [self.showPeaks, self.updateAfter,
                self.makeDifferenceList, self.showStructure]
    bottomButtons = ButtonList(guiFrame, commands=commands,
                               grid=(row,0), texts=texts, tipTexts=tipTexts)
    self.showPeakButton = bottomButtons.buttons[0]
    
    self.updatePeakLists()
    self.updateShiftLists()
    
    self.administerNotifiers(self.registerNotify)


  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__','delete','setName'):
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.Experiment', func)
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.DataSource', func)
      notifyFunc(self.updateShiftLists, 'ccp.nmr.Nmr.ShiftList', func)

    for func in ('__init__','delete'):
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.PeakList', func)
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.PeakDimContrib', func)
 
    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.PeakDim', 'setPosition')
    
    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.Shift', 'setValue')
    
  def toggleTab(self, index):
  
    self.updateAfter()
  
  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

  def selectResidue(self, obj, row, col):
  
    self.residue = obj

  def selectRow(self, obj, row, col):
  
    self.currentObject = obj
    self.updateButtons()

  def makeDifferenceList(self):
  
    index = self.tabbedFrame.selected
    resonanceValues = []
    
    if index == 0:    
      data = (self.getPeakListName(self.peakList1),
              self.getPeakListName(self.peakList2))
      details = 'Peak List comparison: %s - %s' % data

      experiments = set([])
      experiments.add(self.peakList1.dataSource.experiment)
      experiments.add(self.peakList2.dataSource.experiment)
      
      for i, (data1, data2) in enumerate(self.peakCompTable.objectList):
        resonance1, values1, peaks1, peakDims1 = data1
        resonance2, values2, peaks2, peakDims2 = data2
        
        if resonance1 and (None not in values1):
          delta = values1[1]-values1[0]
          resonanceValues.append((delta, resonance1, peakDims1, peaks1))     
        
        if resonance2 and (None not in values2):
          delta = values2[1]-values2[0]
          resonanceValues.append((delta, resonance2, peakDims2, peaks2))     
    
    elif index == 1:
      textMatrix = self.shiftCompTable.textMatrix
      data = ('%s:%d' % (self.shiftList1.name, self.shiftList1.serial),
              '%s:%d' % (self.shiftList2.name, self.shiftList2.serial))
      details = 'Shift List comparison: %s - %s' % data
      
      experiments = set([])
      experiments.update(self.shiftList1.experiments)
      experiments.update(self.shiftList2.experiments)
      
      for i, (shift1, shift2) in enumerate(self.shiftCompTable.objectList):
        
        if shift1 and shift2:
          delta = shift1.value - shift2.value     
        
          peakDims = set(shift1.peakDims)
          peakDims.update(shift2.peakDims)
 
          peaks = set(shift1.peaks)
          peaks.update(shift2.peaks)
 
          resonanceValues.append((delta, shift1.resonance, peakDims, peaks))
  
    if resonanceValues:
      sdList = self.nmrProject.newShiftDifferenceList(details=details,
                                              experiments=experiments)
                                              
      sdList.name = 'ShiftDifferences%d' % (sdList.serial)
      newMeasurement = sdList.newShiftDifference
      resonances = set()
      duplicates = []
      
      for value, resonance, peakDims, peaks in resonanceValues:
        if resonance in resonances:
          # Duplicates are theoretically possible
          duplicates.append(resonance)
          continue
      
        resonances.add(resonance)
        newMeasurement(value=value, resonance=resonance,
                       peakDims=peakDims, peaks=peaks)
  
      if duplicates:
        data = ','.join(makeResonanceGuiName(r) for r in duplicates)
        msg =  'Ignoring duplicate entries for resonances: %s', data
        showWarning('Warning', msg, parent=self)
      
      self.guiParent.editMeasurements(measurementList=sdList)
      
  def showStructure(self):
  
    index = self.tabbedFrame.selected
    self.guiParent.viewStructure()
    popup = self.guiParent.popups.get('view_structure')
    popup.update_idletasks()
    if popup.structure:
      atomValues = []
    
      if index == 0:    
        textMatrix = self.peakCompTable.textMatrix
        
        for i, (data1, data2) in enumerate(self.peakCompTable.objectList):
          value = textMatrix[i][-2]
          if value is None:
            continue
          
          r1 = data1[0]
          r2 = data2[0]
          
          atoms = []
          if r1 is not None:
            rs1 = r1.resonanceSet
            
            if rs1:
              atoms.extend(rs1.findFirstAtomSet().atoms)
            
          if r2 is not None:
            rs2 = r2.resonanceSet
            
            if rs2:
              atoms.extend(rs2.findFirstAtomSet().atoms)
           
          if atoms:        
            atomValues.append((value, atoms))     
      
      elif index == 1:
        textMatrix = self.shiftCompTable.textMatrix
        
        for i, (shift1, shift2) in enumerate(self.shiftCompTable.objectList):
          value = textMatrix[i][-1]
          if value is None:
            continue
          
          r1 = shift1.resonance
          r2 = shift2.resonance
          
          rs1 = r1.resonanceSet
          rs2 = r2.resonanceSet
          
          atoms = []
          if rs1:
            atoms.extend(rs1.findFirstAtomSet().atoms)
          if rs2:
            atoms.extend(rs2.findFirstAtomSet().atoms)
          else:
            continue        
          
          atomValues.append((value, atoms))     
      
      popup.displayAtomParamsList(atomValues, size=1.0)       

  def showPeaks(self):
  
    index = self.tabbedFrame.selected
    peaks = set()

    if index == 0:
      objects = self.peakCompTable.currentObjects
      
      for data in objects:
        for datum in data:
         resonance = datum[0]
         if resonance:
           for contrib in resonance.peakDimContribs:
             peak = contrib.peakDim.peak
             
             if peak.peakList is self.peakList1:
               peaks.add(peak)
            
             elif peak.peakList is self.peakList2:
               peaks.add(peak)

    elif index == 1:
      shiftLists = set([self.shiftList1, self.shiftList2])
      for shift1, shift2 in self.shiftCompTable.currentObjects:
        for resonance in (shift1.resonance, shift2.resonance):
           for contrib in resonance.peakDimContribs:
             peak = contrib.peakDim.peak
          
             shiftList = peak.peakList.dataSource.experiment.shiftList
             if shiftList in shiftLists:
               peaks.add(peak)
    
    elif index == 2:
      if self.residue:
        for atom in self.residue.atoms:
          atomSet = atom.atomSet
        
          if atomSet:
            for resonanceSet in atomSet.resonanceSets:
              for resonance in resonanceSet.resonances:
                for contrib in resonance.peakDimContribs:
                  peak = contrib.peakDim.peak
                  peaks.add(peak)
    
    
    if peaks:
      self.guiParent.viewPeaks(list(peaks))      

     
  def getPeakListName(self, peakList):
  
    spectrum = peakList.dataSource
    return '%s:%s:%d' % (spectrum.experiment.name,spectrum.name,peakList.serial)
    
  def getPeakLists(self):
  
    names = []
    peakLists = []
    for experiment in self.nmrProject.sortedExperiments():
      for spectrum in experiment.sortedDataSources():
        if (spectrum.dataType == 'processed') and (spectrum.numDim > 1):
          for peakList in spectrum.sortedPeakLists():
            names.append( [self.getPeakListName(peakList),peakList] )
          
    return names
  
  def changePeakList1(self, peakList):
  
    if peakList is not self.peakList1:
      self.peakList1 = peakList
      self.updateAlignment()
      self.updateAfter()


  def changePeakList2(self, peakList):
  
    if peakList is not self.peakList2:
      self.peakList2 = peakList
      self.updateAlignment()
      self.updateAfter()

  def changeShiftList1(self, shiftList):
  
    if shiftList is not self.shiftList1:
      self.shiftList1 = shiftList
      self.updateAlignment()
      self.updateShiftComp()


  def changeShiftList2(self, shiftList):
  
    if shiftList is not self.shiftList2:
      self.shiftList2 = shiftList
      self.updateAlignment()
      self.updateShiftComp()
  
  def getChainPairs(self):
  
    chainPairs = []
    chains1 = set()
    chains2 = set()
    resonances1 = set()
    resonances2 = set()
    
    if self.peakList1 and self.peakList2:
 
      for peak in self.peakList1.peaks:
        for peakDim in peak.peakDims:
          for contrib in peakDim.peakDimContribs:
            resonances1.add(contrib.resonance)
      
      for peak in self.peakList2.peaks:
        for peakDim in peak.peakDims:
          for contrib in peakDim.peakDimContribs:
            resonances2.add(contrib.resonance)
    
    if not self.alignCheck.get():
      if self.shiftList1 and self.shiftList2:
        for shift in self.shiftList1.measurements:
          resonances1.add(shift.resonance)
 
        for shift in self.shiftList2.measurements:
          resonances2.add(shift.resonance)
 
    for resonance in resonances1:
      resonanceSet = resonance.resonanceSet
      
      if resonanceSet:
        atom = resonanceSet.findFirstAtomSet().findFirstAtom()
        chains1.add(atom.residue.chain)

    for resonance in resonances2:
      resonanceSet = resonance.resonanceSet
      
      if resonanceSet:
        atom = resonanceSet.findFirstAtomSet().findFirstAtom()
        chains2.add(atom.residue.chain)
    
    alignments = []
    for chain1 in chains1:
      for chain2 in chains2:
        if chain1 is chain2:
          continue
          
        mapping, score = getChainResidueMapping(chain1, chain2)
        alignments.append((score, mapping, chain1, chain2)) 
  
    if alignments:
      alignments.sort()
      alignments.reverse()
    
      done = set()
      for score, mapping, chain1, chain2 in alignments:
      
      
        if chain1 in done:
          continue
        if chain2 in done:
          continue
        
        done.add(chain1)
        done.add(chain2)
        
        chainPairs.append((chain1, chain2, mapping))
      
    return chainPairs
  
  def updateAlignment(self, *event):
  
    textMatrix = []
    objectList = []
    colorMatrix = []
    self.alignment = alignment = {}    
    
    if self.alignCheck.get() or self.alignCheck2.get():
      chainPairs = self.getChainPairs()
      
      for chainA, chainB, mapping in chainPairs:
      
        alignDict ={}
        for r1, r2 in mapping:
          if r1 and r2:
            alignDict[r1] = r2
      
        for residue1 in chainA.sortedResidues():
          residue2 = alignDict.get(residue1)
          
          residueText1 = '%d%s' % (residue1.seqCode,getResidueCode(residue1))
          
          if residue2:
            residueText2 = '%d%s' % (residue2.seqCode,getResidueCode(residue2))

            nameDict1 = {}
            nameDict2 = {}
            
            names1 = set()
            for spinSystem in residue1.resonanceGroups:
              for resonance in spinSystem.resonances:
                name = makeResonanceGuiName(resonance, fullName=False)
                names1.add(name)
                nameDict1[name] = resonance
                
            names2 = set()
            for spinSystem in residue2.resonanceGroups:
              for resonance in spinSystem.resonances:
                name = makeResonanceGuiName(resonance, fullName=False)
                names2.add(name)
                nameDict2[name] = resonance
           
            common = names1.intersection(names2)
            common = list(common)
            common.sort()
            
            # So the actual resonance mapping
            
            for name in common:
              resonance1 = nameDict1[name]
              resonance2 = nameDict2[name]
              alignment[resonance1] = resonance2
              alignment[resonance2] = resonance1
            
            commonText = ' '.join(common)
          
          else:
            residueText2 = None
            commonText = None
 
          datum = [residueText1, residueText2, commonText]
          colors = [None,None,None]
 
          textMatrix.append(datum)
          objectList.append((residue1, residue2))
          colorMatrix.append(colors)
    
    
    self.alignmentTable.update(textMatrix=textMatrix,
                               colorMatrix=colorMatrix,
                               objectList=objectList)
    self.updateAfter()
  
  def updateShiftLists(self, *obj):
  
    index1 = 0
    index2 = 0
    names = []
    mLists = self.nmrProject.sortedMeasurementLists()
    shiftLists = [ml for ml in mLists if ml.className == 'ShiftList']
    shiftList1 = self.shiftList1
    shiftList2 = self.shiftList2
    
    if shiftLists:
      names = ['%s [%d]' % (sl.name or '<No name>', sl.serial) for sl in shiftLists]
      
      if shiftList1 not in shiftLists:
        shiftList1 = shiftLists[0]
    
      if shiftList2 not in shiftLists:
        shiftList2 = shiftLists[-1]
    
      index1 = shiftLists.index(shiftList1)
      index2 = shiftLists.index(shiftList2)
    
    else:
      shiftList1 = None
      shiftList2 = None
        
    if (shiftList1 is not self.shiftList1) or (shiftList2 is not self.shiftList2):
      self.shiftList1 = shiftList1
      self.shiftList2 = shiftList2
      self.updateAlignment()
      self.updateAfter()
         
    self.shiftList1Pulldown.setup(names, shiftLists, index1)
    self.shiftList2Pulldown.setup(names, shiftLists, index2)
  
  def updatePeakLists(self, *peakList):
    
    index1 = 0
    index2 = 0
    names = []
    peakLists = []
    peakList1 = self.peakList1
    peakList2 = self.peakList2
    
    peakData = self.getPeakLists()
    if peakData:
      names = [x[0] for x in peakData]
      peakLists = [x[1] for x in peakData]
      
      if peakList1 not in peakLists:
        peakList1 = peakLists[0]
    
      if peakList2 not in peakLists:
        peakList2 = peakLists[0]
    
      index1 = peakLists.index(peakList1)
      index2 = peakLists.index(peakList2)
        
    if (peakList1 is not self.peakList1) or (peakList2 is not self.peakList2):
      self.peakList1 = peakList1
      self.peakList2 = peakList2
      self.updateAlignment()
      self.updateAfter()
         
    self.peakList1Pulldown.setup(names, peakLists, index1)
    self.peakList2Pulldown.setup(names, peakLists, index2)
        
  def updateButtons(self):
  
    if self.currentObject or self.residue:
      self.showPeakButton.enable()
    else:
      self.showPeakButton.disable()
 
  def updateAfter(self, object=None):
    
    index = self.tabbedFrame.selected
    
    if hasattr(object,'className'):
      if index == 0:
        if object.className == 'PeakDim':
          peakList = object.peak.peakList
        elif object.className == 'PeakDimContrib':
          peakList = object.peakDim.peak.peakList
        elif object.className == 'Shift':
          resonance = object.resonance
          for peakDimContrib in resonance.peakDimContribs:
            peakList = peakDimContrib.peakDim.peak.peakList
            if peakList in (self.peakList1, self.peakList2):
              break
          else:
            return
        
        if peakList not in (self.peakList1, self.peakList2):
          return
          
      if index == 1:
        if object.className == 'Shift':
          if object.parentList not in (self.shiftList1, self.shiftList2):
            return

    if self.waiting:
      return
    else:
      self.waiting = True
      
      if index == 0:
        self.after_idle(self.updatePeakComp)
      elif index == 1:
        self.after_idle(self.updateShiftComp)
      else:
        self.waiting = False
    
  def updatePeakComp(self):
   
    self.updateButtons()

    scaleFactor1 = float(self.scaleEntry1.get() or 1.0)
    scaleFactor2 = float(self.scaleEntry2.get() or 1.0)

    atomText1 = self.atomEntry1.get()
    atomText2 = self.atomEntry2.get()
    atomText1 = atomText1.replace(',',' ')
    atomText2 = atomText2.replace(',',' ')

    atomNames1 = set([s.upper() for s in atomText1.split()])
    atomNames2 = set([s.upper() for s in atomText2.split()])

    objectList  = []
    textMatrix  = []
    
    if self.peakList1 and self.peakList2:
      
      onebondDims1 = getOnebondDataDims(self.peakList1.dataSource)
      onebondDims2 = getOnebondDataDims(self.peakList2.dataSource)

      objectList   = []

      resonances1  = set()
      resonances2  = set()

      # Peak dims assigned to each resonance
      peakDimDict1 = {}
      peakDimDict2 = {}

      # Peaks assigned to found resonances
      peaks1 = {}
      peaks2 = {}
      
      assignedOnly = self.assignCheck.get()
      
      for peak in self.peakList1.sortedPeaks():
        for peakDim in peak.sortedPeakDims():
          if not peakDim.dataDimRef:
            continue
        
          for contrib in peakDim.sortedPeakDimContribs():
            resonance = contrib.resonance
            if assignedOnly and not resonance.resonanceSet:
              continue

            if peakDimDict1.get(resonance) is None:
              peakDimDict1[resonance] = []

            if resonance not in resonances1:
              peaks1[resonance] = set()
              if atomNames1:
                if resonance.resonanceSet:
                  for atomSet in resonance.resonanceSet.atomSets:
                    for atom in atomSet.atoms:
                      if atom.name in atomNames1:
                        resonances1.add(resonance)
                        peaks1[resonance].add(peak)
                        break
                else:
                  for assignName in resonance.assignNames:
                    if assignName in atomNames1:
                      resonances1.add(resonance)
                      peaks1[resonance].add(peak)
                      break
              else:
                resonances1.add(resonance)
                peaks1[resonance].add(peak)
                
                      
            else:
              peaks1[resonance].add(peak)

            if resonance not in resonances2:
              peaks2[resonance] = set()
              if atomNames2:
                if resonance.resonanceSet:
                  for atomSet in resonance.resonanceSet.atomSets:
                    for atom in atomSet.atoms:
                      if atom.name in atomNames2:
                        resonances2.add(resonance)
                        peaks2[resonance].add(peak)
                        break
                else:
                  for assignName in resonance.assignNames:
                    if assignName in atomNames2:
                      resonances2.add(resonance)
                      peaks2[resonance].add(peak)
                      break
              else:
                resonances2.add(resonance)
                peaks2[resonance].add(peak)
                
            else:
              peaks2[resonance].add(peak)
       
            peakDimDict1[resonance].append(peakDim)


      for peak in self.peakList2.peaks:
        for peakDim in peak.sortedPeakDims():
          if not peakDim.dataDimRef:
            continue
            
          for contrib in peakDim.peakDimContribs:
            resonance = contrib.resonance
            if assignedOnly and not resonance.resonanceSet:
              continue

            if peakDimDict2.get(resonance) is None:
              peakDimDict2[resonance] = []

            if resonance not in resonances1:
              peaks1[resonance] = set()
              if resonance.resonanceSet:
                for atomSet in resonance.resonanceSet.atomSets:
                  for atom in atomSet.atoms:
                    if atom.name in atomNames1:
                      resonances1.add(resonance)
                      peaks1[resonance].add(peak)
                      break
              else:
                for assignName in resonance.assignNames:
                  if assignName in atomNames1:
                    resonances1.add(resonance)
                    peaks1[resonance].add(peak)
                    break
                      
            else:
              peaks1[resonance].add(peak)

            if resonance not in resonances2:
              peaks2[resonance] = set()
              if resonance.resonanceSet:
                for atomSet in resonance.resonanceSet.atomSets:
                  for atom in atomSet.atoms:
                    if atom.name in atomNames2:
                      resonances2.add(resonance)
                      peaks2[resonance].add(peak)
                      break
              else:
                for assignName in resonance.assignNames:
                  if assignName in atomNames2:
                    resonances2.add(resonance)
                    peaks2[resonance].add(peak)
                    break
            else:
              peaks2[resonance].add(peak)
       
            peakDimDict2[resonance].append(peakDim)

      measuredPeaks1 = set([])
      measuredPeakDims1 = set([])
      measuredPeaks2 = set([])
      measuredPeakDims2 = set([])

      for resonance in peakDimDict1.keys():
        value = 0.0
        
        for peakDim in peakDimDict1[resonance]:
          value += peakDim.realValue
          measuredPeaks1.add(peakDim.peak)
          measuredPeakDims1.add(peakDim)
          
        value /= float(len(peakDimDict1[resonance] or 1.0))
        peakDimDict1[resonance] = value
        

      for resonance in peakDimDict2.keys():
        value = 0.0
        
        for peakDim in peakDimDict2[resonance]:
          value += peakDim.realValue
          measuredPeaks2.add(peakDim.peak)
          measuredPeakDims2.add(peakDim)
         
        value /= float(len(peakDimDict2[resonance] or 1.0))
        peakDimDict2[resonance] = value

      if self.alignment and self.alignCheck.get():
        alignment = self.alignment
        allowed = set(alignment.keys())
        resonances1 = resonances1.intersection(allowed)
        resonances2 = resonances2.intersection(allowed)
        
      else:
        alignment = {}
        for resonance in resonances1:
          alignment[resonance] = resonance
        for resonance in resonances2:
          alignment[resonance] = resonance
        

      paired = set()
      for resonance1 in resonances1:
        valueA = peakDimDict1.get(resonance1)
        valueB = peakDimDict2.get(alignment[resonance1])
        bound1 = resonance1.covalentlyBound
        
        resonanceA = None
        if valueA and valueB:
          resonanceA = resonance1
                
        for resonance2 in resonances2:
          
          # peaks assigned to both resonances on row
          p2 = peaks2.get(resonance2)
          if not p2:
            continue
          
          matchPeak = peaks1[resonance1].intersection(p2)
          if not matchPeak:
            continue

          skip = False
          if bound1 and matchPeak:
            if resonance2 not in bound1:
            
              # If not bound must be in same spin
              # system or residue
              if resonance1.resonanceGroup is not resonance2.resonanceGroup:
                residue1 = getResonanceResidue(resonance1)
                residue2 = getResonanceResidue(resonance2)
                if residue1 is not residue2:
                  continue
                
              # resonances have a co-assigned peak
              # but are not bound
              dataDim1 = None
              dataDim2 = None
              for peakDim in matchPeak.pop().peakDims:
                dataDimRef = peakDim.dataDimRef
                
                if not dataDimRef:
                  continue
                  
                for contrib in peakDim.peakDimContribs:
                  if contrib.resonance is resonance1:
                    dataDim1 = dataDimRef.dataDim
                  if contrib.resonance is resonance2:
                    dataDim2 = dataDimRef.dataDim
              
              if dataDim1 and dataDim2:
                # found spec dims on co-assigned peak
              
                for dataDims in onebondDims1:
                  if (dataDim1 in dataDims) and (dataDim2 in dataDims):
                    # co-assigned resonances are both in
                    # bound dims
                    skip = True
 
                for dataDims in onebondDims2:
                  if (dataDim1 in dataDims) and (dataDim2 in dataDims):
                    skip = True
                      
          if skip:
            continue

          valueC = peakDimDict1.get(resonance2)
          valueD = peakDimDict2.get(alignment[resonance2])

          resonanceB = None
          if valueC and valueD:
            resonanceB = resonance2
        
          if resonanceB and resonanceA:
            dataList = [(resonanceA, (valueA, valueB), frozenset(measuredPeaks1), frozenset(measuredPeakDims1)),
                        (resonanceB, (valueC, valueD), frozenset(measuredPeaks2), frozenset(measuredPeakDims2))]
            objectList.append(dataList)
            paired.add(resonanceA)
            paired.add(resonanceB)
        
        
      for resonance in resonances1:
        if resonance not in paired:
          valueA = peakDimDict1.get(resonance)
          valueB = peakDimDict2.get(alignment[resonance])
          if valueA and valueB:
            dataList = [(resonance, (valueA,valueB), frozenset(measuredPeaks1), frozenset(measuredPeakDims1)),
                        (None, (None, None), (), ())]
            objectList.append(dataList)

      for resonance in resonances2:
        if resonance not in paired:
          valueC = peakDimDict1.get(resonance)
          valueD = peakDimDict2.get(alignment[resonance])
          if valueC and valueD:
            dataList = [(None, (None, None), (), ()),
                        (resonance, (valueC,valueD), frozenset(measuredPeaks2), frozenset(measuredPeakDims2))]
            objectList.append(dataList)
         
            
    for data1, data2 in objectList:

      datum  = []
      nameA  = None
      nameB  = None
      valueA = None
      valueB = None
      valueC = None
      valueD = None
      delta1 = None
      delta2 = None
      residue1 = None
      residue2 = None
      
      shiftDist = None
      shiftSum  = None

      resonance1, values1, peaks1, peakDims1 = data1
      resonance2, values2, peaks2, peakDims2 = data2
      
      residues1 = set()
      residues2 = set()
      spinSystems = set()
      
      if resonance1:
        nameA = makeResonanceGuiName(resonance1, fullName=False)
        valueA, valueB = values1
        if (valueA is not None) and (valueB is not None):
          delta1 = valueB - valueA
          
        residues1.add(getResonanceResidue(resonance1))
        residues2.add(getResonanceResidue(alignment[resonance1]))
        spinSystems.add(resonance1.resonanceGroup)

      if resonance2:
        nameB = makeResonanceGuiName(resonance2, fullName=False)
        valueC, valueD = values2
        if (valueC is not None) and (valueD is not None):
          delta2 = valueD - valueC

        residues1.add(getResonanceResidue(resonance2))
        residues2.add(getResonanceResidue(alignment[resonance2]))
        spinSystems.add(resonance2.resonanceGroup)
      
      if (delta1 is not None) and (delta2 is not None):
        d1 = delta1*scaleFactor1
        d2 = delta2*scaleFactor2
        shiftDist = sqrt( (d1*d1) + (d2*d2) )
        shiftSum  = abs(d1) + abs(d2)
      #elif delta1 is not None:
      #  shiftDist = delta1*scaleFactor1

      #elif delta2 is not None:
      #  shiftDist = delta2*scaleFactor2
      
      residues1 = residues1.difference(NONE_SET)
      residues2 = residues2.difference(NONE_SET)
      spinSystems = spinSystems.difference(NONE_SET)
      
      seqNum = None
      if len(residues1) == 1:
        seqNum = tuple(residues1)[0].seqCode
        
      residueNames1 = '/'.join(['%d%s' % (r.seqCode, getResidueCode(r)) for r in residues1])
      residueNames2 = '/'.join(['%d%s' % (r.seqCode, getResidueCode(r)) for r in residues2])
      
      if residueNames1 != residueNames2:
        residueText = '%s:%s' % (residueNames1, residueNames2)
      elif not residues1:
        serials = ['{%s}' % ss.serial for ss in spinSystems]
        serials.sort()
        residueText = ','.join(serials)
      else:
        residueText = residueNames1
      
      datum = [residueText,
               nameA, valueA, valueB, delta1,
               nameB, valueC, valueD, delta2,
               shiftSum, shiftDist, seqNum]
      
      textMatrix.append( datum )
   
    self.peakCompTable.update(textMatrix=textMatrix,
                               objectList=objectList)
    
    self.waiting = False

  def updateShiftComp(self, event=None):
  
    textMatrix = []
    objectList = []
    textMatrixAppend = textMatrix.append
    objectListAppend = objectList.append
    alignment = self.alignment
    useAlignment = self.alignCheck2.get()

    atomText = self.atomEntryS.get()
    atomText = atomText.replace(',',' ')
    atomNames = set([s.upper() for s in atomText.split()])
    
    if self.shiftList1 and self.shiftList2:
      data = []
      shiftList = self.shiftList2
                       
      for shift1 in self.shiftList1.measurements:
        resonance1 = shift1.resonance
        name = makeResonanceGuiName(resonance1, fullName=False)
        
        if atomNames and (name.upper() not in atomNames):
          continue
                
        if alignment and useAlignment:
          resonance2 = alignment.get(resonance1)
          
          if not resonance2:
            continue
        
          resonanceSet = resonance1.resonanceSet
          if resonanceSet:
            residue1 = resonanceSet.findFirstAtomSet().findFirstAtom().residue
            resText1 =  '%d%s' % (residue1.seqCode, getResidueCode(residue1))
            sortKey = residue1.seqCode
          else:
            resText1 = None
            sortKey = 'z'

          resonanceSet = resonance2.resonanceSet
          if resonanceSet:
            residue2 = resonanceSet.findFirstAtomSet().findFirstAtom().residue
            resText2 =  '%d%s' % (residue2.seqCode, getResidueCode(residue2))
          else:
            resText2 = None   
        
        else:
          resonance2 = resonance1
          resonanceSet = resonance1.resonanceSet
          if resonanceSet:
            residue1 = resonanceSet.findFirstAtomSet().findFirstAtom().residue
            resText1 =  '%d%s' % (residue1.seqCode, getResidueCode(residue1))
            sortKey = residue1.seqCode
          else:
            resText1 = None
            sortKey = 'z'
          resText2 = resText1
        
        shift2 = resonance2.findFirstShift(parentList=shiftList)
        
        if not shift2:
          continue
        
        ppm1 = shift1.value
        ppm2 = shift2.value
        
        datum = [resonance1.serial, resonance2.serial, name,
                 resText1, resText2, ppm1, ppm2, ppm1-ppm2]        
        
        data.append((sortKey, name, datum, (shift1, shift2)))
      
      data.sort()
      
      for sortKey, name, datum, shifts in data:
        textMatrixAppend(datum)
        objectListAppend(shifts)
                     
    self.shiftCompTable.update(textMatrix=textMatrix,
                               objectList=objectList)

    self.waiting = False

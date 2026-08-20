
"""
======================COPYRIGHT/LICENSE START==========================

ViewQualityReports.py: Part of the CcpNmr Analysis program

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

# Option show all, or only flagged


from ccpnmr.analysis.popups.BasePopup  import BasePopup

from memops.gui.Button          import Button
from memops.gui.ButtonList      import ButtonList, UtilityButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.Entry           import Entry
from memops.gui.Frame           import Frame
from memops.gui.Label           import Label
from memops.gui.PulldownList    import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.TabbedFrame     import TabbedFrame
from memops.gui.MessageReporter import showWarning

from ccpnmr.analysis.core.QualityControlBasic import analyseAssignmentCompleteness, analyseNoeAssignments
from ccpnmr.analysis.core.QualityControlBasic import analyseChemicalShifts, analysePeaks
from ccpnmr.analysis.core.AssignmentBasic import getShiftLists, getBoundResonances, getPeakDimResonances, getAtomResonances

from ccpnmr.analysis.core.ExperimentBasic  import getThroughSpacePeakLists, getEquivalentDataDims, getOnebondExpDimRefs, getPrimaryExpDimRef

from ccpnmr.analysis.core.PeakBasic  import findSymmetryPeaks

from ccpnmr.analysis.core.MoleculeBasic import getRandomCoilShift

def qualityTestMacro(argServer):
  
  popup = ViewQualityReportsPopup(argServer.parent)

  popup.open()

ppmSdThresholds    = {'1H':0.015,'15N':0.06,'13C':0.08}
ppmDeltaThresholds = {'1H':0.035,'15N':0.16,'13C':0.18}

atomColorDict = {'H':'#a0d0a0','N':'#a0a0d0',
                 'C':'#d0d0a0','O':'#d0a0a0',
                 'P':'#d0a0d0','S':'#d0c090',
                 'Se':'#e0b0a0','F':'#d0a0d0'}

class ViewQualityReportsPopup(BasePopup):
  """
  **Analyses of Potential Errors and Completeness in Assignment**
  
  This popup presents a series of tables that enables the user to get a sense of
  how complete the assignment is within a given project and whether there are
  any assignments which look like they might be mistakes. This system is not
  designed to say that something is definitely missing or definitely a mistake,
  but rather tries to focus attention to any places where there may potentially
  be a problem. Often something that is highlighted in red as "unusual" merely
  needs checking, to ensure that the user has a good explanation. The layout
  consists of four tabs, each relating to a different aspect of assignment.

  **Atom Assign Status**
  
  The first tab contains various categories of atom, from the selected molecular
  system, like "backbone" or "Lys" (for amino acid type lysine). Each category
  has an indication of how many, out of the total available for that type, have
  been assigned and what chemical shift range they cover. This kind of
  information is useful to track the assignment progress of a project and to
  give a final summary. By default the table represents the atoms in the whole
  of the selected molecular system, but the user may specify only a subset by
  typing the appropriate chain letters and residue ranges into the "Residue
  Selection" field. Only atoms with element type
  "H", "C", "N", "P", "F" are considered.

  Explanation of categories, in terms of atom name:

  - (Backbone)N+H: "H" or "N"

  - Backbone+H+HA: "H", "HA", "HA2", "HA3", "C", "CA" or "N"

  - Backbone: "C", "CA" or "N"

  - Side Chain = not (Backbone+H+HA)

  **Through-Space Status**
  
  The second tab contains a table that gives the counts of various classes of
  through-space connection, e.g. for NOESY experiments. Counts are given on a
  per-chain and per-residue basis. In this way the user can gauge the 
  information density that has been or may be used during structure calculation.
  By default the table shows data for all peak lists that have and experiment
  type that is deemed to be "through-space", like NOESY or solid-state
  dipole-dipole. Data may be displayed for only one of these peak lists by
  selecting from the "Peak List" pulldown menu.

  **Resonances**

  The third table considers all of the resonances that carry chemical shifts in
  a specified shift list. The purpose of this table is to indicate potential
  problems with assignment from the perspective of a resonance, rather than
  peaks, which are covered in the next tab. Three basic category of error
  are highlighted (and appear as red cells):
  
  
  - Each resonance is checked for which others it has covalent bond type
    connectivity to, considering both atom assignments and the peak dimensions
    it is assigned to. Typical errors include hydrogens with more than one bound
    partner and carbons with more than four. Often errors are the result or
    (redundant) resonance duplication.

  - The chemical shift of the resonance is compared to the known shift
    distribution derived from BioMagResBank data. Outliers are coloured yellow,
    orange or red in the "Type Score" column according to how unusual
    the chemical shift is. Even red warnings may actually represent correct
    atomic assignments, but the user should be able to conform the assignment.

  - The distribution of peak dimension position information, from which the
    chemical shifts are derived, is used to fill the "SD" (standard deviation)
    and "Max Peak U+0394" columns. The former gives and indication of the width 
    of spread in the shift measurement, and is weighted by spectrum dimension,
    like all shift measurements in Analysis. The later represents the biggest
    difference between the chemical shift value (an average) and peak dimension
    position; and thus is useful for finding outliers, especially where the
    number of measurements is large and the standard deviation is small.


  The user may eliminate rows from the table which do not carry any warnings
  by selecting the "Show only alerts" option.  
  
  **Peaks**
  
  For a selected peak list, the forth table displays assignment and intensity
  problems on a per-peak basis. Two basic kinds of check are made for each peak:
  
  
  - For each peak any resonances assigned to its dimensions are analysed to
    highlight potential errors. Some of the checks mirror those found in the
    "Resonances" tab, e.g. looking for inappropriate covalent bond context.
    Checks include: missing assignments, like on one side of a covalent pair
    for an HSQC peak, and seemingly impossible assignments, like having an
    HSQC peak assigned to two resonances from different residues.

  - The intensity (height and volume) measures for each peak are compared to the
    average for the whole peak list and significant outliers are highlighted.
    For magnitude comparisons peak heights and volumes are compared in
    log-space, and red cells indicate three standard deviations from the sample
    mean of the logarithm of the intensity. Warnings may indicate the presence
    of problems, like unintentional diagonal peaks, but some very intense or
    very weak peaks may be expected. A check is also made if the peak sign seems
    inappropriate e.g. peaks of opposite sign in a NOESY experiment.


  The user may eliminate rows from the table which do not carry any warnings
  by selecting the "Show only alerts" option.
  

  **Caveats & Tips**

  Because the analyses are relatively slow to perform the tables will not
  automatically update every time there is an assignment change. If significant
  changes have been made and the tables need to be updated the user can click
  the [Update Table] button to force a refresh.

  The [Show Peaks] and [Show Resonances] buttons can be used to display all of
  the peaks and resonances, in tabular form, from which the information in the
  selected rows are derived. These options are not available for the first tab.

  Lots of red warnings in the "Resonances" tab, especially in in the "Bound"
  column, may be an indication of a problem with the specification of one or
  more experiment types, rather than a direct problem with the peak assignments.
  Knowledge of which experimental dimensions represent 'onebond' transfers, by
  virtue of peak assignment, indicates which atomic assignments are in a
  covalent bond context. The most common mistakes in this regard are the
  incorrect choice of experiment type or an incorrect choice of dimension order
  (i.e. how the real experiment dimensions map to the CCPN reference). Both of
  these problems may be fixed by considering the "Experiment Types" tab in the
  Experiments_ popup window. Here, most dimension order problems are fixed by
  editing the "Transfers To Dim" column in the lower left table, using the
  user's knowledge of how the dimensions are connected/correlated.

  .. _Experiments: EditExperimentPopup.html

  """

  def __init__(self, parent, *args, **kw):

    self.guiParent  = parent
    self.molSystemA = None
    self.molSystemB = None
    self.molSystemD = None
    self.shiftListA = None
    self.shiftListC = None
    self.peakListB  = None
    self.peakListD  = None
    self.object     = None
    
    BasePopup.__init__(self, parent=parent, title='Assignment : Quality Reports', **kw)
    
  def open(self):
  
    self.updateMolSystems()
    self.updateShiftLists()
    self.updatePeakLists()
  
    BasePopup.open(self)

  def close(self):
  
    BasePopup.close(self)

 
  def body(self, guiFrame):
  

    self.geometry('600x700')

    guiFrame.expandGrid(0,0)
    
    tipTexts = ['',
                '',
                '',
                '']
    options = ['Atom Assign Status','Through-Space Status',
               'Resonances','Peaks']
               
    tabbedFrame = TabbedFrame(guiFrame, options=options, callback=self.selectTab,
                              grid=(0,0), tipTexts=tipTexts)
    frameA, frameB, frameC, frameD = tabbedFrame.frames
    frameA.grid_columnconfigure(0, weight=1)
    frameB.grid_columnconfigure(0, weight=1)
    frameC.grid_columnconfigure(0, weight=1)
    frameD.grid_columnconfigure(0, weight=1)
    self.tabbedFrame = tabbedFrame
   
    #
    # Atom 
    #

    frameA.grid_rowconfigure(2, weight=1)
    
    frame = Frame(frameA, grid=(0,0),sticky='ew')
    frame.grid_columnconfigure(4, weight=1)
    
    label = Label(frame, text='Mol System:', grid=(0,0))
    tipText = 'Selects which molecular system (group of chains) to display assignment status information for'
    self.molSystemPulldownA = PulldownList(frame, grid=(0,1), tipText=tipText,
                                           callback=self.changeMolSystemA)
    
    label = Label(frame, text='Shift List:', grid=(0,2))
    tipText = 'Selects which shift list is used as the source of chemical shift information'
    self.shiftListPulldownA = PulldownList(frame, grid=(0,3), tipText=tipText,
                                           callback=self.changeShiftListA)

    label = Label(frame, text='Residue Selection:', grid=(0,5))
    tipText = 'Sets which chain and range of residues to show assignment status for, e.g. "2-10, 13-77" or "A1-99, B20-120"'
    self.residueEntry = Entry(frame, text='', grid=(0,6), tipText=tipText,
                              returnCallback=self.reportAssignStatus)
    
    text = 'Include water exchangeable atoms'
    tipText = 'Whether to include water exchangeable atoms in the statistics'
    self.waterExchangeableCheckButton = CheckButton(frame, grid=(1,0), text=text, tipText=tipText,
                                                    callback=self.refresh)
    
    tipTexts = ['The name of a category into which atoms of the selected chain/residue have been grouped',
                'The total number of atoms (or equivalent sets for methyl groups) that are available for assignment within the category',
                'The number of atoms (or equivalence groups) that have actually been assigned within the category',
                'The percentage of available atoms (or equivalence groups) hat have been assigned within the category',
                'The lowest chemical shift value of the assigned atoms within in the category',
                'The average chemical shift value of the assigned atoms within in the category',
                'The highest chemical shift value of the assigned atoms within in the category']
                
    headingList = ['Category','Available','Assigned','% Assigned',
                   'Min\nShift','Mean\nShift','Max\nShift']
    self.atomMatrix = ScrolledMatrix(frameA, headingList=headingList, grid=(2,0),
                                     multiSelect=True, callback=self.selectObject,
                                     tipTexts=tipTexts)
    
   
    #
    # NOE
    #
   
    frameB.grid_rowconfigure(1, weight=1)
     
    frame = Frame(frameB, grid=(0,0),sticky='ew')
    frame.grid_columnconfigure(3, weight=1)
   
    label = Label(frame, text='Mol System:', grid=(0,0))
    tipText = 'Selects which molecular system (group of chains) to show through-space connectivity information for'
    self.molSystemPulldownB = PulldownList(frame, callback=self.changeMolSystemB,
                                           grid=(0,1), tipText=tipText)
    
    label = Label(frame, text='Peak List:', grid=(0,2))
    tipText = 'Selects which through-space peak list(s) to show connectivity information for'
    self.peakListPulldownB = PulldownList(frame, callback=self.changePeakListB,
                                          grid=(0,3), tipText=tipText)
     
    tipTexts = ['The identity of the chain or residue that the data in the row pertains',
                'The total number of though-space connections for the chain/residue in the selected peak list(s)',
                'The number of connections that lie within a single residue; do not cross between different residues',
                'The number of connections that cross from one residue to another',
                'The number of connections between residues that are sequential neighbours',
                'The number of connections between residues that are within four sequence positions of one another',
                'The number of connections between residues that are not sequential neighbours but are within four sequence positions of one another',
                'The number of connections between residues that are five or more sequence positions apart',
                'The number of connections that lie within a single molecular chain',
                'The number of connections that cross between different molecular chains',
                'For a given residue, a list of sequence numbers for all other residues that are connected']
    headingList = ['Seq\nElement','Total','Intra\nResidue','Inter\nResidue',
                   'Sequential','Short\nRange','Short\nNon-seq','Long\nRange',
                   'Intra\nChain','Inter\nChain','Contacted Residues']

    justifyList = ['center'] * 10
    justifyList.append('left')

    self.noeMatrix = ScrolledMatrix(frameB, headingList=headingList,
                                    justifyList=justifyList,
                                    multiSelect=True, callback=self.selectObject,
                                    grid=(1,0), tipTexts=tipTexts)
   
    #
    # Resonances
    #
   
    frameC.grid_rowconfigure(1, weight=1)
     
    frame = Frame(frameC)
    frame.grid(row=0, column=0,sticky='ew')
    frame.grid_columnconfigure(2, weight=1)
    
    label = Label(frame, text='Shift List:', grid=(0,0))
    tipText = 'Selects which shift list to derive chemical shift values from in the resonance analysis'
    self.shiftListPulldownC = PulldownList(frame, grid=(0,1), tipText=tipText,
                                           callback=self.changeShiftListC)
    
    label = Label(frame, text='Show only alerts:', grid=(0,3))
    tipText = 'Sets whether to show only rows in the table for resonances that carry warnings'
    self.resAlertCheck = CheckButton(frame, selected=False, grid=(0,4),
                                     callback=self.reportChemShifts,
                                     tipText=tipText)

    tipTexts = ['The serial number of the resonance',
                'The nuclear isotope type of the resonance',
                'The assignment annotation for the resonance; just its serial number in square brackets if unassigned',
                'The assignment annotations of all the other resonances that are considered covalently bound; either by virtue of atom assignments or "onebond" linked experimental dimensions',
                'When atom assigned, the average chemical shift for the atom type given in BioMagResBank data',
                'When atom assigned, the sequence-adjusted random coil chemical shift value for the atom type',
                'When atom assigned, a score indicating how well the resonance chemical shift fits with the known distribution of shifts for the atom type',
                'The chemical shift value of the resonance, in the selected shift list',
                'The standard deviation of the chemical shift value; over the peak dimensions to which the resonance is assigned',
                'The maximum difference in the averaged chemical shift value and a peak dimension assigned to the resonance',
                'The number of peak dimensions the resonance is assigned to']
    headingList = ['#','Iso.','Resonance','Bound','BMRB\nMean',
                   'Random\nCoil','Type\nScore','Value','SD',
                   u'Max Peak\n\u0394','Dim\nContribs']
    self.resonanceMatrix = ScrolledMatrix(frameC, headingList=headingList,
                                          multiSelect=True, callback=self.selectObject,
                                          tipTexts=tipTexts, grid=(1,0))
   
    
    #
    # Peaks
    #
    
    
    frameD.grid_rowconfigure(1, weight=1)
    frameD.grid_rowconfigure(3, weight=1)
     
    frame = Frame(frameD)
    frame.grid(row=0, column=0,sticky='ew')
    frame.grid_columnconfigure(2, weight=1)
   
    label = Label(frame, text='Peak List:', grid=(0,0))
    tipText = 'Selects which peak list to show peak analyses for'
    self.peakListPulldownD = PulldownList(frame, grid=(0,1), tipText=tipText,
                                          callback=self.changePeakListD)

    label = Label(frame, text='Show only alerts:', grid=(0,3))
    tipText = 'Sets whether to show only rows in the table for peaks that carry warnings'
    self.peakAlertCheck = CheckButton(frame, grid=(0,4), selected=False,
                                      callback=self.reportPeaks,
                                      tipText=tipText)


    tipTexts = ['The serial number of the peak in the selected peak list',
                'The resonance assignment annotation of the peak',
                'Any error that relate to the peak\'s resonance assignment',
                'The position of the peak, in the units of the experimental axes; usually ppm',
                'The integration volume of the peak; colour indicates deviation from the peak list average on a logarithmic scale',
                'The standard deviation of the peak integration volume on a logarithmic scale; highlights strong/weak outliers (which could be genuine)',
                'The height intensity of the peak; colour indicates deviation from the peak list average on a logarithmic scale',
                'The standard deviation of the height intensity  on a logarithmic scale; highlights strong/weak outliers (which could be genuine) and inappropriate peak sign',
                'An indication of whether or not there is a symmetry peak in the same location and with the same assignment',
                'The nearest other peak',
                'The distance to the nearest other peak, an indication of overlap']
    headingList = ['#','Assignment','Assign\nErrors',
                   'Location','Volume','Log Vol\nDev',
                   'Height','Log Height\nDev',
                   'Symmetry\nPeaks', 'Nearest\nPeak', 'Nearest\nDistance']
    
    self.peakMatrixTipTexts = tipTexts
    self.peakMatrix = ScrolledMatrix(frameD, headingList=headingList,
                                     multiSelect=True, callback=self.selectObject,
                                     tipTexts=tipTexts, grid=(1,0))
    
    frame = Frame(frameD, grid=(2,0))
    frame.grid_columnconfigure(2, weight=1)
   
    label = Label(frame, text='Mol System:', grid=(0,0))
    tipText = 'Selects which molecular system (group of chains) to display assignment status information for'
    self.molSystemPulldownD = PulldownList(frame, grid=(0,1), tipText=tipText,
                                           callback=self.changeMolSystemD)

    tipTexts = ['The row number',
                'The chain code',] + 7*['']
    headingList = ['#','Chain','Seq Code','Ccp Code',
                   'Atom','Resonances','Shifts',
                   'Delta', 'Stdev']
    
    self.peakAtomMatrix = ScrolledMatrix(frameD, headingList=headingList,
                                         multiSelect=True, callback=self.selectObject,
                                         tipTexts=tipTexts, grid=(3,0))
    
    
    #
    # Main
    #

    sideButtons = UtilityButtonList(tabbedFrame.sideFrame, helpUrl=self.help_url,
                                    grid=(0,0), sticky='e')
                                    
  
    tipTexts = ['Manually update the currently viewed table in light of assignment changes etc.',
                'Show a table of the peaks that relate to the selected rows, either directly or by virtue of assignment (depending on the current table)',
                'Show a table of the peaks and symmetry peaks that relate to the selected rows',
                'Show a table of resonances that relate to the selected rows, either directly or by virtue of assignment (depending on the current table)']
    texts = ['Update Table','Show Peaks','Show Peaks and Symmetry Peaks','Show Resonances']
    commands = [self.refresh, self.showObjects, self.showSymmetryPeaks, self.showResonances]
    self.bottomButtons = ButtonList(guiFrame, texts=texts, tipTexts=tipTexts,
                                    commands=commands, grid=(1,0))
    self.showButton = self.bottomButtons.buttons[1]
    self.showButton2 = self.bottomButtons.buttons[2]
    
    self.updateMolSystems()
    self.updateShiftLists()
    self.updatePeakLists()
    self.reportAssignStatus()
    
    self.administerNotifiers(self.registerNotify)

  def administerNotifiers(self, notifyFunc):

    for func in ('__init__', 'delete','setName'):
      for clazz in ('ccp.nmr.Nmr.DataSource', 'ccp.nmr.Nmr.Experiment', 'ccp.nmr.Nmr.PeakList'):
        notifyFunc(self.updatePeakLists, clazz, func)
      
    for func in ('__init__', 'delete'):
      notifyFunc(self.updateShiftLists, 'ccp.nmr.Nmr.ShiftList', func)

    for func in ('__init__', 'delete'):
      for clazz in ('ccp.molecule.MolSystem.MolSystem', 'ccp.molecule.MolSystem.Chain'):
        notifyFunc(self.updateMolSystems, clazz, func)
  
  def refresh(self, *junk):
  
    self.selectTab(self.tabbedFrame.selected)
        
  def selectTab(self, i):
    
    commands = [self.reportAssignStatus,
                self.reportNoe,
                self.reportChemShifts,
                self.reportPeaksAndAtoms]
                
    commands[i]()
  
  def showResonances(self):
  
    if type(self.object) == type(''):
      return

    resonances = set()
    shiftList = None
    
    if self.object:
      if self.object.className == 'Peak':
        shiftList = self.peakListD.dataSource.experiment.shiftList
        for peak in self.peakMatrix.currentObjects:
          for peakDim in peak.peakDims:
            for contrib in peakDim.peakDimContribs:
              if contrib.peakDimComponent:
                continue
              
              resonances.add(contrib.resonance)
        
      elif self.object.className == 'Atom':
        for atom in self.peakAtomMatrix.currentObjects:
          resonances.update(getAtomResonances(atom))

      elif self.object.className == 'Shift':
        shiftList = self.shiftListC
        for shift in self.resonanceMatrix.currentObjects:
          resonances.add(shift.resonance)
      
      elif self.object.className == 'Residue':
        shiftList = self.peakListB.dataSource.experiment.shiftList
        for residue in self.noeMatrix.currentObjects:
          for atom in residue.atoms:
            atomSet = atom.atomSet
            
            if atomSet:
              for resonanceSet in atomSet.resonanceSets:
                resonances.update(resonanceSet.resonances)
      
      elif self.object.className == 'Chain':
      
        shiftList = self.peakListB.dataSource.experiment.shiftList
        for residue in self.object.residues:
          for atom in residue.atoms:
            atomSet = atom.atomSet
            
            if atomSet:
              for resonanceSet in atomSet.resonanceSets:
                resonances.update(resonanceSet.resonances)
 
      if resonances:
        self.guiParent.viewSelectedResonances(resonances, shiftList)
    
  def getRelatedPeaks(self):

    peaks = None
    object = self.object
    if object and type(object) != type(''):
    
      residues = []
    
      className = object.className
      if className == 'Peak':
        peaks = self.peakMatrix.currentObjects
        
      elif className == 'Shift':
        peaks = set()
        shiftList = self.shiftListC
        for shift in self.resonanceMatrix.currentObjects:
          for contrib in shift.resonance.peakDimContribs:
            peak = contrib.peakDim.peak
            #peaks.add(peak)
            
            experiment = peak.peakList.dataSource.experiment
            if experiment.shiftList is shiftList:
              peaks.add(peak)
            
      
        peaks = list(peaks)
      
      elif className == 'Residue':
      
        residues = self.noeMatrix.currentObjects
      
      elif className == 'Chain':
      
        residues = object.residues

      if residues:
      
        noesyPeakList = {}
        for peakList in getThroughSpacePeakLists(self.project):
          noesyPeakList[peakList] = True

        peaks = set()
        for residue in residues:
          for atom in residue.atoms:
            if atom.atomSet:
              for resonanceSet in atom.atomSet.resonanceSets:
                for resonance in resonanceSet.resonances:
                  for contrib in resonance.peakDimContribs:
                    peak = contrib.peakDim.peak
                    if noesyPeakList.get(peak.peakList):
                       peaks.add(peak)

        peaks = list(peaks)

    return peaks

  def showObjects(self):
  
    peaks = self.getRelatedPeaks()
    if peaks:
      self.guiParent.viewPeaks(peaks)

  def showSymmetryPeaks(self):
  
    peaks = self.getRelatedPeaks()
    if peaks:
      allPeaks = []
      for peak in peaks:
        allPeaks.append(peak)
        symmetryPeaks = findSymmetryPeaks(peak)
        allPeaks.extend(symmetryPeaks)
      self.guiParent.viewPeaks(allPeaks)

  def selectObject(self, object, row, col):
  
    self.object = object

  def getResidueSelection(self):
    
    from re import match
    
    if not self.molSystemA:
      return 

    if not self.molSystemA.chains:
      return 
      
    text = self.residueEntry.get().strip()
    
    resRanges = text.split(',') 
    
    if not resRanges:
      return
        
    residues = []
    for resRange in resRanges:
      resRange = resRange.strip()
      
      matchObj = match('([a-zA-Z_]*)\s*(-*\d+)\s*-\s*([a-zA-Z_]*)\s*(-*\d+)', resRange)
      
      if matchObj:
        chain1 = matchObj.group(1)
        res1 = int(matchObj.group(2))
        chain2 = matchObj.group(3)
        res2 = int(matchObj.group(4))
      
      else:
        matchObj = match('([a-zA-Z_]*)\s*(-*\d+)', resRange)
          
        if matchObj:
          chain1 = matchObj.group(1)
          res1 = int(matchObj.group(2))
          chain2 = None
          res2 = res1
        
        else:
          msg = 'Cannot interpret residue range "%s"' % resRange
          showWarning('Failure',msg, parent=self)
          continue
        
      if chain2 and (chain1 != chain2):
        msg = 'Residue range "%s" has ambiguous chain specification' % resRange
        showWarning('Failure',msg, parent=self)
        continue
      
      if chain1:
        chain = self.molSystemA.findFirstChain(code=chain1)
        if not chain:
          msg = 'Could not find chain "%s" in residue range specification' % chain1
          showWarning('Failure',msg, parent=self)
          continue
      
      else:
        chain = self.molSystemA.sortedChains()[0]
        
      if res1 > res2:
        res1, res2 = res2, res1
      
      residue1 =  chain.findFirstResidue(seqCode=res1) 
      if not residue1:
        msg = "Couldn't find residue %d in range specification" % res1
        showWarning('Failure',msg, parent=self)
        continue
      
      residue2 =  chain.findFirstResidue(seqCode=res2) 
      if not residue2:
        msg = "Couldn't find residue %d in range specification" % res2
        showWarning('Failure',msg, parent=self)
        continue
      
      for residue in chain.sortedResidues():
        if res1 <= residue.seqCode <= res2:
          residues.append(residue)
    
    return residues      

  def reportAssignStatus(self, *event):
  
    self.showButton.disable()
    self.showButton.config(background='#D8D8D8')  
    self.showButton2.disable()
    self.showButton2.config(background='#D8D8D8')  
   
    textMatrix = []
    objectList = []
    colorMatrix = []
    
    if self.molSystemA and self.shiftListA:
      residueSelection = self.getResidueSelection()
      excludeWaterExchangeable = not self.waterExchangeableCheckButton.get()
      data = analyseAssignmentCompleteness(self.molSystemA,
                                           self.shiftListA,
                                           residueSelection,
                                           excludeWaterExchangeable)
            
      for datum in data:
      
        objectList.append(datum[0])
        textMatrix.append(datum)
        colorMatrix.append([None] * 7)

    self.atomMatrix.update(textMatrix=textMatrix,
                           objectList=objectList,
                           colorMatrix=colorMatrix)  

  def reportNoe(self):
  
    self.showButton.enable()
    self.showButton.config(background='#B0B0FF')
    self.showButton2.enable()
    self.showButton2.config(background='#B0B0FF')  
    
    textMatrix = []
    objectList = []
    colorMatrix = []

    if self.molSystemB:
    
      if self.peakListB:
        peakLists = [self.peakListB,]
      else:
        peakLists = None  
    
      data = analyseNoeAssignments(self.molSystemB, peakLists)
    
      for residue, datum in data:
      
        objectList.append(residue)
        textMatrix.append(datum)
        colorMatrix.append([None] * 11)

    self.noeMatrix.update(textMatrix=textMatrix,objectList=objectList,
                          colorMatrix=colorMatrix)  

  def reportChemShifts(self, *event):
  
    self.showButton.enable()
    self.showButton.config(background='#B0B0FF')
    self.showButton2.enable()
    self.showButton2.config(background='#B0B0FF')  
    #self.peakListPulldown.disable()

    #  T1, shifts -  Row per resonance - Cols: shift, shiftSD, maxPeak delta, SD per spectrum ++
    textMatrix  = []
    objectList  = []
    colorMatrix   = []
                   
    onlyAlerts = self.resAlertCheck.get()
    
    if self.shiftListC:
    
      data = analyseChemicalShifts(self.shiftListC)
      
      for shift, datum in data:
        num, isotope, name, boundResonances, mean, \
         coil, prob, value, sd, peakDelta, \
         nContribs, isDuplicate, boundWarn = datum
         
        color = [None] * 11
        alert = False
      
        t1 = ppmSdThresholds.get(isotope, 0.1)
        t2 = ppmDeltaThresholds.get(isotope,0.15)
        
        isotope = shift.resonance.isotope
        
        if isotope:
          color[0] = color[1] = atomColorDict.get(isotope.chemElement.symbol)
        
        if isDuplicate:
          datum[2] += ' repeat'
          color[2] = '#FF0000'
          alert = True
        
        if boundWarn:
          color[3] = '#FF0000'
          alert = True
        
        if prob is not None:
          if prob < 0.001:
            color[7] = '#FF0000'
            color[6] = '#FF0000'
            alert = True
          elif prob < 0.01:
            color[6] = '#d0b8a0'
          elif prob < 0.05:
            color[6] = '#d0d0a0'
            
        if sd > t1:
          alert = True
          color[8] = '#FF0000'
      
        if peakDelta > t2:
          alert = True
          color[9] = '#FF0000'
        
        if onlyAlerts and not alert:
          continue
          
        objectList.append(shift)
        textMatrix.append(datum[:-2])
        colorMatrix.append(color)
 
    self.resonanceMatrix.update(textMatrix=textMatrix,
                                objectList=objectList,
                                colorMatrix=colorMatrix)  

  def commonResonances(self, peakDim1, peakDim2):

    # say it is the same if any resonance is in both (a bit liberal)

    resonances1 = getPeakDimResonances(peakDim1)
    resonances2 = getPeakDimResonances(peakDim2)

    return resonances1 & resonances2

  def boundResonancesOfCorrectIsotopeCode(self, peak, dataDim, isotopeCode):

    peakDim = peak.findFirstPeakDim(dataDim=dataDim)
    resonances = getPeakDimResonances(peakDim)
    resonanceSet = set()
    for resonance in resonances:
      boundResonances = getBoundResonances(resonance)
      resonanceSet.update([boundResonance for boundResonance in boundResonances if isotopeCode == boundResonance.isotopeCode])

    return resonanceSet

  def determineSymmetryPeakColor(self, dataDimPairs, peak, symmetryPeaks):

    if not dataDimPairs:
      return None

    pairs = []

    for dataDim1, dataDim2, isotopeCode1, isotopeCode2 in dataDimPairs:
      if isotopeCode1:
        boundResonances = self.boundResonancesOfCorrectIsotopeCode(peak, dataDim2, isotopeCode1)
      else: # isotopeCode2
        boundResonances = self.boundResonancesOfCorrectIsotopeCode(peak, dataDim1, isotopeCode2)

      if boundResonances: # it is worth checking for symmetry peaks
        peakDim1 = peak.findFirstPeakDim(dataDim=dataDim1)
        peakDim2 = peak.findFirstPeakDim(dataDim=dataDim2)
        pairs.append((peakDim1, peakDim2, dataDim1, dataDim2))


    if not pairs:
      return None

    if not symmetryPeaks:
      return '#FF0000'

    for peak2 in symmetryPeaks:
      for peakDim1, peakDim2, dataDim1, dataDim2 in pairs:
        peakDim3 = peak2.findFirstPeakDim(dataDim=dataDim1)
        peakDim4 = peak2.findFirstPeakDim(dataDim=dataDim2)
        if not self.commonResonances(peakDim1, peakDim4) or not self.commonResonances(peakDim2, peakDim3):
          break
      else:
        #return '#a0FFa0'
        return None

    return '#FFa0a0'

  def getOnebondIsotopeCode(self, dataDim):

    expDim = dataDim.expDim
    expDimRef = getPrimaryExpDimRef(expDim)
    expDimRefPairs = getOnebondExpDimRefs(expDim.experiment)

    for expDimRef1, expDimRef2 in expDimRefPairs:
      if expDimRef1 is expDimRef:
        return expDimRef2.isotopeCodes[0]
      elif expDimRef2 is expDimRef:
        return expDimRef1.isotopeCodes[0]

    return None
    
  def reportPeaksAndAtoms(self, *event):

    self.reportPeaks()
    self.reportAtoms()
  
  def reportPeaks(self, *event):
  
    self.showButton.enable()
    self.showButton.config(background='#B0B0FF')
    self.showButton2.enable()
    self.showButton2.config(background='#B0B0FF')  
    #self.peakListPulldown.enable()
    
    #   Row per peak - Cols: height vs PL mean, vol vs PL mean,  sign check, F1 shift delta, F2 shift delta ++
    textMatrix  = []
    objectList  = []
    colorMatrix = []

    headingList = ['#','Assignment','Assign\nErrors',
                   'Location','Volume','Log Vol\nDev',
                   'Height','Log Height\nDev',
                   'Symmetry\nPeaks', 'Nearest\nPeak', 'Nearest\nDistance']
    
    tipTexts = self.peakMatrixTipTexts[:]
                   
    onlyAlerts = self.peakAlertCheck.get()

    if self.peakListD:
      for i in range(self.peakListD.dataSource.numDim):
        headingList.append(u'F%d Max\nShift \u0394' % (i+1))
        tip = "For dimension %d, the maximum difference between an assigned resonance's (averaged) chemical shift value and the peak position"
        tipTexts.append(tip % (i+1))

      equivDataDimPairs = getEquivalentDataDims(self.peakListD.dataSource)
      dataDimPairs = []
      for dataDim1, dataDim2 in equivDataDimPairs:
        isotopeCode1 = self.getOnebondIsotopeCode(dataDim1)
        isotopeCode2 = self.getOnebondIsotopeCode(dataDim2)
        # no point checking if both fixed or if neither fixed
        if (isotopeCode1 and not isotopeCode2) or (isotopeCode2 and not isotopeCode1):
          dataDimPairs.append((dataDim1, dataDim2, isotopeCode1, isotopeCode2))
    
      posProp, volSd, htSd, data = analysePeaks(self.peakListD)

      for peak, datum in data:
        colors = []
        serial,annotation,assignErrors,location,volume,volDelta,height,htDelta,maxDeltas, nearestPeak, nearestDistance, symmetryPeaks = datum
        alert = False

        colors = [None] * 11
        maxDeltaData = []
        for maxDelta, isotope in maxDeltas:
          maxDeltaData.append(maxDelta)
          
          th = ppmDeltaThresholds.get(isotope,0.15)
          if maxDelta and (maxDelta > th):
            colors.append('#FF0000')
          else:
            colors.append(None)
        
        if volSd:
          if volDelta is None:
            nVolSd = 999
          else:
            nVolSd = volDelta/volSd
        else:
          nVolSd = 0
          
        if htSd:
          if htDelta is None:
            nHtSd = 999  
          else:
            nHtSd = htDelta/htSd
        else:
          nHtSd = 0  
        
        if assignErrors:
          alert = True
          colors[1] = '#FF0000'
          colors[2] = '#FF0000'
  
        if nVolSd > 3:
          alert = True
          colors[4] = '#FF0000'
          colors[5] = '#FF0000'
        elif nVolSd > 2:
          colors[4] = '#FFa0a0'
          colors[5] = '#FFa0a0'
        
        if nHtSd > 3:
          alert = True
          colors[6] = '#FF0000'
          colors[7] = '#FF0000'
        
        elif nHtSd > 2:
          colors[6] = '#FFa0a0'
          colors[7] = '#FFa0a0'
  
        if posProp > 0.75:
          if height < 0.0:
            alert = True
            colors[6] = '#FF0000'
            if volume < 0.0:
              colors[4] = '#FF0000'
       
        elif posProp < 0.25:
          if height > 0.0:
            alert = True
            colors[6] = '#FF0000'
            if volume > 0.0:
              colors[4] = '#FF0000'
         
        colors[8] = self.determineSymmetryPeakColor(dataDimPairs, peak, symmetryPeaks)
        if colors[8]:
          alert = True

        colors[9] = colors[10] = None
        if nearestPeak:
          if nearestDistance < 0.5:
            colors[10] = '#FF0000'
            alert = True
          elif nearestDistance < 1.0:
            colors[10] = '#FFa0a0'
            alert = True

        if onlyAlerts and not alert:
          continue

        datum = [serial,annotation,
                 ','.join(assignErrors),location,
                 volume, nVolSd, height, nHtSd]

        datum.append(', '.join([str(pk.serial) for pk in symmetryPeaks]))
        datum.extend([nearestPeak and nearestPeak.serial, nearestDistance])

        datum.extend(maxDeltaData)

        objectList.append(peak)
        textMatrix.append(datum)
        colorMatrix.append(colors)

    self.peakMatrix.update(textMatrix=textMatrix,
                           objectList=objectList,
                           headingList=headingList,
                           colorMatrix=colorMatrix,
                           tipTexts=tipTexts)  

  def reportAtoms(self, *event):

    molSystem = self.molSystemD

    if molSystem:
      chains = molSystem.sortedChains()
      shiftList = self.peakListD.dataSource.experiment.shiftList
    else:
      chains = []

    textMatrix = []
    objectList = []
    n = 0
    for chain in chains:
      for residue in chain.sortedResidues():
        for atom in residue.sortedAtoms():
          if atom.name[0] not in ('H', 'C', 'N'):
            continue
          chemAtom = atom.chemAtom
          resonances = sorted(list(getAtomResonances(atom)))
          resonanceText = ', '.join(['%d' % resonance.serial for resonance in resonances])
          if shiftList:
            values = []
            deltas = []
            randomCoilShift = getRandomCoilShift(chemAtom)
            for resonance in resonances:
              shift = resonance.findFirstShift(parentList=shiftList)
              if shift:
                values.append('%.3f' % shift.value)
                deltas.append('%.3f' % (shift.value - randomCoilShift))
              else:
                values.append(None)
                deltas.append(None)
            shiftText = ', '.join(['%s' % value for value in values])
            deltaText = ', '.join(['%s' % delta for delta in deltas])
          else:
            shiftText = ''
            deltaText = ''

          datum = [n+1, chain.code, residue.seqCode, residue.ccpCode,
                   atom.name, resonanceText, shiftText, deltaText, 0]
          textMatrix.append(datum)
          objectList.append(atom)
  
    self.peakAtomMatrix.update(textMatrix=textMatrix,
                               objectList=objectList)

  def changeMolSystemA(self, molSystem):
  
    if molSystem is not self.molSystemA:
      self.molSystemA = molSystem
      resRanges = []
      for chain in molSystem.sortedChains():
        residues = chain.sortedResidues()
        if not residues:
          continue
          
        if len(residues) > 1:
          data = (chain.code,residues[0].seqCode,residues[-1].seqCode)
          resRanges.append( '%s%d-%d' %  data)
        else:
          data = (chain.code,residues[0].seqCode)
          resRanges.append( '%s%d' % (data) )
          
      self.residueEntry.set(','.join(resRanges))
      self.refresh()

  def changeMolSystemB(self, molSystem):
  
    if molSystem is not self.molSystemB:
      self.molSystemB = molSystem
      self.refresh()
  
  def changeMolSystemD(self, molSystem):
  
    if molSystem is not self.molSystemD:
      self.molSystemD = molSystem
      self.refresh()
  
  def changeShiftListA(self, shiftList):
  
    if shiftList is not self.shiftListA:
      self.shiftListA = shiftList
      self.refresh()
  
  def changeShiftListC(self, shiftList):
  
    if shiftList is not self.shiftListC:
      self.shiftListC = shiftList
      self.refresh()
  
  def changePeakListB(self, peakList):
  
    if peakList is not self.peakListB:
      self.peakListB = peakList
      self.refresh()
      
  def changePeakListD(self, peakList):
  
    if peakList is not self.peakListD:
      self.peakListD = peakList
      self.refresh()
  
  def getPeakLists(self):
  
    peakLists = []
    
    for experiment in self.nmrProject.sortedExperiments():
      for spectrum in experiment.sortedDataSources():
        if spectrum.dataType == 'processed':
          for peakList in spectrum.sortedPeakLists():
            peakLists.append(peakList)
    
    return peakLists      
  
  
  def getShiftLists(self):
  
    return getShiftLists(self.nmrProject)

   
  def getMolSystems(self):
    
    return self.project.sortedMolSystems()
 
 
  def updateMolSystems(self, *obj):
  
    molSystems = self.getMolSystems()
    indexA = 0
    indexB = 0
    names  = [ms.code for ms in molSystems]
    molSystemA = self.molSystemA
    molSystemB = self.molSystemB
    
    if names:
      if molSystemA not in molSystems:
        molSystemA = molSystems[0]
        
      if molSystemB not in molSystems:
        molSystemB = molSystems[0]
      
      indexA = molSystems.index(molSystemA)
      indexB = molSystems.index(molSystemB)
    
    if molSystemA is not self.molSystemA:
      self.changeMolSystemA(molSystemA)

    if molSystemB is not self.molSystemB:
      self.changeMolSystemB(molSystemB)
     
    self.molSystemPulldownA.setup(names, molSystems, indexA)
    self.molSystemPulldownB.setup(names, molSystems, indexB)

    self.updateMolSystemD()

  def updateMolSystemD(self, *obj):
  
    indexD = 0
    molSystemD = self.molSystemD

    if self.peakListD:

      experiment = self.peakListD.dataSource.experiment
      molSystems = experiment.sortedMolSystems()

      names  = [ms.code for ms in molSystems]

    else:

      molSystems = []
      names = []

    if names:
      if molSystemD not in molSystems:
        molSystemD = molSystems[0]
      
      indexD = molSystems.index(molSystemD)
    
    if molSystemD is not self.molSystemD:
      self.changeMolSystemD(molSystemD)

    self.molSystemPulldownD.setup(names, molSystems, indexD)
    
  def updateShiftLists(self, *obj):
  
    shiftLists = self.getShiftLists()
    indexA = 0
    indexC = 0
    names  = ['Shift List %d' % sl.serial for sl in shiftLists]
    if names:
      if self.shiftListA not in shiftLists:
        self.shiftListA = shiftLists[0]
        
      if self.shiftListC not in shiftLists:
        self.shiftListC = shiftLists[0]
      
      indexA = shiftLists.index(self.shiftListA)  
      indexC = shiftLists.index(self.shiftListC)  
    
    self.shiftListPulldownA.setup(names, shiftLists, indexA)
    self.shiftListPulldownC.setup(names, shiftLists, indexC)
    
  def updatePeakLists(self, *obj):
  
    indexB = 0
    peakLists = getThroughSpacePeakLists(self.project)
    names  = ['%s:%s:%d' % (pl.dataSource.experiment.name, \
              pl.dataSource.name, pl.serial) for pl in peakLists]

    peakLists.insert(0, None)
    names.insert(0, '<All>')

    if names:
      if self.peakListB not in peakLists:
        self.peakListB = peakLists[0]
        
      indexB = peakLists.index(self.peakListB)  
    
    self.peakListPulldownB.setup(names, peakLists, indexB)


    indexD = 0
    peakLists = self.getPeakLists()
    names  = ['%s:%s:%d' % (pl.dataSource.experiment.name, \
              pl.dataSource.name, pl.serial) for pl in peakLists]
    
    if names:
      if self.peakListD not in peakLists:
        self.peakListD = peakLists[0]
        
      indexD = peakLists.index(self.peakListD)  
    
    self.peakListPulldownD.setup(names, peakLists, indexD)
    
    self.updateMolSystemD()
    
  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

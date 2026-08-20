
"""
======================COPYRIGHT/LICENSE START==========================

CopyAssignments.py: Part of the CcpNmr Analysis program

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

from math import sqrt

from memops.gui.MessageReporter import showOkCancel, showWarning

from ccpnmr.analysis.popups.BasePopup     import BasePopup
from ccpnmr.analysis.core.AssignmentBasic import propagatePeakAssignments, assignResToDim, getResidueResonances
from ccpnmr.analysis.core.MoleculeBasic   import getChainResidueMapping, duplicateResidueAssignments, getResidueCode
from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes
from ccpnmr.analysis.core.MarkBasic       import createPeakMark
from ccpnmr.analysis.core.PeakBasic       import arePeaksAssignedSame, getIsotopeWeightedTolerances, findShiftDistPeakMatches
from ccpnmr.analysis.core.WindowBasic     import isSpectrumInWindowPane, getDataDimAxisMapping
from ccpnmr.analysis.core.WindowBasic     import getWindowPaneName, getActiveWindows

from memops.gui.ButtonList      import ButtonList, UtilityButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Frame           import Frame
from memops.gui.LabelFrame      import LabelFrame
from memops.gui.LabelDivider    import LabelDivider
from memops.gui.Label           import Label
from memops.gui.MessageReporter import showOkCancel
from memops.gui.ProgressBar     import ProgressBar
from memops.gui.PulldownList    import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix 
from memops.gui.TabbedFrame     import TabbedFrame

# "From" spin systems
#
# Select spin systems table
# Select molSystem pulldown
# Select chain pulldown
# Pick peaks [] Pick non-extrema peaks (at ref position) []
# Pick tolerances
#
# Follow in window pulldown 
#
#  - Grid (always?) C,H,N - Z axis?
#  - Z is bound to X or Y (change strip dir)
#  - X,Z of window must match bound of expt
#  - Irrelevant for 2D  
#  - Consider existing code

COPY_OPERATIONS = ["Don't Move",'Move to Target','Duplicate']

class CopyAssignmentsPopup(BasePopup):

  """
  **Copy Resonance Assignments to Different Peak Lists or Molecules**
  
  This popup has two separate but related purposes: firstly, to transfer
  resonance assignments from one peak list to another peak list where the peaks
  are in similar, but not identical, positions; secondly, to move or duplicate
  resonances that are assigned to atoms in one molecular chain to another of
  similar sequence.
  
  **Copying Assignments Between Peak Lists**
  
  The general idea of the first tab is to select two peak lists, one of which 
  acts as the "source" where resonance assignments are copied from and another
  to act as the "target" where assignments are copied to. In general both the
  source and targets peak lists will be from experiments with the same number
  and kinds of axis. However, it is also possible to copy assignments *to* a
  target peak list that is of *lower* dimensionality than the source. For
  example you can copy 3D 15N HSQC-NOESY (H,N,H dimensions) assignments to
  15N-HSQC (H,N dimensions) or 2D NOESY (H,H dimensions), and in these cases
  only resonance assignments that go on dimensions common to source and target
  will be transferred.

  The "Source Peaks" table lists all the peaks that assignments may be copied
  from, together with an indication the number of peaks that are close and the
  one that is closest. The "Target Peaks" table is updated for an individual
  source peak when you click on a row of the "Source Peaks" table; to show you
  details of all the peaks that match in terms of position.
  
  The matching of peak positions is made according to a chemical shift distance
  measure and only possibilities that lie within the distance threshold are 
  considered. The chemical shift distance is calculated by taking the difference
  in peak positions for each dimension, dividing each difference by the scale
  factor for the isotope that appears on that dimension, squaring the
  differences and taking the square root of their summation. The isotope scale
  factors, which you can set in the interface mean that differences in
  dissimilar dimensions can be compared. By default the distances can be thought
  of in terms as the "1H" equivalent.

  If suitable matches are found, resonance assignments are copied between peaks
  either on an individual basis, by selecting the required target in the lower
  table and using [Assign Selected Target], or *en masse* by using the other
  buttons which process all of the peaks in the list; copying assignments if there
  is only a single matching target within the distance threshold, or
  to the  closest matching target. When using the process-all functions it is
  common to start with a very strict/short threshold radius, assign some peaks
  and then increase the threshold to consider the poorer matching ones.

  **Copying Assignments Between Molecular Chains**

  The second tab is used to move resonances' atomic assignments to a different
  molecule/chain, while maintaining assignments to peaks. When copying
  assignments between molecular chains, i.e. the entities  with residues and
  atoms that you assign, the general principle is that you choose one chain as a
  source, to get assignments from, and another as a target, to transfer
  assignments to. Chains may have identical sequences, for example when copying
  assignments within a homodimer, or reasonably different sequences. In the
  latter case a pairwise sequence alignment is used to determine the initial
  mapping between source and target residues.

  With the chain selection setup, the next task is usually to consider how
  moving assignments will affect the peaks that the source chain resonances are
  assigned to. In this regard there are three options and they are all specified
  on a per-experiment basis (and thus affect that experiment's spectra and
  peaks). The first option is to leave an experiment alone ("don't move") so
  that its peaks are not affected at all. The second option is to move to the
  target, whereby all resonance assignments on the peaks are moved to the other
  chain. The third option is to duplicate assignments on the peaks so that they
  are assigned to both the source and the target chains; this makes each peak
  assigned to double the number of resonances. If a resonance is moved entirely
  to the target chain (no peaks left at source) then the resonance's atomic
  assignment is simply pointed to the different chain. If a resonance remains
  partly assigned to the original chain, then a new resonance is made and this
  is the one assigned to the target chain; the old resonance remains on the
  source chain. 

  Where residues don't match exactly between the two chains, any resonances
  assigned to atoms with no direct equivalent (nothing of the same name) will
  still be copied across; their spin system will become the target residue, they
  retain their original atom types, but they will not have full atomic
  assignments. If the destination residue does not appear correctly in the lower
  table, the user may change the residue-residue mapping by double clicking in
  the "Destination" column.
  """

  def __init__(self, parent, *args, **kw):

    self.guiParent  = parent
    self.follow     = True
    self.restrict   = True
    self.waitingPeak = False
    self.overwrite  = True
    self.showCopied = False

    self.considerAliased = True
    self.sourcePeak      = None
    self.targetPeak      = None
    self.sourcePeakList  = None
    self.targetPeakList  = None
    self.windowPane = None

    self.sourceChain = None
    self.targetChain = None
    self.waitingChain = False
    self.guiParent   = parent
    self.origMapping = {}
    self.mapping     = {}
    self.residue     = None
    self.residueDict = {}
    self.expOperations = {}

    BasePopup.__init__(self, parent=parent, title='Assignment : Copy Assignments', **kw)

  def open(self):
  
    BasePopup.open(self)
    self.updatePeakSource()    
    self.updateSourceChains()
    self.updateExperimentChains()
 
  def body(self, guiFrame):

    self.geometry('600x700')

    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)
    
    options = ['Between Peak Lists','Between Molecule Chains']
    tabbedFrame = TabbedFrame(guiFrame, options=options)
    tabbedFrame.grid(row=0, column=0, sticky='nsew')
    frameA, frameB = tabbedFrame.frames

    #
    #
    # Peak Lists
    #
    #
    
    frameA.grid_columnconfigure(0, weight=1)
    frameA.grid_rowconfigure(5, weight=1)
    frameA.grid_rowconfigure(7, weight=1)
    
    row = 0
    div = LabelDivider(frameA, text='Options')
    div.grid(row=row, column=0, sticky='ew')
    
    row += 1
    frame = Frame(frameA)
    frame.grid(row=row, column=0, sticky='ew')
    frame.grid_columnconfigure(3, weight=1)
    
    label = Label(frame, text='Source Peak List:', grid=(0,0))
    tipText = 'Selects the peak list to copy assignments from'
    self.sourcePulldown = PulldownList(frame, tipText=tipText, grid=(0,1),
                                       callback=self.changeSourcePeakList)

    label = Label(frame, text='Target Peak List:', grid=(0,2))
    tipText = 'Selects the peak list which will be assigned'
    self.targetPulldown = PulldownList(frame, grid=(0,3), tipText=tipText, 
                                       callback=self.changeTargetPeakList)

    label = Label(frame, text='Show Already Copied?', grid=(1,0))
    tipText = 'Whether to show source peaks in the table whose assignments are already copied'
    self.followSelect = CheckButton(frame, callback=self.setShowCopied,
                                    grid=(1,1), selected=False, tipText=tipText)
    
    label = Label(frame, text='Overwrite Assignments?', grid=(1,2))
    tipText = 'Whether to overwrite any existing assignments in the target peak list'
    self.overwriteSelect = CheckButton(frame, callback=self.setOverwrite,
                                       grid=(1,3), selected=True, tipText=tipText)

    label = Label(frame, text='Follow Peaks?', grid=(2,0))
    tipText = 'Whether to follow the location of peaks in the spectrum window when clicking on a peak row'
    self.followSelect = CheckButton(frame, callback=self.setFollow,
                                    grid=(2,1), selected=True, tipText=tipText)
    
    label = Label(frame, text='Only Good Matches?', grid=(2,2))
    tipText = 'When selected only source peaks with at least one good target match will be shown, otherwise all source peaks are shown'
    self.restrictSelect = CheckButton(frame, callback=self.setRestrict,
                                      grid=(2,3), selected=True, tipText=tipText)

    label = Label(frame, text='Consider Aliased?', grid=(3,0))
    tipText = 'Whether target peaks could be aliased; position matching can add whole numbers of sweep widths'
    self.aliasedSelect = CheckButton(frame, callback=self.setConsiderAliased,
                                     grid=(3,1), selected=True, tipText=tipText)

    label = Label(frame, text='Follow Window:', grid=(4,0))
    tipText = 'Selects which spectrum window will be used to navigate to source & target peak positions'
    self.windowPulldown = PulldownList(frame, grid=(4,1), tipText=tipText, 
                                       callback=self.setWindowPane)
                                       
    #self.windowPulldown.inactivate()

    label = Label(frame, text='Dist. Threshold:', grid=(4,2))
    tipText = 'The ppm search radius to match peak positions. Note spectrum dimensions will be weighted by the isotope scale factors'
    self.thresholdEntry = FloatEntry(frame, text=0.08, width=7, tipText=tipText, 
                                     returnCallback=self.changeThreshold, grid=(4,3))

    row +=1
    div = LabelDivider(frameA, text='Scale Factors')
    div.grid(row=row, column=0, sticky='ew')
    
    row += 1
    frame = Frame(frameA)
    frame.grid(row=row, column=0, sticky='ew')
    frame.grid_columnconfigure(5, weight=1)

    label = Label(frame, text='1H:', grid=(0,0))
    tipText = 'The scaling factor used to weight ppm distances in 1H dimensions; used in position radius search'
    self.scaleEntry1H = FloatEntry(frame, text=1.0, width=7, tipText=tipText,
                                   returnCallback=self.changeScaleEntry, grid=(0,1))
 
    label = Label(frame, text='15N:', grid=(0,2))
    tipText = 'The scaling factor used to weight ppm distances in 15N dimensions; used in position radius search'
    self.scaleEntry15N = FloatEntry(frame, text=5.0, width=7, tipText=tipText,
                                    returnCallback=self.changeScaleEntry, grid=(0,3))

    label = Label(frame, text='13C:', grid=(0,4))
    tipText = 'The scaling factor used to weight ppm distances in 13C dimensions; used in position radius search'
    self.scaleEntry13C = FloatEntry(frame, text=10.0, width=7, tipText=tipText,
                                    returnCallback=self.changeScaleEntry, grid=(0,5))

    row +=1
    div = LabelDivider(frameA, text='Source Peaks', grid=(row,0))
    
    row += 1
    tipTexts = ['Serial number of the source peak',
                'Assignment of source peak; which will be copied',
                'The number of target peaks that match the source peak within the ppm distance',
                'The smallest, isotope weighted, ppm distance to a target peak',
                'Assignment annotation of closest matching target peak']
    headingList = ['#','Assignment','Num.\nMatches','Closest\nDistance.','Best\nMatch']
    self.sourcePeakMatrix = ScrolledMatrix(frameA, headingList=headingList, tipTexts=tipTexts,
                                           callback=self.selectSourcePeak, grid=(row,0))
    
    row +=1
    div = LabelDivider(frameA, text='Target Peaks', grid=(row,0))

    row +=1
    tipTexts = ['Serial number of target peak possibility',
                'Assignment of the target peak, which may be overwritten',
                'Isotope weighted (by scale factor per dimension) ppm distance between source and target']
    headingList = ['#','Assignment','Distance']
    self.targetPeakMatrix = ScrolledMatrix(frameA, headingList=headingList, tipTexts=tipTexts,
                                           callback=self.selectTargetPeak, grid=(row,0))
    
    row +=1
    tipTexts = ['Copy assignments from the selected source peak to the selected target',
                'Process all source peaks and copy assignments to those that match a single target within the search radius',
                'Process all source peaks and copy assignments to the closest matching target peak (if there is one in the search radius)']
    texts = ['Assign Selected\nTarget','Assign All\nSingly Matched','Assign All\nTo Closest']
    commands = [self.assignSelected, self.assignSingles, self.assignClosest]
    self.peakButtons = ButtonList(frameA, commands=commands, texts=texts,
                                  grid=(row,0), tipTexts=tipTexts)
 
    #
    #
    # Chains
    #
    #
    
    frameB.grid_columnconfigure(0, weight=1)
    frameB.grid_rowconfigure(3, weight=1)
    frameB.grid_rowconfigure(5, weight=1)
   
    row = 0
    div = LabelDivider(frameB, text='Options', grid=(row,0))
    
    row += 1
    frame = Frame(frameB)
    frame.grid(row=row, column=0, sticky='ew')
    frame.grid_columnconfigure(5, weight=1)
    frame.grid_rowconfigure(0, weight=1)

    label = Label(frame, text='Source Chain:', grid=(0,0), sticky='e')
    tipText = 'Selects the molecular chain from which atom assignments will be copied'
    self.sourceChainPulldown = PulldownList(frame, grid=(0,1), tipText=tipText,
                                            callback=self.changeSourceChain)

    label = Label(frame, text='Target Chain:', grid=(0,2))
    tipText = 'Selects the molecular chain which will receive new atom assignments'
    self.targetChainPulldown = PulldownList(frame, grid=(0,3), tipText=tipText,
                                            callback=self.changeTargetChain)

    row += 1
    div = LabelDivider(frameB, text='Peak Assignment Transfers', grid=(row,0))

    row += 1
    self.exptChainPulldown = PulldownList(self, callback=self.setExptChain,
                                          objects=(0,1,2), texts=COPY_OPERATIONS)
    tipTexts = ['Serial number of experiment to consider',
                'Name of experiment to consider',
                'Sets whether atom assignments for an experiment will be transferred entirely, duplicated on both chains or left unaltered']
    headingList = ['#','Experiment','Operation']
    editWidgets      = [None, None, self.exptChainPulldown]
    editGetCallbacks = [None, None, self.getOperation]
    editSetCallbacks = [None, None, self.setExptChain]
    self.exptChainMatrix = ScrolledMatrix(frameB, headingList=headingList,
                                          multiSelect=True, highlightType=1,
                                          editSetCallbacks=editSetCallbacks,
                                          editGetCallbacks=editGetCallbacks,
                                          editWidgets=editWidgets,
                                          tipTexts=tipTexts,
                                          callback=self.selectExperimentChain,
                                          grid=(row,0), gridSpan=(1,2))
    
    row += 1
    frame = Frame(frameB, grid=(row,0), sticky='ew')
    frame.expandGrid(0,1)
    
    label = Label(frame, text='Set operation\nfor selected:', grid=(0,0), sticky='e')
    tipTexts = ['Sets the selected experiments so that their assignments will be left unaltered',
                'Sets the selected experiments so that their assignments will be transferred entirely to the target chain',
                'Sets the selected experiments so that their assignments will be duplicated for both source and target chains']
    commands = [self.setSelectedExperimentsSource,
                self.setSelectedExperimentsTarget,
                self.setSelectedExperimentsDuplicate]
    self.exptChainButtons = ButtonList(frame, texts=COPY_OPERATIONS, tipTexts=tipTexts,
                                       commands=commands, grid=(0,1), expands=True)
    
    row +=1
    div = LabelDivider(frameB, text='Residue Mapping', grid=(row,0))

    self.targetResiduePulldown = PulldownList(self, callback=self.setTargetResidue)
    tipTexts = ['The residue from which assignments may be copied',
                'The number of resonance assignments currently on the source residue\'s atoms',
                'The residue in the target chain to which assignments may be copied',
                'The number of resonance assignments currently on the target residue\'s atoms']
    headingList = ['Source','Source\nResonances','Destination','Destination\nResonances']
    editWidgets      = [None, None,self.targetResiduePulldown]
    editGetCallbacks = [None, None,self.getTargetResidue]
    editSetCallbacks = [None, None,self.setTargetResidue]
    self.residueMatrix = ScrolledMatrix(frameB, headingList=headingList,
                                        multiSelect=True,
                                        tipTexts=tipTexts,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks, 
                                        editWidgets=editWidgets,
                                        callback=self.selectResidue,
                                        grid=(row,0))
                                        
    
    row +=1
    tipTexts = ['Go ahead and copy assignments between chains using the selected options',
                'Enable copying of assignments for the selected residues',
                'Disable copying of assignments for the selected residues']
    texts    = ['Copy Assignments!', 'Enable Selected','Disable Selected']
    commands = [self.copyChainAssignments, self.enableSelected,self.disableSelected]
    self.residueButtons = ButtonList(frameB, texts=texts, commands=commands, 
                                     grid=(row,0), tipTexts=tipTexts)
    self.residueButtons.buttons[0].config(bg='#A0FFA0')

    self.chainButtons = UtilityButtonList(tabbedFrame.sideFrame, helpUrl=self.help_url,
                                          grid=(0,0), sticky='e')

        
    self.updatePeakSource()    
    self.updateSourceChains()
    self.updateExperimentChains()
    
    self.administerNotifiers(self.registerNotify)

  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete'):
      for clazz in ('ccp.molecule.MolSystem.Chain',):
        notifyFunc(self.updateChains, clazz, func)

    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.Experiment',):
        notifyFunc(self.updateExperimentChains, clazz, func)

    for func in ('__init__', 'delete', 'setName'):
      for clazz in ('ccp.nmr.Nmr.PeakList', 'ccp.nmr.Nmr.DataSource',
                    'ccp.nmr.Nmr.Experiment'):
        notifyFunc(self.updatePeakSource, clazz, func)
        notifyFunc(self.updateTarget, clazz, func)
 
    for func in ('__init__', 'delete',):
      notifyFunc(self.updatePeaksAfter, 'ccp.nmr.Nmr.Peak', func)

    for func in ('setAnnotation',):
      notifyFunc(self.updatePeakDim, 'ccp.nmr.Nmr.PeakDim', func)
    
  def destroy(self):
 
    self.administerNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)
  
  
  def updateResidueButtons(self):
  
    self.residueButtons.buttons[0].disable()
    self.residueButtons.buttons[1].disable()
    self.residueButtons.buttons[2].disable()
    self.chainButtons.buttons[0].disable()
  
    if self.sourceChain and self.targetChain:
      self.residueButtons.buttons[0].enable()
      if self.mapping:
        self.chainButtons.buttons[0].enable()
      
      if self.residueMatrix.currentObjects:
       self.residueButtons.buttons[1].enable()
       self.residueButtons.buttons[2].enable()
  
  # Chain functions
  
  def getChainName(self, chain):
  
    return '%s:%s' % (chain.molSystem.code,chain.code)

  def getExperimentChainKey(self, experiment, chain):

    return '%d:%s:%s' % (experiment.serial,chain.molSystem.code,chain.code)

  def updateExperimentChains(self, dummy=None):

    textMatrix  = []
    objectList  = []
    colorMatrix = []

    if self.sourceChain and self.targetChain:
      for experiment in self.nmrProject.sortedExperiments():
        for molSystem in experiment.molSystems:
          for chain in molSystem.sortedChains():
            if chain is self.sourceChain:
              key = self.getExperimentChainKey(experiment, chain)
              operation = self.expOperations.get(key, 1)
              self.expOperations[key] = operation
              
              datum = [experiment.serial,
                       experiment.name,
                       COPY_OPERATIONS[operation]]
 
              colors = [None,None,None]
              if operation == 2:
                colors[2] = '#F08080'
              elif operation == 1:
                colors[2] = '#A0D0A0'
 
              textMatrix.append(datum)
              objectList.append([experiment,chain,key])
              colorMatrix.append(colors)
 

    self.exptChainMatrix.update(textMatrix=textMatrix,
                                colorMatrix=colorMatrix,
                                objectList=objectList)
    self.updateResidueButtons()

  def setSelectedExperimentsSource(self):
  
    if self.sourceChain:
      for experiment, chain, key in self.exptChainMatrix.currentObjects:
        self.expOperations[key] = 0
 
    self.updateExperimentChains()

  def setSelectedExperimentsTarget(self):

    if self.targetChain:
      for experiment, chain, key in self.exptChainMatrix.currentObjects:
        self.expOperations[key] = 1

    self.updateExperimentChains()

  def setSelectedExperimentsDuplicate(self):

    if self.sourceChain and self.targetChain:
      for experiment, chain, key in self.exptChainMatrix.currentObjects:
        self.expOperations[key] = 2

    self.updateExperimentChains()
   
  def getOperation(self, object):
  
    if self.sourceChain and self.targetChain:
      experiment, chain, key = object
      operation = self.expOperations.get(key, 1)
      self.exptChainPulldown.setIndex(operation)
    
  def setExptChain(self, null):
  
    operation = self.exptChainPulldown.getObject()
    
    if self.experimentChain and self.sourceChain and self.targetChain:
      key = self.experimentChain[-1]
      self.expOperations[key] = operation
 
    self.updateExperimentChains()

  def selectExperimentChain(self, obj, row, col):

    self.experimentChain = obj

  def updateChains(self, chain):
    
    self.updateSourceChains()
    self.updateTargetChains()
    self.updateExperimentChains()

  def updateSourceChains(self):
  
    names = []
    index = -1
  
    chains = self.getSourceChains()
    chain  = self.sourceChain
    
    if chains:
      names = ['%s:%s' % (ch.molSystem.code,ch.code) for ch in chains ]
      if chain not in chains:
        chain = chains[0]
        index = 0
        
      else:
        index = chains.index(chain)  
      
    else:
      chain = None
      self.residue = None
      
    if chain is not self.sourceChain:
      self.sourceChain = chain
      self.targetChain = None
      self.updateTargetChains()
      self.getResidueMapping()
  
    self.sourceChainPulldown.setup(names,chains,index)
  
  def updateTargetChains(self):


    names = []
    index = -1
  
    chains = self.getTargetChains()
    chain  = self.targetChain
    
    if chains:
      names = ['%s:%s' % (ch.molSystem.code,ch.code) for ch in chains ]
      if chain not in chains:
        chain = chains[0]
        index = 0
        
      else:
        index = chains.index(chain)  
      
    else:
      chain = None
      
    if chain is not self.targetChain:
      self.targetChain = chain
      self.residueDict = {}
      if chain:
        for residue in chain.residues:
          key = self.getResidueKey(residue)
          self.residueDict[key] = residue
        
      self.getResidueMapping()
  
    self.targetChainPulldown.setup(names,chains,index)
   
  def changeSourceChain(self, chain):
  
    if chain is not self.sourceChain:
      self.sourceChain = chain
      self.updateTargetChains()
      self.getResidueMapping()
      self.updateExperimentChains()
  
  def changeTargetChain(self, chain):
  
    if chain is not self.targetChain:
      self.targetChain = chain
      self.residueDict = {}
      for residue in chain.residues:
        key = self.getResidueKey(residue)
        self.residueDict[key] = residue
      self.getResidueMapping()
      self.updateExperimentChains()
       
  
  def getResidueKey(self, residue):
    
    return '%d%s' % (residue.seqCode,getResidueCode(residue))
  
  def getSourceChains(self):
  
    chains = []
    for molSystem in self.project.sortedMolSystems():
      for chain in molSystem.sortedChains():
        if len(chain.residues) > 1:
          chains.append(chain)

    return chains    
  
  def getTargetChains(self):
  
    chains = []
    if self.sourceChain:
      molTypes = {}
      
      for molResidue in self.sourceChain.molecule.molResidues:
        molType = molResidue.molType
        molTypes[molType] = True
   
      for molSystem in self.project.sortedMolSystems():
        for chain in molSystem.sortedChains():
          if chain is self.sourceChain:
            continue
            
          if len(chain.residues) < 2:
            continue  
            
          for molResidue in chain.molecule.sortedMolResidues():
            if molTypes.get(molResidue.molType):
              chains.append(chain)
              break
    
    return chains

  def selectResidue(self, residue, row, col):
  
    self.residue = residue
    self.updateResidueButtons()

  def enableSelected(self):
  
    if self.residue:
      residues = self.residueMatrix.currentObjects
      dict = {}
      for residue1 in self.sourceChain.residues:
        dict[self.mapping.get(residue1)] = residue1
        
      for residue in residues:
        if self.mapping.get(residue) is None:
          origResidue = self.origMapping.get(residue)
          if origResidue and (dict.get(origResidue) is None):
            self.mapping[residue] = origResidue
  
      self.updateResMappingAfter()
  
  def disableSelected(self):
  
    if self.residue:
      residues = self.residueMatrix.currentObjects
      for residue in residues:
        self.mapping[residue] = None
 
      self.updateResMappingAfter()
 
  def getResidueMapping(self):
  
    self.mapping = {}
    if self.sourceChain and self.targetChain:
      mapping, score = getChainResidueMapping(self.sourceChain, self.targetChain)
      
      for residueA, residueB in mapping:
        if residueA is not None:
          self.mapping[residueA] = residueB
    
    self.origMapping = self.mapping.copy()
    self.updateResMappingAfter()
 
  def copyChainAssignments(self):
  
    if self.mapping and showOkCancel('Confirm','OK to continue?', parent=self):
     
      experimentChains = {}
      for experiment, chain, key in self.exptChainMatrix.objectList:
        operation = self.expOperations.get(key, 1)
        chains = [self.sourceChain, self.targetChain, None]
        experimentChains[experiment] = chains[operation]
          
      progressBar = ProgressBar(self, text="Transferring residue assignments", total=len(self.mapping.keys()))

      for residueA in self.mapping:
        if residueA is not None:
          residueB = self.mapping[residueA]
          if residueB is not None:
            duplicateResidueAssignments(residueA, residueB, experimentChains=experimentChains)
            
        progressBar.increment()
      
      progressBar.destroy()
      self.updateResMappingAfter()
  
  def getTargetResidue(self, residue):
  
   if self.mapping:
     names = ['<None>']
     residues = [None, ]  
     
     residue2 = self.mapping.get(residue)
     
     dict = {}
     for residue1 in self.sourceChain.residues:
       dict[self.mapping.get(residue1)] = residue1
     
     if residue2:
       names.append(self.getResidueKey(residue2))
       residues.append(residue2)
       index = 1
     else:
       index = 0
     
     for residue1 in self.targetChain.residues:
       if not dict.has_key(residue1):
         names.append(self.getResidueKey(residue1))
         residues.append(residue1)
          
     self.targetResiduePulldown.setup(names, residues, index)
    
  def setTargetResidue(self, null):
  
    residue = self.targetResiduePulldown.getObject()
    
    if self.mapping and self.residue:
      self.mapping[self.residue] = residue
      self.updateResMappingAfter()

  def updateResMappingAfter(self, obj=None):
 
    if obj:
      pass
         
    if self.waitingChain:
      return
    else:
      self.waitingChain = True
      self.after_idle(self.updateResMapping)

  def getResonanceInfo(self, residue):
  
    resonances = getResidueResonances(residue)
    
    return len(resonances)

  def updateResMapping(self):
  
    if self.residue and (self.residue.chain is not self.sourceChain):
      self.residue = None
  
    textMatrix = []
    objectList = []
    colorMatrix = []

    if self.sourceChain:
      for residue in self.sourceChain.sortedResidues():
        targetKey = None
        target = self.mapping.get(residue)
        info2  = None
        if target:
          targetKey = self.getResidueKey(target)
          info2 = self.getResonanceInfo(target)
      
        info = self.getResonanceInfo(residue)
        
        datum = []
        datum.append(self.getResidueKey(residue))
        datum.append(info)
        datum.append(targetKey)
        datum.append(info2)
        
        if info:
          if info2 is None:
            colors = [None,None,None]
            
          elif info2 > 0:  
            colors = ['#f08080',None,'#f08080']
            
          else:
            colors = ['#a0a0d0',None,'#a0d0a0']
          
        else:    
          if info2:
            colors = [None,None,None]
          
          else:
            colors = [None,None,None]
        
        textMatrix.append(datum)
        objectList.append(residue)
        colorMatrix.append(colors)
    

    self.residueMatrix.update(textMatrix=textMatrix,
                              colorMatrix=colorMatrix,
                              objectList=objectList)
    self.updateResidueButtons()
    self.waitingChain = False
  
    # Peak List functions
    
  def changeThreshold(self, *event):                                                                
  
    #if self.restrict:
    
    self.updatePeaksAfter()


  def changeScaleEntry(self, *event):
  
    #if self.restrict:
    
    self.updatePeaksAfter()
    

  def setShowCopied(self, boolean):
  
    self.showCopied = boolean
    self.updatePeaksAfter()


  def setOverwrite(self, boolean):
  
    self.overwrite = boolean
    self.updatePeaksAfter()


  def setRestrict(self, boolean):
  
    self.restrict = boolean
    self.updatePeaksAfter()

  def setConsiderAliased(self, boolean):
  
    self.considerAliased = boolean
    self.updatePeaksAfter()

  def setFollow(self, boolean):
  
    self.follow = boolean
    
    #if self.follow:
    #  self.windowPulldown.activate()
    #else:
    #  self.windowPulldown.inactivate()


  def getPeakLists(self, isotopes=None):
  
    peakLists = []
    
    if isotopes is None:
      for experiment in self.nmrProject.sortedExperiments():
        for spectrum in experiment.sortedDataSources():
          for peakList in spectrum.sortedPeakLists():
             name = '%s:%s:%d' % (experiment.name,spectrum.name,peakList.serial)
             peakLists.append( (name, peakList) )
    
    else:
      if None in isotopes:
        isotopes.remove(None)
      
      for experiment in self.nmrProject.sortedExperiments():
        for spectrum in experiment.sortedDataSources():
          isotopes2 = getSpectrumIsotopes(spectrum)
          if None in isotopes2:
            isotopes2.remove(None)
          
          isotopes3 = isotopes2[:]
          for isotope in isotopes:
            if isotope in isotopes3:
              isotopes3.remove(isotope)
 
          if not isotopes3: # were able to match all target dims
            for peakList in spectrum.sortedPeakLists():
              name = '%s:%s:%d' % (experiment.name,spectrum.name,peakList.serial)
              peakLists.append( (name, peakList) )
   
    peakLists.sort()
    return peakLists
   
  def getDimMapping(self):
 
    windowPane = self.windowPane

    dimMapping = {}
    if windowPane and self.sourcePeakList and self.targetPeakList:
      sourceSpectrum = self.sourcePeakList.dataSource
      targetSpectrum = self.targetPeakList.dataSource
      mapping1 = getDataDimAxisMapping(sourceSpectrum, windowPane)
      mapping2 = getDataDimAxisMapping(targetSpectrum, windowPane)

      for axisLabel in mapping1.keys():
        dataDim1 = mapping1[axisLabel]
        dataDim2 = mapping2.get(axisLabel)

        if dataDim2:
          j = dataDim2.dim-1
        else:
          j = None

	dimMapping[dataDim1.dim-1] = j 

      """
      # TBD: this only works if have one SampledDataDim and is a bit of a hack
      if sourceSpectrum.numDim < targetSpectrum.numDim:
        sampledDataDim = targetSpectrum.findFirstDataDim(className='SampledDataDim')
        if sampledDataDim:
          dim = sampledDataDim.dim
          if dim < targetSpectrum.numDim:
            for key in dimMapping.keys():
              value = dimMapping[key]
              if value >= dim:
                dimMapping[key] = value-1
      """
      
    return dimMapping

  def focusOnPeak(self, peak):
  
    windowPane = self.windowPane
    
    if peak and windowPane:
      createPeakMark(peak, lineWidth=2.0)
      windowFrame = windowPane.getWindowFrame()
      windowFrame.gotoPeak(peak)  


  def updateWindows(self):
  
    index = 0
    panes = []
    names = []
    pane = self.windowPane
    
    if self.sourcePeakList:
      spectrum = self.sourcePeakList.dataSource
      activeWindows = getActiveWindows(self.sourcePeakList.root)
      for window in activeWindows:
        windowPanes = window.sortedSpectrumWindowPanes()
        
        for windowPane in windowPanes:
          if not isSpectrumInWindowPane(windowPane, spectrum):
            continue
          panes.append(windowPane)
          names.append(getWindowPaneName(windowPane))  
    
      if panes:
        if pane not in panes:
          pane = panes[0]
        
        index = panes.index(pane)   
        
      else:
        pane = None
    
      if self.windowPane is not pane:
        self.setWindowPane(pane)
    
    self.windowPulldown.setup(names, panes, index)
  
  def setWindowPane(self, pane):
  
    if self.windowPane is not pane:
      self.windowPane = pane

  def updatePeakSource(self, *opt):
 
    index     = -1
    peakList  = self.sourcePeakList
    dataList  = self.getPeakLists()
    names     = [x[0] for x in dataList]
    peakLists = [x[1] for x in dataList]
    
    if peakLists:
      if peakList in peakLists:
        index = peakLists.index(peakList)
        
      else:
        peakList = peakLists[0]
        index = 0
        
    else:
      peakList = None

    if peakList is not self.sourcePeakList:
      self.sourcePeakList = peakList
      self.sourcePeak = None
      self.updateTarget()
      self.updatePeaksAfter()
      self.updateWindows()
      
    self.sourcePulldown.setup(names, peakLists, index)


  def updateTarget(self, *opt):
    # target should always be appropriate to source

    index    = 0
    names    = []
    peakList = self.targetPeakList
    peakLists = []
    
    if self.sourcePeakList:
      isotopes = getSpectrumIsotopes(self.sourcePeakList.dataSource)
      dataList = self.getPeakLists(isotopes=isotopes)
      if dataList:
        names     = [x[0] for x in dataList]
        peakLists = [x[1] for x in dataList]
        
        if peakList not in peakLists:
          peakList = peakLists[-1]
          
        index = peakLists.index(peakList)
        
    else:
      peakList = None

    if peakList is not self.targetPeakList:
      self.targetPeakList = peakList
      self.targetPeak = None
      self.updatePeaksAfter()

    self.targetPulldown.setup(names, peakLists, index)


  def changeSourcePeakList(self, peakList):
  
    if peakList is not self.sourcePeakList:
      self.sourcePeakList = peakList
      self.sourcePeak = None
      self.targetPeak = None
      self.updateTarget()
      self.updatePeaksAfter()
      self.updateWindows()


  def changeTargetPeakList(self, peakList):
  
    if peakList is not self.targetPeakList:
      self.targetPeakList = peakList
      self.targetPeak = None
      self.updatePeaksAfter()


  def updatePeakDim(self, peakDim):
  
    self.updatePeaksAfter(peak=peakDim.peak)


  def getScaleFactorDict(self):
  
    scaleFactorDict = {}
    scaleFactorDict['1H']  = float(self.scaleEntry1H.get()  or 1.0)
    scaleFactorDict['15N'] = float(self.scaleEntry15N.get() or 1.0)
    scaleFactorDict['13C'] = float(self.scaleEntry13C.get() or 1.0)
  
    return scaleFactorDict
  
  def updatePeaksAfter(self, peak=None):
  
    if peak and peak.peakList not in (self.targetPeakList, self.sourcePeakList):
      return
  
    if self.waitingPeak:
      return
    else:
      self.waitingPeak = True
      self.after_idle(self.updatePeaks)


  def updatePeaks(self):
  
    textMatrix = []
    objectList = []
    threshold  = self.thresholdEntry.get()
    
    if self.sourcePeak and (self.sourcePeak.peakList is not self.sourcePeakList):
      self.sourcePeak = None
    
    if self.sourcePeakList and self.targetPeakList and \
         (self.sourcePeakList is not self.targetPeakList):
      scaleFactors = self.getScaleFactorDict()
      dimMapping = self.getDimMapping()

      #print "dimMapping A ", dimMapping
      #print self.sourcePeakList, self.targetPeakList, threshold, scaleFactors

      for peak in self.sourcePeakList.peaks:
        
        matches  = findShiftDistPeakMatches(peak, self.targetPeakList,
                                            threshold, scaleFactors, dimMapping=dimMapping,
                                            considerAliased=self.considerAliased)
        bestDist = None
        bestPeak = None
        barf     = False
        for dist, targetPeak in matches:
          if (bestDist is None) or (dist<bestDist):
            bestDist = dist
            bestPeak = targetPeak
            
          if not self.showCopied:
            if arePeaksAssignedSame(peak, targetPeak):
              # NB arePeaksAssignedSame returns true only if assignments
              # match exactly for common dims
	      barf = True
              break
        
        if barf:
          continue
 
        if self.restrict:
          if bestDist is None:
            # effect of restrict being true is that peaks with no
            # match are excluded
	    continue
          elif bestDist > 1.0:
            #this should never be true for getPeakMatches output
	    continue
  
        datum = []
        datum.append(peak.serial)
        datum.append(' '.join([pd.annotation or '-' for pd in peak.sortedPeakDims()]))
        datum.append(len(matches))
        datum.append(bestDist)
        
        if bestPeak:
          datum.append(' '.join([pd.annotation or '-' for pd in bestPeak.sortedPeakDims()]))
        else:
          datum.append('')
        
        textMatrix.append(datum)
        objectList.append(peak)
    
    self.sourcePeakMatrix.update(objectList=objectList, textMatrix=textMatrix)
    self.updateTargetPeaks()
    self.waitingPeak = False


  def updateTargetPeaks(self):
  
    textMatrix = []
    objectList = []
    threshold  = self.thresholdEntry.get()
    
    
    if self.sourcePeak and self.sourcePeak in self.sourcePeakMatrix.objectList:
      scaleFactors = self.getScaleFactorDict()
      dimMapping = self.getDimMapping()
      #print "dimMapping B ", dimMapping
      #print self.sourcePeakList, self.targetPeakList, threshold, scaleFactors
      matches = findShiftDistPeakMatches(self.sourcePeak, self.targetPeakList,
                                         threshold, scaleFactors, dimMapping=dimMapping,
                                         considerAliased=self.considerAliased)
      for dist, peak in matches:
        
        datum = []
        datum.append(peak.serial)
        datum.append(' '.join([pd.annotation or '-' for pd in peak.sortedPeakDims()]))
        datum.append(dist)
        
        textMatrix.append(datum)
        objectList.append(peak)

    if self.targetPeak and (self.targetPeak.peakList is not self.targetPeakList):
      self.targetPeak = None
    elif self.targetPeak not in objectList:
      self.targetPeak = None

    self.targetPeakMatrix.update(objectList=objectList, textMatrix=textMatrix)
    self.updatePeakButtons()
    
  
  def updatePeakButtons(self):
  
    if self.targetPeak:
      self.peakButtons.buttons[0].enable()

    else:
      self.peakButtons.buttons[0].disable()

    if self.sourcePeakList and self.targetPeakList \
     and self.sourcePeakList.peaks and self.targetPeakList.peaks:
      self.peakButtons.buttons[1].enable()
      self.peakButtons.buttons[2].enable()

    else:
      self.peakButtons.buttons[1].disable()
      self.peakButtons.buttons[2].disable()
    
    #if self.follow:
    #  self.windowPulldown.activate()
    #
    #else:
    #  self.windowPulldown.inactivate()
      

  def selectSourcePeak(self, peak, row, col):

    self.sourcePeak = peak
    self.updateTargetPeaks()
    
    if self.follow:
      self.focusOnPeak(peak)
    

  def selectTargetPeak(self, peak, row, col):

    self.targetPeak = peak
    self.updatePeakButtons()
    
    if self.follow:
      self.focusOnPeak(peak)

  def assignChecks(self):
  
    shiftList = self.targetPeakList.dataSource.experiment.shiftList
    if self.sourcePeakList.dataSource.experiment.shiftList is shiftList:
      msg = 'Source and target use the same shiftList. '
      msg += 'Assignment will use tolerances of the target spectrum. Continue?'
      if not showOkCancel('Warning', msg, parent=self):
        return False
    
    if self.overwrite:
      msg = 'This function will overwrite assignments in the target peak list'
      if not showOkCancel('Warning', msg, parent=self ):
        return False
    
    if not self.sourcePeakMatrix.objectList:
      showWarning('Warning', 'No peaks in source table.', parent=self)
      return False
    
    return True
    
    

  def assignSelected(self):
  
    if self.targetPeak and self.sourcePeak:
    
      threshold = self.thresholdEntry.get()
      
      tolerances, null = getIsotopeWeightedTolerances(self.targetPeak, threshold,
                                                      scaleFactorDict=self.getScaleFactorDict())
      propagatePeakAssignments([self.targetPeak,], refPeak=self.sourcePeak,
                               cleanNonRef=self.overwrite,
                               tolerances=tolerances)
      
 
  def assignSingles(self):
  
     if self.assignChecks() and self.sourcePeakList and self.targetPeakList:
    
      targets = set([])
      threshold    = self.thresholdEntry.get()
      scaleFactors = self.getScaleFactorDict()
      dimMapping   = self.getDimMapping()
      peak         = self.sourcePeakMatrix.objectList[0]
      tolerances, null = getIsotopeWeightedTolerances(peak, threshold,
                                                      scaleFactorDict=scaleFactors)
       
      for peak in self.sourcePeakMatrix.objectList:
        matches  = findShiftDistPeakMatches(peak, self.targetPeakList,
                                            threshold, scaleFactors,dimMapping=dimMapping,
                                            considerAliased=self.considerAliased)

        if len(matches) == 1:
          dist, targetPeak = matches[0]
          if targetPeak in targets:
            # If assigned this target on this run
            # First assignment may clean, but second gets merged
            cleanNonRef = False
          else:
            cleanNonRef = self.overwrite
          targets.add(targetPeak)
          
          propagatePeakAssignments([targetPeak,], refPeak=peak,
                                   cleanNonRef=cleanNonRef,
                                   tolerances=tolerances)
 
  
  def assignClosest(self):
    """Descrn: Assign closest matching target peak for all source peaks.
               Useful for estimating "minimal shift distance" changes.
    """
    if self.assignChecks() and self.sourcePeakList and self.targetPeakList:
    
      targets = set([])
      threshold    = self.thresholdEntry.get()
      scaleFactors = self.getScaleFactorDict()
      dimMapping   = self.getDimMapping()
      peak         = self.sourcePeakMatrix.objectList[0]
      tolerances, null = getIsotopeWeightedTolerances(peak, threshold,
                                                      scaleFactorDict=scaleFactors)
      noMatches  = []
     
      for peak in self.sourcePeakMatrix.objectList:
        matches  = findShiftDistPeakMatches(peak, self.targetPeakList,
                                            threshold, scaleFactors,dimMapping=dimMapping,
                                            considerAliased=self.considerAliased)
        
        if not matches:
          msg = 'No peak matches source peak %s within scaled distance threshold %3g'
          showWarning('Warning', msg % (peak.serial, threshold), parent=self)
          noMatches.append(peak)
          
        else:
          dist, targetPeak = matches[0]
          if targetPeak in targets:
            # If assigned this target on this run
            # First assignment may clean, but second gets merged
            cleanNonRef = False
          else:
            cleanNonRef = self.overwrite
            
          targets.add(targetPeak)
          propagatePeakAssignments([targetPeak,], refPeak=peak,
                                   cleanNonRef=cleanNonRef,
                                   tolerances=tolerances)
  


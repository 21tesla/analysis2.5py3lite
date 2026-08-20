
"""
======================COPYRIGHT/LICENSE START==========================

EditAssignment.py: Part of the CcpNmr Analysis program

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

from memops.universal.BlockData import cumulativeArray, arrayOfIndex

from memops.general import Implementation

from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.DataEntry import askString
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.LabelDivider import LabelDivider
from memops.gui.MessageReporter import showWarning, showYesNo, showOkCancel
from memops.gui.MessageReporter import showMulti
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix

from ccp.api.nmr import Nmr
from ccpnmr.api import Analysis

from ccpnmr.analysis.core import CouplingBasic

from ccpnmr.analysis.core.ChemicalShiftBasic import getChemAtomNmrRef, lookupAtomProbability
from ccpnmr.analysis.core.ExperimentBasic import getThroughSpaceDataDims, getOnebondExpDimRefs
from ccpnmr.analysis.core import AssignmentBasic
from ccpnmr.analysis.core.MoleculeBasic import areResonancesBound, getResidueCode
from ccpnmr.analysis.core.PeakBasic import pickPeak
from ccpnmr.analysis.core.UnitConverter import ppm2pnt
from ccpnmr.analysis.core.Util import getAnalysisDataDim, getAnalysisPeakList
from ccpnmr.analysis.frames.ResonanceFrame import ResonanceFrame
from ccpnmr.analysis.popups.BasePopup import BasePopup

PPM_MEASUREMENT_TYPES = set(['Shift','MQShift'])

class EditAssignmentPopup(BasePopup):

  """
  **Assign Resonances To Peak Dimensions**
  
  The purpose of this popup window is to control the assignment of resonances to
  a specific, picked, peak in a spectrum. Such resonances might already be
  associated with a specific atom or they could be anonymous, and only
  represented by a resonance number e.g. "[123]".
  
  The general layout of the window is such that you have two tables for each
  dimension of the cross peak, i.e. the dimensions are represented as separate
  rows. The first column of tables on the left indicate the resonances which 
  *are* currently assigned to the peak dimensions. The second column of tables
  on the right indicate the resonances that *might be* assigned to the peak
  dimensions; these are usually ordered according to chemical shift closeness
  to the peak position in that dimension.

  In normal operation to assign the dimension of a peak to a specific resonance
  by clicking on one of the resonance row from the right hand table for the
  specific peak dimension. If the resonance is not already assigned, it will
  appear in the left hand table to indicate that is now assigned. More than one
  resonance may be assigned to a peak dimension (e.g. for overlapping or
  ambiguous signals) by clicking on another resonance for the same dimension.
  Totally new, never seen before, resonances are added to a dimension with the 
  [<New>] option; the new resonance then gains a chemical shift based on this
  peak dimension position.
  
  Clicking on a resonance row in the left hand tables allows you to remove it as
  a peak assignment [Clear Dim Contrib] or change the way that the specific
  resonance is linked to atoms or residue spin systems. The light blue buttons
  toward the bottom will allow you to add specific atomic assignments, atom type
  information or tentative assignments. The light orange buttons allow you to
  remove atomic and residue/spin system assignments. Double clicking on a
  resonance row will open a pulldown menu with some of these options available.

  Various options are given to control which resonance assignment possibilities
  appear in the right hand tables. Some options show more possibilities, i.e. by
  considering wider position-to-shift tolerances or including matches outside of
  the current sweep with. Other options restrict the possibilities, usually to
  make the selection more relevant to the current situation; you may want only 
  assignments within a given molecular system (group of chains), within a single
  residue, appropriate for an isotope labelling scheme or those that obey
  covalent bond relationships of the recorded experiment.

  **Caveats & Tips**

  To change the peak that is being assigned in the window move the mouse over
  the desired peak in a spectrum window and press the "a" key (lower case).
  The experiment:spectrum identity of the peak is indicated at the top left.
  
  Only resonances that are within the chemical shift tolerances (set per
  spectrum dimension) will be shown as assignment possibilities. The chemical
  shift value used for each resonance is taken from the shift list for the
  peak's experiment.
  
  Resonance and spin system names can be set by double clicking on an assigned
  resonance row and selecting the appropriate option from the pulldown menu.
  Resonance and spin system names will be annotated in brackets with serial
  numbers, e.g. "{123:contamination}[456:aliphatic]". Names will not be
  displayed if there is any atom, residue or type information. Accordingly names
  cannot be set here under such circumstances. These names are free text notes
  and have absolutely no formal meaning with regards to assignment or molecules.

  Distances only appear in the tables if you have an appropriate structure set
  at the top and then only if there are through-space transfers in particular
  dimensions of the peak's experiment.
  
  "Assignment Groups" are only used when you wish to specify mutually exclusive
  assignment possibilities, e.g. when you have an NOE peak caused by two signals
  and only certain assignments go together. Resonances with the same group
  number go together and those with different numbers do not mix; so you would
  not get an inappropriate distance restraint for example.
  """

  def __init__(self, parent, *args, **kw):
 
    BasePopup.__init__(self, parent=parent, title="Assignment : Assignment Panel", **kw)
    self.guiParent = parent
    
  def body(self, guiFrame):

    self.geometry("500x500")

    self.peak = None
    self.peakDim = None
    self.contrib = None
    self.structure = None
    self.restrictMolSys = True
    self.showMolSys = False
    self.boundAtoms = True
    self.spectrumDict = {}
    self.aliasing = True
    self.intraResidue = False
    self.refresh = False 
    self.waitingForAtoms = False
    self.waitingForType = False
    self.waitingForTentative = False
    self.doubleTol = False
    self.boundExpDimRefs = {}
    self.labellingScheme = True
    self.minIsoFrac = None
    self.maxDim = 0
    self.table = None

    guiFrame.grid_columnconfigure(0, weight=1)
    row = 0
    peakFrame = Frame(guiFrame, grid=(row, 0), sticky='ew')
    peakFrame.expandGrid(None, 1)
    
    tipText = 'The experiment, spectrum & peak list for the current peak and which shift list it currently uses'
    self.peakLabel = Label(peakFrame, text='Peak:  Shift List:', grid=(0,0), tipText=tipText)

    tipText = 'The serial number of the selected resonance'
    self.resonanceLabel = Label(peakFrame, text='Resonance: ', grid=(0,1), tipText=tipText)
        
    label = Label(peakFrame, text=' Structure:', grid=(0,2), sticky='e') 
    tipText = 'Selects a structure ensemble from which interatomic distances are calculated'
    self.structurePulldown = PulldownList(peakFrame, grid=(0,3), tipText=tipText,
                                          callback=self.setStructure, sticky='e')
   
    utilButtons = UtilityButtonList(peakFrame, doClone=False,
                                    helpUrl=self.help_url, grid=(0,4))
    
    row += 1
    div = LabelDivider(guiFrame, text='Assignment possibilities', grid=(row,0))
    
    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
           
    assnFrame = Frame(guiFrame, grid=(row,0))
    assnFrame.grid_columnconfigure(0, weight=1, minsize=140)
    assnFrame.grid_columnconfigure(2, weight=10)
    assnFrame.grid_rowconfigure(0, weight=1)
    assnFrame.grid_rowconfigure(2, weight=1)
    assnFrame.grid_rowconfigure(4, weight=1)
    assnFrame.grid_rowconfigure(6, weight=1)
    self.assignFrame = assnFrame

    self.groupPulldown = PulldownList(self, callback=self.setGroup)
    self.assignOptPulldown = PulldownList(self, callback=self.setAssignOpt)
  
    self.peakDimPanels = []
    self.newResButtons = []
    self.mainCouplingButtons = []
    self.resonancePanels = []
    for i in range(4):
      self.addDimension(i+1)
  
    self.resonancePanels[0].grid(row=0, column=1, rowspan=2, sticky = 'nsew')
    self.resonancePanels[1].grid(row=2, column=1, rowspan=2, sticky = 'nsew')
    #self.resonancePanels[2].grid(row=4, column=1, rowspan=2, sticky = 'nsew')
    self.peakDimPanels[0].grid(row=0, column=0, sticky = 'nsew')
    self.peakDimPanels[1].grid(row=2, column=0, sticky = 'nsew')
    #self.peakDimPanels[2].grid(row=4, column=0, sticky = 'nsew')
    self.newResButtons[0].grid(row=1, column=0, sticky = 'ew')
    self.newResButtons[1].grid(row=3, column=0, sticky = 'ew')
    #self.newResButtons[2].grid(row=5, column=0, sticky = 'ew')
   
    row += 1
    div = LabelDivider(guiFrame, text='Options', grid=(row,0))
    
    row += 1
    frame = Frame(guiFrame,  grid=(row,0))
    frame.grid_columnconfigure(5, weight=1)
    
    label = Label(frame, text='Aliased Possible', grid=(0,0))
    tipText = 'Whether to the peak dimension position may be aliased; matching chemical shifts +/- a whole number of sweep widths'
    self.aliasSelect = CheckButton(frame, callback=self.setAlias,
                                   selected=True, grid=(0,1), tipText=tipText)
    
    label = Label(frame, text='Restrict Mol System', grid=(0,2))
    tipText = 'Whether to restrict assignment possibilities to only those within a single molecular system (complex)'
    self.molSysSelect = CheckButton(frame, callback=self.setRestrictMolSys,
                                    selected=True, grid=(0,3), tipText=tipText)

    label = Label(frame, text='Correlated Dims', grid=(1,0))
    tipText = 'Whether to filter the number of assignment possibilities by enforcing that for directly bound dimensions only directly bound resonances are shown'
    self.boundAtomsSelect = CheckButton(frame, callback=self.setBoundAtoms,
                                        selected=True, grid=(1,1), tipText=tipText)

    label = Label(frame, text='Double Tolerances', grid=(1,2))
    tipText = 'Whether to double the normal assignment tolerances for matching peak dim positions to resonance chemical shifts'
    self.doubleTolSelect = CheckButton(frame, callback=self.setDoubleTolerance,
                                       selected=False, grid=(1,3), tipText=tipText)

    label = Label(frame, text='Intra-residue', grid=(0,4))
    tipText = 'Whether to restrict assignment possibilities to only assignments within the same residue'
    self.intraResSelect = CheckButton(frame, callback=self.setIntraResidue,
                                      selected=False, grid=(0,5), tipText=tipText)

    label = Label(frame, text='Assignment Groups', grid=(1,4))
    tipText = 'Whether to display information to specify mutually exclusive assignment groups for superimposed/ambiguous assignments'
    self.groupsSelect = CheckButton(frame, callback=self.toggleGroups,
                                    selected=False, grid=(1,5), tipText=tipText)

    label = Label(frame, text='Show Mol System', grid=(2,0))
    tipText = 'Whether to show molecular system in the proposed resonance table'
    self.showMolSysSelect = CheckButton(frame, callback=self.setShowMolSys,
                                        selected=False, grid=(2,1), tipText=tipText)

    frame2 = Frame(frame, grid=(3,0), gridSpan=(1,6))
    frame2.expandGrid(1,4)
    
    label = Label(frame2, text='Isotope Labelling:', grid=(0,0))
    tipText = 'Selects an isotope labelling scheme, if required, to filter the assignment possibilities'
    self.labellingPulldown = PulldownList(frame2, grid=(0,1), tipText=tipText,
                                          callback=self.changeLabellingScheme)

    label = Label(frame2, text='Min Isotope Fraction:', grid=(0,2))
    tipText = 'The lower limit for spin active isotope incorporation at an atom site for a resonance to be considered'
    self.minIsoFracEntry = FloatEntry(frame2, grid=(0,3), width=7, text=0.1,
                                      returnCallback=self.updateAfterEntry,
                                      tipText=tipText)
    self.minIsoFracEntry.bind('<Leave>', self.updateAfterEntry, '+')

    row += 1
    
    commands = [self.assignAtoms,self.assignAtomType,
                self.assignAtomsTentative, self.deassignResonance,
                self.deassignSpinSystem,self.clearContrib]
    tipTexts = ['Assign the resonance selected in a peak dimension table to atoms',
                'Set the atom type of a resonance selected in a peak dimension table',
                'Assign the selected resonance to atoms in a fuzzy or putative manner',
                'Remove atomic assignments from the selected resonance',
                'Remove residue or spin system assignments from the selected resonance',
                'Disconnect the selected resonance from the peak dimension it its currently linked to']
                
    texts    = ['Assign\nResonance','Set\nType',
                'Tentative\nAssign','Deassign\nResonance',
                'Deassign\nSpin System', 'Clear Dim\nContrib']
                
    functionButtons = ButtonList(guiFrame, commands=commands,
                                 texts=texts, grid=(row,0),
                                 tipTexts=tipTexts)
                                 
    buttons = functionButtons.buttons
    self.atomsButton  = buttons[0]
    self.typesButton  = buttons[1]
    self.tentaButton  = buttons[2]
    self.deassButton  = buttons[3]
    self.deassButtonS = buttons[4]
    self.clearButton  = buttons[5]
    self.atomsButton.config(bg='#C0D0E0')
    self.typesButton.config(bg='#C0D0E0')
    self.tentaButton.config(bg='#C0D0E0') 
    self.deassButton.config(bg='#E0D0C0')
    self.deassButtonS.config(bg='#E0D0C0')
    self.clearButton.config(bg='#E0D0C0') 
  

    row += 1
    tipTexts = ['Set all resonances assigned to the peak dimensions to be in the same spin system',
                'Show any atomic connectivities on a graphical structure display. If the peak is assigned the assigned connections are used, otherwise potential possibilities are displayed',
                'Show a table of all the peaks currently assigned to the selected resonance',
                'Merge the selected resonance with the others assigned to the same peak dimension',
                'Show a table of information for the selected resonance, including the shifts of all the peaks it is linked to',
                'Predict peaks based on the possible assignments in each dimension']
    commands = [self.addSpinSystem, self.showStructConnections,
                self.showPeaks, self.mergeResonances, self.resonanceInfo, self.predictPeaks]
    texts    = ['Set Same\nSpin System','Show On\nStructure',
                'Show\nPeaks','Merge\nResonances','Resonance\nInfo', 'Predict\nPeaks']
    bottomButtons = ButtonList(guiFrame, commands=commands, texts=texts,
                               tipTexts=tipTexts, grid=(row,0))
    
    self.ssystButton = bottomButtons.buttons[0]
    self.structButton = bottomButtons.buttons[1]
    self.peaksButton = bottomButtons.buttons[2]
    self.mergeButton = bottomButtons.buttons[3]
    self.infoButton = bottomButtons.buttons[4]
    self.predictButton = bottomButtons.buttons[5]
   
    self.updateStructures() # calls updateAfter - but not always (only of struc new or changed)
    self.updateLabellingSchemes()
  
    self.curateNotifiers(self.registerNotify)
    
    self.bind('<KeyPress-Delete>', self.askClearContrib)
    self.update()
  
  def addDimension(self, dim):
  
    headingList = ['F%d' % dim]
    tipTexts = ['The peak dimension number and position',
                'The (mutual exclusion) assignment group to for the linked resonance']
    editSetCallbacks = [self.setAssignOpt, self.setGroup]
    editGetCallbacks = [self.getAssignOpt, self.getGroup]
    editWidgets      = [self.assignOptPulldown, self.groupPulldown]
    scrolledMatrix = ScrolledMatrix(self.assignFrame, initialRows=2,
                                    headingList=headingList, 
                                    callback=self.setContrib,
                                    passSelfToCallback=True,
                                    editWidgets=editWidgets, 
                                    editSetCallbacks=editSetCallbacks,
                                    editGetCallbacks=editGetCallbacks,
                                    tipTexts=tipTexts) 
                                    
    resonancePanel = ResonanceFrame(self.assignFrame, self)
                                    
    tipText = 'Make a new resonance object and link to this peak dimension'
    newResButton   = Button(self.assignFrame, text='<New>', 
                            command=lambda n=dim:self.newResonance(n),
                            pady=0, relief='groove',
                            bd=2, background='grey82',
                            tipText=tipText)

    tipText = 'Make a new resonance object and link to this coupling/splitting'
    mainCouplingButton = Button(self.assignFrame, text='<Main>', 
                                command=lambda n=dim:self.newResonance(n),
                                pady=0, relief='groove',
                                bd=2, background='grey82',
                                tipText=tipText)
                            
    self.peakDimPanels.append( scrolledMatrix )
    self.newResButtons.append( newResButton )
    self.resonancePanels.append( resonancePanel )
    self.mainCouplingButtons.append( mainCouplingButton )
    self.maxDim += 1
    
  
  def curateNotifiers(self, notifyFunc):
    
    for func in ('__init__','delete'):
      for clazz in ('ccp.nmr.Nmr.PeakDimContrib',
                    'ccp.nmr.Nmr.Resonance',
                    'ccp.nmr.Nmr.ResonanceSet',
                    'ccp.nmr.Nmr.Peak',
                    'ccp.nmr.Nmr.PeakDim',
                    'ccp.nmr.Nmr.PeakContrib',
                    'ccp.nmr.Nmr.PeakIntensity'):
        notifyFunc(self.updateAfter, clazz, func)

    # need this because notifier for resonanceGroup.removeResonance
    # is called when resonanceGroup no longer has resonance in it
    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.PeakDim', 'setAnnotation')
    
    for func in ('setPeakDimContribs','addPeakDimContrib','removePeakDimContrib'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.PeakContrib', func)
    
    for func in ('__init__','delete', 'setName',
                 'addResonance','removeResonance',
                 'setResonances','setResidue','setCcpCode'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)
    
    for func in ('removeAssignName','setAssignNames',
                 'addAssignName','setIsotopeCode',
                 'setNuclGroupType','setAtomSiteType','delete'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.Resonance', func)

    for func in ('delete','__init__','setWeight','setPossibility'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.ResidueProb', func)

    for func in ('delete','__init__'):
      notifyFunc(self.updateStructures, 'ccp.molecule.MolStructure.StructureEnsemble', func)
    
    for func in ('delete','__init__'):
      notifyFunc(self.updateLabellingSchemes, 'ccp.molecule.ChemCompLabel.LabelingScheme', func)
     


  def open(self):
  
    self.updateStructures() # calls updateAfter
    self.updateAfter()
    self.curateNotifiers(self.registerNotify) 
    BasePopup.open(self)

  def newSplitting(self, component):
    # Reduced dimensionality or JCoupling
    
    peak = self.peak
    if peak:
      project = peak.root
      self.peakDim = component.peakDim
      expDimRef = component.dataDimRef.expDimRef
      
      # TBD: Is below OK? How else to know if this is a coupling?
      # RefExpDimRef.coupledIsotopeCodes?
      if expDimRef.unit == 'Hz':
        isotopeA = expDimRef.isotopeCodes[0]
        isotopeB = expDimRef.isotopeCodes[-1]
        resonanceA = AssignmentBasic.newResonance(project, isotopeCode=isotopeA)
        resonanceB = AssignmentBasic.newResonance(project, isotopeCode=isotopeB)
        value = CouplingBasic.getPeakDimComponentSplitting(component)
        
        experiment = expDimRef.expDim.experiment
        jCouplingList = CouplingBasic.getExperimentJCouplingList(experiment)
        
        resonances = set([resonanceA, resonanceB])   
        jCoupling = jCouplingList.newJCoupling(value=value, resonances=resonances)

        self.assign(jCoupling, component)
        
      else:
        isotope = expDimRef.isotopeCodes[0]
        resonance = AssignmentBasic.newResonance(project, isotopeCode=isotope)
        self.assign(resonance, component)
  
  def mainResonanceCouplings(self, dim):
  
    peak = self.peak
    if peak:
      project = peak.root
      self.peakDim = peak.findFirstPeakDim(dim=dim)
      for cluster in peak.peakClusters:
        CouplingBasic.assignPrimaryClusterCoupling(cluster)
        # Alternatively could do only one dim, 
        # Above will do whole cluster
        #  - the whole cluster would get the same assignment in any case

  def newResonance(self, dim, component=None):
    
    peak = self.peak
    if peak:
      project = peak.root
      self.peakDim = peak.findFirstPeakDim(dim=dim)
      
      if component:
        isotope = component.dataDimRef.expDimRef.isotopeCodes[0]
      else:
        isotope = self.peakDim.dataDimRef.expDimRef.isotopeCodes[0]
      
      resonance = AssignmentBasic.newResonance(project, isotopeCode=isotope)
      self.assign(resonance, component)

  def showPeaks(self):
  
    if self.contrib:
      
      peaks = set([])
      if isinstance(self.contrib, Nmr.PeakDimContribN):
        for resonance in self.contrib.resonances:
           for contrib in resonance.peakDimContribs:
             peaks.add(contrib.peakDim.peak)

      else:
        for contrib in self.contrib.resonance.peakDimContribs:
          peaks.add(contrib.peakDim.peak)
 
      self.guiParent.viewPeaks(peaks)

  def mergeResonances(self):

    if self.contrib:
      resonances = []
      contribs = list(self.contrib.peakDim.peakDimContribs)
      for contrib in contribs:
        resonances.append(contrib.resonance)

      if len(resonances) == 1:
        msg = 'More than one resonance must be selected.'
        showWarning('Merge failed', msg, parent=self)
      elif len(resonances) > 1:
        msg = 'Are you sure you want to\nmerge %d resonances?' % len(resonances)
        
        if showOkCancel('Merge resonances', msg, parent=self):
          resonance1 = resonances[0]
          for contrib in contribs[1:]:
            for peakContrib in contrib.peakContribs:
              if len(peakContrib.peakDimContribs) < 2:
                peakContrib.delete()
            contrib.delete() # otherwise we get 2+ contribs linked to the same thing
         
          for resonance2 in resonances[1:]:
            AssignmentBasic.mergeResonances(resonance1, resonance2)

  def resonanceInfo(self):
  
     if self.contrib:
       if self.contrib.peakDimComponent:
         return
     
       self.guiParent.browseResonanceInfo(self.contrib.resonance)
  
  def predictPeaks(self):

    if not self.peak:
      return

    spectrum  = self.peak.peakList.dataSource
    numDim = spectrum.numDim
    shiftList = spectrum.experiment.shiftList
    resonanceDict = {}
    resonancesList = []
    for resonancePanel in self.resonancePanels[:numDim]:
      objectList = resonancePanel.scrolledMatrix.objectList
      if not objectList or objectList[0].className != 'Resonance':
        return

      resonances = []
      for resonance in objectList:
        shift = resonance.findFirstShift(parent=shiftList)
        if not shift:
          continue
        resonances.append(resonance)
        resonanceDict[resonance] = shift.value
      if not resonances:
        return
      
      resonancesList.append(resonances)

    peakList  = spectrum.findFirstPeakList(isSimulated=True, name='AssignPredict')
    if not peakList:
      peakList = spectrum.newPeakList(isSimulated=True, name='AssignPredict')
      analysisPeakList = getAnalysisPeakList(peakList)
      analysisPeakList.symbolStyle = '+'
      analysisPeakList.symbolColor = '#FF0000'
      analysisPeakList.textColor = '#BB0000'

    array = [len(resonances) for resonances in resonancesList]
    n, cumul = cumulativeArray(array)

    assignResToDim = AssignmentBasic.assignResToDim
    r = range(numDim)
    for i in range(n):
      indices = arrayOfIndex(i, cumul)
      resonances = [resonancesList[j][indices[j]] for j in r]
      values = [resonanceDict[resonance] for resonance in resonances]
      peak = pickPeak(peakList, values, unit=shiftList.unit, doFit=False)
      for j, peakDim in enumerate(peak.sortedPeakDims()):
        assignResToDim(peakDim,resonances[j])

  def deassignOptAtom(self):

    resonance = self.contrib.resonance
    AssignmentBasic.deassignResonance(resonance, clearAssignNames=False)
  
  def deassignOptAtomAndType(self):

    resonance = self.contrib.resonance
    AssignmentBasic.deassignResonance(resonance, clearAssignNames=True)
    if resonance.name in resonance.assignNames:
      resonance.name = None
  
  def deassignOptType(self):

    resonance = self.contrib.resonance
    if resonance.name in resonance.assignNames:
      resonance.name = None
    AssignmentBasic.assignResonanceType(resonance, None)
  
  def deassignOptSpinSystem(self):

    resonance = self.contrib.resonance
    spinSystem = resonance.resonanceGroup
    AssignmentBasic.deassignResonance(resonance, clearAssignNames=False)
    AssignmentBasic.removeSpinSystemResonance(spinSystem, resonance)

  def deassignOptCcpCode(self):

    resonance = self.contrib.resonance
    spinSystem = resonance.resonanceGroup
    for residueProb in resonance.resonanceGroup.residueProbs:
      residueProb.delete()
    AssignmentBasic.deassignResonance(resonance, clearAssignNames=False)
    AssignmentBasic.assignSpinSystemResidue(spinSystem, None)
    AssignmentBasic.assignSpinSystemType(spinSystem, None)
  
  def deassignOptResidue(self):

    resonance = self.contrib.resonance
    spinSystem = resonance.resonanceGroup
    AssignmentBasic.deassignResonance(resonance, clearAssignNames=False)
    AssignmentBasic.assignSpinSystemResidue(spinSystem, None)
  
  def deassignOptTentative(self):

    resonance = self.contrib.resonance
    spinSystem = resonance.resonanceGroup
    AssignmentBasic.deassignResonance(resonance, clearAssignNames=False)
    for residueProb in resonance.resonanceGroup.residueProbs:
      residueProb.delete()

  def deassignSpinSystem(self):
 
    if self.contrib:
      resonance = self.contrib.resonance
      spinSystem = resonance.resonanceGroup
      msg = 'What would you like to deassign?'
      texts = []
      objects = []
 
      if spinSystem:
        texts.append('Resonance from spin system')
        objects.append(self.deassignOptSpinSystem)
      
      if spinSystem and spinSystem.residue:
        texts.append('Spin system\'s residue')
        objects.append(self.deassignOptResidue)
        texts.append('Spin system\'s residue and type')
        objects.append(self.deassignOptCcpCode)
      
      elif spinSystem and spinSystem.ccpCode:
        texts.append('Spin system\'s residue type')
        objects.append(self.deassignOptCcpCode)

      elif spinSystem and spinSystem.residueProbs:
        texts.append('Tentative residue')
        objects.append(self.deassignOptTentative)
        texts.append('Tentative residue and type')
        objects.append(self.deassignOptCcpCode)
      
      if len(objects) == 1:
        objects[0]()
        return
      
      if not objects:
        return
      
      texts.append('Cancel')
      objects.append(None)

      func = showMulti('Query', msg, texts, objects, parent=self)
      if func:
        func()
        
  def deassignResonance(self):
  
    if self.contrib:
      resonance = self.contrib.resonance
      msg = 'What would you like to deassign?'
      texts = []
      objects = []
      
      if resonance.resonanceSet:
        texts.append('Atom assignment')
        objects.append(self.deassignOptAtom)
        texts.append('Atom assignment and type')
        objects.append(self.deassignOptAtomAndType)
        
      elif resonance.assignNames:
        texts.append('Atom type')
        objects.append(self.deassignOptType)
      
      if not objects:
        return
      
      if len(objects) == 1:
        objects[0]()
        return
     
      texts.append('Cancel')
      objects.append(None)

      func = showMulti('Query', msg, texts, objects, parent=self)
      if func:
        func()
      
  def changeLabellingScheme(self, scheme):
  
    if scheme is not self.labellingScheme:
      self.labellingScheme = scheme
      self.updateAfter()

  def updateLabellingSchemes(self, obj=None):
  
    index = 0
    schemes = [True, None,] + self.project.sortedLabelingSchemes()
    scheme = self.labellingScheme
    names = ['Automatic from sample', '<None>',] + ['Scheme: %s' % sc.name for sc in schemes[2:]]
  
    if scheme not in schemes:
      scheme = schemes[0]
      
    index = schemes.index(scheme)
      
    if scheme is not self.labellingScheme:
      self.labellingScheme = scheme
      self.updateAfter()
  
    self.labellingPulldown.setup(names, schemes, index)

  def setGroup(self, event):

    peakContrib = self.groupPulldown.getObject()
    contrib = self.contrib
    peakContribs = contrib.peakContribs
    
    if peakContrib is True:
      peakContrib = self.peak.newPeakContrib()
    
    if peakContrib is None:
      peakContribs0 = [pc for pc in self.peak.sortedPeakContribs() if pc.peakDimContribs]
      for peakContrib in peakContribs0:
        if peakContrib not in peakContribs:
          peakContrib.addPeakDimContrib(contrib)
        
    else:
    
      if peakContrib not in peakContribs:
        peakContrib.addPeakDimContrib(contrib)

      for peakContrib0 in peakContribs:
        if peakContrib0 is peakContrib:
          continue
      
        if len(peakContrib0.peakDimContribs) < 2: # Also clears orphans
          peakContrib0.delete()
        else:
          peakContrib0.removePeakDimContrib(contrib)
       
  def getGroup(self, contrib):
  
    index = 0
    names = []
    objects = []
    
    if self.peak:
      peakContribs = [pc for pc in self.peak.sortedPeakContribs() if pc.peakDimContribs]
      if len(peakContribs) > 1:
        names = ['*',]
        objects = [None,] 
             
      for i, peakContrib in enumerate(peakContribs):
        names.append('%d' % (i+1))
        objects.append(peakContrib)

      if contrib.peakContribs == frozenset(peakContribs):
        index = 0
      else:
        index = peakContribs.index(contrib.findFirstPeakContrib()) +1
    
      names.append('<New>')
      objects.append(True)

    # TBD allow arbitrary multiple membership

    self.groupPulldown.setup(names, objects, index)

  def setAssignOpt(self, event):
  
    if not self.contrib:
      return
  
    option = self.assignOptPulldown.getObject()
    resonance = self.contrib.resonance
    spinSystem = resonance.resonanceGroup
  
    if self.table:
      self.table.keyPressEscape()
  
    if option is 1:
      self.after_idle(self.assignAtoms)
    
    elif option is 2:
      self.assignAtomType()
      
    elif option is 3:
      self.assignAtomsTentative()
    
    elif option is 4:
      name = askString('Query', 'Enter Resonance Name:',
                       resonance.name or '', parent=self)

      if name is not None:
        name = name.strip() or None
        resonance.setName(name)
    
    elif spinSystem and (option is 5):
      name = askString('Query', 'Enter Spin System Name:',
                       spinSystem.name or '', parent=self)

      if name is not None:
        name = name.strip() or None
        spinSystem.setName(name)
 
    
  def getAssignOpt(self, contrib):

    index = 0
    names = ['<Select Option>','Assign Atom','Set Type','Tentative Assign']
    options = [0,1,2,3]
 
    resonance = self.contrib.resonance
 
    if (not resonance.assignNames) and (not resonance.resonanceSet):
      names.append('Set Resonance Name')
      options.append(4)
 
    spinSystem = resonance.resonanceGroup
    if spinSystem and (not spinSystem.ccpCode) and (not spinSystem.residue):
      names.append('Set Spin System Name')
      options.append(5)

    self.assignOptPulldown.setup(names, options, index)
  

  def setStructure(self, structure):

    if structure is not self.structure:
      self.structure = structure
      
      if self.peak:
        experiment = self.peak.peakList.dataSource.experiment
        experiment.editAssignmentPopupStructure = structure
        
      self.updateAfter()


  def getStructures(self):

    structures = []
    if self.peak:
      for molSystem in self.peak.peakList.dataSource.experiment.molSystems:
        for structure in molSystem.sortedStructureEnsembles():
          structures.append(structure)
 
    return structures

  def updateStructures(self, *object):
 
    names = []
    index = 0
    structures = self.getStructures()
    
    if self.peak:
      experiment = self.peak.peakList.dataSource.experiment
    else:
      experiment = None  

    if structures:
      if self.structure in structures:
        structure = self.structure
      else:
        structure = structures[0]
        
        if experiment and hasattr(experiment,'editAssignmentPopupStructure'):
          structureB = experiment.editAssignmentPopupStructure

          if structureB in structures:
            structure = structureB         
        
      names = [str(x.ensembleId) for x in structures]
      index = structures.index(structure)
 
    else:
      structure = None

    self.structurePulldown.setup(names, structures, index)

    if structure is not self.structure:
      self.structure = structure
      
      if experiment:
        experiment.editAssignmentPopupStructure = structure
      
      self.updateAfter()

  def haveCommonPeakContrib(self, resonance1, resonance2):
    
    peakDimContribs1 = resonance1.peakDimContribs
    if not peakDimContribs1:
      return True # not sure if this is correct

    peakDimContribs2 = resonance2.peakDimContribs
    if not peakDimContribs2:
      return True # not sure if this is correct
      
    peakContribs1 = set()
    for peakDimContrib in peakDimContribs1:
      peakContribs1 |= peakDimContrib.peakContribs
    if not peakContribs1:
      return True # not sure if this is correct
    
    peakContribs2 = set()
    for peakDimContrib in peakDimContribs2:
      peakContribs2 |= peakDimContrib.peakContribs
    if not peakContribs2:
      return True # not sure if this is correct
      
    return peakContribs1 & peakContribs2

  def showStructConnections(self):

    if self.peak:
      self.guiParent.viewStructure(self.structure)
      popup = self.guiParent.popups['view_structure']
      popup.clearConnections()
      
      spectrum = self.peak.peakList.dataSource
      dims = [dd.dim for dd in getThroughSpaceDataDims(spectrum)]

      if dims:
        assigned = {}
        for peakDim in self.peak.peakDims:
          resonances = []
         
          if peakDim.dataDimRef:
            for contrib in peakDim.peakDimContribs:
              if contrib.resonance.resonanceSet:
                resonances.append(contrib.resonance)
 
          assigned[peakDim.dim] = resonances
 
        dim1, dim2 = dims[:2]
        resonances1 = assigned[dim1]
        resonances2 = assigned[dim2]
        
        if not resonances1:
          resonances1 = self.resonancePanels[dim1-1].resonances

        if not resonances2:
          resonances2 = self.resonancePanels[dim2-1].resonances
      
        for r1 in resonances1:
          if r1:
            for r2 in resonances2:
              if r2 and self.haveCommonPeakContrib(r1, r2):
                popup.showResonancesConnection(r1,r2)

      else:
        # For non-through-space spectra there are no connections, but the atoms are highlighted anyhow.
        popup.showPeakConnection(self.peak)

  def close(self):
  
    self.curateNotifiers(self.unregisterNotify) 
    self.cancelAllWaits()
    BasePopup.close(self)

  def destroy(self):
  
    self.curateNotifiers(self.unregisterNotify) 
    BasePopup.destroy(self)


  def setIntraResidue(self, boolean):
  
    if self.intraResidue != boolean:
      self.updateAfter()
    
    self.intraResidue = boolean


  def setRestrictMolSys(self, boolean):
  
    if self.restrictMolSys != boolean:
      self.updateAfter()
    
    self.restrictMolSys = boolean


  def setShowMolSys(self, boolean):
  
    if self.showMolSys != boolean:
      self.updateAfter()
    
    self.showMolSys = boolean


  def setDoubleTolerance(self, boolean):
  
    if self.doubleTol != boolean:
      self.updateAfter()

      
    self.doubleTol = boolean

  def setBoundAtoms(self, boolean):
  
    if self.boundAtoms != boolean:
      self.updateAfter()

    self.boundAtoms = boolean


  def setAlias(self, boolean):
  
    if self.aliasing != boolean:
      self.updateAfter()
      
    self.aliasing = boolean

  def toggleGroups(self, boolean):

      self.updateAfter()
       
  def setContrib(self, object, row, col, table=None):
  
    self.table = table
    
    if object is None:
      contrib = None
      peakDim = self.peakDim
    elif object.className == 'PeakDim':
      contrib = None
      peakDim = object
    else:
      contrib = object
      peakDim = contrib.peakDim
      
    for j in range(self.maxDim):
      scrolledMatrix = self.peakDimPanels[j]
      if scrolledMatrix.currentObject is not contrib:
        scrolledMatrix.deselectAll()
      
      #if scrolledMatrix is table:
      #  if (col==0) and (contrib is self.contrib):
      #    contrib = None
      #    scrolledMatrix.deselectAll()
 
      #else:
      #  scrolledMatrix.deselectAll()
        
    self.typesButton.config(text='Set\nType')
    self.tentaButton.config(text='Tentative\nAssign')
    if contrib is None:
      self.atomsButton.config(text='Assign')
      self.atomsButton.disable()
      self.deassButton.disable()
      self.deassButtonS.disable()
      self.typesButton.disable()
      self.tentaButton.disable()
      self.clearButton.disable()
      self.mergeButton.config(text='Merge\nResonances')
      self.mergeButton.disable()
      self.infoButton.disable()
      self.peaksButton.disable()
      self.resonanceLabel.set('Resonance: ')
    
    elif contrib and isinstance(contrib, Nmr.PeakDimContribN):
      self.atomsButton.config(text='Assign')
      self.atomsButton.disable()
      self.deassButton.disable()
      self.deassButtonS.disable()
      self.typesButton.disable()
      self.tentaButton.disable()
      self.mergeButton.config(text='Merge\nResonances')
      self.mergeButton.disable()
      
      self.clearButton.enable()
      self.peaksButton.enable()
      
      res1, res2 = contrib.resonances
      self.resonanceLabel.set('Resonances: %d %d' % (res1.serial, res2.serial))
      
    else: #if contrib is not self.contrib:
      res1 = contrib.resonance
      text = AssignmentBasic.makeResonanceGuiName(res1)
       
      self.clearButton.enable()
      self.tentaButton.enable()
      self.typesButton.enable()
      
      if res1.resonanceGroup:
        self.deassButtonS.enable()
      
      if res1.resonanceSet:
        self.deassButton.enable()
      else:
        if res1.assignNames:
          self.deassButton.enable()
        else:
          self.deassButton.disable()
        
      self.atomsButton.enable()
      if res1.resonanceSet:
        self.atomsButton.config(text='Re-assign\n' + text)
      else:
        self.atomsButton.config(text='Assign\n' + text)
        
      self.mergeButton.config(text='Merge F%d\nResonances' % (peakDim.dim))
      self.resonanceLabel.set('Resonance: %d' % res1.serial)
      
      if len(contrib.peakDim.peakDimContribs) > 1:
        self.mergeButton.enable()
      else:
        self.mergeButton.disable()
        
      self.peaksButton.enable()
      self.infoButton.enable()

    self.waitingForAtoms = False
    self.waitingForType  = False
    self.waitingForTentative = False
    self.contrib = contrib
    self.peakDim = peakDim
 
  def addSpinSystem(self):
  
    if self.peak:
      AssignmentBasic.addPeakResonancesToSpinSystem([self.peak,])
 
  def askClearContrib(self, *event):
  
    if self.contrib and self.peakDim:
      name = AssignmentBasic.makeResonanceGuiName(self.contrib.resonance)
      msg = 'Clear assignment contribution of %s?' % name
      if showOkCancel('Confirm', msg , parent=self):
        self.clearContrib()

  def clearContrib(self):
  
    if self.contrib and self.peakDim:
      AssignmentBasic.clearPeakDim(self.peakDim,contrib=self.contrib)

  def assignType(self, atomSetMapping):
  
    if not (self.peakDim and self.contrib):
      return

    if atomSetMapping.mappingType == 'ambiguous':
      msg = 'Cannot set resonance atom type to an ambiguous selection' 
      showWarning('Warning', msg, parent=self)
      return
      
    if isinstance(self.contrib, Nmr.PeakDimContribN): 
      msg = 'Cannot set resonance atom type for multiple resonances' 
      showWarning('Warning', msg, parent=self)
      return
       
    resonance = self.contrib.resonance
    if resonance.isotopeCode:
      elementSymbol = resonance.isotope.chemElement.symbol
      
      if elementSymbol != atomSetMapping.elementSymbol:
        data = (resonance.isotopeCode, atomSetMapping.elementSymbol)
        msg = 'Cannot assign a %s resonance to %s atoms' % data
        showWarning('Failure', msg, parent=self)
        return
    
    if resonance.resonanceSet:
      msg = 'Resonance is already assigned to specific atoms. '
      msg += 'Deassign first if you want to set only its type.'
      showWarning('Failure', msg, parent=self)
      return
    
    atomSets  = atomSetMapping.atomSets
    AssignmentBasic.assignResonanceType(resonance,atomSets)

  def assignResidueType(self, residue):
  
    if not (self.peakDim and self.contrib):
      return

    if isinstance(self.contrib, Nmr.PeakDimContribN): 
      msg = 'Cannot set residue type for multiple resonances' 
      showWarning('Warning', msg, parent=self)
      return

    resonance = self.contrib.resonance
    spinSystem = resonance.resonanceGroup
    
    if not spinSystem:
      spinSystem = self.nmrProject.newResonanceGroup(resonances=[resonance,])

    origResidue = spinSystem.residue
    if origResidue and (origResidue.ccpCode is not residue.ccpCode):
      msg = 'Existing residue assignment will be removed.' 
      if not showOkCancel('Confirm', msg, parent=self):
        return
      AssignmentBasic.assignSpinSystemResidue(spinSystem, None)
      
    for residueProb in spinSystem.residueProbs:
          residueProb.delete()

    AssignmentBasic.assignSpinSystemType(spinSystem, residue.ccpCode,
                                         residue.molType)

  def assignResidueTentative(self, residue):
  
    if not (self.peakDim and self.contrib):
      return

    if isinstance(self.contrib, Nmr.PeakDimContribN): 
      msg = 'Cannot set tentative assignments for multiple resonances' 
      showWarning('Warning', msg, parent=self)
      return
    
    resonance = self.contrib.resonance
    spinSystem = resonance.resonanceGroup
    
    if not spinSystem:
      spinSystem = self.nmrProject.newResonanceGroup(resonances=[resonance,])
    
    origResidue = spinSystem.residue
    if origResidue and (origResidue is not residue):
      msg = 'Existing residue assignment will be removed.' 
      if not showOkCancel('Confirm', msg, parent=self):
        return
      
      
    AssignmentBasic.assignTentativeSpinSystemResidues(spinSystem, [residue,],
                                                      doWarnings=True)
                                                      
  def assignTentative(self, atomSetMapping):
  
    if not (self.peakDim and self.contrib):
      return

    if isinstance(self.contrib, Nmr.PeakDimContribN): 
      msg = 'Cannot set tentative assignments for multiple resonances' 
      showWarning('Warning', msg, parent=self)
      return   
    
    resonance = self.contrib.resonance
    if resonance.resonanceSet:
      msg = 'Continue and replace assignment with a tentative one?' 
      if not showOkCancel('Confirm', msg, parent=self):
        return
    
    atomSets  = atomSetMapping.atomSets
    AssignmentBasic.assignTentativeAtoms(atomSets, resonance)

  def assignResidue(self, residue):
  
    if not (self.peakDim and self.contrib):
      return

    if isinstance(self.contrib, Nmr.PeakDimContribN): 
      msg = 'Cannot set residue assignments for multiple resonances' 
      showWarning('Warning', msg, parent=self)
      return   
    
    resonance = self.contrib.resonance
    spinSystem = resonance.resonanceGroup
    
    if not spinSystem:
      spinSystem = self.nmrProject.newResonanceGroup(resonances=[resonance,])

    origResidue = spinSystem.residue
    if origResidue and (origResidue is not residue):
      msg = 'Existing residue assignment will be changed'
      msg += 'for all resonances in this spin system.' 
      if not showOkCancel('Confirm', msg, parent=self):
        return
    
    AssignmentBasic.assignSpinSystemResidue(spinSystem, residue, warnMerge=True)
    
    
  def assign(self, object, component=None):
    
    peakDim = self.peakDim
    if peakDim is None:
      return

    if peakDim.isDeleted:
      msg = 'Oops. Current peak was deleted.' 
      showWarning('Warning', msg, parent=self)
      return   
            
    peak = peakDim.peak
    resonances = []
    atomSets = []
    mappings = []

    if object is None:
      return
      
    elif isinstance(object, Nmr.Resonance):
      # add contrib to peakDim
      resonances.append(object)
      atomSets.append(None)

    elif isinstance(object, Nmr.JCoupling):
      if component:
        CouplingBasic.assignPeakDimComponentCoupling(component, object)
        self.contrib = None
        return
      
    elif isinstance(object, Analysis.AtomSetMapping):
      atomSetMappings = []
      isotopes = peakDim.dataDimRef.expDimRef.isotopeCodes
      elements = []
      for isotope in peakDim.dataDimRef.expDimRef.isotopeCodes:
        elements.append( re.match('\d+([A-Z]\D*)', isotope).group(1) )
      
      if object.mappingType == 'ambiguous':
 	if not hasattr(object, 'subSets'):
	  raise "Don't know what the ambiguous atomSetMapping represents: need to set atomSetMapping.subSets"
 	for subSet in object.subSets:
          if subSet.mappingType != 'nonstereo':
            atomSetMappings.append(subSet)
      else:
        atomSetMappings.append(object)
      
      for atomSetMapping in atomSetMappings:
      
        for element in elements:
          if element != atomSetMapping.elementSymbol:
            elementString = ','.join(elements)
            msg = 'Cannot assign resonances in %s dimensions to %s atoms'
            showWarning('Assignment failed', msg % (elementString,atomSetMapping.elementSymbol), parent=self)
            return

	atomSets2 = list(atomSetMapping.atomSets)
        
	atomSets.append(atomSets2)
        mappings.append(atomSetMapping)
      
        if len(atomSetMapping.resonanceSerials)<1:
          resonance = None
          resonances.append(None)
        else:
          resonance = list(atomSetMapping.resonances)[0]
          resonances.append(resonance)
          
    assignResToDim = AssignmentBasic.assignResToDim
    newResonance = AssignmentBasic.newResonance
    
    for i in range(len(resonances)):
      resonance = resonances[i]
      atomSets0 = atomSets[i]
           
      if self.contrib: 
        if resonance is None:
          resonance = self.contrib.resonance
          
        elif atomSets0 and (self.contrib.resonance is not resonance):
          
          shiftList = peak.peakList.dataSource.experiment.shiftList
          shift     = resonance.findFirstShift(parentList=shiftList)
          tolerance = getAnalysisDataDim(peakDim.dataDim).assignTolerance
 
          if self.doubleTol:
            tolerance *= 2.0
          
          residue = atomSets0[0].findFirstAtom().residue
          atomStr = ','.join(['%d%s %s' % (residue.seqCode, getResidueCode(residue), ass.name) for ass in atomSets0])
          if shift and (abs(peakDim.realValue-shift.value)>tolerance):
            resStr  = '%d' % resonance.serial
            resStr += ' at %.3f %s' % (shift.value, shiftList.unit)
            msg = 'Atom %s already assigned to another resonance with a different shift:%s'
            showWarning('Failure', msg % (atomStr,resStr))            
          
          else: # Move the assignment onto existing atom-linked resonance
            resonance2 = self.contrib.resonance
            for peakContrib in self.contrib.peakContribs:
              if len(peakContrib.peakDimContribs) < 2:
                peakContrib.delete()
            
            if not self.contrib.isDeleted:
              # Notifier overlap means this could have been 
              # deleted just before getting here
              self.contrib.delete()
            
            existingShift = resonance.findFirstShift(parentList=shiftList)
            
            if component:
              CouplingBasic.assignPeakDimComponentResonance(component, resonance,
                                                            tolerance=tolerance)
            else:
              assignResToDim(peakDim,resonance,
                             tolerance=tolerance,
                             doWarning=True)
          
            if existingShift:
              msg  = 'Resonance already present for atom selection. '
              msg += 'The existing atom-assigned resonance was used for this peak dimension. '
              msg += 'Merge resonance [%d] with the one for %s?' % (resonance2.serial, atomStr)
 
              if showYesNo('Query',msg):
                AssignmentBasic.mergeResonances(resonance, resonance2)
          
            else:
              AssignmentBasic.mergeResonances(resonance, resonance2)
              print 'NOTICE: Atom assignment %s existed in a different shift list for a different resonance.' % atomStr
              print 'NOTICE: Resonance [%d] has been merged with resonance [%d].' % (resonance2.serial, resonance.serial)
          
          return
          
      else:
        if resonance is None:
          project  = peakDim.root
          isotope = peakDim.dataDimRef.expDimRef.isotopeCodes[0]
          resonance = AssignmentBasic.newResonance(project, isotopeCode=isotope)
       
      if self.checkSpinSystem(resonance, atomSets0):
      
        if (not atomSets0) or self.checkBmrbShift(resonance, atomSets0[0], peakDim):
          tolerance = getAnalysisDataDim(peakDim.dataDim).assignTolerance
          
          if self.doubleTol:
            tolerance *= 2.0
          
          if component:
            CouplingBasic.assignPeakDimComponentResonance(component, resonance,
                                                          tolerance=tolerance)
            if atomSets0:
              AssignmentBasic.assignAtomsToRes(atomSets0,resonance)
          
          else:
            AssignmentBasic.assignPeakDim(resonance,peakDim,atomSets0,
                                          self.contrib,tolerance)
 
          self.setContrib(None, 0, 0)

  def checkBmrbShift(self, resonance, atomSet, peakDim):

    if not resonance.shifts:
      return True

    project   = resonance.root
    shiftList = peakDim.peak.peakList.dataSource.experiment.shiftList

    if shiftList:
      shift = resonance.findFirstShift(parentList=shiftList)
    
    if not shift:
      shift = resonance.findFirstShift()

    if shift:
      value = shift.value
    else:
      value = peakDim.realValue  

    residue = atomSet.findFirstAtom().residue
    molType = residue.molResidue.molType
    name    = atomSet.name
    ccpCode = residue.ccpCode
    chemAtomNmrRef = getChemAtomNmrRef(project, name, ccpCode,
                                       molType=molType)
    if chemAtomNmrRef is None:
      name = atomSet.findFirstAtom().name
      chemAtomNmrRef = getChemAtomNmrRef(project, name, ccpCode,
                                         molType=molType)

    if chemAtomNmrRef:
      typeScore  = lookupAtomProbability(project, ccpCode, name, value, molType)
      if (typeScore is not None) and (typeScore < 0.001):
        msg = 'Atom type %s is unlikely for shift %.3f. Continue?' % (name, value)
        if not showOkCancel('Warning',msg, parent=self):
          return False

    return True


  def checkSpinSystem(self, resonance, atomSets):
  
    spinSystem = resonance.resonanceGroup
    if spinSystem and atomSets:
      if not spinSystem.residue:
        return True
      if len(spinSystem.resonances) < 2:
        return True
      
      residue = atomSets[0].findFirstAtom().residue
      if spinSystem.residue is not residue:
        ccpCode = getResidueCode(residue)
        getName = AssignmentBasic.makeResonanceGuiName
        
        resText = ''
        for resonance2 in spinSystem.resonances:
          if resonance2 is not resonance:
            resText += '%s ' % getName(resonance2,fullName=False)
        
        msg = 'Move resonance %s to different residue %d%s?'
        data = (getName(resonance),residue.seqCode,ccpCode)
        if not showOkCancel('Confirm', msg % data, parent=self):
          return False                
        
        molSystems = residue.root.sortedMolSystems()
        chains = []
        for molSystem in molSystems:
          chains.extend(molSystem.sortedChains())
        if len(chains) > 1:
          if len(molSystems) > 1:
            m1 = '%s:' % spinSystem.residue.chain.molSystem.code
            m2 = '%s:' % residue.chain.molSystem.code
          else:
            m1 = m2 = ''
          c1 = '%s:' % spinSystem.residue.chain.code
          c2 = '%s:' % residue.chain.code
        else:
          m1 = m2 = c1 = c2 = ''
     
        msg = 'Re-assign all other %s%s%d%s resonances (%s) to the residue %s%s%d%s?'
        data = (m1, c1, spinSystem.residue.seqCode, getResidueCode(spinSystem.residue),
                resText, m2, c2, residue.seqCode, ccpCode) 
        if not showYesNo('Question', msg % data , self):
          AssignmentBasic.removeSpinSystemResonance(spinSystem, resonance)

    return True

  def chooseAtoms(self, atomSetMapping):
  
    if self.contrib:
      if self.contrib.isDeleted:
        self.setContrib(None, 0, 0)
        return
    
      self.peakDim = peakDim = self.contrib.peakDim
  
      if self.waitingForAtoms:
        self.assign(atomSetMapping)
        self.waitingForAtoms = False
 
      elif self.waitingForType:
        self.assignType(atomSetMapping)
        self.waitingForType = False
 
      elif self.waitingForTentative:
        self.assignTentative(atomSetMapping)
        self.waitingForTentative = False
     
      if self.peakDim is None:
        # Could have been removed by an interrupted notifier
        self.peakDim = peakDim 
      self.setContrib(self.contrib, 0, 0)

  def chooseResidue(self, residue):
  
    if self.contrib:
      if self.contrib.isDeleted:
        self.setContrib(None, 0, 0)
        return
    
      self.peakDim = peakDim = self.contrib.peakDim
  
      if self.waitingForAtoms:
        self.assignResidue(residue)
        self.waitingForAtoms = False
 
      elif self.waitingForType:
        self.assignResidueType(residue)
        self.waitingForType = False
 
      elif self.waitingForTentative:
        self.assignResidueTentative(residue)
        self.waitingForTentative = False
     
      if self.peakDim is None:
        # Could have been removed by an interrupted notifier
        self.peakDim = peakDim 
      self.setContrib(self.contrib, 0, 0)

  def assignAtoms(self, *opt):
  
    if self.waitingForAtoms:
      self.waitingForAtoms = False
      self.setContrib(self.contrib,0,0)
    else:
      if self.contrib is None:
        return

      if isinstance(self.contrib, Nmr.PeakDimContribN):
        return

      self.waitingForType  = False
      self.waitingForTentative = False
      self.waitingForAtoms = True
      self.deassButton.disable()
      self.deassButtonS.disable()
      self.mergeButton.disable()
      self.typesButton.disable()
      self.tentaButton.disable()
      self.atomsButton.config(text='Cancel\nAssign')

      self._openAtomBrowser()
     
  def assignAtomType(self):
  
    if self.waitingForType:
      self.waitingForType = False
      self.setContrib(self.contrib,0,0)
    
    else:
      self.waitingForType  = True
      self.waitingForTentative = False
      self.waitingForAtoms = False
      self.deassButton.disable()
      self.deassButtonS.disable()
      self.mergeButton.disable()
      self.atomsButton.disable()
      self.tentaButton.disable()
      self.typesButton.config(text='Cancel\nSet Type')
      
      self._openAtomBrowser()
     
  def assignAtomsTentative(self):
  
    if self.waitingForTentative:
      self.waitingForTentative = False
      self.setContrib(self.contrib,0,0)
    
    else:
      self.waitingForType = False
      self.waitingForTentative = True
      self.waitingForAtoms = False
      self.deassButton.disable()
      self.deassButtonS.disable()
      self.mergeButton.disable()
      self.typesButton.disable()
      self.atomsButton.disable()
      self.tentaButton.config(text='Cancel\nTentative')
      
      self._openAtomBrowser()

  def _openAtomBrowser(self):
   
    resonance = self.contrib.resonance
    spinSystem = resonance.resonanceGroup
    resonanceSet = resonance.resonanceSet
   
    if resonanceSet:
      atom = resonanceSet.findFirstAtomSet().findFirstAtom()
      chains = [atom.residue.chain, ]
    
    elif spinSystem and spinSystem.residue:
      chains = [spinSystem.residue.chain, ]
          
    else:
      chains = []
      
      for peakDim in self.peak.peakDims:
        for contrib in peakDim.peakDimContribs:
          if isinstance(contrib, Nmr.PeakDimContribN):
            continue
            
          if contrib is not self.contrib:
            resonance2 = contrib.resonance
            resonanceSet2 = resonance.resonanceSet
          
            if resonanceSet2:
              atom = resonanceSet2.findFirstAtomSet().findFirstAtom()
              chain = atom.residue.chain
              chains.append(chain)
              
            elif resonance2.resonanceGroup and resonance2.resonanceGroup.residue:
              chain = resonance2.resonanceGroup.residue.chain
              chains.append(chain)


    popup = self.guiParent.browseAtoms(requestor=self)
    if chains:
      if popup.chain not in chains:
        popup.setChain(chains[0])
      

  def cancelAllWaits(self):
  
    if self.waitingForType:
      self.waitingForType = False
      self.setContrib(self.contrib,0,0)
   
    if self.waitingForAtoms:
      self.waitingForAtoms = False
      self.setContrib(self.contrib,0,0)

    if self.waitingForTentative:
      self.waitingForTentative = False
      self.setContrib(self.contrib,0,0)
        
  def updateAfterEntry(self, event=None):

    if self.refresh:
      return
    
    value = self.minIsoFracEntry.get()
    if value == self.minIsoFrac:
      return 

    self.refresh = True
    self.after_idle(self.update)
          
  def updateAfter(self, object=None):

    if self.refresh:
      return

    if object and self.peak:
      name = object.className
      
      if name == 'Peak':
        if object is self.peak:
          if object.isDeleted:
            self.peak = None
          self.refresh = True
          self.after_idle(self.update)
    
      elif name == 'PeakDimContrib':
        if object.peakDim.peak is self.peak:
          self.refresh = True
          self.after_idle(self.update)
        else:
          return  

      elif name == 'PeakDimContribN':
        if object.peakDim.peak is self.peak:
          self.refresh = True
          self.after_idle(self.update)
        else:
          return  
 
      elif name == 'PeakDim':
        if object.peak is self.peak:
          self.refresh = True
          self.after_idle(self.update)

      elif name == 'PeakContrib':
        if object.peak is self.peak:
          self.refresh = True
          self.after_idle(self.update)
 
      elif (name == 'ResonanceSet') or (name == 'ResonanceGroup'):
        for peakDim in self.peak.peakDims:
          for contrib in peakDim.peakDimContribs:
            if isinstance(contrib, Nmr.PeakDimContribN):
              for resonance in contrib.resonances:
                if resonance in object.resonances:
                  self.refresh = True
                  self.after_idle(self.update)
                  return

            elif contrib.resonance in object.resonances:
              self.refresh = True
              self.after_idle(self.update)
              return
 
      elif name == 'PeakIntensity':
        if object.peak is self.peak:
          self.refresh = True
          self.after_idle(self.update)

      elif name == 'Resonance':
        for contrib in object.peakDimContribs:
          if contrib.peakDim.peak is self.peak:
            self.refresh = True
            self.after_idle(self.update)
            return
        if object.isDeleted:
          for resonancePanel in self.resonancePanels:
            if object in resonancePanel.resonances:
              self.refresh = True
              self.after_idle(self.update)
              return

      elif name == 'ResidueProb':
        resonances = object.resonanceGroup.resonances
      
        for peakDim in self.peak.peakDims:
          for contrib in peakDim.peakDimContribs:
            if contrib.resonance in resonances:
              self.refresh = True
              self.after_idle(self.update)
              return
      
    else:
      self.refresh = True
      self.after_idle(self.update)
    
  def update(self, peak=None):
    
    self.atomsButton.disable()
    self.typesButton.disable()
    self.tentaButton.disable()
    self.deassButton.disable()
    self.deassButtonS.disable()
    self.clearButton.disable()
    self.structButton.disable()
    self.mergeButton.disable()
    self.peaksButton.disable()
    self.ssystButton.disable()
    self.waitingForAtoms = False
    self.waitingForType = False
    self.contrib = None
    self.peakDim = None

    showGroups = self.groupsSelect.get()
    
    self.minIsoFrac = value = self.minIsoFracEntry.get()
    if value < 0.0:
      self.minIsoFracEntry.set(0.0)
    elif value > 1.0:
      self.minIsoFracEntry.set(1.0)
    
    resetScrollbars=False
    if peak:
      if peak.isDeleted:
        self.peak = None
      else:
        if peak is not self.peak:
          self.peak = peak
          resetScrollbars=True
    
    self.resonanceLabel.set('Resonance: ')
    mainResonancePanels = []
    
    peak = self.peak
    if peak:
      self.ssystButton.enable()
      self.updateStructures()
      if self.structure:
        self.structButton.enable()
      
      peakContribs = [pc for pc in peak.sortedPeakContribs() if pc.peakDimContribs] 
      annotation = ''
      spectrum   = peak.peakList.dataSource
      experiment = spectrum.experiment
      peakDims   = [ pd for pd in peak.sortedPeakDims() if pd.dataDimRef ]
      peakDimsB = []
      componentsB = []
      
      if self.restrictMolSys:
        molSystems = experiment.molSystems
      else:
        molSystems = []
      
      self.boundExpDimRefs = {}
      for expDimRef1, expDimRef2 in getOnebondExpDimRefs(experiment):
        if self.boundExpDimRefs.get(expDimRef1) is None:
          self.boundExpDimRefs[expDimRef1] = []

        if self.boundExpDimRefs.get(expDimRef2) is None:
          self.boundExpDimRefs[expDimRef2] = []
        
        self.boundExpDimRefs[expDimRef1].append(expDimRef2)
        self.boundExpDimRefs[expDimRef2].append(expDimRef1)
      
      
      noDisplay = self.boundAtoms or self.intraResidue or self.labellingScheme
      getName = AssignmentBasic.makeResonanceGuiName
      
      for peakContrib in peak.peakContribs:
        if not peakContrib.peakDimContribs:
          peakContrib.delete()

      i = 0 
      for peakDim in peakDims:
      
        if not peakDim.annotation:
          annotation = annotation + '-'
        else:
          annotation = annotation +peakDim.annotation

        if (peakDim.realValue is not None) and (peakDim.realValue != peakDim.value):
          posn = '%8.4f (%.4f)' % (peakDim.realValue, peakDim.value)

        elif (peakDim.realValue is not None):
          # NOTE:ED - ppm.realValue can still be None
          posn  = '%8.4f' % peakDim.value
        else:
          continue

        peakDimComponents = [None,] + peakDim.sortedPeakDimComponents()
        allContribs = peakDim.sortedPeakDimContribs()
        
        for component in peakDimComponents:
          textMatrix  = []
          colorMatrix = []
          contribs = [c for c in allContribs if c.peakDimComponent is component]
        
          isCoupling = False
          if component:
            dataDimRef = component.dataDimRef
            expDimRef = dataDimRef.expDimRef
            isotopesCode = '/'.join(expDimRef.isotopeCodes)
            
            if expDimRef.measurementType in PPM_MEASUREMENT_TYPES:
              ppm = abs(peakDim.value-peakDim.realValue)
              data = (peakDim.dim, component.serial, isotopesCode,
                      ppm, peakDim.value)
              headingList = [ 'F%d.%d %s %8.4f (%.4f)' %  data]
            
            else:
              isCoupling = True
              realPoint = ppm2pnt(peakDim.realValue, peakDim.dataDimRef)
              deltaPoints = abs(peakDim.position-realPoint)
              hz = component.scalingFactor*dataDimRef.pointToValue(deltaPoints)
              data = (peakDim.dim, component.serial, isotopesCode,
                      hz, dataDimRef.expDimRef.unit)
              headingList = [ 'J%d.%d %s %.3f %s' %  data]
            
            for contrib in contribs:
              if hasattr(contrib, 'resonances'):
                resonanceName = ','.join([getName(r) for r in contrib.resonances])
              else:
                resonanceName = getName(contrib.resonance)
                
              textMatrix.append( [resonanceName,] )
              colorMatrix.append( ['#c0d0e0'], )
          
          else:
            isotopesCode = '/'.join(peakDim.dataDimRef.expDimRef.isotopeCodes)
            headingList = [ 'F%d %s %s' % (peakDim.dim,isotopesCode,posn) ]
 
            if showGroups and peakContribs:
              headingList.append('G')
              for contrib in contribs:
                groupContribs = contrib.peakContribs
                if groupContribs:
                  if (groupContribs == frozenset(peakContribs)) \
                      and (len(groupContribs) > 1):
                    group = '*'
                  else:
                    nums = ['%d' % (peakContribs.index(pc)+1) for pc in groupContribs]
                    group = ','.join(nums)
                else:
                  group = None
 
                textMatrix.append( [getName(contrib.resonance),group] )
                colorMatrix.append(['#c0d0e0', None])
 
            else:
              for contrib in contribs:
                textMatrix.append( [getName(contrib.resonance),] )
                colorMatrix.append( ['#c0d0e0'], )
          
          dim = peakDim.dim 
          if i >= self.maxDim:
            self.addDimension(dim)
 
          if isCoupling:
            self.newResButtons[i].config(command=lambda c=component: self.newSplitting(c))
            self.mainCouplingButtons[i].config(command=lambda n=dim: self.mainResonanceCouplings(n))
          else:
            self.newResButtons[i].config(command=lambda n=dim, c=component: self.newResonance(n, c))
            self.mainCouplingButtons[i].config(command=None)
            
 
          tipTexts = ['The peak dimension number, isotope and position',
                      'The (mutual exclusion) assignment group to for the linked resonance']
          self.peakDimPanels[i].update(headingList=headingList,
                                       objectList=contribs,
                                       colorMatrix=colorMatrix,
                                       textMatrix=textMatrix,
                                       tipTexts=tipTexts)
          
 
          if isCoupling:
            self.resonancePanels[i].update(None,
                                           peakDim,
                                           self.aliasing,
                                           self.structure,
                                           component=component,
                                           resetScrollbars=resetScrollbars,
                                           molSystems=molSystems,
                                           doubleTol=self.doubleTol,
                                           showMolSystem=self.showMolSys)

            self.mainCouplingButtons[i].grid(row=2*i+1, column=0, columnspan=1, sticky='nsew') 
            self.newResButtons[i].grid(row=2*i+1, column=1, columnspan=1, sticky='nsew')

          else:
            self.resonancePanels[i].update(None,
                                           peakDim,
                                           self.aliasing,
                                           self.structure,
                                           component=component,
                                           noDisplay=noDisplay,
                                           molSystems=molSystems,
                                           doubleTol=self.doubleTol,
                                           showMolSystem=self.showMolSys)
                                         
            peakDimsB.append(peakDim)
            componentsB.append(component)
            mainResonancePanels.append(self.resonancePanels[i])
            self.newResButtons[i].grid(row=2*i+1, column=0, columnspan=2, sticky='nsew')
            self.mainCouplingButtons[i].grid_forget()
          
          self.peakDimPanels[i].grid(row=2*i, column=0, columnspan=2, sticky='nsew')
          self.assignFrame.grid_rowconfigure(2*i, weight=1)
          self.resonancePanels[i].grid(row=2*i, column=2, rowspan=2, sticky='nsew')
          self.peakDimPanels[i].deselect()
 
          i += 1
      
      if i < self.maxDim:
        for j in range(i,self.maxDim):
          self.peakDimPanels[j].grid_forget()
          self.newResButtons[j].grid_forget()
          self.resonancePanels[j].grid_forget()
          self.mainCouplingButtons[j].grid_forget()
          self.assignFrame.grid_rowconfigure(2*j, weight=0)
     
      data = (experiment.name,spectrum.name,
              peak.peakList.serial,peak.serial,
              experiment.shiftList.serial)
      
      self.peakLabel.set( 'Peak: %s:%s:%d:%d Shift List:%d' %  data)
 
    else:
      self.boundExpDimRefs = {}
      for i in range(2):
        self.assignFrame.grid_rowconfigure(2*i, weight=1)

      for i in range(self.maxDim):
        headingList = ['F%d\n' % (i+1)]
        self.peakDimPanels[i].deselect()
        self.peakDimPanels[i].update(objectList=[], textMatrix=[[],],colorMatrix=None)
        self.resonancePanels[i].update(contrib=None, peakDim=None)
        
      for i in range(2,self.maxDim):
        self.assignFrame.grid_rowconfigure(2*i, weight=0)
        self.peakDimPanels[i].grid_forget()

      self.peakLabel.set( 'Peak:  Shift List:' )
     
    if self.peak and (self.boundAtoms or self.intraResidue or self.labellingScheme):
      N = len(mainResonancePanels)
      
      resonances = []
      for i in range(N):
        resonances.append(mainResonancePanels[i].resonances or [])
 
      # Check for labelling incorporation 
      
      if self.labellingScheme is True:
        if experiment.labeledMixtures:
          labelling = experiment
        else:
          labelling = None
      else:
        labelling = self.labellingScheme  
      
      minFraction = self.minIsoFracEntry.get() or 0.01
      if labelling:
        getFrac = AssignmentBasic.getResonanceLabellingFraction
        
        for i, dimResonances in enumerate(resonances):
          filtered = []
          
          for resonance in dimResonances:
            fraction = getFrac(resonance, labelling)
            
            if fraction >= minFraction:
              filtered.append(resonance)
                  
          resonances[i] = filtered

      # Check bound atoms
      
      if self.boundAtoms:
        restricted = [set() for i in range(N)]
        restrictedDims = set()
        getFrac = AssignmentBasic.getResonancePairLabellingFraction
  
        for i in range(N-1):
          expDimRef1 = peakDimsB[i].dataDimRef.expDimRef
          boundExpDimRefs = self.boundExpDimRefs.get(expDimRef1, [])
 
          for j in range(i+1,N):
            expDimRef2 = peakDimsB[j].dataDimRef.expDimRef

            if expDimRef2 in boundExpDimRefs:
              for resonance0 in resonances[i]:
                for resonance1 in resonances[j]:
                  if resonance0 is resonance1:
                    continue
 
                  if resonance0 and resonance1:
                    if areResonancesBound(resonance0, resonance1):
                      restrictedDims.add(i)
                      restrictedDims.add(j)
                      
                      if labelling:
                        fraction = getFrac(resonance0, resonance1,
                                           labelling)
                        if fraction < minFraction:
                          continue
                    
                      restricted[i].add(resonance0)
                      restricted[j].add(resonance1)
  
        for i in restrictedDims:
          resonances[i] = restricted[i]
      
 
      if self.intraResidue:
        spinSystems = {}
        getBound = AssignmentBasic.getBoundResonances
        
        for i in range(N):
          for resonance in resonances[i]:
            spinSystem = resonance.resonanceGroup
            
            if spinSystem:
              if spinSystems.get(spinSystem) is None:
                spinSystems[spinSystem] = [set() for x in range(N)]
 
              spinSystems[spinSystem][i].add(resonance)
              
              for resonanceB in getBound(resonance):
                spinSystemB = resonanceB.resonanceGroup
                
                if spinSystemB and (spinSystemB is not spinSystem):
                  for j in range(N):
                    if j != i:
                      if resonanceB in resonances[j]:
                        spinSystems[spinSystem][j].add(resonanceB)
                       
 
        for i in range(N):
          resonances[i] = set()
          
        for spinSystem in spinSystems.keys():
          
          j = 0
          for i in range(N):
            if spinSystems[spinSystem][i]:
              j += 1
          
          if j < 2:
            continue
 
          for i in range(N):
            for resonance in spinSystems[spinSystem][i]:
              resonances[i].add(resonance)
 
      
 
      for i in range(N):
        mainResonancePanels[i].update(None, peakDimsB[i],
                                      self.aliasing,
                                      self.structure,
                                      limitedResonances=resonances[i],
                                      resetScrollbars=resetScrollbars,
                                      molSystems=molSystems,
                                      component=componentsB[i],
                                      doubleTol=self.doubleTol,
                                      showMolSystem=self.showMolSys)
       
    self.refresh = False
      

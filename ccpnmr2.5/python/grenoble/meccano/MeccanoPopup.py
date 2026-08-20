
# CCPN:MECCANO GUi Tim Stevens 17 April 2008
# Updates Tim Stevens 6 June 2008: Added link to C functions

from memops.gui.ButtonList      import ButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.Entry           import Entry
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Frame           import Frame
from memops.gui.IntEntry        import IntEntry
from memops.gui.LabelDivider    import LabelDivider
from memops.gui.Label           import Label
from memops.gui.MessageReporter import showWarning
from memops.gui.PulldownList    import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.TabbedFrame     import TabbedFrame

from memops.api import Implementation

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.AssignmentBasic import getResonanceName
from ccpnmr.analysis.core.ConstraintBasic import angleDifference
from memops.editor.Util import createDismissHelpButtonList

from ccp.lib.StructureIo import getStructureFromFile
from ccp.lib.MoleculeQuery import getLinkedResidue

try:
  from grenoble.c import Meccano
except Exception, e:
  ee = Exception('There was a problem importing Meccano module, perhaps the C code was not compiled, exception thrown by system was: %s (please contact ccpn-dev@ccpn.ac.uk for further information)' % e)
  raise ee

import time, os

# TBD
#
# 
# Copyright etc
# 
# Check GSL library installed
# 
# Import from files options? - Load RDCs, ,PhiPsi etc?
#
# Test on test data
#
# Distribution agreement

PHI_ATOM_NAMES = ['C','N','CA','C']
PSI_ATOM_NAMES = ['N','CA','C','M']

DEFAULT_ERRORS = {('H','N'):1.1,
                  ('C','N'):0.13,
                  ('C','CA'):0.21,
                  ('C','H'):0.36,
                  ('CA','HA'):2.2,
                  ('CA','CB'):0.21,
                  ('H','H'):1.0,
                  }

RAMACHANDRAN_DATABASE = os.path.join(os.path.dirname(__file__), 'data','phi_psi_database_loop_glysymm')

PDB_FORMAT = '%-6.6s%5.2d %4.4s%s%s %s%4.1d%s   %8.3f%8.3f%8.3f%6.2f%6.2f'

PI_OVER_180 = 0.017453292519943278

def meccanoPopupMacro(argServer):

  project = argServer.getProject()
  gui = argServer.parent
  
  popup = MeccanoPopup(gui, project)
  popup.open()
  
class MeccanoPopup(BasePopup):

  def __init__(self, parent, project, *args, **kw):
  
    self.alignMedium = None
    self.chain = None
    self.constraint = None
    self.constraintSet = None
    self.molSystem = None
    self.project = project
    self.run = None
    self.shiftList = None
    self.tensor = None

    
    BasePopup.__init__(self, parent=parent, title='MECCANO', *args, **kw)

    self.curateNotifiers(self.registerNotify)

  def body(self, guiFrame):
  
    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)

    options = ['Parameters','Restraints','Alignment Media & Tensors','About Meccano']
    tabbedFrame = TabbedFrame(guiFrame, options=options)
    tabbedFrame.grid(row=0, column=0, sticky='nsew')
    
    frameA, frameB, frameC, frameD = tabbedFrame.frames
    frameA.grid_columnconfigure(1, weight=1)
    frameA.grid_rowconfigure(13, weight=1)
    frameB.grid_columnconfigure(1, weight=1)
    frameB.grid_rowconfigure(1, weight=1)
    frameC.grid_columnconfigure(0, weight=1)
    frameC.grid_rowconfigure(1, weight=1)
    frameD.grid_columnconfigure(0, weight=1)
    frameD.grid_rowconfigure(0, weight=1)
    
    texts = ['Run MECCANO!']
    commands = [self.runMeccano]
    bottomButtons = createDismissHelpButtonList(guiFrame, texts=texts,
                                                commands=commands, expands=True)
    bottomButtons.grid(row=1, column=0, sticky='ew')

    if not Meccano:
      bottomButtons.buttons[0].disable()
  
    # Parameters
        
    row = 0
    label = Label(frameA, text='Calculation Run:')
    label.grid(row=row,column=0,sticky='w')
    self.runPulldown = PulldownList(frameA, callback=self.selectRun)
    self.runPulldown.grid(row=row,column=1,sticky='w')
    
    row += 1    
    label = Label(frameA, text='Shift List (for CO):')
    label.grid(row=row,column=0,sticky='w')
    self.shiftListPulldown = PulldownList(frameA, callback=self.selectShiftList)
    self.shiftListPulldown.grid(row=row,column=1,sticky='w')
           
    row += 1    
    label = Label(frameA, text='Keep Copy of Used Shifts:')
    label.grid(row=row,column=0,sticky='w')
    self.toggleCopyShifts = CheckButton(frameA)
    self.toggleCopyShifts.grid(row=row,column=1,sticky='w')
    self.toggleCopyShifts.set(True)
        
    row += 1    
    label = Label(frameA, text='Molecular System:')
    label.grid(row=row,column=0,sticky='w')
    self.molSystemPulldown = PulldownList(frameA, callback=self.selectMolSystem)
    self.molSystemPulldown.grid(row=row,column=1,sticky='w')
        
    row += 1    
    label = Label(frameA, text='Chain:')
    label.grid(row=row,column=0,sticky='w')
    self.chainPulldown = PulldownList(frameA, callback=self.selectChain)
    self.chainPulldown.grid(row=row,column=1,sticky='w')
    self.chainPulldown.bind('<Leave>', self.updateRunParams) 
        
    row += 1    
    label = Label(frameA, text='First Peptide Plane:')
    label.grid(row=row,column=0,sticky='w')
    self.firstResEntry = IntEntry(frameA, text=None, width=8)
    self.firstResEntry.grid(row=row,column=1,sticky='w')
    self.firstResEntry.bind('<Leave>', self.updateRunParams) 
        
    row += 1    
    label = Label(frameA, text='Last Peptide Plane:')
    label.grid(row=row,column=0,sticky='w')
    self.lastResEntry = IntEntry(frameA, text=None, width=8)
    self.lastResEntry.grid(row=row,column=1,sticky='w')
    self.lastResEntry.bind('<Leave>', self.updateRunParams) 
        
    row += 1    
    label = Label(frameA, text='Max Num Optimisation Steps:')
    label.grid(row=row,column=0,sticky='w')
    self.maxOptStepEntry = IntEntry(frameA, text=500, width=8)
    self.maxOptStepEntry.grid(row=row,column=1,sticky='w')
    self.maxOptStepEntry.bind('<Leave>', self.updateRunParams) 
        
    row += 1    
    label = Label(frameA, text='Num Optimisation Peptide Planes:')
    label.grid(row=row,column=0,sticky='w')
    self.numOptPlaneEntry = IntEntry(frameA, text=2, width=8)
    self.numOptPlaneEntry.grid(row=row,column=1,sticky='w')
    self.numOptPlaneEntry.bind('<Leave>', self.updateRunParams) 
        
    row += 1    
    label = Label(frameA, text='Min Num Optimisation Hits:')
    label.grid(row=row,column=0,sticky='w')
    self.numOptHitsEntry = IntEntry(frameA, text=5, width=8)
    self.numOptHitsEntry.grid(row=row,column=1,sticky='w')
    self.numOptHitsEntry.bind('<Leave>', self.updateRunParams) 

    row += 1    
    label = Label(frameA, text='File Name Prefix:')
    label.grid(row=row,column=0,sticky='w')
    self.pdbFileEntry = Entry(frameA, text='Meccano', width=8)
    self.pdbFileEntry.grid(row=row,column=1,sticky='w')
    self.pdbFileEntry.bind('<Leave>', self.updateRunParams) 
           
    row += 1    
    label = Label(frameA, text='Write Output File (.out):')
    label.grid(row=row,column=0,sticky='w')
    self.toggleWriteOutFile = CheckButton(frameA)
    self.toggleWriteOutFile.grid(row=row,column=1,sticky='w')
    self.toggleWriteOutFile.set(False)
    self.toggleWriteOutFile.bind('<Leave>', self.updateRunParams) 
           
    row += 1    
    label = Label(frameA, text='Write PDB File (.pdb):')
    label.grid(row=row,column=0,sticky='w')
    self.toggleWritePdbFile = CheckButton(frameA)
    self.toggleWritePdbFile.grid(row=row,column=1,sticky='w')
    self.toggleWritePdbFile.set(True)
    self.toggleWritePdbFile.bind('<Leave>', self.updateRunParams) 
    
    if not Meccano:
      row += 1    
      label = Label(frameA, text='The Meccano executable is not available (it needs to be compiled)', fg='red')
      label.grid(row=row,column=0,columnspan=2,sticky='w')

    # Restraints
    
    label = Label(frameB, text='Constraint Set:')
    label.grid(row=0,column=0,sticky='w')
    
    self.constraintSetPulldown = PulldownList(frameB, callback=self.selectConstraintSet)
    self.constraintSetPulldown.grid(row=0,column=1,sticky='w')
    
    self.alignMediumPulldown= PulldownList(self, callback=self.setAlignMedium)
    
    headingList = ['#','List Type','Use?','Alignment\nMedium','Num\nRestraints']
    editWidgets      = [None,None,None,self.alignMediumPulldown,None]
    editGetCallbacks = [None,None,self.toggleUseRestraints,self.getAlignMedium,None]
    editSetCallbacks = [None,None,None,self.setAlignMedium,None]
    self.restraintMatrix = ScrolledMatrix(frameB,
                                          headingList=headingList,
                                          editSetCallbacks=editSetCallbacks,
                                          editGetCallbacks=editGetCallbacks, 
                                          editWidgets=editWidgets,
                                          callback=None,
                                          multiSelect=True)
    self.restraintMatrix.grid(row=1,column=0,columnspan=2,sticky='nsew')
    
    
    # Alignment Media
    
    div = LabelDivider(frameC,text='Alignment Media')
    div.grid(row=0,column=0,sticky='ew')
    
    self.mediumNameEntry = Entry(self, returnCallback=self.setMediumName)
    self.mediumDetailsEntry = Entry(self, returnCallback=self.setMediumDetails)
    
    headingList = ['#','Name','Details','Static Tensor','Dynamic Tensor']
    editWidgets      = [None, self.mediumNameEntry, self.mediumDetailsEntry, None, None]
    editGetCallbacks = [None, self.getMediumName, self.getMediumDetails, None, None]
    editSetCallbacks = [None, self.setMediumName, self.setMediumDetails, None, None]
    self.mediaMatrix = ScrolledMatrix(frameC,
                                      headingList=headingList,
                                      editSetCallbacks=editSetCallbacks,
                                      editGetCallbacks=editGetCallbacks, 
                                      editWidgets=editWidgets,
                                      callback=self.selectAlignMedium,
                                      multiSelect=True)
                                 
    self.mediaMatrix.grid(row=1,column=0,sticky='nsew')
     
    
    texts = ['Add Alignment medium','Remove Alignment Medium']
    commands = [self.addAlignMedium,self.removeAlignMedium]
    buttonList = ButtonList(frameC, texts=texts, commands=commands, expands=True)
    buttonList.grid(row=2,column=0,sticky='nsew')
    
    self.editAxialEntry = FloatEntry(self, returnCallback=self.setAxial)
    self.editRhombicEntry = FloatEntry(self, returnCallback=self.setRhombic)
    self.editAlphaEulerEntry = FloatEntry(self, returnCallback=self.setEulerAlpha)
    self.editBetaEulerEntry = FloatEntry(self, returnCallback=self.setEulerBeta)
    self.editGammaEulerEntry = FloatEntry(self, returnCallback=self.setEulerGamma)
    
    
    div = LabelDivider(frameC,text='Alignment Tensors')
    div.grid(row=3,column=0,sticky='ew')
    
    headingList = ['Type', u'Axial (\u03B6)',u'Rhombic (\u03B7)',
                   u'Euler \u03B1',u'Euler \u03B2',u'Euler \u03B3']
    editWidgets      = [None,self.editAxialEntry,
                        self.editRhombicEntry,self.editAlphaEulerEntry,
                        self.editBetaEulerEntry,self.editGammaEulerEntry]
    editSetCallbacks = [None,self.setAxial,self.setRhombic,
                        self.setEulerAlpha,self.setEulerBeta,self.setEulerGamma]
    editGetCallbacks = [None,self.getAxial,self.getRhombic,
                        self.getEulerAlpha,self.getEulerBeta,self.getEulerGamma]
                   
    self.tensorMatrix = ScrolledMatrix(frameC, maxRows=2,
                                       headingList=headingList,
                                       editSetCallbacks=editSetCallbacks,
                                       editGetCallbacks=editGetCallbacks, 
                                       editWidgets=editWidgets,
                                       callback=self.selectTensor,
                                       multiSelect=True)
                                          
    self.tensorMatrix.grid(row=4,column=0,sticky='nsew')
    
    texts = ['Add Static Tensor','Add Dynamic Tensor','Remove Tensor']
    commands = [self.addStaticTensor,self.addDynamicTensor,self.removeTensor]
    buttonList = ButtonList(frameC,texts=texts, commands=commands, expands=True)
    buttonList.grid(row=5,column=0,sticky='ew')
       
    # About
    
    label = Label(frameD, text='About Meccano...')
    label.grid(row=0,column=0,sticky='w')
  
    #
  
    self.geometry('500x400')

    self.updateShiftLists()
    self.updateMolSystems()
    self.updateResidueRanges()
    self.updateConstraintSets()
    self.updateAlignMedia()
    self.updateRuns()
    
  def close(self):
  
    self.updateRunParams()
    
    BasePopup.close(self)
  
  def destroy(self):
  
    self.updateRunParams()
    self.curateNotifiers(self.unregisterNotify)
    
    BasePopup.destroy(self)
  
  def curateNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete'):
      notifyFunc(self.updateConstraintSetsAfter, 'ccp.nmr.NmrConstraint.NmrConstraintStore', func)
    
    for func in ('__init__', 'delete','setName','setConditionState'):
      for clazz in ('ccp.nmr.NmrConstraint.CsaConstraintList',
                    'ccp.nmr.NmrConstraint.DihedralConstraintList',
                    'ccp.nmr.NmrConstraint.DistanceConstraintList',
                    'ccp.nmr.NmrConstraint.HBondConstraintList',
                    'ccp.nmr.NmrConstraint.JCouplingConstraintList',
                    'ccp.nmr.NmrConstraint.RdcConstraintList'):
        notifyFunc(self.updateConstraintListsAfter, clazz, func)
        
    for func in ('__init__', 'delete',):
      for clazz in ('ccp.nmr.NmrConstraint.CsaConstraint',
                    'ccp.nmr.NmrConstraint.DihedralConstraint',
                    'ccp.nmr.NmrConstraint.DistanceConstraint',
                    'ccp.nmr.NmrConstraint.HBondConstraint',
                    'ccp.nmr.NmrConstraint.JCouplingConstraint',
                    'ccp.nmr.NmrConstraint.RdcConstraint'):
        notifyFunc(self.updateConstraintsAfter, clazz, func)    
    
    for func in ('__init__', 'delete'):
      notifyFunc(self.updateShiftListsAfter,'ccp.nmr.Nmr.ShiftList', func)

    for func in ('__init__', 'delete'):
      notifyFunc(self.updateMolSystemsAfter,'ccp.molecule.MolSystem.MolSystem', func)

    for func in ('__init__', 'delete'):
      notifyFunc(self.updateChainsAfter,'ccp.molecule.MolSystem.Chain', func)

    for func in ('__init__', 'delete','setDynamicAlignment',
                 'setStaticAlignment','setName','setDetails'):
      notifyFunc(self.updateAlignMediaAfter,'ccp.nmr.NmrConstraint.ConditionState', func)

  def updateAlignMediaAfter(self, alignMedium):
  
     if alignMedium.nmrConstraintStore is self.constraintSet:
       self.after_idle(self.updateAlignMedia)
 
       if alignMedium is self.alignMedium:
         self.after_idle(self.updateTensors)

  def updateConstraintSetsAfter(self, constraintSet):
  
     self.after_idle(self.updateConstraintSets)

  def updateShiftListsAfter(self, shiftList):
  
     self.after_idle(self.updateShiftLists)

  def updateMolSystemsAfter(self, molSystem):
  
     self.after_idle(self.updateMolSystems)

  def updateChainsAfter(self, chain):
  
     self.after_idle(self.updateChains)
  
  def updateConstraintsAfter(self, constraint):
  
     if constraint.parent.parent is self.constraintSet:
       self.after_idle(self.updateConstraintLists)
  
  def updateConstraintListsAfter(self, constraintList):
  
     if constraintList.parent is self.constraintSet:
       self.after_idle(self.updateConstraintLists)
    
  def runMeccano(self):
    
    #
    #
    # Input checks first
    #
    #
  
    warning = ''
    if not self.molSystem:
      warning += 'No molecular system selected\n'
    
    if not self.chain:
      warning += 'No chain selected\n'
  
    if not self.constraintSet:
      warning += 'No selected constraint set\n'  
    else:
      constraintLists = [cl for cl in self.constraintSet.constraintLists if cl.useForMeccano]  
      if not constraintLists:
        warning += 'No constraint lists selected for use\n'   
 
    first, last = self.updateResidueRanges()
    if (last-first) < 2:
      warning += 'Too few peptide planes selected\n'         
  
    if warning:
      showWarning('Cannot run MECCANO','Encountered the following problems:\n' + warning)
      return
      
    if not self.run:
      self.run = self.makeSimRun()
      
    self.updateRunParams()
    
    if self.toggleCopyShifts.get() and self.shiftList:
      shiftList = self.run.findFirstOutputMeasurementList(className='ShiftList')
      
      if not shiftList:
        shiftList = self.project.currentNmrProject.newShiftList(name='Meccano Input')
        
      self.run.setOutputMeasurementLists([shiftList,])
    
      shiftDict = {}
      for shift in shiftList.shifts:
        shiftDict[shift.resonance] = shift
      
    
      for shift in self.shiftList.shifts:
        resonance = shift.resonance
        resonanceSet = resonance.resonanceSet
        
        if resonanceSet:
          atom = resonanceSet.findFirstAtomSet().findFirstAtom()
          
          if (atom.name == 'C') and (atom.residue.molResidue.molType == 'protein'):
            shift2 = shiftDict.get(resonance)
            if shift2:
              shift2.value = shift.value
              shift2.error = shift.error
              
            else:
              shiftList.newShift(resonance=resonance, value=shift.value, error=shift.error)  
    
    #
    #
    # Accumulate data from CCPN data model & GUI
    #
    #

    # Sequence

    residues = self.chain.sortedResidues()
    residueDict = {}
    
    seqData = []
    for residue in residues:
      molResidue = residue.molResidue
      
      code1Letter = molResidue.chemComp.code1Letter
      
      if not code1Letter:
        msg = 'Encountered non-standard residue type: %s'
        showWarning('Cannot run MECCANO', msg % residue.ccpCode)
        return
     
      seqCode = residue.seqCode
      seqData.append((seqCode, code1Letter))
      residueDict[seqCode] = residue.chemCompVar.chemComp.code3Letter

    # Media, RDCs & Dihedrals

    rdcLists = []
    dihedralLists = []

    for constraintList in constraintLists:
      if constraintList.className == 'RdcConsraintList':
        if constraintList.conditionState:
          rdcLists.append(constraintList)
      
      elif constraintList.className == 'DihedralConstraintList':
        dihedralLists.append(dihedralLists)
    
    f = PI_OVER_180  
    mediaData = []
    for constraintList in rdcLists:
      medium = constraintList.conditionState
      dynamicTensor = medium.dynamicAlignment
      staticTensor = medium.staticAlignment
    
      if not (dynamicTensor or staticTensor):
        continue
    
      if dynamicTensor:
        dynamicTensorData = ['Dynamic', dynamicTensor.aAxial, dynamicTensor.aRhombic,
                             f*dynamicTensor.alpha, f*dynamicTensor.beta, f*dynamicTensor.gamma]
   
      if staticTensor:
        staticTensorData = ['Static', staticTensor.aAxial, staticTensor.aRhombic,
                            f*staticTensor.alpha, f*staticTensor.beta, f*staticTensor.gamma]
      
      if not dynamicTensor:
        dynamicTensorData = staticTensorData
      
      elif not staticTensor:
        staticTensorData = dynamicTensorData
      
      rdcData = []
      for restraint in constraintList.constraints:
        items = list(restraint.items)
        
        if len(items) != 1:
          continue
          
        resonanceA, resonanceB = [fr.resonance for fr in items[0].resonances] 
        
        resonanceSetA = resonanceA.resonanceSet
        if not resonanceSetA:
          continue
        
        resonanceSetB = resonanceB.resonanceSet
        if not resonanceSetB:
          continue
        
        resNameA = getResonanceName(resonanceA)
        resNameB = getResonanceName(resonanceB)
          
        residueA = resonanceSetA.findFirstAtomSet().findFirstAtom().residue  
        residueB = resonanceSetB.findFirstAtomSet().findFirstAtom().residue  
        
        seqA = residueA.seqCode
        seqB = residueB.seqCode
        
        value = restraint.targetValue
        error = restraint.error
        
        if error is None:
          key = [resNameA,resNameB]
          key.sort()
          sigma = DEFAULT_ERRORS.get(tuple(key), 1.0)
        
        rdcData.append([seqA, resNameA, seqB, resNameB, value, error])
      
      mediaData.append((dynamicTensorData,staticTensorData,rdcData))
      
    oneTurn = 360.0
    dihedralDict = {}
    for constraintList in dihedralLists:
      for restraint in constraintList.constraints:
        items = list(restraint.items)
        
        if len(items) != 1:
          continue
        
        item = items[0]
        resonances = [fr.resonance for fr in item.resonances] 
        
        resonanceSets = [r.resonanceSet for r in resonances]
        
        if None in resonanceSets:
          continue
          
        atoms = [rs.findFirstAtomSet().findFirstAtom() for rs in resonanceSets]  
    
        atomNames = [a.name for a in atoms]
        
        if atomNames == PHI_ATOM_NAMES:
          isPhi = True
        elif atomNames == PSI_ATOM_NAMES:
          isPhi = False
        else:
          continue
    
        residue = atoms[2].residue
        
        if residue.chain is not self.chain:
          continue
        
        if isPhi:
          residuePrev = getLinkedResidue(residue, linkCode='prev')
          
          if not residuePrev:
            continue
            
          atomC0 = residuePrev.findFirstAtom(name='C')
          atomN  = residue.findFirstAtom(name='N')
          atomCa = residue.findFirstAtom(name='CA')
          atomC  = residue.findFirstAtom(name='C')
          atoms2 = [atomC0, atomN, atomCa, atomC]
          
        else:
          residueNext = getLinkedResidue(residue, linkCode='next')
          
          if not residueNext:
            continue
          
          atomN  = residue.findFirstAtom(name='N')
          atomCa = residue.findFirstAtom(name='CA')
          atomC  = residue.findFirstAtom(name='C')
          atomN2 = residueNext.findFirstAtom(name='N')
          atoms2 = [atomN, atomCa, atomC, atomN2]
          
        if atoms != atoms2:
          continue
        
        value = item.targetValue
        error = item.error
        
        if error is None:
          upper = item.upperLimit
          lower = item.lowerLimit
          
          if (upper is not None) and (lower is not None):
            dUpper = angleDifference(value, lower, oneTurn)
            dLower = angleDifference(upper, value, oneTurn)
            error  = max(dUpper, dLower)
            
          elif lower is not None:
            error = angleDifference(value, lower, oneTurn)
            
          elif upper is not None:
            error = angleDifference(upper, value, oneTurn)
            
          else:
            error = 30.0 # Arbitrary, but sensible for TALOS, DANGLE
        
        seqCode = residue.seqCode
        if not dihedralDict.has_key(seqCode):
          dihedralDict[seqCode] = [None, None, None, None] # Phi, Psi, PhiErr, PsiErr
        
        if isPhi:
          dihedralDict[seqCode][0] = value
          dihedralDict[seqCode][2] = error
        else:
          dihedralDict[seqCode][1] = value
          dihedralDict[seqCode][3] = error
          
          
    phipsiData = []
    seqCodes = dihedralDict.keys()
    seqCodes.sort()
    
    for seqCode in seqCodes:
      data = dihedralDict[seqCode]
      
      if None not in data:
        phi, psi, phiErr, psiErr = data
        phipsiData.append((seqCode, phi, psi, phiErr, psiErr))
        
    
    # User options
    
    firstPPlaneFrag = self.firstResEntry.get() or 1
    lastPPlaneFrag  = self.lastResEntry.get() or 1
    ppNbMin         = self.numOptPlaneEntry.get() or 1
    minValueBest    = self.numOptHitsEntry.get() or 2
    maxValueBest    = self.maxOptStepEntry.get() or 5

    strucData = Meccano.runFwd(firstPPlaneFrag, lastPPlaneFrag, ppNbMin,
                               minValueBest, maxValueBest, RAMACHANDRAN_DATABASE,
                               seqData, mediaData, phipsiData)
   
    if strucData:
      fileName = 'CcpnMeccanoPdb%f.pdb' % time.time()
      fileObj = open(fileName, 'w')
 
      ch = self.chain.pdbOneLetterCode.strip()
      if not ch:
        ch = self.chain.code[0].upper()
 
        i = 1
        for atomType, resNb, x, y, z in strucData:
          resType = residueDict.get(resNb, '???')
          line = PDB_FORMAT % ('ATOM',i,'%-3s' % atomType,'',resType,ch,resNb,'',x,y,z,1.0,1.0)
 
          i += 1

      fileObj.close()
      
      ensemble = getStructureFromFile(self.molSystem, fileName)
      
      
      if not self.toggleWritePdbFile.get():
        os.unlink(fileName)
         
      self.run.outputEnsemble = ensemble
      self.run = None
    
    self.updateRuns()

  def getMediumName(self, alignMedium):
  
    self.mediumNameEntry.set(alignMedium.name)
    
  def getMediumDetails(self, alignMedium):
  
    self.mediumDetailsEntry.set(alignMedium.details)
    
  def setMediumName(self, event): 

    value = self.mediumNameEntry.get()
    self.alignMedium.name = value or None
    
  def setMediumDetails(self, event): 

    value = self.mediumDetailsEntry.get()
    self.alignMedium.details = value or None
      
  def setAlignMedium(self, alignMedium):

    if self.constraintSet:
      self.constraintSet.conditionState = alignMedium
    
  def getAlignMedium(self, constraintList):

    media = self.getAlignmentMedia()
    names = [am.name for am in media]
    
    if constraintList.conditionState in media:
      index = media.index(constraintList.conditionState)
    else:
      index = 0
  
    self.alignMediumPulldown.setup(names, media, index)

  def toggleUseRestraints(self, constraintList):
 
    bool = constraintList.useForMeccano
    bool = not bool
 
    if bool and (not constraintList.conditionState) \
     and (constraintList.className == 'RdcConsraintList'):
      msg = 'Cannot use RDC restraint list for Meccano '
      msg += 'unless it is first associated with an amigment medium'
      showWarning('Warning', msg, parent=self)
    else:
      constraintList.useForMeccano = bool
      
    self.updateConstraintLists()
 
  def addStaticTensor(self):
  
    if self.alignMedium:
      tensor = Implementation.SymmTracelessMatrix(aAxial=0.0,aRhombic=0.0,
                                                  alpha=0.0,beta=0.0,
                                                  gamma=0.0)
      

      self.alignMedium.staticAlignment = tensor
 
      self.updateAlignMediaAfter(self.alignMedium)
 
  def addDynamicTensor(self):
  
    if self.alignMedium:
      tensor = Implementation.SymmTracelessMatrix(aAxial=0.0,aRhombic=0.0,
                                                  alpha=0.0,beta=0.0,
                                                  gamma=0.0)
      self.alignMedium.dynamicAlignment = tensor
      
      self.updateAlignMediaAfter(self.alignMedium)
  
  def removeTensor(self):
  
    if self.alignMedium and self.tensor:
      if self.tensor is self.alignMedium.dynamicAlignment:
        self.alignMedium.dynamicAlignment = None
        
      elif self.tensor is self.alignMedium.staticAlignment:
        self.alignMedium.staticAlignment = None
      
      self.updateAlignMediaAfter(self.alignMedium)
        
  def addAlignMedium(self):
  
    if self.constraintSet:
      medium = self.constraintSet.newConditionState()
      medium.name = 'Align Medium %d' % medium.serial   
  
  def removeAlignMedium(self):

    if self.alignMedium:
      self.alignMedium.delete()

  def updateTensor(self, aAxial=None, aRhombic=None, alpha=None, beta=None, gamma=None):
  
    aAxial   = aAxial   or self.tensor.aAxial
    aRhombic = aRhombic or self.tensor.aRhombic
    alpha    = alpha    or self.tensor.alpha
    beta     = beta     or self.tensor.beta
    gamma    = gamma    or self.tensor.gamma
    
    tensor = Implementation.SymmTracelessMatrix(aAxial=aAxial,
                                                aRhombic=aRhombic,
                                                alpha=alpha,beta=beta,
                                                gamma=gamma)
 
    if self.alignMedium:
      if self.tensor is self.alignMedium.dynamicAlignment:
        self.alignMedium.dynamicAlignment = tensor
        
      elif self.tensor is self.alignMedium.staticAlignment:
        self.alignMedium.staticAlignment = tensor
 
    self.tensor = tensor
 
  def setAxial(self, event): 
  
    value = self.editAxialEntry.get()
    self.updateTensor(aAxial=value)
    self.updateTensors()
    
  def setRhombic(self,  event): 
  
    value = self.editRhombicEntry.get()
    self.updateTensor(aRhombic=value)
    self.updateTensors()
   
  def setEulerAlpha(self,  event): 
  
    value = self.editAlphaEulerEntry.get()
    self.updateTensor(alpha=value)
    self.updateTensors()
    
  def setEulerBeta(self,  event): 
  
    value = self.editBetaEulerEntry.get()
    self.updateTensor(beta=value)
    self.updateTensors()
    
  def setEulerGamma(self,  event): 
  
    value = self.editGammaEulerEntry.get()
    self.updateTensor(gamma=value)
    self.updateTensors()

  def getAxial(self, tensor): 
  
    value = tensor.aAxial
    self.editAxialEntry.set(value)
     
  def getRhombic(self, tensor): 
  
    value = tensor.aRhombic
    self.editRhombicEntry.set(value)
     
  def getEulerAlpha(self, tensor): 
  
    value = tensor.alpha
    self.editAlphaEulerEntry.set(value)
     
  def getEulerBeta(self, tensor): 
  
    value = tensor.beta
    self.editBetaEulerEntry.set(value)
     
  def getEulerGamma(self, tensor): 
  
    value = tensor.gamma
    self.editGammaEulerEntry.set(value)
 
  def selectTensor(self, tensor, row, col):
  
    self.tensor = tensor
 
  def selectAlignMedium(self, alignMedium, row, col):
  
    self.alignMedium = alignMedium
    self.updateTensors()
 
  def getAlignmentMedia(self):
  
    if self.constraintSet:
      return self.constraintSet.sortedConditionStates()
    else:
      return []
 
  def updateAlignMedia(self):

    textMatrix = []
    objectList = []
    
    if self.constraintSet:
      objectList = self.getAlignmentMedia()
        
      for conditionState in objectList:
         
         staticTensor = None
         dyamicTensor = None
         tensor = conditionState.dynamicAlignment
         if tensor:
           vals = (tensor.aAxial, tensor.aRhombic, tensor.alpha, tensor.beta, tensor.gamma)
           dyamicTensor = u'\u03B6:%.3f \u03B7:%.3f \u03B1:%.3f \u03B2:%.3f \u03B3:%.3f ' % vals

         tensor = conditionState.staticAlignment
         if tensor:
           vals = (tensor.aAxial, tensor.aRhombic, tensor.alpha, tensor.beta, tensor.gamma)
           staticTensor = u'\u03B6:%.3f \u03B7:%.3f \u03B1:%.3f \u03B2:%.3f \u03B3:%.3f ' % vals
           
         datum = [conditionState.serial,
                  conditionState.name,
                  conditionState.details,
                  dyamicTensor,
                  staticTensor]
         textMatrix.append(datum)

         if dyamicTensor or staticTensor:
           if not self.alignMedium:
             self.alignMedium = conditionState
   
    self.mediaMatrix.update(textMatrix=textMatrix, objectList=objectList)
      
    if self.alignMedium:
      self.mediaMatrix.selectObject(self.alignMedium)
 
  def updateTensors(self):
    
    textMatrix = []
    objectList = []
    conditionState = self.alignMedium
    
    if conditionState:
    
      tensor = conditionState.dynamicAlignment
      if tensor:
        datum = ['Dynamic', tensor.aAxial, tensor.aRhombic,
                 tensor.alpha, tensor.beta, tensor.gamma]
        textMatrix.append(datum)
        objectList.append(tensor)
      
      tensor = conditionState.staticAlignment
      if tensor:
        datum = ['Static', tensor.aAxial, tensor.aRhombic,
                 tensor.alpha, tensor.beta, tensor.gamma]
        textMatrix.append(datum)
        objectList.append(tensor)

    self.tensorMatrix.update(textMatrix=textMatrix, objectList=objectList)  
 
  def getMolSystems(self):
  
    molSystems = []
    
    for molSystem in self.project.sortedMolSystems():
      if molSystem.chains:
        molSystems.append(molSystem)
        
    return molSystems    
     
     
  def updateMolSystems(self, *notifyObj):
  
    index = 0
    names = []
    
    molSystems = self.getMolSystems()
    molSystem = self.molSystem
    
    if molSystems:
      if molSystem not in molSystems:
        molSystem = molSystems[0]
      
      index = molSystems.index(molSystem)
      names = [ms.code for ms in molSystems] 
    
    else:
      self.molSystem = None
    
    if self.molSystem is not molSystem:
      if self.run:
        self.run.molSystem = molSystem
        
      self.molSystem = molSystem
      self.updateChains()
    
    self.molSystemPulldown.setup(texts=names, objects=molSystems, index=index)
    
  
  def selectMolSystem(self, molSystem):
  
    if self.molSystem is not molSystem:
      if self.run:
        self.run.molSystem = molSystem
        
      self.molSystem = molSystem
      self.updateChains()
     
     
  def updateChains(self, *notifyObj):
  
    index = 0
    names = []
    chains = []
    chain = self.chain
    
    if self.molSystem:
      chains = [c for c in self.molSystem.sortedChains() if c.residues]
    
    if chains: 
      if chain not in chains:
        chain = chains[0]
        
      index = chains.index(chain)
      names = [c.code for c in chains]
    
    if chain is not self.chain:
      self.chain = chain
      self.updateResidueRanges()
    
    self.chainPulldown.setup(texts=names, objects=chains, index=index)
     
  def selectChain(self, chain):
  
    if chain is not self.chain:
      self.chain = chain
      self.updateResidueRanges()

  def updateResidueRanges(self):
  
    first = self.firstResEntry.get()
    last = self.lastResEntry.get()
    
    if self.chain:
      residues = self.chain.sortedResidues()
      firstSeq = residues[0].seqCode
      lastSeq = residues[-2].seqCode
      
      if first < firstSeq:
        first = firstSeq
     
      if last == first:
        last = lastSeq
      
      elif last > lastSeq:
        last = lastSeq
        
      if first > last:
        last, first = first, last    
      
  
    self.firstResEntry.set(first)
    self.lastResEntry.set(last)
 
    return first, last

  def getConstraintSets(self):
  
    constraintSets = []
    nmrProject = self.project.currentNmrProject
    
    for constraintSet in nmrProject.sortedNmrConstraintStores():
      for constraintList in constraintSet.constraintLists:
        if constraintList.className not in ('ChemShiftConstraintList','somethingElse'):
          constraintSets.append(constraintSet)
          break
    
    return constraintSets
  
  def updateConstraintSets(self, *notifyObj):

    index = 0
    names = []
    constraintSets = self.getConstraintSets()
    constraintSet = self.constraintSet

    if constraintSets:
      if constraintSet not in constraintSets:
        constraintSet = constraintSets[0]
        
      index = constraintSets.index(constraintSet)
      names = ['%d' % cs.serial for cs in constraintSets]
    
    if constraintSet is not self.constraintSet:
      if self.run:
        self.run.inputConstraintStore = constraintSet
    
      self.constraintSet = constraintSet
      self.updateConstraintLists()    
    
    self.constraintSetPulldown.setup(texts=names, objects=constraintSets, index=index)

  def selectConstraintSet(self, constraintSet):
  
    if self.constraintSet is not constraintSet:
      if self.run:
        self.run.inputConstraintStore = constraintSet
        
      self.constraintSet = constraintSet
      self.updateConstraintLists() 
  
  def getConstraintLists(self):
  
    constraintLists = []
    if self.constraintSet:
      for constraintList in self.constraintSet.sortedConstraintLists():
        if constraintList.className not in ('ChemShiftConstraintList','somethingElse'):
          constraintLists.append(constraintList)
  
    return constraintLists
    
  def updateConstraintLists(self, *notifyObj):
  
    textMatrix = []
    objectList = self.getConstraintLists()
    
    for constraintList in objectList:
      if not hasattr(constraintList, 'useForMeccano'):
        if constraintList.conditionState \
         or (constraintList.className != 'RdcConstraintList'):
          bool = True
        else:
          bool = False  
      
        constraintList.useForMeccano = bool
    
      if constraintList.conditionState:
        alignMedium = constraintList.conditionState.name
      else:
        alignMedium = None
    
      datum = [constraintList.serial,
               constraintList.className[:-14],
               constraintList.useForMeccano and 'Yes' or 'No',
               alignMedium,
               len(constraintList.constraints)]
  
      textMatrix.append(datum)
  
    self.restraintMatrix.update(textMatrix=textMatrix, objectList=objectList)
  
  def selectConstraint(self, obj, row, column):
  
    if self.constraint is not obj:
      self.constraint = obj

  def getSimStore(self):
  
    simStore = self.project.findFirstNmrSimStore(name='meccano')
    if not simStore:
      simStore = self.project.newNmrSimStore(name='meccano')
    
    return simStore

  def getRuns(self):
  
    runs = [None, ]
    
    if self.molSystem and self.constraintSet:
      simStore = self.getSimStore()
      runs += simStore.sortedRuns()
    
    return runs
  
  def updateRuns(self, *notifyObj):
  
    index = 0
    names = ['<New>']
    runs = self.getRuns()
    run = self.run

    if runs:
      if run not in runs:
        run = runs[0]
    
      index = runs.index(run)
      names += [r.serial for r in runs if r]
    
    if run is not self.run:
      self.updateConstraintSets()
      self.updateMolSystems()
      self.updateShiftLists()
    
    self.runPulldown.setup(names, runs, index)
  
  def updateRunParams(self, event=None):
  
    if self.run and self.molSystem and self.constraintSet:
      simRun = self.run
      
      simRun.inputConstraintStore = self.constraintSet
      simRun.molSystem = self.molSystem
      
      if self.shiftList:
        simRun.setInputMeasurementLists([self.shiftList,])                 

      simRun.newRunParameter(code='FirstPepPlane',id=1, 
                             intValue=self.firstResEntry.get() or 0)
      simRun.newRunParameter(code='LastPepPlane' ,id=1, 
                             intValue=self.lastResEntry.get() or 0)
      simRun.newRunParameter(code='MaxOptSteps',  id=1, 
                             intValue=self.maxOptStepEntry.get() or 0)
      simRun.newRunParameter(code='NumOptPlanes', id=1, 
                             intValue=self.numOptPlaneEntry.get() or 0)
      simRun.newRunParameter(code='MinOptHits',   id=1, 
                             intValue=self.numOptHitsEntry.get() or 0)
  
  def makeSimRun(self, template=None):

    simStore = self.getSimStore()
   
    if template:
      molSystem = template.molSystem
      constraintSet = template.inputConstraintStore
      shiftList = template.findFirstInputMeasurementList(className='ShiftList')
      protocol = template.annealProtocol
    else:
      molSystem = self.molSystem
      constraintSet = self.constraintSet
      shiftList = self.shiftList
      protocol = self.annealProtocol
    
    simRun = simStore.newRun(annealProtocol=protocol,
                             molSystem=molSystem,
                             inputConstraintStore=protocol)
    
    if shiftList:
      simRun.addInputMeasurementList(shiftList)                 
    
    if template:
      for param in template.runParameters:
        simRun.newRunParameter(code=param.code,
                               id=param.id,
                               booleanValue=param.booleanValue,
                               floatValue=param.floatValue,
                               intValue=param.intValue,
                               textValue=param.textValue)
    else:
       if self.chain:
         chainCode = self.chain.code
       else:
         chaincode = 'A'
    
       simRun.newRunParameter(code='FirstPepPlane',id=1, 
                              intValue=self.firstResEntry.get() or 0)
       simRun.newRunParameter(code='LastPepPlane' ,id=1, 
                              intValue=self.lastResEntry.get() or 0)
       simRun.newRunParameter(code='MaxOptSteps',  id=1, 
                              intValue=self.maxOptStepEntry.get() or 0)
       simRun.newRunParameter(code='NumOptPlanes', id=1, 
                              intValue=self.numOptPlaneEntry.get() or 0)
       simRun.newRunParameter(code='MinOptHits',   id=1, 
                              intValue=self.numOptHitsEntry.get() or 0)
       simRun.newRunParameter(code='ChainCode',    id=1,
                              textValue=chainCode)

    return simRun
 
  def selectRun(self, simRun):
  
    if self.run is not simRun:
      if simRun:
        if simRun.outputEnsemble:
          msg  = 'Selected run has already been used to generate a structure.'
          msg += 'A new run will be setup based on the selection.'
          showWarning('Warning',msg)
          simRun = self.makeSimRun(template=simRun)
 
        molSystem = simRun.molSystem
        constraintSet = simRun.inputConstraintStore
        shiftList = simRun.findFirstInputMeasurementList(className='ShiftList')
 
        if molSystem and (self.molSystem is not molSystem):
          self.molSystem = molSystem
          self.updateMolSystems()
 
        if constraintSet and (self.constraintSet is not constraintSet):
          self.constraintSet = constraintSet
          self.updateConstraintSets()
    
        if shiftList and (self.shiftList is not shiftList): 
          self.shiftList = shiftList
          self.updateShiftLists()
        
        firstPepPlane = simRun.findFirstrunParameter(code='FirstPepPlane')
        lastPepPlane = simRun.findFirstrunParameter(code='LastPepPlane')
        maxOptSteps = simRun.findFirstrunParameter(code='MaxOptSteps')
        numOptPlanes = simRun.findFirstrunParameter(code='NumOptPlanes')
        minOptHits = simRun.findFirstrunParameter(code='MinOptHits')
        chainCode = simRun.findFirstrunParameter(code='ChainCode')

        
        if firstPepPlane is not None:
          self.firstResEntry.set(firstPepPlane.intValue or 0)
        
        if lastPepPlane is not None:
          self.lastResEntry.set(lastPepPlane.intValue or 0)
          
        if maxOptSteps is not None:
          self.maxOptStepEntry.set(maxOptSteps.intValue or 0)
            
        if numOptPlanes is not None:
          self.numOptPlaneEntry.set(numOptPlanes.intValue or 0)
          
        if minOptHits is not None:
          self.numOptHitsEntry.set(minOptHits.intValue or 0)
           
        if chainCode is not None:
          chainCode = chainCode.textValue or 'A'
          if self.molSystem:
            chain = self.molSystem.findFirsChain(code=chainCode)
 
            if chain:
              self.selectChain(chain) 
                 
      self.run = simRun

  def updateShiftLists(self, *notifyObj):
    
    index = 0
    names = []
    nmrProject = self.project.currentNmrProject
    
    shiftLists = nmrProject.findAllMeasurementLists(className='ShiftList')
    shiftLists = [(sl.serial, sl) for sl in shiftLists]
    shiftLists.sort()
    shiftLists = [x[1] for x in shiftLists]
    
    shiftList = self.shiftList

    if shiftLists:
      if shiftList not in shiftLists:
        shiftList = shiftLists[0]
    
      index = shiftLists.index(shiftList)
      names = ['%d'% sl.serial for sl in shiftLists]
    
    if shiftList is not self.shiftList:
      if self.run:
        self.run.setInputMeasurementLists([shiftList])
  
    self.shiftListPulldown.setup(texts=names, objects=shiftLists, index=index)

  def selectShiftList(self, shiftList):
  
    if shiftList is not self.shiftList:
      if self.run:
        self.run.setInputMeasurementLists([shiftList])
      
      self.shiftList = shiftList                   
        
      

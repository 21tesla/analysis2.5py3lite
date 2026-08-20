"""
======================COPYRIGHT/LICENSE START==========================

HcloudsMdPopup.py: Part of the CcpNmr Clouds program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
import os
import Tkinter

from memops.general import Implementation
from ccp.api.nmr import Nmr

from memops.editor.BasePopup    import BasePopup
from memops.gui.ButtonList      import ButtonList
from memops.gui.IntEntry        import IntEntry
from memops.gui.Entry           import Entry
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Label           import Label
from memops.gui.LabelFrame      import LabelFrame
from memops.gui.MessageReporter import showInfo
from memops.gui.PulldownMenu    import PulldownMenu
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.Util            import createDismissHelpButtonList

from ccpnmr.analysis.core.AssignmentBasic import newResonance, assignAtomsToRes
from ccpnmr.analysis.core.ExperimentBasic import getSpectraByType

from ccpnmr.clouds.HydrogenDynamics import *

class HcloudsMdPopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.guiParent = parent
    self.project   = parent.getProject()
    self.waiting   = 0
    
    self.constraintSet = None
    self.constrLists = [None] * 4
    self.numClouds    = 100
    self.filePrefix   = 't_intra_'
    self.cloudsFiles  = []
    self.adcAtomTypes = 'HN'
    
    # step num, initial temp, final temp, cooling steps, MD steps, MD tau, rep scale
    self.coolingScheme = []
    self.coolingScheme.append([ 1,    1,    1,  3,  500, 0.001, 0])
    self.coolingScheme.append([ 2,80000, 4000, 19, 1000, 0.001, 0])
    self.coolingScheme.append([ 3, 4000,    1,  5,  500, 0.001, 0])
    self.coolingScheme.append([ 4,15000,    1,  3, 1000, 0.001, 0])
    self.coolingScheme.append([ 5,    1,    1,  5,  500, 0.001, 0])
    self.coolingScheme.append([ 6, 8000,    1,  3, 1000, 0.001, 0])
    self.coolingScheme.append([ 7,    1,    1,  5,  500, 0.001, 0])
    self.coolingScheme.append([ 8, 3000,   25, 60, 2500, 0.001, 1])
    self.coolingScheme.append([ 9,   25,   25,  1, 7500, 0.001, 1])
    self.coolingScheme.append([10,   10,   10,  1, 7500, 0.001, 1])
    self.coolingScheme.append([11, 0.01, 0.01,  1, 7500,0.0005, 1])

    self.coolingStep = None

    BasePopup.__init__(self, parent, title="Hydrogen Cloud Molecular Dynamics", **kw)
  
  def body(self, guiFrame):
        
    self.mdInitTempEntry  = FloatEntry(self,text='',returnCallback=self.setMdInitTemp )
    self.mdFinTempEntry   = FloatEntry(self,text='',returnCallback=self.setMdFinTemp  )
    self.mdCoolStepsEntry = IntEntry  (self,text='',returnCallback=self.setMdCoolSteps)
    self.mdSimStepsEntry  = IntEntry  (self,text='',returnCallback=self.setMdSimSteps )
    self.mdTauEntry       = FloatEntry(self,text='',returnCallback=self.setMdTau      )
    self.mdRepScaleEntry  = FloatEntry(self,text='',returnCallback=self.setMdRepScale )
    
    guiFrame.grid_columnconfigure(0, weight=1)

    row = 0
    guiFrame.grid_rowconfigure(row, weight=1)
    frame = LabelFrame(guiFrame, text='Input constraints')
    frame.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame.grid_columnconfigure(2, weight=1)

    srow = 0
    label = Label(frame, text='Constraint set:')
    label.grid(row=srow, column=0, sticky=Tkinter.W)
    self.constraintSetPulldown = PulldownMenu(frame, callback=self.changeConstraintSet, selected_index=0, do_initial_callback=0)
    self.constraintSetPulldown.grid(row=srow, column=1,sticky=Tkinter.W )

    srow += 1
    label = Label(frame, text='Dist constraint list 1:')
    label.grid(row=srow, column=0, sticky=Tkinter.W)
    self.distance1Pulldown = PulldownMenu(frame, callback=self.changeDistance1ConstraintList, selected_index=0, do_initial_callback=0)
    self.distance1Pulldown.grid(row=srow, column=1,sticky=Tkinter.W )
    self.numConstr1Label = Label(frame, text='Constraints: 0')
    self.numConstr1Label.grid(row=srow, column=2, sticky=Tkinter.W)

    srow += 1
    label = Label(frame, text='Dist constraint list 2:')
    label.grid(row=srow, column=0, sticky=Tkinter.W)
    self.distance2Pulldown = PulldownMenu(frame, callback=self.changeDistance2ConstraintList, selected_index=0, do_initial_callback=0)
    self.distance2Pulldown.grid(row=srow, column=1,sticky=Tkinter.W )
    self.numConstr2Label = Label(frame, text='Constraints: 0')
    self.numConstr2Label.grid(row=srow, column=2, sticky=Tkinter.W)

    srow += 1
    label = Label(frame, text='Dist constraint list 3:')
    label.grid(row=srow, column=0, sticky=Tkinter.W)
    self.distance3Pulldown = PulldownMenu(frame, callback=self.changeDistance3ConstraintList, selected_index=0, do_initial_callback=0)
    self.distance3Pulldown.grid(row=srow, column=1,sticky=Tkinter.W )
    self.numConstr3Label = Label(frame, text='Constraints: 0')
    self.numConstr3Label.grid(row=srow, column=2, sticky=Tkinter.W)

    srow += 1
    label = Label(frame, text='Dist constraint list 4:')
    label.grid(row=srow, column=0, sticky=Tkinter.W)
    self.distance4Pulldown = PulldownMenu(frame, callback=self.changeDistance4ConstraintList, selected_index=0, do_initial_callback=0)
    self.distance4Pulldown.grid(row=srow, column=1,sticky=Tkinter.W )
    self.numConstr4Label = Label(frame, text='Constraints: 0')
    self.numConstr4Label.grid(row=srow, column=2, sticky=Tkinter.W)

    row += 1
    frame0 = LabelFrame(guiFrame, text='Cooling scheme')
    frame0.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame0.grid_columnconfigure(1, weight=1)

    f0row = 0
    frame0.grid_rowconfigure(f0row, weight=1)
    colHeadings      = ['Step','Initial\nTemp.','Final\nTemp.','Cooling\nSteps','MD Steps','MD Tau','Rep.\nScale']
    editWidgets      = [None,self.mdInitTempEntry,self.mdFinTempEntry,self.mdCoolStepsEntry,self.mdSimStepsEntry,self.mdTauEntry,self.mdRepScaleEntry]
    editGetCallbacks = [None,self.getMdInitTemp,  self.getMdFinTemp,  self.getMdCoolSteps,  self.getMdSimSteps,  self.getMdTau,  self.getMdRepScale  ]
    editSetCallbacks = [None,self.setMdInitTemp,  self.setMdFinTemp,  self.setMdCoolSteps,  self.setMdSimSteps,  self.setMdTau,  self.setMdRepScale  ]
    self.coolingSchemeMatrix = ScrolledMatrix(frame0, editSetCallbacks=editSetCallbacks, editGetCallbacks=editGetCallbacks, editWidgets=editWidgets,
                                              maxRows=9, initialRows=12, headingList=colHeadings, callback=self.selectCoolingStep, objectList=self.coolingScheme, textMatrix=self.coolingScheme)
    self.coolingSchemeMatrix.grid(row=f0row, column=0, columnspan=4, sticky=Tkinter.NSEW)

    f0row +=1
    texts    = ['Move earlier','Move later','Add step','Remove step']
    commands = [self.moveStepEarlier,self.moveStepLater,self.addCoolingStep,self.removeCoolingStep]
    self.coolingSchemeButtons = ButtonList(frame0,expands=1,commands=commands,texts=texts)
    self.coolingSchemeButtons.grid(row=f0row, column=0, columnspan=4, sticky=Tkinter.EW)

    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    frame1 = LabelFrame(guiFrame, text='Dynamics control')
    frame1.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame1.grid_columnconfigure(1, weight=1)

    f1row =0
    label20 = Label(frame1, text='Number of clouds:')
    label20.grid(row=f1row, column=0, sticky=Tkinter.NW)
    self.numCloudsEntry = IntEntry(frame1, text=1, returnCallback=self.setNumClouds, width=10)
    self.numCloudsEntry.grid(row=f1row, column=1, sticky=Tkinter.NW)
    label21 =Label(frame1, text='Cloud file prefix:')
    label21.grid(row=f1row, column=2, sticky=Tkinter.NW)
    self.filePrefixEntry = Entry(frame1, text='t_intra_', returnCallback=self.setFilePrefix, width=10)
    self.filePrefixEntry.grid(row=f1row, column=3, sticky=Tkinter.NW)

    f1row +=1
    texts    = ['Start molecular dynamics','Show dynamics progress']
    commands = [self.startMd,self.showMdProgress]
    self.mdButtons = ButtonList(frame1, expands=1, commands=commands, texts=texts)
    self.mdButtons.grid(row=f1row, column=0, columnspan=4, sticky=Tkinter.NSEW)

    row += 1
    self.bottomButtons = createDismissHelpButtonList(guiFrame,expands=0,help_url=None)
    self.bottomButtons.grid(row=row, column=0, sticky=Tkinter.EW)
    
    self.update()    

    for func in ('__init__', 'delete', 'setName'):
      for clazz in ('ccp.nmr.NmrConstraint.DistanceConstraintList',):
        Implementation.registerNotify(self.updateConstraintLists,clazz,func)
    for func in ('__init__', 'delete'):
      Implementation.registerNotify(self.updateConstraintSets,'ccp.nmr.NmrConstraint.NmrConstraintStore', func)
    
  def getContraintSetNames(self):
  
    names = []
    constraintSets = self.project.currentNmrProject.nmrConstraintStores
    for set in constraintSets:
      names.append('%d' % set.serial)
      
    return names
            
  def changeConstraintSet(self, i, name):    

    project = self.project
    if project.currentNmrProject.nmrConstraintStores:
      constraintSet = project.currentNmrProject.sortedNmrConstraintStores()[i]
    else:
      constraintSet = None
    
    if constraintSet is not self.constraintSet: 
      self.constraintSet  = constraintSet
    
    self.updateConstraintLists()
    self.update()
 
  def updateConstraintLists(self, *opt):
  
    constrListData  = self.getConstraintLists()
    constrListNames = [x[0] for x in constrListData]
    constraintLists = [x[1] for x in constrListData]
    # copes with self.constraintSet being None
    
    if constrListNames:
      i = 0
      for constraintList in self.constrLists:
        if constraintList not in constraintLists:
          if i == 0:
            self.constrLists[i] = constraintLists[0]
          else:
            self.constrLists[i] = None
        i += 1
      
      constraintLists.append(None)
      constrListNames.append('<None>')
      self.distance1Pulldown.setup(constrListNames,constraintLists.index(self.constrLists[0]))
      self.distance2Pulldown.setup(constrListNames,constraintLists.index(self.constrLists[1]))
      self.distance3Pulldown.setup(constrListNames,constraintLists.index(self.constrLists[2]))
      self.distance4Pulldown.setup(constrListNames,constraintLists.index(self.constrLists[3]))
        
    else:
      self.constrLists = [None] * 4
      self.distance1Pulldown.setup([],-1)
      self.distance2Pulldown.setup([],-1)
      self.distance3Pulldown.setup([],-1)
      self.distance4Pulldown.setup([],-1)
 
  def updateConstraintSets(self, *opt):
  
    project    = self.project
    constraintSets = list(project.currentNmrProject.nmrConstraintStores)
    if constraintSets:
      constraintSetNames = self.getContraintSetNames()
      
      # set defaults
      if self.constraintSet not in constraintSets:
        self.constraintSet = constraintSets[0]
        
      if self.constraintSet:
        j = 0
        for constraintList in self.constrLists:
          if constraintList and (constraintList.nmrConstraintStore is not self.constraintSet):
            if self.constraintSet.constraintLists and j==0:
              self.constrLists[j] = self.constraintSet.constraintLists[0]
            else:
              self.constrLists[j] = None
          j += 1

      else:
        self.constrLists = [None] * 4
          
      i =  constraintSets.index(self.constraintSet)
      self.constraintSetPulldown.setup(constraintSetNames,i)

    else:
      self.constraintSet  = None
      self.constrLists = [None] * 4
      self.constraintSetPulldown.setup([],-1)

  def getConstraintListName(self, constraintList):
    
    if constraintList.name:
      listName = ':%s' % (constraintList.name) 
    else:
      listName = ''
      
    name = '%d:%d:%s%s' % (constraintList.nmrConstraintStore.serial,constraintList.serial,constraintList.className[:-14],listName)
    return name


  def getConstraintLists(self):
  
    constraintLists = []
    if self.constraintSet:
      for constraintList in self.constraintSet.constraintLists:
        if constraintList.className == 'DistanceConstraintList':
          name = self.getConstraintListName(constraintList)
          constraintLists.append( [name, constraintList] )
    
    return constraintLists

    
  def changeDistance1ConstraintList(self, i, name): 
    self.changeDistanceConstraintList(i, name, 0)
    
  def changeDistance2ConstraintList(self, i, name): 
    self.changeDistanceConstraintList(i, name, 1)
    
  def changeDistance3ConstraintList(self, i, name): 
    self.changeDistanceConstraintList(i, name, 2)
    
  def changeDistance4ConstraintList(self, i, name): 
    self.changeDistanceConstraintList(i, name, 3)

  def changeDistanceConstraintList(self, i, name, listNum): 
    
    project = self.project
    constraintLists = self.getConstraintLists()
        
    if constraintLists and (i < 4):
      self.constrLists[listNum] = constraintLists[i][1]
    else:
      self.constrLists[listNum] = None
    
    self.update()

  def startMd(self):
  
    self.setNumClouds()
    self.setFilePrefix()
    if ( (self.constrLists != [None] * 4)
        and (self.numClouds > 0) and self.filePrefix):
      
      resDict = {}
      for resonance in self.guiParent.project.currentNmrProject.resonances:
        resDict[resonance.serial] = resonance
      
      constraints = []
      
      constraintStore = None
      for dcl in self.constrLists:
        if dcl:
          constraintStore = dcl.nmrConstraintStore
          constraints.extend( list(dcl.constraints) )

      resonances = []
      for constraint in constraints:
        for item in constraint.items:
          for fixedResonance in item.resonances:
            if fixedResonance.resonanceSerial is None:
              resonance = newResonance(self.guiParent.project, isotopeCode=fixedResonance.isotopeCode)
              resonance.setName(fixedResonance.name)
              fixedResonance.setResonanceSerial(resonance.serial)
              resDict[resonance.serial] = resonance
              if fixedResonance.resonanceSet:
                atomSets = list(fixedResonance.resonanceSet.atomSets)
                assignAtomsToRes(atomSets, resonance)
            
            if resDict.get(fixedResonance.resonanceSerial) is not None:
              resonances.append( resDict[fixedResonance.resonanceSerial] )
              resDict[fixedResonance.resonanceSerial] = None

      resonances, intraConstraintList = self.makeIntraConstraints(resonances, constraintStore)
      constraints.extend(  list(intraConstraintList.constraints) )
      resonances, interConstraintList = self.makeInterConstraints(resonances, constraintStore)
      constraints.extend(  list(interConstraintList.constraints) )
    
      startMdProcess(self.numClouds, constraints,
                     resonances, self.coolingScheme, self.filePrefix)
   
      #structGen = self.distanceConstraintList.structureGeneration

      serials = []
      for resonance in resonances:
        serials.append(resonance.serial)
      clouds = []
      for i in range(self.numClouds):
        clouds.append( '%s%3.3d.pdb' % (self.filePrefix,i) )
      self.guiParent.application.setValues(constraintStore, 'clouds', values=clouds)
      self.guiParent.application.setValues(constraintStore, 'cloudsResonances', values=serials)

      # do better than this check for creation
  
  def makeInterConstraints(self, resonances, constraintSet):

    from ccpnmr.analysis.core.ConstraintBasic import getFixedResonance
    from ccpnmr.analysis.core.AssignmentBasic import findConnectedSpinSystem
    
    project = constraintSet.root
    constraintList = constraintSet.newDistanceConstraintList()
    constraintList.name = 'Seq connections'

    resDict = {}
    spinSystemDict = {}
    for resonance in resonances:
      resDict[resonance] = True
      spinSystem = resonance.resonanceGroup
      if spinSystem:
        spinSystemDict[spinSystem] = None

    spinSystems = spinSystemDict.keys()
    
    for spinSystem in spinSystems:
      nextSpinSystem = findConnectedSpinSystem(spinSystem, delta=1)
      if nextSpinSystem:
        ca = spinSystem.newAtoms.get('CA')
        c  = spinSystem.newAtoms.get('C')
        n  = nextSpinSystem.newAtoms.get('N')
        
        
        if ca and c and n:
          if resDict.get(ca) is None:
            resonances.append(ca)
          if resDict.get(c) is None:
            resonances.append(c)
          if resDict.get(n) is None:
            resonances.append(n)
            
          c_n = constraintList.newDistanceConstraint(weight=1.0, origData=1.0, targetValue=1.32,
                                                     upperLimit=1.35, lowerLimit=1.29, error=0.06)
 
          # below based on angle constraints
          ca_n = constraintList.newDistanceConstraint(weight=1.0, origData=1.0, targetValue=2.415,
                                                      upperLimit=2.513, lowerLimit=2.316, error=0.197)
                                                      
          frCa = getFixedResonance(constraintSet,ca)
          frC  = getFixedResonance(constraintSet,c)
          frN  = getFixedResonance(constraintSet,n)
 
          item = ca_n.newDistanceConstraintItem(resonances=[frCa,frN])
          item = c_n.newDistanceConstraintItem(resonances=[frC,frN])

    return resonances, constraintList

  def makeIntraConstraints(self, resonances, constraintSet):
    
    from ccpnmr.analysis.core.ConstraintBasic import getFixedResonance
    
    project = constraintSet.root
    constraintList = constraintSet.newDistanceConstraintList()
    constraintList.name = 'Backbone intra'
    
    dict = {}
    for resonance in resonances:
      if resonance.resonanceGroup and resonance.assignNames:
        dict[resonance] = 1
        resonance.resonanceGroup.newAtoms = {}
    
    for resonance in dict.keys():
      ss = resonance.resonanceGroup
      if resonance.assignNames[0] == 'H':
        for resonance2 in ss.resonances:
          if dict.get(resonance2) and resonance2.assignNames[0][:2] == 'HA':
          
            ca = ss.newAtoms.get('CA')
            if ca is None:
              ca =  project.newResonance(isotopeCode='13C', assignNames= ['CA',])
              resonances.append(ca)

            c = ss.newAtoms.get('C')
            if c is None:
              c =  project.newResonance(isotopeCode='13C', assignNames= ['C',])
              resonances.append(c)
             
            n = ss.newAtoms.get('N')
            if n is None:
              n = project.newResonance(isotopeCode='15N', assignNames= ['N',])
              resonances.append(n)  
            
            ss.newAtoms['C']  = c
            ss.newAtoms['CA'] = ca
            ss.newAtoms['N']  = n
 
            frCa = getFixedResonance(constraintSet,ca)
            frC  = getFixedResonance(constraintSet,c)
            frN  = getFixedResonance(constraintSet,n)
            frH  = getFixedResonance(constraintSet,resonance)
            frha  = getFixedResonance(constraintSet,resonance2)
            
            h_n  = constraintList.newDistanceConstraint(weight=1.0, origData=1.0, targetValue=1.01,
                                                        upperLimit=1.1, lowerLimit=0.9,  error=0.2)
            n_ca = constraintList.newDistanceConstraint(weight=1.0, origData=1.0, targetValue=1.465,
                                                        upperLimit=1.50, lowerLimit=1.43, error=0.07)
            ca_c = constraintList.newDistanceConstraint(weight=1.0, origData=1.0, targetValue=1.525,
                                                        upperLimit=1.55, lowerLimit=1.50, error=0.05)
                                                        
            # below based on angle constraints                                            
            n_c = constraintList.newDistanceConstraint(weight=1.0, origData=1.0, targetValue=2.464,
                                                        upperLimit=2.572, lowerLimit=2.356, error=0.216)
                                                       
            item = h_n.newDistanceConstraintItem(resonances=[frH,frN])
            item = n_ca.newDistanceConstraintItem(resonances=[frN,frCa])
            item = ca_c.newDistanceConstraintItem(resonances=[frCa,frC])
            item = n_c.newDistanceConstraintItem(resonances=[frN,frC])
            
            if ss.newAtoms.get('CAHA') is None:
              ca_ha = constraintList.newDistanceConstraint(weight=1.0, origData=1.0, targetValue=1.09,
                                                         upperLimit=1.19, lowerLimit=0.99, error=0.2)
              item = ca_ha.newDistanceConstraintItem(resonances=[frC,frha])
              ss.newAtoms['CAHA']   = 1
      
      if resonance.assignNames[0][:2] == 'HA':
        for resonance2 in ss.resonances:
          if dict.get(resonance2) and resonance2.assignNames[0][:2] == 'HB':
            
            ca = ss.newAtoms.get('CA')
            if ca is None:
              ca =  project.newResonance(isotopeCode='13C', assignNames= ['CA',])
              resonances.append(ca)

            cb = ss.newAtoms.get('CB')
            if cb is None:
              cb =  project.newResonance(isotopeCode='13C', assignNames= ['CB',])
              resonances.append(cb)
 
            ss.newAtoms['CA']   = cb
            ss.newAtoms['CB']   = ca
            ss.newAtoms['CAHA'] = 1
 
            frCA = getFixedResonance(constraintSet,ca)
            frCB = getFixedResonance(constraintSet,cb)
            frHA = getFixedResonance(constraintSet,resonance)
            frHB = getFixedResonance(constraintSet,resonance2)
            
            c_b = constraintList.newDistanceConstraint(weight=1.0, origData=1.0, targetValue=1.09,
                                                       upperLimit=1.19, lowerLimit=0.99, error=0.2)
            c_c = constraintList.newDistanceConstraint(weight=1.0, origData=1.0, targetValue=1.53,
                                                       upperLimit=1.56, lowerLimit=1.51, error=0.05)
                                                       
            item = c_b.newDistanceConstraintItem(resonances=[frCB,frHB])
            item = c_c.newDistanceConstraintItem(resonances=[frCA,frCB])

            if ss.newAtoms.get('CAHA') is None:
              c_a = constraintList.newDistanceConstraint(weight=1.0, origData=1.0, targetValue=1.09,
                                                         upperLimit=1.19, lowerLimit=0.99, error=0.2)
              item = c_a.newDistanceConstraintItem(resonances=[frCA,frHA])
              ss.newAtoms['CAHA']   = 1


    return resonances, constraintList
   
  def showMdProgress(self):

    n = 0
    m = self.numClouds
    for i in range(m):
      pdbFileName = '%s%3.3d.pdb' % (self.filePrefix,i)
      if os.path.exists(pdbFileName):
        n += 1
      
    p = n*100/float(m)
    text = 'Done %d of %d clouds (%1.2f)%%' % (n,m,p)
    showInfo('MD Progress',text)

  def setFilePrefix(self, text=None):
  
    if not text:
      text = self.filePrefixEntry.get()
    
    if text:
      self.filePrefix = text
  
  def setNumClouds(self, n=None, *event):

   if not n:
     n = self.numCloudsEntry.get()
   
   if n:
     self.numClouds = int(n)

  def update(self):
    
    self.updateConstraintSets()    
    self.updateConstraintLists()    
 
    if (self.constrLists != [None] * 4) and self.coolingScheme:
      self.mdButtons.buttons[0].enable()
      self.mdButtons.buttons[1].enable()
    else:
      self.mdButtons.buttons[0].disable()
      self.mdButtons.buttons[1].disable()

    if self.constrLists[0]:
      self.numConstr1Label.set('Constraints: ' + str(len(self.constrLists[0].constraints)))
    else:
      self.numConstr1Label.set('Constraints: 0')

    if self.constrLists[1]:
      self.numConstr2Label.set('Constraints: ' + str(len(self.constrLists[1].constraints)))
    else:
      self.numConstr2Label.set('Constraints: 0')

    if self.constrLists[2]:
      self.numConstr3Label.set('Constraints: ' + str(len(self.constrLists[2].constraints)))
    else:
      self.numConstr3Label.set('Constraints: 0')

    if self.constrLists[3]:
      self.numConstr4Label.set('Constraints: ' + str(len(self.constrLists[3].constraints)))
    else:
      self.numConstr4Label.set('Constraints: 0')


  def getMdInitTemp(self, coolingStep):

    self.mdInitTempEntry.set(coolingStep[1])

  def getMdFinTemp(self, coolingStep): 

    self.mdFinTempEntry.set(coolingStep[2])

  def getMdCoolSteps(self, coolingStep):

    self.mdCoolStepsEntry.set(coolingStep[3])

  def getMdSimSteps(self, coolingStep):

    self.mdSimStepsEntry.set(coolingStep[4])

  def getMdTau(self, coolingStep): 

    self.mdTauEntry.set(coolingStep[5])

  def getMdRepScale(self, coolingStep):

    self.mdRepScaleEntry.set(coolingStep[6])  

  def setMdInitTemp(self, event):

    value = self.mdInitTempEntry.get()
    if value is not None:
      self.coolingStep[1] = value
      
    self.updateCoolingScheme()

  def setMdFinTemp(self, event):  

    value = self.mdFinTempEntry.get()
    if value is not None:
      self.coolingStep[2] = value
      
    self.updateCoolingScheme()

  def setMdCoolSteps(self, event):

    value = self.mdCoolStepsEntry.get()
    if value is not None:
      self.coolingStep[3] = value
      
    self.updateCoolingScheme()

  def setMdSimSteps(self, event): 

    value = self.mdSimStepsEntry.get()
    if value is not None:
      self.coolingStep[4] = value
      
    self.updateCoolingScheme()

  def setMdTau(self, event):      

    value = self.mdTauEntry.get()
    if value is not None:
      self.coolingStep[5] = value
      
    self.updateCoolingScheme()

  def setMdRepScale(self, event): 

    value = self.mdRepScaleEntry.get()
    if value is not None:
      self.coolingStep[6] = value
      
    self.updateCoolingScheme()

  def selectCoolingStep(self, object, row, col):
  
    self.coolingStep = object

  def moveStepEarlier(self):
  
    if self.coolingStep:
      i = self.coolingStep[0] - 1
      if i > 0:
        coolingStep = self.coolingScheme[i-1]
        coolingStep[0] = i+1
        self.coolingStep[0] = i
        self.coolingScheme[i-1] = self.coolingStep
        self.coolingScheme[i]   = coolingStep

        self.updateCoolingScheme()
        self.coolingSchemeMatrix.hilightObject(self.coolingStep)
  
  def moveStepLater(self):

    if self.coolingStep:
      i = self.coolingStep[0] -1
      if i < len(self.coolingScheme) - 1:
        coolingStep = self.coolingScheme[i+1]
        coolingStep[0] = i+1
        self.coolingStep[0] = i+2
        self.coolingScheme[i+1] = self.coolingStep
        self.coolingScheme[i]   = coolingStep

        self.updateCoolingScheme()
        self.coolingSchemeMatrix.hilightObject(self.coolingStep)

  def addCoolingStep(self):
  
    i = len(self.coolingScheme) + 1
    datum = [i, 3000,100, 10, 2500, 0.001, 1]
  
    self.coolingScheme.append(datum)
    self.updateCoolingScheme()
  
  def removeCoolingStep(self):

    if self.coolingStep:
      coolingScheme = []
      i = 0
      for coolingStep in self.coolingScheme:
        if coolingStep is not self.coolingStep:
          i  += 1
          coolingStep[0] = i
          coolingScheme.append( coolingStep )
      
      self.coolingScheme = coolingScheme
      self.updateCoolingScheme()

  def updateCoolingScheme(self):
    
    objectList = self.coolingScheme
    textMatrix = self.coolingScheme
    self.coolingSchemeMatrix.update(objectList=objectList,textMatrix=textMatrix)

  def updateMidgeParams(self):
    
    data = [self.specFreq,self.maxIter,self.mixTime,self.corrTime,self.leakRate,self.maxIntens]

    self.midgeParamsMatrix.update(textMatrix=[data,])

  def destroy(self):

    for func in ('__init__', 'delete', 'setName'):
      for clazz in ('ccp.nmr.NmrConstraint.DistanceConstraintList',):
        Implementation.unregisterNotify(self.updateConstraintLists,clazz,func)
    for func in ('__init__', 'delete'):
      Implementation.unregisterNotify(self.updateConstraintSets,'ccp.nmr.NmrConstraint.NmrConstraintStore', func)

    BasePopup.destroy(self)

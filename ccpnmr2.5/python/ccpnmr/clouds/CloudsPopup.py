"""
======================COPYRIGHT/LICENSE START==========================

CloudsPopup.py: Part of the CcpNmr Clouds program

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

from memops.api import Implementation
from ccp.api.nmr import Nmr, NmrConstraint

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

from ccpnmr.analysis.core.ConstraintBasic import makeNmrConstraintStore
from ccpnmr.analysis.core.ExperimentBasic import getSpectraByType
from ccpnmr.analysis.core.StructureBasic  import getAtomSetsDistance

from ccpnmr.clouds.ResonanceIdentification import getCloudsResonanceList, makeNoeAdcs, constrainSpinSystems
from ccpnmr.clouds.NoeRelaxation           import optimiseRelaxation, disambiguateNoesyPeaks
from ccpnmr.clouds.HydrogenDynamics        import *

class CloudsPopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.guiParent = parent
    self.project   = parent.getProject()
    self.waiting   = 0
    self.specFreq  = 800.13
    self.maxIter   = 50
    self.mixTime   = 60
    self.corrTime  = 11.5
    self.leakRate  = 2.0
    self.peakListDict  = {}
    self.noesyPeakList = None
    self.tocsyPeakList = None  
    self.noesy3dPeakList = None
    self.hsqcPeakList  = None
    self.maxIntens     = 37000000 
    
    self.resonances = None
    self.origResonances = None
    self.noesyPeaks = None
    self.distanceConstraintList = None
    self.antiDistConstraintList = None
    self.numClouds    = 100
    self.filePrefix   = 'cloud_'
    self.cloudsFiles  = []
    self.adcAtomTypes = 'HN'
    self.structure    = None
    
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

    BasePopup.__init__(self, parent, title="Resonance Clouds Analysis", **kw)
  
  def body(self, guiFrame):

    self.specFreqEntry = IntEntry  (self,text=self.specFreq, width=8,returnCallback = self.setSpecFreq)
    self.maxIterEntry  = IntEntry  (self,text=self.maxIter,  width=8,returnCallback = self.setMaxIter )
    self.mixTimeEntry  = FloatEntry(self,text=self.mixTime,  width=8,returnCallback = self.setMixTime )
    self.corrTimeEntry = FloatEntry(self,text=self.corrTime, width=8,returnCallback = self.setCorrTime)
    self.leakRateEntry = FloatEntry(self,text=self.leakRate, width=8,returnCallback = self.setLeakRate)
    self.maxIntensEntry = IntEntry(self,text=self.maxIntens, width=8,returnCallback = self.setMaxIntens)
        
    self.mdInitTempEntry  = FloatEntry(self,text='',returnCallback=self.setMdInitTemp )
    self.mdFinTempEntry   = FloatEntry(self,text='',returnCallback=self.setMdFinTemp  )
    self.mdCoolStepsEntry = IntEntry  (self,text='',returnCallback=self.setMdCoolSteps)
    self.mdSimStepsEntry  = IntEntry  (self,text='',returnCallback=self.setMdSimSteps )
    self.mdTauEntry       = FloatEntry(self,text='',returnCallback=self.setMdTau      )
    self.mdRepScaleEntry  = FloatEntry(self,text='',returnCallback=self.setMdRepScale )
    
    guiFrame.grid_columnconfigure(0, weight=1)

    row = 0
    frame0 = LabelFrame(guiFrame, text='Setup peak lists')
    frame0.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame0.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame0.grid_columnconfigure(1, weight=1)
    
    f0row = 0
    label00 = Label(frame0, text='1H-1H NOESY spectrum')
    label00.grid(row=f0row, column=0, sticky=Tkinter.NW)
    self.noesyPulldown = PulldownMenu(frame0, entries=self.getNoesys(), callback=self.setNoesy, selected_index=0, do_initial_callback=0)
    self.noesyPulldown.grid(row=f0row, column=1,sticky=Tkinter.NW )
  
    f0row += 1
    label01 = Label(frame0, text='15N HSQC spectrum')
    label01.grid(row=f0row, column=0, sticky=Tkinter.NW)
    self.hsqcPulldown = PulldownMenu(frame0, entries=self.getHsqcs(), callback=self.setHsqc, selected_index=0, do_initial_callback=0)
    self.hsqcPulldown.grid(row=f0row, column=1,sticky=Tkinter.NW )

    f0row += 1
    label02 = Label(frame0, text='15N HSQC TOCSY spectrum')
    label02.grid(row=f0row, column=0, sticky=Tkinter.NW)
    self.tocsyPulldown = PulldownMenu(frame0, entries=self.getTocsys(), callback=self.setTocsy, selected_index=0, do_initial_callback=0)
    self.tocsyPulldown.grid(row=f0row, column=1,sticky=Tkinter.NW )

    f0row += 1
    label02 = Label(frame0, text='15N HSQC NOESY spectrum')
    label02.grid(row=f0row, column=0, sticky=Tkinter.NW)
    self.noesy3dPulldown = PulldownMenu(frame0, entries=self.getNoesy3ds(), callback=self.setNoesy3d, selected_index=0, do_initial_callback=0)
    self.noesy3dPulldown.grid(row=f0row, column=1,sticky=Tkinter.NW )

    f0row += 1
    texts    = ['Setup resonances & peaks','Show Peaks','Show resonances']
    commands = [self.setupResonances,self.showPeaks,self.showResonances]
    self.setupButtons = ButtonList(frame0,expands=1,texts=texts,commands=commands)
    self.setupButtons.grid(row=f0row, column=0, columnspan=2, sticky=Tkinter.NSEW)

    f0row += 1
    self.label03a = Label(frame0, text='Resonances found: 0')
    self.label03a.grid(row=f0row, column=0, sticky=Tkinter.NW)
    self.label03b = Label(frame0, text='NOESY peaks found: 0')
    self.label03b.grid(row=f0row, column=1, sticky=Tkinter.NW)
  
    row += 1
    frame1 = LabelFrame(guiFrame, text='Calculate distance constraints')
    frame1.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame1.grid_columnconfigure(3, weight=1)
    
    f1row = 0
    frame1.grid_rowconfigure(f1row, weight=1)
    data = [self.specFreq,self.maxIter,self.mixTime,self.corrTime,self.leakRate,self.maxIntens]
    colHeadings      = ['Spectrometer\nfrequency','Max\niterations','Mixing\ntime (ms)','Correl.\ntime (ns)','Leak\nrate','Max\nintensity']
    editWidgets      = [self.specFreqEntry,self.maxIterEntry,self.mixTimeEntry,self.corrTimeEntry,self.leakRateEntry,self.maxIntensEntry,]
    editGetCallbacks = [self.getSpecFreq,  self.getMaxIter,  self.getMixTime,  self.getCorrTime,  self.getLeakRate,  self.getMaxIntens,  ]
    editSetCallbacks = [self.setSpecFreq,  self.setMaxIter,  self.setMixTime,  self.setCorrTime,  self.setLeakRate,  self.setMaxIntens,  ]
    self.midgeParamsMatrix = ScrolledMatrix(frame1, editSetCallbacks=editSetCallbacks, editGetCallbacks=editGetCallbacks, editWidgets=editWidgets,
                                            maxRows=1, initialCols=5, headingList=colHeadings, callback=None, objectList=['None',], textMatrix=[data,])
    self.midgeParamsMatrix.grid(row=f1row, column=0, columnspan=4, sticky=Tkinter.NSEW)

    f1row += 1
    label10 = Label(frame1, text='Benchmark structure')
    label10.grid(row=f1row, column=0, sticky=Tkinter.NW)
    self.structurePulldown = PulldownMenu(frame1, entries=self.getStructures(), callback=self.setStructure, selected_index=0, do_initial_callback=0)
    self.structurePulldown.grid(row=f1row, column=1,sticky=Tkinter.NW )

    label11 = Label(frame1, text='ADC atom types:')
    label11.grid(row=f1row, column=2, sticky=Tkinter.NW)
    self.adcAtomsPulldown = PulldownMenu(frame1, entries=self.getAdcAtomTypes(), callback=self.setAdcAtomTypes, selected_index=0, do_initial_callback=0)
    self.adcAtomsPulldown.grid(row=f1row, column=3,sticky=Tkinter.NW )

    f1row += 1
    texts    = ['Calculate distances','Show distance\nconstraints','Show anti-distance\nconstraints']
    commands = [self.calculateDistances,self.showConstraints,self.showAntiConstraints]
    self.midgeButtons = ButtonList(frame1,expands=1,texts=texts,commands=commands)
    self.midgeButtons.grid(row=f1row, column=0, columnspan=4, sticky=Tkinter.NSEW)

    f1row += 1
    self.distConstrLabel = Label(frame1, text='Distance constraints:')
    self.distConstrLabel.grid(row=f1row, column=0, columnspan=2, sticky=Tkinter.NW)
    self.antiConstrLabel = Label(frame1, text='Anti-distance constraints:')
    self.antiConstrLabel.grid(row=f1row, column=2, columnspan=2, sticky=Tkinter.NW)
  
    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    frame2 = LabelFrame(guiFrame, text='Proton cloud molecular dynamics')
    frame2.grid(row=row, column=0, sticky=Tkinter.NSEW)
    frame2.grid_columnconfigure(1, weight=1)

    f2row = 0
    frame2.grid_rowconfigure(f2row, weight=1)
    data = [self.specFreq,self.maxIter,self.mixTime,self.corrTime,self.leakRate]
    colHeadings      = ['Step','Initial temp.','Final temp.','Cooling steps','MD steps','MD tau','Rep. scale']
    editWidgets      = [None,self.mdInitTempEntry,self.mdFinTempEntry,self.mdCoolStepsEntry,self.mdSimStepsEntry,self.mdTauEntry,self.mdRepScaleEntry]
    editGetCallbacks = [None,self.getMdInitTemp,  self.getMdFinTemp,  self.getMdCoolSteps,  self.getMdSimSteps,  self.getMdTau,  self.getMdRepScale  ]
    editSetCallbacks = [None,self.setMdInitTemp,  self.setMdFinTemp,  self.setMdCoolSteps,  self.setMdSimSteps,  self.setMdTau,  self.setMdRepScale  ]
    self.coolingSchemeMatrix = ScrolledMatrix(frame2, editSetCallbacks=editSetCallbacks, editGetCallbacks=editGetCallbacks, editWidgets=editWidgets,
                                              maxRows=9, initialRows=12, headingList=colHeadings, callback=self.selectCoolingStep, objectList=self.coolingScheme, textMatrix=self.coolingScheme)
    self.coolingSchemeMatrix.grid(row=f2row, column=0, columnspan=4, sticky=Tkinter.NSEW)

    f2row +=1
    texts    = ['Move earlier','Move later','Add step','Remove step']
    commands = [self.moveStepEarlier,self.moveStepLater,self.addCoolingStep,self.removeCoolingStep]
    self.coolingSchemeButtons = ButtonList(frame2,expands=1,commands=commands,texts=texts)
    self.coolingSchemeButtons.grid(row=f2row, column=0, columnspan=4, sticky=Tkinter.EW)

    f2row +=1
    label20 = Label(frame2, text='Number of clouds:')
    label20.grid(row=f2row, column=0, sticky=Tkinter.NW)
    self.numCloudsEntry = FloatEntry(frame2, text=100, returnCallback=self.setNumClouds, width=10)
    self.numCloudsEntry.grid(row=f2row, column=1, sticky=Tkinter.NW)
    label21 =Label(frame2, text='Cloud file prefix:')
    label21.grid(row=f2row, column=2, sticky=Tkinter.NW)
    self.filePrefixEntry = Entry(frame2, text='cloud_', returnCallback=self.setFilePrefix, width=10)
    self.filePrefixEntry.grid(row=f2row, column=3, sticky=Tkinter.NW)

    f2row +=1
    texts    = ['Start molecular dynamics','Show dynamics progress']
    commands = [self.startMd,self.showMdProgress]
    self.mdButtons = ButtonList(frame2, expands=1, commands=commands, texts=texts)
    self.mdButtons.grid(row=f2row, column=0, columnspan=4, sticky=Tkinter.NSEW)

    row += 1
    self.bottomButtons = createDismissHelpButtonList(guiFrame,expands=0,help_url=None)
    self.bottomButtons.grid(row=row, column=0, sticky=Tkinter.EW)
    
    self.setButtonStates()    
    
  def getStructures(self):
  
    names = ['<None>',]
    for molSystem in self.project.sortedMolSystems():
      for structure in molSystem.sortedStructureEnsembles():
        names.append('%s:%d' % (molSystem.name,structure.ensembleId) )
        
    return names
  
  def setStructure(self, index, name=None):
  
    if index < 1:
      self.structure = None
    else:
      structures = []
      for molSystem in self.project.molSystems:
        for structure in molSystem.structureEnsembles:
          structures.append( structure )
 
      self.structure = structures[index-1]
    
  def getAdcAtomTypes(self):
    
    return ['HN','HN HA','HN HA HB']
  
  def setAdcAtomTypes(self, index, name=None):
  
    if name is None:
      name = self.adcAtomsPulldown.getSelected()
  
    self.adcAtomTypes = name
  
  def startMd(self):
  
    self.setNumClouds()
    self.setFilePrefix()
    if (self.distanceConstraintList
        and self.antiDistConstraintList
        and (self.numClouds > 0) and self.filePrefix):
    
      
      resDict = {}
      for resonance in self.guiParent.project.currentNmrProject.resonances:
        resDict[resonance.serial] = resonance

      resonances = []
      for constraint in self.distanceConstraintList.constraints:
        for item in constraint.items:
          for fixedResonance in item.resonances:
            if resDict.get(fixedResonance.resonanceSerial) is not None:
              resonances.append( resDict[fixedResonance.resonanceSerial] )
              resDict[fixedResonance.resonanceSerial] = None
    
      startMdProcess(self.numClouds, self.distanceConstraintList,
                     resonances, self.coolingScheme, self.filePrefix)
   
      #structGen = self.distanceConstraintList.structureGeneration

      serials = []
      for resonance in resonances:
        serials.append(resonance.serial)
      clouds = []
      for i in range(self.numClouds):
        clouds.append( '%s%3.3d.pdb' % (self.filePrefix,i) )
      self.guiParent.application.setValues(self.distanceConstraintList.nmrConstraintStore, 'clouds', values=clouds)
      self.guiParent.application.setValues(self.distanceConstraintList.nmrConstraintStore, 'cloudsResonances', values=serials)

      # do better than this check for creation
   
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
   
  def calculateDistances(self):

    # setup normalisation factor intensityMax

    # self.maxIter
    # what if failure ?
    
    resDict = {}
    for resonance in self.project.currentNmrProject.resonances:
      resDict[resonance.serial] = resonance
    
    self.resonances  = self.origResonances
    intensityFactors = [1.0 for x in range(len(self.resonances))]
    
    
    # optimiseRelaxation will remove unconstrained resonances
    self.distanceConstraintList = optimiseRelaxation(self.resonances,self.noesyPeaks,intensityMax=self.maxIntens,intensityFactors=intensityFactors,
                                                     tmix=self.mixTime,sf=self.specFreq,tcor=self.corrTime,rleak=self.leakRate)

    constrainSpinSystems(self.distanceConstraintList)
    # for testing calculate distances from structure overrides any resonances: uses assigned ones
    #(self.distanceConstraintList, self.resonances) = self.cheatForTesting()
    #self.antiDistConstraintList = self.distanceConstraintList
    protonNumbs = { 'CH3': 3, 'Haro': 2, 'HN': 1, 'H': 1}
    
    PI = 3.1415926535897931
    GH = 2.6752e4
    HB = 1.05459e-27
    CONST = GH*GH*GH*GH * HB*HB
    tc = 1.0e-9 * self.corrTime
    wh = 2.0 * PI * self.specFreq * 1.0e6
    j0 = CONST * tc
    j1 = CONST * tc / (1.0 + wh*wh*tc*tc)
    j2 = CONST * tc / (1.0 + 4.0*wh*wh*tc*tc)
    #jself = 6.0*j2 + 3.0*j1 + j0
    jcross = 6.0*j2 - j0
    
    
    if self.distanceConstraintList:
      constraintStore = self.distanceConstraintList.nmrConstraintStore
      
      
      dict = {'HN':['H'],
              'HN HA':['H','HA','HA1','HA2'],
              'HN HA HB':['H','HA','HA1','HA2','HB','HB2','HB3']}
      
      self.antiDistConstraintList = makeNoeAdcs(self.resonances,self.noesyPeakList.dataSource,constraintStore,allowedAtomTypes=dict[self.adcAtomTypes])
    
      if self.structure:
      
        N = len(self.resonances)
        sigmas = [[] for i in range(N) ]
        for i in range(N):
          sigmas[i] = [0.0 for j in range(N)]
          
        for constraint in self.distanceConstraintList.constraints:
          resonances = list(constraint,findFirstItem().resonances)
          
          ri = resDict[resonances[0].resonanceSerial]
          rj = resDict[resonances[1].resonanceSerial]
          i = self.resonances.index(ri)
          j = self.resonances.index(rj)
          atomSets1 = list(ri.resonanceSet.atomSets)
          atomSets2 = list(rj.resonanceSet.atomSets)
          if atomSets1 == atomSets2:
            ass = list(atomSets1)
            atomSets1 = [ass[0],]
            atomSets2 = [ass[-1],]
            
          distance = getAtomSetsDistance(atomSets1, atomSets2, self.structure)
          r     = distance * 1e-8
          nhs   = protonNumbs[rj.name]
          sigma = 0.1*jcross*nhs/(r**6)
          sigmas[i][j] = sigma
          
          constraint.setDetails('Known Dist: %4.3f' % (distance))
          #for constraint in self.antiDistConstraintList.constraints:
          #  atomSets1 = list(resonances[0].resonanceSet.atomSets)
          #  atomSets2 = list(resonances[1].resonanceSet.atomSets)
          #  distance = getAtomSetsDistance(atomSets1, atomSets2, self.structure)
          #  constraint.setDetails('Known Dist: %4.3f' % (distance))
       
        fp = open('sigmas.out', 'w')
        for i in range(N-1):
          for j in range(i+1,N):
            if sigmas[i][j] != 0.0:
              fp.write('%3.1d  %3.1d   %9.2e\n' % (i,j,sigmas[i][j]))
          #fp.write('\n')
        fp.close()

    self.setButtonStates()
  
  def cheatForTesting(self, atomSelection='H'):
    """ Makes a perfect cloud from a structure. """
  
    project   = self.project
    structure = self.guiParent.argumentServer.getStructure()
 
    constraintStore = makeNmrConstraintStore(project)
    distConstraintList  = NmrConstraint.DistanceConstraintList(constraintStore)
 
    chain = structure.findFirstCoodChain()
    structureGeneration.hydrogenResonances = []

    molSystem     = structure.molSystem
    atomSets      = []
    resonances    = []
    i =0
    for resonance in project.currentNmrProject.resonances:
 
      if resonance.isotopeCode == '1H':
 
        if resonance.resonanceSet:
 
          atomSet = resonance.resonanceSet.findFirstAtomSet()
          atom = atomSet.findFirstAtom()
          seqId = atom.residue.seqId
          if (seqId < 9) or (seqId > 78):
            continue
 
          if atom.residue.chain.molSystem is molSystem:

            if atomSelection == 'methyl':
              if len(atomSet.atoms) == 3:
                if atom.residue.ccpCode not in ('Ala','Val','Ile','Leu','Thr','Met'):
                   continue
              elif atom.name != 'H':
                 continue
 
            elif atomSelection == 'amide':
              if atom.name != 'H':
                continue
 
            if atom.name == 'H':
              resonance.name = 'HN'
            else:
              resonance.name = 'H'
 
            resonances.append(resonance)
            atomSets.append(list(resonance.resonanceSet.atomSets))
            i += 1
 
    print "Found %d atomSets" % (len(atomSets))
    weight    = 1
    adcWeight = 1
    constrDict = {}
    N = len(atomSets)
    for i in range(N-1):
      atomSets0 = atomSets[i]
      residue0 = atomSets0[0].findFirstAtom().residue.seqId
      print "R", residue0
 
      for j in range(i+1,N):
        if j == i:
          continue
        atomSets1 = atomSets[j]

        dist = getAtomSetsDistance(atomSets0,atomSets1,structure)
        if not dist:
          continue
 
        if dist < 5.5:
          fixedResonance0 = getFixedResonance(constraintStore,resonances[i])
          fixedResonance1 = getFixedResonance(constraintStore,resonances[j])
          constrDict[i] = 1
          constrDict[j] = 1
          constraint = NmrConstraint.DistanceConstraint(distConstraintList, weight=weight, targetValue=dist, upperLimit=dist+(dist/10), lowerLimit=dist-(dist/10), error=dist/5)
          item = NmrConstraint.DistanceConstraintItem(constraint, resonances=[fixedResonance0,fixedResonance1])
 
        elif (atomSets1[0].findFirstAtom().name == 'H') and (atomSets0[0].findFirstAtom().name == 'H') and (dist > 7):
          #else:
          fixedResonance0 = getFixedResonance(constraintStore,resonances[i])
          fixedResonance1 = getFixedResonance(constraintStore,resonances[j])
          constrDict[i] = 1
          constrDict[j] = 1
          constraint = NmrConstraint.DistanceConstraint(distConstraintList, weight=adcWeight, targetValue=75, upperLimit=175, lowerLimit=5.0, error=94.5)
          item = NmrConstraint.DistanceConstraintItem(constraint, resonances=[fixedResonance0,fixedResonance1])
     
    return (distConstraintList, resonances)
  
  def showConstraints(self):
  
    if self.distanceConstraintList:
      self.guiParent.browseConstraints(constraintList=self.distanceConstraintList)
  
  def showAntiConstraints(self):
  
    if self.antiDistConstraintList:
      self.guiParent.browseConstraints(constraintList=self.antiDistConstraintList)

  def showPeaks(self):
  
    self.guiParent.viewPeaks(peaks=self.noesyPeaks)
  
  def showResonances(self):
    
    pass
    #self.guiParent.viewResonances(resonances=self.resonances)

  def setupResonances(self):

    if self.noesyPeakList and self.noesy3dPeakList and self.tocsyPeakList and self.hsqcPeakList:
    
      disambiguateNoesyPeaks(self.noesyPeakList,self.noesy3dPeakList,self.tocsyPeakList,self.hsqcPeakList)
    
      (self.origResonances,self.noesyPeaks,null) = getCloudsResonanceList(self.guiParent.argumentServer,
                                                                      hsqcPeakList=self.hsqcPeakList,
                                                                      tocsy3dPeakList=self.tocsyPeakList,
                                                                      noesy2dPeakList=self.noesyPeakList )
      self.setButtonStates()

  def setButtonStates(self):
     
    if self.origResonances:
      self.label03a.set('Resonances found: %d' % (len(self.origResonances)))
      
    if self.noesyPeaks:
      self.label03b.set('NOESY peaks found: %d' % (len(self.noesyPeaks)))
    
    if self.noesyPeakList and self.tocsyPeakList and self.hsqcPeakList:
      self.setupButtons.buttons[0].enable()
    else:
      self.setupButtons.buttons[0].disable()
    
    if self.noesyPeaks:
      self.setupButtons.buttons[1].enable()
    else:
      self.setupButtons.buttons[1].disable()

    if self.origResonances:
      self.setupButtons.buttons[2].enable()
    else:
      self.setupButtons.buttons[2].disable()
      
    if self.noesyPeaks and self.origResonances:
      self.midgeButtons.buttons[0].enable()
    else:    
      self.midgeButtons.buttons[0].disable()

    if self.distanceConstraintList:
      self.midgeButtons.buttons[1].enable()
      self.distConstrLabel.set('Distance constraints: %d' % len(self.distanceConstraintList.constraints))
    else:    
      self.distConstrLabel.set('Distance constraints:')
      self.midgeButtons.buttons[1].disable()
    
    if self.antiDistConstraintList:
      self.antiConstrLabel.set('Anti-distance constraints: %d' % len(self.antiDistConstraintList.constraints))
      self.midgeButtons.buttons[2].enable()
    else:    
      self.antiConstrLabel.set('Anti-distance constraints:')
      self.midgeButtons.buttons[2].disable()

    if (self.distanceConstraintList
       and self.antiDistConstraintList
       and (self.numClouds > 0)
       and self.filePrefix):
      self.mdButtons.buttons[0].enable()
      self.mdButtons.buttons[1].enable()
    else:
      self.mdButtons.buttons[0].disable()
      self.mdButtons.buttons[1].disable()

  def getNoesys(self):
  
    names = []
    spectra = getSpectraByType(self.project, '2dNOESY')
    for spectrum in spectra:
      for peakList in spectrum.peakLists:
        name = '%s:%s:%s' % (spectrum.experiment.name,spectrum.name,peakList.serial)
        names.append( name )
        self.peakListDict[name] = peakList
        if not self.noesyPeakList:
          self.noesyPeakList = peakList
    
    return names
  
  def setNoesy(self, index, name=None):

    if not name:
      name = self.noesyPulldown.getSelected()
      
    self.noesyPeakList = self.peakListDict[name]
    self.setButtonStates()

  def getTocsys(self):

    names = []
    spectra = getSpectraByType(self.project, '3dTOCSY')
    for spectrum in spectra:
      for peakList in spectrum.peakLists:
        name = '%s:%s:%s' % (spectrum.experiment.name,spectrum.name,peakList.serial)
        names.append( name )
        self.peakListDict[name] = peakList
        if not self.tocsyPeakList:
          self.tocsyPeakList = peakList
    
    return names
  
  def getNoesy3ds(self):

    names = []
    spectra = getSpectraByType(self.project, '3dNOESY')
    for spectrum in spectra:
      for peakList in spectrum.peakLists:
        name = '%s:%s:%s' % (spectrum.experiment.name,spectrum.name,peakList.serial)
        names.append( name )
        self.peakListDict[name] = peakList
        if not self.noesy3dPeakList:
          self.noesy3dPeakList = peakList
    
    return names
  
  def setTocsy(self, index, name=None):

    if not name:
      name = self.tocsyPulldown.getSelected()
      
    self.tocsyPeakList = self.peakListDict[name]
    self.setButtonStates()

  def setNoesy3d(self, index, name=None):

    if not name:
      name = self.noesy3dPulldown.getSelected()
      
    self.noesy3dPeakList = self.peakListDict[name]
    self.setButtonStates()

  def getHsqcs(self):

    names = []
    spectra = getSpectraByType(self.project, 'HSQC')
    for spectrum in spectra:
      for peakList in spectrum.peakLists:
        name = '%s:%s:%s' % (spectrum.experiment.name,spectrum.name,peakList.serial)
        names.append( name )
        self.peakListDict[name] = peakList
        if not self.hsqcPeakList:
          self.hsqcPeakList = peakList
    
    return names

  def setHsqc(self, index, name=None):

    if not name:
      name = self.hsqcPulldown.getSelected()
      
    self.hsqcPeakList  = self.peakListDict[name]
    self.setButtonStates()

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

  def getSpecFreq(self, obj):

    self.specFreqEntry.set(self.specFreq)

  def getMaxIter(self, obj):

    self.maxIterEntry.set(self.maxIter)

  def getMixTime(self, obj):

    self.mixTimeEntry.set(self.mixTime)

  def getCorrTime(self, obj):

    self.corrTimeEntry.set(self.corrTime)

  def getLeakRate(self, obj):

    self.leakRateEntry.set(self.leakRate)

  def getMaxIntens(self, obj):

    self.maxIntensEntry.set(self.maxIntens)

  def setSpecFreq(self, event):

    value = self.specFreqEntry.get()
    if value is not None:
      self.specFreq = value
      
    self.updateMidgeParams()

  def setMaxIter(self, event):

    value = self.maxIterEntry.get()
    if value is not None:
      self.maxIter = value
      
    self.updateMidgeParams()

  def setMixTime(self, event):

    value = self.mixTimeEntry.get()
    if value is not None:
      self.mixTime = value
      
    self.updateMidgeParams()

  def setCorrTime(self, event):

    value = self.corrTimeEntry.get()
    if value is not None:
      self.corrTime = value
      
    self.updateMidgeParams()

  def setLeakRate(self, event):

    value = self.leakRateEntry.get()
    if value is not None:
      self.leakRate = value
      
    self.updateMidgeParams()
 
  def setMaxIntens(self, event):

    value = self.maxIntensEntry.get()
    if value is not None:
      self.maxIntens = value
      
    self.updateMidgeParams()

  def destroy(self):

    BasePopup.destroy(self)

"""
======================COPYRIGHT/LICENSE START==========================

MidgePopup.py: Part of the CcpNmr Clouds program

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
from ccp.api.nmr import Nmr

from memops.editor.BasePopup    import BasePopup
from memops.gui.ButtonList      import ButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.IntEntry        import IntEntry
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Label           import Label
from memops.gui.LabelFrame      import LabelFrame
from memops.gui.PulldownMenu    import PulldownMenu
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.Util            import createDismissHelpButtonList

from ccpnmr.analysis.core.ExperimentBasic import getSpectraByType, getThroughSpacePeakLists
from ccpnmr.analysis.core.StructureBasic  import getAtomSetsDistance

from ccpnmr.clouds.ResonanceIdentification import getCloudsResonanceList, makeNoeAdcs, constrainSpinSystems
from ccpnmr.clouds.NoeRelaxation           import optimiseRelaxation, disambiguateNoesyPeaks
from ccpnmr.clouds.NoeMatrix               import getNoeMatrixFromPeaks

class MidgePopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.guiParent = parent
    self.project   = parent.getProject()
    self.waiting   = 0
    self.specFreq  = 800.13
    self.maxIter   = 15
    self.mixTime   = 60
    self.corrTime  = 11.5
    self.leakRate  = 2.0
    self.ratioHD   = 0.9
    self.peakListDict    = {}
    self.peakListDict3d  = {}
    self.noesyPeakList   = None
    self.noesy3dPeakList = None
    self.carbonLabel     = 0
    self.nitrogenLabel   = 1
    self.noesyPeakList1  = None
    self.noesyPeakList2  = None
    self.noesyPeakList3  = None
    self.noesyPeakList3d = None
    
    self.resonances = None
    self.noesyPeaks = None
    self.distanceConstraintList = None
    self.antiDistConstraintList = None
    self.adcAtomTypes = None
    self.structure    = None
    
    BasePopup.__init__(self, parent, title="Relaxation Matrix Optimisation", **kw)
  
  def body(self, guiFrame):

    self.specFreqEntry  = IntEntry  (self,text=self.specFreq, width=8,returnCallback = self.setSpecFreq)
    self.maxIterEntry   = IntEntry  (self,text=self.maxIter,  width=8,returnCallback = self.setMaxIter )
    self.mixTimeEntry   = FloatEntry(self,text=self.mixTime,  width=8,returnCallback = self.setMixTime )
    self.corrTimeEntry  = FloatEntry(self,text=self.corrTime, width=8,returnCallback = self.setCorrTime)
    self.leakRateEntry  = FloatEntry(self,text=self.leakRate, width=8,returnCallback = self.setLeakRate)
        
    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(1, weight=1)

    row = 0
    labelFrame0 = LabelFrame(guiFrame, text='Input data')
    labelFrame0.grid(row=row, column=0, sticky=Tkinter.NSEW)
    labelFrame0.grid_columnconfigure(3, weight=1)
    
    label = Label(labelFrame0, text='Assigned NOESY spectrum')
    label.grid(row=0, column=0, sticky=Tkinter.NW)
    self.noesyPulldown = PulldownMenu(labelFrame0, entries=self.getNoesys(), callback=self.setNoesy, selected_index=0, do_initial_callback=0)
    self.noesyPulldown.grid(row=0, column=1,sticky=Tkinter.NW )
    
    label = Label(labelFrame0, text='H/D ratio: ')
    label.grid(row=0, column=2, sticky=Tkinter.NW)
    self.ratioHDEntry = FloatEntry(labelFrame0,text=self.ratioHD,  width=6)
    self.ratioHDEntry.grid(row=0, column=3, sticky=Tkinter.NW)
 
    label = Label(labelFrame0, text='NOESY spectrum 1:')
    label.grid(row=1, column=0, sticky=Tkinter.NW)
    self.tmix1Pulldown = PulldownMenu(labelFrame0, entries=self.getNoesys(), callback=self.setNoesy1, selected_index=-0, do_initial_callback=0)
    self.tmix1Pulldown.grid(row=1, column=1,sticky=Tkinter.NW )
    label = Label(labelFrame0, text='Tmix (ms): ')
    label.grid(row=1, column=2, sticky=Tkinter.NW)
    self.tmix1Entry = FloatEntry(labelFrame0,text=60,  width=6)
    self.tmix1Entry.grid(row=1, column=3, sticky=Tkinter.NW)
    
    label = Label(labelFrame0, text='NOESY spectrum 2:')
    label.grid(row=2, column=0, sticky=Tkinter.NW)
    self.tmix2Pulldown = PulldownMenu(labelFrame0, entries=self.getNoesys(), callback=self.setNoesy2, selected_index=0, do_initial_callback=0)
    self.tmix2Pulldown.grid(row=2, column=1,sticky=Tkinter.NW )
    label = Label(labelFrame0, text='Tmix (ms): ')
    label.grid(row=2, column=2, sticky=Tkinter.NW)
    self.tmix2Entry = FloatEntry(labelFrame0,text=120,  width=6)
    self.tmix2Entry.grid(row=2, column=3, sticky=Tkinter.NW)
    
    label = Label(labelFrame0, text='NOESY spectrum 3:')
    label.grid(row=3, column=0, sticky=Tkinter.NW)
    self.tmix3Pulldown = PulldownMenu(labelFrame0, entries=self.getNoesys(), callback=self.setNoesy3, selected_index=0, do_initial_callback=0)
    self.tmix3Pulldown.grid(row=3, column=1,sticky=Tkinter.NW )
    label = Label(labelFrame0, text='Tmix (ms): ')
    label.grid(row=3, column=2, sticky=Tkinter.NW)
    self.tmix3Entry = FloatEntry(labelFrame0,text=200,  width=6)
    self.tmix3Entry.grid(row=3, column=3, sticky=Tkinter.NW)
    
    label = Label(labelFrame0, text='3D NOESY:')
    label.grid(row=4, column=0, sticky=Tkinter.NW)
    self.noesy3dPulldown = PulldownMenu(labelFrame0, entries=self.getNoesys3d(), callback=self.setNoesy3d, selected_index=0, do_initial_callback=0)
    self.noesy3dPulldown.grid(row=4, column=1,sticky=Tkinter.NW )

    
    label10 = Label(labelFrame0, text='Num peaks:')
    label10.grid(row=5, column=0, sticky=Tkinter.NW)
    self.numPeaksLabel = Label(labelFrame0, text='0')
    self.numPeaksLabel.grid(row=5, column=1, sticky=Tkinter.NW)

    label11 = Label(labelFrame0, text='Num resonances:')
    label11.grid(row=5, column=2, sticky=Tkinter.NW)
    self.numResonancesLabel = Label(labelFrame0, text='0')
    self.numResonancesLabel.grid(row=5, column=3, sticky=Tkinter.NW)
    
    row += 1
    labelFrame1 = LabelFrame(guiFrame, text='Parameters')
    labelFrame1.grid(row=row, column=0, sticky=Tkinter.NSEW)
    labelFrame1.grid_columnconfigure(3, weight=1)
    
    label = Label(labelFrame1, text='15N labelled sample:')
    label.grid(row=0, column=0, sticky=Tkinter.NW)
    self.nitrogenSelect = CheckButton(labelFrame1, callback=self.setNitrogenLabel)
    self.nitrogenSelect.grid(row=0, column=1, sticky=Tkinter.W)
    self.nitrogenSelect.set(1)
  
    label = Label(labelFrame1, text='13C labelled sample:')
    label.grid(row=0, column=2, sticky=Tkinter.NW)
    self.carbonSelect = CheckButton(labelFrame1, callback=self.setCarbonLabel)
    self.carbonSelect.grid(row=0, column=3, sticky=Tkinter.W)
    self.carbonSelect.set(0)

    labelFrame1.grid_rowconfigure(1, weight=1)
    data = [self.specFreq,self.maxIter,self.mixTime,self.corrTime,self.leakRate]
    colHeadings      = ['Spectrometer\nfrequency','Max\niterations','Mixing\ntime (ms)','Correl.\ntime (ns)','Leak\nrate']
    editWidgets      = [self.specFreqEntry,self.maxIterEntry,self.mixTimeEntry,self.corrTimeEntry,self.leakRateEntry,]
    editGetCallbacks = [self.getSpecFreq,  self.getMaxIter,  self.getMixTime,  self.getCorrTime,  self.getLeakRate,  ]
    editSetCallbacks = [self.setSpecFreq,  self.setMaxIter,  self.setMixTime,  self.setCorrTime,  self.setLeakRate,  ]
    self.midgeParamsMatrix = ScrolledMatrix(labelFrame1, editSetCallbacks=editSetCallbacks, editGetCallbacks=editGetCallbacks, editWidgets=editWidgets,
                                            maxRows=1, initialCols=5, headingList=colHeadings, callback=None, objectList=['None',], textMatrix=[data,])
    self.midgeParamsMatrix.grid(row=1, column=0, columnspan=4, sticky=Tkinter.NSEW)

    label10 = Label(labelFrame1, text='Benchmark structure')
    label10.grid(row=2, column=0, sticky=Tkinter.NW)
    self.structurePulldown = PulldownMenu(labelFrame1, entries=self.getStructures(), callback=self.setStructure, selected_index=0, do_initial_callback=0)
    self.structurePulldown.grid(row=2, column=1,sticky=Tkinter.NW )

    label11 = Label(labelFrame1, text='ADC atom types:')
    label11.grid(row=2, column=2, sticky=Tkinter.NW)
    self.adcAtomsPulldown = PulldownMenu(labelFrame1, entries=self.getAdcAtomTypes(), callback=self.setAdcAtomTypes, selected_index=0, do_initial_callback=0)
    self.adcAtomsPulldown.grid(row=2, column=3,sticky=Tkinter.NW )


    row += 1
    labelFrame2 = LabelFrame(guiFrame, text='Output')
    labelFrame2.grid(row=row, column=0, sticky=Tkinter.NSEW)
    labelFrame2.grid_columnconfigure(3, weight=1)
    
    label20 = Label(labelFrame2, text='Distance constraints:')
    label20.grid(row=0, column=0, sticky=Tkinter.NW)
    self.distConstrLabel = Label(labelFrame2, text='0')
    self.distConstrLabel.grid(row=0, column=1, sticky=Tkinter.NW)

    label21 = Label(labelFrame2, text='Anti-distance constraints:')
    label21.grid(row=0, column=2, sticky=Tkinter.NW)
    self.antiConstrLabel = Label(labelFrame2, text='0')
    self.antiConstrLabel.grid(row=0, column=3, sticky=Tkinter.NW)

    texts    = ['Calculate distances','Show distance\nconstraints','Show anti-distance\nconstraints']
    commands = [self.calculateDistances,self.showConstraints,self.showAntiConstraints]
    self.midgeButtons = ButtonList(labelFrame2,expands=1,texts=texts,commands=commands)
    self.midgeButtons.grid(row=1, column=0, columnspan=4, sticky=Tkinter.NSEW)
  
    row += 1
    self.bottomButtons = createDismissHelpButtonList(guiFrame,expands=0,help_url=None)
    self.bottomButtons.grid(row=row, column=0, columnspan=4, sticky=Tkinter.EW)
    
    self.getPeaks()
    self.getResonances()
    self.update()    
    
    self.geometry('600x400')

  def setCarbonLabel(self, boolean):
 
    self.carbonLabel = boolean

  def setNitrogenLabel(self, boolean):

    self.nitrogenLabel = boolean

  def update(self):
     
    if self.resonances and ((self.noesyPeaks and self.noesyPeakList1 and self.noesyPeakList2 and self.noesyPeakList3) or self.noesyPeakList3d):
      self.midgeButtons.buttons[0].enable()
    else:
      self.midgeButtons.buttons[0].disable()
     
    if self.distanceConstraintList:
      self.distConstrLabel.set(str(len(self.distanceConstraintList.constraints)))
      self.midgeButtons.buttons[1].enable()
    else:
      self.distConstrLabel.set('')
      self.midgeButtons.buttons[1].disable()
    
    if self.antiDistConstraintList:
      self.antiConstrLabel.set(str(len(self.antiDistConstraintList.constraints)))
      self.midgeButtons.buttons[2].enable()
    else:
      self.antiConstrLabel.set('')
      self.midgeButtons.buttons[2].disable()
     
    if self.resonances:
      self.numResonancesLabel.set(str(len(self.resonances)))
    else:
      self.numResonancesLabel.set('')
     
    if self.noesyPeaks:
      self.numPeaksLabel.set(str(len(self.noesyPeaks)))
    else:
      self.numPeaksLabel.set('')
      
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
    
    return ['<None>','HN','HN HA','HN HA HB']
  
  def setAdcAtomTypes(self, index, name=None):
  
    if name is None:
      name = self.adcAtomsPulldown.getSelected()
    
    if name == '<None>':
      name = None
    
    self.adcAtomTypes = name
  
  def getResonances(self):
  
    resonanceDict = {}
    if self.noesyPeaks:
      for peak in self.noesyPeaks:
        for peakDim in peak.peakDims:
          for contrib in peakDim.peakDimContribs:
            resonanceDict[contrib.resonance] = 1
            # TBD: Set resonance.name for typing
         
    self.resonances = resonanceDict.keys()
  
  def getPeaks(self):
     
    if self.noesyPeakList:
      self.noesyPeaks = self.noesyPeakList.sortedPeaks()
     
  def calculateDistances(self):

    resonances = list(self.resonances)
    
    resDict = {}
    for resonance in resonances:
      resDict[resonance.serial] = resonance
    
    ratioHD = self.ratioHDEntry.get() or self.ratioHD
    
    tmix1 = self.tmix1Entry.get() or 60
    tmix2 = self.tmix2Entry.get() or 120
    tmix3 = self.tmix3Entry.get() or 200
    
    data = [(tmix1,self.noesyPeakList1), (tmix2,self.noesyPeakList2), (tmix3,self.noesyPeakList3)]
    data.sort()
    
    mixingTimes = [x[0] for x in data]
    peakLists   = [x[1] for x in data]
    
    # get a clean, symmetric and normalised NOE matrix
    noeMatrix = getNoeMatrixFromPeaks(self.noesyPeaks, resonances, peakLists, mixingTimes, ratioHD=ratioHD, analysis=self.guiParent)
    
    # optimiseRelaxation will remove unconstrained resonances
    self.distanceConstraintList, resonances = optimiseRelaxation(resonances,noeMatrix,self.mixTime,self.specFreq,
                                                                 self.corrTime,self.leakRate,self.carbonLabel,
                                                                 self.nitrogenLabel,maxIter=self.maxIter)


    #constrainSpinSystems(self.distanceConstraintList)
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
    
    
    if self.distanceConstraintList and self.noesyPeakList:
      constraintHead = self.distanceConstraintList.nmrConstraintStore
      
      if self.adcAtomTypes:
        adcDict = {'HN':['H'],
                   'HN HA':['H','HA','HA1','HA2'],
                   'HN HA HB':['H','HA','HA1','HA2','HB','HB2','HB3']}
 
        allowedAtomTypes= adcDict[self.adcAtomTypes]
 
        print "Making ADCs"
        self.antiDistConstraintList = makeNoeAdcs(resonances[:],self.noesyPeakList.dataSource,constraintHead,allowedAtomTypes=allowedAtomTypes)
        print "Done ADCs"
 
    
      if self.structure:
      
        N = len(self.resonances)
        sigmas = [[] for i in range(N) ]
        for i in range(N):
          sigmas[i] = [0.0 for j in range(N)]
          
        for constraint in self.distanceConstraintList.constraints:
          item = constraint.findFirstItem()
          resonances = list(item.resonances)
          
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
          
          constraint.setOrigData(distance)
       
    self.update()
  
  
  def showConstraints(self):
  
    if self.distanceConstraintList:
      self.guiParent.browseConstraints(constraintList=self.distanceConstraintList)
  
  def showAntiConstraints(self):
  
    if self.antiDistConstraintList:
      self.guiParent.browseConstraints(constraintList=self.antiDistConstraintList)

  def getNoesys3d(self):

    peakLists = getThroughSpacePeakLists(self.project)
  
    names = ['<None>',]
    for peakList in peakLists:
      spectrum = peakList.dataSource
      if spectrum.numDim != 3:
        continue

      name = '%s:%s:%s' % (spectrum.experiment.name,spectrum.name,peakList.serial)
      names.append( name )
      self.peakListDict3d[name] = peakList
      if not self.noesyPeakList:
        self.noesyPeakList = peakList
    
    return names

  def getNoesys(self):
    
    peakLists = getThroughSpacePeakLists(self.project)
    
    names = ['<None>',]
    for peakList in peakLists:
      spectrum = peakList.dataSource
      name = '%s:%s:%s' % (spectrum.experiment.name,spectrum.name,peakList.serial)
      names.append( name )
      self.peakListDict[name] = peakList
      
      if not self.noesyPeakList:
        self.noesyPeakList = peakList
    
    return names
 
  def setNoesy(self, index, name=None):

    if not name:
      name = self.noesyPulldown.getSelected()
    
    if name == '<None>':
      self.noesyPeakList = None
    
    else:  
      self.noesyPeakList = self.peakListDict[name]
    
    self.getPeaks()
    self.getResonances()
    self.update()
  
  def setNoesy1(self, index, name=None):

    if not name:
      name = self.tmix1Pulldown.getSelected()

    if name != '<None>':
      self.noesyPeakList1 = self.peakListDict[name]
    else:
      self.noesyPeakList1 = None
    
    self.update()
   
  def setNoesy2(self, index, name=None):

    if not name:
      name = self.tmix2Pulldown.getSelected()
      
    if name != '<None>':
      self.noesyPeakList2 = self.peakListDict[name]
    else:
      self.noesyPeakList2 = None
    
    self.update()
    
  def setNoesy3(self, index, name=None):

    if not name:
      name = self.tmix3Pulldown.getSelected()
 
    if name != '<None>':
      self.noesyPeakList3 = self.peakListDict[name]
    else:
      self.noesyPeakList3 = None
    
    self.update()
  
  def setNoesy3d(self, index, name=None):

    if not name:
      name = self.noesy3dPulldown.getSelected()
 
    if name != '<None>':
      self.noesyPeakList3d = self.peakListDict3d[name]
      self.noesyPeaks      = self.noesyPeakList3d.sortedPeaks()
    
    else:
      self.noesyPeakList3d = None
      self.noesyPeaks      = []
 
    self.getResonances()
    self.update()
  
  def updateMidgeParams(self):
    
    data = [self.specFreq,self.maxIter,self.mixTime,self.corrTime,self.leakRate]

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
 
  def destroy(self):

    BasePopup.destroy(self)

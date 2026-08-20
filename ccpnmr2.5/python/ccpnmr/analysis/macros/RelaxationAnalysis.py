
"""
======================COPYRIGHT/LICENSE START==========================

ScrolledGraph.py: <write function here>

Copyright (C) 2009 Tim Stevens (University of Cambridge)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

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

#T1, T2, HetNOE graph + errs

from math import sqrt, exp
from random import random, randint

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.AssignmentBasic import makeResonanceGuiName

from memops.gui.Frame import Frame
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Label import Label
from memops.gui.LabelDivider import LabelDivider
from memops.gui.LabelFrame import LabelFrame
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.ButtonList import ButtonList, UtilityButtonList
import Tkinter

MS_UNIT_MULTIPLIERS = {'s':1000.0,'ms':1.0,'us':1e-3,'ns':1e-6,}

T1 = u'T\u2081'
T2 = u'T\u2082'
S2 = u'S\u00B2'
Te = u'\u03C4e'

def relaxationAnalysisMacro(argServer):

  popup = RelaxationAnalysisPopup(argServer.parent)
  popup.open()

class RelaxationAnalysisPopup(BasePopup):

  def __init__(self, parent, *args, **kw):
    
    self.guiParent = parent
    self.t1List = None
    self.t2List = None
    self.noeList = None
    self.waiting = False
   
    BasePopup.__init__(self, parent, title="Relaxation Analyses", **kw)

    periodicTable = self.project.currentChemElementStore
    
    nitrogen = periodicTable.findFirstChemElement(elementSymbol='N')
    hydrogen = periodicTable.findFirstChemElement(elementSymbol='H')

  def body(self, guiFrame):

    self.geometry('700x840')
    
    guiFrame.expandGrid(1,0)
    
    # Top frame
    
    frame = Frame(guiFrame, grid=(0,0))
    frame.expandGrid(None,8)
    
    label = Label(frame, text=u' %s List:' % T1, grid=(0,0))
    self.t1Pulldown = PulldownList(frame, callback=self.selectT1List, grid=(0,1))
    
    label = Label(frame, text=u'  %s List:' % T2, grid=(0,2))
    self.t2Pulldown = PulldownList(frame, callback=self.selectT2List, grid=(0,3))
    
    label = Label(frame, text='  NOE List:', grid=(0,4))
    self.noePulldown = PulldownList(frame, callback=self.selectNoeList, grid=(0,5))
 
    label = Label(frame, text='  Spectrometer Freq (MHz):', grid=(0,6))
    self.sfEntry = FloatEntry(frame, grid=(0,7), text=600.00, width=6)

    UtilityButtonList(frame, grid=(0,9), sticky='e')
    
    # Tabs
    
    options = [u'%s vs %s Scatter' % (T1, T2),
               u'%s,%s & NOE Graphs' % (T1, T2) ,
               u'%s/%s Graph' % (T1, T2),
               u'%s Estimate Graph' % S2,
               'Options']
    self.tabbedFrame = TabbedFrame(guiFrame, options=options,
                                   callback=self.toggleTab, grid=(1,0))
    frameA, frameD, frameC, frameE, frameB = self.tabbedFrame.frames

    # T1 vs T2 Graph

    frameA.expandGrid(0,0)
 
    self.t1t2Graph = T1VersesT2Plot(frameA, grid=(0,0))

    # T1 & T2 Graph

    frameD.expandGrid(0,0)
    frameD.expandGrid(1,0)
    frameD.expandGrid(2,0)
 
    self.t1Graph = MeasurementPlot(frameD, T1, grid=(0,0))
    self.t2Graph = MeasurementPlot(frameD, T2, grid=(1,0))
    self.noeGraph = NoePlot(frameD, 'NOE', grid=(2,0))

    # T1 over T2 Graph

    frameC.expandGrid(0,0)
 
    self.t1t2GraphB = T1OverT2Plot(frameC, grid=(0,0))
    
    # Order params graph
    
    frameE.expandGrid(0,0)
    frameE.expandGrid(1,0)
    frameE.expandGrid(2,0)
    
    self.s2Graph = ScrolledGraph(frameE, title=u'%s vs Residue Sequence' % S2,
                                 xLabel='Residue number', yLabel=S2,
                                 width=500, height=150, graphType='histogram', 
                                 xGrid=True, yGrid=False, grid=(0,0),
                                 dataColors=['#0000A0','#808000'],
                                 dataNames=['Isotropic',])

    self.teGraph = ScrolledGraph(frameE, title=u'%s vs Residue Sequence' % Te,
                                 xLabel='Residue number', yLabel=u'%s (ps)' % Te,
                                 width=500, height=150, graphType='histogram',
                                 xGrid=True, yGrid=False, grid=(1,0),
                                 dataColors=['#008000',])
                                 
    self.rexGraph = ScrolledGraph(frameE, title=u'%s vs Residue Sequence' % 'Rex',
                                  xLabel='Residue number', yLabel='Rex',
                                  width=500, height=150, graphType='histogram',
                                  xGrid=True, yGrid=False, grid=(2,0),
                                  dataColors=['#900000',])
    
    # Options
    
    frameB.expandGrid(4,2)
    
    frame = LabelFrame(frameB, text='Physical Params', grid=(0,0))
    frame.expandGrid(None,3)
   
    label = Label(frame, text=u'N-H bond length (\u00C5)', grid=(0,0))
    self.lenNhEntry = FloatEntry(frame, text=1.015, grid=(0,1), width=8)

    label = Label(frame, text=u'Internal correlation\ntime, \u03C4e (ps)', grid=(1,0))
    self.ictEntry = FloatEntry(frame, text=50.0, grid=(1,1), width=8)
    
    label = Label(frame, text=u'15N Chemical Shift\nAnisotopy,\u0394N (ppm)',
                  grid=(2,0))
    self.csaNEntry = FloatEntry(frame, text=-160.0, grid=(2,1), width=8)
    
    frame = LabelFrame(frameB, text=u'%s vs %s Scatter' % (T1, T2), grid=(1,0))
    
    label = Label(frame, text='Max cluster difference (ms):', grid=(0,0))
    self.clusterDictEntry = FloatEntry(frame, text=20.0, grid=(0,1), width=8)    

    label = Label(frame, text='Min cluster size:', grid=(1,0))
    self.clusterSizeEntry = FloatEntry(frame, text=5, grid=(1,1), width=8)    

    label = Label(frame, text=u'Min graph %s (ms):' % T1, grid=(2,0))
    self.minT1Entry = FloatEntry(frame, text=300.0, grid=(2,1), width=8)    

    label = Label(frame, text=u'Max graph %s (ms):' % T1, grid=(3,0))
    self.maxT1Entry = FloatEntry(frame, text=1000.0, grid=(3,1), width=8)    

    label = Label(frame, text=u'Min graph %s (ms):' % T2, grid=(4,0))
    self.minT2Entry = FloatEntry(frame, text=0.0, grid=(4,1), width=8)    

    label = Label(frame, text=u'Max graph %s (ms):' % T2, grid=(5,0))
    self.maxT2Entry = FloatEntry(frame, text=600.0, grid=(5,1), width=8)    
    
    frame = LabelFrame(frameB, text=u'%s Contours' % S2, grid=(0,1))
    frame.expandGrid(4,3)

    label = Label(frame, text='(Order Parameter Lines)',
                  grid=(0,0), gridSpan=(1,2))

    label = Label(frame, text='Min value:', grid=(1,0))
    self.minS2Entry = FloatEntry(frame, text=0.3, grid=(1,1), width=8)
    
    label = Label(frame, text='Max value:', grid=(2,0))
    self.maxS2Entry = FloatEntry(frame, text=1.0, grid=(2,1), width=8)
    
    label = Label(frame, text='Step:', grid=(3,0))
    self.stepS2Entry = FloatEntry(frame, text=0.1, grid=(3,1), width=8)

    frame = LabelFrame(frameB, text=u'\u03C4m Contours', grid=(1,1))
    frame.expandGrid(4,3)

    label = Label(frame, text='(Rotational Correlation Time Lines)',
                  grid=(0,0), gridSpan=(1,2))
                  
    label = Label(frame, text='Min value (ns):', grid=(1,0))
    self.minRctEntry = FloatEntry(frame, text=5.0, grid=(1,1), width=8)
    
    label = Label(frame, text='Max value (ns):', grid=(2,0))
    self.maxRctEntry = FloatEntry(frame, text=14.0, grid=(2,1), width=8)
    
    label = Label(frame, text='Step (ns):', grid=(3,0))
    self.stepRctEntry = FloatEntry(frame, text=1.0, grid=(3,1), width=8)

    # Bottom frame
 
    texts = [u'Show %s Table' % T1,u'Show %s Table' % T2, u'Estimate %s' % S2]
    commands = [self.showT1List, self.showT2List, self.mc]
    buttonList = ButtonList(guiFrame, grid=(2,0), texts=texts, commands=commands)

    # Update

    self.updateRelaxationLists()

    self.drawAfter()

    self.administerNotifiers(self.registerNotify)
  
  def open(self):
  
    self.updateRelaxationLists()
    self.drawAfter()
  
    BasePopup.open(self)
  
  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
    
    BasePopup.destroy(self)
  
  def getSanitisedParams(self):
  
     cDist = self.clusterDictEntry.get() or 0.0
     lenNh = max(min(self.lenNhEntry.get() or 1.015, 1.2), 0.8)
     ict = self.ictEntry.get() or 0.0
     csaN = self.csaNEntry.get() or -180.0
     
     if csaN > 0:
       csaN *= -1.0
     
     minS2 = max(min(self.minS2Entry.get() or 0.0, 0.9), 0.0)
     maxS2 = max(min(self.maxS2Entry.get() or 1.0, 1.0), 0.1)
     
     if minS2 > maxS2:
       maxS2, minS2 = minS2, maxS2
     
     stepS2 = max(self.stepS2Entry.get() or 0.1, (maxS2-minS2)/100)
     
     minRct = self.minRctEntry.get()
     maxRct = self.maxRctEntry.get()
     
     if minRct > maxRct:
       maxRct, minRct = minRct, maxRct
     
     stepRct = max(self.stepRctEntry.get() or 0.1, (maxRct-minRct)/100)
     
     cluster = max(abs(self.clusterSizeEntry.get() or 10), 1)
     
     minT1 = self.minT1Entry.get() or 300.0
     maxT1 = self.maxT1Entry.get() or 1000.0
     minT2 = self.minT2Entry.get() or 0.0
     maxT2 = self.maxT2Entry.get() or 600.0
   
     
     if minT1 < 0:
       minT1 = 0.0

     if minT2 < 0:
       minT2 = 0.0


     if maxT1 < 0:
       maxT1 = 0.0

     if maxT2 < 0:
       maxT2 = 0.0


     if minT1 > maxT1:
       minT1, maxT1 = maxT1, minT1
       
     if minT2 > maxT2:
       minT2, maxT2 = maxT2, minT2


     if minT1 == maxT1:
       maxT1 += 100
     if minT2 == maxT2:
       maxT2 += 100
     
     
     self.clusterSizeEntry.set(cluster)
     self.minT1Entry.set(minT1)
     self.maxT1Entry.set(maxT1)
     self.minT2Entry.set(minT2)
     self.maxT2Entry.set(maxT2)
     self.clusterDictEntry.set(cDist)
     self.lenNhEntry.set(lenNh)
     self.ictEntry.set(ict)
     self.csaNEntry.set(csaN)
     self.minS2Entry.set(minS2)
     self.maxS2Entry.set(maxS2)
     self.stepS2Entry.set(stepS2)
     self.minRctEntry.set(minRct)
     self.maxRctEntry.set(maxRct)
     self.stepRctEntry.set(stepRct)
     
     paramsS2 = (minS2, maxS2, stepS2)
     paramsRct = (minRct, maxRct, stepRct)
     tBounds = (minT1, maxT1, minT2, maxT2)

     return cluster, cDist, lenNh, ict, csaN,  paramsS2, paramsRct, tBounds

  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete','setName'):
      for clazz in ('ccp.nmr.Nmr.T1List','ccp.nmr.Nmr.T2List'):
        notifyFunc(self.updateRelaxationLists,clazz, func)

    for func in ('__init__', 'delete', 'setValue'):
      notifyFunc(self.updateMeasurementAfter,'ccp.nmr.Nmr.T1', func)
      notifyFunc(self.updateMeasurementAfter,'ccp.nmr.Nmr.T2', func)

  def showT1List(self):
  
    if self.t1List:
      self.guiParent.editMeasurements(self.t1List)

  def showT2List(self):
  
    if self.t2List:
      self.guiParent.editMeasurements(self.t2List)

  def toggleTab(self, index):
  
    if index == 0:
      self.drawAfter()

  def updateMeasurementAfter(self, measurement):
  
    if measurement.parentList in (self.t1List, self.t2List):
      self.drawAfter()
  
  def updateRelaxationLists(self, obj=None):
  
    mLists = self.nmrProject.sortedMeasurementLists()
    

    # T1
    t1List = self.t1List
    t1Lists = [ml for ml in mLists if ml.className == 'T1List']

    index = 0
    names = []
    
    if t1Lists:
      if t1List not in t1Lists:
        t1List = t1Lists[0]
      
      index = t1Lists.index(t1List)
      names = [ml.name or 'List %d' % ml.serial for ml in t1Lists]
    
    else:
      t1List = None
      
    self.selectT1List(t1List)
    self.t1Pulldown.setup(names, t1Lists, index)


    # T2
    t2List = self.t2List
    t2Lists = [ml for ml in mLists if ml.className == 'T2List']
    
    index = 0
    names = []
    
    if t2Lists:
      if t2List not in t2Lists:
        t2List = t2Lists[0]
      
      index = t2Lists.index(t2List)
      names = [ml.name or 'List %d' % ml.serial for ml in t2Lists]
    
    else:
      t2List = None
      
    self.selectT2List(t2List)
    self.t2Pulldown.setup(names, t2Lists, index)
    
    
    # NOE
    noeList = self.noeList
    noeLists = [ml for ml in mLists if ml.className == 'NoeList']
    
    index = 0
    names = []
    
    if noeLists:
      if noeList not in noeLists:
        noeList = noeLists[0]
      
      index = noeLists.index(noeList)
      names = [ml.name or 'List %d' % ml.serial for ml in noeLists]
    
    else:
      noeList = None
      
    self.selectNoeList(noeList)
    self.noePulldown.setup(names, noeLists, index)


  def drawAfter(self):
  
    if self.waiting:
      return
      
    self.waiting = True
    self.after_idle(self.draw)

  def selectT1List(self, t1List):
  
    if t1List is not self.t1List:
      self.t1List = t1List
      self.sfEntry.set(t1List.sf)
      values = [m.value for m in t1List.measurements]
      values.sort()
      minVal = self.minT1Entry.get()
      maxVal = self.maxT1Entry.get()
      minVal0 = values[0]*MS_UNIT_MULTIPLIERS[t1List.unit]
      maxVal0 = values[-1]*MS_UNIT_MULTIPLIERS[t1List.unit]

      if minVal0 < minVal:
        self.minT1Entry.set(max(minVal0 - 100, 0))
      
      if maxVal0 > maxVal:
        self.maxT1Entry.set(maxVal0 + 100)
      
      self.drawAfter()
    
  
  def selectT2List(self, t2List):
  
    if t2List is not self.t2List:
      self.t2List = t2List
      self.sfEntry.set(t2List.sf)
      values = [m.value for m in t2List.measurements]
      values.sort()
      values.sort()
      minVal = self.minT2Entry.get()
      maxVal = self.maxT2Entry.get()
      minVal0 = values[0]*MS_UNIT_MULTIPLIERS[t2List.unit]
      maxVal0 = values[-1]*MS_UNIT_MULTIPLIERS[t2List.unit]

      if minVal0 < minVal:
        self.minT2Entry.set(max(minVal0 - 50, 0))
      
      if maxVal0 > maxVal:
        self.maxT2Entry.set(maxVal0 + 50)
      
      self.drawAfter()
  
  
  def selectNoeList(self, noeList):
  
    if noeList is not self.noeList:
      self.noeList = noeList
      self.sfEntry.set(noeList.sf)
      
      #self.drawAfter()

  def draw(self):
  
    params = self.getSanitisedParams()
    
    cSize, cDist, lenNh, ict, csaN, paramsS2, paramsRct, tBounds = params
 
    sf = abs(self.sfEntry.get() or 500.010)
    
    self.t1t2Graph.update(self.t1List, self.t2List, cSize, cDist, sf, lenNh,
                          ict, csaN, paramsS2, paramsRct, tBounds)
                          
    self.t1t2GraphB.update(self.t1List, self.t2List)
    
    self.t1Graph.update(self.t1List)
    self.t2Graph.update(self.t2List)
    self.noeGraph.update(self.noeList)
      
    self.waiting = False


  def mc(self):
  
    sf = abs(self.sfEntry.get() or 500.010)
    
    if not self.t1List and self.t2List:
      return
    
    chains = set([])
    t1t2Points = {}
    resonancesT1 = {}
    resonancesNoe = {}
    residueResonances = {}

    t1Unit = MS_UNIT_MULTIPLIERS[self.t1List.unit]
    t2Unit = MS_UNIT_MULTIPLIERS[self.t2List.unit]
    
    
    if self.noeList:
      for noe in self.noeList.measurements:
        for resonance in noe.resonances:
          resonancesNoe[resonance] = noe
    
    for t1 in self.t1List.measurements:
      resonancesT1[t1.resonance] = t1
       
    for t2 in self.t2List.measurements:
      resonance = t2.resonance
      if not resonance.resonanceSet:
        continue
      
      residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
      
      t1 = resonancesT1.get(resonance)
      noe = resonancesNoe.get(resonance)
      
      if t1 and noe:
        t1t2Points[residue] = 1e-3*t1Unit*t1.value, 1e-3*t2Unit*t2.value, noe.value
        chains.add(residue.chain)
        residueResonances[residue] = resonance
      
      if t1:
        t1t2Points[residue] = 1e-3*t1Unit*t1.value, 1e-3*t2Unit*t2.value, None
        chains.add(residue.chain)
        residueResonances[residue] = resonance
     
    for chain in chains:
      residues = chain.sortedResidues()
      t1Values = []
      t2Values = []
      noeValues = []
      n = len(residues)
      n2 = n-1
      jj = range(n)
      s2Best = [0.0] * n
      teBest = [0.0] * n
      rexBest = [0.0] * n
      
      print ''
      
      for residue in residues:
        t1,t2, noe = t1t2Points.get(residue, (None, None, None))
        t1Values.append(t1)
        t2Values.append(t2)
        noeValues.append(noe)
      
      t1s = [v for v in t1Values if v]
      t2s = [v for v in t2Values if v]
      noes = [v for v in noeValues if v]
      
      m = len(t1s)
      
      t1 = sum(t1s)/float(m)
      t2 = sum(t2s)/float(m)
      noe = sum(noes)/float(m)
      
      t12 = [v/t2s[i] for i, v in enumerate(t1s) if v]
      t12m = sum(t12)/float(len(t12))
      deltas = [abs(v-t12m)/t12m for v in t12]
      w = [1.0/(v*v) for v in deltas]
      
      t1 = sum([ t1s[j]*w[j] for j in range(m)])/sum(w)
      t2 = sum([ t2s[j]*w[j] for j in range(m)])/sum(w)
      
      if noes:
        noe = sum([ noes[j]*w[j] for j in range(m)])/sum(w)
      else:
        noe = None

      i, ensemble = self.fitT1T2(t1, t2, noe, sf, tmFix=None, teFix=None, s2Fix=None,
                                 tmMin=1e-9, tmMax=100e-9, teMin=1e-12, teMax=1e-10,
                                 rexFix=0.0)
      ensemble.sort()
      score, s2, te0, tm0, rex, t1t, t2t, noet = ensemble[0]
 
      data = (s2, te0*1e12, tm0*1e9, rex, t1, t1t, t2, t2t, noe or 0.0, noet,  score, i)
      print 'Mean A S2:%5.3f Te:%5.1f Tm:%5.3f Rex:%5.3f T1:%5.3f %5.3f T2:%5.3f %5.3f NOE:%5.3f %5.3f %e %6d' % data
      
      rexCheck = 999 # 1.40 * t1/t2
       
      # Prior is mean, then region 
   
      for j in jj:
        t1 = t1Values[j]
        t2 = t2Values[j]
        noe = noeValues[j]
        if t1 is None:
          continue
         
        residue = residues[j]
        print '%3d%s' %  (residue.seqCode, residue.ccpCode),
        
        i, ensemble = self.fitT1T2(t1, t2, noe, sf, tmFix=tm0, teFix=None, s2Fix=None,
                                   tmMin=tm0*0.1, tmMax=tm0*5, teMin=te0/100, teMax=te0*20,
                                   rexFix=0.0, rexMax=15.0)
        ensemble.sort()
        score, s2, te, tm, rex, t1t, t2t, noet = ensemble[0]
        
        if s2 > 0.995:
          i, ensemble = self.fitT1T2(t1, t2, noe, sf, tmFix=tm0, teFix=None, s2Fix=None,
                                     teMin=te/10, teMax=te*10,
                                     rexFix=None, rexMax=15.0)
          ensemble.sort()
          score, s2, te, tm, rex, t1t, t2t, noet = ensemble[0]
         
          
        teBest[j] = te
        s2Best[j] = s2
        rexBest[j] = rex
 
        data = (s2, te*1e12, tm*1e9, rex, t1, t1t, t2, t2t, noe or 0.0, noet, score, i)
        print 'S2:%5.3f Te:%5.1f Tm:%5.3f Rex:%5.3f T1:%5.3f %5.3f T2:%5.3f %5.3f NOE:%5.3f %5.3f %e %6d' % data

      dataSet1 = []
      dataSet2 = []
      dataSet3 = []
      for j in jj:
        residue = residues[j]
        resonance = residueResonances.get(residue)
        
        if not resonance:
          continue
        
        dataSet1.append((residue.seqCode, s2Best[j]))
        dataSet2.append((residue.seqCode, teBest[j]*1e9))
        dataSet3.append((residue.seqCode, rexBest[j]))
        
      self.s2Graph.update([dataSet1,])
      self.teGraph.update([dataSet2,])
      self.rexGraph.update([dataSet3,])
      
      self.tabbedFrame.select(3)

  def fitT1T2(self, t1, t2, noe, sf, tmFix=None, teFix=None, s2Fix=None,
              tmMin=1e-9, tmMax=100e-9, teMin=1e-12, teMax=1e-9,
              rexFix=None, rexMax=15.0, niter=10000):
    
    # Could have multiple sf
    
    rNH = 1.015
    csaN = 160*1e-6 # 160 ppm
    omegaH = sf * 2 * PI * 1e6 # Rad/s
    omegaN = omegaH * GAMMA_N/GAMMA_H
    gammaHN = GAMMA_H/GAMMA_N
     
    A = REDUCED_PERM_VACUUM * REDUCED_PLANK * GAMMA_N * GAMMA_H * 1e30 # Cubic Angstrom adjust
    A /= rNH**3.0
    A = A*A/4.0
    C = omegaN * omegaN * csaN * csaN 
    C /= 3.0
    
    # Init params
    if s2Fix is None:
      s2Start = [x*0.05 for x in range(1,20)]
    else:
      s2Start = [s2Fix,]
    
    if teFix is None:
      teStart = [x*1e-12 for x in [10,30,70,100,300,700,1000]]
    else:
      teStart = [teFix,]
    
    if tmFix is None:
      tmStart = [x*1e-9 for x in range(1,30)]
    else:
      tmStart = [tmFix,]
    
    if rexFix is None:
      rexStart = [float(x) for x in range(int(rexMax))]
    else:
      rexStart = [rexFix,]
    
    # Init ensemble of solutions
    ensemble = []
    for s2 in s2Start:
      for te in teStart:
        for tm in tmStart:
          for rex in rexStart:
             ensemble.append((None, s2, te, tm, rex))
   
    
    # Init scores
    for k, (score, s2, te, tm, rex) in enumerate(ensemble):
      jH   = getSpectralDensity(s2, tm, te, omegaH)
      jN   = getSpectralDensity(s2, tm, te, omegaN)
      jHpN = getSpectralDensity(s2, tm, te, omegaH+omegaN)
      jHmN = getSpectralDensity(s2, tm, te, omegaH-omegaN)
      j0   = getSpectralDensity(s2, tm, te, 0.0)
 
      r1 = A*( (3*jN) + (6*jHpN) + jHmN ) + C*jN
      r2 = 0.5*A*( (4*j0) + (3*jN) + (6*jHpN) + (6*jH) + jHmN ) + C*( 2*j0/3.0 + 0.5*jN ) + rex
 
      t1p = 1.0/r1
      t2p = 1.0/r2
 
 
      d1 = (t1p-t1)/t1
      d2 = (t2p-t2)/t2
 
      score = (d1*d1) + (d2*d2)
      
      if noe is None:
        noep = 0.0
        
      else:
        noep = 1.0 + ( A * gammaHN * t1 * ((6*jHpN)-jHmN) )
        dn = (noep-noe)/2.0
        score += (dn*dn) 
      
      ensemble[k] = (score, s2, te, tm, rex, t1p, t2p, noep)
      
      # Matrix
      # 
      # Pred    Jh        
      # R1
      # R2
      # NOE
      #   r = 1.0/rct + 1.0/ict 
      #   t = 1.0/r
      #   j = (s2*rct) / (1.0 + w*w*rct*rct)
      #   j += ((1.0-s2)*t) / (1.0 + w*w*t*t)
      #       
      #   return j*0.4
    
    ensemble.sort()
    ensemble = ensemble[:10]
    ensembleSize = len(ensemble)
         
    bestScore = 1e99
    
    for i in xrange(niter):
 
      f = i/float(niter)
      f = exp(-10.0*f)
      # Mutate
 
      ensemble.sort()
      prevScore, s2, te, tm, rex, t1p, t2p, noep = ensemble[-1] #Biggest is worst
 
      if ensemble[0][0] < 1e-10:
        break
 
      if not s2Fix:
        #d = ((random() + 0.618 - 1.0) * f) + 1.0
        d = ((random() - 0.382) * f) + 1.0
        s2 = max(0.0, min(1.0, s2*d))
 
      if not tmFix:
        d = ((random() - 0.382) * f) + 1.0
        tm = max(tmMin, min(tmMax, tm*d))
      
      if not teFix:
        d = ((random() - 0.382) * f) + 1.0 
        te = max(teMin, min(teMax, te*d))

      d = ((random() - 0.382) * f) + 1.0 
      rex = max(0.0, min(rexMax, rex*d))
 
      jH   = getSpectralDensity(s2, tm, te, omegaH)
      jN   = getSpectralDensity(s2, tm, te, omegaN)
      jHpN = getSpectralDensity(s2, tm, te, omegaH+omegaN)
      jHmN = getSpectralDensity(s2, tm, te, omegaH-omegaN)
      j0   = getSpectralDensity(s2, tm, te, 0.0)

      r1 = A*( (3*jN) + (6*jHpN) + jHmN ) + C*jN
      r2 = 0.5*A*( (4*j0) + (3*jN) + (6*jHpN) + (6*jH) + jHmN ) + C*( 2*j0/3.0 + 0.5*jN ) + rex
 
      t1p = 1.0/r1
      t2p = 1.0/r2
 
      d1 = (t1p-t1)/t1
      d2 = (t2p-t2)/t2
 
      score = (d1*d1) + (d2*d2)
      if noe is None:
        noep = 0.0
      else:
        noep = 1.0 + ( A * gammaHN * t1 * ((6*jHpN)-jHmN) )
        dn = (noep-noe)/2.0
        score += (dn*dn) 
      
      ratio = exp(prevScore-score)
                      
      if ratio > 1.0: # random():
        ensemble[-1] = (score, s2, te, tm, rex, t1p, t2p, noep)
      else:
        k = randint(0,ensembleSize-1)
        score, s2, te, tm, rex, t1p, t2p, noep = ensemble[k]
        ensemble[-1] = (score, s2, te, tm, rex, t1p, t2p, noep)
    
      
      if score < bestScore:
        bestScore = score
    
    #print bestScore, ', '.join(['%.3e' % x[0] for x in ensemble])
        
    return i, ensemble
    
CLUSTER_COLORS = ('#800000','#000080',
                  '#008000','#808000',
                  '#800080','#008080',
                  '#808080','#000000',
                  '#804000','#004080')
                   
CLUSTER_SYMBOLS = ('circle', 'square', 'triangle')

REDUCED_PLANK = 1.05457162853e-34
PI = 3.14159265358979323846
GAMMA_H = 42.576 * 1e6 * 2 * PI # Hz/T
GAMMA_N = -4.3156 * 1e6 * 2 * PI # Hz/T

REDUCED_PERM_VACUUM = 1e-7

class T1VersesT2Plot(ScrolledGraph):

  def __init__(self, parent, contourColorA='#FFC0B0', contourColorB='#B0C0FF', *args, **kw):
   
    kw['title'] =  u'%s vs %s Relaxation Analysis' % (T1, T2)
    kw['xLabel'] = u'%s (ms)' % T1
    kw['yLabel'] = u'%s (ms)' % T2
    kw['width'] = 500
    kw['height'] = 500
    kw['xGrid'] = False
    kw['yGrid'] = False
    kw['graphType'] = 'scatter'
    
    self.minT1 = 0.0
    self.maxT1 = 1000.0
    self.minT2 = 0.0
    self.maxT2 = 1000.0
    self.contourColorA = contourColorA
    self.contourColorB = contourColorB
    self.orderParamLines = []
    self.rotCorrTimeLines = []
    self.contourItems = set([])
    self.outlierLabels = []
      
    ScrolledGraph.__init__(self, parent, **kw)

  def calcContourLines(self, sf, lenNh, ict, csaN, paramsS2, paramsRct):

    rctLines = []
    s2Lines = []

    minS2, maxS2, stepS2 = paramsS2
    minRct, maxRct, stepRct = paramsRct
    
    omegaH = sf * 2 * PI * 1e6 # Rad/s
    omegaN = omegaH * GAMMA_N/GAMMA_H
    csaN /= 1e6 # Was  ppm
     
    A = REDUCED_PERM_VACUUM * REDUCED_PLANK * GAMMA_N * GAMMA_H * 1e30 # Cubic Angstrom adjust
    A /= lenNh**3.0
    A = A*A/4.0
    C = omegaN * omegaN * csaN * csaN 
    C /= 3.0
    
    rct = maxRct
    ict *= 1e-12 # Seconds
    
    while rct >= minRct:
      line = []
      s2 = 50
 
      while s2 > 0.1:
        rctB = rct * 1e-9 # Seconds
        jH = getSpectralDensity(s2, rctB, ict, omegaH)
        jN = getSpectralDensity(s2, rctB, ict, omegaN)
        jHpN = getSpectralDensity(s2, rctB, ict, omegaH+omegaN)
        jHmN = getSpectralDensity(s2, rctB, ict, omegaH-omegaN)
        j0 = getSpectralDensity(s2, rctB, ict, 0.0)
 
        r1 = A*( (3*jN) + (6*jHpN) + jHmN ) + C*jN
        r2 = 0.5*A*( (4*j0) + (3*jN) + (6*jHpN) + (6*jH) + jHmN ) + C*( 2*j0/3.0 + 0.5*jN )
        t1 = 1e3/r1 # ms
        t2 = 1e3/r2 # ms
        
        line.append((t1,t2))
        
        
        s2 /= 1.4
        
      rctLines.append((rct, line))
      rct -= stepRct
    
    s2 = maxS2
    while s2 >= minS2:
      line = []
      
      rct = 0.05
      while rct < 50:
        rctB = rct * 1e-9 # Nanoseconds
        jH = getSpectralDensity(s2, rctB, ict, omegaH)
        jN = getSpectralDensity(s2, rctB, ict, omegaN)
        jHpN = getSpectralDensity(s2, rctB, ict, omegaH+omegaN)
        jHmN = getSpectralDensity(s2, rctB, ict, omegaH-omegaN)
        j0 = getSpectralDensity(s2, rctB, ict, 0.0)
 
        r1 = A*( (3*jN) + (6*jHpN) + jHmN ) + C*jN
        r2 = 0.5*A*( (4*j0) + (3*jN) + (6*jHpN) + (6*jH) + jHmN ) + C*( 2*j0/3.0 + 0.5*jN )

        t1 = 1e3/r1 # ms
        t2 = 1e3/r2 # ms
        
        line.append((t1,t2))
            
        rct += 0.1
      
      s2Lines.append((s2, line))
      s2 -= stepS2
         
    self.rotCorrTimeLines = rctLines
    self.orderParamLines = s2Lines
      
  def update(self, t1List, t2List, cSize, cDist, sf, lenNh, ict,
             csaN, paramsS2, paramsRct, graphBounds):
    
    self.calcContourLines(sf, lenNh, ict, csaN, paramsS2, paramsRct)
    
    self.minT1, self.maxT1, self.minT2, self.maxT2 = graphBounds
    
    dataSets = []
    dataColors = []
    symbols = []
    outliers = []
    nColors = len(CLUSTER_COLORS)
    nSymbols = len(CLUSTER_SYMBOLS)
    cDist2 = cDist**2
    
    if t1List and t2List:
      t1t2Points = []
      resonancesT1 = {}

      t1Unit = MS_UNIT_MULTIPLIERS[t1List.unit]
      t2Unit = MS_UNIT_MULTIPLIERS[t2List.unit]
      
      for t1 in t1List.measurements:
        resonancesT1[t1.resonance] = t1
         
      for t2 in t2List.measurements:
        resonance = t2.resonance
        t1 = resonancesT1.get(resonance)
        
        if t1:
          xyr = (t1Unit*t1.value,t2Unit*t2.value, resonance)
          t1t2Points.append(xyr)
      
      # Do clustering
      clusters = [[pt,] for pt in t1t2Points]
      
      nMerge = 1
      while nMerge:
        nMerge = 0
        
        for i, clusterA in enumerate(clusters[:-1]):
          for j in range(i+1,len(clusters)):
            clusterB = clusters[j]
            
            for x1, y1, r1 in clusterA:
              for x2, y2, r2 in clusterB:
                dx = x2-x1
                dy = y2-y1
                
                if (dx*dx) + (dy*dy) <= cDist2:
                  clusterA += clusterB
                  clusters[j] = []
                  nMerge += 1
                  break
              
              else:
                continue
              break      
        
        clusters = [c for c in clusters if c]
      
      self.outlierLabels = []
      for c in clusters:
        if len(c) < cSize:
          for x,y,r in c:
            outliers.append((x,y))
            label =  makeResonanceGuiName(r)
            self.outlierLabels.append((label,x,y))
        else:
          data = [xyr[:2] for xyr in c]
          dataSets.append(data)
               
      dataSets.append(outliers)
      
      for i, dataSet in enumerate(dataSets):
        dataSet.sort()
        dataColors.append(CLUSTER_COLORS[i % nColors])
        symbols.append(CLUSTER_SYMBOLS[i % nSymbols])
    
    if outliers:
      symbols[-1] = 'star'
      
    #else:
    #  dataSets.append([(0,0),(1000,1000)])

    
    ScrolledGraph.update(self, dataSets, dataColors, symbols=symbols)


  def configMenu(self):

    # Cut-down from superclass

    zoom_items = [ \
      { 'kind': 'command', 'label': 'In',
      'command': lambda event: self.setZoom(self.zoom*1.2) },
      { 'kind': 'command', 'label': 'Out',
      'command': lambda event: self.setZoom(self.zoom/1.2) },
      { 'kind': 'command', 'label': 'Reset',
      'command': lambda event: self.setZoom(1.0) },
    ]

    symbolsize_items = [ \
      { 'kind': 'command', 'label': '+',
      'command':lambda event: self.setSymbolSize( self.symbolSize+1 ) },
      { 'kind': 'command', 'label': '-',
      'command':lambda event: self.setSymbolSize( max(self.symbolSize-1,1) ) },
      { 'kind': 'command', 'label': 'Small',
      'command':lambda event: self.setSymbolSize(3) },
      { 'kind': 'command', 'label': 'Normal',
      'command':lambda event: self.setSymbolSize(5) },
      { 'kind': 'command', 'label': 'Big',
      'command':lambda event: self.setSymbolSize(7) },
    ]

    _items = [ \
      { 'kind': 'command', 'label': '', 'command':'' },
      { 'kind': 'command', 'label': '', 'command':'' },
      { 'kind': 'command', 'label': '', 'command':'' },
    ]

    items = [
             { 'kind': 'cascade', 'label': 'Zoom',
             'submenu': zoom_items },
             { 'kind': 'cascade', 'label': 'Symbol size',
             'submenu': symbolsize_items },
             { 'kind': 'command', 'label': 'Set Title',
             'command' : self.setGraphTitle },
             { 'kind': 'command', 'label': 'Print to file',
             'command' : self.scrolledCanvas.printCanvas },
            ]
 
    self.scrolledCanvas.menu.setMenuItems(items)

  def draw(self):
  
    self.drawContours()
  
    ScrolledGraph.draw(self)
  
    plotRegion = self.getPlotRegion()

    minX = self.minT1
    maxX = self.maxT1
    minY = self.minT2
    maxY = self.maxT2

    dataRegion = [minX, minY, maxX, maxY]
    
    deltaXplot = plotRegion[2] - plotRegion[0]
    deltaYplot = plotRegion[3] - plotRegion[1] 
    deltaXdata = dataRegion[2] - dataRegion[0]
    deltaYdata = dataRegion[3] - dataRegion[1]
  
    ppvX = deltaXplot/float(deltaXdata)
    ppvY = deltaYplot/float(deltaYdata)
    createText = self.canvas.create_text

    for label, x, y in self.outlierLabels:
      x0 = (x - dataRegion[0])*ppvX
      y0 = (y - dataRegion[1])*ppvY
      x0 += plotRegion[0]
      y0 += plotRegion[1]
      
      item = createText(x0+8,y0, text=label,
                        fill='#000000',anchor='w')
      self.contourItems.add(item)
      
  
  def getMinMaxValues(self, dataSets):
  
    minX,maxX,minY,maxY,maxN,minDx,minDy = ScrolledGraph.getMinMaxValues(self, dataSets)
    
    minX = self.minT1
    maxX = self.maxT1
    minY = self.minT2
    maxY = self.maxT2
    
    return minX,maxX,minY,maxY,maxN,minDx,minDy
  
  def drawContours(self, linewidth=1):
    
    delete = self.canvas.delete
    for item in self.contourItems:
      delete(item)
      
    plotRegion = self.getPlotRegion()

    minX = self.minT1
    maxX = self.maxT1
    minY = self.minT2
    maxY = self.maxT2

    dataRegion = [minX, minY, maxX, maxY]
    color = self.contourColorA
    
    deltaXplot = plotRegion[2] - plotRegion[0]
    deltaYplot = plotRegion[3] - plotRegion[1] 
    deltaXdata = dataRegion[2] - dataRegion[0]
    deltaYdata = dataRegion[3] - dataRegion[1]
  
    ppvX = deltaXplot/float(deltaXdata)
    ppvY = deltaYplot/float(deltaYdata)
  
    s2Labels = []
    rctLabelsR = []
    rctLabelsT = []
  
    s = 0
    for lines in self.orderParamLines, self.rotCorrTimeLines:
      for value, points in lines:
        label = '%1.1f' % value
        if s == 0:
          s2t1t2 = points[:]
          s2t1t2.sort()
          t1, t2 = s2t1t2[0]
          if (minX < t1 < maxX) and (minY < t2 < maxY): 
            s2Labels.append((label, t1, t2)) 
        
        coords = []
 
        for i, (x,y) in enumerate(points):
          
          if i:
            xp, yp = points[i-1]
            dy = y-yp
            dx = x-xp
            g = dy/dx
            
            if (x < maxX) and (xp > maxX):
              y = y + (g*(maxX-x))
              x = maxX
            
              if s == 1 and (y < maxY):
                rctLabelsR.append((label, x,y))
            
            elif (xp < maxX) and (x > maxX):
              y = yp + (g*(maxX-xp))
              x = maxX

              if s == 1 and (y < maxY):
                rctLabelsR.append((label, x,y))
            
            elif (x > minX) and (xp < minX):
              y = y + (g*(minX-x))
              x = minX
            
            elif (xp > minX) and (x < minX):
              y = yp + (g*(minX-xp))
              x = minX
            
            if (y < maxY) and (yp > maxY):
              x = x + ((maxY-y)/g)
              y = maxY

              if s == 1 and (x < maxX):
                rctLabelsT.append((label, x,y))
            
            elif (yp < maxY) and (y > maxY):
              x = xp + ((maxY-yp)/g)
              y = maxY

              if s == 1 and (x < maxX):
                rctLabelsT.append((label, x,y))
            
            elif (y > minY) and (yp < minY):
              x = x + ((minY-y)/g)
              y = minY
            
            elif (yp > minY) and (y < minY):
              x = xp + ((minY-yp)/g)
              y = minY

          if not (minX <= x <= maxX):
            continue 
          if not (minY <= y <= maxY):
            continue 
 
          x0 = (x - dataRegion[0])*ppvX
          y0 = (y - dataRegion[1])*ppvY
          x0 += plotRegion[0]
          y0 += plotRegion[1]
 
          coords.append((x0,y0,None,None))

        items = self.drawLines(coords, color, linewidth)
        self.contourItems.update(items)
        
      color = self.contourColorB
      s = 1
    
       
    createText = self.canvas.create_text
    
    for label, x, y in s2Labels:
      x0 = (x - dataRegion[0])*ppvX
      y0 = (y - dataRegion[1])*ppvY
      x0 += plotRegion[0]
      y0 += plotRegion[1]
      
      item = createText(x0+2,y0, text=label,
                        fill=self.contourColorA,anchor='w')
      self.contourItems.add(item)
      #self.canvas.lift(item)
    
    for label, x, y in rctLabelsR:
      x0 = (x - dataRegion[0])*ppvX
      y0 = (y - dataRegion[1])*ppvY
      x0 += plotRegion[0]
      y0 += plotRegion[1]
      
      item = createText(x0+2,y0, text=label,
                        fill=self.contourColorB,anchor='w')
      self.contourItems.add(item)
    
    for label, x, y in rctLabelsT:
      x0 = (x - dataRegion[0])*ppvX
      y0 = (y - dataRegion[1])*ppvY
      x0 += plotRegion[0]
      y0 += plotRegion[1]
      
      item = createText(x0,y0+2, text=label,
                        fill=self.contourColorB,anchor='s')
      self.contourItems.add(item)


class T1OverT2Plot(ScrolledGraph):

  def __init__(self, parent, *args, **kw):
   
    kw['title'] = u'%s/%s vs Residue Sequence' % (T1, T2)
    kw['xLabel'] = 'Residue number'
    kw['yLabel'] = u'%s/%s' % (T1, T2)
    kw['width'] = 500
    kw['height'] = 500
    kw['xGrid'] = True
    kw['yGrid'] = False
    kw['graphType'] = 'histogram'
    
    self.minRes = 0
    self.maxRes = 100
    self.minVal = 0.0
    self.maxVal = 10.0
      
    ScrolledGraph.__init__(self, parent, **kw)
    
  def update(self, t1List, t2List):
    
    nColors = len(CLUSTER_COLORS)

    dataSets = []
    dataColors = []
    dataNames = []
    seqNums = set()
    values = set()
    
    if t1List and t2List:
      t1t2Points = {}
      resonancesT1 = {}

      t1Unit = MS_UNIT_MULTIPLIERS[t1List.unit]
      t2Unit = MS_UNIT_MULTIPLIERS[t2List.unit]
      
      for t1 in t1List.measurements:
        resonancesT1[t1.resonance] = t1
         
      for t2 in t2List.measurements:
        resonance = t2.resonance
        
        if not resonance.resonanceSet:
          continue
        
        t1 = resonancesT1.get(resonance)
        
        if t1:
          residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
          value = (t1Unit*t1.value) / (t2Unit*t2.value)
          chain = residue.chain
          
          if chain not in t1t2Points:
            t1t2Points[chain] = []
          
          seqCode = residue.seqCode
          t1t2Points[chain].append((seqCode, value))
          seqNums.add(seqCode)
          values.add(value)
       
      chains = t1t2Points.keys()
      molSystems = set([c.molSystem for c in chains])
      
      if len(molSystems) > 1:
        chainLabels = [('Chain %s:%s' % (c.molSystem.code, c.code),  c) for c in chains]
        chainLabels.sort()
        
      else:
        chainLabels = [('Chain %s' % c.code, c) for c in chains]
        chainLabels.sort()
        
      for label, chain in chainLabels:
        data = t1t2Points[chain]
        data.sort()
        dataSets.append(data)
        dataNames.append(label)
      
      for i, dataSet in enumerate(dataSets):
        dataColors.append(CLUSTER_COLORS[i % nColors])
    
      res = [c.sortedResidues() for c in chains]
      
      if seqNums:
        self.minRes = min(seqNums) - 1
        self.maxRes = max(seqNums) + 1
      
      if values:
        self.minVal = 0.0 # min(values)
        self.maxVal = int(max(values) + 1) 
    
    if len(dataNames) < 2:
      dataNames = None
       
    ScrolledGraph.update(self, dataSets, dataColors, dataNames=dataNames)

  def getMinMaxValues(self, dataSets):
  
    minX,maxX,minY,maxY,maxN,minDx,minDy = ScrolledGraph.getMinMaxValues(self, dataSets)
    
    minX = self.minRes
    maxX = self.maxRes
    minY = self.minVal
    maxY = self.maxVal
    
    return minX,maxX,minY,maxY,maxN,minDx,minDy
  
class MeasurementPlot(ScrolledGraph):

  def __init__(self, parent, name, *args, **kw):
   
    kw['title'] = u'%s vs Residue Sequence' % name
    kw['xLabel'] = 'Residue number'
    kw['yLabel'] = u'%s (ms)' % name
    kw['width'] = 500
    kw['height'] = 150
    kw['xGrid'] = True
    kw['yGrid'] = False
    kw['graphType'] = 'histogram'
    
    self.minRes = 0
    self.maxRes = 100
    self.minVal = 0.0
    self.maxVal = 1.0
      
    ScrolledGraph.__init__(self, parent, **kw)
    
  def update(self, mList):
    
    nColors = len(CLUSTER_COLORS)

    dataSets = []
    dataColors = []
    dataNames = []
    seqNums = set()
    values = set()
    
    if mList:
      points = {}
      mUnit = MS_UNIT_MULTIPLIERS[mList.unit]
          
      for measurement in mList.measurements:
        resonance = measurement.resonance
        
        if not resonance.resonanceSet:
          continue
        
        residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
        value = mUnit*measurement.value
        chain = residue.chain
        
        if chain not in points:
          points[chain] = []
        
        seqCode = residue.seqCode
        points[chain].append((seqCode, value))
        seqNums.add(seqCode)
        values.add(value)
       
      chains = points.keys()
      molSystems = set([c.molSystem for c in chains])
      
      if len(molSystems) > 1:
        chainLabels = [('Chain %s:%s' % (c.molSystem.code, c.code),  c) for c in chains]
        chainLabels.sort()
        
      else:
        chainLabels = [('Chain %s' % c.code, c) for c in chains]
        chainLabels.sort()
        
      for label, chain in chainLabels:
        data = points[chain]
        data.sort()
        dataSets.append(data)
        dataNames.append(label)
      
      for i, dataSet in enumerate(dataSets):
        dataColors.append(CLUSTER_COLORS[i % nColors])
    
      res = [c.sortedResidues() for c in chains]
      
      if seqNums:
        self.minRes = min(seqNums) - 1
        self.maxRes = max(seqNums) + 1
      
      if values:
        self.minVal = 0.0 # min(values)
        self.maxVal = int(max(values) + 1) 
    
    if len(dataNames) < 2:
      dataNames = None
       
    ScrolledGraph.update(self, dataSets, dataColors, dataNames=dataNames)

  def getMinMaxValues(self, dataSets):
  
    minX,maxX,minY,maxY,maxN,minDx,minDy = ScrolledGraph.getMinMaxValues(self, dataSets)
    
    minX = self.minRes
    maxX = self.maxRes
    minY = self.minVal
    maxY = self.maxVal
    
    return minX,maxX,minY,maxY,maxN,minDx,minDy
  
class NoePlot(ScrolledGraph):

  def __init__(self, parent, name, *args, **kw):
   
    kw['title'] = '%s vs Residue Sequence' % name
    kw['xLabel'] = 'Residue number'
    kw['yLabel'] = '%s' % name
    kw['width'] = 500
    kw['height'] = 150
    kw['xGrid'] = True
    kw['yGrid'] = False
    kw['graphType'] = 'histogram'
    
    self.minRes = 0
    self.maxRes = 100
    self.minVal = -1.0
    self.maxVal = 1.0
      
    ScrolledGraph.__init__(self, parent, **kw)
    
  def update(self, mList):
    
    nColors = len(CLUSTER_COLORS)

    dataSets = []
    dataColors = []
    dataNames = []
    seqNums = set()
    values = set()
    
    if mList:
      points = {}
           
      for measurement in mList.measurements:
        resonance = list(measurement.resonances)[0]
        
        if not resonance.resonanceSet:
          continue
        
        residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
        value = measurement.value
        chain = residue.chain
        
        if chain not in points:
          points[chain] = []
        
        seqCode = residue.seqCode
        points[chain].append((seqCode, value))
        seqNums.add(seqCode)
        values.add(value)
       
      chains = points.keys()
      molSystems = set([c.molSystem for c in chains])
      
      if len(molSystems) > 1:
        chainLabels = [('Chain %s:%s' % (c.molSystem.code, c.code),  c) for c in chains]
        chainLabels.sort()
        
      else:
        chainLabels = [('Chain %s' % c.code, c) for c in chains]
        chainLabels.sort()
        
      for label, chain in chainLabels:
        data = points[chain]
        data.sort()
        dataSets.append(data)
        dataNames.append(label)
      
      for i, dataSet in enumerate(dataSets):
        dataColors.append(CLUSTER_COLORS[i % nColors])
    
      res = [c.sortedResidues() for c in chains]
      
      if seqNums:
        self.minRes = min(seqNums) - 1
        self.maxRes = max(seqNums) + 1
      
      if values:
        self.minVal = int(min(values))
        self.maxVal = int(max(values) + 1) 
    
    if len(dataNames) < 2:
      dataNames = None
       
    ScrolledGraph.update(self, dataSets, dataColors, dataNames=dataNames)

  def getMinMaxValues(self, dataSets):
  
    minX,maxX,minY,maxY,maxN,minDx,minDy = ScrolledGraph.getMinMaxValues(self, dataSets)
    
    minX = self.minRes
    maxX = self.maxRes
    #minY = self.minVal
    #maxY = self.maxVal
    
    return minX,maxX,minY,maxY,maxN,minDx,minDy
      
def getSpectralDensity(s2, rct, ict, w):
  # Isotropic

  r = 1.0/rct + 1.0/ict 
  t = 1.0/r
  j = (s2*rct) / (1.0 + w*w*rct*rct)
  j += ((1.0-s2)*t) / (1.0 + w*w*t*t)
      
  return j*0.4

"""    
def getSpectralDensityAxial(s2, rct, ict, w):

  r = 1.0/rct + 1.0/ict 
  t = 1.0/r
  
  j  = (a1*t1) / (1.0 + w*w*t1*t1)
  j += (a2*t2) / (1.0 + w*w*t2*t2)
  j += (a3*t3) / (1.0 + w*w*t3*t3)
  j *= s2
  j += ((1.0-s2)*t) / (1.0 + w*w*t*t)
      
  return j*0.4
"""

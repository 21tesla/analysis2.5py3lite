
"""
======================COPYRIGHT/LICENSE START==========================

LinkSideChains.py: Part of the CcpNmr Analysis program

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
# * * * * * * * * * * * * * * TBD * * * * * * * * * * * * * * *
#
# Type spin system as we go?
# Orthogonal window & force windows open?
# Match threshold? Or max num?
# Persistence on settings
#
# Merged Phe CD1/CD2 etc - only one resonance row
# Prochirals options and finding existing peaks and redundant shifts
# Don't add peaks if no intensity
# No strips at all when one atomSite missing

from math import sqrt, exp
from numpy import matrix, zeros, multiply, linalg, identity, ones

from ccpnmr.analysis.core.AssignmentBasic  import getResonanceName, assignResToDim, assignAtomsToRes
from ccpnmr.analysis.popups.BasePopup        import BasePopup
from ccpnmr.analysis.core.ExperimentBasic  import getSpectrumIsotopes, getOnebondDataDims, getDataDimIsotopes
from ccpnmr.analysis.core.ExperimentBasic  import getPrimaryDataDimRef
from ccpnmr.analysis.core.MarkBasic        import createRuler
from ccpnmr.analysis.core.MoleculeBasic    import greekSortAtomNames, getResidueCode
from ccpnmr.analysis.core.PeakBasic        import findPeaks, searchPeaks, pickPeak, snapPeak
from ccpnmr.analysis.core.UnitConverter    import pnt2ppm
from ccpnmr.analysis.core.Util             import getAnalysisPeakList
from ccpnmr.analysis.core.WindowBasic import getWindowPaneName, findOrthogonalWindows, getDataDimAxisMapping
from ccpnmr.analysis.core.WindowBasic import getSpinSystemWindowShifts, displayStrips
from ccpnmr.analysis.core.ChemicalShiftBasic import getAtomProbability, getChemAtomNmrRef

from memops.gui.ButtonList      import ButtonList, UtilityButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.IntEntry        import IntEntry
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Frame           import Frame
from memops.gui.Label           import Label
from memops.gui.LabelFrame      import LabelFrame
from memops.gui.MessageReporter import showWarning, showOkCancel
from memops.gui.ProgressBar     import ProgressBar
from memops.gui.PulldownMenu    import PulldownMenu
from memops.gui.ScrolledMatrix  import ScrolledMatrix


sideChainSortOrder = {'A':0,'B':1,'G':2,'D':3,'E':4,'H':5,'Z':6}
defaultTolerances  = {'1H':0.03,'13C':0.3,'15N':0.6}

root2pi     = 2.506628274631
gausWidth   = 0.3 # Gaussian width relative to tolerance window
twoGaus     = gausWidth*gausWidth*2.0
gausRoot2pi = gausWidth*root2pi
      
def testPopup(argServer):

  popup = LinkSideChainsPopup(argServer.parent)
  popup.open()


def excludeRegion(newRegions, oldRegion, dim, minVal, maxVal):
  reg = list(oldRegion[dim])
  reg.sort()
  start, end = reg
  
  if start < minVal:
  
    if end < minVal:
      newRegions.append(list(oldRegion))        
    
    else:
      region2 = list(oldRegion)
      region2[dim] = (start, minVal)
      newRegions.append(region2)
 
      if end > maxVal:
        region2 = list(oldRegion)
        region2[dim] = (maxVal, end)
        newRegions.append(region2)
  
  elif start < maxVal:
    if end > maxVal:
      region2 = list(oldRegion)
      region2[dim] = (minVal, end)
      newRegions.append(region2)
  
  else: 
    newRegions.append(list(oldRegion))        
  
def getSideChainAtomSets(spinSystem):

  # ************************* SETUP ATOMS **************************

  residue    = spinSystem.residue
  ccpCode    = spinSystem.ccpCode
  project    = spinSystem.root
  linking    = None
  descriptor = None
  molType    = 'protein'
  
  if residue:
    molResidue = residue.molResidue
    linking    = molResidue.linking
    descriptor = molResidue.descriptor
    molType    = molResidue.molType

  chemComp = project.findFirstChemComp(ccpCode=ccpCode, molType=molType)
    
  if linking and descriptor:
    chemCompVar = chemComp.findFirstChemCompVar(linking=linking, descriptor=descriptor)
    
  else:
    chemCompVar = chemComp.findFirstChemCompVar(isDefaultVar=True)
    
    if not chemCompVar:
      chemCompVar = chemComp.findFirstChemCompVar(linking='middle')
  
  if not chemCompVar:
    # Warn
    return

  from ccpnmr.analysis.core.MoleculeBasic import greekSortAtomNames
  
  sideChainAtomSets = set()
  for chemAtom in chemCompVar.chemAtoms:
    if chemAtom.elementSymbol == 'H':
      if not chemAtom.waterExchangeable:
        chemAtoms2 = list(chemAtom.findFirstChemBond().chemAtoms)
        chemAtoms2.remove(chemAtom)
        chemAtom2 = chemAtoms2[0]
        
        if chemAtom2.elementSymbol == 'C':
          chemAtomSet  = chemAtom.chemAtomSet

          if chemAtomSet and chemAtomSet.isEquivalent:
            atom = chemAtomSet
          else:
            atom = chemAtom
          
          sideChainAtomSets.add( (atom.name, atom, chemAtom2) )
       
  sideChainAtomSets = greekSortAtomNames(list(sideChainAtomSets))
  sideChainAtomSets = [(x[1], x[2]) for x in sideChainAtomSets]   
     
  return sideChainAtomSets
  

class LinkSideChainsPopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.waiting    = False
    self.wPane1    = None
    self.wPane2    = None
    self.spectrum   = None
    self.tPeakList  = None
    self.peakList   = None
    self.rulers     = []
    self.spinSystem = None
    self.shiftList  = None
    self.assignment = None
    
    self.dataDims    = []
    self.isotopes    = []
    self.tolerances  = []
    self.tolDict     = {}
    self.knownRoots  = []
    self.tempPeaks   = []
    self.rootDims    = []
    self.nonRootDims = []
    self.specWidths  = []  
    self.fullRegion  = []
    self.diagDims    = []
    self.tempPeaks   = []
    
    BasePopup.__init__(self, parent=parent, title='Link Side Chain Resonances', **kw)

    self.color      = '#808080'

    return

    notify = self.registerNotify
    for func in ('__init__', 'delete','setName'):
      for clazz in ('ccp.nmr.Nmr.DataSource', 'ccp.nmr.Nmr.Experiment'):
        notify(self.updateSpectraAfter, clazz, func)
      
      notify(self.updateMainWindows, 'ccpnmr.Analysis.SpectrumWindow', func)

    for func in ('__init__', 'delete',):
      notify(self.updatePeakListsAfter, 'ccp.nmr.Nmr.PeakList', func)
 
    for func in ('__init__', 'delete','setValue'):
      notify(self.updateShiftsAfter, 'ccp.nmr.Nmr.Shift', func)
    
    for func in ('__init__', 'delete','setResidue',
                 'setResonances','addResonance',
                 'removeResonance', 'setName'):
      notify(self.updateAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)
 
    for func in ('setResonanceSet',):
      notify(self.updateResonancesAfter, 'ccp.nmr.Nmr.Resonance', func)
    
 
  def open(self):
  
    BasePopup.open(self)
    self.updateSpectra()
    self.updateMainWindows()
    self.updateOrthogWindows()
    self.updateAfter()


  def close(self):
  
    BasePopup.close(self)


  def destroy(self):

    self.clearTempPeaks()

    notify = self.unregisterNotify
    for func in ('__init__', 'delete','setName'):
      for clazz in ('ccp.nmr.Nmr.DataSource', 'ccp.nmr.Nmr.Experiment'):
        notify(self.updateSpectraAfter, clazz, func)
      
      notify(self.updateMainWindows, 'ccpnmr.Analysis.SpectrumWindow', func)

    for func in ('__init__', 'delete',):
      notify(self.updatePeakListsAfter, 'ccp.nmr.Nmr.PeakList', func)
 
    for func in ('__init__', 'delete','setValue'):
      notify(self.updateShiftsAfter, 'ccp.nmr.Nmr.Shift', func)
    
    for func in ('__init__', 'delete','setResidue','setResonances','addResonance','removeResonance'):
      notify(self.updateAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)
 
    for func in ('setResonanceSet',):
      notify(self.updateResonancesAfter, 'ccp.nmr.Nmr.Resonance', func)
  
    BasePopup.destroy(self)


  def body(self, guiFrame):

    guiFrame.grid_columnconfigure(0, weight=1)

    row   = 0
    frame = LabelFrame(guiFrame, text='Options')
    frame.grid(row=row, column=0, sticky='ew')
    frame.grid_columnconfigure(6, weight=1)
    
    label = Label(frame, text='Spectrum:')
    label.grid(row=0, column=0, sticky='w')
    self.spectrumPulldown = PulldownMenu(frame, callback=self.changeSpectrum,
                                         do_initial_callback=False)
    self.spectrumPulldown.grid(row=0, column=1, sticky='w')
    
    label = Label(frame, text='Peak List')
    label.grid(row=0, column=2, sticky='w')
    self.peakListPulldown = PulldownMenu(frame, callback=self.changePeakList,
                                         do_initial_callback=False)
    self.peakListPulldown.grid(row=0, column=3, sticky='w')
    
    label = Label(frame, text='Window:')
    label.grid(row=0, column=4, sticky='w')
    self.windowPulldown = PulldownMenu(frame, callback=self.changeMainWindow,
                                       do_initial_callback=False)
    self.windowPulldown.grid(row=0, column=5, sticky='w')
    
    #label = Label(frame, text='Window 2:')
    #label.grid(row=1, column=2, sticky='w')
    self.orthogWindowPulldown = PulldownMenu(frame, callback=self.changeOrthogWindow,
                                             do_initial_callback=False)
    #self.orthogWindowPulldown.grid(row=1, column=3, sticky='w')

         
    label = Label(frame, text='Min Water:')
    label.grid(row=2, column=0, sticky='w')
    self.minWaterEntry = FloatEntry(frame, width=8, text=4.76)
    self.minWaterEntry.grid(row=2, column=1, sticky='w')

    label = Label(frame, text='Max Water:')
    label.grid(row=2, column=2, sticky='w')
    self.maxWaterEntry = FloatEntry(frame, width=8, text=5.00)
    self.maxWaterEntry.grid(row=2, column=3, sticky='w')
     
    label = Label(frame, text='Diagonal width:')
    label.grid(row=2, column=4, sticky='w')
    self.diagEntry = FloatEntry(frame, width=8, text=0.15)
    self.diagEntry.grid(row=2, column=5, sticky='w') 

    frame0 = Frame(frame) 
    frame0.grid(row=3, column=0, columnspan=6, sticky='w')
    frame0.grid_columnconfigure(4, weight=1)

    label = Label(frame0, text='Tolerance 1H:')
    label.grid(row=1, column=0, sticky='w')
    self.tolHydrogenEntry = FloatEntry(frame0, width=8, text=defaultTolerances['1H'])
    self.tolHydrogenEntry.grid(row=1, column=1, sticky='w')
    label = Label(frame0, text='Tolerance 13C:')
    label.grid(row=1, column=2, sticky='w')
    self.tolCarbonEntry = FloatEntry(frame0, width=8, text=defaultTolerances['13C'])
    self.tolCarbonEntry.grid(row=1, column=3, sticky='w')
    #label = Label(frame0, text='Tolerance 15N:')
    #label.grid(row=1, column=4, sticky='w')
    self.tolNitrogenEntry = FloatEntry(frame0, width=8, text=defaultTolerances['15N'])
    #self.tolNitrogenEntry.grid(row=1, column=5, sticky='w')
       
    row  += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    frame = LabelFrame(guiFrame, text='Backbone Spin Systems')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)

    headingList = ['#','Assignment','Side chain C','Side chain H']
    justifyList = ['center','center','left','left']
    editWidgets      = [None, None, None, None,]
    editGetCallbacks = [None, None, None, None,]
    editSetCallbacks = [None, None, None, None,]
    self.spinSystemMatrix = ScrolledMatrix(frame, headingList=headingList,
                                           editSetCallbacks=editSetCallbacks,
                                           editGetCallbacks=editGetCallbacks,
                                           editWidgets=editWidgets,
                                           justifyList=justifyList,
                                           callback=self.selectSpinSystem)
    
    self.spinSystemMatrix.grid(row=0, column=0, sticky='nsew')
    
    commands = [self.prevSpinSystem,self.nextSpinSystem]
    texts = ['Previous Spin System','Next Spin System']
    self.spinSystemButtons = ButtonList(frame, texts=texts, commands=commands, expands=True)
    self.spinSystemButtons.grid(row=1, column=0, sticky='ew')
    
    row  += 1
    frame = LabelFrame(guiFrame, text='Possible Assignments')
    frame.grid(row=row, column=0, sticky='ew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)

    headingList      = ['Rank','Score','Shifts']
    self.matchMatrix = ScrolledMatrix(frame, headingList=headingList,
                                      callback=self.selectMatch)
                                 
    self.matchMatrix.grid(row=0, column=0, sticky='nsew')

    row  += 1
    frame = LabelFrame(guiFrame, text='Side Chain Assignment')
    frame.grid(row=row, column=0, sticky='ew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)

    headingList = ['Assignment','1H Shift','13C shift','Serials']
    self.resonanceMatrix = ScrolledMatrix(frame, headingList=headingList,
                                          multiSelect=True, 
                                          callback=self.selectResonance)
                                          
    self.resonanceMatrix.grid(row=0, column=0, sticky='nsew')



    row +=1
    commands = [self.pickKnownLocations,self.updateMatchLocations,
                self.commitAssignments,self.commitAllAssignments]
    texts    = ['Add Known\nPeak Locations','Search\nStrips',
                'Commit Selected\nAssignments','Commit All\nAssignments']
    self.bottomButtons = UtilityButtonList(guiFrame, helpUrl=self.help_url,
                                           commands=commands, texts=texts)
    self.bottomButtons.grid(row=row, column=0, sticky='ew')

    self.updateSpectra()
    self.updateMainWindows()
    self.updateOrthogWindows()
    self.updateAll()

  def updateSpectraAfter(self, obj=None):

   if obj.className == 'Experiment':
     spectra = obj.dataSources
   else:
     spectra = [obj,]
     
   for spectrum in spectra:
     isotopes = getSpectrumIsotopes(spectrum)
    
     numH = isotopes.count('1H')
     numC = isotopes.count('13C')
    
     if (numH + numC) > 2:
       self.after_idle(self.updateSpectra)
       return

  def selectMatch(self, obj, row, col):

    self.assignment = obj
    self.updateButtons()
    self.displayMatchLocations()
    #self.clearTempPeaks()

  def setTempPeakList(self):
  
    self.clearTempPeaks()
  
    peakList = None
    
    if self.spectrum:
      peakList = self.spectrum.findFirstPeakList(isSimulated=True, name='SideChainTemp')
    
      if not peakList:
        peakList = self.spectrum.newPeakList(isSimulated=True, name='SideChainTemp')
        analysisPeakList = getAnalysisPeakList(peakList)
        analysisPeakList.symbolColor = '#FF0000'
        analysisPeakList.textColor = '#C00000'
        analysisPeakList.symbolStyle = '+'  
        
    self.tPeakList = peakList
  
    return peakList
  
  def getPeakLists(self):
 
    peakLists = []
    
    if self.spectrum:
      for peakList in self.spectrum.sortedPeakLists():
        if not peakList.isSimulated:
          peakLists.append(peakList)
    
    return peakLists
 
  def getSpectra(self):
  
    spectra = set()
    for experiment in self.nmrProject.experiments:
      for spectrum in experiment.dataSources:
        isotopes = getSpectrumIsotopes(spectrum)
        
        numH = isotopes.count('1H')
        numC = isotopes.count('13C')
        
        if (numH + numC) > 2:
          spectra.add(spectrum)
  
    return list(spectra)
  
  def updateWindowsAfter(self, window):
  
    if self.spectrum:
      for windowPane in window.spectrumWindowPanes:
        for view in windowPane.spectrumWindowViews:
          spectrum = view.analysisSpectrum.dataSource
 
          if spectrum is self.spectrum:
            self.after_idle(self.updateMainWindows)
            return
 
  def getMainWindows(self):
  
    windows = set()
  
    if self.spectrum:
      for window in self.analysisProject.spectrumWindows:
        #if window is self.wPane2:
        #  continue
          
        for windowPane in window.spectrumWindowPanes:
          for view in windowPane.spectrumWindowViews:
            spectrum = view.analysisSpectrum.dataSource
            
            if self.spectrum is spectrum:
              windows.add(windowPane)
  
    return list(windows)
  
  def getOrthogWindows(self):
  
    windows = set()

    if self.spectrum:
      for window in self.analysisProject.spectrumWindows:
        if window is self.wPane1:
          continue
      
        for windowPane in window.spectrumWindowPanes:
          for view in windowPane.spectrumWindowViews:
            if view.isPosVisible:
              try:
                spectrum = view.analysisSpectrum.dataSource
              except:
                continue
 
              if self.spectrum is spectrum:
                windows.add(windowPane)
  
    return list(windows)
  
  def updateMainWindows(self):
  
    panes  = self.getMainWindows()
    index    = -1
    names    = []
    pane   = self.wPane1
 
    if panes:
      names = [getWindowPaneName(p) for p in panes]
 
      if pane not in panes:
        pane = panes[0]
        
      index = panes.index(pane)  
 
    else:
      pane = None
 
    if pane is not self.wPane1:
      self.wPane1 = pane
      self.displayKnownLocations()
      
 
    self.windowPulldown.setup(names, index) 
     
  def changeMainWindow(self, index, name):

    panes = self.getMainWindows()
    if panes:
      pane = panes[index]
      
      if pane is not self.wPane1:
        self.wPane1 = pane
        self.displayKnownLocations()
        
  
  def updateOrthogWindows(self):
    
    panes = self.getOrthogWindows()
    index = -1
    names = []
    pane = self.wPane2
 
    if panes:
      names = [getWindowPaneName(p) for p in panes]
 
      if pane not in panes:
        pane = panes[0]
        
      index = panes.index(pane)  
 
    else:
      pane = None
 
    if pane is not self.wPane2:
      self.wPane2 = pane
      self.displayKnownLocations()
 
    self.orthogWindowPulldown.setup(names, index) 
  
  def changeOrthogWindow(self, index, name):

    panes = self.getOthogWindows()
    if panes:
      pane = panes[index]
      
      if pane is not self.wPane2:
        self.wPane2 = pane
        self.displayKnownLocations()
 
  def updateSpectra(self):

    spectra  = self.getSpectra()
    index    = -1
    names    = []
    spectrum = self.spectrum
 
    if spectra:
      names = ['%s:%s' % (s.experiment.name,s.name) for s in spectra]
 
      if spectrum not in spectra:
        spectrum = spectra[0]
        
      index = spectra.index(spectrum)  
 
    else:
      spectrum = None
 
    if spectrum is not self.spectrum:
      self.spectrum = spectrum
      self.shiftList = spectrum.experiment.shiftList
      self.setTempPeakList()
      self.updateSpectrum()
      self.updatePeakLists()
 
    self.spectrumPulldown.setup(names, index)     
  
  def changeSpectrum(self, index, name):
  
    spectra = self.getSpectra()
    
    if spectra:
      spectrum = spectra[index]
      
      if spectrum is not self.spectrum:
        self.spectrum  = spectrum
        self.shiftList = spectrum.experiment.shiftList
        self.setTempPeakList()
        self.updateSpectrum()
        self.updatePeakLists()
  
  
  def updatePeakListsAfter(self, peakList):
  
    if peakList.dataSource is self.spectrum:
      self.updatePeakLists()
  
  
  def updatePeakLists(self):

    peakLists = self.getPeakLists()
    index    = -1
    names    = []
    peakList = self.peakList
 
    if peakLists:
      names = ['%d' % pl.serial for pl in peakLists]
 
      if peakList not in peakLists:
        peakList = peakLists[0]
        
      index = peakLists.index(peakList)  
 
    else:
      peakList = None
 
    if peakList is not self.peakList:
      self.peakList = peakList
      self.displayKnownLocations()
 
    self.peakListPulldown.setup(names, index)     
 
 
  def changePeakList(self, index, name):
  
    peakLists = self.getPeakLists()
    
    if peakLists:
      peakList = peakLists[index]
      
      if peakList is not self.peakList:
         self.peakList = peakList
         self.displayKnownLocations()
  
  def updateResonancesAfter(self, resonance):
  
    if self.waiting:
      return
    
    if self.spinSystem and (resonance.resonanceGroup is self.spinSystem):
      
      self.after_idle(self.updateResonances)
      
         
  def updateShiftsAfter(self, shift):
  
    if self.waiting:
      return

    resonance = shift.resonance
    if self.spinSystem and (resonance.resonanceGroup is self.spinSystem):
      if self.spectrum and (shift.parentList is not self.shiftList):
        return
   
      self.after_idle(self.updateResonances)
      
  
  def updateAfter(self, spinSystem=None):

    if self.waiting:
      return
      
    elif spinSystem:
      residue = spinSystem.residue
      
      if not residue:
        return

      if self.spectrum:
        if residue.chain.molSystem not in self.spectrum.experiment.molSystems:
          return     
    
    self.waiting = True
    self.after_idle(self.updateAll)
  
  def updateAll(self):

    self.updateSpinSystems()
    self.updateResonances()
    self.updateButtons()
    self.waiting = False
  
  def getSpinSystems(self):
  
    spinSystems = []
    molSystems  = []
    if self.spectrum:
      molSystems = self.spectrum.experiment.molSystems
  
    for spinSystem in self.nmrProject.resonanceGroups:
      residue    = spinSystem.residue
      resonances = spinSystem.resonances
      
      if residue and molSystems:
        if residue.chain.molSystem not in molSystems:
          continue
      
      for resonance in resonances:
        if resonance.isotopeCode == '1H':
          if resonance.covalentlyBound:
 
            bound = list(resonance.covalentlyBound)[0]
 
            if bound.isotopeCode == '13C':
 
              if residue:
                ident = '%s%5.5d' % (residue.chain.code, residue.seqCode)
              else:
                ident = '%d' % spinSystem.serial
 
              spinSystems.append((ident, spinSystem))
              break
  
    spinSystems.sort()
    
    return [s[1] for s in spinSystems]
  
  def selectSpinSystem(self, obj, row, col):
  
    if obj is not self.spinSystem:
      self.spinSystem = obj
      self.assignment = None
      self.matchMatrix.update(textMatrix=[],objectList=[])
      self.updateResonances()
      self.updateButtons()
      self.clearTempPeaks()
      self.displayKnownLocations()
  
  def updateSpinSystems(self):

    textMatrix = []
    objectList = []
    for spinSystem in self.getSpinSystems():
      
      assignStr = ''
      residue = spinSystem.residue
      if residue:
        chain = residue.chain
        
        if len(chain.molSystem.chains) > 1:
          assignStr = '%s %d%s' % (chain.code, residue.seqCode, getResidueCode(residue))
        else:
          assignStr = '%d%s' % (residue.seqCode, getResidueCode(residue))
         
      elif spinSystem.ccpCode:
        assignStr = getResidueCode(spinSystem)
      
      carbons   = spinSystem.findAllResonances(isotopeCode='13C')
      hydrogens = spinSystem.findAllResonances(isotopeCode='1H')

      carbonNames = []
      for resonance in carbons:
         if 'C' not in resonance.assignNames:
           carbonNames.append(getResonanceName(resonance))
           
      carbonNames = greekSortAtomNames(carbonNames)

      hydrogenNames = []
      for resonance in hydrogens:
         if 'H' not in resonance.assignNames:
           hydrogenNames.append(getResonanceName(resonance))
           
      hydrogenNames = greekSortAtomNames(hydrogenNames)

      datum = [spinSystem.serial,
               assignStr,
               ' '.join(carbonNames),
               ' '.join(hydrogenNames),
               ]
 
      textMatrix.append(datum)
      objectList.append(spinSystem)
     
    self.spinSystemMatrix.update(textMatrix=textMatrix,objectList=objectList)
    
  def selectResonance(self, obj, row, col):
  
    #self.resonanceAssignment = obj
    pass
  
  def updateResonances(self):

    return

    textMatrix = []
    objectList = []
    colorMatrix = []
    
    if self.spinSystem:
    
      resonanceDict = {}
      prochiralDict = {}
      for resonance in self.spinSystem.resonances:
        resonanceSet = resonance.resonanceSet
        
        if resonanceSet:
          atomSets = list(resonanceSet.atomSets)
          
          index = 0
          if len(atomSets) > 1:
            resonances = list(resonanceSet.resonances)
            if len(resonances) > 1:
              index = resonances.index(resonance)
            #  if index == 0:
            #    atomSets.reverse()
            #  
            #prochiralDict[resonance] = '/'.join([as.name for as in atomSets])  
          
          chemAtom = atomSets[index].findFirstAtom().chemAtom
          chemAtomSet = chemAtom.chemAtomSet
          
          if chemAtomSet and chemAtomSet.isEquivalent:
            resonanceDict[chemAtomSet] = resonance
          else:
            resonanceDict[chemAtom] = resonance
      
      assignDict = {}
      if self.assignment:
        for atomSets, ppms, peaks in self.assignment:
          atomH, atomC = atomSets
          ppmH, ppmC   = ppms
          
          assignDict[atomH] = ppmH
          assignDict[atomC] = ppmC
            
      for atomH, atomC in getSideChainAtomSets(self.spinSystem):
 
        colors = [None] * 4
 
        resonanceH = resonanceDict.get(atomH)
        resonanceC = resonanceDict.get(atomC)
 
        if prochiralDict.get(resonanceH):
          nameH = prochiralDict[resonanceH] or '?'
        else:
          nameH = atomH.name or '?'

        if prochiralDict.get(resonanceC):
          nameC = prochiralDict[resonanceC] or '?'
        else:
          nameC = atomC.name or '?'
 
        serials = []
        
        shiftH = None
        if resonanceH:
          shift = resonanceH.findFirstShift(parentList=self.shiftList)
          serials.append('[%d]' % resonanceH.serial)
        
          if shift:
            shiftH = shift.value
        
        else:
          shiftH = assignDict.get(atomH)
          #serials.append(None)
          
          if shiftH is not None:
            colors[1] = '#B0FFB0'
        
        shiftC = None
        if resonanceC:
          shift = resonanceC.findFirstShift(parentList=self.shiftList)
          serials.append('[%d]' % resonanceC.serial)
        
          if shift:
            shiftC = shift.value

        else:
          shiftC = assignDict.get(atomC)
          #serials.append(None)

          if shiftC is not None:
            colors[2] = '#B0FFB0'

        datum = ['%s %s' % (nameH, nameC),
                 shiftH,
                 shiftC,
                 ','.join(serials)
                ]
 
        textMatrix.append(datum)
        objectList.append((atomH, atomC, resonanceH, resonanceC, shiftH, shiftC))
        colorMatrix.append(colors)

    self.resonanceMatrix.update(textMatrix=textMatrix,objectList=objectList, colorMatrix=colorMatrix)
  
  def getTolerances(self):
  
    self.tolerances = []
    for isotopes in self.isotopes:
      isotope = isotopes.copy().pop()
    
      tol = 0.1
      if isotope == '13C':
        tol = self.tolCarbonEntry.get() or defaultTolerances[isotope]
      elif isotope == '1H':
        tol = self.tolHydrogenEntry.get() or defaultTolerances[isotope]
      elif isotope == '15N':
        tol = self.tolNitrogenEntry.get() or defaultTolerances[isotope]
    
      self.tolDict[isotope] = tol
      self.tolerances.append(tol)
 
  def updateButtons(self):
    
    buttons1 = self.spinSystemButtons.buttons
    if self.spinSystemMatrix.objectList:
      for button in buttons1:
        button.enable()
    else:
      for button in buttons1:
        button.disable()
    
    buttons2 = self.bottomButtons.buttons[:4]
    if self.spinSystem:
      for button in buttons2[:2]:
        button.enable()
      
      if self.assignment:
        buttons2[2].enable()
        buttons2[3].enable()
      else:
        buttons2[2].disable()  
        buttons2[3].enable()
        
    else:
      for button in buttons2:
        button.disable()
          

  def updateSpectrum(self):
  
    dataDims = self.spectrum.sortedDataDims()
    
    self.dataDims    = dataDims
    self.isotopes    = []
    self.rootDims    = []
    self.nonRootDims = []
    self.specWidths  = []  
    self.fullRegion  = []
    self.diagDims    = []
    
    boundDataDims  = getOnebondDataDims(self.spectrum)
    if not boundDataDims:
      showWarning('Invalid experiment type','Selected spectrum has no bound dimensions')
      return
    
    for dataDim in dataDims:
      self.isotopes.append(getDataDimIsotopes(dataDim))
    
    checkedDims = []
    for dim1, dim2 in boundDataDims:
      if '1H' not in self.isotopes[dim1.dim-1]:
        if '1H' not in self.isotopes[dim2.dim-1]:
          continue
          
      checkedDims.append((dim1, dim2))

    isRoot = {}
    for dim1, dim2 in checkedDims:

      if len(checkedDims) == 1:
        isRoot[dim1] = True
        isRoot[dim2] = True
          
      else: # 4D
        for dim1, dim2 in boundDataDims:
          if (dim1.expDimRef.isAcquisition) or (dim2.expDimRef.isAcquisition):
            isRoot[dim1] = True
            isRoot[dim2] = True
            break
    
    for dataDim in dataDims:
      i = dataDim.dim-1
      dataDimRef  = getPrimaryDataDimRef(dataDim)
      expDimRef   = dataDimRef.expDimRef
      minPointPpm = pnt2ppm(dataDim.numPoints,dataDimRef)
      maxPointPpm = pnt2ppm(0,dataDimRef)
      
      self.specWidths.append(maxPointPpm-minPointPpm)
      
      valueMin   = expDimRef.minAliasedFreq
      if valueMin is None:
        valueMin   = minPointPpm

      valueMax   = expDimRef.maxAliasedFreq
      if valueMax is None:
        valueMax   = maxPointPpm

      self.fullRegion.append([valueMin,valueMax])    
     
      if isRoot.get(dataDim):
        self.rootDims.append(i)      
      else:
        self.nonRootDims.append(i)

    for dim1 in self.nonRootDims:
      isotopes = self.isotopes[dim1]
      
      for dim2 in self.rootDims:
        
        if isotopes == self.isotopes[dim2]:
          self.diagDims.append((dim2, dim1))
          break
       

  def displayKnownLocations(self):
    
    return
    
    self.clearRulers()
    self.getTolerances()
     
    if self.wPane1 and self.spinSystem and self.dataDims:
    
      for view in self.wPane1.spectrumWindowViews:
        try:
          if view.analysisSpectrum.dataSource is self.spectrum:
            if not view.isPosVisible:
              view.isPosVisible = True
        except:
          continue
      
      ccpCode  = self.spinSystem.ccpCode
      project  = self.spinSystem.root 
    
      shifts, orthoShifts = getSpinSystemWindowShifts([self.spinSystem,], self.wPane1,
                                                      shiftList=self.shiftList)
            
      stripAxis = self.wPane1.spectrumWindow.stripAxis
      sortList = []
      for shift in shifts:
        resonance = shift[stripAxis].resonance
        if resonance.assignNames:
          letter = resonance.assignNames[0][1]
          order = sideChainSortOrder.get(letter, 7)
        else:
          order = 7
       
        sortList.append((order,shift))
      sortList.sort()
       
      orthoPositions = None
      #if splitIntoCells and orthoShifts:
      #  orthoPositions = [shift.value for shift in orthoShifts]
      #  orthoPositions.sort()
      
      shifts          = [x[1] for x in sortList]
      axisMapping     = getDataDimAxisMapping(self.spectrum, self.wPane1)
 
      for shiftDict in shifts:
        for key in shiftDict:
          shift = shiftDict[key]
          shiftDict[key] = shift.value
 
      # Make strips
      displayStrips(self.parent, shifts,
                    orthoPositions=orthoPositions,
                    windowPane=self.wPane1)
      self.update_idletasks()
      displayStrips(self.parent, shifts,
                    orthoPositions=orthoPositions,
                    windowPane=self.wPane1)

  def getKnownRoots(self, spinSystem, sideChainAtomSets):
  
    rootResonances = []
    knownRoots     = []
    
    for resonance in spinSystem.resonances:
      if resonance.isotopeCode == '1H':
        resonanceSet = resonance.resonanceSet
        
        if resonanceSet:
          bound = resonance.findFirstCovalentlyBound()
   
          if bound and bound.isotopeCode == '13C' and bound.resonanceSet:
            rootResonances.append((resonance, bound)) 
   
    dims = [None, None]   
    dataDims = self.dataDims
    for i in self.rootDims:
      dataDim = dataDims[i]
      if '1H' in self.isotopes[i]:
        dims[0] = dataDim
      else: 
        dims[1] = dataDim
    
    for resonances in rootResonances:
      knownLocation = {}
      resonanceDict = {}
      atomSite = []
      
      j = 0
      for resonance in resonances:
        shift = resonance.findFirstShift(parentList=self.shiftList)
        
        if not shift:
          break
      
        dataDim = dims[j]
        dim     = dataDim.dim-1
        resonanceSet = resonance.resonanceSet
        index = list(resonanceSet.resonances).index(resonance)
        atom  = list(resonanceSet.atomSets)[index].findFirstAtom()
 
        chemAtom = atom.chemAtom
        if chemAtom.chemAtomSet and chemAtom.chemAtomSet.isEquivalent:
          atomKey = chemAtom.chemAtomSet
        else:
          atomKey = chemAtom     
        
        atomSite.append(atomKey)
        
        ppms = []
        for contrib in resonance.peakDimContribs:
          peakDim = contrib.peakDim
          if peakDim.dataDimRef.dataDim is dataDim:
            ppms.append(peakDim.value)
           
        if ppms:
          ppm = sum(ppms)/float(len(ppms))
        else:
          ppm = shift.value
        
        error = self.tolerances[dim]
        knownLocation[dim] = (ppm, ppm-error, ppm+error)
        resonanceDict[dim] = resonance
        
        j += 1
        
      else:
        knownRoots.append((knownLocation, resonanceDict, tuple(atomSite) ))
     
    return knownRoots

  def pickKnownLocations(self, sideChainAtomSets=None, peakList=None):

    self.clearRulers()
    self.getTolerances()
    self.clearTempPeaks()
    self.assignment = None
    self.matchMatrix.update(textMatrix=[],objectList=[])
     
    if not self.spinSystem:
      return
    
    if not sideChainAtomSets:
      sideChainAtomSets = getSideChainAtomSets(self.spinSystem)
 
    if not peakList:
      peakList = self.peakList
    
    rootDims   = self.rootDims
    knownRoots = self.getKnownRoots(self.spinSystem, sideChainAtomSets)
    if not knownRoots:
      return
     
    atomSitePeaks = {}
    
    if self.spectrum and self.spinSystem:
    
      minWater  = self.minWaterEntry.get() or 0.0
      maxWater  = self.maxWaterEntry.get() or 0.0
      diagTol   = self.diagEntry.get()     or 0.0
      ccpCode   = self.spinSystem.ccpCode
      numDim    = self.spectrum.numDim
      shiftList = self.shiftList
      shiftDict = {}
      resonDict = {}
      
      nonRootDims = self.nonRootDims
      diagDims    = self.diagDims
      isotopes    = self.isotopes
      
      fixedLocations = {}
      partLocations  = {}
      
      fullRegions = []
      
      unassignedResonances = {} 
         
      for resonance in self.spinSystem.resonances:
        shift = resonance.findFirstShift(parentList=shiftList)
        
        if not shift:
          continue
      
        values = []
        for contrib in resonance.peakDimContribs:
          peakDim = contrib.peakDim
          
          if peakDim.peak.peakList is self.peakList:
            values.append(peakDim.value)
      
        if values:
          shiftValue = sum(values)/float(len(values))
        else:
          shiftValue = shift.value   
      
        resonanceSet = resonance.resonanceSet
        if resonanceSet:
          index = list(resonanceSet.resonances).index(resonance)
          atom  = list(resonanceSet.atomSets)[index].findFirstAtom()
          
          chemAtom = atom.chemAtom
          if chemAtom.chemAtomSet and chemAtom.chemAtomSet.isEquivalent:
            atomKey = chemAtom.chemAtomSet
          else:
            atomKey = chemAtom
    
          shiftDict[atomKey] = shiftValue
          resonDict[atomKey] = resonanceSet.resonances
          
        elif resonance.isotopeCode == '1H':
          uShifts = [shiftValue]
          bound   = resonance.findFirstCovalentlyBound()
          
          if bound:
            shiftX  = bound.findFirstShift(parentList=shiftList)

            if shiftX:
              valuesX = []
 
              for contrib in bound.peakDimContribs:
                peakDim = contrib.peakDim
 
                if peakDim.peak.peakList is self.peakList:
                  valuesX.append(peakDim.value)
 
              if valuesX:
                shiftValueX = sum(valuesX)/float(len(valuesX))
              else:
                shiftValueX = shiftX.value
            
              uShifts.append(shiftValueX) 
        
          unassignedResonances[resonance] = uShifts
    
      for atomKey in sideChainAtomSets:

        atomSitePeaks[atomKey] = []
        ppms = []
        for i in nonRootDims:
          if '1H'in self.isotopes[i]:
            ppms.append(shiftDict.get(atomKey[0]))
          else:
            ppms.append(shiftDict.get(atomKey[1]))
        
        if None in ppms:
          partLocations[atomKey] = ppms
 
        elif ppms:
          fixedLocations[atomKey] = ppms
    
      for pickLocation, rootResonances, knownAtomSite in knownRoots:
            
        for atomKey in sideChainAtomSets:
          fixedLocation = fixedLocations.get(atomKey)
          
          if fixedLocation is None:
            # Try to pick something near BMRB shift
            # - pick in region about root and BMRB optimum
            # - very generous in indirect, but remove all but best
            region = [None] * numDim
            for i in rootDims:
              region[i] = pickLocation[i][1:]
                          
            for i in nonRootDims:
              if '1H' in isotopes[i]:
                atom = atomKey[0]
              else:
                atom = atomKey[1]

              # atom is actually a chemAtom or equivalent chemAtomSet
              
              ppmRange = None
              nmrRef   = getChemAtomNmrRef(self.project, atom.name, ccpCode, molType='protein')
              
              if nmrRef:
                distribution  = nmrRef.distribution
                refPoint      = nmrRef.refPoint
                refValue      = nmrRef.refValue
                valuePerPoint = nmrRef.valuePerPoint
                n = len(distribution)
                
                minPt  = 0
                maxPt  = n-1
                
                for v in distribution:
                  if v < 0.0001:
                    minPt += 1
                  else:
                    break  
                
                for v in distribution[::-1]:
                  if v < 0.0001:
                    maxPt -= 1
                  else:
                    break  
                 
                maxPpm = ((maxPt - refPoint) * valuePerPoint) + refValue + 0.5
                minPpm = ((minPt - refPoint) * valuePerPoint) + refValue - 0.5
                region[i] = (minPpm, maxPpm)
            
            if None not in region:              
              # chop out water stripe from pick region
              regions = []
              for dataDim in self.dataDims:
                if dataDim.expDim.isAcquisition:
                  dim = dataDim.dim-1
                  if '1H' in isotopes[dim]:
                    excludeRegion(regions, region, dim, minWater, maxWater)
                            
              # chop out diagonal from pick region
              for dim1, dim2 in diagDims: # root, nonRoot
                ppm    = pickLocation[dim1][0]
                minPpm = ppm-diagTol
                maxPpm = ppm+diagTol
                
                regions2 = []
                for region2 in regions:
                  excludeRegion(regions2, region2, dim2, minPpm, maxPpm)
 
                regions = regions2
              
              doneUnassigned = {}
              peaks0 = []
              # Look for existing peaks in proper list assigned to the root
              for peak in searchPeaks([self.peakList,], region2):
                peakDims = peak.sortedPeakDims()
                
                for i in rootDims:
                  peakDim = peakDims[i]
                  resonances = [contrib.resonance for contrib in peakDim.peakDimContribs]
                  if rootResonances[peakDim.dim-1] not in resonances:
                    break
                    
                else:
                  for k in nonRootDims:
                    if '1H' in isotopes[k]:
                      for resonance in [contrib.resonance for contrib in peakDims[k].peakDimContribs]:
                        if unassignedResonances.get(resonance):
                          doneUnassigned[resonance] = True
                 
                  peaks0.append(peak)
              
              # Pick peaks at positions where we have sufficient unassigned resonances
              # - these must be bound pairs for 4Ds
              for resonance in unassignedResonances:
                if not doneUnassigned.get(resonance):
                  ppms = unassignedResonances[resonance]
                  
                  if len(ppms) >= len(nonRootDims):
                    position = [None] * numDim
 
                    for i in rootDims:
                      position[i] = pickLocation[i][0] # value, min, max
 
                    for i in nonRootDims:
                      if '1H' in isotopes[i]:
                        ppm = ppms[0]
                      else:  
                        ppm = ppms[1]
                      position[i] = ppm
                    
                    region3 = [(p-0.001,p+0.001) for p in position]
                    existingPeaks = searchPeaks([peakList,], region3)
                    
                    if existingPeaks:
                      peak = existingPeaks[0]
                    else:  
                      peak = pickPeak(peakList, position, unit='ppm', doFit=True)
                      #snapPeak(peak)
                    
                    peaks0.append(peak)
                  
              for region2 in regions:
                peaks0.extend( searchPeaks([peakList,], region2) )
                peaks0.extend( findPeaks(peakList, region2) )
              
              peakScores = []
              for peak in peaks0:
                peakDims = peak.sortedPeakDims()
              
                score = 0.0
                
                for i in nonRootDims:
                  if '1H' in isotopes[i]:
                    atom = atomKey[0]
                  else:
                    atom = atomKey[1]
                 
                  ppm = peakDims[i].value
                  sc  = getAtomProbability(ccpCode, atom.name, ppm)
                  score += sc*sc
                  
                peakScores.append((score, peak))

              if peakScores:
                peakScores.sort()
                for score, peak in peakScores[:-1]:
                  if peak.peakList is not self.peakList:
                    peak.delete()
 
                peak = peakScores[-1][1]
                peakDims = peak.sortedPeakDims()
                atomSitePeaks[atomKey].append(peak)

                for i in rootDims:
                  assignResToDim(peakDims[i], rootResonances[i])
 
              del peakScores
            
          else:
            peak = None
            assignment = [None] * numDim
            
            for i in rootDims:
              assignment[i] = (rootResonances[i],)
          
            j = 0
            for i in nonRootDims:
              if '1H' in isotopes[i]:
                atom = atomKey[0]
              else:
                atom = atomKey[1] 
                            
              assignment[i] = resonDict.get(atom)  
              j += 1
          
          
            # Try to find existing peak of required assignment
            for resonance in rootResonances.values():
              for contrib in resonance.peakDimContribs:
                peak2 = contrib.peakDim.peak
                
                if peak2.peakList is not peakList:
                  continue
                 
                assignment0 = []
                for peakDim in peak2.sortedPeakDims():
                  resonances = [contrib.resonance for contrib in peakDim.peakDimContribs]
                  
                  if resonances:
                    assignment0.append(tuple(resonances))
    
                if assignment0 == assignment:
                  peak = peak2
                  break
  
              else:
                continue
                
              break

            
            if not peak:
              position = [None] * numDim
 
              for i in rootDims:
                position[i] = pickLocation[i][0] # value, min, max
 
              j = 0
              for i in nonRootDims:
                position[i] = fixedLocation[j]
                j += 1
 
              region3 = [(p-0.001,p+0.001) for p in position]
              existingPeaks = searchPeaks([peakList,], region3)
 
              if existingPeaks:
                peak = existingPeaks[0]
              else:
                peak = pickPeak(peakList, position, unit='ppm', doFit=True)
                #snapPeak(peak)
 
              peakDims = peak.sortedPeakDims()
              for i in rootDims:
                assignResToDim(peakDims[i], rootResonances[i])

              for i in nonRootDims:
                if '1H' in isotopes[i]:
                  atom = atomKey[0]
                else:
                  atom = atomKey[1]
  
                nonRootResonances = resonDict.get(atom)
                if nonRootResonances:
                  peakDim = peakDims[i]
 
                  for resonance in nonRootResonances:
                    assignResToDim(peakDim, resonance)

            
            atomSitePeaks[atomKey].append(peak)  
               
    return atomSitePeaks

  def nextSpinSystem(self):
    
    spinSystems = self.spinSystemsMatrix.objectList
    
    index = self.spinSystemsMatrix.currentCell[0]
    if index is None:
      index = -1
 
    index += 1  
    if index >= len(spinSystems):
      index = 0  

    self.peakMatrix.selectObject(spinSystems[index])

  def prevSpinSystem(self):
    
    spinSystems = self.spinSystemsMatrix.objectList
    
    index = self.spinSystemsMatrix.currentCell[0]
    if index is None:
      index = 1
      
    index -= 1  
    if index < 0:
      index = len(spinSystems)-1  
 
    self.peakMatrix.selectObject(spinSystems[index])

  def mergeIndirectValues(self, atomSitePeaks):
    
    # list of list of indirect values (could be 4D)
    
    # TBD: Merge only peaks of equivalent atom sites
    
    # positions are non-rootDim ppms

    
    mergedPositions = {}
    tolerances2 = [self.tolerances[dim] for dim in self.nonRootDims]  
    
    for atomSite in atomSitePeaks:
      peaks = atomSitePeaks[atomSite]
      positions = []
      
      for peak in peaks:
        peakDims = peak.sortedPeakDims()
        position = []
        
        for i in self.nonRootDims:
          position.append(peakDims[i].value)
      
        positions.append(position)
           
      M = len(self.nonRootDims)
      N = len(positions)
      if N == 1:
        mergedPositions[atomSite] = positions[0]
 
      elif N > 1:
        # Take the mean position of the largest cluster, or the most likely
      
        clusters = [[p,] for p in positions]
        
        for i in range(N-1):
          clusterA = clusters[i]
          
          for j in range(i+1,N):
            clusterB = clusters[j]
            
            merge = True
            for posA in clusterA:
              for posB in clusterB:
                for k in range(M):
                  delta = abs(posA[k]-posB[k])
 
                  if delta > tolerances2[k]/10.0:
                    break
                else:
                  break
              else:
                continue
              break  
            else:
              merge = False
      
            if merge:
              clusters[i] = clusterA + clusterB
              clusters[j] = []
      
        sortClusters = [(len(c), c) for c in clusters]
        sortClusters.sort()
        
        mean = [0.0] * M
        num  = 0.0
        for position in sortClusters[-1][1]:
          for k in range(M):
            mean[k] += position[k]
          num += 1.0
          
        for k in range(M):
          mean[k] /= num
        
        mergedPositions[atomSite] = mean
             
    return mergedPositions
    
  def updateMatchLocations(self):

    self.clearRulers()
    self.getTolerances()
    self.clearTempPeaks()
    
    if self.spinSystem and self.spectrum:
      
      # # # # # # # # SETUP CONSTANTS # # # # # # # #
      ccpCode     = self.spinSystem.ccpCode
      project     = self.spinSystem.root
      rootDims    = self.rootDims
      nonRootDims = self.nonRootDims
      isotopes    = self.isotopes
      tolerances  = self.tolerances
      shiftList   = self.shiftList
       
      # Get the assignable locations along the side chain
      sideChainAtomSets = getSideChainAtomSets(self.spinSystem)
      
      # Get the peaks that are or could be responsible for each atom site 
      atomSitePeaks  = self.pickKnownLocations(sideChainAtomSets, peakList=self.tPeakList)
      for atomSite in atomSitePeaks:
        self.tempPeaks.extend(atomSitePeaks[atomSite])
 
      # Merge results from all roots into one list of indirect positions
      # for all atom sites
      indirectPositions = self.mergeIndirectValues(atomSitePeaks)
  
      # Make window rulers of indirect search locations
      if self.wPane1:
        axisMapping = getDataDimAxisMapping(self.spectrum, self.wPane1)
        dimMapping = {}
        for axis in axisMapping:
          dimMapping[axisMapping[axis]] = axis
       
        panelTypes = []
        for dim in nonRootDims:
          axis = dimMapping[self.dataDims[dim]]
          axisPanel = self.wPane1.findFirstAxisPanel(label=axis)
          panelTypes.append(axisPanel.panelType)
       
        for atomSite in indirectPositions:
          ppms = indirectPositions[atomSite]
          
          for i in range(len(nonRootDims)):
            ruler = createRuler(ppms[i], panelTypes[i], dashLength=2, gapLength=2,
                                color=self.color, remove=False)
            self.rulers.append(ruler)
  
      # Extend spin system by picking peaks in potential
      # additional carbon locations
      # Uses indirect H locations to search through all 13C planes
      newRootPeakLocations = self.searchCarbonPlanes(indirectPositions) 
      
      # Cluster peaks at similar 13C planes to provide
      # estimate of potential root locations
      rawRoots = self.clusterRootLocations(newRootPeakLocations)
      for position, peaks in rawRoots:
        self.tempPeaks.extend(peaks)
      
      # Expand roots to include aliased positions
      # and trim according to atom types 
      newPositions = self.unaliasRootLocations(rawRoots, indirectPositions, sideChainAtomSets)    
            
      # # # # # FILL IN MATRICES FOR SIDE CHAIN ROUTE SEARCH # # # # # 
      fixedPositions = {}
      for resonance in self.spinSystem.resonances:
        shift = resonance.findFirstShift(parentList=shiftList)
        resonanceSet = resonance.resonanceSet
        
        if resonanceSet and shift:
          atomSets = resonanceSet.sortedAtomSets()
          index = 0
          
          if len(atomSets) > 1:
            resonances = list(resonanceSet.resonances)
            if len(resonances) > 1:
              index = resonances.index(resonance)
      
          chemAtom = atomSets[index].findFirstAtom().chemAtom
          chemAtomSet = chemAtom.chemAtomSet
          
          if chemAtomSet and chemAtomSet.isEquivalent:
            fixedPositions[chemAtomSet] = shift.value
          else:
            fixedPositions[chemAtom] = shift.value

      
      linkThreshold  = 0.00001
      matchThreshold = 0.00001
      
      I = len(newPositions)
      J = len(sideChainAtomSets)
      R = range(len(rootDims))
      
      objects           = range(I)
      locations         = range(J)
      voidLocations     = {}
      linkScoreMatrix   = [[0.0]*I for i in objects]
      matchScoreMatrix  = [[0.0]*J for i in objects]
      contextVoidMatrix = [[None]*J for i in objects]
      atomConnMatrix    = [[0.0]*J for i in locations]
      
      for i in range(J-1):
        iAtomH, iAtomC = sideChainAtomSets[i]
        bonds = list(iAtomC.chemBonds)
        
        for j in range(i+1, J):
          jAtomH, jAtomC = sideChainAtomSets[j]
          if iAtomH is jAtomH:
            continue
          
          atomConnMatrix[i][j] = 0.7
          atomConnMatrix[j][i] = 0.7
        
          if iAtomC is jAtomC: # Prochirals share carbon
            atomConnMatrix[i][j] = 1.0
            atomConnMatrix[j][i] = 1.0
            
          else:
            for bond in bonds:
              if jAtomC in bond.chemAtoms:
                atomConnMatrix[i][j] = 1.0
                atomConnMatrix[j][i] = 1.0
      
      for i, (rootPos, nonRootPositions, peaks, score) in enumerate(newPositions):
        #score /= bestMatchListScore
        
        #print 'PK', i, ['%.3f' % x for x in rootPos], '%.4f' % score

        for j, atomSets in enumerate(sideChainAtomSets):
          # Work with atomSets or chemAtomSets?
          
          matchScore   = 0.0  # From BMRB match and peak closeness
          contextVoids = contextVoidMatrix[i][j] or []

          fixedPosition = fixedPositions.get(atomSets)
          
          for k in R:
            fixedPosition = fixedPositions.get(atomSets[k])
            
            if fixedPosition:
              delta = abs(rootPos[k] - fixedPosition)
              
              if delta > tolerances[rootDims[k]]:
                matchScore = 0.0 # Impossible
                break

          if not matchScore:
            for k in R:
              ppm         = rootPos[k]
              atomSet     = atomSets[k]
              atomSetName = atomSet.name
              element     = atomSetName[0]
              score       = getAtomProbability(ccpCode, atomSetName, ppm)
 
              if score < matchThreshold: # e.g. if one of an H-C pair is impossible, both no good
                matchScore = 0.0
                break
              else:
                matchScore += score
 
              if j > 0:
                atomSets2 = sideChainAtomSets[j-1]
 
                if atomSet in atomSets2: # second of a prochiral
 
                  for k, (rootPos2, nonRootPositions2, peaks2, score2) in enumerate(newPositions):
                    if k == i:
                      continue
                      
                    for l, dim in enumerate(rootDims):
                      if '13C' in isotopes[dim]:
                        diff = abs(rootPos[l]-rootPos2[l])
                        
                        if diff > tolerances[dim]:
                          contextVoids.append(k)

          contextVoidMatrix[i][j] = contextVoids
          matchScoreMatrix[i][j]  = matchScore
          
          # - With object A at location I this object B is not
          #   allowed to follow at location I+1
          # = With rootPos A at atomSite I the rootPos B is not
          #   allowed to follow at atomSite I+1 unless its 13C
          #   shift is similar
          
        for k, (rootPos2, nonRootPositions2, score2, peaks2) in enumerate(newPositions):
          if k == i:
            continue
            
          #print 'LINK', ','.join(['%.3f' % x for x in rootPos]), ','.join(['%.3f' % x for x in rootPos2])
          #print 'NR', 

          linkScore = 0.0      # From shared peak positions, or NOEs...
 
          for nonRootPos in nonRootPositions:
          
            closest = None
            for nonRootPos2 in nonRootPositions2:
              delta2 = 0.0
 
              for l, dim in enumerate(nonRootDims):
                diff = (nonRootPos[l]-nonRootPos2[l])/tolerances[dim]
                delta2 += diff*diff
 
              if (closest is None) or (delta2 < closest[0]):
                closest = (delta2,  nonRootPos2)
 
            if closest:
              delta2,  nonRootPos2 = closest
              e = delta2 /twoGaus
              q = exp(-e)/gausRoot2pi
              linkScore += q
            
              #print '%.3f' % nonRootPos[0], '%.3f' % closest[1][0],
          
          #print ''
          #print linkScore
          
          linkScoreMatrix[i][k] = linkScore

      
      print I, J # Check size of problem

      # # # # # # # # # FIND AND SCORE WHOLE SIDE CHAIN ROUTES # # # # # # # # # #
            
      """
      Matrix of link scores - Matches between indirect carbon positions for each CH
                            - Size i,i is CH roots vs CH roots
 
      Matrix of match scores - Matches between CH chem shifts and BMRB distributions
                             - Size i,j is CH roots vs side chain atom sets
 
      Matrix of contextual - For each CH chem shift and side chain atom set
      voids                - When a prochiral CH is present, then the next in the CH2
                             group must follow, identified by same carbon shift
                             
                             
      C is a transformation in resonance space that 
      1) Gets you from this spin system to the i-1 spin system
      2) Gets you from this CH resonance to the whole side chain of CH resonances
      
      R is a transformation in atom space that does the equivalent of C
      1) Gets you from this residue to the i-1 residue
      2) Gets you from this CH atom to the whole side chain of CH atoms
      
      A is the assignment matrix to map from resonance space to atom space
      The initial values in A are based on chemical shifts
      1) 
      2) 
      
      """
      import Tkinter
      from memops.gui.ScrolledDensityMatrix import ScrolledDensityMatrix
      
      root = Tkinter.Tk()
      root.geometry('800x850')
      root.grid_columnconfigure(0, weight=1)
      root.grid_columnconfigure(1, weight=1)
      root.grid_columnconfigure(2, weight=1)
      root.grid_rowconfigure(0, weight=1)
      root.grid_rowconfigure(1, weight=1)
      root.grid_rowconfigure(2, weight=1)
      
      C = matrix(linkScoreMatrix)
      
      #A = 0.01 * ones((I,J+1))
      #A[:I,:J] = matrix(matchScoreMatrix)
      A = matrix(matchScoreMatrix)
      
      #R = ones((J+1,J+1)) - identity(J+1)
      #R[:J,:J] = matrix(atomConnMatrix)
      R =  matrix(atomConnMatrix)
      
      C = matrix(C)/C.sum()
      A = matrix(A)/A.sum()
      R = matrix(R)/R.sum()
      
      
      null = ['-',] 
      atomLabels = [a[0].name for a in sideChainAtomSets]  + null
      numLabels = ['%d' % (i+1) for i in range(I)]
      
      graph = ScrolledDensityMatrix(root, matrix=C.tolist(), boxSize=12, maxVal=None,
                                    xLabels=numLabels, yLabels=numLabels, borderWidth=1,
                                    borderColor='grey', zoom=1.0,
                                    title='Links')
      graph.grid(row=0, column=0, sticky='nsew')
      
      graph = ScrolledDensityMatrix(root, matrix=A.tolist(), boxSize=12, maxVal=None,
                                    xLabels=numLabels, yLabels=atomLabels, borderWidth=1,
                                    borderColor='grey', zoom=1.0,
                                    title='Matches')
      graph.grid(row=0, column=1, sticky='nsew')
     
      graph = ScrolledDensityMatrix(root, matrix=R.tolist(), boxSize=12, maxVal=None,
                                    xLabels=atomLabels, yLabels=atomLabels, borderWidth=1,
                                    borderColor='grey', zoom=1.0,
                                    title='Atom Connections')
      graph.grid(row=0, column=2, sticky='nsew')
      
      
      A1 = C*A
      graph = ScrolledDensityMatrix(root, matrix=A1.tolist(), boxSize=12, maxVal=None,
                                    xLabels=numLabels, yLabels=atomLabels, borderWidth=1,
                                    borderColor='grey', zoom=1.0,
                                    title='Iteration')
      graph.grid(row=1, column=0, sticky='nsew')
      
      
      A2 = A*R
      graph = ScrolledDensityMatrix(root, matrix=A2.tolist(), boxSize=12, maxVal=None,
                                    xLabels=numLabels, yLabels=atomLabels, borderWidth=1,
                                    borderColor='grey', zoom=1.0,
                                    title='Iteration')
      graph.grid(row=1, column=1, sticky='nsew')
 
      
      A3 = iteration(C,A,R)
      graph = ScrolledDensityMatrix(root, matrix=A3.tolist(), boxSize=12, maxVal=None,
                                    xLabels=numLabels, yLabels=atomLabels, borderWidth=1,
                                    borderColor='grey', zoom=1.0,
                                    title='Iteration')
      graph.grid(row=1, column=2, sticky='nsew')

      
      A4 = iteration(C,A3,R)
      graph = ScrolledDensityMatrix(root, matrix=A4.tolist(), boxSize=12, maxVal=None,
                                    xLabels=numLabels, yLabels=atomLabels, borderWidth=1,
                                    borderColor='grey', zoom=1.0,
                                    title='Iteration')
      graph.grid(row=2, column=0, sticky='nsew')
      

      A5 = iteration(C,A4,R)
      graph = ScrolledDensityMatrix(root, matrix=A5.tolist(), boxSize=12, maxVal=None,
                                    xLabels=numLabels, yLabels=atomLabels, borderWidth=1,
                                    borderColor='grey', zoom=1.0,
                                    title='Iteration')
      graph.grid(row=2, column=1, sticky='nsew')
 
      A6 = iteration(C,A5,R)
      graph = ScrolledDensityMatrix(root, matrix=A6.tolist(), boxSize=12, maxVal=None,
                                    xLabels=numLabels, yLabels=atomLabels, borderWidth=1,
                                    borderColor='grey', zoom=1.0,
                                    title='Iteration')
      graph.grid(row=2, column=2, sticky='nsew')
      
      return
      
      assignments = dynamicBranchPrune(I, J,
                                       linkScoreMatrix, linkThreshold,
                                       matchScoreMatrix, matchThreshold,
                                       contextVoidMatrix,
                                       isLinearSystem=False)


      #print assignments
      #print locations
      #print sideChainAtomSets
      #print newPositions

      textMatrix = []
      objectList = []

      c = 0
      for score, assignment in assignments:
        c += 1      
        assignTexts = []
        assignmentObj = []
        
        for i, j in enumerate(assignment):
          atomSets = sideChainAtomSets[i]
          ppms, null, peaks, sc = newPositions[j]
          ppm = ','.join(['%.2f' % x for x in ppms])

          assignText =  '%s:%s' % (atomSets[0].name, ppm)
          assignTexts.append(assignText)
          assignmentObj.append( (atomSets, ppms, peaks) )

        datum = [c, score, ' '.join(assignTexts)]
        textMatrix.append(datum)
        objectList.append(assignmentObj)

      self.matchMatrix.update(textMatrix=textMatrix, objectList=objectList)
    
  def displayMatchLocations(self):

    return 
    
    if self.assignment and self.spectrum and self.wPane1:
    
      optCache = self.analysisProject.doDetailAnnotations
      self.analysisProject.doDetailAnnotations = True
      
      N = range(len(self.rootDims))
      locations = []
      residue = self.spinSystem.residue
      
      for atomSets, ppms, peaks in self.assignment:
        dimMapping  = getDataDimAxisMapping(self.spectrum, self.wPane1)
        axisMapping = {}
        location    = {} 
        
        for axis in dimMapping:
          axisMapping[dimMapping[axis].dim-1] = axis
                  
        j = 0
        for dim in self.rootDims:
          location[axisMapping[dim]] = ppms[j]
          j  += 1
        
        locations.append(location)
      
        atomLabel = '%d:%s?' % (residue.seqCode,','.join([a.name for a in atomSets]))
 
        for peak in peaks:
          peak.details = atomLabel
 
      # Make strips, location['x'] = ppm
      displayStrips(self.parent, locations,
                    orthoPositions=None,
                    windowPane=self.wPane1)
      self.update_idletasks()
      displayStrips(self.parent, locations,
                    orthoPositions=None,
                    windowPane=self.wPane1)

      self.updateResonances()
 
      self.analysisProject.doDetailAnnotations = optCache
 
  def searchCarbonPlanes(self, indirectPositions):
    
    M = len(self.rootDims)
    peaksLocations = []

    tolerances = self.tolerances

    for atomSite in indirectPositions:
      location = indirectPositions[atomSite]
      region1  = list(self.fullRegion)
      
      j = 0
      for dim1, dim2 in self.diagDims: # root, nonRoot
        ppm = location[j]
        tol = tolerances[dim1]
        region1[dim1] = (ppm-tol,ppm+tol)
        j += 1
      
      for atomSite2 in indirectPositions:
        if atomSite == atomSite2:
          continue
          
        location2 =  indirectPositions[atomSite2]
        region2   = list(region1)
        
        i = 0
        for dataDim in self.nonRootDims:
          ppm = location2[i]
          tol = tolerances[j]/2.0
          region2[j] = (ppm-tol,ppm+tol)
          i += 1
          
        peaks2 = set()
        peaks2.update( searchPeaks([self.tPeakList,], region2) )
        peaks2.update( findPeaks(self.tPeakList, region2) )
        
        peaksLocations.append( (location2, peaks2) )

    return peaksLocations  
    
  def clusterRootLocations(self, newRootPeakLocations):
    
    rootDims    = self.rootDims
    nonRootDims = self.nonRootDims
    newRoots    = []
    tolerances  = self.tolerances
    
    done = {}
    for location2, peaks2 in newRootPeakLocations:
      for peak in peaks2:
        peakDims = peak.sortedPeakDims()
        #print "PEAK", ['%.3f' % x.value for x in peakDims], peak.serial
        if done.get(peak):
          continue
        else:
          done[peak] = True  
        
        position = []
        delta = 0.0
 
        i = 0
        for j in nonRootDims:
          d       = (peakDims[j].value - location2[i])/tolerances[j]
          delta  += d*d
          i += 1
 
        for i in rootDims:
          position.append(peakDims[i].value)
 
        newRoots.append((1.0, position, [peak,], delta))
        #print 'PRE-MERGED', ['%.3f' % x for x in position]
    
    # # # # # MERGE REPEAT FOUND H-C LOCATIONS # # # # # 
    
    N = len(newRoots)
    M = len(rootDims)
    if N > 1:
    
      i = 0
      while i < N-1:
        p1, position, peaks, delta = newRoots[i]
      
        j = i+1
        while j < N:
          p2, position2, peaks2, delta2 = newRoots[j]
          
          l = 0
          for k in range(M):
            tol  = tolerances[rootDims[k]]
            delta = abs(position[k]-position2[k])
                            
            if delta > tol:
              break
              
          else: # Merge clusters
            pT = p1+p2
            f1 = p1/pT
            f2 = p2/pT
            for l in range(M):
              position[l] = (f1*position[l]) + (f2*position2[l]) 

            N -= 1
            p1 = pT
            delta += delta2
            peaks.extend(peaks2)
            newRoots[i] = (p1, position, peaks, delta)
            newRoots.pop(j)
            continue
          
          j += 1
        i += 1
        
    testRoots = []
    for pT, position, peaks, delta in newRoots:
      #dev = delta/(pT) # Can use to select according to closeness to original indirect dim search ppm
      testRoots.append((position, peaks))
 
      #print 'ROOTMERGE', ['%.3f' % x for x in position], ','.join(['%d' % p.serial for p in peaks])
        
    return testRoots
        
  def unaliasRootLocations(self, rawRoots, indirectPositions, sideChainAtomSets):
    
    nonRootDims = self.nonRootDims
    rootDims    = self.rootDims
    diagDims    = self.diagDims
    tols        = self.tolerances
    fullRegion  = self.fullRegion
    specWidths  = self.specWidths
    project     = self.project
    ccpCode     = self.spinSystem.ccpCode
    
    R = range(len(rootDims))
    S = range(len(nonRootDims))
    matchThreshold = 0.00001
    doneRoots     = {}
    rootPositions = []
    
    bestMatchListScore = 0.0
    for position, peaks in rawRoots:   
      
      position2 = []
      for dim1, dim2 in diagDims: # root, nonRoot
        index = rootDims.index(dim1)
        position2.append(position[index])
         
      nonRootPositions = [position2, ]
      for peak in peaks:
        position2 = []
        peakDims = peak.sortedPeakDims() 

        for dim in nonRootDims:
          ppm1 = round(peakDims[dim].value,4)
          position2.append(ppm1)
 
        nonRootPositions.append(position2)
 
      matchListScore = 0.0
      for position2 in nonRootPositions:
        closest = None
        i = 0
        
        for dim in nonRootDims:
          ppm1 = position2[i]
          
          # Score by closeness to master list of indirect side chain locations
          for atomSite in indirectPositions:
            location = indirectPositions[atomSite]
            ppm2 = location[i]
            diff = (ppm1-ppm2)/tols[dim]
            delta2 = diff*diff
 
            if (closest is None) or (diff<closest[0]):
              closest = (delta2, location)
 
          i += 1
 
        e = closest[0]/twoGaus
        matchListScore += exp(-e)/gausRoot2pi
      
      if matchListScore > bestMatchListScore:
        bestMatchListScore = matchListScore
      
      #rootPositions.append( (position, nonRootPositions) )
      
      # # # # # # Deal with aliased positions
      altPositions = []
      i = 0
      for dim in rootDims:
        minVal, maxVal = fullRegion[dim]
        ppm  = position[i]
        sw   = specWidths[dim]
        pos2 = []
        
        ppm0 = ppm # include orig location
        while ppm0>minVal:
          pos2.append(round(ppm0,4))
          ppm0 -= sw  

        ppm0 = ppm+sw
        while ppm0<maxVal:
          pos2.append(round(ppm0,4))
          ppm0 += sw  
    
        altPositions.append(pos2)
    
        i += 1
    
      positions = [[],]
      for alternatives in altPositions:
        positions2 = []
      
        for position0 in positions:
          for alternative in alternatives:
            position2 = list(position0)
            position2.append(alternative)
            positions2.append(tuple(position2))
    
        positions = positions2
    
      # # # # # Check possible given atom types
      for position0 in positions:
        if doneRoots.get(position0):
          continue
    
        use = False
        for atomSets in sideChainAtomSets:
        
          matchScore = 0.0
          for k in R:
            ppm         = position0[k]
            atomSet     = atomSets[k]
            atomSetName = atomSet.name
            element     = atomSetName[0]
            score       = getAtomProbability(ccpCode, atomSetName, ppm)
 
            if score < matchThreshold: # e.g. if one of an H-C pair is impossible, both no good
              matchScore = 0.0
              break
            else:
              matchScore += score
    
          if matchScore:
            use = True
            break
            
        if use:
          doneRoots[position0] = True
          rootPositions.append( (position0, nonRootPositions, peaks, matchListScore) )
          #print 'ALIAS', ['%f' % x for x in position0], '%.3f' % matchListScore
  
    return rootPositions
          
  def clearRulers(self):
  
    for ruler in self.rulers:
      if not ruler.isDeleted:
        ruler.delete()
    
    self.rulers = []
      
  def commitAssignments(self):
  
    self.getTolerances()
  
    if self.assignment:
      for atomAssign in self.resonanceMatrix.currentObjects:
        self.commitSingleAssignment(atomAssign)
    
    self.pickKnownLocations()
         
  def commitAllAssignments(self):
  
    self.getTolerances()

    if self.assignment:
      for atomAssign in self.resonanceMatrix.objectList:
        self.commitSingleAssignment(atomAssign)
    
    self.pickKnownLocations()
    
  def commitSingleAssignment(self, atomAssign):
  
  
    atomH, atomC, resonanceH, resonanceC, shiftH, shiftC = atomAssign
  
    resonanceH = self.commitResonanceAssignment(atomH.name, resonanceH, '1H', shiftH)
    resonanceC = self.commitResonanceAssignment(atomC.name, resonanceC, '13C', shiftC)
    
    #resonanceH.setCovalentlyBound((resonanceC,))
    #if resonanceH not in resonanceC.covalentlyBound:
    #  resonanceC.addCovalentlyBound(resonanceH)
    
    #print "Commit", atomH.name, atomC.name, shiftH, shiftC, resonanceH.serial, resonanceC.serial
    
    self.updateAfter() # Update spin system table and selected atom assignments

  def commitResonanceAssignment(self, atomSetName, resonance, isotope, ppm):
  
    if resonance is None:
      for resonance0 in self.spinSystem.resonances:
        resonanceSet = resonance0.resonanceSet
        
        if resonanceSet:
          atomSet = resonanceSet.findFirstAtomSet()
           
          if atomSet.name == atomSetName:
            resonance = resonance0
            break
          
    if resonance is None:
      for resonance0 in self.spinSystem.resonances:
        if not resonance0.resonanceSet:
          if resonance0.isotopeCode == isotope:
            shift = resonance0.findFirstShift(parentList=self.shiftList)
          
            if shift and (abs(ppm-shift.value) < self.tolDict.get(isotope, 0.01)):
              resonance = resonance0
              break
            
      if not resonance:
        resonance = self.nmrProject.newResonance(isotopeCode=isotope)
    
    if not resonance.resonanceSet:
      atomSets = set()

      for atom in self.spinSystem.residue.atoms:
        atomSet = atom.atomSet
        
        if atomSet and (atomSet.name == atomSetName):
          atomSets.add(atomSet)
  
      if atomSets:
        assignAtomsToRes(list(atomSets), resonance)    
    
    shift = resonance.findFirstShift(parentList=self.shiftList)
    if not shift:
      shift = self.shiftList.newShift(resonance=resonance, value=ppm)
    
    return resonance
  
  def clearTempPeaks(self):
  
    for peak in self.tempPeaks:
      if not peak.isDeleted:
        if peak.peakList is self.tPeakList:
          peak.delete()
    
    self.tempPeaks = []    

  
def dynamicBranchPrune(I, J, linkScoreMatrix, linkThreshold, matchScoreMatrix, matchThreshold,
                       voidContextMatrix=None, voidLocations=None, ensembleSize=100, magic=10,
                       verbose=True, isLinearSystem=True):


  if not voidContextMatrix:
    voidContextMatrix = [[[],] * range(J) for x in range(I)]
                        
  baseScore  = linkThreshold + matchThreshold
  worstScore = baseScore * J * 2
  bestRoute  = []
  
  iNodes     = []
  neighbours = []
    
  i2 = I
  k = 0
  for j in range(J):
    iList = []

    for i in range(I):
      matchScore = matchScoreMatrix[i][j]
      
      if matchScore > matchThreshold:
        iNode = (k,j,i)
        k += 1
        iList.append(iNode)

    if not iList:
      iNode = (k,j,i2)
      k  += 1
      i2 += 1
      iList.append(iNode)      

    iNodes.append(iList)


  for j in range(J):
    for i in range(i2,J+1):
      iNode = (k,j,i)
      iNodes[j].append(iNode)      
      k  += 1

  I2 = max(i2,J+1)
   
  neighbours = [[] for x in range(k)]
  for j in range(J-1):
    j2 = j+1
    
    for ka, ja, ia in iNodes[j]: # All objects at this pos
      if ia >= I: # Null node
        for kb, jb, ib in iNodes[j2]: # All objs at next pos
          if ia == ib:
            continue
            
          neighbours[ka].append((kb, jb, ib))          
          
      else:
        for kb, jb, ib in iNodes[j2]: # All objs at next pos
          if ia == ib:
            continue

          if ib >= I:
            neighbours[ka].append((kb, jb, ib))

          else:
            linkScore = linkScoreMatrix[ia][ib]
  
            if linkScore > linkThreshold:
              neighbours[ka].append((kb, jb, ib))

  matchScoreMatrix2  = [[matchThreshold] * J for x in range(I2)]
  linkScoreMatrix2   = [[linkThreshold] * I2 for x in range(I2)]
  voidContextMatrix2 = [[[],] * J for x in range(I2)]

  for i in range(I):
    for i2 in range(I):
      linkScoreMatrix2[i][i2] = linkScoreMatrix[i][i2]
    for j in range(J):
      matchScoreMatrix2[i][j]  = matchScoreMatrix[i][j]
      voidContextMatrix2[i][j] = voidContextMatrix[i][j]

  matchScoreMatrix  = matchScoreMatrix2 
  linkScoreMatrix   = linkScoreMatrix2  
  voidContextMatrix = voidContextMatrix2

  route     = []
  score     = []
  routes    = []
  openList  = iNodes[0]
  scores    = {}
  bestScore = worstScore

  #c = 0
  while openList:
    
    #c += 1
    # 
    # if c > 1000000:
    #   break
    
    node = openList.pop(0)
    k, j, i = node

    route = route[:j]
    score = score[:j]
    
    matchScore  = matchScoreMatrix[i][j]
    
    if route:
      contextVoid = voidContextMatrix[i][j]
      i2 = route[-1]
      if contextVoid and (i2 in contextVoid):
        continue
        
      if isLinearSystem:
        linkScore = linkScoreMatrix[i2][i]
      else:
        linkScore = sum([linkScoreMatrix[i2][i] for i2 in route])
        # Adding up all link permutations - fair in this context
            
    else:
      linkScore = linkThreshold
    
    route.append(i)
    score.append(matchScore+linkScore)
    
    key = frozenset(route[-magic:]+[(j,i),])

    #totalScore = (J-i) * baseScore
    totalScore = sum(score) + (J-j) * baseScore
    
    prevScore = scores.get(key, worstScore)
      
    if totalScore > prevScore:
      scores[key] = totalScore
      
      if len(scores) > 20000:
        scores.popitem()
      
      if totalScore > bestScore:
        bestScore = totalScore
        routes.append((totalScore,list(route)))
        #print route, neighbours[k]
      
    else:
      continue
    
    nextNodes = neighbours[k]

    if j < J-1:
      nextNodes = [n for n in nextNodes if n[2] not in route]
      openList  = nextNodes + openList
      
    else:
      routes.append((totalScore,list(route)))
    
  routes.sort()
  routes.reverse()
      
  return routes
 
 

def iteration(C,A,R):

  m, n = A.shape
  # shape = ss = rows, rr = cols
  #R2 = R*R
  
  Ri = linalg.inv(R)
  
  CA = C*A
  AR = A*R
  
  B1 = multiply(CA,AR)
  F = B1*Ri
  
  # Remove negative
  F = F + abs(F)
  F = multiply(A,F)
  F /= F.sum()
  
  """
  l = []
  for i in xrange(m):
    for j in xrange(n):
      l.append((F[i,j],(i,j)))
  
  l.sort()
  #l.reverse()
  
  while l:
    v, (i,j) = l.pop()
    
    if v:
 
      s = F[:,j].sum()
      if s:
        F[:,j]/s

      s = F[i].sum()
      if s:
        F[i]/s

  # On spin system
  #n, m = F.shape
 
  for s in xrange(m):
    s3 = F[s].sum()
    if s3:
      F[s] /= s3
  
  # On residue
  F = F.T
  #n, m = F.shape
  for s in xrange(n):
    s3 = F[s].sum()
    if s3:
      F[s] /= s3
  F = F.T
  
      
  """
  
  
  return F



"""
======================COPYRIGHT/LICENSE START==========================

CalcCouplings.py: Part of the CcpNmr Analysis program

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
from ccpnmr.analysis.core.CouplingBasic import assignPrimaryClusterCoupling, getClusterCoupling, deleteCluster, \
                          MULTIPLET_PEAK_DICT, getCoupledExpDims, EXPT_MULTIPLETS,\
                          makeMultipletPeakCluster, getCouplingAnnotations, getPeakListsPeakClusters, \
                          COUPLING_TYPES, getSpectrumMultipletPattern, setSpectrumMultipletPattern

from ccpnmr.analysis.core.WindowBasic import getDataDimAxisMapping, getWindowPaneName

from ccpnmr.analysis.popups.BasePopup import BasePopup

from memops.gui.ButtonList          import ButtonList, UtilityButtonList
from memops.gui.IntEntry            import IntEntry
from memops.gui.Label               import Label
from memops.gui.LabelFrame          import LabelFrame
from memops.gui.MessageReporter     import showWarning, showOkCancel
from memops.gui.PartitionedSelector import PartitionedSelector
from memops.gui.PulldownList        import PulldownList
from memops.gui.ScrolledMatrix      import ScrolledMatrix


def testCalcCouplingPopupMacro(argServer):

  popup = CalcCouplingPopup(argServer.parent)
  popup.open()

# Pattern prediction - needed?
#     previous pattern can be readily determined for the spectra of existing peak clusters
#     check for uniformity of measurement, significant absolute size, max num of peaks
#     per spectrum and sign. Can score a pattern by laying out a grid based on coupling value
#     peaks should be close to intersections and have the right sign
#     can check several patterns and several clusters
#     limit to patterns comes from num peaks per spec
#     Warn if type is changed but clusters exist

# Auto function to predict clusters from picked peaks - COSY picker will be an alternative
#   - needs expected coupling ranges in addition to multiplet pattern
#     Min and max Hz for each coupled dimension
#     Search goes through all peaks checking for the Hz distance and peak sign
#     Builds on initial matches to complete the pattern (just direction and sign)
#     Can use same grid method as pattern predictor, use average grid with tolerances initially
#     Then try grid sizes within bounds.
#     Check equivalency of absolute intensity
#   - will inform of ambiguity of cluster peaks - to be sorted manually

# Essential
# 1, Tabbed popup
#   Spectra Tab:
#     List of coupled spectra
#     - expected multiplet pattern for each
#     - uses defaults for each spectrum type
#     - coupling measurement list for each c.f. shift list
#
#   Couplings Tab:
#     Select peakList
#     Select windowPane for finding peaks
#     List of clusters and associated couplings
#     Make J coupling list func
#
#   RDC Alignment Tab:
#     Set isotropic/unordered sample
#     Add/remove alignment media and associated peakList(s) - Two for IPAP

#   Calc RDCs Tab:
#     Set isotropic sample
#     Set aligned sample
#     - Effectively gives peak-cluster list for each
#     - Check cluster integrity
#     - Connect across the different lists with centre resonance assignments,
#       or failing that use centre location, within some tolerance
#     Table of subtracted couplings per resonance/centre
#     - list residual Hz
#     - resonance assignment for the coupling,
#     - main peak assignment
#     Create RDC constraint list from the residual data
# 2, Macro and window menu functions
# 3, On screen annotation
# 4, Add more refExperiment data

# Maybe
# 5, Auto clustering
#   - specify min/max Hz ranges per dimension
#   - already have multiplet types
#   - only do the easy ones with unique matches
# 6, Pattern prediction




class CalcCouplingPopup(BasePopup):

  def __init__(self, parent, *args, **kw):
  
    self.waiting   = False
    self.peakList  = None
    self.peakLists = []
    self.windowPane = None
    self.cluster   = None
    self.clusters  = {}
    self.multiplet = None
    self.guiParent = parent
    
    BasePopup.__init__(self, parent=parent, title='Measure Couplings', **kw)
  
    self.updateWindows()
    self.updateAfter()

  def body(self, guiFrame):

    guiFrame.grid_columnconfigure(0, weight=1)
    
    row = 0
    frame = LabelFrame(guiFrame, text='Options')
    frame.grid(row=row, column=0, sticky='ew')
    frame.grid_columnconfigure(1, weight=1)
    frame.grid_rowconfigure(1, weight=1)
    
    label = Label(frame, text='Window:')
    label.grid(row=0, column=0, sticky='nw')
    self.windowPulldown = PulldownList(frame, callback=self.changeWindow)
    self.windowPulldown.grid(row=0, column=1, sticky='nw')

    label = Label(frame, text='Multiplet Pattern:')
    label.grid(row=1, column=0, columnspan=2, sticky='nw')
    self.multipletButtons = PartitionedSelector(frame, callback=self.setMultiplet,
                                                radio=True)
    self.multipletButtons.grid(row=2, column=0, columnspan=2, sticky='ew')
  
    row += 1
    frame = LabelFrame(guiFrame, text='Active Peak Lists')
    frame.grid(row=row, column=0, sticky='ew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)
    
    headingList = ['Experiment','Spectrum','List','Coupled Dims','Experiment Type']
    self.peakListMatrix = ScrolledMatrix(frame,
                                         headingList=headingList,
                                         callback=None,
                                         multiSelect=False)
    self.peakListMatrix.grid(row=0, column=0, sticky = 'nsew')

    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    frame = LabelFrame(guiFrame, text= 'Multiplet Peak Clusters')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)
    
    headingList = ['#','Main\nAssignment','Num\nPeaks',
                   'Coupling\nAssignment','Value','Value']
    self.clusterMatrix = ScrolledMatrix(frame,
                                        headingList=headingList,
                                        callback=self.selectCluster,
                                        multiSelect=True)
    self.clusterMatrix.grid(row=0, column=0, sticky = 'nsew')

    row += 1
    texts    = ['Cluster Selected\nPeaks','Assign\nCouplings',
                'List\nPeaks','Find\nPeaks','Delete\nClusters']
    commands = [self.clusterSelectedPeaks, self.assignCouplings,
                self.showPeaks, self.findPeaks,self.deleteClusters]
    self.bottomButtons = UtilityButtonList(guiFrame, texts=texts, expands=True,
                                           commands=commands, helpUrl=self.help_url)
    self.bottomButtons.grid(row=row, column=0, sticky = 'ew')

    self.administerNotifiers(self.registerNotify)

  def administerNotifiers(self, notifyFunc):

    for clazz in ('ccp.nmr.Nmr.DataSource', 'ccp.nmr.Nmr.Experiment'):
      notifyFunc(self.updatePeakListsAfter, clazz, 'setName')
   
    for func in ('__init__', 'delete','setName'):
      notifyFunc(self.updateWindows, 'ccpnmr.Analysis.SpectrumWindow', func)

    for func in ('__init__', 'delete'):
      notifyFunc(self.updatePeakListsAfter, 'ccp.nmr.Nmr.PeakList', func)
 
    for func in ('__init__', 'delete','addPeak','setPeaks','removePeak','setAnnotation'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.PeakCluster', func)
 
  def deleteClusters(self):
  
    clusters = self.clusterMatrix.currentObjects
    for cluster in clusters:
      deleteCluster(cluster)
 
 
  def findPeaks(self):
  
    if self.windowPane and self.cluster:
      peaks = list(self.cluster.peaks)
      if not peaks:
        return
      
      spectrum = peaks[0].peakList.dataSource
      analysisSpectrum = spectrum.analysisSpectrum
      view = self.windowPane.findFirstSpectrumWindowView(analysisSpectrum=analysisSpectrum)
      windowFrame = self.windowPane.getWindowFrame()
     
      position = {}
      n = float(len(peaks))
     
      for axisMapping in view.axisMappings:
        dim = axisMapping.analysisDataDim.dataDim.dim
        xyz = axisMapping.label
      
        mean = 0.0
        for peak in peaks:
          peakDim = peak.findFirstPeakDim(dim=dim)
          mean += peakDim.realValue
      
        mean /= n
        position[xyz] = mean
        
      windowFrame.gotoPosition(position)  
 
  def updatePeakListsAfter(self, object):
  
    if object.className == 'Experiment':
      for spectrum in object.dataSources:
        for peakList in spectrum.peakLists:
          if peakList in self.peakLists:
            self.updatePeakLists()
    
    elif object.className == 'DataSource':
      for peakList in object.peakLists:
        if peakList in self.peakLists:
          self.updatePeakLists()
  
    else:
      if object in self.peakLists:
        self.updatePeakLists(object)

  def open(self):
  
    BasePopup.open(self)

    self.updateWindows()
    self.updateAfter()
  
  def showPeaks(self):
  
    peaks = {}  
    for peakCluster in self.clusterMatrix.currentObjects:
      for peak in peakCluster.peaks:
        peaks[peak] = True
    
    if peaks:
      self.parent.viewPeaks(peaks.keys())
  
  def getWindows(self):
    
    windowPanes = []
    
    for window in self.analysisProject.sortedSpectrumWindows():
      for windowPane in window.sortedSpectrumWindowPanes():
        for view in windowPane.spectrumWindowViews:
          if view.isPosVisible or view.isNegVisible:
            spectrum = view.analysisSpectrum.dataSource

            if self.getCoupledExpDims(spectrum.experiment):
              windowPanes.append(windowPane)
              break
      
    return windowPanes
    
  def getCoupledExpDims(self, experiment):
    """ Descrn: List the dimensions (ExpDims) of an experiment which carry couplings.
        Inputs: Nmr.Experiment
        Output: List of Nmr.ExpDims
    """
    COUPLING_TYPES = ('JCoupling', 'Rdc', 'DipolarCoupling') 

    expDims = []
    refExperiment = experiment.refExperiment
 
    if refExperiment:
      for expDim in experiment.expDims:
        refExpDim = expDim.refExpDim
 
        if not refExpDim: # Could be sampled dim
          continue
 
        for refExpDimRef in refExpDim.refExpDimRefs:
          if refExpDimRef.coupledIsotopeCodes:
            expDims.append(expDim)
            break
    
    if not expDims:
      for expDim in experiment.expDims:
        for expDimRef in expDim.expDimRefs:
          if expDimRef.measurementType in COUPLING_TYPES:
            expDims.append(expDim)
            break
 
    return expDims
            
  def changeWindow(self, windowPane):
  
    if windowPane is not self.windowPane:
      self.windowPane = windowPane
      self.updatePeakLists()
  
  def updateWindows(self, *object):
  
    panes = self.getWindows()
    index   = 0
    names   = []
    pane  = self.windowPane
 
    if panes:
      names  = [getWindowPaneName(wp) for wp in panes]
     
      if pane not in panes:
        pane = panes[0]
        
      index = panes.index(pane) 
 
    if pane is not self.windowPane:
      self.windowPane = pane
      self.updatePeakLists()
 
    self.windowPulldown.setup(names, panes, index)
  
  def getPeakLists(self):
  
    peakLists = []
    
    if self.windowPane:
      for view in self.windowPane.spectrumWindowViews:
        try:
          spectrum = view.analysisSpectrum.dataSource
        except:
          continue
        
        if spectrum.peakLists:
          peakList = spectrum.activePeakList  
        
          if not peakList:
            peakList = spectrum.findFirstPeakList()
          
          peakLists.append(peakList)
    
    return peakLists
  
  def updatePeakLists(self, peakList=None):
  
    objectList = self.getPeakLists()
    textMatrix = []
    spectra = {}
    
    for peakList0 in objectList:
      spectrum = peakList0.dataSource
      experiment = spectrum.experiment
      spectra[spectrum] = True
      refExperiment = experiment.refExperiment
      coupledDims   = []
      if refExperiment:
        experimentType = refExperiment.name
        
        for expDim in experiment.expDims:
          refExpDim = expDim.refExpDim
          isotopes  = {}
          
          for refExpDimRef in refExpDim.refExpDimRefs:
            for isotope in refExpDimRef.coupledIsotopeCodes:
              isotopes[isotope] = True
          
          if isotopes:
            isoString = ','.join(isotopes)
            coupledDims.append('%d (%s)' % (expDim.dim, isoString))
          
      else:
        experimentType = None
        
      if not coupledDims:
        coupledDims = None
      else:
        coupledDims = ' '.join(coupledDims)
    
      datum = [experiment.name,
               spectrum.name,
               peakList0.serial,
               coupledDims,
               experimentType
              ]
      textMatrix.append(datum)
    
    self.peakLists = objectList
    self.peakListMatrix.update(objectList=objectList, textMatrix=textMatrix)
    
    self.multiplet = None
    for spectrum in spectra.keys():
      multiplet = self.getSpectrumMultiplet(spectrum)
      
      if multiplet:
        self.multiplet = MULTIPLET_PEAK_DICT.get(multiplet)
        break
    
    self.updateMultipletButtons()
    self.updateAfter()
  
  def setMultiplet(self, pattern):
  
    for peakList in self.peakLists:
      spectrum = peakList.dataSource
      
      for name in MULTIPLET_PEAK_DICT.keys():
        if MULTIPLET_PEAK_DICT[name] == pattern:
          setSpectrumMultipletPattern(spectrum, name)
          break
  
    self.multiplet = pattern
  
  def setSpectrumMultiplet(self, spectrum, pattern):
  
    prevPattern = getSpectrumMultipletPattern(spectrum)
    
    if prevPattern and (prevPattern != pattern):
      if showOkCancel('Query',
                      'Really change multiplet pattern?', parent=self):
  
        setSpectrumMultipletPattern(spectrum, pattern)
  
  def getSpectrumMultiplet(self, spectrum):
  
    name = getSpectrumMultipletPattern(spectrum)
    
    if not name:
      name = self.predictSpectrumMultiplet(spectrum)
      setSpectrumMultipletPattern(spectrum, name)
      
    return name
    
  def predictSpectrumMultiplet(self, spectrum):
  
    name = None
    
    refExperiment = spectrum.experiment.refExperiment
    if refExperiment:
      name = EXPT_MULTIPLETS.get(refExperiment.name)
      
    return name
  
  def updateMultipletButtons(self):
    
    import re
    
    labels   = []
    objects  = []
    colors   = []
    selected = None
    
    if self.windowPane:
      coupledAxes = []
      spectra = {}
      for peakList in self.peakLists:
        spectra[peakList.dataSource] = True
        multiplet = self.getSpectrumMultiplet(peakList.dataSource)
        selected = MULTIPLET_PEAK_DICT[multiplet]
        
      for spectrum in spectra.keys():
        mapping = getDataDimAxisMapping(spectrum, self.windowPane)
 
        for axisLabel in ('x','y'):
          dataDim = mapping[axisLabel]
 
          for expDimRef in dataDim.expDim.expDimRefs:
            if expDimRef.refExpDimRef and expDimRef.refExpDimRef.coupledIsotopes:
              if axisLabel not in coupledAxes:
                coupledAxes.append(axisLabel)
                break
                            
    else:
      coupledAxes = ['x','y']
        
    for name in MULTIPLET_PEAK_DICT.keys():
    
      pattern = MULTIPLET_PEAK_DICT[name]

      if type(pattern[0]) in (type(()), type([])):
        if 'y' not in coupledAxes:
          continue
        
        if len(pattern[0]) > 1: # x & y
          if 'x' not in coupledAxes:
            continue
            
        else: # y only
          if 'x' in coupledAxes:
            continue
          
          
      else: # x only
        if 'x' not in coupledAxes:
          continue

        if 'y' in coupledAxes:
          continue
    
      label = re.sub('\|','\n',name)
    
      objects.append(pattern)
      labels.append(label)
      colors.append('#8080FF')
            
    self.multipletButtons.update(objects=objects, colors=colors,
                                 labels=labels)
   
    self.multipletButtons.setSelected([selected,])
   
   
  def selectPeakList(self, object, row, col):
  
    self.peakList = object
  
    self.updateButtons()

  def selectCluster(self, object, row, col):
  
    self.cluster = object
  
    self.updateButtons()
  
  def assignCouplings(self):
 
    for cluster in self.clusterMatrix.currentObjects:
      assignPrimaryClusterCoupling(cluster)
  
  def assignAllCouplings(self):
  
    for cluster in self.clusters.keys():
      assignPrimaryClusterCoupling(cluster)
  
  def clusterSelectedPeaks(self):
  
    peaks0 = self.parent.currentPeaks
    peaks = []
    for peak in peaks0:
      if peak.peakList in self.peakLists:
        peaks.append(peak)
    
    if not peaks:
      showWarning('Cluster Failure',
                  'No peaks selected from active peak lists', parent=self)
      return
    
    if not self.multiplet:
      showWarning('Cluster Failure',
                  'No multiplet pattern selected', parent=self)
      return
    
    cluster = makeMultipletPeakCluster(peaks, self.multiplet, self.windowPane)
      
    if cluster:
      self.clusterMatrix.selectObject(cluster)
 
  def getActiveClusterCouplings(self):
  
  
    clusters = {}
    if self.peakLists:
      dims = {}
      
      for peakList in self.peakLists:
        experiment = peakList.dataSource.experiment
        
        for expDim in getCoupledExpDims(experiment):
          dims[expDim.dim] = True
      
      dims = dims.keys()
        
      for cluster in getPeakListsPeakClusters(self.peakLists):
        values = []
        
        for dim in dims: 
          values.append( getClusterCoupling(cluster, dim) )
          # Allow Nones? 
          
        clusters[cluster] = values
  
    self.clusters = clusters
    
    return clusters
  
  def updateButtons(self):
  
    pass
  
  def updateAfter(self, cluster=None):
  
    if self.waiting:
      return
    
    if cluster:
      for peak in cluster.peaks:
        if peak.peakList in self.peakLists:
          self.waiting = True
          self.after_idle(self.update)
          break
    
    else:
      self.waiting = True
      self.after_idle(self.update)
 
  
  def update(self):
  
    clusters = self.getActiveClusterCouplings()
  
    textMatrix = []
    objectList = []
    
    headingList = ['#','Main\nAssignment','Num\nPeaks','Coupling\nAssignment']
     
    if clusters:
      cluster0 = clusters.keys()[0]
      i = 1
      
      for value in clusters[cluster0]:
        headingList.append('F%d Value\n(Hz)' % i)
        i += 1
      
    else:
      headingList.append('Value')
      
    for cluster in clusters.keys():
 
      couplingAssign = ','.join(getCouplingAnnotations(cluster))
    
      datum = [cluster.serial,
               cluster.annotation,
               len(cluster.peaks),
               couplingAssign]
               
      values = clusters[cluster]
      for value in values:
        datum.append(value)
        
      textMatrix.append(datum)
      objectList.append(cluster)

    self.clusterMatrix.update(headingList=headingList,
                              textMatrix=textMatrix,
                              objectList=objectList)
  
    self.updateButtons()
    self.waiting  = False
  
  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
  
    BasePopup.destroy(self)
    
  
  
  

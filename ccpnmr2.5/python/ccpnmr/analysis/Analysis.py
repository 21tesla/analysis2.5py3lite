
"""
======================COPYRIGHT/LICENSE START==========================

Analysis.py: Part of the CcpNmr Analysis program

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
import os, time, re, sys, traceback, cPickle, shutil

from memops.general import Implementation

from memops.general.Application import Application
from memops.general.Io import backupProject
# from memops.general import Version

from memops.universal.Io import joinPath
from memops.universal.Util import isWindowsOS
from memops.universal.Region1D import Region1D

from memops.gui.MessageReporter import showYesNo
from memops.gui.Color import hexToHsb, hexToRgb

from memops.api.Implementation import Url

from ccp.api.nmr import Nmr

from ccp.general.Io import getDataSourceFileName
from ccp.util.ShapeUtil import create1dShape

# from ccpnmr.api.AnalysisProfile import Macro

from ccpnmr.analysis.Citation import analysisReference
from ccpnmr.analysis.macros.ArgumentServer import ArgumentServer
from ccpnmr.analysis.core.BlockUtil import getBlockFile, getShapeBlockFile
from ccpnmr.analysis.core import AssignmentBasic
from ccpnmr.analysis.core import CouplingBasic
from ccpnmr.analysis import Copyright
from ccpnmr.analysis.core import ExperimentBasic
from ccpnmr.analysis.core import MarkBasic
from ccpnmr.analysis.core import MoleculeBasic
from ccpnmr.analysis.core import PeakBasic
from ccpnmr.analysis.core import WindowBasic
from ccpnmr.analysis.core import UnitConverter
from ccpnmr.analysis.core import Util

try:
  # NOTE:ED - quick test to find which environment variables are missing
  # print
  # print 'ENVIRONMENT VARIABLES'
  # for k, v in sorted(os.environ.items()):
  #   if ';' in v:
  #     print '>>>', k
  #     vList = v.split(';')
  #     for vl in vList:
  #       print '      ', vl
  #   else:
  #     print '>>>', k, v
  # print
  from memops.c.MemCache import MemCache
  import memops.c.StoreFile as StoreFile
  from ccpnmr.c.ContourFile import ContourFile, StoredContourFile
  from ccpnmr.c.PeakList import PeakList as CPeakList
  from ccpnmr.c.PeakCluster import PeakCluster as CPeakCluster
  from ccpnmr.c.SliceFile import SliceFile
  from ccpnmr.c.WinPeakList import WinPeakList as CWinPeakList

except Exception, e:
  print 'Error, the Analysis module will not work, something is wrong with the C code.'
  print 'Exception:', e
  print 'Will continue without Analysis C functionality'
  MemCache = StoreFile = ContourFile = StoredContourFile = CPeakList = CPeakCluster = SliceFile = CWinPeakList = None

X_Y = set(['x', 'y'])
MULTIPLET = 'multiplet'
MULTIPLET_ENUM = 0
LOCAL_HELP_DOC_DIR =  '/python/ccpnmr/analysis/doc/build/html/'

class Analysis:

  def __init__(self, cache_size):

    self.cache_size = cache_size

    name = Copyright.suite + '_' + Copyright.program
    alternativeNames = ('CCPNMR_Analysis2',
                        'CCPNMR_Analysis') # previous names
    self.application = Application(name=name, alternativeNames=alternativeNames)
    self.project = None
    self.mem_cache = None
    self.shiftList = None
    self.currentPeak = None
    self.currentPeaks = []
    self.currentPeakLists = []
    self.currentSpectra = []
    self.currentSpectrum = None
    self.currentWindow = None
    self.currentWindowPane = None
    self.currentWindowPopup = None
    self.currentWindowFrame = None
    self.currentCanvas = None
    self.currentEvent = None
    self.currentPosition = None
    self.currentRegion = None
    self.changedAxisPanel = None

    # below are for the NewWindow popup so 
    # that know about visibility before windows created
    self.visibleSpectra = {}
    self.toolbarSpectra = {}

    self.showYesNo = showYesNo

    version = Copyright.version
    self.versionInfo = '%s %s Version %s.%s.%s Release %s (%s)' % (Copyright.suite, Copyright.program, version.major, version.minor, version.level, version.release, Copyright.copyright)
    if Copyright.build:
      self.versionInfo += ('\n\t\tWindows Build number %s' % (Copyright.build))

    self.printCommandLineInfo()

    # below are done in initProject
    ###AssignmentBasic.setupAssignmentNotifiers()
    ###CouplingBasic.setupCouplingNotifiers()

  def after(self, delay, func, *args):

    raise Exception('need to override in subclass')

  def after_idle(self, func, *args):

    raise Exception('need to override in subclass')

  def after_cancel(self, eventId):

    raise Exception('need to override in subclass')

  def printCommandLineInfo(self):
    
    version = Copyright.version
    print self.versionInfo
    if version.timestamp:
      print 'Distribution created %s' % version.timestamp
    
    print analysisReference

  def checkAxisPanels(self):

    # check that orthogonal axisPanels have
    # number of regions matching x axisPanel
    for window in self.analysisProject.spectrumWindows:
      for windowPane in window.spectrumWindowPanes:
        axisPanels = windowPane.sortedAxisPanels()
        axisPanel  = windowPane.findFirstAxisPanel(label=window.stripAxis)
        n = len(axisPanel.axisRegions)

        for axisPanel in axisPanels[2:]:
          if not axisPanel.axisType.isSampled:
            axisRegions = axisPanel.sortedAxisRegions()
            m = len(axisRegions)
            
            if m < n:
              self.addAxisRegions(axisPanel, n-m)
            
            elif m > n:
              for axisRegion in axisRegions[n:]:
                axisRegion.delete()


  def addAxisRegions(self, axisPanel, n):

    axisRegions = axisPanel.sortedAxisRegions()
    #region = axisRegions[-1].region
    #size = axisRegions[-1].size
    region = axisRegions[0].region
    size = axisRegions[0].size
    for i in range(n):
      Util.addAxisPanelRegion(axisPanel, region=region, size=size)

  def addedAxisRegion(self, axisRegion):

    # if have added extra stripAxis axisRegion
    # then add orthogonal ones
    axisPanel = axisRegion.axisPanel
    #print 'Analysis.addedAxisRegion', axisPanel.label, axisPanel.spectrumWindow.stripAxis
    
    if axisPanel.label == axisPanel.spectrumWindowPane.spectrumWindow.stripAxis:
      axisPanels = axisPanel.spectrumWindowPane.sortedAxisPanels()
      for axisPanel in axisPanels[2:]:
        ###if Util.isIsotopeAxisType(axisPanel.axisType):
        if not axisPanel.axisType.isSampled:
          self.addAxisRegions(axisPanel, 1)

  def deletedAxisRegion(self, axisRegion):

    pass
    ''' this now handled in WindowDraw
    # if have deleted extra x axisRegion then delete orthogonal ones
    axisPanel = axisRegion.axisPanel
    if (axisPanel.label == 'x'):
      axisPanels = axisPanel.spectrumWindowPane.sortedAxisPanels()[2:]
      for axisPanel in axisPanels:
        axisPanel.axisRegions[-1].delete()
    '''

  def initTopObjects(self, project):

    if not project.currentInstrumentStore:
      instrumentStores = project.sortedInstrumentStores()
      if instrumentStores:
        project.currentInstrumentStore = instrumentStores[0]
      else:
        # this automatically sets project.currentInstrumentStore
        project.newInstrumentStore(name=project.name)

    if not project.currentNmrProject:
      nmrProjects = project.sortedNmrProjects()
      if nmrProjects:
        project.currentNmrProject = nmrProjects[0]
      else:
        # this automatically sets project.currentNmrProject
        project.newNmrProject(name=project.name)
        
    if not project.currentAnalysisProject:
      analysisProjects = project.sortedAnalysisProjects()
      if analysisProjects:
        project.currentAnalysisProject = analysisProjects[0]
      else:
        # this automatically sets project.currentAnalysisProject
        analysisProject = project.newAnalysisProject(name=project.name, nmrProject=project.currentNmrProject)
        analysisProject.maxMarks = analysisProject.maxRulers = 10

    if not project.currentAnalysisProfile:
      analysisProfiles = project.sortedAnalysisProfiles()
      if analysisProfiles:
        project.currentAnalysisProfile = analysisProfiles[0]
      else:
        # this automatically sets project.currentAnalysisProfile
        project.newAnalysisProfile(name=project.name)
    
    self.nmrProject = project.currentNmrProject
    self.analysisProject = project.currentAnalysisProject
    self.analysisProfile = project.currentAnalysisProfile
    self.analysisProject.isThreadingAllowed = False
    #self.analysisProject.doOneLetterAnnotations = False
    #self.analysisProject.doCompressAnnotations = True
    #self.analysisProject.doLeadingNumberAnnotations = True


  def initProject(self, project):

    self.project = project
    self.mem_cache = MemCache and MemCache(self.cache_size)
    
    if not hasattr(self, 'argumentServer'):
      self.argumentServer = ArgumentServer(self)
      self.argServer = self.argumentServer

    project.application = self.application

    if os.environ.has_key('USER'):
      user = os.environ['USER']
      user = user.replace(' ', '_')
      try:
        project.currentUserId = user
      except:
        pass

    self.backup_alarm_id = None

    self.visibleSpectra = {}
    self.toolbarSpectra = {}

    self.initTopObjects(project)
    self.checkAxisPanels()
    
    Util.getSoftware(project) # just so that it is initialised always
    Util.defaultColors(project)
    Util.defaultAxisTypes(project)
    Util.defaultWindowGroup(project)

    notify = self.registerNotify
    notify(self.initPeakList, 'ccp.nmr.Nmr.PeakList', '__init__')

    self.initSpectra()

    # Re-average shifts - done centrally instead of doing it in initPeakLists
    for shiftList in set(x.shiftList for x in project.currentNmrProject.experiments):
      AssignmentBasic.updateAllShifts(shiftList)
    
    if not Util.isTopObjectAnalysisInitialised(project.currentNmrProject):
      self.initResonances()
      self.initSpinSystems()
      self.initShiftLists()

    self.initPeakClusters()
    self.initProfile()
    self.initSpectrumViews()

    if not Util.isTopObjectAnalysisInitialised(project.currentAnalysisProject):
      self.initAtomSetMappings()

    self.initMolSystems()
    self.currentPeaks = []
    self.currentPeak = None
    self.currentPeakDim = None
    self.currentSpectrum = None
    self.currentSpectra = []
    self.currentPeakLists = []

    self.setupBackup()

    self.axisRegionGroups = []

    self.curateNotifiers(notify)

    AssignmentBasic.setupAssignmentNotifiers()
    CouplingBasic.setupCouplingNotifiers()
    #self.project.checkAllValid()

  def registerNotify(self, notifyFunc, classname, funcname, keyword = None):

    self.application.registerNotify(notifyFunc, classname, funcname, keyword)

  def unregisterNotify(self, notifyFunc, classname, funcname, keyword = None):

    self.application.unregisterNotify(notifyFunc, classname, funcname, keyword)

  def curateNotifiers(self, notify):

    appDataNotify = notify

    ### TBD v2: probably need more
    
    notify(self.initProfile, 'memops.Implementation.MemopsRoot', 'setCurrentAnalysisProfile')
    for func in ('__init__','delete','setGuiName'):
      notify(self.cacheResidueCode,'ccpnmr.AnalysisProfile.ResidueProfile', func)
    
    for func in ('setNumPoints', 'setNumPointsOrig', 'setIsComplex', 'setPointOffset'):
      notify(self.updatedDataDimDetails, 'ccp.nmr.Nmr.FreqDataDim', func)
    for func in ('setNumPoints', 'setIsComplex'):
      notify(self.updatedDataDimDetails, 'ccp.nmr.Nmr.SampledDataDim', func)
    for func in ('setIsBigEndian', 'setNByte', 'setNumberType', 'setHeaderSize', 'setHasBlockPadding', 'setBlockSizes', 'setPath', 'setNumPoints', 'setDataUrl'):
      notify(self.updatedNumericMatrixDetails, 'ccp.general.DataLocation.BlockedBinaryMatrix', func)
    for func in ('setNumberType', 'setHeaderSize', 'setPath', 'setNumPoints', 'setDataUrl'):
      notify(self.updatedNumericMatrixDetails, 'ccp.general.DataLocation.ShapeMatrix', func)
    notify(self.updatedDataUrl, 'ccp.general.DataLocation.DataUrl', 'setUrl')

    notify(self.initMolSystemChain, 'ccp.molecule.MolSystem.Chain', '__init__')
    notify(self.deletedMolSystemChain, 'ccp.molecule.MolSystem.Chain', 'delete')
    notify(self.refreshPeakResNumberAfter, 'ccp.molecule.MolSystem.Residue', 'setSeqCode')
    for func in ('__init__','setResonances','addResonance','removeResonance'):
      notify(self.initResonanceSet, 'ccp.nmr.Nmr.ResonanceSet', func)
    notify(self.deletedResonanceSet, 'ccp.nmr.Nmr.ResonanceSet', 'delete')
    notify(self.initExperiment, 'ccp.nmr.Nmr.Experiment', '__init__')
    notify(self.initSpectrum, 'ccp.nmr.Nmr.DataSource', '__init__')
    notify(self.initAxisPanel, 'ccpnmr.Analysis.AxisPanel', '__init__')
    notify(self.deletedAxisPanel, 'ccpnmr.Analysis.AxisPanel', 'delete')
    notify(self.deletedPeakList, 'ccp.nmr.Nmr.PeakList', 'delete')
    notify(self.updateChemShifts, 'ccp.nmr.Nmr.DataDimRef', '')
    notify(self.initWinPeakList, 'ccpnmr.Analysis.WindowPeakList', '__init__')
    notify(self.deletedWinPeakList, 'ccpnmr.Analysis.WindowPeakList', 'delete')
    notify(self.initSpectrumWindow, 'ccpnmr.Analysis.SpectrumWindow', '__init__')
    notify(self.initSpectrumWindowView, 'ccpnmr.Analysis.SpectrumWindowView', '__init__')
    notify(self.initSpectrumWindowGroup, 'ccpnmr.Analysis.SpectrumWindowGroup', '__init__')
    notify(self.deletedSpectrumWindowGroup, 'ccpnmr.Analysis.SpectrumWindowGroup', 'delete')
    notify(self.setCPeaksText, 'ccpnmr.Analysis.AnalysisProject', 'setDoMinimalAnnotations')

    notify(self.setupCPeakCluster, 'ccp.nmr.Nmr.PeakCluster', '__init__')
    notify(self.removeCPeakCluster, 'ccp.nmr.Nmr.PeakCluster', 'delete')
    notify(self.setCPeakClusterPeaks, 'ccp.nmr.Nmr.PeakCluster', 'addPeak')
    notify(self.setCPeakClusterPeaks, 'ccp.nmr.Nmr.PeakCluster', 'removePeak')
    notify(self.setCPeakClusterPeaks, 'ccp.nmr.Nmr.PeakCluster', 'setPeaks')
    notify(self.setCPeakClusterText, 'ccp.nmr.Nmr.PeakCluster', 'setAnnotation')
    notify(self.setupCPeak, 'ccp.nmr.Nmr.Peak', '__init__')
    notify(self.removeCPeak, 'ccp.nmr.Nmr.Peak', 'delete')
    notify(self.setCPeakText, 'ccp.nmr.Nmr.Peak', 'setAnnotation')
    notify(self.setCPeakDimText, 'ccp.nmr.Nmr.PeakDim', 'setAnnotation')
    notify(self.setCPeakPosition, 'ccp.nmr.Nmr.PeakDim', 'setPosition')
    notify(self.setCPeakNumAliasing, 'ccp.nmr.Nmr.PeakDim', 'setNumAliasing')
    notify(self.setCPeakDimLineWidth, 'ccp.nmr.Nmr.PeakDim', 'setLineWidth')
    for func in ('__init__', 'setValue'):
      notify(self.setCPeakSize, 'ccp.nmr.Nmr.PeakIntensity', func)

    for clazz in ('DataSource', 'PeakList', 'Peak'):
      for func in ('add', 'remove', 'set'):
        notify(self.changedAppData, 'ccp.nmr.Nmr.%s' % clazz, '%sApplicationData' % func)

    notify(self.symbolDrawnPeakList, 'ccpnmr.Analysis.WindowPeakList', 'setIsSymbolDrawn')
    notify(self.textDrawnPeakList, 'ccpnmr.Analysis.WindowPeakList', 'setIsAnnotationDrawn')
    notify(MarkBasic.setPeakMarkColor, 'ccp.nmr.Nmr.PeakList', 'setPeakColor')
    notify(self.setupPeakAnnotation, 'ccp.nmr.Nmr.Peak', 'setFigOfMerit')
    notify(self.setupPeakAnnotation, 'ccp.nmr.Nmr.Peak', 'setDetails')

    notify(self.updateStripRegion, 'ccpnmr.Analysis.AxisRegion', 'setRegion')
    notify(self.updateTiedAxesAfter, 'ccpnmr.Analysis.AxisRegion', 'setRegion')
    notify(self.addedAxisRegion, 'ccpnmr.Analysis.AxisRegion', '__init__')
    notify(self.deletedAxisRegion, 'ccpnmr.Analysis.AxisRegion', 'delete')
    notify(self.initStoredContour, 'ccpnmr.Analysis.StoredContour', '__init__')
    notify(self.deletedStoredContour, 'ccpnmr.Analysis.StoredContour', 'delete')

    notify(self.changedPeakTextOffset, 'ccp.nmr.Nmr.Peak', 'setAppDataValue', keyword='textOffset')
    notify(self.changedPeakDimTextOffset, 'ccp.nmr.Nmr.PeakDim', 'setAppDataValue', keyword='textOffset')

    notify(self.changedPeakListSymbolColor, 'ccpnmr.Analysis.AnalysisPeakList', 'setSymbolColor')
    notify(self.changedPeakListSymbolStyle, 'ccpnmr.Analysis.AnalysisPeakList', 'setSymbolStyle')

    notify(self.changedUsePeakArrow, 'ccpnmr.Analysis.AnalysisSpectrum', 'setUsePeakArrow')

    notify(self.changedExpMeasurementType, 'ccp.nmr.Nmr.ExpDimRef', 'setMeasurementType')

    notify(self.newSpectrumWindowPaneAfter, 'ccpnmr.Analysis.SpectrumWindowPane', '__init__')
    notify(self.newSpectrumWindowViewAfter, 'ccpnmr.Analysis.SpectrumWindowView', '__init__')

  # in theory below could be done with __del__ but cannot
  # guarantee when garbage collection will happen
  def destroy(self, closingProject=False):

    self.visibleSpectra = {}
    self.toolbarSpectra = {}

    denotify = self.unregisterNotify
    denotify(self.initPeakList, 'ccp.nmr.Nmr.PeakList', '__init__')

    self.curateNotifiers(denotify)

    AssignmentBasic.removeAssignmentNotifiers()
    CouplingBasic.removeCouplingNotifiers()

    try:
      self.setBackupOff()
    except:
      pass

    self.deleteCMemory()

    if self.project:
      if closingProject:
        self.initProject(None)
      else:
        self.project = None
      self.nmrProject = self.analysisProject = None

  def deleteCMemory(self):

    if not self.project:
      return

    for expt in self.nmrProject.experiments:
      for spectrum in expt.dataSources:
        # NOTE:ED - causes fatal crash
        # self.deleteBlockFile(spectrum)

        for peakList in spectrum.peakLists:
          self.removeCPeakList(peakList)

    for window in self.analysisProject.spectrumWindows:
      for windowPane in window.spectrumWindowPanes:
        for view in windowPane.spectrumWindowViews:
          self.deleteContourSliceFiles(view)
          for winPeakList in view.windowPeakLists:
            self.removeCWinPeakList(winPeakList)

    removeCPeakCluster = self.removeCPeakCluster
    for peakCluster in self.nmrProject.peakClusters:
      if peakCluster.clusterType == MULTIPLET:
        removeCPeakCluster(peakCluster)

    if self.mem_cache: # true if C code working
      del self.mem_cache
      self.mem_cache = None

  def closeProject(self):

    self.argServer.destroy()
    self.argumentServer = self.argServer = None
    # NOTE:ED testing destruction error
    Analysis.destroy(self, closingProject=True)
    # self.destroy(closingProject=True)
    sys.stdout.flush()

  def updateStripRegion(self, axisRegion):

    axisRegionGroup = axisRegion.axisRegionGroup

    if not axisRegionGroup:
      return

    if axisRegionGroup in self.axisRegionGroups:
      return

    self.axisRegionGroups.append(axisRegionGroup)

    # seem to need after_idle, otherwise when zooming out at edge
    # of universe in one direction region in other direction
    # actually changes, so aspect ratio gets mucked up
    # (must be because this function is called by notifier)
    self.after_idle(lambda: self.doUpdateStripRegion(axisRegion))

  def doUpdateStripRegion(self, axisRegion):

    axisRegionGroup = axisRegion.axisRegionGroup
    newRegion = axisRegion.region
    try:
      for otherAxisRegion in axisRegionGroup.axisRegions:
        if otherAxisRegion is not axisRegion:
          origRegion = otherAxisRegion.region
          otherAxisRegion.region = newRegion

          # WindowPopup seems to be updating x and y axes
          # of own window ok  but not of other windows,
          # so do it here explicitly (not brilliant)
          axisPanel = otherAxisRegion.axisPanel
          label = axisPanel.label
          if label in X_Y:
            windowPane = axisPanel.spectrumWindowPane
            
            if label == 'x':
              otherLabel = 'y'
            else:
              otherLabel = 'x'
            
            otherAxisRegion = windowPane.findFirstAxisPanel(label=otherLabel).findFirstAxisRegion()
            r = Region1D(*(otherAxisRegion.region))
            s = abs(float(newRegion[1]-newRegion[0]) / (origRegion[1]-origRegion[0]))
            r.zoom(s)
            otherAxisRegion.region = (r[0], r[1])
    
    finally:
      self.axisRegionGroups.remove(axisRegionGroup)

  def updateTiedAxesAfter(self, axisRegion):

    axisPanel = axisRegion.axisPanel
    if len(axisPanel.axisRegions) != 1:
      return

    if self.changedAxisPanel:
      return

    if len(axisPanel.axisRegions) == 1:
      self.unregisterNotify(self.updateTiedAxesAfter, 'ccpnmr.Analysis.AxisRegion', 'setRegion')
      self.changedAxisPanel = axisPanel
      self.after_idle(self.updateTiedAxes)

  def updateTiedAxes(self):

    data = []
    axisPanel = self.changedAxisPanel
    axisRegion = axisPanel.findFirstAxisRegion()
    tiedAxisPanels = WindowBasic.getTiedAxisPanels(axisPanel)
    for tiedAxisPanel in tiedAxisPanels:
      if tiedAxisPanel != self.changedAxisPanel:
        data.append((axisRegion, tiedAxisPanel.findFirstAxisRegion()))

    for axisRegion, tiedAxisRegion in data:
      if axisRegion.axisPanel.label in ('x', 'y') and tiedAxisRegion.axisPanel.label in ('x', 'y'):
        # change entire region
        tiedAxisRegion.region = axisRegion.region
      else:
        # change only center point, not extent
        (r1, r2) = tiedAxisRegion.region
        dr = 0.5 * (r2 - r1)
        (r1, r2) = axisRegion.region
        c = 0.5 * (r1 + r2)
        tiedAxisRegion.region = (c - dr, c + dr)

    # if you don't do after_idle then it gets into an infinite loop
    self.after_idle(lambda: self.registerNotify(self.updateTiedAxesAfter, 'ccpnmr.Analysis.AxisRegion', 'setRegion'))
    self.changedAxisPanel = None
      
  def setupBackup(self):

    if not self.project:
      return

    if self.analysisProject.doAutoBackup:
      self.setBackupOn()
    else:
      self.setBackupOff()

  def setBackupOn(self):

    self.setBackupOff()

    repository = self.project.findFirstRepository(name='backup')
    if repository:
      directory = repository.url.path
      freq = self.analysisProject.autoBackupFreq
      
      if freq > 0:
        freq = 60000 * freq # convert from minutes to msec
        self.backup_alarm_id = self.after(freq, self.doBackupAuto)

  def setBackupOff(self):

    if self.backup_alarm_id is not None:
      self.after_cancel(self.backup_alarm_id)
      self.backup_alarm_id = None

  def doBackupAuto(self):

    if self.analysisProject.doAutoBackup:
      self.doBackup()
      freq = self.analysisProject.autoBackupFreq

      if freq > 0:
        freq = 60000 * freq # convert from minutes to msec
        self.backup_alarm_id = self.after(freq, self.doBackupAuto)

  def doBackup(self):

    try:
      backupProject(self.project)
      
    except Implementation.ApiError, e:
      print 'Backup error %s: %s' % (time.ctime(time.time()), e.error_msg)

  def setupPeakAnnotation(self, peak):

    PeakBasic.makePeakAnnotation(peak)

  def refreshPeakResNumberAfter(self, residue):

    self.after_idle( lambda r=residue: self.refreshPeakResNumber(r) )

  def refreshPeakResNumber(self, residue):

    resonances = set()
    chain = residue.chain
    for residue2 in chain.residues:
      for atom in residue2.atoms:
        atomSet = atom.atomSet
      
        if atomSet and atomSet.resonanceSets:
          for resonanceSet in atomSet.resonanceSets:
            for resonance in resonanceSet.resonances:
              resonances.add(resonance)

    peaks = set()
    for resonance in resonances:
      for contrib in resonance.peakDimContribs:
        peaks.add(contrib.peakDim.peak)

    for peak in peaks:
      PeakBasic.makePeakAnnotation(peak)

  def updateChemShifts(self, *object):

    if self.project:
      for shiftList in self.nmrProject.findAllMeasurementLists(className='ShiftList'):
        AssignmentBasic.updateAllShifts(shiftList)

  def setCPeakClusterPeaks(self, peakCluster):
 
    peaks = peakCluster.peaks
    cPeakCluster = hasattr(peakCluster, 'cPeakCluster') and peakCluster.cPeakCluster
    
    if cPeakCluster:
      if peaks:
        cPeaks = [peak.cPeak for peak in peaks]
        cPeakCluster.setPeaks(cPeaks)
      else:
        del peakCluster.cPeakCluster

  def setCPeakClusterText(self, peakCluster):
 
    if not hasattr(peakCluster, 'cPeakCluster'):
      self.setupCPeakCluster(peakCluster)
 
    if not (hasattr(peakCluster, 'cPeakCluster') and peakCluster.cPeakCluster):
      return

    cPeakCluster = peakCluster.cPeakCluster
    cPeakCluster.setText(peakCluster.annotation or '-')
    
    for i, text in enumerate(CouplingBasic.getCouplingAnnotations(peakCluster)):
      cPeakCluster.setDimText(i, text)

  def setCPeaksText(self, analysisProject):
    
    setCPeakText = self.setCPeakText
    nmrProject = self.nmrProject
    for experiment in nmrProject.experiments:
      for dataSource in experiment.dataSources:
        for peakList in dataSource.peakLists:
          for peak in peakList.peaks:
            setCPeakText(peak)

  def setCPeakText(self, peak):
 
    if not (hasattr(peak, 'cPeak') and peak.cPeak):
      return

    analysisProject = self.analysisProject
    if analysisProject.doMinimalAnnotations:
      text = PeakBasic.getSimplePeakAnnotation(peak, doChain=analysisProject.doChainAnnotations)
      
    else: # analysisProject.doCompressAnnotations:
      if len(peak.annotation or '') < 5:
        PeakBasic.makePeakAnnotation(peak)
    
      text = PeakBasic.getPeakAnnotation(peak, noPeakDimAnnotationChar='',
                                         joinChar='', doPeakDims=False)
    
    #text = PeakBasic.getPeakAnnotation(peak, noPeakDimAnnotationChar='-', joinChar=' ')
    
    peak.cPeak.setText(text)
 
  def setCPeakDimText(self, peakDim):

    self.setCPeakText(peakDim.peak)

  def setCPeakDimTextOffset(self, peakDim, offset):

    peak = peakDim.peak
    if not (hasattr(peak, 'cPeak') and peak.cPeak):
      return

    if offset is None:
      peak.cPeak.resetTextOffset(peakDim.dim-1)
    else:
      peak.cPeak.setTextOffset(peakDim.dim-1, offset)

  def setCPeakTextOffset(self, peak, offset):

    if not (hasattr(peak, 'cPeak') and peak.cPeak):
      return

    if offset is None:
      peak.cPeak.resetTextOffset(-1)
    else:
      peak.cPeak.setTextOffset(-1, offset)

  def setCPeakPosition(self, peakDim=None, peak=None, peakDims=None):

    if not peak:
      peak = peakDim.peak

    if not (hasattr(peak, 'cPeak') and peak.cPeak):
      return

    if not peakDims:
      peakDims = peak.sortedPeakDims()

    #if (len(peakDims) != peak.peakList.dataSource.numDim):
    #  return

    position = [ peakDim.position for peakDim in peakDims ]
    
    if None not in position:
      peak.cPeak.setPosition(position)
      
      if peakDim and peakDim.dataDim.className != 'FreqDataDim':
        AssignmentBasic.makePeakDimAnnotation(peakDim)
        self.setCPeakDimText(peakDim)

  def setCPeakNumAliasing(self, peakDim=None, peak=None, peakDims=None):

    if not peak:
      peak = peakDim.peak

    if not (hasattr(peak, 'cPeak') and peak.cPeak):
      return

    if not peakDims:
      peakDims = peak.sortedPeakDims()

    #if (len(peakDims) != peak.peakList.dataSource.numDim):
    #  return

    numAliasing = [ peakDim.numAliasing for peakDim in peakDims ]
    if None not in numAliasing:
      peak.cPeak.setNumAliasing(numAliasing)

  def setCPeakSize(self, peakIntensity):

    intensityType = peakIntensity.intensityType
    peak = peakIntensity.peak
    if not (hasattr(peak, 'cPeak') and peak.cPeak):
      return

    spectrum = peakIntensity.peak.peakList.dataSource
    scale = spectrum.scale
    # dividing by scale means that C world is working in data values
    # as given in the file rather than actual values after scaling

    if intensityType == 'height':
      peak.cPeak.setIntensity(peakIntensity.value/scale)
    
    elif intensityType == 'volume':
      peak.cPeak.setVolume(peakIntensity.value/scale)

  def setCPeakDimLineWidth(self, peakDim):

    peak = peakDim.peak
    if not hasattr(peak, 'cPeak'):
      return

    lineWidth = peakDim.lineWidth
    if lineWidth is None:
      lineWidth = 0

    peak.cPeak.setLineWidth(peakDim.dim-1, lineWidth)

  def setCPeakLineWidth(self, peak):

    for peakDim in peak.peakDims:
      self.setCPeakDimLineWidth(peakDim)

  def setupCPeak(self, peak):

    if not self.mem_cache:
      return

    # TBD: need to keep own copy of cPeakList since when peak.peakList
    # is being deleted then peak is deleted but peak.peakList
    # no longer valid by the time the notifier is called
    
    if not hasattr(peak, 'cPeakList'): # should not exist but protect against this
      peak.cPeakList = peak.peakList.cPeakList
    
    if not hasattr(peak, 'cPeak'): # should not exist but protect against this
      peak.cPeak = peak.cPeakList.addPeak()
    
    peakDims = peak.sortedPeakDims()
    self.setCPeakText(peak)
    self.setCPeakPosition(peak=peak, peakDims=peakDims)
    self.setCPeakNumAliasing(peak=peak, peakDims=peakDims)
    self.setCPeakLineWidth(peak)

    for peakIntensity in peak.peakIntensities:
      self.setCPeakSize(peakIntensity)

    offset = Util.getPeakTextOffset(peak)
    if offset is not None:
      self.setCPeakTextOffset(peak, offset)

    # make sure that peakDim.dataDimRef is set
    # TBD: perhaps this code should be in initPeakList()
    # but we already have peakDims here so take the opportunity
    # or should we be using PeakBasic.setupPeakDataDimRef()??
    for peakDim in peakDims:
      if not peakDim.dataDimRef:
        dataDim = peakDim.dataDim
        if dataDim.className == 'FreqDataDim':
          peakDim.dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)

  def removeCPeak(self, peak):

    #print 'removeCPeak1'
    if hasattr(peak,'cPeakList') and hasattr(peak,'cPeak'):
      try:
        if peak.cPeak.getIsSelected():
          try:
            self.currentPeaks.remove(peak)
          except:
            print 'C <-> Python peak mismatch. Selected C peak was not in currentPeaks.'

      except:
        pass

      peak.cPeakList.removePeak(peak.cPeak)

    if peak is self.currentPeak:
      self.currentPeak = None

    #print 'removeCPeak2'

  def setupCPeakCluster(self, peakCluster):
    
    if peakCluster.clusterType != MULTIPLET:
      return
    
    peaks = list(peakCluster.peaks)
    if not peaks:
      return
    
    if not self.mem_cache:
      return

    nDim = len(peaks[0].peakDims) 
    if not hasattr(peakCluster, 'cPeakCluster'):
      peakCluster.cPeakCluster = cPeakCluster = CPeakCluster(nDim, MULTIPLET_ENUM)
    else:
      cPeakCluster = peakCluster.cPeakCluster
    
    for peak in peaks:
      if hasattr(peak, 'cPeak'):
        cPeakCluster.addPeak(peak.cPeak)
    
    self.setCPeakClusterText(peakCluster)

  def setupCPeakList(self, peakList, initialising=False):

    if hasattr(peakList, 'cPeakList'): # should not exist but protect against this
      return

    if not self.mem_cache:
      return

    spectrum = peakList.dataSource
    
    # TBD: below assumes all dataDims are freqDataDims
    npoints = [ dataDim.numPoints for dataDim in spectrum.sortedDataDims() ]
    
    peakList.cPeakList = CPeakList(npoints)
    peaks = peakList.sortedPeaks()
    for peak in peaks:
      self.setupCPeak(peak)
      if initialising:
        for peakDim in peak.peakDims:
          offset = Util.getPeakDimTextOffset(peakDim)
          if offset is not None:
            self.setCPeakDimTextOffset(peakDim, offset)

    self.colorPeakList(peakList)
    self.symbolPeakList(peakList)

  def removeCPeakCluster(self, peakCluster):

    if hasattr(peakCluster,'cPeakCluster'):
      del peakCluster.cPeakCluster

  def removeCPeakList(self, peakList):

    #print 'removeCPeakList'
    for peak in peakList.peaks:
      if hasattr(peak, 'cPeakList'):
        del peak.cPeakList
        
    if hasattr(peakList,'cPeakList'):
      del peakList.cPeakList
      
    for peak in peakList.peaks:
      if hasattr(peak, 'cPeak'):
        del peak.cPeak

  def setupCWinPeakList(self, windowPeakList):

    peakList = windowPeakList.analysisPeakList.peakList
    if not hasattr(peakList, 'cPeakList'):
      return

    view = windowPeakList.spectrumWindowView
    
    valueAxis = WindowBasic.windowPaneHasValueAxis(view.spectrumWindowPane)
    windowPeakList.cWinPeakList = CWinPeakList(peakList.cPeakList, valueAxis)

    self.symbolDrawnPeakList(windowPeakList)
    self.textDrawnPeakList(windowPeakList)
    self.textPointerDrawnPeakList(windowPeakList)

  def removeCWinPeakList(self, winPeakList):

    #print 'removeCWinPeakList'
    if hasattr(winPeakList,'cWinPeakList'):
      del winPeakList.cWinPeakList

  def initMolSystems(self):

    for molSystem in self.project.molSystems:
      if not  Util.isTopObjectAnalysisInitialised(molSystem):
        for chain in molSystem.chains:
          self.initMolSystemChain(chain)

        Util.setTopObjectAnalysisSaveTime(molSystem)

  def pinitMolSystemChain(self, chain):

    import profile, pstats

    self.chainn = chain
    profile.run('top.initMolSystemChain(top.chainn)','initChainOut.txt')

  def initMolSystemChain(self, chain):

    analysisProject = self.analysisProject
    if (not analysisProject.chainMappings) or \
          (chain.residues and not hasattr(chain.findFirstResidue(), 'residueMapping')):
      
      aromaticEquiv = True
      
      """
      for residue in chain.residues:
        if residue.ccpCode in ('Phe','Tyr','Ptr'):
          if not residue.findFirstAtom().atomSet:
            msg = 'Are chain %s Phe/Tyr (Hd1,Hd2) and (He1,He2) atom sets '
            msg += 'equivalent due to rotation?'
            if self.showYesNo('Question',msg % (chain.code)):
              aromaticEquiv = True
            break
      """
      
      atomSetMappings = []
      getMapping = MoleculeBasic.getResidueMapping
      
      if self.argServer and self.argServer.inGui:
        progressBar = None
        if not analysisProject.findFirstChainMapping(chain=chain):
          text = 'Making residue mappings for %s chain %s'
          data = (chain.molSystem.code, chain.code)
          progressBar = self.argServer.getProgressBar(text=text % data,
                                                      total=len(chain.residues))


        for residue in chain.residues:
          residueMapping = getMapping(residue, aromaticsEquivalent=aromaticEquiv)
          atomSetMappings.extend( residueMapping.atomSetMappings )

          if progressBar:
            progressBar.increment()

      else:
        for residue in chain.residues:
          msg = "Making Atom Sets and Mappings for residue %s %s %d"
          print msg % (chain.code,residue.ccpCode,residue.seqCode)
          residueMapping = getMapping(residue, aromaticsEquivalent=aromaticEquiv)
          atomSetMappings.extend( residueMapping.atomSetMappings )

      atomSetDict = {}
      for atomSet in self.nmrProject.atomSets:
        atomSetDict[atomSet.serial] = atomSet

      for atomSetMapping in atomSetMappings:
        serials  = atomSetMapping.atomSetSerials
        atomSets = []
        for serial in serials:
          atomSet = atomSetDict.get(serial)
          
          if atomSet:
            atomSets.append(atomSet)
        
        if atomSets:
          AssignmentBasic.updateAtomSetMapping(atomSetMapping,atomSets=atomSets)
        else:
          atomSetMapping.delete()

  def deletedMolSystemChain(self, chain):

    msCode = chain.molSystem.code
    code   = chain.code
    findMappings = self.analysisProject.findAllChainMappings

    for chainMapping in findMappings(molSystemCode=msCode, chainCode=code):
      chainMapping.delete()

  def initAtomSetMappings(self):
    
    AssignmentBasic.initAtomSetMappings(self.analysisProject)


  def initShiftLists(self):

    shiftList = self.nmrProject.findFirstMeasurementList(className='ShiftList')
    if not shiftList:
      shiftList = ExperimentBasic.newShiftList(self.project, unit='ppm')

    self.shiftList = shiftList

  def initExperiment(self, experiment):

    ExperimentBasic.initExpTransfers(experiment, overwrite=False)
    if not self.shiftList:
      self.initShiftLists()

    if not experiment.findFirstExpDim(isAcquisition=True):
      if experiment.refExperiment:
        ExperimentBasic.setRefExperiment(experiment, experiment.refExperiment)

    if hasattr(experiment, 'shiftList'):
      if experiment.shiftList:
        return

    experiment.setShiftList(self.shiftList)

  def changedExpMeasurementType(self, expDimRef):
    
    # Need this for identifying DQ experiments, which need
    # DQ axis after load time 
    
    experiment = expDimRef.expDim.experiment
    
    spectra = ExperimentBasic.getExperimentSpectra(experiment)
    
    for spectrum in spectra:
      analysisSpectrum = Util.getAnalysisSpectrum(spectrum)
    
      for view in analysisSpectrum.spectrumWindowViews:
        view.delete()
    
      if hasattr(spectrum, 'isFinished') and spectrum.isFinished:
        self.checkCreateSpectrumViews(spectrum)

  def newSpectrumWindowPaneAfter(self, windowPane):

    self.after_idle(lambda: self.newSpectrumWindowPane(windowPane))
    
  def newSpectrumWindowPane(self, windowPane):

    WindowBasic.updateWindowPaneSliceRange(windowPane)

  def newSpectrumWindowViewAfter(self, view):

    self.after_idle(lambda: self.newSpectrumWindowView(view))
    
  def newSpectrumWindowView(self, view):

    windowPane = view.spectrumWindowPane
    if len(windowPane.spectrumWindowViews) == 1: # this is only view in pane
      WindowBasic.updateWindowPaneSliceRangeByView(view)

  def initSpectra(self):

    experiments = self.nmrProject.experiments
    for experiment in experiments:
      if not Util.isTopObjectAnalysisInitialised(experiment.nmrProject):
        self.initExperiment(experiment)
      
      for expDim in experiment.expDims:
        for expDimRef in expDim.expDimRefs:
          if expDimRef.measurementType == 'shift':
            expDimRef.measurementType = 'Shift'
          elif expDimRef.measurementType == 'JCoupling':
            expDimRef.measurementType = 'Jcoupling'
      
      spectra = ExperimentBasic.getExperimentSpectra(experiment)
      
      for spectrum in spectra:
        self.initAnalysisSpectrum(spectrum)
        self.initBlockFile(spectrum)
        self.initStoredContours(spectrum.analysisSpectrum)
        self.initSpectrumPeakList(spectrum)
        self.checkSpectrumContourRange(spectrum)
        
        for peakList in spectrum.peakLists:
          self.initPeakList(peakList, initialising=True)
          
        self.checkCreateSpectrumViews(spectrum)

  def checkSpectrumContourRange(self, spectrum):

    if self.analysisProject.contourToUnaliased:
      unitConvert = UnitConverter.unit_converter
    
      for dataDim in spectrum.dataDims:
        if (dataDim.className == 'FreqDataDim') and dataDim.dataDimRefs:
          dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
          expDimRef  = dataDimRef.expDimRef

          unit = expDimRef.unit
          if not unit or unit == 'None':
            unit = expDimRef.unit = 'ppm'
          valRange = [unitConvert[('point', unit)](1.0, dataDimRef),
                      unitConvert[('point', unit)](1.0+dataDim.numPoints, dataDimRef) ]
          valRange.sort()

          if expDimRef.minAliasedFreq is None:
            expDimRef.minAliasedFreq = valRange[0]

          if expDimRef.maxAliasedFreq is None:
            expDimRef.maxAliasedFreq = valRange[1]

  def initBlockFile(self, spectrum, updateContourLevels=False):

    self.deleteBlockFile(spectrum)

    writeable = hasattr(spectrum, 'writeable') and spectrum.writeable
    spectrum.block_file = None
    
    if not self.mem_cache:
      return

    if ExperimentBasic.isShapeSpectrum(spectrum):
      block_file = getShapeBlockFile(spectrum)
      if not block_file:
        return
      spectrum.block_file = block_file
    
    else:
      block_file = getBlockFile(spectrum,
                                mem_cache=self.mem_cache,
                                writeable=writeable)
      if not block_file:
        return

      spectrum.block_file = block_file
      try:
        block_file.open()
        if isWindowsOS(): # work around for memory trample - not complete fix 
          block_file.close()
      
      except:
        fileName = getDataSourceFileName(spectrum)
        msg = 'Warning: spectrum "%s": data file "%s" cannot be opened'
        print msg % (spectrum.name, fileName)
        self.deleteBlockFile(spectrum)

    if updateContourLevels:
      Util.defaultContourLevels(spectrum,
                                updateContourLevels=updateContourLevels)
      analysisSpectrum = spectrum.analysisSpectrum
      Util.updateSpectrumLevelParams(analysisSpectrum,
                                     analysisSpectrum.posLevels,
                                     analysisSpectrum.negLevels)

  def deleteBlockFile(self, spectrum):

    if not spectrum:
      return

    try:
      if hasattr(spectrum, 'block_file'):
        del spectrum.block_file
    except:
      pass

    spectrum.block_file = None

  def initStoreFile(self, storedContour):

    fileName = storedContour.fullPath
    dims = storedContour.dims
    xdim = dims[0] - 1
    ydim = dims[1] - 1

    spectrum = storedContour.analysisSpectrum.dataSource
    ndim = spectrum.numDim
    # TBD: this doesn't work for shape data
    blockSize = spectrum.dataStore.blockSizes

    storeFile = StoreFile and StoreFile.StoreFile(fileName, ndim, xdim, ydim, blockSize)

    return storeFile

  def initStoredContour(self, storedContour, doViews = True):

    try:
      storedContour.storeFile = self.initStoreFile(storedContour)
    except StoreFile.error, e:
      fileName = storedContour.fullPath
      msg = 'Contour file "%s": %s: delete corresponding stored contour?'
      
      if self.showYesNo('Stored contours', msg % (fileName, str(e))):
        try:
          storedContour.delete()
        except:
          pass
      
      return

    if doViews:
      analysisSpectrum = storedContour.analysisSpectrum
      for window in self.analysisProject.spectrumWindows:
        for windowPane in window.spectrumWindowPanes:
          view = windowPane.findFirstSpectrumWindowView(analysisSpectrum=analysisSpectrum)
          if view:
            self.setupViewStoredContourFile(storedContour, view)

  def deletedStoredContour(self, storedContour):

    try:
      del storedContour.storeFile
    except:
      pass

  def initStoredContours(self, analysisSpectrum):

    for storedContour in analysisSpectrum.storedContours:
      self.initStoredContour(storedContour, doViews=False)

  def initSpectrumViews(self):

    windows = self.analysisProject.spectrumWindows
    badViews = []
    
    for window in windows:
      for windowPane in window.spectrumWindowPanes:
        self.initSpectrumWindowPane(windowPane)
        
        views = windowPane.spectrumWindowViews
 
        for view in views:
          # first check if there is a WindowPeakList for every PeakList
          try:
            spectrum = view.analysisSpectrum.dataSource
          except:
            spectrum = None
 
          if spectrum:
            peakLists = list(spectrum.peakLists)
            # do not bother with those that already exist
            for winPeakList in view.windowPeakLists:
              try:
                peakList = winPeakList.analysisPeakList.peakList
                peakLists.remove(peakList)
              except:
                winPeakList.delete()
 
            # create a WindowPeakList for the remainder
            for peakList in peakLists:
              view.newWindowPeakList(analysisPeakList=peakList.analysisPeakList)

            # tell window that view axis mapping changed (to get slice initialised)
            self.changedViewAxisMapping(view)
          else:
            badViews.append(view)

    for view in badViews:
      try:
        view.delete()
      except:
        pass # a bit stuffed at this point, have view with no spectrum

  def initAnalysisSpectrum(self, spectrum):

    # TBD: it sets up contour levels in the below but
    # this is called before block_file is set up so
    # it should probably be removed
    if not spectrum.analysisSpectrum:
      Util.getAnalysisSpectrum(spectrum)
    
    analysisSpectrum = spectrum.analysisSpectrum
    repository = analysisSpectrum.root.findFirstRepository(name='userData')
    repositoryPath = repository.url.dataLocation
    contourDir = analysisSpectrum.contourDir
    if contourDir:
      oldContourPath = contourDir.dataLocation
      if not oldContourPath.startswith(repositoryPath): # slightly dangerous
        newContourPath = joinPath(repositoryPath, 'storedContours')
        if not os.path.exists(newContourPath) and os.path.exists(oldContourPath):
          shutil.copytree(oldContourPath, newContourPath)
        analysisSpectrum.contourDir = Url(path=newContourPath)

    else:
      path = joinPath(repositoryPath, 'storedContours')
      contourDir = analysisSpectrum.contourDir = Url(path=path)

    if not spectrum.activePeakList:

      for peakListB in spectrum.sortedPeakLists():
        # not sure how it can be deleted but sometimes it happens
        if peakListB.peaks and not peakListB.isDeleted:
          peakList = peakListB
          break
      else:
        peakList = spectrum.findFirstPeakList(isDeleted=False)
    
      spectrum.activePeakList = peakList  # could be None
 
  def initAnalysisPeakList(self, peakList):

    if not peakList.analysisPeakList:
      analysisPeakList = Util.getAnalysisPeakList(peakList)

  def initAnalysisDataDim(self, dataDim):

    if not dataDim.analysisDataDim:
      analysisDataDim = Util.getAnalysisDataDim(dataDim)
    
  def initSpectrum(self, spectrum):

    spectrum.isFinished = False

    if not ExperimentBasic.isSpectrum(spectrum):
      return

    self.initAnalysisSpectrum(spectrum)

    # use after_idle because required spectrum children not created when notify called

    self.checkSpectrumContourRange(spectrum)
    # below separated out into separate function below
    #self.parent.after_idle(lambda: self.initSpectrumPeakList(spectrum))
    #self.parent.after_idle(lambda: self.initBlockFile(spectrum))
    #self.parent.after_idle(lambda: self.checkCreateSpectrumViews(spectrum))

  # below has to be called explicitly
  # this is so that user opening spectrum can change details first
  def finishInitSpectrum(self, spectrum):

    if not ExperimentBasic.isSpectrum(spectrum):
      return
    
    self.initAnalysisSpectrum(spectrum)

    self.initBlockFile(spectrum, updateContourLevels=True)
    self.initStoredContours(spectrum.analysisSpectrum)
    self.initSpectrumPeakList(spectrum)
    self.checkCreateSpectrumViews(spectrum)

    spectrum.isFinished = True

  def deletedResonanceSet(self, resonanceSet):

    self.initResonanceSet(resonanceSet)

  def initResonanceSet(self, resonanceSet):

    residue = resonanceSet.findFirstAtomSet().findFirstAtom().residue
    residueMapping = MoleculeBasic.getResidueMapping(residue)

    serials = [atomSet.serial for atomSet in resonanceSet.atomSets]
    serials.sort()

    for atomSetMapping in residueMapping.atomSetMappings:
      serials2 = list(atomSetMapping.atomSetSerials)
      serials2.sort()

      if serials2 == serials:
        # unambiguous
        AssignmentBasic.updateAtomSetMapping(atomSetMapping, resonanceSet.atomSets)
        # if the resonance set is just deleted the atomSet to resonanceSet link
        # will be gone; resonances will not be found; serials will be set to []

      elif (atomSetMapping.mappingType == 'ambiguous') and (len(serials) == 1):
        if serials[0] in serials:
          AssignmentBasic.updateAtomSetMapping(atomSetMapping)

  def initResonances(self):

    resonances = self.nmrProject.resonances
    for resonance in resonances:
      # create spin systems
      # make assignNames
      AssignmentBasic.initResonance(resonance, doMerge=False)
      for shift in resonance.shifts:
        AssignmentBasic.setQuickShiftList(shift)

  def initSpinSystems(self):

    return

  def initSpectrumPeakList(self, spectrum):

    if not spectrum.peakLists:
      spectrum.newPeakList(details='Default list')
     
    if not spectrum.activePeakList:
      spectrum.activePeakList = spectrum.sortedPeakLists()[0]

  def updatedDataDimDetails(self, dataDim):

    self.updatedSpectrumFileDetails(dataDim.dataSource)

  def updatedNumericMatrixDetails(self, numericMatrix):

    for dataSource in numericMatrix.nmrDataSources:
      self.updatedSpectrumFileDetails(dataSource)

  def updatedSpectrumFileDetails(self, spectrum):

    if not ExperimentBasic.isSpectrum(spectrum):
      return

    self.initBlockFile(spectrum)

    views = WindowBasic.getSpectrumViews(spectrum)
    for view in views:
      self.initContourSliceFiles(view)

  def updatedDataUrl(self, dataUrl):

    spectra = ExperimentBasic.getSpectra(self.project)
    for spectrum in spectra:
      if spectrum.dataStore and spectrum.dataStore.dataUrl is dataUrl:
        self.updatedSpectrumFileDetails(spectrum)

  def updateValueAxisTypeRegion(self, spectrum, axisType):

    if axisType.name != 'value':
      return

    if spectrum.numDim != 1:
      return

    first = (0,)
    last = (spectrum.findFirstDataDim().numPoints,)

    if hasattr(spectrum, 'block_file') and spectrum.block_file:
      block_file = spectrum.block_file
    else:
      return

    w0 = block_file.minValue(first, last)
    w1 = block_file.maxValue(first, last)

    (r0, r1) = axisType.region

    if w0 < r0:
      r0 = w0 - 0.05*(w1-w0)  
      
    if w1 > r1:
      r1 = w1 + 0.05*(w1-w0)  

    axisType.region = (r0, r1)

  def checkCreateSpectrumViews(self, spectrum):

    createdView = None
    windows = self.analysisProject.spectrumWindows
    for window in windows:
      for windowPane in window.spectrumWindowPanes:
        createdView = self.checkCreateSpectrumView(windowPane, spectrum) or createdView

    if not createdView:
      # there is no matching window already so create one
      axisTypes = spectrum.numDim * [None]
      for dataDim in spectrum.sortedDataDims():
        axisType = Util.findAxisTypeMatch(dataDim)
        if not axisType:
          break
        #print 'checkCreateSpectrumViews', dataDim.dim, axisType.name
        axisTypes[dataDim.dim-1] = axisType

      if None not in axisTypes:
        name = WindowBasic.defaultWindowName(self.project)
        # NOTE: the below relies on last dimension in pseudo-2D experiment being the sampled dim
        if spectrum.numDim == 1 or (spectrum.numDim == 2 and axisType.name == 'sampled'):
          # add value axis
          axisType = self.analysisProject.findFirstAxisType(name='value')
          self.updateValueAxisTypeRegion(spectrum, axisType)
          axisTypes[1:1] = [axisType]
        window = WindowBasic.createSpectrumWindow(self.project, name, [axisTypes,],
                                                       spectrum=spectrum)
        self.checkCreateSpectrumView(window.findFirstSpectrumWindowPane(), spectrum)

  def initAxisPanel(self, axisPanel):
    
    # use after_idle because all required information not set up when notify called
    self.after_idle(lambda: self.checkCreateSpectraViews(axisPanel.spectrumWindowPane))
    
  def deletedAxisPanel(self, axisPanel):

    try:
      if axisPanel.spectrumWindowPane:
        self.checkDeleteSpectraViews(axisPanel.spectrumWindowPane)
    
    except:
      pass

  def checkCreateSpectraViews(self, windowPane):

    spectra = ExperimentBasic.getSpectra(self.project)
    for view in windowPane.spectrumWindowViews:
      # remove ones already mapped
      spectra.remove(view.analysisSpectrum.dataSource)

    for spectrum in spectra:
      self.checkCreateSpectrumView(windowPane, spectrum)

    hasValueAxis = WindowBasic.windowPaneHasValueAxis(windowPane)
    views = windowPane.spectrumWindowViews
    if views and hasattr(windowPane, 'hasDefaultRegions'):
      del windowPane.hasDefaultRegions
      for axisPanel in windowPane.sortedAxisPanels():
        if axisPanel.axisType.isSampled:
          continue
        label = axisPanel.label
        if hasValueAxis and label == 'y':
          continue
        axisUnit = axisPanel.axisUnit
        converter = UnitConverter.unit_converter[('point', axisUnit.unit)]
        rmin = rmax = None
        for view in views:
          analysisSpectrum = view.analysisSpectrum
          spectrum = analysisSpectrum.dataSource
          axisMapping = view.findFirstAxisMapping(label=label)
          if not axisMapping:
            continue
          analysisDataDim = axisMapping.analysisDataDim
          dataDim = analysisDataDim.dataDim
          numPoints = dataDim.numPoints

          if label in X_Y:
            r0 = 0.5
            r1 = dataDim.numPoints + 0.5
          else:
            p = int(dataDim.numPoints / 2)
            r0 = p - 0.5
            r1 = p + 0.5

          dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
          region = [ converter(r, dataDimRef) for r in (r0, r1) ]
          (r0, r1) = Util.checkSwapRegion(region, axisUnit)
          if rmin is not None:
            rmin = min(rmin, r0)
            rmax = max(rmax, r1)
          else:
            rmin = r0
            rmax = r1
          if label not in X_Y:
            break  # just use first spectrum to set orthogonal region
        if rmin is not None:
          d = rmax - rmin
          dmax = 20.0  # arbitrary
          if d > dmax:
            c = 0.5 * (rmin + rmax)
            rmin = c - 0.5 * dmax
            rmax = c + 0.5 * dmax
          for axisRegion in axisPanel.axisRegions:
            axisRegion.region = (rmin, rmax)

  # axisPanelMapping maps dataDim to axisPanel of window
  def createSpectrumView(self, windowPane, spectrum, axisPanelMapping):

    isVisible = self.visibleSpectra.get(spectrum, True)

    if isVisible is 1: # This deliberate: Map spectrum to first window view only.
      isVisible = True
      setVisibleFalse = True
    
    else:
      setVisibleFalse = False

    isInToolbar = self.toolbarSpectra.get(spectrum, True)

    if spectrum.numDim == 1:
      axisType = windowPane.topObject.findFirstAxisType(name='value')
      self.updateValueAxisTypeRegion(spectrum, axisType)

    view = Util.createSpectrumWindowView(windowPane, spectrum, axisPanelMapping,
                                         isPosVisible=isVisible, isNegVisible=isVisible,
                                         isInToolbar=isInToolbar, isSliceVisible=isVisible)
    if view:
      self.changedViewAxisMapping(view)
      self.setWindowRegions(view)

    # need to set this here, not above, because initSpectrumWindowView() uses
    # self.visibleSpectra to determine value of winPeakList.isSymbolDrawn
    if setVisibleFalse:
      self.visibleSpectra[spectrum] = False

    return view

  def checkCreateSpectrumView(self, windowPane, spectrum):

    if WindowBasic.isSpectrumInWindowPane(windowPane, spectrum):
      return True

    """ Ah, also need ordinary ND-reconstructed spectrum, so do that first
    # for shape spectra do not try to put into some other existing window
    if ExperimentBasic.isShapeSpectrum(spectrum):
      return False
    """

    try:
      # sometimes can get exception thrown because called too soon
      # (probably because idle not working) so try this ploy
      mapping = self.determineMapping(spectrum, windowPane)
    except:
      self.after_idle(lambda: self.checkCreateSpectrumView(windowPane, spectrum))
      return

    if mapping:
      self.createSpectrumView(windowPane, spectrum, mapping)
      return True
      
    else:
      return False

  def checkDeleteSpectraViews(self, windowPane):

    labels = [ ap.label for ap in windowPane.axisPanels if not ap.isDeleted ]

    for view in windowPane.spectrumWindowViews:
      for axisMapping in view.axisMappings:
        if axisMapping.label not in labels:
          try:
            dataDim = axisMapping.analysisDataDim.dataDim
          except:
            dataDim = None
          axisMapping.delete()
          if not dataDim or not self.findNewMapping(view, dataDim):
            view.delete()
          break

  def changedViewAxisMapping(self, view):

    self.initContourSliceFiles(view)

  def initContourSliceFiles(self, view):

    self.deleteContourSliceFiles(view)

    if not WindowBasic.windowPaneHasValueAxis(view.spectrumWindowPane):
      self.initContourFile(view)
      self.initStoredContourFiles(view)

    # this is needed for hasValueAxis windows because main window uses this
    self.initSliceFile(view)

  def deleteContourSliceFiles(self, view):

    if hasattr(view, 'contourFile') and view.contourFile:
      del view.contourFile

    if hasattr(view, 'storedContourFiles'):
      del view.storedContourFiles

    if hasattr(view, 'sliceFile'):
      labels = view.sliceFile.keys()
      for label in labels:
        del view.sliceFile[label]

    view.sliceFile = {}

  def initContourFile(self, view):

    view.contourFile = None

    dataSource = view.analysisSpectrum.dataSource
    if dataSource.numDim < 2:
      return

    xMapping = view.findFirstAxisMapping(label='x')
    if not xMapping:
      return

    yMapping = view.findFirstAxisMapping(label='y')
    if not yMapping:
      return
    
    xdim = xMapping.analysisDataDim.dataDim.dim - 1
    ydim = yMapping.analysisDataDim.dataDim.dim - 1

    if not isWindowsOS():
      # NOTE:ED - this is needed because there are threading issues with Windows
      if hasattr(dataSource, 'block_file') and dataSource.block_file:
        view.contourFile = ContourFile(xdim, ydim, dataSource.block_file, self.mem_cache)

  def initStoredContourFiles(self, view):

    view.storedContourFiles = []

    analysisSpectrum = view.analysisSpectrum
    for storedContour in analysisSpectrum.storedContours:
      self.setupViewStoredContourFile(storedContour, view)

  def setupViewStoredContourFile(self, storedContour, view):

    if not hasattr(storedContour, 'storeFile'):  # should always be the case
      return

    axisMappings = view.axisMappings
    dimList = [(am.label, am.analysisDataDim.dataDim.dim) for am in axisMappings]
    
    xdim = [ dim for label, dim in dimList if label == 'x' ][0]
    ydim = [ dim for label, dim in dimList if label == 'y' ][0]
    view_dims = (xdim, ydim)

    dims = tuple(storedContour.dims)

    if dims == view_dims:
      transposed = 0
    elif (dims[1], dims[0]) == view_dims:
      transposed = 1
    else: # incompatible contour file, so ignore
      return

    storeFile = storedContour.storeFile
    if self.mem_cache:
      storedContourFile = StoredContourFile(storeFile, self.mem_cache, transposed)
      view.storedContourFiles.append(storedContourFile)

  def initSliceFile(self, view):

    for label in X_Y:
      axisMapping = view.findFirstAxisMapping(label=label)
      if axisMapping:
        self.initAxisMapping(axisMapping)

  def initAxisMapping(self, axisMapping):

    view = axisMapping.spectrumWindowView
    dataSource = view.analysisSpectrum.dataSource
    if not hasattr(dataSource, 'block_file') or not dataSource.block_file:
      return

    block_file = dataSource.block_file
    if axisMapping.label == 'x':
      try:
        xdim = axisMapping.analysisDataDim.dataDim.dim - 1
        # first argument = 1 for horizontal, 0 for vertical
        if isWindowsOS():
          view._createSliceFileObjectX = (xdim, block_file, self.mem_cache)
        else:
          view.sliceFile['x'] = SliceFile(1, xdim, block_file, self.mem_cache)
      except:
        pass
    elif axisMapping.label == 'y' and not WindowBasic.windowPaneHasValueAxis(view.spectrumWindowPane):
      try:
        ydim = axisMapping.analysisDataDim.dataDim.dim - 1
        # first argument = 1 for horizontal, 0 for vertical
        if isWindowsOS():
          view._createSliceFileObjectY = (ydim, block_file, self.mem_cache)
        else:
          view.sliceFile['y'] = SliceFile(0, ydim, block_file, self.mem_cache)
      except:
        pass

  # below does not seem to be doing anything any more
  def setWindowRegions(self, spectrumView):

    # set regions if they have not been set already

    for axisMapping in spectrumView.axisMappings:

      label = axisMapping.label
      axisPanel = spectrumView.spectrumWindowPane.findFirstAxisPanel(label=label)

      if (not hasattr(axisPanel, 'firstMapping')):
        continue

      del axisPanel.firstMapping

      """ commented out 12 Jan 2006 while testing sampled dataDims
      # NBNB obsolete - uses findFirstDataDimRef
      try:
        dataDim = axisMapping.analysisDataDim.dataDim
      except:
        continue

      axisUnit = axisPanel.axisUnit
      converter = UnitConverter.unit_converter[('point', axisUnit.unit)]

      if (label in X_Y):
        r0 = 0.5
        r1 = dataDim.numPoints + 0.5
      else:
        p = int(dataDim.numPoints / 2)
        r0 = p - 0.5
        r1 = p + 0.5

      # TBD: more general dataDimRefs
      region = [ converter(r, dataDim.findFirstDataDimRef()) for r in (r0, r1) ]
      region = Util.checkSwapRegion(region, axisUnit)

      #for axisRegion in axisPanel.axisRegions:
      #  axisRegion.region = region
      """

  # TBD: lots of assumptions in determineMapping()
  def determineMapping(self, spectrum, windowPane):

    mapping = {}
    panels = list(windowPane.sortedAxisPanels())
    hasValueAxis = WindowBasic.windowPaneHasValueAxis(windowPane)
    if hasValueAxis:
      del panels[1]
    
    for dataDim in spectrum.sortedDataDims():

      for panel in panels:
        if Util.haveTypeMatch(panel, dataDim):
          break
          
      else:
        #print 'determineMapping failed at dataDim', dataDim.dim
        return None # no mapping possible

      mapping[dataDim] = panel
      panels.remove(panel)

    if (spectrum.numDim > 1):
      # in > 1D require both x and y mapped
      if panels and (panels[0].label in X_Y): # always want x,y mapped
        return None
        
    elif not hasValueAxis or (panel.label not in X_Y): # in 1D require x or y mapped
      return None

    return mapping

  def findNewMapping(self, view, dataDim):

    labels = [ axisMapping.label for axisMapping in view.axisMappings ]
    panels = [panel for panel in view.spectrumWindowPane.sortedAxisPanels() if not panel.isDeleted]
    if WindowBasic.windowPaneHasValueAxis(view.spectrumWindowPane):
      del panels[1]
      
    remaining = []
    for panel in panels:
      if panel.label not in labels:
        remaining.append(panel)

    for panel in remaining:
      if Util.haveTypeMatch(panel, dataDim):
        break
    else:
      return False

    analysisDataDim = Util.getAnalysisDataDim(dataDim)
    view.newAxisMapping(label=panel.label, analysisDataDim=analysisDataDim)
    self.changedViewAxisMapping(view)

    return True

  def initPeakClusters(self):
  
    setupCPeakCluster = self.setupCPeakCluster
    for peakCluster in self.nmrProject.peakClusters:
      if peakCluster.clusterType == MULTIPLET:
        setupCPeakCluster(peakCluster)

  def initPeakList(self, peakList, initialising=False):

    project = peakList.root
    analysisProject = self.analysisProject
    spectrum = peakList.dataSource

    self.initAnalysisPeakList(peakList)

    # Moved up from below because
    # AssignmentBasic.makePeakDimAnnotation(peakDim)
    # needs cPeaks for notifier to work
    self.setupCPeakList(peakList, initialising)

    if not Util.isTopObjectAnalysisInitialised(peakList.topObject):

      fitMethod = Util.getMethod(project, task='fit peak center',
                                 procedure='peak.cPeak.fitCenter',
                                 parameters=(('method', 'parabolic'),))
      manMethod = Util.getMethod(project, task='set peak intensity',
                                 procedure='manual')
      volMethod = Util.getMethod(project, task='fit peak volume',
                                 procedure='peak.cPeak.fitVolume',
                                 parameters=(('method', 'gaussian3'),))
      hytMethod = Util.getMethod(project,task='find peak height',
                                 procedure='peak.cPeak.getIntensity')

      if not ExperimentBasic.isSpectrum(spectrum):
        return

      peaks = peakList.peaks

      if peaks:
        dimRegionDict = {}
        peak = peakList.findFirstPeak()
        for peakDim in peak.peakDims:
          dataDimRef = peakDim.dataDimRef
          if dataDimRef:
            dimRegionDict[dataDimRef] = (peakDim.value,peakDim.value)

      molSystems = {}

      for peak in peaks:
        # make sure peak dims are annotated if assigned
        #if not peak.annotation:
        #  PeakBasic.makePeakAnnotation(peak)

        if peak.figOfMerit == 0.0 and not initialising:
          # If initialising we average shifts elsewhere. Do not know why the figOfMerit
          AssignmentBasic.updatePeakShifts(peak)

        for peakDim in peak.peakDims:
          dataDimRef = peakDim.dataDimRef

          if dataDimRef and dimRegionDict.get(dataDimRef):
            (minPos,maxPos) = dimRegionDict[dataDimRef]
            if peakDim.value > maxPos:
              dimRegionDict[dataDimRef] = (minPos,peakDim.value)
            elif peakDim.value < minPos:
              dimRegionDict[dataDimRef] = (peakDim.value,maxPos)

            for contrib in peakDim.peakDimContribs:
              resonanceSet = contrib.resonance.resonanceSet

              if resonanceSet:
                atom = resonanceSet.findFirstAtomSet().findFirstAtom()
                molSystems[atom.topObject] = True

          if not peakDim.annotation:
            if len(peakDim.peakDimContribs) > 0:
              AssignmentBasic.makePeakDimAnnotation(peakDim)
            elif not peakDim.dataDimRef:
              AssignmentBasic.makePeakDimAnnotation(peakDim)

        """ reinstate if I check the software?
        # get rid of extra methods
        if peak.fitMethod is not fitMethod:
          m = peak.fitMethod
          peak.fitMethod = fitMethod
          if m:
            m.delete()

        if not peak.findFirstPeakIntensity(method=hytMethod):
          height =  peak.findFirstPeakIntensity(intensityType = 'height')
          if height and height.method is not manMethod:
            peakIntensity = peak.newPeakIntensity(method=hytMethod, intensityType='height',
                                                  value=height.value, error=height.error)
            if len(height.method.peakIntensities) == 1:
              height.method.delete()
            else:
              height.delete()

        if not peak.findFirstPeakIntensity(method=volMethod):
          volume =  peak.findFirstPeakIntensity(intensityType = 'volume')
          if volume and volume.method is not manMethod:
            peakIntensity = peak.newPeakIntensity(method=volMethod, intensityType='volume',
                                                  value=volume.value, error=volume.error)
            if len(volume.method.peakIntensities) == 1:
              volume.method.delete()
            else:
              volume.delete()
        """

      expMolSystems = spectrum.experiment.molSystems
      for molSystem in molSystems.keys():
        if molSystem not in expMolSystems:
          spectrum.experiment.addMolSystem(molSystem)

      if analysisProject.contourToUnaliased:
        peak = peakList.findFirstPeak()
        if peak:
          for peakDim in peak.peakDims:
            if peakDim.dataDimRef: # only exists for freqDataDim
              (minPos,maxPos) = dimRegionDict[peakDim.dataDimRef]
              maxAliasedFreq  = peakDim.dataDimRef.expDimRef.maxAliasedFreq
              minAliasedFreq  = peakDim.dataDimRef.expDimRef.minAliasedFreq

              if (minAliasedFreq is None) or (minPos < minAliasedFreq):
                peakDim.dataDimRef.expDimRef.minAliasedFreq = minPos

              if (maxAliasedFreq is None) or (maxPos > maxAliasedFreq):
                peakDim.dataDimRef.expDimRef.maxAliasedFreq = maxPos

    windows = analysisProject.spectrumWindows
    for window in windows:
      for windowPane in window.spectrumWindowPanes:
        view = WindowBasic.getSpectrumWindowView(windowPane, spectrum)
        if view: # exists only if spectrum is mapped onto window
          analysisPeakList = peakList.analysisPeakList
          windowPeakList = view.findFirstWindowPeakList(analysisPeakList=analysisPeakList)
 
          if windowPeakList:
            self.initWinPeakList(windowPeakList)
          else:
            if view.isPosVisible or view.isNegVisible:
              isDrawn = True
            else:
              isDrawn = False

            windowPeakList = view.newWindowPeakList(analysisPeakList=analysisPeakList,
                                                    isSymbolDrawn=isDrawn,
                                                    isAnnotationDrawn=isDrawn)
            if initialising: # otherwise called by notifier
              self.initWinPeakList(windowPeakList)


  def deletedPeakCluster(self, peakCluster):
    
    self.removeCPeakCluster(peakCluster)
  
  def deletedPeakList(self, peakList):

    analysisPeakList = peakList.analysisPeakList
    
    # TBD: Do we need this? Gone automatically due to link? 
    
    if analysisPeakList:
      for winPeakList in analysisPeakList.windowPeakLists:
        if not winPeakList.isDeleted:
          winPeakList.delete()

    self.removeCPeakList(peakList)

  def initWinPeakList(self, winPeakList):

    self.setupCWinPeakList(winPeakList)

  def deletedWinPeakList(self, winPeakList):

    self.removeCWinPeakList(winPeakList)

  def initSpectrumWindow(self, spectrumWindow):

    group = self.analysisProject.activeWindowGroup
    if group:
      group.addSpectrumWindow(spectrumWindow)

  def initSpectrumWindowPane(self, windowPane, windowFrame=None):
    """
    Expects to be overwritten in AnalysisPopup
    """      
    pass

  def initSpectrumWindowView(self, view):

    analysisSpectrum = view.analysisSpectrum

    for analysisPeakList in analysisSpectrum.analysisPeakLists:
      isDrawn = self.visibleSpectra.get(analysisSpectrum.dataSource, True)
      view.newWindowPeakList(analysisPeakList=analysisPeakList,
                             isSymbolDrawn=isDrawn,
                             isAnnotationDrawn=isDrawn)

  def initSpectrumWindowGroup(self, spectrumWindowGroup):

    analysisProject = self.analysisProject
    groups = analysisProject.spectrumWindowGroups
    
    if len(groups) == 1: # hence groups[0] == spectrumWindowGroup
      analysisProject.activeWindowGroup = spectrumWindowGroup

  def deletedSpectrumWindowGroup(self, spectrumWindowGroup):

    analysisProject = self.analysisProject
    
    try:
      group = analysisProject.activeWindowGroup
    except:
      group = None
    
    if (not group) or group is spectrumWindowGroup:
      windowGroups = analysisProject.sortedSpectrumWindowGroups()
      
      if windowGroups:
        group = windowGroups[0]
      else:
        group = None
      
      analysisProject.activeWindowGroup = group

  def changedPeakTextOffset(self, owner):

    if owner.isDeleted:
      return

    appData = owner.findFirstApplicationData(application=self.application.name,
                                             keyword='textOffset')
    if appData:
      value = appData.value
    else:
      value = None

    self.setCPeakTextOffset(owner, value)

  def changedPeakDimTextOffset(self, owner):

    if owner.isDeleted:
      return

    appData = owner.findFirstApplicationData(application=self.application.name,
                                             keyword='textOffset')
    if appData:
      value = appData.value
    else:
      value = None

    self.setCPeakDimTextOffset(owner, value)

  def changedAppData(self, owner):
    
    # NBNMB TBD completely wrong and partially reworked. For later
    return

    appData = owner.applicationData
    if appData.application != owner.root.application.name:
      return

    keyword = appData.keyword
    
    if keyword in ('projectAutoBackup', 'projectAutoBackupFreq', 'projectBackupDirectory'):
      self.setupBackup()
    
    elif isinstance(owner, Nmr.PeakDim):
      if keyword == 'textOffset':
        self.setCPeakDimTextOffset(owner, appData.value)
    
    elif isinstance(owner, Nmr.Peak):
      if keyword == 'textOffset':
        self.setCPeakTextOffset(owner, appData.value)
    
    elif isinstance(owner, Nmr.PeakList):
      if keyword == 'color':
        self.colorPeakList(owner)
      elif keyword == 'symbol':
        self.symbolPeakList(owner)
    
    elif isinstance(owner, Nmr.DataSource):
      if keyword == 'isSpectrumPointerShown':
        for window in self.analysisProject.spectrumWindows:
          for windowPane in window.spectrumWindowPanes:
            for view in windowPane.spectrumWindowViews:
              if view.analysisSpectrum.dataSource is owner:
                for winPeakList in view.windowPeakLists:
                  self.textPointerDrawnPeakList(winPeakList)
 


  def deletedAppData(self, appData):

    if appData.application != appData.root.application.name:
      return

    owner = appData.owner

    if owner.isDeleted:
      return

    keyword = appData.keyword

    if isinstance(owner, Nmr.PeakDim):
      if keyword == 'textOffset':
        self.setCPeakDimTextOffset(owner, None)
        
    elif isinstance(owner, Nmr.Peak):
      if keyword == 'textOffset':
        self.setCPeakTextOffset(owner, None)

  def changedPeakListSymbolColor(self, analysisPeakList):

    color = analysisPeakList.symbolColor
    analysisPeakList.peakList.cPeakList.setColor(hexToRgb(color))

  def colorPeakList(self, peakList):

    self.changedPeakListSymbolColor(Util.getAnalysisPeakList(peakList))

  def changedPeakListSymbolStyle(self, analysisPeakList):

    symbolList = ['x', '+', 'o', '*']
    try:
      symbolCode = symbolList.index(analysisPeakList.symbolStyle)
    except:
      symbolCode = 0
      
    analysisPeakList.peakList.cPeakList.setSymbol(symbolCode)

  def symbolPeakList(self, peakList):

    self.changedPeakListSymbolStyle(Util.getAnalysisPeakList(peakList))
    
  def symbolDrawnPeakList(self, windowPeakList):

    #print 'symbolDrawnPeakList', windowPeakList.spectrumWindowView.spectrumWindow.name, \
    #      windowPeakList.spectrumWindowView.analysisSpectrum.dataSource.experiment.name, \
    #      windowPeakList.spectrumWindowView.analysisSpectrum.dataSource.name, \
    #      windowPeakList.isSymbolDrawn

    windowPeakList.cWinPeakList.setIsSymbolDrawn(windowPeakList.isSymbolDrawn)
    # turn annotation off if symbol off
    ### 10 Nov 06: no longer do this
    ###if (not windowPeakList.isSymbolDrawn):
    ###  windowPeakList.isAnnotationDrawn = False

  def textDrawnPeakList(self, windowPeakList):

    windowPeakList.cWinPeakList.setIsTextDrawn(windowPeakList.isAnnotationDrawn)
    # turn symbol on if annotation on
    ### 10 Nov 06: no longer do this
    ###if (windowPeakList.isAnnotationDrawn):
    ###  windowPeakList.isSymbolDrawn = True

  def textPointerDrawnPeakList(self, windowPeakList, analysisSpectrum=None):
 
    if not analysisSpectrum:
      analysisSpectrum = windowPeakList.analysisPeakList.analysisSpectrum

    windowPeakList.cWinPeakList.setIsTextPointerDrawn(analysisSpectrum.usePeakArrow)
 
  def changedUsePeakArrow(self, analysisSpectrum):

    for analysisPeakList in analysisSpectrum.analysisPeakLists:
      windowPeakLists = analysisPeakList.windowPeakLists
      for windowPeakList in windowPeakLists:
        self.textPointerDrawnPeakList(windowPeakList, analysisSpectrum)

  def getApplication(self):

    return self.application

  def getProject(self):

    return self.project

  def getExperiments(self):

    if self.project:
      return self.nmrProject.sortedExperiments()
      
    else:
      return []

  def getSpectra(self):

    spectra = []
    experiments = self.getExperiments()
    for experiment in experiments:
      spectra.extend(experiment.sortedDataSources())

    return spectra

  def getMolSystems(self):

    if self.project:
      return self.project.sortedMolSystems()
      
    else:
      return ()

  def getColors(self):

    # TBD: BROKEN
    if self.project:
      profile = self.project.currentAnalysisProfile
      schemes = [cs for cs in profile.colorSchemes if len(cs.colors) == 1]
      colors = [(hexToHsb(cs.colors[0]), (cs.name, cs.colors[0])) for cs in schemes]
      colors.sort()
      
      return [c[1] for c in colors]
    else:
      return ()

  def getColorNames(self):

    colors = self.getColors()
    names = [ x[0] for x in colors ]

    return names

  def getColorSchemes(self):

    if self.project:
      return self.analysisProject.colorSchemes
      
    else:
      return ()

  def getColorSchemeNames(self):

    schemes = self.getColorSchemes()
    names = [ scheme.name for scheme in schemes ]

    return names

  def getPanelTypes(self):

    if self.project:
      return self.analysisProject.sortedPanelTypes()
      
    else:
      return []

  def getAxisTypes(self):

    if self.project:
      axisTypes = self.analysisProject.sortedAxisTypes()
      data1 = []
      data2 = []
      for axisType in axisTypes:
        isotope = axisType.findFirstIsotope()
        
        if isotope:
          data1.append((isotope.massNumber, axisType))
        else:
          data2.append((axisType.name, axisType))
          
      data1.sort()
      data2.sort()
      axisTypes= [x[1] for x in data1] + [x[1] for x in data2]
      return axisTypes
      
    else:
      return []

  def getAxisUnits(self):

    if self.project:
      return self.analysisProject.sortedAxisUnits()
      
    else:
      return []

  def getMeasurementTypes(self):

    # TBD: more general
    return ['shift', 'temperature', 'time']

  def getWindows(self):

    if self.project:
      return self.analysisProject.sortedSpectrumWindows()
      
    else:
      return []

  def getActiveWindows(self):

    if self.project:
      group = self.analysisProject.activeWindowGroup
      if group:
        return group.sortedSpectrumWindows()

    return []

  # can be overridden in superclass, but must exist here since called from initProfile
  def initMacros(self):

    pass

  def initProfile(self, obj=None):
    # Called at load and when changing current analysisProfile
    
    profile = self.analysisProfile
    
    if not profile: # is this ever true?
      return
    
    self.initMacros()

    # Change spectrum colors
    
    
    
    # Change background colours
  
  
  
    # Setup the quick lookup residue user codes dict from residue profiles
    # & backward compatibility to convert appData into ResidueProfiles
    
    if self.analysisProject and not profile.residueProfiles:
      application = self.project.application
      data = application.getValue(self.analysisProject,
                                  keyword='customResidueCodesDict',
                                  defaultValue=None)
  
      if data:
        resDict = cPickle.loads(data)
      else:
        resDict = MoleculeBasic.userResidueCodesDict
      
      ffcc = self.project.findFirstChemComp
      nrp = profile.newResidueProfile
      for molType in resDict.keys():
        for ccpCode in resDict[molType].keys():
          if not ffcc(molType=molType, ccpCode=ccpCode):
            continue
          
          guiName = resDict[molType][ccpCode]
          nrp(molType=molType, ccpCode=ccpCode, guiName=guiName)
     
      
    resDict = {}
    for resProfile in profile.residueProfiles:
      molType = resProfile.molType
      if not resDict.has_key(molType):
        resDict[molType] = {}
    
      resDict[molType][resProfile.ccpCode] = resProfile.guiName
     
    MoleculeBasic.setUserResidueCodesDict(resDict)

  def cacheResidueCode(self, residueProfile):
   
    resDict = MoleculeBasic.userResidueCodesDict
   
    if residueProfile.isDeleted:
      del resDict[residueProfile.molType][residueProfile.ccpCode]
      
    else:
      molType = residueProfile.molType
      if not resDict.has_key(molType):
        resDict[molType] = {}
        
      resDict[molType][residueProfile.ccpCode] = residueProfile.guiName

  def create1dShape(self, dataDim):

    newSpectrum = create1dShape(dataDim, applicationName=self.application.name)
    self.finishInitSpectrum(newSpectrum)

    return newSpectrum


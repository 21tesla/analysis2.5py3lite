LICENSE = """
======================COPYRIGHT/LICENSE START==========================

Util.py: Part of the CcpNmr Analysis program

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
import os, stat
import sys

from time import time

from memops.format.xml.Util import findTopObjectPath

from memops.universal.Io import getPythonDirectory, joinPath
from memops.universal.Region1D import Region1D
from memops.universal.Util import isWindowsOS, isMacOS

from memops.general.Implementation import ApiError
from memops.general import Version

try:
  from memops.gui.MessageReporter import showError
except ImportError:
  from memops.universal.MessageReporter import showError

from memops.api import Implementation

import memops.gui.Color as Color

from ccp.api.nmr.Nmr import FreqDataDim, SampledDataDim

from ccp.general.Command import Command

from ccp.util import Software

from ccpnmr.analysis import Copyright
from ccpnmr.analysis.core.ExperimentBasic import getPrimaryDataDimRef
from ccpnmr.analysis.core.ExperimentBasic import getSpectra, getNoiseEstimate
from ccpnmr.analysis.core.UnitConverter import unit_converter

toleranceDefaults = {'1H':0.05,'13C':0.2,'15N':0.4}
spectrum_contour_separator = ','

###spectrum_shortcuts = map(str, range(1, 10)) + ['0']
if isMacOS():
  nmax = 8
else:
  nmax = 12
spectrum_shortcuts = [ 'F%d' % x for x in range(1, nmax+1) ]

shift_key_modifier = 'Shift-'
ctrl_key_modifier = 'Ctrl-'
alt_key_modifier = 'Alt-'

graphics_handlers = []
try:
  import memops.c.TkHandler
  graphics_handlers.append('Tk')
except:
  pass
try:
  import memops.c.GlHandler
  graphics_handlers.append('OpenGL')
except:
  pass


def getHueSortedColorSchemes(analysisProfile):
  """Get the colour shemes of an AnalysisProfile in average hue order"""

  schemes = analysisProfile.colorSchemes
  
  sortList = [(getAverageHue(cs), cs) for cs in schemes]
  
  sortList.sort()
  
  return [x[1] for x in sortList]
 

def getAverageHue(colorScheme):
  """Get the average colour hue for a color scheme"""
  
  # Default black
  R = 0
  G = 0
  B = 0
  
  colors = colorScheme.colors
  
  for hexColor in colors:
    r, g, b = Color.hexToRgb(hexColor)
 
    R += r
    G += g
    B += b
  
  if colors:
    n = float(len(colors))
    R /= n
    G /= n
    B /= n   
  
  h, s, v = Color.rgbToHsb(R, G, B)
  
  return h

def getAnalysisProfile(project):
  """Retrieve or create an AnalysisProfile for a memops root"""

  analysisProfile = project.currentAnalysisProfile
  
  if analysisProfile:
    # Back compatibility
    # Cleanup old installation macros
    
    #for macro in analysisProfile.macros:
    #  print macro.path
    pass
  
  else:
    analysisProfile = project.newAnalysisProfile(name='standard')

  return analysisProfile

def getAnalysisProject(project):
  """Retrieve or create an AnalysisProject for a memops root"""

  analysisProject = project.currentAnalysisProject

  if not analysisProject:
    
    nmrProject = project.currentNmrProject
    
    if not nmrProject:
      nmrProject = project.newNmrProject(name=project.name)
    
    analysisProject = project.newAnalysisProject(name=project.name,
                                                 nmrProject=nmrProject)

  return analysisProject

def getAnalysisSpectrum(spectrum):
  """Retrieve or create an AnalysisSpectrum for an Nmr.DataSource"""

  analysisSpectrum = spectrum.analysisSpectrum
  
  if (not analysisSpectrum 
      or not (analysisSpectrum.posLevels or analysisSpectrum.negLevels)):
    analysisProject = getAnalysisProject(spectrum.root)
    
    if not analysisSpectrum:
      analysisSpectrum = analysisProject.newAnalysisSpectrum(dataSource=spectrum)
    
    setupAnalysisSpectrum(analysisSpectrum)
    
  return analysisSpectrum

def setupAnalysisSpectrum(analysisSpectrum):

  #analysisSpectrum = getAnalysisSpectrum(spectrum)
  spectrum = analysisSpectrum.dataSource

  if not (analysisSpectrum.posLevels or 
          analysisSpectrum.negLevels or 
          analysisSpectrum.posColors or
          analysisSpectrum.negColors):

    analysisProfile = getAnalysisProfile(spectrum.root)
  
    levels = getSpectrumContourLevels(spectrum)
    posLevels = [l for l in levels if l >= 0]
    negLevels = [l for l in levels if l < 0]
    
    analysisSpectrum.posLevels = posLevels
    analysisSpectrum.negLevels = negLevels

    if levels:
      defaultBase = min([abs(l) for l in levels]) or 1000.0
    else:
      defaultBase = 1000.0

    numberLevels = max(len(posLevels), len(negLevels), 1)

    posLevels.sort()
    if posLevels and len(posLevels) > 1:
      levelChanger = posLevels[1] / posLevels[0]
    else:
      levelChanger = 1.5

    application = spectrum.root.application
    baseLevel = application.getValue(spectrum,  keyword='baseLevel', defaultValue=defaultBase, deleteAppData=True)
    numberLevels = application.getValue(spectrum,  keyword='numberLevels', defaultValue=numberLevels, deleteAppData=True)
    levelChanger = application.getValue(spectrum,  keyword='levelChanger', defaultValue=levelChanger, deleteAppData=True)
    changeMode = application.getValue(spectrum,  keyword='levelChangeMode', defaultValue=0, deleteAppData=True) and 'add' or 'multiply'
   
    analysisSpectrum.autoBaseLevel    = baseLevel   
    analysisSpectrum.autoNumLevels    = numberLevels
    analysisSpectrum.autoLevelChanger = levelChanger
    analysisSpectrum.autoLevelMode    = changeMode  

    # reset the above to something sensible based on levels
    updateSpectrumLevelParams(analysisSpectrum, posLevels, negLevels)

    schemeName = getSpectrumPosContourColor(spectrum)
    scheme = analysisProfile.findFirstColorScheme(name=schemeName)
    
    if scheme:
      analysisSpectrum.posColors = scheme.colors
    else:
      analysisSpectrum.posColors = ['#000080',]

    schemeName = getSpectrumNegContourColor(spectrum)
    scheme = analysisProfile.findFirstColorScheme(name=schemeName)

    if scheme:
      analysisSpectrum.negColors = scheme.colors
    else:
      analysisSpectrum.negColors = ['#808000',]

    priority = getSpectrumContourPriority(spectrum)
    if priority == 'precalculated':
      bool = True
    else:
      bool = False
    analysisSpectrum.usePrecalculated = True
    
    bool = getSpectrumPointerVisibility(spectrum)
    analysisSpectrum.usePeakArrow = bool
    
    bool = getSpectrumBoxVisibility(spectrum)
    analysisSpectrum.useBoundingBox = bool
    
    rank = application.getValue(spectrum,  keyword='order', defaultValue=0, deleteAppData=True)
    analysisSpectrum.rank = rank
 
    shortcut = getSpectrumShortcut(spectrum)
    analysisSpectrum.shortcut = shortcut
    
    font = getSpectrumFont(spectrum)
    analysisSpectrum.font = '%s %s' % font
    
    if not spectrum.activePeakList:
      getSpectrumActivePeakList(spectrum)

  return analysisSpectrum

def getAnalysisDataDim(dataDim):
  """Retrieve or create an AnalysisDataDim for an Nmr.DataDim"""

  analysisDataDim = dataDim.analysisDataDim

  if not analysisDataDim:
    analysisSpectrum = getAnalysisSpectrum(dataDim.dataSource)
    analysisDataDim = analysisSpectrum.newAnalysisDataDim(dataDim=dataDim)

    if isinstance(dataDim, FreqDataDim):
 
      tol = getDataDimAssignmentTolerance(dataDim)
      analysisDataDim.assignTolerance = tol

      tol = getDataDimNoeTolerance(dataDim)
      analysisDataDim.noeTolerance = tol

      wt = getDataDimShiftWeight(dataDim)
      analysisDataDim.chemShiftWeight = wt

      application = dataDim.root.application
      
      minLw = application.getValue(dataDim, keyword='peakFindmin_lw',
                                   defaultValue = 0.0)
      analysisDataDim.peakFindMinLineWidth = minLw

      boxwidth = application.getValue(dataDim, keyword='peakFindboxwidth',
                                      defaultValue = 1.0)
      analysisDataDim.peakFindBoxWidth = boxwidth

    else:
 
      plane = getSampledDataDimReferencePlane(dataDim)
      analysisDataDim.refSamplePlane = plane
      

  return analysisDataDim


def getAnalysisPeakList(peakList):

  analysisPeakList = peakList.analysisPeakList

  if not analysisPeakList:
    analysisSpectrum = getAnalysisSpectrum(peakList.dataSource)
    analysisProfile = analysisSpectrum.root.currentAnalysisProfile
    if analysisProfile and analysisProfile.bgColor in (Color.black.hex, Color.darkgrey.hex):
      color = Color.white.hex
    else:
      color = Color.black.hex
    symbol = getPeakListSymbol(peakList)
    noeIntensityType = getSpectrumNoeIntensityType(peakList.dataSource)
    analysisPeakList = analysisSpectrum.newAnalysisPeakList(peakList=peakList, symbolColor=color, textColor=color, symbolStyle=symbol, noeIntensityType=noeIntensityType)
    
  return analysisPeakList


# TBD: look at again if default value changes in data model
def getFormatConverterThreading(analysisProject):

  if isWindowsOS():
    return False
  
  # Added Rasmus 5/4/10 to allow starting FC without a project
  # This effectively sets the default threading to False 
  # (which should be the less errror-prone)
  if analysisProject is None:
    return False
  
  keyword = 'formatConverterThreading'
  application = analysisProject.root.application

  result = application.getValue(analysisProject, keyword=keyword)
  if result is not None:
    application.setValue(analysisProject, keyword=keyword,
                         value=None)
    analysisProject.isThreadingAllowed = result
 
  return analysisProject.isThreadingAllowed

def setFormatConverterThreading(analysisProject, threading):

  if isWindowsOS():
    threading = False

  analysisProject.isThreadingAllowed = threading



def getBackgroundContrast(project):

  analysisProfile = project.currentAnalysisProfile
  color = analysisProfile.bgColor
  
  return Color.inverseGrey(color)
  
def getIsotopeExclusion(isotope):

  defaultValue = 0.0

  project = isotope.root
  analysisProject = project.currentAnalysisProject
  axisType = analysisProject.findFirstAxisType(isotopes=(isotope,))
  if not axisType:
    return defaultValue

  return axisType.diagonalExclusion

def setIsotopeExclusion(isotope, exclusion):

  project = isotope.root
  analysisProject = project.currentAnalysisProject
  axisType = analysisProject.findFirstAxisType(isotopes=(isotope,))
  if not axisType:
    return

  axisType.diagonalExclusion = exclusion

def getShortcutSpectrum(project, shortcut):

  application = project.application
  #print 'getShortcutSpectrum1', shortcut, type(shortcut)
  spectra = getSpectra(project)
  for spectrum in spectra:
    #print 'getShortcutSpectrum2', spectrum.name
    #value = application.getValue(spectrum, keyword='shortcut')
    value = spectrum.analysisSpectrum.shortcut
    #print 'getShortcutSpectrum3', value, type(value)
    if value == shortcut:
      return spectrum
    # make alt and ctrl synonymous for now
    if value and shortcut:  # current usage always has shortcut not None but play safe
      if value.startswith(alt_key_modifier):
        if shortcut.startswith(ctrl_key_modifier):
          if value[len(alt_key_modifier):] == shortcut[len(ctrl_key_modifier):]:
            return spectrum
      elif value.startswith(ctrl_key_modifier):
        if shortcut.startswith(alt_key_modifier):
          if value[len(ctrl_key_modifier):] == shortcut[len(alt_key_modifier):]:
            return spectrum

  return None

spectrumContourPriorities = [ 'precalculated', 'on-the-fly', 'both' ]

def getSpectrumContourPriority(spectrum):

  application = spectrum.root.application
  value = application.getValue(spectrum,  keyword='contourPriority',
                               defaultValue=spectrumContourPriorities[0],
                               deleteAppData=True)

  # add this in because names changed, so protect against that
  if value not in spectrumContourPriorities:
    value = spectrumContourPriorities[0]

  return value


def getSpectrumFont(spectrum):

  application = spectrum.root.application
  name = application.getValue(spectrum,  keyword='fontName',
                              defaultValue='Helvetica', deleteAppData=True)
  size = application.getValue(spectrum,  keyword='fontSize',
                              defaultValue=10, deleteAppData=True)

  return (name, size)

def getSpectrumShortcut(spectrum):

  project = spectrum.root
  application = project.application
  shortcut = application.getValue(spectrum,  keyword='shortcut',
                                  deleteAppData=True)

  try:
    shortcut = int(shortcut)
    if shortcut == 0:
      shortcut = 10
    shortcut = 'F%d' % shortcut
    
  except:
    pass

  return shortcut


def getNamedSpectrumView(window, name):

  expName, specName = stringToExperimentSpectrum(name)

  for pane in window.spectrumWindowPanes:
    for view in pane.spectrumWindowViews:
      analysisSpectrum = view.analysisSpectrum
      dataSource = analysisSpectrum.dataSource
      
      if dataSource.name == specName \
       and dataSource.experiment.name == expName:
 
        return view

  return None

def getSpectrumActivePeakList(spectrum):

  peakList = spectrum.activePeakList

  if not peakList:
    # v1 methodology
    application = spectrum.root.application
    serial = application.getValue(spectrum,  keyword='activePeakList',
                                  deleteAppData=True)
    peakList = spectrum.findFirstPeakList(serial=serial)

    if not peakList:
      peakList = spectrum.findFirstPeakList()
      for peakList0 in spectrum.sortedPeakLists():
        if peakList0.findFirstPeak():
          peakList = peakList0
          break

    if peakList:
      spectrum.activePeakList = peakList

  return peakList

def getSpectrumNoeIntensityType(spectrum):

  project     = spectrum.root
  application = project.application
  value = application.getValue(spectrum,  keyword='noeIntensityType',
                               defaultValue='volume', deleteAppData=True)

  return value

def getSampledDataDimReferencePlane(dataDim):

  project     = dataDim.root
  application = project.application
  value = application.getValue(dataDim,  keyword='sampleReferencePlane',
                               defaultValue=None)

  return value


def setSpectrumNoeIntensityType(spectrum, intensityType):

  application = spectrum.root.application
  application.setValue(spectrum,  keyword='noeIntensityType',
                        value=intensityType)

NOE_DIST = 'noeDistanceClass'

def getSpectrumNoeDistanceClasses(spectrum):

  distanceClasses = []
  project = spectrum.root

  if hasattr(project, 'application'):
    application = project.application

    # Old method for back compatibility
    for i in range(50): # Suitably large
      value = application.getValue(spectrum,  keyword='%s%d' % (NOE_DIST,1),
                                   defaultValue=None, deleteAppData=True)
      
      if value:
        distanceClasses.append( [float(x) for x in value.split(':')] )
      else:
        break                              
    
    # New method
    for datum in spectrum.findAllApplicationData(application=application.name,
                                                 keyword=NOE_DIST):
      distanceClasses.append( [float(x) for x in datum.value.split(':')] )


  return distanceClasses

def setSpectrumNoeDistanceClasses(spectrum, distanceClasses):

  project = spectrum.root
  application = project.application

  appData = spectrum.findAllApplicationData(application=application.name,
                                            keyword=NOE_DIST)

  for appDatum in appData:
    spectrum.removeApplicationData(appDatum)
  
  for datum in distanceClasses:
    value = ':'.join([str(x) for x in datum])
    appDatum = Implementation.AppDataString(application=application.name,
                                            keyword=NOE_DIST, value=value)
    spectrum.addApplicationData(appDatum)                                      

def setSpectrumNoeResidueRanges(spectrum, ranges):

  keyword = 'noeResidueRange'
  project = spectrum.root
  application = project.application
  appData = spectrum.findAllApplicationData(application=application.name,
                                            keyword=keyword)

  for appDatum in appData:
    spectrum.removeApplicationData(appDatum)

  for dims, chain, startRes, endRes in ranges:
    dimNums = []
    for dataDim in dims:
      dataDims = spectrum.sortedDataDims()
      if dataDim in dataDims:
        dimNums.append(str(dataDim.dim-1))

    data = [chain.molSystem.code, chain.code,
            'seqId'+str(startRes.seqId),
            'seqId'+str(endRes.seqId)]
    data.extend(dimNums)
    value = '%'.join(data)
 
    appDatum = Implementation.AppDataString(application=application.name,
                                            keyword=keyword, value=value)
    spectrum.addApplicationData(appDatum)                                      
                                            

def getSpectrumNoesResidueRanges(spectrum):

  keyword = 'noeResidueRange'
  ranges = []
  project = spectrum.root
  application = project.application
  aName = application.name
  
  # Old method for back compatibility
  values = []
  for i in range(50): # Arbitrarily large
    value = application.getValue(spectrum,  keyword='%s%d' % (keyword,i),
                                 defaultValue=None, deleteAppData=True)
    if value:
      values.append(value)
    else:
      break                              
  
  # New method
  appData = list(spectrum.findAllApplicationData(application=aName, keyword=keyword))
  values.extend([appDatum.value for appDatum in appData])
  
  sortedDataDims = spectrum.sortedDataDims()                       
  for datum in values:
    data = datum.split('%')
    (msCode, chainCode, startVal, endVal) = data[:4]

    dataDims  = [sortedDataDims[int(x)] for x in data[4:]]
    molSystem = None
    chain     = None

    useSeqIds = False
    if startVal[:5] == 'seqId':
      useSeqIds = True
      start     = int(startVal[5:])
      end       = int(endVal[5:])

    else:
      start     = int(startVal)
      end       = int(endVal)

    molSystem = project.findFirstMolSystem(code=msCode)
    if molSystem:
      chain = molSystem.findFirstChain(code=chainCode)

    if chain:
      if useSeqIds:
        startRes = chain.findFirstResidue(seqId=start)
        endRes   = chain.findFirstResidue(seqId=end)
      else:
        startRes = chain.findFirstResidue(seqCode=start)
        endRes   = chain.findFirstResidue(seqCode=end)

      if not startRes:
        startRes = chain.findFirstResidue()

      if not endRes:
        endRes = chain.sortedResidues()[-1]

      ranges.append([dataDims, chain, startRes, endRes])

  #print "getSpectrumNoesResidueRanges", ranges

  return ranges

def setSpectrumNoeShiftRanges(spectrum, ranges):

  keyword = 'noeShiftRange'
  project = spectrum.root
  application = project.application
  appData = spectrum.findAllApplicationData(application=application.name,
                                            keyword=keyword)
  for appDatum in appData:
    spectrum.removeApplicationData(appDatum)

  for dataDim, isotopes, start, end in ranges:
    dimNum = dataDim.dim - 1
    isotopeCodes =  ','.join(isotopes)
    data = [str(dimNum), isotopeCodes, str(start), str(end)]
    value = '%'.join(data)
    
    appDatum = Implementation.AppDataString(application=application.name,
                                            keyword=keyword, value=value)
    spectrum.addApplicationData(appDatum)                                      


def getSpectrumNoesShiftRanges(spectrum):

  keyword = 'noeShiftRange'
  ranges = []
  project = spectrum.root
  application = project.application
  aName = application.name
  
  # Old method for back compatibility
  values = []
  for i in range(50): # Arbitrarily large
    value = application.getValue(spectrum,  keyword='%s%d' % (keyword,i),
                               defaultValue=None, deleteAppData=True)
    if value:
      values.append(value)
    else:
      break                              

  # New method
  appData = list(spectrum.findAllApplicationData(application=aName, keyword=keyword))
  values.extend([appDatum.value for appDatum in appData])
  
  for datum in values:
    data = datum.split('%')
    (dimNum, isotopeCodes, startVal, endVal) = data[:4]

    isotopes = isotopeCodes.split(',')
    dataDim = spectrum.findFirstDataDim(dim=int(dimNum)+1)
    start   = float(startVal)
    end     = float(endVal)

    ranges.append([dataDim, isotopes, start, end])

  return ranges

def getSpectrumPosContourColor(spectrum):

  colors = ('teal','purple','green','navy','red','skyblue','maroon','pink')

  project = spectrum.root
  application = project.application
  initialColor = colors[ (int(spectrum.experiment.serial) % 8)-1 ]
  return application.getValue(spectrum,  keyword='posContourColor',
                              defaultValue=initialColor, deleteAppData=True)

def getSpectrumNegContourColor(spectrum):

  colors = ('magenta','orange','mauve','yellow',
            'teal','orange','cyan','midgreen')

  project = spectrum.root
  application = project.application
  initialColor = colors[ (int(spectrum.experiment.serial) % 8)-1 ]
  return application.getValue(spectrum,  keyword='negContourColor',
                              defaultValue=initialColor, deleteAppData=True)

def getSpectrumSliceColor(spectrum):

  name = getSpectrumPosContourColor(spectrum)
  getScheme = spectrum.root.currentAnalysisProject.findFirstColorScheme
  colorScheme = getScheme(name=name)
  if (colorScheme):
    defaultValue = colorScheme.colorNames[0]
  else:
    defaultValue = getBackgroundContrast(spectrum.root)

  application = spectrum.root.application
  return application.getValue(spectrum,  keyword='sliceColor',
                              defaultValue=defaultValue, deleteAppData=True)

def getSpectrumBoxVisibility(spectrum):

  application = spectrum.root.application
  return application.getValue(spectrum, keyword='isSpectrumBoxShown',
                              defaultValue=True, deleteAppData=True)

def getSpectrumPointerVisibility(spectrum):

  application = spectrum.root.application
  return application.getValue(spectrum, keyword='isSpectrumPointerShown',
                              defaultValue=True, deleteAppData=True)

def updateSpectrumLevelParams(analysisSpectrum, posLevels, negLevels):
  
  baseLevel    = analysisSpectrum.autoBaseLevel   
  numberLevels = analysisSpectrum.autoNumLevels   
  levelChanger = analysisSpectrum.autoLevelChanger
  changeMode   = analysisSpectrum.autoLevelMode   

  if posLevels:
    posLevels = list(posLevels)
    posLevels.sort()
    baseLevel = posLevels[0]
    numberLevels = len(posLevels)
    if numberLevels > 1:
      if changeMode == 'multiply':
        levelChanger = posLevels[1] / baseLevel
      else:
        levelChanger = posLevels[1] - baseLevel
        
  elif negLevels:
    negLevels = list(negLevels)
    negLevels.sort()
    baseLevel = negLevels[-1]
    numberLevels = len(negLevels)
    if numberLevels > 1:
      if changeMode == 'multiply':
        levelChanger = negLevels[-2] / baseLevel
      else:
        levelChanger = negLevels[-2] - baseLevel

  analysisSpectrum.autoBaseLevel    = abs(baseLevel)
  analysisSpectrum.autoNumLevels    = numberLevels
  analysisSpectrum.autoLevelChanger = levelChanger

def getSpectrumContourLevels(spectrum):

  analysisSpectrum = spectrum.analysisSpectrum
  levels = list(analysisSpectrum.posLevels + analysisSpectrum.negLevels)

  if not levels:
    # below is v1 way of doing contour levels
    # problem with using array of floats instead of string containing array
    # is that former means notifier gets updated for each level being added
    # (and removed) so get multiple updates which is disconcerting
    if not hasattr(spectrum.root, 'application'):
      return []

    application = spectrum.root.application
    #levels = application.getValues(spectrum,  keyword='contourLevel')
    levels = application.getValue(spectrum,  keyword='contourLevel',
                                  deleteAppData=True)
    if levels:
      levels = [float(level.strip()) for level in levels.split(spectrum_contour_separator) if level.strip()]
    else:
      levels = []

  return levels

def defaultWindowGroup(project):

  analysisProject = project.currentAnalysisProject
  group = analysisProject.activeWindowGroup
  if (not group):
    group = analysisProject.findFirstSpectrumWindowGroup()
    if group is None:
      group = analysisProject.newSpectrumWindowGroup(name='default_group')
    analysisProject.activeWindowGroup = group

def getPeakListSymbol(peakList):

  application = peakList.root.application
  symbol = application.getValue(peakList,  keyword='symbol',
                              defaultValue='x')

  if symbol == 'cross':
    symbol = 'x'

  if symbol == 'plus':
    symbol = '+'

  return symbol

def getPeakListColor(peakList):

  analysisPeakList = peakList.analysisPeakList
  if analysisPeakList:
    return analysisPeakList.symbolColor
  else:
    return Color.black.hex

def getDataDimAssignmentTolerance(dataDim):

  isotope = getPrimaryDataDimRef(dataDim).expDimRef.isotopeCodes[0]
  default = toleranceDefaults.get(isotope, 0.1)

  """
  if tolerance is None:
    axisType = Util.findAxisTypeMatch(dataDim)
    if (axisType is None): # should not be True
      tol = 0.1 # arbitrary
    else:
      tol = axisType.peakSize
  """

  application = dataDim.root.application
  tolerance = application.getValue(dataDim,  keyword='assignmentTolerance',
                                   defaultValue=default)

  return tolerance

def getDataDimNoeTolerance(dataDim):

  if hasattr(dataDim,'noeTolerance' ):
    return dataDim.noeTolerance

  application = dataDim.root.application
  default     = getDataDimAssignmentTolerance(dataDim)
  tolerance   = application.getValue(dataDim,  keyword='noeTolerance',
                                     defaultValue=default)
  dataDim.noeTolerance = tolerance

  return tolerance


def getDataDimShiftWeight(dataDim):

  if hasattr(dataDim,'shiftWeight' ):
    return dataDim.shiftWeight

  application = dataDim.root.application
  weight = application.getValue(dataDim,  keyword='shiftWeight',
                                defaultValue=1.0)
  dataDim.shiftWeight = weight

  return weight


def getPeakDimTextOffset(peakDim):

  application = peakDim.root.application
  textOffset = application.getValue(peakDim, keyword='textOffset')

  return textOffset

def setPeakDimTextOffset(peakDim, textOffset):

  application = peakDim.root.application
  application.setValue(peakDim, keyword='textOffset', value=textOffset)

  return textOffset


def getPeakTextOffset(peak):

  application = peak.root.application
  textOffset = application.getValue(peak, keyword='textOffset')

  return textOffset

def setPeakTextOffset(peak, textOffset):

  application = peak.root.application
  application.setValue(peak, keyword='textOffset', value=textOffset)

  return textOffset


def isTopObjectAnalysisInitialised(topObject, tolerance=50):
  """tolerance is max difference between app data save time and file mod time (in seconds)
     RHF 5 June 2014. Increased tolerance from 10 to 50s, as discrepancy can be at least 28s
  """

  repository = topObject.findFirstActiveRepository()
  if not repository:
    return False

  file = findTopObjectPath(repository.url.path, topObject)

  if not os.path.exists(file):  # it should exist
    return False

  modTime = os.stat(file)[stat.ST_MTIME]

  saveTime = getTopObjectAnalysisSaveTime(topObject)

  if saveTime:
    if abs(modTime-saveTime) < tolerance:
      return True

    else:
      return False

  else:
    return False


def getTopObjectAnalysisSaveTime(topObject):

  application = topObject.root.application
  timeStamp  = application.getValue(topObject, keyword='CcpNmrAnalysisTimeStamp')

  return timeStamp


def setTopObjectAnalysisSaveTime(topObject):

  timeStamp = int(time())

  application = topObject.root.application
  application.setValue(topObject, keyword='CcpNmrAnalysisTimeStamp',
                        value=timeStamp)

  #print "Set app time stamp", timeStamp, topObject

def findPanelType(project, axisType, axisPanels):

  analysisProject = project.currentAnalysisProject

  panelTypes = [ panelType for panelType in analysisProject.panelTypes \
                   if panelType.axisType == axisType ]
  n = len(panelTypes) + 1

  for axisPanel in axisPanels:
    if (axisPanel.panelType in panelTypes):
      panelTypes.remove(axisPanel.panelType)

  if panelTypes:
    panelType = panelTypes[0]
    
  else:
    symbols = [ isotope.chemElement.symbol for isotope in axisType.isotopes ]
    names = [ panelType.name for panelType in analysisProject.panelTypes ]
    name = ''.join(symbols) + str(n)
    while (name in names):
      n = n + 1
      name = ''.join(symbols) + str(n)
    panelType = analysisProject.newPanelType(name=name, axisType=axisType)

  return panelType

def createAxisPanel(windowPane, axisType, region = None):

  # TBD: below a bit dangerous
  axisPanels = windowPane.sortedAxisPanels()[2:]
  n = 0
  for axisPanel in axisPanels:
    m = int(axisPanel.label[1:])
    n = max(n, m)
    
  label = 'z' + str(n+1)
  panelType = findPanelType(windowPane.root, axisType, windowPane.sortedAxisPanels())
  if axisType.axisUnits:
    axisUnit = axisType.findFirstAxisUnit(unit='ppm')
  else:
    axisUnit = None
    
  panel = windowPane.newAxisPanel(label=label, panelType=panelType,
                                  axisUnit=axisUnit)
                            
  if region is None:
    (r0, r1) = axisType.region
    r = Region1D(r0, r1)
    r.zoom(0.01) # just make it some smallish fraction of the world
    region = (r[0], r[1])
    
  axisRegion = panel.findFirstAxisRegion()
  axisRegion.region = region
  windowPane.newSlicePanel(label=label)
  panel.isVisible = True
  panel.firstMapping = True

  return panel

# axisPanelMapping maps dataDim to axisPanel of window
def createSpectrumWindowView(windowPane, spectrum, axisPanelMapping,
                             isPosVisible=True, isNegVisible=True,
                             isInToolbar=True, isSliceVisible=True):

  analysisSpectrum = getAnalysisSpectrum(spectrum)

  try:
    view = windowPane.newSpectrumWindowView(analysisSpectrum=analysisSpectrum,
                                            isPosVisible=isPosVisible,
                                            isNegVisible=isNegVisible,
                                            isInToolbar=isInToolbar,
                                            isSliceVisible=isSliceVisible)
  except:
    return None # trying to create a second view with same key

  for dataDim in spectrum.dataDims:
    analysisDataDim = getAnalysisDataDim(dataDim)
  
    axisPanel = axisPanelMapping[dataDim]
    view.newAxisMapping(label=axisPanel.label, analysisDataDim=analysisDataDim)

  return view

def addAxisPanelRegion(axisPanel, region, size):

  #print 'addAxisPanelRegion', region, size
  axisPanel.newAxisRegion(region=region, size=size)

expt_spectrum_separator = ': '

def stringToExperimentSpectrum(s):

  try:
    (experiment_name, spectrum_name) = s.split(expt_spectrum_separator)
  except:
    experiment_name = spectrum_name = ''

  return (experiment_name, spectrum_name)

def stringFromExperimentSpectrum(experiment_name, spectrum_name):

  return experiment_name + expt_spectrum_separator + spectrum_name

def defaultProject():

  project = Implementation.MemopsRoot(name='defaultProject')

  return project

def defaultColors(project):

  color_dict = {}
  for c in Color.standardColors:
    color_dict[c.name] = c.hex

  # default color schemes
  colorSchemeTable = ( \
    ('orangeshade' , ('#FFE0C0','#FFC890','#FFB060','#FF9830','#FF8000',
                      '#E17100','#C26100','#A35200','#854200','#663300')),
    ('rainbow'     , ('#FF00FF','#FF0080','#FF0000','#FF8000','#FFFF00','#80FF00',
                      '#00FF00','#00FF80','#00FFFF','#0080FF','#0000FF','#8000FF')),
    ('cyanshade'   , ('#00FFFF','#00ECEC','#00D8D8','#00C4C4','#00B0B0',
                      '#009C9C','#008888','#007474','#006060','#004C4C')),
    ('greyshade'   , ('#CCCCCC','#BBBBBB','#AAAAAA','#999999','#888888',
                      '#777777','#666666','#555555','#444444','#333333')),
    ('wimbledon'   , ('#008000','#1C8E00','#389C00','#55AB00','#71B900',
                      '#8EC700','#AAD500','#C7E300','#E3F100','#FFFF00')),
    ('yellowshade' , ('#FFFF99','#FFFF4C','#FFFF00','#E7E700','#CFCF00',
                      '#B6B600','#9E9E00','#868600','#6D6D00','#555500')),
    ('redshade'    , ('#FFC0C0','#FF9A9A','#FF7373','#FF4D4D','#FF2626',
                      '#FF0000','#D90000','#B30000','#8C0000','#660000')),
    ('purpleshade' , ('#E6CCFF','#D399F0','#C066E0','#AC33D0','#9900C0',
                      '#8500AC','#700097','#5C0082','#47006E','#330059')),
    ('toothpaste'  , ('#C0C0FF','#9A9AFF','#7373FF','#4D4DFF','#2626FF',
                      '#0000FF','#0040FF','#0080FF','#00C0FF','#00FFFF')),
    ('cmy'         , ('#00FFFF','#33CCFF','#6699FF','#9966FF','#CC33FF',
                      '#FF00FF','#FF33CC','#FF6699','#FF9966','#FFCC33','#FFFF00')),
    ('steel'       , ('#C0C0C0','#ABABB9','#9595B2','#8080AB','#6B6BA4',
                      '#55559D','#404095','#2A2A8E','#151587','#000080')),
    ('rgb'         , ('#FF0000','#CC1900','#993300','#664D00','#336600',
                      '#008000','#006633','#004D66','#003399','#0019CC','#0000FF')),
    ('tropicana'   , ('#FFFF00','#FFE30E','#FFC71C','#FFAA2A','#FF8E39',
                      '#FF7147','#FF5555','#FF3863','#FF1C72','#FF0080')),
    ('sunset'      , ('#FFC0C0','#FF9A9A','#FF7373','#FF4D4D','#FF2626',
                      '#FF0000','#FF4000','#FF8000','#FFC000','#FFFF00')),
    ('magma'       , ('#000000','#400000','#800000','#C00000','#FF0000',
                      '#FF3300','#FF6600','#FF9900','#FFCC00','#FFFF00')),
    ('holly'       , ('#80FF80','#66E666','#4DCD4D','#33B333','#199A19',
                      '#008000','#FF0000','#D50000','#AB0000','#800000')),
    ('greenshade'  , ('#99FF99','#73F073','#4CE04C','#26D026','#00C000',
                      '#00AE00','#009C00','#008A00','#007800','#006600')),
    ('glacier'     , ('#000000','#000040','#000080','#0000C0','#0000FF',
                      '#2626FF','#4D4DFF','#7373FF','#9A9AFF','#C0C0FF')),
    ('monarchy'    , ('#C0C0FF','#6060FF','#0000FF','#3300CC','#660099',
                      '#990066','#CC0033','#FF0000','#C00000','#800000')),
    ('blueshade'   , ('#C0C0FF','#9A9AFF','#7373FF','#4D4DFF','#2626FF',
                      '#0000FF','#0000D9','#0000B3','#00008C','#000066')),
    ('contrast'    , ('#FF0000','#008000','#0000FF','#FFFF00','#FF00FF',
                      '#00FFFF')),
  )
  
  analysisProfile = project.currentAnalysisProfile
  if not analysisProfile.colorSchemes:

    for color in Color.standardColors:
      analysisProfile.newColorScheme(name=color.name, colors=[color.hex,])

    for scheme in colorSchemeTable:
      (name, colorValues) = scheme
      analysisProfile.newColorScheme(name=name, colors=colorValues)

# different than other defaults because no Symbol class


default_axis_units_table = ( \
    # unit, isBackwards
    ('ppm', True),
    ('Hz', True),
    ('K', False),
    ('C', False),
    ('sec', False),
    ('arbitrary', False),
  )

def defaultAxisUnits(project):

  analysisProject = project.currentAnalysisProject
  for entry in default_axis_units_table:
    (unit, isBackwards) = entry
    if (not analysisProject.findFirstAxisUnit(unit=unit)):
      # axis unit not defined yet
      analysisProject.newAxisUnit(unit=unit, isBackwards=isBackwards)

default_axis_types_table = ( \
    # name, region, isSampled, numDecimals, peakSize, units, isotopeCodes, measurementType
    ('1H', (-4.0, 20.0), False, 3, 0.02, ('ppm', 'Hz'), ('1H',), 'Shift'),
    ('13C', (0.0, 250.0), False, 2, 0.2, ('ppm', 'Hz'), ('13C',), 'Shift'),
    ('14N', (80.0, 155.0), False, 2, 0.2, ('ppm', 'Hz'), ('14N',), 'Shift'),
    ('15N', (80.0, 155.0), False, 2, 0.2, ('ppm', 'Hz'), ('15N',), 'Shift'),
    ('2H', (-4.0, 20.0), False, 3, 0.1, ('ppm', 'Hz'), ('2H',), 'Shift'),
    ('19F', (-250.0, 100.0), False, 2, 0.2, ('ppm', 'Hz'), ('19F',), 'Shift'),
    ('31P', (-50.0, 100.0), False, 2, 0.2, ('ppm', 'Hz'), ('31P',), 'Shift'),
    ('79Br', (-4000.0, 2000.0), False, 2, 0.2, ('ppm', 'Hz'), ('79Br',), 'Shift'),
    ('T', (273.0, 333.0), True, 1, 1.0, ('K', 'C'), (), 'temperature'),
    ('t', (-1.0e6, 1.0e6), True, 0, 1.0, ('sec',), (), 'time'),
    ('sampled', (0.0, 10.0), True, 2, 1.0, (), (), 'sampled'),
    ('value', (-1.0e6, 1.0e6), False, 0, 1.0, ('arbitrary',), (), 'none'),
    ('DQ13C', (-50.0, 400.0), False, 2, 0.2, ('ppm', 'Hz'), ('13C',), 'MQShift'),
  )

def defaultAxisTypes(project):

  defaultAxisUnits(project)

  analysisProject = project.currentAnalysisProject
  for entry in default_axis_types_table:
    (name, region, isSampled, numDecimals, peakSize, units, isotopeCodes, measurementType) = entry
    if (not analysisProject.findFirstAxisType(name=name)):
      # axis type not defined yet
      t = analysisProject.newAxisType(name=name, region=region, isSampled=isSampled,
                            numDecimals=numDecimals, peakSize=peakSize,
                            isotopeCodes=isotopeCodes, measurementType=measurementType)
      for unit in units:
        axis_unit = analysisProject.findFirstAxisUnit(unit=unit)
        t.addAxisUnit(axis_unit)

def isIsotopeAxisType(axisType):

  # a bit of a hack this
  return axisType.name in ('1H', '2H', '13C', '15N', '19F', '31P')

default_isotopes_table = ( \
      ('H', 1),
      ('C', 13),
      ('N', 15),
      ('H', 2),
      ('F', 19),
      ('P', 31),
  )

def getNmrIsotopes(project):

  isotopes = []
  for entry in default_isotopes_table:
    (symbol, mass) = entry
    nucleus = project.currentChemElementStore.findFirstChemElement(symbol=symbol)
    if nucleus:
      isotope = nucleus.findFirstIsotope(massNumber=mass)
      if isotope:
        isotopes.append(isotope)

  return isotopes

def fitViewInWorld(view, world):

  #print 'fitViewInWorld', view, world

  if (world[0] < world[1]):
    if (view[1] <= world[0] or view[0] >= world[1]):
      view[0] = world[0]
      view[1] = world[1]
    else:
      view[0] = max(view[0], world[0])
      view[1] = min(view[1], world[1])
  else:
    if (view[1] >= world[0] or view[0] <= world[1]):
      view[0] = world[0]
      view[1] = world[1]
    else:
      view[0] = min(view[0], world[0])
      view[1] = max(view[1], world[1])

def expandWorldToView(view, world):

  #print 'expandWorldToView', view, world

  if (world[0] < world[1]):
    world[0] = max(view[0], world[0])
    world[1] = min(view[1], world[1])
  else:
    world[0] = min(view[0], world[0])
    world[1] = max(view[1], world[1])

# convert position from one unit to another
def convertPosition(position, dataDimRef, fromUnit = 'point', toUnit = 'point',
                    relative = False):

  if (fromUnit == toUnit):
    return position

  converter = unit_converter[(fromUnit, toUnit)]
  position = converter(position, dataDimRef)

  if (relative):
    position = position - converter(0.0, dataDimRef)

  return position

# convert region to points (starting from 0)
def convertRegion(region, axisUnit, dataDim):

  if isinstance(dataDim, FreqDataDim):

    region = checkSwapRegion(region, axisUnit)

    unit = axisUnit.unit
    converter = unit_converter[(unit, 'point')]
    # TBD: more general dataDimRefs
    # -1 is because points start from 1 but want region to start from 0
    region = [ converter(x, getPrimaryDataDimRef(dataDim))-1 for x in region ]

  # TBD: this assumes region defined in planes counting from 1
  elif isinstance(dataDim, SampledDataDim):

    region = [x - 1 for x in region]

  return region

def haveTypeMatch(axisPanel, dataDim):

  axisType = axisPanel.panelType.axisType
  if isinstance(dataDim, FreqDataDim):
    #TBD: look at below again
    #print 'haveTypeMatch', dataDim.dim, dataDim.findFirstDataDimRef().expDimRef.isotopes, axisType.isotopes
    
    if dataDim.dataDimRefs:
      expDimRef = getPrimaryDataDimRef(dataDim).expDimRef
      
      if expDimRef.measurementType != axisType.measurementType:
        # E.g. don't match DQ to normal axes
        return False
    
      if expDimRef.isotopes == axisType.isotopes:
        return True
        
  elif isinstance(dataDim, SampledDataDim):
    return axisType.name in ('T', 't', 'sampled')

  return False

def findAxisTypeMatch(dataDim):

  project = dataDim.root
  analysisProject = project.currentAnalysisProject
  findAxisType = analysisProject.findFirstAxisType
  if (isinstance(dataDim, FreqDataDim)):
    #TBD: look at below again
    
    expDimRef = getPrimaryDataDimRef(dataDim).expDimRef
    axisType = findAxisType(isotopeCodes=expDimRef.isotopeCodes,
                            measurementType=expDimRef.measurementType)
                            
  elif (isinstance(dataDim, SampledDataDim) and \
        (dataDim.conditionVaried == 'temperature')):
    axisType = findAxisType(name='T')
    
  # TBD: below is short-term hack
  elif (isinstance(dataDim, SampledDataDim) and \
        (not dataDim.conditionVaried)):
    axisType = findAxisType(name='sampled')
    
  else:
    axisType = None

  return axisType


def changeSpectrumContourLevels(analysisSpectrum, levelChanger, changeMode='multiply'):
  """Descrn: Change the contour levels of a spectrum by a given factor
     Inputs: Analysis.AnalysisSpectrum, Float, Word (multiply or add)
     Output: None
  """
 
  pos = analysisSpectrum.posLevels
  neg = analysisSpectrum.negLevels
  
  if changeMode == 'multiply':
    pos = [ levelChanger*level for level in pos ]
    neg = [ levelChanger*level for level in neg ]
 
  else:
    pos = [ levelChanger+abs(level) for level in pos ]
    neg = [ -1 * (levelChanger+abs(level)) for level in neg ]
    
  if pos and min(pos) > 0:
    analysisSpectrum.posLevels = pos

  if neg and max(neg) < 0:
    analysisSpectrum.negLevels = neg

# set up contour levels for spectrum if none currently exist
def defaultContourLevels(spectrum, updateContourLevels = False):

  analysisSpectrum = getAnalysisSpectrum(spectrum)
  analysisProject = analysisSpectrum.analysisProject
  
  posLevels = analysisSpectrum.posLevels 
  negLevels = analysisSpectrum.negLevels  
  
  if not updateContourLevels and (posLevels or negLevels): # already have contour levels
    return

  if not hasattr(spectrum, 'block_file') or not spectrum.block_file:
    return

  v = 3 * getNoiseEstimate(spectrum) / analysisProject.globalContourScale
  nLevels = 10
  multiplier = 1.5
  w = [v]
  for n in range(nLevels-1):
    v = v * multiplier
    w.append(v)

  #print "defaultContourLevels", getNoiseEstimate(spectrum), analysisProject.globalContourScale, w

  analysisSpectrum.posLevels = w
  analysisSpectrum.negLevels = [-v,]

def greaterOrEqualEntry(entry, entries):

  for e in entries:
    if (entry <= e):
      break
  else:
    e = entries[-1]

  return e

# get software if it exists, otherwise create
def getSoftware(project):

  name = Copyright.suite + ' ' + Copyright.program
  version = Copyright.version
  version = '%s.%s' % (version.major, version.minor)

  return Software.getSoftware(project, name, version, Copyright.vendorName, Copyright.vendorAddress, Copyright.vendorWebAddress, Copyright.details)

# get method if it exists, otherwise create
# two methods the same if software, procedure and parameters are the same
def getMethod(project, task, procedure = None, parameters = None, details = None):

  software = getSoftware(project)

  return Software.getMethod(software, task, procedure, parameters, details)

# always unpack region into (r0, r1) in case region is not
# tuple but something fancier (with __getitem__ defined)
def checkSwapRegion(region, axisUnit):

  (r0, r1) = region
  if (axisUnit and axisUnit.isBackwards):
    (r0, r1) = (r1, r0)

  return [r0, r1]

def printStackTrace(fp = None):

  import sys
  import traceback
  if (not fp):
    fp = sys.stdout
  try:
    raise Exception()
  except:
    t = traceback.extract_stack()[:-1]
    m = len(t)
    for x in t:
     m = m - 1
     (file, line, func, code) = x
     n = file.rfind('.')
     if (n >= 0):
       file = file[:n]
     n = file.rfind('/')
     if (n >= 0):
       file = file[n+1:]
     fp.write('%3d: %6d %-20s %-20s %-s\n' % (m, line, file, func, code))


def runMacro(macro,argumentServer):

  path = macro.path
  if not os.path.isabs(path):  # assumed to be relative to our python directory
    path = joinPath(getPythonDirectory(), path)

  sys.path.append(path)
  try:
    command = Command(argumentServer,macro.name,macro.module,macro.function)
  finally:
    del sys.path[-1]
  command.run()

def reloadMacro(macro,argumentServer):

  path = macro.path
  if not os.path.isabs(path):  # assumed to be relative to our python directory
    path = joinPath(getPythonDirectory(), path)

  sys.path.append(path)
  try:
    command = Command(argumentServer,macro.name,macro.module,macro.function)
    command.reload()
  finally:
    del sys.path[-1]

# this is so that C world know what dims to skip for peak manipulations
# it does not take aliasing into consideration (use getDimWrapped() for that)
def getDimChecked(spectrum):

  if not hasattr(spectrum, 'dimChecked'):
    dimChecked = spectrum.numDim * [1]
    for dataDim in spectrum.dataDims:
      if not isinstance(dataDim, FreqDataDim):
        dimChecked[dataDim.dim-1] = 0
    spectrum.dimChecked = dimChecked

  return spectrum.dimChecked

# this is so that C world know what dims are wrapped
def getDimWrapped(spectrum):

  dimWrapped = []
  for dataDim in spectrum.sortedDataDims():
    if isinstance(dataDim, FreqDataDim) and dataDim.numPoints == dataDim.numPointsOrig:
      f = 1
    else:
      f = 0
    dimWrapped.append(f)

  return dimWrapped

def calcPointsPerPixel(view, axisLabel):

  axisPanel = view.spectrumWindowPane.findFirstAxisPanel(label=axisLabel)
  axisMapping = view.findFirstAxisMapping(label=axisLabel)
  if not axisPanel or not axisMapping:
    return 1.0

  dataDim = axisMapping.analysisDataDim.dataDim
  dataDimRef = getPrimaryDataDimRef(dataDim)
  if not dataDimRef:
    return 1.0

  axisRegion = axisPanel.findFirstAxisRegion()
  region = axisRegion.region
  d = region[1] - region[0]
  d = convertPosition(d, dataDimRef, fromUnit=axisPanel.axisUnit.unit, relative=True)

  return abs(d) / axisRegion.size


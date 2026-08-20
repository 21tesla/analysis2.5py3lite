
"""
======================COPYRIGHT/LICENSE START==========================

WindowDraw.py: Part of the CcpNmr Analysis program

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
import math

from memops.universal.Pdf import Pdf
from memops.universal.PostScript import PostScript
from memops.universal.PrintHandler import PrintHandler
from memops.universal.Region1D import Region1D
from memops.universal.Region2D import Region2D
from memops.universal.Util import cumulativeProductArray, arrayOfIndex, formatDecimals
 
from ccp.api.nmr import Nmr
from ccpnmr.analysis.core import ExperimentBasic
from ccpnmr.analysis.core.PeakBasic import getPeakDimPosition, getPeakAnnotation
from ccpnmr.analysis.core.PeakFindParams import getPeakFindParams, getPeakFindBoxwidth
from ccpnmr.analysis.core.PeakDrawUtil import peak_draw_methods, getCPeakDrawMethod
from ccpnmr.analysis.core.Util import getAnalysisSpectrum, getAnalysisProfile
from ccpnmr.analysis.core.Util import checkSwapRegion, calcPointsPerPixel, \
                                      fitViewInWorld, convertPosition, convertRegion, \
                                      getDimWrapped, isIsotopeAxisType

from ccpnmr.analysis.core.WindowBasic import windowPaneHasValueAxis, getSpectrumWindowView
from ccpnmr.analysis.core.WindowBasic import getAxisRegionRegion, getXAxisOffset, getValueAxisOffset
from memops.universal.Util import isWindowsOS

try:
  from memops.gui.Color import hexToRgb, inverseGrey, inverseRgb
except ImportError:
  pass # code below might fail, but not much can do, above functions should not be in gui

try:
  import ccpnmr.c.ContourStyle as ContourStyle
  import ccpnmr.c.ContourLevels as ContourLevels
except Exception, e:
  print 'Error, the WindowDraw module will not work, something is wrong with the C code.'
  ContourStyle = ContourLevels = None
  print 'Exception:', e
  print 'Will continue without Analysis window drawing functionality'

no_peak_text = 'No peak text'

class WindowDraw:

  extraPad = 0.5  # for peak drawing selection in orthogonal planes

  def __init__(self, parent, windowPane):
 
    self.parent = parent
    self.windowPane = windowPane

    self.hasValueAxis = windowPaneHasValueAxis(windowPane)
    self.deltaPositions = None
 
  # can be overridden in subclass
  def lift(self):

    pass

  # can be overridden in subclass
  def after_idle(self, func):

    func()

  def getNRows(self):

    return len(self.windowPane.findFirstAxisPanel(label='y').axisRegions)

  def getNCols(self):

    return len(self.windowPane.findFirstAxisPanel(label='x').axisRegions)

  def getAspectRatio(self):

    return self.windowPane.aspectRatio

  def gotoPeak(self, peak, row=None, col=None):

    windowPane = self.windowPane
    spectrum = peak.peakList.dataSource
    view = getSpectrumWindowView(windowPane, spectrum)

    if view:
      if windowPane.spectrumWindow.stripAxis == 'x':
        axisPanel   = self.windowPane.findFirstAxisPanel(label='x')
        axisRegions = axisPanel.sortedAxisRegions()
        cols = []
        if col is None:
          if len(axisRegions) > 1:
            i = 0
            for axisRegion in axisRegions:
              if axisRegion.isActive:
                cols.append(i)
              i += 1
 
          if not cols:
            cols = [0,]

        else:
          cols = [col,]

        if row is None:
          row = 0

        for axisMapping in view.axisMappings:
          axisPanel = self.windowPane.findFirstAxisPanel(label=axisMapping.label)
          axisRegions = axisPanel.sortedAxisRegions()
          peakDim = peak.findFirstPeakDim(dataDim=axisMapping.analysisDataDim.dataDim)
          if peakDim.dataDimRef:
            unit = axisPanel.axisUnit.unit
          else:
            unit = 'point'
          p = getPeakDimPosition(peakDim, toUnit=unit)
            
          for col in cols:
        
            if axisMapping.label == 'x':
              axisRegion = axisRegions[col]
            elif axisMapping.label == 'y':
              axisRegion = axisRegions[row]
            else:
              axisRegion = axisRegions[col]
            
            region = axisRegion.region
          
            (r0, r1) = region
            d = 0.5 * (r1 - r0)
            region = (p-d, p+d)
            axisRegion.region = region

      else: # windowPane.spectrumWindow.stripAxis == 'y':
        axisPanel   = self.windowPane.findFirstAxisPanel(label='y')
        axisRegions = axisPanel.sortedAxisRegions()
        rows = []
        if row is None:
          if len(axisRegions) > 1:
            i = 0
            for axisRegion in axisRegions:
              if axisRegion.isActive:
                rows.append(i)
              i += 1
 
          if not rows:
            rows = [0,]

        else:
          rows = [row,]

        if col is None:
          col = 0

        for axisMapping in view.axisMappings:
          axisPanel = self.windowPane.findFirstAxisPanel(label=axisMapping.label)
          axisRegions = axisPanel.sortedAxisRegions()
          peakDim = peak.findFirstPeakDim(dataDim=axisMapping.analysisDataDim.dataDim)
          if peakDim.dataDimRef:
            unit = axisPanel.axisUnit.unit
          else:
            unit = 'point'
          p = getPeakDimPosition(peakDim, toUnit=unit)
            
          for row in rows:
        
            if axisMapping.label == 'x':
              axisRegion = axisRegions[col]
            elif axisMapping.label == 'y':
              axisRegion = axisRegions[row]
            else:
              axisRegion = axisRegions[row]
            
            region = axisRegion.region
          
            (r0, r1) = region
            d = 0.5 * (r1 - r0)
            region = (p-d, p+d)
            axisRegion.region = region

      self.lift()

  # Note: position is assumed to be in same units as axisPanels in window
  #       it is a dictionary keyed on label
  #       axes with labels as keys have region changed so that position is at center
  def gotoPosition(self, position, row=None, col=None, doLift=True):

    windowPane = self.windowPane
    
    if windowPane.spectrumWindow.stripAxis == 'x':
      self.gotoXStripPosition(windowPane, position, row, col, doLift)
    else:
      self.gotoYStripPosition(windowPane, position, row, col, doLift)

  def gotoXStripPosition(self, windowPane, position, row=None, col=None, doLift=True):

    # setup column for strips; same for x, z1, z2...
    axisPanel = windowPane.findFirstAxisPanel(label='x')
    axisRegions = axisPanel.sortedAxisRegions()
    cols = []
    if col is None:
      if len(axisRegions) > 1:
        i = 0
        for axisRegion in axisRegions:
          if axisRegion.isActive:
            cols.append(i)
          i += 1
      
      if not cols:
        cols = [0,]

    else:
      cols = [col,]
                
    for axisRegion in axisRegions:
      axisRegion.isActive = False
      
    for col in cols:  
      axisRegions[col].isActive = True

    for col in cols:
      for label in position.keys():

        axisPanel = windowPane.findFirstAxisPanel(label=label)
        axisRegions = axisPanel.sortedAxisRegions()

        if label != 'y':
          axisRegion = axisRegions[col]

        else:
          if row is None:
            row = 0
            if len(axisRegions) > 1:
              i = 0
              for axisRegion in axisRegions[1:]:
                i += 1
                if axisRegion.isActive:
                  row = i
                  break
 
          axisRegion = axisRegions[row]
 
        region = axisRegion.region

        p = position[label]
        (r0, r1) = region
        d = 0.5 * (r1 - r0)
        region = (p-d, p+d)

        # print 'mid1 gotoXStripPosition', col, label, position
        if label in ('x', 'y'):
          axisRegion.region = region
          # above automatically updates scrollbar
        else:
          #(r0, r1) = checkSwapRegion(region, axisPanel.axisUnit)
          #axisPanel.region_selector.setViewRegion(r0, r1)
          # above directly updates scrollbar and indirectly sets axisRegion.region
          # should not be updating widgets in WindowDraw!
          axisRegion.region = region
        # print 'mid2 gotoXStripPosition', col, label, position

    if doLift:
      self.lift()

  def gotoYStripPosition(self, windowPane, position, row=None, col=None, doLift=True):

    # setup column for strips; same for x, z1, z2...
    findFirstAxisPanel = windowPane.findFirstAxisPanel
    axisPanel = findFirstAxisPanel(label='y')
    axisRegions = axisPanel.sortedAxisRegions()
    rows = []
    if row is None:
      if len(axisRegions) > 1:
        i = 0
        for axisRegion in axisRegions:
          if axisRegion.isActive:
            rows.append(i)
          i += 1
      
      if not rows:
        rows = [0,]

    else:
      rows = [row,]
                
    for axisRegion in axisRegions:
      axisRegion.isActive = False
      
    for row in rows:  
      axisRegions[row].isActive = True

    for row in rows:
      for label in position.keys():

        axisPanel = findFirstAxisPanel(label=label)
        axisRegions = axisPanel.sortedAxisRegions()

        if label != 'x':
          axisRegion = axisRegions[row]

        else:
          if col is None:
            col = 0
            if len(axisRegions) > 1:
              i = 0
              for axisRegion in axisRegions[1:]:
                i += 1
                if axisRegion.isActive:
                  col = i
                  break
 
          axisRegion = axisRegions[col]
 
        region = axisRegion.region

        p = position[label]
        (r0, r1) = region
        d = 0.5 * (r1 - r0)
        region = (p-d, p+d)

        #print 'mid1 gotoYStripPosition', row, label
        if (label in ('x', 'y')):
          axisRegion.region = region
          # above automatically updates scrollbar
        else:
          #(r0, r1) = checkSwapRegion(region, axisPanel.axisUnit)
          #axisPanel.region_selector.setViewRegion(r0, r1)
          # above directly updates scrollbar and indirectly sets axisRegion.region
          # should not be updating widgets in WindowDraw!
          axisRegion.region = region
        #print 'mid2 gotoYStripPosition', row, label

    if doLift:
      self.lift()

  def getSpectrumViews(self):

    return self.windowPane.sortedSpectrumWindowViews()

  def isViewPosVisible(self, view):

    if hasattr(view, 'printPositive'):
      return view.printPositive
    else:
      return view.isPosVisible

  def isViewNegVisible(self, view):

    if hasattr(view, 'printNegative'):
      return view.printNegative
    else:
      return view.isNegVisible

  def isViewVisible(self, view):

    if self.hasValueAxis:
      return view.isSliceVisible
    else:
      return self.isViewPosVisible(view) or self.isViewNegVisible(view)

  def isWinPeakListDrawn(self, winPeakList):

    if hasattr(winPeakList, 'printPeaks'):
      return winPeakList.printPeaks or winPeakList.printFont != no_peak_text
    else:
      return winPeakList.isSymbolDrawn or winPeakList.isAnnotationDrawn

  def getActiveSpectrumViews(self):

    return [ view for view in self.getSpectrumViews() if self.isViewVisible(view) ]

  def getActiveSliceViews(self):

    return [ view for view in self.getSpectrumViews() if view.isSliceVisible ]

  def getActivePeakLists(self):

    views = self.getActiveSpectrumViews()

    peakLists = []
    for view in views:
      peakList = view.analysisSpectrum.dataSource.activePeakList
      if peakList:
        peakLists.append(peakList)

    return peakLists

  def getMinMaxFreq(self, label, view):

    axisMapping = view.findFirstAxisMapping(label=label)
    axisPanel = axisMapping.axisPanel
    axisType = axisPanel.axisType
    dataDim = axisMapping.analysisDataDim.dataDim
    dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
    expDimRef = dataDimRef.expDimRef
    minAliasedFreq = expDimRef.minAliasedFreq
    maxAliasedFreq = expDimRef.maxAliasedFreq
    if minAliasedFreq is None:
      minAliasedFreq = convertPosition(float(dataDim.numPoints), dataDimRef, toUnit=axisType.findFirstAxisUnit(unit='ppm').unit)
 
    if maxAliasedFreq is None:
      maxAliasedFreq = convertPosition(1.0, dataDimRef, toUnit=axisType.findFirstAxisUnit(unit='ppm').unit)

    return (minAliasedFreq, maxAliasedFreq)

  def zoomToSpectraRegion(self, row, col):
    """ zoom x,y region just enough to see all visible spectra
    """
     
    minX = maxX = minY = maxY = None
    for view in self.windowPane.spectrumWindowViews:
      if self.isViewVisible(view):
        (minFreq, maxFreq) = self.getMinMaxFreq('x', view)
        if minX is None:
          minX = minFreq
          maxX = maxFreq
        else:
          minX = min(minX, minFreq)
          maxX = max(maxX, maxFreq)
        if not self.hasValueAxis:
          (minFreq, maxFreq) = self.getMinMaxFreq('y', view)
          if minY is None:
            minY = minFreq
            maxY = maxFreq
          else:
            minY = min(minY, minFreq)
            maxY = max(maxY, maxFreq)
          
    if minX is not None and maxX > minX and  minY is not None and maxY > minY:
      aspectRatio = self.windowPane.aspectRatio
      dx = maxX - minX
      dy = maxY - minY
      s = dx * aspectRatio
      if s < dy:
        dx = dy / aspectRatio
        minX -= 0.5 * dx
        maxX += 0.5 * dx
      else:
        dy = dx * aspectRatio
        minY -= 0.5 * dy
        maxY += 0.5 * dy

    if minX is not None and maxX > minX:
      axisPanel = self.windowPane.findFirstAxisPanel(label='x')
      axisRegion = axisPanel.sortedAxisRegions()[col]
      axisRegion.region = (minX, maxX)

    if minY is not None and maxY > minY:
      axisPanel = self.windowPane.findFirstAxisPanel(label='y')
      axisRegion = axisPanel.sortedAxisRegions()[row]
      axisRegion.region = (minY, maxY)

  # if returnDict is True, returns dictionary keyed on label
  # otherwise returns array in same order as self.windowPane.sortedAxisPanels()
  def findPosition(self, a, b, n, returnDict = True):

    axisPanels = self.windowPane.sortedAxisPanels()[2:]

    calcAxisMidpoint = self.calcAxisMidpoint
    if returnDict:

      position = {}
      position['x'] = a
      position['y'] = b

      for axisPanel in axisPanels:
        axisRegion = axisPanel.sortedAxisRegions()[n]
        position[axisPanel.label] = calcAxisMidpoint(axisRegion)

    else:

      position = [a, b]
      for axisPanel in axisPanels:
        axisRegion = axisPanel.sortedAxisRegions()[n]
        position.append(calcAxisMidpoint(axisRegion))

    return position

  # returns dictionary keyed on label
  # have position in (x, y) and region in other dimensions
  def findPositionRegion(self, a, b, n):

    position_region = {}
    position_region['x'] = a
    position_region['y'] = b

    for axisPanel in self.windowPane.sortedAxisPanels()[2:]:
      axisRegions = axisPanel.sortedAxisRegions()
      if axisPanel.axisType.isSampled:
        # TBD: look at again
        axisRegion = axisRegions[0]
      else:
        axisRegion = axisRegions[n]
      position_region[axisPanel.label] = axisRegion.region

    #print 'findPositionRegion', a, b, n, position_region

    return position_region

  def determinePosition(self, view, position):

    spectrum = view.analysisSpectrum.dataSource
    spectrum_position = spectrum.numDim * [0]
    spectrum_tile = spectrum.numDim * [0]

    for axisMapping in view.axisMappings:
      dataDim = axisMapping.analysisDataDim.dataDim
      dim = dataDim.dim - 1
      axisPanel = view.spectrumWindowPane.findFirstAxisPanel(label=axisMapping.label)
      if isinstance(dataDim, Nmr.FreqDataDim):
        p = convertPosition(position[axisMapping.label],
                            ExperimentBasic.getPrimaryDataDimRef(dataDim), 
                            fromUnit=axisPanel.axisUnit.unit)
        n = dataDim.numPointsOrig
        (tile, p) = divmod(p-1, n)
        spectrum_tile[dim] = int(tile)
        spectrum_position[dim] = p + 1.0
      else:
        spectrum_tile[dim] = 0
        spectrum_position[dim] = position[axisMapping.label]

    return (spectrum_position, spectrum_tile)

  def findNearbyPeak(self, peakList, region, xdim, ydim, xscale, yscale):

    if not (hasattr(peakList, 'cPeakList') and peakList.cPeakList):
      return None

    try:
      first = [ r[0] for r in region ]
      last  = [ r[1] for r in region ]
    except:
      # this can happen if region not set up correctly
      return None
    
    #print 'findNearbyPeak', self.windowPane.name, xdim, ydim, xscale, yscale, first, last
    allow_aliasing = getDimWrapped(peakList.dataSource)

    #peakInd = peakList.cPeakList.nearestPeak(xdim, ydim, xscale, yscale, first, last, allow_aliasing)
    (peakInd, d2Min) = peakList.cPeakList.nearestPeak(xdim, ydim, xscale, yscale, first, last)
    
    if peakInd >= 0:
      #return peakList.pickPeak(peakInd)
      return (peakList.sortedPeaks()[peakInd], d2Min)
    else:
      return None

  # returns dictionary keyed on label
  def findRegion(self, a0, b0, a1, b1, n):

    region = {}

    region['x'] = (a0, a1)
    region['y'] = (b0, b1)

    for axisPanel in self.windowPane.sortedAxisPanels()[2:]:
      axisRegions = axisPanel.sortedAxisRegions()
      if axisPanel.axisType.isSampled:
        # TBD: look at again
        axisRegion = axisRegions[0]
      else:
        axisRegion = axisRegions[n]
      region[axisPanel.label] = getAxisRegionRegion(axisRegion)

    return region

  def calcAxisMidpoint(self, axisRegion):

    (z0, z1) = axisRegion.region
    midpoint = 0.5 * (z0 + z1)
    if axisRegion.axisPanel.axisType.isSampled:
      midpoint = max(1.0, midpoint)
    return midpoint

  def getWorldRegion(self):

    # TBD: need to do better: should use panel.axisUnit region
    findFirstAxisPanel = self.windowPane.findFirstAxisPanel
    
    axisType = findFirstAxisPanel(label='x').axisType
    (t0, t1) = checkSwapRegion(axisType.region, axisType.findFirstAxisUnit(unit='ppm'))
    xregion = Region1D(t0, t1)

    axisType = findFirstAxisPanel(label='y').axisType
    (t0, t1) = checkSwapRegion(axisType.region, axisType.findFirstAxisUnit(unit='ppm'))
    yregion = Region1D(t0, t1)

    #print 'getWorldRegion', xregion, yregion

    return Region2D(xregion, yregion)

  def getViewRegions(self, panel):

    regions = []
    for axisRegion in panel.sortedAxisRegions():
      (t0, t1) = checkSwapRegion(axisRegion.region, panel.axisUnit)
      regions.append(Region1D(t0, t1))

    return regions

  def setViewRegions(self, regions, panel):

    for n in range(len(regions)):

      region = regions[n]
      axisRegion = panel.sortedAxisRegions()[n]
      (t0, t1) = checkSwapRegion(region, panel.axisUnit)
      #print 'setViewRegions', region, t0, t1
      axisRegion.region = (t0, t1)

  def getXviewRegions(self, worldRegion):

    xPanel  = self.windowPane.findFirstAxisPanel(label='x')
    regions = self.getViewRegions(xPanel)

    for region in regions:
      fitViewInWorld(region, worldRegion)

    self.setViewRegions(regions, xPanel)

    return regions

  def getYviewRegions(self, worldRegion):

    yPanel  = self.windowPane.findFirstAxisPanel(label='y')
    regions = self.getViewRegions(yPanel)

    for region in regions:
      fitViewInWorld(region, worldRegion)

    self.setViewRegions(regions, yPanel)

    return regions

  # convert position/region to points
  def convertPositionRegion(self, position_region, axisPanel, dataDim):

    label = axisPanel.label
    if label in ('x', 'y'):
      region = (position_region, position_region)
    else:
      region = position_region

    return convertRegion(region, axisPanel.axisUnit, dataDim)

  def determineDimRange(self, view, axisMapping, row, col, axisRegion = None):

    dataDim = axisMapping.analysisDataDim.dataDim
    axisPanel = self.windowPane.findFirstAxisPanel(label=axisMapping.label)

    #print 'determineDimRange', view.analysisSpectrum.dataSource.name, axisMapping.label
    if not axisRegion:
      if axisMapping.label == 'x':
        n = col
      elif axisMapping.label == 'y':
        n = row
      elif self.windowPane.spectrumWindow.stripAxis == 'x':
        n = col
      else:
        n = row

      axisRegion = axisPanel.sortedAxisRegions()[n]

    view_region = getAxisRegionRegion(axisRegion)
    axisType = axisPanel.axisType
    #print 'determineDimRange0', axisPanel.label, dataDim.dim, view_region

    if axisType.isSampled:
      npoints = dataDim.numPoints
      view_region = convertRegion(view_region, axisPanel.axisUnit, dataDim)
      world_region = (0.0, npoints)
      
    else:
      npoints = dataDim.numPointsOrig
      world_region = list(axisType.region)
      # TBD: assumes that axisType.findFirstAxisUnit() is same as min(max)AliasedFreq units
      dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
      expDimRef = dataDimRef.expDimRef
      minAliasedFreq = expDimRef.minAliasedFreq
      maxAliasedFreq = expDimRef.maxAliasedFreq
      if minAliasedFreq is None:
        minAliasedFreq = convertPosition(float(dataDim.numPoints), dataDimRef, toUnit=axisType.findFirstAxisUnit(unit='ppm').unit)
 
      if maxAliasedFreq is None:
        maxAliasedFreq = convertPosition(1.0, dataDimRef, toUnit=axisType.findFirstAxisUnit(unit='ppm').unit)

      world_region[0] = max(world_region[0], minAliasedFreq)
      world_region[1] = min(world_region[1], maxAliasedFreq)

      #print 'determineDimRange1', axisPanel.label, dataDim.dim, view_region, world_region

      if self.hasValueAxis and axisPanel.label == 'x':
        offset = getXAxisOffset(view)
        view_region = [view_region[0]-offset, view_region[1]-offset]

      view_region = convertRegion(view_region, axisPanel.axisUnit, dataDim)
      world_region = convertRegion(world_region, axisType.findFirstAxisUnit(unit='ppm'), dataDim)
      
      #print 'determineDimRange2', axisPanel.label, gdataDim.dim, view_region, world_region

    axisMapping.view_region = view_region
    axisMapping.world_region = world_region

    r0 = max(view_region[0], world_region[0])
    r1 = min(view_region[1], world_region[1])
    axisMapping.tile0 = int(r0) / npoints
    axisMapping.tile1 = int(r1+1) / npoints
    axisMapping.ntiles = axisMapping.tile1 - axisMapping.tile0 + 1

    #print 'determineDimRange3', axisPanel.label, dataDim.dim, r0, r1, \
    #         npoints, axisMapping.tile0, axisMapping.tile1, axisMapping.ntiles

  def determineTileRange(self, axisMapping, tile):

    tile = tile + axisMapping.tile0

    dataDim = axisMapping.analysisDataDim.dataDim

    if isinstance(dataDim, Nmr.FreqDataDim):
      # 23 Mar 09: looks like pointOffset does not belong here
      # so wp is relative to beginning of actual spectrum data
      ###wp = tile * dataDim.numPointsOrig + dataDim.pointOffset
      wp = tile * dataDim.numPointsOrig
    else:
      if tile != 0: # only draw fundamental region for pseudo-ND
        return None
      wp = tile * dataDim.numPoints
    wq = wp + dataDim.numPoints
    (r0, r1) = axisMapping.view_region
    (w0, w1) = axisMapping.world_region
    r0 = max(r0, w0)
    r1 = min(r1, w1)

    if wp > r1:
      return None

    if wq < r0:
      return None

    sp = 0
    if (wp < r0):
      sp = r0 - wp
      wp = r0

    sq =  dataDim.numPoints
    if (wq > r1):
      sq = sq + r1 - wq
      wq = r1

    worldPointRange = (wp, wq)
    spectrumPointRange = (sp, sq)

    #print 'determineTileRange', (sp, sq), (wp, wq)

    return (worldPointRange, spectrumPointRange)

  def drawViewTile(self, handler, object, view, contourLevels, contourStyle,
                   worldPointRanges, spectrumPointRanges, row, components=None):

    self.drawViewTileReal(handler, view, contourLevels, contourStyle,
                          worldPointRanges, spectrumPointRanges, row, components)

  def setupRanges(self, handler, view, worldPointRanges, spectrumPointRanges, row, pad = False):

    spectrum = view.analysisSpectrum.dataSource
    xaxisMapping = view.findFirstAxisMapping(label='x')
    
    xdim = xaxisMapping.analysisDataDim.dataDim.dim
    
    (a0, a1) = spectrumPointRanges[xdim-1]
    if a0 == a1:
      return None

    (x0, x1) = worldPointRanges[xdim-1]
    (r0, r1) = xaxisMapping.view_region
    x0 = (x0 - r0) / (r1 - r0)
    x1 = (x1 - r0) / (r1 - r0)
    if x0 == x1:
      return None

    if self.hasValueAxis:
      axisPanel = self.windowPane.findFirstAxisPanel(label='y')
      ###(b0, b1) = (y0, y1) = axisPanel.axisType.region
      # TBD: look at again, ### lines are wrong and below assumes only one y axis region
      axisRegion = axisPanel.sortedAxisRegions()[row]
      ###(r0, r1) = getAxisRegionRegion(axisRegion)
      (r0, r1) = (y0, y1) = (0.0, 1.0)
      (b0, b1) = getAxisRegionRegion(axisRegion)
      scale = spectrum.scale
      offset = getValueAxisOffset(view)
      b0 = b0 / scale - offset
      b1 = b1 / scale - offset
    else:
      yaxisMapping = view.findFirstAxisMapping(label='y')
      ydim = yaxisMapping.analysisDataDim.dataDim.dim
     
      (b0, b1) = spectrumPointRanges[ydim-1]
      (y0, y1) = worldPointRanges[ydim-1]
      (r0, r1) = yaxisMapping.view_region

    if b0 == b1:
      return None

    y0 = (y0 - r0) / (r1 - r0)
    y1 = (y1 - r0) / (r1 - r0)
    if y0 == y1:
      return None

    if not handler.makeCurrent():
      return None

    handler.mapRanges(x0, y0, x1, y1, a0, b0, a1, b1)

    if pad:
      params = getPeakFindParams(self.windowPane.root)
      thickness = params['thickness']

    ndim = spectrum.numDim
    firstFloat = ndim * [0]
    lastFloat = ndim * [0]
    for axisMapping in view.axisMappings:
      if not pad or axisMapping.label in ('x', 'y'):
        t = 0
      else:
        t = thickness + self.extraPad  # extra padding so that interpolated peaks shown when picked
      dim = axisMapping.analysisDataDim.dataDim.dim - 1
      firstFloat[dim] = spectrumPointRanges[dim][0] - t
      lastFloat[dim] = spectrumPointRanges[dim][1] + t

    if self.hasValueAxis:
      firstFloat.append(r0)
      lastFloat.append(r1)

    return (firstFloat, lastFloat)

  def drawViewTileReal(self, handler, view, contourLevels, contourStyle,
                       worldPointRanges, spectrumPointRanges, row, components=None):

    #print 'drawViewTileReal0'

    t = self.setupRanges(handler, view, worldPointRanges, spectrumPointRanges, row)
    if t is None:
      return
    (firstFloat, lastFloat) = t

    analysisSpectrum = view.analysisSpectrum
    spectrum = analysisSpectrum.dataSource
    ndim = spectrum.numDim
    #firstInt = [ int(math.ceil(x)) for x in firstFloat ]
    #lastInt = [ int(math.ceil(x)) for x in lastFloat ]
    #firstInt = [ int(math.floor(x)) for x in firstFloat ]
    #lastInt = [ 1+int(math.ceil(x)) for x in lastFloat ]
    firstInt = ndim * [0]
    lastInt = ndim * [0]
    for axisMapping in view.axisMappings:
      dataDim = axisMapping.analysisDataDim.dataDim
      dim = dataDim.dim - 1
      if axisMapping.label in ('x', 'y'):
        firstInt[dim] = int(math.floor(firstFloat[dim]))
        lastInt[dim] = 1 + int(math.ceil(lastFloat[dim]))
        
      else:
        firstInt[dim] = int(math.ceil(firstFloat[dim]))
        lastInt[dim] = 1 + int(math.floor(lastFloat[dim]))
        if lastInt[dim] == firstInt[dim]:
          p = int(firstFloat[dim] + 0.5)
          firstInt[dim] = p
          lastInt[dim] = p + 1
          
    for n in range(ndim):
      lastInt[n] = min(lastInt[n], spectrum.sortedDataDims()[n].numPoints)

    #print 'drawViewTileReal1', self.windowPane.name, firstInt, lastInt
    #print 'drawViewTileReal2', self.windowPane.name, firstFloat, lastFloat

    for n in range(ndim):
      if firstInt[n] >= lastInt[n]:
        break
        
    else:
      # some handlers (e.g. GL) are inherently C, some (e.g. PS)
      # have C handler attached, handler required here is C one
      if hasattr(handler, 'cHandler'):
        handler = handler.cHandler

      if self.hasValueAxis:
        if analysisSpectrum.posColors:
          n = len(analysisSpectrum.posColors) // 2
          color = analysisSpectrum.posColors[n]
        else:
          color = analysisSpectrum.sliceColor
        handler.setColor(hexToRgb(color))
        if components: # TBD: this is not needed but will allow non-C update
          # TBD: below is short-term hack
          if hasattr(analysisSpectrum, 'sliceColors'):
            sliceColors = analysisSpectrum.sliceColors
          else:
            sliceColors = (analysisSpectrum.sliceColor,)
          n = len(sliceColors)
          for (i, component) in enumerate(components):
            handler.setColor(hexToRgb(sliceColors[i%n]))
            view.sliceFile['x'].drawAll(handler, firstInt, lastInt, (component,))
        else:
          view.sliceFile['x'].drawAll(handler, firstInt, lastInt, None)
        
      else:
        usePrecalculated = analysisSpectrum.usePrecalculated

        contourFile = view.contourFile
        storedContourFiles = view.storedContourFiles

        if not view.contourFile and isWindowsOS():
          # NOTE:ED - this is needed because there are threading issues with Windows

          from ccpnmr.c.ContourFile import ContourFile
          from ccpnmr.analysis.core.BlockUtil import getBlockFile

          dataSource = view.analysisSpectrum.dataSource
          if dataSource.numDim >= 2:
            xMapping = view.findFirstAxisMapping(label='x')
            yMapping = view.findFirstAxisMapping(label='y')
            if xMapping and yMapping:
              xdim = xMapping.analysisDataDim.dataDim.dim - 1
              ydim = yMapping.analysisDataDim.dataDim.dim - 1
              if hasattr(dataSource, 'block_file') and dataSource.block_file:
                mem_cache = self.windowPopup.top.mem_cache

                writeable = hasattr(dataSource, 'writeable') and dataSource.writeable
                if not hasattr(self, '_tempContourBlockFile'):
                  self._tempContourBlockFile = getBlockFile(dataSource,
                                                            mem_cache=mem_cache,
                                                            writeable=writeable)
                if self._tempContourBlockFile:
                  contourFile = view.contourFile = ContourFile(xdim, ydim, self._tempContourBlockFile, mem_cache)

        if contourFile and (not storedContourFiles or not usePrecalculated):
          if components: # TBD: this is not needed but will allow non-C update
            contourFile.draw(handler, firstInt, lastInt, contourLevels, contourStyle, components)
          else:
            contourFile.draw(handler, firstInt, lastInt, contourLevels, contourStyle)

        if not contourFile or usePrecalculated:
          for storedContourFile in storedContourFiles:
            # TBD: components
            storedContourFile.draw(handler, firstInt, lastInt, contourLevels, contourStyle)

  def drawViewTilePeaks(self, handler, view, xscale, yscale, xdim, ydim,
                        worldPointRanges, spectrumPointRanges, row,
                        center, thickness, tile,
                        drawMethod, intensityMax, volumeMax, xpix, ypix):

    #print 'drawViewTilePeaks0', self.windowPane.name, center, thickness, tile
    #print 'drawViewTilePeaks1', self.windowPane.name, xscale, yscale, xdim, ydim
    t = self.setupRanges(handler, view, worldPointRanges, spectrumPointRanges, row, pad=True)
    if t is None:
      return
    (firstFloat, lastFloat) = t

    analysisSpectrum = view.analysisSpectrum
    spectrum = analysisSpectrum.dataSource
    font = analysisSpectrum.font
    (spectrumFontName, spectrumFontSize) = font.split()[:2]

    if hasattr(handler, 'cHandler'):
      handler = handler.cHandler

    profile = getAnalysisProfile(view.root)
    profileBgColor = hexToRgb(profile.bgColor)
    parent = self.parent
    while not hasattr(parent, 'setupCWinPeakList'):
      parent = parent.parent
    setupCWinPeakList = parent.setupCWinPeakList
    
    analysisPeakLists = set()
    for winPeakList in view.windowPeakLists:
      #print 'drawViewTilePeaks2:', winPeakList.analysisPeakList.peakList.serial, firstFloat, lastFloat, center, thickness
      if not hasattr(winPeakList, 'cWinPeakList'):
        # should not be here but looks like timing can sometimes
        # cause draw to happen before winPeakList is setup
        setupCWinPeakList(winPeakList)
        # this could be called if exiting program so check that ok
        if not hasattr(winPeakList, 'cWinPeakList'):
          return
	
      if self.isWinPeakListDrawn(winPeakList):
        if hasattr(winPeakList, 'printFont'):
          if winPeakList.printPeaks:
            winPeakList.cWinPeakList.setIsSymbolDrawn(1)
          else:
            winPeakList.cWinPeakList.setIsSymbolDrawn(0)
          font = winPeakList.printFont
          if font == no_peak_text:
            winPeakList.cWinPeakList.setIsTextDrawn(0)
          else:
            winPeakList.cWinPeakList.setIsTextDrawn(1)
            (name, size) = font.split()[:2]
            self.setHandlerFont(handler, name, size)
          analysisPeakList = winPeakList.analysisPeakList
          if analysisPeakList not in analysisPeakLists:
            color = analysisPeakList.symbolColor              # printColor
            if color and color[0] != '#': # assume it is a name
              scheme = profile.findFirstColorScheme(name=color)
              if scheme:
                color = scheme.colors[0]
            analysisPeakList.peakList.cPeakList.setColor(hexToRgb(color))
            analysisPeakLists.add(analysisPeakList)
          bgColor = (1.0, 1.0, 1.0)
        else:
          bgColor = profileBgColor
          self.setHandlerFont(handler, spectrumFontName, spectrumFontSize)

        # NOTE:ED - sometimes the font is defined as a strokeFont which gives a freeGlut error
        if hasattr(handler, 'hasFont'):
          thisFont = handler.hasFont()
          if thisFont and thisFont in (('Roman', 1), ('MonoRoman', 1)):
            # CHECK THIS CONDITION
            self.setHandlerFont(handler, spectrumFontName, spectrumFontSize)

        winPeakList.cWinPeakList.drawPeaks(handler, xdim, ydim, xpix, ypix,
                                 xscale, yscale, firstFloat, lastFloat,
                                 drawMethod, intensityMax, volumeMax,
                                 bgColor, center, thickness, tile)
        if hasattr(winPeakList, 'printFont'):
          winPeakList.cWinPeakList.setIsSymbolDrawn(winPeakList.isSymbolDrawn)
          winPeakList.cWinPeakList.setIsTextDrawn(winPeakList.isAnnotationDrawn)
      #print 'drawViewTilePeaks3'
    
    for analysisPeakList in analysisPeakLists:
      color = analysisPeakList.symbolColor
      analysisPeakList.peakList.cPeakList.setColor(hexToRgb(color))

  def calcPeakScale(self, axisPanel, axisMapping, size, fromUnit, toUnit='point'):

    #print 'calcPeakScale', axisPanel.label, size, fromUnit

    dataDimRef = ExperimentBasic.getPrimaryDataDimRef(axisMapping.analysisDataDim.dataDim)
    p = convertPosition(size, dataDimRef, fromUnit=fromUnit, toUnit=toUnit, relative=True)

    return abs(p)

  def getPeakScale(self, view, xAxisRegion, yAxisRegion, xDataDim=None, yDataDim=None):

    if self.hasValueAxis:
      return (1, 1)

    #spectrum = view.analysisSpectrum.dataSource
    xaxisPanel = self.windowPane.findFirstAxisPanel(label='x')
    yaxisPanel = self.windowPane.findFirstAxisPanel(label='y')

    xaxisMapping = view.findFirstAxisMapping(label='x')
    yaxisMapping = view.findFirstAxisMapping(label='y')
    
    analysisProject = view.topObject
    drawMethod = analysisProject.peakDrawMethod
    toUnit = 'point'
    if drawMethod == peak_draw_methods[0]: # uniform in pixels
      p = analysisProject.peakPixelSize
      (r0, r1) = xAxisRegion.region
      px = p * abs(r1 - r0) / float(xAxisRegion.size)
      (r0, r1) = yAxisRegion.region
      py = p * abs(r1 - r0) / float(yAxisRegion.size)

      xFromUnit = xaxisPanel.axisUnit.unit
      yFromUnit = yaxisPanel.axisUnit.unit
    elif drawMethod == peak_draw_methods[4]: # line width
      px = py = 2.0 # changed from 1.0 on 7 Oct 2014 because was showing double line width
      xFromUnit = yFromUnit = 'point'
      toUnit = 'Hz'
    elif drawMethod == peak_draw_methods[5]: # box width
      # added 0.5 factor on 7 Oct 2014 because was showing double box width
      px = 0.5 * getPeakFindBoxwidth(xDataDim) if xDataDim else 1.0
      py = 0.5 * getPeakFindBoxwidth(yDataDim) if yDataDim else 1.0
      xFromUnit = yFromUnit = 'point'
    else:
      (r0, r1) = xAxisRegion.region
      px = xaxisPanel.axisType.peakSize
      sx = px * float(xAxisRegion.size) / abs(r1 - r0)
      (r0, r1) = yAxisRegion.region
      py = yaxisPanel.axisType.peakSize
      sy = py * float(yAxisRegion.size) / abs(r1 - r0)

      # scale so that peak symbol is square
      if (sx > sy):
        py = py * sx / sy
      else:
        px = px * sy / sx

      xFromUnit = yFromUnit = 'ppm' # peak size always specified in ppm

    xscale = self.calcPeakScale(xaxisPanel, xaxisMapping, px, xFromUnit, toUnit)
    yscale = self.calcPeakScale(yaxisPanel, yaxisMapping, py, yFromUnit, toUnit)

    if drawMethod == peak_draw_methods[4]: # line width
      xscale = 1.0 / xscale
      yscale = 1.0 / yscale

    return (xscale, yscale)

  def getSampleRegions(self):

    axisPanels = self.windowPane.sortedAxisPanels()[2:]
    sampleRegions = []
    for axisPanel in axisPanels:
      if axisPanel.axisType.isSampled:
        n = len(axisPanel.axisRegions)
      else:
        n = 1
      sampleRegions.append(n)

    return sampleRegions

  def viewNotReadyYet(self, view):

    # check if axisMappings not ready yet

    if not view.axisMappings:
      return True

    for label in ('x', 'y'): # should not need this check but seem to get exceptions when axisMapping is None in y
      if self.hasValueAxis and label == 'y':
        continue
      if not view.findFirstAxisMapping(label=label):
        return True

    if len(view.axisMappings) != view.analysisSpectrum.dataSource.numDim:
      return True

    return False

  def drawView(self, handler, object, view, row, col):

    #print 'drawViewA', self.windowPane.name

    # # What is below doing? TJS  (if both of these are not set then nothing can draw)
    # if not self.hasValueAxis and \
    #     (not hasattr(view, 'contourFile') or not view.contourFile) and \
    #     (not hasattr(view, 'storedContourFiles') or not view.storedContourFiles):
    #   return

    if self.viewNotReadyYet(view):
      return

    analysisSpectrum = view.analysisSpectrum
    spectrum = analysisSpectrum.dataSource

    block_file = hasattr(spectrum, 'block_file') and spectrum.block_file
    storedContourFiles = hasattr(view, 'storedContourFiles') and view.storedContourFiles
    if not block_file and not storedContourFiles:
      return

    # TBD: more general cases
    """
    for dataDim in spectrum.dataDims:
      assert isinstance(dataDim, Nmr.FreqDataDim), 'Spectrum %s:%s has non FreqDataDim %d' % (spectrum.experiment.name, spectrum.name, dim+1)
      assert len(dataDim.dataDimRefs) == 1, 'Spectrum %s:%s dim %d should have exactly 1 dataDimRef, has %d' % (spectrum.experiment.name, spectrum.name, 1, len(dataDim.dataDimRefs))
    """

    project = spectrum.root
    scale = spectrum.scale / self.windowPane.spectrumWindow.analysisProject.globalContourScale

    if self.hasValueAxis:
      if not view.isSliceVisible:
        return
      contourLevels = contourStyle = None
      
    else:
      posLevels = [ abs(l/scale) for l in analysisSpectrum.posLevels ]
      negLevels = [-abs(l/scale) for l in analysisSpectrum.negLevels ]

      levels = []
      if self.isViewPosVisible(view):
        levels += posLevels
        
      if self.isViewNegVisible(view):
        levels += negLevels

      if not levels:
        return

      contourLevels = ContourLevels.ContourLevels(levels)

      # TBD: switch positive and negative if isAliased
      #print 'drawViewB', self.windowPane.name
      
      if hasattr(analysisSpectrum, 'printPositiveColors'):
        posColors = analysisSpectrum.printPositiveColors
      else:
        posColors = analysisSpectrum.posColors
      if hasattr(analysisSpectrum, 'printNegativeColors'):
        negColors = analysisSpectrum.printNegativeColors
      else:
        negColors = analysisSpectrum.negColors
      
      rgbPos = [hexToRgb(c) for c in posColors]
      rgbNeg = [hexToRgb(c) for c in negColors]
      
      contourStyle = ContourStyle.ContourStyle(rgbPos, rgbNeg, 0, 0)

    # TBD: this is a bit of a hack until in data model
    if hasattr(analysisSpectrum, 'components'):
      components = analysisSpectrum.components
    else:
      components = None

    # determine view region in points and tiles
    xdim = ydim = -1
    ndim = spectrum.numDim
    ntiles_array = ndim * [0]

    sampleRegions = self.getSampleRegions()
    (nsampleRegions, cumSampleRegions) = cumulativeProductArray(sampleRegions)

    for sample in range(nsampleRegions):
      sampleRegion = arrayOfIndex(sample, cumSampleRegions)
      dd = {}
      for k in range(len(sampleRegion)):
        dd['z%d' % (k+1)] = sampleRegion[k]

      for analysisDataDim in analysisSpectrum.analysisDataDims:
        dataDim = analysisDataDim.dataDim
        dim = dataDim.dim - 1
        axisMapping = view.findFirstAxisMapping(analysisDataDim=analysisDataDim)
        label = axisMapping.label
        axisPanel = self.windowPane.findFirstAxisPanel(label=label)
        if axisPanel.axisType.isSampled:
          axisRegion = axisPanel.sortedAxisRegions()[dd[label]]
        else:
          axisRegion = None
        self.determineDimRange(view, axisMapping, row, col, axisRegion)
        ntiles_array[dim] = axisMapping.ntiles

      #print 'drawViewC', self.windowPane.name, ntiles_array
      (ntiles, cum_tiles) = cumulativeProductArray(ntiles_array)
      #print 'drawViewD', self.windowPane.name, ntiles, cum_tiles
    
      drawViewTile = self.drawViewTile
      analysisSpectrum = view.analysisSpectrum
      for tile in range(ntiles):
        worldPointRanges = ndim * [0]
        spectrumPointRanges = ndim * [0]
        tile_array = arrayOfIndex(tile, cum_tiles)
        #print 'drawView0', self.windowPane.name, tile, tile_array
      
        for analysisDataDim in analysisSpectrum.analysisDataDims:
          dataDim = analysisDataDim.dataDim
          dim = dataDim.dim - 1
          axisMapping = view.findFirstAxisMapping(analysisDataDim=analysisDataDim)
          result = self.determineTileRange(axisMapping, tile_array[dim])
          if not result:
            break # this tile no good
          
          (worldPointRanges[dim], spectrumPointRanges[dim]) = result
          #print 'drawViewA0', self.windowPane.name, dim, worldPointRanges[dim], spectrumPointRanges[dim]
        
        else: # all dims ok so this tile good
          #print 'drawView2: about to draw tile:',  self.windowPane.name, str(tile_array)
          # TBD: below a kludge, see if can do better
          if (not object or \
              (hasattr(object, 'doubleBuffer') and object.doubleBuffer)):
            drawViewTile(handler, object, view, contourLevels, contourStyle,
                         worldPointRanges, spectrumPointRanges, row, components)
          else:
            self.after_idle(lambda contourLevels=contourLevels, contourStyle=contourStyle,
                            worldPointRanges=worldPointRanges, spectrumPointRanges=spectrumPointRanges: \
                            drawViewTile(handler, object, view, contourLevels, contourStyle,
                                         worldPointRanges, spectrumPointRanges, row, components))
      #print 'drawView3', self.windowPane.name
    #print 'drawView4', self.windowPane.name

  def findFundamentalRegion(self, view, label, r0, r1):

    axisMapping = view.findFirstAxisMapping(label=label)
    dataDim = axisMapping.analysisDataDim.dataDim
    assert isinstance(dataDim, Nmr.FreqDataDim), 'Need FreqDataDim for box'

    dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
    axisPanel = view.spectrumWindowPane.findFirstAxisPanel(label=axisMapping.label)
    toUnit = axisPanel.axisUnit.unit
    # 23 Mar 09: below was just plain wrong, but it looks like point 1 is
    # relative to the actual data, so point -pointOffset is where the left
    # side of the box should be, and the right side should be at numPointsOrig-pointOffset
    # (the whole box is numPointsOrig wide)
    ###r = convertPosition(dataDim.pointOffset+0.5, dataDimRef, toUnit=toUnit)
    r = convertPosition(-dataDim.pointOffset+0.5, dataDimRef, toUnit=toUnit)
    t0 = float(r - r0) / (r1 - r0)
    ###r = convertPosition(dataDim.numPointsOrig+0.5, dataDimRef, toUnit=toUnit)
    r = convertPosition(dataDim.numPointsOrig-dataDim.pointOffset+0.5, dataDimRef, toUnit=toUnit)
    t1 = float(r - r0) / (r1 - r0)

    return (t0, t1)

  def drawViewBox(self, handler, view, x0, x1, y0, y1):

    if self.hasValueAxis:
      return

    if self.viewNotReadyYet(view):
      return

    analysisSpectrum = view.analysisSpectrum

    if not handler.makeCurrent():
      return

    handler.mapRanges(0, 0, 1, 1, 0, 0, 1, 1)

    if hasattr(handler, 'cHandler'):
      handler = handler.cHandler

    colors = analysisSpectrum.posColors
    color = colors[int(len(colors)/2)]
    handler.setColor(hexToRgb(color))

    (a0, a1) = self.findFundamentalRegion(view, 'x', x0, x1)
    (b0, b1) = self.findFundamentalRegion(view, 'y', y0, y1)
    handler.drawDashBox(a0, b0, a1, b1)

  def setHandlerFont(self, handler, name, size):

    # TBD: remove below when Font is ComplexDataType in data model
    size = int(size)

    if isinstance(handler, PrintHandler) \
     and (isinstance(handler.outputHandler, Pdf) \
     or isinstance(handler.outputHandler, PostScript)):
      # need to call C world setFont instead of Python world one
      size = max(1, int(handler.output_scaling * size + 0.5))
      #size = max(1, int(max(handler.output_scaling) * size + 0.5))
      handler.cHandler.setFont(name, size)
      
    else:
      handler.setFont(name, size)

  def drawViewPeaks(self, handler, view, row, col):

    if self.viewNotReadyYet(view):
      return

    analysisSpectrum = view.analysisSpectrum
    try:
      spectrum = analysisSpectrum.dataSource
      if spectrum.isDeleted:
        return
        
    except:
      return # spectrum has been deleted

    xAxisPanel = self.windowPane.findFirstAxisPanel(label='x')
    yAxisPanel = self.windowPane.findFirstAxisPanel(label='y')
    xAxisRegion = xAxisPanel.sortedAxisRegions()[col]
    yAxisRegion = yAxisPanel.sortedAxisRegions()[row]

    project = analysisSpectrum.root
    analysisProject = analysisSpectrum.analysisProject
    drawMethod = getCPeakDrawMethod(project)
    intensityMax = analysisProject.peakIntensityScale
    volumeMax = analysisProject.peakVolumeScale

    xaxisMapping = view.findFirstAxisMapping(label='x')
    xDataDim = xaxisMapping.analysisDataDim.dataDim
    xdim = xDataDim.dim - 1
    if self.hasValueAxis:
      yDataDim = None
      ydim = -1
      
    else:
      yaxisMapping = view.findFirstAxisMapping(label='y')
      yDataDim = yaxisMapping.analysisDataDim.dataDim
      ydim = yDataDim.dim - 1

    (xscale, yscale) = self.getPeakScale(view, xAxisRegion, yAxisRegion, xDataDim, yDataDim)

    xpix = calcPointsPerPixel(view, 'x')
    ypix = calcPointsPerPixel(view, 'y')

    analysisSpectrum = view.analysisSpectrum
    ndim = spectrum.numDim
    center = ndim * [0]
    thickness = ndim * [0]
    ntiles_array = ndim * [0]
    tile_base = ndim * [0]

    sampleRegions = self.getSampleRegions()
    (nsampleRegions, cumSampleRegions) = cumulativeProductArray(sampleRegions)

    for sample in range(nsampleRegions):
      sampleRegion = arrayOfIndex(sample, cumSampleRegions)
      dd = {}
      for k in range(len(sampleRegion)):
        dd['z%d' % (k+1)] = sampleRegion[k]

      # determine view region in points and tiles
      for analysisDataDim in analysisSpectrum.analysisDataDims:
        dataDim = analysisDataDim.dataDim
        dim = dataDim.dim - 1
        axisMapping = view.findFirstAxisMapping(analysisDataDim=analysisDataDim)
        label = axisMapping.label
        axisPanel = self.windowPane.findFirstAxisPanel(label=label)
        if axisPanel.axisType.isSampled:
          axisRegion = axisPanel.sortedAxisRegions()[dd[label]]
        else:
          axisRegion = None
        self.determineDimRange(view, axisMapping, row, col, axisRegion)
        ntiles_array[dim] = axisMapping.ntiles
        tile_base[dim] = axisMapping.tile0
        r = axisMapping.view_region
        center[dim] = 1.0 + 0.5 * (r[0] + r[1])
        thickness[dim] = 0.5 * (r[1] - r[0])

      # cannot do this here because in Postscript needs to be done after setting up range
      ##self.setHandlerFont(handler, name, size)

      (ntiles, cum_tiles) = cumulativeProductArray(ntiles_array)
    
      for tile in range(ntiles):
        worldPointRanges = ndim * [0]
        spectrumPointRanges = ndim * [0]
        tile_array = arrayOfIndex(tile, cum_tiles)
        #print 'drawViewPeaks0', self.windowPane.name, tile, tile_array
      
        for analysisDataDim in analysisSpectrum.analysisDataDims:
          dataDim = analysisDataDim.dataDim
          dim = dataDim.dim - 1
          axisMapping = view.findFirstAxisMapping(analysisDataDim=analysisDataDim)
          result = self.determineTileRange(axisMapping, tile_array[dim])
          if not result:
            break # this tile no good
          
          (worldPointRanges[dim], spectrumPointRanges[dim]) = result
        
        else: # all dims ok so this tile good
          alias_tile = ndim * [0]
          for dim in range(ndim):
            alias_tile[dim] = tile_base[dim] + tile_array[dim]
          self.drawViewTilePeaks(handler, view, xscale, yscale, xdim, ydim,
                                 worldPointRanges, spectrumPointRanges, row,
                                 center, thickness, alias_tile,
                                 drawMethod, intensityMax, volumeMax, xpix, ypix)

  ###def drawPeakClusters(self, handler, row, col):
  ###
  ###  pass

  def drawMarkDimRuler(self, handler, object, axisLabel, dashLength, gapLength,
                       x0, x1, y0, y1):

    # TBD: bit of a hack this, scale dashLength and gapLength for PS output
    if hasattr(handler, 'cHandler'):
      scale = max(1, 10/dashLength)
      dashLength = scale * dashLength
      gapLength = scale * gapLength

    if axisLabel == 'x':
      x = float(object.position - x0) / (x1 - x0)
      #print 'drawMarkDimRuler x', x0, x1, x
      handler.drawDashLine(x, 0, x, 1, dashLength, gapLength)
      
    elif axisLabel == 'y':
      y = float(object.position - y0) / (y1 - y0)
      #print 'drawMarkDimRuler y', y0, y1, y
      handler.drawDashLine(0, y, 1, y, dashLength, gapLength)

  def drawMark(self, handler, mark, x0, x1, y0, y1):

    worldPanels = self.windowPane.sortedAxisPanels()[:2]

    color = mark.color
    handler.setColor(hexToRgb(color))
    handler.setLineWidth(mark.lineWidth)
    dashLength = mark.dashLength
    gapLength = mark.gapLength

    drawMarkDimRuler = self.drawMarkDimRuler
    for markDim in mark.markDims:
      for worldPanel in worldPanels:
        if markDim.axisType == worldPanel.axisType:
          drawMarkDimRuler(handler, markDim, worldPanel.label,
                           dashLength, gapLength, x0, x1, y0, y1)
      
  def drawMarks(self, handler, x0, x1, y0, y1):

    if not handler.makeCurrent():
      return

    handler.mapRanges(0, 0, 1, 1, 0, 0, 1, 1)
    drawMark = self.drawMark
    
    for mark in self.windowPane.topObject.sortedMarks():
      drawMark(handler, mark, x0, x1, y0, y1)

    # handler.resetLineWidth() # reset to default
    # NOTE:ED - quickest way to fix the bug without re-compiling
    #           error caused because tkhandler DEFAULT_LINEWIDTH is 0.0, should be 1.0 under conda
    handler.setLineWidth(1.0)

  def drawRuler(self, handler, ruler, x0, x1, y0, y1):

    worldPanels = self.windowPane.sortedAxisPanels()[:2]

    color = ruler.color
    handler.setColor(hexToRgb(color))
    handler.setLineWidth(ruler.lineWidth)
    dashLength = ruler.dashLength
    gapLength = ruler.gapLength

    drawMarkDimRuler = self.drawMarkDimRuler
    for worldPanel in worldPanels:
      if (ruler.panelType == worldPanel.panelType):
        drawMarkDimRuler(handler, ruler, worldPanel.label,
                         dashLength, gapLength, x0, x1, y0, y1)
      
  def drawRulers(self, handler, x0, x1, y0, y1):

    if not handler.makeCurrent():
      return

    handler.mapRanges(0, 0, 1, 1, 0, 0, 1, 1)
    drawRuler = self.drawRuler
    
    for ruler in self.windowPane.topObject.sortedRulers():
      drawRuler(handler, ruler, x0, x1, y0, y1)

    # handler.resetLineWidth() # reset to default
    handler.setLineWidth(1.0)

  def drawDeltaMarker(self, handler, x0, x1, y0, y1, color):

    if not self.deltaPositions:
      return

    if not handler.makeCurrent():
      return

    handler.mapRanges(0, 0, 1, 1, 0, 0, 1, 1)
    handler.setColor(hexToRgb(color))

    (x, y) = self.deltaPositions
    x = float(x - x0) / (x1 - x0)
    handler.drawLine(x, 0, x, 1)

    y = float(y - y0) / (y1 - y0)
    handler.drawLine(0, y, 1, y)

  def drawDiagonal(self, handler, x0, x1, y0, y1, color, isDashed=False):

    if not handler.makeCurrent():
      return

    handler.mapRanges(0, 0, 1, 1, 0, 0, 1, 1)
    handler.setColor(hexToRgb(color))

    yy0 = float(x0 - y0) / (y1 - y0)
    yy1 = float(x1 - y0) / (y1 - y0)
    #print 'drawDiagonal', x0, x1, y0, y1, yy0, yy1
    if isDashed:
      handler.drawDashLine(0, yy0, 1, yy1, 2, 2)
    else:
      handler.drawLine(0, yy0, 1, yy1)

  def drawZeroLine(self, handler, y0, y1, color):

    if not handler.makeCurrent():
      return

    handler.mapRanges(0, 0, 1, 1, 0, 0, 1, 1)
    handler.setColor(hexToRgb(color))

    y = -float(y0) / (y1-y0)
    #print 'drawZeroLine', y0, y1, y
    handler.drawLine(0, y, 1, y)

  def drawLowestContourLine(self, handler, view, y0, y1, haveHigh, haveLow):

    if not handler.makeCurrent():
      return

    analysisSpectrum = view.analysisSpectrum
    spectrum = analysisSpectrum.dataSource
    analysisProject = analysisSpectrum.analysisProject

    ##scale = spectrum.scale / analysisProject.globalContourScale
    # instead of scaling contour line down by spectrum.scale, scale spectrum drawing up
    scale = 1.0 / analysisProject.globalContourScale
    
    posLevels = [ l/scale for l in analysisSpectrum.posLevels ]
    negLevels = [ l/scale for l in analysisSpectrum.negLevels ]
    
    if not (posLevels or negLevels):
      return

    handler.mapRanges(0, 0, 1, 1, 0, 0, 1, 1)
    if analysisSpectrum.posColors:
      n = len(analysisSpectrum.posColors) // 2
      color = analysisSpectrum.posColors[n]
    else:
      color = analysisSpectrum.sliceColor
    handler.setColor(hexToRgb(color))

    dashLength = 1
    gapLength = 3

    if haveHigh:
      if posLevels:
        y = min(posLevels)
        y = float(y-y0) / (y1-y0)
        #handler.drawLine(0, y, 1, y)
        handler.drawDashLine(0, y, 1, y, dashLength, gapLength)

    if haveLow:
      if negLevels:
        y = max(negLevels)
        y = float(y-y0) / (y1-y0)
        #handler.drawLine(0, y, 1, y)
        handler.drawDashLine(0, y, 1, y, dashLength, gapLength)

  def drawCanvasLabel(self, handler, strip, x0, x1, y0, y1, color):

    if not handler.makeCurrent():
      return

    windowPane = self.windowPane
    if windowPane.spectrumWindow.stripAxis == 'x':
      axisPanel   = windowPane.findFirstAxisPanel(label='x')
    else:
      axisPanel   = windowPane.findFirstAxisPanel(label='y')
    axisRegion = axisPanel.sortedAxisRegions()[strip]
    if axisRegion.isActive:
      ss = '*'
    else:
      ss = ''

    handler.mapRanges(0, 0, 1, 1, 0, 0, 1, 1)
    self.setHandlerFont(handler, name='Helvetica', size=12)
    handler.setColor( hexToRgb(color) )
    handler.drawText('%d%s' % (strip+1,ss), 0.02, 0.99, 0, 1) # top left corner

  def drawCanvasMidpoint(self, handler, strip, x0, x1, y0, y1, color):

    if not handler.makeCurrent():
      return

    handler.mapRanges(0, 0, 1, 1, 0, 0, 1, 1)
    self.setHandlerFont(handler, name='Helvetica', size=12)
    handler.setColor( hexToRgb(color) )
    texts = []
    for axisPanel in self.windowPane.sortedAxisPanels()[2:]:
      axisType = axisPanel.axisType
      if axisType and not axisType.isSampled:
        #region = axisPanel.sortedAxisRegions()[strip].region
        axisRegion = axisPanel.sortedAxisRegions()[strip]
        region = getAxisRegionRegion(axisRegion)
        midpoint = 0.5 * (region[0] + region[1])
        text = formatDecimals(midpoint, decimals=axisType.numDecimals)
        texts.append(text)
        
    if texts:
      text = ','.join(texts)
      handler.drawText(text, 0.02, 0.01, 0, 0) # bottom left corner

  def findAxisRegion(self, axisPanel, n):

    axisRegion = axisPanel.sortedAxisRegions()[n]
    (t0, t1) = getAxisRegionRegion(axisRegion)
    if axisPanel.axisUnit and axisPanel.axisUnit.isBackwards:
      (t0, t1) = (t1, t0)
 
    return (t0, t1)

  # used to sort view list by order
  def compareViewOrder(self, view1, view2):

    analysisSpectrum1 = view1.analysisSpectrum
    analysisSpectrum2 = view2.analysisSpectrum

    c = cmp(analysisSpectrum1.rank,analysisSpectrum2.rank)
    if c != 0:
      return c

    spectrum1 = analysisSpectrum1.dataSource
    spectrum2 = analysisSpectrum2.dataSource

    expt1 = spectrum1.experiment
    expt2 = spectrum2.experiment

    c = cmp(expt1.name, expt2.name)
    if c != 0:
      return c

    return cmp(spectrum1.name, spectrum2.name)

  def doCanvas(self, handler, object = None, row = 0, col = 0):

    windowPane = self.windowPane
    window = windowPane.spectrumWindow
    allViews = list(self.getSpectrumViews())
    allViews.sort(self.compareViewOrder)
    allViews.reverse()

    # print 'doCanvas1'
    for view in allViews:
      #if (view.analysisSpectrum.dataSource.numDim >= 2):
      self.drawView(handler, object, view, row, col)

    # print 'doCanvas2'
    for view in allViews:
      #if (view.analysisSpectrum.dataSource.numDim >= 2):
      self.drawViewPeaks(handler, view, row, col)

    ###self.drawPeakClusters(handler, row, col)

    # print 'doCanvas3'
    xPanel = windowPane.findFirstAxisPanel(label='x')
    yPanel = windowPane.findFirstAxisPanel(label='y')
    (x0, x1) = self.findAxisRegion(xPanel, col)
    (y0, y1) = self.findAxisRegion(yPanel, row)

    if not self.hasValueAxis:
      for view in allViews:
        analysisSpectrum = view.analysisSpectrum
        spectrum = analysisSpectrum.dataSource
        if spectrum.numDim >= 2:
          if self.isViewVisible(view) and analysisSpectrum.useBoundingBox:
            self.drawViewBox(handler, view, x0, x1, y0, y1)

    # print 'doCanvas4'
    self.drawMarks(handler, x0, x1, y0, y1)
    self.drawRulers(handler, x0, x1, y0, y1)

    project = window.root
    profile = getAnalysisProfile(project)
    color = inverseGrey(profile.bgColor)

    self.drawDeltaMarker(handler, x0, x1, y0, y1, color)

    xaxisType = xPanel.axisType
    yaxisType = yPanel.axisType

    #print 'doCanvas5'
    if xaxisType == yaxisType:
      self.drawDiagonal(handler, x0, x1, y0, y1, color)

      # pseudo-diagonals
      # TBD: assume for now that have ppm
      if xPanel.axisUnit.unit == yPanel.axisUnit.unit == 'ppm':
        for view in allViews:
          if view.isPosVisible or view.isNegVisible:
            analysisSpectrum = view.analysisSpectrum
            spectrum = analysisSpectrum.dataSource
            experiment = spectrum.experiment
            spinningRate = experiment.spinningRate
            if spinningRate:
              dataDim = view.findFirstAxisMapping(label='x').analysisDataDim.dataDim
              dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
              expDimRef = dataDimRef.expDimRef
              spinningRate /= expDimRef.sf  # assumes y expDimRef would give the same
              nmin = int((y1-x0) // spinningRate)
              nmax = - int((x1-y0) // spinningRate)
              for n in range(nmin, nmax+1):
                if n:  # n = 0 is normal diagonal
                  self.drawDiagonal(handler, x0+n*spinningRate, x1+n*spinningRate, y0, y1, color, isDashed=True)

    # extra multiple-quantum diagonals
    if xaxisType.isotopeCodes == yaxisType.isotopeCodes:
      if xaxisType.measurementType == 'MQShift' and yaxisType.measurementType == 'Shift':
        self.drawDiagonal(handler, x0, x1, 2*y0, 2*y1, color)
      elif xaxisType.measurementType == 'Shift' and yaxisType.measurementType == 'MQShift':
        self.drawDiagonal(handler, 2*x0, 2*x1, y0, y1, color)

    #print 'doCanvas6'
    if self.hasValueAxis and window.isZeroLineShown:
      self.drawZeroLine(handler, y0, y1, color)

    if self.hasValueAxis:
      # a bit of a hack this, only want lowest contour line on screen, not on paper
      if self.__class__.__name__ != 'WindowDraw':
        peakFindParams = getPeakFindParams(project)
        haveHigh = peakFindParams['haveHigh']
        haveLow = peakFindParams['haveLow']
        for view in allViews:
          if view.isSliceVisible and view.isContourLineVisible:
            self.drawLowestContourLine(handler, view, y0, y1, haveHigh, haveLow)
    #print 'doCanvas7'

    strip = -1
    if window.stripAxis == 'x':
      if len(xPanel.axisRegions) > 1:
        strip = col
    else:
      if len(yPanel.axisRegions) > 1:
        strip = row

    if strip >= 0:
      if window.isCanvasLabelShown or window.isCanvasMidpointShown:
        if window.isCanvasLabelShown:
          self.drawCanvasLabel(handler, strip, x0, x1, y0, y1, color)
        if window.isCanvasMidpointShown:
          self.drawCanvasMidpoint(handler, strip, x0, x1, y0, y1, color)

  def drawCanvas(self, handler, object, row, col):
 
    self.drawCanvasReal(handler, object, row, col)

  def drawCanvasReal(self, handler, object, row, col):

    self.doCanvas(handler, object, row, col)

  def drawRow(self, handlers, object, row):

    ncols = self.getNCols()
    #print 'WindowDraw: drawRow1', handlers, row, ncols
    drawCanvas = self.drawCanvas
    
    for i in range(ncols):
      drawCanvas(handlers[row][i], object, row, i)

    #print 'WindowDraw: drawRow2'

  def drawCol(self, handlers, object, col):

    nrows = self.getNRows()
    drawCanvas = self.drawCanvas
    for j in range(nrows):
      drawCanvas(handlers[j][col], object, j, col)

  # object is there in case need extra information
  # handlers is what contains the drawing funtionality
  def drawAll(self, handlers, object = None):

    nrows = self.getNRows()
    drawRow = self.drawRow
    for j in range(nrows):
      drawRow(handlers, object, j)

  def deleteCol(self, col = -1):

    axisPanels = self.windowPane.sortedAxisPanels()
    axisPanels[0].sortedAxisRegions()[col].delete()
    # remove corresponding z axisRegions
    for axisPanel in axisPanels[2:]:
      ###if self.windowPane.spectrumWindow.stripAxis == 'x':
      if self.windowPane.spectrumWindow.stripAxis == 'x' and isIsotopeAxisType(axisPanel.axisType):
        axisPanel.sortedAxisRegions()[col].delete()

  def deleteRow(self, row=-1):

    axisPanels = self.windowPane.sortedAxisPanels()
    axisPanels[1].sortedAxisRegions()[row].delete()
    # remove corresponding z axisRegions
    for axisPanel in axisPanels[2:]:
      ###if self.windowPane.spectrumWindow.stripAxis == 'y':
      if self.windowPane.spectrumWindow.stripAxis == 'y' and isIsotopeAxisType(axisPanel.axisType):
        axisPanel.sortedAxisRegions()[row].delete()

  def drawViewSlicePeaks(self, slice, view, thickness):

    if self.viewNotReadyYet(view):
      return

    thickness = float(thickness)
    windowPane = self.windowPane
    analysisProject = windowPane.topObject
    handler = slice.handler
    
    analysisSpectrum = view.analysisSpectrum
    (name, size) = analysisSpectrum.font.split()[:2]
    self.setHandlerFont(handler, name, size)
    
    (r0, r1) = slice.view_region
    handler.mapRanges(0, 0, 1, 1, r0, 0, r1, 1)
    # note that r0 > r1

    unit = windowPane.findFirstAxisPanel(label='x').axisUnit.unit
    axisMapping = view.findFirstAxisMapping(label='x')
    dataDim = axisMapping.analysisDataDim.dataDim
    dim = dataDim.dim
    dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
    numPointsOrig = dataDim.numPointsOrig
    
    parent = self.parent
    while parent and not hasattr(parent, 'currentPeaks'):
      parent = parent.parent
    if parent:
      currentPeaks = parent.currentPeaks
    else:
      currentPeaks = []
    selColor = hexToRgb(inverseRgb(getAnalysisProfile(view.root).bgColor))
    
    for winPeakList in view.windowPeakLists:
      if self.isWinPeakListDrawn(winPeakList):
        analysisPeakList = winPeakList.analysisPeakList
        peakList = analysisPeakList.peakList
        rgb1 = hexToRgb(analysisPeakList.symbolColor)
        rgb2 = hexToRgb(analysisPeakList.textColor)
        
        for peak in peakList.peaks:
          if peak in currentPeaks:
            handler.setColor(selColor)
            
          else:
            handler.setColor(rgb1)

          peakDim = peak.findFirstPeakDim(dim=dim)
          position = peakDim.position + peakDim.numAliasing*numPointsOrig
          position = convertPosition(peakDim.position, dataDimRef, toUnit=unit)
          if position > r1 and position < r0:
            if winPeakList.isAnnotationDrawn:
              text = getPeakAnnotation(peak, noPeakDimAnnotationChar='-')
              # draw text vertically
              n = len(text)
              r = size / thickness  # approximate height of one character
              s = n * r
            else:
              s = 0

            if winPeakList.isSymbolDrawn:
              handler.drawLine(position, s, position, 1)
            
            if winPeakList.isAnnotationDrawn:
              handler.setColor(rgb2)
              for i in range(n):
                handler.drawText(text[i], position, s-i*r, 0.5, 1.0)


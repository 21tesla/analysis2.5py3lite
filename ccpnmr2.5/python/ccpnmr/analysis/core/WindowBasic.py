
LICENSE = """
======================COPYRIGHT/LICENSE START==========================

WindowBasic.py: Part of the CcpNmr Analysis program

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

from ccp.api.nmr import Nmr

from ccpnmr.analysis.core import ExperimentBasic
from ccpnmr.analysis.core.PeakBasic import getPeakDimPosition, uniteResonancePeakPositions
from ccpnmr.analysis.core import Util

try:
  from memops.gui.MessageReporter import showWarning
except ImportError:
  from memops.universal.MessageReporter import showWarning

from memops.universal.Region1D import Region1D

VALUE_AXIS_NAME = 'value'

def getWindowPaneName(windowPane):
  """
  Get a name which identifies a window pane and
  its parent spectrum window
  
  .. describe:: Input
  
  Analysis.SpectrumWindowPane
  
  .. describe: Output
  
  Word
  """
  
  window = windowPane.spectrumWindow
  if len(window.spectrumWindowPanes) > 1:
    return '%s:%d' % (window.name, windowPane.serial)
  else:
    return window.name

def resetPeakLabelPositions(peaks, windowPane):
  """
  Reset the peak annotation labels in the x-y plane
  of a window for a given list of peaks.
  
  .. describe:: Input
  
  List of nmr.Nmr,peaks, Analysis.SpectrumWindow
  
  .. describe: Output
  
  None
  """

  setPeakOffset = Util.setPeakTextOffset
  setPeakDimOffset = Util.setPeakDimTextOffset
  
  for peak in peaks:
    view = getSpectrumWindowView(windowPane, peak.peakList.dataSource)
    if not view:
      continue
      
    setPeakOffset( peak, None )
    for label in ('x', 'y'):
      dataDim = view.findFirstAxisMapping(label=label).analysisDataDim.dataDim
      peakDim = peak.findFirstPeakDim(dim=dataDim.dim)
      setPeakDimOffset( peakDim, None )

  
def refinePeakLabelPositions(peaks, windowPane, analysis):
  """
  Automatically tidy the peak annotation labels in the x-y plane
  of a window for a given list of peaks.
  
  .. describe:: Input
  
  List of nmr.Nmr,peaks, Analysis.SpectrumWindow, Analysis GUI
  
  .. describe: Output
  
  None
  """
  from math import acos, sqrt, cos, sin
  
  offsX = []
  offsY = []
  posnX = []
  posnY = []
  dimsX = []
  dimsY = []
  facsX = []
  facsY = []
  
  
  
  for peak in peaks:
    view = getSpectrumWindowView(windowPane, peak.peakList.dataSource)
    if not view:
      continue
 
    if len(peak.peakDims) < 2:
      continue
 
    ppp  = []
    dims = []
    for label in ('x', 'y'):
      dataDim = view.findFirstAxisMapping(label=label).analysisDataDim.dataDim
      ppp.append( Util.calcPointsPerPixel(view, label) )
      dims.append( dataDim.dim )
      
    peakDimX = peak.findFirstPeakDim(dim=dims[0])
    peakDimY = peak.findFirstPeakDim(dim=dims[1])
    
    eX = ppp[0] * 20 # * pixesls -> points
    eY = ppp[1] * 20
        
    offsX.append(-eX) # Arbitrary start point
    offsY.append(-eY)
    
    posnX.append(peakDimX.position)
    posnY.append(peakDimY.position)
  
    dimsX.append(peakDimX)
    dimsY.append(peakDimY)

    facsX.append(eX)
    facsY.append(eY)
  
  N    = len(offsX)
  pi   = 3.14159265358979
  step = pi/20.0
  Etot = 1000.0
  Epre = 0.0
  
  setOff = Util.setPeakDimTextOffset
    
  c = 0
  while abs(Etot-Epre) > 0.0:
    
    if c > 100:
      break
    
    Epre = Etot
    Etot = 0.0   
    for i in range(N):
      
      oX = offsX[i]
      oY = offsY[i]
      
      fX = facsX[i]
      fY = facsY[i]
                 
      x0 = oX/fX
      y0 = oY/fY
      
      p0X = posnX[i]
      p0Y = posnY[i]
      
      setOff( dimsX[i], oX )
      setOff( dimsY[i], oY )

      Epk = None
      xB  = None
      yB  = None
 
      r = sqrt((x0*x0) + (y0*y0))
  
      a0 = acos(y0/r)
      
      if x0 < 0:
        a0 = -a0
      
      a2 = a0 - step
      a3 = a0 + step
 
      x1 = -x0
      y1 = -y0

      x2 = -x0
      y2 =  y0
      
      x3 =  x0
      y3 = -y0

      x4 = r * sin(a2) 
      y4 = r * cos(a2) 
      
      x5 = r * sin(a3)
      y5 = r * cos(a3)
   
      for x,y  in ((x0,y0),(x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5)):
 
        pX = p0X + (x * fX)
        pY = p0Y + (y * fY)
        
        E = 0.0 
        for j in range(N):
          if i == j:
            continue
 
          dx = fX * (pX-posnX[j]-offsX[j])
          dy = fY * (pY-posnY[j]-offsY[j])
          r2 = (dx * dx) + (dy * dy)
          E += 1.0/r2

          dx = fX * (pX-posnX[j])
          dy = fY * (pY-posnY[j])
          r2 = (dx * dx) + (dy * dy)
          E += 1.0/r2
        
        if (Epk is None) or (E <= Epk):
          Epk = E
          xB  = x
          yB  = y
          
      Etot += Epk
      offsX[i] = xB  * fX
      offsY[i] = yB  * fY     
    
    c += 1
    analysis.update_idletasks()

# convenience function
def getWindowFrameCanvas(argServer, canvas = None):
  """
  Get the popup, and canvas if None passed in, for this argServer.
  
  .. describe:: Input
  
  ArgumentServer, WindowCanvas
  
  .. describe: Output
  
  (WindowFrame, WindowCanvas)
  """

  if not canvas:
    canvas = argServer.getCurrentCanvas()

  try:
    windowFrame = canvas.parent.parent # slightly dangerous code, isolated here
  except:
    windowFrame = None

  return (windowFrame, canvas)

# convenience function
def getWindowFrameCanvasLocation(argServer, event = None):
  """
  Get the popup, canvas and event location, for this argServer.
  
  .. describe:: Input
  
  ArgumentServer, 'e'vent
  
  .. describe: Output
  
  (WindowFrame, WindowCanvas)
  """

  if event:
    canvas = event.widget
  else:
    canvas = None

  (windowFrame, canvas) = getWindowFrameCanvas(argServer, canvas)

  if not event:
    event = argServer.getCurrentEvent()
    if (event.widget != canvas):
      event = None

  if event:
    # x = event.x
    # y = event.y

    # NOTE:ED - replace the above event.<pos> because doesn't work on MacOS (OSbug)
    px, py = event.widget.winfo_pointerxy()
    rx, ry = (event.widget.winfo_rootx(), event.widget.winfo_rooty())
    x, y = (px - rx, py - ry)
  else:
    x = y = None

  return (windowFrame, canvas, x, y)


def getSpectrumViews(spectrum):
  """
  Get a list of all the window views of a spectrum.
  
  .. describe:: Input
  
  Nmr.DataSource
  
  .. describe: Output
  
  List of Analysis.SpectrumWindowViews
  """

  analysisSpectrum = Util.getAnalysisSpectrum(spectrum)

  return analysisSpectrum.sortedSpectrumWindowViews()


def isSpectrumInWindowPane(spectrumWindowPane, spectrum):
  """
  Determine if a spectrum is displayable in a given window.
  
  .. describe:: Input
  
  Analysis.SpectrumWindowPane, Nmr.DataSource
  
  .. describe: Output
  
  Boolean
  """
  analysisSpectrum = Util.getAnalysisSpectrum(spectrum)

  view = analysisSpectrum.findFirstSpectrumWindowView(spectrumWindowPane=spectrumWindowPane)

  return (view is not None)

def togglePeakList(window, peakList):
  """
  Toggle whether the picked peaks of a peak list are
  displayed in a spectrum window.
  
  .. describe:: Input
  
  Analysis.SpectrumWindow, Nmr.PeakList
  
  .. describe: Output
  
  None
  """
    
  winPeakList = None
  
  for windowPane in window.spectrumWindowPanes:
    view = getSpectrumWindowView(windowPane, peakList.dataSource)
    
    if view and view.isPosVisible:
      analysisPeakList = peakList.analysisPeakList
      winPeakList = view.findFirstWindowPeakList(analysisPeakList=analysisPeakList)
 
      if winPeakList:
        state = not winPeakList.isSymbolDrawn
        setWindowPeakListState(winPeakList, state)


def setWindowPeakListState(winPeakList, state):
  """
  Set the view state for a peak list in a window.
  This will set the drawing state of the annotation
  and the crosspeak symbol.
  
  .. describe:: Input
  
  Analysis.WindowPeakList, Boolean
  
  .. describe: Output
  
  None
  """
  
  winPeakList.state = state
  winPeakList.isSymbolDrawn = state
  winPeakList.isAnnotationDrawn = state


def toggleSpectrum(window, shortcut=None, spectrum=None, isGlobal=False):
  """
  Toggle whether a spectrum is displayed in a spectrum window.
  
  .. describe:: Input
  
  Analysis.SpectrumWindow, Int, Nmr.DataSource, Boolean
  
  .. describe: Output
  
  None
  """

  assert shortcut or spectrum
  # TJS: Currently spectrum overides shortcut

  # TBD: figure out how to deal with shift key (for negative contours)
  # event.keysym in this case is what appears above key (e.g. exclam)
  # as is event.char
  # but event.keycode is same as without shift key (e.g. 10 for '1', etc.)
  # can one rely that keycode is always 10 for '1', etc.?

  if spectrum is None:
    spectrum = Util.getShortcutSpectrum(window.root, shortcut)
  
  if spectrum:
    analysisSpectrum = Util.getAnalysisSpectrum(spectrum)
    refView = None
    
    for windowPane in window.spectrumWindowPanes:
      getView = windowPane.findFirstSpectrumWindowView
      refView = getView(analysisSpectrum=analysisSpectrum)
      if refView:
        break
 
    if not refView:
      return
 
    views = []
    if isGlobal:
      views = analysisSpectrum.spectrumWindowViews
  
    else:
      for windowPane in window.spectrumWindowPanes:
        getView = windowPane.findFirstSpectrumWindowView
        view = getView(analysisSpectrum=analysisSpectrum)
        
        if view:
          views.append(view)
 
    if windowPaneHasValueAxis(refView.spectrumWindowPane):
      state = not refView.isSliceVisible
    else:
      state = not refView.isPosVisible
 
    for view in views:
      if windowPaneHasValueAxis(view.spectrumWindowPane):
        view.isSliceVisible = state
      else:
        view.isPosVisible = view.isNegVisible = view.isSliceVisible = state
 
      for winPeakList in view.windowPeakLists:
        if state and hasattr(winPeakList, 'state'):
          winPeakList.isSymbolDrawn = winPeakList.isAnnotationDrawn = winPeakList.state
        else:
          winPeakList.isSymbolDrawn = winPeakList.isAnnotationDrawn = state
 

  return spectrum

def toggleWindowSpectra(window, isGlobal=False):
  """
  Turn all spectra off in a spectrum window (all windows if isGlobal is True)
  if any on (in that window) or all spectra on if all of them are off.
  
  .. describe:: Input
  
  Analysis.SpectrumWindow, Boolean
  
  .. describe: Output
  
  None
  """

  views = []
  for windowPane in window.spectrumWindowPanes:
    views2 = [view for view in windowPane.spectrumWindowViews if view.isInToolbar]
    views.extend(views2)
  
  state = True
  for view in views:
    if windowPaneHasValueAxis(view.spectrumWindowPane):
      if view.isSliceVisible:
        state = False
        break
    else:
      if view.isPosVisible or view.isNegVisible:
        state = False
        break

  if isGlobal:
    for window2 in window.parent.spectrumWindows:
      if window2 is not window:
        for windowPane in window2.spectrumWindowPanes:
          views2 = [view for view in windowPane.spectrumWindowViews if view.isInToolbar]
          views.extend(views2)
  
  for view in views:
    if windowPaneHasValueAxis(view.spectrumWindowPane):
      view.isSliceVisible = state
    else:
      view.isPosVisible = view.isNegVisible = view.isSliceVisible = state
 
    for winPeakList in view.windowPeakLists:
      winPeakList.isSymbolDrawn = winPeakList.isAnnotationDrawn = state
 
def turnWindowSpectraOn(window, isGlobal=False):
  """
  Turn all spectra on in a spectrum window (all windows if isGlobal is True).
  
  .. describe:: Input
  
  Analysis.SpectrumWindow, Boolean
  
  .. describe: Output
  
  None
  """

  views = []
  for windowPane in window.spectrumWindowPanes:
    views2 = [view for view in windowPane.spectrumWindowViews if view.isInToolbar]
    views.extend(views2)
  
  if isGlobal:
    for window2 in window.parent.spectrumWindows:
      if window2 is not window:
        for windowPane in window2.spectrumWindowPanes:
          views2 = [view for view in windowPane.spectrumWindowViews if view.isInToolbar]
          views.extend(views2)
  
  for view in views:
    if windowPaneHasValueAxis(view.spectrumWindowPane):
      view.isSliceVisible = True
    else:
      view.isPosVisible = view.isNegVisible = view.isSliceVisible = True
 
    for winPeakList in view.windowPeakLists:
      winPeakList.isSymbolDrawn = winPeakList.isAnnotationDrawn = True
 
def turnWindowSpectraOff(window, isGlobal=False):
  """
  Turn all spectra off in a spectrum window (all windows if isGlobal is True).
  
  .. describe:: Input
  
  Analysis.SpectrumWindow, Boolean
  
  .. describe: Output
  
  None
  """

  views = []
  for windowPane in window.spectrumWindowPanes:
    views2 = [view for view in windowPane.spectrumWindowViews if view.isInToolbar]
    views.extend(views2)
  
  if isGlobal:
    for window2 in window.parent.spectrumWindows:
      if window2 is not window:
        for windowPane in window2.spectrumWindowPanes:
          views2 = [view for view in windowPane.spectrumWindowViews if view.isInToolbar]
          views.extend(views2)
  
  for view in views:
    if windowPaneHasValueAxis(view.spectrumWindowPane):
      view.isSliceVisible = False
    else:
      view.isPosVisible = view.isNegVisible = view.isSliceVisible = False
 
    for winPeakList in view.windowPeakLists:
      winPeakList.isSymbolDrawn = winPeakList.isAnnotationDrawn = False
 
def getSpectrumWindowView(windowPane, spectrum):
  """
  Get the SpectrumWindowView object for a given spectrum in a
  given window pane.
  
  .. describe:: Input
  
  Analysis.SpectrumWindowPane, Nmr.DataSource
  
  .. describe: Output
  
  Analysis.SpectrumWindowView
  """

  if not windowPane or not spectrum:
    return None

  analysisSpectrum = Util.getAnalysisSpectrum(spectrum)
  
  view = windowPane.findFirstSpectrumWindowView(analysisSpectrum=analysisSpectrum)
   
  return view
    
def getDefaultAspectRatio(axisTypes):

  weight = {'1H': 1.0,
            '13C': 10.0,
            'DQ13C': 10.0,
            '15N': 8.0}

  if len(axisTypes) >= 2:
    return weight.get(axisTypes[1].name, 1.0) / weight.get(axisTypes[0].name, 1.0)
  else:
    return 1.0

def defaultWindowName(project):
  """
  Determine a default window name.
  
  .. describe:: Input
  
  Project
  
  .. describe: Output
  
  Window name
  """
  
  names = [ window.name for window in project.currentAnalysisProject.spectrumWindows ]
  n = len(names) + 1
  name = 'window%d' % (n,)
  
  while name in names:
    n += 1
    name = 'window%d' % (n,)

  return name

def copySpectrumWindowRegions(oldWindow, newWindow):
  """
  Copy regions from one window to another (pane and region numbers must match or else)
  
  .. describe:: Input
  
  Analysis.SpectrumWindow, Analysis.SpectrumWindow
  
  .. describe: Output
  
  None
  """

  newPanes = newWindow.sortedSpectrumWindowPanes()
  oldPanes = oldWindow.sortedSpectrumWindowPanes()

  for i, oldPane in enumerate(oldPanes):
    newPane = newPanes[i]
  
    oldAxisPanels = oldPane.sortedAxisPanels()
    newAxisPanels = newPane.sortedAxisPanels()

    for j, oldAxisPanel in enumerate(oldAxisPanels):
      newAxisPanel = newAxisPanels[j]

      oldAxisRegions = oldAxisPanel.sortedAxisRegions()
      newAxisRegions = newAxisPanel.sortedAxisRegions()

      for k, oldAxisRegion in enumerate(oldAxisRegions):
        newAxisRegion = newAxisRegions[k]
        for attr in ('isActive', 'overrideRegion', 'region', 'size'):
          setattr(newAxisRegion, attr, getattr(oldAxisRegion, attr))

# note, should probably also call copySpectrumWindowViewProperties() after an idle
# and also copySpectrumWindowRegions()
def cloneSpectrumWindow(window, name = None):
  """
  Clone an existing spectrum window to make a new one
  
  .. describe:: Input
  
  Analysis.SpectrumWindow
  
  .. describe: Output
  
  Analysis.SpectrumWindow
  """

  project = window.root
  if not name:
    name = defaultWindowName(project)
  
  axisTypes = []
  regions = []
  sizes = []
  for windowPane in window.spectrumWindowPanes:
    axisPanels = windowPane.sortedAxisPanels()
    rgs = []
    szs = []
    for axisPanel in axisPanels:
      axisRegions = axisPanel.sortedAxisRegions()
      rgs.append([axisRegion.region for axisRegion in axisRegions])
      if axisPanel.label in ('x', 'y'):
        szs.append([axisRegion.size for axisRegion in axisRegions])
    regions.append(rgs)
    sizes.append(szs)
    axisTypes.append( [ axisPanel.axisType for axisPanel in axisPanels ] )
 
  newWindow = createSpectrumWindow(project, name, axisTypes,
                                   sizes=sizes, regions=regions)

  for attr in ('isXTickShown', 'isYTickShown', 'useOverrideRegion',
               'isCanvasLabelShown', 'isCanvasMidpointShown', 'stripAxis',
               'useMultiplePeakLists', 'isXSliceDrawn', 'isYSliceDrawn'):
    setattr(newWindow, attr, getattr(window, attr))

  newPanes = newWindow.sortedSpectrumWindowPanes()
  oldPanes = window.sortedSpectrumWindowPanes()

  for j, windowPane in enumerate(newPanes):
    oldPane = oldPanes[j]

    for attr in ('aspectRatio', 'sliceRange'):
      setattr(windowPane, attr, getattr(oldPane, attr))
  
    for axisPanel in windowPane.axisPanels:
      windowPanel = oldPane.findFirstAxisPanel(label=axisPanel.label)
      for attr in ('isVisible', 'thickness', 'axisUnit', 'panelType'):
        setattr(axisPanel, attr, getattr(windowPanel, attr))

    for slicePanel in windowPane.slicePanels:
      panel = oldPane.findFirstSlicePanel(label=slicePanel.label)
      for attr in ('isVisible', 'thickness'):
        setattr(slicePanel, attr, getattr(panel, attr))

    newAxisPanels = windowPane.sortedAxisPanels()
    oldAxisPanels = oldPane.sortedAxisPanels()
    for n in range(len(oldAxisPanels)):
      newAxisPanel = newAxisPanels[n]
      newAxisRegions = newAxisPanel.sortedAxisRegions()
      axisPanel = oldAxisPanels[n]
      axisRegions = axisPanel.sortedAxisRegions()
      for i in range(len(axisRegions)):
        axisRegion = axisRegions[i]
        if n < 2:
          """ below gets trampled by later calls so put in copySpectrumWindowRegions
          if i == 0:
            newAxisRegion = newAxisRegions[0]
            for attr in ('region', 'size', 'isActive', 'overrideRegion'):
              setattr(newAxisRegion, attr, getattr(axisRegion, attr))
          else:
"""
          if i > 0:
            newAxisPanel.newAxisRegion(region=axisRegion.region, size=axisRegion.size,
                isActive=axisRegion.isActive, overrideRegion=axisRegion.overrideRegion)
        """ below gets trampled by later calls so put in copySpectrumWindowRegions
        else:
          # the orthogonal axisRegions created by notifier on x/y
          # axisRegion creation in Analysis.py
          newAxisRegion = newAxisRegions[i]
          for attr in ('isActive', 'overrideRegion', 'region', 'size'):
            setattr(newAxisRegion, attr, getattr(axisRegion, attr))
"""

  return newWindow

def copySpectrumWindowViewProperties(oldWindow, newWindow):
  """
  Copy view visibility properties from one window to another
  
  .. describe:: Input
  
  Analysis.SpectrumWindow, Analysis.SpectrumWindow
  
  .. describe: Output
  
  None
  """

  newPanes = newWindow.sortedSpectrumWindowPanes()
  oldPanes = oldWindow.sortedSpectrumWindowPanes()
 
  for j, windowPane in enumerate(newPanes):
  
    # Assumes panes are in same order
    if j < len(oldPanes):
      oldPane = oldPanes[j]
    else:
      break
        
    findOldView = oldPane.findFirstSpectrumWindowView

    for newView in windowPane.spectrumWindowViews:
 
      oldView = findOldView(analysisSpectrum=newView.analysisSpectrum)
      if oldView: # should be the case
        isPosVis = newView.isPosVisible = oldView.isPosVisible
        isNegVis = newView.isNegVisible = oldView.isNegVisible
        newView.isSliceVisible = oldView.isSliceVisible
        newView.isInToolbar = oldView.isInToolbar

        for winPeakList in newView.windowPeakLists:
          oldWinPeakList = oldView.findFirstWindowPeakList(analysisPeakList=winPeakList.analysisPeakList)
          if oldWinPeakList:
            winPeakList.isSymbolDrawn = oldWinPeakList.isSymbolDrawn
            winPeakList.isAnnotationDrawn = oldWinPeakList.isAnnotationDrawn
          elif not (isPosVis or isNegVis):
            winPeakList.isSymbolDrawn = False
            winPeakList.isAnnotationDrawn = False

      newAxisMappings = newView.sortedAxisMappings()
      oldAxisMappings = oldView.sortedAxisMappings()
      needAxisMappingChange = False
      
      for i in range(len(oldAxisMappings)):
        if newAxisMappings[i].analysisDataDim is not \
           oldAxisMappings[i].analysisDataDim:
          needAxisMappingChange = True
          break

      if needAxisMappingChange:
        makeMapping = newView.newAxisMapping
        
        for axisMapping in newAxisMappings:
          axisMapping.delete()
          
        for oldAxisMapping in oldAxisMappings:
          makeMapping(label=oldAxisMapping.label,
                      analysisDataDim=oldAxisMapping.analysisDataDim)

def createSpectrumWindow(project, name, axisTypes, spectrum=None,
                        sizes=None, ncols=1, nrows=1, regions=None, **kw):
  """
  Make a new spectrum window
  
  .. describe:: Input
  
  Project, String,List of List of Analysis.AxisTypes,
  Nmr.DataSource, RgbColor, List of 2-List of List of Floats, Int, Int,
  List of List of List of Floats
  
  .. describe: Output
  
  Analysis.SpectrumWindow
  """

  analysisProject = Util.getAnalysisProject(project)

  # TBD: this is a bit of a hack
  # problem is that we only want frequency axis types
  # in x, y or value axis for y
  
  overideRatio = kw.get('aspectRatio')
  if overideRatio:
    # Can't be set on an Analysis spectrum any longer
    del kw['aspectRatio']
  
  window = analysisProject.newSpectrumWindow(name=name, **kw)
  
  aspectRatios = []
  for p, paneAxisTypes in enumerate(axisTypes):
    xyInds = []
    for n, axisType in enumerate(paneAxisTypes):
      if axisType.isotopeCodes or (axisType.name == VALUE_AXIS_NAME):
        xyInds.append(n)
        if len(xyInds) == 2:
          break

    # add arbitrary axisType for < 2D spectrum case
    n = len(xyInds)
    if n < 2:
      axisType = analysisProject.findFirstAxisType(name='H')
      # H is arbitary, should always exist
      m = len(axisTypes)
      for i in range(n, 2):
        xyInds.append(m+i)
        paneAxisTypes.append(axisType)

    xyAxisTypes = [ paneAxisTypes[xyInds[n]] for n in range(2) ]
    if overideRatio:
      aspectRatio = overideRatio
    
    else:
      if VALUE_AXIS_NAME in [at.name for at in xyAxisTypes]:
        aspectRatio = 1
      else:
        aspectRatio = getDefaultAspectRatio(xyAxisTypes)

    if p < 26:
      name = chr(ord('A') + p)
    else:
      name = chr(ord('A') - 1 + p // 26) + chr(ord('A') + p % 26)
    
    windowPane = window.newSpectrumWindowPane(aspectRatio=aspectRatio,
                                        name=name)
 
    sortedAxisPanels = windowPane.sortedAxisPanels()
 
    if regions:
      hadRegions = True
      paneRegions = regions[p]
    else:
      windowPane.hasDefaultRegions = True # so that later can set region sensibly
      hadRegions = False
      paneRegions = 2 * [0]
      for n in range(2):
        axisType = xyAxisTypes[n]
        # for 1D windows try and select sensible region
        if n == 1 and axisType.name == 'value' and \
            spectrum and hasattr(spectrum, 'block_file') and \
            spectrum.block_file:
          (minValue, maxValue) = ExperimentBasic.getMinMaxValues(spectrum)
          delta = maxValue - minValue
          scale = 0.2
          minValue -= scale * delta
          maxValue += scale * delta
          paneRegions[n] = (minValue, maxValue)
        else:
          # TBD: because of introduction of xyInds,
          # is below still ok?
          panelType = Util.findPanelType(project, axisType,
                                         sortedAxisPanels[:n])
          paneRegions[n] = panelType.axisType.region

    if sizes:
      size = sizes[p]
      
    else:
      size = (int(500/ncols), int(500/nrows))

    for n in range(2):
      m = xyInds[n]
      axisType = paneAxisTypes[m]
      axisPanel = sortedAxisPanels[n]
      axisPanel.panelType = Util.findPanelType(project, axisType,
                                               sortedAxisPanels[:n])
      axisPanel.axisUnit = axisType.findFirstAxisUnit(unit='ppm')
      axisPanel.isVisible = True
      axisPanel.firstMapping = True
 
      # at this point there is only one axisRegion, others get added below
      axisRegion = axisPanel.findFirstAxisRegion()
      if hadRegions:
        region = paneRegions[m][0]
      else:
        region = paneRegions[n]
      axisRegion.region = region

      if sizes:
        sz = size[n][0]
      else:
        sz = size[n]
      axisRegion.size = sz
 
    # below was when only had slice panel to show assignments, now also do on canvas
    # so there is not so much point making this slicePanel visible by default
    ##axisPanel = sortedAxisPanels[1]
    ##if axisPanel.axisType.name == VALUE_AXIS_NAME:
    ##  axisPanel = sortedAxisPanels[0]
    ##  slicePanel = axisPanel.slicePanel
    ##  slicePanel.isVisible = True
    ##  slicePanel.thickness = 100

    m = 0
    for n in range(len(paneAxisTypes)):
      axisType = paneAxisTypes[n]
 
      # this needs to be before continue so that m is incremented correctly
      if axisType.name != VALUE_AXIS_NAME:
        m += 1

      if n in xyInds:
        continue
 
      if hadRegions:
        region = paneRegions[n][0]
 
      else:
        region = None
 
        if spectrum is not None:
          dataDim = spectrum.findFirstDataDim(dim=m)
 
          if dataDim.className == 'FreqDataDim':
            dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
            axisUnit = axisType.findFirstAxisUnit(unit='ppm')
            p = 0.5 * dataDim.numPoints
            r0 = Util.convertPosition(p-0.5, dataDimRef, toUnit=axisUnit.unit)
            r1 = Util.convertPosition(p+0.5, dataDimRef, toUnit=axisUnit.unit)
 
            if axisUnit.isBackwards:
              region = (r1, r0)
            else:
              region = (r0, r1)
          
          else:
            region = (0.5, 1.5) # show first plane only be default
        elif axisType.name == 'sampled':
          region = (0.5, 1.5) # show first plane only be default
 
      Util.createAxisPanel(windowPane, axisType, region=region)

    nstrips = (ncols-1, nrows-1)
    for n in range(2):
      axisPanel  = sortedAxisPanels[n]
      if hadRegions:
        rgns = paneRegions[n][1:]
      else:
        axisRegion = axisPanel.findFirstAxisRegion()
        rgns = nstrips[n] * (axisRegion.region,)
      if sizes:
        szs = size[n][1:]
      else:
        szs = nstrips[n] * (size[n],)
      for i in range(nstrips[n]):
        axisPanel.newAxisRegion(region=rgns[i], size=szs[i])

  if spectrum:
    Util.defaultContourLevels(spectrum)

  return window

def isActiveWindow(window):
  """
  Determine if a window is active within its project
  
  .. describe:: Input
  
  Analysis.SpectrumWindow
  
  .. describe: Output
  
  Boolean
  """

  group = window.analysisProject.activeWindowGroup
  if group and window in group.spectrumWindows:
    return True
    
  else:
    return False

def getActiveWindows(project):
  """
  Get the currently active windows in a project
  
  .. describe:: Input
  
  Implementation.Project
  
  .. describe: Output
  
  List of Analysis.SpectrumWindows
  """
  
  analysisProject = project.currentAnalysisProject

  #return [ w for w in analysisProject.sortedSpectrumWindows() if isActiveWindow(w) ]
  return sorted([ w for w in analysisProject.spectrumWindows if isActiveWindow(w) ], key=lambda window: window.name)

def getEquivalentWindows(queryWindow):
  """
  Get a list of windows with the same axis types as the query window
  
  .. describe:: Input
  
  Analysis.SpectrumWindow
  
  .. describe: Output
  
  List of Analysis.SpectrumWindows
  """

  analysisProject = queryWindow.topObject
  queryTypes = getWindowAxisTypes(queryWindow)
  windows    = []
  for window in analysisProject.sortedSpectrumWindows():
    if getWindowAxisTypes(window) == queryTypes:
      windows.append(window)
    
  return windows

def getWindowAxisTypes(window):
  """
  Get a list of names representing the axis types of a window.
  
  .. describe:: Input
  
  Analysis.SpectrumWindow
  
  .. describe: Output
  
  List of Strings (Analysis.AxisType.name)
  """
 
  types = []
  
  for windowPane in window.sortedSpectrumWindowPanes():
    paneTypes = []
    
    for axisPanel in windowPane.sortedAxisPanels():

      if axisPanel.axisType:
        paneTypes.append(axisPanel.axisType.name)
 
      else:
        paneTypes.append(None)
    
    types.append(paneTypes)
   
  return types


def windowPaneHasValueAxis(windowPane):
  """
  Determine if a window has a value axis
  
  .. describe:: Input
  
  Analysis.SpectrumWindow
  
  .. describe: Output
  
  Boolean
  """

  # TBD: more general?
  try:
    if windowPane.findFirstAxisPanel(label='y').axisType.name == VALUE_AXIS_NAME:
      return True
        
    else:
      return False
      
  except:
    return False

def getViewMinMaxValues(view):
  """
  Get the minimum and maximum spectrum data values in this view.
  Returns (None, None) if there is no data file for this spectrum.
  
  .. describe:: Input
  
  SpectrumWindowView
  
  .. describe: Output
  
  (minValue, maxValue)
  """

  import math

  minValue = maxValue = None

  analysisSpectrum = view.analysisSpectrum
  if analysisSpectrum:
    spectrum = analysisSpectrum.dataSource
    if spectrum and hasattr(spectrum, 'block_file') and spectrum.block_file:
      windowPane = view.spectrumWindowPane
      block_file = spectrum.block_file
      ndim = spectrum.numDim
      first = ndim * [0]
      last = ndim * [0]
      for axisMapping in view.axisMappings:
        dataDim = axisMapping.analysisDataDim.dataDim
        """
        dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
        axisPanel = windowPane.findFirstAxisPanel(label=axisMapping.label)
        axisRegion = axisPanel.findFirstAxisRegion()
        axisUnit = axisPanel.axisUnit
        region = axisRegion.region
        unit = axisUnit.unit
        (r0, r1) = [Util.convertPosition(r, dataDimRef, fromUnit=unit) for r in region]
        if axisUnit.isBackwards:
          (r0, r1) = (r1, r0)
"""
        dim = dataDim.dim
        numPoints = dataDim.numPoints
        """
        numPointsOrig = dataDim.numPointsOrig
        pointOffset = dataDim.pointOffset
        r0 = int(math.floor(r0))
        r1 = int(math.ceil(r1))
        if (r0 // numPointsOrig) != (r1 // numPointsOrig): # complicated case, so just use everything (so to do it properly would need to consider multiple (first, last) regions
          r0 = 0
          r1 = numPoints
        else:
          r0 = max(r0 % numPointsOrig - pointOffset, 0)
          r1 = min(r1 % numPointsOrig - pointOffset, numPoints)
"""
        c = numPoints // 2
        n = max(1, numPoints // 10)
        r0 = c
        r1 = c + n
        
        first[dim-1] = r0
        last[dim-1] = r1

      v = block_file.maxValue(first, last)
      if maxValue is None or v > maxValue:
        maxValue = v
      v = block_file.minValue(first, last)
      if minValue is None or v < minValue:
        minValue = v

  return (minValue, maxValue)

def updateWindowPaneSliceRange(windowPane):
  """
  Update the windowPane sliceRange based on the window spectrum values 

  .. describe:: Input
  
  SpectrumWindowPane
  
  .. describe: Output
  
  None
  """

  minValue = maxValue = None
  for view in windowPane.spectrumWindowViews:
    try:
      (minV, maxV) = getViewMinMaxValues(view)
    except:
      continue
    if maxV is not None and (maxValue is None or maxValue < maxV):
      maxValue = maxV
    if minV is not None and (minValue is None or minValue > minV):
      minValue = minV

  if maxValue is not None and minValue is not None:
    d = 0.1 * (maxValue - minValue)
    maxValue += d
    minValue -= d
    windowPane.sliceRange = (minValue, maxValue)

def updateWindowPaneSliceRangeByView(view):
  """
  Update the view.windowPane sliceRange based on the spectrum values 

  .. describe:: Input
  
  SpectrumWindowView
  
  .. describe: Output
  
  None
  """

  (minValue, maxValue) = getViewMinMaxValues(view)
  if maxValue is not None and minValue is not None:
    d = 0.1 * (maxValue - minValue)
    maxValue += d
    minValue -= d
    view.spectrumWindowPane.sliceRange = (minValue, maxValue)


def autoPeakLabelPos(argServer):
  """
  Automatically tidy the annotation labels of selected peaks.
  
  .. describe:: Input
  
  ArgumentServer
  
  .. describe: Output
  
  None
  """

  peaks  = argServer.getCurrentPeaks()
  windowPane = argServer.getCurrentWindowPane()
  
  if windowPane and peaks:
    refinePeakLabelPositions(peaks, windowPane, argServer.parent)

def resetPeakLabelsPos(argServer):
  """
  Automatically tidy the annotation labels of selected peaks.
  
  .. describe:: Input
  
  ArgumentServer
  
  .. describe: Output
  
  None
  """

  peaks  = argServer.getCurrentPeaks()
  windowPane = argServer.getCurrentWindowPane()
  
  if windowPane and peaks:
    resetPeakLabelPositions(peaks, windowPane)

def zoomWindowOut(argServer = None, canvas = None, factor = 2.0):
  """
  Zoom out window by factor (default 2).
  
  .. describe:: Input
  
  ArgumentServer, WindowCanvas, factor
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas) = getWindowFrameCanvas(argServer, canvas)

  if windowFrame and canvas:
    windowFrame.zoom(canvas, factor)

def zoomWindowIn(argServer = None, canvas = None, factor = 2.0):
  """
  Zoom in window by factor (default 2).
  
  .. describe:: Input
  
  ArgumentServer, WindowCanvas, factor
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas) = getWindowFrameCanvas(argServer, canvas)

  if (windowFrame and canvas):
    windowFrame.zoom(canvas, 1.0/factor)

def zoomWindowToSpectra(argServer = None, canvas = None):
  """
  Zoom window to just contain the visible spectra regions.
  
  .. describe:: Input
  
  ArgumentServer, WindowCanvas
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas) = getWindowFrameCanvas(argServer, canvas)

  if windowFrame and canvas:
    windowFrame.zoomToSpectra(canvas)

# convenience function
def translationFraction(argServer, fraction):
  """
  Translate window by fraction, with up/down choice determined by pan mode
  
  .. describe:: Input
  
  ArgumentServer, fraction
  
  .. describe: Output
  
  None
  """

  panMode = argServer.getAnalysisProfile().panView
  if (panMode):
    t = -fraction
  else:
    t = fraction

  return t

def translateWindowUp(argServer=None, canvas=None, fraction=0.25):
  """
  Translate window up by fraction (default 0.25)
  
  .. describe:: Input
  
  ArgumentServer, WindowCanvas, fraction
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas) = getWindowFrameCanvas(argServer, canvas)

  if windowFrame and canvas:
    t = - translationFraction(argServer, fraction)
    windowFrame.translate(canvas, 0, t)

def translateWindowDown(argServer=None, canvas=None, fraction=0.25):
  """
  Translate window down by fraction (default 0.25)
  
  .. describe:: Input
  
  ArgumentServer, WindowCanvas, fraction
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas) = getWindowFrameCanvas(argServer, canvas)

  if windowFrame and canvas:
    t = translationFraction(argServer, fraction)
    windowFrame.translate(canvas, 0, t)

def translateWindowRight(argServer=None, canvas=None, fraction=0.25):
  """
  Translate window right by fraction (default 0.25)
  
  .. describe:: Input
  
  ArgumentServer, WindowCanvas, fraction
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas) = getWindowFrameCanvas(argServer, canvas)

  if windowFrame and canvas:
    t = - translationFraction(argServer, fraction)
    windowFrame.translate(canvas, t, 0)

def translateWindowLeft(argServer=None, canvas=None, fraction=0.25):
  """
  Translate window left by fraction (default 0.25)
  
  .. describe:: Input
  
  ArgumentServer, WindowCanvas, fraction
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas) = getWindowFrameCanvas(argServer, canvas)

  if windowFrame and canvas:
    t = translationFraction(argServer, fraction)
    windowFrame.translate(canvas, t, 0)

def zoomSliceRangeDown(argServer=None, windowPane=None, factor=2.0):
  """
  Zoom slice range by factor (default 2).
  
  .. describe:: Input
  
  ArgumentServer, Analysis.WindowPane, float
  
  .. describe: Output
  
  None
  """

  if not windowPane:
    windowPane = argServer.getCurrentWindowPane()

  if windowPane:
  
    if argServer:
      canvas = argServer.getCurrentCanvas()
    else:
      canvas = None

    if hasattr(windowPane, 'getWindowFrame'):
      windowFrame = windowPane.getWindowFrame()
      windowFrame.changeSliceRange(1.0/factor, canvas)

def zoomSliceRangeUp(argServer=None, windowPane=None, factor=2.0):
  """
  Zoom slice range by factor (default 2).
  
  .. describe:: Input
  
  ArgumentServer, Analysis.WindowPane, float
  
  .. describe: Output
  
  None
  """

  if not windowPane:
    windowPane = argServer.getCurrentWindowPane()

  if windowPane:
  
    if argServer:
      canvas = argServer.getCurrentCanvas()
    else:
      canvas = None

    if hasattr(windowPane, 'getWindowFrame'):
      windowFrame = windowPane.getWindowFrame()
      windowFrame.changeSliceRange(factor, canvas)

def deleteSelected(argServer=None, analysisPopup=None):
  """
  Delete selected objects.
  
  .. describe:: Input
  
  ArgumentServer, AnalysisPopup
  
  .. describe: Output
  
  None
  """

  if not analysisPopup:
    if argServer.inGui:
      analysisPopup = argServer.parent

  if analysisPopup:
    if argServer:
      originator = argServer.getCurrentWindowPopup()
    else:
      originator = None
    analysisPopup.queryDeleteSelected(originator)

def assignPeak(argServer=None, event=None):
  """
  Assign peak at location (if any nearby).
  
  .. describe:: Input
  
  ArgumentServer, Tk.Event
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas, x, y) = getWindowFrameCanvasLocation(argServer, event)

  if (windowFrame and canvas and x is not None and y is not None):
    windowFrame.assignPeakAtLocation(canvas, x, y)

def showSelectedPeaks(argServer=None, analysisPopup=None):
  """
  Show selected peaks.
  
  .. describe:: Input
  
  ArgumentServer, AnalysisPopup
  
  .. describe: Output
  
  None
  """

  if not analysisPopup:
    if argServer.inGui:
      analysisPopup = argServer.parent

  if analysisPopup:
    analysisPopup.viewSelectedPeaks()

def centerAxisRegion(axisRegion, center):
  """
  Move and axis region so that it is centred on the input value
  
  .. describe:: Input
  
  Analysis.AsisRegion, Float
  
  .. describe: Output
  
  None
  """

  (r0, r1) = axisRegion.region
  d = center - 0.5*(r1+r0)
  axisRegion.region = (d+r0, d+r1)

def centerWindow(argServer=None, event=None):
  """
  Center window at location.
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas, x, y) = getWindowFrameCanvasLocation(argServer, event)

  if (windowFrame and canvas and x is not None and y is not None):
    windowFrame.centerAtLocation(canvas, x, y)

def createMark(argServer=None, event=None):
  """
  Create mark at location.
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas, x, y) = getWindowFrameCanvasLocation(argServer, event)

  if (windowFrame and canvas and x is not None and y is not None):
    windowFrame.doCreateMark(canvas, x, y)

def createSidebandMarks(argServer=None, event=None):
  """
  Create mark at location and sideband copies.
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas, x, y) = getWindowFrameCanvasLocation(argServer, event)

  if (windowFrame and canvas and x is not None and y is not None):
    windowFrame.doCreateMark(canvas, x, y, doSidebands=True)

def createHorizontalRuler(argServer=None, event=None):
  """
  Create horizontal ruler at location.
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas, x, y) = getWindowFrameCanvasLocation(argServer, event)

  if (windowFrame and canvas and x is not None and y is not None):
    windowFrame.doCreateHorizontalRuler(canvas, x, y)

def createSidebandHorizontalRulers(argServer=None, event=None):
  """
  Create horizontal ruler at location and sideband copies.
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas, x, y) = getWindowFrameCanvasLocation(argServer, event)

  if (windowFrame and canvas and x is not None and y is not None):
    windowFrame.doCreateHorizontalRuler(canvas, x, y, doSidebands=True)

def createVerticalRuler(argServer=None, event=None):
  """
  Create vertical ruler at location.
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas, x, y) = getWindowFrameCanvasLocation(argServer, event)

  if (windowFrame and canvas and x is not None and y is not None):
    windowFrame.doCreateVerticalRuler(canvas, x, y)

def createSidebandVerticalRulers(argServer=None, event=None):
  """
  Create vertical ruler at location and sideband copies.
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas, x, y) = getWindowFrameCanvasLocation(argServer, event)

  if (windowFrame and canvas and x is not None and y is not None):
    windowFrame.doCreateVerticalRuler(canvas, x, y, doSidebands=True)

def movePeak(argServer=None, event=None):
  """
  Move peak to location.
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas, x, y) = getWindowFrameCanvasLocation(argServer, event)

  if (windowFrame and canvas and x is not None and y is not None):
    windowFrame.moveSelectedPeak(canvas, x, y)

def movePeakAnnotation(argServer=None, event=None):
  """
  Move peak to location.
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas, x, y) = getWindowFrameCanvasLocation(argServer, event)

  if (windowFrame and canvas and x is not None and y is not None):
    windowFrame.moveSelectedPeakAnnotation(canvas, x, y)

def snapPeaks(argServer=None, windowPane=None):
  """
  Move peak to location.
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  if argServer and not windowPane:
    windowPane = argServer.getCurrentWindowPane()

  if windowPane and hasattr(windowPane, 'getWindowFrame'):
    windowFrame = windowPane.getWindowFrame()
    windowFrame.snapSelectedPeaks()

def popupMenu(argServer=None, event=None):
  """
  Popup right mouse menu
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  if not event and argServer:
    event = argServer.getCurrentEvent()

  if not event:
    return

  (windowFrame, canvas) = getWindowFrameCanvas(argServer, event.widget)

  if windowFrame and canvas:
    windowFrame.popupMenu(event)

def unitePeakPositions(argServer=None, event=None, allowMultiplePeaks=False):
  """
  Move peak to current cursor location and coordinate positions
  of other peaks in same peak list assigned to same resonances.
  If allowMultiplePeaks then allow multiple peaks to be moved.
  
  .. describe:: Input
  
  ArgumentServer, Nmr.Peaks, Bool
  
  .. describe: Output
  
  None
  """

  if not event and argServer:
    event = argServer.getCurrentEvent()

  if not event:
    return

  peaks = argServer.getCurrentPeaks()
    
  if peaks:
    windowFrame, canvas = getWindowFrameCanvas(argServer, event.widget)
    if len(peaks) > 1 and not allowMultiplePeaks:
      msg = 'Only one peak should be selected to coordinate positions.'
      showWarning('Multiple peaks', msg, parent=windowFrame)
      return
      
    for peak in peaks:
      uniteResonancePeakPositions(peak)
    
def unitePeakPositionsMulti(argServer=None, event=None):
  """
  Move peak to current cursor location and coordinate positions
  of other peaks in same peak list assigned to same resonances.
  
  .. describe:: Input
  
  ArgumentServer, Nmr.Peaks
  
  .. describe: Output
  
  None
  """

  unitePeakPositions(argServer, event, allowMultiplePeaks=True)

def orthogScrollLeft(argServer=None, canvas=None):
  """
  Scroll to left in orthogonal (z1) dimension in window.
  
  .. describe:: Input
  
  ArgumentServer, WindowCanvas
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas) = getWindowFrameCanvas(argServer, canvas)

  if windowFrame and canvas:
    windowFrame.orthogScrollLeft(canvas)

def orthogScrollRight(argServer=None, canvas=None):
  """
  Scroll to right in orthogonal (z1) dimension in window.
  
  .. describe:: Input
  
  ArgumentServer, WindowCanvas
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas) = getWindowFrameCanvas(argServer, canvas)

  if windowFrame and canvas:
    windowFrame.orthogScrollRight(canvas)

def ncontoursIncr(argServer=None, windowPopup=None):
  """
  .. describe:: Input
    
  .. describe: Output
   
  """

  if argServer and not windowPopup:
    windowPopup = argServer.getCurrentWindowPopup()

  if windowPopup:
    windowPopup.ncontoursIncr()

def ncontoursDecr(argServer=None, windowPopup=None):
  """
  .. describe:: Input
    
  .. describe: Output
   
  """

  if argServer and not windowPopup:
    windowPopup = argServer.getCurrentWindowPopup()

  if windowPopup:
    windowPopup.ncontoursDecr()

def contourLevelsUp(argServer=None, windowPopup=None):
  """
  .. describe:: Input
    
  .. describe: Output
   
  """

  if argServer and not windowPopup:
    windowPopup = argServer.getCurrentWindowPopup()

  if windowPopup:
    windowPopup.contourLevelsUp()

def contourLevelsDown(argServer=None, windowPopup=None):
  """
  .. describe:: Input
    
  .. describe: Output
   
  """

  if argServer and not windowPopup:
    windowPopup = argServer.getCurrentWindowPopup()

  if windowPopup:
    windowPopup.contourLevelsDown()

def startPositionDelta(argServer=None, event=None):
  """
  Create marker at location and start to display delta for x, y positions.
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas, x, y) = getWindowFrameCanvasLocation(argServer, event)

  if (windowFrame and canvas and x is not None and y is not None):
    windowFrame.startPositionDelta(canvas, x, y)

def endPositionDelta(argServer=None, event=None):
  """
  Stop using marker to display delta for x, y positions.
  
  .. describe:: Input
  
  ArgumentServer, TkEvent
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas, x, y) = getWindowFrameCanvasLocation(argServer, event)

  if windowFrame:
    windowFrame.endPositionDelta()

def resetWindowSliceRange(argServer=None):
  """
  Reset the window sliceRange based on the spectra data values
  
  .. describe:: Input
  
  ArgumentServer
  
  .. describe: Output
  
  None
  """

  windowPane = argServer.getCurrentWindowPane()

  if windowPane:
    updateWindowPaneSliceRange(windowPane)
    
def calculateNoise(argServer=None, canvas=None):
  """
  Calculate noise of spectra visible in window
  
  .. describe:: Input
  
  ArgumentServer
  
  .. describe: Output
  
  None
  """

  (windowFrame, canvas) = getWindowFrameCanvas(argServer, canvas)

  if windowFrame and canvas:
    windowFrame.calculateNoise(canvas)

def toggleSpectra(argServer=None):
  """
  Turn all spectrum in window off (if any on) or on (if all off)
  
  .. describe:: Input
  
  ArgumentServer
  
  .. describe: Output
  
  None
  """

  windowPane = argServer.getCurrentWindowPane()

  if windowPane:
    toggleWindowSpectra(windowPane.spectrumWindow)

def turnSpectraOn(argServer=None):
  """
  Turn all spectrum in window on
  
  .. describe:: Input
  
  ArgumentServer
  
  .. describe: Output
  
  None
  """

  windowPane = argServer.getCurrentWindowPane()

  if windowPane:
    turnWindowSpectraOn(windowPane.spectrumWindow)

def turnSpectraOff(argServer=None):
  """
  Turn all spectrum in window off
  
  .. describe:: Input
  
  ArgumentServer
  
  .. describe: Output
  
  None
  """

  windowPane = argServer.getCurrentWindowPane()

  if windowPane:
    turnWindowSpectraOff(windowPane.spectrumWindow)

def initWindowMacros(project):
  """
  Initialise window macros, adding default ones that not already included
  
  .. describe:: Input
  
  Implementation.Project (the project)
  
  .. describe: Output
  
  None
  """

  analysisProfile = Util.getAnalysisProfile(project)

  module = __name__
  n = module.rfind('.')
  if n >= 0:
    path = module[:n].replace('.', '/')
    module = module[n+1:]
  else:
    path = ''
  ordering = 1

  macros = analysisProfile.macros
  shortcuts = [ macro.shortcut for macro in macros ]
  functions = [ macro.function for macro in macros if macro.module == module ]

  for (function, shortcut) in ( \
         ('zoomWindowOut',         'Prior'),
         ('zoomWindowIn',          'Next'),
         ('zoomWindowToSpectra',   'Z'),
         ('translateWindowUp',     'Up'),
         ('translateWindowDown',   'Down'),
         ('translateWindowLeft',   'Left'),
         ('translateWindowRight',  'Right'),
         ('zoomSliceRangeDown',    'Home'),
         ('zoomSliceRangeUp',      'End'),
         ('deleteSelected',        'Delete'),
         ('assignPeak',            'a'),
         ('showSelectedPeaks',     's'),
         ('centerWindow',          'c'),
         ('createMark',            'm'),
         ('createSidebandMarks',   'M'),
         ('createHorizontalRuler', 'h'),
         ('createSidebandHorizontalRulers',   'H'),
         ('createVerticalRuler',   'v'),
         ('createSidebandVerticalRulers',   'V'),
         ('movePeak',              'p'),
         ('movePeakAnnotation',    'q'),
         ('snapPeaks',             'P'),
         ('popupMenu',             'u'),
         ('unitePeakPositions',    'l'),
         ('unitePeakPositionsMulti','L'),
         ('orthogScrollLeft',      'j'),
         ('orthogScrollRight',     'k'),
         ('ncontoursIncr',         'i'),
         ('ncontoursDecr',         'o'),
         ('contourLevelsUp',       'e'),
         ('contourLevelsDown',     'r'),
         ('autoPeakLabelPos',      'w'),
         ('resetPeakLabelsPos',    'W'),
         ('saveProjectShortcut',   'S'),
         ('startPositionDelta',    'd'),
         ('endPositionDelta',      'D'),
         ('resetWindowSliceRange', 'z'),
         ('calculateNoise',        'N'),
         ('toggleSpectra',         't'),
         ('turnSpectraOn',         'f'),
         ('turnSpectraOff',        'g'),
       ):

    if function not in functions:
      if shortcut in shortcuts:
        shortcut = None
      macro = analysisProfile.newMacro(name=function, path=path, module=module,
                                       function=function, ordering=ordering,
                                       shortcut=shortcut)

def saveProjectShortcut(argServer):
  """
  Save the modified parts of the current project
  
  .. describe:: Input
  
  Implementation.Project
  
  .. describe: Output
  
  None
  """

  project = argServer.getProject()
  if project:
    project.saveModified()

def zoomToShowPeaks(peaks, windowPane, row=None, col=None):
  """
  Zoom in a given row and column a given window so that all of the
  input peaks are visible.
  
  .. describe:: Input
  
  List of Nmr.Peaks, Analysis.SpectrumWindowPane,
  Int (AxisPanel.axisRegion index)
  
  .. describe: Output
  
  None
  """

  window = windowPane.spectrumWindow
  stripAxis = window.stripAxis
  if stripAxis == 'x':
    if row is None:
      row = 0
      
    if col is None:
      stripPanel = windowPane.findFirstAxisPanel(label=stripAxis)
      activeRegion = stripPanel.findFirstAxisRegion(isActive=True)
      
      if activeRegion:
        col = stripPanel.sortedAxisRegions().index(activeRegion)
      else:
        col = 0
  
  else:
    if col is None:
      col = 0
  
    if row is None:
      stripPanel = windowPane.findFirstAxisPanel(label=stripAxis)
      activeRegion = stripPanel.findFirstAxisRegion(isActive=True)
      
      if activeRegion:
        row = stripPanel.sortedAxisRegions().index(activeRegion)
      else:
        row = 0  

  region = []
  spectra = set()
  for view in windowPane.spectrumWindowViews:
    spectra.add(view.analysisSpectrum.dataSource)

  minDelta = []
  for peak in peaks:
    if peak.peakList.dataSource in spectra:
      dimMapping = getPeakDimAxisMapping(peak, windowPane)
      axes = dimMapping.keys()
      axes.sort()
    
      for i, axis in enumerate(axes):
        peakDim = dimMapping[axis]
 
        if peakDim.dataDimRef:
          ppm = peakDim.value
          tol = Util.getAnalysisDataDim(peakDim.dataDim).assignTolerance

          if i >= len(region):
            region.append([ppm,ppm])
            minDelta.append(tol)
          else:
            minPpm, maxPpm = region[i]
            if ppm < minPpm:
              region[i][0] = ppm
            elif ppm > maxPpm:
              region[i][1] = ppm
            if tol >  minDelta[i]:
              minDelta[i] = tol
 
        #else: # wb104 1 Jun 2011: not sure why we ever did this and it breaks enumerate in following loop
        #  region.append(None)
  
  if not region:
    # No selected peaks are visible in this window
    return

  for i, (ppmA, ppmB) in enumerate(region):
    delta = max(minDelta[i], ppmB-ppmA)
    region[i] = ppmA-delta, ppmB+delta

  if len(region) < 2:
    axisRegion = windowPane.findFirstAxisPanel(label='y').findFirstAxisRegion()
    region.append(axisRegion.region)
 
  dx = region[0][1]-region[0][0]
  dy = region[1][1]-region[1][0]

  axisPanels = windowPane.sortedAxisPanels()
  xRegions = axisPanels[0].sortedAxisRegions()
  yRegions = axisPanels[1].sortedAxisRegions()
 
  r = xRegions[col].region
  dx2 = abs(r[1]-r[0])
  r = yRegions[row].region
  dy2 = abs(r[1]-r[0])

  windowRatio = dx2/dy2
  regionRatio = dx/dy

  if regionRatio < windowRatio:
    # region is thinner in x: do not stretch x to fit: recalc x range
    mid       = (region[0][1]+region[0][0])/2.0
    halfDx    = (dy * windowRatio) * 0.5
    region[0] = (mid-halfDx, mid+halfDx)

  else:
    # region is shorter in y: recalc y range
    mid       = (region[1][1]+region[1][0])/2.0
    halfDy    = (dx / windowRatio) * 0.5
    region[1] = (mid-halfDy, mid+halfDy)

  for i, region1d in enumerate(region):
    if region1d is not None:
      if i == 0:
        strip = col
      elif i == 1:
        strip = row
      else:
        if stripAxis == 'x':
          strip = col
        else:
          strip = row
 
      axisPanels[i].sortedAxisRegions()[strip].region = region1d

def getRepeatingAxes(windowPane):
  """
  For a spectrum window pane, gets a list of pairs of axis panels
  that are of the same type.
  
  .. describe:: Input
  
  Analysis.SpectrumWindowPane
  
  .. describe: Output
  
  List of Tuple of Analysis.AxisPanels (pairs of same type)
  """

  axisPairs = []
  
  axisPanels = windowPane.sortedAxisPanels()
 
  N = len(axisPanels)
  if N >= 2:
 
    for i in range(N-1):
      axisPanel1 = axisPanels[i]
 
      if axisPanel1.axisUnit:
        for j in range(i+1,N):
          axisPanel2 = axisPanels[j]
 
          if axisPanel2.axisUnit:
            if (axisPanel1.axisUnit.unit == axisPanel2.axisUnit.unit) \
             and  (axisPanel1.axisType == axisPanel2.axisType):
              axisPairs.append( (axisPanel1,axisPanel2) )
 
  return axisPairs
  

def findWindowTransposes(queryWindowPane, xyz):
  """
  Given a position in the input query window, find all equivalent
  points in the window which represent a swapping of coordinates for
  similarly types axes. Returns a list of mappings for each target
  window, comprising a mapping name (e.g. for GUIs) the target
  spectrum window and the target position.
  
  .. describe:: Input
  
  Analysis.SpectrumWindowPane, List of Floats (position)
  
  .. describe: Output
  
  List of Lists [String (Mapping Name), Analysis.SpectrumWindow,
  List of Floats (position)]
  """

  axisNames = ['x','y','z1','z2','z3','z4']
  
  axes = getRepeatingAxes(queryWindowPane) 
  if not axes:
    return []
  
  mappings  = []
  axisPanels = queryWindowPane.sortedAxisPanels()
    
  n = len(axisPanels)
  if n >= 2:
    coords = xyz[:n]
 
    for axisPanel1, axisPanel2 in axes:
      i = axisPanels.index(axisPanel1)
      j = axisPanels.index(axisPanel2)
      name1 = axisNames[i]
      name2 = axisNames[j]
      pos1 = coords[i]
      pos2 = coords[j]
 
      data = (axisPanel1.axisType.name, name1, name2, pos2, pos1)
      label = '%s %s-%s: %.3f %.3f' % data
      xyz2 = list(coords)
      xyz2[i] = coords[j]
      xyz2[j] = coords[i]
 
      posDict = {}
      for k, coord in enumerate(xyz2):
        posDict[axisNames[k]] = coord
 
      mappings.append( [label, posDict] )
  
  return mappings

def findOrthogonalWindows(queryWindowPane, positions, excludeQuery=True, minDims=3):
  """
  Given positions in the input query window, find all equivalent
  points in all spectrum windows which have either a) a Z dimension
  that matches any one of the dimensions of the query window or b)
  two dimemsions that match the query window. Returns a list of
  mappings for each target window, comprising a mapping name (e.g.
  for GUIs) the target spectrum window and the target position.
  There is boolean option, to say whether the query window mapping
  is to be excluded form the output, this is usually true for
  searches originating from windows, but may not be so  if the query
  comes from a table. The minDims optional argument determines the
  minimum number of dimensions that a target window must have to be
  considered.
  
  .. describe:: Input
  
  Analysis.SpectrumWindowPane, List of List of Floats (position),
  Boolean, Int
  
  .. describe: Output
  
  List of Lists [String (Mapping Name), Analysis.SpectrumWindow,
  List of Floats (position)]
  """

  N = len(queryWindowPane.axisPanels)
  if N < 2:
    return []
  
  analysisProject = queryWindowPane.topObject
  axisUnits = []
  axisTypes = []
  axisNames = ['x','y','z1','z2','z3','z4']
  coords    = [xyz[:N]for xyz in positions]
  mappings  = []
  labels    = {}

  windowZPlanes = []
  
  for axisPanel in queryWindowPane.sortedAxisPanels():
    if axisPanel.axisUnit:
      axisUnits.append(axisPanel.axisUnit.unit)
      axisTypes.append(axisPanel.axisType)
  
  N = len(axisTypes)
  
  # recursive routine to find all the different matching axis combinations
  def traceRoute(N, axisTypes2, axisUnits2, i, paths, path):
    if i < N:
      match = False
      for j in range(len(axisTypes2)):
        if j not in [x[1] for x in path]:
          if (axisTypes[i] is axisTypes2[j]) and (axisUnits[i] is axisUnits2[j]):
            traceRoute(N, axisTypes2, axisUnits2, i+1, paths, path[:] + [(i,j),])
            match = True
          
      if not match:
        traceRoute(N, axisTypes2, axisUnits2, i+1, paths, path)
      
    elif path:
      paths.append(tuple(path))
      
  # setup labels to show relation to query window
  for i in range(N):
    if axisTypes.count(axisTypes[i]) > 1:
      labels[i] = '(F%d)' % (i+1)
  
  windowsByName = [(w.name,w) for w in analysisProject.spectrumWindows]
  windowsByName.sort()
  windowsByName = [x[1] for x in windowsByName]
  
  # find windows with matching dims
  windowPanes = []
  for window in windowsByName:
    for windowPane in window.sortedSpectrumWindowPanes():
      axisPanels = windowPane.sortedAxisPanels()
      
      M = len(axisPanels)
      if M < minDims:
        continue
 
      axisTypes2 = [ ap.axisType for ap in axisPanels if ap.axisUnit ]
      axisUnits2 = [ ap.axisUnit.unit for ap in axisPanels if ap.axisUnit ]
 
      paths = []
      traceRoute(N, axisTypes2, axisUnits2, 0, paths, [])
 
      paths = [p for p in paths if len(p) > 1]
 
      if paths:
        windowPanes.append(windowPane)
        mappings.append(paths) 
 
  # create the mapping names and find the mapped positions
  for p, windowPane in enumerate(windowPanes):
    paths  = mappings[p]
    window = windowPane.spectrumWindow
    
    axisPanels = windowPane.sortedAxisPanels()
    for path in paths:
      if excludeQuery and (windowPane is queryWindowPane):
        for i, j in path:
          if i != j:
            break
        else: # no transformation
          continue
          
      positions  = [{} for xyz in coords]
      axisLabels = ['-' for ap in axisPanels]
 
      for i, j in path:
        label = ''
        if len(paths) > 1:
          label = labels.get(i) or ''
        axisLabels[j] = label + axisPanels[j].axisType.name

        for k in range(len(coords)):
          position = coords[k][i]
          if position is not None:
            positions[k][axisNames[j]] = position
      data = (' '.join(axisLabels), getWindowPaneName(windowPane))
      mappingName = '%s in %s' % data
        
      windowZPlanes.append( [mappingName, windowPane, positions] )

  return windowZPlanes

def findDataDimAxisMapping(spectrumWindowView, dataDim):

  """
  Get the AxisMapping object for a given dataDim in a given spectrumWindowView.
  
  .. describe:: Input
  
  Analysis.SpectrumWindowView, Nmr.AbstractDataDim
  
  .. describe: Output
  
  Analysis.AxisMapping
  """

  for axisMapping in spectrumWindowView.axisMappings:
    if axisMapping.analysisDataDim.dataDim is dataDim:
      return axisMapping

  return None

def getDataDimAxisMapping(spectrum, windowPane):
  """
  Find a mapping between axes in a spectrum window pane and the data
  dimensions of a spectrum: Works initially with dataDim matches.
  
  .. describe:: Input
  
  Nmr.DataSource, Analysis.SpectrumWindow
  
  .. describe: Output
  
  Dictionary of Nmr.DataDims keyed by 'x','y','z1' etc
  """
  
  mapping  = {}
  
  analysisSpectrum = Util.getAnalysisSpectrum(spectrum)
  
  view = analysisSpectrum.findFirstSpectrumWindowView(spectrumWindowPane=windowPane)

  if view:
    for axisMapping in view.axisMappings:
      dataDim = axisMapping.analysisDataDim.dataDim
      if dataDim:
        mapping[axisMapping.label] = dataDim
  
  else:
    view = analysisSpectrum.findFirstSpectrumWindowView()
    
    if view:
      panelTypes = {}
      for analysisDataDim in analysisSpectrum.analysisDataDims:
        axisMapping = view.findFirstAxisMapping(analysisDataDim=analysisDataDim)
        if axisMapping:
          axisPanel = windowPane.findFirstAxisPanel(label=axisMapping.label)
          if axisPanel:
            panelTypes[axisPanel.panelType] = analysisDataDim.dataDim
  
      for axisPanel in windowPane.axisPanels:
        dataDim = panelTypes.get(axisPanel.panelType)
        if dataDim:
          mapping[axisPanel.label] = dataDim
  
    else:
      dataDims = list(spectrum.dataDims)
      for axisPanel in windowPane.axisPanels:
        axisType = axisPanel.axisType
        isotopes = axisType.isotopeCodes
        
        for dataDim in dataDims:
          if isinstance(dataDim, Nmr.FreqDataDim):
            dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
            expDimRef =  dataDimRef.expDimRef
            
            if dataDimRef and (expDimRef.isotopeCodes == isotopes) \
               and (expDimRef.measurementType == axisType.measurementType):
              mapping[axisPanel.label] = dataDim
              dataDims.remove(dataDim)
              break
  
  return mapping

def getPeakDimAxisMapping(peak, windowPane):
  """
  Find a mapping between axes in a spectrum window pane
  and the dimensions of a peak:
  Works initially with dataDim matches.
  
  .. describe:: Input
  
  Nmr.Peak, Analysis.SpectrumWindowPane
  
  .. describe: Output
  
  Dictionary of Nmr.PeakDims keyed by 'x','y','z1' etc
  """
  
  dataDimMapping = getDataDimAxisMapping(peak.peakList.dataSource, windowPane)
  
  mapping  = {}
  for key in dataDimMapping.keys():
    dataDim = dataDimMapping[key]
    mapping[key] = peak.findFirstPeakDim(dim=dataDim.dim)
  
  return mapping


def getPeakWindowPosition(peak, windowPane, default=None, useDefault=True):
  """
  Find the position of a peak in a given window's coordinates
  Copes with peakDims ordering being different to window axis order.
  If peaks do not match a window axis a default position may be used otherwise
  the default position is the current window location. If useDefault is False
  None will be entered instead of a default location.
  
  .. describe:: Input
  
  Nmr.Peak, Analysis.SpectrumWindowPane
  
  .. describe: Output
  
  List of Floats
  """

  dimMapping = getPeakDimAxisMapping(peak, windowPane)
  
  i = 0
  position = []
  for axisPanel in windowPane.sortedAxisPanels():
    peakDim = dimMapping.get(axisPanel.label)
    if peakDim:
      if axisPanel.axisUnit:
        value = getPeakDimPosition(peakDim, toUnit=axisPanel.axisUnit.unit)
      else:
        value = peakDim.position
      position.append( value )

    elif not useDefault:
      position.append( None )

    elif default is not None:
      position.append( default[i] )

    else:
      axisRegion = axisPanel.findFirstAxisRegion()
      for axisRegion0 in axisPanel.sortedAxisRegions():
        if axisRegion0.isActive:
          axisRegion = axisRegion0
          break
      
      (r0, r1) = axisRegion.region
      position.append( 0.5 * (r1 - r0) )
      
    i += 1
  
  return position
  
def displayPeakStrips(analysis, peaks, windowPane):
  """
  Create a number of strips at peak locations along a given strip axis.
  The ortho value of the strips will be the mean peak position for
  the orthogonal axis. A new window may be created if needed. 
  
  .. describe:: Input
  
  Analysis (Tkinter top level), List of Nmr.Peaks, Analysis.SpectrumWindow
  
  .. describe: Output
  
  Analysis.SpectrumWindow
  """

  if windowPane:
    stripAxis = windowPane.spectrumWindow.stripAxis
  else:
    stripAxis = 'x'  

  orthoLabel = 'y'
  if stripAxis == 'y':
    orthoLabel = 'x'

  if not (peaks and windowPane):
    return


  if len(peaks) > 15:
    peaks = peaks[:15]

  
  orthoPos = 0.0
  orthoNum = 0.0
  positions = []
  
  for peak in peaks:
    dimMapping = getPeakDimAxisMapping(peak, windowPane)
    peakDim = dimMapping.get(orthoLabel)
    
    if peakDim:
      orthoPos += peakDim.realValue
      orthoNum += 1.0

    position = {}
    for axisPanel in windowPane.axisPanels:
      label   = axisPanel.label
      peakDim = dimMapping.get(label)
      if peakDim:
        position[label] = peakDim.realValue
 
    positions.append(position)

  spectrum = peaks[0].peakList.dataSource
 
  if orthoNum:
    orthoPos /= orthoNum
    for position in positions:
       position[orthoLabel] = orthoPos
  
  displayStrips(analysis, positions, spectrum=spectrum, windowPane=windowPane)
  

def swapStripAxis(analysis, window):
  """
  Swap vertical strips for horizontal ones and vice versa. Preserves number
  active index and z positions
  
  .. describe:: Input
  
  Analysis (Tkinter top level), Analysis.SpectrumWindow
  
  .. describe: Output
  
  None
  """

  stripAxis = window.stripAxis
  if stripAxis == 'x':
    orthoAxis = 'y'
  else:
    orthoAxis = 'x'  
  
  window.stripAxis = orthoAxis 

  popup = analysis.getWindowPopup(window.name, doOpen=True)

  windowFrames = popup.windowFrames
  for p, windowPane in enumerate(window.sortedSpectrumWindowPanes()):
    windowFrame = windowFrames[p]
    stripPanel = windowPane.findFirstAxisPanel(label=stripAxis)
    orthoPanel = windowPane.findFirstAxisPanel(label=orthoAxis)

    axisRegionsS = stripPanel.sortedAxisRegions()
    axisRegionsD = orthoPanel.sortedAxisRegions()
 
    activeS = [r.isActive for r in axisRegionsS]
    activeD = [r.isActive for r in axisRegionsD]
 
    nSource = len(axisRegionsS)
    nDest   = len(axisRegionsD)
 
    if p == 0:
      # TBD: Assumes we only want to do this once
      # do delete before add
      if nSource > nDest:
        if stripAxis == 'x':
          windowFrame.setColsTo(nDest)
          windowFrame.setRowsTo(nSource)
        else:
          windowFrame.setRowsTo(nDest)
          windowFrame.setColsTo(nSource)
 
      elif nDest > nSource:
        if stripAxis == 'x':
          windowFrame.setRowsTo(nSource)
          windowFrame.setColsTo(nDest)
        else:
          windowFrame.setColsTo(nSource)
          windowFrame.setRowsTo(nDest)

    axisRegionsS = stripPanel.sortedAxisRegions()
    axisRegionsD = orthoPanel.sortedAxisRegions()

    for i in range(nDest):
      axisRegionsS[i].isActive = activeD[i]

    for i in range(nSource):
      axisRegionsD[i].isActive = activeS[i]
  
  ###window.stripAxis = orthoAxis 

def displayStrips(analysis, positions, orthoPositions=None,
                  spectrum=None, windowPane=None):
  """
  Create a number of strips at specified positions in a new or
  existing spectrum window panel. bPositions are specified in
  dictionary form. e.g. {'x': 7.43, 'y': 4.21, 'z1': 121.96}
  Orthogonal positions may be specified to further divide strips 
  nto cells. A spectrum is specified so that the appropriate
  type of window can be made if needed.
  
  .. describe:: Input
  
  Analysis (Tkinter top level), List of Dicts (axis:location),
  List of Floats, Nmr,DataSource, Analysis.SpectrumWindowPanel
  
  .. describe: Output
  
  Analysis.SpectrumWindow
  """
  assert spectrum or windowPane

  if windowPane:
    stripAxis = windowPane.spectrumWindow.stripAxis
  else:
    stripAxis = 'x'
  
  if stripAxis == 'x':
    orthoAxis = 'y'
    
  else:
    orthoAxis = 'x'  
    numStrips = 1

  if orthoPositions is not None:
    numDivs = len(orthoPositions)
  elif windowPane:
    numDivs = len(windowPane.findFirstAxisPanel(label=orthoAxis).axisRegions)
  else:
    numDivs = 1
    
  numStrips = len(positions)
  
  project = analysis.getProject()
  # analysis is the guiParent/top level object

  turnedOff = False
  try:

    if windowPane is None:
      # make a new window with the desired strips
      analysisProject = project.currentAnalysisProject

      i = len(analysisProject.spectrumWindows)
      name = 'window%s' % i
      while analysisProject.findFirstSpectrumWindow(name = name):
        i += 1
        name = 'window%s' % i

      axisTypes = []
      for i in range(spectrum.numDim):
        isotopeCodes = spectrum.sortedDataDims()[i].expDim.findFirstExpDimRef().isotopeCodes
        axisType = analysisProject.findFirstAxisType(isotopeCodes = isotopeCodes)
        axisTypes.append(axisType)

      window = createSpectrumWindow(project, name, [axisTypes,], 
                                    ncols=numStrips, nrows=numDivs)
      popup = analysis.getWindowPopup(window.name)
      windowPane = window.findFirstSpectrumWindowPane()
      windowFrame = popup.windowFrames[0]
      
      turnedOff = popup.turnDrawRequestsOff()

    else:
      # add or remove strips from an existing window
      window = windowPane.spectrumWindow
      axisRegions = windowPane.findFirstAxisPanel(label=stripAxis).sortedAxisRegions()
      nStrips = len(axisRegions)
      popup = analysis.getWindowPopup(window.name, doOpen=True)
      windowFrame = windowPane.getWindowFrame()
      turnedOff = popup.turnDrawRequestsOff()

      if stripAxis == 'x':
        windowFrame.setColsTo(numStrips)
        windowFrame.setRowsTo(numDivs)
      else:
        windowFrame.setColsTo(numDivs)
        windowFrame.setRowsTo(numStrips)

    for num in range(numStrips):
      # locate strips and any cells
      position = positions[num]
    
      if orthoPositions:
        for div in range(numDivs):
          position[orthoAxis] = orthoPositions[div]
          if stripAxis == 'x':
            row=div
            col=num
          else:
            row=num
            col=div
        
          #analysis.gotoPosition(window_name=window.name, position=position, row=row, col=col, doOpen=False)
          windowFrame.gotoPosition(position, row=row, col=col, doLift=False)
   
      else:
        if stripAxis == 'x':
          #analysis.gotoPosition(window_name=window.name, position=position, col=num, doOpen=False)
          windowFrame.gotoPosition(position, col=num, doLift=False)
        else:
          #analysis.gotoPosition(window_name=window.name, position=position, row=num, doOpen=False)
          windowFrame.gotoPosition(position, row=num, doLift=False)

  finally:

    if turnedOff:
      popup.turnDrawRequestsOn()

    popup.open()

  return windowPane

def getSpinSystemWindowShifts(spinSystems, windowPane, shiftList=None):
  """
  Get the shifts (linked to primary/z axes) visible in a given
  window pane that are represented by the resonances of a spin
  system. Also returns a list of orthogonal axis shifts.
  
  .. describe:: Input
  
  List of Nmr.ResonanceGroups, Analysis.SpectrumWindowPane,
  Nmr.ShiftList
  
  .. describe: Output
  
  List of Dict of Word:Nmr.Shift (dimName:position eg 'z1':Shift),
  List of Nmr.Shifts
  """
  from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes
  from ccpnmr.analysis.core.AssignmentBasic import getDataDimFullShiftRange, getBoundResonances
  from ccpnmr.analysis.core.MoleculeBasic   import areResonancesBound
  
  stripAxis = windowPane.spectrumWindow.stripAxis
  sortedAxisPanels = windowPane.sortedAxisPanels()
  xIsotopes = sortedAxisPanels[0].axisType.isotopeCodes
  yIsotopes = sortedAxisPanels[1].axisType.isotopeCodes
  zIsotopeList = [ ap.axisType.isotopeCodes for ap in sortedAxisPanels[2:] ]
  
  spectra = []
  activeViews = [view for view in windowPane.spectrumWindowViews if view.isPosVisible]
  for view in activeViews:
    spectra.append(view.analysisSpectrum.dataSource)
  
  # below does not work because spinSystems might be empty
  #nmrProject = spinSystems[0].topObject
  # so do this instead
  nmrProject = windowPane.topObject.nmrProject
  
  ranges = []
  if not shiftList:
    if spectra:
      for spectrum in spectra:
        if spectrum.experiment.shiftList:
          shiftList = spectrum.experiment.shiftList
          break

    else:  
      shiftList = nmrProject.findFirstMeasurementList(className='ShiftList')
  
  for spectrum in spectra[:1]:
    dimMapping = getDataDimAxisMapping(spectrum, windowPane)
    axisNames = dimMapping.keys()
    axisNames.sort()
    
    for key in axisNames:
      dataDim = dimMapping[key]
      ranges.append( getDataDimFullShiftRange(dataDim) )

  for spectrum in spectra[1:]:
    dimMapping = getDataDimAxisMapping(spectrum, windowPane)
    axisNames = dimMapping.keys()
    axisNames.sort()

    i = 0
    for key in axisNames:
      dataDim = dimMapping[key]
      minVal, maxVal = getDataDimFullShiftRange(dataDim)
      if i < len(ranges):
        if ranges[i][0] > minVal:
          ranges[i][0] = minVal
 
        if ranges[i][1] < maxVal:
          ranges[i][1] = maxVal
      
      else:
        ranges.append([minVal, maxVal])
      
      i += 1
 
  if not ranges:
    ranges = [axisPanel.axisType.region for axisPanel in sortedAxisPanels]
 
  sortShifts = []
  xShifts = []
  yShifts = []
  zListShifts = [ [] for axisPanel in sortedAxisPanels[2:] ]

  for j, spinSystem in enumerate(spinSystems):
    resonanceDict = {}

    for resonance in spinSystem.resonances:
      resonanceDict[resonance] = True

      for resonance2 in getBoundResonances(resonance):
        resonanceDict[resonance2] = True

    for resonance in resonanceDict.keys():
      shift = resonance.findFirstShift(parentList=shiftList)
      if shift:
        if resonance.isotopeCode in xIsotopes:
          if (shift.value >= ranges[0][0]) and (shift.value <= ranges[0][1]):
            if stripAxis == 'x':
              sortShifts.append( (j, shift) )
            xShifts.append( shift )

        if resonance.isotopeCode in yIsotopes:
          if (shift.value >= ranges[1][0]) and (shift.value <= ranges[1][1]):
            if stripAxis == 'y':
              sortShifts.append( (j, shift) )
            yShifts.append( shift )

        for i, zIsotopes in enumerate(zIsotopeList):
          # Could have a 2D mapped to a 3D window
          if i+2 >= len(ranges):
            break

          if resonance.isotopeCode in zIsotopes:
            if (shift.value >= ranges[i+2][0]) \
               and (shift.value <= ranges[i+2][1]):
              zListShifts[i].append( shift )

  if stripAxis == 'x':
    orthoShifts = yShifts
  else:
    orthoShifts = xShifts

  shifts    = []
  for j, shift in sortShifts:
    shiftDict = {stripAxis:shift}

    foundZ = False
    for i, zShifts in enumerate(zListShifts):
      for zShift in zShifts:
        if areResonancesBound(shift.resonance,zShift.resonance):
          shiftDict['z%d' % (i+1,)] = zShift
          foundZ = True
          break

    if not foundZ:
      for orthoShift in orthoShifts:
        # Maybe:
        #if not areResonancesBound(orthoShift.resonance,shift.resonance):
        #  continue

        for i, zShifts in enumerate(zListShifts):
          for zShift in zShifts:
            if areResonancesBound(orthoShift.resonance,zShift.resonance):
              shiftDict['z%d' % (i+1,)] = zShift
              foundZ = True
              break

    if zIsotopeList and not foundZ:
      continue

    shifts.append((j,shift.value, shiftDict))

  shifts.sort()
  shifts.reverse()
  shifts = [x[2] for x in shifts]

  return shifts, orthoShifts

def gotoResonancePosition(analysis, resonance, windowPane, shiftList, row=None, col=None):
  """
  Go to the position of a given resonance (according to the input
  shiftlist) in a given window. Note currently only goes to first
  available axis not all.
  
  .. describe:: Input
  
  Analysis (Tkinter top level), Nmr.Resonance,
  Analysis.SpectrumWindowPane, Nmr.ShiftList
  
  .. describe: Output
  
  None
  """

  shift = resonance.findFirstShift(parentList=shiftList)

  if shift:
    axisPanels = []
    for axisPanel in windowPane.axisPanels:
      if axisPanel.axisType:
        if resonance.isotopeCode in axisPanel.axisType.isotopeCodes:
          axisPanels.append(axisPanel)
          #break
 
    if len(axisPanels) > 1:
      stripAxis = windowPane.spectrumWindow.stripAxis
      
      axisLabels = []
      for axisPanel in axisPanels:
        label = axisPanel.label
        
        if (label in 'xy') and (label != stripAxis):
          continue
        
        else:
          axisLabels.append(label) 
   
    else:   
      axisLabels = [ap.label for ap in axisPanels]
 
    position = {}
    for axis in axisLabels:
      position[axis] = shift.value

    if hasattr(windowPane, 'getWindowFrame'):
      windowFrame = windowPane.getWindowFrame()
      windowFrame.gotoPosition(position, row=row, col=col)

def displaySpinSystemStrips(analysis, spinSystems, windowPane,
                            shiftList=None, splitIntoCells=False):
  """
  For a given window pane, display the strips appropriate for a
  given spin system. Option to further divide the strips
  along the y axis according to spin system shifts.
  
  .. describe:: Input
  
  Analysis (Tkinter top level), Nmr.ResonanceGroup,
  Analysis.SpectrumWindowPane, Nmr.ShiftList, Boolean
  
  .. describe: Output
  
  None
  """
  
  shifts, orthoShifts = getSpinSystemWindowShifts(spinSystems, windowPane,
                                                  shiftList=shiftList)
  
  if not shifts:
    serials = ','.join([str(ss.serial) for ss in spinSystems])
    msg = 'No possible strips in window %s for spin systems [%s]'
    data = (windowPane.spectrumWindow.name, serials)
    showWarning('Spin System Strips Failure', msg % data)
    return
      
  orthoPositions = None
  if splitIntoCells and orthoShifts:
    orthoPositions = [shift.value for shift in orthoShifts]
    orthoPositions.sort()
  
  for shiftDict in shifts:
    for key in shiftDict.keys():
      shiftDict[key] = shiftDict[key].value
  
  stripAxis = windowPane.spectrumWindow.stripAxis
  stripPanel = windowPane.findFirstAxisPanel(label=stripAxis)
  activeRegion = stripPanel.findFirstAxisRegion(isActive=True)
      
  displayStrips(analysis, shifts, orthoPositions=orthoPositions, windowPane=windowPane)

  # Restore originally active
  for region in stripPanel.axisRegions:
    if region is activeRegion:
      region.isActive = True
    else:
      region.isActive = False
      

def getAxisRegionRegion(axisRegion):
  """
  .. describe:: Input
    
  .. describe: Output
   
  """

  # note: inOverrideMode is a runtime thing
  # it is set if both window.useOverrideRegion (an API attribute) and we are
  # in a mode (e.g. printing) where this applies
  window = axisRegion.axisPanel.spectrumWindowPane.spectrumWindow
  if hasattr(window, 'inOverrideMode') and window.inOverrideMode \
      and axisRegion.overrideRegion:
    region = axisRegion.overrideRegion
  else:
    region = axisRegion.region

  return region

def getPossibleTiedAxisPanels(axisPanel):

  window = axisPanel.spectrumWindowPane.spectrumWindow
  axisType = axisPanel.axisType
  analysisProject = window.analysisProject

  data = []
  for tiedWindow in analysisProject.spectrumWindows:
    if tiedWindow is window:
      continue
    tiedWindowPane = tiedWindow.findFirstSpectrumWindowPane()
    if tiedWindowPane: # should always be the case
      for tiedAxisPanel in tiedWindowPane.axisPanels:
        if tiedAxisPanel.axisType is axisType:
          data.append(('%s:%s' % (tiedWindow.name, tiedAxisPanel.label), tiedAxisPanel))

  data.sort()
  axisPanels = [d[1] for d in data]

  return axisPanels

TIED_AXIS = 'tiedAxis'

def getTiedAxisPanels(axisPanel):

  window = axisPanel.spectrumWindowPane.spectrumWindow
  analysisProject = window.analysisProject
  project = analysisProject.root
  application = project.application
  tiedAxisNames = application.getValues(axisPanel, keyword=TIED_AXIS)
  axisPanels = []
  for tiedAxisName in tiedAxisNames:
    fields = tiedAxisName.split(':')
    if len(fields) == 2:
      (tiedWindowName, tiedLabel) = fields
      tiedWindow = analysisProject.findFirstSpectrumWindow(name=tiedWindowName)
      if tiedWindow:
        tiedWindowPane = tiedWindow.findFirstSpectrumWindowPane()
        if tiedWindowPane: # should always be True
          tiedAxisPanel = tiedWindowPane.findFirstAxisPanel(label=tiedLabel)
          if tiedAxisPanel:
            axisPanels.append(tiedAxisPanel)

  return axisPanels

def setTiedAxisPanels(axisPanel, tiedAxisPanels):

  values = ['%s:%s' % (tiedAxisPanel.spectrumWindowPane.spectrumWindow.name, tiedAxisPanel.label) for tiedAxisPanel in tiedAxisPanels]
  project = axisPanel.root
  application = project.application
  application.setValues(axisPanel, TIED_AXIS, values)

VALUE_AXIS_OFFSET = 'valueAxisOffset'

def getValueAxisOffset(view):
  
  application = view.root.application
  
  return application.getValue(view, keyword=VALUE_AXIS_OFFSET) or 0.0
  
def setValueAxisOffset(view, value):
  
  application = view.root.application
 
  application.setValue(view, keyword=VALUE_AXIS_OFFSET, value=value)
  
X_AXIS_OFFSET = 'xAxisOffset'

def getXAxisOffset(view):
  
  application = view.root.application
  
  return application.getValue(view, keyword=X_AXIS_OFFSET) or 0.0
  
def setXAxisOffset(view, value):
  
  application = view.root.application
 
  application.setValue(view, keyword=X_AXIS_OFFSET, value=value)

  

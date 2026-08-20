
"""
======================COPYRIGHT/LICENSE START==========================

WindowFrame.py: Part of the CcpNmr Analysis program

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
import time
import sys
import Tkinter

from memops.general import Implementation

from memops.universal.Region1D import Region1D
from memops.universal.Util import formatDecimals, isWindowsOS, useWheelMouse, OSButton, isMacOS

from memops.gui.ButtonScrollbar     import ButtonScrollbar
from memops.gui.Color               import hexToRgb, inverseGrey, hexXor, hexInvert
from memops.gui.DataEntry           import askString
from memops.gui.Frame               import Frame
from memops.gui.MessageReporter     import showError, showOkCancel, showWarning, showYesNo
from memops.gui.RegionSelector      import RegionSelector
from memops.gui.Scrollbar           import top_left_trough_mode, bottom_right_trough_mode
from memops.gui.ScrolledWindow      import ScrolledWindow, no_key_state, \
                                           shift_key_state, ctrl_key_state, alt_key_state

from ccpnmr.analysis.core.ExperimentBasic import getPrimaryDataDimRef, calculateNoiseInBox, getRegionStats

try:
  import ccpnmr.c.ContourFile as ContourFile
  import ccpnmr.c.ContourStyle as ContourStyle
  import ccpnmr.c.ContourLevels as ContourLevels
  import ccpnmr.c.PeakList as PeakList
  import ccpnmr.c.SliceFile as SliceFile
  import ccpnmr.c.WinPeakList as WinPeakList
except Exception, e:
  print 'Error, the WindowFrame module will not work, something is wrong with the C code.'
  print 'Exception:', e
  print 'Will continue without Analysis window drawing functionality'
  ContourFile = ContourStyle = ContourLevels = PeakList = SliceFile = WinPeakList = None

try:
  import memops.c.GlHandler as GlHandler
except:
  GlHandler = None
try:
  import memops.c.TkHandler as TkHandler
except:
  TkHandler = None
if (not GlHandler and not TkHandler):
  print 'Error, the WindowFrame module will not work, something is wrong with both the GlHandler and TkHandler C code.'

from ccpnmr.analysis.core.AssignmentBasic import clearSeqSpinSystemLinks, addPeakResonancesToSeqSpinSystems, propagatePeakAssignments, addPeakResonancesToSpinSystem
from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes, getDataDimRefFullRange
from ccpnmr.analysis.core.MarkBasic import createPeakMark, createNonPeakMark, createRuler, \
                                      removeMarks, removeRulers
from ccpnmr.analysis.core.PeakBasic import addPeak, addPeakToSelected, getPeakDimPosition, \
                                      removePeakFromSelected, findMatchingPeaks, \
                                      setPeakPosition, findPeaks, searchPeaks, removePeaksAssignment, \
                                      makePeakAnnotation, makeIntermediatePeak, propagatePeakUnaliasing, \
                                      findPositionMatches, snapPeaks, findSymmetryPeaks, \
                                      uniteResonancePeakPositions, PEAK_FIT_METHODS, fitPeaks, \
                                      copyPeaksToPeakList
                                      
from ccpnmr.analysis.core.PeakFindParams import getPeakFindParams
from ccpnmr.analysis.popups.SetReference import SetReferencePopup
from ccpnmr.analysis.popups.TranslatePeak import TranslatePeakPopup
from ccpnmr.analysis.core import Util
from ccpnmr.analysis.core.WindowBasic import displayPeakStrips, getEquivalentWindows, displayStrips, \
                                        isSpectrumInWindowPane, getSpectrumWindowView, setWindowPeakListState, \
                                        createSpectrumWindow, centerAxisRegion, getPeakWindowPosition, \
                                        zoomToShowPeaks, findOrthogonalWindows, findWindowTransposes, \
                                        refinePeakLabelPositions, resetPeakLabelPositions, \
                                        getDataDimAxisMapping, getWindowPaneName, VALUE_AXIS_OFFSET, X_AXIS_OFFSET

from ccpnmr.analysis.core.WindowDraw  import WindowDraw

# maximum time do double buffer before switch to single buffer
doubleBufferTimeout = 0.5 # seconds

slice_height_entries = [50, 100, 150, 200]

all_label = 'All'

PEAK_FIT_TYPES = ('singly', 'together')

# internal class
class WindowTimeoutException:
  pass

class FakeSlice:
  pass

# public class
class WindowFrame(Frame, WindowDraw):

  def __init__(self, parent, windowPane):

    self.windowPopup = parent.parent
    self.topPopup = self.windowPopup.parent

    self.drawCanvasCount = 0
    self.waitDraw     = False
    self.waitPeak     = False
    self.waitResize   = False
    self.keypressTime = None
    self.updatingAspectRatio = False
    self.crosshairLocations = []
    self.fontChanged = False
    self.canvasXs = None
    self.canvasYs = None
    self.keysym = None
    self.draggedBox = None
    # below now done in WindowDraw
    #self.deltaPositions = None

    profile = Util.getAnalysisProfile(windowPane.root)
    application = profile.root.application
    expand_y_axis_tick_decimals = application.getValue(profile, 'yAxisDecimalsExpand', defaultValue=True)

    if TkHandler and not GlHandler:
      handlerClass = TkHandler.TkHandler
      useGl = False
    elif GlHandler and not TkHandler:
      handlerClass = GlHandler.GlHandler
      useGl = True
    elif GlHandler and TkHandler:
      graphicsHandler = profile.graphicsHandler
      
      # NOTE:ED - MacOS only works with Tk
      #           Windows only works with OpenGL
      #             not checked linux yet.
      thisOS = sys.platform[:3].lower()
      if thisOS == 'win':
        handlerClass = GlHandler.GlHandler
        useGl = True
      elif thisOS == 'dar':
        handlerClass = TkHandler.TkHandler
        useGl = False
      else:
        # linux/other version here - use graphics handler option
        if (graphicsHandler == 'OpenGL'):
          handlerClass = GlHandler.GlHandler
          useGl = True
        else:
          handlerClass = TkHandler.TkHandler
          useGl = False

    else:
      useGl = False  # arbitrary
      handlerClass = None

    if useGl:
      if hasattr(self.topPopup, 'glDirect'):
        handlerArgs = (self.topPopup.glDirect,)
      else:
        handlerArgs = ()
      handlerExc = GlHandler.error
      # below was for when we had GL_EQUIV
      #handlerXorColor = (0.0, 0.0, 0.0)
      # below is for GL_XOR
      handlerXorColor = (1.0, 1.0, 1.0)
    else:
      handlerArgs = ()
      handlerExc = TkHandler and TkHandler.error
      handlerXorColor = (1.0, 1.0, 1.0)

    self.handlerClass    = handlerClass
    self.handlerArgs     = handlerArgs
    self.handlerExc      = handlerExc
    self.handlerXorColor = handlerXorColor # for xor mode

    xname = windowPane.findFirstAxisPanel(label='x').axisType.name
    yname = windowPane.findFirstAxisPanel(label='y').axisType.name
    self.xyName = '(%s,%s)' % (xname,yname)

    #print 'handler info', handlerClass, handlerXorColor
    
    # WindowDraw.__init__ sets self.windowPane
    WindowDraw.__init__(self, parent, windowPane)
    Frame.__init__(self, parent)


    self.grid_rowconfigure(0, weight=1)
    self.grid_columnconfigure(0, weight=1)

    # Scrolled Window

    world_region = self.getWorldRegion()
    xview_regions = self.getXviewRegions(world_region.x_region)
    yview_regions = self.getYviewRegions(world_region.y_region)

    xpanel = windowPane.findFirstAxisPanel(label='x')
    ypanel = windowPane.findFirstAxisPanel(label='y')

    ncols = len(xpanel.axisRegions)
    nrows = len(ypanel.axisRegions)
    width = [ axisRegion.size for axisRegion in xpanel.sortedAxisRegions() ]
    height = [ axisRegion.size for axisRegion in ypanel.sortedAxisRegions() ]
    xslice = xpanel.slicePanel
    yslice = ypanel.slicePanel
    show_xslices = xslice.isVisible
    show_yslices = yslice.isVisible
    xslice_size = Util.greaterOrEqualEntry(xslice.thickness, slice_height_entries)
    yslice_size = Util.greaterOrEqualEntry(yslice.thickness, slice_height_entries)
    show_xscrollbar = xpanel.isVisible
    show_yscrollbar = ypanel.isVisible
    self.scrolled_window = ScrolledWindow(self, ncols, nrows,
                             width=width, height=height,
                             world_region=world_region,
                             xview_regions=xview_regions,
                             yview_regions=yview_regions,
                             configure_func=self.resize,
                             expose_func=self.expose,
                             motion_func=self.motion,
                             leave_func=self.leave,
                             slice_configure_func=self.sliceResize,
                             slice_expose_func=self.sliceExpose,
                             slice_motion_func=self.sliceMotion,
                             slice_leave_func=self.leave,
                             view_set_func=self.viewSet,
                             key_func=self.keypress,
                             slice_key_func=self.sliceKeypress,
                             select_single_func=self.selectSingle,
                             select_multi_func=self.selectMulti,
                             drag_func=self.dragBox,
                             show_xslices=show_xslices,
                             show_yslices=show_yslices,
                             xslice_size=xslice_size,
                             yslice_size=yslice_size,
                             get_aspect_ratio = self.getAspectRatio,
                             get_geometry=self.windowPopup.getGeometry,
                             set_geometry=self.windowPopup.setGeometry,
                             has_value_axis=self.hasValueAxis,
                             show_xscrollbar=show_xscrollbar,
                             show_yscrollbar=show_yscrollbar,
                             expand_y_axis_tick_decimals=expand_y_axis_tick_decimals)

    # TBD: not sure if this is what we want here
    self.scrolled_window.windowFrame = self # needed for WindowBasic macros

    for axisPanel in windowPane.sortedAxisPanels()[2:]:
      self.createRegionSelector(axisPanel)

    if not self.hasValueAxis:
      for slicePanel in windowPane.sortedSlicePanels()[2:]:
        self.createSliceCanvas(slicePanel)

    self.gridAll()

    self.setMenuItems()
    if not self.hasValueAxis:
      self.setSliceMenuItems()
    self.menuPeak = None
    
    scrolledWindow = self.scrolled_window
    scrolledWindow.canvasBind('<Double-1>', self.clickActivateStrip)
    #self.bind('<FocusIn>', self.focusIn)
    
    scrolledWindow.translateBind(button=OSButton(2))
    scrolledWindow.zoomBind(button=OSButton(2), state=shift_key_state)
    scrolledWindow.regionBind(button=OSButton(2), state=ctrl_key_state)
    #scrolledWindow.selectSingleBind(button=OSButton(1))
    scrolledWindow.selectSingleBind(button=OSButton(1), state=ctrl_key_state)
    scrolledWindow.selectMultiBind(button=OSButton(1), state=no_key_state)
    scrolledWindow.selectMultiBind(button=OSButton(1), state=shift_key_state)
    scrolledWindow.selectMultiBind(button=OSButton(1), state=shift_key_state+ctrl_key_state)
    scrolledWindow.selectMultiBind(button=OSButton(1), state=ctrl_key_state)

    # NOTE:ED - change the menu binding for MacOS
    scrolledWindow.menuBind(button=OSButton(3), menu_items=self.menu_items,
                            update_func=self.updateMenuState)

    if useWheelMouse():
      self.windowPopup.bind('<MouseWheel>', self.windowsZoom)
      #self.windowPopup.bind('<KeyPress>', self.keypress)
    else:
      scrolledWindow.canvasBind('<Button-4>', self.zoomIn)
      scrolledWindow.canvasBind('<Button-5>', self.zoomOut)

    if not self.hasValueAxis:
      scrolledWindow.sliceMenuBind(button=OSButton(3), menu_items=self.slice_menu_items,
                                   update_func=self.updateSliceMenuState)

    self.topPopup.initSpectrumWindowPane(windowPane, self)
    self.curateNotifiers(self.windowPopup.registerNotify)

  # overrides WindowDraw
  def lift(self):

    self.windowPopup.lift()

  def updateAll(self):

    self.setHandlerBackground()
    self.drawAllAfter()

  def togglePeaks(self, winPeakList):

    if winPeakList.analysisPeakList.peakList:
      if self.windowPopup.peakListCheck.get():
        setWindowPeakListState(winPeakList, not winPeakList.isSymbolDrawn)

      else:
        setWindowPeakListState(winPeakList, not winPeakList.isSymbolDrawn)

        winPeakLists = self.windowPopup.peaksToggleSelector.objects
        for winPeakList1 in winPeakLists:
          if winPeakList1 is not winPeakList:
            setWindowPeakListState(winPeakList1, 0)

  def toggleStrips(self, object):

    stripAxis = self.windowPane.spectrumWindow.stripAxis
    axisRegions = self.windowPane.findFirstAxisPanel(label=stripAxis).sortedAxisRegions()

    if object == 'all':
      for axisRegion in axisRegions:
        axisRegion.isActive = True
    else:
      for axisRegion in axisRegions:
        axisRegion.isActive = False
      axisRegion = object
      axisRegion.isActive = True

      num = axisRegions.index(axisRegion)
      self.updateRegionSelectors(num)

  def getActiveStrip(self):

    stripAxis = self.windowPane.spectrumWindow.stripAxis
    axisRegions = self.windowPane.findFirstAxisPanel(label=stripAxis).sortedAxisRegions()
    for n in range(len(axisRegions)):
      axisRegion = axisRegions[n]
      if (axisRegion.isActive):
        if (n == 0):
          break
        else:
          return n
    else:
      return None

    for m in range(1, len(axisRegions)):
      axisRegion = axisRegions[m]
      if (not axisRegion.isActive):
        return n

    return 'all'

  def clickActivateStrip(self, event):
  
    scrolledWindow = self.scrolled_window
  
    canvas = event.widget
  
    (row, col) = scrolledWindow.getCanvasRowCol(canvas)

    if self.windowPane.spectrumWindow.stripAxis == 'x':
      currentStrip = col
    else:
      currentStrip = row
 
    self.activateStrip(currentStrip)

  def activateStrip(self, num):

    stripAxis = self.windowPane.spectrumWindow.stripAxis
    axisRegions = self.windowPane.findFirstAxisPanel(label=stripAxis).sortedAxisRegions()
    N = len(axisRegions)
    if num >= N:
      num = 0

    for i in range(N):
      axisRegion = axisRegions[i]
      if i == num:
        axisRegion.isActive = True
      else:
        axisRegion.isActive = False

    self.updateRegionSelectors(num)

  def updateRegionSelectors(self, num):

    for axisPanel in self.windowPane.sortedAxisPanels()[2:]:
      if axisPanel.axisType.isSampled:
        self.setPseudoState(axisPanel)
      else:
        axisRegion = axisPanel.sortedAxisRegions()[num]
        (r0, r1) = Util.checkSwapRegion(axisRegion.region, axisPanel.axisUnit)
        axisPanel.region_selector.setViewRegion(r0, r1, do_callback=True)

  def unpostMenu(self, *event):

    self.scrolled_window.menu.unpost()
    if not self.hasValueAxis:
      self.scrolled_window.sliceMenu.unpost()

  def moveStripToWindow(self, sourceNum, windowPane=None, targetNum=0):

    label = self.windowPane.spectrumWindow.stripAxis
    if label == 'x':
      doCols = True
    else:
      doCols = False

    # TBD: get rid of shuffleAxisRegions
    if windowPane is self.windowPane:
      self.shuffleAxisRegions(sourceNum=sourceNum, targetNum=targetNum, doCols=doCols)
      return

    strip = self.windowPane.findFirstAxisPanel(label=label).sortedAxisRegions()[sourceNum]

    if windowPane:
      windowFrame = windowPane.getWindowFrame()

      n = len(windowPane.findFirstAxisPanel(label=label).axisRegions)
      if targetNum >= n:
        targetNum = n

      if doCols:
        windowFrame.addCol(col=targetNum)
      else:
        windowFrame.addRow(row=targetNum)

    else:
      analysisProject = self.topPopup.analysisProject
      i = len(analysisProject.spectrumWindows)
      name = 'window%s' % i
      while analysisProject.findFirstSpectrumWindow(name = name):
        i += 1
        name = 'window%s' % i

      axisTypes = [ap.axisType for ap in self.windowPane.sortedAxisPanels()]

      window = createSpectrumWindow(self.windowPopup.project, name, [axisTypes,], ncols=1, nrows=1)
      windowPane = window.findFirstSpectrumWindowPane()

    if windowPane is self.windowPane:
      if targetNum < sourceNum:
        sourceNum += 1

    sortedAxisPanelsS = self.windowPane.sortedAxisPanels()
    sortedAxisPanelsT = windowPane.sortedAxisPanels()
    for i in range(len(sortedAxisPanelsS)):
      if sourceNum < len(sortedAxisPanelsS[i].axisRegions):
        if targetNum < len(sortedAxisPanelsT[i].axisRegions):
          sourceStrip = sortedAxisPanelsS[i].sortedAxisRegions()[sourceNum]
          targetStrip = sortedAxisPanelsT[i].sortedAxisRegions()[targetNum]
          targetStrip.region = tuple(sourceStrip.region)
          targetStrip.size   = sourceStrip.size

    if doCols:
      self.deleteCol(col=sourceNum)
    else:
      self.deleteRow(row=sourceNum)

  def getMoveStripItems(self, num):

    items = [{ 'kind': 'command', 'label': 'New window', 'shortcut': 'N',
               'command': lambda event, s=num:self.moveStripToWindow(s) }]
    windows = getEquivalentWindows(self.windowPane.spectrumWindow)

    stripAxis = self.windowPane.spectrumWindow.stripAxis
    if stripAxis == 'x':
      doCols = True
    else:
      doCols = False

    for window in windows:
      for windowPane in window.sortedSpectrumWindowPanes():
        axisRegions = windowPane.findFirstAxisPanel(label=stripAxis).sortedAxisRegions()

        for i in range(len(axisRegions)+1):

          if (window is self.windowPane) and ( (i == num) or (i == len(axisRegions)) ):
            continue

          cmd = lambda event, w=windowPane, s=num, t=i: self.moveStripToWindow(s,w,t)
          label = '%s strip %d' % (getWindowPaneName(windowPane),i+1)
          item = { 'kind': 'command', 'label': label, 'shortcut': '%d' % (i+1),
                   'command': cmd}
          items.append(item)


    return items

  # TBD: dangerous way this done here, using explicit indices into submenus
  def updateMenuState(self):

    scrolled_window = self.scrolled_window
    windowMenu = scrolled_window.menu
    windowSubmenuDict = windowMenu.submenuDict
    event  = windowMenu.menu_event
    canvas = event.widget
    topPopup = self.topPopup
    windowPane = self.windowPane
    currentPeaks = topPopup.currentPeaks
    
    (x, y, propX, propY) = self.calcWorldCoord(canvas, event.x, event.y)
    (row, col) = scrolled_window.getCanvasRowCol(canvas)
    #print 'updateMenuState', row, col
    xAxisRegion = windowPane.findFirstAxisPanel(label='x').sortedAxisRegions()[col]
    yAxisRegion = windowPane.findFirstAxisPanel(label='y').sortedAxisRegions()[row]

    # View option
    if not self.hasValueAxis:
      view_menu = windowSubmenuDict['View']

      if self.getSpectrumViews():
        state = Tkinter.NORMAL
      else:
        state = Tkinter.DISABLED

      view_menu.entryconfig(2, state=state)

    stripAxis = self.windowPane.spectrumWindow.stripAxis
    if stripAxis == 'x':
      stripNum  = col
      otherAxis = 'y'
    else:
      otherAxis = 'x'
      stripNum  = row

    # Strip option
    cursorPosition = self.findPosition(x,y,stripNum,returnDict=False)
    allPositions   = [cursorPosition]

    if not currentPeaks:
      state = Tkinter.DISABLED
    else:
      state = Tkinter.NORMAL
      for peak in currentPeaks:
        spectrum = peak.peakList.dataSource
        if isSpectrumInWindowPane(windowPane, spectrum):
          position = getPeakWindowPosition(peak, windowPane, default=cursorPosition)
          allPositions.append( position )

    # Setup orthogonal planes for strip and general navigation
    window = windowPane.spectrumWindow
    transposes     = findWindowTransposes(windowPane, cursorPosition)
    windowZplanes  = findOrthogonalWindows(windowPane, allPositions, minDims=2)
    peak_strip_items = [ { 'kind': 'command', 'label': 'Current Window', 'shortcut': 'C',
                           'command': self.showPeakStrips, 'state': state} ]

    if windowZplanes and currentPeaks:
      for (planeName, pane, positions) in windowZplanes:
        cmd = lambda event, w=pane, p=positions[1:]: self.makeOrthogonalStrips(w, p)
        peak_strip_items.append({'kind': 'command', 'label': planeName, 'state': state,
                                'command': cmd} )

    strip_items = []

    if len(windowPane.findFirstAxisPanel(label=stripAxis).axisRegions) > 1:
      move_strip_items = self.getMoveStripItems(stripNum)

      if stripAxis == 'x':
         strip_items = [ \
           { 'kind': 'command', 'label': 'Add New Strip', 'shortcut': 'N',
             'tipText': 'Add a new vertical strip sub-division to the window',
             'command': lambda event: self.addCol(col=stripNum+1)  },
           { 'kind': 'command', 'label': 'Delete Strip', 'shortcut': 'D',
             'tipText': 'Delete the strip currently under the cursor',
             'command': lambda event: self.deleteCol(col=stripNum) },
           { 'kind': 'command', 'label': 'Make Active', 'shortcut': 'A',
             'tipText': 'Make the strip under the cursor the active one; for Z-axis navigation, moving, deletion etc.',
             'command': lambda event: self.activateStrip(stripNum) },
           { 'kind': 'command', 'label': 'Move Left', 'shortcut': 'L',
             'tipText': 'Shuffle the strip under the cursor one position to the left',
             'command': lambda event: self.moveStripsLeft(stripNum) },
           { 'kind': 'command', 'label': 'Move Right', 'shortcut': 'R',
             'tipText': 'Shuffle the strip under the cursor one position to the right',
             'command': lambda event: self.moveStripsRight(stripNum) },
           { 'kind': 'cascade', 'label': 'Peak Location Strips', 'shortcut': 'P',
             'tipText': 'Make vertical strips based on the selected peak positions in a window with matching axes',
             'submenu': peak_strip_items, 'state': state },
           { 'kind': 'cascade', 'label': 'Change window', 'shortcut': 'C',
             'tipText': 'Move the vertical strip under the cursor to a different window with matching axes',
             'submenu': move_strip_items },
           { 'kind': 'command', 'label': 'Add horizontal separator', 'shortcut': 'h',
             'tipText': 'Split the Y-axis of the window to give a separate region; X- and Z-axes are still coordinated, obeying any strips',
             'command': lambda event: self.addRow() },
           { 'kind': 'command', 'label': 'Switch to horizontal', 'shortcut': 'S',
             'tipText': 'Change from vertical window strips to horizontal strips (some positional information may be lost)',
             'command': lambda event: self.windowPopup.toggleStripDir() },
         ]

      else:
        strip_items = [ \
          { 'kind': 'command', 'label': 'Add New Strip', 'shortcut': 'N',
            'tipText': 'Add a new horizontal strip sub-division to the window',
            'command': lambda event: self.addRow(row=stripNum+1)  },
          { 'kind': 'command', 'label': 'Delete Strip', 'shortcut': 'D',
            'tipText': 'Delete the strip currently under the curso',
            'command': lambda event: self.deleteRow(row=stripNum) },
          { 'kind': 'command', 'label': 'Make Active', 'shortcut': 'A',
            'tipText': 'Make the strip under the cursor the active one; for Z-axis navigation, moving, deletion etc.',
            'command': lambda event: self.activateStrip(stripNum) },
          { 'kind': 'command', 'label': 'Move Up', 'shortcut': 'U',
            'tipText': 'Shuffle the strip under the cursor one position upwards',
            'command': lambda event: self.moveStripsLeft(stripNum) },
          { 'kind': 'command', 'label': 'Move Down', 'shortcut': 'D',
            'tipText': 'Shuffle the strip under the cursor one position downwards',
            'command': lambda event: self.moveStripsRight(stripNum) },
          { 'kind': 'cascade', 'label': 'Peak Location Strips', 'shortcut': 'P',
            'tipText': 'Make horizontal strips based on the selected peak positions in a window with matching axes',
            'submenu': peak_strip_items, 'state': state },
          { 'kind': 'cascade', 'label': 'Change window', 'shortcut': 'C',
            'tipText': 'Move the horizontal strip under the cursor to a different window with matching axes',
            'submenu': move_strip_items },
          { 'kind': 'command', 'label': 'Add vertical separator', 'shortcut': 'v',
            'tipText': 'Split the X-axis of the window to give a separate region; Y- and Z-axes are still coordinated, obeying any strips',
            'command': lambda event: self.addCol() },
          { 'kind': 'command', 'label': 'Switch to vertical', 'shortcut': 'S',
            'tipText': 'Change from horizontal window strips to vertical strips (some positional information may be lost)',
            'command': lambda event: self.windowPopup.toggleStripDir() },
        ]


    else:
      if stripAxis == 'x':
          strip_items = [ \
           { 'kind': 'command', 'label': 'Add new strip',
             'tipText': 'Add a new vertical strip sub-division to the window',
             'command': lambda event: self.addCol(), 'shortcut': 'A'  },
           { 'kind': 'cascade', 'label': 'Peak Location Strips', 'shortcut': 'P',
             'tipText': 'Make vertical strips based on the selected peak positions in a window with matching axes',
             'submenu': peak_strip_items, 'state': state },
           { 'kind': 'command', 'label': 'Add horizontal separator', 'shortcut': 'h',
             'tipText': 'Split the Y-axis of the window to give a separate region; X- and Z-axes are still coordinated',
             'command':  lambda event: self.addRow()},
           { 'kind': 'command', 'label': 'Switch to horizontal', 'shortcut': 'S',
             'tipText': 'Switch from adding vertical window strips to adding horizontal strips',
             'command': lambda event: self.windowPopup.toggleStripDir() },
         ]
      else:
         strip_items = [ \
           { 'kind': 'command', 'label': 'Add new strip',
             'tipText': 'Add a new horizontal strip sub-division to the window',
             'command': lambda event: self.addRow(), 'shortcut': 'A'  },
           { 'kind': 'cascade', 'label': 'Peak Location Strips', 'shortcut': 'P',
             'tipText': 'Make horizontal strips based on the selected peak positions in a window with matching axes',
             'submenu': peak_strip_items, 'state': state },
           { 'kind': 'command', 'label': 'Add vertical separator', 'shortcut': 'v',
             'tipText': 'Split the X-axis of the window to give a separate region; Y- and Z-axes are still coordinated',
             'command': lambda event: self.addCol()},
           { 'kind': 'command', 'label': 'Switch to vertical', 'shortcut': 'S',
             'tipText': 'Switch from adding horizontal window strips to adding vertical strips',
             'command': lambda event: self.windowPopup.toggleStripDir() },
         ]

    if len(windowPane.findFirstAxisPanel(label=otherAxis).axisRegions) > 1:
      if stripAxis == 'x':
        strip_items.append( { 'kind': 'command', 'label': 'Delete horizontal separator',
            'tipText': 'Remove the Y-axis sub-division currently under the cursor; does not effect any vertical strips',
            'command': lambda event: self.deleteRow(row=row), 'shortcut': 'o' } )
        strip_items.append( { 'kind': 'command', 'label': 'Delete all separators',
            'tipText': 'Remove all Y-axis window sub-divisions; does not effect any vertical strips',
            'command': lambda event: self.deleteSeparators(), 'shortcut': 'a' } )
      else:
        strip_items.append( { 'kind': 'command', 'label': 'Delete vertical separator',
            'tipText': 'Remove the X-axis sub-division currently under the cursor; does not effect any horizontal strips',
            'command': lambda event: self.deleteCol(col=col), 'shortcut': 'o' } )
        strip_items.append( { 'kind': 'command', 'label': 'Delete all separators',
            'tipText': 'Remove all X-axis window sub-divisions; does not effect any horizontal strips ',
            'command': lambda event: self.deleteSeparators(), 'shortcut': 'a' } )
      

    stripMenu = windowSubmenuDict['Strip']
    stripMenu.setMenuItems(strip_items)

    # peak and assign options
    position_region = self.findPositionRegion(x, y, stripNum)

    peak = self.findNearestPeak(position_region, xAxisRegion, yAxisRegion)
    find_items  = []
    seqss_items = []
    propsister_items = []
    if peak:
      label = 'Assign %s:%s:%s:%s' % (peak.peakList.dataSource.experiment.name,
                                      peak.peakList.dataSource.name,
                                      peak.peakList.serial,
                                      peak.serial)
      state = Tkinter.NORMAL
      i = 1
      for peakDim in peak.sortedPeakDims():
        if not peakDim.dataDimRef:
          continue

        expDimRef = peakDim.dataDimRef.expDimRef
        isotopesCode = ','.join( [ic for ic in expDimRef.isotopeCodes]  )

        spec_find_items = []
        for experiment in peakDim.topObject.sortedExperiments():
          for spectrum in experiment.sortedDataSources():
            isotopes = getSpectrumIsotopes(spectrum)
            for isotope in expDimRef.isotopeCodes:
              if isotope in isotopes:
                specText = '%s:%s' % (experiment.name,spectrum.name)
                spec_find_items.append({'kind': 'command', 'label':specText,
                                        'tipText': 'Search for peaks with a similar dimension %d position in this spectrum' % i,
                                        'command': lambda event, pd=peakDim, sp=spectrum: self.findClosePeaks(pd,sp)})
                break

        if len(spec_find_items) > 1:
          spec_find_items.append({'kind': 'command', 'label':'All spectra', 'shortcut': 'A',
                                  'tipText': 'Find any peaks close to the selected position in all spectra listed above',
                                  'command': lambda event, pd=peakDim: self.findClosePeaks(pd)})

        text = 'F%d: %s %7.3f %s' % (i,isotopesCode,getPeakDimPosition(peakDim,toUnit=expDimRef.unit),expDimRef.unit)
        find_items.append({'kind': 'cascade', 'label':text,
                           'tipText': 'Find other peaks that have a similar dimension %d position to the peak at the cursor' % i,
                           'shortcut': '%d' % i, 'submenu': spec_find_items})
        i += 1

      seqOffsets = []
      N = len(peak.peakDims)
      if N == 2:
        seqOffsets = [(None, -1),(None,1)]

      elif N == 3:
        seqOffsets = [(None,-1,None),(None,1,None),
                      (None,None,-1),(None,None,1),
                      (-1,None,None),(1,None,None)]

      elif N == 4:
        seqOffsets = [(-1,-1,None,None),(1,1,None,None),
                      (None,-1,-1,None),(None,1,1,None)]

      elif N == 5:
        seqOffsets = [(-1,-1,None,None,None),(1,1,None,None,None),
                      (-1,-1,-1,None,None),  (1,1,1,None,None)]

      for seqOffset in seqOffsets:
        text = []
        tips = []
        f = 1
        for j in seqOffset:
          if j is None:
            text.append('F%d 0' % (f) )
            tips.append('the dimension %d assignment to have no sequence offset' % (f))
          elif j > 0:
            text.append('F%d +%d' % (f,j) )
            tips.append('the dimension %d assignment to be %d sequence position later' % (f,j))
          else:
            text.append('F%d %d' % (f,j) )
            tips.append('the dimension %d assignment to be %d sequence position earlier' % (f,j))

          f += 1
        seqss_items.append({'kind': 'command', 'label':','.join(text),
                            'tipText': 'Set ' + ','.join(tips),
                            'command': lambda event, s=seqOffset: self.makeSeqSpinSystemConnections(s)})

      # propagate assignments to sister peak lists
      peakList0 = peak.peakList
      for peakList in peakList0.dataSource.sortedPeakLists():
        if peakList is not peakList0:
          label2 = '%s:%s:%s' % (peakList.dataSource.experiment.name,
                                 peakList.dataSource.name,
                                 peakList.serial)
          tipText = 'Propagate assignments to %s' % label2
          propsister_items.append({'kind': 'command', 'label':label2,
                            'tipText': tipText,
                            'command': lambda event, pl=peakList: self.propagatePeakAssignmentsToSister(pl)})
          
    else:
      label = 'Assign peak'
      state = Tkinter.DISABLED

    self.menuPeak = peak
    assn_menu = windowSubmenuDict['Assign']
    assn_menu.entryconfig(0, state=state, label=label)
    assn_menu.entryconfig(3, state=state)
    assn_menu.entryconfig(6, state=state)
    assn_menu.entryconfig(7, state=state)
    if currentPeaks:
      assn_menu.entryconfig(8, state=Tkinter.NORMAL)
    else:
      assn_menu.entryconfig(8, state=state)
    seqss_menu = assn_menu.submenuDict['Set sequential spin systems']
    seqss_menu.setMenuItems(seqss_items)
    propsister_menu = assn_menu.submenuDict['Propagate to sister peak lists']
    propsister_menu.setMenuItems(propsister_items)

    peak_menu = windowSubmenuDict['Peak']
    peak_menu.entryconfig(2, state=Tkinter.NORMAL)
    for n in (3, 6, 12, 13):
      peak_menu.entryconfig(n, state=state)
    if currentPeaks: # Fit peaks item
      peak_menu.entryconfig(4, state=Tkinter.NORMAL)
      peak_menu.entryconfig(14, state=Tkinter.NORMAL)
      peak_menu.entryconfig(15, state=Tkinter.NORMAL)
    else:
      peak_menu.entryconfig(4, state=state)
      peak_menu.entryconfig(14, state=state)
      peak_menu.entryconfig(15, state=state)

    peak_nav_menu = windowSubmenuDict['Locate Peaks']
    peak_nav_menu.entryconfig(2, state=Tkinter.NORMAL)
    for n in (1, 2, 3, 6):
      peak_nav_menu.entryconfig(n, state=state)

    find_menu = peak_nav_menu.submenuDict['Find close peaks']
    find_menu.setMenuItems(find_items)

    state = Tkinter.DISABLED
    if peak:
      for peakDim in peak.peakDims:
        for contrib in peakDim.peakDimContribs:
          if hasattr(contrib, 'resonance'):
            if contrib.resonance.resonanceGroup:
              state = Tkinter.NORMAL
              break
 
          else:
            for resonance in contrib.resonances:
              if contrib.resonance.resonanceGroup:
                state = Tkinter.NORMAL
                break
            else:
              continue
            break  
               
        else:
          continue
        break

    assn_menu.entryconfig(3, state=state)

    state = Tkinter.DISABLED
    match_items = []
    peakLists   = self.getActivePeakLists()
    if peakLists and currentPeaks:
      state = Tkinter.NORMAL
      peak = currentPeaks[0]

      for peakList in peakLists:
        spectrum = peakList.dataSource
        specName = '%s:%s' % (spectrum.experiment.name,spectrum.name)
        specId = (specName,peakList.serial)
        dim_items = []
        # TBD: Automate this to give only options for non-aligned peakDims
        for i in range(spectrum.numDim):
          text = 'F%d' % (i+1)
          dim_items.append({'kind': 'command', 'label':text, 'shortcut': '%d' % (i+1),
                            'tipText': 'Match selected peaks to the stated spectrum using dimension %d' % (i+1),
                            'command': lambda event, p=currentPeaks, i=i,
                             pl=peakList: self.findPositionMatches(p,i,pl)})

        match_items.append({'kind': 'cascade', 'label':'In %s:%d' % specId,
                            'tipText': 'Match the selected peaks to those in spectrum %s, peak list %d ' % specId,
                            'submenu': dim_items})

      dim_items = []
      for i in range(len(peak.peakDims)):
        text = 'F%d' % (i+1)
        dim_items.append({'kind': 'command', 'label':text,
                          'tipText': 'Match selected peaks to all spectra displayed in the window using dimension %d' % (i+1),
                          'command': lambda event, p=currentPeaks, i=i: self.findPositionMatches(p,i)})
      match_items.append({'kind': 'cascade', 'label':'In active spectra',
                          'tipText': 'Match the selected peaks to all spectra displayed in the window',
                          'submenu': dim_items, 'shortcut': 'I' })

    fit_items = []
    if len(currentPeaks) > 1:
      for fitType in PEAK_FIT_TYPES:
        for fitMethod in PEAK_FIT_METHODS[1:]:
          fit_items.append({'kind': 'command', 'label':'%s, %s' % (fitMethod, fitType),
            'tipText': 'Fit the selected peaks %s using the %s method' % (fitType, fitMethod),
            'command': lambda event, fitMethod=fitMethod, fitType=fitType: self.fitPeaks(fitMethod, fitType), 'shortcut': fitMethod[0]})
    elif len(currentPeaks) == 1:
      for fitMethod in PEAK_FIT_METHODS[1:]:
        fit_items.append({'kind': 'command', 'label':'%s' % fitMethod,
        'command': lambda event, fitMethod=fitMethod: self.fitPeaks(fitMethod), 'shortcut': fitMethod[0]})

    fit_menu = peak_menu.submenuDict['Fit selected']
    fit_menu.setMenuItems(fit_items)
    
    move_items = []
    if currentPeaks:
      peak = currentPeaks[0]
      peakList = peak.peakList
      dataSource = peakList.dataSource
      peakLists = [pl for pl in dataSource.sortedPeakLists() if pl is not peakList]
      for peakList in peakLists:
        label = str(peakList.serial)
        details = peakList.details
        if details:
          nmax = 30
          label += ' (%s' % details[:nmax]
          if len(details) > nmax:
            label += '...'
          label += ')'
        shortcut = str(peakList.serial)[0]
        move_items.append({'kind': 'command', 'label':label,
          'command': lambda event, peakList=peakList: self.movePeaksToPeakList(peakList), 'shortcut': shortcut})

    move_menu = peak_menu.submenuDict['Move to peak list']
    move_menu.setMenuItems(move_items)
    
    match_menu = peak_nav_menu.submenuDict['Match multiple peaks']
    match_menu.setMenuItems(match_items)

    assn_menu.entryconfig(1, state=state)
    assn_menu.entryconfig(2, state=state)
    assn_menu.entryconfig(4, state=state)
    for n in (0, 2, 5, 9, 10, 11):
      peak_menu.entryconfig(n, state=state)
      
    for n in (0, 4, 5):
      peak_nav_menu.entryconfig(n, state=state)
    
    state = Tkinter.DISABLED
    if peak and currentPeaks:
      state = Tkinter.NORMAL
    peak_menu.entryconfig(7, state=state)
    peak_menu.entryconfig(13, state=state)

    # Navigate option
    goto_items = []

    if windowZplanes:
      for (planeName, pane, positions) in windowZplanes:
        cmd = lambda event, w=pane, position=positions[0]: \
                      self.gotoOrthogonalPlane(w, position)
        goto_items.append({'kind': 'command', 'label': planeName,
                           'tipText': 'Navigate to the equivalent position of the cursor in a different window, mapping the current axes to the target axes',
                           'command': cmd} )

      state = Tkinter.DISABLED
      if currentPeaks:
        state = Tkinter.NORMAL
      goto_items.append({'kind': 'cascade', 'label': 'Peak Location Strips', 'shortcut': 'P',
                         'tipText': 'Make strips based on the selected peak positions in a window with matching axes',
                         'submenu': peak_strip_items, 'state': state })

      # Below is just too confusing
      #if len(windowZplanes) > 1:
      #  cmd = lambda event, wzps=windowZplanes: self.gotoOrthogonalAllPlanes(wzps)
      #  goto_items.append({'kind': 'command','label': 'All window positions', 
      #                     'tipText': '',
      #                     'shortcut': 'A','command': cmd} )

    if transposes:
      if windowZplanes:
        goto_items.append({ 'kind': 'separator'})

      if len(transposes) > 1:
        transpose_items = []
        for label, position in transposes:
          cmd = lambda event, pos=position, col=col, row=row: \
                        self.gotoPosition(pos, row, col)
          transpose_items.append({'kind': 'command', 'label': label,
                                  'tipText': 'Navigate to the stated symmetry related homonuclear position, e.g. swapping two 13C or 1H locations',
                                  'command': cmd })
        goto_items.append({'kind': 'cascade', 'label':'Transpose',
                           'tipText': 'Navigate to symmetry related positions',
                           'submenu': transpose_items, 'shortcut': 'T' })

      else:
        label, position = transposes[0]
        cmd = lambda event, pos=position, col=col, row=row: \
                      self.gotoPosition(pos, row, col)
        goto_items.append({'kind': 'command', 'label': 'Transpose '+label, 'shortcut': 'T',
                           'tipText': 'Navigate to the stated symmetry related homonuclear position, e.g. swapping two 13C or 1H locations',
                           'command': cmd})


    goto_menu = windowSubmenuDict['Navigate']
    goto_menu.setMenuItems(goto_items)


    # mark & rulers options
    mark_menu = windowSubmenuDict['Markers']
    analysisProject = window.analysisProject

    state = Tkinter.DISABLED
    if analysisProject.marks:
      state = Tkinter.NORMAL
    mark_menu.entryconfig(3, state=state) # marks

    state = Tkinter.DISABLED
    if analysisProject.rulers:
      state = Tkinter.NORMAL
    mark_menu.entryconfig(4, state=state) # rulers


    # macros options
    argServer   = topPopup.argumentServer
    macro_menu  = windowSubmenuDict['Macros']
    macro_items = [{'kind': 'command', 'label': 'Reload macros',
                    'tipText': 'Reload the Python macro scripts currently listed in the window menu',
                    'command': self.reloadMouseMacros, 'shortcut': 'R' }]
    for macro in self.windowPopup.analysisProfile.macros:
      if macro.isInMouseMenu:
        macro_items.append( {'kind': 'command', 'label': macro.name,
                             'tipText': 'Launch the macro called "%s"; any positional information comes from the cursor location' % macro.name,
                             'command': lambda event, m=macro: Util.runMacro(m,argServer)} )
    macro_menu.setMenuItems(macro_items)

    state = Tkinter.DISABLED
    if len(macro_items) > 1:
      state = Tkinter.NORMAL
    macro_menu.entryconfig(0, state=state)

    self.setCurrentObjects(event)

    topPopup.currentPeakLists = []
    views = self.getActiveSpectrumViews()
    for view in views:
      for winPeakList in view.windowPeakLists:
        if winPeakList.isSymbolDrawn:
          try:
            peakList = winPeakList.analysisPeakList.peakList
          except:
            peakList = None
          if (peakList and not peakList.isDeleted):
            topPopup.currentPeakLists.append(peakList)

    topPopup.currentPosition = (x,y)

  def setCurrentObjects(self, event):

    canvas = event.widget
    (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
    views = self.getActiveSpectrumViews()

    topPopup = self.topPopup
    topPopup.currentPeak   = self.menuPeak
    topPopup.currentWindow = self.windowPane.spectrumWindow
    topPopup.currentWindowPane = self.windowPane

    if self.windowPane.spectrumWindow.stripAxis == 'x':
      topPopup.currentStrip = col
    else:
      topPopup.currentStrip = row

    topPopup.currentWindowPopup  = self.windowPopup
    topPopup.currentWindowFrame  = self
    topPopup.currentSpectra   = [view.analysisSpectrum.dataSource for view in views]
    topPopup.currentCanvas    = canvas
    topPopup.currentEvent     = event

  def reloadMouseMacros(self, *event):

    profile = self.windowPane.root.currentAnalysisProfile
    for macro in profile.macros:
      if macro.isInMouseMenu:
        Util.reloadMacro(macro, self.topPopup.argumentServer)

  # TBD: more general
  def setMenuItems(self):


    zoom_items = [ \
      { 'kind': 'command', 'label': 'Out by 2',
        'tipText': 'Change the viewed region of the contour display to see twice as much on the X and Y axes',
        'command': lambda event: self.zoom(event.widget, 2.0), 'shortcut': 'O' },
      { 'kind': 'command', 'label': 'In by 2',
        'tipText': 'Change the viewed region of the contour display to see half as much on the X and Y axes',
        'command': lambda event: self.zoom(event.widget, 0.5), 'shortcut': 'I' },
      { 'kind': 'command', 'label': 'Zoom to selected peaks',
        'tipText': 'Center the viewed region of the contour display to cover any peaks selected in the window',
        'command': self.zoomToSelectedPeaks, 'shortcut': 'Z' },
      { 'kind': 'command', 'label': 'Zoom to spin system peaks', 'shortcut': 's',
        'tipText': 'Change the window view to see all peaks assigned to the same, primary spin system as the peak under the cursor',
        'command': self.zoomToSpinSystem },
    ]


    mark_items = [ \
      { 'kind': 'command', 'label': 'Add Mark',
        'tipText': 'Add a multi-dimensional cross mark at the current cursor location',
        'command': self.createMark, 'shortcut': 'M' },
      { 'kind': 'command', 'label': 'Add Horizontal ruler',
        'tipText': 'Add a single horizontal line at the current cursor location',
        'command': self.createHorizontalRuler, 'shortcut': 'H' },
      { 'kind': 'command', 'label': 'Add Vertical ruler',
        'tipText': 'Add a single vertical line at the current cursor location',
        'command': self.createVerticalRuler, 'shortcut': 'V' },
      { 'kind': 'command', 'label': 'Delete marks',
        'tipText': 'Delete all multi-dimensional cross marks in all windows',
        'command': self.clearMarks, 'shortcut': 'k' },
      { 'kind': 'command', 'label': 'Delete rulers',
        'tipText': 'Delete all horizontal and vertical ruler (1D) lines in all windows',
        'command': self.clearRulers, 'shortcut': 'l' },
      { 'kind': 'separator'},
      { 'kind': 'command', 'label': 'Options',
        'tipText': 'Open a popup window to alter the display settings for cross marks and ruler lines',
        'command': lambda event: self.topPopup.editMarks(), 'shortcut': 'O' },
    ]


    strip_items = [ \
      { 'kind': 'command', 'label': 'Add Vertical strip',
        'tipText': 'Add a new vertical strip/sub-division to the window',
        'command': lambda: self.addCol(), 'shortcut': 'V' },
      { 'kind': 'command', 'label': 'Add Horizontal strip',
        'tipText': 'Add a new Horizontal strip/sub-division to the window',
        'command': lambda event: self.addRow(), 'shortcut': 'H'},
      { 'kind': 'command', 'label': 'Delete vertical strip',
        'tipText': 'Delete the vertical strip currently under the cursor',
        'command': lambda event: self.deleteCol(), 'shortcut': 'v'},
      { 'kind': 'command', 'label': 'Delete horizontal strip',
        'tipText': 'Delete the horizontal strip currently under the cursor',
        'command': lambda event: self.deleteRow(), 'shortcut': 'h'},
      { 'kind': 'command', 'label': 'Switch direction', 'shortcut': 'S',
        'tipText': 'Swap between vertical window strips and horizontal strips',
        'command': lambda event: self.windowPopup.toggleStripDir() },
    ]

    contour_items = []
    macro_items = [] # created properly when menu pops up
    find_items  = [] # created properly when menu pops up
    match_items = [] # created properly when menu pops up
    seqss_items = [] # created properly when menu pops up
    propsister_items = [] # created properly when menu pops up
    merit_items = []
    merit_items.append({'kind': 'command', 'label':'0.0',
        'tipText': 'Set the figure-of-merit value for the peak under the cursor to 0.0; a "bad" peak',
        'command': lambda event, m=0.0: self.setPeakMerit(m), 'shortcut': '0'})
    merit_items.append({'kind': 'command', 'label':'0.5',
        'tipText': 'Set the figure-of-merit value for the peak under the cursor to 0.5; a "mediocre" peak',
        'command': lambda event, m=0.5: self.setPeakMerit(m), 'shortcut': '5'})
    merit_items.append({'kind': 'command', 'label':'1.0',
        'tipText': 'Set the figure-of-merit value for the peak under the cursor to 1.0; a "good" peak',
        'command': lambda event, m=1.0: self.setPeakMerit(m), 'shortcut': '1'})

    assn_items = [ \
      { 'kind': 'command', 'label': 'Assign peak',
        'tipText': 'Open a popup window to assign the dimensions of the peak under the cursor to resonances',
        'command': self.assignPeak, 'shortcut': 'A' },
      { 'kind': 'command', 'label': 'Propagate assignments',
        'tipText': 'Spread resonance assignments from the peak under the cursor to all selected peaks; where shifts match peak dimension positions',
        'command': self.propagatePeakAssignments, 'shortcut': 'P' },
      { 'kind': 'cascade', 'label': 'Propagate to sister peak lists',
        'tipText': 'Spread resonance assignments from the peak under the cursor to peaks in sister peak lists; where shifts match peak dimension positions',
        'submenu': propsister_items, 'shortcut': 'g' },
      { 'kind': 'command', 'label': 'Add to spin system',
        'tipText': 'For the peak under the cursor, put all resonance assignments into the same spin system (resonance group)',
        'command': self.addPeakSpinSystem, 'shortcut': 's' },
      { 'kind': 'command', 'label': 'Predict Spin System Type',
        'tipText': 'Predict the residue type of the primary spin system the peak under the cursor is assigned to',
        'command': self.predictSpinSystemType, 'shortcut': 'T' },
      { 'kind': 'command', 'label': 'Deassign selected peaks',
        'tipText': 'Clear all peak dimension assignments (to resonances) for all selected peaks (potentially in many windows)',
        'command': self.clearPeakAssignments, 'shortcut': 'D' },
      { 'kind': 'cascade', 'label': 'Set sequential spin systems',
        'tipText': 'Set the sequential links between resonances assigned to the dimensions of the peak under the cursor',
        'submenu': seqss_items, 'shortcut': 'q' },
      { 'kind': 'command', 'label': 'Clear all sequential links',
        'tipText': 'Remove all sequential links between resonances assigned to the dimensions of the peak under the cursor',
        'command': self.clearSeqLinks, 'shortcut': 'C' },
      { 'kind': 'command', 'label': 'Unite resonance positions',
        'tipText': 'Based on the cursor peak having representative shift values for its assigned resonances, align all peak dimensions assigned to the same resonances to the same value(s)',
        'command': self.uniteResonancePositions, 'shortcut': 'n' },
    ]

    peak_symm_items = [ \
      { 'kind': 'command', 'label': 'Use strips',
        'tipText': 'Display the results of finding symmetry related (e.g. return) peaks using strips in the window',
        'command': self.stripSymmetryRelatedPeaks, 'shortcut': 's' },
      { 'kind': 'command', 'label': 'Use tables',
        'tipText': 'Display the results of finding symmetry related (e.g. return) peaks in a table',
        'command': self.findSymmetryRelatedPeaks, 'shortcut': 't' },
    ]

    fit_items = [] # created properly when menu pops up
    move_items = []

    peak_items = [ \
      { 'kind': 'command', 'label': 'Selection Table',
        'tipText': 'Show a table of the peaks currently selected in spectrum windows',
        'command': self.showSelectedPeaks, 'shortcut': 'S' },
      { 'kind': 'command', 'label': 'Add new peak(s)',
        'tipText': 'Make a new peak at the current cursor location; makes peaks in all spectra displayed in the window',
        'command': self.createPeak, 'shortcut': 'A' },
      { 'kind': 'command', 'label': 'Delete selected',
        'tipText': 'Delete the selected peaks, which may be in more than one window',
        'command': self.deleteSelectedPeaks, 'shortcut': 'D' },
      { 'kind': 'command', 'label': 'Set details',
        'tipText': 'Set a textual comment for the peak under the cursor',
        'command': self.editPeakDetails, 'shortcut': 'd' },
      { 'kind': 'cascade', 'label': 'Set merit',
        'tipText': 'Set the figure-of-merit value for the peak under the cursor, e.g. "bad" peaks may be excluded from certain operations',
        'submenu': merit_items, 'shortcut': 'm' },
      { 'kind': 'command', 'label': 'Snap selected peaks',
        'tipText': 'Re-center the selected peaks to their nearest spectrum intensity extremum, if within the normal peak finding tolerances',
        'command': self.snapSelectedPeaks, 'shortcut': 't' },
      { 'kind': 'command', 'label': 'Unalias peak',
        'tipText': 'Move the peak under the cursor to its real ppm value by adding or removing spectrum sweep widths, to unfold or unalias',
        'command': self.unaliasPeak, 'shortcut': 'U' },
      { 'kind': 'command', 'label': 'Unaliasing propagate',
        'tipText': 'Set the selected peaks to be in the same spectrum (sweep width based) ppm range as the peak under the cursor',
        'command': self.propagatePeakUnaliasing, 'shortcut': 'g' },
      { 'kind': 'separator'},
      { 'kind': 'command', 'label': 'Label Auto Arrange',
        'tipText': 'For the selected peaks, automatically adjust the positions of the textual peak labels to avoid overlap; this operation is zoom-level dependent',
        'command': self.autoArrangePeakLabels, 'shortcut': 'L' },
      { 'kind': 'command', 'label': 'Label Position Reset',
        'tipText': 'For the selected peaks, reset all textual peak labels to their default/unmoved positions',
        'command': self.resetPeakLabels, 'shortcut': 'l' },
      { 'kind': 'command', 'label': 'Make intermediate peak',
        'tipText': 'Make a new peak in the geometric center of the selected peaks; useful for COSY etc.',
        'command': self.makeIntermediatePeak, 'shortcut': 'i' },
      { 'kind': 'command', 'label': 'Structure connections',
        'tipText': 'Show any non-onebond connections described by the resonance assignments of the cursor peak in a graphical structure display',
        'command': self.showStructConnections, 'shortcut': 'c' },
      { 'kind': 'command', 'label': 'Re-reference to this peak',
        'tipText': 'Re-reference the points to ppm relationship for the spectra of the selected peaks by aligning peak positions; cursor peak is the reference',
        'command': self.translatePeak, 'shortcut': 'k' },
      { 'kind': 'cascade', 'label': 'Fit selected',
        'tipText': 'Fit selected peaks singly or together',
        'submenu': fit_items, 'shortcut': 'F' },
      { 'kind': 'cascade', 'label': 'Move to peak list',
        'tipText': 'Move selected peaks to a sister peak list',
        'submenu': move_items, 'shortcut': 'e' },
    ]
    peak_nav_items = [ \
      { 'kind': 'cascade', 'label': 'Match multiple peaks',
        'tipText': 'For the selected peaks, find other peaks with similar positions; in a specified spectrum and dimension',
        'submenu': match_items, 'shortcut': 'p' },
      { 'kind': 'command', 'label': 'Center z planes on peak',
        'tipText': 'Move the viewed region of the window so the peak under the cursor is in the center of the screen-orthogonal (Z) dimensions',
        'command': self.peakCenterOrthogonalPlanes, 'shortcut': 'C' },
      { 'kind': 'cascade', 'label': 'Find close peaks',
        'tipText': 'Find peaks in spectra that are close to the peak under the cursor',
        'submenu': find_items, 'shortcut': 'F' },
      { 'kind': 'cascade', 'label': 'Find symmetry related',
        'tipText': 'Find peaks that are related by homonuclear symmetry to the cursor peak, e.g. NOE return peaks with swapped 1H locations',
        'submenu': peak_symm_items, 'shortcut': 'r' },
      { 'kind': 'command', 'label': 'Strip plot selected peaks',
        'tipText': 'Make strips in the current window based on the locations of the selected peaks',
        'command': self.showPeakStrips, 'shortcut': 'p'  },
      { 'kind': 'command', 'label': 'Zoom to selected peaks',
        'tipText': 'Center the viewed region of the contour display to cover any peaks selected in the window',
        'command': self.zoomToSelectedPeaks, 'shortcut': 'Z' },
      { 'kind': 'command', 'label': 'Zoom to spin system peaks',
        'tipText': 'Change the window view to see all peaks assigned to the same, primary spin system as the peak under the cursor',
        'command': self.zoomToSpinSystem, 'shortcut': 'y' },
    ]

    contour_items = [ \
      { 'kind': 'command', 'label': 'Up by 2',
        'tipText': 'Increase the contour base level ,for all displayed spectra, by a factor of two; move away from zero',
        'command': lambda event: self.changeScale(2.0), 'shortcut': 'U' },
      { 'kind': 'command', 'label': 'Down by 2',
        'tipText': 'Decrease the contour base level, for all displayed spectra, by half; move away towards zero',
        'command': lambda event: self.changeScale(0.5), 'shortcut': 'D' },
      { 'kind': 'command', 'label': 'General...',
        'tipText': 'Open a popup window to give fine control over spectrum contour levels',
        'command': lambda event: self.windowPopup.generalContours(), 'shortcut': 'G' },
    ]

    #split_items = [ \
    #  { 'kind': 'command', 'label': 'Horizontally', 'command': self.splitWindowHorizontal },
    #  { 'kind': 'command', 'label': 'Vertically', 'command': self.splitWindowVertical },
    #]

    view_items = [ \
      { 'kind': 'cascade', 'label': 'Zoom',
        'tipText': 'Change the viewed region of the spectrum contour display',
        'submenu' : zoom_items, 'shortcut': 'Z' },
      { 'kind': 'command', 'label': 'Center here',
        'tipText': 'Center the window so that its middle is at the current cursor location',
        'command' : self.center, 'shortcut': 'h'},
      { 'kind': 'cascade', 'label': 'Contour levels',
        'tipText': 'Adjust the levels of the contours in the window',
        'submenu' : contour_items, 'shortcut': 'C' },
    ]

    slice_items = []
    for t in (50,100,150,200,250):
      item = {'kind': 'command', 'label': '%d' % t,
              'tipText': '',
              'command': lambda event, thick=t: self.changeSliceThickness(thick) }
      slice_items.append(item)
      
 
    scrollbar_items = [{'kind': 'command', 'label': 'All',
                        'tipText': 'Show navigation scrollbars on all axes of the current window',
                        'command': lambda event: self.setScrollbarVisibility(['x','y','z']) }]
    if self.windowPane.findFirstAxisPanel(label='z1'):  # >= 3D window
      scrollbar_items.append(
                       {'kind': 'command', 'label': 'Z only',
                        'tipText': 'Show navigation scrollbars for only the screen-orthogonal (Z) axes; removes X and Y',
                        'command': lambda event: self.setScrollbarVisibility(['z'])})
    scrollbar_items.append(
                       {'kind': 'command', 'label': 'None',
                        'tipText': 'Remove all navigation scrollbars on all axes of the current window',
                        'command': lambda event: self.setScrollbarVisibility([])})

    crosshair_items = [{'kind': 'command', 'label': 'X & Y',
                        'tipText': 'For both X and Y axes, superimpose 1D slice/transect traces, from the cursor location, on the crosshairs',
                        'command': lambda event: self.setCrosshairTraces(['x','y']) },
                       {'kind': 'command', 'label': 'X',
                        'tipText': 'For only the X axis, superimpose 1D slice/transect traces, from the cursor location, on the crosshairs',
                        'command': lambda event: self.setCrosshairTraces(['x'])},
                       {'kind': 'command', 'label': 'Y',
                        'tipText': 'For only the Y axis, superimpose 1D slice/transect traces, from the cursor location, on the crosshairs',
                        'command': lambda event: self.setCrosshairTraces(['y'])},
                       {'kind': 'command', 'label': 'None',
                        'tipText': 'Remove all 1D slice/transect traces superimposed on crosshairs',
                        'command': lambda event: self.setCrosshairTraces([])},
                      ]    
    
    cloneTip = 'Make a near identical copy of the current window'
    deleteTip = 'Delete the current window'
    printTip = 'Export a PostScript, EPS of PDF rendering of the window for printing etc.'
    propertiesTip = 'Open a popup window to give fine control over spectrum window setting'
    scrollTip = 'Control which window axes, if any, have navigation scrollbars'

    if self.hasValueAxis:
      window_items = [
                      {'kind': 'command', 'label': 'Clone',
                       'tipText': cloneTip,
                       'command' : self.cloneWindow, 'shortcut': 'C' },
                      {'kind': 'command', 'label': 'Delete',
                       'tipText': deleteTip,
                       'command' : self.deleteWindow, 'shortcut': 'D' },
                      {'kind': 'command', 'label': 'Print',
                       'tipText': printTip,
                       'command' : self.printWindow, 'shortcut': 'P' },
                      {'kind': 'command', 'label': 'Window properties',
                       'tipText': propertiesTip,
                       'command' : self.editWindow, 'shortcut': 'W' },
                      {'kind': 'cascade', 'label': 'Scrollbars',
                       'tipText': scrollTip,
                       'submenu' : scrollbar_items, 'shortcut': 'S' },
                     ]
    
    else:
      window_items = [
                      {'kind': 'command', 'label': 'Clone',
                       'tipText': cloneTip,
                       'command' : self.cloneWindow, 'shortcut': 'C' },
                      {'kind': 'command', 'label': 'Delete',
                       'tipText': deleteTip,
                       'command' : self.deleteWindow, 'shortcut': 'D' },
                      {'kind': 'command', 'label': 'Print',
                       'tipText': printTip,
                       'command' : self.printWindow, 'shortcut': 'P' },
                      {'kind': 'command', 'label': 'Window properties',
                       'tipText': propertiesTip,
                       'command' : self.editWindow, 'shortcut': 'W' },
                      {'kind': 'cascade', 'label': 'Scrollbars',
                       'tipText': scrollTip,
                       'submenu' : scrollbar_items, 'shortcut': 'S' },
                      {'kind': 'cascade', 'label': 'Crosshair traces',
                       'tipText': 'Control whether to superimpose 1D slice/transect traces on the crosshairs',
                       'submenu' : crosshair_items, 'shortcut': 'h' },
                      {'kind': 'command', 'label': 'Add side traces',
                       'tipText': 'Add displays of one-dimensional spectrum slices/transects in panels at the edges of X and Y axes',
                       'command' : self.addAxisSlices, 'shortcut': 'A' },
                      {'kind': 'command', 'label': 'Remove side traces',
                       'tipText': 'Remove the one-dimensional spectrum slices/transects panels from the window edges ',
                       'command' : self.removeAxisSlices, 'shortcut': 'R' },
                      {'kind': 'cascade', 'label': 'Side trace size',
                       'tipText': 'Sets how many pixels the one-dimensional spectrum slice/transect panels use for their intensity axis',
                       'submenu' : slice_items, 'shortcut': 'z' },
                     ]

    goto_items = [] # created properly when menu pops up
    self.menu_items = [ \
      { 'kind': 'cascade', 'label': 'Assign',
        'tipText': 'Options relating to the assignment of spectrum peaks to resonances',
        'submenu': assn_items, 'shortcut': 'A' },
      { 'kind': 'cascade', 'label': 'Peak',
        'tipText': 'Options administering peaks and their properties',
        'submenu': peak_items, 'shortcut': 'P' },
      { 'kind': 'cascade', 'label': 'Locate Peaks',
        'tipText': 'Options for finding or navigating to peak positions,',
        'submenu': peak_nav_items, 'shortcut': 'L' },
      { 'kind': 'cascade', 'label': 'View',
        'tipText': 'Options to change the viewed region of contours',
        'submenu': view_items, 'shortcut': 'V' },
      { 'kind': 'cascade', 'label': 'Navigate',
        'tipText': 'Options to navigate to equivalent positions in different windows',
        'submenu': goto_items, 'shortcut': 'N' },
      { 'kind': 'cascade', 'label': 'Strip',
        'tipText': 'Options to administer strips and orthogonal sub-divisions of windows',
        'submenu': strip_items, 'shortcut': 'S' },
      { 'kind': 'cascade', 'label': 'Markers',
        'tipText': 'Options to add and remove cross marks and ruler lines',
        'submenu': mark_items, 'shortcut': 'M' },
      { 'kind': 'cascade', 'label': 'Window',
        'tipText': 'Options for general window properties and one-dimensional slice displays',
        'submenu': window_items, 'shortcut': 'W' },
      { 'kind': 'separator'},
      { 'kind': 'cascade', 'label': 'Macros',
        'tipText': 'Python macro scripts that may be executed from the window menu',
        'submenu': macro_items, 'shortcut': 'c' },
      #{ 'kind': 'cascade', 'label': 'Split screen', 'submenu': split_items },
    ]


  def changeSliceThickness(self, value):

    if self.windowPane:
      for axisPanel in self.windowPane.axisPanels:
        if axisPanel.label in ('x','y'):
          axisPanel.slicePanel.thickness = value

  def addAxisSlices(self, *event):

    if self.windowPane:
      for axisPanel in self.windowPane.axisPanels:
        if axisPanel.label in ('x','y'):
          axisPanel.slicePanel.isVisible = True

  def removeAxisSlices(self, *event):

    if self.windowPane:
      for axisPanel in self.windowPane.axisPanels:
        if axisPanel.label in ('x','y'):
          axisPanel.slicePanel.isVisible = False

  def cloneWindow(self, *event):

    if self.windowPane:
      self.topPopup.cloneWindow(self.windowPane.spectrumWindow)

  def deleteWindow(self, *event):

    if self.windowPane:
      if showYesNo('Delete window', 'Do you really want to delete this window?', parent=self):
        self.windowPane.spectrumWindow.delete()

  def printWindow(self, *event):

    self.topPopup.printWindow(self.windowPane.spectrumWindow)

  # TBD: dangerous way this done here, using explicit indices into submenus
  def updateSliceMenuState(self):

    slice_menu = self.scrolled_window.sliceMenu

    views = self.getActiveSliceViews()
    if len(views) == 1:
      state = Tkinter.NORMAL
    else:
      state = Tkinter.DISABLED

    slice_menu.entryconfig(0, state=state) # reference

  # TBD: more general
  def setSliceMenuItems(self):

    self.slice_menu_items = [ \
      { 'kind': 'command', 'label': 'Reference...',
        'tipText': 'Use the cursor position in the side panel to set spectrum point to ppm referencing',
        'command' : self.reference, 'shortcut': 'R' },
    ]

  def assignPeak(self, *event):

    self.topPopup.assignmentPanel()
    if self.menuPeak:
      self.topPopup.popups['edit_assignment'].update(self.menuPeak)

  def editPeakDetails(self, event):

    if self.menuPeak:
      analysisProject = self.menuPeak.root.currentAnalysisProject

      doDetails = analysisProject.doDetailAnnotations
      doMerit = analysisProject.doMeritAnnotations
      
      analysisProject.doDetailAnnotations = True
      analysisProject.doMeritAnnotations = False

      details = askString('Peak Details','Enter comment:',self.menuPeak.details, parent=self) or None
      
      # Need a clear out first; to force a redraw
      # Otherwise merit can override if details unchanged
      self.menuPeak.setDetails(None)
      self.menuPeak.setDetails(details)

      analysisProject.doDetailAnnotations = doDetails
      analysisProject.doMeritAnnotations = doMerit
      
  def unaliasPeak(self, event):

    if self.menuPeak:
      self.topPopup.editPeakAliasing(peak=self.menuPeak)

  def setPeakMerit(self, merit):

    ##if self.menuPeak:
    ##  self.menuPeak.figOfMerit = min(1,max(0,merit))
    ##  makePeakAnnotation(self.menuPeak)

    peaks = self.topPopup.currentPeaks
    if not peaks:
      return

    merit = min(1,max(0,merit))

    if len(peaks) > 1:
      if not showYesNo('Set merit', 'Do you really want to set merit to %.1f for %d peaks?' % (merit, len(peaks)), parent=self):
        return

    for peak in peaks:
      peak.figOfMerit = merit
      makePeakAnnotation(peak)

  def fitPeaks(self, fitMethod, fitType=None):

    peaks = self.topPopup.selected_objects
    if not peaks:
      return

    if not fitType:
      fitType = PEAK_FIT_TYPES[0]

    if fitType == PEAK_FIT_TYPES[0]:
      for peak in peaks:
        fitPeaks([peak], fitMethod)
    else:
      fitPeaks(peaks, fitMethod)

  def movePeaksToPeakList(self, peakList):
    
    peaks = self.topPopup.selected_objects
    if not peaks:
      return
      
    copyPeaksToPeakList(peaks, peakList)
    self.topPopup.deleteSelected()
    
  def deleteSelectedPeaks(self, *event):

    self.topPopup.queryDeleteSelected(self)

  def showSelectedPeaks(self, *event):

    self.topPopup.viewSelectedPeaks()

  def findClosePeaks(self, peakDim, spectrum=None):

    nmrProject = peakDim.topObject
    expDimRefs = peakDim.dataDim.expDim.expDimRefs
    peaks      = []

    if spectrum:
      peakList = spectrum.activePeakList
      peaks.extend( findMatchingPeaks(peakList, peakDim) )

    else:
      for experiment in nmrProject.experiments:
        for spectrum1 in experiment.dataSources:
          peakList = spectrum1.activePeakList
          peaks.extend( findMatchingPeaks(peakList, peakDim) )

    # tbd with isotope checks
    peaks2 = []
    ppm = peakDim.value
    for peak in peaks:
      peakDims = peak.sortedPeakDims()
      minD = abs(peakDims[0].value-ppm)
      for peakDim2 in peakDims[1:]:
        d = abs(peakDim2.value-ppm)
        if d < minD:
          minD = d
      peaks2.append( (minD, peak) )

    peaks2.sort()
    peaks = [x[1] for x in peaks2]

    self.topPopup.viewPeaks(peaks)

  def stripSymmetryRelatedPeaks(self, *event):

    if self.menuPeak:
      peakLists = self.getActivePeakLists() or None
      peaks = findSymmetryPeaks(self.menuPeak, peakLists=peakLists)
      if peaks:
        if len(peaks) > 9:
          peaks = peaks[:9]

        createPeakMark(self.menuPeak, lineWidth=2)
        peaks = [self.menuPeak,] + peaks
        displayPeakStrips(self.topPopup, peaks, self.windowPane)
        self.update_idletasks()
        displayPeakStrips(self.topPopup, peaks, self.windowPane)
      else:
        showWarning('Warning','None found', parent=self)


  def findSymmetryRelatedPeaks(self, *event):

    if self.menuPeak:
      peakLists = self.getActivePeakLists() or None
      peaks = findSymmetryPeaks(self.menuPeak, peakLists=peakLists)
      if peaks:
        self.topPopup.viewPeaks(peaks)

      else:
        showWarning('Warning','None found', parent=self)



  def propagatePeakAssignments(self, event):

    peaks = self.topPopup.currentPeaks
    if peaks:
      propagatePeakAssignments(peaks, warnUnalias=True)

  def propagatePeakAssignmentsToSister(self, peakList):

    peaks = list(self.topPopup.currentPeaks)
    if peaks:
      peaks.extend(peakList.sortedPeaks())
      propagatePeakAssignments(peaks, warnUnalias=True)

  def propagatePeakUnaliasing(self, event):

    refPeak = self.menuPeak
    peaks   = self.topPopup.currentPeaks
    if refPeak and peaks:
      message = 'Propagate unaliasing from peak\nto all selected peaks?\n'

      i = 0
      for peakDim in refPeak.sortedPeakDims():
        i += 1
        message += 'Dimension %d: Num aliasing %d\n(Ref peak position %f)\n' % (i,peakDim.numAliasing, peakDim.value)

      if showOkCancel('Confirm',message, parent=self):
        propagatePeakUnaliasing(refPeak, peaks)

  def uniteResonancePositions(self, event):
  
    peaks = self.topPopup.currentPeaks
    if len(peaks) > 1:
      if not showOkCancel('Confirm','Do you really want to unite resonance positions for the %d selected peaks? This will align other peaks which carry the same assignments.' % len(peaks), parent=self):
        return
    if not peaks and self.menuPeak:
      peaks = [self.menuPeak]
    for peak in peaks: 
      uniteResonancePeakPositions(peak)

  def addPeakSpinSystem(self, event):

    peaks = self.topPopup.currentPeaks
    if peaks:
      addPeakResonancesToSpinSystem(peaks)

  def clearPeakAssignments(self, event):

    peaks = self.topPopup.currentPeaks

    if peaks:
      if showOkCancel('Confirm','Clear assignments for %d selected peaks?' % len(peaks), parent=self):
        removePeaksAssignment(peaks)

  def clearSeqLinks(self, event):

    peaks = self.topPopup.currentPeaks

    if peaks:
      if showOkCancel('Confirm','Clear all sequential spin system connections for %d selected peaks?' % len(peaks), parent=self):
        spinSystems = []
        for peak in peaks:
          for peakDim in peak.peakDims:
            for contrib in peakDim.peakDimContribs:
              spinSystem = contrib.resonance.resonanceGroup
              if spinSystem and (spinSystem not in spinSystems):
                spinSystems.append(spinSystem)

        for spinSystem in spinSystems:
          clearSeqSpinSystemLinks(spinSystem, delta=None)

  def createMark(self, event):

    self.doCreateMark(event.widget, event.x, event.y)

  def createHorizontalRuler(self, event):

    self.doCreateHorizontalRuler(event.widget, event.x, event.y)

  def createVerticalRuler(self, event):

    self.doCreateVerticalRuler(event.widget, event.x, event.y)

  def createPeak(self, event):

    canvas = event.widget
    (a, b, x, y) = self.calcWorldCoord(canvas, event.x, event.y)
    (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
    self.createNewPeak(a, b, row, col)

  def clearMarks(self, *event):

    removeMarks(self.windowPane.root, removeAll=True)

  def clearRulers(self, *event):

    removeRulers(self.windowPane.root, removeAll=True)

  def changeScale(self, factor):

    views = self.getActiveSpectrumViews()
    for view in views:
      analysisSpectrum = view.analysisSpectrum
      changeMode = analysisSpectrum.autoLevelMode
      
      if changeMode == 'add':
        if factor <= 0.0:
          continue
        factor = math.log(factor)
        
      Util.changeSpectrumContourLevels(analysisSpectrum, factor, changeMode)

  def center(self, event):

    self.centerAtLocation(event.widget, event.x, event.y)

  def centerAtLocation(self, canvas, x, y):

    (a, b, s, t) = self.calcWorldCoord(canvas, x, y)
    #print 'center', a, b, s, t

    s = s - 0.5
    if (not self.hasValueAxis):
      t = t - 0.5
    else:
      t = 0
    #self.scrolled_window.translate(canvas, s, t)
    self.translate(canvas, s, t)

  def translate(self, canvas, sx, sy):

    #f = 0.5
    #self.scrolled_window.translate(canvas, f*sx, f*sy)
    self.scrolled_window.translate(canvas, sx, sy)

  def windowsZoom(self, event):

    #print "windowsZoom", event

    delta  = event.delta
    canvas = event.widget
    state  = event.state

    if state & 4: # Control
      self.scrollZPlane(canvas, 'z1', delta)

    elif state & 1: # Shift
      self.scrollZPlane(canvas, 'z2', delta)

    ## djo35 - in Windows 7 state seems to be constantly 8
    #elif state & 8: # Alt
    #  self.scrollZPlane(canvas, 'z1', delta)

    else:
      if delta > 0:
        self.scrolled_window.zoom(canvas, 0.8)
      elif delta < 0:
        self.scrolled_window.zoom(canvas, 1.2)

  def zoomIn(self, event):

    state  = event.state
    canvas = event.widget
    if state & 4: # Control
      self.scrollZPlane(canvas, 'z1', -1)

    elif state & 1: # Shift
      self.scrollZPlane(canvas, 'z2', -1)

    elif state & 8: # Alt
      self.scrollZPlane(canvas, 'z1', -1)

    else:
      self.scrolled_window.zoom(canvas, 0.8)

  def zoomOut(self, event):

    state  = event.state
    canvas = event.widget
    if state & 4: # Control
      self.scrollZPlane(canvas, 'z1', 1)

    elif state & 1: # Shift
      self.scrollZPlane(canvas, 'z2', 1)

    elif state & 8: # Alt
      self.scrollZPlane(canvas, 'z1', 1)

    else:
      self.scrolled_window.zoom(canvas, 1.2)

  def scrollZPlane(self, canvas, label, step):

    try:
      # it's possible some other widget has focus, so an exception can be thrown
      (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
    except:
      return

    windowPane = self.windowPane
    axisPanelZ = windowPane.findFirstAxisPanel(label=label)

    if axisPanelZ:
      stripAxis = windowPane.spectrumWindow.stripAxis
      axisPanel = windowPane.findFirstAxisPanel(label=stripAxis)

      if stripAxis == 'x':
        strip = col
      else:
        strip = row

      axisRegionsZ = axisPanelZ.sortedAxisRegions()

      if strip < len(axisRegionsZ):
        self.orthogScroll(axisRegionsZ[strip], step)


  def zoom(self, canvas, scale):

    self.scrolled_window.zoom(canvas, scale)

  def zoomToSpectra(self, canvas):
    """ zoom x,y region just enough to see all visible spectra
    """
    (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
    self.zoomToSpectraRegion(row, col)

  def orthogScroll(self, axisRegion, step):

    axisPanel = axisRegion.axisPanel
    (w0, w1) = self.getOrthogonalWorldRegion(axisPanel)
    if axisPanel.axisUnit and axisPanel.axisUnit.isBackwards:
      step = -step
      (w0, w1) = (w1, w0)
    (r0, r1) = axisRegion.region
    d = r1 - r0
    #print 'orthogScroll1', step, d, (r0, r1), (w0, w1)
    r0 = r0 + step*d
    r1 = r0 + d
    #print 'orthogScroll2', (r0, r1)

    if (r0 < w0):
      r0 = w0
      r1 = w0 + d
    elif (r1 > w1):
      r0 = w1 - d
      r1 = w1

    axisRegion.region = (r0, r1)

  def orthogScrollRight(self, canvas):

    axisPanels = self.windowPane.sortedAxisPanels()
    if len(axisPanels) > 2:
      (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
      axisPanel = axisPanels[0]
      axisRegion = axisPanel.sortedAxisRegions()[col]
      axisPanel = axisPanels[2]
      if axisRegion.isActive and hasattr(axisPanel, 'region_selector'):
        axisPanel.region_selector.moveScrollbar(bottom_right_trough_mode)
      else:
        axisRegion = axisPanel.sortedAxisRegions()[col]
        self.orthogScroll(axisRegion, 1)

  def orthogScrollLeft(self, canvas):

    axisPanels = self.windowPane.sortedAxisPanels()
    if len(axisPanels) > 2:
      (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
      axisPanel = axisPanels[0]
      axisRegion = axisPanel.sortedAxisRegions()[col]
      axisPanel = axisPanels[2]
      if axisRegion.isActive and hasattr(axisPanel, 'region_selector'):
        axisPanel.region_selector.moveScrollbar(top_left_trough_mode)
      else:
        axisRegion = axisPanel.sortedAxisRegions()[col]
        self.orthogScroll(axisRegion, -1)

  def calculateNoise(self, canvas):

    # TBD: 1D windows??
    if self.hasValueAxis:
      return

    draggedBox = self.draggedBox
    if draggedBox:
      (a0, b0, a1, b1) = draggedBox
      xmin = min(a0, a1)
      xmax = max(a0, a1)
      xregion = (xmin, xmax)
      ymin = min(b0, b1)
      ymax = max(b0, b1)
      yregion = (ymin, ymax)
      (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
      windowPane = self.windowPane
      for view in windowPane.spectrumWindowViews:
        if view.isPosVisible or view.isNegVisible:
          analysisSpectrum = view.analysisSpectrum
          spectrum = analysisSpectrum.dataSource
          if not hasattr(spectrum, 'block_file') or not spectrum.block_file:
            continue
          ndim = spectrum.numDim
          boxMin = ndim*[0]
          boxMax = ndim*[0]
          for axisMapping in view.axisMappings:
            label = axisMapping.label
            axisPanel = self.windowPane.findFirstAxisPanel(label=label)
            if label == 'x':
              region = xregion
            elif label == 'y':
              region = yregion
            else:
              if self.windowPane.spectrumWindow.stripAxis == 'x':
                n = col
              else:
                n = row
              axisRegion = axisPanel.sortedAxisRegions()[n]
              region = axisRegion.region
            analysisDataDim = axisMapping.analysisDataDim
            dataDim = analysisDataDim.dataDim
            (region0, region1) = Util.convertRegion(region, axisPanel.axisUnit, dataDim)
            region0 = int(math.floor(region0))
            region1 = int(math.ceil(region1))
            if region0 < 0 or region1 > dataDim.numPoints:
              print 'Outside fundamental region for spectrum %s:%s and noise calculation only works inside that' %(spectrum.experiment.name, spectrum.name)
              break
            boxMin[dataDim.dim-1] = region0
            boxMax[dataDim.dim-1] = region1
          else:
            calculateNoiseInBox(spectrum, boxMin, boxMax)

  def reference(self, event):

    views = self.getActiveSpectrumViews()
    if (len(views) != 1):
      return

    view = views[0]

    slice = event.widget
    if (slice.orient == Tkinter.HORIZONTAL):
      w = slice.slice_width
      x = event.x
      s = (x + 0.5) / w
      axisPanel = self.windowPane.findFirstAxisPanel(label='x')
    else:
      h = slice.slice_height
      y = h - 1 - event.y
      s = (y + 0.5) / h
      axisPanel = self.windowPane.findFirstAxisPanel(label='y')

    axisMapping = view.findFirstAxisMapping(label=axisPanel.label)
    if (not axisMapping):
      return

    dataDim = axisMapping.analysisDataDim.dataDim
    (t0, t1) = slice.view_region
    t = t0 * (1 - s) + t1 * s
    t = Util.convertPosition(t, getPrimaryDataDimRef(dataDim), fromUnit=axisPanel.axisUnit.unit)
    popup = SetReferencePopup(self, dataDim, t)
    popup.destroy()

  def editWindow(self, *event):

    self.topPopup.editWindow(self.windowPane)

  def isInAllowedRegion(self, dataDim, label, region):

    if label in ('x', 'y'):
      r0 = r1 = region
    else:
      (r0, r1) = region

    if dataDim.className == 'FreqDataDim':
      dataDimRef = getPrimaryDataDimRef(dataDim)
      expDimRef = dataDimRef.expDimRef
      minAliasedFreq = expDimRef.minAliasedFreq
      maxAliasedFreq = expDimRef.maxAliasedFreq
      if (minAliasedFreq is None):
        minFreq = Util.convertPosition(dataDim.numPoints, dataDimRef, toUnit='ppm')
      else:
        minFreq = minAliasedFreq
      if r1 < minFreq:
        return False
      if (maxAliasedFreq is None):
        maxFreq = Util.convertPosition(1.0, dataDimRef, toUnit='ppm')
      else:
        maxFreq = maxAliasedFreq
      if r0 > maxFreq:
        return False
    else:
      # TBD: is this correct?
      if r1 < 0:
        return False
      if r0 > dataDim.numPoints:
        return False

    return True

  def findNearestPeak(self, position_region, xAxisRegion, yAxisRegion):

    peakMin = None
    activeViews = self.getActiveSpectrumViews()
    d2Min = None
    for view in activeViews:
      result = self.findNearbyViewPeak(view, position_region, xAxisRegion, yAxisRegion)
      if result:
        (peak, d2) = result
        if (d2Min is None) or (d2 < d2Min):
          peakMin = peak
          d2Min = d2

    return peakMin

  def findNearbyViewPeak(self, view, position_region, xAxisRegion, yAxisRegion):

    spectrum = view.analysisSpectrum.dataSource
    peakList = spectrum.activePeakList
    
    if not peakList:
      return None
      
    analysisPeakList = peakList.analysisPeakList
    winPeakList = view.findFirstWindowPeakList(analysisPeakList=analysisPeakList)
    
    if not winPeakList:
      return None
      
    if not winPeakList.isSymbolDrawn:
      return None
      
    (xscale, yscale) = self.getPeakScale(view, xAxisRegion, yAxisRegion)

    spectrum_region = spectrum.numDim * [0]
    params = getPeakFindParams(self.windowPane.root)
    thickness = params['thickness'] + self.extraPad
    if self.hasValueAxis:
      ydim = -1
    for axisMapping in view.axisMappings:
      dataDim = axisMapping.analysisDataDim.dataDim
      dim = dataDim.dim - 1
      label = axisMapping.label
      if (not self.isInAllowedRegion(dataDim, label, position_region[label])):
        return None
      axisPanel = self.windowPane.findFirstAxisPanel(label=label)
      if label == 'x':
        xdim = dim
        tt = self.nearbyThickness(dataDim, axisPanel)
      elif label == 'y':
        ydim = dim
        tt = self.nearbyThickness(dataDim, axisPanel)
      else:
        tt = thickness
      r = position_region[label]
      (r0, r1) = self.convertPositionRegion(r, axisPanel, dataDim)
      spectrum_region[dim] = [ r0-tt, r1+tt ]

    #print 'findNearbyViewPeak', spectrum.name, spectrum_region
    result = self.findNearbyPeak(peakList, spectrum_region, xdim, ydim, xscale, yscale)

    return result

  def nearbyThickness(self, dataDim, axisPanel):

    axisRegion = axisPanel.findFirstAxisRegion()
    if self.hasValueAxis and axisPanel.label == 'y':
      (r0, r1) = axisRegion.region
    else:
      (r0, r1) = Util.convertRegion(axisRegion.region, axisPanel.axisUnit, dataDim)

    pixels = 7
    thickness = abs((r1-r0) * pixels / axisRegion.size)

    return thickness

  def selectSingle(self, canvas, a, b, x, y, button = 1,
                   state = no_key_state, newSelection = True, event = None):

    self.unpostMenu()

    if (button != 1):
      return

    topPopup = self.topPopup
    topPopup.startSelection()

    if newSelection:
      topPopup.clearSelected()

    (row, col) = self.scrolled_window.getCanvasRowCol(canvas)

    if self.windowPane.spectrumWindow.stripAxis == 'x':
      stripNum = col
    else:
      stripNum = row

    if state != ctrl_key_state:
      getPanel = self.windowPane.findFirstAxisPanel

      xAxisRegion = getPanel(label='x').sortedAxisRegions()[col]
      yAxisRegion = getPanel(label='y').sortedAxisRegions()[row]

      position_region = self.findPositionRegion(a, b, stripNum)

      #print 'selectSingle1', position_region
      activeViews = self.getActiveSpectrumViews()
      for view in activeViews:
        result = self.findNearbyViewPeak(view, position_region, xAxisRegion, yAxisRegion)
        if result:
          (peak, d2Min) = result
          if (not newSelection and topPopup.isSelected(peak)):
            removePeakFromSelected(peak, topPopup)
          else:
            self.menuPeak = peak
            addPeakToSelected(peak, parent=topPopup)
        else:
          self.menuPeak = None

    else:
      self.menuPeak = self.createNewPeak(a, b, row, col)

    self.setCurrentObjects(event)
    self.topPopup.endSelection()
    self.drawCanvas(canvas)

  def createNewPeak(self, a, b, row, col):

    if self.windowPane.spectrumWindow.stripAxis == 'x':
      n = col
    else:
      n = row

    position = self.findPosition(a, b, n)
    activeViews = self.getActiveSpectrumViews()
    windowPane = self.windowPane
    
    for view in activeViews:
      spectrum = view.analysisSpectrum.dataSource
      peakList = spectrum.activePeakList
      if not peakList:
        continue
        
      analysisPeakList = peakList.analysisPeakList
      winPeakList = view.findFirstWindowPeakList(analysisPeakList=analysisPeakList)
      if not winPeakList:
        continue
      if not winPeakList.isSymbolDrawn:
        continue
        
      # Check we're not picking outside contoured region
      axisMapping = getDataDimAxisMapping(spectrum, windowPane)
      for axis in position.keys():
        if axisMapping.has_key(axis): # might not if lower dimensional spectrum
          dataDim = axisMapping[axis]
          if dataDim.className == 'FreqDataDim':
            ppm = position[axis]
        
            for dataDimRef in dataDim.dataDimRefs:
              ppmMin, ppmMax = getDataDimRefFullRange(dataDimRef)

              if ppmMin < ppm < ppmMax:
                break
          
            else:
              break
          else:
            pnt = position[axis]
            if pnt > dataDim.numPoints:
              break
      else:
        specPosition, specTile = self.determinePosition(view, position)
        dataDim = spectrum.findFirstDataDim(className='SampledDataDim')
        if dataDim: # Pseudo 3Ds etc.
          axisMapping = view.findFirstAxisMapping(analysisDataDim=dataDim.analysisDataDim)
          axisPanel   = self.windowPane.findFirstAxisPanel(label=axisMapping.label)
          axisRegion  = axisPanel.findFirstAxisRegion()
          region = axisRegion.region
 
          points = range(int(region[0]),int(region[1])+1)
          index = dataDim.dim -1
 
          for point in points:
            specPosition[index] = point
            addPeak(peakList, position=specPosition, tile=specTile,
                    parent=self.topPopup, doFit=False)
 
        else:
          #print 'createNewPeak', position, spectrum_position, spectrum_tile
          addPeak(peakList, position=specPosition, tile=specTile,
                 parent=self.topPopup, doFit=False)
 
 

  def examineRegion(self, canvas, a0, b0, a1, b1, x0, y0, x1, y1, state=no_key_state):

    # TBD: do more
    #print 'examineRegion1', a0, b0, a1, b1, x0, y0, x1, y1

    (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
    if self.windowPane.spectrumWindow.stripAxis == 'x':
      stripNum = col
    else:
      stripNum = row

    region = self.findRegion(min(a0,a1), min(b0,b1), max(a0,a1), max(b0,b1), n=stripNum)
    #print 'examineRegion2', region

    if (state == shift_key_state or state == no_key_state):

      if state == no_key_state:
        self.topPopup.clearSelected()

      params = getPeakFindParams(self.windowPane.root)
      thickness = params['thickness'] + self.extraPad

      # TBD: check more, i.e. not just peaks
      views = self.getActiveSpectrumViews()
      for view in views:
        spectrum = view.analysisSpectrum.dataSource
        spectrum_region = spectrum.numDim * [0]
        spectrum_thickness = spectrum.numDim * [0]
        for axisMapping in view.axisMappings:
          label = axisMapping.label
          if self.hasValueAxis and label == 'y':
            continue
          if (label in ('x', 'y')):
            t = 0
          else:
            t = thickness
          dim = axisMapping.analysisDataDim.dataDim.dim - 1
          spectrum_region[dim] = region[label]
          spectrum_thickness[dim] = t
          
        peakList = spectrum.activePeakList
        if peakList:
          analysisPeakList = peakList.analysisPeakList
          winPeakList = view.findFirstWindowPeakList(analysisPeakList=analysisPeakList)
          if winPeakList and winPeakList.isSymbolDrawn:
            peakLists = [ peakList ]
            searchPeaks(peakLists, region=spectrum_region, parent=self.topPopup,
                        thickness=spectrum_thickness)

    elif (state == ctrl_key_state):

      views = self.getActiveSpectrumViews()
      for view in views:
        spectrum = view.analysisSpectrum.dataSource
        spectrum_region = spectrum.numDim * [0]
        for axisMapping in view.axisMappings:
          label = axisMapping.label
          if self.hasValueAxis and label == 'y':
            continue
          dim = axisMapping.analysisDataDim.dataDim.dim - 1
          spectrum_region[dim] = region[label]
        result = getRegionStats(spectrum, spectrum_region)
        if result:
          n, avg, std = result
          print('spectrum=%s:%s, npts=%d, sum=%.3e, avg=%.3e, std=%.3e' % (spectrum.experiment.name, spectrum.name, n, n*avg, avg, std))
        else:
          print('stats could not be obtained (spectrum file not available)')

    elif (state == (shift_key_state + ctrl_key_state)):

      params = getPeakFindParams(self.windowPane.root)
      thickness = params['thickness']

      views = self.getActiveSpectrumViews()
      for view in views:
        spectrum = view.analysisSpectrum.dataSource
        spectrum_region = spectrum.numDim * [0]
        spectrum_thickness = spectrum.numDim * [0]
        peakList = spectrum.activePeakList
        if (not peakList):
          continue
        winPeakList = view.findFirstWindowPeakList(analysisPeakList=peakList.analysisPeakList)
        if (not winPeakList):
          continue
        if not winPeakList.isSymbolDrawn:
          continue
        for axisMapping in view.axisMappings:
          label = axisMapping.label
          if self.hasValueAxis and label == 'y':
            continue
          dataDim = axisMapping.analysisDataDim.dataDim
          if dataDim.className == 'SampledDataDim':
            (r0, r1) = region[label]
            r1 = min(r1, dataDim.numPoints)
            if r1 < r0:
              break
            region[label] = (r0, r1)
          if (label in ('x', 'y')):
            t = 0
          else:
            t = thickness
          dim = dataDim.dim - 1
          spectrum_region[dim] = region[label]
          spectrum_thickness[dim] = t
        else:
          findPeaks(peakList, region=spectrum_region, parent=self.topPopup,
                    thickness=spectrum_thickness)

    # TBD: why was below here?
    #canvas.handler.mapRanges(x0, y0, x1, y1, a0, b0, a1, b1)

    # TBD: below might not be needed now because parent redraws all
    # but clearXor seems necessary so leave both for now
    #self.drawCanvas(canvas)
    if canvas.handler:
      canvas.handler.clearXor()

  def selectMulti(self, canvas, a0, b0, a1, b1, x0, y0, x1, y1,
                  button = 1, state = no_key_state, event = None):

    self.draggedBox = None

    self.unpostMenu()
    
    #print 'selectMulti', state
    #print 'selectMulti', a0, b0, a1, b1, x0, y0, x1, y1
    self.topPopup.currentRegion = (a0, b0, a1, b1)
    self.topPopup.startSelection()

    if ((x0 == x1) and (y0 == y1)):
      if state == shift_key_state:
        newSelection = False
      else:
        newSelection = True
      self.selectSingle(canvas, a0, b0, x0, y0, state=state, newSelection=newSelection, event=event)
    else:
      self.examineRegion(canvas, a0, b0, a1, b1, x0, y0, x1, y1, state=state)

    if hasattr(canvas, 'prev_cross'):
      del canvas.prev_cross
      
    if hasattr(canvas, 'prev_box'):
      del canvas.prev_box

    self.setCurrentObjects(event)
    self.topPopup.endSelection(redraw=True)
    self.drawCanvas(canvas)

  def dragBox(self, canvas, a0, b0, a1, b1, x0, y0, x1, y1,
              button = 1, state = no_key_state, event = None):

    #print 'dragBox', a0, b0, a1, b1, x0, y0, x1, y1

    self.draggedBox = (a0, b0, a1, b1)

    handler = canvas.handler
    if not handler:
      return

    handler.makeCurrent()
    handler.startXor()

    # TBD: more general color
    color = (0.3, 0.5, 0.7)
    handler.setColor(color)

    if hasattr(canvas, 'doubleBuffer') and not canvas.doubleBuffer:
      doubleBuffer = False
    else:
      doubleBuffer = True

    if not doubleBuffer and hasattr(canvas, 'prev_box'):
      (prev_a0, prev_b0, prev_a1, prev_b1, prev_x0, prev_y0, prev_x1, prev_y1) = canvas.prev_box
      if self.hasValueAxis:
        handler.mapRanges(prev_x0, 0.0, prev_x1, 1.0, prev_a0, 0.0, prev_a1, 1.0)
        handler.drawXorBox(prev_a0, 0.0, prev_a1, 1.0)
      else:
        handler.mapRanges(prev_x0, prev_y0, prev_x1, prev_y1, prev_a0, prev_b0, prev_a1, prev_b1)
        handler.drawXorBox(prev_a0, prev_b0, prev_a1, prev_b1)

    if self.hasValueAxis:
      handler.mapRanges(x0, 0.0, x1, 1.0, a0, 0.0, a1, 1.0)
      handler.drawXorBox(a0, 0.0, a1, 1.0)
    else:
      handler.mapRanges(x0, y0, x1, y1, a0, b0, a1, b1)
      handler.drawXorBox(a0, b0, a1, b1)

    handler.finishXor()

    if not doubleBuffer:
      handler.flush()
      canvas.prev_box = (a0, b0, a1, b1, x0, y0, x1, y1)

  # TBD: checks in "changed" functions below are not very clever
  # (or very useful) but at some point something better could perhaps be done
  def changedSpectrum(self, spectrum):

    if isSpectrumInWindowPane(self.windowPane, spectrum):
      self.drawAllAfter()

  def changedPeakList(self, peakList):

    if isSpectrumInWindowPane(self.windowPane, peakList.dataSource):
      self.drawAllAfter()

  def changedWinPeakList(self, winPeakList):

    if self.windowPane == winPeakList.spectrumWindowView.spectrumWindowPane:
      self.windowPopup.setPeaksSelectorAfter()
      self.drawAllAfter()

  def changedPeak(self, peak):

    if isSpectrumInWindowPane(self.windowPane, peak.peakList.dataSource):
      self.drawAllAfter()

  def changedPeakDim(self, peakDim):

    if isSpectrumInWindowPane(self.windowPane, peakDim.peak.peakList.dataSource):
      self.drawAllAfter()

  def changedMark(self, *markOrMarkDim):

    #self.after_idle(self.drawAll)
    self.drawAllAfter()

  def changedMaxMarks(self, analysisProject):

    removeMarks(self.windowPane.root)

  def changedRuler(self, *ruler):

    #self.after_idle(self.drawAll)
    self.drawAllAfter()

  def changedMaxRulers(self, analysisProject):

    removeRulers(self.windowPane.root)

  def changedStoredContour(self, storedContour):

    try:
      if isSpectrumInWindowPane(self.windowPane, storedContour.dataSource):
        self.drawAllAfter()
    except:
      pass

  def changedAppData(self, appData):

    if (appData.application != self.windowPane.root.application.name):
      return

    keyword = appData.keyword
    if keyword not in ('textOffset',):
      return

    # TBD: is this needed any more??
    if self.handlerClass and self.handlerClass == TkHandler.TkHandler:
      # for some reason in Tk the font gets changed but
      # the old font is still used until another draw is done
      if keyword in ('fontName', 'fontSize'):
        self.fontChanged = True

    self.drawAllAfter()

  def initAxisMapping(self, axisMapping):

    view = axisMapping.spectrumWindowView
    analysisSpectrum = view.analysisSpectrum
    spectrum = analysisSpectrum.dataSource
    if len(view.axisMappings) == spectrum.numDim:
      self.changedReferencing()

  def changedReferencing(self, *extra):

    self.setMinThickness()
    self.updateOrthogonalWorlds()
    self.windowPopup.updatePositionView()
    self.drawAllAfter()

  def changedColor(self, *color):

    self.setHandlerBackground()
    #self.after_idle(self.drawAll)
    self.drawAllAfter()

  def initSpectrumWindowView(self, spectrumWindowView):

    if not spectrumWindowView.spectrumWindowPane == self.windowPane:
      return

    self.after_idle(self.setMinThickness)

    self.after_idle(lambda spectrumWindowView=spectrumWindowView: \
                    self.changedSpectrumWindowView(spectrumWindowView))
    self.after_idle(self.updateOrthogonalWorlds)

  def setMinThickness(self):

    for axisPanel in self.windowPane.sortedAxisPanels()[2:]:
      if hasattr(axisPanel, 'region_selector') and not axisPanel.axisType.isSampled:
        min_thickness = self.minAxisPanelThickness(axisPanel)
        axisPanel.region_selector.setMinThickness(min_thickness)

  def deleteSpectrumWindowView(self, spectrumWindowView):

    self.windowPopup.setSpectrumSelector()
    
    if not self.windowPane.isDeleted:
      self.after_idle(self.drawAll)
      self.after_idle(self.setMinThickness)

  def changedSpectrumWindowView(self, spectrumWindowView):

    self.windowPopup.setSpectrumSelector()
    if spectrumWindowView.spectrumWindowPane == self.windowPane:
      self.drawAllAfter()

  def changedSpectrumWindowViewSlice(self, spectrumWindowView):

    if spectrumWindowView.spectrumWindowPane == self.windowPane:
      if self.hasValueAxis:
        self.windowPopup.setSpectrumSelector()
        self.drawAll()
      else:
        self.drawAllSlices()

  def changedAxisTypeRegion(self, *axisType):

    world_region  = self.getWorldRegion()
    xview_regions = self.getXviewRegions(world_region.x_region)
    yview_regions = self.getYviewRegions(world_region.y_region)

    self.scrolled_window.setWorldRegion(world_region, xview_regions, yview_regions)
    self.drawAllAfter()

  def changedViewRegion(self, axisRegion):

    # TBD: not sure if below is good enough
    panel = axisRegion.axisPanel
    #print 'changedViewRegion', panel.label, axisRegion.region, panel.spectrumWindow.name, self.windowPane.name
    if panel.spectrumWindowPane == self.windowPane:
      world_region = self.getWorldRegion()
      #print 'changedViewRegion0', self.windowPane.name, axisRegion.region
      view_region = Util.checkSwapRegion(axisRegion.region, panel.axisUnit)
      #print 'changedViewRegion1', self.windowPane.name, view_region, axisRegion.region
      axisRegions = panel.sortedAxisRegions()
      n = axisRegions.index(axisRegion)
      #print 'changedViewRegion2', self.windowPane.name, n, panel.label
      if (panel.label == 'x'):
        Util.fitViewInWorld(view_region, world_region.x_region)
        #self.scrolled_window.setXviewRegion(view_region, n)
        (a0, a1) = view_region
        self.scrolled_window.setXviewRegion(Region1D(a0, a1), n)
        if not self.updatingAspectRatio:
          self.after_idle(lambda: self.updateAspectRatio(doDraw=False))
        self.drawCol(n)
        self.drawSliceCol(n)
      else:
        Util.fitViewInWorld(view_region, world_region.y_region)
        #self.scrolled_window.setYviewRegion(view_region, n)
        (b0, b1) = view_region
        self.scrolled_window.setYviewRegion(Region1D(b0, b1), n)
        if not self.updatingAspectRatio:
          self.after_idle(lambda: self.updateAspectRatio(doDraw=False))
        self.drawRow(n)
        self.drawSliceRow(n)

  def changedAxisSize(self, axisRegion):

    panel = axisRegion.axisPanel
    label = panel.label
    if (label not in ('x', 'y')):
      return

    (w, h, x, y) = self.windowPopup.getGeometry()
    size = axisRegion.size
    if (panel.spectrumWindowPane == self.windowPane):
      if (label == 'x'):
        col = panel.sortedAxisRegions().index(axisRegion)
        ww = self.scrolled_window.canvases[0][col].winfo_width()
        #self.scrolled_window.zeroColWeights(col=col)
        #self.scrolled_window.gridAll(col=col)
        #self.scrolled_window.zeroWeights(col=col)
        #self.setGeometry(w+size-ww, h, x, y)
        #for j in range(self.scrolled_window.nrows):
        #  self.scrolled_window.canvases[j][col].config(width=size)
        #self.scrolled_window.canvases[0][col].config(width=size)
        #self.after_idle(self.scrolled_window.oneColWeights)
        #self.after_idle(self.scrolled_window.gridAll)
        #print 'changedAxisSize1', panel.sortedAxisRegions().index(axisRegion), w, ww, size, w+size-ww
      elif (label == 'y'):
        row = panel.sortedAxisRegions().index(axisRegion)
        hh = self.scrolled_window.canvases[row][0].winfo_height()
        #self.scrolled_window.zeroRowWeights(row=row)
        #self.scrolled_window.gridAll(row=row)
        #self.scrolled_window.zeroWeights(row=row)
        #self.setGeometry(w, h+size-hh, x, y)
        #for i in range(self.scrolled_window.ncols):
        #  self.scrolled_window.canvases[row][i].config(height=size)
        #self.scrolled_window.canvases[row][0].config(height=size)
        #self.after_idle(self.scrolled_window.oneRowWeights)
        #self.after_idle(self.scrolled_window.gridAll)
        #print 'changedAxisSize2', panel.sortedAxisRegions().index(axisRegion), h, hh, size, h+size-hh
      #self.drawAll()

  def changedOrthogonalRegion(self, axisRegion):

    windowPane = self.windowPane
    axisPanel = axisRegion.axisPanel
    if axisPanel.spectrumWindowPane == windowPane:
      self.drawAllAfter()

      if axisPanel.axisType.isSampled:
        self.setPseudoState(axisPanel)
      else:
        axisRegions = axisPanel.sortedAxisRegions()
        if axisRegion not in axisRegions:
          return
        n = axisRegions.index(axisRegion)
        # Can we use axisRegion.serial?

        # check if corresponding x axis region is active, if so change selector
        stripAxis = windowPane.spectrumWindow.stripAxis
        stripAxisPanel = windowPane.findFirstAxisPanel(label=stripAxis)
        stripAxisRegions = stripAxisPanel.sortedAxisRegions()
        if n >= len(stripAxisRegions):
          return
        stripAxisRegion = stripAxisRegions[n]

        if len(axisRegions) == 1:
          stripAxisRegion.isActive = True

        if stripAxisRegion.isActive:
          (r0, r1) = Util.checkSwapRegion(axisRegion.region, axisPanel.axisUnit)
          axisPanel.region_selector.setViewRegion(r0, r1, do_callback=True)

  def changedAxisRegion(self, axisRegion):

    axisPanel = axisRegion.axisPanel
    #print 'WindowPopup.changedAxisRegion', axisPanel.label, axisRegion.region, axisPanel.sortedAxisRegions().index(axisRegion)

    if axisPanel.spectrumWindowPane != self.windowPane:
      return

    if axisPanel.label in ('x', 'y'):
      self.changedViewRegion(axisRegion)
    else:
      self.changedOrthogonalRegion(axisRegion)

  def changedIsActive(self, axisRegion):

    if axisRegion.axisPanel.spectrumWindowPane != self.windowPane:
      return

    self.drawAllAfter()

  # some axisRegion has been created or deleted
  def changedAxisRegions(self, axisRegion):

    axisPanel = axisRegion.axisPanel
    windowPane = axisPanel.spectrumWindowPane
    if windowPane != self.windowPane:
      return

    self.windowPopup.setStripsSelectorAfter()

  def changedPanelVisibility(self, panel):

    if panel.spectrumWindowPane is self.windowPane:
      if panel.label in ('x', 'y'):
        #print 'changedPanel', panel.label
        if panel.label == 'x':
          self.scrolled_window.setIsShownXScrollbar(panel.isVisible)
        else: # panel.label == 'y'
          self.scrolled_window.setIsShownYScrollbar(panel.isVisible)
      else:
        self.gridAll()

  def changedSliceVisibility(self, panel):

    if panel.spectrumWindowPane is self.windowPane:
      if panel.label == 'x':
        self.scrolled_window.setIsShownXSlice(panel.isVisible)
      elif panel.label == 'y':
        self.scrolled_window.setIsShownYSlice(panel.isVisible)
      else:
        self.gridAll()
      if (panel.isVisible):
        self.drawAllSlices()

  def changedSliceThickness(self, panel):

    if panel.spectrumWindowPane is self.windowPane:
      if panel.label == 'x':
        t = Util.greaterOrEqualEntry(panel.thickness, slice_height_entries)
        self.scrolled_window.setSizeXSlice(t)
      elif panel.label == 'y':
        t = Util.greaterOrEqualEntry(panel.thickness, slice_height_entries)
        self.scrolled_window.setSizeYSlice(t)
      #draw automatically happens
      if panel.isVisible:
        self.drawAllSlices()

  def minAxisPanelThickness(self, axisPanel):

    label = axisPanel.label

    min_thickness = None
    if axisPanel.axisUnit:
      unit = axisPanel.axisUnit.unit
      for view in self.getSpectrumViews():
        axisMapping = view.findFirstAxisMapping(label=label)
        if axisMapping:
          t = Util.convertPosition(1.0, getPrimaryDataDimRef(axisMapping.analysisDataDim.dataDim),
                                   toUnit=unit, relative=True)
          t = abs(t)
          if min_thickness is None:
            min_thickness = t
          else:
            min_thickness = min(min_thickness, t)

    if min_thickness is None:
      min_thickness = 1.0

    return min_thickness

  def getOrthogonalWorldRegion(self, axisPanel):

    # TBD: some assumptions in here that working with ppm...

    # TBD: this looks like it belongs in WindowDraw

    axisType = axisPanel.axisType
    label = axisPanel.label
    possViews = self.getSpectrumViews()
    views = []
    # try to protect against deleted or messed up views
    for view in possViews:
      if not view.isDeleted:
        try:
          axisMapping = view.findFirstAxisMapping(label=label)
          dataDim = axisMapping.analysisDataDim.dataDim
          if axisMapping and dataDim:
            views.append(view)
        except:
          pass

    # special case
    if axisType.isSampled:
      if hasattr(axisPanel, 'region_selector'):
        region_selector = axisPanel.region_selector
        region = (0, max(region_selector.numButtons, 1))
      else:
        region = (0, 10)
      if views:
        r1 = None
        for view in views:
          axisMapping = view.findFirstAxisMapping(label=label)
          dataDim = axisMapping.analysisDataDim.dataDim
          n = dataDim.numPoints
          if r1 is None:
            r1 = n
          else:
            r1 = max(r1, n)
        if r1 is not None:
          region = (0.5, r1+0.5)
      return region

    if not views:
      # if no views just use the default axisType world region
      region = Util.checkSwapRegion(axisType.region, axisType.findFirstAxisUnit(unit='ppm'))
      return region

    r0 = None
    r1 = None

    for view in views:
      axisMapping = view.findFirstAxisMapping(label=label)
      dataDim = axisMapping.analysisDataDim.dataDim
      dataDimRef = getPrimaryDataDimRef(dataDim)
      expDimRef = dataDimRef.expDimRef

      minAliasedFreq = expDimRef.minAliasedFreq
      if minAliasedFreq is None:
        minAliasedFreq = Util.convertPosition(float(dataDim.numPoints), dataDimRef, toUnit=axisType.findFirstAxisUnit(unit='ppm').unit)
      if r0 is None:
        r0 = minAliasedFreq
      else:
        r0 = min(r0, minAliasedFreq)

      maxAliasedFreq = expDimRef.maxAliasedFreq
      if maxAliasedFreq is None:
        maxAliasedFreq = Util.convertPosition(1.0, dataDimRef, toUnit=axisType.findFirstAxisUnit(unit='ppm').unit)
      if r1 is None:
        r1 = maxAliasedFreq
      else:
        r1 = max(r1, maxAliasedFreq)

    (r0, r1) = Util.checkSwapRegion((r0, r1), axisType.findFirstAxisUnit(unit='ppm'))

    return Region1D(r0, r1)

  def updateOrthogonalWorlds(self):

    for axisPanel in self.windowPane.sortedAxisPanels()[2:]:
      if hasattr(axisPanel, 'region_selector'):
        world_region = self.getOrthogonalWorldRegion(axisPanel)
        axisType = axisPanel.axisType
        if axisType.isSampled:
          axisPanel.region_selector.setNumButtons(int(world_region[1]))
        else:
          axisPanel.region_selector.setWorldRegion(world_region)

  def getPseudoState(self, axisPanel):

    axisRegions = axisPanel.sortedAxisRegions()
    (w0, w1) = self.getOrthogonalWorldRegion(axisPanel)
    nstates = int(w1)
    state = nstates * [False]
    flag = True
    for axisRegion in axisRegions:
      (v0, v1) = axisRegion.region
      v0 = max(0, int(v0)-1)
      v1 = min(int(v1), nstates)
      if v0 < v1:
        flag = False
        for v in range(v0, v1):
          state[v] = True

    if flag: # want at least one on
      state[0] = True

    return state

  def setPseudoState(self, axisPanel):

    state = self.getPseudoState(axisPanel)
    axisPanel.region_selector.setState(state)

  def createRegionSelector(self, axisPanel):

    if axisPanel.spectrumWindowPane != self.windowPane:
      return

    self.deleteRegionSelector(axisPanel)

    axisType = axisPanel.axisType
    (w0, w1) = world_region = self.getOrthogonalWorldRegion(axisPanel)
    callback = lambda view_region, axisPanel=axisPanel: \
                   self.regionSelectorCallback(view_region, axisPanel)

    if axisType.isSampled:
      numButtons = int(w1)
      #state = w0*[False] + (w1-w0+1)*[True]
      state = self.getPseudoState(axisPanel)
      axisPanel.region_selector = ButtonScrollbar(self, numButtons=numButtons,
                                                  state=state, label='Plane',
                                                  callback=callback)
    else:

      #if True:
      label = axisType.name
      (v0, v1) = Util.checkSwapRegion(axisPanel.findFirstAxisRegion().region, axisPanel.axisUnit)
      if (w0 < w1):
        if v0 < w0:
          d = v1 - v0
          v0 = w0
          v1 = min(v0+d, w1)
        elif v1 > w1:
          d = v1 - v0
          v1 = w1
          v0 = max(v1-d, w0)
      else:
        if v0 > w0:
          d = v1 - v0
          v0 = w0
          v1 = max(v0+d, w1)
        elif v1 < w1:
          d = v1 - v0
          v1 = w1
          v0 = min(v1-d, w0)
      view_region = Region1D(v0, v1)

      min_thickness = self.minAxisPanelThickness(axisPanel)
      #print 'createRegionSelector', label, world_region, view_region
      axisPanel.region_selector = RegionSelector(self, label=label,
          callback=callback, text_decimals=axisType.numDecimals,
          world_region=world_region, view_region=view_region,
          min_thickness=min_thickness)

  def deleteRegionSelector(self, axisPanel):

    if (hasattr(axisPanel, 'region_selector')):
      axisPanel.region_selector.destroy()
      del axisPanel.region_selector

  def createSliceCanvas(self, slicePanel):

    if slicePanel.spectrumWindowPane != self.windowPane:
      return

    self.deleteSliceCanvas(slicePanel)
    # TBD: more

  def deleteSliceCanvas(self, slicePanel):

    if hasattr(slicePanel, 'slice_canvas'):
      slicePanel.slice_canvas.destroy()
      del slicePanel.slice_canvas

  def gridAll(self):

    row = 0
    self.scrolled_window.grid(row=row, column=0, columnspan=2, sticky='nsew')
    row += 1

    for axisPanel in self.windowPane.sortedAxisPanels()[2:]:
      if axisPanel.isVisible:
        axisPanel.region_selector.grid(row=row, column=0, columnspan=2, sticky='ew')
        row += 1
        
      else:
        axisPanel.region_selector.grid_forget()

    # TBD: uncomment this (after modifying)
    #for slicePanel in self.windowPane.slicePanels[2:]:
    #  if (slicePanel.isVisible):
    #    slicePanel.region_selector.grid(row=row, column=0, sticky='ew')
    #    row = row + 1
    #  else:
    #    slicePanel.region_selector.grid_forget()

  def regionSelectorCallback(self, viewRegion, axisPanel):

    if not hasattr(axisPanel, 'region_selector'):
      return

    stripAxis = self.windowPane.spectrumWindow.stripAxis
    stripAxisRegions = self.windowPane.findFirstAxisPanel(label=stripAxis).sortedAxisRegions()
    if len(stripAxisRegions) == 1:
      stripAxisRegions[0].isActive = True

    if axisPanel.axisType.isSampled:
      numPlanes = len(viewRegion)
      ends   = []
      starts = [viewRegion[0]+1,]
      for i in range(1,numPlanes):
        this = viewRegion[i]
        prev = viewRegion[i-1]
        if this != prev+1:
          ends.append(prev+1)
          starts.append(this+1)

      ends.append(viewRegion[-1]+1)

      regions = [(starts[i],ends[i]) for i in range(len(starts))]

      axisRegions = axisPanel.sortedAxisRegions()
      numR = len(regions)
      numA = len(axisRegions)

      if numA > numR:
        for i in range(numR,numA):
          axisRegions[i].delete()

      for i in range(numR):
        if i < numA:
          axisRegions[i].region = regions[i]
        else:
          axisRegion = axisPanel.newAxisRegion(region=regions[i])

      axisRegions = axisPanel.sortedAxisRegions()

    else:
      viewRegion = Util.checkSwapRegion(viewRegion, axisPanel.axisUnit)
      axisRegions = axisPanel.sortedAxisRegions()

      for n in range(len(stripAxisRegions)):
        if (stripAxisRegions[n].isActive):
          try:
            axisRegions[n].region = viewRegion
          except:
            pass  # can fail when orthogonal axis being created

    #self.drawAll()

  def destroy(self):

    self.curateNotifiers(self.windowPopup.unregisterNotify)

    self.deleteHandlers()

    Frame.destroy(self)

  def curateNotifiers(self, notify):

    # AnalysisSpectrum
    
    clazz = 'ccpnmr.Analysis.AnalysisSpectrum'
    callback = self.changedAnalysisSpectrum
    for func in ('set', 'add', 'remove'):
      for attr in ('NegLevels', 'PosLevels', 'PosColors', 'NegColors'):
        notify(callback, clazz, func+attr)

    for attr in ('Rank', 'UseBoundingBox', 'SliceColor', 'UsePeakArrow', 'UsePrecalculated'):
      notify(callback, clazz, 'set'+attr)

    # AnalysisPeakList
    
    clazz = 'ccpnmr.Analysis.AnalysisPeakList'
    callback = self.changedAnalysisPeakList
    for attr in ('SymbolStyle', 'SymbolColor', 'TextColor'):
      notify(callback, clazz, 'set'+attr)
    
    #
    
    for func in ('setIsXSliceDrawn', 'setIsYSliceDrawn'):
      notify(self.drawAllAfter, 'ccpnmr.Analysis.SpectrumWindow', func)
    
    notify(self.changedAxisTypeRegion, 'ccpnmr.Analysis.AxisType', 'setRegion')
    notify(self.drawAllAfter, 'ccpnmr.Analysis.AxisType', 'setPeakSize')
    notify(self.changedAspectRatio, 'ccpnmr.Analysis.SpectrumWindowPane', 'setAspectRatio')
    notify(self.changedSliceRange, 'ccpnmr.Analysis.SpectrumWindowPane', 'setSliceRange')
    for func in ('setIsCanvasLabelShown', 'setIsCanvasMidpointShown'):
      notify(self.drawAllAfter, 'ccpnmr.Analysis.SpectrumWindow', func)
    notify(self.changedAxisRegion, 'ccpnmr.Analysis.AxisRegion', 'setRegion')
    notify(self.changedIsActive, 'ccpnmr.Analysis.AxisRegion', 'setIsActive')
    for func in ('__init__', 'delete'):
      notify(self.changedAxisRegions, 'ccpnmr.Analysis.AxisRegion', func)
    notify(self.changedAxisSize, 'ccpnmr.Analysis.AxisRegion', 'setSize')
    notify(self.initSpectrumWindowView, 'ccpnmr.Analysis.SpectrumWindowView', '__init__')
    notify(self.changedSpectrumWindowView, 'ccpnmr.Analysis.SpectrumWindowView', 'setIsPosVisible')
    notify(self.changedSpectrumWindowView, 'ccpnmr.Analysis.SpectrumWindowView', 'setIsNegVisible')
    notify(self.changedSpectrumWindowViewSlice, 'ccpnmr.Analysis.SpectrumWindowView', 'setIsSliceVisible')
    notify(self.changedSpectrumWindowViewSlice, 'ccpnmr.Analysis.SpectrumWindowView', 'setIsContourLineVisible')
    notify(self.deleteSpectrumWindowView, 'ccpnmr.Analysis.SpectrumWindowView', 'delete')
    notify(self.changedPanelVisibility, 'ccpnmr.Analysis.AxisPanel', 'setIsVisible')
    notify(self.initAxisPanel, 'ccpnmr.Analysis.AxisPanel', '__init__')
    notify(self.deleteAxisPanel, 'ccpnmr.Analysis.AxisPanel', 'delete')
    notify(self.createSliceCanvas, 'ccpnmr.Analysis.SlicePanel', '__init__')
    notify(self.changedSliceVisibility, 'ccpnmr.Analysis.SlicePanel', 'setIsVisible')
    notify(self.changedSliceThickness, 'ccpnmr.Analysis.SlicePanel', 'setThickness')
    notify(self.changedSpectrum, 'ccp.nmr.Nmr.DataSource', 'setScale')
    notify(self.changedPeakList, 'ccp.nmr.Nmr.PeakList', 'setPeakColor')
    notify(self.changedPeakList, 'ccp.nmr.Nmr.PeakList', 'setPeakSymbol')
    notify(self.changedWinPeakList, 'ccpnmr.Analysis.WindowPeakList', 'setIsSymbolDrawn')
    notify(self.changedWinPeakList, 'ccpnmr.Analysis.WindowPeakList', 'setIsAnnotationDrawn')
    notify(self.initPeak, 'ccp.nmr.Nmr.Peak', '__init__')
    notify(self.deletePeak, 'ccp.nmr.Nmr.Peak', 'delete')
    notify(self.changedPeak, 'ccp.nmr.Nmr.Peak', 'setAnnotation')
    notify(self.changedPeakDim, 'ccp.nmr.Nmr.PeakDim', 'setPosition')
    notify(self.changedPeakDim, 'ccp.nmr.Nmr.PeakDim', 'setAnnotation')
    notify(self.changedPeakDim, 'ccp.nmr.Nmr.PeakDim', 'setLineWidth')
    
    for func in ('setPeakDrawMethod', 'setPeakPixelSize', 'setPeakIntensityScale',
                 'setPeakVolumeScale', 'setPeakFindThickness', 'setGlobalContourScale',
                 'setDoMinimalAnnotations'):
      notify(self.drawAllAfter, 'ccpnmr.Analysis.AnalysisProject', func)
      
    notify(self.changedMaxMarks, 'ccpnmr.Analysis.AnalysisProject', 'setMaxMarks')
    notify(self.changedMaxRulers, 'ccpnmr.Analysis.AnalysisProject', 'setMaxRulers')
    notify(self.changedMark, 'ccpnmr.Analysis.Mark', '')
    notify(self.changedMark, 'ccpnmr.Analysis.MarkDim', '')
    notify(self.changedMark, 'ccpnmr.Analysis.MarkDim', '__init__')
    notify(self.changedMark, 'ccpnmr.Analysis.MarkDim', 'delete')
    notify(self.changedRuler, 'ccpnmr.Analysis.Ruler', '')
    notify(self.changedRuler, 'ccpnmr.Analysis.Ruler', '__init__')
    notify(self.changedRuler, 'ccpnmr.Analysis.Ruler', 'delete')
    notify(self.changedColor, 'ccpnmr.AnalysisProfile.ColorScheme', 'setColors')
    notify(self.changedBackground, 'ccpnmr.AnalysisProfile.AnalysisProfile', 'setBgColor')
    notify(self.changedReferencing, 'ccp.nmr.Nmr.ExpDimRef', '')
    notify(self.changedReferencing, 'ccp.nmr.Nmr.DataDimRef', '')
    notify(self.changedReferencing, 'ccp.nmr.Nmr.FreqDataDim', '')
    notify(self.initAxisMapping, 'ccpnmr.Analysis.AxisMapping', '__init__')

    for func in ('setIsBigEndian', 'setNByte', 'setNumberType',
                 'setHeaderSize', 'setHasBlockPadding', 'setBlockSizes',
                 'setPath', 'setNumPoints', 'setDataUrl'):
      notify(self.drawAllAfter, 'ccp.general.DataLocation.BlockedBinaryMatrix', func)
    
    for func in ('setNumberType', 'setHeaderSize', 'setPath', 'setNumPoints', 'setDataUrl'):
      notify(self.drawAllAfter, 'ccp.general.DataLocation.ShapeMatrix', func)
    notify(self.drawAllAfter, 'ccp.general.DataLocation.DataUrl', 'setUrl')

    notify(self.drawAllAfter, 'ccp.nmr.Nmr.PeakDim', 'setAppDataValue', keyword='textOffset')
    notify(self.drawAllAfter, 'ccpnmr.Analysis.SpectrumWindowView', 'setAppDataValue', keyword=VALUE_AXIS_OFFSET)
    notify(self.drawAllAfter, 'ccpnmr.Analysis.SpectrumWindowView', 'setAppDataValue', keyword=X_AXIS_OFFSET)

    ### TBD: v2
    ###for func in ('__init__', 'delete', 'setValue'):
    ###  for t in ('Int', 'Float', 'String', 'Boolean'):
    ###    clazz = 'memops.Implementation.AppData' + t
    ###    notify(self.changedAppData, clazz, func)
    for func in ('__init__', 'delete'):
      notify(self.changedStoredContour, 'ccpnmr.Analysis.StoredContour', func)

  # def _returnWidget(self, handlerReturn):
  #   return handlerReturn
  #
  # def _getWidgetHandler(self, handlerClass, widget, msg_queue, args):
  #   """Get the handler from the c code - separate thread test
  #   """
  #   widgetHandler = None
  #   try:
  #     widgetHandler = handlerClass and handlerClass(widget, *args)
  #
  #   except Exception as e:
  #     msg_queue.put(lambda val=str(e): self._returnWidget(val))
  #   else:
  #     msg_queue.put(lambda val=widgetHandler: self._returnWidget(val))
  #
  # def setupWidgetHandler(self, widget, isCanvas = True):
  #
  #   import threading
  #   from Queue import Queue
  #
  #   if not hasattr(widget, 'handler'):
  #
  #     msg_queue = Queue()
  #     thread_count = threading.active_count()
  #     t = threading.Thread(target=self._getWidgetHandler, args=(self.handlerClass, widget, msg_queue, self.handlerArgs))
  #     t.start()
  #     result = None
  #     while threading.active_count() > thread_count:
  #       event = msg_queue.get()
  #       result = event()
  #
  #     if isinstance(result, str):
  #       widget.handler = None
  #       print(result)
  #
  #     else:
  #       widget.handler = result
  #       if result:
  #         # widget.handler = self.handlerClass and self.handlerClass(widget, *self.handlerArgs)
  #         if (isCanvas):
  #           (j, i) = self.scrolled_window.getCanvasRowCol(widget)
  #           widget.xview = self.windowPane.findFirstAxisPanel(label='x').sortedAxisRegions()[i]
  #           widget.yview = self.windowPane.findFirstAxisPanel(label='y').sortedAxisRegions()[j]
  #           self.setWidgetHandlerBackground(widget)

  def setupWidgetHandler(self, widget, isCanvas = True):

    if not hasattr(widget, 'handler'):
      #print 'setupWidgetHandler'
      widget.handler = self.handlerClass and self.handlerClass(widget, *self.handlerArgs)
      if (isCanvas):
        (j, i) = self.scrolled_window.getCanvasRowCol(widget)
        widget.xview = self.windowPane.findFirstAxisPanel(label='x').sortedAxisRegions()[i]
        widget.yview = self.windowPane.findFirstAxisPanel(label='y').sortedAxisRegions()[j]
        self.setWidgetHandlerBackground(widget)

  def setHandlers(self):

    scrolledWindow = self.scrolled_window
    handlerClass   = self.handlerClass
    handlerArgs    = self.handlerArgs
    rows           = range(scrolledWindow.nrows)
    cols           = range(scrolledWindow.ncols)
    findFirstAxisPanel = self.windowPane.findFirstAxisPanel
    canvases = scrolledWindow.canvases
    xRegions = findFirstAxisPanel(label='x').sortedAxisRegions()
    yRegions = findFirstAxisPanel(label='y').sortedAxisRegions()
    
    for j in rows:
      for i in cols:
        canvas = canvases[j][i]
        canvas.handler = handlerClass and handlerClass(canvas, *handlerArgs)
        canvas.xview = xRegions[i]
        canvas.yview = yRegions[j]

    if self.hasValueAxis:
      return

    for i in cols:
      slice = scrolledWindow.xslices[i]
      slice.handler = handlerClass and handlerClass(slice, *handlerArgs)

    for j in rows:
      slice = scrolledWindow.yslices[j]
      slice.handler = handlerClass and handlerClass(slice, *handlerArgs)

  def deleteHandlers(self):

    # force C world clean-up

    #print 'deleteHandlers start'
    #for j in range(self.scrolled_window.nrows):
    #  for i in range(self.scrolled_window.ncols):
    #    canvas = self.scrolled_window.canvases[j][i]
    #    canvas.handler = None

    if not hasattr(self, 'scrolled_window'):
      return

    for j in range(self.scrolled_window.nrows):
      if j >= len(self.scrolled_window.canvases):
        print "Cleanup called on non-existent canvases"
        continue

      for i in range(self.scrolled_window.ncols):
        if i >= len(self.scrolled_window.canvases[j]):
          print "Cleanup called on non-existent canvases"
          continue

        canvas = self.scrolled_window.canvases[j][i]
        try:
          canvas.handler.delete()
        except:
          pass
        canvas.handler = None

    if (self.hasValueAxis):
      return

    for i in range(self.scrolled_window.ncols):
      slice = self.scrolled_window.xslices[i]
      try:
        slice.handler.delete()
      except:
        pass
      slice.handler = None

    for j in range(self.scrolled_window.nrows):
      slice = self.scrolled_window.yslices[j]
      try:
        slice.handler.delete()
      except:
        pass
      slice.handler = None

    #print 'deleteHandlers end'

  def setHandlersRow(self, row=-1):
  
    scrolledWindow = self.scrolled_window
    if row == -1:
      row = scrolledWindow.nrows - 1

    findFirstAxisPanel = self.windowPane.findFirstAxisPanel
    xRegions = findFirstAxisPanel(label='x').sortedAxisRegions()
    yRegions = findFirstAxisPanel(label='y').sortedAxisRegions()
    canvases = scrolledWindow.canvases
    
    for i in range(scrolledWindow.ncols):
      canvas = canvases[row][i]
      if not hasattr(canvas, 'handler'):
        canvas.handler = self.handlerClass and self.handlerClass(canvas, *self.handlerArgs)
      canvas.xview = xRegions[i]
      canvas.yview = yRegions[row]

    if self.hasValueAxis:
      return

    slice = scrolledWindow.yslices[row]
    if not hasattr(slice, 'handler'):
      slice.handler = self.handlerClass and self.handlerClass(slice, *self.handlerArgs)

  def setHandlersCol(self, col=-1):

    scrolledWindow = self.scrolled_window
    if col == -1:
      col = scrolledWindow.ncols - 1

    findFirstAxisPanel = self.windowPane.findFirstAxisPanel
    xRegions = findFirstAxisPanel(label='x').sortedAxisRegions()
    yRegions = findFirstAxisPanel(label='y').sortedAxisRegions()
    canvases = scrolledWindow.canvases
    
    #print 'setHandlersCol', col
    for j in range(scrolledWindow.nrows):
      canvas = canvases[j][col]
      if not hasattr(canvas, 'handler'):
        canvas.handler = self.handlerClass and self.handlerClass(canvas, *self.handlerArgs)
      canvas.xview = xRegions[col]
      canvas.yview = yRegions[j]

    if self.hasValueAxis:
      return

    slice = scrolledWindow.xslices[col]
    if not hasattr(slice, 'handler'):
      slice.handler = self.handlerClass and self.handlerClass(slice, *self.handlerArgs)

  def setWidgetHandlerBackground(self, widget):

    color = self.windowPopup.analysisProfile.bgColor
    if widget.handler:
      widget.handler.setBackground(hexToRgb(color))

  def changedBackground(self, *extra):

    self.setHandlerBackground()
    self.drawAllAfter()

  def setHandlerBackground(self):

    color = hexToRgb(self.windowPopup.analysisProfile.bgColor)
    scrolledWindow = self.scrolled_window
    rows = range(scrolledWindow.nrows)
    cols = range(scrolledWindow.ncols)
    setupWidgetHandler =  self.setupWidgetHandler
    canvases = scrolledWindow.canvases
    
    #print 'setHandlerBackground', self.windowPane.name, color

    self.scrolled_window.setBackground(self.windowPopup.analysisProfile.bgColor)

    for j in rows:
      for i in cols:
        canvas = canvases[j][i]
        setupWidgetHandler(canvas, isCanvas=True)
        if canvas.handler:
          canvas.handler.setBackground(color)

    if self.hasValueAxis:
      return

    for i in cols:
      slice = scrolledWindow.xslices[i]
      setupWidgetHandler(slice, isCanvas=False)
      if slice.handler:
        slice.handler.setBackground(color)

    for j in rows:
      slice = scrolledWindow.yslices[j]
      setupWidgetHandler(slice, isCanvas=False)
      if slice.handler:
        slice.handler.setBackground(color)

  def changedAnalysisSpectrum(self, analysisSpec):
  
    for view in analysisSpec.spectrumWindowViews:
      if view.spectrumWindowPane is self.windowPane:
        self.windowPopup.setSpectrumSelector()
        if self.isViewVisible(view):
          self.drawAllAfter()
        return
  
  def changedAnalysisPeakList(self, analysisPeakList):
  
    windowPane = self.windowPane
    for winPeakList in analysisPeakList.windowPeakLists:
      if winPeakList.spectrumWindowView.spectrumWindowPane is windowPane:
        if self.isWinPeakListDrawn(winPeakList):
          self.drawAllAfter()
        return

  def changedAspectRatio(self, windowPane, doDraw = True):

    if windowPane != self.windowPane:
      return

    self.updateAspectRatio(doDraw=doDraw)

  def updateAspectRatio(self, doDraw = True):

    if self.hasValueAxis:
      return

    scrolled_window = self.scrolled_window
    canvases = scrolled_window.canvases
    self.updatingAspectRatio = True
    checkAspectRatio = self.checkAspectRatio
    
    try:
      for j in range(scrolled_window.nrows):
        checkAspectRatio(canvases[j][0])

      for i in range(1,scrolled_window.ncols):
        checkAspectRatio(canvases[0][i])

    finally:
      self.updatingAspectRatio = False

    if doDraw:
      self.drawAllAfter()

  def changedSliceRange(self, windowPane):

    if windowPane is not self.windowPane:
      return

    if self.hasValueAxis:

      axisPanel = windowPane.findFirstAxisPanel(label='y')
      # TBD: this is not very good: just set some random axisRegion to sliceRange
      axisRegion = axisPanel.findFirstAxisRegion()
      axisRegion.region = windowPane.sliceRange

    else:
      self.drawAllSlices()

      windowPane = self.windowPane
      window = windowPane.spectrumWindow
      if (window.isXSliceDrawn or window.isYSliceDrawn) \
          and self.canvasXs is not None \
          and self.canvasYs is not None:
        scrolledWindow = self.scrolled_window
        nCols = scrolledWindow.ncols
        nRows = scrolledWindow.nrows
        canvases = scrolledWindow.canvases
        drawCanvasCrosshairs = self.drawCanvasCrosshairs
        canvasXs = self.canvasXs
        canvasYs = self.canvasYs
        for j in range(nRows):
          for i in range(nCols):
            drawCanvasCrosshairs(canvases[j][i], canvasXs, canvasYs, self)

  def setMaxExtent(self, canvas):

    worldRegion = canvas.parent.world_region
    (x0, x1) = worldRegion.x_region
    (y0, y1) = worldRegion.y_region

    if not self.hasValueAxis:
      r = self.windowPane.aspectRatio
      w = float(canvas.canvas_width)
      h = float(canvas.canvas_height)

      rr = abs(((y1-y0)/h) / ((x1-x0)/w))
      if (rr > r):
        d = y0 + y1
        e = r * (x1-x0) * h / w
        y0 = 0.5 * (d - e)
        y1 = 0.5 * (d + e)
      else:
        d = x0 + x1
        e = (y1-y0) * w / (h * r)
        x0 = 0.5 * (d - e)
        x1 = 0.5 * (d + e)

    canvas.xmax_extent = abs(x1-x0)
    canvas.ymax_extent = abs(y1-y0)

    #print 'setWorldRegion0', canvas.parent.world_region.x_region, canvas.xmax_extent
    #print 'setWorldRegion1', canvas.parent.world_region.y_region, canvas.ymax_extent

  def checkAspectRatio(self, canvas):

    if self.hasValueAxis:
      self.setMaxExtent(canvas)
      return True

    windowPane = self.windowPane
    window = windowPane.spectrumWindow
    r = windowPane.aspectRatio
    try:
      w = float(canvas.canvas_width)
      h = float(canvas.canvas_height)
    except:
      return

    (x0, x1) = canvas.xview_region
    (y0, y1) = canvas.yview_region
    #print 'checkAspectRatio0', x0, x1, y0, y1, w, h, r

    (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
    xRegions = windowPane.findFirstAxisPanel(label='x').sortedAxisRegions()
    yRegions = windowPane.findFirstAxisPanel(label='y').sortedAxisRegions()
    xgroup = xRegions[col].axisRegionGroup
    ygroup = yRegions[row].axisRegionGroup

    if xgroup is not None and ygroup is not None: # ignore aspect ratio
      return False

    #print 'WindowPopup.checkAspectRatio0', row, col, x0, x1, y0, y1, xgroup, ygroup, self.scrolled_window.nrows, self.scrolled_window.ncols

    if window.stripAxis == 'y': # self.((not ygroup) and (self.scrolled_window.ncols == 1)): # resize y direction

      d = y0 + y1
      e = r * (x1-x0) * h / w
      yy0 = 0.5 * (d - e)
      yy1 = 0.5 * (d + e)
      if (yy0 == y0) and (yy1 == y1):
        return False

      #print 'checkAspectRatio1', yy0, yy1, row, col
      canvas.yview_region.set(yy0, yy1)
      axisRegion = yRegions[row]
      #print 'WindowPopup.checkAspectRatio1', (min(yy0, yy1), max(yy0, yy1))
      axisRegion.region = (min(yy0, yy1), max(yy0, yy1))

    else: # resize x direction

      d = x0 + x1
      e = (y1-y0) * w / (r * h)
      xx0 = 0.5 * (d - e)
      xx1 = 0.5 * (d + e)
      if (xx0 == x0) and (xx1 == x1):
        return False

      #print 'checkAspectRatio2', xx0, xx1, row, col
      canvas.xview_region.set(xx0, xx1)
      axisRegion = xRegions[col]
      #print 'WindowPopup.checkAspectRatio2', (min(xx0, xx1), max(xx0, xx1))
      axisRegion.region = (min(xx0, xx1), max(xx0, xx1))

    self.setMaxExtent(canvas)
    # not sure why after_idle required but otherwise get blank
    #canvas.after_idle(lambda self=self, canvas=canvas: self.scrolled_window.updateView(canvas))
    self.scrolled_window.updateView(canvas)

    return True

  def focusIn(self, event):

    if self is not event.widget:
      return
      
    #print 'focusIn', self.windowPane.name
    #self.scrolled_window.oneWeights()

  def resize(self, event):

    if self.waitResize:
      return

    #print 'in resize', self.windowPane.name
    canvas = event.widget
    width  = event.width
    height = event.height
    findFirstAxisPanel  = self.windowPane.findFirstAxisPanel
    
    (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
    #print 'resize', row, col, width, height, canvas.winfo_width(), canvas.winfo_height()
    findFirstAxisPanel(label='x').sortedAxisRegions()[col].size = width
    findFirstAxisPanel(label='y').sortedAxisRegions()[row].size = height

    self.setupWidgetHandler(canvas, isCanvas=True)

    canvas.canvas_width = width
    canvas.canvas_height = height

    if canvas.handler:
      canvas.handler.resize(width, height)
    redrawn = self.checkAspectRatio(canvas)

    ### do you always get expose events after resize in any case??
    if (not redrawn)  and (not self.waitDraw):
      # not sure why after_idle required but otherwise get blank
      #canvas.after_idle(lambda self=self, canvas=canvas: self.drawCanvas(canvas))
      self.drawCanvas(canvas)

  def expose(self, event):

    #print 'in expose', self.windowPane.name
    #print 'in expose', self.windowPane.name, event.x, event.y, event.width, event.height
    #print '***', event.state, event.num, '***'

    if not self.waitDraw:
      canvas = event.widget
      self.drawCanvas(canvas)
      
    #canvas.after_idle(lambda: self.drawCanvas(canvas))

  def showMotion(self, x, y, canvas = None):

    windowPane = self.windowPane
    stripAxis = windowPane.spectrumWindow.stripAxis
    getCanvasRowCol = self.scrolled_window.getCanvasRowCol
    
    windowPane.findFirstAxisPanel(label='x').pointLocation = x
    windowPane.findFirstAxisPanel(label='y').pointLocation = y
    for axisPanel in windowPane.sortedAxisPanels()[2:]:
      if canvas:
        try:
          (row, col) = getCanvasRowCol(canvas)
        except:
          return  # TBD: is this the right strategy??
        axisRegions = axisPanel.sortedAxisRegions()
        if axisPanel.axisType.isSampled:
          # TBD: does the below make the most sense???
          axisRegion = axisRegions[0]
        elif stripAxis == 'x':
          axisRegion = axisRegions[col]
        else:
          axisRegion = axisRegions[row]
        axisPanel.pointLocation = self.calcAxisMidpoint(axisRegion)
      else:
        axisPanel.pointLocation = None

    self.setLocationLabel(x, y)

    n = len(windowPane.axisPanels)
    typeLocation = n * [0]
    typeLocation = [ (ap.panelType, ap.pointLocation) for ap in windowPane.sortedAxisPanels() ]

    # extra multiple-quantum crosshairs

    if x is not None and y is not None: # if over 1D slice pane than can be None
      xaxisPanel = windowPane.findFirstAxisPanel(label='x')
      xaxisType = xaxisPanel.axisType
      yaxisPanel = windowPane.findFirstAxisPanel(label='y')
      yaxisType = yaxisPanel.axisType
      if xaxisType.isotopeCodes == yaxisType.isotopeCodes:
        if xaxisType.measurementType == 'MQShift' and yaxisType.measurementType == 'Shift':
          typeLocation.append((yaxisPanel.panelType, xaxisPanel.pointLocation-yaxisPanel.pointLocation))
        elif xaxisType.measurementType == 'Shift' and yaxisType.measurementType == 'MQShift':
          typeLocation.append((xaxisPanel.panelType, yaxisPanel.pointLocation-xaxisPanel.pointLocation))

    self.topPopup.drawCrosshairs(typeLocation, self)

  def motion(self, event):

    canvas = event.widget
    
    if not (hasattr(canvas, 'canvas_width') and 
            hasattr(canvas, 'canvas_height')):
      print ('WARNING, canvas called when lacking width or height')
      return
    
    w = canvas.canvas_width
    h = canvas.canvas_height
    x = event.x
    y = h - 1 - event.y
    (x0, x1) = canvas.xview_region
    (y0, y1) = canvas.yview_region
    r = (x + 0.5) / w
    s = (y + 0.5) / h
    x = x0 * (1 - r) + x1 * r
    y = y0 * (1 - s) + y1 * s

    self.showMotion(x, y, canvas)
    if not self.hasValueAxis:
      self.drawAllSlices()

  def setLocationLabel(self, x = None, y = None):

    nonetext = 'undef'
    windowPopup = self.windowPopup
    analysisProject = windowPopup.window.analysisProject
    unit = windowPopup.getPositionUnit()
    # TBD: usage of below assumes window is using ppm
    if unit == windowPopup.UNIT_HZ:
      (sfx, sfy) = windowPopup.getWindowXySf()
      view = windowPopup.getPositionView()
    else:
      sfx = sfy = 1.0
      view = None
    if x is None:
      xtext = nonetext
    else:
      # TBD: below doesn't really make sense because independent of unit
      xdecimals = self.windowPane.findFirstAxisPanel(label='x').axisType.numDecimals
      if self.deltaPositions:
        x -= self.deltaPositions[0]
      x *= sfx
      xtext = formatDecimals(x, decimals=xdecimals)

    if y is None:
      ytext = nonetext
    else:
      ydecimals = self.windowPane.findFirstAxisPanel(label='y').axisType.numDecimals
      if self.deltaPositions:
        y -= self.deltaPositions[1]
      y *= sfy
      ytext = formatDecimals(y, decimals=ydecimals)

    if view:
      tt = '(sf=%.0f,%.0f)' % (sfx, sfy)
    else:
      tt = self.xyName
    text = '%6s, %6s %s' % (xtext,ytext,tt)
    windowPopup.setLocationLabel(text)

  def drawCanvasCrosshairs(self, canvas, xs, ys, originatingWindowFrame, color=(1, 0, 0)):

    scrolled_window = self.scrolled_window
    if scrolled_window.isDeleting:
      return

    # this can be called before widget has been set up properly
    if (not hasattr(canvas, 'handler')):
      return

    if self.waitDraw:
      return

    handler = canvas.handler
    if not handler:
      return  # this can happen if no C code

    handler.makeCurrent()

    if hasattr(canvas, 'doubleBuffer') and not canvas.doubleBuffer:
      doubleBuffer = False
    else:
      doubleBuffer = True

    self.canvasXs = xs
    self.canvasYs = ys
    windowPane = self.windowPane 
    window = windowPane.spectrumWindow

    (x0, x1) = canvas.xview_region
    (y0, y1) = canvas.yview_region

    handler.mapRanges(0, 0, 1, 1, 0, 0, 1, 1)

    handler.startXor()
    handler.setColor(self.handlerXorColor)
    # above in xor mode means that below will appear in inverse colors

    # undo previous crosshairs if in single buffer mode

    if not doubleBuffer and hasattr(canvas, 'prev_cross'):
      # TBD: will this work for extra multiple-quantum crosshairs??
      (prev_xs, prev_x0, prev_x1, prev_ys, prev_y0, prev_y1) = canvas.prev_cross

      if prev_x0 == x0 and prev_x1 == x1 and prev_y0 == y0 and prev_y1 == y1:
        for x in prev_xs:
          x = float(x - prev_x0) / (prev_x1 - prev_x0)
          handler.drawLine(x, 0, x, 1)

        for y in prev_ys:
          y = float(y - prev_y0) / (prev_y1 - prev_y0)
          handler.drawLine(0, y, 1, y)

    # ordinary crosshairs

    # NOTE:ED - use the foreground colour here, but check that Xor works
    # wrong place
    # project = window.root
    # profile = Util.getAnalysisProfile(project)
    # color = inverseGrey(profile.bgColor)
    # handler.setColor(hexToRgb(color))
    fg = self.windowPopup.analysisProfile.fgColor
    # handler.setColor(hexToRgb(fg))
    bg = self.windowPopup.analysisProfile.bgColor
    if isMacOS():
      # only Tk
      handler.setColor(hexToRgb(fg))

    else:
      handler.setColor(hexToRgb(hexXor(fg, bg)))

    for x in xs:
      x = float(x - x0) / (x1 - x0)
      handler.drawLine(x, 0, x, 1)

    for y in ys:
      y = float(y - y0) / (y1 - y0)
      handler.drawLine(0, y, 1, y)

    xaxisPanel = windowPane.findFirstAxisPanel(label='x')
    xaxisType = xaxisPanel.axisType
    yaxisPanel = windowPane.findFirstAxisPanel(label='y')
    yaxisType = yaxisPanel.axisType

    if originatingWindowFrame == self:
      pass
    elif xs:
      xaxisPanel.pointLocation = xs[0]
    else:
      xaxisPanel.pointLocation = None

    if originatingWindowFrame == self:
      pass
    elif ys:
      yaxisPanel.pointLocation = ys[0]
    else:
      yaxisPanel.pointLocation = None

    for axisPanel in windowPane.sortedAxisPanels()[2:]:
      if (canvas):
        (row, col) = scrolled_window.getCanvasRowCol(canvas)
        axisRegions = axisPanel.sortedAxisRegions()
        if axisPanel.axisType.isSampled:
          axisRegion = axisRegions[0]
        elif window.stripAxis == 'x':
          axisRegion = axisRegions[col]
        else:
          axisRegion = axisRegions[row]
        axisPanel.pointLocation = self.calcAxisMidpoint(axisRegion)
      else:
        axisPanel.pointLocation = None

    (row, col) = scrolled_window.getCanvasRowCol(canvas)
    axisUnit = canvas.xview.axisPanel.axisType.findFirstAxisUnit(unit='ppm')
    (b0, b1) = xsliceRange = ysliceRange = windowPane.sliceRange
    unit = axisUnit.unit

    convertPosition = Util.convertPosition
    for view in self.getActiveSliceViews():
      if window.isXSliceDrawn:
        if yaxisPanel.pointLocation:
          a1 = (b1 - b0) * float(y1 - yaxisPanel.pointLocation) / (y1 - y0)
          a0 = a1 + b0 - b1
          ysliceRange = (a0, a1)
        axisMapping = view.findFirstAxisMapping(label='x')
        if not axisMapping: # can happen if change dimMapping
          continue
        ##if not hasattr(axisMapping, 'view_region'):
        ##  continue
        self.determineDimRange(view, axisMapping, row, col)
        dataDimRef = getPrimaryDataDimRef(axisMapping.analysisDataDim.dataDim)
        (r0, r1) = axisMapping.view_region
        # above in points, starting from 0, not 1; convert to ppm
        r0 = convertPosition(r0+1, dataDimRef, toUnit=unit)
        r1 = convertPosition(r1+1, dataDimRef, toUnit=unit)
        s = FakeSlice()
        s.col = col
        s.label = 'x'
        s.handler = handler
        s.orient = 'horizontal'
        s.view_region = (r0, r1)
        self.drawViewSlice(s, view, ysliceRange)

      if window.isYSliceDrawn:
        if xaxisPanel.pointLocation:
          a1 = (b1 - b0) * float(xaxisPanel.pointLocation - x0) / (x1 - x0)
          a0 = a1 + b0 - b1
          xsliceRange = (a0, a1)
        axisMapping = view.findFirstAxisMapping(label='y')
        if not axisMapping: # can happen if change dimMapping
          continue
        ##if not hasattr(axisMapping, 'view_region'):
        ##  continue
        self.determineDimRange(view, axisMapping, row, col)
        dataDimRef = getPrimaryDataDimRef(axisMapping.analysisDataDim.dataDim)
        (r0, r1) = axisMapping.view_region
        # above in points, starting from 0, not 1; convert to ppm
        r0 = convertPosition(r0+1, dataDimRef, toUnit=unit)
        r1 = convertPosition(r1+1, dataDimRef, toUnit=unit)
        s = FakeSlice()
        s.row = row
        s.label = 'y'
        s.handler = handler
        s.orient = 'vertical'
        s.view_region = (r0, r1)
        self.drawViewSlice(s, view, xsliceRange)

    # below needed for OpenGL implementation
    # otherwise xor of crosshairs messes up
    ###handler.mapRanges(0, 0, 1, 1, 0, 0, 1, 1)

    handler.finishXor()

    if not doubleBuffer:
      handler.flush()
      canvas.prev_cross = (xs, x0, x1, ys, y0, y1)

  def drawCrosshairs(self, typeLocation, originatingWindowFrame):

    scrolledWindow = self.scrolled_window
    windowPane = self.windowPane

    if scrolledWindow.isDeleting:
      return

    if self.windowPopup.state() != 'normal':
      return

    axisPanels = windowPane.sortedAxisPanels()
    xs = []
    ys = []
    for n in (0,1):
      axisPanel = axisPanels[n]
      
      for panelType, pointLocation in typeLocation:
        if panelType is axisPanel.panelType:
          if n == 0:
            if pointLocation is not None:
              xs.append(float(pointLocation))
          else:
            if pointLocation is not None:
              ys.append(float(pointLocation))

    self.crosshairLocations = (xs, ys)

    if originatingWindowFrame != self:

      if xs:
        x = xs[0]
      else:
        x = None

      if ys:
        y = ys[0]
      else:
        y = None

      self.setLocationLabel(x, y)

      # extra multiple-quantum crosshairs
      xaxisPanel = windowPane.findFirstAxisPanel(label='x')
      xaxisType = xaxisPanel.axisType
      yaxisPanel = windowPane.findFirstAxisPanel(label='y')
      yaxisType = yaxisPanel.axisType
      if xaxisType.isotopeCodes == yaxisType.isotopeCodes:
        if xaxisType.measurementType == 'MQShift' and yaxisType.measurementType == 'Shift':
          if len(ys) == 2:
            x = ys[0] + ys[1]
            xs.append(x)
        elif xaxisType.measurementType == 'Shift' and yaxisType.measurementType == 'MQShift':
          if len(xs) == 2:
            y = xs[0] + xs[1]
            ys.append(y)

    if not xs and not ys:
      return
    
    rows = range(scrolledWindow.nrows)
    cols = range(scrolledWindow.ncols)
    
    drawCanvasCrosshairs = self.drawCanvasCrosshairs
    canvases = scrolledWindow.canvases
    for j in rows:
      for i in cols:
        drawCanvasCrosshairs(canvases[j][i], xs, ys, originatingWindowFrame)

    for i in cols:
      ticks = scrolledWindow.xticks[i]
      ticks.drawCrosshairs(xs)

    for j in rows:
      ticks = scrolledWindow.yticks[j]
      ticks.drawCrosshairs(ys)

  def leave(self, event):

    axisPanels = self.windowPane.sortedAxisPanels()
    for i in range(2):
      axisPanels[i].pointLocation = None

    self.topPopup.endCrosshair()
    self.crosshairLocations = []

  def endCanvasCrosshair(self, canvas):

    if self.scrolled_window.isDeleting:
      return

    # this can be called before widget has been set up properly
    if not hasattr(canvas, 'handler'):
      return

    handler = canvas.handler
    if not handler:
      return
    handler.makeCurrent()
    handler.clearXor()

  def endSliceCrosshair(self, slice):

    if self.scrolled_window.isDeleting:
      return

    # this can be called before widget has been set up properly
    if not hasattr(slice, 'handler'):
      return

    handler = slice.handler
    if not handler:
      return
    handler.makeCurrent()
    handler.clearXor()

  def endCrosshair(self):

    scrolledWindow = self.scrolled_window
    if scrolledWindow.isDeleting:
      return

    if self.windowPopup.state() != 'normal':
      return

    rows = range(scrolledWindow.nrows)
    cols = range(scrolledWindow.ncols)
    canvases = scrolledWindow.canvases

    for j in rows:
      for i in cols:
        self.endCanvasCrosshair(canvases[j][i])

    for i in cols:
      ticks = scrolledWindow.xticks[i]
      ticks.clearCrosshairs()

    for j in rows:
      ticks = scrolledWindow.yticks[j]
      ticks.clearCrosshairs()

    if self.hasValueAxis:
      return

    endSliceCrosshair = self.endSliceCrosshair
    for i in cols:
      endSliceCrosshair(scrolledWindow.xslices[i])

    for j in rows:
      endSliceCrosshair(scrolledWindow.yslices[j])

  def sliceResize(self, event):

    #print 'in sliceResize'
    slice = event.widget
    width = event.width
    height = event.height

    self.setupWidgetHandler(slice, isCanvas=False)

    slice.slice_width = width
    slice.slice_height = height
    if slice.handler:
      slice.handler.resize(width, height)

  def sliceExpose(self, event):

    #print 'in sliceExpose'
    slice = event.widget
    self.drawSlice(slice)

  def sliceMotion(self, event):

    slice = event.widget

    if (slice.orient == Tkinter.HORIZONTAL):
      w = slice.slice_width
      x = event.x
      (x0, x1) = slice.view_region
      r = (x + 0.5) / w
      x = x0 * (1 - r) + x1 * r
      y = None
    else:
      h = slice.slice_height
      y = h - 1 - event.y
      (y0, y1) = slice.view_region
      s = (y + 0.5) / h
      y = y0 * (1 - s) + y1 * s
      x = None

    self.showMotion(x, y)
    #self.drawSliceCrosshair(slice, x=x, y=y)
    self.drawAllSlices()

  def viewSet(self, canvas, xview_region, yview_region):

    #print 'viewSet1', self.windowPane.name, xview_region[0], xview_region[1], yview_region[0], yview_region[1]

    # first time viewSet called setHandlers has not yet been called so xview, yview not set
    if (hasattr(canvas, 'xview')):

      #print 'viewSet2', self.windowPane.name
      axisUnit = canvas.xview.axisPanel.axisType.findFirstAxisUnit(unit='ppm')
      canvas.xview.region = Util.checkSwapRegion(xview_region, axisUnit)
      #print 'viewSet3', canvas.xview.region[0], canvas.xview.region[1]

      axisUnit = canvas.yview.axisPanel.axisType.findFirstAxisUnit(unit='ppm')
      canvas.yview.region = Util.checkSwapRegion(yview_region, axisUnit)
      #print 'viewSet4', canvas.yview.region[0], canvas.yview.region[1]

    else:

      pass # TBD: anything one should do?
      #print 'viewSet5', self.windowPane.name

    #print 'viewSet6', self.windowPane.name
    # TBD: is this needed, and when?
    #self.drawCanvas(canvas)
    #print 'viewSet7', self.windowPane.name

  def doCreateHorizontalRuler(self, canvas, x, y, doSidebands=False):

    (a, b, x, y) = self.calcWorldCoord(canvas, x, y)
    createRuler(b, self.windowPane.findFirstAxisPanel(label='y').panelType)
    #print 'doCreateHorizontalRuler', len(self.windowPane.root.rulers)

    if doSidebands and not self.hasValueAxis:
      # extra sideband rulers
      # TBD: assume for now that have ppm
      windowPane = self.windowPane 
      yaxisPanel = windowPane.findFirstAxisPanel(label='y')
      axisUnit = yaxisPanel.axisUnit
      if axisUnit.unit == 'ppm':
        (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
        if windowPane.spectrumWindow.stripAxis == 'x':
          strip = col
        else:
          strip = row
        positions = self.findPosition(a, b, strip, returnDict=False)
        views = self.getActiveSpectrumViews()
        for view in views:
          analysisSpectrum = view.analysisSpectrum
          spectrum = analysisSpectrum.dataSource
          experiment = spectrum.experiment
          spinningRate = experiment.spinningRate
          if spinningRate:
            dataDim = view.findFirstAxisMapping(label='y').analysisDataDim.dataDim
            dataDimRef = getPrimaryDataDimRef(dataDim)
            expDimRef = dataDimRef.expDimRef
            posn = positions[1]
            spinningRate /= expDimRef.sf
            minAliasedFreq = expDimRef.minAliasedFreq
            maxAliasedFreq = expDimRef.maxAliasedFreq
            if minAliasedFreq is None:
              minAliasedFreq = Util.convertPosition(float(dataDim.numPoints), dataDimRef, toUnit=axisUnit.unit)
            if maxAliasedFreq is None:
              maxAliasedFreq = Util.convertPosition(1.0, dataDimRef, toUnit=axisUnit.unit)
            maxRange = int(math.floor((maxAliasedFreq-posn)/spinningRate))
            minRange = int(math.ceil((minAliasedFreq-posn)/spinningRate))
            yrng = range(minRange, maxRange+1)
            for j in yrng:
              posn = positions[1] + j*spinningRate
              if j:  # skip if 0 because then ordinary ruler
                createRuler(posn, yaxisPanel.panelType, remove=False)

  def doCreateVerticalRuler(self, canvas, x, y, doSidebands=False):

    (a, b, x, y) = self.calcWorldCoord(canvas, x, y)
    createRuler(a, self.windowPane.findFirstAxisPanel(label='x').panelType)

    if doSidebands:
      # extra sideband rulers
      # TBD: assume for now that have ppm
      windowPane = self.windowPane 
      xaxisPanel = windowPane.findFirstAxisPanel(label='x')
      axisUnit = xaxisPanel.axisUnit
      if axisUnit.unit == 'ppm':
        (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
        if windowPane.spectrumWindow.stripAxis == 'x':
          strip = col
        else:
          strip = row
        positions = self.findPosition(a, b, strip, returnDict=False)
        views = self.getActiveSpectrumViews()
        for view in views:
          analysisSpectrum = view.analysisSpectrum
          spectrum = analysisSpectrum.dataSource
          experiment = spectrum.experiment
          spinningRate = experiment.spinningRate
          if spinningRate:
            dataDim = view.findFirstAxisMapping(label='x').analysisDataDim.dataDim
            dataDimRef = getPrimaryDataDimRef(dataDim)
            expDimRef = dataDimRef.expDimRef
            posn = positions[0]
            spinningRate /= expDimRef.sf
            minAliasedFreq = expDimRef.minAliasedFreq
            maxAliasedFreq = expDimRef.maxAliasedFreq
            if minAliasedFreq is None:
              minAliasedFreq = Util.convertPosition(float(dataDim.numPoints), dataDimRef, toUnit=axisUnit.unit)
            if maxAliasedFreq is None:
              maxAliasedFreq = Util.convertPosition(1.0, dataDimRef, toUnit=axisUnit.unit)
            maxRange = int(math.floor((maxAliasedFreq-posn)/spinningRate))
            minRange = int(math.ceil((minAliasedFreq-posn)/spinningRate))
            xrng = range(minRange, maxRange+1)
            for j in xrng:
              posn = positions[0] + j*spinningRate
              if j:  # skip if 0 because then ordinary ruler
                createRuler(posn, xaxisPanel.panelType, remove=False)

  def doCreateMark(self, canvas, x, y, doSidebands=False):
    
    getPanel = self.windowPane.findFirstAxisPanel
    (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
    xAxisRegion = getPanel(label='x').sortedAxisRegions()[col]
    yAxisRegion = getPanel(label='y').sortedAxisRegions()[row]
    (a, b, x, y) = self.calcWorldCoord(canvas, x, y)

    if self.windowPane.spectrumWindow.stripAxis == 'x':
      strip = col
    else:
      strip = row

    position_region = self.findPositionRegion(a, b, strip)

    views = self.getActiveSpectrumViews()
    peak  = None
    ydim  = -1

    params = getPeakFindParams(self.windowPane.root)
    thickness = params['thickness'] + self.extraPad

    for view in views:
      spectrum         = view.analysisSpectrum.dataSource
      (xscale, yscale) = self.getPeakScale(view, xAxisRegion, yAxisRegion)
      spectrum_region  = spectrum.numDim * [0]
      axisPanels       = spectrum.numDim * [0]
      
      for axisMapping in view.axisMappings:
        dataDim = axisMapping.analysisDataDim.dataDim
        dim     = dataDim.dim - 1
        label   = axisMapping.label
        
        if not self.isInAllowedRegion(dataDim, label, position_region[label]):
          break
          
        r = position_region[label]
        if label == 'x':
          xdim = dim
        elif label == 'y':
          ydim = dim
        else:
          (r0, r1) = r
          r = (r0-thickness, r1+thickness)
          
        axisPanels[dim] = getPanel(label=label)
        spectrum_region[dim] = self.convertPositionRegion(r, axisPanels[dim], dataDim)
        
      else:
        peakList = spectrum.activePeakList
        
        if peakList:
          result = self.findNearbyPeak(peakList, spectrum_region, xdim, ydim, xscale, yscale)
          
          if result:
            (peak, d2Min) = result
            break

    # extra multiple-quantum mark

    windowPane = self.windowPane 
    xaxisPanel = windowPane.findFirstAxisPanel(label='x')
    xaxisType = xaxisPanel.axisType
    yaxisPanel = windowPane.findFirstAxisPanel(label='y')
    yaxisType = yaxisPanel.axisType

    remove = True
    if xaxisType.isotopeCodes == yaxisType.isotopeCodes:
      analysisProject = windowPane.topObject
      if xaxisType.measurementType == 'MQShift' and yaxisType.measurementType == 'Shift':
        positions = self.findPosition(a, b, strip, returnDict=False)
        positions[1] = positions[0] - positions[1]
        axisTypes = [axisPanel.axisType for axisPanel in windowPane.sortedAxisPanels()]
        remove = False
        createNonPeakMark(positions, axisTypes, remove=remove)

      elif xaxisType.measurementType == 'Shift' and yaxisType.measurementType == 'MQShift':
        positions = self.findPosition(a, b, strip, returnDict=False)
        positions[0] = positions[1] - positions[0]
        axisTypes = [axisPanel.axisType for axisPanel in windowPane.sortedAxisPanels()]
        remove = False
        createNonPeakMark(positions, axisTypes, remove=remove)

      elif xaxisType.measurementType == 'Shift' and yaxisType.measurementType == 'Shift':
        # Dangerous: relies on MQ axis types defined in core.Util.defaultAxisTypes
        name = 'DQ' + xaxisType.isotopeCodes[0]
        axisType = analysisProject.findFirstAxisType(name=name)
        if axisType:
          # create mark for MQ axes
          positions = self.findPosition(a, b, strip, returnDict=False)
          position = positions[0] + positions[1]
          positions = (position,)
          axisTypes = (axisType,)
          remove = False
          createNonPeakMark(positions, axisTypes, remove=remove)

    if doSidebands:
      # extra sideband marks
      # TBD: assume for now that have ppm
      axisUnit = xaxisPanel.axisUnit
      if self.hasValueAxis:
        if axisUnit.unit == 'ppm':
          axisTypes = [axisPanel.axisType for axisPanel in windowPane.sortedAxisPanels()]
          positions = self.findPosition(a, b, strip, returnDict=False)
          for view in views:
            analysisSpectrum = view.analysisSpectrum
            spectrum = analysisSpectrum.dataSource
            experiment = spectrum.experiment
            spinningRate = experiment.spinningRate
            if spinningRate:
              dataDim = view.findFirstAxisMapping(label='x').analysisDataDim.dataDim
              dataDimRef = getPrimaryDataDimRef(dataDim)
              expDimRef = dataDimRef.expDimRef
              posn = positions[0]
              spinningRate /= expDimRef.sf
              minAliasedFreq = expDimRef.minAliasedFreq
              maxAliasedFreq = expDimRef.maxAliasedFreq
              if minAliasedFreq is None:
                minAliasedFreq = Util.convertPosition(float(dataDim.numPoints), dataDimRef, toUnit=axisUnit.unit)
              if maxAliasedFreq is None:
                maxAliasedFreq = Util.convertPosition(1.0, dataDimRef, toUnit=axisUnit.unit)
              maxRange = int(math.floor((maxAliasedFreq-posn)/spinningRate))
              minRange = int(math.ceil((minAliasedFreq-posn)/spinningRate))
              xrng = range(minRange, maxRange+1)
              posns = positions[:]
              for j in xrng:
                if j:  # skip if 0 because then ordinary mark
                  posns[0] = positions[0] + j*spinningRate
                  remove = False
                  createNonPeakMark(posns, axisTypes, remove=remove)
      else: # not self.hasValueAxis:
        if axisUnit.unit == yaxisPanel.axisUnit.unit == 'ppm':
          axisTypes = [axisPanel.axisType for axisPanel in windowPane.sortedAxisPanels()]
          positions = self.findPosition(a, b, strip, returnDict=False)
          for view in views:
            analysisSpectrum = view.analysisSpectrum
            spectrum = analysisSpectrum.dataSource
            experiment = spectrum.experiment
            spinningRate = experiment.spinningRate
            if spinningRate:
              for label in ('x', 'y'):
                dataDim = view.findFirstAxisMapping(label=label).analysisDataDim.dataDim
                dataDimRef = getPrimaryDataDimRef(dataDim)
                expDimRef = dataDimRef.expDimRef
                if label == 'x':
                  posn = positions[0]
                  spinningRate /= expDimRef.sf  # assumes y expDimRef would give the same
                else:
                  posn = positions[1]
                minAliasedFreq = expDimRef.minAliasedFreq
                maxAliasedFreq = expDimRef.maxAliasedFreq
                if minAliasedFreq is None:
                  minAliasedFreq = Util.convertPosition(float(dataDim.numPoints), dataDimRef, toUnit=axisUnit.unit)
                if maxAliasedFreq is None:
                  maxAliasedFreq = Util.convertPosition(1.0, dataDimRef, toUnit=axisUnit.unit)
                maxRange = int(math.floor((maxAliasedFreq-posn)/spinningRate))
                minRange = int(math.ceil((minAliasedFreq-posn)/spinningRate))
                rng = range(minRange, maxRange+1)
                if label == 'x':
                  xrng = rng
                else:
                  yrng = rng
              posns = positions[:]
              for j in xrng:
                posns[0] = positions[0] + j*spinningRate
                for k in yrng:
                  posns[1] = positions[1] + k*spinningRate
                  if j and k:  # skip if both 0 because then ordinary mark
                    remove = False
                    createNonPeakMark(posns, axisTypes, remove=remove)

    # ordinary mark

    if peak:
      axisTypeDict = {}
      if xaxisType.measurementType == 'MQShift':
        dataSource = peak.peakList.dataSource
        analysisSpectrum = dataSource.analysisSpectrum
        view = windowPane.findFirstSpectrumWindowView(analysisSpectrum=analysisSpectrum)
        axisMapping = view.findFirstAxisMapping(label='x')
        dataDim = axisMapping.analysisDataDim.dataDim
        for peakDim in peak.peakDims:
          if peakDim.dataDimRef.dataDim == dataDim:
            axisTypeDict[peakDim] = xaxisType
      elif yaxisType.measurementType == 'MQShift':
        dataSource = peak.peakList.dataSource
        analysisSpectrum = dataSource.analysisSpectrum
        view = windowPane.findFirstSpectrumWindowView(analysisSpectrum=analysisSpectrum)
        axisMapping = view.findFirstAxisMapping(label='y')
        dataDim = axisMapping.analysisDataDim.dataDim
        for peakDim in peak.peakDims:
          if peakDim.dataDimRef.dataDim == dataDim:
            axisTypeDict[peakDim] = yaxisType
        
      createPeakMark(peak, remove=remove, axisTypeDict=axisTypeDict)
    else:
      positions = self.findPosition(a, b, strip, returnDict=False)
      axisTypes = [axisPanel.axisType for axisPanel in windowPane.sortedAxisPanels()]
      createNonPeakMark(positions, axisTypes, remove=remove)

  def moveSelectedPeak(self, canvas, x, y):

    peaks = self.topPopup.currentPeaks
    if (len(peaks) > 1):
      showError('Multiple peaks', 'Multiple peaks selected, can only move one at a time.', parent=self)
      return
    elif (not len(peaks)):
      return

    peak = peaks[0]
    self.movePeak(peak, canvas, x, y)

  def moveSelectedPeakAnnotation(self, canvas, x, y):

    peaks = self.topPopup.currentPeaks
    if (len(peaks) > 1):
      showError('Multiple peaks', 'Multiple peaks selected, can only move one annotation at a time.', parent=self)
      return
    elif (not len(peaks)):
      return

    peak = peaks[0]
    self.movePeakAnnotation(peak, canvas, x, y)

  def startPositionDelta(self, canvas, x, y):
    
    getPanel = self.windowPane.findFirstAxisPanel
    (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
    xAxisRegion = getPanel(label='x').sortedAxisRegions()[col]
    yAxisRegion = getPanel(label='y').sortedAxisRegions()[row]
    (a, b, x, y) = self.calcWorldCoord(canvas, x, y)

    if self.windowPane.spectrumWindow.stripAxis == 'x':
      strip = col
    else:
      strip = row

    positions = self.findPosition(a, b, strip, returnDict=False)[:2]
    self.deltaPositions = positions
    self.drawAllAfter()

  def endPositionDelta(self):

    self.deltaPositions = None
    self.drawAllAfter()
    
  def movePeak(self, peak, canvas, x, y):

    spectrum = peak.peakList.dataSource
    view = getSpectrumWindowView(self.windowPane, spectrum)

    if (not view):
      showError('Peak not in window', 'Selected peak not in this window.', parent=self)
      return

    (a, b, x, y) = self.calcWorldCoord(canvas, x, y)
    (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
    if self.windowPane.spectrumWindow.stripAxis == 'x':
      position = self.findPosition(a, b, col)
    else:
      position = self.findPosition(a, b, row)

    (spectrum_position, spectrum_tile) = self.determinePosition(view, position)

    setPeakPosition(peak, spectrum_position, spectrum_tile)

  def movePeakAnnotation(self, peak, canvas, x, y):

    spectrum = peak.peakList.dataSource
    view = getSpectrumWindowView(self.windowPane, spectrum)

    if (not view):
      showError('Peak not in window', 'Selected peak not in this window.', parent=self)
      return

    (a, b, x, y) = self.calcWorldCoord(canvas, x, y)
    (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
    if self.windowPane.spectrumWindow.stripAxis == 'x':
      position = self.findPosition(a, b, col)
    else:
      position = self.findPosition(a, b, row)

    (position, tile) = self.determinePosition(view, position)

    if self.hasValueAxis:
      labels = ('x',)
      peakIntensity = peak.findFirstPeakIntensity(intensityType='height')
      if peakIntensity:
        offset = (b - peakIntensity.value) / spectrum.scale
        Util.setPeakTextOffset(peak, offset)
    else:
      labels = ('x','y')
    for label in labels:
      dataDim = view.findFirstAxisMapping(label=label).analysisDataDim.dataDim
      dim = dataDim.dim
      peakDim = peak.findFirstPeakDim(dim=dim)
      dataDimRef = peakDim.dataDimRef
      if dataDimRef:
        npts = dataDimRef.dataDim.numPointsOrig
        dim = dim - 1
        offset = position[dim] - peakDim.position
        offset = offset % npts
        if offset > npts/2: # choose closest offset
          offset = offset - npts

        Util.setPeakDimTextOffset(peakDim, offset)

  def snapSelectedPeaks(self, *event):

    peaks = self.topPopup.currentPeaks
    snapPeaks(peaks)

  def assignPeakAtLocation(self, canvas, x, y):

    (x, y, propX, propY) = self.calcWorldCoord(canvas, x, y)
    (row, col)  = self.scrolled_window.getCanvasRowCol(canvas)
    getPanel = self.windowPane.findFirstAxisPanel
    xAxisRegion = getPanel(label='x').sortedAxisRegions()[col]
    yAxisRegion = getPanel(label='y').sortedAxisRegions()[row]

    if self.windowPane.spectrumWindow.stripAxis == 'x':
      region = self.findPositionRegion(x, y, col)
    else:
      region = self.findPositionRegion(x, y, row)

    peak = self.findNearestPeak(region, xAxisRegion, yAxisRegion)
    if peak:
      self.menuPeak = peak
      self.assignPeak()

  # Hack for Windows because ScrolledWindow focus pops window to front
  # so need to put bind on WindowPopup rather than on ScrolledWindow.
  # This hack mangles it so that can pretend event was on ScrolledWindow.
  # It seems that on Windows the event is relative to the top-left of self
  # (so just under the title bar).
  # (On OSX it is relative to the top-left of the scrolled_window canvases.
  # 13 Mar 2012: decide that the popping up of the focus is not that bad
  # And also, without the focus on enter() in ScrolledWindow then just
  # clicking in the z Entry in the RegionSelector means that can no longer
  # get any events to the normal window, they all go to the RegionSelector
  def modifyKeyEvent(self, event):

    #print 'modifyKeyEvent1', self.winfo_rootx(), self.winfo_rooty()
    #print 'modifyKeyEvent2', event.x, event.y

    x = self.winfo_rootx() + event.x
    y = self.winfo_rooty() + event.y

    #print 'modifyKeyEvent3', x, y

    canvases = self.scrolled_window.canvases

    ncols = len(canvases[0])
    if ncols == 1:
      col = 0
      canvas = canvases[0][0]
      x = x - canvas.winfo_rootx()
    else:
      for col in range(ncols-1):
        canvas = canvases[0][col+1]
        #print 'modifyKeyEvent: in x loop:', col, x, canvas.winfo_rootx()
        if x < canvas.winfo_rootx():
          canvas = canvases[0][col]
          x = x - canvas.winfo_rootx()
          break
      else:
        canvas = canvases[0][ncols-1]
        x = x - canvas.winfo_rootx()

    nrows = len(canvases)
    if nrows == 1:
      row = 0
      canvas = canvases[0][0]
      y = y - canvas.winfo_rooty()
    else:
      for row in range(nrows-1, 0, -1):
        canvas = canvases[row-1][0]
        #print 'modifyKeyEvent: in y loop:', col, y, canvas.winfo_rooty()
        if y < canvas.winfo_rooty():
          canvas = canvases[row][0]
          y = y - canvas.winfo_rooty()
          break
      else:
        canvas = canvases[0][0]
        y = y - canvas.winfo_rooty()

    #print 'modifyKeyEvent4', canvas.winfo_rootx(), canvas.winfo_rooty()
    #print 'modifyKeyEvent5', x, y
    #print 'modifyKeyEvent6', row, col

    event.widget = canvases[row][col]
    event.x = x
    # this adjusts the root origin for Windows - will break when float windows introduced
    # TBD: will also break if exact buttonList implementation changes
    dy = self.windowPopup.buttonList.parent.parent.master.winfo_height()
    if self.windowPopup.activeButtonFrame:
      dy += self.windowPopup.activeButtonFrame.winfo_height()
    event.y = y - dy
 
    
  def keypress(self, event):

    if self.keypressTime is not None:
      time  = event.time
      delta = time - self.keypressTime
      if delta < 100: # Fix for really fast key repeat speeds slowing things down
        return
       
    self.keypressTime = event.time

    if self.waitDraw:
      return

    # 13 Mar 2012: see comment above modifyKeyEvent() definition
    #if isWindowsOS():
    #  self.modifyKeyEvent(event)

    keysym = event.keysym
    if keysym.startswith('Control_') or keysym.startswith('Shift_'):
      return

    if keysym[0] == 'F' and len(keysym) > 1:  # F keys
      try:
        int(keysym[1:])  # check that F followed by number
        state = event.state & 255
        if state & shift_key_state:
          keysym = Util.shift_key_modifier + keysym
        if state & alt_key_state:
          keysym = Util.alt_key_modifier + keysym
        if state & ctrl_key_state:
          keysym = Util.ctrl_key_modifier + keysym
      except:
        pass
      
    areShortcutsTwoChars = self.windowPopup.analysisProfile.twoCharShortcuts
    if areShortcutsTwoChars:
      if keysym == 'Escape':
        self.keysym = None
        return
      if self.keysym:
        keysym = self.keysym + keysym
        self.keysym = None

    #print 'keypress1', keysym
    spectrum = self.topPopup.toggleSpectrum(self.windowPane.spectrumWindow, shortcut=keysym)
    #print 'keypress1A', spectrum
    if not spectrum: # this means was not spectrum shortcut, so try other macros
      macro = self.windowPopup.analysisProfile.findFirstMacro(shortcut=keysym)
      #print 'keypress2', macro
      if macro:
        self.setCurrentObjects(event)
        Util.runMacro(macro, self.topPopup.argumentServer)
      elif areShortcutsTwoChars:
        if len(keysym) == 1:
          self.keysym = keysym
        else:
          self.keysym = None

  def sliceKeypress(self, event):

    #print 'sliceKeypress', event.keysym
    if event.keysym == 'Home':
      self.changeSliceRange(0.5, event.widget)
    elif event.keysym == 'End':
      self.changeSliceRange(2.0, event.widget)

  def popupMenu(self, event):

    self.scrolled_window.popupMenu(event)

  def changeSliceRange(self, scale, canvas = None):

    if self.hasValueAxis:
      if canvas and hasattr(canvas, 'ymax_extent'):
        ymax_extent = canvas.ymax_extent
      else:
        ymax_extent = None
      yAxisPanel = self.windowPane.findFirstAxisPanel(label='y')
      ymin = ymax = None
      for axisRegion in yAxisPanel.sortedAxisRegions():
        (r0, r1) = axisRegion.region
        if ymax_extent and r1 > r0:
          maxScale = ymax_extent / (r1-r0)
          scale = min(scale, maxScale)
        r0 = scale * r0
        r1 = scale * r1
        axisRegion.region = (r0, r1)
        if ymin is None:
          ymin = r0
          ymax = r1
        else:
          ymin = min(ymin, r0)
          ymax = max(ymax, r1)
      self.windowPane.sliceRange = (ymin, ymax)
    else:
      (r0, r1) = self.windowPane.sliceRange

      """ changed 20 Nov 2006 to scale relative to 0
      a = 0.5 * (r0 + r1)
      b = 0.5 * scale * (r1 - r0)

      self.windowPane.sliceRange = (a - b, a + b)
      """
      self.windowPane.sliceRange = (scale * r0, scale * r1)

  # overrides WindowDraw version
  def drawViewTile(self, handler, canvas, view, contourLevels, contourStyle,
                   worldPointRanges, spectrumPointRanges, row, components=None):

    #if (self.windowPane.name == 'w'):
    #  print 'entering drawViewTile', self.windowPane.name

    if canvas.doubleBuffer:

      #print 'drawViewTile double', self.windowPane.name

      if not isWindowsOS():
        t = time.time()
        if ((t - canvas.initTime) > doubleBufferTimeout):
          raise WindowTimeoutException()

      # exceptions caught higher up
      self.drawViewTileReal(handler, view, contourLevels, contourStyle,
                            worldPointRanges, spectrumPointRanges, row, components)

    else:

      #print 'drawViewTile single', self.windowPane.name

      # exceptions caught here because called with after_idle
      try:
        #print 'drawViewTile0', self.windowPane.name
        self.drawViewTileReal(handler, view, contourLevels, contourStyle,
                              worldPointRanges, spectrumPointRanges, row)
        #print 'drawViewTile1', self.windowPane.name
      except Implementation.ApiError, e:
        print 'Drawing canvas tile error:', e.error_msg
      except self.handlerExc, e:
        print 'Drawing canvas tile handler error:', e
      except ContourFile.error, e:
        print 'Drawing canvas tile ContourFile error:', e
      except SliceFile.error, e:
        print 'Drawing canvas tile SliceFile error:', e
      except PeakList.error, e:
        print 'Drawing canvas tile PeakList error:', e
      except WinPeakList.error, e:
        print 'Drawing canvas tile WinPeakList error:', e
      except:
        #print 'Unknown canvas tile error'
        pass

  # overrides WindowDraw version
  def drawCanvas(self, canvas, row = None, col = None):

    if self.scrolled_window.isDeleting:
      return

    if self.windowPopup.state() != 'normal':
      return

    # if self.waitDraw:
    #   return

    # below does not work because get splatted by Tcl/Tk if do nothing
    # if (hasattr(canvas, 'beingDrawn') and canvas.beingDrawn):
    #   return

    if row is None or col is None:
      try:
        (row, col) = self.scrolled_window.getCanvasRowCol(canvas)
      except:
        return

    self.drawCanvasCount = self.drawCanvasCount + 1

    #canvas.beingDrawn = True
    if not hasattr(canvas, 'drawCount'):
      canvas.drawCount = 0

    #print 'drawCanvas1', self.windowPane.name, canvas.drawCount
    canvas.drawCount = canvas.drawCount + 1
    self.after_idle(lambda: self.drawCanvasReal(canvas, row, col))
    #print 'drawCanvas2', self.windowPane.name

  ##def getCanvasState(self, canvas, row, col):
  ##
  ##  state = {}
  ##  state['width'] = canvas.canvas_width
  ##  state['height'] = canvas.canvas_height
  ##
  ##  return state

  # overrides WindowDraw version
  def drawCanvasReal(self, canvas, row, col):

    # print 'drawCanvasReal0', self.windowPane.name

    canvas.drawCount = canvas.drawCount - 1
    #print 'drawCanvasReal1', self.windowPane.name, canvas.drawCount
    if canvas.drawCount:
      return # only draw when there is only one request remaining
    #print 'drawCanvasReal2', self.windowPane.name

    if not hasattr(canvas, 'canvas_width'):
      #canvas.beingDrawn = False
      return

    if not hasattr(canvas, 'xview_region'):
      #canvas.beingDrawn = False
      return

    if self.waitDraw: # In the middle of processing another draw
      return

    #print 'drawCanvasReal3', self.windowPane.name, canvas.canvas_width, canvas.canvas_height

    ## optimisation: do not draw again if state is the same
    #state = self.getCanvasState(canvas, row, col)
    #if (hasattr(canvas, 'state')):
    #  if state == canvas.state:
    #    return
    #canvas.state = state

    # safety (check that canvas still what you think it is)
    try:
      (rr, cc) = self.scrolled_window.getCanvasRowCol(canvas)
    except:
      return

    if row != rr or col != cc:
      return

    #print 'drawCanvasReal4', self.windowPane.name, canvas.canvas_width, canvas.canvas_height

    self.setupWidgetHandler(canvas, isCanvas=True)

    handler = canvas.handler
    if not handler:
      return
    
    try:

      canvas.doubleBuffer = True
      handler.setIsDoubleBuffer(1)

      canvas.initTime = time.time()

      handler.makeCurrent()
      handler.startBack()
      handler.expose(0, 0, canvas.canvas_width, canvas.canvas_height)
      #handler.clearXor()

      self.doCanvas(handler, canvas, row, col)
      #####handler.flush()

      # doing an after_idle probably messes things up since can
      # have multiple draws before a single swapBuffer is called
      ##canvas.after_idle(canvas.handler.swapBuffers)
      handler.endBack()
      handler.swapBuffers()
      #canvas.beingDrawn = False

    except WindowTimeoutException:

      #print 'drawCanvasReal: about to single buffer doCanvas', self.windowPane.name
      canvas.doubleBuffer = False
      handler.setIsDoubleBuffer(0)

      handler.makeCurrent()
      handler.clearXor()
      handler.startFront()
      handler.expose(0, 0, canvas.canvas_width, canvas.canvas_height)

      self.doCanvas(handler, canvas, row, col)
      handler.flush()

    except Implementation.ApiError, e:
      print 'Drawing real canvas error:', e.error_msg
      #canvas.beingDrawn = False

    except self.handlerExc, e:
      print 'Drawing real canvas handler error:', e

    except ContourFile.error, e:
      print 'Drawing real canvas ContourFile error:', e

    except SliceFile.error, e:
      print 'Drawing real canvas SliceFile error:', e

    except ContourLevels.error, e:
      print 'Drawing real canvas ContourLevels error:', e

    except ContourStyle.error, e:
      print 'Drawing real canvas ContourStyle error:', e

    except PeakList.error, e:
      print 'Drawing real canvas PeakList error:', e

    except WinPeakList.error, e:
      print 'Drawing real canvas WinPeakList error:', e

    except:
      #canvas.beingDrawn = False
      #print 'Unknown drawing real canvas error'
      raise
      pass

    #print 'drawCanvas2', self.windowPane.name

  # overrides WindowDraw version
  def drawRow(self, row):

    if self.windowPane.isDeleted:
      return

    #print 'WindowPopup: drawRow1'

    ncols = self.getNCols()
    canvases = self.scrolled_window.canvases
    for i in range(ncols):
      canvas = canvases[row][i]
      self.drawCanvas(canvas, row, i)

    #print 'WindowPopup: drawRow2'

  # overrides WindowDraw version
  def drawCol(self, col):

    if self.windowPane.isDeleted:
      return

    nrows = self.getNRows()
    canvases = self.scrolled_window.canvases
    for j in range(nrows):
      canvas = canvases[j][col]
      self.drawCanvas(canvas, j, col)

  def drawAllAfter(self, *extra):

    if self.waitDraw:
      return

    self.waitDraw = True
    self.after_idle(self.drawAll)

  # overrides WindowDraw version
  def drawAll(self, *extra):

    if self.windowPane.isDeleted:
      self.waitDraw = False
      return

    # exceptions caught here because called with after_idle
    try:
      w = self.scrolled_window
      for j in range(w.nrows):
        self.drawRow(j)

      self.drawAllSlices()
    except Implementation.ApiError, e:
      print 'Drawing all error:', e.error_msg
    except self.handlerExc, e:
      print 'Drawing all handler error:', e
    except ContourFile.error, e:
      print 'Drawing all ContourFile error:', e
    except SliceFile.error, e:
      print 'Drawing all SliceFile error:', e
    except PeakList.error, e:
      print 'Drawing all PeakList error:', e
    except WinPeakList.error, e:
      print 'Drawing all WinPeakList error:', e
    except ContourLevels.error, e:
      print 'Drawing all ContourLevels error:', e
    except ContourStyle.error, e:
      print 'Drawing all ContourStyle error:', e
    except:
      #print 'Unknown drawing all error'
      pass

    #print 'WindowPopup: drawAll2'
    self.waitPeak = False
    self.waitDraw = False

    if self.fontChanged:
      self.fontChanged = False
      self.drawAllAfter()

  def drawViewSlice(self, slice, view, sliceRange = None):

    from ccpnmr.c.SliceFile import SliceFile
    from ccpnmr.analysis.core.BlockUtil import getBlockFile

    analysisProject = view.topObject

    if isWindowsOS():
      if hasattr(view, '_createSliceFileObjectX'):
        xdim, _, mem_cache = view._createSliceFileObjectX
        dataSource = view.analysisSpectrum.dataSource
        writeable = hasattr(dataSource, 'writeable') and dataSource.writeable
        self._tempblock_fileX = getBlockFile(dataSource,
                                  mem_cache=mem_cache,
                                  writeable=writeable)
        if self._tempblock_fileX:
          view.sliceFile['x'] = SliceFile(1, xdim, self._tempblock_fileX, mem_cache)
        del view._createSliceFileObjectX

      if hasattr(view, '_createSliceFileObjectY'):
        ydim, _, mem_cache = view._createSliceFileObjectY
        dataSource = view.analysisSpectrum.dataSource
        writeable = hasattr(dataSource, 'writeable') and dataSource.writeable
        self._tempblock_fileY = getBlockFile(dataSource,
                                  mem_cache=mem_cache,
                                  writeable=writeable)
        if self._tempblock_fileY:
          view.sliceFile['y'] = SliceFile(0, ydim, self._tempblock_fileY, mem_cache)
        del view._createSliceFileObjectY

    #print 'drawViewSliceA', view.analysisSpectrum.dataSource.name
    label = slice.label
    sliceFile = view.sliceFile.get(label)
    if (not sliceFile):
      return

    #print 'drawViewSliceB', view.analysisSpectrum.dataSource.name
    handler = slice.handler
    if not handler:
      return

    # TBD: do again when view.sliceColor in data model
    analysisSpectrum = view.analysisSpectrum
    color = analysisSpectrum.sliceColor
 
    if isinstance(slice, FakeSlice):
      bg = self.windowPopup.analysisProfile.bgColor
      if isMacOS():
        # only Tk
        handler.setColor(hexToRgb(color))

      else:
        handler.setColor(hexToRgb(hexXor(color, bg)))

    else:
      handler.setColor(hexToRgb(color))

    #print 'drawViewSliceC', view.analysisSpectrum.dataSource.name
    position = analysisSpectrum.dataSource.numDim * [0]

    for axisMapping in view.axisMappings:
      axisPanel = self.windowPane.findFirstAxisPanel(label=axisMapping.label)
      #print 'drawViewSliceC1', self.windowPane.name, axisPanel, label, hasattr(axisPanel, 'pointLocation') and axisPanel.pointLocation
      if (axisPanel.label == label): # position irrelevant for this dim
        continue
      if (not hasattr(axisPanel, 'pointLocation')):
        return
      p = axisPanel.pointLocation
      if axisPanel.label not in ('x', 'y'):
        m = None
        stripAxis = self.windowPane.spectrumWindow.stripAxis
        orient = slice.orient[0].lower()
        if (stripAxis == 'x') and (orient == 'h'):
          if isinstance(slice, FakeSlice):
            m = slice.col
          else:
            m = self.scrolled_window.getSliceCol(slice)
        elif (stripAxis == 'y') and (orient == 'v'):
          if isinstance(slice, FakeSlice):
            m = slice.row
          else:
            m = self.scrolled_window.getSliceRow(slice)
        if m is not None:
          axisRegion = axisPanel.sortedAxisRegions()[m]
          p = self.calcAxisMidpoint(axisRegion)
      if (p is None):
        return
      dataDim = axisMapping.analysisDataDim.dataDim
      if dataDim.className == 'FreqDataDim':
        # -1 below to get points starting from 0
        p = Util.convertPosition(p, getPrimaryDataDimRef(dataDim),
                            fromUnit=axisPanel.axisUnit.unit) - 1
        n = dataDim.numPointsOrig
        p = position[dataDim.dim-1] = p % n
        if p >= dataDim.numPoints:
          return
      else:
        axisRegion = axisPanel.findFirstAxisRegion()
        p = position[dataDim.dim-1] = int(axisRegion.region[0] - 1 + 0.0001)
        if p >= dataDim.numPoints:
          return

    #print 'drawViewSlice1', label, position, view.analysisSpectrum.dataSource.name

    slicePanel = self.windowPane.findFirstSlicePanel(label=label)
    y0 = 0.0
    y1 = 1.0
    if not sliceRange:
      sliceRange = self.windowPane.sliceRange
    (b0, b1) = sliceRange
    scale = view.analysisSpectrum.dataSource.scale
    (b0, b1) = (b0/scale, b1/scale)

    axisMapping = view.findFirstAxisMapping(label=label)
    dataDim = axisMapping.analysisDataDim.dataDim
    axisPanel = self.windowPane.findFirstAxisPanel(label=label)
    #print 'drawViewSlice1A', slice.view_region, dataDim.dim
    (t0, t1) = Util.convertRegion(slice.view_region, axisPanel.axisUnit, dataDim)
    if (t0 > t1):
      (t0, t1) = (t1, t0)

    dataDimRef = getPrimaryDataDimRef(dataDim)
    expDimRef = dataDimRef.expDimRef
    minAliasedFreq = expDimRef.minAliasedFreq
    maxAliasedFreq = expDimRef.maxAliasedFreq
    unit = axisPanel.axisType.findFirstAxisUnit(unit='ppm').unit
    if (maxAliasedFreq is None):
      minFreqPts = 0.0
    else:
      minFreqPts = Util.convertPosition(maxAliasedFreq, dataDimRef, unit) - 1.0
    if (minAliasedFreq is None):
      maxFreqPts = float(dataDim.numPoints) - 1.0
    else:
      maxFreqPts = Util.convertPosition(minAliasedFreq, dataDimRef, unit) - 1.0

    n = dataDim.numPointsOrig
    tile0 = int(math.floor(max(t0, minFreqPts) / n))
    tile1 = int(math.floor(min(t1, maxFreqPts) / n))
    #print 'drawViewSlice2', label, t0, t1, n, tile0, tile1

    for tile in range(tile0, tile1+1):
      o = tile * n
      if (tile == tile0):
        s0 = t0
      else:
        s0 = tile * n
      if (tile == tile1):
        s1 = t1
      else:
        s1 = (tile+1) * n
      #print 'drawViewSlice3', label, tile, s0, s1, t0, t1
      s0 = max(s0, minFreqPts)
      s1 = min(s1, maxFreqPts)
      #print 'drawViewSlice4', label, tile, s0, s1, t0, t1
      if (s1 <= s0):
        continue

      x0 = (s0-t0) / (t1-t0)
      x1 = (s1-t0) / (t1-t0)
      # TBD: here set s0 and s1 to lie in freq range of spectrum
      a0 = float(s0 - o)
      a1 = float(s1 - o)
      first = int(math.floor(a0))
      last = min(dataDim.numPoints, int(math.ceil(a1+1)))
      if first >= last:
        continue

      #print 'drawViewSlice5', label, x0, y0, x1, y1, a0, b0, a1, b1
      if (slice.orient == Tkinter.HORIZONTAL):
        handler.mapRanges(x0, y0, x1, y1, a0, b0, a1, b1)
      else:
        # note b1, b0 swap, due to wanting positive values to point left, not right
        handler.mapRanges(y0, x0, y1, x1, b1, a0, b0, a1)
      #print 'drawViewSlice6', self.windowPane.name, label, first, last, position
      sliceFile.draw(handler, first, last, position)
      #print 'drawViewSlice7', label

  def doSlice(self, slice):

    for view in self.getSpectrumViews():
      #print 'doSlice', view.analysisSpectrum.dataSource.name, view.isSliceVisible
      if view.isSliceVisible:
        if self.hasValueAxis:
          self.drawViewSlicePeaks(slice, view, slice.winfo_height())
        else:
          self.drawViewSlice(slice, view)
        if slice.handler:
          slice.handler.flush()

    if not self.hasValueAxis:
      self.drawSliceAxis(slice)

  def drawSlice(self, slice):

    if self.windowPopup.state() != 'normal':
      return

    if (not hasattr(slice, 'drawCount')):
      slice.drawCount = 0

    slice.drawCount = slice.drawCount + 1
    self.after_idle(lambda: self.drawSliceReal(slice))

  def drawSliceReal(self, slice):

    slice.drawCount = slice.drawCount - 1
    if slice.drawCount:
      return # only draw when there is only one request remaining

    if not hasattr(slice, 'slice_width'):
      return

    if not hasattr(slice, 'view_region'):
      return

    #print 'drawSliceReal', slice.slice_width, slice.slice_height
    self.setupWidgetHandler(slice, isCanvas=False)

    try:

      handler = slice.handler
      if not handler:
        return

      handler.makeCurrent()
      handler.startBack()
      handler.expose(0, 0, slice.slice_width, slice.slice_height)
      handler.clearXor()

      self.doSlice(slice)

      handler.endBack()
      handler.swapBuffers()

    except Implementation.ApiError, e:
      print 'Drawing slice error:', e.error_msg

    except self.handlerExc, e:
      print 'Drawing slice handler error:', e

    except SliceFile.error, e:
      print 'Drawing slice SliceFile error:', e

    except:
      #print 'Unknown drawing slice error'
      pass

    try:
      if (slice.label == 'x'):
        #self.drawSliceCrosshair(slice, x=self.windowPane.findFirstAxisPanel(label='x').pointLocation)
        self.drawSliceCrosshairs(slice, xs=self.crosshairLocations[0])
      else:
        #self.drawSliceCrosshair(slice, y=self.windowPane.findFirstAxisPanel(label='y').pointLocation)
        self.drawSliceCrosshairs(slice, ys=self.crosshairLocations[1])
    except:
      pass

  def drawSliceCol(self, col):

    self.drawSlice(self.scrolled_window.xslices[col])

  def drawSliceRow(self, row):

    self.drawSlice(self.scrolled_window.yslices[row])

  def drawAllSlices(self):

    slicePanels = self.windowPane.sortedSlicePanels()

    if slicePanels[0].isVisible:
      for i in range(self.scrolled_window.ncols):
        self.drawSliceCol(i)

    if slicePanels[1].isVisible:
      for j in range(self.scrolled_window.nrows):
        self.drawSliceRow(j)

  def drawSliceAxis(self, slice):

    handler = slice.handler
    if not handler:
      return

    color = hexToRgb(inverseGrey(self.windowPopup.analysisProfile.bgColor))
    handler.setColor(color)

    slicePanel = self.windowPane.findFirstSlicePanel(label=slice.label)
    y0 = 0.0
    y1 = 1.0
    (b0, b1) = self.windowPane.sliceRange
    x0 = 0.0
    x1 = 1.0
    a0 = 0.0
    a1 = 1.0
    if (slice.orient == Tkinter.HORIZONTAL):
      handler.mapRanges(x0, y0, x1, y1, a0, b0, a1, b1)
      handler.drawLine(0.0, 0.0, 1.0, 0.0)
    else:
      handler.mapRanges(y1, x0, y0, x1, b0, a0, b1, a1)
      handler.drawLine(0.0, 0.0, 0.0, 1.0)

  def drawSliceCrosshairs(self, slice, xs = None, ys = None):

    if (self.scrolled_window.isDeleting):
      return

    # this can be called before widget has been set up properly
    if (not hasattr(slice, 'handler')):
      return

    handler = slice.handler
    if not handler:
      return

    if (not self.topPopup.isCrosshairVisible()):
      return

    if xs is None:
      xs = []

    if ys is None:
      ys = []

    handler.makeCurrent()
    handler.startXor()
    handler.setBlack()
    # above in xor mode means that below will appear in inverse colors
    handler.mapRanges(0, 0, 1, 1, 0, 0, 1, 1)

    for x in xs:
      (x0, x1) = slice.view_region
      x = float(x - x0) / (x1 - x0)
      handler.drawLine(x, 0, x, 1)

    for y in ys:
      (y0, y1) = slice.view_region
      y = float(y - y0) / (y1 - y0)
      handler.drawLine(0, y, 1, y)

    handler.finishXor()

  def initAxisPanel(self, axisPanel):

    if axisPanel.spectrumWindowPane != self.windowPane:
      return

    self.createRegionSelector(axisPanel)

  def deleteAxisPanel(self, axisPanel):

    if axisPanel.spectrumWindowPane != self.windowPane:
      return

    self.deleteRegionSelector(axisPanel)

  # TBD: could make initPeak and drawPeak more intelligent
  def initPeak(self, peak):

    if self.waitPeak:
      return

    self.waitPeak = True
    self.after_idle(self.drawAll)

  def deletePeak(self, peak):

    if self.waitPeak:
      return

    self.waitPeak = True
    self.after_idle(self.drawAll)

  def splitWindowHorizontal(self, event):

    pass

  def splitWindowVertical(self, event):

    canvas = event.widget
    col = self.scrolled_window.getCanvasRowCol(canvas)[1]
    self.scrolled_window.splitWindowVertical(event)
    new_canvas = self.scrolled_window.canvases[0][col+1] # row 0 arbitrary
    (r0, r1) = new_canvas.xview_region
    if (r0 > r1):
      (r0, r1) = (r1, r0)
    Util.addAxisPanelRegion(self.windowPane.findFirstAxisPanel(label='x'), region=(r0, r1),
                 size=new_canvas.winfo_width())
    self.setHandlersCol(col+1)
    self.setHandlerBackground()

  def addRow(self, row=-1, centerPositionDict = None):

    #print 'addRow', row
    if centerPositionDict is None:
      centerPositionDict = {}

    self.scrolled_window.addRow()

    canvas = self.scrolled_window.canvases[row][0] # col 0 arbitrary
    (r0, r1) = canvas.yview_region
    if (r0 > r1):
      (r0, r1) = (r1, r0)

    yAxisPanel = self.windowPane.findFirstAxisPanel(label='y')

    Util.addAxisPanelRegion(yAxisPanel, region=(r0, r1),
                            size=canvas.winfo_height())

    self.setHandlersRow()
    self.setHandlerBackground()
    self.shuffleAxisRegions(sourceNum=row, doCols=False)

    if (row != -1):
      if (row == 0):
        row = 1
      axisPanels = self.windowPane.sortedAxisPanels()
      for i in range(len(axisPanels)):
        if (i == 0):
          continue
        axisRegions = axisPanels[i].sortedAxisRegions()
        axisRegions[row].region = axisRegions[row-1].region

    if len(self.windowPane.findFirstAxisPanel(label='y').axisRegions) > 1:
      self.windowPopup.stripsButtons.buttons[1].enable()

    if row == -1:
      self.activateStrip(len(yAxisPanel.axisRegions)-1)
    else:
      self.activateStrip(row)

    for label in centerPositionDict.keys():
      axisPanel = self.windowPane.findFirstAxisPanel(label=label)
      axisRegion = axisPanel.sortedAxisRegions()[row]
      centerAxisRegion(axisRegion, centerPositionDict[label])


  def deleteRow(self, row=-1):

    #print 'deleteRow'
    self.update_idletasks()
    N = len(self.windowPane.findFirstAxisPanel(label='y').axisRegions)
    if N > 1:
      active_row = self.getActiveStrip()
      if (row == -1) and (active_row != None) and (active_row != 'all') and self.windowPane.spectrumWindow.stripAxis == 'y':
        row = active_row

      self.scrolled_window.deleteRow(row)
      WindowDraw.deleteRow(self, row)
      if (active_row == row):
        self.activateStrip(min(row, N-2)) # new active row

      if N < 3:
        self.windowPopup.stripsButtons.buttons[1].disable()

  def setRowsTo(self, nrows = 1):

    self.update_idletasks() # see if this helps crashing stop

    old_nrows = self.scrolled_window.nrows
    if old_nrows == nrows:
      return

    turnedOff = self.turnDrawRequestsOff()
    self.waitResize = True

    try:

      active_row = self.getActiveStrip()
      newActive = (nrows > old_nrows) or (active_row >= nrows)

      self.scrolled_window.setRowsTo(nrows)

      if old_nrows < nrows:  # add rows
        canvas = self.scrolled_window.canvases[-1][0] # col 0 arbitrary
        size = canvas.winfo_height()
        (r0, r1) = canvas.yview_region
        if r0 > r1:
          (r0, r1) = (r1, r0)
        axisPanel = self.windowPane.findFirstAxisPanel(label='y')
        for i in range(old_nrows, nrows):
          Util.addAxisPanelRegion(axisPanel, region=(r0, r1), size=size)
          self.setHandlersRow(i)
        self.setHandlerBackground()
      else:  # delete rows
        for i in range(nrows, old_nrows):
          WindowDraw.deleteRow(self)

      if newActive:
        self.activateStrip(nrows-1) # new active col

      if nrows == 1:
        self.windowPopup.stripsButtons.buttons[1].disable()
      elif old_nrows == 1:
        self.windowPopup.stripsButtons.buttons[1].enable()

    finally:

      self.waitResize = False
      if turnedOff: self.turnDrawRequestsOn()

  def shuffleAxisRegions(self, sourceNum=-1, targetNum=-1, doCols=True):

    if sourceNum == -1:
      return

    axisPanels = self.windowPane.sortedAxisPanels()
    if doCols:
      n = len(axisPanels[0].axisRegions)
    else:
      n = len(axisPanels[1].axisRegions)

    if targetNum == -1:
      targetNum = n-1
    elif targetNum >= n:
      targetNum = n-1
      

    if sourceNum == targetNum:
      return

    if sourceNum < targetNum:
      d = 1
    else:
      d = -1

    turnedOff = self.turnDrawRequestsOff()

    try:

      for i in range(len(axisPanels)):
        if (i==1) and doCols: # do not need to do y
          continue
        elif (i==0) and not doCols:# do not need to do x
          continue

        axisRegions = axisPanels[i].sortedAxisRegions()
        region = axisRegions[sourceNum].region
        size = axisRegions[sourceNum].size
        isActive = axisRegions[sourceNum].isActive
        n = len(axisRegions)
        for j in range(sourceNum, targetNum, d):
          axisRegions[j].region = axisRegions[j+d].region
          axisRegions[j].size = axisRegions[j+d].size
          axisRegions[j].isActive = axisRegions[j+d].isActive
          
        axisRegions[targetNum].region = region
        axisRegions[targetNum].size = size
        axisRegions[targetNum].isActive = isActive
 
    finally:

      if turnedOff: self.turnDrawRequestsOn()

  def addStrip(self, num=-1, centerPositionDict=None):

    if self.windowPane.spectrumWindow.stripAxis == 'x':
      self.addCol(col=num, centerPositionDict=centerPositionDict)
    else:
      self.addRow(row=num, centerPositionDict=centerPositionDict)

  def deleteStrip(self, num=-1):

    if self.windowPane.spectrumWindow.stripAxis == 'x':
      self.deleteCol(col=num)
    else:
      self.deleteRow(row=num)

  def addCol(self, col=-1, centerPositionDict=None):

    #print 'addCol1', col
    if (centerPositionDict is None):
      centerPositionDict = {}

    #self.scrolled_window.oneWeights()

    # need to append everything as far as scrolled_window concerned
    # otherwise handlers and regions get confused
    # then need to do shuffle on axis regions to make up for this
    ###self.scrolled_window.addCol(col)
    self.scrolled_window.addCol()

    canvas = self.scrolled_window.canvases[0][col] # row 0 arbitrary
    (r0, r1) = canvas.xview_region
    if r0 > r1:
      (r0, r1) = (r1, r0)
    #print 'WindowPopup.addCol', (r0, r1)
    Util.addAxisPanelRegion(self.windowPane.findFirstAxisPanel(label='x'), region=(r0, r1),
                 size=canvas.winfo_width())

    # see above note
    ###self.setHandlersCol(col)
    self.setHandlersCol()
    self.setHandlerBackground()
    self.shuffleAxisRegions(sourceNum=col, doCols=True)

    if col != -1: # reset region to be that of strip where event set off
      if col == 0:
        col = 1
      axisPanels = self.windowPane.sortedAxisPanels()
      for i in range(len(axisPanels)):
        if (i == 1): # do not need to do y
          continue
        axisRegions = axisPanels[i].sortedAxisRegions()
        axisRegions[col].region = axisRegions[col-1].region

    if len(self.windowPane.findFirstAxisPanel(label='x').axisRegions) > 1:
      self.windowPopup.stripsButtons.buttons[1].enable()

    if col == -1:
      self.activateStrip(len(self.windowPane.findFirstAxisPanel(label='x').axisRegions)-1)
    else:
      self.activateStrip(col)

    for label in centerPositionDict.keys():
      axisPanel = self.windowPane.findFirstAxisPanel(label=label)
      axisRegion = axisPanel.sortedAxisRegions()[col]
      centerAxisRegion(axisRegion, centerPositionDict[label])

  def deleteCol(self, col = -1):

    #print 'deleteCol', col
    self.update_idletasks() # see if this helps crashing stop
    #self.scrolled_window.oneWeights()
    N = len(self.windowPane.findFirstAxisPanel(label='x').axisRegions)
    if N > 1:
      active_col = self.getActiveStrip()
      if (col == -1) and (active_col != None) and (active_col != 'all') and self.windowPane.spectrumWindow.stripAxis == 'x':
        col = active_col

      self.scrolled_window.deleteCol(col)
      WindowDraw.deleteCol(self, col)
      if (active_col == col):
        self.activateStrip(min(col, N-2)) # new active col

      if N < 3:
        self.windowPopup.stripsButtons.buttons[1].disable()

  def setColsTo(self, ncols = 1):

    self.update_idletasks() # see if this helps crashing stop

    old_ncols = self.scrolled_window.ncols
    if old_ncols == ncols:
      return

    turnedOff = self.turnDrawRequestsOff()
    self.waitResize = True

    try:

      active_col = self.getActiveStrip()
      newActive = (ncols > old_ncols) or (active_col >= ncols)

      self.scrolled_window.setColsTo(ncols)

      if old_ncols < ncols:
        canvas = self.scrolled_window.canvases[0][-1] # row 0 arbitrary
        size = canvas.winfo_width()
        (r0, r1) = canvas.xview_region
        if r0 > r1:
          (r0, r1) = (r1, r0)
        axisPanel = self.windowPane.findFirstAxisPanel(label='x')
        for i in range(old_ncols, ncols):
          Util.addAxisPanelRegion(axisPanel, region=(r0, r1), size=size)
          self.setHandlersCol(i)
        self.setHandlerBackground()
      else:
        for i in range(ncols, old_ncols):
          WindowDraw.deleteCol(self)

      if newActive:
        self.activateStrip(ncols-1) # new active col

      if ncols == 1:
        self.windowPopup.stripsButtons.buttons[1].disable()
      elif old_ncols == 1:
        self.windowPopup.stripsButtons.buttons[1].enable()

    finally:

      self.waitResize = False
      if turnedOff: self.turnDrawRequestsOn()

  def deleteSeparators(self):
  
    self.update_idletasks() # see if this helps crashing stop
    
    if self.windowPane.spectrumWindow.stripAxis == 'x':
      self.setRowsTo(1)
    else:
      self.setColsTo(1)
  
  def deleteStrips(self, num=None):

    self.update_idletasks() # see if this helps crashing stop
    if self.windowPane.spectrumWindow.stripAxis == 'x':
      axisPanel = self.windowPane.findFirstAxisPanel(label='x')
      n = len(axisPanel.axisRegions)
      
      if (num is None) or (num >= n):
        self.setColsTo(1)
      else: 
        self.setColsTo(n-num)

    else:
      axisPanel = self.windowPane.findFirstAxisPanel(label='y')
      n = len(axisPanel.axisRegions)
      
      if (num is None) or (num >= n):
        self.setRowsTo(1)
      else: 
        self.setRowsTo(n-num)


  def moveStripsLeft(self, num=None):

    self.moveStrips(delta=-1, num=num)

  def moveStripsRight(self, num=None):

    self.moveStrips(delta=1, num=num)

  def moveStrips(self, delta, num=None):

    N = len(self.windowPane.findFirstAxisPanel(label=self.windowPane.spectrumWindow.stripAxis).axisRegions)

    if num is None:
      num = self.getActiveStrip()
      if num is None:
        return

    if num == 'all':
      if delta > 0:
        while delta > 0:
          self.moveStripToWindow(N-1,self.windowPane,0)
          delta -= 1

      if delta < 0:
        while delta < 0:
          self.moveStripToWindow(0,self.windowPane,N-1)
          delta += 1


    else:
      num2 = num + delta
      if num2 >= N:
        num2 = 0

      self.moveStripToWindow(num,self.windowPane,num2)


  def findPositionMatches(self, peaks, dimNum, peakList=None):

    activePeakLists = []
    windowingPeakLists = []

    # only want peaks represented in the current window
    for view in self.getSpectrumViews():
      for winPeakList in view.windowPeakLists:
        try:
          peakList0 = winPeakList.analysisPeakList.peakList
        except:
          peakList0 = None
        if (peakList0 and not peakList0.isDeleted and dimNum < peakList0.dataSource.numDim):
          activePeakLists.append(peakList0)

    for peak in peaks:
      if peak.peakList not in activePeakLists:
        peaks.remove(peak)

    for peak in peaks:
      if peak.peakList not in windowingPeakLists:
        windowingPeakLists.append(peak.peakList)

    if peakList is None:
      searchPeakLists = activePeakLists
    else:
      searchPeakLists = [peakList,]

    dimDict = {}
    peakLists = set([p.peakList for p in peaks])
    peakLists.update(searchPeakLists)
    
    for peakList in peakLists:
      if dimNum < peakList.dataSource.numDim:
        dimDict[peakList] = dimNum+1

    groups = findPositionMatches(peaks, dimDict, searchPeakLists)

    self.topPopup.viewPeakGroups(groups, dimNum, windowingPeakLists, peaks)

  def makeIntermediatePeak(self, *event):

    peaks = self.topPopup.currentPeaks
    if peaks:
      makeIntermediatePeak(peaks)

  # TBD: functions which call this need updating to windowPane
  def gotoOrthogonalPlane(self, windowPane, position):

    windowPane.getWindowFrame().gotoPosition(position=position)


  def makeOrthogonalStrips(self, windowPane, positions):

    N = len(positions)
    if N > 10 and not showOkCancel('Warning','%d positions selected. Really make %d strips in window %s' % (N,N,getWindowPaneName(windowPane)), parent=self):
      return

    label = 'y'
    if self.windowPane.spectrumWindow.stripAxis == 'y':
      label = 'x'

    M = 0.0
    yPos = 0.0
    for xyz in positions:
      if xyz.get(label):
        yPos += xyz.get(label)
        M += 1.0

    if M > 0.0:
      yPos /= M
      for xyz in positions:
        xyz[label] = yPos

    if positions:
      displayStrips(self.topPopup, positions, orthoPositions=None,
                    spectrum=None, windowPane=windowPane)


  def predictSpinSystemType(self, *event):

    if self.menuPeak:
      spinSystem = None
      for peakDim in self.menuPeak.peakDims:
        for contrib in peakDim.peakDimContribs:
          if contrib.resonance.resonanceGroup:
            spinSystem = contrib.resonance.resonanceGroup
            break
        else:
          continue
        break

      if spinSystem:
        shiftList = self.menuPeak.peakList.dataSource.experiment.shiftList
        self.topPopup.typeSpinSystem(spinSystem=spinSystem, shiftList=shiftList)

  def makeSeqSpinSystemConnections(self, seqOffsets, *event):

    peaks = []
    if self.menuPeak and self.topPopup.currentPeaks:
      for peak in self.topPopup.currentPeaks:
        if peak.peakList is self.menuPeak.peakList:
          peaks.append(peak)

    elif self.menuPeak:
      peaks = [self.menuPeak,]

    for peak in peaks:
      addPeakResonancesToSeqSpinSystems(peak, seqOffsets)

  def showPeakStrips(self, *event):

    peaks = self.topPopup.currentPeaks
    if peaks:
      pane = self.windowPane
      axis = pane.spectrumWindow.stripAxis
      axisPanel = pane.findFirstAxisPanel(label=axis)
      
      sortList = []
      
      for peak in peaks:
        analysisSpectrum = peak.peakList.dataSource.analysisSpectrum  
        view = pane.findFirstSpectrumWindowView(analysisSpectrum=analysisSpectrum)
        
        if view:
          axisMapping = view.findFirstAxisMapping(label=axis)
 
          if axisMapping:
            dim = axisMapping.analysisDataDim.dataDim.dim
          else:
            dim = 1
        else:
          dim = 1
          
        peakDim = peak.findFirstPeakDim(dim=dim)
        sortList.append((peakDim.value, peak))
     
      sortList.sort()
      
      if axisPanel.axisUnit.isBackwards:
        sortList.reverse()
      
      peaks = [x[1] for x in sortList]
     
    elif self.menuPeak:
      peaks = [self.menuPeak,]

    if peaks:
      displayPeakStrips(self.topPopup, peaks, self.windowPane)


  def showStructConnections(self, *event):

    peaks = []
    if self.topPopup.currentPeaks:
      peaks = self.topPopup.currentPeaks
    elif self.menuPeak:
      peaks = [self.menuPeak,]

    if peaks:
      self.topPopup.viewStructure()
      popup = self.topPopup.popups['view_structure']
      popup.clearConnections()
      for peak in peaks:
        popup.showPeakConnection(peak)
      popup.updateAfter()

  def gotoOrthogonalAllPlanes(self, windowZplanes):

    for planeName, windowPane, positions in windowZplanes:
       self.gotoOrthogonalPlane(windowPane, positions[0])

  def peakCenterOrthogonalPlanes(self, event):

    canvas = event.widget
    (row, col) = self.scrolled_window.getCanvasRowCol(canvas)

    position = {}
    if self.menuPeak:
      peak = self.menuPeak

      for view in self.getSpectrumViews():
        if view.analysisSpectrum.dataSource is peak.peakList.dataSource:
          for axisMapping in view.axisMappings:
            if axisMapping.label not in ('x','y'):
              for peakDim in peak.peakDims:
                if peakDim.dataDim is axisMapping.analysisDataDim.dataDim:
                  dataDimRef = peakDim.dataDimRef
                  if dataDimRef:
                    position[axisMapping.label] = getPeakDimPosition(peakDim, toUnit=dataDimRef.expDimRef.unit)
                    break
                  else:
                    position[axisMapping.label] = peakDim.position

      self.gotoPosition(position, row=row, col=col)

  def translatePeak(self, event):

    peaks = self.topPopup.currentPeaks
    peak = self.menuPeak # this is supposed to be the referencePeak

    if len(peaks) != 2:
      showError('Error','Need to have two peaks selected (have %d)' % len(peaks), parent=self)
      return

    if peak:
      if peak not in peaks:
        showError('Error','Two peaks selected and menu peak is not one of them', parent=self)
        return

      if peak == peaks[1]:
        peaks = (peaks[1], peaks[0])

    (referencePeak, translatePeak) = peaks

    referenceSpectrum = referencePeak.peakList.dataSource
    translateSpectrum = translatePeak.peakList.dataSource

    if referenceSpectrum == translateSpectrum:
      showError('Error','Two peaks selected are in same spectrum', parent=self)
      return

    getView = self.windowPane.findFirstSpectrumWindowView
    referenceView = getView(analysisSpectrum=referenceSpectrum.analysisSpectrum)
    translateView = getView(analysisSpectrum=translateSpectrum.analysisSpectrum)

    popup = TranslatePeakPopup(self, referencePeak, translatePeak)
    popup.destroy()

  def zoomToSpinSystem(self, event):

    peak = self.menuPeak
    if peak:
      canvas = event.widget
      row, col = self.scrolled_window.getCanvasRowCol(canvas)

      specDict = {}
      for view in self.windowPane.spectrumWindowViews:
        specDict[view.analysisSpectrum.dataSource] = 1

      ssDict = {}
      for peakDim in peak.peakDims:
        for contrib in peakDim.peakDimContribs:
          ss = contrib.resonance.resonanceGroup
          if ss:
            ssDict[ss] = ssDict.get(ss, 0) + 1

      if ssDict:
        ssBest = ssDict.keys()[0]
        for ss in ssDict.keys():
          if ssDict[ss] > ssDict[ssBest]:
            ssBest = ss

        peaks = []
        for resonance in ssBest.resonances:
          for contrib in resonance.peakDimContribs:
            peak = contrib.peakDim.peak
            if specDict.get(peak.peakList.dataSource):
              peaks.append(peak)

        if peaks:
          zoomToShowPeaks(peaks, self.windowPane, row, col)

  def zoomToSelectedPeaks(self, event):

    peaks  = self.topPopup.currentPeaks
    if peaks:
      canvas = event.widget
      row, col    = self.scrolled_window.getCanvasRowCol(canvas)
      zoomToShowPeaks(peaks, self.windowPane, row, col)

  def autoArrangePeakLabels(self, event):

    if self.waitDraw:
      return

    peaks  = self.topPopup.currentPeaks
    if peaks:
      refinePeakLabelPositions(peaks, self.windowPane, self.topPopup)

  def resetPeakLabels(self, event):

    peaks  = self.topPopup.currentPeaks
    if peaks:
      resetPeakLabelPositions(peaks, self.windowPane)
      
  def setScrollbarVisibility(self, axes):
  
    for axisPanel in self.windowPane.axisPanels:
      
      for axis in axes:
        if axis in axisPanel.label: # e.g. 'z' in 'z1','z2'
          axisPanel.isVisible = True
          break
      
      else: 
        axisPanel.isVisible = False    
        
  def setCrosshairTraces(self, axes):
  
    window = self.windowPane.spectrumWindow
    window.isXSliceDrawn = 'x' in axes
    window.isYSliceDrawn = 'y' in axes

  def calcWorldCoord(self, canvas, x, y):

    # input: x, y are in pixels
    # output: (a, b, s, t) where a, b are in ppm and s, t are proportion of width/height

    return self.scrolled_window.calcWorldCoord(canvas, x, y)

  #
  # turnDrawRequestsOff / turnDrawRequestsOn is for cases where you change a lot of
  # data model attributes but don't want a re-drawing to occur until at the end
  # turnDrawRequestsOff returns False if someone else has already turned drawing
  # requests off, otherwise it returns True
  # You should only call turnDrawRequestsOn if turnDrawRequestsOff returned True
  # Code should look like:
  #
  # turnedOff = self.turnDrawRequestsOff()
  # try:
  #   ... (lots of code with API model changes)
  # finally:
  #   if turnedOff: self.turnDrawRequestsOn()
  #

  def turnDrawRequestsOff(self):

    if self.waitDraw:
      return False

    self.waitDraw = True

    return True

  def turnDrawRequestsOn(self, doDraw=True):

    self.waitDraw = False

    if doDraw:
      self.drawAll()

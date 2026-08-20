
"""
======================COPYRIGHT/LICENSE START==========================

EditWindow.py: Part of the CcpNmr Analysis program

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
from memops.general import Implementation

from memops.universal.Util import formatFloat

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.DataEntry import askString
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.LabelDivider import LabelDivider
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError, showYesNo
from memops.gui.MultiWidget import MultiWidget
from memops.gui.PulldownList import PulldownList
from memops.gui.RadioButtons import RadioButtons
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.ToggleLabel import ToggleLabel

from ccpnmr.analysis.core import Util
from ccpnmr.analysis.core import WindowBasic
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.popups.ChangeAxisMapping import ChangeAxisMappingPopup
from ccpnmr.analysis.popups.CreateAxis import CreateAxisPopup

SLICE_HEIGHT_VALS = [50, 100, 150, 200]
ACTIVE_OPTIONS = ['All windows', 'Active windows']

def hasVisibleSlicePanel(axisPanel):

  if axisPanel.label not in ('x', 'y'):
    return False

  if WindowBasic.windowPaneHasValueAxis(axisPanel.spectrumWindowPane):
    return False

  return True

class WindowListScrolledMatrix(ScrolledMatrix):

  def doEditMarkExtraRules(self, windowPane, row, col):

    # DANGEROUS: make sure keeps in step with columns below
    if col not in (3, 4, 8):
      return True

    hasValueAxis = WindowBasic.windowPaneHasValueAxis(windowPane)
    if col in (3, 4):
      return not hasValueAxis
    else: # col = 8
      return hasValueAxis

class AxisListScrolledMatrix(ScrolledMatrix):

  def doEditMarkExtraRules(self, axisPanel, row, col):

    # DANGEROUS: make sure keeps in step with columns below
    if col < 4:
      return True

    if col == 6:
      return axisPanel.axisType.name in ('1H', '13C', '15N', '2H', '19F', '31P', '79Br')

    return hasVisibleSlicePanel(axisPanel)

class EditWindowPopup(BasePopup):
  """
  **Edit Spectrum Window Parameters**
  
  This popup window controls all of the settings that relate to the spectrum
  display windows used in Analysis. All of the spectrum windows within the
  project are listed and all of the details which control how they show spectra
  and peak lists are given, and in many cases the parameters can be edited. The
  settings available in this popup are often also, and sometimes more commonly,
  manipulated elsewhere. For example, the toggle buttons at the top of spectrum
  windows are usually used to control spectrum and peak list visibility.
  However, this popup window provides an authoritative list of all relevant
  settings, including a few that are not accessible from elsewhere.

  The layout of the popup is divided into three tabs; the first lists the main
  properties of the spectrum windows and their axis, the second controls how
  spectra and peak lists are mapped into the windows and displayed, and the
  third allows the user to group windows together so that different sets may be
  reserved for different kinds of operation. 

  **Windows & Axes**
  
  The first tab is divided into two sections. The upper table lists the spectrum
  windows which have been defined in the current CCPN project or, if the "Active
  windows" options is set, those in the active window group. The settings given
  in the table, which may all be edited, relate to general properties that
  govern the whole window display. Typically user may change the name of a
  window and adjust the aspect ratio; all spectrum window displays in Analysis
  preserve aspect ratio when zooming, so peaks remain the same aesthetic shape.
  The crosshair trace settings control whether a 1D intensity slice is
  superimposed upon the cursor location, although this is normally toggled on or
  off via the right-click mouse menu of the particular spectrum window. For both
  the spectrum traces that appear at the crosshair location and those that may
  be added to side panels (at the window's edge) the "Slice Range" column
  dictates the scale that is used to show the intensity axis. In normal
  operation the user controls this with the <PgUp> and <PgDn> keys, rather than
  via this table. The last two columns of the table indicate whether strip
  numbers and Z axis positions are superimposed on the spectrum display (in the
  corners).

  The lower table lists the kinds of axes are present in the selected window;
  the selected row in the upper table. Here the user can control whether a
  scrollbar to move window contents appears or not, e.g. to save some screen
  space. The "Panel Type" column is used to make sub-groups of window axes that
  otherwise have the same isotope type. The panel type definitions are set
  elsewhere, in the Axes_ popup. The panel type is only relevant for dictating 
  where 1D ruler lines go in different windows; a ruler is only displayed on
  axes with the same panel type, thus the user can distinguish amide 1H axes
  (panel type "H1") and other 1H axes (panel type "H2"). The last two columns
  adjust the depth and presence of the "Side Trace". When visible, a side trace
  shows a 1D intensity graph at the edge of the window axis. Naturally, this
  only applies for the X & Y screen dimensions.

  **Spectrum & Peak List Mappings**
  
  The second tab controls how spectra and peak lists are displayed within the
  spectrum window selected at the top. It should be noted that spectrum contour levels and
  colours are set elsewhere, via the main Spectra_ popup. Likewise, the peak
  list symbols and colours are changed in the `Peak Lists`_ popup. These
  peaklist and spectrum settings are set on a project-wide basis, and thus are
  not set for an individual display window. 

  The upper table lists all of the spectra that may be displayed in the selected
  spectrum window. Here the user can control: whether a spectrum is listed in
  the  toolbar at the top of the window, where spectra are switched on and off;
  whether the positive contours are drawn; whether the negative contours are
  drawn; and whether the 1D intensity slices, when used, include the spectrum.
  The last "Dim. Mapping" column is important if the user wishes to change how
  the data dimensions of a spectrum are mapped to the X, y & Z axes of the
  spectrum window, especially given that Analysis cannot always guess this
  correctly when there are two axis with the same isotope type. If, for example,
  the user has a 3D HC_cH.TOCSY spectrum and the first 1H data dimension (from
  the spectrum file) is displayed on the X axis, when it should really be on the
  Z axis, then the user clicks on the "Dim. Mapping" column. In the resulting
  mini-popup, changing the "Dimension 1" pulldown menu from "x" to "z1" will
  correct the axis mapping and effectively rotate the spectrum data as far as
  the window view is concerned. It should be noted that multiple windows may
  need to be adjusted in this way for any given spectrum; only the user can
  really decide what is correct, especially where data axes have been swapped.

  The lower "Mapped Peak Lists" table shows all of the peak lists that may be
  displayed in the selected spectrum window. In the display sense, a peak list
  consists of a list of peak symbols (i.e. positioned crosses) and assignment
  annotations. The two columns of the table that can be edited dictate whether a
  peak list's peak symbols and annotations are drawn, thus giving total control.
  However, these settings are usually controlled via the toolbar at the top of a
  spectrum window, and normally the peak symbols and textual annotations are both
  switched on or off at the same time. 

  **Window Groups**
 
  Defining window groups is a mechanism by which different selections of
  spectrum windows may be made so that only some are visible at the same time.
  There will only be one active group of spectrum windows at a given time, and
  only those in this active group can be opened and will be shown in the main
  Analysis menu. The user can create as many named window groupings as is
  required and can add any window to any group. Clicking on a group in the upper
  table lists all of the windows in the project and whether they are included in
  the selected group or not. To add a window to a group, or remove it from a
  group, the user simply toggles the "In Group?" column. The selected group may
  be set as the active one by clicking on the central [Make Active] button.

  .. _Axes: EditAxisPanelPopup.html
  .. _Spectra: EditSpectrumPopup.html
  .. _`Peak Lists`: EditPeakListsPopup.html
  """

  def __init__(self, parent, *args, **kw):

    self.windowPane = None
    self.axisPanel = None
    self.groupWindow = None
    self.windowGroup = None
    self.mapWindowPane = None
    self.waitingM = False
    self.specWindowView = None
    self.winPeakList = None

    BasePopup.__init__(self, parent=parent, title='Window : Windows', **kw)

  def body(self, guiFrame):

    self.geometry('600x700')  

    guiFrame.expandGrid(0,0)
   
    tipTexts = ['Shows a table of spectrum contour windows and details of their axes',
                'Sets how spectra & peak lists are drawn in contour windows, including axis mappings; note colours are set elsewhere in peak list & spectrum tables',
                'Allows spectrum windows to be grouped into sets so that only some are displayed at a given time']
    options = ['Windows & Axes', 'Spectrum & Peak List Mappings', 'Window Groups']
      
    tabbedFrame = TabbedFrame(guiFrame, options=options,
                              grid=(0,0), tipTexts=tipTexts)

    self.tabbedFrame = tabbedFrame
    frameA, frameB, frameC, = tabbedFrame.frames

    #
    # Windows & Axes
    #
    
    frameA.grid_columnconfigure(0, weight=1)
    frameA.grid_rowconfigure(2, weight=1)
    frameA.grid_rowconfigure(5, weight=1)

    row = 0
    frame = Frame(frameA, grid=(row,0))
    label = Label(frame, text=' Windows:', grid=(0,0))
    
    tipTexts = ['For the below table, show all of the spectrum windows in the current project, irrespective of the active group',
                'For the below table, show only spectrum windows that are in the currently active window group']
    self.which_windows = RadioButtons(frame, entries=ACTIVE_OPTIONS, tipTexts=tipTexts,
                                      select_callback=self.updateWindowList)
    self.which_windows.grid(row=0, column=1, sticky='w')

    row += 1
    div = LabelDivider(frameA, text ="Edit Windows", grid=(row,0))
    
    row += 1
    self.windowTable = None
    tipTexts = ['Row number',
                'Short textual name to identify the window',
                'The screen (X-Y) pixel:unit aspect ratio of the window, defaults according to isotope; "8.0" means Y axis shows 8.0 units for every 1.0 on X',
                'Whether to show a 1D slice of the specta at the cursor position, along the horizontal screen axis',
                'Whether to show a 1D slice of the specta at the cursor position, along the vertical screen axis',
                'The upper and lower intensity bounds for the 1D cursor slice display; used to scale the magnitude of the drawn amplitude',
                'Whether to display strip numbers, in the upper left of each window subdivision',
                'Whether to display the Z-axis positions of each window strip, in the lower left corner of each window subdivision',
                ] # 'Whether to display the zero line in 1D windows']
    headings = ('#', 'Name', 'Aspect\nRatio',
                'X Crosshair\nTrace', 'Y Crosshair\nTrace',
                'Slice Intensity\nRange', 'Strip\nNumbers',
                'Strip Z\nLocations') #, 'Show\nZero Line')
    self.windowNameEntry = Entry(self, returnCallback=self.setWindowName, width=10)
    self.aspectWidget = FloatEntry(self, returnCallback=self.setAspectRatio, width=5)
    self.sliceWidget = FloatEntry(self, returnCallback=self.setSliceRange, isArray=True, width=24)
    
    editWidgets = [None, self.windowNameEntry, self.aspectWidget,
                   None, None, self.sliceWidget, None, None] #, None]
                   
    editGetCallbacks = [None, self.getWindowName, self.getAspectRatio,
                        self.toggleXSliceVisible, self.toggleYSliceVisible,
                        self.getSliceRange, self.toggleLabelVisible,
                        self.toggleMidpointVisible] #, self.toggleZeroLineVisible]
                         
    editSetCallbacks = [None, self.setWindowName, self.setAspectRatio,
                        None, None, self.setSliceRange, None, None] #, None]
                         
    self.windowTable = WindowListScrolledMatrix(frameA, headingList=headings,
                                                callback=self.selectWindowPane,
                                                multiSelect=True,
                                                editWidgets=editWidgets,
                                                editGetCallbacks=editGetCallbacks,
                                                editSetCallbacks=editSetCallbacks,
                                                deleteFunc=self.deleteWindow,
                                                grid=(row,0), tipTexts=tipTexts)
  
    row += 1
    tipTexts = ['Open and display the selected spectrum window',
                'Open a tool to create a new spectrum window; specifying names and axes etc. ',
                'Make a copy of the currently selected window, which may then be independently manipulated',
                'Delete the currently selected spectrum window; does not affect the loaded spectra']
    texts = [ 'Open', 'Create New', 'Clone', 'Delete' ]
    commands = [self.openWindow, self.createWindow,
                self.cloneWindow, self.deleteWindow]
    self.window_buttons = ButtonList(frameA, texts=texts, commands=commands,
                                     grid=(row,0), tipTexts=tipTexts)
    
    row += 1
    div = LabelDivider(frameA, text ="Window Axes")
    div.grid(row=row,column=0,sticky='ew')

    row += 1
    tipTexts = ['A short label for identifying the axis within its window, e.g. "x", "y", "z1", "z2"',
                'A name for the kind of thing represented along this axis; usually serves to subdivide frequency axes into different isotope types',
                'For a given kind of axis, which sub-type is currently set; used to differentiate axes with the same isotope ',
                'Whether the scrollbar, which shifts the spectrum view position, is displayed for the window axis',
                'For an X or Y axis, whether a 1D slice trace (at the cursor position) is displayed alongside the main contour display',
                'For an X or Y axis, sets the number of pixels use for the amplitude display of 1D slice traces (using cursor position)',
                'Axes in other windows that are tied (synchronized) to this axis']
    headings = ('Axis\nLabel', 'Axis\nType', 'Panel\nType',
                'Scrollbar\nVisible', 'Side Trace\nVisible',
                'Side Trace\nHeight', 'Tied Axes in\nOther Windows')
                
    self.panelTypeWidget = PulldownList(self, callback=self.setAxisPanelType)
                                         
    self.sliceThicknessWidget = PulldownList(self, callback=self.setSliceThickness,
                                             texts=map(str, SLICE_HEIGHT_VALS),
                                             objects=SLICE_HEIGHT_VALS)
                                             
    self.tiedAxesWidget = MultiWidget(self, CheckButton, callback=self.setTiedAxes, minRows=0, useImages=False)

    editWidgets = [None, None, self.panelTypeWidget,
                   None, None, self.sliceThicknessWidget, self.tiedAxesWidget]
                   
    editGetCallbacks = [None, None, self.getPanelType, self.toggleScrollbarVisible,
                        self.toggleSliceVisible, self.getThickness, self.getTiedAxes]
                         
    editSetCallbacks = [None, None, self.setPanelType,
                        None, None, self.setThickness, self.setTiedAxes]
                        
    self.axisTable = AxisListScrolledMatrix(frameA, headingList=headings,
                                            initialRows=4,
                                            callback=self.selectAxis,
                                            editWidgets=editWidgets,
                                            editGetCallbacks=editGetCallbacks,
                                            editSetCallbacks=editSetCallbacks,
                                            deleteFunc=self.deleteAxis,
                                            grid=(row,0), tipTexts=tipTexts)

    row += 1
    tipTexts = ['For the currently selected spectrum window, add a new axis',]
                #'Delete the selected spectrum window axis']
    texts = ['New Axis',] # 'Delete']
    commands = [ self.createAxis, ] # self.deleteAxis ]
    self.axis_buttons = ButtonList(frameA, texts=texts, tipTexts=tipTexts,
                                   commands=commands, grid=(row,0))
  
    #
    # Mappings
    #
    
    frameB.grid_columnconfigure(1, weight=1)
    
    row = 0
    frame = Frame(frameB, grid=(row, 0))
    
    label = Label(frame, text='Spectrum Window: ', grid=(0,0))

    tipText = 'Selects the spectrum window to display the mappable spectra & peak lists for'
    self.mapWindowPulldown = PulldownList(frame, callback=self.changeMapWindow,
                                          grid=(0,1), tipText=tipText)

    texts = ('All On', 'All Off', 'Invert')
    commands = (lambda: self.setAllSpectraVisibility(True), lambda: self.setAllSpectraVisibility(False), self.invertSpectraVisibility)
    tipTexts = ('Turn on all spectra (contours, etc.) and peak lists for this window',
                'Turn off all spectra (contours, etc.) and peak lists for this window',
                'Invert visibility of spectra (contours, etc.) and peak lists for this window')
    buttons = ButtonList(frameB, texts=texts, commands=commands, tipTexts=tipTexts, grid=(row, 1), sticky='e')
    
    row += 1
    div = LabelDivider(frameB, text='Mapped Spectra', grid=(row,0), gridSpan=(1,2))

    row += 1
    frameB.grid_rowconfigure(row, weight=1)
    tipTexts = ['The name of the experiment record for the window-mappable spectrum',
                'The name of the window-mappable spectrum',
                'Whether the spectrum appears with its own toggle button in the toolbar at the top of the spectrum window',
                'Whether the positive contour levels are shown in the window for the spectrum',
                'Whether the negative contour levels are shown in the window for the spectrum',
                'Whether 1D slices are visible for the spectrum, if any slice displays are active for the window',
                'Whether contour line in 1D window is visible for the spectrum; if on this helps indicate where peaks would be picked',
                'Sets how the spectrum dimensions map to corresponding window axes; for a given kind, enables swapping of spectrum dimensions relative to the window',
                'How much the spectrum is offset in the x axis', 'How much the spectrum is offset in the y axis']
    headings = ('Experiment', 'Spectrum', 'Spectrum In\nToolbar',
                'Pos. Contours\nVisible', 'Neg. Contours\nVisible',
                'Slice\nVisible', 'Contour Line\nVisible', 'Dim. Mapping', 'X Offset', 'Y Offset')
                
    self.xAxisOffsetWidget = FloatEntry(self, returnCallback=self.setXAxisOffset, width=5)
    self.valueAxisOffsetWidget = FloatEntry(self, returnCallback=self.setValueAxisOffset, width=5)
    
    editWidgets = editSetCallbacks = 10 * [None] # changed later depending on whether 1D or ND window
    editGetCallbacks = [None, None, self.toggleInToolbar,
                        self.togglePosContoursVisible, self.toggleNegContoursVisible,
                        self.toggleViewSliceVisible, self.toggleContourLineVisible,
                        self.changeViewMapping, self.getXAxisOffset, self.getValueAxisOffset ]
    self.spectraMatrix = ScrolledMatrix(frameB, headingList=headings,
                                        callback=self.selectSpectrum,
                                        initialRows=5, editWidgets=editWidgets,
                                        editGetCallbacks=editGetCallbacks,
                                        editSetCallbacks=editSetCallbacks,
                                        grid=(row,0), gridSpan=(1,2), tipTexts=tipTexts)

    row += 1
    tipTexts = ['Opens a table allowing editing of per-spectrum (rather than window) parameters',]
    texts = ['Edit Spectra']
    commands = [ self.editSpectrum ]
    self.spectrumButtons = ButtonList(frameB, texts=texts, commands=commands,
                                      grid=(row,0), gridSpan=(1,2), tipTexts=tipTexts)

    row += 1
    div = LabelDivider(frameB, text='Mapped Peak Lists', grid=(row,0), gridSpan=(1,2))
    
    row += 1
    frameB.grid_rowconfigure(row, weight=1)
    
    tipTexts = ['The name of the experiment record for the window-mappable peak list',
                'The name of the window-mappable spectrum that contains the peak list',
                'The serial number of the peak list, which is displayable in the spectrum window',
                'Whether the peak crosses/symbols are drawn for the peak list in the selected window; to indicate the position of extrema/resonances',
                'Whether the peak annotation texts are drawn for the peak list in the selected window']
    headings = ('Experiment', 'Spectrum', 'Peak list',
                'Symbol\nDrawn?', 'Annotation\nDrawn?')
    editWidgets = 5 * [None]
    editGetCallbacks = [None, None, None,
                        self.togglePeakSymbol,
                        self.togglePeakAnnotation ]
    editSetCallbacks = 5 * [None]
    
    self.peakListMatrix = ScrolledMatrix(frameB, headingList=headings,
                                         callback=self.selectPeakList,
                                         initialRows=5, editWidgets=editWidgets,
                                         editGetCallbacks=editGetCallbacks,
                                         editSetCallbacks=editSetCallbacks,
                                         grid=(row,0), gridSpan=(1,2),
                                         tipTexts=tipTexts)

    row += 1
    tipTexts = ['Show a table of the peaks within the selected peak list',]
    texts = ['Edit Peak List',]
    commands = [self.editPeakList,]
    self.peakListButtons = ButtonList(frameB, texts=texts, commands=commands,
                                       grid=(row,0), gridSpan=(1,2), tipTexts=tipTexts)
    
    #
    # Groups
    #
    
    frameC.grid_columnconfigure(0, weight=1)
    frameC.grid_rowconfigure(1, weight=1) 
    frameC.grid_rowconfigure(4, weight=1)
    
    row = 0
    div = LabelDivider(frameC, text ="Window Groups", grid=(row,0))
    
    row += 1
    self.windowGroupNameWidget = Entry(self, width=10,
                                 returnCallback=self.setGroupName)
                                 
    tipTexts = ['Row number',
                'The name for the group of spectrum windows',
                'The number of windows contained within the group']
    headings = ('#', 'Group', 'Number of windows')
    editWidgets = [ None, self.windowGroupNameWidget, None ]
    editGetCallbacks = [ None, self.getGroupName, None ]
    editSetCallbacks = [ None, self.setGroupName, None ]
    
    self.windowGroupTable = ScrolledMatrix(frameC, headingList=headings,
                                     initialRows=8, callback=self.selectGroup,
                                     editWidgets=editWidgets, multiSelect=True,
                                     editGetCallbacks=editGetCallbacks,
                                     editSetCallbacks=editSetCallbacks,
                                     deleteFunc=self.deleteGroup,
                                     grid=(row,0), tipTexts=tipTexts)

    row += 1
    tipTexts = ['Sets the selected window group as the active one; only windows within this group will be displayable (in menus etc.)',
                'Make a new grouping of spectrum windows, which can be used to quickly swap between different displayed sets',
                'Delete the selected grouping from the CCPN project; does not otherwise affect the contained windows']
    texts = [ 'Make Active', 'New Group', 'Delete Selected' ]
    commands = [ self.openGroup, self.createGroup, self.deleteGroup ]
    self.windowGroup_buttons = ButtonList(frameC, texts=texts, commands=commands,
                                          grid=(row,0), tipTexts=tipTexts)
    
    row += 1
    div = LabelDivider(frameC, text ="Group Membership", grid=(row,0))
    
    row += 1
    tipTexts = ['Row number',
                'The name of the spectrum window, which may or may not be in the selected group',
                'Whether the spectrum window is part of the selected grouping; can be toggled to include & exclude']
    headings = ('#', 'Window', 'In Group?')
    self.windowNameEntry2 = Entry(self, returnCallback=self.setWindowName2, width=10)
    editWidgets = [ None, self.windowNameEntry2, None ]
    editGetCallbacks = [ None, self.getWindowName2, self.getInGroup ]
    editSetCallbacks = [ None, self.setWindowName2, None ]
    self.groupWindowsTable = ScrolledMatrix(frameC, headingList=headings,
                                            initialRows=8,
                                            callback=self.selectWindow,
                                            editWidgets=editWidgets,
                                            editGetCallbacks=editGetCallbacks,
                                            editSetCallbacks=editSetCallbacks,
                                            grid=(row,0), gridSpan=(1,2),
                                            tipTexts=tipTexts)
  
    #
    # Main
    #
  
    buttons = UtilityButtonList(tabbedFrame.sideFrame, helpUrl=self.help_url,
                                grid=(0,0), sticky='e')

    self.updateWindowList()
    self.setAxisState()
    self.updateGroupList()
    self.setGroupState()
    self.updateMapWindowPulldown()
    self.updateMappingsAfter()

    self.administerNotifiers(self.registerNotify)

  def administerNotifiers(self, notifyFunc):

    # classes
    
    spectrumWindow = 'ccpnmr.Analysis.SpectrumWindow'
    spectrumWindowPane = 'ccpnmr.Analysis.SpectrumWindowPane'
    axisPanel = 'ccpnmr.Analysis.AxisPanel'
    axisRegion = 'ccpnmr.Analysis.AxisRegion'
    windowGroup = 'ccpnmr.Analysis.SpectrumWindowGroup'
    slicePanel = 'ccpnmr.Analysis.SlicePanel'
    
    # Windows tab

    for func in ('__init__', 'delete'):
      notifyFunc(self.delayedUpdateWindowList, spectrumWindow, func)
      notifyFunc(self.delayedUpdateAxisList, axisPanel, func)
      notifyFunc(self.delayedUpdateRegionList, axisRegion, func)
      notifyFunc(self.delayedUpdateRegionList, axisPanel, func)

    notifyFunc(self.updateWindowList, spectrumWindow, '')
    notifyFunc(self.updateWindowList, spectrumWindowPane, '')
    notifyFunc(self.delayedUpdateRegionList, axisRegion, 'setRegion')
    notifyFunc(self.delayedUpdateRegionList, axisRegion, 'setSize')
    notifyFunc(self.updateAxisList, axisPanel, '')
    notifyFunc(self.updateAxisList, slicePanel, '')

    for func in ('addSpectrumWindowGroup', 'removeSpectrumWindowGroup'):
      notifyFunc(self.updateAllLists, spectrumWindow, func)

    for func in ('delete', 'addSpectrumWindow', 'removeSpectrumWindow'):
      notifyFunc(self.updateAllLists, windowGroup, func)

    # Mappings tab
    
    for func in ('__init__', 'delete', ''):
      self.registerNotify(self.updateMappingsAfter, 'ccpnmr.Analysis.SpectrumWindowView', func)
      notifyFunc(self.updateMappingsAfter, 'ccpnmr.Analysis.WindowPeakList', func)

    notifyFunc(self.updateMappingsAfter, 'ccpnmr.Analysis.SpectrumWindow', '')
    notifyFunc(self.updateMappingsAfter, 'ccp.nmr.Nmr.DataSource', 'setName')
    notifyFunc(self.updateMappingsAfter, 'ccp.nmr.Nmr.Experiment', 'setName')
    notifyFunc(self.updateMappingsAfter, 'ccpnmr.Analysis.AxisMapping', 'setDim')
    
    # Groups tab

    for func in ('__init__', 'delete', 'setName',
                 'addSpectrumWindow', 'removeSpectrumWindow'):
      notifyFunc(self.updateGroupList, windowGroup, func)

    for func in ('delete', 'addSpectrumWindow', 'removeSpectrumWindow'):
      notifyFunc(self.updateGroupWindows, windowGroup, func)

    for func in ('addSpectrumWindowGroup', 'removeSpectrumWindowGroup'):
      notifyFunc(self.updateGroupList, spectrumWindow, func)
      notifyFunc(self.updateGroupWindows, spectrumWindow, func)

    for func in ('__init__', 'delete', 'setName'):
      notifyFunc(self.updateGroupWindowsAfter, spectrumWindow, func)
      notifyFunc(self.updateGroupWindowsAfter, spectrumWindowPane, func)

    notifyFunc(self.updateGroupList, 'ccpnmr.Analysis.AnalysisProject', 'setActiveWindowGroup')

    #notifyFunc(self.changedAppData, 'memops.Implementation.AppDataBoolean', 'setValue')

  def open(self):
  
    self.updateWindowList()
    self.setAxisState()
    self.updateGroupList()
    self.setGroupState()
    self.updateMapWindowPulldown()
    self.updateMappingsAfter()
    BasePopup.open(self)
 
  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)

  def updateAllLists(self, *extra):

    self.updateWindowList()
    self.updateAxisList()

  def delayedUpdateWindowList(self, *extra):

    self.after_idle(self.updateWindowList)

  def getWindowPanes(self):
  
    which = self.which_windows.get()
    if which == ACTIVE_OPTIONS[0]:
      windows = self.analysisProject.spectrumWindows
    else:
      windows = self.parent.getActiveWindows()
    windows = sorted(windows, key=lambda window: window.name)
    
    panes = []
    for window in windows:
      panes.extend(window.sortedSpectrumWindowPanes())

    return panes

  def updateWindowList(self, *extra):

    self.updateMapWindowPulldown()

    windowPanes = self.getWindowPanes()
      
    textMatrix = []
    
    for n, windowPane in enumerate(windowPanes):
      window = windowPane.spectrumWindow
      
      if WindowBasic.windowPaneHasValueAxis(windowPane):
        xSlice = ySlice = ''
      else:
        xSlice = window.isXSliceDrawn and 'Yes' or 'No'
        ySlice = window.isYSliceDrawn and 'Yes' or 'No'
      
      sliceMin = formatFloat(windowPane.sliceRange[0])
      sliceMax = formatFloat(windowPane.sliceRange[1])
      
      datum = [n+1, WindowBasic.getWindowPaneName(windowPane),
              formatFloat(windowPane.aspectRatio),
              xSlice, ySlice,
              sliceMin + ',' + sliceMax,
              window.isCanvasLabelShown and 'Yes' or 'No',
              window.isCanvasMidpointShown and 'Yes' or 'No']

      """ not needed for now
      if WindowBasic.windowPaneHasValueAxis(windowPane):
        datum.append(window.isZeroLineShown and 'Yes' or 'No')
      else:
        datum.append(None)
"""
              
      textMatrix.append(datum)

    self.windowTable.update(objectList=windowPanes,
                           textMatrix=textMatrix)

    self.setWindowState(self.windowPane)

  def setWindowState(self, windowPane=None):

    if windowPane and WindowBasic.isActiveWindow(windowPane.spectrumWindow):
      state = 'normal'
    else:
      state = 'disabled'
 
    button = self.window_buttons.buttons[0]
    button.config(state=state)

    if windowPane:
      state = 'normal'
    else:
      state = 'disabled'
 
    for n in (2, 3):
      button = self.window_buttons.buttons[n]
      button.config(state=state)

  def selectWindowPane(self, windowPane, row, col):

    self.windowPane = windowPane
    self.initWindowPane(windowPane)

  def initWindowPane(self, windowPane=None):

    self.setWindowState(windowPane)
    self.setAxisState()
    #self.setRegionState()
    self.updateAxisList()

  def getWindowName(self, windowPane):
 
    self.windowNameEntry.set(windowPane.spectrumWindow.name)

  def setWindowName(self, *extra):
 
    name = self.windowNameEntry.get()
    if self.windowPane.spectrumWindow.name == name:
      return
    
    if name in [window.name for window in self.analysisProject.spectrumWindows ]:
      showError('Setting window name', 'Window name "%s" already used' % name, parent=self)
      return

    try:
      self.windowPane.spectrumWindow.name = name
      
    except Implementation.ApiError, e:
      showError('Setting window name', e.error_msg, parent=self)

  def selectWindow(self, window, row, col):

    self.window = window

  def getWindowName2(self, spectrumWindow):
 
    self.windowNameEntry2.set(spectrumWindow.name)

  def setWindowName2(self, *extra):
 
    name = self.windowNameEntry2.get()
    if self.window.name == name:
      return
    
    if name in [window.name for window in self.analysisProject.spectrumWindows ]:
      showError('Setting window name', 'Window name "%s" already used' % name, parent=self)
      return

    try:
      self.window.name = name
      
    except Implementation.ApiError, e:
      showError('Setting window name', e.error_msg, parent=self)

  def getAspectRatio(self, windowPane):
 
    self.aspectWidget.set(windowPane.aspectRatio)

  def setAspectRatio(self, *extra):
 
    windowPane = self.windowPane
    
    try:
      aspectRatio = self.aspectWidget.get()
      if (aspectRatio <= 0):
        raise ValueError
      windowPane.aspectRatio = aspectRatio
      
    except ValueError:
      msg = 'Aspect ratio must be positive float'
      showError('Setting aspect ratio', msg, parent=self)
      
    except Implementation.ApiError, e:
      showError('Setting aspect ratio', e.error_msg, parent=self)

  def getSliceRange(self, windowPane):

    self.sliceWidget.set(windowPane.sliceRange)

  def setSliceRange(self, *extra):
 
    windowPane = self.windowPane
    try:
      sliceRange = self.sliceWidget.get()
      if len(sliceRange) != 2 or (sliceRange[0] >= sliceRange[1]):
        raise ValueError
        
      windowPane.sliceRange = sliceRange
      
    except ValueError:
      msg = 'Slice range must be two floats with first less than second'
      showError('Setting slice range', msg, parent=self)
   
    except Implementation.ApiError, e:
      showError('Setting slice range', e.error_msg, parent=self)

  def toggleXSliceVisible(self, windowPane):

    window = windowPane.spectrumWindow
    window.isXSliceDrawn = not window.isXSliceDrawn
    
  def toggleYSliceVisible(self, windowPane):

    window = windowPane.spectrumWindow
    window.isYSliceDrawn = not window.isYSliceDrawn
    
  def toggleLabelVisible(self, windowPane):

    window = windowPane.spectrumWindow
    window.isCanvasLabelShown = not window.isCanvasLabelShown
    
  def toggleMidpointVisible(self, windowPane):

    window = windowPane.spectrumWindow
    window.isCanvasMidpointShown = not window.isCanvasMidpointShown
    
  """ not needed for now
  def toggleZeroLineVisible(self, windowPane):

    window = windowPane.spectrumWindow
    window.isZeroLineShown = not window.isZeroLineShown
"""
    
  def openWindow(self, *extra):

    windowPane = self.windowPane

    if windowPane:
      self.parent.openWindow(windowPane.spectrumWindow)

  def createWindow(self):

    self.parent.newWindow()

  def cloneWindow(self):

    windowPane = self.windowPane

    if windowPane:
      window = windowPane.spectrumWindow
      self.parent.cloneWindow(window)

  def deleteWindow(self, *event):

    windowPanes = self.windowTable.currentObjects

    if len(windowPanes) == 1:
      windowPane = windowPanes[0]
      name = WindowBasic.getWindowPaneName(windowPane)
      msg = 'Are you sure you want to delete window "%s"?' % name
    
    else:
      names =  ', '.join([WindowBasic.getWindowPaneName(p) for p in windowPanes])
      msg = 'Are you sure you want to delete'
      msg += '%d windows [ %s]?' % (len(windowPanes),names)
            
    if showYesNo('Confirm', msg, parent=self):
      for windowPane in windowPanes:
        window = windowPane.spectrumWindow
        
        if len(window.spectrumWindowPanes) == 1:
          windowPane.delete()
          window.delete()
        
        else:
          windowPane.delete()
        
      self.updateWindowList()
      self.initWindowPane()


  def delayedUpdateAxisList(self, *extra):

    self.after_idle(self.updateAxisList)

  def updateAxisList(self, *extra):

    windowPane = self.windowPane

    if windowPane:
      axisPanels = [ap for ap in windowPane.sortedAxisPanels() if not ap.isDeleted]
    else:
      axisPanels = []

    textMatrix = []
    for axisPanel in axisPanels:
      # axisPanel.axisUnit.unit
      text = [axisPanel.label,
              axisPanel.axisType.name,
              axisPanel.panelType.name,
              axisPanel.isVisible and 'Yes' or 'No']
              
      if hasVisibleSlicePanel(axisPanel):
        text.append(axisPanel.slicePanel.isVisible and 'Yes' or 'No')
        thickness = Util.greaterOrEqualEntry(axisPanel.slicePanel.thickness,
                                             SLICE_HEIGHT_VALS)
        text.append(thickness)
      
      else:
        text.append(None)
        text.append(None)
        
      tiedAxisPanels = WindowBasic.getTiedAxisPanels(axisPanel)
      tiedNames = self._getTiedNames(tiedAxisPanels)
      text.append(','.join(tiedNames))

      textMatrix.append(text)

    self.axisTable.update(objectList=axisPanels,
                          textMatrix=textMatrix)

  def setAxisState(self, axisPanel=None):

    windowPane = self.windowPane
    buttons = self.axis_buttons.buttons
    
    if windowPane:
      buttons[0].config(state='normal')
    else:
      buttons[0].config(state='disabled')

    #if windowPane and axisPanel and axisPanel.label not in ('x', 'y'):
    #  buttons[1].config(state='normal')
    #else:
    #  buttons[1].config(state='disabled')

  def getPanelTypes(self, axisPanel=None):

    if not axisPanel:
      axisPanel = self.axisPanel

    if not axisPanel:
      return []

    panelTypes = self.parent.getPanelTypes()
    panelTypes = [ panelType for panelType in panelTypes \
                   if panelType.axisType == axisPanel.axisType ]

    return panelTypes

  def selectAxis(self, axisPanel, row, col):

    self.axisPanel = axisPanel
    self.setAxisState(axisPanel)

  def getPanelType(self, axisPanel):

    panelTypes = self.getPanelTypes(axisPanel)
    entries = [panelType.name for panelType in panelTypes]
    ind = panelTypes.index(axisPanel.panelType)
    self.panelTypeWidget.setup(entries, panelTypes, ind)

  def setPanelType(self, axisPanel):

    panelType = self.panelTypeWidget.getObject()
    self.setAxisPanelType(panelType)

  def setAxisPanelType(self, panelType):

    if panelType:
      axisPanel = self.axisPanel
      axisPanel.panelType = panelType

  def getThickness(self, axisPanel):

    thickness = Util.greaterOrEqualEntry(axisPanel.slicePanel.thickness,
                                         SLICE_HEIGHT_VALS)
    self.sliceThicknessWidget.set(thickness)

  def setThickness(self, axisPanel):

    thickness = self.sliceThicknessWidget.getObject()
    self.setSliceThickness(thickness) # -1 arbitrary

  def _getTiedNames(self, axisPanels):

    return ['%s:%s' % (axisPanel.spectrumWindowPane.spectrumWindow.name, axisPanel.label) for axisPanel in axisPanels]

  def getTiedAxes(self, axisPanel):

    axisPanels = WindowBasic.getPossibleTiedAxisPanels(axisPanel)
    options = self._getTiedNames(axisPanels)
    tiedAxisPanels = WindowBasic.getTiedAxisPanels(axisPanel)
    values = [(axisPanel in tiedAxisPanels and True or False) for axisPanel in axisPanels]
    
    self.tiedAxesWidget.set(values=values, options=options)
 
  def setTiedAxes(self, obj):

    if obj is not None:
      chosen = self.tiedAxesWidget.get()
      axisPanels = WindowBasic.getPossibleTiedAxisPanels(self.axisPanel)
      tiedAxisPanels = [axisPanel for n, axisPanel in enumerate(axisPanels) if chosen[n]]
      WindowBasic.setTiedAxisPanels(self.axisPanel, tiedAxisPanels)

    self.axisTable.keyPressEscape()
    
  def setSliceThickness(self, thickness):

    thickness = Util.greaterOrEqualEntry(thickness, SLICE_HEIGHT_VALS)
    # TBD: does axisPanel.slicePanel always exists here?
    axisPanel = self.axisPanel
    if (axisPanel):
      axisPanel.slicePanel.thickness = thickness

  def toggleScrollbarVisible(self, axisPanel):

    axisPanel.isVisible = not axisPanel.isVisible
    
  def toggleSliceVisible(self, axisPanel):

    axisPanel.slicePanel.isVisible = not axisPanel.slicePanel.isVisible

  def createAxis(self):

    popup = CreateAxisPopup(self)

    axisType = popup.axisType
    popup.destroy()

    if axisType and self.windowPane:
      Util.createAxisPanel(self.windowPane, axisType)

  def deleteAxis(self, *event):

    axisPanel = self.axisPanel

    if (axisPanel and not axisPanel.isDeleted):
      self.checkDeleteAxisPanel(axisPanel)

  def checkDeleteAxisPanel(self, axisPanel):

    msg = 'Are you sure you want to delete axis "%s"?' % axisPanel.label
    if showYesNo('Delete axis', msg, parent=self):
      axisPanel.slicePanel.delete()
      axisPanel.delete()
      self.updateAxisList()
      self.setAxisState()
      #self.updateRegionList()
      #self.setRegionState()
      

  def delayedUpdateRegionList(self, *extra):

    pass
    #self.after_idle(self.updateRegionList)

  def setRegionState(self, axisRegion=None):

    windowPane = self.windowPane
    if windowPane:
      state = 'normal'
    else:
      state = 'disabled'

    for n in range(2):
      self.region_buttons.buttons[n].config(state=state)

    if axisRegion:
      axisPanel = axisRegion.axisPanel
      axisRegions = axisPanel.sortedAxisRegions()
      
      if axisPanel.label in ('x', 'y'):
        if (len(axisRegions) > 1) and (axisRegions[-1] == axisRegion):
          # TBD: for now can only delete last row/col and can never delete if only one
          state = 'normal'
        else:
          state = 'disabled'
      else:
        state = 'normal'
    else:
      state = 'disabled'

    self.region_buttons.buttons[2].config(state=state)

  def selectRegion(self, axisRegion, row, col):

    self.setRegionState(axisRegion)

  def addRow(self, *extra):

    self.windowPane.getWindowFrame().addRow()

  def addCol(self, *extra):

    self.windowPane.getWindowFrame().addCol()

  def deleteRegion(self, *extra):

    window = self.windowPane.spectrumWindow
    window_popup = self.parent.getWindowPopup(window.name)
    
    if window_popup:
      axisRegion = self.region_list.currentObject
      
      if axisRegion:
        axisPanel = axisRegion.axisPanel
        label = axisPanel.label
        
        if label == 'x':
          window_popup.deleteCol()
          
        elif label == 'y':
          window_popup.deleteRow()
          
        else: # TBD: for now assume only have one region in orthogonal regions
          self.checkDeleteAxisPanel(axisPanel)

  def setWindow(self, window):

    self.windowTable.selectObject(window)

  def changedAppData(self, appData):

    if not self.windowPane:
      return

    window = self.windowPane.spectrumWindow
    if window.root.application.name != appData.application:
      return

    if appData.keyword not in ('isCanvasLabelShown', 'isCanvasMidpointShown'):
      return

    self.updateWindowList()
  
  # Mapping functions

  def updateMapWindowPulldown(self):
  
    index = 0
    windowPanes = self.getWindowPanes()
    names = [WindowBasic.getWindowPaneName(wp) for wp in windowPanes]
    
    if windowPanes:
      if self.mapWindowPane not in windowPanes :
        self.mapWindowPane = windowPanes[0]
        
      index = windowPanes.index(self.mapWindowPane)  
      
    else:
      self.mapWindowPane = None
    
    self.mapWindowPulldown.setup(names, windowPanes, index)
  
  
  def editSpectrum(self):
  
    if self.specWindowView:
      self.parent.editSpectrum(self.specWindowView.analysisSpectrum.dataSource)
  
  def editPeakList(self):

    if self.winPeakList:
      self.parent.editPeakList(self.winPeakList.analysisPeakList.peakList)

  def togglePeakSymbol(self, winPeakList):

    winPeakList.isSymbolDrawn = not winPeakList.isSymbolDrawn

  def togglePeakAnnotation(self, winPeakList):

    winPeakList.isAnnotationDrawn = not winPeakList.isAnnotationDrawn

  def toggleInToolbar(self, view):

    view.isInToolbar = not view.isInToolbar

  def togglePosContoursVisible(self, view):

    view.isPosVisible = not view.isPosVisible
    
  def toggleNegContoursVisible(self, view):

    view.isNegVisible = not view.isNegVisible
    
  def toggleViewSliceVisible(self, view):

    view.isSliceVisible = not view.isSliceVisible

  def toggleContourLineVisible(self, view):

    view.isContourLineVisible = not view.isContourLineVisible

  def changeViewMapping(self, view):

    popup = ChangeAxisMappingPopup(self, view)
    redraw = popup.redraw
    popup.destroy()
    # slightly messy this, do not want to use notifiers
    # because have several AxisMappings deleted and re-created
    
    if redraw:
      self.parent.changedViewAxisMapping(view)
      self.updateMappingsAfter()
      window = self.mapWindowPane.spectrumWindow
      window_popup = self.parent.getWindowPopup(window.name)
      
      if window_popup:
        window_popup.drawAll()

  def getXAxisOffset(self, view):

    application = self.project.application
    value = WindowBasic.getXAxisOffset(view)
    self.xAxisOffsetWidget.set(value)

  def setXAxisOffset(self, obj):
    
    value = self.xAxisOffsetWidget.get()

    view = self.spectraMatrix.currentObject
    if view: # should always be true but play safe
    
      try:
        WindowBasic.setXAxisOffset(view, value)
      
      except Implementation.ApiError, e:
        showError('Setting x axis offset', e.error_msg, parent=self)
 
  def getValueAxisOffset(self, view):

    application = self.project.application
    value = WindowBasic.getValueAxisOffset(view)
    self.valueAxisOffsetWidget.set(value)

  def setValueAxisOffset(self, obj):
    
    value = self.valueAxisOffsetWidget.get()

    view = self.spectraMatrix.currentObject
    if view: # should always be true but play safe
    
      try:
        WindowBasic.setValueAxisOffset(view, value)
      
      except Implementation.ApiError, e:
        showError('Setting y axis offset', e.error_msg, parent=self)
 
  def selectSpectrum(self, specWindowView, row, col):
  
    self.specWindowView = specWindowView
    
    if self.specWindowView:
      self.spectrumButtons.buttons[0].enable()
    else:
      self.spectrumButtons.buttons[0].disable()

  def selectPeakList(self, winPeakList, row, col):

    self.winPeakList = winPeakList

    if self.winPeakList:
      self.peakListButtons.buttons[0].enable()
    else:
      self.peakListButtons.buttons[0].disable()

  def changeMapWindow(self, windowPane):
  
    self.mapWindowPane = windowPane
    self.updateMappingsAfter()

  def setAllSpectraVisibility(self, isVisible):
    
    windowPane = self.mapWindowPane
    if not windowPane:
      return
      
    for view in windowPane.spectrumWindowViews:
      view.isNegVisible = view.isPosVisible = view.isSliceVisible = isVisible
      for winPeakList in view.windowPeakLists:
        winPeakList.isAnnotationDrawn = winPeakList.isSymbolDrawn = isVisible
      
  def invertSpectraVisibility(self):
  
    windowPane = self.mapWindowPane
    if not windowPane:
      return
      
    for view in windowPane.spectrumWindowViews:
      view.isNegVisible = not view.isNegVisible
      view.isPosVisible = not view.isPosVisible
      view.isSliceVisible = not view.isSliceVisible
      for winPeakList in view.windowPeakLists:
        winPeakList.isAnnotationDrawn = not winPeakList.isAnnotationDrawn
        winPeakList.isSymbolDrawn = not winPeakList.isSymbolDrawn
    
  def updateMappingsAfter(self, *obj):

    if self.waitingM:
      return
    else:
      self.waitingM = True
      self.after_idle(self.updateMappings)

  def dimMappingString(self, view):

    ndim = view.analysisSpectrum.dataSource.numDim
    texts = ndim * ['']
    for axisMapping in view.axisMappings:
      dim = axisMapping.analysisDataDim.dataDim.dim
      texts[dim-1] = 'D%s=%s' % (dim, axisMapping.label)
      
    text = ','.join(texts)

    return text
  
  def updateMappings(self):

    if self.specWindowView:
      self.spectrumButtons.buttons[0].enable()
    else:
      self.spectrumButtons.buttons[0].disable()
    
    if self.winPeakList:
      self.peakListButtons.buttons[0].enable()
    else:
      self.peakListButtons.buttons[0].disable()


    if self.mapWindowPane:
      views = self.mapWindowPane.spectrumWindowViews
      hasValueAxis = WindowBasic.windowPaneHasValueAxis(self.mapWindowPane)
   
    else:
      views = []
      hasValueAxis = False  # arbitrary

    tipTexts1 = ['The name of the experiment record for the window-mappable spectrum',
                'The name of the window-mappable spectrum',
                'Whether the spectrum appears with its own toggle button in the toolbar at the top of the spectrum window']
 
    tipTexts2 = ['Whether the positive contour levels are shown in the window for the spectrum',
                'Whether the negative contour levels are shown in the window for the spectrum']

    tipTexts3 = ['Whether 1D slices are visible for the spectrum, if any slice displays are active for the window']

    tipTexts4 = ['Whether contour line in 1D window is visible for the spectrum; if on this helps indicate where peaks would be picked']

    tipTexts5 = ['Sets how the spectrum dimensions map to corresponding window axes; for a given kind, enables swapping of spectrum dimensions relative to the window']
    tipTexts6 = ['How much the spectrum is offset in the x axis', 'How much the spectrum is offset in the y axis']
    headings1 = ('Experiment', 'Spectrum', 'Spectrum In\nToolbar')
    headings2 = ('Pos. Contours\nVisible', 'Neg. Contours\nVisible')
    headings3 = ('Slice\nVisible',)
    headings4 = ('Contour Line\nVisible',)
    headings5 = ('Dim. Mapping',)
    headings6 = ('X Offset', 'Y Offset',)
    editGetCallbacks1 = [None, None, self.toggleInToolbar]
    editGetCallbacks2 = [self.togglePosContoursVisible, self.toggleNegContoursVisible]
    editGetCallbacks3 = [self.toggleViewSliceVisible]
    editGetCallbacks4 = [self.toggleContourLineVisible]
    editGetCallbacks5 = [self.changeViewMapping]
    editGetCallbacks6 = [self.getXAxisOffset, self.getValueAxisOffset]

    if hasValueAxis:
      tipTexts = tipTexts1 + tipTexts3 + tipTexts4 + tipTexts5 + tipTexts6
      headings = headings1 + headings3 + headings4 + headings5 + headings6
      editGetCallbacks = editGetCallbacks1 + editGetCallbacks3 + editGetCallbacks4 + editGetCallbacks5 + editGetCallbacks6
      n = len(headings) - 2
      editSetCallbacks = n * [None] + [self.setXAxisOffset, self.setValueAxisOffset]
      editWidgets = n * [None] + [self.xAxisOffsetWidget, self.valueAxisOffsetWidget]
    else:
      tipTexts = tipTexts1 + tipTexts2 + tipTexts3 + tipTexts5
      headings = headings1 + headings2 + headings3 + headings5
      editGetCallbacks = editGetCallbacks1 + editGetCallbacks2 + editGetCallbacks3 + editGetCallbacks5
      n = len(headings)
      editWidgets = editSetCallbacks = n * [None]
    
    textMatrix = []
    objectList = []
    for view in views:
      try:
        spectrum = view.analysisSpectrum.dataSource
      
        text = [spectrum.experiment.name,
                spectrum.name,
                view.isInToolbar and 'Yes' or 'No']
                
        if not hasValueAxis:
          text.append(view.isPosVisible and 'Yes' or 'No')
          text.append(view.isNegVisible and 'Yes' or 'No')
          
        text.append(view.isSliceVisible and 'Yes' or 'No')

        if hasValueAxis:
          text.append(view.isContourLineVisible and 'Yes' or 'No')

        text.append(self.dimMappingString(view))
        
        if hasValueAxis:
          text.append(WindowBasic.getXAxisOffset(view))
          text.append(WindowBasic.getValueAxisOffset(view))
          
        textMatrix.append(text)
        objectList.append(view)
      
      except:
        pass

    self.spectraMatrix.update(objectList=objectList,
                              textMatrix=textMatrix,
                              headingList = headings,
                              tipTexts=tipTexts,
                              editWidgets=editWidgets,
                              editGetCallbacks=editGetCallbacks,
                              editSetCallbacks=editSetCallbacks)
    
    textMatrix = []
    objectList = []
    for view in views:
      winPeakLists = view.windowPeakLists
      for winPeakList in winPeakLists:
        try:
          peakList = winPeakList.analysisPeakList.peakList
          text = [view.analysisSpectrum.dataSource.experiment.name,
                  view.analysisSpectrum.dataSource.name,
                  peakList.serial,
                  winPeakList.isSymbolDrawn and 'Yes' or 'No',
                  winPeakList.isAnnotationDrawn and 'Yes' or 'No']
 
          textMatrix.append(text)
          objectList.append(winPeakList)
        except:
          pass

    self.peakListMatrix.update(objectList=objectList,
                              textMatrix=textMatrix)
    
    self.waitingM = False
    
  # Window group functions

  def updateGroupList(self, notifyObj=None):

    if notifyObj and (self.which_windows.get() == ACTIVE_OPTIONS[1]):
      self.updateWindowList()

    activeGroup = self.analysisProject.activeWindowGroup
    groups = self.analysisProject.sortedSpectrumWindowGroups()
    
    textMatrix = []
    for i, group in enumerate(groups):

      name = group.name
      if group is activeGroup:
        name = name + ' (active)'

      text = [i+1, name, len(group.spectrumWindows)]
      textMatrix.append(text)

    self.windowGroupTable.update(objectList=groups,
                                 textMatrix=textMatrix)

  def selectGroup(self, group, row, col):

    self.windowGroup = group
    self.setGroupState(group)
    self.updateGroupWindows()

  def getGroupName(self, group):

    self.windowGroupNameWidget.set(group.name)

  def setGroupName(self, *extra):

    name = self.windowGroupNameWidget.get()
    group = self.windowGroup
    
    try:
      group.name = name
      
    except Implementation.ApiError, e:
      showError('Changing group name', e.error_msg, parent=self)

  def setGroupState(self, group=None):
 
    if group and (self.analysisProject.activeWindowGroup != group):
      state = 'normal'
    else:
      state = 'disabled'
 
    self.windowGroup_buttons.buttons[0].config(state=state)

    if group:
      state = 'normal'
    else:
      state = 'disabled'
 
    self.windowGroup_buttons.buttons[2].config(state=state)

  def openGroup(self):

    group = self.windowGroup
    if group:
    
      self.parent.openWindowGroup(group)
      self.setGroupState(group)

  def createGroup(self):

    #msg = 'What is name of new window group?'
    #name = askString('Window group name', msg, parent=self)
    #if name:
    #project = self.project
    
    try:
      group = self.analysisProject.newSpectrumWindowGroup(name='New Group')
      group.name = 'New Group %s' % group.serial
      self.windowGroupTable.selectObject(group)
      
    except Implementation.ApiError, e:
      showError('Window group creation error', e.error_msg, parent=self)

  def deleteGroup(self, *event):

    groups = self.windowGroupTable.currentObjects

    if len(groups) == 1:
      msg = 'Are you sure you want to delete group "%s"?'
      
      if showYesNo('Delete group', msg  % groups[0].name, parent=self):
        groups[0].delete()
        
    elif len(groups) > 1:
      names = ' '.join([g.name for g in groups])
      msg =  'Are you sure you want to delete %d groups [%s]?'
      
      if showYesNo('Delete group', msg % (len(groups),names), parent=self):
        for g in groups:
          g.delete()


  def updateGroupWindowsAfter(self, *event):

    self.after_idle(self.updateGroupWindows) # Give axes a moment to be made


  def updateGroupWindows(self, *extra):

    windows = self.analysisProject.sortedSpectrumWindows()
    group = self.windowGroup

    textMatrix = []
    for i, window in enumerate(windows):
      
      if not group:
        inGroup = None
      elif window in group.spectrumWindows:
        inGroup = 'Yes'
      else:
        inGroup = 'No'
        
      text = [i+1, window.name, inGroup]
      
      textMatrix.append(text)

    self.groupWindowsTable.update(objectList=windows,
                                  textMatrix=textMatrix)


  def getInGroup(self, window):

    group = self.windowGroup
    if group:
    
      if window in group.spectrumWindows:
        group.removeSpectrumWindow(window)
        
      else:
        group.addSpectrumWindow(window)


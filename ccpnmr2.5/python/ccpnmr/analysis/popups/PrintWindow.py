
"""
======================COPYRIGHT/LICENSE START==========================

PrintWindow.py: Part of the CcpNmr Analysis program

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
from memops.gui.Button import Button
from memops.gui.ButtonList import UtilityButtonList, ButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.FloatEntry import FloatEntry
from memops.gui.FontList import FontList, printNames
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.LabelDivider import LabelDivider
from memops.gui.PrintFrame import PrintFrame, spacing_choices, tick_length_choices
from memops.gui.PulldownList import PulldownList
from memops.gui.MessageReporter import showError, showInfo
from memops.gui.RadioButtons import RadioButtons
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.WindowBasic import getWindowPaneName, getAxisRegionRegion, windowPaneHasValueAxis
from ccpnmr.analysis.core.WindowDraw import WindowDraw, no_peak_text

import ccpnmr.analysis.core.PrintBasic as PrintBasic

USE_ENTRIES = ('Use Min and Max', 'Use Center and Width')

MIN_MAX_HEADINGS = [ 'Axis', 'Region\nNumber',
                 'Override Region\nMin (ppm)',
                 'Override Region\nMax (ppm)' ]

CENTER_WIDTH_HEADINGS = [ 'Axis', 'Region\nNumber',
                 'Override Region\nCenter (ppm)',
                 'Override Region\nWidth (ppm)' ]

class RegionScrolledMatrix(ScrolledMatrix):

  def doEditMarkExtraRules(self, axisRegion, row, col):

    window = axisRegion.axisPanel.spectrumWindowPane.spectrumWindow
    if window.useOverrideRegion:
      return True
    else:
      return False

class PrintWindowPopup(BasePopup):

  """
  **Print Window to Output File**

  The purpose of this dialog is to allow the saving of the drawing of
  one of the spectrum windows to a file, in one of the following formats:
  PostScript (PS), Encapsulated PostScript (EPS) or Portable Document
  Format (PDF).

  The one window that is being printed out is specified at the top.
  There are four tabs.  The first one, Options, is the most important.
  In particular, it is used to specify the File name.  At its simplest
  to print out a window you just need to specify the File name, and
  then click "Save Print File".  But it is likely you will at the very
  least want to change some of the settings in the Options tab.

  You can specify a Title and a label for the X axis and/or Y axis.

  This tab is also used to specify the Paper size (default A4), the
  Orientation of the paper (default Portrait), whether the printout
  is Color or Black and white (the Style, default Color), and what
  the Format is (PS, EPS or PDF, default PS).

  The ticks for the rulers can be chosen to be Inside or Outside the
  main frame and can be in any combination of the Top, Bottom, Left
  or Right side of the main frame.  The Tick Font includes the option
  of not printing the tick labels at all.  The Tick spacing between
  the major and minor ticks can be set automatically (so the program
  determines it) or manually.  For the latter the user has to specify
  the Major and Minor spacings in terms of the unit of the display
  (normally ppm), and also the number of decimal places for the Tick
  labels.

  How the main frame fits into the paper is determined by the Scaling
  option.  The Percentage option just means that the main frame is
  scaled by that amount relative to the biggest size it could be and
  still fit on the paper.  The remaining options are if you want to
  specify the cms or inches per unit (normally ppm), or the inverse
  of these.  In this case it could be the case that the main frame
  actually exceeds the size of the paper.

  You can also include the Time and Date and/or the File Name in the
  printout, and you can specify the font used for this.  The same
  font is used for the Title and X and Y axis labels, except that
  the Title font is 6 pts bigger.

  Finally, you can also set the linewidth, in points (the default is 0.1).

  The other three tabs provide fine tuning of what is output.  In many
  cases they can be ignored.

  The Spectra tab lets you choose settings for which of the window's
  spectra are drawn, in terms of both the positive and negative contours.
  This is independent of what is actually drawn on the screen.  But you
  need to check the "Use below settings when printing" checkbutton if you
  want the values specified here to be used rather than the screen settings.
  The values in the table are initially set to be the screen values but
  afterwards can only be changed manually.  Clicking on the "Reset Selected"
  button changes the values in the table to the current screen ones.

  The Peak Lists tab is similar, except it applies to the peak lists in
  the window rather than the spectra.  Again, if you want the values in
  the table to be used then you need to check the "Use below settings
  when printing" checkbutton.

  The Region tab is for specifying the region for the x, y and orthogonal
  axes, if you do not want to use the regions as seen in the window on the
  screen.  Again, you have to check "Use override region when printing"
  checkbutton to actually have the values in the table be used in the
  printout.  By default, the override region is set to the current
  window region if it is not set already, otherwise it is left to the
  previous value unless you click the "Set Region from Window" or the
  "Set Width from Window" or the "Set Center from Window" buttons.
  And the override region can be specified either using the min and max
  values of the region, or the center and width.
"""

  def __init__(self, parent, *args, **kw):

    self.waiting = False
    
    title='Window : Print Window'
    BasePopup.__init__(self, parent=parent, title=title, **kw)

  def body(self, guiFrame):

    self.geometry('600x350')

    project = self.project
    analysisProject = self.analysisProject

    guiFrame.grid_columnconfigure(1, weight=1)

    row = 0
    label = Label(guiFrame, text=' Window:', grid=(0,0))
    self.windowPulldown = PulldownList(guiFrame, grid=(0,1),
                                      tipText='The window that will be printed out',
                                      callback=self.selectWindow)
    texts = [ 'Save Print File' ]
    tipTexts = [ 'Save the printout to the specified file' ]
    commands = [ self.saveFile ]
    buttons = UtilityButtonList(guiFrame, helpUrl=self.help_url, grid=(row,2),
                                commands=commands, texts=texts, tipTexts=tipTexts)
    self.buttons = buttons
    buttons.buttons[0].config(bg='#B0FFB0')
    
    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    options = ['Options', 'Spectra', 'Peak Lists', 'Region']
    tipTexts = ['Optional settings for spectra', 'Optional settings for peak lists', 'Optional settings for the region']
    tabbedFrame = TabbedFrame(guiFrame, options=options, tipTexts=tipTexts)
    tabbedFrame.grid(row=row, column=0, columnspan=3, sticky='nsew')
    self.tabbedFrame = tabbedFrame

    optionFrame, spectrumFrame, peakListFrame, regionFrame = tabbedFrame.frames

    optionFrame.expandGrid(0, 0)
    getOption = lambda key, defaultValue: PrintBasic.getPrintOption(analysisProject, key, defaultValue)
    setOption = lambda key, value: PrintBasic.setPrintOption(analysisProject, key, value)
    self.printFrame = PrintFrame(optionFrame, getOption=getOption,
                                 grid=(0,0), gridSpan=(1,1),
                                 setOption=setOption, haveTicks=True,
                                 doOutlineBox=False)

    spectrumFrame.expandGrid(0, 0)
    frame = Frame(spectrumFrame, grid=(0,0), gridSpan=(1,1))
    frame.expandGrid(1,0)

    self.overrideSpectrum = CheckButton(frame,
       text='Use below settings when printing',
       tipText='Use below settings when printing instead of the window values',
       grid=(0,0), sticky='w')

    tipText = 'Change the settings of the selected spectra back to their window values'
    button = Button(frame, text='Reset Selected', tipText=tipText,
                    command=self.resetSelected, grid=(0,1), sticky='e')

    headings = ['Spectrum', 'Pos. Contours\nDrawn', 'Neg. Contours\nDrawn']
    tipTexts = ['Spectrum in window', 'Whether the positive contours should be drawn', 'Whether the negative contours should be drawn']
    editWidgets      = [ None, None, None]
    editGetCallbacks = [ None, self.togglePos, self.toggleNeg]
    editSetCallbacks = [ None, None, None]
    self.spectrumTable = ScrolledMatrix(frame, headingList=headings,
                                        tipTexts=tipTexts,
                                        multiSelect=True,
                                        editWidgets=editWidgets,
                                        editGetCallbacks=editGetCallbacks,
                                        editSetCallbacks=editSetCallbacks,
                                        grid=(1,0), gridSpan=(1,2))

    peakListFrame.expandGrid(0, 0)
    frame = Frame(peakListFrame, grid=(0,0), gridSpan=(1,3))
    frame.expandGrid(1,0)

    self.overridePeakList = CheckButton(frame,
       text='Use below settings when printing',
       tipText='Use below settings when printing instead of the window values',
       grid=(0,0))

    tipText = 'Change the settings of the selected peak lists back to their window values'
    button = Button(frame, text='Reset Selected', tipText=tipText,
                    command=self.resetSelected, grid=(0,1), sticky='e')

    headings = [ 'Peak List', 'Symbols Drawn', 'Peak Font']
    self.fontMenu = FontList(self, mode='Print', extraTexts=[no_peak_text])
    editWidgets      = [ None, None, self.fontMenu]
    editGetCallbacks = [ None, self.togglePeaks, self.getPeakFont ]
    editSetCallbacks = [ None, None, self.setPeakFont ]
    self.peakListTable = ScrolledMatrix(frame, headingList=headings,
                                        multiSelect=True,
                                        editWidgets=editWidgets,
                                        editGetCallbacks=editGetCallbacks,
                                        editSetCallbacks=editSetCallbacks,
                                        grid=(1,0), gridSpan=(1,2))

    regionFrame.expandGrid(0, 0)
    frame = Frame(regionFrame, grid=(0,0), gridSpan=(1,3))
    frame.expandGrid(3,0)
    tipText = 'Use the specified override region when printing rather than the window values'
    self.overrideButton = CheckButton(frame, text='Use override region when printing',
                                      tipText=tipText,
                                      callback=self.toggledOverride, grid=(0,0))

    tipTexts = ('Use min and max to specify override region', 'Use center and width to specify override region')
    self.use_entry = USE_ENTRIES[0]
    self.useButtons = RadioButtons(frame, entries=USE_ENTRIES,
                                      tipTexts=tipTexts,
                                      select_callback=self.changedUseEntry,
                                      grid=(1,0))

    texts = ('Set Region from Window', 'Set Center from Window', 'Set Width from Window')
    tipTexts = ('Set the override region to be the current window region',
                'Set the center of the override region to be the center of the current window region',
                'Set the width of the override region to be the width of the current window region')
    commands = (self.setRegionFromWindow, self.setCenterFromWindow, self.setWidthFromWindow)
    self.setRegionButton = ButtonList(frame, texts=texts,
                                      tipTexts=tipTexts,
                                      commands=commands, grid=(2,0))

    self.minRegionWidget = FloatEntry(self, returnCallback=self.setMinRegion, width=10)
    self.maxRegionWidget = FloatEntry(self, returnCallback=self.setMaxRegion, width=10)
    headings = MIN_MAX_HEADINGS
    editWidgets      = [ None, None, self.minRegionWidget, self.maxRegionWidget ]
    editGetCallbacks = [ None, None, self.getMinRegion,    self.getMaxRegion ]
    editSetCallbacks = [ None, None, self.setMinRegion,    self.setMaxRegion ]
    self.regionTable = RegionScrolledMatrix(frame, headingList=headings,
                                            editWidgets=editWidgets,
                                            editGetCallbacks=editGetCallbacks,
                                            editSetCallbacks=editSetCallbacks,
                                            grid=(3,0))

    self.updateWindows()
    self.updateAfter()

    self.administerNotifiers(self.registerNotify)

  def open(self):
    # NOTE:ED - subclass the open so that the window pulldown gets updated on subsequent open events
    super(PrintWindowPopup, self).open()
    self.updateWindows()

  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete', 'setOverrideRegion'):
      notifyFunc(self.updateAfter, 'ccpnmr.Analysis.AxisRegion', func)

    notifyFunc(self.updateAfter, 'ccpnmr.Analysis.SpectrumWindow', 'setUseOverrideRegion')

    for func in ('__init__', 'delete', 'setName'):
      notifyFunc(self.updateWindows, 'ccpnmr.Analysis.SpectrumWindow', func)
      notifyFunc(self.updateWindows, 'ccpnmr.Analysis.SpectrumWindowPane', func)

  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

  def changedUseEntry(self, entry):

    self.use_entry = entry
    self.updateRegionTable()

  def resetSpectrum(self):

    spectrumWindowViews = self.spectrumTable.currentObjects
    for spectrumWindowView in spectrumWindowViews:
      PrintBasic.setPrintOption(spectrumWindowView, 'PositiveOn', spectrumWindowView.isPosVisible)
      PrintBasic.setPrintOption(spectrumWindowView, 'NegativeOn', spectrumWindowView.isNegVisible)

    self.updateAfter()

  def resetPeakList(self):

    windowPeakLists = self.peakListTable.currentObjects
    for windowPeakList in windowPeakLists:
      PrintBasic.setPrintOption(windowPeakList, 'PeaksOn', windowPeakList.isSymbolDrawn)
      PrintBasic.setPrintOption(windowPeakList, 'PeakFont', windowPeakList.spectrumWindowView.analysisSpectrum.font)

    self.updateAfter()

  def resetSelected(self):

    n = self.tabbedFrame.selected
    if n == 0:
      self.resetSpectrum()
    elif n == 1:
      self.resetPeakList()

  def togglePos(self, spectrumWindowView):

    PrintBasic.setPrintOption(spectrumWindowView, 'PositiveOn',
      not PrintBasic.getPrintOption(spectrumWindowView, 'PositiveOn', defaultValue=spectrumWindowView.isPosVisible))

    self.updateAfter()

  def toggleNeg(self, spectrumWindowView):

    PrintBasic.setPrintOption(spectrumWindowView, 'NegativeOn',
      not PrintBasic.getPrintOption(spectrumWindowView, 'NegativeOn', defaultValue=spectrumWindowView.isNegVisible))

    self.updateAfter()

  def togglePeaks(self, windowPeakList):

    PrintBasic.setPrintOption(windowPeakList, 'PeaksOn',
      not PrintBasic.getPrintOption(windowPeakList, 'PeaksOn', defaultValue=windowPeakList.isSymbolDrawn))

    self.updateAfter()

  def getPeakFont(self, windowPeakList):

    if windowPeakList.isAnnotationDrawn:
      default = windowPeakList.analysisPeakList.analysisSpectrum.font
    else:
      default = no_peak_text
    font = PrintBasic.getPrintOption(windowPeakList, 'PeakFont', defaultValue=default)

    self.fontMenu.set(font)

  def setPeakFont(self, windowPeakList):

    PrintBasic.setPrintOption(windowPeakList, 'PeakFont', self.fontMenu.getText())

    self.updateAfter()

  def updateWindows(self, obj=None):

    index = 0
    windowPane = self.windowPulldown.getObject()
    windowPanes = []
    names = []
    
    if not windowPane:
      application = self.project.application
      name = application.getValue(self.analysisProject, keyword='printWindowWindow')
      window = self.analysisProject.findFirstSpectrumWindow(name=name)
      if window:
        windowPane = window.findFirstSpectrumWindowPane()
     

    for window in self.analysisProject.sortedSpectrumWindows():
      for windowPane0 in window.sortedSpectrumWindowPanes():
        windowPanes.append(windowPane0)
        names.append(getWindowPaneName(windowPane0))
    
    if windowPanes:
      if windowPane not in windowPanes:
        windowPane = windowPanes[0]
      
      index = windowPanes.index(windowPane)  
    
    else:
      windowPane = None
      
    self.selectWindow(windowPane)
        
    self.windowPulldown.setup(names, windowPanes, index)

  def saveFile(self):

    windowPane = self.windowPulldown.getObject()

    if not windowPane:
      return

    axisPanels = windowPane.sortedAxisPanels()
    aspectRatio = self.getPrintAspectRatio()
    pixelWidth = self.totalSize(axisPanels[0])
    pixelHeight = aspectRatio*self.totalSize(axisPanels[1])
    unitWidth = self.totalOverrideRegion(axisPanels[0])
    unitHeight = self.totalOverrideRegion(axisPanels[1])

    printFrame = self.printFrame

    isOverrideSpectrumSelected = self.overrideSpectrum.isSelected()
    if isOverrideSpectrumSelected:
      spectrumWindowViews = self.spectrumTable.objectList
      # alternatively, spectrumWindowViews = windowPane.spectrumWindowViews
      for spectrumWindowView in spectrumWindowViews:
        spectrumWindowView.printPositive = PrintBasic.getPrintOption(spectrumWindowView, 'PositiveOn', spectrumWindowView.isPosVisible)
        spectrumWindowView.printNegative = PrintBasic.getPrintOption(spectrumWindowView, 'NegativeOn', spectrumWindowView.isNegVisible)

    isOverridePeakListSelected = self.overridePeakList.isSelected()
    if isOverridePeakListSelected:
      windowPeakLists = self.peakListTable.objectList
      for windowPeakList in windowPeakLists:
        windowPeakList.printPeaks = PrintBasic.getPrintOption(windowPeakList, 'PeaksOn', windowPeakList.isSymbolDrawn)
        if windowPeakList.isAnnotationDrawn:
          default = windowPeakList.analysisPeakList.analysisSpectrum.font
        else:
          default = no_peak_text
        windowPeakList.printFont = PrintBasic.getPrintOption(windowPeakList, 'PeakFont', default)

    xrr = axisPanels[0].findFirstAxisRegion().region
    dxx = abs(xrr[0]-xrr[1])
    yrr = axisPanels[1].findFirstAxisRegion().region
    dyy = abs(yrr[0]-yrr[1])
    try:
      outputHandler = printFrame.getOutputHandler(pixelWidth, pixelHeight, unitWidth, unitHeight, fonts=printNames)
      if not outputHandler:
        return
      
      analysisProject = self.analysisProject
      major_minor_dict = {}
      spacing_choice = PrintBasic.getPrintOption(analysisProject, 'SpacingChoice', spacing_choices[0])
      if spacing_choice != spacing_choices[0]:
        for attr in ('XMajor', 'XMinor', 'XDecimal',
                     'YMajor', 'YMinor', 'YDecimal',):
          val = PrintBasic.getPrintOption(analysisProject, attr, None)
          if val is not None:
            major_minor_dict[attr] = val

      tick_length_choice = PrintBasic.getPrintOption(analysisProject, 'TickLengthChoice', tick_length_choices[0])
      if tick_length_choice != tick_length_choices[0]:
        for attr in ('TickMajor', 'TickMinor'):
          val = PrintBasic.getPrintOption(analysisProject, attr, None)
          if val is not None:
            major_minor_dict[attr] = val

      windowDraw = WindowDraw(self.parent, windowPane)

      PrintBasic.printWindow(windowDraw,
                             outputHandler,
                             printFrame.tick_location,
                             printFrame.tick_placement,
                             aspectRatio,
                             printFrame.tick_font,
                             major_minor_dict)
      
      msg = 'Saved to file "%s"' % printFrame.file_name                        
      showInfo('Success', msg, parent=self)
      
    except IOError, e:
      showError('IO Error', str(e), parent=self)

    if isOverrideSpectrumSelected:
      for spectrumWindowView in spectrumWindowViews:
        del spectrumWindowView.printPositive
        del spectrumWindowView.printNegative

    if isOverridePeakListSelected:
      for windowPeakList in windowPeakLists:
        del windowPeakList.printPeaks
        del windowPeakList.printFont

  def getPrintAspectRatio(self):

    windowPane = self.windowPulldown.getObject()
    window = windowPane.spectrumWindow
    
    if windowPaneHasValueAxis(windowPane):
      aspectRatio = self.printFrame.getAspectRatio()
      
    elif window.useOverrideRegion:
      axisPanel = windowPane.findFirstAxisPanel(label='x')
      x = self.totalOverrideRegion(axisPanel) / float(self.totalRegion(axisPanel))
      axisPanel = windowPane.findFirstAxisPanel(label='y')
      y = self.totalOverrideRegion(axisPanel) / float(self.totalRegion(axisPanel))
      ####aspectRatio = y / (window.aspectRatio*x)
      aspectRatio = y / x
    
    else:
      aspectRatio = 1

    return aspectRatio

  def totalRegion(self, axisPanel):

    r = 0
    for axisRegion in axisPanel.axisRegions:
      (r0, r1) = axisRegion.region
      r += r1 - r0

    return abs(r)

  def totalOverrideRegion(self, axisPanel):

    r = 0
    for axisRegion in axisPanel.axisRegions:
      if axisRegion.overrideRegion:
        (r0, r1) = axisRegion.overrideRegion
      else:
        (r0, r1) = axisRegion.region
      r = r + r1 - r0

    return abs(r)

  def totalSize(self, axisPanel):

    size = 0
    for axisRegion in axisPanel.axisRegions:
      size = size + axisRegion.size

    return size

  def setWindow(self, windowPane):

    self.windowPulldown.setSelected(windowPane)

  def updateAfter(self, *extra):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)

  def update(self):

    windowPane = self.windowPulldown.getObject()
    
    if windowPane:
      self.overrideButton.set(windowPane.spectrumWindow.useOverrideRegion)
      self.updateSpectrumTable()
      self.updatePeakListTable()
      self.updateRegionTable()
    
    self.waiting = False

  def selectWindow(self, windowPane):

    if windowPane:
      # TBD: change this code when printWindowName in data model
      
      application = self.project.application
      window = windowPane.spectrumWindow
      self.printFrame.setOptionValues()
      application.setValue(self.analysisProject,
                           keyword='printWindowWindow',
                           value=window.name)
                           
      self.updateSpectrumTable()
      self.updatePeakListTable()
      self.updateRegionTable()

  def updateSpectrumTable(self):

    textMatrix = []
    spectrumWindowViews = []

    windowPane = self.windowPulldown.getObject()
    if windowPane:
      # Note that windowPane.sortedSpectrumWindowViews() does not sort
      # consistently from one instantation to the next so do it by hand
      spectrumWindowViews = windowPane.spectrumWindowViews
      tt = [(x.analysisSpectrum.dataSource, x) for x in spectrumWindowViews]
      tt.sort()
      spectrumWindowViews = [x[1] for x in tt]
      for spectrumWindowView in spectrumWindowViews:
        analysisSpectrum = spectrumWindowView.analysisSpectrum
        spectrum = analysisSpectrum.dataSource
        experiment = spectrum.experiment
        
        text = []
        text.append('%s:%s' % (experiment.name, spectrum.name))
        positiveOn = PrintBasic.getPrintOption(spectrumWindowView, 'PositiveOn', spectrumWindowView.isPosVisible)
        text.append(positiveOn and 'Yes' or 'No')
        negativeOn = PrintBasic.getPrintOption(spectrumWindowView, 'NegativeOn', spectrumWindowView.isNegVisible)
        text.append(negativeOn and 'Yes' or 'No')
        textMatrix.append(text)

    self.spectrumTable.update(objectList=spectrumWindowViews,
                            textMatrix=textMatrix)

  def updatePeakListTable(self):

    textMatrix = []
    windowPeakLists = []

    windowPane = self.windowPulldown.getObject()
    if windowPane:
      # Note that windowPane.sortedSpectrumWindowViews() does not sort
      # consistently from one instantation to the next nor does
      # spectrumWindowView.sortedWindowPeakLists() so do it by hand
      for spectrumWindowView in windowPane.spectrumWindowViews:
        windowPeakLists.extend(spectrumWindowView.windowPeakLists)
      tt = [(x.analysisPeakList.peakList, x) for x in windowPeakLists]
      tt.sort()
      windowPeakLists = [x[1] for x in tt]

      for windowPeakList in windowPeakLists:
        analysisPeakList = windowPeakList.analysisPeakList
        peakList = analysisPeakList.peakList
        spectrum = peakList.dataSource
        experiment = spectrum.experiment
        
        text = []
        text.append('%s:%s:%s' % (experiment.name, spectrum.name, peakList.serial))
        peaksOn = PrintBasic.getPrintOption(windowPeakList, 'PeaksOn', windowPeakList.isSymbolDrawn)
        text.append(peaksOn and 'Yes' or 'No')
        if windowPeakList.isAnnotationDrawn:
          default = analysisPeakList.analysisSpectrum.font
        else:
          default = no_peak_text
        text.append(PrintBasic.getPrintOption(windowPeakList, 'PeakFont', default))
        textMatrix.append(text)

    self.peakListTable.update(objectList=windowPeakLists,
                              textMatrix=textMatrix)

  def updateRegionTable(self):

    if self.use_entry == USE_ENTRIES[0]:
      headings = MIN_MAX_HEADINGS
      isCenterWidth = False
    else:
      headings = CENTER_WIDTH_HEADINGS
      isCenterWidth = True

    textMatrix = []
    axisRegions = []
    
    windowPane = self.windowPulldown.getObject()
    if windowPane:
      window = windowPane.spectrumWindow
      
      for axisPanel in windowPane.sortedAxisPanels():
        if axisPanel.axisType.isSampled:
          continue
 
        n = 1
        for axisRegion in axisPanel.sortedAxisRegions():
          if axisRegion.overrideRegion:
            (r0, r1) = axisRegion.overrideRegion
            if isCenterWidth:
              (r0, r1) = (0.5*(r0+r1), r1-r0)
 
          elif window.useOverrideRegion:
            (r0, r1) = axisRegion.region
            if isCenterWidth:
              (r0, r1) = (0.5*(r0+r1), r1-r0)
 
          else:
            r0 = r1 = ''
 
          text = [ axisPanel.label, n, r0, r1 ]
          textMatrix.append(text)
          axisRegions.append(axisRegion)
          n += 1

    self.regionTable.update(headingList=headings,
                            objectList=axisRegions,
                            textMatrix=textMatrix)

  def getMinRegion(self, axisRegion):

    if axisRegion.overrideRegion:
      if self.use_entry == USE_ENTRIES[0]:
        r = axisRegion.overrideRegion[0]
      else:
        r = 0.5 * (axisRegion.overrideRegion[0] + axisRegion.overrideRegion[1])
    else:
      r = ''
    self.minRegionWidget.set(r)

  def setMinRegion(self, *extra):

    axisRegion = self.getAxisRegion()
    if not axisRegion:
      return

    r = self.minRegionWidget.get()
    if axisRegion.overrideRegion:
      (r0, r1) = axisRegion.overrideRegion
    else:
      (r0, r1) = axisRegion.region

    if self.use_entry == USE_ENTRIES[0]:
      if r is not None and r < r1:
        axisRegion.overrideRegion = (r, r1)
    else:
      if r is not None:
        center = r
        halfwidth = 0.5 * (r1-r0)
        axisRegion.overrideRegion = (center-halfwidth, center+halfwidth)

  def getMaxRegion(self, axisRegion):

    if axisRegion.overrideRegion:
      if self.use_entry == USE_ENTRIES[0]:
        r = axisRegion.overrideRegion[1]
      else:
        r = axisRegion.overrideRegion[1] - axisRegion.overrideRegion[0]
    else:
      r = ''
    self.maxRegionWidget.set(r)

  def setMaxRegion(self, *extra):

    axisRegion = self.getAxisRegion()
    if not axisRegion:
      return

    r = self.maxRegionWidget.get()
    if axisRegion.overrideRegion:
      (r0, r1) = axisRegion.overrideRegion
    else:
      (r0, r1) = axisRegion.region

    if self.use_entry == USE_ENTRIES[0]:
      if r is not None and r > r0:
        axisRegion.overrideRegion = (r0, r)
    else:
      if r is not None and r > 0:
        halfwidth = 0.5 * r
        center = 0.5 * (r0+r1)
        axisRegion.overrideRegion = (center-halfwidth, center+halfwidth)

  def getAxisRegion(self):

    return self.regionTable.currentObject

  def toggledOverride(self, isSelected):

    windowPane = self.windowPulldown.getObject()

    if windowPane:
      window = windowPane.spectrumWindow
      window.useOverrideRegion = isSelected
      self.updateAfter()

  def setRegionFromWindow(self):

    windowPane = self.windowPulldown.getObject()
    if windowPane:
      window = windowPane.spectrumWindow
      
      for axisPanel in windowPane.sortedAxisPanels():
        if axisPanel.axisType.isSampled:
          continue
 
        for axisRegion in axisPanel.sortedAxisRegions():
          axisRegion.overrideRegion = axisRegion.region

      self.updateRegionTable()

  def setCenterFromWindow(self):

    windowPane = self.windowPulldown.getObject()
    if windowPane:
      window = windowPane.spectrumWindow
      
      for axisPanel in windowPane.sortedAxisPanels():
        if axisPanel.axisType.isSampled:
          continue
 
        for axisRegion in axisPanel.sortedAxisRegions():
          overrideRegion = axisRegion.overrideRegion
          if overrideRegion:
            (r0, r1) = overrideRegion
            halfwidth = 0.5*(r1 - r0)
            (r0, r1) = axisRegion.region
            center = 0.5*(r0 + r1)
            axisRegion.overrideRegion = (center-halfwidth, center+halfwidth)
          else:
            axisRegion.overrideRegion = axisRegion.region

      self.updateRegionTable()

  def setWidthFromWindow(self):

    windowPane = self.windowPulldown.getObject()
    if windowPane:
      window = windowPane.spectrumWindow
      
      for axisPanel in windowPane.sortedAxisPanels():
        if axisPanel.axisType.isSampled:
          continue
 
        for axisRegion in axisPanel.sortedAxisRegions():
          overrideRegion = axisRegion.overrideRegion
          if overrideRegion:
            (r0, r1) = overrideRegion
            center = 0.5*(r0 + r1)
            (r0, r1) = axisRegion.region
            halfwidth = 0.5*(r1 - r0)
            axisRegion.overrideRegion = (center-halfwidth, center+halfwidth)
          else:
            axisRegion.overrideRegion = axisRegion.region

      self.updateRegionTable()


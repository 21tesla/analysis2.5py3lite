
LICENSE = """
======================COPYRIGHT/LICENSE START==========================

PrintBasic.py: Part of the CcpNmr Analysis program

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
from memops.universal.PrintHandler import PrintHandler
from memops.universal.PrintTicks import PrintTicks, Top, Bottom, Left, Right

from memops.general.Implementation import ApiError

from ccpnmr.analysis.core.WindowBasic import getAxisRegionRegion, windowPaneHasValueAxis

# windowDraw supposed to be object of type WindowDraw
# but all it requires (so far) is attribute window (an API
# SpectrumWindow object) and a function drawAll(handlers, object=None)
def printWindow(windowDraw, outputHandler, tick_location, tick_placement, aspectRatio, tick_font=None, major_minor_dict=None):

  # TBD: below just considers main window

  # first draw spectrum window
  windowPane = windowDraw.windowPane
  window = windowPane.spectrumWindow
  window.inOverrideMode = window.useOverrideRegion
  try:
    handlers = setupHandlers(windowPane, outputHandler, aspectRatio)
    outputHandler.save()
    if outputHandler.cHandler:
      windowDraw.drawAll(handlers)

  finally:
    del window.inOverrideMode

  drawInternalBorders(windowPane, outputHandler, aspectRatio)
  outputHandler.restore()

  printTicks(outputHandler, handlers, tick_location, tick_placement, tick_font, major_minor_dict)

  outputHandler.close()

def setupHandlers(windowPane, outputHandler, aspectRatio):

  handlers = []
  xAxisPanel = windowPane.findFirstAxisPanel(label='x')
  yAxisPanel = windowPane.findFirstAxisPanel(label='y')
  y = 0
  yRegions = yAxisPanel.sortedAxisRegions()
  yRegions.reverse()  # because need to flip y axis
  for yAxisRegion in yRegions:
    (y0, y1) = getAxisRegionRegion(yAxisRegion)
    if not windowPaneHasValueAxis(windowPane) and yAxisPanel.axisUnit.isBackwards:
      (y0, y1) = (y1, y0)
 
    ySize = aspectRatio*yAxisRegion.size
    row_handlers = []
    handlers.append(row_handlers)
    x = 0
    for xAxisRegion in xAxisPanel.sortedAxisRegions():
      (x0, x1) = getAxisRegionRegion(xAxisRegion)
      if (xAxisPanel.axisUnit.isBackwards):
        (x0, x1) = (x1, x0)
      xSize = xAxisRegion.size
      handler = PrintHandler(outputHandler, xRegion=(x0, x1), yRegion=(y0, y1),
                             xSize=xSize, ySize=ySize, x=x, y=y)
      row_handlers.append(handler)
      x = x + xSize
    y = y + ySize

  handlers.reverse() # because need to flip y axis

  return handlers
    
def printTicks(outputHandler, handlers, tick_location, tick_placement, tick_font, major_minor_dict=None):

  plot_size = outputHandler.plot_size
  outputHandler.newRange(0, 0, plot_size[0], plot_size[1])

  tickMajor = major_minor_dict.get('TickMajor')
  tickMinor = major_minor_dict.get('TickMinor')

  if Bottom in tick_placement:
    deltaMajor = major_minor_dict.get('XMajor')
    deltaMinor = major_minor_dict.get('XMinor')
    numberDecimals = major_minor_dict.get('XDecimal')
    for handler in handlers[0]:
      PrintTicks(handler, plot_size, tick_location, Bottom, deltaMajor, deltaMinor, numberDecimals, tick_font, tickMajor, tickMinor)

  if Top in tick_placement:
    deltaMajor = major_minor_dict.get('XMajor')
    deltaMinor = major_minor_dict.get('XMinor')
    numberDecimals = major_minor_dict.get('XDecimal')
    for handler in handlers[-1]:
      PrintTicks(handler, plot_size, tick_location, Top, deltaMajor, deltaMinor, numberDecimals, tick_font, tickMajor, tickMinor)

  if Left in tick_placement:
    deltaMajor = major_minor_dict.get('YMajor')
    deltaMinor = major_minor_dict.get('YMinor')
    numberDecimals = major_minor_dict.get('YDecimal')
    for row_handlers in handlers:
      handler = row_handlers[0]
      PrintTicks(handler, plot_size, tick_location, Left, deltaMajor, deltaMinor, numberDecimals, tick_font, tickMajor, tickMinor)

  if Right in tick_placement:
    deltaMajor = major_minor_dict.get('YMajor')
    deltaMinor = major_minor_dict.get('YMinor')
    numberDecimals = major_minor_dict.get('YDecimal')
    for row_handlers in handlers:
      handler = row_handlers[-1]
      PrintTicks(handler, plot_size, tick_location, Right, deltaMajor, deltaMinor, numberDecimals, tick_font, tickMajor, tickMinor)

def drawInternalBorders(windowPane, outputHandler, aspectRatio):

  outputHandler.newRange(0, 0, 1, 1)
  outputHandler.setColor((0.0, 0.0, 0.0))

  xAxisPanel = windowPane.findFirstAxisPanel(label='x')
  xAxisRegions = xAxisPanel.sortedAxisRegions()
  w = outputHandler.width
  yAxisPanel = windowPane.findFirstAxisPanel(label='y')
  yAxisRegions = yAxisPanel.sortedAxisRegions()
  h = outputHandler.height

  x = 0
  for xAxisRegion in xAxisRegions:
    width = xAxisRegion.size
    r0 = float(x) / w
    r1 = float(width) / w
    y = 0
    for yAxisRegion in yAxisRegions:
      height = aspectRatio*yAxisRegion.size
      s0 = float(y) / h
      s1 = float(height) / h
      outputHandler.drawRectangle((r0, s0, r1, s1))
      y += height
    x += width

def getPrintOption(appObject, key, defaultValue):

  name = 'printWin' + key
  if hasattr(appObject, name):
    value = getattr(appObject, name)
    if key in ('FileName', 'OtherHeight', 'OtherWidth', 'TickOutside') and value is None:
      value = defaultValue
  elif key in ('PositiveColors', 'NegativeColors'):
    value = appObject.root.application.getValues(appObject, keyword=name, defaultValues=defaultValue)
  else:
    value = appObject.root.application.getValue(appObject, keyword=name, defaultValue=defaultValue)

  return value

def setPrintOption(appObject, key, value):

  name = 'printWin' + key
  # Empty string not allowed
  if value == '':
    value = None

  if hasattr(appObject, name):
    try:
      setattr(appObject, name, value)
    except ApiError, e:
      if key != 'TickPlacement':
        print 'Warning: print option "%s" not stored in data model as "%s" because of exception: %s' % (key, value, str(e))

  elif key in ('PositiveColors', 'NegativeColors'):
    appObject.root.application.setValues(appObject, keyword=name, values=value)
  else:
    appObject.root.application.setValue(appObject, keyword=name, value=value)


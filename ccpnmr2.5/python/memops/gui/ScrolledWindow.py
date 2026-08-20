
"""
======================COPYRIGHT/LICENSE START==========================

ScrolledWindow.py: <write function here>

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)
- PDBe website (http://www.ebi.ac.uk/pdbe/)

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

import Tkinter

from memops.universal.Region1D import Region1D
from memops.universal.Region2D import Region2D
from memops.universal.Ticks import Ticks
from memops.universal.Util import isArray

from memops.gui.Canvas import Canvas
import memops.gui.Color as Color
from memops.gui.Frame import Frame
from memops.gui.Menu import Menu
from memops.gui.Scrollbar import Scrollbar
from memops.gui.Separator import Separator
from memops.universal.Util import OSButton

no_key_state = 0
shift_key_state = 1
ctrl_key_state = 4
alt_key_state = 8

# internal class
class WindowScrollbar(Scrollbar):

  def __init__(self, parent, *args, **kw):
 
    self.parent = parent

    self.old_lo = self.old_hi = None

    kw['callback'] = self.view
    kw['background'] = '#c9c9c9'
    Scrollbar.__init__(self, parent, *args, **kw)

  # color is a tuple
  # TBD: does nothing currently
  def setBackground(self, color):

    #self.config(troughcolor=color)

    color = Color.inverseGrey(color)
    #self.config(background=color)
    #self.config(activebackground=color)
    
  def setRegion(self, view_region):

    if (self.orient == Tkinter.HORIZONTAL):
      (lo, hi) = self.parent.xviewRegionToFraction(view_region)
      #print 'setRegion x', lo, hi
    else:
      (lo, hi) = self.parent.yviewRegionToFraction(view_region)
      (lo, hi) = (1.0-hi, 1.0-lo)
      #print 'setRegion y', lo, hi

    self.set(lo, hi)

  def view(self, lo, hi):

    if (self.orient == Tkinter.VERTICAL):
      (lo, hi) = (1.0 - hi, 1.0 - lo)
      
    # a bit dangerous but otherwise get infinite recursion
    # perhaps need to have add some tolerance
    if ((lo == self.old_lo) and (hi == self.old_hi)):
      return
    #print 'view1', lo, hi, self.old_lo, self.old_hi
    self.old_lo = lo
    self.old_hi = hi

    if (self.orient == Tkinter.HORIZONTAL):
      self.parent.xScroll(self, lo, hi)
    else:
      self.parent.yScroll(self, lo, hi)

# internal class
class TickCanvas(Canvas, Ticks):

  def __init__(self, parent, orient, size, tickScale=None,
               expandYAxisDecimals=True, *args, **kw):

    self.orient = orient
    self.crosshairs = []
    self.expandYAxisDecimals = expandYAxisDecimals

    defaultSize = 10 # gets overridden by other widgets in container
    # but need to do something otherwise canvas asks for too much space

    if (orient == Tkinter.HORIZONTAL):
      kw['width'] = defaultSize
      kw['height'] = size
    else:
      size = max(48, size)
      kw['width'] = size
      kw['height'] = defaultSize

    Canvas.__init__(self, parent, *args, **kw)

    if (orient == Tkinter.HORIZONTAL):
      # width arbitrary, will be changed in updateTicks()
      Ticks.__init__(self, width=defaultSize, height=size, isHorizontal=True, tickScale=tickScale)
    else:
      # height arbitrary, will be changed in updateTicks()
      Ticks.__init__(self, width=size, height=defaultSize, isHorizontal=False, tickScale=tickScale)

    self.bind('<Expose>', self.updateTicks)

  # color is a tuple
  def setBackground(self, color):

    self.config(background=color)
    self.fill = Color.inverseGrey(color)
    self.itemconfig('all', fill=self.fill)

  # overrides Ticks updateTicks()
  def updateTicks(self, *event):

    if (not self.region):
      return

    self.delete('all') # clear all previous items
    self.width = self.winfo_width()
    self.height = self.winfo_height()
    if (self.orient == Tkinter.HORIZONTAL):
      length = self.width
    else:
      length = self.height
    self.setMaxMajorMinor(length)

    Ticks.updateTicks(self)

    # 30 Mar 2010: below is hack to get the vertical tick canvas to
    # resize itself in width as you get more and more decimal places
    if not self.isHorizontal and self.majorTexts and self.expandYAxisDecimals:
      maxlen = max([len(x[0]) for x in self.majorTexts])
      #maxlen = max(3, maxlen)
      #width = 40 + (maxlen-3) * 8
      # 19 Aug 2010: change to below so that get less jiggling
      maxlen = max(5, maxlen)
      width = 48 + (maxlen-5) * 8
      self.config(width=width)
      self.width = width

  # overrides Ticks convertCoord()
  def convertCoord(self, r, r0, r1, d, haveInversion = False):

    c = Ticks.convertCoord(self, r, r0, r1, d, haveInversion)

    return int(c)

  # implements Ticks drawLine()
  def drawLine(self, x0, y0, x1, y1):

    return self.create_line(x0, y0, x1, y1, fill=self.fill)

  # implements Ticks drawText()
  def drawText(self, text, x, y, a, b):

    if (a == 0.0):
      if (b == 0.0):
        anchor = Tkinter.SW
      elif (b == 1.0):
        anchor = Tkinter.NW
      else:
        anchor = Tkinter.W
    elif (a == 1.0):
      if (b == 0.0):
        anchor = Tkinter.SE
      elif (b == 1.0):
        anchor = Tkinter.NE
      else:
        anchor = Tkinter.E
    else:
      if (b == 0.0):
        anchor = Tkinter.S
      elif (b == 1.0):
        anchor = Tkinter.N
      else:
        anchor = Tkinter.CENTER

    #print 'drawText', x, y, a, b, text, anchor
    self.create_text(x, y, text=text, anchor=anchor, font=self.font, fill=self.fill)

  def clearCrosshairs(self):

    for crosshair in self.crosshairs:
      self.delete(crosshair)

    self.crosshairs = []

  def drawCrosshairs(self, ts):

    if not self.region:
      return

    self.clearCrosshairs()
    horizontal = Tkinter.HORIZONTAL
    drawLine = self.drawLine
    
    w = self.winfo_width()
    h = self.winfo_height()
    (t0, t1) = self.region
    for t in ts:
      if self.orient == horizontal:
        x = self.convertCoord(t, t0, t1, w, False)
        self.crosshairs.append(drawLine(x, 0, x, h))
      else:
        y = self.convertCoord(t, t0, t1, h, True)
        self.crosshairs.append(drawLine(0, y, w, y))

# public class
# see note for WindowCanvas about superclass
#class SliceCanvas(Canvas):
class SliceCanvas(Frame):

  background = Color.white.hex

  def __init__(self, parent, orient, label, size, *args, **kw):

    #if (orient == Tkinter.HORIZONTAL):
    #  kw['height'] = size
    #else:
    #  kw['width'] = size

    #apply(Canvas.__init__, (self, parent) + args, kw)
    apply(Frame.__init__, (self, parent) + args, kw)

    self.orient = orient
    self.label = label
    self.size = None

    self.setSize(size)

  # need this because Tk does not seem to clean up bindings
  # and so a right-mouse click then causes an X11 crash
  def destroy(self):

    #print 'SliceCanvas destroy1', self.bind()
    for name in self.parent.slice_event_handlers.keys():
      #print 'SliceCanvas destroy2', name
      self.bind(name)

    Frame.destroy(self)
    #print 'SliceCanvas destroy3'

  def setRegion(self, view_region):

    self.view_region = view_region
    #if (self.orient == Tkinter.HORIZONTAL):
    #  pass
    #else:
    #  pass

  def setSize(self, size):

    if (size != self.size):

      self.size = size

      if (self.orient == Tkinter.HORIZONTAL):
        self.config(height=size)
      else:
        self.config(width=size)

  # color is a tuple
  def setBackground(self, color):

    self.background = color
    self.config(background=self.background)

#internal class
class CornerCanvas(Canvas):

  # TBD: does not work with splitting yet

  resize_button = 1
  move_button = OSButton(2)

  resize_delta = 2

  def __init__(self, parent, *args, **kw):

    apply(Canvas.__init__, (self, parent) + args, kw)

    self.xMultiplier = 1 # arbitrary
    self.yMultiplier = 1 # arbitrary

    # arbitrary sizes to start
    self.xline = self.create_line(0, 0, 1, 0)
    self.yline = self.create_line(0, 0, 0, 1)
    self.view_rects = []
    for j in range(self.parent.nrows):
      self.view_rects.append([])
      for i in range(self.parent.ncols):
        rect = self.create_rectangle(0, 0, 1, 1)
        self.view_rects[j].append(rect)
      
    self.bind('<Configure>', self.redraw)
    self.bind('<ButtonPress-%s>' % self.resize_button, self.markResize)
    self.bind('<B%s-Motion>' % self.resize_button, self.doResize)
    self.bind('<ButtonRelease-%s>' % self.resize_button, self.doResize)
    self.bind('<ButtonPress-%s>' % self.move_button, self.markMove)
    self.bind('<B%s-Motion>' % self.move_button, self.doMove)
    self.bind('<ButtonRelease-%s>' % self.move_button, self.doMove)

    self.state = no_key_state
    self.resizeCanvas = None
    self.moveCanvas = None

  def okDelta(self, u, up, v, v0, v1):

    if (abs(u-up) > self.resize_delta):
      return False

    if ((v0 <= v) and (v <= v1)):
      return True

    return False

  def findResizeCanvas(self, x, y):

    for j in range(self.parent.nrows):
      for i in range(self.parent.ncols):
        (x0, y0, x1, y1) = map(int, map(round, self.coords(self.view_rects[j][i])))
        if (self.okDelta(x, x0, y, y0, y1)):
          self.resizeType = 0
          self.resizeDirn = 'x'
        elif (self.okDelta(x, x1, y, y0, y1)):
          self.resizeType = 1
          self.resizeDirn = 'x'
        elif (self.okDelta(y, y0, x, x0, x1)):
          self.resizeType = 1 # note (because y from top in canvas)
          self.resizeDirn = 'y'
        elif (self.okDelta(y, y1, x, x0, x1)):
          self.resizeType = 0 # note
          self.resizeDirn = 'y'
        else:
          continue
        self.resizeCanvas = self.parent.canvases[j][i]
        return

    self.resizeCanvas = None

  def markResize(self, event):

    x = self.x = event.x
    y = self.y = event.y

    self.findResizeCanvas(x, y)

    if not self.resizeCanvas:
      return

  def calcXMultiplier(self, canvas):

    (w0, w1) = self.parent.world_region.x_region
    (v0, v1) = canvas.xview_region
    if v0 == v1:
      self.xMultiplier = 1 # arbitrary
    else:
      self.xMultiplier = (w1 - w0) / (v1 - v0)

  def calcYMultiplier(self, canvas):

    (w0, w1) = self.parent.world_region.y_region
    (v0, v1) = canvas.yview_region
    if v0 == v1:
      self.yMultiplier = 1 # arbitrary
    else:
      self.yMultiplier = (w1 - w0) / (v1 - v0)

  def doResize(self, event):

    if (not self.resizeCanvas):
      return

    if self.resizeDirn == 'x':
      w = self.winfo_width()
      dx = event.x - self.x
      self.calcXMultiplier(self.resizeCanvas)
      s = self.xMultiplier * dx / float(w)
      #print 'doResize1', s, self.resizeType
      self.parent.xResize(self.resizeCanvas, s, self.resizeType)
    else:
      h = self.winfo_height()
      dy = event.y - self.y
      self.calcYMultiplier(self.resizeCanvas)
      t = - self.yMultiplier * dy / float(h)
      #print 'doResize2', t, self.resizeType
      self.parent.yResize(self.resizeCanvas, t, self.resizeType)

    self.x = event.x
    self.y = event.y

  def findMoveCanvas(self, x, y):

    if ((self.parent.nrows == 1) and (self.parent.ncols == 1)):
      return self.parent.canvases[0][0]

    for j in range(self.parent.nrows):
      for i in range(self.parent.ncols):
        (x0, y0, x1, y1) = map(int, map(round, self.coords(self.view_rects[j][i])))
        if ((x0 <= x) and (x <= x1) and (y0 <= y) and (y <= y1)):
          return self.parent.canvases[j][i]

    return None

  def markMove(self, event):

    x = self.x = event.x
    y = self.y = event.y
    self.state = event.state

    moveCanvas = self.moveCanvas = self.findMoveCanvas(x, y)

    if moveCanvas:
      return

    self.calcXMultiplier(moveCanvas)
    self.calcYMultiplier(moveCanvas)

  def doMove(self, event):
  
    moveCanvas = self.moveCanvas
    if not moveCanvas:
      return
      
    parent = self.parent
    if self.state == shift_key_state:

      dy = event.y - self.y
      h = self.winfo_height()
      t = dy / float(h)
      scale = math.exp(-parent.speed*t)

      parent.zoom(moveCanvas, scale)

    else:

      w = self.winfo_width()
      h = self.winfo_height()
      dx = event.x - self.x
      dy = event.y - self.y
      s = self.xMultiplier * dx / float(w)
      t = - self.yMultiplier * dy / float(h)

      parent.translate(moveCanvas, s, t)

    self.x = event.x
    self.y = event.y

  def redraw(self, *event):

    w = self.winfo_width()
    h = self.winfo_height()
    parent = self.parent
    
    coords =self.coords 
    coords(self.xline, 0, 0, w, 0)
    coords(self.yline, 0, 0, 0, h)

    (s0, t0, s1, t1) = parent.world_region
    rows = range(parent.nrows)
    cols = range(parent.ncols)
    canvases = self.parent.canvases
    view_rects = self.view_rects
    
    for j in rows:
      (v0, v1) = canvases[j][0].yview_region
      y0 = h - int(round(h * (v0 - t0) / (t1 - t0)))
      y1 = h - int(round(h * (v1 - t0) / (t1 - t0)))
      
      for i in cols:
        (v0, v1) = canvases[j][i].xview_region
        x0 = int(round(w * (v0 - s0) / (s1 - s0)))
        x1 = int(round(w * (v1 - s0) / (s1 - s0)))
        coords(view_rects[j][i], x0, y0, x1, y1)

  # color is a tuple
  def setBackground(self, color):

    parent = self.parent
    itemconfig = self.itemconfig
    self.config(background=color)
    outline = fill = Color.inverseGrey(color)
    itemconfig(self.xline, fill=fill)
    itemconfig(self.yline, fill=fill)
    rows = range(parent.nrows)
    cols = range(parent.ncols)
    view_rects = self.view_rects

    for j in rows:
      for i in cols:
        itemconfig(view_rects[j][i], outline=outline)

# internal class
# should perhaps be subclass of Canvas but that does not
# seem to recognise keyboard events
#class WindowCanvas(Canvas):
class WindowCanvas(Frame):

  def __init__(self, parent, xview_region, yview_region, *args, **kw):

    #apply(Canvas.__init__, (self, parent) + args, kw)
    apply(Frame.__init__, (self, parent) + args, kw)

    self.xview_region = xview_region
    self.yview_region = yview_region

  # need this because Tk does not seem to clean up bindings
  # and so a right-mouse click then causes an X11 crash
  def destroy(self):

    #print 'WindowCanvas destroy1', self.bind()
    for name in self.parent.event_handlers.keys():
      #print 'WindowCanvas destroy2', name
      self.bind(name)

    Frame.destroy(self)
    #print 'WindowCanvas destroy3'

  def setXviewRegion(self, xview_region):

    self.xview_region = xview_region
    self.updateView()

  def setYviewRegion(self, yview_region):

    self.yview_region = yview_region
    self.updateView()

  def updateView(self):

    #print 'WindowCanvas updateView', self.parent.getCanvasRowCol(self), self.xview_region, self.yview_region
    if (self.parent.view_set_func):
      self.parent.view_set_func(self, self.xview_region, self.yview_region)

scrollName = 'scroll'
tickName = 'tick'
sliceName = 'slice'

# world_region and x/yview_region defined in whatever units drawing picture in

# public class
class ScrolledWindow(Frame):

  background = Color.white.hex

  def __init__(self, parent, ncols = 1, nrows = 1,
               world_region = None, # initial world size
               xview_regions = None, # initial x view size(s), Region1D
               yview_regions = None, # initial y view size(s)
               width = None, height = None, # initial width/height in pixels, can be tuples
               show_xcrosshair = True, show_ycrosshair = True,
               show_xticks = True, show_yticks = True,
               xtick_size = 30, ytick_size = 40,
               show_xslices = False, show_yslices = False,
               xslice_size = 100, yslice_size = 100,
               show_xscrollbar = True, show_yscrollbar = True,
               view_set_func = None, key_func = None,
               configure_func = None, expose_func = None, motion_func = None,
               enter_func = None, leave_func = None,
               slice_configure_func = None, slice_expose_func = None, slice_motion_func = None,
               slice_enter_func = None, slice_leave_func = None, slice_key_func = None,
               select_single_func = None, select_multi_func = None, drag_func = None,
               units_scroll = 0.1, pages_scroll = 1.0, # fraction for scrolling
               scroll_borderwidth = 1,
               separator_button = 1, separator_thickness = 4,
               get_aspect_ratio = None, get_geometry = None, set_geometry = None,
               has_value_axis = False, min_region_event_delta = 20,
               expand_y_axis_tick_decimals = True, *args, **kw):

    self.isDeleting = False

    if (world_region is None):
      world_region = Region2D()

    if (xview_regions):
      if (isArray(xview_regions)):
        assert ncols == len(xview_regions), 'ncols = %d, len(xview_regions) = %d' % (ncols, len(xview_regions))
      else:
        xview_region = xview_regions
        xview_regions = ncols * [0]
        for i in range(ncols):
          xview_regions[i] = xview_region.copy()
    else:
      xview_regions = ncols * [0]
      for i in range(ncols):
        xview_regions[i] = Region1D()

    if (yview_regions):
      if (isArray(yview_regions)):
        assert nrows == len(yview_regions), 'nrows = %d, len(yview_regions) = %d' % (nrows, len(yview_regions))
      else:
        yview_region = yview_regions
        yview_regions = nrows * [0]
        for i in range(nrows):
          yview_regions[i] = yview_region.copy()
    else:
      yview_regions = nrows * [0]
      for i in range(nrows):
        yview_regions[i] = Region1D()

    if (not width):
      width = int(600 / ncols)

    if (not height):
      height = int(600 / nrows)

    if (isArray(width)):
      assert ncols == len(width), 'ncols = %d, len(width) = %d' % (ncols, len(width))
      width = list(width)
    else:
      width = ncols * [width]

    if (isArray(height)):
      assert nrows == len(height), 'nrows = %d, len(height) = %d' % (nrows, len(height))
      height = list(height)
    else:
      height = nrows * [height]

    self.nrows = nrows
    self.ncols = ncols

    #self.xpanel_names = [ scrollName ]
    #self.ypanel_names = [ scrollName ]
    self.xpanel_names = [ tickName, sliceName, scrollName ]
    self.ypanel_names = [ tickName, sliceName, scrollName ]
    #self.xpanel_names = [ sliceName, tickName, scrollName ]
    #self.ypanel_names = [ sliceName, tickName, scrollName ]
    self.xpanel_dict = {}
    self.ypanel_dict = {}

    self.show_xcrosshair = show_xcrosshair
    self.show_ycrosshair = show_ycrosshair

    self.view_set_func = view_set_func

    if (has_value_axis):
      ytick_size = ytick_size + 15
    self.xtick_size = xtick_size
    self.ytick_size = ytick_size
    self.xslice_size = xslice_size
    self.yslice_size = yslice_size

    self.units_scroll = units_scroll
    self.pages_scroll = pages_scroll
    self.scroll_borderwidth = scroll_borderwidth

    self.separator_button = separator_button
    self.separator_thickness = separator_thickness

    self.select_single_func = select_single_func
    self.select_multi_func = select_multi_func
    self.drag_func = drag_func

    self.has_value_axis = has_value_axis
    self.min_region_event_delta = min_region_event_delta

    self.expand_y_axis_tick_decimals = expand_y_axis_tick_decimals

    self.x = None
    self.y = None
    self.state = no_key_state
    self.pressFuncs = {}
    self.motionFuncs = {}
    self.releaseFuncs = {}

    self.event_handlers = {}
    self.slice_event_handlers = {}

    self.get_aspect_ratio = get_aspect_ratio
    self.get_geometry = get_geometry
    self.set_geometry = set_geometry

    apply(Frame.__init__, (self, parent) + args, kw)
 
    #self.grid_rowconfigure(0, weight=1)
    #self.grid_columnconfigure(0, weight=1)

    # below no good because the empty lists all refer to the same list
    # self.canvases = nrows * [[]]
    self.canvases = []
    for j in range(nrows):
      self.canvases.append([])

    for j in range(nrows):
      for i in range(ncols):
        c = WindowCanvas(self, xview_regions[i], yview_regions[j],
                         width=width[i], height=height[j],
                         background=self.background)
        # background overwritten by GL (if GL used)

        self.canvases[j].append(c)

        #c.bind('<Motion>', lambda event, j=j, i=i: self.motion(j, i, event))
        #c.bind('<Configure>', lambda event, j=j, i=i: self.configure(j, i, event))
        #c.bind('<Expose>', lambda event, j=j, i=i: self.expose(j, i, event))
        #c.bind('<Leave>', lambda event, j=j, i=i: self.leave(j, i, event))

    self.xseparators = []
    for i in range(ncols-1):
      self.xseparators.append(Separator(self, orient=Tkinter.VERTICAL,
                                        width=self.separator_thickness))

    self.yseparators = []
    for j in range(nrows-1):
      self.yseparators.append(Separator(self, orient=Tkinter.HORIZONTAL,
                                        height=self.separator_thickness))

    self.xscrollbars = []
    for i in range(ncols):
      s = WindowScrollbar(self, orient=Tkinter.HORIZONTAL,
                          units_scroll=units_scroll, pages_scroll=pages_scroll,
                          borderwidth=scroll_borderwidth)
      self.xscrollbars.append(s)

    self.xpanel_dict[scrollName] = (self.xscrollbars, show_xscrollbar)

    self.yscrollbars = []
    for j in range(nrows):
      s = WindowScrollbar(self, orient=Tkinter.VERTICAL,
                          units_scroll=units_scroll, pages_scroll=pages_scroll,
                          borderwidth=scroll_borderwidth)
      self.yscrollbars.append(s)

    self.ypanel_dict[scrollName] = (self.yscrollbars, show_yscrollbar)

    self.xticks = []
    for i in range(ncols):
      c = TickCanvas(self, Tkinter.HORIZONTAL, xtick_size)
      self.xticks.append(c)

    self.xpanel_dict[tickName] = (self.xticks, show_xticks)

    self.yticks = []
    for j in range(nrows):
      c = TickCanvas(self, Tkinter.VERTICAL, ytick_size, tickScale=xtick_size,
                     expandYAxisDecimals=expand_y_axis_tick_decimals)
      self.yticks.append(c)

    self.ypanel_dict[tickName] = (self.yticks, show_yticks)

    self.xslices = []
    for i in range(ncols):
      c = SliceCanvas(self, Tkinter.HORIZONTAL, 'x', xslice_size)
      self.xslices.append(c)

    self.xpanel_dict[sliceName] = (self.xslices, show_xslices)

    self.yslices = []
    for j in range(nrows):
      c = SliceCanvas(self, Tkinter.VERTICAL, 'y', yslice_size)
      self.yslices.append(c)

    self.ypanel_dict[sliceName] = (self.yslices, show_yslices)

    self.corner_canvas = CornerCanvas(self, width=1, height=1)

    self.update_idletasks()
    self.gridAll()

    self.separatorBind()

    self.enter_func = enter_func
    self.canvasBind('<Enter>', self.enter)

    self.leave_func = leave_func
    self.canvasBind('<Leave>', self.leave)

    self.slice_enter_func = slice_enter_func
    self.sliceBind('<Enter>', self.slice_enter)

    self.slice_leave_func = leave_func
    self.sliceBind('<Leave>', self.slice_leave)

    if (configure_func):
      self.canvasBind('<Configure>', configure_func)

    if (expose_func):
      self.canvasBind('<Expose>', expose_func)

    if (motion_func):
      self.canvasBind('<Motion>', motion_func)

    if (key_func):
      self.canvasBind('<KeyPress>', key_func)

    if (slice_configure_func):
      self.sliceBind('<Configure>', slice_configure_func)

    if (slice_expose_func):
      self.sliceBind('<Expose>', slice_expose_func)

    if (slice_motion_func):
      self.sliceBind('<Motion>', slice_motion_func)

    if (slice_key_func):
      self.sliceBind('<KeyPress>', slice_key_func)

    self.setWorldRegion(world_region, xview_regions, yview_regions)

    self.menu = Menu(self, tearoff=0, include_event=True)
    self.sliceMenu = Menu(self, tearoff=0, include_event=True)

  def setExpandYAxisTickDecimals(self, expand_y_axis_tick_decimals):

    self.expand_y_axis_tick_decimals = expand_y_axis_tick_decimals
    for tick in self.yticks:
      tick.expand_y_axis_tick_decimals = expand_y_axis_tick_decimals

  # color is a tuple
  def setBackground(self, color):

    self.background = color
    self.config(background=self.background)

    for i in range(self.ncols):
      slice = self.xslices[i]
      slice.setBackground(color)
      ticks = self.xticks[i]
      ticks.setBackground(color)
      scroll = self.xscrollbars[i]
      scroll.setBackground(color)
 
    for j in range(self.nrows):
      slice = self.yslices[j]
      slice.setBackground(color)
      ticks = self.yticks[j]
      ticks.setBackground(color)
      scroll = self.yscrollbars[j]
      scroll.setBackground(color)
 
    self.corner_canvas.setBackground(color)

    color = Color.inverseGrey(color)
 
    for i in range(self.ncols-1):
      separator = self.xseparators[i]
      separator.setColor(color)
 
    for j in range(self.nrows-1):
      separator = self.yseparators[j]
      separator.setColor(color)

  def separatorBind(self):

    for i in range(self.ncols-1):
      self.separatorBindCol(i)

    for j in range(self.nrows-1):
      self.separatorBindRow(j)

  def separatorBindCol(self, col):

    b = self.separator_button
    '''
    separator = self.xseparators[col]
    separator.canvas.bind('<ButtonPress-%s>' % b,
                          lambda event: self.markXSeparator(event, separator))
    separator.canvas.bind('<ButtonRelease-%s>' % b,
                          lambda event: self.moveXSeparator(event, separator))
    '''

  def separatorBindRow(self, row):

    b = self.separator_button
    '''
    separator = self.yseparators[row]
    separator.canvas.bind('<ButtonPress-%s>' % b,
                          lambda event: self.markYSeparator(event, separator))
    separator.canvas.bind('<ButtonRelease-%s>' % b,
                          lambda event: self.moveYSeparator(event, separator))
    '''

  def translateBind(self, button, state = no_key_state):

    self.setupBind(button, self.markLocation, self.translateEvent, self.unmarkLocation, state=state)

  def zoomBind(self, button, state = no_key_state):

    self.setupBind(button, self.markLocation, self.zoomEvent, self.unmarkLocation, state=state)

  def regionBind(self, button, state = no_key_state):

    self.setupBind(button, self.markLocation, self.dragBox, self.regionEvent, state=state)

  def selectSingleBind(self, button, state = no_key_state):

    self.setupBind(button, self.markLocation,
                   releaseFunc=self.selectLocation, state=state)

  def selectMultiBind(self, button, state = no_key_state):

    self.setupBind(button, self.markLocation, self.dragBox,
                   self.selectRegion, state=state)

  def menuBind(self, button, menu_items = None, update_func = None):

    if (menu_items is None):
      menu_items = []

    self.menu.setMenuItems(menu_items)
    self.setupBind(button, self.menu.popupMenu)

    if (update_func):
      self.menu.config(postcommand=update_func)

  def sliceMenuBind(self, button, menu_items = None, update_func = None):

    if (menu_items is None):
      menu_items = []

    self.sliceMenu.setMenuItems(menu_items)
    s = '<ButtonPress-' + str(button) + '>'
    self.sliceBind(s, self.sliceMenu.popupMenu)

    if (update_func):
      self.sliceMenu.config(postcommand=update_func)

  def pressFunc(self, event, button, state):

    #self.menu.unpost()

    #state = event.state
    key = (button, state)
    #print 'pressFunc', key
    if (not self.pressFuncs.has_key(key)):
      state = no_key_state
      key = (button, state)
    self.state = state
    func = self.pressFuncs[key]
    if (func):
      event.state = state
      func(event)

  def motionFunc(self, event, button, state):

    key = (button, self.state)
    #print 'motionFunc', key
    func = self.motionFuncs[key]
    if (func):
      event.state = self.state
      func(event)

  def releaseFunc(self, event, button, state):

    key = (button, self.state)
    #print 'releaseFunc', key
    func = self.releaseFuncs[key]
    if (func):
      event.state = self.state
      func(event)

  def setupBind(self, button, pressFunc = '', motionFunc = '',
                releaseFunc = '', state = no_key_state):

    t = []
    if (state & shift_key_state):
      t.append('Shift')
    if (state & ctrl_key_state):
      t.append('Control')
    t = '-'.join(t)
    if (t):
      t = t + '-'

    self.pressFuncs[(button, state)] = pressFunc
    self.motionFuncs[(button, state)] = motionFunc
    self.releaseFuncs[(button, state)] = releaseFunc

    # sButton = OSButton(button)

    s = '<%sButtonPress-%s>' % (t, button)
    cmd = lambda event, button=button, state=state: self.pressFunc(event, button, state)
    self.canvasBind(s, cmd)

    s = '<%sB%s-Motion>' % (t, button)
    cmd = lambda event, button=button, state=state: self.motionFunc(event, button, state)
    self.canvasBind(s, cmd)

    s = '<%sButtonRelease-%s>' % (t, button)
    cmd = lambda event, button=button, state=state: self.releaseFunc(event, button, state)
    self.canvasBind(s, cmd)

  def sliceBind(self, name, handler = ''):

    self.slice_event_handlers[name] = handler

    for i in range(self.ncols):
      self.xslices[i].bind(name, handler)

    for j in range(self.nrows):
      self.yslices[j].bind(name, handler)

  def canvasBind(self, name, handler = ''):

    self.event_handlers[name] = handler

    for j in range(self.nrows):
      for i in range(self.ncols):
        self.canvases[j][i].bind(name, handler)

  def newRowBind(self, row):

    for name in self.event_handlers.keys():
      handler = self.event_handlers[name]
      for i in range(self.ncols):
        self.canvases[row][i].bind(name, handler)

    for name in self.slice_event_handlers.keys():
      handler = self.slice_event_handlers[name]
      self.yslices[row].bind(name, handler)

    self.separatorBindRow(row-1)

  def newColBind(self, col):

    for name in self.event_handlers.keys():
      handler = self.event_handlers[name]
      for j in range(self.nrows):
        self.canvases[j][col].bind(name, handler)

    for name in self.slice_event_handlers.keys():
      handler = self.slice_event_handlers[name]
      self.xslices[col].bind(name, handler)

    self.separatorBindCol(col-1)

  def setIsShownXScrollbar(self, isShown):

    (panel, shown) = self.xpanel_dict[scrollName]

    if (shown != isShown):
      self.xpanel_dict[scrollName] = (panel, isShown)
      self.gridAll()

  def setIsShownYScrollbar(self, isShown):

    (panel, shown) = self.ypanel_dict[scrollName]

    if (shown != isShown):
      self.ypanel_dict[scrollName] = (panel, isShown)
      self.gridAll()

  def setIsShownXSlice(self, isShown):

    (panel, shown) = self.xpanel_dict[sliceName]

    if (shown != isShown):
      self.xpanel_dict[sliceName] = (panel, isShown)
      self.gridAll()

  def setIsShownYSlice(self, isShown):

    (panel, shown) = self.ypanel_dict[sliceName]

    if (shown != isShown):
      self.ypanel_dict[sliceName] = (panel, isShown)
      self.gridAll()

  '''
  def zeroRowWeights(self, row):

    for j in range(self.nrows):
      if (row != j):
        self.grid_rowconfigure(2*j, weight=0, minsize=1)

    for j in range(self.nrows):
      for i in range(self.ncols):
        self.canvases[j][i].grid(row=2*j, column=2*i, sticky=Tkinter.NSEW)

  def zeroColWeights(self, col):

    for i in range(self.ncols):
      if (col != i):
        self.grid_columnconfigure(2*i, weight=0, minsize=1)

    for j in range(self.nrows):
      for i in range(self.ncols):
        self.canvases[j][i].grid(row=2*j, column=2*i, sticky=Tkinter.NSEW)

  def oneRowWeights(self):

    for j in range(self.nrows):
      self.grid_rowconfigure(2*j, weight=1, minsize=1)

    for j in range(self.nrows):
      for i in range(self.ncols):
        self.canvases[j][i].grid(row=2*j, column=2*i, sticky=Tkinter.NSEW)

  def oneColWeights(self):

    for i in range(self.ncols):
      self.grid_columnconfigure(2*i, weight=1, minsize=1)

    for j in range(self.nrows):
      for i in range(self.ncols):
        self.canvases[j][i].grid(row=2*j, column=2*i, sticky=Tkinter.NSEW)

  def zeroWeights(self, row = None, col = None):

    #print 'zeroWeights', row, col

    for j in range(self.nrows):
      if (row == j):
        self.grid_rowconfigure(2*j, weight=1)
      else:
        self.grid_rowconfigure(2*j, weight=0)

    for i in range(self.ncols):
      if (col == i):
        self.grid_columnconfigure(2*i, weight=1)
      else:
        self.grid_columnconfigure(2*i, weight=0)

  def oneWeights(self):

    #print 'oneWeights'

    for j in range(self.nrows):
      self.grid_rowconfigure(2*j, weight=1)

    for i in range(self.ncols):
      self.grid_columnconfigure(2*i, weight=1)
  '''

  #def gridAll(self, row=None, col=None):
  def gridAll(self):

    h = self.winfo_height()
    w = self.winfo_width()
    grid_columnconfigure = self.grid_columnconfigure
    grid_rowconfigure = self.grid_rowconfigure
    
    ncols = self.ncols
    nrows = self.nrows
    canvases = self.canvases
    
    rowspan = 2*nrows - 1
    for name in self.xpanel_names:
      (p, shown) = self.xpanel_dict[name]
      if (shown):
        rowspan = rowspan + 1

    columnspan = 2*ncols - 1
    for name in self.ypanel_names:
      (p, shown) = self.ypanel_dict[name]
      if (shown):
        columnspan = columnspan + 1

    rows = range(nrows)
    cols = range(ncols)
    for j in range(rowspan):
      grid_rowconfigure(j, weight=0, minsize=1)
    #for j in range(2*nrows, rowspan):
    #  self.grid_rowconfigure(j, weight=0, minsize=1)
    #for j in range(nrows-1):
    #  self.grid_rowconfigure(2*j+1, weight=0)

    for i in range(columnspan):
      grid_columnconfigure(i, weight=0, minsize=1)
    #for i in range(2*ncols, columnspan):
    #  self.grid_columnconfigure(i, weight=0, minsize=1)
    #for i in range(ncols-1):
    #  self.grid_columnconfigure(2*i+1, weight=0)

    for j in rows:
      grid_rowconfigure(2*j, weight=1, minsize=1)
      #if (row is None or row == j):
      #  self.grid_rowconfigure(2*j, weight=1, minsize=1)
      #self.grid_rowconfigure(2*j, minsize=1)
      #self.grid_rowconfigure(2*j, weight=1)

    for i in cols:
      grid_columnconfigure(2*i, weight=1, minsize=1)
      #if (col is None or col == i):
      #  self.grid_columnconfigure(2*i, weight=1, minsize=1)
      #self.grid_columnconfigure(2*i, minsize=1)
      #self.grid_columnconfigure(2*i, weight=1)

    for j in rows:
      for i in cols:
        canvases[j][i].grid(row=2*j, column=2*i, sticky=Tkinter.NSEW)

    for i in range(ncols-1):
      self.xseparators[i].grid(row=0, column=2*i+1, rowspan=rowspan, sticky=Tkinter.NS)
      w = w - self.xseparators[i].winfo_width()

    for j in range(nrows-1):
      self.yseparators[j].grid(row=2*j+1, column=0, columnspan=columnspan, sticky=Tkinter.EW)
      h = h - self.yseparators[j].winfo_height()

    row = 2*nrows - 1
    for name in self.xpanel_names:
      (p, shown) = self.xpanel_dict[name]
      ###if (self.has_value_axis and (name == sliceName)):
      ###  shown = False
      if (shown):
        if (isArray(p)):
          for i in cols:
            p[i].grid(row=row, column=2*i, sticky=Tkinter.EW)
          h = h - p[0].winfo_height()
        else:
          p.grid(row=row, column=0, columnspan=2*ncols-1, sticky=Tkinter.EW)
          h = h - p.winfo_height()
        row = row + 1
      else:
        if (isArray(p)):
          for i in cols:
            p[i].grid_forget()
        else:
          p.grid_forget()

    col = 2*ncols - 1
    for name in self.ypanel_names:
      (p, shown) = self.ypanel_dict[name]
      if (self.has_value_axis and (name == sliceName)):
        shown = False
      if (shown):
        if (isArray(p)):
          for j in rows:
            p[j].grid(row=2*j, column=col, sticky=Tkinter.NS)
          w = w - p[0].winfo_width()
        else:
          p.grid(row=0, column=col, rowspan=2*nrows-1, sticky=Tkinter.NS)
          w = w - p.winfo_width()
        col = col + 1
      else:
        if (isArray(p)):
          for j in rows:
            p[j].grid_forget()
        else:
          p.grid_forget()

    h = h / nrows
    w = w / ncols
    for j in rows:
      for i in cols:
        canvases[j][i].config(width=w, height=h)

    rowspan = row - 2*nrows + 1
    colspan = col - 2*ncols + 1
    if ((rowspan > 0) and (colspan > 0)):
      self.corner_canvas.grid(row=2*nrows-1, column=2*ncols-1,
                              rowspan=rowspan, columnspan=colspan, sticky=Tkinter.NSEW)
    else:
      self.corner_canvas.grid_forget()

  def viewRegionToFraction(self, v0, v1, w0, w1):

    d = float(w1 - w0)
    lo = (v0 - w0) / d
    hi = (v1 - w0) / d

    #print 'viewRegionToFraction', v0, v1, w0, w1, lo, hi

    return (lo, hi)

  def xviewRegionToFraction(self, xview_region):

    (x0, y0, x1, y1) = self.world_region
    (a0, a1) = xview_region

    return self.viewRegionToFraction(a0, a1, x0, x1)

  def yviewRegionToFraction(self, yview_region):

    (x0, y0, x1, y1) = self.world_region
    (b0, b1) = yview_region

    return self.viewRegionToFraction(b0, b1, y0, y1)

  def fractionToViewRegion(self, lo, hi, w0, w1):

    v0 = (w1 - w0) * lo + w0
    v1 = (w1 - w0) * hi + w0

    return Region1D(v0, v1)

  def fractionToXviewRegion(self, lo, hi):

    (x0, y0, x1, y1) = self.world_region

    return self.fractionToViewRegion(lo, hi, x0, x1)

  def fractionToYviewRegion(self, lo, hi):

    (x0, y0, x1, y1) = self.world_region

    return self.fractionToViewRegion(lo, hi, y0, y1)

  def xScroll(self, scrollbar, lo, hi):

    xview_region = self.fractionToXviewRegion(lo, hi)
    col = self.xscrollbars.index(scrollbar)

    #print 'xScroll1', lo, hi, col, xview_region
    for j in range(self.nrows):
      #print 'xScroll2', j, col
      self.canvases[j][col].setXviewRegion(xview_region)

    for name in self.xpanel_names:
      (panel, shown) = self.xpanel_dict[name]
      if (isArray(panel)):
        #print 'xScroll3', name, len(panel)
        panel[col].setRegion(xview_region)

    self.corner_canvas.redraw()

  def yScroll(self, scrollbar, lo, hi):

    yview_region = self.fractionToYviewRegion(lo, hi)
    #print 'yScroll', lo, hi, yview_region
    row = self.yscrollbars.index(scrollbar)

    for i in range(self.ncols):
      self.canvases[row][i].setYviewRegion(yview_region)

    for name in self.ypanel_names:
      (panel, shown) = self.ypanel_dict[name]
      if (isArray(panel)):
        panel[row].setRegion(yview_region)

    self.corner_canvas.redraw()

  def enter(self, event):

    #print 'enter', id(event.widget)
    event.widget.focus()
    #event.widget.haveFocus = True

    if (self.enter_func):
      self.enter_func(event)

  def leave(self, event):

    #print 'leave', id(event.widget)
    #event.widget.haveFocus = False

    if (self.leave_func):
      self.leave_func(event)

  def slice_enter(self, event):

    #print 'sliceEnter', id(event.widget)
    event.widget.focus()

    if (self.slice_enter_func):
      self.slice_enter_func(event)

  def slice_leave(self, event):

    if (self.slice_leave_func):
      self.slice_leave_func(event)

  def xResize(self, canvas, s, resize_type):

    if (self.get_aspect_ratio):
      r = self.get_aspect_ratio()
      if r is None:
        other_region = None
      else:
        other_region = canvas.yview_region
    else:
      r = other_region = None
    if (hasattr(canvas, 'ymax_extent')):
      ymax_extent = canvas.ymax_extent
    else:
      ymax_extent = None
    canvas.xview_region.resizeFraction(s, resize_type, r, other_region,
                                       self.world_region.y_region, ymax_extent)
    self.updateView(canvas)

  def yResize(self, canvas, t, resizeType):

    if (self.get_aspect_ratio):
      r = self.get_aspect_ratio()
      if r is None:
        other_region = None
      else:
        r = 1.0 / r
        other_region = canvas.xview_region
    else:
      r = other_region = None
    if (hasattr(canvas, 'xmax_extent')):
      xmax_extent = canvas.xmax_extent
    else:
      xmax_extent = None
    canvas.yview_region.resizeFraction(t, resizeType, r, other_region,
                                       self.world_region.x_region, xmax_extent)
    self.updateView(canvas)

  def markLocation(self, event):

    #print 'markLocation', event.state, event.x, event.y
    self.x = event.x
    self.y = event.y

  def unmarkLocation(self, event):

    self.x = None
    self.y = None

  def translateEvent(self, event):

    if self.x is None:
      self.markLocation(event)

    canvas = event.widget
    #print 'translateEvent1', self.getCanvasRowCol(canvas)

    w = canvas.winfo_width()
    h = canvas.winfo_height()
    dx = event.x - self.x
    dy = event.y - self.y
    s = -dx / float(w)
    t = dy / float(h)
    #print 'translateEvent2', s, t
    
    # TBD: Trial only 
    self.update_idletasks()
    self.markLocation(event)
    self.translate(canvas, s, t)

  def translate(self, canvas, s, t):

    world_region = self.world_region
    canvas.xview_region.translateFraction(s, world_region.x_region)
    canvas.yview_region.translateFraction(t, world_region.y_region)

    self.updateView(canvas)

  speed = 4 # TODO: does this need to be customisable?

  def zoomEvent(self, event):

    if self.y is None:
      self.markLocation(event)

    dy = event.y - self.y
    h = event.widget.winfo_height()
    t = dy / float(h)
    scale = math.exp(self.speed*t)
    
    # TBD: Trial only 
    self.update_idletasks()
    self.markLocation(event)
    self.zoom(event.widget, scale)

  def zoom(self, canvas, scale):

    if (hasattr(canvas, 'xmax_extent')):
      xmax_extent = canvas.xmax_extent
      ymax_extent = canvas.ymax_extent
    else:
      xmax_extent = None
      ymax_extent = None

    ''' old code
    # cannot just do below: all canvases need to be scaled because of aspect ratio
    #print 'ScrolledWindow zoom1', canvas.xview_region
    canvas.xview_region.zoom(scale, self.world_region.x_region, xmax_extent)
    #print 'ScrolledWindow zoom2', canvas.yview_region
    canvas.yview_region.zoom(scale, self.world_region.y_region, ymax_extent)
    #print 'ScrolledWindow zoom3', canvas.xview_region, canvas.yview_region

    self.updateView(canvas)
    '''
    
    world_region = self.world_region
    updateView = self.updateView
    canvases = self.canvases
    
    for i in range(self.ncols):
      canvases[0][i].xview_region.zoom(scale, world_region.x_region, xmax_extent)
      updateView(canvases[0][i])

    if not self.has_value_axis:
      for j in range(self.nrows):
        canvases[j][0].yview_region.zoom(scale, world_region.y_region, ymax_extent)
        updateView(canvases[j][0])

  def regionEvent(self, event):

    canvas = event.widget
    (row, col) = self.getCanvasRowCol(canvas)

    w = canvas.winfo_width()
    a0 = min(self.x, event.x)
    a1 = 1 + max(self.x, event.x)
    if (a1-a0) < self.min_region_event_delta:
      return

    h = canvas.winfo_height()
    b0 = h - 1 - self.y
    b1 = h - 1 - event.y
    (b0, b1) = (min(b0, b1), max(b0, b1))
    b1 = b1 + 1
    if (b1-b0) < self.min_region_event_delta:
      return

    # TBD below assumes that aspect ratio fixed
    if (a1-a0)*h > (b1-b0)*w:
      d = (0.5 * h * (a1-a0)) / w
      c = 0.5 * (b0+b1)
      b0 = c - d
      b1 = c + d
      if b0 < 0:
        b1 = b1 - b0
        b0 = 0
      elif b1 > h:
        b0 = h + b0 - b1
        b1 = h
      b0 = int(b0)
      b1 = int(b1)
    else:
      d = (0.5 * w * (b1-b0)) / h
      c = 0.5 * (a0+a1)
      a0 = c - d
      a1 = c + d
      if a0 < 0:
        a1 = a1 - a0
        a0 = 0
      elif a1 > w:
        a0 = w + a0 - a1
        a1 = w
      a0 = int(a0)
      a1 = int(a1)

    (x0, x1) = canvas.xview_region
    d = x1 - x0
    x1 = x0 + (d*a1)/w
    x0 = x0 + (d*a0)/w
    canvas.xview_region.set(x0, x1)
    self.updateView(canvas)

    (y0, y1) = canvas.yview_region
    d = y1 - y0
    y1 = y0 + (d*b1)/h
    y0 = y0 + (d*b0)/h
    canvas.yview_region.set(y0, y1)
    self.updateView(canvas)

    self.unmarkLocation(event)

  def selectLocation(self, event):

    if ((event.x == self.x) and (event.y == self.y) and self.select_single_func):
      (a, b, x, y) = self.calcWorldCoord(event.widget, event.x, event.y)
      button = event.num
      state = event.state & 255
      self.select_single_func(event.widget, a, b, x, y, button, state, event=event)

  def selectRegion(self, event):

    if (self.select_multi_func and self.x is not None and self.y is not None):
      (a0, b0, x0, y0) = self.calcWorldCoord(event.widget, self.x, self.y)
      (a1, b1, x1, y1) = self.calcWorldCoord(event.widget, event.x, event.y)
      button = event.num
      state = event.state & 255
      self.select_multi_func(event.widget, a0, b0, a1, b1, x0, y0, x1, y1, button, state, event=event)
      self.x = None
      self.y = None

  def dragBox(self, event):

    if self.x is None:
      self.markLocation(event)

    if (self.drag_func and self.x is not None and self.y is not None):
      widget = event.widget
      (a0, b0, x0, y0) = self.calcWorldCoord(widget, self.x, self.y)
      (a1, b1, x1, y1) = self.calcWorldCoord(widget, event.x, event.y)
      button = event.num
      state = event.state & 255
      self.drag_func(widget, a0, b0, a1, b1, x0, y0, x1, y1, button, state, event=event)

  # input: (x, y) = pixel coordinates (relative to view)
  # output: (a, b) = world coordinates
  # output: (x, y) = fraction of coordinates relative to screen
  def calcWorldCoord(self, canvas, x, y):

    x = (float(x)+0.5) / float(canvas.winfo_width())
    (a0, a1) = canvas.xview_region
    a = a0 + x*(a1-a0)

    y = 1.0 - (float(y)+0.5) / float(canvas.winfo_height())
    (b0, b1) = canvas.yview_region
    b = b0 + y*(b1-b0)

    return (a, b, x, y)

  def getSliceRow(self, slice):

    try:
      row = self.yslices.index(slice)
    except:
      raise ValueError, 'slice not in row list' # should never be here

    return row

  def getSliceCol(self, slice):

    try:
      col = self.xslices.index(slice)
    except:
      raise ValueError, 'slice not in col list' # should never be here

    return col

  def getCanvasRowCol(self, canvas):

    for row in range(self.nrows):
      try:
        col = self.canvases[row].index(canvas)
        break
      except:
        pass
    else:
      raise ValueError, 'canvas not in list' # should never be here

    #print 'row, col', row, col
    return (row, col)

  def updateView(self, canvas):

    (row, col) = self.getCanvasRowCol(canvas)

    #print 'updateView1', row, col
    for name in self.xpanel_names:
      (panel, shown) = self.xpanel_dict[name]
      if (isArray(panel)):
        #print 'updateView2', name, id(panel[col])
        panel[col].setRegion(canvas.xview_region)

    for name in self.ypanel_names:
      (panel, shown) = self.ypanel_dict[name]
      if (isArray(panel)):
        #print 'updateView3', name, id(panel[row])
        panel[row].setRegion(canvas.yview_region)

    self.canvases[row][col].updateView()
    #for j in range(self.nrows):
    #  #print 'updateView4', j, id(self.canvases[j][col])
    #  self.canvases[j][col].updateView()

    #for i in range(self.ncols):
    #  if (i != col):
    #    #print 'updateView5', i, id(self.canvases[row][i])
    #    self.canvases[row][i].updateView()

    self.corner_canvas.redraw()

  def setWorldRegion(self, world_region, xview_regions, yview_regions):

    ncols = self.ncols
    nrows = self.nrows

    assert ncols == len(xview_regions), 'ncols = %d, len(xview_regions) = %d' % (ncols, len(xview_regions))
    assert nrows == len(yview_regions), 'nrows = %d, len(yview_regions) = %d' % (nrows, len(yview_regions))

    self.world_region = world_region

    (x0, y0, x1, y1) = world_region

    assert x0 != x1, 'x0 = x1 = %f' % x0
    assert y0 != y1, 'y0 = y1 = %f' % y0

    rows = range(self.nrows)
    cols = range(self.ncols)

    for i in cols:

      (a0, a1) = xview_regions[i]
      if (x0 < x1):
        assert x0 <= a0 and a0 < a1 and a1 <= x1, 'x0 = %f, a0 = %f, a1 = %f, x1 = %f' % (x0, a0, a1, x1)
      else:
        assert x0 >= a0 and a0 > a1 and a1 >= x1, 'x0 = %f, a0 = %f, a1 = %f, x1 = %f' % (x0, a0, a1, x1)

    for j in rows:

      (b0, b1) = yview_regions[j]
      if (y0 < y1):
        assert y0 <= b0 and b0 < b1 and b1 <= y1, 'y0 = %f, b0 = %f, b1 = %f, y1 = %f' % (y0, b0, b1, y1)
      else:
        assert y0 >= b0 and b0 > b1 and b1 >= y1, 'y0 = %f, b0 = %f, b1 = %f, y1 = %f' % (y0, b0, b1, y1)

    #print 'setWorldRegion', x0, y0, x1, y1, a0, a1, b0, b1

    for name in self.xpanel_names:
      (panel, shown) = self.xpanel_dict[name]
      if (isArray(panel)):
        for i in cols:
          panel[i].setRegion(xview_regions[i])

    for name in self.ypanel_names:
      (panel, shown) = self.ypanel_dict[name]
      if (isArray(panel)):
        for j in rows:
          panel[j].setRegion(yview_regions[j])

    for j in rows:
      for i in cols:
        self.canvases[j][i].updateView()

    self.corner_canvas.redraw()

  def setXviewRegion(self, xview_region, col):

    (x0, y0, x1, y1) = self.world_region

    (a0, a1) = xview_region
    if (x0 < x1):
      assert x0 <= a0 and a0 < a1 and a1 <= x1, 'x0 = %f, a0 = %f, a1 = %f, x1 = %f' % (x0, a0, a1, x1)
    else:
      assert x0 >= a0 and a0 > a1 and a1 >= x1, 'x0 = %f, a0 = %f, a1 = %f, x1 = %f' % (x0, a0, a1, x1)

    for name in self.xpanel_names:
      (panel, shown) = self.xpanel_dict[name]
      if (isArray(panel)):
        panel[col].setRegion(xview_region)

    for j in range(self.nrows):
      #self.canvases[j][col].updateView()
      self.canvases[j][col].setXviewRegion(xview_region)

  def setYviewRegion(self, yview_region, row):

    (x0, y0, x1, y1) = self.world_region

    (b0, b1) = yview_region

    if (y0 < y1):
      assert y0 <= b0 and b0 < b1 and b1 <= y1, 'y0 = %f, b0 = %f, b1 = %f, y1 = %f' % (y0, b0, b1, y1)
    else:
      assert y0 >= b0 and b0 > b1 and b1 >= y1, 'y0 = %f, b0 = %f, b1 = %f, y1 = %f' % (y0, b0, b1, y1)

    for name in self.ypanel_names:
      (panel, shown) = self.ypanel_dict[name]
      if (isArray(panel)):
        panel[row].setRegion(yview_region)

    for i in range(self.ncols):
      #self.canvases[row][i].updateView()
      self.canvases[row][i].setYviewRegion(yview_region)

  def setXviewWidth(self, width, col):

    #self.parent.grid_propagate(False)
    #self.grid_propagate(False)

    #for j in range(self.nrows):
    #  for i in range(self.ncols):
    #    self.canvases[j][i].grid_propagate(False)

    #print 'setXviewWidth1', width, col
    # below does not work for scrollbars since width means thickness
    for name in self.xpanel_names:
      (panel, shown) = self.xpanel_dict[name]
      if (isArray(panel)):
        #print 'setXviewWidth2', name, panel[col].winfo_width(), panel[col].winfo_height()
        #panel[col].config(width=width)
        #print 'setXviewWidth3', name, panel[col].winfo_width(), panel[col].winfo_height()
        pass # TBD: for now

    for j in range(self.nrows):
      #print 'setXviewWidth4', j, col, self.canvases[j][col].winfo_width(), width
      self.canvases[j][col].config(width=width)
      #print 'setXviewWidth5', j, col, self.canvases[j][col].winfo_width()

    #for j in range(self.nrows):
    #  for i in range(self.ncols):
    #    self.canvases[j][i].grid_propagate(True)

    #self.grid_propagate(True)
    #self.parent.grid_propagate(True)

  def setSizeXSlice(self, size):

    (self.xslices, show_xslices) = self.xpanel_dict[sliceName]

    self.xslice_size = size
    for slice in self.xslices:
      slice.setSize(size)

  def setSizeYSlice(self, size):

    (self.yslices, show_yslices) = self.ypanel_dict[sliceName]

    self.yslice_size = size
    for slice in self.yslices:
      slice.setSize(size)

  def markXSeparator(self, event, separator):

    separator.x = event.x

  def markYSeparator(self, event, separator):

    separator.y = event.y

  def moveXSeparator(self, event, separator):

    x = event.x
    col = self.xseparators.index(separator)

    #print 'moveXSeparator1', col, x
    dx = x - separator.x
    separator.x = x

    w0 = self.canvases[0][col].winfo_width()
    w1 = self.canvases[0][col+1].winfo_width()
    #print 'moveXSeparator2', w0, w1
    self.setXviewWidth(w0+dx, col)
    self.setXviewWidth(w1-dx, col+1)

    #w0 = self.canvases[0][col].winfo_width()
    #w1 = self.canvases[0][col+1].winfo_width()
    #print 'moveXSeparator3', w0, w1

  def moveYSeparator(self, event, separator):

    pass

  def splitWindowHorizontal(self, event):

    pass

  def splitWindowVertical(self, event):

    canvas = event.widget
    x = event.x
    w = canvas.winfo_width()
    col = self.getCanvasRowCol(canvas)[1]
    (v0, v1) = canvas.xview_region
    r = float(x) / float(w-1)
    v = r*(v1-v0) + v0
    xview_region0 = Region1D(v0, v)
    xview_region1 = Region1D(v, v1)
    w0 = x - 1
    w1 = w - x - 1

    for j in range(self.nrows):
      rect = self.corner_canvas.create_rectangle(0, 0, 1, 1)
      self.corner_canvas.view_rects[j][col+1:col+1] = [rect]

    for j in range(self.nrows):
      cc = self.canvases[j][col]
      c = WindowCanvas(self, xview_region1, cc.yview_region,
                         width=w1, height=cc.winfo_height(),
                         background=self.background)
      # background overwritten by GL (if GL used)
      self.canvases[j][col+1:col+1] = [c]

    s = Separator(self, orient=Tkinter.VERTICAL, width=self.separator_thickness)
    self.xseparators[col:col] = [s]

    s = WindowScrollbar(self, orient=Tkinter.HORIZONTAL,
                        units_scroll=self.units_scroll, pages_scroll=self.pages_scroll,
                        borderwidth=self.scroll_borderwidth)
    self.xscrollbars[col+1:col+1] = [s]

    c = TickCanvas(self, Tkinter.HORIZONTAL, self.xtick_size)
    self.xticks[col+1:col+1] = [c]

    c = SliceCanvas(self, Tkinter.HORIZONTAL, 'x', self.xslice_size)
    self.xslices[col+1:col+1] = [c]

    self.newColBind(col+1)

    self.setXviewWidth(w0, col)
    self.setXviewRegion(xview_region0, col)

    self.setXviewWidth(w1, col+1)
    self.setXviewRegion(xview_region1, col+1)

    self.ncols = self.ncols + 1

    self.gridAll()

  def setRowHeight(self, row, h):

    canvas = self.canvases[row][0]
    h_orig = canvas.winfo_height()
    scale = h / float(h_orig)
    yview_region = canvas.yview_region
    if (hasattr(canvas, 'ymax_extent')):
      ymax_extent = canvas.ymax_extent
    else:
      ymax_extent = None
    yview_region.zoom(scale, self.world_region.y_region, ymax_extent)
    self.setYviewRegion(yview_region, row)

    for i in range(self.ncols):
      self.canvases[row][i].config(height=h)

  def addRow(self, doGrid = True):

    h = self.canvases[-1][0].winfo_height()
    """
    h = 0
    for j in range(self.nrows):
      h = h + self.canvases[j][0].winfo_height()
    h = max(1, int(h / (self.nrows+1)))
    for j in range(self.nrows):
      self.setRowHeight(j, h)
    """

    rects = []
    self.corner_canvas.view_rects.append(rects)
    for i in range(self.ncols):
      rect = self.corner_canvas.create_rectangle(0, 0, 1, 1)
      rects.append(rect)

    c = self.canvases[-1][0]
    yview_region = c.yview_region.copy()

    cs = []
    self.canvases.append(cs)
    for i in range(self.ncols):
      c = self.canvases[0][i]
      xview_region = c.xview_region.copy()
      w = c.winfo_width()
      c = WindowCanvas(self, xview_region, yview_region,
                       width=w, height=h, background=self.background)
      # background overwritten by GL (if GL used)
      cs.append(c)

    s = Separator(self, orient=Tkinter.HORIZONTAL, height=self.separator_thickness)
    self.yseparators.append(s)

    s = WindowScrollbar(self, orient=Tkinter.VERTICAL,
                        units_scroll=self.units_scroll, pages_scroll=self.pages_scroll,
                        borderwidth=self.scroll_borderwidth)
    self.yscrollbars.append(s)

    c = TickCanvas(self, Tkinter.VERTICAL, self.ytick_size,
                   expandYAxisDecimals=self.expand_y_axis_tick_decimals)
    self.yticks.append(c)

    c = SliceCanvas(self, Tkinter.VERTICAL, 'y', self.yslice_size)
    self.yslices.append(c)

    self.newRowBind(self.nrows)

    #if (self.get_geometry and self.set_geometry):
    #  (w, hh, x, y) = self.get_geometry()
    #  self.set_geometry(w, hh+h+self.separator_thickness, x, y)

    self.nrows = self.nrows + 1

    self.gridAll()

    if doGrid:
      self.gridAll()

  def deleteRow(self, row=-1, doGrid = True):

    if (self.nrows <= 1):                       
      return                                    

    self.isDeleting = True                      

    rects = self.corner_canvas.view_rects[row]  
    for rect in rects:                          
      self.corner_canvas.delete(rect)           
    del self.corner_canvas.view_rects[row]      

    cs = self.canvases[row]                     
    del self.canvases[row]                      
    for c in cs:                                
      if hasattr(c, 'handler'):                 
        del c.handler                           
      c.destroy()                               

    if (row == -1):
      seprow = -1
    else:
      seprow = max(row - 1, 0)
      
    s = self.yseparators[seprow]                
    del self.yseparators[seprow]                
    s.destroy()                                 
 
    s = self.yscrollbars[row]                   
    del self.yscrollbars[row]                   
    s.destroy()                                 

    c = self.yticks[row]                        
    del self.yticks[row]                        
    c.destroy()                                 
 
    c = self.yslices[row]                       
    del self.yslices[row]                       
    if hasattr(c, 'handler'):                   
      del c.handler                             
    c.destroy()                                 
 
    self.nrows = self.nrows - 1                 
 
    self.isDeleting = False                     
 
    if doGrid:
      self.gridAll()

  def setColWidth(self, col, w):

    canvas = self.canvases[0][col]
    w_orig = canvas.winfo_width()
    scale = w / float(w_orig)
    xview_region = canvas.xview_region
    if (hasattr(canvas, 'xmax_extent')):
      xmax_extent = canvas.xmax_extent
    else:
      xmax_extent = None
    xview_region.zoom(scale, self.world_region.x_region, xmax_extent)
    self.setXviewRegion(xview_region, col)

    for j in range(self.nrows):
      self.canvases[j][col].config(width=w)

  def addCol(self, col = -1, doGrid = True):

    #w = self.canvases[0][-1].winfo_width()
    w = 0
    for i in range(self.ncols):
      w = w + self.canvases[0][i].winfo_width()
    w = max(1, int(w / (self.ncols+1)))
    for i in range(self.ncols):
      self.setColWidth(i, w)

    for j in range(self.nrows):
      rect = self.corner_canvas.create_rectangle(0, 0, 1, 1)
      if (col == -1):
        self.corner_canvas.view_rects[j].append(rect)
      else:
        self.corner_canvas.view_rects[j][col:col] = [rect]

    c = self.canvases[0][-1]
    xview_region = c.xview_region.copy()

    for j in range(self.nrows):
      c = self.canvases[j][0]
      yview_region = c.yview_region.copy()
      h = c.winfo_height()
      c = WindowCanvas(self, xview_region, yview_region,
                       width=w, height=h, background=self.background)
      # background overwritten by GL (if GL used)
      if (col == -1):
        self.canvases[j].append(c)
      else:
        self.canvases[j][col:col] = [c]

    s = Separator(self, orient=Tkinter.VERTICAL, width=self.separator_thickness)
    self.xseparators.append(s)

    s = WindowScrollbar(self, orient=Tkinter.HORIZONTAL,
                        units_scroll=self.units_scroll, pages_scroll=self.pages_scroll,
                        borderwidth=self.scroll_borderwidth)
    if (col == -1):
      self.xscrollbars.append(s)
    else:
      self.xscrollbars[col:col] = [s]

    c = TickCanvas(self, Tkinter.HORIZONTAL, self.xtick_size)
    if (col == -1):
      self.xticks.append(c)
    else:
      self.xticks[col:col] = [c]

    c = SliceCanvas(self, Tkinter.HORIZONTAL, 'x', self.xslice_size)
    if (col == -1):
      self.xslices.append(c)
    else:
      self.xslices[col:col] = [c]

    #self.newColBind(self.ncols)
    self.newColBind(col)

    #if (self.get_geometry and self.set_geometry):
    #  (ww, h, x, y) = self.get_geometry()
    #  self.set_geometry(ww+w+self.separator_thickness, h, x, y)

    self.ncols = self.ncols + 1

    if doGrid:
      self.gridAll()

  def deleteCol(self, col = -1, doGrid = True):

    if (self.ncols <= 1):
      return

    self.isDeleting = True

    for j in range(self.nrows):
      rects = self.corner_canvas.view_rects[j]
      self.corner_canvas.delete(rects[col])
      del rects[col]

    for j in range(self.nrows):
      c = self.canvases[j][col]
      del self.canvases[j][col]
      if hasattr(c, 'handler'):
        del c.handler
      c.destroy()

    if (col == -1):
      sepcol = -1
    else:
      sepcol = max(col - 1, 0)
    s = self.xseparators[sepcol]
    del self.xseparators[sepcol]
    s.destroy()

    s = self.xscrollbars[col]
    del self.xscrollbars[col]
    s.destroy()

    c = self.xticks[col]
    del self.xticks[col]
    c.destroy()

    c = self.xslices[col]
    del self.xslices[col]
    if hasattr(c, 'handler'):
      del c.handler
    c.destroy()

    self.ncols = self.ncols - 1

    self.isDeleting = False

    if doGrid:
      self.gridAll()

  def popupMenu(self, event):

    self.menu.popupMenu(event)
    self.menu.focus_set()

  def setRowsTo(self, nrows = 1, doGrid = True):

    while self.nrows > nrows:
      self.deleteRow(doGrid=False) 

    while self.nrows < nrows:
      self.addRow(doGrid=False)

    if doGrid:
      self.gridAll()

  def setColsTo(self, ncols = 1, doGrid = True):

    while self.ncols > ncols:
      self.deleteCol(doGrid=False)

    while self.ncols < ncols:
      self.addCol(doGrid=False)

    if doGrid:
      self.gridAll()

if (__name__ == '__main__'):

  import sys
  from memops.gui.Button import Button

  highlightthickness = 0
  borderwidth = 0
  n = 0
  colors = ('red', 'green', 'blue')

  # does not work now that WindowCanvas is Frame rather than Canvas
  def addLine():

    global n

    c = s.canvases[0][0]
    w = c.winfo_width()
    h = c.winfo_height()
    col = colors[n % len(colors)]
    p = highlightthickness + borderwidth
    q = -1 - highlightthickness - borderwidth
    c.create_line(n+p, n+p, w+q-n, n+p, fill=col)
    c.create_line(w+q-n, n+p, w+q-n, h+q-n, fill=col)
    c.create_line(w+q-n, h+q-n, n+p, h+q-n, fill=col)
    c.create_line(n+p, h+q-n, n+p, n+p, fill=col)
    n = n + 1

  def printWidths():

    ws = []
    for i in range(s.ncols):
      c = s.canvases[0][i]
      w = c.winfo_width()
      ws.append(w)
    print ws
        
  def printHeights():

    hs = []
    for j in range(s.nrows):
      c = s.canvases[j][0]
      h = c.winfo_height()
      hs.append(h)
    print hs

  def change(dx):
    w0 = s.canvases[0][0].winfo_width()
    w1 = s.canvases[0][1].winfo_width()
    for i in range(s.nrows):
      s.canvases[i][0].config(width=w0+dx)
      s.canvases[i][1].config(width=w1-dx)
    
  def changePlus():
    change(5)

  def changeMinus():
    change(-5)

  def configure(event):
    print 'configure:', event.width, event.height, s.getCanvasRowCol(event.widget)

  def expose(event):
    #print 'expose:', event.x, event.y, event.width, event.height
    pass

  #def keypress(canvas, event):
  def keypress(event):
    print 'keypress'

  def view_set(canvas, xview, yview):
    #print 'view_set:', xview.t0, xview.t1, yview.t0, yview.t1
    pass

  def test(n, event):
    print 'test:', n, event.x, event.y

  root = Tkinter.Tk()

  xr = Region1D(-0.5, 1.5)
  yr = Region1D(-1.5, 1.5)
  wr = Region2D(xr, yr)

  root.minsize(1, 1)
  s = ScrolledWindow(root, ncols=2, nrows=3, world_region=wr,
                     configure_func = configure,
                     expose_func = expose,
                     key_func = keypress,
                     view_set_func = view_set)
  s.pack(expand=Tkinter.YES, fill=Tkinter.BOTH)

  s.translateBind(2)

  subitems = [ \
    { 'kind': 'command', 'label': 'subtest2', 'command': lambda event: test(2, event) },
    ]

  menu_items = [ \
    { 'kind': 'command', 'label': 'test1', 'command': lambda event: test(1, event) },
    { 'kind': 'separator' },
    { 'kind': 'cascade', 'label': 'test2', 'submenu': subitems },
    ]

  s.menuBind(3, menu_items)

  f = Frame(root)
  f.pack()

  #Button(f, text='add line', command=addLine).pack(side=Tkinter.LEFT)
  Button(f, text='widths', command=printWidths).pack(side=Tkinter.LEFT)
  Button(f, text='heights', command=printHeights).pack(side=Tkinter.LEFT)
  Button(f, text='+ 5', command=changePlus).pack(side=Tkinter.LEFT)
  Button(f, text='- 5', command=changeMinus).pack(side=Tkinter.LEFT)
  Button(f, text='quit', command=sys.exit).pack(side=Tkinter.LEFT)

  root.mainloop()

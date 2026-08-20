
"""
======================COPYRIGHT/LICENSE START==========================

Scrollbar.py: <write function here>

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
# not a subclass of Tkinter.Scrollbar
# but with approximately same functionality
# but allows resizing with middle button
# as well as usual scrolling with left button

import Tkinter

from memops.gui.Canvas import Canvas
from memops.gui.Color import scaleColor, invertColor
from memops.gui.Frame import Frame
from memops.universal.Util import buttonClick, buttonMotion


shift_key_state = 1

min_delta = 1.0e-5

min_bubble_size = 8

top_left_arrow_mode = 0
top_left_trough_mode = 1
bubble_mode = 2
bottom_right_trough_mode = 3
bottom_right_arrow_mode = 4

class Scrollbar(Frame):

  def __init__(self, parent,
               # arguments below are Tkinter.Scrollbar options
               orient = Tkinter.HORIZONTAL,
               width = 15, borderwidth = 1,
               background = None, troughcolor = None,
               repeatdelay = 300, # msecs before repeat action on arrows and troughs
               repeatinterval = 100, # msecs between repeat action on arrows and troughs
               # arguments below are not Tkinter.Scrollbar options
               # callback is approximately same as Tkinter.Scrollbar command
               # option but has different (more convenient) arguments
               allow_resize = False, callback = None, relief='sunken',
               show_text = False, text_color = '#000000', text_func = None,
               units_scroll = 0.1, pages_scroll = 1.0,
               min_thickness = None, *args, **kw):

    assert orient in (Tkinter.HORIZONTAL, Tkinter.VERTICAL)
    if (show_text):
      assert text_func is not None

    self.relief = relief
    self.orient = orient
    self.callback = callback
    self.borderwidth = borderwidth
    self.repeatdelay = repeatdelay
    self.repeatinterval = repeatinterval
    self.show_text = show_text
    self.text_func = text_func
    self.units_scroll = units_scroll
    self.pages_scroll = pages_scroll
    self.min_thickness = min_thickness
    
    self.mode = 0
    self.lo = 0.0
    self.hi = 1.0
    self.first_pass = True

    self.haveShift = False

    Frame.__init__(self, parent, *args, **kw)

    if (background is None):
      bg = parent.cget('bg')
    else:
      bg = background

    if (troughcolor is None):
      bg_trough = scaleColor(self, bg, 0.9)
    else:
      bg_trough = troughcolor

    # below are for shadows
    
    if relief == 'sunken':
      color1 = scaleColor(self, bg_trough, 0.7)
      color2 = scaleColor(self, bg_trough, 1.2)
      color3 = scaleColor(self, bg, 1.2)
      color4 = scaleColor(self, bg, 0.7)
      
    elif relief == 'raised':
      color1  = scaleColor(self, bg_trough, 1.2)
      color2 = scaleColor(self, bg_trough, 0.7)
      color3 = scaleColor(self, bg, 0.7)
      color4 = scaleColor(self, bg, 1.2)
      
    else:
      color1 = color2 = scaleColor(self, bg_trough, 0.7)
      color3 = color4 = scaleColor(self, bg, 1.2)

    self.brColor = color4
    self.tlColor = color3
    
    defaultSize = 10 # gets overriden by other widgets in container
    # but need to do something otherwise canvas asks for too much space

    kw2 = {}
    if (orient == Tkinter.HORIZONTAL):
      kw2['width'] = defaultSize
      kw2['height'] = width
    else:
      kw2['width'] = width
      kw2['height'] = defaultSize

    self.canvas = c = Canvas(self, bg=bg, **kw2)
    # grid does not seem to work because once canvas size specified
    # in any way it never seems to be resized
    c.pack(side=Tkinter.TOP, expand=Tkinter.YES, fill=Tkinter.BOTH)

    self.tlOuter = c.create_polygon(0, 0, 0, 0, 0, 0,0, 0,fill=color1, outline='')
    self.brOuter = c.create_polygon(0, 0, 0, 0, 0, 0,0, 0,fill=color2, outline='')
    self.bgOuter = c.create_rectangle(0, 0, 0, 0, fill=bg_trough, outline='')


    # self.arrowXTl,Br are for shadows of the arrows
    # self.arrowXFill is for triangle infill
    # arrow and rectangle boundaries which are colored the
    # same as (and so in some sense are part of) the trough
    # the main part of the arrow does not need its own polygon
    # because it is the same color as the main rectangle
    self.arrow1Tl = c.create_polygon(0, 0, 0, 0, 0, 0, 0, 0,
                                    fill=color3, outline='')
    self.arrow1Br = c.create_polygon(0, 0, 0, 0, 0, 0, 0, 0,
                                    fill=color4, outline='')
    self.arrow1Fill = c.create_polygon(0, 0, 0, 0, 0, 0,
                                    fill=bg, outline='')
                                    
    self.arrow2Tl = c.create_polygon(0, 0, 0, 0, 0, 0, 0, 0,
                                    fill=color3, outline='')
    self.arrow2Br = c.create_polygon(0, 0, 0, 0, 0, 0, 0, 0,
                                    fill=color4, outline='')
    self.arrow2Fill = c.create_polygon(0, 0, 0, 0, 0, 0,
                                    fill=bg, outline='')

    # slider is the bit you grab in the middle
    # created with two triangles for the shadowsed border
    # with rectangle panel in middle
    
    self.tlSlider = c.create_polygon(0, 0, 0, 0, 0, 0,fill=color3, outline='')
    self.brSlider = c.create_polygon(0, 0, 0, 0, 0, 0,fill=color4, outline='')
    self.bgSlider = c.create_rectangle(0, 0, 0, 0, fill=bg, outline='')
    

    if (show_text):
      self.text0 = c.create_text(0, 0, fill=text_color)
      self.text1 = c.create_text(0, 0, fill=text_color)
      if (orient == Tkinter.HORIZONTAL):
        c.itemconfig(self.text0, anchor=Tkinter.W)
        c.itemconfig(self.text1, anchor=Tkinter.E)
      else:
        c.itemconfig(self.text0, anchor=Tkinter.S)
        c.itemconfig(self.text1, anchor=Tkinter.N)

    c.bind('<Configure>', self.configure)
    c.bind('<Button-1>', self.buttonPressMove)
    c.bind('<B1-Motion>', self.buttonMotionMove)
    c.bind('<ButtonRelease-1>', self.buttonReleaseMove) 

    if (allow_resize):
      # c.bind('<Button-2>', self.buttonPressResize)
      # c.bind('<B2-Motion>', self.resize)
      c.bind(buttonClick(2), self.buttonPressResize)
      c.bind(buttonMotion(2), self.resize)

  def updateScrollbar(self):

    self.updateText()
    self.updateBubble()
 
  def setArrow1State(self, isNormal):
 
    bg_light = self.tlColor
    bg_dark  = self.brColor
    itemconfig = self.canvas.itemconfig
    if isNormal:
      itemconfig(self.arrow1Tl, fill=bg_light)
      itemconfig(self.arrow1Br, fill=bg_dark )
    else:
      itemconfig(self.arrow1Tl, fill=bg_dark )
      itemconfig(self.arrow1Br, fill=bg_light)

  def setArrow2State(self, isNormal):
 
    bg_light = self.tlColor
    bg_dark  = self.brColor
    itemconfig = self.canvas.itemconfig
    if isNormal:
      itemconfig(self.arrow2Tl, fill=bg_light)
      itemconfig(self.arrow2Br, fill=bg_dark)
    else:
      itemconfig(self.arrow2Tl, fill=bg_dark)
      itemconfig(self.arrow2Br, fill=bg_light)

  def setBubbleState(self, isNormal):
 
    tlColor = self.tlColor
    brColor = self.brColor
    itemconfig = self.canvas.itemconfig
    
    if isNormal:
      itemconfig(self.tlSlider, fill=tlColor)
      itemconfig(self.brSlider, fill=brColor)
    else:
      itemconfig(self.tlSlider, fill=brColor)
      itemconfig(self.brSlider, fill=tlColor)

  def getLocation(self, event):

    if (self.orient == Tkinter.HORIZONTAL):
      h = self.canvas.winfo_height()
      s = event.x - h
    else:
      w = self.canvas.winfo_width()
      s = event.y - w
 
    return s

  def getLength(self):

    w = self.canvas.winfo_width()
    h = self.canvas.winfo_height()
    if (self.orient == Tkinter.HORIZONTAL):
      d = w - 2*h
    else:
      d = h - 2*w

    d = max(1, d) # w = h = 1 early on

    return d

  def getBubble(self):

    d = self.getLength()
    s0 = int(self.lo * d)
    s1 = int(self.hi * d)

    if self.relief == 'flat':
      s0 += 1
      s1 -= 1

    dims = [self.canvas.winfo_width(),self.canvas.winfo_height()]

    # do not let bubble get too small
    limit = max(dims)
    limit -= 2 * min(dims)
    size = s1 - s0
    if size < min_bubble_size:
      s1 = s0 + min_bubble_size
      
      if s1 > limit:
        s1 = limit
        s0 = s1-min_bubble_size
      
    return (s0, s1)

  def updateBubble(self):
 
    w = int(self.canvas.winfo_width())
    h = int(self.canvas.winfo_height())
    b = self.borderwidth

    (s0, s1) = self.getBubble()
 
    if (self.orient == Tkinter.HORIZONTAL):
      s0 = s0 + h
      s1 = s1 + h
      limit = w/3
    else:
      s0 = s0 + w
      s1 = s1 + w
      limit = h/3
 
    #print 'updateBubble1', s0, s1, w, h, b
 
    c = self.canvas
    coords = c.coords
    if (self.orient == Tkinter.HORIZONTAL):
      # top left border triangle
      coords(self.tlOuter, 0, 0, w, 0, 2, h-2, 0, h, )
      # bottom right border triangle
      coords(self.brOuter, 0, h, 2, h-2, w, 0, w, h)
    
    else:
      # top left border triangle
      coords(self.tlOuter, 0, 0, w, 0, w-2, 2, 0, h)
      # bottom right border triangle
      coords(self.brOuter, w, h, w, 0, w-2, 2, 0, h)

    # back panel (trough)
    coords(self.bgOuter, b, b, w-b, h-b)
    
    if (self.orient == Tkinter.HORIZONTAL):
      #w2 = min(limit, w)
    
      x0 = b
      x2 = h
      x1 = (x0+x2)/2
      x3 = w-h
      x5 = w-1
      x4 = (x3+x5)/2
      
      y0 = b-1
      y2 = h-1
      y1 = (y0+y2)/2
      
      # left arrow
          
      # top left
      coords(self.arrow1Tl, x0, y1, x2, y0, x1, y1, x0, y1)
      # bottom right
      coords(self.arrow1Br, x0, y1, x1, y1, x2, y0, x2, y2)
      # infill
      coords(self.arrow1Fill, x0+b, y1, x2-b, y0+b, x2-b, y2-b)
      
      # right arrow
      
      # top left
      coords(self.arrow2Tl, x3, y0, x5, y1, x4, y1, x3, y2)
      # bottom right
      coords(self.arrow2Br, x4, y1, x5, y1, x3, y2, x4, y1)
      # infill
      coords(self.arrow2Fill, x3+b, y0+b, x5-b, y1, x3+b, y2-b)
      
    else:
      #w2 = min(limit, w)
    
      x0 = b
      x2 = w-x0
      x1 = (x0+x2)/2
      y0 = b-1
      y2 = w
      y1 = (y0+y2)/2
      y3 = h-w
      y5 = h-b-1
      y4 = (y3+y5)/2
      
      # upper arrow
          
      # top left
      coords(self.arrow1Tl, x0, y2, x1, y0, x1, y1, x0, y2)
      # bottom right
      coords(self.arrow1Br, x1, y0, x2, y2, x0, y2, x1, y1)
      # infill
      coords(self.arrow1Fill, x0+b, y2-b, x1, y0+b, x2-b, y2-b)
      
      # lower arrow
      
      # top left
      coords(self.arrow2Tl, x0, y3, x2, y3, x1, y4, x1, y5)
      # bottom right
      coords(self.arrow2Br, x1, y4, x2, y3, x1, y5, x1, y4, )
      # infill
      coords(self.arrow2Fill, x0+b, y3+b, x2-b, y3+b, x1, y5-b)
    
    if self.orient == Tkinter.VERTICAL:
      x0 = b
      y0 = s0
      x1 = w-x0
      y1 = s1
    else:
      x0 = s0
      y0 = b
      x1 = s1
      y1 = h-y0
  
    coords(self.tlSlider, x0, y0, x1, y0, x0, y1)
    coords(self.brSlider, x0, y1-1, x1-1, y0, x1, y1)
    coords(self.bgSlider, x0+b, y0+b, x1-b, y1-b)
    
  def updateText(self):
 
    if (not self.show_text):
      return

    # TBD: no guarantee for now that text fits inside canvas
 
    (t0, t1) = self.text_func()
    #print 'updateText', t0, t1
    self.canvas.itemconfig(self.text0, text=t0)
    self.canvas.itemconfig(self.text1, text=t1)
 
    w = self.canvas.winfo_width()
    h = self.canvas.winfo_height()
    if (self.orient == Tkinter.HORIZONTAL):
      x = h + 1
      y = h / 2
      self.canvas.coords(self.text0, x, y)
      x = w - h - 1
      self.canvas.coords(self.text1, x, y)
    else:
      x = w / 2
      y = h - w - 1
      self.canvas.coords(self.text0, x, y)
      y = w + 1
      self.canvas.coords(self.text1, x, y)

  def configure(self, event):
 
    if (self.first_pass):
      self.first_pass = False
      self.update_idletasks()

    #self.canvas.config(width=event.width)
    #self.canvas.config(height=event.height)
    self.updateScrollbar()

  def buttonPressMove(self, event):
 
    #print 'buttonPressMove1', event.state

    s = self.getLocation(event)
    d = self.getLength()

    if (s < 0 or s >= d):
      if (self.min_thickness):
        w = self.hi - self.lo
        if (w > 0):
          f = max(self.units_scroll, self.min_thickness/w)
        else:
          f = self.units_scroll
      else:
        f = self.units_scroll

    (s0, s1) = self.getBubble()
    self.start = s

    if (s < 0): # top (left) arrow click
      self.mode = top_left_arrow_mode # top (left) arrow click
      self.setArrow1State(False)
      self.fraction = -f
    elif (s < s0):
      self.mode = top_left_trough_mode # top (left) trough click
      self.fraction = -self.pages_scroll
    elif (s < s1):
      self.setBubbleState(False)
      self.mode = bubble_mode # bubble click
    elif (s < d):
      self.mode = bottom_right_trough_mode # bottom (right) trough click
      self.fraction = self.pages_scroll
    else:
      self.setArrow2State(False)
      self.mode = bottom_right_arrow_mode # bottom (right) arrow click
      self.fraction = f
 
    if (self.mode != bubble_mode):
      self.continueMove = True
      self.delay = self.repeatdelay
      self.moveFraction()
 
    #print 'buttonPressMove2', self.mode

  def buttonMotionMove(self, event):
 
    if (self.mode == bubble_mode):
      self.moveContinuous(event)

  def buttonReleaseMove(self, event):
 
    #print 'buttonReleaseMove', event.state, self.mode

    self.continueMove = False

    if (self.mode == top_left_arrow_mode):
      self.setArrow1State(True)
    elif (self.mode == bubble_mode):
      self.setBubbleState(True)
    elif (self.mode == bottom_right_arrow_mode):
      self.setArrow2State(True)

  def moveFraction(self):

    if (not self.continueMove):
      return

    d = self.fraction * (self.hi - self.lo)
    #print 'moveFraction', self.lo, self.hi, self.fraction, self.lo+d, self.hi+d
    self.set(self.lo+d, self.hi+d, do_callback=True)

    if (self.continueMove): # should be True
      (s0, s1) = self.getBubble()
      if (self.mode < bubble_mode):
        if (self.start >= s0):
          return
      else: # (self.mode > bubble_mode)
        if (self.start <= s1):
          return
      self.after(self.delay, self.moveFraction)
      self.delay = self.repeatinterval

  def moveContinuous(self, event):
 
    s = self.getLocation(event)
    d = self.getLength()
    ds = s - self.start
    r = float(ds) / d
    self.start = s
    self.set(self.lo+r, self.hi+r, do_callback=True)
 
  def buttonPressResize(self, event):
 
    #print 'buttonPressResize', event.state

    s = self.getLocation(event)
    #print 'buttonPressResize', event.x, event.y, s, self.lo, self.hi
    d = self.getLength()
    r = float(s) / d
 
    if (r < 0.5*(self.lo+self.hi)):
      self.mode = 0
    else:
      self.mode = 1

    if event.state & shift_key_state:
      self.haveShift = True
    else:
      self.haveShift = False
 
    self.start = s
 
  def resize(self, event):
 
    s = self.getLocation(event)
    d = self.getLength()
    ds = s - self.start
    r = float(ds) / d
 
    lo = self.lo
    hi = self.hi
    if self.haveShift:
      d = 0.25*(hi-lo)
    else:
      d = 0.5*(hi-lo)
    #print 'resize', event.x, event.y, self.mode, s, ds, lo, hi, r, d
    # protect against hi < lo
    if (self.mode == 0):
      if (d < min_delta and r > 0):
        return
      r = min(r, d)
      lo = max(0.0, lo + r)
      if self.haveShift:
        hi = min(1.0, hi - r)
    else:
      if (d < min_delta and r < 0):
        return
      r = -r
      r = min(r, d)
      hi = min(1.0, hi - r)
      if self.haveShift:
        lo = max(0.0, lo + r)
 
    self.start = s
    self.set(lo, hi, do_callback=True)
    #print 'resize2', lo, hi

  def get(self):

    return (self.lo, self.hi)

  def set(self, lo, hi, do_callback=False):

    assert lo < hi and (hi-lo) <= 1.0, 'lo = %s, hi = %s' % (lo, hi)

    #print 'set1', lo, hi, self.lo, self.hi

    if (lo < 0):
      hi = hi - lo
      lo = 0.0
    elif (hi > 1.0):
      lo = lo - hi + 1.0
      hi = 1.0

    #print 'set2', lo, hi, self.lo, self.hi

    if ((abs(lo - self.lo) < min_delta) and \
        (abs(hi - self.hi) < min_delta)):
      return

    if (self.min_thickness and (hi-lo) < 0.1*self.min_thickness):
      # here have new small thickness, but only do something if existing position not minimal
      # or if expanding rather than contracting region (can happen that min_thickness jumps up)
      #if ((self.hi - self.lo) > (self.min_thickness + min_delta)):
      if (((self.hi - self.lo) > (self.min_thickness + min_delta))
          or ((hi - lo) > (self.hi - self.lo + min_delta))):
        a = 0.5 * (lo + hi)
        b = 0.5 * self.min_thickness
        #print 'set3', lo, hi, self.min_thickness, a-b, a+b
        lo = a - b
        hi = a + b
        if (lo < 0):
          hi = hi - lo
          lo = 0
        elif (hi > 1):
          lo = lo - hi + 1
          hi = 1
      else:
        return
 
    #print 'set4', lo, hi
    self.lo = lo
    self.hi = hi

    self.updateScrollbar()

    if do_callback and self.callback:
      self.callback(lo, hi)

  def moveScrollbar(self, mode):

    if (mode == top_left_arrow_mode): # top (left) arrow click
      fraction = -self.units_scroll
    elif (mode == top_left_trough_mode): # top (left) trough click
      fraction = -self.pages_scroll
    elif (mode == bubble_mode): # bubble click
      pass
    elif (mode == bottom_right_trough_mode): # bottom (right) trough click
      fraction = self.pages_scroll
    elif (mode == bottom_right_arrow_mode): # bottom (right) arrow click
      fraction = self.units_scroll
    else:
      pass

    d = fraction * (self.hi - self.lo)
    self.set(self.lo+d, self.hi+d, do_callback=True)

if __name__ == '__main__':
 
  from memops.universal.Util import formatDecimals

  s = None

  def callback(lo, hi):
    print 'callback', lo, hi
 
  def getText():

    (lo, hi) = s.get()
    t0 = formatDecimals(lo, decimals=3)
    t1 = formatDecimals(hi, decimals=3)

    return (t0, t1)

  root = Tkinter.Tk()
  root.grid_rowconfigure(0, weight=1)
  root.grid_columnconfigure(0, weight=1)
 
  s = Scrollbar(root,
                #orient=Tkinter.VERTICAL,
                show_text = True, text_func = getText,
                width = 30,
                borderwidth = 1,
                allow_resize = True,
                callback=callback)
  s.set(0.2, 0.6)
  #s.pack(fill=Tkinter.BOTH, expand=Tkinter.YES)
  s.grid(row=0, column=0, sticky=Tkinter.NSEW)
  s.grid_rowconfigure(0, weight=1)
  s.grid_columnconfigure(0, weight=1)
 
  root.mainloop()


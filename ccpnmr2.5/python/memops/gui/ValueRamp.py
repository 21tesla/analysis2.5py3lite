
"""
======================COPYRIGHT/LICENSE START==========================

ValueRamp.py: <write function here>

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
import Tkinter


from memops.gui.Canvas import Canvas
from memops.gui.Frame import Frame

class ValueRamp(Frame):

  def __init__(self, parent, callback, speed = 1.5, delay = 50,
               fill = 'grey73', outline = 'grey90', hilight='grey90',
               width = 144, height = 24, docKey=None, tipText=None, *args, **kw):

    self.callback = callback
    self.speed   = speed
    self.delay   = delay
    self.fill    = fill
    self.outline = outline
    self.hilight = hilight

    Frame.__init__(self, parent=parent, docKey=docKey, tipText=tipText, createToolTip=True, *args, **kw)

    # if canvas width/height fixed then resizing has no effect, it seems
    #self.grid_rowconfigure(0, weight=1)
    #self.grid_columnconfigure(0, weight=1)
    font = '-schumacher-clean-bold-r-normal--14-140-75-75-c-70-iso646.1991-irv'
    self.canvas = c = Canvas(self, width=width, height=height, relief='flat', borderwidth=0)
    c.bind('<Configure>', self.resize)
    c.bind('<ButtonPress>', self.press)
    c.bind('<Motion>', self.motion)
    c.bind('<Enter>', self.mouseEnter)
    c.bind('<Leave>', self.mouseLeave)
    c.bind('<ButtonRelease>', self.release)
    self.bg      = c.create_rectangle(1, 1, 1, 1, fill='red', outline='grey50')
    self.boxes0   = []
    self.boxes1   = []
    self.boxes2   = []
    self.strikes  = []
    self.dashes   = []

    for n in range(24):
      box1    = c.create_rectangle(0, 0, 0, 0, fill=outline, outline=outline)
      box2    = c.create_rectangle(0, 0, 0, 0, fill=outline, outline=outline)
      box0    = c.create_rectangle(0, 0, 0, 0, fill=outline, outline=outline)
      self.boxes0.append(box0) 
      self.boxes1.append(box1) 
      self.boxes2.append(box2) 
     
    for n in range(8):
      dash    = c.create_line(0,0,0,0,fill='black')
      self.dashes.append( dash )

    for n in range(4):
      strike = c.create_line(0,0,0,0,fill='black')
      self.strikes.append( strike  )
     
    c.grid(row=0, column=0, sticky=Tkinter.NSEW)

    self.continue_callback = False
    self.multiplier = 1.0

  def mouseOver(self):
    
    self.canvas.config(background = 'grey50')
    self.canvas.itemconfigure(self.bg, outline='white')
    for n in range(24):
      self.canvas.itemconfigure(self.boxes0[n],fill=self.fill, outline=self.fill)
      self.canvas.itemconfigure(self.boxes1[n],fill=self.fill, outline=self.fill)
      self.canvas.itemconfigure(self.boxes2[n],fill=self.fill, outline=self.fill)

  def mouseEnter(self, *opt):

    self.mouseOver()
    event = opt[0]
    x = self.canvas.canvasx(event.x)
    y = self.canvas.canvasy(event.y)
    obj0 = self.canvas.find('closest',x,y)[0]
    if obj0 in self.boxes0:
      i = self.boxes0.index(obj0)
      obj1 = self.boxes1[i]
      obj2 = self.boxes2[i]
      self.canvas.itemconfigure(obj0,fill=self.hilight, outline=self.hilight)
      self.canvas.itemconfigure(obj1,outline='grey50')
      self.canvas.itemconfigure(obj2,outline='white')
  
  def mouseLeave(self, *opt):

   self.canvas.config(background ='white')
   self.canvas.itemconfigure(self.bg, outline='grey50')
   for n in range(24):
      self.canvas.itemconfigure(self.boxes0[n],fill=self.outline, outline=self.outline)
      self.canvas.itemconfigure(self.boxes1[n],fill=self.outline, outline=self.outline)
      self.canvas.itemconfigure(self.boxes2[n],fill=self.outline, outline=self.outline)

  def resize(self, event):

    q = (4,3,2,1,1,2,3,4)
    self.update_idletasks()
    w = self.canvas.winfo_width()-2
    h = self.canvas.winfo_height() 
    self.canvas.coords(self.bg, 1, 1, self.canvas.winfo_width()-1, h-1)

    for n in range(8):
      x1 = int((n*w)/8)
      x2 = int(((n+1)*w)/8)
      x3 = int(((n+0.5)*w)/8)
      dw = q[n]
      h2 = h/2
      self.canvas.coords(self.dashes[n], x3-dw, h2, x3+dw+1, h2)
      if n >= 4:
        self.canvas.coords(self.strikes[n-4], x3, h2-dw, x3, h2+dw+1)
 
    for n in range(24):
      x1 = 1+(n*w/24)
      x2 = 1+((n+1)*w/24)
      self.canvas.coords(self.boxes0[n],  x1+1,   2, x2-2, h-3 )
      self.canvas.coords(self.boxes1[n],  x1+0,   1, x2-1, h-2 )
      self.canvas.coords(self.boxes2[n],  x1+1,   2, x2-1, h-2 )

  def calcMultiplier(self, x):

    w = self.canvas.winfo_width()
    r = float(x) / w
    r = 2*r - 1 # puts in range -1 to 1
    self.multiplier = pow(self.speed, r)

  def press(self, event):

    self.calcMultiplier(event.x)
    self.continue_callback = True
    self.after(self.delay, self.doCallback)
 
  def motion(self, event):

    self.mouseOver()
    if (self.continue_callback):
      self.calcMultiplier(event.x)
    x = self.canvas.canvasx(event.x)
    y = self.canvas.canvasy(event.y)
    obj0 = self.canvas.find('closest',x,y)[0]
    if obj0 in self.boxes0:
      i = self.boxes0.index(obj0)
      obj1 = self.boxes1[i]
      obj2 = self.boxes2[i]
      self.canvas.itemconfigure(obj0,fill=self.hilight, outline=self.hilight)
      self.canvas.itemconfigure(obj1,outline='grey50')
      self.canvas.itemconfigure(obj2,outline='white')

  def release(self, event):

    self.continue_callback = False

  def doCallback(self):

    self.update_idletasks()

    if (self.continue_callback):
      if (self.multiplier != 1):
        self.callback(self.multiplier)
      self.after(self.delay, self.doCallback)

if (__name__ == '__main__'):

  def callback(multiplier):
    print 'callback', multiplier

  r = Tkinter.Tk()

  f = ValueRamp(r, callback=callback, tipText='Click to increase or decrease value')
  f.pack()

  r.mainloop()

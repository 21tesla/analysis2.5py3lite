
"""
======================COPYRIGHT/LICENSE START==========================

TabbedeFrame.py: <write function here>

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

# Creates a system of embedded frames, each with a selectable toggle tab at the top

import Tkinter

from memops.gui.Color import hexRepr
from memops.gui.Canvas import Canvas
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.ToolTip import ToolTip
from memops.universal.Util import buttonClick

triangleOutline='#602000'
triangleFill='#B05848'

h0  = 30 # Top panel height
bw  = 2  # Inner border width
th  = 1  # Bead inner width
de  = 4  # Width to inner edge of bead
bw2 = 2*bw
pad = 2  # Top exclusion

class TabbedFrame(Frame):

  def __init__(self, parent, options, selected=0, toggleOff=False, relief='raised',
               callback=None, font='Helvetica 10', docKeys=None, tipTexts=None, *args, **kw):
     
    Frame.__init__(self, parent, *args, **kw)
        
    self.selected = selected
    self.options  = options
    self.tipTexts = tipTexts
    self.docKeys  = docKeys
    self.font       = font
    self.callback   = callback
    self.bg         = self.cget('bg')
    self.toggleOff  = toggleOff
    self.inactColor = '#D0C8C0'
    self._frameCache = []
    self.sideFrame = None
    
    self.toolTip = ToolTip(self, text='')
    
    if relief == 'sunken':
      self.edge1Color ='grey60'
      self.inlayColor ='#D0B0A0'
      self.edge2Color ='white'
    else:
      self.edge1Color ='white'
      self.inlayColor ='#D0B0A0'
      self.edge2Color ='grey60'
    
    self.canvas = Tkinter.Canvas(self, highlightbackground=self.bg,
                                 background=self.bg, width=1, height=1)
    self.canvas.grid()
    
    self.setTabs(options, selected, tipTexts, docKeys)
    
    self.canvas.bind('<Button-1>', self.mouseClick)

    # self.canvas.bind('<Button-2>', self.mouseClick)
    # self.canvas.bind('<Button-3>', self.mouseClick)
    self.canvas.bind(buttonClick(2), self.mouseClick)
    self.canvas.bind(buttonClick(3), self.mouseClick)

    self.bind('<Configure>', self.redraw)
  
  def setTabs(self, options, selected=None, tipTexts=None, docKeys=None):
  
    self.selected = selected or self.selected
    self.options  = options
    nOptions = len(options)
    self.numTabs  = len(options)
     
    self.texts  = []
    self.frames = []
    self.boxes  = []
    self.markers = []
    self.optionDict = {}
    
    if self.selected >= nOptions:
      self.selected = 0
    
    if not docKeys:
      docKeys = options[:]
      
    while len(docKeys) < nOptions:
      docKeys.append(None)
    
    tipTexts = tipTexts or self.tipTexts or [None] * nOptions
    while len(tipTexts) < nOptions:
      tipTexts.append(None)
        
    rect = self.canvas.create_rectangle
    poly = self.canvas.create_polygon
    ctext = self.canvas.create_text
    tagBind = self.canvas.tag_bind
    
    self.canvas.delete('all')
    for frame in self._frameCache:
      frame.grid_forget()
      
    if self.sideFrame:  
      self.sideFrame.grid_forget()
    
    self.borders = []
    for i in range(5):
      b1 = rect(-1,-1,-1,-1,width=0,fill=self.edge1Color)
      b2 = rect(-1,-1,-1,-1,width=0,fill=self.edge2Color)
      b3 = rect(-1,-1,-1,-1,width=0,fill=self.inlayColor)
      self.borders.append((b1,b2,b3))
         
    for i in range(self.numTabs):
      box1 = rect(-1,-1,-1,-1,width=0,fill=self.edge1Color)
      box2 = rect(-1,-1,-1,-1,width=0,fill='')
      box3 = rect(-1,-1,-1,-1,width=0,fill='')
      mark = poly(-1,-1,-1,-1,-1,-1,fill=triangleFill,outline=triangleOutline)
      text = ctext(0,0,text=self.options[i])
      
      self.boxes.append((box1,box2,box3))
      self.markers.append(mark)
      self.texts.append(text)
            
      for item in (box1,box2,box3,text,mark):
        tagBind(item, '<Enter>', self.mouseOver, add='+')
        tagBind(item, '<Leave>', self.mouseLeave, add='+')
        self.optionDict[item] = i
      
      frame = Frame(self,background=self.bg, docKey=docKeys[i],
                    tipText=tipTexts[i], createToolTip=False)
 
                    
      self.frames.append(frame)
                      
    bbox = self.canvas.bbox(text)
    self.minHeight = (bbox[3]-bbox[1]) + de+de+bw2
    
    self.sideFrame = Tkinter.Frame(self,background=self.bg)
    self.sideFrame.grid_columnconfigure(0, weight=1)
    self._frameCache = self.frames
    
    self.redraw()
  
  def mouseLeave(self, event):
  
    self.toolTip.deactivate()
  
  def mouseOver(self, event):
  
    self.toolTip.deactivate()
  
    x = self.canvas.canvasx(event.x)
    y = self.canvas.canvasy(event.y)
    item = self.canvas.find('closest',x,y)[0]
    index = self.optionDict.get(item)
    
    if index is not None:
      tipText = self.frames[index].tipText
         
      if tipText:
        self.toolTip.text = tipText
        self.toolTip.activate()
  
  def setFont(self, font):
  
    self.font = font
    self.redraw()
        
  def redraw(self, event=None):
   
    if event:
      w = event.width
      h = event.height
      self.config(width=w, height=h)
    else:
      w = int(self.winfo_width())
      h = int(self.winfo_height())

    self.after_idle(lambda: self.draw(w, h))
     
  def draw(self, w, h):
    
    if self.selected is None:
      minH = h0+th+2+pad
    else:
      minH = h0+th+2+pad + self.minHeight
    
    h = max(minH,h)
    w = max(de,w)
    
    c  = self.canvas
    #c.configure(width=w, height=h)
    cc = c.coords
    ci = c.itemconfigure
    bb = c.bbox
    cl = c.lift
    
    font = self.font
    options = self.options

    c.place(x=0, y=0, anchor='nw', width=w, height=h+2, bordermode='ignore')
    
    numChars = 0.0
    for option in self.options:
      numChars += float(len(option))

    
    l4 = de
    l3 = l4-1
    l2 = l3-th
    l1 = l2-1
    
    t1 = h0
    t2 = t1+1
    t3 = t2+th
    t4 = t3+1
    
    r1 = w-de
    r2 = r1+1
    r3 = r2+th
    r4 = r3+1
    
    b1 = h-de
    b2 = b1+1
    b3 = b2+th
    b4 = b3+1
    
    selectedX = None
    
    x = l4
    for i in range(self.numTabs):
      frame = self.frames[i]
      text  = self.texts[i]
      mark  = self.markers[i]
     
      ci(text, text=options[i], font=font)
      bbox = bb(text)
      bwidth = bbox[2]-bbox[0]
      
      width = min(bwidth+16, int((r1-l4)*len(options[i])/numChars)+16 )
      
      x2    = x+width
      x1    = x2+1
      xm    = x+(width/2)
      cc(text, xm, 0.6*h0)
      box1,box2,box3 = self.boxes[i]
      
      chars = options[i]
      while chars and (bwidth >= width):
        chars = chars[:-1]
        ci(text, text=chars)
        bbox = bb(text)
        bwidth = bbox[2]-bbox[0]
            
      if i == self.selected:
        selectedX = (x,x2)
        xf = bw+de
        yf = h0+bw+th+2
        frame.place(x=xf, y=yf, anchor='nw', width=w-(xf+xf),
                    height=max(0,h-(yf+xf)), bordermode='ignore')
        ci(text, fill='black')
        cc(mark, x+3, pad+3, x+11, pad+3, x+3, pad+11)
        ci(mark, fill=triangleFill,outline=triangleOutline)
        ci(box1, fill=self.edge1Color)
        ci(box2, fill=self.edge2Color)
        ci(box3, fill=self.inlayColor)
        cc(box1, x,   pad,   x1,   t4)
        cc(box2, x+1, pad+1, x1,   t4)
        cc(box3, x+1, pad+1, x1-1, t3)
        
      else:
        frame.place_forget()
        ci(text, fill='grey40')
        cc(mark, x+3, pad+3, x+11, pad+3, x+3, pad+11)
        ci(mark, fill='grey97',outline='grey65')
        ci(box1, fill=self.edge2Color)
        ci(box2, fill=self.edge2Color)
        ci(box3, fill=self.inactColor)
        cc(box1, x,   pad+1, x1,   t4)
        cc(box2, x+1, pad+2, x1,   t4)
        cc(box3, x+1, pad+2, x1-1, t3)

      x = x2
    
    x += bw
    self.sideFrame.place(x=x, y=0, anchor='nw', bordermode='ignore',
                         width=w-x, height=t1)
    x += int(self.sideFrame.winfo_width())
    
    if self.selected is None:
      for i in range(5):
        e1, e2, e3 = self.borders[i]
        cc(e1, -1, -1, -1, -1)
        cc(e2, -1, -1, -1, -1)
        cc(e3, -1, -1, -1, -1)
      self.config(height=minH)
      return  
    
    else:
      self.config(height=b4)
    
    e1, e2, e3 = self.borders[0]
    cc(e1, l1, t1, l4, b4)
    cc(e2, l2, t2, l4, b4)
    cc(e3, l2, t2, l3, b3)

    e1, e2, e3 = self.borders[1]
    cc(e1, l4, b1, r4, b4)
    cc(e2, l3, b2, r4, b4)
    cc(e3, l3, b2, r3, b3)
    
    e1, e2, e3 = self.borders[2]
    cc(e1, r1, t1, r4, b2)
    cc(e2, r2, t2, r4, b2)
    cc(e3, r2, t2, r3, b2)

    if selectedX:
      x1, x2 = selectedX
      e1, e2, e3 = self.borders[3]
      cc(e1, l3, t1, x1+1, t4)
      cc(e2, l3, t2, x1+3, t4)
      cc(e3, l3, t2, x1+3, t3)
      cl(e1); cl(e2); cl(e3)
      
      e1, e2, e3 = self.borders[4]
      cc(e1, x2,   t1, r2, t4)
      cc(e2, x2-1, t2, r2, t4)
      cc(e3, x2-2, t2, r2, t3)
      cl(e1); cl(e2); cl(e3)
         
    else:
      for i in (3,4):
        e1, e2, e3 = self.borders[i]
        cc(e1, l3, t1, r2, t4)
        cc(e2, l3, t2, r2, t4)
        cc(e3, l3, t2, r2, t3)
        cl(e1); cl(e2); cl(e3)

  def mouseClick(self, event):

    self.configure(cursor="watch")
    self.update_idletasks()
    self.toolTip.deactivate()
    
    x = self.canvas.canvasx(event.x)
    y = self.canvas.canvasy(event.y)
    obj = self.canvas.find('closest',x,y)[0]
    
    selected = self.optionDict.get(obj)

    if selected is not None:
      if self.toggleOff and (selected == self.selected):
        self.selected = None
        if self.callback:
          self.callback(None)
 
      else:
        if selected is not self.selected:
          self.selected = selected
          if self.callback:
            self.callback(selected)
          

    self.redraw()
    
    self.configure(cursor="")

  def select(self, n):
  
    if self.toggleOff and (n == self.selected):
      self.selected = None
      if self.callback:
        self.callback(n)
 
    else:
      if n is not self.selected:
        self.selected = n
        if self.callback:
          self.callback(n)

    self.redraw()
    

if __name__ == '__main__':

  from ScrolledMatrix import ScrolledMatrix
  from LabelFrame import LabelFrame
  from Button import Button

  root = Tkinter.Tk()
  
  root.grid_columnconfigure(0, weight=1)
  root.grid_rowconfigure(1, weight=1)
  root.grid_rowconfigure(2, weight=2)
  
  label = Label(root, text='Label outside tabbed frame')
  label.grid(row=0, column=0, sticky='ew')
  
  tabbedFrame = TabbedFrame(root, options=['Option A','Option B','Option C','Option D'],
                            toggleOff=True, selected=2, tipTexts=['Tip A','Tip B','Tip C','Tip D'])
  tabbedFrame.grid(row=1, column=0, sticky='nsew')
  
  frames = tabbedFrame.frames

  l1 = Label(frames[0], text='Inside Frame1')
  l1.grid()

  l2 = Label(frames[1], text='Inside Frame2')
  l2.grid()
  
  frames[2].grid_columnconfigure(0, weight=1)
  frames[2].grid_rowconfigure(0, weight=1)
  l3 = ScrolledMatrix(frames[2], headingList=['Apple','Orange','Banana'])
  l3.grid(sticky='nsew')
  
  l4 = Label(frames[3], text='Inside Frame3')
  l4.grid()
  
  def changeTabs():
    tabbedFrame.setTabs(['Option 1','Option 2','Option 3'])
    
    frames = tabbedFrame.frames

    l1 = Label(frames[0], text='Totally new frame 1')
    l1.grid()

    l2 = Label(frames[1], text='Totally new frame 2')
    l2.grid()
   
    l3 = Label(frames[2], text='Totally new frame 3')
    l3.grid()
  
  button = Button(root, text='Change tabs', command=changeTabs, grid=(2,0), sticky='ew')
  
  frame2 = LabelFrame(root, text='Text Label')
  frame2.grid(row=3, column=0, sticky='nsew')
  frame2.grid_columnconfigure(0, weight=1)
  frame2.grid_rowconfigure(0, weight=1)

  l5 = ScrolledMatrix(frame2, headingList=['Peregrin','Meriadoc','Samwise'])
  l5.grid(row=0, column=0, sticky='nsew')
  #l5 = Label(frame2, text='LowerFrame')
  #l5.grid(row=0, column=0, sticky='nsew')

  root.mainloop()

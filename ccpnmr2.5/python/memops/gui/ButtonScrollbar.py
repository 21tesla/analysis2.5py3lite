
"""
======================COPYRIGHT/LICENSE START==========================

ButtonScrollbar.py: <write function here>

Copyright (C) 2006 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

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
from memops.gui.Color  import scaleColor
from memops.gui.Frame  import Frame
from memops.gui.Label  import Label
from memops.universal.Util import buttonClick, buttonRelease, buttonMotion

TOOL_TIP = """Left-click (& drag): Slide selection.
Middle-click (& drag): Toggle button on/off."""

class ButtonScrollbar(Frame):

  def __init__(self, parent,numButtons,state=None,orient='h',background=None,
               selectColor='#C0C0FF',size=15,callback=None,label=None,
               grid=None, *args,**kw):
   
    self.numButtons  = numButtons
    self.orient      = orient
    self.selectColor = selectColor
    self.size        = size
    self.callback    = callback
    self.state       = state or [False] * numButtons
    self.canvasDict  = {}
    self.canvasNums  = {}
    self.initialN    = None
    self.initialS    = None

    if (background is None):
      self.bg = parent.cget('bg')
    else:
      self.bg = background

    Frame.__init__(self, parent, grid=grid, *args, **kw)

    self.bgLight  = scaleColor(self, self.bg, 1.5)
    self.bgDark   = scaleColor(self, self.bg, 0.4)
    self.bgLight2 = scaleColor(self, self.bg, 1.2)
    self.bgDark2  = scaleColor(self, self.bg, 0.7)
    
    if label:
      self.label = Label(self, text=label, tipText=TOOL_TIP)
      self.label.grid(row=0, column=0, sticky=Tkinter.EW)
      if orient in (Tkinter.HORIZONTAL,'horizontal','h'):
        c_row = 0
        c_col = 1
      else:
        c_row = 1
        c_col = 0
    else:
      self.label =None

    if orient in (Tkinter.HORIZONTAL,'horizontal','h'):
      self.orient = 'h'
      self.canvas = c = Canvas(self, bg=self.bg, height=size)
      self.grid_columnconfigure(c_col, weight=1)
    else:
      self.orient = 'v'
      self.canvas = c = Canvas(self, bg=self.bg, width=size)
      self.grid_rowconfigure(c_row, weight=1)
   
    c.grid(row=c_row,column=c_col,sticky='nsew')
    
    c.bind('<Configure>',        self.resizeCanvas)
    c.bind('<Button-1>',        self.slideState)
    c.bind('<B1-Motion>',        self.slideState)
    c.bind('<ButtonRelease-1>', self.slideStateDone)

    # c.bind('<Button-2>',        self.clickState)
    # c.bind('<B2-Motion>',        self.clickState)
    # c.bind('<ButtonRelease-2>', self.clickStateDone)
    c.bind(buttonClick(2),        self.clickState)
    c.bind(buttonMotion(2),        self.clickState)
    c.bind(buttonRelease(2), self.clickStateDone)

  def resizeCanvas(self, event):

    self.updateButtons()

  def clickState(self, event):

    if not self.initialN:
      self.initialN = i = self.getButtonNum(event)
      state = not self.state[i]
      if (not state) and (self.state.count(True) == 1):
        state = True

      self.state[i] = state
      self.initialS = state
      self.updateButtons()
      
    else:
      i = self.getButtonNum(event)
      if i < 0:
        return
      if i >= self.numButtons:
        return
      
      if i < self.initialN:
        for j in range(i,self.initialN):
          if (not self.initialS) and (self.state.count(True) == 1):
            continue
          
          self.state[j] = self.initialS
        self.updateButtons()
            
      elif i > self.initialN:
        for j in range(self.initialN+1,i+1):
          if (not self.initialS) and (self.state.count(True) == 1):
            continue
            
          self.state[j] = self.initialS
        self.updateButtons()

    if self.callback:
      self.callback(self.getSelected())

  def slideState(self, event):

    n = self.getButtonNum(event)
    if self.initialN is None:
      if self.state[n]:
        self.initialN = n
      else:
        minDiff = None
        for i in range(self.numButtons):
          if self.state[i]:
            diff = i-n
            if (minDiff is None) or (abs(diff) < abs(minDiff)):
              minDiff = diff
              
        if minDiff:
          if minDiff < 0:
            self.state = [self.state[-1],] + self.state[:-1]
          else:
            self.state = self.state[1:] + [self.state[0],]
          self.updateButtons()
                 
    else:
      delta = self.initialN-n 
      self.state = self.state[delta:] + self.state[:delta]
      self.updateButtons()
      self.initialN = n

    if self.callback:
      self.callback(self.getSelected())

  def slideStateDone(self, event):
    
    self.initialN = None

  def clickStateDone(self, event):
    
    self.initialN = None
    self.initialS = None
   
  def extendSelection(self, event):

    #print "extend", event.x, event.y

    if self.callback:
      self.callback(self.getSelected())

  def getButtonNum(self, event):
    
    item = self.canvas.find('closest',event.x,event.y)[0]
    n = self.canvasNums[item]
       
    return n
 
  def setNumButtons(self, i):

    for j in range(i, self.numButtons):
      if self.canvasDict.get(j):
        hilight, lolight, body = self.canvasDict[j]
        self.canvas.delete(hilight)
        self.canvas.delete(lolight)
        self.canvas.delete(body)
        del self.canvasDict[j]
      
    self.numButtons = i
    while len(self.state) < i:
      self.state.append(False) 

    while i > len(self.state):
      self.state.pop(-1)

    self.updateButtons()

  def getState(self):
  
    return self.state

  def setState(self, booleans):

    for i in range(self.numButtons):
      if i < len(booleans):
        if booleans[i]:
          self.state[i] = True
        else:
          self.state[i] = False

    self.updateButtons()
  
  def drawButton(self, i, x1, y1, x2, y2):
    
    if self.state[i]:
      color1 = self.selectColor
      color2 = self.bgLight
      color3 = self.bgDark
    else:
      color1 = self.bg
      color2 = self.bgDark2
      color3 = self.bgLight2
    
    if self.canvasDict.get(i):
      hilight, lolight, body = self.canvasDict[i]
    
      self.canvas.coords(hilight, x1,y1,x2,y2)
      self.canvas.itemconfigure(hilight, outline=color2,fill=color2)
      self.canvas.coords(lolight, x1+1,y1+1,x2,y2)
      self.canvas.itemconfigure(lolight, outline=color3,fill=color3)
      self.canvas.coords(body, x1+1,y1+1,x2-1,y2-1)
      self.canvas.itemconfigure(body, outline=color1,fill=color1)
    
    else:
      hilight = self.canvas.create_rectangle(x1  ,y1  ,x2  ,y2  ,outline=color2,fill=color2)
      lolight = self.canvas.create_rectangle(x1+1,y1+1,x2  ,y2  ,outline=color3,fill=color3)
      body    = self.canvas.create_rectangle(x1+1,y1+1,x2-1,y2-1,outline=color1,fill=color1)

      self.canvasDict[i] = (hilight, lolight, body)
      self.canvasNums[hilight] = i
      self.canvasNums[lolight] = i
      self.canvasNums[body] = i
    
  def updateButtons(self):  

    if self.orient == 'h':
      y1 = 0
      y2 = self.size
      cw = self.canvas.winfo_width()
      w  = cw/self.numButtons
      x1 = 0
      sp = cw - (w * self.numButtons)

      for i in range(self.numButtons):
        x2 = x1+w-1
        if sp:
          x2 += 1
          sp -= 1
          
        self.drawButton(i, x1, y1, x2, y2)
        x1 = x2+1

    else:
      x1 = 0
      x2 = self.size
      ch = self.canvas.winfo_height() 
      h  = ch/self.numButtons
      y1 = 0
      spe = ch - (h * self.numButtons)

      for i in range(self.numButtons):
        y2 = y1+h-1
        if sp:
          y2 += 1
          sp -= 1
          
        self.drawButton(i, x1, y1, x2, y2)
        y1 = y2+1

  def getSelected(self):

    return [ n for n in range(len(self.state)) if self.state[n] ]

if __name__ == '__main__':
 
  def callback(selected):
    print 'callback', selected

  root = Tkinter.Tk()
  root.grid_rowconfigure(0, weight=1)
  root.grid_columnconfigure(0, weight=1)

  sb = ButtonScrollbar(root, 200, background='#C0C0C0',
                       callback=callback, label='label')
  sb.grid(row=0, column=0, sticky=Tkinter.EW)
  state = [0] * 200
  state[77] = 1
  sb.setState(state)
  
  #sb.setNumButtons(100)
  
  root.mainloop()

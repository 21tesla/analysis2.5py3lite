
"""
======================COPYRIGHT/LICENSE START==========================

ScrolledDensityMatrix.py: <write function here>

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
import Tkinter, re

from math import log
log10 = log(10)

from memops.gui.Canvas import Canvas
from memops.gui.Frame import Frame
from memops.gui.ScrolledCanvas import ScrolledCanvas
from memops.gui.Color import hexRepr

class ScrolledDensityMatrix(ScrolledCanvas):

  def __init__(self, parent,  matrix=None, boxSize=8, maxVal=None,
               title=None, xLabels=[], yLabels=[], borderWidth=0.1,
               canvasOrigin=None, labelAxes=True, barPlot=False,
               borderColor='grey', zoom=1.0, doLegend=True,
               xAxisLabel=None, yAxisLabel=None,
               *args, **kw):

    ScrolledCanvas.__init__(self, parent, width=250, **kw)
    self.normFont = ("Helvetica", 8)
    self.boldFont = ("Helvetica", 10, "bold")
    
    self.origin = canvasOrigin or (0,0)
    self.labelAxes = labelAxes
    self.border = None
    self.indent = 50
    self.matrix = matrix
    self.borderWidth = borderWidth
    self.borderColor = borderColor
    self.zoom = zoom
    self.boxSize = boxSize
    self.maxVal = maxVal
    self.minVal = 0
    if not self.maxVal:
      if self.matrix is not None:
        self.setMatrixMaxVal()
      else:
        self.maxVal = 1.0
 
    self.hilightedItems = []
    
    self.barPlot = barPlot
    self.title  = title
    self.legend = []
    self.miscItems = []
    self.canvasDict = {}
    self.boxes = []
    self.xGrid = []
    self.yGrid = []
    self.xLabels = xLabels or []
    self.yLabels = yLabels or []
    self.xAxisLabel = xAxisLabel
    self.yAxisLabel = yAxisLabel
    self.xText = None
    self.yText = None
    self.xBg   = None
    self.yBg   = None
    self.doLegend = doLegend
    self.draw()

    self.canvas.bind('<Button-1>', self.mouseClick)
    self.canvas.bind('<Double-1>', self.mouseDoubleClick)
    self.canvas.bind('<Motion>',   self.mouseEnter)
    self.canvas.bind('<Enter>',    self.mouseEnter)
    self.canvas.bind('<Leave>',    self.mouseLeave)

  def update(self, matrix):
  
    self.matrix = matrix
    self.setMatrixMaxVal()
    self.draw()
    self.update_idletasks()

  def mouseEnter(self, event):
   
    self.removeMenu()

    if self.matrix is None:
      return
  
    c =self.canvas
    c.focus()
    x = c.canvasx(event.x)
    y = c.canvasy(event.y)
    
    f = (self.boxSize*self.zoom)
  
    self.mouseLeave()
    items = c.find('overlapping',x-f,y-f,x+f,y+f)
    for item in items:
      if self.canvasDict.get(item):
        if self.canvasDict[item][0][:4] == 'grid':
          self.hilightedItems.append(item)
          c.itemconfigure(item, fill='#5050b0')
          
    item = c.find('closest',x,y)[0]
    if self.canvasDict.get(item) and (self.canvasDict[item][0] == 'box'):
      textX = self.canvasDict[item][1]
      textY = self.canvasDict[item][2]
      if self.canvasDict[item][1] < len(self.xLabels):
        textX = self.xLabels[self.canvasDict[item][1]]
      if self.canvasDict[item][2] < len(self.yLabels):
        textY = self.yLabels[self.canvasDict[item][2]]

      x1 = f+c.canvasx(self.origin[0])
      y1 = y
      x2 = x
      y2 = f+c.canvasy(self.origin[1])
    
      if self.xText:
        c.itemconfigure(self.xText,text=textX)
        c.itemconfigure(self.yText,text=textY)
        c.coords(self.xText,x2,y2)
        c.coords(self.yText,x1,y1)
        (x1,y1,x2,y2) = c.bbox(self.xText)
        c.coords(self.xBg,x1-4,y1-4,x2+4,y2+4)
        (x1,y1,x2,y2) = c.bbox(self.yText)
        c.coords(self.yBg,x1-4,y1-4,x2+4,y2+4)
        c.itemconfigure(self.xBg,fill='#a0a0ff',outline='#5050b0')
        c.itemconfigure(self.yBg,fill='#a0a0ff',outline='#5050b0')
        c.lift(self.xBg)
        c.lift(self.yBg)
        c.lift(self.xText)
        c.lift(self.yText)
        
      else:
        self.xText = c.create_text(x2,y2,font=self.normFont,text=textX, anchor='nw', fill='darkBlue')
        self.yText = c.create_text(x1,y1,font=self.normFont,text=textY, anchor='sw', fill='darkBlue')

        (x1,y1,x2,y2) = c.bbox(self.xText)
        self.xBg = c.create_rectangle(x1-4,y1-4,x2+4,y2+4,fill='#a0a0ff',outline='#5050b0',width=2)
        (x1,y1,x2,y2) = c.bbox(self.yText)
        self.yBg = c.create_rectangle(x1-4,y1-4,x2+4,y2+4,fill='#a0a0ff',outline='#5050b0',width=2)

        c.lift(self.xText)
        c.lift(self.yText)
          
  
  def mouseLeave(self, *event):
  
    if self.xText:
      self.canvas.itemconfigure(self.xText,text='')
      self.canvas.itemconfigure(self.yText,text='')
      self.canvas.itemconfigure(self.xBg,fill='white',outline='white')
      self.canvas.itemconfigure(self.yBg,fill='white',outline='white')
      self.canvas.coords(self.xBg,0,0,0,0)
      self.canvas.coords(self.yBg,0,0,0,0)
    
    for item in self.hilightedItems:
      self.canvas.itemconfigure(item, fill=self.borderColor)
      
    self.hilightedItems = []
  
  def mouseClick(self, event):
  
    x = self.canvas.canvasx(event.x)
    y = self.canvas.canvasy(event.y)
    item = self.canvas.find('closest',x,y)
    if self.canvasDict.get(item) and (self.canvasDict[item][0] == 'box'):
      print ">", self.canvasDict[item][2], self.canvasDict[item][1]
  
  def mouseDoubleClick(self, event):
    
    x = self.canvas.canvasx(event.x)
    y = self.canvas.canvasy(event.y)
    item = self.canvas.find('closest',x,y)
    if self.canvasDict.get(item) and (self.canvasDict[item][0] == 'box'):
      print ">", self.canvasDict[item][2], self.canvasDict[item][1]
    
  def draw(self):
   
    if self.matrix is None:
      return
    
    c = self.canvas
    canvasDict = self.canvasDict
    matrix     = self.matrix
    boxes      = self.boxes
    legend     = self.legend
    w = max(3,int(self.boxSize * self.zoom))
    h = max(3,int(self.boxSize * self.zoom))
    h2 = h/2.0
    w2 = w/2.0
    N = len(matrix)
    M = len(matrix[0])
    
    createRect = c.create_rectangle
    createText = c.create_text
    cCoords = c.coords
    cConfig = c.itemconfigure
    xGrid = self.xGrid
    yGrid = self.yGrid
    barPlot = self.barPlot
    
    for item in self.miscItems:
      c.delete(item)
      
    if self.border:
      cCoords(self.border, 0,0, (N*w)+(2*self.indent), (M*h)+(2*self.indent))
    else:
      self.border = createRect(0,0, self.origin[0]+(N*w)+(2*self.indent),
                               self.origin[1]+(M*h)+(2*self.indent),
                               outline='white', width=1)

    if self.title:
      midWay = self.origin[0] + ((N*w)+(2*self.indent))/2.0
      item = createText(midWay,self.origin[1]+self.indent-h,
                        font=self.boldFont,text=self.title,
                        anchor='s', fill='black')
      self.miscItems.append(item)                  
    
    if self.xAxisLabel:
      midWay = self.origin[0] + ((N*w)+(2*self.indent))/2.0
      item = createText(midWay, self.origin[1]+(M*h)+self.indent+h+h+h,
                        font=self.boldFont, text=self.xAxisLabel,
                        anchor='s', fill='black')
      self.miscItems.append(item)                  
      
    if self.yAxisLabel:
      item = createText(self.indent-h2,self.origin[1]+self.indent-h2,
                        font=self.boldFont, text=self.yAxisLabel,
                        anchor='s', fill='black')
      self.miscItems.append(item)                  
      
    if len(xGrid) > N+1:
      for i in range(N+1,len(xGrid)):
        c.delete(xGrid[i])
      
    if len(yGrid) > M+1:
      for j in range(M+1,len(yGrid)):
        c.delete(yGrid[j])

    x = self.indent + self.origin[0]
    y = self.indent + self.origin[1]
    for j in range(M+1):
      
      if j >= len(yGrid):
        yGrid.append(None)
    
      if yGrid[j]:
        c.coords(yGrid[j],x,y,x+(N*h),y)
      else:
        yGrid[j] = self.makeLine(x,y,x+(N*h),y)
        canvasDict[yGrid[j]][0] = 'gridY'
        canvasDict[yGrid[j]][1] = j
      y += h

    y = self.indent + self.origin[1]
    for i in range(N+1):

      if i >= len(xGrid):
        xGrid.append(None)
        
      if xGrid[i]:
        c.coords(xGrid[i],x,y,x,y+(M*w))
      else:
        xGrid[i] = self.makeLine(x,y,x,y+(M*w))
        canvasDict[xGrid[i]][0] = 'gridX'
        canvasDict[xGrid[i]][1]  = i
      x += w

    x = self.indent + self.origin[0]
    y = self.indent + self.origin[1]
    makeBox = self.makeBox
    setBoxCoords = self.setBoxCoords
    configureBox = self.configureBox
    deleteBox    = self.deleteBox
    maxVal       = self.maxVal
    labelAxes    = self.labelAxes
    yLabels      = self.yLabels
    xLabels      = self.xLabels
    
    if labelAxes:
      for j in range(M):
        yText = createText(x-w2, y+h2, font=self.normFont,
                           text=yLabels[j], anchor='e',
                           fill='black')
        self.miscItems.append(yText)
        y += h
    
    dv1   = (maxVal-self.minVal)/20.0
    order = int(round((log(dv1)/log10)))
    dv    = round(dv1,-order) + 10**order
    start = (round(maxVal/dv)) * dv
    end   = min(0.0, (round(self.minVal/dv)) * dv)
    maxDiff = max(abs(start),abs(end))
    
    for i in range(N):
      y = self.indent + self.origin[1]
      if i >= len(boxes):
        boxes.append([])
      
      for j in range(M):
        if j >= len(boxes[i]):
          boxes[i].append(None)
      
      
        if barPlot:
          v =  matrix[i][j]

          if v < 0:
            color = '#A00000'
            f = -min(1.0,v/maxDiff)
            yA = y
            yB = y+(f*h)
          elif v >= 0:
            color = '#0000A0'
            f = min(1.0,v/maxDiff)
            yB = y+h
            yA = yB - (f*h)

          item = boxes[i][j]
          if item is None:
            boxes[i][j] = makeBox(c, x, yA, i, j, width=w,
                                  height=yB-yA, fillColor=color)
          else:
            #setBoxCoords(c,i,j,x,y,w,h) # Substituted for speed
            cCoords(item, x, yA, x+w, yB)

            #configureBox(c,i,j,color) # Substituted for speed
            cConfig(item,fill=color)
            
        else: # Normal colour density
        
          #idx = 1-(min(matrix[i][j],maxVal)/float(maxVal))
          #color = hexRepr(idx,idx,idx)
          v =  matrix[i][j]
 
          if v < 0:
            idx   = -min(1.0,v/maxDiff)
            color = hexRepr(1-(0.5*idx),1-idx,1-idx)
 
          elif v > 0:
            idx   = min(1.0,v/maxDiff)
            color = hexRepr(1-idx,1-idx,1-(0.5*idx))
 
          else:
            color = hexRepr(1,1,1)
 
          item = boxes[i][j]
          if item is None:
            boxes[i][j] = makeBox(c, x, y, i, j, width=w,
                                  height=h, fillColor=color)
          else:
            #setBoxCoords(c,i,j,x,y,w,h) # Substituted for speed
            cCoords(item, x, y, x+w, y+h)
 
            #configureBox(c,i,j,color) # Substituted for speed
            cConfig(item,fill=color)

        y += h
        
      if labelAxes:
        xText = createText(x+w2, y+h2, font=self.normFont,
                           text=xLabels[i], anchor='n',
                           fill='black')
        self.miscItems.append(xText)
        
      x += w

      # cleanup end of rows
      if len(boxes[i]) > M:
        for j in range(M,len(boxes[i])):
          if boxes[i][j]:
            deleteBox(i,j)
        boxes[i] = boxes[i][0:M]
    
    # clean up unused rows  
    if len(boxes) > N:
      for i in range(N,len(boxes)):
        for j in range(len(boxes[i])):
          if boxes[i][j]:
            deleteBox(i,j)
      boxes = boxes[0:N]

    if not self.doLegend:
      return

    # Create legend
    x += w+2
    y  = self.indent + self.origin[1]

    borderColor = self.borderColor
    
    if -5 < order < 6:
      textFormat  = '%%.%df' % max(1,-order)
    else:
      textFormat = '%.2e'

    # TBD should have text bounding box here
    sz = max([w,h])
    B = int((start-end)/dv)+1
    v = start
    for i in range(B):
      v = start - (i*dv)
      text = '?'
    
      if barPlot:
        if v < 0:
          color = '#A00000'
          f = -min(1.0,v/maxDiff)
          yA = y
          yB = y+(f*sz)
          text  = textFormat % v
        elif v >= 0:
          color = '#0000A0'
          f = min(1.0,v/maxDiff)
          yB = y+sz
          yA = yB - (f*sz)
          text  = ' ' + textFormat % v
          
        if len(legend) < B:
          itemA = createRect(x,yA,x+sz,yB,fill=color, outline=borderColor, width=1)
          itemB = createText(x+sz+sz/2, y, font=self.normFont,
                             text=text, anchor='w', fill='black')
          itemC = createRect(x,yA,x+sz,yB, outline=borderColor, width=1)
          legend.append((itemA, itemB, itemC))
        #else:
 
        itemA, itemB, itemC = legend[i]
        cCoords(itemA,x,yA,x+sz,yB)
        cCoords(itemB,x+sz+sz/2,y+sz/2)
        cCoords(itemC,x,y,x+sz,y+sz)
        cConfig(itemA, fill=color, outline=borderColor)
        cConfig(itemB, text=text, font=self.normFont)
        cConfig(itemC, outline=borderColor)
      
      else:
    
        if v < 0:
          idx   = -v/maxDiff
          color = hexRepr(1-(0.5*idx),1-idx,1-idx)
          text  = textFormat % v
 
        elif v >= 0:
          idx   = v/maxDiff
          color = hexRepr(1-idx,1-idx,1-(0.5*idx))
          text  = ' ' + textFormat % v
 
        if len(legend) < B:
          itemA = createRect(x,y,x+sz,y+sz,fill=color,
                             outline=borderColor, width=1)
          itemB = createText(x+sz+sz/2, y, font=self.normFont,
                             text=text, anchor='w', fill='black')
          legend.append((itemA, itemB))
        #else:
 
        itemA, itemB = legend[i]
        cCoords(itemA,x,y,x+sz,y+sz)
        cCoords(itemB,x+sz+sz/2,y+sz/2)
        cConfig(itemA, fill=color, outline=borderColor)
        cConfig(itemB, text=text, font=self.normFont)

      y += sz
    
    while len(legend) > B:
      items = legend.pop()
      
      for item in items:
        c.delete(item)
    
    self.legend = legend    

  def makeLine(self,x1,y1,x2,y2):
  
    item = self.canvas.create_line(x1,y1,x2,y2,fill=self.borderColor,width=self.borderWidth)
    obj = ['line',0]
    self.canvasDict[item] = obj
    
    return item
      
  def setZoom(self, zoom=1.0):
    
    if zoom:
      self.zoom = max(zoom, 0.1)
      self.draw()
      self.refresh()
    
  def setMatrixMaxVal(self, val=None):
  
    matrix = self.matrix
    m = range(len(matrix))
    n = range(len(matrix[0]))
    
    if val is None:
      maxVal = matrix[0][0]
      minVal = matrix[0][0]
      for i in m:
        for j in n:
          if matrix[i][j] > maxVal:
            maxVal = matrix[i][j]
          if matrix[i][j] < minVal:
            minVal = matrix[i][j]
    else:
      maxVal = val
      minVal = 0
      
    if maxVal == 0:
      maxVal = 1
      minVal = 0
    
    if (-1*minVal) > maxVal:
      maxVal = -1*minVal 
      
    self.maxVal = maxVal
    self.minVal = minVal

  def makeBox(self, canvas, x, y, i, j, fillColor='white', width=8, height=8):
  
    item = canvas.create_rectangle(x, y, x+width, y+height,
                                   fill=fillColor, width=0,
                                   outline=fillColor)
    obj = ['box',i,j,x + (width/2),y + (height/2)]
    self.canvasDict[item] = obj
    canvas.lower(item)
    return item
  
  def setBoxCoords(self,canvas,i,j,x,y,width,height):

    item = self.boxes[i][j]
    x2 = x + width
    y2 = y + height
  
    canvas.coords(item, x, y, x2, y2)

  def configureBox(self,canvas,i,j,color,borderWidth=None,borderColor=None):
  
    item = self.boxes[i][j]
    borderWidth = borderWidth or self.borderWidth
    borderColor = borderColor or self.borderColor
    
    canvas.itemconfigure(item,fill=color,outline=borderColor,width=borderWidth)
  
  def callback(self):
  
    print "Callback", self

  def deleteBox(self,i,j):
  
    item = self.boxes[i][j]
    self.canvas.delete(item)
    if self.canvasDict.get(item):
      del self.canvasDict[item]
    
if __name__ == '__main__':

  import Tkinter

  root = Tkinter.Tk()
  root.grid_columnconfigure(0, weight=1)
  root.grid_rowconfigure(0, weight=1)

  xx = [x/10.0 for x in range(1,11)]
  yy = [x*x for x in xx]

  xLabels = ['%.2f' % x for x in xx]
  yLabels = ['%.2f' % y for y in yy]

  matrix = []
  
  for x in xx:
    matrix.append([])
    for y in yy:
      matrix[-1].append(x*y - 0.4)

  matrix = ScrolledDensityMatrix(root,  matrix=matrix, title='Test Plot', barPlot=False,
                                boxSize=24,  xLabels=xLabels, yLabels=yLabels,
                                xAxisLabel='Numbers', yAxisLabel='Squares')

  matrix.grid(sticky='nsew')
  
  root.mainloop()

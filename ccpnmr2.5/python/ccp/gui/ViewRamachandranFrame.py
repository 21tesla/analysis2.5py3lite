
"""
======================COPYRIGHT/LICENSE START==========================

ViewRamachandranFrame.py: 

Copyright (C) 2005 Tim Stevens, Nicole Cheung (University of Cambridge)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in license/LGPL.license
 
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

- email: ccpn@bioc.cam.ac.uk

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

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

"""
import math
import tkFont
import cPickle

from memops.gui.ScrolledCanvas  import ScrolledCanvas
from memops.gui.Label           import Label
from memops.gui.Color           import scaleColor
from memops.gui.MessageReporter import showWarning
from memops.universal.Io        import getTopDirectory, joinPath
from memops.universal.Util      import buttonMotion

UNICODE_PHI = u'\u03A6'
UNICODE_PSI = u'\u03A8'

class ViewRamachandranFrame(ScrolledCanvas):

  def __init__(self, parent, callback=None, rightCallback=None, width=420, height=440, bgColor='#FFFFFF',
               showCoords=True, xTicks=True, yTicks=True, xLabel=UNICODE_PHI, yLabel=UNICODE_PSI, titleText='',
               nullColor='#F0F0F0', enterCallback=None, leaveCallback=None, defaultPlot=True, *args, **kw):
    
    self.parent         = parent
    ScrolledCanvas.__init__(self, parent, resizeCallback=self.drawCanvas,
                            width=width, height=height )
    
    self.bitmapDir = joinPath(getTopDirectory(),'data','ccpnmr','ramachandranPlot')
    self.canvas.config(bg=bgColor)
    self.rightCallback  = rightCallback  # handle for change of chosen quadrant
    self.callback       = callback       # handle for estimation selection
    self.enterCallback  = enterCallback  # handle for mouse enter an estimation
    self.leaveCallback  = leaveCallback  # handle for mouse leave an estimation
    self.defaultPlot    = defaultPlot
    
    self.showCoords = showCoords
    self.xTicks = xTicks
    self.yTicks = yTicks
    self.xLabel = xLabel
    self.yLabel = yLabel
    
    self.matrix         = []    # bin matrix forming the plot
    self.canvasDict     = {}    # aaType : a list of color distribution
    self.intensityDict  = {}    # aaType : a list of color distribution
    self.currentObjDict = {}    # object : estimation circle on plot
    self.currentCirDict = {}    # circle on plot : object
    self.residueLabels  = {}
    self.currentObjects  = []
    self.selectedObjects = []  # a list of selected objects
    
    self.binSize        = 10        # 10 degree for each bin
    self.cirTag         = "est"     # arbitrary tag of estimation circles
    self.maxInt         = 255       # max number of color used
    self.selectColor    = "#00C000"   # color of a selected estimation
    self.unselectColor  = "#000000"  # color of a unselected estimation
    self.nullColor      = nullColor
    self.bgColor        = bgColor
    self.objectColors   = []
    self.titleText      = titleText
    self.title          = None
    self.phi            = None
    self.psi            = None
    self.sd1            = None
    self.sd2            = None
    self.phiPsiList     = []
    self.labels         = []
    self.crosshair      = []
    self.plotCoords     = []
    
    self.aminoAcid      = None
    
    self.binWidth    = 10 # Arbitrary deafult - reset dynamically
    self.binHeight   = 10 # Arbitrary deafult - reset dynamically
    self.cirRadius   = 4
    self.axisObjects = []
    self.zoom = 1.0
    
    # add handle function to mouse events on estimation objects
    self.canvas.bind("<Button-1>", self.circleClick)
    
    #if self.enterCallback:
    #  self.canvas.tag_bind(self.cirTag, "<Enter>", self.circleEnter, "+")
      
    #if self.leaveCallback:
    #  self.canvas.tag_bind(self.cirTag, "<Leave>", self.circleLeave, "+")

    # self.canvas.bind('<B2-Motion>', self.mouseMotion)
    self.canvas.bind(buttonMotion(2), self.mouseMotion)

    self.canvas.bind('<Button-4>', self.zoomIn)
    self.canvas.bind('<Button-5>', self.zoomOut)
    self.canvas.bind('<Motion>',   self.mouseEnter)
    self.canvas.bind('<Enter>',    self.mouseEnter)
    self.canvas.bind('<Leave>',    self.mouseLeave)
    
    # add handle function to double mouse click on Rama plot
    # # self.canvas.bind("<Button-3>", self.rightClick)
    #self.canvas.bind(buttonClick(3), self.rightClick)

    self.drawCanvas()
    
  def mouseEnter(self, event):

    if not self.plotCoords:
      return
    
    color = '#000000'
      
    x0,y0,x1,y1 = self.plotCoords

    canvas  = self.canvas
    cCoords = canvas.coords
    cLine   = canvas.create_line
    cText   = canvas.create_text
    
    crosshair = self.crosshair
    canvas.focus_force()
    self.removeMenu()
    
    x = canvas.canvasx(event.x)
    y = canvas.canvasy(event.y)
        
    if (x < x0) or (x > x1) or (y < y0) or (y > y1):
      self.mouseLeave()
    else:
      getAngle = self.getAngleFromPlot
      phi = getAngle(x, angleType='phi')
      psi = getAngle(y, angleType='psi')
      
      if self.showCoords:
        coordStr = '[ %3.3f, %3.3f ]' % (phi,psi)
	coordXtxt = x0+(0.9*(x1-x0))
        coordYtxt = y0+y1
      else:
        coordStr = ''
	coordXtxt = 0
        coordYtxt = 0

      if crosshair:
        (itemV,itemH,itemT) = crosshair
        cCoords(itemV,x,y0,x,y1)
        cCoords(itemH,x0,y,x1,y)
        cCoords(itemT,coordXtxt,coordYtxt)
        canvas.itemconfig(itemT,text=coordStr)
      else:
        itemV = cLine(x,y0,x,y1,fill=color,width=1)
        itemH = cLine(x0,y,x1,y,fill=color,width=1)
        itemT = cText(coordXtxt,coordYtxt, text=coordStr, fill=color)
        self.crosshair = (itemV,itemH,itemT)
      
  def mouseLeave(self, *event):

    if self.crosshair:
      for item in self.crosshair:
        self.canvas.delete(item)
        
      self.crosshair = []

  def configMenu(self):

    items = [
            { 'kind': 'command', 'label': 'Print to file', 'command' : self.printCanvas,
              'tipText': 'Save the plot graphics as a PostScript file' },
            { 'kind': 'command', 'label': 'Reset Zoom',    'command' : self.resetZoom,
              'tipText': 'Go back to the original plot size' },
            ]
      
    self.menu.setMenuItems(items)

  def mouseMotion(self, event):


    if event.state & 1:
      # shift - zoom
      self.removeMenu()
 
      dx = event.x - self.initialX
      dy = self.initialY - event.y
 
      if dx == 0 and dy ==0:
        return
 
      if dy < 0:
        self.zoomOut()
      else:
        self.zoomIn()
    
    else:
      ScrolledCanvas.mouseScroll(self, event)
    
       
  def zoomIn(self, *event):
  
    self.zoom *= 1.1
    self.zoom = min(5, self.zoom)
    self.drawCanvas()
  
  
  def zoomOut(self, *event):

    self.zoom *= 0.9
    self.zoom = max(0.5, self.zoom)
    self.drawCanvas()


  def resetZoom(self, *event):

    self.zoom = 1.0
    self.drawCanvas()


  def drawCanvas(self, width=None, height=None):

    canvas = self.canvas
    createText = canvas.create_text
    createRect = canvas.create_rectangle
    
    if not width: # Zero and None
      width  = int(canvas.cget('width')) 
    width  = max(100, width)

    if not height: 
      height = int(canvas.cget('height'))
    height = max(100, height)

    size = min(width,height)

    width  = self.zoom * size
    height = self.zoom * size
         
    pad = 5
    x0  = y0 = pad*8
    x1  = width  - x0
    y1  = height - y0
    
    bold  = 'Helvetica 14 bold'
    big  = 'Helvetica 21'
    
    for item in self.axisObjects:
      canvas.delete(item)
      
    self.axisObjects =[]
    
    # Y (psi) axis 
    yM  = (y1+y0)/2
    yM1 = (yM+y0)/2
    yM2 = (y1+yM)/2
    
    if self.yTicks:
      it0 = createText(x0-pad, y0, anchor='e', text="180")
      it1 = createText(x0-pad, yM1,anchor='e', text="90")
      it2 = createText(x0-pad, yM, anchor='e', text="0")
      it3 = createText(x0-pad, yM2,anchor='e', text="-90")
      it4 = createText(x0-pad, y1, anchor='e', text="-180")
      self.axisObjects.extend([it0,it1,it2,it3,it4])
    
    if self.yLabel: 
      it5 = createText(x0-15,  yM, anchor='e', text=self.yLabel, font=big)
      self.axisObjects.append(it5)
    
    # X (phi) axis 
    xM  = (x1+x0)/2
    xM1 = (xM+x0)/2
    xM2 = (x1+xM)/2
    
    if self.xTicks:
      it6 = createText(x0, y1+pad,anchor='n', text="-180")
      it7 = createText(xM1,y1+pad,anchor='n', text="-90")
      it8 = createText(xM, y1+pad,anchor='n', text="0")
      it9 = createText(xM2, y1+pad,anchor='n', text="90")
      itA = createText(x1, y1+pad,anchor='n', text="180")
      self.axisObjects.extend([it6,it7,it8,it9,itA])
    
    if self.xLabel:
      itB = createText(xM, y1+15, anchor='n', text=self.xLabel, font=big)
      self.axisObjects.append(itB)
       
    itC = createRect(x0,y0,x1,y1, outline='#000000')

    itD = createRect(0,0,x1+pad,y1+pad, outline=self.canvas.cget('bg'))
    
    # Title
    self.title = createText(xM, y0-(4*pad), text=self.titleText, font=bold)

    self.axisObjects += [itC,itD,self.title]
    

    # Draw the plot itself
    self.plotCoords = (x0,y0,x1,y1)
    self.drawRamachandran(x0,y0,x1,y1) # left,top,bottom, right

    self.refresh()

  def drawRamachandran(self,x0,y0,x1,y1):
  
    width  = x1-x0
    height = y1-y0
  
    matrix = self.matrix
    canvas = self.canvas
    nBins  = 360 / self.binSize   
    binWidth  = self.binWidth  = width/float(nBins)
    binHeight = self.binHeight = height/float(nBins)
    #self.cirRadius = min(5, self.binWidth/4.0)
    rangeNbins = range(nBins)
    coordFunc  = canvas.coords

    # draw bin matrix column by column
    if not matrix:
      createFunc = canvas.create_rectangle
      for x in rangeNbins:
        row = []
        bx0 = x0 + (x*binWidth)
        bx1 = bx0 + binWidth
 
        for y in rangeNbins:
          by0 = y0 + (y*binHeight)
          by1 = by0 + binHeight
          bin = createFunc(bx0,by0,bx1,by1,fill="#FFFFFF", outline="grey80")
          row.append(bin)
 
        matrix.append(row)

    else: # refresh matrix given canvas size
      for x in rangeNbins:
        bx0 = x0 + float(x*binWidth)
        bx1 = bx0 + binWidth
        row = [coordFunc(matrix[x][y],bx0,y0 + (y*binHeight),bx1,y0 + (y*binHeight) + binHeight) for y in rangeNbins]


    if self.defaultPlot:
      self.showAminoPlot()
    
    self.updateObjects()

    """
    # create crosshair at (0,0)
    mid = nBins / 2 * binLength
    self.phi = self.canvas.create_line(mid, 0, mid, mid*2, fill="yellow", width=2)
    self.psi = self.canvas.create_line(0, mid, mid*2, mid, fill="yellow", width=2)
    """

  def setAminoAcid(self, ccpCode):
  
    if self.aminoAcid != ccpCode:
      self.aminoAcid = ccpCode
      self.showAminoPlot()

  def scaleColorQuick(self, color, value, dict={}):
  
    if dict.get(color) is None:
      dict[color] = {}
    
    d = dict[color]  
    if d.get(value) is None:
      color2 = scaleColor(self, color, value)
      d[value] = color2
      return color2
      
    else:
      return d[value]

  def showAminoPlot(self):
    
    matrix   = self.matrix
    binWidth = self.binWidth
    canvas   = self.canvas
    itemconf = canvas.itemconfigure
    scaleColorQuick = self.scaleColorQuick

    if self.aminoAcid:
      colors = self.canvasDict.get(self.aminoAcid)
       
      if not colors:
        # fileName is the name of a PGM file in P2 format
        fileName = joinPath(self.bitmapDir,'log_%s.pgm' % self.aminoAcid.upper())
        colors = self.readPGM(fileName)
        self.canvasDict[self.aminoAcid] = colors

      if not colors:
        showWarning('Warning','Ramachandran plot data for %s does not exist! Using default' % self.aminoAcid)
        fileName = joinPath(self.bitmapDir,'log_ALL.pgm')
        colors = self.readPGM(fileName)
        self.canvasDict[self.aminoAcid] = colors

    else:
      colors = self.canvasDict.get('ALL')
      
      if not colors:
        fileName = joinPath(self.bitmapDir,'log_ALL.pgm')
        colors = self.readPGM(fileName)
        self.canvasDict['ALL'] = colors


    for i in range(len(colors)):
      for j in range(len(colors[0])):
        item  = matrix[i][j]
        color = colors[i][j]
        itemconf(item, fill=color)
        if binWidth < 7:
          itemconf(item, outline=color)
        elif binWidth < 12:
          itemconf(item, outline=scaleColorQuick(color,0.9))
        else:
          itemconf(item, outline=scaleColorQuick(color,0.8))

    itemconf(self.title, text=self.aminoAcid or 'All Residue Types')

  def getIntensityValue(self, phi, psi):
  
    intensities = self.intensityDict.get(self.aminoAcid)
    if not intensities:
      if self.aminoAcid:
        fileName = joinPath(self.bitmapDir,'log_%s.pgm' % self.aminoAcid.upper())
      else:
        fileName = joinPath(self.bitmapDir,'log_ALL.pgm')
 
      self.readPGM(fileName)
      intensities = self.intensityDict[self.aminoAcid]
      
    i = int((phi+180)/self.binSize)
    j = int((psi-180)/-self.binSize)
            
    return intensities[i][j]
    

  def readPGM(self, fileName):
    # return a list of populations specified by pgm file

    binSize = self.binSize

    file  = open(fileName, "r")
    # magic number "P2"
    # comment
    # dims 
    # max intensity
    # real intensity info
    
    lines = file.readlines()
    dims  = lines[2].split()
    lines = lines[4:]
    
    # initialize intensities[i][j] to 0
    x = int(dims[0]) / binSize
    y = int(dims[1]) / binSize
    intensities = [[0] * y for i in range(x)]

    # Here the intensities are read row by row corresponding to the matrix.
    i      = 0
    j      = 0
    x      = 0
    y      = 0
    #maxVal = 0
    total  = 0.0
    for line in lines:
      if i >= int(dims[0]):
        i = 0
        j = j + 1
        y = j / binSize
      x = i / binSize
      value = int(line.strip())
      intensities[x][y] += value
      total += value
      #if intensities[x][y] > maxVal:
      #  maxVal = intensities[x][y]
        
      i += 1

    #maxVal = float(maxVal)

    # convert relative intensity into color representation
    getPlotColor = self.getPlotColor
    maxInt = self.maxInt
    colors = []
    for i in range(len(intensities)):
      row  = []
      for j in range(len(intensities[0])):
        v = intensities[i][j] / total

        if v < 6.564e-5: # PROCKECK
          color = '#FFFFFF'
        elif v < 0.000821: # PROCKECK
          color = '#D0D0D0'
        else:
          color = getPlotColor(v, maxInt)
        row.append(color)
        
        intensities[i][j] = v

      colors.append(row)

    self.intensityDict[self.aminoAcid] = intensities

    return colors


  def getPlotColor(self, i, maxInt):
 
    j = i**0.5

    red   = 0.6 + ( 0.4 * j )
    green = 0.6 * (1-j)
    blue  = green

    return "#%02x%02x%02x" % (red*maxInt, green*maxInt, blue*maxInt)


  def updateCrosshair(self, phi, phiSD, psi, psiSD):
    # phi and psi must lie between -180 and 180
    
    self.canvas.delete(self.phi)
    self.canvas.delete(self.psi)
    self.canvas.delete(self.sd1)
    self.canvas.delete(self.sd2)

    # crosshairs
    x = self.getCoordFromAngle(phi, "phi") 
    y = self.getCoordFromAngle(psi, "psi")
    length = self.binLength * len(self.matrix)
    self.phi = self.canvas.create_line(x, 0, x, length, fill="yellow", width=2)
    self.psi = self.canvas.create_line(0, y, length, y, fill="yellow", width=2)

    # oval of 1 standard deviation
    x0 = self.getCoordFromAngle(phi-phiSD, "phi")
    x1 = self.getCoordFromAngle(phi+phiSD, "phi")
    y0 = self.getCoordFromAngle(psi-psiSD, "psi")
    y1 = self.getCoordFromAngle(psi+psiSD, "psi")
    self.sd1 = self.canvas.create_oval(x0, y0, x1, y1, outline="yellow", width=2)

    # oval of 2 standard deviation
    x0 = self.getCoordFromAngle(phi-phiSD*2, "phi")
    x1 = self.getCoordFromAngle(phi+phiSD*2, "phi")
    y0 = self.getCoordFromAngle(psi-psiSD*2, "psi")
    y1 = self.getCoordFromAngle(psi+psiSD*2, "psi")
    self.sd2 = self.canvas.create_oval(x0, y0, x1, y1, outline="yellow", width=2)


  def updateObjects(self, phiPsiAccList=None, objectList=None, labels=None, colors=None):
 
    # phiPsiAccList is a list of tuples (phi, psi, accept_state)
    # objectList is a list of estimation objects

    canvas = self.canvas
    scaleColorQuick = self.scaleColorQuick

    if objectList is not None:
      self.currentObjects = objectList
    else:
      objectList = self.currentObjects

    if labels is not None:
      self.labels = labels
    else:
      labels = self.labels  

    if phiPsiAccList is not None:
      self.phiPsiList = phiPsiAccList
    else:
      phiPsiAccList = self.phiPsiList

    if colors is not None:
      self.objectColors = colors
    else:
      colors = self.objectColors

    N = len(phiPsiAccList)

    n = N - len(objectList)
    if n > 0:
      objectList += [None] * n

    n = N - len(labels)
    if n > 0:
      labels += [None] * n

    n = N - len(colors)
    if n > 0:
      colors += [None] * n

    # delete existing estimation circles
    delete = canvas.delete
    residueLabels = self.residueLabels
    for text, index in residueLabels.values():
      delete(text)
    
    residueLabels = self.residueLabels = {}

    # draw new estimation circles
    r = self.cirRadius
    cirTag = self.cirTag
    currentCirDict    = self.currentCirDict  = {}
    currentObjDict    = self.currentObjDict
    selectedObjects   = self.selectedObjects = []
    getCoordFromAngle = self.getCoordFromAngle
    
    for i in range(N):
      phi, psi, accept = phiPsiAccList[i]
      if (phi is None) or (psi is None):
        continue
 
      obj = objectList[i]
      color  = colors[i] or '#808080'
      outline = scaleColorQuick(color,0.8)
 
      phi = int(getCoordFromAngle(phi, "phi"))
      psi = int(getCoordFromAngle(psi, "psi"))
      
      circle = currentObjDict.get(obj)
      if circle:
        canvas.coords(circle, phi-r, psi-r, phi+r, psi+r)
        canvas.itemconfigure(circle, fill=color, outline=outline)
      else:
        circle = canvas.create_oval(phi-r, psi-r, phi+r, psi+r, fill=color,
                                    outline=outline)
        currentObjDict[obj] = circle
        
      currentCirDict[circle] = obj

      if accept:
        selectedObjects.append(obj)
        label = labels[i]
      else:
        label = ''
      
      text = canvas.create_text(phi+r, psi+r, text=label, fill=outline, anchor='sw')
      residueLabels[obj] = (text, i)

    objs = currentObjDict.keys()
    for obj in objs:
      cir = currentObjDict[obj]
      if currentCirDict.get(cir) is None:
        canvas.delete(cir)
        del currentObjDict[obj]

    lift = canvas.lift
    for label, index in residueLabels.values():
      lift(label)
        
  def rejectObject(self, obj):
  
    circle = self.currentObjDict.get(obj)
    if circle:
      self.canvas.itemconfigure(circle, fill=self.unselectColor)
 
      if obj in self.selectedObjects:
        self.selectedObjects.remove(obj)


  def acceptObject(self, obj):
  
    circle = self.currentObjDict.get(obj)
    if circle:
      self.canvas.itemconfigure(circle, fill=self.selectColor)
 
      if obj not in self.selectedObjects:
        self.selectedObjects.append(obj)


  def recalculate(self, q):
    # re-evaluate mean phi, mean psi and members of selection
    # according to the quadrant specified

    # 4 steps:
    # 1. ignore those lying outside q
    newSelect  = []
    phiPsiList = []
    
    for object in self.currentObjDict:
      coords = self.canvas.coords(self.currentObjDict[object])
      phi    = self.getAngleFromPlot(coords[0] + self.cirRadius, "phi")
      psi    = self.getAngleFromPlot(coords[1] + self.cirRadius, "psi")

      if (phi < 0):
        if ((psi > -90) and (psi < 90)):
          if (q == 2):
            newSelect.append(object)
            phiPsiList.append((phi, psi))
        else:
          if (q == 1):
            newSelect.append(object)
            phiPsiList.append((phi, psi))
      else:
        if ((psi > -90) and (psi < 90)):
          if (q == 4):
            newSelect.append(object)
            phiPsiList.append((phi, psi))
        else:
          if (q == 3):
            newSelect.append(object)
            phiPsiList.append((phi, psi))

    if (len(newSelect) <= 0):
      print "No estimation in selected quadrant!"
      return
    
    # 2. calculate crude mean and sd
    phiMean, phiSD, psiMean, psiSD = find_mean_sd(phiPsiList)

    # 3. throw away those outside (2 * sd)
    phiThres = 2 * phiSD
    psiThres = 2 * psiSD
    for object in newSelect:
      index    = newSelect.index(object)
      phi, psi = phiPsiList[index]
      if ((math.fabs(phi-phiMean) > phiThres) or
          (math.fabs(psi-psiMean) > psiThres)):
        newSelect.remove(object)
        phiPsiList.remove((phi,psi))
    
    # 4. calculate true mean and sd
    phiMean, phiSD, psiMean, psiSD = find_mean_sd(phiPsiList)

    
    # make necessary updates
    self.updateCrosshair(phiMean, phiSD, psiMean, psiSD)
    for object in self.currentObjDict:
      if object in newSelect:
        self.acceptObject(object)
      else:
        self.rejectObject(object)
    

    return (phiMean, phiSD, psiMean, psiSD)
    

  def getCoordFromAngle(self, angle, angleType):
    # transform an angle between -180 and 180 to a coordinate on Rama plot
    # according to its type (phi or psi)
    # type must either be "phi" or "psi"
    
    x0, y0, x1, y1 = self.plotCoords
    if angleType == "phi":
      x = (angle + 180) / float(self.binSize) * self.binWidth
      coord = x0 + x
      
    else:  
      y = (180 - angle) / float(self.binSize) * self.binHeight
      coord = y0 + y

    return coord

  def getAngleFromPlot(self, coord, angleType):
  
    # transform a coordinate in Rama plot to an angle between -180 and 180

    x0, y0, x1, y1 = self.plotCoords
    if angleType == "phi":
      x = coord - x0
      angle = (x / float(self.binWidth) * self.binSize) - 180
      
    else:
      y = coord - y0
      angle = 180 - (y / float(self.binHeight) * self.binSize)

    return angle

  def circleClick(self, event):

    # find the actual estimation object pointed by cursor
    x = self.canvas.canvasx(event.x)
    y = self.canvas.canvasy(event.y)
    idents = self.canvas.find_overlapping(x-1,y-1,x+1,y+1)
    
    for ident in idents:
      obj = self.currentCirDict.get(ident)
    
      if obj:
        (text, index) = self.residueLabels[obj] 
        phi, psi, accepted = self.phiPsiList[index]
        accepted = not accepted
      
        if accepted:
          self.canvas.itemconfigure(text, text=self.labels[index])
        
        else:
          self.canvas.itemconfigure(text, text='')
        
        self.phiPsiList[index] = phi, psi, accepted
        if self.callback:
          self.callback(obj)
          
    # change circle color
    """currentColor = self.canvas.itemcget(cirID, "fill")
    if currentColor == self.selectColor:
      self.canvas.itemconfigure(cirID, fill=self.unselectColor)
      self.selectedObjects.remove(object)

    else:
      self.canvas.itemconfigure(cirID, fill=self.selectColor)
      self.selectedObjects.append(object)"""


  def circleEnter(self, event):
  
    cirID  = self.canvas.find_closest(event.x, event.y)[0]
    object = self.currentCirDict.get(cirID)
    text   = self.residueDict.get(object)

    #if text:
      

    if object:
      self.enterCallback(object)
    

  def circleLeave(self, event):
  
    self.leaveCallback()


  def rightClick(self, event):
  
    x = self.getAngleFromPlot(event.x, "phi")
    y = self.getAngleFromPlot(event.y, "psi")
    mean_sd_tuple = None

    if x < 0:                   
      if ((y > -90) and (y < 90)):
        mean_sd_tuple = self.recalculate(2)
        q = 2
        
      else:
        mean_sd_tuple = self.recalculate(1)
        q = 1
        
    else:                         
      if ((y > -90) and (y < 90)):
        mean_sd_tuple = self.recalculate(4)
        q = 4
        
      else:
        mean_sd_tuple = self.recalculate(3)
        q = 3

    if self.rightCallback and mean_sd_tuple:
      self.rightCallback(q, mean_sd_tuple)


if (__name__ == "__main__"):
  
  import Tkinter
  
  root = Tkinter.Tk()
  root.grid_columnconfigure(0, weight=1)
  root.grid_rowconfigure(0, weight=1)
  
  rama = ViewRamachandranFrame(root)
  rama.grid(row=0, column=0, sticky='nsew')
  rama.setAminoAcid('Gly')
  
  phiPsi = [(-60,145,1),(-66,0,1),(-55,-45,1),(90,-90,0),(-120,150,1),]
  objList = [1,2,3,4,5]
  
  rama.updateObjects(phiPsi,objectList=objList)
  
  root.mainloop()

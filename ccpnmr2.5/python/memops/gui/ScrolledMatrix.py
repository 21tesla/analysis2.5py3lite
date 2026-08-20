
"""
======================COPYRIGHT/LICENSE START==========================

ScrolledMatrix.py: <write function here>

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
import re, sys
from datetime import datetime

import Tkinter

from memops.universal.Util import useWheelMouse, buttonPress, buttonMotion

from memops.gui.Base import getPopup
from memops.gui.BasePopup import BasePopup
from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList
from memops.gui.Canvas import Canvas
from memops.gui.Color import scaleColor, invertColorRgb, invertColor
from memops.gui.Entry import Entry
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.Menu import Menu
from memops.gui.FontMenu import FontMenu
from memops.gui.PulldownList import PulldownList
from memops.gui.RadioButtons import RadioButtons
from memops.gui.Scrollbar import Scrollbar
from memops.gui.ToolTip import ToolTip

float_re = re.compile('\s*(-?\d+\.*\d*|\d*\.*\d+)')
numStr_re = re.compile('(\D*)(-?\d+\.*\d*|\d*\.\d+)(.*)')

datetimeType = type(datetime.today())
floatType    = type(1.0)
integerType  = type(1)
stringType   = type('')
unicodeType  = type(u'')
months = (None,'Jan','Feb','Mar','Apr','May','Jun',
          'Jul','Aug','Sep','Oct','Nov','Dec') # Months start at 1
          
HEADING_CELL_ITEMS = set(['headingText','headingRect','headingSort','minimise','filter'])

class ScrolledMatrix(Frame):

  def __init__(self, parent, initialCols=None, initialRows=None, minCellWidth=3,
               minCellHeight=1, maxRows=None, multiSelect=False,
               headingBg='#A8A8A8', headingFg='black', headingOutline='#787878',
               cellBg='grey82', cellFg='black', cellOutline='grey75' ,
               specialBg='darkBlue', specialFg='white', selectBg='#D0B0A0', 
               selectFg='#501000', selectOutline='#908070', headingList=None,
               objectList=None, justifyList=None, textMatrix=None, 
               colorMatrix=None, editWidgets=None, editSetCallbacks=None,
               editGetCallbacks=None, callback=None, doubleCallback=None, 
               sorting=True, font='Helvetica 10', highlightType='object',
               deleteFunc=None, passSelfToCallback=False, docKey=None, tipTexts=None,
               *args, **kw):

    self.parent = parent

    if textMatrix is None:
      textMatrix  = [[],]
    if colorMatrix is None:
      colorMatrix = [[],]
    if  headingList is None:
      headingList = []
    if objectList is None:
      objectList  = []
    if justifyList is None:
      justifyList  = ['center' for x in range(len(headingList))]
    if editWidgets is None:
      editWidgets = []
    if editSetCallbacks is None:
      editSetCallbacks = []
    if editGetCallbacks is None:
      editGetCallbacks = []
    if not initialCols:
      initialCols = len(headingList)
    if initialRows:
      if maxRows:
        initialRows = min(initialRows,maxRows)
    else:
      if maxRows:
        initialRows = min(max(5, len(objectList)), maxRows)
      else:
        initialRows = min(max(5, len(objectList)), 10)
      
    assert len(editWidgets) == len(editSetCallbacks) == len(editGetCallbacks)
    
    if not deleteFunc:
      deleteFunc = self.deselectAll

    self.headingBg      = headingBg
    self.headingFg      = headingFg
    self.headingOutline = headingOutline
    self.cellBg         = cellBg
    self.cellFg         = cellFg
    self.cellOutline    = cellOutline
    self.specialBg      = specialBg
    self.specialFg      = specialFg
    self.selectBg       = selectBg
    self.selectFg       = selectFg
    self.selectOutline  = selectOutline
    self.rangeBg        = scaleColor(parent,selectBg,1.15)
    self.rangeFg        = scaleColor(parent,selectFg,1.15)
    self.rangeOutline   = scaleColor(parent,selectOutline,1.15)
    self.headingList    = list(headingList)
    self.objectList     = list(objectList)
    self.justifyList    = list(justifyList)
    self.editWidgets    = editWidgets
    self.editSetCallbacks  = editSetCallbacks
    self.editGetCallbacks  = editGetCallbacks
    self.textMatrix     = textMatrix
    self.colorMatrix    = colorMatrix
    self.callback       = callback
    self.doubleCallback = doubleCallback
    self.sorting        = sorting
    self.highlightType  = highlightType
    self.blank          = ''
    self.initialCols    = initialCols
    self.initialRows    = initialRows
    self.maxRows        = maxRows
    self.multiSelect    = multiSelect
    self.passSelfToCallback  = passSelfToCallback
    
    self.scrollbarWidth = 15
    
    self.font = font
    self.smallFont = 'Courier 10'
    self.smallBoldFont = 'Courier 10 bold'
    
    self.vOffset = 0
    self.hOffset = 0

    self.cellDict = {}
    self.cellTextMatrix  = [[],]
    self.cellRectMatrix  = [[],]
    self.cellMarkMatrix  = [[],]
    self.textHeading = []
    self.rectHeading = []
    self.minColBox = []
    self.minColVer = []
    self.minColHor = []
    #self.srtColBox = []
    self.filterMks = []
    self.srtColArr = []
    self.toggledCols = []
    self.filteredCols =  []
    self.currentObject  = None
    self.currentObjects = []
    self.currentCell   = [0,0]
    self.currentRange  = []
    self.lastSortLine  = None
    self.sortDown = True
    self.mouseCell = ['','','']
    self.tableExportSelection = {}

    self.editWidget = None
    self.editObject = None
    self.editWidgetItem = None
    self.editSetCallback = None

    self.headHeight = 28

    Frame.__init__(self, parent, docKey=docKey, *args, **kw)

    self.setupTipTexts(tipTexts)
    
    self.visibleRows = initialRows
    self.visibleCols = initialCols

    self.toolTip = ToolTip(self, text='')
    
    self.canvas = Canvas(self, relief='flat', borderwidth=0, background = self.cget('background'),width=50, height=(self.initialRows+1) * 30 )
    self.canvas.pack()
      
    self.bind('<Configure>',self.refreshSizeAfter)
    
    self.canvas.bind('<Button-1>', self.mouseClick)

    if useWheelMouse():
      self.canvas.bind('<MouseWheel>', self.windowsOsScroll)
    else:
      self.canvas.bind('<Button-4>', self.scrollUp)
      self.canvas.bind('<Button-5>', self.scrollDown)

    self.canvas.bind('<Double-1>', self.mouseDoubleClick)
    self.canvas.bind('<Motion>',   self.mouseEnter)
    self.canvas.bind('<Enter>',    self.mouseEnter)
    self.canvas.bind('<Leave>',    self.mouseLeave)
    self.canvas.bind('<f>', self.selectFilter)
    self.canvas.bind('<a>', self.selectAll)
    self.canvas.bind('<Control-a>', self.selectAll)
    self.canvas.bind('<p>', self.printCanvas)
    self.canvas.bind('<KeyPress-Prior>', self.keyPressPrior)
    self.canvas.bind('<KeyPress-Next>',  self.keyPressNext )
    self.canvas.bind('<KeyPress-Up>',    self.keyPressUp   )
    self.canvas.bind('<KeyPress-Down>',  self.keyPressDown )
    self.canvas.bind('<KeyPress-Right>', self.keyPressRight)
    self.canvas.bind('<KeyPress-Left>',  self.keyPressLeft )
    self.canvas.bind('<KeyPress-Home>',  self.keyPressHome )
    self.canvas.bind('<KeyPress-End>',   self.keyPressEnd  )
    self.canvas.bind('<KeyPress-Delete>', deleteFunc)
    
    self.filterCanvas = Canvas(self, relief='flat', borderwidth=0,
                               background = self.cget('background'),
                               width=self.scrollbarWidth,
                               height=self.headHeight)
    self.filterCanvas.bind('<Button-1>', self.selectFilter)
    self.filterCanvas.bind('<Enter>',    self.mouseEnterfilterCanvas)
    self.filterCanvas.bind('<Leave>',    self.mouseLeavefilterCanvas)
    self.filterCanvas.bind('f', self.selectFilter)
    self.filterCanvas.bind('p', self.printCanvas)
    self.filterBox   = None
    self.cornerCanvas = Canvas(self, relief='flat', borderwidth=0,
                               background = self.cget('background'),
                               width=self.scrollbarWidth,
                               height=self.scrollbarWidth)
    self.filterPopup = None
    self.graphPopup  = None
       
    self.horizScrollbar = Scrollbar(self, orient=Tkinter.HORIZONTAL, background=cellBg,
                                    troughcolor=headingBg, width=self.scrollbarWidth,
                                    borderwidth=1, callback=self.moveHScrollbar)
    self.vertScrollbar  = Scrollbar(self, orient=Tkinter.VERTICAL, width=self.scrollbarWidth,
                                    borderwidth=1, callback=self.moveVScrollbar, background=cellBg)

    infoItems = [{ 'kind': 'command', 'label':'','command': None },
                 { 'kind': 'command', 'label':'','command': None },
                 { 'kind': 'command', 'label':'','command': None }]

    self.menu = Menu(self, tearoff=False)
    menu_items = [ \
      { 'kind': 'command', 'label': 'Filter', 'command': self.selectFilter },
      { 'kind': 'command', 'label': 'Export', 'command': self.popupExport  },
      { 'kind': 'cascade', 'label': 'Graph',  'submenu': [] },
      { 'kind': 'command', 'label': 'Print',  'command': self.printCanvas  },
      { 'kind': 'cascade', 'label': 'Table info','submenu': infoItems  },
    ]

    self.menu.setMenuItems(menu_items)
    fontMenu = FontMenu(self.menu, self.setFont, sizes=(8,10,12), doBold=False, doBoldItalic=False, tearoff=False)
    self.menu.add_cascade(label = 'Font', shortcut='F', menu=fontMenu)
    self.setGraphMenu()

    # self.canvas.bind('<ButtonPress-3>', self.popupMenu)
    # self.canvas.bind('<B3-Motion>', self.doNothing)
    self.canvas.bind(buttonPress(3), self.popupMenu)
    self.canvas.bind(buttonMotion(3), self.doNothing)

    self.refreshScrollbars()
    self.drawCanvas()
  
  def popupMenu(self, event):
  
    self.setInfoMenu()
    self.menu.popupMenu(event)

  def setupTipTexts(self, tipTexts=None):

    objectDoc = self.objectDoc
    headingList = self.headingList
    if not tipTexts and objectDoc:
      tipTexts = []
      headings = [heading.replace('\n', ' ') for heading in headingList]
      countDict = {}
      for heading in headings:
        count = countDict.get(heading, -1) + 1
        columnDoc = objectDoc.get((heading, count))
        countDict[heading] = count
        if columnDoc:
          tipText = columnDoc.documentation
        else:
          tipText = None
        tipTexts.append(tipText)

    self.tipTexts = tipTexts or [None] * len(headingList)
    
    while len(self.tipTexts) < len(headingList):
      self.tipTexts.append(None)
    
  def windowsOsScroll(self, event):
  
    delta = event.delta
    if delta > 0:
      self.scrollUp(event)
    elif delta < 0:
      self.scrollDown(event)
   
  
  def scrollUp(self, event):
  
    numObjs = len(self.objectList)
    numRows = self.visibleRows
  
    if numObjs <= numRows:
      return

    delta   = numRows/3 or 1
    vOffset = max(0, self.vOffset-delta)

    if vOffset != self.vOffset:
      self.vOffset = vOffset
      if vOffset < 1:
        if numRows < numObjs:
          self.moveVScrollbar(0,0,False)
      else:
        prop  = float(vOffset)/numObjs
        prop2 = float(vOffset+numRows)/numObjs
        self.moveVScrollbar(prop,prop2,False)

      self.drawCanvas()


  def scrollDown(self, event):
  
    numObjs = len(self.objectList)
    numRows = self.visibleRows
 
    if numObjs <= numRows:
      return
  
    delta   = numRows/3 or 1
    vOffset = min(numObjs-numRows, self.vOffset+delta)

    if vOffset != self.vOffset:
      self.vOffset = vOffset
      if self.vOffset < 1:
        if numRows < numObjs:
          self.moveVScrollbar(0,0,False)
      else:
        prop  = float(self.vOffset)/numObjs
        prop2 = float(self.vOffset+numRows)/numObjs
        self.moveVScrollbar(prop,prop2,False)

      self.drawCanvas()

      
  def doNothing(self, event):
  
    pass
    
  def keyPressPrior(self, *event):

    if self.currentObject is not None:
      (row, col) = self.currentCell
      i = max(0, row-self.visibleRows)
      if i >= 0:
        self.selectNthObject(i)
  
  def keyPressNext(self, *event):

    if self.currentObject is not None:
      (row, col) = self.currentCell
      i = min(len(self.objectList)-1, row+self.visibleRows)
      if i < len(self.objectList):
        self.selectNthObject(i)
  
  def keyPressUp(self,*event):
    
    if self.currentObject is not None:
      (row, col) = self.currentCell
      i = row-1
      if i >= 0:
        self.selectNthObject(i)
  
  def keyPressDown(self, *event):

    if self.currentObject is not None:
      (row, col) = self.currentCell
      i = row+1
      if i < len(self.objectList):
        self.selectNthObject(i)
  
  def keyPressLeft(self, *event):
    pass
    
  def keyPressRight(self, *event):
    pass
  
  def keyPressHome(self, *event):
  
    if self.currentObject is not None:
      self.selectNthObject(0)
  
  def keyPressEnd(self, *event):
  
    if self.currentObject is not None:
      self.selectNthObject(len(self.objectList)-1)

  def printCanvas(self, *event):
  
    from memops.gui.FileSelect import FileType
    from memops.gui.FileSelectPopup import FileSelectPopup

    fileTypes = [  FileType('PostScript', ['*.ps']), FileType('All', ['*'])]
    fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
               title = 'Print canvas to file', dismiss_text = 'Cancel',
               selected_file_must_exist = False)

    fileName = fileSelectPopup.getFile()
        
    self.canvas.postscript(colormode='color',file=fileName)
  
  def setInfoMenu(self):
  
    (row, col) = self.currentCell
    #screenRow  = row-self.visibleRows
    #screenCol  = col-self.visibleCols
    
    if row is None:
      row = '-'
    else:
      row = str(row)  
    
    if col is None:
      col = '-'
    else:
      col = str(col)  
     
    submenu = self.menu.submenuDict['Table info']
    submenu.entryconfig(0, label='Total rows: %d' % len(self.objectList))
    submenu.entryconfig(1, label='Selected rows: %d' %  len(self.currentObjects))
    submenu.entryconfig(2, label='Current row, col: %s,%s' % (row, col))
         
  def setGraphMenu(self):
  
    options = self.getGraphMenuOptions()
    
    if not options:
      self.menu.entryconfig(2, state=Tkinter.DISABLED)

    else:
      self.menu.entryconfig(2, state=Tkinter.NORMAL)
     
  def getGraphMenuOptions(self):
  
    submenu = self.menu.submenuDict['Graph']
    items = []
    if not self.textMatrix:
      return 
    
    heads = ['Row Number',]
    if self.headingList:
      for i, heading in enumerate(self.headingList):
        heading = heading.replace('\n',' ').strip()
        heads.append(heading or 'Column %d' % (i+1))
    
    N     = len(self.textMatrix[0])
    opts  = [0,]
    types = (integerType, floatType, datetimeType)
    
    for col in range(N):
      colType = None
      for row in range(len(self.textMatrix)):
        if self.textMatrix[row][col] is not None:
          colType = type(self.textMatrix[row][col])
          if colType in types:
            opts.append(col+1)
            break
          elif colType is stringType:
            match = re.match(numStr_re, self.textMatrix[row][col])
            if match:
              opts.append(col+1)
              break
 
    dict = {}
    for opt in opts:
      name = 'X:%s' % (heads[opt])
      yItems = []
      for optY in opts:
        if opt is optY:
          continue
        nameY = 'Y:%s' % (heads[optY])
        yItems.append( {'kind': 'command','label': nameY,'command': lambda x=opt, y=optY: self.selectGraph(x, y)} )
      items.append( {'kind': 'cascade','label': name,'submenu': yItems} )
      
    submenu.setMenuItems(items)
    return items
       
  def setFont(self, font):
  
    self.font = font
    #self.refreshSize()
    self.drawCanvas()
  
  def popupExport(self, *event):

    # put here because otherwise get circular imports
    from memops.gui.TableExportPopup import TableExportPopup
    popup = TableExportPopup(self, headings=self.headingList, entries=self.textMatrix)
    popup.destroy()

  def refreshfilterCanvas(self, height):
  
    height += 1
    self.vertScrollbar.canvas.config(height=height)
    fc = self.filterCanvas
    width = fc.winfo_width()
    fc.config(height=height+2)
    fc.delete('[All]')
    self.filterBox = fc.create_rectangle(0,1,width-1,height,fill=self.headingBg, outline=self.headingOutline)
    fc.create_text((width-1)/2,height/2+2,text='?', font=self.font, fill=self.headingFg, justify='center')

  def mouseEnterfilterCanvas(self, *event):

    self.filterCanvas.focus()
    fill    = scaleColor(self,self.headingBg,1.1,)
    outline = scaleColor(self,self.headingOutline,1.1,)
    self.filterCanvas.itemconfigure(self.filterBox,fill=fill, outline=outline)

  def mouseLeavefilterCanvas(self, *event):

    self.filterCanvas.itemconfigure(self.filterBox,fill=self.headingBg, outline=self.headingOutline)

  def selectFilter(self, *event):
  
    if self.filterPopup:
      self.filterPopup.open()
    else:
      self.filterPopup = FilterMatrixPopup(getPopup(self),self)

  def selectGraph(self, xAxis, yAxis):
  
    if self.graphPopup:
      self.graphPopup.open()
      self.graphPopup.draw(xAxis, yAxis)
    else:
      self.graphPopup = GraphMatrixPopup(getPopup(self),self, xAxis, yAxis)
 
  def filterCol(self, col):
  
    col = col+self.hOffset
    self.selectFilter()
    self.filterPopup.setColumn(col)
 
  def keyPressEscape(self, *event):
  
    self.menu.popdownMenu()
    self.editSetCallback = None
    self.finishEditCell()
 
  def finishEditCell(self):

    if self.editWidget:
      self.canvas.focus_force()
      self.canvas.delete(self.editWidgetItem)
      self.editWidget = None
      if self.editSetCallback:
        self.editSetCallback(self.editObject)
      self.editObject = None

  def getEditWidget(self, i):

    return self.editWidgets[i]

  def editCell(self, row, col):
  
    if col >= len(self.editWidgets):
       return
    
    objRow = row + self.vOffset
    objCol = col + self.hOffset
    if not self.doEditMarkExtraRules(self.objectList[objRow], objRow, objCol):
      return
    
    if self.editWidget:
      self.keyPressEscape()
       
    self.editWidget = self.getEditWidget(objCol)
    self.editObject = self.currentObject
    self.editSetCallback = self.editSetCallbacks[objCol]

    editGetCallback = self.editGetCallbacks[objCol]
    if editGetCallback:
      editGetCallback(self.currentObject)  

    if not self.editWidget:
      return

    rectItem = self.cellRectMatrix[row][col]
    bbox = self.canvas.bbox(rectItem)
    x = 0.5*(bbox[0]+bbox[2])
    y = 0.5*(bbox[1]+bbox[3])
    
    self.editWidget.configure(borderwidth=1, relief='solid')
    self.editWidgetItem = self.canvas.create_window(x, y, window=self.editWidget, anchor='center') 
    self.update_idletasks()

    if not self.editWidget:
      return
    
    w = self.canvas.winfo_width()
    h = self.canvas.winfo_height()
    w2 = self.editWidget.winfo_width()/2.0
    h2 = self.editWidget.winfo_height()/2.0
        
    x0 = x-w2
    x1 = x+w2
    y0 = y-h2
    y1 = y+h2
    
    move = False    
    if y1 > h:
      y0 -= (y1-h)
      y1 = h
      move = True

    if x1 > w:
      x0 -= (x1-w)
      x1 = w
      move = True
    
    if x0 < 0:
      x1 -= x0
      x0 = 0
      move = True
      
    if y0 < 0:
      y1 -= y0
      y0 = 0
      move = True
    
    x2 = 0.5*(x0+x1)
    y2 = 0.5*(y0+y1)

    if move:
      self.canvas.coords(self.editWidgetItem, x2, y2)
    
    self.editWidget.focus_force()
    self.editWidget.bind('<KeyPress-Escape>', self.keyPressEscape)
    self.canvas.lift(self.editWidgetItem)

  def selectDoubleCell(self, row,col):

    objRow = row + self.vOffset
    objCol = col + self.hOffset
    self.currentCell = [objRow, objCol]
    if (len(self.objectList) > 0) and (objRow < len(self.objectList)):
      self.currentObject = self.objectList[objRow]
      if self.doubleCallback:
        self.doubleCallback(self.currentObject, objRow, objCol)
    else:
      self.currentObject = None

  def mouseDoubleClick(self, event):

    self.toolTip.deactivate()
    if self.editWidget:
      self.finishEditCell()

    x = self.canvas.canvasx(event.x)
    y = self.canvas.canvasy(event.y)
    obj = self.canvas.find('closest',x,y)[0]
    cell = self.cellDict.get(obj)
    if cell:
      (row,col,typ) = cell
      if self.multiSelect and (event.state & 4 or event.state & 1):
          if event.state & 4:
            self.toggleInRange(row,col)
          elif event.state & 1:
            self.selectRange(row,col)
      else:
        if self.toggledCols[col]:
          return
        if row < len(self.objectList):
          if self.doubleCallback:
            self.selectDoubleCell(row,col)
          if typ not in HEADING_CELL_ITEMS:
            #self.selectCell(row,col)
            if ((self.editWidgets and self.editSetCallbacks) or self.editGetCallbacks) and self.currentObject is not None:
              self.editCell(row,col)
          
 
  def mouseClick(self, event):
    
    self.menu.popdownMenu(event)
    self.toolTip.deactivate()

    if self.editWidget:
      self.finishEditCell()

    x = self.canvas.canvasx(event.x)
    y = self.canvas.canvasy(event.y)
    obj = self.canvas.find('closest',x,y)[0]
    cell = self.cellDict.get(obj)
    if cell:
      (row,col,typ) = cell
      if typ in ('headingText','headingRect','headingSort'):
        if self.sorting and not self.toggledCols[col]:
           self.sortLine(col)
      elif typ == 'minimise':
        self.toggleCol(col)
      elif typ == 'filter':
        self.filterCol(col)
      elif not self.toggledCols[col]:
        if row < len(self.objectList):
          if self.multiSelect and (event.state & 4 or event.state & 1):
            if event.state & 4:
              self.toggleInRange(row,col)
            elif event.state & 1:
              self.selectRange(row,col)
          else:
            self.selectCell(row,col)

  def mouseEnter(self, event):
    
    self.menu.popdownMenu()  
    canvas = self.canvas
    
    if not self.editWidgetItem:
      canvas.focus_force()

    x = canvas.canvasx(event.x)
    y = canvas.canvasy(event.y)
    t = canvas.find('closest',x,y)
    obj =  t[0]
    if obj == self.mouseCell[0]:
      return
        
    cell = self.cellDict.get(obj)
    self.toolTip.deactivate()
    
    
    if cell:
      (row,col,typ) = cell
      
      if typ in HEADING_CELL_ITEMS:
        tipText = self.tipTexts[col]
        
        if tipText:
          self.toolTip.text = tipText
          self.toolTip.activate()
      else:
        self.toolTip.deactivate()
        
      if not self.objectList:
        return
      
      if typ == 'cellText' or typ == 'headingText' or typ == 'cellMark':
        bbox = canvas.bbox(obj)
        objs = canvas.find('overlapping',bbox[0],bbox[1],bbox[2],bbox[3])
        
              
        for obj2 in objs:
          cell = self.cellDict.get(obj2)
          if cell and obj2:
            (row,col,typ) = cell
            if typ == 'cellRect' or typ == 'headingRect':
              obj = obj2
          if obj == self.mouseCell[0]:
            return
            
      (row,col,typ) = self.cellDict.get(obj)
      if typ == 'cellRect' or typ == 'headingRect':
        self.mouseLeave()
        fill1    = canvas.itemcget(obj,'fill')
        outline1 = canvas.itemcget(obj,'outline')
        self.mouseCell = [obj, fill1, outline1]
        fill2    = scaleColor(self,fill1,1.1,)
        outline2 = scaleColor(self,outline1,1.1)
        self.canvas.itemconfigure(obj,fill=fill2, outline=outline2)
   
  def mouseLeave(self, *event):
    
    if self.mouseCell[0]:
      (obj0, fill0, outline0) = self.mouseCell
      self.canvas.itemconfigure(obj0,fill=fill0, outline=outline0)
      if event:
        self.mouseCell = ['','','']
        self.toolTip.deactivate()
    
    else: # Left whole canvas
      self.toolTip.deactivate()

  def toggleCol(self, col):
  
    col = col+self.hOffset
    sum = 0
    for i in self.toggledCols:
      sum += i
     
    if (sum >= len(self.toggledCols)-1) and not self.toggledCols[col]:
      return
    
    self.toggledCols[col] = not self.toggledCols[col]
    self.mouseLeave()
    #self.refreshSize()
    self.drawCanvas()
    
  def formatText(self,text):
  
    textType = type(text)
  
    if text is None:
      text = ''
  
    elif textType == floatType:
      if text == 0:
        text = '0'
      elif abs(text) > 999999 or abs(text) < 0.01:
        text = '%5.3e' % text
      else:
        text = '%5.3f' % text
    
    elif textType is datetimeType:
      text = '%2d %s %d %2d:%2.2d' % (text.day, months[text.month], text.year, text.hour, text.minute)
    
    elif textType is unicodeType:
      text = text.encode('utf-8')
    
    elif textType != stringType:
      text = str(text)
      
    if text and text[0:1] == '@':
      text = ''

    return text

  def refreshSizeAfter(self, event, scrollRefresh=True):
  
    self.after_idle(lambda: self.refreshSize(event, scrollRefresh))

  def refreshSize(self, event=None, scrollRefresh=True):
    
    self.keyPressEscape()
    pad = 2
    if event:
      width  = event.width
      height = event.height
    else:
      width  = int(self.winfo_width())
      height = int(self.winfo_height())

    height -= self.scrollbarWidth + 1
    width  -= self.scrollbarWidth + 1

    lastCol = 1
    cellHeight = 1
    headHeight = 1
    extraRowList = []
    rowWidth = 0
    textRows = len(self.textMatrix)
    columns = range(self.hOffset,len(self.headingList))
    unseen  = range(self.hOffset)
    unseen.reverse()
    for x in unseen:
      columns.append(x)

    textItem = self.canvas.create_text(0,0,text='A', font=self.font)
    coords   = self.canvas.bbox(textItem)
    self.canvas.delete(textItem)
    
    cellHeight = coords[3] - coords[1] + (3*pad)
    visibleRows = None
    
    for j in columns:
      j2 = j - self.hOffset

      headHeight0 = 1
      if ( j2 < len(self.textHeading) ) and (j2 >= 0):
        if self.textHeading[j2]:
          coords = self.canvas.bbox(self.textHeading[j2])
          headHeight0 = coords[3] - coords[1] + (3*pad)
          headHeight = max(headHeight, headHeight0)

    for j in columns:
      if j >= len(self.toggledCols):
        self.toggledCols.append( 0 )
         
      y = headHeight
      
      if self.toggledCols[j]:
        cellWidth = 4*pad
        row = int((height-headHeight)/cellHeight)
      
      else:
        row = self.vOffset
        while y < height:
          if row < textRows:
            if j < len(self.textMatrix[row]):
              tryText = self.textMatrix[row][j]
              textType = type(tryText)
              if textType in (stringType, unicodeType):
                numLines = tryText.count('\n') + 1
                y += numLines*cellHeight
              else:
                y += cellHeight
            else:
              y += cellHeight
          else:
            y += cellHeight
          
          row += 1
        row -= 1
            

        text = self.headingList[j].strip() or 'AAA'
        textItem = self.canvas.create_text(0,0,text=text, font=self.font)
        coords   = self.canvas.bbox(textItem)
        cellWidth = coords[2] - coords[0] + (4*pad)
        self.canvas.delete(textItem)

      numRows = row-self.vOffset
      if (visibleRows is None) or (numRows < visibleRows):
        visibleRows = numRows
 
      rowWidth += cellWidth

      if rowWidth < width:
        if j >= self.hOffset:
          lastCol = j
        else:
          self.hOffset = j
          
      else:
        break     
          
    self.visibleCols = max(min(lastCol+1-self.hOffset, len(self.headingList) ), 1)
    
    if self.maxRows:
      self.visibleRows = max(min(self.maxRows,visibleRows), 1)
    else:
      self.visibleRows = max(visibleRows, 1)
      
    self.headHeight = headHeight
         
    if scrollRefresh:
      self.refreshScrollbars()
    
    self.drawCanvas()
    
  def drawCanvas(self):
    
    blank          = self.blank
    cellBg         = self.cellBg
    cellDict       = self.cellDict
    cellFg         = self.cellFg
    cellMarkMatrix = self.cellMarkMatrix
    cellOutline    = self.cellOutline
    cellRectMatrix = self.cellRectMatrix
    cellTextMatrix = self.cellTextMatrix
    colorCell      = self.colorCell  
    colorMatrix    = self.colorMatrix
    doEditMark     = self.doEditMark
    editWidgets    = self.editWidgets
    filteredCols   = self.filteredCols
    filterMks      = self.filterMks
    font           = self.font
    formatText     = self.formatText
    headingBg      = self.headingBg
    headingFg      = self.headingFg
    headingList    = self.headingList
    headingOutline = self.headingOutline
    justifyList    = self.justifyList
    lastSortLine   = self.lastSortLine
    minColBox      = self.minColBox
    minColHor      = self.minColHor
    minColVer      = self.minColVer
    rectHeading    = self.rectHeading
    selectBg       = self.selectBg
    selectOutline  = self.selectOutline
    smallBoldFont  = self.smallBoldFont
    smallFont      = self.smallFont
    sortDown       = self.sortDown
    srtColArr      = self.srtColArr
    #srtColBox      = self.srtColBox
    textHeading    = self.textHeading
    textMatrix     = self.textMatrix
    toggledCols    = self.toggledCols
    visibleCols    = self.visibleCols
    visibleRows    = self.visibleRows
    scrollbarWidth = self.scrollbarWidth

    canvas = self.canvas
    coords = canvas.coords
    itconf = canvas.itemconfigure
    lift   = canvas.lift
    lower  = canvas.lower
    cbbox  = canvas.bbox
    createRect = canvas.create_rectangle
    createLine = canvas.create_line
    createText = canvas.create_text
    createPoly = canvas.create_polygon
    deleteCanvasItems = self.deleteCanvasItems

    canvasWidth  = int(self.winfo_width()) - 1
    canvasHeight = int(self.winfo_height()) -1

    if len(headingList) > visibleCols: 
      canvasHeight -= scrollbarWidth
 
    if len(textMatrix) > visibleRows:
      canvasWidth  -= scrollbarWidth + 1
            

    headingLine = 1
    if len(headingList) == 0:
      headingLine = 0
    
    N  = self.visibleCols
    ho = self.hOffset
    vo = self.vOffset

    colWidths = []
    pad  = 2
    pad2 = 2*pad
    y    = 1
    if headingLine:
      x = 1
      for i in range(N):
        colWidths.append(0)
        text = ''
        if i+ho >= len(toggledCols):
          toggledCols.append(0)
        if i+ho >= len(filteredCols):
          filteredCols.append(0)
        if i+ho >= len(editWidgets):
          editWidgets.append(None)
        if i+ho >= len(justifyList):
          justifyList.append('center')
        if toggledCols[i+ho]:
          text = ''
        elif i+ho < len(headingList):
          try:
            text = headingList[i+ho]
          except:
            print 'Missing heading in ScrolledMatrix input'
            
        justify = justifyList[i+ho]
        if (i >= len(textHeading)) or not textHeading[i]:
          textItem = createText(x+pad,y+pad,text=text, font=font, fill=headingFg, justify='center')
          bbox = cbbox(textItem)
          x2 = bbox[2] + pad2 + pad2
          y2 = bbox[3] + pad2
          rectItem = createRect(x,y,x2,y2,fill=headingBg, outline=headingOutline)
          minimBox = createRect(x+7,y,x,y+7,fill=headingBg, outline=headingOutline)
          minimVer = createLine(x+7,y+3,x,y+3,fill=headingOutline)
          minimHor = createLine(x+3,y,x+3,y+7,fill=headingOutline)
          #sortsBox = createRect(x2-7,y-7,x2,y2,fill=headingBg, outline=headingOutline)
          filterMk = createText(x2-7,y-7,text='?', font=smallFont, fill=headingOutline, justify='center')
          sortsArr = createPoly(x+7,y-7,x,y2-3,x+7,y2,fill=selectBg, outline=selectOutline)
          if i >= len(textHeading):
            textHeading.append( textItem )
            rectHeading.append( rectItem )
            minColBox.append( minimBox )
            minColVer.append( minimVer )
            minColHor.append( minimHor )
            #srtColBox.append( sortsBox )
            filterMks.append( filterMk )
            srtColArr.append( sortsArr )
          else:
            textHeading[i] = textItem
            rectHeading[i] = rectItem
            minColBox[i] = minimBox
            minColVer[i] = minimVer
            minColHor[i] = minimHor
            #srtColBox[i] = sortsBox
            filterMks[i] = filterMk 
            srtColArr[i] = sortsArr
          cellDict[textItem] = (0,i,'headingText')
          cellDict[rectItem] = (0,i,'headingRect')
          cellDict[minimBox] = (0,i,'minimise')
          cellDict[minimVer] = (0,i,'minimise')
          cellDict[minimHor] = (0,i,'minimise')
          #cellDict[sortsBox] = (0,i,'heading')
          cellDict[filterMk] = (0,i,'filter')
          cellDict[sortsArr] = (0,i,'headingSort')

        textItem = textHeading[i]
        coords(textItem, x+pad, y+pad)
        itconf(textItem, fill=headingFg, text=text, font=font)
        bbox = cbbox(textItem)
        if toggledCols[i+ho]:
          colWidth = 3*pad
        else:
          colWidth = bbox[2] - bbox[0] + (2*pad)
        if colWidth > colWidths[i]:
          colWidths[i] = colWidth
        x = x + colWidth + pad

      rowHeight = 0
      for i in range(N):
        bbox   = cbbox(textHeading[i])
        height = bbox[3] - bbox[1]
        if height > rowHeight:
          rowHeight = height
      rowHeight += 2*pad
      y = y + rowHeight + pad
            
    for i in range(visibleRows):
      if i >= len(cellTextMatrix):
        cellTextMatrix.append( [] )
        cellRectMatrix.append( [] )
        cellMarkMatrix.append( [] )
        
      x = 1
      rowHeight = 0
      for j in range(N):
        if j >= len(colWidths):
          colWidths.append(0)
 
        text    = None
        color   = cellBg
        justify = justifyList[j+ho]

        if (justify == 'left'):
          anchor = Tkinter.W
        elif (justify == 'center'):
          anchor = Tkinter.CENTER
        else:
          anchor = Tkinter.E

        if toggledCols[j+ho]:
          text = ''
        else:
          if i+vo < len(textMatrix):
            if j < len(textMatrix[i+vo]):
              try:
                text = textMatrix[i+vo][j+ho]
              except:
                print 'Missing text in ScrolledMatrix input'
              try:
                color = colorMatrix[i+vo][j+ho]
              except:
                pass
 
          if text is not None:
            text = formatText(text);
 
        if (j >= len(cellTextMatrix[i])) or not cellTextMatrix[i][j]:
          textItem = createText(x+pad,y+pad,text=text, font=font, fill=cellFg, justify=justify, anchor=anchor)
          bbox = cbbox(textItem)
          w = bbox[2] - bbox[0] + pad2
          h = bbox[3] - bbox[1] + pad2
          rowHeight = max(h, rowHeight)
          rectItem = createRect(x,y,x+w,y+h,fill=color, outline=cellOutline)
          markItem = createLine(x,y,x,y,x,y,x,y,x,y,x,y,x,y,x,y, fill=headingBg)
          if j >= len(cellTextMatrix[i]):
            cellTextMatrix[i].append( textItem )
            cellRectMatrix[i].append( rectItem )
            cellMarkMatrix[i].append( markItem )
          else:
            cellTextMatrix[i][j] = textItem 
            cellRectMatrix[i][j] = rectItem
            cellMarkMatrix[i][j] = markItem
          cellDict[textItem] = (i,j,'cellText')
          cellDict[rectItem] = (i,j,'cellRect')
          cellDict[markItem] = (i,j,'cellMark')
        
        else:       
          textItem = cellTextMatrix[i][j]
          coords(textItem, x+pad, y+pad)
        
        if text is None:
          itconf(textItem, fill=headingFg, font=font, text=blank,
                 justify=justify, anchor=anchor)
        else:
          itconf(textItem, fill=headingFg, font=font, text=text,
                 justify=justify, anchor=anchor)
          
        bbox = cbbox(textItem)
        colWidth  = bbox[2] - bbox[0] + pad2
        rowHeight = max(rowHeight, bbox[3] - bbox[1] + pad2)
        
        if colWidth > colWidths[j]:
          colWidths[j] = colWidth
        
        x = x + colWidth + pad
      y = y + rowHeight + pad
    
    totalX = 0
    M = 0
    for i in range(N):
      totalX += colWidths[i] + pad
      if toggledCols[i+ho]:
        M += 1
    L = N-M
    spareX = canvasWidth - totalX
    extraX = max(0, int(spareX/L) )
    errorX = spareX % L
    for i in range(N):
      if not toggledCols[i+ho]:
        colWidths[i] = colWidths[i] + extraX
        
    i=0
    while errorX > 0:
      if not toggledCols[i+ho]:
        colWidths[i] += 1
        errorX -= 1
      i+=1
      i = i % N
        
    #extraY = max(0, int( (int(self.innerFrame.winfo_height()) - y)/visibleRows ) )
    extraY = max(0, int((canvasHeight - y)/visibleRows) )
    #errorY = (int(self.innerFrame.winfo_height()) - y) % visibleRows
    errorY = (canvasHeight - y) % visibleRows
    # final headings
    
    y = 1
    if headingLine:
      x = 0
      rowHeight = 0
      for i in range(N):
        bbox   = cbbox(textHeading[i])
        height = bbox[3] - bbox[1]
        if height > rowHeight:
          rowHeight = height
      rowHeight += pad2

      x = 1
      for i in range(N):
        rectItem = rectHeading[i]
        textItem = textHeading[i]
        minimBox = minColBox[i]
        minimVer = minColVer[i]
        minimHor = minColHor[i]
        #sortsBox = #srtColBox[i]
        filterMk = filterMks[i]
        sortsArr = srtColArr[i]

        w2 = int(colWidths[i]/2)
        h2 = int(rowHeight/2) + 1
        x2 = x+colWidths[i]
        y2 = y+rowHeight
        coords(textItem, x+w2,y+h2)
        coords(rectItem, x,   y,   x2,  y2)
        coords(minimBox, x+6,y,   x,  y+6)
        coords(minimVer, x+3,y+1, x+3,y+6)
        coords(minimHor, x+6,y+3, x,  y+3)
        coords(filterMk, x2-4, y+6)
        #coords(#sortsBox, x2-6,y2-6,x2,  y2)
        if i+ho == lastSortLine:
          lift(sortsArr)
          if sortDown:
            coords(sortsArr, x+6,y2-6,x,  y2-6, x+3,y2)
          else:
            coords(sortsArr, x+6,y2,  x+3,y2-6, x,  y2)
        else:
          lower(sortsArr)
          coords(sortsArr, x+6,y2-6,x,  y2-3, x+6,  y2)
        lift(minimBox)
        if toggledCols[i+ho]:
          itconf(minimBox, fill=selectBg, outline=selectOutline)
          itconf(minimVer, fill=selectOutline)
          itconf(minimHor, fill=selectOutline)
          lift(minimVer)
          lower(filterMk)
        else:
          itconf(minimBox, fill=headingBg, outline=headingOutline)
          itconf(minimVer, fill=headingOutline)
          itconf(minimHor, fill=headingOutline)
          lower(minimVer)
          lift(filterMk)
        if filteredCols[i+ho]:
          itconf(filterMk, fill=selectOutline, text='!', font=smallBoldFont)
        else:
          itconf(filterMk, fill=headingOutline, text='?', font=smallFont)
        lift(minimHor)
        lift(textItem)
        x = x + colWidths[i] + pad
          
      if len(textHeading) > N:
        for i in range(N,len(textHeading)):
          deleteCanvasItems([textHeading[i],rectHeading[i],filterMks[i],minColBox[i],minColVer[i],minColHor[i],srtColArr[i]])
          (textHeading[i],rectHeading[i],filterMks[i],minColBox[i],minColVer[i],minColHor[i],srtColArr[i]) = (None, None, None, None, None, None, None)
      y = y + rowHeight + pad
    else:
      for i in range(len(textHeading)):
        deleteCanvasItems([textHeading[i],rectHeading[i],filterMks[i],minColBox[i],minColVer[i],minColHor[i],srtColArr[i]])
        (textHeading[i],rectHeading[i],filterMks[i],minColBox[i],minColVer[i],minColHor[i],srtColArr[i]) = (None, None, None, None, None, None, None)
        
    self.refreshfilterCanvas(rowHeight)
    
    headHeight = y
    # final cells
    for i in range(visibleRows):
      x = 1
      rowHeight = 0
      for j in range(N):
        bbox   = cbbox(cellTextMatrix[i][j])
        height = bbox[3] - bbox[1]
        if height > rowHeight:
          rowHeight = height
      rowHeight += (2*pad) + extraY
      if errorY > 0:
        rowHeight +=1
        errorY -=1

      x = 1
      for j in range(N):
        rectItem = cellRectMatrix[i][j]
        textItem = cellTextMatrix[i][j]
        markItem = cellMarkMatrix[i][j]

        w2 = int(colWidths[j]/2) + pad
        h2 = int(rowHeight/2) + 1
        x2 = x+colWidths[j]
        y2 = y+rowHeight
       
        justify = justifyList[j+ho]
        if (justify == 'left'):
          xOffset = pad2
        elif (justify == 'center'):
          xOffset = w2
        else:
          xOffset = colWidths[j] - pad2
        coords(textItem, x+xOffset, y+h2)
        coords(rectItem, x, y, x2,y2)
        if toggledCols[j+ho]:
          coords(markItem, 3,3,3,3,3,3,3,3,3,3,3,3,3,3)
          lower(markItem)
        elif doEditMark(i+vo, j+ho):
          coords(markItem, x+6,y2-2,x+2,y2-2,x+2,y2-5,x+4,y2-5,x+2,y2-5,x+2,y2-8,x+6,y2-8)
          lift(markItem)
        else:
          coords(markItem, 3,3,3,3,3,3,3,3,3,3,3,3,3,3)
          lower(markItem)
        itconf(markItem, fill=headingOutline)
        itconf(textItem, fill=cellFg)
        itconf(rectItem, fill=cellBg, outline=cellOutline)
        lift(textItem)
        x = x + colWidths[j] + pad

        color = cellBg
        try:
          color = colorMatrix[i+vo][j+ho]
          colorCell(i,j,color)
        except:
          pass
     
      if len(cellTextMatrix[i]) > visibleCols:
        for j in range(visibleCols,len(cellTextMatrix[i])):
          deleteCanvasItems([cellTextMatrix[i][j],cellRectMatrix[i][j],cellMarkMatrix[i][j]])
          cellTextMatrix[i][j] = None
          cellRectMatrix[i][j] = None
          cellMarkMatrix[i][j] = None
      
      y = y + rowHeight + pad

    if len(cellTextMatrix) > visibleRows:
      for i in range(visibleRows,len(cellTextMatrix)):
        for j in range(len(cellTextMatrix[i])):
          deleteCanvasItems([cellTextMatrix[i][j],cellRectMatrix[i][j],cellMarkMatrix[i][j]])
          cellTextMatrix[i][j] = None
          cellRectMatrix[i][j] = None
          cellMarkMatrix[i][j] = None
        
    # sort object and highlights
    if self.currentObject is not None:
      self.hilightObject(self.currentObject, focus=False)
        
    if self.multiSelect:
      self.hilightRange()
  
    self.canvas.place(x=0,y=0,width=canvasWidth,height=canvasHeight)

  
  def doEditMark(self, row, col):
  
    if row < len(self.textMatrix) and col < len(self.textMatrix[row]):
      if col < len(self.editGetCallbacks) and self.editGetCallbacks[col]:
        #if self.textMatrix[row][col] is not None:
        return self.doEditMarkExtraRules(self.objectList[row], row, col)
    
    return False

  def doEditMarkExtraRules(self, obj, row, col):
  
    # dummy function to be overwritten
    return True
  
  def deleteCanvasItems(self,items):
  
    c = self.canvas
    d = self.cellDict
    for item in items:
      c.delete(item)
      d[item] = None
   
  def sortLine(self, col, noUpdate=None):
  
    if len(self.textMatrix) < 2:
      return

    if not self.sorting:
      return

    objCol = col + self.hOffset
    matrix = self.textMatrix
    colors = self.colorMatrix
    array = []
    
    if objCol >= len(matrix[0]):
      return
    
    for i, textRow in enumerate(matrix):
      sortKey = textRow[objCol]
      if type(sortKey) is stringType:
        match = re.match(numStr_re,sortKey)
        if match:
          sortKey = (match.group(1),float(match.group(2)),match.group(3))
        else:
          sortKey = (sortKey, -1.0e20, '')
      
      array.append( (sortKey, i) )
     
    if noUpdate:
      array.sort()
      indices = [x[1] for x in array]
      if not self.sortDown:
        indices.reverse()
    elif self.lastSortLine == objCol:
      indices = [x[1] for x in array]
      indices.reverse()
      self.sortDown = not self.sortDown
    else:
      array.sort()
      indices = [x[1] for x in array]
      self.sortDown = True
    
    nRows = len(matrix)
    newMatrix = [None] * nRows
    newColors = [None] * nRows
    newList = [None] * nRows
    newRange = []
    currentObjects = self.currentObjects
    objectList = self.objectList 
    
    for i in range(nRows):
      j = indices[i]
      newMatrix[i] = matrix[j]

      if j < len(colors):
        newColors[i] = colors[j]
      
      object = objectList[j]
      newList[i] = object
      
      if object in currentObjects:
        newRange.append(i)

    obj = self.currentObject
    
    self.currentObjects = [newList[x] for x in newRange]
    self.objectList = newList
    self.textMatrix = newMatrix
    self.colorMatrix = newColors
    self.currentRange = newRange
    self.lastSortLine = objCol

    if not noUpdate:
      if obj is not None:
        self.hilightObject(obj, focus=True)
      self.drawCanvas()
  
  def selectNthObject(self, n, focus=True):
  
    (row,col) = self.hilightObject(None, n, focus=focus)
    self.selectCell(row, col)

  def selectObject(self, obj):
  
    (row,col) = self.hilightObject(obj)
    self.selectCell(row, col)

  def hilightObject(self, obj, i=None, focus=True):
  
    if i is None:
      if obj not in self.objectList:
        return [0, 0]
      else:
        i = self.objectList.index(obj)
 
    self.currentObject = self.objectList[i]
    col = self.currentCell[1]
    self.currentCell = [i,col]

    if focus:
      self.hOffset = 0
      vOffset = self.vOffset
      row = self.offsetRow(i)
      if vOffset != self.vOffset:
        if self.vOffset < 1:
          if self.visibleRows < len(self.objectList):
            self.moveVScrollbar(0,0,False)
        else:
          prop  = float(self.vOffset)/len(self.objectList)
          prop2 = float(self.vOffset+self.visibleRows)/len(self.objectList)
          self.moveVScrollbar(prop,prop2,False)

        row = self.currentCell[0] - self.vOffset
        col = self.currentCell[1] - self.hOffset
        self.drawCanvas()
        
    else:
      row = i - self.vOffset
 
    if self.currentObject is not None:
      if (row >= 0) and (row < self.visibleRows):
        if (col >= 0) and (col < self.visibleCols):
          if type(self.highlightType) is integerType:
            col = self.highlightType
          self.hilightCell(row, col)

    return [row,col]
                                                
  def hilightRange(self):
  
    self.clearHighlights()
  
    if self.currentObjects:
      if self.currentObject not in self.currentObjects:
        self.currentRange.sort()
        objRow = self.currentRange[0]
        self.currentObject = self.objectList[objRow]
        self.currentCell   = [objRow - self.vOffset, 0]
        if self.callback:
          if self.passSelfToCallback:
            self.callback(self.currentObject, objRow,self.currentCell[1], self)
          else:
            self.callback(self.currentObject, objRow,self.currentCell[1])

      for objRow in self.currentRange:
        row = objRow - self.vOffset
        if (row >= 0) and (row < self.visibleRows):
          if self.currentObject is self.objectList[objRow]:
            for col in range(self.visibleCols):
              self.hilightSelectItems(row,col)
            
          else:
            for col in range(self.visibleCols):
              self.hilightInRangeItems(row,col)

  def hilightInRangeItems(self, row, col):
  
    color  = self.getVisibleColor(row, col)
    if color:
      self.inverseColorCell(row,col,color)
      return
    
    color  = self.getVisibleColor(row, col)
    canvas = self.canvas
    textItem = self.cellTextMatrix[row][col]
    rectItem = self.cellRectMatrix[row][col]
    markItem = self.cellMarkMatrix[row][col]
    canvas.itemconfigure(rectItem, fill=self.rangeBg, outline=self.rangeOutline)
    canvas.itemconfigure(textItem, fill=self.rangeFg)
    canvas.itemconfigure(markItem, fill=self.rangeOutline)
    if self.mouseCell[0] == rectItem:
      self.mouseCell[1] = self.rangeBg
      self.mouseCell[2] = self.rangeOutline

  def hilightSelectItems(self, row, col):
  
    color  = self.getVisibleColor(row, col)
    if color:
      self.inverseColorCell(row,col,color)
      return
    
    canvas = self.canvas
    textItem = self.cellTextMatrix[row][col]
    rectItem = self.cellRectMatrix[row][col]
    markItem = self.cellMarkMatrix[row][col]
    canvas.itemconfigure(rectItem, fill=self.selectBg, outline=self.selectOutline)
    canvas.itemconfigure(textItem, fill=self.selectFg)
    canvas.itemconfigure(markItem, fill=self.selectOutline)
    if self.mouseCell[0] == rectItem:
      self.mouseCell[1] = self.selectBg
      self.mouseCell[2] = self.selectOutline

  def hilightSpecialItems(self, row, col):
  
    color  = self.getVisibleColor(row, col)
    if color:
      self.inverseColorCell(row,col,color)
      return
    
    canvas = self.canvas
    textItem = self.cellTextMatrix[row][col]
    rectItem = self.cellRectMatrix[row][col]
    markItem = self.cellMarkMatrix[row][col]
    canvas.itemconfigure(rectItem, fill=self.specialBg, outline=self.specialFg)
    canvas.itemconfigure(textItem, fill=self.specialFg)
    canvas.itemconfigure(markItem, fill=self.specialFg)
    if self.mouseCell[0] == rectItem:
      self.mouseCell[1] = self.specialBg
      self.mouseCell[2] = self.specialFg

  def hilightCell(self, row, col):              
          
    if self.highlightType is not None:                       

      self.clearHighlights()
      if self.highlightType in ('object','cross'):
        for j in range(self.visibleCols):
          self.hilightSelectItems(row,j)
          #if self.multiSelect:
          #  self.hilightSelectItems(row,j)
          #else:
          #  self.hilightSelectItems(row,j)
               
      if (col >= 0) and (col < self.visibleCols): 
        if self.highlightType == 'cross':
          for i in range(self.visibleRows):
            self.hilightSelectItems(i,col)
            
        if self.highlightType in ('cell','cross'):
          self.hilightSpecialItems(row,col)

        elif type(self.highlightType) is integerType:
          self.hilightSpecialItems(row,self.highlightType)
             
  def clearHighlights(self):
  
    for i in range(self.visibleRows):
       for j in range(self.visibleCols):
         color = self.getVisibleColor(i,j)
         if color:
           self.colorCell(i,j,color)
         else:
           self.uncolorCell(i,j)

  def getVisibleColor(self, row, col):
  
    try:
      color = self.colorMatrix[row+self.vOffset][col+self.hOffset]
    except:
      color = None

    return color

  def deselectAll(self, *event):
    
    if self.editWidget:
      self.finishEditCell()
      
    self.currentObject  = None
    self.currentObjects = []
    self.currentCell    = [None,None]
    self.currentRange   = []
    self.mouseCell      = ['','','']
    self.clearHighlights()

  def deselect(self, object=None, i=None):

    if not object and (i is not None):
      object = self.objectList[i]
    
    elif object and i is None:
      i = self.objectList.index(object)
        
    if object:
      if i == self.currentCell[0]:
        self.currentObject  = None
        self.currentCell    = [None,None]
        self.mouseCell      = ['','','']
        self.currentRange.remove(i)
        self.currentObjects.remove(object)
     
      elif i in self.currentRange:
        self.currentRange.remove(i)
        self.currentObjects.remove(object)
    
      self.hilightRange()

    else:
      self.deselectAll()

  def selectAll(self, *event):
    
    if not self.multiSelect:
      return
    
    indices = range(len(self.objectList))
    self.currentObjects = [self.objectList[i] for i in indices]
    self.currentRange   = [i for i in indices]
    
    if not self.currentObject:
      objRow = max(0, min(self.vOffset, len(self.objectList)-1))
      objCol = self.hOffset
      self.currentObject  = self.objectList[objRow]
      self.currentCell    = [objRow,objCol]
     

    self.hilightRange()
 
  def selectObjects(self, objects):
  
    if not self.multiSelect:
      return
    
    objectList = self.objectList
    indexDict = dict([(obj,i) for i, obj in enumerate(objectList)])
    
    indices = []
    for object in objects:
      index = indexDict.get(object)
      
      if index is not None:
        indices.append(index)
    
    indices.sort()
    self.currentObjects = [objectList[i] for i in indices]
    self.currentRange = indices
    
    if self.currentObjects and (self.currentObject not in objects):
      objRow = indices[0]
      objCol = self.hOffset
      self.currentObject  = self.currentObjects[0]
      self.currentCell    = [objRow,objCol]

    self.hilightRange()
    
 
  def selectRange(self,row,col):
    
    objRow = max(0, min(row + self.vOffset, len(self.objectList)-1))
    objCol = col + self.hOffset

    if self.currentRange:
      self.currentRange.sort()
      if objRow < self.currentRange[0]:
        start = objRow
        end   = self.currentRange[-1]
      elif objRow > self.currentRange[-1]:
        start = self.currentRange[0]
        end   = objRow
      else:
        start = self.currentRange[0]
        end   = self.currentRange[-1]
      
    elif self.currentObject is not None:
      objRow0 = self.currentCell[0]
      if objRow0 > objRow:
        start = objRow
        end   = objRow0
      else:
        start = objRow0
        end   = objRow
    else:
      start = objRow
      end   = objRow
    
    extent = min(len(self.objectList)-1, end)
    
    
    self.currentObjects = self.objectList[start: end+1]
    self.currentRange = range(start, end+1)
        
    if self.objectList:
      self.currentObject  = self.objectList[objRow]
      self.currentCell    = [objRow, objCol]
      if self.callback:
        if self.passSelfToCallback:
          self.callback(self.currentObject, objRow,objCol, self)
        else:
          self.callback(self.currentObject, objRow,objCol)

    self.hilightRange()

  def toggleInRange(self,row,col):
    
    objRow = row + self.vOffset
    objCol = col + self.hOffset
    
    if (len(self.objectList) > 0) and (objRow < len(self.objectList)):
      object = self.objectList[objRow]
      if object in self.currentObjects:
        self.currentRange.remove(objRow)
        self.currentObjects.remove(object)
      else:
        self.currentRange.append(objRow)
        self.currentRange.sort()
        self.currentObjects = [self.objectList[x] for x in self.currentRange]
        self.currentObject  = self.objectList[objRow]
        self.currentCell    = [objRow, objCol]
      if self.callback:
        if self.passSelfToCallback:
          self.callback(self.currentObject, objRow,objCol, self)
        else:
          self.callback(self.currentObject, objRow,objCol)
          
      self.hilightRange()

  def selectCell(self, row, col):
   
    objRow = row + self.vOffset
    objCol = col + self.hOffset
    self.currentCell = [objRow, objCol]
    if (len(self.objectList) > 0) and (objRow < len(self.objectList)):
      self.currentObject  = self.objectList[objRow]
      self.currentObjects = [self.currentObject]
      self.currentRange   = [objRow]
      if self.callback:
        if self.passSelfToCallback:
          self.callback(self.currentObject, objRow, objCol, self)
        else:
          self.callback(self.currentObject, objRow, objCol)

    else:
      self.currentObject  = None
      self.currentObjects = []
      self.currentRange   = []
      self.mouseCell      = ['','','']

    if self.currentObject is not None:
      self.hilightCell(row, col)
      
  def offsetRow(self, i):
  
    N = len(self.objectList)

    if i < self.vOffset:
      self.vOffset = min(i, (N-self.visibleRows) )
    elif i > (self.vOffset + self.visibleRows - 1):
      self.vOffset = max( (i-self.visibleRows + 1) ,0)
    row = i - self.vOffset
    return row
  
  def offsetCol(self, i):

    N = len(self.headingList)
  
    if i < self.hOffset:
      self.hOffset = min(i, (N-self.visibleCols) )
    elif i > (self.hOffset + self.visibleCols - 1):
      self.hOffset = max(0, (i-self.visibleCols + 1) )
    col = i - self.hOffset
    return col
  
  def colorCell(self, row, col, color):
    
    itemconf = self.canvas.itemconfigure
    textItem = self.cellTextMatrix[row][col]
    rectItem = self.cellRectMatrix[row][col]
    markItem = self.cellMarkMatrix[row][col]
    invColor = invertColor(self, color)
    
    itemconf(rectItem, fill=color, outline=self.cellOutline)
    itemconf(textItem, fill=invColor)
    itemconf(markItem, fill=invColor)
    if self.mouseCell[0] == rectItem:
      self.mouseCell[1] = color
      self.mouseCell[2] = self.cellOutline

  def inverseColorCell(self, row, col, color):
  
    #colors = [color, self.cellOutline, self.cellFg, self.headingOutline]
    #for i in range(4):
    #  colors[i] = invertColorRgb(self, colors[i])
    
    itemconf = self.canvas.itemconfigure
    invCol   = invertColorRgb(self, color)
    contrast = invertColor(self, invCol)
    colors   = [invCol, self.cellOutline, contrast, contrast]
     
    textItem = self.cellTextMatrix[row][col]
    rectItem = self.cellRectMatrix[row][col]
    markItem = self.cellMarkMatrix[row][col]
    itemconf(rectItem, fill=colors[0], outline=colors[1])
    itemconf(textItem, fill=colors[2])
    itemconf(markItem, fill=colors[3])
    if self.mouseCell[0] == rectItem:
      self.mouseCell[1] = colors[0]
      self.mouseCell[2] = colors[1]

  def uncolorCell(self, row, col):
     
    itemconf = self.canvas.itemconfigure
    textItem = self.cellTextMatrix[row][col]
    rectItem = self.cellRectMatrix[row][col]
    markItem = self.cellMarkMatrix[row][col]
    itemconf(rectItem, fill=self.cellBg, outline=self.cellOutline)
    itemconf(textItem, fill=self.cellFg)
    itemconf(markItem, fill=self.headingOutline)
    if self.mouseCell[0] == rectItem:
      self.mouseCell[1] = self.cellBg
      self.mouseCell[2] = self.cellOutline
     
  def moveHScrollbar(self, prop, prop2, refresh=True):
  
    self.keyPressEscape()
    N = len(self.headingList)
      
    gaps = prop+1-prop2

    if gaps == 0 or N == 0:
      self.hOffset = 0
      self.horizScrollbar.set(0,1)
    else:
      self.hOffset = max(0, int(round(prop*N)))
      frac = min(1.0, float(self.visibleCols)/N )
      gaps = 1-frac
      start, end = min(gaps, max(0, prop)), max(frac, min(1.0, prop + frac))
      self.horizScrollbar.set(start, end)
           
    if refresh:
      self.drawCanvas()

  def moveVScrollbar(self, prop, prop2, refresh=True):

    self.keyPressEscape()
    N = len(self.objectList)
    
    gaps = prop+1-prop2
        
    if gaps == 0 or N == 0:
      self.vOffset = 0
      self.vertScrollbar.set(0,1)
    else:
      self.vOffset = max(0, int(round(prop*N)))
      self.refreshSizeAfter(None, False)
      frac = min(1.0, float(self.visibleRows)/N )
      gaps = 1-frac
      start, end = min(gaps, max(0, prop)), max(frac, min(1.0, prop + frac))
      self.vertScrollbar.set(start, end)
        
    if refresh:
      self.drawCanvas()
        

  def update(self, headingList=None, objectList=None, colorMatrix=None, textMatrix=None,
            justifyList=None, editWidgets=None, editSetCallbacks=None, editGetCallbacks=None,
            tipTexts=None, resetScrollbars=True):
  
    self.keyPressEscape()
      
    if justifyList is not None:
      self.justifyList = list(justifyList)
    
    if headingList is not None:
      if self.graphPopup and self.textMatrix and self.textMatrix[0]:
        if self.graphPopup.xAxis and self.graphPopup.yAxis:
          xHeading = self.headingList[self.graphPopup.xAxis-1]
          yHeading = self.headingList[self.graphPopup.yAxis-1]
          if (xHeading in headingList) and (yHeading in headingList):
            self.graphPopup.xAxis = list(headingList).index(xHeading) + 1
            self.graphPopup.yAxis = list(headingList).index(yHeading) + 1
          else:
            self.graphPopup.xAxis = None
            self.graphPopup.yAxis = None
 
      if (self.lastSortLine is not None) and (self.lastSortLine < len(self.headingList)):
        refHeading = self.headingList[self.lastSortLine]
        if refHeading in headingList:
          self.lastSortLine = list(headingList).index(refHeading)
          
      else:
        self.lastSortLine = None
      
      self.headingList = list(headingList)
      if self.filterPopup:
        self.filterPopup.updateColumns()
      
      n = len(headingList)
      m = len(self.justifyList)
      if n > m:
        for x in range(m,n):
          self.justifyList.append('center')
          self.toggledCols.append(False)
          self.filteredCols.append(False)

    self.currentRange   = []
    if objectList is None:
      self.currentObject  = None
      self.mouseCell = ['','','']
      self.currentObjects = []
    else:
      self.objectList = list(objectList)
      
      if objectList and (isinstance(objectList[0], list) or isinstance(objectList[-1], list)):
        for i, obj in enumerate(self.objectList):
          if isinstance(obj, list):
            self.objectList[i] = tuple(obj)
      
      if self.currentObject not in self.objectList:
        self.currentObject = None
        self.mouseCell = ['','','']
      for object in list(self.currentObjects):
        if object not in self.objectList:
          self.currentObjects.remove(object)
        else:
          self.currentRange.append( self.objectList.index(object) )

    if editWidgets is not None:
      self.editWidgets = editWidgets

    if editSetCallbacks is not None:
      self.editSetCallbacks = editSetCallbacks

    if editGetCallbacks is not None:
      self.editGetCallbacks = editGetCallbacks

    if colorMatrix is not None:
      self.colorMatrix = colorMatrix
  
    if textMatrix is not None:
      self.textMatrix = textMatrix
      for i in range(len(self.filteredCols)):
        self.filteredCols[i] = False
        self.toggledCols[i] = False
        
      if len(self.colorMatrix) > len(self.textMatrix):
        self.colorMatrix = self.colorMatrix[:len(self.textMatrix)]
    
    if (tipTexts is not None) or headingList:
      self.setupTipTexts(tipTexts)
    
    if self.lastSortLine is not None:
      self.sortLine(self.lastSortLine, noUpdate=True)
      
    if self.filterPopup:
      self.filterPopup.textMatrix  = self.textMatrix
      self.filterPopup.colorMatrix = self.colorMatrix
      self.filterPopup.objectList  = self.objectList
          
    if self.graphPopup and self.textMatrix and self.textMatrix[0]:
      self.graphPopup.draw()
    
    limit = len(self.objectList) - self.visibleRows
    if self.vOffset > limit:
      self.vOffset = max(0, limit)

    limit = len(self.headingList) - self.visibleCols
    if self.hOffset > limit:
      self.hOffset = max(0, limit)
          
    self.setGraphMenu()
    self.refreshSize()
    
  def refreshScrollbars(self, reset=False):

    width  = int(self.winfo_width())-1
    height = int(self.winfo_height())-1

    headingLine = 1
    if len(self.headingList) == 0:
      headingLine = 0

    N = len(self.headingList)
    
    useHoriz = False
    useVert  = False
    size = self.scrollbarWidth
    
    if N <= self.visibleCols:
      self.hOffset = 0
      self.horizScrollbar.set(0,1)
    else:
      if reset:
        prop = 0.0
      else:  
        prop = self.hOffset/float(N)
      frac = min(1.0, float(self.visibleCols)/N )
      gaps = 1-frac
      start, end = min(gaps, max(0, prop)), max(frac, min(1.0, prop + frac))
      self.horizScrollbar.set(start, end)
      useHoriz = True

       
    N = len(self.objectList)
    
    if N <= self.visibleRows:
      self.vOffset = 0
      self.vertScrollbar.set(0,1)
    else:
      if reset:
        prop = 0.0
      else:
        prop = self.vOffset/float(N)
        
      frac = min(1.0, float(self.visibleRows)/N )
      gaps = 1-frac
      start, end = min(gaps, max(0, prop)), max(frac, min(1.0, prop + frac))
      self.vertScrollbar.set(start, end)
      useVert = True

    if useHoriz:
      x2 = width-size
      y2 = height-size
      if useVert:
        self.horizScrollbar.place(x=1,y=y2,width=x2-2,height=size)
        self.vertScrollbar.place(x=x2,y=self.headHeight+1,width=size,height=y2-self.headHeight-2)
        self.filterCanvas.place(x=x2,y=0,width=size,height=self.headHeight)
        self.cornerCanvas.place(x=x2,y=y2,width=size,height=size)
        
        w2 = size-1
        self.cornerCanvas.delete('[All]')
        self.cornerCanvas.create_rectangle(0,0,w2,w2,fill=self.headingBg, outline=self.headingOutline)
        
      else:
        self.horizScrollbar.place(x=1,y=y2,width=width-2,height=size)
        self.vertScrollbar.place_forget()
        self.filterCanvas.place_forget()
        self.cornerCanvas.place_forget()
    
    else:
      self.horizScrollbar.place_forget()    
      self.cornerCanvas.place_forget()
      
      if useVert:
        self.vertScrollbar.place(x=width-size,y=self.headHeight+1,width=size,height=height-self.headHeight-2)
        self.filterCanvas.place(x=width-size,y=0,width=size,height=self.headHeight)
      else:
        self.vertScrollbar.place_forget()
        self.filterCanvas.place_forget()
    
    self.update_idletasks()
      
    
  def destroy(self):
  
    if self.filterPopup:
      self.filterPopup.destroy()
      
    if self.graphPopup:
      self.graphPopup.destroy()
    
    Frame.destroy(self)

class GraphMatrixPopup(BasePopup):


  def __init__(self, parent, scrolledMatrix, xAxis, yAxis):

    self.scrolledMatrix = scrolledMatrix
    self.xAxis          = xAxis or 0
    self.yAxis          = yAxis or 0
    BasePopup.__init__(self, parent=parent, title="Columns Graph",transient=False, borderwidth=6)

    parent.protocol("WM_DELETE_WINDOW", self.close)

  def body(self, guiFrame):
    
    from ScrolledGraph import ScrolledGraph
    from ButtonList    import UtilityButtonList
 
    guiFrame.grid_rowconfigure(0, weight = 1)
    guiFrame.grid_columnconfigure(0, weight = 1)

    graph = self.scrolledGraph = ScrolledGraph(guiFrame, showCoords=True, grid=(0,0))

    texts = ['Zoom In','Zoom Out','Reset Zoom',
             'Line','Scatter','Bar Chart','Histogram']
    commands = [lambda : graph.setZoom(graph.zoom*1.2),
                lambda : graph.setZoom(graph.zoom/1.2),
                lambda : graph.setZoom(1.0),
                lambda : graph.setGraphType('line'),
                lambda : graph.setGraphType('scatter'),
                lambda : graph.setGraphType('bar'),
                lambda : graph.setGraphType('histogram')]
    bottomButtons = UtilityButtonList(guiFrame, texts=texts, doHelp=False, doClone=False,
                                      commands=commands, expands=True)
    bottomButtons.grid(row=1, column=0, sticky=Tkinter.EW)

    self.draw(self.xAxis,self.yAxis)

  def draw(self, xAxis=None, yAxis=None):
  
    if xAxis is not None:
      self.xAxis = xAxis
      
    if yAxis is not None:
      self.yAxis = yAxis
  
    dataSet = []
    if (self.yAxis is None) or (self.yAxis is None):
      xLabel = ''
      yLabel = ''
      dataSet = [[0,0],[0,0]]
    
    else:
      for i in range(len(self.scrolledMatrix.textMatrix)):
        row = self.scrolledMatrix.textMatrix[i]
        if not row:
          continue

        if self.xAxis > 0:
          x =  row[self.xAxis-1]
          if type(x) is stringType:
            match = re.match(numStr_re,x)
            if match:
              x = float(match.group(2))
            else:
              x = None  
              continue
          
        else:
          x = i+1

        if self.yAxis > 0:
          y =  row[self.yAxis-1]
          if type(y) is stringType:
            match = re.match(numStr_re,y)
            if match:
              y = float(match.group(2))
            else:
              y = None  
              continue
              
        else:
          y = i+1
 
        dataSet.append( [x,y] )

      if self.xAxis > 0:
        xLabel  = self.scrolledMatrix.headingList[self.xAxis-1]
        xLabel  = re.sub('\n',' ', xLabel)
      else:
        xLabel = 'Row'
 
      if self.yAxis > 0:
        yLabel  = self.scrolledMatrix.headingList[self.yAxis-1]
        yLabel  = re.sub('\n',' ', yLabel)
      else:
        yLabel = 'Row'
      
    title   = 'Scrolled Matrix Graph - %s vs %s' % (xLabel,yLabel)
    
    self.scrolledGraph.update(dataSets=[dataSet,],xLabel=xLabel,yLabel=yLabel,title=title)

SEARCH_MODES = [ 'Literal','Case Sensitive Literal','Regular Expression' ]

class FilterMatrixPopup(BasePopup):

  def __init__(self, parent, scrolledMatrix):

    self.scrolledMatrix = scrolledMatrix
                       
    BasePopup.__init__(self, parent=parent, title="Filter Options",
                       transient=False, borderwidth=6)
    

    parent.protocol("WM_DELETE_WINDOW", self.close)
    self.geometry("+150+150")
    
  def body(self, guiFrame):
    
    guiFrame.grid_rowconfigure(1, weight = 1)
    guiFrame.grid_columnconfigure(0, weight = 0)
    guiFrame.grid_columnconfigure(1, weight = 0)
    guiFrame.grid_columnconfigure(2, weight = 1)
  
    objEntries = [ 'All', 'Selected' ]
    columnNames = self.getColumnNames()

    self.textMatrix  = self.scrolledMatrix.textMatrix
    self.colorMatrix = self.scrolledMatrix.colorMatrix
    self.objectList  = self.scrolledMatrix.objectList

    self.mode = SEARCH_MODES[0]
    self.objsToBeFiltered = None
    self.column = '[All]'

    label        = Label(guiFrame, text='Filter column:', grid=(0,0))
    self.pulldownMenu = PulldownList(guiFrame, texts=columnNames,
                                     callback=self.changeColumn, grid=(0,1))
    self.entry = Entry(guiFrame, text='', returnCallback=self.filterInclude)
    modeButtons  = RadioButtons(guiFrame, SEARCH_MODES, self.selectMode, borderwidth=0, relief='flat')
    
    texts = ['Close','Reset','Filter\nInclude','Filter\nExclude']
    commands = [self.close, self.unfilterMatrix,
                self.filterInclude, self.filterExclude]
    buttons = ButtonList(guiFrame, texts=texts, commands=commands,
                         grid=(3,0), gridSpan=(1,3))
    
    objFrame     = Frame(guiFrame, grid=(0,2))
    objLabel     = Label(objFrame, text='Objects to filter:', grid=(0,0))
    objButtons   = RadioButtons(objFrame, objEntries,
                                self.filteredObjs, borderwidth=0,
                                relief='flat', grid=(0,1))

    self.entry.grid(row=1, column=0, columnspan=3, sticky=Tkinter.NSEW)
    modeButtons.grid(row=2, column=0, columnspan=3, sticky=Tkinter.EW)
  
  def updateColumns(self):
  
    columnNames = self.getColumnNames()
    index = self.pulldownMenu.getSelectedIndex()
    if index >= len(columnNames):
      index = 0
    
    self.pulldownMenu.setup(columnNames, columnNames, index)
  
  def setColumn(self, col):
  
    self.column = col
    self.pulldownMenu.setIndex(col+1)
    
  def changeColumn(self, name):
  
    if name == '[All]':
      self.column = '[All]'
    else:
      self.column = self.pulldownMenu.getSelectedIndex() - 1
      
  def unfilterMatrix(self):
  
    for i in range(len(self.scrolledMatrix.filteredCols)):
      self.scrolledMatrix.filteredCols[i] = 0
    
    self.scrolledMatrix.update(textMatrix  = self.textMatrix,
                               colorMatrix = self.colorMatrix,
                               objectList  = self.objectList)

  def filterInclude(self, *event):
  
    self.filterMatrix(True)
  
  def filterExclude(self, *event):
  
    self.filterMatrix(False)
     
  def filterMatrix(self, includeMatches=True):
  
    if not self.objectList:
      return
  
    string = self.entry.get()
    if not string:
      return
    
    col = self.column
    flag = re.S
    
    def exclude(a,b,c):
      return not re.search(a,b,c)
    
    if includeMatches:
      find = re.search
    else:
      find = exclude
    
    if self.mode != SEARCH_MODES[2]:
      string  = re.escape(string)
    
    if self.mode == SEARCH_MODES[0]:
      flag = re.I
    
    scrolledMatrix = self.scrolledMatrix
    textMatrix     = self.textMatrix
    colorMatrix    = self.colorMatrix

    scrolledMatrix.textMatrix  = []
    scrolledMatrix.colorMatrix = []
    scrolledMatrix.objectList  = []

    appendText = scrolledMatrix.textMatrix.append
    appendObject = scrolledMatrix.objectList.append
    appendColor = scrolledMatrix.colorMatrix.append
    formatText = scrolledMatrix.formatText
    
    currentObjects = set(scrolledMatrix.currentObjects)
    selfObjectList = self.objectList  
    objsToBeFiltered = self.objsToBeFiltered

    newObjs  = []
    newRange = []
    for row, textRow in enumerate(textMatrix):
      object = selfObjectList[row]

      if objsToBeFiltered == 'Selected':
        if object not in currentObjects:
          continue

      if col == '[All]':
        nMatch = 0
        for datum in textRow:
          datum = formatText(datum)
          match = find(string,datum,flag)

          if match:
            nMatch += 1
          
            if includeMatches:
              break
         
        if ((not includeMatches) and (nMatch == len(textRow))) or \
           (includeMatches and nMatch):
          
          appendText( textRow )
          if row < len(colorMatrix):
            appendColor(colorMatrix[row] )
          appendObject( object )

          if object in currentObjects:
            newObjs.append( object )
            newRange.append( len(scrolledMatrix.objectList)-1 )

      else:
        datum = formatText(textRow[col])
        match = find(string,datum,flag)

        if match:
          appendText( textRow )
          if row < len(colorMatrix):
            appendColor(colorMatrix[row] )
          appendObject( object )

          if object in currentObjects:
            newObjs.append( object )
            newRange.append( len(scrolledMatrix.objectList)-1 )

    if len(self.objectList) != len(scrolledMatrix.objectList):
      if col == '[All]':
        for i in range(len(scrolledMatrix.filteredCols)):
          scrolledMatrix.filteredCols[i] = 1
      else:
        for i in range(len(scrolledMatrix.filteredCols)):
          scrolledMatrix.filteredCols[i] = 0
        scrolledMatrix.filteredCols[col] = 1
    
    scrolledMatrix.objectSet = set(scrolledMatrix.objectList)     
    scrolledMatrix.currentObjects = newObjs
    scrolledMatrix.currentRange   = newRange
    scrolledMatrix.keyPressEscape()
    scrolledMatrix.refreshScrollbars()
    scrolledMatrix.drawCanvas()

  def selectMode(self, mode):

    self.mode = mode

  def filteredObjs(self, objsToBeFiltered):

    self.objsToBeFiltered = objsToBeFiltered

  def getColumnNames(self):
  
    columnList = ['[All]']
    for name in self.scrolledMatrix.headingList:
      name = re.sub('\s+',' ',name)
      columnList.append( name )
    
    return columnList

class FilterMatrix:

  def __init__(self, parent, scrolledMatrix, searchString,
               column=None, mode=SEARCH_MODES[0], selectedObjectsOnly=False):

    self.scrolledMatrix = scrolledMatrix
    self.searchString   = searchString
    
    self.column      = column #  None or a number - If None all cols are filtered
    self.mode        = mode   # 'literal', 'regexp' or if anything else case sensitive literal
    self.selected    = selectedObjectsOnly
    self.textMatrix  = self.scrolledMatrix.textMatrix
    self.colorMatrix = self.scrolledMatrix.colorMatrix
    self.objectList  = self.scrolledMatrix.objectList
   
    
  def remove(self):
  
    for i in range(len(self.scrolledMatrix.filteredCols)):
      self.scrolledMatrix.filteredCols[i] = 0
    
    self.scrolledMatrix.update(textMatrix  = self.textMatrix,
                               colorMatrix = self.colorMatrix,
                               objectList  = self.objectList)
    
  def apply(self, includeMatches=True):
  
    col    = self.column
    string = self.searchString
    flag   = re.S
    
    if includeMatches:
      find = re.search
    else:
      find = lambda: not re.search 
    
    if self.mode != SEARCH_MODES[2]:
      string  = re.escape(string)
    
    if self.mode == SEARCH_MODES[0]:
      flag = re.I
    
    scrolledMatrix = self.scrolledMatrix
    textMatrix     = self.textMatrix
    colorMatrix    = self.colorMatrix
    
    scrolledMatrix.textMatrix  = []
    scrolledMatrix.colorMatrix = []
    scrolledMatrix.objectList  = []
    
    appendText = scrolledMatrix.textMatrix.append
    appendObject = scrolledMatrix.objectList.append
    appendColor = scrolledMatrix.colorMatrix.append
    
    currentObjects = set(scrolledMatrix.currentObjects)
    
    newObjs  = []
    newRange = []
    for row in range(len(textMatrix)):
      object = self.objectList[row]
      if self.selected:
        if object not in currentObjects:
          continue
          
      if col is None:
        for col2 in range(len(textMatrix[row])):
          datum = textMatrix[row][col2]
          match = find(string,str(datum),flag)
          if match:
            appendText( textMatrix[row] )
            if row < len(colorMatrix):
              appendColor(colorMatrix[row] )
            appendObject( object )
            if object in currentObjects:
              newObjs.append( object )
              newRange.append( len(scrolledMatrix.objectList)-1 )
            break
      else:
        datum = textMatrix[row][col]
        match = find(string,str(datum),flag)
        if match:
          appendText( textMatrix[row] )
          if row < len(colorMatrix):
            appendColor(colorMatrix[row] )
          appendObject( object )
          if object in currentObjects:
            newObjs.append( object )
            newRange.append( len(scrolledMatrix.objectList)-1 )

    if len(self.objectList) != len(scrolledMatrix.objectList):
      if col is None:
        for i in range(len(scrolledMatrix.filteredCols)):
          scrolledMatrix.filteredCols[i] = 1
      else:
        for i in range(len(scrolledMatrix.filteredCols)):
          scrolledMatrix.filteredCols[i] = 0
        scrolledMatrix.filteredCols[col] = 1
          
    scrolledMatrix.objectSet = set(scrolledMatrix.objectList)
    scrolledMatrix.currentObjects = newObjs
    scrolledMatrix.currentRange   = newRange
    scrolledMatrix.keyPressEscape()
    scrolledMatrix.refreshScrollbars()
    scrolledMatrix.drawCanvas()

if __name__ == '__main__':

  root = Tkinter.Tk()
  root.grid_rowconfigure(0, weight=1)
  root.grid_columnconfigure(0, weight=1)
  root.geometry('600x200')
  
  headingList = ['Col A','Col B','Col C','Col D']
  
  table = ScrolledMatrix(root,headingList=headingList,multiSelect=True)
  table.grid(row=0, column=0, sticky='nsew')
  
  textMatrix = []
  objectList = []

  textMatrix.append([1,'One',1.00,u'\u03B1'])
  objectList.append(1)
  
  textMatrix.append([2,'Two',4.00,u'\u03B2'])
  objectList.append(2)
  
  textMatrix.append([3,'Three',9.00, u'\u03B3'])
  objectList.append(3)

  textMatrix.append([4,'Four',16.00, u'\u03B4'])
  objectList.append(4)

  textMatrix.append([5,'Five',25.00, u'\u03B5'])
  objectList.append(5)

  textMatrix.append([6,'Six',36.00,  u'\u03B6'])
  objectList.append(6)

  textMatrix.append([7,'Seven',49.00,u'\u03B7'])
  objectList.append(7)

  textMatrix.append([8,'Eight',64.00,u'\u03B8'])
  objectList.append(8)

  textMatrix.append([9,'Nine',81.00,u'\u03B9'])
  objectList.append(9)

  textMatrix.append([10,'Ten',100.00,u'\u03BA'])
  objectList.append(10)


  tipTexts = ['Tip A','Tip B','Tip C','Tip D']
  table.update(objectList=objectList, textMatrix=textMatrix, tipTexts=tipTexts)

  root.mainloop()

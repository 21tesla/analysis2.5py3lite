import Tkinter

from math import sqrt

from memops.gui.Canvas import Canvas
from memops.gui.Menu import Menu
from memops.gui.Frame import Frame
from memops.gui.MessageReporter import showOkCancel
from memops.gui import Color
from memops.universal.Util import buttonPress


CARDINAL_ANCHORS = ('N','E','S','W')

ANCHORS_GRIDDING = {'N':(0,0,1,0,'ew'),'E':(0,1,0,0,'ns'),
                    'S':(1,0,0,0,'ew'),'W':(0,0,0,1,'ns'),}

SHAPES = ('CIRCLE','SQUARE','DIAMOND','HEXAGON','OCTAGON')  
CIRCLE, SQUARE, DIAMOND, HEXAGON, OCTAGON = SHAPES

MOVE_STEP = 10
CORNER_SIZE = 8

ACTIVE = 'normal'
INACTIVE = 'disabled'


# Interacting with a data model: init(), delete()
# - could pass in initFunc, deleteFunc
# - could use subclasses
# - both nodes and links

# Persistence - Data Model XML? Regenerate state.

# Little 'E' for editable

# Adaptive rectangle

# Rounded corners

# Arrows. Dotted lines.

# Lock & unlock

# Line type change - merge 1-to-many

# Shadow

# Fixed areas and primatives - Canvas but easy 
# - Zones, Text box

# Cut-and-paste

# Unified colors for links

# Collision aviodance?

# Auto-arrange

# Text choppin...

# Zoom?

# Fame superclass for selection panel


from memops.universal.Util import isWindowsOS, useWheelMouse, isMacOS

class LinkChart(Frame):

  def __init__(self, parent, scrollX=True, scrollY=True, *args, **kw):
  
    self.parent = parent
    self.scrollX = scrollX
    self.scrollY = scrollY
    
    self.prototypes = []
    self.nodes = []
    self.node = None # Selected
    self.links = []
    self.link = None # Selected
    self.objDict = {} # Map canvas to nodes, links
    self.cornerDict = {}
    self.editWidget = None
    self.mouseOver = None
    self.inMotion = None
    self.inMotion2 = None
    self.inMotion3 = None
    self.dragRegion = None
    self._resizeBusy = False

    self.selectedNodes = []
    self.selectedLinks = []
    self.highlightNodes = []
    self.highlightLinks = []
 
    kw['relief'] = 'flat'
    kw['borderwidth'] = 0
    
    Frame.__init__(self, parent, *args, **kw)

    self.canvas = canvas = Canvas(self, relief='flat', borderwidth=0)
    self.canvas.configure(xscrollincrement=2, yscrollincrement=2, bg=self.cget('bg'))
    self.canvas.grid(row=0, column=0, sticky='nsew')
    
    if scrollX:
      self.horizScrollbar = Tkinter.Scrollbar(self, bg=self.cget('bg'),
                                              command=self.canvas.xview,
                                              orient=Tkinter.HORIZONTAL,
                                              borderwidth=1)
      self.canvas.configure(xscrollcommand=self.horizScrollbar.set)
    
    if scrollY:  
      self.vertScrollbar  = Tkinter.Scrollbar(self, bg=self.cget('bg'),
                                              command=self.canvas.yview,
                                              orient=Tkinter.VERTICAL,
                                              borderwidth=1)
      self.canvas.configure(yscrollcommand=self.vertScrollbar.set)
 
    canvas.bind('<Button-1>', self._mouseSingleClick)

    if useWheelMouse():
      canvas.bind('<MouseWheel>', self._windowsOsScroll)
    else:
      canvas.bind('<Button-4>', self.scrollUp)
      canvas.bind('<Button-5>', self.scrollDown)

    canvas.bind('<Double-1>', self._mouseDoubleClick)
    canvas.bind('<B1-Motion>', self._mouseDrag)
    canvas.bind('<ButtonRelease-1>', self._mouseRelease)                                 
    canvas.bind('<Motion>', self._mouseEnter)
    canvas.bind('<Enter>', self._mouseEnter)
    canvas.bind('<Leave>', self._mouseLeave)
    canvas.bind('<KeyPress-Up>', self._moveUp)
    canvas.bind('<KeyPress-Down>', self._moveDown)
    canvas.bind('<KeyPress-Right>', self._moveRight)
    canvas.bind('<KeyPress-Left>', self._moveLeft)
    canvas.bind('<KeyPress-Delete>',self.deleteSelected)
    self.bind('<Configure>', self._resizeAfter)
    
    self.menu = Menu(self, tearoff=False)
    
    alignItems = [{ 'kind': 'command', 'label':'Y axis','command': self._alignHorizontal },
                  { 'kind': 'command', 'label':'X axis','command': self._alignVertical }]

    distributeItems = [{ 'kind': 'command', 'label':'Horizontal','command': self._distributeHorizontal },
                       { 'kind': 'command', 'label':'Verical','command': self._distributeVertical }]
    
    menuItems = [
                 {'kind': 'cascade', 'label': 'Align',
                  'submenu': alignItems }, #0
                 {'kind': 'cascade', 'label': 'Distribute',
                  'submenu': distributeItems }, #1
                 {'kind': 'command', 'label': 'Link Selected', 
                  'command': self.linkSelected },#2
                 {'kind': 'command', 'label': 'Reset Links', 
                  'command': self._resetLinks },#3
                 {'kind': 'separator'}, 
                 {'kind': 'command', 'label': 'Delete Links', #5
                  'command': self.deleteLinks },
                 {'kind': 'command', 'label': 'Delete Items', #6
                  'command': self.deleteNodes}, 
                 {'kind': 'separator'}, 
                 {'kind': 'command', 'label': 'Save PostScript File', 
                  'command': self.printCanvas},
                 ]
                 
    self.menuItems = menuItems              
    self.menu.setMenuItems(menuItems)
    
    # canvas.bind('<ButtonPress-3>', self._popupMenu)
    canvas.bind(buttonPress(3), self._popupMenu)

    self._drawAfter()

  def _resizeAfter(self, event):
  
    if self._resizeBusy:
      return
          
    self._resizeBusy = True
    self._popdownMenu(event)
    self.after_idle( lambda:self._refreshScrollbars(event) )

  def _refreshScrollbars(self, event):
  
    bbox = self.canvas.bbox('all')
    
    if not bbox:
      self._resizeBusy = False
      return
    
    cWidth  = event.width
    cHeight = event.height
    
    vertScrollbar = self.vertScrollbar
    horizScrollbar = self.horizScrollbar
    
    xSpace = bbox[2]+bbox[0] # Symmetric space
    ySpace = bbox[3]+bbox[1]
    
    horizHeight = int(horizScrollbar.winfo_height())
    vertWidth = int(vertScrollbar.winfo_width())
    
    if cHeight < ySpace+horizHeight:
      vertScrollbar.grid(row=0, column=1, sticky='ns')
      cWidth -= vertWidth
    else:
      vertScrollbar.grid_forget()

    if cWidth < xSpace+vertWidth:
      horizScrollbar.grid(row=1, column=0, sticky='ew')
      cHeight -= horizHeight
    else:
      horizScrollbar.grid_forget()
    
    x1 = 0
    y1 = 0
    x2 = max(cWidth, xSpace)
    y2 = max(cHeight, ySpace)
        
    self.canvas.configure(width=cWidth,height=cHeight)
    self.canvas.configure( scrollregion = (x1,y1,x2,y2) )
    self.update_idletasks()
    self._resizeBusy = False

  def _popupMenu(self, event):
  
    self._finishEdit()
    self._configMenu()
    self.menu.popupMenu(event)

  def _popdownMenu(self, event):
  
    self._finishEdit()
    self.menu.unpost()
  
  def _configMenu(self):
  

    if self.selectedLinks:
      self.menu.entryconfig(5, state=ACTIVE)
      self.menu.entryconfig(3, state=ACTIVE)
    else:
      self.menu.entryconfig(5, state=INACTIVE)
      self.menu.entryconfig(3, state=INACTIVE)
        
    
    if self.selectedNodes:
      self.menu.entryconfig(2, state=ACTIVE)
      self.menu.entryconfig(6, state=ACTIVE)
    else:
      self.menu.entryconfig(2, state=INACTIVE)
      self.menu.entryconfig(6, state=INACTIVE)
         
    if self.selectedNodes or self.selectedLinks:
      self.menu.entryconfig(0, state=ACTIVE)
      self.menu.entryconfig(1, state=ACTIVE)
    else:
      self.menu.entryconfig(0, state=INACTIVE)
      self.menu.entryconfig(1, state=INACTIVE)
  
  def _resetLinks(self):
  
    for link in self.selectedLinks:
      link.resetCoords()
  
  def _alignSelected(self, vertical=False, horizontal=False):
  
    obj = self.node or self.link
    
    if not obj:
      return 
  
    x0, y0 = obj.coords  
    
    objs = self.selectedNodes + self.selectedLinks
    n = len(objs)
    
    if n < 2:
      return

    if horizontal:
      for obj in objs:
        obj.move(obj.coords[0],y0)
    
    if vertical:
      for obj in objs:
        obj.move(x0,obj.coords[1])
    

  def _distributeSelected(self, vertical=False, horizontal=False):
  
    xCoords = []
    yCoords = []
    
    objs = self.selectedNodes + self.selectedLinks
    n = len(objs)
    
    if n < 3:
      return
    
    for obj in objs:
      x, y = obj.coords
      xCoords.append((x, obj))
      yCoords.append((y, obj))
    
    if horizontal:
      xCoords.sort()
      xA = xCoords[0][0]
      dx = (xCoords[-1][0]-xA)/float(n-1)
      
      for i, (x0, obj) in enumerate(xCoords):
        y = obj.coords[1]
        x = xA + (i*dx)
        obj.move(x,y)
    
    if vertical:
      yCoords.sort()
      yA = yCoords[0][0]
      dy = (yCoords[-1][0]-yA)/float(n-1)

      for i, (y0, obj) in enumerate(yCoords):
        x = obj.coords[0]
        y = yA + (i*dy)
        obj.move(x,y)
        
  def _alignHorizontal(self):
  
    self._alignSelected(horizontal=True)
  
  def _alignVertical(self):

    self._alignSelected(vertical=True)

  def _distributeHorizontal(self):
  
    self._distributeSelected(horizontal=True)
  
  def _distributeVertical(self):
  
    self._distributeSelected(vertical=True)
  
  def linkSelected(self):
  
    nodes = self.selectedNodes
    n = len(nodes)
    if n > 1:
      
      for i in range(n-1):
        nodeA = nodes[i]
      
        for j in range(i+1,n):
          nodeB = nodes[j]
        
          for link in nodeA.links:
            if nodeB in link.nodes:
              #print link, link.nodes, nodeA, nodeB
              break
        
          else:
            link = Link(self, nodeA, nodeB, text='*')
            #print link
  
  def deleteNodes(self):
  
    if self.selectedNodes:
      msg = 'Really delete %d selected nodes?' % len(self.selectedNodes)
      if showOkCancel('Confirm', msg, parent=self):
        for node in self.selectedNodes[:]:
          node.delete()
        self.selectedNodes = []
        self.node = None 

 
  def deleteLinks(self):
    
    if self.selectedLinks:
      msg = 'Really delete %d selected links?' % len(self.selectedLinks)
      if showOkCancel('Confirm', msg, parent=self):
        for link in self.selectedLinks[:]:
          link.delete()
                
        self.selectedLinks = []
        self.link = None

  def deleteSelected(self, event=None):
    
    objs = self.selectedLinks + self.selectedNodes
    if objs:
      n1 = len(self.selectedNodes)
      n2 = len(self.selectedLinks)
      msg = 'Really delete %d selected objects?\n' % (n1+n2)
      msg += '(%d nodes, %d links)' % (n1, n2)
      if showOkCancel('Confirm', msg, parent=self):
        for obj in objs:
           obj.delete()
                
        self.selectedLinks = []
        self.link = None
        self.selectedNodes = []
        self.node = None 
 
       
  def printCanvas(self, *event):
  
    from memops.gui.FileSelect import FileType
    from memops.gui.FileSelectPopup import FileSelectPopup

    fileTypes = [  FileType('PostScript', ['*.ps']), FileType('All', ['*'])]
    fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
               title = 'Print canvas to file', dismiss_text = 'Cancel',
               selected_file_must_exist = False)

    fileName = fileSelectPopup.getFile()
        
    self.postscript(colormode='color',file=fileName)

  def _clearSelection(self):
  
    for node in self.selectedNodes:
      node.isSelected = False

    for link in self.selectedLinks:
      link.isSelected = False
    
    self.selectedNodes = []
    self.selectedLinks = []    
    self.link = None
    self.node = None
    
  def _clearHighlights(self):

    for node in self.highlightNodes:
      node.highlight = False
      
    for link in self.highlightLinks:
      link.highlight = False
    
    self.mouseOver = None
    self.highlightNodes = []
    self.highlightLinks = []
    
  def _mouseSingleClick(self, event):
 
    self._mouseClick(event)

  def _mouseDoubleClick(self, event):
 
    self._mouseClick(event, True)
  
  def _mouseClick(self, event, edit=False):
 
 
    #self.menu.popdownMenu(event)
 
    self._finishEdit()

    canvas = event.widget
    canvas.focus_force()
  
    x = canvas.canvasx(event.x)
    y = canvas.canvasy(event.y)
    items = canvas.find('overlapping',x-1,y-1,x+1,y+1)
    
     
    #item = self.find('closest',x,y)[0]
    if not items:
      self._clearSelection()
      self._clearHighlights()
      self.draw()
      return
    
    objs = [self.objDict.get(item) for item in items if self.objDict.get(item)] 
    
    for obj in objs:
      if isinstance(obj, Node):        
        if not (event.state & 4 or event.state & 1):
          self._clearSelection()
   
        obj.isSelected = True
        self.node = obj
        if obj not in self.selectedNodes:
          self.selectedNodes.append(obj)
          
        break 
    
    else:
      self.node = None
      
      for obj in objs:
        if isinstance(obj, Link):
          if not (event.state & 4 or event.state & 1):
            self._clearSelection()
           
          obj.isSelected = True
          self.link = obj
          
          if obj not in self.selectedLinks:
            self.selectedLinks.append(obj)
          
          break
          
      else:
        self.link = None
    
      
    self._drawAfter()
    if not obj:
      return
    
    if edit and obj.editWidget:
      obj.edit()
      self.editWidget = obj.editWidget
  
  def _clearDragRegion(self, canvas):
  
    if self.dragRegion:
      canvas.delete(self.dragRegion[0])
      self.dragRegion = None
      
  def _mouseEnter(self, event):
    
    self.menu.popdownMenu(event)
    canvas = event.widget
    canvas.focus_force()
    
    x = canvas.canvasx(event.x)
    y = canvas.canvasy(event.y)
    items = canvas.find('overlapping',x-1,y-1,x+1,y+1)
    
    #item = self.find('closest',x,y)[0]
    if not items:
      self._clearHighlights()
      self.draw()
      return
    
    obj = self.objDict.get(items[0])

    self._clearHighlights()
    
    if obj is None:
      self.draw()
      return    
      
    if obj is self.mouseOver:
      if self.inMotion2:
        for item in items:
          if self.cornerDict.has_key(item):
            self.inMotion2 = None
            break 
            
      return 
    
    objs = [self.objDict.get(item) for item in items if self.objDict.get(item)] 
    
    for obj in objs:
      if isinstance(obj, Node):
        self.highlightNodes = [obj,]
        break 
    
    else:
      for obj in objs:
        if isinstance(obj, Link):
          self.highlightLinks = [obj,]
    
    obj.highlight = True
    
    self.mouseOver = obj
    self.draw()
    
  def _mouseLeave(self, event):
    
    #self._mouseRelease(event)    
    #self._finishEdit()
    self._clearHighlights()
    self.draw()
  
  def _mouseDrag(self, event):
 
    self.menu.popdownMenu(event)
 
    canvas = event.widget
    canvas.focus_force()
    x = canvas.canvasx(event.x)
    y = canvas.canvasy(event.y)
    
    #width = int(canvas.winfo_width())
    #height = int (canvas.winfo_height())
    #
    #if not (0 < x < width) or not (0 < y < height):
    #  self._mouseRelease(event)
    #  return
    
    inMotion = self.inMotion
    
    objs = self.selectedNodes + self.selectedLinks
    if objs:
      items = canvas.find('overlapping',x-1,y-1,x+1,y+1)
      
      if not inMotion:
        for item in items:
          if self.inMotion2 or self.cornerDict.has_key(item):
            obj = self.node
 
            x1, y1 = obj.coords # offsets
            dx = abs(x-x1)
            dy = abs(y-y1)
 
            size = max(dx, dy)+1

            obj.width = size
            obj.draw()
 
            self.inMotion2 = True
 
            return
      
      if inMotion:
        offsets = inMotion # Get offsets

      else:
        offsets = []
        for obj in objs:
          x1, y1 = obj.coords # offsets
          dx = x-x1 # Offset relative to centre
          dy = y-y1
          
          if isinstance(obj, Node) and obj.isPrototype:
            obj = obj.copy()
          
          offsets.append((dx,dy,obj))
 
        self.inMotion = offsets  # Store click offset and Node

      #print event, event.widget._root()

      for dx, dy, obj in offsets:
        x2 = x-dx
        y2 = y-dy

        obj.move(x2,y2)

      x0, y0 = event.x_root, event.y_root
      under = self.winfo_containing(x0, y0)
      
      if isinstance(under, Canvas):
        under = under.parent
      
      if isinstance(under, LinkChart) and (under is not self):
        for dx, dy, obj in offsets:
          if not isinstance(obj, Node):
            continue
        
          original = obj.original
          if original and original.isPrototype:
            xT = x0-under.winfo_rootx()-dx
            yT = y0-under.winfo_rooty()-dy
            
            if original.spawn:
              original.spawn.move(xT, yT)
            else:
              original.spawn = original.copy(under, coords=(xT,yT))

      else:
        for dx, dy, obj in offsets:
          if not isinstance(obj, Node):
            continue

          original = obj.original
          if original and original.spawn:
            original.spawn.delete()
            original.spawn = None

    else:
      if self.inMotion3:
        x0, y0, = self.inMotion3

      else:
        #self._mouseRelease(event)
        x0 = x
        y0 = y
        self.inMotion3 = (x, y)
        
      if self.dragRegion:
        item = self.dragRegion[0]
        canvas.coords(item,x0,y0,x,y)
      else:
        item = canvas.create_rectangle(x0,y0,x,y,outline='#008000',
                                       fill='#00E000', stipple='gray25')
      
      self.dragRegion = (item,x0,y0,x,y)
  
  def _mouseRelease(self, event):

    canvas = event.widget
    
    if self.dragRegion:
      item, x1, y1, x2, y2 = self.dragRegion
      
      xCoords = [x1, x2]
      yCoords = [y1, y2]
      xCoords.sort()
      yCoords.sort()
      
      x1, x2 = xCoords
      y1, y2 = yCoords
      
      self.node = None
      self.selectedNodes = []
      for node in self.nodes:
        x0, y0 = node.coords
        
        if (x1 < x0 < x2) and (y1 < y0 < y2):
          self.selectedNodes.append(node)
          node.isSelected = True
          self.node = node
        else:
          node.isSelected = False

      self.link = None
      self.selectedLinks = []
      for link in self.links:
        x0, y0 = link.coords
        
        if (x1 < x0 < x2) and (y1 < y0 < y2):
          self.selectedLinks.append(link)
          link.isSelected = True
          self.link = link
        else:
          link.isSelected = False

      self._clearDragRegion(canvas)
      
      self.draw()

    if self.inMotion:
      for dx, dy, obj in self.inMotion:
        if isinstance(obj, Node) and obj.original and obj.original.isPrototype:
          if obj.original.spawn:
            obj.original.spawn.original = None 
            obj.original.spawn = None
          obj.delete()

    self.inMotion = None
    self.inMotion2 = None
    self.inMotion3 = None
 
  def _finishEdit(self):
    
    if self.editWidget:
      self.focus_force()
      self.editWidget.commit()
      self.editWidget.hide()
      self.editWidget = None

  def _moveUp(self, event):
  
    self._moveSelected(event, 0, -1)
  
  def _moveDown(self, event):
  
    self._moveSelected(event, 0, 1)
  
  def _moveRight(self, event):
  
    self._moveSelected(event, 1, 0)
  
  def _moveLeft(self, event):
  
    self._moveSelected(event, -1, 0)

  def _moveSelected(self, event, dx, dy):
    
    if not (event.state & 4 or event.state & 1):
      dx *= MOVE_STEP
      dy *= MOVE_STEP  
  
    for node in self.selectedNodes:
      x0, y0 = node.coords
      node.move(x0+dx, y0+dy)
      
    for link in self.selectedLinks:
      x0, y0 = link.coords
      link.move(x0+dx, y0+dy)

  def scrollUp(self, event=None):
    
    pass

  def scrollDown(self, event):
    
    pass

  def _windowsOsScroll(self, event=None):
  
    delta = event.delta
    
    if delta > 0:
      self.scrollUp(event)
      
    elif delta < 0:
      self.scrollDown(event)
  
  def _drawAfter(self):
  
    self.after_idle(self.draw)
  
  def draw(self):
  
    for link in self.links:
      link.draw()
  
    for node in self.nodes:
      node.draw()
       
      
class Node:

  def __init__(self, linkChart, title='', updateFunc=None, width=150,
               isSelected=False, shape=SQUARE, coords=None, object=None,
               color='#8080F0',  editWidget=None,
               highlight=False, resizeable=False, original=None,
               titleFont='Helvetica 10', mainFont='Helvetica 10 italic'):
  
    self.parent = linkChart
    self.object = object
    self.title = title
    self.updateFunc = updateFunc
    self.width = width/2.0
    self.color = color
    self.shape = shape
    self.isPrototype = False
    self.coords = coords or (0,0) # Centre
    self.canvasItems = []
    self.links = []
    self.editWidget = editWidget
    self.isSelected = isSelected
    self.highlight = highlight
    self.resizeable = resizeable
    self.original = original
    self.titleFont = titleFont
    self.mainFont = mainFont
    self.text = None
    
    self._updateText()
    
    if isSelected:
      linkChart.selectedNodes.append(self)
      
    linkChart.nodes.append(self)

  def _updateText(self):
  
    if self.updateFunc:
      self.text = self.updateFunc(self.object)
    else:
      self.text = ''  

    self.draw()

  def copy(self, linkChart=None, coords=None):
  
    if not linkChart:
      linkChart = self.parent
  
    if not coords:
      coords = self.coords
  
    copy = Node(linkChart, self.title, self.updateFunc,
                self.width*2.0, self.isSelected, self.shape,
                coords, self.object, self.color, 
                self.editWidget, self.highlight, self.resizeable, 
                self, self.titleFont, self.mainFont)
    
    return copy

  def delete(self):
   
    parent = self.parent
    canvas = parent.canvas
  
    if self.editWidget:
      self.editWidget.hide()
   
    for link in list(self.links):
      link.delete()
   
    if self in parent.selectedNodes:
      parent.selectedNodes.remove(self)
   
    parent.nodes.remove(self)
    
    objDict = parent.objDict
    cornerDict = parent.cornerDict
    
    for item in self.canvasItems:
      del objDict[item]
      canvas.delete(item)
      
      if cornerDict.get(item):
        del cornerDict[item]
      
    del self
   
  def move(self, x1, y1):
  
  
    if self.editWidget:
      self.editWidget.hide()

    if self.original and self.original.isPrototype:
      x, y = x1, y1
    else:
      width = int(self.parent.winfo_width())
      height = int(self.parent.winfo_height())
      x = max(0, min(width, x1))
      y = max(0, min(height, y1))

    self.coords = x, y
    
    for link in self.links:
      link.draw()
    
    self.draw()

  def _getBoxCoords(self, cx, cy):
  
    w = self.width
    return cx-w, cy-w, cx+w, cy+w

  def _getOutlineCoords(self, thickness, width, cx, cy):
    
    w = width
    
    if self.shape == CIRCLE:
      d = sqrt(2*w*w)
    
      x3 = cx-d
      y3 = cy-d
      x4 = cx+d
      y4 = cy+d
    
      p = 5
      outlineCoords = x3-p, y3-p, x4+p, y4+p
      
      p -= thickness
      inlineCoords = x3-p, y3-p, x4+p, y4+p
    
    elif self.shape == SQUARE:
      x1, y1, x2, y2 = self._getBoxCoords(cx, cy)
      
      p = 12
      x1 -= p
      x2 += p 
      y1 -= p
      y2 += p
      
      outlineCoords = x1, y1, x2, y2
      
      t = thickness
      inlineCoords = x1+t, y1+t, x2-t, y2-t
  
    elif self.shape == DIAMOND:
      x1 = cx
      x2 = cx+(w*2)
      x3 = cx
      x4 = cx-(w*2)
      y1 = cy-(w*2)
      y2 = cy
      y3 = cy+(w*2)
      y4 = cy
      
      p = 5
      outlineCoords = x1, y1-p, x2+p, y2, x3, y3+p, x4-p, y4

      p -= thickness
      inlineCoords = x1, y1-p, x2+p, y2, x3, y3+p, x4-p, y4
  
    elif self.shape == HEXAGON:
      w2 = w+8
  
      dx = 1.5*w2*0.8660254037844386
      dy = 1.5*w2
      dy2 = dy/2.0
      
      x1 = cx
      x2 = cx+dx
      x3 = cx+dx
      x4 = cx
      x5 = cx-dx
      x6 = cx-dx
      y1 = cy-dy
      y2 = cy-dy2
      y3 = cy+dy2
      y4 = cy+dy
      y5 = cy+dy2
      y6 = cy-dy2
      
      outlineCoords = x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6

      p = thickness
      inlineCoords = x1, y1+p, x2-p, y2, x3-p, y3, x4, y4-p, x5+p, y5, x6+p, y6
      
    elif self.shape == OCTAGON:
      w2 = w+8
      d = 1.41421356237309503*w2
      
      x1 = cx-d
      x2 = cx-w2
      x3 = cx
      x4 = cx+w2
      x5 = cx+d
      y1 = cy-d
      y2 = cy-w2
      y3 = cy
      y4 = cy+w2
      y5 = cy+d
  
      outlineCoords = (x3, y1, x4, y2, x5, y3, x4, y4,
                       x3, y5, x2, y4, x1, y3, x2, y2)
      
      p = thickness
      inlineCoords = (x3, y1+p, x4-p, y2+p, x5-p, y3, x4-p, y4-p,
                      x3, y5-p, x2+p, y4-p, x1+p, y3, x2+p, y2+p)
      
    return outlineCoords, inlineCoords

  def draw(self, coords=None):
  
    if coords:
      cx, cy = coords
    else:
      cx, cy = self.coords
    
    color = self.color
    canvas = self.parent.canvas
    cCoords = canvas.coords
    cLift = canvas.lift

    rgb = Color.hexToRgb(color)
    
    if self.isSelected:
      fillColor = Color.inverseRgb(color)
      thickness = 1
      
    elif self.highlight:
      fillColor = Color.scaleColor(canvas, color, 1.2)
      thickness = 1
      
    else:
      fillColor = color
      thickness = 1

    textColor = Color.inverseGrey(fillColor)
    colorA = Color.scaleColor(canvas, color, 0.40)
    colorB = Color.scaleColor(canvas, color, 1.4)
    colorC = Color.scaleColor(canvas, fillColor, 0.90)

    
    outlineCoords, inlineCoords = self._getOutlineCoords(thickness, self.width,
                                                         cx, cy)
    
    x1, y1, x2, y2 = self._getBoxCoords(cx, cy)
    
    if self.resizeable:
      c = CORNER_SIZE
    else:
      c = 0  
    
    if self.canvasItems:
    
      cConfig = canvas.itemconfigure
    
      outline, inline, box, corner, title, text = self.canvasItems
      
      cConfig(outline, outline=colorA, fill=colorA, width=thickness)
      cConfig(inline, outline=colorB, fill=color, width=thickness)
      cConfig(box, outline=colorC, fill=fillColor, width=thickness)
      cConfig(corner, outline=colorC, fill=colorC)
      cConfig(title, text=self.title, fill=textColor, font=self.titleFont)
      cConfig(text, text=self.text, fill=textColor, font=self.mainFont)

      cCoords(outline, *outlineCoords)
      cCoords(inline, *inlineCoords)
      cCoords(box, x1, y1, x2, y2)
      cCoords(corner, x2, y2-c, x2, y2, x2-c, y2)
      
    else:
      
      if self.shape == CIRCLE:
        create = canvas.create_oval
        
      elif self.shape == SQUARE:
        create = canvas.create_rectangle
        
      else:
        create = canvas.create_polygon
      
      outline = create(width=thickness, outline=colorA,
                       fill=colorA, *outlineCoords)

      inline = create(width=thickness, outline=colorB,
                      fill=color, *inlineCoords)

      box = canvas.create_rectangle(x1, y1, x2, y2, width=thickness,
                                    outline=colorC, fill=fillColor)
                                    
      corner = canvas.create_polygon(x2, y2-c, x2, y2, x2-c, y2,
                                     outline=colorC, fill=colorC)

      title = canvas.create_text(0, 0, text=self.title, font=self.titleFont,
                                 fill=textColor, anchor='nw')
      
      text = canvas.create_text(0, 0, text=self.text, font=self.mainFont,
                                fill=textColor, anchor='nw')

      self.canvasItems = outline, inline, box, corner, title, text

    t = 2.0*thickness
    xA, yA, xB, yB = canvas.bbox(title)
    cCoords(title, x1+t, y1+t)
    cCoords(text, x1+t, y1+t+t+(yB-yA))
      
    cLift(outline)  
    cLift(inline)  
    cLift(box)  
    cLift(corner)
    cLift(title)  
    cLift(text)  
       
    for item in self.canvasItems:
      self.parent.objDict[item] = self
    
    if self.resizeable:
      self.parent.cornerDict[corner] = self

  def select(self):  

    self.isSelected = True
    
    # Affect links?
    
  def edit(self):
  
    if self.editWidget:
      self.editWidget.update(self.object)
      x, y = self.coords # Fix this
      self.editWidget.display(self, x, y)

class PrototypeNode(Node):

  def __init__(self, *args, **kw):
  
    Node.__init__(self, *args, **kw)
    self.isPrototype= True
    self.spawn = None
  
  def move(self, x, y):
    pass
    
  def delete(self):
    pass

  def edit(self):
    pass
  

class Link:

  def __init__(self, linkChart, nodeA, nodeB, text='',
               color='#000000', object=None, editWidget=None,
               isSelected=False, highlight=False, coords=None,
               textOffset=None, curved=True):
  
    self.parent = linkChart
    self.nodes = [nodeA, nodeB]
    self.nodeA = nodeA
    self.nodeB = nodeB
    self.text = text
    self.color = color
    self.object = object
    self.editWidget = editWidget
    self.coords = coords
    self.textOffset = textOffset or (10, 10)
    self.canvasItems = []
    self.isSelected = isSelected
    self.highlight = highlight
    self.curved = curved

    if not coords:
      self.resetCoords()
    
    if isSelected:
      linkChart.selectedLinks.append(self)
    
    linkChart.links.append(self)
    nodeA.links.append(self)  
    nodeB.links.append(self)  
  
  def resetCoords(self):
  
    x1, y1 = self.nodeA.coords
    x3, y3 = self.nodeB.coords
    x2, y2 = (x1+x3)/2.0, (y1+y3)/2.0

    self.coords = x2, y2
    self.textOffset = 10, 10 
    
  def delete(self):
  
    if self.editWidget:
      self.editWidget.hide()
      
    for node in self.nodes:
      node.links.remove(self)

    if self in self.parent.selectedLinks:
      self.parent.selectedLinks.remove(self)

    self.parent.links.remove(self)
    
    for item in self.canvasItems:
      del self.parent.objDict[item]
      self.parent.canvas.delete(item)
      
    del self  
    
  def move(self, x1, y1):
  
    if self.nodeA.isPrototype or self.nodeB.isPrototype:
      return 
      
    if self.editWidget:
      self.editWidget.hide()
    
    width = int(self.parent.winfo_width())
    height = int(self.parent.winfo_height())
    x = max(0, min(width, x1))
    y = max(0, min(height, y1))
    
    self.coords = (x, y)
    self.draw()
    #self.nodeA.draw()
    #self.nodeB.draw()
    
  def draw(self):
    
    
    if self.isSelected:
      lineColor = '#800000'
      fillColor='#FFFF80'
      thickness=2
      
    elif self.highlight:
      lineColor = self.color
      fillColor='#8080B0'
      thickness=1
      
    else:
      lineColor = self.color
      fillColor = '#FFFFFF'
      thickness = 1
    
    canvas = self.parent.canvas
    t = thickness
    
    x1, y1 = self.nodeA.coords
    x0, y0 = self.coords
    x3, y3 = self.nodeB.coords
    xT, yT = self.textOffset


    if self.curved:
      smooth = 1
      xM, yM = (x1+x3)/2.0, (y1+y3)/2.0
      dx, dy = x0-xM, y0-yM
      x2, y2 = xM+(2*dx), yM+(2*dy)
    
    else:
      smooth = 0
      x2, y2 = x0, y0
    
    d = 3
     
    if self.canvasItems:
    
      line, box, text = self.canvasItems
      
      canvas.itemconfigure(line, fill=lineColor,
                           width=thickness, smooth=smooth)
      canvas.coords(line, x1, y1, x2, y2, x3, y3)
      
      canvas.itemconfigure(text, text=self.text, fill=lineColor)
      canvas.coords(text, x0, y0)
      
      xA, yA, xB, yB = canvas.bbox(text)
      canvas.itemconfigure(box, outline=lineColor, fill=fillColor)
      canvas.coords(box, xA, yA, xB, yB)
      
    else:
      line = canvas.create_line(x1, y1, x2, y2, x3, y3, smooth=smooth,
                                fill=lineColor, width=thickness, splines=25)
      text = canvas.create_text(x0, y0, text=self.text,
                                fill=lineColor, anchor='center')
      
      xA, yA, xB, yB = canvas.bbox(text)
      box = canvas.create_rectangle(xA, yA, xB, yB,
                                    outline=lineColor, fill=fillColor)
      
      self.canvasItems = [line, box, text]
    
    
    canvas.lower(text)
    canvas.lower(box)
    canvas.lower(line)
    
    for item in self.canvasItems:
      self.parent.objDict[item] = self
      
      
  def select(self):  

    self.isSelected = True
    
    self.draw()
    
  def edit(self):
  
    if self.editWidget:
      x,y = self.getCenter()
    
      self.editWidget.update(self.object)
      self.editWidget.display(x, y)

  def getCenter(self):
  
    # Straight line only at the mo...
    
    x1, y1, x2, y2 = self.coords

    return (x1+x2)/2.0, (y1+y2)/2.0


class EditWidget:

  def __init__(self, tkWidget, canvas, setupCallback=None, commitCallback=None):
  
    self.tkWidget = tkWidget
    self.canvas = canvas
    self.setupCallback = setupCallback
    self.commitCallback = commitCallback
    self.canvasItem = None
    self.object = None
    self.node = None
  
  def display(self, node, x, y):
  
    make = self.canvas.create_window
  
    self.node = node
    self.tkWidget.configure(borderwidth=1, relief='solid')
    self.canvasItem = make(x, y, window=self.tkWidget, anchor='center')    
    self.tkWidget.focus_force()
    self.tkWidget.bind('<KeyPress-Escape>', self.hide)
    self.canvas.lift(self.tkWidget)
  
  def hide(self):
  
    if self.canvasItem:
      self.canvas.delete(self.canvasItem)
      self.canvasItem = None
      self.node = None
  
  def update(self, object):
  
    self.object = object
    
    if self.setupCallback:
      self.setValue(self.setupCallback(object))
  
  def commit(self):
  
    if self.node and self.commitCallback:
      value = self.getValue()
      self.commitCallback(self.object, value)
      self.node._updateText()

  def getValue(self):
    # Overwrite for widgets without .get()
  
    return self.tkWidget.get()
    
  def setValue(self, value):
    # Overwrite for widgets without .set()
  
    self.tkWidget.set(value)
      
      
class PulldownEditWidget(EditWidget):

  def getValue(self):
  
    return self.tkWidget.getObject()
    
  def setValue(self, value):
  
    names, objects, index = value
    self.tkWidget.setup(names, objects, index)
    
    
if __name__ == '__main__':

      
  import Tkinter
  
  from memops.gui.Entry import Entry
  from memops.gui.LabelFrame import LabelFrame
  from memops.gui.ScrolledFrame import ScrolledFrame
  

  root = Tkinter.Tk()
  root.grid_columnconfigure(0, weight=1)
  root.grid_rowconfigure(1, weight=1)
  root.geometry('900x600')
  
  frame = LabelFrame(root, text='Prototypes', grid=(0,0))
  frame.expandGrid(0,0)
    
  linkChartA = LinkChart(frame, background='#8080B0', height=160)
  linkChartA.grid(row=0, column=0, sticky='nsew')

  linkChartB = LinkChart(root)
  linkChartB.grid(row=1, column=0, sticky='nsew')
 
  tkWidget = Entry(root, text='', width=8)
  
  def setupFunc(obj):
    # This function creates widget contents
    # from the underlying node object
    return obj[0]
  
  def commitFunc(obj, value):
    # This function takes a value from an
    # edit widget and applies it to the
    # underlying object
    obj[0] = value
  
  def updateFunc(obj):
    # This function defines a Node's main text
    # based upon any underlying object
    # This function would typically list attributes.
    return 'x = ' + obj[0]


  
  editWidget = EditWidget(tkWidget, linkChartB,
                          setupCallback=setupFunc,
                          commitCallback=commitFunc)
  
  nodes = []
  for x in range(len(SHAPES)):
    y = x+1
  
    if x % 2 == 0:
      resizeable = False
    else:
      resizeable = True
  
    node = Node(linkChartB, title='Node %d' % (y), width=60,
                updateFunc=updateFunc ,resizeable=resizeable,
                shape=SHAPES[x], coords=(y*160, 160),
                editWidget=editWidget, object=['A%d' % y])

    if nodes:
      prev = nodes[-1]
      link = Link(linkChartB, node, prev,
                  text='%s-%s' % (prev.title[-1], node.title[-1]))
     
    nodes.append(node)

    node = PrototypeNode(linkChartA, title='Type %d' % (y), width=60,
                         updateFunc=None, resizeable=False, shape=SHAPES[x],
                         coords=(60+x*120, 80), object=None)

    
  link = Link(linkChartB, nodes[0], nodes[2],
              text='%s-%s' % (1,3), coords=(320, 50))

  linkChartA.draw()
  linkChartB.draw()

  root.mainloop()

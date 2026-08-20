from sys import platform

from os import listdir, path

from Tkinter import PhotoImage

from memops.gui.Scrollbar import Scrollbar
from memops.gui.Canvas import Canvas
from memops.gui.Frame import Frame
from memops.gui.Menu import Menu

from memops.universal.Io import getTopDirectory
from memops.universal.Util import useWheelMouse, buttonPress

# TBD
# EditWidgets
# Right mouse menu: collapse/expand all; print
# Select between small/med/large icons

# isWindowsOS = platform[:3].lower() == 'win'

TOGGLE = 0
ICON = 1
TEXT = 2

SELECT_BG_COLOR = '#D0B0A0'
SELECT_FG_COLOR = '#908070'
LINE_COLOR='#A0A0A0'
NODE_BG_COLOR = '#FFFFFF'
BG_COLOR = '#D4D4D4'

ICON_SIZE_DICT = {'small':16,'medium':22,'large':32}

class Tree(Frame):

  def __init__(self, parent, iconSize='medium', multiSelect=False,
               font='Helvetica 8', doubleCallback=None, *args, **kw):
  
    Frame.__init__(self, parent, *args, **kw)
    
    self.multiSelect = multiSelect
    self.icons = {}
    self.openDict = {}
    self.scrollbarWidth = 15
    self.canvasList = []
    self.canvasDict = {}
    self.canvasLines = []
    self._wait = False
    self.vOffset = 0
    self.nodes = [] # Just a list of nodes
    self.nodeDict = {} # To fetch nodes via objects
    self.font = font
    self.doubleCallback = doubleCallback

    #bg = self.cget('background')

    iconItems = [{ 'kind': 'command', 'label':'Small','command': self.smallIcons },
                 { 'kind': 'command', 'label':'Medium','command': self.medIcons },
                 { 'kind': 'command', 'label':'Large','command': self.largeIcons }]
                 
    self.menu = Menu(self, tearoff=False)
    menu_items = [{'kind': 'cascade',
                   'label': 'Icon size',
                   'submenu': iconItems },
                 ]

    self.menu.setMenuItems(menu_items)

    self.rowHeight = None
    self.visibleRows = 80
  
    self.setIconSize(iconSize, redraw=False)
    
    self.canvas = Canvas(self, relief='flat', borderwidth=0,
                         background=BG_COLOR, xscrollcommand=self._setHScrollbar)

    self.canvas.pack()

    self.xScrollbar  = Scrollbar(self, orient='horizontal', width=self.scrollbarWidth,
                                    borderwidth=1, callback=self._moveHScrollbar,
                                    background=BG_COLOR)

    self.yScrollbar  = Scrollbar(self, orient='vertical', width=self.scrollbarWidth,
                                    borderwidth=1, callback=self._setVScrollbar,
                                    background=BG_COLOR)
    
   
    self.canvas.bind('<Button-1>', self._mouseClick)
    self.canvas.bind('<Double-1>', self._mouseDoubleClick)

    if useWheelMouse():
      self.canvas.bind('<MouseWheel>', self._windowsOsScroll)
    else:
      self.canvas.bind('<Button-4>', self._mouseUp)
      self.canvas.bind('<Button-5>', self._mouseDown)

    self.canvas.bind('<p>', self._printCanvas)
    self.canvas.bind('<KeyPress-Prior>', self._pageUp)
    self.canvas.bind('<KeyPress-Next>',  self._pageDown)
    self.canvas.bind('<KeyPress-Up>',    self._keyUp)
    self.canvas.bind('<KeyPress-Down>',  self._keyDown)
    self.canvas.bind('<KeyPress-Home>',  self._keyHome)
    self.canvas.bind('<KeyPress-End>',   self._keyEnd)
    self.canvas.bind('<Enter>',  self._enter)

    # self.canvas.bind('<ButtonPress-3>', self._popupMenu)
    self.canvas.bind(buttonPress(3), self._popupMenu)

    self.bind('<Configure>',self._changeSizeAfter)
    
  def _popupMenu(self, event):
  
    self.menu.popupMenu(event)

  def smallIcons(self):
  
    self.setIconSize('small')
  
  def medIcons(self):
  
    self.setIconSize('medium')
  
  def largeIcons(self):
  
    self.setIconSize('large')

  def setIconSize(self, iconSize, redraw=True):
  
    size = ICON_SIZE_DICT.get(iconSize, 22)
    if size != self.rowHeight:
      self.rowHeight = size + 1
      self.font = 'Helvetica %d' % (size/2)

      ccpnDir = getTopDirectory()
      iconDir = '%dx%d' % (size,size)
      imageDir = path.join(ccpnDir,'python','memops','gui','graphics',iconDir)
 
      files = [f for f in listdir(imageDir) if f.endswith('.gif')]
 
      for file in files:
        imageFile = path.join(imageDir,file)
        self.icons[file[:-4]] = PhotoImage(file=imageFile)

      if redraw:
        self._changeSize()
      
    return size

  def _setHScrollbar(self, start, end):

    self.xScrollbar.set(float(start), float(end))

  def _moveHScrollbar(self, start, end):
  
    self.canvas.xview('moveto',start)
  
  def _enter(self, event):
  
    self.menu.popdownMenu()
    self.canvas.focus_force()

  def _printCanvas(self, *event):
  
    from memops.gui.FileSelect import FileType
    from memops.gui.FileSelectPopup import FileSelectPopup

    fileTypes = [  FileType('PostScript', ['*.ps']), FileType('All', ['*'])]
    fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
               title = 'Print canvas to file', dismiss_text = 'Cancel',
               selected_file_must_exist = False)

    fileName = fileSelectPopup.getFile()
        
    self.canvas.postscript(colormode='color',file=fileName)

  def _windowsOsScroll(self, event):
  
    delta = event.delta
    if delta > 0:
      self._mouseUp(event)
    elif delta < 0:
      self._mouseDown(event)

  def _keyHome(self, event):

    self.vOffset = 0
    self._drawAfter()

  def _keyEnd(self, event):

    numNodes = len(self.nodes)
    visNodes = self.visibleRows
    self.vOffset = max(0, numNodes-visNodes)
    self._drawAfter()

  def _scroll(self, delta=1):

    height = int(self.winfo_height())
    numNodes = len(self.nodes)
    visNodes = self.visibleRows
  
    if numNodes <= visNodes:
      return

    vOffset = self.vOffset+delta
    vOffset = max(0, vOffset)
    vOffset = min(numNodes-visNodes, vOffset)

    if vOffset != self.vOffset:
      self.vOffset = vOffset
      self._drawAfter()

  def _pageUp(self, event):
  
    self._scroll(-1*self.visibleRows)

  def _pageDown(self, event):
  
    self._scroll(self.visibleRows)

  def _mouseUp(self, event):

    self._scroll(-8)

  def _mouseDown(self, event):
  
    self._scroll(8)
 
  def _keyUp(self, event):

    self._scroll(-1)

  def _keyDown(self, event):
  
    self._scroll(1)

  def _setVScrollbar(self, start, end):

    self.cancelEdits()

    height = int(self.winfo_height())

    gaps = start+1-end

    nRows = float(len(self.nodes))
        
    if gaps == 0 or nRows == 0.0:
      self.yScrollbar.set(0,1)
      self.vOffset = 0
    else:

      # Have to do this so that the bar keeps up with the mouse
      overlap = int(nRows - self.visibleRows)
      size = min(1.0, self.visibleRows/nRows )
      self.yScrollbar.set(min(gaps, max(0, start)), max(size, min(1.0, start + size)))
      self.vOffset = min(max(0, int('%1.1d' % (start/gaps * (overlap+1)) ) ), overlap)
     
      #vOffset = max(0, start*nRows)
      #vOffset = min(nRows-self.visibleRows, vOffset)
      #self.vOffset = int(vOffset)

    self._drawAfter()

  def _updateScrollBars(self):
   
    allRows = len(self.nodes)
    visRows = self.visibleRows
        
    if allRows < visRows:
      self.vOffset = 0
      self.yScrollbar.set(0,1)
      self.yScrollbar.place_forget()
    
    else:
      size = self.scrollbarWidth
      width  = int(self.winfo_width())-1
      height = int(self.winfo_height())-1
      start = float(self.vOffset)/allRows
      end = float(self.vOffset+visRows)/allRows

      self.yScrollbar.set(start, end)
      self.yScrollbar.place(x=width-size+1,y=1,width=size,height=height)
              
    self.update_idletasks()    

  def _mouseClick(self, event):

    self.cancelEdits()
    
    x = self.canvas.canvasx(event.x)
    y = self.canvas.canvasy(event.y)
    obj = self.canvas.find('closest',x,y)[0]
    node, typ = self.canvasDict.get(obj, (None, None))
    
    if node:
      if typ is TOGGLE:
        node.toggle()
        
      else:
        if self.multiSelect and (event.state & 4 or event.state & 1):
          if event.state & 4:
            node.isSelected = not node.isSelected
            
          elif event.state & 1:
            indices = [self.nodes.index(node),]
            for i, node2 in enumerate(self.nodes):
              if node2.isSelected:
                indices.append(i)
        
            for node2 in self.nodes[min(indices):max(indices)+1]:
              node2.isSelected = True
        
        else:
          for node2 in self.nodes:
            if node2.isSelected:
              node2.isSelected = False
        
          node.isSelected = True
      
      self._drawAfter() 
      
  def _mouseDoubleClick(self, event):

    self.cancelEdits()
    
    x = self.canvas.canvasx(event.x)
    y = self.canvas.canvasy(event.y)
    obj = self.canvas.find('closest',x,y)[0]
    node, typ = self.canvasDict.get(obj, (None, None))
 
    if node:
      if typ is TOGGLE:
        node.toggle()
 
      elif self.doubleCallback:
        self.doubleCallback(node)
        
    
    self._drawAfter() 

  def _drawAfter(self):
  
    if self._wait:
      return
    
    self._wait = True
    self.after_idle(self._draw)  

  def _draw(self):
  
    self._updateScrollBars()
  
    delta = self.rowHeight
    canvas = self.canvas    
    cImage = canvas.create_image
    cList = self.canvasList
    cDict = self.canvasDict
    cLines = self.canvasLines
    cConfig = canvas.itemconfigure
    cDel = canvas.delete
    cCoords = canvas.coords
    cText = canvas.create_text
    cRect = canvas.create_rectangle
    cBbox = canvas.bbox
    cLine = canvas.create_line
    cLift = canvas.lift
    cLower = canvas.lower
    cPoly = canvas.create_polygon
    nodes = self.nodes
    nodeDict = self.nodeDict
    visibleRows = self.visibleRows
    icons = self.icons
    font = self.font
    vOffset = self.vOffset
    bg = BG_COLOR # self.cget('bg')

    roots = [n for n in nodes if not n.parent]
    
    if not roots:
      raise Exception('No rooted nodes')

    for item in cLines:
      cDel(item)

    stack = [(node, 0) for node in roots]
    nodes = []
    while stack:
      node, level = stack.pop(0)
      nodes.append((node, level))
      
      if node.isOpen:
        level2 = level+1
        stack = [(c, level2) for c in node.children] + stack
    
    yDict = {}
    pLevel = 0    
    pNode = None
    pad = 2 + (delta/2) 
    bSpace = 4
    nodes = nodes[vOffset:]
    for row, data in enumerate(nodes):
    
      node, level = data

      y = pad + row*delta
      x = pad + level*delta
      x1 = x + delta
      x2 = x1 + delta - bSpace
      
      yDict[node] = y, x1
      
      yP, xP = yDict.get(node.parent, (0, x))
        
      if yP is not None:
        line = cLine(xP, yP, xP, y, fill=LINE_COLOR, width=1)
        cLower(line)
        cLines.append(line)
              
      if row > visibleRows:
        break 
                
      iconImage = icons[node.icon]
      
      if node.isSelected:
        bgColor = SELECT_BG_COLOR
        fgColor = SELECT_FG_COLOR
      else:
        bgColor = bg
        fgColor = bg
      
      if row < len(cList):
        line, rect, tog1, tog2, box, icon, text = cList[row]
        cConfig(rect, fill=bgColor, outline=fgColor)
        cConfig(icon, image=iconImage)
        cConfig(text, text=node.label or '', font=font)
        cCoords(icon, x1, y)
        cCoords(text, x2, y)
        cCoords(line, x, y, x1, y)

      else:
        line = cLine(x, y, x1, y, fill=LINE_COLOR)
        rect = cRect(x,y,x2,y, fill=bgColor, outline=fgColor)
        box  = cRect(-2,-2,-2,-2, fill=NODE_BG_COLOR)
        tog1 = cLine(-2,-2,-2,-2)
        tog2 = cLine(-2,-2,-2,-2)
        icon = cImage(x1, y, image=iconImage)
        text = cText(x2, y, text=node.label or '', justify='left', anchor='w', font=font)
        cList.append( (line, rect, tog1, tog2, box, icon, text) )
      
      if node.callback:
        cCoords(box, x-4, y-4, x+4, y+4)
        
        if node.isOpen:
          cCoords(tog1, -2,-2,-2,-2)
          cCoords(tog2, x-2, y, x+3, y)
          
        else:
          cCoords(tog1, x, y-2, x, y+3)
          cCoords(tog2, x-2, y, x+3, y)
      
      else:
        cCoords(tog1, -2,-2,-2,-2)
        cCoords(tog2, -2,-2,-2,-2)
        cCoords(box, -2,-2,-2,-2)

      bbox = cBbox(text)
      if bbox:
        cCoords(rect, bbox[0]-2, bbox[1], bbox[2]+2, bbox[3])
      
 
      cLift(box)        
      cLift(tog1)
      cLift(tog2)
      
      cDict[tog1] = (node, TOGGLE)
      cDict[tog2] = (node, TOGGLE)
      cDict[box]  = (node, TOGGLE)
      cDict[text] = (node, TEXT)
      cDict[icon] = (node, ICON)
    
    nRows = min(visibleRows, len(nodes))    
    if len(cList) > nRows:
      for row in range(nRows, len(cList)):
        for item in cList[row]:
          if cDict.get(item):
            del cDict[item]
          cDel(item)
      
      self.canvasList = cList[:nRows]    

    self._wait = False

  def setFont(self, font):
    
    self.font = font
    self._drawAfter()

  def _changeSizeAfter(self, event):
  
    self.after_idle(lambda: self._changeSize(event))


  def _changeSize(self, event=None):
    
    self.cancelEdits()
    
    if event:
      width  = event.width
      height = event.height
    else:
      width  = int(self.winfo_width())
      height = int(self.winfo_height())

    self.canvas.config(width=width, height=height)
    self.visibleRows = height/self.rowHeight
    
    bbox = self.canvas.bbox('all')
    
    if bbox: # None at startup: nothing drawn
      canvasWidth = bbox[2]-bbox[0]
      if width < canvasWidth:
        size = self.scrollbarWidth
        end = width/float(canvasWidth)
        self.canvas.config(scrollregion="0 0 %s %s" % (canvasWidth, height) )

        if self.visibleRows < len(self.nodes):
          self.xScrollbar.place(x=0,y=height-size,width=width-size,height=size)
        else:
          self.xScrollbar.place(x=0,y=height-size,width=width,height=size)
        
        self.xScrollbar.set(0, end)
      else:
        self.xScrollbar.place_forget()
    else:
      self.xScrollbar.place_forget()

    self._drawAfter()

  # Public

  def cancelEdits(self):
  
    self.menu.popdownMenu()

  def getSelected(self):
  
    return [node.object for node in self.nodes if node.isSelected]
  
  def select(self, object, expand=True):
  
    node = self.nodeDict.get(object)
    if node is not None:
      node.isSelected = True
    
      if expand:
        parent = node.parent
        
        while parent:
          parent.isOpen = True
          parent = parent.parent

      self._drawAfter()

  def selectAll(self, expand=True):
  
    for node in self.nodes:
      node.isSelected = True
      if expand:
        node.isOpen = expand
  
  def selectDeep(self, object, expand=True):
  
    node = self.nodeDict.get(object)
    if node is not None:
      stack = [node,]
      nodes = []
      
      while stack:
        node2 = stack.pop()
        stack.extend(node2.children)
        nodes.append(node2)
    
      for node2 in nodes:
        if expand:
          node2.isOpen = True
        node2.isSelected = True
      
      self._drawAfter()  
    
    
  def setSelected(self, objects, expand=True):
  
    select = self.select
    for object in objects:
      select(object, expand)

  def expand(self, object):
  
    node = self.nodeDict.get(object)
    if node is not None:
      node.expand()
      self._drawAfter()

  def expandAll(self):
  
    for node in self.nodes:
      node.expand()
    
    self._drawAfter()
  
  def expandDeep(self, object):
  
    node = self.nodeDict.get(object)
    if node is not None:
      stack = [node,]
      nodes = []
      
      while stack:
        node2 = stack.pop()
        stack.extend(node2.children)
        nodes.append(node2)
    
      for node2 in nodes:
        node2.iexpand()
      
      self._drawAfter()  
    
  def contract(self, object):
  
    node = self.nodeDict.get(object)
    if node is not None:
      node.collapse()
      self._drawAfter()

  def contractAll(self):
  
    for node in self.nodes:
      node.collapse()
    
    self._drawAfter()

  def contractDeep(self, object):
  
    node = self.nodeDict.get(object)
    if node is not None:
      stack = [node,]
      nodes = []
      
      while stack:
        node2 = stack.pop()
        stack.extend(node2.children)
        nodes.append(node2)
    
      for node2 in nodes:
        node2.collapse()
      
      self._drawAfter()  

  def update(self, parents, objects, labels=None, icons=None,
             callbacks=None, editWidgets=None):
  
    n = len(parents)
    
    for test in (objects, callbacks, labels, icons, editWidgets):
      if test is not None:
        assert len(test) == n 
  
    self.nodes = nodes = []
    self.nodeDict = nodeDict = {}
    self.openDict = {}
    
    if not labels:
      labels = [None] * n
    
    if not icons:
      icons = [None] * n  

    if not callbacks:
      callbacks = [None] * n  

    if not editWidgets:
      editWidgets = [[] for x in xrange(n)]  

    for i, parent in enumerate(parents):
      object = objects[i]
      
      node = Node(self, parent, object, labels[i], icons[i],
                  callbacks[i], False, False, editWidgets[i])


  def add(self, parent, object, label=None, icon='media-playback-stop',
          callback=None, isOpen=False, isSelected=False, editWidgets=None):
    
    nodeDict = self.nodeDict
    if nodeDict.get(object):
      nodeDict[object].delete()
      
      #raise Exception('Node already present for object %s' % object)
  
    node = Node(self, parent, object, label, icon, callback, 
                isOpen, isSelected, editWidgets)

    
  def remove(self, object):
  
    node = self.nodeDict.get(object)
    if node is not None:
      if self.openDict.get(object):
        del self.openDict[object]
        
      node.delete()

      self._drawAfter()

class Node:

  def __init__(self, tree, parent, object, label=None, icon='list-add', callback=None,
               isOpen=False, isSelected=False, editWidgets=None):

    self.tree = tree

    self.callback = callback
    self.children = []
    self.editWidgets = editWidgets or []
    self.icon = icon
    self.isSelected = isSelected
    self.label = label
    self.object = object
    self.parent = parent

    if parent is not None:
      parent.children.append(self)
    
    self.isOpen = tree.openDict.get(object, isOpen)
    
    if parent in tree.nodes:
      index = tree.nodes.index(parent)
      tree.nodes.insert(index+len(parent.children), self)
    else:  
      tree.nodes.append(self)
    
    tree.nodeDict[object] = self  

  def expand(self):
  
    self.isOpen = True
    if self.callback:
      tree = self.tree
      tree.openDict[self.object] = True
      self.callback(self)
      
      for child in self.children:
        if tree.openDict.get(child.object):
          child.expand()
  
  def collapse(self):
           
    self.isOpen = False 
    if self.callback:
      tree.openDict[self.object] = False
      
      for child in self.children[:]:
        child.delete()
    
  def toggle(self):
  
    if self.isOpen:
      self.collapse()
    else:
      self.expand()
    
  def delete(self):
    
    for node in self.children:
      node.delete()

    parent = self.parent  
    if parent is not None:
      parent.children.remove(self)

    tree = self.tree
    object = self.object

    #if tree.openDict.get(object):
    #  del tree.openDict[object]
    
    del tree.nodeDict[object]
    tree.nodes.remove(self)
    
    del self

class TreeEditWidget:

  def __init__(self):
  
    pass

if __name__ == '__main__':

  # File browse

  import Tkinter, mimetypes
  
  from os import R_OK, access
  
  root = Tkinter.Tk()
  root.grid_columnconfigure(0, weight=1)
  root.grid_rowconfigure(0, weight=1)
  
  
  icon_dict = {'none':'text-x-generic-template',
               'text':'text-x-generic',
               'text/plain':'text-x-generic',
               'text/html':'text-html',
               'text/script':'text-x-script',
               'text/x-python':'text-x-script',
               'text/x-perl':'text-x-script',
               'image':'image-x-generic',
               'application':'applications-other',
               'video':'video-x-generic',
               'message':'internet-mail',
               'audio':'audio-x-generic',
               'chemical':'applications-internet',
               'application/x-tar':'package-x-generic',
               'application/msword':'x-office-document',
               'application/postscript':'x-office-drawing',
               'application/pdf':'x-office-drawing',
               'application/xml':'text-x-script',
               'application/x-python-code':'application-x-executable',}
               
  def getData(file):
    
    dirr, text = path.split(file)
    
    readable = access(file, R_OK)
    
    isDir = False
    if path.isdir(file):
      
      if not readable:
        icon = 'emblem-unreadable'
      elif path.islink(file):
        isDir = True
        icon = 'folder-remote'
      else:
        isDir = True
        icon = 'folder'
      text += '/'  
        
    elif path.isfile(file):
      
      if not readable:
        icon = 'emblem-unreadable'
      #elif path.islink(file):
      #  icon = 'text-html'
      else:
        mimeType = mimetypes.guess_type(file)[0]
        
        if mimeType:
          if not icon_dict.get(mimeType):
            print mimeType
        
          icon = icon_dict.get(mimeType) \
                  or icon_dict.get(mimeType.split('/')[0]) \
                  or icon_dict['none']
        else:
          icon = icon_dict['none']
    
    elif not readable:
      icon = 'emblem-unreadable'
      
    else:
      icon = 'emblem-system'
    
    return text, icon, isDir

  def doubleClick(node):
    data = (node.isOpen, node.label, node.object, len(node.children))
    print 'Callback on node -  open:%s label:%s, object:%s, children:%d' % data

  
  def openDir(node):
    
    tree = node.tree
    filePath = node.object
    files = listdir(filePath)
        
    for file in files:
       
       file = path.join(filePath, file)
       text, icon, isDir = getData(file)
       parent, null = path.split(file)
  
       if isDir and access(file, R_OK):
         callback = openDir
       else:
         callback = None
           
       tree.add(node, file, text, icon, callback)
  
  tree = Tree(root, multiSelect=True, grid=(0,0), doubleCallback=doubleClick)

  dirr = '/home'

  files = listdir(dirr)   
  parents = [None] * len(files)
  icons = []
  objects = []
  texts = []
  callbacks = []
  for file in files:
    file = path.join(dirr, file)
    
    text, icon, isDir = getData(file)
    icons.append(icon)
    texts.append(text)    
    objects.append(file)
    
    if isDir:
      callbacks.append(openDir)
    else:
      callbacks.append(None)
    
  tree.update(parents, objects, texts, icons, callbacks)

  root.mainloop()

  


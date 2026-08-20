
"""
======================COPYRIGHT/LICENSE START==========================

PulldownMenu.py: <write function here>

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
# old code so probably of dubious coding style at best

# TBD check if list is really long and then go sideways
# 
#  root = Tkinter.Tk()
#  root.withdraw()
#
#  popup = Popup(parent)
#  (w,h) = popup.maxsize()
#  popup.geometry('%dx%d+0+0' % (w,h-50))

import types

import Tkinter


from memops.gui.Canvas import Canvas
from memops.gui.Frame import Frame
from memops.gui.ScrolledFrame import ScrolledFrame
from memops.gui.Label import Label
from memops.gui.Menu import Menu
from memops.gui.Color import scaleColor

No_entries_label = '<None>'

class PulldownMenu(Frame):

  # if using indentation then should use list as stack only
  # in other words only delete or insert at end

  # entries entries must have attribute text_attr or if this is
  # None then assume that entries entries themselves provide name

  def __init__(self, parent, callback = None, entries = None, text_attr = None, colors=None,
               outline='#602000', fill='#B05848', label_color = '#501000',
               selected_index=-1, indent='', extra='', sticky='w',
               do_initial_callback=True, force_callback=False,
               docKey=None, tipText=None, *args, **kw):

    if (entries is None):
      entries = []

    Frame.__init__(self, parent, sticky=sticky, docKey=docKey, tipText=tipText,
                   createToolTip=True, *args, **kw)

    self.text_attr = text_attr

    #indentation and extra string for popup information
    self.indent = indent
    self.extra = extra

    self.callback = callback
    self.do_initial_callback = do_initial_callback
    self.force_callback = force_callback

    frame = Frame(self)
    #frame = parent
    self.frame   = frame
    self.bg      = frame.cget('bg')
    self.outline = outline
    self.fill    = fill
    self.label_color = label_color

    self.menu = Menu(parent, tearoff=0, borderwidth=1, activeborderwidth=1)
    self.menu.images = []
    self.label = Label(self, foreground = label_color)

    self.first_pass = True
    #self.entries = entries
    self.entries = []
    self.colors = []
    self.setup(entries, selected_index, colors=colors)

    s = 10
    self.canvas = Canvas(self, width=s, height=s, background=self.bg)

    self.label.bind("<Button-1>", self.labelPopup)
    # below does not work for some reason
    #self.bind("<Button-1>", self.labelPopup)

    # below does not work any more since can have submenus and
    # no way to reach those with below since involves leaving menu
    #self.menu.bind("<Leave>", self.popdown)
    self.menu.bind("<Leave>", self.leave)
    self.menu.bind("<Button-1>", self.buttonPress)
    self.poppedUpSubmenu = False

    self.canvas.bind("<Button-1>", self.canvasPopup)
    self.canvas.bind("<Configure>", self.resizeCallback)

    self.label.grid(row=0, column=0, sticky=Tkinter.W)
    self.canvas.grid(row=0, column=1, sticky=Tkinter.E, padx=2)

  def leave(self, event):

    if (not self.poppedUpSubmenu):
      x = event.x
      y = event.y
      x1 = self.menu.winfo_width()
      y1 = self.menu.winfo_height()
      
      if (x < 0) or (y<0) or (x>=x1) or (y>=y1):
        self.popdown()

  def buttonPress(self, event):

    if (event.widget != self.menu):
      return

    for n in range(len(self.entries)-1, -1, -1):
      if event.y >= self.menu.yposition(n):
        break
    else:
      return

    item = self.entries[n]
    if (isinstance(item, {}.__class__) and item['kind'] == 'cascade'):
      self.poppedUpSubmenu = True
    else:
      self.poppedUpSubmenu = False

  def resizeCallback(self, *event):

    w = self.canvas.winfo_width() / 2
    h = self.canvas.winfo_height() / 2
    self.canvas.delete('all')
    self.canvas.create_polygon(1, 1, w, 2*h-1, 2*w-1, 1, fill=self.fill, outline=self.outline)

  def getEntry(self, entry_index, entries = None):

    if not entries:
      entries = self.entries

    if type(entry_index) == types.TupleType:
      e = entry_index[0]
      if len(entry_index) > 1:
        items = entries[e]['submenu']
        entry = self.getEntry(entry_index[1:], items)
      else:
        entry = entries[e]['label']
    elif (entry_index >= 0):
      entry = entries[entry_index]
    else:
      entry = None

    return entry

  def getEntryText(self, entry_index, entries = None):
    
    if not entries:
      entries = self.entries

    if type(entry_index) in (types.TupleType, types.ListType):
      e = entry_index[0]
      if len(entry_index) > 1:
        items = entries[e]['submenu']
        text = self.getEntryText(entry_index[1:], items)
      else:
        text = entries[e]['label']
    elif (entry_index != -1):
      #print 'getEntryText', entry_index, entries
      entry = entries[entry_index]
      if (self.text_attr):
        text = getattr(entry, self.text_attr)
      else:
        text = entry
    else:
      text = ''

    # TBD: for now text returned should not be unicode
    text = str(text)
    # below did not work for one person (Chinese characters??) but above did
    ###text = text.encode('utf-8')

    return text

  def setup(self, entries, selected_index, first_pass=True, colors=None):

    #print 'setup', entries, selected_index
    if first_pass is not None:
      self.first_pass = first_pass
      
    self.clearMenuItems()
    if (entries):
      e = list(entries)
    else:
      e = [No_entries_label]

    self.entries = e
    self.colors  = colors or [None] * len(e)
    
    if type(selected_index) in (types.TupleType, types.ListType):
      if (not selected_index or selected_index[0] >= len(entries)):
        selected_index = -1
    else:
      if (selected_index >= len(entries)):
        selected_index = -1

    if (selected_index == -1):
      selected_index = 0

    self.selected_index = -1 # changed in setSelectedIndex()
    self.setMenuItems()
    self.setSelectedIndex(selected_index)

  def substituteCascadeCallbacks(self, item, selected_index):

    if (item['kind'] == 'cascade'):
      n = 0
      for subitem in item['submenu']:
        self.substituteCascadeCallbacks(subitem, selected_index+[n])
        n = n + 1
    else:
      callback = item.get('command')
      if (callback):
        label = item['label']
        item['command'] = lambda :self.cascadeCallback(callback, label, tuple(selected_index))

  def cascadeCallback(self, callback, label, selected_index):

    self.selected_index = selected_index
    self.setLabel(selected_index)
    callback(selected_index, label)

  def setMenuItems(self):

    self.menu.images = []
    for n in range(len(self.entries)):
      columnbreak = 0
      if n and n % 20 == 0:
        columnbreak = 1
      item  = self.entries[n]
      color = None
      if n < len(self.colors):
        color = self.colors[n]
      
      if isinstance(item, {}.__class__):
        if (item['kind'] == 'cascade'):
          self.substituteCascadeCallbacks(item, [n])
      else:
        t = (self.indent * n) + self.extra + str(self.getEntryText(n))
        command = lambda n=n: self.setSelectedIndex(n)
        if color:
          image = self.makeColorTile(color)
          item = { 'kind': 'command', 'accelerator': t, 'command': command,
                   'image': image, 'columnbreak': columnbreak}
        else:
          item = { 'kind': 'command', 'label': t, 'command': command, 'columnbreak': columnbreak }
      self.menu.addMenuItem(item)

  def clearMenuItems(self):

    if self.entries:
      self.menu.delete(0, len(self.entries))

  # below assumes entry_index is integer not tuple
  def insert(self, entry_index, entry, make_selected = False, color=None):

    if (self.haveNoEntries()):
      if (entry_index != 0):
        return
    elif ((entry_index < 0) or (entry_index >= len(self.entries))):
      return

    self.clearMenuItems()
    self.checkForNone()
    self.entries.insert(entry_index, entry)
    self.colors.insert(entry_index, color)
    self.setMenuItems()

    if (make_selected or (self.selected_index == -1)):
      self.setSelectedIndex(entry_index, force_callback = True)
    else:
      self.setSelectedIndex(self.selected_index)

  # below assumes entry is integer not tuple
  def append(self, entry, make_selected = False, color=None):

    self.clearMenuItems()
    self.checkForNone()
    self.entries.append(entry)
    self.colors.append(color)
    self.setMenuItems()

    if (make_selected or (self.selected_index == -1)):
      self.setSelectedIndex(len(self.entries)-1, force_callback = True)
    else:
      self.setSelectedIndex(self.selected_index)

  def checkForNone(self):

    if (self.haveNoEntries()):
      self.entries = []
      self.colors  = []
      self.menu.delete(0)

  def haveNoEntries(self):

    if (len(self.entries) == 1 and self.entries[0] == No_entries_label):
      return True
    else:
      return False

  def baseIndex(self, ind):

    if type(ind) == types.TupleType:
      ind = ind[0]

    return ind

  # below assumes entry_index is integer not tuple
  def delete(self, entry_index, n = 1):

    if (self.haveNoEntries()):
      return

    if ((entry_index < 0) or (entry_index >= len(self.entries))):
      return

    self.clearMenuItems()

    n1 = entry_index
    n2 = n1 + n

    m = self.baseIndex(self.selected_index)

    del self.entries[n1:n2]
    del self.colors[n1:n2]

    if (not self.entries):
      self.entries = [No_entries_label]
      self.colors = [None]

    self.setMenuItems()

    if ((m >= n1) and (m < n2)):
      if (n1 == 0):
        selected_index = 0
      else:
        selected_index = n1-1
      self.setSelectedIndex(selected_index, force_callback = True)

  def replace(self, entries, selected_index = -1, colors=None):

    self.menu.delete(0, len(self.entries))
    self.setup(entries, selected_index, first_pass=False, colors=colors)

  def popdown(self, *event):

    self.menu.unpost()

  def setLabel(self, selected_index):

    text = self.getEntryText(selected_index)
    if (not text):
      text = self.extra # arbitrary

    self.label.set(text=text)

  def labelPopup(self, event):

    self.poppedUpSubmenu = False

    s = self.baseIndex(self.selected_index)
    x = event.x_root - event.x + 2
    y = event.y_root - event.y - \
        max(0, s) * (self.label.winfo_height() + 1)
    self.menu.post(x, y)

  def canvasPopup(self, event):

    s = self.baseIndex(self.selected_index)
    x = event.x_root - event.x -self.label.winfo_width() + 2
    y = event.y_root - event.y - \
        max(0, s) * (self.label.winfo_height() + 1) + \
        self.canvas.winfo_height() - self.label.winfo_height()
    self.menu.post(x, y)

  def getSelectedIndex(self):

    return self.selected_index

  def setSelectedIndex(self, selected_index, force_callback = False):

    if (force_callback or self.force_callback or \
        (selected_index != self.selected_index)):

      #print 'setSelectedIndex', self.entries, selected_index, self.selected_index, self.first_pass
      self.selected_index = selected_index
      self.setLabel(selected_index)

      if (self.callback):
        if (self.haveNoEntries()):
          if (not self.first_pass):
            self.callback(-1, None)
        elif (selected_index != -1):
          if (not self.first_pass or self.do_initial_callback):
            self.callback(selected_index, self.entries[selected_index])

    self.first_pass = False

  def get(self):

    return self.getSelected()

  def getSelected(self):

    ind = self.getSelectedIndex()

    if (ind == -1):
      return None
    else:
      return self.getEntry(ind)

  def set(self, selected):

    self.setSelected(selected)

  def findEntryIndex(self, entry, entries = None):

    #print 'findEntryIndex1', entry, entries
    if (entries):
      noEntries = True
    else:
      noEntries = False
      entries = self.entries
      try:
        ind = entries.index(entry)
        return ind
      except:
        pass

    #print 'findEntryIndex2', entry, entries
    for n in range(len(entries)):
      e = entries[n]
      #print 'findEntryIndex3', entry, e
      if e == entry:
        if noEntries:
          return [n]
        else:
          return n
      elif type(e) == types.DictType:
        if (len(entries) > 1):
          try:
            ind = self.findEntryIndex(entry, e['submenu'])
            return [n] + ind
          except:
            pass
        else:
          if (entry == e['label']):
            return [n]

    raise Exception('unknown entry "%s"' % entry)

  # selects first item found in entries which matches
  def setSelected(self, selected):

    #print 'setSelected1', selected, type(selected), self.entries
    try:
      selected_index = self.findEntryIndex(selected)
    except:
      #selected_index = -1
      return

    #print 'setSelected2', selected_index
    self.setSelectedIndex(selected_index)

  def makeColorTile(self, color):
  
    image = Tkinter.PhotoImage()
    self.menu.images.append(image)
    
    if type(color) == type([]):
      colors = [ scaleColor(self.menu, c, 1.0) for c in color ] 
    else:
      colors = [ scaleColor(self.menu, color, 1.0), ]

    cols = max(8, len(colors))
    
    for x in range(cols):
      i = x % len(colors)
      c = colors[i]
      
      for y in range(16):
        image.put('{%s %s}' % (c,c), to=(2*x,y))
    
    return image

if __name__ == '__main__':
 
  from memops.gui.Button import Button

  import sys

  n = 0
  popup_select = None

  def callback(n, text):
    print 'callback: n =', n, 'text =', text

  def callback2(inds, text):
    print 'callback2: inds =', inds, 'text =', text

  def add():
    global n, popup_select
    n = n + 1
    text = str(n+1)
    popup_select.append(text, make_selected=True)

  def getSelected():
    global popup_select
    selected = popup_select.getSelected()
    print 'getSelected', selected

  def getSelectedInd():
    global popup_select
    selected = popup_select.getSelectedIndex()
    print 'getSelectedInd', selected

  def remove():
    global popup_select
    popup_select.delete(0)

  root = Tkinter.Tk()

  subitems = []
  for label in ('dim1', 'dim2', 'dim3'):
    subitems.append({'kind': 'command', 'label': label, 'command': callback2})
  item = { 'kind': 'cascade', 'label': 'cascade', 'submenu': subitems }
  entries = ["aaaa", "bbbbb", "cccccc", "dddddddeeeeffffgggghhhh", "asdfawer", item]
  colors = [None, None, 'red',['red','green'], ['#000040','#000080','#0000A0','#0000FF','#0040FF','#0080FF','#00A0FF','#00C0FF','#00E0FF','#00FFFF','#40FFFF','#80FFFF','#A0FFFF','#D0FFFF','#FFFFFF'], None]
  selected_index = 2

  popup_select = PulldownMenu(root, callback=callback, entries=entries,
                              colors=colors, selected_index=selected_index,
                              tipText='Explanation')
  popup_select.grid(row=0, column=0, columnspan=3)

  button = Button(root, text='add', command=add)
  button.grid(row=1, column=0)

  button = Button(root, text='remove', command=remove)
  button.grid(row=1, column=1)

  button = Button(root, text='getSelected', command=getSelected)
  button.grid(row=1, column=2)

  button = Button(root, text='getSelectedInd', command=getSelectedInd)
  button.grid(row=1, column=3)

  button = Button(root, text='quit', command=sys.exit)
  button.grid(row=1, column=4)

  root.mainloop()

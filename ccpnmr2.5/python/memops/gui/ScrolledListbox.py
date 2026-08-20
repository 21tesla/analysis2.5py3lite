
"""
======================COPYRIGHT/LICENSE START==========================

ScrolledListbox.py: <write function here>

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
import Tkinter, tkSimpleDialog

from memops.gui.ButtonList import ButtonList

from memops.gui.Frame import Frame

class ScrolledListbox(Frame):

  def __init__(self, parent, initial_list = None,
               width = 60, height = 5,
               xscroll = True, yscroll = True, addDeleteButtons = False,
               selectmode = Tkinter.BROWSE, exportselection = 0,
               select_callback = None, double_callback = None,
               list_background = 'white',
               *args, **kw):

    if (initial_list is None):
      initial_list = []

    Frame.__init__(self, parent, *args, **kw)

    self.selectmode = selectmode

    self.grid_rowconfigure(0, weight=1)
    self.grid_columnconfigure(0, weight=1)
    listbox = Tkinter.Listbox(self, width=width, height=height,
                      background = list_background,
                      selectmode=selectmode, exportselection=exportselection)
    listbox.grid(row=0, column=0, sticky=Tkinter.NSEW)

    if (xscroll):
      xscrollbar = Tkinter.Scrollbar(self, orient=Tkinter.HORIZONTAL)
      xscrollbar.config(command=listbox.xview)
      listbox.config(xscrollcommand=xscrollbar.set)
      xscrollbar.grid(row=1, column=0, sticky=Tkinter.EW)

    if (yscroll):
      yscrollbar = Tkinter.Scrollbar(self, orient=Tkinter.VERTICAL)
      yscrollbar.config(command=listbox.yview)
      listbox.config(yscrollcommand=yscrollbar.set)
      yscrollbar.grid(row=0, column=1, sticky=Tkinter.NS)

    if addDeleteButtons:
      texts = ['Add item','Delete selected']
      commands = [self.addItem,self.deleteItems]
      buttons = ButtonList(self,texts = texts, commands = commands)
      buttons.grid(row =1, columnspan = 2, sticky = Tkinter.EW)

    # bind frame, not listbox, because listbox with focus has
    # activation which means in particular that get underlining
    self.bind('<Enter>', self.enterCallback)
    self.bind('<KeyPress>', self.keypressCallback)

    self.listbox = listbox

    self.setItems(initial_list)

    self.setSelectCallback(select_callback)
    self.setDoubleCallback(double_callback)

    self.size = self.listbox.size # otherwise get Frame size called

  def moveSelectedItemUp(self, *event):

    if (self.selectmode != Tkinter.SINGLE):
      return

    ind = self.getSelected()

    if (ind):

      s = self.get(ind)
      self.delete(ind)
      self.insert(ind-1, s)
      self.select_set(ind-1)
      self.see(ind-1)

  def moveSelectedItemDown(self, *event):

    if (self.selectmode != Tkinter.SINGLE):
      return

    n = self.size()
    ind = self.getSelected()

    if ((ind != None) and (ind < n-1)):

      s = self.get(ind)
      self.delete(ind)
      self.insert(ind+1, s)
      self.select_set(ind+1)
      self.see(ind+1)

  def addItem(self):
  
    newItem = tkSimpleDialog.askstring('Enter new value','New value:',parent = self)
    
    if newItem:
      newItem = newItem.strip()
      self.append(newItem)

  def deleteItems(self):
  
    selectedItems = self.getSelectedItems()
    if type(selectedItems) != type([]):
      selectedItems = [selectedItems]
    
    unSelectedItems = []

    for item in self.listbox.get(0, Tkinter.END):
      if item not in selectedItems:
        unSelectedItems.append(item)
    
    self.setItems(unSelectedItems)
    
  def getAllItems(self):
  
    return self.listbox.get(0, Tkinter.END)

  def getSelected(self):

    inds = [ int(sel) for sel in self.listbox.curselection() ]

    if (self.selectmode == Tkinter.SINGLE):
      if (inds):
        return inds[0]
      else:
        return None
    else:
      return inds

  def getSelectedItems(self):

    items = [ self.get(sel) for sel in self.listbox.curselection() ]

    if (self.selectmode == Tkinter.SINGLE):
      if (items):
        return items[0]
      else:
        return None
    else:
      return items

  def setSelectedItems(self, items):

    entries = list(self.getItems())

    self.listbox.select_clear(0, Tkinter.END)

    if (self.selectmode == Tkinter.SINGLE):
      if (items):
        items = [items]

    for item in items:
      try:
        ind = entries.index(item)
        self.select_set(ind)
      except:
        pass
      
  def getItems(self):

    return self.get(0, Tkinter.END)

  def setItems(self, entries):

    self.clear()
    for e in entries:
      self.append(e)

  def append(self, x):

    self.listbox.insert(Tkinter.END, x)

  def clear(self):

    self.listbox.delete(0, Tkinter.END)

  def enterCallback(self, event):

    self.focus()

  def keypressCallback(self, event):

    #print dir(event), str(event.keysym), str(event.keycode)

    if (event.keysym in [ 'Up', 'Down' ]):

      selected = self.curselection()

      if (selected):

        ind = int(selected[0]) # tough if more than one selected

        if (event.keysym == 'Up'):

          if (ind != 0):
            ind = ind - 1

        else: # (event.keysym == 'Down')

          size = self.size()
          if (ind != (size-1)):
            ind = ind + 1

        self.select_clear(0, Tkinter.END)
        self.select_set(ind)
        self.see(ind)
        
        if (self.select_callback):
          self.select_callback(event)

  def setSelectCallback(self, select_callback):

    self.select_callback = select_callback

    if (select_callback):
      self.listbox.bind('<ButtonRelease-1>', select_callback)

  def setDoubleCallback(self, double_callback):

    self.double_callback = double_callback

    if (double_callback):
      self.listbox.bind('<Double-1>', double_callback)

  def __getattr__(self, name):

    # dispatch everything not defined by ScrolledListbox to listbox widget

    try:
      return getattr(self.__dict__['listbox'], name)
    except:
      raise AttributeError, "%s instance has no attribute '%s'" % (self.__class__.__name__, name)

if (__name__ == '__main__'):

  def myCallback(event):

    print 'myCallback:', listbox.getSelected(), listbox.getSelectedItems()

  def myCallback2(event):

    print 'myCallback2:', listbox.getSelected()
    #listbox.moveSelectedItemUp()
    listbox.moveSelectedItemDown()

  #
  # Use for list selection
  #

  root = Tkinter.Tk()

  listbox = ScrolledListbox(root, width=60, height=10,
              xscroll=False, selectmode=Tkinter.SINGLE,
              initial_list = ('init1', 'init2'),
              select_callback=myCallback, double_callback=myCallback2)

  for n in [ 'abc', '123', 'edf', 456 ]:
    listbox.append(str(n))

  #listbox.grid(sticky=Tkinter.NSEW)
  listbox.pack(expand=1, fill=Tkinter.BOTH)

  root.mainloop()
  
  
  #
  # Use for list editing
  #

  root = Tkinter.Tk()

  listbox = ScrolledListbox(root, width=60, height=10,
              xscroll=False, selectmode=Tkinter.SINGLE,
              initial_list = ('init1', 'init2'),
              select_callback=myCallback, double_callback=myCallback2,addDeleteButtons = True)

  for n in [ 'abc', '123', 'edf', 456 ]:
    listbox.append(str(n))

  #listbox.grid(sticky=Tkinter.NSEW)
  listbox.pack(expand=1, fill=Tkinter.BOTH)

  root.mainloop()
  
  listbox.getAllItems()

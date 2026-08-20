
"""
======================COPYRIGHT/LICENSE START==========================

RadioButtons.py: <write function here>

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

from memops.gui.Frame import Frame
from memops.gui.RadioButton import RadioButton

# option_dict keys and defaults:
#   borderwidth         1               # borderwidth for frame
#   padding             1               # padding for frame
#   selectcolor         None            # selection color for button (the diamond)
#   select_relief       Tkinter.SUNKEN  # relief for frame when button selected
#   unselect_relief     Tkinter.RAISED  # relief for frame when button unselected
#   select_fg           None            # fg for button when selected
#   unselect_fg         None            # fg for button when unselected
#   select_bg           None            # bg for button when selected
#   unselect_bg         None            # bg for button when unselected

class RadioButtons(Frame):

  def __init__(self, parent, entries = None, select_callback = None,
               direction = Tkinter.HORIZONTAL, option_dict = None, min_text_width=0,
               selected_index = 0, docKeys=None, tipTexts=None, *args, **kw):

    if option_dict is None:
      option_dict = {'selectcolor':'#D0B0A0'}

    if entries is None:
      entries = []

    if tipTexts and (len(tipTexts) == 1):
      # for length 1 tooltips just put on frame, rather than individual checkButtons
      tipText=tipTexts[0]
      if docKeys:
        docKey=docKeys[0]
      else:
        docKey = None
      tipTexts = None
      docKeys = None
    else:
      tipText = None
      docKey = None

    Frame.__init__(self, parent, docKey=docKey, tipText=tipText, createToolTip=True, *args, **kw)

    self.font = kw.get('font')
    self.var = Tkinter.IntVar()
    self.var.set(selected_index)
    self.direction = direction
    self.select_callback = select_callback
    self.min_text_width = min_text_width

    self.borderwidth = option_dict.get('borderwidth', 1)
    self.padding = option_dict.get('padding', 1)
    self.selectcolor = option_dict.get('selectcolor')
    self.select_relief = option_dict.get('select_relief', Tkinter.FLAT)
    self.unselect_relief = option_dict.get('unselect_relief', Tkinter.FLAT)

    self.select_dict = {}
    self.unselect_dict = {}
    for key in option_dict.keys():
      if (key[-7:] == '_relief'):
        continue
      elif (key[:7] == 'select_'):
        self.select_dict[key[7:]] = option_dict[key]
      elif (key[:9] == 'unselect_'):
        self.unselect_dict[key[9:]] = option_dict[key]
    
    self.entries = None
    self.setEntries(entries, docKeys, tipTexts)

  def setEntries(self, entries, docKeys=None, tipTexts=None):

    if entries == self.entries:
      return

    nentries = len(entries)
    if not docKeys:
      docKeys = [None] * nentries

    if not tipTexts:
      tipTexts = [None] * nentries

    if self.entries:
      for b in self.buttons:
        b.destroy()

    self.entries = entries
    self.buttons = nentries * [None]
    self.frames = nentries * [None]

    for n, entry in enumerate(entries):

      frame = Frame(self, borderwidth=self.borderwidth)

      entryText = entry
      m =  self.min_text_width - len(entryText)
      if m > 0:
        entryText = entryText + m * ' '
      
      bgColor = self.cget('bg')
      callback = lambda entry=entry, callback=self.selectCallback: callback(entry)
      button = RadioButton(frame, bg=bgColor, text=entryText, callback=callback,
                                   highlightbackground=bgColor, highlightcolor=bgColor,
                                   activebackground='#D0B0A0', activeforeground='white',
                                   variable=self.var, font=self.font, value=n,
                                   docKey=docKeys[n], tipText=tipTexts[n])

      if self.selectcolor:
        button.config(selectcolor=self.selectcolor)

      if self.direction == Tkinter.HORIZONTAL:
        frame.grid(row = 0, pady=0, column = n, sticky=Tkinter.NSEW, padx=self.padding)
        self.grid_columnconfigure(n, weight=1)
      else:
        frame.grid(row = n, padx=0, sticky=Tkinter.NSEW, pady=self.padding)
        self.grid_rowconfigure(n, weight=1)
 
      frame.grid_rowconfigure(0, weight=1)
      frame.grid_columnconfigure(0, weight=1)
      button.grid(sticky='nw')

      self.buttons[n] = button
      self.frames[n] = frame

    self.setOptions()

  def selectCallback(self, entry):

    self.setOptions()

    if self.select_callback:
      self.select_callback(entry)

  def get(self):

    n = self.getIndex()

    return self.entries[n]

  def getIndex(self):

    return self.var.get()

  def set(self, entry):

    n = self.entries.index(entry)
    self.setIndex(n)

  def setIndex(self, n):

    if n != self.getIndex():
      self.var.set(n)
      self.selectCallback(self.entries[n])

  def setOptions(self):

    n = self.getIndex()
    for m in range(len(self.frames)):
      if (m == n):
        relief = self.select_relief
        dict = self.select_dict
      else:
        relief = self.unselect_relief
        dict = self.unselect_dict
      self.frames[m].config(relief=relief)
      fg = dict.get('fg')
      if (fg):
        self.buttons[m].config(fg=fg)
      bg = dict.get('bg')
      if (bg):
        self.buttons[m].config(bg=bg)

if __name__ == '__main__':

  from memops.gui.Button import Button

  def get_me():
    print 'get_me:', r.get(), r.getIndex()

  def set_me():
    c = r.getIndex()
    c = (c + 1) % 3
    print 'set_me:', c
    r.setIndex(c)

  def my_callback(text):
    print 'my_callback:', text

  root = Tkinter.Tk()

  r = RadioButtons(root, ['one', 'two', 'three'], select_callback=my_callback,
                   tipTexts=['Tip A','Tip B','Tip C'])
  r.grid()
  b = Button(root, text='get me', command=get_me)
  b.grid()
  b = Button(root, text='set me', command=set_me)
  b.grid()

  root.mainloop()

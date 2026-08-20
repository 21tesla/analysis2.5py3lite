
"""
======================COPYRIGHT/LICENSE START==========================

LabeledEntries.py: <write function here>

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
import types

import Tkinter

from memops.universal.Util import checkListIsSet

from memops.gui.Entry import Entry
from memops.gui.Frame import Frame
from memops.gui.Label import Label

# this creates a table with left column a bunch of labels
# and right column a bunch of entries
# the separator is appended to the label, and should not be changed
# the entry texts can be initialised, or after initialisation
# set or get, via a list or via a dictionary keyed on the label
# the entry texts can be frozen, so only settable at constructor time
# in which case the entry becomes disabled
# the entry text can also be hidden using show, keyed on label

class LabeledEntries(Frame):

  def __init__(self, parent, labels = None, entries = None, separator = ': ',
               label_width = 20, entry_width = 60, label_anchor = Tkinter.E,
               show = None, frozen = None, returnCallback = None,
               put_in_row = False, *args, **kw):

    self.separator = separator
    self.label_width = label_width
    self.entry_width = entry_width
    self.label_anchor = label_anchor
    self.put_in_row = put_in_row
    self.returnCallback = returnCallback

    if (entries is None):
      entries = len(labels) * ['']

    nentries = len(entries)

    if (labels is None):
      labels = nentries * ['']

    assert len(labels) == nentries, '%d vs %d' % (len(labels), nentries)

    if (show is None):
      show = nentries * ['']

    if (frozen is None):
      frozen = nentries * [False]

    Frame.__init__(self, parent, *args, **kw)

    self.labels = []
    self.entries = []

    for n in range(nentries):
      self.addLabelEntry(n, labels[n], entries[n], show[n], frozen[n], doGrid=False)

    self.gridAll()

  def addLabelEntry(self, n, label = '', entry = '', show = '', frozen = False,
                    doGrid = True):

    text = label + self.separator
    l = Tkinter.Label(self, text=text, width=self.label_width, anchor=self.label_anchor)
    self.labels[n:n] = [l]

    e = Entry(self, text=entry, show=show, width=self.entry_width,
              returnCallback=lambda n=n: self.returnCallback(n))
    if (frozen):
      e.config(bg='gray', state=Tkinter.DISABLED)
      e.bind('<Button>', lambda event, e=e: e.focus())
    else:
      e.config(bg='white')
    self.entries[n:n] = [e]

    if (doGrid):
      self.gridAll()

  def gridAll(self):

    nentries = len(self.entries)

    row = 0
    col = 0
    for n in range(nentries):

      if (not self.put_in_row):
        col = 0

      l = self.labels[n]
      if (l):
        l.grid(row=row, column=col, sticky=Tkinter.EW)
        col = col + 1
      else:
        l.grid_forget()

      self.grid_columnconfigure(col, weight=1)

      e = self.entries[n]
      e.grid(row=row, column=col, sticky=Tkinter.EW)
      self.entries.append(e)
      col = col + 1

      if (not self.put_in_row):
        row = row + 1

  def getLabel(self, n):

    text = self.labels[n].get()
    m = text.find(self.separator)
    if (m >= 0):
      text = text[:m]

    return text

  def getLabels(self):

    labels = []

    for n in range(len(self.labels)):
      labels.append(self.getLabel(n))

    return labels

  def setLabel(self, n, text):

    text = text + self.separator
    self.labels[n].set(text)

  def getEntry(self, n):

    e = self.entries[n]

    return e.get()

  def getEntries(self):

    entries = []

    for n in range(len(self.entries)):
      entries.append(self.getEntry(n))

    return entries

  def setEntry(self, n, text):

    e = self.entries[n]

    # for frozen (DISABLED) entries this is just silently ignored
    e.set(text)

if __name__ == '__main__':
  
  root = Tkinter.Tk()

  """
  labels = ['user', 'password', 'age', 'height', 'weight']
  entries = ['wb104', '', '', '', '']
  frozen = [True, False, False, False, False]
  show = ['', '*', '', '', '']
  put_in_row = False
  label_width = 10
  entry_width = 20
  separator = ': '
  """

  labels = ['', '', '', '', '15N']
  entries = None
  frozen = None
  show = None
  put_in_row = True
  label_width = 5
  entry_width = 10
  separator = ''

  root.grid_columnconfigure(0, weight=1)
  w = LabeledEntries(root, labels, entries=entries, show=show, frozen=frozen,
                     label_width=label_width, entry_width=entry_width,
                     put_in_row=put_in_row, separator=separator)
  w.grid(sticky=Tkinter.EW)

  root.mainloop()

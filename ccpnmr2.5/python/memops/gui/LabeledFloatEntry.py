
"""
======================COPYRIGHT/LICENSE START==========================

LabeledFloatEntry.py: <write function here>

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


from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.Label import Label

class LabeledFloatEntry(Frame):

  def __init__(self, parent, label, entry = '', separator = ': ',
               label_width = 20, entry_width = 60,
               label_anchor=Tkinter.E, show = '', isArray = False,
               returnCallback = None, *args, **kw):

    apply(Frame.__init__, (self, parent) + args, kw)

    self.grid_columnconfigure(1, weight=1)

    self.separator = separator

    text = label + separator
    self.label = Label(self, text=text, width=label_width, anchor=label_anchor)
    self.label.grid(row=0, column=0, sticky=Tkinter.EW)

    self.entry = FloatEntry(self, show=show, width=entry_width, isArray=isArray,
                            returnCallback=returnCallback)
    self.entry.grid(row=0, column=1, sticky=Tkinter.EW)
    self.setEntry(entry)

  def getLabel(self):

    text = self.label.get()
    n = text.find(self.separator)
    if (n >= 0):
      text = text[:n]

    return text

  def setLabel(self, text = ''):

    text = text + self.separator
    self.label.set(text)

  def getEntry(self):

    return self.entry.get()

  def setEntry(self, value = ''):

    self.entry.set(value)

if __name__ == '__main__':
  
  root = Tkinter.Tk()

  entry1 = LabeledFloatEntry(root, 'x', '1.5')
  entry1.grid()

  entry2 = LabeledFloatEntry(root, 'list', (3.5, 6.6), isArray=True)
  entry2.grid()

  root.mainloop()

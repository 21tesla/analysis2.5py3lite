
"""
======================COPYRIGHT/LICENSE START==========================

LabeledScale.py: <write function here>

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
from memops.gui.Label import Label

class LabeledScale(Frame):

  def __init__(self, parent, values, initial_index = 0,
               orient = Tkinter.HORIZONTAL, scale_length = 200,
               label_width = 5, label_format = '%2.1f',
               set_callback = None,
               *args, **kw):

    assert len(values) > 0

    self.label_format = label_format
    self.values = tuple(values)
    self.set_callback = set_callback

    apply(Frame.__init__, (self, parent) + args, kw)

    if (orient == Tkinter.HORIZONTAL):
      self.grid_columnconfigure(1, weight=1)
      sticky = Tkinter.EW
      row = 0
      col = 1
    else:
      self.grid_rowconfigure(1, weight=1)
      sticky = Tkinter.NS
      row = 1
      col = 0

    self.label = Label(self, width=label_width)
    self.label.grid(row=0, column=0, sticky=sticky)

    self.scale = Tkinter.Scale(self, orient=orient, length=scale_length,
                               from_=0, to=len(values)-1, label='',
                               showvalue=0, command=self.setCallback)
    self.scale.grid(row=row, column=col, sticky=sticky)

    self.set(initial_index)

  def get(self):

    return self.scale.get()

  def getValue(self):

    return self.values[self.get()]

  def set(self, index):

    assert index >= 0 and index < len(self.values), 'index = %s, len(self.values) = %s' % (index, len(self.values))
 
    self.scale.set(index)
    #self.setLabel() # seems that automatically sets label via setCallback

  def setLabel(self):

    index = self.get()
    text = self.label_format % self.values[index]
    self.label.set(text)

  def setCallback(self, ignore):

    self.setLabel()
    if (self.set_callback):
      self.set_callback(self.get())

if __name__ == '__main__':
 
  import sys
  from memops.gui.Button import Button

  def getValue():

    print 'index = %s, value = %s' % (labeled_scale.get(), labeled_scale.getValue())

  def setCallback(index):

    print 'index = %s' % index

  root = Tkinter.Tk()

  labeled_scale = LabeledScale(root, values= (12.3, 34.5, 192, 12.3), set_callback=setCallback)
  labeled_scale.pack(side=Tkinter.TOP, expand=Tkinter.YES, fill=Tkinter.X)

  button = Button(root, text='get value', command=getValue)
  button.pack(side=Tkinter.TOP)

  button = Button(root, text='quit', command=sys.exit)
  button.pack(side=Tkinter.TOP)

  root.mainloop()

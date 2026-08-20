
"""
======================COPYRIGHT/LICENSE START==========================

ToggleArrow.py: <write function here>

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


from memops.gui.Canvas import Canvas
from memops.gui.Frame import Frame

class ToggleArrow(Frame):

  def __init__(self, parent, arrowSize = 11, outline='#5050b0', fill='#a0a0ff',
               isArrowClosed = True, callback = None, *args, **kw):
 
    self.arrowSize = arrowSize
    self.isArrowClosed = isArrowClosed
    self.callback = callback

    apply(Frame.__init__, (self, parent) + args, kw)

    bg = kw.get('bg')
    if (not bg):
      bg = kw.get('background')
    if (not bg):
      bg = parent.cget('bg')
    # below does not work for some reason, instead explicitly create rectangle below
    #self.canvas = Canvas(self, width=arrowSize, height=arrowSize, bg=bg)
    self.canvas = Canvas(self, width=arrowSize+3, height=arrowSize+1)
    self.canvas.grid(row=0, column=0)

    # need +/-1 otherwise get line along borders
    self.canvas.create_rectangle(-1, -1, arrowSize+3, arrowSize+1, fill=bg)

    # arbitrary coords
    self.arrow = self.canvas.create_polygon(0, 0, 1, 0, 1, 1, fill=fill, outline=outline)
    self.drawArrow()

    self.canvas.bind('<Button-1>', self.flipState)

  def drawArrow(self):

    s = self.arrowSize
    if (self.isArrowClosed):
      self.canvas.coords(self.arrow, 2, 0, s+2, s/2, 2, s)
    else:
      self.canvas.coords(self.arrow, 2, 0, s+2, 0, (s/2)+2, s)

  def flipState(self, *event):

    self.isArrowClosed = not self.isArrowClosed
    if (self.callback):
      self.callback(self.isArrowClosed)

    self.drawArrow()

if (__name__ == '__main__'):

  from Label import Label

  frame = None

  def callback(isClosed):

    if (isClosed):
      frame.grid_forget()
    else:
      frame.grid(row=1, column=0, columnspan=2, sticky=Tkinter.NSEW)

  root = Tkinter.Tk()
 
  root.grid_rowconfigure(1, weight=1)
  root.grid_columnconfigure(1, weight=1)

  arrow = ToggleArrow(root, callback=callback, arrowSize=20)
  arrow.grid(row=0, column=0, sticky=Tkinter.W)
  label = Label(root, text='Test frame:')
  label.grid(row=0, column=1, sticky=Tkinter.W)

  frame = Frame(root, bg='green', width=300, height=300)

  root.mainloop()


"""
======================COPYRIGHT/LICENSE START==========================

CrossLine.py: Creates cross lines canvas

Copyright (C) 2005 Wim Vranken (European Bioinformatics Institute)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license
 
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

- contact Wim Vranken (wim@ebi.ac.uk)
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

from memops.gui.Color import hexRepr
from memops.gui.Canvas import Canvas
from memops.gui.Frame import Frame

class CrossLine(Frame):

  def __init__(self, parent, color = 'black', side = 'both', canvas_bg = 'lightgrey', *args, **kw):

    self.color = color
    self.side = side

    apply(Frame.__init__, (self, parent) + args, kw)

    #self.grid_rowconfigure(0, weight=1)
    #self.grid_columnconfigure(1, weight=1)

    fill = Tkinter.X
    w = kw.get('width', 1)
    z = h = kw.get('height', 1)
    #self.canvas.config(height=h)
    self.canvas = Canvas(self, width=w, height=h, background = canvas_bg)

    self.cross = {'upper': [], 'lower': [], 'both': []}
    
    #
    # This is setup only... coordinates do not matter here
    #
    
    if self.side in ['both','upper']:
      self.cross['upper'].append(self.canvas.create_line(0, 0, 0, 0, fill=color))
      self.cross['upper'].append(self.canvas.create_line(0, 0, 0, 0, fill=color))

    self.cross['both'].append(self.canvas.create_line(0, 0, 0, 0, fill=color))

    if self.side in ['both','lower']:
      self.cross['lower'].append(self.canvas.create_line(0, 0, 0, 0, fill=color))
      self.cross['lower'].append(self.canvas.create_line(0, 0, 0, 0, fill=color))
    
    self.canvas.pack(expand=Tkinter.YES, fill=fill)

    #self.canvas.bind('<Configure>', self.resizeCallback)
    self.bind('<Configure>', self.resizeCallback)

  def resizeCallback(self, event):

    self.update_idletasks()

    #w = self.canvas.winfo_width()
    #h = self.canvas.winfo_height()
    w = self.winfo_width()
    h = self.winfo_height()
    
    if self.side == 'both':
      mh = h/2
    elif self.side == 'lower':
      mh = 0
    else:
      mh = h

    if self.side in ['both','upper']:
      self.canvas.coords(self.cross['upper'][0], 0, 0, w/3, mh)
      self.canvas.coords(self.cross['upper'][1], w*2/3, mh, w, 0)

    self.canvas.coords(self.cross['both'][0], w/3, mh, w*2/3, mh)
    
    if self.side in ['both','lower']:
      self.canvas.coords(self.cross['lower'][0], 0, h, w/3, mh)
      self.canvas.coords(self.cross['lower'][1], w*2/3, mh, w, h)

  # color is a tuple
  def setColor(self, color):

    (r, g, b) = color
    self.canvas.itemconfig(self.item, fill=hexRepr(r, g, b))

if __name__ == '__main__':

  def printEvent(event):
    print event.x, event.y

  from memops.gui.Label import Label

  root = Tkinter.Tk()

  root.grid_rowconfigure(0, weight=1)
  root.grid_columnconfigure(0, weight=1)

  sep = CrossLine(root, color='black', width = 20, height = 50, side = 'both')
  sep.grid(row=0, column=0, sticky=Tkinter.EW)

  root.mainloop()

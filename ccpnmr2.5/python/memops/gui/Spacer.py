
"""
======================COPYRIGHT/LICENSE START==========================

Spacer.py: <write function here>

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

from memops.gui.Color import hexRepr
from memops.gui.Canvas import Canvas
from memops.gui.Frame import Frame

class Spacer(Frame):

  def __init__(self, parent, relief='raised', *args, **kw):

    apply(Frame.__init__, (self, parent) + args, kw)

    #self.grid_rowconfigure(0, weight=1)
    #self.grid_columnconfigure(1, weight=1)

    self.bg      = self.cget('bg')
    self.relief  = relief
    self.grid_columnconfigure(0, weight=1)
    self.canvas  = Canvas(self, background=self.bg, width=1, height=6)
    self.canvas.grid(row=0, column=0, sticky=Tkinter.NSEW)

    #self.canvas.bind('<Configure>', self.resizeCallback)
    self.grid(sticky=Tkinter.EW)
    self.canvasL1 = []
    self.canvasL2 = []
    self.canvasL3 = []
    self.refresh  = 0
    self.event    = None
    self.bind('<Configure>', self.resizeAfter)
    
  def resizeAfter(self, event):
    
    self.event = event
    if self.refresh:
      return
    else:
      self.refresh = 1
      self.after_idle(lambda: self.resize(self.event))    
    
  def resize(self, event):
    
    m = 6
    w = event.width
    h = event.height
    N = int(h/m)
    
    self.canvas.config(width=w, height=h)
    if self.relief  == 'sunken':
      fill1 = 'grey65'
      fill2 = 'grey95'
    else:
      fill1 = 'grey95'
      fill2 = 'grey65'
    
    for i in range(N):
      y = 1 + (i/float(N))*h
      if i >= len(self.canvasL1):
        rect1 = self.canvas.create_rectangle(4, y+0, w-5, y+2,width=0,fill=fill1)
        rect2 = self.canvas.create_rectangle(5, y+1, w-4, y+3,width=0,fill=fill2)
        rect3 = self.canvas.create_rectangle(5, y+1, w-5, y+2,width=0,fill=self.bg)
        self.canvasL1.append(rect1)
        self.canvasL2.append(rect2)
        self.canvasL3.append(rect3)
      else: 
        self.canvas.coords(self.canvasL1[i], 4, y+0, w-5, y+2)
        self.canvas.coords(self.canvasL2[i], 5, y+1, w-4, y+3)
        self.canvas.coords(self.canvasL3[i], 5, y+1, w-5, y+2)

    if N > len(self.canvasL1):
      for i in range(N, len(self.canvasL1)):
        self.canvas.delete(self.canvasL1[i])
        self.canvas.delete(self.canvasL2[i])
        self.canvas.delete(self.canvasL3[i])
    
    self.refresh = 0

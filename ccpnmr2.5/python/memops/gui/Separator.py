
"""
======================COPYRIGHT/LICENSE START==========================

Separator.py: <write function here>

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

from memops.gui.Color  import hexRepr
from memops.gui.Canvas import Canvas
from memops.gui.Frame  import Frame

class Separator(Frame):
  
  refresh  = 0

  def __init__(self, parent, orient = Tkinter.HORIZONTAL, relief='solid', color='black',
               bgColor = 'grey', *args, **kw):

    if (orient == Tkinter.HORIZONTAL):
      size = kw.get('height', 1)
    else:
      size = kw.get('width', 1)
      
    if relief =='solid':
      kw['relief'] = 'flat'
    else:
      kw['relief'] = relief
    
    
    kw['borderwidth'] = (size-1)/2
    
    apply(Frame.__init__, (self, parent) + args, kw)

    self.config(bg=color)
    self.relief = relief

    if self.relief in ('solid','flat'):
      color2 = bgColor
    else:
      color2 = color
 
    self.grid_columnconfigure(0, weight=1)
    self.grid_rowconfigure(0, weight=1)
    self.innerFrame = Frame(self, borderwidth=0, bg=color2)
    self.innerFrame.grid(sticky = Tkinter.NSEW)     
    self.innerFrame.grid_columnconfigure(0, weight=1)
    self.innerFrame.grid_rowconfigure(0, weight=1)
 
  # color is a tuple
  def setColor(self, color, bgColor='grey'):

    if type(color) is type(()) :
      (r,g,b) = color
      color   = hexRepr(r,g,b)
      
    if type(bgColor) is type(()) :
      (r,g,b) = bgColor
      bgColor = hexRepr(r,g,b)
    
    self.config(bg = color)
    
    if self.relief in ('solid','flat'):
      self.innerFrame.config(bg = bgColor)
    else:
      self.innerFrame.config(bg = color)

if __name__ == '__main__':

  def printEvent(event):
    print event.x, event.y

  from Label import Label

  root = Tkinter.Tk()

  root.grid_rowconfigure(0, weight=1)
  root.grid_columnconfigure(0, weight=1)

  nrows = 5
  ncols = 4

  for j in range(nrows):
    for i in range(ncols):
      label = Label(root, text='label number '+str(j)+str(i))
      label.grid(row=2*j, column=2*i)

  for j in range(1,nrows):
    sep = Separator(root, orient=Tkinter.HORIZONTAL, color='red')
    sep.canvas.bind('<ButtonPress-1>', printEvent)
    sep.grid(row=2*j-1, column=0, columnspan=2*ncols-1, sticky=Tkinter.EW)

  for i in range(1,ncols):
    sep = Separator(root, orient=Tkinter.VERTICAL, color='green')
    sep.canvas.bind('<ButtonPress-1>', printEvent)
    sep.grid(row=0, column=2*i-1, rowspan=2*nrows-1, sticky=Tkinter.NS)

  root.mainloop()

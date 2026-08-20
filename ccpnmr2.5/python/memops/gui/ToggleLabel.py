
"""
======================COPYRIGHT/LICENSE START==========================

ToggleLabel.py: <write function here>

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
 
 
from memops.gui.Button import Button
from memops.gui.Frame import Frame
from memops.gui.ToggleArrow import ToggleArrow
 
class ToggleLabel(Frame):
 
  def __init__(self, parent, arrowSize = 12, outline='#5050b0', fill='#a0a0ff', font=None,
               isArrowClosed = True, callback = None, text = '', *args, **kw):

    apply(Frame.__init__, (self, parent) + args, kw)

    self.grid_rowconfigure(0, weight=1)
    self.grid_columnconfigure(1, weight=1)

    self.arrow = ToggleArrow(self, arrowSize=arrowSize, outline=outline, fill=fill,
                             isArrowClosed=isArrowClosed, callback=callback)
    self.arrow.grid(row=0, column=0, sticky=Tkinter.W)
    if font:
      self.label = Button(self, relief='flat',borderwidth=0, padx=1, pady=1, text=text, font=font, command=self.arrow.flipState)
    else:
      self.label = Button(self, relief='flat',borderwidth=0, padx=1, pady=1, text=text, command=self.arrow.flipState)
    self.label.grid(row=0, column=1, sticky=Tkinter.W)

  def arrowOn(self):
  
    if self.arrow.isArrowClosed:
      self.arrow.isArrowClosed = 0
      self.arrow.drawArrow()
  
  def arrowOff(self):
  
    if not self.arrow.isArrowClosed:
      self.arrow.isArrowClosed = 1
      self.arrow.drawArrow()

  def __getattr__(self, name):
 
    # dispatch everything not defined by ToggleLabel to arrow widget
 
    try:
      return getattr(self.__dict__['arrow'], name)
    except:
      raise AttributeError, "%s instance has no attribute '%s'" % (self.__class__.__name__, name)

if (__name__ == '__main__'):
 
  frame = None
 
  def callback(isClosed):
 
    if (isClosed):
      frame.grid_forget()
    else:
      frame.grid(row=1, column=0, sticky=Tkinter.NSEW)
 
  root = Tkinter.Tk()
 
  root.grid_rowconfigure(1, weight=1)
  root.grid_columnconfigure(0, weight=1)
 
  toggle_label = ToggleLabel(root, callback=callback, text='Test frame:',
                             arrowSize=20)
  toggle_label.grid(row=0, column=0, sticky=Tkinter.W)
 
  frame = Frame(root, bg='green', width=300, height=300)
 
  root.mainloop()

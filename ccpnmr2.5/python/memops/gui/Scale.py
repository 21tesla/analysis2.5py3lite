
"""
======================COPYRIGHT/LICENSE START==========================

Scale.py: <write function here>

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

from memops.gui.Base import Base
from memops.gui.Color import getIntRgb

# scale seems to be brain dead in Tkinter, for some reason need DoubleVar
# or IntVar to change scale
# if you have a variable to hand you can pass it in using kw['variable']

class Scale(Tkinter.Scale, Base):

  def __init__(self, parent, value=None, oddNumbers=False, docKey=None, tipText=None, *args, **kw):

    self.var = kw.get('variable', Tkinter.DoubleVar())
    kw['variable'] = self.var
    kw['bg'] = 'grey90'
    
    if oddNumbers:
      self._userCallback = kw.get('command')
      kw['command'] = self.odd
 
    Tkinter.Scale.__init__(self, parent, *args, **kw)
    Base.__init__(self, docKey=docKey, tipText=tipText)

    self.parent = parent

    if value is not None:
      self.set(value)

    self.determineFgs()

  def determineFgs(self):

    fg = self.cget('fg')
    bg = self.cget('bg')

    (r1, b1, g1) = getIntRgb(self, fg)
    (r2, b2, g2) = getIntRgb(self, bg)

    r = (r1 + r2) / 2
    g = (g1 + g2) / 2
    b = (b1 + b2) / 2

    self.enableFg = fg
    self.disableFg = '#%02x%02x%02x' % (r, g, b)

  def get(self):

    return self.var.get()

  def set(self, value):

    self.var.set(value)

  def disable(self):

    self.config(fg=self.disableFg, state=Tkinter.DISABLED)

  def enable(self):

    self.config(fg=self.enableFg, state=Tkinter.NORMAL)

  def odd(self, arg):

    value = self.var.get()
    value = 2*int(value/2)+1
    self.set(value)
    if self._userCallback:
      self._userCallback(str(value))


if (__name__ == '__main__'):

  from memops.gui.Button import Button

  d = 10
  scale = None

  def func():
    global d
    v = scale.get()
    if (v > (float(scale['to'])-d)):
      d = -d
    elif (v < (float(scale['from'])-d)):
      d = -d
    scale.set(v+d)

  root = Tkinter.Tk()
 
  def func(args):
    print args
 
  scale = Scale(root, orient=Tkinter.HORIZONTAL, value=50.0,
                from_=1.0, to=99, oddNumbers=True, command=func,
                tipText='Some advice')
  
  scale.grid()
  
  scale.set(3)
  
  button = Button(root, text='hit me', command=func)
  button.grid()
  button = Button(root, text='Enable', command=scale.enable)
  button.grid()
  button = Button(root, text='Disable', command=scale.disable)
  button.grid()

  root.mainloop()


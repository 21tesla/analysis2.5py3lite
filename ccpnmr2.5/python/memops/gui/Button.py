
"""
======================COPYRIGHT/LICENSE START==========================

Button.py: <write function here>

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

DEFAULTBACKGROUND = 'grey82'
DEFAULTMIDCOLOR = 'grey90'


# this extends the Tkinter.Button class with a few extra functions
class Button(Tkinter.Button, Base):

  def __init__(self, parent, grid=None, gridSpan=(1,1), sticky='w', 
               docKey=None, tipText=None, *args, **kw):

    if not docKey:
      docKey = kw.get('text')

    self.font = kw.get('font')
    if not kw.has_key('bg'):
      kw['bg'] = DEFAULTBACKGROUND

    if not kw.get('activebackground'):
      kw['activebackground'] = '#D0B0A0'

    if not kw.get('highlightthickness'):
      kw['highlightthickness'] = 1

    if not kw.get('highlightbackground'):
      kw['highlightbackground'] = DEFAULTMIDCOLOR

    if not kw.get('highlightcolor'):
      kw['highlightcolor'] = DEFAULTMIDCOLOR

    if kw.has_key('command'):
      self.buttonCommand = kw['command']
    else:
      self.buttonCommand = None

    Tkinter.Button.__init__(self, parent, *args, **kw)
    Base.__init__(self, docKey=docKey, tipText=tipText)

    if grid is not None:
      row, col = grid
      rowSpan, colSpan = gridSpan
      self.grid(row=row, column=col, rowspan=rowSpan,
                columnspan=colSpan, sticky=sticky)

    self.parent = parent
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

  def disable(self):

    self.config(fg=self.disableFg, state=Tkinter.DISABLED)

  def enable(self):

    self.config(fg=self.enableFg, state=Tkinter.NORMAL)

  def setText(self, text):

    self.config(text=text)

  def setState(self, state):

    if (state.lower() == Tkinter.DISABLED):
      self.disable()
    else:
      self.enable()

  def configure(self, **options):

    if options.has_key('command'):
      cmd = self.cget('command')
      if cmd and cmd in self._tclCommands: # o/w get memory leak until button destroyed
        self.deletecommand(str(cmd))
      self.buttonCommand = options['command']

    return Tkinter.Button.config(self, **options)

  #config = configure - Taken out for now, as it messes analysis up
  def config(self, **options):

    return self.configure(**options)

if __name__ == '__main__':

  root = Tkinter.Tk()
  
  def click():
    print "Clicked"
    root.destroy()
  
  #b = Button(root, text='+ - + -\n- + - +\n+ - + -\n- + - +', command=click)
  b = Button(root, text='Click Me', command=click,
             activebackground='red',
             activeforeground='white',
             tipText='Click for action')
  b.pack()
  
  root.mainloop()

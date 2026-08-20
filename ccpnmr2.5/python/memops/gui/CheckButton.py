
"""
======================COPYRIGHT/LICENSE START==========================

CheckButton.py: <write function here>

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

class CheckButton(Tkinter.Checkbutton, Base):

  def __init__(self, parent, text ='', selected=False,
               callback=None, grid=None, gridSpan=(1,1),
               sticky='w', docKey=None, tipText=None, *args, **kw):
   
    if not docKey:
      docKey = text

    self.font = kw.get('font')
    if not kw.has_key('bg'):
      kw['bg'] = 'grey90'
      
    if not kw.has_key('highlightbackground'):
      kw['highlightbackground'] = kw.get('bg') or 'grey90'
      
    if not kw.has_key('highlightcolor'):
      kw['highlightcolor'] = kw.get('bg') or 'grey90'
 
    if not kw.has_key('activebackground'):
      kw['activebackground'] = '#D0B0A0'
 
    if not kw.has_key('selectcolor'):
      kw['selectcolor'] = '#D0B0A0'
      
    Tkinter.Checkbutton.__init__(self, parent, *args, **kw)
    Base.__init__(self, docKey=docKey, tipText=tipText)

    self.callback = callback
    self.var = Tkinter.IntVar()
    self.config(text=text, variable=self.var)

    if (callback):
      self.config(command=self.doCallback)

    self.setSelected(selected)
 
    if grid is not None:
      row, col = grid
      rowSpan, colSpan = gridSpan
      self.grid(row=row, column=col, rowspan=rowSpan,
                columnspan=colSpan, sticky=sticky)

  def doCallback(self):

    self.callback(self.getSelected())

  def isSelected(self):

    return self.getSelected()

  def get(self):

    return self.getSelected()

  def getSelected(self):

    if self.var.get():
      return True
    else:
      return False

  def set(self, selected):

    self.setSelected(selected)

  def setSelected(self, selected):

    if selected:
      self.var.set(1)
    else:
      self.var.set(0)

  def toggle(self):

    self.setSelected(not self.getSelected())

if __name__ == '__main__':

  from memops.gui.Button import Button

  def get_me():
    print 'get_me:', c.getSelected()

  def toggle_me():
    c.toggle()

  def my_callback(selected):
    print 'my_callback:', selected

  root = Tkinter.Tk()

  c = CheckButton(root, 'check button', callback=my_callback,
                  tipText='What this does')
  c.grid()
  b = Button(root, text='get me', command=get_me)
  b.grid()
  b = Button(root, text='toggle me', command=toggle_me)
  b.grid()

  root.mainloop()

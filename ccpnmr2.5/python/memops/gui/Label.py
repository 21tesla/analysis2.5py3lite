
"""
======================COPYRIGHT/LICENSE START==========================

Label.py: <write function here>

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

# labels seem to be brain dead in Tkinter, for some reason need StringVar to change label
# if you have a StringVar to hand you can pass it in using kw['textvariable']
# if you want the label to be initialised to something use kw['text']

class Label(Tkinter.Label, Base):

  def __init__(self, parent, grid=None, gridSpan=(1,1), sticky='w',
              docKey=None, tipText=None, *args, **kw):

    if not docKey:
      docKey = kw.get('text')

    self.font = kw.get('font')
    self.var = kw.get('textvariable', Tkinter.StringVar())
    kw['textvariable'] = self.var
    kw['bg'] = 'grey90'
    
    Tkinter.Label.__init__(self, parent, *args, **kw)
    Base.__init__(self, docKey=docKey, tipText=tipText)

    self.parent = parent

    text = kw.get('text', '')

    self.set(text)
    
    if grid is not None:
      row, col = grid
      rowSpan, colSpan = gridSpan
      self.grid(row=row, column=col, rowspan=rowSpan,
                columnspan=colSpan, sticky=sticky)
    
  def get(self):

    return self.var.get()

  def set(self, text = ''):

    self.var.set(text)

if (__name__ == '__main__'):

  from memops.gui.Button import Button

  msg = 'hello world'
  count = 0

  def func():

    global count

    count = count + 1
    label.set(msg + ' ' + str(count))
    print label.get()

  root = Tkinter.Tk()
 
  label = Label(root, text='hello world', relief='raised')
  label.grid()
  button = Button(root, text='hit me', command=func)
  button.grid()

  root.mainloop()


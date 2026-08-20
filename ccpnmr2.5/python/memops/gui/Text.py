
"""
======================COPYRIGHT/LICENSE START==========================

Text.py: <write function here>

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

class Text(Tkinter.Text, Base):

  def __init__(self, parent, text='', grid=None, gridSpan=(1,1),
               sticky='w',docKey=None, tipText=None, *args, **kw):

    if (not kw.has_key('bg') and not kw.has_key('background')):
      kw['bg'] = 'white'

    Tkinter.Text.__init__(self, parent,  *args, **kw)
    Base.__init__(self, docKey=docKey, tipText=tipText)

    self.parent = parent
    self.append(text)
    
    if grid is not None:
      row, col = grid
      rowSpan, colSpan = gridSpan
      self.grid(row=row, column=col, rowspan=rowSpan,
                columnspan=colSpan, sticky=sticky)

  def append(self, text):

    self.insert(Tkinter.END, str(text))

  def clear(self):

    self.delete('0.0', Tkinter.END)

  def getText(self):

    # TBD: str below is used to convert unicode to string
    # (in Tcl/Tk 8.5 it looks like everything is unicode)
    ###return str(self.get('0.0', Tkinter.END))
    return self.get('0.0', Tkinter.END).encode('utf-8')

  def setText(self, text):

    self.clear()
    self.append(text)

  def clearSelection(self):

    self.tag_remove(Tkinter.SEL, '1.0', Tkinter.END)

  def setSelection(self, start, end = Tkinter.END, clear = True):

    if (clear):
      self.clearSelection()

    self.tag_add(Tkinter.SEL, start, end)

if __name__ == '__main__':

  from memops.gui.Button import Button
  from memops.gui.Frame import Frame
  from memops.gui.Util import createQuitButton

  s = 7

  def addText():

    global s

    text.append(s)
    s = s + 10

    print 'text = ' + text.getText()

  def setText():

    global s

    text.setText(s)
    s = s + 10

    print 'text = ' + text.getText()

  def clearText():

    text.clear()

  root = Tkinter.Tk()

  frame = Frame(root)
  frame.pack(side=Tkinter.TOP)

  button = Button(frame, text='add text', command=addText)
  button.pack(side=Tkinter.LEFT)

  button = Button(frame, text='set text', command=setText)
  button.pack(side=Tkinter.LEFT)

  button = Button(frame, text='clear text', command=clearText)
  button.pack(side=Tkinter.LEFT)

  button = createQuitButton(frame)
  button.pack(side=Tkinter.LEFT)

  text = Text(root, width=20, height=10, tipText='Helpful tip')
  text.pack(side=Tkinter.TOP, expand=Tkinter.YES, fill=Tkinter.BOTH)

  root.mainloop()

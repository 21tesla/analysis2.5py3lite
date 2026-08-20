
"""
======================COPYRIGHT/LICENSE START==========================

ScrolledText.py: <write function here>

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


from memops.gui.Frame import Frame
from memops.gui.Text import Text

class ScrolledText(Frame):

  def __init__(self, parent, text = '',
               width = 60, height = 10,
               xscroll = True, yscroll = True,
               *args, **kw):

    apply(Frame.__init__, (self, parent) + args, kw)

    self.grid_rowconfigure(0, weight=1)
    self.grid_columnconfigure(0, weight=1)

    self.text_area = Text(self, width=width, height=height)
    self.text_area.grid(row=0, column=0, sticky=Tkinter.NSEW)

    if xscroll:
      self.text_area.config(wrap=Tkinter.NONE)
      xscrollbar = Tkinter.Scrollbar(self, orient=Tkinter.HORIZONTAL, bd=1)
      xscrollbar.config(command=self.text_area.xview)
      self.text_area.config(xscrollcommand=xscrollbar.set)
      xscrollbar.grid(row=1, column=0, sticky=Tkinter.EW)
    else:
      self.text_area.config(wrap=Tkinter.WORD)

    if yscroll:
      yscrollbar = Tkinter.Scrollbar(self, orient=Tkinter.VERTICAL, bd=1)
      yscrollbar.config(command=self.text_area.yview)
      self.text_area.config(yscrollcommand=yscrollbar.set)
      yscrollbar.grid(row=0, column=1, sticky=Tkinter.NS)

    self.text_area.append(text)

    self.size = self.text_area.size # otherwise get Frame size called

  def __getattr__(self, name):

    # dispatch everything not defined by ScrolledText to text widget

    try:
      return getattr(self.__dict__['text_area'], name)
    except:
      raise AttributeError, "ScrolledText instance has no attribute '%s'" % name

  def setState(self, state = Tkinter.NORMAL):

    assert state in (Tkinter.NORMAL, Tkinter.DISABLED), 'state = %s' % state

    self.text_area.config(state=state)

  def setCursor(self, cursor = ''):

    self.text_area.config(cursor=cursor)
  
  def setText(self, text=''):
    self.text_area.setText(text)

if (__name__ == '__main__'):

  from memops.gui.Button import Button
  from memops.gui.Util import createQuitButton

  def getText():

    print 'getText:', text.getText()

  root = Tkinter.Tk()

  frame = Frame(root)
  frame.pack(side=Tkinter.TOP)

  button = Button(frame, text='get text', command = getText)
  button.pack(side=Tkinter.LEFT)

  button = createQuitButton(frame)
  button.pack(side=Tkinter.LEFT)

  text = ScrolledText(root, width=60, height=10,
                      #text='in the beginning')
                      text='in the beginning', xscroll=False)
  text.pack(side=Tkinter.TOP, expand=Tkinter.YES, fill=Tkinter.BOTH)

  root.mainloop()

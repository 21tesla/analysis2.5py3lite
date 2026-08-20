
"""
======================COPYRIGHT/LICENSE START==========================

BooleanPulldownMenu.py: <write function here>

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

from memops.gui.PulldownMenu import PulldownMenu


class BooleanPulldownMenu(PulldownMenu):

  booleanToIndex = { True: 0, False: 1 }
  indexToBoolean = { 0: True, 1: False }

  def __init__(self, parent, callback = None, entries = None,
               selected = True, *args, **kw):

    self.booleanCallback = callback

    if (entries):
      assert len(entries) == 2
      kw['entries'] = entries
    else:
      kw['entries'] = ['True', 'False']

    kw['selected_index'] = self.booleanToIndex[selected]
    kw['callback'] = self.doCallback

    apply(PulldownMenu.__init__, (self, parent) + args, kw)

  def doCallback(self, selected_index, selected):

    if (self.booleanCallback):
      if (selected == None):
        self.booleanCallback(None)
      else:
        self.booleanCallback(self.indexToBoolean[selected_index])

  def setSelected(self, value):

    PulldownMenu.setSelectedIndex(self, self.booleanToIndex[value])

  def getSelected(self):

    selected_index = PulldownMenu.getSelectedIndex(self)

    return self.indexToBoolean[selected_index]

if __name__ == '__main__':
 
  import sys
  import Tkinter
  from memops.gui.Button import Button

  popup_select = None

  def callback(text):
    global popup_select
    if (popup_select):
      print 'callback: selected =', popup_select.getSelected()

  def toggle():
    global popup_select
    if (popup_select.getSelected()):
      popup_select.setSelected(False)
    else:
      popup_select.setSelected(True)

  root = Tkinter.Tk()

  popup_select = BooleanPulldownMenu(root, callback=callback, selected=False)
  popup_select.grid(row=0, column=0, columnspan=2)

  button = Button(root, text='toggle extra state', command=toggle)
  button.grid(row=1, column=0)

  button = Button(root, text='quit', command=sys.exit)
  button.grid(row=1, column=1)

  root.mainloop()

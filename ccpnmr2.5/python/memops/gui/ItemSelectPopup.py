
"""
======================COPYRIGHT/LICENSE START==========================

ItemSelectPopup.py: <write function here>

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
import os

import Tkinter


from memops.gui.BasePopup import BasePopup
from memops.gui.Label import Label
from memops.gui.PulldownList import PulldownList
from memops.gui.ButtonList import ButtonList

class ItemSelectPopup(BasePopup):

  def __init__(self, parent, entries, label='', message='',
               select_text='Select', default=0, *args, **kw):
 
    self.entries = entries
    self.label = label
    self.message = message
    self.select_text = select_text
    self.default = default
    self.item = None

    kw['title'] = 'Select Item' 
    kw['transient'] = True
    kw['modal'] = True
    BasePopup.__init__(self, parent=parent, *args, **kw)

  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(1, weight=1)

    row = 0
    if self.message:
      label = Label(guiFrame, text=self.message, gridSpan=(1,2), grid=(row, 0))
      row += 1

    if self.label:
      label = Label(guiFrame, text=self.label, grid=(row,0))

    self.itemMenu = PulldownList(guiFrame, texts=self.entries, objects=self.entries,
                                 index=self.default, grid=(row,1))

    row += 1
    texts = [ self.select_text, 'Cancel' ]
    commands = [ self.ok, self.cancel]
    buttons = ButtonList(guiFrame, texts=texts, commands=commands,
                         gridSpan=(1,2), grid=(row, 0))

  def cancel(self):
    
    self.destroy()
    
    return None

  def apply(self):

    self.item = self.itemMenu.getText()

    return True

if __name__ == '__main__':

  root = Tkinter.Tk()

  popup = ItemSelectPopup(root, entries=['a', 'b', 'c', 'd'], label='Entries: ',
                          message='Select one of the entries')
  print 'item =', popup.item
  popup.destroy()

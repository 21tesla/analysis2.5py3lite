
"""
======================COPYRIGHT/LICENSE START==========================

NamingSystemPopup.py: GUI window for selecting naming system

Copyright (C) 2005 Wim Vranken (European Bioinformatics Institute)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license
 
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

- contact Wim Vranken (wim@ebi.ac.uk)
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

from memops.universal.Io import joinPath

from memops.gui.Label import Label
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.PulldownMenu import PulldownMenu

from ccpnmr.format.general.Io import getHelpUrlDir

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

class NamingSystemPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'NamingSystem.html')

  def __init__(self, parent, namingSysList = None, namingSysDict = None):

    # Constructor doesn't do much except call body
    # The parent is self.parent (parent of the popup)
   
    self.namingSys = None
    
    self.namingSysList = namingSysList
    self.namingSysDict = namingSysDict

    # modal = true means that it won't continue unless this one returns value
    TemporaryBasePopup.__init__(self, parent=parent, title='Choose namingSys', modal=False, transient=True)

  def body(self, master):
    
    #
    # Popup gui
    #

    # First row
    row = 0
    label = Label(master, text= "Choose a naming system:")
    label.grid(row=row, column=0, sticky=Tkinter.W)

    row = row + 1
    label = Label(master, text= "(% matches in brackets)")
    label.grid(row=row, column=0, sticky=Tkinter.W)

    row = row + 1
    self.menu = PulldownMenu(master, entries = self.namingSysList)
    self.menu.grid(row=row, column=0, sticky=Tkinter.EW)
   
    row = row + 1
    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, dismiss_text = 'Exit', help_url=self.help_url)
    buttons.grid(row=row, column=0)

  def apply(self):
    
    self.namingSysName = self.namingSysDict[self.menu.getSelected()]
    
    return True

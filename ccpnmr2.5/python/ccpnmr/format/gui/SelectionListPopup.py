
"""
======================COPYRIGHT/LICENSE START==========================

SelectionListPopup.py: GUI window for setting up generic selection list

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

from ccpnmr.format.general.Io import getHelpUrlDir

from memops.universal.Io import joinPath

from memops.gui.PulldownList import PulldownList
from memops.gui.Label import Label
from memops.gui.Util import createDismissHelpButtonList, createHelpButtonList

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

class SelectionListPopup(TemporaryBasePopup):
 
  def __init__(self, parent, selectionList, title = 'Select', text = 'Select', topText = None, dismissText = None, selected = None, selectionDict = None, urlFile = None, dismissButton = True, modal = False):
  
    self.selectionList = selectionList
    self.selectionDict = selectionDict
    
    self.text = text
    
    self.dismissButton = dismissButton
    
    if dismissButton:
      if dismissText:
        self.dismissText = dismissText
      else:
        self.dismissText = 'dismiss'
      
    self.topText = topText
    self.isSelected = None

    if not selected:
      self.selectedIndex = 0
    else:
      self.selectedIndex = self.selectionList.index(selected)
      
    if urlFile:
      self.help_url = joinPath(getHelpUrlDir(),urlFile + '.html')
    else:
      self.help_url = None
          
    TemporaryBasePopup.__init__(self,parent = parent, title = title, modal = modal, transient=True)
 
  def body(self, master):
    
    #
    # Popup window
    #

    row = 0
    
    if self.topText:
      label = Label(master, text= self.topText)
      label.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.EW)
    
      row = row + 1

    label = Label(master, text= self.text)
    label.grid(row=row, column=0, sticky=Tkinter.EW)

    self.menu = PulldownList(master, texts = self.selectionList, index = self.selectedIndex)
    self.menu.grid(row=row, column=1, sticky=Tkinter.E, ipadx = 20)
 
    row = row + 1
    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    
    if self.dismissButton:
      buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, dismiss_text = self.dismissText, help_url=self.help_url)
    else:
      buttons = createHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)

    buttons.grid(row=row, column=0, columnspan = 3)
   

  def apply(self):
        
    self.isSelected = self.menu.getText()
    
    if self.selectionDict:
      self.selection = self.selectionDict[self.isSelected]
    else:
      self.selection = self.isSelected

    return True

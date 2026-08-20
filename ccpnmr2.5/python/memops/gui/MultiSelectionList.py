  
"""
======================COPYRIGHT/LICENSE START==========================

MultiSelectionList.py: GUI window for setting up generic multi selection list

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

from memops.universal.Io import joinPath

from memops.gui.ScrolledListbox import ScrolledListbox
from memops.gui.Label import Label
from memops.gui.Util import createDismissHelpButtonList

from memops.gui.BasePopup import BasePopup

class MultiSelectionListSkeleton:
 
  def __init__(self, parent, selectionList, selectionDict = {}, selectedItems = [], title = 'MultiSelect', text = 'MultiSelect', endText = None, dismissText = 'dismiss', urlFile = None, modal = False):
  
    self.selectionList = selectionList
    self.selectionDict = selectionDict
    self.selectedItems = selectedItems
    self.text = text
    self.endText = endText
    self.dismissText = dismissText
    self.isSelectedList = None
    
    self.setUrlFile(urlFile)
    
    self.startBasePopup(parent,title,modal,True)
    
  def startBasePopup(self,parent,title,modal,True):
  
    pass
     
  def setUrlFile(self,urlFile):
  
    self.help_url = None
 
  def body(self, master):
    
    #
    # Popup window
    #

    row = 0
    label = Label(master, text= self.text)
    label.grid(row=row, column=0, sticky=Tkinter.EW)

    row = row + 1
    self.multiSelect = ScrolledListbox(master,
                                       width = 50,
	                                      height = 5,
	                                      selectmode = Tkinter.MULTIPLE,
                                       initial_list = self.selectionList)
    self.multiSelect.grid(row=row, column=0, sticky=Tkinter.E, ipadx = 20)
    master.grid_columnconfigure(0, weight=1)
    master.grid_rowconfigure(row, weight=1)
    
    if self.selectedItems:
      self.multiSelect.setSelectedItems(self.selectedItems)
 
    if self.endText:
      row = row + 1
      label = Label(master, text= self.endText)
      label.grid(row=row, column=0, sticky=Tkinter.EW)
 
    row = row + 1
    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, dismiss_text = self.dismissText, help_url=self.help_url)
    buttons.grid(row=row, column=0)
   

  def apply(self):
    
    self.isSelectedList = []
    
    for selectedItem in self.multiSelect.getSelectedItems():
    
      if self.selectionDict.has_key(selectedItem):
        self.isSelectedList.append(self.selectionDict[selectedItem])
      else:
        self.isSelectedList.append(selectedItem)
        
    return True


class MultiSelectionListPopup(MultiSelectionListSkeleton,BasePopup):

  def startBasePopup(self,parent,title,modal,True):
  
    BasePopup.__init__(self,parent = parent, title = title, modal = modal, transient=True)

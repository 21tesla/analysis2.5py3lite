
"""
======================COPYRIGHT/LICENSE START==========================

AtomSelectPopup.py: GUI window for selecting atoms to link to a resonance

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
import string

from ccpnmr.format.general.Io import getHelpUrlDir

from memops.universal.Io import joinPath

from memops.gui.Label import Label
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.PulldownMenu import PulldownMenu

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

class AtomSelectPopup(TemporaryBasePopup):

  def __init__(self, parent, matchName, chemAtomList, chemAtomDict, propagateList, ccpCode, shifts, matchType = 'Resonance', matchInfoPopup = None, selectedIndex = 4, title = 'Select resonance atom match', helpUrl = 'AtomSelect.html', headerLines = None, modal = False):

    self.help_url = joinPath(getHelpUrlDir(),helpUrl)

    # Constructor doesn't do much except call body
    # The parent is self.parent (parent of the popup)
   
    self.chemAtomOrSets = None
    self.propagate = None
    
    self.matchName = matchName
    self.matchType = matchType
    self.matchInfoPopup = matchInfoPopup
    self.chemAtomList = chemAtomList
    self.chemAtomDict = chemAtomDict
    self.propagateList = propagateList
    self.ccpCode = ccpCode
    self.shifts = shifts
    self.guiParent = parent
    self.headerLines = headerLines
    
    self.selectedIndex = selectedIndex

    # modal = true means that it won't continue unless this one returns value
    TemporaryBasePopup.__init__(self, parent=parent, title=title, modal=modal, transient=True)

  def showAllAtoms(self):

    self.chemAtomList = self.chemAtomDict.keys()
    
    self.chemAtomList.sort()

    self.menu.replace(self.chemAtomList)

  def body(self, master):    
        
    #
    # Popup window
    #

    row = 0

    if self.headerLines:
      for headerLine in self.headerLines:
        label = Label(master, text= headerLine)
        label.grid(row=row, column=0, columnspan = 3, sticky=Tkinter.EW)
        row += 1

    label = Label(master, text= "Residue type '%s'" % self.ccpCode)
    label.grid(row=row, column=0, columnspan = 3, sticky=Tkinter.EW)

    row = row + 1
    
    matchNameString = "[%s]" % self.matchName
    
    if self.matchInfoPopup:
      label = Label(master, text= "%s:" % self.matchType)
      label.grid(row=row, column=0, sticky=Tkinter.E)
      button = Tkinter.Button(master, text= matchNameString, command = lambda : self.matchInfoPopup(self,self.matchName), fg = 'blue')
      button.grid(row=row, column=1, columnspan = 2)
    else:
      label = Label(master, text= "%s:" %self.matchType)
      label.grid(row=row, column=0, sticky=Tkinter.E)
      label = Label(master, text= matchNameString, fg = 'blue')
      label.grid(row=row, column=1, columnspan = 2)
    
    if self.shifts:
      row = row + 1
      label = Label(master, text= "(shifts: %s ppm)" % self.shifts)
      label.grid(row=row, column=1, columnspan = 2, sticky=Tkinter.EW)

    row = row + 1
    
    label = Label(master, text= "Pick atom:")
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.menu = PulldownMenu(master, entries = self.chemAtomList)
    self.menu.grid(row=row, column=1,  columnspan = 2, ipadx = 20)
    
    row = row + 1
    
    label = Label(master, text= "Propagate\nmapping to:")
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.propagateMenu = PulldownMenu(master, entries = self.propagateList, selected_index = self.selectedIndex)
    self.propagateMenu.grid(row=row, column=1,  columnspan = 2)

    row = row + 1
    texts = [ 'OK' , 'Show all atoms']
   
    commands = [ self.ok, self.showAllAtoms]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, dismiss_text = 'Do not link', dismiss_cmd = self.doNotLink, help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = 3)

  def apply(self):
  
    selectionText = self.menu.getSelected()
    
    if selectionText == '<None>':
      return False
    
    self.chemAtomOrSets = self.chemAtomDict[selectionText]
    
    self.propagate = self.propagateMenu.getSelectedIndex()
    
    return True 
  
  def doNotLink(self):
  
    self.propagate = self.propagateMenu.getSelectedIndex()

    self.close()
    

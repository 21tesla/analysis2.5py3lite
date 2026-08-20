
"""
======================COPYRIGHT/LICENSE START==========================

ResonanceGroupPopup.py: GUI window for grouping resonances together

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

from memops.universal.Io import joinPath

from memops.gui.BasePopup import BasePopup
from memops.gui.Label import Label
from memops.gui.Util import createHelpButtonList
from memops.gui.CheckButton import CheckButton

from ccpnmr.format.general.Io import getHelpUrlDir

class ResonanceGroupPopup(BasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'ResonanceGroup.html')

  def __init__(self, parent, resonanceNameList, resonanceLineDict, ccpCode, title = 'Group resonances'):
   
    self.resNameGroups = None
     
    self.resonanceNameList = resonanceNameList
    self.resonanceLineDict = resonanceLineDict
    self.checkButtons = []
    self.ccpCode = ccpCode
        
    BasePopup.__init__(self, parent=parent, title=title, modal=True, transient=True)

  def body(self, master):
      
    master.grid_columnconfigure(0, weight = 1)
    master.grid_columnconfigure(1, weight = 1)
  
    self.geometry('600x400')

    #
    # Popup window
    #

    row = 0
    label = Label(master, text= "Residue type '%s'" % self.ccpCode)
    label.grid(row=row, column=0, sticky=Tkinter.EW)
      
    row = row + 1
    label = Label(master, text= "Please select resonances that belong together (remaining ones are also grouped)")
    label.grid(row=row, column=0, sticky=Tkinter.EW)
     
    for resName in self.resonanceNameList:

      row = row + 1
      
      label = Label(master, text= "Resonance '%s'" % resName)
      label.grid(row=row, column=0, sticky=Tkinter.EW)     
      
      self.checkButtons.append(CheckButton(master))
      self.checkButtons[-1].grid(row=row, column=1, sticky=Tkinter.EW)     

    row = row + 1
    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, column=0)

    for i in range(row):
      master.grid_rowconfigure(i, weight = 1)

  def apply(self):
    
    self.resNameGroups = []
    self.resNameGroups.append([])
    self.resNameGroups.append([])
    
    for i in range(0,len(self.checkButtons)):
      
      checkButton = self.checkButtons[i]
      
      if checkButton.isSelected():
        
        groupIndex = 0
        
      else:
      
        groupIndex = 1

      self.resNameGroups[groupIndex].append(self.resonanceLineDict[i])
        
    return True
  

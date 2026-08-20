
"""
======================COPYRIGHT/LICENSE START==========================

ChemCompSelectionPopup.py: GUI window to handle chemComp selection

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

import os

from memops.universal.Io import joinPath

from memops.gui.Util import createDismissHelpButtonList
from memops.gui.MessageReporter import showInfo

from ccpnmr.format.general.Io import getHelpUrlDir

from ccp.gui.ChemCompFrame import ChemCompFrame

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

class ChemCompSelectionPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'ChemCompSelection.html')

  def __init__(self, parent, project, **keywds):
  
    self.project = project
    self.keywds = keywds
    
    origCode = self.keywds['origCode']
    del self.keywds['origCode']

    TemporaryBasePopup.__init__(self, parent=parent, title="Project '%s': " % project.name + "ChemComp selection for code '%s'" % origCode, modal=False, transient=False, width = 200)

  def body(self, master):
    
    self.frame = ChemCompFrame(master,self.project,**self.keywds)
    self.frame.pack(side=Tkinter.TOP, expand=Tkinter.YES, fill=Tkinter.BOTH)
    
    texts = ['OK']
    commands = [self.ok]
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, dismiss_text = 'Cancel', help_url=self.help_url)
    
    buttons.pack() 

  def apply(self):
    
    if self.frame.getSelectedChemComp():
    
      return True
      
    else:
      
      showInfo('No selection','Please select a chemCompVar from the bottom half of the window, or press Cancel to quit.')
      
      return False

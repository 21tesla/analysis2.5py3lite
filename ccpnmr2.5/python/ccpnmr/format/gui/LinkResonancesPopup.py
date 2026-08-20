
"""
======================COPYRIGHT/LICENSE START==========================

LinkResonancesPopup.py: GUI window for setting parameters to start linkResonances script

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

from ccpnmr.format.general.Io import getHelpUrlDir

from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.Label import Label
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.MessageReporter import showInfo
from memops.gui.CheckButton import CheckButton

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

class LinkResonancesPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'LinkResonances.html')
 
  def __init__(self, parent, formatNamesList, formatNamesDict, title = 'Link resonances setup'):
  
    self.guiParent = parent
    self.formatNamesList = formatNamesList
    self.formatNamesDict = formatNamesDict
    self.status = False

    TemporaryBasePopup.__init__(self,parent = parent, title = title, modal = False, transient=False)
 
  def body(self, master):      
    
    row = 0

    label = Label(master, text = 'File format reference:')
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.menu = PulldownMenu(master, entries = self.formatNamesList)
    self.menu.grid(row=row, column=1, sticky=Tkinter.E, ipadx = 20)
        
    row = row + 1
        
    label = Label(master, text = 'Try to link unrecognized atom names')
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.linkRes = CheckButton(master, selected = True)
    self.linkRes.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1
        
    label = Label(master, text = 'Try IUPAC names (as backup)')
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.iupacNames = CheckButton(master)
    self.iupacNames.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1

    label = Label(master, text = 'Use ambiguous name information')
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.useAmb = CheckButton(master)
    self.useAmb.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1

    label = Label(master, text = 'Complete stereospecific assignment')
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.allStereo = CheckButton(master)
    self.allStereo.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1

    label = Label(master, text = 'Force shift merges')
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.shiftMerge = CheckButton(master)
    self.shiftMerge.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1

    label = Label(master, text = 'Status other atom for all single prochiral atoms')
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.allSingleProchiral = PulldownMenu(master, entries = ['Unknown','Same information','Always ignore'])
    self.allSingleProchiral.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1

    label = Label(master, text = 'Status other atom for all possibly equivalent single atoms')
    label.grid(row=row, column=0, sticky=Tkinter.E)
        
    self.allSinglePossEquiv = PulldownMenu(master, entries = ['Unknown','Always equivalent','Always ignore'])
    self.allSinglePossEquiv.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1

    label = Label(master, text = 'Automatically connect ambiguous resonances to stereospecific ones')
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.connStereo = CheckButton(master)
    self.connStereo.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1
        
    label = Label(master, text = 'Make methyl protons equivalent')
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.methylsEquivalent = CheckButton(master, selected = True)
    self.methylsEquivalent.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1

    label = Label(master, text = 'Use minimal number of popups')
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.minimalPopups = CheckButton(master)
    self.minimalPopups.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1
        
    label = Label(master, text = 'Verbose output')
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.verbose = CheckButton(master, selected = True)
    self.verbose.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1

    label = Label(master, text = 'Warning: it is recommended you save your project first,\nin case linkResonances is interrupted (this might corrupt the data).', fg = 'red')
    label.grid(row=row, column=0, columnspan = 2, ipady = 15, sticky=Tkinter.EW)

    row = row + 1

    texts = [ 'Link resonances to atoms' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = 3)
  
      
  def apply(self):

    self.assignFormat = self.formatNamesDict[self.menu.getSelected()]

    self.globalStereoAssign = self.allStereo.isSelected()
    self.useIupacMatching = self.iupacNames.isSelected()
    self.useLinkResonancePopup = self.linkRes.isSelected()
    self.minimalPrompts = self.minimalPopups.isSelected()
    self.verbose = self.verbose.isSelected()
    self.forceShiftMerge = self.shiftMerge.isSelected()
    self.useAmbiguity = self.useAmb.isSelected()
    self.autoConnectStereo = self.connStereo.isSelected()
    self.makeMethylsEquivalent = self.methylsEquivalent.isSelected()
    
    self.setSingleProchiral = self.allSingleProchiral.getSelected()
    self.setSinglePossEquiv = self.allSinglePossEquiv.getSelected()
    
    self.status = True
    
    return True

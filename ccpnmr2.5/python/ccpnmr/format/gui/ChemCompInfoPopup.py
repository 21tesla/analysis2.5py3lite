
"""
======================COPYRIGHT/LICENSE START==========================

ChemCompInfoPopup.py: GUI window for defining chemComp information

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

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup
from memops.gui.Label import Label
from memops.gui.Entry import Entry
from memops.gui.CheckButton import CheckButton
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.Util import createHelpButtonList
from memops.gui.CheckButton import CheckButton

from ccpnmr.format.general.Io import getHelpUrlDir

class ChemCompInfoPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'ChemCompInfo.html')

  def __init__(self, parent, formula, bondNumber, chemCompInfo):
   
    self.formula = formula
    self.bondNumber = bondNumber
    self.chemCompInfo = chemCompInfo

    self.updated = 0
    
    self.nonEntryAttributes = {
    
      'molType': (PulldownMenu,['protein','DNA','RNA','carbohydrate','other']),
      'hasStdChirality': (CheckButton,None)
    
    }
        
    TemporaryBasePopup.__init__(self, parent=parent, title='ChemComp creation', modal=False, transient=True)

  def body(self, master):
      
    #
    # Popup window
    #
    
    self.widgets = []

    row = 0
    label = Label(master, text= "ChemComp formula '%s'" % self.formula)
    label.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.EW)
      
    row = row + 1
    label = Label(master, text= "Number of bonds: %d" % self.bondNumber)
    label.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.EW)
    
    #
    # Show relevant attributes...
    #
    
    for chemCompAttrInfo in self.chemCompInfo:
    
      row = row + 1
      
      attrName = chemCompAttrInfo[0]

      label = Label(master, text = attrName)
      label.grid(row=row, column=0, sticky=Tkinter.EW)     
      
      if attrName in self.nonEntryAttributes:
      
        widgetInfo = self.nonEntryAttributes[attrName]
        
        if widgetInfo[0] == PulldownMenu:
        
          self.widgets.append(PulldownMenu(master, entries = widgetInfo[1], selected_index = widgetInfo[1].index(chemCompAttrInfo[1])))

        elif widgetInfo[0] == CheckButton:

          self.widgets.append(CheckButton(master, selected = widgetInfo[1]))
      
      else:
      
        text = chemCompAttrInfo[1]
        if not text:
          text = ''

        self.widgets.append(Entry(master, text = text))
      
      self.widgets[-1].grid(row=row, column=1, sticky=Tkinter.EW)     

    row = row + 1
    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = 2)

  def apply(self):
  
    self.updated = 1
    
    for widgetNum in range(0,len(self.chemCompInfo)):
      
      chemCompAttrInfo = self.chemCompInfo[widgetNum]
    
      attrName = chemCompAttrInfo[0]
      
      if attrName in self.nonEntryAttributes:
      
        widgetInfo = self.nonEntryAttributes[attrName]
        
        if widgetInfo[0] in [PulldownMenu,CheckButton]:
        
          value = self.widgets[widgetNum].getSelected()
      
      else:
      
        value = self.widgets[widgetNum].get()
        
        if value == '':
          value = None
    
      self.chemCompInfo[widgetNum] = (attrName,value)
        
    return True
  

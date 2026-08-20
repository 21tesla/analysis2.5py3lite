
"""
======================COPYRIGHT/LICENSE START==========================

ColumnInfoSelectionPopup.py: GUI window for selecting which information is in which column

Copyright (C) 2006 Wim Vranken (European Bioinformatics Institute)

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
from memops.universal.Constants import dirsep
from memops.universal.Io import joinPath

from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.Label import Label
from memops.gui.Util import createHelpButtonList
from memops.gui.MessageReporter import showError

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

class ColumnInfoSelectionPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'ColumnInfoSelection.html')

  def __init__(self, parent, columnData, columnContents, title = None):
    
    self.columnData = columnData
    self.columnContents = columnContents
    
    self.columnLen = len(self.columnData)
    
    if not title:
      title = "Column information selection"
    print title
    
    TemporaryBasePopup.__init__(self,parent = parent, title = title, modal = False, transient=True)
 
  def body(self, master):

    #
    # Initialize
    #
    
    contentList = []
    
    for (fullName,attrName,func) in self.columnContents:
      contentList.append(fullName)   
                    
    #
    # Popup info
    #
    
    #
    # Header labels
    #

    row = 0
    
    label = Label(master, text= "Column number")
    label.grid(row=row, column=0, sticky=Tkinter.EW)
    
    label = Label(master, text= "Column data range")
    label.grid(row=row, column=1, sticky=Tkinter.EW)
    
    label = Label(master, text= "Column content selection")
    label.grid(row=row, column=2, sticky=Tkinter.EW)

    #
    # Selection per column
    #
    
    self.columnInfoWidgets = []
    
    for colIndex in range(0,self.columnLen):

      row = row + 1

      label = Label(master, text= str(colIndex + 1))
      label.grid(row=row, column=0, sticky=Tkinter.EW)

      label = Label(master, text= "%s - %s" % (str(self.columnData[colIndex][0]), str(self.columnData[colIndex][1])))
      label.grid(row=row, column=1, sticky=Tkinter.EW)

      self.columnInfoWidgets.append(PulldownMenu(master, entries = contentList))
      self.columnInfoWidgets[-1].grid(row=row, column=2, sticky=Tkinter.W, ipadx = 20)
    
    row = row + 1
    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = 3)
   

  def apply(self):
    
    columnInfoLong = []
    
    for colIndex in range(0,self.columnLen):
      
      columnInfoLongName = self.columnInfoWidgets[colIndex].getSelected()
      
      if columnInfoLongName != 'None' and columnInfoLongName in columnInfoLong:
        showError("Unique names","Each column needs a unique content - %s repeated twice." % columnInfoLongName)
        return False
      
      columnInfoLong.append(columnInfoLongName)
      
    if columnInfoLong.count('None') == len(columnInfoLong):
      showError("Need identifiers","Need content for at least one column.")
      return False

    self.columnConvert = []
    
    for i in range(self.columnLen):
      for (longName,attrName,function) in self.columnContents:
        if longName == columnInfoLong[i]:
          if attrName:
            self.columnConvert.append((attrName,function))
          else:
            self.columnConvert.append(None)
          break
    
    return True

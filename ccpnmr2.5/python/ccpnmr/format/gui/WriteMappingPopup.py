
"""
======================COPYRIGHT/LICENSE START==========================

WriteMappingPopup.py: GUI window for setting parameters to write ccpn atom/resonance mapping

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

from memops.gui.BasePopup import BasePopup
from memops.gui.Label import Label
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.ScrolledListbox import ScrolledListbox
from memops.gui.MessageReporter import showInfo, showError

from ccpnmr.format.gui.FormatFilePopup import FormatFilePopup

from ccpnmr.format.general.Io import getHelpUrlDir
from ccpnmr.format.general.Constants import allFormatsDict
from ccpnmr.format.general.Constants import ccpNmr_kw

#
# Get script to write mapping file for resonances -> atoms
#

from ccpnmr.format.process.writeMappingFile import writeMappingFile

class WriteMappingPopup(BasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'WriteMapping.html')

  def __init__(self, parent, project):

    self.project = project
    self.selectedFormats = []
    self.defaultText = 'Select file'
    
    self.formats = allFormatsDict.keys()
    self.formats.sort()
     
    BasePopup.__init__(self, parent=parent, title= "Project '%s': " % project.name +  'Write ccpNmr mapping file', modal=False, transient=True)

  def body(self, master):
      
    row = 0
    
    label = Label(master, text= "CcpNmr resonance -> atom mapping file writer.")
    label.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.W)
    
    row += 1
    
    label = Label(master, text= "File formats:")
    label.grid(row=row, column=0, sticky=Tkinter.W)

    self.formatListBox = ScrolledListbox(master, width = 50, height = 5, selectmode = Tkinter.MULTIPLE,
                                       initial_list = self.formats)
    self.formatListBox.grid(row=row, column=1, sticky=Tkinter.EW)
    
    # Default is all formats selected...
    self.formatListBox.setSelectedItems(self.formatListBox.getItems())

    row += 1
    
    label = Label(master, text= "Mapping output file:")
    label.grid(row=row, column=0, sticky=Tkinter.W)

    self.fileButton = Tkinter.Button(master, text = self.defaultText, command = self.selectFile)
    self.fileButton.grid(row=row, column=1, sticky=Tkinter.W)

    row += 1

    texts = [ 'Write mapping' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, columnspan = 2, column=0)

  def selectFile(self):
   
    fileName = self.fileButton.__getitem__('text')
    
    if fileName == self.defaultText:
      fileName = 'ccpNmr.map'
    
    popup = FormatFilePopup(self, file = fileName, component = 'mapping', format = ccpNmr_kw)
    
    if popup.fileSelected:
    
      self.fileButton.config(text = popup.file)      
      popup.destroy()

  def apply(self):
    
    self.selectedFormats = self.formatListBox.getSelectedItems()
    
    fileName = self.fileButton.__getitem__('text')
    
    if fileName == self.defaultText:
      return False
    
    fileCreated = writeMappingFile(self.project,fileName,originalFormats = self.selectedFormats)
    
    if fileCreated:
      showInfo("Success","Succesfully wrote mapping file")
    else:
      showError("Not written","Error writing file %s. File not written" % fileName)
      return False
          
    return True


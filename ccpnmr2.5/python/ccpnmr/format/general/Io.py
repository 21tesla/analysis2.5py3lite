
"""
======================COPYRIGHT/LICENSE START==========================

Io.py: Defines file interactions for ccpnmr/format and subdirectories

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

import os, string

from memops.api import Implementation

from memops.universal.Io import getTopDirectory
from memops.universal.Io import getPythonDirectory
from memops.universal.Io import joinPath

class TextPipe:

  def __init__(self,textArea, guiWindow = None):
    
    self.textArea = textArea
    self.guiWindow = guiWindow
    
    if self.guiWindow:
      self.textArea.tag_config("warning", foreground="orange")
      self.textArea.tag_config("error", foreground="red")
    
  def write(self,text):
  
    if self.guiWindow:
    
      if self.guiWindow.state() == 'withdrawn':
        self.guiWindow.open()
      else:
        self.guiWindow.focus()
      
      endIndex = self.textArea.index('end')

      if text.find('Warning') != -1 or text.find('WARNING') != -1:
        self.textArea.insert(endIndex,text,('warning'))
      elif text.find('Error') != -1 or text.find('ERROR') != -1:
        self.textArea.insert(endIndex,text,('error'))
      else:
        self.textArea.append(text)
      
    else:
      self.textArea.append(text)
      
    if self.guiWindow:
      import Tkinter
      self.textArea.see(Tkinter.END)

def getHelpUrlDir():

  pythonDir = getPythonDirectory()  

  return 'file:' + joinPath(pythonDir,'ccpnmr','format','doc')

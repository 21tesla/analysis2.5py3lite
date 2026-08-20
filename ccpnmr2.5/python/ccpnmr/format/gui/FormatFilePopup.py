
"""
======================COPYRIGHT/LICENSE START==========================

FormatFilePopup.py: GUI window for selecting file in a particular format

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

import os, string

from memops.universal.Constants import dirsep
from memops.universal.Io import splitPath, normalisePath

from memops.gui.FileSelect import FileSelect, FileType
from memops.gui.MessageReporter import showError
from memops.gui.Util import createDismissHelpButtonList

from memops.gui.BasePopup import BasePopup

from ccpnmr.format.general.Constants import fileTypeDict

class FormatFilePopup(BasePopup):

  def __init__(self, parent, file = None, files = None, component = None, selectionText = 'Select', title = None, format = None, multiSelect = False):
    
    #
    # Also get format (from self?)
    #
    
    self.file = None # Use this to return single file info, for backward compat
    
    if not multiSelect and file:
      self.files = [file]
      
    elif files:
      self.files = files
      
    else:
      self.files = []
      
    if self.files:
      for i in range(len(self.files)):
        self.files[i] = normalisePath(self.files[i])

    self.component = component
    self.fileSelected = False
    self.selectionText = selectionText
    self.format = format
    
    self.multiSelect = multiSelect
    
    if not title:
      title = 'Select %s file' % component

    BasePopup.__init__(self, parent=parent, title = title, modal = True, transient=True)

  def body(self, master):

    master.grid_columnconfigure(0, weight = 1)
    master.grid_rowconfigure(0, weight = 1)
    master.grid_rowconfigure(1, weight = 1)
  
    self.geometry('500x500')

    fileNames = []
    
    # All files have to be in same directory!
    if self.files and self.files[0].find(dirsep) > -1:
      (path,fileName) = splitPath(self.files[0])
      directory = path
      fileNames.append(fileName)
      
      for file in self.files[1:]:
        (path,fileName) = splitPath(file)
        fileNames.append(fileName)
    
    else:
      path = None
      directory = normalisePath(os.path.abspath(os.curdir))
      fileNames = []
      
    #print directory
    #print fileName
    
    file_types = []
    
    if fileTypeDict.has_key(self.component):
    
      if self.format and fileTypeDict[self.component].has_key(self.format):
      
        file_types.extend([FileType(self.format, fileTypeDict[self.component][self.format])])

      if fileTypeDict[self.component].has_key('generic'):
      
        file_types.extend([FileType('generic', fileTypeDict[self.component]['generic'])])
        
    file_types.extend([FileType('All', ['*'])])
      
    self.file_select = FileSelect(master,
                                  file_types=file_types,
	                                directory = directory,
                                  double_callback=self.ok,
                                  multiSelect = self.multiSelect)
    
    if self.multiSelect:
      self.file_select.setFiles(fileNames)
    else:
      if fileNames:
        fileName = fileNames[0]
      else:
        fileName = ''
      self.file_select.setFile(fileName)
	              
    self.file_select.grid(row=0, column=0, sticky=Tkinter.NSEW)

    texts = [ self.selectionText ]
    commands = [ self.ok ]
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands)
    buttons.grid(row=1, column=0, sticky=Tkinter.EW)

  def apply(self):
  
    if not self.multiSelect:
      file = self.file_select.getFile()
      if file:
        files = [file]
      else:
        files = []
        
    else:
      files = self.file_select.getFiles()
    
    if (not files):
      self.files = None
      return False
          
    self.fileSelected = True
    self.files = []
    
    # Handle symbolic links in Linux/unix
    for file in files:
      symbolicLinkIndex = string.find(file,'->')
      if symbolicLinkIndex != -1:
        self.files.append(file[symbolicLinkIndex + 3:])
      else:
        self.files.append(file)
    
    # For single files, set self.file for backward compat
    if not self.multiSelect:
      self.file = self.files[0]
    
    return True

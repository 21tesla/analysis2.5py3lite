
"""
======================COPYRIGHT/LICENSE START==========================

FileSelectPopup.py: <write function here>

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
import os

import Tkinter


from memops.gui.BasePopup import BasePopup
from memops.gui.FileSelect import FileSelect
from memops.gui.MessageReporter import showError
from memops.gui.Util import createDismissHelpButtonList

class FileSelectPopup(BasePopup):

  def __init__(self, parent, file_types = None, directory = None, multiSelect=False,
               title = 'Browse files', prompt = None, show_file = True,
               file = '', dismiss_text = 'Close', extra_dismiss_text = '',
               selected_file_must_exist = False, default_dir = None, *args, **kw):
 
    self.file_types = file_types
    self.directory = directory
    self.prompt = prompt
    self.show_file = show_file
    self.dismiss_text = dismiss_text
    self.extra_dismiss_text = extra_dismiss_text
    self.initial_file = file
    self.selected_file_must_exist = selected_file_must_exist
    self.default_dir = default_dir
    self.multiSelect = multiSelect

    kw['title'] = title
    kw['transient'] = True
    kw['modal'] = True
    BasePopup.__init__(self, parent=parent, *args, **kw)

  def body(self, master):

    self.geometry('600x400')
    master.grid_rowconfigure(0, weight=1)
    master.grid_columnconfigure(0, weight=1)

    self.result = ''
    self.file_select = FileSelect(master, file_types=self.file_types,
                                  directory=self.directory, prompt=self.prompt,
                                  show_file=self.show_file, file=self.initial_file,
                                  multiSelect=self.multiSelect,
                                  default_dir = self.default_dir)
    self.file_select.grid(row=0, column=0, sticky=Tkinter.NSEW)

    texts = [ 'Ok' ]
    commands = [ self.ok ]
    if (self.extra_dismiss_text):
      texts.append(self.extra_dismiss_text)
      commands.append(self.extra)
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands,
                                          dismiss_text=self.dismiss_text)
    buttons.grid(row=1, column=0, sticky=Tkinter.EW)

  def apply(self):

    self.result = ''

    file = self.file_select.getFile()
    if (self.selected_file_must_exist and not os.path.exists(file)):
      showError('No file', 'File does not exist', self)
      return False

    self.result = file

    return True

  def extra(self):

    self.result = None
    self.close()

  def getFile(self):

    if (self.result):
      return self.file_select.getFile()
    else:
      return self.result

  def getDirectory(self):

    if (self.result):
      return self.file_select.getDirectory()
    else:
      return self.result

  def open(self, file = ''):

    if (file):
      self.setFile(file)
    BasePopup.open(self)

  def __getattr__(self, name):
 
    try:
      return getattr(self.__dict__['file_select'], name)
    except:
      raise AttributeError, "%s instance has no attribute '%s'" % (self.__class__.__name__, name)

if __name__ == '__main__':

  from FileSelect import FileType

  root = Tkinter.Tk()

  type1 = FileType("All", ["*"])
  type2 = FileType("Text", ["*.txt"])
  type3 = FileType("Python Source", ["*.py"])

  file_types = [ type1, type2, type3 ]

  popup = FileSelectPopup(root, file_types, dismiss_text='Cancel',
                          extra_dismiss_text='Skip')
  print 'file = "' +  str(popup.getFile()) + '"'
  popup.destroy()

  popup = FileSelectPopup(root, show_file=False)
  print 'file = "' +  popup.getFile() + '"'
  popup.destroy()

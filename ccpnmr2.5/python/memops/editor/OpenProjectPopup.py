
"""
======================COPYRIGHT/LICENSE START==========================

OpenProjectPopup.py: <write function here>

Copyright (C) 2007 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

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
import time
import Tkinter

from memops.general.Implementation import ApiError
from memops.general.Io import loadProject
from memops.format.xml.Util import doesRepositoryContainProject, getPossibleProjectFiles

from memops.universal.Io import joinPath, printExceptionAndTrace

from memops.gui.DataEntry import askDir, askFile
from memops.gui.FileSelect import FileSelect, FileType
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError, showWarning
from memops.gui.Util import createDismissHelpButtonList

from memops.editor.BasePopup import BasePopup

class OpenProjectPopup(BasePopup):

  def __init__(self, parent, title = 'Project : Open', callback = None,
               help_msg = '', help_url = '', load_project = None, *args, **kw):

    self.callback = callback
    self.help_msg = help_msg
    self.help_url = help_url
    self.project = None

    if not load_project:

      askdir = lambda title, prompt, initial_value: askDir(title, prompt,
                      initial_value, parent=self, extra_dismiss_text='Skip')
      askfile = lambda title, prompt, initial_value: askFile(title, prompt,
                         initial_value, parent=self, extra_dismiss_text='Skip')
      load_project = lambda path: loadProject(path, showWarning=showWarning,
                                              askDir=askdir, askFile=askfile)

    self.load_project = load_project

    BasePopup.__init__(self, parent=parent, title=title, *args, **kw)

  def body(self, master):

    self.geometry('600x500')
    master.grid_rowconfigure(1, weight=1)
    master.grid_columnconfigure(0, weight=1)

    label = Label(master, text='Select Project Directory')
    label.grid(row=0, column=0, sticky=Tkinter.W)

    file_types = [ FileType('Project', ['*.xml']), FileType('All', ['*']) ]
    self.file_select = FileSelect(master, file_types=file_types, show_file=False,
                                  double_callback=self.ok, getRowColor=self.getEntryColor,
                                  getExtraCell=self.getProjectFileText,
                                  extraHeadings=('Status',), extraJustifies=('left',))
    self.file_select.grid(row=1, column=0, sticky=Tkinter.NSEW)

    texts = [ 'Open' ]
    commands = [ self.ok ]
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands,
                dismiss_text='Cancel', help_msg=self.help_msg, help_url=self.help_url)
    buttons.grid(row=2, column=0, sticky=Tkinter.EW)

    self.ok_button = buttons.buttons[0]

  def getEntryColor(self, fileName):

    # this function returns background color for entry

    if doesRepositoryContainProject(fileName):
      color = '#C0FFC0'  # green
    elif getPossibleProjectFiles(fileName):
      color = '#FFFFC0'  # yellow
    else:
      color = None  # default

    return color

  def getProjectFileText(self, fileName):

    if doesRepositoryContainProject(fileName):
      tt = time.localtime(os.path.getmtime(fileName))
      tt = '%4d-%02d-%02d %02d:%02d:%02d' % tt[:6]
      data = 'OK (%s)' % tt
    else:
      files =  getPossibleProjectFiles(fileName)
      nfiles = len(files)
      if nfiles == 1:
        data = '%s' % os.path.basename(files[0])
      elif nfiles:
        data = '%d possible' % nfiles
      else:
        data = None
    data = (data,)

    return data

  def disableOk(self):

    self.file_select.double_callback = None
    self.ok_button.config(state=Tkinter.DISABLED)

  def enableOk(self):

    self.file_select.double_callback = self.ok
    self.ok_button.config(state=Tkinter.NORMAL)

  def apply(self):

    # protect against repeated double-clicking

    if (self.project):
      return False

    path = self.file_select.getDirectory()
    if not doesRepositoryContainProject(path) and not getPossibleProjectFiles(path):
      self.file_select.changeDir(path)
      return False

    self.disableOk()

    self.project = None

    if not path:
      showError('No path', 'Need to select path', self)
      self.enableOk()
      return False

    try:
      self.project = self.load_project(path)
    except ApiError, e:
      showError('Loading project', e.error_msg, self)
      self.enableOk()
      return False

    except IOError, e:
      showError('Loading project', str(e), self)
      self.enableOk()
      return False

    except Exception, e:
      showError('Loading project', str(e), self)
      self.enableOk()
      return False

    except:
      showError('Loading project', 'Unknown exception', self)
      self.enableOk()
      return False

    if self.project:
      try:
        if self.callback:
          self.callback(self.project)
        ok = True
      except:
        ok = False
        printExceptionAndTrace()
      self.enableOk()
      return ok
    else:
      self.enableOk()
      return False

if (__name__ == '__main__'):

  r = Tkinter.Tk()

  p = OpenProjectPopup(r, transient=True, modal=True)
  if (p.project):
    print p.project.name
  p.destroy()

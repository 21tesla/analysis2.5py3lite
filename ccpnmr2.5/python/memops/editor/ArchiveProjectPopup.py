
"""
======================COPYRIGHT/LICENSE START==========================

ArchiveProjectPopup.py: <write function here>

Copyright (C) 2010 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

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
- MSD website (http://www.ebi.ac.uk/msd/)

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
import datetime
import os
import Tkinter

from memops.universal.Io import joinPath

from memops.general.Io import getUserDataPath, packageProject

from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FileSelect import FileSelect
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showYesNo
from memops.gui.Util import createDismissHelpButtonList

from memops.editor.BasePopup import BasePopup

class ArchiveProjectPopup(BasePopup):

  """
  **Archive the Project**

  This popup window enables the user to archive the CCPN Project into a gzipped
  tar file.  By default it only includes the *.xml files which reside in the
  project directory.  If the "Include *.xml.bak files" check button is checked
  then it will also include the *.xml.bak files which reside in the project
  directory.  If the "Include data files which are in the project directory"
  check button is checked then it will also include the binary data files
  which reside in the project directory.

  **Caveats & Tips**

  The archive excludes the reference data *.xml files (quite sensibly).
"""

  def __init__(self, parent, project, title = 'Project : Archive', callback = None,
               help_msg = '', help_url = '', dismiss_text = '', *args, **kw):

    self.callback = callback
    self.help_msg = help_msg
    self.help_url = help_url
    self.dismiss_text = dismiss_text
    BasePopup.__init__(self, parent=parent, project=project, title=title, *args, **kw)

  def body(self, guiParent):

    now = datetime.date.today().strftime('%y%m%d')
    filePrefix = '%s_%s' % (self.project.name, now)
    projDir = getUserDataPath(self.project)
    directory = os.path.dirname(projDir)

    guiParent.grid_rowconfigure(0, weight=1)
    guiParent.grid_columnconfigure(1, weight=1)

    row = 0
    label = Label(guiParent, text='Archive File:')
    label.grid(row=row, column=0, sticky=Tkinter.E)
    tipText = 'File name (excluding .tgz ending) for archive'
    self.fileEntry = Entry(guiParent, text=filePrefix, tipText=tipText)
    self.fileEntry.grid(row=row, column=1, sticky=Tkinter.EW)
    label = Label(guiParent, text='.tgz (automatically appended)')
    label.grid(row=row, column=2, sticky=Tkinter.W)

    row = row + 1
    self.backupCheck = CheckButton(guiParent, text='Include *.xml.bak files',
                         tipText='If checked include *.xml.bak files')
    self.backupCheck.grid(row=row, column=1, columnspan=2, sticky=Tkinter.W)

    row = row + 1
    self.dataCheck = CheckButton(guiParent, text='Include data files which are in project directory',
                         tipText='If checked include data files if they are located in project directory')
    self.dataCheck.grid(row=row, column=1, columnspan=2, sticky=Tkinter.W)

    row = row + 1
    labelFrame = LabelFrame(guiParent, text='Archive Location')
    labelFrame.grid(row=row, column=0, columnspan=3, sticky=Tkinter.NSEW)
    labelFrame.grid_rowconfigure(0, weight=1)
    labelFrame.grid_columnconfigure(0, weight=1)

    self.dirSelect = FileSelect(labelFrame, directory=directory,
                                show_file=False)
    self.dirSelect.grid(row=0, column=0, sticky=Tkinter.NSEW)

    guiParent.grid_rowconfigure(row, weight=1)

    row = row + 1
    texts = [ 'Save' ]
    tipTexts = [ 'Create archive file' ]
    commands = [ self.save ]
    buttons = createDismissHelpButtonList(guiParent, texts=texts,
                tipTexts=tipTexts, commands=commands,
                help_msg=self.help_msg, help_url=self.help_url,
                dismiss_text=self.dismiss_text, expands=True)
    buttons.grid(row=row, column=0, columnspan=3, sticky=Tkinter.EW)

  def save(self):

    filePrefix = self.fileEntry.get()
    directory = self.dirSelect.getDirectory()
    filePrefix = joinPath(directory, filePrefix)
    fileName = filePrefix + '.tgz'

    if os.path.exists(fileName):
      if not showYesNo('File exists', 'File "%s" exists, overwrite?' % fileName, parent=self):
        return

    includeBackups = self.backupCheck.isSelected()
    includeData = self.dataCheck.isSelected()
    packageProject(self.project, filePrefix, includeBackups, includeData)


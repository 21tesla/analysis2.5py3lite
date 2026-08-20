
"""
======================COPYRIGHT/LICENSE START==========================

SaveProjectFrame.py: <write function here>

Copyright (C) 2008 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

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

from memops.universal.Io import joinPath

from memops.general import Implementation, Io

from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FileSelect import FileSelect
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showError, showInfo, showYesNo, showWarning, showOkCancel
from memops.gui.Util import createDismissHelpButtonList

from memops.format.xml.Util import doesRepositoryContainProject, getPossibleProjectFiles

class SaveProjectFrame(Frame):

  def __init__(self, guiParent, project, callback = None,
               help_msg='', help_url='', dismiss_text='',
               modal=False, *args, **kw):

    self.project = project
    self.callback = callback
    self.help_msg = help_msg
    self.help_url = help_url
    self.dismiss_text = dismiss_text
    self.modal = modal

    self.did_save = False

    Frame.__init__(self, guiParent, *args, **kw)

    projDir = Io.getUserDataPath(self.project)

    guiParent.grid_columnconfigure(1, weight=1)

    row = 0

    label = Label(guiParent, text='Project Name:') 
    label.grid(row=row, column=0, sticky=Tkinter.E)

    self.proj_name_entry = Entry(guiParent, text=self.project.name, returnCallback=self.updateInfo,
                                 leaveCallback=self.updateInfo,
                                 tipText='The name used for the project save directory')
    self.proj_name_entry.grid(row=row, column=1, sticky=Tkinter.EW)

    row = row + 1
    label = Label(guiParent, text='Project Directory:')
    label.grid(row=row, column=0, sticky=Tkinter.E)

    label = self.proj_dir_label = Label(guiParent, text=projDir)
    label.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1
    label = Label(guiParent, text='Note: Project Directory = Save Location + Project Name')
    label.grid(row=row, column=1, sticky=Tkinter.W)

    text = 'Save binary data with project'
    tipText = 'Copy data files (e.g. for spectra) into new project directory if not already in current project directory: careful, this can take some time'
    row = row + 1
    self.dataCheckButton = CheckButton(guiParent, text=text, tipText=tipText)
    self.dataCheckButton.grid(row=row, column=0, columnspan=2, sticky=Tkinter.W)

    row = row + 1
    guiParent.grid_rowconfigure(row, weight=1)
    labelFrame = LabelFrame(guiParent, text='Save Location')
    labelFrame.grid(row=row, column=0, columnspan=2, sticky=Tkinter.NSEW)
    labelFrame.grid_rowconfigure(0, weight=1)
    labelFrame.grid_columnconfigure(0, weight=1)

    directory = os.path.dirname(projDir)
    self.proj_dir_select = FileSelect(labelFrame, directory=directory,
                             select_dir_callback=self.selectDir,
                             change_dir_callback=self.updateInfo,
                             should_change_dir_callback=self.shouldChangeDir,
                             getRowColor=self.getEntryColor,
                             show_file=False)
    self.proj_dir_select.grid(row=0, column=0, sticky=Tkinter.NSEW)

    row = row + 1
    texts = [ 'Save' ]
    tipTexts = [ 'Save project with specified name in specified directory' ]
    commands = [ self.save ]
    buttons = createDismissHelpButtonList(guiParent, texts=texts,
                tipTexts=tipTexts, commands=commands,
                help_msg=self.help_msg, help_url=self.help_url,
                dismiss_text=self.dismiss_text, expands=True)
    buttons.grid(row=row, column=0, columnspan=2, sticky=Tkinter.EW)

  def save(self):

    projName = self.proj_name_entry.get()
    directory = self.proj_dir_select.getDirectory()
    directory = joinPath(directory, projName)
    if self.isProjectDirectory(directory):
      if not showOkCancel('Overwrite directory', 'Overwrite existing project directory?', parent=self):
        return

    self.updateInfo()

    self.did_save = False
    changeDataLocations = self.dataCheckButton.isSelected()
    done = False
    try:
      done = Io.saveProject(self.project, newPath=directory, newProjectName=projName,
                            createFallback=True, showYesNo=showYesNo,
                            changeDataLocations=changeDataLocations,
                            showWarning=showWarning)

      if done:
        showInfo('Project saved', 'Project saved successfully')
        self.did_save = True
        if self.callback:
          self.callback(self.project)
      elif self.modal:
        return # give another chance
    except Implementation.ApiError, e:
      showError('Save project', e.error_msg)
    except IOError, e:
      showError('Save project', str(e))

    if done:
      self.close()

  def isProjectDirectory(self, directory):

    return (doesRepositoryContainProject(directory) or getPossibleProjectFiles(directory))

  def selectDir(self, directory):

    set_dir = True
    if self.isProjectDirectory(directory):
      baseName = os.path.basename(directory)
      self.proj_name_entry.set(baseName)
      set_dir = False

    self.updateInfo()
    return set_dir

  def shouldChangeDir(self, directory):

    shouldNotChange = self.isProjectDirectory(directory)
    if shouldNotChange:
      showWarning('Warning', 'Cannot save project inside a project directory', parent=self)

    return not shouldNotChange

  def updateInfo(self, *extra):

    projName = self.proj_name_entry.get()
    directory = self.proj_dir_select.getDirectory()
    directory = joinPath(directory, projName)

    self.proj_dir_label.set(directory)

  def getEntryColor(self, fileName):

    # this function returns background color for entry

    if self.isProjectDirectory(fileName):
      color = '#F08080'  # pink
    else:
      color = None  # default

    return color


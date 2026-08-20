
"""
======================COPYRIGHT/LICENSE START==========================

BackupProject.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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

from memops.api.Implementation import Url

from memops.universal.Io import joinPath, normalisePath

from memops.gui.Button import Button
from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Entry import Entry
from memops.gui.IntEntry import IntEntry
from memops.gui.Label import Label
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.MessageReporter import showYesNo, showWarning
from memops.gui.RadioButtons import RadioButtons

from ccpnmr.analysis.popups.BasePopup import BasePopup

# for this to work parent needs defined setBackupOn(), setBackupOff() and doBackup()

class BackupProjectPopup(BasePopup):

  """
  **Create Automatic Project Backup**

  The purpose of this dialog is to allow the user to create backups
  of their project automatically.  The XML files for the backup go
  into a separate directory from the project itself.  If the project
  directory is called PROJECTDIR then the default for the backup
  directory is PROJECTDIR_backup, although that can be changed to
  something else in this dialog.

  Other than the backup directory, the user can specify the frequency
  of the backup, in minutes.

  The "Apply Auto Settings" does not have to be applied if the user
  has changed the on/off setting or if the user has changed the
  Auto-backup frequency and entered a carriage return.

  The "Do Immediate Backup" is in case the user wants to do a backup
  then and there.

  There is no backup for the backup.
  """

  on_off_entries = ['on', 'off']

  def __init__(self, parent, help_msg = '', help_url = '', *args, **kw):

    self.help_msg = help_msg
    self.help_url = help_url

    BasePopup.__init__(self, parent=parent, title='Project : Backup', **kw)

  def body(self, guiFrame):

    self.geometry('500x130')

    guiFrame.grid_columnconfigure(2,weight=1)
    guiFrame.grid_rowconfigure(3,weight=1)

    self.alarm_id = None

    row = 0
    tipText = 'Browse for the directory into which project backups will be made'
    button = Button(guiFrame, text='Directory:', command=self.selectDir,
                    grid=(row, 0), tipText=tipText, sticky='ew')

    repository = self.project.findFirstRepository(name='backup')
    if repository:
      text = repository.url.path
    else: # this is trouble, but should not happen
      text = ''
    
    tipText = 'Enter the name of the directory into which project backups will be made' 
    self.dir_entry = Entry(guiFrame, text=text, tipText=tipText,
                           width=40, returnCallback=self.applyAuto,
			   grid=(row, 1), gridSpan=(1,2), sticky='ew')
			   
    row += 1
    label = Label(guiFrame, text='Auto-backup:', grid=(row, 0))
    if self.analysisProject.doAutoBackup:
      ind = 0
    else:
      ind = 1
    
    tipTexts = ['Toggle the timed automatic backup on',
                'Toggle the timed automatic backup off']  
    self.on_off_buttons = RadioButtons(guiFrame, entries=self.on_off_entries,
                                       tipTexts=tipTexts, select_callback=self.applyAuto,
				        selected_index=ind, grid=(row,1))
 
    row +=1
    tipText = 'The number of minutes to wait before automatic project backups'
    label = Label(guiFrame, text='Auto-backup frequency:',
                  grid=(row,0))
                   
    self.freq_entry = IntEntry(guiFrame, text=self.analysisProject.autoBackupFreq,
                               returnCallback=self.applyAuto, tipText=tipText)
    self.freq_entry.grid(row=row, column=1)
    
    label = Label(guiFrame, text=' (minutes)', grid=(row,2))
 
    row = row + 1
    # Blank for expansion

    row = row + 1
    texts = [ 'Apply Auto Settings', 'Do Immediate Backup' ]
    commands = [ self.applyAuto, self.applyManual ]
    tipTexts = ['Commit the specified settings and commence automated CCPN project backup',
                'Backup the CCPN project now, into the specified backup directory']
    buttons = UtilityButtonList(guiFrame, texts=texts, commands=commands,
                                doClone=False, helpMsg=self.help_msg,
                                helpUrl=self.help_url, tipTexts=tipTexts,
                                grid=(row, 0), gridSpan=(1,3), sticky='nsew')

  def selectDir(self):

    popup = FileSelectPopup(self, show_file=False)
    dir = popup.getDirectory()
    popup.destroy()

    if (dir):
      self.dir_entry.set(dir)
      self.setDir()

  def setDir(self):

    directory = self.dir_entry.get()
    if not directory:
      showWarning('No directory', 'No directory specified, not setting backup directory', parent=self)
      return

    if not os.path.abspath(directory):
      directory = joinPath(os.getcwd(), directory)

    if os.path.exists(directory):
      if not os.path.isdir(directory):
        showWarning('Path not directory', 'Path "%s" exists but is not a directory' % directory, parent=self)
        return
    else:
      if showYesNo('Directory does not exist', 'Directory "%s" does not exist, should it be created?' % directory, parent=self):
        os.mkdir(directory)
      else:
        showWarning('Directory not created', 'Directory "%s" not created so not setting backup directory' % directory, parent=self)
        return

    repository = self.project.findFirstRepository(name='backup')
    if not repository:
      showWarning('No backup repository', 'No backup repository found (something wrong somewhere)', parent=self)
      return

    url = repository.url
    if url.path != directory:
      repository.url = Url(path=normalisePath(directory))

  def setFreq(self):

    freq = self.freq_entry.get()

    if (freq is None):
      freq = 0

    self.analysisProject.autoBackupFreq = freq

  def setInfo(self):

    self.setFreq()
    self.setDir()

    on_off = (self.on_off_buttons.get() == 'on')
    self.analysisProject.doAutoBackup = on_off

  def applyAuto(self, *event):

    self.setInfo()
    on_off = self.on_off_buttons.get()
    if on_off == 'on':
      self.parent.setBackupOn()
    else:
      self.parent.setBackupOff()

  def applyManual(self):

    self.setDir()
    self.parent.doBackup()

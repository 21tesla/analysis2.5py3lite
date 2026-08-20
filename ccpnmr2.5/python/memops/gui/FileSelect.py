
"""
======================COPYRIGHT/LICENSE START==========================

FileSelect.py: <write function here>

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
import Tkinter

import glob
import os
import sys

from datetime import datetime

from memops.universal.Constants import dirsep
from memops.universal.Io import normalisePath, splitPath, joinPath
from memops.universal.Util import isWindowsOS

from memops.gui.ButtonList import ButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.LabeledEntry import LabeledEntry
from memops.gui.MessageReporter import showError
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix

from memops.universal.Io import getTopDirectory

ICON_NAMES = ['go-previous','go-next',
              'go-up', 'go-top', 'go-home',
              'folder-new','view-refresh']
GFX_DIR = os.path.join(getTopDirectory(),'python','memops','gui','graphics','22x22')

defaultFilter = '*'

class FileType:

  def __init__(self, message, filters=None):
    
    if not filters:
      filters = [defaultFilter,]
    
    self.message = message
    self.filters = filters

class FileSelect(Frame):

  def __init__(self, parent, file_types = None, directory = None,
               single_callback = None, double_callback = None,
               select_dir_callback = None, change_dir_callback = None,
               should_change_dir_callback = None,
               prompt = None, show_file = True, file = '', multiSelect=False,
               default_dir = None, getRowColor = None,
               getExtraCell = None, extraHeadings = None, extraJustifies = None,
               displayExtra = True, manualFileFilter=False, extraTipTexts=None, 
               *args, **kw):
    
    if file_types is None:
      file_types = [ FileType("All", ["*"]) ]
    
    if manualFileFilter:
      self.manualFilter = FileType("Manual")
      file_types.append(self.manualFilter)
    else:
      self.manualFilter = None
    
    if directory is None:
      directory = normalisePath(os.getcwd())

    if extraHeadings is None:
      extraHeadings = ()
    else:
      extraHeadings = tuple(extraHeadings)

    if extraJustifies is None:
      extraJustifies = ()
    else:
      extraJustifies = tuple(extraJustifies)

    if extraHeadings or extraJustifies:
      assert len(extraHeadings) == len(extraJustifies)
      assert getExtraCell

    self.extraHeadings = extraHeadings
    self.extraJustifies = extraJustifies
    self.extraTipTexts = extraTipTexts or []
    self.displayExtra = displayExtra

    Frame.__init__(self, parent, *args, **kw)

    self.file_types = file_types
    self.fileType = file_types[0]
    self.single_callback = single_callback
    self.double_callback = double_callback
    self.select_dir_callback = select_dir_callback
    self.change_dir_callback = change_dir_callback
    self.should_change_dir_callback = should_change_dir_callback
    self.show_file = show_file
    self.directory = None
    self.historyBack = []
    self.historyFwd = []
    self.determineDir(directory)
    self.default_dir = default_dir
    self.getRowColor = getRowColor
    self.getExtraCell = getExtraCell

    
    self.grid_columnconfigure(0, weight=1)

    row = 0
    if prompt:
      label = Label(self, text=prompt, grid=(row,0))
      row += 1

    self.grid_rowconfigure(row, weight=1)
    if show_file:
      headings = ('Name', 'Size', 'Date')
      justifies = ('left','right','right')
      tipTexts = ['Name of file, directory or link',
                  'Size of file in bytes',
                  'Date of last modification']
    else:
      headings = ('Directory',)
      justifies = ('left',)
      tipTexts = ['Name of directory']
    
    self.normalHeadings = headings
    self.normalJustifies = justifies
    self.normalTipTexts = tipTexts
    headings = headings + extraHeadings
    justifies = justifies + extraJustifies
    tipTexts += [None] * len(extraHeadings) 

    self.fileList = ScrolledMatrix(self, headingList=headings,
                                   justifyList=justifies,
                                   initialRows=10,
                                   callback=self.singleCallback,
                                   doubleCallback=self.doubleCallback,
                                   tipTexts=tipTexts,
                                   multiSelect=multiSelect, grid=(row,0))

    row += 1
    tipTexts = ['Go to previous location in history',
                'Go forward in location history',
                'Go up one directory level',
                'Go to root directory',
                'Go to home directory',
                'Make a new directory',
                'Refresh directory listing']
    
    texts = ['Back', 'Forward', 'Up', 'Top', 'Home', 'New', 'Refresh']
    commands = [self.backDir, self.fwdDir, self.upDir, self.topDir,
                self.homeDir, self.createDir, self.updateFileList ]

    if self.default_dir:
      texts.append('Default')
      commands.append(self.setDefaultDir)
      tipTexts.append('Set the current directory as the default')
      
    self.icons = []
    for name in ICON_NAMES:
      icon = Tkinter.PhotoImage(file=os.path.join(GFX_DIR,name+'.gif'))
      self.icons.append(icon)
      
    self.buttons = ButtonList(self, texts=texts, commands=commands,
                              images=self.icons, grid=(row,0),
                              tipTexts=tipTexts)


    if show_file:
      row += 1
      self.file_entry = LabeledEntry(self, label='File name',
                              label_width=10, entry_width=40,
                              returnCallback=self.setSelectedFile)
      self.file_entry.grid(row=row, column=0, sticky=Tkinter.EW)
    else:
      self.file_entry = None

    row += 1
    self.directory_entry = LabeledEntry(self, label='Directory',
                                        entry = directory,
                                        label_width=10, entry_width=40,
                                        returnCallback=self.entryDir)
    self.directory_entry.grid(row=row, column=0, sticky=Tkinter.EW)


    row += 1
    subFrame = Frame(self, grid=(row,0))
    subFrame.expandGrid(None,6)

    if show_file:

      label = Label(subFrame, text='File type:', grid=(0,0))

      type_labels = self.determineTypeLabels()
      self.fileType_menu = PulldownList(subFrame, callback=self.typeCallback,
                                        texts=type_labels, objects=self.file_types,
                                        grid=(0,1),
                                        tipText='Restrict listed files to selected suffix')

    label = Label(subFrame, text=' Show hidden:', grid=(0,2))

    self.hidden_checkbutton = CheckButton(subFrame, text='', selected=False,
                                          callback=self.updateFileList, grid=(0,3),
                                          tipText='Show hidden files beginning with "." etc.')

    label = Label(subFrame, text='Dir path:', grid=(0,4))

    self.pathMenu = PulldownList(subFrame, callback=self.dirCallback,
                                 objects=range(len(self.dirs)),
                                 texts=self.dirs, index=len(self.dirs)-1,
                                 indent='  ', prefix=dirsep, grid=(0,5),
                                 tipText='Directory navigation to current location')

    if show_file and self.manualFilter is not None:
      row += 1
      self.manual_filter_entry = LabeledEntry(self, label='Manual Select',
                                          entry = defaultFilter,
                                          label_width=13, entry_width=40,
                                          returnCallback=self.entryFilter,
                                          tipText='Path specification with wildcards to select multiple files')
      self.manual_filter_entry.grid(row=row, column=0, sticky=Tkinter.EW)
    
    
    self.updateFileList()
    self.updateButtons()

    if file:
      self.setFile(file)

  def updateDisplayExtra(self, displayExtra):

    self.displayExtra = displayExtra
    self.updateFileList()

  def isRootDirectory(self, directory = ''):

    if not directory:
      directory = self.directory

    if directory:
      return os.path.dirname(directory) == directory

    return True  # arbitrary

  def updateButtons(self):
  
    buttons = self.buttons.buttons
    isRoot = self.isRootDirectory

    if self.historyBack:
      buttons[0].enable()
    else:  
      buttons[0].disable()
    
    if self.historyFwd:
      buttons[1].enable()
    else:  
      buttons[1].disable()
    
    if self.directory and not isRoot(self.directory):
      buttons[2].enable()
      buttons[3].enable()
    else:  
      buttons[2].disable()
      buttons[3].disable()
   
    homeDir = self.getHomeDir()
    if homeDir and os.path.isdir(homeDir):
      buttons[4].enable()
    else:  
      buttons[4].disable()

  def backDir(self):

    if self.historyBack:
      self.historyFwd.insert(0, self.directory)
      dirName = self.historyBack.pop()
      self.changeDir(dirName, addHistory=False)

  def fwdDir(self):

    if self.historyFwd:
      self.historyBack.append(self.directory)
      dirName = self.historyFwd.pop(0)
      self.changeDir(dirName, addHistory=False)
  
  def topDir(self):
  
    if self.directory:
      isRoot = self.isRootDirectory
      dirName = self.directory
      while not isRoot(dirName):
        dirName = joinPath(dirName, os.pardir)
  
      self.changeDir(dirName)
  
  def homeDir(self):
  
    homeDir = self.getHomeDir()
    if homeDir:
      self.changeDir(homeDir)

  def getHomeDir(self):
  
    return os.environ.get('HOME') or os.environ.get('HOMEPATH')

  def upDir(self):

    self.changeDir(joinPath(self.directory, os.pardir))

  def setDefaultDir(self):

    self.changeDir(self.default_dir)

  def createDir(self):

    from memops.gui.DataEntry import askString
    
    msg = 'Enter new sub-directory name'
    dirName = askString('New directory', msg, parent=self)
    if dirName:
      dirName = joinPath(self.directory, dirName)
      if os.path.exists(dirName):
        msg = 'Directory "%s" already exists' % dirName
        showError('Directory exists', msg, parent=self)
      else:
        os.mkdir(dirName)
        self.updateFileList()

  def determineTypeLabels(self):
    
    type_labels = []
    for t in self.file_types:
      s = t.message + ' (' + ', '.join(t.filters) + ')'
      type_labels.append(s)
    
    return type_labels

  def setFileTypes(self, file_types):
     
    self.file_types = file_types
    if self.fileType not in file_types:
      self.fileType = file_types[0]
    if self.show_file:
      ind = self.file_types.index(self.fileType)
      type_labels = self.determineTypeLabels()
      self.fileType_menu.setup(type_labels, self.file_types, ind)
      self.updateFileList()

  def determineDir(self, directory):

    directory = normalisePath(directory)

    assert os.path.isdir(directory), '"%s" is not a directory' % directory

    self.prev_directory = self.directory

    if not os.path.isabs(directory):
      if self.directory is None: # first time around
        directory = joinPath(normalisePath(os.getcwd()), directory)
      else:
        directory = joinPath(self.directory, directory)

    self.directory = directory
    self.dirs = self.directory.split(dirsep)

  def setSelectedFile(self, *event):

    file = self.file_entry.getEntry()
    try:
      self.fileList.selectObject(file)
    except:
      showError('Some err', 'some err')

  def entryDir(self, *event):

    dirName = self.directory_entry.getEntry()
    try:
      self.changeDir(dirName)
    except:
      showError('Not a directory', '"' + dirName + '" is not a directory.')

  def entryFilter(self, *event):
    
    filterString = self.manual_filter_entry.getEntry()
    filters = [str.strip(x) for x in filterString.split(',')]
    self.manualFilter.filters = filters
    self.fileType = self.manualFilter
    
    self.setFileTypes(self.file_types) # to trigger update after filter change
    self.updateFileList()
    
    self.manual_filter_entry.setEntry(filterString)

  def changeDir(self, dirName, addHistory=True):

    if not os.path.isdir(dirName):
      return

    if self.should_change_dir_callback:
      if not self.should_change_dir_callback(dirName):
        return

    oldDir = self.directory

    self.determineDir(dirName)
    self.updateFileList()
    self.directory_entry.setEntry(self.directory)
    
    if self.change_dir_callback:
      self.change_dir_callback(self.directory)

    nDirs = len(self.dirs)
    
    self.pathMenu.setup(self.dirs, range(nDirs), index=nDirs-1)
    
    if addHistory:
      if oldDir and (self.directory != oldDir):
        self.historyBack.append(oldDir)
      self.historyFwd = []
    
    self.updateButtons()

  def getEntryInfo(self, entry):

    file = joinPath(self.directory, entry)

    if os.path.islink(file):
      # plain arrow: u' \u2192 '
      entry = entry + u' \u21D2 ' + unicode(os.readlink(file), 'utf-8')
      size  = None
      color = '#E0D0C0'
      
    elif os.path.isdir(file):
      if not self.isRootDirectory(entry):
        entry = entry + dirsep
      size = None
      color = '#C0D0C0'
      
    else:
      color = '#C0C0D0'
      
      if self.show_file:
        try:
          size = str(os.path.getsize(file))
        except:
          size = None
      else:
        size = None
        
    try:
      fileTime = self.fileTime(file)
    except:
      fileTime = ''

    return (entry, size, fileTime, color)

  def getDrives(self):

    if isWindowsOS():
      #import win32api
      #drives = win32api.GetLogicalDriveStrings()
      #drives = drives.split('\x00')
      #drives = [normalisePath(drive) for drive in drives if drive]
      #drives.sort()
      # 19 Mar 2012: do not use win32api because not standard
      from ctypes import windll
      import string
      drives = []
      bitmask = windll.kernel32.GetLogicalDrives()
      for letter in string.ascii_uppercase:
        if bitmask & 1:
          drives.append(letter)
        bitmask >>= 1
      drives = [drive+':\\' for drive in drives]
      drives = [normalisePath(drive) for drive in drives]
    else:
      drives = []

    return drives

  def updateFileList(self, *extra):

    try:
      if self.directory:
        directory = self.directory
      else:
        directory = dirsep
      entries = self.getFilterFiles(directory)
    except OSError, e:
      showError('OS Error', str(e))
      if self.prev_directory:
        self.directory = None
        self.changeDir(self.prev_directory)
      return
      
    if self.directory:
      if self.isRootDirectory():
        entries[:0] = [drive for drive in self.getDrives() if drive not in self.directory]
      else:
        entries[:0] = [os.pardir] # add parent directory at front of list

    # DJOD: Not elegant, but robust !?!
    showHidden = self.hidden_checkbutton.getSelected()
    show_ent = []
    for entry in entries:
      if not showHidden:
        if not entry.startswith('.'):
          show_ent.append( entry )
      else:
        show_ent.append( entry )

    # if not show_ent.__contains__('.'):
    #   show_ent.append('.')

    if not show_ent.__contains__('..'):
      show_ent.insert(0, '..')

    entries = show_ent

    textMatrix = []
    colorMatrix = []
    for entryActual in entries:

      (entry, size, fileTime, color) = self.getEntryInfo(entryActual)
      fullEntry = joinPath(directory, entryActual)
      if self.getRowColor:
        color = self.getRowColor(fullEntry)

      if self.displayExtra and self.getExtraCell:
        data = self.getExtraCell(fullEntry)
        headingList = self.normalHeadings + self.extraHeadings
        justifyList = self.normalJustifies + self.extraJustifies
        tipTexts = self.normalTipTexts + self.extraTipTexts
      else:
        data = ()
        headingList = self.normalHeadings
        justifyList = self.normalJustifies
        tipTexts = self.normalTipTexts
        # in case we had sorted on a now removed line :
        self.fileList.lastSortLine = None 
        
      if self.show_file:
        text = [entry, size, fileTime]
      else:
        text = [entry,]
      text = text + list(data)

      textMatrix.append(text)
      numCols = len(headingList)
      colorMatrix.append(numCols*[color])

    self.fileList.update(objectList=entries, textMatrix=textMatrix,
                         colorMatrix=colorMatrix, headingList=headingList,
                         justifyList=justifyList, tipTexts=tipTexts)

  def fileTime(self, file):

    return datetime.fromtimestamp(os.path.getmtime(file))

  def getFilterFiles(self, directory):

    isDir = os.path.isdir

    filterFiles = set()
    for filtr in self.fileType.filters:
      # have to use fullfilter rather than filter because cwd is not directory in general
      fullfilter = joinPath(directory, filtr)
      filterFiles.update(glob.glob(fullfilter))

    # remove directory from file names

    n = len(directory)
    if directory[-1] != dirsep:
      n = n + 1

    filterFiles = list(set([x[n:] for x in filterFiles if not isDir(x)]))
    filterFiles.sort()


    # add in subdirectories
    subDirs = [x for x in os.listdir(directory) if isDir(joinPath(directory, x))]
    subDirs.sort()
    
    return subDirs + filterFiles

  def singleCallback(self, entry, row, col):

    file = joinPath(self.directory, entry)
    if os.path.isdir(file):
      set_dir = not self.show_file
      if self.select_dir_callback:
        set_dir = self.select_dir_callback(file)
      if set_dir:
        self.directory_entry.setEntry(joinPath(self.directory, entry))
    else:
      if self.show_file:
        self.file_entry.setEntry(entry)
      if self.single_callback:
        self.single_callback(file)

  def doubleCallback(self, entry, row, col):

    # TODO: need to worry about symbolic links leading you astray on path
    file = joinPath(self.directory, entry)
    if (self.show_file or not self.double_callback) and os.path.isdir(file):
      self.changeDir(file)
    elif self.double_callback:
      self.double_callback(file)

  def dirCallback(self, dir_index):

    if dir_index < len(self.dirs)-1:
      directory = dirsep.join(self.dirs[:dir_index+1])
      if not directory: # root directory
        directory = dirsep
      self.changeDir(directory)

  def typeCallback(self, fileType):

    if fileType != self.fileType:
      self.fileType = fileType
      self.updateFileList()

  def getFile(self, full_path = True):

    if self.show_file:
      file = self.file_entry.getEntry()
      if file and full_path:
        file = joinPath(self.directory, file)
    else: # ignore full_path
      file = self.directory_entry.getEntry()

    return file

  def getFiles(self, full_path=True):

    fileList = self.fileList.currentObjects

    if self.fileList.multiSelect and fileList:
      if full_path:
        files = [joinPath(self.directory, f) for f in fileList]
      else:
        files = fileList

    else:
      return [self.getFile(),]

    return files

  def getDirectory(self):

    return self.directory_entry.getEntry()

  def setFile(self, file):

    if self.show_file:
      (dir, file) = splitPath(file)
    else:
      dir = file

    self.setDirectory(dir)

    if self.show_file:
      self.file_entry.setEntry(file)

  def setFiles(self, fileNames):
    
    if not self.fileList.multiSelect:
      self.setFile(fileNames[0])
      return
    
    dirA = None
    
    files = []
    for fileName in fileNames:
      dirB, file = splitPath(fileName)

      if dirA is None:
        dirA = dirB

      elif dirA != dirB:
        msg = 'Cannot select files from different directories'
        showError('Set files failure', msg, parent=self)
        return

      if file in self.fileList.objectList:
        files.append(file)
    
    if dirA:
      self.setDirectory(dirA)
    
    self.fileList.selectObjects(files)

  def setDirectory(self, directory):

    if directory and directory != self.directory:
      self.changeDir(directory)

if __name__ == '__main__':

  def f(file):
    print 'double_callback', file

  def g(directory):
    print 'change_dir_callback', directory

  root = Tkinter.Tk()

  type1 = FileType("All", ["*"])
  type2 = FileType("Text", ["*.txt"])
  type3 = FileType("Python Source", ["*.py"])

  file_types = [ type1, type2, type3 ]

  file_select = FileSelect(root, file_types=file_types,
                           double_callback=f, change_dir_callback=g,
                           prompt='Select a file of your choice', multiSelect=True)
  file_select.setFiles(['Button.py','ButtonList.py','Base.py'])
  file_select.pack(side=Tkinter.TOP, expand=Tkinter.YES, fill=Tkinter.BOTH)

  root.mainloop()


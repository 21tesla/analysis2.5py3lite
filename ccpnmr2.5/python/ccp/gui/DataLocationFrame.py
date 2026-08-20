
"""
======================COPYRIGHT/LICENSE START==========================

DataLocationFrame.py: <write function here>

Copyright (C) 2008 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

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

from memops.universal.Io import joinPath, splitPath, converseSplitPath

from memops.api.Implementation import Url

from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Frame import Frame
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.MessageReporter import showOkCancel, showWarning
from memops.gui.Base import getPopup

from ccp.api.general.DataLocation import NumericMatrix
from ccp.api.general.DataLocation import MimeTypeDataStore

from ccp.general.Io import changeDataStoreUrl, setDataStoreFileName

class DataLocationFrame(Frame):

  def __init__(self, parent, project, utilityButtons=False, *args, **kw):

    self.project = project
    self.waiting  = False
    self.dataStore = None
    Frame.__init__(self, parent, *args, **kw)

    self.grid_rowconfigure(0, weight=1)
    self.grid_columnconfigure(0, weight=1)

    tipTexts = ['Row number',
                'Location of the absolute path directory, that may contain multiple spectrum files',
                'Location of spectrum file relative to its absolute path directory',
                'Identity of the spectrum as experiment:spectrum']

    headingList = ['#', 'Absolute Path', 'Relative Path\nand File', 'Spectrum / FID / Image']
    editWidgets = [None, None, None, None]
    editGetCallbacks = [None, self.getDataUrl, self.getDataPath, None]
    editSetCallbacks = [None, None, None, None]

    self.dataStoreTable = ScrolledMatrix(self, headingList=headingList, tipTexts=tipTexts,
                                          callback=self.selectDataStore,
                                          editWidgets=editWidgets, multiSelect=True,
                                          editGetCallbacks=editGetCallbacks,
                                          editSetCallbacks=editSetCallbacks,
                                          ###deleteFunc=self.deleteDataStores,
                                          grid=(0,0))

    texts = ['Propagate\nAbsolute Path',
             'Shift Directory\nTo Absolute Path',
             'Shift Directory\nTo Relative Path']
    tipTexts = ['Copy absolute path from last row selected to other selected spectra',
                'Shift the top directory from the relative path to the absolute path',
                'Shift the bottom directory from the absolute path to the relative path']
    commands = [self.propagatePath,
                self.moveDirToAbsolute,
                self.moveDirToRelative]
    
    isModal = False
    popup = getPopup(parent)
    if popup:
      isModal = popup.modal
    
    if utilityButtons:
      if isModal:
        texts = ['All Done!',] + texts
        tipTexts = ['Close the popup, satisfied that paths are OK'] + tipTexts
        commands = [popup.close,] + commands
      self.buttonList = UtilityButtonList(self, texts=texts, commands=commands,
                                          tipTexts=tipTexts, doClone=False, grid=(1,0))
      if isModal:
        self.buttonList.buttons[0].config(bg='#B0FFB0')
    else:
      self.buttonList = ButtonList(self, texts=texts, commands=commands, grid=(1,0), tipTexts=tipTexts)

    self.updateAfter()
      
  def selectDataStore(self, dataStore, row, col):
  
    self.dataStore = dataStore
  
  def moveDirToAbsolute(self):

    dataStores = self.dataStoreTable.currentObjects
    
    for dataStore in dataStores:
      absPath = dataStore.dataUrl.url.dataLocation
      relPath = dataStore.path
      
      (head, tail) = converseSplitPath(relPath)

      if head:
        if head[0] in ('/', '.'):
          head = head[1:]
        relPathNew = tail
        absPathNew = os.path.join(absPath, head)
 
        if relPathNew != relPath:
          # No notifiers for temporarily invalid path
          ###dataUrl = dataStore.dataUrl
          self.project.__dict__['override'] = True
          # 22 Aug 2014: the below completely messes up the other dataStores
          ###for dataUrlStore in dataUrl.dataStores:
          ###  dataUrlStore.path = converseSplitPath(dataUrlStore.path)[1]
          dataStore.path = relPathNew
          self.project.__dict__['override'] = False
          changeDataStoreUrl(dataStore, absPathNew)
      
  def moveDirToRelative(self):

    dataStores = self.dataStoreTable.currentObjects
    for dataStore in dataStores:
      absPath = dataStore.dataUrl.url.dataLocation
      relPath = dataStore.path
      
      absPaths = splitPath(absPath)

      if absPaths[1]:
        absPathNew = absPaths[0]
        relPathNew = joinPath(absPaths[1], relPath)
 
        if relPathNew != relPath:
          # No notifiers for temporarily invalid path
          ###dataUrl = dataStore.dataUrl
          self.project.__dict__['override'] = True
          # 22 Aug 2014: the below completely messes up the other dataStores
          # everything is handled in changeDataStoreUrl so only need to change dataStore
          ###for dataUrlStore in dataUrl.dataStores:
          ###  dataUrlStore.path = joinPath(absPaths[1], dataUrlStore.path)
          dataStore.path = relPathNew
          self.project.__dict__['override'] = False
          changeDataStoreUrl(dataStore, absPathNew)
  
  def propagatePath(self):
    
    dataStores = self.dataStoreTable.currentObjects
    
    if self.dataStore and len(dataStores) > 1:
      sourcePath = self.dataStore.dataUrl.url.dataLocation
      if os.path.exists(sourcePath):
        directory = sourcePath
      elif os.path.exists(os.path.dirname(sourcePath)):
        directory = os.path.dirname(sourcePath)
      else:
        directory = None
        
      if not directory:
        msg = 'Source directory does not exist'
        showWarning('Failure', msg, parent=self)
        return
      
      msg = 'Transfer path "%s" to selection?' % directory
      if showOkCancel('Confirm', msg, parent=self):
        for dataStore in dataStores:
          if dataStore is self.dataStore:
            continue
          
          if dataStore.dataUrl.url.dataLocation == directory:
            continue
          
          changeDataStoreUrl(dataStore, sourcePath)

  def updateAfter(self, *extra):
    
    if self.waiting:
      return
    
    self.waiting = True
    self.after_idle(self.update)

  def update(self):

    dataStores = []
    for dataLocationStore in self.project.sortedDataLocationStores():
      dataStores.extend(dataLocationStore.sortedDataStores())

    textMatrix = []
    colorMatrix = []
    for i, dataStore in enumerate(dataStores):
      if isinstance(dataStore, NumericMatrix):
        tt = ['%s:%s' % (dataSource.experiment.name, dataSource.name) for dataSource in dataStore.nmrDataSources]
      elif isinstance(dataStore, MimeTypeDataStore):
        tt = ['%s:%s:%s' % (dataSourceImage.dataSource.experiment.name, dataSourceImage.dataSource.name, dataSourceImage.serial) for dataSourceImage in dataStore.nmrDataSourceImages]
      else:
        tt = []

      n = len(tt)
      if n == 0:
        tt = ''
      elif n == 1:
        tt = tt[0]
      elif n == 2:
        pass
      else:
        tt = '[%s,%s,...]' % (tt[0], tt[1])

      datum = [i+1, dataStore.dataUrl.url.dataLocation, dataStore.path, tt]
      textMatrix.append(datum)

      colors = 4*[None]  # number of columns=4
      color = '#FF8080'  # red
      if not os.path.exists(dataStore.dataUrl.url.dataLocation):
        colors[1] = color
      elif not os.path.exists(dataStore.fullPath):
        colors[2] = color
      colorMatrix.append(colors)

    self.dataStoreTable.update(textMatrix=textMatrix, objectList=dataStores,
                               colorMatrix=colorMatrix)

    self.waiting = False

  def getDataUrl(self, *extra):

    if self.dataStore:
      dataStore = self.dataStore
      oldPath = dataStore.dataUrl.url.dataLocation
      if os.path.exists(oldPath):
        directory = oldPath
      elif os.path.exists(os.path.dirname(oldPath)):
        directory = os.path.dirname(oldPath)
      else:
        directory = None
      popup = FileSelectPopup(self, directory=directory, show_file=False)
      newPath = popup.getDirectory()
      popup.destroy()
      if newPath and oldPath != newPath:
        changeDataStoreUrl(dataStore, newPath)

  def getDataPath(self, *extra):

    if self.dataStore:
      dataStore = self.dataStore
      oldPath = dataStore.fullPath
      directory = file = ''
      if os.path.exists(oldPath):
        file = oldPath
      elif os.path.exists(os.path.dirname(oldPath)):
        directory = os.path.dirname(oldPath)
      elif os.path.exists(dataStore.dataUrl.url.dataLocation):
        directory = dataStore.dataUrl.url.dataLocation
      popup = FileSelectPopup(self, file=file, directory=directory)
      newPath = popup.getFile()
      popup.destroy()
      if newPath and oldPath != newPath:
        pp = dataStore.dataUrl.url.dataLocation
        if newPath.startswith(pp):
          newPath = newPath[len(pp)+1:]
          dataStore.path = newPath
        else:
          setDataStoreFileName(dataStore, newPath)

  def administerNotifiers(self, notifyFunc):

    for func in ('__init__', 'delete', 'setDataUrl', 'setPath'):
      for clazz in ('BlockedBinaryMatrix', 'ShapeMatrix', 'MimeTypeDataStore'):
        notifyFunc(self.updateAfter, 'ccp.general.DataLocation.%s' % clazz, func)

    notifyFunc(self.updateAfter, 'ccp.general.DataLocation.DataUrl', 'setUrl')

    for func in ('__init__', 'delete', 'setName'):
      for clazz in ('Experiment', 'DataSource', 'DataSourceImage'):
        notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.%s' % clazz, func)

if __name__ == '__main__':

  import sys

  if len(sys.argv) != 2:
    print 'must specify project directory'
    sys.exit()

  path = sys.argv[1]

  import Tkinter
  from memops.gui.DataEntry import askDir, askFile
  from memops.gui.MessageReporter import showWarning
  from memops.general.Io import loadProject
  from memops.universal.Io import normalisePath

  path = normalisePath(path)
  askdir = lambda title, prompt, initial_value: askDir(title, prompt,
              initial_value, parent=top, extra_dismiss_text='Skip')
  askfile = lambda title, prompt, initial_value: askFile(title, prompt,
              initial_value, parent=top, extra_dismiss_text='Skip')
  project = loadProject(path, showWarning=showWarning,
                        askDir=askdir, askFile=askfile)

  r = Tkinter.Tk()

  popup = DataLocationPopup(r, project)

  r.withdraw()
  r.mainloop()

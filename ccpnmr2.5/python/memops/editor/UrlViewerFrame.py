
"""
======================COPYRIGHT/LICENSE START==========================

UrlViewerFrame.py: <write function here>

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

from memops.api import Implementation as ApiImplementation
from ccp.api.nmr import Nmr

from memops.universal.Io import normalisePath

from memops.general import Implementation, Io

from memops.gui.Base import getPopup
from memops.gui.ButtonList import ButtonList
from memops.gui.DataEntry import askString
from memops.gui.Entry import Entry
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Frame import Frame
from memops.gui.MessageReporter import showYesNo, showError
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.RadioButtons import RadioButtons
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.ToggleLabel import ToggleLabel

unchangedColor = None
changedColor   = '#FF8080'

class UrlViewerFrame(Frame):

  all_used_entries = ('all', 'used')

  def __init__(self, parent, project, url_closed = True,
               show_storages = True, storage_closed = True,
               show_data_locations = True, data_location_closed = True,
               *args, **kw):

    self.project = project
    self.show_storages = show_storages
    self.show_data_locations = show_data_locations

    self.dir_popup = None

    Frame.__init__(self, parent, *args, **kw)

    popup = getPopup(parent)

    self.grid_columnconfigure(0, weight=1)
 
    row = 0
    label = ToggleLabel(self, text='Urls ', isArrowClosed=url_closed,
                        callback=self.urlToggle)
    label.grid(row=row, column=0, sticky=Tkinter.W)
    self.url_radio = RadioButtons(self, entries=self.all_used_entries,
                                  select_callback=self.updateUrlList,
                                  selected_index=1)
    row = row + 1

    self.url_frame = Frame(self)
    self.url_row = row
    row = row + 1
    self.url_frame.grid_rowconfigure(0, weight=1)
    self.url_frame.grid_columnconfigure(0, weight=1)

    self.grid_rowconfigure(row, weight=1)
    headings = ('Number', 'Serial', 'Name', 'Path', 'Number of\nstorages', 'Number of\ndata locations')
    justifies = (Tkinter.CENTER, Tkinter.CENTER, Tkinter.LEFT, Tkinter.LEFT, Tkinter.CENTER, Tkinter.CENTER)
    self.urlNameWidget = Entry(popup, returnCallback=self.setUrlName)
    editWidgets = [ None, None, self.urlNameWidget, None, None, None ]
    editGetCallbacks = [ None, None, self.getUrlName, self.getUrlPath, None, None ]
    editSetCallbacks = [ None, None, self.setUrlName, None, None, None ]
    self.url_list = ScrolledMatrix(self.url_frame, headingList=headings,
                        justifyList=justifies, initialRows=5, callback=self.selectUrl,
                        editWidgets=editWidgets,
                        editGetCallbacks=editGetCallbacks,
                        editSetCallbacks=editSetCallbacks)
    self.url_list.grid(row=0, column=0, sticky=Tkinter.NSEW)

    texts = [ 'Create', 'Delete' ]
    commands = [ self.createUrl, self.deleteUrl ]
    self.url_buttons = ButtonList(self.url_frame, texts=texts, commands=commands)
    self.url_buttons.grid(row=1, column=0, sticky=Tkinter.EW)

    if (show_storages):

 
      label = ToggleLabel(self, text='Storages ', isArrowClosed=storage_closed,
                          callback=self.storageToggle)
      label.grid(row=row, column=0, sticky=Tkinter.W)
      self.storage_radio = RadioButtons(self, entries=self.all_used_entries,
                                        select_callback=self.updateStorageList,
                                        selected_index=1)
      row = row + 1

      self.storage_frame = Frame(self)
      self.storage_row = row
      row = row + 1
      self.storage_frame.grid_rowconfigure(0, weight=1)
      self.storage_frame.grid_columnconfigure(0, weight=1)

      headings = ('Number', 'Type', 'Package', 'Url serial', 'Storage serial',
                  'Path', 'Is modifiable', 'Is stored', 'Is loaded', 'Is modified')
      justifies = (Tkinter.CENTER, Tkinter.CENTER, Tkinter.LEFT, Tkinter.CENTER, Tkinter.CENTER,
                   Tkinter.LEFT, Tkinter.CENTER, Tkinter.CENTER, Tkinter.CENTER, Tkinter.CENTER)
      self.storageUrlSerialWidget = PulldownMenu(popup, callback=self.setStorageUrlIndSerial)
      self.storagePathWidget = Entry(popup, returnCallback=self.setStoragePath)
      editWidgets = [ None, None, None, self.storageUrlSerialWidget, None,
                      self.storagePathWidget, None, None, None, None ]
      editGetCallbacks = [ None, None, None, self.getStorageUrlSerial, None,
                           self.getStoragePath, self.toggleIsModifiable, None, None, None ]
      editSetCallbacks = [ None, None, None, self.setStorageUrlSerial, None,
                           self.setStoragePath, None, None, None, None ]
      self.storage_list = ScrolledMatrix(self.storage_frame, headingList=headings,
                            justifyList=justifies, initialRows=8,
                            editWidgets=editWidgets,
                            editGetCallbacks=editGetCallbacks,
                            editSetCallbacks=editSetCallbacks, highlightType=None)
      self.storage_list.grid(row=0, column=0, sticky=Tkinter.NSEW)

      texts = [ 'Refresh' ]
      commands = [ self.updateStorageList ]
      self.storage_buttons = ButtonList(self.storage_frame, texts=texts, commands=commands)
      self.storage_buttons.grid(row=1, column=0, sticky=Tkinter.EW)

    if (show_data_locations):
 
      label = ToggleLabel(self, text='Data locations ', isArrowClosed=data_location_closed,
                          callback=self.dataLocationToggle)
      label.grid(row=row, column=0, sticky=Tkinter.W)
      row = row + 1

      self.data_location_frame = Frame(self)
      self.data_location_row = row
      row = row + 1
      self.data_location_frame.grid_rowconfigure(0, weight=1)
      self.data_location_frame.grid_columnconfigure(0, weight=1)

      headings = ('Number', 'Url serial', 'Path')
      justifies = (Tkinter.CENTER, Tkinter.CENTER, Tkinter.LEFT)
      self.dataLocationUrlSerialWidget = PulldownMenu(popup, callback=self.setDataLocationUrlIndSerial)
      self.dataLocationPathWidget = Entry(popup, returnCallback=self.setDataLocationPath)
      editWidgets = [ None, self.dataLocationUrlSerialWidget, self.dataLocationPathWidget ]
      editGetCallbacks = [ None, self.getDataLocationUrlSerial, self.getDataLocationPath ]
      editSetCallbacks = [ None, self.setDataLocationUrlSerial, self.setDataLocationPath ]
      self.data_location_list = ScrolledMatrix(self.data_location_frame,
                                  headingList=headings, justifyList=justifies,
                                  initialRows=3,
                                  editWidgets=editWidgets,
                                  editGetCallbacks=editGetCallbacks,
                                  editSetCallbacks=editSetCallbacks, highlightType=None)
      self.data_location_list.grid(row=0, column=0, sticky=Tkinter.NSEW)

    self.updateUrlList()
    self.urlToggle(url_closed)
    if (show_storages):
      self.updateStorageList()
      self.storageToggle(storage_closed)
    if (show_data_locations):
      self.updateDataLocationList()
      self.dataLocationToggle(data_location_closed)
    self.setUrlState()
    self.setUrlEntries()

    for func in ('__init__', 'delete', ''):
      Implementation.registerNotify(self.updateUrlList, 'memops.Implementation.Url', func)

    for func in ('__init__', 'delete', 'setUrl'):
      Implementation.registerNotify(self.updateUrlList, 'memops.Implementation.NormalStorage', func)
      Implementation.registerNotify(self.updateUrlList, 'memops.Implementation.ContentStorage', func)

    for func in ('',):
      Implementation.registerNotify(self.updateUrlList, 'memops.Implementation.Project', func)

    for func in ('__init__', 'delete'):
      Implementation.registerNotify(self.setUrlEntries, 'memops.Implementation.Url', func)

    if (self.show_storages):
      for func in ('__init__', 'delete', ''):
        Implementation.registerNotify(self.updateStorageList, 'memops.Implementation.NormalStorage', func)
        Implementation.registerNotify(self.updateStorageList, 'memops.Implementation.ContentStorage', func)
      for func in ('',):
        Implementation.registerNotify(self.updateStorageList, 'memops.Implementation.Project', func)

    if (self.show_data_locations):
      for func in ('__init__', 'delete', ''):
        Implementation.registerNotify(self.updateDataLocationList, 'memops.Implementation.DataLocation', func)

  def destroy(self):

    for func in ('__init__', 'delete', ''):
      Implementation.unregisterNotify(self.updateUrlList, 'memops.Implementation.Url', func)

    for func in ('__init__', 'delete', 'setUrl'):
      Implementation.unregisterNotify(self.updateUrlList, 'memops.Implementation.NormalStorage', func)
      Implementation.unregisterNotify(self.updateUrlList, 'memops.Implementation.ContentStorage', func)

    for func in ('',):
      Implementation.unregisterNotify(self.updateUrlList, 'memops.Implementation.Project', func)

    for func in ('__init__', 'delete'):
      Implementation.unregisterNotify(self.setUrlEntries, 'memops.Implementation.Url', func)

    if (self.show_storages):
      for func in ('__init__', 'delete', ''):
        Implementation.unregisterNotify(self.updateStorageList, 'memops.Implementation.NormalStorage', func)
        Implementation.unregisterNotify(self.updateStorageList, 'memops.Implementation.ContentStorage', func)
      for func in ('',):
        Implementation.unregisterNotify(self.updateStorageList, 'memops.Implementation.Project', func)

    if (self.show_data_locations):
      for func in ('__init__', 'delete', ''):
        Implementation.unregisterNotify(self.updateDataLocationList, 'memops.Implementation.DataLocation', func)

    if (self.dir_popup):
      self.dir_popup.destroy()

    Frame.destroy(self)

  """
  def touchStorage(self, storage):

    if (isinstance(storage, ApiImplementation.ContentStorage)):
      if not storage.headObject: # for some reason the storage is an orphan
        return

    if (storage != self.project):
      if (storage.isStored and not storage.isLoaded):
        storage.load()
    if (storage.isStored):
      # storage.isLoaded should be True here (since would have been loaded above)
      # except for project which has isLoaded always True
      storage.isModifiable = True
      storage.touch() # sets isModified to True
  """

  def urlToggle(self, isClosed):
 
    if (isClosed):
      self.url_radio.grid_forget()
      self.grid_rowconfigure(self.url_row, weight=0)
      self.url_frame.grid_forget()
    else:
      self.url_radio.grid(row=self.url_row-1, column=1, sticky=Tkinter.W)
      self.grid_rowconfigure(self.url_row, weight=1)
      self.url_frame.grid(row=self.url_row, column=0, columnspan=2, sticky=Tkinter.NSEW)

  def updateUrlList(self, *extra):

    ind = self.url_radio.getIndex()
    urls = self.project.urls
    objectList = []
    textMatrix = []
    colorMatrix = []
    n = 0
    for url in urls:
      if ((ind == 0) or Io.isUrlUsed(url)):
        n = n + 1
        objectList.append(url)
        text = []
        colors = []
        text.append(n)
        colors.append(unchangedColor)
        text.append(url.serial)
        colors.append(unchangedColor)
        text.append(url.name)
        colors.append(unchangedColor)
        (path, color) = self.getPathColor(url)
        text.append(path)
        colors.append(color)
        text.append(len(url.storages))
        colors.append(unchangedColor)
        text.append(len(url.dataLocations))
        colors.append(unchangedColor)
        textMatrix.append(text)
        colorMatrix.append(colors)
    self.url_list.update(objectList=objectList, textMatrix=textMatrix, colorMatrix=colorMatrix)

  def selectUrl(self, url, row, col):

    self.setUrlState(url)

  def getUrlName(self, url):

    self.urlNameWidget.set(url.name)

  def setUrlName(self, *extra):

    name = self.urlNameWidget.get()
    url = Io.getIoUrl()
    try:
      url.name = name
    except Implementation.ApiError, e:
      showError('Changing url name', e.error_msg)

  def getUrlPath(self, url):

    if (self.dir_popup):
      self.dir_popup.open()
    else:
      file_types = [ FileType('All', ['*']) ]
      self.dir_popup = FileSelectPopup(self, title='Select directory',
                                       file_types=file_types, show_file=False)
    
    dir = self.dir_popup.getDirectory()
    if (dir and dir != url.path):
      Io.setIoPath(url, dir)
    else:
      Io.clearIoPath(url)
    self.updateUrlList()

  def setUrlState(self, url = None):

    if (url and not url.storages and not url.dataLocations):
      state = Tkinter.NORMAL
    else:
      state = Tkinter.DISABLED

    self.url_buttons.buttons[1].config(state=state)

  def createUrl(self):

    name = askString('URL name', 'What is the name of the new URL?', parent=self)
    if (name):
      try:
        ApiImplementation.Url(self.project, name=name, path=normalisePath(os.getcwd()))
      except Implementation.ApiError, e:
        showError('Creating url', e.error_msg)

  def deleteUrl(self):

    url = self.url_list.currentObject
 
    if (url):
      if (showYesNo('Delete URL',
            'Are you sure you want to delete URL "%s"?' % url.name)):
        url.delete()

  def getUrl(self):

    url = self.url_list.currentObject

    return url

  def setUrlEntries(self, *extra):

    entries = [ str(url.serial) for url in self.project.urls ]
    if (self.show_storages):
      self.storageUrlSerialWidget.replace(entries)
    if (self.show_data_locations):
      self.dataLocationUrlSerialWidget.replace(entries)

  def refresh(self):

    self.updateUrlList()
    self.updateStorageList()
    self.updateDataLocationList()

  def storageToggle(self, isClosed):
 
    if (isClosed):
      self.storage_radio.grid_forget()
      self.grid_rowconfigure(self.storage_row, weight=0)
      self.storage_frame.grid_forget()
    else:
      self.storage_radio.grid(row=self.storage_row-1, column=1, sticky=Tkinter.W)
      self.grid_rowconfigure(self.storage_row, weight=1)
      self.storage_frame.grid(row=self.storage_row, column=0, columnspan=2, sticky=Tkinter.NSEW)

  def updateStorageList(self, *extra):

    ind = self.storage_radio.getIndex()
    storages = (self.project,) + self.project.storages
    objectList = []
    textMatrix = []
    colorMatrix = []
    n = 0
    for storage in storages:
      if ((ind == 0) or not Io.isStorageEmpty(storage)):
        n = n + 1
        objectList.append(storage)
        if (storage == self.project):
          t = 'project'
          serial = None
        elif (isinstance(storage, ApiImplementation.NormalStorage)):
          t = 'normal'
          serial = storage.serial
        elif (isinstance(storage, ApiImplementation.ContentStorage)):
          t = 'content'
          serial = storage.serial
        else:
          t = ''
        text = []
        colors = []
        text.append(n)
        colors.append(unchangedColor)
        text.append(t)
        colors.append(unchangedColor)
        text.append(storage.package)
        colors.append(unchangedColor)
        (url, color) = self.getUrlColor(storage)
        text.append(url.serial)
        colors.append(color)
        text.append(serial)
        colors.append(unchangedColor)
        (path, color) = self.getPathColor(storage)
        text.append(path)
        colors.append(color)
        text.append(self.getYesNo(storage.isModifiable))
        colors.append(unchangedColor)
        text.append(self.getYesNo(storage.isStored))
        colors.append(unchangedColor)
        text.append(self.getYesNo(storage.isLoaded))
        colors.append(unchangedColor)
        text.append(self.getYesNo(storage.isModified))
        colors.append(unchangedColor)
        textMatrix.append(text)
        colorMatrix.append(colors)
    self.storage_list.update(objectList=objectList, textMatrix=textMatrix, colorMatrix=colorMatrix)

  def getYesNo(self, value):

    if (value):
      return 'yes'
    else:
      return 'no'

  def getStorageUrlSerial(self, storage):

    url = Io.getIoUrl(storage)
    self.storageUrlSerialWidget.setSelected(str(url.serial))

  def setStorageUrlSerial(self, storage):

    serial = self.storageUrlSerialWidget.getSelected()
    self.setStorageUrlIndSerial(-1, serial) # -1 arbitrary

  def setStorageUrlIndSerial(self, serial_ind, serial):

    storage = self.getStorage()
    if (storage):
      serial = int(serial)
      url = self.project.findFirstUrl(serial=serial)
      if (url and storage.url != url):
        Io.setIoUrl(storage, url)
      else:
        Io.clearIoUrl(storage)
      self.updateStorageList()

  def getStoragePath(self, storage):

    path = Io.getIoPath(storage)
    self.storagePathWidget.set(path)

  def setStoragePath(self, *extra):

    storage = self.getStorage()
    if (storage):
      path = self.storagePathWidget.get()
      if (path and storage.path != path):
        Io.setIoPath(storage, path)
      else:
        Io.clearIoPath(storage)
      self.updateStorageList()

  def toggleIsModifiable(self, storage):

    try:
      storage.isModifiable = not storage.isModifiable
    except Implementation.ApiError, e:
      showError('Changing storage isModifiable', e.error_msg)

  def getStorage(self):

    storage = self.storage_list.currentObject

    return storage

  def dataLocationToggle(self, isClosed):
 
    if (isClosed):
      self.grid_rowconfigure(self.data_location_row, weight=0)
      self.data_location_frame.grid_forget()
    else:
      self.grid_rowconfigure(self.data_location_row, weight=1)
      self.data_location_frame.grid(row=self.data_location_row, column=0, columnspan=2, sticky=Tkinter.NSEW)

  def updateDataLocationList(self, *extra):

    data_locations = self.project.dataLocations
    textMatrix = []
    colorMatrix = []
    n = 0
    for data_location in data_locations:
      n = n + 1
      text = []
      colors = []
      text.append(n)
      colors.append(unchangedColor)
      (url, color) = self.getUrlColor(data_location)
      text.append(url.serial)
      colors.append(color)
      (path, color) = self.getPathColor(data_location)
      text.append(path)
      colors.append(color)
      textMatrix.append(text)
      colorMatrix.append(colors)
    self.data_location_list.update(objectList=data_locations, textMatrix=textMatrix, colorMatrix=colorMatrix)

  def getDataLocationUrlSerial(self, data_location):

    url = Io.getIoUrl(data_location)
    self.dataLocationUrlSerialWidget.setSelected(str(url.serial))

  def setDataLocationUrlSerial(self, data_location):

    serial = self.dataLocationUrlSerialWidget.getSelected()
    self.setDataLocationUrlIndSerial(-1, serial) # -1 arbitrary

  def setDataLocationUrlIndSerial(self, serial_ind, serial):

    data_location = self.getDataLocation()
    if (data_location):
      serial = int(serial)
      url = self.project.findFirstUrl(serial=serial)
      if (url and data_location.url != url):
        Io.setIoUrl(data_location, url)
      else:
        Io.clearIoUrl(data_location)
      self.updateDataLocationList()

  def getDataLocationPath(self, data_location):

    path = Io.getIoPath(data_location)
    self.dataLocationPathWidget.set(path)

  def setDataLocationPath(self, *extra):

    data_location = self.getDataLocation()
    if data_location:
      path = self.dataLocationPathWidget.get()
      if (path and data_location.path != path):
        Io.setIoPath(data_location, path)
      else:
        Io.clearIoPath(data_location)
      self.updateDataLocationList()

  def getDataLocation(self):

    data_location = self.data_location_list.currentObject

    return data_location

  def getPathColor(self, obj):

    path = Io.getIoPath(obj)

    if path == obj.path:
      color = unchangedColor
    else:
      color = changedColor

    return (path, color)

  def getUrlColor(self, obj):

    url = Io.getIoUrl(obj)

    if url == obj.url:
      color = unchangedColor
    else:
      color = changedColor

    return (url, color)

if (__name__ == '__main__'):

  p = ApiImplementation.MemopsRoot(name='testMemopsRoot')
  nmrProject = Nmr.NmrProject(p, name = p.name)
  r = Tkinter.Tk()
  r.grid_rowconfigure(0, weight=1)
  r.grid_columnconfigure(0, weight=1)
  r.top = r # needed so that getPopup() works
  f = UrlViewerFrame(r, p)
  f.grid(sticky=Tkinter.NSEW)
  r.mainloop()

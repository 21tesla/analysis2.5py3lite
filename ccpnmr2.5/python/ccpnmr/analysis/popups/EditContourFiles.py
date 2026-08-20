
"""
======================COPYRIGHT/LICENSE START==========================

EditContourFiles.py: Part of the CcpNmr Analysis program

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

from memops.universal.Io import joinPath

from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Entry import Entry
from memops.gui.MessageReporter import showYesNo, showError
from memops.gui.ScrolledMatrix import ScrolledMatrix

from ccpnmr.analysis.popups.BasePopup import BasePopup

class EditContourFilesPopup(BasePopup):

  """
  **Manage Contour Files in Project**

  The purpose of this dialog is to allow the user to manage contour files.
  A contour file stores contours calculated at specified contour levels in
  a specified region of a spectrum, as an alternative to calculating contours
  on the fly.  You cannot pick peaks using contour files but after peaks have
  been picked, contour files are a reasonable alternative to contours for large
  3D and 4D data sets, since the spectrum data does not need to be read from disk.

  This dialog lists the contour files that are in the project.  To create a
  new contour file click on the "Create New File" button and to add a contour
  file that already exists on disk (e.g. from another project) then click on
  the "Find Existing File" button.

  When a contour file is created for a spectrum then by default it is used to
  display contours rather than having them calculated on the fly.  This can be
  changed in the main `Spectra`_ dialog in the Display Options tab.

  See also: `Creating Contour Files`_, `Add Existing Contour Files`_.

  .. _`Spectra`: EditSpectrumPopup.html
  .. _`Creating Contour Files`: CreateContourFilePopup.html
  .. _`Add Existing Contour Files`: AddContourFilePopup.html
"""

  def __init__(self, parent, *args, **kw):

    self.waiting = False
    self.contourFileDir = None

    BasePopup.__init__(self, parent=parent, title='Spectrum Contour Files', **kw)

  def body(self, master):

    self.geometry('650x200')

    master.grid_columnconfigure(0, weight=1)

    row = 0
    master.grid_rowconfigure(row, weight=1)
    self.storedContourTable = None
    headings = ('#', 'Experiment', 'Spectrum', 'Dims', 'Directory', 'File')
    ###self.urlWidget = Entry(self, returnCallback=self.setUrl, width=10)
    ###self.pathWidget   = Entry(self, returnCallback=self.setPath, width=10)
    ###editWidgets = [ None, None, None, None, self.urlWidget, self.pathWidget ]
    ###editGetCallbacks = [ None, None, None, None, self.getUrl, self.getPath ]
    ###editSetCallbacks = [ None, None, None, None, self.setUrl, self.setPath ]
    tipTexts = ['Row number',
                'The experiment that contains the spectrum to which the contour file relates',
                'The spectrum to which the contour file relates',
                'The projected spectrum dimensions which are represented by the contour file, i.e. to form the X-Y plane of the screen',
                'Directory in which the contour file is saved',
                'The name of the file in which contours for thei spectrum projection are saved']
    editWidgets = 6 * [None]
    editGetCallbacks = 6 * [None]
    editSetCallbacks = 6 * [None]
    self.storedContourTable = ScrolledMatrix(master, headingList=headings,
                                        callback=self.setButtonState,
                                        editWidgets=editWidgets, multiSelect=True,
                                        editGetCallbacks=editGetCallbacks,
                                        editSetCallbacks=editSetCallbacks,
                                        deleteFunc=self.deleteContourFile,
                                        tipTexts=tipTexts)
    self.storedContourTable.grid(row=row, column=0, sticky='nsew')
    row = row + 1

    tipTexts = ['Locate an exiting file on disk and use that as a spectrum contour file',
                'Create a new spectrum contour file; requires that the spectrum data matrix is accessible',
                'Remove the selected contour file from the CCPN project and optionally delete the file from disk']
    texts = [ 'Find Existing File', 'Create New File', 'Delete' ]
    commands = [ self.addContourFile, self.createContourFile, self.deleteContourFile ]
    self.buttons = UtilityButtonList(master, texts=texts, doClone=False,
                                     commands=commands, helpUrl=self.help_url)
    self.buttons.grid(row=row, column=0, sticky='ew')

    self.updateAfter()

    for clazz in ('ccp.nmr.Nmr.Experiment', 'ccp.nmr.Nmr.DataSource'):
      self.registerNotify(self.updateAfter, clazz, 'setName')

    for func in ('__init__', 'delete', ''):
      self.registerNotify(self.updateAfter, 'ccpnmr.Analysis.StoredContour', func)

    ###self.registerNotify(self.updateAfter, 'memops.Implementation.Url', 'setPath')

  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)
    
  def destroy(self):

    for clazz in ('ccp.nmr.Nmr.Experiment', 'ccp.nmr.Nmr.DataSource'):
      self.registerNotify(self.updateAfter, clazz, 'setName')

    for func in ('__init__', 'delete', ''):
      self.unregisterNotify(self.updateAfter, 'ccpnmr.Analysis.StoredContour', func)

    ###self.unregisterNotify(self.updateAfter, 'memops.Implementation.Url', 'setPath')

    BasePopup.destroy(self)

  def updateAfter(self, *extra):
 
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)

  def update(self):

    self.updateStoredContourTable()
    self.setButtonState()
    self.waiting = False

  def updateStoredContourTable(self, *extra):

    analysisSpectra = self.analysisProject.analysisSpectra
      
    objectList = []
    textMatrix = []
    n = 0
    for analysisSpectrum in analysisSpectra:
      spectrum = analysisSpectrum.dataSource
      for storedContour in analysisSpectrum.storedContours:
        if storedContour.isDeleted:
          continue

        n = n + 1
        dims = ', '.join(['%d' % dim for dim in storedContour.dims])

        text = [n, spectrum.experiment.name, spectrum.name,
                dims, analysisSpectrum.contourDir.dataLocation, storedContour.path]
  
        objectList.append(storedContour)
        textMatrix.append(text)

    self.storedContourTable.update(objectList=objectList, textMatrix=textMatrix)

  """
  def getUrl(self, storedContour):

    self.urlWidget.set(storedContour.url.path)

  def setUrl(self, *extra):

    storedContour = self.getStoredContour()
    if not storedContour:
      return
    
    path = self.urlWidget.get()
    if not path:
      return

    url = storedContour.url
    project = self.project
    analysisProject = self.analysisProject
    if len(analysisProject.findAllStoredContours(url=url)) > 1 or \
        project.findAllStorages(url=url) or project.findAllDataLocations(url=url):
      # create a new url with new path
      url = project.newUrl(name='storedContourUrl', path=path)
    else:
      # just change url path since (hopefully) nobody else is using it
      try:
        url.path = path
      except Implementation.ApiError, e:
        showError('Setting url path', e.error_msg, parent=self)

  def getPath(self, storedContour):

    self.pathWidget.set(storedContour.path)

  def setPath(self, *extra):

    storedContour = self.getStoredContour()
    if not storedContour:
      return
    
    path = self.pathWidget.get()
    if not path:
      return

    try:
      storedContour.path = path
    except Implementation.ApiError, e:
      showError('Setting stored contour path', e.error_msg, parent=self)
"""

  def setButtonState(self, *extra):

    if self.storedContourTable.objectList:
      state = 'normal'
    else:
      state = 'disabled'
 
    self.buttons.buttons[2].config(state=state)

  def getSpectrum(self):

    storedContour = self.getStoredContour()
    if storedContour:
      try:
        spectrum = storedContour.dataSource
      except:
        spectrum = None
    else:
      spectrum = None

    return spectrum

  def getStoredContour(self):

    if self.storedContourTable:
      storedContour = self.storedContourTable.currentObject
      if storedContour and storedContour.isDeleted:
        storedContour = None
    else:
      storedContour = None

    return storedContour

  def getStoredContours(self):

    if (self.storedContourTable):
      storedContours = self.storedContourTable.currentObjects

    return storedContours

  def addContourFile(self, *extra):

    spectrum = self.getSpectrum()
    self.parent.addSpectrumContourFile(spectrum)

  def createContourFile(self, *extra):

    spectrum = self.getSpectrum()
    self.parent.createSpectrumContourFile(spectrum)

  def deleteContourFile(self, *extra):

    storedContours = self.getStoredContours()
    if not storedContours:
      return

    if len(storedContours) == 1:
      s = t = ''
    else:
      s = ' %d ' % len(storedContours)
      t = 's'
    if (showYesNo('Delete stored contour%s from project' % t,
            'Are you sure you want to delete selected %sstored contour%s from the project?' % (s, t), parent=self)):
      if (showYesNo('Delete underlying contour file%s on disk' % t,
            'Do you also want to delete underlying contour file%s on disk (normally yes)?' % t, parent=self)):
        for storedContour in storedContours:
          path = storedContour.fullPath
          try:
            os.remove(path)
          except:
            print 'Warning: could not remove %s' % path
      for storedContour in storedContours:
        storedContour.delete()


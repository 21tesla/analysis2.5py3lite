"""
======================COPYRIGHT/LICENSE START==========================

AddContourFile.py: Part of the CcpNmr Analysis program

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
import shutil

from memops.universal.Io import normalisePath

from memops.gui.Button import Button
from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Entry import Entry
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError, showInfo
from memops.gui.PulldownList import PulldownList

from ccp.api.general.DataLocation import BlockedBinaryMatrix

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.ContourStore import getStoredContourHeader, createStoredContour

class AddContourFilePopup(BasePopup):

  """
  **Add Existing Contour Files to Project**

  The purpose of this dialog is to allow the user to add pre-existing
  contour files to the project.  Contour files only depend on the spectrum
  data so the same contour files can be used across multiple projects,
  and that is the reason this dialog might be used.

  See also: `Spectrum Contour Files`_, `Creating Contour Files`_.

  .. _`Spectrum Contour Files`: EditContourFilesPopup.html
  .. _`Creating Contour Files`: CreateContourFilePopup.html
"""

  def __init__(self, parent, *args, **kw):

    self.spectrum = None

    BasePopup.__init__(self, parent=parent, title='Add existing contour file', **kw)

  def body(self, master):

    self.geometry('600x130')
    master.grid_columnconfigure(1, weight=1)
    for n in range(5):
      master.grid_rowconfigure(n, weight=1)

    row = 0
    label = Label(master, text='Spectrum: ')
    label.grid(row=row, column=0, sticky='e')
    tipText = 'The spectrum for which the contour file is being added'
    self.expt_spectrum = PulldownList(master, callback=self.updateContourDir,
                                      tipText=tipText)
    self.expt_spectrum.grid(row=row, column=1, sticky='w')

    row = row + 1
    tipText = 'The location of the directory where contour files are stored on disk'
    label = Label(master, text='Contour dir: ')
    label.grid(row=row, column=0, sticky='e')
    self.dir_label = Label(master, text='', tipText=tipText)
    self.dir_label.grid(row=row, column=1, sticky='w')

    row = row + 1
    label = Label(master, text='(file will be copied into Contour dir if it is not already in there)')
    label.grid(row=row, column=1, sticky='w')

    row = row + 1
    tipText = 'Browse for a file store contour data'
    button = Button(master, text='File name: ', command=self.selectFile, tipText=tipText)
    button.grid(row=row, column=0, sticky='e')
    tipText = 'Enter the name of the file to store contour data'
    self.file_entry = Entry(master, tipText=tipText)
    self.file_entry.grid(row=row, column=1, sticky='ew')

    row = row + 1
    texts = [ 'Add File' ]
    commands = [ self.addFile ]
    tipTexts = ['Use the selected contour file in the current project, copying it to the contour directory if required',]
    self.buttons = UtilityButtonList(master, texts=texts, doClone=False, tipTexts=tipTexts,
                                     commands=commands, helpUrl=self.help_url)
    self.buttons.grid(row=row, column=0, columnspan=2, sticky='ew')

    self.curateNotifiers(self.registerNotify)
    self.updateSpectrum()

  def destroy(self):

    self.curateNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

  def curateNotifiers(self, notifyFunc):

    for clazz in ('Experiment', 'DataSource'):
      for func in ('__init__', 'delete', 'setName'):
        notifyFunc(self.updateNotifier, 'ccp.nmr.Nmr.%s' % clazz, func)

  def updateSpectrum(self, spectrum=None):

    if not spectrum:
      spectrum = self.spectrum

    spectra = self.parent.getSpectra()
    if spectra:
      if spectrum not in spectra:
        spectrum = spectra[0]
      index = spectra.index(spectrum)
      names = ['%s:%s' % (x.experiment.name, x.name) for x in spectra]
    else:
      index = 0
      names = []

    self.expt_spectrum.setup(names, spectra, index)

    self.updateContourDir(spectrum)

  def updateNotifier(self, *extra):

    self.updateSpectrum()

  def updateContourDir(self, spectrum):

    if spectrum is self.spectrum:
      return

    self.spectrum = spectrum

    if spectrum:
      path = spectrum.analysisSpectrum.contourDir.dataLocation
    else:
      path = ''
    self.dir_label.set(path)

  def selectFile(self):

    spectrum = self.spectrum
    if spectrum:
      directory = spectrum.analysisSpectrum.contourDir.dataLocation
    else:
      directory = os.getcwd()
    popup = FileSelectPopup(self, directory=directory)
    fileName = popup.getFile()
    popup.destroy()
    if fileName:
      self.file_entry.set(fileName)

  def addFile(self):

    spectrum = self.spectrum
    if not spectrum:
      return

    dataStore = spectrum.dataStore
    if not dataStore:
      showError('No dataStore', 'Spectrum does not have associated dataStore', parent=self)
      return

    if not isinstance(dataStore, BlockedBinaryMatrix):
      showError('No blockedBinaryMatrix', 'Spectrum dataStore is not a blockedBinaryMatrix', parent=self)
      return

    if not dataStore.blockSizes:
      showError('No blockSize', 'Spectrum dataStore does not have blockSize set', parent=self)
      return

    blockSize = list(dataStore.blockSizes)

    fileName = self.file_entry.get()
    if not fileName:
      showError('No filename', 'No filename given', parent=self)
      return

    fileName = normalisePath(fileName, makeAbsolute=True)

    contourDir = spectrum.analysisSpectrum.contourDir.dataLocation
    if fileName.startswith(contourDir):
      path = fileName[len(contourDir)+1:]
    else:
      path = os.path.basename(fileName)
      if not os.path.exists(contourDir):
        os.makedirs(contourDir)
      print 'Copying %s to %s' % (fileName, contourDir)
      shutil.copy(fileName, contourDir)

    try:
      header = getStoredContourHeader(fileName)
    except Exception, e:
      showError('File error', str(e), parent=self)
      return

    if header['ndim'] != spectrum.numDim:
      showError('Number of dimensions', 'Number of dimensions in file (%d) does not match spectrum (%d)' % (header['ndim'], spectrum.numDim), parent=self)
      return

    dataDims = spectrum.sortedDataDims()
    npoints = [x.numPoints for x in dataDims]
    if header['npoints'] != npoints:
      showError('Number of points', 'Number of points in file (%s) does not match spectrum (%s)' % (header['npoints'], npoints), parent=self)
      return

    if header['blockSize'] != blockSize:
      showError('Block size', 'Block size in file (%s) does not match spectrum (%s)' % (header['blockSize'], blockSize), parent=self)
      return

    createStoredContour(spectrum, path, header['xdim'], header['ydim'])

    showInfo('Added file', 'Successfully added stored contour file', parent=self)

  def setSpectrum(self, spectrum):

    self.expt_spectrum.set(spectrum, doCallback=True)


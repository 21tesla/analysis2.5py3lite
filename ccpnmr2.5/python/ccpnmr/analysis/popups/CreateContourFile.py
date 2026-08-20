"""
======================COPYRIGHT/LICENSE START==========================

CreateContourFile.py: Part of the CcpNmr Analysis program

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
import re

from memops.universal.Io import normalisePath, joinPath

from memops.api.Implementation import Url

from memops.gui.Button import Button
from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Entry import Entry
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError, showYesNo, showInfo
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.ContourStore import saveSpectrumContours, createStoredContour
from ccpnmr.analysis.core import UnitConverter
from ccpnmr.analysis.core import ExperimentBasic

class CreateContourFilePopup(BasePopup):

  """
  **Create Contour Files for Project**

  The purpose of this dialog is to allow the user to create new contour
  files for the project.  Not only the spectrum has to be specified, but
  also the two dimensions (X, Y) that are going to be contoured.  The
  table specifies the conditions that specify the region being contoured,
  in terms of what is included and excluded.

  For now, only one condition (which is an "include" condition) is allowed.
  This also means that the "Add Condition" and "Delete Condition" buttons
  are disabled.

  The contour directory is specified by the program, using the project
  directory but the file name can be specified by the user.

  The contour levels used are the ones that are current for the spectrum,
  as set in the `Contour Levels dialog`_.

  See also: `Spectrum Contour Files`_, `Add Existing Contour Files`_.

  .. _`Contour Levels dialog`: EditContourLevelsPopup.html
  .. _`Spectrum Contour Files`: EditContourFilesPopup.html
  .. _`Add Existing Contour Files`: AddContourFilePopup.html
"""

  def __init__(self, parent, *args, **kw):

    self.spectrumConditions = []
    self.col = None
    self.spectrum = None

    BasePopup.__init__(self, parent=parent, title='New Contour File', **kw)

  def body(self, master):

    self.geometry('650x200')

    master.grid_columnconfigure(1, weight=1)

    row = 0
    frame = Frame(master, grid=(row, 0), gridSpan=(1,2))
    
    label = Label(frame, text='Spectrum: ', grid=(row, 0))
    tipText = 'Selects the experiment and spectrum for which to make a contour file'
    self.expt_spectrum = PulldownList(frame, grid=(row, 1),
                                      tipText=tipText, callback=self.update)

    label = Label(frame, text=' (X-dim, Y-dim): ', grid=(row, 2))
    tipText = 'Selects which dimensions (projection) of the spectrum form the X and Y axes of the contour file'
    self.dim_menu = PulldownList(frame, grid=(row, 3), tipText=tipText,
                                 callback=self.updateFile)

    row = row + 1
    master.grid_rowconfigure(row, weight=1)
    #### for now not editable because only allow one include region
    ###self.conditionMenu = PulldownMenu(self, entries=('include', 'exclude'),
    ###                       callback=self.selectedCondition, do_initial_callback=False)
    self.regionEntry = FloatEntry(self, text='', returnCallback=self.setRegion, width=10)
    tipTexts = ['Whether to include or exclude the region in the contour file',
                '',
                '',
                '']
    headingList = ['Condition','','','','']
    editSetCallbacks = [None]
    editGetCallbacks = [None]
    ###editWidgets = [self.conditionMenu]
    editWidgets = [None]
    self.conditionTable = ScrolledMatrix(master, initialRows=6,
                                         grid=(row, 0), gridSpan=(1,2),
                                         tipTexts=tipTexts,
                                         headingList=headingList,
                                         callback=self.selectCell,
                                         editWidgets=editWidgets,
                                         editGetCallbacks=editGetCallbacks,
                                         editSetCallbacks=editSetCallbacks)

    # TBD: make directory editable
    row = row + 1
    label = Label(master, text='Contour dir: ', grid=(row, 0), sticky='e')
    tipText = 'The directory location on disk into which contour files are saved'
    self.dir_label = Label(master, text='', grid=(row, 1), tipText=tipText)

    row = row + 1
    label = Label(master, text='File name: ', grid=(row, 0), sticky='e')
    tipText = 'Sets the name of the contour file to save to disk'
    self.file_entry = Entry(master, grid=(row, 1), tipText=tipText)

    ##row = row + 1
    ##label = Label(master, text='Contour levels: ')
    ##label.grid(row=row, column=0, sticky='e')
    ##self.levels_entry = FloatEntry(master, isArray=True, returnCallback=self.saveLevels)
    ##self.levels_entry.grid(row=row, column=1, sticky='ew')

    row = row + 1
    tipTexts = ['Add a new contour region for including in or excluding from the contour file ',
                'Remove the selected contour region from the table',
                'Make the specified contour file using the input settings & regions']
    texts = [ 'Add Condition', 'Delete Condition', 'Contour and Save' ]
    commands = [ self.addCondition, self.deleteCondition, self.contourAndSaveSpectrum ]
    self.buttons = UtilityButtonList(master, texts=texts, commands=commands,
                                     doClone=False, helpUrl=self.help_url,
                                     grid=(row, 0), gridSpan=(1,2), tipTexts=tipTexts)

    ### TBD disabled for now
    for n in range(2):
      self.buttons.buttons[n].config(state='disabled')

    self.curateNotifiers(self.registerNotify)
    self.updateSpectrum()
    self.updateDimMenu()

  def destroy(self):

    self.curateNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

  def curateNotifiers(self, notifyFunc):

    for clazz in ('Experiment', 'DataSource'):
      for func in ('__init__', 'delete', 'setName'):
        notifyFunc(self.updateNotifier, 'ccp.nmr.Nmr.%s' % clazz, func)

  def update(self, spectrum):

    if spectrum is self.spectrum:
      return

    self.spectrum = spectrum

    self.updateDimMenu()
    self.updateConditionTable()
    self.updateFile()
    ##self.updateLevels()

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

    self.update(spectrum)

  def updateNotifier(self, *extra):

    self.updateSpectrum()

  def updateFile(self, *extra):

    spectrum = self.spectrum
    if not spectrum:
      return

    analysisSpectrum = spectrum.analysisSpectrum
    if not analysisSpectrum:
      return

    dims = self.dim_menu.getText()
    dims = [int(dim) for dim in dims.split(', ')]
    (xdim, ydim) = dims
    storedContour = analysisSpectrum.findFirstStoredContour(dims=dims)
    if not storedContour:
      storedContour = analysisSpectrum.findFirstStoredContour()
    if storedContour:
      fileName = storedContour.path
    else:
      fileName = '%s_%s_%d_%d.cnt' % (spectrum.experiment.name, spectrum.name, xdim, ydim)

    path = self.getContourDir()
    self.dir_label.set(path)
    self.file_entry.set(fileName)

  def getContourDir(self):

    spectrum = self.spectrum
    if not spectrum:
      return ''

    analysisSpectrum = spectrum.analysisSpectrum
    if not analysisSpectrum:
      return ''

    url = analysisSpectrum.contourDir
    if url: # should be the case since set in Analysis.py
      path = url.dataLocation
    else:
      path = ''

    return path

  """
  def saveFile(self, *extra):

    spectrum = self.spectrum
    if not spectrum:
      return

    analysisSpectrum = spectrum.analysisSpectrum
    if not analysisSpectrum:
      return

    fileName = self.file_entry.get()
    if fileName:
      
      application = self.project.application
      application.setValue(spectrum, keyword='contourFileName', value=fileName)
"""
    
  def updateDimMenu(self):

    entries = []
    spectrum = self.spectrum
    if spectrum:
      ndim = spectrum.numDim
      for xdim in range(1, ndim):
        dataDim = spectrum.findFirstDataDim(dim=xdim)
        if dataDim.className != 'SampledDataDim':
          for ydim in range(xdim+1, ndim+1):
            dataDim = spectrum.findFirstDataDim(dim=ydim)
            if dataDim.className != 'SampledDataDim':
              entries.append('%d, %d' % (xdim, ydim))
    self.dim_menu.setup(entries, objects=None, index=0)

  def updateConditionTable(self, *extra):

    spectrum = self.spectrum
    if spectrum:
      ndim = spectrum.numDim
    else:
      ndim = 0
    tipTexts = ['Whether to include or exclude the region in the contour file',]
    headingList = ['Condition']
    ###editWidgets = [self.conditionMenu] + 2*ndim*[self.regionEntry]
    ###editGetCallbacks = [self.getCondition]
    ###editSetCallbacks = [self.setCondition]
    editWidgets = [None] + 2*ndim*[self.regionEntry]
    editGetCallbacks = [None]
    editSetCallbacks = [None]
    for i in range(ndim):
      dim = i + 1
      headingList.extend(['Dim %d Min' % dim, 'Dim %d Max' % dim])
      editGetCallbacks.append(lambda obj, i=i: self.getRegionMin(obj, i))
      editGetCallbacks.append(lambda obj, i=i: self.getRegionMax(obj, i))
      editSetCallbacks.append(lambda obj, i=i: self.setRegionMin(obj, i))
      editSetCallbacks.append(lambda obj, i=i: self.setRegionMax(obj, i))
      tipTexts.append('Lower bound of dimension %d region to include or exclude' % dim)
      tipTexts.append('Upper bound of dimension %d region to include or exclude' % dim)
      
    objectList = []
    textMatrix = []
    self.spectrumConditions = self.getSpectrumConditions()
    for spectrumCondition in self.spectrumConditions:
      (condition, region) = spectrumCondition
      objectList.append(spectrumCondition)
      textRow = [condition]
      for i in range(ndim):
        r = region[i]
        if r:
          rMin = r[0]
          rMax = r[1]
        else:
          rMin = rMax = None
        textRow.append(rMin)
        textRow.append(rMax)
      textMatrix.append(textRow)

    self.conditionTable.update(objectList=objectList, textMatrix=textMatrix,
                               headingList=headingList,
                               editSetCallbacks=editSetCallbacks,
                               editGetCallbacks=editGetCallbacks,
                               editWidgets=editWidgets)

  """
  def updateLevels(self):

    spectrum = self.getSpectrum()
    if spectrum:
      analysisSpectrum = spectrum.analysisSpectrum
      levels = analysisSpectrum.posLevels + analysisSpectrum.negLevels
      scale = spectrum.scale / self.analysisProject.globalContourScale
      levels = [ level/scale for level in levels ]
    else:
      levels = []

    self.levels_entry.set(levels)

  def saveLevels(self):

    spectrum = self.spectrum
    if not spectrum:
      return
"""

  def getRegionMin(self, spectrumCondition, dim):
   
    #print 'getRegionMin'
    (condition, region) = spectrumCondition
    self.regionEntry.set(region[dim][0])

  def setRegionMin(self, spectrumCondition, dim):
   
    #print 'setRegionMin'

    (condition, region) = spectrumCondition
    (r0, r1) = region[dim]
    r = self.regionEntry.get()
    region[dim] = (r, r1)
    self.setSpectrumConditions()
    self.updateConditionTable()

  def getRegionMax(self, spectrumCondition, dim):
   
    #print 'getRegionMax'
    (condition, region) = spectrumCondition
    self.regionEntry.set(region[dim][1])

  def setRegionMax(self, spectrumCondition, dim):

    #print 'setRegionMax'
    (condition, region) = spectrumCondition
    (r0, r1) = region[dim]
    r = self.regionEntry.get()
    region[dim] = (r0, r)
    self.setSpectrumConditions()
    self.updateConditionTable()

  def setRegion(self, *extra):

    spectrumCondition = self.getCurrentObject()
    if spectrumCondition and self.col is not None:
      col = self.col - 1 # -1 because of condition
      dim = col / 2
      if col % 2: # max
        self.setRegionMax(spectrumCondition, dim)
      else: # min
        self.setRegionMin(spectrumCondition, dim)

  """ not needed for now
  def getCondition(self, spectrumCondition):
   
    #print 'getCondition'
    (condition, region) = spectrumCondition
    self.conditionMenu.set(condition)

  def setCondition(self, spectrumCondition):

    #print 'setCondition'
    spectrumCondition[0] = self.conditionMenu.get()
    self.setSpectrumConditions()

  def selectedCondition(self, ind, condition):

    spectrumCondition = self.getCurrentObject()

    if spectrumCondition is not None:
      self.setCondition(spectrumCondition)
"""

  def getSpectrumConditions(self):

    spectrum = self.spectrum
    if not spectrum:
      return []

    application = self.project.application
    spectrumConditions = application.getValue(spectrum, keyword='contourFileConditions')
    if spectrumConditions:
      spectrumConditions = eval(spectrumConditions)
    else:
      region = []
      for i in range(spectrum.numDim):
        dim = i + 1
        region.append(self.getWholeRegion(spectrum, dim))
      spectrumCondition = ['include', region]
      spectrumConditions = [spectrumCondition]
      self.setSpectrumConditions(spectrumConditions)

    return spectrumConditions

  def setSpectrumConditions(self, spectrumConditions = None):

    spectrum = self.spectrum

    if spectrum:
      if not spectrumConditions:
        spectrumConditions = self.spectrumConditions

      application = self.project.application
      application.setValue(spectrum, keyword='contourFileConditions', value=str(spectrumConditions))

  def getDimMin(self, spectrum, dim):

    dataDim = spectrum.findFirstDataDim(dim=dim)
    if dataDim.className == 'SampledDataDim':
      r = 1.0
    else:
      converter = UnitConverter.pnt2ppm
      dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
      r = converter(float(dataDim.numPoints), dataDimRef)

    return r

  def getDimMax(self, spectrum, dim):

    dataDim = spectrum.findFirstDataDim(dim=dim)
    if dataDim.className == 'SampledDataDim':
      r = float(dataDim.numPoints)
    else:
      converter = UnitConverter.pnt2ppm
      dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
      r = converter(1.0, dataDimRef)

    return r

  def getWholeRegion(self, spectrum, dim):

    rMin = self.getDimMin(spectrum, dim)
    rMax = self.getDimMax(spectrum, dim)

    return (rMin, rMax)

  def getCurrentObject(self):

    # sometimes row highlighting stops so obj passed into selectCell might no longer be current
    # so instead use direct interrogation of scrolledMatrix
    obj = self.conditionTable.currentObject

    return obj

  def selectCell(self, obj, row, col):

    # see note about obj in getCurrentObject()
    self.col = col

  def addCondition(self):

    spectrum = self.spectrum
    if spectrum:
      ndim = spectrum.numDim
      spectrumCondition = ['exclude', ndim*[(None,None)]]
      self.spectrumConditions.append(spectrumCondition)
      self.setSpectrumConditions()
      self.updateConditionTable()

  def deleteCondition(self):

    spectrumCondition = self.getCurrentObject()
    if spectrumCondition:
      self.spectrumConditions.remove(spectrumCondition)
      self.setSpectrumConditions()
      self.updateConditionTable()

  def continueIfFileNameUsed(self, fileName):

    result = True

    storedContoursToDelete = []
    for analysisSpectrum in self.analysisProject.analysisSpectra:
      for storedContour in analysisSpectrum.storedContours:
        if storedContour.fullPath == fileName:
          storedContoursToDelete.append(storedContour)

    if storedContoursToDelete:
      if len(storedContoursToDelete) == 1:
        s = ''
        t = 's'
      else:
        s = 's'
        t = ''
      result = showYesNo('File used', 'Stored contour%s already use%s this fileName, and will be deleted: continue?' % (s, t), parent=self)
      if result:
        for storedContour in storedContoursToDelete:
          try:
            storedContour.delete()
          except:
            pass

    return result

  def contourAndSaveSpectrum(self):

    spectrum = self.spectrum
    if not spectrum:
      return

    if not self.spectrumConditions:
      return

    ###self.saveFile()
    fileName = self.file_entry.get()
    if not fileName:
      showError('No filename', 'No filename given', parent=self)
      return

    contourDir = self.getContourDir()

    fullPath = joinPath(contourDir, fileName)
    if not self.continueIfFileNameUsed(fullPath):
      return

    directory = os.path.dirname(fullPath)
    if not os.path.exists(directory):
      os.makedirs(directory)

    dims = self.dim_menu.getText()
    (xdim, ydim) = [int(dim) for dim in dims.split(', ')]

    ##self.saveLevels()
    ##levels = self.levels_entry.get()
    ##if not levels:
    ##  showError('No levels', 'No contour levels given', parent=self)
    ##  return
    analysisSpectrum = spectrum.analysisSpectrum
    levels = analysisSpectrum.posLevels + analysisSpectrum.negLevels
    scale = spectrum.scale / self.analysisProject.globalContourScale
    levels = [ level/scale for level in levels ]

    spectrumCondition = self.spectrumConditions[0]
    (condition, region) = spectrumCondition

    ndim = spectrum.numDim
    firstInt = ndim * [0]
    lastInt = ndim * [0]
    for i in range(ndim):
      try:
        (firstInt[i], lastInt[i]) = self.convertToPoints(spectrum, i, region[i])
      except Exception, e:
        showError('Invalid region', str(e), parent=self)

    try:
      #print 'about to saveSpectrumContours', fullPath, xdim, ydim, levels, firstInt, lastInt
      saveSpectrumContours(spectrum, fullPath, xdim, ydim, levels, firstInt, lastInt,
                         mem_cache=self.parent.mem_cache)
    except Exception, e:
      showError('Save error', str(e), parent=self)
      return

    createStoredContour(spectrum, fileName, xdim, ydim)

    showInfo('Created file', 'Successfully created stored contour file', parent=self)

  def convertToPoints(self, spectrum, n, region):

    (rMin, rMax) = region
    dim = n + 1
    dataDim = spectrum.findFirstDataDim(dim=dim)

    if dataDim.className != 'SampledDataDim':
      converter = UnitConverter.ppm2pnt
      dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
      (rMin, rMax) = (converter(rMax, dataDimRef), converter(rMin, dataDimRef))
      rMin = max(1, rMin)
      rMax = min(dataDim.numPoints, rMax)
    rMin = rMin - 1
    rMin = int(rMin+0.5)
    rMax = int(rMax+0.5)

    if rMin >= rMax:
      raise Exception('in dimension %d have invalid region' % dim)

    return (rMin, rMax)

  def setSpectrum(self, spectrum):

    self.expt_spectrum.set(spectrum, doCallback=True)


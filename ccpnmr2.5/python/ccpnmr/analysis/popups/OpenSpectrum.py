
"""
======================COPYRIGHT/LICENSE START==========================

OpenSpectrum.py: Part of the CcpNmr Analysis program

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

from memops.general import Implementation

from memops.universal import Io as uniIo

from memops.gui.ButtonList import UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.Frame import Frame
from memops.gui.FloatEntry import FloatEntry
from memops.gui.FileSelect import FileSelect, FileType
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showError, showWarning, showOkCancel
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix

from ccp.api.nmr import Nmr

from ccpnmr.analysis.core.AssignmentBasic import getShiftLists
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.popups.BrukerPseudo import BrukerPseudoPopup
from ccpnmr.analysis.popups.EditSpectrum import EditSpectrumPopup
from ccpnmr.analysis.core.ExperimentBasic import setAllIsotopeCodes, getRefExperiments, setRefExperiment, setExperimentShiftList, newShiftList, isSpectrumBigEndian, getIsSpectrumBigEndian, setIsSpectrumBigEndian
from ccpnmr.analysis.core.ExperimentBasic import isReferencingIncorrect
from ccpnmr.analysis.popups.NmrPipePseudo import NmrPipePseudoPopup
from ccpnmr.analysis.core.Util import defaultProject, getSoftware

# external data formats understood
from ccp.format.spectra.params.AzaraParams import AzaraParams
from ccp.format.spectra.params.BrukerParams import BrukerParams
from ccp.format.spectra.params.FelixParams import FelixParams
from ccp.format.spectra.params.NmrPipeParams import NmrPipeParams
from ccp.format.spectra.params.NmrViewParams import NmrViewParams
from ccp.format.spectra.params.UcsfParams import UcsfParams
from ccp.format.spectra.params.VarianParams import VarianParams
from ccp.format.spectra.params.FactorisedParams import FactorisedParams

# external data format information
file_formats = ('Azara', 'Bruker', 'Felix', 'NMRPipe', 'NMRView', 'UCSF', 
                'Varian', 'Factorised')

file_type_dict = { \
  'Azara':   FileType('Azara', ['*.par', '*.par.*']),
  'Bruker':  FileType('Bruker', ['*/*/*/procs']),
  'Felix':   FileType('Felix', ['*.mat']),
  'NMRPipe': FileType('NMRPipe', ['*.ft*']),
  'NmrView': FileType('NmrView', ['*.nv', '*.par']),
  'UCSF':    FileType('UCSF', ['*.ucsf']),
  'Varian':  FileType('Varian', ['procpar']),
  'Factorised':    FileType('Factorised', ['*.usf3', '*.xml']),
}

params_class_dict = { \
  'Azara':   AzaraParams,
  'Bruker':  BrukerParams,
  'Felix':   FelixParams,
  'NMRPipe': NmrPipeParams,
  'NMRView': NmrViewParams,
  'UCSF':    UcsfParams,
  'Varian':  VarianParams,
  'Factorised':    FactorisedParams,
}

details_file_dict = { \
  'Bruker':  'title',
  'Varian':  'text',
}

WINDOW_OPTS = ['All','First','None']

# 7 Jun 2012: RowObject introduced because ScrolledMatrix cannot have list
# objects and tuple objects cannot be modified, so use proper objects instead

class RowObject:

  def __init__(self, exptName, specName, fileName, window, shiftListName):

    self.exptName = exptName
    self.specName = specName
    self.fileName = fileName
    self.window = window
    self.shiftListName = shiftListName

class OpenSpectrumPopup(BasePopup):
  r"""
  **Locate Spectrum Data for Use in CCPN Project**
  
  This popup window enables the user to locate spectrum data within a file
  system and associate the files (typically binary) with an experiment and
  spectrum name so that it may be visualised and accessed within the current
  CCPN project. Spectra of many different origins and file formats may be
  loaded, which currently includes Bruker, Varian, Felix, NMRPipe, NmrView,
  SPARKY/UCSF, Azara and the factorised shape format "USF3". Depending upon the
  file format of the spectrum, data loaded the user may be required to either
  select a parameter file which then refers to the actual spectrum intensity
  data; this is true for Bruker "procs" and AZARA ".par" files, or alternatively
  a spectrum data file itself that contains referencing information; this is
  the case for SPARKY/UCSF, NmrView and NMRPipe spectra.

  The layout of the popup involved two sections; the upper of which is for
  navigating to and selecting the spectrum or parameter files within the
  file-system, and the lower is for specifying how each spectrum is loaded into
  the CCPN project. It should be noted that when spectrum parameters are read
  the first time, the relevant information is copied into the CCPN project,
  where it may be adjusted independently of the original file information. No
  copies of the spectrum intensity data are made, the CCPN project merely refers
  to the spectrum data on disk, although the data file for a loaded spectrum may
  subsequently be moved or replaced.

  In normal operation the user first selects the kind of spectrum file format
  that will be loaded via the upper "File format" pulldown menu and then checks
  that the "File type" pulldown (toward the bottom of the file browser) is set
  to detect the appropriate kinds of filename; if a helpful file name filter is
  not available the user can add one via the "Manual Select" field, taking care
  to add any wild-card symbols, like the asterisk in "\*.ft3". Next the spectrum
  data or parameter files, appropriate to the selected format, are located by
  navigating within the file-system browser. When the required spectrum files are
  visible the user selects one *or more* to load. Multiple file selections may
  be made using left-click with <Ctrl> (toggle selection) or <Shift> (select
  range). It should be noted that when selecting Bruker files, when using the
  standard Bruker directory structure, the user only needs to navigate to the
  numbered spectrum directory; by default the "procs" file two levels down is
  searched for, e.g. "123/pdata/1/procs" is shown in the directory containing
  the "123" directory.

  When spectrum or parameter files are selected in the file table, the lower
  "Spectra To Open" table is filled to reflect the selection. The user should
  then be mindful of the settings within this table and may choose to edit
  various things by double-clicking on the appropriate cell. Typically the user
  just adjusts the name of the independent "Experiment" and "Spectrum" records.
  These names are usually concatenated like "expName:specName" in CCPN graphical
  displays so there is no need to repeat a name in both fields; this only takes
  up more space. The Experiment, which is a record of *what was done
  experimentally*, commonly has a short name like "HNCA" or "HSQC_298K" so the
  user readily knows how to interpret the experimental data. The Spectrum, which
  is a record of *the data that was collected*, commonly has a short name to
  identify the spectrum number or file name. An Experiment record may contain
  several Spectrum records, so the spectrum's name need minimally only identify
  it amongst others from the same experiment. The Shift List value may be
  changed if the user knows that the experiment represents a distinct set of
  conditions, with different spectrum peak/resonance positions, to existing or
  other experiments being entered. Each shift list will be curated separately,
  to give separate chemical shift values for assignments made under different
  conditions (even when relating to the same atoms). The shift list that an
  experiment uses may also be changed at any time after loading.

  When all spectra and options are specified the [Open Spectrum] button will
  load the relevant data into the CCPN project. If the "Skip verification
  dialogs" option is set it is assumed that all of the spectrum point to
  frequency referencing information, and any data file references, are correct.
  Otherwise, the user will be prompted to confirm the file details and
  referencing information for each spectrum in turn. Finally, after loading the
  user is asked to set the type of NMR experiment, in terms of general
  magnetisation transfer pathway, that was performed.

  **Caveats & Tips**

  If the name of an Experiment that is *already within the CCPN project* is
  used, then the loaded spectrum will (assuming it is compatible) be entered
  under that existing experiment record; no new experiment entity will be
  defined. The user may legitimately use this feature to load several spectra
  that relate to the same experiment; typically where spectra are different
  projections. To facilitate this the "Use shared experiment" option can be
  selected.

  Although experiments and spectra may be renamed after loading, a spectrum
  record may not be placed under a different experiment once created; deletion
  and re-loading is the only mans of achieving this, and care must be taken in
  transferring any assignments.

  """
  
  def __init__(self, parent, *args, **kw):

    self.experiment     = None
    self.currentObject  = None
    #self.currentObjects = [] # looks obsolete
    
    BasePopup.__init__(self, parent=parent, title='Experiment : Open Spectra', **kw)

  def open(self):
  
    self.message()
    BasePopup.open(self)

  def body(self, guiFrame):
    
    self.fileSelect = None
    names, objects = self.getShiftLists()
    self.shiftListPulldown = PulldownList(self, callback=self.setShiftList,
                                          texts=names, objects=objects)
    self.windowPulldown    = PulldownList(self, texts=WINDOW_OPTS,
                                          callback=self.setWindow)
    self.experimentEntry = Entry(self, width=16, returnCallback=self.setExperiment)
    self.spectrumEntry   = Entry(self, width=16, returnCallback=self.setSpectrum)
 

    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(1, weight=1)
    
    leftFrame = LabelFrame(guiFrame, text='File Selection')
    leftFrame.grid(row=0, column=0, sticky='nsew')
    leftFrame.grid_columnconfigure(3, weight=1)

    row = 0

    label = Label(leftFrame, text='File format:')
    label.grid(row=row, column=0, sticky='w')
    tipText = 'Selects which kind of spectrum file is being loaded; what its data matrix format is'
    self.formatPulldown = PulldownList(leftFrame, callback=self.chooseFormat,
                                       texts=file_formats, tipText=tipText, grid=(row,1))
    
    self.detailsLabel = Label(leftFrame, text='Show details:')
    tipText = 'Whether to show an annotation that describes the spectrum in the file selection; currently only uses comment fields from Bruker spectra'
    self.detailsSelect = CheckButton(leftFrame, selected=False,
                                     callback=self.showDetails,
                                     tipText=tipText)
    self.titleRow = row
    self.detailsSelected = False

    row = row + 1
    leftFrame.grid_rowconfigure(row, weight=1)
    file_types = [ FileType('All', ['*']) ]
    self.fileSelect = FileSelect(leftFrame, multiSelect=True, file_types=file_types, single_callback=self.chooseFiles,
        extraHeadings=('Details',), extraJustifies=('left',), displayExtra=False, getExtraCell=self.getDetails,
        manualFileFilter=True)
    self.fileSelect.grid(row=row, column=0, columnspan=6, sticky='nsew')
    
    rightFrame = LabelFrame(guiFrame, text='Spectra To Open')
    rightFrame.grid(row=1, column=0, sticky='nsew')
    rightFrame.grid_columnconfigure(3, weight=1)

    row = 0
    label = Label(rightFrame, text='Skip verification dialogs:', grid=(row,0))
    tipText = 'Whether to allow the user to check file interpretation and referencing information before the spectrum is loaded'
    self.verifySelect = CheckButton(rightFrame, selected=False,
                                    grid=(row,1), tipText=tipText)
    
    label = Label(rightFrame, text='Use shared experiment:', grid=(row,2))
    tipText = 'When selecting multiple spectrum files, whether the loaded spectra will all belong to (derive from) the same experiment; useful for projection spectra etc.'
    self.sharedExpSelect = CheckButton(rightFrame, selected=False, tipText=tipText,
                                       callback=self.useShared, grid=(row,3))
    
    row = row + 1
    rightFrame.grid_rowconfigure(row, weight=1)
    tipTexts = ['A short textual name for the experiment record that the loaded spectrum will belong to; may be a new experiment or the name of an existing one',
                'A short textual name to identify the spectrum within its experiment; typically a few characters or spectrum number, rather than a repeat of the experiment name',
                'The location of the file, relative to the current directory, that the spectrum data will be loaded from',
                'Sets which window or windows the spectrum will initially appear within once loaded',
                'Sets which shift list the experiment (and hence loaded spectrum) will use to curate chemical shift information; can be changed after load time']
    headingList      = ['Experiment','Spectrum','File','Windows','Shift List']
    editWidgets      = [self.experimentEntry,self.spectrumEntry,None,self.windowPulldown, self.shiftListPulldown]
    editGetCallbacks = [self.getExperiment,  self.getSpectrum,  None,self.getWindow,      self.getShiftList]
    editSetCallbacks = [self.setExperiment,  self.setSpectrum,  None,self.setWindow,      self.setShiftList]
    self.scrolledMatrix = ScrolledMatrix(rightFrame, headingList=headingList,
                                         callback=self.selectCell,
                                         editWidgets=editWidgets, multiSelect=True,
                                         editGetCallbacks=editGetCallbacks,
                                         editSetCallbacks=editSetCallbacks,
                                         tipTexts=tipTexts,
                                         grid=(row,0), gridSpan=(1,4))

    row = row + 1
    tipTexts = ['Load spectrum or spectra into the CCPN project using the selected file(s)',]
    texts = ['Open Spectrum']
    commands = [self.openSpectra]
    bottomButtons = UtilityButtonList(guiFrame, texts=texts, tipTexts=tipTexts, doClone=False,
                                      commands=commands, helpUrl=self.help_url)
    bottomButtons.grid(row=row, column=0, columnspan=1, sticky='ew')
    self.openButton = bottomButtons.buttons[0]

    self.chooseFormat('Azara')
    self.message()

  def message(self):

    if not self.project or len(self.nmrProject.experiments) < 1:
      pass
      #self.parent.ticker.setMessage('Choose spectrum files to open....     ')

  def showDetails(self, isSelected):

    self.detailsSelected = isSelected
    self.fileSelect.updateDisplayExtra(isSelected)
    # below is so that when Details column is toggled on it will actually
    # be seen without having to use the scrollbar
    self.fileSelect.fileList.refreshSize()

  def useShared(self, isSelected):
    
    self.chooseFiles(forceUpdate=True)
    
    #if isSelected:
      
      #objects = self.scrolledMatrix.objectList
      #if len(objects) > 1:
      #  self.currentObject = objects[0]
      #  text = objects[0][0]
      #  self.chooseFiles()
      #  for oo in objects[1:]:
      #    oo[0] = text
      #  if self.project:
      #    self.experiment = self.nmrProject.findFirstExperiment(name=text)
      #self.update()

  def gridDetails(self, bool):

    if bool:
      self.detailsLabel.grid(row=self.titleRow, column=2, sticky='w')
      self.detailsSelect.grid(row=self.titleRow, column=3, sticky='w')
      self.fileSelect.updateDisplayExtra(self.detailsSelected)
    else:
      self.detailsLabel.grid_forget()
      self.detailsSelect.grid_forget()
      self.fileSelect.updateDisplayExtra(False)
    
  def openSpectra(self):
    
    noVerify = self.verifySelect.getSelected()
    
    # tracks if 'add to existing experiment' has already ben OK'ed
    self.okExpSet = set()
    
    directory = self.fileSelect.getDirectory()
    spectra = []
    specIndex = 0
    for obj in self.scrolledMatrix.objectList:
      fileName = uniIo.joinPath(directory, obj.fileName)
      spectrum = self.openSpectrum(obj.exptName, obj.specName, fileName, obj.window, 
                                   obj.shiftListName)
      specIndex += 1
      
      if (spectrum):
        # check endianness if we are not verifying
        spectra.append(spectrum)
        if noVerify:
          isBigEndian = isSpectrumBigEndian(spectrum) # according to data in file
          if isBigEndian is not None:
            isBigEndianCurr = getIsSpectrumBigEndian(spectrum) # according to data model
            setIsSpectrumBigEndian(spectrum, isBigEndian)
            if isBigEndian != isBigEndianCurr:
              if isBigEndian:
                s = 'big'
              else:
                s = 'little'
              print 'WARNING: swapped endianess of spectrum to %s endian' % s
    #
    del self.okExpSet
    
    
    if noVerify and len(spectra) > 1 and self.sharedExpSelect.getSelected():
      # if we are using a shared experiment and not verifying,
      # set referencing to match first spectrum for all
    
      # get reference spectrum and set up data structure
      # use most recent pre-existing spectrum, otherwise first new one
      refSpec = spectra[0]
      for spec in spectra[0].experiment.sortedDataSources():
        if spec in spectra:
          break
        else:
          refSpec = spec
      
      ddrLists = {}
      refDdrs = []
      for dataDim in refSpec.sortedDataDims():
        for ddr in dataDim.dataDimRefs:
          ddrLists[ddr.expDimRef] = []
          refDdrs.append(ddr)
      
      # get dataDimRefs, store by ExpDimRef, 
      # checking that all spectra have data dim refs for same set of xdr
      nTotal = len(ddrLists)
      for spec in spectra:
        nFound = 0
        for dataDim in spec.sortedDataDims():
          for ddr in dataDim.dataDimRefs:
            xdr = ddr.expDimRef
            ll = ddrLists.get(xdr)
            if ll is None:
              # something did not match - do nothing
              break
            else:
              ll.append(ddr)
              nFound += 1
        else:
          if nFound == nTotal:
            # we are OK. Do next spectrum
            continue
        # something did not match - do nothing
        break
      else:
 
        # all spectra matched. Now reset O1 references to match reference
        if refSpec is spectra[0]:
          startAt = 1
        else:
          startAt = 0
 
        for refDdr in refDdrs:
          dataDim = refDdr.dataDim
          centrePoint = dataDim.numPointsOrig/2 - dataDim.pointOffset + 1
          refValue = refDdr.pointToValue(centrePoint)
 
          xdr = refDdr.expDimRef
          for ddr in ddrLists[xdr][startAt:]:
            dataDim = ddr.dataDim
            centrePoint = dataDim.numPointsOrig/2 - dataDim.pointOffset + 1
            ddr.refPoint = centrePoint
            ddr.refValue = refValue
            
    # set refExperiment if there is only one possibility
    experiments = []
    ignoreSet = set()
    showPopup = False
    for spectrum in spectra:
      experiment = spectrum.experiment
      if experiment not in ignoreSet:
        ignoreSet.add(experiment)
        if not experiment.refExperiment:
          experiments.append(spectrum.experiment)
          if noVerify:
            resetCategory = False
            if not hasattr(experiment, 'category'):
              if (hasattr(experiment,'pulProgName') and
                  hasattr(experiment,'pulProgType')):
                # this is first time we get here, and we have external name and source 
                # use external source to set fullType
                experiment.category = 'use external'
                resetCategory = True
            refExperiments = getRefExperiments(experiment)
            if resetCategory and not refExperiments:
              # no refExperiments match external source. 
              # unset 'use external' category
              del experiment.category
              
            if len(refExperiments) == 1:
              # only one possibility, just set it
              setRefExperiment(experiment, refExperiments[0])
            
            # wb104: 20 Oct 2014: do not popup Experiment types dialog if noVerify
            #else:
            #  showPopup = True
    
    # Pop up refExperiment verification
    if experiments and (showPopup or not noVerify):
      self.parent.initRefExperiments(experiments)
    
    # set up internal Analysis data
    for spectrum in spectra:
      self.parent.finishInitSpectrum(spectrum)
      print 'finished opening spectrum', spectrum.experiment.name, spectrum.name

  def chooseFiles(self, forceUpdate=False, *file):
  
    directory = self.fileSelect.getDirectory()

    fileNames = self.fileSelect.fileList.currentObjects
    fullFileNames1 = [uniIo.joinPath(directory, x) for x in fileNames]

    fullFileNames2 = [x.fileName for x in self.scrolledMatrix.objectList]
    fullFileNames2 = [uniIo.joinPath(directory, x) for x in fullFileNames2]

    if fullFileNames1 == fullFileNames2 and not forceUpdate:
      return

    objectList  = []
    textMatrix  = []
    
    format = self.formatPulldown.getText()

    shiftListName = self.getShiftLists()[0][0]
    windowOpt = WINDOW_OPTS[1]
    oneUp = os.path.dirname
    
    if format == 'Bruker':
      if self.sharedExpSelect.getSelected():
        nameTemplate='Bruker_%d'
        next = self.getNextExpNum(nfiles=len(fileNames), 
                                   nameTemplate=nameTemplate)
        
        exptName = nameTemplate % (next)
        for i,fileName in enumerate(fileNames):
          fullFileName = fullFileNames1[i]
          specName = os.path.basename(oneUp(oneUp(oneUp(fullFileName))))
          datum = (exptName, specName, fileName, windowOpt, shiftListName)
          dataObj = RowObject(*datum)
          textMatrix.append(datum)
          objectList.append(dataObj)
      
      else:
        for i,fileName in enumerate(fileNames):
          fullFileName = fullFileNames1[i]
          try: # below should not fail
            ss1 = oneUp(fullFileName)
            specName = os.path.basename(ss1)
            ss2 = os.path.basename(oneUp(oneUp(ss1)))
            exptName = 'Bruker_' + ss2
          except: # just put in something
            ss = os.path.basename(fullFileName)
            exptName = 'Bruker_' + ss
            specName = ss
   
          datum = (exptName, specName, fileName, windowOpt, shiftListName)
          dataObj = RowObject(*datum)
          textMatrix.append(datum)
          objectList.append(dataObj)
    
    else:
      next = self.getNextExpNum(nfiles=len(fileNames))
      if self.sharedExpSelect.getSelected():
        exptName = 'Expt_%d' % (next)
        for i,fileName in enumerate(fileNames):
          specName = re.sub('\.\w+$','',fileName)

          datum = (exptName, specName, fileName, windowOpt, shiftListName)
          dataObj = RowObject(*datum)
          textMatrix.append(datum)
          objectList.append(dataObj)

      else:
        for i,fileName in enumerate(fileNames):
          exptName = 'Expt_%d' % (next+i)
          specName = re.sub('\.\w+$','',fileName)

          datum = (exptName, specName, fileName, windowOpt, shiftListName)
          dataObj = RowObject(*datum)
          textMatrix.append(datum)
          objectList.append(dataObj)

    if len(fileNames) > 1:
      self.openButton.config(text='Open Spectra')
    
    else:
      self.openButton.config(text='Open Spectrum')

    self.scrolledMatrix.update(objectList=objectList, textMatrix=textMatrix)
  
  def getNextExpNum(self, nfiles=0, nameTemplate = 'Expt_%d'):
    """ get suitable free integer to use for exp names 
    """
    next = 1
    if self.project:
      nmrProject = self.nmrProject
      ii = len(nmrProject.experiments)
 
      # find first exp number that is not taken
      # NBNB TBD could consider expname = specname, specname = proc dir
      next = ii + 1
      if nfiles:
        while ii < next + nfiles:
          ii += 1
          if nmrProject.findFirstExperiment(name=nameTemplate % ii):
            next = ii + 1
    #
    return next
    
  def getDetails(self, fullfile):

    details = ''
    if os.path.isfile(fullfile):
      format = self.formatPulldown.getText()
      detailsDir = os.path.dirname(fullfile)
      detailsFile = uniIo.joinPath(detailsDir, details_file_dict[format])
      if os.path.exists(detailsFile):
        fp = open(detailsFile)
        details = fp.read().strip().replace('\n', ' ').replace('\r', ' ')
        fp.close()

    return (details,)

  def update(self):
  
    objectList = self.scrolledMatrix.objectList
    textMatrix = [(obj.exptName, obj.specName, obj.fileName, obj.window, obj.shiftListName) for obj in objectList]
   
    self.scrolledMatrix.update(objectList=objectList, textMatrix=textMatrix)
 
  def selectCell(self, obj, row, col):
  
    self.currentObject = obj
    
    if self.project:
      self.experiment = self.nmrProject.findFirstExperiment(name=obj.exptName)
    else:
      self.experiment = None
    
  def getWindow(self, obj):
    
    if obj:  
      self.windowPulldown.set(obj.window)


  def setWindow(self, opt):

    if isinstance(opt, RowObject):
      self.currentObject.window = opt.window
    else:
      self.currentObject.window = opt
    self.update()  

 
  def setShiftList(self, obj=None):

    if self.project:
      project = self.project
      shiftList = self.shiftListPulldown.getObject()
      if shiftList is None:
        shiftList = newShiftList(project, unit='ppm')
      
      if self.experiment and shiftList and (shiftList is not self.experiment.shiftList):
        setExperimentShiftList(self.experiment, shiftList)
      
      self.currentObject.shiftListName = shiftList.name

    self.update()

  def getShiftList(self, object):

    names, shiftLists = self.getShiftLists()
    if names:
      self.shiftListPulldown.setup(names, shiftLists, 0)
      if self.experiment and self.experiment.shiftList:
        name = self.experiment.shiftList.name
      else:
        name = object.shiftListName
        
      if name is not None:    
        self.shiftListPulldown.set(name)
  
  def getShiftLists(self):
  
    if self.project:
      names = []
      objects = getShiftLists(self.nmrProject)
      for shiftList in objects:
        if not shiftList.name:
          shiftList.name = 'ShiftList %d' % shiftList.serial
        names.append( shiftList.name )
      
      objects.append(None)
      names.append('<New>')
      
    else:
      objects = [None,]
      names = ['ShiftList 1',]
     
    return names, objects

  def chooseFormat(self, format):

    if format in ('Bruker', 'Varian'):
      self.gridDetails(True)
    else:
      self.gridDetails(False)

    file_types = []
    file_type = file_type_dict.get(format)
    if (file_type):
      file_types.extend([ file_type ])
    file_types.append(FileType('All', ['*']))
    file_types.append(self.fileSelect.manualFilter)
    self.fileSelect.setFileTypes(file_types)

  def getSpectrum(self, obj):

    if obj:
      self.spectrumEntry.set(obj.specName)

  def setSpectrum(self, *event):

    if self.currentObject:
      text = self.spectrumEntry.get()
      if text and text != ' ':
        for data in self.scrolledMatrix.objectList:
          if data is self.currentObject:
            continue
          if (data.specName == text) and (data.exptName == self.currentObject.exptName):
            showWarning('Repeated name','Spectrum name (%s) already in use for experiment (%s)' % (data.specName,data.exptName), parent=self)
            return
          elif (self.experiment) and (self.experiment.findFirstDataSource(name=text)):
            showWarning('Repeated name','Spectrum name (%s) already in use for experiment (%s)' % (data.specName,data.exptName), parent=self)
            return
        
        self.currentObject.specName = text
      self.update()

  def getExperiment(self, obj):

    if obj:
      self.experimentEntry.set(obj.exptName)

  def setExperiment(self, *event):

    if self.currentObject:
      text = self.experimentEntry.get()
      if text and text != ' ':
        if self.sharedExpSelect.getSelected():
          # share one experiment for all rows
          for oo in self.scrolledMatrix.objectList:
            oo.exptName = text
        else:
          #separate experiments
          self.currentObject.exptName = text
        if self.project:
          self.experiment = self.nmrProject.findFirstExperiment(name=text)
      self.update()

  def updateShiftLists(self):
  
    if self.project:
      name = self.expt_entry.get()
      e = self.nmrProject.findFirstExperiment(name=name)
    
    else:
      e = None

    names, objects = self.getShiftLists()
    if e and e.shiftList:
      index = objects.index(e.shiftList)
    else:
      index = 0
       
    self.shiftListPulldown.setup(names, objects, index)

  def openSpectrum(self, exptName, specName, file, windowOpt=WINDOW_OPTS[2], 
                   shiftListName='<New>', extraData=None):
    
    # input checks
    if not file:
      showError('No file', 'Need to enter file', parent=self)
      return None

    if not exptName:
      showError('Experiment', 'Need to enter experiment name', parent=self)
      return None

    if not specName:
      showError('Spectrum', 'Need to enter spectrum name', parent=self)
      return None

    # get or set up project
    project = self.project
    if not project:
      self.project = project = defaultProject()
      self.parent.initProject(project)
      self.nmrProject = self.parent.nmrProject
      self.analysisProject = self.parent.analysisProject
      #Default ShiftList with name 'ShiftList 1' created
    
    # set up shift list
    if shiftListName =='<New>':
      shiftList = None
    else:
      shiftList = self.nmrProject.findFirstMeasurementList(className='ShiftList',
                                                           name=shiftListName)

    # read params

    format = self.formatPulldown.getText()
    clazz = params_class_dict[format]
    try:
      params = clazz(file, extraData=extraData)
    except Implementation.ApiError, e:
      showError('Reading params file', 'Fatal error: ' + e.error_msg, 
                parent=self)
      return None

    dim = params.pseudoDataDim()
    if dim is not None:
      if format == 'NMRPipe':
        popup = NmrPipePseudoPopup(self, params, dim, file)
        popup.destroy()
      elif format == 'Bruker':
        popup = BrukerPseudoPopup(self, params, dim)
        popup.destroy()

    # get or set up experiment

    experiment = self.nmrProject.findFirstExperiment(name=exptName)
    if experiment:
      expIsNew = False
      if experiment.findFirstDataSource(name=specName):
        showError('Duplicate name', 
                  'Duplicate spectrum name "%s" in experiment %s' 
                  % (specName, experiment.name), parent=self)
        return None
      elif (experiment.dataSources and experiment not in self.okExpSet):
        if showOkCancel('Multiple Spectra Warning', 'Really put multiple '
                            'spectra into existing experiment %s?' 
                            % experiment.name, parent=self):
          self.okExpSet.add(experiment)
        else:
          return
      
    else:
      expIsNew = True
      try:
        # Will also work for shiftList == None
        experiment = Nmr.Experiment(self.nmrProject, name=exptName, 
                                    numDim=params.ndim, shiftList=shiftList)
          
      except Implementation.ApiError, experiment:
        showError('Experiment', experiment.error_msg, parent=self)
        return None
    
    
    # set shiftList
    if shiftListName == '<New>':
      shiftList = newShiftList(project)
      
    if shiftList and  experiment.shiftList is not shiftList:
        setExperimentShiftList(experiment, shiftList)
    
    # create spectrum
    try:
      spectrum = params.createDataSource(experiment, specName)
      
    except Implementation.ApiError, exc:
      showError('Spectrum', exc.error_msg, parent=self)
      raise
      return None
    
    
    if windowOpt == WINDOW_OPTS[1]: # Put in first window only
      self.parent.visibleSpectra[spectrum] = 1
      
    elif windowOpt == WINDOW_OPTS[2]: # Do not put in any windows
      self.parent.visibleSpectra[spectrum] = False

    #self.close()
    
    
    if format == 'Bruker' and expIsNew:
      # Check referencing and possibly reset.
      if self.verifySelect.getSelected():
        incorrect, text = isReferencingIncorrect(spectrum, fixErrors=True)
        if incorrect:
          print """\nWARNING: Referencing was reset
         The following problem(s) were found and fixed:
%s\n""" % text
      
      else:
        incorrect, text = isReferencingIncorrect(spectrum)
        if (incorrect
            and showOkCancel('Reset Referencing?', text, parent=self)):
          isReferencingIncorrect(spectrum, fixErrors=True, printChecks=False)

    if self.verifySelect.getSelected():
      # skipping verification - ensure isotopeCodes are set somehow
      setAllIsotopeCodes(experiment)
      
    else:
      # do verification
      
      useReducedDim = False
      if extraData:
        if extraData.get('displayNames') or extraData.get('scalingFactors'):
          useReducedDim = True
    
      popup = EditSpectrumPopup(self.parent, transient=False, modal=True, 
                                useReducedDim=useReducedDim, spectrum=spectrum)
      cancelled = popup.cancelled
      popup.destroy()
      
      if cancelled:
        spectrum.delete()
        if expIsNew:
          if experiment.findFirstDataSource() is None:
            # should be unnecessary to ask, but no harm
            experiment.delete()
        return None
    
    # handle extra parameters:
    if params.format == 'Factorised':
      # factorised (ShapeMatrix) data - read peaks
      software = getSoftware(spectrum.root)
      params.parsePeaks(spectrum, software=software)
    
    pulProgName = params.pulProgName
    if pulProgName is not None:
      if experiment.refExperiment is None:
        if not (hasattr(experiment, 'pulProgName') 
                and experiment.pulProgName is not None):
          experiment.pulProgName = pulProgName
          experiment.pulProgType = params.pulProgType
    
    print 'Spectrum successfully opened'

    return spectrum

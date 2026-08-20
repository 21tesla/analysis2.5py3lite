
"""
======================COPYRIGHT/LICENSE START==========================

BrukerPseudo.py: Part of the CcpNmr Analysis program

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
import re

from memops.general.Implementation import ApiError

from memops.gui.Button import Button
from memops.gui.ButtonList import UtilityButtonList
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError
#from memops.gui.MultiWidget import MultiWidget
from memops.gui.PulldownList import PulldownList
from memops.gui.RadioButtons import RadioButtons

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.ExperimentBasic import getSampledDimExperiments, getExperimentSampledDim

PSEUDO_ENTRIES = ('Is Pseudo-%dD Expt', 'Is Not Pseudo-%dD Expt')

class BrukerPseudoPopup(BasePopup):

  def __init__(self, parent, params, dim, *args, **kw):

    self.dim = dim
    self.params = params

    BasePopup.__init__(self, parent=parent, title='Bruker Pseudo Data', modal=True, **kw)

  def body(self, master):

    pseudoExpts = getSampledDimExperiments(self.parent.nmrProject)

    master.rowconfigure(0, weight=1)
    master.rowconfigure(1, weight=1)
    master.columnconfigure(0, weight=1)

    tipTexts = ['The experiment is pseudo-N dimensional, with a sampled axis',
                'The experiment is the regular kind with only NMR frequency axes']
    self.pseudoEntries = [x % len(self.params.npts) for x in PSEUDO_ENTRIES]            
    self.pseudoButton = RadioButtons(master, entries=self.pseudoEntries,
                                     select_callback=self.changedPseudoMode,
                                     grid=(0,0), sticky='nw', tipTexts=tipTexts)

    frame = self.pseudoFrame = Frame(master)
    self.pseudoFrame.grid(row=1, column=0, sticky='nsew')

    row = 0
    if pseudoExpts:
      tipText = 'Select from existing pseudo nD experiments to copy sampled axis values from'
      texts = [x.name for x in pseudoExpts]
      label = Label(frame, text='Existing pseudo expts: ')
      label.grid(row=row, column=0, sticky='e')
      self.pseudoList = PulldownList(frame, texts=texts, objects=pseudoExpts, tipText=tipText)
      self.pseudoList.grid(row=row, column=1, sticky='w')
      tipText = 'Transfer the sampled axis values from the existing experiment to the new one'
      Button(frame, text='Copy values down', command=self.copyValues,
             tipText=tipText, grid=(row,2))
      
      row += 1

    npts = self.params.npts[self.dim]
    tipText = 'Number of data points (planes) along sampled axis'
    label = Label(frame, text='Number of points: ')
    label.grid(row=row, column=0, sticky='e')
    self.nptsEntry = IntEntry(frame, text=npts, tipText=tipText, width=8, grid=(row,1))
     
    tipText = 'Load the values for the sampled axis from a text file containing a list of numeric values'
    Button(frame, text='Load File', command=self.loadValues,
           tipText=tipText, grid=(row,2), sticky='ew')
     
    row += 1

    tipText = 'The values (e.g. T1, T2) corresponding to each data point (plane) along sampled axis'
    label = Label(frame, text='Point values: ')
    label.grid(row=row, column=0, sticky='e')
    self.valueEntry = FloatEntry(frame, isArray=True, tipText=tipText)
    #minRows = self.params.npts[self.dim]
    #self.valueEntry = MultiWidget(frame, FloatEntry, callback=None, minRows=minRows, maxRows=None,
    #                              options=None, values=[], useImages=False)
    self.valueEntry.grid(row=row, column=1, columnspan=2, sticky='ew')
    row += 1

    label = Label(frame, text='(requires comma-separated list, of length number of points)')
    label.grid(row=row, column=1, columnspan=2, sticky='w')
    row += 1

    for n in range(row):
      frame.rowconfigure(n, weight=1)
    frame.columnconfigure(1, weight=1)

    buttons = UtilityButtonList(master, closeText='Ok', closeCmd=self.updateParams,
                                helpUrl=self.help_url)
    buttons.grid(row=row, column=0, sticky='ew')

  def loadValues(self):
  
    directory = self.parent.fileSelect.getDirectory()
    fileSelectPopup = FileSelectPopup(self, title='Select Sampled Data File',
                                      dismiss_text='Cancel',
                                      selected_file_must_exist=True,
                                      multiSelect=False,
                                      directory=directory)
 

    fileName = fileSelectPopup.file_select.getFile()
    
    fileObj = open(fileName, 'rU')
    
    data = ''
    line = fileObj.readline()
    while line:
      data += line
      line = fileObj.readline()
    
    data = re.sub(',\s+', ',', data)
    data = re.sub('\s+', ',', data)
    data = re.sub(',,', ',', data)
    data = re.sub('[^0-9,.\-+eE]', '', data)
      
    self.valueEntry.set(data)

  def copyValues(self):

    expt = self.pseudoList.getObject()
    if expt:
      dataDim = getExperimentSampledDim(expt)
      values = dataDim.pointValues
      self.nptsEntry.set(len(values))
      self.valueEntry.set(values)

  def updateParams(self):

    params = self.params
    if self.pseudoButton.get() == self.pseudoEntries[0]:
      npts = self.nptsEntry.get()
      params.npts[self.dim] = npts
      values = self.valueEntry.get()
      try:
        params.setSampledDim(self.dim, values)
      except ApiError, e:
        showError('Set Sampled Dim', e.error_msg, parent=self)
        return
    else:
      params.fixNullDims()

    self.close()

  def changedPseudoMode(self, option):

    if option == self.pseudoEntries[0]:
      self.pseudoFrame.grid(row=1, column=0, sticky='nsew')
    else:
      self.pseudoFrame.grid_forget()


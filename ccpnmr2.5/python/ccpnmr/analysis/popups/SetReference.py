
"""
======================COPYRIGHT/LICENSE START==========================

SetReference.py: Part of the CcpNmr Analysis program

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
import Tkinter

from memops.universal.Util import formatFloat

from memops.gui.ButtonList import UtilityButtonList
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError
from memops.gui.ScrolledMatrix import ScrolledMatrix

from ccp.api.nmr import Nmr
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.ExperimentBasic import isSpectrum, getPrimaryDataDimRef
from ccpnmr.analysis.core.Util import convertPosition

yes_string = 'yes'
no_string = 'no'

class SetReferencePopup(BasePopup):

  def __init__(self, parent, dataDim, position):

    self.dataDim  = dataDim
    self.position = position

    BasePopup.__init__(self, parent=parent, title='Set spectrum referencing',
                       modal=True, transient=True)

  def body(self, guiFrame):

    dataDim = self.dataDim
    dataDimRef = getPrimaryDataDimRef(dataDim)
    spectrum = dataDim.dataSource

    guiFrame.grid_columnconfigure(0, weight=1)

    row = 0
    label = Label(guiFrame, text='Experiment:', grid=(row,0), sticky='e')
    tipText = 'The name of the experiment that contains the spectrum to be re-referenced'
    label = Label(guiFrame, text=spectrum.experiment.name,
                  grid=(row,1), tipText=tipText)

    row += 1
    label = Label(guiFrame, text='Spectrum:', grid=(row,0), sticky='e')
    tipText = 'The name of the spectrum to be re-referenced'
    label = Label(guiFrame, text=spectrum.name,
                  grid=(row,1), tipText=tipText)

    row += 1
    label = Label(guiFrame, text='Dimension:', grid=(row,0), sticky='e')
    tipText = 'The number of the spectrum dimension that will be re-referenced'
    label = Label(guiFrame, text=self.dataDimText(dataDim),
                  grid=(row,1), tipText=tipText)

    row += 1
    label = Label(guiFrame, text='Chosen reference point:', grid=(row,0), sticky='e')
    text = formatFloat(self.position, places=6)
    tipText = 'The chosen position in the spectrum dimension (where the mouse was clicked) in data matrix points'
    label = Label(guiFrame, text=text,
                  grid=(row,1), tipText=tipText)

    row += 1
    label = Label(guiFrame, text='Current reference ppm at this point:', grid=(row,0), sticky='e')
    self.current_refppm = convertPosition(self.position, dataDimRef, toUnit='ppm')
    text = formatFloat(self.current_refppm, places=6)
    tipText = 'The ppm value currently associated with this points location'
    label = Label(guiFrame, text=text,
                  grid=(row,1), tipText=tipText)

    row += 1
    label = Label(guiFrame, text='New reference ppm at this point:', grid=(row,0), sticky='e')
    tipText = 'Inputs a new value for the exact ppm position of the selected point; consequently re-referencing the spectrum dimension'
    self.refppm_entry = FloatEntry(guiFrame, width=30, grid=(row,1),
                                  sticky='ew', tipText=tipText)
    self.refppm_entry.set(0.0) # default

    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    tipTexts = ['The name of experiments within the project containing similar spectrum dimensions',
                'The name of the spectrum within the project that may be re-referenced',
                'The spectrum dimension to which the re-referencing may apply',
                'The current reference point for the spectrum dimension in its spectrum data matrix',
                'The current reference ppm value for the spectrum dimension; the value at the reference point',
                'Sets whether the spectrum dimension will be re-referenced using the above reference point and ppm value']
    headings = ('Experiment', 'Spectrum', 'Dimension', 'Reference\npoint',
                'Reference\nppm', 'Change\nreferencing')
    editWidgets = editSetCallbacks = 6 * [ None ]
    editGetCallbacks = [ None, None, None, None, None, self.toggleChangeRef ]
    self.datadim_list = ScrolledMatrix(guiFrame, headingList=headings,
                          initialRows=4, editWidgets=editWidgets,
                          editGetCallbacks=editGetCallbacks,
                          editSetCallbacks=editSetCallbacks,
                          grid=(row,0), gridSpan=(1,2), tipTexts=tipTexts)
    tipTexts = ['Re-reference the selected spectrum dimensions using the stated ppm value',]
    row += 1
    texts = [ 'Commit' ]
    commands = [ self.ok ]
    buttons = UtilityButtonList(guiFrame, texts=texts, commands=commands, doClone=False,
                                closeText='Cancel', helpUrl=self.help_url,
                                grid=(row,0), gridSpan=(1,2), tipTexts=tipTexts)

    self.updateDataDimList()

  def dataDimText(self, dataDim):

    dataDimRef = getPrimaryDataDimRef(dataDim)
    
    if dataDimRef:
      isotope = '/'.join(dataDimRef.expDimRef.isotopeCodes)
    else:
      isotope = 'None'
      
    text = str(dataDim.dim) + ' (%s)' % isotope

    return text

  def updateDataDimList(self, *extra):

    textMatrix = []
    self.dataDims = []
    isotopeCodes = getPrimaryDataDimRef(self.dataDim).expDimRef.isotopeCodes
    project = self.dataDim.root
    for experiment in self.nmrProject.experiments:
      for spectrum in experiment.dataSources:
        if (isSpectrum(spectrum)):
          for dataDim in spectrum.sortedDataDims():
            if (isinstance(dataDim, Nmr.FreqDataDim)):
              expDimRef = getPrimaryDataDimRef(dataDim).expDimRef
              if (expDimRef.isotopeCodes == isotopeCodes):
                if (not hasattr(dataDim, 'change_ref')):
                  if (dataDim == self.dataDim):
                    dataDim.change_ref = True
                  else:
                    dataDim.change_ref = False
                text = []
                text.append(experiment.name)
                text.append(spectrum.name)
                text.append(self.dataDimText(dataDim))
                dataDimRef = getPrimaryDataDimRef(dataDim)
                text.append(formatFloat(dataDimRef.refPoint, places=6))
                text.append(formatFloat(dataDimRef.refValue, places=6))
                if (dataDim.change_ref):
                  text.append(yes_string)
                else:
                  text.append(no_string)
                textMatrix.append(text)
                self.dataDims.append(dataDim)

    self.datadim_list.update(objectList=self.dataDims, textMatrix=textMatrix)

  def toggleChangeRef(self, dataDim):

    dataDim.change_ref = not dataDim.change_ref
    self.updateDataDimList()

  def apply(self):

    refppm = self.refppm_entry.get()

    if (refppm is None):
      showError('No reference', 'Must specify reference ppm.', parent=self)
      return False

    delta_refppm = refppm - self.current_refppm

    for dataDim in self.dataDims:
      if (dataDim.change_ref):
        dataDimRef = getPrimaryDataDimRef(dataDim)
        dataDimRef.refValue = dataDimRef.refValue + delta_refppm
      del dataDim.change_ref

    return True

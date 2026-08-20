
"""
======================COPYRIGHT/LICENSE START==========================

ChangeAxisMapping.py: Part of the CcpNmr Analysis program

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
import string

import Tkinter


from ccpnmr.api import Analysis

from memops.gui.Label import Label
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ButtonList import UtilityButtonList

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.Util import haveTypeMatch
from ccpnmr.analysis.core.WindowBasic import findDataDimAxisMapping, getWindowPaneName

class ChangeAxisMappingPopup(BasePopup):

  def __init__(self, parent, spectrumWindowView, **kw):

    self.spectrumWindowView = spectrumWindowView

    BasePopup.__init__(self, parent=parent, title='Change axis mappings',
                       modal=True, transient=False, **kw)

  def body(self, master):

    master.grid_columnconfigure(1, weight=1)
    master.grid_rowconfigure(0, weight=1)

    view = self.spectrumWindowView
    spectrum = view.analysisSpectrum.dataSource
    dataDims = spectrum.sortedDataDims()
    
    self.ndim = spectrum.numDim
    self.prev_label = self.ndim * [ None ]
    self.redraw = False

    row = 0
    label = Label(master, text='Window:')
    label.grid(row=row, column=0, sticky='w')
    
    name = getWindowPaneName(self.spectrumWindowView.spectrumWindowPane)
    tipText = 'The name of the window where the dimension mapping is changing '
    label = Label(master, text=name, tipText=tipText)
    label.grid(row=row, column=1, sticky='w')

    row += 1
    label = Label(master, text='Experiment:')
    label.grid(row=row, column=0, sticky='w')
    tipText = 'The NMR experiment to change dimension mapping for'
    label = Label(master, text=spectrum.experiment.name, tipText=tipText)
    label.grid(row=row, column=1, sticky='w')

    row += 1
    label = Label(master, text='Spectrum:')
    label.grid(row=row, column=0, sticky='w')
    tipText = 'The spectrum data to change the dimension mapping for'
    label = Label(master, text=spectrum.name, tipText=tipText)
    label.grid(row=row, column=1, sticky='w')

    self.axis_lists = self.ndim * [0]
    for n in range(self.ndim):
      row += 1
      callback = lambda entry_index, entry, dim=n: self.callback(dim, entry)
      dataDim = dataDims[n]
      axisMapping = findDataDimAxisMapping(view, dataDim)
      entries = self.validAxisLabels(dataDim)
      selected = axisMapping.label
      if selected not in entries:
        label = Label(master, text='Dimension %d:\n(sampled)' % (n+1,))
        label.grid(row=row, column=0, sticky='w')
      else:
        label = Label(master, text='Dimension %d:' % (n+1,))
        label.grid(row=row, column=0, sticky='w')
        selected_index = entries.index(selected)
        tipText = 'Selects which numbered spectrum dimension goes with which window axis'
        self.axis_lists[n] = PulldownMenu(master, callback=callback,
                                          entries=entries, tipText=tipText,
                                          selected_index=selected_index)
        self.axis_lists[n].grid(row=row, column=1, sticky='w')

    row += 1
    texts = commands = []
    tipTexts = ['Commit any changes to the spectrum dimension to window axis mapping and close the popup',]
    texts = [ ' Commit ',]
    commands = [ self.ok,]
    buttons = UtilityButtonList(master, texts=texts, tipTexts=tipTexts,
                                commands=commands, doClone=False,
                                helpUrl=self.help_url)
    buttons.grid(row=row, column=0, columnspan=2, sticky='ew')
    
  def validAxisLabels(self, dataDim):

    axisPanels = self.spectrumWindowView.spectrumWindowPane.sortedAxisPanels()
    labels = [ axisPanel.label for axisPanel in axisPanels if haveTypeMatch(axisPanel, dataDim) ]

    return labels

  def callback(self, n, entry):

    prev_entry = self.prev_label[n]
    self.prev_label[n] = entry

    if (not prev_entry):
      return

    if (entry == prev_entry):
      return

    for m in range(self.ndim):
      if (m == n):
        continue
      if (entry == self.axis_lists[m].getSelected()):
        self.prev_label[m] = prev_entry
        self.axis_lists[m].setSelected(prev_entry)
        break

  def apply(self):

    view = self.spectrumWindowView
    for axisMapping in view.axisMappings:
      axisMapping.delete()

    spectrum = view.analysisSpectrum.dataSource
    dataDims = spectrum.sortedDataDims()
    for n in range(self.ndim):
      dataDim = dataDims[n]
      label = self.axis_lists[n].getSelected()
      view.newAxisMapping(label=label, analysisDataDim=dataDim.analysisDataDim)

    self.redraw = True

    return True

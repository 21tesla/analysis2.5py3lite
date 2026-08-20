
"""
======================COPYRIGHT/LICENSE START==========================

ExptSpectrumRows.py: Part of the CcpNmr Analysis program

Copyright (C) 2005 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in ../../../license/CCPN.license.

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

from memops.general import Implementation

from memops.gui.Label import Label
from memops.gui.Frame import Frame

from ccpnmr.analysis.core.ExperimentBasic import getExperimentSpectra
from ccpnmr.analysis.frames.ExperimentList import ExperimentList
from ccpnmr.analysis.frames.SpectrumList import SpectrumList

notify_funcs = ('__init__', 'delete', 'setName')

class ExptSpectrumRows(Frame):

  def __init__(self, parent, analysis, orientation='horizontal', callback = None, *args, **kw):

    self.analysis = analysis
    self.callback = callback

    apply(Frame.__init__, (self, parent) + args, kw)

    label = Label(self, text='Experiment:')
    label.grid(row=0, column=0, sticky='ne')
    self.expt_list = ExperimentList(self, self.analysis.getExperiments,
                                    callback=self.setSpectra)
    self.expt_list.grid(row=0, column=1, sticky='nw')
 
    if orientation in ['horizontal','h','H',Tkinter.HORIZONTAL]:
      row = 0
      col = 2
    else:
      row = 1
      col = 0
    
    label = Label(self, text='Spectrum:')
    label.grid(row=row, column=col, sticky='ne')
    self.spectrum_list = SpectrumList(self, self.getSpectra,
                                      callback=self.setSpectrumProperties)
    self.spectrum_list.grid(row=row, column=col+1, sticky='nw')

    for func in notify_funcs:
      Implementation.registerNotify(self.setExperiments, 'ccp.nmr.Nmr.Experiment', func)
      Implementation.registerNotify(self.setSpectra, 'ccp.nmr.Nmr.DataSource', func)

  def destroy(self):

    self.spectrum_list.destroy()
    self.expt_list.destroy()

    for func in notify_funcs:
      Implementation.unregisterNotify(self.setExperiments, 'ccp.nmr.Nmr.Experiment', func)
      Implementation.unregisterNotify(self.setSpectra, 'ccp.nmr.Nmr.DataSource', func)

  def setSpectra(self, *extra):
 
    self.spectrum_list.setSpectra()

  def setSpectrumProperties(self, *extra):

    spectrum = self.getSpectrum()

    if (self.callback):
      self.callback(spectrum)
 
  def setExperiments(self, *extra):

    self.expt_list.setExperiments()

  def getExperiment(self):
 
    ind = self.expt_list.getSelectedIndex()
    if (ind >= 0):
      project = self.analysis.getProject()
      try:
        experiment = project.currentNmrProject.sortedExperiments()[ind]
      except IndexError: # if no experiments for some reason ind = 0 rather than -1
        experiment = None
    else:
      experiment = None
 
    return experiment
 
  def getSpectra(self):
 
    experiment = self.getExperiment()
    if (experiment):
      return getExperimentSpectra(experiment)
    else:
      return []
 
  def getSpectrum(self):

    spectra = self.getSpectra()
    ind = self.spectrum_list.getSelectedIndex()
    if (ind >= 0 and (ind < len(spectra))):
      spectrum = spectra[ind]
    else:
      spectrum = None
 
    return spectrum

  def setSpectrum(self, spectrum):

    if (not spectrum):
      return

    experiments = list(self.analysis.getExperiments())
    expt_index = experiments.index(spectrum.experiment)
    spectra = getExperimentSpectra(spectrum.experiment)
    spectrum_index = spectra.index(spectrum)
    
    self.expt_list.setSelectedIndex(expt_index)
    self.spectrum_list.setSelectedIndex(spectrum_index)

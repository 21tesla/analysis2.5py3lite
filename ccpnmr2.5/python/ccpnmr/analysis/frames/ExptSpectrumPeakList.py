
"""
======================COPYRIGHT/LICENSE START==========================

ExptSpectrumPeakList.py: Part of the CcpNmr Analysis program

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
from memops.general import Implementation

from memops.gui.Label import Label
from memops.gui.Frame import Frame

from ccpnmr.analysis.core.ExperimentBasic import getExperimentSpectra
from ccpnmr.analysis.frames.ExperimentList import ExperimentList
from ccpnmr.analysis.frames.PeakListList import PeakListList
from ccpnmr.analysis.frames.SpectrumList import SpectrumList

notify_funcs = ('__init__', 'delete', 'setName')

class ExptSpectrumPeakList(Frame):

  def __init__(self, parent, analysis, callback = None, *args, **kw):

    self.analysis = analysis
    self.callback = callback

    apply(Frame.__init__, (self, parent) + args, kw)

    label = Label(self, text='Experiment:')
    label.grid(row=0, column=0, sticky='ne')
    self.expt_list = ExperimentList(self, self.analysis.getExperiments,
                                    callback=self.setSpectra)
    self.expt_list.grid(row=0, column=1, sticky='nw')
 
    label = Label(self, text='Spectrum:')
    label.grid(row=0, column=2, sticky='ne')
    self.spectrum_list = SpectrumList(self, self.getSpectra,
                                      callback=self.setPeakLists)
    self.spectrum_list.grid(row=0, column=3, sticky='nw')
 
    label = Label(self, text='PeakList:')
    label.grid(row=0, column=4, sticky='ne')
    self.peak_list_list = PeakListList(self, self.getPeakLists,
                                      callback=self.doCallback)
    self.peak_list_list.grid(row=0, column=5, sticky='nw')

    for func in notify_funcs:
      Implementation.registerNotify(self.setExperiments, 'ccp.nmr.Nmr.Experiment', func)
      Implementation.registerNotify(self.setSpectra, 'ccp.nmr.Nmr.DataSource', func)
      if func != 'setName':
        Implementation.registerNotify(self.setPeakLists, 'ccp.nmr.Nmr.PeakList', func)

  def destroy(self):

    self.peak_list_list.destroy()
    self.spectrum_list.destroy()
    self.expt_list.destroy()

    for func in notify_funcs:
      Implementation.unregisterNotify(self.setExperiments, 'ccp.nmr.Nmr.Experiment', func)
      Implementation.unregisterNotify(self.setSpectra, 'ccp.nmr.Nmr.DataSource', func)
      if func != 'setName':
        Implementation.unregisterNotify(self.setPeakLists, 'ccp.nmr.Nmr.PeakList', func)

  def doCallback(self, *extra):

    peakList = self.getPeakList()

    if (self.callback):
      self.callback(peakList)
 
  def setExperiments(self, *extra):

    self.expt_list.setExperiments()

  def setSpectra(self, *extra):
 
    self.spectrum_list.setSpectra()

  def setPeakLists(self, *extra):

    self.peak_list_list.setPeakLists()

  def getExperiment(self):
 
    ind = self.expt_list.getSelectedIndex()
    if (ind >= 0):
      project = self.analysis.getProject()
      try:
        experiment = project.currentNmrProject.experiments[ind]
      except: # if no experiments for some reason ind = 0 rather than -1
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
 
  def getPeakLists(self):

    spectrum = self.getSpectrum()
    if (spectrum):
      return spectrum.peakLists
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

  def getPeakList(self):

    peakLists = self.getPeakLists()
    ind = self.peak_list_list.getSelectedIndex()
    if (ind >= 0 and (ind < len(peakLists))):
      peakList = peakLists[ind]
    else:
      peakList = None
 
    return peakList

  def setPeakList(self, peakList):

    if (not peakList):
      return

    spectrum = peakList.spectrum
    experiments = list(self.analysis.getExperiments())
    expt_index = experiments.index(spectrum.experiment)
    spectra = getExperimentSpectra(spectrum.experiment)
    spectrum_index = spectra.index(spectrum)

    self.expt_list.setSelectedIndex(expt_index)
    self.spectrum_list.setSelectedIndex(spectrum_index)
    self.peak_list_list.setSelected(peakList.serial)

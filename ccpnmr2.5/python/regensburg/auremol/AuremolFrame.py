import Tkinter, os

from memops.gui.Frame import Frame
from memops.gui.FloatEntry import FloatEntry
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.ButtonList import ButtonList
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError
from memops.gui.FloatEntry import FloatEntry
from memops.gui.IntEntry import IntEntry
from memops.gui.PulldownList import PulldownList
from memops.gui.LabelDivider import LabelDivider
from memops.gui.BasePopup import BasePopup

from memops.universal.Io import getTopDirectory

from regensburg.auremol.findAuremolPeaks import findAuremolPeaksThreshold, findAuremolPeaksAdaptive

PEAK_MODES = ['Positive & Negative','Positive only','Negative only']

class AuremolFrame(Frame):

  def __init__(self, guiParent, ccpnProject=None, **kw):

    self.guiParent = guiParent
    self.project = ccpnProject
    self.spectrum = None
    self.peakMode = 0
    
    if ccpnProject:
      self.nmrProject = ccpnProject.currentNmrProject
    else:
      self.nmrProject = None
    
    Frame.__init__(self, guiParent, **kw)
  
    self.expandGrid(0,0)

    options = ['Peak Picking',] #,'About Auremol' 'NOE assignment','Homology Modelling',]
    self.tabbedFrame = TabbedFrame(self, options=options)
    self.tabbedFrame.grid(row=0,column=0,sticky='nsew')
    frameA = self.tabbedFrame.frames[0]
    
    #frameC.grid_columnconfigure(0, weight=1)
    #frameC.grid_rowconfigure(0, weight=1)
    #frameD.grid_columnconfigure(0, weight=1)
    #frameD.grid_rowconfigure(0, weight=1)
    
    #
    # Frame A
    #
    frameA.expandGrid(2,0)
    frameA.expandGrid(3,0)
    frameA.expandGrid(4,0)
    frameA.expandGrid(5,0)
    
    
    frame = Frame(frameA, grid=(0,0))
    frame.expandGrid(0,4)
    
    label = Label(frame, text='Spectrum:', grid=(0,0))
    self.spectrumPulldown = PulldownList(frame, self.changeSpectrum, grid=(0,1))
    
    label = Label(frame, text='  Use Peak Sign:', grid=(0,2))
    self.peakModePulldown = PulldownList(frame, self.changePeakMode, texts=PEAK_MODES,
                                         objects=[0,1,2], grid=(0,3))
    
    
    frame = Frame(frameA, grid=(1,0))
    frame.expandGrid(0,4)
    
    label = Label(frame, text='Integration Depth (Relative to max):', grid=(1,0))
    self.segLevelEntry = FloatEntry(frame, text=0.1, grid=(1,1), width=8)
    
    label = Label(frame, text='Threshold (Threshold only):', grid=(1,3))
    self.thresholdEntry = IntEntry(frame, text=100000, grid=(1,4), width=8)
    
    label = Label(frame, text='Keep Peaks (Adaptive only):', grid=(1,5))
    self.keepPeakEntry = IntEntry(frame, text=4000, grid=(1,6), width=8)
    
    texts = ['Threshold\nPeak Pick','Adaptive\nPeak Pick']
    commands = [self.pickThreshold, self.pickAdaptive]
    self.buttons = ButtonList(frameA, texts=texts, commands=commands,
                              grid=(2,0),  sticky='NSEW')
    
    frame = Frame(frameA, grid=(3,0))
    frame.expandGrid(0,0)
    frame = Frame(frameA, grid=(4,0))
    frame.expandGrid(0,0)
    frame = Frame(frameA, grid=(5,0))
    frame.expandGrid(0,0)
     
    #
    # About
    """
    frameB.expandGrid(4,0)

    label = Label(frameB, text='References', font='Helvetica 12 bold')
    label.grid(row=0, column=0, sticky='w')
    
    text = 
    * Gronwald W, Brunner K, Kirchhofer R, Nasser A, Trenner J, Ganslmeier B,
    Riepl H, Ried A, Scheiber J, Elsner R, Neidig K-P, Kalbitzer HR
    AUREMOL, a New Program for the Automated Structure Elucidation of Biological Macromolecules
    Bruker Reports 2004; 154/155: 11-14

    * Ried A, Gronwald W, Trenner JM, Brunner K, Neidig KP, Kalbitzer HR
    Improved simulation of NOESY spectra by RELAX-JT2 including effects of J-coupling,
    transverse relaxation and chemical shift anisotrophy
    J Biomol NMR. 2004 Oct;30(2):121-31

    * Gronwald W, Moussa S, Elsner R, Jung A, Ganslmeier B, Trenner J, Kremer W, Neidig KP, Kalbitzer HR
    Automated assignment of NOESY NMR spectra using a knowledge based method (KNOWNOE)
    J Biomol NMR. 2002 Aug;23(4):271-87

    * Gronwald W, Kirchhofer R, Gorler A, Kremer W, Ganslmeier B, Neidig KP, Kalbitzer HR
    RFAC, a program for automated NMR R-factor estimation
    J Biomol NMR. 2000 Jun;17(2):137-51
    
    label = Label(frameB, text=text)
    label.grid(row=1, column=0, sticky='w')
    """
   
    #
    # Frame C
    #

    
    #
    # Frame D
    #

  
    self.updateAll()
  
  def getEntryData(self):
  
    segLevel = self.segLevelEntry.get() or 0.001
    threshold = self.thresholdEntry.get() or 100000
    maxPeaks = self.keepPeakEntry.get() or 1
  
    segLevel = min(1.0, segLevel)
  
    self.segLevelEntry.set(segLevel)
    self.thresholdEntry.set(threshold)
    self.keepPeakEntry.set(maxPeaks)
  
    return segLevel, threshold, maxPeaks
  
  def pickThreshold(self):
  
    if self.spectrum:
      segLevel, threshold, maxPeaks = self.getEntryData()
      try:
        findAuremolPeaksThreshold(spectrum=self.spectrum, mode=self.peakMode,
                                  useAutoThreshold=0, threshold=threshold, seglevel=segLevel)
      except Exception, e:
        showError('pickThreshold', str(e), parent=self)
  
  def pickAdaptive(self):
  
    if self.spectrum:
      segLevel, threshold, maxPeaks = self.getEntryData()
      try:
        findAuremolPeaksAdaptive(spectrum=self.spectrum, mode=self.peakMode,
                                 number=maxPeaks, seglevel=segLevel)
      except Exception, e:
        showError('pickAdaptive', str(e), parent=self)
      
  def getSpectra(self):
  
    spectra = []
    
    if self.nmrProject:
      for experiment in self.nmrProject.sortedExperiments():
        for spectrum in experiment.sortedDataSources():
          if spectrum.dataType == 'processed':
            spectra.append(spectrum)
  
    return spectra
  
  def updateSpectra(self):
  
    index = 0
    names = []
    spectrum = self.spectrum
    spectra = self.getSpectra()
    
    if spectra:
      names = ['%s:%s' % (s.experiment.name, s.name) for s in spectra]
      
      if spectrum not in spectra:
        spectrum = spectra[0]
        
      index = spectra.index(spectrum)  
    
    else:
      spectrum = None
      
    self.changeSpectrum(spectrum)
    self.spectrumPulldown.setup(names, spectra, index)  
  
  def changeSpectrum(self, spectrum):
  
    if self.spectrum is not spectrum:
      self.spectrum = spectrum
      # Any other dependent updates in future
  
  def changePeakMode(self, peakMode):
  
    self.peakMode = peakMode
  
  def open(self):
  
    self.updateAll()
  
    Frame.open(self)
  
  def updateAll(self, project=None):
  
    if project:
      self.project = project
      self.nmrProject = project.currentNmrProject
    
      if not self.nmrProject:
        self.nmrProject = project.newNmrProject(name=project.name)
    
    self.updateSpectra()

  def destroy(self):
  
    Frame.destroy(self)  


"""
======================COPYRIGHT/LICENSE START==========================

EditNoeClasses.py: Part of the CcpNmr Analysis program

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
from memops.general import Implementation

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes, getThroughSpacePeakLists
from ccpnmr.analysis.core.ConstraintBasic import getIntensityDistanceTable
from ccpnmr.analysis.core.Util import setSpectrumNoeDistanceClasses

from memops.gui.ButtonList import UtilityButtonList
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Label import Label
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ScrolledMatrix import ScrolledMatrix 

class EditNoeClassesPopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.guiParent = parent
    BasePopup.__init__(self, parent=parent, title='NOE Distance Classes', **kw)

  def body(self, guiFrame):
  
    self.noeClassChoice = None
    self.spectrum = None
    self.intensEntry = FloatEntry(self, returnCallback=self.setIntens, width=5)
    self.targetEntry = FloatEntry(self, returnCallback=self.setTarget, width=5)
    self.minEntry    = FloatEntry(self, returnCallback=self.setMin,    width=5)
    self.maxEntry    = FloatEntry(self, returnCallback=self.setMax,    width=5)
   
    row = 0

    label = Label(guiFrame, text='Spectrum: ', grid=(row,0))
    tipText = ''
    self.spectrumPulldown = PulldownMenu(guiFrame,self.changeSpectrum, grid=(row,1))

    row +=1

    guiFrame.expandGrid(row, 1)

    tipTexts = ['Lower bound of this intensity category. Values are relative to reference intensity.',
                'Target restraint distance for this category',
                'Lower bound distance for this category',
                'Upper bound distance for this category']
    headingList = ['Min. NOE\nIntensity','Target\nDist','Min\nDist','Max\nDist']
    editWidgets = [self.intensEntry,self.targetEntry,self.minEntry,self.maxEntry]
    editGetCallbacks = [self.getIntens,self.getTarget,self.getMin,self.getMax]
    editSetCallbacks = [self.setIntens,self.setTarget,self.setMin,self.setMax]
    
    self.noeClassMatrix = ScrolledMatrix(guiFrame,
                                         headingList=headingList,
                                         callback=self.selectClass,
                                         tipTexts=tipTexts,
                                         editWidgets=editWidgets,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         deleteFunc=self.deleteClass,
                                         grid=(row,0), gridSpan=(1,2))
                                         

    row +=1

    tipTexts = ['Add a new distance restraint category',
                'Deleted selected restraint categor']
    texts = ['Add Class','Delete Class']
    commands = [self.addClass,self.deleteClass]
    self.bottomButtons = UtilityButtonList(guiFrame, doClone=False, grid=(row,0),
                                           gridSpan=(1,2), tipTexts=tipTexts,
                                           commands=commands, texts=texts)

    for func in ('__init__','delete','setName'):
      self.registerNotify(self.updateSpectra, 'ccp.nmr.Nmr.Experiment', func)
      self.registerNotify(self.updateSpectra, 'ccp.nmr.Nmr.DataSource', func)

    self.updateSpectra()
    self.update()

  def open(self):
  
    self.updateSpectra()
    self.update()
    BasePopup.open(self)

  def updateSpectra(self, *opt):
    
    spectra = self.getSpectra()
    if not spectra:
      return
    
    names = [self.getSpectrumName(x) for x in spectra]
    if (not self.spectrum) or (self.spectrum not in spectra):
      self.spectrum = spectra[0]
    
    self.spectrumPulldown.setup(names, names.index(self.getSpectrumName(self.spectrum)) )
    
    self.update()

  def changeSpectrum(self, i, name):
  
    self.spectrum = self.getSpectra()[i]
    self.update()

  def getSpectrumName(self,spectrum):
  
    name = '%s:%s' % (spectrum.experiment.name,spectrum.name)
    return name
  
  def getSpectra(self):
  
    spectra = set()
    peakLists = getThroughSpacePeakLists(self.nmrProject)
    
    for peakList in peakLists:
      spectra.add(peakList.dataSource)
 
    spectra = list(spectra)
    spectra.sort()
 
    return spectra
    
  def selectClass(self, noeClass, row, col):
  
    if noeClass:
      self.noeClassChoice = (row, noeClass)
      
    if len(self.noeClassMatrix.objectList) > 1:
      self.bottomButtons.buttons[1].enable()
    else:
      self.bottomButtons.buttons[1].disable()

  def addClass(self):
  
    if self.spectrum:
      noeClass = [0.0,6.0,0.0,6.0]
 
      noeClasses = getIntensityDistanceTable(self.spectrum)
      noeClasses.append(noeClass)
      setSpectrumNoeDistanceClasses(self.spectrum, noeClasses)

      self.update()
  
  def deleteClass(self, *event):
    
    if self.spectrum:
      noeClasses = getIntensityDistanceTable(self.spectrum)

      if self.noeClassChoice and (self.noeClassChoice[1] in noeClasses):
        if len(noeClasses) > 1:
          (i,noeClass) = self.noeClassChoice
          noeClasses.remove(noeClass)
          self.noeClassChoice = None
          setSpectrumNoeDistanceClasses(self.spectrum, noeClasses)
          self.update()
    
  def setIntens(self, event):
  
    if self.noeClassChoice:
      val = self.intensEntry.get() or 0.0
      self.noeClassChoice[1][0] = val
      
    self.updateClass()
  
  def getIntens(self, row):
  
    if row:
      self.intensEntry.set(row[0])
  
  def setTarget(self, event):
  
    if self.noeClassChoice:
      val = self.targetEntry.get() or 0.0
      self.noeClassChoice[1][1] = val
      
    self.updateClass()
  
  def getTarget(self, row):
  
    if row:
      self.targetEntry.set(row[1])
  
  def setMin(self, event):
  
    if self.noeClassChoice:
      val = self.minEntry.get() or 0.0
      self.noeClassChoice[1][2] = val
      
    self.updateClass()
  
  def getMin(self, row):
  
    if row:
      self.minEntry.set(row[2])
  
  def setMax(self, event):
  
    if self.noeClassChoice:
      val = self.maxEntry.get() or 0.0
      self.noeClassChoice[1][3] = val
      
    self.updateClass()
  
  def getMax(self, row):
  
    if row:
      self.maxEntry.set(row[3])
    
  def getClasses(self):
  
    noeClasses = []
    if self.spectrum:
      noeClasses = getIntensityDistanceTable(self.spectrum)

    if noeClasses:
      for i in range(len(noeClasses)):
        (intens,target,minimum,maximum) = noeClasses[i]

        if minimum > maximum:
          (minimum,maximum) = (maximum,minimum)
        minimum = min(target, minimum)
        maximum = max(target, maximum)
        intens  = max(intens, 0.0)
        
        noeClasses[i] = [intens,target,minimum,maximum]
      noeClasses.sort()
      noeClasses.reverse()
    
    else:
      noeClasses = []
      if self.spectrum:
        # default
        noeClasses = getIntensityDistanceTable(self.spectrum)
      
    return noeClasses
    
  def updateClass(self):
  
    if self.spectrum and self.noeClassChoice:
      (i, noeClass) = self.noeClassChoice
      noeClasses = getIntensityDistanceTable(self.spectrum)
      noeClasses[i] = noeClass
      setSpectrumNoeDistanceClasses(self.spectrum, noeClasses)
      self.update()
    
  def update(self):


    textMatrix = []
    objectList = self.getClasses()
    
    if self.spectrum:
      if self.noeClassChoice and (len(objectList) > 1):
        self.bottomButtons.buttons[1].enable()
      else:
        self.bottomButtons.buttons[1].disable()
      self.bottomButtons.buttons[0].enable()
    else:
      self.bottomButtons.buttons[0].disable()
      self.bottomButtons.buttons[1].disable()
      
    for (intens,target,minimum,maximum) in objectList:
      datum = []
      datum.append(intens)
      datum.append(target)
      datum.append(minimum)
      datum.append(maximum)
      textMatrix.append(datum)
    
    self.noeClassMatrix.update(objectList=objectList,textMatrix=textMatrix)
    
    if self.spectrum:
      setSpectrumNoeDistanceClasses(self.spectrum,objectList)
  
  def destroy(self):

    for func in ('__init__','delete','setName'):
      self.unregisterNotify(self.updateSpectra, 'ccp.nmr.Nmr.Experiment', func)
      self.unregisterNotify(self.updateSpectra, 'ccp.nmr.Nmr.DataSource', func)

    BasePopup.destroy(self)

   

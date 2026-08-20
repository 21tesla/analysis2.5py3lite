
"""
======================COPYRIGHT/LICENSE START==========================

NoeDistParamsFrame.py: Part of the CcpNmr Analysis program

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

from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix 

from ccpnmr.analysis.core.ConstraintBasic import getMeanPeakIntensity, getIntensityDistanceTable
from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes, getThroughSpacePeakLists
from ccpnmr.analysis.core.Util import getAnalysisPeakList, setSpectrumNoeDistanceClasses

noeDistanceFunctions = ['intensity^-1/6','intensity^-1/4',
                        'intensity^-1/3','Distance bins']

ANGSTROM = u'\u00C5'

class RowObject:
  
  def __init__(self, intensity, targetDist, minDist, maxDist):
    self.intensity = intensity
    self.targetDist = targetDist
    self.minDist = minDist
    self.maxDist = maxDist

class NoeDistParamsFrame(Frame):

  def __init__(self, parent, nmrProject, peakList=None,
               intensityTypeCallback=None, distanceFunctionCallback=None,
               *args, **kw):

    self.nmrProject = nmrProject
    self.peakList = peakList
    self.analysisPeakList = None
    self.parent = parent
    self.intensityTypeCallback = intensityTypeCallback
    self.distanceFunctionCallback = distanceFunctionCallback

    Frame.__init__(self, parent=parent, **kw)
    self.grid_rowconfigure(1, weight=1)
    self.grid_columnconfigure(0, weight=1)

    labelFrame = LabelFrame(self, text='Restraint Distance Params')
    labelFrame.grid_columnconfigure(5, weight=1)
    labelFrame.grid(row=0, column=0, sticky='ew')

    label = Label(labelFrame, text='Intensity type: ')
    label.grid(row=0, column=0, sticky='ne')
    entries = ['height','volume']
    self.intensityType = entries[0]
    self.intensityTypePulldown = PulldownList(labelFrame, callback=self.changeIntensityType,
                                             texts=entries,
                                             tipText='The kind of intensity to estimate distances with')
    self.intensityTypePulldown.grid(row=0, column=1, sticky='nw')
    
    label = Label(labelFrame, text='Distance function: ')
    label.grid(row=0, column=2, sticky='ne')
    self.distFunction = noeDistanceFunctions[0]
    self.distFunctionPulldown = PulldownList(labelFrame, callback=self.changeDistFunction,
                                            texts=noeDistanceFunctions,
                                            tipText='How to derive restraint distances from peak intensity')
    self.distFunctionPulldown.grid(row=0, column=3, sticky='nw')

    label = Label(labelFrame, text='Only assigned\nref intensities:')
    label.grid(row=0, column=4,  sticky='ne')
    self.assignedCheck = CheckButton(labelFrame, callback=self.toggleUseAssigned,
                                     tipText='Whether to consider only assigned peaks when calculating average, reference intensity')
    self.assignedCheck.grid(row=0, column=5, sticky='nw')
    self.assignedCheck.set(0)

    label = Label(labelFrame, text='Ref intensity: ')
    label.grid(row=1, column=0, sticky='ne')
    self.scaleEntry = FloatEntry(labelFrame, text=1.0,
                                 returnCallback=None, width=8,
                                 tipText='The peak intensity that corresponds to the reference distance. Defaults to peak list mean')
    self.scaleEntry.grid(row=1, column=1, sticky='nw')

    self.label1 = Label(labelFrame, text='Ref distance (%s): ' % ANGSTROM)
    self.refDistEntry = FloatEntry(labelFrame, text=3.2,
                                   returnCallback=self.setRefDist, width=8,
                                   tipText='The distance corresponding to the reference peak intensity; to perform calibration')
    self.refDistEntry.bind('<Leave>' , self.setRefDist, '+')

    self.label2 = Label(labelFrame, text='Lower frac error: ')
    self.lowerErrorEntry = FloatEntry(labelFrame, text=0.20,
                                      returnCallback=None, width=8,
                                      tipText='The fractional lower restraint bound, defaults to 20% closer')
 
    self.label3 = Label(labelFrame, text='Upper frac error: ')
    self.upperErrorEntry = FloatEntry(labelFrame, text=0.20, returnCallback=None, width=8,
                                      tipText='The fractional upper restraint bound bound, defaults to 20% further')

    self.label4 = Label(labelFrame, text='Lower dist limit (%s): ' % ANGSTROM)
    self.lowerLimEntry = FloatEntry(labelFrame, text=1.72, returnCallback=None, width=8,
                                    tipText='The absolute lower distance limit, overriding fractional error')

    self.label5 = Label(labelFrame, text='Upper dist limit (%s): ' % ANGSTROM)
    self.upperLimEntry = FloatEntry(labelFrame, text=8.0, returnCallback=None, width=8,
                                    tipText='The absolute upper distance limit, overriding fractional error')
  
    self.binsFrame = EditNoeClassesFrame(self, nmrProject)
    
    self.gridNoeParamsEntries()
  
    if self.peakList:
      self.update(self.peakList)
  
  def setRefDist(self, *event):
      
    if self.analysisPeakList:
      self.analysisPeakList.noeRefDistance = self.refDistEntry.get() or 3.2
  
  def gridNoeBinsButton(self):
  
    self.binsFrame.grid(row=1, column=0, sticky='nsew')
    
    self.label1.grid_forget()
    self.refDistEntry.grid_forget()
    self.label2.grid_forget()
    self.lowerErrorEntry.grid_forget()
    self.label3.grid_forget()
    self.upperErrorEntry.grid_forget()
    self.label4.grid_forget()
    self.lowerLimEntry.grid_forget()
    self.label5.grid_forget()
    self.upperLimEntry.grid_forget()

  
  def gridNoeParamsEntries(self):
  
    self.binsFrame.grid_forget()

    self.label1.grid(row=2, column=0, sticky='ne')
    self.refDistEntry.grid(row=2, column=1, sticky='nw')

    self.label2.grid(row=1, column=2, sticky='ne')
    self.lowerErrorEntry.grid(row=1, column=3, stick='nw')
 
    self.label3.grid(row=2, column=2, sticky='ne')
    self.upperErrorEntry.grid(row=2, column=3, stick='nw')

    self.label4.grid(row=1, column=4, sticky='ne')
    self.lowerLimEntry.grid(row=1, column=5, stick='nw')

    self.label5.grid(row=2, column=4, sticky='ne')
    self.upperLimEntry.grid(row=2, column=5, stick='nw')
  
  def toggleUseAssigned(self, boolean):

    self.updateBaseline()
  
  def getScale(self):
  
    return self.scaleEntry.get() or 0.0
    
  def setDistanceClasses(self):
  
    if hasattr(self.parent.parent, 'guiParent'):
      self.parent.parent.guiParent.editNoeClasses()
    else:
      self.parent.parent.parent.editNoeClasses()

  def changeDistFunction(self, name):

    if self.distFunction != name:
      self.distFunction = name

      if name == noeDistanceFunctions[-1]:
        self.gridNoeBinsButton()
 
      else:
        self.gridNoeParamsEntries()
        
      self.update(self.peakList)

      if self.distanceFunctionCallback:
        self.distanceFunctionCallback(name)

  def changeIntensityType(self, intensityType):
 
    self.intensityType = intensityType

    if self.analysisPeakList:
      self.analysisPeakList.noeIntensityType = intensityType

    self.updateBaseline()

    if self.intensityTypeCallback:
      self.intensityTypeCallback(intensityType)

  def updateIntensityTypes(self, *obj):
  
    if self.analysisPeakList:
      self.intensityType = self.analysisPeakList.noeIntensityType
      self.intensityTypePulldown.set(self.intensityType)
  
  def updateBaseline(self):
      
    onlyAssigned = self.assignedCheck.get()
    peaks = []
    
    if self.peakList:
      if onlyAssigned:
        for peak in self.peakList.peaks:
          for peakDim in peak.peakDims:
            if peakDim.peakDimContribs:
              peaks.append(peak)
              break
 
      else:
        peaks = self.peakList.peaks
    
    
    if self.peakList and peaks:
      value = getMeanPeakIntensity(peaks, intensityType=self.intensityType)
      self.scaleEntry.set(value)
      self.analysisPeakList.noeRefIntensity = value or None
    else:
      self.scaleEntry.set(0.0)

  
  def getParams(self):

    if self.distFunction == noeDistanceFunctions[-1]:
      return None
    
    refDist  = self.refDistEntry.get() or 3.2
    negError = self.lowerErrorEntry.get() or 0.0
    posError = self.upperErrorEntry.get() or 0.0
    absMin   = self.lowerLimEntry.get() or 0.0
    absMax   = self.upperLimEntry.get() or 8.0
    
    if self.distFunction == noeDistanceFunctions[2]:
      power = 3.0
    elif self.distFunction  == noeDistanceFunctions[1]:
      power = 4.0
    else:
      power = 6.0  
    
    return refDist, negError, posError, absMin, absMax, power

  def update(self, peakList=None):
 
    self.peakList = peakList
    #self.updateDistFunctions()
    
    if self.peakList:
      self.binsFrame.changeSpectrum(peakList.dataSource)
      self.analysisPeakList = getAnalysisPeakList(self.peakList)
      self.refDistEntry.set(self.analysisPeakList.noeRefDistance or 3.2)
      
    else:
      self.analysisPeakList = None
    
    self.updateIntensityTypes()
    self.updateBaseline()
      
class EditNoeClassesFrame(LabelFrame):

  def __init__(self, parent, nmrProject, *args, **kw):

    self.parent = parent
    self.nmrProject = nmrProject
    LabelFrame.__init__(self, parent=parent, text='Binned Distance Classes', **kw)
  
    self.noeClassChoice = None
    self.spectrum = None
    self.intensEntry = FloatEntry(self.parent, returnCallback=self.setIntens, width=5)
    self.targetEntry = FloatEntry(self.parent, returnCallback=self.setTarget, width=5)
    self.minEntry    = FloatEntry(self.parent, returnCallback=self.setMin,    width=5)
    self.maxEntry    = FloatEntry(self.parent, returnCallback=self.setMax,    width=5)
    
    self.grid_columnconfigure(0, weight=1)
   
    #label = Label(self, text='Spectrum: ')
    #label.grid(row=row,column=0, sticky='nw')
    #self.spectrumPulldown = PulldownList(self,callback=self.changeSpectrum)
    #self.spectrumPulldown.grid(row=row, column=1, sticky='nw')
    #row +=1

    row = 0

    self.grid_rowconfigure(row, weight=1)
    tipTexts = ['Lower bound of this intensity category. Values are relative to reference intensity.',
                'Target restraint distance for this category',
                'Lower bound distance for this category',
                'Upper bound distance for this category']
    headingList = ['Min. Peak\nIntensity','Target\nDist',
                   'Min\nDist','Max\nDist']
    editWidgets = [self.intensEntry,self.targetEntry,
                   self.minEntry,self.maxEntry]
    editGetCallbacks = [self.getIntens,self.getTarget,
                        self.getMin,self.getMax]
    editSetCallbacks = [self.setIntens,self.setTarget,
                        self.setMin,self.setMax]
    
    self.noeClassMatrix = ScrolledMatrix(self, grid=(row,0), gridSpan=(None, 2),
                                         headingList=headingList,
                                         callback=self.selectClass,
                                         tipTexts=tipTexts,
                                         editWidgets=editWidgets,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         deleteFunc=self.deleteClass)
                                         

    row +=1

    tipTexts = ['Add a new distance restraint category',
                'Deleted selected restraint category',
                'Copy restraint categories from selected spectrum']
    texts = ['Add Class','Delete Class','Copy From:']
    commands = [self.addClass,self.deleteClass,self.copyFromSpec]
    self.bottomButtons = ButtonList(self, tipTexts=tipTexts, grid=(row,0),
                                    commands=commands, texts=texts)
    self.copySpecPulldown = PulldownList(self, callback=None, grid=(row,1),
                                         tipText='Select spectrum to copy NOE classes from')

    self.updateSpectra()
    self.update()

  def open(self):
  
    self.updateSpectra()
    self.update()

  def updateSpectra(self, *opt):
    
    spectra = self.getSpectra()
    if self.spectrum not in spectra:
      self.spectrum = None
    
    index = 0
    names = []
    if self.spectrum:
      spectra.remove(self.spectrum)
      names = [self.getSpectrumName(s) for s in spectra]

    self.copySpecPulldown.setup(names, spectra, index)
    
    self.update()
    
  def copyFromSpec(self):
  
    if self.spectrum:
      spectrum2 = self.copySpecPulldown.getObject()
      
      if spectrum2:
        noeClasses = getIntensityDistanceTable(spectrum2)
        setSpectrumNoeDistanceClasses(self.spectrum, noeClasses)
 
        self.update()
  
  def changeSpectrum(self, spectrum):
  
    self.spectrum = spectrum
    self.updateSpectra()
    self.update()

  def getSpectrumName(self,spectrum):
  
    name = '%s:%s' % (spectrum.experiment.name,spectrum.name)
    return name
  
  def getSpectra(self):
  
    spectra = set()
    peakLists = getThroughSpacePeakLists(self.nmrProject.root)
    
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

      if self.noeClassChoice:
        choice = self.noeClassChoice[1]
        noeClass = [choice.intensity, choice.targetDist, choice.minDist, choice.maxDist]
        if noeClass in noeClasses:
          noeClasses.remove(noeClass)
          self.noeClassChoice = None
          setSpectrumNoeDistanceClasses(self.spectrum, noeClasses)
          self.update()
    
  def setIntens(self, event):
  
    if self.noeClassChoice:
      val = self.intensEntry.get() or 0.0
      self.noeClassChoice[1].intensity = val
      
    self.updateClass()
  
  def getIntens(self, row):
  
    if row:
      self.intensEntry.set(row.intensity)
  
  def setTarget(self, event):
  
    if self.noeClassChoice:
      val = self.targetEntry.get() or 0.0
      self.noeClassChoice[1].targetDist = val
      
    self.updateClass()
  
  def getTarget(self, row):
  
    if row:
      self.targetEntry.set(row.targetDist)
  
  def setMin(self, event):
  
    if self.noeClassChoice:
      val = self.minEntry.get() or 0.0
      self.noeClassChoice[1].minDist = val
      
    self.updateClass()
  
  def getMin(self, row):
    
    if row:
      self.minEntry.set(row.minDist)
  
  def setMax(self, event):
  
    if self.noeClassChoice:
      val = self.maxEntry.get() or 0.0
      self.noeClassChoice[1].maxDist = val
      
    self.updateClass()
  
  def getMax(self, row):
  
    if row:
      self.maxEntry.set(row.maxDist)
    
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
      noeClasses[i] = (noeClass.intensity, noeClass.targetDist, noeClass.minDist, noeClass.maxDist)
      setSpectrumNoeDistanceClasses(self.spectrum, noeClasses)
      self.update()
    
  def update(self):


    textMatrix = []
    classList = self.getClasses()
    objectList = []
    
    if self.spectrum:
      if self.noeClassChoice and (len(objectList) > 1):
        self.bottomButtons.buttons[1].enable()
      else:
        self.bottomButtons.buttons[1].disable()
      self.bottomButtons.buttons[0].enable()
    else:
      self.bottomButtons.buttons[0].disable()
      self.bottomButtons.buttons[1].disable()
      
    for (intens,target,minimum,maximum) in classList:
      datum = []
      datum.append(intens)
      datum.append(target)
      datum.append(minimum)
      datum.append(maximum)
      textMatrix.append(datum)
      objectList.append(RowObject(intens,target,minimum,maximum))
    
    self.noeClassMatrix.update(objectList=objectList,textMatrix=textMatrix)
    
    if self.spectrum:
      setSpectrumNoeDistanceClasses(self.spectrum,classList)
  
  def destroy(self):

    LabelFrame.destroy(self)

   

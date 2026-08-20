"""
======================COPYRIGHT/LICENSE START==========================

EditExperimentSeries.py: Part of the CcpNmr Analysis program

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

from memops.gui import Color
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.DataEntry import askInteger
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.LabelDivider import LabelDivider
from memops.gui.Label import Label
from memops.gui.MessageReporter import showOkCancel, showWarning
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.DataAnalysisBasic import getExperimentConditionSet, getNmrExpSeriesSampleConditions
from ccpnmr.analysis.core.ExperimentBasic import getSampledDimExperiments, getExperimentSampledDim
from ccpnmr.analysis.core.Util import getAnalysisSpectrum, getHueSortedColorSchemes, getAnalysisDataDim

#'concentration':['M','mM','uM','nM','pM','fM'],
#                      'spin-lock field':['Hz'],
#                      'spin-lock offset':['Hz'],

CONDITION_UNITS_DICT ={
                      'temperature':['K','C','F'],
                      'pressure':['MPa','kPa','Pa','atm','bar','mbar','ksi','psi','torr','at','mmHg'],
                      'pH':['pH'],
                      'concentration':['M','mM','uM','nM','pM','fM'],
                      'ligand conc.':['M','mM','uM','nM','pM','fM'],
                      'ligand ratio':['None'],
                      'initial temperature':['K','C','F'],
                      'final temperature':['K','C','F'],
                      'delay time':['ms','s','us','ns','mins','hours','days'],
                      'mixing time':['ms','s','us','ns','mins','hours','days'],
                      'Time':['ms','s','us','ns','mins','hours','days'],
                      'num delays':['None'],
                      'gradient strength':['%','G/cm'],
                      'pulsing frequency':['Hz','kHz','MHz','GHz','THz','mHz','uHz'],
                      'frequency':['Hz','kHz','MHz','GHz','THz','mHz','uHz'],
                      }


# TBD
# 'Blank experiment'?
# Chnaging condition type more easily

class EditExperimentSeriesPopup(BasePopup):
  """
  **Setup Experiment Series for Chemical Shift and Intensity Changes**
  
  The purpose of this popup is to setup ordered groups of experiments that are
  related by the variation in some condition or parameter, but otherwise the
  kind of experiment being run is the same. For example the user could setup
  experiments for a temperature series, ligand binding titration or relaxation
  rate measurement.

  The layout is divided into two tables. The upper table shows all of the series
  that are known to the current project and the lower table shows all of the
  experiments or planes that comprise the selected series. These series relate
  to both groups of separate experiments (and hence spectra) and also single
  experiment where one dimension is not for NMR frequency but rather a "sampled"
  dimension and thus effectively combines many experiments (e.g. for different
  T1 values) as different planes. A stack of effectively 2D experiments combined
  in this manner is typically referred to as pseudo-3D. - The experiment is 3D
  but there are only two NMR dimensions.

  Series which are stacks of planes in a single experiment entity are
  automatically detected and entered into the table once they are loaded. Series
  that are constructed from multiple, separate experiments however must be setup
  by the user. To setup a new series of this kind [Add Series] makes a new,
  empty NMR series. Next the user should change the "Parameter varied" column to
  specify what type of thing is varied between the different experiments. The
  user then adds experiments into this series with the [Add Series Point]
  function at the bottom. Once the points have been added to the series the name
  of the experiment for each point may be changed. Initially, arbitrary
  experiments appear for the series points, so these almost always have to be
  adjusted.

  Once a stacked-plane experiment or series of experiments is setup, the user
  next sets (or checks) the value of the parameter associated with each point in
  the lower table. When loading stacked-plane experiments these values may come
  though automatically, if they are present in the spectrum header or parameter
  file. Given a completed NMR series specification the user may then extract the
  relevant data from the series using one of the analysis tools like `Follow
  Intensity Changes`_ or `Follow Shift Changes`_.

  **Caveats & Tips**

  Make sure the "Parameter Varied" for a given NMR series is appropriate to the
  type of analysis being performed. Many tools that extract T1 or Kd
  measurements for example look for specific types of series.

  The "Set/Unset Ref Plane" function is only used in certain kinds of series
  such as those that use trains of CPMG pulses. 

  .. _`Follow Intensity Changes`: CalcRatesPopup.html
  .. _`Follow Shift Changes`: FollowShiftChangesPopup.html

  """

  def __init__(self, parent, *args, **kw):

    self.guiParent   = parent
    self.expSeries   = None
    self.conditionPoint   = None
    self.waiting     = 0
  
    BasePopup.__init__(self, parent, title="Experiment : NMR Series", **kw)

  def body(self, guiFrame):

    self.geometry("500x500")

    self.nameEntry  = Entry(self, text='', returnCallback=self.setName,    width=12)
    self.detailsEntry = Entry(self, text='', returnCallback=self.setDetails, width=16)
    self.valueEntry = FloatEntry(self, text='', returnCallback=self.setValue, width=10)
    self.errorEntry = FloatEntry(self, text='', returnCallback=self.setError, width=8)
    
    self.conditionNamesPulldown = PulldownList(self, callback=self.setConditionName,
                                               texts=self.getConditionNames())
    self.unitPulldown = PulldownList(self, callback=self.setUnit,
                                     texts=self.getUnits())
    self.experimentPulldown = PulldownList(self, callback=self.setExperiment)

    guiFrame.grid_columnconfigure(0, weight=1)

    row = 0
    frame = Frame(guiFrame, grid=(row, 0))
    frame.expandGrid(None,0)
    div = LabelDivider(frame, text='Current Series', grid=(0, 0))
    utilButtons = UtilityButtonList(frame, helpUrl=self.help_url, grid=(0,1))

    row += 1
    frame0 = Frame(guiFrame, grid=(row, 0))
    frame0.expandGrid(0,0)
    tipTexts = ['The serial number of the experiment series, but left blank if the series as actually a pseudo-nD experiment (with a sampled non-frequency axis)',
                'The name of the experiment series, which may be a single pseudo-nD experiment',
                'The number of separate experiments (and hence spectra) present in the series',
                'The kind of quantity that varies for different experiments/planes within the NMR series, e.g. delay time, temperature, ligand concentration etc.',
                'The number of separate points, each with a separate experiment/plane and parameter value, in the series']
    headingList      = ['#','Name','Experiments','Parameter\nVaried','Num\nPoints']
    editWidgets      = [None, self.nameEntry, None, self.conditionNamesPulldown, None]
    editGetCallbacks = [None, self.getName,   None, self.getConditionName, None]
    editSetCallbacks = [None, self.setName,   None, self.setConditionName, None]
    self.seriesMatrix = ScrolledMatrix(frame0, tipTexts=tipTexts,
                                       editSetCallbacks=editSetCallbacks,
                                       editGetCallbacks=editGetCallbacks,
                                       editWidgets=editWidgets,
                                       headingList=headingList,
                                       callback=self.selectExpSeries,
                                       deleteFunc=self.deleteExpSeries,
                                       grid=(0,0), gridSpan=(None, 3))

    tipTexts = ['Make a new, blank NMR series specification in the CCPN project',
                'Delete the selected NMR series from the project, although any component experiments remain. Note you cannot delete pseudo-nD series; delete the actual experiment instead',
                'Colour the spectrum contours for each experiment in the selected series (not pseudo-nD) using a specified scheme']
    texts    = ['Add Series','Delete Series',
                'Auto Colour Spectra']
    commands = [self.addExpSeries,self.deleteExpSeries,
                self.autoColorSpectra]
                
    self.seriesButtons = ButtonList(frame0, texts=texts, commands=commands,
                                    grid=(1,0), tipTexts=tipTexts)

    label = Label(frame0, text='Scheme:', grid=(1,1))
    
    tipText = 'Selects which colour scheme to apply to the contours of (separate) experiments within an NMR series'
    self.colorSchemePulldown = PulldownList(frame0, grid=(1,2), tipText=tipText)

    row += 1
    div = LabelDivider(guiFrame, text='Experimental Parameters & Conditions', grid=(row, 0))

    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    frame1 = Frame(guiFrame, grid=(row, 0))
    frame1.expandGrid(0,0)
    tipTexts = ['The kind of experimental parameter that is being used to define the NMR series',
                'The experiment that corresponds to the specified parameter value; can be edited from an arbitrary initial experiment',
                'The numeric value of the parameter (condition) that relates to the experiment or point in the NMR series',
                'The estimated error in value of the condition',
                'The measurement unit in which the value of the condition is represented']
    headingList      = ['Parameter','Experiment','Value','Error','Unit']
    editWidgets      = [None,self.experimentPulldown,self.valueEntry,self.errorEntry, self.unitPulldown]
    editGetCallbacks = [None,self.getExperiment,     self.getValue,  self.getError,   self.getUnit]
    editSetCallbacks = [None,self.setExperiment,     self.setValue,  self.setError,   self.setUnit]
    self.conditionPointsMatrix = ScrolledMatrix(frame1, grid=(0,0), tipTexts=tipTexts,
                                                editSetCallbacks=editSetCallbacks,
                                                editGetCallbacks=editGetCallbacks,
                                                editWidgets=editWidgets,
                                                headingList=headingList,
                                                callback=self.selectConditionPoint,
                                                deleteFunc=self.deleteConditionPoint)
    
    self.conditionPointsMatrix.doEditMarkExtraRules = self.conditionTableShow 
    tipTexts = ['Add a new point to the NMR series with an associated parameter value and experiment',
                'Remove the selected point from the series, including any associated parameter value',
                'For appropriate kinds of NMR series, set or unset a point as representing the plane to use as a reference']
    texts    = ['Add Series Point','Delete Series Point','Set/Unset Ref Plane']
    commands = [self.addConditionPoint,self.deleteConditionPoint,self.setSampledReferencePlane]
    self.conditionPointsButtons = ButtonList(frame1, texts=texts, commands=commands,
                                             tipTexts=tipTexts, grid=(1,0))
    
    self.updateAfter()
    self.updateColorSchemes()

    self.administerNotifiers(self.registerNotify)
    
  def administerNotifiers(self, notifyFunc):

    #for func in ('__init__', 'delete','setName'):
    for func in ('__init__', 'delete','setName','setConditionNames',
                 'addConditionName','removeConditionName'):
      notifyFunc(self.updateAfter,'ccp.nmr.Nmr.NmrExpSeries', func)

    for func in ('__init__', 'delete','setName'):
      notifyFunc(self.updateExperiments,'ccp.nmr.Nmr.Experiment', func)

    for func in ('__init__', 'delete'):
      notifyFunc(self.updateDataDim,'ccp.nmr.Nmr.SampledDataDim', func)

    for func in ('__init__', 'delete','setCondition','setUnit','setValue','setError'):
      notifyFunc(self.updateConditionsAfter,'ccp.nmr.Nmr.SampleCondition', func)

    for func in ('__init__', 'delete','setCondition'):
      notifyFunc(self.updateAfter,'ccp.nmr.Nmr.SampleCondition', func)
      
    for func in ('setConditionVaried', 'setPointErrors', 
                 'addPointError', 'removePointError', 
                 'setPointValues','addPointValue', 
                 'removePointValue','setUnit'): 
      notifyFunc(self.updateAfter,'ccp.nmr.Nmr.SampledDataDim', func)    

    for func in ('__init__', 'delete'):
      notifyFunc(self.updateColorSchemes,'ccpnmr.AnalysisProfile.ColorScheme', func)  

  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)
  
  def updateColorSchemes(self, scheme=None):
    
    index = 0
    prevScheme = self.colorSchemePulldown.getObject()
    
    schemes = getHueSortedColorSchemes(self.analysisProfile)
    schemes = [s for s in schemes if len(s.colors) > 1]
    colors = [list(s.colors) for s in schemes]
    
    if schemes:
      names = [s.name for s in schemes]
      
      if prevScheme in schemes:
        index = schemes.index(prevScheme)
    
    else:
      names  = []
  
    self.colorSchemePulldown.setup(names, schemes, index, colors)
  
  def autoColorSpectra(self):
  
    if self.expSeries and (self.expSeries.className != 'Experiment'):
      scheme = self.colorSchemePulldown.getObject()
      
      if scheme:
        colors = scheme.colors
      else:
        colors = ['#FF0000','#00FF00','#0000FF']  
      
      cdict = getNmrExpSeriesSampleConditions(self.expSeries)
      conditionName = list(self.expSeries.conditionNames)[0]
      expList = []
       
      for sampleCondition in cdict.get(conditionName, []):
        expList.append( (sampleCondition.value, sampleCondition.parent.experiments) )
 
      expList.sort()
      m = len(expList)-1.0
      c = len(colors)-1
       
      for i, (value, experiments) in enumerate(expList):
        p = c*i/m
        j = int(p)
        
        r1, g1, b1 = Color.hexToRgb(colors[j])
        r2, g2, b2 = Color.hexToRgb(colors[min(c,j+1)])
        
        f2 = p-j
        f1 = 1.0-f2
      
        r = (r1*f1)+(r2*f2)
        g = (g1*f1)+(g2*f2)
        b = (b1*f1)+(b2*f2)

        hexColor = Color.hexRepr(r,g,b)
       
        for experiment in experiments:
          for spectrum in experiment.dataSources:
            if spectrum.dataType == 'processed':
              analysisSpec = getAnalysisSpectrum(spectrum)
              
              if analysisSpec.posColors:
                analysisSpec.posColors = [hexColor,]
              elif analysisSpec.negColors:
                analysisSpec.negColors = [hexColor,]
  
  def getUnusedExperiments(self):
  
    sampleExperiments = getSampledDimExperiments(self.nmrProject)
    
    experiments = []
    for experiment in self.nmrProject.sortedExperiments():
      if experiment in sampleExperiments:
        continue
    
      if self.expSeries and (self.expSeries.className != 'Experiment'):
        if experiment in self.expSeries.experiments:
          continue
            
      experiments.append(experiment)
    
    return experiments
  
  def conditionTableShow(self, object, row, col):
  
    if type(object) is type(()):
      dataDim, index = object
      refPlane = dataDim.analysisDataDim.refSamplePlane
      if refPlane == index:
        return False
      if col == 1:
        return False
      
    return True
  
  def setSampledReferencePlane(self):
  
    if self.expSeries and (self.expSeries.className == 'Experiment'):
      if self.conditionPoint:
        dataDim, point = self.conditionPoint
        analysisDataDim = dataDim.analysisDataDim
        refPoint = analysisDataDim.refSamplePlane
        
        if refPoint == point:
          analysisDataDim.refSamplePlane = None
        else:
          analysisDataDim.refSamplePlane = point
    
        self.updateAfter()
    
  def checkAddSampleCondition(self, experiment):

    conditionSet = getExperimentConditionSet(experiment)
    conditionName = self.expSeries.conditionNames[0]
    condDict = getNmrExpSeriesSampleConditions(self.expSeries)

    sampleConditions = condDict.get(conditionName, [])
    units = [sc.unit for sc in sampleConditions if sc.unit]

    if units:
      sortList = list(set(units))
      sortList.sort(key = lambda u:units.count(u))
      unit = sortList[-1]
    else:
      unit = CONDITION_UNITS_DICT[conditionName][0]

    condition = conditionSet.findFirstSampleCondition(condition=conditionName)
    if not condition:
      condition = conditionSet.newSampleCondition(condition=conditionName,
                                                  unit=unit, value=0.0, error=0.0)
                                                  
  def addConditionPoint(self):
  
    if self.expSeries and (self.expSeries.className != 'Experiment'):
      
      experiments = self.getUnusedExperiments()
      if not experiments:
        showWarning('Warning','No experiments available', parent=self)
        return
      
      experiment = experiments[0]
      
      if experiment not in self.expSeries.experiments:
        self.expSeries.addExperiment(experiment)
      
      self.checkAddSampleCondition(experiment)
      self.updateAfter()

  def deleteConditionPoint(self, *event):
  
    if self.conditionPoint and (self.expSeries.className != 'Experiment'):
      if showOkCancel('Confirm','Really delete series point?', parent=self):
        conditionSet = self.conditionPoint.sampleConditionSet 
        
        experiments = [e for e in conditionSet.experiments if e in self.expSeries.experiments]
        
        for experiment in experiments:
          self.expSeries.removeExperiment(experiment)

        for experiment in experiments:
          for expSeries in experiment.nmrExpSeries:
             if self.conditionPoint.condition in self.expSeries.conditionNames:
               break
          else:
            continue
          break
          
        else:    
          self.conditionPoint.delete()
          
        self.conditionPoint = None
           

  def selectConditionPoint(self, object, row, col):
  
    if object:
      self.conditionPoint = object
      self.updateButtons()

  def selectExpSeries(self, object, row, col):

    if object:
      self.expSeries = object
      self.checkExperimentConditionsConsistent()
      self.updateConditions()

  def checkExperimentConditionsConsistent(self):

    if self.expSeries.className != 'Experiment':
      for experiment in self.expSeries.experiments:
        self.checkAddSampleCondition(experiment)

  def getUnits(self):
  
    units = []
    if self.expSeries:
      if self.expSeries.className == 'Experiment':
        conditionName = getExperimentSampledDim(self.expSeries).conditionVaried
     
      else:
        conditionName = self.expSeries.conditionNames[0]
        
      units = CONDITION_UNITS_DICT.get(conditionName)
      if not units:
        units = ['?',]
    
    return units


  def getUnit(self, sampleCondition):

    index = -1
    units = self.getUnits()
    if units:
      if sampleCondition:
        if type(sampleCondition) is type(()):
          dataDim, index = sampleCondition
          unit = dataDim.unit
        else:
          unit = sampleCondition.unit
      
      if unit not in units:
        unit = units[0]
      
      index = units.index(unit)
    
    self.unitPulldown.setup(units,units,index)
  
  def setUnit(self, obj):
  
    name = self.unitPulldown.getObject()
      
    if self.conditionPoint:
      if type(self.conditionPoint) is type(()):
        dataDim, index = self.conditionPoint
        dataDim.setUnit(name)
      
      else:
        self.conditionPoint.setUnit(name)

  def getConditionNames(self):
    
    if self.expSeries and (self.expSeries.className == 'Experiment'):
      names = ['delay time','mixing time','num delays','pulsing frequency','gradient strength']
    else:
      names = CONDITION_UNITS_DICT.keys()
      names.sort()
      
    return names

  def setConditionName(self, obj):
  
    name = self.conditionNamesPulldown.getObject()
  
    if self.expSeries:
      if self.expSeries.className == 'Experiment':
        dataDim = getExperimentSampledDim(self.expSeries)
        dataDim.conditionVaried = name
        
      else:
        self.expSeries.setConditionNames([name,])

  def getConditionName(self, expSeries):
  
    index = 0
    names = self.getConditionNames()
    
    if names:
    
      if expSeries:
        if expSeries.className == 'Experiment':
          name  = getExperimentSampledDim(expSeries).conditionVaried
        else:
          name  = expSeries.conditionNames[0]
        
        if name:  
          index = names.index(name)
      
      else:
        index = 0
    
    self.conditionNamesPulldown.setup(names,names,index)

  def getName(self, expSeries):
    
    if expSeries :
      self.nameEntry.set(expSeries.name)
 
  def setName(self, event):

    text = self.nameEntry.get()
    if text and text != ' ':
      self.expSeries.setName( text )

  def getValue(self, conditionPoint):
    
    if conditionPoint:
      if type(self.conditionPoint) is type(()):
        dataDim, index = conditionPoint
        value = dataDim.pointValues[index]
        
      else:
        value = conditionPoint.value
        
      self.valueEntry.set(value)
 
  def setValue(self, event):

    value = self.valueEntry.get()
    if value is not None:
      if type(self.conditionPoint) is type(()):
        dataDim, index = self.conditionPoint
        values = list(dataDim.pointValues)
        values[index] = value
        dataDim.setPointValues(values)
        
      else:
        self.conditionPoint.setValue( value )

  def getError(self, conditionPoint):
    
    if conditionPoint:
      if type(self.conditionPoint) is type(()):
        dataDim, index = conditionPoint
        
        if index < len(dataDim.pointErrors):
          error = dataDim.pointValues[index]
        else:
          error = 0.0
      
      else:
        error = conditionPoint.error
      
      self.errorEntry.set(error)
 
  def setError(self, event):

    value = self.errorEntry.get()
    if value is not None:
      if type(self.conditionPoint) is type(()):
        dataDim, index = self.conditionPoint
        pointErrors = dataDim.pointErrors
        if pointErrors:
          values = list(pointErrors)
        else:
          values = [0.0] * dataDim.numPoints
 
        values[index] = value
        dataDim.setPointErrors(values)
        
      else:  
        self.conditionPoint.setError( value )
  
  def getDetails(self, expSeries):
    
    if expSeries :
      self.detailsEntry.set(expSeries.details)
 
  def setDetails(self, event):

    text = self.detailsEntry.get()
    if text and text != ' ':
      self.expSeries.setDetails( text )

  def addExpSeries(self):
    
    expSeries = self.nmrProject.newNmrExpSeries(conditionNames=['delay time',])

  def deleteExpSeries(self, *event):
    
    if self.expSeries and (self.expSeries.className != 'Experiment'):
      if showOkCancel('Confirm','Really delete series?', parent=self):
        self.expSeries.delete()
        self.expSeries = None
        self.conditionPoint = None

  def getExperiments(self):
        
    if self.expSeries and (self.expSeries.className == 'Experiment'):
      return [self.expSeries,]
    
    else:
      return self.nmrProject.sortedExperiments()

  def getExperiment(self, sampleCondition):

    index = 0
    names = []
    
    if self.conditionPoint and (type(self.conditionPoint) != type(())):
      index = 0
      experiment = self.conditionPoint.parent.findFirstExperiment()
      experiments = self.getUnusedExperiments()
      
      name  = experiment.name
      names = [name] + [e.name for e in experiments]
      experiments = [experiment,] + experiments
    
    self.experimentPulldown.setup(names,experiments,index)
    
  def setExperiment(self, obj):
    
    experiment = self.experimentPulldown.getObject()
       
    if self.conditionPoint and (type(self.conditionPoint) != type(())):
      conditionSet = getExperimentConditionSet(experiment)
      
      if conditionSet is not self.conditionPoint.parent:
      
        if experiment not in self.expSeries.experiments:
          self.expSeries.addExperiment(experiment)
      
        unit = self.conditionPoint.unit
        if not unit:
          unit = CONDITION_UNITS_DICT[self.expSeries.conditionNames[0]]
        
        condition = self.conditionPoint.condition
      
        experiments = set(self.expSeries.experiments)
        for experiment0 in self.conditionPoint.parent.experiments:
          if experiment0 in experiments:
            experiments.remove(experiment0)
        self.expSeries.experiments = experiments
      
        value = self.conditionPoint.value
        error = self.conditionPoint.error
        self.conditionPoint.delete()
        
        self.conditionPoint = conditionSet.findFirstSampleCondition(condition=condition)
        if self.conditionPoint:
          self.conditionPoint.unit  = unit 
          self.updateAfter()          
        else:
          self.conditionPoint = conditionSet.newSampleCondition(condition=condition,unit=unit,
                                                                value=value,error=error)

  def updateDataDim(self, sampledDataDim):
  
    experiment = sampledDataDim.dataSource.experiment
    self.updateExperiments(experiment)

  def updateExperiments(self, experiment):
    
    experiments = self.getExperiments()
    names = [e.name for e in experiments]
    self.experimentPulldown.setup(names, experiments,0)
    
    if getExperimentSampledDim(experiment):
      self.updateAfter()
  
    elif self.expSeries:
      if self.expSeries.className == 'Experiment':
        if experiment is self.expSeries:
          self.updateConditionsAfter()

      elif experiment in self.expSeries.experiments:
        self.updateConditionsAfter()
 
 
  def updateConditionsAfter(self, sampleCondition=None):


    if self.waitingConditions:
      return
    
    if sampleCondition:
      experiments = sampleCondition.sampleConditionSet.experiments
      for experiment in experiments:
         if self.expSeries.className == 'Experiment':
           if experiment is self.expSeries:
             self.waitingConditions = True
             self.after_idle(self.updateConditions)
             break
      
         elif experiment in self.expSeries.experiments:
          self.waitingConditions = True
          self.after_idle(self.updateConditions)
          break

    else:
      self.waitingConditions = True
      self.after_idle(self.updateConditions)
      
  def updateConditions(self):
  
    self.updateButtons()
    objectList = []
    textMatrix = []
    colorMatrix = []
    nCols = len(self.conditionPointsMatrix.headingList)
    defaultColors = [None] * nCols
    if self.expSeries:
      if self.expSeries.className == 'Experiment':
        dataDim     = getExperimentSampledDim(self.expSeries)
        analysisDataDim = getAnalysisDataDim(dataDim)
        conditionVaried =  dataDim.conditionVaried
        expName     = self.expSeries.name
        unit        = dataDim.unit
        pointValues = dataDim.pointValues
        pointErrors = dataDim.pointErrors
        refPlane    = analysisDataDim.refSamplePlane
        
        for i in range(dataDim.numPoints):
          if i < len(pointErrors):
            error = pointErrors[i]
          else:
            error = None  
          
          pointText = ':%3d' % (i+1)
          if i == refPlane:
            datum      = ['* Ref Plane *',
                          expName+pointText,
                          None,
                          None,
                          None]
            colorMatrix.append(['#f08080'] * nCols)
          else:
            datum      = [conditionVaried ,
                          expName+pointText,
                          pointValues[i],
                          error,
                          unit]
            colorMatrix.append(defaultColors)
 

          textMatrix.append(datum)
          objectList.append((dataDim, i))
             
      else:
        condDict = getNmrExpSeriesSampleConditions(self.expSeries)
        conditionNames = self.expSeries.conditionNames
        
        for conditionName in conditionNames:
          for sampleCondition in condDict.get(conditionName, []):
 
            datum      = [sampleCondition.condition,
                          ' '.join([e.name for e in sampleCondition.parent.experiments]),
                          sampleCondition.value,
                          sampleCondition.error,
                          sampleCondition.unit]
 
            textMatrix.append(datum)
            objectList.append(sampleCondition)
            colorMatrix.append(defaultColors)

    self.conditionPointsMatrix.update(objectList=objectList,
                                      colorMatrix=colorMatrix,
                                      textMatrix=textMatrix)
    
    self.waitingConditions = 0

  def updateAfter(self, object=None):
  
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
 
  def updateButtons(self):
  
    if self.expSeries is None:
      self.seriesButtons.buttons[1].disable()
      self.seriesButtons.buttons[2].disable()
      self.conditionPointsButtons.buttons[0].disable()
      self.conditionPointsButtons.buttons[1].disable()
      self.conditionPointsButtons.buttons[2].disable()

    elif self.expSeries.className == 'Experiment':
      self.seriesButtons.buttons[1].disable()
      self.seriesButtons.buttons[2].disable()
      self.conditionPointsButtons.buttons[0].disable()
      self.conditionPointsButtons.buttons[1].disable()
      self.conditionPointsButtons.buttons[2].enable()
    
    else:
      self.seriesButtons.buttons[1].enable()
      self.seriesButtons.buttons[2].enable()
      self.conditionPointsButtons.buttons[0].enable()
      self.conditionPointsButtons.buttons[2].disable()
      
      if self.conditionPoint is None:
        self.conditionPointsButtons.buttons[1].disable()
      else:
        self.conditionPointsButtons.buttons[1].enable()
 
 
  def update(self):
    
    self.updateButtons()
    
    objectList = []
    textMatrix = []
    for experiment in getSampledDimExperiments(self.nmrProject):
      
      getExperimentConditionSet(experiment)
      sampledDim = getExperimentSampledDim(experiment)
      datum      = [None,
                    experiment.name,
                    1,
                    sampledDim.conditionVaried,
                    sampledDim.numPoints]
 
      textMatrix.append(datum)
      objectList.append(experiment)    
    
    for expSeries in self.nmrProject.sortedNmrExpSeries():

      experiments   = expSeries.experiments
      conditionSets = len([e.sampleConditionSet for e in experiments if e.sampleConditionSet])
      
      datum      = [expSeries.serial,
                    expSeries.name or ' ',
                    len(experiments),
                    ','.join(expSeries.conditionNames),
                    conditionSets]
 
      textMatrix.append(datum)
      objectList.append(expSeries)
       
    self.seriesMatrix.update(objectList=objectList, textMatrix=textMatrix)
    self.updateConditions()

    self.waiting = False
    
  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)


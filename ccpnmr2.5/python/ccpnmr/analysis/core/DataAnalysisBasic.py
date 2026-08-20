
LICENSE = """
======================COPYRIGHT/LICENSE START==========================

DataAnalysisBasic.py: Part of the CcpNmr Analysis program

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
import math

from ccpnmr.analysis.core.AssignmentBasic import propagatePeakAssignments, assignResToDim, isPeakAssigned
from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes, getSpectrumNoise, getExperimentSampledDim
from ccpnmr.analysis.core.PeakBasic import findSameAssignmentPeaks, findClosePeaks, getPeakDimPpm
from ccpnmr.analysis.core.PeakBasic import pickPeak, arePeaksAssignedSame
from ccpnmr.analysis.core.Util import getAnalysisDataDim

from memops.general.Implementation import ApiError

try:
  from memops.gui.MessageReporter import showWarning
except ImportError:
  from memops.universal.MessageReporter import showWarning

from ccp.general.Constants import chemShiftRefRatios

try:
  from memops.c.FitMethod import fitData
  from memops.c.FitMethod import bootstrapData
  from memops.c.FitMethod import runFit
  from memops.c.FitMethod import error as FitMethodError
except Exception, e:
  print 'Error, the DataAnalysisBasic module fitting functionality will not work, something is wrong with the C code.'
  print 'Exception:', e
  print 'Will continue without Analysis C fitting functionality'

class DataFitting:
  """
  Object to store data prior to the fitting of a function and hold the 
  data relating to the fit afterwards
  """

  def __init__(self, refObject=None, objects=None, fitFunction=1, 
               noiseLevel=1.0, dataX=None, dataY=None, guiParent=None,
               errorsX=None, errorsY=None, fitErrorFunction='covariance'):
    """
    .refObject - Some reference object, e.g. Peak from which assignments are derived
    .objects - The objects the data are extracted from, used for recalculation
    .fitFunction - The index of the fitting function, starts at 1
    .dataX - The X axis data to fit
    .dataY - The Y axis data to fit
    .errorsX - Errors for the X axis data
    .errorsY - Errors for the Y axis data
    .noiseLevel - Noise level for chi square analysis
    .fitErrorFunction - The name of the fitting error function
    .className - Attribute for type testing - type testing often appears in Analysis
    .parameters - The fit parameters
    .parameterErrors - The estimated errors in each of the fit parameters
    .derivedData - Some arbitrary data derived from the fit, e.g. relaxation rate or half life
    .derivedDataErrors - Errors in the above arbitrary data derived from the fit
    .fitError - The goodness of fit of the selected function to the input data
    .fittedY - The Y axis values of the best fitting function
    .guiParent - The grphical parent object to display warnings relative to 
    """
    
    self.refObject        = refObject
    self.objects          = objects or []
    self.fitFunction      = fitFunction
    self.fitErrorFunction = fitErrorFunction
    self.dataX            = dataX or []
    self.dataY            = dataY or []
    self.errorsX          = errorsX or []
    self.errorsY          = errorsY or []
    self.noiseLevel       = noiseLevel
    self.className        = 'FittedData'
    self.guiParent        = guiParent
    
    self.parameters        = []
    self.parameterErrors   = []
    self.derivedData       = None
    self.derivedDataErrors = None
    self.fitError          = None
    self.fittedY           = []
  
  def fit(self):
    """Calculate the fit given the specified input data.
    Returns True if fit went ok, False otherwise
    """
    
    self.resetFit()
   
    error = ''
  
    if not self.dataX:
      error += 'No X axis data in fit'
      
    if not self.dataY:
      error += 'No Y axis data in fit'      
    
    if not error:
      try:
        result = functionFitData(self)
        if result:
          return True
        else:
          return False
      except Exception, e:
        #showWarning('Function fit failure', e, parent=self.guiParent)
        print e
        raise
       
    else:
      #showWarning('Function fit failure',error, parent=self.guiParent)
      print error
      raise Exception(error)

    return False
  
  def resetFit(self):
    """Reset the ouput fit attributes to blank or null values
    """

    self.parameters        = []
    self.parameterErrors   = []
    self.derivedData       = None
    self.derivedDataErrors = None
    self.fitError          = None
    self.fittedY           = []
  
  
  def resetAll(self):
    """Reset all attributes, i.e. including input data, to blank or null values
    """
  
    self.refObject     = None
    self.objects       = []
    self.dataX         = []
    self.dataY         = []
    self.errorsX       = []
    self.errorsY       = []

    self.parameters        = []
    self.parameterErrors   = []
    self.derivedData       = None
    self.derivedDataErrors = None
    self.fitError          = None
    self.fittedY           = []
   
  
  def delete(self):
    """Delete the FittedData and clean up
    """
  
    del self
  

def getExperimentConditionSet(experiment):
  """Get or make a new SampleConditionSet for an NMR experiment
\n.. describe:: Input\n\nNmr.Experiment
\n.. describe:: Output\n\nNmr.SampleConditionSet
  """

  conditionSet = experiment.sampleConditionSet
  if conditionSet is None:
    conditionSet = experiment.topObject.newSampleConditionSet(experiments=[experiment,])

  return conditionSet

def getNmrExpSeriesSampleConditions(nmrExpSeries):
  """Get the Sample conditions that relate to an NmrExpSeries
\n.. describe:: Input\n\nNmr.NmrExpSeries
\n.. describe:: Output\n\nDict of SampleCondition Type : List of Nmr.SampleConditions
  """
  
  conditionDict = {}
  sampleConditionTypes = nmrExpSeries.conditionNames
  for experiment in nmrExpSeries.experiments:
    sampleConditionSet = experiment.sampleConditionSet
    if sampleConditionSet:
      for sampleCondition in sampleConditionSet.sampleConditions:
        condition = sampleCondition.condition
        if condition in sampleConditionTypes:
          conditionDict[condition] = conditionDict.get(condition, []) + [sampleCondition,]

  return conditionDict

def makeDataList(name, unit, values, errors=None, peaks=None, peakDims=None, measurements=None):
  """Makes an NMR data list e.g. to record chemical shift changes during a titration.
             Assumes for the moment that one derivation gives one datum.
\n.. describe:: Input\n\nWord, Word, List of Floats, List of Floats, List of List of Nmr.Peaks
             List of List of Nmr.PeakDims, List of List of Nmr.AstractMeasuremenrs 
\n.. describe:: Output\n\nNmr.DataList
  """

  N = len(values)
  assert len(peaks) == N

  if errors is None:
    errors = [0 for x in range(N)]

  nmrProject  = peaks[0].topObject
  dataList    = nmrProject.newDataList(name=name, unit=unit)
  
  if not errors:
    errors = [None] * N
  if not peaks:
    peaks = [None] * N
  if not peakDims:
    peakDims = [None] * N
  if not measurements:
    measurements = [None] * N

  for i in range(N):
    derivation = dataList.newDataDerivation(inputMeasurements=measurements[i],
                                            peakDims=peakDims[i],peaks=peaks[i])
    try:
      resonances = measurements[i].resonances
    except:
      resonances = [measurements[i].resonance,]
      
    datum = derivation.newDatum(resonances=resonances, value=values[i], error=errors[i]) # figOfMerit

  return dataList

def makeRatesList(rateType, experimentOrSeries, specificType, peakGroups, values, errors=None):
  """
  Makes an NMR rates list of a specified type (T1, T2 etc) with input values.
  If an experiment is passed in, rather than a series, the experiment has a sampled dimension
  
  .. describe:: Input

  String ('T1','T2','H/D exchange' or 'H/D protection'),
  Nmr.NmrExpSeries or Nmr.Experiment, String (Nmr.Measurement.coherenceType/protectionType)
  List of List of Nmr.Peaks (rate series), List of Floats (Nmr.AbstractMeasurement.values),
  List of Floats (Nmr.AbstractMeasurement.errors)
  
  .. describe:: Output

  Nmr.AbstractMeasurementList (depends on rateType)
  """

  N = len(values)
  assert len(peakGroups) == N

  if errors is None:
    errors = [0 for x in range(N)]

  nmrProject = experimentOrSeries.topObject
  
  if experimentOrSeries.className == 'Experiment':
    dataDim = getExperimentSampledDim(experimentOrSeries)
    experiments = [experimentOrSeries,]
    if not dataDim:
      msg = 'Experiment %s has no sampled dimension'
      showWarning('Failure',msg % experimentOrSeries.name)
      return
    if dataDim.conditionVaried != 'delay time':
      msg = 'Experiment %s has no delay time dimension'
      showWarning('Failure', msg % experimentOrSeries.name)
      return
      
    unit = dataDim.unit  
      
  else:
    conditions = getNmrExpSeriesSampleConditions(experimentOrSeries).get('delay time')
    experiments = list(experimentOrSeries.experiments)
    if not conditions:
      msg = 'NmrExpSeries %d has no associated delay time conditions'
      showWarning('Failure', msg % experimentOrSeries.serial)
      return
      
    unit    = conditions[0].unit
 
  if not unit:
    msg = 'No unit set for the experiment series time points. Please set and retry.'
    showWarning('Failure', msg)
    return
    
  
  sf      = None
  expDims = experiments[0].expDims
  
  for expDim in expDims:
    expDimRef = expDim.findFirstExpDimRef()
    if expDimRef and expDimRef.isotopeCodes:
      isotope = expDimRef.isotopeCodes[0]
      if chemShiftRefRatios.get(isotope) is None:
        continue
      sf = expDimRef.sf/chemShiftRefRatios[isotope]
      break

  if rateType in ('T1','T2','T1rho'):
    isotope0 = '15N'
  else:
    isotope0 = '1H'
    
  resonanceCheck = {}
  resonances  = []
  peakGroups2 = []
  values2     = []
  errors2     = []

  for i in range(N):
    if values[i] is None:
      continue
  
    peak      = peakGroups[i][0]
    resonance = None
    peakDim   = peak.sortedPeakDims()[0]
    for peakDim2 in peak.sortedPeakDims():
      expDimRef = peakDim2.dataDim.expDim.findFirstExpDimRef()
      if expDimRef and (isotope0 in expDimRef.isotopeCodes):
        peakDim = peakDim2
        break
        
    if peakDim.peakDimContribs:
      resonance = peakDim.findFirstPeakDimContrib().resonance

    if resonance is None:
      contrib = assignResToDim(peakDim)
      if contrib:
        resonance = contrib.resonance
      else:
        continue
    
    if resonanceCheck.get(resonance):
      print "Warning: Resonance %s repeated in peak groups for %s list" % (resonance,rateType)
      continue
    
    peaks = [p for p in peakGroups[i] if not p.isDeleted]  
    peakGroups2.append(peaks)
    values2.append(values[i])
    errors2.append(errors[i]) 
    resonances.append(resonance)
    resonanceCheck[resonance] = 1

  peakGroups = peakGroups2
  values     = values2
  errors     = errors2
   
  N = len(values)
  measurementList = None
  if rateType == 'T1':
    measurementList = nmrProject.newT1List(unit=unit, sf=sf,
                                           coherenceType=specificType,
                                           experiments=experiments)
    for i in range(N):
      measurement = measurementList.newT1(value=values[i], error=errors[i],
                                          resonance=resonances[i], peaks=peakGroups[i])
        
  elif rateType == 'T1rho':
    measurementList = nmrProject.newT1RhoList(unit=unit, sf=sf,
                                              coherenceType=specificType,
                                              experiments=experiments)
    # wb104: 5 Apr 2011: spelling typo in data model
    # so put code in that will also work when fixed
    if hasattr(measurementList, 'newT1Rho'):
      newObjFunc = measurementList.newT1Rho  # it is this
    else:
      newObjFunc = measurementList.newT1rho  # it should be this

    for i in range(N):
      measurement = newObjFunc(value=values[i], error=errors[i],
                                             resonance=resonances[i],
                                             peaks=peakGroups[i])
     
  elif rateType == 'T2':
    measurementList = nmrProject.newT2List(unit=unit, sf=sf,
                                           coherenceType=specificType,
                                           experiments=experiments)
    for i in range(N):
      measurement = measurementList.newT2(value=values[i], error=errors[i],
                                          resonance=resonances[i],
                                          peaks=peakGroups[i])
    
  elif rateType == 'H/D exchange':
    measurementList = nmrProject.newHExchRateList(unit=unit, experiments=experiments)
    for i in range(N):
      measurement = measurementList.newHExchRate(value=values[i], error=errors[i],
                                                 resonance=resonances[i],
                                                 peaks=peakGroups[i])
    
  elif rateType == 'H/D protection':
    measurementList = nmrProject.newHExchProtectionList(unit=unit,
                                                        protectionType=specificType,
                                                        experiments=experiments)
    for i in range(N):
      measurement = measurementList.newHExchProtection(value=values[i], error=errors[i],
                                                       resonance=resonances[i],
                                                       peaks=peakGroups[i])

  return measurementList

def isResidueInRange(residue, residueRanges, dataDim):
  """Determine if a residue is in a residue range
\n.. describe:: Input\n\nNmr.Residue, List of (Nmr.DataDim, Nmr.Chain, Integer
             (Nmr.Residue.seqCode - first), Integer
             (Nmr.Residue.seqCode - last)), Nmr.DataDim
\n.. describe:: Output\n\nBoolean
  """

  for (dataDims, chain, start, end) in residueRanges:
    if dataDim in dataDims:
      if residue.chain is chain:
        if residue.seqCode >= start:
          if residue.seqCode <= end:
            return True
  
  return False
  
  
def pairHnHaPeaks(peakList, residueRanges=None):
  """Generates pairs of Hn and Ha peaks from a peak list
             which may be restricted to optional residue ranges. 
\n.. describe:: Input\n\nNmr.PeakList. List of (Nmr.DataDim, Nmr.Chain,
             Integer (Nmr.Residue.seqCode - first)
\n.. describe:: Output\n\nList of (Nmr.Peak, Nmr.Peak)
  """
  # HnHa peaks must be asigned

  if not residueRanges:
    residueRanges = None
    
  pairs = []
  
  dimH  = 0
  dimHx = 1
  dimN  = 2
  
  amides = {}
  alphas = {}
  
  for peak in peakList.peaks:
    peakDims = peak.sortedPeakDims()
    if isPeakAssigned(peak, fully=1):
      resonanceN  = peakDims[dimN].findFirstPeakDimContrib().resonance
      resonanceH  = peakDims[dimH].findFirstPeakDimContrib().resonance
      resonanceHx = peakDims[dimHx].findFirstPeakDimContrib().resonance
      
      if not (resonanceN.resonanceSet):
        continue
      if not (resonanceH.resonaceSet):
        continue
      if not(resonanceHx.resonanceSet):
        continue         
      
      atomN  = resonanceN.resonanceSet.findFirstAtomSet().findFirstAtom()
      atomH  = resonanceH.resonanceSet.findFirstAtomSet().findFirstAtom()
      atomHx = resonanceHx.resonanceSet.findFirstAtomSet().findFirstAtom()
      
      residue = atomN.residue
      inRange = isResidueInRange(residue, residueRanges, peakDims[dimN].dataDim)
      if residueRanges and not inRange:
        continue

      if residue is not atomH.residue:
        continue
      if residue is not atomHx.residue:
        continue
      
      if (atomN.name == 'N') and (atomH.name == 'H'):
        if atomHx.name == 'H':
          amides[residue] = peak
	
	elif atomHx.name == 'HA':
	  alphas[residue] = peak

  for residue in amides.keys():
    if alphas.get(residue) is not None:
      pairs.append( [ amides[residue],alphas[residue] ])

  return pairs

def getPeakSampledDimIntensities(peaks, experiment, intensityType='volume'):
  """
  Get the intensities (of specified type) of the input peaks and
  list them with their corresponding  sampled condition values.
  E.g. List peak intensities with the sampled T1/T2 time.
  Often used to generate graphs.
  Note the input experiment must have a sampled data dim.
  
  .. describe:: Input
  
  List of Nmr.Peaks, Nmr.Experiment, String (Nmr.PeakIntensity.intensityType)
  
  .. describe:: Output
  
  Tuple of (List of Floats (Nmr.PeakIntensity.value), List of Floats (Nmr.ExpSeriesCondition.value))
  """

  dataDim = getExperimentSampledDim(experiment)
  if not dataDim:
    showWarning('Failure','Experiment %s has mo sampled data dimensions' % experiment.name)
    return
  
  referencePoint = getAnalysisDataDim(dataDim).refSamplePlane
  refIntensity = None
    
  pointValues = dataDim.pointValues
  pointErrors = dataDim.pointErrors
  nValues = len(pointValues)
  dim = dataDim.dim

  xData   = []
  yData   = []
  xWidths = []
  yWidths = []
  okPeaks = []
  missing = 0
  for peak in peaks:
    yVal = None
    xVal = None
    yErr = None
    xErr = None
    
    intensity = peak.findFirstPeakIntensity(intensityType=intensityType)
    peakDim = peak.findFirstPeakDim(dim=dim)
    index = int(round(peakDim.position))-1
    if index == referencePoint:
      refIntensity = intensity.value
      continue    

    if intensity:
      yVal = intensity.value
      yErr = intensity.error 
    else:
      missing += 1
      
    xVal    = pointValues[index]
    xErr    = 0.0
    
    if not (0 <= index < nValues):
      continue
    
    if index < len(pointErrors):
      xErr = pointErrors[index] or 0.0
  
    if (yVal is not None) and (xVal is not None):
      yData.append( yVal )
      xData.append( xVal )
      yWidths.append( yErr )
      xWidths.append( xErr )
      okPeaks.append( peak )
  
  n = 0.0
  mean = 0.0
  for yErr in yWidths:
    if yErr is not None:
      n += 1.0
      mean += yErr
  
  if n > 0.0:
    mean /= n
    for i in range(len(yWidths)):
      if yWidths[i] is None:
        yWidths[i] = mean
    
  else:
    for i in range(len(yWidths)):
      yWidths[i] = getSpectrumNoise(okPeaks[i].peakList.dataSource)
  
  if missing > 0:
    print "Warning: %d peak intensities (type %s) missing" % (missing,intensityType)
      
  return (xData,yData, xWidths, yWidths, refIntensity)


def matchSampledExperimentPeaks(refPeakList, experiment, tolerances,
                                pickPeaks=True, pickNonMaxima=False,
                                doAssign=True, progressBar=None,
                                noiseThreshold=None):
  """
  Matches the peaks of the experiment (with a sampled dimension) to
  a reference peak list according to the input tolerances. Generates
  a list of peaks that are grouped together by having overlapping
  positions. Useful precursor for determining rates form the peak
  intensities of groups. Optional arguments to specify whether new
  peaks can be picked, whether at non-maximal contours and whether to
  assign matched peaks - all relative to the reference peak list. A
  progress bar popup can be passed in. Peaks will not be searched for
  below the noiseThreshold, but may be picked at the ref position.
  
  .. describe:: Input
  
  Nmr.PeakList, Nmr.Experiment,
  List of Floats (tolerances),
  Boolean, Boolean, Boolean, Boolean 
        
  .. describe:: Output
  
  List of List of Nmr.Peaks
  """
    
  from ccpnmr.analysis.core.PeakBasic import getClosestPeak
  
  refSpec = refPeakList.dataSource
  
  refIsotopes = getSpectrumIsotopes(refSpec)
  groups = []
  
  dataDim = getExperimentSampledDim(experiment)
  if not dataDim:
    showWarning('Group Peaks Failure',
                'Experiment %s has no sampled dims' % (experiment.name))
    return
  
  numSamples = dataDim.numPoints
  dim = dataDim.dim
  for spectrum in experiment.dataSources:
    if spectrum.dataType == 'processed':
      peakList = spectrum.activePeakList or spectrum.findFirstPeakList()
      break

  if progressBar:
    progressBar.total = len(refPeakList.peaks)
  
  groupDict = {}
  refPeaks = refPeakList.peaks
  for peak in refPeaks:
    if groupDict.get(peak) is True:
      continue
  
    group = [peak]
    if doAssign:
      for peakDim in peak.peakDims:
        if not peakDim.peakDimContribs:
          assignResToDim(peakDim)
    
    samplePeaks = []
    for i in range(numSamples):
      samplePeaks.append([])
    
    for peak2 in findSameAssignmentPeaks(peak, peakList):
      peakDim = peak2.findFirstPeakDim(dim=dim)
      i = int(peakDim.position)-1
      samplePeaks[i].append(peak2)
     
    if [] in samplePeaks:
      for peak2 in findClosePeaks(peak, peakList,
                                  tolerances=tolerances,
                                  pickNewPeaks=pickPeaks,
                                  noiseThreshold=noiseThreshold):
        if (not doAssign) and (not pickNonMaxima):
          if isPeakAssigned(peak2):
            continue
        
        peakDim = peak2.findFirstPeakDim(dim=dim)
        i = int(peakDim.position)-1
        samplePeaks[i].append(peak2)
   
    i = 1 
    for peaks in samplePeaks:
      peak2 = None
      
      if peaks:
        if len(peaks) == 1:
          peak2 = peaks[0]
          
        else:
          for peak3 in peaks:
            if arePeaksAssignedSame(peak,peak3):
              peak2 = peak3
              break
 
          if not peak2:
            peak2 = getClosestPeak(peak,peaks,tolerances)
       
      else:
        if pickNonMaxima and pickPeaks:
          # Both below functions default to points for sampled dims
          expDims = experiment.sortedExpDims()
          position = [getPeakDimPpm(x) for x in peak.sortedPeakDims()]
          for j in range(len(expDims)):
            if not expDims[j].expDimRefs:
              break
          else:
            j = len(expDims) - 1
          if len(position) < len(expDims):
            position.insert(j, i)
          else:
            position[j] = i
            
          peak2 = pickPeak(peakList, position, unit='ppm')
   
      if peak2:
        group.append(peak2)
        groupDict[peak2] = True
    
      i +=1
      
    if group and doAssign:
      propagatePeakAssignments(group, refPeak=peak)
  
    groups.append( group )
    
    if progressBar:
      progressBar.increment()
  
  return groups


def matchSeriesPeaks(refPeakList, expSeries, tolerances, pickPeaks=True,
                     pickNonMaxima=False, doAssign=True,
                     progressBar=None, noiseThreshold=None):
  """
  Matches the peaks of the experiments in an NMR experiment series to
  a reference peak list according to the input tolerances. Generates
  a list of peaks that are grouped together by having overlapping
  positions. Useful precursor for determining rates from the peak
  intensities of groups. Optional arguments to specify whether new
  peaks can be picked, whether at non-maximal contours and whether to
  assign matched peaks - all relative to the reference peak list. A
  progress bar popup can be passed in.
  
  .. describe:: Input
  
  Nmr.PeakList, Nmr.NmrExpSeries, List of Floats (tolerances),
  Boolean, Boolean, Boolean 
                    
  .. describe:: Output
  
  List of List of Nmr.Peaks
  """
  
  from ccpnmr.analysis.core.PeakBasic import getClosestPeak
  refSpec = refPeakList.dataSource
  N = refSpec.numDim
  assert len(tolerances) == N
  
  dimMappings = []
  groups = []
  targetPeakLists = []
  experiments = list(expSeries.experiments)
  refIsotopes = getSpectrumIsotopes(refSpec)
  for experiment in experiments:
    for spectrum in experiment.dataSources:
      if spectrum.dataType != 'processed':
        continue
        
      isotopes = getSpectrumIsotopes(spectrum)
      dimMapping = None
 
      if isotopes == refIsotopes:
        peakList = spectrum.activePeakList or spectrum.findFirstPeakList()
      
      else:
        isotopes.reverse()
        
        if isotopes == refIsotopes:
          dimMapping = {} 
          for i in range(N):
            dimMapping[i] = N-i-1
          
          peakList = spectrum.activePeakList or spectrum.findFirstPeakList()
        else:
          continue
        
      if peakList:
        targetPeakLists.append((peakList, dimMapping))
      
  
  if progressBar:
    progressBar.total = len(refPeakList.peaks)
  
  for peak in refPeakList.peaks:
  
    group = [peak]
    if doAssign:
      for peakDim in peak.peakDims:
        if not peakDim.peakDimContribs:
          assignResToDim(peakDim)
   
    for peakList, dimMapping in targetPeakLists:
      if peakList is refPeakList:
        group.append(peak)
        continue
    
      peaks = findSameAssignmentPeaks(peak, peakList)
      
      if not peaks:
        
        if (not doAssign) and (not pickNonMaxima):
          peaks = []
          for peak3 in findClosePeaks(peak, peakList,
                                      tolerances=tolerances,
                                      pickNewPeaks=pickPeaks,
                                      dimMapping=dimMapping,
                                      noiseThreshold=noiseThreshold):
            if not isPeakAssigned(peak3):
              peaks.append(peak3)
                                         
        else:
          peaks =  findClosePeaks(peak, peakList,
                                  tolerances=tolerances,
                                  pickNewPeaks=pickPeaks,
                                  dimMapping=dimMapping,
                                  noiseThreshold=noiseThreshold)
           
      peak2 = None
      if peaks:
        if len(peaks) == 1:
          peak2 = peaks[0]
        else:
          sameAssn = []
          for peak3 in peaks:
            if arePeaksAssignedSame(peak,peak3):
              sameAssn.append(peak3)
              
          if sameAssn:
            peak2 = getClosestPeak(peak, sameAssn, tolerances, dimMapping)
              
          else:
            peak2 = getClosestPeak(peak, peaks, tolerances, dimMapping)
              
      elif pickNonMaxima and pickPeaks:
        refPeakDims = peak.sortedPeakDims()
        
        if dimMapping:
          indices = [dimMapping[x] for x in range(len(refPeakDims))]
          position = [getPeakDimPpm(refPeakDims[x]) for x in indices]
        else:
          position = [getPeakDimPpm(x) for x in refPeakDims] 
          
        peak2 = pickPeak(peakList, position, unit='ppm')
   
      if peak2:
        group.append(peak2)
          
    if group and doAssign:
      propagatePeakAssignments(group, refPeak=peak)
  
    groups.append( group )
    
    if progressBar:
      progressBar.increment()
  
  return groups


def getPeakSeriesIntensities(peaks, expSeries, intensityType='volume', sampleConditionType='delay time'):
  """Get the intensities (of specified type) of the input peaks and list them with their corresponding
             NMR experiment series condition values. E.g. List peak intensities with the pH of their experiment.
             Often used to generate graphs.
             If an experiment, rather than an NmrExpSeries is passed in the experiment has a sampled data dim.
\n.. describe:: Input\n\nList of Nmr.Peaks, Nmr.NmrExpSeries, String (Nmr.PeakIntensity.intensityType), String
\n.. describe:: Output\n\nTuple of (List of Floats (Nmr.PeakIntensity.value), List of Floats (Nmr.ExpSeriesCondition.value))
  """
  xData   = []
  yData   = []
  xWidths = []
  yWidths = []

  conditions = getNmrExpSeriesSampleConditions(expSeries).get(sampleConditionType)
  if not conditions:
    data = (sampleConditionType, expSeries.serial)
    showWarning('Failure','No sample conditions of type "%s" in NmrExpSeries %d' % data)
    return (xData, yData, xWidths, yWidths)
 
  okPeaks = []
  missing = 0
  for peak in peaks:
    yVal = None
    xVal = None
    yErr = None
    xErr = None
    intensity = peak.findFirstPeakIntensity(intensityType=intensityType)
    if intensity:
      yVal = intensity.value
      yErr = intensity.error 
    else:
      missing += 1
      
    experiment = peak.peakList.dataSource.experiment
    for condition in conditions:
      if experiment in condition.sampleConditionSet.experiments:
        xVal = condition.value
        xErr = condition.error or 0.0
        break
  
    if (yVal is not None) and (xVal is not None):
      yData.append( yVal )
      xData.append( xVal )
      yWidths.append( yErr )
      xWidths.append( xErr )
      okPeaks.append( peak )
  
  n = 0.0
  mean = 0.0
  for yErr in yWidths:
    if yErr is not None:
      n += 1.0
      mean += yErr
  
  if n > 0.0:
    mean /= n
    for i in range(len(yWidths)):
      if yWidths[i] is None:
        yWidths[i] = mean
    
  else:
    for i in range(len(yWidths)):
      yWidths[i] = getSpectrumNoise(okPeaks[i].peakList.dataSource)
  
  if missing > 0:
    print "Warning: %d peak intensities (type %s) missing" % (missing,intensityType)
      
  return (xData, yData, xWidths, yWidths)
  
def getFitErrorInfo():

  data = [ ('covariance', 'Calculates covariance (valid if have normally distributed errors)'),
           ('bootstrap',  'Bootstrap method (samples data with replacement, i.e. by allowing points to be repeated)'),
           ('jiggling',   'Jiggles data (samples data by moving points in both x and y)'),
         ]

  return data

def getFitErrorMethods():

  fitErrorMethods = [x[0] for x in getFitErrorInfo()]

  return fitErrorMethods

def getFitErrorMethod(analysisProject):

  fitErrorMethods = getFitErrorMethods()
  defaultErrorMethod = fitErrorMethods[0]

  application = analysisProject.root.application
  value = application.getValue(analysisProject,  keyword='fitErrorMethod', defaultValue=defaultErrorMethod)

  if value not in fitErrorMethods:
    value = defaultErrorMethod
  
  return value

def setFitErrorMethod(analysisProject, value):

  fitErrorMethods = getFitErrorMethods()
  if value not in fitErrorMethods:
    return

  application = analysisProject.root.application
  application.setValue(analysisProject,  keyword='fitErrorMethod', value=value)

  return value

def getFitMethodInfo():
  """Get the data: equation name, number of parameters, formatted formula string 
             from a numbered fitting function & which patameter, if any, is a "rate".
\n.. describe:: Input\n\nInteger
\n.. describe:: Output\n\nList of Tuple(String, Int, Tuple(String, Tuple(Ints)) )
  """
  
  data = [('Ax + B',                                 2, ('%4fx + %4f', (0,1))                                       ),
          ('log (A exp(-Bx))',                       2, ('log (%4f exp(-%4fx)))', (0,1))                            ),
          ('A exp(-Bx)',                             2, ('%4f exp(-%4fx)', (0,1))                                   ),
          ('A exp(-Bx) + C',                         3, ('%4f exp(-%4fx) + %4f', (0,1,2))                           ),
          ('A (1-sin(BX)/Bx) + C',                   3, ('%4f (1-sin(%4fx)/%4fx) + %4f', (0,1,1,2))                 ),
          ('Ax / (1+Ax)',                            1, ('%4fx / (1 + %4fx)', (0,0))                                ),
          ('A((B+4x-sqrt((B+4x)^2-(4x)^2))/4x - C)', 3, ('%4f((%4f+x-sqrt((%4f+x)^2-(4x)^2))/4x - %4f)', (0,1,1,2)) ),
          ('A(1/2 - exp(-Bx))',                      2, ('%4f(1/2- exp(-%4fx))', (0,1)) ),
          ('A(B+x-sqrt((B+x)^2-4x))', 2, ('%4f(%4f+x-sqrt((%4f+x)^2-4x))', (0,1,1)) ),
          ('A exp(-Bx^2)', 2, ('%4f exp(-%4fx^2)', (0,1)) ),
          ('A cos(Bx)', 2, ('%4f cos(%4fx)', (0,1)) ),
          ('CPMG fast, kAB=kBA', 3, ('CPMG fast, kAB=kBA', ()) ),
          ('CPMG slow, kAB=kBA', 3, ('CPMG slow, kAB=kBA', ()) ),
          ('CPMG fast', 4, ('CPMG fast', ()) ),
          ('CPMG slow', 4, ('CPMG slow', ()) ),
          ('C - A exp(-Bx)',     3, ('%4f - %4f exp(-%4fx)', (2,0,1)) ),
          ('A(1 - exp(-Bx))',    2, ('%4f(1 - exp(-%4fx)', (0,1)) ),
         ]
 
  return data
    
def functionFitData(dataFitting, numIterations=1000):
  """
  Tries to fit a function (type of which is specified by method) to X and Y
  values. Noise level can be added to estimate the error rate for the overall
  fit. A separate error is given for the fit parameters, i.e. for rate
  constants.
  
  .. describe:: Input
  
  ccpnmr.analysis.dataFitting object (non CCPN API), Float
  
  .. describe:: Output
  
  ccpnmr.analysis.dataFitting object (non CCPN API)
  """
  
  method = dataFitting.fitFunction
  if method <1:
    return None
  
  errorMethod = dataFitting.fitErrorFunction

  noise  = dataFitting.noiseLevel
  x      = dataFitting.dataX
  y      = dataFitting.dataY
  xW     = dataFitting.errorsX
  yW     = dataFitting.errorsY
  nIter  = numIterations
  
  
  if method == 10: # A exp(-B x^2)
    x = [v*v for v in x]
    method = 3 # A exp(-B x)
  
  try:
    if errorMethod == 'covariance':
      (params, paramsDev, yFit, chiSq) = runFit(method,noise,x,y)
    elif errorMethod == 'bootstrap':
      (params, paramsDev, yFit, chiSq) = bootstrapData(method,nIter,noise,x,y)
    else:
      (params, paramsDev, yFit, chiSq) = fitData(method,nIter,noise,x,y,xW,yW)
  except FitMethodError, e:
    msg  = 'Problem with fitting: method = %s, nIter = %s, '
    msg += 'noise = %s, x = %s, y = %s, xW = %s, yW = %s: %s'
    raise ApiError(msg % (method, nIter, noise, x, y, xW, yW, e))
  
  dataFitting.parameters      = params
  dataFitting.parameterErrors = paramsDev
  dataFitting.fittedY         = yFit
  dataFitting.fitError        = math.sqrt(abs(chiSq))
  
  return dataFitting

def matchHnoePeaks(assignPeakList, refPeakList, satPeakList, tolerancesA, tolerancesB,
                   pickNewPeaks=True, doAssignments=False):
  """Match reference and saturation peak list peaks to an assigned (potentially) peak list
             using input toleranes. Used to calaulate the NOE from peak intensity ratios.
\n.. describe:: Input\n\nNmr.PeakList, Nmr.PeakList, Nmr.PeakList, List of Floats, List of Floats
\n.. describe:: Output\n\nList of (Nmr.Peak, Nmr.Peak)
  """
  
  from ccpnmr.analysis.core.PeakBasic import getClosestPeak
  
  # kludge below: deals with two very closely
  # picked peaks in the Assign spectrum
  donePairs = {}
  peakPairs = []
 
  for peak in assignPeakList.peaks:
  
    if doAssignments:
      for peakDim in peak.peakDims:
        if not peakDim.peakDimContribs:
          assignResToDim(peakDim)

    peakDims = peak.sortedPeakDims()
    ppm0 = peakDims[0].value
    ppm1 = peakDims[1].value

    refPeaks = findSameAssignmentPeaks(peak, refPeakList)
    if not refPeaks:
      refPeaks = findClosePeaks(peak, refPeakList,
                                tolerances=tolerancesA,
                                pickNewPeaks=pickNewPeaks)
 
    if pickNewPeaks and not refPeaks:
      position = [ppm0,ppm1]
      refPeak = pickPeak(refPeakList, position, unit='ppm')
    elif not refPeaks:
      continue
    else:
      refPeak = refPeaks[0]
      if len(refPeaks) > 1:
        refPeak = getClosestPeak(peak,refPeaks,tolerancesA)

    (tol0,tol1) = tolerancesB

    satPeaks   = findSameAssignmentPeaks(peak, satPeakList)
    if not satPeaks:
      satPeaks = findClosePeaks(peak, satPeakList,
                                tolerances=tolerancesB,
                                pickNewPeaks=pickNewPeaks)
      
    if pickNewPeaks and not satPeaks:
      position = [ppm0,ppm1]
      satPeak = pickPeak(satPeakList, position, unit='ppm')
    elif not satPeaks:
      continue
    else:
      satPeak = satPeaks[0]
      if len(satPeaks) > 1:
        satPeak = getClosestPeak(peak,satPeaks,tolerancesB)
          
    if doAssignments:
      propagatePeakAssignments([refPeak,satPeak], refPeak=peak)

    if not donePairs.get( (satPeak,refPeak) ):
      peakPairs.append( (satPeak,refPeak) )
      donePairs[(satPeak,refPeak)] = 1

  return peakPairs


def getPeakMatchRegion(peak, tolerances):
  """Legacy code copy only. Function now in PeakBasic.
\n.. describe:: Input\n\nNmr.Peak, List of Floats (tolerances for each Nmr.PeakDim)
\n.. describe:: Output\n\nList or (Float, Float)
  """

  print 'Warning: Deprecated use of DataAnalysisBasic getPeakMatchRegion()\n Use PeakBasic getPeakMatchRegion() instead'
  from ccpnmr.analysis.core.PeakBasic import getPeakMatchRegion
  return getPeakMatchRegion(peak, tolerances)


def getClosestPeak(peak,peaks,tolerances):
  """Legacy code copy only. Function now in PeakBasic.
\n.. describe:: Input\n\nNmr.Peak, List of Nmr.Peaks, List of Floats
\n.. describe:: Output\n\nNmr.Peak
  """

  print 'Warning: Deprecated use of DataAnalysisBasic getClosestPeak()\n Use PeakBasic getClosestPeak() instead'
  from ccpnmr.analysis.core.PeakBasic import getClosestPeak
  return getClosestPeak(peak,peaks,tolerances)


def getOverlapScore(peakA,peakB,tolerances):
  """Legacy code copy only. Function now in PeakBasic.
\n.. describe:: Input\n\nNmr.Peak, Nmr.Peak, List of Floats (tolerances for each Nmr.PeakDim)
\n.. describe:: Output\n\nFloat
  """
   
  print 'Warning: Deprecated use of DataAnalysisBasic getOverlapScore()\n Use PeakBasic getPeaksOverlapScore() instead'
  from ccpnmr.analysis.core.PeakBasic import getPeaksOverlapScore
  return getPeaksOverlapScore(peakA,peakB,tolerances)

def calcT2List(t1List, t1rhoList, spinLock):
  """Make a new T2 list based on a T1 list and a T1rho list

     .. describe:: Input

     Nmr.T1List, Nmr.T1rhoList, float

     .. describe:: Output

     Nmr.T2List
  """
  
  spinLock2 = spinLock ** 2

  nmrProject = t1List.nmrProject
  sf = t1List.sf

  if abs(t1rhoList.sf - sf) > 0.001 * sf:
    msg = 'T1 sf = %.1f, T1rho sf = %.1f. Both T1 and T1rho list '
    msg += 'should have the same spectrometer frequency'
    showWarning('Error', msg % (sf, t1rhoList.sf))
    return
  
  t1Experiments  = set(t1List.experiments)
  
  if not t1Experiments:
    t1Experiments = set()
    for t1 in t1List.measurements:
      for peak in t1.peaks:
        t1Experiments.add(peak.peakList.dataSource.experiment)
    
  t1rhoExperiments = set(t1rhoList.experiments)

  if not t1rhoExperiments:
    t1rhoExperiments = set()
    for t1r in t1rhoList.measurements:
      for peak in t1r.peaks:
        t1rhoExperiments.add(peak.peakList.dataSource.experiment)
  
  experiments = t1Experiments.union(t1rhoExperiments)
  experiments = [e for e in experiments if e.shiftList]
  
  shiftLists = [expt.shiftList for expt in experiments]
  if not shiftLists:
    msg = 'Could not determine compatible shiftList for T1 and T1rho lists'
    showWarning('Error', msg)
    return
  
  spectra = [e.findFirstDataSource() for e in experiments if e.dataSources]
  if not spectra:
    msg = 'No spectrum to derive reference info from'
    showWarning('Error', msg)
    return

  experiment = experiments[0]
  spectrum = spectra[0]
  shiftList = shiftLists[0]

  t2Lists = nmrProject.findAllMeasurementLists(className='T2List')
  
  i = len(t2Lists)+1
  format = 'T2 List %d'
  name = format % (i)
  
  while nmrProject.findFirstMeasurementList(name=name):
    i += 1
    name = format % (i)
  
  details = 'Calculated from T1 and T1rho'
  t2List = nmrProject.newT2List(sf=sf, details=details, name=name, unit=t1List.unit)
      
  # the sf above is for 1H and here we want sf for the other isotope
  for dataDim in spectrum.dataDims:
    for dataDimRef in dataDim.dataDimRefs:
      expDimRef = dataDimRef.expDimRef
      if '1H' not in expDimRef.isotopeCodes:
        sf = expDimRef.sf
        break
    
    else:
      continue
    break    
  
  # calculate frequency of midpoint of spectrum
  # take into account fact that spectrum might have been chopped
  point = 1 + 0.5*dataDim.numPointsOrig - dataDim.pointOffset
  ref = dataDimRef.refValue - (point-dataDimRef.refPoint)*dataDimRef.valuePerPoint
    
  for t1 in t1List.measurements:
    resonance = t1.resonance
    t1rho = t1rhoList.findFirstMeasurement(resonance=resonance)
    
    if t1rho:
      shift = resonance.findFirstShift(parentList=shiftList)
      
      if shift:
        # now calculate t2
        offset = (shift.value - ref) * sf
        norm = math.sqrt(offset*offset + spinLock2)
        sin_t = offset / norm
        cos_t = spinLock / norm
        sin_t2 = sin_t * sin_t
        cos_t2 = cos_t * cos_t
        r2 = (1.0/t1rho.value - sin_t2/t1.value) / cos_t2
        value = 1.0 / r2
        if t1rho.error > 0 and t1.error > 0:
          error = cos_t2 / math.sqrt((1.0/t1rho.error)**2 + (sin_t2/t1.error)**2)
        else:
          error = 0.0
        t2 = t2List.newT2(value=value, error=error, resonance=resonance)
  
  return t2List

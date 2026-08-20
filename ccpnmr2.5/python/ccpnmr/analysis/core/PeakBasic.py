
LICENSE = """
======================COPYRIGHT/LICENSE START==========================

PeakBasic.py: Part of the CcpNmr Analysis program

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
import operator
# from Python 2.3 can use sum() instead of reduce(operator.add, )

from memops.universal.Util import cumulativeProductArray, arrayOfIndex

from memops.general.Implementation import ApiError
from memops.general.Util import copySubTree
from ccpnmr.analysis.core.AssignmentBasic import assignResToDim, makePeakDimAnnotation, clearPeakDim
from ccpnmr.analysis.core.AssignmentBasic import getPeakDimAtomTuple, isShiftInRange
from ccpnmr.analysis.core.ExperimentBasic import getEquivalentDataDims, getSpectrumIsotopes
from ccpnmr.analysis.core.ExperimentBasic import getDataDimRefFullRange, getThroughSpaceDataDims
from ccpnmr.analysis.core.ExperimentBasic import ppmDataDimBoundedRegion, ptsDataDimBoundedRegion
from ccpnmr.analysis.core.ExperimentBasic import getPrimaryDataDimRef, getOnebondDataDims
from ccpnmr.analysis.core.BoxIntegral     import truncatedBoxIntegral
from ccpnmr.analysis.core.MergeObjects    import mergeObjects
from ccpnmr.analysis.core.MoleculeBasic   import areAtomsTocsyLinked, getResidueObservableAtoms
from ccpnmr.analysis.core.MoleculeBasic   import getNumConnectingBonds, getLinkedResidue, DEFAULT_ISOTOPES
from ccpnmr.analysis.core.PeakFindParams  import getPeakFindParams, getPeakFindMinLinewidth
from ccpnmr.analysis.core.PeakFindParams  import getPeakFindBoxwidth
                                            
from ccpnmr.analysis.core.UnitConverter import unit_converter, pnt2ppm, ppm2pnt
from ccpnmr.analysis.core.Util import getAnalysisDataDim, getAnalysisSpectrum
from ccpnmr.analysis.core.Util import getMethod, convertPosition, convertRegion, \
                                      getIsotopeExclusion, \
                                      getDimChecked, getDimWrapped

from ccp.util.LabeledMolecule import getSchemeAtomPairFractions, getSchemeAtomFractions
from ccp.util.LabeledMolecule import getExperimentAtomPairFractions, getExperimentAtomFractions
from ccp.util.NmrExpPrototype import longRangeTransfers

try:
  from memops.gui.MessageReporter import showOkCancel, showWarning
  from memops.gui.DataEntry import askInteger
except ImportError:
  from memops.universal.MessageReporter import showOkCancel, showWarning
  from memops.universal.DataEntry import askInteger

from ccp.api.nmr.Nmr import SampledDataDim, PeakDimContribN

PEAK_FIT_METHODS = [
  'parabolic',
  'gaussian',
  'lorentzian',
]

PEAK_VOLUME_METHODS = [
  'box sum',
  'truncated box sum',
  'peak fit',
]

def uniteResonancePeakPositions(peak):
  """
  Coordinate all peaks in the peak list of the input peak so that
  all peak dimensions that are assigned to the same resonances as the 
  input peak have the same peak dimension position.

  .. describe:: Input
  
  ccp.nmr.Nmr.Peak
  .. describe:: Output
  
  ccp.nmr.Nmr.Peaks
  """

  peakList = peak.peakList
  peakDims = peak.sortedPeakDims()
  peaks    = set()
  
  resPositions = {}
  for peakDim in peakDims:
    resonances = set()
    
    for contrib in peakDim.peakDimContribs:
      resonance = contrib.resonance
      resonances.add(resonance)
      
      for contrib2 in resonance.peakDimContribs:
        peak2 = contrib2.peakDim.peak
      
        if peak2.peakList is peakList:
          peaks.add(contrib2.peakDim.peak)
  
    if resonances:
      resonances = tuple(resonances)
    
      if resPositions.get(resonances) is None:
        resPositions[resonances] = []
 
      resPositions[resonances].append(peakDim.realValue)

  for resonances in resPositions:
    values = resPositions[resonances]
    mean = sum(values) / float(len(values))
    resPositions[resonances] = mean

  for peak2 in peaks:
    if peak2 is peak:
      continue
  
    for peakDim in peak2.sortedPeakDims():
      resonances = set()
 
      for contrib in peakDim.peakDimContribs:
        resonance = contrib.resonance
        resonances.add(resonance)
 
      position = resPositions.get(tuple(resonances))
      if position is not None:
        peakDim.value = position
 
  return list(peaks)

def makeTransposePeakList(peakList):
  """Make a transpose peak list, where equivalent dimension positions are swapped,
             using peak positions form an existing peak list. The new peak list is created 
             in the same spectrum as the input. Assignments will be copied, and this will
             not affect the shift values. 

  .. describe:: Input
  
  ccp.nmr.Nmr.PeakList
  .. describe:: Output
  
  ccp.nmr.Nmr.PeakList
  """
  
  from ccpnmr.analysis.core.AssignmentBasic import propagatePeakAssignments
  
  spectrum    = peakList.dataSource
  
  tolerances  = [getAnalysisDataDim(dd).assignTolerance for dd in spectrum.sortedDataDims()]
  dimPairList = getEquivalentDataDims(spectrum)

  if not dimPairList:
    showWarning('Transpose Peak List Failure',
                'Spectrum has no equivalent dimension pairs')
    return

  dimPairs = {}
  for dataDimA, dataDimB in dimPairList:
    dimPairs[dataDimA] = dataDimB.dim
    dimPairs[dataDimB] = dataDimA.dim
  
  unit = spectrum.experiment.shiftList.unit
  newPeakList = spectrum.newPeakList()
  newPeakList.details = 'Synthetic transpose peak list'

  for peak in peakList.peaks:
    position = []

    for peakDim in peak.sortedPeakDims():
      dataDim  = peakDim.dataDim
      equivDim = dimPairs.get(dataDim)
      
      if equivDim is not None:
        peakDim2 = peak.findFirstPeakDim(dim=equivDim)
        position.append(peakDim2.value)
      
      elif peakDim.dataDimRef:
        position.append(peakDim.value)
        
      else:
        position.append(peakDim.position)

    newPeak = pickPeak(newPeakList, position, unit=unit)
    newPeak.details = peak.details

    propagatePeakAssignments([newPeak,], refPeak=peak, tolerances=tolerances)

  return newPeakList


def getIsotopeWeightedTolerances(peak, tolerance, scaleFactorDict=None):
  """Weight a shift distance tolerance according to the isotopes
             that correspond to a peak dimension, with the option to
             pass in the isotope scale factors.
             Also returns scale factors used for each dimension.

  .. describe:: Input
  
  Nmr.Peak,  Float (tolerance), Dict of Resonance.isotopeCode:Float (isotope: scale factor)
  .. describe:: Output
  
  List of Floats (tolerences), List of Floats (scale factors)
  """

  if not scaleFactorDict:
    scaleFactorDict = {'1H':1.0 ,'13C':10.0 ,'15N':5.0}

  tolerances   = []
  scaleFactors = []
  
  for peakDim in peak.sortedPeakDims():
    dataDimRef = peakDim.dataDimRef
    
    if dataDimRef:
      expDimRef  = dataDimRef.expDimRef
      factor     = 0.0
 
      i = 0.0
      for isotope in expDimRef.isotopeCodes:
        factor += scaleFactorDict.get(isotope, 1.0)
        i += 1.0
 
      if i:
        factor /= i
      else:
        factor = 1.0
 
      # the tolerence for a given dim is the shift
      # distance threshold * the isotope scale factor
 
      tolerances.append(tolerance*factor)
      scaleFactors.append(factor)

    else:
      tolerances.append(None)
      scaleFactors.append(1.0)
          

  return tolerances, scaleFactors
  
  
def findShiftDistPeakMatches(peak, peakList, threshold, scaleFactorDict=None,
                             considerAliased=True, dimMapping=None):
  """Find peaks within a given shift distance threshold of a peak and
             rank them by distance. A different threshold is uesd for each peak
             dim accoding to a dictonary of scale factors keyed by isotope.
             Option to specify dim mapping for source to target dims
             (starts at 0)

  .. describe:: Input
  
  Nmr.Peak, Nmr.PeakList, Floats (threshold),
             Dict of ExpdimRef.isotopeCode:Float, Dict of Int:Int
  .. describe:: Output
  
  List of 2-Lists [Float, Nmr.Peak] - (distance, peak)
  """

  sqrt = math.sqrt

  if not scaleFactorDict:
    scaleFactorDict = {'1H':1.0 ,'13C':10.0 ,'15N':5.0}

  matches = []
  
  tolerances, scaleFactors = getIsotopeWeightedTolerances(peak, threshold, scaleFactorDict)
      
  peaks = findClosePeaks(peak, peakList, tolerances,
                         considerAliased=considerAliased,
                         dimMapping=dimMapping)
  for peakB in peaks:
    value = getPeaksOverlapScore(peak,peakB,scaleFactors,
                                 considerAliased=considerAliased,
                                 dimMapping=dimMapping)

    dist = sqrt(value)
    if dist <= threshold:
      matches.append([dist, peakB])
  
  matches.sort()
  return matches


def getPeakMatchRegion(peak, tolerances, considerAliased=False):
  """Get a region (list of min max bound for each dimension) around a
             peak position with bounds set by input tolerances. The region will
             be twice the tolerance width and centered on the peakDim position.
             Option to also consider aliased/unaliased positions

  .. describe:: Input
  
  Nmr.Peak, List of Floats (tolerances for each Nmr.PeakDim), Boolean
  .. describe:: Output
  
  List of List of (Float, Float)
  """

  peakDims = peak.sortedPeakDims()
  regions = [[],]
  for i in range(len(peakDims)):
    peakDim =  peakDims[i]
    dataDimRef = peakDim.dataDimRef
    if dataDimRef:
      
      if considerAliased:
        ranges = [getDataDimRefFullRange(dataDimRef),]
        points = getAliasedPeakDimPositions(peakDim, ranges)
        
        newRegions = []
        
        for region in regions:
          for point in points:
            ppm = pnt2ppm(point,dataDimRef)
          
            newRegion = list(region)
            newRegion.append((ppm-tolerances[i],ppm+tolerances[i]))
            newRegions.append(newRegion)
     
        regions = newRegions
     
      else:  
        ppm = peakDim.value
        
        for region in regions:
          region.append((ppm-tolerances[i],ppm+tolerances[i]))

    else:
      dataDim = peakDim.dataDim
      for region in regions:
        region.append((1, dataDim.numPoints)) # Whole sampled region

  return regions

def findSameAssignmentPeaks(peak, peakList=None, onlyFullyAssigned=True):
  """
  Find a peak in a peak list (or all peak lists if none is specified)
  with the same assignment as the input peak, if it exists. Option to
  specifiy that the funtion only works for fully assigned peaks.
  
  .. describe:: Input
  
  Nmr.Peak, Nmr.PeakList

  .. describe:: Output

  List of Nmr.Peaks
  """
  
  # WARNING: 22 Jan 2015: this code did not work for diagonal peaks because it found any peak with one dim matching rather than two
  # Check for diagonal peaks but only code it for 2D for now (which should be the only case this function is being used for)
  
  i = 0
  resonances = [] 
  peakDims   = [pd for pd in peak.sortedPeakDims() if pd.dataDimRef]
   
  for peakDim in peakDims:
    dimResonances = []
    for contrib in peakDim.peakDimContribs:
      dimResonances.append(contrib.resonance)
    if dimResonances:
      i += 1
    resonances.append(dimResonances)

  if len(peakDims) == 2 and set(resonances[0]).intersection(set(resonances[1])):
    isDiagonal = True
  else:
    isDiagonal = False
    
  if i == 0:
    return []
  elif onlyFullyAssigned and (i != len(peakDims)):
    return []
  
  candidates = set([])
  
  if isDiagonal:
    resonances = set(resonances[0])
    for resonance in resonances:
      for contrib in resonance.peakDimContribs:
        peak2 = contrib.peakDim.peak
        if (not peakList) or (peak2.peakList is peakList):
          peakDims = [pd for pd in peak2.sortedPeakDims() if pd.dataDimRef]
          for peakDim in peakDims:
            if peakDim is not contrib.peakDim:
              dimResonances = [contrib2.resonance for contrib2 in peakDim.peakDimContribs]
              if resonance not in dimResonances:
                break
          else:
            candidates.add(peak2)
      
  else:
    for n, dimResonances in enumerate(resonances):
      peaks = set([])
    
      for resonance in dimResonances:
        for contrib in resonance.peakDimContribs:
          peak2 = contrib.peakDim.peak
          if (not peakList) or (peak2.peakList is peakList):
            peaks.add( contrib.peakDim.peak)
            
      if n > 0:
        candidates = candidates.intersection(peaks)
      else:
        candidates = peaks
        
      if not candidates:
        break

  return list(candidates)

def findClosePeaks(peak, peakList, tolerances=None, pickNewPeaks=False,
                   considerAliased=False, dimMapping=None, noiseThreshold=None):
  """
  
  Find peaks in a peak list within position tolerances of the input peak.
  Default tolerances are the assignment tolerances for each dimension
  The peak dimension distances are weighted according to the tolerances.
  Option to pick new peaks if no existing peaks are found.
  Option to look for peaks at aliased positions.
  Optional dict can be passed in to say how peak dims map together (list
  indices start at 0, key is for input peak).
  New peaks will not be picked below the noiseThreshold.
  
  .. describe:: Input
  
  Nmr.Peak, Nmr.PeakList, List of Floats, Boolean, Boolean, Boolean
  
  .. describe:: Output
  
  List of Nmr.Peaks
  """

  if not tolerances:
    tolerances = []
    for peakDim in peak.sortedPeakDims():
      dataDimRef = peakDim.dataDimRef
      if dataDimRef:
        analysisDataDim = getAnalysisDataDim(dataDimRef.dataDim)
        tolerances.append(analysisDataDim.assignTolerance)
      else:
        tolerances.append(None)  

  regions = getPeakMatchRegion(peak, tolerances, considerAliased=considerAliased)
  peaks   = set()
  
  if dimMapping:
    regions2 = []
    N = range(len(tolerances))

    for region in regions:
      region2 = []
      dataDims = peakList.dataSource.sortedDataDims()
      for (i, dataDim) in enumerate(dataDims):
        if dataDim.className == 'FreqDataDim':
          dataDimRef = getPrimaryDataDimRef(dataDim)
          rr = getDataDimRefFullRange(dataDimRef)
        else:
          rr = (1, dataDim.numPoints)
        region2.append(rr)
      for i in N:
        j = dimMapping.get(i)

        if j is not None:
          region2[j] = region[i]

      regions2.append(region2)      
 
    regions = regions2

  for region in regions:
    peaks.update( searchPeaks([peakList],region,considerAliased=considerAliased) )
 
    if pickNewPeaks:
      peaks.update( findPeaks(peakList, region, intensityThreshold=noiseThreshold) )

  if peak in peaks:
    peaks.remove(peak)
  
  return list(peaks)

def getPeaksOverlapScore(peakA,peakB,tolerances=None, considerAliased=False,dimMapping=None):
  """Get a score for how well overlapped two peaks are. Scores for each peak dimension are normalised with respect to
             input dimension tolerances (often DataDim assignment tolerance)
             Default tolerances are the assignment tolerances for each dimension
             Option to score peaks when there may be an aliasing mismatch.
             Optional dict can be passed in to say how peak dims map together (list
             indices start at 0).

  .. describe:: Input
  
  Nmr.Peak, Nmr.Peak, List of Floats (tolerances for each Nmr.PeakDim), Boolean, Dict of Int:Int
  .. describe:: Output
  
  Float
  """
  
  if not tolerances:
    tolerances = []
    for peakDim in peakA.sortedPeakDims():
      dataDimRef = peakDim.dataDimRef
      if dataDimRef:
        analysisDataDim = getAnalysisDataDim(dataDimRef.dataDim)
        tolerances.append(analysisDataDim.assignTolerance)
      else:
        tolerances.append(None)

  score = 0.0
  pdA = peakA.sortedPeakDims()
  pdB = peakB.sortedPeakDims()

  N = len(pdA)
  if not dimMapping:
    dimMapping = {}

    for i in range(N):
      dimMapping[i] = i

  for i in range(N):
    j = dimMapping.get(i)

    if j is None:
      continue 

    tolerance = tolerances[i]
    peakDimA  = pdA[i]
    peakDimB  = pdB[j]
    
    dataDimRefA = peakDimA.dataDimRef
    dataDimRefB = peakDimB.dataDimRef
    if not (dataDimRefA and dataDimRefB):
      continue
    
    if tolerance is not None:
      ppmB = peakDimB.value
      diff = abs(peakDimA.value-ppmB)
      
      if considerAliased:
        fullRange  = getDataDimRefFullRange(dataDimRefA)
        
        for point in getAliasedPeakDimPositions(peakDimA, [fullRange,]):
          ppm = pnt2ppm(point, dataDimRefA)
          d2  = abs(ppm-ppmB)
          
          if d2 < diff:
            diff = d2
      
      diff /= tolerance
      
      score += diff*diff
    
  return score


def getClosestPeak(peak,peaks,tolerances,dimMapping=None):
  """
  Find the closest peak to the input peak from a list of peaks.
  The peak dimensions are weighted according to input tolerances.
  Optional dict can be passed in to say how peak dims map together (list
  indices start at 0).
  
  .. describe:: Input
  
  Nmr.Peak, List of Nmr.Peaks, List of Floats, Dict of Int:Int
  
  .. describe:: Output
  
  Nmr.Peak
  """

  score = peakMin = None

  for peak2 in peaks:
    if peak2 is peak:
      continue
    s = getPeaksOverlapScore(peak,peak2,tolerances,dimMapping=dimMapping)
    if not peakMin or s < score:
      score = s
      peakMin = peak2

  return peakMin
   
def findSymmetryPeaks(peak, tolerances=None, peakLists=None, considerAliased=True):
  """Find the symmetry related transpose peaks (e.g. NOE return) 
             within tolerances given an input peak.
             Option to search a given set of peak lists the
             default is otherwise the input peaks list.

  .. describe:: Input
  
  List of Nmr.Peak
  .. describe:: Output
  
  List of Nmr.Peaks
  """
  from AssignmentBasic import getPeakDimFullShiftRange

  dimPairs = {}
  peakDims = peak.sortedPeakDims()
  N = len(peakDims)
  
  if not peakLists:
    peakLists = [peak.peakList,]
  
  if not tolerances:
    tolerances = [getAnalysisDataDim(pd.dataDimRef.dataDim).assignTolerance for pd in peakDims]
  
  for i in range(N-1):
    peakDim  = peakDims[i]
    dataDimRef = peakDim.dataDimRef
    
    if dataDimRef:
      isotopes = list(peakDim.dataDimRef.expDimRef.isotopeCodes)
      isotopes.sort()
      
      for j in range(i+1,N):
        peakDim2 = peakDims[j]
        dataDimRef2 = peakDim2.dataDimRef
        
        if dataDimRef2:
          isotopes2 = list(dataDimRef2.expDimRef.isotopeCodes)
          isotopes2.sort()
          if isotopes == isotopes2:
            dimPairs[peakDim] = peakDim2
            dimPairs[peakDim2] = peakDim

  peaks = []
  
  if dimPairs:
    region = []
    center = []
    i = 0
    for peakDim in peakDims:
      tol = tolerances[i]
 
      pairedDim = dimPairs.get(peakDim)
      if pairedDim:
        pos = pairedDim.value
        dimRange = [pos-tol,pos+tol]
        center.append(pos)
 
      else:
        dimRange = getPeakDimFullShiftRange(peakDim)
        center.append(None)
 
      region.append(dimRange)
      i += 1

    peaks = searchPeaks(peakLists, region, considerAliased=considerAliased)

    if peaks:
      rankedPeaks = []
      for peakB in peaks:
        score = 0.0
        for i in range(N):
          pos = center[i]
          if pos is not None:
            d = abs( peakB.sortedPeakDims()[i].value - pos )/tolerances[i]
            score += d*d
      
        rankedPeaks.append((score, peakB))
  
      rankedPeaks.sort()
      peaks = [x[1] for x in rankedPeaks]
  
  return peaks

def findPositionMatches(peaks, dimDict, peakLists=None, matchTolerance=None,
                        separateSpinSystems=False):
  """
  Search the target peaklist or the peak lists of the input peaks for other
  groups of peaks which have matches to the positions of the input peaks in the
  given dimension. The search dimension is given but the input dictionary, keyed
  by peak list. Returns a ranked list of group positions in order of strength of
  match. Option to pass in a tolerance for the match dimension. Option to
  separate groups according to spinSystem assignment.
  
  .. describe:: Input
  
  List of Nmr.Peaks, Dict of (NmrPeakList:Int) Int, Nmr.PeakLists,
  Float, Boolean
  
  .. describe:: Output
  
  Tuple of (Float, List of Floats (positions), List of Nmr.Peaks)
  """
  from AssignmentBasic import getPeakDimFullShiftRange, getDataDimFullShiftRange
  
  
  if not peakLists:
    peakLists = []
    
  queries = []
  mRegions = []
  peakListDict = {}
  
  for peakList in peakLists:
    peakListDict[peakList] = 1
  
  if len(peaks) < 1:
    return
  
  K = len(peaks[0].peakDims)
  
  for peak in peaks:
    if len(peak.peakDims) != K:
      continue
  
    dimNum = dimDict[peak.peakList] -1
    for d, peakDim in enumerate(peak.sortedPeakDims()):
      dataDimRef = peakDim.dataDimRef
      
      if dataDimRef:
        if d == dimNum:
          if matchTolerance:
            tol = matchTolerance
          else:
            tol = getAnalysisDataDim(dataDimRef.dataDim).assignTolerance
          value = peakDim.value
          queries.append(value)
          mRegions.append([value-tol,value+tol])
 
  fullRegions = {}
  for peakList in peakLists:
    
    region = []
    spectrum = peakList.dataSource
    for dataDim in spectrum.sortedDataDims():
      region.append(getDataDimFullShiftRange(dataDim))
  
    fullRegions[peakList] = region
  
  # should get the initial peak => pos control for testing
  N = len(mRegions)
  resultGroups = []
  spinSystems = {}
  c = 0
  for mRegion in mRegions:
    
    resultGroup = []
    for peakList in peakLists:
      region = fullRegions[peakList][:]
      region[dimDict[peakList]-1] = mRegion
      resultGroup += searchPeaks([peakList,], region)
    
    resultGroups.append(resultGroup)
    for peak in resultGroup:
      refPos = []
      tols = []
      spinSystem = set()
      dimNum = dimDict[peak.peakList] -1
      
      for d, peakDim in enumerate(peak.sortedPeakDims()):
        dataDimRef = peakDim.dataDimRef
      
        if dataDimRef:
          analysisDataDim = getAnalysisDataDim(dataDimRef.dataDim)
        
          if d != dimNum:
            refPos.append(peakDim.value)
            tols.append(analysisDataDim.assignTolerance)
            
            for contrib in peakDim.peakDimContribs:
              # None is OK
              spinSystem.add(contrib.resonance.resonanceGroup)
                
          else:
            peak.qtol = analysisDataDim.assignTolerance
            peak.qdim = peakDim

      peak.column = c
      peak.refPos = refPos
      peak.tols = tols
      peak.columnList = [peak,]
      
      spinSystems[peak] = spinSystem
      
      c += 1

  Km1 = K-1
  for i in range(N-1):
    for peakI in resultGroups[i]:
      for j in range(i,N):
        for peakJ in resultGroups[j]:
          if peakI.column == peakJ.column:
            continue

          if separateSpinSystems and (spinSystems[peakI] != spinSystems[peakJ]):
            continue
        
          for k in range(Km1):
            if abs(peakI.refPos[k] - peakJ.refPos[k]) > peakJ.tols[k]:
              break
          
          else:
            for peak in peakJ.columnList:
              for k in range(Km1):
                if abs(peakI.refPos[k] - peak.refPos[k]) < peakI.tols[k]:
                  if peak not in peakI.columnList:
                    peakI.columnList.append(peak)
            
            peakJ.columnList = peakI.columnList
            peakJ.column = peakI.column

  done = {}
  columnGroups = []
  for group in resultGroups:
    for peak in group:
      if done.get(peak.column) is None:
        columnGroups.append( peak.columnList )
        done[peak.column] = 1
  
  sqrt = math.sqrt      
  groups = []
  for group in columnGroups:
  
    N = float(len(group))
    if N < 1:
      continue
    
    score = 0
    sums  = []
    for i in range(Km1):
      sums.append(0.0)
    
    for peak in group:
      
      for i in range(Km1):
        sums[i] += peak.refPos[i]
    
      qtol = peak.qtol
      for value in queries:
        delta = abs(value-peak.qdim.value) 
        if delta < qtol:
          d = (qtol-delta)/qtol
          score += d*d

    position = []
    for i in range(Km1):
      position.append(sums[i]/N)
    
    score = sqrt(score)
    groups.append( (score,position,group) )

  groups.sort()
  groups.reverse()
  
  return groups

def doPeaksOverlap(peakA, peakB):
  """Determine if two peaks overlap accordint to their data dim
             assignment tolerance.

  .. describe:: Input
  
  Nmr.Peak, Nmr.Peak
  .. describe:: Output
  
  Boolean
  """

  N = len(peakA.peakDims)
  M = 0
  
  for i in range(N):
    peakDimA = peakA.sortedPeakDims()[i]
    peakDimB = peakB.sortedPeakDims()[i]
    tol = getAnalysisDataDim(peakDimA.dataDimRef.dataDim).assignTolerance
    if abs(peakDimA.value-peakDimB.value) < tol:
      M += 1
    
  if M == N:
    return True
  else:
    return False

def getAliasedPeakDimPositions(peakDim, shiftRanges, returnPpms=False):
  """Give all the aliased/unaliased positions of a peakDim either in a
             specified shift range or the full range for the dimension type.
             Units for the shift ranges are ppm. Note this function uses the
             actual spectrum locations, rather than coupling center. Use the
             AssignmentBasic version to consider the center of multiplets.

  .. describe:: Input
  
  Nmr.PeakDim, List of (Tuples of Floats (MinShift,MaxShift) ), Boolean 
  .. describe:: Output
  
  List of Floats (Nmr.PeakDim.positions)
  """

  dataDimRef = peakDim.dataDimRef
  positions = []

  if dataDimRef:
    sw = peakDim.dataDimRef.dataDim.numPointsOrig
    peakDimPos = peakDim.position + (peakDim.numAliasing*sw)


    points = peakDimPos
    while isShiftInRange( pnt2ppm(points,dataDimRef), shiftRanges):
      if returnPpms:
        positions.append(pnt2ppm(points,dataDimRef))
      else:
        positions.append( points )
      points -= sw

    points = peakDimPos+sw
    while isShiftInRange( pnt2ppm(points,dataDimRef), shiftRanges):
      if returnPpms:
        positions.append(pnt2ppm(points,dataDimRef))
      else:
        positions.append( points )
      points += sw
  
  else:
    positions.append(peakDim.position)

  return positions

def makePeakAnnotation(peak):
  """Sets the annotation for a peak that indicates peak molSystem, merit and details
             (not the annotation of peakDims)

  .. describe:: Input
  
  Nmr.Peak
  .. describe:: Output
  
  None
  """
  assert peak
  project = peak.root
  analysisProject = project.currentAnalysisProject
  
  symbol     = ''
  detail     = ''
  molSysCode = ''
  
  if analysisProject.doMeritAnnotations:
    if peak.figOfMerit > 0.66:
      symbol = analysisProject.meritAnnotationGood or ''
    elif peak.figOfMerit < 0.34:
      symbol = analysisProject.meritAnnotationBad or '!'
    else:
      symbol = analysisProject.meritAnnotationMediocre or '~'
    symbol += ' '
 
  if analysisProject.doDetailAnnotations and (peak.details is not None):
    detail =  peak.details

  if analysisProject.doMolSysAnnotations:
    for peakDim in peak.peakDims:
      molSystem = None
      for contrib in peakDim.peakDimContribs:
        if contrib.resonance.resonanceSet:
          atom = contrib.resonance.resonanceSet.findFirstAtomSet().findFirstAtom()
          molSystem = atom.topObject
          break
      if molSystem:
        molSysCode = '%s:' % (molSystem.code)
        break
  
  compress = True # analysisProject.doCompressAnnotations
  if compress:
    
    peakDimAnnos = [getPeakDimAtomTuple(pd) for pd in peak.sortedPeakDims()]
    
    chain0, residue0, atom0 = peakDimAnnos[0]

    compress = True
    if len(peakDimAnnos) > 1:
      for chain, residue, atom in peakDimAnnos[1:]:
        if atom:
          if (chain != chain0) or (residue != residue0):
            compress = False
        
    if compress:
      assignAnno = '%s%s' % (chain0,residue0) + ','.join(text[2] or '-' for text in peakDimAnnos)
    else:
      assignAnno = ' '.join('%s%s%s' % texts for texts in peakDimAnnos)
 
    peakAnnotation = '%s%s%s %s ' % (symbol,molSysCode,assignAnno,detail)  
  
  else:
    peakAnnotation = '%s%s %s ' % (symbol,molSysCode,detail)
  
  peak.setAnnotation(peakAnnotation[:80])  


def getSpinSystemLinksLabel(peak):
  """Makes a setring describing the spin system connectivities of a peak.

  .. describe:: Input
  
  Nmr.Peak
  .. describe:: Output
  
  Line
  """
  
  seqLinks = []
  i = 1
  for peakDim in peak.sortedPeakDims():
    offsets = {}
   
    for contrib in peakDim.peakDimContribs:
      spinSystem = contrib.resonance.resonanceGroup
      if spinSystem:
      
        for link in spinSystem.findAllResonanceGroupProbs(linkType='sequential',
                                                          isSelected=True):
          delta = str(link.sequenceOffset)
          if delta[0] != '-':
            delta = '+' + delta
 
          offsets[delta] = None
 
        for link in spinSystem.findAllFromResonanceGroups(linkType='sequential',
                                                          isSelected=True):
          delta = str(-1 * link.sequenceOffset)
          if delta[0] != '-':
            delta = '+' + delta
 
          offsets[delta] = None
    
    if offsets:
      seqLinks.append( 'F%d:%s' % (i,','.join(offsets.keys())) )
        
    i += 1  

  if seqLinks:
    seqLinkText = 'Link %s' % (' '.join(seqLinks))
  else:
    seqLinkText = None

  return seqLinkText


def refreshPeakAnnotations(project):
  """Refreshes all of the peak and peakDim annotations for all peaks in a project

  .. describe:: Input
  
  Project
  .. describe:: Output
  
  None
  """
  assert project
  
  for resonance in project.currentNmrProject.resonances:
    if hasattr(resonance, 'guiName'):
      del resonance.guiName
    if hasattr(resonance, 'label'):
      del resonance.label
    
  for experiment in project.currentNmrProject.experiments:
    for dataSource in experiment.dataSources:
      for peakList in dataSource.peakLists:
        for peak in peakList.peaks:
          makePeakAnnotation(peak)
          
          for peakDim in peak.peakDims:
            makePeakDimAnnotation(peakDim)

def getDataDimMapping(dataSource, targetDataSource):
  """Get a mapping of the mapping of the dimensions of one
             spectrum onto another, where they match.

  .. describe:: Input
  
  Nmr.DataSource, Nmr.DataSource
  .. describe:: Output
  
  Dict of Int:Int (DataDim.dim:DataDim.dim)
  """

  mapping = {}
  
  dataDimsA = dataSource.dataDims
  dataDimsB = list(targetDataSource.dataDims)
  for dataDimA in dataDimsA:
    keysA = []
    
    if isinstance(dataDimA, SampledDataDim):
      continue
    
    for dataDimRef in dataDimA.dataDimRefs:
      expDimRef = dataDimRef.expDimRef
      isotopeCodes = expDimRef.isotopeCodes
      measurementType = expDimRef.measurementType.lower()
      #isAcquisition = expDimRef.expDim.isAcquisition
      
      #keysA.append( (isotopeCodes, measurementType, isAcquisition) )
      keysA.append( (isotopeCodes, measurementType) )
    
    for dataDimB in dataDimsB:
      keysB = []
      
      if isinstance(dataDimB, SampledDataDim):
        continue
      
      for dataDimRef in dataDimB.dataDimRefs:
        expDimRef = dataDimRef.expDimRef
        isotopeCodesB = expDimRef.isotopeCodes
        measurementTypeB = expDimRef.measurementType.lower()
        #isAcquisitionB = expDimRef.expDim.isAcquisition
      
        keysB.append( (isotopeCodesB, measurementTypeB) )

      if keysA == keysB:
        dataDimsB.remove(dataDimB)
        mapping[dataDimA.dim] = dataDimB.dim 
        break

  return mapping

def getDataDimRefMapping(dataSource, targetDataSource):
  """Get a mapping of the dataDimRefs for a data source to the equivalent
             dataDimRefs in a target data source. Assumes spectra with same types
             and number of dims.

  .. describe:: Input
  
  Nmr.DataSource, Nmr.DataSource
  .. describe:: Output
  
  Dict of Nmr.DataDimRef:Nmr.DataDimRef
  """

  objectMap = {}
  for i in range(dataSource.numDim):
    dataDimA = dataSource.sortedDataDims()[i]
    dataDimB = targetDataSource.sortedDataDims()[i]
    
    if isinstance(dataDimA, SampledDataDim):
      continue
    if isinstance(dataDimB, SampledDataDim):
      continue
        
    for dataDimRefA in dataDimA.dataDimRefs:
      expDimRefA      = dataDimRefA.expDimRef
      isotopeCodes    = expDimRefA.isotopeCodes
      measurementType = expDimRefA.measurementType
      isAcquisition   = expDimRefA.expDim.isAcquisition
      dataDimRefB     = None
      
      # Could consider ScalingFactors
      
      for dataDimRef in dataDimB.dataDimRefs:
        expDimRefB = dataDimRef.expDimRef
        if expDimRefB.isotopeCodes == isotopeCodes:
          if expDimRefB.measurementType == measurementType:
            if expDimRefB.expDim.isAcquisition == isAcquisition:
              dataDimRefB = dataDimRef
              break
      
      if dataDimRefB:
        objectMap[dataDimRefA] = dataDimRefB

  return objectMap

def copyPeakListNew(peakList, destDataSource, peakListParams=None):
  """
  Copies an entire peak list to another (or the same) DataSource)
  There are no appropriateness tests (so far)
  NB this routine *will* handle subpeaks correctly
  NB spectra must have same referencing otherwise peak points will
  not make sense 
             
  .. describe:: Input
  
  Nmr.PeakList, Nmr.DataSource, dictionary
  peakListParams is a dictionary of parameters that will be passed
  to the copy of the peaklist (for e.g. changing its name)
  
  .. describe:: Output

  Nmr.PeakList (new peak list)
  """
  
  objectMap  = getDataDimRefMapping(peakList.dataSource, destDataSource)
   
  # Maybe put appropriateness tests here
  
  newPeakList = copySubTree(peakList, destDataSource, 
                            topObjectParameters=peakListParams,
                            objectMap=objectMap)

  return newPeakList


def copyPeakList(peakList, destPeakList, progressBar=None, rePick=False, peaks=None):
  """
  Copies peaks from one peak list to another if appropriate, i.e dimension isotopes match
  Has the option to re-pick peaks, rather than doing a straight copy,
  NB this routine will *not* handle subpeaks correctly.
  Optional argument to copy only a subset of peaks from the source peak list.
  
  .. describe:: Input
  
  Nmr.PeakList, Nmr.PeakList, ProgressBar, Boolean, List of Nmr.Peaks
  
  .. describe:: Output
  
  None
  """
  assert peakList, destPeakList

  # are the peakLists equivalent
  
  if progressBar:
    guiParent = progressBar.parent
  else:
    guiParent = None
  
  spec1 = peakList.dataSource
  spec2 = destPeakList.dataSource
  
  if peaks:
    checkedPeaks = [p for p in peaks if p.peakList is peakList]
    
    if len(checkedPeaks) != len(peaks):
      msg = 'Some of the peaks specified were not in the source peak list; these were ignored'
      showWarning('Warning', msg, parent=guiParent)
    
    peaks = checkedPeaks
    
  else:
    peaks = peakList.sortedPeaks()
  
  dimMapping = getDataDimMapping(spec1, spec2)
  
  dimMappingRev = {}
  for dim1 in dimMapping:
    dim2 = dimMapping[dim1]
    dimMappingRev[dim2] = dim1

  for dataDim in spec1.dataDims:
    if not isinstance(dataDim, SampledDataDim):
      if dataDim.dim not in dimMapping:
        msg = "Could not map dimension %d of source's spectrum" % dataDim.dim
        showWarning('Error', msg, parent=guiParent)
        return

  for dataDim in spec2.dataDims:
    if not isinstance(dataDim, SampledDataDim):
      if dataDim.dim not in dimMappingRev:
        msg = "Could not map dimension %d of target's spectrum" % dataDim.dim
        showWarning('Error', msg, parent=guiParent)
        return
  
  nDim1 = spec1.numDim
  nDim2 = spec2.numDim
  
  if not rePick and (nDim1 != nDim2):
    msg  = 'Source and destination spectra must have'
    msg += 'same number of dimensions for cloning'
    showWarning('Error', msg, parent=guiParent)
    return

  samplePlane = None
  sampledDim  = None
  if nDim1 == nDim2:
    sampledDim = spec2.findFirstDataDim(className='SampledDataDim')
  
  elif nDim1 < nDim2:
    sampledDim = spec2.findFirstDataDim(className='SampledDataDim')
    
    if sampledDim:
      samplePlane = 0
      msg = 'Destination plane 1-%d or 0 => all' % sampledDim.numPoints
      samplePlane = askInteger('Enter Number', msg, 0, parent=guiParent)

      if (samplePlane is None) or (samplePlane>sampledDim.numPoints):
        showWarning('Error', 'Invalid plane selection', parent=guiParent)
        return

  elif nDim1 > nDim2:
    sampledDim = spec1.findFirstDataDim(className='SampledDataDim')
    
    if sampledDim:
      samplePlane = 0
      msg = 'Source plane 1-%d or 0 => all' % sampledDim.numPoints
      samplePlane = askInteger('Enter Number', msg,0, parent=guiParent)
  
      if (samplePlane is None) or (samplePlane>sampledDim.numPoints):
        showWarning('Error', 'Invalid plane selection', parent=guiParent)
        return
  
  if rePick:
    for peak in peaks:
      # Positions for new peak
      position = [None] * nDim2
      peakDims = peak.sortedPeakDims()

      for peakDim in peakDims:
        if peakDim.dataDimRef:
          # Normal, mapped freq dims
          i = dimMapping[peakDim.dim]-1
          position[i] = peakDim.value
        
        elif samplePlane is None:
          # Source and destination both have sampled dim
          plane = peakDim.position
          
          if plane > sampledDim.numPoints:
            # Source not in destination range
            break
          
          position[sampledDim.dim-1] = plane
      
        elif samplePlane and (peakDim.position != samplePlane):
          # Source only sampled
          # Peak is only used if plane choice matches
          break
      
      else:
        if samplePlane and (sampledDim.dataSource is spec2):
          # Destination uses chosen sampled dim
          position[sampledDim.dim-1] = samplePlane
 
        if (samplePlane == 0) and (sampledDim.dataSource is spec2):
          # Destination sampled, pick in all planes
          d = sampledDim.dim-1
          newPeaks = []
          for samplePoint in range(1,sampledDim.numPoints+1):
            position[d] = samplePoint
            newPeak = pickPeak(destPeakList, position, 'ppm', doFit=False)
            newPeaks.append(newPeak)

        else:
          newPeaks = [pickPeak(destPeakList, position, 'ppm', doFit=False),]
        
        for newPeak in newPeaks:
          checkContourRegion(newPeak)
          newPeakDims = newPeak.sortedPeakDims()
 
          for peakDim in peakDims:
            i = dimMapping.get(peakDim.dim)
            
            if i is None:
              # Can't assign to sampled
              continue
          
            for contrib in peakDim.peakDimContribs:
              resonance = contrib.resonance
              assignResToDim(newPeakDims[i-1], resonance, tolerance=1.0) 

          # add the details here
          newPeak.details = peak.details

      if progressBar:
        progressBar.increment()
        
  else:
    
    failedPeaks = []
    dataDimRefMapping = getDataDimRefMapping(spec1, spec2)
    
    for peak in peaks:
      try:
        newPeak = copySubTree(peak,destPeakList,objectMap=dataDimRefMapping)
      except Exception:
        failedPeaks.append(peak)
        
      if progressBar:
        progressBar.increment()
    
    if failedPeaks:
      showWarning('Error','%s peaks could not be copied' % len(failedPeaks))

def copyPeaksToPeakList(peaks, peakList):
  """
  Copy peaks from one peakList to another if they have the same spectrum.
  
  .. describe:: Input
  
  Nmr.Peak, Nmr.PeakList
  
  .. describe:: Output
  
  None
  """
  
  for peak in peaks:
    dataSource = peak.peakList.dataSource
    if dataSource is peakList.dataSource:
      newPeak = copySubTree(peak, peakList)
  
def _havePeakNearPosition(values, tolerances, peaks):

  for peak in peaks:
    for i, peakDim in enumerate(peak.sortedPeakDims()):
      if abs(peakDim.value - values[i]) > tolerances[i]:
        break
    else:
      return peak

  return None

def _copyPeakInfo(peakList, peak):

  peak2 = peakList.newPeak()
  peakDims2 = peak2.sortedPeakDims()
  for i, peakDim in enumerate(peak.sortedPeakDims()):
    for attr in ('boxWidth', 'decayRate', 'decayRateError', 'lineWidth', 'numAliasing', 'phase', 'phaseError', 'value', 'valueError', 'dataDimRef'):
      setattr(peakDims2[i], attr, getattr(peakDim, attr))
  
  for peakIntensity in peak.sortedPeakIntensities():
    peakIntensity2 = peak2.newPeakIntensity(method=peakIntensity.method)
    for attr in ('error', 'intensityType', 'value'):
      setattr(peakIntensity2, attr, getattr(peakIntensity, attr))

  peak2.details = peak.details

  return peak2

def subtractPeakLists(peakList1, peakList2, progressBar=None):
  """
  Subtracts peaks in peakList2 from peaks in peakList1, based on position,
  and puts those in a new peakList3.  Assumes a common spectrum for now.
  
  .. describe:: Input
  
  Nmr.PeakList, Nmr.PeakList, ProgressBar
  
  .. describe:: Output
  
  Nmr.PeakList
  """

  spectrum = peakList1.dataSource

  assert spectrum is peakList2.dataSource, 'For now requires both peak lists to be in same spectrum'

  dataDims = spectrum.sortedDataDims()
  tolerances = [getAnalysisDataDim(dataDim).assignTolerance for dataDim in dataDims]
  
  peaks2 = peakList2.peaks
  peakList3 = spectrum.newPeakList()

  for peak1 in peakList1.sortedPeaks():
    values1 = [peakDim.value for peakDim in peak1.sortedPeakDims()]
    if not _havePeakNearPosition(values1, tolerances, peaks2):
      _copyPeakInfo(peakList3, peak1)

def makeIntermediatePeak(peaks):
  """Make a peak which is intermediate in position to
             the input peaks

  .. describe:: Input
  
  List of Nmr.Peaks
  
  .. describe:: Output
  
  Nmr.Peak
  """

  if len(peaks) < 2:
    return

  N = 0
  
  peakList = peaks[0].peakList
  M = len(peaks[0].peakDims)
  position = [0] * M
  for peak in peaks:
    if peak.peakList is not peakList:
      showWarning('Error',
                  'Selected peak from different peak list. Skipping.')
      continue

    for i in range(M):
      position[i] += peak.sortedPeakDims()[i].position
    
    N +=1

  if N < 1:
    return

  for i in range(M):
    position[i] /= N
   
  return pickPeak(peakList, position)


def mergePeaks(peaks):
  """Merge two peaks together. * Function needs work/checking * Not really used in Analysis currently

  .. describe:: Input
  
  List of Nmr.Peaks
  .. describe:: Output
  
  Nmr.Peak
  """

  if len(peaks) < 2:
    return
  else:
    if not showOkCancel('Confirm', 'Merge %d peaks?' % len(peaks)):
      return

  targetPeak = peaks[0]
  peakList = targetPeak.peakList
  for peak in peaks[1:]:
    if peak.peakList is not peakList:
      showWarning('Error','Merging peak from different peak list. Skipping.')
      continue
    
    # Children to be merged?
    # peakContribs 
    # peakDims
    # peakIntensities
    mergeObjects(peak, targetPeak)


def deletePeak(peaks):
  """Delete given peaks with a warning

  .. describe:: Input
  
  List of Nmr.Peaks
  .. describe:: Output
  
  Boolean (success)
  """

  numPeaks = len(peaks) 
  warning = "Delete"

  if numPeaks > 1:
    warning = "%s %d crosspeaks?" % ( warning, len(peaks) )
    
  elif numPeaks == 1:
    warning = warning+" crosspeak?"
    
  else:
    return False

  if showOkCancel('Delete Peak', warning):
    for peak in peaks:
      if not peak.isDeleted:
        peak.delete()
    return True
    
  else:
    return False
    

def aliasedPeakDimPosition(peakDim):
  """Give the aliased points position of a given peak dimension

  .. describe:: Input
  
  Nmr.PeakDim
  .. describe:: Output
  
  Float (Nmr.PeakDim.position points) 
  """
  
  dataDimRef = peakDim.dataDimRef
  position   = peakDim.position
  
  if dataDimRef:
    position += peakDim.numAliasing*dataDimRef.dataDim.numPointsOrig
  
  return position 


def movePeak(peak, posn):
  """Move a peak to a given position. Automatically calculates any peak aliasing

  .. describe:: Input
  
  Nmr.Peak, List of Floats (Nmr.PeakDim.postion in points)
  .. describe:: Output
  
  None 
  """

  for i, peakDim in enumerate(peak.sortedPeakDims()):
    dataDimRef  = peakDim.dataDimRef
    
    if dataDimRef: # Not sampled
      peakDim.numAliasing = findNumAliasing(dataDimRef, posn[i])
      peakDim.position = posn[i] - peakDim.numAliasing*peakDim.dataDimRef.dataDim.numPointsOrig

  checkContourRegion(peak)
  setupPeakCoupling(peak)


def setPeakDimNumAliasing(peakDim, numAliasing, doWarning=True):
  """Set the number of peakDim aliasings to set a peakDim position and
             put it at its non-folded, correct ppm

  .. describe:: Input
  
  Nmr.PeakDim, Integer (peakDim.numAliasing) 
  .. describe:: Output
  
  None
  """
   
  dataDimRef = peakDim.dataDimRef
  if not dataDimRef:
    return
  
  expDimRef = dataDimRef.expDimRef
  origNumAliasing = peakDim.numAliasing 
  numAliasing = int(numAliasing)
  peakDim.setNumAliasing( numAliasing )
  
  minAliasedFreq, maxAliasedFreq = getDataDimRefFullRange(dataDimRef)
  
  ppm = peakDim.value
  if (ppm < minAliasedFreq) or (ppm > maxAliasedFreq):
  
    msg  = 'Value puts peak outside current spectrum bounds'
    msg += '(%.2f - %.2f %s). Continue and extend the bounds to %.2f?'
    data = (minAliasedFreq,maxAliasedFreq,expDimRef.unit, peakDim.value)
    if not doWarning or showOkCancel('Query', msg % data):

      analysisProject = peakDim.root.currentAnalysisProject
      if analysisProject and analysisProject.contourToUnaliased:
        expDimRef.minAliasedFreq = minAliasedFreq
        expDimRef.maxAliasedFreq = maxAliasedFreq
    else:
      peakDim.setNumAliasing( origNumAliasing )
      return
  
  resonances = []
  if peakDim.peakDimContribs:
    shiftList = expDimRef.expDim.experiment.shiftList
    for contrib in peakDim.peakDimContribs:
      resonances.append( contrib.resonance )
    
    for resonance in resonances:
      for contrib in resonance.peakDimContribs:
        peakDim2 = contrib.peakDim
        if peakDim2 is peakDim:
          continue
          
        position = unit_converter[(shiftList.unit,'point')](ppm,peakDim2.dataDimRef)
        numPointsOrig = float(peakDim2.dataDimRef.dataDim.numPointsOrig)
        n = int(round((position-aliasedPeakDimPosition(peakDim2))/numPointsOrig))
        
        if abs(n) > 0:
          for contrib2 in peakDim2.peakDimContribs:
            if contrib2.resonance not in resonances:
              clearPeakDim(peakDim2,contrib=contrib2)
              
          peakDim2.setNumAliasing( n + peakDim2.numAliasing )
          checkContourRegion(peakDim2.peak)

  checkContourRegion(peakDim.peak)


def propagatePeakUnaliasing(refPeak, peaks):
  """Gives the same positional unaliasing (as the input reference peak)
             to a group of peaks.

  .. describe:: Input
  
  Nmr.Peak, List of Nmr.Peaks
  .. describe:: Output
  
  None
  """
  
  if len (peaks) <  2:
    return
  
  refSpec = refPeak.peakList.dataSource
  N = len(refPeak.peakDims)
  for peak in peaks:
    if peak.peakList.dataSource is not refSpec:
      showWarning('Failure',
                  'Cannot set unaliasing for peaks from different spectra')
      return
  
  numAliasings = []
  for peakDim in refPeak.sortedPeakDims():
    numAliasings.append( peakDim.numAliasing )
      
  for peak in peaks:
    peakDims = peak.sortedPeakDims()
    for i in range(N):
      setPeakDimNumAliasing(peakDims[i], numAliasings[i])
      
      
def findMatchingPeaks(peakList, peakDim):
  """Find peaks within a peak list that are within the data dim assignment tolerance

  .. describe:: Input
  
  Nmr.PeakList, Nmr.PeakDim
  .. describe:: Output
  
  List of Nmr.Peaks 
  """

  peaks = []
  assert peakList and peakDim
  
  shiftList = peakList.dataSource.experiment.shiftList

  if shiftList:
    unit = shiftList.unit
  else:
    unit = 'ppm'
  
  dataDimRef = peakDim.dataDimRef
  if not dataDimRef:
    return
  
  isotopeCodes = dataDimRef.expDimRef.isotopeCodes
  dataDim   = dataDimRef.dataDim
  tolerance = getAnalysisDataDim(dataDim).assignTolerance
  if unit == 'point':
    positionA = aliasedPeakDimPosition(peakDim)
  else:
    positionA = unit_converter[('point',unit)](aliasedPeakDimPosition(peakDim),peakDim.dataDimRef)
    
  # find dimensions with matching isotopeCodes
  dims = []
  for i in range(peakList.dataSource.numDim):
    dataDim = peakList.dataSource.sortedDataDims()[i]
    if getPrimaryDataDimRef(dataDim).expDimRef.isotopeCodes == isotopeCodes:
      dims.append(i)
  
  # find peaks close in appropriate dimensions
  for peak in peakList.peaks:
    for i in dims:
      peakDim2 = peak.sortedPeakDims()[i]
      p2 = peakDim2.position + (peakDim2.numAliasing*peakDim2.dataDimRef.dataDim.numPointsOrig)
      if unit == 'point':
        positionB = p2
      else:
        positionB = unit_converter[('point',unit)](p2,peakDim2.dataDimRef)
      if abs(positionA - positionB) < tolerance:
        peaks.append(peak)
   
  return peaks

def findNumAliasing(dataDimRef, p):
  """Find the number if aliased spectral widths at a given point for a reference data dimension 

  .. describe:: Input
  
  Nmr.DataDimRef, Float
  .. describe:: Output
  
  Int (number of aliasings)
  """

  o = dataDimRef.dataDim.numPointsOrig
  #(n, r) = divmod(p-1, o)
  n = divmod(p-1, o)[0]
  return int(n)


def removePeaksAssignment(peaks):
  """Clear all peak dimension contribs for a list of peaks

  .. describe:: Input
  
  List of Nmr.Peaks
  .. describe:: Output
  
  None
  """

  for peak in peaks:
    for peakDim in peak.peakDims:
      clearPeakDim(peakDim)


def pickPeak(peakList, position, unit='point', doFit=True, figOfMerit=1.0, serial=None):
  """Pick a peak at a given position in given units in a peak list. 

  .. describe:: Input
  
  Nmr.PeakList, List of Floats (Nmr.PeakDim.positions), Word (Nmr.ShiftList.unit), Boolean
  .. describe:: Output
  
  Nmr.Peak
  """

  spectrum = peakList.dataSource
  if serial:
    peak = peakList.newPeak(figOfMerit=figOfMerit, serial=serial)
  else:
    peak = peakList.newPeak(figOfMerit=figOfMerit)

  dataDims = spectrum.sortedDataDims()
  peakDims = peak.sortedPeakDims()

  for i, peakDim in enumerate(peakDims):
    dataDim = dataDims[i] 
    
    if isinstance(dataDim, SampledDataDim):
      dataDimRef = None
    else:
      dataDimRef = getPrimaryDataDimRef(dataDim)

    if dataDimRef:
      if unit == 'point':
        p = position[i]
      else:
        p = dataDimRef.valueToPoint(position[i])
 
      peakDim.numAliasing = findNumAliasing(dataDimRef, p)
      peakDim.position = p - peakDim.numAliasing * dataDim.numPointsOrig

    else:
      peakDim.position = position[i]
    
  setupPeak(peak, doFit=doFit)
  checkContourRegion(peak)

  return peak


# TBD: peak routines assume position ordered as dim, and peakDim in same order as spectrum.dataDim
def addPeak(peakList, position, tile=None, parent=None, doFit=True, doOtherFit=True):
  """Add a peak to a peak list at a position in a specified GUI tile location

  .. describe:: Input
  
  Nmr.PeakList, List of Floats (Nmr.PeakDim.positions),
             List of Ints (Nmr.PeakDim.numAliasing), Tkinter.widget (GUI parent)
             Bool (whether position fit), Bool (whether height/linewidth fit)
  .. describe:: Output
  
  Nmr.Peak 
  """
 
  if not tile:
    tile = peakList.dataSource.numDim*[0]
 
  peak = peakList.newPeak()
  setPeakPosition(peak, position, tile)
  setupPeak(peak, doFit=doFit, doOtherFit=True)
  addPeakToSelected(peak, parent)
  checkContourRegion(peak)

  return peak


def arePeaksAssignedSame(peakA, peakB):
  """
  Determine whether peaks carry the same resonance assignments.
  PeakDimContrib order unimportant. If one peak has a lower dimensionality
  than the other then the best overlapping mapping of dims is considered:
  e.g. 3D and 2D NOESY peaks have same asignment if common 1H dims carry same
  resonances. Reciprocal peaks; e.g. with resonances (r1, r2) and (r2, r1)
  respectively on their dims, are considered to be assigned the same.

  .. describe:: Input
  
  Nmr.Peak, Nmr.Peak
  
  .. describe:: Output
  
  Boolean
  """
 
  peakDimsA = [pd for pd in peakA.sortedPeakDims() if pd.dataDimRef]
  peakDimsB = [pd for pd in peakB.sortedPeakDims() if pd.dataDimRef]
 
  if len(peakDimsA) > len(peakDimsB):
    peakDimsA, peakDimsB = peakDimsB, peakDimsA

  serialListB = []
  for peakDim in peakDimsB:
    serialsB = [ x.resonance.serial for x in peakDim.peakDimContribs]
    serialsB.sort()
    serialListB.append(serialsB)
  
  for peakDim in peakDimsA:
    serialsA = [ x.resonance.serial for x in peakDim.peakDimContribs ]
    serialsA.sort()
    for serialsB in serialListB:
      if serialsA == serialsB:
        serialListB.remove(serialsB)
        break

    else: # Found no dim that carries same resonances
      return False

  return True


def getPeakDimPpm(peakDim):
  """**Often depricated given peakDim.value** Get the position of a peak dimension in PPM

  .. describe:: Input
  
  Nmr.PeakDim
  .. describe:: Output
  
  Float (ppm)
  """

  dataDimRef = peakDim.dataDimRef
  if dataDimRef:
    return pnt2ppm( aliasedPeakDimPosition(peakDim), peakDim.dataDimRef )

  else:
    return peakDim.position


def checkContourRegion(peak):
  """Ensure the contours of a peaks spectrum extend to cover the input peak's position.

  .. describe:: Input
  
  Nmr.Peak
  .. describe:: Output
  
  None
  """

  analysisProject = peak.root.currentAnalysisProject
  if analysisProject and analysisProject.contourToUnaliased:
    for peakDim in peak.peakDims:
 
      if peakDim.dataDimRef:
        maxAliasedFreq = peakDim.dataDimRef.expDimRef.maxAliasedFreq
        minAliasedFreq = peakDim.dataDimRef.expDimRef.minAliasedFreq
        position       = peakDim.value
        if (minAliasedFreq is not None) and (position < minAliasedFreq):
          peakDim.dataDimRef.expDimRef.minAliasedFreq = position

        if (maxAliasedFreq is not None) and (position > maxAliasedFreq):
          peakDim.dataDimRef.expDimRef.maxAliasedFreq = position


def setPeakPosition(peak, position, tile):
  """Define the position of a peak and its aliasing tile

  .. describe:: Input
  
  Nmr.Peak, List of Floats (Nmr.PeakDim.position), List of Ints (Nmr.PeakDim.numAliasing)
  .. describe:: Output
  
  None
  """
 
  spectrum = peak.peakList.dataSource
  peakDims = peak.sortedPeakDims()
  for dim, peakDim in enumerate(peakDims):
    peakDim.numAliasing = tile[dim]
    if peakDim.position != position[dim]:
      peakDim.position = position[dim]

  checkContourRegion(peak)


def getPeakDimPosition(peakDim, toUnit = 'point'):
  """Give the position of peak dimension in specified units

  .. describe:: Input
  
  Nmr.PeakDim, Word (ShiftList.unit)
  .. describe:: Output
  
  Float
  """

  dataDimRef = peakDim.dataDimRef
  
  if not dataDimRef:
    # this might not have been set if peak imported via Format Converter (for example)
    dataDim = peakDim.dataDim
    if dataDim.className == 'FreqDataDim':
      peakDim.dataDimRef = dataDimRef = getPrimaryDataDimRef(peakDim.dataDim)
     
  if dataDimRef:
    #dataDim = dataDimRef.dataDim
    # TBD: looks like dataDim.pointOffset should not need be added here
    p = aliasedPeakDimPosition(peakDim)
    if (toUnit != 'point'):
      p = convertPosition(p, dataDimRef, toUnit=toUnit)
  
  else:
    if toUnit == 'point':
      p = peakDim.position
    else:
      p = peakDim.dataDim.pointValues[int(peakDim.position)-1]

  return p


def setupPeak(peak, doFit=True, doOtherFit=True):
  """Initialise reference data dimensions, peak position and intensities for a new peak.

  .. describe:: Input
  
  Nmr.Peak (new)
  Bool (whether position fit), Bool (whether height/linewidth fit)
  .. describe:: Output
  
  None
  """

  setupPeakDataDimRef(peak)
  setupPeakFit(peak, doFit=doFit, doOtherFit=doOtherFit)
  if doOtherFit:
    setupPeakVolume(peak)
  setupPeakCoupling(peak)

def setupPeakCoupling(peak):
  """
  Setup the underlying real peakDim positional values for MQ spectra etc.
  
  .. describe:: Input
  
  Nmr.Peak 
  
  .. describe:: Output
  
  None
  """
  
  for peakDim in peak.peakDims:
    dataDimRef = peakDim.dataDimRef
    if dataDimRef and (dataDimRef.expDimRef.measurementType == 'MQShift'):
      spectrum = peak.peakList.dataSource
      ppm = None
      
      for dataDimA in spectrum.dataDims:
        dataDimRefA =  getPrimaryDataDimRef(dataDimA)
        expDimRefA = dataDimRefA.expDimRef
        # TBD: will only use first match
        if (expDimRefA.measurementType == 'Shift') and \
          (expDimRefA.isotopes == dataDimRef.expDimRef.isotopes):
          expDimRefs = set([expDimRefA, dataDimRef.expDimRef])
          expTransfer = spectrum.experiment.findFirstExpTransfer(expDimRefs=expDimRefs)
          
          if expTransfer and expTransfer.transferType == 'Jcoupling':
            peakDimA = peak.findFirstPeakDim(dim=dataDimA.dim)
            ppm = peakDim.value - peakDimA.value
          
      if ppm is not None:
        peakDim.realValue = ppm 


def setupPeakDataDimRef(peak):
  """Gives a peak a default data dim reference if it has none

  .. describe:: Input
  
  Nmr.Peak 
  .. describe:: Output
  
  None
  """

  for peakDim in peak.peakDims:
    dataDimRef = peakDim.dataDimRef
    
    if (not dataDimRef) and (peakDim.dataDim.className == 'FreqDataDim'):
      dataDimRef = getPrimaryDataDimRef(peakDim.dataDim)
      
      if dataDimRef is not None:
        # TBD: more general NBNB may be OK now? Rasmus Oct 08
        peakDim.dataDimRef = dataDimRef
        
    # Cope with MQ axes
    if dataDimRef:
      expDimRef = dataDimRef.expDimRef
      
      if (expDimRef.measurementType == 'MQShift') and (expDimRef.refExpDimRef):
        for scalingFactor in expDimRef.refExpDimRef.validScalingFactors:
          # TBD: Int is horrid!!
          peakDim.newPeakDimComponent(dataDimRef=dataDimRef, \
                                      scalingFactor=int(scalingFactor))
        
 
def setupPeakFit(peak, fitMethod=None, doFit=True, doOtherFit=True):
  """Fit a newly created Nmr.Peak to a position in the spectral data specified by its dataSource

  .. describe:: Input
  
  Nmr.Peak (new)
  Optionally String describing the fitMethod

  .. describe:: Output
  
  None
  """
  
  if (not doFit) and (not doOtherFit):
    return

  spectrum = peak.peakList.dataSource
  if not hasattr(spectrum, 'block_file'):
    return
  block_file = spectrum.block_file
  if not block_file:
    return

  if not fitMethod:
    fitMethod = PEAK_FIT_METHODS[0]

  if fitMethod not in PEAK_FIT_METHODS:
    return

  dimDone = getDimChecked(spectrum)
  if fitMethod == PEAK_FIT_METHODS[0]:
    if doFit:
      method = 0
      center = peak.cPeak.fitCenter(method, block_file, dimDone)
      if center is not None:
        dim = 0
        # TBD: should it be dim or peakDim.dim-1 below?
        for peakDim in peak.sortedPeakDims():
          peakDim.position = int(math.floor(peakDim.position+0.5)) + center[dim]
          dim = dim + 1
          
        # not very elegant but ought to be what we want
        peak.fitMethod = getMethod(peak.root, task='fit peak',
                                   procedure='%s peak fit' % fitMethod,
                                   parameters=(('method', fitMethod),))
    if doOtherFit:
      setupPeakHeight(peak)
      setupPeakLinewidth(peak)
  else:
    fitPeaks([peak], fitMethod, updatePosition=doFit)

    # not very elegant but ought to be what we want
    peak.fitMethod = getMethod(peak.root, task='fit peak',
                               procedure='%s peak fit' % fitMethod,
                               parameters=(('method', fitMethod),))

def fitPeaks(peaks, fitMethod, updatePosition=True):
  """ Fit peaks using the specified fitMethod.
      Input: peaks, all assumed to be from the same peakList.
             fitMethod (String)
      Output: None
  """
  try:
    from ccpnmr.c.PeakList import fitPeaksInRegion
  except Exception, e:
    print 'Could not import fitPeaksInRegion'
    print e
    return

  fittedPeaks = []
  cPeaks = []
  first = None
  last = None
  peakList = None
  for peak in peaks:
    if not hasattr(peak, 'cPeak'):
      continue
    cPeak = peak.cPeak
    if not cPeak:
      continue
    if peakList:
      if peak.peakList is not peakList:
        raise Exception('Mismatched peakLists: %s and %s' % (peakList, peak.peakList))
    else:
      peakList = peak.peakList
      spectrum = peakList.dataSource
      if not hasattr(spectrum, 'block_file') or not spectrum.block_file:
        return
      block_file = spectrum.block_file

    fittedPeaks.append(peak)
    cPeaks.append(cPeak)

    if not first:
      ndim = len(peak.peakDims)
      first = ndim * [None]
      last = ndim * [None]
      dimDone = ndim * [None]

    for peakDim in peak.sortedPeakDims():
      dim = peakDim.dim - 1
      dataDim = peakDim.dataDim
      if dataDim.className == 'FreqDataDim':
        boxWidth = peakDim.boxWidth
        if not boxWidth: # if None (or even if 0)
          boxWidth = getPeakFindBoxwidth(dataDim)
        halfBoxWidth = max(1, int(boxWidth/2))
        position = int(math.floor(peakDim.position)) - 1
        f = max(position - halfBoxWidth, 0)
        l = min(position + halfBoxWidth + 1, dataDim.numPoints)
        if first[dim] is None:
          first[dim] = f
          last[dim] = l
          dimDone[dim] = 1

        else:
          first[dim] = min(first[dim], f)
          last[dim] = max(last[dim], l)
      else:
        position = int(math.floor(peakDim.position+0.5)) - 1
        first[dim] = position
        last[dim] = first[dim] + 1
        dimDone[dim] = 0
          
  if not cPeaks:
    return

  apiFitMethod = getMethod(fittedPeaks[0].root, task='fit peak',
                           procedure='%s peak fit' % fitMethod,
                           parameters=(('method', fitMethod),))
  method = PEAK_FIT_METHODS.index(fitMethod) - 1
  result = fitPeaksInRegion(method, cPeaks, block_file, first, last, dimDone)
  #for i, (height, position, lineWidth, heightDev, positionDev, lineWidthDev) in enumerate(result):
  for i, (height, position, lineWidth) in enumerate(result):
    peak = fittedPeaks[i]
    spectrum = peak.peakList.dataSource
    height *= spectrum.scale
    peak.fitMethod = apiFitMethod
    setManualPeakIntensity(peak, height)
    for i, peakDim in enumerate(peak.sortedPeakDims()):
      if dimDone[i]:
        if updatePosition:
          peakDim.position = position[i]
        dataDimRef = peakDim.dataDimRef 
        peakDim.lineWidth = abs(convertPosition(lineWidth[i], dataDimRef, toUnit='Hz', relative=True))
 
def setupPeakHeight(peak):
  """Determine a peaks height from the spectral data specified by its dataSource

  .. describe:: Input
  
  Nmr.Peak
  .. describe:: Output
  
  None
  """

  spectrum = peak.peakList.dataSource
  if (not hasattr(spectrum, 'block_file')) or (not spectrum.block_file):
    return

  prevHt = peak.findFirstPeakIntensity(intensityType='height')
  if prevHt:
    prevHt.delete()
    
  height = peak.cPeak.getIntensity(spectrum.block_file)
  height *= spectrum.scale
  method = getMethod(peak.root, task='find peak height',
                     procedure='peak.cPeak.getIntensity')
  peak.newPeakIntensity(intensityType='height', value=height, method=method)
 
 
def setManualPeakIntensity(peak, value, intensityType='height'):
  """Manually define a peaks volume/height

  .. describe:: Input
  
  Nmr.Peak, Float (Nmr.PeakIntensity.value), String (Nmr.PeakIntensity.intensityType)
  .. describe:: Output
  
  None
  """
 
  if (value is not None):
    method = getMethod(peak.root, task='set peak intensity', procedure='manual %s' % intensityType)
    peakIntensity = peak.findFirstPeakIntensity(intensityType=intensityType)
    if peakIntensity and (peakIntensity.method is not method):
      peakIntensity.delete()
    
    peakIntensity = peak.findFirstPeakIntensity(intensityType=intensityType, method=method)
    if peakIntensity:
      if peakIntensity.value != value:
        peakIntensity.setValue(value)
    else:
      peak.newPeakIntensity(intensityType=intensityType, value=value, method=method)
 
 
def getPeakVolume(peak):
  """Gives a peak's volume

  .. describe:: Input
  
  Nmr.Peak
  .. describe:: Output
  
  Float (Nmr.PeakIntensity.value)
  """

  if not peak:
    return
    
  peakIntensity = peak.findFirstPeakIntensity(intensityType='volume')

  if peakIntensity:
    return peakIntensity.value


def getPeakHeight(peak):
  """Gives a peak's height

  .. describe:: Input
  
  Nmr.Peak
  .. describe:: Output
  
  Float (Nmr.PeakIntensity.value)
  """

  if not peak:
    return
    
  peakIntensity = peak.findFirstPeakIntensity(intensityType='height')

  if peakIntensity:
    return peakIntensity.value


def getPeakIntensity(peak, intensityType):
  """Gives a peak's intensity value according to the input intensity type

  .. describe:: Input
  
  Nmr.Peak
  .. describe:: Output
  
  Float (Nmr.PeakIntensity.value)
  """
 
  if not peak:
    return

  peakIntensity = peak.findFirstPeakIntensity(intensityType=intensityType)

  if peakIntensity:
    return peakIntensity.value
  
  
def findPeakBoxValues(peak):
  """Estimates the value of a peak's box width

  .. describe:: Input
  
  Nmr.Peak
  .. describe:: Output
  
  List of List of Floats (min,max points positions for dim), List of Floats (Nmr.PeakDim.BoxWidth in points), List of Floats (centre position)
  """

  spectrum = peak.peakList.dataSource
  block_file = spectrum.block_file
  if (not block_file):
    return None

  # TBD: look at again when boxWidth has units in data model
  peakDims = peak.sortedPeakDims()
  boxMin = []
  boxMax = []
  boxSize = []
  center = []
  for peakDim in peakDims:
    dataDim = peakDim.dataDim
    boxWidth = peakDim.boxWidth
    if not boxWidth: # if None (or even if 0)
      boxWidth = getPeakFindBoxwidth(dataDim)
    halfBoxWidth = max(1, int(boxWidth/2))
    position = int(0.5 + peakDim.position - 1) # round to nearest integer
    if dataDim.className == 'FreqDataDim':
      a = max(0, position-halfBoxWidth)
      b = min(dataDim.numPoints, position+halfBoxWidth+1)
    else:
      a = position
      b = a + 1
    boxMin.append(a)
    boxMax.append(b)
    boxSize.append(b-a)
    center.append(position-a)

  #print 'findPeakBoxValues1', boxMin, boxMax
  values = block_file.getValues(boxMin, boxMax)
  #print 'findPeakBoxValues2', len(values)

  return (values, boxSize, center)


def findPeakVolume(peak, volumeMethod=None):
  """Finds the volume of a peak using the volume method passed in, or if not set then the one stored in application specific data.

  .. describe:: Input
  
  Nmr.Peak
  .. describe:: Output
  
  Float
  """

  spectrum = peak.peakList.dataSource
  block_file = spectrum.block_file
  if not block_file:
    return None

  if not volumeMethod:
    volumeMethod = getPeakFindParams(spectrum.root)['volumeMethod']
  if volumeMethod in ('parabolic fit', 'peak fit'):
    fitMethod = peak.fitMethod
    if not fitMethod:
      return

    parameter = fitMethod.findFirstParameter(name='method')
    if not parameter:
      return

    method = parameter.value
    if method == 'parabolic':
      method = 0 # only one so far for C world
      dimDone = getDimChecked(spectrum)
      volume = peak.cPeak.fitVolume(method, block_file, dimDone)
      #print 'findPeakVolume1', volume

    else:
      heightIntensity = peak.findFirstPeakIntensity(intensityType='height')
      if not heightIntensity:
        return None
      volume = heightIntensity.value
      if method == 'gaussian':
        mult = math.sqrt(math.pi/math.log(16))
      else:
        mult = math.pi / 2
      for peakDim in peak.peakDims:
        if peakDim.dataDim.className != 'FreqDataDim':
          continue
        lineWidth = peakDim.lineWidth
        if not lineWidth:
          return None
        volume *= mult * peakDim.lineWidth
        
  else:
    (values, boxSize, center) = findPeakBoxValues(peak)
    if volumeMethod == 'truncated box sum':
      volume = truncatedBoxIntegral(values, boxSize, center)
      #print 'findPeakVolume2', volume
    elif volumeMethod == 'box sum':
      #volume = reduce(operator.add, values)
      volume = sum(values)
      #print 'findPeakVolume3', volume
    else:
      raise ApiError('unknown volumeMethod "%s" in findPeakVolume()' % volumeMethod)

  volume *= spectrum.scale

  return volume


def setupPeakVolume(peak, volumeMethod=None):
  """Determine a peaks volume from the spectral data specified by its dataSource

  .. describe:: Input
  
  Nmr.Peak
  .. describe:: Output
  
  None
  """

  spectrum = peak.peakList.dataSource
  if (not hasattr(spectrum, 'block_file')) or (not spectrum.block_file):
    return

  if not volumeMethod:
    volumeMethod = getPeakFindParams(spectrum.root)['volumeMethod']
  volume = findPeakVolume(peak, volumeMethod)

  if volume is not None:
    prevVol = peak.findFirstPeakIntensity(intensityType='volume')
    if prevVol:
      prevVol.delete()
    method = getMethod(peak.root, task='fit peak volume',
                       procedure='peak.cPeak.fitVolume',
                       parameters=(('method', volumeMethod),))
    peak.newPeakIntensity(intensityType='volume', value=volume, method=method)
 
 
def setupPeakLinewidth(peak):
  """Fit linewidth of Nmr.Peak

  .. describe:: Input
  
  Nmr.Peak
  .. describe:: Output
  
  None
  """
  
  spectrum = peak.peakList.dataSource
  if (not hasattr(spectrum, 'block_file')) or (not spectrum.block_file):
    return

  dimDone = getDimChecked(spectrum)
  linewidth = peak.cPeak.fitLinewidth(spectrum.block_file, dimDone)
  # TBD: should it be dim or peakDim.dim-1 below?

  for peakDim in peak.peakDims:
    dataDim = peakDim.dataDim
    dim = dataDim.dim - 1
    dataDimRef = peakDim.dataDimRef 
    if dimDone[dim] and dataDimRef:
      # TBD: should this always be in Hz??
      peakDim.lineWidth = abs(convertPosition(linewidth[dim], dataDimRef, toUnit='Hz', relative=True))



def addPeakToSelected(peak, parent=None):
  """Add a peak to the GUI current peak selection

  .. describe:: Input
  
  Nmr.Peak, Tkinter.widget 
  .. describe:: Output
  
  None
  """

  if parent:
    parent.addSelected(peak)


def removePeakFromSelected(peak, parent=None):
  """Remove a peak to the GUI current peak selection

  .. describe:: Input
  
  Nmr.Peak, Tkinter.widget 
  .. describe:: Output
  
  None
  """

  if parent:
    parent.removeSelected(peak)


# below assumes region is in points, with region[0] < region[1]
# below assumes that peakDim.dataDimRef is defined
def isPeakInRegion(peak, region, acceptAlias=True):
  """Determine whether a peak is within a given region optionally considering aliasing

  .. describe:: Input
  
  Nmr.Peak, List of List of Floats (min,max poitions), Boolean
  .. describe:: Output
  
  Boolean
  """

  n = 0
  for peakDim in peak.sortedPeakDims():
    p = peakDim.position - 1
    (r0, r1) = region[n]
    #print 'isPeakInRegion', n, peak.serial, peakDim.dim, p, r0, r1
    npoints = peakDim.dataDim.numPointsOrig
    if (acceptAlias):
      m = int(math.floor(float(p - r0) / npoints))
      p = p - m * npoints
      #assert p >= r0 and p <= r0+npoints, 'p = %f, r0 = %f, npoints = %d' % (p, r0, npoints)
      if (p > r1):
        return False
    else:
      p = p + peakDim.numAliasing*npoints
      if ((p < r0) or (p > r1)):
        return False
    n = n + 1
 
  #print 'isPeakInRegion returning True', peak.serial

  return True


# internal function
# input is region defined in units specified by axisUnit
def findDataDimRegions(dataDim, region, axisUnit, thickness, addOneToUpperPoint=True):
  """Define a region of a data dimension in terms of a  given axis unit with extra thickness if required

  .. describe:: Input
  
  Nmr.FreqDataDim, List of List of Floats (min,max poitions), Word (ShiftList.unit), List of Floats, Bool
  .. describe:: Output
  
  List of Tuple of Floats min,max poitions), Float (min value), Float (max value)
  """

  # find contiguous regions: there can be 0, 1 or 2
 
  (p0, p1) = convertRegion(region, axisUnit, dataDim)
  #print 'findDataDimRegions0', dataDim.dim, region, p0, p1
  # 15 Sep 2009: subtract and add 0.5 because the grid point maximum
  # can be outside the search bin but the interpolated maximum inside
  # the bin, and without the -+ 0.5 those maxima will be missed
  p0 = int(math.ceil(p0-0.5) - thickness)
  p1 = int(math.floor(p1+0.5) + thickness)
  if addOneToUpperPoint:
    p1 += 1
  #print 'findDataDimRegions0A', dataDim.dim, region, p0, p1
 
  if dataDim.className == 'FreqDataDim':
    # TBD: assumes that axisType.axisUnits[0] is same as min(max)AliasedFreq units
    dataDimRef =  getPrimaryDataDimRef(dataDim)
    expDimRef = dataDimRef.expDimRef
    minAliasedFreq = expDimRef.minAliasedFreq
    maxAliasedFreq = expDimRef.maxAliasedFreq
    # note p0, p1 in points so maxAlised for min and minAliased for max
    if (maxAliasedFreq is None): # restrict to fundamental region
      p0 = max(p0, 0)
    if (minAliasedFreq is None): # restrict to fundamental region
      p1 = min(p1, dataDim.numPoints-1) # TBD: is -1 correct?
    p1 = max(p0, p1) # safety
 
    min_point = p0
    max_point = p1
 
    n = dataDim.numPointsOrig
    # below is wrong, a should be counted for data set as is,
    # not data set as it would have been if it hadn't had an offset
    ###a = dataDim.pointOffset
    a = 0
    b = a + dataDim.numPoints
    #print 'findDataDimRegions1', dataDim.dim, p0, p1, a, b, n
 
    d = p1 - p0
    p0 = p0 % n
    p1 = p0 + d
    #print 'findDataDimRegions2', dataDim.dim, p0, p1
 
    if (p0 <= a):
      if (p1 <= a):
        regions = []
      else:
        regions = [ (0, min(b,p1)-a) ]
    elif (p0 >= b):
      if (p1 <= n+a):
        regions = []
      else:
        regions = [ (0, min(b,p1-n)-a) ]
    else:
      if (p1 <= n+a):
        regions = [ (p0-a, min(b,p1)-a) ]
      elif (p1 >= n+p0):
        regions = [ (0, b-a) ]
      else:
        regions = [ (0, p1-n-a), (p0-a, b-a) ]
 
  else:

    #n = dataDim.numPoints - 1
    n = dataDim.numPoints
    p0 = min(max(0, p0), n)
    p1 = min(max(0, p1), n)

    min_point = p0
    max_point = p1
    regions = [ (p0, p1) ]
 
  #print 'findDataDimRegions3', dataDim.dim, regions
 
  return (regions, min_point, max_point)


# TBD: region assumed to be in ppm
def findPeaks(peakList, region, parent=None, thickness=None,
              excludedRegions=None, intensityThreshold=None):
  """
  Pick new peaks in a given region (units: ppm) of a peak list.
  Thickness (units: points) can be specified for orthogonal dimensions
  to include extra planes (must be of length ndim).
  ExcludedRegions (units: points) can be specified to exclude a list of
  regions.  Each list is an ndim list of 2-tuples, with each 2-tuple the
  (min, max) of that dimension.
  Peaks are added to the GUI current selection.
  
  .. describe:: Input
  
  Nmr.PeakList, List of List of Floats (min,max positions) in ppm, 'w'idget,
  List of Floats, List of List of 2-tuples of Floats in points
  
  .. describe:: Output
  
  List of Nmr.Peaks
  """

  #print 'findPeaks0', peakList.serial, region, thickness

  spectrum = peakList.dataSource
  block_file = spectrum.block_file

  if not block_file:
    return []

  ndim = spectrum.numDim
  project = spectrum.root
  analysisSpectrum = getAnalysisSpectrum(spectrum)
  analysisProject = analysisSpectrum.analysisProject
  axisUnit = analysisProject.findFirstAxisUnit(unit='ppm')

  params = getPeakFindParams(project)
  scale = params['scale']
  buf = params['buffer']
  drop_factor = params['drop']
  nonadjacent = params['nonadjacent']
  have_high = params['haveHigh']
  have_low = params['haveLow']

  if (thickness is None):
    thickness = ndim * [0]

  if (excludedRegions is None):
    excludedRegions = []

  dimChecked = getDimChecked(spectrum)

  scale = spectrum.scale / (scale * analysisProject.globalContourScale)

  high = 0
  if (have_high):
    if intensityThreshold is None:
      posLevels = analysisSpectrum.posLevels
      if posLevels:
        high = min(posLevels) / scale
      else:
        have_high = False
    else:
      high = abs(intensityThreshold)

  low = 0
  if have_low:
    if intensityThreshold is None:
      negLevels = analysisSpectrum.negLevels
      if negLevels:
        low = max(negLevels) / scale
      else:
        have_low = False
    else:
      low = -abs(intensityThreshold)
    
  buff = ndim * [0]
  for i in range(ndim):
    dataDim = spectrum.findFirstDataDim(dim=i+1)
    if dataDim.className == 'FreqDataDim':
      buff[i] = buf
    else:
      buff[i] = 0

  min_linewidth = ndim * [0]
  for dataDim in spectrum.dataDims:
    if dataDim.className == 'FreqDataDim':
      min_lw = getPeakFindMinLinewidth(dataDim)
      min_lw = convertPosition(min_lw,  getPrimaryDataDimRef(dataDim),
                               fromUnit='Hz', relative=True)
    else:
      min_lw = 0
    min_linewidth[dataDim.dim-1] = abs(min_lw)

  regions = ndim * [0]
  min_point = ndim * [0]
  max_point = ndim * [0]

  dataDims = spectrum.sortedDataDims()
  rdim = 0
  missingDim = (len(region) < len(dataDims))
  for dataDim in dataDims:
    dim = dataDim.dim - 1
    if dataDim.className == 'FreqDataDim':
      (r0, r1) = ppmDataDimBoundedRegion(region[rdim], dataDim)
      rdim += 1

      (regions[dim], min_point[dim], max_point[dim]) = \
           findDataDimRegions(dataDim, (r0, r1), axisUnit, thickness[dim])
    else:
      npoints = dataDim.numPoints
      if missingDim:
        r0, r1 = (0, npoints) # whole sampled region
      else:
        (r0, r1) = region[rdim]
        r0 = int(r0-1)
        r1 = int(r1)
        rdim += 1
      regions[dim] = [(r0, r1)]
      min_point[dim] = min(max(r0, 0), npoints)
      max_point[dim] = min(max(r1, 0), npoints)

  #print 'findPeaks0A', regions

  peaks = []

  tile_min = ndim * [0]
  tile_max = ndim * [0]
  ntiles_array = ndim * [0]
  tile_offset = ndim * [0]
  tile_abs = ndim * [0]

  for dataDim in spectrum.dataDims:
    dim = dataDim.dim - 1
    if dataDim.className == 'FreqDataDim':
      n = dataDim.numPointsOrig
    else:
      n = dataDim.numPoints
      
    t0 = int(math.floor(min_point[dim] / n))
    t1 = int(math.floor(max_point[dim] / n))
    tile_min[dim] = t0
    tile_max[dim] = t1

  diagonalExclusions = getDiagonalExclusions(spectrum)

  first = ndim * [0]
  last = ndim * [0]
  (nregions, cumulative) = cumulativeProductArray([len(reg) for reg in regions])
  for r in range(nregions):
    a = arrayOfIndex(r, cumulative)
    reg = [ regions[i][a[i]] for i in range(len(a)) ]
    for dataDim in spectrum.dataDims:
      dim = dataDim.dim - 1
      #first[dim] = int(reg[dim][0])
      #last[dim] = int(reg[dim][1])
      first[dim] = reg[dim][0]
      last[dim] = reg[dim][1]

    #print 'findPeaks1', first, last, have_high, have_low, high, low, buffer, nonadjacent, drop_factor, min_linewidth
    positions = peakList.cPeakList.findPeaks(first, last,
                  block_file, have_high=have_high, have_low=have_low,
                  high=high, low=low, buffer=buff, nonadjacent=nonadjacent,
                  drop_factor=drop_factor, min_linewidth=min_linewidth,
                  diagonal_exclusions=diagonalExclusions, excluded_regions=excludedRegions,
                  dim_checked=dimChecked)
    #print 'findPeaks2', positions

    for position in positions:
      for dataDim in spectrum.dataDims:
        dim = dataDim.dim - 1
        if dataDim.className == 'FreqDataDim':
          n = dataDim.numPointsOrig
        else:
          n = dataDim.numPoints

        t = tile_min[dim]
        min_t = None
        while (t <= tile_max[dim]):
          p = position[dim] + t*n  # actual position with aliasing accounted for
          if (p >= min_point[dim]):
            min_t = t
            break
          t = t + 1
        if (min_t is None):
          break

        t = tile_max[dim]
        max_t = None
        while (t >= min_t):
          p = position[dim] + t*n  # actual position with aliasing accounted for
          if (p <= max_point[dim]):
            max_t = t
            break
          t = t - 1
        if (max_t is None):
          break

        ntiles_array[dim] = max_t - min_t + 1
        tile_offset[dim] = min_t
        #print 'findPeaks3', dim, min_point[dim], max_point[dim], position[dim], n, min_t, max_t

      (ntiles, cum_tiles) = cumulativeProductArray(ntiles_array)
      for tile_ind in range(ntiles):
        tile_rel = arrayOfIndex(tile_ind, cum_tiles)
        for dim in range(ndim):
          tile_abs[dim] = tile_rel[dim] + tile_offset[dim]
        #print 'findPeaks4', tile_ind, tile_rel, tile_abs, position
        peak = addPeak(peakList, position, tile=tile_abs, parent=parent)
        peaks.append(peak)
      
  return peaks

def getDiagonalExclusions(spectrum):
  """Get the diagonal exclusions for a given spectrum (used in C world).

  .. describe:: Input
  
  Nmr.DataSource
  .. describe:: Output
  
  List of diagonal exclusions, each is a 5-tuple:
             (dim1, dim2, dim1 scale, dim2 scale, dim1 to dim2 translation)
  """

  project = spectrum.root
  isotopes = {}
  for dataDim in spectrum.dataDims:
    if isinstance(dataDim, SampledDataDim):
      continue
  
    dataDimRef = getPrimaryDataDimRef(dataDim) # Filter for main TBD
    
    if dataDimRef:
      isotope = dataDimRef.expDimRef.isotopes[0]
      isotopes[dataDim] = isotope

  diagonalExclusions = []
  for dataDim1 in isotopes.keys():
    isotope1 = isotopes[dataDim1]
    exclusion = getIsotopeExclusion(isotope1)
    if exclusion == 0:
      continue
    for dataDim2 in isotopes.keys():
      if dataDim1.dim < dataDim2.dim:
        isotope2 = isotopes[dataDim2]
        if isotope1 == isotope2:
          diagonalExclusion = getDiagonalExclusion(dataDim1, dataDim2, exclusion)
          diagonalExclusions.append(diagonalExclusion)

  return diagonalExclusions


def getDiagonalExclusion(dataDim1, dataDim2, exclusion):
  """Get the diagonal exclusion for (dataDim1, dataDim2)

  .. describe:: Input
  
  Nmr.AbstractDataDim, Nmr.AbstractDataDim
  .. describe:: Output
  
  The diagonal exclusion as a 5-tuple:
             (dim1, dim2, dim1 scale, dim2 scale, dim1 to dim2 translation)
  """

  (a1, b1) = getExclusionInfo(dataDim1)
  (a2, b2) = getExclusionInfo(dataDim2)

  return (dataDim1.dim-1, dataDim2.dim-1, a1, a2, b1-b2, exclusion)


def getExclusionInfo(dataDim):
  """Get the exclusion information for dataDim

  .. describe:: Input
  
  Nmr.AbstractDataDim
  .. describe:: Output
  
  The exclusion information as a 2-tuple:
             (scale, translation)
  """

  npoints = dataDim.numPoints
  dataDimRef = getPrimaryDataDimRef(dataDim)
  sw = dataDim.spectralWidth
  sf = dataDimRef.expDimRef.sf
  refpt = dataDimRef.refPoint
  refppm = dataDimRef.refValue

  a = - sw / (npoints * sf)
  b = - refpt * a + refppm

  return (a, b)


# TBD: region assumed to be in ppm
def searchPeaks(peakLists, region, parent=None, thickness=None, considerAliased=True):
  """
  Find existing peaks within a given region of a peak list. Add peaks to the
  current GUI selection. Thickness (units: points) can be specified for
  orthogonal dimensions to include extra planes (must be of length ndim or is
  ignored for a given peakList).
  
  .. describe:: Input
  
  Nmr.PeakList, List of List of Floats (min,max positions), 'w'idget, List of Floats
  
  .. describe:: Output
  
  List of Nmr.Peaks
  """

  #print 'searchPeaks1', region, parent
  all_peaks = []
  for peakList in peakLists:
    if not (hasattr(peakList, 'cPeakList') and peakList.cPeakList):
      continue
    spectrum = peakList.dataSource
    project  = spectrum.root
    axisUnit = project.currentAnalysisProject.findFirstAxisUnit(unit='ppm')
    ndim = spectrum.numDim
    
    if (thickness is None or len(thickness) != ndim):
      thickness = ndim * [0]

    first = ndim * [0]
    last  = ndim * [0]
    dataDims = spectrum.sortedDataDims()
    rdim = 0
    missingDim = (len(region) < len(dataDims))
    for dataDim in dataDims:
      dim = dataDim.dim - 1 
      if isinstance(dataDim, SampledDataDim):
        if missingDim:
          r0, r1 = (1, dataDim.numPoints) # whole sampled region
        else:
          (r0, r1) = region[rdim]
          rdim += 1
        r0 = r0 - 1
        r1 = r1 - 1
      else:
        r = ppmDataDimBoundedRegion(region[rdim], dataDim)
        (r0, r1) = convertRegion(r, axisUnit, dataDim)
        rdim += 1
        
      r0 = r0 - thickness[dim]
      r1 = r1 + thickness[dim]
      r = (r0, r1)
      (r0, r1) = ptsDataDimBoundedRegion(r, dataDim)
      if r0 == r1:
        break
      first[dim] = r0
      last[dim] = r1
    else:
      #print 'searchPeaks2', first, last
      allow_aliasing = getDimWrapped(spectrum)
      for d in range(ndim):
        allow_aliasing[d] &= considerAliased
      peakIndices = peakList.cPeakList.searchPeaks(first, last, allow_aliasing)
      #print 'searchPeaks3', peakIndices
      peaks = peakList.sortedPeaks()
      for ind in peakIndices:
        peak = peaks[ind]
        addPeakToSelected(peak, parent=parent)
        all_peaks.append(peak)

  return all_peaks


def snapPeaks(peaks, doFit=True):
  """Snap peaks to the local maximum
             (optionally do a fit for the center point)

  .. describe:: Input
  
  List of Nmr.Peaks, Boolean
  .. describe:: Output
  
  None
  """

  for peak in peaks:
    snapPeak(peak, doFit)


def snapPeak(peak, doFit=True):
  """Snap peak to the local maximum
             (optionally do a fit for the center point)

  .. describe:: Input
  
  List of Nmr.Peaks, Boolean
  .. describe:: Output
  
  None
  """

  spectrum = peak.peakList.dataSource
  block_file = spectrum.block_file
  if not block_file:
    return

  project = spectrum.root
  ndim = spectrum.numDim

  position = [ peakDim.position-1 for peakDim in peak.sortedPeakDims() ]
  value = block_file.getValue(position)
  if value > 0:
    have_high = True
    have_low = False
    high = value
    low = 0
  else:
    have_high = False
    have_low = True
    high = 0
    low = value

  params = getPeakFindParams(project)
  #scale = params['scale']
  buf = params['buffer']
  drop_factor = params['drop']
  nonadjacent = params['nonadjacent']

  first = ndim * [0]
  last = ndim * [0]
  buff = ndim * [0]

  for i in range(ndim):
    dataDim = spectrum.findFirstDataDim(dim=i+1)
    peakDim = peak.findFirstPeakDim(dim=dataDim.dim)
    if dataDim.className == 'FreqDataDim':
      boxWidth = peakDim.boxWidth
      if not boxWidth: # if None (or even if 0)
        boxWidth = getPeakFindBoxwidth(dataDim)
      halfBoxWidth = dataDim.numPoints / 100  # a bit of a hack
      halfBoxWidth = max(halfBoxWidth, 1, int(boxWidth/2))
      first[i] = max(0, int(math.floor(position[i]-halfBoxWidth)))
      last[i] = min(dataDim.numPoints, int(math.ceil(position[i]+1+halfBoxWidth)))
      buff[i] = buf
    else:
      first[i] = max(0, int(math.floor(position[i]+0.5)))
      last[i] = min(dataDim.numPoints, first[i]+1)
      buff[i] = 0

  min_linewidth = ndim * [0]
  for dataDim in spectrum.dataDims:
    if dataDim.className == 'FreqDataDim':
      min_lw = getPeakFindMinLinewidth(dataDim)
      min_lw = convertPosition(min_lw, getPrimaryDataDimRef(dataDim), fromUnit='Hz',
                               relative=True)
    else:
      min_lw = 0
    min_linewidth[dataDim.dim-1] = abs(min_lw)

  diagonalExclusions = getDiagonalExclusions(spectrum)
  excludedRegions = []
  dimChecked = getDimChecked(spectrum)

  peakNumber = list(peak.peakList.sortedPeaks()).index(peak)
  positions = peak.peakList.cPeakList.findPeaks(first, last,
                  block_file, have_high=have_high, have_low=have_low,
                  high=high, low=low, buffer=buff, nonadjacent=nonadjacent,
                  drop_factor=drop_factor, min_linewidth=min_linewidth,
                  diagonal_exclusions=diagonalExclusions, excluded_regions=excludedRegions,
                  dim_checked=dimChecked, ignore_peak=peakNumber)

  if len(positions) == 1:
    position = positions[0]
    tile = [ peakDim.numAliasing for peakDim in peak.sortedPeakDims() ]
    setPeakPosition(peak, position, tile)
    setupPeak(peak, doFit=doFit)
    checkContourRegion(peak)
    
    
def makePeakListFromShifts(spectrum, useUnassigned=True, progressBar=None,
                           shiftList=None, molSystem=None, bondLimit=6,
                           residueLimit=1, labelling=None, labellingThreshold=0.1):
  """
  Make an artificial peak list using shift intersections from a shift list
  Boolean option to consider only shifts with atom assigned resonances. Option
  to filter on a specific molSystem. Option to limit through-space transfers to a
  limited number of bonds and a given residue range.
  Optional labelling scheme/mixture and threshold to filter according to isotopomers.
 
  .. describe:: Input
  
  Nmr.DataSource, Boolean, memops.gui.ProgressBar, Nmr.ShiftList
  MolSystem.MolSystem, Int, 2-Tuple of Ints,
  ChemComLabel.LabelingScheme or True (automatic from experiment MolLabel),
  Float

  .. describe:: Output
  
  None
  """
  from ccpnmr.analysis.core.MoleculeBasic import areResonancesBound

  peakList = None
  experiment = spectrum.experiment
  refExperiment = experiment.refExperiment
  refType = refExperiment.name
  
  if not shiftList:
    shiftList = experiment.shiftList
  
  observableAtoms = {}

  if not (shiftList and spectrum):
    return  

  if labelling is True:
    labelling = spectrum.experiment
    
    if labelling.labeledMixtures:
      getLabelAtomFractions = getExperimentAtomFractions
      getLabelAtomPairFractions = getExperimentAtomPairFractions
    else:
      labelling = None  
  
  elif labelling and labelling.className == 'LabelingScheme':
    getLabelAtomFractions = getSchemeAtomFractions
    getLabelAtomPairFractions = getSchemeAtomPairFractions
        
  N = spectrum.numDim
  isotopes = []
  fullRegion = []
  transferDims = {}

  atomSiteDims = {}
  dimAtomSites = []
  
  for dim, dataDim in enumerate(spectrum.sortedDataDims()):
    atomSites = set()
    expDim = dataDim.expDim
    for expDimRef in expDim.expDimRefs:
      if not expDimRef.refExpDimRef:
        continue
    
      measurement = expDimRef.refExpDimRef.expMeasurement
      for atomSite in measurement.atomSites:
        if atomSite not in atomSiteDims:
          atomSiteDims[atomSite] = []
        
        atomSiteDims[atomSite].append(dim)
        atomSites.add(atomSite)
        
    dimAtomSites.append(atomSites)
    
  for dim, dataDim in enumerate(spectrum.sortedDataDims()):
    dataDimRef = getPrimaryDataDimRef(dataDim)
    expDimRef  = dataDimRef.expDimRef
    fullRange  = getDataDimRefFullRange(dataDimRef)
    fullRegion.append(fullRange)
    isotopes.append( list(expDimRef.isotopeCodes) )
    expTransfers = expDimRef.expTransfers
    
    if expTransfers:         
      for expTransfer in expTransfers:
        transferType = expTransfer.transferType
        expDimRefs2 = list(expTransfer.expDimRefs)
        expDimRefs2.remove(expDimRef)
        dim2 = expDimRefs2[0].expDim.dim-1
        transferDims[frozenset([dim,dim2])]= transferType
     
    else:
      refExpDimRef = expDimRef.refExpDimRef
      
      expDimRef2 = None
      atomSites2 = set()
      for atomSite in refExpDimRef.expMeasurement.atomSites:
        for expTransfer in atomSite.expTransfers:
          if expTransfer.transferType not in ('onebond','Jcoupling'):
            continue
        
          atomSites2.update(expTransfer.atomSites)
      
      for atomSite in atomSites2:
        if atomSite in atomSiteDims:
          continue
          
        for expTransfer in atomSite.expTransfers:
          if expTransfer.transferType not in ('onebond','Jcoupling'):
            continue
 
          for atomSite2 in expTransfer.atomSites:
            if atomSite2 not in atomSiteDims:
              continue
              
            for dim2 in atomSiteDims[atomSite2]:
              if dim2 == dim:
                continue
              
              transferDims[frozenset([dim,dim2])]= 'twostep'
  
  dims = range(N)
  dimResonances = []
  for i in dims:
    dimResonances.append([])

  measurements = shiftList.measurements
  tick = len(measurements)//100
  tick = max(1, tick)
  if progressBar:
    progressBar.total = 100
    progressBar.set(0)
    progressBar.setText('Filtering shifts')
  
  minDimShifts = []
  maxDimShifts = []
  for i in dims:
    minShifts = [aSite.minShift for aSite in dimAtomSites[i]]
    maxShifts = [aSite.maxShift for aSite in dimAtomSites[i]]
    
    if minShifts:
      minDimShifts.append(min(minShifts))
    else:
      minDimShifts.append(None)
    
    if maxShifts:
      maxDimShifts.append(max(maxShifts))
    else:
      maxDimShifts.append(None)
  
  for j, shift in enumerate(measurements):
    if progressBar and (j % tick == 0):
      progressBar.increment()
    
    ppm = shift.value
    resonance = shift.resonance
    isotopeCode = resonance.isotopeCode
            
    atoms = []
    resonanceSet = resonance.resonanceSet
    
    if resonanceSet:
      for atomSet in resonanceSet.atomSets:
        if labelling:
          for atom in atomSet.atoms:
            frac = getLabelAtomFractions(labelling, atom).get(resonance.isotopeCode)
          
            if frac > labellingThreshold:
              atoms.append(atom)
        
        else:  
          atoms.extend(atomSet.atoms)      

    elif not useUnassigned:
      continue
    
    if molSystem:
      if atoms:
        molSystem2 = atoms[0].topObject
        if molSystem2 is not molSystem:
          continue
        
      else:
        continue
    
    if labelling and not atoms:
      continue
    
    for i in dims:
      if isotopeCode in isotopes[i]:
        if (ppm > fullRegion[i][0]) and (ppm < fullRegion[i][1]):
          if (minDimShifts[i] is not None) and (ppm < minDimShifts[i]):
            continue

          if (maxDimShifts[i] is not None) and (ppm > maxDimShifts[i]):
            continue
          
          if useUnassigned or resonanceSet:
            if atoms:
              name = atoms[0].name
          
            dimResonances[i].append((resonance, frozenset(atoms)))

  correlations = {}
  for pair in transferDims:
    dimA, dimB = pair
    transferType = transferDims[pair]
    
    usedA = set()
    usedB = set()
    for resonanceA, atomsA in dimResonances[dimA]:
      keyA = (dimA, resonanceA)
    
      for resonanceB, atomsB in dimResonances[dimB]:
        keyB = (dimB, resonanceB)
        correlate = False
 
        if (transferType == 'onebond') and areResonancesBound(resonanceA,resonanceB):
          if labelling:
            if atomsA and atomsB:
              isotopes = (resonanceA.isotopeCode, resonanceB.isotopeCode)

              atomPairs = []
              for atomA in atomsA:
                for atomB in atomsB:
                  atomPairs.append((atomA, atomB))

              for atomA, atomB in atomPairs:
                fracDict = getLabelAtomPairFractions(labelling, atomA, atomB)
                frac = fracDict.get(isotopes)

                if (frac is not None) and (frac >= labellingThreshold):
                  break

              else:
                # No labelled possibility; skip
                continue

            else:
              # No possible pairs; skip
              continue

          correlate = True
          
        for atomB in atomsB:
          for atomA in atomsA:
            if labelling:
              isotopes = (resonanceA.isotopeCode, resonanceB.isotopeCode)
              fracDict = getLabelAtomPairFractions(labelling, atomA, atomB)
              frac = fracDict.get(isotopes)

              if (frac is None) or (frac < labellingThreshold):
                continue
                               
            if transferType == 'relayed':
              if areAtomsTocsyLinked(atomA, atomB):
                correlate = True
                break
                
            elif transferType in longRangeTransfers:
              
              if atomB is atomA:
                # Not diagonal
                break
            
              if atomA.residue is atomB.residue:
                if bondLimit:
                  numBonds = getNumConnectingBonds(atomA, atomB)
                  if numBonds and (numBonds <= bondLimit):
                    correlate = True
                    break
                
                else:
                  correlate = True
                  break
                
            
              elif residueLimit:
                prevA = getLinkedResidue(atomA.residue, linkCode='prev')
                prevB = getLinkedResidue(atomB.residue, linkCode='prev')
                if (prevA and (atomB.residue is prevA)) or \
                   (prevB and (atomA.residue is prevB)):
                  if bondLimit:
                    numBonds = getNumConnectingBonds(atomA, atomB)
                    if numBonds and (numBonds <= bondLimit):
                      correlate = True
                      break
                                    
                  else:
                    correlate = True
                    break
                  
                if residueLimit == 2:
                  if prevA:
                    prevA2 = getLinkedResidue(prevA, linkCode='prev')
                  else:
                    prevA2 = None
                    
                  if prevB:
                    prevB2 = getLinkedResidue(prevB, linkCode='prev')
                  else:
                    prevB2 = None
                    
                  if (prevA2 and (atomB.residue is prevA2)) or \
                     (prevB2 and (atomA.residue is prevB2)):
                    if bondLimit:
                      numBonds = getNumConnectingBonds(atomA, atomB)
                      if numBonds and (numBonds <= bondLimit):
                        correlate = True
                        break
                    else:
                      correlate = True
                      break

           
            elif transferType == 'twostep':
              if refType in ('H[N[co[CA]]]', 'H[N[co[{CA|ca[C]}]]]'):
                
                if refType == 'H[N[co[CA]]]':
                  atomNames = ('CA',)
                else:
                  atomNames = ('CA','CB',)  
                
                if (atomB.name in atomNames) and (atomA.name == 'N'):
                  prev = getLinkedResidue(atomA.residue, linkCode='prev')
                  if prev and (atomB.residue is prev):
                    correlate = True
                    break
                    
                elif (atomA.name in atomNames) and (atomB.name == 'N'):
                  prev = getLinkedResidue(atomB.residue, linkCode='prev')
                  if prev and (atomA.residue is prev):
                    correlate = True
                    break
             
              elif refType == 'H[N[ca[CO]]]':
                atomNames = ('C',)
  
                if (atomB.name in atomNames) and (atomA.name == 'N'):
                  if atomB.residue is atomA.residue:
                    correlate = True
                    break
                
                  prev = getLinkedResidue(atomA.residue, linkCode='prev')
                  if prev and (atomB.residue is prev):
                    correlate = True
                    break
                    
                elif (atomA.name in atomNames) and (atomB.name == 'N'):
                  if atomB.residue is atomA.residue:
                    correlate = True
                    break
                    
                  prev = getLinkedResidue(atomB.residue, linkCode='prev')
                  if prev and (atomA.residue is prev):
                    correlate = True
                    break

            elif transferType in ('Jcoupling','Jmultibond'):
              if refType in ('H[N[{CA|ca[Cali]}]]',
                             'h{CA|Cca}NH','H[N[CA]]'):
                
                if refType == 'H[N[CA]]':
                  atomNames = ('CA',)
                else:
                  atomNames = ('CA','CB',)  
  
                if (atomB.name in atomNames) and (atomA.name == 'N'):
                  if atomB.residue is atomA.residue:
                    correlate = True
                    break
                
                  prev = getLinkedResidue(atomA.residue, linkCode='prev')
                  if prev and (atomB.residue is prev):
                    correlate = True
                    break
                    
                elif (atomA.name in atomNames) and (atomB.name == 'N'):
                  if atomB.residue is atomA.residue:
                    correlate = True
                    break
                    
                  prev = getLinkedResidue(atomB.residue, linkCode='prev')
                  if prev and (atomA.residue is prev):
                    correlate = True
                    break
  
                    
              elif refType in ('H{[N]+[HA]}', 'H[N[ca[HA]]]'):
                if atomB.name in ('HA','HA2','HA3') and (atomA.name == 'H'):
                  if atomB.residue is atomA.residue:
                    correlate = True
                    break

                if atomA.name in ('HA','HA2','HA3') and (atomB.name == 'H'):
                  if atomB.residue is atomA.residue:
                    correlate = True
                    break
               
              else:
                residue = atomA.residue
                observable = observableAtoms.get(residue)
                if observable is None:
                  observable = getResidueObservableAtoms(residue,
                                                         refExperiment)
                  observableAtoms[residue] = observable
  
                if atomB in observable:
                  correlate = True
                  break
                                     
          else:
            continue
          break
          
        if correlate:
          usedA.add( (resonanceA, atomsA) )
          usedB.add( (resonanceB, atomsB) )
          
          if keyA in correlations:
            correlations[keyA].append(keyB)
          else:
            correlations[keyA] = [keyB]

          if keyB in correlations:
            correlations[keyB].append(keyA)
          else:
            correlations[keyB] = [keyA]
          
    dimResonances[dimA] = usedA
    dimResonances[dimB] = usedB
  
  stack = [ set([key,]) for key in correlations if key[0] == 0]
  intersections = set()
  
  while stack:
    intersection = stack.pop()
    dimNums = set([x[0] for x in intersection])
    
    for keyA in intersection:
      for keyB in correlations[keyA]:
        if keyB in intersection:
          continue
        
        if keyB[0] in dimNums:
          continue
         
        intersection2 = set(intersection)
        intersection2.add(keyB)
        
        if len(intersection2) == N:
          intersections.add(frozenset(intersection2))
          
        else:
          stack.append(intersection2)
   
  if intersections:
    I = len(intersections)
    peakList = spectrum.newPeakList(isSimulated=True)
    peakList.details = 'Synthetic peak list made with shift list "%s"' % shiftList.name

    if progressBar:
      progressBar.total = I
      progressBar.set(0)
      progressBar.setText("Making %d peaks" % I)
      progressBar.open()

    unit = shiftList.unit
    for intersection in intersections:
      resonances = list(intersection)
      resonances.sort()
      resonances = [x[1] for x in resonances]
    
      if progressBar:
        progressBar.increment()
 
      position = []
      figOfMerit = 1.0
      for resonance in resonances:
        shift = resonance.findFirstShift(parentList=shiftList)
        figOfMerit *= shift.figOfMerit
        position.append(shift.value)

      peak = pickPeak(peakList, position, unit=unit, figOfMerit=figOfMerit, doFit=False)
      peakDims = peak.sortedPeakDims()
 
      for dim, resonance in enumerate(resonances):
        assignResToDim(peakDims[dim],resonance,tolerance=10.0)

  
  else:
    msg = 'Experiment type not supported or '
    msg += 'cannot find resonance intersections to make any peaks'
    showWarning('Failure', msg)
    
  if progressBar:
    progressBar.close()
  
  return peakList   
      
def structurePredictNoePeakList(structure, spectrum, distThreshold=5.0, progressBar=None,
                                labelling=None, labellingThreshold=0.1, model=None,
                                minHeight=None, shiftList=None):
  """
  Predict the positions of NOESY peaks given a shift list and a structure
  ensemble. Produces a peak list which peaks that would represent distances less
  than a given threshold. Optional labelling scheme/mixture to filter according
  to isotopomers. Optional structure model to use in ensemble, otherwise all
  models are used -  in this case the average distance must be bwlow threshold.
  Option to specify a minumum hight value so that peaks are not picked below a
  threshold.
  
  .. describe:: Input
  
  MolStructure.StructureEnsemble, Nmr.DataSource, Float, ProgressBar (Tk widget),
  ChemCompLabel.LabelingScheme or rue (automatic from experiment MolLabel),
  Float, MolStructure.Model, Float
  
  .. describe:: Output
  
  Nmr.PeakList
  """

  from ccpnmr.analysis.core.StructureBasic  import getAtomSetsDistance
  from ccpnmr.analysis.core.AssignmentBasic import getAtomSetShifts
  from ccpnmr.analysis.core.ExperimentBasic import getOnebondDataDims
  from ccpnmr.analysis.core.Util import getAnalysisPeakList

  project   = structure.root
  
  if not shiftList:
    shiftList = spectrum.experiment.shiftList
  
  if not shiftList:
    shiftList = project.currentNmrProject.findFirstMeasurementList(className='ShiftList')

  if labelling is True:
    labelling = spectrum.experiment

    if labelling.labeledMixtures:
      getLabelAtomFractions = getExperimentAtomFractions
      getLabelAtomPairFractions = getExperimentAtomPairFractions
    else:
      labelling = None  
 
  elif labelling and labelling.className == 'LabelingScheme':
    getLabelAtomFractions = getSchemeAtomFractions
    getLabelAtomPairFractions = getSchemeAtomPairFractions

  # Works for 4Ds
  bonded   = {}
  dataDims = getOnebondDataDims(spectrum)
  for dataDim1, dataDim2 in dataDims:
    expDim1 = dataDim1.expDim
    expDim2 = dataDim2.expDim
  
    if (not expDim1.expDimRefs) or (not expDim2.expDimRefs):
      continue
    
    isotopes1 = set()
    isotopes2 = set()
    
    for expDimRef in expDim1.expDimRefs:
      if expDimRef.measurementType in ('Shift','shift'):
        if len(expDimRef.isotopeCodes) == 1:
          isotopes1.add(expDimRef.isotopeCodes[0])
    
    for expDimRef in expDim2.expDimRefs:
      if expDimRef.measurementType in ('Shift','shift'):
        if len(expDimRef.isotopeCodes) == 1:
          isotopes2.add(expDimRef.isotopeCodes[0])

    bonded[dataDim1] = (dataDim2, isotopes2)
    bonded[dataDim2] = (dataDim1, isotopes1)
  
  hDataDims = getThroughSpaceDataDims(spectrum)
  peakList = spectrum.newPeakList(isSimulated=True)
  dData = (structure.molSystem.code, structure.ensembleId, distThreshold)
  peakList.details = 'Synthetic through-space peaks for structure %s:%d within %.3f' % dData
  analysisPeakList = getAnalysisPeakList(peakList)
  
  color = '#FF0000'
  if analysisPeakList.symbolColor == color:
    color = '#0000FF'
    
  analysisPeakList.symbolStyle = '+' 
  analysisPeakList.symbolColor = color
  analysisPeakList.textColor = color

  if progressBar:
    i = 0
    for coordChain in structure.coordChains:
      for coordResidue in coordChain.residues:
        i += 1
  
    progressBar.setText('Finding atom coordinates')
    progressBar.total = i or 1
    progressBar.set(0)
    progressBar.open()
  
  from re import match
    
  tsAtoms = set()
  for dataDim in hDataDims:
    dataDimRef = getPrimaryDataDimRef(dataDim)
    for ic in dataDimRef.expDimRef.isotopeCodes:
      aName = match('\d+(\D+)',ic).group(1)
      tsAtoms.add(aName)
      
  atomSetDict = {}
  closeAtomSets = {}
  for coordChain in structure.coordChains:
    for coordResidue in coordChain.residues:
      for coordAtom in coordResidue.atoms:
        if coordAtom.name[0] in tsAtoms:
          atom = coordAtom.atom
          if atom:
            atomSet = atom.atomSet
            
            if atomSet:
              atomSetDict[coordAtom] = atomSet
              closeAtomSets[atomSet] = set([])
            
      if progressBar:
        progressBar.increment()
  
  atomSetPairs = []        
  N = len(closeAtomSets)
  if N == 0:
    showWarning('Failure', 'Cannot find any atom sets for structure')
    return
  
  if model:
    models = [model,]
  else:
    models = structure.models   
  
  distThreshold2 = (distThreshold+1)**2
  for model0 in models:
    sortList = [(c.x, c.y, c.z, atomSetDict[c.atom]) for c in model0.coords if c.atom in atomSetDict]
    sortList.sort()
    
    for i, (x, y, z, atomSet1) in enumerate(sortList[:-1]):
      for x2, y2 ,z2, atomSet2 in sortList[i+1:]:
        dx = x2-x
        if dx > distThreshold+1:
          break
      
        dy = y2-y
        if abs(dy) > distThreshold+1:
          continue
        
        dz = z2-z
        if abs(dz) > distThreshold+1:
          continue
      
        if (dx*dx) + (dy*dy) + (dz*dz) > distThreshold2:
          continue
      
        closeAtomSets[atomSet1].add(atomSet2)
        closeAtomSets[atomSet2].add(atomSet1)
  
  if progressBar:
    progressBar.setText('Finding close atoms amongst %d' % N)
    progressBar.total = N
    progressBar.set(0)
    progressBar.open()
    
  for atomSetI in closeAtomSets.keys():
  
    if progressBar:
      progressBar.increment()
    
    for atomSetJ in closeAtomSets[atomSetI]:
    
      dist = getAtomSetsDistance([atomSetI,],[atomSetJ,],structure,
                                 model=model, method='noe')
      if dist <= distThreshold:
        
        if labelling:
          for atomI in atomSetI.atoms:
	    elemI = atomI.chemAtom.elementSymbol
            isotopeI = DEFAULT_ISOTOPES[elemI]
	    for atomJ in atomSetJ.atoms:
              elemJ = atomJ.chemAtom.elementSymbol
              isotopeJ = DEFAULT_ISOTOPES[elemJ]
              fracDict = getLabelAtomPairFractions(labelling, atomI, atomJ)
              frac = fracDict.get((isotopeI,isotopeJ)) 
              if (frac is not None) and (frac >= labellingThreshold):
                atomSetPairs.append( (atomSetI, atomSetJ, dist) )
                break
            else:
              continue
            break  
           
        else:
          atomSetPairs.append( (atomSetI, atomSetJ, dist) )

  N = len(atomSetPairs)
  if N == 0:
    showWarning('Failure', 'Cannot find any relevent atom pairs in structure')
    return
 
  if hasattr(spectrum, 'block_file') and spectrum.block_file:
    blockFileGetValue = spectrum.block_file.getValue
  else:
    blockFileGetValue= None
 
  
  if progressBar:
    #progressBar.setText('Making %d synthetic peaks' % N)
    progressBar.setText('Making synthetic peaks')
    progressBar.total = N
    progressBar.set(0)
    progressBar.open()

  ranges = {}
  sortedDataDims = spectrum.sortedDataDims()
  for dataDim in sortedDataDims:
    dataDimRef = getPrimaryDataDimRef(dataDim)
    ranges[dataDim] = getDataDimRefFullRange(dataDimRef)

  duplicates = {}
  dataDimRefs = [getPrimaryDataDimRef(dd) for dd in sortedDataDims]
  N = spectrum.numDim
  for atomSet1, atomSet2, dist in atomSetPairs:
    if progressBar:
      progressBar.increment()
    
    shifts1 = getAtomSetShifts(atomSet1, shiftList)
    shifts2 = getAtomSetShifts(atomSet2, shiftList)
    for shift1 in shifts1:
      resonance1 = shift1.resonance
      bound1 = None
      if resonance1.covalentlyBound:
        bound1 = resonance1.findFirstCovalentlyBound().findFirstShift(parentList=shiftList)
     
      
      for shift2 in shifts2:
        if shift1 is shift2:
          continue
      
        resonance2 = shift2.resonance
        bound2 = None
        if resonance2.covalentlyBound:
          bound2 = resonance2.findFirstCovalentlyBound().findFirstShift(parentList=shiftList)
       
        for i in range(2):
          shiftDict = {}
          hDim1 = hDataDims[i]
          hDim2 = hDataDims[1-i]
          shiftDict[hDim1] = shift1
          shiftDict[hDim2] = shift2
 
          xDim1, isotopesX1 = bonded.get(hDim1, (None, None))
          xDim2, isotopesX2 = bonded.get(hDim2, (None, None))
          
          if xDim1 and bound1:
            if bound1.resonance.isotopeCode not in isotopesX1:
              continue
          
            shiftDict[xDim1] = bound1

            if labelling:
              # If we have a labelling scheme, check bound atom is suitably enriched
              resonanceSet1 = bound1.resonance.resonanceSet
              
              if resonanceSet1:
                atom1 = resonanceSet1.findFirstAtomSet().findFirstAtom()
                frac  = getLabelAtomFractions(labelling, atom1).get(bound1.resonance.isotopeCode)
                if frac < labellingThreshold:
                  continue
                  

          if xDim2 and bound2:
            if bound2.resonance.isotopeCode not in isotopesX2:
              continue
              
            shiftDict[xDim2] = bound2

            if labelling:
              # If we have a labelling scheme, check bound atom is suitably enriched
              resonanceSet2 = bound2.resonance.resonanceSet
              
              if resonanceSet2:
                atom2 = resonanceSet2.findFirstAtomSet().findFirstAtom()
                frac  = getLabelAtomFractions(labelling, atom2).get(bound2.resonance.isotopeCode)
                if frac < labellingThreshold:
                  continue
 
          outOfRange = False
          for dataDim in shiftDict.keys():
            minAliasedFreq = ranges[dataDim][0]
            if shiftDict[dataDim].value < minAliasedFreq:
              outOfRange = True
 
            maxAliasedFreq = ranges[dataDim][1]
            if shiftDict[dataDim].value > maxAliasedFreq:
              outOfRange = True
 
          if outOfRange:
            continue
 
          position = []
          for dataDim in spectrum.sortedDataDims():
            shift = shiftDict.get(dataDim)
            if shift:
              position.append(shift.value)
            else:
              break

          else:
            
            value = None
            if blockFileGetValue and (minHeight is not None):
              points = [ppm2pnt(ppm, dataDimRefs[k]) for k, ppm in enumerate(position)]
              
              for j, point in enumerate(points):
                dataDim = sortedDataDims[j]
                p = (point-1) % dataDim.numPointsOrig
                
                if p >= dataDim.numPoints:
                  points[j] = None
                  
                else:
                  points[j] = p  
              
              if None in points:
                continue
              
              value = blockFileGetValue(points)
              if value < minHeight:
                continue
             
          
            key = ':'.join(['%.3f' % x for x in position])
            
            # e.g. redundant prochiral intersections go to the same peak
            peak = duplicates.get(key)
            
            if not peak:
              peak = pickPeak(peakList, position, unit=shiftList.unit, doFit=False)
              peak.setAnnotation('%.2f: ' % dist)
              duplicates[key] = peak 
 
            for j, peakDim in enumerate(peak.sortedPeakDims()):
              assignResToDim(peakDim, shiftDict[sortedDataDims[j]].resonance, tolerance=10.0)
 
          
  if progressBar:
    progressBar.close()
  
  return peakList        

def translateSpectrumUsingPeaks(referencePeak, translatePeak, dimMapping):
  """Translate spectrum referencing using a translatePeak from that
             spectrum and a referencePeak from another spectrum.  Use a dimMapping
             to specify what dimensions in the translate spectrum are mapped
             to what dimensions in the reference spectrum, and only these
             dimensions are translated.

  .. describe:: Input
  
  ccp.nmr.Nmr.Peak, ccp.nmr.Nmr.Peak, dict: int->int
  .. describe:: Output
  
  None
  """

  #referenceSpectrum = referencePeak.peakList.dataSource
  translateSpectrum = translatePeak.peakList.dataSource

  s = '%s:%s' % (translateSpectrum.experiment.name, translateSpectrum.name)
  for trnDim in dimMapping.keys():
    refDim = dimMapping[trnDim]
    #refDataDim = referenceSpectrum.findFirstDataDim(dim=refDim)
    refPeakDim = referencePeak.findFirstPeakDim(dim=refDim)
    trnDataDim = translateSpectrum.findFirstDataDim(dim=trnDim)
    trnPeakDim = translatePeak.findFirstPeakDim(dim=trnDim)
    trnDataDimRef = getPrimaryDataDimRef(trnDataDim)
    
    msg = 'translating %s dim %d dataDimRef refValue from %4.3f to ' % (s, trnDim, trnDataDimRef.refValue)
    trnDataDimRef.refValue += refPeakDim.value - trnPeakDim.value
    
    print '%s%4.3f' % (msg, trnDataDimRef.refValue)

def getPeakAnnotation(peak, noPeakAnnotationChar='', noPeakDimAnnotationChar='', joinChar='', doPeakDims=True):
  """
  Get annotation text for the peak.
  
  .. describe:: Input
  
  ccp.nmr.Nmr.Peak, String, String, String
  
  .. describe:: Output
  
  String
  """

  text = peak.annotation or noPeakAnnotationChar

  if doPeakDims:
    texts = []
    for peakDim in peak.sortedPeakDims():
      texts.append(peakDim.annotation or noPeakDimAnnotationChar)

    text = text + joinChar.join(texts)

  return text

def getSimplePeakAnnotation(peak, doChain = False):

  aList = []
  for peakDim in peak.sortedPeakDims():
    aText = ''
    bList = []
    chainCodes = {}

    for contrib in peakDim.peakDimContribs:
      bText = '?'
      resonance = contrib.resonance
      spinSystem = resonance.resonanceGroup

      if resonance.resonanceSet:
        residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
        code1Letter = residue.molResidue.chemComp.code1Letter or 'X'
        bText = '%d%s' % (residue.seqCode,code1Letter)

        if doChain:
          chainCodes[residue.chain.code] = True

      elif spinSystem and spinSystem.residue:
        residue = resonance.resonanceGroup.residue
        code1Letter = residue.molResidue.chemComp.code1Letter or 'X'
        bText = '%d%s' % (residue.seqCode,code1Letter)

        if doChain:
          chainCodes[residue.chain.code] = True

      elif spinSystem and spinSystem.residueProbs:
        resTexts = []
        resSeqs = []
        resCodes = set()
 
        for residueProb in spinSystem.residueProbs:
          if not residueProb.weight:
            continue
          
          residue = residueProb.possibility
          seq = residue.seqCode
          
          if doChain:
            chainCodes[residue.chain.code] = True
        
          resCode = residue.chemCompVar.chemComp.code1Letter or '?'
          resText = '%d?%s' % (seq, resCode)

          resTexts.append(resText)
          resSeqs.append('%d?' % seq)
          resCodes.add(resCode)
 
        if len(resCodes) == 1:
          bText = '/'.join(resSeqs) + resCodes.pop()
        else:
          bText = '/'.join(resTexts)

      elif spinSystem:
        bText = '{%d}' % spinSystem.serial

      bList.append(bText)

    if bList:
      for bText in bList:
        if bText == bList[0]:
          aText = bList[0]
        else:
          aText = '*'
          break

      if doChain:
        codes = chainCodes.keys()
        codes.sort()
        aText = '/'.join(codes) + aText
      aList.append(aText)

  simpleText = ''
  if aList:
    simpleText = aList[0]

  for aText in aList:
    if aText != aList[0]:
      simpleText = ''
      for aText in aList:
        simpleText += '%s' % aText
      break

  return simpleText

def getPeakSeqCodes(peak):

  seqCodes = set()
  for peakDim in peak.peakDims:
    for peakDimContrib in peakDim.peakDimContribs:
      if isinstance(peakDimContrib, PeakDimContribN):
        continue
      resonance = peakDimContrib.resonance
      resonanceSet = resonance.resonanceSet
      if resonanceSet:
        for atomSet in resonanceSet.atomSets:
          for atom in atomSet.atoms:
            residue = atom.residue
            seqCodes.add(residue.seqCode)

  return sorted(seqCodes)

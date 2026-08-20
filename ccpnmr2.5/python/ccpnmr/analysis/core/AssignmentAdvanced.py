
"""
======================COPYRIGHT/LICENSE START==========================

AssignmentAdvanced.py: Part of the CcpNmr Analysis program

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

from ccpnmr.analysis.core.AssignmentBasic import  assignResToDim, findConnectedSpinSystem, makeSeqSpinSystemLink
from ccpnmr.analysis.core.AssignmentBasic import  setResonanceTypeFromRefExp, findSpinSystem, addSpinSystemResonance
from ccpnmr.analysis.core import ExperimentBasic

from ccpnmr.analysis.core.PeakBasic import findPeaks, searchPeaks, pickPeak
from ccpnmr.analysis.core.UnitConverter import ppm2pnt, pnt2ppm
from ccpnmr.analysis.core.Util import getAnalysisDataDim

try:
  from memops.gui.MessageReporter import showWarning
except ImportError:
  from memops.universal.MessageReporter import showWarning

from math import sqrt

intraSpinSystemExperiments = ['H[{N|C}]_H.TOCSY','H[N]_H.TOCSY','H[C]_H.TOCSY',
                              'H[C[C]]','HCCH','HC_CH.TOCSY','HC_cNH.TOCSY','HNCAHA',
                              'H[N[CA[HA]]]','H[N[{CA|ca[C]}[H]]]','H{[N]+[HA]}','H[N[HB]]']

# intra HNCACO, HNCA/ HNCA/CB 

def assignSpecNonRootResonances(rootPeaks, targetPeakList, tolerances=None, diagTolerance=0.5,
                                waterMinPpm=4.88,waterMaxPpm=4.92, progressBar=None,
                                assignType=False, refExpsCO=None):
  """Descrn: 
     Inputs: 
     Output: 
  """
  peaks, dimMapping = pickAssignSpecFromRoot(rootPeaks, targetPeakList, tolerances=tolerances,
                                             progressBar=progressBar, pickNew=False)
  
  spectrum = targetPeakList.dataSource
  targetDataDim = None
  rootDataDims = set()
  for dataDim in spectrum.dataDims:
    if dataDim.dim in dimMapping.values():
      rootDataDims.add(dataDim.dim)
    else:
      targetDataDim = dataDim
  
  if targetDataDim is None:
    msg = 'No non-root dimension found for %s:%s'
    showWarning('Failure', msg % (spectrum.experiment.name, spectrum.name))
  
  resonances = assignDimNewResonances(peaks, targetDataDim.dim,
                                      diagTolerance=diagTolerance,
                                      waterMinPpm=waterMinPpm,
                                      waterMaxPpm=waterMaxPpm,
                                      assignType=assignType)
  if not refExpsCO:
    refExps, refExpsCO = ExperimentBasic.getSeqAssignRefExperiments(targetPeakList.root)
  
  refExperiment = targetPeakList.dataSource.experiment.refExperiment
  
  if refExperiment in refExpsCO:
    resonances2 = set([r for r in resonances if not r.resonanceGroup])
    
    rootSpinSystems = set()
    for peak in peaks:
      for peakDim in peak.peakDims:
        if peakDim.dim in rootDataDims:
          for contrib in peakDim.peakDimContribs:
            spinSystem = contrib.resonance.resonanceGroup

            if spinSystem:
              rootSpinSystems.add(spinSystem)
   
    if len(rootSpinSystems) == 1:
      rootSpinSystem = rootSpinSystems.pop()
      prevSpinSystem = findConnectedSpinSystem(rootSpinSystem, delta=-1)
      
      if not prevSpinSystem:
        prevSpinSystem = rootSpinSystem.nmrProject.newResonanceGroup()
        makeSeqSpinSystemLink(prevSpinSystem, rootSpinSystem, delta=1)
      
      resonances2.update(prevSpinSystem.resonances)
      prevSpinSystem.setResonances(resonances2)
      
def assignDimNewResonances(peaks, dimNum, diagTolerance=0.5, waterMinPpm=4.88,
                           waterMaxPpm=4.92,assignType=False):
  """Descrn: 
     Inputs: 
     Output: 
  """
    
  resonances = []
  for peak in peaks:
    peakDims = peak.sortedPeakDims()
    peakDim  = peak.findFirstPeakDim(dim=dimNum)
    ppm      = peakDim.value
    
    # Check for 1H diagonal peaks and water range
    if '1H' in peakDim.dataDimRef.expDimRef.isotopeCodes:
      if (ppm >= waterMinPpm) and (ppm <= waterMaxPpm):
        continue

      ignore = False
      for peakDim0 in peak.peakDims:
        if peakDim0 is not peakDim:
          if '1H' in peakDim0.dataDimRef.expDimRef.isotopeCodes:
            if abs( ppm - peakDim0.value) < diagTolerance:
              ignore = True
              break
      
      if ignore:
        continue  
      
    refExperiment = peak.peakList.dataSource.experiment.refExperiment
    
    if len(peakDim.peakDimContribs) < 1:
      resonance = assignResToDim(peakDim).resonance
      resonances.append(resonance)
      
      if refExperiment and (refExperiment.nmrExpPrototype.name in intraSpinSystemExperiments):
        contrib0 = None
        for peakDim0 in peak.peakDims:
          if (peakDim0 is not peakDim) and peakDim0.peakDimContribs:
            contrib0 = peakDim0.findFirstPeakDimContrib()
        
        if contrib0:
          resonance0 = contrib0.resonance
          if resonance0.resonanceGroup:
            addSpinSystemResonance(resonance0.resonanceGroup, resonance)

      if assignType:
        refExpDimRef = peakDim.dataDimRef.expDimRef.refExpDimRef
        if refExpDimRef:
          setResonanceTypeFromRefExp(resonance, refExpDimRef)
          
    else:
      resonance = peakDim.findFirstPeakDimContrib().resonance

      if assignType:
        refExpDimRef = peakDim.dataDimRef.expDimRef.refExpDimRef
        if refExpDimRef:
          setResonanceTypeFromRefExp(resonance, refExpDimRef)

      resonances.append(resonance)

  return resonances


def pickAssignSpecFromRoot(rootPeaks, targetPeakList, tolerances=None,
                           progressBar=None, pickNew=True, diagTolerance=None,
                           waterExclusion=None):
  """Descrn: Pick peaks in a spectrum based upon a root peak list
             and assign picked peaks to root resonances.
     Inputs: List of Nmr.Peaks, Nmr.PeakList, List of Floats (target dim order),
             memops.gui.ProgressBar, Boolean, Float or None,
             2-Tuple of Floats (min, max) or None
     Output: Nmr.Peaks, Dict of Int:Int (dim mapping root:target)
  """

  # TBD: perhaps a list of peaks rather than a set should be passed in
  rootPeakList = tuple(rootPeaks)[0].peakList
  project      = rootPeakList.root

  spectrum0 = rootPeakList.dataSource
  spectrum1 = targetPeakList.dataSource
  
  if not tolerances:
    tolerances = []
    for dataDim in spectrum1.sortedDataDims():
      tolerances.append(getAnalysisDataDim(dataDim).assignTolerance)
  
  expDimRefs0 = ExperimentBasic.getOnebondExpDimRefs(spectrum0.experiment)
  expDimRefs1 = ExperimentBasic.getOnebondExpDimRefs(spectrum1.experiment)
  
  if not expDimRefs0:
    msg = 'Spectrum %s:%s has no directly bonded dims.'
    showWarning('Failure', msg % (spectrum0.experiment.name, spectrum0.name))
    return

  if not expDimRefs1:
    msg = 'Spectrum %s:%s has no directly bonded dims.'
    showWarning('Failure', msg % (spectrum1.experiment.name, spectrum1.name))
    return
  
  # could really do with a generic get matching dataDims function
  
  # TBD try PeakBasic.getDataDimMapping()
  
  dimMapping = {}
  
  for expDimRef0, expDimRef1 in expDimRefs0:
    dataDim0  = spectrum0.findFirstDataDim(expDim=expDimRef0.expDim)
    dataDim1  = spectrum0.findFirstDataDim(expDim=expDimRef1.expDim)
    isotopes0 = expDimRef0.isotopeCodes
    isotopes1 = expDimRef1.isotopeCodes
    if  '1H' in isotopes0 or '1H' in isotopes1:
      break
  
  boundDims = set([dataDim0.dim, dataDim1.dim])
  
  for expDimRef2, expDimRef3 in expDimRefs1:
    
    isotopes2 = expDimRef2.isotopeCodes
    isotopes3 = expDimRef3.isotopeCodes
    dataDim2  = spectrum1.findFirstDataDim(expDim=expDimRef2.expDim)
    dataDim3  = spectrum1.findFirstDataDim(expDim=expDimRef3.expDim)
        
    if not (dataDim2 and dataDim3):
      continue
    
    if (isotopes0 == isotopes2) and (isotopes1 == isotopes3):
      dimMapping[dataDim0.dim] = dataDim2.dim
      dimMapping[dataDim1.dim] = dataDim3.dim
      break
    
    elif (isotopes0 == isotopes3) and (isotopes1 == isotopes2):
      dimMapping[dataDim0.dim] = dataDim3.dim
      dimMapping[dataDim1.dim] = dataDim2.dim
      break

  if not dimMapping:
    msg = 'Could not find equivalent dims in spectra %s:%s & %s:%s'
    data = (spectrum0.experiment.name, spectrum0.name,
            spectrum1.experiment.name, spectrum1.name)
    showWarning('Failure', msg % data)
    return
    
  if progressBar:
    progressBar.total = len(rootPeaks)
    progressBar.set(0)
    progressBar.open()
     
  fullRegion = []
  getDimRef = ExperimentBasic.getPrimaryDataDimRef
  targetDims = targetPeakList.dataSource.sortedDataDims()
  targetDims = [(getDimRef(dd), dd) for dd in targetDims]
  targetIsotopes = [','.join(ddr.expDimRef.isotopeCodes) for (ddr, dd) in targetDims]

  for dataDimRef, dataDim in targetDims:
    valueMax   = pnt2ppm(1,dataDimRef)
    valueMin   = pnt2ppm(dataDim.numPoints,dataDimRef)
    fullRegion.append([valueMin,valueMax])

  foundPeaks = []
  for peak in rootPeaks:
    
    peakRegion = list(fullRegion)
    rootIsotopes = []
    rootMapDims = set()
    values = []
    
    boundPeakDims = [pd for pd in peak.sortedPeakDims() if pd.dim in boundDims]
    
    for i, peakDim in enumerate(boundPeakDims):
      error = tolerances[i]
      value = peakDim.value
      j = dimMapping[peakDim.dim]-1
      peakRegion[j] = (value-error, value+error)
      
      rootIsotopes.append(targetIsotopes[j])
      rootMapDims.add(j)
      values.append(value)

    
    regions = [peakRegion,]
    if diagTolerance:
      for i, isotopes in enumerate(targetIsotopes):
        if (i not in rootMapDims) and (isotopes in rootIsotopes):
          value = values[rootIsotopes.index(isotopes)]
          
          # Chop in two at diagonal
          regions2 = []
          for region in regions:
            regionA = list(region)
            regionB = list(region)
            valueMin, valueMax = region[i]
            regionA[i] = (valueMin,value-diagTolerance)
            regionB[i] = (value+diagTolerance,valueMax)
            regions2.append(regionA)
            regions2.append(regionB)
          
          regions = regions2
           
    if waterExclusion:
      waterMin, waterMax = waterExclusion
      for i, isotopes in enumerate(targetIsotopes):
        if (i not in rootMapDims) and (isotopes in rootIsotopes):
          regions2 = []
    
          # Chop at water
          for region in regions:
            valueMin, valueMax = region[i]
            
            
            if valueMin < waterMin < valueMax:
              regionA = list(region)
              regionA[i] = (valueMin, waterMin)
              regions2.append(regionA)
              
              if valueMin < waterMax < valueMax:
                regionB = list(region)
                regionB[i] = (waterMax, valueMax)
                regions2.append(regionB)
              
            elif valueMin < waterMax < valueMax:
              regionB = list(region)
              regionB[i] = (waterMax, valueMax)
              regions2.append(regionB)
           
            else:
              regions2.append(region)
          
          regions = regions2
            
    
    peaks = []
    for region in regions:
      peaks.extend( searchPeaks([targetPeakList,], region) )
      
      if pickNew:
        peaks.extend( findPeaks(targetPeakList, region) )
    
    # possible to grow the region here if no peaks are found
     
    for peak2 in peaks:
      for intens in peak.peakIntensities:
        if str(intens.value) == 'inf':
          peaks.remove(peak2)
          peak2.delete()
          break
               
    foundPeaks.extend(peaks)
    
    for i, peakDim in enumerate(boundPeakDims):
      dim2 = dimMapping[peakDim.dim]
      resonances = [c.resonance for c in peakDim.peakDimContribs]
      tolerance = tolerances[i]
      
      for peak2 in peaks:
        peakDim2 = peak2.findFirstPeakDim(dim=dim2)
        if not peakDim2.peakDimContribs:
          for resonance in resonances:
            assignResToDim(peakDim2, resonance, tolerance=10*tolerance,
                           doWarning=False)
     
    if progressBar:
      progressBar.increment()
    
  if progressBar:
    progressBar.close()
  
  return foundPeaks, dimMapping


def initialiseAmideExpts(argServer,hsqc=None, tocsy=None, noesy=None):
  
  func = ExperimentBasic.getPrimaryDataDimRef
  
  xDim = 1

  peakSize = (0.09,0.49)

  ratio = peakSize[0]/ peakSize[1]
  tol   = peakSize[0]/5.0

  minPpmH =  6.0
  maxPpmH = 9.84

  waterRegion = [4.90,4.92]

  #tocsy = argServer.getSpectrum()
  
  #nhsqcPl = getBlankPeakList(nhsqc)
  #tocsyPl = getBlankPeakList(tocsy)
  #noesyPl = getBlankPeakList(noesy)
  
  noesyPl = argServer.getPeakList()
  noesy   = noesyPl.dataSource

  tocsyPl = argServer.getPeakList()
  tocsy   = noesyPl.dataSource
  
  dataDims = noesy.sortedDataDims()
  
  if not noesyPl.peaks:
    print "Picking new NOE peaks"
    wholeRegion   = [[pnt2ppm(dd.numPointsOrig,func(dd)), pnt2ppm(0,func(dd))] 
                     for dd in dataDims]
    excludeRegion = [[0,dd.numPointsOrig] for dd in dataDims]
 
    dataDimRef = func(dataDims[xDim])
    excludeRegion[xDim] = [ppm2pnt(waterRegion[0], dataDimRef),
                           ppm2pnt(waterRegion[1], dataDimRef)]
 
    findPeaks(noesyPl, wholeRegion, argServer.parent, [1,1,1], )

  if not tocsyPl.peaks:
    print "Picking new TOCSY peaks"
    wholeRegion   = [[pnt2ppm(dd.numPointsOrig,func(dd)), pnt2ppm(0,func(dd))] 
                     for dd in dataDims]
    excludeRegion = [[0,dd.numPointsOrig] for dd in dataDims]
 
    dataDimRef = func(dataDims[xDim])
    excludeRegion[xDim] = [ppm2pnt(waterRegion[0], dataDimRef),
                           ppm2pnt(waterRegion[1], dataDimRef)]
 
    findPeaks(tocsyPl, wholeRegion, argServer.parent, [1,1,1], )
 
  nhsqcPl = argServer.getPeakList()
  nhsqc   = nhsqcPl.dataSource
  
  noise = ExperimentBasic.getNoiseEstimate(nhsqc) * 2.0
 
  dataDims = nhsqc.sortedDataDims()
  dd0  = dataDims[0]
  dd1  = dataDims[1]
  
  ddr0 = ExperimentBasic.getPrimaryDataDimRef(dd0)
  ddr1 = ExperimentBasic.getPrimaryDataDimRef(dd1)
  
  print "Initial NOESY filter"
  amides = []
  
  allPeaks = list(noesyPl.peaks)
  allPeaks.extend(tocsyPl.peaks)
  
  for peak in allPeaks:
    peakDims = peak.sortedPeakDims()
    ppm0 = peakDims[0].value
    ppm2 = peakDims[2].value
    
    pnt0 = ppm2pnt(ppm0, ddr0)
    pnt1 = ppm2pnt(ppm2, ddr1)
    
    if ppm0 < minPpmH:
      peak.delete()
      continue
    
    if ppm0 > maxPpmH:
      peak.delete()
      continue
    
    if (pnt0-1 < 0) or (pnt0 > dd0.numPointsOrig):
      peak.delete()
      continue
    
    if (pnt1-1 < 0) or (pnt1 > dd1.numPointsOrig):
      peak.delete()
      continue
      
    height1 = nhsqc.block_file.getValue( (pnt0-1,pnt1) )
    height2 = nhsqc.block_file.getValue( (pnt0,pnt1-1) )
    height3 = nhsqc.block_file.getValue( (pnt0-1,pnt1-1) )
    height4 = nhsqc.block_file.getValue( (pnt0,pnt1) )
    
    if height1 < noise:
      peak.delete()
      continue
    if height2 < noise:
      peak.delete()
      continue
    if height3 < noise:
      peak.delete()
      continue
    if height4 < noise:
      peak.delete()
      continue

    if peak.peakList is noesyPl:
      amides.append( (peak, ppm0,ppm2) )
    
    peak.ppmH = ppm0
    peak.ppmN = ppm2

  print "Cluster %d amides" % len(amides)
  cluster = {}
  for i in range(len(amides)-1):
    if i and i % 100 == 0: 
      print i
    peak1, ppm0, ppm1 = amides[i]
    if cluster.get(peak1) is None:
      cluster[peak1] = [peak1]
  
    for j in range(i+1, len(amides)):
      peak2, ppm2, ppm3 = amides[j]
      if peak1 is peak2:
        continue
        
      if cluster.get(peak2) is None:
        cluster[peak2] = [peak2]
      
      if cluster[peak1] == cluster[peak2]:
        continue
              
      deltaH = ppm2-ppm0
      deltaN = ratio * (ppm3-ppm1)
      delta  = sqrt( (deltaH*deltaH) + (deltaN*deltaN))
    
      if (delta <=tol):
        cluster[peak1].extend(cluster[peak2])
        
        for peak3 in cluster[peak2]:
          cluster[peak3] = cluster[peak1]

  print "Remove isolated peaks"
  clusters2 = {}
  for peak in cluster.keys():
    c = cluster[peak]
    if len(c) < 2:
      peak.delete()
      del c
    else:
      clusters2[c[0]] = c
  
  clusters = clusters2.values()

  ss = {}
  centres = []
  print "Check for overlapped clusters"
  for peaks in clusters:
  
    p = peaks[0]
    print 'CLUSTER', p.ppmH, p.ppmN 
  
    for n in range(1):
      print 'Iteration', n
      f =  0.30
      tolF = tol * f

      cluster = {}
      cluster2 = {}
      M = 0
      for i in range(len(peaks)-1):

        peak1 = peaks[i]
        if cluster.get(peak1) is None:
          cluster[peak1] = [peak1]
          cluster2[peak1] = M
          M += 1
 
        for j in range(i+1, len(peaks)):
          peak2 = peaks[j]
          if peak1 is peak2:
            continue
 
          if cluster.get(peak2) is None:
            cluster[peak2] = [peak2]
            cluster2[peak2] = M
            M += 1
 
          if cluster2[peak1] == cluster2[peak2]:
            continue
 
          deltaH = peak1.ppmH-peak2.ppmH
          deltaN = ratio * (peak1.ppmN-peak2.ppmN)
          delta  = sqrt( (deltaH*deltaH) + (deltaN*deltaN))

          if delta <= tolF:
            cluster[peak1].extend(cluster[peak2])
            
            for peak3 in cluster[peak2]:
              cluster[peak3] = cluster[peak1]
              cluster2[peak3] = cluster2[peak1]
      
      cluster3 = []
      for i in range(M):
        cluster3.append([])

      for peak in peaks:
        cluster3[cluster2[peak]].append(peak)

      
      print '  F %3f' % (f),
      for i in range(M):
        N = float(len(cluster3[i]))
        if N > 1.0:
          aveH = 0.0
          aveN = 0.0
          for peak in cluster3[i]:
            aveH += peak.ppmH
            aveN += peak.ppmN
          
          aveH /= N
          aveN /= N
          hsqcPeak = pickPeak(nhsqcPl, (aveH, aveN), unit='ppm')
          
          centres.append( [hsqcPeak, aveH, aveN] )
          ss[hsqcPeak] = []
          print len(cluster3[i]), 
          
  print "Assign 15N HSQC"
  #assignAllNewResonances(peaks=nhsqcPl.peaks)
  #assignSpinSystemPerPeak(peaks=nhsqcPl.peaks)
  
  print "Assign NOESY & TOCSY" 
  for peak in allPeaks:
    minDist = tol
    best = centres[0][0]
    for hsqcPeak, aveH, aveN in centres:
      deltaH = peak.ppmH-aveH
      deltaN = ratio * (peak.ppmN-aveN)
      delta  = sqrt( (deltaH*deltaH) + (deltaN*deltaN))
      
      if delta < minDist:
        minDist = delta
        best =  hsqcPeak
      
    ss[best].append(peak)
    
  for hsqcPeak in ss.keys():
    peaks = ss[hsqcPeak]
    #propagatePeakAssignments(peaks, refPeak=hsqcPeak)
  
      
def assignSpinSystemPerPeak(argServer=None, peaks=None):
  assert argServer or peaks
  if not peaks:
    peaks = argServer.getCurrentPeaks()

  for peak in peaks:
    spinSystem = None
    resonances = []
    peakDims = peak.sortedPeakDims()
    for peakDim in peakDims:
      if len(peakDim.peakDimContribs) == 1:
        resonance = peakDim.peakDimContribs[0].resonance
        resonances.append( resonance )
    
    for resonance in resonances:
      if resonance.resonanceGroup:
        spinSystem = resonance.resonanceGroup
        
    for resonance in resonances:
      if spinSystem is None:
        spinSystem = findSpinSystem(resonance)
      addSpinSystemResonance(spinSystem, resonance)


def assignAllNewResonances(argServer=None, peaks=None):
  assert argServer or peaks
  if not peaks:
    peaks = argServer.getCurrentPeaks()
  # e.g. Assign all of an unassigned HSQC to different Hn and N resonances 

  for peak in peaks:
    for peakDim in peak.peakDims:
      if len(peakDim.peakDimContribs) < 1:
        assignResToDim(peakDim)

def getBlankPeakList(spectrum):

  peakList = None
  
  if spectrum.peakLists:
    for peakList0 in spectrum.peakLists:
      if len(peakList0.peaks) == 0:
        peakList = peakList0
        break
    
  if not peakList:
    peakList = spectrum.newPeakList()
    
    
  return peakList  

def possiblePeakAssigments(peak, shiftRanges=None, tolerances=None,
                           aliasing=True, findAssigned=False):
  """ Get list of resonance that make possible assignments for a peak
  NB could be made more sophisticated
  NB does  not (yet) deal with TOCSY or intraresidue
  """
  
  dataSource = peak.peakList.dataSource
  numDim = dataSource.numDim
  if shiftRanges is None:
    shiftRanges = [None] * numDim
  if tolerances is None:
    tolerances = estimateAssignmentTolerances(dataSource)
    tolerances = [None] * numDim
  
  
  matchResDict = {}
  for ii,peakDim in enumerate(peak.sortedPeakDims()):
    ll = findMatchingPeakDimShifts(peakDim, shiftRanges[ii], 
                                                   tolerances[ii], aliasing, 
                                                   findAssigned)
    matchResDict[peakDim.dim] = set(x.resonance for x in ll)
  
  
  boundDimNos = []
  for dim1,dim2 in ExperimentBasic.getOnebondDataDims(dataSource):
   boundDimNos.append((dim1.dim, dim2.dim))
  
  
  # NBNB Multiple runs to guard against one dimension being bound to two others.
  changing = True
  while changing:
    changing = False
    for dimNo,resonances in matchResDict.items():
      for dimNo1, dimNo2 in boundDimNos:
        if dimNo == dimNo1:
          set1 = set()
          set2 = set()
          for res1 in matchResDict[dimNo1]:
            for res2 in matchResDict[dimNo2]:
              if areResonancesBound(res1, res2):
                set1.add(res1)
                set2.add(res2)
          if (len(set1) < len(matchResDict[dimNo1]) or
              len(set2) < len(matchResDict[dimNo2])):
            changing = True
            matchResDict[dimNo1] = set1
            matchResDict[dimNo2] = set2
    
  #
  return [val for key,val in sorted(matchResDict.items())]
    
  
  


"""
======================COPYRIGHT/LICENSE START==========================

NoeRelaxation.py: Part of the CcpNmr Clouds program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
import math
from ccpnmr.c.Midge import Midge
from ccpnmr.analysis.core.ConstraintBasic import makeNmrConstraintStore, getFixedResonance
from ccpnmr.analysis.core.PeakBasic import getPeakDimPpm, findPeaks, searchPeaks
from ccpnmr.analysis.core.UnitConverter import ppm2pnt
from ccpnmr.analysis.core.AssignmentBasic import clearPeakDim, findMatchingShifts, assignResToDim, getResonanceAtomTuple
from ccp.api.nmr import Nmr, NmrConstraint

def matchPeaks(peak1, peak2, tolerances):

  n = len(tolerances)
  score = 0
  
  peakDims1 = peak1.sortedPeakDims()
  peakDims2 = peak2.sortedPeakDims()
  
  for i in range(n):
    ppm1 = getPeakDimPpm(peakDims1[i])
    ppm2 = getPeakDimPpm(peakDims2[i])
    d = abs(ppm1-ppm2)
    if d > tolerances[i]:
      return 0
    else:
      s = d/tolerances[i]
      score += s*s

  return score

def matchPeakDims(peakDims1, peakDims2, tolerances):

  n = len(tolerances)
  score = 0
  
  for i in range(n):
    ppm1 = getPeakDimPpm(peakDims1[i])
    ppm2 = getPeakDimPpm(peakDims2[i])
    d = abs(ppm1-ppm2)
    if d > tolerances[i]:
      return 0
    else:
      s = d/tolerances[i]
      score += s*s

  return score
  
def disambiguateNoesyPeaks(noesy2dPeakList,noesy3dPeakList,tocsy3dPeakList,hsqcPeakList):

  # assign HSQC resonances and spin systems
  # assign pick 3d tocsy F1 F3
  # assign tocsy F2
  # assign pick 3d noesy F1 F3

  tolerances = [0.04, 0.04, 0.4]
  t1 = tolerances[0]/2.0
  t2 = tolerances[1]/2.0

  spinSystems = []
  for peak in hsqcPeakList.peaks:
    for peakDim in peak.peakDims:
      if peakDim.dataDim.expDim.isAcquisition:
        resonance = peakDim.findFirstpeakDimContrib().resonance
        spinSystems.append(resonance.resonanceGroup)
        break
  
  ii = 0
  for spinSystem in spinSystems:
    resonances = spinSystem.resonances
  
    intraPeaks = []
    interPeaks = []
    anchors = []
    t3peaks = []
    n3peaks = []
    matched = []
    
    # find 3d tocsy and 3d noesy peaks for spin system
    for resonance in resonances:
      for contrib in resonance.peakDimContribs:
        peak = contrib.peakDim.peak
	
	if peak.peakList is noesy3dPeakList:
	  n3peaks.append(peak) 
	elif peak.peakList is tocsy3dPeakList:
	  t3peaks.append(peak)
	  
    # get intra list
    for pn in n3peaks:
      match = None
      score = 0
      for pt in t3peaks:
        s = matchPeaks(pt,pn,tolerances)
	if s and (s>score):
	  score = s
	  match = pt
	  
      if match and match.sortedPeakDims()[1].peakDimContribs:
        assignResToDim(pn.sortedPeakDims()[1],match.sortedPeakDims()[1].findFirstpeakDimContrib().resonance)
        intraPeaks.append(pn)
    
    # get unambiguous inter list
    for pn in n3peaks:
      if pn not in intraPeaks:
        peakDim = pn.sortedPeakDims()[1]
        if peakDim.peakDimContribs and peakDim.findFirstpeakDimContrib().resonance:
          interPeaks.append(pn)
        
	ppm = getPeakDimPpm(peakDim)
        position  = peakDim.position + (peakDim.numAliasing*peakDim.dataDimRef.dataDim.numPointsOrig)
	shifts = findMatchingShifts(peakDim.dataDimRef,position,tolerance=t1)
                
	if shifts and len(shifts) == 1:
	  assignResToDim(pn.sortedPeakDims()[1],shifts[0].resonance)
	  interPeaks.append(pn)
        else:
 	  shifts = findMatchingShifts(peakDim.dataDimRef,position,tolerance=t1)
          if not shifts:
            pass
            #assignResToDim(peakDim)
	    #interPeaks.append(pn)
            
          elif len(shifts) == 1:
	    assignResToDim(peakDim,shifts[0].resonance)
	    interPeaks.append(pn)
            
            
    peaks3d = intraPeaks
    peaks3d.extend(interPeaks)
    
    ii += 1
    print "SS", ii, spinSystem
    
    # pick assign 2d equivalents to intra and unambiguous 3d and find anchor point peaks
    for p3 in peaks3d:
      matched = {}
      f = 0.0
      n = 0
      pd1 = p3.sortedPeakDims()[0]
      pd2 = p3.sortedPeakDims()[1]
      ppm1 = getPeakDimPpm(pd1)
      ppm2 = getPeakDimPpm(pd2)
      
      region = [[ppm1-t1,ppm1+t1],[ppm2-t2,ppm2+t2]]
      n2peaks = searchPeaks([noesy2dPeakList,], region)
      n2peaks.extend( findPeaks(noesy2dPeakList, region))
      match = None
      score = 0
      
      if len(n2peaks) == 1:
        match = n2peaks[0]
      else:
        for p2 in n2peaks:
          s = matchPeakDims(p3.peakDims,p2.peakDims,tolerances[:2])
          if s and (s>score):
            score = s
            match = p2
          
      if match:
        anchors.append( (match,p3) )
        f += match.findFirstPeakIntensity().value/p3.findFirstPeakIntensity().value
        n += 1
        peakDimsM = match.sortedPeakDims()
        
        clearPeakDim(peakDimsM[0])
        clearPeakDim(peakDimsM[1])
        assignResToDim(peakDimsM[0],pd1.findFirstpeakDimContrib().resonance)
        assignResToDim(peakDimsM[1],pd2.findFirstpeakDimContrib().resonance)
        matched[p3] = 1

      region = [[ppm2-t1,ppm2+t1],[ppm1-t2,ppm1+t2]]
      n2peaks = searchPeaks([noesy2dPeakList,], region)
      n2peaks.extend(findPeaks(noesy2dPeakList, region))
      match = None
      score = 0
      
      if len(n2peaks) == 1:
        match = n2peaks[0]
      else:
        for p2 in n2peaks:
          s = matchPeakDims(p3.peakDims,p2.peakDims,tolerances[:2])
          if s and (s>score):
            score = s
            match = p2
          
      if match:
        peakDimsM = match.sortedPeakDims()
        anchors.append( (match,p3) )
        f += match.findFirstPeakIntensity().value/p3.findFirstPeakIntensity().value
        n += 1
        clearPeakDim(peakDimsM[0])
        clearPeakDim(peakDimsM[1])
        assignResToDim(peakDimsM[0],pd1.findFirstpeakDimContrib().resonance)
        assignResToDim(peakDimsM[1],pd2.findFirstpeakDimContrib().resonance)
        matched[p3] = 1
        
    if matched:
         
      f /= n
      
      # use the peaks that intially failed in the 2d picking 
           
      for p3 in peaks3d:
        if not matched.get(p3):
          pd1 = p3.sortedPeakDims()[0]
          pd2 = p3.sortedPeakDims()[1]
        
          r1 = pd1.findFirstpeakDimContrib().resonance
          r2 = pd2.findFirstpeakDimContrib().resonance
        
          ppm1 = getPeakDimPpm(pd1)
          ppm2 = getPeakDimPpm(pd2)
        
          region = [[ppm1-tolerances[0],ppm1+tolerances[0]],[ppm2-tolerances[1],ppm2+tolerances[1]]]
          peaks = searchPeaks([noesy2dPeakList,], region)
          peaks.extend( findPeaks(noesy2dPeakList, region) )
          for p2 in peaks:
            p2.annotation = '*'
            p2.findFirstPeakIntensity().value = p3.findFirstPeakIntensity().value * f
            assignResToDim(p2.sortedPeakDims()[0],r1)
            assignResToDim(p2.sortedPeakDims()[1],r2)
          
          region = [[ppm2-tolerances[0],ppm2+tolerances[0]],[ppm1-tolerances[1],ppm1+tolerances[1]]]
          peaks = searchPeaks([noesy2dPeakList,], region)
          for p2 in peaks:
            p2.annotation = '*'
            p2.findFirstPeakIntensity().value = p3.findFirstPeakIntensity().value * f
            assignResToDim(p2.sortedPeakDims()[0],r2)
            assignResToDim(p2.sortedPeakDims()[1],r1)
  
  return noesy2dPeakList
  resonances = list(noesy2dPeakList.root.resonances)
  N = len(resonances)
  
  print "Final 2d matches for %d resonances" % N
  
  for i in range(N-1):
    print "R", i
    for j in range(i+1,N):
      ppm1 = resonances[i].findFirstShift().value
      ppm2 = resonances[j].findFirstShift().value
      region = [[ppm1-t1,ppm1+t1],[ppm2-t2,ppm2+t2]]
      
      peaks = searchPeaks([noesy2dPeakList,], region)
      if not peaks:
        peaks = findPeaks(noesy2dPeakList, region)
        
      if len(peaks) == 1:
        for peakDim in peaks[0].peakDims:
          shifts = findMatchingShifts(peakDim.dataDimRef,position,tolerance=t1)
          if len(shifts) == 1:
            assignResToDim(peakDim,shifts[0].resonance)
        
  return noesy2dPeakList

def assignAmbigPeaks(peaks, resonances):

  tolerance = 0.015
  shifts = []
  
  N - len(resonances)
  for resonance in resonances:
    shifts.append(resonance.findFirstShift().value)
      
  for peak in peaks:
    for peakDim in peak.peakDims:
      ppm = pnt2ppm(peakDim.position,peakDim.dataDimRef)
      candidates = []
      for i in range(N):
        if abs(ppm-shift) <= tolerance:
          candidates.append(i)

      if len(candidates) == 1:
        resonance = resonances[candidates[0]]
        assignResToDim(peakDim,resonance)

def optimiseRelaxation(resonances,amat,tmix=60,sf=500,tcor=3,rleak=2,C13=1,N15=1,minConnectivity=1,maxIter=50):

  n = len(resonances) # number of H 'atoms'
  
  assert n == len(amat)
  assert n == len(amat[0])
 
  atom_types  = { 'CH3': 1, 'Haro': 3, 'HN': 4 }
  protonNumbs = { 'CH3': 3, 'Haro': 2, 'HN': 1, 'H': 1}

  # make blank matrices
  nhs   = n * [0]
  types = n * [0]
  atoms = n * [0]

  for i in range(n):
    typ      = resonances[i].name
    atoms[i] = typ
    nhs[i]   = protonNumbs[typ]
    types[i] = atom_types.get(typ, 0)
  
  symmetriseNoeMatrix(amat)
  removeOutliers(amat)  

  exclude={}
  for a in range(len(amat)):
    nonZero = 0
    for b in range(len(amat[a])):
      if amat[a][b] > 0:
        nonZero += 1
    if nonZero < minConnectivity:
      exclude[a] = 1
      
  amatNew=[]
  for a in range(len(amat)):
    if exclude.get(a) is None:
      new = []
      for b in range(len(amat[a])):
        if exclude.get(b) is None:
          new.append(amat[a][b])
      amatNew.append(new)
      
  for a in range(len(amat)-1,-1,-1):
    if exclude.get(a):
      resonances.pop(a)
      atoms.pop(a)
      nhs.pop(a)
      types.pop(a)

  print len(atoms), len(nhs), len(types), len(amat), len(amatNew)

  n = len(resonances) # number of atoms
  rmat  = n * [0]
  for i in range(n):
    rmat[i] = n * [0]

  amat = amatNew      
  #print amat
  num = 0
  for i in range(n):
    for j in range(n):
      if amat[i][j] != 0.0:
        num += 1
  
  fp = open('noe01.inp', 'w')
  fp.write('%6.1d\n' % (num))
  for i in range(n):
    amat[i][i] = nhs[i]
    for j in range(n):
      if amat[i][j] != 0.0:
        fp.write('%3.1d  %3.1d   %9.2e\n' % (i+1,j+1,amat[i][j]))
    #fp.write('\n')
  fp.close()

  fp = open('symbols.inp', 'w')
  for i in range(n):
    fp.write('%3.1d %-4s %d\n' % (i+1,resonances[i].name,nhs[i]))
  fp.close()
  

  m = Midge(nhs, types)
  err = m.run(amat, rmat, maxIter, sf, tmix, tcor, rleak, N15, C13)
  print 'error = %3.2e' % err

  """
  fp = open('amat.out', 'w')
  for i in range(n):
    for j in range(n):
      fp.write('%8.3e ' % amat[i][j])
    fp.write('\n')
  fp.close()

  fp = open('rmat.out', 'w')
  for i in range(n):
    for j in range(n):
      fp.write('%8.3e ' % rmat[i][j])
    fp.write('\n')
  fp.close()
  """

  PI = math.pi
  GH = 2.6752e4
  if N15:
    GN = -2.7126e3
  else:
    GN = 1.9340e3
  HB = 1.05459e-27
  CONST = GH*GH*GH*GH * HB*HB
  tc = 1.0e-9 * tcor
  wh = 2.0 * PI * sf * 1.0e6
  wn = wh * GN / GH
  j0 = CONST * tc
  j1 = CONST * tc / (1.0 + wh*wh*tc*tc)
  j2 = CONST * tc / (1.0 + 4.0*wh*wh*tc*tc)
  #jself = 6.0*j2 + 3.0*j1 + j0
  jcross = 6.0*j2 - j0

  project             = resonances[0].root
  constraintHead      = project.newNmrConstraintStore(nmrProject=resonances[0].topObject)
  distConstraintList  = NmrConstraint.DistanceConstraintList(constraintHead)
  for i in range(n-1):
    for j in range(i+1, n):
      if ((amatNew[i][j] != 0) and (rmat[i][j] != 0.0)):
        dist = 1.0e8 * ((10*abs(rmat[i][j]/(jcross*nhs[j])))**(-1.0/6.0))
        if dist > 6.0:
          continue
        fixedResonanceI = getFixedResonance(constraintHead,resonances[i])
        fixedResonanceJ = getFixedResonance(constraintHead,resonances[j])
        #fp.write('%3d %3d %3d %4s %3d %4s  %5.2f\n' % (i+1, j+1, i+1, atoms[i], j+1, atoms[j], d))
        constraint = NmrConstraint.DistanceConstraint(distConstraintList, weight=1, targetValue=dist, upperLimit=dist+(dist/2.0), lowerLimit=dist-(dist/2.0), error=dist/5)
        item = NmrConstraint.DistanceConstraintItem(constraint, resonances=[fixedResonanceI,fixedResonanceJ])
        
  print "Midge for CcpNmr Done"
  return distConstraintList, resonances

def getResonancesFromPeaks(peaks): 

  n = len(peaks[0].peakDims)
  resonances = []
  resonanceDict = {}
  
  for peak in peaks:
    peakDims = peak.sortedPeakDims()
    for i in range(n):
      for contrib in peakDims[i].peakDimContribs:
        resonance = contrib.resonance
        if resonance.resonanceSet:
          if resonanceDict.get(resonance) is None:
            resonances.append(resonance)
            resonanceDict[resonance] = 1

  return resonances

def removeOutliers(matrix):

  routes = []
  cluster = {}
  clusters = []
  N = len(matrix)
  M = len(matrix[0])
  for i in range(N):
    routes.append([])
    for j in range(M):
      if j != i:
        if matrix[i][j] > 0:
          routes[i].append(j)

  c = 0
  for i in range(N-1):
    
    if cluster.get(i) is None:
      cluster[i] = c
      members = routes[i]
    else:
      continue
      
    for j in members:
      cluster[j] = c
      for k in routes[j]:
        if cluster.get(k) is None:
          cluster[k] = c
          members.append(k)
     
    clusters.append( members )
    c += 1
           
  if c > 1:
    best = 0
    maxLen = len(clusters[best])
    for i in range(1,c):
      tryLen = len(clusters[i])
      if tryLen > maxLen:
        maxLen = tryLen
        best = i
    
    for i in range(N):
      if (cluster.get(i) is None) or (cluster[i] != best):
        for j in range(M):
          matrix[j][i] = 0
          matrix[i][j] = 0
         
  return matrix         
        
def symmetriseNoeMatrix(matrix):

  N = len(matrix)
  M = len(matrix[0])
  for i in range(N-1):
    for j in range(i+1,M):
      if not matrix[i][j]:
        matrix[i][j] = matrix[j][i]
      elif not matrix[j][i]:
        matrix[j][i] = matrix[i][j]
      else:
        v = matrix[i][j] + matrix[j][i]
        v /= 2
        
        matrix[i][j] = v
        matrix[j][i] = v

  return matrix

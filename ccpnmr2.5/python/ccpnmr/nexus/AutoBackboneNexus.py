
"""
======================COPYRIGHT/LICENSE START==========================

AutoBackboneNexus.py: Part of the CcpNmr Nexus program

Copyright (C) 2003-2010 Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the author: tjs23@cam.ac.uk
=======================================================================
===========================REFERENCE START=============================
===========================REFERENCE END===============================

"""               
                         
from memops.gui.ScrolledDensityMatrix import ScrolledDensityMatrix
from memops.gui.ScrolledGraph import ScrolledGraph
from ccpnmr.analysis.core.ExperimentBasic import getSpectraByType, getSpectrumIsotopes
from ccpnmr.analysis.core.ChemicalShiftBasic import getShiftsChainProbabilities
from ccpnmr.analysis.core.MoleculeBasic import getLinkedResidue
from ccpnmr.analysis.core.AssignmentBasic import findConnectedSpinSystem
from ccpnmr.nexus.NexusBasic import getSpinSystemInterIntraResonances

import Tkinter

from math import exp, log, sqrt

root2pi = 2.506628274631

from numpy import matrix, zeros, multiply, where, linalg, ones
from random import randint, seed, random

import time, numpy

seed(time.time())

PROLYL = ('Pro','Hyp')

def matrixRefine(C,A,R, matrixIter, progressBar):

  if progressBar:
    progressBar.setText('Initial Optimisation')
    progressBar.total = 100
    progressBar.set(0)
    progressBar.open()
    progressBar.parent.update_idletasks()

  m, n = A.shape
  # shape = ss = rows, rr = cols
  
  Rt = R.T
  Ct = C.T
  B = A
  
  j = matrixIter/100
  for i in range(matrixIter):

    # i-1, i+1
    F = multiply((C*B*Rt) + (Ct*B*R), B)
    
    # On spin system
    for s in xrange(m):
      s3 = F[s].sum()
      if s3:
        F[s] /= s3
    
    # On residue
    F = F.T
    for s in xrange(n):
      s3 = F[s].sum()
      if s3:
        F[s] /= s3
    
    B = F.T
    #B /= F.sum()
    
    if (i%j == 0) and i and progressBar:
      progressBar.increment()
    
    #graph.title = title + ' - Iteration %d of %d' % (i+1, iterations)
    #graph.update(B.tolist())


  if progressBar:
    progressBar.close()

  return B

def autoBackboneNexus(chain, spinSystems, shiftList,
                      peakLists, iterations=25,
                      varianceC=0.09, varianceH=0.04,
                      keepExisting=False, progressBar=None):

  #ss2 = [(ss.residue.seqCode, ss) for ss in spinSystems if ss.residue]
  #ss2.sort()
  #spinSystems = [x[1] for x in ss2]

  A, R, spinSystems, residues = getInitialAssignMatrix(chain, spinSystems, peakLists,
                                                       shiftList, keepExisting,
                                                       progressBar)
  if not iterations:
    return
  
  matrixIter = 400
  
  m, n = A.shape
      
  cMatrix = getConnectivityMatrix(spinSystems, peakLists, sqrt(varianceC),
                                  sqrt(varianceH), progressBar)

  C = matrix(cMatrix)
  C /= C.sum()

  #root = Tkinter.Tk()
  #availHeight = int(root.winfo_screenheight())
  #root.geometry('%dx%d' % (availHeight,availHeight))
  #root.grid_columnconfigure(0, weight=1)
  #root.grid_rowconfigure(0, weight=1)

  #xLabels = [(i+1) for i in range(m)]
  #yLabels = ['%3d %s' % (r.seqCode, r.chemCompVar.chemComp.code1Letter) for r in residues]
  #title = 'Assignment Matrix Refinement'
  #boxSize = (availHeight - 100.0)/len(yLabels)
  #graph = ScrolledDensityMatrix(root, matrix=M.tolist(), boxSize=boxSize, maxVal=None,
  #                              xLabels=xLabels, yLabels=yLabels, borderWidth=1,
  #                              borderColor='grey', zoom=1.0, title=title,
  #                              xAxisLabel='Spin Systems', yAxisLabel='Residues')

  #graph.grid(row=0, column=0, sticky='nsew')

  #return {}

  # Nexus matrix method

  B = matrixRefine(C,A,R,matrixIter,progressBar)

  # Resolve

  maxS = B.max()
  assign = {}
  scores = []
  for s in xrange(m):
    for r in xrange(n):
      scores.append((A[s,r]*B[s,r],s,r))
  
  scores.sort()
  scores2 = []
  
  assign = {}
  while scores:
    v,s,r = scores.pop()
    assign[r] = s 
    i = 0
    scores2 = scores[:]
    scores2.reverse()
    for v2, s2, r2 in scores2:
      if i > 1:
        break
      if v2/maxS < 0.001:
        break
    
      if r2 == r:
        i += 1
    
    scores = [x for x in scores if (x[1] != s) and (x[2] != r)]
  
  # Monte carlo
  
  A /= A.max()
  B /= B.max()
  C /= C.max()

  m, n = A.shape
  for r in xrange(n):
    if r not in assign:
      assign[r] = None
  
  sc = 0.0
  for r in assign:
    s = assign[r]
    if s is None:
      continue
    
    sc += A[s,r]
    r2 = r-1
    s2 = assign.get(r2)
    
    if s2 is not None:
      sc += C[s2,s]
  
  e = 10
  ensemble = []
  for i in xrange(e):
    ensemble.append((sc, assign.copy()))

  if progressBar:
    progressBar.setText('Refinement')
    progressBar.total = iterations
    progressBar.set(0)
    progressBar.open()
    progressBar.parent.update_idletasks()

  n2 = n-1
  nR = range(n)
  scBest = sc
  assignBest = assign.copy()
  
  multi = 1000
  i =  multi * iterations
  I = float(i)
  t = i/e
  w = 5
  indices = range(n)
  
  if nR > w:
   for r in nR[:-w]:
     for r1 in range(r,r+w):
       s = assignBest[r1]
       if s is None:
         break
       
       sc3 = 0.0
     
       s2 = assignBest.get(r1-1)
       if s2 is None:
         break

       s3 = assignBest.get(r1+1)
       if s3 is None:
         break

       sc3 = C[s2,s] + C[s,s3]
       if sc3 <= 1.0:
         break
     
     else:
       for r1 in range(r,r+w):  
         if r1 in indices:
           indices.remove(r1)
           n2 -=1
 
  
  while i:
    ensemble.sort()
    
    if i and (i % t == 0) and (len(ensemble) > 2):
      ensemble.pop(0)
          
    if i and (i % multi == 0):
      if i < I/2:
        indices = range(n)
        n2 = n-1 
                         
      if progressBar:
        progressBar.increment()
    
    if ensemble[-1][0] > scBest:
      scBest = ensemble[-1][0]
      assignBest = ensemble[-1][1].copy()
      print 'Best Score:', scBest
    
    randint(1,len(ensemble))
    scW, assignW = ensemble[0]
    sc, assign = ensemble[randint(1,len(ensemble)-1)]
  
    rA = indices[randint(0,n2)]    
    rB = rA 
    while rA == rB:
      rB = indices[randint(0,n2)]
    
    assignB = assign.copy() 
    sB = assignB[rB]
    assignB[rB] = assignB[rA]
    assignB[rA] = sB
    
    sc2 = 0.0
    for r in nR:
      s = assignB[r]
      if s is None:
        continue
      
      sc2 += A[s,r]
    
      r2 = r-1
      s2 = assignB.get(r2)
    
      if s2 is not None:
        sc2 += C[s2,s]
            
    accept = scW-sc2
    if accept <= 0:
      ensemble[0] = (sc2, assignB)
      
    elif random() < (0.00001/(accept*accept)):
      ensemble[0] = (sc2, assignB)
      
    else:
      ensemble[0] = (sc, assign.copy())
      i -= 1

  rAssign = {}
  
  if progressBar:
    progressBar.close()
  
  for r in assignBest:
    s = assignBest[r]
     
    if s is not None:
      ss = set([s,])
      residue = residues[r]
      rAssign[residue] = [1.0, [spinSystems[s],]]
  
      for sc2, assign2 in ensemble:
        s3 = assign2[r]
      
        if (s3 is not None) and (s3 not in ss):
          ss.add(s3)
      
          if len(ss) < 3:
            rAssign[residue][1].append(spinSystems[s3])

      #rAssign[residue][0] /= len(ss)

  # Scale by connectivity score
  conns = [1.0] * n
  for i in range(n):
    residue = residues[i]
    prev = getLinkedResidue(residue, 'prev')
    next = getLinkedResidue(residue, 'next')

    connPrev = 0.5
    connNext = 0.5
    
    if prev and (prev.ccpCode in PROLYL):
      connPrev = 0.8

    else:
      spinSystemsA = rAssign.get(prev)
      spinSystemsB = rAssign.get(residue)
 
      if spinSystemsA and spinSystemsB:
        spinSystemA = spinSystemsA[1][0]
        spinSystemB = spinSystemsB[1][0]
        connPrev = getSpinSystemLikelihood(spinSystemA,
                                           spinSystemB,
                                           peakLists)
    
    if next and (next.ccpCode in PROLYL):
      connNext = 0.8
      
    else:
      spinSystemsB = rAssign.get(residue)
      spinSystemsC = rAssign.get(next)
 
      if spinSystemsB and spinSystemsC:
        spinSystemB = spinSystemsB[1][0]
        spinSystemC = spinSystemsC[1][0]
        connNext = getSpinSystemLikelihood(spinSystemB,
                                           spinSystemC,
                                           peakLists)
         
    conns[i] = connPrev*connNext
     
  
  for i, residue in enumerate(residues):
    if rAssign.get(residue):
      #rAssign[residue][0] = conns[i]
      rAssign[residue][0] *= conns[i]
  
  return rAssign

def getInitialAssignMatrix(chain, spinSystems, peakLists, shiftList,
                           keepExisting=False, progressBar=None):

  if progressBar:
    progressBar.setText('Preparing Spin System Typing')
    progressBar.total = len(spinSystems)
    progressBar.set(0)
    progressBar.open()
    progressBar.parent.update_idletasks()

  peakListDict = {}
  for peakList in peakLists:
    peakListDict[peakList] = True

  allResidues = chain.sortedResidues()
  nResidues = len(allResidues)
  
  prolyl = []
  residues = []
  prevResidueDict = {}
  
  rMatrix = []
  for i, residue in enumerate(allResidues):
    row = [0.0] * nResidues
    
    if i+1 < nResidues:
      row[i+1] = 1.0
    
    rMatrix.append(row)
  
    prevResidueDict[residue] = getLinkedResidue(residue, 'prev')
  
    if residue.ccpCode in PROLYL:
      prolyl.append(i)
    else:
      residues.append(residue)
     
  prolyl.reverse()
  for row in rMatrix:
    for i in prolyl:
      del row[i]
  
  for i in prolyl:
    del rMatrix[i]
  
  nResidues = len(residues)
  nSpinSystems = len(spinSystems)
  
  #t0 = time.time()
  
  aMatrix = [[0] * nResidues for x in range(nSpinSystems)]

  for i, spinSystem in enumerate(spinSystems):
    #if hasattr(spinSystem, 'codeScoreDict'):
    #  del spinSystem.codeScoreDict
    
    if keepExisting:
      residueA = spinSystem.residue
      ccpCodeA = spinSystem.ccpCode
      
      if residueA and (residueA.chain is chain):
        for j, residue in enumerate(residues):
          if residue is residueA:
            score = 1.0
          else:
            score = 0.0
        
          aMatrix[i][j] = score
        continue
      
      elif ccpCodeA:
        for j, residue in enumerate(residues):
          if residue.ccpCode == ccpCodeA:
            score = 1.0
          else:
            score = 0.0
        
          aMatrix[i][j] = score
        continue
      
    shifts = []
    for resonance in spinSystem.resonances:
      for contrib in resonance.peakDimContribs:
        peak = contrib.peakDim.peak
        
        if peakListDict.get(peak.peakList):
          shift = resonance.findFirstShift(parentList=shiftList)
 
          if shift:
            shifts.append(shift)
            break
     
    scores = getShiftsChainProbabilities(shifts, chain)
    
    prevSpinSystem = findConnectedSpinSystem(spinSystem, delta=-1)
    if spinSystem.residue and not prevSpinSystem:
      residueB = getLinkedResidue(spinSystem.residue, linkCode='prev')
      nmrProject = spinSystem.topObject
      prevSpinSystem = nmrProject.findFirstResonanceGroup(residue=residueB)
     
    prevScores = None
    if prevSpinSystem:
      prevShifts = []
      for resonance in prevSpinSystem.resonances:
        for contrib in resonance.peakDimContribs:
          peak = contrib.peakDim.peak
        
          if peakListDict.get(peak.peakList):
            shift = resonance.findFirstShift(parentList=shiftList)
 
            if shift:
              prevShifts.append(shift)
              break
      
      prevSpinSystem.shifts = prevShifts
      prevScores = getShiftsChainProbabilities(prevShifts, chain) 
 
    for j, residue in enumerate(residues):
      aMatrix[i][j] = scores[residue.ccpCode] or 1e-150
      prevResidue = prevResidueDict.get(residue)
      if prevResidue and prevScores:
        aMatrix[i][j] *= prevScores[prevResidue.ccpCode] or 1e-150
      else:
        aMatrix[i][j] *= aMatrix[i][j]
 
    if progressBar:
      progressBar.increment()

  #print "Time", time.time() - t0
    
  A = matrix(aMatrix)
  R = matrix(rMatrix)
  
  # Normalise per Spin System
  for s in xrange(A.shape[0]):
    A[s] /= A[s].max()
  
  # Normalise per Residue
  A = A.T
  for s in xrange(A.shape[0]):
    A[s] /= A[s].max()
  A = A.T
  
  A /= A.max()
  
  if progressBar:
    progressBar.close()
  
  return A, R, spinSystems, residues

  
def getConnectivityMatrix(spinSystems, peakLists, varianceC=0.09,
                          varianceH=0.04, progressBar=None):


  if progressBar:
    progressBar.setText('Preparing Spin System Connectivity')
    progressBar.total = len(spinSystems)
    progressBar.set(0)
    progressBar.open()
    progressBar.parent.update_idletasks()

  nSpinSystems = len(spinSystems)
  
  cMatrix = [[0] * nSpinSystems for x in range(nSpinSystems)]

  for i, spinSystemA in enumerate(spinSystems):
    for j, spinSystemB in enumerate(spinSystems):
      p = getSpinSystemLikelihood(spinSystemA, spinSystemB, peakLists,
                                  varianceC, varianceH)
      
      cMatrix[i][j] = p
  
    s = sum(cMatrix[i])
    
    if s:
      cMatrix[i] = [x/s for x in cMatrix[i]]
  
    if progressBar:
      progressBar.increment()
  
  if progressBar:
    progressBar.close()
  
  return cMatrix
 
  
def getSpinSystemLikelihood(spinSystemA, spinSystemB, peakLists,
                            sigmaC=0.3, sigmaH=0.2, limit=3.0, cache={}):

  data = getSpinSystemInterIntraResonances(spinSystemA, peakLists)
  # pRes, iRes, pPpm, iPpm, pTypes, iTypes = data 
  intraIsotopes = [r.isotopeCode for r in data[1]]
  intra = data[3]
  intraTypes = data[5]

  data = getSpinSystemInterIntraResonances(spinSystemB, peakLists)
  interIsotopes = [r.isotopeCode for r in data[0]]
  inter = data[2]
  interTypes = data[4]
  
  ss2C      = sigmaC*sigmaC*2
  #sroot2piC = sigmaC*root2pi
  ss2H      = sigmaH*sigmaH*2
  #sroot2piH = sigmaH*root2pi
  
  p = None
  
  for i, ppm0 in enumerate(inter):
    interType = interTypes[i]
    interIsotope = interIsotopes[i]
    if interIsotope == '15N':
      continue
    
    if interIsotope == '1H':
      ss2 = ss2H
      #sroot2pi = sroot2piH
      sMax = limit*sigmaH
    else:
      ss2 = ss2C
      #sroot2pi = sroot2piC
      sMax = limit*sigmaC
  
    for j, ppm1 in enumerate(intra):
      intraType = intraTypes[j]
    
      if interIsotope != intraIsotopes[j]:
        continue
    
      if intraType:
        if intraType != interType:
          continue

      d = ppm0-ppm1
      e = (d*d)/ss2
      q = exp(-e) # /sroot2pi
 
      if abs(d) > sMax:
        q = 0.0
        
        if not intraType:
          continue

      if p is None:
        p = q
      else:
        p *= q         
  
      if p == 0.0:
        break
    
    if p == 0.0:
      break
  
  return p or 0.0
  

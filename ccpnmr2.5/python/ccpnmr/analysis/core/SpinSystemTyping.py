
"""
======================COPYRIGHT/LICENSE START==========================

SpinSystemTyping.py: Part of the CcpNmr Analysis program

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
import sys

from math import sqrt, log, exp
from random import randint, random
from ccpnmr.analysis.core.ChemicalShiftBasic import getShiftsResidueProbability, getAtomProbability

olc = {'Ala':'A','Cys':'C','Asp':'D','Glu':'E','Phe':'F','Gly':'G',
       'His':'H','Ile':'I','Lys':'K','Leu':'L','Met':'M','Asn':'N',
       'Gln':'Q','Arg':'R','Ser':'S','Thr':'T','Val':'V','Trp':'W',
       'Tyr':'Y','Asx':'B','Glx':'Z','Hyp':'O','Pro':'P',None:'?'
       }

notVisible = {'Lys':["H''","HZ1","HZ2","HZ3"],
              'Arg':["H''","HH11","HH12","HH21","HH22","CZ"],
              'Ser':["H''","HG",],
              'Cys':["H''","HG",],
              'His':["H''","HD1","CG"],
              'Thr':["H''","HG1",],
              'Tyr':["H''","HH","CG"],
              'Phe':["H''","CG"],
              'Glu':["H''","HE2",],
              'Asp':["H''","HD2",],
              'Trp':["H''","CG","CD2","CE2"],
             }
       
atomNameDict = {}

class Blank:

  def __init__(self, project, serial):

    self.root    = project
    self.project = project
    self.serial  = serial
    self.seqId   = serial
    self.ccpCode = None
    self.residue = None
    self.className = 'blank'
    self.sstTypes = []
    self.resonances = []
    self.shifts     = []

def getSpinSystemTypes(residues, spinSystems, preserveAssign, isotopes=('1H','13C'),
                       shiftList=None, graph=None, numSteps=100000,
                       numBest=20, progressBar=None):

  residues = list(residues)

  print "Typing spin systems"
  project  = residues[0].root

  N = len(residues)
  M = len(spinSystems)

  for spinSystem in spinSystems:
    spinSystem.codeScoreDict = {}
    spinSystem.shifts = getSpinSystemShifts(spinSystem, shiftList, isotopes)

  if M > N:
    j = 1
    for i in range(N, M):
      residues.append(Blank(project, j+residues[-1].seqId))
      j += 1
  
  elif N > M:
    j = 1
    for i in range(M, N):
      spinSystems.append(Blank(project, j+spinSystems[-1].serial))
      j += 1

  atomTypes = []
  for isotope in isotopes:
    atomTypes.append(isotope[-1])

  ccpCodes0 = set()
  ccpCodeCount = {}
  for residue in residues:
    ccpCode = residue.ccpCode
    exclude = notVisible.get(residue.ccpCode) or ["H''",]
    
    if atomNameDict.get(ccpCode) is None:
      atomNames = getAtomNames(residue, atomTypes)
      if atomNames:
        atomNameDict[ccpCode] = atomNames

    if not ccpCode:
      continue
    
    ccpCodeCount[ccpCode] = ccpCodeCount.get(ccpCode, 0) + 1
    ccpCodes0.add(ccpCode)

  if progressBar:
    progressBar.total = len(spinSystems)
    progressBar.setText('Preparing chemical shift lookup')
    progressBar.set(0)

  scoreList = []

  typeScores = {}
  startingTypes = {}
  for ss in spinSystems:
    
    typeScores[ss] = {}
    startingTypes[ss] = None
    if not hasattr(ss, 'codeScoreDict'):
      ss.codeScoreDict = {}

    for ccpCode in ccpCodes0:
      score = getSpinSystemScore(ccpCode, ss)
      scoreList.append((score, ss, ccpCode))

    if progressBar:
      progressBar.increment()

  scoreList.sort()
  while scoreList:

    score, ss, ccpCode = scoreList.pop()
    startingTypes[ss] = ccpCode
    
    for i in range(len(scoreList)-1,-1,-1):
      if scoreList[i][1] is ss:
        scoreList.pop(i)

    if ccpCodeCount[ccpCode] < 2:
      for i in range(len(scoreList)-1,-1,-1):
        if scoreList[i][2] == ccpCode:
          scoreList.pop(i)

    else:
     ccpCodeCount[ccpCode] -= 1

  if progressBar:
    progressBar.total = 100
    progressBar.setText('Monte Carlo search')
    progressBar.set(0)
      

  ensemble = searchPosterior(spinSystems, startingTypes, graph=graph,
                             numBest=numBest, numSteps=numSteps,
                             preserveAssign=preserveAssign,
                             progressBar=progressBar)

  ensemble.reverse()
  cc0 = ensemble[0][1]
  
  for j in range(10):
    
    p, cc = ensemble[j]
    seq = ''.join([olc.get(cc[ss]) or '?' for ss in spinSystems])
    
    print p, seq
    
    for ss in spinSystems:   
      ccpCode = cc[ss]
      typeScores[ss][ccpCode] = typeScores[ss].get(ccpCode, 0) + 1
      
  ccpCodes = atomNameDict.keys()
  for ss in spinSystems:
    for ccpCode in ccpCodes:
      if typeScores[ss].get(ccpCode) is not None:
        score = getSpinSystemScore(ccpCode, ss) # possible - in ensemble
      else:
        score = -30.0 # very unlikely indeed - not in ensemble
      
      typeScores[ss][ccpCode] = score

  print "Done typing spin systems"
  
  return typeScores, cc0
 
def searchPosterior(spinSystems, cc0, graph=None, numSteps=100000, numBest=20,
                    preserveAssign=0, v=False, progressBar=None):

  from time import time

  udtA = time() - 1.51
  udtB = udtA - 1.51

  if progressBar:
    progressBar.open()
    progressBar.parent.update_idletasks()
    progressBar.total = 100
    progressBar.set(0)

  pBest = 0.0

  ensemble = []
  p0 = generatePosterior(cc0, spinSystems, preserveAssign)
  for j in range(numBest):
    ensemble.append( (p0, cc0.copy()) )

  
  q = numSteps/100

  i = 0
  iBest = 0
  while i < numSteps:
    i += 1
    num = 1  
      
    ensemble.sort()  
    p0, cc0 = ensemble[0]
      
    cc = getNewClassifications(cc0, num)
    p  = generatePosterior(cc, spinSystems, preserveAssign)
    accept = 0
    
    if p >= p0:
      accept = 1
       
    else:
      d = p/p0
      r = random()
     
      if r < d:
        accept = 1

    if accept:
      iBest = i
      passes = 0
      ensemble[0] = (p, cc)
      ensemble.sort()
      
      if ensemble[-1][0] > pBest:
        pBest = ensemble[-1][0]
        udt2 = time()
 
        if graph and (udt2-udtA > 1.5):
          udtA = udt2
          updateGraph(graph, spinSystems, ensemble[-1][1])
      
        if progressBar and (udt2-udtB > 0.2):
          udtB = udt2
          
          seqs = []
          for p1, cc1 in ensemble:
            seq = ''.join([ olc.get(cc1[ss]) or '?' for ss in spinSystems ])
            seqs.append(seq)
    
          progressBar.updateSequences(i, seqs)

        
      if v:
        known = ''.join([olc.get(ss.ccpCode) or '?' for ss in spinSystems])
        data = []
        for ss in spinSystems:
          shifts = []
          
        
        scores = ''.join([str(10+int(getSpinSystemScore( cc0[ss],ss ))) for ss in spinSystems if cc0[ss]])
        correct = ''
 
        for j in range(len(known)):
          if known[j] == current[j]:
            correct += known[j]
          else:
            correct += '-'
 
        print i, p0
        print scores
        print current
        print correct
        print known

    else:
      k = randint(1,numBest-1)
      p, cc = ensemble[k]
      ensemble[0] = (p, cc.copy())

    if progressBar and (i % q == 0):
      progressBar.increment()
      
  pBest, ccBest = ensemble[-1]
  if graph:
    updateGraph(graph, spinSystems, ccBest)
    
  if progressBar:
    seqs = []
    for p1, cc1 in ensemble:
      seq = ''.join([ olc.get(cc1[ss]) or '?' for ss in spinSystems ])
      seqs.append(seq)
    
    progressBar.updateSequences(iBest, seqs)
  
  if v:
    known = ''.join([olc.get(ss.ccpCode) or '?' for ss in spinSystems])
    current = ''.join([olc.get(ccBest[ss]) or '?'  for ss in spinSystems])
    correct = ''
    scores = ''.join([str(10+int(getSpinSystemScore( ccBest[ss],ss ))) for ss in spinSystems if ccBest[ss]])

    for j in range(len(known)):
      if known[j] == current[j]:
        correct += known[j]
      else:
        correct += '-'

    print i, pBest
    print scores
    print current
    print correct
    print known
  
  return ensemble

def getSpinSystemShifts(ss, shiftList, isotopes):
  
  shifts = []
  for resonance in ss.resonances:
    if resonance.isotopeCode in isotopes:
      if shiftList:
        shift = resonance.findFirstShift(parentList=shiftList)
      elif resonance.shifts:
        shift = resonance.findFirstShift()
      else:
        shift = None
      
      if shift:
        shifts.append(shift)

  return shifts
  

def updateGraph(graph, spinSystems, cc):

  dataSet = []
  for ss in spinSystems:
    ccpCode = cc[ss]
    if ccpCode:
      dataSet.append( [ss.serial, getSpinSystemScore(ccpCode, ss) or 0.0] )

  graph.update(dataSets=[dataSet,])
  #graph.parent.update_idletasks()

def generatePosterior(classifications, spinSystems, preserveAssign=0):

  likelihood = generateLikelihood(classifications, spinSystems, preserveAssign)
  prior = generatePrior(classifications)
  
  posterior = likelihood * prior

  return posterior

def generateLikelihood(cc, spinSystems, preserveAssign):

  T = len(spinSystems);
  p = 0.0
  
  for ss in spinSystems:
    if preserveAssign and ss.ccpCode and (ss.ccpCode not in ['Asx','Glx']):
      if cc[ss] == ss.ccpCode:
        p += 0.0
      else:
        p += -30.0
        
    elif cc[ss] is None:
    
      T -= 1
    else:
      p += getSpinSystemScore(cc[ss], ss)
  
  if T > 0:
    p /= float(T)
  else:
    p = 0.0
  
  return exp(p)

def generatePrior(classifications):

  # Flat prior - what else?
  # Sequence entropy?
  
  return 1.0

def getNewClassifications(cc, N):

  sss = cc.keys()
  L   = len(sss)-1
  
  new = {}
  for ss in sss:
    new[ss] = cc[ss]
  
  for i in range(N):
    j = randint(0,L)
    k = j
    while new[sss[k]] == new[sss[j]]:
      k = randint(0,L)
    
    swap        = new[sss[k]] 
    new[sss[k]] = new[sss[j]]
    new[sss[j]] = swap
    
  return new

def getAtomShiftScore(project, atomName, ccpCode, shift, normalise=True):

  symbol = shift.resonance.isotope.chemElement.symbol
  if atomName[:len(symbol)] != symbol:
    return 0.0

  return getAtomProbability(project, ccpCode, atomName, shift.value, normalise=normalise)

def getSpinSystemScore(ccpCode, spinSystem):

  if not hasattr(spinSystem, 'codeScoreDict'):
    spinSystem.codeScoreDict = {}

  elif spinSystem.codeScoreDict.get(ccpCode):
    return spinSystem.codeScoreDict[ccpCode]
  
  p = getShiftsResidueProbability(spinSystem.shifts, ccpCode) or 1e-12
  
  p = log(p)
  spinSystem.codeScoreDict[ccpCode] = p

  return p

def getAtomNames(residue, atomTypes=('H','C','N')):

  atomNames = []
  exclude   = notVisible.get(residue.ccpCode) or ["H''",]

  if (residue.className != 'blank') and (residue.molResidue.linking == 'middle'):
    for atom in residue.atoms:
      if atom.name[0] in atomTypes:
        if atom.name in exclude:
          continue
 
        if atom.atomSet:
          if len(atom.atomSet.atoms) > 1:
            if atom is atom.atomSet.findFirstAtom():
              atomNames.append(atom.name)
          else:
            atomNames.append(atom.name)
        else:
          atomNames.append(atom.name)
 
  return atomNames

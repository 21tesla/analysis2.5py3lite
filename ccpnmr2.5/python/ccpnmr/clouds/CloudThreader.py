"""
======================COPYRIGHT/LICENSE START==========================

CloudThreader.py: Part of the CcpNmr Clouds program

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
import cPickle
from os.path import exists, isfile, isdir
from os import listdir, path
from math import sqrt, log, exp
from random import randint, random
from memops.general.Io import loadProject
from memops.gui.MessageReporter import showWarning

from ccpnmr.analysis.core.SpinSystemTyping import getSpinSystemTypes
from ccpnmr.analysis.core.ChemicalShiftBasic import getSpinSystemChainProbabilities
from ccpnmr.analysis.core.AssignmentBasic  import makeResonanceGuiName
from ccpnmr.clouds.CloudBasic import code1LetterDict

#from ccpnmr.analysis.core.StructureBasic import makeStructureDictFromRoughPdb
olc = code1LetterDict

logDict = {}
nmrRefStoreDict = {}
chemAtomShiftDict = {}
global J
J = 0
global inter
global intra
inter = None

def scoreAssignment(chain, clouds, spinSystems=None, graph=None):

  global inter
  project  = chain.root
  if not spinSystems:
    spinSystems = list(project.currentNmrProject.resonanceGroups)
  
  assignment0 = {}
  for ss in spinSystems:
    ss.resonanceCoords = {}
    if ss.residue:
      assignment0[ss.residue] = ss
    
  residueScores = {}
  if inter is None:
    intra, inter = getDistDistributions()
  
  for residue in chain.residues:
    ss = assignment0.get(residue)
    ll = getResidueLikelihood(residue, ss, assignment0, clouds, inter) or 0.0
    residueScores[residue] = ll
  
  if graph:
    updateGraph(graph, chain.residues, assignment0, 0)
    
  return 1.0, residueScores, assignment0

def cloudThreader(chain, clouds, spinSystems=None, nSteps=3000, graph=None, preserveAssign=0, progressBar=None):

  project  = chain.root
  residues = list(chain.residues)

  if not spinSystems:
    spinSystems = list(project.currentNmrProject.resonanceGroups)

  for residue in residues:
    residue.spinSystemScoreDict = {}

  for ss in spinSystems:
    ss.resonanceCoords = {}
    
  score, scores, assignment = searchPosterior(project, residues, spinSystems,
                                              clouds, nSteps, graph=graph,
                                              preserveAssign=preserveAssign,
                                              progressBar=progressBar )

  return score, scores, assignment


def searchPosterior(project, residues, spinSystems, clouds, nSteps=3000, graph=None, verbose=0, preserveAssign=0, progressBar=None):

  chain = residues[0].chain
  for spinSystem in spinSystems:
    for resonance in spinSystem.resonances:
      if not resonance.resonanceSet:
        pass # resonance.delete()
      else:
        resonance.name == 'r'

  sequence = ''
  for residue in chain.residues:
    sequence += olc.get(residue.ccpCode, 'X')

  for residue in residues:
    residue.atomDict = {}
    for atom in residue.atoms:
      residue.atomDict[atom.name] = atom

  print "Generating intial assignment"
  typeScores, assignment0 = getInitialAssignment(residues, spinSystems, preserveAssign, graph=graph)
  print "Reading distance distributions"
  intra, inter = getDistDistributions()

  p0 = generatePosterior(typeScores, residues, assignment0, clouds, intra, inter, 0, preserveAssign)
  p  = p0

  
  if progressBar:
    progressBar.parent.update_idletasks()
    progressBar.total = 100
    progressBar.label.set('Searching cloud')
    progressBar.set(0)
    progressBar.open()

  q = max(1, nSteps/100)
  
  if graph:
    updateGraph(graph, chain.residues, assignment0, 0)

  print "Generate posterior"
  for i in range(nSteps):

    """if i < 20:
      n = 20
    elif i < 30:
      n = 15
    elif i < 50:
      n = 10
    elif i < 80:
      n = 5
    elif i < 100:
      n = 1
    else:"""
    n = 1

    assignment = swapSpinSystem(typeScores, assignment0, clouds, inter, n, i)  
    p = generatePosterior(typeScores, residues, assignment, clouds, intra, inter, i, preserveAssign)

    print i, n , p, p0
    success = 0
    if p > p0:
      success = 1
    
    elif p < p0:
      delta = p0 - p
      r = random()
      
      if i > 500:
        e = exp(delta*-5e8)
        if e > r:
          print e, r
          success = 1
   
      
    if success: 
      p0 = p
      print p0
      assignment0 = assignment
      
      if graph:
        updateGraph(graph, chain.residues, assignment0, i)
      
      if verbose:
        found = 0
        assign = ''
        typed  = ''
        foundCodes = ''
        scores = ''
        for residue in chain.residues:
          ss = assignment0[residue]
          if hasattr(residue, 'likelihood'):
            ll = residue.likelihood
          else:
            ll = 0.0
          score = int(ll+10)
          if score > 9:
            scores += '*'
          else:
            scores += str(score)
          if ss.residue:
            if ss.residue is residue:
              found +=1
              foundCodes += olc.get(residue.ccpCode, 'X')
              #if ll == -10.0:
              #  print residue.seqCode, residue.ccpCode
              #  getResidueLikelihood(residue, ss, assignment0, cloud, inter, v=1)
              
            else:
              foundCodes += '-'

            if ss.residue.ccpCode == residue.ccpCode:
              typed += olc.get(residue.ccpCode, 'X')
            else:
              typed += '-'

            assign += olc.get(ss.residue.ccpCode, 'X')
          else:
            assign += '-'
            foundCodes += '-'
            typed += '-'

        print i, p, p0
        print 'Found: %d' % found
        print scores
        print foundCodes
        print typed
        print sequence
        print assign
    
    if progressBar and (i % q == 0):
      progressBar.increment()
    
  residueScores = {}
  for residue in chain.residues:
    ss = assignment0[residue]
    ll = getResidueLikelihood(residue, ss, assignment0, clouds, inter) or 0.0
    residueScores[residue] = ll
   
  return p0, residueScores, assignment0

def updateGraph(graph, residues, assignment, step):

  title  = 'Clouds Threader Fitting Scores - step %d' % step
  xLabel = 'Residue'
  yLabel = 'Score'

  dataSet = []
  dataSet2 = []
  for residue in residues:
    ss = assignment[residue]
    residue2 = ss.residue
    if residue2 is None:
      check = -5.0
    elif residue2 is residue:
      check = 0.0
    else:
      check = -10.0
  
    if hasattr(residue, 'likelihood'):
      ll = residue.likelihood
    else:
      ll = 0.0
    dataSet.append(  [residue.seqCode, check] )
    dataSet2.append( [residue.seqCode, ll] )

  graph.update(dataSets=[dataSet,dataSet2],xLabel=xLabel,yLabel=yLabel,title=title)
  graph.parent.update_idletasks()

def printResidueScores(typeScores, residues, assignment, clouds, interDistrib):

  for residue in residues:
    spinSystem = assignment[residue]
    l = getResidueLikelihood(residue, spinSystem, assignment, clouds, interDistrib) or 0.0
    s = getSpinSystemScore(typeScores, residue, spinSystem) or 0.0
    print "%d %s %f %f" % (residue.seqCode,residue.ccpCode, s, l)
 
 
def generatePosterior(typeScores, residues, assignment, clouds, intra, inter, i, preserveAssign=0):

  
  prior = generatePrior(typeScores, residues, assignment)
  likelihood = generateLikelihood(residues, assignment, clouds, intra, inter, preserveAssign)

  p = prior * likelihood 
  # no normalisation (yet)
  
  return p

def getSpinSystemScore(typeScores, residue, spinSystem):
  
  score = typeScores[spinSystem].get(residue.ccpCode) or -10.0
  
  return score
  
  if p > 0.0:
    score = Log(p)
  else:
    score = -10.0

  return score

def generatePrior(typeScores, residues, assignment):

  # how well do the shifts of the resonances in a spin system match a residue
  # could the shift possibly be in the spin system
  # no conditional dependance in the intial instance

  score = 0.0
  N = 0

  for residue in residues:
    spinSystem = assignment[residue]
    score2 = getSpinSystemScore(typeScores, residue, spinSystem)
    
    N += 1
    score += score2

  if N > 0:
    return exp(score/len(residues))

  else:
    return 0.0


def getResidueLikelihood(residue1, spinSystem1, assignment, clouds, interDistrib, v=0, floor=0.0025):

  chain   = residue1.chain
  amide1  = getResonanceCoords(clouds, spinSystem1, 'H') or getResonanceCoords(clouds, spinSystem1, 'HD')
  alpha1  = getResonanceCoords(clouds, spinSystem1, 'HA')
  betas1  = getResonanceCoords(clouds, spinSystem1, 'HB')
  coords1 = [amide1,alpha1,betas1]
  coords2 = None 
  
  if not hasattr(chain, 'lookupSeqID'):
    chain.lookupSeqID = {}
    for residue in chain.residues:
      chain.lookupSeqID[residue.seqId] = residue
  lookupSeqID = chain.lookupSeqID
  
  r2 = []
  
  #i = 0
  i = 2
  N = 0
  p = 0.0
  badPrev = 0
  #for delta in (-3,-2,-1,1,2,3):
  for delta in (-1,1):
    p2 = 0.0
    N2 = 0
    residue2 = lookupSeqID.get(residue1.seqId+delta)
    if residue2 is None:
      i += 1
      continue
  
    spinSystem2 = assignment.get(residue2)
    if spinSystem2 is None:
      i += 1
      continue
  
    r2.append( '%d%s' % (residue2.seqCode, residue2.ccpCode))
  
    dist = []
    
    amide2  = getResonanceCoords(clouds, spinSystem2, 'H') or getResonanceCoords(clouds, spinSystem2, 'HD')
    alpha2  = getResonanceCoords(clouds, spinSystem2, 'HA')
    betas2  = getResonanceCoords(clouds, spinSystem2, 'HB')
    coords2 = [amide2,alpha2,betas2]
    
    for coord1 in coords1:
      for coord2 in coords2:
        if coord1 and coord2:
          dist.append(getEnsembleCoordsDist(coord1,coord2))
        else:
          dist.append(None)

    for j in range(9):
      if dist[j] is not None:
        key = str(round(dist[j],1))  
        q   = interDistrib[i][j].get(key, 0.0)
          
        if q:
          q =  max(q, floor)
        else:
          q = floor

 	p2 += log(q)
 	N2 += 1
                    
    if N2 and p2 > -9.9:
      p2 /= float(N2)
      N2 = 1
      
    N += N2
    p += p2
    i += 1
  
  
  
  if N == 0:
    if v:
      print "barf"
      print r2
      print coords1
      print coords2
      
    
    residue1.likelihood = -10.0
    return -10.0
  
  out = p/float(N)  
  
  if v:
    print "ok", out
    print r2
    print coords1
    print coords2
  
  residue1.likelihood = out 
  return out


def Log(q):
  
  key = str(q)
  if logDict.get(key) is not None:
    return logDict[key]
  else:
    v = log(q)
    logDict[key] = v
    return v


def generateLikelihood(residues, assignment, clouds, intraDistrib, interDistrib, preserveAssign=0):

  p = 0.0
  N = 0
  for residue1 in residues:
      
    spinSystem1 = assignment[residue1]
    if preserveAssign and spinSystem1.residue:
      if spinSystem1.residue is residue1:
        p2 = 0.0
      else:
        p2 = -10.0
    else:
      p2 = getResidueLikelihood(residue1, spinSystem1, assignment, clouds, interDistrib)
    
    p += p2
    N += 1
    
  if N == 0:
    return 0.0
 
  return exp(p/len(residues))

def swapSpinSystem(typeScores, assignment, clouds, interDistrib, n, iteration):
  
  """
  if iteration < 500:
    threshold = -1.3
  else:
  """
  threshold = -5.9
 # threshold  = -3.0
  threshold2 = -4.3
  newAssignment = {}
  residues = assignment.keys()
  for residue in residues:
    newAssignment[residue] = assignment[residue]
  
  N = len(residues)
  
  for i in range(n):
    j = randint(0,N-1)
    spinSystem1 = newAssignment[residues[j]]
    
    h = 0
    while h<25 and getSpinSystemScore(typeScores, residues[j],spinSystem1) > threshold and getResidueLikelihood(residues[j],spinSystem1, newAssignment, clouds, interDistrib) > threshold2:
      j = randint(0,N-1)
      spinSystem1 = newAssignment[residues[j]]
      h += 1

    k = j
 
    h = 0
    while (k == j) and h<25:
      k = randint(0,N-1)
      
      spinSystem2 = newAssignment[residues[k]]
      if k != j:
        if getSpinSystemScore(typeScores, residues[k],spinSystem1) < threshold:
          if getSpinSystemScore(typeScores, residues[j],spinSystem2) < threshold:
            k = j
        
      
        if (k != j) and (iteration <500):
          s0  = getResidueLikelihood(residues[j],spinSystem1, newAssignment, clouds, interDistrib)
          s0 += getResidueLikelihood(residues[k],spinSystem2, newAssignment, clouds, interDistrib)
          s1  = getResidueLikelihood(residues[k],spinSystem1, newAssignment, clouds, interDistrib)
          s1 += getResidueLikelihood(residues[j],spinSystem2, newAssignment, clouds, interDistrib)
          if s1 <= s0:
            k = j
      
      h += 1
    
    #spinSystem2 = newAssignment[residues[k]]
    newAssignment[residues[j]] = spinSystem2
    newAssignment[residues[k]] = spinSystem1

  return newAssignment

def getInitialAssignment(residues, spinSystems, preserveAssign, graph=None, progressBar=None):

  print "Gen initial assignment"
  chain      = residues[0].chain
  shiftList = chain.root.findFirstNmrMeasurementList(className='ShiftList')
  codeList   = {}
  typeScores = {}
  for ss in spinSystems:
    print '  ', ss
    scores = getSpinSystemChainProbabilities(ss, chain, dshiftList)
    
    typeScores[ss] = {}
    for ccpCode in scores.keys():
      typeScores[ss][ccpCode] = scores[ccpCode]
  
  typeScores, bestGuess = getSpinSystemTypes(residues, spinSystems,
                                             preserveAssign, ('1H','13C'),
                                             shiftList=shiftList,
                                             graph=graph,
                                             progressBar=progressBar)
  
  sss = list(spinSystems)
  assignment = {}
  for residue in residues:
    i  = randint(1,len(sss)) - 1
    n  = 0
    ss = sss[i]
    while (n < len(sss)) and (typeScores[ss][residue.ccpCode] < -1.3):
      i  = randint(1,len(sss)) - 1
      n += 1
    
    if typeScores[ss][residue.ccpCode] < -1.3:
      n  = 0
      while (n < len(sss)) and (typeScores[ss][residue.ccpCode] < -2.3):
        i  = randint(1,len(sss)) - 1
        n += 1
    
    assignment[residue] = sss.pop(i) 
  
  """
  for ss in sss:
    if ss.residue:
      assignment[ss.residue] = ss
  
  # typeScores[ss][ccpCode] = score
  
  # take most likely neighbouring spinSystems
  """
  
  return typeScores, assignment


def fitSideChains(scores, assignment, clouds, threshold, shiftList=None):

  # assignment[residue] = spinSystem

  # Get prior form ato chemical shift match
  
  # Search with likelihood from clouds
  
  print 'Calculating distance distributions'
  
  intraDistrib = readIntraDistribution('intraDistribs001.txt')
  
  print 'Monte Carlo...'

  # Get only the good assignments
  assign = {}
  for residue in assignment.keys():
    score = scores.get(residue)
    if score is not None:
      if score > threshold:
        assign[residue] = assignment[residue]

  residues = assign.keys()
  project  = residues[0].root
  
  if shiftList is None:
    shiftList = project.currentNmrProject.findFirstMeasurementList(className='ShiftList')
  
  for residue in residues:
    spinSystem = assign[residue]
    shifts = []
    for resonance in spinSystem.resonances:
      shift = resonance.findFirstShift(parentList=shiftList)
      if shift:
        shifts.append(shift)
    
    if shifts:
      distrib = intraDistrib.get(residue.ccpCode) or {}
      fitSideChain(shifts, residue, clouds, distrib)

def fitSideChain(shifts, residue, clouds, distrib, useAssignNames=False):
  
  ccpCode    = residue.ccpCode
  atomScores = {}
  project = residue.root

  from ccpnmr.analysis.core.SpinSystemTyping import getAtomNames, getAtomShiftScore

  atomSets = {}
  for shift in shifts:
    
    """scores[shift] = None
    
    name = None
    if shift.resonance.assignNames:
      name = shift.resonance.assignNames
      
    elif shift.resonance.name:
      name = shift.resonance.name  
      
    name = shift.resonance.name
    
    if name == 'HN':
      name = 'H'"""
 
    done = {}
    atomScores[shift] = {}

    for atom in residue.atoms:
      if not atom.chemAtom.waterExchangeable:
        atomSet = atom.atomSet
        if atomSet and (done.get(atomSet) is None):
          atomSetName       = atomSet.name
          atomName          = atom.name
          done[atomSet]     = True
          atomSets[atomSet] = True
 
          if atom.chemAtom.elementSymbol == shift.resonance.isotope.chemElement.symbol:
            score = getAtomShiftScore(project, atomSetName, ccpCode, shift)
            if not score:
              score = getAtomShiftScore(project, atomName, ccpCode, shift)
            
            if score:
              atomScores[shift][atomSetName] = log(score)
            else:
              atomScores[shift][atomSetName] = None
 
          else:
            atomScores[shift][atomSetName] = None
  
  atomNames0 = [(ass.atoms[0].name, ass.name) for ass in atomSets.keys()] # getAtomNames(residue)
    
  best    = -20000.0
  bestM   = None
  passes  = 0
  J       = 0  
  i       = 0
  
  print "Fitting residue %d %s" % (residue.seqCode, residue.ccpCode)
  #print atomNames0
  
  while passes < (len(atomNames0)*len(atomNames0)):
  
    i += 1
    if bestM is None:
      bestM = {}
      shifts1 = list(shifts)
      
      for atomName, atomSetName in atomNames0:
        if shifts1:
          j = randint(0,len(shifts1)-1)
          bestM[atomName] = shifts1.pop(j)
        else:
          bestM[atomName] = None
      
      mapping = bestM.copy()
      
    else:
      mapping = bestM.copy()
      
      if len(atomNames0) < 2:
        break
      
      j = J
      J += 1
      if J >= len(atomNames0):
        J = 0
      
      k = j
      while k == j:
        k = randint(0, len(atomNames0)-1)    

      swap = mapping[atomNames0[j][0]]
      mapping[atomNames0[j][0]] = mapping[atomNames0[k][0]]
      mapping[atomNames0[k][0]] = swap
       
    prior = 0.0
    for atomName, atomSetName in atomNames0:
      shift = mapping[atomName]
      if shift:
        s = atomScores[shift][atomSetName]
        if s:
          prior += s
        else:
          prior += -50.0
          
    likelihood = 0.1 * getAtomCloudScores(mapping, ccpCode, clouds, distrib)
          
    score = prior + likelihood

    if score > best:
      #print i, score, likelihood, prior
      passes = 0
      best   = score
      bestM  = mapping
      
    else:
      passes += 1
 
  # do some provisional assignment
  if bestM:
    for atomName in bestM.keys():
      shift = bestM[atomName]
      if shift:
        resonance = shift.resonance
        resonance.setName(atomName)
        
        atomNames4 = []
        if resonance.resonanceSet:
          for atomSet in resonance.resonanceSet.atomSets:
            for atom4 in atomSet.atoms:
              atomNames4.append(atom4.name) 
        
        if atomNames4 and (atomName not in atomNames4):
          #print i, best
          print 'Error %5.5s %3.3f - %s' % (atomName, shift.value, makeResonanceGuiName(resonance))
          for atom4 in residue.atoms:
            if not atom4.chemAtom.waterExchangeable:
              print '     %3.3f %5.5s %.3f' % (shift.value, atom4.atomSet.name, atomScores[shift][atom4.atomSet.name] or -1000.0)
 
  return bestM

def getAtomCloudScores(mapping, ccpCode, clouds, distrib):

  score = 0.0
  atomNames = mapping.keys()
  
  for atomName1 in atomNames:
    shift1  = mapping[atomName1]
    if shift1:
      ensemble1 = []
      for cloud in clouds:
        coords = []
        coord = getResonanceCoord(cloud, shift1.resonance)
        if coord:
          coords.append(coord)

        if coords:
          coord = averageCoord(coords)
          ensemble1.append(coord)

      for atomName2 in atomNames:
        if atomName1 != atomName2:
          shift2 = mapping[atomName2]
          if shift2:
            ensemble2 = []
            for cloud in clouds:
              coords = []
              coord = getResonanceCoord(cloud, shift2.resonance)
              if coord:
                coords.append(coord)

              if coords:
                coord = averageCoord(coords)
                ensemble2.append(coord)
 
            if ensemble1 and ensemble2:
              dist = getEnsembleCoordsDist(ensemble1,ensemble2)
              key  = str(round(dist,1))
              val  = distrib[atomName1][atomName2].get(key)
              if val:
                score += log(val)
              else:
                score += -10.0

  return score


def getDistDistributions(dataDirName='/home/tjs234/ccpn/recoord/', intraFileName='intraDistribs001.txt', interFileName='interDistribs001.txt'):

  if not( exists(intraFileName) and isfile(intraFileName) and exists(interFileName) and isfile(interFileName)):
    calcDistDistributions(dataDirName, intraFileName, interFileName)

  inter = readInterDistribution(interFileName)
  intra = {} #readIntraDistribution(intraFileName)

  return intra, inter


def makeStructureDictFromXml(path, fileName):

  print "Reading %s" % fileName
  project = XmlIO.loadProject(fileName, showWarning=showWarning)

  molSystem = project.findFirstMolSystem()
  structure = molSystem.findFirstStructureEnsemble()
  print "Got structure %s" % structure

  model = structure.findFirstModel()

  # set chain and residue dicts
  dict = {}
  for chain in structure.coordChains:
    dd1 = dict[chain.code] = {}
    for residue in chain.residues:
      dd2 = dd1[residue.seqId] = {}
      dd2['resName'] = residue.residue.ccpCode
      dd2['atoms'] = {}
  
  # fill in atom coordinates
  offset = 0
  coordinates = model.coordinates
  for atom in structure.orderedAtoms:
    next = offset + 3
    residue = atom.residue
    dict[residue.chain.code][residue.seqId]['atoms'][atom.name] = (
     list(coordinates[offset:next])
    )
    offset = next
  
  del project
  
  return dict   

def calcDistDistributions(dirName, intraFileName='intraDistribs001.txt', interFileName='interDistribs001.txt'):

  intra = {}
  data  = []
  posnLabels = ['i-3','i-2','i-1','i+1','i+2','i+3','other']
  atomLabels = ['amide','alpha','betas']
  for i in range(7):
    data.append([])
    for j in range(9):
      data[i].append({})
  
  files = listdir(dirName)

  for fileName in files:
    path = '%s%s/%s.xml' % (dirName, fileName, fileName)
    path2 = '%s%s/' % (dirName, fileName)
    if not exists(path):
      continue
  
    dict =  makeStructureDictFromXml(path2, path)
    #model = dict.keys()[0]

    for chain in dict.keys():
      chainDict = dict[chain]
      for residue1 in chainDict.keys():
        ccpCode = chainDict[residue1]['resName']
        atoms   = chainDict[residue1]['atoms'].keys()
        if intra.get(ccpCode) is None:
          intra[ccpCode] = {}
          for atom1 in atoms:
            intra[ccpCode][atom1] = {}
            for atom2 in atoms:
              intra[ccpCode][atom1][atom2] = {}
          
        amide1 = getAmideCoord(chainDict[residue1]['atoms'])
        alpha1 = getAlphaCoord(chainDict[residue1]['atoms'])
        betas1 = getBetasCoord(chainDict[residue1]['atoms'])
        coords1 = [amide1,alpha1,betas1]

        if not (amide1 and alpha1):
          continue

        for residue2 in chainDict.keys():
          delta = residue2 - residue1
          if delta == 0:
            continue
          
          amide2 = getAmideCoord(chainDict[residue2]['atoms'])
          alpha2 = getAlphaCoord(chainDict[residue2]['atoms'])
          betas2 = getBetasCoord(chainDict[residue2]['atoms'])
          coords2 = [amide2,alpha2,betas2]

          if not (amide1 and alpha2):
            continue

          if abs(delta) > 3:
            i = 6
          elif delta < 0:
            i = delta + 3
          else:
            i = delta + 2

          for j in range(9):
            coord1 = coords1[j/3]
            coord2 = coords2[j%3]
            if coord1 and coord2:
              dist = getCoordsDist(coord1, coord2)
              key  = str(round(dist, 1))
              data[i][j][key] = data[i][j].get(key, 0) + 1
          
        for atom1 in atoms:
          coord1 = chainDict[residue1]['atoms'][atom1]
          for atom2 in atoms:
            if atom1 == atom2:
              continue
            
            coord2 = chainDict[residue1]['atoms'][atom2]
            dist = getCoordsDist(coord1, coord2)
            key = str(round(dist, 1))
            if not intra[ccpCode].get(atom1):
              intra[ccpCode][atom1] = {}
              
            if not intra[ccpCode][atom1].get(atom2):
              intra[ccpCode][atom1][atom2] = {}
              
            intra[ccpCode][atom1][atom2][key] = intra[ccpCode][atom1][atom2].get(key, 0) + 1

  interFile = open(interFileName, 'w')

  for i in range(7):
    posnLabel = posnLabels[i]  
    for j in range(9):
      atomLabel1 = atomLabels[j/3]
      atomLabel2 = atomLabels[j%3]

      distr = normaliseDistribution(data[i][j])
      data2 = ' '.join([str(x) for x in distr])

      interFile.write('%s %s %s %s\n' % (posnLabel, atomLabel1, atomLabel2, data2))

  for ccpCode in intra.keys():
    for atom1 in intra[ccpCode].keys():
      for atom2 in intra[ccpCode][atom1].keys():
        distr = intra[ccpCode][atom1][atom2]
        N = 0.0
        for key in distr.keys():
          N += float(distr[key])
        
        for key in distr.keys():
          distr[key] /= N
          
        intra[ccpCode][atom1][atom2] = distr

  intraFile = open(intraFileName, 'w')
  cPickle.dump(intra, intraFile)
  
  intraFile.close()
  interFile.close()
  
def getCoordsDist(coord1, coord2):

  dx = coord1[0] - coord2[0]
  dy = coord1[1] - coord2[1]
  dz = coord1[2] - coord2[2]
  dist = sqrt((dx*dx) + (dy*dy) + (dz*dz))

  return dist

def readIntraDistribution(fileName):

  
  file = open(fileName, 'r')
  dict = cPickle.load(file)

  return dict


def readInterDistribution(fileName):

  data = []
  for i in range(7):
    data.append([])
    for j in range(9):
      data[i].append(None)
  

  file = open(fileName)

  line = file.readline()

  i = 0
  while line:
    array = line.split()
    distr = array[3:]
    j = i/9
    k = i%9
    dict = {}
    
    for l in range(200):
      m = str(l/10.0)
      dict[m] = float(distr[l])

    data[j][k] = dict
    line = file.readline()
    i += 1

  return data


def normaliseDistribution(dict):


  distribution = []

  N = 0
  for i in range(200):
    k = str(i/10.0)
    N += dict.get(k, 0)
    
  for i in range(200):
    k = str(i/10.0)
    
    if N > 0:
      distribution.append(dict.get(k, 0.0)/float(N))
    else:
      distribution.append(0.0)

  return distribution


def getAlphaCoord(dict):

  coord = None

  if dict.get('HA'):
    coord = dict['HA']
  else:
    coords = []
    for atom in dict.keys():
      if atom[:2] == 'HA':
        coords.append(dict[atom])

    if coords:
      coord = averageCoord(coords)
  
  return coord


def getAmideCoord(dict):

  coord = None

  if dict.get('H'):
    coord = dict['H']
  elif dict.get('HN'):
    coord = dict['HN']
  else:
    coords = []
    for atom in dict.keys():
      if atom[:2] == 'HD':
        coords.append(dict[atom])

    if coords:
      coord = averageCoord(coords)
  
  return coord


def getBetasCoord(dict):    
  
  coord = None

  if dict.get('HB'):
    coord = dict['HB']
  else:
    coords = []
    for atom in dict.keys():
      if atom[:2] == 'HB':
        coords.append(dict[atom])

    if coords:
      coord = averageCoord(coords)

  return coord


def averageCoord(coordList):

  sum = [0.0, 0.0, 0.0]
  for coord in coordList:
    for i in range(3):
      sum[i] += coord[i]

  N = float(len(coordList))
  for i in range(3):
    sum[i] /= N

  return sum


def getEnsembleCoordsDist(ensemble1, ensemble2):

  n    = len(ensemble1) 
  dist = 0.0

  for j in range(n):
    coords1 = ensemble1[j]
    coords2 = ensemble2[j]
    sum = 0.0
    for i in range(3):
      diff = coords2[i]-coords1[i]
      sum += diff * diff
    
    dist += sqrt(sum)
  dist /= float(n)
    
  return sqrt(sum)


def getResonanceCoords(clouds, spinSystem, atomType):

  # find spin system in cloud
  # find atom type within spin system
  # get coords for each cloud of an FOC

  if not hasattr(spinSystem, 'resonanceCoords'):
    spinSystem.resonanceCoords = {}
    
  if spinSystem.resonanceCoords.get(atomType):
    return spinSystem.resonanceCoords[atomType]

  resonances = []
  for resonance in spinSystem.resonances:
    if resonance.name == atomType:
      resonances.append(resonance)
    elif atomType in resonance.assignNames:
      resonances.append(resonance)
    elif resonance.name and resonance.name[:2] == atomType:
      resonances.append(resonance)

  ensemble = []
  for cloud in clouds:
    coords = []
    for resonance in resonances:
      coord = getResonanceCoord(cloud, resonance)
      if coord:
        coords.append(coord)

    if coords:
      coord = averageCoord(coords)
      ensemble.append(coord)

  spinSystem.resonanceCoords[atomType] = ensemble

  return ensemble
  
def getResonanceCoord(cloud, resonance):

  # not all spinSystem resonances are represented in a cloud
  return cloud.get(resonance)

if __name__ == '__main__':
  calcDistDistributions('/home/tjs234/ccpn/recoord/')
  
  

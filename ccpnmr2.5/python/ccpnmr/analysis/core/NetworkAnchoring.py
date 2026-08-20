
"""
======================COPYRIGHT/LICENSE START==========================

NetworkAnchoring.py: Part of the CcpNmr Analysis program

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

from math import log, sqrt

from ccpnmr.analysis.core.ExperimentBasic import getDataDimIsotopes, getOnebondDataDims
from ccpnmr.analysis.core.ExperimentBasic import findSpectrumDimsByIsotope, getThroughSpaceDataDims
from ccpnmr.analysis.core.MoleculeBasic   import getLinkedResidue, areResonancesBound
from ccpnmr.analysis.core.MoleculeBasic   import DEFAULT_ISOTOPES, getNumConnectingBonds
from ccpnmr.analysis.core.AssignmentBasic import findConnectedSpinSystem, makeResonanceGuiName
from ccpnmr.analysis.core.AssignmentBasic import findMatchingPeakDimShifts, assignResToDim
from ccpnmr.analysis.core.AssignmentBasic import getResonanceLabellingFraction, getResonancePairLabellingFraction
from ccpnmr.analysis.core.ConstraintBasic import getMeanPeakIntensity, getIntensityDistanceTable, getNoeDistance
#from ccpnmr.analysis.core.StructureBasic     import getAtomSetsDistance
#from ccpnmr.analysis.core.PeakBasic          import pickPeak


hTol = 0.04 
nTol = 0.4  
cTol = 0.4  

TOLERANCE_DICT = {'1H':hTol,'15N':nTol,'13C':cTol,'13C,15N':cTol}

def networkAnchorAssign(peakLists, intensityType='height', strictness=2, threshold=1.0,
                        isotopeTolerances=None, assignPeakList=True, constraintSet=None,
                        labelling=None, minLabelFraction=0.1,
                        distParams=None, structure=None, progressBar=None, nexus=None):

  if not peakLists:
    return
 
  if isotopeTolerances:
    TOLERANCE_DICT = isotopeTolerances
  
  distanceFunction = None
  if distParams:  
    distanceFunction = lambda val:getNoeDistance(val, distParams)
     
  project    = peakLists[0].root
  nmrProject = peakLists[0].topObject
  molSystem  = project.findFirstMolSystem()
  shiftList  = peakLists[0].dataSource.experiment.shiftList
  isotopes = set([])

  # Get known network of connectivities inc covalent and NOE
  network = {}
  covalent = {}
  
  # Assign some peaks with uniquely matching shifts
  for peakList in peakLists:
     if labelling is True:
       expLabelling = peakList.dataSource.experiment
     else:
       expLabelling = labelling
     
     network, covalent = getCloseSingleShiftMatches(peakList, network, covalent,
                                                    labelling=expLabelling,
                                                    minLabelFraction=minLabelFraction,
                                                    intensityType=intensityType,
                                                    progressBar=progressBar)
  
  
  # Get existing assigned NOE network
  totalPeaks = 0
  for peakList in peakLists:
    if not peakList.peaks:
      continue
  
    meanIntensity = getMeanPeakIntensity(peakList.peaks, intensityType=intensityType)
      
    spectrum = peakList.dataSource
    distDims = getThroughSpaceDataDims(spectrum)
    for dataDim in distDims:
      isotopes.update(getDataDimIsotopes(dataDim))
    
    hDims = [dd.dim-1 for dd in distDims]
    

    numPeaks = len(peakList.peaks)
    totalPeaks += numPeaks
    info = (spectrum.experiment.name, spectrum.name, peakList.serial, numPeaks)
    
    if progressBar:
      progressBar.setText('Get existing NOE network\nfor %s:%s:%d - %d peaks' % info)
      progressBar.set(0)
      progressBar.total = numPeaks
      progressBar.open()
      progressBar.update_idletasks()
      
    else:
      print 'Get existing NOE network for %s:%s:%d - %d peaks' % info

    if len(hDims) != 2:
      continue
    
    for peak in peakList.peaks:
      if progressBar:
        progressBar.increment()
    
      peakDims  = peak.sortedPeakDims()
      peakDim1  = peakDims[hDims[0]]
      contribs1 = peakDim1.peakDimContribs
      if contribs1:
        peakDim2  = peakDims[hDims[1]]
        contribs2 = peakDim2.peakDimContribs
        if contribs2:
          peakIntensity = peak.findFirstPeakIntensity(intensityType=intensityType)
          if peakIntensity:
            intensity = peakIntensity.value
            intensity /= float(len(contribs1))
            intensity /= float(len(contribs2))
            intensity /= meanIntensity
            
            for contrib1 in contribs1:
              resonance1 = contrib1.resonance
              if network.get(resonance1) is None:
                covalent[resonance1] = {}
                network[resonance1] = {}
                
              intensity2 = intensity
              if resonance1.resonanceSet:
                intensity2 /= float(len(resonance1.resonanceSet.findFirstAtomSet().atoms))
                
              for contrib2 in contribs2:
                resonance2 = contrib2.resonance
                if network.get(resonance2) is None:
                  network[resonance2] = {}
                  covalent[resonance2] = {}
                
                if resonance2.resonanceSet:
                  intensity2 /= float(len(resonance2.resonanceSet.findFirstAtomSet().atoms))
                    
                if not contrib2.peakContribs:
                  if not contrib1.peakContribs:
                    network[resonance1][resonance2] = [intensity2,peak]
                    network[resonance2][resonance1] = [intensity2,peak]
 
                else:
                  for peakContrib in contrib2.peakContribs:
                    if peakContrib in contrib1.peakContribs:
                      network[resonance1][resonance2] = [intensity2,peak]
                      network[resonance2][resonance1] = [intensity2,peak]
                      break
            
          else:
            pass # Should warn
  
  covalentIntensity = 5.0 # Need to optimise this
  
  # Get covalent network
  if progressBar:
    progressBar.setText('Getting covalent network')
    progressBar.set(0)
    progressBar.total = len(nmrProject.resonances)
    progressBar.open()
    progressBar.update_idletasks()
    
  else:  
    print 'Getting covalent network - %d resonances' % len(nmrProject.resonances)

  neighbours = {}
  chemAtomToAtom = {}
  
  c = 0
  for resonance in nmrProject.resonances:
    if progressBar:
      progressBar.increment()
    
    if resonance.isotopeCode not in isotopes:
      continue
      
    resonanceSet = resonance.resonanceSet
    if not resonanceSet:
      continue
      
    if network.get(resonance) is None:
      network[resonance] = {}
      covalent[resonance] = {}
      
    for atomSet in resonanceSet.atomSets:
      for atom in atomSet.atoms:
        residue = atom.residue
        
        if chemAtomToAtom.get(residue) is None:
          chemAtomToAtom[residue] = {}
          for atom2 in residue.atoms:
            chemAtomToAtom[residue][atom2.chemAtom] = atom2
        
        chemAtom = atom.chemAtom
        if chemAtom.waterExchangeable:
          continue
          
        chemAtoms = neighbours.get(chemAtom)
        
        if chemAtoms is None:
          chemAtoms = []
          for atom2 in residue.atoms:
            if atom2 is atom:
              continue
            
            chemAtom2 = atom2.chemAtom
            if DEFAULT_ISOTOPES.get(chemAtom2.elementSymbol) not in isotopes:
              continue
            
            numBonds = getNumConnectingBonds(atom, atom2, limit=6)
            if numBonds < 5:
              chemAtoms.append(chemAtom2)
                    
          neighbours[chemAtom] = chemAtoms
        
        atoms = []
        for chemAtomB in chemAtoms:
          atom2 = chemAtomToAtom[residue].get(chemAtomB)
          if atom2 is not None:
            atoms.append(atom2)
        
        residue2 = getLinkedResidue(residue, 'prev')
        if residue2:
          for atom2 in residue2.atoms:
            chemAtom2 = atom2.chemAtom
            if DEFAULT_ISOTOPES.get(chemAtom2.elementSymbol) not in isotopes:
              continue
            
            numBonds = getNumConnectingBonds(atom, atom2, limit=6)
            if numBonds < 5:
              atoms.append(atom2)
               
        residue2 = getLinkedResidue(residue, 'next')
        if residue2:
          for atom2 in residue2.atoms:
            chemAtom2 = atom2.chemAtom
            if DEFAULT_ISOTOPES.get(chemAtom2.elementSymbol) not in isotopes:
              continue
            
            numBonds = getNumConnectingBonds(atom, atom2, limit=6)
            if numBonds < 5:
              atoms.append(atom2)
         
        for atom2 in atoms:
          atomSet2 = atom2.atomSet
          if atomSet2 and (atomSet2 is not atomSet):
            for resonanceSet2 in atomSet2.resonanceSets:
              for resonance2 in resonanceSet2.resonances:
                if network.get(resonance2) is None:
                  network[resonance2] = {}
                  covalent[resonance2] = {}
                
                if network[resonance].get(resonance2) is None: # Not already in network
                  network[resonance][resonance2] = [covalentIntensity, None]
                  network[resonance2][resonance] = [covalentIntensity, None]
                  covalent[resonance][resonance2] = True
                  covalent[resonance2][resonance] = True
                  c += 1
                 
  #print 'Atom pair network connections %d' % c
  
  c = 0
  for ss in nmrProject.resonanceGroups:
    ss2 = findConnectedSpinSystem(ss, delta=-1)
    if ss2:
      for r1 in ss.resonances:
        if r1.isotopeCode not in isotopes:
          continue
          
        for r2 in ss.resonances:
          if r2.isotopeCode not in isotopes:
            continue
            
          if network.get(r1) is None:
            network[r1]  = {}
            covalent[r1] = {}
          if network.get(r2) is None:
            network[r2]  = {}
            covalent[r2] = {}
          
          if network[r1].get(r2) is None:
            network[r1][r2]  = [covalentIntensity, None]
            network[r2][r1]  = [covalentIntensity, None]
            covalent[r1][r2] = True
            covalent[r2][r1] = True
            c += 1

  #print 'Anonymous intra residue connections %d' % c
  
  done = {}
  iter = 0
  nAssign = 1
  k = 0
  dataSets = []
  while nAssign > 0:
  #while iter < 1:
    data = []

    nAssign = 0
    iter += 1
    
    if progressBar:
      progressBar.setText('Anchoring iteration %d' % iter)
      progressBar.set(0)
      progressBar.total = totalPeaks
      progressBar.open()
      progressBar.update_idletasks()
    else:
      print 'Anchoring iteration %d' % iter
    
    closeResonancesDict = {}    
    
    for peakList in peakLists:
      if not peakList.peaks:
        continue
        
      spectrum = peakList.dataSource
      distDims = getThroughSpaceDataDims(spectrum)
      hDims = [dd.dim-1 for dd in distDims]
      tolerances = getDimTolerances(spectrum)
      bondedDims = getBondedDimsDict(spectrum)
      meanIntensity = getMeanPeakIntensity(peakList.peaks, intensityType=intensityType)

      info = (spectrum.experiment.name, spectrum.name, peakList.serial, len(peakList.peaks))
      #print '  Using %s:%s:%d - %d peaks' % info
      if len(hDims) != 2:
        continue
 
      for peak in peakList.peaks:
        if progressBar:
          progressBar.increment()
        
        if done.get(peak):
          continue
        
        peakDims = peak.sortedPeakDims()
 
        if peakDims[hDims[0]].peakDimContribs:
          if peakDims[hDims[1]].peakDimContribs:
            continue
 
        boundResonances = {}
        
        if closeResonancesDict.get(peak) is None:
          possibles = []
          for dim in hDims:
            peakDim = peakDims[dim]
            if peakDim.peakDimContribs:
              possibles.append([contrib.resonance for contrib in peakDim.peakDimContribs])
              continue
 
            resonances = []
            shifts = findMatchingPeakDimShifts(peakDim, shiftRanges=None,
                                               tolerance=tolerances[dim], aliasing=True,
                                               findAssigned=False)
            dim2 = bondedDims.get(dim)
            peakDim2 = None
            if dim2 is not None:
              peakDim2 = peakDims[dim2]
              shifts2  = findMatchingPeakDimShifts(peakDim2, shiftRanges=None,
                                                  tolerance=tolerances[dim2],
                                                  aliasing=True,findAssigned=False)

              for shift in shifts:
                resonance1 = shift.resonance
                for shift2 in shifts2:
                  resonance2 = shift2.resonance
                  if areResonancesBound(resonance1, resonance2):
                    if labelling:
                      fraction = getResonancePairLabellingFraction(resonance1,
                                                                   resonance2,
                                                                   expLabelling)
                      if fraction < minLabelFraction:
                        continue
 
 
                    resonances.append(resonance1)
                    boundResonances[resonance1] = resonance2
                    break
 
            else:
              for shift in shifts:
                resonance = shift.resonance
 
                if labelling:
                  fraction = getResonanceLabellingFraction(resonance, expLabelling)
 
                  if fraction < minLabelFraction:
                        continue
 
                resonances.append(resonance)
 
            possibles.append(resonances)
          closeResonancesDict[peak] = possibles
          
        else:
          possibles = closeResonancesDict[peak]

        if not possibles[0]:
          continue # warn
        
        if not possibles[1]:
          continue # warn
                
        peakIntensity = peak.findFirstPeakIntensity(intensityType=intensityType)
        if not peakIntensity:
          print 'Peak missing intensity', peak
          continue
       
        else:  
          intensity = peakIntensity.value/meanIntensity
          
        scores = {}
        numbers = {}
 
        bestScore = None
        bestPair  = None
        for resonance1 in possibles[0]:
          if not resonance1.resonanceGroup:
            continue
        
          if network.get(resonance1) is None:
            continue
       
          anchors1 = network[resonance1].keys()
          spinSystem1 = resonance1.resonanceGroup
          
          for resonance2 in possibles[1]:
            if not resonance2.resonanceGroup:
              continue
              
            if network.get(resonance2) is None:
              continue
 
            anchors2 = network[resonance2].keys()
            pair = (resonance1, resonance2)
            spinSystem2 = resonance2.resonanceGroup
                          
            for resonance3 in anchors1:
 
              spinSystem3 = resonance3.resonanceGroup
              if (strictness > 0) and spinSystem3:
                check = False
                if (spinSystem3 is spinSystem1) or (spinSystem3 is spinSystem2):
                  check = True
 
                if (strictness>1) and spinSystem3.residue:
                  if spinSystem1.residue:
                    if abs(spinSystem1.residue.seqCode - spinSystem3.residue.seqCode) < 2:
                      check = True
 
                  if spinSystem2.residue:
                    if abs(spinSystem2.residue.seqCode - spinSystem3.residue.seqCode) < 2:
                      check = True
 
  
              else:
                check = True
                
 
              if check and (resonance3 in anchors2):
                intensityA, peakA = network[resonance1][resonance3]
                intensityB, peakB = network[resonance2][resonance3]
 
                if scores.get(pair) is None:
                  scores[pair] = 0.0
                  numbers[pair] = 0.0

                shift1 = resonance1.findFirstShift(parentList=shiftList)
                shift2 = resonance2.findFirstShift(parentList=shiftList)

                if shift1 and shift2:

                  delta1 = abs(peakDims[hDims[0]].realValue - shift1.value)
                  delta2 = abs(peakDims[hDims[1]].realValue - shift2.value)
                  
                  tol1 = TOLERANCE_DICT[resonance1.isotopeCode]
                  tol2 = TOLERANCE_DICT[resonance2.isotopeCode]

                  match1 = (tol1-delta1)/tol1
                  match2 = (tol2-delta2)/tol2

                  scores[pair] += intensityA*intensityB*match1*match2
                  numbers[pair] += 1
                
        nGood = 0
        for pair in scores.keys():
          resonance1, resonance2 = pair
                  
          score = scores[pair]
 
          if score > threshold:
            nGood += 1
 
            
          if (bestScore is None) or (score > bestScore):
            bestScore = score
            bestPair  = pair
 
        #if bestScore and (nGood < 2):
        for pair in scores.keys():
          resonance1, resonance2 = pair
                      
          if scores[pair] < threshold:
            #if len(scores.keys()) > 1:
            continue
          else:  
            intensity2 = intensity/nGood
            
          bestPair = pair
          bestScore = scores[pair] 
        
          for i in (0,1):
            resonance = bestPair[i]

            if assignPeakList:
              assignResToDim(peakDims[hDims[i]], resonance, doWarning=False)

              bound = boundResonances.get(resonance)
              if bound:
                dim2 = bondedDims.get(hDims[i])
                if dim2 is not None:
                  assignResToDim(peakDims[dim2], bound, doWarning=False)

            if network.get(resonance) is None:
              network[resonance] = {}

          #name1 = makeResonanceGuiName(bestPair[0])
          #name2 = makeResonanceGuiName(bestPair[1])          
          
          if resonance1.resonanceSet:
            intensity2 /= float(len(resonance1.resonanceSet.findFirstAtomSet().atoms))
          if resonance2.resonanceSet:
            intensity2 /= float(len(resonance2.resonanceSet.findFirstAtomSet().atoms))
          
          if labelling:
            intensity2 /= getResonanceLabellingFraction(resonance1, expLabelling)
            intensity2 /= getResonanceLabellingFraction(resonance2, expLabelling)
           
          nAssign += 1
          
          covalent[bestPair[0]][bestPair[1]] = None
          covalent[bestPair[1]][bestPair[0]] = None
          network[bestPair[0]][bestPair[1]]  = [intensity2, peak]
          network[bestPair[1]][bestPair[0]]  = [intensity2, peak]
          done[peak]= True
     
    #print '  Assigned:', nAssign
    dataSets.append(data)
  

  from ccpnmr.analysis.core.ConstraintBasic import getDistancesFromIntensity, getIntensityDistanceTable
  from ccpnmr.analysis.core.ConstraintBasic import makeNmrConstraintStore, getFixedResonance
  
  if constraintSet is None:
    constraintSet = makeNmrConstraintStore(nmrProject)
  
  if not constraintSet:
    return
  
  constraintList = constraintSet.newDistanceConstraintList()
  
  peakConstraints = {}
  doneResonances  = {}
  for resonance1 in network.keys():
    fixedResonance1 = getFixedResonance(constraintSet,resonance1)
    
    for resonance2 in network[resonance1].keys():
      if resonance1 is resonance2: 
        continue
    
      key = [resonance1.serial, resonance2.serial]
      key.sort()
      key = tuple(key)
      
      if doneResonances.get(key):
        continue
      else:
        doneResonances[key] = True 
    
      if covalent.get(resonance1):
        if covalent[resonance1].get(resonance2):
          # J connected are close so what do we do...? 
          #print "Skip", makeResonanceGuiName(resonance1), makeResonanceGuiName(resonance2)
          continue
    
      fixedResonance2 = getFixedResonance(constraintSet,resonance2)
      intensity, peak = network[resonance1][resonance2]
      
      if peak:
        peakList        = peak.peakList
        spectrum        = peakList.dataSource
        experiment      = spectrum.experiment
 
        if not distanceFunction:
          noeDistClasses   = getIntensityDistanceTable(spectrum)
          distanceFunction = lambda val:getDistancesFromIntensity(noeDistClasses,val)
        
        constraint = peakConstraints.get(peak)
        if not constraint:
          constraint = constraintList.newDistanceConstraint(weight=1.0, origData=intensity)
          peakContrib = constraint.newConstraintPeakContrib(experimentSerial=experiment.serial,
                                                            dataSourceSerial=spectrum.serial,
                                                            peakListSerial=peakList.serial,
                                                            peakSerial=peak.serial)
          peakConstraints[peak] = constraint
        else:
          intensity += constraint.origData 
 
        dist, minDist, maxDist = distanceFunction(intensity)
        error = abs(maxDist - minDist)
        
        constraint.origData    = intensity 
        constraint.targetValue = dist
        constraint.upperLimit  = maxDist
        constraint.lowerLimit  = minDist
        constraint.error       = error

 
        item = constraint.newDistanceConstraintItem(resonances=[fixedResonance1,fixedResonance2])

  return constraintList

def getBondedDimsDict(spectrum):

  bondedDims = {}
  for dataDim1, dataDim2 in getOnebondDataDims(spectrum):
    i = dataDim1.dim - 1
    j = dataDim2.dim - 1
    bondedDims[i] = j
    bondedDims[j] = i
  
  return bondedDims

def getDimTolerances(spectrum):

  tolerances = []
  for dataDim in spectrum.sortedDataDims():
    ll = [(x.expDimRef.serial, x.expDimRef) for x in dataDim.dataDimRefs]
    ll.sort()
    for junk,expDimRef in ll:
      if expDimRef.measurementType in ('Shift','shift'):
        isotopes = list(expDimRef.isotopeCodes)
        isotopes.sort()
        key = ','.join(isotopes)
        tolerances.append(TOLERANCE_DICT.get(key, 0.2))
        break
  
  return tolerances

def getNexusDist(nexus, resonanceA, resonanceB):

  ensemble1 = getNexusCoord(nexus, resonanceA)
  ensemble2 = getNexusCoord(nexus, resonanceB)

  dist = None

  if ensemble1 and ensemble2:
    dist = getEnsembleCoordsDist(ensemble1, ensemble2)

  return dist

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
    
  return dist


def getNexusCoord(nexus, resonance):

  ensemble = []
  for cloud in nexus:
    coords = []
    for resonance in [resonance,]:
      coord = getResonanceCoord(cloud, resonance)
      if coord:
        coords.append(coord)

    if coords:
      coord = averageCoord(coords)
      ensemble.append(coord)

  return ensemble

def getResonanceCoord(cloud, resonance):

  return cloud.get(resonance)
  
def averageCoord(coordList):

  sum = [0.0, 0.0, 0.0]
  for coord in coordList:
    for i in range(3):
      sum[i] += coord[i]

  N = float(len(coordList))
  for i in range(3):
    sum[i] /= N

  return sum
          
def assignCloseSingleShiftMatches(peakList):


  spectrum   = peakList.dataSource
  experiment = spectrum.experiment
  shiftList  = experiment.shiftList
  hDims = findSpectrumDimsByIsotope(spectrum, '1H')
  tolerances = getDimTolerances(spectrum)
  bondedDims = getBondedDimsDict(spectrum)

  info = (experiment.name, spectrum.name, peakList.serial, len(peakList.peaks))
  print 'Assigning obvious shift matches for %s:%s:%d - %d peaks' %  info
  
  c = 0
  for peak in peakList.peaks:
    if c and (c%100 == 0):
      print '   %d' % (c,)
    c += 1
    
    peakDims = peak.sortedPeakDims()
    boundResonances = {}
    for dim in hDims:
      peakDim = peakDims[dim]
      if not peakDim.peakDimContribs:
        resonances = []
        shifts = findMatchingPeakDimShifts(peakDim, shiftRanges=None,
                                           tolerance=tolerances[dim], aliasing=True,
                                           findAssigned=False)
        dim2 = bondedDims.get(dim)
        peakDim2 = None
        if dim2 is not None:
          peakDim2 = peakDims[dim2]
          shifts2  = []

          if peakDim2.peakDimContribs:
            for contrib in peakDim2.peakDimContribs:
              shift = contrib.resonance.findFirstShift(parentList=shiftList)
              if shift:
                shifts2.append(shift)

          else:
            shifts2 = findMatchingPeakDimShifts(peakDim2, shiftRanges=None,
                                                tolerance=tolerances[dim2],
                                                aliasing=True,findAssigned=False)
          for shift in shifts:
            resonance1 = shift.resonance
            for shift2 in shifts2:
              resonance2 = shift2.resonance
              if areResonancesBound(resonance1, resonance2):
                resonances.append(resonance1)
                boundResonances[resonance1] = resonance2
                break

        else:
          for shift in shifts:
            resonances.append(shift.resonance)

        if len(resonances) == 1:
          resonance  = resonances[0]
          assignResToDim(peakDim,resonance,doWarning=False)

          resonance2 = boundResonances.get(resonance)
          if peakDim2 and (resonance2 is not None):
             assignResToDim(peakDim2,resonance2,doWarning=False)
            
def getCloseSingleShiftMatches(peakList, network, covalent,
                               labelling=None, minLabelFraction=0.1,
                               intensityType='height', progressBar=None):

  spectrum   = peakList.dataSource
  experiment = spectrum.experiment
  shiftList  = experiment.shiftList
  hDims      = findSpectrumDimsByIsotope(spectrum, '1H')
  tolerances = getDimTolerances(spectrum)
  bondedDims = getBondedDimsDict(spectrum)
  
  meanIntensity = getMeanPeakIntensity(peakList.peaks, intensityType=intensityType)        

  info = (experiment.name, spectrum.name, peakList.serial, len(peakList.peaks))
  
  if progressBar:
    progressBar.setText('Unique shift matches for\n%s:%s:%d - %d peaks' %  info)
    progressBar.set(0)
    progressBar.total = len(peakList.peaks)
    progressBar.open()
    progressBar.update_idletasks()
  else:
    print 'Determining unique shift matches for %s:%s:%d - %d peaks' %  info
  
  for peak in peakList.peaks:

    if progressBar:
      progressBar.increment()
    
    hResonances = []
    
    peakDims = peak.sortedPeakDims()
    boundResonances = {}
    for dim in hDims:
      peakDim = peakDims[dim]
      if peakDim.peakDimContribs:
        hResonances.append([contrib.resonance for contrib in peakDim.peakDimContribs])
      
      else:
        resonances = []
        shifts = findMatchingPeakDimShifts(peakDim, shiftRanges=None,
                                           tolerance=tolerances[dim], aliasing=True,
                                           findAssigned=False)
        dim2 = bondedDims.get(dim)
        peakDim2 = None
        if dim2 is not None:
          peakDim2 = peakDims[dim2]
          shifts2  = []

          if peakDim2.peakDimContribs:
            for contrib in peakDim2.peakDimContribs:
              shift = contrib.resonance.findFirstShift(parentList=shiftList)
              if shift:
                shifts2.append(shift)

          else:
            shifts2 = findMatchingPeakDimShifts(peakDim2, shiftRanges=None,
                                                tolerance=tolerances[dim2],
                                                aliasing=True,findAssigned=False)
          for shift in shifts:
            resonance1 = shift.resonance
            for shift2 in shifts2:
              resonance2 = shift2.resonance
              if areResonancesBound(resonance1, resonance2):
                if labelling:
                  fraction = getResonancePairLabellingFraction(resonance1,
                                                               resonance2,
                                                               labelling)
                  if fraction < minLabelFraction:
                    continue
                    
                resonances.append(resonance1)
                boundResonances[resonance1] = resonance2
                break

        else:
          for shift in shifts:
            resonance = shift.resonance
            
            if labelling:
              fraction = getResonanceLabellingFraction(resonance, labelling)
 
              if fraction < minLabelFraction:
                    continue
 
            resonances.append(resonance)

        if len(resonances) == 1:
          #resonance  = resonances[0]
          #assignResToDim(peakDim,resonance,doWarning=False)

          #resonance2 = boundResonances.get(resonance)
          #if peakDim2 and (resonance2 is not None):
          #   assignResToDim(peakDim2,resonance2,doWarning=False)
    
          hResonances.append(resonances)
      
    if len(hResonances) == 2:        
            
      peakIntensity = peak.findFirstPeakIntensity(intensityType=intensityType)
      if peakIntensity:
        intensity = peakIntensity.value
        intensity /= float(len(hResonances[0]))
        intensity /= float(len(hResonances[1]))
        intensity /= meanIntensity
 
        for resonance1 in hResonances[0]:
          if network.get(resonance1) is None:
            covalent[resonance1] = {}
            network[resonance1] = {}
 
          intensity2 = intensity
          if resonance1.resonanceSet:
            intensity2 /= float(len(resonance1.resonanceSet.findFirstAtomSet().atoms))
 
          for resonance2 in hResonances[1]:
            if network.get(resonance2) is None:
              network[resonance2] = {}
              covalent[resonance2] = {}
 
            if resonance2.resonanceSet:
              intensity2 /= float(len(resonance2.resonanceSet.findFirstAtomSet().atoms))
             
            network[resonance1][resonance2] = [intensity2,peak]
            network[resonance2][resonance1] = [intensity2,peak]
 
  return network, covalent
     
      
      
      
      
      
      
      
      
      
      
      
      
            
            
            
            
            
            
            
            
            

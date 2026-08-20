
"""
======================COPYRIGHT/LICENSE START==========================

NexusBasic.py: Part of the CcpNmr Nexus program

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
from ccpnmr.analysis.core.AssignmentBasic import findConnectedSpinSystem, getBoundResonances
from ccpnmr.analysis.core.AssignmentBasic import addSpinSystemResonance, assignResToDim
from ccpnmr.analysis.core.AssignmentBasic import clearPeakDim, assignResonanceType
from ccpnmr.analysis.core.AssignmentBasic import makeSeqSpinSystemLink

from ccpnmr.analysis.core.Util import getAnalysisDataDim
from ccpnmr.analysis.core.MoleculeBasic import getLinkedResidue
from ccpnmr.analysis.core.PeakBasic import setupPeakHeight
from ccpnmr.analysis.core.ExperimentBasic import getOnebondDataDims, getPrimaryDataDimRef
from ccpnmr.analysis.core.ExperimentBasic import getSeqAssignRefExperiments

CARBONYL_TYPE = set(('H[N[CO]]', 'H[N[ca[CO]]]',
                     'CONH','HNCO','HNcaCO',
                     'H[N[coca][ca[CO]]','H[N_[ca[CO]].onebond)'))
                     
CA_TYPE = set(('H[N[co[CA]]]', 'H[N[CA]]', 'hCcoNH','H[N[CA]].seq',
               'hCAcoNH','CANH','HNCA','HNcoCA','haCANH',
               'H[N[coca][CA]]','H[N_[CA]].onebond' ))

CB_TYPE = set(('H[N[co[c[C]]]]','H[N[ca[CB]]]','H[N[co[ca[CB]]]]'))

CACB_SIGN_TYPE = set(('H[N[co[{CA|ca[C]}]]]',
                      'H[N[{CA|ca[Cali]}]]',
                      'H[N[coca][{CA|ca[Cali]}]',
                      'H[N[{CA|ca[Cali]}]].seq',
                      'H[N_[{CA|ca[Cali]}]].onebond'))

CACB_TYPE = set(('h{CA|Cca}coNH',
                 'h{CA|Cca}NH',
                 '{CA|Cca}coNH',
                 '{CA|Cca}NH',
                 'h{CA|Cca}NH.seq'))

HA_TYPE = set(('HA[ca[N[H]]]','H[N[ca[HA]]]',
               'HAcaNH','HNcaHA',
               'HNcocaH','HcacoNH','H{[N]+[HA]}'))

HB_TYPE = set(('H[N[HB]]',))

HAHB_SIGN_TYPE = set(('H[N[co[{ca[H]|ca[c[H]]}]]]','H[N[co[{ca[H]|cac[H]}]]]',
                      'H[N[{ca[H]|ca[c[H]}]]','H[N[{ca|cac}[H]]]'))

HAHB_TYPE = set(('H{ca|cca}NH','H{ca|cca}coNH',
                 'H{ca|cca}NH.seq','HN{ca|cac}H.seq',
                 'HNco{ca|cac}H','HN{ca|cac}H','H{ca|cca}NH'))

def guessAtomType(refExperiment, ppm, height, v=False):

  expType = refExperiment.name
  atomType = None

  if expType in CARBONYL_TYPE:
    atomType = 'C'
 
  elif expType in CA_TYPE:
    atomType = 'CA'
 
  elif expType in CACB_SIGN_TYPE:
    if height < 0:
      atomType = 'CB'
    else:
      if ppm < 42.0:
        atomType = 'CB'
      elif ppm < 59.0:
        atomType = 'CA'
      elif ppm > 70.0:
        atomType = 'CB'
      else:
        atomType = 'CA'

  elif expType in CB_TYPE:
    atomType = 'CB'
 
  elif expType in CACB_TYPE:
    if ppm < 42.0:
      atomType = 'CB'
    elif ppm < 59.0:
      atomType = 'CA'
    elif ppm > 70.0:
      atomType = 'CB'
      
  elif expType in HA_TYPE:
    atomType = 'HA'

  elif expType in HB_TYPE:
    atomType = 'HB'

  elif atomType in HAHB_TYPE:
    if ppm > 4.6:
      atomType = 'HA'
    elif ppm < 2.5:
      atomType = 'HB' 
  
  elif atomType in HAHB_SIGN_TYPE:
    if height < 0:
      atomType = 'HB'
    else:
      atomType = 'HA'

  return atomType

def isResonanceAmide(resonance):                  

  isotope = resonance.isotopeCode                 
  assignName = '/'.join(resonance.assignNames)    
  bound = list(resonance.covalentlyBound)         
 
  if isotope == '15N':                            
    if assignName and (assignName != 'N'):        
      return False                                
 
    numH = 0                                      
    for resonance2 in bound:                      
      if resonance2.isotopeCode == '1H':          
        numH +=1                                  
 
    if numH == 1:                                 
      return True                                 
    else:                                         
      return False                                
 
  elif isotope == '1H':                           
    for resonance2 in bound:                      
      if resonance2.isotopeCode == '1H':          
        continue                                  
                                                  
      return isResonanceAmide(resonance2)
                                                  
    else:                                         
      return False
 
  else:
    return False

def getAmideSpinSystems(peakLists, minNitrogenPpm=93.0):

  #
  
  spinSystems = set()
  
  for peakList in peakLists:
    shiftList = peakList.dataSource.experiment.shiftList
  
    for peak in peakList.peaks:
      for peakDim in peak.peakDims:
        ppm = peakDim.value
        
        # Check ppm range
        
        if ppm < minNitrogenPpm: # Also removes 1H and most 1C 
          continue
        
        for contrib in peakDim.peakDimContribs:
          resonance = contrib.resonance
          
          # Check isotope
          
          if resonance.isotopeCode != '15N':
            continue
          
          # Check we have an amide
          
          if not isResonanceAmide(resonance):
            continue
          
          # Check not side chain bound
          
          numH = 0
          for bound in getBoundResonances(resonance):
            if bound.isotopeCode == '1H':
              numH +=1
            
          if numH > 1:
            continue  
          
          spinSystem = resonance.resonanceGroup
 
          if spinSystem:
            spinSystems.add(spinSystem)
        
  ss = [(s.serial, s) for s in spinSystems]
  ss.sort()
  
  return [x[1] for x in ss]
  

def getSpinSystemInterIntraResonances(spinSystem, peakLists):

  nmrProject = spinSystem.nmrProject
  shiftLists = []
  
  for peakList in peakLists:
    spectrum = peakList.dataSource
    experiment = spectrum.experiment
    shiftLists.append(experiment.shiftList)
    
  interResonances = []
  intraResonances = []
  interPpms = []
  intraPpms = []
  interTypes = []
  intraTypes = []
  
  for resonance in spinSystem.resonances:
    if isResonanceAmide(resonance):
      continue
  
    if resonance.shifts:
      for contrib in resonance.peakDimContribs:
        peakList2 = contrib.peakDim.peak.peakList
        
        if peakList2 in peakLists:
          assignName = '/'.join(resonance.assignNames)
          intraResonances.append(resonance)
          intraTypes.append(assignName)
          break
 

  prevSpinSystem = findConnectedSpinSystem(spinSystem)
  if not prevSpinSystem:
    if spinSystem.residue:
      residueB = getLinkedResidue(spinSystem.residue, linkCode='prev')
      prevSpinSystem = nmrProject.findFirstResonanceGroup(residue=residueB)
  
  if prevSpinSystem:
    for resonance in prevSpinSystem.resonances:
      if isResonanceAmide(resonance):
        continue
        
      if resonance.shifts:
        for contrib in resonance.peakDimContribs:
          peakList2 = contrib.peakDim.peak.peakList
        
          if peakList2 in peakLists:
            assignName = '/'.join(resonance.assignNames)
            interResonances.append(resonance)
            interTypes.append(assignName)
            break          

  inter = [interResonances,interPpms]
  intra = [intraResonances,intraPpms]
  
  for resonances, ppms in (inter, intra):
    for resonance in resonances:
      ppmSum = {}
      counts = {}
 
      for shiftList in shiftLists:
        ppmSum[shiftList] = 0.0
        counts[shiftList] = 0.0
 
      for shiftList in shiftLists:
        shift = resonance.findFirstShift(parentList=shiftList)
 
        if shift:
          ppmSum[shiftList] += shift.value
          counts[shiftList] += 1.0
 
      ppm = None
      if shiftLists:
        ppm = 0.0
        for shiftList in shiftLists:
          n = counts[shiftList]
          
          if n:
            ppm += ppmSum[shiftList]/n
 
        ppm /= len(shiftLists)
 
      ppms.append(ppm)
 

  return interResonances, intraResonances, interPpms, intraPpms, interTypes, intraTypes

def linkSpinSystemInterIntraResonances(spinSystem, activeLists, tolerances=None):
  
  nmrProject = spinSystem.topObject
  prevSpinSystem = findConnectedSpinSystem(spinSystem)
  allowedRefExps, coRefExps = getSeqAssignRefExperiments(nmrProject.root)
  
  if not prevSpinSystem:
    if spinSystem.residue:
      residueB = getLinkedResidue(spinSystem.residue, linkCode='prev')
      prevSpinSystem = nmrProject.findFirstResonanceGroup(residue=residueB)
  
  if not prevSpinSystem:
    prevSpinSystem = nmrProject.newResonanceGroup()
    makeSeqSpinSystemLink(prevSpinSystem, spinSystem)
  
  peaks = []
  found = {}
  expTypes = {}
  linkDims = {}    
    
  # go though the peaks associated with this spin system
  # within the selected peak lists
  # find the indirect, linking dimension
  for resonance in spinSystem.resonances:
    if resonance.isotopeCode != '15N':
      continue
  
    if not isResonanceAmide(resonance):
      continue
  
    for contrib in resonance.peakDimContribs:
      peakDim = contrib.peakDim
      peak = peakDim.peak
      
      if found.get(peak):
        continue
      
      peakList = peak.peakList
      if peakList not in activeLists:
        continue
      
      setupPeakHeight(peak)  
      intensity = peak.findFirstPeakIntensity(intensityType='height')
      if not intensity:
        continue
  
      spectrum = peakList.dataSource
      expType = expTypes.get(peakList, spectrum.experiment.refExperiment)
      expTypes[peakList] = expType
  
      linkDim = linkDims.get(peakList)
      if linkDim is None:
        boundDict = {}
        for dataDimA, dataDimB in getOnebondDataDims(spectrum):
          boundDict[dataDimA] = dataDimB
          boundDict[dataDimB] = dataDimA
 
        dims = []
        for dataDim in spectrum.dataDims:
          dataDimRef = getPrimaryDataDimRef(dataDim)
        
          if not dataDimRef:
            continue
            
          isotopes = '/'.join(dataDimRef.expDimRef.isotopeCodes)
          dims.append((isotopes, dataDim))
        
        dims.sort()
        for isotopes, dataDim in dims:
          if '13C' == isotopes:
            linkDim = (dataDim.dim, isotopes)
            break
          
          elif '1H' == isotopes:
            dataDimB = boundDict.get(dataDim)
            
            if dataDimB:
              dataDimRefB = getPrimaryDataDimRef(dataDimB)
              
              if '15N' in dataDimRefB.expDimRef.isotopeCodes:
                continue
          
            linkDim = (dataDim.dim, isotopes)
            break
            
        linkDims[peakList] = linkDim
        
      if linkDim is None:
        continue     
      
      if peakDim.dim is linkDim[0]:
        continue
      
      found[peak] = True
      peakDimL = peak.findFirstPeakDim(dim=linkDim[0])
      isotope = linkDim[1]
      
      peaks.append((peak, expType, peakDimL, isotope, intensity.value))
               
  peakTypes = []
  interResonances = {}    
  intraResonances = {}    
  for peak, expType, peakDim, isotope, height in peaks:
  
    ppm = peakDim.value
    atomType = None
    
    if expType in coRefExps:
      isInter = True
    elif 'N[coca]' in expType.name:
      isInter = False 
    elif 'N_' in expType.name:
      isInter = False 
    else:
      isInter = None  
    
    resonance = None
    contrib = peakDim.findFirstPeakDimContrib()
    if contrib:
      resonance = contrib.resonance
      
      if isInter and resonance.assignNames:
        atomType = resonance.assignNames[0]
        interResonances[atomType] = (resonance, ppm, prevSpinSystem)
      
      if (isInter is False) and resonance.assignNames:
        atomType = resonance.assignNames[0]
        intraResonances[atomType] = (resonance, ppm, spinSystem)
          
    if not atomType:
      atomType = guessAtomType(expType, ppm, height)
                        
    peakTypes.append((peakDim, atomType, contrib, isotope, isInter, ppm))        
  
  # Get known inter/intra resonance assignments
  # for each atom type
  for peakDim, atomType, contrib, isotope, isInter, ppm in peakTypes:
    if isInter is None:
      continue
    elif isInter:
      resDict = interResonances
      uniqSpinSystem = prevSpinSystem
    else:
      resDict = intraResonances 
      uniqSpinSystem = spinSystem
    
    if atomType:
      resonance, ppm2, ss = resDict.get(atomType, (None, None, None))
      
      if (resonance is None) and contrib:
        resDict[atomType] = (contrib.resonance, ppm, uniqSpinSystem)
  
  # Make any new assignments for the unambig peaks
  untypedIntra = []
  untypedInter = []
  for i, data in enumerate(peakTypes):
    (peakDim, atomType, contrib, isotope, isInter, ppm) = data

    if isInter is None:
      continue
      
    elif isInter:
      untyped = untypedInter
      resDict = interResonances
      uniqSpinSystem = prevSpinSystem
      
    else:
      resDict = intraResonances 
      untyped = untypedIntra
      uniqSpinSystem = spinSystem

    if atomType:
      # Get any previously used resonances for this atom type
      resonance, ppm2, ss = resDict.get(atomType, (None, None, None))
      
      # If no prev resonance, look at the peak assignment
      if (not resonance) and contrib:
        resonance = contrib.resonance
      
      # Check to ensure that any existing resonance
      # is from the correct, this/prev spin system
      if resonance:
        spinSystem2 = resonance.resonanceGroup
        if spinSystem2 and (spinSystem2 is not uniqSpinSystem):
          resonance = None
         
        
      # If no valid reasonance yet, check named ones
      # in the required spin system
      if not resonance:
        for resonance2 in uniqSpinSystem.resonances:
          if atomType in resonance2.assignNames:
            resonance = resonance2
            break

      # If no valid resonance yet, make a new one
      if not resonance:
        resonance = nmrProject.newResonance(isotopeCode=isotope)
      
      # If peak dim is assigned to the wrong resonance
      # clear the assignment
      if contrib and contrib.resonance is not resonance:
        clearPeakDim(peakDim)
        contrib = None
      
      # Ensure the peak dim is assigned correctly  
      if not contrib:
        assignResToDim(peakDim, resonance)
      
      # Check type of resonance set correctly              
      if not resonance.assignNames:
        assignResonanceType(resonance, assignNames=(atomType,))
      
      # Check spin system set correctly
      if resonance.resonanceGroup is not uniqSpinSystem:
        addSpinSystemResonance(uniqSpinSystem, resonance)
              
      # Store resonance by atom type, so that it can
      # be picked up later
      resDict[atomType] = (resonance, ppm, uniqSpinSystem)
        
    else:
      untyped.append((peakDim, contrib, isotope, isInter, ppm, i)) 
    
  uniqResonances = interResonances.values() +  intraResonances.values()
  
  # Cluster untyped peaks
  interCluster = {}
  intraCluster = ()
  
  interData = (untypedInter, interCluster, prevSpinSystem)
  intraData = (untypedIntra, intraCluster, spinSystem)
   
  for untyped, clusters, uniqSpinSystem in (interData, intraData):
 
    for peakDim, contrib, isotope, isInter, ppm, i in untyped:
 
      if peakDim in clusters:
        continue
 
      if tolerances:
        tolerance = tolerances[isotope]
      else:
        tolerance = getAnalysisDataDim(peakDim.dataDim).assignTolerance
 
      if contrib:
        resonance = contrib.resonance
      else:
        resonance = None
 
      cluster = [peakDim]
      interCluster[peakDim] = True
      for peakDimB, contribB, isotopeB, isInterB, ppmB, i in untyped:
        if isotope != isotopeB:
          continue
 
        if abs(ppmB-ppm) > tolerance:
          continue
 
        if peakDim.peak.peakList == peakDimB.peak.peakList:
          continue
 
        if contribB:
          if not resonance:
            resonance = contribB.resonance
          if resonance is not contribB.resonance:
            clearPeakDim(peakDimB)
            assignResToDim(peakDimB, resonance,
                           tolerance=tolerance)
 
        cluster.append(peakDimB)
        clusters[peakDimB] = True
 
      if not resonance:
        # Try to macth to any typed peaks
        for peakDimB, atomType, contribB, isotopeB, isInterB, ppmB in peakTypes:
          if not contribB:
            continue
 
          if not atomType:
            continue
 
          if isotope != isotopeB:
            continue
 
          if abs(ppmB-ppm) > tolerance:
            continue
 
          if peakDim.peak.peakList == peakDimB.peak.peakList:
            continue
 
          resonance = contribB.resonance
 
          # Set type for this peak
          peakTypes[i]= (peakDim, atomType, contrib, isotope, isInter, ppm)
          break
 
 
      if not resonance:
        resonance = nmrProject.newResonance(isotopeCode=isotope)
 
      for peakDimC in cluster:
        clearPeakDim(peakDimC)
        assignResToDim(peakDimC, resonance,
                       tolerance=tolerance)
 
 
      uniqResonances.append((resonance, ppm, uniqSpinSystem))
 
      if resonance.resonanceGroup is not uniqSpinSystem:
        addSpinSystemResonance(uniqSpinSystem, resonance)

  
  # E.g. Find the HNCA peaks which best match the HNcoCA/iHNCA resonances 
  matchDict = {}
  closestDict = {}
  for peakDim, atomType, contrib, isotope, isInter, ppm in peakTypes:
    if isInter is not None:
      # Not through HNcoCA or iHNCA
      continue
 
    if tolerances:
      tolerance = tolerances[isotope]
    else:
      tolerance = getAnalysisDataDim(peakDim.dataDim).assignTolerance
      
    matches = []  
    shiftList = peakDim.peak.peakList.dataSource.experiment.shiftList
    for resonanceB, ppm2, uniqSpinSystem in uniqResonances:

      if resonanceB.isotopeCode != isotope:
        continue
    
      assignNames = resonanceB.assignNames
      if assignNames and (atomType not in assignNames):
        continue
      
      if not ppm2:
        shift = resonanceB.findFirstShift(parentList=shiftList)
        if not shift:
          continue
        ppm2 = shift.value
	
      delta = abs(ppm-ppm2)
      prevDelta = closestDict.get(peakDim)
      if (prevDelta is None) or (delta < prevDelta):
        closestDict[peakDim] = delta
      
      
      if delta > tolerance:
        continue
      
      matches.append((delta, resonanceB, uniqSpinSystem)) 
    
    # Best match has smallest delta
    if matches:
      matches.sort()   
      delta, resonance, uniqSpinSystem = matches[0]
      
      prevDelta, peakDimB, ss = matchDict.get(resonance, (None, None, None))
      if not peakDimB or (delta < prevDelta):
        matchDict[resonance] = (delta, peakDim, uniqSpinSystem)
  
  uniqResonanceDict = {}
  for resonance in matchDict.keys():
    delta, peakDim, uniqSpinSystem = matchDict[resonance]
    uniqResonanceDict[peakDim] = resonance, uniqSpinSystem
  
  # E.g. go through HNCA peaks and assign to
  # HNcoCA or iHNCA resonances if required
  nonMatched  = {} 
  for peakDim, atomType, contrib, isotope, isInter, ppm in peakTypes:
    if isInter is not None:
      continue
    
    if closestDict.get(peakDim) is None:
      # No inter residue peaks at all
      intensity = peak.findFirstPeakIntensity(intensityType='height')
      
      if intensity:
        if not nonMatched.has_key(atomType):
          nonMatched[atomType] = []
        nonMatched[atomType].append((-abs(intensity.value), peakDim, contrib, isotope))

      continue
      
    match = uniqResonanceDict.get(peakDim)
    if match:
      resonance, uniqSpinSystem = match
      
      if tolerances:
        tolerance = tolerances[isotope]
      else:
        tolerance = getAnalysisDataDim(peakDim.dataDim).assignTolerance

      # untyped through-carbonyl resonance can interit type
      # by matching typed peaks
    
      if atomType and not resonance.assignNames:
        resonance.addAssignName(atomType)
    
      # e.g. the HNcoCA resonance is the one to put on this HNCA peak
      if contrib:
        if contrib.resonance is not resonance:
          clearPeakDim(peakDim)
          assignResToDim(peakDim, resonance,
                         tolerance=tolerance)
      else:
        assignResToDim(peakDim, resonance,
                       tolerance=tolerance)
    
    else:
      # No match to an unambig peak
      # store to work out which is best for each atom type
      # e.g. just in case we have two intra residue CA possibilites
      if not nonMatched.has_key(atomType):
        nonMatched[atomType] = []
    
      nonMatched[atomType].append((closestDict[peakDim], peakDim, contrib, isotope))
  
  
  for atomType in nonMatched.keys():
    
    if (atomType in intraResonances) and (atomType in interResonances):
      # Both resonances already found, peak is spurious
      continue
      
    if (atomType in intraResonances) and (atomType not in interResonances):
      # E.g. there was iHNCA but no HNcoCA
      otherSpinSystem = prevSpinSystem
    
    elif (atomType not in intraResonances) and (atomType in interResonances):
      # E.g. there was HNcoCA but no iHNCA
      otherSpinSystem = spinSystem
    
    else:
      # No info
      continue
    
    peakDimList = nonMatched[atomType]
    peakDimList.sort()
    # Assume that if we have two possible ambiguous peaks
    # for any given atom type, then the one furthest from an
    # unambigous peak is the best to take
    deltaPpm, peakDim, contrib, isotope = peakDimList[-1]

    if tolerances:
      tolerance = tolerances[isotope]
    else:
      tolerance = getAnalysisDataDim(peakDim.dataDim).assignTolerance
    
    # E.g. this HNCA peak matches no HNcoCA/iHNCA resonances
    resonance = None
    if contrib:
      resonance = contrib.resonance
    
    if resonance:
      # Already have an assignment
      spinSystem2 = resonance.resonanceGroup
      if spinSystem2 is None:
        addSpinSystemResonance(otherSpinSystem, resonance)
        
      elif spinSystem2 is not otherSpinSystem:
        resonance = nmrProject.newResonance(isotopeCode=isotope)
        addSpinSystemResonance(otherSpinSystem, resonance)
        clearPeakDim(peakDim)
        contrib = assignResToDim(peakDim, resonance,
                                 tolerance=tolerance)
      
    else:
      # Needs a new assignment
      for resonance2 in otherSpinSystem.resonances:
        if atomType in resonance2.assignNames:
          resonance = resonance2
          break
 
      if not resonance:
        resonance = nmrProject.newResonance(isotopeCode=isotope)
        contrib = assignResToDim(peakDim, resonance,
                                 tolerance=tolerance)
        addSpinSystemResonance(otherSpinSystem, resonance)
    
    if not contrib:
      # Resonance has come from the spin system due to its atomTypes
      contrib = assignResToDim(peakDim, resonance,
                               tolerance=2*tolerance)

    if not contrib:
      # Existing resonance on spin system was too far away
      resonance = nmrProject.newResonance(isotopeCode=isotope)
      contrib = assignResToDim(peakDim, resonance,
        		       tolerance=tolerance)
      addSpinSystemResonance(otherSpinSystem, resonance)
                
    if atomType:
      assignNames =(atomType,)
    else:
      assignNames = []
      
    if not resonance.assignNames:
      assignResonanceType(resonance, assignNames=assignNames)
    

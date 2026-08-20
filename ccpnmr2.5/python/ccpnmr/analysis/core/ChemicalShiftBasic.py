
"""
======================COPYRIGHT/LICENSE START==========================

ChemicalShiftBasic.py: Part of the CcpNmr Analysis program

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
from ccpnmr.analysis.core.ChemicalShiftRef import REFDB_SD_MEAN

from math import exp

# TBD DNA/RNA residue probs
#
# Sanity checks:
#   Only on amide
#   Gly CA - ignore CB
#

ROOT_TWO_PI = 2.506628274631
PROTEIN_MOLTYPE = 'protein'
REF_STORE_DICT = {}
CHEM_ATOM_REF_DICT = {}

def makeRandomCoilShiftList(molSystems):
  """
  Make a synthetic chemical shift list using random coil values,
  adjusting protein backbone values for sequence where approprate.
 
  .. describe:: Input
  
  MolSystem.MolSystem

  .. describe:: Output
  
  Nmr.ShiftList
  """
  
  from ccpnmr.analysis.core.AssignmentBasic import assignAtomsToRes
  from ccpnmr.analysis.core.MoleculeBasic import DEFAULT_ISOTOPES, getRandomCoilShift
  
  project = list(molSystems)[0].root
  done = set()
  nmr = project.currentNmrProject
  shiftList = nmr.newShiftList(isSimulated=True, unit='ppm', name='Random coil')
  resonanceDict = {}
  
  for molSystem in molSystems:
    for chain in molSystem.sortedChains():
      if chain.molecule in done:
        continue
 
      done.add(chain.molecule)
      residues = chain.sortedResidues()
      n = len(residues)
 
      for i, residue in enumerate(residues):
        context = [None] * 5
        for k, j in enumerate(range(i-2, i+3)):
          if 0 <= j < n:
            context[k] = residues[j]
 
        atomSets = set()
        for atom in residue.atoms:
          atomSet = atom.atomSet
 
          if atomSet:
            atomSets.add(atomSet)
 
        for atomSet in atomSets:
          atom = atomSet.findFirstAtom()
          chemAtom = atom.chemAtom
          value = getRandomCoilShift(chemAtom, context)
 
          if value is not None:
            resonanceSet = atomSet.findFirstResonanceSet()
 
            if resonanceSet:
              resonances = list(resonanceSet.resonances)
              index = list(resonanceSet.atomSets).index(atomSet)
              resonance = resonances[min(index, len(resonances)-1)]
 
            else:
              isotope = DEFAULT_ISOTOPES.get(atom.chemAtom.elementSymbol)
 
              if isotope is None:
                continue
 
              resonance = nmr.newResonance(isotopeCode=isotope)
              resonanceSet = assignAtomsToRes((atomSet,), resonance)
 
            resonanceDict[resonance] = value
  
  for resonance in resonanceDict:
    value = resonanceDict[resonance] 
    shift = shiftList.newShift(value=value, resonance=resonance)
  
  return shiftList
            
def chemShiftBasicMacro(argServer):
    
  shiftValues = [8.21, 121.4, 53.9, 38.2]
  elements = ['H','N','C','C']
  atomTypes = None
  
  print (shiftValues, 'Asp', elements, atomTypes)
  
def getAtomProbability(ccpCode, atomName, shiftValue, molType=PROTEIN_MOLTYPE):

  shiftRefs = REFDB_SD_MEAN.get((molType, ccpCode))
    
  if not shiftRefs:
    return 

  stats = shiftRefs.get(atomName)
  if not stats:
    return

  mean, sd, pMissing, bound = stats
  d = shiftValue-mean
  e = d/sd   
  p = exp(-0.5*e*e)/(sd*ROOT_TWO_PI)
 
  return p  

def lookupAtomProbability(project, ccpCode, atomName, shiftValue,
                          molType=PROTEIN_MOLTYPE):
  
  if molType == PROTEIN_MOLTYPE:
    if (ccpCode == 'Gly') and (atomName == 'HA'):
      atomName = 'HA2' 
  
    value = getAtomProbability(ccpCode, atomName, shiftValue)
  
  if value is None:
    chemAtomNmrRef = getChemAtomNmrRef(project, atomName, ccpCode, molType)
 
    if not chemAtomNmrRef:
      return 
 
    value = getAtomRefLikelihood(chemAtomNmrRef, shiftValue)
 
  return value


def getResidueProbability(ppms, ccpCode, elements, atomTypes=None, ppmsBound=None,
                          prior=0.05, molType=PROTEIN_MOLTYPE, cutoff=1e-10):

  # Use refExperiment info
  # Use bound resonances info

  shiftRefs = REFDB_SD_MEAN.get((molType, ccpCode))
    
  if not shiftRefs:
    return None
  
  if not atomTypes:
    atomTypes = [None] * len(ppms)
  
  if not ppmsBound:
    ppmsBound = [None] * len(ppms)
    
  
  atomData = [(x, shiftRefs[x]) for x in shiftRefs.keys()]
  
  data = []
  dataAppend = data.append
  for i, ppm in enumerate(ppms):
    element = elements[i]
    atomType = atomTypes[i] or ()
    ppmB = ppmsBound[i]
    n = 0
    
    
    for j, (atomName, stats) in enumerate(atomData):
      if element != atomName[0]:
        continue 
      
      if atomType and (atomName not in atomType):
        for atomType0 in atomType:
          if (atomType0[-1] == '*') and (atomName == atomType0[:-1]):
            # Useful fior HE*, H* etc
            break
        else:
          continue
      
      mean, sd, pMissing, bound = stats
      d = ppm-mean
      
      if (not atomType) and (abs(d) > 5*sd):
        continue
      
      e = d/sd   
      p = exp(-0.5*e*e)/(sd*ROOT_TWO_PI)
      
      if bound and (ppmB is not None):
        boundData = shiftRefs.get(bound)
	
	if boundData:
          meanB, sdB, pMissingB, boundB = boundData      
	  dB = ppmB-meanB
          eB = dB/sdB   
          pB = exp(-0.5*eB*eB)/(sdB*ROOT_TWO_PI)
      
          p = (p*pB) ** 0.5
      
      if (not atomType) and (p < cutoff):
        continue
      
      
      dataAppend((i,j,p))
      n += 1
    
    if n == 0:
      return 0.0  
  
  groups = [set([node,]) for node in data if node[0] == 0]
  
  while data:
    node = data.pop()
    i, j, p = node
    
    for group in groups[:]:
      for node2 in group:
        i2, j2, p2 = node2
       
        if (i == i2) or (j == j2):
          break
      
      else:
        newGroup = group.copy()
        newGroup.add(node)
        groups.append(newGroup)
 
  probTot = 0.0
  for group in groups:
  
    if len(group) != len(ppms):
      continue
    
    found = set([])
    prob = 1.0
    for i,j, p in group:
      found.add(j)  
      prob *= p
    
    #for k, datum in enumerate(atomData):
    #  atomName, stats = datum
    #  pMissing = stats[2]
    #  
    #  if k in found:
    #    prob *= 1-pMissing
    #  else:
    #    prob *= pMissing
    
    if found:
      probTot += prob
      
  return probTot
  
def getShiftsChainProbabilities(shifts, chain):

  probDict = {}
  getProb = getShiftsResidueProbability
  priors = getChainResTypesPriors(chain)
  
  ccpCodes = set(getChainResidueCodes(chain))
  
  total = 0.0
  for ccpCode, molType in ccpCodes:
    prob = getProb(shifts, ccpCode, priors[ccpCode], molType)
    probDict[ccpCode] = prob
    
    if prob is not None:
      total += prob

  if not total:
    total = 1.0

  for ccpCode, molType in ccpCodes:
    if probDict[ccpCode] is None:
      probDict[ccpCode] = 1.0 / len(chain.residues)
    else:
      probDict[ccpCode] /= total
  
  # Have to do this until have stats at var level
  if 'Cyss' in probDict:
    if 'Cys' in probDict:
      if probDict['Cyss'] > probDict['Cys']:
        probDict['Cys'] = probDict['Cyss']
      
    else:
      probDict['Cys'] = probDict['Cyss']
    
    del probDict['Cyss']
  
  return probDict

def getSpinSystemChainProbabilities(spinSystem, chain, shiftList):

  probDict = {}
  getProb = getSpinSystemResidueProbability
  priors = getChainResTypesPriors(chain)
  
  ccpCodes = set(getChainResidueCodes(chain))
  
  for ccpCode, molType in ccpCodes:
    probDict[ccpCode] = getProb(spinSystem, shiftList, ccpCode,
                                priors[ccpCode], molType)
  
  return probDict

def getChainResidueCodes(chain):

  ccpCodes = []
  for residue in chain.residues:
    ccpCode = residue.ccpCode
    if (ccpCode == 'Cys') and (residue.descriptor == 'link:SG'):
      ccpCode = 'Cyss'

    ccpCodes.append((ccpCode, residue.molType))

  return ccpCodes

def getChainResTypesPriors(chain):

  priors = {}
  
  ccpCodes = [x[0] for x in getChainResidueCodes(chain)]
  n = float(len(ccpCodes))
  
  for ccpCode in set(ccpCodes):
    priors[ccpCode] = ccpCodes.count(ccpCode)/n
  
  return priors

def getShiftsResidueProbability(shifts, ccpCode, prior=0.05, molType=PROTEIN_MOLTYPE):

  ppms = []
  boundPpms = []
  elements = []
  atomTypes = []
  ppmsAppend = ppms.append
  boundPpmsAppend = boundPpms.append
  elementsAppend = elements.append
  atomTypesAppend = atomTypes.append
  
  betaBranch = set(['Val','Ile','Thr'])
  
  for shift in shifts:
    resonance = shift.resonance
    isotope = resonance.isotope
    
    if isotope:
      element = isotope.chemElement.symbol
      
      if element == 'H':
        bound = resonance.findFirstCovalentlyBound()
      else:
        bound = resonance.findFirstCovalentlyBound(isotopeCode='1H')
    
      if bound:
        shift2 = bound.findFirstShift(parentList=shift.parentList)
	if shift2:
	  ppm2 = shift2.value
	else:
	  ppm2 = None
      else:
        ppm2 = None
      
      
      assignNames = resonance.assignNames or set([])
      
      if (not assignNames) and resonance.peakDimContribs:
	refExpDimRefs = set([])

	for contrib in resonance.peakDimContribs:
	  refExpDimRef = contrib.peakDim.dataDimRef.expDimRef.refExpDimRef
  	  if refExpDimRef:
	    refExpDimRefs.add(refExpDimRef)
	  
	for refExpDimRef in refExpDimRefs:
  	  expMeasurement = refExpDimRef.expMeasurement
  	  atomSites = expMeasurement.atomSites
 
  	  for atomSite in atomSites:
  	    name = atomSite.name
 
  	    if name == 'CO':
  	      name == 'C'

  	    elif name in ('H','N',): # Not specific sites
  	      continue
              
  	    elif (name == 'HA') and (ccpCode == 'Gly'):
  	      name = 'HA2'

  	    elif (name == 'HB') and (ccpCode not in betaBranch):
  	      name = 'HB2'
 
  	    elif name in ('C','Cali'):
  	      for expTransfer in atomSite.expTransfers:
  	        if expTransfer.transferType in ('onebond','CP'):
  	          atomSites2 = list(expTransfer.atomSites)
  	          atomSites2.remove(atomSite)
  	          name2 = atomSites2[0].name
 
  	          if (name2 == 'CA') and (ccpCode != 'Gly'):
  	            name = 'CB'
  	            break
  	          elif name2 == 'CO':
  	            name = 'CA'
  	            break

  	      else:
  	        continue

            assignNames.add(name)

      boundPpmsAppend(ppm2)
      ppmsAppend(shift.value)
      elementsAppend(element)
      atomTypesAppend(assignNames)
    
  prob = getResidueProbability(ppms, ccpCode, elements, atomTypes,
                               boundPpms, prior, molType)
  
  return prob

def getSpinSystemResidueProbability(spinSystem, shiftList, ccpCode,
                                    prior=0.05, molType=PROTEIN_MOLTYPE):

  ppms = []
  elements = []
  atomTypes = []
  ppmsAppend = ppms.append
  elementsAppend = elements.append
  atomTypesAppend = atomTypes.append
  
  for resonance in spinSystem.resonances:
    
    isotope = resonance.isotope
    if isotope:
    
      shift = resonance.findFirstShift(parentList=shiftList)
      if shift:
        ppmsAppend(shift.value)
        elementsAppend(isotope.chemElement.symbol)
        atomTypesAppend(resonance.assignNames)
 
    
  prob = getResidueProbability(ppms, ccpCode, elements,
                               atomTypes, prior, molType)
  
  return prob

 

def getChemAtomNmrRef(project, atomName, ccpCode, molType=PROTEIN_MOLTYPE,
                      sourceName='RefDB'):

  atomKey = molType + ccpCode + atomName
  chemAtomNmrRef  = CHEM_ATOM_REF_DICT.get(atomKey)


  if not chemAtomNmrRef:
    key = molType + ccpCode
    getRefStore = project.findFirstNmrReferenceStore
    nmrRefStore = REF_STORE_DICT.get(key, getRefStore(molType=molType,ccpCode=ccpCode)) 
    
    if nmrRefStore:
      REF_STORE_DICT[key] = nmrRefStore
      chemCompNmrRef = nmrRefStore.findFirstChemCompNmrRef(sourceName=sourceName)
      if chemCompNmrRef:
        chemCompVarNmrRef = chemCompNmrRef.findFirstChemCompVarNmrRef(linking='any',descriptor='any')
        
        if chemCompVarNmrRef:
          for chemAtomNmrRef1 in chemCompVarNmrRef.chemAtomNmrRefs:
            if atomName == chemAtomNmrRef1.name:
              chemAtomNmrRef = chemAtomNmrRef1
              CHEM_ATOM_REF_DICT[atomKey] = chemAtomNmrRef
              break
        else:
          return
      else:
        return
    else:
      return
 
  if not chemAtomNmrRef:
    atomName2 = atomName[:-1]
    for chemAtomNmrRef1 in chemCompVarNmrRef.chemAtomNmrRefs:
      if atomName2 == chemAtomNmrRef1.name:
        chemAtomNmrRef = chemAtomNmrRef1
        CHEM_ATOM_REF_DICT[atomKey] = chemAtomNmrRef
        break

  return chemAtomNmrRef

def getAtomRefLikelihood(chemAtomNmrRef, shiftValue):

  distribution  = chemAtomNmrRef.distribution
  refPoint      = chemAtomNmrRef.refPoint
  refValue      = chemAtomNmrRef.refValue
  valuePerPoint = chemAtomNmrRef.valuePerPoint
  N = len(distribution)
  f = 1/3.0
  
  refDelta = shiftValue-refValue
  point    = int(refPoint + (refDelta/valuePerPoint))

  if point < 0:
    return 0.0

  elif point >= N:
    return 0.0

  else:
    #return distribution[int(point)]
    return (distribution[max(point-1,0)] + distribution[point] + distribution[min(point+1,N-1)]) * f


if __name__ == '__main__':

  shiftValues =  [8.868, 112.28, 46.065, 3.982, 4.408]

  elements = ['H','N','C','H','H']
  #atomTypes = [('CB',), ('CA',), ('N',), ('H',)]
  atomTypes = [ ('H',), ('N',), set(['CA']), set(), set()]
  vals = {}
  ccpCodes = ['Ala','Cys','Asp','Glu',
 	      'Phe','Gly','His','Ile',
	      'Lys','Leu','Met','Asn',
 	      'Pro','Gln','Arg','Ser',
 	      'Thr','Val','Trp','Tyr',]
  
  for ccpCode in ccpCodes:
    vals[ccpCode] = getResidueProbability(shiftValues, ccpCode, elements, atomTypes, prior=0.0578512396694,
                                          ppmsBound=[112.28, 8.868, None, None, None])
  

  tot = sum(vals.values())
  for ccpCode in ccpCodes:
    v = vals[ccpCode] / tot
    print ccpCode, '%.3f' % v
 




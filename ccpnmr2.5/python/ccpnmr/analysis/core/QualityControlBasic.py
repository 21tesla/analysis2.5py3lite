
LICENSE = """
======================COPYRIGHT/LICENSE START==========================

QualityControlBasic.py: Part of the CcpNmr Analysis program

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

from math import sqrt, log
from ccpnmr.analysis.core.AssignmentBasic  import getResonanceResidue, makeResonanceGuiName, getResonanceAtomTuple
from ccpnmr.analysis.core.AssignmentBasic  import updateAllShifts, getBoundResonances, getAtomSetShifts
from ccpnmr.analysis.core.ChemicalShiftBasic import getChemAtomNmrRef, lookupAtomProbability
from ccpnmr.analysis.core.PeakBasic        import getPeakVolume, getPeakHeight, findSymmetryPeaks, searchPeaks, getClosestPeak, getPeaksOverlapScore, findClosePeaks
from ccpnmr.analysis.core.ExperimentBasic  import getOnebondDataDims, getThroughSpacePeakLists
from ccpnmr.analysis.core.ExperimentBasic  import getThroughSpaceDataDims, initExpTransfers
from ccpnmr.analysis.core.MoleculeBasic    import areAtomsBound, getResidueCode
from ccpnmr.analysis.core.Util             import getAnalysisDataDim


def analysePeaks(peakList):


  # Row per peak - Cols: height vs PL mean, vol vs PL mean,
  # sign check, F1 shift delta, F2 shift delta ++

  meanVolume = 0.0
  meanHeight = 0.0
  meanVolume2 = 0.0
  meanHeight2 = 0.0
  nVolume = 0
  nHeight = 0

  pos = 0
  tot = 0

  data0 = []
  for peak in peakList.sortedPeaks():
    
    volume = getPeakVolume(peak)
    height = getPeakHeight(peak)
    logVolume = None
    logHeight = None
        
    if volume:
      logVolume = log(abs(volume))
      meanVolume += logVolume
      meanVolume2 += logVolume * logVolume
      nVolume += 1
    else:
      volume = 0.0
    
    if height:
      logHeight = log(abs(height))
      meanHeight += logHeight
      meanHeight2 += logHeight * logHeight
      nHeight += 1
  
      tot += 1
      if height > 0:
        pos += 1
        
    else:
      height = 0.0
  
    data0.append([peak, volume,logVolume, height, logHeight])
  
  
  if nVolume:
    meanVolume /= float(nVolume)
    meanVolume2 /= float(nVolume)
  
  if nHeight:
    meanHeight /= float(nHeight)
    meanHeight2 /= float(nHeight)
  
  heightSd =  sqrt(abs(meanHeight2 - (meanHeight * meanHeight)))
  # abs() due to float point error - v small negs rather than zero 
  volumeSd =  sqrt(abs(meanVolume2 - (meanVolume * meanVolume)))
    
  shiftList = peakList.dataSource.experiment.shiftList
  
  boundDataDims = {}
  for dataDim1, dataDim2 in getOnebondDataDims(peakList.dataSource):
    boundDataDims[dataDim1] = dataDim2
    boundDataDims[dataDim2] = dataDim1

  data = []
  for peak, volume, logVolume, height, logHeight in data0:
  
    deltaLogVolume = None
    deltaLogHeight = None
    
    if volume:
      deltaLogVolume = abs(logVolume - meanVolume)

    if height:
      deltaLogHeight = abs(logHeight - meanHeight)

    maxDeltas = []
    assignErrors = {}
    # Row per peak: - Cols: Causes muli 1bond,
    # causes atom impossible 1bond, Partly assigned onebond,
    # missing/extra assign
    for peakDim in peak.sortedPeakDims():
      
      maxDelta = None
      isotopeCode  = None
      
      if shiftList:
        for contrib in peakDim.peakDimContribs:
          resonance = contrib.resonance
          isotopeCode = resonance.isotopeCode
          shift = resonance.findFirstShift(parentList=shiftList)
          if shift:
            delta = abs(peakDim.realValue - shift.value)
            
            if (maxDelta is None) or (delta > maxDelta):
              maxDelta = delta
      
      maxDeltas.append([maxDelta,isotopeCode])  
      
      for contrib in peakDim.peakDimContribs:
        resonance = contrib.resonance
        # Warn the user if prochirals have similar shifts
        # but only one is linked to the peak.
        resonanceSet = resonance.resonanceSet
        if resonanceSet:
          if len(resonanceSet.atomSets) > 1:
            resonances2 = list(resonanceSet.resonances)
            resonances2.remove(resonance)
            
            for resonance2 in resonances2:
              if peakDim.findFirstPeakDimContrib(resonance=resonance2):
                continue
              
              for contrib2 in resonance2.peakDimContribs:
                if contrib2.peakDim.peak.peakList is peakList:
                  break
              
              else:
                shift = resonance.findFirstShift(parentList=shiftList)
                shift2 = resonance2.findFirstShift(parentList=shiftList)
                if shift and shift2:
                  delta = abs(shift.value-shift2.value)
                  analysisDataDim = getAnalysisDataDim(peakDim.dataDim)
 
                  if analysisDataDim and (delta < analysisDataDim.assignTolerance/5.0):
                    name = makeResonanceGuiName(resonance2, fullName=False)
                    msg = 'Also expected dim %d %s assignment' % (peakDim.dim,name)
                    assignErrors[msg] = True

      boundDataDim = boundDataDims.get(peakDim.dataDim)
      if boundDataDim:
        peakDim2 = peak.findFirstPeakDim(dataDim=boundDataDim)
    
        if peakDim.peakDimContribs:
          if peakDim2.peakDimContribs:
            
            # Count number of atomSets to see if result makes sense
            atomSets1 = set()
            atomSets2 = set()
            for contrib2 in peakDim2.peakDimContribs:
              resonanceSet = contrib2.resonance.resonanceSet
              if resonanceSet:
                for atomSet in resonanceSet.atomSets:
                  atomSets2.add(atomSet)
          
            for contrib in peakDim.peakDimContribs:
              resonance = contrib.resonance
              resonanceSet = resonance.resonanceSet
              name = makeResonanceGuiName(resonance)
              bound = getBoundResonances(resonance, recalculate=True)
              impossibles  = []
              
              if resonanceSet:
                for atomSet in resonanceSet.atomSets:
                  atomSets1.add(atomSet)

              if bound:
        
                nBound = len(bound)
                
                for reson in bound:
                  # Correct for multiatom atomSets (CH3 resonances e,g,)
                  resSet = reson.resonanceSet
                  if resSet:
                    ll = [1]
                    for atomSet in resSet.atomSets:
                      length = len(atomSet.atoms)
                      cas = atomSet.findFirstAtom().chemAtom.chemAtomSet
                      if not cas or cas.isEquivalent is not None:
                        # Test excludes e.g. Tyr and Phe side chains
                        # that otherwise give spurious errors
                        ll.append(length)
                    # Add additional number of atoms for each bound resonance
                    nBound += min(ll) - 1
                    
                if isotopeCode in ('1H','2H','19F'):
                  if nBound > 1:
                    assignErrors['%s mutiple %s bonds' % (name,isotopeCode)] = True
 
                elif resonanceSet:
                  chemAtom = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().chemAtom
                  if nBound > len(chemAtom.chemBonds):
                    assignErrors['%s excessive bonds' % name] = True
                    
                okAtoms = False
                foundBound = False
                for contrib2 in peakDim2.peakDimContribs:
                  resonance2 = contrib2.resonance
                  
                  if resonance2 in bound:
                    resonanceSet2 = resonance2.resonanceSet
                    foundBound = True
                    if resonanceSet and resonanceSet2:
                      if not areResonanceSetAtomsBound(resonanceSet, resonanceSet2):
                        names = [name,makeResonanceGuiName(resonance2)]
                        names.sort()
                        impossibles.append(names)
                        
                      else:
                        okAtoms = True
                        
                if not okAtoms:
                  for names in impossibles:
                    assignErrors['Impossible bond %s-%s' % tuple(names)] = True
                
                if not foundBound:
                  assignErrors['%s no bound partner' % name] = True
            
            if isotopeCode in ('1H','2H','19F'):
              if len(atomSets2) > len(atomSets1):
                assignErrors['Only %s H to bind %s Non-H' % (len(atomSets1), len(atomSets2))] = True
                
          else:
            assignErrors['Dim %d empty' % peakDim2.dataDim.dim] = True
        
        elif peakDim2.peakDimContribs:
          assignErrors['Dim %d empty' % peakDim.dataDim.dim] = True
      

    sortedPeakDims = peak.sortedPeakDims()
    annotation = ' '.join([pd.annotation or '-' for pd in sortedPeakDims])
    
    locations = []
    for pd in sortedPeakDims:
      if pd.value is None:
        if pd.position is None:
          locations.append('Error')
        else:
          locations.append('%.2f' % pd.position)
          
      else:
        locations.append('%.2f' % pd.value)
    
    location   = ' '.join(locations)

    nearestPeak, nearestDistance = findClosestOtherPeak(peak, scale=2.0)
    symmetryPeaks = findSymmetryPeaks(peak)

    datum = [peak.serial,annotation,assignErrors.keys(),location,
             volume,deltaLogVolume,
             height,deltaLogHeight,maxDeltas,
             nearestPeak, nearestDistance, symmetryPeaks]
    data.append([peak, datum])

  posProp = 1.0
  if tot:
    posProp = pos/float(tot)

  return posProp, volumeSd, heightSd, data

def findClosestOtherPeak(peak, scale=1.0):

  peakList = peak.peakList
  spectrum = peakList.dataSource
  analysisSpectrum = spectrum.analysisSpectrum
  tolerances = spectrum.numDim * [0]
  for analysisDataDim in analysisSpectrum.analysisDataDims:
    dim = analysisDataDim.dataDim.dim
    tolerances[dim-1] = analysisDataDim.assignTolerance

  peaks = findClosePeaks(peak, peakList, tolerances)
  peak2 = getClosestPeak(peak, peaks, tolerances)
  if peak2:
    dist = getPeaksOverlapScore(peak, peak2, tolerances)
  else:
    dist = None

  return peak2, dist

def areResonanceSetAtomsBound(resonanceSet1, resonanceSet2):

  for atomSet in resonanceSet1.atomSets:
    for atom in atomSet.atoms:
      residue = atom.residue
      for atomSet2 in resonanceSet2.atomSets:
        for atom2 in atomSet2.atoms:
          if areAtomsBound(atom, atom2):
            return True

  return False

def analyseChemicalShifts(shiftList):
    
    updateAllShifts(shiftList)
    
    # T1, shifts -  Row per resonance - Cols: shift, shiftSD,
    # maxPeak delta, SD per spectrum ++
    data = []
    
    duplicateAssign = {}
    data2 = []
    if shiftList:
      project = shiftList.root
    
      for shift in shiftList.measurements:
        resonance = shift.resonance
        assignTuple = getResonanceAtomTuple(resonance)
        data2.append(['%s%s%8.8s%s' % assignTuple,shift])
        duplicateAssign[assignTuple] = duplicateAssign.get(assignTuple, 0) + 1
    
    data2.sort()
    for name, shift in data2:
      resonance  = shift.resonance
      deltaMax   = None
      nContribs  = 0
      bmrbMean   = None
      randomCoil = None
      typeScore  = None
      
      resonanceSet = resonance.resonanceSet
      if resonanceSet:
        atomSet = resonanceSet.findFirstAtomSet()
        residue = atomSet.findFirstAtom().residue
        ccpCode = residue.ccpCode
        molType = residue.molResidue.molType
        atomName = atomSet.name
        chemAtomNmrRef = getChemAtomNmrRef(project, atomName, ccpCode, molType=molType)
        if chemAtomNmrRef is None:
          atomName = atomSet.findFirstAtom().name
          chemAtomNmrRef = getChemAtomNmrRef(project, atomName, ccpCode, molType=molType)
         
        if chemAtomNmrRef:
        
          # Horrid kludge until we have chem shift ref info per var
          if (ccpCode == 'Cys') and ('link:SG' in residue.chemCompVar.descriptor):
            ccpCode = 'Cyss'
                
          typeScore  = lookupAtomProbability(project, ccpCode, atomName, shift.value, molType)
          bmrbMean   = chemAtomNmrRef.meanValue
          randomCoil = chemAtomNmrRef.randomCoilValue
      
      for contrib in resonance.peakDimContribs:
        peakDim = contrib.peakDim
        shiftList1 = peakDim.peak.peakList.dataSource.experiment.shiftList
        if shiftList1 and (shiftList1 is shiftList):
          nContribs +=1
          delta = abs(peakDim.realValue-shift.value)
          if (deltaMax is None) or (delta > deltaMax):
            deltaMax = delta
            
      boundWarn = False      
      bound = getBoundResonances(resonance, recalculate=True)
      if bound:
        
        nBound = len(bound)
        
        for reson in bound:
          # Correct for multiatom atomSets (CH3 resonances e,g,)
          resSet = reson.resonanceSet
          if resSet:
            ll = [1]
            for atomSet in resSet.atomSets:
              length = len(atomSet.atoms)
              cas = atomSet.findFirstAtom().chemAtom.chemAtomSet
              if not cas or cas.isEquivalent is not None:
                # Test excludes e.g. Tyr and Phe side chains 
                # that otherwise give spurious errors
                ll.append(length)
            # Add additional number of atoms for each bound resonance
            nBound += min(ll) - 1
          
        if resonance.isotopeCode in ('1H','2H','19F'):
          if nBound > 1:
            boundWarn = True
 
        elif resonanceSet:
          chemAtom = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().chemAtom
          if nBound > len(chemAtom.chemBonds):
            boundWarn = True
    
        if resonanceSet and not boundWarn:
          atom = resonanceSet.findFirstAtomSet().findFirstAtom()
          
          for resonance1 in bound:
            resonanceSet1 = resonance1.resonanceSet
            
            if resonanceSet1:
              for atomSet1 in resonanceSet1.atomSets:
                for atom1 in atomSet1.atoms:
                  if areAtomsBound(atom, atom1):
                    break
                    
                else:
                  continue
                break   
              
              else:
                boundWarn = True
    
      assignTuple = getResonanceAtomTuple(resonance)
      sameResBound = []
      otherBound   = []
      residue = getResonanceResidue(resonance)
      for resonance1 in bound:
        residue1 = getResonanceResidue(resonance1)
        if residue1 is residue:
          sameResBound.append(makeResonanceGuiName(resonance1,fullName=False))
        else:
          otherBound.append(makeResonanceGuiName(resonance1))
      
      boundResonances = ''
      if residue:
        ccpCode = getResidueCode(residue)
        resName1 = '%d%s' % (residue.seqCode,ccpCode)
        
        if len(sameResBound) > 1:
          boundResonances += '%s[%s] ' % (resName1,','.join([x for x in sameResBound]))
        elif sameResBound:
          boundResonances += '%s%s ' % (resName1,sameResBound[0])
      
      else:
         boundResonances += ','.join([x for x in sameResBound])
         
      boundResonances += ','.join([x for x in otherBound])
     
      isDuplicate = False
      if duplicateAssign.get(assignTuple, 0) > 1:
        isDuplicate = True
        
      resName = makeResonanceGuiName(resonance)
      datum = [resonance.serial,resonance.isotopeCode,resName,boundResonances,
               bmrbMean,randomCoil,typeScore,shift.value,shift.error,
               deltaMax,nContribs,isDuplicate,boundWarn]
      
      data.append([shift, datum])
      
    return data  

def analyseNoeAssignments(molSystem, peakLists=None):


  # Row per residue - cols are num NOES: total, intra, inter,
  # seq, short, short non-seq, long, inter chain, intra chain, contacted res
  
  project = molSystem.root
  
  if not peakLists:
    peakLists = getThroughSpacePeakLists(project)
    
  total = {}
  intra = {}
  inter = {}
  seqen = {}
  short = {}
  nonSeqShort = {}
  longr = {}
  interCh = {}
  intraCh = {}
  contact = {}
  
  for peakList in peakLists:
    dataDims = getThroughSpaceDataDims(peakList.dataSource)
    
    if not dataDims:
      initExpTransfers(peakList.dataSource.experiment)
      dataDims = getThroughSpaceDataDims(peakList.dataSource)
    
    if not dataDims:
      experiment = peakList.dataSource.experiment
      print "Something is wrong with ExpTransfer setup for experiment %s" % (experiment.name)
      continue
    
    dims = set([dd.dim for dd in dataDims])
  
    for peak in peakList.peaks:
      dict = {}

      peakDims = [pd for pd in peak.peakDims if pd.dim in dims]
      
      dimGrouped = []
      dimUngrouped = []
      for peakDim in peakDims:
        grouped = []
        ungrouped = []
      
        for contrib in peakDim.peakDimContribs:
          resonanceSet = contrib.resonance.resonanceSet
          
          if not resonanceSet:
            continue
          
          if contrib.peakContribs:
            grouped.append((contrib, resonanceSet))
          else:
            ungrouped.append((contrib, resonanceSet))
        
        dimGrouped.append(grouped)
        dimUngrouped.append(ungrouped)
      
      for contribs in dimGrouped:
        for contrib, resonanceSet in contribs:
          residue = resonanceSet.findFirstAtomSet().findFirstAtom().residue

          for peakContrib in contrib.peakContribs:
            if dict.get(peakContrib) is None:
              dict[peakContrib] = []
            dict[peakContrib].append(residue)
      
      if (len(dimUngrouped[0]) < 2) or (len(dimUngrouped[-1]) < 2):
        for contribs in dimUngrouped:
          for contrib, resonanceSet in contribs:
            residue = resonanceSet.findFirstAtomSet().findFirstAtom().residue
        
            if dict.get(None) is None:
              dict[None] = []
            dict[None].append(residue)
  
      for peakContrib in dict.keys():
        residues = dict[peakContrib]
        
        if len(residues) > 1:
          for i in range(len(residues)-1):
            residueI = residues[i]
            chainI   = residueI.chain
            if contact.get(residueI) is None:
              contact[residueI] = {}
              
            for j in range(i+1,len(residues)):
              residueJ = residues[j]
              chainJ   = residueJ.chain
      
              delta = abs(residueI.seqCode-residueJ.seqCode)
              
              if contact.get(residueJ) is None:
                contact[residueJ] = {}
              
              total[chainI] = total.get(chainI, 0) + 1
              if chainI is not chainJ:
                total[chainJ] = total.get(chainJ, 0) + 1
              
              if residueI is residueJ:
                intra[residueI] = intra.get(residueI, 0) + 1
                total[residueI] = total.get(residueI, 0) + 1
                intra[chainI]   = intra.get(chainI, 0) + 1
              else:
                contact[residueJ][residueI] = True
                contact[residueI][residueJ] = True
                total[residueI] = total.get(residueI, 0) + 1
                total[residueJ] = total.get(residueJ, 0) + 1
                inter[residueI] = inter.get(residueI, 0) + 1
                inter[residueJ] = inter.get(residueJ, 0) + 1
                inter[chainI] = inter.get(chainI, 0) + 1
                if chainI is not chainJ:
                  inter[chainJ] = inter.get(chainJ, 0) + 1
              
              if chainI is chainJ:
                intraCh[chainI] = intraCh.get(chainI, 0) + 1
                intraCh[residueI] = intraCh.get(residueI, 0) + 1
                if residueI is not residueJ:
                  intraCh[residueJ] = intraCh.get(residueJ, 0) + 1
 
                if delta < 5:
                
                  if delta > 0:
                    short[chainI] = short.get(chainI, 0) + 1
                    short[residueI] = short.get(residueI, 0) + 1
                    if residueI is not residueJ:
                      short[residueJ] = short.get(residueJ, 0) + 1
 
                    if delta > 1:
                      nonSeqShort[chainI]   = nonSeqShort.get(chainI, 0) + 1
                      nonSeqShort[residueI] = nonSeqShort.get(residueI, 0) + 1
                      nonSeqShort[residueJ] = nonSeqShort.get(residueJ, 0) + 1
 
                    elif delta == 1:
                      seqen[residueI] = seqen.get(residueI, 0) + 1
                      seqen[residueJ] = seqen.get(residueJ, 0) + 1
                      seqen[chainI]   = seqen.get(chainI, 0) + 1
 
                else:
                  longr[chainI]   = longr.get(chainI, 0) + 1
                  longr[residueI] = longr.get(residueI, 0) + 1
                  longr[residueJ] = longr.get(residueJ, 0) + 1
                 
 
              else:
                longr[residueI] = longr.get(residueI, 0) + 1
                longr[residueJ] = longr.get(residueJ, 0) + 1
                longr[chainI]   = longr.get(chainI, 0) + 1
                longr[chainJ]   = longr.get(chainJ, 0) + 1
                
                interCh[residueI] = interCh.get(residueI, 0) + 1
                interCh[residueJ] = interCh.get(residueJ, 0) + 1
                interCh[chainI]   = interCh.get(chainI, 0) + 1
                interCh[chainJ]   = interCh.get(chainJ, 0) + 1

  
  if len(molSystem.chains) > 1:
    residues = [('%s%5.5d' % (r.chain.code,r.seqCode),
                 '%s%d%s'% (r.chain.code,r.seqCode,r.ccpCode),r)  for r in contact.keys()]
  else:
    residues = [('%5.5d' % (r.seqCode,),
                 '%d%s' % (r.seqCode,r.ccpCode),r) for r in contact.keys()]
  
  residues.sort()
  
  data = []
  
  for key, name, residue in residues:
    datum = [name,
             total.get(residue, 0),
             intra.get(residue, 0),
             inter.get(residue, 0),
             seqen.get(residue, 0),
             short.get(residue, 0),
             nonSeqShort.get(residue, 0),
             longr.get(residue, 0),
             intraCh.get(residue, 0),
             interCh.get(residue, 0),
             contact[residue].keys()]
             
    data.append([residue, datum])        

  if len(molSystem.chains) > 1:
    for residue, datum in data:
      residues2 = [('%s%5.5d' % (r.chain.code,r.seqCode),
                    '%s%d'% (r.chain.code,r.seqCode)) for r in datum[-1]]
                    
      residues2.sort()
      datum[-1] = ' '.join([x[1] for x in residues2])
     
 
  else:
    for residue, datum in data:
      residues2 = [('%5.5d' % (r.seqCode),
                   '%d' % (r.seqCode)) for r in datum[-1]]
      residues2.sort()
      datum[-1] = ' '.join([x[1] for x in residues2])
 
  chains = [(chain.code, chain) for chain in molSystem.chains]
  chains.sort()
  chains.reverse()
  
  for code, chain in chains:
    datum = ['Chain %s' % code,]
    datum.append(total.get(chain, 0))
    datum.append(intra.get(chain, 0))
    datum.append(inter.get(chain, 0))
    datum.append(seqen.get(chain, 0))
    datum.append(short.get(chain, 0))
    datum.append(nonSeqShort.get(chain, 0))
    datum.append(longr.get(chain, 0))
    datum.append(intraCh.get(chain, 0))
    datum.append(interCh.get(chain, 0))
    datum.append(None)
    data.insert(0, [None, datum])      
     
  return data

def updateShiftSatsList(array, minPpm, meanPpm, maxPpm):  
            
  if meanPpm:
    if (array[2][0] is None) or (minPpm < array[2][0]):
      array[2][0] = minPpm
    
    if array[2][1] is None:
      array[2][1] = meanPpm
    else:
      array[2][1] += meanPpm
      
    if (array[2][2] is None) or (maxPpm > array[2][2]):
      array[2][2] = maxPpm

    array[1] += 1

def analyseAssignmentCompleteness(molSystem, shiftList, residueSelection=None,
                                  excludeWaterExchangeable=True):

  from ccpnmr.analysis.core.MoleculeBasic import greekSortAtomNames
  # ['Category','Available','Assigned','% Assigned']
  # By residue also
            
  elements  = {}
  amide     = [0,0,[None,None,None]]
  backbone  = [0,0,[None,None,None]]
  backboneX = [0,0,[None,None,None]]
  sideH     = [0,0,[None,None,None]]
  sideX     = [0,0,[None,None,None]]
  riboseH   = [0,0,[None,None,None]]
  riboseX   = [0,0,[None,None,None]]
  sugPosH   = [0,0,[None,None,None]]
  sugPosX   = [0,0,[None,None,None]]
  atomTypes = {}
  residues  = {}
  residueSelection = set(residueSelection or [])
  
  project = molSystem.root
  
  for chain in molSystem.chains:
  
    for residue in chain.sortedResidues():
      if residueSelection and (residue not in residueSelection):
        continue
    
      molType = residue.molResidue.molType
      
      if residues.get(residue.ccpCode) is None:
        residues[residue.ccpCode] = [0,0,[None,None,None]] 
        
      residues[residue.ccpCode][0] += 1
      
      residueAssigned = False
      dict = {}
      
      for atom in residue.atoms:
        
        if atom.chemAtom.waterExchangeable and excludeWaterExchangeable:
          continue
        
        atomSet = atom.atomSet
        if atomSet:
          if dict.get(atomSet) is None:
            dict[atomSet] = True
            
            minPpm  = None
            meanPpm = None
            maxPpm  = None
            
            if atomSet.resonanceSets: # could also check for peak links
              
              shifts  = getAtomSetShifts(atomSet, shiftList=shiftList)
              sum  = 0.0
              for shift in shifts:
                value = shift.value
              
                if (minPpm is None) or (value < minPpm):
                  minPpm = value
                  
                if (maxPpm is None) or (value > maxPpm):
                  maxPpm = value
                  
                sum += value 
              
              if shifts:
                isAssigned = True
                meanPpm = sum/ float(len(shifts))
              else:
                isAssigned = False
              
            else:
              isAssigned = False
              
            element  = atom.chemAtom.elementSymbol
            if element in ('H','C','N','P','F'):
            
              atomType = element
              if len(element) < len(atom.name):
                atomType += atom.name[len(element):][0]
 
              if elements.get(element) is None:
                elements[element] = [0,0,[None,None,None]]
              elements[element][0] += 1
              
 
              if atomTypes.get(atomType) is None:
                atomTypes[atomType] = [0,0,[None,None,None]]
              atomTypes[atomType][0] += 1
              
              if molType == 'protein':
                # wb104: 4 Aug 2014: added HA,HA2,HA3 to the backbone
                if atom.name in ('H','N','C','CA','HA','HA2','HA3'):
                  
                  if element != 'H':
                    backboneX[0] +=1
                    if isAssigned:
                      backboneX[1] += 1
                   
                  backbone[0] +=1
                  if isAssigned:
                    backbone[1] += 1
 
                  if atom.name in ('H','N'):
                    amide[0] += 1
                    if isAssigned:
                      amide[1] += 1
 
                else:
                  if element == 'H':
                    sideH[0] += 1
                    if isAssigned:
                      updateShiftSatsList(sideH, minPpm, meanPpm, maxPpm)
                      
                  else:
                    sideX[0] += 1
                    if isAssigned:
                      updateShiftSatsList(sideX, minPpm, meanPpm, maxPpm)
 
 
              if molType in ('RNA','DNA'):
                if "'" in atom.name:
                  if element == 'H':
                    riboseH[0] += 1
                    if isAssigned:
                      updateShiftSatsList(riboseH, minPpm, meanPpm, maxPpm)
 
                  else:
                    riboseX[0] += 1
                    if isAssigned:
                      updateShiftSatsList(riboseX, minPpm, meanPpm, maxPpm)
 
                else:
                  if element == 'H':
                    sugPosH[0] += 1
                    if isAssigned:
                      updateShiftSatsList(sugPosH, minPpm, meanPpm, maxPpm)
 
                  else:
                    sugPosX[0] += 1
                    if isAssigned:
                      updateShiftSatsList(sugPosX, minPpm, meanPpm, maxPpm) 
 
 
              if isAssigned:
                updateShiftSatsList(elements[element], minPpm, meanPpm, maxPpm)
                updateShiftSatsList(atomTypes[atomType], minPpm, meanPpm, maxPpm)
 
                if not residueAssigned:
                  residues[residue.ccpCode][1] += 1
                  residueAssigned = True
 
  data = []
  
  elementList = elements.keys()
  elementList.sort()
  for e in elementList:
    maxNum, found, shiftData = elements[e]
    minPpm, meanPpm, maxPpm = shiftData
    N = float(maxNum)
    if meanPpm is not None:
      meanPpm /= float(found)
    data.append(['Element ' + e, maxNum, found, 100.0*found/N, minPpm, meanPpm, maxPpm])
  
  classes = [('(Backbone)N+H',amide),
             ('Backbone+H+HA',backbone),
             ('Backbone',backboneX),
             ('Side Chain H',sideH),
             ('Side Chain non-H',sideX),
             ('Ribose H',riboseH),
             ('Ribose non-H',riboseX),
             ('Sugar Pos H',sugPosH),
             ('Sugar Pos non-H',sugPosX)]
  
  for name, stats in classes:
    if stats[0] != 0:
      maxNum, found, shiftData = stats
      minPpm, meanPpm, maxPpm = shiftData
      N = float(maxNum)
      if meanPpm is not None:
        meanPpm /= float(found)
      data.append([name, maxNum, found, 100.0*found/N, minPpm, meanPpm, maxPpm])
  
  atomTypeList = atomTypes.keys()
  atomTypeList = greekSortAtomNames(atomTypeList)
  for atomType in atomTypeList:
    maxNum, found, shiftData = atomTypes[atomType]
    minPpm, meanPpm, maxPpm = shiftData
    N = float(maxNum)
    if meanPpm is not None:
      meanPpm /= float(found)
    
    data.append(['Type ' + atomType, maxNum, found, 100.0*found/N, minPpm, meanPpm, maxPpm])
  
  allMaxNum = 0
  allFound = 0  
  residueList = residues.keys()
  residueList.sort()
  for ccpCode in residueList:
    maxNum, found, shiftData  = residues[ccpCode]
    minPpm, meanPpm, maxPpm = shiftData
    N = float(maxNum)
    if meanPpm is not None:
      meanPpm /= float(found)
    
    allMaxNum += maxNum
    allFound  += found
      
    data.append(['Residue ' + ccpCode, maxNum, found, 100.0*found/N, minPpm, meanPpm, maxPpm])
  
  percentage = 100.0*allFound/float(allMaxNum) if allMaxNum else 0.0
  data.append(['All Residues', allMaxNum, allFound, percentage,None,None,None])
  
  spinSystems = shiftList.nmrProject.resonanceGroups
  maxNum = len(spinSystems)
  N = float(maxNum)
  found = 0
  for spinSystem in spinSystems:
    if spinSystem.residue:
      found += 1
  
  if N:
    data.append(['Spin Systems', maxNum, found, 100.0*found/N, None, None, None])
    
  return data
 


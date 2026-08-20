LICENSE = """
======================COPYRIGHT/LICENSE START==========================

ConstraintBasic.py: Part of the CcpNmr Analysis program

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
import re
import os

from math import sqrt

from ccpnmr.analysis.core.PeakBasic       import pickPeak, getAliasedPeakDimPositions
from ccpnmr.analysis.core.AssignmentBasic import areResonancesProchiral, findMatchingPeakDimShifts
from ccpnmr.analysis.core.AssignmentBasic import getBoundResonances, getOnebondResonance
from ccpnmr.analysis.core.AssignmentBasic import isShiftInRange, assignAtomsToRes, assignResToDim
from ccpnmr.analysis.core.AssignmentBasic import getResonanceLabellingFraction, getResonancePairLabellingFraction
from ccpnmr.analysis.core.ExperimentBasic import findSpectrumDimsByIsotope, getSpectrumIsotopes
from ccpnmr.analysis.core.ExperimentBasic import getOnebondDataDims, getThroughSpaceDataDims
from ccpnmr.analysis.core.ExperimentBasic import getIndirectThroughSpaceIsotopes, getIndirectDataDims
from ccpnmr.analysis.core.MoleculeBasic   import areResonancesBound, getBoundAtoms
from ccpnmr.analysis.core.Util            import getSpectrumNoeDistanceClasses, getAnalysisPeakList
from ccpnmr.analysis.core.UnitConverter   import pnt2ppm, unit_converter
from memops.general.Util             import copySubTree

from ccp.util.NmrExpPrototype import longRangeTransfers

try:
  from memops.gui.MessageReporter import showWarning, showYesNo
  #from memops.gui.DataEntry       import askString
except ImportError:
  from memops.universal.MessageReporter import showWarning, showYesNo
  #from memops.universal.DataEntry       import askString

def updateDistConstraintFromPeakAssign(constraint):
  """
  Update the items/contributions of a distance constraint given
  any new assignments to the peak from which it was originally derived.
  
  .. describe:: Input
  
  List of NmrConstraint.DistanceConstraint
  
  .. describe:: Output
  
  None
  """

  peaks = [contrib.peak for contrib in constraint.peakContribs if contrib.peak]

  resonancePairs = set()
  for peak in peaks:    
    peakDims    = peak.sortedPeakDims()
    expDimRefs  = [pd.dataDimRef.expDimRef for pd in peakDims if pd.dataDimRef]
    experiment  = peak.peakList.dataSource.experiment
    
    for expTransfer in experiment.expTransfers:
      if expTransfer.transferType in longRangeTransfers:
        expDimRefs0 = list(expTransfer.expDimRefs)
        indexA      = expDimRefs.index(expDimRefs0[0])
        indexB      = expDimRefs.index(expDimRefs0[1])
        
        peakDimA = peakDims[indexA]
        peakDimB = peakDims[indexB]
        
        for contribA in peakDimA.peakDimContribs:
          resonanceA = contribA.resonance
          contribsB  = []
          
          for peakContrib in contribA.peakContribs:
            for contribB in peakContrib.peakDimContribs:
              if (contribB.peakDim is peakDimB) and (contribB is not contribA):
                contribsB.append(contribB)
        
          if not contribsB:
            contribsB = peakDimB.peakDimContribs
          
          for contribB in contribsB:
            resonanceB = contribB.resonance
            
            if resonanceA is not resonanceB:
              resonances = [resonanceA,resonanceB]
              resonances.sort()
            
              resonancePairs.add(frozenset(resonances))

  refreshAssign = []
  
  if resonancePairs:  
    constraintSet = constraint.parentList.nmrConstraintStore
  
    for item in constraint.items:
      item.delete()
      
    for resonanceA, resonanceB in resonancePairs:
      fixedResonanceA = getFixedResonance(constraintSet, resonanceA)
      fixedResonanceB = getFixedResonance(constraintSet, resonanceB)
    
      resonanceSetA = resonanceA.resonanceSet
      resonanceSetB = resonanceB.resonanceSet
      fResonanceSetA = fixedResonanceA.resonanceSet
      fResonanceSetB = fixedResonanceB.resonanceSet
    
      if resonanceSetA and fResonanceSetA:
        atomsR = [aSet.atoms for aSet in resonanceSetA.sortedAtomSets()]
        atomsF = [aSet.atoms for aSet in fResonanceSetA.sortedAtomSets()]
        if atomsR != atomsF:
          refreshAssign.append((resonanceA, fixedResonanceA))
      
      elif resonanceSetA and not fResonanceSetA:
        refreshAssign.append((resonanceA, fixedResonanceA))
      
      if resonanceSetB and fResonanceSetB:
        atomsR = [aSet.atoms for aSet in resonanceSetB.sortedAtomSets()]
        atomsF = [aSet.atoms for aSet in fResonanceSetB.sortedAtomSets()]
        if atomsR != atomsF:
          refreshAssign.append((resonanceB, fixedResonanceB))
          
      elif resonanceSetB and not fResonanceSetB:
          refreshAssign.append((resonanceB, fixedResonanceB))
    
      constraint.newDistanceConstraintItem(resonances=[fixedResonanceA,fixedResonanceB])

  if refreshAssign:
    msg = 'Atomic assignments of resonances have changed since this restraint was made. '
    msg += 'Update restraints to relect this?' 
    msg += '(This may affect more than one restraint)'
    
    if showYesNo('Query', msg):
      nmrConstraintStore = constraint.topObject
      
      for resonance, fixedResonance in refreshAssign:
        updateFixedResonance(fixedResonance, resonance)


def updateFixedResonance(fResonance, resonance=None):
  """
  Update the atom assignments of a fixed resonance to reflect a connected
  Nmr.Resonance. Will not remove asslignments if Nmr.Resonance is missing of
  unassigned. A new Nmr.Resonance may optinally be specified.
  
  .. describe:: Input
  
  NmrConstraint.FixedResonance, Nmr.Resonance
  
  .. describe:: Output
  
  None
  """
  
  if not resonance:
    resonance = fResonance.resonance
  
  if resonance:
    resonanceSet = resonance.resonanceSet
    
    if resonance.serial != fResonance.resonanceSerial:
      fResonance.resonanceSerial = resonance.serial
    
    if resonanceSet:
      nmrConstraintStore = fResonance.topObject
      fResonanceSet = fResonance.resonanceSet
      fixedResonances = [fResonance,]
      fResonance2 = None
      
      if fResonanceSet:
        for fResonance0 in fResonanceSet.resonances:
          if fResonance0 is not fResonance:
            fResonance2 = fResonance0
        
        fResonanceSet.delete()
        
      if fResonance2:
        resonance2 = fResonance2.resonance
        
        if (not resonance2) or (not resonance2.resonanceSet) or \
           (resonance2.resonanceSet is resonanceSet):
          fixedResonances.append(fResonance2)
        
        else:
          updateFixedResonance(fResonance2, resonance2)
      
      fixedAtomSets = []
      for atomSet in resonanceSet.atomSets:
        fixedAtomSets.append( getFixedAtomSet(nmrConstraintStore, atomSet.atoms) )
        
      nmrConstraintStore.newFixedResonanceSet(atomSets=fixedAtomSets,
                                              resonances=fixedResonances)
  
    if hasattr(fResonance, 'guiName'):
      del fResonance.guiName
              
def updatePeaksFromConstraints(constraints, replace=False):
  """
  Update the assignments of the peaks that are linked to distance constraints of
  a constraint set, given changes to the constraint items. This function is
  ideal for visualising the NOE disambiguation of procedures like ARIA and
  CYANA. Option to replace any existing assignments; otherwise assignments will
  be added to.

  .. describe:: Input
  
  List of NmrConstraint.DistanceConstraints, Boolean
  
  .. describe:: Output
  
  None
  """

  isotopeWeights = {'1H':1.0,'2H':1.0,'15N':0.15,'13C':0.1}

  boundDimDict = {}

  # find any existing linked peaks
  for constraint in constraints:
    if constraint.className == 'DistanceConstraint':
      peakContribs = constraint.peakContribs
      
      for peakContrib in constraint.peakContribs:
        peak = peakContrib.peak
        
        if peak is None:
          continue
        
        peakList = peak.peakList
        spectrum = peakList.dataSource
        shiftList = spectrum.experiment.shiftList
  
        boundDims = boundDimDict.get(spectrum)
        if boundDims is None:
          boundDims = {}
          
          for dataDim1, dataDim2 in getOnebondDataDims(spectrum):
            boundDims[dataDim1] = dataDim2
            boundDims[dataDim2] = dataDim1
        
          boundDimDict[spectrum] = boundDims
  
        peakDims = peak.sortedPeakDims()
        dimIsotopes = []
        dimBound = []
        dimPpms = []
        for peakDim in peakDims:
          dataDimRef = peakDim.dataDimRef
          dimIsotopes.append(dataDimRef.expDimRef.isotopeCodes)
          dimPpms.append(peakDim.realValue)
          
          boundDim = boundDims.get(dataDimRef.dataDim)
          if boundDim:
            dimBound.append(boundDim.dim-1)
          else:
            dimBound.append(None)
            
          if replace:
            for contrib in peakDim.peakDimContribs:
              for peakContrib in contrib.peakContribs:
                if len(peakContrib.peakDimContribs) < 2:
                  peakContrib.delete()

              contrib.delete()
        
        for item in constraint.items:
          fResonances = list(item.resonances)
          resonance0, resonance1 = [fr.resonance for fr in fResonances]
 
          if not (resonance0 and resonance1):
            continue
 
          shift0 = resonance0.findFirstShift(parentList=shiftList)
          shift1 = resonance1.findFirstShift(parentList=shiftList)
          
          if not (shift0 and shift1):
            continue
 
          isotope0, isotope1 = [fr.isotopeCode for fr in fResonances]
          bound0 = resonance0.findFirstCovalentlyBound()
          bound1 = resonance1.findFirstCovalentlyBound()

          if bound0:
            shiftX0 = bound0.findFirstShift(parentList=shiftList)
          else:
            shiftX0 = None
          
          if bound1:
            shiftX1 = bound1.findFirstShift(parentList=shiftList)
          else:
            shiftX1 = None

          isotopeDims0 = []
          isotopeDims1 = []
          
          for i, isotopes in enumerate(dimIsotopes):
            boundDim = dimBound[i]
            
            if boundDim:
              if (isotope0 in isotopes) and bound0 and (bound0.isotopeCode in dimIsotopes[boundDim]):
                isotopeDims0.append(i)
              if (isotope1 in isotopes) and bound1 and (bound1.isotopeCode in dimIsotopes[boundDim]):
                isotopeDims1.append(i)
            
            else:
              if isotope0 in isotopes:
                isotopeDims0.append(i)
              if isotope1 in isotopes:
                isotopeDims1.append(i)
 
          assign = []
          for i in isotopeDims0:
            iX = dimBound[i]
            delta0 = abs(dimPpms[i]-shift0.value)
            
            if (iX is not None) and (shiftX0 is not None):
              delta0 += isotopeWeights[bound0.isotopeCode] * abs(dimPpms[iX]-shiftX0.value)
              delta0 /= 2.0
            
            for j in isotopeDims1:
              if i == j:
                continue
          
              jX = dimBound[j]
              delta1 = abs(dimPpms[j]-shift1.value)
              
              if (jX is not None) and (shiftX1 is not None):
                delta1 += isotopeWeights[bound1.isotopeCode] * abs(dimPpms[jX]-shiftX1.value)
                delta1 /= 2.0
  
              score = delta0 * delta0 + delta1 * delta1
              assign.append( (score, i, j, iX, jX) )  

          assign.sort()
          
          if assign:
            peakContribs = [peak.newPeakContrib(),]
          
            score, i, j, iX, jX = assign[0]
          
            contrib = peakDims[i].findFirstPeakDimContrib(resonance=resonance0)
            assignResToDim(peakDims[i], resonance0, contrib=contrib,
                           peakContribs=peakContribs)
                           
            contrib = peakDims[j].findFirstPeakDimContrib(resonance=resonance1)
            assignResToDim(peakDims[j], resonance1, contrib=contrib,
                           peakContribs=peakContribs)

            if (iX is not None) and bound0:
              contrib = peakDims[iX].findFirstPeakDimContrib(resonance=bound0)
              assignResToDim(peakDims[iX], bound0, contrib=contrib,
                             peakContribs=peakContribs)
            
            if (jX is not None) and bound1:
              contrib = peakDims[jX].findFirstPeakDimContrib(resonance=bound1)
              assignResToDim(peakDims[jX], bound1, contrib=contrib,
                             peakContribs=peakContribs)
             
  
def makePeaksFromConstraints(constraints, spectrum=None):
  """
  Predict the occurence of a set of NOESY peaks from a list of constraints. This
  function attempts to use sister lists to the input restraints linked peaks
  unless an overring target spectrum is passed in. This function is ideal for
  visualising the NOE disambiguation of procedures like ARIA and CYANA.
  
  .. describe:: Input
  
  List of NmrConstraint.DistanceConstraints, Nmr.DataSource
  
  .. describe:: Output
  
  List of Nmr.Peaks
  """

  isotopeWeights = {'1H':1.0,'2H':1.0,'15N':0.15,'13C':0.1}

  if spectrum:
    shiftList = spectrum.experiment.shiftList
    
    if not shiftList:
      msg  = 'Cannot make peaks from restraints.' 
      msg += ' Shift list is unset for input spectrum' 
      showWarning('Failure', msg)
      return []
    
  spectraDict = {}
  missedConstraints = []
  foundPeaks = []
  constraintPeaks = {}

  # find any existing linked peaks
  for constraint in constraints:
    if constraint.className == 'DistanceConstraint':
      peakContribs = constraint.peakContribs

      peakData = []
      if peakContribs:
        for peakContrib in peakContribs:
          peak0 = peakContrib.peak
          if not peak0: # Peak is now gone
            continue

          peakDims0 = peak0.sortedPeakDims()
          peakList0 = peak0.peakList
          spectrum0 = peakList0.dataSource
          spectraDict[spectrum0] = True
          
          positions = []
          for i in range(spectrum0.numDim):
            numAliasing = peakDims0[i].numAliasing
            positions.append([peakDims0[i].dataDimRef,peakDims0[i].position,numAliasing])
          
          intensities = []
          for peakIntensity in peak0.peakIntensities:
            intensityType = peakIntensity.intensityType
            method = peakIntensity.method
            value  = peakIntensity.value
            error  = peakIntensity.error
            intensities.append([intensityType,method,value,error])

          peakData.append([spectrum0,positions,intensities])
          
      constraintPeaks[constraint] = peakData
  
  if spectrum:
    spectra = [spectrum,]
  else:
    spectra = spectraDict.keys()
  
  peakLists = {}
  for spectrum0 in spectra:
    peakList = spectrum0.findFirstPeakList(isSimulated=True,name='ConstraintsPredict')
    
    if not peakList:
      peakList = spectrum0.newPeakList(isSimulated=True,name='ConstraintsPredict',
                                       details='Constraints Predicted')
      
      analysisPeakList = getAnalysisPeakList(peakList)
      analysisPeakList.symbolStyle = '+'
      analysisPeakList.symbolColor = '#FF0000'
      analysisPeakList.textColor = '#C00000'
    
    peakLists[spectrum0] = peakList
    
  if not peakLists:
    showWarning('Failure',
                'Cannot make peaks from restraints. No input spectra')
    return []
 
  # Setup isotopes and bound dims of target peakLists
  boundDims = {}  
  isotopeDims = {}
  for spectrum0 in spectra:
    for dimA, dimB in getOnebondDataDims(spectrum0):
      boundDims[dimA] = dimB
      boundDims[dimB] = dimA
    
    for dataDim in spectrum0.dataDims:
      for dataDimRef in dataDim.dataDimRefs:
        for isotopeCode in dataDimRef.expDimRef.isotopeCodes:
          if isotopeDims.get(isotopeCode) is None:
            isotopeDims[isotopeCode] = []
          
          if dataDim not in isotopeDims[isotopeCode]:
            isotopeDims[isotopeCode].append(dataDim)
  
  for constraint in constraintPeaks.keys():
    peakData = constraintPeaks[constraint]
    specPeakResonances = []

    # get all resonances, including for bound dims
    for spectrum0 in spectra:
      peakResonances = []
      
      for item in constraint.items:
        fResonances = list(item.resonances)
        resonance0, resonance1 = [fr.resonance for fr in fResonances]
        
        if not (resonance0 and resonance1):
          continue
    
        isotope0, isotope1 = [fr.isotopeCode for fr in fResonances]
        bound0   = resonance0.covalentlyBound
        bound1   = resonance1.covalentlyBound
        
        for dataDim0 in isotopeDims.get(isotope0, []):
          dim0 = dataDim0.dim
          boundDim0 = boundDims.get(dataDim0)
  
          for dataDim1 in isotopeDims.get(isotope1, []):
            if dataDim1.dataSource is not dataDim0.dataSource:
              continue
 
            dim1 = dataDim1.dim

            if dim0 != dim1:
              boundDim1 = boundDims.get(dataDim1)
               
              resonances0 = {dim0:resonance0,dim1:resonance1}
 
              if boundDim0 and boundDim1:
                for resonance2 in bound0:
                  for resonance3 in bound1:
                    resonances = resonances0.copy()
                    resonances[boundDim0.dim] = resonance2
                    resonances[boundDim1.dim] = resonance3
                    peakResonances.append(resonances)
 
              elif boundDim0:
                for resonance2 in bound0:
                  resonances = resonances0.copy()
                  resonances[boundDim0.dim] = resonance2
                  peakResonances.append(resonances)
  
              elif boundDim1:
                for resonance3 in bound1:
                  resonances = resonances0.copy()
                  resonances[boundDim1.dim] = resonance3
                  peakResonances.append(resonances)
  
              else:
                peakResonances.append(resonances0)
                
      specPeakResonances.append(peakResonances)
 
    if not peakData: # No peak contribs
      for spectrum0 in spectra:
        peakResonances = specPeakResonances.pop(0)
        if not peakResonances:
          missedConstraints.append(constraint)
          continue 
         
        # make peaks using resonance positions only
        peaks = []
        peakList  = peakLists[spectrum0]
        shiftList = spectrum0.experiment.shiftList
        dataDims  = spectrum0.sortedDataDims()
        #numDim    = spectrum0.numDim
        dims      = [dataDim.dim for dataDim in dataDims]
        specRanges = []
        isotopes = []

        for dataDim in dataDims:
          minFreq = None
          maxFreq = None
          isotope = []
 
          for dataDimRef in dataDim.dataDimRefs:
            expDimRef = dataDimRef.expDimRef
            valueMin  = expDimRef.minAliasedFreq or pnt2ppm(1,dataDimRef)
            valueMax  = expDimRef.maxAliasedFreq or pnt2ppm(dataDim.numPoints,dataDimRef)

            if (minFreq is None) or (valueMin < minFreq):
              minFreq = valueMin

            if (maxFreq is None) or (valueMax < maxFreq):
              maxFreq = valueMax
 
            isotope.extend(expDimRef.isotopeCodes)
 
          specRanges.append((minFreq,maxFreq))
          isotopes.append(isotope)
 
        for resonances in peakResonances:
          positions = []
 
          for dim in dims:
            resonance = resonances.get(dim)
            position  = None
 
            if resonance:
              if resonance.isotopeCode not in isotopes[dim-1]:
                break
            
              shift = resonance.findFirstShift(parentList=shiftList)

              if shift:
                value = shift.value
 
                if specRanges[dim-1][0] < value < specRanges[dim-1][1]:
                  position = value
 
            positions.append(position)
          
          else:
            if None not in positions:
              peak = pickPeak(peakList, positions, unit='ppm', doFit=True)
              for peakDim in peak.sortedPeakDims():
                resonance = resonances.get(peakDim.dim)
                if resonance:
                  assignResToDim(peakDim, resonance)
 
              peaks.append(peak)
 
            
        foundPeaks.extend(peaks)
          
    else:
      for spectrum0 in spectra:
        peakResonances = specPeakResonances.pop(0)
        
        if not peakResonances:
          missedConstraints.append(constraint)
          continue 
        
        peaks = []
        for spectrum1,positions,intensities in peakData:
          if spectrum1 is not spectrum0:
            continue
          peakList = peakLists[spectrum0]
          peak     = peakList.newPeak()
          peaks.append(peak)
 
          for peakDim in peak.sortedPeakDims():
            peakDim.dataDimRef, peakDim.position, peakDim.numAliasing = positions[peakDim.dim-1]
 
          for intensityType, method, value, error in intensities:
            peak.newPeakIntensity(method=method,intensityType=intensityType,
                                  value=value,error=error)

          shiftList = spectrum0.experiment.shiftList
          #bestScore = None
          #bestResonances = None
 
          for resonances in peakResonances:
            assignments = {}
          
            for peakDim in peak.peakDims:
              resonance = resonances.get(peakDim.dim)
 
              if resonance:
                if resonance.isotopeCode not in peakDim.dataDimRef.expDimRef.isotopeCodes:
                  continue
              
                shift = resonance.findFirstShift(parentList=shiftList)
 
                if shift:
                  delta = abs(shift.value-peakDim.value)*isotopeWeights.get(resonance.isotopeCode, 0.1)
                  if delta > 1.0: # Ought not be used: pathalogical case
                    continue
 
                  assignments[peakDim] = assignments.get(peakDim, []) + [resonance,]
                    
            else:
              for peakDim in assignments:
                boundDim = boundDims.get(peakDim.dataDim)

                for resonance in assignments[peakDim]:
                  
                  if boundDim:
                    peakDim2 = peak.findFirstPeakDim(dim=boundDim.dim)
                    
                    for resonance2 in assignments.get(peakDim2, []):
                      if areResonancesBound(resonance, resonance2):
                        break
                        
                    else:
                      continue    

                    assignResToDim(peakDim, resonance)
                  
                  else:
                    assignResToDim(peakDim, resonance)
 
        foundPeaks.extend(peaks)
 
  if missedConstraints:
    msg = 'Peaks could not be made for %d constraints, due to missing resonances.'
    showWarning('Warning', msg % len(missedConstraints))

  return foundPeaks

def copyFixedResonance(fixedResonance, nmrConstraintStore):
  """
  Make an equivalent fixed resonance in the input constraint store
  
  .. describe:: Input
  
  NmrConstraint.FixedResonance, Nmr.NmrConstraintStore
  
  .. describe:: Output
  
  NmrConstraint.FixedResonance
  """
  
  tags = ('resonanceSerial', 'name', 'details', 'isotopeCode')
  resonanceMap = {}
  
  cSet = nmrConstraintStore
  
  if fixedResonance.resonanceSerial:
    fixedResonance2 = cSet.findFirstFixedResonance(resonanceSerial=fixedResonance.resonanceSerial)

    if fixedResonance2:
      return fixedResonance2
    
    resonanceSet = fixedResonance.resonanceSet

  if resonanceSet:
    atomSets2 = [getFixedAtomSet(cSet, ass.atoms) for ass in resonanceSet.atomSets]
    resonanceSet2 = cSet.findFirstFixedResonanceSet(atomSets=atomSets2)
    
    if resonanceSet2:
      index = list(resonanceSet.resonances).index(fixedResonance)
      fixedResonance2 = list(fixedResonance2.resonances)[index]
    
    else:
      kw = {}
      for tag in tags:
        kw[tag] = getattr(fixedResonance,tag)
      fixedResonance2 =  cSet.newFixedResonance(**kw)
      #fixedResonances2 = [fixedResonance2]
      resonanceMap[fixedResonance] = fixedResonance2

      for fixedResonance1 in resonanceSet.resonances:
        if fixedResonance1 is not fixedResonance:
          kw = {}
          for tag in tags:
            kw[tag] = getattr(fixedResonance1,tag)
          fixedResonance3 = cSet.newFixedResonance(**kw)
          #fixedResonances2.append(fixedResonance3)
          resonanceMap[fixedResonance1] = fixedResonance3

      resonanceSet2 = cSet.newFixedResonanceSet(atomSets=atomSets2,
                                                resonances=resonanceMap.values)

  else:
    kw = {}
    for tag in tags:
      kw[tag] = getattr(fixedResonance,tag)
    fixedResonance2 = cSet.newFixedResonance(**kw)
    resonanceMap[fixedResonance] = fixedResonance2
  
  for fr1, fr2 in resonanceMap.items():
    transferCovalentlyBound(fr1, fr2)
  
  return fixedResonance2


def getFixedAtomSet(nmrConstraintStore, atoms):
  """
  Finds or creates a fixed set of atoms that is used in an NMR
  constraint head object (equivalent to one NmrConstraint file).
  Creating fixed atom sets allows assignments to change but
  old constraints to be preserved.

  .. describe:: Input
  
  Nmr.NmrConstraintStore, List of MolSystem.Atoms
  
  .. describe:: Output
  
  NmrConstraint.FixedAtomSet
  """

  atom = list(atoms)[0]

  if not hasattr(nmrConstraintStore, 'quickAtomSets'):
    nmrConstraintStore.quickAtomSets = {}
 
  fixedAtomSet = nmrConstraintStore.quickAtomSets.get(atoms)

  if not fixedAtomSet:
    fixedAtomSet = nmrConstraintStore.findFirstFixedAtomSet(atoms=atoms)

  if not fixedAtomSet:
    atomSet = atom.atomSet
    if atomSet:
      fixedAtomSet = nmrConstraintStore.newFixedAtomSet(atoms=atomSet.atoms, name=atomSet.name)
    
    else:
      fixedAtomSet2 = atom.findFirstFixedAtomSet(atoms=atoms)
      if fixedAtomSet2:
        fixedAtomSet = nmrConstraintStore.newFixedAtomSet(atoms=fixedAtomSet2.atoms,name=fixedAtomSet2.name)

  nmrConstraintStore.quickAtomSets[atoms] = fixedAtomSet
  
  return fixedAtomSet
       

def getConstraintAtoms(constraint):
  """
  Get the atoms that may be assigned to the constrained resonances
  The outer most list is due to restraint ambiguity, the middle list
  is the list for different resonances and the inner list
  is for equivalent atoms.

  .. describe:: Input

  NmrConstraint.AbstractConstraint

  .. describe:: Output

  List of List of List of MolSystem.Atoms
  """

  atoms = []
  fixedResonances = []
  className = constraint.className
  
  if className == 'DihedralConstraint':
    fixedResonances.append(constraint.resonances)

  elif className in ('ChemShiftConstraint', 'CsaConstraint'):
    fixedResonances.append([constraint.resonance,])

  else:
    for item in constraint.items:
      fixedResonances.append(item.resonances)

  for fixedResonanceList in fixedResonances:
    atomList = []
    for fixedResonance in fixedResonanceList:
      fixedResonanceSet = fixedResonance.resonanceSet
    
      if fixedResonanceSet:
        equivAtoms = {}
        
        for fixedAtomSet in fixedResonanceSet.atomSets:
          for atom in fixedAtomSet.atoms:
            equivAtoms[atom] = True
            
        atomList.append(equivAtoms.keys())    

    if len(atomList) == len(fixedResonanceList):
      atoms.append(atomList)

  return atoms

def getJCouplingsFromConstraint(jCouplingConstraint, jCouplingList):
  """
  The the J couplings associated with a J coupling constraint

  .. describe:: Input
  
  NmrConstraint.JCouplingConstraint, Nmr.JCouplingList
  
  .. describe:: Output
  
  List of Nmr.JCouplings
  """
  #nmrProject = jCouplingConstraint.parentList.nmrConstraintStore.nmrProject

  jCouplings = []
  
  for item in jCouplingConstraint.items:
    resonances = [fr.resonance for fr in item.resonances]
    
    if None not in resonances:
      jCoupling = jCouplingList.findFirstMeasurement(resonances=resonances)
      
      if jCoupling:
        jCouplings.append(jCoupling)

  return jCouplings

def mergeDuplicateConstraints(constraintList):
  """
  Merge redundant constraints in a constraint list (constraining the same
  resonances)
  
  .. describe:: Input
  
  NmrConstraint.ConstraintList
  
  .. describe:: Output
  
  NmrConstraint.ConstraintList
  """

  merged = []
  resonanceDict = {}
  
  className = constraintList.className
  if className == 'DihedralConstraintList':

    # Can this ever occur? Uniqueness of resonances checked?

    for constraint in list(constraintList.constraints):
      key = ' '.join(['%d' % r.serial for r in constraint.resonances])   
      if resonanceDict.get(key) is None:
        resonanceDict[key] = []
      resonanceDict[key].append(constraint)

    for key in resonanceDict.keys():
      constraints = resonanceDict[key]

      if len(constraints) > 1:
        constraint0 = constraints[0]
        merged.append(constraint0)

        for constraint in constraints[1:]:
          for item in constraint.items:
            constraint0.newDihedralConstraintItem(targetValue=item.targetValue,
                                                  upperLimit=item.upperLimit,
                                                  lowerLimit=item.lowerLimit,
                                                  error=item.error)
          constraint.delete()

  else:
    if className == 'ChemShiftConstraintList':
  
      for constraint in constraintList.constraints:
        key = constraint.resonance.serial
        if resonanceDict.get(key) is None:
          resonanceDict[key] = []
        resonanceDict[key].append(constraint)

    else:

      for constraint in list(constraintList.constraints):
        serialList = []
      
        for item in constraint.items:
          subSerials = [r.serial for r in item.resonances]
          subSerials.sort()
          serialList.append(subSerials)

        serialList.sort()

        serials = []
        for subSerials in serialList:
          for serial in subSerials:
            serials.append('%d' % serial) 

        key = ' '.join([s for s in serials])
        if resonanceDict.get(key) is None:
          resonanceDict[key] = []
        resonanceDict[key].append(constraint)
    
    for key in resonanceDict.keys():
      constraints = resonanceDict[key]

      if len(constraints) > 1:
        
        constraint0 = constraints[0]
        upperLimit  = constraint0.upperLimit
        lowerLimit  = constraint0.lowerLimit
        targetValue = constraint0.targetValue
        contribSerials = []
        merged.append(constraint0)
        
        for constraint in constraints[1:]:
       
          if constraint.upperLimit < upperLimit:
            upperLimit = constraint.upperLimit

          if constraint.lowerLimit < lowerLimit:
            lowerLimit = constraint.lowerLimit

          if constraint.targetValue < targetValue:
            targetValue = constraint.targetValue
          
          for contrib in constraint.peakContribs:
            serials = [contrib.experimentSerial, contrib.dataSourceSerial,
                       contrib.peakListSerial, contrib.peakSerial, contrib.peak]
            contribSerials.append(serials)
          
          constraint.delete()

        constraint0.targetValue = targetValue
        constraint0.upperLimit  = upperLimit
        constraint0.lowerLimit  = lowerLimit

        for e,s,pl,p, peak in contribSerials:
          if constraint0.findFirstPeakContrib(peak=peak) is None:
            constraint0.newConstraintPeakContrib(experimentSerial=e,
                                                 dataSourceSerial=s,
                                                 peakListSerial=pl,
                                                 peakSerial=p)

  return merged

def makeStructureGeneration(structure, constraintSet=None, generationType='denovo'):
  """
  Make a structure generation object to group structures together
  
  .. describe:: Input
  
  MolSystem.StructureEnsemble, Nmr.NmrConstraintStore
  
  .. describe:: Output
  
  Nmr.StructureGeneration
  """
  
  if constraintSet:
    nmrProject = constraintSet.nmrProject
  else:
    nmrProject = structure.root.currentNmrProject

  structureGeneration = nmrProject.newStructureGeneration(generationType=generationType,
                                                          structureEnsemble=structure)

  if constraintSet:
    structureGeneration.nmrConstraintStore = constraintSet

  return structureGeneration

def mergeConstraintLists(constraintListA, constraintListB):
  """
  Merge two constraint lists from a constraint set into one.
  
  .. describe:: Input
  
  NmrConstraint.ConstraintList, NmrConstraint.ConstraintList
  
  .. describe:: Output
  
  NmrConstraint.ConstraintList
  """
  
  if constraintListA.className != constraintListB.className:
    showWarning('Failure','Cannot merge constraint lists of different types')
    return

  if constraintListA.nmrConstraintStore is not constraintListB.nmrConstraintStore:
    showWarning('Failure','Cannot merge constraint lists from different parent sets')
    return

  for constraint in constraintListB.constraints:
    copySubTree(constraint, constraintListA,
                objectMap=None,
                maySkipCrosslinks=True)
    constraint.delete()
    
  if constraintListA.name in ('Ambig','Unambig'):
    if constraintListB.name in ('Ambig','Unambig'):
      if constraintListB.name != constraintListA.name:
        constraintListA.name = 'Merged'

  for serial in constraintListB.experimentSerials:
    if serial not in constraintListA.experimentSerials:
      constraintListA.addExperimentSerial(serial)

  constraintListB.delete()
  
  return constraintListA


def splitConstraintListAmbigUnambig(constraintList):
  """
  Spit a constraint list into two lists according to whether
  constraints are ambigous or unambiguous.

  .. describe:: Input
  
  NmrConstraint.ConstraintList
  
  .. describe:: Output
  
  NmrConstraint.ConstraintList, NmrConstraint.ConstraintList

  """

  className = constraintList.className
  if className == 'ChemShiftConstraintList':
    msg = 'ChemShiftConstraintList class does not have ambiguous constraints'
    showWarning('Split Ambig/Unambig Warning', msg)
    return

  from ccp.api.nmr import NmrConstraint

  ListClass = getattr(NmrConstraint, className)

  ambigList         = []
  constraintListB   = None
  nmrConstraintStore = constraintList.nmrConstraintStore
  
  if className == 'DihedralConstraintList':
    # ambiguity is at level of ranges not resonances
    for constraint in constraintList.constraints:
      if len(constraint.items) > 1:
        ambigList.append(constraint)

  else:
    for constraint in constraintList.constraints:
      items = constraint.items
      if len(items) > 1:
        prochiralA = True
        prochiralB = True
        resDictA = {}
        resDictB = {}
        
        for item in items:
          resonances = tuple(item.resonances)
          resDictA[ resonances[0] ] = None
          resDictB[ resonances[1] ] = None
      
        resonancesA = resDictA.keys()
        resonancesB = resDictB.keys()
        
        for resonance in resonancesA[1:]:
          if not areResonancesProchiral(resonancesA[0], resonance):
            prochiralA = False
            break

        if prochiralA:
          for resonance in resonancesB[1:]:
            if not areResonancesProchiral(resonancesB[0], resonance):
              prochiralB = False
              break
     
        if (not prochiralA) or (not prochiralB): # Still true if single resonance
          ambigList.append(constraint)
    

  if ambigList:
  
    unit = constraintList.unit
    experimentSerials  = list(constraintList.experimentSerials)
    measureListSerials = list(constraintList.measureListSerials)
    constraintListB    = ListClass(nmrConstraintStore, unit=unit,
                                   experimentSerials=experimentSerials,
                                   measureListSerials=measureListSerials)
                               
    for constraint in ambigList:
      moveConstraintToList(constraint, constraintListB)

    if (not constraintList.name) or (constraintList.name == 'Ambig'):
      constraintList.name = 'Unambig'
    
    constraintListB.name = 'Ambig'

  return constraintListB

def splitConstraintListViol(constraintList, violationList=None):
  """
  Spit a constraint list into two lists according to whether
  constraints show a violation (in th einput violation list).
  
  .. describe:: Input
  
  NmrConstraint.ConstraintList, NmrConstraint.ViolationList 
  
  .. describe:: Output
  
  NmrConstraint.ConstraintList, NmrConstraint.ConstraintList

  """

  violList  = []
  constraintListB  = None
  nmrConstraintStore = constraintList.nmrConstraintStore

  for constraint in constraintList.constraints:
  
    if violationList:
      violation = constraint.findFirstViolation(violationList=violationList)
    else:
      violation = constraint.findFirstViolation()
    
    if violation:
      violList.append(constraint)
    

  if violList:
    from ccp.api.nmr import NmrConstraint
    ListClass = getattr(NmrConstraint, constraintList.className)
 
    unit = constraintList.unit
    experimentSerials  = list(constraintList.experimentSerials)
    measureListSerials = list(constraintList.measureListSerials)
    constraintListB    = ListClass(nmrConstraintStore, unit=unit,
                                   experimentSerials=experimentSerials,
                                   measureListSerials=measureListSerials)
                               
    for constraint in violList:
      moveConstraintToList(constraint, constraintListB)
    
    constraintListB.name = 'Violated'

  return constraintListB

def moveConstraintToList(constraint, constraintList):
  """
  Move a constraint to a different constraint list and take care of any
violations.
  
  .. describe:: Input
  
  NmrConstraint.Constraint,  NmrConstraint.ConstraintList
  
  .. describe:: Output
  
  NmrConstraint.Constraint

  """

  newConstraint = copySubTree(constraint, constraintList, objectMap=None, maySkipCrosslinks=True)
  
  for violation in constraint.violations:
    violationList = violation.violationList
    violationList.newViolation(violation=violation.violation,
                               calcValue=violation.calcValue,
                               calcValueError=violation.calcValueError,
                               method=violation.method,
                               fractionViolated=violation.fractionViolated,
                               constraint=newConstraint)
                               
  for violation in constraint.violations:
    violation.delete()
  
  constraint.delete()
  
  return newConstraint

def getStructureViolations(constraintList, structure, violationList=None):
  """
  Calculate the violated restraints in a restraint list given in input
  structure. Currently limited to distance and dihedral restraints. Option to
  add violations to a speficied violation list (otherwise a new one will be
  made) 
  
  .. describe:: Input
  
  NmrConstraint.ConstraintList, MolStructure.StructureEnsemble
  
  .. describe:: Output
  
  NmrConstraint.ViolationList
  """

  from ccpnmr.analysis.core.StructureBasic import getAtomSetsDistance, getAtomSetsDihedral
  from math import cos, sin, atan2

  negSixth = -1/6.0

  if constraintList.className not in ('DistanceConstraintList',
                                      'DihedralConstraintList',
                                      'HBondConstraintList'):
    return

  # Violation only if all items violated

  constraintSet = constraintList.nmrConstraintStore
 
  violDict = {}
  if violationList:
    if violationList.nmrConstraintStore is not constraintSet:
      showWarning('Warning','Violation list does not match constraint set')
      return

    for violation in violationList.violations:
      violDict[violation.constraint] = violation

  models = structure.models

  if constraintList.className != 'DihedralConstraintList':
    
    # resolve prochirals
    bestAtomSet = {}
    for constraint in constraintList.constraints:
      targetValue = constraint.targetValue

      if targetValue is None:
        targetValue = (constraint.lowerLimit + constraint.upperLimit)/2.0
        
      for item in constraint.items:
        resonance1, resonance2 = item.resonances
        
        if bestAtomSet.get(resonance1) is None:
          bestAtomSet[resonance1] = {}
        if bestAtomSet.get(resonance2) is None:
          bestAtomSet[resonance2] = {}
        
        resonanceSet1 = resonance1.resonanceSet
        resonanceSet2 = resonance2.resonanceSet
        
        if resonanceSet1:
          if resonanceSet2:
            atomSets1 = resonanceSet1.atomSets
            atomSets2 = resonanceSet2.atomSets
            
            delta  = None
            best1V = None
            best2V = None
            deltaV = None
            for atomSet1 in atomSets1:
              for atomSet2 in atomSets2:
                # Average over whole ensemble
                dist  = getAtomSetsDistance([atomSet1,], [atomSet2,],
                                            structure, method='noe')
                if dist is None:
                  continue
                
                delta = abs(targetValue - dist)
                   
                if dist < constraint.lowerLimit:
                  delta += constraint.lowerLimit - dist
               
                elif dist > constraint.upperLimit:
                  delta += dist - constraint.upperLimit
                                       
                if (deltaV is None) or (delta < deltaV):
                  best1V = atomSet1
                  best2V = atomSet2
                  deltaV = delta
              
            if deltaV is not None:
              atomSet1 = best1V
              atomSet2 = best2V
 
            bestAtomSet[resonance1][atomSet1] = bestAtomSet[resonance1].get(atomSet1, 0) + 1
            bestAtomSet[resonance2][atomSet2] = bestAtomSet[resonance2].get(atomSet2, 0) + 1
    
    resolvedAtomSets = {}
    for resonance in bestAtomSet.keys():
      atomSet   = None
      bestCount = None
      
      for atomSet1 in bestAtomSet[resonance].keys():
        if (atomSet is None) or (bestAtomSet[resonance][atomSet1] > bestCount):
          bestCount = bestAtomSet[resonance][atomSet1]
          atomSet = atomSet1
      
      resolvedAtomSets[resonance] = atomSet
      
      # Make sure prochirals don't both go the same way
      # Can happen in unstructured proteins
      resonanceSet = resonance.resonanceSet
      if resonanceSet:
        resonances = list(resonanceSet.resonances)
        resonances.remove(resonance)
        atomSets = list(resonanceSet.atomSets)
        atomSets.remove(atomSet)
        
        if resonances and atomSets:
          resolvedAtomSets[resonances[0]] = atomSets[0]
        
             
    for constraint in constraintList.constraints:
      targetValue = constraint.targetValue

      if targetValue is None:
        targetValue = (constraint.lowerLimit + constraint.upperLimit)/2.0
            
      totalViols  = 0
      violAmounts = []
      noeSum = 0.0
      noeCount = 0
      
      for model in models:
        dists = []
        for item in constraint.items:
          resonance1, resonance2 = item.resonances
          if resonance1.resonanceSet:
            if resonance2.resonanceSet:
              atomSet1 = resolvedAtomSets[resonance1]
              atomSet2 = resolvedAtomSets[resonance2]
              
              dist  = getAtomSetsDistance([atomSet1,], [atomSet2,], structure,
                                          model=model, method='noe')
                                  
              if not dist:
                continue
 
              dists.append(dist)
                
        if dists:
          noe = 0.0
          for dist in dists:
            noe += dist ** -6.0
          # c.f. similar logic in StructureBasic.getAtomSetsDistance
          noeSum += sqrt(noe/len(dists))
          noeCount += 1
 
          dist2  = noe ** negSixth
          deltaL = constraint.lowerLimit - dist2
          deltaU = dist2 - constraint.upperLimit
          
          if deltaL > 0:
            totalViols += 1 
            violAmounts.append(deltaL)
          
          elif deltaU > 0:
            totalViols += 1 
            violAmounts.append(deltaU)
            
      if totalViols:
        if noeCount > 0:
          noeSum /= noeCount
        calcValue = noeSum ** (-1/3.0)
        N = float(len(violAmounts))
      
        if violationList is None:
          violationList = constraintSet.newViolationList(molStructures=models)
          violDict = {}
 
        aveViolAmount = 0.0
        for violAmount in violAmounts:
          aveViolAmount += violAmount
 
        aveViolAmount = aveViolAmount/N
 
        sumSq = 0.0
        for violAmount in violAmounts:
          diff  = violAmount - aveViolAmount
          sumSq += diff * diff
            
        sd = sqrt(sumSq/N)
 
        fracViolated = totalViols/float(len(models))
        violObj = violDict.get(constraint)
        
        if violObj:
          violObj.violation = aveViolAmount
          violObj.calcValue = calcValue
          violObj.calcValueError = sd
          violObj.fractionViolated=fracViolated
 
        else:
          violationList.newViolation(violation=aveViolAmount, calcValue=calcValue,
                                     calcValueError=sd, fractionViolated=fracViolated,
                                     constraint=constraint)
 

  else: 

    twoPi = 6.2831853071796
    if constraintList.unit in  ('degrees', None):
      inDegrees = True
    else:
      inDegrees = False 

    if inDegrees:
      oneTurn = 360.0
    else:
      oneTurn = twoPi

    for constraint in constraintList.constraints:
      skip = False
      for resonance in constraint.resonances:
        if not resonance.resonanceSet:
          skip = True
          break
      
      if skip:
        continue
              
      atomSets = []
      for resonance in constraint.resonances:
        if resonance.resonanceSet:
          atomSets.append(resonance.resonanceSet.atomSets)

      if len(atomSets) != 4:
        continue

      mean   = 0.0
      delta  = 0.0
      angles = []
      totalViols = 0
      for model in models:
        angle = getAtomSetsDihedral(atomSets, structure,
                                    model=model, inDegrees=inDegrees)
                
        if angle is None:
          continue
          
        angles.append(angle)
        nViols = 0
        minDelta = None
        
        for item in constraint.items:
          angleValue = angle % oneTurn
          upperLimit = item.upperLimit % oneTurn
          lowerLimit = item.lowerLimit % oneTurn
          
          diff1 = angleDifference(angleValue, lowerLimit, oneTurn)
          diff2 = angleDifference(angleValue, upperLimit, oneTurn)
          
          while lowerLimit > upperLimit:
            lowerLimit -= oneTurn

          while angleValue > upperLimit:
            angleValue -= oneTurn
                    
          if angleValue < lowerLimit:
            nViols += 1
            violAngle = min(diff1, diff2) 
            
            if (minDelta is None) or violAngle < minDelta:
              minDelta = violAngle

        if nViols and nViols == len(constraint.items): # if every range item is violated
          totalViols += 1
          delta += minDelta


      if totalViols > 0:
        if violationList is None:
          violationList = constraintSet.newViolationList(molStructures=models)
          violDict = {}

        violAmount = delta/float(totalViols)

        meanCos = 0.0
        meanSin = 0.0
        N = float(len(angles))
        for angle in angles:
          if inDegrees:
            a = twoPi * angle/360.0
          else:
            a = angle

          meanCos += cos(a)
          meanSin += sin(a) 

        meanCos /= N
        meanSin /= N

        mean = atan2(meanSin,meanCos)
        if inDegrees:
          mean *= 360/twoPi

        sigma = 0.0
        for angle in angles:
          delta = angleDifference(mean,angle,oneTurn)
          sigma += (delta*delta)

        sd = sqrt(sigma/N)
        fracViolated = totalViols/float(len(models))
        violObj = violDict.get(constraint)
        
        if violObj:
          violObj.violation = violAmount
          violObj.calcValue = mean 
          violObj.calcValueError = sd
          violObj.fractionViolated=fracViolated 
          
        else:
          violationList.newViolation(violation=violAmount, calcValue=mean,
                                     calcValueError=sd, fractionViolated=fracViolated,
                                     constraint=constraint)

  return violationList

def angleDifference(a1, a2, oneTurn):

  delta = abs(a1- a2) % oneTurn
  
  if delta > (oneTurn/2.0):
    delta = oneTurn - delta

  return delta

def exportAriaTbl(constraints, fileName):
  """
  Exports a constraint list to a file as an ARIA 1.x readable .tbl file (similar
  to CNS format). For proper CNS format use CcpNmr FormmatConverter.
  NBNB will expand restraisnt (A<->C) or B<->D) to (A or B)<->(C or D) 
  
  .. describe:: Input
  
  List of NmrConstraint.DistanceConstraints, String (export file name)
  
  .. describe:: Output
  
  None
  """

  resonanceList = []
  useSegIds = False
  chains = set()

  fileHandle = open(fileName, 'w')
  for c in constraints:
    if not c.peaks:
      continue
    
    p = list(c.peaks)[0]
    shiftList = p.peakList.dataSource.experiment.shiftList
    
    for useIsotopeCode in ('1H', '13C', '15N'):
      hd = findSpectrumDimsByIsotope(p.peakList.dataSource, useIsotopeCode)
      if len(hd) == 2:
        break
    else:
      continue
    
    #hd = findSpectrumDimsByIsotope(p.peakList.dataSource, '1H')
    #if len(hd) != 2:
    
    #  hd = findSpectrumDimsByIsotope(p.peakList.dataSource, '13C')
    #  if len(hd) != 2:
    
    #    hd = findSpectrumDimsByIsotope(p.peakList.dataSource, '15N')
    #    if len(hd) != 2:
    #       continue
 
    ppms = []
    possiblePpms = []
    peakDims = p.sortedPeakDims()
    for n in range(2):
      peakDim = peakDims[hd[n]]
      value = peakDim.value
      ppms.append(value)
      dataDimRef = peakDim.dataDimRef
      if dataDimRef:
        expDimRef = dataDimRef.expDimRef
        unit = expDimRef.unit
        if not unit or unit == 'None':
          unit = 'ppm'
        valRange = [unit_converter[('point', unit)](1.0, dataDimRef),
                    unit_converter[('point', unit)](1.0+dataDimRef.dataDim.numPoints, dataDimRef) ]
        valRange.sort()
        if expDimRef.minAliasedFreq is None:
          minFreq = valRange[0]
        else:
          minFreq = expDimRef.minAliasedFreq
        if expDimRef.maxAliasedFreq is None:
          maxFreq = valRange[1]
        else:
          maxFreq = expDimRef.maxAliasedFreq
        possiblePpms.append(getAliasedPeakDimPositions(peakDim, [(minFreq, maxFreq)], returnPpms=True))
      else:
        possiblePpms.append([value])
    
    possiblePpms1, possiblePpms2 = possiblePpms
    ppm1, ppm2 = ppms

    #pd = p.sortedPeakDims()[hd[0]]
    #ppm1 = pd.value
    #pd = p.sortedPeakDims()[hd[1]]
    #ppm2 = pd.value
    
    # 16 Sep 2014: wb104: removed the resonances0 and resonances1 sets and replaced with a dictionary 
    # the problem is that in a given constraint the resonances in the items don't necessarily pair
    # with resonances in all other items
    # it's very messy to deal with this but hopefully the below does it correctly
    # resonances01Dict has as key a resonance in an item and the value is all other resonances that
    # get paired with it in other items, as a set
    # then we use resonances10Dict to determine which resonances have the same set
    
    #resonances0 = set()
    #resonances1 = set()
    resonances01Dict = {}
    for item in c.sortedItems():
      fixedResonances = list(item.resonances)
      
      #rA, rB = [fr.resonance for fr in fixedResonances]
      
      # Bug fix: If experiment is a H,C,C 3D subset of a 4D NOESY, 
      # and resonances are given as protons (they will be)
      # you must get the shift of the bound carbon for the item matching to work.
      # Rasmus 12/3/13
      resl = [fr.resonance for fr in fixedResonances]
      
      for ii,rr in enumerate(resl):
        if useIsotopeCode != rr.isotopeCode:
          newRes = getOnebondResonance(rr, useIsotopeCode)
          if newRes:
            resl[ii] = newRes
      rA,rB = resl
      # end bug fix
      
      if not (rA and rB):
        continue 
       
      shiftA = rA.findFirstShift(parentList=shiftList)
      shiftB = rB.findFirstShift(parentList=shiftList)
      
      if not (shiftA and shiftB):
        continue 
      
      valueA = shiftA.value
      valueB = shiftB.value
      
      # wb104: 12 Sep 2014
      # problem with below is that it ignores aliasing so might pair up incorrectly
      # replaced it with _bestMatch formalism

      #deltaA1 = valueA-ppm1
      #deltaB2 = valueB-ppm2
      
      #deltaA2 = valueA-ppm2
      #deltaB1 = valueB-ppm1
     
      #distAB = (deltaA1*deltaA1) + (deltaB2*deltaB2)
      #distBA = (deltaA2*deltaA2) + (deltaB1*deltaB1)
      
      def _bestMatch(value1, value2):
        bestDist = None
        for possiblePpm1 in possiblePpms1:
          delta1 = value1 - possiblePpm1
          for possiblePpm2 in possiblePpms2:
            delta2 = value2 - possiblePpm2
            dist = delta1*delta1 + delta2*delta2
            if bestDist is None or dist < bestDist:
              bestDist = dist
        return bestDist

      distAB = _bestMatch(valueA, valueB)
      distBA = _bestMatch(valueB, valueA)

      if distAB is None or distBA is None or distAB > distBA:
        #resonances0.add(fixedResonances[0])
        #resonances1.add(fixedResonances[1])
        resonances01Dict.setdefault(fixedResonances[0], set()).add(fixedResonances[1])
      else:
        #resonances0.add(fixedResonances[1])
        #resonances1.add(fixedResonances[0])
        resonances01Dict.setdefault(fixedResonances[1], set()).add(fixedResonances[0])
   
    # need to check now which resonances0 match with which resonances1
    # need separate assignment for each set of pairings which match
    
    resonances10Dict = {}
    for resonance0 in resonances01Dict:
      resonances1 = frozenset(resonances01Dict[resonance0])
      resonances10Dict.setdefault(resonances1, set()).add(resonance0)
    
    for resonances1 in resonances10Dict:
      resonances0 = resonances10Dict[resonances1]
      
      resonances = []
      for resonance in resonances1:
        seqCode = 0
        resonanceSet = resonance.resonanceSet
        if resonanceSet:
          residue = resonanceSet.findFirstAtomSet().findFirstAtom().residue
          seqCode = residue.seqCode
          chains.add(residue.chain)
        resonances.append([seqCode,resonance])

      resonances.sort()
      resonancesA = [x[1] for x in resonances]
    
      resonances = []
      for resonance in resonances0:
        seqCode = 0
        resonanceSet = resonance.resonanceSet
        if resonanceSet:
          residue = resonanceSet.findFirstAtomSet().findFirstAtom().residue
          seqCode = residue.seqCode
          chains.add(residue.chain)
        resonances.append([seqCode,resonance])
      resonances.sort()
      resonancesB = [x[1] for x in resonances]
    
      resonanceList.append((c, resonancesA, resonancesB, ppm1, ppm2))

  if len(chains) > 1:
    if showYesNo('Query','Use SEGIDs in output?'):
      useSegIds = True
    else:
      useSegIds = False

  for c, resonancesA, resonancesB, ppm1, ppm2 in resonanceList:
    
    for resonance in resonancesA:
      if resonance is resonancesA[0]:
        a = 'assign ( '
      else:
        a = '\t or '
        
      resonanceSet = resonance.resonanceSet
      if resonanceSet:
        atomSets = resonanceSet.sortedAtomSets()
        residue  = atomSets[0].findFirstAtom().residue
         
        if resonance is resonanceSet.sortedResonances()[0]:
          name = atomSets[0].name
        else:
          name = atomSets[-1].name
          
        name = re.sub('\*','#',name)

        if useSegIds:
          b = 'segid %s and resid %d and name %s' % (residue.chain.code, residue.seqCode, name.upper())
        else:
          b = 'resid %d and name %s' % (residue.seqCode, name.upper())
      else:
        b = 'resid ? and name ?'
      fileHandle.write('%s%s\n' % (a,b))
    fileHandle.write('\t )\n')
    
    for resonance in resonancesB:
      if resonance is resonancesB[0]:
        a = '       ( '
      else:
        a = '\t or '
        
      resonanceSet = resonance.resonanceSet
      if resonance.resonanceSet:
        atomSets = resonanceSet.sortedAtomSets()
        residue  = atomSets[0].findFirstAtom().residue
         
        if resonance is resonanceSet.sortedResonances()[0]:
          name = atomSets[0].name
        else:
          name = atomSets[-1].name
          
        name = re.sub('\*','#',name)
          
        if useSegIds:
          b = 'segid %s and resid %d and name %s' % (residue.chain.code, residue.seqCode, name.upper())
        else:
          b = 'resid %d and name %s' % (residue.seqCode, name.upper() )
      else:
        b = 'resid ? and name ?'
      fileHandle.write('%s%s\n' % (a,b))
      
    deltaMin = abs(c.targetValue - c.lowerLimit)
    deltaMax = abs(c.targetValue - c.upperLimit)
    data = (c.targetValue,deltaMin,deltaMax,c.origData or 1.0,c.serial,ppm1,ppm2 )
    fileHandle.write('\t ) %.1f %.1f %.1f volume=%.3f peak=%d ppm1=%.3f ppm2=%.3f\n\n' % data)
  
  fileHandle.close()


def getConstraintStoreResonances(nmrConstraintStore, makeUnassignedResonances=True):
  """
  Makes any missing resonances for assigned fixedResonances in an
  NmrConstraintStore with the option to also make new resonances for unassigned
  fixedResonances. Returns all the resonances linkable from the constraints in
  an NmrConstraintStore object.
  
  .. describe:: Input
  
  Nmr.NmrConstraintStore, Boolean
  
  .. describe:: Output
  
  List of Nmr.Resonances
  """
  
  from ccpnmr.analysis.core.AssignmentBasic import initResonance
  
  nmrProject = nmrConstraintStore.topObject.nmrProject
  atomSetDict = {}
  for atomSet in nmrProject.atomSets:
    atomSetDict[frozenset(atomSet.atoms)] = atomSet

  resonances = []
  for fixedResonance in nmrConstraintStore.fixedResonances:
    if fixedResonance.resonanceSerial:
      resonances.append(fixedResonance.resonance)
    
    elif fixedResonance.resonanceSet:
      atomSets = []
      resonance = None
      resonanceSet = None
      
      for atomSet in fixedResonance.resonanceSet.atomSets:
        atoms = frozenset(atomSet.atoms)
        if atomSetDict.get(atoms) is not None:
          atomSets.append( atomSetDict[atoms] )
      
      for atomSet in atomSets:
        for resonanceSet0 in atomSet.resonanceSets:
          if resonanceSet0.atomSets == set(atomSets):
            resonanceSet = resonanceSet0
      
      if resonanceSet:
        fixedResonances = fixedResonance.resonanceSet.sortedResonances()
        i = fixedResonances.index(fixedResonance)
        if i < len(resonanceSet.resonances):
          resonance = resonanceSet.sortedResonances()[i]
      
        else:
          for resonance0 in resonanceSet.resonances:
            if resonance0.name == fixedResonance.name:
              resonance = resonance0
              break
          else:
            resonance = resonanceSet.findFirstResonance()
           
      if not resonance:
        resonance = nmrProject.newResonance(isotopeCode=fixedResonance.isotopeCode)
        resonance.setName(fixedResonance.name)
        assignAtomsToRes(atomSets, resonance)
      
      fixedResonance.setResonanceSerial(resonance.serial)
      resonances.append(resonance)
      initResonance(resonance)
    
    elif makeUnassignedResonances:
      resonance = nmrProject.newResonance(isotopeCode=fixedResonance.isotopeCode)
      fixedResonance.setResonanceSerial(resonance.serial)
      resonance.setName(fixedResonance.name)
      resonances.append(resonance)  

  return resonances

def makeFixedResonance(nmrConstraintStore,resonance):
  """
  Make a new fixed resonance for an NMR constraint top object based on a normal
  resonance. The fixed resonance will preserve assignment information at the
  time of constraint generation.
  NB should work equally well for fixedResonance input

  .. describe:: Input
  
  Nmr.NmrConstraintStore, Nmr.Resonance
  
  .. describe:: Output
  
  NmrConsraints.FixedResonance
  """

  fixedResonance = nmrConstraintStore.newFixedResonance(resonanceSerial=resonance.serial,
                                                        isotopeCode=resonance.isotopeCode,
                                                        name=resonance.name)
  transferCovalentlyBound(resonance, fixedResonance)
  return fixedResonance

def isValidFixedResonance(fixedResonance):

  resonance = fixedResonance.resonance

  if not resonance:
    return False

  resonanceSet = resonance.resonanceSet
  if not resonanceSet:
    return False

  fixedResonanceSet = fixedResonance.resonanceSet
  if not fixedResonanceSet:
    return False

  atoms = set()
  for atomSet in resonanceSet.atomSets:
    atoms.update(atomSet.atoms)

  for fixedAtomSet in fixedResonanceSet.atomSets:
    fixedAtoms = set(fixedAtomSet.atoms)
    if not fixedAtoms.issubset(atoms):
      return False

  return True

def getFixedResonance(nmrConstraintStore,resonance):
  """
  Find or create a fixed resonance for an NMR constraint top  object equivalent
  to the input normal resonance. The fixed resonance will preserve assignment
  information at the time of constraint generation.
  NB should work equally well for fixedResonance input
  
  .. describe:: Input
  
  Nmr.NmrConstraintStore, Nmr.Resonance
  
  .. describe:: Output
  
  NmrConsraints.FixedResonance
  """
  
  if not resonance:
    return None
  
  if not hasattr(nmrConstraintStore, 'quickResonances'):
    quickDict = {}
    nmrConstraintStore.quickResonances = quickDict
  else:
    quickDict = nmrConstraintStore.quickResonances
    
  serial = resonance.serial
  fixedResonance = quickDict.get(serial)

  if fixedResonance and fixedResonance.isDeleted:
    del quickDict[serial]
    fixedResonance = None
  
  if not fixedResonance:
    fixedResonance = nmrConstraintStore.findFirstFixedResonance(resonanceSerial=serial)

  if fixedResonance and not isValidFixedResonance(fixedResonance):
    # we should not delete but it is no longer validly pointing to resonance
    fixedResonance.resonanceSerial = None
    fixedResonance = None
  
  if not fixedResonance:
    fixedResonance = makeFixedResonance(nmrConstraintStore,resonance)
    fixedResonances = [fixedResonance,]
  
    if resonance.resonanceSet:
      fixedAtomSets = []
      for atomSet in resonance.resonanceSet.atomSets:
        fixedAtomSets.append( getFixedAtomSet(nmrConstraintStore, atomSet.atoms) )

      for resonance2 in resonance.resonanceSet.resonances:
        if resonance2 is not resonance:
          fixedResonance2 = makeFixedResonance(nmrConstraintStore,resonance2)
          quickDict[resonance2.serial] = fixedResonance2
          fixedResonances.append(fixedResonance2)

      nmrConstraintStore.newFixedResonanceSet(atomSets=fixedAtomSets,resonances=fixedResonances)
  
  quickDict[resonance.serial] = fixedResonance
  
  return fixedResonance

def makeNmrConstraintStore(nmrProject):
  """
  Make a new NMR constraint top object for a project which will contain
  constraints and violations.
  
  .. describe:: Input
  
  Nmr.NmrProject
  
  .. describe:: Output
  
  Nmr.NmrConstraintStore
  """

  project = nmrProject.root

  nmrConstraintStore = project.newNmrConstraintStore(nmrProject=nmrProject)
  nmrConstraintStore.quickResonances = {}
  nmrConstraintStore.quickAtomSets   = {}
  
  return nmrConstraintStore

def getIntensityDistanceTable(spectrum):
  """
  For a given spectrum ind or get a default table of NOE intensity to constraint
  distance relations.
  
  .. describe:: Input
  
  Nmr.DataSource
  
  .. describe:: Output
  
  List of Tuples of Floats (relative NOE intensity, target dist, min dist, max dist)
  """
  
  distanceClasses = getSpectrumNoeDistanceClasses(spectrum)
  
  if not distanceClasses:
    distanceClasses =  [(3.5,2.5,0.0,2.5),
                        (1.3,2.8,0.0,2.8),
                        (0.3,4.0,0.0,4.0),
                        (0.1,5.0,0.0,5.0),
                        (0.0,6.0,0.0,6.0)]

  return distanceClasses

def getDistancesFromIntensity(noeClasses, value):
  """
  Get constraining distances appropriate to a given NOE intensity value using a
  distance relation table.
  
  .. describe:: Input
  
  List of Tuples of Floats (relative NOE intensity, target dist, min dist, max dist), Float (NOE intensity/average for peak list)
  
  .. describe:: Output
  
  Float, Float, Float (targetValue, upperLimit, lowerLimit)
  """

  (targetValue, upperLimit, lowerLimit) = (6.0,6.0,0.0)
  for noeClass in noeClasses:
    if value >= noeClass[0]:
      (targetValue, upperLimit, lowerLimit) = noeClass[1:]
      break

  return (targetValue, upperLimit, lowerLimit)

def isResidueInRange(residue, residueRanges, dataDim):
  """
  Determine if a residue is in a residue range table for a given data dim. The
  range table lists residue bounds (first and last residue objects) for a chain
  appropriate to a list of data dims (often bonded).
  
  .. describe:: Input
  
  MolSystem.Residue, List of (List of Nmr.DataDims, MolSystem.Chain,
  MolSystem.Residue, MolSystem.Residue), Nmr.DataDim
  
  .. describe:: Output
  
  Boolean
  """
  
  for (dataDims, chain, startResidue, endResidue) in residueRanges:
    if dataDim in dataDims:
      if residue.chain is chain:
        if residue.seqCode >= startResidue.seqCode:
          if residue.seqCode <= endResidue.seqCode:
            return True
  
  return False

def getMeanPeakIntensity(peaks, intensityType='volume'):
  """
  Calculate the mean of the unsigned intensities of input peaks.
  
  .. describe:: Input
  
  List of Nmr.Peaks, String (Nmr.PeakIntensity.intensityType)
  
  .. describe:: Output
  
  Float
  """
     
  sumV = 0
  n    = 0.0
  for peak in peaks:
    peakDims  = peak.sortedPeakDims()
    intensity = peak.findFirstPeakIntensity(intensityType=intensityType)
    
    if intensity:
      value = abs(intensity.value)
      sumV += value
      n    += 1.0
            
  return sumV/max(n,1.0)


def getNoeDistance(intensity, params):
  """
  Get target, minimum and maximum distance for an NOE intensity Given ISPA
  parameters, erros and limits 
  
  .. describe:: Input
  
  Float, 6 List of Floats 
  
  .. describe:: Output
  
  Float, Float, Float
  """

  refDist, negError, posError, absMin, absMax, power = params
  dist    = refDist / (abs(intensity) ** (1/power))
  dist    = max(min(dist, absMax), absMin)
  minDist = max(absMin,dist-(negError*dist))
  maxDist = min(absMax,dist+(posError*dist))

  return dist, minDist, maxDist

def getDistMinMax(intensityValue, intensityScale, resonances0, resonances1, distanceFunction, normalise=True, labelling=None):

  if normalise:
    isoCorr0 = 0.0
    isoCorr1 = 0.0

    weight0 = 1.0/len(resonances0) 
    weight1 = 1.0/len(resonances1) 

    for resonance, indirect in resonances0:
      resonanceSet = resonance.resonanceSet

      propAtoms = 1.0 

      if indirect:
        propAtoms = 0.0
        for resonanceB in indirect:
          if labelling:
            frac = max(0.1, getResonancePairLabellingFraction(resonance, resonanceB, labelling) or 1.0)
          else:
            frac = 1.0  

          propAtoms += frac

      elif labelling:
        propAtoms *= max(0.1, getResonanceLabellingFraction(resonance, labelling) or 1.0)

      fac = weight0/propAtoms
      isoCorr0 +=  fac

    for resonance, indirect in resonances1:
      resonanceSet = resonance.resonanceSet

      propAtoms = 1.0 

      if indirect:
        propAtoms = 0.0
        for resonanceB in indirect:
          if labelling:
            frac = max(0.1, getResonancePairLabellingFraction(resonance, resonanceB, labelling) or 1.0)
          else:
            frac = 1.0  

          propAtoms += frac

      elif labelling:
        propAtoms *= max(0.1, getResonanceLabellingFraction(resonance, labelling) or 1.0)

      fac = weight1/propAtoms
      isoCorr1 += fac

  else:
    isoCorr0 = 1
    isoCorr1 = 1

  (dist,minDist,maxDist) = distanceFunction(intensityValue*isoCorr0*isoCorr1/intensityScale)

  return (dist,minDist,maxDist)

def findPeakConstraintPairs(peak, residueRanges=None, distDataDims=None, distIndices=None, indirectDims=None):

  spectrum = peak.peakList.dataSource
  experiment = spectrum.experiment
  nmrProject = experiment.nmrProject

  if distDataDims is None:
    distDataDims = getThroughSpaceDataDims(spectrum)

  if distIndices is None:
    distIndices  = [dd.dim-1 for dd in distDataDims]

  if indirectDims is None:
    indirectDims = {}
    for dataDims in getIndirectDataDims(spectrum):
      if  set(dataDims) == set(distDataDims):
        isotopesDict = getIndirectThroughSpaceIsotopes(experiment)

        for dataDim in dataDims:
          expDimRef = dataDim.expDim.sortedExpDimRefs()[0]
          indirectDims[dataDim.dim] = isotopesDict[expDimRef]

  resonances0 = []
  resonances1 = []

  distDim0, distDim1 = distDataDims
  peakDims = peak.sortedPeakDims()

  if peakDims[distIndices[0]].peakDimContribs and peakDims[distIndices[1]].peakDimContribs:

    peakDim0 = peakDims[distIndices[0]]
    for contrib in peakDim0.peakDimContribs:
      resonance = contrib.resonance
      if resonance.resonanceSet:
        if residueRanges:
          residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
          if isResidueInRange(residue, residueRanges, distDim0):
            resonances0.append( (contrib.resonance, []) )
        else:
          resonances0.append( (contrib.resonance, []) )

    if peakDim0.dim in indirectDims and indirectDims[peakDim0.dim][1]:
      isotopeA, isotopeB = indirectDims[peakDim0.dim]
      chemElement = isotopeB.chemElement

      for resonance, indirect in resonances0:

        isotopeCode = '%d%s' % (isotopeB.massNumber, chemElement.symbol)

        # Use getBoundResonances, to get from Cga to Hga* (and not also Hgb*)
        resonancesA = set(x for x in getBoundResonances(resonance, recalculate=True)
                          if x.isotopeCode == isotopeCode
                          and x.resonanceSet)

        # get covalently bound atomSets
        atoms = set()
        for atomSet in resonance.resonanceSet.atomSets:
          atoms.update(getBoundAtoms(atomSet.findFirstAtom()))

        atomSets = set(a.atomSet for a in atoms if a.atomSet and \
                       a.chemAtom.chemElement is chemElement)

        if resonancesA:
          # remove covalently impossible resonances
          resonanceSets = set(y for x in atomSets for y in x.resonanceSets)
          resonancesA = set(x for x in resonancesA 
                            if x.resonanceSet in resonanceSets)

        if not resonancesA:
          # make new resonances to fit covalent atoms.
          for atomSet in atomSets:
            resonanceB = nmrProject.newResonance(isotopeCode=isotopeCode)
            assignAtomsToRes([atomSet,], resonanceB)
            resonancesA.add(resonanceB)

        indirect.extend(resonancesA)

    peakDim1 = peakDims[distIndices[1]]
    for contrib in peakDim1.peakDimContribs:
     resonance = contrib.resonance
     if resonance.resonanceSet and (resonance not in resonances0):
       if residueRanges:
         residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
         if isResidueInRange(residue, residueRanges, distDim1):
           resonances1.append( (contrib.resonance, []) )
       else:
         resonances1.append( (contrib.resonance, []) )

    if peakDim1.dim in indirectDims and indirectDims[peakDim1.dim][1]:
      isotopeA, isotopeB = indirectDims[peakDim1.dim]
      chemElement = isotopeB.chemElement

      for resonance, indirect in resonances1:

        isotopeCode = '%d%s' % (isotopeB.massNumber, chemElement.symbol)

        resonancesA = set(x for x in getBoundResonances(resonance, recalculate=True)
                          if x.isotopeCode == isotopeCode
                          and x.resonanceSet)

        atoms = set()
        for atomSet in resonance.resonanceSet.atomSets:
          atoms.update(getBoundAtoms(atomSet.findFirstAtom()))

        atomSets = set(a.atomSet for a in atoms if a.atomSet and \
                       a.chemAtom.chemElement is chemElement)

        if resonancesA:
          # remove covalently impossible resonances
          resonanceSets = set(y for x in atomSets for y in x.resonanceSets)
          resonancesA = set(x for x in resonancesA 
                            if x.resonanceSet in resonanceSets)

        if not resonancesA:  
          for atomSet in atomSets:
            resonanceB = nmrProject.newResonance(isotopeCode=isotopeCode)
            assignAtomsToRes([atomSet,], resonanceB)
            resonancesA.add(resonanceB)

        indirect.extend(resonancesA)       

  return resonances0, resonances1


def makeDistConstraints(peakList, constraintSet=None, intensityType='volume',
                        distanceFunction=None, normalise=True, labelling=None,
                        params=None, residueRanges=None, minMerit=0.0, scale=None,
                        scaleDict=None):
  """
  Makes a constraint list with constraints based upon the assigned peaks within
  a NOESY peak list. Constraints will be put in a new NMR constraint head object
  if none is specified. An NOE intensity-to-distance function, minimum peak
  merit and allowed residue ranges may be input for the calculation. A
  normalising peak list can be specified to help calibrate the peak intensities
  according to the sensitivity of a corresponding root resonances. The scale
  option is the value by which peak intensities are scaled for NOE
  table/function lookup. Params relate to the generic NOE distance function if
  neither these nor a distance function is specified a lookuptable is used.
  
  .. describe:: Input
  
  Nmr.PeakList (NOESY), Nmr.NmrConstraintStore,
  String (Nmr.PeakIntensity.intensityType)
  Function (to get constraining distances), Boolean
  ChemComLabel.LabelingScheme or True (automatic from experiment MolLabel),
  List of Floats (distance function parameters),
  List of (List of Nmr.DataDims, Nmr.Chain, Integer, Integer),
  Float (Nmr.Peak.figOfMerit), Float (baseline intensity) 
  Dict of Nmr.Peak --> Float (peak-specific scale)
  
  .. describe:: Output
  
  NmrConstraint.DistanceConstraintList
  """
  
  if not scaleDict:
    scaleDict = {}

  spectrum = peakList.dataSource
  #peaks  = peakList.peaks
  experiment = spectrum.experiment
  nmrProject = experiment.nmrProject

  distDataDims = getThroughSpaceDataDims(spectrum)
  distIndices  = [dd.dim-1 for dd in distDataDims]

  if len(distDataDims) != 2:
    msg = 'Experiment appears to not have two through-space linked dimensions. '
    msg += 'Check experiment type and dim-dim transfer settings.' 
    showWarning('Failure', msg)
    return

  distDim0, distDim1 = distDataDims

  if labelling is True:
    labelling = experiment
    
  if not residueRanges:
    residueRanges = None

  #workingPeaks = []
  #for peak in peaks:
  #  if peak.figOfMerit < minMerit:
  #    continue
  #  workingPeaks.append(peak)
  workingPeaks = [x for x in peakList.sortedPeaks() if x.figOfMerit >= minMerit]
  
  mean = getMeanPeakIntensity(workingPeaks, intensityType=intensityType)
  if scale: # neither Zero nor None
    mean = scale
  
  if not mean: 
    msg  = 'Cannot make restraints: peak %s is zero on average.' % intensityType
    msg += ' Maybe intensities are missing or the peak list is empty' 
    showWarning('Failure', msg)
    return 
   
  if not constraintSet:
    constraintSet = makeNmrConstraintStore(peakList.topObject)
  distConstraintList = constraintSet.newDistanceConstraintList()
  distConstraintList.addExperimentSerial(experiment.serial)
  newConstraint = distConstraintList.newDistanceConstraint
  
  if not distanceFunction:
    if params:
      distanceFunction = lambda val:getNoeDistance(val, params)
    else:
      noeDistClasses = getIntensityDistanceTable(spectrum)
      distanceFunction = lambda val:getDistancesFromIntensity(noeDistClasses,val)
  
  # Check for indirect transfers
  indirectDims = {}
  for dataDims in getIndirectDataDims(spectrum):
    if  set(dataDims) == set(distDataDims):
      isotopesDict = getIndirectThroughSpaceIsotopes(experiment)
      
      for dataDim in dataDims:
        expDimRef = dataDim.expDim.sortedExpDimRefs()[0]
        indirectDims[dataDim.dim] = isotopesDict[expDimRef]
  
  
  for peak in workingPeaks:

    if peak.figOfMerit < minMerit:
      continue

    intensity = peak.findFirstPeakIntensity(intensityType=intensityType)
    if not intensity:
      continue
    intensityValue = abs(intensity.value)

    resonances0, resonances1 = findPeakConstraintPairs(peak, residueRanges, distDataDims, distIndices, indirectDims)
    if resonances0 and resonances1:

      peakDims = peak.sortedPeakDims()
      peakDim0 = peakDims[distIndices[0]]
      peakDim1 = peakDims[distIndices[1]]

      # Filter by correlated contributions
      contribFilter = set()
      for peakContrib in peak.peakContribs:
        contribs = peakContrib.peakDimContribs
        for contrib0 in contribs:
          if contrib0.peakDim is peakDim0:
            for contrib1 in contribs:
              if contrib1.peakDim is peakDim1:
                resonancesF = (contrib0.resonance, contrib1.resonance)
                contribFilter.add(resonancesF)

      peakMean = scaleDict.get(peak)
      if peakMean is None:
        if peak in scaleDict: # in this case the peak is in the dict but with a value of None
          continue # and that means that we should skip this peak
        peakMean = mean
      (dist,minDist,maxDist) = getDistMinMax(intensityValue, peakMean, resonances0, resonances1, distanceFunction, normalise=normalise, labelling=labelling)
      error = abs(maxDist - minDist)

      fResonancePairs = set()                                     
      for resonance0, indirect0 in resonances0:
        resonances0 = indirect0 or [resonance0,]
        fixedResonances0 = [getFixedResonance(constraintSet,r) for r in resonances0]

        for resonance1, indirect1 in resonances1:
          if contribFilter:
            if (resonance0, resonance1) not in contribFilter:
              continue

          resonances1A = indirect1 or [resonance1,]
          fixedResonances1 = [getFixedResonance(constraintSet,r) for r in resonances1A]

          for fixedRes0 in fixedResonances0:
            for fixedRes1 in fixedResonances1:
              if fixedRes0 is not fixedRes1:
                fResonancePairs.add(frozenset([fixedRes0, fixedRes1]))

      if fResonancePairs:
        # Otherwise you generate restraints between nothing and nothing

        constraint  = newConstraint(weight=1.0, origData=intensityValue, targetValue=dist,
                                    upperLimit=maxDist, lowerLimit=minDist, error=error)

        constraint.newConstraintPeakContrib(experimentSerial=experiment.serial,
                                            dataSourceSerial=spectrum.serial,
                                            peakListSerial=peakList.serial,
                                            peakSerial=peak.serial)

        for fixedRes0, fixedRes1 in fResonancePairs:
          constraint.newDistanceConstraintItem(resonances=[fixedRes0,fixedRes1])

  return distConstraintList

def getPeakDimTolerance(peakDim, minTol, maxTol, multiplier):
  """
  Finds the shift matching tolerance for a peak dim based on its line width.
  Defaults to the minimum tolerance if no line width is present.
  
  .. describe:: Input
  
  Nmr.PeakDim, Float, Float, Float (Nmr,PeakDim.lineWidth to tolerance factor)
  
  .. describe:: Output
  
  Float (PPM tolerance)
  """

  tolerance = None
  if peakDim.dataDimRef:
    # works in ppm
    tolerance = minTol
    if peakDim.lineWidth:
      zeroVal = pnt2ppm(0,peakDim.dataDimRef)
      width   = pnt2ppm(multiplier*peakDim.lineWidth,peakDim.dataDimRef)
      width  -= zeroVal
      tolerance = min(maxTol, max(minTol, width))
      
  return tolerance

def makeAmbigDistConstraints(peakList, tolerances, chemShiftRanges, constraintSet=None,
                             testOnly=False, labelling=None, minLabelFraction=0.1,
                             distanceFunction=None, residueRanges=None, minMerit=0.0, progressBar=None,
                             intensityType='volume', ignoreDiagonals=True, doAliasing=True,
                             structure=None, maxDist=None,
                             scale=None, params=None, peakCategories=None):
  """
  Makes a constraint list with constraints by matching known shifts
  in given range and within specified tolerances to a NOESY peak
  list. Optional labelling scheme/mixture and labelling threshold to filter
  according to a given set of residue isotopomers.
  Constraints will be put in a new NMR constraint store object
  if none is specified. Peaks are catergorised into various lists. 
  A structure and max distance can be used to filter contributions.
  The scale option is the value by which peak intensities are scaled
  for NOE table/function lookup. Params relate to the generic NOE
  distance function if neither these nor a distance function is
  specified a lookuptable is used
  
  .. describe:: Input
  
  .PeakList (NOESY), List of (Nmr.DataDim, Float, Float, Float) (shift tolerances),
  List of (Nmr.DataDim, String (isotope code), Float, Float) (chem shift ranges)
  Nmr.NmrConstraintStore, Boolean (test only),
  ChemCompLabel.LabelingScheme or True (automatic from experiment MolLabel), Float,
  Function (to get distances from NOEs),  Nmr.PeakIntensity.intensityType,
  List of (List of Nmr.DataDims, Nmr.Chain, Integer, Integer) (residue ranges),
  Float (min Peak.figOfMerit), ProgressBar (Analysis popup),
  Boolean, Boolean, MolStructure.StructureEnsemble, Float
  Float, List of Floats (distance function parameters)
  Dict (for Category Name:Lists of Nmr.Peaks)
  
  .. describe:: Output
  
  NmrConstraint.DistanceConstraintList
  """
  
  from ccpnmr.analysis.core.StructureBasic import getAtomSetsDistance
  
  if peakCategories is None:
    peakCategories = {}
    
  assignedPeaks = peakCategories['Assigned'] = []
  diagonalPeaks = peakCategories['Diagonal'] = []
  unmatchedPeaks = peakCategories['Unmatchable'] = []
  poorMeritPeaks = peakCategories['Poor Merit'] = []
  outOfRangePeaks = peakCategories['Out of range'] = []
  distalPeaks = peakCategories['Too Distal'] = []
      
  #peaks = peakList.peaks
  spectrum = peakList.dataSource
  experiment = spectrum.experiment
  nmrProject = experiment.nmrProject 
  distDataDims = getThroughSpaceDataDims(spectrum)
  distIndices  = [dd.dim-1 for dd in distDataDims]

  if len(distDataDims) != 2:
    return
    
  distDim1, distDim2 = distDataDims
  
  if labelling is True:
    labelling = experiment
    
  if not residueRanges:
    residueRanges = None
  
  bondedDims = {}  
  for dataDim1, dataDim2 in getOnebondDataDims(spectrum):
    bondedDims[dataDim1] = dataDim2
    bondedDims[dataDim2] = dataDim1
    
  if testOnly:
    distConstraintList = None
  else:
    if not constraintSet:
      constraintSet = makeNmrConstraintStore(experiment.topObject)
    distConstraintList = constraintSet.newDistanceConstraintList()
    distConstraintList.addExperimentSerial(experiment.serial)
    newConstraint = distConstraintList.newDistanceConstraint

  tolDict = {}
  for (dataDim,minT,maxT,multi) in tolerances:
    tolDict[dataDim] = (minT,maxT,multi)
  
  chemShiftRangesDict = {} 
  for (dataDim, iso, minShift, maxShift) in chemShiftRanges:
    if chemShiftRangesDict.get(dataDim) is None:
      chemShiftRangesDict[dataDim] = []
      
    chemShiftRangesDict[dataDim].append([minShift, maxShift])
    if tolDict.get(dataDim) is None:
      msg = 'No tolerance set for dataDim %s of dataSource %s' % (dataDim,spectrum)
      raise Exception(msg)
  
  # go through peaks
  # if not assigned in all the Hydrogen dims or H + bonded dim
  
  # Check for indirect transfers
  indirectDims = {}
  for dataDims in getIndirectDataDims(spectrum):
    if set(dataDims) == set(distDataDims):
      isotopesDict = getIndirectThroughSpaceIsotopes(experiment)
      
      for dataDim in dataDims:
        expDimRef = dataDim.expDim.sortedExpDimRefs()[0]
        indirectDims[dataDim.dim] = isotopesDict[expDimRef]

  workingPeaks = []
  for peak in peakList.sortedPeaks():
    # filter out diagonals
    if ignoreDiagonals:
      peakDims = peak.sortedPeakDims()
      peakDim1 = peakDims[distIndices[0]]
      peakDim2 = peakDims[distIndices[1]]
      ppm1 = peakDim1.realValue
      ppm2 = peakDim2.realValue
      
      delta = abs(ppm1-ppm2)
      
      if (delta <= tolDict[distDim1][0] ) or (delta <= tolDict[distDim2][0]):
        dataDimA = bondedDims.get(distDim1)
        dataDimB = bondedDims.get(distDim2)
        
        if dataDimA and dataDimB :
          peakDimA = peak.findFirstPeakDim(dataDim=dataDimA)
          peakDimB = peak.findFirstPeakDim(dataDim=dataDimB)
          ppmA = pnt2ppm(peakDimA.position,peakDimA.dataDimRef)
          ppmB = pnt2ppm(peakDimB.position,peakDimB.dataDimRef)
          
          delta2 = abs(ppmA-ppmB)
          if (delta2 <= tolDict[dataDimA][0] ) or (delta2 <= tolDict[dataDimB][0]):
            diagonalPeaks.append(peak)
            continue
        
        else:
          diagonalPeaks.append(peak)
          continue

    if peak.figOfMerit < minMerit:
      poorMeritPeaks.append(peak)
      continue

    workingPeaks.append(peak)
    
  mean = getMeanPeakIntensity(workingPeaks, intensityType=intensityType)
  if scale: # neither Zero nor None
    mean = scale
  
  if not mean: 
    msg  = 'Cannot make restraints: peak %s is zero on average.' % intensityType
    msg += ' Maybe intensities are missing or the peak list is empty' 
    showWarning('Failure', msg)
    return 
 
  
  if not distanceFunction:
    if params:
      distanceFunction = lambda val:getNoeDistance(val, params)
    else:
      noeDistClasses = getIntensityDistanceTable(spectrum)
      distanceFunction = lambda val:getDistancesFromIntensity(noeDistClasses,val)
  
  for peak in workingPeaks:
    
    if progressBar:
      progressBar.increment()

    intensity = peak.findFirstPeakIntensity(intensityType=intensityType)
    if not intensity:
      continue
    intensityValue = abs(intensity.value)

    outOfRange = 0
       
    unassignedPeakDims = []
    peakDims = peak.sortedPeakDims()
    
    outOfShiftRange = False
    for peakDim in peakDims:
      inRange = isShiftInRange(peakDim.realValue,chemShiftRangesDict[peakDim.dataDim])
      if chemShiftRanges and not inRange:
        outOfShiftRange = True
        break
      

    if outOfShiftRange:
      unmatchedPeaks.append(peak)
      continue
           
    #n = 0
    for i in distIndices:
      peakDim = peakDims[i]
      if not peakDim.peakDimContribs:
        unassignedPeakDims.append( peakDim )
    
     
    # filter out assigned peaks
    if not unassignedPeakDims:
      assignedPeaks.append(peak)
      continue
    
    peakResonances = []
    for i in distIndices:
      resonances = [] 
      peakDim = peakDims[i]
      dataDim = peakDim.dataDim
      
      if peakDim in unassignedPeakDims:
        #isotope    = dataDim.expDim.findFirstExpDimRef().isotopeCodes[0]
        #peakDimPos = peakDim.position + (peakDim.numAliasing*dataDim.numPointsOrig)
        (minT,maxT,multi) = tolDict[dataDim]
        tolerance  = getPeakDimTolerance(peakDim,minT,maxT,multi)
        
        bondedDim = bondedDims.get(peakDim.dataDimRef.dataDim)

        if bondedDim:
          # check that both bonded dim possibilities are within tolerances
          
          shifts = findMatchingPeakDimShifts(peakDim,
                                             chemShiftRangesDict[dataDim],
                                             tolerance=tolerance,
                                             aliasing=doAliasing,
                                             findAssigned=True)
            
          if shifts:
            for peakDim2 in peakDims:
              if peakDim2.dataDimRef.dataDim is bondedDim:
                
                shifts2 = []
                if peakDim2.peakDimContribs:
                  for contrib in peakDim2.peakDimContribs:
                    shift = contrib.resonance.findFirstShift(parentList=experiment.shiftList)
                    if shift:
                      shifts2.append(shift)
                    
                else:
                  dataDim2    = peakDim2.dataDim
                  (minT,maxT,multi) = tolDict[dataDim2]
                  tolerance2  = getPeakDimTolerance(peakDim2,minT,maxT,multi)

                  shifts2 = findMatchingPeakDimShifts(peakDim2,
                                                      chemShiftRangesDict[dataDim2],
                                                      tolerance=tolerance2,
                                                      aliasing=doAliasing,
                                                      findAssigned=True)

                for shift in shifts:
                  resonance = shift.resonance
                
                  for shift2 in shifts2:
                    resonance2 = shift2.resonance
                  
                    if areResonancesBound(resonance, resonance2):
                      if residueRanges:
                        residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
                        if isResidueInRange(residue, residueRanges, dataDim):
                          resonances.append((resonance, resonance2, []))
                        else:
                          outOfRange += 1
                      else:
                        resonances.append((resonance, resonance2, []))
                      break

        else:

          shifts = findMatchingPeakDimShifts(peakDim,
                                             chemShiftRangesDict[dataDim],
                                             tolerance=tolerance,
                                             aliasing=doAliasing,
                                             findAssigned=True)
                      
          for shift in shifts:
            resonance = shift.resonance
          
            if residueRanges:
              residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
              if isResidueInRange(residue, residueRanges, dataDim):
                resonances.append((resonance, None, []))
              else:
                outOfRange += 1
            else:
              resonances.append((resonance, None, []))
        
      else:
        # this dim is assigned
        for contrib in peakDim.peakDimContribs:
          resonance = contrib.resonance
          resonanceSet = resonance.resonanceSet
        
          if resonanceSet:
            if residueRanges:
              residue = resonanceSet.findFirstAtomSet().findFirstAtom().residue
              if isResidueInRange(residue, residueRanges, dataDim):
                resonances.append((resonance, None, []))
              else:
                outOfRange += 1
            else:
              resonances.append((resonance, None, []))
      
      # Deal with indirect transfers
      
      if peakDim.dim in indirectDims and indirectDims[peakDim.dim][1]:
        isotopeA, isotopeB = indirectDims[peakDim.dim]
        chemElement = isotopeB.chemElement
        
        for resonance, bound, indirect in resonances:
          
          isotopeCode = '%d%s' % (isotopeB.massNumber, chemElement.symbol)
          
          # Use getBoundResonance to get from e.g. Cga to Hga* and not Hgb*
          resonancesA = set(x for x in getBoundResonances(resonance, recalculate=True)
                            if x.isotopeCode == isotopeCode
                            and x.resonanceSet)
          
          # get covalently bound atomSts
          atoms = set()
          for atomSet in resonance.resonanceSet.atomSets:
            atoms.update(getBoundAtoms(atomSet.findFirstAtom()))
          
          atomSets = set(a.atomSet for a in atoms if a.atomSet and \
                         a.chemAtom.chemElement is chemElement)
                         
 
          if resonancesA:
            # remove covalently impossible resonances
            resonanceSets = set(y for x in atomSets for y in x.resonanceSets)
            resonancesA = set(x for x in resonancesA 
                              if x.resonanceSet in resonanceSets)
          
          if not resonancesA:
            # make new resonances for covanlently bound atoms
            for atomSet in atomSets:
              resonanceB = nmrProject.newResonance(isotopeCode=isotopeCode)
              assignAtomsToRes([atomSet,], resonanceB)
              resonancesA.add(resonanceB)
          
          indirect.extend(resonancesA)
       
      # Store resonances for this dim
      
      peakResonances.append( resonances )
      
    if peakResonances[0] and peakResonances[1]:
      distal = False
      
      resonancePairs = set()
      for resonance0, bound0, indirect0 in peakResonances[0]:
        resonanceSet0 = resonance0.resonanceSet
        
        for resonance1, bound1, indirect1 in peakResonances[1]:
          if resonance1 is resonance0:
            continue
        
          if labelling:
            if bound0:
              fraction0 = getResonancePairLabellingFraction(resonance0,
                                                            bound0,
                                                            labelling)
            else:
              fraction0 = getResonanceLabellingFraction(resonance0,
                                                        labelling)
            if fraction0 < minLabelFraction:
              continue
                                                           
            if bound1:
              fraction1 = getResonancePairLabellingFraction(resonance1,
                                                            bound1,
                                                            labelling)
            else:
              fraction1 = getResonanceLabellingFraction(resonance1,
                                                        labelling)
                                                        
            if fraction1 < minLabelFraction:
              continue
        
          if structure and resonanceSet0 and (maxDist is not None):
            resonanceSet1 = resonance1.resonanceSet
            
            if resonanceSet1:
              atomSets0 = list(resonanceSet0.atomSets)
              atomSets1 = list(resonanceSet1.atomSets)
              dist = getAtomSetsDistance(atomSets0, atomSets1, structure, method='noe')
              
              if dist > maxDist:
                distal = True
                continue
          
          resonances0 = indirect0 or [resonance0,]
          resonances1 = indirect1 or [resonance1,]
          
          for resonanceA in resonances0:
            for resonanceB in resonances1:
              if resonanceA is not resonanceB:
                resonancePairs.add(frozenset([resonanceA, resonanceB]))
      
      if not resonancePairs:
        unmatchedPeaks.append(peak)
        
        if distal:
          distalPeaks.append(peak)
    
      elif not testOnly:
        resonances0 = [(resonance,indirect) for resonance, bound, indirect in peakResonances[0]]
        resonances1 = [(resonance,indirect) for resonance, bound, indirect in peakResonances[1]]
        (dist,minDistL,maxDistL) = getDistMinMax(intensityValue, mean, resonances0, resonances1, distanceFunction, labelling=labelling)
        error = abs(maxDistL - minDistL)
        constraint  = newConstraint(weight=1.0, origData=intensityValue,
                                    targetValue=dist, upperLimit=maxDistL,
                                    lowerLimit=minDistL, error=error)
        
        constraint.newConstraintPeakContrib(experimentSerial=experiment.serial,
                                            dataSourceSerial=spectrum.serial,
                                            peakListSerial=peakList.serial,
                                            peakSerial=peak.serial)

        for resonance0, resonance1 in resonancePairs:
          fixedResonance0 = getFixedResonance(constraintSet,resonance0)
          fixedResonance1 = getFixedResonance(constraintSet,resonance1)
          constraint.newDistanceConstraintItem(resonances=[fixedResonance0,fixedResonance1])
 
    elif outOfRange > 1:
      outOfRangePeaks.append(peak)
    else:
      unmatchedPeaks.append(peak)

  return distConstraintList

def transferCovalentlyBound(sourceResonance, targetResonance):
  """
  Set targetResonance.covalentlyBound to match sourceResonance.covalentlyBound
  Uses serial for Resonance, and resonanceSerial for FixedResonance
  
  .. describe:: Input
  
  Nmr.Resonance or NmrConstraint.FixedResonance (twice)
  
  .. describe:: Output
  
  None
  """
  
  if sourceResonance.className == 'FixedResonance':
    covalentSerials = set(x.resonanceSerial 
                       for x in sourceResonance.covalentlyBound)
    if None in covalentSerials:
      covalentSerials.remove(None)
  else:
    # normal resonance, make sure we update status
    covalentSerials = set(x.serial for x in getBoundResonances(sourceResonance,
                                                            recalculate=True))
  
  topObj = targetResonance.topObject
  if targetResonance.className == 'FixedResonance':
    boundRes = set(topObj.findFirstFixedResonance(resonanceSerial=x) 
                   for x in covalentSerials)
  else:
    boundRes = set(topObj.findFirstResonance(serial=x) 
                   for x in covalentSerials)
  if None in boundRes:
    boundRes.remove(None)
    
  targetResonance.covalentlyBound = boundRes

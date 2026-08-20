
LICENSE = """
======================COPYRIGHT/LICENSE START==========================

AssignmentBasic.py: Part of the CcpNmr Analysis program

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

from math import sqrt

import re, operator

from memops.api.Implementation import AppDataString
from memops.general import Implementation

try:
  from memops.gui.MessageReporter import showOkCancel, showWarning, showYesNo
except ImportError:
  from memops.universal.MessageReporter import showOkCancel, showWarning, showYesNo

from ccpnmr.analysis.core.MergeObjects import mergeObjects
from ccpnmr.analysis.core.Util import getAnalysisDataDim
from ccpnmr.analysis.core.UnitConverter import pnt2ppm, ppm2pnt, unit_converter

from ccp.util.LabeledMolecule import getIsotopomerSingleAtomFractions, getIsotopomerAtomPairFractions
from ccp.util.LabeledMolecule import singleAtomFractions, atomPairFractions

from ccp.api.nmr import Nmr

def getResonanceLabellingFraction(resonance, labelling):
  """
  Get the fraction of labelling for a given resonance's assignment
  or make a guess if it is atom typed and in a residue typed spin system.
  Can work with a reference isotopomer scheme or a labeled mixture.

  .. describe:: Input
  
  Nmr.Resonance, ChemComLabel.LabelingScheme or LabeledMolecule.LabeledMixture
  
  .. describe:: Output
  
  Float
  """

  fraction = 1.0 # In the absence of any assignment
  
  resonanceSet = resonance.resonanceSet
  
  if labelling.className == 'LabelingScheme':
    if resonanceSet:
      atomSets = list(resonanceSet.atomSets)
      residue = atomSets[0].findFirstAtom().residue
      ccpCode = residue.ccpCode
      molType = residue.molType
      chemCompLabel = labelling.findFirstChemCompLabel(ccpCode=ccpCode,
                                                             molType=molType)
 
      if not chemCompLabel:
        natAbun = resonance.root.findFirstLabelingScheme(name='NatAbun')
 
        if natAbun:
          chemCompLabel = natAbun.findFirstChemCompLabel(ccpCode=ccpCode,
                                                         molType=molType)
 
      if chemCompLabel:
        isotopomers = chemCompLabel.isotopomers
        isotope = resonance.isotopeCode
 
        fractions = []
        for atomSet in atomSets:
          atoms = atomSet.atoms
          atomFrac = 0.0
 
          for atom in atoms:
            subType = atom.chemAtom.subType
 
            fracDict = getIsotopomerSingleAtomFractions(isotopomers,
                                                        atom.name, subType)
            atomFrac += fracDict.get(isotope, 1.0)
 
          atomFrac /= float(len(atoms))
 
          fractions.append(atomFrac)
 
        fraction = max(fractions)
 
    elif resonance.assignNames:
      atomNames = resonance.assignNames
      spinSystem = resonance.resonanceGroup
 
      if spinSystem and spinSystem.ccpCode:
        ccpCode = spinSystem.ccpCode
        molType = spinSystem.molType or 'protein'
        chemCompLabel = labelling.findFirstChemCompLabel(ccpCode=ccpCode,
                                                         molType=molType)
 
        if not chemCompLabel:
          natAbun = resonance.root.findFirstLabelingScheme(name='NatAbun')
 
          if natAbun:
            chemCompLabel = natAbun.findFirstChemCompLabel(ccpCode=ccpCode,
                                                           molType=molType)
 
        if chemCompLabel:
          isotopomers = chemCompLabel.isotopomers
          isotope = resonance.isotopeCode
          fraction = 0.0
 
          for atomName in atomNames:
            fracDict = getIsotopomerSingleAtomFractions(isotopomers,
                                                        atomName, 1)
            fraction += fracDict.get(isotope, 1.0)
 
          fraction /= float(len(atomNames))
  
  else: # get from experiment abelled mixture
    
    if resonanceSet:
      atomSets = list(resonanceSet.atomSets)
      isotope = resonance.isotopeCode
      labelledMixtures = labelling.labeledMixtures
      molResidue = atomSets[0].findFirstAtom().residue.molResidue
      molecule = molResidue.molecule
      resId = molResidue.serial
      
      for mixture in labelledMixtures:
        if mixture.labeledMolecule.molecule is molecule:
          fractions = []
          for atomSet in atomSets:
            atoms = atomSet.atoms
            atomFrac = 0.0
 
            for atom in atoms:
              fracDict = singleAtomFractions(mixture, resId, atom.name)
              atomFrac += fracDict.get(isotope, 1.0)
 
            atomFrac /= float(len(atoms))
 
            fractions.append(atomFrac)
 
          fraction = max(fractions)
          break
          
  return fraction


def getResonancePairLabellingFraction(resonanceA, resonanceB, labelling):
  """
  Get the fraction of a pair of resonances both being labelled
  given a labelling scheme. Considers individual isotopomers if
  the resonances are bound within the same residue.
  Can work with a reference isotopomer scheme or a labeled mixture.

  .. describe:: Input
  
  Nmr.Resonance, Nmr.Resonance,
  ChemComLabel.LabelingScheme or LabeledMolecule.LabeledMixture
  
  .. describe:: Output
  
  Float
  """
  
  from ccpnmr.analysis.core.MoleculeBasic import areResonancesBound
  
  fraction = 1.0 # In the absence of any assignment
  
  resonanceSetA = resonanceA.resonanceSet
  resonanceSetB = resonanceB.resonanceSet
  
  if resonanceSetA and resonanceSetB:
    isotopes = (resonanceA.isotopeCode, resonanceB.isotopeCode)
    atomA = resonanceSetA.findFirstAtomSet().findFirstAtom()
    atomB = resonanceSetB.findFirstAtomSet().findFirstAtom()
    residueA = atomA.residue
    residueB = atomB.residue
    
    if labelling.className == 'LabelingScheme':
      findFirstChemCompLabel = labelling.findFirstChemCompLabel
  
      subTypeA = atomA.chemAtom.subType
      subTypeB = atomB.chemAtom.subType
 
      if residueA is residueB:
        chemCompLabel = findFirstChemCompLabel(ccpCode=residueA.ccpCode,
                                               molType=residueA.molType)
 
        if not chemCompLabel:
          natAbun = resonanceA.root.findFirstLabelingScheme(name='NatAbun')
 
          if natAbun:
            chemCompLabel = natAbun.findFirstChemCompLabel(ccpCode=residueA.ccpCode,
                                                           molType=residueA.molType)
        if not chemCompLabel:
          return 1.0 # Nothing can be done, no isotopomers
 
        isotopomers  = chemCompLabel.isotopomers
 
        fractions = []
        for atomSetA in resonanceSetA.atomSets:
          for atomSetB in resonanceSetB.atomSets:
 
            n = 0.0
            pairFrac = 0.0
            for atomA in atomSetA.atoms:
              nameA = atomA.name
              subTypeA = atomA.chemAtom.subType
 
              for atomB in atomSetB.atoms:
                atomNames = (nameA, atomB.name)
                subTypes  = (subTypeA, atomB.chemAtom.subType)
                pairDict  = getIsotopomerAtomPairFractions(isotopomers, atomNames, subTypes)
                pairFrac += pairDict.get(isotopes, 1.0)
                n += 1.0
 
            pairFrac /= n
            fractions.append(pairFrac)
 
        fraction = max(fractions)

      else: # Assumes filly mixed
        fractionA = getResonanceLabellingFraction(resonanceA, labelling)
        fractionB = getResonanceLabellingFraction(resonanceB, labelling)
        fraction = fractionA * fractionB
     
    else: # Get Labelling mixture from experiment
      molResidueA = residueA.molResidue
      molResidueB = residueB.molResidue
      resIds = (molResidueA.serial, molResidueB.serial)
      labelledMixtures = labelling.labeledMixtures
      
      moleculeA = molResidueA.molecule
      moleculeB = molResidueA.molecule
      
      if moleculeA is moleculeB:
        for mixture in labelledMixtures:
          if mixture.labeledMolecule.molecule is moleculeA:
            fractions = []
            for atomSetA in resonanceSetA.atomSets:
              for atomSetB in resonanceSetB.atomSets:

                n = 0.0
                pairFrac = 0.0
                for atomA in atomSetA.atoms:
                  nameA = atomA.name

                  for atomB in atomSetB.atoms:
                    atomNames = (nameA, atomB.name)
 
                    pairDict  = atomPairFractions(mixture, resIds, atomNames)
                    pairFrac += pairDict.get(isotopes, 1.0)
                    n += 1.0

                pairFrac /= n
                fractions.append(pairFrac)

            fraction = max(fractions)
            break
      else:
        fractionA = getResonanceLabellingFraction(resonanceA, labelling)
        fractionB = getResonanceLabellingFraction(resonanceB, labelling)
        fraction = fractionA * fractionB
     
  else:
    fractionA = getResonanceLabellingFraction(resonanceA, labelling)
    fractionB = getResonanceLabellingFraction(resonanceB, labelling)
    fraction = fractionA * fractionB
      
  return fraction


def getShiftLists(nmrProject):
  """
  Get all the shift lists associated with an NMR project.

  .. describe:: Input
  
  Nmr.NmrProject
  
  .. describe:: Output
  
  Nmr.ShiftList
  """

  if not nmrProject:
    return []

  sls = [(sl.serial, sl) for sl in nmrProject.findAllMeasurementLists(className='ShiftList')]
  sls.sort()
  
  return [x[1] for x in sls]


def getSyntheticShiftLists(nmrProject):
  """
  Get all the synthetic shift lists associated with an NMR project.

  .. describe:: Input
  
  Nmr.NmrProject
  
  .. describe:: Output
  
  Nmr.ShiftList
  """
     
   #NBNB Filtering on isSimulated might be enough, but is it reliable?

  if not nmrProject:
    return []
  
  result = [x for x in nmrProject.sortedMeasurementLists()
            if x.className=='ShiftList'
            and x.findFirstExperiment() is None
            #and x.isSimulated # disabled pending check on how to set it
           ]
  
  return result

def updateResonanceAnnotation(resonance):
  """
  Update the annotations for the peaks assigned to a given resonance

  .. describe:: Input
  
  Nmr.Resonance
  
  .. describe:: Output
  
  None 
  """
  
  if hasattr(resonance,'guiName'):
    del resonance.guiName

  if hasattr(resonance,'label'):
    del resonance.label
  
  makeResonanceGuiName(resonance)
      
  for contrib in resonance.peakDimContribs:
    contrib.peakDim.setAnnotation( makePeakDimAnnotation(contrib.peakDim) )

  if hasattr(resonance, 'onebond'):
    updateResonanceAnnotation(resonance.onebond)


def updateContribPeakDimAnnotation(contrib):
  """
  Updates annotation string for the peak dimension of a given
  contribution according to its assignment status

  .. describe:: Input
  
  Nmr.PeakDimContrib
  
  .. describe:: Output
  
  None
  """
  from ccpnmr.analysis.core.CouplingBasic import updateClusterAssignments

  peakDim = contrib.peakDim
  resonance = contrib.resonance
   
  #if contrib.peakDimComponent:
  for peakCluster in peakDim.peak.sortedPeakClusters():
    updateClusterAssignments(peakCluster, peakDim.peak)
    break # First only

  if not peakDim.isDeleted:
    makePeakDimAnnotation(peakDim)
  
  if not resonance.isDeleted:
    getBoundResonances(resonance, recalculate=True)
    updateResonShift(resonance, peakDim)

def updateResonanceSetAnnotation(resonanceSet):
  """
  Update the annotations for the peaks assigned to a resonanceSet's resonances

  .. describe:: Input
  
  Nmr.ResonanceSet
  
  .. describe:: Output
  
  None 
  """
  
  for resonance in resonanceSet.resonances:
    if not resonance.isDeleted: 
      getBoundResonances(resonance, recalculate=True, contribs=None)
      updateResonanceAnnotation(resonance)

def updateResidueAnnotation(residue):
  """
  Update the annotations for the peaks assigned to a residues resonances

  .. describe:: Input
  
  MolSystem.Residue
  
  .. describe:: Output
  
  None 
  """

  for atom in residue.atoms:
    if atom.atomSet:
      for resonanceSet in atom.atomSet.resonanceSets:
        updateResonanceSetAnnotation(resonanceSet)

  resonances = set()
  for residueProb in residue.residueProbs:
    resonances.update(residueProb.resonanceGroup.resonances)

  for spinSystem in residue.resonanceGroups:
    resonances.update(spinSystem.resonances)

  for resonance in resonances:
    updateResonanceAnnotation(resonance)

def updateResidueProbAnnotation(residueProb):
  """
  Update the annotations for the peaks assigned to a resonances
  within a spin system that has a tentative assignment.

  .. describe:: Input
  
  Nmr.ResidueProb
  
  .. describe:: Output
  
  None 
  """

  for resonance in residueProb.resonanceGroup.resonances:
    updateResonanceAnnotation(resonance)

def updateResonanceGroupAnnotation(spinSystem):
  """
  Update the annotations for the peaks assigned to a resonances
  within a spin system.

  .. describe:: Input
  
  Nmr.ResidueProb
  
  .. describe:: Output
  
  None 
  """

  for resonance in spinSystem.resonances:
    updateResonanceAnnotation(resonance)


def getBoundResonances(resonance, recalculate=False, contribs=None, doWarning=False,
                        recursiveCall=False):
  """
  Find all resonances that have a single bond connection to the input resonance
  Option to recalculate given assignment status (e.g. if something changes)
  Option to specify peakDimContribs to search

  .. describe:: Input
  
  Nmr.Resonance, Boolean, List of Nmr.PeakDimContribs
  
  .. describe:: Output
  
  List of Nmr.Resonances
  """
  
  from ccpnmr.analysis.core.MoleculeBasic import getBoundAtoms
  
  if (not recalculate) and resonance.covalentlyBound:
    return list(resonance.covalentlyBound)
    
  resonances = set() # Linked by bound atoms irrespective of spectra
  pairResonances = set() # prochiral or other pairs that can not be determined imemdiately
  resonanceSet   = resonance.resonanceSet
  
  funnyResonances = set()

  if resonanceSet:
    #residue  = resonanceSet.findFirstAtomSet().findFirstAtom().residue
    atomSets = resonanceSet.atomSets

    for atomSet in atomSets:
      #for atom in atomSet.atoms:
      atom = atomSet.findFirstAtom()
      
      for atom2 in getBoundAtoms(atom):
        atomSet2 = atom2.atomSet
      
        if atomSet2 and atomSet2.resonanceSets:
          
          usePaired = False
          if len(atomSets) > 1:
            chemAtomSet = atom2.chemAtom.chemAtomSet
            if chemAtomSet:
              usePaired = (chemAtomSet.isProchiral or 
                           (chemAtomSet.chemAtomSet and chemAtomSet.chemAtomSet.isProchiral))

          for resonanceSet2 in atomSet2.resonanceSets:
            for resonance2 in resonanceSet2.resonances:
              if resonance2 is resonance: # should not happen
                if resonance not in funnyResonances:
                  print 'WARNING: in getBoundResonances(): resonance %d tried to be linked to itself' % resonance.serial
                  funnyResonances.add(resonance)
              elif usePaired:
                pairResonances.add(resonance2)
              else:
                resonances.add(resonance2)
  
  if not contribs:
    contribs = resonance.peakDimContribs

  expResonances = set()
  foundBothPaired = False
  for contrib in contribs:
    peakDim      = contrib.peakDim
    expDimRef1   = peakDim.dataDimRef.expDimRef
    expTransfers = expDimRef1.expTransfers
    
    for expTransfer in expTransfers:
      if expTransfer.transferType in ('onebond','CP'):
        expDimRef2 = None

        for expDimRef in expTransfer.expDimRefs:
          if expDimRef is not expDimRef1:
            expDimRef2 = expDimRef
            break

        if expDimRef2:
          for peakDim2 in peakDim.peak.peakDims:
            if peakDim2.dataDimRef and (peakDim2.dataDimRef.expDimRef is expDimRef2):
              expBound = set()
 
              for contrib2 in peakDim2.peakDimContribs:
                if (not contrib.peakContribs) and (not contrib2.peakContribs):
                  resonance2 = contrib2.resonance
                  
                  if resonance is not resonance2:
                    expBound.add(resonance2)
 
                else:
                  for peakContrib in contrib.peakContribs:
                    if peakContrib in contrib2.peakContribs:
                      resonance2 = contrib2.resonance
                      
                      if resonance is not resonance2:
                        expBound.add(resonance2)

                      break

              if len(expBound) > 1:
                # Ambiguity
                for bound in expBound:
                  # Leave the covalently bound one
                  if bound in resonances:
                    break

                else:
                  aSet = set(x for x in expBound if x in resonance.covalentlyBound)
                  if aSet and aSet != pairResonances:
                    # Resonances found. Previously linked.
                    # Not the pairResonances. Use them
                    expResonances.update(aSet)
                    
                  else:
                    # check presence of prochiral pairs
                    ll = [x for x in pairResonances if x in expBound]
                    if len(pairResonances) == 2 and len(ll) == 2:
                      foundBothPaired= True
                    elif ll:
                      # found some prochiral pair resonances - use them
                      expResonances.update(ll)
              else:
                expResonances.update(expBound)
  
  if foundBothPaired and not [x for x in expResonances if x in pairResonances]:
    # particular special case. 
    # Resonnce is bound to both prochiral altrnatives but always as a pair.
    
    if recursiveCall:
      # This was called from elsewhere. We could resolve nothing, so send back to caller
      pass
    
    else:
      # call for sister resonances and see
      resons = resonanceSet.sortedResonances()
      newResonances = set()
      if len(resons)> 1:
        # there are sister resonances
        resons.remove(resonance)
        for reson in resons:
          boundResons = getBoundResonances(reson, recalculate=True, contribs=contribs,
                                           doWarning=False, recursiveCall=True)
          ll = [x for x in pairResonances if x not in boundResons]
          if not ll:
            # One sister was bound to both. Incorrect data. Bind to both here too
            newResonances.update(pairResonances)
            break
          elif len(ll) < len(pairResonances):
            # Some resonances were taken. Use the free ones.
            newResonances.update(ll)
      
      if newResonances:
        expResonances.update(newResonances)
      else:
        # No data anywhere to resolve which is which. Match on serials
        pairResonList = list(sorted(pairResonances, key=operator.attrgetter('serial')))
        rr = pairResonList[resonanceSet.sortedResonances().index(resonance)]
        expResonances.add(rr)
        
  
  resonances.update(expResonances)
 
  #if doWarning and (resonance.isotopeCode == '1H') and (len(resonances) > 1):
  #  pass
 
  if resonances:
    resonance.setCovalentlyBound(resonances)
  else:
    resonance.setCovalentlyBound([])
  
  return list(resonances)


def getOnebondResonance(resonance, isotopeCode=None):
  """
  Find any resonance that may have a single bond connection to the input resonance
  Option to specify the isotope type

  .. describe:: Input
  
  Nmr.Resonance, Nmr.Resonance.isotopeCode
  
  .. describe:: Output
  
  Nmr.Resonance
  """
  
  resonances = getBoundResonances(resonance)
  if resonances:
    if isotopeCode:
      for resonance1 in resonances:
        if resonance1.isotopeCode == isotopeCode:
          return resonance1
    
    else:
      return resonances[0]
  
  resonance2 = None
  
  for contrib in resonance.peakDimContribs:
    peakDim      = contrib.peakDim
    expDimRef1   = peakDim.dataDimRef.expDimRef
    expTransfers = expDimRef1.expTransfers

    for expTransfer in expTransfers:
      if expTransfer.transferType in ('onebond','CP'):
        expDimRef2 = None

        for expDimRef in expTransfer.expDimRefs:
          if expDimRef is not expDimRef1:
            expDimRef2 = expDimRef
            break

        if expDimRef2:
          if (not isotopeCode) or (isotopeCode in expDimRef2.isotopeCodes):
            for peakDim2 in peakDim.peak.peakDims:
              if peakDim2.dataDimRef and (peakDim2.dataDimRef.expDimRef is expDimRef2):
                for contrib2 in peakDim2.peakDimContribs:
                  if (not isotopeCode) or (contrib2.resonance.isotopeCode == isotopeCode):
                    resonance2 = contrib2.resonance
 
                break
              
    if resonance2:
      break   

  return resonance2
 
def getOnebondResonances(resonance, isotopeCode=None):
  """
  Find all resonances that may have a single bond connection to the input resonance
  Option to specify the isotope type

  .. describe:: Input

  Nmr.Resonance, Nmr.Resonance.isotopeCode

  .. describe:: Output

  List of Nmr.Resonances
  """

  resonances = getBoundResonances(resonance)
  if resonances:
    if isotopeCode:
      resonances = [res for res in resonances if res.isotopeCode == isotopeCode]

    return resonances

  resonances = set()
  for contrib in resonance.peakDimContribs:
    peakDim      = contrib.peakDim
    expDimRef1   = peakDim.dataDimRef.expDimRef
    expTransfers = expDimRef1.expTransfers

    for expTransfer in expTransfers:
      if expTransfer.transferType in ('onebond','CP'):
        expDimRef2 = None

        for expDimRef in expTransfer.expDimRefs:
          if expDimRef is not expDimRef1:
            expDimRef2 = expDimRef
            break

        if expDimRef2:
          if (not isotopeCode) or (isotopeCode in expDimRef2.isotopeCodes):
            for peakDim2 in peakDim.peak.peakDims:
              if peakDim2.dataDimRef and (peakDim2.dataDimRef.expDimRef is expDimRef2):
                for contrib2 in peakDim2.peakDimContribs:
                  if (not isotopeCode) or (contrib2.resonance.isotopeCode == isotopeCode):
                    resonances.add(contrib2.resonance)

                break

  return list(resonances)



def getOnebondFixedResonance(resonance, isotopeCode=None):
  """
  Find any Fixedesonance that may have a single bond connetion to the input resonance
  Option to specify the isotope type.
  FixedResonance version of getOnebondResonance

  .. describe:: Input
  
  Nmr.Resonance, Nmr.Resonance.isotopeCode
  
  .. describe:: Output
  
  Nmr.Rssonance
  """
  
  resonances = resonance.sortedCovalentlyBound()
  if resonances:
    if isotopeCode:
      for resonance1 in resonances:
        if resonance1.isotopeCode == isotopeCode:
          return resonance1
    
    else:
      return resonances[0]
      
  else:
    return None


def getAmbigProchiralLabel(resonance):
  """
  Deterimine if an ambigous prochiral resonance (non-stereospecifically assigned)
  Has an "a" label or a "b" label. "a" is reserved for the upfield proton and any
  other nulceus bound to it.

  .. describe:: Input
  
  Nmr.Resonance
  
  .. describe:: Output
  
  Character
  """

  letter = ''
  if hasattr(resonance, 'onebond'):
    del resonance.onebond

  resonanceSet = resonance.resonanceSet
  
  if resonanceSet:
    if resonance.isotopeCode == '1H':
      data = []
      for resonance2 in resonanceSet.sortedResonances():
        if resonance2.shifts:
          data.append( ('%f%d' % (resonance2.findFirstShift().value,resonance2.serial),resonance2) )
        else:
          data.append( (resonance2.serial,resonance2) )
 
      data.sort()
      resonances = [x[1] for x in data]
      i = resonances.index(resonance)
      letter = chr(ord('a')+i)

    else:
      resonance2 = getOnebondResonance(resonance, isotopeCode='1H') 
 
      if resonance2 and resonance2.resonanceSet and (len(resonance2.resonanceSet.atomSets) > 1):
        letter = getAmbigProchiralLabel(resonance2)
        resonance2.onebond = resonance
      
      elif (len(resonanceSet.resonances) > 1) and (len(resonanceSet.atomSets) > 1):
        for resonance2 in resonanceSet.resonances:
          if resonance2 is not resonance:
            resonance3 = getOnebondResonance(resonance2)
            if resonance3 and resonance3.resonanceSet and (len(resonance3.resonanceSet.atomSets) > 1):
              letter = 'b'
            break
             
      if not letter:
        data = []
        for resonance2 in resonanceSet.resonances:
          if resonance2.shifts:
            data.append( (resonance2.findFirstShift().value,resonance2) )
          else:
            data.append( (resonance2.serial,resonance2) )
 
        data.sort()
        resonances = [x[1] for x in data]
        i = resonances.index(resonance)
        letter = chr(ord('a')+i)       

  #keyword = 'ambigProchiralLabel'
  #app     = 'Analysis'
  #appData = resonance.findFirstApplicationData(application=app, keyword=keyword)
  #
  #if appData and (appData.value != letter):
  #  appData.delete()
  #  appData = None
  #    
  #if not appData:
  #  AppDataString(resonance,application=app,keyword=keyword, value=letter)
    
  return letter


def swapProchiralResonance(resonance, makeAmbiguous=False):
  """
  If a resonance is assigned to a prociral centre in a stereospecific manner, change
  the assignment to the other prochiral atom set. Automatically swaps over any other
  resonance that may be assigned to the other prochiral option.
  If the makeAmbiguous option is True then the unambiougously assigned resonances
  are changed to ambigous assignments by assigning resonances to both prochiral options. 

  .. describe:: Input
  
  Nmr.Resonance, Boolean
  
  .. describe:: Output
  
  None
  """
  
  from ccpnmr.analysis.core.MoleculeBasic import areAtomsBound

  if resonance.resonanceSet and (len(resonance.resonanceSet.atomSets) == 1):
    atomSet     = resonance.resonanceSet.findFirstAtomSet()
    atom        = atomSet.findFirstAtom()
    chemAtom    = atom.chemAtom
    chemAtomSet = chemAtom.chemAtomSet
    chemAtom0   = None
    
    if chemAtomSet:
      if chemAtomSet.isEquivalent:
        for chemAtomSet0 in chemAtom.chemComp.findAllChemAtomSets(isProchiral=True):
          chemAtomSets = list(chemAtomSet0.chemAtomSets)
          if chemAtomSet in chemAtomSet0.chemAtomSets:
            chemAtomSets.remove(chemAtomSet)
            chemAtom0 = chemAtomSets[0].findFirstChemAtom()
            break
    
      elif chemAtomSet.isProchiral:
        for chemAtom1 in chemAtomSet.chemAtoms:
          if chemAtom1 is not chemAtom:
            chemAtom0 = chemAtom1
            break

    if chemAtom0:
      atom0 = atom.residue.findFirstAtom(name=chemAtom0.name)
      if atom0 and atom0.atomSet:
      
        resonances0 = list(resonance.resonanceSet.resonances)
        resonances1 = []
      
        for resonanceSet in atom0.atomSet.resonanceSets:
          if len(resonanceSet.atomSets) == 1:
            resonances1 = list(resonanceSet.resonances)
            break
      
        resonance.resonanceSet.delete()
        if resonances1:
          resonances1[0].resonanceSet.delete()
 
        if makeAmbiguous:
          for resonance1 in resonances1:
            assignAtomsToRes([atomSet,atom0.atomSet],resonance1)

          for resonance0 in resonances0:
            assignAtomsToRes([atomSet,atom0.atomSet],resonance0)

        else:
 
          for resonance1 in resonances1:
            assignAtomsToRes([atomSet, ],resonance1)

          for resonance0 in resonances0:
            assignAtomsToRes([atom0.atomSet, ],resonance0)
        
        if chemAtomSet.isEquivalent:
          for resonance0 in resonances0:
            resonancesX = getBoundResonances(resonance0, recalculate=True)
            
            for bound in resonancesX:
              if bound.resonanceSet:
                atomB = bound.resonanceSet.findFirstAtomSet().findFirstAtom()
 
                if areAtomsBound(atom, atomB):
                  swapProchiralResonance(bound, makeAmbiguous)
 
                break
        
          else:
            for resonance1 in resonances1:
              resonancesX = getBoundResonances(resonance0, recalculate=True)
            
              for bound in resonancesX:
                if bound.resonanceSet:
                  atomB = bound.resonanceSet.findFirstAtomSet().findFirstAtom()
 
                  if areAtomsBound(atom0, atomB):
                    swapProchiralResonance(bound, makeAmbiguous)
 
                  break
 
 
def getAtomSetShifts(atomSet, shiftList=None):
  """
  Gives the shifts that an atom set may be assigned to in a specified shift list
  or the first shift list in the project if none is specified.

  .. describe:: Input
  
  Nmr.AtomSet, Nmr.ShiftList
  
  .. describe:: Output
  
  List of Nmr.Shifts
  """

  if not shiftList:

    molSystem = atomSet.findFirstAtom().residue.chain.molSystem
    for experiment in molSystem.root.currentNmrProject.sortedExperiments():
      if experiment.shiftList:
        shiftList = experiment.shiftList
        break

    if not shiftList:
      shiftList = atomSet.nmrProject.findFirstMeasurementList(className='ShiftList')
  
  if not shiftList:
    return []

  shifts = []
  for resonanceSet in atomSet.resonanceSets:
    for resonance in resonanceSet.resonances:
      shift = resonance.findFirstShift(parentList = shiftList)
      if shift:
        shifts.append(shift)
  
  return shifts

def getAliasedPeakDimPositions(peakDim, shiftRanges=None, ppm=None):
  """
  Give all the aliased/unaliased positions of a peakDim either in a
  specified shift range or the full range for the dimension type.
  Units for the shift ranges are ppm. Note this function uses the
  peakDim.realValue, i.e. center of couplings for main assignment.
  Use the PeakBasic version if the actual extremum location should
  be used. Optional ppm argument if main position is not the one
  to be considered, e.g. in reduced dimensionality or MQ.

  .. describe:: Input
  
  Nmr.PeakDim, List of (Tuples of Floats (MinShift,MaxShift) )
  
  .. describe:: Output
  
  List of Floats (Nmr.PeakDim.positions)
  """

  if not shiftRanges:
    shiftRanges = [getPeakDimFullShiftRange(peakDim), ]

  positions  = []
  dataDimRef = peakDim.dataDimRef
  
  if dataDimRef:
    sw = dataDimRef.dataDim.numPointsOrig
    
    if ppm is None:
      ppm = peakDim.realValue
    
    points = points0 = ppm2pnt(ppm, dataDimRef)
    while isShiftInRange(ppm, shiftRanges):
      positions.append( points )
      points -= sw
      ppm = pnt2ppm(points, dataDimRef)

    points = points0+sw
    while isShiftInRange(ppm, shiftRanges):
      positions.append( points )
      points += sw
      ppm = pnt2ppm(points, dataDimRef)
  
  else:
    positions.append(peakDim.position)
    

  return positions

def areResonancesProchiral(resonanceA, resonanceB):
  """
  Determine if a pair of resonances are assigned to a
  potentially prochiral pair of atoms.

  .. describe:: Input
  
  Nmr.Resonance, Nmr.Resonance
  
  .. describe:: Output
  
  Boolean
  """

  if resonanceA.resonanceSet:
    if resonanceB.resonanceSet:
      atomA = resonanceA.resonanceSet.findFirstAtomSet().findFirstAtom()
      atomB = resonanceB.resonanceSet.findFirstAtomSet().findFirstAtom()
    
      residueA = atomA.residue
      residueB = atomB.residue
    
      if residueA is residueB:
        chemAtomSetA = atomA.chemAtom.chemAtomSet
        chemAtomSetB = atomB.chemAtom.chemAtomSet

        if chemAtomSetA and (chemAtomSetA is chemAtomSetB):
          return True
      
  return False


def isShiftInRange(shiftValue, shiftRanges):
  """
  Determine whether a chemical shift value is contained within
  one of the input ranges.

  .. describe:: Input
  
  Float (Nmr.Shift.value), List of Tuples of Floats (min, max range bounds)
  
  .. describe:: Output
  
  Boolean
  """

  for (minVal, maxVal) in shiftRanges:
    if shiftValue >= minVal:
      if shiftValue <= maxVal:
        return True
          
  return False
  
def isChainAssigned(chain):
  """
  Determines if an NMR chain has any atoms in its
  residues that are assigned to resonances

  .. describe:: Input
  
  Nmr.Chain
  
  .. describe:: Output
  
  Boolean
  """

  for residue in chain.residues:
    for atom in residue.atoms:
      if atom.atomSet:
        if atom.atomSet.resonanceSets:
          return True
      
  return False

def isResidueAssigned(residue):
  """
  Determines if an NMR residue has any atoms that are assigned to resonances

  .. describe:: Input
  
  MolSystem.Residue
  
  .. describe:: Output
  
  Boolean
  """

  for atom in residue.atoms:
    if atom.atomSet:
      if atom.atomSet.resonanceSets:
        return True
  
  return False

def isAtomAssigned(atom, toPeaks=False):
  """
  Determines if an atom is assigned to any resonances, with
  an optional argument to check if the resonances are linked
  to peaks.

  .. describe:: Input
  
  MolSystem.Atom
  
  .. describe:: Output
  
  Boolean
  """

  atomSet = atom.atomSet
  if atomSet and atomSet.resonanceSets:
    if toPeaks:
      for resonanceSet in atomSet.resonanceSets:
        for resonance in resonanceSet.resonances:
          if resonance.peakDimContribs:
            return True
    else:
      return True
  
  return False

def findFirstAtomShiftInShiftList(atom, shiftList):
  """
  Finds the first shift in a shiftList for the specified atom.

  .. describe:: Input
  
  MolSystem.Atom, Nmr.ShiftList
  
  .. describe:: Output
  
  Nmr.Shift
  """

  atomSet = atom.atomSet
  if atomSet:
    for resonanceSet in atomSet.resonanceSets:
      for resonance in resonanceSet.resonances:
        shift = shiftList.findFirstMeasurement(resonance=resonance)
        if shift:
          return shift

  return None

def isPeakAssigned(peak, fully=True):
  """
  Determines if a peak is assigned, ether fully or partially

  .. describe:: Input
  
  Peak, Boolean
  
  .. describe:: Output
  
  Boolean
  """

  n = 0
  for peakDim in peak.peakDims:
    if len(peakDim.peakDimContribs) > 0:
      n +=1
      
  if n == len(peak.peakDims):
    return True
    
  elif n > 0:
    if fully:
      return False
    else:
      return True
      
  else:
    return False

def initResonance(resonance, doMerge=True):
  """
  Initialise a resonance. If the resonance is assigned it is
  placed in a spin system if appropriate and its assignName is set.
  Optional argument to set whether merging spin systems is allowed.
  Typically this is not allowed upon load, but is allowed when
  assigning within Analysis.

  .. describe:: Input
  
  Nmr.Resonance, Boolean
  
  .. describe:: Output
  
  None
  """
  from ccpnmr.analysis.core.MoleculeBasic import DEFAULT_ISOTOPES
  
  if resonance.isDeleted:
    return
  
  #print 'initResonance', makeResonanceGuiName(resonance)
  resonanceSet = resonance.resonanceSet
  
  if resonance.name == 'r%d' % resonance.serial:
    resonance.setName(None)
  
  if resonanceSet:
    atomSets   = tuple(resonanceSet.atomSets)
    atom       = atomSets[0].findFirstAtom()
    element    = atom.chemAtom.elementSymbol
    if resonance.isotopeCode[-len(element):] != element:
      msg = 'Resonance %d isotope-assignment mismatch: Resetting isotope' 
      print msg % resonance.serial
      
      resonance.isotopeCode = DEFAULT_ISOTOPES.get(element, 'unknown')
      getBoundResonances(resonance, recalculate=True, contribs=None)   
      
    spinSystem = findSpinSystem(resonance)
    if spinSystem:
      if resonance.resonanceGroup: # same as the spin system
        residue = atom.residue
        
        if spinSystem.residue is not residue:
          #assignSpinSystemResidue(spinSystem,residue)
              
          if doMerge:
            ccpCode = residue.molResidue.ccpCode
            spinSystems = list(resonance.nmrProject.findAllResonanceGroups(residue=residue))
            N = len(spinSystems)
            msg = 'There are %d separate %d%s spin systems. Merge together?'
            if (N > 1) and showOkCancel('Confirm', msg % (N,residue.seqCode,ccpCode)):
              assignSpinSystemResidue(spinSystem,residue)  
                
              for spinSystem1 in spinSystems[1:]:
                mergeSpinSystems(spinSystem1,spinSystems[0])
                
            else:
              assignSpinSystemResidue(spinSystem,residue)    
              
          else:  # 24 Nov 2017: with the else, ended up with problem when had mixed DNA/RNA terms
          #elif spinSystem.residue:
            try:  # 27 Nov 2017: make it a bit less dramatic
              assignSpinSystemResidue(spinSystem,residue)
            except Implementation.ApiError, e:
              pass # this can happen when have mixed molType in spinSystem, so just give up
                           
      else:
        addSpinSystemResonance(spinSystem, resonance)

    assignResonanceType(resonance, atomSets)
    
  else:
    updateResonanceAnnotation(resonance)

  if not resonance.shifts:
    if resonance.peakDimContribs:
      for contrib in resonance.peakDimContribs:
        peakDim = contrib.peakDim
        experiment = peakDim.peak.peakList.dataSource.experiment
        
        if experiment.shiftList is None:
          shiftList = resonance.nmrProject.findFirstMeasurementList(className='ShiftList')
          if shiftList is None:
            shiftList = resonance.nmrProject.newShiftList(unit='ppm')
          experiment.setShiftList( shiftList )
        
        updateResonShift(resonance,peakDim)
        
 

def getFixedResonanceName(resonance):
  """
  Give the name for a fixedResonance, based upon a normal resonance

  .. describe:: Input
  
  Nmr.Resonance
  
  .. describe:: Output
  
  String
  """
  
  resonanceSet = resonance.resonanceSet
  if resonanceSet:
    atomSets = tuple(resonanceSet.atomSets)

    if len(atomSets) > 3: # Whole residue selections
      return '*'

    name = makeAtomSetsGuiName(atomSets)
    if len(atomSets) > 1:
      resonances = resonanceSet.sortedResonances()
      i = resonances.index(resonance)
      name = name[:-1] + chr(ord('a')+i)
      if len(atomSets[0].atoms) > 1:
        name = name + '*'

  else:
    name = '[%d]' % resonance.serial

  return name


def getResonanceName(resonance):
  """
  Generate a name for a resonance or fixedResonance based upon its assignment

  .. describe:: Input
  
  Nmr.Resonance or NmrConstraint.FixedResonance
  
  .. describe:: Output
  
  String
  """

  if resonance.className == 'FixedResonance':
    resonance2 = resonance.getResonance()
    if resonance2 is None:
      return getFixedResonanceName(resonance)
    
    else:
      resonanceSet = resonance.resonanceSet
      resonanceSet2 = resonance2.resonanceSet
    
      if not resonanceSet2:
        return getFixedResonanceName(resonance)
      
      elif not resonanceSet:
        resonance = resonance2
    
      else:
        atoms = [aSet.atoms for aSet in resonanceSet.sortedAtomSets()]
        atoms2 = [aSet.atoms for aSet in resonanceSet2.sortedAtomSets()]
          
        if atoms == atoms2:
          resonance = resonance2
        else:
          return getFixedResonanceName(resonance)

  if hasattr(resonance, 'label'):
    return resonance.label

  resonanceSet = resonance.resonanceSet
  
  if resonanceSet:
    atomSets = tuple(resonanceSet.atomSets)
    name = makeAtomSetsGuiName(atomSets)
    if len(atomSets) > 1:
      name = name[:-1] + getAmbigProchiralLabel(resonance)
      if len(atomSets[0].atoms) > 1:
        name = name + '*'
    
  elif resonance.assignNames:
    assignNames = tuple(resonance.assignNames)
    name = '[%d]' % resonance.serial

    N = len(assignNames)
    if N > 1:
      for i in range(N):
        name = assignNames[i] + name
        if i < (N-1):
          name = '|' + name
    else:
      name = resonance.assignNames[0] + name
  
  elif resonance.name and (resonance.name != 'r%d' % resonance.serial):
    name = '[%d:%s]' % (resonance.serial,resonance.name)

  else:
    name = '[%d]' % resonance.serial

  resonance.label = name
              
  return name


def newResonance(project, isotopeCode='unknown'):
  """
  Make a new resonance in a given project. Sets initial isotopeCode.
  Uses the current NMR project.

  .. describe:: Input
  
  Project, Word (Resonance.isotopeCode)
  
  .. describe:: Output
  
  Nmr.Resonance
  """

  resonance = project.currentNmrProject.newResonance(isotopeCode=isotopeCode)
  
  return resonance


def makePeakDimAnnotation(peakDim):
  """
  Makes an annotation string for a given peak dimension according to its assignment status

  .. describe:: Input
  
  Nmr.PeakDim
  
  .. describe:: Output
  
  String (Nmr.PeakDim.annotation)
  """
  
  from ccpnmr.analysis.core.PeakBasic import makePeakAnnotation
  from ccpnmr.analysis.core.CouplingBasic import updateClusterAnnotation

  peakDim.setAnnotation(None)

  chainCode, seqId, name = getPeakDimAtomTuple(peakDim)
    
  annotation = '%s%s%s' % (chainCode,seqId,name)
  
  if len(annotation) > 80:
    annotation = '**TooLong!**'
  
  peakDim.annotation = annotation or None
  
  peak = peakDim.peak
  if peak.root.currentAnalysisProject:
    makePeakAnnotation(peak)

  for peakCluster in peak.peakClusters:
    updateClusterAnnotation(peakCluster)

  return annotation

def getPeakDimAtomTuple(peakDim):
  """
  Give a tupe of string identifiers for a peakDim indicating, chain, residue and atomic assignment

  .. describe:: Input
  
  Nmr.PeakDim
  
  .. describe:: Output
  
  3-Tuple of Words
  """
  from ccpnmr.analysis.core.MoleculeBasic import getResidueCode


  contribs = peakDim.peakDimContribs
  Ncontribs = len(contribs)

  if Ncontribs < 1:
    if not peakDim.dataDimRef:
      dataDim = peakDim.dataDim
      if dataDim.className != 'FreqDataDim':
        index = int(peakDim.position or 1) -1
        value = dataDim.pointValues[index]
        return ('','','%s' % value)

    else:
      return ('','','')
  
  analysisProject = peakDim.root.currentAnalysisProject
  if analysisProject:
    sysAnno = analysisProject.doSpinSystemAnnotations
    assAnno = analysisProject.doAssignmentAnnotations
    chnAnno = analysisProject.doChainAnnotations
    oneAnno = False
    numAnno = True
  else:
    sysAnno = True
    assAnno = True
    chnAnno = False
    oneAnno = False
    numAnno = True

    
  chemAtomSetDict = {}

  residues = []
  atoms  = []
  chains = []
  for contrib in contribs:
    if isinstance(contrib, Nmr.PeakDimContribN):
      continue
  
    resonance = contrib.resonance
    resonanceGroup = resonance.resonanceGroup
    resonanceSet = resonance.resonanceSet
    myChain = ''
    myResidue = ''
    myAtom  = makeResonanceGuiName( resonance, fullName=False )

    if resonanceSet:
      atomSet = resonanceSet.findFirstAtomSet()
      atom = atomSet.findFirstAtom()
      residue = atom.residue
      chemAtomSet = atom.chemAtom.chemAtomSet
      
      key = (residue, chemAtomSet)
      if chemAtomSet and chemAtomSetDict.has_key(key):
        index = chemAtomSetDict[key]
        if index < len(atoms):
          if len(atomSet.atoms) > 1:
            atoms[index] = myAtom[:-2] + '*'
          else:
            atoms[index] = myAtom[:-1] + '*'

        continue
      
      chemAtomSetDict[key] = len(atoms)
      
      if assAnno:
        if chnAnno:
          myChain = residue.chain.code
        
        if oneAnno:
          resCode = residue.chemCompVar.chemComp.code1Letter or '?'
        else:
          resCode = getResidueCode(residue)
          
        if numAnno:
          myResidue = str(residue.seqCode) + resCode
        
        else:  
          myResidue = resCode + str(residue.seqCode)  

      elif sysAnno and  resonanceGroup:
        residue = resonanceGroup.residue

        if residue:
          if chnAnno:
            myChain = residue.chain.code
         
          if oneAnno:
            resCode = residue.chemCompVar.chemComp.code1Letter or '?'
          else:
            resCode = getResidueCode(residue)  
                    
          if numAnno:
            myResidue = '%d%s' % (residue.seqCode, resCode)
          else:
            myResidue = '%s%d' % (resCode, residue.seqCode, )
           
        
        elif resonanceGroup.ccpCode:
          myResidue = '{%d}%s' % (resonanceGroup.serial, resonanceGroup.ccpCode)
        
        elif resonanceGroup.name:
          myResidue = '{%d:%s}' % (resonanceGroup.serial, resonanceGroup.name)
        
        else:
          myResidue = '{%d}' % resonanceGroup.serial

    elif resonanceGroup:       
      if sysAnno:
        
        residue = resonanceGroup.residue
        if residue:
          myResidue = '%d' % residue.seqCode
          
          if chnAnno:
            myChain = residue.chain.code
            
          if oneAnno:
            resCode = residue.chemCompVar.chemComp.code1Letter or '?'
          else:
            resCode = getResidueCode(residue)  
             
          if numAnno:
            myResidue += resCode
          else:
            myResidue = resCode + myResidue
        
        elif resonanceGroup.residueProbs:
          resTexts = []
          resSeqs = []
          resCodes = set()
 
          for residueProb in resonanceGroup.residueProbs:
            if not residueProb.weight:
              continue
            
            residue = residueProb.possibility
            seq = residue.seqCode
         
            if oneAnno:
              resCode = residue.chemCompVar.chemComp.code1Letter or '?'
            else:
              resCode = getResidueCode(residue)  
                    
            if numAnno:
              resText = '%d?%s' % (seq, resCode)
            else:
              resText = '%s%d?' % (resCode, seq)

            resTexts.append(resText)
            resSeqs.append('%d?' % seq)
            resCodes.add(resCode)
 
          if len(resCodes) == 1:
            myResidue = '/'.join(resSeqs) + resCodes.pop()
          else:
            myResidue = '/'.join(resTexts)
        
        elif resonanceGroup.ccpCode:
          myResidue = '{%d}' % resonanceGroup.serial
          
          if numAnno:
            myResidue += resonanceGroup.ccpCode
          else:
            myResidue = resonanceGroup.ccpCode + myResidue
        
        elif resonanceGroup.name:
          myResidue = '{%d:%s}' % (resonanceGroup.serial, resonanceGroup.name)
        
        else:
          myResidue = '{%d}' % resonanceGroup.serial
    
    chains.append( myChain )
    residues.append( myResidue )
    atoms.append ( myAtom  )

  # we always have contribs & atoms at this point
  # - otherwise would return before now
  
  if not chains:
    return ('','','')
  
  chainCode = chains[-1]
  residue  = residues[-1]

  if len(atoms) == 1:
    atom = atoms[0]
    
  else:  
    atom      = atoms[-1]
    atomLen   = len(atom)
    
    if chains[0] != chainCode:
      if len(chains) == 2:
        chainCode = '(%s)' % ('/'.join([c or '?' for c in chains]))
      else:
        chainCode = '?'
 
    if residues[0] != residue:
      if len(residues) == 2:
        residue = '(%s)' % ('/'.join([s or '?' for s in residues]))
      else:
        residue = '{*}'

    for i in range(len(atoms[:-1])):
      if len(atoms[i]) > atomLen:
        atomLen = len(atoms[i])
        atom = atoms[i]
 
    common = 0
    for i in range(atomLen):
      for atom0 in atoms:
        if i >= len(atom0) or ( atom0[i] != atom[i]):
          break
        if atom0[i] == '[':
          common += 1
          break
 
      else:
        common += 1
        continue
 
      break
 
    atom = atoms[0][:common]
    if peakDim.dataDimRef and peakDim.dataDimRef.expDimRef.measurementType == 'MQShift':
      atom = '(%s)' % ( '+'.join(atoms) )
    
    elif atom == atoms[0]:
      pass
      
    elif atom and (atom[-1] == '['):
      common -= 1
      atom = atom[:-1]
      atom += '(%s)' % ('/'.join([x[common:] for x in atoms]) )
    
    elif (common == atomLen-1):
      atom += '*' 
 
    else:
      atom = '(%s)' % ('/'.join(atoms) )

  if not numAnno:
    if chainCode:
      chainCode += ':'
   
  else:
    if oneAnno:
      residue += ' '

  return chainCode, residue, atom 
  
def makeAtomSetsGuiName(atomSets):
  """
  Give the name of a given set of atoms for GUIs

  .. describe:: Input
  
  List of Nmr.AtomSets
  
  .. describe:: Output
  
  Word
  """

  #from ccpnmr.analysis.core.MoleculeBasic import unicodeGreek

  atomNames = []
  for atomSet in atomSets:
    for atom in atomSet.atoms:
      atomNames.append(atom.name)

  tryName = atomNames[0]
  end = ''
  for atomName in atomNames[1:]:
    size = len(atomName)
    i = 0
    n = min(size, len(tryName))
    while i < n and atomName[i] == tryName[i]:
      i +=1
    
    tryName = tryName[0:i]
    end = '*'
        
  name = tryName[:1] + tryName[1:].lower() + end
  
  #if atom0.residue.molType == 'protein':
  #  i = len(atom0.chemAtom.elementSymbol)
  #  
  #  if i < len(name):
  #    letter = name[i]
  #    
  #    if letter.islower():
  #      name = name[:i] + unicodeGreek.get(letter,letter) + name[i+1:]
  #      name = name.encode('utf-8')
  
  return name

            
def makeResonanceGuiName(resonance, fullName=True, doAtoms=None):
  """
  Give the name of a resonance for GUIs. Either full identifier or
  just the atom type name (if any), Option to force the display
  of atom assignements set doAtoms to True/False

  .. describe:: Input
  
  Nmr.Resonance, Boolean, Boolean
  
  .. describe:: Output
  
  Word
  """

  molSysAnno = False
  assAnno    = True
  chnAnno    = False
  
  analysisProject = resonance.root.currentAnalysisProject
  if analysisProject:
    chnAnno    = analysisProject.doChainAnnotations
    molSysAnno = analysisProject.doMolSysAnnotations
    assAnno    = analysisProject.doAssignmentAnnotations
  
  if doAtoms is not None:
    assAnno = doAtoms
    
  if not fullName:
    if not assAnno:
      return '[%d]' % resonance.serial
    else:
      if hasattr(resonance, 'label'):
        return resonance.label
      else:  
        return getResonanceName(resonance)

  if hasattr(resonance, 'guiName'):
    return resonance.guiName

  (molSystem, chain, residue, name) = getResonanceAtomTuple(resonance)

  if not assAnno:
    name = '[%d]' % resonance.serial

  if molSysAnno:
    molSystem = molSystem+':'
  else:
    molSystem = ''
  if not chnAnno:
    chain = ''
  
  guiName = '%s%s%s%s' % (molSystem,chain,residue,name)
  if resonance.details:
    guiName = '%s:%s' % (guiName, resonance.details[:4])
  resonance.guiName = guiName
  resonance.label = name
  
  return guiName
  
def getResonanceAtomTuple(resonance):
  """
  Give a tupe of string identifiers for a resonance
  indicating, chain, residue and atomic assignment

  .. describe:: Input
  
  Nmr.Resonance
  
  .. describe:: Output
  
  Tuple of Words (MolSystem.Chain.code,
  MolSystem.Residue identifier, Nmr.Resonance identifier)
  """
  from ccpnmr.analysis.core.MoleculeBasic import getResidueCode
  
  molSystemCode = ''
  chainCode = ''
  res   = ''
  name  = getResonanceName(resonance)

  if resonance.resonanceSet:
    residue = getResonanceResidue(resonance)
    chain = residue.chain
    chainCode = chain.code
    molSystemCode = chain.molSystem.code
    res = str(residue.seqCode) + getResidueCode(residue)
    
  elif hasattr(resonance, 'resonanceGroup'): # Not FixedResonance
    spinSystem = resonance.resonanceGroup
    
    if spinSystem:
      residue = spinSystem.residue
      ccpCode = spinSystem.ccpCode
      residueProbs = [rp for rp in spinSystem.residueProbs if rp.weight > 0.0]
 
      if residue:
        sysChain = residue.chain
        chainCode = sysChain.code
        molSystemCode = sysChain.molSystem.code
        res = str(residue.seqCode) + getResidueCode(residue)

      elif residueProbs:
        resTexts = []
        resSeqs = []
        resCodes = set()
        
        for residueProb in residueProbs:
          residue = residueProb.possibility
          seq = residue.seqCode
          code = getResidueCode(residue)
          
          resTexts.append('%d?%s' % (seq,code))
          resSeqs.append('%d?' % seq)
          resCodes.add(code)
        
        if len(resCodes) == 1:
          res = '/'.join(resSeqs) + resCodes.pop()
        else:
          res = '/'.join(resTexts)
        
      elif ccpCode:
        res  = '{%d}%s' % (spinSystem.serial,ccpCode)

      elif spinSystem.name:
        res  = '{%d:%s}' % (spinSystem.serial,spinSystem.name)
 
      else:
        res  = '{%d}' % (spinSystem.serial)
   
  return (molSystemCode, chainCode, res, name)


def getDataDimFullShiftRange(dataDim):
  """
  Give the min and max possible chem shifts for a spectrums data dim
  based on set min/max epDim frequencies or full spec with

  .. describe:: Input
  
  Nmr.FreqDataDim
  
  .. describe:: Output
  
  List of Floats 
  """

  from ccpnmr.analysis.core.ExperimentBasic import getPrimaryDataDimRef

  expDimRef = dataDim.expDim.findFirstExpDimRef()
  shiftList = expDimRef.expDim.experiment.shiftList
  unit      = shiftList.unit

  if expDimRef.minAliasedFreq is None:
    if unit == 'point':
      minShift = dataDim.numPointsOrig
    else:
      dataDimRef = getPrimaryDataDimRef(dataDim)
      minShift = unit_converter[('point',unit)](dataDim.numPointsOrig,dataDimRef)

  else:
    minShift = expDimRef.minAliasedFreq

  if expDimRef.maxAliasedFreq is None:
    if unit == 'point':
      maxShift = 0
    else:
      dataDimRef = getPrimaryDataDimRef(dataDim)
      maxShift = unit_converter[('point',unit)](0,dataDimRef)

  else:
    maxShift = expDimRef.maxAliasedFreq
  
  shiftRange = [minShift,maxShift]
  shiftRange.sort()
  
  return shiftRange


def getPeakDimFullShiftRange(peakDim):
  """
  Give the min and max possible chem shifts for a peak dim
  based on set min/max epDim frequencies or full spec with

  .. describe:: Input
  
  Nmr.PeakDim
  
  .. describe:: Output
  
  List of Floats 
  """


  dataDimRef = peakDim.dataDimRef
  dataDim    = peakDim.dataDim
  if not dataDimRef:
    values = dataDim.pointValues
    return [ min(values), max(values) ]
  
  expDimRef = dataDim.expDim.findFirstExpDimRef()
  shiftList = expDimRef.expDim.experiment.shiftList
  unit      = shiftList.unit

  if expDimRef.minAliasedFreq is None:
    if unit == 'point':
      minShift = dataDim.numPointsOrig
    else:
      minShift = unit_converter[('point',unit)](dataDim.numPointsOrig,dataDimRef)

  else:
    minShift = expDimRef.minAliasedFreq

  if expDimRef.maxAliasedFreq is None:
    if unit == 'point':
      maxShift = 0
    else:
      maxShift = unit_converter[('point',unit)](0,dataDimRef)

  else:
    maxShift = expDimRef.maxAliasedFreq
  
  shiftRange = [minShift,maxShift]
  shiftRange.sort()
  
  return shiftRange

def findMatchingPeakDimShifts(peakDim, shiftRanges=None, tolerance=None,
                              aliasing=True, findAssigned=False, ppm=None):
  """
  For peakDim give the NMR shifts within a tolerance
  (which may be specified) of its position in a given set of ranges
  with options if the peakDim position might be aliased
  and to find only shifts with resonances assigned to atoms.
  Obeys molSystem assignments for atom linked peaks - these
  must match the peakDim's experiment molSystems.
  Option to pas in an alternative ppm value if the main peakDim value
  is not to be used (e.g. reduced dimensionality MQ etc).

  .. describe:: Input
  
  Nmr.PeakDim, List of Tuples (Float, Float) Float, Boolean, Boolean, Float
  
  .. describe:: Output
  
  List of Nmr.Shifts
  """
  
  shifts = []
  dataDimRef = peakDim.dataDimRef

  if not dataDimRef:
    return shifts
  
  if ppm is None:
    ppm = peakDim.realValue

    # NOTE:ED - in some cases ppm can still be None
    if ppm is None:
      showWarning('peakDim Failure',
                  'peakDim contains no realValue - possibly wrong experimentType')
      return shifts

  unit = dataDimRef.expDimRef.unit
  dataDim = dataDimRef.dataDim
  peakDimPos = unit_converter[(unit,'point')](ppm, dataDimRef)

  if aliasing:
    positions = getAliasedPeakDimPositions(peakDim, shiftRanges, ppm)
    if peakDimPos not in positions:
      positions.append(peakDimPos)
  else:
    positions  = [peakDimPos, ]
   
  for position in positions:
    shifts.extend(findMatchingShifts(dataDimRef, position,
                                     tolerance=tolerance,
                                     findAssigned=findAssigned))
  
  molSystems = peakDim.peak.peakList.dataSource.experiment.molSystems
  
  if molSystems:
    outShifts = []
    for shift in shifts:
      resonanceSet = shift.resonance.resonanceSet
 
      if resonanceSet:
        atom = resonanceSet.findFirstAtomSet().findFirstAtom()
        molSystem = atom.topObject
 
        if molSystem not in molSystems:
          continue
 
      outShifts.append(shift)
      
  else:
    outShifts = shifts
  
  return outShifts 

def findMatchingShifts(dataDimRef,position,tolerance=None,findAssigned=False):
  """
  For a dataDimRef give the NMR shifts within a tolerance
  (which may be specified) to a given position with an option
  to find only shifts with resonances assigned to atoms

  .. describe:: Input
  
  Nmr.DataDimRef, Float, Float, Boolean
  
  .. describe:: Output
  
  List of Nmr.Shifts
  """


  shiftList = dataDimRef.expDimRef.expDim.experiment.shiftList
  unit      = shiftList.unit
  isotopes  = dataDimRef.expDimRef.isotopeCodes
    
  dataDim   = dataDimRef.dataDim
  if not tolerance:
    tolerance = getAnalysisDataDim(dataDim).assignTolerance
  
  if unit != 'point':
    position = unit_converter[('point',unit)](position,dataDimRef)

  possibilities = []
  
  if not shiftList.measurements:
    return possibilities
  
  if not hasattr(shiftList, 'quickShifts'):
    setQuickShiftList(shiftList.findFirstMeasurement())
 
  # WARNING: this key function has to be same as in setQuickShiftList
  keyMin = int(round(10*(position-tolerance)))
  keyMax = int(round(10*(position+tolerance)))
  keyRange = range(keyMin, keyMax+1)
    
  for isotope in isotopes:
    quickShiftDict = shiftList.quickShifts.get(isotope, {})
    if quickShiftDict:
      for key in keyRange:
        shifts = quickShiftDict.get(key)
        if shifts is not None:
          for shift in shifts:
            if findAssigned and (not shift.resonance.resonanceSet):
              continue
 
            tryResonance = shift.resonance
            if tryResonance.isDeleted:
              continue
  
            difference = abs(position - shift.value)
            if difference < tolerance:
              possibilities.append(shift)
 
  return possibilities

def setAssignmentMolSystem(residue,peakDim=None,resonance=None):
  """
  Sets the MolSystem.Molsystem for all experiments with assignments to a given
  resonance and/or the experiment for a given peakDim
  Displays a warning if an assignment is made outside an experiment's
  list of molSystems.

  .. describe:: Input
  
  MolSystem.Residue, Nmr.PeakDim, Nmr.Resonance
  
  .. describe:: Output
  
  Boolean
  """

  molSystem    = residue.chain.molSystem
  peakDims     = []
  experiments  = {}
  experiments2 = {}
  if peakDim:
    peakDims = [peakDim]
  
  if resonance:
    for peakDimContrib in resonance.peakDimContribs:
      peakDims.append(peakDimContrib.peakDim)
  
  for peakDim2 in peakDims:
    experiment = peakDim2.peak.peakList.dataSource.experiment
    if experiment.molSystems:
      if molSystem not in experiment.molSystems:
        experiments[experiment] =  True

    else:
      experiments2[experiment] = True

  if experiments:
    expNames = ' '.join([e.name for e in experiments])
    msg  = 'Molecular system %s of assignment doesn\'t match\n' % molSystem.code
    msg += 'previous assignments for experiments '
    msg += '[%s]\nAssociate molecular system with experiments?' % expNames
    if showOkCancel('Warning', msg):
      for experiment in experiments.keys():
        experiment.addMolSystem(molSystem)
      for experiment in experiments2.keys():
        experiment.addMolSystem(molSystem)
    else:
      return None
  else:
    for experiment in experiments2.keys():
      if molSystem not in experiment.molSystems:
        experiment.addMolSystem(molSystem)
      
  return True
  
  
def assignPeakDim(resonance, peakDim, atomSets=None, contrib=None, tolerance=None):
  """
  Assigns a given peak dimension to a given resonance,
  via a given peakDim contribution if specified
  and assign the resonance to atom sets if specified

  .. describe:: Input
  
  Nmr.Resonance, Nmr.PeakDim, List of Nmr.AtomSets,
             Nmr.PeakDimContrib, Float
  
  .. describe:: Output
  
  Nmr.PeakDimContrib
  """

  contrib = assignResToDim(peakDim,resonance,contrib,tolerance)
    
  if atomSets:
    assignAtomsToRes(atomSets,resonance)
    
  return contrib


def clearPeakDim(peakDim,contrib=None):
  """
  Clear a peak dimension of any or a specified peakDim contribution
  Recalculates shifts of the resonances involved 

  .. describe:: Input
  
  Nmr.PeakDim, Nmr.PeakDimContrib
  
  .. describe:: Output
  
  None
  """
  
  if contrib:
    if contrib.isDeleted:
      contribs = []
    else:
      contribs = [contrib,]
  else:
    contribs = peakDim.peakDimContribs 
  
  for contrib in contribs:
    contrib.delete()
      
  for peakContrib in peakDim.peak.peakContribs:
     if not peakContrib.peakDimContribs:
       peakContrib.delete()


def clearResonancePeakDimContribs(resonance,peaks=None):
  """
  Remove a resonance's peak dim contributions present in a 
  specified list of peaks or all peaks in project (if peak=None)

  .. describe:: Input
  
  Nmr.Resonance, List of Nmr.Peaks
  
  .. describe:: Output
  
  None
  """

  if not peaks:
    peaks = []

  peakDict = {}
  for peak in peaks:
    peakDict[peak] = True
  
  peakDims = {}  
  for contrib in resonance.peakDimContribs:
    peakDim = contrib.peakDim
    
    if (not peakDict) or peakDict.get(peakDim.peak):
      peakDims[peakDim] = True
      peakContribs = contrib.peakContribs
      contrib.delete()
      
      for peakContrib in peakContribs:
        if not peakContrib.peakDimContribs:
          peakContrib.delete()


def findResonanceSet(resonance,atomSets):
  """
  Find the resonance set, if any, that connects a resonance to given atomSets
  Or that is shared between atomSets

  .. describe:: Input
  
  Nmr.Resonance, List of Nmr.AtomSets
  
  .. describe:: Output
  
  Nmr.ResonanceSet
  """

  atomSets = list(atomSets)
  result = None
  
  resonanceSet = resonance.resonanceSet
  if resonanceSet:
    atomSets.sort()
    atomSets2 = list(resonanceSet.atomSets)
    atomSets2.sort()
    
    if atomSets == atomSets2:
      result = resonanceSet

  if result is None and len(atomSets) > 1:
    for tryResonanceSet in atomSets[0].resonanceSets:
      if tryResonanceSet in atomSets[1].resonanceSets:
        result = tryResonanceSet
 
  return result
  
def assignAtomsToRes(atomSets,resonance,resonanceSet=None):
  """
  Assign a resonance to given atom sets via a resonanceSet (optionally specified)
  Checks ensure that the resonance ends up in the correct spin system

  .. describe:: Input
  
  List of Nmr.AtomSets, Nmr.Resonance, Nmr.ResonanceSet
  
  .. describe:: Output
  
  Nmr.ResonanceSet
  """
  
  atomSets = list(atomSets)
  nmrProject = resonance.nmrProject

  # NOTE:ED - could possibly insert a temporary logger here

  chemAtomSetRef = None
  for atomSet in atomSets:
    atom = atomSet.findFirstAtom()
    if not setAssignmentMolSystem(atom.residue, resonance=resonance):
      return
    
    if len(atomSets) > 1:
      chemAtomSet = atom.chemAtom.chemAtomSet
      if chemAtomSet and chemAtomSet.chemAtomSet:
        chemAtomSet = chemAtomSet.chemAtomSet
      
      if not chemAtomSetRef:
        chemAtomSetRef = chemAtomSet
              
      #if (not chemAtomSet) or (not chemAtomSet.isProchiral) or (chemAtomSetRef is not chemAtomSet):
      if (not chemAtomSet) or (chemAtomSetRef is not chemAtomSet):
        info = [(aS, aS.findFirstAtom().residue) for aS in atomSets]
        data = ['%d%s %s' % (r.seqCode, r.ccpCode, aS.name) for aS, r in info]
        names = ','.join(data)
        msg  = 'Resonance can only be assigned to multiple sets of atoms if '
        msg += 'the sets are prochiral. This is not true for the input %s. ' % names
        msg += 'Assignment will be made to %s only' % data[0]
        showWarning('Assignment failed', msg)
        atomSets = [atomSets[0],]
        break
  
  if resonance.isotopeCode:  
    if resonance.isotopeCode != 'unknown':
      element = re.match('\d+([A-Z]\D*)', resonance.isotopeCode).group(1)
      if element != atom.chemAtom.elementSymbol:
        data = (resonance.isotopeCode,atom.chemAtom.elementSymbol)
        msg  = 'A %s resonance cannot be assigned to %s atoms'
        showWarning('Assignment failed', msg % data)
        return
    
  if not resonanceSet:
    resonanceSet = findResonanceSet(resonance,atomSets)
      
  oldResonanceSet = resonance.resonanceSet
  if oldResonanceSet and (oldResonanceSet is not resonanceSet):
    if len(oldResonanceSet.resonances) == 1:
      oldResonanceSet.delete()
    else:
      oldAtomSets   = list(oldResonanceSet.atomSets)
      oldResonances = list(oldResonanceSet.resonances)  
          
      oldResonances.remove(resonance)
      for atomSet in atomSets:
        if atomSet in oldAtomSets:
          oldAtomSets.remove(atomSet)
      
      oldResonanceSet.delete() 
      # Other half of a now split prochiral
      nmrProject.newResonanceSet(atomSets=oldAtomSets,
                                 resonances=oldResonances)
         
  if resonanceSet:
    resonances = list(resonanceSet.resonances)

    if not resonance in resonances:
      resonanceSet.addResonance(resonance)
      resonances.append(resonance)

      # NOTE:ED atomSets always appears to be of length 2
      # resonances can be 2 or more

      if len(resonances) > len(atomSets):
        residue = atomSets[0].findFirstAtom().residue
        aName = '/'.join([ass.name for ass in atomSets])
        data  = (len(resonances), residue.seqCode, residue.ccpCode, aName)
        msg   = 'There are more resonances (%d) than atoms sets for %d%s %s'
        showWarning('Redundant resonance', msg % data)

    for atomSet in resonanceSet.atomSets:
      if atomSet not in atomSets:
        resonanceSet.delete()
        resonanceSet = nmrProject.newResonanceSet(atomSets=atomSets,resonances=resonances)
        break
      
  else:
    resonanceSet = nmrProject.newResonanceSet(atomSets=atomSets,resonances=[resonance, ])

  initResonance(resonance)
  resonances = getBoundResonances(resonance, recalculate=True)
 
  if len(resonances) == 1:
    from ccpnmr.analysis.core.MoleculeBasic import getBoundAtoms
    
    resonance2 = resonances[0]
    if not resonance2.resonanceSet:

      isotope    = resonance2.isotope
      if isotope:
        element  = isotope.chemElement.symbol
        
        atomDict = {}
        for atomSet in atomSets:
          for atom in getBoundAtoms(atomSet.findFirstAtom()):
            if atom.chemAtom.elementSymbol == element:
              atomDict[atom] = True
 
        atomSetDict = {}
        for atom in atomDict.keys():
          atomSet2 = atom.atomSet
          if atomSet2:
            atomSetDict[atomSet2] = True
 
        atomSets2 = atomSetDict.keys()
        if len(atomSets2) == 1:
          assignAtomsToRes(atomSets2,resonances[0])
 	 
  return resonanceSet

def assignTentativeAtoms(atomSets, resonance, weight=1.0, doWarnings=False):
  """
  Set a tentative/provisional/trial assignment from 
  a resonance to the input atoms. Optionally set weight of the
  tentative assignment.

  .. describe:: Input
  
  List of Nmr.AtomSets, Nmr.Resonance, Float, Boolean
  
  .. describe:: Output
  
  None
  """
  spinSystem = resonance.resonanceGroup or \
               resonance.nmrProject.newResonanceGroup(resonances=[resonance,])
    
  isotope = resonance.isotopeCode
  atomNames = set(resonance.assignNames)
  residues = set()
  
  for atomSet in atomSets:
    atom = atomSet.findFirstAtom()
    residues.add(atom.residue)
    
    element = atom.chemAtom.elementSymbol
    if isotope[-len(element):] == element:
      atomNames.add(atomSet.name)

  resonance.setAssignNames(atomNames)

  assignTentativeSpinSystemResidues(spinSystem, residues, weight, doWarnings)
    
  updateResonanceAnnotation(resonance)

def assignTentativeSpinSystemResidues(spinSystem, residues, weight=1.0, doWarnings=False):
  """
  Set a tentative/provisional/trial residue assignments 
  for a given spin system. Optionally set weight of the
  tentative assignment.

  .. describe:: Input
  
  Nmr.ResonanceGroup, List of MolSystem.Residues, Float, Boolean 
  
  .. describe:: Output
  
  List of Nmr.ResidueProbs
  """

  if spinSystem.residue:
    spinSystem.setResidue(None)
    
    if doWarnings:
      msg = 'Remove sequential spin system links?'

      if getSeqSpinSystemLinks(spinSystem) and showYesNo('Query',msg):
        clearSeqSpinSystemLinks(spinSystem)

  for resonance2 in spinSystem.resonances:
    resonanceSet = resonance2.resonanceSet
    if resonanceSet:
      if len(resonanceSet.resonances) == 1:
        resonanceSet.delete()
      else:
        resonanceSet.removeResonance(resonance2)
  
  residueProbs = []
  for residue in residues:
    residueProb = spinSystem.findFirstResidueProb(possibility=residue) or \
                  spinSystem.newResidueProb(possibility=residue)
    residueProb.weight = weight
    residueProbs.append(residueProb)
    
  return residueProbs

def setResonanceTypeFromRefExp(resonance, refExpDimRef):
  """
  Set a resonance's assign names based on a refExpDimRef
  to which it is assigned. Currently limited to CO, CA, CB, HA.

  .. describe:: Input
  
  Nmr.Resonance, NmrExpPrototype.RefExpDimRef
  
  .. describe:: Output
  
  None
  """

  name = None
  if (not resonance.resonanceSet) and (not resonance.assignNames):
    atomSites = refExpDimRef.expMeasurement.atomSites
    if len(atomSites) == 1:
      if atomSites[0].name in ('CO','CA','CB','HA'):
        name = atomSites[0].name
        
        if name == 'CO':
          name = 'C'
  
  if name:    
    resonance.setAssignNames([name, ])
    if resonance.resonanceGroup:
      residue = resonance.resonanceGroup.residue
      if residue:
        atom = residue.findFirstAtom(name=name)
        if atom and atom.atomSet:
          assignAtomsToRes([atom.atomSet],resonance)
 
    updateResonanceAnnotation(resonance)

def assignResonanceType(resonance, atomSets=None, assignNames=None):
  """
  Sets the assign names for a resonance according to given AtomSets or atom names
  without performing a real assignment. If no atom sets or assign names are given
  then all typing information is removed.

  .. describe:: Input
  
  Nmr.Resonance, List of Nmr.AtomSets or None, List of Words or NOne
  
  .. describe:: Output
  
  None
  """

  if not (atomSets or assignNames):
    if not resonance.resonanceSet:
      resonance.assignNames = []
      updateResonanceAnnotation(resonance)
    return

  isotope = resonance.isotopeCode
  atomNames = []
  
  for atomSet in atomSets or []:
    atom = atomSet.findFirstAtom()
    element = atom.chemAtom.elementSymbol
    
    if isotope[-len(element):] == element:
      name = atomSet.name
  
      if name not in atomNames:
        atomNames.append( name )
  
  if assignNames:
    atomNames.extend(assignNames)
        
  resonance.setAssignNames(atomNames)
  
  if resonance.resonanceSet:
    updateResonanceAnnotation(resonance)
    return
  
  if resonance.resonanceGroup:
    residue = resonance.resonanceGroup.residue
    if residue:
      atomSets = {}
      for atom in residue.atoms:
        atomSet = atom.atomSet
        
        if atomSet and (atomSet.name in atomNames):
          atomSets[atomSet] = True
          
      if atomSets:
        assignAtomsToRes(atomSets.keys(),resonance)
  
  updateResonanceAnnotation(resonance)

def deassignResonance(resonance, clearAssignNames=True):
  """
  Disconnects a resonance from any assigned atomSets.
  Option to remove previous assignment atom names

  .. describe:: Input
  
  Nmr.Resonance
  
  .. describe:: Output
  
  None
  """
  
  if clearAssignNames:
    resonance.setAssignNames([])
    
  resonanceSet = resonance.resonanceSet
  if resonanceSet:
    if len(resonanceSet.resonances) == 1:
      resonanceSet.delete()
    else:
      resonanceSet.removeResonance(resonance)
      updateResonanceAnnotation(resonance)    
      # the resonance is no longer in the
      # resonanceSet's list to be updated
  else:
    updateResonanceAnnotation(resonance)

def deassignSpinSystem(spinSystem, clearResidueType=True):

    assignSpinSystemResidue(spinSystem, None)
    if clearResidueType:
      for residueProb in spinSystem.residueProbs:
        residueProb.delete()
      assignSpinSystemType(spinSystem, None)

def aliasedPeakDimPosition(peakDim):
  """
  Give the position in PPM for a peak dimension

  .. describe:: Input
  
  Nmr.PeakDim
  
  .. describe:: Output
  
  Float (ppm)
  """

  dataDimRef = peakDim.dataDimRef
  position   = peakDim.position
  
  if dataDimRef:
    position += peakDim.numAliasing * dataDimRef.dataDim.numPointsOrig
  
  return position 


def assignResToDim(peakDim, resonance=None, contrib=None,
                   tolerance=None, doWarning=True, peakContribs=None):
  """
  Assign a resonance to a peak dimension, via a specified
  peakDimContrib if needed. Assigns a new resonance if none
  is input. Can also specify a set of peakContribs (grouped
  assigment possibilities) that the assignment relates to.

  .. describe:: Input
  
  Nmr.PeakDim, Nmr.Resonance, Nmr.PeakDimContrib, Float
             Boolean, List of Nmr.PeakContribs
  
  .. describe:: Output
  
  Nmr.PeakDimContrib
  """
  from ccpnmr.analysis.core.PeakBasic import setPeakDimNumAliasing, findNumAliasing
  
  peak = peakDim.peak
  dataDimRef = peakDim.dataDimRef
  if not dataDimRef: # Is sampled
    return
  
  isotopeCode = dataDimRef.expDimRef.isotopeCodes[0]
  if resonance:
    resonanceSet = resonance.resonanceSet
    if resonanceSet:
      residue = resonanceSet.findFirstAtomSet().findFirstAtom().residue
      if not setAssignmentMolSystem(residue, peakDim=peakDim):
        return
  
  peakContribs = set(peakContribs or [])
  if contrib and contrib.peakContribs:
    peakContribs.update(contrib.peakContribs)
  elif (not peakContribs) and (peak.peakContribs):
    # Default => add to all
    peakContribs = set(peak.peakContribs)
    
  if resonance is None:
    resonance = peakDim.topObject.newResonance(isotopeCode=isotopeCode)
  else:
    if resonance.isotopeCode != isotopeCode:
      if resonance.isotopeCode != 'unknown':
        msg = 'Attempt to link %s resonance to %s dimension'
        #raise Exception(msg % (resonance.isotopeCode, isotopeCode))
        print msg % (resonance.isotopeCode, isotopeCode)
        return

    #resonance.isotopeCode = isotopeCode
    for checkContrib in peakDim.peakDimContribs:
      if checkContrib.peakDimComponent:
        continue
    
      if checkContrib.resonance is resonance:
        if peakContribs:
          if peakContribs != set(checkContrib.peakContribs):
            checkContrib.setPeakContribs(peakContribs)
        return checkContrib

  positionB  = peakDim.realValue
  dataDim = dataDimRef.dataDim
  spectrum = dataDim.dataSource
  experiment = spectrum.experiment
  shiftList  = experiment.shiftList

  if doWarning:
    msg = isOnebondMismatched(peakDim, resonance, spectrum)
   
    if msg and not showYesNo('Warning',msg + ' Continue with assignment?'):
      return None

  if not tolerance:
    tolerance = getAnalysisDataDim(dataDim).assignTolerance
  
  if not shiftList:
    nmrProject = peakDim.topObject
    shiftList  = nmrProject.findFirstMeasurementList(className='ShiftList')
    
    if not shiftList:
      shiftList = nmrProject.newShiftList(name='ShiftList 1',
                                          details='Assignment Default',
                                          unit='ppm')
  
    experiment.shiftList = shiftList
  
  unit  = shiftList.unit
  shift = resonance.findFirstShift(parentList=shiftList)

  if shift:    
    positionA = shift.value
    if (positionB is None) or (abs(positionA - positionB) < tolerance):
      if contrib:
        contrib.delete()
        
      contrib = newPeakDimContrib(peakDim, resonance, peakContribs=peakContribs)
      
    else:
      n = findNumAliasing(dataDimRef, unit_converter[(unit,'point')](positionA,dataDimRef))
      if n != peakDim.numAliasing:
        sw = dataDim.numPointsOrig
        position  = peakDim.position + sw*n
        positionC = unit_converter[('point',unit)](position,dataDimRef)

        if abs(positionA - positionC) < tolerance:
          n = findNumAliasing(dataDimRef, position)

          if n != peakDim.numAliasing:
            if doWarning:
              if not showOkCancel('Confirm','Peak dimension will be unaliased'):
                return None
               
            if contrib:
              contrib.delete()
            contrib = newPeakDimContrib(peakDim, resonance,peakContribs=peakContribs)
            setPeakDimNumAliasing(peakDim, n, doWarning=doWarning)
 
          else:
            return None

        else:
          return None

      else:
        return  None
      
  else:
    if contrib:
      contrib.delete()
    contrib = newPeakDimContrib(peakDim, resonance,peakContribs=peakContribs)
 
  return contrib

def addPeakResonances(peaks):
  """
  Add new resonance assignments to all unassigned dimensions
  of the input peaks.

  .. describe:: Input
  
  Nmr.Peaks
  
  .. describe:: Output
  
  Nmr.Resonances
  """
  
  contribs = []
  for peak in peaks:
    for peakDim in peak.peakDims:
      if len(peakDim.peakDimContribs) < 1:
        contrib = assignResToDim(peakDim)
        contribs.append(contrib)
        
  resonances = [c.resonance for c in contribs]
  
  return resonances   

def isOnebondMismatched(peakDim, resonance, spectrum=None):
  """
  Check to make sure a resonance to peakDim assignment obeys
  experimental onbond relationshps. 

  .. describe:: Input
  
  Nmr.PeakDim, Nmr.Resonance, Nmr.DataSource
  
  .. describe:: Output
  
  String (warning message) or None
  """
  
  if peakDim.peakDimContribs:
    return
  
  if resonance.isotopeCode != '1H':
    return
    
  bound = resonance.covalentlyBound
  
  if bound:
    from ccpnmr.analysis.core.ExperimentBasic import getOnebondDataDims
    dim = peakDim.dim
    resonanceSet = resonance.resonanceSet
    
    for dataDim0, dataDim1 in getOnebondDataDims(spectrum):
 
      if dataDim0.dim == dim:
        peakDimB = peakDim.peak.findFirstPeakDim(dim=dataDim1.dim)

      elif dataDim1.dim == dim:
        peakDimB = peakDim.peak.findFirstPeakDim(dim=dataDim0.dim)

      else:
        continue
 
      if len(peakDimB.peakDimContribs) == 1:
        resonanceB = peakDimB.findFirstPeakDimContrib().resonance
        
        if resonanceB in bound:
          continue
 
        for resonanceC in resonanceB.covalentlyBound:
          if resonanceC is resonance:
            continue
          
          elif resonanceSet and (resonanceC.resonanceSet is resonanceSet):
            continue
            
          elif resonanceC.isotopeCode == '1H':
            rText1 = makeResonanceGuiName(resonance)
            rText2 = makeResonanceGuiName(resonanceB)
            rText3 = makeResonanceGuiName(resonanceC)
            msg = 'Dimensions %d and %d are directly bound but ' % (dataDim0.dim, dataDim1.dim)
            msg += 'resonance %s is already bound to %s not %s' % (rText2, rText3, rText1)
            return msg
 

def newPeakDimContrib(peakDim, resonance, peakContribs=None):
  """
  Make a peakDim to resonance assignment via a PeakDimContrib.
  Optional argument to specify which PeakContrib assignment
  possibilities the peakDim assignment relates to.

  .. describe:: Input
  
  Nmr.PeakDim, Nmr.Resonance, List of Nmr.PeakContribs
  
  .. describe:: Output
  
  Nmr.PeakDimContrib
  """

  if peakContribs is None:
    peakContribs = peakDim.peak.peakContribs
  
  existing = peakDim.findFirstPeakDimContrib(resonance=resonance)
  if existing:
    existing.setPeakContribs(peakContribs)
    return existing
  
  if len(peakContribs) < 1:
    return peakDim.newPeakDimContrib(resonance=resonance)

  else:
    return peakDim.newPeakDimContrib(resonance=resonance,
                                     peakContribs=peakContribs)
    

def updateResonShift(resonance, peakDim):
  """
  Update the shift values (and error) for a given resonance
  given assignment to a given peak dimension. The shift list is approriate
  to the peakDim's experiment. Note that shift to peak and shift to peak
  dim links are also set.

  .. describe:: Input
  
  Nmr.Resonance, Nmr.PeakDim
  
  .. describe:: Output
  
  None
  """
  
  peak       = peakDim.peak
  peakList   = peak.peakList
  shiftList  = peakList.dataSource.experiment.shiftList
  
  if not shiftList:
    return
  
  dataDimRef = peakDim.dataDimRef
  if not dataDimRef:
    return
  
  if peakDim.isDeleted:
    shift = resonance.findFirstShift(parentList=shiftList)
    if shift:
      averageShiftValue(shift)
    return  
      
  value = peakDim.realValue # OK, even if split
  if value is None:
    return
     
  shift = resonance.findFirstShift(parentList=shiftList)
  
  if shift is None:
    unit  = dataDimRef.expDimRef.unit or 'ppm' # 1D shapes fix default
    unit2 = shiftList.unit
    
    if unit2 != unit:
      shiftValue = unit_converter[(unit,unit2)](value,dataDimRef)
    else:
      shiftValue = value
      
    shift = shiftList.newShift(value=shiftValue, resonance=resonance)
  
  averageShiftValue(shift)
    
def setQuickShiftList(shift):
  """
  Sets up a list of lists of shifts for a shift list
  accessed by shift.value: quicker than searching all

  .. describe:: Input
  
  Nmr.Shift
  
  .. describe:: Output
  
  None
  """

  shiftList = shift.parentList
  if not hasattr(shiftList,'quickShifts'):
    
    shiftList.quickShifts = {}
    for shift2 in shiftList.measurements:
      isotope = shift2.resonance.isotopeCode
      
      quickShiftDict = shiftList.quickShifts.get(isotope, {})
      if not quickShiftDict:
        shiftList.quickShifts[isotope] = quickShiftDict
      
      # WARNING: this key function has to be same as below and in findMatchingShifts
      key = int(round(10*shift2.value))
      if quickShiftDict.get(key) is None:
        quickShiftDict[key] = []
      quickShiftDict[key].append(shift2)
      shift2.quickShiftIndex = key

  isotope = shift.resonance.isotopeCode
  quickShiftDict = shiftList.quickShifts.get(isotope, {})
  if not quickShiftDict:
    shiftList.quickShifts[isotope] = quickShiftDict
    
  oldKey = hasattr(shift,'quickShiftIndex') and shift.quickShiftIndex

  if shift.isDeleted:
    if oldKey is not False:
      try:
        quickShiftDict[oldKey].remove(shift)
      except:
        pass
      del shift.quickShiftIndex
    return
    
  # WARNING: this key function has to be same as above and in findMatchingShifts
  key = int(round(10*shift.value))
  if oldKey is not False and key == oldKey:
    return

  if oldKey is not False:
    # remove old key
    try:
      quickShiftDict[oldKey].remove(shift)
    except ValueError, e:
      print 'Warning: Quick shift index value failure', e
    except KeyError, e:
      print 'Warning: Quick shift index key failure', e

  if quickShiftDict.get(key) is None:
    quickShiftDict[key] = []

  quickShiftDict[key].append(shift)
  shift.quickShiftIndex = key

def updatePeakDimShifts(peakDim):
  """
  Update the shifts of resonances assigned to
  a peakDim if that peakDim changes

  .. describe:: Input
  
  Nmr.PeakDim
  
  .. describe:: Output
  
  None
  """

  for contrib in peakDim.peakDimContribs:
    if contrib.peakDimComponent:
      continue
    updateResonShift(contrib.resonance,peakDim)
  
def updatePeakShifts(peak):
  """
  Update the shifts of resonances assigned to a peak.
  Required if the peak's figure of merit changes

  .. describe:: Input
  
  Nmr.Peak
  
  .. describe:: Output
  
  None
  """

  for peakDim in peak.peakDims:
    for contrib in peakDim.peakDimContribs:
      if contrib.peakDimComponent:
        continue
      updateResonShift(contrib.resonance,peakDim)       
  
def updateDataDimShifts(analysisDataDim):
  """
  Update the shifts of resonances assigned to a given dim of peaks
  of a given spectrum. Required if the spectrum data dims chem shift
  weighting changes.

  .. describe:: Input
  
  Analysis.AnalysisDataDim  
  
  .. describe:: Output
  
  None
  """

  analysisSpectrum = analysisDataDim.analysisSpectrum
  spectrum = analysisSpectrum.dataSource
  shiftList = spectrum.experiment.shiftList
  
  if not shiftList:
    return
  
  dataDim = analysisDataDim.dataDim
  if dataDim.className != 'FreqDataDim':
    return
  
  dim = dataDim.dim
  resonances = set()
  
  for peakList in spectrum.peakLists:
    for peak in peakList.peaks:
      peakDim = peak.findFirstPeakDim(dim=dim)
      
      for contrib in peakDim.peakDimContribs:
        if contrib.peakDimComponent:
          continue
          
        resonances.add(contrib.resonance)
  
  for resonance in resonances:
    shift = resonance.findFirstShift(parentList=shiftList)

    if shift:
      averageShiftValue(shift)

def setupAssignmentNotifiers():
  """
  Set up the notification calls to automatically update annotations
  and quick shift lookups when assignments change.

  .. describe:: Input
  
  None
  
  .. describe:: Output
  
  None
  """

  assignmentNotifiers(Implementation.registerNotify)
  
def removeAssignmentNotifiers():
  """
  Remove the notification calls to automatically update annotations
  and quick shift lookups when assignments change.

  .. describe:: Input
  
  None
  
  .. describe:: Output
  
  None
  """

  assignmentNotifiers(Implementation.unregisterNotify)


def assignmentNotifiers(notify):
  
  notify(updatePeakDimShifts, 'ccp.nmr.Nmr.PeakDim', 'setPosition')
  notify(updatePeakDimShifts, 'ccp.nmr.Nmr.PeakDim', 'setNumAliasing')
  
  notify(updateContribPeakDimAnnotation, 'ccp.nmr.Nmr.PeakDimContrib', '__init__')
  notify(updateContribPeakDimAnnotation, 'ccp.nmr.Nmr.PeakDimContrib', 'delete')
  
  for func in ('addResonance','removeResonance','setResonances','__init__','delete'):
    notify(updateResonanceSetAnnotation, 'ccp.nmr.Nmr.ResonanceSet', func)
  
  notify(updateResonanceAnnotation, 'ccp.nmr.Nmr.Resonance', 'setName')

  for func in ('setName', 'setResidue', 'setCcpCode'):
    notify(updateResonanceGroupAnnotation, 'ccp.nmr.Nmr.ResonanceGroup')
  
  notify(updateResidueAnnotation, 'ccp.molecule.MolSystem.Residue', 'setSeqCode')
  
  for func in ('__init__','setValue','delete'):
    notify(setQuickShiftList, 'ccp.nmr.Nmr.Shift',func)
  
  notify(updatePeakShifts, 'ccp.nmr.Nmr.Peak', 'setFigOfMerit')

  for func in ('delete','__init__','setWeight','setPossibility'):
    notify(updateResidueProbAnnotation, 'ccp.nmr.Nmr.ResidueProb', func)

  notify(updateDataDimShifts, 'ccpnmr.Analysis.AnalysisDataDim', 'setChemShiftWeight')
    
    
def updateAllShifts(shiftList):
  """
  Recalculate the values for all shifts in a shift list given peak assignments

  .. describe:: Input
  
  Nmr.ShiftList
  
  .. describe:: Output
  
  None
  """
  
  for shift in shiftList.measurements:
    averageShiftValue(shift)
    
    
def averageShiftValue(shift, simulatedPeakScale=0.0001):
  """
  Calculates the value and error for a given shift based upon the
  peaks to which its resonance is assigned. Also sets links to
  appropriate peaks and peak dims. Optional float to scale simulated
  peak contributions (normally so that they have much less influence)

  .. describe:: Input
  
  Nmr.Shift, Float
  
  .. describe:: Output
  
  Float (shift.value)
  """
  from ccpnmr.analysis.core.CouplingBasic import getPeakDimComponentSplitting

  hasApp = hasattr(shift.root, 'application')
  shiftList = shift.parentList
  experiments = shiftList.experiments
  resonance = shift.getResonance()
  peakDimContribs = resonance.peakDimContribs
  sum1  = 0.0
  sum2  = 0.0
  N     = 0.0
  mean  = 0.0
  mean2 = 0.0
  sigma2= 0.0
  sigma = 0.0
  peakDims = []
  peaks = set([])

  for contrib in peakDimContribs:
    peakDim = contrib.peakDim
    peak = peakDim.peak
    
    if peak.figOfMerit == 0.0:
      continue
      
    peakList = peak.peakList
    experiment = peakList.dataSource.experiment
    if experiment not in experiments:
      continue
    
    component = contrib.peakDimComponent
    if component:
      if isinstance(contrib, Nmr.PeakDimContribN):
        continue
      
      expDimRef = component.dataDimRef.expDimRef
      if expDimRef.unit == shiftList.unit:
        # Works for MQ etc
        value = getPeakDimComponentSplitting(component)
      else:
        continue
      
    else:
      value = peakDim.realValue
 
    if hasApp:
      weight = getAnalysisDataDim(peakDim.dataDim).chemShiftWeight
    else:
      weight = 1.0
 
    if peakList.isSimulated:
      weight *= simulatedPeakScale
 
    peakDims.append(peakDim)
    peaks.add(peak)
 
    vw    = value * weight
    sum1 += vw
    sum2 += value * vw
    N    += weight

  if N > 0.0:
    mean  = sum1/N
    mean2 = sum2/N
    sigma2= abs(mean2 - (mean * mean))
    sigma = sqrt(sigma2)
    
  else:
    # resonance has run out of contribs
    # - leave it be, even if orphaned
    # or all dataDimWeighting are zero
    shift.setError( 0.0 )
    return
   
  if shift is not None:
    shift.setValue( mean )
    shift.setError( sigma )
    shift.setPeaks( peaks )
    shift.setPeakDims( peakDims )

  return mean

def splitResonance(resonance):
  """
  Split assignments that were once thought to be from one resonance
  into assignments from two.

  .. describe:: Input
  
  Nmr.Resonance
  
  .. describe:: Output
  
  Nmr.Resonance
  """
  
  nmrProject = resonance.nmrProject
  resonanceB = nmrProject.newResonance(isotopeCode=resonance.isotopeCode)
  
  for contrib in resonance.peakDimContribs:
    peakDim = contrib.peakDim
    peakDim.newPeakDimContrib(resonance=resonanceB)

  return resonanceB
   
def mergeResonances(resonanceB, resonanceA):
  """
  Merge two resonances and their shifts into one

  .. describe:: Input
  
  Nmr.Resonance, Nmr.Resonance
  
  .. describe:: Output
  
  Nmr.Resonance
  """

  from ccpnmr.analysis.core.MoleculeBasic import getResidueMapping
  
  if resonanceB is resonanceA:
    return resonanceA

  if resonanceB.isDeleted:
    return resonanceA

  if resonanceA.isDeleted:
    return resonanceB
  
  removeAssignmentNotifiers()
  
  isotopeA = resonanceA.isotopeCode
  isotopeB = resonanceB.isotopeCode
  
  if isotopeA and isotopeB:
    if isotopeA != isotopeB:
      showWarning('Resonance Merge Failure',
                  'Attempt to merge resonances with different isotope codes')
      setupAssignmentNotifiers()
      return 
  
  mappings = []
  resonanceSet = resonanceB.resonanceSet
  if resonanceSet:
    atomSets = resonanceSet.atomSets
    residue  = resonanceSet.findFirstAtomSet().findFirstAtom().residue
    serials  = [atomSet.serial for atomSet in atomSets]
    serials.sort()
    residueMapping = getResidueMapping(residue)
    for atomSetMapping in residueMapping.atomSetMappings:
      serials2 = list(atomSetMapping.atomSetSerials)
      serials2.sort()
      if serials2 == serials:
        mappings.append([atomSetMapping, atomSets])
  
  # attributes where we have object.resonance
  controlData = {'findFirstMeasurement':('shiftDifferences', 'hExchRates',
                                         'hExchProtections', 'shiftAnisotropies',
                                         't1s', 't1Rhos', 't2s'),
                 'findFirstDerivedData':('pkas',),
                 'findFirstPeakDimContrib':('peakDimContribs',)
                }
  for funcName in controlData:
    for attrName in controlData[funcName]:
      for objectA in list(resonanceA.__dict__.get(attrName)):
        objectB = getattr(objectA.parent, funcName)(resonance=resonanceB)
        if objectB is not None:
          objectA = mergeObjects(objectB, objectA)
  
  # attributes where we have object.resonances
  controlData = {'findFirstMeasurement':('jCouplings',
                                         'noes', 'rdcs', 'dipolarRelaxations'),
                 'findFirstDerivedData':('isotropicS2s', 'spectralDensities',
                                         'datums'),
                 'findFirstPeakDimContribN':('peakDimContribNs',)
                }
  for funcName in controlData:
    for attrName in controlData[funcName]:
      for objectA in list(resonanceA.__dict__.get(attrName)):
        testKey = set(objectA.__dict__['resonances'])
        testKey.remove(resonanceA)
        testKey.add(resonanceB)
        testKey = frozenset(testKey)
        objectB = getattr(objectA.parent, funcName)(resonances=testKey)
    
        if objectB is not None:
          objectA = mergeObjects(objectB, objectA)
  
  resonanceA.setCovalentlyBound([])
  resonanceB.setCovalentlyBound([])
        
  # merge shifts in the same shiftlist
  # NB must be done after other measurements 
  for shiftA in resonanceA.shifts:
    for shiftB in resonanceB.shifts:
      if shiftA.parentList is shiftB.parentList:
        shiftA = mergeObjects(shiftB,shiftA)

  # Get rid of duplicate appData
  for appData in resonanceA.applicationData:
    matchAppData = resonanceB.findFirstApplicationData(application=appData.application,
                                                       keyword=appData.keyword)
    if matchAppData:
      resonanceB.removeApplicationData(matchAppData)
  
  mergeObjects(resonanceB, resonanceA)
  
  # Must be after resonance merge, so that links to peaks are properly set
  for shiftA in resonanceA.shifts:
    averageShiftValue(shiftA)
  
  # Assign names will be merged, but if assigned we only want the correct ones 
  if resonanceA.resonanceSet:
    assignNames = []
    for atomSet in resonanceA.resonanceSet.atomSets:
      assignNames.append( atomSet.name )
      
    resonanceA.setAssignNames(assignNames)  
  
  for atomSetMapping, atomSets in mappings:
    updateAtomSetMapping(atomSetMapping, atomSets)
  
  getBoundResonances(resonanceA, recalculate=True)
  updateResonanceAnnotation(resonanceA)
  
  setupAssignmentNotifiers()
  
  return resonanceA

def updateAtomSetMapping(atomSetMapping, atomSets=None):
  """
  Refresh an AtomSetMapping according to the current assignment status.
  AtomSets optionally passed into increase speed

  .. describe:: Input
  
  Analysis.AtomSetMapping, List of Nmr.AtomSets
  
  .. describe:: Output
  
  None
  """

  if not atomSets:
    atomSets = list(atomSetMapping.atomSets)
  else:
    atomSets = list(atomSets)

  if not atomSets:
    # nothing to be done
    return
  
  resSerials = []
  if atomSetMapping.mappingType == 'ambiguous':
    for atomSet in atomSets:
      for resonanceSet in atomSet.resonanceSets:
        if len(resonanceSet.atomSets) == 1:
          # must be stereo or simple
          resSerials.append(resonanceSet.findFirstResonance().serial)

    if len(resSerials) < 2:
      resSerials = []
  
  elif atomSetMapping.mappingType == 'nonstereo':
    for resonanceSet in atomSets[0].resonanceSets:
      if len(resonanceSet.atomSets) > 1:
        # must be non-stereo
        for resonance in resonanceSet.resonances:
          letter1 = getAmbigProchiralLabel(resonance)
          letter2 = atomSetMapping.name[-1]
          if letter2 == '*':
            if letter1 == atomSetMapping.name[-2]:
              resSerials.append(resonance.serial)
              break
          
          elif letter1 == atomSetMapping.name[-1]:
            resSerials.append(resonance.serial)
            break

  else:
    for resonanceSet in atomSets[0].resonanceSets:
      if len(resonanceSet.atomSets) == 1:
        # must be stereo or simple
        resSerials.append(resonanceSet.findFirstResonance().serial)

  atomSetMapping.setResonanceSerials( resSerials )

def findSpinSystem(resonance):
  """
  Find the spin system in which a resonance resides
  makes a new spin system if none is found

  .. describe:: Input
  
  Nmr.Resonance
  
  .. describe:: Output
  
  Nmr.ResonanceGroup
  """

  resonanceGroup = resonance.resonanceGroup
  resonanceSet   = resonance.resonanceSet
  if resonanceGroup:
    spinSystem = resonanceGroup
    
  elif resonanceSet:
    # find via residue
    residue     = resonanceSet.findFirstAtomSet().findFirstAtom().residue
    molResidue  = residue.molResidue
    molType     = molResidue.molType
    ccpCode     = molResidue.ccpCode
    nmrProject  = resonance.nmrProject
    spinSystem  = None
    spinSystems = nmrProject.findAllResonanceGroups(ccpCode=ccpCode, molType=molType)
    
    for trialSpinSystem in spinSystems:
      if trialSpinSystem.residue is residue:
        spinSystem = trialSpinSystem
        break
        
    if not spinSystem:
      spinSystem = nmrProject.newResonanceGroup()
      
  else:
    spinSystem = resonance.nmrProject.newResonanceGroup()

  return spinSystem

def addSpinSystemResonance(spinSystem, resonance):
  """
  Add a resonance to a given spin system

  .. describe:: Input
  
  Nmr.ResonanceGroup, Nmr.Resonance
  
  .. describe:: Output
  
  List of Nmr.Resonances (belonging to the spin system)
  """

  # when merging spin systems some resonances can end up deleted
  # it's easiest just to put that check here
  if resonance.isDeleted:
    return ()

  residue = getResonanceResidue(resonance)

  if resonance in spinSystem.resonances:
    if residue is spinSystem.residue:
      return spinSystem.resonances
  
  else:
    spinSystem.addResonance(resonance)
  
  if spinSystem.residue and residue:
    if residue != spinSystem.residue:
      assignResonanceResidue(resonance,spinSystem.residue)

  elif residue:
    spinSystem = assignSpinSystemResidue(spinSystem,residue)
    
  elif spinSystem.residue:
    assignResonanceResidue(resonance,spinSystem.residue)
 
  updateResonanceAnnotation(resonance)
    
  return spinSystem.resonances


def removeSpinSystemResonance(spinSystem, resonance):
  """
  Remove a resonance from a given spin system

  .. describe:: Input
  
  Nmr.ResonanceGroup, Nmr.Resonance
  
  .. describe:: Output
  
  None
  """

  if resonance not in spinSystem.resonances:
    return

  spinSystem.removeResonance(resonance)
  updateResonanceAnnotation(resonance)


def mergeSpinSystems(spinSystemB, spinSystemA):
  """
  Merge the resonances from two spin systems into one spin system

  .. describe:: Input
  
  Nmr.ResonanceGroup, Nmr.ResonanceGroup
  
  .. describe:: Output
  
  Nmr.ResonanceGroup
  """

  if spinSystemB is spinSystemA:
    return spinSystemA
  
  if spinSystemB.isDeleted:
    return spinSystemA
  
  if spinSystemA.isDeleted:
    return spinSystemB
  
  residueA = spinSystemA.residue
  if not residueA:
    spinSystemA.setResidue(spinSystemB.residue)
    spinSystemA.setCcpCode(spinSystemB.ccpCode)
    spinSystemA.setMolType(spinSystemB.molType)
    
  if not spinSystemA.ccpCode:
    if residueA:
      spinSystemA.setCcpCode(residueA.ccpCode)
    else:
      spinSystemA.setCcpCode(spinSystemB.ccpCode)
      
  if not spinSystemA.molType:
    if residueA:
      spinSystemA.setMolType(residueA.molResidue.molType)
    else:
      spinSystemA.setMolType(spinSystemB.molType)


  resonanceAssignments = {}
  for resonance in spinSystemA.resonances:
    resonanceSet = resonance.resonanceSet
    
    if resonanceSet:
      atomSets = list(resonanceSet.atomSets)
      atomSets.sort()
      resonanceAssignments[tuple(atomSets)] = resonance
  
    
  mergeList = []
  for resonance in spinSystemB.resonances:
      
    removeSpinSystemResonance(spinSystemB, resonance)
    addSpinSystemResonance(spinSystemA, resonance)
      
    resonanceSet = resonance.resonanceSet
    
    if resonanceSet:
      atomSets = list(resonanceSet.atomSets)
      atomSets.sort()
      resonanceA = resonanceAssignments.get(tuple(atomSets))
      
      if resonanceA:
        n = len(atomSets)
        if n == 1:
          mergeList.append((resonance, resonanceA))
        #elif len(resonanceSet.resonances) == n : # prochiral
  
  for residueProb in list(spinSystemB.residueProbs):
    residue = residueProb.possibility
    
    if not spinSystemA.findFirstResidueProb(possibility=residue):
      spinSystemA.newResidueProb(possibility=residue,
                                 weight=residueProb.weight)
                
  if mergeList:
    residue = spinSystemA.residue
    if residue:
      name = '%d%s' % (residue.seqCode,residue.ccpCode)
    else:
      name = '%s{%d}' % (spinSystemA.ccpCode or '',spinSystemA.serial)
  
    getName = makeResonanceGuiName
    resonanceText = ','.join([getName(r1,fullName=False) for r1,r2 in mergeList])
  
    msg = 'Merge duplicate %s resonances in spin system %s' % (resonanceText, name)
    if showYesNo('Query',msg): 
      for r1, r2 in  mergeList:
        mergeResonances(r1,r2)  

  for link in spinSystemB.findAllResonanceGroupProbs(linkType='sequential',isSelected=True):
    makeSeqSpinSystemLink(spinSystemA, link.possibility, delta=link.sequenceOffset)

  for link in spinSystemB.findAllFromResonanceGroups(linkType='sequential',isSelected=True):
    makeSeqSpinSystemLink(link.fromResonanceGroup, spinSystemA, delta=link.sequenceOffset)
  
  if not spinSystemB.isDeleted:
    # Could be already deleted due to recursive merge
    spinSystemB.delete()

  return spinSystemA

def swapSpinSystemResonances(spinSystemB, spinSystemA):
  """
  Swaps the resonances for as pair of spin systems.
  Spin systen assignments updates according to resonance assignments

  .. describe:: Input
  
  Nmr.ResonanceGroup, Nmr.ResonanceGroup
  
  .. describe:: Output
  
  None
  """

  resonancesA = list(spinSystemA.resonances)
  for resonance in spinSystemA.resonances:
    spinSystemA.removeResonance(resonance)

  resonancesB = list(spinSystemB.resonances)
  for resonance in spinSystemB.resonances:
    spinSystemB.removeResonance(resonance)
 
  molTypeB = spinSystemB.molType
  residueB = spinSystemB.residue
  ccpCodeB = spinSystemB.ccpCode
  molTypeA = spinSystemA.molType
  residueA = spinSystemA.residue
  ccpCodeA = spinSystemA.ccpCode

  spinSystemB.setMolType(molTypeA)
  spinSystemB.setResidue(residueA)
  spinSystemB.setCcpCode(ccpCodeA)
  spinSystemA.setMolType(molTypeB)
  spinSystemA.setResidue(residueB)
  spinSystemA.setCcpCode(ccpCodeB)

  for resonance in resonancesA:
    addSpinSystemResonance(spinSystemB, resonance)   

  for resonance in resonancesB:
    addSpinSystemResonance(spinSystemA, resonance)   

def newSpinSystem(project):
  """
  Make a new spin system for a project

  .. describe:: Input
  
  Project
  
  .. describe:: Output
  
  Nmr.ResonanceGroup (spin system)
  """

  spinSystem = project.currentNmrProject.newResonanceGroup()
  return spinSystem

def getResidueResonances(residue, atomType=None):
  """
  Find the resonances (for a given atom type if specified) assigned to a given residue 

  .. describe:: Input
  
  Nmr.Residue, Word (ChemComp.ChemAtom.elementSymbol)
  
  .. describe:: Output
  
  List of Nmr.Resonances
  """
  
  resonanceDict = {}
  for atom in residue.atoms:
    if atomType and (atom.chemAtom.elementSymbol != atomType):
      continue
   
    if atom.atomSet and atom.atomSet.resonanceSets:
      for resonanceSet in atom.atomSet.resonanceSets:
        for resonance in resonanceSet.resonances:
          resonanceDict[resonance] = True
  
  return resonanceDict.keys()

# NOTE: as of 9 Feb 2009 this function has been copied into ccp/util/Assignment.py
# TBD: remove it from here (part of general tidying up exercise)
def getResonanceResidue(resonance):
  """
  Find the residue, if any, to which a resonance is assigned via atomSets etc

  .. describe:: Input
  
  Nmr.Resonance
  
  .. describe:: Output
  
  Nmr.Residue or None
  """

  residue = None
  if resonance.resonanceSet:
    residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
  
  return residue

def getResonanceMolSystem(resonance):
  """
  Find, if any exists, the molecular system associated with a given resonance
  by virtue of a direct assignment or spin system link.

  .. describe:: Input
  
  Nmr.Resonance
  
  .. describe:: Output
  
  Nmr.MolSystem or None
  """

  molSystem = None
  if resonance.resonanceSet:
    atom = resonance.resonanceSet.findFirstAtomSet().findFirstAtom()
    molSystem = atom.topObject
  
  else:
    spinSystem = resonance.resonanceGroup
    if spinSystem:
      if spinSystem.residue:
        molSystem = spinSystem.residue.chain.molSystem
 
      elif spinSystem.chains:
        molSystems = []
 
        for chain in spinSystem.chains:
          molSystem0 = chain.molSystem
          if molSystem0 not in molSystems:
            molSystems.append( molSystem0 )

        if len(molSystems) == 1:
          molSystem = molSystems[0]

  return molSystem

def assignResonanceResidue(resonance, residue):
  """
  Changes the assignment of a resonance to a given residue if an 
  appropriately named atom set can be found for the resonance's
  assignName.

  .. describe:: Input
  
  Nmr.Resonance, Nmr.Residue
  
  .. describe:: Output
  
  None
  """


  #print 'assignResonanceResidue', makeResonanceGuiName(resonance), residue.seqId
  if residue:
    # we could be setting the residue to None
    if not setAssignmentMolSystem(residue, resonance=resonance):
      return
  
  resonanceSet = resonance.resonanceSet
  if resonanceSet:
    assignNames = set()
    
    for atomSet in resonanceSet.atomSets:
      assignNames.add(atomSet.name)
  
  else:
    assignNames =  set(resonance.assignNames)
    
  assignedResidue = getResonanceResidue(resonance)
  
  if assignNames and (assignedResidue is not residue):
    resonance.setAssignNames(assignNames)
    
    if residue:
      atomSets = []
      isotopeCode = resonance.isotopeCode
      if ('HA' in assignNames) and (residue.ccpCode == 'Gly'):
        assignNames.add('HA2')
        assignNames.add('HA3')
      
      for atom in residue.atoms:
        element = atom.chemAtom.elementSymbol
        
        if isotopeCode[-len(element):] != element:
          continue
      
        atomSet = atom.atomSet
        if atomSet:
          if atomSet.name in assignNames:
            if atomSet not in atomSets:
              atomSets.append(atomSet)
      
      if resonanceSet:
        resonanceSet.delete()
      
      if atomSets:
        # might not always find atom sets
        # e.g. if changing residue types
        assignAtomsToRes(atomSets,resonance,None)
    else:
      if resonanceSet:
        resonanceSet.delete()
        #updateResonanceAnnotation(resonance)

  if not resonance.resonanceGroup:
    nmrProject  = resonance.nmrProject
    spinSystem  = None
    molResidue  = residue.molResidue
    ccpCode     = molResidue.ccpCode
    molType     = molResidue.molType
    spinSystems = nmrProject.findAllResonanceGroups(ccpCode=ccpCode,
                                                    molType=molType)
    
    for extantSpinSystem in spinSystems:
      if extantSpinSystem.residue is residue:
        spinSystem = extantSpinSystem
        break
        
    if not spinSystem:
      spinSystem = nmrProject.newResonanceGroup()
      assignSpinSystemResidue(spinSystem,residue)
      
    if resonance not in spinSystem.resonances:
      spinSystem.addResonance(resonance)
      updateResonanceAnnotation(resonance)
      
  else:
    updateResonanceAnnotation(resonance)


def assignSpinSystemType(spinSystem,ccpCode=None,molType=None):
  """
  Assigns a spin system to a given type of possible residue

  .. describe:: Input
  
  Nmr.ResonanceGroup, Word(Molecule.MolResidue.ccpCode),
  Word(Molecule.MolResidue.molType)
  
  .. describe:: Output
  
  None
  """
  
  # NBNB TBD inconsistent with assignSpinSystemResidue - ResidueTypeProbs not reset.
  
  residue = spinSystem.residue
  if residue and residue.molResidue.ccpCode != ccpCode:
    assignSpinSystemResidue(spinSystem,None)
    spinSystem.setCcpCode(ccpCode)
    
    if molType:
      spinSystem.setMolType(molType)
      
  else:
    spinSystem.setCcpCode(ccpCode)
    
    if molType:
      spinSystem.setMolType(molType)
      
    for resonance in spinSystem.resonances:
      updateResonanceAnnotation(resonance)


def assignSpinSystemResidue(spinSystem,residue=None, warnMerge=False):
  """
  Connect a given spin system to a residue. Resonances in the spin
  system are assigned where possible. Also used to disconnect from
  a residue if Residue is None. Optionally merge spin systems
  assigned to same residue, given warning.

  .. describe:: Input
  
  Nmr.ResonanceGroup, Nmr.Residue, Boolean, Boolean
  
  .. describe:: Output
  
  Nmr.ResonanceGroup
  """
  
  from ccpnmr.analysis.core.MoleculeBasic import getLinkedResidue
  
  for resonanceProb in spinSystem.resonanceProbs:
    # NBNB TBD should this happen??? Rasmus Jan 2010
    resonanceProb.delete()

  if residue:
    #print 'assignSpinSystemResidue', residue.seqId
    nmrProject  = spinSystem.nmrProject
    molResidue  = residue.molResidue
    ccpCode     = molResidue.ccpCode
    molType     = molResidue.molType
    resDict     = {}
    
    for ss in nmrProject.resonanceGroups:
      if ss.residue:
        resDict[ss.residue] = ss
    
    spinSystem.setResidue(None)
    spinSystem.setCcpCode(ccpCode)
    
    duplicateSpinSystem = resDict.get(residue)
    if duplicateSpinSystem:
      resonancesD = duplicateSpinSystem.resonances
      contribsD = set()
      for resonanceD in resonancesD:
        for contribD in resonanceD.peakDimContribs:
          contribsD.add(contribD)
    
      if not contribsD:
        clearSeqSpinSystemLinks(duplicateSpinSystem)
        spinSystem = mergeSpinSystems(duplicateSpinSystem, spinSystem)
    
      elif warnMerge:
        msg = 'There\'s an existing %d %s spin system. OK to merge?'
        
        if showYesNo('Warning', msg % (residue.seqCode,ccpCode)):
          spinSystem = mergeSpinSystems(duplicateSpinSystem, spinSystem)
    
    spinSystem.setResidue(residue)
    spinSystem.setMolType(molType)
    
    for resonance in spinSystem.resonances:
      assignResonanceResidue(resonance,residue)

    residueN = getLinkedResidue(residue, 'prev')
    residueC = getLinkedResidue(residue, 'next')
    
    if residueN:
      connSpinSystemN = findConnectedSpinSystem(spinSystem, delta=-1)
      assnSpinSystemN = resDict.get(residueN)
      if assnSpinSystemN:
        if connSpinSystemN and connSpinSystemN is not assnSpinSystemN:
          clearSeqSpinSystemLinks(connSpinSystemN, delta=1)
          makeSeqSpinSystemLink(assnSpinSystemN, spinSystem, delta=1)
        
      elif connSpinSystemN:
        assignSpinSystemResidue(connSpinSystemN, residueN, warnMerge=0)

    if residueC:
      connSpinSystemC = findConnectedSpinSystem(spinSystem, delta=1)
      assnSpinSystemC = resDict.get(residueC)
      if assnSpinSystemC:
        if connSpinSystemC and connSpinSystemC is not assnSpinSystemC:
          clearSeqSpinSystemLinks(connSpinSystemC, delta=-1)
          makeSeqSpinSystemLink(spinSystem, assnSpinSystemC, delta=1)
        
      elif connSpinSystemC:
        assignSpinSystemResidue(connSpinSystemC, residueC, warnMerge=0)
     
  else:
    spinSystem.setResidue(None)
    if warnMerge:
      msg = 'Remove sequential spin system links?'

      if getSeqSpinSystemLinks(spinSystem) and showYesNo('Query',msg):
        clearSeqSpinSystemLinks(spinSystem)

    for resonance in spinSystem.resonances:
      assignResonanceResidue(resonance,None)
  
  return spinSystem    

def addPeakResonancesToSpinSystem(peaks):
  """
  Add the resonances assigned to given peaks to a spin system.
  Makes a new spin system if none exists for the resonances
  Multiple spin systems will be merged, if present, for resonances upon confirmation

  .. describe:: Input
  
  List of Nmr.Peaks
  
  .. describe:: Output
  
  Nmr.ResonanceGroup
  """
   
  # TBD check experiment type of the peak
  
  if not peaks:
    return
  
  resonances  = []
  for peak in peaks:
    for peakDim in peak.peakDims:
      for contrib in peakDim.peakDimContribs:
        if contrib.peakDimComponent:
          continue
        resonances.append(contrib.resonance)
   
  spinSystems = []
  for resonance in resonances:
    resonanceGroup = resonance.resonanceGroup
    if resonanceGroup and (resonanceGroup not in spinSystems):
      spinSystems.append(resonanceGroup)

  spinSystem = None
  if len(spinSystems) == 1:
    spinSystem = spinSystems[0]
  elif len(spinSystems) > 1:
    msg  = 'There are multiple spin systems for these peaks.\n'
    msg += 'Continue and merge spin systems together?'
    if showOkCancel('Confirm',msg):
      spinSystem = spinSystems[0]
      for spinSystem2 in spinSystems[1:]:
        mergeSpinSystems(spinSystem2,spinSystem)
    else:
      return
  
  if spinSystem is None:
    spinSystem = peaks[0].topObject.newResonanceGroup()

  for resonance in resonances:
    addSpinSystemResonance(spinSystem,resonance)

  return spinSystem

def propagatePeakAssignments(peaks, refPeak=None, cleanNonRef=False,
                             tolerances=None, warnUnalias=False):
  """
  Propogates approprate assignments across a group of peaks
  e.g. to spead F1 F3 assignments to a column of peaks in the same spin system.
  If a refPeak is present it will not be affected but its assignments will
  be used instead of the assignments from the other peaks.
  cleanNonRef Boolean option is used when there is a reference peak so
  that the non-reference peaks may be cleaned prior to propagation.
  Optional argument to pass in the tolerences for each dimension.
  Note peaks can be of different dimensions, but if you pass in 
  tolerances, tolerances map 1:1 to peak.sortedPeakDims()

  .. describe:: Input
  
  List of Nmr.Peaks, Nmr.Peak, Boolean, Dict of Resonance.isotopeCode:Float, Boolean
  
  .. describe:: Output
  
  None
  """

  if refPeak:
    peaksIn = [refPeak, ]
  else:
    peaksIn = peaks
  
  if not tolerances:
    tolerances = []
  
  dimResonances = {}
  resonanceDims = {}
  for peak in peaksIn:
    for i, peakDim in enumerate(peak.sortedPeakDims()):
      dataDim = peakDim.dataDim
      expDimRef = dataDim.expDim.findFirstExpDimRef()
      
      if not expDimRef:
        continue
      
      key = expDimRef.isotopeCodes
      if dimResonances.get(key) is None:
        dimResonances[key] = []
        
      if peakDim.peakDimContribs:
        # could be in different spectra
 
        for contrib in peakDim.peakDimContribs:
          resonance = contrib.resonance
          
          dimResonances[key].append(resonance)
          if resonanceDims.get(resonance) is None:
            resonanceDims[resonance] = []
          
          if i not in resonanceDims[resonance]:
            resonanceDims[resonance].append(i)

  if refPeak and cleanNonRef:
    for peak in peaks:
      if peak is refPeak:
        continue
 
      for peakDim in peak.peakDims:
        clearPeakDim(peakDim)

  shiftRanges = {}
  for peak in peaks:
    if peak is refPeak:
      continue

    for i, peakDim in enumerate(peak.sortedPeakDims()):
      dataDimRef = peakDim.dataDimRef
    
      if dataDimRef:
        dataDim = dataDimRef.dataDim
              
        if dataDim not in shiftRanges:
          shiftMin, shiftMax = getDataDimFullShiftRange(dataDim)
          shiftRanges[dataDim] = (shiftMin, shiftMax)
        else:
          shiftMin, shiftMax = shiftRanges[dataDim]
        
        if i < len(tolerances):
          tolerance = tolerances[i]
        else:
          tolerance = getAnalysisDataDim(dataDim).assignTolerance
        
        key = dataDimRef.expDimRef.isotopeCodes
        pValue = peakDim.realValue

        extantResonances = []
        for contrib in peakDim.peakDimContribs:
          if contrib.peakDimComponent:
            continue
          extantResonances.append(contrib.resonance)
 
        assignResonances = []
        closeResonances  = []
        for resonance in dimResonances[key]:
          if resonance not in extantResonances:
            shiftList = peak.peakList.dataSource.experiment.shiftList
            shift = resonance.findFirstShift(parentList=shiftList)
 
            if shift:
              # Could result in unaliasing the peak

              sValue = shift.value
              # Only assign if within known bounds
              if not (shiftMin < sValue < shiftMax): # Inside, not on edge
                continue
              
              assignResonances.append(resonance)
              
              if abs(sValue-pValue) <= tolerance:
                closeResonances.append(resonance)
 
            elif i in resonanceDims.get(resonance, []):
              # No shift so only propagate across the same dim numbers
              assignResonances.append(resonance)
 
        # Can't have both aliased and unaliased resonances: go for the
        # unaliased/close ppm ones in preference  
        
        if closeResonances:
          for resonance in closeResonances:
            assignResToDim(peakDim, resonance, tolerance=tolerance,
                           doWarning=False)
          
        elif not extantResonances:
          # Don't risk aliasing changes if already assigned
          # warn for aliasing changes
          for resonance in assignResonances:
            assignResToDim(peakDim, resonance, tolerance=tolerance,
                           doWarning=warnUnalias)
        
    
def getPeakDimPpm(peakDim):
  """
  Gives the position in PPM for a given peak dimension

  .. describe:: Input
  
  Nmr.PeakDim
  
  .. describe:: Output
  
  Float (ppm)
  """

  dataDimRef = peakDim.dataDimRef
  if dataDimRef:
    return pnt2ppm( aliasedPeakDimPosition(peakDim), peakDim.dataDimRef )

  else:
    return peakDim.position

def clearSeqSpinSystemLinks(spinSystem, delta=None):
  """
  Remove sequential links for a spin system. A sequence offset may be specified
  otherwise all sequential links are removed.

  .. describe:: Input
  
  Nmr.ResonanceGroup, Int
  
  .. describe:: Output
  
  None
  """

  links = getSeqSpinSystemLinks(spinSystem, delta=delta)
  for link in links:
    link.delete()

def makeSeqSpinSystemLink(spinSystemA, spinSystemB, delta=1):
  """
  Make a sequential link with a give offset between two spin systems.
  Assigning A to B with a delta of 1 assumes B is C-terminal/upstream of A

  .. describe:: Input
  
  Nmr.ResonanceGroup, Nmr.ResonanceGroup, Int
  
  .. describe:: Output
  
  Nmr.ResonanceGroupProb
  """

  if spinSystemA is spinSystemB:
    showWarning('Failure','Attempt to link spin system to itself.')
    return

  residueA = spinSystemA.residue
  residueB = spinSystemB.residue

  if residueA:
    idA = '%d %s' % (residueA.seqCode, residueA.ccpCode)
    residueC = residueA.chain.findFirstResidue(seqId=residueA.seqId + delta)
    if not residueC:
      if delta >= 0:
        d = '+%d' % delta
      else:
        d = '%d' % delta  
      showWarning('Failure','Impossible spin system link attempted: %s to i%s' % (idA,d))
      return

    #idC = '%d %s' % (residueC.seqCode, residueC.ccpCode)
    
    if residueC is not residueB:
      assignSpinSystemResidue(spinSystemB, residueC, warnMerge=False) 
  
  elif residueB:
    idB = '%d %s' % (residueB.seqCode, residueB.ccpCode)
    residueC = residueB.chain.findFirstResidue(seqId = residueB.seqId - delta)
    if not residueC:
      print 'Impossible spin system link attempted: %s to i - %d' % (idB,delta)
    else:
      assignSpinSystemResidue(spinSystemA, residueC, warnMerge=False) 
    
  
  clearSeqSpinSystemLinks(spinSystemA, delta=delta)
  clearSeqSpinSystemLinks(spinSystemB, delta=-delta)
  
  link = spinSystemA.findFirstResonanceGroupProb(linkType='sequential',
                                                 possibility=spinSystemB)
  
  if link:
    link.sequenceOffset = delta
    link.isSelected = True
    
  else:
    link = spinSystemA.newResonanceGroupProb(linkType='sequential', isSelected=True,
                                             possibility=spinSystemB, sequenceOffset=delta)
  
  return link

def getSeqSpinSystemLinks(spinSystem, delta=None):
  """
  Get any sequential spin system links (resonanceGroupProbs).
  An optional sequence offset may be specified.

  .. describe:: Input
  
  Nmr.ResonanceGroup, Int
  
  .. describe:: Output
  
  List of Nmr.ResonanceGroupProbs
  """

  seqLinks = {}
  for link in spinSystem.findAllResonanceGroupProbs(linkType='sequential',isSelected=True):
    if delta is None:
      seqLinks[link] = None
    
    elif link.sequenceOffset == delta:
      seqLinks[link] = None

  for link in spinSystem.findAllFromResonanceGroups(linkType='sequential',isSelected=True):
    if delta is None:
      seqLinks[link] = None
    
    elif link.sequenceOffset == -delta:
      seqLinks[link] = None

  return seqLinks.keys()

def findConnectedSpinSystem(spinSystem, delta=-1):
  """
  Find a spin system sequentially connected to the input one with given sequence offset.

  .. describe:: Input
  
  Nmr.ResonanceGroup, Int
  
  .. describe:: Output
  
  Nmr.ResonanceGroup
  """

  spinSystemB = None
  
  nLinks = getSeqSpinSystemLinks(spinSystem, delta=delta)
  if nLinks:
    nLink = nLinks[0]
    if spinSystem is nLink.parent:
      spinSystemB = nLink.possibility
      
    else:
      spinSystemB = nLink.parent

  return spinSystemB

def findConnectedSpinSystems(spinSystem, delta=None):
  """
  Find spin systems sequentially connected to the input one with given sequence offset.

  .. describe:: Input
  
  Nmr.ResonanceGroup, Int  
  .. describe:: Output
  
  Nmr.ResonanceGroup
  """

  result = []
  
  for link in getSeqSpinSystemLinks(spinSystem, delta=delta):
    if spinSystem is link.parent:
      result.append(link.possibility)
    else:
      result.append(link.parent)
  #
  return result

def addPeakResonancesToSeqSpinSystems(peak, seqOffsets):
  """
  Set the spin systems of a peak's resonances to have sequential
  connectivity according to sequence offsets. Sequence offsets for
  HN(CO)CA would be (None,-1,None), for H,C,N dims

  .. describe:: Input
  
  Nmr.Peak, List of Int (Nmr.ResonanceGroupProb.sequenceOffset)
  
  .. describe:: Output
  
  List of Nmr.ResonanceGroups
  """
    
  assert len(peak.peakDims) == len(seqOffsets)
  assert None in seqOffsets # otherwise no reference point

  spinSystems = []
  resonanceList = []
  for i, peakDim in enumerate(peak.sortedPeakDims()):
    spinSystem = None
    resonances = []
    for contrib in peakDim.peakDimContribs:
      resonance = contrib.resonance
      resonances.append(resonance)
    
      if resonance.resonanceGroup:
        if not spinSystem:
          spinSystem = resonance.resonanceGroup

        elif spinSystem is not resonance.resonanceGroup:
          msg  = 'There are multiple spin systems for peak dimension %d.\n' % (i+1)
          msg += 'Continue and merge spin systems together?'
          if showOkCancel('Confirm', msg):
            mergeSpinSystems(resonance.resonanceGroup,spinSystem)
          else:
            return

    resonanceList.append(resonances)
    spinSystems.append( spinSystem )

  ref = None
  I = 0
  for i, spinSystem in enumerate(spinSystems):
    if spinSystem is not None:
      if seqOffsets[i] is None:
        if ref is None:
          ref = spinSystem
          I = i
          
        else:
          if spinSystem is not ref:
            msg  = 'Dimensions %d and %d have different spin systems.\n' % (I+1,i+1)
            msg += 'Continue and merge spin systems together?'
            if showOkCancel('Confirm', msg):
              mergeSpinSystems(spinSystem, ref)
            else:
              return
           
  if ref is not None:
    for i, seqOffset in enumerate(seqOffsets):
    
      if seqOffset:
        spinSystem = findConnectedSpinSystem(ref, seqOffset)
        if spinSystems[i] is ref:
          if seqOffsets[i] < 0:
            deltaText = '%d' % seqOffset
          else:
            deltaText = '+%d' % seqOffset
          showWarning('Failure','Spin system cannot be both i and i%s (dimension %d)' % (deltaText,i+1))
          continue
          
          
        if spinSystem and spinSystems[i]:
          if spinSystem is not spinSystems[i]:
            if (not spinSystem.residue) or (not spinSystems[i].residue):
              if seqOffsets[i] < 0:
                deltaText = '%d' % seqOffset
              else:
                deltaText = '+%d' % seqOffset
              
              msg =  'There is an i%s spin system already present (dimension %d).\n' % (deltaText, i+1)
              msg += 'Merge spin systems together?'
              if showOkCancel('Confirm', msg):
                spinSystem = mergeSpinSystems(spinSystems[i],spinSystem)
              else:
                spinSystem = None

            elif spinSystem.residue is spinSystems[i].residue:
              name = '%d%s' % (spinSystem.residue.seqCode,spinSystem.residue.ccpCode)
              msg =  'There are multiple spin systems for residue %s.\n?' % name
              msg += 'Merge spin systems together?'
              
              if showOkCancel('Confirm',msg):
                spinSystem = mergeSpinSystems(spinSystems[i],spinSystem)
              else:
                spinSystem = None

            else:
              txt1 = '%d%s' % (spinSystem.residue.seqCode,spinSystem.residue.ccpCode)
              txt2 = '%d%s' % (spinSystems[i].residue.seqCode,spinSystems[i].residue.ccpCode)
              msg  = 'Cannot set spin system for F%d dim' % (i+1)
              msg += 'Offset %d causes conflict between %s and %s' % (seqOffset, txt1, txt2)
              showWarning('Failure',msg)
              return
         
        if resonanceList[i]:
          nmrProject = resonanceList[i][0].nmrProject
          if not spinSystem:
            if spinSystems[i]:
              spinSystem = spinSystems[i]
            else:
              spinSystem = nmrProject.newResonanceGroup()
          
          makeSeqSpinSystemLink(ref, spinSystem, seqOffsets[i])
          
          for resonance in resonanceList[i]:
            if resonance.resonanceGroup is not spinSystem:
              addSpinSystemResonance(spinSystem,resonance)
              

def getSpinSystemResidues(spinSystem):
  """
  Find a residues that spin system is assigned to, or tentatively assigned to.

  .. describe:: Input
  
  Nmr.ResonanceGroup
  
  .. describe:: Output
  
  list of MolSystem.Residue. Contains no duplicates
  """
  
  residue = spinSystem.residue
  if residue is not None:
    result = [residue]
  
  else:
    result = [x.possibility for x in spinSystem.sortedResidueProbs() 
               if x.weight]
  #
  return result
              

def getSpinSystemChemComps(spinSystem):
  """
  Find a ChemComps that spin system is assigned to, or tentatively assigned to.

  .. describe:: Input
  
  Nmr.ResonanceGroup  
  
  .. describe:: Output
  
  list of MolSystem.Residue. Contains no duplicates
  """
  
  chemComp = spinSystem.chemComp
  if chemComp is not None:
    result = [chemComp]
  
  else:
    result = [x.possibility for x in spinSystem.sortedResidueTypeProbs() 
               if x.weight]
  #
  return result


def initAtomSetMappings(analysisProject):
  """init/refresh AtomSetMappings
  """  

  atomSetDict = {}
  for atomSet in analysisProject.nmrProject.atomSets:
    atomSetDict[atomSet.serial] = atomSet

  msDict = {}
  for molSystem in analysisProject.root.molSystems:
    msDict[molSystem.code] = {}
    
    for chain in molSystem.chains:
      msDict[molSystem.code][chain.code] = {}
      
      for residue in chain.residues:
        msDict[molSystem.code][chain.code][residue.seqId] = 1

  for chainMapping in analysisProject.chainMappings:
    msc = chainMapping.molSystemCode
    chc = chainMapping.chainCode
    
    if msDict.get(msc) is None:
      chainMapping.delete()
      continue
    
    if msDict[msc].get(chc) is None:
      chainMapping.delete()
      continue

    for residueMapping in chainMapping.residueMappings:
      if msDict[msc][chc].get(residueMapping.seqId) is None:
        residueMapping.delete()
        continue

      for atomSetMapping in residueMapping.atomSetMappings:
        if atomSetMapping.atomSetSerials:
          atomSets = []

          for serial in atomSetMapping.atomSetSerials:
            atomSet = atomSetDict.get(serial)

            if atomSet:
              atomSets.append(atomSet)

          if not atomSets:
            continue

          updateAtomSetMapping(atomSetMapping, atomSets)


def estimateAssignmentTolerances(dataSource, defPoints=1.0, minTol=0.02):
  """ get a list of assignment tolerances, in dimension order.
  Dimensions that cannot be give a tolerance have toelrance set to 0
  """
  
  tolerances = []
  for ii,dataDim in enumerate(dataSource.sortedDataDims()):
    useTol = 0
    if hasattr(dataDim,'analysisDataDim'):
      # Set tolerance from AnalysisDataDim
      add = dataDim.analysisDataDim
      tol = add.assignTolerance
      if tol != 0.1:
        # 0.1 is the default value, and so not reliable)
        useTol = tol
    
    if not useTol:
      # Set tolerance from point separation
      dataDimRefs = [x for x in dataDim.dataDimRefs 
                     if x.expDimRef.measurementType.lower() == 'shift']
      if len(dataDimRefs) == 1:
        useTol = max(minTol, dataDimRefs[0].valuePerPoint * defPoints)
    tolerances.append(useTol)
  #
  return tolerances

def getPeakDimResonances(peakDim):
  """
  Find all resonances that a peakDim points to.

  .. describe:: Input

  Nmr.PeakDim

  .. describe:: Output

  set of Resonances
  """

  resonances = set([peakDimContrib.resonance for peakDimContrib in peakDim.peakDimContribs])

  return resonances

def getAtomResonances(atom):
  """
  Find all resonances that an atom points to.

  .. describe:: Input

  MolSystem.Atom

  .. describe:: Output

  set of Resonances
  """

  resonances = set()
  atomSet = atom.atomSet
  if atomSet:
    for resonanceSet in atomSet.resonanceSets:
      resonances.update(resonanceSet.resonances)
  return resonances


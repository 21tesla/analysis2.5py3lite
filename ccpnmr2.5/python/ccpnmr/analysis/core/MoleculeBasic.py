LICENSE = """
======================COPYRIGHT/LICENSE START==========================

MoleculeBasic.py: Part of the CcpNmr Analysis program

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

from memops.general.Util import copySubTree

try:
  from memops.gui.MessageReporter import showWarning, showOkCancel, showYesNo
except ImportError:
  from memops.universal.MessageReporter import showWarning, showOkCancel, showYesNo

from ccp.general.ChemCompOverview import chemCompStdDict

from ccpnmr.analysis.core.AssignmentBasic import assignAtomsToRes, assignResonanceResidue, getResidueResonances

from ccp.util.Molecule import makeChain
from ccp.util.NmrExpPrototype import longRangeTransfers
from ccp.util.LabeledMolecule import getIsotopomerSingleAtomFractions, getIsotopomerAtomPairFractions
from ccp.util.LabeledMolecule import singleAtomFractions, atomPairFractions

from ccp.lib.MoleculeAlign import getChainResidueMapping, getSequenceResidueMapping
# The following should not be necessary:
from ccp.lib.MoleculeAlign import findMatchingChains
# The following really are not necessary (and have been renamed):
# ccp.lib.MoleculeAlign import sequenceAlign, substringAlign, dynamicAlign

from ccp.lib.MoleculeModify import copyMolecule

from ccp.lib import MoleculeQuery

from ccp.lib.MoleculeQuery import getNumConnectingBonds, areAtomsTocsyLinked
from ccp.lib.MoleculeQuery import getBoundAtoms, areAtomsBound, areResonancesBound
from ccp.lib.MoleculeQuery import getMolTypeCcpCodes, getChemCompOverview
from ccp.lib.MoleculeQuery import getResidueCode, getLinkedResidue
from ccp.lib.MoleculeQuery import getResidueObservableAtoms, getAtomsTorsion
from ccp.lib.MoleculeQuery import getChemicalShiftBounds, getRandomCoilShift
from ccp.lib.MoleculeQuery import getRandomCoilShiftCorrectionsDict, getChemAtomNmrRef
from ccp.lib.MoleculeQuery import makeResidueLocTag
from ccp.lib.MoleculeQuery import greekSortAtomNames, getUnicodeAtomName
from ccp.lib.MoleculeQuery import _getUnicodeGreek as getUnicodeGreek

from ccp.lib.MoleculeQuery import STANDARD_ISOTOPES, DEFAULT_ISOTOPES
from ccp.lib.MoleculeQuery import PROTEIN_RESIDUE_CLASS_DICT, unicodeGreek

#STANDARD_ISOTOPES = ['1H','13C','15N','31P','2H','29Si','19F','17O', '79Br']
# 1H, 2H, 3He, 13C, 15N, 14N, 19F, 31P, 17O, 10B, 11B, 35Cl,
# 37Cl, 43Ca, 195Pt 6Li, 7Li, 9Be, 21Ne, 23Na, 25Mg, 27Al,
# 29Si, 33S, 39K, 40K,  41K, 45Sc, 47Ti, 49Ti, 50V, 51V,
# 53Cr, 55Mn, 57Fe, 59Co,  61Ni, 63Cu, 65Cu, 67Zn, 69Ga,
# 71Ga, 73Ge, 77Se, 81Br, 87Rb, 87Sr, 95Mo, 109Ag, 113Cd,
# 125Te, 127I, 133Cs, 135Ba, 137Ba, 139La, 183W, 199Hg

#DEFAULT_ISOTOPES = {'H':'1H','C':'13C','N':'15N','P':'31P','Si':'29Si',
#                    'F':'19F','O':'16O', 'Br':'79Br',}


# Should really be derived or modelled attrs
#PROTEIN_RESIDUE_CLASS_DICT = {'Acidic'       :['Asp','Glu'],
#                              'Basic'        :['Arg','Lys','His'],
#                              'Charged'      :['Asp','Glu','Arg','Lys','His'],
#                              'Polar'        :['Asn','Gln','Asp','Glu','Arg','Lys','His','Ser','Thr','Tyr'],
#                              'Non-polar'    :['Ala','Phe','Gly','Ile','Leu','Met','Pro','Val','Trp','Cys'],
#                              'Hydrophilic'  :['Ser','Asp','Glu','Arg','Lys','His','Asp','Glu','Pro','Tyr'],
#                              'Hydrophobic'  :['Phe','Met','Ile','Leu','Val','Cys','Trp','Ala','Thr','Gly'],
#                              'Amide'        :['Asn','Gln'],
#                              'Hydroxyl'     :['Ser','Thr','Tyr'],
#                              'Aromatic'     :['Phe','Ptr','Tyr','Trp'],
#                              'Beta-branched':['Thr','Val','Ile'],
#                              'Small'        :['Cys','Ala','Ser','Asp','Thr','Gly','Asn'],
#                              'Neutral'      :['Ala','Asn','Cys','Gln','Gly','Ile','Leu','Met',
#                                               'Phe','Pro','Ser','Thr','Trp','Tyr','Val'],
#                              'Methyl'       :['Ala','Met','Ile','Leu','Thr','Val'],
#                             }

# X is an unusual base, not ambiguiuty

STEREO_PREFIX = 'stereo_'
CARBOHYDRATE_MOLTYPE = 'carbohydrate'
PROTEIN_MOLTYPE = 'protein'
OTHER_MOLTYPE = 'other'
DNA_MOLTYPE = 'DNA'
RNA_MOLTYPE = 'RNA'
DNARNA_MOLTYPE = 'DNA/RNA'

########################################################################
#
# Functions scheduled for move to ccp/lib/MoleculeQuery
#


"""
userResidueCodesDict = {DNA_MOLTYPE:{'A':'Ade','T':'Thy','G':'Gua','C':'Cyt','U':'Ura'},
                        RNA_MOLTYPE:{'A':'Ade','T':'Thy','G':'Gua','C':'Cyt','U':'Ura','I':'Ino'},
                        PROTEIN_MOLTYPE:{},
                        CARBOHYDRATE_MOLTYPE:{}
                        }
"""
userResidueCodesDict = MoleculeQuery.userResidueCodesDict

def setUserResidueCodesDict(resDict):

  global userResidueCodesDict
  MoleculeQuery.userResidueCodesDict = resDict
  userResidueCodesDict = resDict
      
# NOT Currently Used.
def getLinkAtoms(residueA, linkCode='prev'):

  chain = residueA.chain
  molResidue = residueA.molResidue
  
  molResLinkEndA = molResidue.findFirstMolResLinkEnd(linkCode=linkCode)
  chemAtomA = molResLinkEndA.linkEnd.boundChemAtom

  atomA = residueA.findFirstAtom(chemAtom=chemAtomA)
  atomB = None
  
  if molResLinkEndA and molResLinkEndA.molResLink:
  
    for molResLinkEndB in molResLinkEndA.molResLink.molResLinkEnds:
    
      if molResLinkEndB is not molResLinkEndA:
      
        residueB = chain.findFirstResidue(molResidue=molResLinkEndB.molResidue)
        chemAtomB = molResLinkEndB.linkEnd.boundChemAtom
        atomB = residueB.findFirstAtom(chemAtom=chemAtomB)
        
        break
  
  return atomA, atomB

      
   

def makeGuiName(name, elementSymbol, molType):
  """Convert atom or atomSet name into name for gui: e.g H becomes Hn
  .. describe:: Input
  
  Word (Nmr.AtomSet.name), Word, Word

  .. describe:: Output
  
  Word 
  """
  return MoleculeQuery.makeGuiName(name, elementSymbol)

  
  
  
####################################################################
#
# Assignment-related - consider move to ccp/lib/AssignmentLib
#
  
  
  
def transferChainAssignments(chainA, chainB):
  """Transfer any atomic assignments from one chain to another where possible.
  .. describe:: Input
  
  MolSystem.Chain, MolSystem.Chain

  .. describe:: Output
  
  None
  """

  mapping = getChainResidueMapping(chainA, chainB)
  for residueA, residueB in mapping:
    if residueB:
      resonancesB = getResidueResonances(residueB)
      if resonancesB:
        msg  = 'Destination residue %s%d has assignments. Continue?.'
        data = (residueB.seqCode,residueB.ccpCode)
        if not showOkCancel('Warning', msg % data):
          return
 
  for residueA, residueB in mapping:
    if residueA:
      if residueB is None:
        msg = 'Residue %d%s has no equivalent in destination chain'
        data = (residueA.seqCode,residueA.ccpCode)
        showWarning('Warning', msg % data)
      else:
        transferResidueAssignments(residueA,residueB)


def transferResidueAssignments(residueA,residueB):
  """
  Move any atomic assignments from one residue to another where possible.
  
  .. describe:: Input
  
  MolSystem.Residue, MolSystem.Residue

  .. describe:: Output
  
  None
  """
  
  resonancesA = getResidueResonances(residueA)

  for resonance in resonancesA:
    assignResonanceResidue(resonance, residueB)

def duplicateResidueAssignments(residueA, residueB, experimentChains=None):
  """
  Assign residueB to an equivalent set of resonances as residueA
  Optional dictionary, keyed by experiment to specifiy which chain (should be parent of residueA
  or residueB or None - for both) an experiment's peak assignments should go with.
 
  .. describe:: Input
  
  MolSystem.Residue, MolSystem.Residue, Dict of Nmr.Experiment:MolSystem.Chain (or None)

  .. describe:: Output
  
  None
  """
  
  from ccpnmr.analysis.core.AssignmentBasic import assignResToDim
  
  project = residueA.root
  nmrProject = project.currentNmrProject
  
  atomSetDict = {}
  for atom in residueB.atoms:
    atomSet = atom.atomSet
    if atomSet:
      atomSetDict[atomSet.name] = atomSet
  
  if experimentChains is None:
    experimentChains = {}
  else:
    for experiment in experimentChains:
      chain = experimentChains[experiment]
      molSystems = experiment.molSystems
      if chain and chain.molSystem not in molSystems:
        experiment.addMolSystem(chain.molSystem)
      else:
        if residueA.chain.molSystem not in molSystems:
          experiment.addMolSystem(residueA.chain.molSystem)
        if residueB.chain.molSystem not in molSystems:
          experiment.addMolSystem(residueB.chain.molSystem)
  
  resonancesA = getResidueResonances(residueA)
  resonancesB = getResidueResonances(residueB)

  resonancesDictA = {}
  for resonance in resonancesA:
    atomSetNames = [ass.name for ass in resonance.resonanceSet.atomSets]
    atomSetNames.sort()
    key = tuple(atomSetNames)
    resonancesDictA[key] = resonancesDictA.get(key, []) + [resonance,]

  resonancesDictB = {}
  for resonance in resonancesB:
    atomSetNames = [ass.name for ass in resonance.resonanceSet.atomSets]
    atomSetNames.sort()
    key = tuple(atomSetNames)
    resonancesDictB[key] = resonancesDictB.get(key, []) + [resonance,]

 
  transfers = []
  pureTransfer = True
  for key in resonancesDictA.keys():
    resonances1 = resonancesDictA[key]
    resonances2 = resonancesDictB.get(key, [])
    
    atomSets = []
    for atomSetName in key:
      if atomSetDict.has_key(atomSetName):
        atomSets.append(atomSetDict[atomSetName])
    
    if len(atomSets) == len(key): # If we find all equiv atom sets in destination (not always true if different res type)

      # Duplications of resonances is only done if needed at a per atomSet level
      doDuplicate = False
      contribs = set()
      
      for resonance1 in resonances1:
        contribs.update(resonance1.peakDimContribs)
        
      if contribs:
        for contrib in contribs:
          experiment = contrib.peakDim.peak.peakList.dataSource.experiment
          if experimentChains.get(experiment, residueB.chain) is not residueB.chain: # Not a plain transfer
            doDuplicate = True
            pureTransfer = False
            break
      
      else: # All orphans
        for experiment in experimentChains:
          if experimentChains[experiment] is not residueB.chain:
            doDuplicate = True
            pureTransfer = False
            break
             
               
      if doDuplicate: # Not a plain transfer Need to make new resonances
        while len(resonances2) < len(resonances1):
          resonances2.append(None)
        

        for i, resonance1 in  enumerate(resonances1):
          resonance2  = resonances2[i]
          peakDimContribs = list(resonance1.peakDimContribs)
          
          if peakDimContribs:
            for peakDimContrib in peakDimContribs:
              peakDim    = peakDimContrib.peakDim
              experiment = peakDim.peak.peakList.dataSource.experiment
              chainB     = experimentChains.get(experiment)
 
              if chainB is None: # Duplication on this peak
                if resonance2 is None:
                  resonance2 = nmrProject.newResonance(isotopeCode=resonances1[0].isotopeCode)
                  assignAtomsToRes(atomSets, resonance2)
                  resonances2[i] = resonance2
 
                assignResToDim(peakDim, resonance2, doWarning=False)
 
              elif chainB is residueB.chain: # Transfer this peak
                if resonance2 is None:
                  resonance2 = nmrProject.newResonance(isotopeCode=resonances1[0].isotopeCode)
                  assignAtomsToRes(atomSets, resonance2)
                  resonances2[i] = resonance2
 
                assignResToDim(peakDim, resonance2, doWarning=False)
                for peakContrib in peakDimContrib.peakContribs:
                  if len(peakContrib.peakDimContribs) < 2:
                    peakContrib.delete()
                
                if not peakDimContrib.isDeleted:
                  peakDimContrib.delete()
           
          else:
            if resonance2 is None:
              resonance2 = nmrProject.newResonance(isotopeCode=resonances1[0].isotopeCode)
              assignAtomsToRes(atomSets, resonance2)
              resonances2[i] = resonance2
             
              
      else: 
        # Plain transfer, move existing resonances to different atoms
        transfers.append((resonances1, atomSets))

  # Transfer original spin system to destination residue
  # to keep existing links and attrs (pure transfers only)
  spinSystem = nmrProject.findFirstResonanceGroup(residue=residueA)
  if spinSystem and pureTransfer:
    spinSystem.residue = residueB
    spinSystem.ccpCode = residueB.ccpCode
    spinSystem.chains = [residueB.chain,]
  
  # Clear all resonance to spinSystem links first so
  # that assigniment can pickup the original spin system
  for resonances, atomSets in transfers:
    for resonance in resonances:
      resonance.resonanceGroup = None
    
  for resonances, atomSets in transfers:
    for resonance in resonances:
      assignAtomsToRes(atomSets, resonance)
   
  
  
  
####################################################################
#
# Analysis-specific (uses Analysis data model)
#
  
        
def makeMolSystemLink(residueA,residueB,linkEndA,linkEndB):
  """Make a molSystemLink given two residues and the linkEnds to be joined
  .. describe:: Input
  
  MolSystem.Residue, MolSystem.Residue, ChemComp.LinkEnd, ChemComp.LinkEnd

  .. describe:: Output
  
  MolSystem.MolSystemLink
  """
  molSystem = residueA.chain.molSystem  
  molSysLinkEndA = residueA.newMolSystemLinkEnd(linkCode=linkEndA.linkCode)
  molSysLinkEndB = residueB.newMolSystemLinkEnd(linkCode=linkEndB.linkCode)
  molSysLinkEnds = (molSysLinkEndA,molSysLinkEndB)
  
  try:
    molSystemBond = molSystem.newMolSystemLink(molSystemLinkEnds=molSysLinkEnds) 
  except Exception:
    return

  removeAtoms = []
  for atom in removeAtoms:
    atomSet = atom.atomSet
    if atomSet:
      for resonanceSet in atomSet.resonanceSets:
        resonanceSet.delete()
      atomSet.delete()
    atom.delete()

  for residue in (residueA,residueB):
    residueMapping = getResidueMapping(residue)
    for atomSetMapping in residueMapping.atomSetMappings:
      atomSetMapping.delete()
      
    makeResidueAtomSets(residue)

  return molSystemBond



def makeAtomSet(guiName,atoms,chemAtomSet,mappingType):
  """
  Make atomSet and atomSetMapping for equivalent atoms
  
  .. describe:: Input
  
  Word (AtomSet.name), List of MolSystem.Atoms,
  ChemComp.ChemAtomSet, Word (AtomSetMapping.mappingType)
  
  .. describe:: Output
  
  Nmr.AtomSet
  """
  
  # RHFogh 3/12/09 - refactored to reduce getAtomSet calls
  
  atom0 = list(atoms)[0]
  project = atom0.root
  
  atomSets = [x.atomSet for x in atoms]
  atomSet0 = atomSets[0]
  aSet = set(atomSets)
  if len(aSet) != 1:
    for atomSet in aSet:
      if atomSet and not atomSet.resonanceSets:
        atomSet.delete()
  
  nmrProject = project.currentNmrProject
  
  if atomSet0 is None:
    atomSet = nmrProject.newAtomSet(atoms=atoms)
  else:
    atomSet = atomSet0
 
  residue = atom0.residue
    
  residueMapping = getResidueMapping(residue)
  if not residueMapping.findFirstAtomSetMapping(name=guiName):
    makeAtomSetMapping(residueMapping,guiName,(atomSet,),chemAtomSet,mappingType)

  atomSet.name = guiName
  return atomSet

def getResidueMapping(residue, aromaticsEquivalent=True):
  """Gives the Analysis.ResidueMapping for a residue
             Makes a new one with new AtomSetsMappings if not exists
             Makes a ChainMapping too if needed.
  .. describe:: Input
  
  MolSystem.Residue

  .. describe:: Output
  
  Analysis.ResidueMapping
  """
  
  if hasattr(residue, 'residueMapping'):
    return residue.residueMapping
  
  residueMapping = None
  analysisProject = residue.root.currentAnalysisProject
  chainMapping = analysisProject.findFirstChainMapping(chain=residue.chain)

  if not chainMapping:
    chainMapping = analysisProject.newChainMapping(molSystemCode=residue.chain.molSystem.code,
                                                   chainCode=residue.chain.code)
    chainMapping.residueMappingDict = {}
  else:
    if not hasattr(chainMapping, 'residueMappingDict'):
      chainMapping.residueMappingDict = {}
    
    residueMapping = chainMapping.residueMappingDict.get(residue.seqId)
    
    if not residueMapping:
      residueMapping = chainMapping.findFirstResidueMapping(seqId=residue.seqId)
      chainMapping.residueMappingDict[residue.seqId] = residueMapping
      
  if not residueMapping:
    residueMapping = chainMapping.newResidueMapping(seqId=residue.seqId)
    #makeResidueAtomSets(residue)

  residue.residueMapping = residueMapping 
  
  if not residueMapping.atomSetMappings:
    makeResidueAtomSets(residue, aromaticsEquivalent=aromaticsEquivalent)

  return residueMapping

def makeAtomSetMapping(residueMapping,name,atomSets,chemAtomSet,mappingType,resonances=None):
  """Make atomSetMapping given atomSets and mapping type
  .. describe:: Input
  
  Analysis.ResidueMapping, Word, MolSystem.Residue,
             List of Nmr.AtomSets, ChemComp.ChemAtomSet,
             Word, Word, List of Nmr.Resonances

  .. describe:: Output
  
  Analysis.AtomSetMapping
  """

  atom          = list(atomSets)[0].findFirstAtom()
  elementSymbol = atom.chemAtom.elementSymbol
  serials       = []
  for atomSet in atomSets:
    serials.append(atomSet.serial)

  molType = residueMapping.residue.molResidue.molType    
  guiName = makeGuiName(name, elementSymbol, molType)
      
  atomSetMapping = residueMapping.newAtomSetMapping(name=guiName,mappingType=mappingType,
                                                    atomSetSerials=serials,
                                                    chemAtomSet=chemAtomSet,
                                                    elementSymbol=elementSymbol)
                                                    
  if resonances is not None:
    resSerials = []
    for resonance in resonances:
      resSerials.append(resonance.serial)
    atomSetMapping.setResonanceSerials(resSerials)

  return atomSetMapping

def makeResidueAtomSetsNonEquivalent(residue):
  """Remake a residue's atom sets if they are found to be non-equivalent: e.g.
             if an aromatic ring does not rotate quickly on the NMR timescale
  .. describe:: Input
  
  MolSystem.Residue

  .. describe:: Output
  
  None
  """
  
  residueMapping = getResidueMapping(residue)
  #chain   = residue.chain
  molType = residue.molResidue.molType
  nonequivalent     = {}
  elementSymbolDict = {}
  chemAtomSetDict   = {}
  for atom in residue.atoms:
    chemAtom = atom.chemAtom
    chemAtomSetDict[atom] = chemAtom
    elementSymbol = chemAtom.elementSymbol
    if chemAtom.chemAtomSet:
      chemAtomSet = chemAtom.chemAtomSet
      name = chemAtomSet.name
      if chemAtomSet.isEquivalent is None: # i.e. not False, aromatic rotation
        if nonequivalent.get(name) is None:
          nonequivalent[name] = []
        nonequivalent[name].append(atom)
        elementSymbolDict[name] = chemAtom.elementSymbol
        chemAtomSetDict[name] = chemAtomSet

  for groupName in nonequivalent.keys():
    atoms = nonequivalent[groupName]
    atomSet = atoms[0].atomSet
    if atomSet:
      for atom in atoms[1:]:
        if atom.atomSet is not atoms[0].atomSet:
          return # already seperate
    
    elementSymbol = elementSymbolDict[groupName]
    chemAtomSet   = chemAtomSetDict[groupName]
    resonances = []

    if atomSet:
      for resonanceSet in atomSet.resonanceSets:
        for resonance in resonanceSet.resonances:
          resonances.append(resonance)
    
      atomSet.delete()
      name = makeGuiName(groupName, elementSymbol, molType)
      atomSetMapping = residueMapping.findFirstAtomSetMapping(name=name)
      if atomSetMapping:
        atomSetMapping.delete()

    atomSets = []
    atomSetNames = []
    for atom in atoms:
      name = chemAtomSetDict[atom].name
      atomSet = makeAtomSet(name,(atom,),chemAtomSet,'stereo')
      atomSets.append(atomSet)
      atomSetNames.append(name)

    resonanceSet = None
    for resonance in resonances:
      resonanceSet = assignAtomsToRes(atomSets,resonance,resonanceSet)

    for n, atom in enumerate(atoms):
      name = chemAtomSetDict[atom].name
      name2 = makeNonStereoName(molType, name, n)
      makeGuiMultiAtomSet(residue, name2, atomSetNames,elementSymbol,'nonstereo',chemAtomSet)

    makeGuiMultiAtomSet(residue, groupName, atomSetNames,elementSymbol,'ambiguous',chemAtomSet)

def makeResidueAtomSetsEquivalent(residue):
  """Remake a residue's atom sets if they are found to be equivalent: e.g.
             if an aromatic ring rotates quickly on the NMR timescale
  .. describe:: Input
  
  MolSystem.Residue

  .. describe:: Output
  
  None
  """

  #chain   = residue.chain
  molType = residue.molResidue.molType
  residueMapping = getResidueMapping(residue)
  equivalent        = {}
  elementSymbolDict = {}
  chemAtomSetDict   = {}
  for atom in residue.atoms:
    chemAtom = atom.chemAtom
    chemAtomSetDict[atom] = chemAtom
    elementSymbol = chemAtom.elementSymbol
    if chemAtom.chemAtomSet:
      chemAtomSet = chemAtom.chemAtomSet
      name = chemAtomSet.name
      if chemAtomSet.isEquivalent is None: # i.e. not False, aromatic rotation
        if equivalent.get(name) is None:
          equivalent[name] = []
        equivalent[name].append(atom)
        elementSymbolDict[name] = chemAtom.elementSymbol
        chemAtomSetDict[name] = chemAtomSet
         
  for groupName in equivalent.keys():
    atoms = equivalent[groupName]
    if atoms[0].atomSet:
      for atom in atoms[1:]:
        if atom.atomSet is atoms[0].atomSet:
          return
    
    elementSymbol = elementSymbolDict[groupName]
    chemAtomSet   = chemAtomSetDict[groupName]
    resonances = []
    for atom in atoms:
      # TBD more nested layers?
      # delete ambiguous
      name = makeGuiName(chemAtomSet.name, elementSymbol, molType)
      atomSetMapping = residueMapping.findFirstAtomSetMapping(name=name)
      if atomSetMapping:
        atomSetMapping.delete()

      # delete stereospecific
      name = makeGuiName(atom.chemAtom.name, elementSymbol, molType)
      atomSetMapping = residueMapping.findFirstAtomSetMapping(name=name)
      if atomSetMapping:
        atomSetMapping.delete()
        
      # delete non-stereospecific
      name = makeGuiName(makeNonStereoName(molType, atom.chemAtom.name), elementSymbol, molType)
      atomSetMapping = residueMapping.findFirstAtomSetMapping(name=name)
      if atomSetMapping:
        atomSetMapping.delete()
        
      atomSet = atom.atomSet
      if atomSet:
        for resonanceSet in atomSet.resonanceSets:
          resonances.extend(resonanceSet.resonances)
        atomSet.delete()
   
    elementSymbol = elementSymbolDict[groupName]
    chemAtomSet = chemAtomSetDict[groupName]
    # make single equivalent group
    #makeAtomSet(,groupName,atoms,chemAtomSet,'simple')
    atomSet = makeAtomSet(groupName,atoms,chemAtomSet,'simple')
    for resonance in resonances:
      assignAtomsToRes([atomSet,],resonance)
    

def makeResidueAtomSets(residue, aromaticsEquivalent=True):
  """Make all atomSets and atomSetMappings for a given residue
             Aromatic Phe, Tyr (Hd1,Hd2), (He1,He2) can be made into 
             single equivalent atom sets due to rotation.
  .. describe:: Input
  
  MolSystem.Residue, Boolean

  .. describe:: Output
  
  None
  """
  
  getResidueMapping(residue)
  
  equivalent = {}
  elementSymbolDict = {}
  nonequivalent = {}
  multiSet = {}
  chemAtomSetDict = {}
  inMultiSet = {}
  molType = residue.molResidue.molType
  
  for atom in residue.atoms:  
    chemAtom = atom.chemAtom
    chemAtomSetDict[atom] = chemAtom
    elementSymbol = chemAtom.elementSymbol
    chemAtomSet = chemAtom.chemAtomSet

    if chemAtomSet is None:
      name = chemAtom.name
      makeAtomSet(name,(atom,),None,'simple')
      
    else:
      name = chemAtomSet.name
      elementSymbolDict[name] = elementSymbol
      chemAtomSetDict[name] = chemAtomSet
      if chemAtomSet.isEquivalent:
        if equivalent.get(name) is None:
          equivalent[name] = []
        equivalent[name].append(atom)
        
      elif (chemAtomSet.isEquivalent is None) and atom.atomSet and (len(atom.atomSet.atoms) > 1):
        # aromatic rotation prev set
        if equivalent.get(name) is None:
          equivalent[name] = []
        equivalent[name].append(atom)
           
      elif (chemAtomSet.isEquivalent is None) and (not atom.atomSet) and aromaticsEquivalent:
        # aromatic rotation to be set
        if equivalent.get(name) is None:
          equivalent[name] = []
        equivalent[name].append(atom)
          
      else:
        if nonequivalent.get(name) is None:
          nonequivalent[name] = []
        nonequivalent[name].append(atom)
   
      if chemAtomSet.chemAtomSet is not None:
        multiName = chemAtomSet.chemAtomSet.name
        chemAtomSetDict[multiName] = chemAtomSet.chemAtomSet
        elementSymbolDict[multiName] = elementSymbol
        if multiSet.get(multiName) is None:
          multiSet[multiName] = {}
        multiSet[multiName][name] = 1
        inMultiSet[name] = multiName

  for groupName in equivalent.keys():
    atoms = equivalent[groupName]
    elementSymbol = elementSymbolDict[groupName]
    chemAtomSet = chemAtomSetDict[groupName]
    if len(atoms)==2:
      # not enough atoms for multi sets!
      makeAtomSet(groupName,atoms,chemAtomSet,'simple')
    else:
      if inMultiSet.get(groupName):
        # e.g. for Val Hg1*
        makeAtomSet(groupName,atoms,chemAtomSet,'stereo')
 
      else:
        makeAtomSet(groupName,atoms,chemAtomSet,'simple')

  for groupName in nonequivalent.keys():
    atoms = nonequivalent[groupName]
    elementSymbol = elementSymbolDict[groupName]
    chemAtomSet = chemAtomSetDict[groupName]
    atomSetNames = []
    
    if len(atoms) == 1:
      atom = atoms[0]
      # not enough atoms for prochiral. Corrupt ChemComp
      makeAtomSet(atom.name, atoms, None, 'simple')
      continue
      
    for atom in atoms:
      name = chemAtomSetDict[atom].name
      makeAtomSet(name,(atom,),chemAtomSet,'stereo')
      atomSetNames.append(name)

    for n, atom in enumerate(atoms):
 
      #name = chemAtomSetDict[atom].name
      #name2 = makeNonStereoName(molType, name, n)
      # Shouldn't have to do this if non-equiv groups have paired names
      
      name2 = makeNonStereoName(molType, '%s%d' % (chemAtomSet.name[:-1], n), n)
        
      makeGuiMultiAtomSet(residue, name2, atomSetNames,
                          elementSymbol,'nonstereo',chemAtomSet)

    makeGuiMultiAtomSet(residue, groupName, atomSetNames,
                        elementSymbol,'ambiguous',chemAtomSet)

  for groupName in multiSet.keys():
    atomSetNames  = multiSet[groupName].keys()
    elementSymbol = elementSymbolDict[groupName]
    chemAtomSet   = chemAtomSetDict[groupName]
    if "|" in groupName:
      # we don't do these pseudoatoms in Analysis
      continue

    # e.g. for Val Hga*
    for n, atomSetName in enumerate(atomSetNames):
      name2 = makeNonStereoName(molType, atomSetName, n)
      makeGuiMultiAtomSet(residue, name2, atomSetNames,
                          elementSymbol,'nonstereo',chemAtomSet)
    
    makeGuiMultiAtomSet(residue, groupName, atomSetNames,
                        elementSymbol,'ambiguous',chemAtomSet)

def makeNonStereoName(molType, name, n=None):
  """Convert a sterospecific atom name into a non-stereospecific one for a GUI
  .. describe:: Input
  
  Word, Int (naming offset from start of alphabet)

  .. describe:: Output
  
  Word
  """

  match = re.match('(\w+)(\d|\'+)(\D*)', name)
  
  if not match:
    #print molType, name, n
    return 
  
  
  letters = match.group(1)
  number  = match.group(2)
  prime   = ''
  
  if number == '\'':
    number = 1
    prime = '\''
  elif number == '\'\'':
    number = 2
    prime = '\''
  
  if n is None:
    n = int(number) - 1

  name = letters + prime + chr(ord('a')+n)+ match.group(3)
    
  return name

def makeGuiMultiAtomSet(residue,multiGuiName,guiSetsNames,elementSymbol,mappingType,chemAtomSet):
  """Make atom set mappings for multiple atom set selections
  .. describe:: Input
  
  MolSystem.Residue, Word (Analysis.AtomSetMapping.name),
             List of Words (Analysis.AtomSetMapping.names),
             Word, Word, ChemComp.ChemAtomSet

  .. describe:: Output
  
  Analysis.AtomSetMapping
  """
  
  if "|" in multiGuiName:
    return
  
  residueMapping = getResidueMapping(residue)  
  molType = residue.molResidue.molType
  for guiName in guiSetsNames:
    atomSetMapping = residueMapping.findFirstAtomSetMapping(name=makeGuiName(guiName, elementSymbol, molType))
    if atomSetMapping is None:
      print "Non-existent group error in makeGuiMultiAtomSet for", residue.molResidue.ccpCode, residue.seqCode, guiName
      return
    #atomSet      = atomSetMapping.atomSets[0]
    chemAtomSet1 = atomSetMapping.chemAtomSet
    
    for guiName2 in guiSetsNames:
      atomSetMapping2 = residueMapping.findFirstAtomSetMapping(name=makeGuiName(guiName2, elementSymbol, molType))
      if atomSetMapping2 is None:
        print "Non-existent group error in makeGuiMultiAtomSet for", residue.molResidue.ccpCode
        return
      #atomSet      = atomSetMapping2.atomSets[0]
      chemAtomSet2 = atomSetMapping2.chemAtomSet
      if chemAtomSet2 and chemAtomSet1:
        if chemAtomSet1.isProchiral != chemAtomSet2.isProchiral:
          print "Prochiratity error in makeGuiMultiAtomSet for", residue.molResidue.ccpCode
          return
        if chemAtomSet1.isEquivalent != chemAtomSet2.isEquivalent:
          print "Equivalent error in makeGuiMultiAtomSet for ", residue.molResidue.ccpCode
          return

  atomSets = []
  for guiName in guiSetsNames:
    name0 = makeGuiName(guiName, elementSymbol, molType)
    atomSetSerials = residueMapping.findFirstAtomSetMapping(name=name0).atomSetSerials
    for atom in residue.atoms:
      atomSet = atom.atomSet
      if atomSet:
        if atomSet.serial in atomSetSerials and atomSet not in atomSets:
          atomSets.append(atomSet)
          break
    
  if not residueMapping.findFirstAtomSetMapping(name=multiGuiName):
    atomSetMapping = makeAtomSetMapping(residueMapping, multiGuiName, atomSets,
                                        chemAtomSet, mappingType)
  
  return atomSetMapping

def moveMolSystemChain(chain, molSystem):
  """Moves a chain from one molSystem to another
  .. describe:: Input
  
  MolSystem.Chain, MolSystem.MolSystem

  .. describe:: Output
  
  None
  """

  from MergeObjects import mergeObjects

  if chain.molSystem is molSystem:
    return

  molecule = chain.molecule
  project  = molSystem.root
  residues = chain.sortedResidues()
  seq      = []

  for residue in residues:
    seq.append( residue.ccpCode )

  newChain = makeChain(molSystem,molecule)
  
  for chainMapping in project.currentAnalysisProject.chainMappings:
    if chainMapping.chain is chain:
      chainMapping.delete()
  
  residues1 = residues
  
  for i in range(len(seq)):
    residue1 = residues1[i]
    residue2 = newChain.sortedResidues()[i]
    for atom in residue1.atoms:

      atom2 = residue2.findFirstAtom(name=atom.name)
      if atom2:
        mergeObjects(atom, atom2)
      else:
        print "missing %d %s %s" % (residue2.seqCode,residue2.ccpCode,atom.name) 
        
    mergeObjects(residue1, residue2)
    
  mergeObjects(chain, newChain)
  
  
####################################################################
#
# Unused, deprecated, or obsolete
#

# Obsolete, unused, and deprecated
def findMatchingChain(molSystem, ccpCodes, excludeChains=None, molTypes=None, doWarning=True):
  """* * * Deprecated in light of findMathcingChains() * * * 
             Find the mol system chain that best matches the input ccpCodes
             (like three letter codes). Useful for trying to match structures
             to existing molecular data. Optional argument to specify which
             chains cannot be matched.
  .. describe:: Input
  
  MolSystem.MolSystem, List of Strings (MolSystem.MolResidue.ccpCodes),
             List of MolSystem.Chains

  .. describe:: Output
  
  Tuple of (MolSystem.Chain, Int (index of first matching ccpCode),
             Int (first matching MolSystem.Residue.seqId))
  """
  print ("DEPRECATED, function findMatchingChain should not be used")

  chains = []
  for chain in molSystem.sortedChains():
    if excludeChains and (chain in excludeChains):
      continue
    
    sequence = []
    if not chain.residues:
      continue
    
    chains.append(chain)
    
    residues = []
    for residue in chain.sortedResidues():
      molType = residue.molType
      if molTypes and (molType not in molTypes): 
        break
      
      residues.append( (residue.seqId, residue) )
    
    else:  
      residues.sort()
 
      for residue in residues:
        sequence.append( residue[1].ccpCode )

      len0 = len(sequence)
      len1 = len(ccpCodes)
 
      if len0 < len1:
        continue
 
      elif len0 == len1:
        if sequence == ccpCodes:
          mapping = [(i, residues[i][1]) for i in range(len0)]
          return chain, mapping

        else:
          misMatch = 0
          for i in range(len0):
            if (ccpCodes[i] is not None) and (ccpCodes[i] != sequence[i]):
              misMatch = 1
              break
          if misMatch:
            continue
          else:
            mapping = [(i, residues[i][1]) for i in range(len0)]
            return chain, mapping
 
      else:
        d = len0 - len1
        for x in range(d+1):
          misMatch = 0
          for i in range(len1):
            if (ccpCodes[i] is not None) and (ccpCodes[i] != sequence[i+x]):
              misMatch = 1
              break
          if misMatch:
            continue
          else:
            mapping = [(i, residues[i+x][1]) for i in range(len1)]
            return chain, mapping
  
  scoreList = []
  
  bestMapping = None
  bestScore   = 0
  bestChain   = None

  for chain in chains:
    mapping, score = getSequenceResidueMapping(chain, ccpCodes)
    scoreList.append((score, chain, mapping))
    
  if scoreList:
    scoreList.sort()
    bestScore, bestChain, bestMapping = scoreList[-1]

  chain = None
  if bestScore and ( bestScore/float(len(bestChain.residues)) > 2.0 ):
    chain = bestChain
    
    if chain.molecule.molType in (CARBOHYDRATE_MOLTYPE): # DNA & RNA now hopefully fine to align
      chain = None
    
    if doWarning:
      msg  = 'Residue sequence matches an existing chain, but not exactly. '
      msg += 'Really link to chain %s? (Otherwise a new chain will be made)' % chain.code
      if not showYesNo('Query', msg):
        chain = None      

  return chain, bestMapping 


# DEPRECATED
def newMolSystem(project):
  """Get a new molSystem for a project with a unique code.
  .. describe:: Input
  
  Implementation.Project

  .. describe:: Output
  
  MolSystem.MolSystem
  """

  i = 1
  while project.findFirstMolSystem(code='MS%d' % (i)):
    i += 1
  molSystem = project.newMolSystem(code='MS%d' % (i))
  molSystem.name = molSystem.code
 
  return molSystem
  

# DEPRECATED
def findBoundResonances(resonance):
  """Find any resonances which are assigned to atoms covalently bound
             to the assigned atoms of the input resonance
  .. describe:: Input
  
  Nmr.Resonance

  .. describe:: Output
  
  List of Nmr.Resonances
  """

  from ccpnmr.analysis.core.AssignmentBasic import getBoundResonances
  
  return getBoundResonances(resonance)

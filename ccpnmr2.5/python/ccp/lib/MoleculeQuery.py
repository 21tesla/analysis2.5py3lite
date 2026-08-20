LICENSE = """
======================COPYRIGHT/LICENSE START==========================

MoleculeQuery.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2013 Wayne Boucher, Rasmus Fogh and Tim Stevens (University of Cambridge)

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

try:
  from memops.gui.MessageReporter import showWarning
except:
  from memops.universal.MessageReporter import showWarning

from ccp.general.ChemCompOverview import chemCompStdDict
from ccp.util.LabeledMolecule import getIsotopomerSingleAtomFractions, getIsotopomerAtomPairFractions
from ccp.util.LabeledMolecule import singleAtomFractions, atomPairFractions

from ccp.util.NmrExpPrototype import longRangeTransfers

STEREO_PREFIX = 'stereo_'
CARBOHYDRATE_MOLTYPE = 'carbohydrate'
PROTEIN_MOLTYPE = 'protein'
OTHER_MOLTYPE = 'other'
DNA_MOLTYPE = 'DNA'
RNA_MOLTYPE = 'RNA'
DNARNA_MOLTYPE = 'DNA/RNA'

userResidueCodesDict = {DNA_MOLTYPE:{'A':'Ade','T':'Thy','G':'Gua','C':'Cyt','U':'Ura'},
                        RNA_MOLTYPE:{'A':'Ade','T':'Thy','G':'Gua','C':'Cyt','U':'Ura','I':'Ino'},
                        PROTEIN_MOLTYPE:{},
                        CARBOHYDRATE_MOLTYPE:{}
                        }

STANDARD_ISOTOPES = ['1H','13C','14N','15N','31P','2H','29Si','19F','17O']

DEFAULT_ISOTOPES = {'H':'1H','C':'13C','N':'15N','P':'31P','Si':'29Si',
                    'F':'19F','O':'16O'}

# Should really be derived or modelled attrs
PROTEIN_RESIDUE_CLASS_DICT = {'Acidic'       :['Asp','Glu'],
                              'Basic'        :['Arg','Lys','His'],
                              'Charged'      :['Asp','Glu','Arg','Lys','His'],
                              'Polar'        :['Asn','Gln','Asp','Glu','Arg','Lys','His','Ser','Thr','Tyr'],
                              'Non-polar'    :['Ala','Phe','Gly','Ile','Leu','Met','Pro','Val','Trp','Cys'],
                              'Hydrophilic'  :['Ser','Asp','Glu','Arg','Lys','His','Asp','Glu','Pro','Tyr'],
                              'Hydrophobic'  :['Phe','Met','Ile','Leu','Val','Cys','Trp','Ala','Thr','Gly'],
                              'Amide'        :['Asn','Gln'],
                              'Hydroxyl'     :['Ser','Thr','Tyr'],
                              'Aromatic'     :['Phe','Ptr','Tyr','Trp'],
                              'Beta-branched':['Thr','Val','Ile'],
                              'Small'        :['Cys','Ala','Ser','Asp','Thr','Gly','Asn'],
                              'Neutral'      :['Ala','Asn','Cys','Gln','Gly','Ile','Leu','Met',
                                               'Phe','Pro','Ser','Thr','Trp','Tyr','Val'],
                              'Methyl'       :['Ala','Met','Ile','Leu','Thr','Val'],
                             }

# X is an unusual base, not ambiguiuty

####################################################################
#
# ChemComops and ccpCodes
#

def getMolTypeCcpCodes(molType='all', project=None):
  """Gives ccpCodes for chemComps according to molecule type: e.g. DNA
             Project can be input to search for non-standard types.
  .. describe:: Input
  
  Implementation.Project, String (ChemComp.molType or 'all')

  .. describe:: Output
  
  List of Words (ChemComp.CcpCodes)
  """

  ccpCodes = []
  if molType == 'all':
    molTypes = [PROTEIN_MOLTYPE, DNA_MOLTYPE, RNA_MOLTYPE,
                CARBOHYDRATE_MOLTYPE, OTHER_MOLTYPE]
  else:
    molTypes = [molType,]
    
  for molType in molTypes:
    chemCompDict = getChemCompOverview(molType, project)
    
    if chemCompDict:
      ccpCodes.extend( chemCompDict.keys() )
        
  if ccpCodes:
    ccpCodes.sort()
  
  return ccpCodes

def getChemCompOverview(molType, project=None):
  """Get a dictionary containing details of all available chemical compounds
             for a given molecule type. Project can be input to search for loaded,
             but non standard chem comps.
  .. describe:: Input
  
  Word (Molecule.MolResidue.MolType), Implementation.Project

  .. describe:: Output
  
  Dict of ChemComp.ccpCode:[Word, Word, Line, Word]
             (1-letter Code, 3-letter Code, name, mol formula)
  """
  
  if molType == OTHER_MOLTYPE:
    from ccp.general.ChemCompOverview import chemCompOverview
    chemCompDict = chemCompOverview.get(molType, {}) 
    
  else:
    chemCompDict = chemCompStdDict.get(molType, {}) 
   
  if project: 
    for chemComp in project.findAllChemComps(molType=molType):
      ccpCode  = chemComp.ccpCode
      
      if chemCompDict.get(ccpCode) is None:
        chemCompVar = chemComp.findFirstChemCompVar(linking=None, isDefaultVar=True) \
                      or chemComp.findFirstChemCompVar(isDefaultVar=True) \
                      or chemComp.findFirstChemCompVar()
        
        if chemCompVar:
          molFormula = chemCompVar.formula
        else:
          molFormula = ''  
      
        chemCompDict[ccpCode] = [chemComp.code1Letter,
                                 chemComp.code3Letter,
                                 chemComp.name,
                                 None]   # Added RHF 1/7/10 for bug fix.
      
  return  chemCompDict



####################################################################
#
# Bonds between atoms
#

def getNumConnectingBonds(atom1, atom2, limit=5):
  """
  Get the minimum number of binds that connect two atoms.
  Stops at a specified limit (and returns None if not within it)
  
  .. describe:: Input
  
  MolSystem.Atom, MolSystem.atom, Int

  .. describe:: Output
  
  Int
  """
  
  num = 0
  atoms = set([atom1,])
  
  while atom2 not in atoms:
    if num > limit:
      return None
      
    atoms2 = atoms.copy()
    
    for atom in atoms2:
      atoms.update(getBoundAtoms(atom))

    num += 1

  return num

def areAtomsTocsyLinked(atom1, atom2):
  """
  Determine if two atoms have a connectivity that may be observable in a TOCSY experiment
  
  .. describe:: Input
  
  MolSystem.Atom, MolSystem.atom

  .. describe:: Output
  
  Boolean
  """
   
  if not hasattr(atom1, 'tocsyDict'):
    atom1.tocsyDict = {}
  elif atom1.tocsyDict.has_key(atom2):
    return atom1.tocsyDict[atom2]

  if not hasattr(atom2, 'tocsyDict'):
    atom2.tocsyDict = {}
  elif atom2.tocsyDict.has_key(atom1):
    return atom2.tocsyDict[atom1]
    
  chemAtom1 = atom1.chemAtom
  chemAtom2 = atom2.chemAtom
  element1  = chemAtom1.elementSymbol
  element2  = chemAtom2.elementSymbol
  
  if element1 != element2:
    boolean = False
    
  elif areAtomsBound(atom1, atom2):
    boolean = True
  
  else:
 
    residue1 = atom1.residue
    residue2 = atom2.residue
 
    if residue1 is not residue2:
      boolean = False
      
    else:
      atomsA = set([atom1,])
      boolean = True
      while atom2 not in atomsA:
        atomsB = atomsA.copy()
 
        for atomB in atomsB:
          for atom3 in getBoundAtoms(atomB):
            if atom3.residue is not residue1:
              continue
            
            if element1 == 'H':
              if atom3.chemAtom.elementSymbol != 'H':
                for atom4 in getBoundAtoms(atom3):
                  if atom4.chemAtom.elementSymbol == 'H':
                    break
                else:
                  continue
 
            if atom3.chemAtom.elementSymbol == element1:
              if not hasattr(atom3, 'tocsyDict'):
                atom3.tocsyDict = {}
 
              atom1.tocsyDict[atom3] = True
              atom3.tocsyDict[atom1] = True
 
            atomsA.add(atom3)
 
        if atomsA == atomsB: # Nothing more to add and atom2 not in linked set
          boolean = False
          break

  atom1.tocsyDict[atom2] = boolean
  atom2.tocsyDict[atom1] = boolean
  return boolean
  
  

def getBoundAtoms(atom):
  """Get a list of atoms bound to a given atom..
  .. describe:: Input
  
  MolSystem.Atom

  .. describe:: Output
  
  List of MolSystem.Atoms
  """
  
  if hasattr(atom, 'boundAtoms'):
    return atom.boundAtoms
  
  atoms    = []
  chemAtom = atom.chemAtom
  residue  = atom.residue
  
  chemAtomDict = {}
  for atom2 in residue.atoms:
    # Only atoms specific to ChemCompVar :-)
    chemAtomDict[atom2.chemAtom] = atom2
  
  for chemBond in chemAtom.chemBonds:
    for chemAtom2 in chemBond.chemAtoms:
      if chemAtom2 is not chemAtom:
        atom2 = chemAtomDict.get(chemAtom2)
        if atom2:
          atoms.append(atom2)
  
  linkEnd = residue.chemCompVar.findFirstLinkEnd(boundChemAtom=chemAtom)
  if linkEnd:
    molResLinkEnd = residue.molResidue.findFirstMolResLinkEnd(linkEnd=linkEnd)
    
    if molResLinkEnd:
      molResLink = molResLinkEnd.molResLink
    
      if molResLink:
        for molResLinkEnd2 in molResLink.molResLinkEnds:
          if molResLinkEnd2 is not molResLinkEnd:
            residue2 = residue.chain.findFirstResidue(molResidue=molResLinkEnd2.molResidue)
            
            if residue2:
              chemAtom2 = molResLinkEnd2.linkEnd.boundChemAtom
              atom2 = residue2.findFirstAtom(chemAtom=chemAtom2)
              
              if atom2:
                atoms.append(atom2)
              
            break
  
  atom.boundAtoms = atoms    
  return atoms
  
  
def areAtomsBound(atom1, atom2):
  """Dertemine whether two atoms are bonded together
  .. describe:: Input
  
  MolSystem.Atom, MolSystem.Atom

  .. describe:: Output
  
  Boolean
  """

  if not hasattr(atom1, 'isAtomBound'):
    atom1.isAtomBound = {}
  elif atom1.isAtomBound.has_key(atom2):
    return atom1.isAtomBound[atom2]

  if not hasattr(atom2, 'isAtomBound'):
    atom2.isAtomBound = {}
  elif atom2.isAtomBound.has_key(atom1):
    return atom2.isAtomBound[atom1]

  isBound = False
  
  if atom1 is not atom2:
    residue1 = atom1.residue
    residue2 = atom2.residue
    
    if residue2.chain is residue1.chain:
      if residue2 is not residue1:

        linkEnd1 = residue1.chemCompVar.findFirstLinkEnd(boundChemAtom=atom1.chemAtom)
        if not linkEnd1:
          isBound = False

        else:
          linkEnd2 = residue2.chemCompVar.findFirstLinkEnd(boundChemAtom=atom2.chemAtom)
          if not linkEnd2:
            isBound = False
 
          else:
            molResLinkEnd1 = residue1.molResidue.findFirstMolResLinkEnd(linkEnd=linkEnd1)
            if not molResLinkEnd1:
              isBound = False
 
            else:
              molResLinkEnd2 = residue2.molResidue.findFirstMolResLinkEnd(linkEnd=linkEnd2)
              if not molResLinkEnd2:
                isBound = False
 
              elif molResLinkEnd2 in molResLinkEnd1.molResLink.molResLinkEnds:
                isBound = True
 
              else:
                isBound = False

      else:
        for chemBond in atom1.chemAtom.chemBonds:
          if atom2.chemAtom in chemBond.chemAtoms:
            isBound = True
            break

  atom1.isAtomBound[atom2] = isBound
  atom2.isAtomBound[atom1] = isBound
  
  return isBound


def areResonancesBound(resonance1,resonance2):
  """
  Determine whether two resonances are assigned to directly bonded atoms
  
  .. describe:: Input
  
  Nmr.Resonance, Nmr.Resonance

  .. describe:: Output
  
  Boolean
  """

  if resonance1 is resonance2:
    return False
  
  
  resonanceSet1 = resonance1.resonanceSet
  resonanceSet2 = resonance2.resonanceSet
    
  if resonanceSet1 and resonanceSet2:
    
    
    atomSets1 = resonanceSet1.atomSets
    atomSets2 = resonanceSet2.atomSets
    
    bound1 = resonance1.covalentlyBound
    bound2 = resonance2.covalentlyBound
    
    # Have to look through everything to get the right equiv
    # & prochiral pair - Val CGa HGb: check both atomSets for each
    # Phe Ce* He*: check both atoms (only one atomSet each)
    for atomSet1 in atomSets1:
      for atom1 in atomSet1.atoms:
        for atomSet2 in atomSets2:
          for atom2 in atomSet2.atoms:
            if areAtomsBound(atom1, atom2):
              # Val Cgb - Hgb* can appear bound,
              # so check resonance links
              if (resonance1.isotopeCode == '1H') and ( len(atomSets2) > 1 ):
                if bound1 and resonance2 not in bound1:
                  continue

              elif (resonance2.isotopeCode == '1H') and ( len(atomSets1) > 1 ):
                if bound2 and resonance1 not in bound2:
                  continue
              
              return True

    return False

  from ccpnmr.analysis.core.AssignmentBasic import getBoundResonances
 
  resonances = getBoundResonances(resonance1)
  if resonance2 in resonances:
    return True
  
  else:
    return False  
    
    
    
def getAtomsTorsion(atoms):
  """Get the chemical torsion object equivalent
             to an ordered input list of four atoms
  .. describe:: Input
  
  4-List of MolSystem.Atoms

  .. describe:: Output
  
  ChemComp.ChemTorsion
  """

  assert len(atoms) == 4
  
  residue = list(atoms)[0].residue
  for i in (1,2,3):
    if not areAtomsBound(atoms[i-1],atoms[i]):
      print ('ERROR, Get atom torsion failed: Atoms %s and %s are not bound' % (atoms[i-1],atoms[i]))
      #showWarning('Failure',
      #            'Get atom torsion failed: Atoms %s and %s are not bound' % (atoms[i-1],atoms[i]))
      return
    
    residue1 = atoms[i].residue
    if residue1 is not residue:
      if residue1.seqCode > residue.seqCode:
        residue = residue1
    
  chemCompVar = residue.chemCompVar
  chemAtoms   = []
  linkAtoms   = []
  
  for atom in atoms:
    if atom.residue is residue:
      chemAtoms.append(atom.chemAtom)
    else:
      if linkAtoms:
        linkAtoms.append(chemCompVar.findFirstChemAtom(name='prev_2'))
      else:
        linkAtoms.append(chemCompVar.findFirstChemAtom(name='prev_1'))
  
  chemAtoms = chemAtoms + linkAtoms
  chemAtoms.sort()

  for chemTorsion in chemCompVar.chemTorsions:
    chemAtoms0 = list(chemTorsion.chemAtoms)
    chemAtoms0.sort()
    
    if chemAtoms0 == chemAtoms:
      return chemTorsion



####################################################################
#
# Residues
#


def getResidueCode(obj):
  """Get a text code for a residue/molResidue/resonanceGroup,
             defaults to the ccpCode if a custom code is not present
             in userResidueCodesDict 
  .. describe:: Input
  
  MolSystem.Residue

  .. describe:: Output
  
  Word
  """ 
  
  ccpCode = obj.ccpCode
  molType = obj.molType
    
  if molType is None:
    molType = PROTEIN_MOLTYPE
  
  ccpCodeDict = userResidueCodesDict.get(molType, {})
  residueCode = ccpCodeDict.get(ccpCode, ccpCode)
  
  if molType == CARBOHYDRATE_MOLTYPE:
    descriptor = obj.descriptor
    n = len(STEREO_PREFIX)
    
    if descriptor and (len(descriptor)> n) and descriptor.startswith(STEREO_PREFIX):
      if descriptor[n] == '1':
        residueCode = 'a-' + residueCode
      elif descriptor[n] == '2':
        residueCode = 'b-' + residueCode
 
  return residueCode
  
  
def makeResidueLocTag(residue, chemAtoms):
  """Make unique identifier for a given residue type in a given chain location
  .. describe:: Input
  
  MolSystem.Residue, List of ChemCmp.chemAtoms

  .. describe:: Output
  
  Word
  """

  # check if any atoms have been deleted (to make a molSystemLink)
  tag = residue.ccpCode + residue.linking
  
  if len(chemAtoms) != len(residue.atoms):
    names = []
    for chemAtom in chemAtoms:
      names.append(chemAtom.name)
    for atom in residue.atoms:
      if atom.name in names:
        names.remove(atom.name)
    for name in names:
      tag = tag + '-' + name
      
  return tag
  

def getLinkedResidue(residue, linkCode='prev'):
  """Find a residue, if it exists, that is linked to the
             input residue by a given type of molResidue link.
  .. describe:: Input
  
  MolSystem.Residue

  .. describe:: Output
  
  MolSystem.Residue
  """  

  if not hasattr(residue, 'linkedResidueDict'):
    residue.linkedResidueDict = {}
  else:
    if residue.linkedResidueDict.has_key(linkCode):
      return residue.linkedResidueDict[linkCode]
  
  residue2    = None
  chain       = residue.chain
  molResidue  = residue.molResidue
  linkEnd     = molResidue.findFirstMolResLinkEnd(linkCode=linkCode)
  
  if linkEnd:
    molResLink = linkEnd.molResLink
    if molResLink:
      for linkEnd2 in molResLink.molResLinkEnds:
        if linkEnd2 is not linkEnd:
          residue2 = chain.findFirstResidue(molResidue=linkEnd2.molResidue)
  
  residue.linkedResidueDict[linkCode] = residue2
  return residue2


def getResidueObservableAtoms(residue, refExperiment=None, labelling=None,
                              minFraction=0.1, jCouplingBonds=(1,2,3),
                              usePermissiveShifts=False,
                              chemElements=('H','C','N','F','P')): 
  """
  Determine which atoms of a chem comp varient would give rise to
  observable resonances considering a given reference experiment
  and/or an isotope labelling scheme. Can specify minimum fraction of
  an isotope to consider something observable and the chemical elements which
  you are observing. Boolean option to match database min and max
  chemical shift bounds to atom sites, rather than randon coil shift
  values (default).
  
  .. describe:: Input
  
  MolSystem.Residue, NmrExpPrototype.RefExperiment,
  ChemCompLabel.LabelingScheme or LabeledMolecule.LabeledMixture,
  Float, Boolean, List of Words

  .. describe:: Output
  
  List of ChemComp.ChemAtoms
  """

  if not jCouplingBonds:
    jCouplingBonds = (0,)

  atomSiteDict   = {}
  isotopomerDict = {}
  atomSitesAll   = {}
  
  if refExperiment:
    for atomSite in refExperiment.nmrExpPrototype.atomSites:
      isotope = atomSite.isotopeCode
      if not atomSitesAll.has_key(isotope):
        atomSitesAll[isotope] = []
  
      atomSitesAll[isotope].append(atomSite)
 
  isotopeDict = {}
  
  if atomSitesAll:
    isotopes = atomSitesAll.keys()
  
  else:
    isotopes = []
    for element in chemElements:
      isotope = DEFAULT_ISOTOPES.get(element)
    
      if isotope:
        isotopes.append(isotope)
  
  for isotope in isotopes:
    element = isotope
 
    while element[0] in '0123456789':
      element = element[1:]
 
    isotopeDict[element] = isotope

 
  filteredAtoms = []
  prevResidue = getLinkedResidue(residue, linkCode='prev')
  nextResidue = getLinkedResidue(residue, linkCode='next')
  
  natAbundance = residue.root.findFirstLabelingScheme(name='NatAbun')
  
  #print residue.seqCode, residue.ccpCode, refExperiment.name
  for residue0 in (prevResidue,residue,nextResidue):
    isotopomers = None
    
    if residue0:
      resId = residue0.molResidue.serial
      atoms = residue0.atoms
      
      # Compile isotopomers for this residue
      if labelling and (labelling.className == 'LabelingScheme'):
        chemComp   = residue0.chemCompVar.chemComp
        ccpCode    = chemComp.ccpCode
        molType    = chemComp.molType
        chemCompLabel = labelling.findFirstChemCompLabel(ccpCode=ccpCode,
                                                         molType=molType)

        if not chemCompLabel:
          chemCompLabel = natAbundance.findFirstChemCompLabel(ccpCode=ccpCode,
                                                              molType=molType)


        if chemCompLabel:
          isotopomers  = chemCompLabel.isotopomers
          isotopomerDict[residue0] = isotopomers
            
          
    else:
      atoms = []  
   
    #atoms0 = [] # Those which make it through the filter
    for atom in atoms:
      chemAtom = atom.chemAtom
      isotope  = isotopeDict.get(chemAtom.elementSymbol)
    
      if not isotope:
        continue
    
      if chemAtom.waterExchangeable:
        continue

      if isotopomers:
        fractionDict = getIsotopomerSingleAtomFractions(isotopomers,atom.name,chemAtom.subType)
        # Exclude if no isotope incorporation above threshold
        fraction = fractionDict.get(isotope, minFraction)
        if fraction < minFraction:
          continue
          
      elif labelling:
        fractionDict = singleAtomFractions(labelling, resId, atom.name)
        if not fractionDict:
          continue
        
        fraction = fractionDict.get(isotope, minFraction)
        if fraction < minFraction:
          continue
        
 
      atomSitesIsotope = atomSitesAll.get(isotope)
      if atomSitesIsotope:
        setSize = None
        
        if usePermissiveShifts:
          shifts = getChemicalShiftBounds(chemAtom)
          if not shifts:
            shifts = [getRandomCoilShift(chemAtom),]

        else:
          shifts = [getRandomCoilShift(chemAtom),]
 
        for atomSite in atomSitesIsotope:
                
          maxShift = atomSite.maxShift
          if (maxShift is not None) and (shifts[0] > maxShift):
            continue
 
          minShift = atomSite.minShift
          if (minShift is not None) and (shifts[-1] < minShift):
            continue
        
          if setSize is None:
            setSize     = 1
            chemAtomSet = chemAtom.chemAtomSet
 
            if chemAtomSet:
              setSize = len(chemAtomSet.chemAtoms)
        
          minNumber =  atomSite.minNumber
          if setSize < minNumber:
            continue
 
          maxNumber = atomSite.maxNumber
          if maxNumber and (setSize>maxNumber):
            continue
 
          numberStep = atomSite.numberStep
          if (setSize-minNumber) % numberStep != 0:
            continue
             
          if atomSiteDict.get(atomSite) is None:
            atomSiteDict[atomSite] = []
          atomSiteDict[atomSite].append(atom)
          
          #print 'AS', atomSite.name, atom.name
        
      filteredAtoms.append(atom)
      
  
  if refExperiment:
    #print refExperiment.name
    
    # Atom sites which are possibly visible given dims
    observableAtomSites = {}
    for refExpDim in refExperiment.refExpDims:
      for refExpDimRef in refExpDim.refExpDimRefs:
        for atomSite in refExpDimRef.expMeasurement.atomSites:
          observableAtomSites[atomSite] = True
  
    # Get prototype graph atomSite routes
  
    graphRoutes = []
    for expGraph in refExperiment.nmrExpPrototype.expGraphs:
      expSteps = [(es.stepNumber, es) for es in expGraph.expSteps]
      expSteps.sort()
      routes = []
      stepNumber, expStep = expSteps[0]

      for atomSite in expStep.expMeasurement.atomSites:
        route = [(atomSite,None,stepNumber)]
        routes.append(route)
      
      while True:
        routes2 = []
      
        for route in routes:
          atomSiteA, null, stepA = route[-1]
          #print atomSiteA.name, step
      
          for expTransfer in atomSiteA.expTransfers:
            atomSites = list(expTransfer.atomSites)
            atomSites.remove(atomSiteA)
            atomSiteB = atomSites[0]
 
            if not expTransfer.transferToSelf:
              if atomSiteB is atomSiteA:
                continue
 
            for stepB, expStepB in expSteps:
              if stepA > stepB:
                continue
            
              if atomSiteB in expStepB.expMeasurement.atomSites:
                routes2.append( route[:] + [(atomSiteB,expTransfer,stepB)] )
                #print ['%s %d' % (a[0].name, a[2]) for a in routes2[-1]]
                break
        
        if routes2:
          routes = routes2
        else:
          break  
      
      for route in routes:
        atomRoutes = []
        lastAtomSite = route[-1][0]
      
        for i in range(len(route)-1):
          atomSiteA, null, stepA = route[i]
          atomSiteB, expTransfer, stepB = route[i+1] 
          transferType = expTransfer.transferType
          
          #print stepA, atomSiteA.name, stepB, atomSiteB.name, transferType

          if atomRoutes:
            atomsA = [r[-1][0] for r in atomRoutes]
          else:
            atomsA = atomSiteDict[atomSiteA]
            
          atomRoutes2 = []
          for atomA in atomsA:
            for atomB in atomSiteDict[atomSiteB]:
              if isotopomerDict:
                chemAtomA = atomA.chemAtom
                chemAtomB = atomB.chemAtom
                subTypeA  = chemAtomA.subType
                subTypeB  = chemAtomB.subType
                isotopeA  = isotopeDict[chemAtomA.elementSymbol]
                isotopeB  = isotopeDict[chemAtomB.elementSymbol]
                residueA  = atomA.residue
                residueB  = atomB.residue
 
                if residueA is residueB:
                  isotopomersA = isotopomerDict.get(residueA)
                  atomNames    = (atomA.name, atomB.name)
                  subTypes     = (subTypeA, subTypeB)
                  pairDict     = getIsotopomerAtomPairFractions(isotopomersA, atomNames, subTypes)
                  fraction     = pairDict.get((isotopeA, isotopeB), minFraction)
 
                  if fraction  < minFraction:
                    continue
 
                else:
                  isotopomersA = isotopomerDict.get(residueA)
                  isotopomersB = isotopomerDict.get(residueB)

                  if isotopomersB and isotopomersA:
                    fractionDictA = getIsotopomerSingleAtomFractions(isotopomersA, atomA.name, subTypeA)
                    fractionDictB = getIsotopomerSingleAtomFractions(isotopomersB, atomB.name, subTypeB)
                    fraction = fractionDictA.get(isotopeA, 1.0) * fractionDictB.get(isotopeB, 1.0)
 
                    if fraction < minFraction:
                      continue
                      
              elif labelling:
                chemAtomA = atomA.chemAtom
                chemAtomB = atomB.chemAtom
                isotopeA  = isotopeDict[chemAtomA.elementSymbol]
                isotopeB  = isotopeDict[chemAtomB.elementSymbol]
                residueA = atomA.residue
                residueB = atomB.residue
                molResidueA = residueA.molResidue
                molResidueB = residueB.molResidue
                resIds = (molResidueA.serial, molResidueB.serial)
                atomNames = (atomA.name, atomB.name)
               
                pairDict = atomPairFractions(labelling, resIds, atomNames)
                fraction = pairDict.get((isotopeA, isotopeB), minFraction)
 
                if fraction  < minFraction:
                  continue
                      
              addAtom = False
              if transferType in longRangeTransfers:
                addAtom = True
 
              elif transferType in ('onebond','CP') and areAtomsBound(atomA, atomB):
                addAtom = True
 
              elif transferType == 'TOCSY'and areAtomsTocsyLinked(atomA, atomB):
                addAtom = True
 
              elif transferType == 'Jcoupling':
                numBonds = getNumConnectingBonds(atomA, atomB, limit=max(jCouplingBonds))
                if numBonds in jCouplingBonds:
                  addAtom = True
 
              elif transferType == 'Jmultibond' and not areAtomsBound(atomA, atomB):
                numBonds = getNumConnectingBonds(atomA, atomB, limit=max(jCouplingBonds))
                if numBonds in jCouplingBonds:
                  addAtom = True
 
              if addAtom:
                grown = True
                #print 'AB', atomA.name, atomA.residue.seqCode,'+', atomB.name, atomB.residue.seqCode
                if not atomRoutes:
                  atomRoutes2.append( [(atomA,atomSiteA),(atomB,atomSiteB),] )
                  #print atomA.name, atomB.name
                
                else:
                  for atomRoute in atomRoutes:
                    atomRoutes2.append( atomRoute[:] + [(atomB,atomSiteB),] )
                  #print '+', atomB.name

               
          atomRoutes = []
          for atomRoute in atomRoutes2:
            if atomRoute[-1][1] is lastAtomSite:
              atomRoutes.append(atomRoute)
              
        graphRoutes.append(atomRoutes)

 
    observableAtoms = set()
    for routes in graphRoutes:
      for route in routes:

        for atomB, atomSiteB in route:
          if atomB.residue is residue: # Must have one atom from this residue
            for atomA, atomSiteA in route:
              if observableAtomSites.get(atomSiteA):
                observableAtoms.add(atomA)
            break
   
  else:
    observableAtoms = filteredAtoms
   
  return list(observableAtoms)
  
  
  

####################################################################
#
# Shift - related
#
# NBNB TBD consider if they need to be moved to ChemicalSHiftBasic
#



def getChemicalShiftBounds(chemAtom, threshold=0.001, sourceName='BMRB'):
  """
  Return the min and max chemical shifts for a
  given atom type observed in the databases
  
  .. describe:: Input
  
  ChemComp.ChemAtom, Float, Word

  .. describe:: Output
  
  Float
  """

  key = '%s:%s' % (sourceName, threshold)
  if not hasattr(chemAtom, 'chemicalShiftBounds'):
    chemAtom.chemicalShiftBounds = {}

  else:
    region = chemAtom.chemicalShiftBounds.get(key)
    if region:
      return region
  
  chemComp = chemAtom.chemComp
  ccpCode  = chemComp.ccpCode
  molType  = chemComp.molType
  
  chemAtomNmrRef = getChemAtomNmrRef(chemAtom.root, chemAtom.name, ccpCode,
                                     molType=molType, sourceName=sourceName)

  if chemAtomNmrRef:
    distribution  = chemAtomNmrRef.distribution
    refPoint      = chemAtomNmrRef.refPoint
    refValue      = chemAtomNmrRef.refValue
    valuePerPoint = chemAtomNmrRef.valuePerPoint
    n = len(distribution)
 
    minPt  = 0
    maxPt  = n-1
 
    for v in distribution:
      if v < threshold:
        minPt += 1
      else:
        break
 
    for v in distribution[::-1]:
      if v < threshold:
        maxPt -= 1
      else:
        break
 
    maxPpm = ((maxPt - refPoint) * valuePerPoint) + refValue
    minPpm = ((minPt - refPoint) * valuePerPoint) + refValue
    region  = (minPpm, maxPpm)
    chemAtom.chemicalShiftBounds[key] = region

  return region
       
       
def getRandomCoilShift(chemAtom, context=None, sourceName='BMRB'):
  """
  Get the random coil chemical shift value of a chemAtom
 
  .. describe:: Input
  
  ChemComp.ChemAtom

  .. describe:: Output
  
  Float
  """
  
  value = None
  if not hasattr(chemAtom, 'randomCoilShiftDict'):
    chemAtom.randomCoilShiftDict = {}
    
  elif chemAtom.randomCoilShiftDict.has_key(sourceName):
    value = chemAtom.randomCoilShiftDict[sourceName]
 
    if not context:
      return value
 
  if value is None:
    chemComp = chemAtom.chemComp
    ccpCode  = chemComp.ccpCode
    molType  = chemComp.molType
    chemAtomNmrRef = getChemAtomNmrRef(chemAtom.root, chemAtom.name, ccpCode,
                                       molType=molType, sourceName=sourceName)
    
    if not chemAtomNmrRef and chemAtom.chemAtomSet:
      chemAtomNmrRef = getChemAtomNmrRef(chemAtom.root, chemAtom.chemAtomSet.name, ccpCode,
                                         molType=molType, sourceName=sourceName)

    if chemAtomNmrRef:
      value = chemAtomNmrRef.randomCoilValue
      if value is None:
        value = chemAtomNmrRef.meanValue
 
    chemAtom.randomCoilShiftDict[sourceName] = value

  if context and (value is not None):
    correctionDict = getRandomCoilShiftCorrectionsDict()
    
    atomName = chemAtom.name
    
    if atomName in ('HA2','HA3'):
      atomName = 'HA'
    
    offsets = [2,1,-1,-2]
    indices = [0,1,3,4]
    
    for i in range(4):
      residue = context[indices[i]]
      
      if residue is not None:
        ccpCode = residue.ccpCode
        atomDict = correctionDict[offsets[i]].get(ccpCode, {})
        value += atomDict.get(atomName, 0.0)
    
  return value
  
  
  
  
def getRandomCoilShiftCorrectionsDict():
  """

  Citation

  Schwarzinger, S., Kroon, G. J. A., Foss, T. R., Chung, J., Wright, P. E., Dyson, H. J.
  "Sequence-Dependent Correlation of Random Coil NMR Chemical Shifts", 
  J. Am. Chem. Soc. 123, 2970-2978 (2001)

  Values obtained from a GGXGG sequence pentapeptide. 

  """

  data = """
  Ala    H     H   -0.01  -0.05   0.07  -0.10
  Ala    HA    H   -0.02  -0.03  -0.03   0.00
  Ala    C     C   -0.11  -0.77  -0.07  -0.02
  Ala    CA    C   -0.02  -0.17   0.06   0.01
  Ala    N     N   -0.12  -0.33  -0.57  -0.15
  Asn    H     H   -0.01  -0.03   0.13  -0.07
  Asn    HA    H   -0.01  -0.01  -0.02  -0.01
  Asn    C     C   -0.09  -0.66  -0.10  -0.03
  Asn    CA    C   -0.06  -0.03   0.23   0.01
  Asn    N     N   -0.18  -0.26   0.87  -0.17
  Asp    H     H   -0.02  -0.03   0.14  -0.11
  Asp    HA    H   -0.02  -0.01  -0.02  -0.01
  Asp    C     C   -0.08  -0.58  -0.13  -0.04
  Asp    CA    C   -0.03   0.00   0.25  -0.01
  Asp    N     N   -0.12  -0.20   0.86  -0.29
  Arg    H     H    0.00  -0.02   0.15  -0.06
  Arg    HA    H   -0.02  -0.02  -0.02   0.00
  Arg    C     C   -0.06  -0.49  -0.19  -0.03
  Arg    CA    C    0.00  -0.07  -0.01   0.02
  Arg    N     N   -0.06  -0.14   1.62  -0.06
  Cys    H     H    0.00  -0.02   0.20  -0.07
  Cys    HA    H   -0.01   0.02   0.00   0.00
  Cys    C     C   -0.08  -0.51  -0.28  -0.07
  Cys    CA    C   -0.03  -0.07   0.10  -0.01
  Cys    N     N   -0.06  -0.26   3.07   0.00
  Gln    H     H   -0.01  -0.02   0.15  -0.06
  Gln    HA    H   -0.01  -0.02  -0.01   0.00
  Gln    C     C   -0.05  -0.48  -0.18  -0.03
  Gln    CA    C   -0.02  -0.06   0.04   0.01
  Gln    N     N   -0.06  -0.14   1.62  -0.06
  Glu    H     H   -0.01  -0.03   0.15  -0.07
  Glu    HA    H   -0.02  -0.02  -0.02   0.00
  Glu    C     C   -0.09  -0.48  -0.20  -0.03
  Glu    CA    C   -0.01  -0.08   0.05   0.01
  Glu    N     N   -0.06  -0.20   1.51  -0.12
  Gly    H     H    0.00   0.00   0.00   0.00
  Gly    HA    H    0.00   0.00   0.00   0.00
  Gly    C     C    0.00   0.00   0.00   0.00
  Gly    CA    C    0.00   0.00   0.00   0.00
  Gly    N     N    0.00   0.00   0.00   0.00
  His    H     H   -0.01  -0.04   0.20   0.00
  His    HA    H   -0.03  -0.06   0.01   0.01
  His    C     C   -0.10  -0.65  -0.22  -0.07
  His    CA    C   -0.05  -0.09   0.02   0.01
  His    N     N   -0.12  -0.55   1.68   0.17
  Ile    H     H   -0.01  -0.06   0.17  -0.09
  Ile    HA    H   -0.03  -0.02  -0.02  -0.01
  Ile    C     C   -0.20  -0.58  -0.18  -0.02
  Ile    CA    C   -0.07  -0.20  -0.01   0.02
  Ile    N     N   -0.18  -0.14   4.87   0.00
  Leu    H     H    0.00  -0.03   0.14  -0.08
  Leu    HA    H   -0.04  -0.03  -0.05  -0.01
  Leu    C     C   -0.13  -0.50  -0.13  -0.01
  Leu    CA    C   -0.01  -0.10   0.03   0.02
  Leu    N     N   -0.06  -0.14   1.05  -0.06
  Lys    H     H    0.00  -0.03   0.14  -0.06
  Lys    HA    H   -0.02  -0.02  -0.01   0.00
  Lys    C     C   -0.08  -0.50  -0.18  -0.03
  Lys    CA    C   -0.01  -0.11  -0.02   0.02
  Lys    N     N   -0.06  -0.20   1.57  -0.06
  Met    H     H    0.00  -0.02   0.15  -0.06
  Met    HA    H   -0.02  -0.01  -0.01   0.00
  Met    C     C   -0.08  -0.41  -0.18  -0.02
  Met    CA    C    0.00   0.10  -0.06   0.01
  Met    N     N   -0.06  -0.20   1.57  -0.06
  Phe    H     H   -0.03  -0.12   0.10  -0.37
  Phe    HA    H   -0.06  -0.09  -0.08  -0.04
  Phe    C     C   -0.27  -0.83  -0.25  -0.10
  Phe    CA    C   -0.07  -0.23   0.06   0.01
  Phe    N     N   -0.18  -0.49   2.78  -0.46
  Pro    H     H   -0.04  -0.18   0.19  -0.12
  Pro    HA    H   -0.01   0.11  -0.03  -0.01
  Pro    C     C   -0.47  -2.84  -0.09  -0.02
  Pro    CA    C   -0.22  -2.00   0.02   0.04
  Pro    N     N   -0.18  -0.32   0.87  -0.17
  Ser    H     H    0.00  -0.03   0.16  -0.08
  Ser    HA    H   -0.01   0.02   0.00  -0.01
  Ser    C     C   -0.08  -0.40  -0.15  -0.06
  Ser    CA    C    0.00  -0.08   0.13   0.00
  Ser    N     N   -0.06  -0.03   2.55  -0.17
  Thr    H     H    0.01   0.00   0.14  -0.06
  Thr    HA    H   -0.01   0.05   0.00  -0.01
  Thr    C     C   -0.08  -0.19  -0.13  -0.05
  Thr    CA    C   -0.01  -0.04   0.12   0.00
  Thr    N     N   -0.06  -0.03   2.78  -0.12
  Trp    H     H   -0.08  -0.13   0.04  -0.62
  Trp    HA    H   -0.08  -0.10  -0.15  -0.16
  Trp    C     C   -0.26  -0.85  -0.30  -0.17
  Trp    CA    C   -0.02  -0.17   0.03  -0.08
  Trp    N     N    0.00  -0.26   3.19  -0.64
  Tyr    H     H   -0.04  -0.11   0.09  -0.42
  Tyr    HA    H   -0.05  -0.10  -0.08  -0.04
  Tyr    C     C   -0.28  -0.85  -0.24  -0.13
  Tyr    CA    C   -0.07  -0.22   0.06  -0.01
  Tyr    N     N   -0.24  -0.43   3.01  -0.52
  Val    H     H   -0.01  -0.05   0.17  -0.08
  Val    HA    H   -0.02  -0.01  -0.02  -0.01
  Val    C     C   -0.20  -0.57  -0.18  -0.03
  Val    CA    C   -0.07  -0.21  -0.02   0.01
  Val    N     N   -0.24  -0.14   4.34  -0.06
  """
  
  rcsDict = {}
  offsets = [-2,-1,1,2]
  for o in offsets:
    rcsDict[o] = {}
  
  lines = data.split('\n')
  for line in lines:
    array = line.split()
    if array:
      tlc    = array[0]
      atom   = array[1]
      values = [float(v) for v in array[3:]]
      
      for i in range(4):
        offset = offsets[i]
        if rcsDict[offset].get(tlc) is None:
          rcsDict[offset][tlc] = {}
  
        rcsDict[offset][tlc][atom] = values[i]

  return rcsDict
  
  


def getChemAtomNmrRef(project, atomName, ccpCode, molType=PROTEIN_MOLTYPE, sourceName='BMRB'):
  """
  Retrieve an NMR chemical shift reference atom record
  
  .. describe:: Input
  
  Implementation.Project, Word (ChemAtom.name), Word (ChemComp.ccpCode),
  Word, (chemComp.molType), Word

  .. describe:: Output
  
  Float
  """
  
  #NBNB TBD Merge with similar function in CheicalShiftBasic???
  
  #key = '%s:%s:%s:%s' % (atomName, ccpCode, molType, sourceName)
  #if not hasattr(project, 'chemAtomNmrRefDict'):
  #  project.chemAtomNmrRefDict = {}
  #else:
    

  nmrRefStore = project.findFirstNmrReferenceStore(molType=molType, ccpCode=ccpCode)
  
  chemAtomNmrRef = None
  if nmrRefStore:
    chemCompNmrRef = nmrRefStore.findFirstChemCompNmrRef(sourceName=sourceName)
    
    if chemCompNmrRef:
      chemCompVarNmrRef = chemCompNmrRef.findFirstChemCompVarNmrRef(linking='any',
                                                                    descriptor='any')
      
      if chemCompVarNmrRef:
        for chemAtomNmrRef1 in chemCompVarNmrRef.chemAtomNmrRefs:
          if atomName == chemAtomNmrRef1.name:
            chemAtomNmrRef = chemAtomNmrRef1
            
            break
      else:
        data = (molType,ccpCode)
        msg  = 'Could not load reference NMR data for'
        msg += 'general chem comp variant %s:%s' % data
        showWarning('Warning', msg)
        return
        
    else:
      data = (molType,ccpCode,sourceName)
      showWarning('Warning',
                  'Could not load reference NMR data for %s:%s source=%s' % data)
      return
      
  else:
    data = (molType,ccpCode)
    showWarning('Warning',
                'Could not load reference NMR data for %s:%s' % data)
    return
 
  if not chemAtomNmrRef:
    atomName2 = atomName[:-1]
    for chemAtomNmrRef1 in chemCompVarNmrRef.chemAtomNmrRefs:
      if atomName2 == chemAtomNmrRef1.name:
        chemAtomNmrRef = chemAtomNmrRef1
        break

  
  return chemAtomNmrRef

##########################################################################
#
# Various
#


def makeGuiName(name, elementSymbol):
  """Convert atom or atomSet name into name for gui
  .. describe:: Input
  
  Word (Nmr.AtomSet.name), Word

  .. describe:: Output
  
  Word 
  """
  if not name.upper().startswith(elementSymbol.upper()):
    raise Exception("Atom name %s does not start with element symbol %s"
                    % (name, elementSymbol))
   
  return elementSymbol + name[len(elementSymbol):].lower()  



def _getUnicodeGreek():
  """
   {'a':u'\u03B1','b':u'\u03B2','g':u'\u03B3','d':u'\u03B4','e':u'\u03B5',
    'z':u'\u03B6','h':u'\u03B7','q':u'\u03B8','i':u'\u03B9','k':u'\u03BA',
    'l':u'\u03BB','m':u'\u03BC','n':u'\u03BD','x':u'\u03BE','o':u'\u03BF',
    'p':u'\u03C0','r':u'\u03C1','j':u'\u03C2','s':u'\u03C3','t':u'\u03C4', # j : Other sigma
    'u':u'\u03C5','f':u'\u03C6','c':u'\u03C7','y':u'\u03C8','w':u'\u03C9'}
  """
  dict = {}
  romanLetterOrder = 'ABGDEZHQIKLMNXOPRJSTUFCYW'
  u = 913
  l = 945
  for a in romanLetterOrder:
    dict[a] = unichr(u)
    dict[a.lower()] = unichr(l)
    u += 1
    l += 1

  return dict   
   
unicodeGreek = _getUnicodeGreek()



def getUnicodeAtomName(name, elementSymbol):
  """Return a unicode string for a protein atom name, i.e.
             including any greek characters 
  .. describe:: Input
  
  Word, Word (ChemElement.symbol)

  .. describe:: Output
  
  Unicode Word
  """ 
  
  l = len(elementSymbol)
  
  if l == len(name):
    return name
  
  else:
    letter = name[l]
    
    if letter.islower():
      return name[:l] + unicodeGreek.get(letter,letter) + name[l+1:]
    else:
      return name
      
      

def greekSortAtomNames(dataList, molType=PROTEIN_MOLTYPE):
  """Sorts a list of atom names according to greek/sidechain order when letters are latinised
  .. describe:: Input
  
  List of Strings

  .. describe:: Output
  
  List of Strings (sorted)
  """
  
  sub = re.sub
  
  #order = 'ABGDEZHQIKLMNXOPRSTUFCYW'
  
  sortList = []
  for x in dataList:
    if type(x) in (tuple,list):
      sortName = x[0]
    else:
      sortName = x
    
    sortName = sortName.upper()  
    sortName = sub('(.+\')', 'zzz@\\1', sortName)
    sortName = sub('^(\d)','zz@\\1', sortName)
    sortName = sub('N(\S*)','\'\'@N\\1',sortName)
    sortName = sub('CO','\'\'@CO',sortName)
    sortName = sub('G', 'c', sortName)
    sortName = sub('Z', 'f', sortName)
    sortName = sub('Q', 'hzz', sortName)
    sortName = sub('X', 'nzz', sortName)
    sortName = sub('F', 'v', sortName)
    sortName = sub('C', 'w', sortName)
    sortName = sub('W', 'z', sortName)
    sortName = sortName.upper()  
    
    if molType == PROTEIN_MOLTYPE:
      sortName = sub('Hn', 'H1n', sortName)
      
    sortList.append( (sortName,x) )
    
  sortList.sort()
  dataList = [e[1] for e in sortList]
  
  return dataList

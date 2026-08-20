LICENSE = """
======================COPYRIGHT/LICENSE START==========================

DataConvertLib.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2013 RAsmyus Fogh, Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk, rhf22@cam.ac.uk
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

from ccp.general.ChemCompOverview import chemCompOverview
from ccp.general.Io import getChemComp, getStdChemComps

# priority order of naming systems
namingSystemPriorityOrder = ['PDB','IUPAC','PDB_REMED','BMRB','XPLOR',
                             'CYANA2.1','DIANA', 'GROMOS','MSD','SYBYL',
                             'UCSF','AQUA','DISGEO','DISMAN', 'MOLMOL','MSI']

    
###################################################################
#
### Input data conversion
#
####################################################################


def getBestNamingSystem(residues, atomNamesList):

  """
  Descrn: Determine best naming system from residues and matching atom names
  Inputs: residues: list of residues 
          atomNamesList: list of list of atom names, matching residues 1:1

  Output: name of best namingSystem
  """
  
  # Make dict of chemComp:merged atom name sets
  atomsDict = {}
  for ii, residue in enumerate(residues):
    if residue:
      # NB residue could be None in some contexts
      names = atomNamesList[ii]
      chemComp = residue.molResidue.chemComp
      if chemComp in atomsDict:
        atomsDict[chemComp].update(names)
      else:
        atomsDict[chemComp]= set(names)
  
  # Make parallel lists of chemcomps and their atomnames
  chemComps,atomNames = zip(*atomsDict.items())
  
  # get naming system
  return getBestNamingSystemCC(chemComps, atomNames)


def getBestNamingSystemCC(chemComps, atomNamesList, prefNamingSystemName=None):

  """
  Descrn: Determine best naming system from ChemComps and matching atom names
          Does not check for duplicate ChemComps.
  Inputs: chemComps: list of ChemComps
          atomNamesList: list of collections of atom names, matching chemComps 1:1

  Output: name of best namingSystem
  """
  
  # Weights for a hit in sysname:
  fullHit = 10
  # Weights for a hit in altSysNames:
  altHit = 5
  
  # Calculate namingSystem score from atomName hits
  scoreDict = {}
  for ii,chemComp in enumerate(chemComps):
    atomNames = atomNamesList[ii]
    
    for namingSystem in chemComp.namingSystems:
      ns = namingSystem.name
      score = scoreDict.get(ns, 0)
 
      for atomName in atomNames:
      
        namingSystemR = namingSystem
        done = False
        
        while namingSystemR:
          # First choice plain mapping
          if namingSystemR.findFirstAtomSysName(sysName=atomName):
            score += fullHit
            break

          # Otherwise try alternative sys names
          else:
            for atomSysName0 in namingSystemR.atomSysNames:
              if atomName in atomSysName0.altSysNames:
                score += altHit
                done = True
                break
            
            if done:
              break
              
            else:
              #  no luck, try again with reference naming system
              namingSystemR = namingSystemR.atomReference

      scoreDict[ns] = score

  bestSc = -1
  bestNs = None
  # read naming systems in priority order (in case of ties)
  sortOrder = list(namingSystemPriorityOrder)
  if prefNamingSystemName:
    if prefNamingSystemName in sortOrder:
      sortOrder.remove(prefNamingSystemName)
    sortOrder.insert(0,prefNamingSystemName)
  for ns in sortOrder:
    if ns in scoreDict:
      score = scoreDict.pop(ns)
      if score > bestSc:
        bestSc = score
        bestNs = ns
  
  # read remaining systems, if any (in alphabetical order, for reproducability)
  for ns in sorted(scoreDict):
    score = scoreDict.pop(ns)
    if score > bestSc:
      bestSc = score
      bestNs = ns
  
  #
  return bestNs
  

def getBestChemComp(project, resName, atomNames, molType=None, download=True):
  """Find the best matching ChemComp for the resName,  using atomNames for determining molType 
  .. describe:: Input
  
  Implementation.Project, Word (imported residue name),
             List of Words (imported atom names)

  .. describe:: Output

  Word (Molecule.MolResidue.ccpCodes)
  """
  
  # TODO NBNB Refactor. Now non-proteins also have names of form 'Xyz'
  
  chemComp = None
  
  # get molType
  if not molType:
    molType  = getBestMolType(atomNames)
  
  if molType in ('DNA','RNA','DNA/RNA'):
    if len(resName) == 1:
      if 'PD' in atomNames:
        resName = resName + '11'
  
  # reset character case depending on molType
  if molType == 'protein':
    resName = resName[0] + resName[1:].lower()
  
  chemComp = project.findFirstChemComp(ccpCode=resName.upper(), molType=molType) or \
             project.findFirstChemComp(ccpCode=resName, molType=molType)
  
  if chemComp:
    return chemComp
    
  ccpCodeDict = {}
  chemComps = getStdChemComps(project, molTypes=[molType,])
  
  for chemComp0 in chemComps:
    ccpCodeDict[chemComp0.ccpCode] = chemComp0
    ccpCodeDict[chemComp0.code3Letter] = chemComp0
  
  # Get ChemComp from std dict
  chemComp = ccpCodeDict.get(resName)
  
  # Get ChemComp from std dict using alternative naming systems
  if chemComp is None:
    for chemCompTest in chemComps:
      for namingSystem in chemCompTest.namingSystems:
        for sysName in namingSystem.chemCompSysNames:
          if sysName.sysName == resName:
            #ccpCodeDict[resName] = chemCompTest
            chemComp = chemCompTest
            break
 
      else:
        continue
      break 
    
    # get ChemComp outside std ChemComps, 
    if not chemComp:
      chemComp = getChemComp(project, molType, resName, download=download)
      
    # get ChemComp outside std ChemComp - try with type Other
    if not chemComp and molType != 'other':
      chemComp = getChemComp(project, 'other', resName, download=download)
      if not chemComp:
        resName = resName[0] + resName[1:].lower()
        chemComp = getChemComp(project, 'other', resName)
    
  return chemComp


def getBestMolType(atomNames, ccpCodes=None):
  """Determine the best molecule type (protein, DNA, RNA, carbohydrate or
             nonpolymer) given the input atom names and residue ccpCodes
  .. describe:: Input
  
  List of Words (imported atom names),
             List of Words (Molecule.MolResidue.ccpCodes)

  .. describe:: Output

  Word (Molecule.Molecule.molType)
  """

  molType = 'other'
  
  if ("C3'" in atomNames) and ("C5'" in atomNames) and ("C2" in atomNames):
    molType = 'DNA'
    if "O2'" in atomNames:
      molType = 'RNA'
      
  elif ("C3*" in atomNames) and ("C5*" in atomNames) and ("C2" in atomNames):
    # PDB Naming system different from others
    molType = 'DNA'
    if "O2*" in atomNames:
      molType = 'RNA'
 
  elif 'CA' in atomNames:
    molType = 'protein'

  elif ("C1" in atomNames) and ("C2" in atomNames) and ("C3" in atomNames) and ("C4" in atomNames) \
       and ( ("O2" in atomNames) or ("O3" in atomNames) or ("O4" in atomNames)):
    molType = 'carbohydrate'
  
  return molType
  
  
def findMatchingMolSystemAtom(atomName, residue, namingSystem, excludeAtoms,
                              fixAtomNames=False):
  """
  Find the best matching CCPN atom name in a residue for the input
  atom name in the input naming system.
  Will try other naming systems if the input one doesn't work
  
  .. describe:: Input
  
  Word (imported atom name), MolSystem.Residue,
  Word (ChemComp.NamingSystem.name), List of MolSystem.Atoms

  .. describe:: Output

  Word (MolSystem.Atom.name)
  """
  #nucleic = ('DNA','RNA','DNA/RNA')
  #weirdos = {'O1P':'OP1','O2P':'OP2','C5A':'C7'}
  #if weirdos.get(atomName) and residue.molResidue.molType in nulciec:
  #  atomName = weirdos[atomName]
  
  # If desired change (e.g.) '2HB' to 'HB2'
  if fixAtomNames and atomName[0] in '123':
    # move leading indices to end of name
    atomName = atomName[1:] + atomName[0]
    #print '### newName', atomName
  
  # get list of atomSysNames for preferred NamingSystem and its reference systems
  # First by sysName, then by altSysNames
  atom = None
  chemComp = residue.chemCompVar.chemComp
  namingSystem0 = chemComp.findFirstNamingSystem(name=namingSystem)
  atomSysNames = []
  atomSysNamesAlt = []
  
  usedNamingSystems = set()
  
  if namingSystem0:
    namingSystemR = namingSystem0
    
    while namingSystemR:
      
      usedNamingSystems.add(namingSystemR)
      
      # First choice plain mapping
      for atomSysName in namingSystemR.findAllAtomSysNames(sysName=atomName):
        atomSysNames.append(atomSysName)
          
      # Otherwise try alternative sys names
      for atomSysName0 in namingSystemR.atomSysNames:
        if atomName in atomSysName0.altSysNames:
          atomSysNamesAlt.append(atomSysName0)
          #break
      
      namingSystemR = namingSystemR.atomReference
  
  atomSysNames.extend(atomSysNamesAlt)
  
  moreNamingSystems = [x for x in priorityOrderedNamingSystems(chemComp)
                       if x not in usedNamingSystems]

  # Next chance any naming system plain mapping 
  for namingSystem0 in moreNamingSystems:
    for atomSysName in namingSystem0.findAllAtomSysNames(sysName=atomName):
      atomSysNames.append(atomSysName)

  # Final chance any naming system alt name
  for namingSystem0 in moreNamingSystems:
    for atomSysName in namingSystem0.atomSysNames:
      if atomName in atomSysName.altSysNames:
        atomSysNames.append(atomSysName)
 
  # Find the molSystem atom
  for atomSysName in atomSysNames:
    atom = residue.findFirstAtom(name=atomSysName.atomName)
    if atom and atom not in excludeAtoms:
      return atom
  #   
  #return None
  return residue.findFirstAtom(name=atomName)


def priorityOrderedNamingSystems(chemComp, prefNamingSystemName=None):
  """
  Give naming systems for chemComp in reproducible priority order
  
  .. describe:: Input
  
  ChemComp (ChemComp being named)

  .. describe:: Output

  List of (ChemComp.NamingSystem)
  """
  
  result = []
  
  # set up name:NamingSystem dict
  dd = {}
  for ns in chemComp.namingSystems:
    dd[ns.name] = ns
  
  # Put naming systems on list in priority order
  sortOrder = list(namingSystemPriorityOrder)
  if prefNamingSystemName:
    if prefNamingSystemName in sortOrder:
      sortOrder.remove(prefNamingSystemName)
    sortOrder.insert(0,prefNamingSystemName)
  for name in sortOrder:
    if name in dd:
      result.append(dd.pop(name))
  
  # Put remaining NamingSystems on in alphabetical order
  for name in sorted(dd):
    result.append(dd[name])
  
  #
  return result


  
def getStdResNameMap(chemComps=None):
  """ Generate map of possible string resNames to tuple of (molType, ccpCode) 
  tuples.
  All ccpCodes, code1Letter, code3Letter are used for ChemCompOverview list;
  for passed-in ChemComps also chemCompSysNames are used. 'other' chemcomp
  names are ignored if clashing, unless the ChemComp is passed in.
  NEW code
  
  NBNB TODO change so that D-Ala is prefered over Dal (etc.)

  Input: List of ChemComps

  Output:
   Dictionary {Word:(Word,Word)}
  """

  molTypeOrder = DataMapper.molTypeOrder

  result = {}

  for molType in molTypeOrder:

    tmpResult = {}

    data = copy.deepcopy(chemCompOverview[molType])

    # add data from entered ChemComps
    ccpCodesCC = set()

    # Get stdChemComps first - avoids certain name clash error loading ChemComps
    chemComps = (   [x for x in chemComps if x.className == 'StdChemComp']
                  + [x for x in chemComps if x.className != 'StdChemComp'])
    for chemComp in chemComps:
      if molType == chemComp.molType:

        ccpCode = chemComp.ccpCode
        ccpCodesCC.add(ccpCode)

        # Add ccpCode and main names
        tags = set((ccpCode, ccpCode.upper(), chemComp.code3Letter))
        tag = chemComp.code1Letter
        if tag and chemComp.className == 'StdChemComp':
          tags.add(tag)
        if None in tags:
          tags.remove(None)

        # also sysNames
        for namingSystem in chemComp.namingSystems:
          for sysName in namingSystem.chemCompSysNames:
            tags.add(sysName.sysName)

        # put in mapping
        ccId = (molType, ccpCode)
        for tag in tags:
          prevId = tmpResult.get(tag)
          if prevId is None:
            tmpResult[tag] = ccId
          elif prevId != ccId:
            # Clash. resolve in favour of shortest ccpCode, if any
            # Heuristic- facvours standard bases over phosphorylation states
            if len(prevId[1]) == 1 and len(ccId[1]) > 1:
              continue
            elif len(prevId[1]) > 1 and len(ccId[1]) == 1:
              tmpResult[tag] = ccId
            else:
              raise Exception(
               "Residue name conflict for %s between ChemComps %s and %s"
               % (tag, ccId, prevId))

    # Add data for all chemComps from overview
    for ccpCode,tt in sorted(data.items()):

      if ccpCode in ccpCodesCC:
        # We have this one already, from chemComps. Skip
        continue

      if ccpCode == ccpCode.upper():
        truCode = ccpCode[0].upper() + ccpCode[1:].lower()

        if ccpCode != truCode and truCode in data:
          # We have this one in correct casing elsewhere. Skip it.
          continue

      # get tags to add
      ccId = (molType, ccpCode)
      tags = set((ccpCode,))
      cifCode = tt[1]
      if cifCode:
        truCif = cifCode[0].upper() + cifCode[1:].lower()
        if cifCode not in data and truCif not in data:
          # Skip cifCodes dealt with elsewhere
          tags.add(cifCode)
      if ccpCode.upper() != cifCode:
        tags.add(ccpCode.upper())
      code1L = tt[0]
      if code1L and code1L not in tmpResult:
        tags.add(code1L)

      # put values in result
      for tag in tags:
        prevId = tmpResult.get(tag)
        if prevId is None:
          tmpResult[tag] = ccId
        elif prevId != ccId:
          raise Exception(
           "Residue name conflict for %s between Entries %s and %s"
           % (tag, ccId, prevId), tags)

    # merge tempResult into global result
    for tag, val in tmpResult.items():
      oldval = result.get(tag)
      if oldval is None:
        result[tag] = (val,)
      else:
        # clash between types
        if molType != 'other':
          # Ignore clashes with 'other. NB 'other' must be last in loop.
          result[tag] = oldval + (val,)
  #
  return result
    
      



LICENSE = """
======================COPYRIGHT/LICENSE START==========================

StructureLib.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2013 Rasmus Fogh, Wayne Boucher and Tim Stevens (University of Cambridge)

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
from math import sqrt, cos, sin, atan2

try:
  from memops.gui.MessageReporter import showOkCancel, showYesNo
  from memops.gui.MessageReporter import showError, showWarning
except:
  from memops.universal.MessageReporter import showOkCancel, showYesNo
  from memops.universal.MessageReporter import showError, showWarning

from ccp.lib.MoleculeQuery import getLinkedResidue
from ccp.general.Geometry import calcTorsionAngleRadians, calcTorsionAngleDegrees
                  
BACKBONE_ATOMS = {'protein':('N','C','CA'),
                  'RNA':("OP1","P","O3'","O5'","C3'","C4'","C5'"),
                  'DNA':("OP1","P","O3'","O5'","C3'","C4'","C5'"),
                  'carbohydrate':("C1","C2","C3","C4","C5","C6","C7","C8"),}

TWOPI = 6.2831853071795864

##############################################################
#
### Ensemble handling
#
##############################################################


def makeEmptyEnsembleCopy(ensemble):
  """Make clone of sourceEnsemble with models and data removed
  .. describe:: Input
  
  MolStructure.StructureEmsemble (source)

  .. describe:: Output

  MolStructure.StructureEmsemble (copy)
  """
  
  # create new StructureEnsemble
  newId = max(x.ensembleId for x in ensemble.root.structureEnsembles) + 1
  dd = {}
  for tag in ('atomNamingSystem', 'resNamingSystem', 'softwareName', 'details',
              'molSystem'):
    dd[tag] = getattr(ensemble,tag)
    dd['ensembleId'] = newId
  copyEnsemble = ensemble.root.newStructureEnsemble(**dd)
  
  residueMap = {}
  for atom in ensemble.orderedAtoms:
    
    # get residue
    oldResidue = atom.residue
    residue = residueMap.get(oldResidue)
    if residue is None:
      
      # get or create chain
      chainCode = oldResidue.chain.code
      chain = copyEnsemble.findFirstCoordChain(code=chainCode)
      if chain is None:
        chain = copyEnsemble.newChain(code=chainCode)
      
      # create residue
      dd = {}
      for tag in ('seqId', 'seqCode', 'seqInsertCode'):
        dd[tag] = getattr(oldResidue,tag)
      residue = chain.newResidue(**dd)
      residueMap[oldResidue] = residue
    
    residue.newAtom(name=atom.name, altLocationCode=atom.altLocationCode)
  #
  return copyEnsemble


def copyModelToEnsemble(ensemble, model, identicalEnsembles=False):
  """ Create copy of model within ensemble. All Chains, Residues and Atoms in
  model must match Ensemble, but Model may contain fewer Atoms.
  If identicalEnsembles is True, atom ordering is assumed to be identical.
  
  Input: 
         ensemble: molecule.MolStructure.StructureEnsemble,
         model: molecule.MolStructure.Model
         identicalEnsembles: boolean
         
  """
  
  if model.structureEnsemble is ensemble:
    showWarning('%s already belongs to target ensemble'  % model)
    return

  if model.structureEnsemble.molSystem is not ensemble.molSystem:
    showError('Ensemble Creation Failed','Molecular System mismatch')
  
  if identicalEnsembles:
    # identical ensembles - we trust tha atom ordering is identical.
  
    if model.structureEnsemble.nAtoms != ensemble.nAtoms:
      showError('Ensemble Creation Failed',
                  "nAtoms differs for 'identical' ensembles")
      
    newCoords = model.coordinates
    newOccs = model.occupancies
    newBFacs = model.bFactors
  
  else:
    #map old data to new - NB the atom ordering may be different
    
    newCoords = [0.0] * (3 * ensemble.nAtoms)
    newBFacs = [0.0] *  ensemble.nAtoms
    newOccs = [0.0] * ensemble.nAtoms
 
    oldCoords = model.coordinates
    oldBFacs = model.bFactors
    oldOccs = model.occupancies
 
    for ii,atom in enumerate(model.structureEnsemble.orderedAtoms):
 
      residue = atom.residue
 
      newAtom = None
      newChain = ensemble.findFirstCoordChain(code=residue.chain.code)
      if newChain is not None:
        newRes = newChain.findFirstResidue(seqId=residue.seqId)
        if newRes is not None:
          newAtom = newRes.findFirstAtom(name=atom.name,
                                         altLocationCode=atom.altLocationCode)
        
        # NB we do not compare residue types
        #  - if the model contains non-matching atoms the error will come there
        
      if newAtom is None:
        showError('Ensemble Creation Failed',
                  "No match found for atom %s" % atom)
        return
 
      else:
        # set new values. NB atoms may not be in same order in the two ensembles
        newii = newAtom.index
        newOccs[newii] = oldOccs[ii]
        newBFacs[newii] = oldBFacs[ii]
        newOff = 3 * newii
        oldOff = 3 * ii
        newCoords[newOff:newOff+3] = oldCoords[oldOff:oldOff+3]
  
  newModel = ensemble.newModel(name=model.name, details=model.details)
  newModel.setSubmatrixData('coordinates', newCoords)
  newModel.setSubmatrixData('occupancies', newOccs)
  newModel.setSubmatrixData('bFactors', newBFacs)
  

def makeEnsemble(models, ensemble=None, replaceStructures=True):
  """Combine the input models into a single structure ensemble. If ensenble
             is set, it will be overwritten and extraneous models deleted, 
             otherwise a new Ensemble will be made. If replaceStructures,
             input models will be deleted, as will (now) empty input ensembles.
  .. describe:: Input
  
  List of MolStructure.Models,  MolStructure.SructureEnsemble or None,
  Boolean

  .. describe:: Output

  MolStructure.StructureEmsemble
  """
  if not models:
    return
  
  # select or create target ensemble; set up control info
  if ensemble is None:
    refEnsemble = models[0].structureEnsemble
    ensemble = makeEmptyEnsembleCopy(refEnsemble)
    nStartModels = 0
    extraModels = []
  else:
    refEnsemble = None
    ll = ensemble.sortedModels()
    nStartModels = len(ll)
    extraModels = [x for x in ll if x not in models]
  
  # add models
  try:
    for model in models:
      ens = model.structureEnsemble
      if ens is ensemble:
        continue
      elif ens is refEnsemble:
        copyModelToEnsemble(ensemble, model, identicalEnsembles=True)
      else:
        copyModelToEnsemble(ensemble, model)
        
  except:
    # Model merge failed - undo changes
    
    if refEnsemble is None:
      # We are overwriting an existing ensemble - remove added models
      ll = ensemble.sortedModels()[nStartModels:]
      for model in reversed(ll):
        model.delete()
    
    else:
      # ensemble was created afresh - delete it
      ensemble.delete()
    
    # re-raise original error
    raise
  
  # Model merge succeeded
  
  # Remove extraneous models from result ensemble
  for model in reversed(extraModels):
    model.delete()
 
  # Cleanup old structures if required
  if replaceStructures:
    for model in models :
      ens = model.structureEnsemble
      if ens is not ensemble:
        model.delete()
        if not ens.models:
          ens.delete()
  #
  return ensemble

##############################################################
#
### Coordinate handling, various
#
##############################################################


def alignStructures(structures):
  """Align structures by minimising weighted atomic RMSD.
             All members of any entered ensembles are aligned.

  .. describe:: Input
  
  List of MolStructure.StructureEnsembles

  .. describe:: Output

  List of MolStructure.StructureEnsembles, Float (Fit error), 
             Float (overall RMSD), List of Floats (Atom RMSDs)
  """

  # TBD: Arbitrary atom selection
  # TBD: Non-equivalent atom mappings
  # NBNB consider refactoring to avoid use of Coords
  
  from ccp.util.Validation import storeResidueValidations, storeModelValidations
  from ccpnmr.analysis.core.ValidationBasic import getEnsembleValidationStore
  from ccpnmr.analysis.core.ValidationBasic import ANALYSIS_RMSD_CONTEXT, RMSD_KEYWORDS
  from ccp.c.StructUtil import alignEnsemble

  structures = list(structures)
  
 
  null = structures[0].root.findFirstChemElementStore()

  nModels = sum([len(s.models) for s in structures])

  if nModels < 2:
    return structures, 0.0, [], {}

  numAtoms  = None
  ensemble  = []
  weights   = []
  allCoords = []
  
  k = 0;
  for structure in structures:
    coordLists = _getStructureCoordinates(structure, model=None)
    coordList = coordLists[0]
    
    if numAtoms is None:
      numAtoms = len(coordList)
      
    elif numAtoms != len(coordList):
      raise Exception('Attempt to align structures with different atoms')

    if not weights:
      for coordObj in coordList:
        chemAtom = coordObj.atom.chemAtom
        mass = chemAtom.chemElement.mass
        weights.append(min(14.0,mass))

    for coordObjs in coordLists:
      coords = [[0.0,0.0,0.0]] * numAtoms
      
      for i in xrange(numAtoms):
        coordObj = coordObjs[i]
        coords[i] = [coordObj.x, coordObj.y, coordObj.z]
      
      allCoords.append(coordObjs)
      ensemble.append(coords)
  
  # Call to the C code to do the alignment
  error, atomRmsds, structureRmsds = alignEnsemble(ensemble, weights)
  weights = [x**-2.0 for x in atomRmsds]
  error, atomRmsds, structureRmsds = alignEnsemble(ensemble, weights)
  
  models = []
  j = 0
  for structure in structures:
    models = structure.sortedModels()
    k = j+len(models)
    scores = structureRmsds[j:k]
    validStore = getEnsembleValidationStore(structure,
                                            ANALYSIS_RMSD_CONTEXT,
                                            RMSD_KEYWORDS)
    storeModelValidations(validStore, ANALYSIS_RMSD_CONTEXT,
                          RMSD_KEYWORDS[1], models, scores)  
    j = k
      
  # Fill in the data model coords
  for i, alignVals in enumerate(ensemble):
    coordObjs = allCoords[i]
    
    for j in xrange(numAtoms):
      coordObj = coordObjs[j]
      coordObj.x, coordObj.y, coordObj.z = alignVals[j]
  
  atomRmsdDict = {}
  for i, coord in enumerate(allCoords[0]):
    atomRmsdDict[coord.atom] = atomRmsds[i]
  
  for structure in structures:
    
    # Reset cache
    structure.coordDict = {}

    residues = []
    caRmsds  = []
    cbRmsds  = []
    coRmsds  = []
    hnRmsds  = []
    bbRmsds  = []
    #scRmsds  = []
    aaRmsds  = []
    for chain in structure.coordChains:
      for residue in chain.residues:
        molType = residue.residue.molResidue.molType
        bbAtoms = BACKBONE_ATOMS.get(molType, [])
        
        caRmsd = None
        cbRmsd = None
        coRmsd = None
        hnRmsd = None
        bbRmsd = []
        #scRmsd = []
        aaRmsd = []
        for atom in residue.atoms:
          atomName = atom.name
          rmsd = atomRmsdDict[atom]
          
          if atomName == 'CA':
            caRmsd = rmsd
          elif atomName in 'O':
            coRmsd = rmsd
          elif atomName == 'CB':
            cbRmsd = rmsd
          elif atomName == 'H':
            hnRmsd = rmsd
            
          if atomName in bbAtoms:
            bbRmsd.append(rmsd)
            
          aaRmsd.append(rmsd)
            
        cAlpha = residue.findFirstAtom(name='CA')
        
        if aaRmsd:
          residues.append(residue)
          caRmsds.append(caRmsd)
          cbRmsds.append(cbRmsd)
          coRmsds.append(coRmsd)
          hnRmsds.append(hnRmsd)
          bbRmsds.append(sum(bbRmsd)/(len(bbRmsd) or 1.0))
          aaRmsds.append(sum(aaRmsd)/(len(aaRmsd) or 1.0))
        
    validStore = getEnsembleValidationStore(structure,
                                            ANALYSIS_RMSD_CONTEXT,
                                            RMSD_KEYWORDS)
    
      
     
    storeResidueValidations(validStore, ANALYSIS_RMSD_CONTEXT, 'CA', residues, caRmsds)
    storeResidueValidations(validStore, ANALYSIS_RMSD_CONTEXT, 'O', residues, coRmsds)
    storeResidueValidations(validStore, ANALYSIS_RMSD_CONTEXT, 'CB', residues, cbRmsds)
    storeResidueValidations(validStore, ANALYSIS_RMSD_CONTEXT, 'H', residues, hnRmsds)
    storeResidueValidations(validStore, ANALYSIS_RMSD_CONTEXT, 'backbone', residues, bbRmsds)
    #storeResidueValidations(validStore, ANALYSIS_RMSD_CONTEXT, 'sidechain', residues, scRmsds)
    storeResidueValidations(validStore, ANALYSIS_RMSD_CONTEXT, 'all', residues, aaRmsds)
       
  return structures, error, structureRmsds, atomRmsdDict


def compareEnsembles(structure1, structure2, compareBackboneOnly=False):
  """Compare models in two ensembles by minimising weighted atomic RMSD.
     The alignment is only done in the C world, not in the Python world.
     If compareBackboneOnly then only backbone atoms are used.

  .. describe:: Input
  
  MolStructure.StructureEnsemble, MolStructure.StructureEnsemble, Boolean

  .. describe:: Output

  Float (overall RMSD), Dict: MolSystem.Residue --> RMSD
  """

  from ccp.c.StructUtil import alignEnsemble

  # make sure ChemElements are loaded??
  null = structure1.root.findFirstChemElementStore()

  atomNameDict = compareBackboneOnly and BACKBONE_ATOMS

  numAtoms  = None
  weights   = []
  coordsDict = {}
  coordListDict = {}
  
  allCoordLists = _getMatchedStructureCoordinates(structure1, structure2, 
                                                          atomNameDict)
  
  structures = (structure1, structure2)
  for ii in (0,1):
    structure = structures[ii]
    allCoordList = allCoordLists[ii]
    for jj,model in enumerate(structure.sortedModels()):
      coordList = allCoordList[jj]
 
      # Should bo longer be necessary, but does no harm
      if numAtoms is None:
        numAtoms = len(coordList)
      elif numAtoms != len(coordList):
        raise Exception('Attempt to compare structures with different atoms')

      if not weights:
        for coordObj in coordList:
          chemAtom = coordObj.atom.chemAtom
          mass = chemAtom.chemElement.mass
          weights.append(min(14.0,mass))

      coords = [[0.0,0.0,0.0]] * numAtoms
      for i in xrange(numAtoms):
        coordObj = coordList[i]
        coords[i] = [coordObj.x, coordObj.y, coordObj.z]
      coordsDict[model] = coords
      coordListDict[model] = coordList

  residueRmsdDict = {}
  residueCountDict = {}
  ensemble  = 2*[0]
  for model1 in structure1.sortedModels():
    ensemble[0] = coordsDict[model1]
    allCoord = coordListDict[model1]
    for model2 in structure2.sortedModels():
      ensemble[1] = coordsDict[model2]

      # Call to the C code to do the alignment
      error, atomRmsds, structureRmsds = alignEnsemble(ensemble, weights)
      # NB this breaks if an atomRmsd is 0, which might happen for 2 models only
      #weights = [x**-2.0 for x in atomRmsds]
      # Set to defWeight if atomRmsd is 0
      # Maybe adaptive weights will not work for just two structures. Trynwithout
      #defWeight = 1.0
      #weights = [(x and x**-2.0) or defWeight for x in atomRmsds]
      #print '~~~2', structure1, model1.serial, structure2, model2.serial
      #error, atomRmsds, structureRmsds = alignEnsemble(ensemble, weights)
      
      # atomRmsds are out by sqrt(2) because C code
      # divides by nensembles not nensembles-1
      for i, coord in enumerate(allCoord):
        rmsd = atomRmsds[i]
        residue = coord.atom.residue.residue
        residueRmsdDict[residue] = residueRmsdDict.get(residue, 0) + 2 * rmsd * rmsd
        residueCountDict[residue] = residueCountDict.get(residue, 0) + 1

  rmsd2Tot = 0
  countTot = 0
  for residue in residueRmsdDict:
    rmsd2 = residueRmsdDict[residue]
    rmsd2Tot += rmsd2
    count = residueCountDict[residue]
    countTot += count
    residueRmsdDict[residue] = sqrt(rmsd2 / count)

  if countTot > 0:
    rmsd = sqrt(rmsd2Tot / countTot)

  return rmsd, residueRmsdDict


  
def _getMatchedStructureCoordinates(structure1, structure2, atomNameDict=None):
  """Get a list of coordinate objects from two structures in matching order, skipping unmatching atoms
             Option to extract only atoms of certain names
  .. describe:: Input
  
  MolStructure.StructureEnsemble, MolStructure.StructureEnsemble

  .. describe:: Output

  2 List of List of MolStructure.Coords (1st list per atom 2nd per model)
  
  NBNB consider refactoring to avoid use of Coords
  """
  
  atoms1 = structure1.orderedAtoms
  atoms2 = structure2.orderedAtoms
  
  # Make atomkey:index dict for structure1, skiping unwanted names
  matchDict1 = {}
  for ii,atom in enumerate(atoms1):
    res = atom.residue
    name = atom.name
    if atomNameDict:
      atomNames = atomNameDict.get(atom.residue.residue.molType, [])
      if atom.name not in atomNames:
        continue
    tt = (name, res.seqId, res.chain.code)
    matchDict1[tt] = ii
  
  # Make index(struc1):index(struc2) for atoms present in both
  matchDict2 = {}
  for ii,atom in enumerate(atoms2):
    res = atom.residue
    tt = (atom.name,res.seqId, res.chain.code)
    jj = matchDict1.get(tt)
    if jj is not None:
      matchDict2[jj] = ii
  
  #print '### matchdicts', len(matchDict1), len(matchDict2)
  
  # make coords for structure 1
  models = structure1.sortedModels()
  coordList1 = [[] for n in models]
  for ii in sorted(matchDict2):
    atom = atoms1[ii]
    for kk, model in enumerate(models):
      coord = atom.newCoord(model=model)
      coordList1[kk].append(coord)
  
  # make coords for structure 2
  models = structure2.sortedModels()
  coordList2 = [[] for n in models]
  for dummy,ii in sorted(matchDict2.items()):
    atom = atoms2[ii]
    for kk, model in enumerate(models):
      coord = atom.newCoord(model=model)
      coordList2[kk].append(coord)

  return coordList1, coordList2
  
def _getStructureCoordinates(structure, model=None, atomNameDict=None):
  """Get a list of coordinate objects from a structure in a consistent order.
             Option to select coords from a given model if required.
             Option to extract only atoms of certain names
  .. describe:: Input
  
  MolStructure.StructureEnsemble, MolStructure.Model, Dict {molTypeSTring:[list of atomNames]}

  .. describe:: Output

  List of List of MolStructure.Coords (1st list per atom 2nd per model)
  
  NBNB consider refactoring to avoid use of Coords
  """

  if model:
    models = [model,]
  else:  
    models = structure.sortedModels()

  coordList = [[] for n in models]

  for atom in structure.orderedAtoms:
    atomNames = atomNameDict and atomNameDict.get(atom.residue.residue.molType, [])
    if atomNameDict and atom.name not in atomNames:
      continue
    for i, model in enumerate(models):
      coord = atom.newCoord(model=model)
      coordList[i].append(coord)

  return coordList



def getResiduePhiPsi(residue, inDegrees=True, model=None):
  """Find the Phi and Psi backbone dihedral angles for a residue in a structure.
             Option inDegrees can be set to False to get an angle in radians.
             Option to specify which model of an ensemble to use, otherwise
             all models are considered and the angles are an average.
  .. describe:: Input
  
  MolStructure.Residue, Boolean, MolStructure.Model

  .. describe:: Output

  2-List of Floats (Phi, Psi)
  NBNB consider refactoring to avoid use of Coords
  """

  chain = residue.chain
  sysResidue = residue.residue
  
  if not sysResidue:
    return None, None

  if model:
    models = [model,]
    structure = model.structureEnsemble
  else:
    structure = chain.structureEnsemble
    models = structure.sortedModels() 

  N = float(len(models))

  sysResiduePrev = getLinkedResidue(sysResidue, linkCode='prev')
  sysResidueNext = getLinkedResidue(sysResidue, linkCode='next')

  if not (sysResiduePrev and sysResidueNext):
    return None, None
    
  meanPhi = None
  meanPsi = None
  residuePrev = chain.findFirstResidue(seqId=sysResiduePrev.seqId)
  residueNext = chain.findFirstResidue(seqId=sysResidueNext.seqId)
  
  if residuePrev and residueNext:
    atomC0 = residuePrev.findFirstAtom(name='C')
    atomN  = residue.findFirstAtom(name='N')
    atomCa = residue.findFirstAtom(name='CA')
    atomC  = residue.findFirstAtom(name='C')
    atomN2 = residueNext.findFirstAtom(name='N')
 
    if atomC0 and atomN and atomCa and atomC and atomN2:
      angles = []
      
      for model in models:
        coords = []
        for atom in (atomC0, atomN, atomCa, atomC, atomN2):
          coordObj = atom.newCoord(model=model)
          coords.append( [coordObj.x,coordObj.y,coordObj.z] )
 
        if inDegrees:
          phi = calcTorsionAngleDegrees(coords[0],coords[1],coords[2],coords[3])
          psi = calcTorsionAngleDegrees(coords[1],coords[2],coords[3],coords[4])
 
        else:
          phi = calcTorsionAngleRadians(coords[0],coords[1],coords[2],coords[3])
          psi = calcTorsionAngleRadians(coords[1],coords[2],coords[3],coords[4])
 
        angles.append((phi,psi))
 
      if N == 1.0:
        meanPhi, meanPsi = angles[0]
        
      else:
        sumCosPhi = 0.0
        sumSinPhi = 0.0
        sumCosPsi = 0.0
        sumSinPsi = 0.0
        
        for phi, psi in angles:
        
          if inDegrees:
            phi = TWOPI * phi/360.0
            psi = TWOPI * psi/360.0

          sumCosPhi += cos(phi)
          sumSinPhi += sin(phi) 
          sumCosPsi += cos(psi)
          sumSinPsi += sin(psi) 

        meanPhi = atan2(sumSinPhi/N,sumCosPhi/N)
        meanPsi = atan2(sumSinPsi/N,sumCosPsi/N)
        
        if inDegrees:
          meanPhi *= 360/TWOPI
          meanPsi *= 360/TWOPI
 
  return meanPhi, meanPsi
  

def getAtomSetCoords(atomSet, structure, model=None):
  """Find the coordinates corresponding to an NMR atom set in a given model
  of a given structure structure. The model defaults to an arbitrary one
  if none is specified.
  
  .. describe:: Input
  
  Nmr.AtomSet, MolStructure.StructureEnsemble, MolStructure.Model

  .. describe:: Output

  List of List of Floats (x,y,z for each atom)
  NBNB Maybe consider refactoring to avoid use of Coords
  """

  if not model:
    model = structure.findFirstModel()

  key = '%s:%s:%d:%d' % (atomSet, # Could be real AtomSet or a ConstraintSet one!
                         structure.molSystem.code,
                         structure.ensembleId,
                         model.serial)

  if not hasattr(structure,'coordDict'):
    structure.coordDict = {}
    
  else:
    coordList = structure.coordDict.get(key)
 
    if coordList:
      return coordList

  atom    = atomSet.findFirstAtom()
  residue = atom.residue
  chain   = residue.chain
  coordChain = structure.findFirstCoordChain(code=chain.code)
  if not coordChain:
    #showWarning('Warning', 'Couldn\'t find coordinate chain')
    return []
  
  coordResidue = coordChain.findFirstResidue(seqId=residue.seqId)
  if not coordResidue:
    data = (residue.ccpCode,residue.seqCode)
    msg  = 'Couldn\'t find coordinate residue %s %d' % data
    print msg
    #showWarning('Warning', msg)
    return []

  coordList = []
  findAtom = coordResidue.findFirstAtom

  for atom in atomSet.atoms:
    coordAtom = findAtom(name=atom.name)
      
    if coordAtom:
      coord = coordAtom.newCoord(model=model)
      coordList.append(coord)

  if not coordList:
    data = (chain.code, residue.ccpCode, residue.seqCode, atomSet.name)
    print 'Couldn\'t find coordinate atoms %s %s %d %s' % data
    return []
  
  structure.coordDict[key] = coordList
  
  return coordList

def getAtomSetsDihedral(atomSets, structure, model=None, inDegrees=True):
  """Measure the dihedral angle in a structure between four groups
             of atoms (lists of atom sets). Input atom sets are a list of lists
             to allow for ambigous assignments (e.g. Ser Hb*) here the mean position is used.
             If mo nodel is specified the values are averates over the whole ensemble.
             Option to specify degrees or radians. 
  .. describe:: Input
  
  4 List of List of Nmr.AtomSets or 4 List of Nmr.AtomSets,
             MolStructure.StructureEnsemble, Int, Boolean

  .. describe:: Output

  Float 
  """

  assert len(atomSets) == 4

  if model:
    models = [model,]
  else:
    models = structure.sortedModels()

  for i, atomSet in enumerate(atomSets):
    if hasattr(atomSet, 'className') and (atomSet.className == 'AtomSet'):
      atomSets[i] = [atomSet,]
 
  angles = []
  
  for model in models:
    coords = []
    
    for atomSetList in atomSets:
      coord = [0.0,0.0,0.0]
      n = 0.0
      # average over ambiguous atom sets e.g. Val Hg1* & Hg2*
      for atomSet in atomSetList:
        # average over atoms in atom set e.g. methyl
        for coord0 in getAtomSetCoords(atomSet, structure, model):
          coord[0] += coord0.x
          coord[1] += coord0.y
          coord[2] += coord0.z
          n += 1.0
      
      if n:
        coords.append([v/n for v in coord])

    if coords and len(coords) == 4:
      angle = calcTorsionAngleRadians(coords[0],coords[1],coords[2],coords[3])
      angles.append(angle)

  meanVal = None
  sumCos  = 0.0
  sumSin  = 0.0
  
  N = float(len(angles))
  for a in angles:
    sumCos += cos(a)
    sumSin += sin(a) 

  if N:
    sumCos /= N
    sumSin /= N
    meanVal = atan2(sumSin/N,sumCos/N)

    if inDegrees:
      meanVal *= 360/TWOPI

  return meanVal

def getAtomSetsDistance(atomSets1, atomSets2, structure, model=None, method='noe'):
  """
  Find the distance between two atom sets in a specified structure or ensemble
  of structures. Distances for multi atom atom sets are calculated using the
  NOE sum method by default. A model may be specified if the structure
  ensemble has many models, otherwise all models in the ensemble will be
  considered and an avergate distance is returned. The method can be either
  "noe" for NOE sum, "min" for minimum distance or "max" for maximum distance.
  
  .. describe:: Input
  
  Nmr.AtomSet, Nmr.AtomSet, MolStructure.StructureEnsemble,
  MolStructure.Model or None, Word

  .. describe:: Output

  Float (distance)
  """
  
  if model:
    models = [model,]
  else:
    models = structure.models 

  atomSets1 = set(atomSets1)
  atomSets2 = set(atomSets2)

  assert atomSets1 and atomSets2 and structure
  
  if atomSets1 == atomSets2:
    if len(atomSets1) == 2:
      atomSets0 = list(atomSets1)
      for resonanceSet in atomSets0[0].resonanceSets:
        if resonanceSet in atomSets0[1].resonanceSets: # Prochiral
          atomSets1 = [atomSets0[0],]
          atomSets2 = [atomSets0[1],]
          break
    
    if atomSets1 == atomSets2: # Not prochiral
      return 0.0

  ensembleCoords = []
  ensembleCoordsAppend = ensembleCoords.append
  
  for model in models:
    coordList1 = []
    coordList2 = []
    coordList1Append = coordList1.append
    coordList2Append = coordList2.append
    
    for atomSet in atomSets1:
      for coord in getAtomSetCoords(atomSet, structure, model):
        coordList1Append(coord)
      
    for atomSet in atomSets2:
      for coord in getAtomSetCoords(atomSet, structure, model):
        coordList2Append(coord)

    if coordList1 and coordList2:
      ensembleCoordsAppend((coordList1,coordList2))
 
  if not ensembleCoords:
    return 0.0
  
  numPairs = float(len(atomSets1)*len(atomSets2))  
  noeEnsemble = 0.0
  minDist2 = None
  maxDist2 = None

  n = 0.0
  for coords1, coords2 in ensembleCoords:
    noeSum = 0.0

    for coord1 in coords1:
      x = coord1.x
      y = coord1.y
      z = coord1.z
      
      for coord2 in coords2:
        dx = x - coord2.x
        dy = y - coord2.y
        dz = z - coord2.z
        dist2 = (dx*dx) + (dy*dy) + (dz*dz)
        
        if dist2 > 0:
          
          if method == 'noe':
            noeSum += dist2 ** -3.0
          
          else:
            if (minDist2 is None) or (dist2 < minDist2):
              minDist2 = dist2

            if (maxDist2 is None) or (dist2 > maxDist2):
              maxDist2 = dist2

    noeEnsemble += sqrt(noeSum/numPairs)
    n += 1.0

  if method == 'min':
    return sqrt(minDist2)
  elif method == 'max':
    return sqrt(maxDist2)
  elif (noeEnsemble > 0.0) and n:
    noeEnsemble /= n
    return noeEnsemble ** (-1/3.0)


def checkChemAtomsForDeepSet(deepChemAtomSet,atomSets):
  """ Check if ChemComp.ChemAtomSet is suitable to cover atomSets
  """

  chemAtomsFound = True
                
  deepCasChemAtoms = []
  for tcas in deepChemAtomSet.chemAtomSets:
    deepCasChemAtoms.extend(list(tcas.chemAtoms))

  currentChemAtoms = []
  for tas in atomSets:
    for ta in tas.atoms:
      if ta.chemAtom not in currentChemAtoms:
        currentChemAtoms.append(ta.chemAtom)

  for tca in deepCasChemAtoms:
    if tca not in currentChemAtoms:
      chemAtomsFound = False
      break
      
  return chemAtomsFound

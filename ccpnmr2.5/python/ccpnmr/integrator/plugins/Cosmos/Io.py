
"""
======================COPYRIGHT/LICENSE START==========================

Io.py: code for CCPN data model and code generation framework

Copyright (C) 2012  (CCPN Project)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

"""
from ccpnmr.analysis.core import AssignmentBasic
from ccpnmr.analysis.core import ConstraintBasic
from ccp.lib.StructureLib import makeEmptyEnsembleCopy
from ccpnmr.integrator.core import Util as intUtil

bondTypeMap = {
 'single':1,
 'singleplanar':1,
 'double':2,
 'triple':3,
 'aromatic':2,
}

interfaceData = [
('Distance Restraints', 'Distance', 'R_CONSTRAINTS', 2),
('Chemical Shifts', 'ChemShift', 'CS_VALUES', 1),
('J Couplings', 'JCoupling', 'J_COUPLINGS', 2),
('RDC Couplings', 'Rdc', 'RDC_VALUES', 2),
]

def makeAtomInput(constraintList, atomIdMap):
  """ Make COSMOS data from input constraint list
  atomIdMap maps atoms to atom IDs
  """
  
  result = []
  
  for constraint in constraintList.sortedConstraints():
    
    data = []
    
    assignments = ConstraintBasic.getConstraintAtoms(constraint)
    if len(assignments) == 1:
      # only unambiguous assignments accepted
      resonances = assignments[0]
      
      # number of resonances correct for data type (none is missing)
      for atoms in resonances:
        # NBNB only single atoms accepted. Prochirals? Methyl groups?
        if len(atoms) == 1:
          ii = atomIdMap.get(atoms[0])
          if ii is None:
            break
          else:
            data.append(ii)
        else:
          break
          
      else:
        # We found one atom for each resonance.
        data.append(constraint.targetValue)
        data.append(constraint.error or -1.0)
        result.append(data)
  #
  return result
  
  
'''
def makeShiftData(shiftList, atomIdMap):
  """ Make COSMOS shift data from shift list
  
  NBNB TBD. Just matches shift with first atom that fits
            No proper handling of prochirals
            No proper handling of methyl or other equivalent groups
  NBNB TBD Would be faster to start from Shift and find Atom
  """
  
  result = []
  
  usedShifts = set()
  for atom in atomIdMap:
    atomSet = atom.atomSet
    if atomSet is not None:
      shifts = AssignmentBasic.getAtomSetShifts(atomSet,shiftList)
      if len(shifts) in (1,2):
        for shift in shifts:
          if shift not in usedShifts:
            usedShifts.add(shift)
            result.append((atomIdMap[atom],shift.value))
            break
  #
  return result
'''


def getAtomData(ensemble):
  """ Get Atom data for ensemble
  returns List:atomNames, List:elementNumbers, Dict{MolSystem.atom:Id}:atomIdMap
  """
  
  # get name stubs for residues - done once for speed
  resNames = {}
  for chain in ensemble.coordChains:
    code = chain.code.strip()
    for res in chain.residues:
      seqId = res.seqId
      ccpCode = res.residue.ccpCode
      if code:
        text = '_%s_%s_%s' % (ccpCode,seqId,code)
      else:
        text = '_%s_%s' % (ccpCode,seqId)
      resNames[res] = text
  
  # element number cache
  elemNumberCache = {}
  
  # initialise results
  atomNames = []
  elementNumbers = []
  atomIdMap = {}
  
  # set result
  for ii,atom in enumerate(ensemble.orderedAtoms):
    
    # process atom name
    name = atom.name + resNames[atom.residue]
    if len(name) > 17:
      # Unlikely. Remove ccpCode if name is too long
      ll = name.split('_')
      del ll[1]
      name = '_'.join(ll)
    atomNames.append(name)
    atomIdMap[atom.atom] = ii + 1
    
    # process element number
    elementSymbol = atom.elementSymbol
    elemNumber = elemNumberCache.get(elementSymbol)
    if elemNumber is None:
      # Use cache to save time
      elemNumber = atom.atom.chemAtom.chemElement.atomNumber
      elemNumberCache[elementSymbol] = elemNumber
    elementNumbers.append(elemNumber)
  #
  return atomNames,elementNumbers,atomIdMap


def getBondData(molSystem, atomIdMap):
  """ Get lit of (atomId1, atomI2, bondType) for Cosmos input
  """
  
  result = []
  
  for chain in molSystem.sortedChains():
    for residue in chain.sortedResidues():
      ccv = residue.chemCompVar
      for bond in ccv.sortedChemBonds():
        data = []
        for chemAtom in bond.chemAtoms:
          if chemAtom.className == 'ChemAtom':
            # normal atom
            atom = residue.findFirstAtom(name=chemAtom.name)
            atomId = atomIdMap.get(atom)
            if atomId is None:
              break
            else:
              data.append(atomId)
            
          else:
            # this is a linkAtom. Links are handled below
            break
            
        else:
          data.append(bondTypeMap[bond.bondType])
          result.append(data)
      
    # do links for Molecule
    for link in chain.molecule.molResLinks:
      data = []
      for linkEnd in link.molResLinkEnds:
        residue = chain.findFirstResidue(seqId=linkEnd.molResidue.serial)
        ccLinkEnd = linkEnd.linkEnd
        atom = residue.findFirstAtom(name=ccLinkEnd.boundChemAtom.name)
        atomId = atomIdMap.get(atom)
        if atomId is None:
          break
        else:
          data.append(atomId)
      else:
        chemComp = ccLinkEnd.chemComp
        bond = chemComp.findFirstChemBond(
         chemAtoms=(ccLinkEnd.boundChemAtom,ccLinkEnd.boundLinkAtom))
        data.append(bondTypeMap[bond.bondType])
        result.append(data)
        
  for link in molSystem.molSystemLinks:
    data = []
    for linkEnd in link.molSystemLinkEnds:
      residue = linkEnd.residue
      ccLinkEnd = linkEnd.linkEnd
      atom = residue.findFirstAtom(name=ccLinkEnd.boundChemAtom.name)
      atomId = atomIdMap.get(atom)
      if atomId is None:
        break
      else:
        data.append(atomId)
    else:
      chemComp = ccLinkEnd.chemComp
      bond = chemComp.findFirstChemBond(
       chemAtoms=(ccLinkEnd.boundChemAtom,ccLinkEnd.boundLinkAtom))
      data.append(bondTypeMap[bond.bondType])
      result.append(data)
  #
  return result


def getCosmosData(nmrCalcRun):
  """ Make data{} and options{} dictionaries for COSMOS input
  """
  
  # get options - use all input options that have a 'code' attribute
  options = {}
  for runParameter in nmrCalcRun.findAllRunParameters(ioRole='input'):
    code = runParameter.code
    if code is not None:
      options[code] = intUtil.getParameterValue(runParameter)
  
  # get data
  data = {}
  
  # get coordinates
  coordinates = data['COO'] = []
  ll = nmrCalcRun.findAllData(className='StructureEnsembleData', ioRole='input')
  if len(ll) == 1:
    ensemble = ll.pop().structureEnsemble
    atomNames,elementNumbers,atomIdMap = getAtomData(ensemble)
    # list of name,elemNumber tuples
    atomData = zip(atomNames, elementNumbers) 
    for model in ensemble.sortedModels():
      # Split into x,y,z triplets
      coordTriplets = zip(*[iter(model.coordinates)]*3) 
      coordinates.append([ list(atomData[ii] + xyz)
                          for ii,xyz in enumerate(coordTriplets)])
    
  else:
    raise Exception(
     "Run must have exactly one input structure ensemble, %s found" 
    % len(ll))
  
  # get bond data
  data['BONDS']= getBondData(ensemble.molSystem, atomIdMap)
  
  # get experiment data
  data['EXP'] = expdict = {}
  for label, name, code, numResonances in interfaceData:
    constraintData = nmrCalcRun.findFirstData(code=code)
    if constraintData is not None:
      # There should be exactly one constraintList
      constraintList = constraintData.findFirstConstraintList()
      expdict[code] = makeAtomInput(constraintList, atomIdMap)
  
  return data, options


def mergeCosmosData(nmrCalcRun, cosmosOutput):
  """ Merge Cosmos output data back into NmrCalcRun
  """
  summaryLog = cosmosOutput[0]
  modelData = cosmosOutput[1:]
  
  nmrCalcRun.newRunParameter(name='summaryLog', code='SUMMARY_LOG', 
                             ioRole='output', textValue=summaryLog)
  
  
  
  if modelData:
    
    inputEnsembleData = nmrCalcRun.findFirstData(ioRole='input',
                                             name='structureEnsemble')
    if inputEnsembleData is None:
      raise KeyError("NmrCalc.run has no structureEnsemble input Data object")
    
    newEnsemble = makeEmptyEnsembleCopy(inputEnsembleData.structureEnsemble)
    nmrCalcRun.newStructureEnsembleData(ioRole='output', code='COO',
                                        name='structureEnsemble',
                                        structureEnsemble=newEnsemble)
    
    nAtoms = newEnsemble.nAtoms
    constraintStore = None
    
    for modelDict in modelData:
      
      # Models and coordinates
      # NBNB problem matching output to input models. NB lacks information.
      model = newEnsemble.newModel()
      
      # NBNB TBD some better model identifier would be useful
      modelSerial = model.serial
      
      newData = nmrCalcRun.newStructureEnsembleData(ioRole='output',
                  #groupId=modelSerial, models=(model,),
                  models=(model,), name='model%s' % modelSerial)
      
      # COSMOS log
      modelLog = modelDict.get('LOG')
      if modelLog is not None:
        nmrCalcRun.newRunParameter(name='log%s' % modelSerial,
                                   code='LOG', ioRole='output', 
                                   textValue=modelLog, data=newData)
                                   #textValue=modelLog, groupId=modelSerial)
      
      coords = modelDict.get('COO')
      if coords:
        # NB assumes that atoms are in same order as input
        # and that coordinates are always present.
        coordinates = [None] * (nAtoms * 3)
        offset = 0
        for tt in coords:
          coordinates[offset] = tt[2]
          offset += 1
          coordinates[offset] = tt[3]
          offset += 1
          coordinates[offset] = tt[4]
          offset += 1
        model.setSubmatrixData('coordinates', coordinates)
      
      inputModelData = None
      parameterGroup = None
      # generate output data as ConstraintLists
      newConstraintLists = {}
      for label, name, code, numResonances in interfaceData:
        data = modelDict.get(code)
        if data is not None:
          
          if parameterGroup is None:
            parameterGroup = nmrCalcRun.newParameterGroup(name=newData.name,
                                                          ioRole='output',
                                                          data=(newData,))
          
          if constraintStore is None:
            # Initialise ConstraintStore
            constraintStore = nmrCalcRun.root.newNmrConstraintStore(
                               nmrProject=nmrCalcRun.nmrCalcStore.nmrProject)
          
          if inputModelData is None:
            # Make input model record, to group of input and utput models
            inputModel = inputEnsembleData.structureEnsemble.findFirstModel(
                           serial=modelSerial)
            inputModelData = nmrCalcRun.newStructureEnsembleData(ioRole='input',
                              models=(inputModel,), 
                              parameterGroup=parameterGroup,
                              name=parameterGroup.name)
          
          constraintList = intUtil.constraintListFromData(nmrCalcRun,
                                                          constraintStore,  
                                                          data, name, code, 
                                                          numResonances, 
                                                          parameterGroup)
          if constraintList:
            ll = newConstraintLists.get(name)
            if ll is None:
              ll = newConstraintLists[name] = []
            ll.append(constraintList)
    
    # NBNB REWORK! in standard operation you should generate only ONE 
    # measurementlist from ALL models (average with std dev.)
    
    for name,constraintLists in newConstraintLists:
      intUtil.measurementListFromConstraints(nmrCalcRun, name, constraintLists)
    
    
    
    nmrCalcRun.status = 'completed'
    
  else:
    nmrCalcRun.status = 'failed'

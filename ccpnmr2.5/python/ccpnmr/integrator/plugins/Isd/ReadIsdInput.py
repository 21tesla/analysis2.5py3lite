
"""
======================COPYRIGHT/LICENSE START==========================

ReadIsdInput.py: code for CCPN data model and code generation framework

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

import sys, os

from memops.universal import Io as uniIo
from memops.general import Io as genIo

from ccpnmr.analysis.core import AssignmentBasic
from ccpnmr.analysis.core.MoleculeBasic import getResidueMapping
from ccp.lib import StructureIo

from ccp.format.talos.projectIO import TalosSequenceFile, TalosChemShiftFile

from ccp.general import Constants as genConstants
from ccp.lib.MoleculeModify import makeMolecule

atomNameMap = {
 'HN':'H',
}

if __name__ == '__main__':
  """ Hack for rough creation, given that neither comaptibility nor PDB load works
  """
  
  pdbFile = '6628.pdb'
  talosFile = '6628.tab'
  molName = 'myMol'
  chainCode = 'A'
  workingDir = os.getcwd() 
  
  project = genIo.newProject('IsdInputTest')
  nmrProject = project.newNmrProject(name=molName)
  shiftList = nmrProject.newShiftList(name=molName)
  project.newAnalysisProject(name='dummy', nmrProject=nmrProject)
  
  # read Talos project
  path = uniIo.joinPath(workingDir, talosFile)
  if os.path.isfile(path):
    
    # get sequence data
    seqReader = TalosSequenceFile(path)
    seqReader.read()
    seqData = [(x.code1Letter, x.seqCode, x.seqInsertCode) 
               for x in seqReader.sequences[-1].elements]
    del seqReader
    
    typeMapDict = genConstants.code1LetterToCcpCodeDict['protein']
    ccpCodes = [typeMapDict.get(x[0].upper()) for x in seqData]
    
    if None in ccpCodes:
      ind = ccpCodes.index(None)
      raise Exception("No residue found for code %s in position %s" 
                      % (seqData[ind][0], ind))
    
    # make molecule and update with read-in data
    molecule = makeMolecule(project, 'protein', ccpCodes, molName=molName)
    residues = molecule.sortedMolResidues()
    for ii,tt in enumerate(seqData):
      res = residues[ii]
      res.seqCode = tt[1]
      res.seqInsertCode = tt[2]
      if tt[0] == 'c':
        # special case, disulfide linked CYS
        # change descriptor to disulfide linked
        # depends on knowledge of descriptor, so a bit hacky
        ss = res.descriptor
        if 'prot:H3,HG' in ss:
          ss = ss.replace('prot:H3,HG','prot:H3')
        elif 'prot:HG' in ss:
          ss = ss.replace('prot:HG','')
        res.descriptor = ';'.join((ss,'link:SG'))
    
    # make MolSystem
    molSystem = project.newMolSystem(name=molName, code=molName)
    chain = molSystem.newChain(code=chainCode, molecule=molecule)
    
    # set up standard mappings. Code based on Analysis initMolSystemChain
    atomSetMappings = []
    getMapping = getResidueMapping
    for residue in chain.sortedResidues():
      msg = "Making Atom Sets and Mappings for residue %s %s %d"
      print msg % (chain.code,residue.ccpCode,residue.seqCode)
      residueMapping = getMapping(residue, aromaticsEquivalent=True)
      atomSetMappings.extend( residueMapping.atomSetMappings )
      
    atomSetDict = {}
    for atomSet in nmrProject.atomSets:
      atomSetDict[atomSet.serial] = atomSet

    for atomSetMapping in atomSetMappings:
      serials  = atomSetMapping.atomSetSerials
      atomSets = []
      for serial in serials:
        atomSets.append(atomSetDict[serial])

      AssignmentBasic.updateAtomSetMapping(atomSetMapping,atomSets=atomSets)
   
    # get shift data
    shiftReader = TalosChemShiftFile(path)
    shiftReader.read()
    shiftData = [
     (x.seqCode, x.seqInsertCode, x.resLabel, x.atomName, x.value)
     for x in shiftReader.chemShifts
    ]
    del shiftReader
    
    for tt in shiftData:
      seqCode, seqInsertCode, resLabel, atomName, value = tt
      res = chain.findFirstResidue(seqCode=seqCode, seqInsertCode=seqInsertCode)
      atName = atomNameMap.get(atomName,atomName)
      atom = res.findFirstAtom(name=atName)
      if atom is None:
        raise Exception("No atom found for shift: %s" % str(tt))
        #print "WARNING, no atom found for shift: %s" % str(tt)
      
      if atName.startswith('H'):
        isotope = '1H'
      elif atName.startswith('C'):
        isotope = '13C'
      elif atName.startswith('N'):
        isotope = '15N'
      else:
        raise Exception("No isotope determined for atom %s" % str(tt))
      
      resonance = nmrProject.newResonance(isotopeCode=isotope)
      shiftList.newShift(resonance=resonance, value=value)
      if atom is not None:
        AssignmentBasic.assignAtomsToRes((atom.atomSet,),resonance)
        if atName in ('HA2','HA3'):
          # Might not be the most efficient way, but should hopefully work
          AssignmentBasic.swapProchiralResonance(resonance, makeAmbiguous=True)
    
  else:
    raise Exception("Talos file not found: %s" % path)
  
  
  # read pdb structures
  path = uniIo.joinPath(workingDir, pdbFile)
  if os.path.isfile(path):
    ensemble = StructureIo.getStructureFromFile(molSystem, path)
    
  else:
    raise Exception("PDB file not found: %s" % path)
  
  
  project.saveModified()

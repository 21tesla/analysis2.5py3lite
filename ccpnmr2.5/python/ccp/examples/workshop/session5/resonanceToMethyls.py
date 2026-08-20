#
# Exercise 1.2: connect Nmr.Resonance(s) in the correct way to the following atom(s):
#
# - Ala 3  HB* (DO THIS AS EXAMPLE!)
# - Leu 17 HD1* and HD2*, ambiguous
#

from memops.api import Implementation
from ccp.examples.workshop.session2.makeMolSystem import makeSequence, makeMolecule
from ccp.examples.workshop.session5.resonanceToSingleAtoms import linkResonanceAtomsUnambiguous, printResonanceToAtoms

if __name__ == '__main__':

  ccpnProject = Implementation.MemopsRoot(name = 'example1b')
  
  #
  # Make a Molecule, MolSystem and Chain to work with
  #
  
  sequence = 'WHATSTHISCCPNAPILIKE'
  sequence = makeSequence(sequence)
  molecule = makeMolecule(ccpnProject, sequence, molName='CcpnMol')
  molSystem = ccpnProject.newMolSystem(code='testMolSystem')
  chain = molSystem.newChain(code = 'A', molecule = molecule)

  #
  # Make nmrProject
  #

  nmrProject = ccpnProject.newNmrProject(name = 'test')

  #
  # Ala 3  HB*
  #
  
  #
  # Find the right residue
  #
  
  residue = chain.findFirstResidue(seqId = 3)
  
  #
  # Now have to find 3 atoms for the HB*. This is best done via the ChemComp reference data...
  #
  
  chemCompVar = residue.chemCompVar
  chemAtomSet = chemCompVar.findFirstChemAtomSet(name = 'HB*')
  
  #
  # Now have the chemAtomSet reference data, find the right atoms...
  #
  
  atomList = []
  
  for chemAtom in chemAtomSet.chemAtoms:
   
    atom = residue.findFirstAtom(name = chemAtom.name)
    atomList.append(atom)
  
  #
  # Then create a direct link between a single resonance and the atomSet as before
  #
  
  resonance = nmrProject.newResonance(isotopeCode = '1H')
    
  linkResonanceAtomsUnambiguous(chemAtomSet.name,atomList,resonance)
  
  printResonanceToAtoms(resonance)

  
  
  #
  # Leu 17 HD1* and HD2*, ambiguous
  #

  residue = chain.findFirstResidue(seqId = 17)
    
  chemCompVar = residue.chemCompVar
  
  atomSetList = []
  resonanceList = []
  
  for chemAtomSetName in ('HD1*','HD2*'):
    chemAtomSet = chemCompVar.findFirstChemAtomSet(name = chemAtomSetName)
  
    atomList = []
    
    for chemAtom in chemAtomSet.chemAtoms:

      atom = residue.findFirstAtom(name = chemAtom.name)
      atomList.append(atom)
      
    atomSetList.append(nmrProject.newAtomSet(name = chemAtomSet.name,atoms = atomList))
    resonanceList.append(nmrProject.newResonance(isotopeCode = '1H'))

  nmrProject.newResonanceSet(resonances = resonanceList, atomSets = atomSetList)

  printResonanceToAtoms(resonanceList[0])
  printResonanceToAtoms(resonanceList[1])

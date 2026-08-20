#
# Exercise 1.1: connect Nmr.Resonance(s) in the correct way to the following atom(s):
#
# - Ile 16 HB  (DO THIS AS EXAMPLE!)
# - Trp 1  CZ3
# - Ser 5  HB2 and HB3 stereo assigned
# - Ser 9  HB2 and HB3 ambiguous
#

from ccp.examples.workshop.session2.makeMolSystem import makeSequence, makeMolecule
from memops.api import Implementation

def linkResonanceAtomsUnambiguous(atomSetName,atomList,resonance):

  nmrProject = resonance.parent
  atomSet   = nmrProject.newAtomSet(name = atomSetName, atoms = atomList)
  resSet    = nmrProject.newResonanceSet(resonances = [resonance], atomSets = [atomSet])
  
def printResonanceToAtoms(resonance):
  
  print 'Resonance %d' % resonance.serial,
  
  resSet = resonance.resonanceSet
  
  if resSet:
    print ' -> ResonanceSet ->'
    
    for atomSet in resSet.atomSets:
      refAtom = atomSet.findFirstAtom()
      residue = refAtom.residue
      print "        AtomSet '%s' with atoms %s (chain '%s', residue %d)" % (atomSet.name,str([atom.name for atom in atomSet.sortedAtoms()]),residue.chain.code,residue.seqId)
  
  else:
    print
    
  print

if __name__ == '__main__':

  ccpnProject = Implementation.MemopsRoot(name = 'example1a')
  
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
  # Ile 16 HB
  #
  
  #
  # Find the right atom via the residue (should be a function)
  #
  
  residue = chain.findFirstResidue(seqId = 16)   # Alternatively seqCode, but less dependable
  atom =    residue.findFirstAtom(name = 'HB')

  #
  # Create a new resonance, isotopeCode has to be set, name is optional (but useful)
  #
  
  resonance = nmrProject.newResonance(isotopeCode = '1H', name = 'Ile.15.HB')
  resonance.details = 'Test resonance'
  
  #
  # Set the link from resonance to atom via atomSet and resonanceSet (should be a function)
  # 
  # Note that atomSet.name is not obligatory, but very useful nevertheless.
  #
  
  atomSet   = nmrProject.newAtomSet(name = atom.name, atoms = [atom])
  resSet    = nmrProject.newResonanceSet(resonances = [resonance], atomSets = [atomSet])
  
  printResonanceToAtoms(resonance)
  
  
  #
  # Trp 1 CZ3
  #

  residue = chain.findFirstResidue(seqId = 1)
  atom =    residue.findFirstAtom(name = 'CZ3')

  resonance = nmrProject.newResonance(isotopeCode = '13C')  
  linkResonanceAtomsUnambiguous(atom.name,[atom],resonance)
  
  printResonanceToAtoms(resonance)


  
  #
  # Ser 5  HB2 and HB3 stereo assigned
  #
  
  residue = chain.findFirstResidue(seqId = 5)
  
  for atomName in ('HB2','HB3'):
    
    atom = residue.findFirstAtom(name = atomName)
    resonance = nmrProject.newResonance(isotopeCode = '1H')
    
    linkResonanceAtomsUnambiguous(atom.name,[atom],resonance)
   
    printResonanceToAtoms(resonance)
   
    
  #
  # Ser 9  HB2 and HB3 ambiguous - can be done more elegantly!
  #
  
  residue = chain.findFirstResidue(seqId = 9)
  atom1 = residue.findFirstAtom(name = 'HB2')
  atom2 = residue.findFirstAtom(name = 'HB3')
  
  resonance1 = nmrProject.newResonance(isotopeCode = '1H')
  resonance2 = nmrProject.newResonance(isotopeCode = '1H')
  
  atomSet1   = nmrProject.newAtomSet(name = atom1.name, atoms = [atom1])
  atomSet2   = nmrProject.newAtomSet(name = atom2.name, atoms = [atom2])
  
  resSet    = nmrProject.newResonanceSet(resonances = [resonance1,resonance2], atomSets = [atomSet1,atomSet2])

  printResonanceToAtoms(resonance1)
  printResonanceToAtoms(resonance2)

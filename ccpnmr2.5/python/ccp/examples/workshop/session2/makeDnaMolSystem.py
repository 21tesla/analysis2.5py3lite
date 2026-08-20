# Exercise 2.2: make a DNA molSystem from a DNA molecule.
#
# - Load the project made in exercise 2.1.
#
# - Create a new molSystem with a chain attached to the molecule from
#      the 'reference' molecule project.
#
# - Print out the residue type, position and all the atoms found in
#     each residue.
#
# - Check that the project is valid.
#
# - Hint: Look at class methods for MemopsRoot objects.
#

from memops.general.Io import loadProject

if __name__ == '__main__':

  # Load the project with the DNA molecule.

  project = loadProject('dnaTest')

  # Locate the DNA molecule in the project.

  molecule = project.findFirstMolecule(name = 'myDnaMoleculeName')

  # Make a new molSystem top object.

  keywds = {'code': 'myMolecularSystemCode',
            'name': 'My molecular system name'}
  molSystem = project.newMolSystem(**keywds)

  # Make a chain for this molSystem that is linked to the reference
  #   DNA molecule.

  chain = molSystem.newChain(code = 'A', molecule = molecule)

  # Print some basic information about the chain and molecule.

  print 'Code [%s], name [%s]' % (chain.code, chain.molecule.name)

  # Navigate residues in this chain and print atom names.

  for res in chain.sortedResidues():
    print 'Residue [%s], position [%s]' % (res.ccpCode, res.seqCode + 1)
    print '  Atoms: ',

    for atom in res.sortedAtoms():
      print atom.name,
    print '\n'

  # Hack to find the currentElementStore.

  project.currentChemElementStore = project.findFirstChemElementStore()

  # Check that the project is valid.

  project.checkAllValid(complete = True)

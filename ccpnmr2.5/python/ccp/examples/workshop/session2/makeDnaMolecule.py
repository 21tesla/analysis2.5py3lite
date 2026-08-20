# Exercise 2.1: make a new DNA reference molecule.
#
# - Create a molecule of DNA with the sequence 'ATGGATCATTAG'.
#
# - Hint: use the example of making a protein given in the notes in ex4.py .
#
# - Molecules are top objects in the Molecule package.
#
# - From this top object, molResidues can be made, but it is first necessary
#     to get the reference ChemComp information for each residue type.
#
# - For each residue:
#
#   i) Find the first ChemComp with the relevant one-letter DNA code and
#        make the 'molType' to be 'DNA'.
#
#   ii) Then set the correct linking ('start' for the first residue,
#         'middle', or 'end' for the last residue).
#
#   iii) Use this linking to find the first ChemCompVar with 'isDefaultVar'
#          set to 'True'.
#
#   iv) Next, make a new MolResidue with the ChemComp, linking, residue
#         position and descriptor information (from the ChemCompVar).
#
#   v) Finally, link the residues to each other - see hint below.
#
# - Check that the molecule has been made into a correct, linear molecule
#     using molecule.isStdLinear .
#
# - Also, navigate the molecule object for its name and the molResidue codes.
#
# - Save this molecule in a new CCPN project to XML files.
#
# - Hint: we have provided a sub-routine to do the linking of residues
#     to each other:
#
#   - from ccp.examples.workshop.session2.makeMolSystem import makeLinearLink
#
#   - makeLinearLink(molecule, molResidue, position, linking) .
#

from memops.api.Implementation import MemopsRoot

# makeLinearLink sub-routine to make the correct links between the molResidues.

from ccp.examples.workshop.session2.makeMolSystem import makeLinearLink

if __name__ == '__main__':

  project = MemopsRoot(name = 'dnaTest')

  # Make a new molecule top object.

  molecule =  project.newMolecule(name = 'myDnaMoleculeName')

  # DNA sequence as a string of one-letter codes.

  seq = 'ATGGATCATTAG'

  # Iterate over this sequence.

  for pos in range(0, len(seq) ):

    # Get the residue code.

    ccpCode = seq[pos]

    # Find the DNA ChemComp for this residue code.

    keywds = {'molType': 'DNA',
              'ccpCode': ccpCode}
    chemComp = project.findFirstChemComp(**keywds)

    # Specify the correct linking information for this residue
    # based on its position in the sequence.

    if pos == 0:
      linking = 'start'
    elif pos == len(seq) - 1:
      linking = 'end'
    else:
      linking = 'middle'

    # Find the ChemCompVar (ChemComp variant) based on the
    # linking information.

    keywds2 = {'linking':      linking,
               'isDefaultVar': True}
    chemCompVar = chemComp.findFirstChemCompVar(**keywds2)

    # Make the molResidue from the parent molecule object -
    # and linked to the correct ChemComp.

    keywds3 = {'chemComp':   chemComp,
               'seqCode':    pos,
               'linking':    linking,
               'descriptor': chemCompVar.descriptor}
    molRes = molecule.newMolResidue(**keywds3)

    # Create a link between two consecutive molResidues.

    makeLinearLink(molecule, molRes, pos, linking)

  print 'Linear: [%s]' % molecule.isStdLinear, '\n'

  print 'Molecule: [%s]' % molecule.name, '\n'

  for molRes in molecule.sortedMolResidues():
    print 'MolResidue: [%s]' % molRes.ccpCode

  print

  project.saveModified()

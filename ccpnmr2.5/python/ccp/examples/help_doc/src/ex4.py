from memops.api.Implementation import MemopsRoot

if __name__ == '__main__':

    # CCPN project.
    project = MemopsRoot(name = 'moleculeTest')

    # Create a molecule. This is a reference object - molSystem contains the
    # actual molecule. e.g., in a homodimer, you only create 1 molecule, but
    # the molSystem has 2 chains (see below).
    molecule =  project.newMolecule(name = 'myMoleculeName')

    # Then we have to create the MolResidues. All the reference information
    # for these comes from the 'ChemComps', which must be loaded first.

    # Here is an example protein sequence.
    seq = ['Ala', 'Gly', 'Tyr', 'Glu', 'Leu', 'Gly', 'Ser', 'His', 'Ile']

    for pos in range(0, len(seq) ):
        ccpCode = seq[pos] 

        # Find the reference data chemComp. Note the type of residue needs to be
        # specified from: 'protein', 'DNA', 'RNA' 'carbohydrate', and 'other'.
        keywds = {'molType': 'protein',
                  'ccpCode': ccpCode}
        chemComp = project.findFirstChemComp(**keywds)

        # For linear biopolymers ('protein', 'DNA' and 'RNA) the 'linking'
        # needs to be set, which indicates the position in the linear chain.
        # For other types ('carbohydrate', 'other') linking indicates the
        # pattern of links to other building blocks.
        if pos == 0:
            linking = 'start'
        elif pos == len(seq) - 1:
            linking = 'end'
        else:
            linking = 'middle'

        # Additionally the 'descriptor', which defines the protonation state 
        # of the residue (the variant of the ChemComp: a ChemCompVar), needs 
        # to be found. The default state is selected using the findFirst
        # function in the API.
        keywds2 = {'linking':      linking,
                   'isDefaultVar': True}
        chemCompVar = chemComp.findFirstChemCompVar(**keywds2)

        # Create a molResidue.
        keywds3 = {'chemComp':   chemComp,
                   'seqCode':    pos,
                   'linking':    linking,
                   'descriptor': chemCompVar.descriptor}
        molRes = molecule.newMolResidue(**keywds3)

        # Finally the residue linking information is set.
        if linking in ['middle', 'end']:
            prevLink = molRes.findFirstMolResLinkEnd(linkCode = 'prev')
            prevMolRes = molecule.findFirstMolResidue(serial = pos)

            if prevLink and prevMolRes:
                nextLink = prevMolRes.findFirstMolResLinkEnd(linkCode = 'next')

                if nextLink:
                    mrle = [nextLink, prevLink]
                    molResLink = molecule.newMolResLink(molResLinkEnds = mrle)

    # Now that the Molecule information is set, creating a molecular system
    # with chains is straightforward. All the residue and atom information for
    # these are taken from the reference Molecule and ChemComp information.
    keywds4 = {'code': 'myMolecularSystemCode',
               'name': 'My molecular system name'}
    molSystem = project.newMolSystem(**keywds4)

    # Now create a homodimer - two chains referring to the same Molecule.
    chainA = molSystem.newChain(code = 'A', molecule = molecule)
    chainB = molSystem.newChain(code = 'B', molecule = molecule)

    # Now we can navigate the system. All the reference info on bonds, angles,
    # etc can also be obtained. Refer to the CCPN data model documentation to
    # see which attributes and links you can access from here.
    print 'Code [%s], name [%s]' % (chainA.code, chainA.molecule.name)
    print 'Code [%s], name [%s]' % (chainB.code, chainB.molecule.name)

    for res in chainA.sortedResidues():
        print 'Residue [%s], position [%s]' % (res.ccpCode, res.seqCode + 1)
        print '  Atoms: ',

        for atom in res.sortedAtoms():
            print atom.name,
        print '\n'

    # Check if whole project is valid using a class method. The first line
    # is a hack from earlier APIs, which should not be necessary eventually,
    # when current<Class> objects will be auto-created.
    project.currentChemElementStore = project.findFirstChemElementStore()
    project.checkAllValid(complete = True)

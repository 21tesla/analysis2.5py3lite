from memops.api.Implementation import MemopsRoot

def makeSequence(sequence):

  seqDict = {
    'A': 'Ala',
    'C': 'Cys',
    'D': 'Asp',
    'E': 'Glu',
    'F': 'Phe',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'K': 'Lys',
    'L': 'Leu',
    'M': 'Met',
    'N': 'Asn',
    'P': 'Pro',
    'Q': 'Gln',
    'R': 'Arg',
    'S': 'Ser',
    'T': 'Thr',
    'V': 'Val',
    'W': 'Trp',
    'Y': 'Tyr'}

  tlcSeq = []

  for res in sequence:
    tlcSeq.append(seqDict[res])

  return tlcSeq


def makeLinearLink(molecule, molRes, pos, linking):

  if linking in ['middle', 'end']:
    prevLink = molRes.findFirstMolResLinkEnd(linkCode = 'prev')
    prevMolRes = molecule.findFirstMolResidue(serial = pos)

    if prevLink and prevMolRes:
      nextLink = prevMolRes.findFirstMolResLinkEnd(linkCode = 'next')

      if nextLink:
        mrle = [nextLink, prevLink]
        molResLink = molecule.newMolResLink(molResLinkEnds = mrle)


def makeBranchedLink(molecule, molRes, linkCode, boundMolRes, boundLinkCode):

  prevLink = molRes.findFirstMolResLinkEnd(linkCode = linkCode)
  
  if prevLink:
    nextLink = boundMolRes.findFirstMolResLinkEnd(linkCode = boundLinkCode)
    
    if nextLink:
      molResLink = molecule.newMolResLink(molResLinkEnds = [nextLink, prevLink])


def makeMolecule(project, sequence, molType = 'protein', molName = 'molName'):

  # Create a molecule. This is a reference object - molSystem contains the
  # actual molecule. e.g., in a homodimer, you only create 1 molecule, but
  # the molSystem has 2 chains (see below).
  molecule = project.newMolecule(name = molName)

  # Then we have to create the MolResidues. All the reference information
  # for these comes from the 'ChemComps', which must be loaded first.

  for pos, ccpCode in enumerate(sequence):

    # Find the reference data chemComp. Note the type of residue needs to be
    # specified from: 'protein', 'DNA', 'RNA' 'carbohydrate', and 'other'.
    keywds = {'molType': molType,
              'ccpCode': ccpCode}
    chemComp = project.findFirstChemComp(**keywds)

    # For linear biopolymers ('protein', 'DNA' and 'RNA) the 'linking'
    # needs to be set, which indicates the position in the linear chain.
    # For other types ('carbohydrate', 'other') linking indicates the
    # pattern of links to other building blocks.
    if pos == 0:
      linking = 'start'
    elif pos == len(sequence) - 1:
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

    # Debug print.
    #print pos, ccpCode, chemComp, linking, chemCompVar

    # Create a molResidue.
    keywds3 = {'chemComp':   chemComp,
               'seqCode':    pos,
               'linking':    linking,
               'descriptor': chemCompVar.descriptor}
    molRes = molecule.newMolResidue(**keywds3)

    # Finally the residue linking information is set.
    makeLinearLink(molecule, molRes, pos, linking)

  return molecule


def makeMolSystem(project, molecule, molSysName = 'molSysName',
                  molSysCode = 'molSysCode', chainCodes = 'A'):

  # Now that the Molecule information is set, creating a molecular system
  # with chains is straightforward. All the residue and atom information for
  # these are taken from the reference Molecule and ChemComp information.
  keywds = {'code': molSysCode,
            'name': molSysName}
  molSystem = project.newMolSystem(**keywds)

  # Now create the chains.
  for code in chainCodes:
    molSystem.newChain(code = code, molecule = molecule)

  return molSystem


if __name__ == '__main__':

  # CCPN project.
  project = MemopsRoot(name = 'moleculeTest')

  #           12345678901234567890
  sequence = 'ACDEFGHIKLMNPQRSTVWY'

  # Make a three-letter sequence list from a one-letter sequence string.
  tlcSeq = makeSequence(sequence)

  molecule = makeMolecule(project, tlcSeq)

  print 'Molecule: [%s]' % molecule.name, '\n'

  for molRes in molecule.sortedMolResidues():
    print 'MolResidue: [%s]' % molRes.ccpCode

  print

  molSystem = makeMolSystem(project, molecule, chainCodes = ['A', 'B'])

  print 'MolSystem: [%s]' % molSystem.name, '\n'

  for chain in molSystem.sortedChains():
    print 'Chain: [%s]' % chain.code

  print

  dnaSeq = 'ACGT'

  dnaMol = makeMolecule(project, dnaSeq, molType = 'DNA', molName = 'dnaMol')

  print 'Molecule: [%s]' % dnaMol.name, '\n'

  for molRes in dnaMol.sortedMolResidues():
    print 'MolResidue: [%s]' % molRes.ccpCode

  print

  dnaMolSys = makeMolSystem(project, dnaMol, chainCodes = 'A', molSysCode = 'dnaMolCode')

  print 'MolSystem: [%s]' % dnaMolSys.name, '\n'

  for chain in dnaMolSys.sortedChains():
    print 'Chain: [%s]' % chain.code

  print

  project.currentChemElementStore = project.findFirstChemElementStore()
  project.checkAllValid(complete = True)

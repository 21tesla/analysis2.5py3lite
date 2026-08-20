# Exercise 2.4: make a carbohydrate molecule and molSystem.
#
# - Make a carbohydrate from two sugar residues - 'molType' is 'carbohydrate'.
#
# - The two residues are in a local data directory: '../data'.
#
# - You will need to add a new repository with the above url to the
#     packageLocator with the 'targetName' of 'ccp.molecule.ChemComp'.
#
# - Note that you don't need to make ChemCompVars like with DNA or proteins.
#
# - Also, the linking string is more complicated for the branched carbohydrates.
#
# - The first residue is called 'dglc-hex-1-5:C1_OMe:C2_NAc' and has the
#     following attributes 'linking': 'link:O4' and 'descriptor': 'stereo_2:C1'.
#
# - The second residue is called 'dglc-hex-1-5:C2_NAc' and has the
#     following attributes 'linking': 'link:C1_2' and 'descriptor': 'stereo_2:C1'.
#
# - Hint: we have provided a sub-routine to do the linking of residues
#     to each other:
#
#   - from ccp.examples.workshop.session2.makeMolSystem import makeBranchedLink
#
#   - makeBranchedLink(molecule, molResidue, linkCode,
#                     boundMolResidue, boundLinkCode) .
#     where: 'linkCode': 'C1_2' and 'boundLinkCode': 'O4'
#
# - Finally, navigate the molecule, molSystem, chain and linked atom objects.
#

import os

from memops.api import Implementation
from memops.universal import Io as uniIo

# makeBranchedLink sub-routine to make the linking between two
#   carbohydrate residues.

from ccp.examples.workshop.session2.makeMolSystem import makeBranchedLink

# Sub routine to set the package location of the local 
#   carbohydrate ChemComp XML files.

def resetPackageLocators(project):

  carbDataPath = uniIo.normalisePath(os.path.join(os.path.abspath('..'), 
                                     'data'))

  print 'Location of local carbohydrate ChemComps: [%s]' % carbDataPath

  carbDataUrl = Implementation.Url(path = carbDataPath)

  carbDataRepos = project.newRepository(name = 'carbData',
                                        url = carbDataUrl)

  chemPackLoc = project.findFirstPackageLocator(targetName = 'ccp.molecule.ChemComp')

  chemPackLoc.addRepository(carbDataRepos)


def makeCarbMolecule(project, molName = 'carbMol'):

  resetPackageLocators(project)

  # Find the two carbohydrate ChemComps.

  keywds = {'molType': 'carbohydrate',
            'ccpCode': 'dglc-hex-1-5:C1_OMe:C2_NAc'}
  glcNAcOMe = project.findFirstChemComp(**keywds)

  keywds2 = {'molType': 'carbohydrate',
             'ccpCode': 'dglc-hex-1-5:C2_NAc'}
  glcNAc = project.findFirstChemComp(**keywds2)

  # Make a molecule.

  molecule = project.newMolecule(name = molName)

  resNum = 1

  # Add two molResidues to this molecule linked to the two earlier ChemComps.

  keywds3 = {'seqCode':    resNum,
             'chemComp':   glcNAcOMe,
             'linking':    'link:O4',
             'descriptor': 'stereo_2:C1'}
  molRes1 = molecule.newMolResidue(**keywds3)

  resNum += 1

  keywds4 = {'seqCode':    resNum,
             'chemComp':   glcNAc,
             'linking':    'link:C1_2',
             'descriptor': 'stereo_2:C1'}
  molRes2 = molecule.newMolResidue(**keywds4)

  # Make the correct link between these two molResidues.

  boundMolRes = molecule.findFirstMolResidue(serial = 1)

  makeBranchedLink(molecule, molRes2, 'C1_2', boundMolRes, 'O4')

  return molecule


def makeCarbMolSystem(project, molecule,  molSysName = 'carbMolSys',
                      molSysCode = 'carbCode', chainCodes = 'A'):

  # Create a molecular system from this 'reference' molecule.

  keywds = {'code': molSysCode,
            'name': molSysName}
  molSystem = project.newMolSystem(**keywds)

  # Add a chain to this molSystem attached to the 'reference' molecule.

  for code in chainCodes:
    molSystem.newChain(code = code, molecule = molecule)

  return molSystem


if __name__ == '__main__':

  project = Implementation.MemopsRoot(name = 'carbTest')

  molecule = makeCarbMolecule(project)

  molSystem = makeCarbMolSystem(project, molecule)

  # Navigate the data model - molecules first.

  print "\nMolecules:"

  for mol in project.sortedMolecules():
    print "  %s" % mol.name

  print

  # MolSystem.

  print "Molecular system: [%s] [%s]\n" % (molSystem.code, molSystem.name)

  # Chains, residues and linking info.

  for chain in molSystem.sortedChains():

    print "Chain: [%s]\n" % chain.code
    molecule = chain.molecule

    for res in chain.sortedResidues():

      print "Residue: [%s.%s]" % (res.ccpCode, res.seqCode)
      molRes = res.molResidue

      for mrle in molRes.molResLinkEnds:
        molResLink = mrle.molResLink
        molResLink_linkEnds = list(molResLink.molResLinkEnds)

        otherMrle = molResLink_linkEnds[not molResLink_linkEnds.index(mrle)]
        otherMolRes = otherMrle.parent
        otherRes = chain.findFirstResidue(molResidue = otherMolRes)

        print "  Atom: [%s] linked to Residue: [%s.%s] Atom: [%s]" % (
          mrle.linkCode, otherRes.ccpCode, otherRes.seqCode, otherMrle.linkCode)

      print

  project.saveModified()

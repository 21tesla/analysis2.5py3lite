"""
##################
# Example script #
##################

Shows how to create a molecule and molecular system from scratch.

Level: intermediate

Contact: Wim Vranken <wim@ebi.ac.uk>
"""

#
# Import the Implementation package - this is the root package
#

from memops.api import Implementation

#
# Import the Molecule and MolSystem packages
#

from ccp.api.molecule import Molecule, MolSystem

#
# Python standard stuff
#

import os, shutil

if __name__ == '__main__':

  currentDir = os.path.abspath('.')
  projectDir = os.path.join(currentDir,'local')
  projectName = 'test'

  #
  # Make sure the projectDir exists, and go to that directory...
  #
  
  if not os.path.exists(projectDir):
    os.mkdir(projectDir)
  
  os.chdir(projectDir)
  
  #
  # Remove if already exists
  #
  
  if os.path.exists(projectName):
    shutil.rmtree(projectName)

  #
  # First you always have to create a project (this is the root object)
  #
  
  project = Implementation.MemopsRoot(name = projectName)

  #
  # Now we can create a molecule. This is a reference object: if you have
  # a homodimer, for example, you only create 1 molecule, but (see below) 2 chains.
  #
  # Note that you have to pass in the parent object (the project in this case).
  # The child object (here the Molecule) can never exist without its parent
  #
  
  molecule =  Molecule.Molecule(project,name = 'myMoleculeName')
  
  #
  # Then we have to create the MolResidues. All the reference information
  # for these comes from the 'ChemComps', and you have to load this ChemComp
  # before you can create the MolResidue.
  #
  
  exampleSequence = ['Ala','Gly','Tyr','Glu','Leu','Gly','Ser','His','Ile']
  
  for seqPosition in range(0,len(exampleSequence)):
  
    ccpCode = exampleSequence[seqPosition] 
    
    #
    # Find the reference data chemComp. Note that you also have to
    # specify the type of residue - available are 'protein', 'DNA', 'RNA' 'carbohydrate', and 'other'
    #
    
    chemComp = project.findFirstChemComp(molType = 'protein', ccpCode = ccpCode)
    
    #
    # For linear biopolymers ('protein','DNA' and 'RNA) you have to set the 'linking'. This
    # indicates the position in the linear chain
    # for other types ('carbohydrate', 'other') linking indicates the pattern
    # of links to other building blocks
    #
    
    if seqPosition == 0:
      linking = 'start'
    elif seqPosition == len(exampleSequence) - 1:
      linking = 'end'
    else:
      linking = 'middle'
      
    #
    # You also have to specify the 'descriptor' which defines the protonation
    # state of the residue (the variant of the ChemComp: a ChemCompVar).
    # Here the default state is selected using the findFirst function provided by the API
    #
    
    chemCompVar = chemComp.findFirstChemCompVar(linking = linking, isDefaultVar = True)
        
    #
    # Create a molResidue. You can also set the 
    #

    molRes = Molecule.MolResidue(molecule, chemComp = chemComp,
                                           seqCode = seqPosition,
                                           linking = linking,
                                           descriptor = chemCompVar.descriptor)
    
    #
    # Finally we have to set the residue linking information...
    #

    if linking in ['middle','end']:

      prevLink = molRes.findFirstMolResLinkEnd(linkCode = 'prev')
      prevMolRes = molecule.findFirstMolResidue(serial = seqPosition)

      if prevLink and prevMolRes:

        nextLink = prevMolRes.findFirstMolResLinkEnd(linkCode = 'next')

        if nextLink:

          molResLink = Molecule.MolResLink(molecule,molResLinkEnds = [nextLink,prevLink])

  #
  # Now that the Molecule information is set, creation a molecular system with chains
  # is very straightforward. All the residue and atom information for these chains is
  # taken from the reference Molecule and ChemComp information and set automatically...
  #

  molSystem = MolSystem.MolSystem(project,code = 'myMolecularSystemCode', name = 'My molecular system name')
  
  #
  # Now create a homodimer - two chains referring to the same Molecule
  #
  
  chainA = MolSystem.Chain(molSystem,code = 'A', molecule = molecule)
  chainB = MolSystem.Chain(molSystem,code = 'B', molecule = molecule)
  
  #
  # Now we can navigate the system... you can also go back to all the reference information to
  # get info on bonds, angles, ... . Refer to the CCPN data model documentation to see which
  # attributes and links you can access from here. For example for residues, see:
  #
  # http://www.ccpn.ac.uk/ccpnmr/ccpnmr/python/ccp/api/doc/MolSystem/Residue/index.html
  #
  # (or even better the local html files that come with the distribution)
  #
  
  for residue in chainA.residues:
    print residue.ccpCode, residue.seqCode
    print "  ",
    for atom in residue.sortedAtoms():
      print atom.name,
    print
    
  print chainA.code, chainA.molecule.name
  print chainB.code, chainB.molecule.name

  #
  # Check if whole project is valid
  #
  
  # This line is a hack for early API - should not be necessary eventually
  project.currentChemElementStore = project.findFirstChemElementStore()
  
  project.checkAllValid(complete = True)

  #
  # Write out the project.
  #
  
  project.saveModified()
  

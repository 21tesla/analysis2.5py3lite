import random

from ccp.examples.workshop.session2.makeMolSystem import makeSequence, makeMolecule

def main(root):

  # make Molecule
  sequence = 'WHATSTHISCCPNAPILIKE'
  sequence = makeSequence(sequence)
  molecule = makeMolecule(root, sequence, molName='CcpnMol')
  
  # create MolSystem
  molSystem = root.newMolSystem(code='testMolSystem')

  # create MolSystem.Chains
  for code in ('A', 'B'):
    molSystem.newChain(code=code, molecule=molecule)

  # temporary bodge
  if not root.currentChemElementStore:
    root.currentChemElementStore = root.findFirstChemElementStore()

  # print some information out
  chain = molSystem.findFirstChain(code='A')
  print 'chain empiricalFormula:', chain.empiricalFormula
  print 'chain formalCharge:', chain.formalCharge
  print 'chain molecularMass:', chain.molecularMass
  print 'molSystem molecularMass:', molSystem.molecularMass
 
  print 'sequence:', [residue.ccpCode for residue in chain.sortedResidues()]

  # create StructureEnsemble
  structureEnsemble = root.newStructureEnsemble(ensembleId=23, molSystem=molSystem)

  # create MolStructure.Chains
  atoms = []
  for code in ('A', 'B'):
    chain = structureEnsemble.newChain(code=code)
    molSysChain = chain.chain

    # create MolStructure.Residues
    for seqId in range(3, 7):
      molSysResidue = molSysChain.findFirstResidue(seqId=seqId)
      residue = chain.newResidue(seqId=seqId, seqCode=molSysResidue.seqCode)

      # create MolStructure.Atoms
      for molSysAtom in molSysResidue.atoms:
        chemAtom = molSysAtom.chemAtom
        if chemAtom.elementSymbol != 'H':
          atom = residue.newAtom(name=molSysAtom.name)
          atoms.append(atom)

  # create Model
  model = structureEnsemble.newModel()
  
  # create MolStructure coordinates
  coordinates = [random.random() for dummy in range(len(atoms)*3)]
  model.setSubmatrixData('coordinates', coordinates)

  #root.saveModified()

if __name__ == '__main__':

  import sys

  if len(sys.argv) == 1:
    from memops.api.Implementation import MemopsRoot
    root = MemopsRoot(name='testNmr')
  elif len(sys.argv) == 2:
    from memops.general.Io import loadProject
    repositoryPath = sys.argv[1]
    root = loadProject(path=repositoryPath)

  main(root)


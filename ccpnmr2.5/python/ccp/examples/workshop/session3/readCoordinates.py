import glob
import os

from ccpnmr.format.converters.PdbFormat import PdbFormat
#from ccpnmr.format.converters.PseudoPdbFormat import PseudoPdbFormat

def main(root, directory = None):

  if not directory:
    directory = os.getcwd()

  fileNames = glob.glob(os.path.join(directory, '*.pdb'))
  print 'fileNames', fileNames

  # create PdbFormat
  pdbFormat = PdbFormat(root, verbose=1)
  #pdbFormat = PseudoPdbFormat(root, verbose=1)

  pdbFormat.readCoordinates(fileNames, minimalPrompts=1, linkAtoms=0)
  #pdbFormat.readCoordinates(fileNames, minimalPrompts=1, linkAtoms=0, forceNamingSystemName='PDB_REMED')

  #root.saveModified()

  print 'structureEnsembles', len(root.structureEnsembles)
  print 'molSystems %d' % len(root.molSystems)

  molSystem = root.findFirstMolSystem()
  print 'chains %d' % len(molSystem.chains)  

  for chain in molSystem.sortedChains():
    n = 0
    for residue in chain.residues:
      n += len(residue.atoms)

    print 'chain %s %d %d' % (chain.code, len(chain.residues), n)

  structureEnsemble = root.findFirstStructureEnsemble()
  print 'models %d' % len(structureEnsemble.models)

if __name__ == '__main__':

  import sys

  if len(sys.argv) == 1:
    from memops.api.Implementation import MemopsRoot
    root = MemopsRoot(name='testNmr')
  elif len(sys.argv) == 2:
    from memops.general.Io import loadProject
    repositoryPath = sys.argv[1]
    root = loadProject(path=repositoryPath)

  directory = '../data/pdb'

  main(root, directory=directory)


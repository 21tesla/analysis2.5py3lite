import os, glob
import Tkinter

# Another way to import the whole of the Implementation package.
from memops.api import Implementation

from ccpnmr.format.converters.PdbFormat import PdbFormat
from ccpnmr.format.converters.CnsFormat import CnsFormat

if __name__ == '__main__':

    # CCPN project.
    project = Implementation.MemopsRoot(name = 'cns2Pdb')

    guiRoot = Tkinter.Tk()

    # Create the format converter CNS object linked to the CCPN project.
    cnsObj = CnsFormat(project, guiRoot, verbose = 1)

    # Remember the current directory.
    cwd = os.getcwd()

    # Location of the CNS files.
    cnsDir = '../data/cns'
    os.chdir(cnsDir)

    cnsList = glob.glob('cns_*.pdb')
    os.chdir(cwd)

    # List with path of the CNS files.
    cnsFiles = []

    for cnsName in sorted(cnsList):
        cnsFile = os.path.join(cnsDir, cnsName)
        cnsFiles.append(cnsFile)

    # Need a structureGeneration object to hold the structures.
    # The parent object is an nmrProject.
    nmrProject = project.newNmrProject(name = project.name)
    strucGen   = nmrProject.newStructureGeneration()

    # Now read in the structures from the CNS files.
    # This will also call readSequence() to create molecule
    # and molSystem objects.
    cnsObj.readCoordinates(cnsFiles,
                           strucGen = strucGen,
                           minimalPrompts = 1,
                           linkAtoms = 0)

    # Create the format converter PDB output object linked to the CCPN project.
    pdbObj = PdbFormat(project, guiRoot, verbose = 1)

    curDir = os.path.abspath('../data')
    pdbDir = os.path.join(curDir, 'pdb')

    # Make sure the output directory exists.
    if not os.path.exists(pdbDir):
        os.mkdir(pdbDir)

    pdbName = 'cns_out.pdb'
    pdbFile = os.path.join(pdbDir, pdbName)

    strucList = list(strucGen.structureEnsemble.sortedModels() )

    # Write the PDB file and catch any exceptions raised.
    try:
        pdbObj.writeCoordinates(pdbFile,
                                structures = strucList,
                                minimalPrompts = 1)
    except Exception, exc:
        print "Can't write coordinates [%s]" % exc

    # Save the CCPN project as XML files using the method saveModified().
    project.saveModified()

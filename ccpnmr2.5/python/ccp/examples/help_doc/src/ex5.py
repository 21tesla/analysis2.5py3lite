# Standard python library.
import os

# Get Tkinter for popups.
import Tkinter

# Read the whole of the Implementation package.
import memops.api.Implementation as Implementation

# Get FastaFormat class for format conversion.
from ccpnmr.format.converters.FastaFormat import FastaFormat

if __name__ == '__main__':

    # CCPN project.
    project = Implementation.MemopsRoot(name = 'readFasta')

    # Create a Tkinter object.
    guiRoot = Tkinter.Tk()

    # Create a FastaFormat object and associate it with the CCPN project and
    # Tkinter object.
    fastaObj = FastaFormat(project, guiRoot)

    # Give the location of the sequence file.
    seqDir  = '../data/seq'
    seqFile = os.path.join(seqDir, 'fasta.seq')

    # Call the readSequence method from the fastaFormat instance
    # and pass in the name of the file to read.
    fastaObj.readSequence(seqFile, minimalPrompts = 1)

    # Check to see if molSystem and molecule objects are in the CCPN project.
    print '\nProject object: [%s]' % project
    print "'list' of MolSystem objects: [%s]" % project.sortedMolSystems()
    print "'list' of Molecule objects: [%s]" % project.sortedMolecules()

    # Select first available molecule.
    molecule = project.findFirstMolecule()
    print 'Molecule type of first molecule: [%s]' % molecule.molType

    # Check that the residues have also been made.
    res1Type = molecule.findFirstMolResidue(seqCode = 1).ccpCode
    print 'Residue 1 type: [%s]\n' % res1Type

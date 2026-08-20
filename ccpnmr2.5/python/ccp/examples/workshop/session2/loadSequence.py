# Exercise 2.3: Load a fasta sequence from disk.
#
# - Use FormatConverter sub-routines to load a sequence located in
#     '../data/seq/fasta.seq'.
#
# - Use an instance of the FastaFormat class located in
#     ccpnmr.format.converters.FastaFormat, and then the readSequence()
#     method for that object.
#
# - Print out the molSystem and molecule information and navigate some
#     of the residue information.
#
# - Hint: create a Tk object from the Tkinter module when making the
#     FastaFormat instance and also create the FastaFormat object with
#     a newly made CCPN project.
#

import os

# Get Tkinter for popups.

import Tkinter

import memops.api.Implementation as Implementation

# Get FastaFormat class for format conversion.

from ccpnmr.format.converters.FastaFormat import FastaFormat

if __name__ == '__main__':

  project = Implementation.MemopsRoot(name = 'readFasta')

  # Create a Tkinter object.

  guiRoot = Tkinter.Tk()

  # Create a FastaFormat object and associate it with the CCPN project and
  #   Tkinter object.

  fastaObj = FastaFormat(project, guiRoot)

  # Give the location of the sequence file.

  seqDir  = os.path.join(os.path.abspath('..'), 'data', 'seq')

  seqFile = os.path.join(seqDir, 'fasta.seq')

  # Call the readSequence() method from the fastaFormat instance
  #   and pass in the name of the file to read.

  fastaObj.readSequence(seqFile, minimalPrompts = 1)

  # Check to see if molSystem and molecule objects are in the CCPN project.

  print '\nProject object: [%s]' % project

  print "\n'list' of MolSystem objects: [%s]" % project.sortedMolSystems()

  print "\n'list' of Molecule objects: [%s]" % project.sortedMolecules()

  # Select first available molecule.

  molecule = project.findFirstMolecule()

  print '\nMolecule type of first molecule: [%s]' % molecule.molType

  # Check that the residues have also been made.

  res1Type = molecule.findFirstMolResidue(seqCode = 1).ccpCode

  print '\nResidue 1 type: [%s]\n' % res1Type

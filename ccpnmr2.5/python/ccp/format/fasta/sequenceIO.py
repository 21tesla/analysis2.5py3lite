#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for Fasta sequence files

Copyright (C) 2005 Wim Vranken (European Bioinformatics Institute)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)
- PDBe website (http://www.ebi.ac.uk/pdbe/)

- contact Wim Vranken (wim@ebi.ac.uk)
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================
"""

import os, string

# Import general settings
from memops.universal.Util import returnInt
from memops.universal.Io import getTopDirectory

from ccp.format.fasta.generalIO import FastaGenericFile

from ccp.format.general.formatIO import Sequence, SequenceElement

#####################
# Class definitions #
#####################
      
class FastaSequenceFile(FastaGenericFile):
  # Information on file level

  def initialize(self):
    self.sequences = []

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s sequence file %s" % (self.format,self.name)

    sequenceRead = 0
    sequenceLine = ""

    fin = open(self.name, 'rU')

    lineErrors = []
    validLines = 0
    
    # Read, look for first line
    line = fin.readline()

    while line:
      
      validLines += 1
      cols = line.split()

      if self.patt['emptyline'].search(line):
        validLines -= 1
        pass

      elif line[0] == '>':

        # Designates line with comment info

        # If sequenceLine present, then parse this info (one sequence already read)

        if sequenceLine:
          self.sequences[-1].parseSequence(sequenceLine)

        sequenceRead = 1
        sequenceLine = ""
        self.sequences.append(FastaSequence(molName = line[1:]))

      elif sequenceRead:
        sequenceLine += cols[0]

      else:
        lineErrors.append(line)

      line = fin.readline()

    fin.close()

    if sequenceLine:
      self.sequences[-1].parseSequence(sequenceLine)
      
    #
    # Check
    #
    
    fileReadOk = True
    if len(lineErrors) > min(5,validLines * 0.5):
      fileReadOk = False
      self.sequences = []
    elif lineErrors:
      print "  Bad fasta format lines:%s" % self.newline
      for lineError in lineErrors:
        print lineError
      
    return fileReadOk

  def write(self,columnLength = 60, verbose = 0):

    if verbose == 1:
      print "Writing fasta sequence file %s" % self.name

    fout = open(self.name,'w')

    for sequence in self.sequences:

      fout.write("> %s" % sequence.molName)
      fout.write(self.newline)

      sequenceCode1Letter = ""
      for sequenceElement in sequence.elements:
        sequenceCode1Letter += sequenceElement.code1Letter

      for i in range(0,len(sequenceCode1Letter),columnLength):
        fout.write(sequenceCode1Letter[i:i+columnLength])
        fout.write(self.newline)

      fout.write(self.newline)

    fout.close()

class FastaSequence(Sequence):

  def parseSequence(self,sequenceLine):

    #
    # Parse sequence
    #
  
    seqCode = 1

    for code1Letter in sequenceLine:
      self.elements.append(FastaSequenceElement(seqCode,code1Letter))
      seqCode += 1

class FastaSequenceElement(SequenceElement):

  def setResidueCode(self,*args):
  
    code1Letter = args[0]
    
    if code1Letter:
      self.code1Letter = string.upper(code1Letter)
    else:
      self.code1Letter = 'X'

###################
# Main of program #
###################

if __name__ == "__main__":  
                                                      
  files = ['../reference/fasta/exampleMulti.seq']
  
  for file in files:
    
    file = os.path.join(getTopDirectory(), file)
    
    sequenceFile = FastaSequenceFile(file)

    sequenceFile.read(verbose = 1)

    for seq in sequenceFile.sequences:
      print seq.molName
      for seqEl in seq.elements:
        print seqEl.seqCode, seqEl.code1Letter

    sequenceFile.name = 'local/seq.testout'

    sequenceFile.write(verbose = 1)

    # Writing tested 25/04/2003 - all OK (Wim)

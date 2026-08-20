#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for Monte sequence export files

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

import string

# Import general functions
from memops.universal.Util import returnInt
from ccp.format.monte.generalIO import MonteGenericFile
from ccp.format.general.formatIO import Sequence, SequenceElement

#####################
# Class definitions #
#####################
      
class MonteSequenceFile(MonteGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.sequences = []

  def read(self,verbose = False):

    if verbose:
      print "Reading %s sequence file %s" % (self.format,self.name)

    self.sequences.append(MonteSequence())

    seqCode = 1
    badCodes = 0

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      cols = line.split()

      if len(cols) == 0 or self.patt['%sComment' % self.format].search(line):
        pass
      
      else:
      
        #
        # Assuming three-letter codes, with multiple ones on same line possible
        #

        for col in cols:
          
          tmpCode = col.strip()
          
          if len(tmpCode) == 3:
            self.sequences[-1].elements.append(MonteSequenceElement(seqCode,col))
            seqCode += 1
          else:
            badCodes += 1

      line = fin.readline()

    fin.close()
    
    #
    # Validity check
    #
    
    fileReadOk = True
    
    if seqCode == 1 or badCodes > seqCode * 0.1:
      self.sequences = []
      fileReadOk = False
      
    return fileReadOk

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing %s sequence file %s" % (self.format,self.name)

    if len(self.sequences) > 1:
      print "Warning: multiple sequences - writing to same file."        

    fout = open(self.name,'w')
    
    stepLen = 4

    for sequence in self.sequences:

      #
      # Write three letter codes (do 4 per line - see stepLen)
      #
      
      seqLen = len(sequence.elements)

      for i in range(0,seqLen,stepLen):
        for j in range(i,min(i+stepLen,seqLen)):
          residue = sequence.elements[j]

          fout.write("%3s " % (string.upper(residue.code3Letter)))
        fout.write(self.newline)

    fout.close()
    
MonteSequence = Sequence
MonteSequenceElement = SequenceElement

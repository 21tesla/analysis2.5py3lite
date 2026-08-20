"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for nmrView sequence file

Copyright (C) 2005-2009 Wim Vranken (European Bioinformatics Institute)

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

import os

# Import general functions
from memops.universal.Util import returnInt
from ccp.format.nmrView.generalIO import NmrViewGenericFile
from ccp.format.general.formatIO import Sequence, SequenceElement

#####################
# Class definitions #
#####################
      
class NmrViewSequenceFile(NmrViewGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.sequences = []

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading nmrView sequence file %s" % self.name

    self.sequences.append(NmrViewSequence())

    seqCode = 1
    
    lineErrors = []
    validLines = 0

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      cols = line.split()

      if len(cols) == 0 or self.patt['hash'].search(line):
        pass

      elif len(cols) > 2:
        lineErrors.append(line)
      
      else:
        validColumn = True
        
        if len(cols) == 2:
          newSeqCode = returnInt(cols[1],default=None,verbose=False)
        
          if newSeqCode != None:
            seqCode = newSeqCode
          else:
            validColumn = False
        
        if validColumn:
          validLines += 1
          self.sequences[-1].elements.append(NmrViewSequenceElement(seqCode,cols[0]))
          seqCode += 1
        else:
          lineErrors.append(line)

      line = fin.readline()

    fin.close()
      
    #
    # Check
    #
    
    if len(lineErrors) > min(5,validLines * 0.5):
      self.sequences = []
      print "  Bad %s format lines:%s" % (self.format,self.newline)
      for lineError in lineErrors:
        print lineError

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing nmrView sequence file %s" % self.name

    if len(self.sequences) > 1:
      print "Warning: multiple sequences - writing to same file."        

    fout = open(self.name,'w')

    for sequence in self.sequences:

      #
      # Write three letter codes (one per line)
      #

      oldSeqCode = 0

      for residue in sequence.elements:

        if oldSeqCode != residue.seqCode - 1:
          seqCodeString = " %d" % residue.seqCode
        else:
          seqCodeString = ""

        fout.write("%3s%s" % (residue.code3Letter.lower(),seqCodeString))
        fout.write(self.newline)

        oldSeqCode = residue.seqCode

    fout.close()

#
# Casting here for imports in ccpnmr.format.converters
#
 
NmrViewSequence = Sequence
NmrViewSequenceElement = SequenceElement

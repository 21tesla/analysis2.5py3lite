#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for autoAssign sequence from project file

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

import os,string

# Import general functions
from memops.universal.Util import returnInt
from ccp.format.autoAssign.generalIO import AutoAssignGenericFile
from memops.universal.Io import getTopDirectory
from ccp.format.general.formatIO import Sequence, SequenceElement

#####################
# Class definitions #
#####################
      
class AutoAssignSequenceFile(AutoAssignGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.sequences = []

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading autoAssign sequence file %s" % self.name

    molName = None

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      cols = line.split()

      if len(cols) == 0 or self.patt['hash'].search(line):
        pass

      else:
      
        if cols[0] == self.proteinNameTag:
          self.sequences.append(AutoAssignSequence(molName = cols[1]))
            
        elif cols[0] == self.sequenceTag:
          seqCode = returnInt(cols[1])
          seqString = cols[2]
          
          while (seqString[-1] != self.sequenceEnd):
            line = fin.readline()
            seqString += string.strip(line)
            
          seqString = seqString[:-1]
          
          #
          # Make the sequence
          #
          
          for code1Letter in seqString:

            self.sequences[-1].elements.append(AutoAssignSequenceElement(seqCode,code1Letter))
            seqCode += 1

      line = fin.readline()

    fin.close()

  def readFromShifts(self,shiftFile, verbose = 0):
  
    if verbose == 1:
      print "Extracting %s sequence from chemical shift file %s" % (self.format,shiftFile.name)

    (pathname,shiftFileName) = os.path.split(shiftFile.name)
    self.sequences.append(AutoAssignSequence(molName = shiftFileName))
    
    for (seqCode,code1Letter) in shiftFile.seqCodes:
      self.sequences[-1].elements.append(AutoAssignSequenceElement(seqCode,code1Letter))

#
# Casting here for imports in ccpnmr.format.converters
#
 
AutoAssignSequence = Sequence

class AutoAssignSequenceElement(SequenceElement):

  def setResidueCode(self,*args):
  
    code1Letter = args[0]
    self.code1Letter = string.upper(code1Letter)

###################
# Main of program #
###################

if __name__ == "__main__":  
                                                      
  files = ['../reference/autoassign/fas_table.aat']

  
  for inFile in files:
    
    sequenceFile = AutoAssignSequenceFile(os.path.join(getTopDirectory(),inFile))

    sequenceFile.read(verbose = 1)
    
    for seq in sequenceFile.sequences:
      print seq.molName
      for seqel in seq.elements:
        print seqel.seqCode, seqel.code1Letter
    

"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for PIPP sequence file

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
from ccp.format.pipp.generalIO import PippGenericFile
from ccp.format.general.formatIO import Sequence, SequenceElement
from ccp.format.general.Constants import defaultSeqInsertCode

#####################
# Class definitions #
#####################
      
class PippSequenceFile(PippGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.sequences = []

  def read(self,verbose = 0):
  
    #
    # TODO TODO not sure if this is correct!!! Check!
    #

    if verbose == 1:
      print "Reading %s sequence file %s" % (self.format,self.name)

    self.sequences.append(PippSequence())

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

      elif len(cols) <= 2:
        validLines += 1
        if len(cols) == 2:
          seqCode = returnInt(cols[1])

        self.sequences[-1].elements.append(PippSequenceElement(seqCode,cols[0]))
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

  def readFromShifts(self,shiftFile, verbose = 0):
  
    if verbose == 1:
      print "Extracting %s sequence from chemical shift file %s" % (self.format,shiftFile.name)

    seqCode = ""
    seqInsertCode = defaultSeqInsertCode
    segId = -1
    residueName = ""

    self.sequences.append(PippSequence(chainCode = shiftFile.defaultMolCode))
    
    for shift in shiftFile.chemShifts:

      if seqCode != shift.seqCode or seqInsertCode != shift.seqInsertCode:
      
        #
        # TODO: Insert empty ones if gaps for this format!?!?!
        #
        
        if seqCode != "" and shift.seqCode > seqCode + 1:
          for i in range(1,shift.seqCode - seqCode):
            fullSeqCode = str(seqCode + i) + shift.seqInsertCode
            self.sequences[-1].elements.append(PippSequenceElement(fullSeqCode,'XXX'))

        #
        # New residue/item
        # 
        
        seqCode = shift.seqCode
        seqInsertCode = shift.seqInsertCode

        fullSeqCode = str(seqCode) + seqInsertCode
        
        self.sequences[-1].elements.append(PippSequenceElement(fullSeqCode,shift.resLabel))
    
#
# Casting here for imports in ccpnmr.format.converters
#
 
PippSequence = Sequence
PippSequenceElement = SequenceElement

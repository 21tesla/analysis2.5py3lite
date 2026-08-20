"""
======================COPYRIGHT/LICENSE START==========================

assignmentIO.py: I/O for Mars spin system assignments

Copyright (C) 2007 Wim Vranken (European Bioinformatics Institute)

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

from ccp.format.mars.generalIO import MarsGenericFile

from memops.universal.Util import returnInt

#####################
# Class definitions #
#####################
      
class MarsAssignmentFile(MarsGenericFile):
  
  def initialize(self):
  
    self.assignments = []
    
  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s assignment file %s" % (self.format,self.name)

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      cols = line.split()

      if self.patt['emptyline'].search(line) or self.patt['hash'].search(line):
        pass
        
      if len(cols) == 3:
        (ccpCode,seqId) = cols[0].split(self.residueSep)
        (void,spinSystemId) = cols[1].split(self.residueSep)
        reliability = cols[2][1]
        
        self.assignments.append(MarsAssignment(returnInt(spinSystemId),returnInt(seqId),ccpCode,reliability))

      line = fin.readline()

    fin.close()

class MarsAssignment:

  def __init__(self,spinSystemId,seqId,ccpCode,reliability):

    self.spinSystemId = spinSystemId
    self.seqId = seqId
    self.ccpCode = ccpCode
    self.reliability = reliability
    

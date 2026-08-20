"""
======================COPYRIGHT/LICENSE START==========================

projectIO.py: I/O for Talos project file (with sequence and shifts)

Copyright (C) 2005-2010 Wim Vranken (European Bioinformatics Institute)

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

# Import general functions
from memops.universal.Util import returnFloat, returnInt
from ccp.format.talos.generalIO import TalosGenericFile

from ccp.format.talos.chemShiftsIO import TalosChemShiftFile
from ccp.format.talos.sequenceIO import TalosSequenceFile

#####################
# Class definitions #
#####################

class TalosProjectFile(TalosGenericFile):

  # Information on file level

  def initialize(self):

    self.chemShiftFiles = []
    self.sequenceFile = None

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading Talos project file %s" % self.name

    #
    # Read the sequence from the project file
    #
    
    self.sequenceFile = TalosSequenceFile(self.name)
    self.sequenceFile.read()
    
    #
    # Read the chemical shifts from the project file
    #

    self.chemShiftFiles.append(TalosChemShiftFile(self.name))
    self.chemShiftFiles[-1].read()

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing Talos project file %s" % self.name

    fout = open(self.name,'w')

    #
    # Write out header
    #

    fout.write("REMARK File written by CcpNmrFormat converter.")
    fout.write(self.newline * 2)
    
    #
    # Read the sequence from the project file
    #
    
    self.sequenceFile.write(use_fout = fout)
    
    #
    # Read the chemical shifts from the project file
    #
    
    for chemShiftFile in self.chemShiftFiles:
      chemShiftFile.write(use_fout = fout)
      
    fout.close()

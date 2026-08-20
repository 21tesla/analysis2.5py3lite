"""
======================COPYRIGHT/LICENSE START==========================

projectIO.py: I/O for Cosmos project file (with coordinates and shifts)

Copyright (C) 2008 Wim Vranken (European Bioinformatics Institute)

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
from ccp.format.cosmos.generalIO import CosmosGenericFile

from ccp.format.cosmos.chemShiftsIO import CosmosChemShiftFile
from ccp.format.cosmos.coordinatesIO import CosmosCoordinateFile

#####################
# Class definitions #
#####################

class CosmosProjectFile(CosmosGenericFile):

  # Information on file level

  def initialize(self):

    # For coo file
    self.chemShiftFile = None
    self.sequenceFile = None
    self.constraintFile = None

    # For cod file
    self.rcOptMolecule = 'ALL'
    self.rcOptCutoff = -1.0
    self.namesOptFitTo = 'EXACT'
    self.groupMembers = [('^C.*',1),('^H.*',2),('^[CN][A_].*',3)]
    
  def write(self, version='003', verbose=False, fileType='COO'):

    if verbose == 1:
      print "Writing %s .%s project file %s" % (self.format,fileType,self.name)

    fout = open(self.name,'w')

    if fileType == 'COO':
    
      #
      # Write out header
      #

      fout.write("$COO%s" % version)
      fout.write(self.newline)

      fout.write("REM File written by CcpNmrFormat converter.")
      fout.write(self.newline)

      #
      # Write the coordinates from the project file
      #

      self.coordinateFile.write(use_fout = fout)

      #
      # Write the chemical shifts from the project file if available
      #
      
      if self.chemShiftFile:

        self.chemShiftFile.write(use_fout = fout)
    
    elif fileType == 'COD':
    
      #
      # Write out header and data options
      #

      fout.write("$COD%s" % version)
      fout.write(self.newline)

      fout.write("REMARK File written by CcpNmrFormat converter.")
      fout.write(self.newline)
      
      # Only constraints within a molecule: IN_MOL, or all constraints: ALL (default)
      fout.write("RC_OPT_MOLECULE %s" % self.rcOptMolecule)
      fout.write(self.newline)
      
      # Cutoff to skip constraints that do not fit structure, set to -1.0 to ignore cutoff
      fout.write("RC_OPT_CUTOFF %.1f" % self.rcOptCutoff)
      fout.write(self.newline)
      
      # Regular expression setting; leave to EXACT
      fout.write("NAMES_OPT_FIT_TO %s" % self.namesOptFitTo)
      fout.write(self.newline)

      # End of data options
      fout.write("DATA_OPT_END")
      fout.write(self.newline)
      
      #
      # Write the group members info, not relevant for now
      #

      fout.write("GROUP_MEMBERS %d" % len(self.groupMembers))
      fout.write(self.newline)
      
      for (groupRegExp,groupNumber) in self.groupMembers:
        fout.write("%s %d%s" % (groupRegExp,groupNumber,self.newline))
        
      
      #
      # Write the constraints
      #
      
      if self.constraintFile:
        self.constraintFile.write(use_fout = fout,rcOptCutoff=self.rcOptCutoff)
 
      #
      # Write the chemical shifts
      #

      if self.chemShiftFile:
        self.chemShiftFile.write(use_fout = fout,useSerial=False)
    
    
    fout.write("END")
    fout.write(self.newline)
    
    fout.close()

"""
======================COPYRIGHT/LICENSE START==========================

distanceConstraintsIO.py: I/O for Gromacs distance constraint files

Copyright (C) 2012 Wim Vranken (Vrije Universiteit Brussel)

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

from memops.universal.Io import getTopDirectory

from memops.universal.Util import returnInt, returnFloat

from ccp.format.gromacs.generalIO import GromacsGenericFile
from ccp.format.gromacs.generalIO import GromacsConstraintItem
from ccp.format.gromacs.generalIO import GromacsConstraintMember

#####################
# Class definitions #
#####################

class GromacsDistanceConstraintFile(GromacsGenericFile):

  def initialize(self):
  
    self.constraints = []
    self.constraintElements = 2
    
    self.constraintFileType = 'distance_restraints'
 
  """
  def read(self,verbose = False):

    if verbose:
    
      print "Reading %s distance constraint list %s" % (self.format,self.name)
    
    fin = open(self.name, 'rU')
     
  """
      
  def write(self,verbose = False):

    """
    Note: in example the column widths don't seem to matter much; as long as they're separated! See local/ dir for explanation
    """
    
    if not self.constraints:
      return

    if verbose: 
      print "Writing %s distance constraint list %s" % (self.format,self.name)
    
    fout = open(self.name,'w')

    #
    # Write out distance constraints
    #
    
    fout.write("[ {} ]{}".format(self.constraintFileType,self.newline))
    fout.write(";   ai     aj  type index type'  low    up1    up2  fac" + self.newline)
    
    for constraint in self.constraints:
      # Note: items should be expanded so that one atom member on every side!
      for item in constraint.items:

        for i in range(self.constraintElements):
        
          member = item.members[i]          
          fout.write("{:7d}".format(member.atomSerial))

        fout.write("  1 {:7d}  1  {:6.3f} {:6.3f} {:6.3f} 1.0".format(constraint.serial,constraint.lowerLimit,constraint.upperLimit1,constraint.upperLimit2))
        fout.write(self.newline)


class GromacsDistanceConstraint:

  def __init__(self,serial):
    
    self.serial = serial
    self.items = []
  
  def setDistances(self,lowerLimit,upperLimit1,upperLimit2):
    
    self.lowerLimit = lowerLimit
    self.upperLimit1 = upperLimit1
    self.upperLimit2 = upperLimit2
    
    

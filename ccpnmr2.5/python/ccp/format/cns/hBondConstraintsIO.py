#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

hBondConstraintsIO.py: I/O for XPLOR/CNS h bond constraint files

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

import os

from memops.universal.Io import getTopDirectory

from ccp.format.cns.distanceConstraintsIO import CnsDistanceConstraintFile
from ccp.format.cns.distanceConstraintsIO import CnsDistanceConstraint

#####################
# Class definitions #
#####################

class CnsHBondConstraintFile(CnsDistanceConstraintFile):

  def checkLinePattern(self,line):

    return 1

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading cns H bonds constraint list %s" % self.name

    return self.readGeneric(CnsHBondConstraint)
    
  def printWriteText(self):
  
      print "Writing cns H bonds constraint list %s" % self.name

class CnsHBondConstraint(CnsDistanceConstraint):

  def setOtherData(self,*args):

    pass
  
###################
# Main of program #
###################

if __name__ == "__main__":

  # NO REAL EXAMPLE FILE AVAILABLE!
  files = ['../reference/aria/ambig.restraints']
  
  for file in files:
    
    file = os.path.join(getTopDirectory(), file)
    
    constraintFile = CnsHBondConstraintFile(file)

    constraintFile.read(verbose = 1)
  
    for constraint in constraintFile.constraints:
      print constraint.Id, constraint.origId,
    
      print constraint.targetDist, constraint.minusDist, constraint.plusDist,
    
      for item in constraint.items:
        for member in item.members:
          print member.seqCode, member.atomName,
        print "|",
      
      print
  
    constraintFile.name = 'local/hbond.testout'
  
    constraintFile.write(verbose = 1)
  

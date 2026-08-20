"""
======================COPYRIGHT/LICENSE START==========================

distanceConstraintsIO.py: I/O for AMBER distance constraint files

Copyright (C) 2010 Wim Vranken (European Bioinformatics Institute)

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

from ccp.format.amber.generalIO import AmberGenericFile
from ccp.format.amber.generalIO import AmberGenericConstraint

from ccp.format.amber.coordinatesIO import AmberCoordinateFile

#####################
# Class definitions #
#####################

class AmberDistanceConstraintFile(AmberGenericFile):

  infoType = 'distance'

  def initialize(self):
  
    self.constraints = []
    self.constraintElements = 2

    self.constraintType = None
    
    self.constraintFileType = 'distance_restraints'
    
  def read(self, coordFile=None, verbose=False):

    if verbose == 1:
    
      print "Reading %s %s constraint list %s" % (self.format,self.infoType,self.name)
    
    constraintInfoList = self.readConstraints()
    
    coordinateFile = AmberCoordinateFile(coordFile)
    coordinateFile.read(maxNum = 1) # Only need one model!
      
    for i in range(len(constraintInfoList)):
     
      constraintInfo = constraintInfoList[i]

      atomSerials = constraintInfo['iat']
      coordinateAtoms = []
      
      for atomSerial in atomSerials:
        coordinateAtoms.append(coordinateFile.serialToCoord[atomSerial])
      
      dihedralConstraint = AmberDistanceConstraint(self,i+1)
      self.constraints.append(dihedralConstraint) 
      
      dihedralConstraint.addItem(coordinateAtoms)
      dihedralConstraint.addDistance(constraintInfo)
    
class AmberDistanceConstraint(AmberGenericConstraint):
  
  def addDistance(self,constraintInfo):
  
    self.lowerDist = constraintInfo['r2'][0]
    self.upperDist = constraintInfo['r3'][0]

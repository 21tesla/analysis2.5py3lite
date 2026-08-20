#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

dihedralConstraintsIO.py: I/O for Discover dihedral constraints file

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
from memops.universal.Util import returnFloat
from memops.universal.Util import returnInt
from ccp.format.discover.generalIO import DiscoverGenericFile, DiscoverGenericConstraint
from memops.universal.Io import getTopDirectory

from ccp.format.general.Util import getSeqAndInsertCode
from ccp.format.general.Constants import defaultMolCode

#####################
# Class definitions #
#####################
    
class DiscoverDihedralConstraintFile(DiscoverGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.constraints = []
    self.constraintElements = 4

    self.constraintType = None

  def read(self,verbose = 0):

    if verbose == 1:
    
      print "Reading %s distance constraint list %s" % (self.format,self.name)
    
    fin = open(self.name, 'rU')
     
    #
    # Start reading...
    #
    
    line = fin.readline()
        
    constraintNum = 0
    
    while line:

      if self.patt['emptyline'].search(line) or self.patt['exclamation'].search(line):
        line = fin.readline()
        continue
      
      elif self.patt['hash'].search(line):
        
        (void,self.constraintType) = string.split(line,'#')
        line = fin.readline()
        continue

      cols = string.split(line)

      if cols[0] != '+':
        constraintNum += 1
      
        constraint = DiscoverDihedralConstraint(self,constraintNum)
        self.constraints.append(constraint)
        
        constraint.addInfo(cols[self.constraintElements:])
        constraint.addItem(cols)
      
      else:
        constraint.addItem(cols[1:])
      
      line = fin.readline()

class DiscoverDihedralConstraint(DiscoverGenericConstraint):
   
  def addInfo(self,cols):
  
    self.JValue = returnFloat(cols[0])
    
    self.lowerAngle = returnFloat(cols[5])
    self.upperAngle = returnFloat(cols[6])
  
###################
# Main of program #
###################

if __name__ == "__main__":  
                                                      
  files = ['../../reference/discover/1ad7.dihed']
  
  for inFile in files:
    
    constraintFile = DiscoverDihedralConstraintFile(os.path.join(getTopDirectory(),inFile))

    constraintFile.read(verbose = 1)
    
    for constraint in constraintFile.constraints:
      print constraint.Id,
      print constraint.upperAngle, constraint.lowerAngle
      for constraintItem in constraint.items:
        for constraintMember in constraintItem.members:
          print constraintMember.chainCode, constraintMember.seqCode, constraintMember.atomName,
        print

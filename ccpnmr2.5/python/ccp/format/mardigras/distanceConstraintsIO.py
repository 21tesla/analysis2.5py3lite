"""
======================COPYRIGHT/LICENSE START==========================

distanceConstraintsIO.py: I/O for Mardigras distance constraint files (very old AMBER format)

Copyright (C) 2009 Wim Vranken (European Bioinformatics Institute)

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

from memops.universal.Util import returnInt, returnFloat

from ccp.format.mardigras.generalIO import MardigrasGenericFile
from ccp.format.mardigras.generalIO import MardigrasConstraintItem
from ccp.format.mardigras.generalIO import MardigrasConstraintMember

#####################
# Class definitions #
#####################

class MardigrasDistanceConstraintFile(MardigrasGenericFile):

  def initialize(self):
  
    self.constraints = []
    self.constraintElements = 2
    
    self.constraintFileType = 'distance_restraints'
    
    self.addColumnNames = ['index','restrnum']
    
  def read(self,verbose = 0):

    if verbose == 1:
    
      print "Reading %s distance constraint list %s" % (self.format,self.name)
    
    fin = open(self.name, 'rU')

    constraintNum = 1

    line = fin.readline()
        
    while line:

      if self.patt['emptyline'].search(line) or self.patt['hash'].search(line):
        line = fin.readline()
        continue

      values = line.split()
      
      # Ignore first line, expected to be:
      #
      #  atom1  residue1  atom2  residue2  Lower Bound  Upper Bound
      #
      
      if values[0] != 'atom1':
        
        constraint = MardigrasDistanceConstraint(self,constraintNum)
        self.constraints.append(constraint)
        
        constraint.addItem(values[:4])
        constraint.setDistances(values[-2:])
        
        constraintNum += 1
      
      line = fin.readline()
      
  """  
  def write(self,verbose = 0):
    
    if not self.constraints:
      return

    if verbose == 1:
 
      print "Writing %s distance constraint list %s" % (self.format,self.name)
    
    
    fout = open(self.name,'w')
    
    fout.close()
    
  """
    
class MardigrasDistanceConstraint:

  def __init__(self,parent,Id):
    
    self.parent = parent
    self.Id = returnInt(Id)
    
    self.items = []
   
  def setDistances(self,values):
    
    self.lowerDist = returnFloat(values[0])
    self.upperDist = returnFloat(values[1])

  def addItem(self,values):
  
    item = MardigrasConstraintItem()
    self.items.append(item)
    
    for i in range(2):
    
      member = MardigrasConstraintMember()
      item.members.append(member)
      
      member.setInfo(values[i*2 + 1],values[i*2])

"""
======================COPYRIGHT/LICENSE START==========================

distanceConstraintsIO.py: I/O for Cosmos distance constraints

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

import os

from memops.universal.Util import returnInt, returnFloat

from ccp.format.cosmos.generalIO import CosmosGenericFile

#####################
# Class definitions #
#####################

class CosmosDistanceConstraintFile(CosmosGenericFile):

  def initialize(self):
  
    self.constraints = []
    self.constraintElements = 2
    
    self.constraintFileType = 'distance_restraints'
    
  def write(self,verbose=False,use_fout=None,rcOptCutoff=-1.0):
    
    if not self.constraints:
      return

    if verbose == 1:
 
      print "Writing %s distance constraint list %s" % (self.format,self.name)
    

    #
    # Write out distance constraints
    #
    
    fileLines = []
    
    for constraint in self.constraints:
      
      # Collating all assignments - there is no OR switching?
      assignments = [[],[]]
      
      for item in constraint.items:      
        for i in range(0,self.constraintElements):
        
          member = item.members[i]
          
          atomOrSetName = member.atomName
         
          assignment = "%s_%s_%s_%d" % (atomOrSetName,member.residueName,member.chainCode,member.seqCode) 
          
          if assignment not in assignments[i]:
            assignments[i].append(assignment)
                  
      if assignments[0] and assignments[1]:
        fileLine = ""
        for assignmentList in assignments:

          if len(assignmentList) > 1:
            fileLine += "("

          fileLine+= "|".join(assignmentList)

          if len(assignmentList) > 1:
            fileLine += ")"
          
          fileLine += " "
          
        if constraint.upperLimit:
          distance = constraint.upperLimit
        elif rcOptCutoff > 0:
          distance = constraint.targetDistance + rcOptCutoff

        fileLine += "{:.1f}(1.0)".format(distance) #,rcOptCutoff)
        
        fileLines.append(fileLine)
  
    #
    # Now write file itself
    #
    
    totalConstraints = len(fileLines)     
          
    if use_fout:
      fout = use_fout
    else:
      fout = open(self.name,'w')

    fout.write("R_CONSTRAINTS %d" % totalConstraints)
    fout.write(self.newline)

    for fileLine in fileLines:
      fout.write(fileLine)
      fout.write(self.newline)

    if not use_fout:
      fout.close()    

class CosmosDistanceConstraint:

  def __init__(self,Id):
    
    self.Id = returnInt(Id)
    
    self.items = []
      
  def setDistanceData(self,upperLimit,targetDistance):
    
    self.upperLimit = upperLimit
    self.targetDistance = targetDistance


#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

hBondConstraintsIO.py: I/O for Diana/Dyana h bond constraint files

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

from ccp.format.dyana.distanceConstraintsIO import DyanaDistanceConstraintFile
from ccp.format.dyana.distanceConstraintsIO import DyanaDistanceConstraint

from ccp.format.dyana.generalIO import DyanaConstraintItem, DyanaConstraintMember

#####################
# Class definitions #
#####################

class DyanaHBondConstraintFile(DyanaDistanceConstraintFile):

  defaultHBondDistance = 2.4
  
  def read(self,verbose = 0):

    if verbose == 1:
      fileText = self.name
      
      if self.lowerDistanceFile:
        fileText += ' and %s.' % self.lowerDistanceFile
        
      print "Reading %s H bonds constraint list %s" % (self.format,fileText)
      
    # First try to read in as hbonds type file
    if not self.readHBondSpecificFile(self.name):
      self.readGeneric(DyanaDistanceConstraint)

  def readHBondSpecificFile(self,fileName):
    
    """
File format is:

hbond H 3 O 68
hbond H 68 O 3
hbond H 65 O 5
hbond H 7 O 63
hbond H 63 O 7
hbond H 4 O 42
hbond H 42 O 4
hbond H 6 O 40    
    """    
    
    returnStatus = False

    #
    # First check whether this is a custom hbond file without distances, read in if so
    #
    
    fin = open(fileName, 'rU')
    lines = fin.readlines()
    fin.close()
    
    self.restrNum = 0
    
    for line in lines:
      cols = line.split()
      if cols[0] == 'hbond':
        returnStatus = True
        
        self.restrNum += 1
        
        # TODO: should probably also get the connected heavy atom, and set this to 3.3 angstrom!
        constraint = DyanaHBondConstraint(self.restrNum)
        self.constraints.append(constraint)

        constraint.setDistanceData(self.defaultHBondDistance)

        constraint.items.append(DyanaConstraintItem())
         
        for i in range(1,5,2):
          atomName = cols[i]
          seqCode = cols[i+1]
        
          constraint.items[-1].members.append(DyanaConstraintMember(self.defaultMolCode,seqCode,None,atomName))
    
    return returnStatus
   
  def write(self,verbose = 0):

    if verbose == 1:
      fileText = self.name
      
      if self.lowerDistanceFile:
        fileText += ' and %s.' % self.lowerDistanceFile

      print "Writing %s H bonds constraint list %s" % (self.format,fileText)

    # Writer is same as for distance constraints...

    self.writeGeneric(noAmbiguous = 1)    
    
class DyanaHBondConstraint(DyanaDistanceConstraint):

  pass
  
###################
# Main of program #
###################

if __name__ == "__main__":

  # NO REAL EXAMPLE FILE AVAILABLE!
  files = ['../reference/diana/u1c.upl']
  
  for file in files:
    
    file = os.path.join(getTopDirectory(), file)
    
    constraintFile = DyanaHBondConstraintFile(file)

    constraintFile.read(verbose = 1)
  
    for constraint in constraintFile.constraints:
      print constraint.Id,

      print constraint.upperDist,
      #print constraint.peakNum, constraint.peakVol, constraint.ppms
      for item in constraint.items:
        for member in item.members:
          print member.seqCode, member.atomName,
        print "|",

      print

    constraintFile.name = 'local/hbonds.testout'

    constraintFile.write(verbose = 1)


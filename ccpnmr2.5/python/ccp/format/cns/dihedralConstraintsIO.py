"""
======================COPYRIGHT/LICENSE START==========================

dihedralConstraintsIO.py: I/O for XPLOR/CNS dihedral constraint files

Copyright (C) 2005-2009 Wim Vranken (European Bioinformatics Institute)

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

from memops.universal.Util import returnInt
from memops.universal.Util import returnFloat

from ccp.format.cns.generalIO import CnsGenericFile

#####################
# Class definitions #
#####################

class CnsDihedralConstraintFile(CnsGenericFile):

  def initialize(self):

    self.constraints = []
    
    self.constraintElements = 4

  def checkLinePattern(self,line):

    return 1

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading cns dihedral constraint list %s" % self.name

    return self.readGeneric(CnsDihedralConstraint)


  def write(self,verbose = 0):

    #
    # Output format is (example following):
    #
    # ASSIGN (resid     4 and name     N) (resid     4 and name    CA)
    #        (resid     4 and name    CB) (resid     4 and name   CG1)  1.00 180.00 20.00 2
    #
    # ...

    if verbose == 1:
      print "Writing cns dihedral constraint list %s" % self.name

    fout = open(self.name,'w')

    for constraint in self.constraints:

      fout.write(" ASSIGN ")

      if len(constraint.items) > 1:
        print "  Warning: multiple items for CNS dihedral restraint - only using first one"

      item = constraint.items[0] 

      for j in range(0,len(item.members)):

        member = item.members[j]

        if member.chainCode != self.defaultMolCode:
          # TODO: this is not correct! Contact Jens!!
          segIdString = "segid \"%4s\" and " % member.chainCode
        else:
          segIdString = ""

        if j == 2:
          startSpace = '        '
        else:
          startSpace = ''

        fout.write("%s %s%sresid %-4d and name %-4s%s" % (startSpace,
    		                  '(',
    		                  segIdString,
                    	                  member.seqCode,
                    	                  member.atomName,
                    	                  ')'))
        if j == 1:
          fout.write(self.newline)

        if j == 3:
          fout.write(" %5.2f %6.2f %5.2f %d" % (constraint.energyCst,
    		        constraint.targetAngle,
    		        constraint.devAngle,
                    	        constraint.funcExp))
          fout.write(self.newline)

      fout.write(self.newline)

    fout.close()


class CnsDihedralConstraint:

  def __init__(self,Id,origId,patt=None,format=None):
    
    self.Id = returnInt(Id)
    
    if origId != None:
      self.origId = returnInt(origId)
    else:
      self.origId = None
    
    self.items = []
    
    self.patt = patt
    self.format = format
    
  def setAngleData(self,targetAngle,devAngle,energyCst,funcExp):
  
    self.targetAngle = returnFloat(targetAngle)
    self.devAngle = returnFloat(devAngle)
    self.energyCst = returnFloat(energyCst)
    self.funcExp = returnInt(funcExp)
        
  def getSpecificInfo(self,assiLine): 

    #
    # Get distances and peaknum, ...
    #
  
    anglePatt = self.patt[self.format + 'RestrAngles'].search(assiLine)
    
    lineRead = False
    
    if anglePatt:
      
      energyCst = anglePatt.group(1)
      targetAngle = anglePatt.group(2)
      devAngle = anglePatt.group(3)
      funcExp = anglePatt.group(4)
  
      self.setAngleData(targetAngle,devAngle,energyCst,funcExp)
      
      lineRead = True
    
    return lineRead
  

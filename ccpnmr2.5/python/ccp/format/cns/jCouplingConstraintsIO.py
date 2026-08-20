#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

jCouplingConstraintsIO.py: I/O for XPLOR/CNS j coupling constraint files

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

import os, string

from memops.universal.Io import getTopDirectory

from memops.universal.Util import returnInt
from memops.universal.Util import returnFloat, returnFloats

from ccp.format.cns.generalIO import CnsGenericFile

#####################
# Class definitions #
#####################

class CnsJCouplingConstraintFile(CnsGenericFile):

  def initialize(self):

    self.constraints = []
    
    self.constraintElements = 4

  def checkLinePattern(self,line):
  
    if string.find(line,'coefficients') != -1:
      cols = string.split(line)
      
      self.coeffs = returnFloats(cols[1:])
      
      return -1

    return 1

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading cns j coupling constraint list %s" % self.name

    return self.readGeneric(CnsJCouplingConstraint)

  def write(self,verbose = 0):

    #
    # Output format is (example following):
    #
    # ASSIGN (resid     4 and name     N) (resid     4 and name    CA)
    #        (resid     4 and name    CB) (resid     4 and name   CG1)  1.00 180.00 20.00 2
    #
    # ...

    if verbose == 1:
      print "Writing cns j coupling constraint list %s" % self.name

    fout = open(self.name,'w')

    fout.write("coefficients")

    for coeff in self.coeffs:
      fout.write(" %.2f" % coeff)

    fout.write(self.newline * 2)

    for constraint in self.constraints:

      fout.write(" ASSIGN ")

      for i in range(0,len(constraint.items)): # Should be only 1

        item = constraint.items[i] 

        for j in range(0,len(item.members)):

          member = item.members[j]

          if member.chainCode != self.defaultMolCode:
            # TODO: this is not correct!?! Contact Jens!!
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
            fout.write(" %7.2f  %7.2f" % (constraint.value,
                                          constraint.error))
            fout.write(self.newline)

      fout.write(self.newline)

    fout.close()

class CnsJCouplingConstraint:

  def __init__(self,Id,origId,patt = None,format = None):
    
    self.Id = returnInt(Id)
    
    if origId != None:
      self.origId = returnInt(origId)
    else:
      self.origId = None
    
    self.items = []
    
    self.patt = patt
    self.format = format
    
  def setJCouplingData(self,value,error):
  
    self.value = returnFloat(value)
    self.error = returnFloat(error)
        
  def getSpecificInfo(self,assiLine): 

    #
    # Get distances and peaknum, ...
    #
  
    couplingPatt = self.patt[self.format + 'RestrCoupling'].search(assiLine)
  
    lineRead = False
    
    if couplingPatt:
      
      value = couplingPatt.group(1)
      error = couplingPatt.group(2)
  
      self.setJCouplingData(value,error)
      
      lineRead = True
      
    return lineRead
     
###################
# Main of program #
###################

if __name__ == "__main__":

  files = ['../reference/aria/phi.restraints']
  
  for file in files:
    
    file = os.path.join(getTopDirectory(), file)
    
    constraintFile = CnsJCouplingConstraintFile(file)

    constraintFile.read(verbose = 1)

    for constraint in constraintFile.constraints:
      print constraint.Id, constraint.origId,

      print constraint.value, constraint.error,

      for item in constraint.items:
        for member in item.members:
          print member.seqCode, member.atomName,
        print "|",

      print

    constraintFile.name = 'local/jcoup.testout'

    constraintFile.write(verbose = 1)


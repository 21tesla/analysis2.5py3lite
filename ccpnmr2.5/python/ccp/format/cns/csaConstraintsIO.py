"""
======================COPYRIGHT/LICENSE START==========================

csaConstraintsIO.py: I/O for XPLOR/CNS CSA constraint files

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

from memops.universal.Util import returnInt
from memops.universal.Util import returnFloat

from ccp.format.cns.generalIO import CnsGenericFile

#####################
# Class definitions #
#####################

class CnsCsaConstraintFile(CnsGenericFile):

  def initialize(self):

    self.constraints = []
    
    self.constraintElements = 5

  def checkLinePattern(self,line):
  
    return True

  def read(self,verbose = False):

    if verbose:
      print "Reading %s CSA constraint list %s" % (self.format,self.name)

    return self.readGeneric(CnsCsaConstraint)

  def write(self,includeWeight = False, verbose = 0):

    #
    # Output format is (example following):
    #
    # assign ( resid 501  and name OO   )
    #  ( resid 501  and name Z   )
    #  ( resid 501  and name X   )
    #  ( resid 501  and name Y   )
    #  ( resid   4  and name CA   )   -0.0400  0.15000
    # ...

    if verbose == 1:
      print "Writing %s CSA constraint list %s" % (self.format,self.name)

    fout = open(self.name,'w')

    for constraint in self.constraints:

      fout.write(" ASSIGN ")

      for i in range(0,len(constraint.items)): # Should be only 1

        item = constraint.items[i] 

        for j in range(0,len(item.members)):

          member = item.members[j]

          if member.chainCode != self.defaultMolCode:
            # TODO: this is not correct! Contact Jens!!
            segIdString = "segid \"%4s\" and " % member.chainCode
          else:
            segIdString = ""

          if j >= 1:
            startSpace = '        '
          else:
            startSpace = ''

          fout.write("%s %s%sresid %-4d and name %-4s%s" % (startSpace,
                                        	'(',
                                        	segIdString,
			member.seqCode,
			member.atomName,
			')'))
          if j < 4:
            fout.write(self.newline)

          if j == 4:
            fout.write(" %7.4f %7.4f" % (constraint.value,
                                         constraint.error,
                                         ))
            
            if includeWeight:
              fout.write(" %7.4f" % constraint.weight)
            
            fout.write(self.newline)

      fout.write(self.newline)

    fout.close()

class CnsCsaConstraint:

  def __init__(self,Id,origId,patt=None,format=None):
    
    self.Id = returnInt(Id)
    
    if origId != None:
      self.origId = returnInt(origId)
    else:
      self.origId = None
    
    self.items = []
    
    self.patt = patt
    self.format = format
    
  def setCsaData(self,value,error,weight):
  
    self.value = returnFloat(value)
    self.error = returnFloat(error)
    self.weight = returnFloat(weight)
        
  def getSpecificInfo(self,assiLine): 

    #
    # Get distances and peaknum, ...
    #
  
    csaPatt = self.patt[self.format + 'RestrCsa'].search(assiLine)
  
    lineRead = False
    
    if csaPatt:

      value = csaPatt.group(1)
      error = csaPatt.group(2)
      
      if csaPatt.group(3):
        weight = csaPatt.group(3)
      else:
        weight = 1.0

      self.setCsaData(value,error,weight)
    
      lineRead = True
      
    return lineRead
  

#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

rdcConstraintsIO.py: I/O for Module program rdc constraint files

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
Ionides and Ernest D. Laue. The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Accepted by Proteins (2004).

===========================REFERENCE END===============================
"""

import os, string

from memops.universal.Util import returnInt, returnFloat
from memops.universal.Io import getTopDirectory

from ccp.format.module.generalIO import ModuleGenericFile
from ccp.format.module.generalIO import ModuleConstraintItem
from ccp.format.module.generalIO import ModuleConstraintMember

#####################
# Class definitions #
#####################

class ModuleRdcConstraintFile(ModuleGenericFile):

  def initialize(self):

    self.constraints = []

    self.constraintElements = 2

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s rdc constraint list %s" % (self.format,self.name)

    rdcId = 0
    chainCode = self.defaultMolCode

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      cols = line.split()

      if len(cols) == 0 or self.patt['hash'].search(line):
        pass

      else:

        rdcId += 1

        keywds = {}

        if len(cols) > 5:
          keywds['error'] = cols[5]
          if len(cols) > 6:
            keywds['energyCst'] = cols[6]

        self.constraints.append(ModuleRdcConstraint(rdcId))
        self.constraints[-1].setRdcData(cols[4],**keywds)
        self.constraints[-1].setAtomMembers(chainCode,cols[0],cols[1],cols[2],cols[3])

      line = fin.readline()

  def write(self,verbose = 0):

    #
    # Output format is (example following):
    #
    # 98    n  98  hn    -0.000    0.5    1.00
    #
    # ...

    if verbose == 1:
      print "Writing %s rdc constraint list %s" % (self.format,self.name)

    fout = open(self.name,'w')

    for constraint in self.constraints:

      for item in constraint.items:
        for member in item.members:
          atomName = member.atomName
          seqCode = member.seqCode

          fout.write("%3d%4s " % (seqCode,atomName))
        
        if constraint.error:
          error = constraint.error
        else:
          error = 0.5

        if constraint.energyCst:
          energyCst = constraint.energyCst
        else:
          energyCst = 1.00
        
        fout.write("%9.3f %6.1f %7.2f" % (constraint.value,error,energyCst))
        fout.write(self.newline)

    fout.close()

class ModuleRdcConstraint:

  def __init__(self,Id):
    
    self.Id = returnInt(Id)
    self.items = []
    
  def setRdcData(self,value,error = 0.00, energyCst = 1.0):
  
    self.error = returnFloat(error)
    self.value = returnFloat(value)
    self.energyCst = returnFloat(energyCst)
      
  def setAtomMembers(self,chainCode,seqCode1,atomName1,seqCode2,atomName2): 

    self.items.append(ModuleConstraintItem())
     
    self.items[-1].members.append(ModuleConstraintMember(chainCode,returnInt(seqCode1),string.upper(atomName1)))
    self.items[-1].members.append(ModuleConstraintMember(chainCode,returnInt(seqCode2),string.upper(atomName2)))
 
###################
# Main of program #
###################

if __name__ == "__main__":

  files = ['../../reference/ccpNmr/stuartSoton/cam_rdc.mod']
  
  for file in files:
    
    file = os.path.join(getTopDirectory(), file)
    
    constraintFile = ModuleRdcConstraintFile(file)

    constraintFile.read(verbose = 1)

    for constraint in constraintFile.constraints:
      print constraint.Id,

      print constraint.value, constraint.error,

      for item in constraint.items:
        for member in item.members:

          print member.seqCode, member.atomName,

        print "|",

      print
    
    constraintFile.name = 'local/rdc.testout'

    constraintFile.write(verbose = 1)


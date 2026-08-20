"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for Cosmo chemical shift file

Copyright (C) 2008 Wim Vranken (European Bioinformatics Institute)

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
from memops.universal.Util import returnFloat, returnFloats
from memops.universal.Util import returnInt
from ccp.format.cosmos.generalIO import CosmosGenericFile
from memops.universal.Io import getTopDirectory

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################
    
class CosmosChemShiftFile(CosmosGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.chemShifts = []
    self.seqNameSep = '.'

  #def read(self,verbose = 0):
    
  #  pass
    
  def write(self, useSerial=True, verbose=False, use_fout=None): #, useResLabel=False

    if verbose == 1:
      print "Writing %s chemical shift list %s" % (self.format,self.name)
                 
    #
    # Store chem shifts lines in list for writing
    #

    fileLines = []

    for chemShift in self.chemShifts:
    
      if useSerial:
        if chemShift.atomSerial != None:
          # Only write if coordinate available!
          fileLines.append("%d %.2f" % (chemShift.atomSerial,chemShift.value))
      else:
        
        fileLine = ""
        
        if len(chemShift.assignments) > 1:
          fileLine = "("

        fileLine+= "|".join(chemShift.assignments)

        if len(chemShift.assignments) > 1:
          fileLine += ")"

        fileLine += " %.2f" % chemShift.value
        
        fileLines.append(fileLine)
    
    #
    # Now write file itself
    #
    
    totalShifts = len(fileLines)     
          
    if use_fout:
      fout = use_fout
    else:
      fout = open(self.name,'w')

    fout.write("CS_VALUES %d" % totalShifts)
    fout.write(self.newline)

    for fileLine in fileLines:
      fout.write(fileLine)
      fout.write(self.newline)

    #if not useSerial:
    #  fout.write("END")
    #  fout.write(self.newline)
      
    if not use_fout:
      fout.close()    

class CosmosChemShift:

  def __init__(self,value,atomSerial = None):
  
    self.value = returnFloat(value)
    self.atomSerial = atomSerial
    
    self.assignments = []
    
  def addAssignment(self,atomName,residueName,chainCode,seqCode):
    
    assignment = "%s_%s_%s_%d" % (atomName.strip(),residueName.strip().upper(),chainCode.upper().strip(),seqCode)
    
    self.assignments.append(assignment)

    #(self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    #self.molCode = defaultMolCode


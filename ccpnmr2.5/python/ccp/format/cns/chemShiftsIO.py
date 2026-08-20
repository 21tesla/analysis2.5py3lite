"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for CNS/ARIA chemical shift file

Copyright (C) 2006 Wim Vranken (European Bioinformatics Institute)

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
from ccp.format.cns.generalIO import CnsGenericFile
from memops.universal.Io import getTopDirectory

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################
    
class CnsChemShiftFile(CnsGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.chemShifts = []

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s chemical shift list %s" % (self.format,self.name)

    fin = open(self.name, 'rU')
    lines = fin.readlines()
    fin.close()
    
    numLines = len(lines)
    errorMessages = []

    for line in lines:

      if self.patt['exclamation'].search(line) or self.patt['emptyline'].search(line):
        pass
      
      else:
        
        chemShiftStoreSearch = self.patt[self.format + 'ChemShiftStore'].search(line)
        seqCodeSearch = self.patt[self.format + 'SeqCode'].search(line)
        atomNameSearch = self.patt[self.format + 'AtomName'].search(line)
        
        if seqCodeSearch and atomNameSearch and chemShiftStoreSearch:
          seqCode = seqCodeSearch.group(1)
          atomName = atomNameSearch.group(1)
          chemShift = chemShiftStoreSearch.group(1)

          self.chemShifts.append(CnsChemShift(chemShift,atomName,seqCode,self.defaultMolCode))
        
        else:
          errorMessages.append("  Error: could not read the following %s chemical shift file line:" % (self.format,))
          errorMessages.append(">>>",line)
    
    #
    # Check if valid file.
    #
    
    fileReadOk = True
    
    if len(errorMessages) > numLines * 0.1:
      print "  Error: file not valid - too many errors."
      fileReadOk = False
    else:
      for errorMessage in errorMessages:
        print errorMessage
    
    return fileReadOk
  
  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing %s chemical shift list %s" % (self.format,self.name)


    fout = open(self.name,'w')

    #
    # Write out chem shifts
    #

    for chemShift in self.chemShifts:

      fout.write("do ( store1 = %.3f )  ( resid %d and name %s )" % (
                      
                       chemShift.value,
                      chemShift.seqCode,
                    	 chemShift.atomName))
                       
      fout.write(self.newline)

class CnsChemShift:

  def __init__(self,value,atomName,seqCode,defaultMolCode):
  
    self.value = returnFloat(value)
    self.atomName = atomName
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.molCode = defaultMolCode

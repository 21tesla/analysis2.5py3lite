"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for XEasy chemical shift files

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

import os

# Import general functions
from memops.universal.Util import returnFloat, returnFloats
from memops.universal.Util import returnInt
from ccp.format.xeasy.generalIO import XEasyGenericFile
from memops.universal.Io import getTopDirectory

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class XEasyChemShiftFile(XEasyGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.chemShifts = []

  def read(self,verbose=False,readAll=False):

    if verbose == 1:
      print "Reading xeasy chemical shift list %s" % self.name

    fin = open(self.name, 'rU')

    line = fin.readline()
    # Read rest file
    while line:

      if not self.patt['hash'].search(line):
      
        # Get rid of trailing hash info if present
        if line.count('#'):
          line = line[:line.index('#')]
      
        cols = line.split()

        if len(cols) == 5:
          atomSerial = cols[0]
          chemShift = cols[1]
          chemShiftError = cols[2]
          atomName = cols[3]
          seqCode = cols[4]

          # TODO: could keep these in and remove at importChemShifts stage...
          if chemShift != '999.000':
            chemShift = XEasyChemShift(self,atomSerial,chemShift,chemShiftError,atomName,seqCode,verbose=False)
            
            if not None in (chemShift.atomSerial,chemShift.value,chemShift.valueError):
              self.chemShifts.append(chemShift)
              
          elif readAll:
            # In case need all atom serials...
            chemShift = XEasyChemShift(self,atomSerial,None,None,atomName,seqCode,verbose=False)  
            self.chemShifts.append(chemShift)          

      line = fin.readline()

    fin.close()
    
    #
    # Crash if nothing read - wrong kind of file
    #
    
    if not self.chemShifts:
      
      raise self.FileParseError("Invalid chemical shift file - no values were read.")

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing xeasy chemical shift list %s" % self.name


    fout = open(self.name,'w')

    #
    # Write three letter codes (one per line), with seqIndex
    #

    for chemShift in self.chemShifts:

      fout.write("%4d %7.3f %5.3f %-5s%4d" % (chemShift.atomSerial,
                                              chemShift.value,
                                              chemShift.valueError,
                                              chemShift.atomName,
                                              chemShift.seqCode))
      fout.write(self.newline)

class XEasyChemShift:

  def __init__(self,parent,atomSerial,value,valueError,atomName,seqCode,verbose=False):
  
    self.atomSerial = returnInt(atomSerial,default=None,verbose=verbose)
    self.value = returnFloat(value,default=None,verbose=verbose)
    self.valueError = returnFloat(valueError,default=None,verbose=verbose)
    self.atomName = atomName
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.molCode = parent.defaultMolCode

###################
# Main of program #
###################

if __name__ == "__main__":  
                                                                                                            
  files = [['../reference/xeasy/3d.1.chemShifts','local/cs.test1']
          ]

  for (inFile,outFile) in files:
    
    chemShiftFile = XEasyChemShiftFile(os.path.join(getTopDirectory(),inFile))

    chemShiftFile.read(verbose = 1)

    for chemShift in chemShiftFile.chemShifts:
      print chemShift.seqCode, chemShift.atomName, chemShift.value


    chemShiftFile.name = outFile

    chemShiftFile.write(verbose = 1)

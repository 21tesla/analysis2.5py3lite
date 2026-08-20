#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for nmrView chemical shift file

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

import os,string

# Import general functions
from memops.universal.Util import returnFloat, returnFloats
from memops.universal.Util import returnInt
from ccp.format.nmrView.generalIO import NmrViewGenericFile
from memops.universal.Io import getTopDirectory

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################
    
class NmrViewChemShiftFile(NmrViewGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.chemShifts = []
    self.seqNameSep = '.'

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading NmrView chemical shift list %s" % self.name

    fin = open(self.name, 'rU')

    line = fin.readline()
    # Read rest file
    while line:

      if not self.patt['hash'].search(line):   # TODO: is this or ! the nmrView comment line?
        cols = line.split()

        if len(cols) != 0:
          (seqCode,atomName) = string.split(cols[0],self.seqNameSep)
          chemShift = cols[1]
          unknownCode = cols[2]

          self.chemShifts.append(NmrViewChemShift(chemShift,atomName,seqCode,self.defaultMolCode,unknownCode = unknownCode))

      line = fin.readline()

    fin.close()

  def write(self,verbose = 0):

    if verbose:
      print "Writing NmrView chemical shift list %s" % self.name


    fout = open(self.name,'w')

    #
    # Write out chem shifts
    #

    for chemShift in self.chemShifts:

      fout.write("%3d%s%-5s %8.3f %d" % (chemShift.seqCode,
                    	 self.seqNameSep,
                    	 chemShift.atomName,
                    	 chemShift.value,
                    	 chemShift.unknownCode))
      fout.write(self.newline)

class NmrViewChemShift:

  def __init__(self,value,atomName,seqCode,defaultMolCode,unknownCode = 0):
  
    self.value = returnFloat(value)
    self.atomName = atomName
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.molCode = defaultMolCode
    self.unknownCode = returnInt(unknownCode)    # TODO: what is this one for?

###################
# Main of program #
###################

if __name__ == "__main__":  
                                                      
  files = [['../reference/nmrView/list.chemicalShifts','local/cs.test1']
          ]

  
  for (inFile,outFile) in files:
    
    chemShiftFile = NmrViewChemShiftFile(os.path.join(getTopDirectory(),inFile))

    chemShiftFile.read(verbose = 1)
    
    for chemShift in chemShiftFile.chemShifts:
      print chemShift.seqCode, chemShift.atomName, chemShift.value


    chemShiftFile.name = outFile
  
    chemShiftFile.write(verbose = 1)

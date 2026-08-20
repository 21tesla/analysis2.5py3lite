#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for Pronto chemical shift files

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

# Import general functions
from memops.universal.Util import returnFloat, returnFloats
from memops.universal.Util import returnInt
from ccp.format.pronto.generalIO import ProntoGenericFile
from memops.universal.Io import getTopDirectory

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class ProntoChemShiftFile(ProntoGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.chemShifts = []
    
  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading Pronto chemical shift list %s" % self.name

    fin = open(self.name, 'rU')

    firstAtom = 'HN'
    secondAtom = 'HA'

    line = fin.readline()

    while line:

      if self.patt['hash'].search(line) or self.patt['emptyline'].search(line):
        line = fin.readline()
        continue


      cols = line.split()

      #
      # Header line
      #  

      if cols[0] == 'Spin':
        firstAtom = cols[2]
        secondAtom = cols[3]

      #
      # Assume all else is line with values on it
      #

      else:
        
        if not self.patt['emptyline'].search(line[0:3]):
        
          spinSys = cols.pop(0)
          spinSys = spinSys[:-1]

          #
          # Read spin systems... don't handle for conversion to data model though.
          #

          if self.patt['emptyline'].search(line[len(spinSys)+1:len(spinSys)+9]):
            resLabel = None
            seqCode = None
            spinSystemId = returnInt(spinSys)

          else:
            (resLabel,seqCode) = cols.pop(0).split(self.resLabelSep)
            spinSystemId = None

          #
          # Is there a firstAtom value listed? Assuming only one value only...
          #

          if not self.patt['emptyline'].search(line[12:19]):
            value = cols.pop(0)
            if not value.count('N/A'):
              self.chemShifts.append(ProntoChemShift(value,firstAtom,seqCode,spinSystemId,resLabel,self.defaultMolCode))

          #
          # Is there a secondAtom value listed? Can be multiple!
          #

          if not self.patt['emptyline'].search(line[20:34]):
            values = line[20:34].split(self.valueSep)

            if len(values) == 1:
              value = cols.pop(0)
              if not value.count('N/A'):
                self.chemShifts.append(ProntoChemShift(value,secondAtom,seqCode,spinSystemId,resLabel,self.defaultMolCode))

            else:
              for i in range(0,len(values)):
                cols.pop(0)
                value = values[i]
                if not value.count('N/A'):
                  self.chemShifts.append(ProntoChemShift(value,secondAtom + str(i+2),seqCode,spinSystemId,resLabel,self.defaultMolCode))

        #
        # Other atoms
        #

        otherAtomsString = ' '.join(cols)

        if otherAtomsString != '':
        
          otherAtoms = otherAtomsString.split(self.valueSep)

          for otherAtom in otherAtoms:
            if otherAtom:
            
              (atomName,value) = otherAtom.split()
              atomName = atomName[:-1]

              if not value.count('N/A'):
                self.chemShifts.append(ProntoChemShift(value,atomName,seqCode,spinSystemId,resLabel,self.defaultMolCode))


      line = fin.readline()

    fin.close()

  def write(self,verbose = 0):
    
    print "NOT WORKING YET"
    return

    if verbose == 1:
      print "Writing Pronto chemical shift list %s" % self.name


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

class ProntoChemShift:

  def __init__(self,value,atomName,seqCode,spinSystemId,resLabel,defaultMolCode):
  
    self.value = returnFloat(value)
    self.atomName = atomName
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.molCode = defaultMolCode
    self.resLabel = resLabel
    
    self.spinSystemId = spinSystemId


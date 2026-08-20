#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

writeBmrbChemShiftDeposition.py: Writes a BMRB chemical shift deposition file

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


from ccpnmr.format.general.Constants import newline

from ccpnmr.format.general.Util import getResName

# Get format classes
from ccpnmr.format.converters.NmrStarFormat import NmrStarFormat

#############
# Functions #
#############

def writeBmrbChemShiftDeposition(guiParent,chain,chemShiftList,outFile):

  proj = chain.root

  format = NmrStarFormat(proj,guiParent,verbose = 1)

  fout = open(outFile,'w')

  resSeq = newline
  
  for residue in chain.sortedResidues():
    
    code1Letter = residue.molResidue.chemComp.getCode1Letter()
    
    if code1Letter:
      resSeq += code1Letter
    else:
      resSeq += 'X'
    
    chopOff = len(resSeq) / 41.0
    
    if chopOff == int(chopOff):
      resSeq += newline
    
  
  fout.write("_Mol_residue_sequence" + newline)
  fout.write(";")
  fout.write(resSeq + newline)
  fout.write(";" + (newline * 2))
  
  fout.write("""
loop_
        _Residue_seq_code
        _Residue_author_seq_code
        _Residue_label""")
  
  fout.write(newline * 2)  
  
  count = 0
  
  for residue in chain.sortedResidues():
  
    count += 1
    
    fout.write("%-3d %-3d %3s" % (residue.seqId, residue.seqCode, residue.ccpCode))
    
    if count / 5.0 == int(count / 5.0):
      fout.write(newline)
    else:
      fout.write("     ")
    
  fout.write("stop_")
  fout.write(newline)  
  
  fout.write("""
###################################################################
#      Chemical Shift Ambiguity Code Definitions                  #
#                                                                 #
#    Codes            Definition                                  #
#                                                                 #
#      1             Unique                                       #
#      2             Ambiguity of geminal atoms or geminal methyl #
#                         proton groups                           #
#      3             Aromatic atoms on opposite sides of the ring #
#                        (e.g. Tyr HE1 and HE2 protons)           #
#      4             Intraresidue ambiguities (e.g. Lys HG and    #
#                         HD protons)                             #
#      5             Interresidue ambiguities (Lys 12 vs. Lys 27) #
#      9             Ambiguous, specific ambiguity not defined    #
#                                                                 #
###################################################################


# INSTRUCTIONS
# 1) Replace the @-signs with appropriate values.
# 2) Text comments concerning the assignments can be
#    supplied in the full deposition.
# 3) Feel free to add or delete rows to the table as needed.
#    The row numbers (_Atom_shift_assign_ID values) will be
#    re-assigned to sequential values by BMRB


# The atom table chosen for this sequence is:

loop_
     _Atom_shift_assign_ID
     _Residue_seq_code
     _Residue_label
     _Atom_name
     _Atom_type
     _Chem_shift_value
     _Chem_shift_value_error
     _Chem_shift_ambiguity_code

#
#Atom   Residue
#shift  Seq     Residue Atom    Atom    Shift/  Error/  Ambiguity
#assign code    Label   Name    Type    ppm     ppm     Code
#---------------------------------------------------------------
#""")
  
  fout.write(newline * 2)
  
  # Using code originally written for XEasy output

  format.getFormatNamingSystemName = getFormatNamingSystemName
  format.writeShifts(outFile,measurementList = chemShiftList,chains = [chain],noWrite = True,compressResonances = False)

  counter = 1
  for chemShift in format.measurementFile.chemShifts:
  
    if chemShift.atomName[-1] == '*':
      atomName = chemShift.atomName[:-1]
    else:
      atomName = chemShift.atomName
       
    fout.write("%-6d   %-5d  %-4s    %-5s   %-5s   %-7.3f %-7.3f %d" % (

	   counter,
	   chemShift.seqCode,
	   chemShift.resLabel,
	   atomName,
	   chemShift.atomType,
	   chemShift.value,
	   chemShift.valueError,
	   chemShift.ambCode
	)
                  )
                  
    counter += 1

    fout.write(newline)

  fout.write(newline)
  fout.write("stop_")
  fout.write(newline*2)
  
  fout.close()
  
  return True

def getFormatNamingSystemName():

  return 'IUPAC'

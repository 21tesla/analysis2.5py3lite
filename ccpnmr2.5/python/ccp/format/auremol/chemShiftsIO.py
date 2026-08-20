#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for auremol chemical shift values

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
from ccp.format.auremol.generalIO import AuremolFile, AuremolMetaFile, AuremolCompFile

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################
    
class AuremolChemShiftFile(AuremolFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.chemShifts = []

  def read(self,metaFile = None, compFile = None, verbose = 0):
  
    compFiles = {}

    if os.path.exists(self.name):
    
      if not metaFile:
    
        metaFile = AuremolMetaFile(self.name)
        metaFile.read(verbose = verbose)
      
      if compFile and os.path.exists(compFile) and len(metaFile.compounds) == 1:
        
        # This only works if one single compound...
        
        compFile = AuremolCompFile(compFile)
        compFile.read(verbose = verbose)
        compFiles[metaFile.compounds[0].label] = compFile
        
      else:
      
        for compound in metaFile.compounds:
          if os.path.exists(compound.compFile):
            compFile = AuremolCompFile(compound.compFile)
            compFile.read(verbose = verbose)
            compFiles[compound.label] = compFile
    
    #
    # Now set the shifts... only works if compfiles available!!
    #
    
    if compFiles:
      
      for compound in metaFile.compounds:
      
        label = compound.label
        compFile = compFiles[label]
        
        for shift in compound.shifts:
          atomKey = shift.atomKey
          
          if compFile.atomDict.has_key(atomKey):

            atomName = compFile.atomDict[atomKey].atomName
            resLabel = compFile.atomDict[atomKey].resLabel

            self.chemShifts.append(AuremolChemShift(shift.value,
                                                    shift.error,
                                                    shift.ambCode,
                                                    shift.atomSetName,
                                                    atomName,
                                                    shift.seqCode,
                                                    label))
          else:
          
            print "  Error: missing atomKey %s for compound %s!" % (atomKey,label)

class AuremolChemShift:

  def __init__(self,value,error,ambCode,atomSetName,atomName,seqCode,molCode):
  
    self.value = value
    self.error = error
    self.ambCode = ambCode
    self.atomSetName = atomSetName
    self.atomName = atomName
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.molCode = molCode


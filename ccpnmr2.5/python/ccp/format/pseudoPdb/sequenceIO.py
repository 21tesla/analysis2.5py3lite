"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for PseudoPDB sequence (coordinate) files (PDB-like but with no header and undefined atom naming)

Copyright (C) 2006-2009 Wim Vranken (European Bioinformatics Institute)

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

import string

from memops.universal.Util import returnInt

from ccp.format.general.Constants import defaultMolCode, defaultSeqInsertCode

from ccp.format.pseudoPdb.coordinatesIO import PseudoPdbCoordinateFile
from ccp.format.pseudoPdb.generalIO import PseudoPdbGenericFile
from ccp.format.general.formatIO import Sequence, SequenceElement

#####################
# Class definitions #
#####################
      
class PseudoPdbSequenceFile(PseudoPdbGenericFile):

  def initialize(self):
  
    self.sequences = []

  def read(self, coordinateFile = None, ignoreResNames = None, verbose = 0):
    
    if not coordinateFile:
      coordinateFile = PseudoPdbCoordinateFile(self.name)
      coordinateFile.read(maxNum = 1)

    seqCode = "XXXXX"
    seqInsertCode = defaultSeqInsertCode
    chainId = -9999
    residueName = ""
    
    modelNums = coordinateFile.modelCoordinates.keys()
    modelNums.sort()

    for coordinate in coordinateFile.modelCoordinates[modelNums[0]]:

      if chainId != coordinate.chainId:

        #
        # New sequence
        #
        
        compoundName = None
        
        for coordChain in coordinateFile.chains:
          if coordChain.chainId == coordinate.chainId and coordChain.compoundName:
            compoundName = coordChain.compoundName

        self.sequences.append(PseudoPdbSequence(chainCode = coordinate.chainId))
        chainId = coordinate.chainId
        seqCode = 'XXXXX'

      if seqCode != coordinate.seqCode or seqInsertCode != coordinate.insertionCode:

        #
        # New residue/item
        # 

        seqCode = coordinate.seqCode
        seqInsertCode = coordinate.insertionCode
        residueName = coordinate.resName

        self.sequences[-1].elements.append(PseudoPdbSequenceElement(str(seqCode) + seqInsertCode,residueName,compoundName = compoundName))
      
      #
      # Keep track of atom names...
      #
      
      self.sequences[-1].elements[-1].addAtomName(coordinate.atomName)
  
class PseudoPdbSequence(Sequence):

  def setFormatSpecific(self,*args,**keywds):
  
    if not self.molName and self.chainCode:
      self.molName = self.chainCode

    self.compoundName = keywds.get('compoundName')
    #self.compoundName = None
    #if keywds.has_key('compoundName'):
    #  self.compoundName = compoundName
  
PseudoPdbSequenceElement = SequenceElement

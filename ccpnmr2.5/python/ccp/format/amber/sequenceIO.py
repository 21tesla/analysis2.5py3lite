"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for AMBER sequence (coordinate) files

Copyright (C) 2010 Wim Vranken (European Bioinformatics Institute)

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

import os, string

from memops.universal.Util import returnInt
from memops.universal.Io import getTopDirectory

from ccp.format.general.Constants import defaultMolCode, defaultSeqInsertCode, bioPolymerCodes

from ccp.format.amber.coordinatesIO import AmberCoordinateFile
from ccp.format.amber.generalIO import AmberGenericFile
from ccp.format.general.formatIO import Sequence, SequenceElement

#####################
# Class definitions #
#####################
      
class AmberSequenceFile(AmberGenericFile):

  def initialize(self):
  
    self.sequences = []
    
    # Do some moltype mappings for standard stuff
    self.molTypeMappings = {'protein': {}}
    
    for resName in bioPolymerCodes['protein'][1]:
      self.molTypeMappings['protein'][resName] = resName
      
    for molType in ('DNA','RNA'):
      self.molTypeMappings[molType] = {}
      for resName in 'ACGTU':
        if molType == 'DNA':
          resName = 'D' + resName
          ccpResName = resName
        elif molType == 'RNA':
          ccpResName = resName
          resName = 'R' + resName
        
        for position in ('','3','5'):
          posResName = resName + position
          
          self.molTypeMappings[molType][posResName] = ccpResName
        
  def read(self, coordinateFile = None, verbose = 0):
    
    if not coordinateFile:
      coordinateFile = AmberCoordinateFile(self.name)
      coordinateFile.read(maxNum = 1)

    seqInsertCode = defaultSeqInsertCode

    modelNums = coordinateFile.modelCoordinates.keys()
    modelNums.sort()
    
    coordinates = coordinateFile.modelCoordinates[modelNums[0]][:]
    
    coordChains = coordinateFile.chains
    
    for coordChain in coordChains:
      
      #
      # New chain
      #
      
      self.sequences.append(AmberSequence(chainCode = coordChain.chainId))
        
      seqCode = None

      for coordinate in coordinates[:]:
      
        if coordinate.chainId == coordChain.chainId:

          #
          # New residue
          #

          if seqCode != coordinate.seqCode or seqInsertCode != coordinate.insertionCode:

            self.sequences[-1].elements.append(AmberSequenceElement(str(coordinate.seqCode) + coordinate.insertionCode,coordinate.resName,self))

            seqCode = coordinate.seqCode
            seqInsertCode = coordinate.insertionCode

          #
          # Keep track of atom names...
          #

          self.sequences[-1].elements[-1].addAtomName(coordinate.atomName)
          
          #
          # Clean up coordinate list to speed up reading (if multiple chains at least)
          #
          
          coordinates.pop(coordinates.index(coordinate))

class AmberSequence(Sequence):

  def setFormatSpecific(self,*args,**keywds):
  
    if not self.molName and self.chainCode:
      self.molName = self.chainCode

class AmberSequenceElement(SequenceElement):

  def setResidueCode(self,*args):
  
    parent = args[1]
    
    residueType = None
    code3Letter = args[0].upper()
    
    # Look if standard polymer, and reset code3Letter if so.
    for molType in parent.molTypeMappings.keys():
      if code3Letter in parent.molTypeMappings[molType]:
        residueType = molType
        code3Letter = parent.molTypeMappings[molType][code3Letter]
        break
    
    self.residueType = residueType
    self.code3Letter = code3Letter
    self.origResName = args[0]
  

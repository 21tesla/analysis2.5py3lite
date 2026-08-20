"""
======================COPYRIGHT/LICENSE START==========================

MMCIFFormat.py: Contains functions specific to mmCIF conversions.

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

import copy

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccp.format.general.Constants import defaultSeqInsertCode

#
# Add some information to IOkeywords...
#
#IOkeywords = copy.deepcopy(IOkeywords)

class MMCIFFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'mmCif'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    
    self.getCoordinates = self.getCoordinatesGeneric

    #self.getChemComps = self.getChemCompsGeneric

  #
  # Deviations from generic import stuff
  #
      
  def getSequenceSetFormatSpecificReadKeywds(self):  
   
    readKeywds = {}

    if self.ignoreResNames:
      readKeywds['ignoreResNames'] = self.ignoreResNames
    
    return readKeywds

  def getCoordinatesSetFormatSpecificReadKeywds(self):
   
    readKeywds = {}
    
    if self.ignoreResNames:
      readKeywds['ignoreResNames'] = self.ignoreResNames
    
    return readKeywds
    
  def getCoordinatesSetSequenceFile(self):
    
    self.setSequenceFileClass()
    self.sequenceFile = self.SequenceFileClass(self.fileName)

    readKeywds = self.getSequenceSetFormatSpecificReadKeywds()
    self.sequenceFile.read(mmCif = self.coordinateFile.mmCif, **readKeywds)

  #
  # Code different from DataFormat
  #

  def duplicateChain(self):
  
    #
    # All chain info should be contained in mmCIF files...
    #
    
    return False

  def getSpecificPdbOneLetterCode(self,currentChainCode):
    
    pdbOneLetterCode = currentChainCode
    
    if self.sequence:
      pdbOneLetterCode = self.sequence.chainCode
      
    # Currently not possible to be more than one letter, so reset to nothing.
    # This restriction will likely change in the future.
    if len(pdbOneLetterCode) > 1:
      pdbOneLetterCode = ' '
  
    return pdbOneLetterCode

  def getPresetChainMapping(self,chainList):
  
    mappingChainDict = {}
    
    for chain in chainList:
    
      mappingChainDict[chain] = (chain.code,1)
    
    return mappingChainDict

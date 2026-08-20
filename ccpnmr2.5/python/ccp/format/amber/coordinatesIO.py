"""
======================COPYRIGHT/LICENSE START==========================

coordinatesIO.py: I/O for AMBER coordinate files

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

# Import general functions
from memops.universal.Util import returnInt, returnFloat
from ccp.format.amber.generalIO import AmberGenericFile

from ccp.format.general.Util import getRegularExpressions
from ccp.format.general.Constants import defaultMolCode, defaultSeqInsertCode

#####################
# Class definitions #
#####################

class AmberCoordinateFile(AmberGenericFile):

  def initialize(self):
  
    self.modelCoordinates = {}
    
    self.patt.update(getRegularExpressions('pdb'))
    
    self.chains = []
    self.serialToCoord = {}
    
    self.chainIds = 'ABCDEFGHIJKLMNO' #PQRSTUVWXYZ' # Limited to 26!

  def read(self, maxNum=999, verbose=False):
  
    """
    Reads in PDB-style coordinate files exported from AMBER only.
    """
    
    if verbose:
      print "  Reading %s coordinate file %s..." % (self.format,self.name)

    #
    # Read file
    #

    fin = open(self.name, 'rU')
    lines = fin.readlines()
    fin.close()
    
    #
    # Get data out. Only one model possible?
    #
    
    modelNum = 1
    chainIdIndex = 0
    
    self.modelCoordinates[modelNum] = []
    self.chains.append(AmberChain(self.getChainId(chainIdIndex)))
    
    lastLineChainAdd = False
    
    for line in lines:

      if line.count("ATOM"):

        serial = line[4:11]
        atomName = line[11:17]
        resName = line[17:20]        
        seqCode = line[20:26]
        
        x = line[30:38]
        y = line[38:46]
        z = line[46:54]
        
        # This is hardset to distinguish between chains
        chainId = self.getChainId(chainIdIndex)
        
        coordinate = AmberCoordinate(serial,chainId,atomName,resName,seqCode,x,y,z)
        
        self.modelCoordinates[modelNum].append(coordinate)
        
        if modelNum == 1:
          self.serialToCoord[coordinate.serial] = coordinate
        
        lastLineChainAdd = False
            
      elif line.count("TER"):
        
        #
        # New chain - only for first model
        #
        
        chainIdIndex += 1
        if modelNum == 1:
          self.chains.append(AmberChain(self.getChainId(chainIdIndex)))
          lastLineChainAdd = True

      elif line.count("END"):
        
        #
        # New model
        #
        
        modelNum += 1
        
        if modelNum > maxNum:
          modelNum -= 1
          break
          
        chainIdIndex = 0
        self.modelCoordinates[modelNum] = []
        
        if lastLineChainAdd:
          self.chains.pop(-1)
          lastLineChainAdd = False
     
    #
    # Clean up
    #
        
    if lastLineChainAdd:
      self.chains.pop(-1)
      
    if not self.modelCoordinates[modelNum]:
      del self.modelCoordinates[modelNum]

  def getChainId(self,chainIdIndex):
  
    chainIdsLen = len(self.chainIds)
  
    if chainIdIndex < chainIdsLen:
      chainId = self.chainIds[chainIdIndex]
    else:
      index1 = (chainIdIndex / chainIdsLen) - 1
      index2 = chainIdIndex % chainIdsLen
      
      chainId = self.chainIds[index1] + self.chainIds[index2]
      
    return chainId      

class AmberChain:

  def __init__(self,chainId):

    self.chainId = chainId
    #self.refChainId = defaultMolCode

class AmberCoordinate:

  def __init__(self,serial,chainId,atomName,resName,seqCode,x,y,z, insertionCode = defaultSeqInsertCode):
  
    self.serial = returnInt(serial)
    self.atomName = atomName.strip()
    self.resName = resName.strip()
    self.seqCode = returnInt(seqCode)

    self.x = returnFloat(x)
    self.y = returnFloat(y)
    self.z = returnFloat(z)
    self.insertionCode = insertionCode

    self.chainId = chainId
    #self.refChainId = defaultMolCode
      
    # Use refChainId - is defaultMolCode. Then set chainID from A-Z, ..

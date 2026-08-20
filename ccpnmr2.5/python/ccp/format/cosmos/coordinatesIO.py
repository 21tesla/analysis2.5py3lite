"""
======================COPYRIGHT/LICENSE START==========================

coordinatesIO.py: I/O for Cosmos coordinate files

Copyright (C) 2008 Wim Vranken (European Bioinformatics Institute)

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
from ccp.format.cosmos.generalIO import CosmosGenericFile

from ccp.format.general.Constants import defaultMolCode, defaultSeqInsertCode

#####################
# Class definitions #
#####################

class CosmosCoordinateFile(CosmosGenericFile):

  def initialize(self):
  
    self.modelCoordinates = {}
    
    self.chains = []

  #def read(self,maxNum = 999):

  #  pass

  def write(self,version = '003',verbose = False, use_fout = None):

    if verbose == 1:
      print "Writing %s coordinate file %s" % (self.format,self.name)
     
    if use_fout:
      fout = use_fout
    else:
      fout = open(self.name,'w')
    
    atomLineFormat = "%s %d %7.5f %7.5f %7.5f"
    
    if not use_fout:
      fout.write("$COO%s" % version)
      fout.write(self.newline)

    fout.write("CELL  1.0000  1.0000  1.0000  90.000  90.000  90.000  0.0000  0.0000  0.0000")
    fout.write(self.newline)
    
    fout.write("SPGROUP P1")
    fout.write(self.newline)

    # TODO can handle multiple models?
    # TODO can handle multiple chains?
    modelNums = self.modelCoordinates.keys()
    modelNums.sort()
    
    if len(modelNums) > 1:
      print "  Warning: only writing out first model!"
    
    for modelNum in modelNums[:1]:
        
      #oldChainId = self.modelCoordinates[modelNum][0].chainId
      
      fout.write("ATOMS %d" % len(self.modelCoordinates[modelNum]))
      fout.write(self.newline)
      
      for coord in self.modelCoordinates[modelNum]:
      
        #if oldChainId != coord.chainId:
        #  oldChainId = coord.chainId
        #  fout.write("TER" + self.newline)
      
        atomId = "%s_%s" % (coord.atomName,coord.resName)
        if coord.chainCode:
          atomId += '_%s' % coord.chainCode
        atomId += '_%s' % coord.seqCode
      
        fout.write(atomLineFormat % (atomId,coord.atomNumber,coord.x,coord.y,coord.z))
        
        fout.write(self.newline)
                
    if not use_fout:
      fout.close()    

class CosmosChain:

  def __init__(self):
    
    # Not supported
    self.chainId = defaultMolCode

class CosmosCoordinate:

  def __init__(self,serial,chainCode,seqCode,atomName,atomNumber,resName,x,y,z):
  
    self.serial = returnInt(serial)
    self.atomName = atomName.strip()
    self.atomNumber = returnInt(atomNumber)
    self.resName = resName.strip().upper()
    self.seqCode = seqCode
    self.chainCode = chainCode.strip()

    self.x = returnFloat(x)
    self.y = returnFloat(y)
    self.z = returnFloat(z)
    
    # Not supported
    self.insertionCode = defaultSeqInsertCode
    self.chainId = defaultMolCode

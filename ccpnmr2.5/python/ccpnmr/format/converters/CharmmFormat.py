#!/usr/bin/pythonc

"""
======================COPYRIGHT/LICENSE START==========================

CharmmFormat.py: Contains functions specific to CHARMM conversions.

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

import copy

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccp.format.general.Constants import defaultSeqInsertCode

from ccpnmr.format.general.Constants import distanceConstraintDefaultLowerLimit
from ccpnmr.format.general.Constants import defaultMolCode

#
# Add some information to IOkeywords...
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['writeRdcConstraints']['refFrameSegId'] = (defaultMolCode,False,'The segment Id for the reference frame.')
IOkeywords['writeRdcConstraints']['refFrameResId'] = (999,False,'The residue code for the reference frame.')

class CharmmFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'charmm'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric

    self.getCoordinates = self.getCoordinatesGeneric
    self.createCoordinateFile = self.createCoordinateFileGeneric 

  #
  # Deviations from generic import stuff
  #
  
  
  #
  # Functions different to default functions in DataFormat
  #

  def duplicateChain(self):
  
    #
    # All chain info should be contained in coordinate file...
    #
    
    return 0   

  def getPresetChainMapping(self,chainList):
  
    return self.getMultiChainFormatPresetChainMapping(chainList)    

  def setRawCoordinate(self):
  
    if not self.coordinateFile.modelCoordinates.has_key(self.modelId):

      self.coordinateFile.modelCoordinates[self.modelId] = []
    
    chemCompVar = self.residue.chemCompVar
    chemComp = chemCompVar.chemComp
    cifCode = self.findChemCompVarSysName('CIF',chemCompVar)
    
    if self.useOriginalData and self.origResLabel != None:
      resName = self.origResLabel
    elif cifCode:
      resName = cifCode
    elif chemComp.code3Letter:
      resName = chemComp.code3Letter
    else:
      resName = chemComp.ccpCode
      
    
    if self.useOriginalData and self.origChainCode != None:
      chainCode = self.origChainCode
    else:
      chainCode = self.exportChainCode
  
    if self.useOriginalData and self.origSeqCode != None:
      seqCode = self.origSeqCode
    else:
      seqCode = self.seqCode
 
    if self.useOriginalData and self.origSeqInsertCode != None:
      seqInsertCode = self.origSeqInsertCode
    else:
      seqInsertCode = defaultSeqInsertCode
 
    if self.useOriginalData and self.origAtomName != None:
      atomName = self.origAtomName
    else:
      atomName = self.atomName
 
    #
    # TODO: setting chainId but is really segId... does charmm use the chainID position?
    #

    modelCoordinate = self.coordinatesIO.CharmmCoordinate(self.coordinateSerial,
                                                  atomName,
                                                  resName,
                                                  chainCode,
                                                  seqCode,
                                                  self.residue.seqId,
                                                  self.x,
                                                  self.y,
                                                  self.z)
      
    self.coordinateFile.modelCoordinates[self.modelId].append(modelCoordinate)

  def getExportChainCode(self,exportChainCode,chain):
    
    # Returns original chain code - this format uses segId so OK if longer than 1 char
    
    return exportChainCode

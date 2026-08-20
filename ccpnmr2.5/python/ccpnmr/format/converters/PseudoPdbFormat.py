#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

PseudoPdbFormat.py: Contains functions specific to PseudoPdb conversions (PDB-like files without the proper header and with possibly different atom naming)

Copyright (C) 2006 Wim Vranken (European Bioinformatics Institute)

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

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readCoordinates']['swapFirstNumberAtom'] = (False,False,'Click if you want to change atom names that start with a number (e.g. 2HB -> HB2).')

class PseudoPdbFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'pseudoPdb'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    
    self.getCoordinates = self.getCoordinatesGeneric
    self.createCoordinateFile = self.createCoordinateFileGeneric 

  #
  # Deviations from generic import stuff
  #

  def getCoordinatesSetFormatSpecificReadKeywds(self):
    
    readKeywds = {}
    
    if self.swapFirstNumberAtom:
      readKeywds['swapFirstNumberAtom'] = self.swapFirstNumberAtom
    
    if self.ignoreResNames:
      readKeywds['ignoreResNames'] = self.ignoreResNames
    
    return readKeywds

  def getSequenceSetFormatSpecificReadKeywds(self):  
   
    readKeywds = {}

    if self.ignoreResNames:
      readKeywds['ignoreResNames'] = self.ignoreResNames
    
    return readKeywds
  
  #
  # Functions different to default functions in DataFormat
  #

  def setRawCoordinate(self):
  
    if not self.coordinateFile.modelCoordinates.has_key(self.modelId):

      self.coordinateFile.modelCoordinates[self.modelId] = []
    
    chemCompVar = self.residue.chemCompVar
    chemComp = chemCompVar.chemComp
    
    preferredResName = self.findChemCompVarSysName(self.namingSystemName,chemCompVar)
    cifCode = self.findChemCompVarSysName('CIF',chemCompVar)
    
    if self.useOriginalData and self.origResLabel != None:
      resName = self.origResLabel
    elif preferredResName:
      resName = preferredResName
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
 
    modelCoordinate = self.coordinatesIO.PseudoPdbCoordinate(
                    
                    self.coordinateFile,
                    self.coordinateSerial,
                    atomName,
                    resName,
                    chainCode,
                    '',
                    seqCode,
                    self.x,
                    self.y,
                    self.z,
                    '',
                    insertionCode = seqInsertCode)
      
    self.coordinateFile.modelCoordinates[self.modelId].append(modelCoordinate)

  def getPresetChainMapping(self,chainList):
  
    return self.getMultiChainFormatPresetChainMapping(chainList)
    

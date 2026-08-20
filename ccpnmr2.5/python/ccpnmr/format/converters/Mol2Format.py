#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

Mol2Format.py: Contains functions specific to TRIPOS .mol2 format conversions.

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

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

class Mol2Format(DataFormat):

  def setFormat(self):
  
    self.format = 'mol2'
    
    self.bondTypeRawToCcpn = {'1': 'single',
                              '2': 'double',
                              '3': 'triple',
                              'am':'singleplanar',
                              'ar':'aromatic'}
    
    self.bondTypeCcpnToRaw = {}
    for rawBondType in self.bondTypeRawToCcpn.keys():
      self.bondTypeCcpnToRaw[self.bondTypeRawToCcpn[rawBondType]] = rawBondType
    
    # This dict is other way around from previous - bit confusing!
    self.molTypeCcpnToRaw = {'protein':      'PROTEIN',
                             'DNA':          'NUCLEIC_ACID',
                             'RNA':          'NUCLEIC_ACID',
                             'DNA/RNA':      'NUCLEIC_ACID',
                             'carbohydrate': 'SACCHARIDE',
                             'other':        'SMALL'}

    self.molTypeRawToCcpn = {}
    for molType in self.molTypeCcpnToRaw.keys():
      self.molTypeRawToCcpn[self.molTypeCcpnToRaw[molType]] = molType

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
#    self.createSequenceFile = self.createSequenceFileGeneric

    self.getChemComps = self.getChemCompsGeneric
    self.createChemCompFile = self.createChemCompFileGeneric
    
  #
  # Deviations from generic import stuff
  #
      

  #
  # Code different from DataFormat
  #
    
  def convertRawBondType(self,rawBondType):
    
    if self.bondTypeRawToCcpn.has_key(rawBondType):
      bondType = self.bondTypeRawToCcpn[rawBondType]
    else:
      bondType = 'void'
    
    return bondType
  
  def convertRawMolType(self,rawMolType):
    
    if self.molTypeRawToCcpn.has_key(rawMolType):
      molType = self.molTypeRawToCcpn[rawMolType]
    else:
      molType = 'other'
    
    return molType

  def setRawChemComp(self):
    
    name = self.chemCompVar.chemComp.ccpCode.upper()
    molType = self.chemCompVar.chemComp.molType
    
    # Not handling 'BIOPOLYMER' raw type in below
    if self.molTypeCcpnToRaw.has_key(molType):
      rawMolType = self.molTypeCcpnToRaw[molType]
    else:
      rawMolType = 'SMALL'
      
    self.rawChemComp = self.chemCompsIO.Mol2ChemComp(self.chemCompFile,name,[],[],{},molType = rawMolType,overallCharge=self.chemCompVar.formalCharge)
    self.chemCompFile.chemComps.append(self.rawChemComp)
    
  def setRawChemCompAtom(self,atomSerial,rawAtomName,x,y,z,rawAtomType,charge,chemCompCode):
      
    atomInfo = {'x': x,'y': y,'z': z,'atomType': rawAtomType,'charge':charge,'substName':chemCompCode.upper(),'substId': 1, 'atomStatus':'****'}
      
    self.rawChemComp.atoms.append(self.chemCompsIO.Mol2Atom(self.rawChemComp,atomSerial,rawAtomName,atomInfo))
   

  def setRawChemCompBond(self,atomSerial1,atomSerial2,bondType):
       
    if self.bondTypeCcpnToRaw.has_key(bondType):
      rawBondType = self.bondTypeCcpnToRaw[bondType]
    else:
      # Default to unknown!
      rawBondType = 'un'
      
    self.rawChemComp.bonds.append(self.chemCompsIO.Mol2Bond(self.rawChemComp,atomSerial1,atomSerial2,rawBondType,''))
        
   

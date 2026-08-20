#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

chemCompsIO.py: I/O for PDB format hetero 'chemical compound' group reading.

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

import os, string

# Import general functions
from memops.universal.Util import returnInt, returnFloat
from memops.universal.Io import getTopDirectory

# Import header reader
from ccp.format.pdb.coordinatesIO import PdbCoordinateFile
from ccp.format.pdb.generalIO import PdbGenericFile

from ccp.format.general.Constants import defaultSeqInsertCode
from ccp.format.general.Util import standardNucleusName

#####################
# Class definitions #
#####################

class PdbChemCompFile(PdbGenericFile):

  def initialize(self):
  
    self.chemComps = []

  def read(self,pdbFile = None,verbose = 0):

    if not pdbFile:
    
      #
      # First read coordinate file (if not passed in)
      #

      pdbFile = PdbCoordinateFile(self.name)
      pdbFile.read(maxNum = 1, readHeader = 0)
      
    
    #
    # Assuming one chemComp only!
    #
    
    self.chemComps.append(PdbChemComp(self))
    
    firstModel = pdbFile.modelCoordinates.keys()[0]
    bonds = []
    coorHandled = {}
    
    for coor in pdbFile.modelCoordinates[firstModel]:
      
      self.chemComps[-1].atoms.append(
      
            PdbAtom(self.chemComps[-1],
                    coor.serial,
                    coor.atomName,
                    coor.x,
                    coor.y,
                    coor.z,
                    coor.atomType))

      for bondType in coor.bonds.keys():
        for bondedCoor in coor.bonds[bondType]:
          if bondedCoor not in coorHandled.keys():
            bonds.append((bondType,coor,bondedCoor))

      coorHandled[coor] = self.chemComps[-1].atoms[-1]
    
    for (bondType,coor1,coor2) in bonds:
    
      atom1 = coorHandled[coor1]
      atom2 = coorHandled[coor2]
       
      self.chemComps[-1].bonds.append(PdbBond(self.chemComps[-1],atom1,atom2,bondType))


class PdbChemComp:

  def __init__(self,parent):

    self.parent = parent

    self.atoms = []
    self.bonds = []

class PdbAtom:

  def __init__(self,parent,atomSerial,name,x,y,z,atomType):
  
    self.parent = parent
    
    self.serial = returnInt(atomSerial)
    self.name = string.strip(name)
    
    self.x = returnFloat(x,verbose = 0)
    self.y = returnFloat(y,verbose = 0)
    self.z = returnFloat(z,verbose = 0)
    
    self.atomType = string.strip(atomType)
    if not self.atomType:
      for char in self.name:
        if char not in '0123456789':
          self.atomType = char
          break

class PdbBond:

  def __init__(self,parent,atom1,atom2,bondType):
  
    self.parent = parent
        
    self.atom1 = atom1
    self.atomSerial1 = parent.atoms.index(atom1) + 1
    self.atom2 = atom2
    self.atomSerial2 = parent.atoms.index(atom2) + 1
    
    self.bondType = 'single' # TODO: can we do better here?

###################
# Main of program #
###################

if __name__ == "__main__":  

  files = [
           '../../reference/pdb/cyclic.pdb',
           '../../reference/pdb/mcys.pdb'
          ]

  
  for inFile in files:
    
    pdbFile = PdbChemCompFile(os.path.join(getTopDirectory(),inFile))

    pdbFile.read(verbose = 1)
    
    for chemComp in pdbFile.chemComps:
      for atom in chemComp.atoms:
        print atom.atomType,
      print
      for bond in chemComp.bonds:
        print bond.bondType
        
    #(path,baseName) = os.path.split(inFile)
    
    #pdbFile.name = ("local/%s" % baseName)
    
    #pdbFile.write(verbose = 1)

"""
======================COPYRIGHT/LICENSE START==========================

AquaFormat.py: Contains functions specific to Aqua conversions.

Copyright (C) 2007 Wim Vranken (European Bioinformatics Institute)

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

import string, copy

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccp.general.Util import getResonancesFromPairwiseConstraintItem

from ccpnmr.format.general.Constants import distanceConstraintDefaultLowerLimit
from ccpnmr.format.general.Util import getNameInfo

class AquaFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'aqua'

  def setGenericImports(self):
    
    self.getConstraints = self.getConstraintsGeneric
    self.createConstraintFile = self.createConstraintFileGeneric
    
  #
  # Functions different to default functions in DataFormat
  #

  def setRawDistanceConstraint(self):
  
    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraint.serial))

    self.rawConstraint = self.constraintFile.constraints[-1]
    
    lowerLimit = self.constraint.lowerLimit
    
    if lowerLimit < 0.0:
      lowerLimit = distanceConstraintDefaultLowerLimit
      
    self.rawConstraint.setLowerDistance(lowerLimit)    
    self.rawConstraint.setUpperDistance(self.constraint.upperLimit)
    
  def setRawDistanceConstraintItem(self):
    
    self.rawConstraint.items.append(self.rawConstraintItemClass())
    self.rawConstraintItem = self.rawConstraint.items[-1]
    
  def setRawDistanceConstraintItemMembers(self):
    
    itemResonances = getResonancesFromPairwiseConstraintItem(self.item)
    
    for i in range(0,2):
      
      resonance = itemResonances[i]

      resLabel = self.getResonanceResLabel(resonance)
      
      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])
      
      self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass(chainCode,seqCode,resLabel,atomName))

  def getResonanceResLabel(self,resonance):
  
   resLabel = None
  
   chemCompVar = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue.chemCompVar
   namingSystem = chemCompVar.chemComp.findFirstNamingSystem(name = self.namingSystemName)
   
   chemCompSysName = chemCompVar.findFirstSpecificSysName(namingSystem = namingSystem)
   
   if not chemCompSysName:
     chemCompSysName = chemCompVar.findFirstChemCompSysName(namingSystem = namingSystem)
     
   if chemCompSysName:
     resLabel = chemCompSysName.sysName
   else:
     resLabel = chemCompVar.chemComp.ccpCode.upper()
   
   return resLabel

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)

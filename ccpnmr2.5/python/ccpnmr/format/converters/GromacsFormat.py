"""
======================COPYRIGHT/LICENSE START==========================

GromacsFormat.py: Contains functions specific to Gromacs conversions.

Copyright (C) 2012 Wim Vranken (Vrije Universiteit Brussel)

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

- contact Wim Vranken (wvranken@vub.ac.be)
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

# Set self.individualAtoms to True, always


import copy

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccp.general.Util import getResonancesFromPairwiseConstraintItem

from ccpnmr.format.general.Constants import atomSerial_kw
from ccpnmr.format.general.Constants import distanceConstraintDefaultLowerLimit

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['writeDistanceConstraints']['individualAtoms'] = (True,False,"Always has to be True for GROMACS.")
IOkeywords['writeDistanceConstraints']['atomSerialFormat'] = ('gromos',False,"Format from where the atom serials for constraint output come from.")
IOkeywords['writeDistanceConstraints']['upperLimitStartLinear'] = (0.1,False,"Add this distance (in nm!!) to upperlimit to define point where penalty function becomes linear.")

class GromacsFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'gromacs'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    #self.getConstraints = self.getConstraintsGeneric
    self.createConstraintFile = self.createConstraintFileGeneric
    
  #
  # Functions different to default functions in DataFormat
  #
  
  def createConstraintFileFormatSpecific(self):
  
    self.resonanceToAtomSerial = {}
    
    # Hard-force this setting, overwritten with above copy. approach (hacky!)
    self.individualAtoms = True

  def setRawDistanceConstraint(self):

    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraint.serial))

    self.rawConstraint = self.constraintFile.constraints[-1]
    
    # Distances are in nm in GROMACS! Divide angstrom by 10...
    lowerLimit = self.constraint.lowerLimit / 10.0
    
    if lowerLimit < 0.0:
      lowerLimit = distanceConstraintDefaultLowerLimit / 10.0
    
    upperLimit1 = self.constraint.upperLimit / 10.0
    upperLimit2 = upperLimit1 + self.upperLimitStartLinear
    
    self.rawConstraint.setDistances(lowerLimit,upperLimit1,upperLimit2)
        
    self.atomSerialCombs = []
    
  def setRawDistanceConstraintItem(self):
  
    # Dealt with at member stage!
    pass
    
  def setRawDistanceConstraintItemMembers(self):
    
    itemResonances = getResonancesFromPairwiseConstraintItem(self.item)
    
    atomSerialList = []
    
    for i in range(2):
    
      resonance = itemResonances[i]
      
      if resonance in self.resonanceToAtomSerial.keys():
        atomSerials = self.resonanceToAtomSerial[resonance]
      
      else:
      
        atomSerials = []
        resonanceToAtoms = self.resonanceToAtoms[resonance]
        
        for resonanceToAtom in resonanceToAtoms:
        
          atom = resonanceToAtom.getAtom()
          atomSerial = self.getAtomSerial(atom)
          
          if atomSerial == None:
            print "  Error: no atom serial for atom '{}.{}.{}, aborting this constraint.'".format(atom.residue.chain.code,atom.residue.seqCode,atom.name)
            atomSerials = []
            break
          
          atomSerials.append(atomSerial)
      
      atomSerialList.append(atomSerials)
      atomSerialList[-1].sort()
    
    #
    # Now make all combinations...
    #

    for atomSerial1 in atomSerialList[0]:
      for atomSerial2 in atomSerialList[1]:
      
        atomSerialComb = (atomSerial1,atomSerial2)
        
        if atomSerialComb in self.atomSerialCombs:
          continue
        else:
          self.atomSerialCombs.append(atomSerialComb)
      
        self.rawConstraint.items.append(self.rawConstraintItemClass())
        self.rawConstraint.items[-1].members.append(self.rawConstraintItemMemberClass(atomSerial1))
        self.rawConstraint.items[-1].members.append(self.rawConstraintItemMemberClass(atomSerial2))
    
  def getAtomSerial(self,atom):
  
    atomSerial = None
    
    atomSerialAppData = atom.findFirstApplicationData(application=self.atomSerialFormat,keyword=atomSerial_kw)

    if atomSerialAppData:
      atomSerial = int(atomSerialAppData.value)
      
    return atomSerial

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)

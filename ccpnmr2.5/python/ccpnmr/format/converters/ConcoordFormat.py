"""
======================COPYRIGHT/LICENSE START==========================

ConcoordFormat.py: Contains functions specific to Concoord conversions.

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

import string

from memops.universal.Io import joinPath, splitPath 

from ccp.general.Util import findChemAtomSysName, findChemAtomSetSysName
from ccp.general.Util import getResonancesFromPairwiseConstraintItem

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccpnmr.format.general.Util import getNameInfo, getResName

from ccpnmr.format.general.Constants import allResidueAtoms_kw

class ConcoordFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'concoord'

  def setGenericImports(self):
    
    self.getConstraints = self.getConstraintsGeneric
    self.createConstraintFile = self.createConstraintFileGeneric

  #
  # Functions different to default functions in DataFormat
  #

  def setRawDistanceConstraint(self):
  
    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraintFile,self.constraint.serial))

    self.rawConstraint = self.constraintFile.constraints[-1]
    
    self.rawConstraint.upperDist = self.constraint.upperLimit
    self.rawConstraint.lowerDist = self.constraint.lowerLimit
    
    self.originalItems = []
  
  def setRawDistanceConstraintInfo(self):

    #
    # Check whether a pseudo correction applies... only necessary if more than one item
    # after reorganization (need to combine stereospecific HB1/HB2 as well...)
    #
    
    if len(self.rawConstraint.items) > 1:
      if len(self.originalItems) > 1:
        checkOriginalItems = []
        checkAtomNames = [[],[]]
      
        for item in self.originalItems:
          checkOriginalItem = []
          
          for i in range(0,2):
            atomName = item[i][2]
            
            if atomName not in checkAtomNames[i]:
              checkAtomNames[i].append(atomName)
            
            if atomName[-1] in '123':
              atomName = atomName[:-1] + '#'
            
            checkOriginalItem.append([item[i][0],item[i][1],atomName])

          if checkOriginalItems.count(checkOriginalItem) == 0:
            checkOriginalItems.append(checkOriginalItem)
        
        #
        # Make sure it is a pseudo group...
        #
        
        if len(checkOriginalItems) == 1:
          for i in range(0,2):
            if len(checkAtomNames[i]) == 1:
              checkOriginalItems[0][i][2] = checkAtomNames[i][0]
        
      else:      
        checkOriginalItems = self.originalItems
            
      #
      # Now set the correction...
      #
      
      if len(checkOriginalItems) == 1:
        pseudoCorrection = 0.0

        for i in range(0,2):
          chain = checkOriginalItems[0][i][0]
          seqId = checkOriginalItems[0][i][1]
          atomName = checkOriginalItems[0][i][2]
          
          if atomName[-1] in '#*%+':
            
            searchAtomName = atomName[:-1] + '*'
            residue = chain.findFirstResidue(seqId = seqId)
            resLabel = residue.molResidue.chemComp.ccpCode
            
            chemCompVar = residue.molResidue.chemCompVar
            
            chemAtomSet = chemCompVar.findFirstChemAtomSet(name = searchAtomName)
            
            if chemAtomSet:
              pseudoCorrection += chemAtomSet.distCorr
                
            else:
              print "  Error: no chemAtomSet with name %s for ccpCode %s." % (searchAtomName,resLabel)
          
        self.rawConstraint.upperPseudoDist = self.constraint.upperLimit + pseudoCorrection
        
      else:
        self.rawConstraint.upperPseudoDist = 999.000
    
    else:
      self.rawConstraint.upperPseudoDist = self.constraint.upperLimit
    
  def setRawDistanceConstraintItem(self):
    
    #
    # Handling almost everything on item member level - need to decompress some restraints!
    #
    
    self.originalItems.append(self.ccpInfo)
    
  def setRawDistanceConstraintItemMembers(self):
  
    #
    # Expand the information... is necessary for Concoord (only handles individual atoms)
    #

    atomNameLists = []
    
    for i in range(0,2):
    
      atomNameLists.append([])
    
      resonanceToAtom = self.ccpInfo[i][3]
      chemAtomOrigName = self.ccpInfo[i][2]
      
      if resonanceToAtom.chemAtomSet:
        chemAtomOrSet = resonanceToAtom.chemAtomSet
      else:
        chemAtomOrSet = resonanceToAtom.chemAtom
      
      if chemAtomOrSet:

        chemCompVar = resonanceToAtom.getResidue().chemCompVar
        namingSystem = chemCompVar.chemComp.findFirstNamingSystem(name = self.namingSystemName)
        if namingSystem:
          chemAtomSysName = findChemAtomSysName(namingSystem,{'sysName': chemAtomOrigName})
        else:
          chemAtomSysName = None
        
        if chemAtomSysName:
          chemAtom = chemCompVar.findFirstChemAtom(name = chemAtomSysName.atomName, subType = chemAtomSysName.atomSubType)
        else:
          chemAtom = chemCompVar.findFirstChemAtom(name = chemAtomOrigName)
        
        if not chemAtom:
        
          chemAtomSetSysName = findChemAtomSetSysName(namingSystem,{'sysName': chemAtomOrigName})
          
          if chemAtomSetSysName:
            chemAtomSet = chemCompVar.findFirstChemAtomSet(name = chemAtomSetSysName.atomName, subType = chemAtomSetSysName.atomSubType)
          else:
            chemAtomSet = chemCompVar.findFirstChemAtomSet(name = chemAtomOrigName)
          
          if chemAtomSet:
          
            (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])

            if chemAtomSet.chemAtomSets:
              chemAtomSets = list(chemAtomSet.chemAtomSets)
            else:
              chemAtomSets = [chemAtomSet]
              
            for chemAtomSet in chemAtomSets:

              chemAtomSysNames = self.getChemAtomSysNamesFromSet(chemAtomSet, chemCompVar, findIupac = 1)

              for chemAtomSysName in chemAtomSysNames:

                if chemAtomSysName:
                  atomNameLists[-1].append(getResName(chainCode,seqCode,chemAtomSysName,seqInsertCode = seqInsertCode))
                else:
                  print "  Error: problems with decompressing Concoord name %s - sysNames missing for %s" % (self.resSetNames[i],chemAtomSet.name)
                
          else:
            print "  Error: can't decompress Concoord name %s" % self.resSetNames[i]
       
      #
      # If no names use default
      #
      
      if not atomNameLists[-1]:
        atomNameLists[-1].append(self.resSetNames[i])
    
    
    #
    # Now write them one by one...
    #
    
    itemResonances = getResonancesFromPairwiseConstraintItem(self.item)
    
    for resSetName0 in atomNameLists[0]:
    
      # TODO: bit risky this one! What if original names (eg GLU-) are being used?!
      resonance = itemResonances[0]
      resLabel0 = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue.molResidue.chemComp.ccpCode.upper()
      
      (chainCode0,seqCode0,spinSystemId0,seqInsertCode0,atomName0) = getNameInfo(resSetName0)

      if atomName0 == allResidueAtoms_kw:
        print "  Error: not handling residue level constraints for %s." % self.format
        continue
      
      for resSetName1 in atomNameLists[1]:

        # TODO: bit risky this one! What if original names (eg GLU-) are being used?!
        resonance = itemResonances[1]
        resLabel1 = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue.molResidue.chemComp.ccpCode.upper()

        (chainCode1,seqCode1,spinSystemId1,seqInsertCode1,atomName1) = getNameInfo(resSetName1)
         
        if atomName1 == allResidueAtoms_kw:
          print "  Error: not handling residue level constraints for %s." % self.format
          continue

        #
        # Now set all items and members...
        #

        self.rawConstraint.items.append(self.rawConstraintItemClass())
        self.rawConstraintItem = self.rawConstraint.items[-1]
      
        self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass())
        self.rawConstraintItem.members[-1].setInfo(chainCode0,seqCode0,resLabel0,atomName0)
        self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass())
        self.rawConstraintItem.members[-1].setInfo(chainCode1,seqCode1,resLabel1,atomName1)

  def getPresetChainMapping(self,chainList):
  
    return self.getMultiChainFormatPresetChainMapping(chainList)

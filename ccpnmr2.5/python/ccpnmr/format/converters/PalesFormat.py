#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

ModuleFormat.py: Contains functions specific to PALES program conversions.

Copyright (C) 2010 Wim Vranken (European Bioinformatics Institute)
                   Brian Smith (University of Glasgow)

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
- MSD website (http://www.ebi.ac.uk/msd/)

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

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccp.format.general.Constants import defaultSeqInsertCode

from ccpnmr.format.general.Constants import distanceConstraintDefaultLowerLimit
from ccpnmr.format.general.Util import getNameInfo

from ccp.general.Util import getResonancesFromPairwiseConstraintItem

class PalesFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'pales'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getConstraints = self.getConstraintsGeneric
    self.createConstraintFile = self.createConstraintFileGeneric

  #
  # Functions different to default functions in DataFormat
  #
  
  def setRawRdcConstraint(self):
  
    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraint.serial))

    self.rawConstraint = self.constraintFile.constraints[-1]
    
    #
    # Have to get some sensible values out
    #
    
    lowerLimit = self.constraint.lowerLimit    
    upperLimit = self.constraint.upperLimit        
    targetValue = self.constraint.targetValue
    if targetValue == None:
      targetValue = (lowerLimit + upperLimit) / 2
        
    self.rawConstraint.setRdcData(targetValue,error = self.constraint.error)
   
  def setRawRdcConstraintItem(self):
    
    self.rawConstraint.items.append(self.rawConstraintItemClass())
    self.rawConstraintItem = self.rawConstraint.items[-1]
    
  def setRawRdcConstraintItemMembers(self):
    
    itemResonances = getResonancesFromPairwiseConstraintItem(self.item)

    for i in range(0,2):
            
      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])
      
      resLabel = self.getResonanceResLabel(itemResonances[i])

      self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass(chainCode,seqCode,resLabel,atomName))

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)

  def getResonanceResLabel(self,resonance):
   # copied from DyanaFormat.py (also in AquaFormat.py - scope for common function in Util?)
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

  def createConstraintFileFormatSpecific(self):
    
    # Here bit tricky because first have to get number of chains, and format
    # doesn't really handle it. I'd do this by having a specific
    # chain keyword for writing RDC constraints for PALES, if it's not set
    # and there's only one chain, all good, otherwise throw an
    # error or ask the user for the relevant chain, then filter out
    # constraints (that's all in place in case it's needed, just ask). Wim

    # I do something similar already in DataFormat.py, but it happens at a later stage in the workflow.
    chains = []
    for constraint in self.constraintList.constraints:
      for item in constraint.sortedItems():
        for fixedResonance in item.resonances:
          if fixedResonance.resonanceSet:
            refAtom = fixedResonance.resonanceSet.findFirstAtomSet().findFirstAtom()
	    if refAtom.residue.chain not in chains:
              chains.append(refAtom.residue.chain)
    
    seqStrings = []
    for chain in chains:
       # Filter out most ligands, hopefully.
       if len(chain.residues) > 1:
         seqStrings.append(chain.molecule.seqString)
    
    if len(seqStrings) > 1:
      print "Warning: multiple sequences present, picking first one"
    
    print 'in', self
    self.writeKeywds['oneLetterSequence'] = seqStrings[0]

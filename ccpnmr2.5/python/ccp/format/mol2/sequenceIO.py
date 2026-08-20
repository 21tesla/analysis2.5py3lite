#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for Tripos Mol2 sequence (coordinate) file

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

import os, string

# Import general functions
from memops.universal.Util import returnInt
from memops.universal.Io import getTopDirectory

# Import header reader
from ccp.format.mol2.generalIO import Mol2GenericFile
from ccp.format.general.formatIO import Sequence, SequenceElement

from ccp.format.general.Constants import bioPolymerCodes, chainCodeString
from ccp.format.general.Constants import defaultSeqInsertCode, bondTypes

#####################
# Class definitions #
#####################  

class Mol2SequenceFile(Mol2GenericFile):
  """
  Information on file level
  """
  def initialize(self):

    self.sequences = []
    
    self.bioPolymerCodes = bioPolymerCodes
    
  def read(self, mol2File = None, verbose = 0):
    
    if not mol2File:
      self.readGeneric()
    else:
      self.infoDict = mol2File.infoDict
    
    if verbose == 1:
      print "Getting sequence/molecule information."
    
    #
    # TODO here do things very similar to PDB.
    #
    # IN Mol2Format, have to load the chemcomp if there's no match
    #                try to match if 'similar' residue (e.g. PHE1, TYR1, ...)
    #
    
    
    
    # OR should I try to be smart, see if there is a 'protein' backbone, and
    # then use modified amino acids instead?
    
    #
    # Molecule information
    #
    
    molName = self.infoDict['MOLECULE']['molName']
    molType = self.infoDict['MOLECULE']['molType']
    
    if self.infoDict['MOLECULE'].has_key('molComment'):
      molDetails = self.infoDict['MOLECULE']['molComment']
    else:
      molDetails = None
    
    totalAtoms = len(self.infoDict['ATOM'])
    
    if totalAtoms != self.infoDict['MOLECULE']['numAtoms']:
      print "  Warning: possible atom numbers mismatch. %d reported, %d counted." % (self.infoDict['MOLECULE']['numAtoms'],totalAtoms)
    
    #
    # Get the chains (?) Or not worth it?
    #
    # TODO: if multiple chains, have to set bonds between them as well!! Currently NOT done!
    #
    
    chainIds = {}
    
    for substInfo in self.infoDict['SUBSTRUCTURE']:
      
      if substInfo.has_key('chain'):
        chainId = substInfo['chain']
      else:
        chainId = self.defaultMolCode
      
      if not chainId in chainIds.keys():
        chainIds[chainId] = []
      
      chainIds[chainId].append(substInfo)
    
    #
    # Create molecules...
    #
    
    chainIdList =  chainIds.keys()
    chainIdList.sort()
    
    for chainId in chainIdList:
    
      self.sequences.append(Mol2Sequence(molName = molName)) # TODO multiple molNames possible?!?
      
      currentSequence = self.sequences[-1]
      currentSequence.setDetails(molDetails)
      currentSequence.setChainCode(chainId)
      
      if molType == 'PROTEIN':
        molType = 'protein'
      elif molType == 'NUCLEIC_ACID':
        molType = 'DNA/RNA'
      elif molType == 'SACCHARIDE':
        molType = 'carbohydrate'
      else:
        molType = 'other'
      
      # TODO: make decision: use protein backbone if available and add other groups as Het (other molecule)
      # groups, or make modified amino acids that one can connect (on the fly?) or make one giant molecule?
      
      # The PDB philosophy is to make the backbone if one can... but how can I match this to
      # the modified amino acids that I have info for? A translation dictionary that they can
      # set themselves? But what about the atom mapping then? And all the extra chemComps?
      
      # Probably the best solution is to make a 'local' chemComp directory with the information as in the file -
      # there's enough info from these mol2 files to do that.
      
      allSubstAtoms = {}

      for substInfo in chainIds[chainId]:
      
        seqCode = substInfo['substId'] # This should do...
        resName = substInfo['substName']
        
        (substAtoms,atomNames,atomIds) = self.findSubstAtoms(totalAtoms,substInfo['rootAtomId'],resName)
        
        # Try to set if it is protein...
        if resName[:3] in bioPolymerCodes['protein'][1]:
          residueType = 'protein'
        elif resName[:3] in bioPolymerCodes['DNA'][1] + bioPolymerCodes['RNA'][1]:
          residueType = 'DNA/RNA'
        else:
          hasProteinBbAtoms = True
          for bbatom in ('C','N','O','CA'):
            if bbatom not in atomNames:
              hasProteinBbAtoms = False
              break
          
          if hasProteinBbAtoms:
            residueType = 'protein'
          else:
            residueType = molType
        
        if residueType == 'DNA/RNA':
          for atomName in atomNames:
            if atomName.count('2HO') or atomName.count('HO2'):
              residueType = 'RNA'
              break
            elif atomName.count('2H2'):
              residueType = 'DNA'
              break

        # TODO WARNING: resName is NOT a 3-letter code!!!
        currentSequence.elements.append(Mol2SequenceElement(seqCode,resName,residueType,atomNames))
        
        #
        # Track atoms for inter-substructure bonds...
        #
        
        allSubstAtoms[currentSequence.elements[-1]] = atomIds
    
      #
      # Track inter-substructure bonds
      #
      
      totalBonds = len(self.infoDict['BOND'])
      circularInfo = []

      for i in range(totalBonds):
        bondInfo = self.infoDict['BOND'][i]
        
        atomId1 = bondInfo['atom1Id']
        atomId2 = bondInfo['atom2Id']
        
        substs = [None,None]
        
        for i in range(0,2):
          atomId = (atomId1,atomId2)[i]
          for seqEl in allSubstAtoms.keys():
            if atomId in allSubstAtoms[seqEl].keys():
              substs[i] = seqEl
              break
        
        if substs[0] != substs[1]:
          substs[0].setBond(bondInfo['bondType'],allSubstAtoms[substs[0]][atomId1]['atomName'],substs[1],allSubstAtoms[substs[1]][atomId2]['atomName'])
          linkPair = [substs[0].seqCode,substs[1].seqCode]
          linkPair.sort()
          if not linkPair in circularInfo:
            circularInfo.append(linkPair)
      
      #
      # Check if circular
      #
      
      if circularInfo:
      
        circularInfo.sort()
        
        while (self.removeUnconnectedLinkPairs(circularInfo)):
          pass
          
        if circularInfo:
          currentSequence.setCircular()

  def removeUnconnectedLinkPairs(self,circularInfo):

    infoPopped = False

    # First remove terminal ones (one side no connection)
    for linkPair in circularInfo[:]:
      connection = [0,0]
      for otherLinkPair in circularInfo:
        if linkPair == otherLinkPair:
          continue
        for i in range(0,2):
          if linkPair[i] in otherLinkPair:
            connection[i] += 1

      if 0 in connection:
        circularInfo.pop(circularInfo.index(linkPair))
        infoPopped = True

    return infoPopped
    
class Mol2Sequence(Sequence):

  def setCircular(self):
    self.isCircular = True

  def setDetails(self,value):
    self.details = value

  def setChainCode(self,value):
    self.chainCode = value
    
class Mol2SequenceElement(SequenceElement):

  def setFormatSpecific(self,*args):
    
    self.residueType = args[1]
    self.atomNames = args[2]

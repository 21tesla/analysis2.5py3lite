"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for mmCIF sequence (coordinate) file

Copyright (C) 2010 Wim Vranken (European Bioinformatics Institute)

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
from memops.universal.Util import returnInt

# Import header reader
from ccp.format.mmCif.generalIO import MMCIFFile
from ccp.format.general.formatIO import Sequence, SequenceElement

#from ccp.format.general.Constants import bioPolymerCodes, chainCodeString
from ccp.format.general.Constants import defaultSeqInsertCode

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################  
      
class MMCIFSequenceFile(MMCIFFile):
  """
  Information on file level
  """
  def initialize(self):
  
    # NOTE: todo this overwrites initialize from generalIo.py class!

    self.sequences = []

    self.chainCodeToSequence = {}
    
    self.seqElDict = {}
    
    # This dictionary tracks the old to the new chain/sequence codes, if any, for molecules
    self.residueMappings = {}
    
  def read(self, mmCif=None, ignoreResNames=None, version=None, verbose=False):

    # Not necessary at the moment
    self.version = version
    self.ignoreResNames = ignoreResNames
        
    #
    # Read mmCIF file and set up information...
    #
    
    if not mmCif:
      self.readGeneric(verbose=verbose)
    else:
      self.mmCif = mmCif
    
    sequenceInfo = self.mmCif.getSequenceInfo()
    bondInfo = self.mmCif.getBondInfo()
    self.chemCompInfo = self.mmCif.getChemCompInfo()
    
    self.coordinateInfo = self.mmCif.getCoordinateInfo()
    models = self.coordinateInfo.keys()
    models.sort()
    self.coordinateRefModel = models[0]
   
    entityIds = sequenceInfo.keys()
    entityIds.sort()
    
    #
    # First defined polymers...
    #   
    
    for entityId in entityIds[:]:
    
      (moleculeType,moleculeName,polymerType,chainCodes,residueInfo) = sequenceInfo[entityId]
      
      # Note: I'm always treating protein and DNA/RNA polymers as stand-alone chains!
      if moleculeType == 'polymer' and polymerType:
        
        entityIds.pop(entityIds.index(entityId))
        
        chainCodes = sequenceInfo[entityId][-1].keys()
        chainCodes.sort()
        
        if polymerType.count('polypep'):
          molType = 'protein'
        elif polymerType.count('deoxy'):
          molType = 'DNA'
        elif polymerType.count('polyribo'):
          molType = 'RNA'
        else:
          print "  Warning: unknown mmCIF molecular type %s." % (polymerType)
          molType = None
                
        for chainCode in chainCodes:
        
          currentSequence = self.addSequence(moleculeName,chainCode,moleculeType)
          
          # TODO: what is originalChainCode here? The PDB??!?!
          
          seqElInfo = sequenceInfo[entityId][-1][chainCode]
          
          for seqEl in seqElInfo:
          
            (cifSeqCode,cifResLabel, pdbResLabel, pdbChainCode, pdbSeqCode, pdbInsertionCode, authResLabel, authSeqCode, isHetero) = seqEl
          
            self.setSequenceElement(currentSequence,pdbChainCode,pdbSeqCode,pdbSeqCode,cifResLabel,molType,pdbInsertionCode)

    #
    # Now the rest, also check for bonds in between to join up elements with different CIF chain codes.
    #   
    
    otherChains = []
    
    for entityId in entityIds:
    
      (moleculeType,moleculeName,polymerType,chainCodes,residueInfo) = sequenceInfo[entityId]
      
      chainCodes = sequenceInfo[entityId][-1].keys()
      chainCodes.sort()
      
      for chainCode in chainCodes:
      
        seqElInfo = sequenceInfo[entityId][-1][chainCode]
      
        #
        # Check if *covalently* bonded to anything...
        #
        
        connectedChains = []

        if bondInfo.has_key('covale'):

          for seqEl in seqElInfo:

            (uniqueSeqCode,cifResLabel, pdbResLabel, pdbChainCode, pdbSeqCode, pdbInsertionCode, authResLabel, authSeqCode, isHetero) = seqEl

            for covBondInfo in bondInfo['covale']:
              
              isBonded = False
              for k in range(2):
                if covBondInfo[k][0] == chainCode and covBondInfo[k][4] == pdbSeqCode:
                  isBonded = True
                  otherIndex = not(k)
                  break
              
              if isBonded:
                (otherChainCode, otherResLabel, otherSeqId, otherAuthChainCode, otherAuthSeqId, otherAtomName, otherPdbInsertionCode) = covBondInfo[otherIndex]
                
                # Only get chain, deal with actual bonds later...
                connectedChains.append(otherChainCode)

        #
        # Check if needs to be connected to existing set of chain info
        #
        
        currentChainInfo = None
        
        for otherChainInfo in otherChains:
          isConnected = False
          if chainCode in otherChainInfo['connectedChains']:
            isConnected = True
            
          # If really far apart, doesn't work...
          for otherChainCode in connectedChains:
            if otherChainCode in otherChainInfo['connectedChains']:
              isConnected = True
              connectedChains.append(chainCode)
              break
          
          if isConnected:
            print "Adding info to existing set of chains"
            currentChainInfo = otherChainInfo
            break
            
        if not currentChainInfo:
          otherChains.append({'connectedChains': [chainCode], 'seqElInfo': [], 'currentSequence': None, 'moleculeName': moleculeName, 'moleculeType': moleculeType})
          currentChainInfo = otherChains[-1]
      
        #
        # If bonded only to existing polymer (as from self.chainCodeToSequence chainCode) then connect to that directly
        # and bail
        #
        
        sequenceAdded = False
        
        for connectedChainCode in connectedChains:
          if connectedChainCode in self.chainCodeToSequence.keys():
          
            currentChainInfo['currentSequence'] = self.chainCodeToSequence[connectedChainCode]
            break
            
        #
        # If already earmarked with a sequence, add it to that one...
        #
        
        if currentChainInfo['currentSequence']:
        
          # Note have to reset sequence codes if going to be part of same molecule!!
          #lastSeqCode = currentChainInfo['currentSequence'].getLastSeqCode()

          for i in range(len(seqElInfo)):
            (uniqueSeqCode,cifResLabel, pdbResLabel, pdbChainCode, pdbSeqCode, pdbInsertionCode, authResLabel, authSeqCode, isHetero) = seqElInfo[i]

            self.setSequenceElement(currentSequence,pdbChainCode,pdbSeqCode,pdbSeqCode,cifResLabel,None,pdbInsertionCode) #lastSeqCode + i

        #
        # Otherwise just add, will deal with creation further down.
        #
        
        else:
        
          for connectedChainCode in connectedChains:
            if not connectedChainCode in currentChainInfo['connectedChains']:
              currentChainInfo['connectedChains'].append(connectedChainCode)
          
          currentChainInfo['seqElInfo'].extend(seqElInfo)


    #
    # Now create rest of sequences...
    #
    
    for otherChainInfo in otherChains:
    
      # Don't do if already set!!
      if not otherChainInfo['currentSequence']:
      
        chainCodes = otherChainInfo['connectedChains']
        chainCodes.sort()

        currentSequence = self.addSequence(otherChainInfo['moleculeName'],chainCodes[0],otherChainInfo['moleculeType'],alternateChainCodes=chainCodes[1:])
        
        seqElInfo = otherChainInfo['seqElInfo']
        seqElInfo.sort()

        for seqEl in seqElInfo:

          (uniqueSeqCode,cifResLabel, pdbResLabel, pdbChainCode, pdbSeqCode, pdbInsertionCode, authResLabel, authSeqCode, isHetero) = seqEl

          self.setSequenceElement(currentSequence,pdbChainCode,pdbSeqCode,pdbSeqCode,cifResLabel,None,pdbInsertionCode)

    #
    # Set the bonds on the sequence element level!
    #
    
    mmCifBondTypes = bondInfo.keys()
    mmCifBondTypes.sort()
    
    for mmCifBondType in mmCifBondTypes:
      
      if mmCifBondType == 'disulf':
        bondType = 'disulfide'
      elif mmCifBondType == 'covale':
        bondType = 'covalent' 
      elif mmCifBondType == 'hydrog':
        bondType = 'hydrogen' 
      # TODO: not sure what type is here.
      #elif mmCifBondType == 'salt':
      #  bondType = 'salt'
      else:
        bondType = 'link' 

      #
      # Loop over the bonds
      #

      for covBondInfo in bondInfo[mmCifBondType]:
        seqEls = []
        for k in range(2):
          (cifChainCode, cifResLabel, cifSeqId, authChainCode, authSeqId, atomName, pdbInsertionCode) = covBondInfo[k]
          
          seqKey = (authChainCode,authSeqId,cifResLabel)

          if self.seqElDict.has_key(seqKey):
            seqEl = (self.seqElDict[seqKey],atomName)
          else:
            print seqKey, self.seqElDict.keys()
            print "  Error: could not find sequence element for '%s.%s.%s' to set bond..." % (cifChainCode,cifSeqId,atomName)
            seqEl = None
            
          seqEls.append(seqEl)
          
        #
        # Set the bond, but ONLY if between different sequence elements.
        # Internal chemComp bonds should be covered by E-MSD info!
        #

        if not None in seqEls and seqEls[0] != seqEls[1]:
          for k in range(2):
            atomName = seqEls[k][1]
            seqEl = seqEls[k][0]
            otherAtomName = seqEls[not k][1]
            otherSeqEl = seqEls[not k][0]

            seqEl.setBond(bondType,atomName,otherSeqEl,otherAtomName)

    #
    # Check if circular protein present...
    #

      
    """
   
      #
      # Check if circular protein...
      #

      if len(proteinCircularCoords) == 2 and proteinCircularCoords[0].bonds.has_key('covalent'):
        for bondedCoord in proteinCircularCoords[0].bonds['covalent']:
          if bondedCoord == proteinCircularCoords[1]:
            currentSequence.setCircular()
            break

    """

    
    #
    # TODO: Set secondary structure information? First have to get it from mmCIF!
    #

    """

    #
    # Set the secondary structure information on the sequence and sequence element level...
    #
    
    for secStrucType in self.pdbFile.seqStrucInfo.keys():
      for secStrucInfo in self.pdbFile.seqStrucInfo[secStrucType]:
        
        (serial,ssId,initResInfo,endResInfo,specificInfo) = secStrucInfo
        
        # Assuming sec struc element always within same chain!
        chainCode = initResInfo[1]
        
        sequence = None
        for tSeq in self.sequences:
          if tSeq.chainCode == chainCode:
            sequence = tSeq
            break
        
        if sequence:
          
          seqStrucTypeText = secStrucType.lower()
          
          seqEls = []
          inSecStruc = False
          
          for seqEl in sequence.elements:
            if seqEl.seqCode == initResInfo[2] and seqEl.insertionCode == initResInfo[3]:
              inSecStruc = True
            
            if inSecStruc:
              seqEl.setSecStrucInfo(seqStrucTypeText,serial,specificInfo)
              seqEls.append(seqEl)
          
            if seqEl.seqCode == endResInfo[2] and seqEl.insertionCode == endResInfo[3]:
              break
              
          sequence.setSecondaryStructure(seqStrucTypeText,serial,specificInfo,seqEls)

        else:
          print "  Warning: could not find chain code '%s' for setting secondary structure info." % chainCode

      """
         
  def getResidueAtomNames(self,pdbChainCode,pdbSeqCode):          

    coordinatesInfo = self.coordinateInfo[self.coordinateRefModel]
    
    atomNames = []
    foundCoord = False

    for coordinateInfo in coordinatesInfo[:]:
    
      if coordinateInfo[15] == pdbChainCode and coordinateInfo[13] == pdbSeqCode:
        atomNames.append(coordinateInfo[2])
        foundCoord = True
      elif foundCoord:
        # Assuming that are in order, so once residue finished don't bother looking for more atom names.
        break


    return atomNames
          
  def addSequence(self,moleculeName,chainCode,moleculeType,alternateChainCodes=None):
  
    if chainCode in self.chainCodeToSequence.keys():
      print "  Warning: trying alternate chain code for molecule %s - '%s' already exists..." % (moleculeName,chainCode)
      
      if alternateChainCodes:
        for tmpChainCode in alternateChainCodes:
          if not tmpChainCode in self.chainCodeToSequence.keys():
            chainCode = tmpChainCode
            break
        
    currentSequence = MMCIFSequence(molName=moleculeName, chainCode=chainCode, moleculeType=moleculeType)

    self.chainCodeToSequence[chainCode] = currentSequence

    self.sequences.append(currentSequence)
    
    return currentSequence

  def setSequenceElement(self,currentSequence,pdbChainCode,pdbSeqCode,seqCode,origResLabel,molType,insertionCode):

    resLabel = self.convertResName(origResLabel)

    # Ignore if not necessary
    if self.ignoreResNames and (resLabel in self.ignoreResNames or origResLabel in self.ignoreResNames):
      return

    atomNames = self.getResidueAtomNames(pdbChainCode,pdbSeqCode)
    
    # Set the molType if not yet available
    if not molType:
      molType = 'other'
    
      if self.chemCompInfo.has_key(resLabel):
        mmCifMolType = self.chemCompInfo[resLabel][0]
        if mmCifMolType.count('sacc'):
          # TODO: for the time being, just use 'other' for this
          molType = 'other'
        elif mmCifMolType.count('pept'):
          molType = 'protein'
        elif mmCifMolType.lower().count('non-polymer'):
          molType = 'other'
        else:
          print "  Warning: did not interpret mmCIF molecular type %s..." % mmCifMolType      
      
    # Make sure this is not None
    if not insertionCode:
      insertionCode = defaultSeqInsertCode
    
    currentSequence.elements.append(MMCIFSequenceElement(seqCode,resLabel,molType,insertionCode,atomNames))
    
    seqElKey = (pdbChainCode,pdbSeqCode,resLabel)
    
    self.seqElDict[seqElKey] = currentSequence.elements[-1]
    
    # Set mapping from original to new information
    self.residueMappings[(pdbChainCode,pdbSeqCode)] = (currentSequence.chainCode,seqCode)
          
class MMCIFSequence(Sequence):

  def setCircular(self):
    self.isCircular = True
    
  def getLastSeqCode(self):
  
    lastSeqCode = None
    for seqEl in self.elements:
      if not lastSeqCode or lastSeqCode < seqEl.seqCode:
        lastSeqCode = seqEl.seqCode
    
    return lastSeqCode
    
  def setFormatSpecific(self,*args,**keywds):
  
    self.molType = keywds['moleculeType']

  #def setSecondaryStructure(self,secStrucType,secStrucSerial,specificInfo,seqEls):
  #  
  #  if not self.secStrucInfo.has_key(secStrucType):
  #    self.secStrucInfo[secStrucType] = {}
  #  
  #  self.secStrucInfo[secStrucType][secStrucSerial] = (specificInfo,tuple(seqEls))
 
class MMCIFSequenceElement(SequenceElement):

  def setFormatSpecific(self,*args):
    
    self.residueType = args[1]
    self.insertionCode = args[2]
    self.atomNames = args[3]

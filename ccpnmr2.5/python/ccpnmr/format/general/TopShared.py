
"""
======================COPYRIGHT/LICENSE START==========================

TopShared.py: Shared functions for DataFormat and linkResonances

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

from ccpnmr.format.general.Util import createSelection, getNameInfo
from memops.general.Util import returnMemopsLine

from ccp.api.molecule import MolSystem
from ccpnmr.api import NameMapping

from ccp.format.general.Constants import defaultSeqInsertCode

class TopShared:

  """
  
  This is a class that is shared by converters.DataModel and 
  process.linkResonances.
  
  Several functions are defined here that occur in both. This
  is because the resonance-atom linking for Nmr data and the
  atom-atom linking for coordinate files uses similar functions
  during setup
  
  """
  
  def __init__(self):
    
    #
    # These are all defined by DataFormat or linkResonances
    #
    
    self.project = None
    self.guiParent = None
    self.multiDialog = None
    self.format = None
    self.formatChains = None
    self.forceDefaultChainMapping = None
    self.minimalPrompts = None
    self.messageReporter = None
    self.chainMappingChain = None
    self.allowPopups = True
  
  def selectMolSystem(self, resonanceList = None):
    
    #
    # If minimalPrompts, just select any available molSystem
    #
    
    if self.minimalPrompts and len(self.project.molSystems) == 1:
      self.molSystem = self.project.findFirstMolSystem()
    elif not hasattr(self,'molSystem'):
      self.molSystem = None
      
    #
    # If resonance info is given, try to determine from to-atom link
    #
    
    if resonanceList:
      for resonance in resonanceList:
        resonanceSet = resonance.resonanceSet
        if resonanceSet:
          refAtom = resonanceSet.findFirstAtomSet().findFirstAtom()
          self.molSystem = refAtom.residue.chain.molSystem
          break
        
    #
    # Try to select from list
    #
    
    if not self.molSystem:
          
      (selectionList,selectionDict) = createSelection(self.project.sortedMolSystems())    
      
      if selectionList:
  
        while (not self.molSystem):
  
          interaction = self.multiDialog.SelectionList(
  
                               self.guiParent,
                               selectionList,
                               selectionDict = selectionDict,
                               title = "Project '%s': " % self.project.name + 'Select molecular system',
                               text = 'Existing molecular system codes:',
                               urlFile = 'SelectMolSystem')
      
          if interaction.isSelected:
            self.molSystem = interaction.selection
    
    #
    # If nothing, create new
    #
    
    if not self.molSystem:
    
      molSysExists = 1
      molSysCode = 'My_molecular_system'
      
      while (molSysExists):
    
        molSysCode = self.dataEntry.askString("Molecular system name","Code for new molecular system:",molSysCode,self.guiParent)
        
        if self.project.findFirstMolSystem(name = molSysCode):
          self.messageReporter.showError("Error","Name for molecular system already exists - select another.",parent = self.guiParent)
        
        elif ' ' in molSysCode:
          self.messageReporter.showError("Error","Name for molecular system cannot have spaces.",parent = self.guiParent)
        
        else:
          molSysExists = 0
          
      self.molSystem = MolSystem.MolSystem(self.project,code = molSysCode, name = molSysCode)

  
  def getNameMapping(self,namingSystemName = None, isOriginalImport = False, resetMapping = False):

    #
    # Initialize... reset nameMapping and set namingSystemName
    #
  
    self.nameMapping = None
    
    if not namingSystemName:
       
      namingSystemName = 'Unknown'
      
    if not self.project.currentNameMappingStore:
      nameMappingStore = self.project.findFirstNameMappingStore(name = self.project.name)
      if not nameMappingStore:
        nameMappingStore = self.project.newNameMappingStore(name = self.project.name)
      self.project.currentNameMappingStore = nameMappingStore
      
    #
    # See if there are any matches
    #
    
    nameMappings = self.project.currentNameMappingStore.findAllNameMappings(formatName = self.format, namingSystem = namingSystemName)

    if nameMappings:

      if len(nameMappings) > 1:
      
        print "  Warning: multiple %s %s nameMappings available - picking first one." % (self.format,namingSystemName)

      self.nameMapping = list(nameMappings)[0]

      # TODO: need selection here! Is highly improbable... is this the key?

    else:
          
      if namingSystemName != 'Unknown':
      
        #
        # If nothing found could be 'cos there's an Unknown naming system - use that
        # one and rename
        #
        
        self.nameMapping = self.project.currentNameMappingStore.findFirstNameMapping(formatName = self.format, namingSystem = 'Unknown')
        
        if self.nameMapping:
          self.nameMapping.namingSystem = namingSystemName
          self.nameMapping.isOriginalImport = isOriginalImport
    
    #
    # If resetting mapping, delete existing one...
    #

    if self.nameMapping and resetMapping:
      self.nameMapping.delete()
      self.nameMapping = None
    
    #
    # If there was nothing at all, then create...
    #
  
    if not self.nameMapping:
      self.nameMapping = self.project.currentNameMappingStore.newNameMapping(formatName = self.format, namingSystem = namingSystemName,isOriginalImport = isOriginalImport)

  def setChainLinks(self,forceChainMappings = None,resetMapping = False,verbose = 0):
  
    self.status = True
    
    #
    # Don't bother if no chains...
    #

    if not self.molSystem.chains:
      self.messageReporter.showError("Error","No chains present in molecular system, cannot link CCPN to external format information.",parent = self.guiParent)
      self.status = False
      return self.status
      
    #
    # Delete if available and resetting mapping
    #
    
    if self.nameMapping and resetMapping:
      self.nameMapping.delete()
      self.nameMapping = None
    
    #
    # Set name mapping (if not yet available)
    #
    
    if not self.nameMapping:

      #
      # This also includes selection of correct name mapping (namingSys)
      # TODO: always set isOriginalImport to True?
      #
      
      if hasattr(self,'namingSystemName'):
        namingSystemName = self.namingSystemName
      else:
        namingSystemName = None

      self.getNameMapping(namingSystemName = namingSystemName, isOriginalImport = True, resetMapping = resetMapping)

    #
    # Force chain mapping, this is always tried even if mappings available.
    # Everything has to match or won't work.
    #

    if forceChainMappings:
    
      forceMapping = 1
      
      for forceChainMapping in forceChainMappings:
      
        chainCode = forceChainMapping[0]
        chain = self.molSystem.findFirstChain(code = chainCode)
        
        if not chain:
          print "  Error: forceChainMappings chain code '%s' is not present in molSystem! Aborting..." % chainCode
          forceMapping = 0
          break
        
      if forceMapping:
        
        # TODO: formatLastSeqCode necessary here? Currently not implemented....
        
        for forceChainMapping in forceChainMappings:
      
          chainCode = forceChainMapping[0]
          chain = self.molSystem.findFirstChain(code = chainCode)
          
          formatChainCode = forceChainMapping[1]
          firstSeqId = forceChainMapping[2]
          offset = forceChainMapping[3]
          
          if len(forceChainMapping) == 5:
            seqInsertCode = forceChainMapping[4]
          else:
            seqInsertCode = defaultSeqInsertCode
          
          self.setChainMappingDicts(chain,formatChainCode,firstSeqId = firstSeqId, offset = offset, seqInsertCode = seqInsertCode)

        return self.status
        
    #
    # Check for chain mapping - take from name mapping if available *and* complete
    #
    
    self.setMappingFromChainMappings = False
    
    if self.nameMapping.chainMappings:
    
      deleteFormatChains = []
      setChainMappings = []
      doSet = 1

      for chainMapping in self.nameMapping.chainMappings:
        
        #print "CHAIN MAPPING", chainMapping.chainCode, chainMapping.formatChainCode, chainMapping.formatFirstSeqCode
        
        chain = self.molSystem.findFirstChain(code = chainMapping.chainCode)
                
        if chain:
        
          #
          # Additional check to see whether these codes make sense or not...
          #
          
          residueMatches = 1
          totalResidues = 1
          
          if hasattr(self,'resUnassignedNames'):
            residueMatches = 0
            totalResidues = 0
            residueList = []
            
            formatChainCode = chainMapping.formatChainCode
            formatFirstSeqCode = chainMapping.formatFirstSeqCode
            formatSeqInsertCode = chainMapping.formatSeqInsertCode

            if chainMapping.formatLastSeqCode != None:
              formatLastSeqCode = chainMapping.formatLastSeqCode
            else:
              formatLastSeqCode = formatFirstSeqCode + chain.molecule.seqLength

            for resName in self.resUnassignedNames.keys():
              (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resName,verbose = 0)

              if (seqCode,seqInsertCode) not in residueList:
                residueList.append((seqCode,seqInsertCode))
                totalResidues += 1
              
                if formatChainCode == chainCode and formatSeqInsertCode == seqInsertCode and seqCode >= formatFirstSeqCode and seqCode <= formatLastSeqCode:
                  residueMatches += 1

          if residueMatches == totalResidues:
          
            setChainMappings.append((chain,chainMapping.formatChainCode,chainMapping))
            if not chainMapping.formatChainCode in deleteFormatChains:
              deleteFormatChains.append(chainMapping.formatChainCode)
          
          else:
            
            doSet = 0
            break
          
        else:
        
          doSet = 0
          print "  Error invalid chain codes for this molecular system - doing manual selection"
          self.formatChainCodeDict = {}
          self.ccpChainDict = {}
          break
      
      #
      # Now check whether all is good and set dicts if so...
      #
      
      if doSet:
        
        self.setMappingFromChainMappings = True
 
        for (chain,formatChainCode,chainMapping) in setChainMappings:
          self.setChainMappingDicts(chain,formatChainCode,chainMapping = chainMapping,verbose = verbose)
        
        for formatChainCode in deleteFormatChains:
          #print formatChainCode,
          if self.formatChains.has_key(formatChainCode):
            #print " DELETED",
            del self.formatChains[formatChainCode]
          #print

      #
      # If no format chains left, mapping is set so continue
      #

      if chain and not self.formatChains:
        return self.status
        
    #
    # Only go ahead if one chain in project and one chain reported,
    # otherwise USER INTERACTION to link up!
    #
    
    linkChains = 1
 
    if len(self.molSystem.chains) == 1 and len(self.formatChains.keys()) == 1:
      
      #
      # Do some initialization
      #
      
      formatChainCode = (self.formatChains.keys())[0]
      chain = self.molSystem.sortedChains()[0]
      
      #
      # Now check for the first seqCode in the data model chain that matches the
      # first seqCode of the format.
      #
      
      if self.formatChains[formatChainCode]:

        (seqCode,seqInsertCode) = self.formatChains[formatChainCode][0]
        
        residues = chain.sortedResidues()

        for resPos in range(len(residues)):      
          residue = residues[resPos]

          if residue.seqCode == seqCode and residue.seqInsertCode == seqInsertCode:

            #
            # Additional check: is the number of residues left big enough?
            #

            if len(chain.residues) - resPos >= len(self.formatChains[formatChainCode]):

              #
              # And does the last sequence code match as well?
              #

              residue = chain.sortedResidues()[resPos + len(self.formatChains[formatChainCode]) - 1]
              (seqCode,seqInsertCode) = self.formatChains[formatChainCode][-1]

              if residue.seqCode == seqCode and residue.seqInsertCode == seqInsertCode:

                self.setSeqAndInsertCodeToChainMapping(formatChainCode,self.formatChains[formatChainCode],chain,resPos,verbose = verbose)

                linkChains = 0
                break


    if linkChains == 1:
    
      #
      # User interaction
      # Link chains to chains that were reported for new resonances
      #
      
      interaction = self.multiDialog.ChainLink(self.guiParent,self.formatChains,self.molSystem,forceDefaultChainMapping = self.forceDefaultChainMapping, doInteraction = self.allowPopups)
      
      # 
      # Chain mapping returned?
      #
      
      if hasattr(interaction,'chainDict') and interaction.chainDict != {}:
        
        for chain in interaction.chainDict:
          
          (formatChainCode,firstSeqId,seqAndInsertCodesList) = interaction.chainDict[chain]

          #
          # Quick check to see if possibly seqCodes have to be reset... could be more
          # sophisticated.
          #
          
          nmrRes = chain.findFirstResidue(seqId = firstSeqId)
          
          firstResidueSeqId = chain.sortedResidues()[0].seqId
          
          if self.allowPopups and not self.minimalPrompts and nmrRes.seqCode != seqAndInsertCodesList[0][0]:
            
            if self.messageReporter.showYesNo("Reset sequence codes","Do you want to reset the sequence codes for chain '%s' in the data model to reflect the original ones?" % chain.code):

              #
              # Reset seqCode for nmrResidues
              #
              
              for seqCodeIndex in range(0,len(seqAndInsertCodesList)):
                (seqCode,seqInsertCode) = seqAndInsertCodesList[seqCodeIndex]
            
                nmrRes = chain.sortedResidues()[firstSeqId - firstResidueSeqId + seqCodeIndex]
                nmrRes.seqCode = seqCode
                nmrRes.seqInsertCode = seqInsertCode
          
          #
          # Setting chain mapping in formatChainCodeDict
          #
          
          self.setSeqAndInsertCodeToChainMapping(formatChainCode,seqAndInsertCodesList,chain,firstSeqId-firstResidueSeqId,verbose = verbose)

      else:
        # Only really bomb out if nothing set, otherwise try to at least map what we can...
        if not self.setMappingFromChainMappings:
          self.status = False
        else:
          print "  Missing mappings for format chains: %s. Continuing because other mappings were set." % str(self.formatChains.keys())

      interaction.destroy()

    #
    # Setting self.formatChainCodeDict to {} is a way to invalidate this routine! ??
    #

    if self.formatChainCodeDict == {}:
      self.status = False
    
    #
    # Return status
    #

    return self.status
    
  def setSeqAndInsertCodeToChainMapping(self,formatChainCode,seqAndInsertCodesList,chain,seqIdOffset,verbose = 0):

    #
    # Set chain mapping in formatChainCodeDict and setting residues seqCodes
    #

    oldSeqInsertCode = None
    oldSeqCode = seqAndInsertCodesList[0][0]
    resPos = None
    
    chainMappingDictSetupList = []

    for formatIndex in range(0,len(seqAndInsertCodesList)):
      
      (seqCode,seqInsertCode) = seqAndInsertCodesList[formatIndex]
      
      #
      # Track changes in the seqCode so resPos remains valid...
      #
      
      seqCodeJump = seqCode - oldSeqCode

      #
      # In case the insertion code changes AND no seqCodeJump, then assume
      # that seqId is one up from before
      #
      
      if resPos == None:
        resPos = seqIdOffset
      elif seqCodeJump == 0 and oldSeqInsertCode != seqInsertCode:
        resPos += 1
      else:
        resPos += seqCodeJump
            
      #
      # Stop when running out of residues...
      #

      if resPos >= len(chain.residues):
        break

      residue = chain.sortedResidues()[resPos]
      
      #
      # Only add new mapping when seqInsertCode changes
      #
      # This means that if there is a real sequence code 'jump' linkResonances has to
      # be run twice...
      #

      if oldSeqInsertCode != seqInsertCode:
        
        if formatIndex != 0:
          chainMappingDictSetupList[-1][3] = seqAndInsertCodesList[formatIndex-1][0]
        
        firstSeqId = residue.seqId
        offset = seqCode - firstSeqId
        
        chainMappingDictSetupList.append([firstSeqId,seqInsertCode,offset,None])
        
      oldSeqCode = seqCode
      oldSeqInsertCode = seqInsertCode
    
    # Add last seqCode to last in setup list
    chainMappingDictSetupList[-1][3] = seqAndInsertCodesList[-1][0]

    for (firstSeqId,seqInsertCode,offset,formatLastSeqCode) in chainMappingDictSetupList:
      self.setChainMappingDicts(chain,formatChainCode,firstSeqId = firstSeqId, seqInsertCode = seqInsertCode, offset = offset,formatLastSeqCode = formatLastSeqCode,verbose = verbose)      
 
  def getChainMappingsSeqCodeInfo(self,chainMappings):
  
    #
    # Gets detailed mapping info for seqCode handling (slightly
    # different from getChainMappingsSeqIdInfo in linkResonances!!)
    #
    
    specificInfo = None
    if not hasattr(self,'specialSeqCodeMapping') or not self.specialSeqCodeMapping:
      specificInfo = self.getSpecificChainMappingsSeqCodeInfo(chainMappings)
    
    if specificInfo:
      (seqCodeToIdMapping,seqCodeToChainMapping,seqCodes) = specificInfo
      
    else:   
      seqCodeToIdMapping = {}
      seqCodes = []
      seqCodeToChainMapping = {}
  
      for chainMapping in chainMappings:
        firstSeqCodeKey = (chainMapping.formatFirstSeqCode,chainMapping.formatSeqInsertCode)
        seqCodeToIdMapping[firstSeqCodeKey] = chainMapping.firstSeqId
        seqCodeToChainMapping[firstSeqCodeKey] = chainMapping
        if chainMapping.formatLastSeqCode != None:
          lastSeqCodeKey = (chainMapping.formatLastSeqCode,chainMapping.formatSeqInsertCode)
        else:
          lastSeqCodeKey = (None,chainMapping.formatSeqInsertCode)
        seqCodes.append((firstSeqCodeKey,lastSeqCodeKey))    
      
      seqCodes.sort()
      seqCodes.reverse()
    
    return (seqCodeToIdMapping,seqCodeToChainMapping,seqCodes)

  def getSpecificChainMappingsSeqCodeInfo(self,chainMappings):
    
    return None

  def setChainMappingDicts(self,chain,formatChainCode,chainMapping = None, firstSeqId = 1,seqInsertCode = defaultSeqInsertCode,offset = 0,formatLastSeqCode = None, verbose = 0):
  
    if not chainMapping:

      #
      # Set chainmapping if not given
      #
      
      keywds = {'chainCode': chain.code,
                'formatChainCode': returnMemopsLine(formatChainCode),
                'firstSeqId': firstSeqId,
                'formatFirstSeqCode': firstSeqId + offset,
                'formatSeqInsertCode': seqInsertCode}
      
      if formatLastSeqCode != None:
        keywds['formatLastSeqCode'] = formatLastSeqCode
        
      chainMapping = self.nameMapping.newChainMapping(**keywds)

    #
    # Possible to have two chains linked to the same resChain (from file!)
    #
    # Also possible to have 'jumps' in sequence that are handled by two
    # different chain mappings for the same chain (06/2004)
    #

    if not self.formatChainCodeDict.has_key(formatChainCode):
      self.formatChainCodeDict[formatChainCode] = []
    self.formatChainCodeDict[formatChainCode].append(chainMapping)    

    if not self.ccpChainDict.has_key(chain):
      self.ccpChainDict[chain] = []
      
    self.ccpChainDict[chain].append(chainMapping)
    self.chainMappingChain[chainMapping] = chain

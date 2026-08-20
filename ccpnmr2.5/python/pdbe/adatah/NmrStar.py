#
# Contains constants, functions required to handle NMR-STAR files
#

from pdbe.adatah.Util import duplicateResonances
from pdbe.adatah.Constants import chemCompArchiveDataDir

#
# Constants
#

from ccpnmr.format.general.Constants import tagSep
from ccp.format.general.Constants import defaultMolCode

#
# General NMR-STAR handler class
#

class NmrStarHandler:
      
  def readNmrStarFile(self,nmrStarFile,chemCompArchiveDataDir = chemCompArchiveDataDir, formatKey = 'NmrStar', components = None, version = '3.0',**otherKeywds):

    #
    # Read project from nmrStar format
    #

    localKeywds = {}

    if self.presets.has_key('readProject') and self.presets['readProject'].has_key('keywds'):
      localKeywds = self.presets['readProject']['keywds']

    if components:
      localKeywds['components'] = components
      
    if otherKeywds:
      localKeywds.update(otherKeywds)

    self.formatObjectDict[formatKey].readProject(nmrStarFile,
                                                 entry = self.entry,
                                                 molSystem = self.molSystem,
                                                 strucGen = self.strucGen,
                                                 useOriginalChainCode = True,
                                                 linkAtoms = False,
                                                 version = version,
                                                 chemCompPath = chemCompArchiveDataDir,
                                                 minimalPrompts = 1,
                                                 **localKeywds)

    #
    # In case reading constraints and is a multimer where resonances have to be copied, do it now...
    #

    if self.presets.has_key('duplicateResonances'):

      nmrConstraintStore = self.strucGen.nmrConstraintStore
      duplicateResonances(nmrConstraintStore,'nmrStar',self.presets['duplicateResonances'])

  def exportOriginalNmrStar(self,nmrStarFile,entry = None):
    
    # This is really obsolete - don't use except for 2.1 files
    
    if not entry:
      entry = self.ccpnProject.currentNmrEntryStore.sortedEntries()[0]

    #self.formatObjectDict['NmrStar'].writeProject(nmrStarFile,entry = entry,useOriginalData = 1,minimalPrompts = 1,exportAll = 1, headerComment = originalHeaderComment) 
    self.formatObjectDict['NmrStar'].writeProject(nmrStarFile,entry = entry,useOriginalData = 1,minimalPrompts = 1,exportAll = 1) 


  #
  # TODO TODO THIS ALSO IN PdbRefDB!!! NEED TO CLEAN THIS UP AND GENERALISE!!
  # THIS CODE is updated - NOT so in PdbRefDb!
  #
  
  def setBmrbNmrStarMapping(self,bmrbNmrStarFile):
  
    print "# Matching NMR-STAR file to molSystem... "
      
    #
    # Get started with the import...
    #

    from ccp.format.nmrStar.projectIO import NmrStarProjectFile

    forceChainMappings = []

    nmrStarFile = NmrStarProjectFile(bmrbNmrStarFile)
    nmrStarFile.read(verbose = 0)
    
    tagSep = nmrStarFile.tagSep
    
    copyMapping = None # Hack for linking RDC values...

    sequenceInfo = {}
    
    #
    # Get chains for chemical shift and coupling constant info - only map those
    #
    
    for (fileType,measurementType) in (('chemShiftFiles','chemShifts'),('jCouplingFiles','jCouplingValues'),('rdcFiles','rdcValues')):

      for valuesFile in getattr(nmrStarFile,fileType):
        for value in getattr(valuesFile,measurementType):
        
          if measurementType == 'chemShifts':
            suffixList = ['']
          elif measurementType in ('jCouplingValues','rdcValues'):
            suffixList = ['1','2']
          
          for suffix in suffixList:
            molCode = getattr(value,"molCode%s" % suffix)            
            seqCode = getattr(value,"seqCode%s" % suffix)

            if hasattr(value,"resLabel%s" % suffix):
              resLabel = getattr(value,"resLabel%s" % suffix)
            else:
              resLabel = None

            if molCode.count(tagSep):
              molCode = molCode.replace(tagSep,'_')
            if not sequenceInfo.has_key(molCode):
              sequenceInfo[molCode] = {}
            if not sequenceInfo[molCode].has_key(seqCode):
              sequenceInfo[molCode][seqCode] = resLabel
            elif sequenceInfo[molCode][seqCode] != resLabel:
              print "  ERROR: mismatch in resLabel for molCode %s, seqCode %d (is %s and %s)! Check this!!" % (molCode,seqCode,resLabel,sequenceInfo[molCode][seqCode])

    resLabelList = {}

    for molCode in sequenceInfo.keys():

      seqCodeList = sequenceInfo[molCode].keys()
      seqCodeList.sort()

      resLabelList[molCode] = []
      
      if seqCodeList[0] != None and seqCodeList[-1] != None:

        for seqCode in range(seqCodeList[0],seqCodeList[-1]):
          
          resLabel = None
          if sequenceInfo[molCode].has_key(seqCode):
            resLabel = sequenceInfo[molCode][seqCode]
            if resLabel and len(resLabel) > 3:
              resLabel = resLabel[:3]
              print "  Warning: truncating BMRB code '%s' to '%s'!" % (sequenceInfo[molCode][seqCode],resLabel)

          resLabelList[molCode].append((resLabel,seqCode))

      if not resLabelList[molCode]:
        del(resLabelList[molCode])

    if not resLabelList:
      print "\n*** No residue label list created for NMR-STAR file, no mapping possible ***\n" 
      return
    
    #
    # Hack for RDCs - molecule codes very often missing, so don't get handled like they should.
    #
    
    if len(resLabelList) == 2:
      molCodes = resLabelList.keys()
      molCodes.sort()
      if molCodes.count(defaultMolCode):
        otherMolCode = molCodes[not molCodes.index(defaultMolCode)]
        sameChain = True
        for (resLabel,seqCode) in resLabelList[defaultMolCode]:
          for (otherResLabel,otherSeqCode) in resLabelList[otherMolCode]:
            if seqCode == otherSeqCode:
              if resLabel and otherResLabel and resLabel != otherResLabel:
                sameChain = False
              break

          if not sameChain:
            break
    
        if sameChain:
          copyMapping = otherMolCode
          
    #
    # Now check if a mapping was given, and set up linkResonances info
    #
  
    if self.presets.has_key('linkResonances'):
      if self.presets['linkResonances'].has_key('keywds'):
        keywds = self.presets['linkResonances']['keywds']
        if 'forceChainMappings' in keywds.keys() or 'forceDefaultChainMapping' in keywds.keys():
        
          if 'forceChainMappings' in keywds.keys():
          
            # Check and track residues that don't match between the sequences.
            mismatchedResidues = []

            for (chainCode,formatChainCode,startSeqId,seqOffset) in self.presets['linkResonances']['keywds']['forceChainMappings']:
                    
              chain = self.molSystem.findFirstChain(code = chainCode)

              if chain and resLabelList.has_key(formatChainCode):

                resSeqIdLen = len(chain.residues)
                
                # Find corresponding BMRB bit
                startSeqCode = startSeqId + seqOffset
                startResLabelIndex = 0
                resLabelListLen = len(resLabelList[formatChainCode])
                for resLabelIndex in range(resLabelListLen):
                  if resLabelList[formatChainCode][resLabelIndex][1] == startSeqCode:
                    startResLabelIndex = resLabelIndex
                    break
                
                # Now see if residue codes match
                for indexOffset in range(resSeqIdLen - startSeqId + 1):
                
                  if startResLabelIndex + indexOffset == resLabelListLen:
                    for endResidueIndex in range(indexOffset,resSeqIdLen-startSeqId + 1):
                      mismatchedResidues.append(chain.findFirstResidue(seqId = startSeqId + endResidueIndex))
                    break
                
                  residue = chain.findFirstResidue(seqId = startSeqId + indexOffset)
                  chemComp = residue.molResidue.chemComp

                  cifNamingSystem = chemComp.findFirstNamingSystem(name = 'CIF')
                  cifCode = None
                  if cifNamingSystem:
                    cifCode = cifNamingSystem.findFirstChemCompSysName().sysName
                  
                  alternateCode = residue.ccpCode.upper()
                  
                  formatResLabel = resLabelList[formatChainCode][startResLabelIndex + indexOffset][0]
                  if alternateCode != formatResLabel and (not cifCode or cifCode != formatResLabel):
                    mismatchedResidues.append(residue)
          
            if mismatchedResidues:
              self.presets['linkResonances']['keywds']['ignoreResidues'] = mismatchedResidues
              print "\n*** Ignoring %d residues that do not match between PDB and BMRB files ***\n" % len(mismatchedResidues)
              
          print "  Exiting mapping..."
          return True

    else:
      self.presets['linkResonances'] = {}

    if not self.presets['linkResonances'].has_key('keywds'):
      self.presets['linkResonances']['keywds'] = {}

      
    #
    # If alignment information was given, try to use that.
    #
    #print self.alignmentInfo

    noMatches = {}
    
    if hasattr(self,'alignmentInfo') and self.alignmentInfo:
    
      pdbChainCodes = self.alignmentInfo.keys()
      pdbChainCodes.sort()
    
      for pdbChainCode in pdbChainCodes:
        bmrbChainCode = self.alignmentInfo[pdbChainCode][0]

        if bmrbChainCode.count(tagSep):
          bmrbChainCode = bmrbChainCode.replace(tagSep,'_')

        chain = self.molSystem.findFirstChain(code=pdbChainCode)
        seqOffset = 0
        otherDiffs = []
        for (pdbIndex,pdbStatus,bmrbIndex,bmrbStatus) in self.alignmentInfo[pdbChainCode][1]:
          # BMRB insertion compared to PDB
          if not pdbStatus:
            # Residues in BMRB sequence before PDB sequence starts
            if pdbIndex == 1:
              seqOffset += 1
            else:
              otherDiffs.append((pdbIndex,1))
          # BMRB deletion compared to PDB
          elif not bmrbStatus:
            # Residues in PDB sequence before BMRB sequence starts
            if bmrbIndex == 1:
              seqOffset -= 1
            else:
              otherDiffs.append((pdbIndex,-1))
          else:
            if not noMatches.has_key((pdbChainCode,bmrbChainCode)):
              noMatches[(pdbChainCode,bmrbChainCode)] = []
            residue = chain.findFirstResidue(seqId = pdbIndex)
            resLabelInfo = (bmrbStatus,bmrbIndex)
            noMatches[(pdbChainCode,bmrbChainCode)].append((residue,resLabelInfo))
          
        forceChainMappings.append((pdbChainCode,bmrbChainCode,1,seqOffset))
        
        for (resIndex,addOffset) in otherDiffs:
          seqOffset += addOffset
          forceChainMappings.append((pdbChainCode,bmrbChainCode,resIndex,seqOffset))
              
      print "\n*** Chain mappings set by alignment information ***\n"

    #
    # If no mappings given, try to find a match
    #
    
    if not forceChainMappings:
    
      for chain in self.molSystem.sortedChains():

        resSeqIds = [residue.seqId for residue in chain.sortedResidues()]

        resSeqIdLen = len(resSeqIds)

        matchFound = False

        for seqPos in range(resSeqIdLen):

          seqId = resSeqIds[seqPos]
          residue = chain.findFirstResidue(seqId = seqId)
                
          pdbResLabels = [residue.ccpCode.upper()]
                  
          cifNamingSys = residue.chemCompVar.chemComp.findFirstNamingSystem(name = 'CIF')
          if cifNamingSys and cifNamingSys.mainChemCompSysName:
            cifCode = cifNamingSys.mainChemCompSysName.sysName
            if not cifCode in pdbResLabels:
              pdbResLabels.append(cifCode)
              
          #print seqPos, pdbResLabels

          for molCode in resLabelList.keys():

            resLabelListLen = len(resLabelList[molCode])

            shortestSequenceLen = min(resSeqIdLen,resLabelListLen)
            
            #
            # Also loop over start residues for the beginning of the BMRB sequence. Start at first match, see what comes out...
            #
            
            resLabelIndexMax = int(resLabelListLen * 0.2)
            # Fix for very short sequences
            if resLabelIndexMax < 3 and resSeqIdLen > 5:
              if resLabelListLen > 1:
                resLabelIndexMax = 2
              else:
                resLabelIndexMax = resLabelListLen
            

            for resLabelIndexStart in range(resLabelIndexMax):
              if resLabelList[molCode][resLabelIndexStart][0] in pdbResLabels:
              
                #print "START MATCH", seqPos, resLabelIndexStart, pdbResLabels, resLabelList[molCode][resLabelIndexStart][0]
                #print
              
                matches = 1
                noMatchesCurrent = []

                # indexOffset just starts to count at 1, is relevant for both CCPN side and BMRB side sequence.
                for indexOffset in range(1,resLabelListLen - resLabelIndexStart):

                  resLabelIndex = resLabelIndexStart + indexOffset
                  otherSeqPos = seqPos + indexOffset

                  if otherSeqPos >= resSeqIdLen:
                    break

                  otherResidue = chain.findFirstResidue(seqId = resSeqIds[otherSeqPos])
                  
                  #
                  # Use CIF code if available - sometimes better for mapping.
                  #
                  
                  otherResLabels = [otherResidue.ccpCode.upper()]
                  
                  cifNamingSys = otherResidue.chemCompVar.chemComp.findFirstNamingSystem(name = 'CIF')
                  if cifNamingSys:
                    for ccsn in cifNamingSys.findAllChemCompSysNames():
                      if not ccsn.sysName in otherResLabels:
                        otherResLabels.append(ccsn.sysName)

                  if otherResidue and (resLabelList[molCode][resLabelIndex][0] in otherResLabels or not resLabelList[molCode][resLabelIndex][0]):
                    matches += 1
                  else:
                    #print otherResLabels, resLabelList[molCode][resLabelIndex]
                    noMatchesCurrent.append((otherResidue,resLabelList[molCode][resLabelIndex]))

                #
                # Cutoff is 95% match... but take seqPos and resLabelIndexStart into account!
                #
                
                if seqPos <= resLabelIndexStart:
                  startPos = max(resLabelIndexStart,seqPos)
                else:
                  startPos = min(resLabelIndexStart,seqPos)
                
                #if 1:#matches > 70:
                #  print matches, shortestSequenceLen, pdbResLabels, resLabelIndexStart, seqPos
                #  print  resSeqIdLen,resLabelListLen            
                #  print noMatchesCurrent
                
                
                # Also allow for small differences if sequence is long enough, or if everything matches
                if matches >= (shortestSequenceLen - startPos) * 0.95 or (shortestSequenceLen > 6 and shortestSequenceLen - matches <= 2) or (not noMatchesCurrent and resLabelIndexStart == 0):
                  forceChainMappings.append([chain.code,molCode,seqId,resLabelList[molCode][resLabelIndexStart][1] - seqId])
                  matchFound = True
                  break

            # Only look for one matching chain
            if matchFound:
              break

          if matchFound:
            if noMatchesCurrent:
              noMatches[(chain.code,molCode)] = noMatchesCurrent
            break
   
    #
    # Set automatically determined mappings, if any
    #

    if forceChainMappings:
      
      if copyMapping:
        for chainMapping in forceChainMappings:
          if chainMapping[1] == copyMapping:
            newChainMapping = chainMapping[:1] + [defaultMolCode] + chainMapping[2:]
            forceChainMappings.append(newChainMapping)
        
    
      self.presets['linkResonances']['keywds']['forceChainMappings'] = forceChainMappings
      print "\n*** Setting chain mapping automatically to: %s ***\n" % str(forceChainMappings)

      if noMatches:
      
        residueList = []
        
        for (chainCode,molCode) in noMatches.keys():
          print "      WARNING: Mismatches in sequence mapping between CCPN chain '%s' and BMRB chain '%s':" % (chainCode,molCode)

          for (residue,resLabelInfo) in noMatches[(chainCode,molCode)]:
            if residue:
              residueInfo = "%d.%s" % (residue.seqId,residue.ccpCode)
              residueList.append(residue)
            else:
              residueInfo = "None"

            print "          - %s  <-> %d.%s" % (residueInfo,resLabelInfo[1],resLabelInfo[0].capitalize())

          print
        print
        
        # Also ignore residues that do not match.
        if residueList:
          self.presets['linkResonances']['keywds']['ignoreResidues'] = residueList
          print "\n*** Ignoring residues that do not match between PDB and BMRB files ***\n"

      returnValue = True

    else:
      print "\n*** Automatic chain mapping did not work ***\n"
      print "DEBUG INFO:"
      print resLabelList
      
      returnValue = False
      
    return returnValue

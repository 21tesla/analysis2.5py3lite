from ccp.general.Util import createOneLetterSequence

from ccp.general.Constants import standardResidueCcpCodes, code1LetterToCcpCodeDict

from ccpnmr.format.general.Util import updateResonanceNamesDict, getNameInfo

from pdbe.adatah.Util import AlignNeedlemanWunsch, getAlignmentInfo

class SequenceCompareError(StandardError):

  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

class SequenceCompare:

  def __init__(self):
  
    self.setCompareList()

  def setCompareList(self):

    """
    The compareList determines how individual amino acids or nucleotides are identified
    from the atom names in a measurement or restraint list
    (in case no three or one letter codes are available!)
    """

    HA_list = ('HA1','HA2','HA3','HA*','HA#','HA+','HA%','1HA','2HA','QA')
    HB_list = ('HB1','HB2','HB3','1HB','2HB','3HB','QB','HB*','HB#','HB%','HB+')
    HG1_list = ('HG1*','HG1#','HG1+','HG1%','CG1','1HG1','2HG1','3HG1','QG1')
    HG2_list = ('HG2*','HG2#','HG2+','HG2%','CG2','1HG2','2HG2','3HG2','QG2')
    HD1_list = ('HD1*','HD1#','HD1+','HD1%','QD1')
    HD2_list = ('HD2*','HD2#','HD2+','HD2%','QD2')

    # The first set of tuples are atom names that HAVE to appear, the second set of tuples is atom names that SHOULD NOT appear
    self.compareList = {

      ('protein',): [

          ('g',(HA_list,)                ,(HB_list,)),
          ('l',(HD1_list, HD2_list)      ,(HG1_list,HG2_list)),
          ('v',(HG1_list, HG2_list)      ,(HD1_list,HD2_list)),
          ('i',(HD1_list, HG2_list)      ,None),
          ('q',(('HE21','HE22'),HB_list) ,None),
          ('w',(('HZ2','HZ3','HH2'),)    ,None),
          ('t',(HG2_list,('HB',))        ,None),
          ('n',(HD2_list,HB_list)        ,(HD1_list,HG1_list,HG2_list))

                 ],

      ('DNA','RNA'): [

          ('t',(('1H5M','2H5M','3H5M','H51','H52','H53','H71','H72','H73'),), None),
          ('c',(('1H4','2H4','H41','H42'),)                                 , None),
          ('a',(('1H6','2H6','H61','H62','H2'),)                            , None),
          ('g',(('1H2','2H2','H21','H22','H1'),)                            , None)

                      ]

    }


  def getFormatFileInformation(self,resonances,importFormatName):
 
    """
    This function will pull the original chain/residue/atom names from the resonances
    that were created on file import (and input to this function).
    
    Created are:
    
    self.formatFileChainDict = {'chainCode': [(seqCode1,seqInsertCode1),(seqCode2,seqInsertCode2), ...]}
    self.formatFileResidueDict = {'chainCode': {(seqCode,seqInsertCode): [residueLabel, [atomName1, atomName2, ...]]}}
    
    self.formatFileHasResLabels  True/False, indicates whether residue labels are available.
    
    """
 
    self.formatFileChainDict = {}
    self.formatFileResidueDict = {}

    self.formatFileHasResLabels = False

    for resonance in resonances:

      #
      # Now set up resonance stuff... for either normal or fixed resonances
      #

      resonanceLinked = False
      resonanceSet = resonance.resonanceSet
      
      # This should never be true, keep anyway.
      if resonanceSet:
        atomSets = resonanceSet.atomSets
        for atomSet in atomSets:
          if atomSet.atoms:
            resonanceLinked = True
            break

      #
      # Only handle if not linked
      #

      if not resonanceLinked:

        resNames = updateResonanceNamesDict(resonance,{},importFormatName)

        resLabel = None
        resLabelAppData = resonance.findFirstApplicationData(application = importFormatName, keyword = 'origResLabel')
        if resLabelAppData:
          resLabel = resLabelAppData.value

        for resName in resNames:

          # 
          # Check if the name can be decomposed or not (no use if not!)
          #

          (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resName,verbose = 0)

          if seqCode != None:

            # Keep track of reported chainCodes for linking later on
            if not self.formatFileChainDict.has_key(chainCode):
              self.formatFileChainDict[chainCode] = []
              self.formatFileResidueDict[chainCode] = {}

            seqCodeKey = (seqCode,seqInsertCode)

            # Also keep track of reported seqCodes and seqInsertCodes
            if seqCodeKey not in self.formatFileChainDict[chainCode]:
              self.formatFileChainDict[chainCode].append(seqCodeKey)
              self.formatFileResidueDict[chainCode][seqCodeKey] = [resLabel,[]] 
              if resLabel:
                self.formatFileHasResLabels = True

            if atomName not in self.formatFileResidueDict[chainCode][seqCodeKey][1]:
              self.formatFileResidueDict[chainCode][seqCodeKey][1].append(atomName)

            if not self.formatFileResidueDict[chainCode][seqCodeKey][0] and resLabel:
              self.formatFileResidueDict[chainCode][seqCodeKey][0] = resLabel
              self.formatFileHasResLabels = True


  def createCcpnChainInformation(self,chains):
 
    """
    Takes a list of CCPN chains as input.
    
    Creates a list with an element per chain:
    
    (molType,chainCode,seqTuple,seqString)
    
    molType     The CCPN molecular type
    chainCode   The CCPN chain.code
    seqTuple    A tuple with three letter codes for ligands and non-linear polymers.
                Is uppercased for comparison purposes
    seqString   The sequence in one-letter codes for linear polymers
    
    
    TODO: not adapted for branched carbohydrates.
    """
  
    self.ccpnChainInfo = []
    
    for chain in chains:
    
      molType = chain.molecule.molType
      seqString = None    
      seqTuple = None

      seqList = []
      resTypeList = []
      for residue in chain.sortedResidues():
        seqList.append((residue.seqId,residue.ccpCode))
        resTypeList.append(residue.molType)
        
      if len(seqList) == 1:
        seqTuple = (seqList[0][-1].upper(),)
      elif molType == 'other':
        # Just make sure it's not 'other' because of some unknown residues...
        if resTypeList.count('other') > len(seqList) * 0.8:
         seqTuple = tuple([seqEl[-1].upper() for seqEl in seqList])

      # Only do this for polymers... 
      if not seqTuple:

        if not molType:
          for tmpMolType in ('protein','DNA','RNA'):
            if resTypeList.count(tmpMolType) > len(seqList) * 0.8:
              molType = tmpMolType
              break
        
        seqString = createOneLetterSequence(molType,seqList)
        
      self.ccpnChainInfo.append((molType,chain.code,seqTuple,seqString))

  def createFormatFileChainInformation(self):

    formatChainCodes = self.formatFileResidueDict.keys()
    formatChainCodes.sort()
    
    if not formatChainCodes:
      raise SequenceCompareError('No format chain codes in imported NMR data file, cannot continue!')

    self.formatFileChainInfo = []
    self.compareToUncertain = True
    
    for formatChainCode in formatChainCodes:

      formatSeqList = []

      seqCodeKeys = self.formatFileChainDict[formatChainCode]
      seqCodeKeys.sort()

      oldSeqCode = None

      for seqCodeKey in seqCodeKeys:

        if oldSeqCode == None:
          oldSeqCode = seqCodeKey[0]
          # Try to catch sequences that have other residues before the first reported one...
          #for i in range(1, seqCodeKey[0]):
          #  formatSeqList.append((i,'',[]))

        if seqCodeKey[0] > oldSeqCode + 1:
          for i in range(oldSeqCode + 1, seqCodeKey[0]):
            formatSeqList.append((i,'Xxx',[]))

        formatSeqList.append((seqCodeKey[0],self.formatFileResidueDict[formatChainCode][seqCodeKey][0],self.formatFileResidueDict[formatChainCode][seqCodeKey][1]))

        oldSeqCode = seqCodeKey[0]
        
      #
      # Create the sequence tuple and one-letter string for comparison to the CCPN chains
      #
      formatSeqTuple = tuple()
      formatSeqString = ""
      
      if self.formatFileHasResLabels:
        
        molType = None
        isOneLetterCodes = False
        
        for formatSeqInfo in formatSeqList:
          
          if formatSeqInfo[1]:
            if formatSeqInfo[1].capitalize() in standardResidueCcpCodes['protein']:
              molType = 'protein'
              break
            elif formatSeqInfo[1] == 'U':
              molType = 'RNA'
            elif formatSeqInfo[1].capitalize() in code1LetterToCcpCodeDict['protein']:
              molType = 'protein'
              isOneLetterCodes = True
              break
        
        if not molType and formatSeqList[0][1].capitalize in standardResidueCcpCodes['DNA']:
          molType = 'DNA'
        
        if molType:
          if not isOneLetterCodes:
            formatSeqString = createOneLetterSequence(molType,formatSeqList)
          else:
            formatSeqString = ''.join([formatSeqInfo[1].upper() for formatSeqInfo in formatSeqList])
          
          formatSeqTuple = tuple([str(seqEl[1]).upper() for seqEl in formatSeqList])
          self.compareToUncertain = False

      else:
        
        #
        # Also need to determine molecular type for multiple incoming chains, linear polymer at least
        #
        
        if len(formatSeqList) == 1:
          formatSeqString = 'X'
        
        else:

          for molTypes in self.compareList.keys():
          
            tmpFormatSeqString = ""
            residueCompareList = self.compareList[molTypes]

            for (seqCode,resLabel,atomList) in formatSeqList:

              codesList = []

              for (code,compareAtomLists,refuseAtomLists) in residueCompareList:

                matches = 0

                for compareAtomNameList in compareAtomLists:
                  for compareAtomName in compareAtomNameList:
                    if compareAtomName in atomList:
                      matches += 1
                      break

                if refuseAtomLists:
                  for refuseAtomNameList in refuseAtomLists:
                    for refuseAtomName in refuseAtomNameList:
                      if refuseAtomName in atomList:
                        matches = -1
                        break

                if matches == len(compareAtomLists):
                  codesList.append(code)

              if len(codesList) == 1:
                tmpFormatSeqString += codesList[0]
              else:
                tmpFormatSeqString += 'X'
                
            #
            # Pick the one with least 'X'
            #
            
            if not formatSeqString or formatSeqString.count('X') > tmpFormatSeqString.count('X'):
              formatSeqString = tmpFormatSeqString

      self.formatFileChainInfo.append((formatChainCode,formatSeqTuple,formatSeqString,seqCodeKeys[0][0]))

  def compareFormatFileToCcpnInfo(self):
    
    self.forceChainMappings = []
    self.seqAlignmentsUsed = {}

    #
    # First do direct comparison, or if ligand connect directly
    #
    
    # Make local copies in case full info needed elsewhere
    ccpnChainInfo = self.ccpnChainInfo[:]
    formatFileChainInfo = self.formatFileChainInfo[:]
    
    ccpnSingleResidue = []
    formatSingleResidue = []
    
    for ccpnChainItem in self.ccpnChainInfo:
      for formatFileChainItem in self.formatFileChainInfo:
        
        (molType,chainCode,seqTuple,seqString) = ccpnChainItem
        (formatChainCode,formatSeqTuple,formatSeqString,formatFirstSeqCode) = formatFileChainItem
        
        if seqTuple == formatSeqTuple:
          
          seqOffset = formatFirstSeqCode - 1
          self.forceChainMappings.append((chainCode,formatChainCode,1,seqOffset))
          
          ccpnChainInfo.pop(ccpnChainInfo.index(ccpnChainItem))          
          formatFileChainInfo.pop(formatFileChainInfo.index(formatFileChainItem))
          
        else:        
          if seqTuple and len(seqTuple) == 1:
            ccpnSingleResidue.append(ccpnChainItem)
          elif len(formatSeqTuple) == 1:
            formatSingleResidue.append(formatFileChainItem)
          
    #
    # If only one single residue, connect up here
    # 
    
    if len(ccpnSingleResidue) == 1 and len(formatSingleResidue) == 1:
        
      ccpnChainItem = ccpnSingleResidue[0]
      formatFileChainItem = formatSingleResidue[0]
        
      (molType,chainCode,seqTuple,seqString) = ccpnChainItem
      (formatChainCode,formatSeqTuple,formatSeqString,formatFirstSeqCode) = formatFileChainItem

      seqOffset = formatFirstSeqCode - 1
      self.forceChainMappings.append((chainCode,formatChainCode,1,seqOffset))

      ccpnChainInfo.pop(ccpnChainInfo.index(ccpnChainItem))          
      formatFileChainInfo.pop(formatFileChainInfo.index(formatFileChainItem))
    
    #
    # Full sequence comparison for linear polymers for anything that's left
    #
    
    for ccpnChainItem in ccpnChainInfo:
    
      bestAlignInfo = None
      
      (molType,chainCode,seqTuple,seqString) = ccpnChainItem
    
      for formatFileChainItem in formatFileChainInfo:
        
        (formatChainCode,formatSeqTuple,formatSeqString,formatFirstSeqCode) = formatFileChainItem
                
        if formatSeqString and seqString:

          sequences = [formatSeqString,seqString]

          alignInfo = self.getAlignInfo(sequences)

          seqAlignInfo = alignInfo[-1]
            
          if not bestAlignInfo:
            bestAlignInfo = (alignInfo,seqAlignInfo,formatFileChainItem)
          elif bestAlignInfo[0][0] < alignInfo[0]:
            bestAlignInfo = (alignInfo,seqAlignInfo,formatFileChainItem)
        
      #
      # Use best match, but only if positive
      #

      if bestAlignInfo and bestAlignInfo[0][0] > 0:

        (alignInfo,seqAlignInfo,formatFileChainItem) = bestAlignInfo

        (molType,chainCode,seqTuple,seqString) = ccpnChainItem
        (formatChainCode,formatSeqTuple,formatSeqString,formatFirstSeqCode) = formatFileChainItem

        seqDifferences = getAlignmentInfo(seqAlignInfo)

        # TODO: some consistency checking? No holes allowed!!!
        #print seqDifferences

        otherDiffs = []
        #noMatches = []

        seqOffset = formatFirstSeqCode - 1
        formatIndexStart = formatFirstSeqCode  # Make sure to start numbering correctly for format info!

        for (ccpnIndex,ccpnStatus,formatIndex,formatStatus) in seqDifferences:
          # CCPN sequence insertion compared to format
          if not ccpnStatus:
            # Residues in format sequence before CCPN sequence starts
            if ccpnIndex == 1:
              seqOffset += 1
            #else:
            #  otherDiffs.append((ccpnIndex,1))
          # Format sequence deletion compared to CCPN
          elif not formatStatus:
            # Residues in CCPN sequence before format sequence starts
            if formatIndex == 1:
              seqOffset -= 1
            #else:
            #  otherDiffs.append((ccpnIndex,-1))
          #else:
          #  resLabelInfo = (formatStatus,formatIndex)
          #  residue = chain.findFirstResidue(seqId = ccpnIndex)
          #  noMatches.append((residue,resLabelInfo))

        #
        # Set this match
        #

        self.forceChainMappings.append((chainCode,formatChainCode,1,seqOffset))

        for (resIndex,addOffset) in otherDiffs:
          seqOffset += addOffset
          self.forceChainMappings.append((chainCode,formatChainCode,resIndex,seqOffset))
        
        #
        # Track sequence alignments that were used in case necessary for other apps
        #
        
        self.seqAlignmentsUsed[chainCode] = seqAlignInfo
        
        # Remove from comparison list - but ONLY if similar amount of residues!
        seqFraction = len(seqString) * 1.0 / len(formatSeqString)
        if 0.9 < seqFraction < 1.1:
          formatFileChainInfo.pop(formatFileChainInfo.index(formatFileChainItem))
            
    return self.forceChainMappings

  def getAlignInfo(self,sequences):
  
    seqSwapped = False
    if len(sequences[0]) < len(sequences[1]):
      sequences.reverse()
      seqSwapped = True

    align = AlignNeedlemanWunsch(sequences[0],sequences[1],verbosity=1,compareToUncertain=self.compareToUncertain)
    #align.make_graph()
    #align.searchPaths()
    #align.printScores()
    alignInfo = align.getBestMatchInfo(verbosity=1)

    #(score, total, degeneracy, matches, mismatches, gaps, seqAlignInfo) = alignInfo

    seqAlignInfo = alignInfo[-1]
    if seqSwapped:
      # Make sure CCPN based sequence is first, always
      seqAlignInfo = (seqAlignInfo[1],seqAlignInfo[0])
      alignInfo = alignInfo[:-1] + (seqAlignInfo,)
    
    return alignInfo

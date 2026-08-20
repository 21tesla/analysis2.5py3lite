"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for PDB sequence (coordinate) file

Copyright (C) 2005-2007 Wim Vranken (European Bioinformatics Institute)

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

# Import general functions
from memops.universal.Util import returnInt

# Import header reader
from ccp.format.pdb.coordinatesIO import PdbCoordinateFile
from ccp.format.pdb.generalIO import PdbGenericFile
from ccp.format.general.formatIO import Sequence, SequenceElement

from ccp.format.general.Constants import bioPolymerCodes, chainCodeString
from ccp.format.general.Constants import defaultSeqInsertCode, bondTypes

from ccp.format.general.Util import getSeqAndInsertCode

from ccp.general.ChemCompOverview import chemCompStdDict

#####################
# Class definitions #
#####################  

def addToCurChain(curChain,pdbCoordinate):

  curChain[0].append(pdbCoordinate.resName)
  curChain[1].append(pdbCoordinate.seqCode)
  curChain[2].append(pdbCoordinate.insertionCode)
  curChain[3].append([])
  curChain[4].append(pdbCoordinate.hetFlag)
  
  #print "ADDED", curChain[0][-1], curChain[1][-1], curChain[4][-1]
  return curChain[3][-1]
  
def insertInCurChain(curChain,pdbCoordinate):

  returnCoordList = None

  for seqIndex in range(len(curChain[0])):
    seqCode = curChain[1][seqIndex]
    insertionCode = curChain[2][seqIndex]
    clist = curChain[3][seqIndex]

    # Insert at first position where smaller - should be sequentially numbered!
    if seqCode > pdbCoordinate.seqCode or (seqCode == pdbCoordinate.seqCode and pdbCoordinate.insertionCode != defaultSeqInsertCode and insertionCode > pdbCoordinate.insertionCode):
      curChain[0].insert(seqIndex,pdbCoordinate.resName)
      curChain[1].insert(seqIndex,pdbCoordinate.seqCode)
      curChain[2].insert(seqIndex,pdbCoordinate.insertionCode)
      curChain[3].insert(seqIndex,[])
      curChain[4].insert(seqIndex,pdbCoordinate.hetFlag)
      returnCoordList = curChain[3][seqIndex]
      #print "INSERTED", seqIndex, curChain[0][seqIndex], curChain[1][seqIndex], curChain[4][seqIndex]
      break
      
    # Ignore if exists!
    elif seqCode == pdbCoordinate.seqCode and insertionCode == pdbCoordinate.insertionCode:
      returnCoordList = clist
      break
  
  if returnCoordList == None:
    returnCoordList = addToCurChain(curChain,pdbCoordinate)
    
  return returnCoordList
      
def getChainId(chainId,chainIdsUsed,chainCodeString):
  
  for tChainId in chainCodeString:
    if tChainId not in chainIdsUsed.keys():
      chainId = tChainId
      break

  return chainId

class PdbSequenceFile(PdbGenericFile):
  """
  Information on file level
  """
  def initialize(self):

    self.sequences = []
    
    self.bioPolymerCodes = bioPolymerCodes
    
  def read(self, pdbFile = None, ignoreResNames = None, version = '3.0', verbose = 0):

    if not pdbFile:
    
      #
      # First read coordinate file (if not passed in)
      #

      pdbFile = PdbCoordinateFile(self.name)
      pdbFile.read(maxNum = 1,ignoreResNames=ignoreResNames)
      
    self.version = pdbFile.version
    self.pdbFile = pdbFile
    
    if not self.version:
      self.version = version
    
    #
    # Get out missing residue info - important for matching to SEQRES info
    #
    
    fakeModelCoordinates = {}
    
    if self.pdbFile.remarks.has_key('465') and self.version[0] == '3':
      
      try:
        remarkText = self.pdbFile.remarks['465']
        remarkLines = remarkText.split("\n")

        dataLine = False
        for remarkLine in remarkLines:
          if remarkLine.count("RES C") or remarkLine.count("SSEQI "):
            dataLine = True
          elif dataLine:
            columns = remarkLine.split()
            if columns:
              if len(columns) == 4:
                (modelKeyString,resName,chainId,seqAndInsertCode) = columns
              else:
                (resName,chainId,seqAndInsertCode) = columns
                modelKeyString = '1'
                
              (seqCode,insertionCode) = getSeqAndInsertCode(seqAndInsertCode)
              
              modelKey = int(modelKeyString)
              
              if not fakeModelCoordinates.has_key(modelKey):
                fakeModelCoordinates[modelKey] = {}
              if not fakeModelCoordinates[modelKey].has_key(chainId):
                fakeModelCoordinates[modelKey][chainId] = []

              fakeModelCoordinates[modelKey][chainId].append(FakeModelCoordinate(chainId,resName,seqCode,insertionCode))

      except:
        print "  Warning: trying to use remark 465 for missing residue info, but failed."
            
    # TODO: Missing coordinate info is in remark '470'!!

    #
    # Set some reference info
    #
    
    if self.version in ('2.5',) or self.version[0] == '3':
      rnaAtoms = ("HO2'",)
      dnaAtoms = ("H2''",)
    
    else:
      rnaAtoms = ('2HO','2HO*','HO2*')
      dnaAtoms = ('2H2','2H2*','H2**')
    
    #
    # Reorganize coordinate information...
    #
    
    chainIdInfo = {}    
    oldChainId = None
    
    modelCoordKeys = pdbFile.modelCoordinates.keys()
    modelCoordKeys.sort()
    modelCoordRefKey = modelCoordKeys[0]
    
    coordToSeq = {}
    bondedCoords = {}
    
    independentHetGroups = []
    
    modelCoordinates = pdbFile.modelCoordinates[modelCoordRefKey]
    
    if fakeModelCoordinates.has_key(modelCoordRefKey):
      # Have to insert the model coordinates at the right position, more or less...
      for pdbCoordinate in modelCoordinates:
        if pdbCoordinate.chainId in fakeModelCoordinates[modelCoordRefKey].keys():
          pdbCoordinateIndex = modelCoordinates.index(pdbCoordinate)
          
          modelCoordinates = modelCoordinates[:pdbCoordinateIndex] + fakeModelCoordinates[modelCoordRefKey][pdbCoordinate.chainId] + modelCoordinates[pdbCoordinateIndex:]
          
          del(fakeModelCoordinates[modelCoordRefKey][pdbCoordinate.chainId])

    for pdbCoordinate in modelCoordinates:
    
      if pdbCoordinate.chainId != oldChainId:

        oldChainId = pdbCoordinate.chainId
        oldSeqCode = pdbCoordinate.seqCode
        oldInsertionCode = pdbCoordinate.insertionCode

        #print "NEW", oldSeqCode, pdbCoordinate.seqCode, self.version

        # Now also have to make provision for HETATM stuff defined at end of file... this is version 3.0 and higher only
        # Can be sure one chain code defines one molecule, so just use that one if already available...
        if self.version[0] == '3' and chainIdInfo.has_key(pdbCoordinate.chainId):
          if not self.isIndependentHetGroup(pdbCoordinate):
            curChain = chainIdInfo[pdbCoordinate.chainId][-1]
            curCoordList = insertInCurChain(curChain,pdbCoordinate)
          else:
            (curChain,curCoordList) = self.createNewChain(pdbCoordinate,chainIdInfo,oldChainId,1)
            independentHetGroups.append(curChain)
        
        else:
          (curChain,curCoordList) = self.createNewChain(pdbCoordinate,chainIdInfo,oldChainId,2)
          
          if self.isIndependentHetGroup(pdbCoordinate):
            independentHetGroups.append(curChain)
                
      elif oldSeqCode == pdbCoordinate.seqCode and oldInsertionCode == pdbCoordinate.insertionCode:
      
        pass
        
      elif (oldSeqCode == pdbCoordinate.seqCode + 1 or oldSeqCode == pdbCoordinate.seqCode - 1) or \
           (oldSeqCode == pdbCoordinate.seqCode and oldInsertionCode != pdbCoordinate.insertionCode):
        #print "SEQUENTIAL", oldSeqCode, pdbCoordinate.seqCode, self.version
        if not self.isIndependentHetGroup(pdbCoordinate):
          # For safety better to insert - will add if necessary.
          curCoordList = insertInCurChain(curChain,pdbCoordinate)
        else:
          (curChain,curCoordList) = self.createNewChain(pdbCoordinate,chainIdInfo,oldChainId,3)
          independentHetGroups.append(curChain)

        oldSeqCode = pdbCoordinate.seqCode
        oldInsertionCode = pdbCoordinate.insertionCode
            
      else:
        #print "JUMP",oldSeqCode, pdbCoordinate.seqCode, self.version
        # Sequence jump, same chain code
        
        oldSeqCode = pdbCoordinate.seqCode
        oldInsertionCode = pdbCoordinate.insertionCode
        
        # Here also make sure that only use existing chain is not an independent HET group, 
        if self.version[0] == '3':
          if not self.isIndependentHetGroup(pdbCoordinate) and chainIdInfo[pdbCoordinate.chainId][-1] not in independentHetGroups:
            curChain = chainIdInfo[pdbCoordinate.chainId][-1]
            curCoordList = insertInCurChain(curChain,pdbCoordinate)
          else:
            (curChain,curCoordList) = self.createNewChain(pdbCoordinate,chainIdInfo,oldChainId,4)
            independentHetGroups.append(curChain)
          
        else:
          (curChain,curCoordList) = self.createNewChain(pdbCoordinate,chainIdInfo,oldChainId,5)

      if not isinstance(pdbCoordinate,FakeModelCoordinate):
        curCoordList.append(pdbCoordinate)
      
      #
      # Also check if there's any linkings to set molTypeBonds later on...
      #
      
      if pdbCoordinate.bonds:
        for bondType in bondTypes:
          if pdbCoordinate.bonds.has_key(bondType):
            if not bondedCoords.has_key(bondType):
              bondedCoords[bondType] = []
              
            for otherPdbCoordinate in pdbCoordinate.bonds[bondType]:
              coordKey = [pdbCoordinate,otherPdbCoordinate]
              coordKey.sort()
              coordKey = tuple(coordKey)
              
              if not coordKey in bondedCoords[bondType]:
                bondedCoords[bondType].append(coordKey)

    chainIds = chainIdInfo.keys()
    
    """ 
    for chainId in chainIds:
      print "CHAIN '%s'" % chainId
      for curChain in chainIdInfo[chainId]:
        print curChain[0]
        print curChain[1]
        print curChain[2]
    
    """ 
    
    #
    # Create molecules...
    #
    
    chainIdsUsed = {}
    
    #print pdbFile.chains

    for pdbChain in pdbFile.chains:
    
      chainId = pdbChain.chainId
      
      if not chainIdInfo.has_key(chainId):
        print "  Error: unrecognized %s code '%s' for polymer... ignored." % (self.format,chainId)
        print chainIdInfo.keys()
        continue

      self.sequences.append(PdbSequence(molName = pdbChain.molName))
      currentSequence = self.sequences[-1]

      currentSequence.setCommonNames(pdbChain.synonyms)
      currentSequence.setDetails(pdbChain.details)
            
      pdbSeq = []

      checkNucleicType = 0
      resNames = pdbChain.resNames[:]
      
      # This is a hack to find out if the sequence numbering is negative - i.e. 1bfw. Is this valid PDB though?
      resNamesString = '_'.join(resNames)      

      lastSeqCode = None
      lastInsertionCode = None
      
      #
      # Try to find a match...
      #
      # TODO: this will fail if residues in SEQRES are missing in coordinates!!!
      #                      with old style PDB files where residues are missing (but are they there then?)
      
      while (1):
     
        matchFound = 0
      
        for matchChain in chainIdInfo[chainId]:
        
          #
          # Check if sequence reversed... only for reasonable amount of residues.
          #
                  
          if len(matchChain[0]) > 5 and len(resNames) > 5:
            matchChainList = []
            for mce in matchChain[0]:
              if mce != None:
                matchChainList.append(mce)
                
            matchChainString = '_'.join(matchChainList)
            
            if not matchChainString.count(resNamesString):
              resNames.reverse()
              resNamesString = '_'.join(resNames)
              
              if not matchChainString.count(resNamesString):
                resNames.reverse()
              else:
                print "WARNING: reversing residue names. Sequence reversed in PDB file!"
          
          #
          # Also need to check if there's het residues in the matchchain...
          #
          
          endHetResidues = 0

          for matchIndex in range(len(matchChain[0])-1,0,-1):
            matchHetFlag = matchChain[4][matchIndex]
            
            if matchHetFlag:
              if self.version[0] != '3' or (resNames and matchChain[0][matchIndex] not in resNames):
                endHetResidues += 1
            else:
              break
          #print endHetResidues
          
          beginHetResidues = 0
          
          # Also track het groups at beginning - but only if they do not appear in sequence any more!
          # In version 3.0 and higher, ACE type groups are not listed any more in SEQRES...
          if self.version[0] == '3':
            for matchIndex in range(len(matchChain[0])):
              if matchChain[4][matchIndex] and matchChain[0][matchIndex] not in resNames:
                beginHetResidues += 1
              else:
                break        
          
          #
          # Now start matching...
          #
         
          #print matchChain[0]
          #print resNames
      
          for matchIndex in range(beginHetResidues,len(matchChain[0])):
          
            matchResName = matchChain[0][matchIndex]

            if resNames and matchResName == resNames[0]:
              
              #
              # Try to catch bad one-residue matches...
              #
              
              if len(resNames) == 1 and len(resNames) != len(matchChain[0]):
                continue
                
              matchStart = matchIndex
              matchNum = 1
              
              # Always take highest number, then check if indexes still OK...
              resNamesLen = len(resNames)
              matchChainLen = len(matchChain[0])
              
              matchAddLen = matchChainLen-matchIndex
              
              #print resNamesLen, matchChainLen

              matchAddIndexCorrection = 0
              for matchAddIndex in range(1,matchAddLen):
              
                #print "START",matchAddIndex, matchAddIndexCorrection
                #print 
                if (matchAddIndex + matchAddIndexCorrection >= resNamesLen) or (matchIndex + matchAddIndex >= matchChainLen):
                  #print "BREAK", matchAddIndex + matchAddIndexCorrection, resNamesLen, matchIndex + matchAddIndex, matchChainLen
                  break

                if matchChain[0][matchIndex + matchAddIndex] == resNames[matchAddIndex + matchAddIndexCorrection]:
                  matchNum +=1
                # Special case for HET groups - often not in SEQRES
                elif matchChain[4][matchIndex + matchAddIndex] and matchChainLen > matchIndex + matchAddIndex + 1:
                  #print "HETS",matchIndex + matchAddIndex, matchChain[4][matchIndex + matchAddIndex]
                  #print "HETS",matchChain[0][matchIndex + matchAddIndex + 1], resNames[matchAddIndex + matchAddIndexCorrection]
                  if matchChain[0][matchIndex + matchAddIndex + 1] == resNames[matchAddIndex + matchAddIndexCorrection]:
                    matchAddIndexCorrection -= 1
                    matchNum += 1
                    #print matchAddIndexCorrection
                  else:
                    break
                else:
                  break
                  
              #print matchNum, len(matchChain[0]), len(resNames), endHetResidues  
              #print matchChain[0]
              #print resNames
              #sys.exit()

              if matchNum == len(matchChain[0]) or matchNum >= len(resNames) or matchNum == (len(matchChain[0]) - endHetResidues):

                matchFound = 1
                
                matchCountStart = matchIndex
                matchCountEnd = matchIndex + matchNum

                if self.version[0] == '3':
                  matchCountStart -= beginHetResidues
                  matchCountEnd   += endHetResidues

                for matchCount in range(matchCountStart,matchCountEnd):
                  #print beginHetResidues, endHetResidues, matchIndex, matchNum, matchCountStart, matchCountEnd, matchCount, len(matchChain[0])
                  resName = matchChain[0][matchCount]
                  seqCode = matchChain[1][matchCount]

                  insertionCode = matchChain[2][matchCount]
                  pdbCoords = matchChain[3][matchCount]
                  
                  #
                  # Renumber previous residues if necessary...
                  #

                  if matchCount == matchIndex and pdbSeq and pdbSeq[-1][1] == None:
                    for j in range(-1,-len(pdbSeq)-1,-1):
                      if pdbSeq[j][1] == None:
                        pdbSeq[j][1] = seqCode + j
                        pdbSeq[j][2] = insertionCode
                        
                  #
                  # Set DNA/RNA if possible
                  #

                  if resName.strip() in self.bioPolymerCodes['DNA'][0] + self.bioPolymerCodes['RNA'][0]:
                    checkNucleicType = 1
                  else:
                    checkNucleicType = 0

                  pdbSeq.append([resName,seqCode,insertionCode,'',[],[]])
                  
                  #
                  # Here get rid of resNames. For remediated files, this has to be a bit more sophisticated
                  # because HET groups do not appear in the sequence when at start/end (e.g. NH2).
                  #
                  #print len(resNames), matchIndex, matchNum, endHetResidues, beginHetResidues, matchCount
                  # TODO Wim: removed -endHetResidues from equation - possibly wrong, but seems to work for version 3 files..
                  if resNames and (self.version[0] != '3' or ((matchIndex + matchNum ) >= matchCount >= (matchIndex-beginHetResidues))):
                    resNames.pop(0)

                  #elif resNames:
                  #  print beginHetResidues, endHetResidues
                  #  print matchIndex + matchNum, matchCount, matchIndex, resNames[0], 'not popped'
                  
                  #
                  # Keep track of missing residues
                  #
                  
                  lastSeqCode = seqCode
                  lastInsertionCode = insertionCode

                  for pdbCoord in pdbCoords:

                    pdbSeq[-1][4].append(pdbCoord.atomName)
                    pdbSeq[-1][5].append(pdbCoord)

                    if checkNucleicType:

                      if pdbCoord.atomName in rnaAtoms:
                        pdbSeq[-1][3] = 'RNA'
                        #print "  Set %d to RNA" % seqCode
                      elif pdbCoord.atomName in dnaAtoms and not pdbSeq[-1][3]:
                        pdbSeq[-1][3] = 'DNA'
                        #print "  Set %d to DNA" % seqCode
                        
                  
                  # Do not reuse info once it's set! (new code!)
                  matchChain[0][matchCount] = matchChain[1][matchCount] = matchChain[2][matchCount] = matchChain[3][matchCount] = None

                break
            
        if resNames:
          #print resNames
          if not matchFound:
            if lastSeqCode != None:
              # This will fail if sequence numbered in reverse...
              lastSeqCode += 1
            #print resNames[0], lastSeqCode
            pdbSeq.append([resNames[0],lastSeqCode,defaultSeqInsertCode,'',[],[]])
            resNames.pop(0)
            
        if not resNames:
          break
      
      #
      # Error if nothing matches...
      #
      
      if not pdbSeq:
        print "  Error: no sequence match found in %s coordinates for chain %s!" % (self.format,chainId)
        continue
      
      #print pdbChain
      #for pdbEl in pdbSeq:
      #  print pdbEl
      #print
      
      #
      # Keep track of codes that are used so there is no duplication...
      #
      
      resetChainId = False
      
      if not chainIdsUsed.has_key(chainId):
        chainIdsUsed[chainId] = []
      
      seqIdList = []
      for pdbSeqEl in pdbSeq:
        #print pdbSeqEl
        seqId = "%d%s" % (pdbSeqEl[1],pdbSeqEl[2])
        if not seqId in chainIdsUsed[chainId]:
          seqIdList.append(seqId)
        else:
          seqIdList = []
          break
          
      if seqIdList:
        chainIdsUsed[chainId].extend(seqIdList)
        
      else:
        chainId = getChainId(pdbChain.chainId,chainIdsUsed,chainCodeString)
        resetChainId = True
      
      currentSequence.setChainCode(chainId)
      currentSequence.setOriginalChainCode(pdbChain.chainId)

      #
      # Set some additional info...
      #
      
      proteinCircularCoords = []
      
      for i in range(0,len(pdbSeq)):

        code3Letter = pdbSeq[i][0]

        seqCode = pdbSeq[i][1]
        insertionCode = pdbSeq[i][2]
        residueType = pdbSeq[i][3]
        atomNames = pdbSeq[i][4]
        pdbCoords = pdbSeq[i][5]
          
        #
        # Check if this element is actually a hetero group...
        #
        
        addedHetGroup = False
                  
        for pdbHetGroup in pdbFile.hetGroups[:]:

          hetGroupId = None
          tempAtomNames = []
          tempPdbCoords = []

          if pdbHetGroup.chainId == pdbChain.chainId and pdbHetGroup.seqNum == seqCode:

            for pdbCoordinate in pdbFile.modelCoordinates[modelCoordRefKey]:

              if pdbCoordinate.hetFlag != 1:
                continue

              tResName = pdbCoordinate.resName
              tSeqCode = pdbCoordinate.seqCode
              tInsertionCode = pdbCoordinate.insertionCode

              if tSeqCode == seqCode and tResName == code3Letter and tInsertionCode == insertionCode:

                if not hetGroupId:
                  hetGroupId = (tResName,tSeqCode,tInsertionCode)
                  tempAtomNames.append(pdbCoordinate.atomName)
                  tempPdbCoords.append(pdbCoordinate)
                elif hetGroupId == (tResName,tSeqCode,tInsertionCode):
                  tempAtomNames.append(pdbCoordinate.atomName)
                  tempPdbCoords.append(pdbCoordinate)
                else:
                  hetGroupId = None
                  tempAtomNames = None
                  break

          if hetGroupId:

            atomNames = tempAtomNames
            pdbCoords = tempPdbCoords
            pdbFile.hetGroups.pop(pdbFile.hetGroups.index(pdbHetGroup))
            addedHetGroup = True
            #print hetGroupId, pdbHetGroup.seqNum, pdbHetGroup.Id
            break

        if addedHetGroup:       
          print "  Warning: adding HET residue %s in position %d (index %d) of chain %s!" % (code3Letter,seqCode,len(currentSequence.elements),currentSequence.chainCode)

        #
        # Add this element
        #

        currentSequence.elements.append(PdbSequenceElement(seqCode,code3Letter,residueType,insertionCode,atomNames))
        
        for pdbCoord in pdbCoords:
          coordToSeq[pdbCoord] = currentSequence.elements[-1]
        
        #
        # Track for circular proteins...
        #
        
        if residueType == '':
          if i == 0:
            for pdbCoord in pdbCoords:
              if pdbCoord.atomName == 'N':
                proteinCircularCoords.append(pdbCoord)
                break
          
          elif i == len(pdbChain.resNames) - 1:
            for pdbCoord in pdbCoords:
              if pdbCoord.atomName == 'C':
                proteinCircularCoords.append(pdbCoord)
                break
      
      #
      # Check if circular protein...
      #

      if len(proteinCircularCoords) == 2 and proteinCircularCoords[0].bonds.has_key('covalent'):
        for bondedCoord in proteinCircularCoords[0].bonds['covalent']:
          if bondedCoord == proteinCircularCoords[1]:
            currentSequence.setCircular()
            break

    
    #
    # Hetgroup handling - also try to group by chain (?).
    # Question is whether this would always be relevant... so checking if bound to each other
    #
    # TODO TODO: ALSO HAVE TO CREATE MOLECULAR BONDS IF THERE ARE 'link' bonds!!
    #
    
    pdbHetChains = {}
    pdbHetChainIds = []
    
    for pdbHetGroup in pdbFile.hetGroups:

      chainId = pdbHetGroup.chainId
      seqCode = pdbHetGroup.seqNum
      iCode = pdbHetGroup.iCode
      pdbHetChainIndex = -1
      
      #print chainId,seqCode, iCode, pdbHetGroup.Id
      if not pdbHetChains.has_key(chainId):
        pdbHetChains[chainId] = []
        pdbHetChainIds.append(chainId)        
      
      #
      # Sort according to seqNum, if possible
      #

      for i in range(len(pdbHetChains[chainId])):
        if pdbHetGroup.seqNum < pdbHetChains[chainId][i][0].seqNum:
          pdbHetChains[chainId].insert(i,[pdbHetGroup,[]])
          pdbHetChainIndex = i
          break

      if pdbHetChainIndex == -1:
        pdbHetChains[chainId].append([pdbHetGroup,[]])

      #
      # Add coordinates
      #      

      for modelCoordKey in modelCoordKeys:
        for pdbCoordinate in pdbFile.modelCoordinates[modelCoordKey]:
          if pdbCoordinate.chainId == chainId:

            if pdbCoordinate.hetFlag != 1:
              continue

            if seqCode == pdbCoordinate.seqCode and iCode == pdbCoordinate.insertionCode:
              pdbHetChains[chainId][pdbHetChainIndex][1].append(pdbCoordinate)      
              #print pdbCoordinate

        if pdbHetChains[chainId][pdbHetChainIndex][1]:
          break

    pdbHetChainIds.sort()
    
    doNotLinkHetGroupIds = ('HOH',)

    # Only look for matches within one chain...
    for pdbHetChainId in pdbHetChainIds:
    
      #pdbHetGroupsLeft = pdbHetChains[pdbHetChainId][:] # Need this??!
      pdbHetGroups = []
      
      #print pdbHetChainId
      for (pdbHetGroup,pdbCoords) in pdbHetChains[pdbHetChainId]:
        #print "DOING",pdbHetGroup.Id,pdbHetGroup.seqNum,pdbHetGroup.iCode
                
        pdbLinkedHetGroups = []
        matchFound = 0
        for pdbCoord in pdbCoords:
          for bondType in ('link','covalent'):
            if pdbCoord.bonds.has_key(bondType):
              #print "  LINK %s" % bondType
              # Do NOT match HOH - gets too complicated to create 'chains'.
              # Will still create bonds...
              if pdbHetGroup.Id in doNotLinkHetGroupIds:
                #print "skip"
                continue
              for otherPdbCoord in pdbCoord.bonds[bondType]:
              
                for (otherPdbHetGroup,otherPdbCoords) in pdbHetChains[pdbHetChainId]:
                  
                  # Intrachemcomp bonding possible as well
                  if otherPdbHetGroup == pdbHetGroup:
                    continue
                
                  # Do NOT match HOH - gets too complicated to create 'chains'.
                  # Will still create bonds...
                  if otherPdbHetGroup.Id in doNotLinkHetGroupIds:
                    #print "skipother"
                    continue
                    
                  if otherPdbCoord in otherPdbCoords:
                    #print "  MATCH", otherPdbHetGroup.Id,otherPdbHetGroup.seqNum

                    matchFound = 1
                    if not [otherPdbHetGroup,otherPdbCoords] in pdbLinkedHetGroups:
                      pdbLinkedHetGroups.append([otherPdbHetGroup,otherPdbCoords])
                      #print "     added"
                    break
                    
              if matchFound:
                break
        
        if (not matchFound or not pdbHetGroups) and (not pdbHetGroup.Id in doNotLinkHetGroupIds):
          #print "new added to hetgroups (1)", pdbHetGroup.Id, pdbHetGroup.seqNum
          pdbHetGroups.append([[pdbHetGroup,pdbCoords]])
          if [pdbHetGroup,pdbCoords] in pdbLinkedHetGroups:
            pdbLinkedHetGroups.pop(pdbLinkedHetGroups.index([pdbHetGroup,pdbCoords]))
          pdbHetGroups[-1].extend(pdbLinkedHetGroups)
          #pdbHetGroupsLeft.pop(pdbHetGroupsLeft.index((pdbHetGroup,pdbCoords)))
        
        else:
          hetGroupExists = False
          
          for pdbHetGroupSet in pdbHetGroups:
            for pdbLinkedHetGroupInfo in pdbLinkedHetGroups:
              if pdbLinkedHetGroupInfo in pdbHetGroupSet and [pdbHetGroup,pdbCoords] not in pdbHetGroupSet:
                #print "  ADD", pdbHetGroup.Id
                pdbHetGroupSet.append([pdbHetGroup,pdbCoords])
                
            if [pdbHetGroup,pdbCoords] in pdbHetGroupSet:
              hetGroupExists = True
          
          if not hetGroupExists:
            #print "new added to hetgroups (2)", pdbHetGroup.Id, pdbHetGroup.seqNum
            pdbHetGroups.append([[pdbHetGroup,pdbCoords]])
            if [pdbHetGroup,pdbCoords] in pdbLinkedHetGroups:
              pdbLinkedHetGroups.pop(pdbLinkedHetGroups.index([pdbHetGroup,pdbCoords]))
            pdbHetGroups[-1].extend(pdbLinkedHetGroups)

      #
      # Remap some data if required
      #
      
      mapChainIds = {}
      
      for pdbHetGroupSet in pdbHetGroups:
      
        pdbHetGroupDict = {}
        pdbHetGroupList = []
        
        for (pdbHetGroup,pdbCoords) in pdbHetGroupSet:
          pdbHetGroupId = (pdbHetGroup.seqNum,pdbHetGroup.iCode)
          pdbHetGroupList.append(pdbHetGroupId)
          pdbHetGroupDict[pdbHetGroupId] = (pdbHetGroup,pdbCoords)
          
        pdbHetGroupList.sort()
        
        #print len(pdbHetGroupSet)
        #for (pdbHetGroup,pdbCoords) in pdbHetGroupSet:
        #  print pdbHetGroup.Id,pdbHetGroup.seqNum, pdbHetGroup.iCode, "   ",
        #print
        
        for pdbHetGroupId in pdbHetGroupList:
        
          (pdbHetGroup,pdbCoords) = pdbHetGroupDict[pdbHetGroupId]

          if pdbHetGroupId == pdbHetGroupList[0]:
            self.sequences.append(PdbSequence(molName = pdbHetGroup.chemicalName))
            currentSequence = self.sequences[-1]

            currentSequence.setCommonNames(pdbHetGroup.synonyms)
            currentSequence.setDetails(pdbHetGroup.description)

            #
            # Keep track of codes that are used so there is no duplication...
            #
            
            chainId = pdbHetGroup.chainId
            resetChainId = False
            
            if mapChainIds.has_key(chainId):
              chainId = mapChainIds[chainId]

            if not chainIdsUsed.has_key(chainId):
              chainIdsUsed[chainId] = []

            seqIdList = []
            for tempPdbHetGroupId in pdbHetGroupList:
              (tempPdbHetGroup,tempPdbCoords) = pdbHetGroupDict[tempPdbHetGroupId]
              
              seqId = None
              
              for pdbCoordinate in pdbCoords:
                resName = pdbCoordinate.resName
                
                if resName == tempPdbHetGroup.Id:

                  seqCode = pdbCoordinate.seqCode
                  insertionCode = pdbCoordinate.insertionCode
            
                  seqId = "%d%s" % (seqCode,insertionCode)
                  break
                  
              if seqId != None and not seqId in chainIdsUsed[chainId]:
                seqIdList.append(seqId)
              else:
                seqIdList = []
                break
                
            originalChainId = chainId

            if seqIdList:
              chainIdsUsed[chainId].extend(seqIdList)

            else:
              oldChainId = chainId
              chainId = getChainId(chainId,chainIdsUsed,chainCodeString)
              mapChainIds[oldChainId] = chainId
              chainIdsUsed[chainId] = []
              resetChainId = True

            currentSequence.setChainCode(chainId)
            currentSequence.setOriginalChainCode(originalChainId)

          if resetChainId:

            pdbHetGroup.refChainId = chainId
            seqCode = pdbHetGroup.seqNum
            iCode = pdbHetGroup.iCode
            for modelKey in pdbFile.modelCoordinates.keys():
              for pdbCoordinate in pdbFile.modelCoordinates[modelKey]:
                if pdbCoordinate.hetFlag == 1 and pdbCoordinate.chainId == pdbHetGroup.chainId:
                  if seqCode == pdbCoordinate.seqCode and iCode == pdbCoordinate.insertionCode:
                    pdbCoordinate.refChainId = chainId

          #self.seqNum = returnInt(seqNum) # is reset from seqCode in coordinates...
          #self.numHetAtoms = returnInt(numHetAtoms)

          #
          # Assuming that order in which residues appear is same as in SEQRES in PDB file header,
          # but start point can be different.
          #

          atomNames = []
          
          seqCode = None
          insertionCode = None

          for pdbCoordinate in pdbCoords:

            resName = pdbCoordinate.resName

            if resName == pdbHetGroup.Id:

              seqCode = pdbCoordinate.seqCode
              insertionCode = pdbCoordinate.insertionCode
              atomNames.append(pdbCoordinate.atomName)

          #
          # Add this element
          #
          
          if seqCode != None:
          
            residueType = 'other'
            
            # CCPN carbohydrate check...
            if pdbCoords:
              for tCcpCode in chemCompStdDict['carbohydrate'].keys():
                # TODO is it worth it here to make a dictionary of the cifCodes?
                if pdbCoords[0].resName == chemCompStdDict['carbohydrate'][tCcpCode][1]:
                  residueType = 'carbohydrate'

            currentSequence.elements.append(PdbSequenceElement(seqCode,pdbHetGroup.Id,residueType,insertionCode,atomNames))
        
            for pdbCoord in pdbCoords:
              coordToSeq[pdbCoord] = currentSequence.elements[-1]
            
          elif not ignoreResNames or not pdbHetGroup.Id in ignoreResNames:
          
            print "  Error: could not find sequence elements for hetGroup '%s'! Not included in sequence..." % pdbHetGroup.Id

    #
    # Set the bonds on the sequence element level!
    #

    for bondType in bondTypes:
      if bondedCoords.has_key(bondType):
        for pdbCoords in bondedCoords[bondType]:
          seqEls = []
          for pdbCoord in pdbCoords:
            if coordToSeq.has_key(pdbCoord):
              seqEl = coordToSeq[pdbCoord]
            else:
              print "  Error: could not find sequence element for coordinate %d ('%s.%d.%s') to set bond..." % (pdbCoord.serial,pdbCoord.chainId,pdbCoord.seqCode,pdbCoord.atomName)
              seqEl = None
            seqEls.append(seqEl)
          
          #
          # Set the bond, but ONLY if between different sequence elements.
          # Internal chemComp bonds should be covered by E-MSD info!
          #
          
          if not None in seqEls and seqEls[0] != seqEls[1]:
            for i in range(0,2):
              atomName = pdbCoords[i].atomName
              seqEl = seqEls[i]
              otherAtomName = pdbCoords[not i].atomName
              otherSeqEl = seqEls[not i]
              
              seqEl.setBond(bondType,atomName,otherSeqEl,otherAtomName)

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

    #
    # Clean up...
    #
    
    for i in range(len(self.sequences)-1,0,-1):
      sequence = self.sequences[i]
      if not sequence.elements:
        self.sequences.pop(i)
        print "Removed empty sequence '%s', '%s'" % (sequence.molName, sequence.chainCode)


  def createNewChain(self,pdbCoordinate,chainIdInfo,oldChainId,location):
  
    if not chainIdInfo.has_key(oldChainId):
      chainIdInfo[oldChainId] = []

    #print "New chain (%d), hetFlag %d" % (location,pdbCoordinate.hetFlag)
    
    # Doublecheck that doesn't exist already if single chemComp! Just in case... is sometimes a problem
    moleculeExists = False
    for tempChain in chainIdInfo[oldChainId]:
      if len(tempChain[1]) == 1:
        if tempChain[1][0] == pdbCoordinate.seqCode and tempChain[2][0] == pdbCoordinate.insertionCode:
          moleculeExists = True
          curChain = tempChain
          curCoordList = tempChain[3][0]
          #print "EXISTING FOUND"
          break
    
    if not moleculeExists:
      chainIdInfo[oldChainId].append([[],[],[],[],[]])
      curChain = chainIdInfo[oldChainId][-1]

      curCoordList = addToCurChain(curChain,pdbCoordinate)
    
    return (curChain,curCoordList)
    
  def isIndependentHetGroup(self,pdbCoordinate):
  
    isIndependent = False
    
    # Means a TER record was read from the previous chain, this is no good...
    if pdbCoordinate.chainTerminated:
      isIndependent = True
      
    elif pdbCoordinate.hetFlag:
      pdbFile = pdbCoordinate.parent
      
      partOfChain = False
      for pdbChain in pdbFile.chains:
        if pdbCoordinate.resName in pdbChain.resNames:
          partOfChain = True
          break

      if not partOfChain:
        isIndependent = True
        #print pdbCoordinate.chainId, pdbCoordinate.seqCode

        for link in pdbFile.links:
          for (atomName,resName,chainId,seqCode,iCode) in link:
            #print "  ",chainId, seqCode
            if pdbCoordinate.chainId == chainId and pdbCoordinate.seqCode == seqCode:
              isIndependent = False
              #print 'FOUND'
              break
          
    return isIndependent
    
class PdbSequence(Sequence):

  def setCircular(self):
    self.isCircular = True

  def setCommonNames(self,value):
    self.commonNames = value

  def setDetails(self,value):
    self.details = value

  def setChainCode(self,value):
    self.chainCode = value
    
  def setOriginalChainCode(self,value):
    self.originalChainCode = value
    
  def setSecondaryStructure(self,secStrucType,secStrucSerial,specificInfo,seqEls):
    
    if not self.secStrucInfo.has_key(secStrucType):
      self.secStrucInfo[secStrucType] = {}
    
    self.secStrucInfo[secStrucType][secStrucSerial] = (specificInfo,tuple(seqEls))
 
class PdbSequenceElement(SequenceElement):

  def setFormatSpecific(self,*args):
    
    self.residueType = args[1]
    self.insertionCode = args[2]
    self.atomNames = args[3]
    
class FakeModelCoordinate:

  def __init__(self,chainId,resName,seqCode,insertionCode):
  
    self.chainId = chainId
    self.resName = resName
    self.seqCode = seqCode
    self.insertionCode = insertionCode
    
    # This info never available for het groups?
    self.hetFlag = 0
    self.bonds = None
    self.atomName = None
    
    self.chainTerminated = False

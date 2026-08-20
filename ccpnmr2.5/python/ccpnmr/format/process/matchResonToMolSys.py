"""
======================COPYRIGHT/LICENSE START==========================

matchResonToMolSys.py: Tries to match resonance info to molecular system info

Copyright (C) 2008 Wim Vranken (European Bioinformatics Institute)

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

from ccp.api.nmr import Nmr, NmrConstraint
from ccp.api.molecule import MolSystem

from ccpnmr.format.general.Util import getNameInfo, updateResonanceNamesDict
from memops.universal.Util import drawBox

import os, sys

comparisonList = (
    
  (('GLY',),('HA1','HA2','HA3','HA*','HA#','HA+','HA%','1HA','2HA','QA'),10.0,2.0),
  (('VAL','ILE','THR'),('HG2*','HG2#','HG2+','HG2%','CG2','1HG2','2HG2','3HG2','QG2'),7.0,1.5),
  (('GLN',),('HE21','HE22'),20.0,2.0),
  (('LEU','ILE'),('HD1*','HD1#','HD1+','HD1%','QD1'),5.0,1.0),
  (('LEU','ASN'),('HD2*','HD2#','HD2+','HD2%','QD2'),7.0,1.0),
  (('TRP',),('HZ2','HZ3','HH2'),6.0,1.0),
  (('THR',),('OG1',),20.0,2.0),
  (('ASP',),('OD2',),20.0,2.0),
  (('ASP','ASN'),('OD1',),10.0,1.0),
  (('GLU',),('OE2',),20.0,2.0),
  (('GLU','GLN'),('OE1',),10.0,1.0),
  (('T',),('1H5M','2H5M','3H5M','H51','H52','H53','H71','H72','H73'),5.0,1.0),
  (('C',),('1H4','2H4','H41','H42'),8.0,1.0),
  (('A',),('1H6','2H6','H61','H62','H2'),6.0,1.0),
  (('G',),('1H2','2H2','H21','H22','H1'),6.0,1.0)
                  )   

def doMatchResidues(bestMatch,comparisonList,comparisonListIndex,residues,resResidues,resResidueDictForChain,fixMatch,test):
  
  (residueCodes,atomNames,matchScore,scoreMultiplier) = comparisonList[comparisonListIndex]
    
  (bestFit,correctFit) = matchResidues(residueCodes,residues,atomNames,resResidues,resResidueDictForChain,matchScore = matchScore,test=test)
  
  if fixMatch:
    correctFit = (correctFit[0],correctFit[1],correctFit[2] * scoreMultiplier)
  
  if bestFit[2] > bestMatch[2]:
    bestMatch = bestFit
  
  if correctFit[0] == None or correctFit[1] == None:
    comparisonListIndex += 1
    
    if len(comparisonList) > comparisonListIndex:
      (correctFit,bestMatch) = doMatchResidues(bestMatch,comparisonList,comparisonListIndex,residues,resResidues,resResidueDictForChain,fixMatch,test)
  
  return (correctFit,bestMatch)

def matchResonToMolSys(resonances,molSystem,assignFormat = None,test = False, matchMode = False):

  #
  # Some quick validity checks...
  #

  if not resonances or (type(resonances) not in (type(('',)),type([]))) or \
      not (isinstance(resonances[0],Nmr.Resonance) or isinstance(resonances[0],NmrConstraint.FixedResonance)):
   print "  Invalid resonance list input - aborting."
   return None

  #
  # Get info from the resonances
  #

  resChainDict = {}
  resResidueDict = {}

  for resonance in resonances:

    #
    # Now set up resonance stuff... for either normal or fixed resonances
    #

    resonanceLinked = 0
    resonanceSet = resonance.resonanceSet

    if resonanceSet:
      atomSets = resonanceSet.atomSets
      for atomSet in atomSets:
        atoms = atomSet.atoms
        if atoms:
          resonanceLinked = 1
          break

    #
    # Only handle if not linked - useless otherwise (???)
    #

    if not resonanceLinked:

      resNames = updateResonanceNamesDict(resonance,{},assignFormat)
      
      resLabel = None
      resLabelAppData = resonance.findFirstApplicationData(application = assignFormat, keyword = 'origResLabel')
      if resLabelAppData:
        resLabel = resLabelAppData.value

      for resName in resNames:

        # 
        # Check if the name can be decomposed or not (no use if not!)
        #

        (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resName,verbose = 0)

        if seqCode != None:

          # Keep track of reported chainCodes for linking later on
          if not resChainDict.has_key(chainCode):
            resChainDict[chainCode] = []
            resResidueDict[chainCode] = {}

          seqCodeKey = (seqCode,seqInsertCode)

          # Also keep track of reported seqCodes and seqInsertCodes
          if seqCodeKey not in resChainDict[chainCode]:
            resChainDict[chainCode].append(seqCodeKey)
            resResidueDict[chainCode][seqCodeKey] = [resLabel,[]] 

          if atomName not in resResidueDict[chainCode][seqCodeKey][1]:
            resResidueDict[chainCode][seqCodeKey][1].append(atomName)
            
          if not resResidueDict[chainCode][seqCodeKey][0] and resLabel:
            resResidueDict[chainCode][seqCodeKey][0] = resLabel


  if test:
    print resResidueDict
    
  return matchToMolSys(resChainDict,resResidueDict,molSystem,test = test, matchMode = matchMode)

def matchCoordAtomsToMolSys(coords,molSystem,test = False, matchMode = False):
  
  #
  # Get info from the coords from the coordinate file...
  #

  resChainDict = {}
  resResidueDict = {}

  for coord in coords:
    chainCode = coord.refChainId
    seqCode = coord.seqCode
    seqInsertCode = coord.insertionCode
    atomName = coord.atomName
    
    resLabel = None
    if hasattr(coord,'resName'):
      resLabel = coord.resName

    if seqCode != None:

      # Keep track of reported chainCodes for linking later on
      if not resChainDict.has_key(chainCode):
        resChainDict[chainCode] = []
        resResidueDict[chainCode] = {}

      seqCodeKey = (seqCode,seqInsertCode)

      # Also keep track of reported seqCodes and seqInsertCodes
      if seqCodeKey not in resChainDict[chainCode]:
        resChainDict[chainCode].append(seqCodeKey)
        resResidueDict[chainCode][seqCodeKey] = [resLabel,[]]

      if atomName not in resResidueDict[chainCode][seqCodeKey][1]:
        resResidueDict[chainCode][seqCodeKey][1].append(atomName)

      if not resResidueDict[chainCode][seqCodeKey][0] and resLabel:
        resResidueDict[chainCode][seqCodeKey][0] = resLabel
  
  return matchToMolSys(resChainDict,resResidueDict,molSystem,test = test, matchMode = matchMode)
  
def matchToMolSys(resChainDict,resResidueDict,molSystem,test = 0, matchMode = 0):

  if not molSystem or not isinstance(molSystem,MolSystem.MolSystem):
   print "  Invalid molSystem input - aborting."
   return None

  #
  # Get info from the molSystem (number chains, ...)
  #

  chainDict = {}
  residueDict = {}

  for chain in molSystem.chains:
   chainDict[chain.code] = chain
   residueDict[chain.code] = {}
   for residue in chain.residues:
     residueDict[chain.code][residue.seqId] = residue

  #
  # Do some sorting...
  #

  for resChainCode in resChainDict.keys():
    resChainDict[resChainCode].sort()

    if len(resChainDict[resChainCode]) > 1:
      # Make sequential...
      seqCodeIndex = 1
      oldSeqCodeKey = resChainDict[resChainCode][0]
      seqCodeKey = None
    
      while seqCodeKey != resChainDict[resChainCode][-1]:
        seqCodeKey = resChainDict[resChainCode][seqCodeIndex]
        if seqCodeKey[0] == (oldSeqCodeKey[0] + 1) or (seqCodeKey[0] == oldSeqCodeKey[0] and seqCodeKey[1] != oldSeqCodeKey[1]):
          pass
        else:
          for i in range(1,seqCodeKey[0]-oldSeqCodeKey[0]):
            resChainDict[resChainCode].insert(seqCodeIndex,(oldSeqCodeKey[0] + i,oldSeqCodeKey[1]))
            seqCodeIndex += 1

        oldSeqCodeKey = seqCodeKey
        seqCodeIndex += 1
    
    for seqCodeKey in resChainDict[resChainCode]:
      if resResidueDict[resChainCode].has_key(seqCodeKey):
        resResidueDict[resChainCode][seqCodeKey][1].sort()
      else:
        resResidueDict[resChainCode][seqCodeKey] = [None,[]]
    
    if test:    
      print resChainCode
      print resChainDict[resChainCode]
      #for seqCodeKey in resChainDict[resChainCode]:
      #  print seqCodeKey, resResidueDict[resChainCode][seqCodeKey]

  #
  # Now try to match...
  #
  
  chainMappingDict = {}
  residueMappingDict = {}
  
  #
  # Chains first...
  #
  
  if len(chainDict) == len(resChainDict):
    if test:
      print "SAME AMOUNT OF CHAINS on molSystem and Resonance level!"
    unmappedCodes = []
    resChainCodes = list(resChainDict.keys())
    resChainCodes.sort()
    
    for chainCode in chainDict.keys():
      if resChainDict.has_key(chainCode):
        chainMappingDict[chainCode] = [chainCode]
        resChainCodes.pop(resChainCodes.index(chainCode))
      else:
        unmappedCodes.append(chainCode)
    
    if len(chainMappingDict) != len(chainDict):
      if len(chainDict) == 1:
        chainMappingDict[chainDict.keys()[0]] = [resChainCodes[0]]
      elif len(unmappedCodes) == 1:
        chainMappingDict[unmappedCodes[0]] = [resChainCodes[0]]
      elif not complexChainMapping(chainDict,resChainDict,resResidueDict,unmappedCodes,resChainCodes,chainMappingDict,test=test,matchMode = matchMode):
        print "  No complex mapping for equal number of chains - aborting."
        printMissingInfo(resResidueDict)
        return None
    
  else:
  
    if test:
      print "DIFFERENT AMOUNT OF CHAINS on molSystem and Resonance level!"

    unmappedCodes = chainDict.keys()
    unmappedCodes.sort()
    resChainCodes = resChainDict.keys()
    resChainCodes.sort()
  
    if len(chainDict) == 1:
      chainMappingDict[chainDict.keys()[0]] = list(resChainDict.keys())
    elif len(resChainDict) == 1:
      for chainCode in chainDict.keys():
        chainMappingDict[chainCode] = [resChainDict.keys()[0]]
    elif not complexChainMapping(chainDict,resChainDict,resResidueDict,unmappedCodes,resChainCodes,chainMappingDict,test=test,matchMode = matchMode):
      print "  No complex mapping for mismatching number of chains - aborting."
      printMissingInfo(resResidueDict)
      return None
    
  #
  # Straight residue map...
  #
  #[' ',' ',1,0]],   # Chain mapping: [chainCode, formatChainCode, firstSeqId, offset  
  #
  
  chainCodes = chainMappingDict.keys()
  chainCodes.sort()

  for chainCode in chainCodes:
    
    chain = chainDict[chainCode]
    residues = chain.sortedResidues()
    
    if test:
      print "CHAIN INFO"
      print chainCode
      for residue in residues:
        print "  %s.%d" % (residue.ccpCode,residue.seqId)
      print
    
    resChainCodes = chainMappingDict[chainCode]
    resChainsMapped = []
    
    for resChainCode in resChainCodes:

      resResidues = resChainDict[resChainCode]
      
      if not resResidues:
        continue
      
      if test:
        print "TRYING TO MATCH TO %s" % resChainCode
        for resResidue in resResidues:
          print resResidue,
        print
        print

      for seqAttr in ('seqId','seqCode'):

        if getattr(residues[0],seqAttr) == resResidues[0][0] and residues[0].seqInsertCode == resResidues[0][1] and \
           getattr(residues[-1],seqAttr) == resResidues[-1][0] and residues[-1].seqInsertCode == resResidues[-1][1]:
          
          offset = resResidues[0][0] - 1
          residueMappingDict[(chainCode,resChainCode)] = (1,offset)
          resChainsMapped.append(resChainCode)

          if test:
            print "  HARD set offset to %d" % offset          

          break
             
    if len(resChainsMapped) != len(resChainCodes):
      matches = []
      for resChainCode in resChainCodes:
        if resChainCode not in resChainsMapped:
        
          #
          # Try to match on specific characteristics...
          #
          
          #
          #  TODO: track most reasonable mapping, and track for 'match' mode!!!
          #
          
          resResidues = resChainDict[resChainCode]
          
          if test:
            print "'%s','%s'" % (chainCode,resChainCode)
          
          (match,bestMatch) = doMatchResidues((None,None,-999),comparisonList,0,residues,resResidues,resResidueDict[resChainCode],0,test)
          (resId,offset,score) = match
                            
          #
          # Try matching on special atom name...
          #
          
          if matchMode:
            if bestMatch[0] != None and bestMatch[1] != None:
              matchText = "start %d, offset %d, score %d" % (bestMatch[0],bestMatch[1],bestMatch[2])
            else:
              matchText = 'None'
              
            print "  Best match for CCPN code '%s' to code '%s': %s" % (chainCode,resChainCode,matchText)

          if len(residues) == 1:
            residue = residues[0]
            atomNames = [atom.name for atom in residue.sortedAtoms()]
            
            for resResidue in resResidues[:]:
              
              if resResidueDict[resChainCode][resResidue][1]:
                fullMatch = True

                for atomName in resResidueDict[resChainCode][resResidue][1]: 

                  if atomName not in atomNames:
                    fullMatch = False
                    break
              
              else:
                fullMatch = False

              if fullMatch:
                if test:
                  print " Matched on atomname, 1 residue: %s-%s" % (chainCode,resChainCode)

                resId = 1
                offset = resResidue[0] - resId
                resResidues.pop(resResidues.index(resResidue))
                break 

          if resId != None and offset != None:
            residueMappingDict[(chainCode,resChainCode)] = (resId,offset)
            resChainDict[resChainCode] = resResidues
            resChainsMapped.append(resChainCode)
            if (resId,offset) not in matches:
              matches.append((resId,offset))
      
      #
      # Set the info based on match (if possible)
      #
      
      if len(matches) == 1:
        (resId,offset) = matches[0]
        for resChainCode in resChainCodes:
          if resChainCode not in resChainsMapped:
            resResidues = resChainDict[resChainCode]
            
            if resResidues:
              newResId = resResidues[0][0] - offset
              residueMappingDict[(chainCode,resChainCode)] = (newResId,offset)
  
  finalChainMappings = []
  resMapDictKeys = residueMappingDict.keys()
  resMapDictKeys.sort()
  for (chainCode,resChainCode) in resMapDictKeys:
    (seqId,offset) = residueMappingDict[(chainCode,resChainCode)]
    finalChainMappings.append([chainCode,resChainCode,seqId,offset])
    if chainDict.has_key(chainCode):
      del chainDict[chainCode]
    if resChainDict.has_key(resChainCode):
      del resChainDict[resChainCode]
  
  if chainDict:
    chainCodes = chainDict.keys()
    chainCodes.sort()
    print "  Warning: Did not map CCPN chain codes %s" % (chainCodes)
  if resChainDict:
    resChainCodes = resChainDict.keys()
    resChainCodes.sort()
    print "  Warning: Did not map resonance chain codes %s" % (resChainCodes)
    printMissingInfo(resResidueDict)
  
  print "  Final mapping: %s" % str(finalChainMappings)
  
  return finalChainMappings

def printMissingInfo(resResidueDict):
  
  print drawBox("Suggested mappings, based on atom names for original restraint data",indent = "  ")
  
  resChainCodes = resResidueDict.keys()
  resChainCodes.sort()
  
  for resChainCode in resChainCodes:
    
    print "    Original chain code '%s'. Higher scores below are better." % resChainCode
    print
    
    resSeqCodes = resResidueDict[resChainCode].keys()
    resSeqCodes.sort()
    
    for resSeqCode in resSeqCodes:
      
      scoreTexts = []

      if resResidueDict[resChainCode][resSeqCode][0]:
        
        scoreTexts.append("%s (from file)" % resResidueDict[resChainCode][resSeqCode][0])
      
      else:
      
        atomNames = resResidueDict[resChainCode][resSeqCode][1]
        
        scoreList = {}
        
        for atomName in atomNames:
          for (ccpCodes,atomNameMatch,score,multiplier) in comparisonList:
            if atomName in atomNameMatch:
              for ccpCode in ccpCodes:
                if not scoreList.has_key(ccpCode):
                  scoreList[ccpCode] = score
                else:
                  scoreList[ccpCode] *= multiplier
        
        
        values = scoreList.values()
        values.sort()
        values.reverse()
        
        ccpCodes = scoreList.keys()
        ccpCodes.sort()
              
        for value in values:
          for ccpCode in ccpCodes:
            if scoreList[ccpCode] == value:
              scoreTexts.append("%s (%d)" % (ccpCode,value))
              ccpCodes.pop(ccpCodes.index(ccpCode))
              break
            
      if scoreTexts:
        print "      %-4s%s:" % (str(resSeqCode[0]),resSeqCode[1]),
        print ', '.join(scoreTexts)
    
    print
  
  print
  
  return None

def setDict(curDict,curKey,curValue):

  if not curDict.has_key(curKey):
    curDict[curKey] = []
  curDict[curKey].append(curValue)

def complexChainMapping(chainDict,resChainDict,resResidueDict,unmappedCodes,resChainCodes,chainMappingDict,test =0, matchMode = 0):
  
  returnStatus = True
  
  unmappedCodes.sort()

  #
  # Check which type...
  #

  molecules = []
  for chainCode in chainDict.keys():
    molecule = chainDict[chainCode].molecule
    if molecule not in molecules:
      molecules.append(molecule)

  if test:
    print "Molecules: ",
    for mol in molecules:
      print " %s (%d)," % (mol.name,len(mol.molResidues)),
    print
    
  #
  # Homomultimer
  #

  if len(molecules) == 1:
    if len(unmappedCodes) == len(resChainCodes):
      for i in range(0,len(unmappedCodes)):
        setDict(chainMappingDict,unmappedCodes[i],resChainCodes[0])
        resChainCodes.pop(0)

    else:
      # Try to group the resChainCodes...
      if test:
        print "RESCHAINGROUPING"
      resChainGroups = []
      for resChainCode in resChainCodes:
        resChainResidues = resChainDict[resChainCode]
        rcrFound = 0
        for i in range(0,len(resChainGroups)):
          (trcr,rccs) = resChainGroups[i]
          fullResMatch = 0
          if len(resChainResidues) == len(trcr):
            fullResMatch = 1
            for j in range(len(resChainResidues)):
              if resChainResidues[j] != trcr[j]:
                fullResMatch = 0
                break
          if fullResMatch:
            rcrFound = 1
            resChainGroups[i][1].append(resChainCode)
        if not rcrFound:
          resChainGroups.append((resChainResidues,[resChainCode]))
          
      for (rcr,rccs) in resChainGroups:
        # Assuming this means it's refering to the real codes...
        if len(unmappedCodes) == len(rccs):
          rccs.sort()
          for i in range(0,len(unmappedCodes)):
            setDict(chainMappingDict,unmappedCodes[i],rccs[i])
            resChainCodes.pop(resChainCodes.index(rccs[i]))

      if test:
        print resChainGroups
        print chainMappingDict
  
  #
  # Do other routines if not matched...
  #
          
  if resChainCodes:
  
    seqDiffs = {}
    doneResChainCodes = []

    #
    # First catch ZN, FE, ...
    #

    for resChainCode in resChainCodes:

      resChainLen = len(resChainDict[resChainCode])
      
      if resChainLen == 1:

        if test:
          print "SINGLEMETAL",resChainCode, resChainLen
          
        resResidue = resChainDict[resChainCode][0]

        for chainCode in unmappedCodes:
    
          chainLen = len(chainDict[chainCode].residues)
          
          if test:
            print "SINGLEMETAL2",chainCode, chainLen
      
          if chainLen == 1:
            residue = chainDict[chainCode].findFirstResidue()
            atomNames = [atom.name for atom in residue.sortedAtoms()]
            
            fullMatch = True

            for atomName in resResidueDict[resChainCode][resResidue][1]: 
            
              if atomName not in atomNames:
                fullMatch = False
                break
            
            if fullMatch:
              if test:
                print " Matched on atomname, 1 residue: %s-%s" % (chainCode,resChainCode)
              doneResChainCodes.append(resChainCode)
              unmappedCodes.pop(unmappedCodes.index(chainCode))
              setDict(chainMappingDict,chainCode,resChainCode)
              break
        
    #
    # Try to do rest...
    #

    if test:
      print "SEQDIFFING"

    for chainCode in unmappedCodes:
    
      chainLen = len(chainDict[chainCode].residues)
      seqDiffs[chainCode] = [999999,[]]
      
      for resChainCode in resChainCodes:
      
        resChainLen = len(resChainDict[resChainCode])
        seqDiff = abs(chainLen-resChainLen)
        
        if seqDiffs[chainCode][0] > seqDiff:
          for oldResChainCode in seqDiffs[chainCode][1]:
            if oldResChainCode in doneResChainCodes:
              doneResChainCodes.pop(doneResChainCodes.index(oldResChainCode))
          seqDiffs[chainCode][0] = seqDiff
          seqDiffs[chainCode][1] = [resChainCode]
          
          if resChainCode not in doneResChainCodes:
            doneResChainCodes.append(resChainCode)
          
        elif seqDiffs[chainCode][0] == seqDiff:
          seqDiffs[chainCode][1].append(resChainCode)
          if resChainCode not in doneResChainCodes:
            doneResChainCodes.append(resChainCode)
    
    if test:
      print "ORIGINAL seqDiff info:"
      print seqDiffs
    
    #
    # If there's problems, reset some stuff...
    #
    
    doneResChainCodes.sort()
    if doneResChainCodes != resChainCodes:
      if test:
        print "PROBLEM - resetting seqDiff info"
      for chainCode in unmappedCodes:
        if seqDiffs[chainCode][0] > 5:
          seqDiffs[chainCode][1] = resChainCodes[:]

        for otherChainCode in unmappedCodes:
          if otherChainCode == chainCode:
            continue
          if seqDiffs[chainCode][1] == seqDiffs[otherChainCode][1]:
            seqDiffs[chainCode][1] = resChainCodes[:]
            seqDiffs[otherChainCode][1] = resChainCodes[:]
            break
          
    if test:
      print unmappedCodes
      print seqDiffs

    # TODO: have a variable that decides whether to set offset if anything left???!
    #       should not set anything if in this else statement...

    molDict = {}

    for chainCode in unmappedCodes:
      molecule = chainDict[chainCode].molecule
      (seqDiff,resChainCodes) = seqDiffs[chainCode]
      if not molDict.has_key(molecule):
        molDict[molecule] = ([],[])
      molDict[molecule][0].append(chainCode)
      for resChainCode in resChainCodes:
        if resChainCode not in molDict[molecule][1]:
          molDict[molecule][1].append(resChainCode)

    if test:
      print "Moldict", molDict
      
    # TODO: should include check to see if resChainCode appears in different molecules?
    for molecule in molDict.keys():
    
      (tempChainCodes,tempResChainCodes) = molDict[molecule]
      if test:
        print "doing moldict",tempChainCodes, tempResChainCodes
      

      #
      # Try to catch ZN, FE, ... if part of original chain...
      #

      chainLen = len(molecule.molResidues)

      if chainLen == 1:
        molResidue = molecule.sortedMolResidues()[0]
        atomNames = [atom.name for atom in molResidue.chemCompVar.chemAtoms]

        if test:
          print molecule.name, atomNames

        for tempChainCode in tempChainCodes[:]:
          for tempResChainCode in tempResChainCodes:
            if test:
              print "TRYATOMMATCH",tempChainCode, tempResChainCode
            for resResidue in resChainDict[tempResChainCode]:

              fullMatch = True

              for atomName in resResidueDict[tempResChainCode][resResidue][1]: 

                if atomName not in atomNames:
                  fullMatch = False
                  break

              if fullMatch:
                if test:
                  print " Matched on atomname, 1 residue: %s-%s" % (tempChainCode,tempResChainCode)

                setDict(chainMappingDict,tempChainCode,tempResChainCode)
                
                if tempChainCode in unmappedCodes:
                  unmappedCodes.pop(unmappedCodes.index(tempChainCode))
              
                if tempChainCode in tempChainCodes:
                  tempChainCodes.pop(tempChainCodes.index(tempChainCode))

                break

      #
      # If 1-1 match, set immediately
      #
      
      if len(tempChainCodes) == len(tempResChainCodes):
        for i in range(0,len(tempChainCodes)):
          setDict(chainMappingDict,tempChainCodes[i],tempResChainCodes[i])
          if resChainCodes.count(tempResChainCodes[i]):
            resChainCodes.pop(resChainCodes.index(tempResChainCodes[i]))
            
      #
      # Otherwise try to figure out best match based on atom names...
      #

      elif len(tempChainCodes) == 1:
      
        residues = chainDict[tempChainCodes[0]].sortedResidues()
        resChainMatches = []
        bestScore = 0

        for tempResChainCode in tempResChainCodes:
          if test:
            print "tempResChainCode check: '%s'" % tempResChainCode
          resResidues = resChainDict[tempResChainCode][:]
                    
          (match,bestMatch) = doMatchResidues((None,None,-999),comparisonList,0,residues,resResidues,resResidueDict[tempResChainCode],1,test)
          (resId,offset,score) = match

          if resId != None and offset != None:
            if test:
              print tempResChainCode, score, bestScore
            if score > bestScore:
              resChainMatches = [tempResChainCode]
              bestScore = score
            elif score == bestScore:
              resChainMatches.append(tempResChainCode)

          if matchMode:
            if bestMatch[0] != None and bestMatch[1] != None:
              matchText = "start %d, offset %d, score %d" % (bestMatch[0],bestMatch[1],bestMatch[2])
            else:
              matchText = 'None'
              
            print "  Best match for CCPN code '%s' to code '%s': %s" % (tempChainCodes[0],tempResChainCode,matchText)

        if test:
          print tempChainCodes[0], resChainMatches

        if len(resChainMatches) == 1:
          setDict(chainMappingDict,tempChainCodes[0],resChainMatches[0])
          unmappedCodes.pop(unmappedCodes.index(tempChainCodes[0]))
          if resChainCodes.count(resChainMatches[0]):
            resChainCodes.pop(resChainCodes.index(resChainMatches[0]))

    # If there's only 1-1 left, map that!
    if len(unmappedCodes) == 1 and len(resChainCodes) == 1:
      setDict(chainMappingDict,unmappedCodes[0],resChainCodes[0])
                
  if not chainMappingDict:  
    print " Cannot map any chains... aborting"
    if test:
      print resChainCodes
    returnStatus = False
    
  elif test:
    print chainMappingDict
  
  return returnStatus
  
def matchResidues(resType,residues,atomNames,resResidues,resAtomDict,matchScore = 10,test = 0):

  if test:
    print resType

  moleculeMap = []
  resonanceMap = []
  
  for residue in residues:
    mapCode = 0
    if residue.ccpCode.upper() in resType:
      mapCode = 1
    moleculeMap.append(mapCode)
  
  for resResidue in resResidues:
    
    if resAtomDict[resResidue][0] and resAtomDict[resResidue][0] in resType:
      # Direct match
      mapCode = 10
    
    else:
      resAtomNames = resAtomDict[resResidue][1]
      
      mapCode = -1
      
      if resAtomNames:
        mapCode += 1
        for atomName in atomNames:
          if atomName in resAtomNames:
            mapCode += 1
    
    resonanceMap.append(mapCode)
  
  if test:  
    print moleculeMap
    print resonanceMap
  
  bestMatch = 0
  bestMatchLen = 0
  bestResStart = None

  correctFit = (None,None,0)
  bestFit = (None,None,0)
  
  maxSeqLen = max(len(residues),len(resResidues))
  moleculeMapLen = len(moleculeMap)
  resonanceMapLen = len(resonanceMap)
  
  for shift in range(-maxSeqLen,maxSeqLen):
    
    match = 0
    firstFit = 0
    offset = 0
    matchLen = 0
    resStart = 0
    
    for i in range(0,maxSeqLen):
      j = i + shift
      
      if i >= moleculeMapLen:
        continue
      if j < 0 or j >= resonanceMapLen:
        continue

      molMapValue = moleculeMap[i]
      resMapValue = resonanceMap[j]

      matchLen += 1
       
      # Ignore if no atom names!!
      if resMapValue == -1:
        continue

      if not firstFit:
        firstFit = i + 1
        offset = resResidues[j][0] - firstFit
        resStart = j

      if molMapValue != resMapValue:
        if molMapValue and resMapValue:
          match += matchScore * resMapValue
        elif molMapValue and not resMapValue:
          match -= (matchScore / 2)
        elif not molMapValue and resMapValue:
          match -= matchScore * resMapValue
        else:
          match -= 1            
      else:
        if molMapValue:
          match += matchScore
        else:
          match += 1
    
    # TODO this a problem for homodimers?!
    if match > bestMatch:
      bestMatch = match
      bestFit = (firstFit,offset,bestMatch)
      bestMatchLen = matchLen
      bestResStart = resStart
  
  if test:
    print "  endMatch %5d %-20s %5d %5d (%5d,%5d)" % (bestMatch,str(bestFit[:2]),bestMatchLen,bestResStart,moleculeMapLen,resonanceMapLen)
  
  # Reset if senseless match...
  if bestMatch > min(moleculeMapLen,resonanceMapLen):
    if test:
      print "  Resetting, bestMatch %d, from mininum of %d, %d (map lengths for molecule and resonance)" % (bestMatch,moleculeMapLen,resonanceMapLen)
    
    correctFit = bestFit
    resResiduesCopy = resResidues[:]
    for i in range(0,len(resResiduesCopy)):
      if bestResStart <= i < bestResStart+bestMatchLen:
        resResidues.pop(resResidues.index(resResiduesCopy[i]))
    
  return (bestFit,correctFit)
    

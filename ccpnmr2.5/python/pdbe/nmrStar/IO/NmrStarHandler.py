#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

NmrStarHandler.py: Contains functions specific to NmrStar conversions.

Copyright (C) 2007 Chris Penkett (European Bioinformatics Institute)

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

- contact Chris Penkett (penkett@ebi.ac.uk)
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

import traceback, sys, os, time, string, copy

from pdbe.nmrStar.IO import nmrStarDict

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccpnmr.format.general.Util import getResName
from ccpnmr.format.general.Util import getNameInfo
from ccpnmr.format.general.Util import getApplResNames
from ccpnmr.format.general.Constants import assign_kw


from memops.universal.Util import returnInt
from memops.general.Util import returnMemopsLine

from ccp.format.nmrStar.generalIO import NmrStarFile, SaveFrame
from ccp.format.nmrStar.util import getNmrStarValue

from ccp.format.general.Constants import defaultSeqInsertCode
from ccp.format.general.Util import getSeqAndInsertCode

from ccpnmr.format.general.Constants import tagSep, defaultMolCode

from ccpnmr.format.general.Util import getOriginalData, setOriginalData, getOriginalDataValue

import memops.api.Implementation as Implementation
import ccp.api.nmr.Nmr as Nmr
import ccp.api.molecule.Molecule as Molecule
import ccp.api.molecule.ChemComp as ChemComp

# Specific for nmrStar!
import ccp.api.nmr.NmrEntry as NmrEntry

# CHANGE START
#
# Import generic reader stuff
# 

#from pdbe.nmrStar.IO import generalIO
#from pdbe.nmrStar.IO import read
#from pdbe.nmrStar.IO import readStarDict
#from pdbe.nmrStar.IO import nmrStarDict
# CHANGE END

#
# Additional IOkeywords definitions
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['writeProject']['entry'] = (None,True,'The entry to be exported to the NMR-STAR file.')
IOkeywords['writeProject']['useOriginalData'] = (False,False,'Use original data for export (if available).')
IOkeywords['writeProject']['exportAll'] = (False,False,'Export all information (if available).')

#
# TODO: nmrStar is a special case... still handling on same level
# as other formats for now.
# This means that many nmrStar scripts are superfluous though...
#

class NmrStarFormat(DataFormat):

  ccpn2NmrStarMolTypes = {
        
        'protein': 'polypeptide(L)',
        'DNA':     'polydeoxyribonucleotide',
        'RNA':     'polyribonucleotide',
        'DNA/RNA': 'DNA/RNA hybrid',
        'carbohydrate': 'carbohydrates'
                        
                        }
                        
  chainIdToCcpn = {}
  
  # TODO SET THIS CORRECTLY!
  version = '3.0'

  def setFormat(self):
  
    self.format = 'nmrStar'
    self.IOkeywords = IOkeywords

  def getSequence(self):
  
    # TODO HERE: have to figure out what to do if project file read...
  
    #self.sequenceFile = self.generalIO.NmrStarFile(self.fileName)
    #self.sequenceFile.read()
    
    #
    # Some special gynmastics to make single access to sequences...
    #
    
    #self.sequenceFile.sequences = []
    
    #for sequenceFile in self.sequenceFile.sequenceFiles:
    #  for sequence in sequenceFile.sequences:
    #    self.sequenceFile.sequences.append(sequence)

    if not self.file:
      self.file = generalIO.NmrStarFile(self.fileName)
      self.file.setCategory('sequence') # This will limit the saveframes that are read
      self.file.readGeneric()
    
    componentDict = read.importData(self.file,readStarDict,nmrStarDict)
    
    # TODO: check if more than one found! Shouldn't happen...
    self.sequenceFile = componentDict[definitions.NmrStarSequenceFile][0] 
    
    self.setCcpnMolTypes()
    
    if self.verbose == 1:
      print "Reading sequence from %s file %s" % (self.formatLabel,self.fileName)


  #
  # Functions different to default functions in DataFormat
  #
  
  def createMolecule(self,molName,createMoleculeInfo):
        
    molecule = Molecule.Molecule(self.project, name = molName)
    
    if hasattr(self.sequenceFile,'commonName'):
      molecule.longName = self.sequenceFile.commonName
        
    # TODO; abbreviation doesn't fit in data model...
    
    if hasattr(self.sequenceFile,'otherNames'):
      # TODO: this not correct... how are multiple names stored in star?
      for commonName in self.sequenceFile.otherNames:
        molecule.addCommonName(returnMemopsLine(commonName))
   
    return molecule
    
  def setChainInfo(self,chain):
  
    #
    # Use this to track ID later if mapping to this chain from within the NMR-STAR file...
    #
    
    setOriginalData(self.format,chain,self.sequence,'Id')

  def getChainOrRefChainId(self,coordOrChain,isChain = False):
  
    chainCode = None
    
    if isChain:
      for chain in self.molSystem.sortedChains():
        starChainId = getOriginalDataValue(self.format,chain,'Id')

        if starChainId == coordOrChain.chainId:
          chainCode = chain.code
          self.chainIdToCcpn[starChainId] = chainCode
          coordOrChain.refChainId = chainCode 
   
    else:
      chainCode = self.chainIdToCcpn[coordOrChain.entityAssemblyId]
    
    return chainCode
   
  def setCcpnMolTypes(self):
      
    #
    # Reset original polymerType!
    #

    for sequence in self.sequenceFile.sequences:
      if hasattr(sequence,'polymerType'):
        if sequence.polymerType in self.ccpn2NmrStarMolTypes.values():
          for molType in self.ccpn2NmrStarMolTypes.keys():
            if self.ccpn2NmrStarMolTypes[molType] == sequence.polymerType:
              sequence.polymerType = molType
        # Hack, otherwise too much hassle
        elif sequence.polymerType == 'polypeptide(D)':
          sequence.polymerType = 'protein'

  def getShiftAmbiguityCode(self,chemShiftValue,resonanceToAtom,shiftList):
        
    ambCode = 1
    
    chemAtom = resonanceToAtom.chemAtom
    
    # Check prochiral and aromatic chemAtoms
    if chemAtom:
      chemAtomSet = chemAtom.chemAtomSet
      if chemAtomSet:
        # Get info for 'deep' chemAtomSets
        if chemAtomSet.isEquivalent and chemAtomSet.chemAtomSet:
          chemAtomSet = chemAtomSet.chemAtomSet
      
        if chemAtomSet.isProchiral and resonanceToAtom.otherGroupResonances:
          refResonance = resonanceToAtom.resonance
      
          for resonance in resonanceToAtom.otherGroupResonances:          
          
            shift = resonance.findFirstShift(parentList = shiftList)
        
            if shift and shift.value != chemShiftValue:
              ambCode = 2
              break
          
        elif chemAtomSet.isProchiral == False and chemAtomSet.isEquivalent == None:
          # Should be aromatic!
          ambCode = 3
          
    return ambCode

  def getPresetChainMapping(self,chainList):
  
    mappingChainDict = {}
    
    for chain in chainList:
    
      mappingChainDict[chain] = (chain.code,1)
    
    return mappingChainDict


class NmrStarFullReaderFile(NmrStarFile):

  def setGeneric(self):
  
    self.format = 'nmrStar'
    self.defaultMolCode = defaultMolCode
    self.version = '3.0' # is default
    self.tagSep = '.'
    self.tagStart = '_'
  
  def initialize(self):
  
    self.components = ['all']
    self.setComponents()

  # TODO HOW CAN I BUILD THIS IN?
  #def setCategory(self,category):
  
  #  self.components = [] # TODO: or append to existing? No... should be addCategory...
    
  #  for sfCatName in readStarDict.sfList:
  #    if readStarDict.sfDict[sfCatName].has_key('category') and readStarDict.sfDict[sfCatName]['category'] == category:
  #      self.components.append(sfCatName)
        
  #  self.setComponents()
      
  def setSfDict(self,component):
  
    #
    # Keep track of component order
    #

    self.componentList = []

    #
    # Generic R/W only allows 3.0 (and later?)
    #

    ###################################
    #                                 #
    # START OF VERSION 3.0 STUFF!!    #
    #                                 #
    ###################################

    if self.version == '3.0':
      
      """
      
      NOTE: all Jurgen's stuff will be added to the dictionary!
      
      """
      
      for sfName in nmrStarDict.sfList:

        if component == sfName or component == 'all':
 
          self.componentList.append(sfName)
          self.sfDict[sfName] = nmrStarDict.sfDict[sfName]

    else:
    
      print "ERROR UNRECOGNIZED VERSION '%s'" % (str(self.version))

  
  def setupSaveFrame(self,saveFrameName,title):
  
    keywds = {}
  
    if not self.sfDict.has_key(saveFrameName):
      print "  Warning: saveframe name %s not in reference data!" % (saveFrameName)

    else:
      keywds['prefix'] = self.getPrefix(saveFrameName)

    if not self.sfs.has_key(saveFrameName):
      self.sfs[saveFrameName] = []
      
    self.sfs[saveFrameName].append(SaveFrame(title,saveFrameName,**keywds))
    
    return self.sfs[saveFrameName][-1]

  def getPrefix(self,saveFrameName):
  
    prefix = self.tagStart + self.sfDict[saveFrameName]['name']
  
    return prefix
    
  #
  # Redefinitions for read stuff
  #

  def setVersion(self,versionHits):
  
    if versionHits['3.0'] > versionHits['2.1.1'] and self.version != '3.0':
     
      print "  Warning: setting nmrStar version to 3.0 for reading."
      self.version = '3.0'
         
    elif versionHits['2.1.1'] > versionHits['3.0'] and self.version != '2.1.1':
     
      print "  Error: no full version 2.1.1 nmrStar reader available - use version 3.0 files."
      return False
    
    return True


  def setReadTags(self,saveFrame,saveFrameDict,tagtable,convError):

    for tagName in saveFrameDict['tagNames']:

      returnValueFunc = saveFrameDict['tags'][tagName][1]

      if self.version == '3.0': # TODO THIS HAS TO BE MODIFIED!! Or is irrelevant

        matchTagName = self.tagStart + saveFrameDict['name'] + self.tagSep + tagName

      starValue = getNmrStarValue(tagtable,matchTagName)

      if returnValueFunc:
        value = returnValueFunc(starValue)
      else:
        value = starValue

      saveFrame.setTag(matchTagName,value,error = convError.getString())

  def setTableTags(self,saveFrame,saveFrameDict,tagtable,convError):

    for tableName in saveFrameDict['tableNames']:

      matchTableName = self.tagStart + tableName

      tableExists = 0

      if self.version == '3.0':

        (tempTableName,tagName) = string.split(tagtable.tagnames[0],self.tagSep)

        if matchTableName == tempTableName:
          tableExists = 1

      if tableExists:

        table = saveFrame.setupTable(matchTableName)
        currentTags = []

        for tagName in saveFrameDict['tables'][tableName]['tagNames']:

          matchTagName = tagName

          if self.version == '3.0':

            matchTagName = matchTableName + self.tagSep + tagName

          if matchTagName in tagtable.tagnames:

            returnValueFunc = saveFrameDict['tables'][tableName]['tags'][tagName][1]
            currentTags.append([matchTagName,returnValueFunc])

            table.tags[matchTagName] = []
            table.tagErrors[matchTagName] = []
            table.tagNames.append(matchTagName)

        for i in range(0,len(tagtable.tagvalues[0])):

          for (matchTagName,returnValueFunc) in currentTags:

            starValue = getNmrStarValue(tagtable,matchTagName,i)

            if returnValueFunc:
              value = returnValueFunc(starValue)
            else:
              value = starValue

            table.setTag(matchTagName,value,error = convError.getString())

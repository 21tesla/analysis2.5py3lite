#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

PdbFormat.py: Contains functions specific to PDB conversions.

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

import copy

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccp.format.general.Constants import defaultSeqInsertCode

#
# Add some information to IOkeywords...
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readCoordinates']['mergeFunctionInfo'] = (None,False,"Tuple with function to merge in external info and file with source info.")
IOkeywords['writeCoordinates']['addPdbHeader'] = (False,False,"Add a header with fake crystal information for wwPDB deposition purposes.")

class PdbFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'pdb'
    self.IOkeywords = IOkeywords


  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    
    self.getCoordinates = self.getCoordinatesGeneric
    self.createCoordinateFile = self.createCoordinateFileGeneric 

    self.getChemComps = self.getChemCompsGeneric

  #
  # Deviations from generic import stuff
  #
      
  def getSequenceSetFormatSpecificReadKeywds(self):  
   
    readKeywds = {}

    if self.ignoreResNames:
      readKeywds['ignoreResNames'] = self.ignoreResNames
    
    return readKeywds

  def getCoordinatesSetFormatSpecificReadKeywds(self):
   
    readKeywds = {}
    
    if self.ignoreResNames:
      readKeywds['ignoreResNames'] = self.ignoreResNames
    
    return readKeywds
    
  def getCoordinatesFormatSpecific(self):
   
    if self.mergeFunctionInfo:
      self.mergeFunctionInfo[0](self.coordinateFile,self.mergeFunctionInfo[1])

  def getCoordinatesSetSequenceFile(self):
    
    self.setSequenceFileClass()
    self.sequenceFile = self.SequenceFileClass(self.fileName)

    readKeywds = self.getSequenceSetFormatSpecificReadKeywds()
    self.sequenceFile.read(pdbFile = self.coordinateFile, **readKeywds)

  def getPeopleAndCitations(self):
  
    # TODO HERE: have to figure out what to do if project file read...
  
    self.peopleAndCitationsFile = self.peopleAndCitationsIO.PdbPeopleAndCitationsFile(self.fileName)
    self.peopleAndCitationsFile.read()

    if self.verbose == 1:
      print "Reading people and citation info from %s file %s" % (self.formatLabel,self.fileName)

  #
  # Code different from DataFormat
  #

  def getHeader(self):

    self.headerFile = self.generalIO.PdbFile(self.fileName)
    self.headerFile.readHeader()

    if self.verbose == 1:
      print "Reading header info from %s file %s" % (self.formatLabel,self.fileName)

  def forceNonEmptyChainCode(self):
    
    returnValue = False
    
    if self.sequenceFile.version >= '3.0':
      returnValue = True
    
    return returnValue

  def getCode3Letter(self,sequenceElement):
    
    code3Letter = sequenceElement.code3Letter.strip().upper()
    
    # TODO THIS IS A HACK TO CONVERT +T, +C, .. (signifying 'a modified nucleotide') to T, C... .
    # This should be cleared up in the actual PDB files eventually
    if code3Letter[0] == '+' and code3Letter[1] in 'CGTUAI':
      code3Letter = code3Letter[1:] 

    return code3Letter

  def duplicateChain(self):
  
    #
    # All chain info should be contained in PDB files...
    #
    
    return 0   

  def getSpecificPdbOneLetterCode(self,currentChainCode):
    
    pdbOneLetterCode = currentChainCode
    
    if self.sequence:
      pdbOneLetterCode = self.sequence.chainCode
      
    # Currently not possible to be more than one letter, so reset to nothing.
    # This restriction will likely change in the future.
    if len(pdbOneLetterCode) > 1:
      pdbOneLetterCode = ' '
  
    return pdbOneLetterCode

  def setCoordChainList(self):
  
    self.coordChainList = self.coordinateFile.chains + self.coordinateFile.hetGroups

  def createCoordinateFileFormatSpecific(self):
    
    if self.addPdbHeader:
      self.writeKeywds['addPdbHeader'] = self.addPdbHeader

  def resetCoordinateSerial(self,coordinateSerial):
    
    # Can use this to reset the coordinate serial used in the output file on a per-model basis
    if self.version not in ('2.2',):
      coordinateSerial = 0
    
    return coordinateSerial
         
  def setRawCoordinate(self):
  
    if not self.coordinateFile.modelCoordinates.has_key(self.modelId):

      self.coordinateFile.modelCoordinates[self.modelId] = []
   
    #
    # Get correct code (use cifCode) and set hetFlag
    #
    
    hetFlag = 0
    
    chemCompVar = self.residue.chemCompVar
    chemComp = chemCompVar.chemComp
    cifCode = self.findChemCompVarSysName('CIF',chemCompVar)
    
    if chemComp.className == 'NonStdChemComp':
      hetFlag = 1
    
    if self.namingSystemName in ('PDB','PDB_REMED') and cifCode:
      resName = cifCode
    elif self.namingSystemName in ('PDB','PDB_REMED') and chemComp.code3Letter:
      resName = chemComp.code3Letter
    else:
      resName = chemComp.ccpCode.upper()

    if self.useOriginalData and self.origChainCode != None:
      chainCode = self.origChainCode
    else:
      chainCode = self.exportChainCode
  
    if self.useOriginalData and self.origSeqCode != None:
      seqCode = self.origSeqCode
    else:
      seqCode = self.seqCode
 
    if self.useOriginalData and self.origSeqInsertCode != None:
      seqInsertCode = self.origSeqInsertCode
    else:
      seqInsertCode = defaultSeqInsertCode
 
    if self.useOriginalData and self.origAtomName != None and not self.forceNamingSystemName:
      atomName = self.origAtomName
    else:
      atomName = self.atomName

    modelCoordinate = self.coordinatesIO.PdbCoordinate(self.coordinateFile,
        self.coordinateSerial,atomName,resName,chainCode,seqCode,
        self.x,self.y,self.z,self.occupancy,self.bFactor,
        self.atomType,hetFlag,insertionCode=seqInsertCode)
      
    self.coordinateFile.modelCoordinates[self.modelId].append(modelCoordinate)

  def getPresetChainMapping(self,chainList):
  
    mappingChainDict = {}
    
    for chain in chainList:
    
      mappingChainDict[chain] = (chain.code,1)
    
    return mappingChainDict

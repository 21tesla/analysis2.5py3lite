"""
======================COPYRIGHT/LICENSE START==========================

coordinatesIO.py: I/O for nmrStar coordinate saveframe

Copyright (C) 2005-2008 Wim Vranken (European Bioinformatics Institute)

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

import os

from ccp.format.nmrStar.generalIO import NmrStarFile
from ccp.format.nmrStar.generalIO import NmrStarGenericFile

from ccp.format.general.Constants import defaultSeqInsertCode
from ccp.format.general.Util import getSeqAndInsertCode

from memops.universal.Util import returnInt

#####################
# Class definitions #
#####################

class NmrStarFile(NmrStarFile):

  def initialize(self, version = '2.1.1'):
  
    self.coordinateFiles = []
    self.files = self.coordinateFiles

    if not self.version:
      self.version = version
    
    self.saveFrameName = 'conformer_family_coord_set'
    self.DataClassFile = NmrStarCoordinateFile
    
    self.components = [self.saveFrameName]
    self.setComponents()

  def read(self,verbose = 0):
    
    self.readComponent(verbose = verbose) 

class NmrStarCoordinateFile(NmrStarGenericFile):

  def initialize(self,parent,saveFrame = None,maxModelNum = 999):

    self.saveFrame = saveFrame
    
    self.parent = parent
    self.version = parent.version

    self.chains = []
    self.hetGroups = []
    self.modelCoordinates = {}
    
    self.code = None
    self.details = None

    # Set maximum number of models to read. Won't work with project level reading...    
    if not hasattr(self,'maxModelNum'):
      self.maxModelNum = maxModelNum

    self.ssBonds = []

    if self.saveFrame:
      self.parseSaveFrame()

  def parseSaveFrame(self):
       
    if not self.checkVersion():
      return
      
    #
    # Here do parsing and reading...
    #

    tableTagNames = {}
    chainIds = []
    
    if self.version == '2.1.1':

      return
      
    elif self.version == '3.0':
    
      tableName = '_Atom_site'
      tableTagNames['modelId'] = 'Model_ID'
      tableTagNames['Id'] = 'ID'
      tableTagNames['entityAssemblyId'] = 'Label_entity_assembly_ID'
      tableTagNames['entityId'] = 'Label_entity_ID'
      tableTagNames['indexId'] = 'Label_comp_index_ID'
      tableTagNames['compId'] = 'Label_comp_ID'
      tableTagNames['atomId'] = 'Label_atom_ID'
      tableTagNames['authorChainId'] = 'Auth_segment_code'
      tableTagNames['authorSeqCode'] = 'Auth_seq_ID'
      tableTagNames['authorResName'] = 'Auth_comp_ID'
      tableTagNames['authorAtomName'] = 'Auth_atom_ID'
      tableTagNames['atomType'] = 'Type_symbol'
      tableTagNames['x'] = 'Cartn_x'
      tableTagNames['y'] = 'Cartn_y'
      tableTagNames['z'] = 'Cartn_z'

    elif self.version == '3.1':
    
      tableName = '_Atom_site'
      tableTagNames['modelId'] = 'Model_ID'
      tableTagNames['Id'] = 'ID'
      tableTagNames['entityAssemblyId'] = 'Label_entity_assembly_ID'
      tableTagNames['entityId'] = 'Label_entity_ID'
      tableTagNames['indexId'] = 'Label_comp_index_ID'
      tableTagNames['compId'] = 'Label_comp_ID'
      tableTagNames['atomId'] = 'Label_atom_ID'
      tableTagNames['atomType'] = 'Type_symbol'
      tableTagNames['x'] = 'Cartn_x'
      tableTagNames['y'] = 'Cartn_y'
      tableTagNames['z'] = 'Cartn_z'
      tableTagNames['occupancy'] = 'Occupancy'
      tableTagNames['authorChainId'] = 'Auth_asym_ID'
      tableTagNames['authorSeqCode'] = 'Auth_seq_ID'
      tableTagNames['authorResName'] = 'Auth_comp_ID'
      tableTagNames['authorAtomName'] = 'Auth_atom_ID'
      tableTagNames['pdbChainId'] = 'PDB_strand_ID'
      tableTagNames['pdbSeqCode'] = 'PDB_residue_no'
      tableTagNames['pdbSeqInsertCode'] = 'PDB_ins_code'
      tableTagNames['pdbResName'] = 'PDB_residue_name'
      tableTagNames['pdbAtomName'] = 'PDB_atom_name'
      tableTagNames['pdbxChainId'] = 'PDBX_label_asym_ID'

    if self.saveFrame.tags.has_key('Details'):
      self.details = self.saveFrame.tags['Details']

    if self.saveFrame.tables.has_key(tableName):

      tableTags = self.saveFrame.tables[tableName].tags

      for i in range(0,len(tableTags[tableTagNames['Id']])):

        Id = tableTags[tableTagNames['entityAssemblyId']][i]
        modelId = tableTags[tableTagNames['modelId']][i]
        cifCode = tableTags[tableTagNames['compId']][i]
        indexId = tableTags[tableTagNames['indexId']][i]
        
        if modelId > self.maxModelNum:
          continue

        if not self.modelCoordinates.has_key(modelId):

          self.modelCoordinates[modelId] = []
        
        # Waters are handled differently - need to distinguish by sequence code as well.
        if cifCode == 'HOH':
          chainId =  (Id,indexId)
          isWaterChain = True
        else:
          chainId = Id
          isWaterChain = False
        
        if not chainId in chainIds:
          
          self.chains.append(NmrStarChain(Id,indexId, isWaterChain = isWaterChain))
          chainIds.append(chainId)

        self.modelCoordinates[modelId].append(NmrStarCoordinate())     
        self.modelCoordinates[modelId][-1].setData(tableTags,i,tableTagNames)
        self.modelCoordinates[modelId][-1].setDefaultAttrs()

      tableName1 = '_Conformer_family_refinement'
      self.setStructMethods(tableName1)

      tableName2 = '_Conformer_family_software'
      self.setStructSoftwares(tableName2)

class NmrStarChain:

  def __init__(self,Id,firstIndexId,isWaterChain = False):

    self.chainId = Id
    self.firstIndexId = firstIndexId
    
    self.isWaterChain = isWaterChain
    
class NmrStarCoordinate:

  def __init__(self):
  
    self.Id = None
    self.modelId = None
    
    self.entityAssemblyId = None
    self.entityId = None
    self.indexId = None
    self.compId = None
    self.atomId = None
    
    self.atomName = None
    self.resName = None
    self.chainCode = None
    self.seqCode = None
    self.insertionCode = defaultSeqInsertCode
    
    ##################
    # All this stuff is handled separately in readCoordinates!!!
    self.authorChainId = None
    self.authorSeqCode = None
    self.authorResName = None
    self.authorAtomName = None
    ##################
    
    ##################
    # Necessary for PDB info transfer!!!
    self.pdbChainId = None
    self.pdbSeqCode = None
    self.pdbSeqInsertCode = None
    self.pdbResName = None
    self.pdbAtomName = None
    ##################

    self.origData = None
    
    self.x = None
    self.y = None
    self.z = None

    self.bFactor = None
    self.occupancy = 1.0
    self.atomType = None

  def setData(self,tableTags,i,tableTagNames):
  
    for attrName in tableTagNames.keys():
      if tableTags.has_key(tableTagNames[attrName]):
        setattr(self,attrName,tableTags[tableTagNames[attrName]][i])
      
  def setDefaultAttrs(self):
  
    self.atomName = self.atomId
    self.seqCode = self.indexId
    self.chainCode = self.entityAssemblyId
    self.resName = self.compId
    

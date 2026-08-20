"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for nmrStar sequence information saveframes

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

from ccp.format.nmrStar.generalIO import NmrStarFile, NmrStarGenericFile

#from ccp.format.nmrStar.projectIO import MoleculeDbLink

from constants import bmrbCodeToCcpCode

from memops.universal.Util import returnInt

#####################
# Class definitions #
#####################

class NmrStarFile(NmrStarFile):

  def initialize(self, version = '2.1.1'):
  
    self.sequenceFiles = []    
    
    self.files = self.sequenceFiles

    if not self.version:
      self.version = version
          
    self.DataClassFile = NmrStarSequenceFile
    
    self.setComponents()

  def setComponentList(self):
  
    otherSfs = []

    if self.version in ('3.0','3.1'):
      self.saveFrameName = 'assembly'
      otherSfs = ['entity','chem_comp']
    else:
      self.saveFrameName = 'monomeric_polymer'
      otherSfs = ['molecular_system']

    self.components = [self.saveFrameName] + otherSfs    

  def read(self,verbose = 0):

    self.readComponent(verbose = verbose)
    
class NmrStarSequenceFile(NmrStarGenericFile):

  """
  Information on file level.
  """
  
  def initialize(self,parent,saveFrame = None):

    self.sequences = []
    self.dbLinks = []
    self.bioFunctions = []
    
    self.saveFrame = saveFrame
    
    self.molSystemCode = None
    self.details = None
    self.EC_number = None
    
    self.parent = parent
    self.version = parent.version

    if self.saveFrame:
      self.parseSaveFrame()

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    if self.version in ('3.0','3.1'): 

      self.molSystemCode = self.saveFrame.tags['Name']

      if self.saveFrame.tags.has_key('Details'):
        self.details = self.saveFrame.tags['Details']

      if self.saveFrame.tags.has_key('Enzyme_commission_number'):
        self.EC_number = self.saveFrame.tags['Enzyme_commission_number']

      if self.saveFrame.tables.has_key('_Entity_assembly'):

        #
        # First look for original sequence codes if available from coordinate section. Hack for reading combined NMR-STAR files
        #
        
        if self.version in ('3.1',):
          
          self.origSeqCodeBySeqId = {}
          self.origSeqInsertCodeBySeqId = {}
          self.atomNamesBySeqId = {}
          
          if self.parent.sfs.has_key('conformer_family_coord_set'):
            for coordSf in self.parent.sfs['conformer_family_coord_set']:
              if coordSf.tables.has_key('_Atom_site'):
                coordTableTags = coordSf.tables['_Atom_site'].tags
                for i in range(len(coordTableTags['Model_ID'])):
                  chainId  = coordTableTags['Label_entity_assembly_ID'][i]
                  seqId    = coordTableTags['Label_comp_index_ID'][i]
                  atomName = coordTableTags['Label_atom_ID'][i]
                  
                  if chainId and seqId:
                    residueKey = (chainId,seqId)
                    if not self.origSeqCodeBySeqId.has_key(residueKey):
                      if coordTableTags.has_key('Auth_seq_ID'):
                        seqCode = coordTableTags['Auth_seq_ID'][i]
                      elif coordTableTags.has_key('PDB_residue_no'):
                        seqCode = coordTableTags['PDB_residue_no'][i]
                      
                      if seqCode != None:
                        self.origSeqCodeBySeqId[residueKey] = returnInt(seqCode)

                    if not self.origSeqInsertCodeBySeqId.has_key(residueKey) and coordTableTags.has_key('PDB_ins_code'):
                      seqInsertCode = coordTableTags['PDB_ins_code'][i]
                      if seqInsertCode:
                        self.origSeqInsertCodeBySeqId[residueKey] = seqInsertCode
                    
                    if atomName:
                      if not self.atomNamesBySeqId.has_key(residueKey):
                        self.atomNamesBySeqId[residueKey] = []
                    
                      if not atomName in self.atomNamesBySeqId[residueKey]:
                        self.atomNamesBySeqId[residueKey].append(atomName)

        #
        # Now start to read the sequence info
        #
  
        tableTags = self.saveFrame.tables['_Entity_assembly'].tags

        for i in range(len(tableTags['Entity_ID'])):

          entityName = tableTags['Entity_label'][i]
          entityAssemblyId =  tableTags['ID'][i]
          
          if self.version in ('3.1',):
            originalChainCode = tableTags['Asym_ID'][i]
            entityId = tableTags['Entity_ID'][i]
            details = tableTags['Details'][i]

            if tableTags.has_key('Role'):
              role = tableTags['Role'][i]
            else:
              role = None

            if tableTags.has_key('Physical_state'):
              physState = tableTags['Physical_state'][i]
            else:
              physState = None

            if tableTags.has_key('Conformational_isomer'):
              confIsomer = tableTags['Conformational_isomer'][i]
            else:
              confIsomer = None

            if tableTags.has_key('Chemical_exchange_state'):
              chemExchState = tableTags['Chemical_exchange_state'][i]
            else:
              chemExchState = None

            if tableTags.has_key('Magnetic_equivalence_group_code'):
              magnEquivCode = tableTags['Magnetic_equivalence_group_code'][i]
            else:
              magnEquivCode = None

          else:
            originalChainCode = None
            entityId = None
          
          for entitySaveFrame in self.parent.sfs['entity']:
            entityFound = False
            if self.version in ('3.0','3.1') and entitySaveFrame.tags['ID'] == entityId:
              entityFound = True
            elif entitySaveFrame.title == entityName[1:]:
              entityFound = True
              
            if entityFound:

              #waterMoleculeInfo = self.getWaterMoleculeInfo(entityAssemblyId,entityId)
              
              #if waterMoleculeInfo:
              #  for compIndexId in waterMoleculeInfo:
              #    self.sequences.append(NmrStarSequence(entitySaveFrame.title,self,entitySaveFrame,Id = entityAssemblyId, entityId = entityId, originalChainCode = originalChainCode, waterMoleculeCompIndexId = compIndexId, details = details, role = role, physState = physState, confIsomer = confIsomer, chemExchState = chemExchState, magnEquivCode = magnEquivCode))
                  
              #else:
                self.sequences.append(NmrStarSequence(entitySaveFrame.title,self,entitySaveFrame,Id = entityAssemblyId, entityId = entityId, originalChainCode = originalChainCode, details = details, role = role, physState = physState, confIsomer = confIsomer, chemExchState = chemExchState, magnEquivCode = magnEquivCode))
                break

        #
        # Set molecular bonds
        #
        
        tableName = '_Bond'
        
        if self.saveFrame.tables.has_key(tableName):
          tableTags = self.saveFrame.tables[tableName].tags
          for i in range(len(tableTags['Type'])):
            
            bondType = tableTags['Type'][i]
            # TODO currently not handled!
            ###bondOrder = tableTags['Order'][i]
            
            sequenceElementData = []
            
            for j in range(1,3):
              for sequence in self.sequences:
                if sequence.Id == tableTags['Entity_assembly_ID_%d' % j][i]:
                  for sequenceElement in sequence.elements:
                    if sequenceElement.seqCode == tableTags['Comp_index_ID_%d' % j][i]:
                      sequenceElementData.append((sequenceElement,tableTags['Atom_ID_%d' % j][i]))
                  
                  break
            
            if len(sequenceElementData) == 2:
              sequenceElementData[0][0].setBond(bondType,sequenceElementData[0][1],sequenceElementData[1][0],sequenceElementData[1][1])
              sequenceElementData[1][0].setBond(bondType,sequenceElementData[1][1],sequenceElementData[0][0],sequenceElementData[0][1])
            else:
              print "  Warning: could not set bond information from NMR-STAR file, information missing."

        
        tableName = '_Assembly_db_link'

        tagNames = ('Database_code','Accession_code','Entry_experimental_method','Entry_structure_resolution')

        if self.saveFrame.tables.has_key(tableName):
          tableTags = self.saveFrame.tables[tableName].tags

          for i in range(0,len(tableTags['Database_code'])):

            infoList = []
            for tagName in tagNames:
              infoList.append(tableTags[tagName][i])           

            self.dbLinks.append(MolSystemDbLink(self,*infoList))

        tableName = '_Assembly_bio_function'

        tagNames = ('Biological_function',)

        if self.saveFrame.tables.has_key(tableName):
          tableTags = self.saveFrame.tables[tableName].tags

          for i in range(0,len(tableTags['Biological_function'])):

            infoList = []
            for tagName in tagNames:
              infoList.append(tableTags[tagName][i])           

            self.bioFunctions.append(MolSystemBioFunction(self,*infoList))

    elif self.version == '2.1.1':
      
      # Handled completely differently from above! Also possible that molecular_system sf not present, will deal with this.

      moleculeName = self.saveFrame.title
      
      if self.parent.sfs.has_key('molecular_system'):    
        molSysSaveFrame = self.parent.sfs['molecular_system'][0]
        tableTags = molSysSaveFrame.tables['_Mol_system_component_name'].tags

        for i in range(len(tableTags['_Mol_system_component_name'])):
          
          chainName = tableTags['_Mol_system_component_name'][i]
          moleculeLabel = tableTags['_Mol_label'][i]
          
          if moleculeLabel[1:] == moleculeName:
            self.sequences.append(NmrStarSequence(self.name,self,self.saveFrame, chainName = chainName))
      
      else:
        self.sequences.append(NmrStarSequence(self.name,self,self.saveFrame, chainName = self.defaultMolCode))
        
  def getWaterMoleculeInfo(self, entityAssemblyId, entityId):
    
    #
    # Special treatment for water molecules
    #
    
    waterMoleculeInfo = []
    
    tableName = '_PDBX_nonpoly_scheme'    
    if self.version in ('3.1',) and self.saveFrame.tables.has_key(tableName):

      tableTags = self.saveFrame.tables[tableName].tags

      for i in range(len(tableTags['Entity_ID'])):

        if entityAssemblyId == tableTags['Entity_assembly_ID'][i] and entityId == tableTags['Entity_ID'][i]:
          
          if tableTags['Comp_ID'][i] == 'HOH':
            waterMoleculeInfo.append(tableTags['Comp_index_ID'][i])
    
    return waterMoleculeInfo
          
class NmrStarSequence:

  def __init__(self,name,parent,sequenceSaveFrame,Id = None,entityId = None, originalChainCode = None, waterMoleculeCompIndexId = None, details = None, role = None, physState = None, confIsomer = None, chemExchState = None, magnEquivCode = None, chainName = None):
  
    self.elements = []
    self.dbLinks = []
    self.bioFunctions = []
    
    self.isCircular = None
    
    self.saveFrame = sequenceSaveFrame

    self.parent = parent
    self.version = self.parent.version
    
    self.Id = Id
    self.entityId = entityId
    
    # 2.1.1 only
    self.chainName = chainName
    
    self.originalChainCode = originalChainCode
    self.waterMoleculeCompIndexId = waterMoleculeCompIndexId
    self.details = details
    self.role = role

    #print 'DATA: [%s] [%s] [%s]' % (physState, confIsomer, chemExchState)
    self.physState = physState
    self.confIsomer = confIsomer
    self.chemExchState = chemExchState

    self.magnEquivCode = magnEquivCode

    tagNames = {}
    tableTagNames = {}
    
    if self.version == '2.1.1':
      self.molName = self.formatMoleculeName(sequenceSaveFrame.title)
      
      tagNames['molType'] = '_Mol_type'
      tagNames['polymerType'] = '_Mol_polymer_class'
      tagNames['commonName'] = '_Name_common'
      tagNames['otherNames'] = '_Name_variant'
      tagNames['abbrev'] = '_Abbreviation_common'
      tagNames['molMass'] = '_Molecular_mass'

      tableName = '_Residue_seq_code'
      
      seqCodeTag = seqIdTag = '_Residue_seq_code'
      authorSeqCodeTag = '_Residue_author_seq_code'
      
      if ( self.saveFrame.tables.has_key(tableName) and
           self.saveFrame.tables[tableName].tags.has_key(authorSeqCodeTag)):
        
        validCode = 1
        sequenceTableTags = self.saveFrame.tables[tableName].tags
        for i in range(len(sequenceTableTags[authorSeqCodeTag])):
        
          if sequenceTableTags[authorSeqCodeTag][i] == None:
            validCode = 0
            break
            
        if validCode:
          seqCodeTag = authorSeqCodeTag
      
      tableTagNames['seqId'] = seqIdTag # This one to track the sequential numbering used in chem shift section
      tableTagNames['seqCode'] = seqCodeTag
      tableTagNames['code3Letter'] = '_Residue_label'
    
    elif self.version == '3.0':
      self.molName = self.formatMoleculeName(name)
    
      tagNames['molType'] = 'Type'
      tagNames['polymerType'] = 'Pol_type'
      #tagNames['commonName'] = '_Name_common'
      #tagNames['otherNames'] = '_Name_variant'
      #tagNames['abbrev'] = '_Abbreviation_common'
      #tagNames['molMass'] = '_Molecular_mass'
      
      """
      From RCSB:
      
      Data items in the ENTITY_POLY_SEQ category specify the sequence
      of monomers in a polymer. Allowance is made for the possibility
      of microheterogeneity in a sample by allowing a given sequence
      number to be correlated with more than one monomer id - the
      corresponding ATOM_SITE entries should reflect this
      heterogeneity.      
      """

      # TODO: just use entity_poly_seq for now - should be OK...
      tableName = '_Entity_poly_seq' 
      tableTagNames['seqCode'] = 'Comp_index_num'
      tableTagNames['code3Letter'] = 'Mon_ID'
    
    elif self.version == '3.1':
    
      if sequenceSaveFrame.tags['Name'] and sequenceSaveFrame.tags['Name'] != '.':
        self.molName = self.formatMoleculeName(sequenceSaveFrame.tags['Name'])
      else:
        self.molName = self.formatMoleculeName(name)
    
      tagNames['molType'] = 'Type'
      tagNames['polymerType'] = 'Polymer_type'
      tagNames['Mutation'] = 'Mutation'
      tagNames['EC_number'] = 'EC_number'
      tagNames['Fragment'] = 'Fragment'
      tagNames['Src_method'] = 'Src_method'
      #tagNames['commonName'] = '_Name_common'
      #tagNames['otherNames'] = '_Name_variant'
      #tagNames['abbrev'] = '_Abbreviation_common'
      #tagNames['molMass'] = '_Molecular_mass'

      tagNames['details'] = 'Details'
      tagNames['seqDetails'] = 'Polymer_author_seq_details'
      tagNames['authDefSeq'] = 'Polymer_author_defined_seq'

      tagNames['ambConfStates'] = 'Ambiguous_conformational_states'
      tagNames['ambChemCompSites'] = 'Ambiguous_chem_comp_sites'

      if sequenceSaveFrame.tags['Type'] in ('non-polymer','water'):
        tableName = '_Entity_comp_index' 
        tableTagNames['seqCode'] = 'ID'
        tableTagNames['code3Letter'] = 'Comp_ID'
      else:
        tableName = '_Entity_comp_index' 
        tableTagNames['seqCode'] = 'ID'
        tableTagNames['code3Letter'] = 'Comp_ID'
        tableTagNames['authorSeqCode'] = 'Auth_seq_ID'

        # Eldon thinks we can use the _Entity_comp_index table instead
        # and then also get the author residue numbers.

        #tableName = '_Entity_poly_seq' 
        #tableTagNames['seqCode'] = 'Comp_index_ID'
        #tableTagNames['code3Letter'] = 'Mon_ID'

    #
    # Set tag names
    #
    
    for attrName in tagNames.keys(): 
      tagName = tagNames[attrName]
      if sequenceSaveFrame.tags.has_key(tagName):
        if sequenceSaveFrame.tags[tagName] != None:
          setattr(self,attrName,sequenceSaveFrame.tags[tagName])
        else:
          setattr(self,attrName,None)

    #
    # Set table tags (sequence elements)
    #
    
    resettingSeqCodes = False

    if self.saveFrame.tables.has_key(tableName):

      sequenceTableTags = self.saveFrame.tables[tableName].tags
      curSeqInsertCode = None
      
      for i in range(len(sequenceTableTags[tableTagNames['seqCode']])):
        
        if self.version == '3.1' and sequenceSaveFrame.tags['Type'] == 'water' and waterMoleculeCompIndexId:
          
          # TODO: HERE SHOULD BE FINDING CODE FROM COORDINATE SECTION! Use self.parent.origSeqCodeBySeqId
          if sequenceTableTags[tableTagNames['seqCode']][i] != waterMoleculeCompIndexId:
            continue
        
        # This is an on-the-fly conversion...
        bmrbCode = sequenceTableTags[tableTagNames['code3Letter']][i]
        ccpMapping = None
        
        if bmrbCode in bmrbCodeToCcpCode.keys():
          ccpMapping = bmrbCodeToCcpCode[bmrbCode]
          
          if ccpMapping[0]:
            self.elements.append(NmrStarSequenceElement(sequenceTableTags,i,tableTagNames))
            self.elements[-1].code3Letter = ccpMapping[0][1]
            self.elements[-1].residueType = ccpMapping[0][0] # This is usually not there!
        
        self.elements.append(NmrStarSequenceElement(sequenceTableTags,i,tableTagNames))
        
        #
        # Reset seqCode if can get it from coordinate section
        #
        
        if self.version == '3.1':
          
          residueKey = (self.Id,sequenceTableTags[tableTagNames['seqCode']][i])
          
          #
          # Now set the sequence code...
          #
          
          if self.parent.origSeqCodeBySeqId:
            #sks = self.parent.origSeqCodeBySeqId.keys()
            #sks.sort()
            #  print sk, self.parent.origSeqCodeBySeqId[sk]
            #for sk in sks:
            #sys.exit()
            if self.parent.origSeqCodeBySeqId.has_key(residueKey):
    
              #
              # Set sequence insertion code specifically for this residue
              #
              
              if self.parent.origSeqInsertCodeBySeqId and self.parent.origSeqInsertCodeBySeqId.has_key(residueKey):
                self.elements[-1].insertionCode = self.parent.origSeqInsertCodeBySeqId[residueKey]
                curSeqInsertCode = self.parent.origSeqInsertCodeBySeqId[residueKey]
              else:
                curSeqInsertCode = None


              resettingSeqCodes = True
              self.elements[-1].seqCode = self.parent.origSeqCodeBySeqId[residueKey]
              #print "RESET ORIG for %s to %d" % (residueKey,self.elements[-1].seqCode)
            
              # Also try and reset previous sequence codes based on coordinate section
              # This is to prevent overlapping sequence codes, and makes code work better in general.
              # Only do this if no sequence insertion code was used in the previous residue...
              for i in range(1,99999):
                if len(self.elements) <= i:
                  break
                prevResidueKey = (residueKey[0],residueKey[1] - i)                
                if not self.parent.origSeqCodeBySeqId.has_key(prevResidueKey):
                  self.elements[-(i+1)].seqCode = self.elements[-1].seqCode - i
                  #print "RESET RANGE", i, residueKey[1]-i, self.elements[-(i+1)].seqCode
                else:
                  break
          
            elif resettingSeqCodes:          
              prevSeqCode = self.elements[-2].seqCode
              self.elements[-1].seqCode = prevSeqCode + 1
              #print "RESET PREV", (prevSeqCode + 1)
          
            if curSeqInsertCode:
              self.elements[-1].insertionCode = curSeqInsertCode

          #
          # Also set the atom names
          #
          
          if self.parent.atomNamesBySeqId and self.parent.atomNamesBySeqId.has_key(residueKey):
            self.elements[-1].setAtomNames(self.parent.atomNamesBySeqId[residueKey])
          
        # Set mapping if required.
        if ccpMapping:
          self.elements[-1].code3Letter = ccpMapping[1][1]
          self.elements[-1].residueType = ccpMapping[1][0]
          
          if ccpMapping[2]:
            self.elements.append(NmrStarSequenceElement(sequenceTableTags,i,tableTagNames))
            self.elements[-1].code3Letter = ccpMapping[2][1]
            self.elements[-1].residueType = ccpMapping[2][0] # This is usually not there!
    
    """
    print "SEQUENCE"  
    for seqEl in self.elements:
      print seqEl.seqCode,
      if hasattr(seqEl, 'insertionCode'):
        print seqEl.insertionCode
      else:
        print
    """

    #
    # Get molecule database links.
    #
    
    tableName = '_Entity_db_link'
    tagNames = ('Database_code','Accession_code','Entry_mol_name','Entry_details')

    if self.saveFrame.tables.has_key(tableName):
      tableTags = self.saveFrame.tables[tableName].tags
  
      for i in range(0,len(tableTags['Database_code'])):
  
        infoList = []
        for tagName in tagNames:
          infoList.append(tableTags[tagName][i])           
            
        self.dbLinks.append(MoleculeDbLink(self,*infoList))

    #
    # Get molecule database links.
    #

    tableName = '_Entity_biological_function'
    tagNames = ('Biological_function',)

    if self.saveFrame.tables.has_key(tableName):
      tableTags = self.saveFrame.tables[tableName].tags
  
      for i in range(0,len(tableTags['Biological_function'])):
  
        infoList = []
        for tagName in tagNames:
          infoList.append(tableTags[tagName][i])           
            
        self.bioFunctions.append(MoleculeBioFunction(self,*infoList))

  def formatMoleculeName(self,molName):
    
    # Have to make sure that this is not too long, only 255 chars allowed in CCPN. Truncating to
    # 250 so some extra place for adding chars in sequence read pipeline
    
    return molName[:65]

class NmrStarSequenceElement:
  """
  Information for an element in the sequence (residue)
  """
  def __init__(self,sequenceTableTags,i,tableTagNames):
  
    self.bonds = {}
    self.formatCode = None
    
    self.atomNames = []
    
    self.secStrucInfo = {}
    
    for attrName in tableTagNames.keys():
      if sequenceTableTags.has_key(tableTagNames[attrName]) and len(sequenceTableTags[tableTagNames[attrName]]) > i and sequenceTableTags[tableTagNames[attrName]][i] is not None:
        setattr(self,attrName,sequenceTableTags[tableTagNames[attrName]][i])
      else:
        setattr(self,attrName,None)

  def setAtomNames(self,atomNames):
    
    self.atomNames = atomNames

  def setBond(self,bondType,atomName,bondedSeqEl,bondedAtomName):
  
    if not self.bonds.has_key(bondType):
      self.bonds[bondType] = {}
    
    if not self.bonds[bondType].has_key(atomName):
      self.bonds[bondType][atomName] = []
    
    if not (bondedSeqEl,bondedAtomName) in self.bonds[bondType][atomName]:
      self.bonds[bondType][atomName].append((bondedSeqEl,bondedAtomName))
      print "  Found %s bond from %s.%s - %s.%s" % (bondType,self.seqCode,atomName,bondedSeqEl.seqCode,bondedAtomName)

# TODO - move these to projectIO or generalIO - it complained earlier due to circular references I think?

class MoleculeDbLink:

  def __init__(self,parent,dbName,dbAcc,dbMolName,dbDetails):

    self.parent = parent
    self.dbName = dbName
    self.dbAcc = dbAcc
    self.dbMolName = dbMolName
    self.dbDetails = dbDetails

class MolSystemDbLink:

  def __init__(self,parent,dbName,dbAcc,exptMethod,structResolution):

    self.parent = parent
    self.dbName = dbName
    self.dbAcc = dbAcc
    self.exptMethod = exptMethod
    self.structResolution = structResolution

class MoleculeBioFunction:

  def __init__(self,parent,bioFunc):

    self.parent  = parent
    self.bioFunc = bioFunc

class MolSystemBioFunction:

  def __init__(self,parent,bioFunc):

    self.parent  = parent
    self.bioFunc = bioFunc

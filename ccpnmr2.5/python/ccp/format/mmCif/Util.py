import sys

from ccp.format.mmCif.sans.handlers import ErrorHandler, ContentHandler

from memops.universal.Util import returnInt, returnFloat

class DataDictionary_mmCIF( ContentHandler, ErrorHandler ):

  separator = "."
  
  def __init__(self):
  
    self.dataContent = {}
    
  def comment( self, line, text ) :
    #print "Comment:", text, "in line", line
    return False
  def startData( self, line, name ) :
    #print "Start data block", name, "in line", line
    return False
  def endData( self, line, name ) :
    #print "End data block", name, "in line", line
    pass
  def startSaveFrame( self, line, name ) :
    #print "Start saveframe", name, "in line", line
    return False
  def endSaveFrame( self, line, name ) :
    #print "End saveframe", name, "in line", line
    return False
  def startLoop( self, line ) :
    #print "Start loop in line", line
    return False
  def endLoop( self, line ) :
    #print "End loop in line", line
    return False
    
  def data( self, tag, tagline, val, valline, delim, inloop ):
    
    (mainName,tagName) = tag.split(self.separator)

    if not self.dataContent.has_key(mainName):
      self.dataContent[mainName] = {}

    if not self.dataContent[mainName].has_key(tagName):
      self.dataContent[mainName][tagName] = []
    
    if val in ('.','?'):
      val = None
    
    self.dataContent[mainName][tagName].append(val)

    return False

  def error( self, line, msg ) :
    print "mmCIF parse error in line", line, ":", msg
    return True
    
  def getPdbCode(self):
  
    pdbCode = None
  
    if self.dataContent.has_key('_database_2'):
    
      for i in range(len(self.dataContent['_database_2']['database_id'])):
        if self.dataContent['_database_2']['database_id'][i] == 'PDB':
          pdbCode = self.dataContent['_database_2']['database_code'][i]
          
    return pdbCode

  def getSequenceInfo(self):
  
    """
    
    Note that need to reconstruct *molecules* on the basis of cif chain codes and bonds between these molecule bits.
      
    """
    
    sequenceInfo = {}

    if self.dataContent.has_key('_entity'):
    
      for i in range(len(self.dataContent['_entity']['id'])):
        entityId = returnInt(self.dataContent['_entity']['id'][i])
        moleculeType = self.dataContent['_entity']['type'][i]
        moleculeName = self.dataContent['_entity']['pdbx_description'][i].strip()
        
        sequenceInfo[entityId] = [moleculeType,moleculeName, None, None, {}]

        if moleculeType == 'polymer':

          mainName = '_entity_poly'
          for j in range(len(self.dataContent[mainName]['entity_id'])):
            curEntityId = returnInt(self.dataContent[mainName]['entity_id'][j])
            if entityId == curEntityId:
              moleculePolymerType = self.dataContent[mainName]['type'][j]
              chainCodes = self.dataContent[mainName]['pdbx_strand_id'][j]
              
              sequenceInfo[entityId][2] = moleculePolymerType
              sequenceInfo[entityId][3] = chainCodes
              
          mainName = '_pdbx_poly_seq_scheme'
          for j in range(len(self.dataContent[mainName]['entity_id'])):
            curEntityId = returnInt(self.dataContent[mainName]['entity_id'][j])
            if entityId == curEntityId:
              cifChainCode = self.dataContent[mainName]['asym_id'][j]
              cifResLabel = self.dataContent[mainName]['mon_id'][j]
              cifSeqId = returnInt(self.dataContent[mainName]['seq_id'][j])
              pdbSeqCode = returnInt(self.dataContent[mainName]['pdb_seq_num'][j])
              authSeqCode = returnInt(self.dataContent[mainName]['auth_seq_num'][j])
              pdbResLabel = self.dataContent[mainName]['pdb_mon_id'][j]
              authResLabel = self.dataContent[mainName]['auth_mon_id'][j]
              pdbChainCode = self.dataContent[mainName]['pdb_strand_id'][j]
              pdbInsertionCode = self.dataContent[mainName]['pdb_ins_code'][j]
              isHetero = self.dataContent[mainName]['hetero'][j]
              
              if not sequenceInfo[entityId][-1].has_key(cifChainCode):
                sequenceInfo[entityId][-1][cifChainCode] = []

              sequenceInfo[entityId][-1][cifChainCode].append((cifSeqId,cifResLabel, pdbResLabel, pdbChainCode, pdbSeqCode, pdbInsertionCode, authResLabel, authSeqCode, isHetero))

          """
          mainName = '_entity_poly_seq'
          for j in range(len(self.dataContent[mainName]['entity_id'])):
            curEntityId = returnInt(self.dataContent[mainName]['entity_id'][j])
            if entityId == curEntityId:
              seqCode = returnInt(self.dataContent[mainName]['num'][j])
              resLabel = self.dataContent[mainName]['mon_id'][j]
              print  "     ", seqCode, resLabel
          """

        mainName = '_pdbx_nonpoly_scheme'
        if self.dataContent.has_key(mainName):
          for j in range(len(self.dataContent[mainName]['entity_id'])):
            curEntityId = returnInt(self.dataContent[mainName]['entity_id'][j])
            if entityId == curEntityId:
              cifChainCode = self.dataContent[mainName]['asym_id'][j]
              cifResLabel = self.dataContent[mainName]['mon_id'][j]
              pdbSeqCode = returnInt(self.dataContent[mainName]['pdb_seq_num'][j])
              authSeqCode = returnInt(self.dataContent[mainName]['auth_seq_num'][j])
              pdbResLabel = self.dataContent[mainName]['pdb_mon_id'][j]
              authResLabel = self.dataContent[mainName]['auth_mon_id'][j]
              pdbChainCode = self.dataContent[mainName]['pdb_strand_id'][j]
              pdbInsertionCode = self.dataContent[mainName]['pdb_ins_code'][j]
              
              if not sequenceInfo[entityId][-1].has_key(cifChainCode):
                sequenceInfo[entityId][-1][cifChainCode] = []

              sequenceInfo[entityId][-1][cifChainCode].append((pdbSeqCode,cifResLabel, pdbResLabel, pdbChainCode, pdbSeqCode, pdbInsertionCode, authResLabel, authSeqCode, True))
  
    return sequenceInfo

  def getBondInfo(self):

    bondInfo = {}
      
    # Bonds
    mainName = '_struct_conn'
    if self.dataContent.has_key(mainName):
      for j in range(len(self.dataContent[mainName]['id'])):
        bondType = self.dataContent[mainName]['conn_type_id'][j]

        if not bondInfo.has_key(bondType):
          bondInfo[bondType] = []

        bondInfo[bondType].append([])

        for k in range(1,3):
          cifChainCode = self.dataContent[mainName]['ptnr%d_label_asym_id' % k][j]
          cifResLabel = self.dataContent[mainName]['ptnr%d_label_comp_id' % k][j]
          cifSeqId = self.dataContent[mainName]['ptnr%d_label_seq_id' % k][j]
          if cifSeqId != None:
            cifSeqId = returnInt(cifSeqId)
          authSeqId = self.dataContent[mainName]['ptnr%d_auth_seq_id' % k][j]
          if authSeqId != None:
            authSeqId = returnInt(authSeqId)
          authChainCode = self.dataContent[mainName]['ptnr%d_auth_asym_id' % k][j]

          pdbInsertionCode = self.dataContent[mainName]['pdbx_ptnr%d_PDB_ins_code' % k][j]
          atomName = self.dataContent[mainName]['ptnr%d_label_atom_id' % k][j]

          bondInfo[bondType][-1].append((cifChainCode, cifResLabel, cifSeqId, authChainCode, authSeqId, atomName, pdbInsertionCode))

    return bondInfo
    
  def getChemCompInfo(self):
  
    chemCompInfo = {}
     
    # Bonds
    mainName = '_chem_comp'
    if self.dataContent.has_key(mainName):
      for j in range(len(self.dataContent[mainName]['id'])):
        cifResLabel = self.dataContent[mainName]['id'][j]
        cifResidueType = self.dataContent[mainName]['type'][j]
        nonStandardFlag = self.dataContent[mainName]['mon_nstd_flag'][j]
        name = self.dataContent[mainName]['name'][j]

        chemCompInfo[cifResLabel] = (cifResidueType, nonStandardFlag, name)

    return chemCompInfo
    
  def getCoordinateInfo(self):
  
    coordinateInfo = {}
  
    mainName = '_atom_site'
    if self.dataContent.has_key(mainName):
      for i in range(len(self.dataContent[mainName]['id'])):
      
        atomClass = self.dataContent[mainName]['group_PDB'][i]

        serial = returnInt(self.dataContent[mainName]['id'][i])
        elementType = self.dataContent[mainName]['type_symbol'][i]
        atomName = self.dataContent[mainName]['label_atom_id'][i]
        
        entityId = returnInt(self.dataContent[mainName]['label_entity_id'][i])

        cifResLabel = self.dataContent[mainName]['label_comp_id'][i]
        cifChainCode = self.dataContent[mainName]['label_asym_id'][i]
        cifSeqId = self.dataContent[mainName]['label_seq_id'][i]
        if cifSeqId != None:
          cifSeqId = returnInt(cifSeqId)

        pdbInsertionCode = self.dataContent[mainName]['pdbx_PDB_ins_code'][i]
        
        x = returnFloat(self.dataContent[mainName]['Cartn_x'][i])
        y = returnFloat(self.dataContent[mainName]['Cartn_y'][i])
        z = returnFloat(self.dataContent[mainName]['Cartn_z'][i])
        occupancy = returnFloat(self.dataContent[mainName]['occupancy'][i])
        Bfactor = returnFloat(self.dataContent[mainName]['B_iso_or_equiv'][i])
        
        authSeqCode = self.dataContent[mainName]['auth_seq_id'][i]
        if authSeqCode != None:
          authSeqCode = returnInt(authSeqCode)
          
        authResLabel = self.dataContent[mainName]['auth_comp_id'][i]
        authChainCode = self.dataContent[mainName]['auth_asym_id'][i]
        authAtomName = self.dataContent[mainName]['auth_atom_id'][i]
        
        # Currently not used.
        formalCharge = self.dataContent[mainName]['pdbx_formal_charge'][i]
        if formalCharge != None:
          formalCharge = returnInt(formalCharge)
        
        # Set the model
        model = returnInt(self.dataContent[mainName]['pdbx_PDB_model_num'][i])

        if not coordinateInfo.has_key(model):
          coordinateInfo[model] = []
          
        coordinateInfo[model].append((serial, elementType, atomName, cifResLabel, cifChainCode, cifSeqId, entityId, pdbInsertionCode, x, y, z, occupancy, Bfactor, authSeqCode, authResLabel, authChainCode, authAtomName, atomClass))

    return coordinateInfo
      
"""
List of mainNames:

['_chem_comp', '_pdbx_nmr_exptl_sample_conditions', '_database_PDB_rev_record', '_audit_author', '_diffrn', '_pdbx_version', '_atom_sites', '_atom_type', '_entity', '_entity_poly', '_struct_ref_seq', '_citation_author', '_pdbx_prerelease_seq', '_citation', '_pdbx_nmr_ensemble', '_database_PDB_matrix', '_pdbx_nmr_refine', '_pdbx_validate_torsion', '_exptl_crystal', '_exptl', '_struct_conn', '_pdbx_nmr_spectrometer', '_pdbx_nmr_software', '_struct_sheet_order', '_pdbx_poly_seq_scheme', '_database_PDB_rev', '_pdbx_struct_sheet_hbond', '_audit_conform', '_struct_keywords', '_entity_keywords', '_struct_asym', '_entity_poly_seq', '_pdbx_nmr_exptl', '_struct_biol', '_struct_sheet_range', '_pdbx_database_related', '_atom_site', '_pdbx_database_status', '_diffrn_radiation', '_entity_src_gen', '_pdbx_nmr_representative', '_struct', '_struct_conf_type', '_struct_sheet', '_struct_ref', '_pdbx_nmr_sample_details', '_database_2', '_struct_conn_type', '_diffrn_radiation_wavelength', '_entry']

"""

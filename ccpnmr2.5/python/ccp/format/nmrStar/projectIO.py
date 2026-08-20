"""
======================COPYRIGHT/LICENSE START==========================

projectIO.py: I/O for nmrStar project file

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

from ccp.format.nmrStar.sequenceIO import NmrStarSequenceFile
from ccp.format.nmrStar.chemShiftsIO import NmrStarChemShiftFile
from ccp.format.nmrStar.orderParamIO import NmrStarOrderParamFile
from ccp.format.nmrStar.hExchRateIO import NmrStarHExchRateFile
from ccp.format.nmrStar.hExchProtectionIO import NmrStarHExchProtectionFile
from ccp.format.nmrStar.t1RelaxIO import NmrStarT1RelaxFile
from ccp.format.nmrStar.t1RhoRelaxIO import NmrStarT1RhoRelaxFile
from ccp.format.nmrStar.t2RelaxIO import NmrStarT2RelaxFile
from ccp.format.nmrStar.hetNoeIO import NmrStarHetNoeFile
from ccp.format.nmrStar.jCouplingIO import NmrStarJCouplingFile
from ccp.format.nmrStar.rdcIO import NmrStarRdcFile
from ccp.format.nmrStar.distanceConstraintsIO import NmrStarDistanceConstraintFile
from ccp.format.nmrStar.dihedralConstraintsIO import NmrStarDihedralConstraintFile
from ccp.format.nmrStar.rdcConstraintsIO import NmrStarRdcConstraintFile
from ccp.format.nmrStar.coordinatesIO import NmrStarCoordinateFile
from ccp.format.nmrStar.peopleAndCitationsIO import NmrStarPeopleAndCitationsFile

from ccp.format.nmrStar.generalIO import NmrStarFile

from memops.universal.Util import returnFloat

#####################
# Class definitions #
#####################

class NmrStarProjectFile(NmrStarFile):

  # Information on file level

  def initialize(self, version = '2.1.1'):

    # These are fully FC supported (in DataFormat.py)
    self.distanceConstraintFiles = []
    self.hBondConstraintFiles = []
    self.dihedralConstraintFiles = []
    self.rdcConstraintFiles = []
    self.chemShiftFiles = []
    self.orderParamFiles = []
    self.hExchRateFiles = []
    self.hExchProtectionFiles = []
    self.t1RelaxFiles = []
    self.t1RhoRelaxFiles = []
    self.t2RelaxFiles = []
    self.hetNoeFiles = []
    self.jCouplingFiles = []
    self.rdcFiles = []
    self.sequenceFiles = []
    self.coordinateFiles = []
    self.peopleAndCitationsFiles = []
    
    # These are only available for 3.1 and on NMR-STAR project level reading
    self.entryInformation = None
    self.studyFiles = []
    self.sampleFiles = []
    self.sampleConditionFiles = []
    self.softwareFiles = []
    self.methodFiles = []
    self.confStatFiles = []
    self.instrumentFiles = []
    self.probeFiles = []
    self.nmrExperimentFiles = []
    self.chemShiftRefFiles = []
    self.sourceFiles = []

    if not self.version:
      self.version = version

    self.components = ['all']
    
    self.setComponents()

  def read(self,verbose = 0,**keywds):
  
    #
    # Read is now version 2.1.1 and 3.0
    #   
    
    if not self.checkVersion():
      return False
      
    self.readGeneric(verbose = verbose)

    if not self.readStatusCheck():
      return False

    MRFileBlockName = '_BMRB_dev_Position_block_in_MR_file'
    MRFileSaveFrame = 'MR_file_characteristics'
    cyanaLibUsed = 0
   
    #
    # Distance and H-bond constraints
    #  

    saveFrameNames = ['distance_constraints' ,'general_distance_constraints']

    for saveFrameName in saveFrameNames:
      if self.sfs.has_key(saveFrameName):
  
        for saveFrame in self.sfs[saveFrameName]:
  
          isHBond = 0
  
          #
          # Check whether it's an H bond constraint bit
          #
          
          if self.version == '2.1.1':
  
            if MRFileBlockName in saveFrame.tags.keys():
              if self.sfs.has_key(MRFileSaveFrame):
                loop = self.sfs[MRFileSaveFrame][0].tables[MRFileBlockName]
                loopIndex = loop.tags[MRFileBlockName].index(saveFrame.tags[MRFileBlockName])
  
                if loop.tags['_BMRB_dev_Subtype'][loopIndex] == "hydrogen bond":
                  isHBond = 1
  
          elif self.version == '3.0':
  
            if saveFrame.tags.has_key('Subtype') and saveFrame.tags['Subtype'] == 'hydrogen bond':
              isHBond = 1
  
          elif self.version == '3.1':
  
            if saveFrame.tags.has_key('Constraint_type') and saveFrame.tags['Constraint_type'] == 'hydrogen bond':
              isHBond = 1
  
          constraintFile = NmrStarDistanceConstraintFile(self.name,self,saveFrame)
  
          if not isHBond:
            self.distanceConstraintFiles.append(constraintFile)
          else:
            self.hBondConstraintFiles.append(constraintFile)

    #
    # Dihedral constraints
    #

    saveFrameName = 'torsion_angle_constraints' 

    if self.sfs.has_key(saveFrameName):

      for saveFrame in self.sfs[saveFrameName]:

        constraintFile = NmrStarDihedralConstraintFile(self.name,self,saveFrame)
        self.dihedralConstraintFiles.append(constraintFile)

        if constraintFile.cyanaLibUsed:
          cyanaLibUsed = 1


    #
    # RDC constraints
    #
    
    if self.version == '3.1':
      saveFrameName = 'RDC_constraints' 
    #else:
    #  saveFrameName = 'residual_dipolar_couplings' 

    if self.sfs.has_key(saveFrameName):

      for saveFrame in self.sfs[saveFrameName]:

        constraintFile = NmrStarRdcConstraintFile(self.name,self,saveFrame)
        self.rdcConstraintFiles.append(constraintFile)

        if constraintFile.cyanaLibUsed:
          cyanaLibUsed = 1

    if cyanaLibUsed:
      print "Using CYANA library - courtesy of Peter Guentert."

    #
    # Sequence
    #
    
    saveFrameName = None
    if self.version == '2.1.1':
      saveFrameName = 'monomeric_polymer' # Using this because assembly saveframe not always present...
    elif self.version[0] == '3':
      saveFrameName = 'assembly' 

    if self.sfs.has_key(saveFrameName):

      for saveFrame in self.sfs[saveFrameName]:
        self.sequenceFiles.append(NmrStarSequenceFile(self.name,self,saveFrame))

    #
    # Coordinates
    #
    
    saveFrameName = None
    if self.version == '2.1.1':
      saveFrameName = 'conformer_family_coord_set' # TODO does this exist? 
    elif self.version[0] == '3':
      saveFrameName = 'conformer_family_coord_set' 

    if self.sfs.has_key(saveFrameName):
      
      localKeywds = {}
      if keywds.has_key('maxNum'):
        localKeywds['maxModelNum'] = keywds['maxNum']

      for saveFrame in self.sfs[saveFrameName]:
        self.coordinateFiles.append(NmrStarCoordinateFile(self.name,self,saveFrame,**localKeywds))

    #
    # Chemical shifts
    #

    saveFrameName = 'assigned_chemical_shifts' 

    if self.sfs.has_key(saveFrameName):

      for saveFrame in self.sfs[saveFrameName]:

        chemShiftFile = NmrStarChemShiftFile(self.name,self,saveFrame)

        tagNameDict = {'Details':            'details',
                       'Chem_shift_1H_err':  'chemShiftErr1H',
                       'Chem_shift_13C_err': 'chemShiftErr13C',
                       'Chem_shift_15N_err': 'chemShiftErr15N',
                       'Chem_shift_31P_err': 'chemShiftErr31P',
                       'Chem_shift_2H_err':  'chemShiftErr2H',
                       'Chem_shift_19F_err': 'chemShiftErr19F'}

        for tagName in tagNameDict.keys():
          if saveFrame.tags.has_key(tagName):
            setattr(chemShiftFile, tagNameDict[tagName], saveFrame.tags[tagName])
          else:
            setattr(chemShiftFile, tagNameDict[tagName], None)

        self.chemShiftFiles.append(chemShiftFile)

    #
    # J couplings
    #

    saveFrameName = 'coupling_constants'
 
    if self.sfs.has_key(saveFrameName):

      for saveFrame in self.sfs[saveFrameName]:

        jCouplingFile = NmrStarJCouplingFile(self.name,self,saveFrame)
        self.jCouplingFiles.append(jCouplingFile)

    #
    # RDCs
    #
    
    saveFrameName = None
    if self.version == '2.1.1':
      saveFrameName = 'residual_dipolar_couplings'
    elif self.version[0] == '3':
      saveFrameName = 'RDCs' 

    if self.sfs.has_key(saveFrameName):

      for saveFrame in self.sfs[saveFrameName]:

        rdcFile = NmrStarRdcFile(self.name,self,saveFrame)
        self.rdcFiles.append(rdcFile)
    
    #
    # Order parameters
    #
    
    saveFrameName = None
    if self.version == '2.1.1':
      saveFrameName = 'S2_parameters'
    elif self.version[0] == '3':
      saveFrameName = 'order_parameters' 

    if self.sfs.has_key(saveFrameName):

      for saveFrame in self.sfs[saveFrameName]:

        orderParamFile = NmrStarOrderParamFile(self.name,self,saveFrame)
        self.orderParamFiles.append(orderParamFile)

    #
    # People and citations
    #
    
    if self.version[0] == '3':
      for saveFrameName in ('entry_information','citations'):
        if self.sfs.has_key(saveFrameName):
          if saveFrameName == 'entry_information':
            saveFrame = self.sfs[saveFrameName][0] # Only ever one...
            self.peopleAndCitationsFiles.append(NmrStarPeopleAndCitationsFile(self.name,self,saveFrame,'contactPersons'))
            self.peopleAndCitationsFiles.append(NmrStarPeopleAndCitationsFile(self.name,self,saveFrame,'authors'))
          else:
            for saveFrame in self.sfs[saveFrameName]:
              self.peopleAndCitationsFiles.append(NmrStarPeopleAndCitationsFile(self.name,self,saveFrame,None))

    #
    # Project level information, only in 3.1 (and higher, eventually)
    #
    
    if self.version == '3.1':
      
      # This is more or less in order of importance!
      
      # Additional entry level information
      saveFrameName = 'entry_information'
      if self.sfs.has_key(saveFrameName):
        # Strictly speaking this is not necessary, but it is for hacking the correct NMR-STAR version
        # from two entry_information saveframes (which is illegal) in the NRG/remediated wwPDB project 
        for saveFrame in self.sfs[saveFrameName]:
          self.entryInformation = NmrStarEntryInformation(self.name,self,saveFrame)

      # Studies
      saveFrameName = 'study_list'
      if self.sfs.has_key(saveFrameName):
        for saveFrame in self.sfs[saveFrameName]:
          self.studyFiles.append(NmrStarStudyFile(self.name,self,saveFrame))

      # Sample conditions
      saveFrameName = 'sample_conditions'
      if self.sfs.has_key(saveFrameName):
        for saveFrame in self.sfs[saveFrameName]:
          self.sampleConditionFiles.append(NmrStarSampleConditionsFile(self.name,self,saveFrame))
      
      # Chemical shift referencing
      saveFrameName = 'chem_shift_reference'
      if self.sfs.has_key(saveFrameName):
        for saveFrame in self.sfs[saveFrameName]:
          self.chemShiftRefFiles.append(NmrStarChemShiftRefFile(self.name,self,saveFrame))
      
      # Sample information
      saveFrameName = 'sample'
      if self.sfs.has_key(saveFrameName):
        for saveFrame in self.sfs[saveFrameName]:
          self.sampleFiles.append(NmrStarSampleFile(self.name,self,saveFrame))

      # NMR experiments
      saveFrameName = 'experiment_list'
      if self.sfs.has_key(saveFrameName):
        for saveFrame in self.sfs[saveFrameName]:
          self.nmrExperimentFiles.append(NmrStarExperimentFile(self.name,self,saveFrame))
      
      # Instruments
      saveFrameName = 'NMR_spectrometer_list'
      if self.sfs.has_key(saveFrameName):
        for saveFrame in self.sfs[saveFrameName]:
          self.instrumentFiles.append(NmrStarInstrumentFile(self.name,self,saveFrame))

      # Probes
      saveFrameName = 'NMR_spectrometer_probe'
      if self.sfs.has_key(saveFrameName):
        for saveFrame in self.sfs[saveFrameName]:
          self.probeFiles.append(NmrStarProbeFile(self.name,self,saveFrame))
            
      # Software
      saveFrameName = 'software'
      if self.sfs.has_key(saveFrameName):
        for saveFrame in self.sfs[saveFrameName]:
          self.softwareFiles.append(NmrStarSoftwareFile(self.name,self,saveFrame))

      # Methods
      saveFrameName = 'method'
      if self.sfs.has_key(saveFrameName):
        for saveFrame in self.sfs[saveFrameName]:
          self.methodFiles.append(NmrStarMethodFile(self.name,self,saveFrame))

      # Ensemble stats
      saveFrameName = 'conformer_statistics'
      if self.sfs.has_key(saveFrameName):
        for saveFrame in self.sfs[saveFrameName]:
          self.confStatFiles.append(NmrStarConfStatFile(self.name,self,saveFrame))

      # Natural/experimental source
      for saveFrameName in ('natural_source','experimental_source'):
        if self.sfs.has_key(saveFrameName):
          for saveFrame in self.sfs[saveFrameName]:
            self.sourceFiles.append(NmrStarSourceFile(self.name,self,saveFrame,saveFrameName))

      # H-Exchange Rates
      saveFrameName = 'H_exch_rates'

      if self.sfs.has_key(saveFrameName):

        for saveFrame in self.sfs[saveFrameName]:

          hExchRateFile = NmrStarHExchRateFile(self.name,self,saveFrame)
          self.hExchRateFiles.append(hExchRateFile)

      # H-Protection Values
      saveFrameName = 'H_exch_protection_factors'

      if self.sfs.has_key(saveFrameName):

        for saveFrame in self.sfs[saveFrameName]:

          hExchProtectionFile = NmrStarHExchProtectionFile(self.name,self,saveFrame)
          self.hExchProtectionFiles.append(hExchProtectionFile)

      # T1 relaxation
      saveFrameName = 'heteronucl_T1_relaxation'

      if self.sfs.has_key(saveFrameName):

        for saveFrame in self.sfs[saveFrameName]:

          t1RelaxFile = NmrStarT1RelaxFile(self.name,self,saveFrame)
          self.t1RelaxFiles.append(t1RelaxFile)

      # T1Rho relaxation
      saveFrameName = 'heteronucl_T1rho_relaxation'

      if self.sfs.has_key(saveFrameName):

        for saveFrame in self.sfs[saveFrameName]:

          t1RhoRelaxFile = NmrStarT1RhoRelaxFile(self.name,self,saveFrame)
          self.t1RhoRelaxFiles.append(t1RhoRelaxFile)

      # T2 relaxation
      saveFrameName = 'heteronucl_T2_relaxation'

      if self.sfs.has_key(saveFrameName):

        for saveFrame in self.sfs[saveFrameName]:

          t2RelaxFile = NmrStarT2RelaxFile(self.name,self,saveFrame)
          self.t2RelaxFiles.append(t2RelaxFile)

      # NOEs
      saveFrameName = 'heteronucl_NOEs'

      if self.sfs.has_key(saveFrameName):

        for saveFrame in self.sfs[saveFrameName]:

          hetNoeFile = NmrStarHetNoeFile(self.name,self,saveFrame)
          self.hetNoeFiles.append(hetNoeFile)

    return True

  def readStatusCheck(self):
  
    #
    # Dummy function that can be used elsewhere to customize the nmrStar project
    # reading
    #
    
    return True

  def write(self,verbose = 0):
     
    #
    # This is taken over by NmrStarExport code for 3.1 and up. Only useful for 2.1.1
    #
    
    self.writeGeneric(verbose = verbose)

#
# Classes specific for NMR-STAR project level reading, some might be split out to IO.py files if become
# part of generic DataFormat data component handling.
#
# Note that the way this is set up is 'overkill', but necessary in case we need to fit this in into
# the DataFormat.py structure.
#

from ccp.format.nmrStar.generalIO import NmrStarGenericFile

class NmrStarProjectDataComponent(NmrStarGenericFile):

  """
  Information on file level.
  """
  
  def genericInit(self,parent,saveFrame):
    
    self.saveFrame = saveFrame
    
    self.parent = parent
    self.version = parent.version

    if self.saveFrame:
      self.parseSaveFrame()

  def checkVersion(self):
    
    if self.version not in ('3.1',):
      print "  Cannot parse nmrStar distance constraints for version %s" % self.version
      return False

    return True        

#
# Entry level information
#

class NmrStarEntryInformation(NmrStarProjectDataComponent):

  def initialize(self,parent,saveFrame):

    self.entry = None
    self.entryId = None
    self.entrySrcs = []
    self.relatedEntries = []
    self.structKeywds = []
    self.releaseInfo = []
    self.strucGenomInfo = []
    
    self.genericInit(parent,saveFrame)

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    self.entryId = self.saveFrame.tags['ID']

    infoList = []

    tagNames = ['Title', 'Details', 'Experimental_method', 'Experimental_method_subtype',
                'Special_processing_instructions', 'Version_type', 'Submission_date',
                'Accession_date', 'Origination', 'NMR_STAR_version', 'Original_NMR_STAR_version',
                'Dep_release_code_coordinates', 'Dep_release_code_nmr_constraints',
                'Dep_release_code_nmr_exptl', 'Dep_release_code_sequence',
                'CASP_target', 'Update_BMRB_accession_code', 'Replace_BMRB_accession_code',
                'Update_PDB_accession_code', 'Replace_PDB_accession_code',
                'BMRB_update_details', 'PDB_update_details', 'Release_request',
                'PDB_deposit_site', 'PDB_process_site', 'BMRB_deposit_site',
                'BMRB_process_site', 'RCSB_annotator', 'Assigned_BMRB_ID',
                'Assigned_BMRB_deposition_code', 'Assigned_PDB_ID', 'Assigned_PDB_deposition_code',
                'Date_nmr_constraints', 'Recvd_author_approval', 'PDB_date_submitted',
                'Author_release_status_code', 'Author_approval_type']
    
    for tagName in tagNames:
      # Possible that tag is not there (if it was empty)... take this into account
      if self.saveFrame.tags.has_key(tagName):
        value = self.saveFrame.tags[tagName]
      else:
        value = None
        
      infoList.append(value)
    
    #
    # Note: this is a NASTY HACK to get correct NMR_STAR_VERSION out of 'joint' coordinate/restraint file in NRG
    # and remediated wwPDB project. Should not be relevant for any other file.
    #
  
    if hasattr(self.parent,'entryInformation') and self.parent.entryInformation:
      if self.parent.entryInformation.entry.nmrStarVers:
        infoList[tagNames.index('NMR_STAR_version')] = self.parent.entryInformation.entry.nmrStarVers
    
    #
    # END nasty hack
    #

    self.entry = Entry(self,*infoList)

    self.tableName1 = '_Entry_src'
    tagNames = ('Organization_full_name',)

    if self.saveFrame.tables.has_key(self.tableName1):
      tableTags = self.saveFrame.tables[self.tableName1].tags
    
      for i in range(0,len(tableTags['Organization_full_name'])):
    
        infoList = []
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])
        
        self.entrySrcs.append(EntrySrc(self,*infoList))

    self.tableName2 = '_Related_entries'
    tagNames = ('Database_name', 'Database_accession_code', 'Relationship')

    if self.saveFrame.tables.has_key(self.tableName2):
      tableTags = self.saveFrame.tables[self.tableName2].tags
    
      for i in range(0,len(tableTags['Database_name'])):
    
        infoList = []
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])
        
        self.relatedEntries.append(RelatedEntry(self,*infoList))

    self.tableName3 = '_Struct_keywords'
    tagNames = ('Entry_ID','Keywords','Text')

    if self.saveFrame.tables.has_key(self.tableName3):
      tableTags = self.saveFrame.tables[self.tableName3].tags
    
      for i in range(0,len(tableTags['Entry_ID'])):
    
        infoList = []
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])

        self.structKeywds.append(StructKeyWord(self,*infoList))

    self.tableName4 = '_Release'
    tagNames = ('Entry_ID','Release_number','Date','Submission_date','Type','Author','Detail')

    if self.saveFrame.tables.has_key(self.tableName4):
      tableTags = self.saveFrame.tables[self.tableName4].tags
    
      for i in range(0,len(tableTags['Entry_ID'])):
    
        infoList = []
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])

        self.releaseInfo.append(Release(self,*infoList))

    self.tableName5 = '_SG_project'
    tagNames = ('Project_name','Full_name_of_center')

    if self.saveFrame.tables.has_key(self.tableName5):
      tableTags = self.saveFrame.tables[self.tableName5].tags
    
      for i in range(0,len(tableTags['Project_name'])):
    
        infoList = []
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])

        self.strucGenomInfo.append(StrucGenomics(self,*infoList))

class Entry:

  def __init__(self,parent,title,details,exptMethod,exptMethodSub,specProcInstr,versionType,subDate,accDate,
               origin,nmrStarVers,origNmrStarVers,depRelCoord,depRelConstr,depRelNmr,depRelSeq,
               caspTarget,updBmrbCode,replBmrbCode,updPdbCode,replPdbCode,bmrbUpdDet,pdbUpdDet,relReq,
               pdbDepSite,pdbProcSite,bmrbDepSite,bmrbProcSite,rcsbAnn,assBmrbId,assBmrbCode,assPdbId,assPdbCode,
               dateNmrConstr, recAuthAppr, pdbDateRec, authRelStatCode, authApprType):

    self.parent = parent

    self.title = title
    self.details = details

    self.exptMethod = exptMethod
    self.exptMethodSub = exptMethodSub
    self.specProcInstr = specProcInstr

    self.versionType = versionType
    self.subDate = subDate
    self.accDate = accDate
    self.origin = origin
    self.nmrStarVers = nmrStarVers
    self.origNmrStarVers = origNmrStarVers
    self.depRelCoord = depRelCoord
    self.depRelConstr = depRelConstr
    self.depRelNmr = depRelNmr
    self.depRelSeq = depRelSeq
    self.caspTarget = caspTarget
    self.updBmrbCode = updBmrbCode
    self.replBmrbCode = replBmrbCode
    self.updPdbCode = updPdbCode
    self.replPdbCode = replPdbCode
    self.bmrbUpdDet = bmrbUpdDet
    self.pdbUpdDet = pdbUpdDet
    self.relReq = relReq
    self.pdbDepSite = pdbDepSite
    self.pdbProcSite = pdbProcSite
    self.bmrbDepSite = bmrbDepSite
    self.bmrbProcSite = bmrbProcSite
    self.rcsbAnn = rcsbAnn
    self.assBmrbId = assBmrbId
    self.assBmrbCode = assBmrbCode
    self.assPdbId = assPdbId
    self.assPdbCode = assPdbCode
    self.dateNmrConstr = dateNmrConstr
    self.recAuthAppr = recAuthAppr
    self.pdbDateRec = pdbDateRec
    self.authRelStatCode = authRelStatCode
    self.authApprType = authApprType

class EntrySrc:

  def __init__(self,parent,orgFullName):

    self.parent = parent
    self.orgFullName = orgFullName

class RelatedEntry:

  def __init__(self,parent,dbName,dbAcc,reln):

    self.parent = parent
    self.dbName = dbName
    self.dbAcc = dbAcc
    self.reln = reln

class StructKeyWord:

  def __init__(self,parent,entryId,keywd,text):

    self.parent = parent
    self.entryId = entryId
    self.keywd = keywd
    self.text = text

class Release:

  def __init__(self,parent,entryId,relNo,date,subDateRel,relType,author,detail):

    self.parent = parent
    self.entryId = entryId
    self.relNo = relNo
    self.date = date
    self.subDateRel = subDateRel
    self.relType = relType
    self.author = author
    self.detail = detail

class StrucGenomics:

  def __init__(self,parent,projName,centerFullName):

    self.parent = parent
    self.projName = projName
    self.centerFullName = centerFullName

#
# Studies information
#

class NmrStarStudyFile(NmrStarProjectDataComponent):

  def initialize(self,parent,saveFrame):
    
    self.studyLists = []
    self.studies = []
    self.studyKeywds = []
    self.studyId = None
    
    self.genericInit(parent,saveFrame)

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    self.studyId = self.saveFrame.tags['ID']

    infoList = []

    infoList.append(self.saveFrame.title) # Safer to use title here - sf_framecode may be missing (Wim)

    self.studyLists.append(StudyList(self,*infoList))

    self.tableName1 = '_Study'
    tagNames = ('ID','Name','Type','Details')

    if self.saveFrame.tables.has_key(self.tableName1):
      tableTags = self.saveFrame.tables[self.tableName1].tags
    
      for i in range(0,len(tableTags['ID'])):
    
        infoList = []
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])
        
        self.studies.append(Study(self,*infoList))

    self.tableName2 = '_Study_keyword'
    tagNames = ('Study_ID','Keyword')

    if self.saveFrame.tables.has_key(self.tableName2):
      tableTags = self.saveFrame.tables[self.tableName2].tags
    
      for i in range(0,len(tableTags['Study_ID'])):
    
        infoList = []
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])           

        self.studyKeywds.append(StudyKeyWord(self,*infoList))


class StudyList:
  
  def __init__(self,parent,name):
    
    # Warning: need to do some retyping here, data type not always correct in original NMR-STAR
    
    self.parent = parent
    self.name = name

class Study:

  def __init__(self,parent,Id,name,type,details):

    self.parent = parent
    self.Id = Id
    self.name = name
    self.type = type
    self.details = details

class StudyKeyWord:

  def __init__(self,parent,studyId,keywd):

    self.parent = parent
    self.studyId = studyId
    self.keywd = keywd

#
# Sample conditions information
#

class NmrStarSampleConditionsFile(NmrStarProjectDataComponent):

  def initialize(self,parent,saveFrame):
    
    # These two are specific for the type of information we're reading
    self.sampleConditions = []
    self.sampleConditionSetId      = None
    self.sampleConditionSetName    = None
    self.sampleConditionSetDetails = None

    self.genericInit(parent,saveFrame)
    self.sampleConditionSetId = None
    
    self.genericInit(parent,saveFrame) # This always has to be here

  def parseSaveFrame(self):

    if not self.checkVersion():  # This always has to be here
      return
    
    self.sampleConditionSetId   = self.saveFrame.tags['ID']
    self.sampleConditionSetName = self.saveFrame.title

    if self.saveFrame.tags.has_key('Details'):
      self.sampleConditionSetDetails = self.saveFrame.tags['Details']

    #if self.sampleConditionSetId in self.sampleConditions:
    #  return    
        
    self.tableName = '_Sample_condition_variable'      # This has to be the table 'prefix' before the . as it appears in the NMR-STAR file
    tagNames = ('Type','Val','Val_err','Val_units')
    if self.saveFrame.tables.has_key(self.tableName):
      tableTags = self.saveFrame.tables[self.tableName].tags
    
      for i in range(len(tableTags['Val'])):    # Just pick any obligatory tag from the table to loop over the rows
    
        infoList = []
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])           
        
        self.sampleConditions.append(SampleCondition(self,*infoList)) # So self is the 'parent' of the SampleCondition data class

# TODO: this might be better kept in ccp.general.formatIO?
class SampleCondition:
  
  def __init__(self,parent,type,value,error,units):
    
    # Warning: need to do some retyping here, data type not always correct in original NMR-STAR

    self.parent = parent
    self.type = type

    self.value = None

    # Exceptions
    if type == 'pressure' and value == 'ambient':
      self.value = 1.0
    elif value:
      self.value = returnFloat(value)

    self.error = None

    if error:
      self.error = returnFloat(error)

    self.units = units
  
#
# Chemical shift referencing information
#

class NmrStarChemShiftRefFile(NmrStarProjectDataComponent):

  def initialize(self,parent,saveFrame):
    
    self.chemShiftReferences = []
    self.chemShiftRefId = None

    self.genericInit(parent,saveFrame)

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    self.chemShiftRefId = self.saveFrame.tags['ID']

    infoList = []

    tagNames = ('Proton_shifts_flag','Carbon_shifts_flag','Nitrogen_shifts_flag','Phosphorus_shifts_flag','Other_shifts_flag')
    
    for tagName in tagNames:
      if self.saveFrame.tags.has_key(tagName):
        value = self.saveFrame.tags[tagName]
      else:
        value = None
        
      infoList.append(value)

    self.chemShiftRef = ChemShiftRefList(self,*infoList)

    self.tableName1 = '_Chem_shift_ref'
    tagNames = ('Atom_type','Atom_isotope_number','Mol_common_name','Atom_group','Chem_shift_units','Chem_shift_val','Ref_method','Ref_type','Indirect_shift_ratio','External_ref_loc','External_ref_sample_geometry','External_ref_axis')

    if self.saveFrame.tables.has_key(self.tableName1):
      tableTags = self.saveFrame.tables[self.tableName1].tags
    
      for i in range(0,len(tableTags['Atom_type'])):
    
        infoList = []
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])
        
        self.chemShiftReferences.append(ChemShiftRef(self,*infoList))

class ChemShiftRefList:

  def __init__(self,parent,protonFlag,carbonFlag,nitrogenFlag,phosphorusFlag,otherFlag):

    self.parent = parent
    self.protonFlag = protonFlag
    self.carbonFlag = carbonFlag
    self.nitrogenFlag = nitrogenFlag
    self.phosphorusFlag = phosphorusFlag
    self.otherFlag = otherFlag

class ChemShiftRef:

  def __init__(self,parent,atomType,isotopeNum,molName,atomGroup,chemShiftUnits,chemShiftVal,refMeth,refType,indirectShiftRatio,location = None,geometry = None, axis = None):

    self.parent = parent
    self.atomType = atomType
    self.isotopeNum = isotopeNum
    self.molName = molName
    self.atomGroup = atomGroup
    self.chemShiftUnits = chemShiftUnits
    self.chemShiftVal = chemShiftVal
    self.refMeth = refMeth
    self.refType = refType
    self.indirectShiftRatio = indirectShiftRatio
    self.location = location
    self.geometry = geometry
    self.axis = axis

#
# Sample information
#

class NmrStarSampleFile(NmrStarProjectDataComponent):

  def initialize(self,parent,saveFrame):
    
    self.samples = []
    self.sampleComponents = []
    self.sampleId = None
    self.sampleDetails = None

    self.genericInit(parent,saveFrame)

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    self.sampleId = self.saveFrame.tags['ID']

    infoList = []

    infoList.append(self.saveFrame.title)

    tagNames = ('Type','Details','Solvent_system','Aggregate_sample_number')
    
    for tagName in tagNames:
      if self.saveFrame.tags.has_key(tagName):
        value = self.saveFrame.tags[tagName]
      else:
        value = None
        
      infoList.append(value)

    self.samples.append(Sample(self,*infoList))

    self.tableName = '_Sample_component'      # This has to be the table 'prefix' before the . as it appears in the NMR-STAR file
    tagNames = ('ID','Mol_common_name','Isotopic_labeling','Assembly_ID','Entity_ID','Concentration_val','Concentration_val_min','Concentration_val_max','Concentration_val_units','Concentration_val_err')
    if self.saveFrame.tables.has_key(self.tableName):
      tableTags = self.saveFrame.tables[self.tableName].tags
    
      for i in range(0,len(tableTags['Mol_common_name'])):    # Just pick any obligatory tag from the table to loop over the rows
    
        infoList = []
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])           

        self.sampleComponents.append(SampleComponent(self,*infoList)) # So self is the 'parent' of the SampleCondition data class


class Sample:
  
  def __init__(self,parent,name,sampleState,details,solventSys,aggrSamNo):
    
    # Warning: need to do some retyping here, data type not always correct in original NMR-STAR
    
    self.parent = parent
    self.name = name

    self.sampleState = sampleState
    self.details = details
    self.solventSys = solventSys
    self.aggrSamNo = aggrSamNo

class SampleComponent:

  def __init__(self,parent,Id,name,label,molSysId,molId,conc,conc_min,conc_max,conc_unit,conc_err):

    self.parent = parent
    self.Id = Id
    self.name = name
    self.label = label
    self.molSysId = molSysId
    self.molId = molId

    self.conc = conc
    self.conc_min = conc_min
    self.conc_max = conc_max
    self.conc_unit = conc_unit
    self.conc_err = conc_err


#
# NMR experiment information
#

class NmrStarExperimentFile(NmrStarProjectDataComponent):

  def initialize(self,parent,saveFrame):
    
    self.experiments = []
    self.experimentListId = None
    
    self.genericInit(parent,saveFrame)

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    self.experimentListId = self.saveFrame.tags['ID']

    tagNameDict = {'Details': 'experimentListDetails'}

    for tagName in tagNameDict.keys():
      if self.saveFrame.tags.has_key(tagName):
        setattr(self, tagNameDict[tagName], self.saveFrame.tags[tagName])
      else:
        setattr(self, tagNameDict[tagName], None)

    self.tableName = '_Experiment'

    tagNames = ('ID','Name','Raw_data_flag','Sample_ID','Sample_label','Sample_state','Sample_condition_list_ID','Sample_condition_list_label','NMR_spectrometer_ID','NMR_spectrometer_label','NMR_spectrometer_probe_ID','NMR_spectrometer_probe_label')
    if self.saveFrame.tables.has_key(self.tableName):
      tableTags = self.saveFrame.tables[self.tableName].tags
    
      for i in range(0,len(tableTags['ID'])):
    
        infoList = []
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])           

        self.experiments.append(Experiment(self,*infoList))

class Experiment:

  def __init__(self,parent,Id,name,raw_data_flag,sampleId,sampleLabel,sampleState,sampleCondListId,sampleCondListLabel,specId,specLabel,probeId,probeLabel):

    self.parent = parent
    self.Id = Id
    self.name = name
    self.raw_data_flag = raw_data_flag
    self.sampleId = sampleId
    self.sampleLabel = sampleLabel
    self.sampleState = sampleState
    self.sampleCondListId = sampleCondListId
    self.sampleCondListLabel = sampleCondListLabel
    self.specId = specId
    self.specLabel = specLabel
    self.probeId = probeId
    self.probeLabel = probeLabel

#
# NMR spectrometer information
#

class NmrStarInstrumentFile(NmrStarProjectDataComponent):

  def initialize(self,parent,saveFrame):
    
    self.instruments = []
    self.specId = None
    self.specListDetails = None

    self.genericInit(parent,saveFrame)

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    self.specId = self.saveFrame.tags['ID']

    if self.saveFrame.tags.has_key('Details'):
      self.specListDetails = self.saveFrame.tags['Details']

    self.tableName = '_NMR_spectrometer_view'

    tagNames = ('ID','Name','Manufacturer','Model','Field_strength','Serial_number','Details')
    if self.saveFrame.tables.has_key(self.tableName):
      tableTags = self.saveFrame.tables[self.tableName].tags
    
      for i in range(0,len(tableTags['ID'])):
    
        infoList = []
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])

        self.instruments.append(NmrSpectrometer(self,*infoList))

class NmrSpectrometer:

  def __init__(self,parent,Id,name,manufacturer,model,field_strength,serialNumber,details):

    self.parent = parent
    self.Id = Id
    self.name = name
    self.manufacturer = manufacturer
    self.model = model
    self.field_strength = field_strength
    self.serialNumber = serialNumber
    self.details = details

#
# NMR Probe information
#

class NmrStarProbeFile(NmrStarProjectDataComponent):

  def initialize(self,parent,saveFrame):

    self.probes = []
    self.probeId = None

    self.genericInit(parent,saveFrame)

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    self.probeId = self.saveFrame.tags['ID']

    infoList = []

    tagNames = ('ID','Manufacturer','Model','Serial_number','Details','Diameter')

    for tagName in tagNames:
      if not tagName in self.saveFrame.tags:
        self.saveFrame.tags[tagName] = None
      infoList.append(self.saveFrame.tags[tagName])

    self.tableName = '_NMR_probe'

    if self.saveFrame.tables.has_key(self.tableName):
      tableTags = self.saveFrame.tables[self.tableName].tags

      if 'Type' in tableTags:
        infoList.append(tableTags['Type'][0])
      else:
        infoList.append(None)

    else:
      infoList.append(None)

    self.probes.append(NmrProbe(self,*infoList))

class NmrProbe:

  def __init__(self,parent,Id,manufacturer,model,serialNumber,details,diameter,pType):

    self.parent = parent
    self.Id = Id
    self.manufacturer = manufacturer
    self.model = model
    self.serialNumber = serialNumber
    self.details = details
    self.diameter = diameter
    self.pType = pType

#
# Software information
#

class NmrStarSoftwareFile(NmrStarProjectDataComponent):

  def initialize(self,parent,saveFrame):
    
    self.software = []
    self.softwareVendors = []
    self.softwareId = None
    
    self.genericInit(parent,saveFrame)

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    self.softwareId = self.saveFrame.tags['ID']

    infoList = []

    tagNames = ('Name','Version','Details')

    for tagName in tagNames:
      if not tagName in self.saveFrame.tags:
        self.saveFrame.tags[tagName] = None
      infoList.append(self.saveFrame.tags[tagName])

    self.tableName = '_Vendor'

    tagNames = ('Name','Address','Electronic_address')
    if self.saveFrame.tables.has_key(self.tableName):
      tableTags = self.saveFrame.tables[self.tableName].tags
    
      for tagName in tagNames:
        if tagName not in tableTags:
          infoList.append(None)
        else:
          infoList.append(tableTags[tagName][0])

    else:
      infoList.append(None)
      infoList.append(None)
      infoList.append(None)

    self.tableName = '_Task'

    taskString = None

    tagNames = ('Task',)
    if self.saveFrame.tables.has_key(self.tableName):
      tableTags = self.saveFrame.tables[self.tableName].tags

      taskList = []
      for i in range(0,len(tableTags['Task'])):
    
        for tagName in tagNames:
          if tagName not in tableTags:
            taskList.append(None)
          else:
            taskList.append(tableTags[tagName][i])

        taskString = ', '.join(taskList)

    infoList.append(taskString)

    self.software.append(Software(self,*infoList))

class Software:

  def __init__(self,parent,name,version,details,venName,venAddress,venEAddress,tasks):

    self.parent = parent
    self.name = name
    self.version = version
    self.details = details

    self.venName = venName
    self.venAddress = venAddress
    self.venEAddress = venEAddress
    self.tasks = tasks

#
# Method information
#

class NmrStarMethodFile(NmrStarProjectDataComponent):

  def initialize(self,parent,saveFrame):
    
    self.method = []
    self.methodId = None
    
    self.genericInit(parent,saveFrame)

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    self.methodId = self.saveFrame.tags['ID']

    infoList = []

    tagNames = ('Details',)

    for tagName in tagNames:
      if not tagName in self.saveFrame.tags:
        self.saveFrame.tags[tagName] = None
      infoList.append(self.saveFrame.tags[tagName])

    self.method.append(Method(self,*infoList))

class Method:

  def __init__(self,parent,details):

    self.parent = parent
    self.details = details

#
# Ensemble stats information
#

class NmrStarConfStatFile(NmrStarProjectDataComponent):

  def initialize(self,parent,saveFrame):
    
    self.confStat = []
    self.confStatId = None
    
    self.genericInit(parent,saveFrame)

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    self.confStatId = self.saveFrame.tags['ID']

    infoList = []

    tagNames = ('Conf_family_coord_set_ID,','Details','Conformer_calculated_total_num','Conformer_submitted_total_num',
                'Conformer_selection_criteria','Representative_conformer','Rep_conformer_selection_criteria')

    for tagName in tagNames:
      if not tagName in self.saveFrame.tags:
        self.saveFrame.tags[tagName] = None
      infoList.append(self.saveFrame.tags[tagName])

    self.confStat.append(ConfStat(self,*infoList))

class ConfStat:

  def __init__(self,parent,ensId,details,calculated,submitted,criteria,representative,repr_criteria):

    self.parent = parent
    self.ensId = ensId
    self.details = details
    self.calculated = calculated
    self.submitted = submitted
    self.criteria = criteria
    self.representative = representative
    self.repr_criteria = repr_criteria

#
# Natural/experimental information
#

class NmrStarSourceFile(NmrStarProjectDataComponent):

  def initialize(self,parent,saveFrame,saveFrameName):

    self.saveFrameName = saveFrameName
    
    if self.saveFrameName == 'natural_source':
      self.source = 'natural'
    elif self.saveFrameName == 'experimental_source':
      self.source = 'experimental'
    else:
      self.source = None
      
    self.sources = []
    self.sourceId = None
    
    self.genericInit(parent,saveFrame)

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    self.sourceId = self.saveFrame.tags['ID']

    if self.source == 'natural':
      self.tableName = '_Entity_natural_src'
      tagNames = ('ID','Entity_ID','Entity_chimera_segment_ID','Type',None,'NCBI_taxonomy_ID','Organism_name_scientific','Organism_name_common','Strain','Variant','Organ','Tissue','Cell_line','Cell_type','ATCC_number',None,'Plasmid','Details','Species','Genus','Superkingdom','Kingdom',None,'Gene_mnemonic','Organelle')

    elif self.source == 'experimental':
      self.tableName = '_Entity_experimental_src'
      tagNames = ('ID','Entity_ID','Entity_chimera_segment_ID',None,'Production_method','Host_org_NCBI_taxonomy_ID','Host_org_scientific_name','Host_org_name_common','Host_org_strain','Host_org_variant','Host_org_organ','Host_org_tissue','Host_org_cell_line','Host_org_cell_type','Host_org_ATCC_number','Vector_type','Vector_name','Details','Host_org_species','Host_org_genus',None,None,'PDBview_host_org_vector_name',None,'Host_org_organelle')

    if self.saveFrame.tables.has_key(self.tableName):
      tableTags = self.saveFrame.tables[self.tableName].tags

      for i in range(0,len(tableTags['ID'])):
    
        infoList = [self.source]
        for tagName in tagNames:
          if tagName not in tableTags:
            infoList.append(None)
          else:
            infoList.append(tableTags[tagName][i])

        self.sources.append(Source(self,*infoList))

class Source:

  def __init__(self,parent,sourceType,Id,entId,entChimId,orgType,prodMethod,ncbiTaxId,orgSci,orgName,strain,variant,organ,tissue,cell_line,cell_type,atccNo,vecType,plasmid,details,species,genus,superKingdom,kingdom,pdbVectorName,geneMnemonic,organelle):

    self.parent = parent
    self.sourceType = sourceType
    self.Id = Id
    self.entId = entId
    self.entChimId = entChimId
    self.orgType = orgType
    self.prodMethod = prodMethod
    self.ncbiTaxId = ncbiTaxId
    self.orgSci = orgSci
    self.orgName = orgName
    self.strain = strain
    self.variant = variant
    self.organ = organ
    self.tissue = tissue
    self.cell_line = cell_line
    self.cell_type = cell_type
    self.atccNo = atccNo
    self.vecType = vecType
    self.plasmid = plasmid
    self.details = details
    self.species = species
    self.genus = genus
    self.superKingdom = superKingdom
    self.kingdom = kingdom
    self.pdbVectorName = pdbVectorName
    self.geneMnemonic = geneMnemonic
    self.organelle = organelle

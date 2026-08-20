"""
======================COPYRIGHT/LICENSE START==========================

generalIO.py: General I/O information and code for nmrStar files.

Copyright (C) 2005-2010 Wim Vranken (European Bioinformatics Institute)

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

import os, sys, string

from memops.universal.Io import getTopDirectory

from ccp.format.general.formatIO import FormatFile

from ccp.format.nmrStar.util import returnStarInt
from ccp.format.nmrStar.util import returnStarFloat
from ccp.format.nmrStar.util import returnStarString
from ccp.format.nmrStar.util import getNmrStarValue

# Import Jurgen Doreleijers' Star reader/writer stuff
import ccp.format.nmrStar.bmrb.File as nmrStar

from ccp.format.general.Constants import defaultMolCode

"""

Notes for rewriting this script for more generic nmrStar treatment:


IF ANY OF BELOW INCLUDES ACTUAL CCPNMR STUFF THEN DO NOT PUT IT HERE!!

- separate reference data classes from 'real' ones, put link in 'real' ones to
  reference data
  
- make loops to set reoccuring tags (e.g. atom_name_1, ...)

- for reference data make table and tag classes!!! Use these (also keep ..Names list for
  quick string checking)

- add 'attribute' names to sfDict (or tag class)

- Set attribute for ref SF class: isCcpMapped (if false, store as application data linked to Entry!)

- Set attribute for ref SF tables: isCcpMapped (if false, store as application data to
  whatever parent object is given during interpretation. Could in principle set this as well...)

- Set attribute for ref SF tags? Is more difficult to set... also give parent object as for tables?
  could then check for appldata on that level... .



- Kinds of tags:

   - not mapped, use if exists as application data
   
   - mapped with direct mapping
   - mapped with function mapping
   - mapped with no direct or function mapping (comes from function at other tag, e.g. chainlabel from atom)
   - each/some of the 'mapped' could also exist as application data? Need flag to know what to use...

- example: Label_atom_ID has function 'getAtomLabels' which ALSO sets the
          Label_entity_assembly_ID, Label_entity_ID, Label_comp_index_ID and Label_comp_ID
          (these tags then have no function or direct mapping but are labelled as mapped)        



Could do setup of saveframes, tables, tags as reference classes (also keep data in there?),
then have functions linked to sets of tags to handle multiple ones (e.g. entity_assembly, entity, ...)

Organize functions for input as datamodel classes, output as saveframe/... ?


"""

#####################
# Class definitions #
#####################

class NmrStarGenericFile(FormatFile):

  def setGeneric(self):
  
    self.format = 'nmrStar'
    self.defaultMolCode = defaultMolCode
    self.tagSep = '.'

    self.measureExpts = []
    self.measureSofts = []

    self.structSofts = []
    self.structMeths = []

    if not self.version:
      self.version = '2.1.1' # is default
      
    # This is a hack so can modify an input file before reading. Useful for chemical shift files
    # where often only the loop_ is reported by programs (e.g. CARA).
    self.text = ""
  
  def checkVersion(self):

    if self.version not in ['2.1.1','3.0','3.1']:

      print "  Cannot parse nmrStar distance constraints for version %s" % self.version
      return False
      
    return True

  def setMeasureExperiments(self,tableName):
    
    if self.saveFrame.tables.has_key(tableName):
      exptTableTags = self.saveFrame.tables[tableName].tags
      numExpts = len(exptTableTags['Experiment_ID'])

      for i in range(0,numExpts):

        self.measureExpts.append(NmrStarMeasureExpt(self))
        self.measureExpts[-1].setData(exptTableTags,i)

  def setMeasureSoftwares(self,tableName):

    if self.saveFrame.tables.has_key(tableName):
      softwareTableTags = self.saveFrame.tables[tableName].tags
      numSoft = len(softwareTableTags['Software_ID'])

      for i in range(0,numSoft):

        self.measureSofts.append(NmrStarMeasureSoft(self))
        self.measureSofts[-1].setData(softwareTableTags,i)

  def setStructSoftwares(self,tableName):

    if self.saveFrame.tables.has_key(tableName):
      softwareTableTags = self.saveFrame.tables[tableName].tags
      numSoft = len(softwareTableTags['Software_ID'])

      for i in range(0,numSoft):

        self.structSofts.append(NmrStarStructSoft(self))
        self.structSofts[-1].setData(softwareTableTags,i)

  def setStructMethods(self,tableName):

    if self.saveFrame.tables.has_key(tableName):
      methodTableTags = self.saveFrame.tables[tableName].tags
      numMeth = len(methodTableTags['Entry_ID'])

      for i in range(0,numMeth):

        self.structMeths.append(NmrStarStructMeth(self))
        self.structMeths[-1].setData(methodTableTags,i)

class NmrStarMeasureExpt:

  def __init__(self,parent):

    self.parent = parent

  def setData(self,exptTableTags,i):

    assignList = [['exptId',      'Experiment_ID', None],
                  ['exptName',    'Experiment_name', None],
                  ['sampleId',    'Sample_ID', None],
                  ['sampleState', 'Sample_state', None]
                  ]

    for (attrName,tagName,default) in assignList:

      if exptTableTags.has_key(tagName):
        if exptTableTags[tagName][i] != None:
          setattr(self,attrName,exptTableTags[tagName][i])
        else:
          setattr(self,attrName,default)

class NmrStarMeasureSoft:

  def __init__(self,parent):

    self.parent = parent

  def setData(self,softTableTags,i):

    assignList = [['softId',      'Software_ID', None],
                  ['softName',    'Software_label', None],
                  ['methodId',    'Method_ID', None],
                  ['methodLabel', 'Method_label', None]
                  ]

    for (attrName,tagName,default) in assignList:

      if softTableTags.has_key(tagName):
        if softTableTags[tagName][i] != None:
          setattr(self,attrName,softTableTags[tagName][i])
        else:
          setattr(self,attrName,default)

class NmrStarStructSoft:

  def __init__(self,parent):

    self.parent = parent

  def setData(self,softTableTags,i):

    assignList = [['softId',      'Software_ID', None],
                  ['softName',    'Software_label', None],
                  ['methodId',    'Method_ID', None],
                  ['methodLabel', 'Method_label', None]
                  ]

    for (attrName,tagName,default) in assignList:

      if softTableTags.has_key(tagName):
        if softTableTags[tagName][i] != None:
          setattr(self,attrName,softTableTags[tagName][i])
        else:
          setattr(self,attrName,default)

class NmrStarStructMeth:

  def __init__(self,parent):

    self.parent = parent

  def setData(self,methTableTags,i):

    assignList = [['methodLabel',   'Refine_method', None],
                  ['methodDetails', 'Refine_details', None]
                  ]

    for (attrName,tagName,default) in assignList:

      if methTableTags.has_key(tagName):
        if methTableTags[tagName][i] != None:
          setattr(self,attrName,methTableTags[tagName][i])
        else:
          setattr(self,attrName,default)

class NmrStarConstraintFile(NmrStarGenericFile):

  def parseCommentsLoop(self):
  
    self.comments = []
      
    #
    # Parse the comments loop (not really necessary but good example)
    #
    
    tableTagNames = {}
    tableName = None
    
    if self.version == '2.1.1':

      tableName = '_Constraint_comment_ID'
      tableTagNames['commentID'] = tableName
      
    elif self.version in ('3.0','3.1'):  
      
      commentTableText = 'constraint_comment_org'
        
      tableName = self.saveFrame.findTableByEndText(commentTableText)
      tableTagNames['commentID'] = 'ID'

    if tableName and self.saveFrame.tables.has_key(tableName):

      constraintTableTags = self.saveFrame.tables[tableName].tags
      
      constraintTableTagKeys = constraintTableTags.keys()
      
      self.commentTags = constraintTableTagKeys

      for i in range(0,len(constraintTableTags[tableTagNames['commentID']])):
        
        commentTags = []
        
        for tagKey in constraintTableTagKeys:
          commentTags.append(constraintTableTags[tagKey][i])
        
        self.comments.append(commentTags)

class NmrStarFile(NmrStarGenericFile):

  def initialize(self):
  
    self.components = ['all']
    
  def setComponentList(self):
  
    pass
    
  def setComponents(self):
        
    self.setComponentList()

    # This should become list of stuff
    self.sfs = {}
    
    #
    # nmrStar DEFINITIONS (to get right value type)
    #
    #  Could be extended to ccpn directly... check this later
    #

    self.sfDict = {}
    
    for component in self.components:
      
      self.setSfDict(component)
      
  def setSfDict(self,component):
  
    #
    # Keep track of component order
    #

    self.componentList = []

    if self.version == '2.1.1':
  
      value = 0

      #
      # This loop occurs multiple times: defined here once
      #

      constraintCommentLoop = [

            '_Constraint_comment_ID', [        

                ['_Constraint_comment_ID',           lambda x = value: returnStarInt(x)],
                ['_Constraint_comment',              lambda x = value: returnStarString(x)],
                ['_Constraint_comment_begin_line',   lambda x = value: returnStarInt(x)],
                ['_Constraint_comment_begin_column', lambda x = value: returnStarInt(x)],
                ['_Constraint_comment_end_line',     lambda x = value: returnStarInt(x)],
                ['_Constraint_comment_end_column',   lambda x = value: returnStarInt(x)]
              ]
            ]

      #
      # Entry information
      #

      sfName = 'entry_information'

      if component == sfName or component == 'all':
 
        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Entry_title',          lambda x = value: x],
            #['_Version_type',         lambda x = value: returnStarString(x)],
            ['_NMR_STAR_version',     lambda x = value: returnStarString(x)],
            ['_Submission_date',      lambda x = value: returnStarString(x)],
            ['_Accession_date',       lambda x = value: returnStarString(x)],
            ['_Experimental_method',  lambda x = value: returnStarString(x)]
            
          ]

        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_Author_ordinal', [        

                ['_Author_ordinal',           lambda x = value: returnStarInt(x)],
                ['_Author_family_name',       lambda x = value: returnStarString(x)],
                ['_Author_given_name',        lambda x = value: returnStarString(x)],
                ['_Author_middle_initials',   lambda x = value: returnStarString(x)],
                ['_Author_family_title',      lambda x = value: returnStarString(x)]
              ]
            ],

            ['_Related_BMRB_accession_number', [        

                ['_Related_BMRB_accession_number',   lambda x = value: returnStarInt(x)],
                ['_Relationship',                    lambda x = value: returnStarString(x)],
              ]
            ]

          ]


      #
      # Entry citation
      #

      sfName = 'entry_citation'

      if component == sfName or component == 'all':
 
        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Citation_title',    lambda x = value: x],
            ['_Citation_status',   lambda x = value: returnStarString(x)],
            ['_Citation_type',     lambda x = value: returnStarString(x)]
            
          ]

        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [
            # Can add author info here if required...
            ['_Keyword', [        
                ['_Keyword',  lambda x = value: returnStarString(x)]
              ]
            ]
          ]

      #
      # Constraint file definitions Jurgen MR format
      #

      sfName = 'MR_file_characteristics'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Saveframe_category',             lambda x = value: x],
            ['_BMRB_dev_PDB_id',                lambda x = value: returnStarString(x)]
          ]


        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_BMRB_dev_Position_block_in_MR_file', [        

                ['_BMRB_dev_Position_block_in_MR_file',  lambda x = value: returnStarInt(x)],
                ['_BMRB_dev_Program',                    lambda x = value: returnStarString(x)],
                ['_BMRB_dev_Type',                       lambda x = value: returnStarString(x)],
                ['_BMRB_dev_Subtype',                    lambda x = value: returnStarString(x)],
                ['_BMRB_dev_Format',                     lambda x = value: returnStarString(x)]
              ]
            ]
          ]

      #
      # Comments Jurgen MR format
      #

      sfName = 'MR_file_comments'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Saveframe_category',                 lambda x = value: returnStarString(x)],
            ['_BMRB_dev_Position_block_in_MR_file', lambda x = value: returnStarInt(x)],
            ['_BMRB_dev_Comment',                   lambda x = value: returnStarString(x)],

          ]
        
      #
      # Molecular system info
      #
      # Warning: not complete! Quick hack to get info out of 2.1.1 files...
      #

      sfName = 'molecular_system'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Mol_system_name',                  lambda x = value: x],
            ['_Abbreviation_common',              lambda x = value: x],
            ['_Enzyme_commission_number',         lambda x = value: x],
            ['_System_physical_state',            lambda x = value: x],
            ['_System_oligomer_state',            lambda x = value: x],
            ['_System_paramagnetic',              lambda x = value: x],
            ['_System_thiol_state',               lambda x = value: x]

          ]


        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_Mol_system_component_name', [        

                ['_Mol_system_component_name',   lambda x = value: x],
                ['_Mol_label',                   lambda x = value: x]

              ]
            ],

            ['_Biological_function', [        

                ['_Biological_function',                          lambda x = value: x]
      
              ]
            ]
          ]

      #
      # Sequence for a monomeric polymer
      #

      sfName = 'monomeric_polymer'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Mol_type',                       lambda x = value: x],
            ['_Mol_polymer_class',              lambda x = value: x],
            ['_Name_common',                    lambda x = value: x],
            ['_Name_variant',                   lambda x = value: x],
            ['_Name_common',                    lambda x = value: x],
            ['_Abbreviation_common',            lambda x = value: x],
            ['_Molecular_mass',                 lambda x = value: returnStarFloat(x)],
            ['_Mol_thiol_state',                lambda x = value: x],
            ['_Residue_count',                  lambda x = value: returnStarInt(x)],
            ['_Mol_residue_sequence',           lambda x = value: x]
          ]


        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_Residue_seq_code', [        

                ['_Residue_seq_code',           lambda x = value: returnStarInt(x)],
                ['_Residue_author_seq_code',    lambda x = value: returnStarInt(x)],
                ['_Residue_label',              lambda x = value: x]
              ]
            ],

            ['_Database_name', [        

                ['_Database_name',                          lambda x = value: x],
                ['_Database_accession_code',                lambda x = value: x],
                ['_Database_entry_mol_name',                lambda x = value: x],
                ['_Sequence_query_to_submitted_percentage', lambda x = value: returnStarFloat(x)],
                ['_Sequence_subject_length',                lambda x = value: returnStarInt(x)],
                ['_Sequence_identity',                      lambda x = value: returnStarFloat(x)],
                ['_Sequence_positive',                      lambda x = value: returnStarFloat(x)],
                ['_Sequence_homology_expectation_value',    lambda x = value: returnStarFloat(x)]
      
              ]
            ]
          ]

      #
      # Ligands
      #

      sfName = 'ligand'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Mol_type',                       lambda x = value: x],
            ['_Name_common',                    lambda x = value: x],
            ['_Abbreviation_common',            lambda x = value: x],
            ['_Name_IUPAC',                     lambda x = value: x],
            ['_BMRB_code',                      lambda x = value: x],
            ['_PDB_code',                       lambda x = value: x],
            ['_Mol_empirical_formula',          lambda x = value: x],
            ['_Mol_paramagnetic',               lambda x = value: x],
            ['_Molecular_mass',                 lambda x = value: returnStarFloat(x)],
            ['_Mol_charge',                     lambda x = value: returnStarInt(x)],
            ['_Details',                        lambda x = value: x]
          ]

      #
      # Assigned chemical shifts
      #

      sfName = 'assigned_chemical_shifts'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Saveframe_category',             lambda x = value: x],
            ['_Mol_system_component_name',      lambda x = value: x],
            ['_Sample_conditions_label',        lambda x = value: x],
            ['_Chem_shift_reference_set_label', lambda x = value: x]
          ]


        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_Sample_label', [        

                ['_Sample_label',      lambda x = value: x]
                
              ]
            ],

            ['_Atom_shift_assign_ID', [        

                ['_Atom_shift_assign_ID',      lambda x = value: returnStarInt(x)],
                ['_Residue_author_seq_code',   lambda x = value: returnStarInt(x)],
                ['_Residue_seq_code',          lambda x = value: returnStarInt(x)],
                ['_Residue_label',             lambda x = value: x],
                ['_Atom_name',                 lambda x = value: x],
                ['_Atom_type',                 lambda x = value: x],
                ['_Chem_shift_value',          lambda x = value: returnStarFloat(x)],
                ['_Chem_shift_value_error',    lambda x = value: returnStarFloat(x)],
                ['_Chem_shift_ambiguity_code', lambda x = value: returnStarInt(x)],
                ['_Chem_shift_ambiguity_type', lambda x = value: returnStarInt(x)]    # This also sometimes occurs...
              ]
            ]
          ]

      #
      # J couplings
      #

      sfName = 'coupling_constants'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Saveframe_category',             lambda x = value: x],
            ['_Mol_system_component_name',      lambda x = value: x],
            ['_Sample_conditions_label',        lambda x = value: x],
          ]


        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_Sample_label', [        

                ['_Sample_label',      lambda x = value: x]
                
              ]
            ],

            ['_Coupling_constant_code', [        

                ['_Coupling_constant_code',    lambda x = value: x],
                ['_Atom_one_residue_seq_code', lambda x = value: returnStarInt(x)],
                ['_Atom_one_residue_label',    lambda x = value: x],
                ['_Atom_one_name',             lambda x = value: x],
                ['_Atom_two_residue_seq_code', lambda x = value: returnStarInt(x)],
                ['_Atom_two_residue_label',    lambda x = value: x],
                ['_Atom_two_name',             lambda x = value: x],
                ['_Coupling_constant_value',       lambda x = value: returnStarFloat(x)],
                ['_Coupling_constant_value_error', lambda x = value: returnStarFloat(x)],
                ['_Coupling_constant_min_value',   lambda x = value: returnStarFloat(x)],
                ['_Coupling_constant_max_value',   lambda x = value: returnStarFloat(x)],
              ]
            ]
          ]


      #
      # Order parameters
      #

      sfName = 'S2_parameters'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Saveframe_category',             lambda x = value: x],
            ['_Mol_system_component_name',      lambda x = value: x],
            ['_Sample_conditions_label',        lambda x = value: x],
          ]


        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_Sample_label', [        

                ['_Sample_label',      lambda x = value: x]
                
              ]
            ],

            ['_Software_label', [        

                ['_Software_label',      lambda x = value: x]
                
              ]
            ],

            ['_Residue_seq_code', [        

                ['_Residue_seq_code', lambda x = value: returnStarInt(x)],
                ['_Residue_label',    lambda x = value: x],
                ['_Atom_name',        lambda x = value: x],
                ['_Model_fit',        lambda x = value: x],
                ['_S2_value',              lambda x = value: returnStarFloat(x)],
                ['_S2_value_fit_error',    lambda x = value: returnStarFloat(x)],
                ['_Tau_e_value',           lambda x = value: returnStarFloat(x)],
                ['_Tau_e_value_fit_error', lambda x = value: returnStarFloat(x)],
                ['_S2f_value',             lambda x = value: returnStarFloat(x)],
                ['_S2f_value_fit_error',   lambda x = value: returnStarFloat(x)],
                ['_S2s_value',             lambda x = value: returnStarFloat(x)],
                ['_S2s_value_fit_error',   lambda x = value: returnStarFloat(x)],
                ['_Tau_s_value',           lambda x = value: returnStarFloat(x)],
                ['_Tau_s_value_fit_error', lambda x = value: returnStarFloat(x)],
                ['_S2H_value',             lambda x = value: returnStarFloat(x)],
                ['_S2H_value_fit_error',   lambda x = value: returnStarFloat(x)],
                ['_S2N_value',             lambda x = value: returnStarFloat(x)],
                ['_S2N_value_fit_error',   lambda x = value: returnStarFloat(x)],
              ]
            ]
          ]

      #
      # Sample conditions
      #

      sfName = 'sample_conditions'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [
            ['_Saveframe_category',             lambda x = value: x]
          ]


        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_Variable_type', [        

                ['_Variable_type',        lambda x = value: x],
                ['_Variable_value',       lambda x = value: returnStarFloat(x)],
                ['_Variable_value_error', lambda x = value: returnStarFloat(x)],
                ['_Variable_value_units', lambda x = value: x]
                
              ]
            ]
          ]

      #
      # Sample 
      #

      sfName = 'sample'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [
            ['_Saveframe_category',             lambda x = value: x],
            ['_Sample_type',                    lambda x = value: x]
          ]


        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_Mol_label', [        

                ['_Mol_label',                 lambda x = value: x],
                ['_Concentration_value',       lambda x = value: returnStarFloat(x)],
                ['_Concentration_value_units', lambda x = value: x],
                ['_Concentration_min_value',   lambda x = value: returnStarFloat(x)],
                ['_Concentration_max_value',   lambda x = value: returnStarFloat(x)],
                ['_Isotopic_labeling',         lambda x = value: x]
                
              ]
            ]
          ]

      #
      # Chemical shift reference
      #

      sfName = 'chemical_shift_reference'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Saveframe_category',  lambda x = value: x],
            ['_Details',             lambda x = value: x],
         ]

        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_Mol_common_name', [        

                ['_Mol_common_name',                    lambda x = value: x],
                ['_Atom_type',                          lambda x = value: x],
                ['_Atom_isotope_number',                lambda x = value: returnStarInt(x)],
                ['_Atom_group',                         lambda x = value: x],
                ['_Chem_shift_units',                   lambda x = value: x],
                ['_Chem_shift_value',                   lambda x = value: returnStarFloat(x)],
                ['_Reference_method',                   lambda x = value: x],
                ['_Reference_type',                     lambda x = value: x],
                ['_External_reference_sample_geometry', lambda x = value: x],
                ['_External_reference_location',        lambda x = value: x],
                ['_External_reference_axis',            lambda x = value: x],
                ['_Indirect_shift_ratio',               lambda x = value: returnStarFloat(x)]
              ]
            ]
          ]

      #
      # Distance constraints (Jurgen format)
      #

      sfName = 'distance_constraints'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Saveframe_category',             lambda x = value: x],
            ['_BMRB_dev_Position_block_in_MR_file', lambda x = value: returnStarInt(x)]
          ]


        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_Constraint_ID', [        

                ['_Constraint_ID',                   lambda x = value: returnStarInt(x)],
                ['_Constraint_tree_node_ID',         lambda x = value: returnStarInt(x)],
                ['_Constraint_tree_logic_operation', lambda x = value: returnStarString(x,strip=1)],
                ['_Constraint_tree_down_node_ID',    lambda x = value: returnStarInt(x)],
                ['_Constraint_tree_right_node_ID',   lambda x = value: returnStarInt(x)]
              ]
            ],

            ['_Constraint_tree_node_member_constraint_ID', [        

                ['_Constraint_tree_node_member_constraint_ID', lambda x = value: returnStarInt(x)],
                ['_Constraint_tree_node_member_node_ID',       lambda x = value: returnStarInt(x)],
                ['_Constraint_tree_node_member_ID',            lambda x = value: returnStarInt(x)],
                ['_Contribution_fractional_value',             lambda x = value: returnStarString(x)],
                ['_Mol_system_component_code',                 lambda x = value: returnStarString(x,strip=1)],
                ['_Residue_seq_code',                          lambda x = value: returnStarInt(x)],
                ['_Atom_name',                                 lambda x = value: returnStarString(x,strip=1)]
              ]
            ],

            ['_Distance_constraint_ID', [        

                ['_Distance_constraint_ID',           lambda x = value: returnStarInt(x)],
                ['_Distance_constraint_tree_node_ID', lambda x = value: returnStarInt(x)],
                ['_Distance_value',                   lambda x = value: returnStarFloat(x)],
                ['_Distance_lower_bound_value',       lambda x = value: returnStarFloat(x)],
                ['_Distance_upper_bound_value',       lambda x = value: returnStarFloat(x)]
              ]
            ],

            constraintCommentLoop
          ]

      #
      # Dihedral constraints (Jurgen format or old?)
      #

      sfName = 'torsion_angle_constraints'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Saveframe_category',                 lambda x = value: x],
            ['_BMRB_dev_Position_block_in_MR_file', lambda x = value: returnStarInt(x)],
            ['_BMRB_dev_Program',                   lambda x = value: x],
            ['_BMRB_dev_Type',                      lambda x = value: x],
            ['_BMRB_dev_Subtype',                   lambda x = value: x],
            ['_BMRB_dev_Format',                    lambda x = value: x]
          ]


        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_Constraint_ID', [        

                ['_Constraint_ID',                    lambda x = value: returnStarInt(x)],
                ['_Angle_name',                       lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_one_mol_system_component_ID', lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_one_residue_seq_code',        lambda x = value: returnStarInt(x)],
                ['_Atom_one_residue_label',           lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_one_atom_name',             lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_two_mol_system_component_ID', lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_two_residue_seq_code',        lambda x = value: returnStarInt(x)],
                ['_Atom_two_residue_label',           lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_two_atom_name',             lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_three_mol_system_component_ID',lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_three_residue_seq_code',      lambda x = value: returnStarInt(x)],
                ['_Atom_three_residue_label',         lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_three_atom_name',             lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_four_mol_system_component_ID',lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_four_residue_seq_code',       lambda x = value: returnStarInt(x)],
                ['_Atom_four_residue_label',          lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_four_atom_name',            lambda x = value: returnStarString(x,strip=1)],
                ['_Angle_upper_bound_value',          lambda x = value: returnStarFloat(x)],
                ['_Angle_lower_bound_value',          lambda x = value: returnStarFloat(x)],
                ['_Force_constant_value',             lambda x = value: returnStarFloat(x)],
                ['_Potential_function_exponent',      lambda x = value: returnStarInt(x)]

              ]
            ],

            constraintCommentLoop
          ]

      #
      # Residual dipolar coupling values
      #

      sfName = 'residual_dipolar_couplings'

      if component == sfName or component == 'all':

        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['_Saveframe_category',             lambda x = value: x],
            ['_Sample_conditions_label',        lambda x = value: x],
            ['_Sample_conditions',              lambda x = value: x],   # Warning; this is not correct but often occurs, so hacking it in...
            ['_Spectrometer_frequency_1H',      lambda x = value: returnStarFloat(x)],
          ]


        # Warning: key here is tagName by which table is identified!
        self.sfDict[sfName].tables = [

            ['_Sample_label', [        

                ['_Sample_label',      lambda x = value: x]
                
              ]
            ],

            ['_Residual_dipolar_coupling_ID', [        

                ['_Residual_dipolar_coupling_ID',     lambda x = value: returnStarInt(x)],
                ['_Residual_dipolar_coupling_code',   lambda x = value: returnStarString(x)],
                ['_Atom_one_mol_system_component_name',lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_one_residue_seq_code',        lambda x = value: returnStarInt(x)],
                ['_Atom_one_residue_label',           lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_one_atom_name',               lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_two_mol_system_component_name',lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_two_residue_seq_code',        lambda x = value: returnStarInt(x)],
                ['_Atom_two_residue_label',           lambda x = value: returnStarString(x,strip=1)],
                ['_Atom_two_atom_name',               lambda x = value: returnStarString(x,strip=1)],
                ['_Residual_dipolar_coupling_value',  lambda x = value: returnStarFloat(x)],
                ['_Residual_dipolar_coupling_error',  lambda x = value: returnStarFloat(x)],
                ['_Residual_dipolar_coupling_value_error',  lambda x = value: returnStarFloat(x)],
                ['_Residual_dipolar_coupling_lower_bound_value',lambda x = value: returnStarFloat(x)],
                ['_Residual_dipolar_coupling_upper_bound_value',lambda x = value: returnStarFloat(x)]

              ]
            ],

            constraintCommentLoop
          ]

    ###################################
    #                                 #
    # START OF VERSION 3.0 STUFF!!    #
    #                                 #
    ###################################

    elif self.version == '3.0':

      value = 0
      
      #
      # TODO: set this up from v3.0 dictionary!!
      #
      # Also: automatically reocurring tags in loops
      # should be handled from sf level (and not defined in dict,
      # or defined by a function when writing!!!)
      #
    
      #
      # Constraint file definitions version 3.0
      #

      sfName = 'file_characteristics'

      if component == sfName or component == 'all':
 
        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['Sf_category',          lambda x = value: x],
            ['ID',                   lambda x = value: returnStarInt(x)],
            ['NMR_STAR_version',     lambda x = value: returnStarString(x)],
            ['PDB_ID',               lambda x = value: returnStarString(x)],
            ['MR_file_name',         lambda x = value: returnStarString(x)],
            ['Conversion_date',      lambda x = value: returnStarString(x)]
            
          ]


        # Warning: key here is prefix in table tag names!
        self.sfDict[sfName].tables = [

            ['_File_characteristic', [        

                ['File_characteristics_ID',  lambda x = value: returnStarInt(x)],
                ['MR_file_block_position',  lambda x = value: returnStarInt(x)],
                ['Program',    lambda x = value: returnStarString(x)],
                ['Type',    lambda x = value: returnStarString(x)],
                ['Subtype',    lambda x = value: returnStarString(x)],
                ['Format',    lambda x = value: returnStarString(x)]
              ]
            ]
          ]

      #
      # Molecular assembly definitions
      #

      sfName = 'assembly'

      if component == sfName or component == 'all':
 
        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['Sf_category',   lambda x = value: x],
            ['ID',            lambda x = value: returnStarInt(x)],
            ['Name',          lambda x = value: returnStarString(x)]
            
          ]

        self.sfDict[sfName].tables = [

            ['_Entity_assembly', [        

                ['Assembly_ID',             lambda x = value: returnStarInt(x)],
                ['ID',    lambda x = value: returnStarInt(x)],
                ['Entity_ID',    lambda x = value: returnStarInt(x)],
                ['Entity_assembly_name',    lambda x = value: returnStarString(x)],
                ['Entity_label',    lambda x = value: returnStarString(x)]
              ]
            ]
          ]


      #
      # Entity information
      #

      sfName = 'entity'

      if component == sfName or component == 'all':
 
        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['Sf_category',   lambda x = value: x],
            ['ID',            lambda x = value: returnStarInt(x)],
            ['Type',          lambda x = value: x],
            ['Pol_type',      lambda x = value: x],
            ['Seq_length',    lambda x = value: returnStarInt(x)],
            ['Seq',           lambda x = value: x]
          ]

        self.sfDict[sfName].tables = [
        
            #
            # Describes moiety
            #

            ['_Entity_comp_index', [        

                ['Entity_ID', lambda x = value: returnStarInt(x)],
                ['Num',       lambda x = value: returnStarInt(x)],
                ['Comp_ID',   lambda x = value: x]
              ]
            ],
       
            #
            # Describes sequence
            #

            ['_Entity_poly_seq', [        

                ['Entity_ID',       lambda x = value: returnStarInt(x)],
                ['Num',             lambda x = value: returnStarInt(x)],
                ['Comp_index_num',  lambda x = value: returnStarInt(x)],
                ['Mon_ID',          lambda x = value: x]
              ]
            ]

          ]

      #
      # ChemComp information
      #

      sfName = 'chem_comp'

      if component == sfName or component == 'all':
 
        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)
        
        #
        # NOTE that IDs do NOT have to be integers!!
        #

        self.sfDict[sfName].tags = [

            ['Sf_category',   lambda x = value: x],
            ['ID',            lambda x = value: x], # is mmCif code
            ['Type',          lambda x = value: x]  # would generally be non-polymer
          ]

        self.sfDict[sfName].tables = [
        
            #
            # Describes moiety
            #

            ['_Chem_comp_atom', [        

                ['Comp_ID',       lambda x = value: x],
                ['Atom_ID',       lambda x = value: x],
                ['Type_symbol',   lambda x = value: x]
              ]
            ]
          ]

      #
      # Distance constraints
      #

      sfName = 'distance_constraints'

      if component == sfName or component == 'all':
 
        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)
        self.sfDict[sfName].preTagShort = '_Dist_'

        self.sfDict[sfName].tags = [

            ['Sf_category',            lambda x = value: x],
            ['ID',                     lambda x = value: returnStarInt(x)],
            ['MR_file_block_position', lambda x = value: returnStarInt(x)],
            ['Program',                lambda x = value: x],
            ['Type',                   lambda x = value: x],
            ['Subtype',                lambda x = value: x],
            ['Format',                 lambda x = value: x]
          ]

        self.sfDict[sfName].tables = [

            ['_Dist_constraint_tree', [        

                ['Constraints_ID',  lambda x = value: returnStarInt(x)],
                ['ID',              lambda x = value: returnStarInt(x)],
                ['Node_ID',         lambda x = value: returnStarInt(x)],
                ['Down_node_ID',    lambda x = value: returnStarInt(x)],
                ['Right_node_ID',   lambda x = value: returnStarInt(x)],
                ['Logic_operation', lambda x = value: returnStarString(x,strip=1)]
              ]
            ],

            ['_Dist_constraint', [        

                ['Constraints_ID',                 lambda x = value: returnStarInt(x)],
                ['Dist_constraint_tree_ID',        lambda x = value: returnStarInt(x)],
                ['Tree_node_member_node_ID',       lambda x = value: returnStarInt(x)],
                ['Contribution_fractional_val',    lambda x = value: returnStarString(x)],
                ['Constraint_tree_node_member_ID', lambda x = value: returnStarInt(x)],
                ['Label_entity_assembly_ID',       lambda x = value: x],
                ['Label_entity_ID',                lambda x = value: x],
                ['Label_comp_index_ID',            lambda x = value: returnStarInt(x)],
                ['Label_comp_ID',                  lambda x = value: returnStarString(x)],
                ['Label_atom_ID',                  lambda x = value: returnStarString(x)],
                ['Auth_segment_code',              lambda x = value: returnStarString(x,strip=1)],
                ['Auth_seq_ID',                    lambda x = value: returnStarString(x)],
                ['Auth_comp_ID',                   lambda x = value: returnStarString(x,strip=1)],
                ['Auth_atom_ID',                   lambda x = value: returnStarString(x,strip=1)]
              ]
            ],

            ['_Dist_constraint_value', [        

                ['Constraints_ID',            lambda x = value: returnStarInt(x)],
                ['Constraint_ID',             lambda x = value: returnStarInt(x)],
                ['Tree_node_ID',              lambda x = value: returnStarInt(x)],
                ['Source_experiment_ID',      lambda x = value: x],
                ['Spectral_peak_ID',          lambda x = value: returnStarInt(x)],
                ['Intensity_val',             lambda x = value: returnStarFloat(x)],
                ['Intensity_lower_val_err',   lambda x = value: returnStarFloat(x)],
                ['Intensity_upper_val_err',   lambda x = value: returnStarFloat(x)],
                ['Distance_val',              lambda x = value: returnStarFloat(x)],
                ['Distance_lower_bound_val',  lambda x = value: returnStarFloat(x)],
                ['Distance_upper_bound_val',  lambda x = value: returnStarFloat(x)],
                ['Weight',                    lambda x = value: returnStarFloat(x)],
                ['Spectral_peak_ppm_1',       lambda x = value: returnStarFloat(x)],
                ['Spectral_peak_ppm_2',       lambda x = value: returnStarFloat(x)]
              ]
            ],

            self.getConstraintCommentLoop(preTag = self.sfDict[sfName].preTagShort),
            
            self.getConstraintParseErrorLoop(preTag = self.sfDict[sfName].preTagShort),
            
            self.getConstraintParseFileLoop(preTag = self.sfDict[sfName].preTagShort),
            
            self.getConstraintConversionErrorLoop(preTag = self.sfDict[sfName].preTagShort)

          ]

      #
      # Dihedral constraints
      #

      sfName = 'torsion_angle_constraints'

      if component == sfName or component == 'all':
 
        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)
        self.sfDict[sfName].preTagShort = '_TA_'

        self.sfDict[sfName].tags = [

            ['Sf_category',            lambda x = value: x],
            ['ID',                     lambda x = value: returnStarInt(x)],
            ['MR_file_block_position', lambda x = value: returnStarInt(x)],
            ['Program',                lambda x = value: x],
            ['Type',                   lambda x = value: x],
            ['Subtype',                lambda x = value: x],
            ['Format',                 lambda x = value: x]
          ]

        self.sfDict[sfName].tables = [

            ['_Torsion_angle_constraint', [        

                ['Constraints_ID',               lambda x = value: returnStarInt(x)],
                ['ID',                           lambda x = value: returnStarInt(x)],
                ['Torsion_angle_name',           lambda x = value: returnStarString(x)],
                
                ['Label_entity_assembly_ID_1',   lambda x = value: x],
                ['Label_entity_ID_1',            lambda x = value: x],
                ['Label_comp_index_ID_1',        lambda x = value: returnStarInt(x)],
                ['Label_comp_ID_1',              lambda x = value: returnStarString(x)],
                ['Label_atom_ID_1',              lambda x = value: returnStarString(x)],
                ['Label_entity_assembly_ID_2',   lambda x = value: x],
                ['Label_entity_ID_2',            lambda x = value: x],
                ['Label_comp_index_ID_2',        lambda x = value: returnStarInt(x)],
                ['Label_comp_ID_2',              lambda x = value: returnStarString(x)],
                ['Label_atom_ID_2',              lambda x = value: returnStarString(x)],
                ['Label_entity_assembly_ID_3',   lambda x = value: x],
                ['Label_entity_ID_3',            lambda x = value: x],
                ['Label_comp_index_ID_3',        lambda x = value: returnStarInt(x)],
                ['Label_comp_ID_3',              lambda x = value: returnStarString(x)],
                ['Label_atom_ID_3',              lambda x = value: returnStarString(x)],
                ['Label_entity_assembly_ID_4',   lambda x = value: x],
                ['Label_entity_ID_4',            lambda x = value: x],
                ['Label_comp_index_ID_4',        lambda x = value: returnStarInt(x)],
                ['Label_comp_ID_4',              lambda x = value: returnStarString(x)],
                ['Label_atom_ID_4',              lambda x = value: returnStarString(x)],

                ['Auth_segment_code_1',            lambda x = value: returnStarString(x,strip=1)],
                ['Auth_seq_ID_1',                  lambda x = value: returnStarString(x)],
                ['Auth_comp_ID_1',                 lambda x = value: returnStarString(x,strip=1)],
                ['Auth_atom_ID_1',                 lambda x = value: returnStarString(x,strip=1)],
                ['Auth_segment_code_2',            lambda x = value: returnStarString(x,strip=1)],
                ['Auth_seq_ID_2',                  lambda x = value: returnStarString(x)],
                ['Auth_comp_ID_2',                 lambda x = value: returnStarString(x,strip=1)],
                ['Auth_atom_ID_2',                 lambda x = value: returnStarString(x,strip=1)],
                ['Auth_segment_code_3',            lambda x = value: returnStarString(x,strip=1)],
                ['Auth_seq_ID_3',                  lambda x = value: returnStarString(x)],
                ['Auth_comp_ID_3',                 lambda x = value: returnStarString(x,strip=1)],
                ['Auth_atom_ID_3',                 lambda x = value: returnStarString(x,strip=1)],
                ['Auth_segment_code_4',            lambda x = value: returnStarString(x,strip=1)],
                ['Auth_seq_ID_4',                  lambda x = value: returnStarString(x)],
                ['Auth_comp_ID_4',                 lambda x = value: returnStarString(x,strip=1)],
                ['Auth_atom_ID_4',                 lambda x = value: returnStarString(x,strip=1)],

                ['Angle_upper_bound_val',          lambda x = value: returnStarFloat(x)],
                ['Angle_lower_bound_val',          lambda x = value: returnStarFloat(x)],
                ['Force_constant_value',           lambda x = value: returnStarFloat(x)],
                ['Potential_function_exponent',    lambda x = value: returnStarFloat(x)]
              ]
            ],

            ['_Dist_constraint', [        

                ['Constraints_ID',                 lambda x = value: returnStarInt(x)],
                ['Dist_constraint_tree_ID',        lambda x = value: returnStarInt(x)],
                ['Tree_node_member_constraint_ID', lambda x = value: returnStarInt(x)],
                ['Tree_node_member_node_ID',       lambda x = value: returnStarInt(x)],
                ['Contribution_fractional_value',  lambda x = value: returnStarString(x)],
                ['Constraint_tree_node_member_ID', lambda x = value: returnStarInt(x)],
                ['Label_entity_assembly_ID',       lambda x = value: x],
                ['Label_entity_ID',                lambda x = value: x],
                ['Label_comp_index_ID',            lambda x = value: x],
                ['Label_comp_ID',                  lambda x = value: returnStarString(x)],
                ['Label_atom_ID',                  lambda x = value: returnStarString(x)],
                ['Auth_segment_code',              lambda x = value: x],
                ['Auth_seq_ID',                    lambda x = value: returnStarString(x)],
                ['Auth_comp_ID',                   lambda x = value: returnStarString(x)],
                ['Auth_atom_ID',                   lambda x = value: returnStarString(x)]
              ]
            ],


            self.getConstraintCommentLoop(preTag = self.sfDict[sfName].preTagShort),
            
            self.getConstraintParseErrorLoop(preTag = self.sfDict[sfName].preTagShort),
            
            self.getConstraintParseFileLoop(preTag = self.sfDict[sfName].preTagShort),
            
            self.getConstraintConversionErrorLoop(preTag = self.sfDict[sfName].preTagShort)

          ]


      #
      # RDC constraints
      #

      sfName = 'residual_dipolar_couplings'

      if component == sfName or component == 'all':
 
        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName, prefix = '_RDC_constraints')
        self.sfDict[sfName].preTagShort = '_RDC_'

        self.sfDict[sfName].tags = [

            ['Sf_category',            lambda x = value: x],
            ['ID',                     lambda x = value: returnStarInt(x)],
            ['MR_file_block_position', lambda x = value: returnStarInt(x)],
            ['Program',                lambda x = value: x],
            ['Type',                   lambda x = value: x],
            ['Subtype',                lambda x = value: x],
            ['Format',                 lambda x = value: x]
          ]

        self.sfDict[sfName].tables = [

            ['_RDC_constraint', [        

                ['Constraints_ID',               lambda x = value: returnStarInt(x)],
                ['ID',                           lambda x = value: returnStarInt(x)],
                ['Label_entity_assembly_ID_1',   lambda x = value: x],
                ['Label_entity_ID_1',            lambda x = value: x],
                ['Label_comp_index_ID_1',        lambda x = value: returnStarInt(x)],
                ['Label_comp_ID_1',              lambda x = value: returnStarString(x)],
                ['Label_atom_ID_1',              lambda x = value: returnStarString(x)],
                ['Label_entity_assembly_ID_2',   lambda x = value: x],
                ['Label_entity_ID_2',            lambda x = value: x],
                ['Label_comp_index_ID_2',        lambda x = value: returnStarInt(x)],
                ['Label_comp_ID_2',              lambda x = value: returnStarString(x)],
                ['Label_atom_ID_2',              lambda x = value: returnStarString(x)],

                ['Auth_segment_code_1',            lambda x = value: returnStarString(x,strip=1)],
                ['Auth_seq_ID_1',                  lambda x = value: returnStarString(x)],
                ['Auth_comp_ID_1',                 lambda x = value: returnStarString(x,strip=1)],
                ['Auth_atom_ID_1',                 lambda x = value: returnStarString(x,strip=1)],
                ['Auth_segment_code_2',            lambda x = value: returnStarString(x,strip=1)],
                ['Auth_seq_ID_2',                  lambda x = value: returnStarString(x)],
                ['Auth_comp_ID_2',                 lambda x = value: returnStarString(x,strip=1)],
                ['Auth_atom_ID_2',                 lambda x = value: returnStarString(x,strip=1)],

                ['RDC_val',          lambda x = value: returnStarFloat(x)],
                ['RDC_lower_bound',  lambda x = value: returnStarFloat(x)],
                ['RDC_upper_bound',  lambda x = value: returnStarFloat(x)],
                ['RDC_val_err',      lambda x = value: returnStarFloat(x)]
              ]
            ],


            self.getConstraintCommentLoop(preTag = self.sfDict[sfName].preTagShort),
            
            self.getConstraintParseErrorLoop(preTag = self.sfDict[sfName].preTagShort),
            
            self.getConstraintParseFileLoop(preTag = self.sfDict[sfName].preTagShort),
            
            self.getConstraintConversionErrorLoop(preTag = self.sfDict[sfName].preTagShort)
            
          ]

      #
      # Comments
      #

      sfName = 'MR_file_comment'

      if component == sfName or component == 'all':
 
        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['Sf_category',            lambda x = value: x],
            ['ID',                     lambda x = value: returnStarInt(x)],
            ['MR_file_block_position', lambda x = value: returnStarInt(x)],
            ['Program',                lambda x = value: x],
            ['Type',                   lambda x = value: x],
            ['Subtype',                lambda x = value: x],
            ['Format',                 lambda x = value: x],
            ['Comment',                lambda x = value: returnStarString(x)]

          ]

      #
      # Coordinates
      #

      sfName = 'conformer_family_coord_set'

      if component == sfName or component == 'all':
 
        self.componentList.append(sfName)
        self.sfDict[sfName] = SaveFrame("",sfName)

        self.sfDict[sfName].tags = [

            ['Sf_category',            lambda x = value: x],
            ['ID',                     lambda x = value: returnStarInt(x)],
            ['Details',                lambda x = value: returnStarString(x)]

          ]

        self.sfDict[sfName].tables = [
        
            #
            # Describes coordinates
            #

            ['_Atom_site', [        

                ['Conformer_family_coord_set_ID',   lambda x = value: returnStarInt(x)],
                ['Model_ID',                        lambda x = value: returnStarInt(x)],
                ['ID',                              lambda x = value: returnStarInt(x)],
                
                ['Label_entity_assembly_ID',        lambda x = value: returnStarInt(x)],
                ['Label_entity_ID',                 lambda x = value: returnStarInt(x)],
                ['Label_comp_index_ID',             lambda x = value: returnStarInt(x)],
                ['Label_comp_ID',                   lambda x = value: returnStarString(x)],
                ['Label_atom_ID',                   lambda x = value: returnStarString(x)],

                ['Auth_segment_code',               lambda x = value: returnStarString(x)],
                ['Auth_seq_ID',                     lambda x = value: returnStarString(x)],
                ['Auth_comp_ID',                    lambda x = value: returnStarString(x)],
                ['Auth_atom_ID',                    lambda x = value: returnStarString(x)],

                ['Type_symbol',                     lambda x = value: returnStarString(x)],
                ['Cartn_x',                         lambda x = value: returnStarFloat(x)],
                ['Cartn_y',                         lambda x = value: returnStarFloat(x)],
                ['Cartn_z',                         lambda x = value: returnStarFloat(x)],

                ['PDB_extracted_Bfactor_col',       lambda x = value: returnStarFloat(x)]
               
              ]
            ]
          ]

    #############################################
    #                                           #
    # START OF TEMPORARY VERSION 3.1 STUFF!!    #
    #                                           #
    #############################################

    elif self.version == '3.1':
      
      #
      # In 3.1, just read everything based on nmrStarDict...
      #
      
      value = 0
      
      try:
        from pdbe.nmrStar.IO.nmrStarDict import sfDict as refSfDict
      except:
        raise self.FormatFileError("No NMR-STAR dictionary available for version 3.1!")
      
      for sfName in refSfDict.keys():
        
        #if sfName == 'assigned_chemical_shifts':
        #  debugMode = True
        #else:
        #  debugMode = False
        
        if component == sfName or component == 'all':
          self.componentList.append(sfName)
          self.sfDict[sfName] = SaveFrame("",sfName,prefix = "_%s" % refSfDict[sfName]['name'])
          
          #
          # Set tags - TODO this is very confusing - .tags is a dictionary when initialising SaveFrame. Why hacked?
          #
          
          self.sfDict[sfName].tags = []
          for tagName in refSfDict[sfName]['tagNames']:
            self.sfDict[sfName].tags.append([tagName,refSfDict[sfName]['tags'][tagName][1]])
            
          #
          # Set tables - see tags comment!
          #
          
          self.sfDict[sfName].tables = []
          
          if refSfDict[sfName].has_key('tableNames'):
           for tableName in refSfDict[sfName]['tableNames']:
             self.sfDict[sfName].tables.append(['_%s' % tableName,[]])
             
             for tagName in refSfDict[sfName]['tables'][tableName]['tagNames']:
               self.sfDict[sfName].tables[-1][-1].append([tagName,refSfDict[sfName]['tables'][tableName]['tags'][tagName][1]])

    else:
    
      print "ERROR UNRECOGNIZED VERSION '%s'" % (str(self.version))

  def getConstraintCommentLoop(self,preTag = 'Dist_'):
  
    #
    # This loop occurs multiple times: defined here once
    #
    
    value = 0

    if self.version == '3.0':

      constraintCommentLoop = [

            preTag + 'constraint_comment_org', [        

                ['Constraints_ID',       lambda x = value: returnStarInt(x)],
                ['ID',                   lambda x = value: returnStarInt(x)],
                ['Comment_begin_line',   lambda x = value: returnStarInt(x)],
                ['Comment_begin_column', lambda x = value: returnStarInt(x)],
                ['Comment_end_line',     lambda x = value: returnStarInt(x)],
                ['Comment_end_column',   lambda x = value: returnStarInt(x)],
                ['Comment',              lambda x = value: returnStarString(x)]
              ]
            ]
    
    elif self.version == '3.1':
    
      constraintCommentLoop = [

            preTag + 'constraint_comment_org', [        

                ['ID',                   lambda x = value: returnStarInt(x)],
                ['Comment_text',         lambda x = value: returnStarString(x)],    # *
                ['Comment_begin_line',   lambda x = value: returnStarInt(x)],
                ['Comment_begin_column', lambda x = value: returnStarInt(x)],
                ['Comment_end_line',     lambda x = value: returnStarInt(x)],
                ['Comment_end_column',   lambda x = value: returnStarInt(x)],
              ]
            ]
    
    return constraintCommentLoop

  def getConstraintParseErrorLoop(self,preTag = 'Dist_'):
  
    #
    # This loop occurs multiple times: defined here once
    #
    
    value = 0

    if self.version == '3.0':

      constraintParseErrorLoop = [

            preTag + 'constraint_org_file_parse_err', [        

                ['Constraints_ID',   lambda x = value: returnStarInt(x)],
                ['ID',               lambda x = value: returnStarInt(x)],
                ['Begin_line',       lambda x = value: returnStarInt(x)],
                ['Begin_column',     lambda x = value: returnStarInt(x)],
                ['End_line',         lambda x = value: returnStarInt(x)],
                ['End_column',       lambda x = value: returnStarInt(x)],
                ['Content',          lambda x = value: returnStarString(x)]
              ]
            ]
    
    elif self.version == '3.1':
    
      constraintParseErrorLoop = [

            preTag + 'constraint_parse_err', [        

                ['ID',               lambda x = value: returnStarInt(x)],
                ['Content',          lambda x = value: returnStarString(x)],
                ['Begin_line',       lambda x = value: returnStarInt(x)],
                ['Begin_column',     lambda x = value: returnStarInt(x)],
                ['End_line',         lambda x = value: returnStarInt(x)],
                ['End_column',       lambda x = value: returnStarInt(x)],
              ]
            ]
    
    return constraintParseErrorLoop

  def getConstraintParseFileLoop(self,preTag = 'Dist_'):
  
    #
    # This loop occurs multiple times: defined here once - TODO 3.1
    #
    
    value = 0

    constraintParseFileLoop = [

          preTag + 'constraint_parse_file', [        

              ['Constraints_ID',   lambda x = value: returnStarInt(x)],
              ['ID',               lambda x = value: returnStarInt(x)],
              ['Name',             lambda x = value: returnStarString(x)]
            ]
          ]
    
    return constraintParseFileLoop

  def getConstraintConversionErrorLoop(self,preTag = 'Dist_'):
  
    #
    # This loop occurs multiple times: defined here once - TODO 3.1
    #
    
    value = 0

    constraintConversionErrorLoop = [

          preTag + 'constraint_parse_file_conv_err', [        

              ['Constraints_ID',           lambda x = value: returnStarInt(x)],
              ['ID',                       lambda x = value: returnStarInt(x)],
              ['Parse_file_ID',            lambda x = value: returnStarInt(x)],
              ['Parse_file_sf_label',      lambda x = value: returnStarString(x)],
              ['Parse_file_constraint_ID', lambda x = value: returnStarInt(x)],
              ['Conv_error_type',          lambda x = value: returnStarInt(x)],
              ['Conv_error_note',          lambda x = value: returnStarString(x)]
            ]
          ]
    
    return constraintConversionErrorLoop

  ###################################
  #                                 #
  # END OF VERSION 3 STUFF!!        #
  #                                 #
  ###################################  
  
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
  
    return self.sfDict[saveFrameName].prefix
    
  def setTags(self,valueDict):
  
    for attr in valueDict.keys():
      value = valueDict[attr]
      setattr(self,attr,self.tags[attr][1](value))

  def readComponent(self,verbose = 0):
  
    #
    # WARNING: this is a silly extra level... currently upheld because
    #          fits more easily with code in importNmrChemShifts.
    #

    self.readGeneric(verbose = verbose)

    cyanaLibUsed = 0
    
    if self.sfs.has_key(self.saveFrameName):

      #
      # Reading all saveframes with dihedral constraint data
      #

      for saveFrame in self.sfs[self.saveFrameName]:

        dataClassFile = self.DataClassFile(self.name,self,saveFrame)
        self.files.append(dataClassFile)

        if hasattr(dataClassFile,'cyanaLibUsed') and dataClassFile.cyanaLibUsed:
          cyanaLibUsed = 1

    if not self.files:
      print "  No %s data found in nmrStar file %s" % (self.saveFrameName,self.name)

    if cyanaLibUsed:
      print "Using CYANA library - courtesy of Peter Guentert."

  def readGeneric(self, verbose = 0):

    origNmrStarFile = nmrStar.File(verbosity = 2, filename = self.name)

    print "Python nmrStar reader courtesy of Jurgen Doreleijers (BMRB) - with modifications added."

    # Read star file in one gulp
    if origNmrStarFile.read(text = self.text):
      print "  Error reading nmrStar file %s" % self.name
      return
      
    #
    # Try to determine version...
    #

    versionHits = {'2.1.1': 0, '3.1': 0} # '3.0': 0,

    for origSaveFrame in origNmrStarFile.datanodes:

      for tagtable in origSaveFrame.tagtables:

        if tagtable.tagnames.count('_Saveframe_category'):
          versionHits['2.1.1'] += 1
        
        else:
          for i in range(len(tagtable.tagnames)):
            
            components = tagtable.tagnames[i].split(self.tagSep)

            if len(components) == 2:
              (sfName,tagName) = components

              if tagName == 'Sf_category':

                versionHits['3.1'] += 1
                break
             
    #
    # Reset if necessary
    #
    
    if not self.setVersion(versionHits):
      return
    
    #
    # Go through the data. Keep track of print statements in 
    # conversionerror class... 
    #
    
    convError = ConversionError()

    for origSaveFrame in origNmrStarFile.datanodes:

      saveFrameCat = None

      for tagtable in origSaveFrame.tagtables:

        try:
          
          tagIndex = -1
          
          if self.version == '2.1.1':
            tagIndex = tagtable.tagnames.index('_Saveframe_category')

          elif self.version in ('3.0','3.1'):

            for i in range(0,len(tagtable.tagnames)):
              
              (sfName,tagName) = tagtable.tagnames[i].split(self.tagSep)

              if tagName == 'Sf_category':
              
                tagIndex = i
                break

          if tagIndex != -1:
            saveFrameCat = tagtable.tagvalues[tagIndex][0]
            break
          
        except:
          pass


      if self.sfDict.has_key(saveFrameCat):

        saveFrame = self.setupSaveFrame(saveFrameCat,origSaveFrame.title)
        saveFrameDict = self.sfDict[saveFrameCat]
        
        #
        # FROM HERE ON: rerout print statements to convError object!
        #
        #
        # NOTE: if want to see error messages, best place is in setTag for SaveFrame/Table!
        #
        
        origStdout = sys.stdout
        sys.stdout = convError

        for tagtable in origSaveFrame.tagtables:

          #
          # Here do tags
          #

          if tagtable.free != None:
          
            self.setReadTags(saveFrame,saveFrameDict,tagtable,convError)
             
          #
          # Here do tables
          #

          else:
          
            self.setTableTags(saveFrame,saveFrameDict,tagtable,convError)
       
        sys.stdout = origStdout

    del(origNmrStarFile)

  def setVersion(self,versionHits):
  
    if versionHits['3.1'] > versionHits['2.1.1'] and self.version not in ('3.0','3.1'):
     
      print "  Warning: setting nmrStar version to 3.1 for reading."
      self.version = '3.1'
      self.setComponents()
         
    elif versionHits['2.1.1'] > versionHits['3.1'] and self.version != '2.1.1':
     
      print "  Warning: setting nmrStar version to 2.1.1 for reading."
      self.version = '2.1.1'
      self.setComponents()
      
    return True

  def setReadTags(self,saveFrame,saveFrameDict,tagtable,convError):

    for tagInfo in saveFrameDict.tags:

      (tagName,returnValue) = tagInfo
              
      matchTagName = tagName

      if self.version in ('3.0','3.1'):

        matchTagName = saveFrameDict.prefix + self.tagSep + tagName

      # TODO: Or could make set saveFrameDict in self, pass tagtable, tagName to
      # setTag...

      starValue = getNmrStarValue(tagtable,matchTagName)
      
      if starValue:
        value = returnValue(starValue)
      else:
        #sys.__stdout__.write(" %s, %s\n" % (matchTagName,str(tagtable.tagnames)))
        value = None

      saveFrame.setTag(tagName,value,error = convError.getString())

  def setTableTags(self,saveFrame,saveFrameDict,tagtable,convError):

    for (tableName,tags) in saveFrameDict.tables:

      tableExists = 0

      if self.version == '2.1.1' and tableName in tagtable.tagnames:
        tableExists = 1

      elif self.version in ('3.0','3.1'):

        (tempTableName,tagName) = tagtable.tagnames[0].split(self.tagSep)

        if tableName == tempTableName:
          tableExists = 1

      if tableExists:

        table = saveFrame.setupTable(tableName)
        currentTags = []

        for tagInfo in tags:

          (tagName,returnValue) = tagInfo

          matchTagName = tagName

          if self.version in ('3.0','3.1'):

            matchTagName = tableName + self.tagSep + tagName                  

          if matchTagName in tagtable.tagnames:

            currentTags.append([tagName,returnValue])

            table.tags[tagName] = []
            table.tagErrors[tagName] = []
            table.tagNames.append(tagName)

        for i in range(len(tagtable.tagvalues[0])):

          for tagInfo in currentTags:

            (tagName,returnValue) = tagInfo

            matchTagName = tagName

            if self.version in ('3.0','3.1'):

              matchTagName = tableName + self.tagSep + tagName                  

            starValue = getNmrStarValue(tagtable,matchTagName,i)
            value = returnValue(starValue)

            table.setTag(tagName,value,error = convError.getString())

  def writeGeneric(self, verbose = 0, title = None, topComment = None):

    if not title:
      title = self.name

    outputNmrStarFile = nmrStar.File(verbosity = 2, flavor = 'NMR-STAR', datanodes = [], title = title, filename = self.name)

    #
    # Jurgen credits
    #
    
    print "Python nmrStar writer courtesy of Jurgen Doreleijers (BMRB) - with modifications added."

    #
    # Write out the saveframes, use ordered list...
    #

    for saveFrameName in self.componentList:

      if not self.sfs.has_key(saveFrameName):

        continue

      for saveFrame in self.sfs[saveFrameName]:

        outputNmrStarFile.datanodes.append(nmrStar.SaveFrame(title = saveFrame.title, comment = saveFrame.comment))
        outputSaveFrame = outputNmrStarFile.datanodes[-1]

        #
        # Here do tags
        #

        outputSaveFrame.tagtables.append(nmrStar.TagTable(free = 1,tagnames = [], tagvalues = []))
        outputTagTable = outputSaveFrame.tagtables[-1]

        for tagName in saveFrame.tagNames:
          value = saveFrame.tags[tagName]

          if value != None and value != '':
            value = str(value)

            outputTagTable.tagnames.append(tagName)
            outputTagTable.tagvalues.append([value])

        #
        # Here do tables
        #

        for tableName in saveFrame.tableNames:

          table = saveFrame.tables[tableName]
        
          if not table.tagNames:
            if verbose:
              print "  ERROR: nmrStar output for table %s... no tags available" % tableName
            continue

          outputSaveFrame.tagtables.append(nmrStar.TagTable(free = None,tagnames = [], tagvalues = []))
          outputTagTable = outputSaveFrame.tagtables[-1]

          for tagName in table.tagNames:

            valFlag = False

            for value in table.tags[tagName]:

              if value is not None:
                valFlag = True
                break

            if not valFlag and verbose:
              print '  Warning: tag %s is always None.' % tagName

          for tagName in table.tagNames:
            outputTagTable.tagnames.append(tagName)
            outputTagTable.tagvalues.append([])
            outputTagValue = outputTagTable.tagvalues[-1]

            for value in table.tags[tagName]:

              if value == None or value == '':
                value = '.'
              else:
                value = str(value)

              outputTagValue.append(value)
              
    #
    # Add a top comment, if necessary
    #
    
    if topComment:
      outputNmrStarFile.datanodes[0].comment = topComment

    #
    # Write star file in one go
    #
    
    if outputNmrStarFile.write():
      print "  Error writing nmrStar file %s" % self.name
      return 0

    else:
      return 1

class SaveFrame:

  def __init__(self,title,name,prefix = None):
    
    self.tags = {}
    self.tagErrors = {}
    self.tagNames = []
    self.tables = {}
    self.tableNames = []
    self.title = title
    
    self.comment = ''
    
    self.name = name
    
    if prefix:
    
      self.prefix = prefix
    
    else:

      if name[0] not in string.ascii_uppercase:
        self.prefix = '_' + name.capitalize()
      else:
        self.prefix = '_' + name

  def setupTable(self,tableName):
    
    self.tables[tableName] = Table()
    self.tableNames.append(tableName)
    
    return self.tables[tableName]
    
  def setTag(self,tagName,value,error=None):
    
    if value != None and value != '':
      self.tags[tagName] = value
      self.tagErrors[tagName] = [error]
      
      #if error:
      #  sys.__stdout__.write("%s,%s\n" % (tagName,error))
      
      # Allow for overwrites...
      if tagName not in self.tagNames:
        self.tagNames.append(tagName)

  def findTableByEndText(self,tableEndText):
  
    tableName = None
    
    for tmpTableName in self.tables.keys():
      if tmpTableName.endswith(tableEndText):
        tableName = tmpTableName
        break
        
    return tableName

class Table:

  def __init__(self):
  
    self.tags = {}
    self.tagNames = []
    self.tagErrors = {}

  def setTag(self,tagName,value,error=None):
    
    if not self.tags.has_key(tagName):
      self.tagNames.append(tagName)
      self.tags[tagName] = []
      self.tagErrors[tagName] = []
      
    self.tags[tagName].append(value)
    self.tagErrors[tagName].append(error)

    #if error:
    #  sys.__stdout__.write("%s,%s\n" % (tagName,error))
    
class ConversionError:

  def __init__(self):
  
    self.string = ""
  
  def write(self,error):
    
    error = error.replace(os.linesep,'')
    error.strip()
    
    if error:
      self.string += error + '|'
      
  def getString(self):

    value = self.string[:-1]
    
    self.string = ""

    return value


class GenericConstraint:

  def setErrors(self,errorTags,tableTagNames,i):
    
    for tableTagName in tableTagNames.keys():
    
      if tableTagName == 'position' or not errorTags.has_key(tableTagNames[tableTagName]):
        continue
        
      if tableTagNames[tableTagName].find('%s') == -1:
        
        error = errorTags[tableTagNames[tableTagName]][i]
          
        if error:

          self.errors.append(error)
            
      else:
      
        for position in tableTagNames['position']:
        
          error = errorTags[tableTagNames[tableTagName] % position][i]

          if error:

            self.errors.append(error)    

###################
# Main of program #
###################

if __name__ == "__main__":  
                                                      
  files = [#'../rawdata/8/1/1/info.general/bmr5106.str',
           #'../rawdata/2/2/1/info.general/bmr5307.str',
           #'../reference/nmrStar/bmr5620.str',
           #'../reference/ccpNmr/aartUtrecht/1bsh/restraints.star'#,
           #'../reference/ccpNmr/aartUtrecht/1d8v/restraints.star',
           #'../reference/ccpNmr/jurgenBmrb/1jwe/restraints.star'
           '../reference/ccpNmr/jurgenBmrb/1ao9/1ao9.str'
           ]
  
  for file in files:
    
    file = os.path.join(getTopDirectory(), file)
    
    nmrStarFile = NmrStarFile(file)

    nmrStarFile.readGeneric(verbose = 1)
    
    for key in nmrStarFile.sfs.keys():
    
      print key

      for Sf in nmrStarFile.sfs[key]:

        for tagName in Sf.tags.keys():

          print tagName, Sf.tags[tagName]
          print tagName, Sf.tagErrors[tagName]

        for tableName in Sf.tables.keys():

          print tableName

          for tagName in Sf.tables[tableName].tags.keys():

            print "  ", tagName #,Sf.tables[tableName].tags[tagName]
            
            for error in Sf.tables[tableName].tagErrors[tagName]:
              if error:
                print error

    #nmrStarFile.name = 'local/fulltest.star'
    #nmrStarFile.writeGeneric()

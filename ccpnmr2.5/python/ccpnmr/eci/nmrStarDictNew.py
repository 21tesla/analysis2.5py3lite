#
# Python dictionary with nmrStar information
# Generated on Wed Oct 17 12:16:20 2007
#
# Values per field are [name of tag,function to get data,foreign key name,obligatory value]
#

from ccp.format.nmrStar.util import returnStarInt, returnStarFloat
from ccp.format.nmrStar.util import returnStarString, returnStarDateTime, returnStarAtCode
from ccp.format.nmrStar.util import returnStarCode, returnStarLine, returnStarName, returnStarIdName
from ccp.format.nmrStar.util import returnStarYesNo, returnStarFaxPhoneEmail, returnStarLabel

value = None

sfDict = {

  'entry_interview': {

    'name': 'Entry_interview',
    'saveFrameCode': 'entry_interview',

    'tags': {

      'Sf_category': ['entry_interview',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'PDB_deposition': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'BMRB_deposition': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'View_mode': [None,lambda x = value: returnStarYesNo(x,length = 15),None,False],
      'Structural_genomics': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Ligands': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Non_standard_residues': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Constraints': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Assigned_chem_shifts': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Coupling_constants': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Chem_shift_anisotropy': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Heteronucl_NOEs': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Heteronucl_T1_relaxation': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Heteronucl_T2_relaxation': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Heteronucl_T1rho_relaxation': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Order_parameters': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Residual_dipolar_couplings': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'H_exchange_rate': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'H_exchange_protection_factors': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Spectral_peak_lists': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Dipole_dipole_couplings': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Quadrupolar_couplings': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Homonucl_NOEs': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Dipole_dipole_relaxation': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'DD_cross_correlation': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Dipole_CSA_cross_correlation': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'PKa_value_data_set': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'D_H_fractionation_factors': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Theoretical_chem_shifts': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Spectral_density_values': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Timedomain_data': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Molecular_interactions': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Secondary_structure_orientations': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Metabolite_coordinates': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Mass_spec_data': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Other_kind_of_data': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'PDB_deposition', 'BMRB_deposition', 'View_mode', 'Structural_genomics', 'Ligands', 'Non_standard_residues', 'Constraints', 'Assigned_chem_shifts', 'Coupling_constants', 'Chem_shift_anisotropy', 'Heteronucl_NOEs', 'Heteronucl_T1_relaxation', 'Heteronucl_T2_relaxation', 'Heteronucl_T1rho_relaxation', 'Order_parameters', 'Residual_dipolar_couplings', 'H_exchange_rate', 'H_exchange_protection_factors', 'Spectral_peak_lists', 'Dipole_dipole_couplings', 'Quadrupolar_couplings', 'Homonucl_NOEs', 'Dipole_dipole_relaxation', 'DD_cross_correlation', 'Dipole_CSA_cross_correlation', 'PKa_value_data_set', 'D_H_fractionation_factors', 'Theoretical_chem_shifts', 'Spectral_density_values', 'Timedomain_data', 'Molecular_interactions', 'Secondary_structure_orientations', 'Metabolite_coordinates', 'Mass_spec_data', 'Other_kind_of_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    },

  'deposited_data_files': {

    'name': 'Deposited_data_files',
    'saveFrameCode': 'deposited_data_files',

    'tags': {

      'Sf_category': ['deposited_data_files',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Atomic_coordinate_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Atomic_coordinate_file_syntax': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Constraint_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Constraint_file_syntax': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Precheck_flag': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Validate_flag': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Atomic_coordinate_file_name', 'Atomic_coordinate_file_syntax', 'Constraint_file_name', 'Constraint_file_syntax', 'Precheck_flag', 'Validate_flag'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Upload_data': {

        'tags': {

          'Data_file_ID': [None,returnStarInt,None,True],
          'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Data_file_Sf_category': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
          'Data_file_syntax': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Deposited_data_files_ID': [None,returnStarInt,'Deposited_data_files.ID',True],

                },

        'tagNames': ['Data_file_ID', 'Data_file_name', 'Data_file_Sf_category', 'Data_file_syntax', 'Entry_ID', 'Deposited_data_files_ID'],
        'sourcePrimaryKeys': ['Data_file_ID', 'Entry_ID', 'Deposited_data_files_ID'],

            }

        },

    'tableNames': ['Upload_data']

    },

  'study_list': {

    'name': 'Study_list',

    'tags': {

      'Sf_category': ['study_list',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Study': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Name': [None,returnStarString,None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Study_list_ID': [None,returnStarInt,'Study_list.ID',True],

                },

        'tagNames': ['ID', 'Name', 'Type', 'Details', 'Entry_ID', 'Study_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Study_list_ID'],

            },

      'Study_keyword': {

        'tags': {

          'Study_ID': [None,returnStarInt,'Study.ID',True],
          'Keyword': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Study_list_ID': [None,returnStarInt,'Study_list.ID',True],

                },

        'tagNames': ['Study_ID', 'Keyword', 'Entry_ID', 'Study_list_ID'],
        'sourcePrimaryKeys': ['Study_ID', 'Keyword', 'Entry_ID', 'Study_list_ID'],

            },

      'Study_entry_list': {

        'tags': {

          'Study_ID': [None,returnStarInt,'Study.ID',True],
          'BMRB_accession_code': [None,lambda x = value: returnStarCode(x,length = 12),None,True],
          'BMRB_entry_description': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Study_list_ID': [None,returnStarInt,'Study_list.ID',True],

                },

        'tagNames': ['Study_ID', 'BMRB_accession_code', 'BMRB_entry_description', 'Details', 'Entry_ID', 'Study_list_ID'],
        'sourcePrimaryKeys': ['Study_ID', 'BMRB_accession_code', 'Entry_ID', 'Study_list_ID'],

            }

        },

    'tableNames': ['Study', 'Study_keyword', 'Study_entry_list']

    },

  'entry_information': {

    'name': 'Entry',
    'saveFrameCode': 'entry_information',

    'tags': {

      'Sf_category': ['entry_information',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'ID': [None,lambda x = value: returnStarCode(x,length = 12),None,True],
      'Title': [None,returnStarString,None,True],
      'Version_type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
      'Submission_date': [None,returnStarDateTime,None,True],
      'Accession_date': [None,returnStarDateTime,None,True],
      'Last_release_date': [None,returnStarDateTime,None,False],
      'Original_release_date': [None,returnStarDateTime,None,False],
      'Origination': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
      'NMR_STAR_version': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
      'Original_NMR_STAR_version': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Experimental_method': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Experimental_method_subtype': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Dep_release_code_coordinates': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Dep_release_code_nmr_constraints': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Dep_release_code_nmr_exptl': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Dep_release_code_sequence': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'CASP_target': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Details': [None,returnStarString,None,False],
      'Special_processing_instructions': [None,returnStarString,None,False],
      'Update_BMRB_accession_code': [None,returnStarInt,None,False],
      'Replace_BMRB_accession_code': [None,returnStarInt,None,False],
      'Update_PDB_accession_code': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'Replace_PDB_accession_code': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'BMRB_update_details': [None,returnStarString,None,False],
      'PDB_update_details': [None,returnStarString,None,False],
      'Release_request': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Release_date_request': [None,returnStarDateTime,None,False],
      'Release_date_justification': [None,returnStarString,None,False],
      'Status_code': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
      'Recvd_deposit_form': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Date_deposition_form': [None,returnStarDateTime,None,False],
      'Recvd_coordinates': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Date_coordinates': [None,returnStarDateTime,None,False],
      'Recvd_nmr_constraints': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Date_nmr_constraints': [None,returnStarDateTime,None,False],
      'Recvd_manuscript': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Date_manuscript': [None,returnStarDateTime,None,False],
      'Recvd_author_approval': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Date_author_approval': [None,returnStarDateTime,None,False],
      'Recvd_initial_deposition_date': [None,returnStarDateTime,None,False],
      'PDB_date_submitted': [None,returnStarDateTime,None,False],
      'Author_release_status_code': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
      'Date_of_PDB_release': [None,returnStarDateTime,None,False],
      'Date_hold_coordinates': [None,returnStarDateTime,None,False],
      'Date_hold_nmr_constraints': [None,returnStarDateTime,None,False],
      'PDB_deposit_site': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'PDB_process_site': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'BMRB_deposit_site': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'BMRB_process_site': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'BMRB_annotator': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'BMRB_internal_directory_name': [None,lambda x = value: returnStarLine(x,length = 255),None,False],
      'RCSB_annotator': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Author_approval_type': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Assigned_BMRB_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
      'Assigned_BMRB_deposition_code': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'Assigned_PDB_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
      'Assigned_PDB_deposition_code': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'Assigned_restart_ID': [None,lambda x = value: returnStarLine(x,length = 255),None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'ID', 'Title', 'Version_type', 'Submission_date', 'Accession_date', 'Last_release_date', 'Original_release_date', 'Origination', 'NMR_STAR_version', 'Original_NMR_STAR_version', 'Experimental_method', 'Experimental_method_subtype', 'Dep_release_code_coordinates', 'Dep_release_code_nmr_constraints', 'Dep_release_code_nmr_exptl', 'Dep_release_code_sequence', 'CASP_target', 'Details', 'Special_processing_instructions', 'Update_BMRB_accession_code', 'Replace_BMRB_accession_code', 'Update_PDB_accession_code', 'Replace_PDB_accession_code', 'BMRB_update_details', 'PDB_update_details', 'Release_request', 'Release_date_request', 'Release_date_justification', 'Status_code', 'Recvd_deposit_form', 'Date_deposition_form', 'Recvd_coordinates', 'Date_coordinates', 'Recvd_nmr_constraints', 'Date_nmr_constraints', 'Recvd_manuscript', 'Date_manuscript', 'Recvd_author_approval', 'Date_author_approval', 'Recvd_initial_deposition_date', 'PDB_date_submitted', 'Author_release_status_code', 'Date_of_PDB_release', 'Date_hold_coordinates', 'Date_hold_nmr_constraints', 'PDB_deposit_site', 'PDB_process_site', 'BMRB_deposit_site', 'BMRB_process_site', 'BMRB_annotator', 'BMRB_internal_directory_name', 'RCSB_annotator', 'Author_approval_type', 'Assigned_BMRB_ID', 'Assigned_BMRB_deposition_code', 'Assigned_PDB_ID', 'Assigned_PDB_deposition_code', 'Assigned_restart_ID'],
    'sourcePrimaryKeys': ['ID'],

    'tables': {

      'Entry_proc_cycle': {

        'tags': {

          'Cycle_ID': [None,returnStarInt,None,True],
          'Date_begin_cycle': [None,returnStarDateTime,None,False],
          'Date_end_cycle': [None,returnStarDateTime,None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

                },

        'tagNames': ['Cycle_ID', 'Date_begin_cycle', 'Date_end_cycle', 'Details', 'Entry_ID'],
        'sourcePrimaryKeys': ['Cycle_ID', 'Entry_ID'],

            },

      'Entry_prerelease_seq': {

        'tags': {

          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Entity_label': [None,lambda x = value: returnStarLine(x,length = 127),'Entity.Sf_framecode',False],
          'Seq_one_letter_code': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

                },

        'tagNames': ['Entity_ID', 'Entity_label', 'Seq_one_letter_code', 'Entry_ID'],
        'sourcePrimaryKeys': ['Entity_ID', 'Entry_ID'],

            },

      'Contact_person': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Email_address': [None,lambda x = value: returnStarFaxPhoneEmail(x,length = 127),None,True],
          'Name_salutation': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Given_name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Family_name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Middle_initials': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
          'Family_title': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
          'Department_and_institution': [None,lambda x = value: returnStarString(x,length = 255),None,False],
          'Mailing_address': [None,lambda x = value: returnStarString(x,length = 255),None,False],
          'Address_1': [None,lambda x = value: returnStarString(x,length = 127),None,True],
          'Address_2': [None,lambda x = value: returnStarString(x,length = 127),None,False],
          'Address_3': [None,lambda x = value: returnStarString(x,length = 127),None,False],
          'City': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'State_province': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Country': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Postal_code': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Phone_number': [None,lambda x = value: returnStarFaxPhoneEmail(x,length = 31),None,True],
          'FAX_number': [None,lambda x = value: returnStarFaxPhoneEmail(x,length = 31),None,False],
          'Role': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Organization_type': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

                },

        'tagNames': ['ID', 'Email_address', 'Name_salutation', 'Given_name', 'Family_name', 'Middle_initials', 'Family_title', 'Department_and_institution', 'Mailing_address', 'Address_1', 'Address_2', 'Address_3', 'City', 'State_province', 'Country', 'Postal_code', 'Phone_number', 'FAX_number', 'Role', 'Organization_type', 'Entry_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID'],

            },

      'Entry_author': {

        'tags': {

          'Ordinal': [None,returnStarInt,None,True],
          'Given_name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Family_name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'First_initial': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
          'Middle_initials': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
          'Family_title': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

                },

        'tagNames': ['Ordinal', 'Given_name', 'Family_name', 'First_initial', 'Middle_initials', 'Family_title', 'Entry_ID'],
        'sourcePrimaryKeys': ['Ordinal', 'Entry_ID'],

            },

      'SG_project': {

        'tags': {

          'SG_project_ID': [None,returnStarInt,None,True],
          'Project_name': [None,lambda x = value: returnStarString(x,length = 127),None,True],
          'Full_name_of_center': [None,lambda x = value: returnStarString(x,length = 127),None,False],
          'Initial_of_center': [None,lambda x = value: returnStarString(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

                },

        'tagNames': ['SG_project_ID', 'Project_name', 'Full_name_of_center', 'Initial_of_center', 'Entry_ID'],
        'sourcePrimaryKeys': ['SG_project_ID', 'Entry_ID'],

            },

      'Entry_src': {

        'tags': {

          'Project_name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Organization_full_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Organization_initials': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

                },

        'tagNames': ['Project_name', 'Organization_full_name', 'Organization_initials', 'Entry_ID'],
        'sourcePrimaryKeys': ['Project_name', 'Organization_initials', 'Entry_ID'],

            },

      'Struct_keywords': {

        'tags': {

          'Keywords': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Text': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

                },

        'tagNames': ['Keywords', 'Text', 'Entry_ID'],
        'sourcePrimaryKeys': ['Keywords', 'Entry_ID'],

            },

      'Data_set': {

        'tags': {

          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Count': [None,returnStarInt,None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

                },

        'tagNames': ['Type', 'Count', 'Entry_ID'],
        'sourcePrimaryKeys': ['Type', 'Entry_ID'],

            },

      'Datum': {

        'tags': {

          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Count': [None,returnStarInt,None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

                },

        'tagNames': ['Type', 'Count', 'Entry_ID'],
        'sourcePrimaryKeys': ['Type', 'Entry_ID'],

            },

      'Release': {

        'tags': {

          'Release_number': [None,returnStarInt,None,True],
          'Date': [None,returnStarDateTime,None,True],
          'Submission_date': [None,returnStarDateTime,None,False],
          'Type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Author': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Detail': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

                },

        'tagNames': ['Release_number', 'Date', 'Submission_date', 'Type', 'Author', 'Detail', 'Entry_ID'],
        'sourcePrimaryKeys': ['Release_number', 'Entry_ID'],

            },

      'Related_entries': {

        'tags': {

          'Database_name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Database_accession_code': [None,lambda x = value: returnStarLine(x,length = 12),None,True],
          'Relationship': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

                },

        'tagNames': ['Database_name', 'Database_accession_code', 'Relationship', 'Entry_ID'],
        'sourcePrimaryKeys': ['Database_name', 'Database_accession_code', 'Entry_ID'],

            }

        },

    'tableNames': ['Entry_proc_cycle', 'Entry_prerelease_seq', 'Contact_person', 'Entry_author', 'SG_project', 'Entry_src', 'Struct_keywords', 'Data_set', 'Datum', 'Release', 'Related_entries']

    },

  'citations': {

    'name': 'Citation',

    'tags': {

      'Sf_category': ['citations',lambda x = value: returnStarLine(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Class': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
      'CAS_abstract_code': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'MEDLINE_UI_code': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'PubMed_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
      'Full_citation': [None,returnStarString,None,False],
      'Title': [None,returnStarString,None,False],
      'Status': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Type': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Journal_abbrev': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Journal_name_full': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Journal_volume': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
      'Journal_issue': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
      'Journal_ASTM': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Journal_ISSN': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Journal_CSD': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Book_title': [None,lambda x = value: returnStarString(x,length = 255),None,False],
      'Book_chapter_title': [None,lambda x = value: returnStarString(x,length = 255),None,False],
      'Book_volume': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
      'Book_series': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
      'Book_publisher': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Book_publisher_city': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Book_ISBN': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Conference_title': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Conference_site': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Conference_state_province': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Conference_country': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Conference_start_date': [None,returnStarDateTime,None,False],
      'Conference_end_date': [None,returnStarDateTime,None,False],
      'Conference_abstract_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Thesis_institution': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Thesis_institution_city': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Thesis_institution_country': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'WWW_URL': [None,returnStarString,None,False],
      'Page_first': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
      'Page_last': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
      'Year': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Class', 'CAS_abstract_code', 'MEDLINE_UI_code', 'PubMed_ID', 'Full_citation', 'Title', 'Status', 'Type', 'Journal_abbrev', 'Journal_name_full', 'Journal_volume', 'Journal_issue', 'Journal_ASTM', 'Journal_ISSN', 'Journal_CSD', 'Book_title', 'Book_chapter_title', 'Book_volume', 'Book_series', 'Book_publisher', 'Book_publisher_city', 'Book_ISBN', 'Conference_title', 'Conference_site', 'Conference_state_province', 'Conference_country', 'Conference_start_date', 'Conference_end_date', 'Conference_abstract_number', 'Thesis_institution', 'Thesis_institution_city', 'Thesis_institution_country', 'WWW_URL', 'Page_first', 'Page_last', 'Year', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Citation_author': {

        'tags': {

          'Ordinal': [None,returnStarInt,None,True],
          'Given_name': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Family_name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'First_initial': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
          'Middle_initials': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
          'Family_title': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Citation_ID': [None,returnStarInt,'Citation.ID',True],

                },

        'tagNames': ['Ordinal', 'Given_name', 'Family_name', 'First_initial', 'Middle_initials', 'Family_title', 'Entry_ID', 'Citation_ID'],
        'sourcePrimaryKeys': ['Ordinal', 'Entry_ID', 'Citation_ID'],

            },

      'Citation_keyword': {

        'tags': {

          'Keyword': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Citation_ID': [None,returnStarInt,'Citation.ID',True],

                },

        'tagNames': ['Keyword', 'Entry_ID', 'Citation_ID'],
        'sourcePrimaryKeys': ['Keyword', 'Entry_ID', 'Citation_ID'],

            },

      'Citation_editor': {

        'tags': {

          'Ordinal': [None,returnStarInt,None,True],
          'Given_name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Family_name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'First_initial': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
          'Middle_initials': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
          'Family_title': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Citation_ID': [None,returnStarInt,'Citation.ID',True],

                },

        'tagNames': ['Ordinal', 'Given_name', 'Family_name', 'First_initial', 'Middle_initials', 'Family_title', 'Entry_ID', 'Citation_ID'],
        'sourcePrimaryKeys': ['Ordinal', 'Entry_ID', 'Citation_ID'],

            }

        },

    'tableNames': ['Citation_author', 'Citation_keyword', 'Citation_editor']

    },

  'assembly': {

    'name': 'Assembly',

    'tags': {

      'Sf_category': ['assembly',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'BMRB_code': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
      'Number_of_components': [None,returnStarInt,None,False],
      'Organic_ligands': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Metal_ions': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Non_standard_bonds': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Ambiguous_conformational_states': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Ambiguous_chem_comp_sites': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Molecules_in_chemical_exchange': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Paramagnetic': [None,lambda x = value: returnStarYesNo(x,length = 3),None,True],
      'Thiol_state': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Molecular_mass': [None,returnStarFloat,None,False],
      'Enzyme_commission_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],
      'DB_query_date': [None,returnStarDateTime,None,False],
      'DB_query_revised_last_date': [None,returnStarDateTime,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Name', 'BMRB_code', 'Number_of_components', 'Organic_ligands', 'Metal_ions', 'Non_standard_bonds', 'Ambiguous_conformational_states', 'Ambiguous_chem_comp_sites', 'Molecules_in_chemical_exchange', 'Paramagnetic', 'Thiol_state', 'Molecular_mass', 'Enzyme_commission_number', 'Details', 'DB_query_date', 'DB_query_revised_last_date'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Assembly_type': {

        'tags': {

          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['Type', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['Type', 'Entry_ID', 'Assembly_ID'],

            },

      'Entity_assembly': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Entity_assembly_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Entity_label': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',False],
          'Asym_ID': [None,lambda x = value: returnStarCode(x,length = 31),'Struct_asym.ID',False],
          'Experimental_data_reported': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
          'Physical_state': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Conformational_isomer': [None,lambda x = value: returnStarYesNo(x,length = 3),None,True],
          'Chemical_exchange_state': [None,lambda x = value: returnStarYesNo(x,length = 3),None,True],
          'Magnetic_equivalence_group_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Role': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['ID', 'Entity_assembly_name', 'Entity_ID', 'Entity_label', 'Asym_ID', 'Experimental_data_reported', 'Physical_state', 'Conformational_isomer', 'Chemical_exchange_state', 'Magnetic_equivalence_group_code', 'Role', 'Details', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['ID', 'Entity_assembly_name', 'Entry_ID', 'Assembly_ID'],

            },

      'Bond': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Order': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_assembly_name_1': [None,lambda x = value: returnStarLine(x,length = 127),'Entity_assembly.Entity_assembly_name',False],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_assembly_name_2': [None,lambda x = value: returnStarLine(x,length = 127),'Entity_assembly.Entity_assembly_name',False],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Auth_entity_assembly_name_1': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Auth_entity_assembly_name_2': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['ID', 'Type', 'Order', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_assembly_name_1', 'Entity_ID_1', 'Comp_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Atom_ID_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_assembly_name_2', 'Entity_ID_2', 'Comp_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Atom_ID_2', 'Auth_entity_assembly_ID_1', 'Auth_entity_assembly_name_1', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_entity_assembly_ID_2', 'Auth_entity_assembly_name_2', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Assembly_ID'],

            },

      'Deleted_atom': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_assembly_name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Entity_label': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_assembly_name', 'Entity_ID', 'Entity_label', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Comp_label', 'Atom_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Atom_type', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Assembly_ID'],

            },

      'Struct_asym': {

        'tags': {

          'ID': [None,lambda x = value: returnStarCode(x,length = 12),None,True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'PDBX_blank_PDB_chainid_flag': [None,lambda x = value: returnStarLine(x,length = 3),None,False],
          'PDBX_modified': [None,lambda x = value: returnStarLine(x,length = 3),None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['ID', 'Entity_ID', 'PDBX_blank_PDB_chainid_flag', 'PDBX_modified', 'Details', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Assembly_ID'],

            },

      'Assembly_db_link': {

        'tags': {

          'Author_supplied': [None,lambda x = value: returnStarYesNo(x,length = 12),None,True],
          'Database_code': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Accession_code': [None,lambda x = value: returnStarLine(x,length = 15),None,True],
          'Entry_mol_code': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_mol_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_experimental_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_structure_resolution': [None,returnStarFloat,None,False],
          'Entry_relation_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['Author_supplied', 'Database_code', 'Accession_code', 'Entry_mol_code', 'Entry_mol_name', 'Entry_experimental_method', 'Entry_structure_resolution', 'Entry_relation_type', 'Entry_details', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['Database_code', 'Accession_code', 'Entry_ID', 'Assembly_ID'],

            },

      'Assembly_common_name': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['Name', 'Type', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['Name', 'Entry_ID', 'Assembly_ID'],

            },

      'Assembly_systematic_name': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Naming_system': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['Name', 'Naming_system', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['Naming_system', 'Entry_ID', 'Assembly_ID'],

            },

      'Assembly_interaction': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Mol_interaction_type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['ID', 'Entity_assembly_ID_1', 'Entity_assembly_ID_2', 'Mol_interaction_type', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Assembly_ID'],

            },

      'Chem_comp_assembly': {

        'tags': {

          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',False],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'Assembly_chem_comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Comp_ID', 'Assembly_chem_comp_ID', 'Seq_ID', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['Entity_assembly_ID', 'Comp_index_ID', 'Entry_ID', 'Assembly_ID'],

            },

      'PDBX_poly_seq_scheme': {

        'tags': {

          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Asym_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Atom_site.Label_asym_ID',False],
          'Entity_ID': [None,returnStarInt,'Entity.ID',False],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Seq_ID': [None,returnStarInt,'Entity_poly_seq.Num',True],
          'Mon_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entity_poly_seq.Mon_ID',False],
          'PDB_seq_num': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_num': [None,lambda x = value: returnStarCode(x,length = 12),'Atom_site.Auth_seq_ID',False],
          'PDB_mon_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_mon_ID': [None,lambda x = value: returnStarCode(x,length = 15),'Atom_site.Auth_comp_ID',False],
          'PDB_strand_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'PDB_ins_code': [None,lambda x = value: returnStarLine(x,length = 12),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['Entity_assembly_ID', 'Asym_ID', 'Entity_ID', 'Comp_index_ID', 'Comp_ID', 'Comp_label', 'Seq_ID', 'Mon_ID', 'PDB_seq_num', 'Auth_seq_num', 'PDB_mon_ID', 'Auth_mon_ID', 'PDB_strand_ID', 'PDB_ins_code', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['Entity_assembly_ID', 'Comp_index_ID', 'Seq_ID', 'Entry_ID', 'Assembly_ID'],

            },

      'PDBX_nonpoly_scheme': {

        'tags': {

          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Asym_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Atom_site.Label_asym_ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Mon_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Atom_site.Label_comp_ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'PDB_seq_num': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_num': [None,lambda x = value: returnStarCode(x,length = 12),'Atom_site.Auth_seq_ID',False],
          'PDB_mon_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_mon_ID': [None,lambda x = value: returnStarCode(x,length = 15),'Atom_site.Auth_comp_ID',False],
          'PDB_strand_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'PDB_ins_code': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['Entity_assembly_ID', 'Asym_ID', 'Entity_ID', 'Mon_ID', 'Comp_index_ID', 'Comp_ID', 'PDB_seq_num', 'Auth_seq_num', 'PDB_mon_ID', 'Auth_mon_ID', 'PDB_strand_ID', 'PDB_ins_code', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['Entity_assembly_ID', 'Entry_ID', 'Assembly_ID'],

            },

      'Atom_type': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Symbol': [None,lambda x = value: returnStarCode(x,length = 3),None,False],
          'Atomic_number': [None,returnStarInt,None,False],
          'Isotope_number': [None,returnStarInt,None,False],
          'Oxidation_number': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Paramagnetic': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
          'Electron_configuration': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Unpaired_electron_number': [None,returnStarInt,None,False],
          'Atomic_mass': [None,returnStarFloat,None,False],
          'Van_der_Vaals_radii': [None,returnStarFloat,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['ID', 'Symbol', 'Atomic_number', 'Isotope_number', 'Oxidation_number', 'Paramagnetic', 'Electron_configuration', 'Unpaired_electron_number', 'Atomic_mass', 'Van_der_Vaals_radii', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Assembly_ID'],

            },

      'Atom': {

        'tags': {

          'Assembly_atom_ID': [None,returnStarInt,None,True],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_type_ID': [None,returnStarInt,'Atom_type.ID',False],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),None,True],
          'Type_symbol': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'PDB_one_letter_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'PDB_strand_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'PDB_ins_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'PDB_asym_ID': [None,lambda x = value: returnStarCode(x,length = 15),'Atom_site.Label_asym_ID',False],
          'PDB_seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'PDB_comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'PDB_group': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'PDB_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'PDB_atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_type_ID', 'Atom_ID', 'Type_symbol', 'PDB_one_letter_code', 'PDB_strand_ID', 'PDB_ins_code', 'PDB_asym_ID', 'PDB_seq_ID', 'PDB_comp_ID', 'PDB_group', 'PDB_atom_ID', 'PDB_atom_type', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['Assembly_atom_ID', 'Entry_ID', 'Assembly_ID'],

            },

      'Assembly_bio_function': {

        'tags': {

          'Biological_function': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['Biological_function', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['Biological_function', 'Entry_ID', 'Assembly_ID'],

            },

      'Angle': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Angle_name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Entity_label_1': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_label_1': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Entity_label_2': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_label_2': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Assembly_atom_ID_3': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_3': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_3': [None,returnStarInt,'Entity.ID',True],
          'Entity_label_3': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Comp_ID_3': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_label_3': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Comp_index_ID_3': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_3': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atom_ID_3': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_3': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['ID', 'Angle_name', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Entity_label_1', 'Comp_ID_1', 'Comp_label_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Entity_label_2', 'Comp_ID_2', 'Comp_label_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Assembly_atom_ID_3', 'Entity_assembly_ID_3', 'Entity_ID_3', 'Entity_label_3', 'Comp_ID_3', 'Comp_label_3', 'Comp_index_ID_3', 'Seq_ID_3', 'Atom_ID_3', 'Atom_type_3', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Assembly_ID'],

            },

      'Torsion_angle': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Torsion_angle_name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Entity_label_1': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_label_1': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Entity_label_2': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_label_2': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Assembly_atom_ID_3': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_3': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_3': [None,returnStarInt,'Entity.ID',True],
          'Entity_label_3': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Comp_ID_3': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_label_3': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Comp_index_ID_3': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_3': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atom_ID_3': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_3': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Assembly_atom_ID_4': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_4': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_4': [None,returnStarInt,'Entity.ID',True],
          'Entity_label_4': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Comp_ID_4': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_label_4': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Comp_index_ID_4': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_4': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atom_ID_4': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_4': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['ID', 'Torsion_angle_name', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Entity_label_1', 'Comp_ID_1', 'Comp_label_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Entity_label_2', 'Comp_ID_2', 'Comp_label_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Assembly_atom_ID_3', 'Entity_assembly_ID_3', 'Entity_ID_3', 'Entity_label_3', 'Comp_ID_3', 'Comp_label_3', 'Comp_index_ID_3', 'Seq_ID_3', 'Atom_ID_3', 'Atom_type_3', 'Assembly_atom_ID_4', 'Entity_assembly_ID_4', 'Entity_ID_4', 'Entity_label_4', 'Comp_ID_4', 'Comp_label_4', 'Comp_index_ID_4', 'Seq_ID_4', 'Atom_ID_4', 'Atom_type_4', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Assembly_ID'],

            },

      'Assembly_segment': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Entity_label': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',False],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['ID', 'Entity_assembly_ID', 'Entity_ID', 'Entity_label', 'Comp_index_ID', 'Comp_ID', 'Comp_label', 'Seq_ID', 'Atom_ID', 'Assembly_atom_ID', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['ID'],

            },

      'Assembly_segment_description': {

        'tags': {

          'Assembly_segment_ID': [None,returnStarInt,'Assembly_segment.ID',True],
          'Code': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['Assembly_segment_ID', 'Code', 'Details', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['Entry_ID', 'Assembly_ID'],

            },

      'Assembly_keyword': {

        'tags': {

          'Keyword': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['Keyword', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['Keyword', 'Entry_ID', 'Assembly_ID'],

            },

      'Assembly_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'Assembly_ID'],
        'sourcePrimaryKeys': ['Citation_ID', 'Entry_ID', 'Assembly_ID'],

            }

        },

    'tableNames': ['Assembly_type', 'Entity_assembly', 'Bond', 'Deleted_atom', 'Struct_asym', 'Assembly_db_link', 'Assembly_common_name', 'Assembly_systematic_name', 'Assembly_interaction', 'Chem_comp_assembly', 'PDBX_poly_seq_scheme', 'PDBX_nonpoly_scheme', 'Atom_type', 'Atom', 'Assembly_bio_function', 'Angle', 'Torsion_angle', 'Assembly_segment', 'Assembly_segment_description', 'Assembly_keyword', 'Assembly_citation']

    },

  'assembly_annotation': {

    'name': 'Assembly_annotation_list',

    'tags': {

      'Sf_category': ['assembly_annotation',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Source': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Source', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Author_annotation': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',True],
          'Assembly_subsystem_ID': [None,returnStarInt,'Assembly_subsystem.ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Comp_index_ID_start': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Comp_index_ID_end': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Seq_ID_start': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Seq_ID_end': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Annotation_code': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_annotation_list_ID': [None,returnStarInt,'Assembly_annotation_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_ID', 'Assembly_subsystem_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Comp_index_ID_start', 'Comp_index_ID_end', 'Seq_ID', 'Seq_ID_start', 'Seq_ID_end', 'Comp_ID', 'Atom_ID', 'Assembly_atom_ID', 'Annotation_code', 'Entry_ID', 'Assembly_annotation_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Assembly_annotation_list_ID'],

            }

        },

    'tableNames': ['Author_annotation']

    },

  'assembly_subsystems': {

    'name': 'Assembly_subsystem',

    'tags': {

      'Sf_category': ['assembly_subsystems',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'CAS_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'CAS_registry_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Enzyme_commission_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Molecular_mass': [None,returnStarFloat,None,False],
      'Details': [None,returnStarString,None,False],
      'DB_query_date': [None,returnStarDateTime,None,False],
      'DB_last_query_revised_last_date': [None,returnStarDateTime,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Name', 'CAS_name', 'CAS_registry_number', 'Enzyme_commission_number', 'Molecular_mass', 'Details', 'DB_query_date', 'DB_last_query_revised_last_date'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Subsystem_common_name': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_subsystem_ID': [None,returnStarInt,'Assembly_subsystem.ID',True],

                },

        'tagNames': ['Name', 'Entry_ID', 'Assembly_subsystem_ID'],
        'sourcePrimaryKeys': ['Name', 'Entry_ID', 'Assembly_subsystem_ID'],

            },

      'Subsystem_type': {

        'tags': {

          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_subsystem_ID': [None,returnStarInt,'Assembly_subsystem.ID',True],

                },

        'tagNames': ['Type', 'Entry_ID', 'Assembly_subsystem_ID'],
        'sourcePrimaryKeys': ['Type', 'Entry_ID', 'Assembly_subsystem_ID'],

            },

      'Subsystem_component': {

        'tags': {

          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Entity_label': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_subsystem_ID': [None,returnStarInt,'Assembly_subsystem.ID',True],

                },

        'tagNames': ['Entity_assembly_ID', 'Entity_ID', 'Entity_label', 'Entry_ID', 'Assembly_subsystem_ID'],
        'sourcePrimaryKeys': ['Entity_assembly_ID', 'Entry_ID', 'Assembly_subsystem_ID'],

            },

      'Subsystem_keyword': {

        'tags': {

          'Keyword': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_subsystem_ID': [None,returnStarInt,'Assembly_subsystem.ID',True],

                },

        'tagNames': ['Keyword', 'Entry_ID', 'Assembly_subsystem_ID'],
        'sourcePrimaryKeys': ['Keyword', 'Entry_ID', 'Assembly_subsystem_ID'],

            },

      'Subsystem_biological_function': {

        'tags': {

          'Biological_function': [None,lambda x = value: returnStarString(x,length = 255),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_subsystem_ID': [None,returnStarInt,'Assembly_subsystem.ID',True],

                },

        'tagNames': ['Biological_function', 'Entry_ID', 'Assembly_subsystem_ID'],
        'sourcePrimaryKeys': ['Biological_function', 'Entry_ID', 'Assembly_subsystem_ID'],

            },

      'Subsystem_db_link': {

        'tags': {

          'Author_supplied': [None,lambda x = value: returnStarYesNo(x,length = 12),None,False],
          'Database_code': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Accession_code': [None,lambda x = value: returnStarLine(x,length = 15),None,True],
          'Entry_mol_code': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_mol_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_experimental_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_relation_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_subsystem_ID': [None,returnStarInt,'Assembly_subsystem.ID',False],

                },

        'tagNames': ['Author_supplied', 'Database_code', 'Accession_code', 'Entry_mol_code', 'Entry_mol_name', 'Entry_experimental_method', 'Entry_relation_type', 'Entry_details', 'Entry_ID', 'Assembly_subsystem_ID'],
        'sourcePrimaryKeys': ['Database_code', 'Accession_code'],

            },

      'Subsystem_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assembly_subsystem_list_ID': [None,returnStarInt,'Assembly_subsystem.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'Assembly_subsystem_list_ID'],
        'sourcePrimaryKeys': ['Citation_ID', 'Entry_ID', 'Assembly_subsystem_list_ID'],

            }

        },

    'tableNames': ['Subsystem_common_name', 'Subsystem_type', 'Subsystem_component', 'Subsystem_keyword', 'Subsystem_biological_function', 'Subsystem_db_link', 'Subsystem_citation']

    },

  'entity': {

    'name': 'Entity',

    'tags': {

      'Sf_category': ['entity',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'BMRB_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
      'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Type': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Polymer_common_type': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Polymer_type': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Polymer_type_details': [None,returnStarString,None,False],
      'Polymer_strand_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
      'Polymer_seq_one_letter_code_can': [None,returnStarString,None,False],
      'Polymer_seq_one_letter_code': [None,returnStarString,None,False],
      'Target_identifier': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Polymer_author_defined_seq': [None,returnStarString,None,False],
      'Polymer_author_seq_details': [None,returnStarString,None,False],
      'Ambiguous_conformational_states': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Ambiguous_chem_comp_sites': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Nstd_monomer': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Nstd_chirality': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Nstd_linkage': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Nonpolymer_comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
      'Nonpolymer_comp_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Number_of_monomers': [None,returnStarInt,None,False],
      'Number_of_nonpolymer_components': [None,returnStarInt,None,False],
      'Paramagnetic': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Thiol_state': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Src_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Parent_entity_ID': [None,returnStarInt,'Entity.ID',False],
      'Fragment': [None,returnStarString,None,False],
      'Mutation': [None,returnStarString,None,False],
      'EC_number': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Calc_isoelectric_point': [None,returnStarFloat,None,False],
      'Formula_weight': [None,returnStarFloat,None,False],
      'Formula_weight_exptl': [None,returnStarFloat,None,False],
      'Formula_weight_exptl_meth': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],
      'DB_query_date': [None,returnStarDateTime,None,False],
      'DB_query_revised_last_date': [None,returnStarDateTime,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'BMRB_code', 'Name', 'Type', 'Polymer_common_type', 'Polymer_type', 'Polymer_type_details', 'Polymer_strand_ID', 'Polymer_seq_one_letter_code_can', 'Polymer_seq_one_letter_code', 'Target_identifier', 'Polymer_author_defined_seq', 'Polymer_author_seq_details', 'Ambiguous_conformational_states', 'Ambiguous_chem_comp_sites', 'Nstd_monomer', 'Nstd_chirality', 'Nstd_linkage', 'Nonpolymer_comp_ID', 'Nonpolymer_comp_label', 'Number_of_monomers', 'Number_of_nonpolymer_components', 'Paramagnetic', 'Thiol_state', 'Src_method', 'Parent_entity_ID', 'Fragment', 'Mutation', 'EC_number', 'Calc_isoelectric_point', 'Formula_weight', 'Formula_weight_exptl', 'Formula_weight_exptl_meth', 'Details', 'DB_query_date', 'DB_query_revised_last_date'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Entity_db_link': {

        'tags': {

          'Author_supplied': [None,lambda x = value: returnStarYesNo(x,length = 12),None,False],
          'Database_code': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Accession_code': [None,lambda x = value: returnStarLine(x,length = 15),None,True],
          'Entry_mol_code': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_mol_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_experimental_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_structure_resolution': [None,returnStarFloat,None,False],
          'Entry_relation_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_details': [None,returnStarString,None,False],
          'Chimera_segment': [None,returnStarInt,None,False],
          'Seq_query_to_submitted_percent': [None,returnStarFloat,None,False],
          'Seq_subject_length': [None,returnStarInt,None,False],
          'Seq_identity': [None,returnStarFloat,None,False],
          'Seq_positive': [None,returnStarFloat,None,False],
          'Seq_homology_expectation_val': [None,returnStarFloat,None,False],
          'Seq_align_begin': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Seq_align_end': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Seq_difference_details': [None,returnStarString,None,False],
          'Seq_alignment_details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],

                },

        'tagNames': ['Author_supplied', 'Database_code', 'Accession_code', 'Entry_mol_code', 'Entry_mol_name', 'Entry_experimental_method', 'Entry_structure_resolution', 'Entry_relation_type', 'Entry_details', 'Chimera_segment', 'Seq_query_to_submitted_percent', 'Seq_subject_length', 'Seq_identity', 'Seq_positive', 'Seq_homology_expectation_val', 'Seq_align_begin', 'Seq_align_end', 'Seq_difference_details', 'Seq_alignment_details', 'Entry_ID', 'Entity_ID'],
        'sourcePrimaryKeys': ['Database_code', 'Accession_code', 'Entry_ID', 'Entity_ID'],

            },

      'Entity_biological_function': {

        'tags': {

          'Biological_function': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],

                },

        'tagNames': ['Biological_function', 'Entry_ID', 'Entity_ID'],
        'sourcePrimaryKeys': ['Biological_function', 'Entry_ID', 'Entity_ID'],

            },

      'Entity_common_name': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],

                },

        'tagNames': ['Name', 'Type', 'Entry_ID', 'Entity_ID'],
        'sourcePrimaryKeys': ['Name', 'Entry_ID', 'Entity_ID'],

            },

      'Entity_systematic_name': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Naming_system': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],

                },

        'tagNames': ['Name', 'Naming_system', 'Entry_ID', 'Entity_ID'],
        'sourcePrimaryKeys': ['Naming_system', 'Entry_ID', 'Entity_ID'],

            },

      'Entity_keyword': {

        'tags': {

          'Keyword': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],

                },

        'tagNames': ['Keyword', 'Entry_ID', 'Entity_ID'],
        'sourcePrimaryKeys': ['Keyword', 'Entry_ID', 'Entity_ID'],

            },

      'Entity_comp_index': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],

                },

        'tagNames': ['ID', 'Auth_seq_ID', 'Comp_ID', 'Comp_label', 'Entry_ID', 'Entity_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Entity_ID'],

            },

      'Entity_poly_seq': {

        'tags': {

          'Hetero': [None,lambda x = value: returnStarCode(x,length = 3),None,False],
          'Mon_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Num': [None,returnStarInt,None,True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],

                },

        'tagNames': ['Hetero', 'Mon_ID', 'Num', 'Comp_index_ID', 'Entry_ID', 'Entity_ID'],
        'sourcePrimaryKeys': ['Mon_ID', 'Num', 'Entry_ID', 'Entity_ID'],

            },

      'Entity_chimera_segment': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Comp_index_ID_begin': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Comp_index_ID_end': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_begin': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Seq_ID_end': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],

                },

        'tagNames': ['ID', 'Comp_index_ID_begin', 'Comp_index_ID_end', 'Seq_ID_begin', 'Seq_ID_end', 'Details', 'Entry_ID', 'Entity_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Entity_ID'],

            },

      'Entity_comp_index_alt': {

        'tags': {

          'Entity_comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],

                },

        'tagNames': ['Entity_comp_index_ID', 'Auth_seq_ID', 'Comp_ID', 'Comp_label', 'Entry_ID', 'Entity_ID'],
        'sourcePrimaryKeys': ['Entity_comp_index_ID', 'Comp_ID', 'Entry_ID', 'Entity_ID'],

            },

      'Entity_bond': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Value_order': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],

                },

        'tagNames': ['ID', 'Type', 'Value_order', 'Comp_index_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Comp_index_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Details', 'Entry_ID', 'Entity_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Entity_ID'],

            },

      'Entity_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'Entity_ID'],
        'sourcePrimaryKeys': ['Citation_ID', 'Entry_ID', 'Entity_ID'],

            }

        },

    'tableNames': ['Entity_db_link', 'Entity_biological_function', 'Entity_common_name', 'Entity_systematic_name', 'Entity_keyword', 'Entity_comp_index', 'Entity_poly_seq', 'Entity_chimera_segment', 'Entity_comp_index_alt', 'Entity_bond', 'Entity_citation']

    },

  'natural_source': {

    'name': 'Entity_natural_src_list',
    'saveFrameCode': 'natural_source',

    'tags': {

      'Sf_category': ['natural_source',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Entity_natural_src': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Entity_label': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Entity_chimera_segment_ID': [None,returnStarInt,'Entity_chimera_segment.ID',False],
          'NCBI_taxonomy_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Common': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
          'Organism_name_scientific': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Organism_name_common': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Organism_acronym': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'ICTVdb_decimal_code': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Superkingdom': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Kingdom': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Genus': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Species': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Strain': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Variant': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Subvariant': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Organ': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Tissue': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Tissue_fraction': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Cell_line': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Cell_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'ATCC_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Organelle': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Cellular_location': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Fragment': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Fraction': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Secretion': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Plasmid': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Plasmid_details': [None,returnStarString,None,False],
          'Gene_mnemonic': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Dev_stage': [None,returnStarString,None,False],
          'Details': [None,returnStarString,None,False],
          'Citation_ID': [None,returnStarInt,'Citation.ID',False],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_natural_src_list_ID': [None,returnStarInt,'Entity_natural_src_list.ID',True],

                },

        'tagNames': ['ID', 'Entity_ID', 'Entity_label', 'Entity_chimera_segment_ID', 'NCBI_taxonomy_ID', 'Type', 'Common', 'Organism_name_scientific', 'Organism_name_common', 'Organism_acronym', 'ICTVdb_decimal_code', 'Superkingdom', 'Kingdom', 'Genus', 'Species', 'Strain', 'Variant', 'Subvariant', 'Organ', 'Tissue', 'Tissue_fraction', 'Cell_line', 'Cell_type', 'ATCC_number', 'Organelle', 'Cellular_location', 'Fragment', 'Fraction', 'Secretion', 'Plasmid', 'Plasmid_details', 'Gene_mnemonic', 'Dev_stage', 'Details', 'Citation_ID', 'Citation_label', 'Entry_ID', 'Entity_natural_src_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Entity_natural_src_list_ID'],

            },

      'Natural_source_db': {

        'tags': {

          'Entity_natural_src_ID': [None,returnStarInt,'Entity_natural_src.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Entity_label': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Entity_chimera_segment_ID': [None,returnStarInt,'Entity_chimera_segment.ID',False],
          'Database_code': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Database_type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_code': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'ORF_code': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Gene_locus_code': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Gene_cDNA_code': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_natural_src_list_ID': [None,returnStarInt,'Entity_natural_src_list.ID',True],

                },

        'tagNames': ['Entity_natural_src_ID', 'Entity_ID', 'Entity_label', 'Entity_chimera_segment_ID', 'Database_code', 'Database_type', 'Entry_code', 'Entry_type', 'ORF_code', 'Gene_locus_code', 'Gene_cDNA_code', 'Entry_ID', 'Entity_natural_src_list_ID'],
        'sourcePrimaryKeys': ['Entity_natural_src_ID', 'Database_code', 'Entry_code', 'Entry_ID', 'Entity_natural_src_list_ID'],

            }

        },

    'tableNames': ['Entity_natural_src', 'Natural_source_db']

    },

  'experimental_source': {

    'name': 'Entity_experimental_src_list',
    'saveFrameCode': 'experimental_source',

    'tags': {

      'Sf_category': ['experimental_source',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Entity_experimental_src': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Entity_label': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Entity_chimera_segment_ID': [None,returnStarInt,'Entity_chimera_segment.ID',False],
          'Production_method': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Host_org_scientific_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_name_common': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_details': [None,returnStarString,None,False],
          'Host_org_NCBI_taxonomy_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Host_org_genus': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_species': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_strain': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_variant': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_subvariant': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_organ': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_tissue': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_tissue_fraction': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_cell_line': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_cell_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_cellular_location': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_organelle': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_gene': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_culture_collection': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_ATCC_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Vector_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'PDBview_host_org_vector_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'PDBview_plasmid_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Vector_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Vector_details': [None,returnStarString,None,False],
          'Vendor_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Host_org_dev_stage': [None,returnStarString,None,False],
          'Details': [None,returnStarString,None,False],
          'Citation_ID': [None,returnStarInt,'Citation.ID',False],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_experimental_src_list_ID': [None,returnStarInt,'Entity_experimental_src_list.ID',True],

                },

        'tagNames': ['ID', 'Entity_ID', 'Entity_label', 'Entity_chimera_segment_ID', 'Production_method', 'Host_org_scientific_name', 'Host_org_name_common', 'Host_org_details', 'Host_org_NCBI_taxonomy_ID', 'Host_org_genus', 'Host_org_species', 'Host_org_strain', 'Host_org_variant', 'Host_org_subvariant', 'Host_org_organ', 'Host_org_tissue', 'Host_org_tissue_fraction', 'Host_org_cell_line', 'Host_org_cell_type', 'Host_org_cellular_location', 'Host_org_organelle', 'Host_org_gene', 'Host_org_culture_collection', 'Host_org_ATCC_number', 'Vector_type', 'PDBview_host_org_vector_name', 'PDBview_plasmid_name', 'Vector_name', 'Vector_details', 'Vendor_name', 'Host_org_dev_stage', 'Details', 'Citation_ID', 'Citation_label', 'Entry_ID', 'Entity_experimental_src_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Entity_experimental_src_list_ID'],

            }

        },

    'tableNames': ['Entity_experimental_src']

    },

  'chem_comp': {

    'name': 'Chem_comp',

    'tags': {

      'Sf_category': ['chem_comp',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,lambda x = value: returnStarCode(x,length = 12),None,True],
      'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
      'BMRB_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
      'PDB_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
      'PDBx_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'PDBx_ambiguous_flag': [None,lambda x = value: returnStarCode(x,length = 3),None,False],
      'PubChem_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
      'InCHi_code': [None,lambda x = value: returnStarString(x,length = 1024),None,False],
      'PDB_NSTD_flag': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Std_deriv_one_letter_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
      'Std_deriv_three_letter_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
      'Std_deriv_BMRB_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
      'Std_deriv_PDB_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
      'Std_deriv_chem_comp_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Formal_charge': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
      'Paramagnetic': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Aromatic': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'Empirical_formula': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Formula_weight': [None,returnStarFloat,None,False],
      'Formula_mono_iso_wt_nat': [None,returnStarFloat,None,False],
      'Formula_mono_iso_wt_13C': [None,returnStarFloat,None,False],
      'Formula_mono_iso_wt_15N': [None,returnStarFloat,None,False],
      'Formula_mono_iso_wt_13C_15N': [None,returnStarFloat,None,False],
      'Image_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Image_file_format': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Topo_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Topo_file_format': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Struct_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Struct_file_format': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Stereochem_param_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Stereochem_param_file_format': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Vendor': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Vendor_product_code': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],
      'DB_query_date': [None,returnStarDateTime,None,False],
      'DB_last_query_revised_last_date': [None,returnStarDateTime,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Name', 'Type', 'BMRB_code', 'PDB_code', 'PDBx_type', 'PDBx_ambiguous_flag', 'PubChem_code', 'InCHi_code', 'PDB_NSTD_flag', 'Std_deriv_one_letter_code', 'Std_deriv_three_letter_code', 'Std_deriv_BMRB_code', 'Std_deriv_PDB_code', 'Std_deriv_chem_comp_name', 'Formal_charge', 'Paramagnetic', 'Aromatic', 'Empirical_formula', 'Formula_weight', 'Formula_mono_iso_wt_nat', 'Formula_mono_iso_wt_13C', 'Formula_mono_iso_wt_15N', 'Formula_mono_iso_wt_13C_15N', 'Image_file_name', 'Image_file_format', 'Topo_file_name', 'Topo_file_format', 'Struct_file_name', 'Struct_file_format', 'Stereochem_param_file_name', 'Stereochem_param_file_format', 'Vendor', 'Vendor_product_code', 'Details', 'DB_query_date', 'DB_last_query_revised_last_date'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Chem_comp_common_name': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 1024),None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['Name', 'Type', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['Name', 'Entry_ID', 'Comp_ID'],

            },

      'Chem_comp_descriptor': {

        'tags': {

          'Ordinal': [None,returnStarInt,None,True],
          'Descriptor': [None,lambda x = value: returnStarString(x,length = 1024),None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Program': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Program_version': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['Ordinal', 'Descriptor', 'Type', 'Program', 'Program_version', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['Ordinal', 'Entry_ID', 'Comp_ID'],

            },

      'Chem_comp_identifier': {

        'tags': {

          'Ordinal': [None,returnStarInt,None,True],
          'Identifier': [None,lambda x = value: returnStarString(x,length = 1024),None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Program': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Program_version': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_id': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['Ordinal', 'Identifier', 'Type', 'Program', 'Program_version', 'Entry_ID', 'Comp_id'],
        'sourcePrimaryKeys': ['Ordinal', 'Entry_ID', 'Comp_id'],

            },

      'Chem_comp_systematic_name': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 1024),None,True],
          'Naming_system': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['Name', 'Naming_system', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['Naming_system', 'Entry_ID', 'Comp_ID'],

            },

      'Chem_comp_SMILES': {

        'tags': {

          'Type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'String': [None,lambda x = value: returnStarString(x,length = 1024),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['Type', 'String', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['Type', 'String', 'Entry_ID', 'Comp_ID'],

            },

      'Chem_comp_keyword': {

        'tags': {

          'Keyword': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['Keyword', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['Keyword', 'Entry_ID', 'Comp_ID'],

            },

      'Characteristic': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Chemical_group': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,True],
          'Source': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['ID', 'Name', 'Atom_ID', 'Chemical_group', 'Val', 'Val_err', 'Source', 'Citation_ID', 'Citation_label', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Comp_ID'],

            },

      'Chem_comp_atom': {

        'tags': {

          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),None,True],
          'PDB_atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),None,False],
          'Alt_atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Type_symbol': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Isotope_number': [None,returnStarInt,None,False],
          'Chirality': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Charge': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Partial_charge': [None,returnStarFloat,None,False],
          'Oxidation_number': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Unpaired_electron_number': [None,lambda x = value: returnStarString(x,length = 3),None,False],
          'PDBx_aromatic_flag': [None,lambda x = value: returnStarCode(x,length = 3),None,False],
          'PDBx_leaving_atom_flag': [None,lambda x = value: returnStarCode(x,length = 3),None,False],
          'Substruct_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Ionizable': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
          '2D_drawing_coord_x': [None,returnStarFloat,None,False],
          '2D_drawing_coord_y': [None,returnStarFloat,None,False],
          'Model_Cartn_x': [None,returnStarFloat,None,False],
          'Model_Cartn_y': [None,returnStarFloat,None,False],
          'Model_Cartn_z': [None,returnStarFloat,None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['Atom_ID', 'PDB_atom_ID', 'Alt_atom_ID', 'Auth_atom_ID', 'Type_symbol', 'Isotope_number', 'Chirality', 'Charge', 'Partial_charge', 'Oxidation_number', 'Unpaired_electron_number', 'PDBx_aromatic_flag', 'PDBx_leaving_atom_flag', 'Substruct_code', 'Ionizable', '2D_drawing_coord_x', '2D_drawing_coord_y', 'Model_Cartn_x', 'Model_Cartn_y', 'Model_Cartn_z', 'Details', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['Atom_ID', 'Entry_ID', 'Comp_ID'],

            },

      'Atom_nomenclature': {

        'tags': {

          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_name': [None,lambda x = value: returnStarLine(x,length = 15),None,True],
          'Naming_system': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['Atom_ID', 'Atom_name', 'Naming_system', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['Atom_ID', 'Entry_ID', 'Comp_ID'],

            },

      'Chem_comp_bond': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Value_order': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'PDB_atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),None,False],
          'PDB_atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['ID', 'Type', 'Value_order', 'Atom_ID_1', 'Atom_ID_2', 'PDB_atom_ID_1', 'PDB_atom_ID_2', 'Details', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Comp_ID'],

            },

      'Chem_comp_tor': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_ID_3': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_ID_4': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['ID', 'Atom_ID_1', 'Atom_ID_2', 'Atom_ID_3', 'Atom_ID_4', 'Details', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Comp_ID'],

            },

      'Chem_comp_angle': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_ID_3': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['ID', 'Atom_ID_1', 'Atom_ID_2', 'Atom_ID_3', 'Details', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Comp_ID'],

            },

      'Chem_comp_db_link': {

        'tags': {

          'Author_supplied': [None,lambda x = value: returnStarYesNo(x,length = 12),None,False],
          'Database_code': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Accession_code': [None,lambda x = value: returnStarLine(x,length = 15),None,True],
          'Accession_code_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_mol_code': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_mol_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_experimental_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_relation_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['Author_supplied', 'Database_code', 'Accession_code', 'Accession_code_type', 'Entry_mol_code', 'Entry_mol_name', 'Entry_experimental_method', 'Entry_relation_type', 'Entry_details', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['Database_code', 'Accession_code', 'Entry_ID', 'Comp_ID'],

            },

      'Chem_comp_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'Comp_ID'],
        'sourcePrimaryKeys': ['Citation_ID', 'Entry_ID', 'Comp_ID'],

            }

        },

    'tableNames': ['Chem_comp_common_name', 'Chem_comp_descriptor', 'Chem_comp_identifier', 'Chem_comp_systematic_name', 'Chem_comp_SMILES', 'Chem_comp_keyword', 'Characteristic', 'Chem_comp_atom', 'Atom_nomenclature', 'Chem_comp_bond', 'Chem_comp_tor', 'Chem_comp_angle', 'Chem_comp_db_link', 'Chem_comp_citation']

    },

  'sample': {

    'name': 'Sample',

    'tags': {

      'Sf_category': ['sample',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Sub_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],
      'Aggregate_sample_number': [None,returnStarInt,None,False],
      'Solvent_system': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Preparation_date': [None,returnStarDateTime,None,False],
      'Preparation_expiration_date': [None,returnStarDateTime,None,False],
      'Polycrystallization_protocol': [None,returnStarString,None,False],
      'Single_crystal_protocol': [None,returnStarString,None,False],
      'Crystal_grow_apparatus': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Crystal_grow_atmosphere': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Crystal_grow_details': [None,returnStarString,None,False],
      'Crystal_grow_method': [None,returnStarString,None,False],
      'Crystal_grow_method_cit_ID': [None,returnStarInt,'Citation.ID',False],
      'Crystal_grow_pH': [None,returnStarFloat,None,False],
      'Crystal_grow_pH_range': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Crystal_grow_pressure': [None,returnStarFloat,None,False],
      'Crystal_grow_pressure_esd': [None,returnStarFloat,None,False],
      'Crystal_grow_seeding': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Crystal_grow_seeding_cit_ID': [None,returnStarInt,'Citation.ID',False],
      'Crystal_grow_temp': [None,returnStarFloat,None,False],
      'Crystal_grow_temp_details': [None,returnStarString,None,False],
      'Crystal_grow_temp_esd': [None,returnStarFloat,None,False],
      'Crystal_grow_time': [None,returnStarFloat,None,False],
      'Oriented_sample_prep_protocol': [None,returnStarString,None,False],
      'Lyophilization_cryo_protectant': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Storage_protocol': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Type', 'Sub_type', 'Details', 'Aggregate_sample_number', 'Solvent_system', 'Preparation_date', 'Preparation_expiration_date', 'Polycrystallization_protocol', 'Single_crystal_protocol', 'Crystal_grow_apparatus', 'Crystal_grow_atmosphere', 'Crystal_grow_details', 'Crystal_grow_method', 'Crystal_grow_method_cit_ID', 'Crystal_grow_pH', 'Crystal_grow_pH_range', 'Crystal_grow_pressure', 'Crystal_grow_pressure_esd', 'Crystal_grow_seeding', 'Crystal_grow_seeding_cit_ID', 'Crystal_grow_temp', 'Crystal_grow_temp_details', 'Crystal_grow_temp_esd', 'Crystal_grow_time', 'Oriented_sample_prep_protocol', 'Lyophilization_cryo_protectant', 'Storage_protocol'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Sample_component': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Mol_common_name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Isotopic_labeling': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Assembly_ID': [None,returnStarInt,'Assembly.ID',False],
          'Assembly_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entity_ID': [None,returnStarInt,'Entity.ID',False],
          'Entity_label': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',False],
          'Product_ID': [None,returnStarInt,None,False],
          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Concentration_val': [None,lambda x = value: returnStarString(x,length = 31),None,False],
          'Concentration_val_min': [None,returnStarFloat,None,False],
          'Concentration_val_max': [None,returnStarFloat,None,False],
          'Concentration_val_units': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Concentration_val_err': [None,returnStarFloat,None,False],
          'Vendor': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Vendor_product_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Vendor_product_code': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],

                },

        'tagNames': ['ID', 'Mol_common_name', 'Isotopic_labeling', 'Assembly_ID', 'Assembly_label', 'Entity_ID', 'Entity_label', 'Product_ID', 'Type', 'Concentration_val', 'Concentration_val_min', 'Concentration_val_max', 'Concentration_val_units', 'Concentration_val_err', 'Vendor', 'Vendor_product_name', 'Vendor_product_code', 'Entry_ID', 'Sample_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Sample_ID'],

            },

      'Sample_component_atom_isotope': {

        'tags': {

          'Sample_component_ID': [None,returnStarInt,'Sample_component.ID',True],
          'Mol_common_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Entity_label': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Comp_isotope_label_code': [None,lambda x = value: returnStarLine(x,length = 12),None,True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarLine(x,length = 12),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,True],
          'Label_pct': [None,returnStarFloat,None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],

                },

        'tagNames': ['Sample_component_ID', 'Mol_common_name', 'Assembly_atom_ID', 'Entity_ID', 'Entity_label', 'Comp_ID', 'Comp_isotope_label_code', 'Comp_index_ID', 'Seq_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Label_pct', 'Entry_ID', 'Sample_ID'],
        'sourcePrimaryKeys': ['Sample_component_ID', 'Entity_ID', 'Comp_index_ID', 'Atom_ID', 'Atom_isotope_number', 'Entry_ID', 'Sample_ID'],

            },

      'Sample_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'Sample_ID'],
        'sourcePrimaryKeys': ['Citation_ID', 'Entry_ID', 'Sample_ID'],

            }

        },

    'tableNames': ['Sample_component', 'Sample_component_atom_isotope', 'Sample_citation']

    },

  'sample_conditions': {

    'name': 'Sample_condition_list',

    'tags': {

      'Sf_category': ['sample_conditions',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Sample_condition_variable': {

        'tags': {

          'Type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Val': [None,lambda x = value: returnStarString(x,length = 31),None,False],
          'Val_err': [None,lambda x = value: returnStarString(x,length = 31),None,False],
          'Val_units': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],

                },

        'tagNames': ['Type', 'Val', 'Val_err', 'Val_units', 'Entry_ID', 'Sample_condition_list_ID'],
        'sourcePrimaryKeys': ['Type'],

            },

      'Sample_condition_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'Sample_condition_list_ID'],
        'sourcePrimaryKeys': ['Entry_ID', 'Sample_condition_list_ID'],

            }

        },

    'tableNames': ['Sample_condition_variable', 'Sample_condition_citation']

    },

  'molecule_purity': {

    'name': 'Entity_purity_list',
    'saveFrameCode': 'sample_mol_purity',

    'tags': {

      'Sf_category': ['molecule_purity',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Entity_purity': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Sample_ID': [None,returnStarInt,'Sample.ID',False],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Entity_label': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
          'Val': [None,returnStarFloat,None,True],
          'Val_units': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Measurement_method': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_purity_list_ID': [None,returnStarInt,'Entity_purity_list.ID',True],

                },

        'tagNames': ['ID', 'Sample_ID', 'Sample_label', 'Entity_ID', 'Entity_label', 'Val', 'Val_units', 'Measurement_method', 'Details', 'Entry_ID', 'Entity_purity_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Entity_purity_list_ID'],

            },

      'Entity_purity_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Entity_purity_list_ID': [None,returnStarInt,'Entity_purity_list.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'Entity_purity_list_ID'],
        'sourcePrimaryKeys': ['Citation_ID', 'Entry_ID', 'Entity_purity_list_ID'],

            }

        },

    'tableNames': ['Entity_purity', 'Entity_purity_citation']

    },

  'software': {

    'name': 'Software',

    'tags': {

      'Sf_category': ['software',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Version': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Name', 'Version', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Vendor': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Address': [None,returnStarString,None,False],
          'Electronic_address': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Software_ID': [None,returnStarInt,'Software.ID',True],

                },

        'tagNames': ['Name', 'Address', 'Electronic_address', 'Entry_ID', 'Software_ID'],
        'sourcePrimaryKeys': ['Name', 'Entry_ID', 'Software_ID'],

            },

      'Task': {

        'tags': {

          'Task': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Software_ID': [None,returnStarInt,'Software.ID',True],

                },

        'tagNames': ['Task', 'Entry_ID', 'Software_ID'],
        'sourcePrimaryKeys': ['Task', 'Entry_ID', 'Software_ID'],

            },

      'Software_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Software_ID': [None,returnStarInt,'Software.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'Software_ID'],
        'sourcePrimaryKeys': ['Citation_ID', 'Entry_ID', 'Software_ID'],

            }

        },

    'tableNames': ['Vendor', 'Task', 'Software_citation']

    },

  'method': {

    'name': 'Method',

    'tags': {

      'Sf_category': ['method',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Derivation_type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Details': [None,returnStarString,None,False],
      'Computer_ID': [None,returnStarInt,'Computer.ID',False],
      'Computer_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Derivation_type', 'Details', 'Computer_ID', 'Computer_label'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Method_file': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Text_format': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Text': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Method_ID': [None,returnStarInt,'Method.ID',True],

                },

        'tagNames': ['Name', 'Text_format', 'Text', 'Entry_ID', 'Method_ID'],
        'sourcePrimaryKeys': ['Name', 'Entry_ID', 'Method_ID'],

            },

      'Method_param': {

        'tags': {

          'File_name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Text_format': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Text': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Method_ID': [None,returnStarInt,'Method.ID',True],

                },

        'tagNames': ['File_name', 'Text_format', 'Text', 'Entry_ID', 'Method_ID'],
        'sourcePrimaryKeys': ['File_name', 'Entry_ID', 'Method_ID'],

            },

      'Method_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Method_ID': [None,returnStarInt,'Method.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'Method_ID'],
        'sourcePrimaryKeys': ['Citation_ID', 'Entry_ID', 'Method_ID'],

            }

        },

    'tableNames': ['Method_file', 'Method_param', 'Method_citation']

    },

  'Mass_spectrometer': {

    'name': 'Mass_spectrometer',

    'tags': {

      'Sf_category': ['Mass_spectrometer',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarString(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Manufacturer': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Model': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Type': [None,lambda x = value: returnStarString(x,length = 127),None,True],
      'Serial_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Mass_ref_introduction_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Manufacturer', 'Model', 'Type', 'Serial_number', 'Mass_ref_introduction_method', 'Details'],
    'sourcePrimaryKeys': ['ID'],

    'tables': {

      'Mass_spectrometer_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Mass_spectrometer_ID': [None,returnStarInt,'Mass_spectrometer.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'Mass_spectrometer_ID'],
        'sourcePrimaryKeys': ['Citation_ID'],

            }

        },

    'tableNames': ['Mass_spectrometer_citation']

    },

  'Mass_spectrometer_list': {

    'name': 'Mass_spectrometer_list',
    'saveFrameCode': 'Mass_spectrometer_list',

    'tags': {

      'Sf_category': ['Mass_spectrometer_list',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarString(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID'],
    'sourcePrimaryKeys': ['ID'],

    'tables': {

      'Mass_spectrometer_view': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Name': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Manufacturer': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Model': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Type': [None,lambda x = value: returnStarString(x,length = 127),None,True],
          'Serial_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Mass_ref_introduction_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Details': [None,returnStarString,None,False],
          'Citation_ID': [None,returnStarInt,'Citation.ID',False],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',False],
          'Mass_spectrometer_list_ID': [None,returnStarInt,'Mass_spectrometer_list.ID',False],

                },

        'tagNames': ['ID', 'Name', 'Manufacturer', 'Model', 'Type', 'Serial_number', 'Mass_ref_introduction_method', 'Details', 'Citation_ID', 'Citation_label', 'Entry_ID', 'Mass_spectrometer_list_ID'],
        'sourcePrimaryKeys': ['ID'],

            }

        },

    'tableNames': ['Mass_spectrometer_view']

    },

  'Mass_spec_ref_compd': {

    'name': 'Mass_spec_ref_compd_set',

    'tags': {

      'Sf_category': ['Mass_spec_ref_compd',lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sf_framecode': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,True],
      'ID': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Mass_spec_ref_compd': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Mono_isotopic_mass': [None,returnStarFloat,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',False],
          'Mass_spec_ref_compd_set_ID': [None,returnStarInt,'Mass_spec_ref_compd_set.ID',False],

                },

        'tagNames': ['ID', 'Name', 'Mono_isotopic_mass', 'Entry_ID', 'Mass_spec_ref_compd_set_ID'],
        'sourcePrimaryKeys': [],

            }

        },

    'tableNames': ['Mass_spec_ref_compd']

    },

  'chromatographic_system': {

    'name': 'Chromatographic_system',

    'tags': {

      'Sf_category': ['chromatographic_system',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarString(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Pump_manufacturer': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Pump_model': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Pump_serial_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Autosampler_manufacturer': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Autosampler_model': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Autosampler_serial_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Col_compartment_manufacturer': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Col_compartment_model': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Col_compartment_serial_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,lambda x = value: returnStarLine(x,length = 127),None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Type', 'Pump_manufacturer', 'Pump_model', 'Pump_serial_number', 'Autosampler_manufacturer', 'Autosampler_model', 'Autosampler_serial_number', 'Col_compartment_manufacturer', 'Col_compartment_model', 'Col_compartment_serial_number', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    },

  'chromatographic_system': {

    'name': 'Chromatographic_system',

    'tags': {

      'Sf_category': ['chromatographic_system',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarString(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Vendor': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Packing_material': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Packing_material_pore_size': [None,returnStarFloat,None,False],
      'Width': [None,returnStarFloat,None,False],
      'Length': [None,returnStarFloat,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Vendor', 'Type', 'Packing_material', 'Packing_material_pore_size', 'Width', 'Length'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    },

  'NMR_spectrometer': {

    'name': 'NMR_spectrometer',

    'tags': {

      'Sf_category': ['NMR_spectrometer',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],
      'Manufacturer': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Model': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Serial_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Field_strength': [None,lambda x = value: returnStarString(x,length = 127),None,True],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Details', 'Manufacturer', 'Model', 'Serial_number', 'Field_strength'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'NMR_spectrometer_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'NMR_spectrometer_ID': [None,returnStarInt,'NMR_spectrometer.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'NMR_spectrometer_ID'],
        'sourcePrimaryKeys': ['Citation_ID', 'Entry_ID', 'NMR_spectrometer_ID'],

            }

        },

    'tableNames': ['NMR_spectrometer_citation']

    },

  'NMR_spectrometer_list': {

    'name': 'NMR_spectrometer_list',
    'saveFrameCode': 'NMR_spectrometer_list',

    'tags': {

      'Sf_category': ['NMR_spectrometer_list',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'NMR_spectrometer_view': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Manufacturer': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Model': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Serial_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Field_strength': [None,lambda x = value: returnStarString(x,length = 127),None,True],
          'Details': [None,returnStarString,None,False],
          'Citation_ID': [None,returnStarInt,'Citation.ID',False],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'NMR_spectrometer_list_ID': [None,returnStarInt,'NMR_spectrometer_list.ID',True],

                },

        'tagNames': ['ID', 'Name', 'Manufacturer', 'Model', 'Serial_number', 'Field_strength', 'Details', 'Citation_ID', 'Citation_label', 'Entry_ID', 'NMR_spectrometer_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Name', 'Entry_ID', 'NMR_spectrometer_list_ID'],

            }

        },

    'tableNames': ['NMR_spectrometer_view']

    },

  'NMR_spectrometer_probe': {

    'name': 'NMR_spectrometer_probe',

    'tags': {

      'Sf_category': ['NMR_spectrometer_probe',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],
      'Manufacturer': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Model': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Serial_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Diameter': [None,lambda x = value: returnStarString(x,length = 127),None,False],
      'Rotor_length': [None,lambda x = value: returnStarString(x,length = 127),None,False],
      'Rotor_composition': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Internal_volume': [None,returnStarFloat,None,False],
      'Spacer_present': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Details', 'Manufacturer', 'Model', 'Serial_number', 'Diameter', 'Rotor_length', 'Rotor_composition', 'Internal_volume', 'Spacer_present'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'NMR_probe': {

        'tags': {

          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'NMR_spectrometer_probe_ID': [None,returnStarInt,'NMR_spectrometer_probe.ID',True],

                },

        'tagNames': ['Type', 'Entry_ID', 'NMR_spectrometer_probe_ID'],
        'sourcePrimaryKeys': ['Type', 'Entry_ID', 'NMR_spectrometer_probe_ID'],

            },

      'NMR_spectrometer_probe_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'NMR_spectrometer_probe_ID': [None,returnStarInt,'NMR_spectrometer_probe.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'NMR_spectrometer_probe_ID'],
        'sourcePrimaryKeys': ['Citation_ID', 'Entry_ID', 'NMR_spectrometer_probe_ID'],

            }

        },

    'tableNames': ['NMR_probe', 'NMR_spectrometer_probe_citation']

    },

  'experiment_list': {

    'name': 'Experiment_list',
    'saveFrameCode': 'experiment_list',

    'tags': {

      'Sf_category': ['experiment_list',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Experiment': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Raw_data_flag': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
          'NMR_spec_expt_ID': [None,returnStarInt,'NMR_spec_expt.ID',False],
          'NMR_spec_expt_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',False],
          'Sample_label': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Sample_volume': [None,returnStarFloat,None,False],
          'Sample_volume_units': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',False],
          'Sample_condition_list_label': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Sample_spinning_rate': [None,returnStarFloat,None,False],
          'Sample_angle': [None,returnStarFloat,None,False],
          'NMR_tube_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'NMR_spectrometer_ID': [None,returnStarInt,'NMR_spectrometer.ID',False],
          'NMR_spectrometer_label': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'NMR_spectrometer_probe_ID': [None,returnStarInt,'NMR_spectrometer_probe.ID',False],
          'NMR_spectrometer_probe_label': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'NMR_spectral_processing_ID': [None,returnStarInt,'NMR_spectral_processing.ID',False],
          'NMR_spectral_processing_label': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Experiment_list_ID': [None,returnStarInt,'Experiment_list.ID',True],

                },

        'tagNames': ['ID', 'Name', 'Raw_data_flag', 'NMR_spec_expt_ID', 'NMR_spec_expt_label', 'Sample_ID', 'Sample_label', 'Sample_state', 'Sample_volume', 'Sample_volume_units', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Sample_spinning_rate', 'Sample_angle', 'NMR_tube_type', 'NMR_spectrometer_ID', 'NMR_spectrometer_label', 'NMR_spectrometer_probe_ID', 'NMR_spectrometer_probe_label', 'NMR_spectral_processing_ID', 'NMR_spectral_processing_label', 'Entry_ID', 'Experiment_list_ID'],
        'sourcePrimaryKeys': ['ID'],

            },

      'Experiment_file': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Directory_path': [None,returnStarString,None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Experiment_list_ID': [None,returnStarInt,'Experiment_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Name', 'Type', 'Directory_path', 'Details', 'Entry_ID', 'Experiment_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID'],

            }

        },

    'tableNames': ['Experiment', 'Experiment_file']

    },

  'NMR_spectrometer_expt': {

    'name': 'NMR_spec_expt',

    'tags': {

      'Sf_category': ['NMR_spectrometer_expt',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_volume': [None,returnStarFloat,None,False],
      'Sample_volume_units': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'NMR_tube_type': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Sample_spinning_rate': [None,returnStarFloat,None,False],
      'Sample_angle': [None,returnStarFloat,None,False],
      'NMR_spectrometer_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'NMR_spectrometer_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'NMR_spectrometer_probe_ID': [None,returnStarInt,'NMR_spectrometer_probe.ID',True],
      'NMR_spectrometer_probe_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Carrier_freq_switch_time': [None,lambda x = value: returnStarString(x,length = 127),None,False],
      'Software_ID': [None,returnStarInt,'Software.ID',True],
      'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Method_ID': [None,returnStarInt,'Method.ID',False],
      'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Pulse_seq_accession_BMRB_code': [None,returnStarInt,None,False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Name', 'Type', 'Sample_volume', 'Sample_volume_units', 'NMR_tube_type', 'Sample_spinning_rate', 'Sample_angle', 'NMR_spectrometer_ID', 'NMR_spectrometer_label', 'NMR_spectrometer_probe_ID', 'NMR_spectrometer_probe_label', 'Carrier_freq_switch_time', 'Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Pulse_seq_accession_BMRB_code', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'NMR_expt_systematic_name': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 255),None,False],
          'Naming_system': [None,lambda x = value: returnStarLine(x,length = 255),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'NMR_spec_expt_ID': [None,returnStarInt,'NMR_spec_expt.ID',True],

                },

        'tagNames': ['Name', 'Naming_system', 'Entry_ID', 'NMR_spec_expt_ID'],
        'sourcePrimaryKeys': ['Entry_ID', 'NMR_spec_expt_ID'],

            },

      'NMR_experiment_file': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Directory_path': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Byte_order': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Bytes_per_data_point': [None,returnStarInt,None,False],
          'File_header_size': [None,returnStarInt,None,False],
          'Record_header_size': [None,returnStarInt,None,False],
          'Record_trailer_size': [None,returnStarInt,None,False],
          'Compression_algorithm': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'NMR_spec_expt_ID': [None,returnStarInt,'NMR_spec_expt.ID',True],

                },

        'tagNames': ['Name', 'Type', 'Directory_path', 'Byte_order', 'Bytes_per_data_point', 'File_header_size', 'Record_header_size', 'Record_trailer_size', 'Compression_algorithm', 'Details', 'Entry_ID', 'NMR_spec_expt_ID'],
        'sourcePrimaryKeys': ['Name', 'Entry_ID', 'NMR_spec_expt_ID'],

            },

      'Spectral_acq_param': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Acquisition_dimension_ID': [None,returnStarInt,None,True],
          'Name': [None,lambda x = value: returnStarLine(x,length = 15),None,True],
          'Val': [None,returnStarFloat,None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'NMR_spec_expt_ID': [None,returnStarInt,'NMR_spec_expt.ID',True],

                },

        'tagNames': ['ID', 'Acquisition_dimension_ID', 'Name', 'Val', 'Entry_ID', 'NMR_spec_expt_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'NMR_spec_expt_ID'],

            },

      'Recoupling_pulse_sequence': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Nucleus': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Time_period': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'NMR_spec_expt_ID': [None,returnStarInt,'NMR_spec_expt.ID',True],

                },

        'tagNames': ['Name', 'Type', 'Nucleus', 'Time_period', 'Entry_ID', 'NMR_spec_expt_ID'],
        'sourcePrimaryKeys': ['Name', 'Time_period', 'Entry_ID', 'NMR_spec_expt_ID'],

            },

      'Decoupling_pulse_sequence': {

        'tags': {

          'Name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Time_period': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'NMR_spec_expt_ID': [None,returnStarInt,'NMR_spec_expt.ID',True],

                },

        'tagNames': ['Name', 'Time_period', 'Entry_ID', 'NMR_spec_expt_ID'],
        'sourcePrimaryKeys': ['Name', 'Time_period', 'Entry_ID', 'NMR_spec_expt_ID'],

            },

      'NMR_experiment_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'NMR_spec_expt_ID': [None,returnStarInt,'NMR_spec_expt.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'NMR_spec_expt_ID'],
        'sourcePrimaryKeys': ['Citation_ID', 'Entry_ID', 'NMR_spec_expt_ID'],

            }

        },

    'tableNames': ['NMR_expt_systematic_name', 'NMR_experiment_file', 'Spectral_acq_param', 'Recoupling_pulse_sequence', 'Decoupling_pulse_sequence', 'NMR_experiment_citation']

    },

  'NMR_spectral_processing': {

    'name': 'NMR_spectral_processing',

    'tags': {

      'Sf_category': ['NMR_spectral_processing',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'NMR_spectral_proc_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'NMR_spectral_processing_ID': [None,returnStarInt,'NMR_spectral_processing.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'NMR_spectral_processing_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'NMR_spectral_processing_ID'],

            },

      'Spectral_processing_param': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Processing_dimension_ID': [None,returnStarInt,None,False],
          'Name': [None,returnStarFloat,None,True],
          'Val': [None,returnStarFloat,None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'NMR_spectral_processing_ID': [None,returnStarInt,'NMR_spectral_processing.ID',True],

                },

        'tagNames': ['ID', 'Processing_dimension_ID', 'Name', 'Val', 'Entry_ID', 'NMR_spectral_processing_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'NMR_spectral_processing_ID'],

            }

        },

    'tableNames': ['NMR_spectral_proc_software', 'Spectral_processing_param']

    },

  'computer': {

    'name': 'Computer',

    'tags': {

      'Sf_category': ['computer',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],
      'Manufacturer': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Model': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Operating_system': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Operating_system_version': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Hardware_code': [None,lambda x = value: returnStarLine(x,length = 127),None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Details', 'Manufacturer', 'Model', 'Operating_system', 'Operating_system_version', 'Hardware_code'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Computer_citation': {

        'tags': {

          'Citation_ID': [None,returnStarInt,'Citation.ID',True],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Computer_ID': [None,returnStarInt,'Computer.ID',True],

                },

        'tagNames': ['Citation_ID', 'Citation_label', 'Entry_ID', 'Computer_ID'],
        'sourcePrimaryKeys': ['Citation_ID', 'Entry_ID', 'Computer_ID'],

            }

        },

    'tableNames': ['Computer_citation']

    },

  'chem_shift_reference': {

    'name': 'Chem_shift_reference',

    'tags': {

      'Sf_category': ['chem_shift_reference',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Proton_shifts_flag': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Carbon_shifts_flag': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Nitrogen_shifts_flag': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Phosphorus_shifts_flag': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Other_shifts_flag': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Proton_shifts_flag', 'Carbon_shifts_flag', 'Nitrogen_shifts_flag', 'Phosphorus_shifts_flag', 'Other_shifts_flag', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Chem_shift_ref': {

        'tags': {

          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,True],
          'Mol_common_name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Atom_group': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Concentration_val': [None,returnStarFloat,None,False],
          'Concentration_units': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Solvent': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Rank': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Chem_shift_units': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Chem_shift_val': [None,returnStarFloat,None,True],
          'Ref_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Ref_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Indirect_shift_ratio': [None,returnStarFloat,None,False],
          'External_ref_loc': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'External_ref_sample_geometry': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'External_ref_axis': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Indirect_shift_ratio_cit_ID': [None,returnStarInt,'Citation.ID',False],
          'Indirect_shift_ratio_cit_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Ref_correction_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Correction_val': [None,returnStarFloat,None,False],
          'Correction_val_cit_ID': [None,returnStarInt,'Citation.ID',False],
          'Correction_val_cit_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Chem_shift_reference_ID': [None,returnStarInt,'Chem_shift_reference.ID',True],

                },

        'tagNames': ['Atom_type', 'Atom_isotope_number', 'Mol_common_name', 'Atom_group', 'Concentration_val', 'Concentration_units', 'Solvent', 'Rank', 'Chem_shift_units', 'Chem_shift_val', 'Ref_method', 'Ref_type', 'Indirect_shift_ratio', 'External_ref_loc', 'External_ref_sample_geometry', 'External_ref_axis', 'Indirect_shift_ratio_cit_ID', 'Indirect_shift_ratio_cit_label', 'Ref_correction_type', 'Correction_val', 'Correction_val_cit_ID', 'Correction_val_cit_label', 'Entry_ID', 'Chem_shift_reference_ID'],
        'sourcePrimaryKeys': ['Atom_type', 'Atom_isotope_number', 'Mol_common_name', 'Entry_ID', 'Chem_shift_reference_ID'],

            }

        },

    'tableNames': ['Chem_shift_ref']

    },

  'assigned_chemical_shifts': {

    'name': 'Assigned_chem_shift_list',

    'tags': {

      'Sf_category': ['assigned_chemical_shifts',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Chem_shift_reference_ID': [None,returnStarInt,'Chem_shift_reference.ID',True],
      'Chem_shift_reference_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Chem_shift_1H_err': [None,returnStarFloat,None,False],
      'Chem_shift_13C_err': [None,returnStarFloat,None,False],
      'Chem_shift_15N_err': [None,returnStarFloat,None,False],
      'Chem_shift_31P_err': [None,returnStarFloat,None,False],
      'Chem_shift_2H_err': [None,returnStarFloat,None,False],
      'Chem_shift_19F_err': [None,returnStarFloat,None,False],
      'Error_derivation_method': [None,returnStarString,None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Chem_shift_reference_ID', 'Chem_shift_reference_label', 'Chem_shift_1H_err', 'Chem_shift_13C_err', 'Chem_shift_15N_err', 'Chem_shift_31P_err', 'Chem_shift_2H_err', 'Chem_shift_19F_err', 'Error_derivation_method', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Chem_shift_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',False],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assigned_chem_shift_list_ID': [None,returnStarInt,'Assigned_chem_shift_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Assigned_chem_shift_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID'],

            },

      'Systematic_chem_shift_offset': {

        'tags': {

          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Val': [None,lambda x = value: returnStarString(x,length = 15),None,True],
          'Val_err': [None,lambda x = value: returnStarString(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assigned_chem_shift_list_ID': [None,returnStarInt,'Assigned_chem_shift_list.ID',True],

                },

        'tagNames': ['Type', 'Atom_type', 'Atom_isotope_number', 'Val', 'Val_err', 'Entry_ID', 'Assigned_chem_shift_list_ID'],
        'sourcePrimaryKeys': ['Type', 'Entry_ID', 'Assigned_chem_shift_list_ID'],

            },

      'Chem_shift_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assigned_chem_shift_list_ID': [None,returnStarInt,'Assigned_chem_shift_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Assigned_chem_shift_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Assigned_chem_shift_list_ID'],

            },

      'Atom_chem_shift': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',False],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,True],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,False],
          'Assign_fig_of_merit': [None,returnStarFloat,None,False],
          'Ambiguity_code': [None,lambda x = value: returnStarString(x,length = 3),None,False],
          'Occupancy': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assigned_chem_shift_list_ID': [None,returnStarInt,'Assigned_chem_shift_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Val', 'Val_err', 'Assign_fig_of_merit', 'Ambiguity_code', 'Occupancy', 'Resonance_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Details', 'Entry_ID', 'Assigned_chem_shift_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Assigned_chem_shift_list_ID'],

            },

      'Ambiguous_atom_chem_shift': {

        'tags': {

          'Ambiguous_shift_set_ID': [None,returnStarInt,None,True],
          'Atom_chem_shift_ID': [None,returnStarInt,'Atom_chem_shift.ID',True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Assigned_chem_shift_list_ID': [None,returnStarInt,'Assigned_chem_shift_list.ID',True],

                },

        'tagNames': ['Ambiguous_shift_set_ID', 'Atom_chem_shift_ID', 'Entry_ID', 'Assigned_chem_shift_list_ID'],
        'sourcePrimaryKeys': ['Ambiguous_shift_set_ID', 'Entry_ID', 'Assigned_chem_shift_list_ID'],

            }

        },

    'tableNames': ['Chem_shift_experiment', 'Systematic_chem_shift_offset', 'Chem_shift_software', 'Atom_chem_shift', 'Ambiguous_atom_chem_shift']

    },

  'coupling_constants': {

    'name': 'Coupling_constant_list',

    'tags': {

      'Sf_category': ['coupling_constants',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Spectrometer_frequency_1H': [None,returnStarFloat,None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Spectrometer_frequency_1H', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Coupling_constant_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',False],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Coupling_constant_list_ID': [None,returnStarInt,'Coupling_constant_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Coupling_constant_list_ID'],
        'sourcePrimaryKeys': ['Entry_ID', 'Coupling_constant_list_ID'],

            },

      'Coupling_constant_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Coupling_constant_list_ID': [None,returnStarInt,'Coupling_constant_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Coupling_constant_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Coupling_constant_list_ID'],

            },

      'Coupling_constant': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Code': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',False],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_1': [None,returnStarInt,None,False],
          'Ambiguity_code_1': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',False],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_2': [None,returnStarInt,None,False],
          'Ambiguity_code_2': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
          'Val': [None,returnStarFloat,None,False],
          'Val_min': [None,returnStarFloat,None,False],
          'Val_max': [None,returnStarFloat,None,False],
          'Val_err': [None,returnStarFloat,None,False],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Coupling_constant_list_ID': [None,returnStarInt,'Coupling_constant_list.ID',True],

                },

        'tagNames': ['ID', 'Code', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Atom_isotope_number_1', 'Ambiguity_code_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Atom_isotope_number_2', 'Ambiguity_code_2', 'Val', 'Val_min', 'Val_max', 'Val_err', 'Resonance_ID_1', 'Resonance_ID_2', 'Auth_entity_assembly_ID_1', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_entity_assembly_ID_2', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Details', 'Entry_ID', 'Coupling_constant_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Coupling_constant_list_ID'],

            }

        },

    'tableNames': ['Coupling_constant_experiment', 'Coupling_constant_software', 'Coupling_constant']

    },

  'spectral_peak_list': {

    'name': 'Spectral_peak_list',

    'tags': {

      'Sf_category': ['spectral_peak_list',lambda x = value: returnStarCode(x,length = 31),None,False],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_ID': [None,returnStarInt,'Sample.ID',True],
      'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
      'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Number_of_spectral_dimensions': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_ID', 'Sample_label', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Experiment_ID', 'Experiment_name', 'Number_of_spectral_dimensions', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Spectral_dim': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,True],
          'Spectral_region': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
          'Magnetization_linkage_ID': [None,returnStarInt,None,False],
          'Sweep_width': [None,returnStarFloat,None,False],
          'Encoding_code': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Encoded_source_dimension_ID': [None,returnStarInt,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_peak_list_ID': [None,returnStarInt,'Spectral_peak_list.ID',True],

                },

        'tagNames': ['ID', 'Atom_type', 'Atom_isotope_number', 'Spectral_region', 'Magnetization_linkage_ID', 'Sweep_width', 'Encoding_code', 'Encoded_source_dimension_ID', 'Entry_ID', 'Spectral_peak_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Atom_type', 'Spectral_region', 'Entry_ID', 'Spectral_peak_list_ID'],

            },

      'Spectral_peak_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_peak_list_ID': [None,returnStarInt,'Spectral_peak_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Spectral_peak_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Spectral_peak_list_ID'],

            },

      'Peak': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Figure_of_merit': [None,returnStarFloat,None,False],
          'Details': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_peak_list_ID': [None,returnStarInt,'Spectral_peak_list.ID',True],

                },

        'tagNames': ['ID', 'Figure_of_merit', 'Details', 'Entry_ID', 'Spectral_peak_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Spectral_peak_list_ID'],

            },

      'Peak_general_char': {

        'tags': {

          'Peak_ID': [None,returnStarInt,'Peak.ID',True],
          'Intensity_val': [None,returnStarFloat,None,True],
          'Intensity_val_err': [None,returnStarFloat,None,False],
          'Measurement_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_peak_list_ID': [None,returnStarInt,'Spectral_peak_list.ID',True],

                },

        'tagNames': ['Peak_ID', 'Intensity_val', 'Intensity_val_err', 'Measurement_method', 'Entry_ID', 'Spectral_peak_list_ID'],
        'sourcePrimaryKeys': ['Peak_ID', 'Intensity_val', 'Entry_ID', 'Spectral_peak_list_ID'],

            },

      'Peak_char': {

        'tags': {

          'Peak_ID': [None,returnStarInt,'Peak.ID',True],
          'Peak_contribution_ID': [None,returnStarInt,'Peak_contribution.ID',False],
          'Spectral_dim_ID': [None,returnStarInt,'Spectral_dim.ID',True],
          'Chem_shift_val': [None,returnStarFloat,None,True],
          'Chem_shift_val_err': [None,returnStarFloat,None,False],
          'Line_width_val': [None,returnStarFloat,None,False],
          'Line_width_val_err': [None,returnStarFloat,None,False],
          'Phase_val': [None,returnStarFloat,None,False],
          'Phase_val_err': [None,returnStarFloat,None,False],
          'Decay_rate_val': [None,returnStarFloat,None,False],
          'Decay_rate_val_err': [None,returnStarFloat,None,False],
          'Coupling_pattern': [None,lambda x = value: returnStarLine(x,length = 15),None,False],
          'Bounding_box_upper_val': [None,returnStarFloat,None,False],
          'Bounding_box_lower_val': [None,returnStarFloat,None,False],
          'Bounding_box_range_val': [None,returnStarFloat,None,False],
          'Details': [None,returnStarString,None,False],
          'Derivation_method_ID': [None,returnStarInt,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_peak_list_ID': [None,returnStarInt,'Spectral_peak_list.ID',True],

                },

        'tagNames': ['Peak_ID', 'Peak_contribution_ID', 'Spectral_dim_ID', 'Chem_shift_val', 'Chem_shift_val_err', 'Line_width_val', 'Line_width_val_err', 'Phase_val', 'Phase_val_err', 'Decay_rate_val', 'Decay_rate_val_err', 'Coupling_pattern', 'Bounding_box_upper_val', 'Bounding_box_lower_val', 'Bounding_box_range_val', 'Details', 'Derivation_method_ID', 'Entry_ID', 'Spectral_peak_list_ID'],
        'sourcePrimaryKeys': ['Peak_ID', 'Spectral_dim_ID', 'Entry_ID', 'Spectral_peak_list_ID'],

            },

      'Peak_contribution': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Contribution_fractional_val': [None,returnStarFloat,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',False],
          'Spectral_peak_list_ID': [None,returnStarInt,'Spectral_peak_list.ID',False],

                },

        'tagNames': ['ID', 'Contribution_fractional_val', 'Entry_ID', 'Spectral_peak_list_ID'],
        'sourcePrimaryKeys': ['ID'],

            },

      'Assigned_peak_chem_shift': {

        'tags': {

          'Peak_ID': [None,returnStarInt,'Peak.ID',True],
          'Spectral_dim_ID': [None,returnStarInt,'Spectral_dim.ID',True],
          'Set_ID': [None,returnStarInt,None,False],
          'Magnetization_linkage_ID': [None,returnStarInt,None,False],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Val': [None,returnStarFloat,None,False],
          'Figure_of_merit': [None,returnStarFloat,None,False],
          'Assigned_chem_shift_list_ID': [None,returnStarInt,'Assigned_chem_shift_list.ID',False],
          'Atom_chem_shift_ID': [None,returnStarInt,'Atom_chem_shift.ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',False],
          'Entity_ID': [None,returnStarInt,'Entity.ID',False],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_peak_list_ID': [None,returnStarInt,'Spectral_peak_list.ID',True],

                },

        'tagNames': ['Peak_ID', 'Spectral_dim_ID', 'Set_ID', 'Magnetization_linkage_ID', 'Assembly_atom_ID', 'Val', 'Figure_of_merit', 'Assigned_chem_shift_list_ID', 'Atom_chem_shift_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Comp_ID', 'Atom_ID', 'Resonance_ID', 'Details', 'Entry_ID', 'Spectral_peak_list_ID'],
        'sourcePrimaryKeys': ['Peak_ID', 'Spectral_dim_ID', 'Entry_ID', 'Spectral_peak_list_ID'],

            },

      'Spectral_transition': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Figure_of_merit': [None,returnStarFloat,None,False],
          'Details': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_peak_list_ID': [None,returnStarInt,'Spectral_peak_list.ID',True],

                },

        'tagNames': ['ID', 'Figure_of_merit', 'Details', 'Entry_ID', 'Spectral_peak_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Spectral_peak_list_ID'],

            },

      'Spectral_transition_general_char': {

        'tags': {

          'Spectral_transition_ID': [None,returnStarInt,'Peak.ID',True],
          'Intensity_val': [None,returnStarFloat,None,True],
          'Intensity_val_err': [None,returnStarFloat,None,False],
          'Measurement_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_peak_list_ID': [None,returnStarInt,'Spectral_peak_list.ID',True],

                },

        'tagNames': ['Spectral_transition_ID', 'Intensity_val', 'Intensity_val_err', 'Measurement_method', 'Entry_ID', 'Spectral_peak_list_ID'],
        'sourcePrimaryKeys': ['Spectral_transition_ID', 'Intensity_val', 'Entry_ID', 'Spectral_peak_list_ID'],

            },

      'Spectral_transition_char': {

        'tags': {

          'Spectral_transition_ID': [None,returnStarInt,'Spectral_transition.ID',True],
          'Spectral_dim_ID': [None,returnStarInt,'Spectral_dim.ID',True],
          'Spectral_transition_contrib_ID': [None,returnStarInt,'Spectral_transition_contrib.ID',False],
          'Chem_shift_val': [None,returnStarFloat,None,True],
          'Chem_shift_val_err': [None,returnStarFloat,None,False],
          'Line_width_val': [None,returnStarFloat,None,False],
          'Line_width_val_err': [None,returnStarFloat,None,False],
          'Phase_val': [None,returnStarFloat,None,False],
          'Phase_val_err': [None,returnStarFloat,None,False],
          'Decay_rate_val': [None,returnStarFloat,None,False],
          'Decay_rate_val_err': [None,returnStarFloat,None,False],
          'Bounding_box_upper_val': [None,returnStarFloat,None,False],
          'Bounding_box_lower_val': [None,returnStarFloat,None,False],
          'Bounding_box_width_val': [None,returnStarFloat,None,False],
          'Derivation_method_ID': [None,returnStarInt,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_peak_list_ID': [None,returnStarInt,'Spectral_peak_list.ID',True],

                },

        'tagNames': ['Spectral_transition_ID', 'Spectral_dim_ID', 'Spectral_transition_contrib_ID', 'Chem_shift_val', 'Chem_shift_val_err', 'Line_width_val', 'Line_width_val_err', 'Phase_val', 'Phase_val_err', 'Decay_rate_val', 'Decay_rate_val_err', 'Bounding_box_upper_val', 'Bounding_box_lower_val', 'Bounding_box_width_val', 'Derivation_method_ID', 'Entry_ID', 'Spectral_peak_list_ID'],
        'sourcePrimaryKeys': ['Spectral_transition_ID', 'Spectral_dim_ID', 'Entry_ID', 'Spectral_peak_list_ID'],

            },

      'Spectral_transition_contrib': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Contribution_fractional_val': [None,returnStarFloat,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',False],
          'Spectral_peak_list_ID': [None,returnStarInt,'Spectral_peak_list.ID',False],

                },

        'tagNames': ['ID', 'Contribution_fractional_val', 'Entry_ID', 'Spectral_peak_list_ID'],
        'sourcePrimaryKeys': ['ID'],

            },

      'Assigned_spectral_transition': {

        'tags': {

          'Spectral_transition_ID': [None,returnStarInt,'Spectral_transition.ID',True],
          'Spectral_dim_ID': [None,returnStarInt,'Spectral_dim.ID',True],
          'Peak_ID': [None,returnStarInt,'Peak.ID',False],
          'Figure_of_merit': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_peak_list_ID': [None,returnStarInt,'Spectral_peak_list.ID',True],

                },

        'tagNames': ['Spectral_transition_ID', 'Spectral_dim_ID', 'Peak_ID', 'Figure_of_merit', 'Resonance_ID', 'Details', 'Entry_ID', 'Spectral_peak_list_ID'],
        'sourcePrimaryKeys': ['Spectral_transition_ID', 'Spectral_dim_ID', 'Entry_ID', 'Spectral_peak_list_ID'],

            }

        },

    'tableNames': ['Spectral_dim', 'Spectral_peak_software', 'Peak', 'Peak_general_char', 'Peak_char', 'Peak_contribution', 'Assigned_peak_chem_shift', 'Spectral_transition', 'Spectral_transition_general_char', 'Spectral_transition_char', 'Spectral_transition_contrib', 'Assigned_spectral_transition']

    },

  'resonance_linker': {

    'name': 'Resonance_linker_list',

    'tags': {

      'Sf_category': ['resonance_linker',lambda x = value: returnStarCode(x,length = 31),None,False],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Resonance': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Name': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Resonance_set_ID': [None,returnStarInt,'Resonance_assignment.Resonance_set_ID',False],
          'Spin_system_ID': [None,returnStarInt,'Spin_system.ID',False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Resonance_linker_list_ID': [None,returnStarInt,'Resonance_linker_list.ID',True],

                },

        'tagNames': ['ID', 'Name', 'Resonance_set_ID', 'Spin_system_ID', 'Entry_ID', 'Resonance_linker_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Resonance_linker_list_ID'],

            },

      'Resonance_covalent_link': {

        'tags': {

          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Resonance_linker_list_ID': [None,returnStarInt,None,False],

                },

        'tagNames': ['Resonance_ID_1', 'Resonance_ID_2', 'Entry_ID', 'Resonance_linker_list_ID'],
        'sourcePrimaryKeys': [],

            },

      'Resonance_assignment': {

        'tags': {

          'Resonance_set_ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',False],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Atom_set_ID': [None,returnStarInt,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Resonance_linker_list_ID': [None,returnStarInt,'Resonance_linker_list.ID',True],

                },

        'tagNames': ['Resonance_set_ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Atom_set_ID', 'Entry_ID', 'Resonance_linker_list_ID'],
        'sourcePrimaryKeys': ['Resonance_set_ID', 'Entity_assembly_ID', 'Comp_index_ID', 'Atom_ID', 'Entry_ID', 'Resonance_linker_list_ID'],

            },

      'Spin_system': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',False],
          'Entity_ID': [None,returnStarInt,'Entity.ID',False],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Resonance_linker_list_ID': [None,returnStarInt,'Resonance_linker_list.ID',True],

                },

        'tagNames': ['ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Comp_ID', 'Entry_ID', 'Resonance_linker_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Resonance_linker_list_ID'],

            },

      'Spin_system_link': {

        'tags': {

          'From_spin_system_ID': [None,returnStarInt,None,False],
          'To_spin_system_ID': [None,returnStarInt,None,False],
          'Offset': [None,returnStarInt,None,False],
          'Type': [None,lambda x = value: returnStarString(x,length = 31),None,False],
          'Selected': [None,lambda x = value: returnStarString(x,length = 3),None,False],
          'Entry_ID': [None,lambda x = value: returnStarString(x,length = 12),None,False],
          'Resonance_linker_list_ID': [None,returnStarInt,None,False],

                },

        'tagNames': ['From_spin_system_ID', 'To_spin_system_ID', 'Offset', 'Type', 'Selected', 'Entry_ID', 'Resonance_linker_list_ID'],
        'sourcePrimaryKeys': [],

            }

        },

    'tableNames': ['Resonance', 'Resonance_covalent_link', 'Resonance_assignment', 'Spin_system', 'Spin_system_link']

    },

  'chem_shift_isotope_effect': {

    'name': 'Chem_shift_isotope_effect_list',

    'tags': {

      'Sf_category': ['chem_shift_isotope_effect',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Isotope_effect_type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Isotope_effect_val_units': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Isotope_effect_type', 'Isotope_effect_val_units', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Isotope_effect_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Chem_shift_isotope_effect_list_ID': [None,returnStarInt,'Chem_shift_isotope_effect_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Chem_shift_isotope_effect_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Chem_shift_isotope_effect_list_ID'],

            },

      'Isotope_effect_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',True],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Chem_shift_isotope_effect_list_ID': [None,returnStarInt,'Chem_shift_isotope_effect_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Chem_shift_isotope_effect_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Chem_shift_isotope_effect_list_ID'],

            },

      'Isotope_effect': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Label_sample_ID_1': [None,returnStarInt,'Sample.ID',True],
          'Label_pattern_ID_1': [None,returnStarInt,'Isotope_label_pattern.ID',False],
          'Label_1_ID_chem_shift_val': [None,returnStarFloat,None,False],
          'Label_1_ID_chem_shift_val_err': [None,returnStarFloat,None,False],
          'Label_sample_ID_2': [None,returnStarInt,'Sample.ID',True],
          'Label_pattern_ID_2': [None,returnStarInt,'Isotope_label_pattern.ID',False],
          'Label_2_ID_chem_shift_val': [None,returnStarFloat,None,False],
          'Label_2_ID_chem_shift_val_err': [None,returnStarFloat,None,False],
          'Chem_shift_val': [None,returnStarFloat,None,False],
          'Chem_shift_val_err': [None,returnStarFloat,None,False],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Chem_shift_isotope_effect_list_ID': [None,returnStarInt,'Chem_shift_isotope_effect_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Label_sample_ID_1', 'Label_pattern_ID_1', 'Label_1_ID_chem_shift_val', 'Label_1_ID_chem_shift_val_err', 'Label_sample_ID_2', 'Label_pattern_ID_2', 'Label_2_ID_chem_shift_val', 'Label_2_ID_chem_shift_val_err', 'Chem_shift_val', 'Chem_shift_val_err', 'Resonance_ID_1', 'Resonance_ID_2', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Chem_shift_isotope_effect_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Chem_shift_isotope_effect_list_ID'],

            },

      'Isotope_label_pattern': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Chem_shift_isotope_effect_list_ID': [None,returnStarInt,'Chem_shift_isotope_effect_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Chem_shift_isotope_effect_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Chem_shift_isotope_effect_list_ID'],

            }

        },

    'tableNames': ['Isotope_effect_experiment', 'Isotope_effect_software', 'Isotope_effect', 'Isotope_label_pattern']

    },

  'chem_shift_interaction_diff': {

    'name': 'Mol_interaction_chem_shift_diff',

    'tags': {

      'Sf_category': ['chem_shift_interaction_diff',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Chem_shift_ref_set_ID': [None,returnStarInt,'Chem_shift_reference.ID',True],
      'Chem_shift_ref_set_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Chem_shift_ref_set_ID', 'Chem_shift_ref_set_label', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Mol_interaction_diff_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Mol_interaction_chem_shift_diff_ID': [None,returnStarInt,'Mol_interaction_chem_shift_diff.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Mol_interaction_chem_shift_diff_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Mol_interaction_chem_shift_diff_ID'],

            },

      'Mol_interaction_diff_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',True],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Mol_interaction_chem_shift_diff_ID': [None,returnStarInt,'Mol_interaction_chem_shift_diff.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Mol_interaction_chem_shift_diff_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Mol_interaction_chem_shift_diff_ID'],

            },

      'Mol_interaction_chem_shift': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Chem_shift_val': [None,returnStarFloat,None,False],
          'Chem_shift_val_err': [None,returnStarFloat,None,False],
          'Chem_shift_diff_val': [None,returnStarFloat,None,False],
          'Chem_shift_diff_val_err': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Mol_interaction_chem_shift_diff_ID': [None,returnStarInt,'Mol_interaction_chem_shift_diff.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Chem_shift_val', 'Chem_shift_val_err', 'Chem_shift_diff_val', 'Chem_shift_diff_val_err', 'Resonance_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Mol_interaction_chem_shift_diff_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Mol_interaction_chem_shift_diff_ID'],

            }

        },

    'tableNames': ['Mol_interaction_diff_experiment', 'Mol_interaction_diff_software', 'Mol_interaction_chem_shift']

    },

  'chem_shift_anisotropy': {

    'name': 'Chem_shift_anisotropy',

    'tags': {

      'Sf_category': ['chem_shift_anisotropy',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Spectrometer_frequency_1H': [None,returnStarFloat,None,False],
      'Val_units': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Spectrometer_frequency_1H', 'Val_units', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'CS_anisotropy_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Chem_shift_anisotropy_ID': [None,returnStarInt,'Chem_shift_anisotropy.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Chem_shift_anisotropy_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Chem_shift_anisotropy_ID'],

            },

      'CS_anisotropy_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',True],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Chem_shift_anisotropy_ID': [None,returnStarInt,'Chem_shift_anisotropy.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Chem_shift_anisotropy_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Chem_shift_anisotropy_ID'],

            },

      'CS_anisotropy': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,False],
          'Principal_value_sigma_11_val': [None,returnStarFloat,None,False],
          'Principal_value_sigma_22_val': [None,returnStarFloat,None,False],
          'Principal_value_sigma_33_val': [None,returnStarFloat,None,False],
          'Principal_Euler_angle_alpha_val': [None,returnStarFloat,None,False],
          'Principal_Euler_angle_beta_val': [None,returnStarFloat,None,False],
          'Principal_Euler_angle_gamma_val': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Chem_shift_anisotropy_ID': [None,returnStarInt,'Chem_shift_anisotropy.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Val', 'Val_err', 'Principal_value_sigma_11_val', 'Principal_value_sigma_22_val', 'Principal_value_sigma_33_val', 'Principal_Euler_angle_alpha_val', 'Principal_Euler_angle_beta_val', 'Principal_Euler_angle_gamma_val', 'Resonance_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Chem_shift_anisotropy_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Chem_shift_anisotropy_ID'],

            }

        },

    'tableNames': ['CS_anisotropy_experiment', 'CS_anisotropy_software', 'CS_anisotropy']

    },

  'chem_shifts_calc_type': {

    'name': 'Chem_shifts_calc_type',

    'tags': {

      'Sf_category': ['chem_shifts_calc_type',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Calculation_level': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Quantum_mechanical_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Quantum_mechanical_theory_level': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Quantum_mechanical_basis_set': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Chem_shift_nucleus': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Modeled_sample_cond_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Modeled_sample_cond_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Chem_shift_reference_ID': [None,returnStarInt,'Chem_shift_reference.ID',True],
      'Chem_shift_reference_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Calculation_level', 'Quantum_mechanical_method', 'Quantum_mechanical_theory_level', 'Quantum_mechanical_basis_set', 'Chem_shift_nucleus', 'Modeled_sample_cond_list_ID', 'Modeled_sample_cond_list_label', 'Chem_shift_reference_ID', 'Chem_shift_reference_label', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Chem_shifts_calc_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Chem_shifts_calc_type_ID': [None,returnStarInt,'Chem_shifts_calc_type.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Chem_shifts_calc_type_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Chem_shifts_calc_type_ID'],

            }

        },

    'tableNames': ['Chem_shifts_calc_software']

    },

  'shielding_tensors': {

    'name': 'Shielding_tensor_list',

    'tags': {

      'Sf_category': ['shielding_tensors',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Shielding_tensor': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Shieldings_calc_type_ID': [None,returnStarInt,'Chem_shifts_calc_type.ID',True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Isotropic_comp_1_1_val': [None,returnStarFloat,None,False],
          'Isotropic_comp_2_2_val': [None,returnStarFloat,None,False],
          'Isotropic_comp_3_3_val': [None,returnStarFloat,None,False],
          'Anti_sym_comp_1_2_val': [None,returnStarFloat,None,False],
          'Anti_sym_comp_1_3_val': [None,returnStarFloat,None,False],
          'Anti_sym_comp_2_3_val': [None,returnStarFloat,None,False],
          'Sym_traceless_comp_1_1_val': [None,returnStarFloat,None,False],
          'Sym_traceless_comp_1_2_val': [None,returnStarFloat,None,False],
          'Sym_traceless_comp_1_3_val': [None,returnStarFloat,None,False],
          'Sym_traceless_comp_2_2_val': [None,returnStarFloat,None,False],
          'Sym_traceless_comp_2_3_val': [None,returnStarFloat,None,False],
          'Reduceable_matrix_1_1_val': [None,returnStarFloat,None,False],
          'Reduceable_matrix_1_2_val': [None,returnStarFloat,None,False],
          'Reduceable_matrix_1_3_val': [None,returnStarFloat,None,False],
          'Reduceable_matrix_2_1_val': [None,returnStarFloat,None,False],
          'Reduceable_matrix_2_2_val': [None,returnStarFloat,None,False],
          'Reduceable_matrix_2_3_val': [None,returnStarFloat,None,False],
          'Reduceable_matrix_3_1_val': [None,returnStarFloat,None,False],
          'Reduceable_matrix_3_2_val': [None,returnStarFloat,None,False],
          'Reduceable_matrix_3_3_val': [None,returnStarFloat,None,False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Shielding_tensor_list_ID': [None,returnStarInt,'Shielding_tensor_list.ID',True],

                },

        'tagNames': ['ID', 'Shieldings_calc_type_ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Isotropic_comp_1_1_val', 'Isotropic_comp_2_2_val', 'Isotropic_comp_3_3_val', 'Anti_sym_comp_1_2_val', 'Anti_sym_comp_1_3_val', 'Anti_sym_comp_2_3_val', 'Sym_traceless_comp_1_1_val', 'Sym_traceless_comp_1_2_val', 'Sym_traceless_comp_1_3_val', 'Sym_traceless_comp_2_2_val', 'Sym_traceless_comp_2_3_val', 'Reduceable_matrix_1_1_val', 'Reduceable_matrix_1_2_val', 'Reduceable_matrix_1_3_val', 'Reduceable_matrix_2_1_val', 'Reduceable_matrix_2_2_val', 'Reduceable_matrix_2_3_val', 'Reduceable_matrix_3_1_val', 'Reduceable_matrix_3_2_val', 'Reduceable_matrix_3_3_val', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Shielding_tensor_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Shielding_tensor_list_ID'],

            }

        },

    'tableNames': ['Shielding_tensor']

    },

  'theoretical_chem_shifts': {

    'name': 'Theoretical_chem_shift_list',

    'tags': {

      'Sf_category': ['theoretical_chem_shifts',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Chem_shifts_calc_type_ID': [None,returnStarInt,'Chem_shifts_calc_type.ID',True],
      'Chem_shifts_calc_type_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Model_atomic_coordinates_ID': [None,returnStarInt,'Representative_conformer.ID',True],
      'Model_atomic_coordinates_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Chem_shielding_tensor_list_ID': [None,returnStarInt,'Shielding_tensor_list.ID',True],
      'Chem_shielding_tensor_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Fermi_contact_spin_density_units': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Chem_shift_1H_err': [None,returnStarFloat,None,False],
      'Chem_shift_2H_err': [None,returnStarFloat,None,False],
      'Chem_shift_13C_err': [None,returnStarFloat,None,False],
      'Chem_shift_15N_err': [None,returnStarFloat,None,False],
      'Chem_shift_19F_err': [None,returnStarFloat,None,False],
      'Chem_shift_31P_err': [None,returnStarFloat,None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Chem_shifts_calc_type_ID', 'Chem_shifts_calc_type_label', 'Model_atomic_coordinates_ID', 'Model_atomic_coordinates_label', 'Chem_shielding_tensor_list_ID', 'Chem_shielding_tensor_list_label', 'Fermi_contact_spin_density_units', 'Chem_shift_1H_err', 'Chem_shift_2H_err', 'Chem_shift_13C_err', 'Chem_shift_15N_err', 'Chem_shift_19F_err', 'Chem_shift_31P_err', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Theoretical_chem_shift': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Fermi_contact_spin_density': [None,lambda x = value: returnStarString(x,length = 127),None,False],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Theoretical_chem_shift_list_ID': [None,returnStarInt,'Theoretical_chem_shift_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Fermi_contact_spin_density', 'Val', 'Val_err', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Theoretical_chem_shift_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Theoretical_chem_shift_list_ID'],

            }

        },

    'tableNames': ['Theoretical_chem_shift']

    },

  'RDCs': {

    'name': 'RDC_list',

    'tags': {

      'Sf_category': ['RDCs',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Spectrometer_frequency_1H': [None,returnStarFloat,None,True],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Spectrometer_frequency_1H', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'RDC_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'RDC_list_ID': [None,returnStarInt,'RDC_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'RDC_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'RDC_list_ID'],

            },

      'RDC_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'RDC_list_ID': [None,returnStarInt,'RDC_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'RDC_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'RDC_list_ID'],

            },

      'RDC': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'RDC_code': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_1': [None,returnStarInt,None,False],
          'Ambiguity_code_1': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_2': [None,returnStarInt,None,False],
          'Ambiguity_code_2': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
          'Val': [None,returnStarFloat,None,False],
          'Val_min': [None,returnStarFloat,None,False],
          'Val_max': [None,returnStarFloat,None,False],
          'Val_err': [None,returnStarFloat,None,False],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'RDC_list_ID': [None,returnStarInt,'RDC_list.ID',True],

                },

        'tagNames': ['ID', 'RDC_code', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Atom_isotope_number_1', 'Ambiguity_code_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Atom_isotope_number_2', 'Ambiguity_code_2', 'Val', 'Val_min', 'Val_max', 'Val_err', 'Resonance_ID_1', 'Resonance_ID_2', 'Auth_entity_assembly_ID_1', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_entity_assembly_ID_2', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Entry_ID', 'RDC_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'RDC_list_ID'],

            }

        },

    'tableNames': ['RDC_experiment', 'RDC_software', 'RDC']

    },

  'dipolar_couplings': {

    'name': 'Dipolar_coupling_list',

    'tags': {

      'Sf_category': ['dipolar_couplings',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Spectrometer_frequency_1H': [None,returnStarFloat,None,True],
      'Scaling_factor': [None,returnStarFloat,None,False],
      'Fitting_procedure': [None,returnStarString,None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Spectrometer_frequency_1H', 'Scaling_factor', 'Fitting_procedure', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Dipolar_coupling_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Dipolar_coupling_list_ID': [None,returnStarInt,'Dipolar_coupling_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Dipolar_coupling_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Dipolar_coupling_list_ID'],

            },

      'Dipolar_coupling_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Dipolar_coupling_list_ID': [None,returnStarInt,'Dipolar_coupling_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Dipolar_coupling_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Dipolar_coupling_list_ID'],

            },

      'Dipolar_coupling': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Dipolar_coupling_code': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_1': [None,returnStarInt,None,False],
          'Ambiguity_code_1': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_2': [None,returnStarInt,None,False],
          'Ambiguity_code_2': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
          'Val': [None,returnStarFloat,None,False],
          'Val_min': [None,returnStarFloat,None,False],
          'Val_max': [None,returnStarFloat,None,False],
          'Val_err': [None,returnStarFloat,None,False],
          'Principal_Euler_angle_alpha_val': [None,returnStarFloat,None,False],
          'Principal_Euler_angle_beta_val': [None,returnStarFloat,None,False],
          'Principal_Euler_angle_gamma_val': [None,returnStarFloat,None,False],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Dipolar_coupling_list_ID': [None,returnStarInt,'Dipolar_coupling_list.ID',True],

                },

        'tagNames': ['ID', 'Dipolar_coupling_code', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Atom_isotope_number_1', 'Ambiguity_code_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Atom_isotope_number_2', 'Ambiguity_code_2', 'Val', 'Val_min', 'Val_max', 'Val_err', 'Principal_Euler_angle_alpha_val', 'Principal_Euler_angle_beta_val', 'Principal_Euler_angle_gamma_val', 'Resonance_ID_1', 'Resonance_ID_2', 'Auth_entity_assembly_ID_1', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_entity_assembly_ID_2', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Entry_ID', 'Dipolar_coupling_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Dipolar_coupling_list_ID'],

            }

        },

    'tableNames': ['Dipolar_coupling_experiment', 'Dipolar_coupling_software', 'Dipolar_coupling']

    },

  'spectral_density_values': {

    'name': 'Spectral_density_list',

    'tags': {

      'Sf_category': ['spectral_density_values',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Spectral_density_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_density_list_ID': [None,returnStarInt,'Spectral_density_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Spectral_density_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Spectral_density_list_ID'],

            },

      'Spectral_density_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',True],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_density_list_ID': [None,returnStarInt,'Spectral_density_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Spectral_density_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Spectral_density_list_ID'],

            },

      'Spectral_density': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'W_zero_val': [None,returnStarFloat,None,False],
          'W_zero_val_err': [None,returnStarFloat,None,False],
          'W_1H_val': [None,returnStarFloat,None,False],
          'W_1H_val_err': [None,returnStarFloat,None,False],
          'W_13C_val': [None,returnStarFloat,None,False],
          'W_13C_val_err': [None,returnStarFloat,None,False],
          'W_15N_val': [None,returnStarFloat,None,False],
          'W_15N_val_err': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Spectral_density_list_ID': [None,returnStarInt,'Spectral_density_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'W_zero_val', 'W_zero_val_err', 'W_1H_val', 'W_1H_val_err', 'W_13C_val', 'W_13C_val_err', 'W_15N_val', 'W_15N_val_err', 'Resonance_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Spectral_density_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Spectral_density_list_ID'],

            }

        },

    'tableNames': ['Spectral_density_experiment', 'Spectral_density_software', 'Spectral_density']

    },

  'other_data_types': {

    'name': 'Other_data_type_list',

    'tags': {

      'Sf_category': ['other_data_types',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Definition': [None,returnStarString,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Name', 'Definition', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Other_data_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Other_data_type_list_ID': [None,returnStarInt,'Other_data_type_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Other_data_type_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Other_data_type_list_ID'],

            },

      'Other_data_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Other_data_type_list_ID': [None,returnStarInt,'Other_data_type_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Other_data_type_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Other_data_type_list_ID'],

            },

      'Other_data': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Other_data_type_list_ID': [None,returnStarInt,'Other_data_type_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Val', 'Val_err', 'Resonance_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Other_data_type_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Other_data_type_list_ID'],

            }

        },

    'tableNames': ['Other_data_experiment', 'Other_data_software', 'Other_data']

    },

  'H_exch_rates': {

    'name': 'H_exch_rate_list',

    'tags': {

      'Sf_category': ['H_exch_rates',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Val_units': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Val_units', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'H_exch_rate_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'H_exch_rate_list_ID': [None,returnStarInt,'H_exch_rate_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'H_exch_rate_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'H_exch_rate_list_ID'],

            },

      'H_exch_rate_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'H_exch_rate_list_ID': [None,returnStarInt,'H_exch_rate_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'H_exch_rate_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'H_exch_rate_list_ID'],

            },

      'H_exch_rate': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Val': [None,returnStarFloat,None,False],
          'Val_min': [None,returnStarFloat,None,False],
          'Val_max': [None,returnStarFloat,None,False],
          'Val_err': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'H_exch_rate_list_ID': [None,returnStarInt,'H_exch_rate_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Val', 'Val_min', 'Val_max', 'Val_err', 'Resonance_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'H_exch_rate_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'H_exch_rate_list_ID'],

            }

        },

    'tableNames': ['H_exch_rate_experiment', 'H_exch_rate_software', 'H_exch_rate']

    },

  'H_exch_protection_factors': {

    'name': 'H_exch_protection_factor_list',

    'tags': {

      'Sf_category': ['H_exch_protection_factors',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Std_values_source_cit_ID': [None,returnStarInt,'Citation.ID',False],
      'Std_values_source_cit_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Std_values_source_cit_ID', 'Std_values_source_cit_label', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'H_exch_protection_fact_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'H_exch_protection_factor_list_ID': [None,returnStarInt,'H_exch_protection_factor_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'H_exch_protection_factor_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'H_exch_protection_factor_list_ID'],

            },

      'H_exch_protect_fact_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'H_exch_protection_factor_list_ID': [None,returnStarInt,'H_exch_protection_factor_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'H_exch_protection_factor_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'H_exch_protection_factor_list_ID'],

            },

      'H_exch_protection_factor': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Calculated_intrinsic_rate': [None,returnStarFloat,None,False],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'H_exch_protection_factor_list_ID': [None,returnStarInt,'H_exch_protection_factor_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Calculated_intrinsic_rate', 'Val', 'Val_err', 'Resonance_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'H_exch_protection_factor_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'H_exch_protection_factor_list_ID'],

            }

        },

    'tableNames': ['H_exch_protection_fact_experiment', 'H_exch_protect_fact_software', 'H_exch_protection_factor']

    },

  'homonucl_NOEs': {

    'name': 'Homonucl_NOE_list',

    'tags': {

      'Sf_category': ['homonucl_NOEs',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Homonuclear_NOE_val_type': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'NOE_ref_val': [None,returnStarFloat,None,False],
      'NOE_ref_description': [None,returnStarString,None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Homonuclear_NOE_val_type', 'NOE_ref_val', 'NOE_ref_description', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Homonucl_NOE_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Homonucl_NOE_list_ID': [None,returnStarInt,'Homonucl_NOE_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Homonucl_NOE_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Homonucl_NOE_list_ID'],

            },

      'Homonucl_NOE_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Homonucl_NOE_list_ID': [None,returnStarInt,'Homonucl_NOE_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Homonucl_NOE_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Homonucl_NOE_list_ID'],

            },

      'Homonucl_NOE': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_1': [None,returnStarInt,None,False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_2': [None,returnStarInt,None,False],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,False],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Homonucl_NOE_list_ID': [None,returnStarInt,'Homonucl_NOE_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Atom_isotope_number_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Atom_isotope_number_2', 'Val', 'Val_err', 'Resonance_ID_1', 'Resonance_ID_2', 'Auth_entity_assembly_ID_1', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_entity_assembly_ID_2', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Entry_ID', 'Homonucl_NOE_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Homonucl_NOE_list_ID'],

            }

        },

    'tableNames': ['Homonucl_NOE_experiment', 'Homonucl_NOE_software', 'Homonucl_NOE']

    },

  'heteronucl_NOEs': {

    'name': 'Heteronucl_NOE_list',

    'tags': {

      'Sf_category': ['heteronucl_NOEs',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Spectrometer_frequency_1H': [None,returnStarFloat,None,True],
      'Heteronuclear_NOE_val_type': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'NOE_ref_val': [None,returnStarFloat,None,False],
      'NOE_ref_description': [None,returnStarString,None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Spectrometer_frequency_1H', 'Heteronuclear_NOE_val_type', 'NOE_ref_val', 'NOE_ref_description', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Heteronucl_NOE_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Heteronucl_NOE_list_ID': [None,returnStarInt,'Heteronucl_NOE_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Heteronucl_NOE_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Heteronucl_NOE_list_ID'],

            },

      'Heteronucl_NOE_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Heteronucl_NOE_list_ID': [None,returnStarInt,'Heteronucl_NOE_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Heteronucl_NOE_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Heteronucl_NOE_list_ID'],

            },

      'Heteronucl_NOE': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_1': [None,returnStarInt,None,False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_2': [None,returnStarInt,None,False],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,False],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Heteronucl_NOE_list_ID': [None,returnStarInt,'Heteronucl_NOE_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Atom_isotope_number_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Atom_isotope_number_2', 'Val', 'Val_err', 'Resonance_ID_1', 'Resonance_ID_2', 'Auth_entity_assembly_ID_1', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_entity_assembly_ID_2', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Entry_ID', 'Heteronucl_NOE_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Heteronucl_NOE_list_ID'],

            }

        },

    'tableNames': ['Heteronucl_NOE_experiment', 'Heteronucl_NOE_software', 'Heteronucl_NOE']

    },

  'heteronucl_T1_relaxation': {

    'name': 'Heteronucl_T1_list',

    'tags': {

      'Sf_category': ['heteronucl_T1_relaxation',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Spectrometer_frequency_1H': [None,returnStarFloat,None,True],
      'T1_coherence_type': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'T1_val_units': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Spectrometer_frequency_1H', 'T1_coherence_type', 'T1_val_units', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Heteronucl_T1_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Heteronucl_T1_list_ID': [None,returnStarInt,'Heteronucl_T1_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Heteronucl_T1_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Heteronucl_T1_list_ID'],

            },

      'Heteronucl_T1_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Heteronucl_T1_list_ID': [None,returnStarInt,'Heteronucl_T1_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Heteronucl_T1_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Heteronucl_T1_list_ID'],

            },

      'T1': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Heteronucl_T1_list_ID': [None,returnStarInt,'Heteronucl_T1_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Val', 'Val_err', 'Resonance_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Heteronucl_T1_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Heteronucl_T1_list_ID'],

            }

        },

    'tableNames': ['Heteronucl_T1_experiment', 'Heteronucl_T1_software', 'T1']

    },

  'heteronucl_T1rho_relaxation': {

    'name': 'Heteronucl_T1rho_list',

    'tags': {

      'Sf_category': ['heteronucl_T1rho_relaxation',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Spectrometer_frequency_1H': [None,returnStarFloat,None,True],
      'T1rho_coherence_type': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'T1rho_val_units': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'Rex_units': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Spectrometer_frequency_1H', 'T1rho_coherence_type', 'T1rho_val_units', 'Rex_units', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Heteronucl_T1rho_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Heteronucl_T1rho_list_ID': [None,returnStarInt,'Heteronucl_T1rho_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Heteronucl_T1rho_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Heteronucl_T1rho_list_ID'],

            },

      'Heteronucl_T1rho_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Heteronucl_T1rho_list_ID': [None,returnStarInt,'Heteronucl_T1rho_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Heteronucl_T1rho_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Heteronucl_T1rho_list_ID'],

            },

      'T1rho': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'T1rho_val': [None,returnStarFloat,None,True],
          'T1rho_val_err': [None,returnStarFloat,None,False],
          'Rex_val': [None,returnStarFloat,None,False],
          'Rex_val_err': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Heteronucl_T1rho_list_ID': [None,returnStarInt,'Heteronucl_T1rho_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'T1rho_val', 'T1rho_val_err', 'Rex_val', 'Rex_val_err', 'Resonance_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Heteronucl_T1rho_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Heteronucl_T1rho_list_ID'],

            }

        },

    'tableNames': ['Heteronucl_T1rho_experiment', 'Heteronucl_T1rho_software', 'T1rho']

    },

  'heteronucl_T2_relaxation': {

    'name': 'Heteronucl_T2_list',

    'tags': {

      'Sf_category': ['heteronucl_T2_relaxation',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Spectrometer_frequency_1H': [None,returnStarFloat,None,True],
      'T2_coherence_type': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'T2_val_units': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'Rex_units': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Spectrometer_frequency_1H', 'T2_coherence_type', 'T2_val_units', 'Rex_units', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Heteronucl_T2_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Heteronucl_T2_list_ID': [None,returnStarInt,'Heteronucl_T2_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Heteronucl_T2_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Heteronucl_T2_list_ID'],

            },

      'Heteronucl_T2_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Heteronucl_T2_list_ID': [None,returnStarInt,'Heteronucl_T2_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Heteronucl_T2_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Heteronucl_T2_list_ID'],

            },

      'T2': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'T2_val': [None,returnStarFloat,None,True],
          'T2_val_err': [None,returnStarFloat,None,False],
          'Rex_val': [None,returnStarFloat,None,False],
          'Rex_err': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Heteronucl_T2_list_ID': [None,returnStarInt,'Heteronucl_T2_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'T2_val', 'T2_val_err', 'Rex_val', 'Rex_err', 'Resonance_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Heteronucl_T2_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Heteronucl_T2_list_ID'],

            }

        },

    'tableNames': ['Heteronucl_T2_experiment', 'Heteronucl_T2_software', 'T2']

    },

  'dipole_dipole_relaxation': {

    'name': 'Dipole_dipole_relax_list',

    'tags': {

      'Sf_category': ['dipole_dipole_relaxation',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Spectrometer_frequency_1H': [None,returnStarFloat,None,True],
      'Val_units': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Spectrometer_frequency_1H', 'Val_units', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Dipole_dipole_relax_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Dipole_dipole_relax_list_ID': [None,returnStarInt,'Dipole_dipole_relax_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Dipole_dipole_relax_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Dipole_dipole_relax_list_ID'],

            },

      'Dipole_dipole_relax_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Dipole_dipole_relax_list_ID': [None,returnStarInt,'Dipole_dipole_relax_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Dipole_dipole_relax_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Dipole_dipole_relax_list_ID'],

            },

      'Dipole_dipole_relax': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_1': [None,returnStarInt,None,False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_2': [None,returnStarInt,None,False],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,False],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Dipole_dipole_relax_list_ID': [None,returnStarInt,'Dipole_dipole_relax_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Atom_isotope_number_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Atom_isotope_number_2', 'Val', 'Val_err', 'Resonance_ID_1', 'Resonance_ID_2', 'Auth_entity_assembly_ID_1', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_entity_assembly_ID_2', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Entry_ID', 'Dipole_dipole_relax_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Dipole_dipole_relax_list_ID'],

            }

        },

    'tableNames': ['Dipole_dipole_relax_experiment', 'Dipole_dipole_relax_software', 'Dipole_dipole_relax']

    },

  'dipole_dipole_cross_correlations': {

    'name': 'Cross_correlation_DD_list',

    'tags': {

      'Sf_category': ['dipole_dipole_cross_correlations',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Spectrometer_frequency_1H': [None,returnStarFloat,None,True],
      'Val_units': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Spectrometer_frequency_1H', 'Val_units', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Cross_correlation_DD_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Cross_correlation_list_ID': [None,returnStarInt,'Cross_correlation_DD_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Cross_correlation_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Cross_correlation_list_ID'],

            },

      'Cross_correlation_DD_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Cross_correlation_list_ID': [None,returnStarInt,'Cross_correlation_DD_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Cross_correlation_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Cross_correlation_list_ID'],

            },

      'Cross_correlation_DD': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Dipole_1_entry_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Dipole_1_entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',False],
          'Dipole_1_entity_ID_1': [None,returnStarInt,'Entity.ID',False],
          'Dipole_1_comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',False],
          'Dipole_1_seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Dipole_1_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'Dipole_1_atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
          'Dipole_1_atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_1_atom_isotope_number_1': [None,returnStarInt,None,False],
          'Dipole_1_entry_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Dipole_1_entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',False],
          'Dipole_1_entity_ID_2': [None,returnStarInt,'Entity.ID',False],
          'Dipole_1_comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',False],
          'Dipole_1_seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Dipole_1_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'Dipole_1_atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
          'Dipole_1_atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_1_atom_isotope_number_2': [None,returnStarInt,None,False],
          'Dipole_2_entry_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Dipole_2_entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',False],
          'Dipole_2_entity_ID_1': [None,returnStarInt,'Entity.ID',False],
          'Dipole_2_comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',False],
          'Dipole_2_seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Dipole_2_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'Dipole_2_atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
          'Dipole_2_atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_2_atom_isotope_number_1': [None,returnStarInt,None,False],
          'Dipole_2_entry_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Dipole_2_entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',False],
          'Dipole_2_entity_ID_2': [None,returnStarInt,'Entity.ID',False],
          'Dipole_2_chem_comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',False],
          'Dipole_2_seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Dipole_2_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'Dipole_2_atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
          'Dipole_2_atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_2_atom_isotope_number_2': [None,returnStarInt,None,False],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Dipole_1_auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_1_auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_1_auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_1_auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_1_auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_1_auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_1_auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_1_auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_2_auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_2_auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_2_auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_2_auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_2_auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_2_auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_2_auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_2_auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Cross_correlation_list_ID': [None,returnStarInt,'Cross_correlation_DD_list.ID',True],

                },

        'tagNames': ['ID', 'Dipole_1_entry_atom_ID_1', 'Dipole_1_entity_assembly_ID_1', 'Dipole_1_entity_ID_1', 'Dipole_1_comp_index_ID_1', 'Dipole_1_seq_ID_1', 'Dipole_1_comp_ID_1', 'Dipole_1_atom_ID_1', 'Dipole_1_atom_type_1', 'Dipole_1_atom_isotope_number_1', 'Dipole_1_entry_atom_ID_2', 'Dipole_1_entity_assembly_ID_2', 'Dipole_1_entity_ID_2', 'Dipole_1_comp_index_ID_2', 'Dipole_1_seq_ID_2', 'Dipole_1_comp_ID_2', 'Dipole_1_atom_ID_2', 'Dipole_1_atom_type_2', 'Dipole_1_atom_isotope_number_2', 'Dipole_2_entry_atom_ID_1', 'Dipole_2_entity_assembly_ID_1', 'Dipole_2_entity_ID_1', 'Dipole_2_comp_index_ID_1', 'Dipole_2_seq_ID_1', 'Dipole_2_comp_ID_1', 'Dipole_2_atom_ID_1', 'Dipole_2_atom_type_1', 'Dipole_2_atom_isotope_number_1', 'Dipole_2_entry_atom_ID_2', 'Dipole_2_entity_assembly_ID_2', 'Dipole_2_entity_ID_2', 'Dipole_2_chem_comp_index_ID_2', 'Dipole_2_seq_ID_2', 'Dipole_2_comp_ID_2', 'Dipole_2_atom_ID_2', 'Dipole_2_atom_type_2', 'Dipole_2_atom_isotope_number_2', 'Val', 'Val_err', 'Resonance_ID', 'Dipole_1_auth_entity_assembly_ID_1', 'Dipole_1_auth_seq_ID_1', 'Dipole_1_auth_comp_ID_1', 'Dipole_1_auth_atom_ID_1', 'Dipole_1_auth_entity_assembly_ID_2', 'Dipole_1_auth_seq_ID_2', 'Dipole_1_auth_comp_ID_2', 'Dipole_1_auth_atom_ID_2', 'Dipole_2_auth_entity_assembly_ID_1', 'Dipole_2_auth_seq_ID_1', 'Dipole_2_auth_comp_ID_1', 'Dipole_2_auth_atom_ID_1', 'Dipole_2_auth_entity_assembly_ID_2', 'Dipole_2_auth_seq_ID_2', 'Dipole_2_auth_comp_ID_2', 'Dipole_2_auth_atom_ID_2', 'Entry_ID', 'Cross_correlation_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Cross_correlation_list_ID'],

            }

        },

    'tableNames': ['Cross_correlation_DD_experiment', 'Cross_correlation_DD_software', 'Cross_correlation_DD']

    },

  'dipole_CSA_cross_correlations': {

    'name': 'Cross_correlation_D_CSA_list',

    'tags': {

      'Sf_category': ['dipole_CSA_cross_correlations',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Spectrometer_frequency_1H': [None,returnStarFloat,None,True],
      'Val_units': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Spectrometer_frequency_1H', 'Val_units', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Cross_correlation_D_CSA_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Cross_correlation_list_ID': [None,returnStarInt,'Cross_correlation_D_CSA_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Cross_correlation_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Cross_correlation_list_ID'],

            },

      'Cross_correlation_D_CSA_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Cross_correlation_list_ID': [None,returnStarInt,'Cross_correlation_D_CSA_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Cross_correlation_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Cross_correlation_list_ID'],

            },

      'Cross_correlation_D_CSA': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Dipole_entry_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Dipole_entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',False],
          'Dipole_entity_ID_1': [None,returnStarInt,'Entity.ID',False],
          'Dipole_comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',False],
          'Dipole_seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Dipole_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'Dipole_atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
          'Dipole_atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_atom_isotope_number_1': [None,returnStarInt,None,False],
          'Dipole_entry_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Dipole_entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',False],
          'Dipole_entity_ID_2': [None,returnStarInt,'Entity.ID',False],
          'Dipole_comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',False],
          'Dipole_seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Dipole_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'Dipole_atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
          'Dipole_atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_atom_isotope_number_2': [None,returnStarInt,None,False],
          'CSA_entry_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'CSA_entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',False],
          'CSA_entity_ID_1': [None,returnStarInt,'Entity.ID',False],
          'CSA_comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',False],
          'CSA_seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'CSA_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'CSA_atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
          'CSA_atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'CSA_atom_isotope_number_1': [None,returnStarInt,None,False],
          'CSA_entry_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'CSA_entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',False],
          'CSA_entity_ID_2': [None,returnStarInt,'Entity.ID',False],
          'CSA_comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',False],
          'CSA_seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'CSA_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
          'CSA_atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
          'CSA_atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'CSA_atom_isotope_number_2': [None,returnStarInt,None,False],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,False],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Dipole_auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Dipole_auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'CSA_auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'CSA_auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'CSA_auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'CSA_auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'CSA_auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'CSA_auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'CSA_auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'CSA_auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Cross_correlation_list_ID': [None,returnStarInt,'Cross_correlation_D_CSA_list.ID',True],

                },

        'tagNames': ['ID', 'Dipole_entry_atom_ID_1', 'Dipole_entity_assembly_ID_1', 'Dipole_entity_ID_1', 'Dipole_comp_index_ID_1', 'Dipole_seq_ID_1', 'Dipole_comp_ID_1', 'Dipole_atom_ID_1', 'Dipole_atom_type_1', 'Dipole_atom_isotope_number_1', 'Dipole_entry_atom_ID_2', 'Dipole_entity_assembly_ID_2', 'Dipole_entity_ID_2', 'Dipole_comp_index_ID_2', 'Dipole_seq_ID_2', 'Dipole_comp_ID_2', 'Dipole_atom_ID_2', 'Dipole_atom_type_2', 'Dipole_atom_isotope_number_2', 'CSA_entry_atom_ID_1', 'CSA_entity_assembly_ID_1', 'CSA_entity_ID_1', 'CSA_comp_index_ID_1', 'CSA_seq_ID_1', 'CSA_comp_ID_1', 'CSA_atom_ID_1', 'CSA_atom_type_1', 'CSA_atom_isotope_number_1', 'CSA_entry_atom_ID_2', 'CSA_entity_assembly_ID_2', 'CSA_entity_ID_2', 'CSA_comp_index_ID_2', 'CSA_seq_ID_2', 'CSA_comp_ID_2', 'CSA_atom_ID_2', 'CSA_atom_type_2', 'CSA_atom_isotope_number_2', 'Val', 'Val_err', 'Resonance_ID_1', 'Resonance_ID_2', 'Dipole_auth_entity_assembly_ID_1', 'Dipole_auth_seq_ID_1', 'Dipole_auth_comp_ID_1', 'Dipole_auth_atom_ID_1', 'Dipole_auth_entity_assembly_ID_2', 'Dipole_auth_seq_ID_2', 'Dipole_auth_comp_ID_2', 'Dipole_auth_atom_ID_2', 'CSA_auth_entity_assembly_ID_1', 'CSA_auth_seq_ID_1', 'CSA_auth_comp_ID_1', 'CSA_auth_atom_ID_1', 'CSA_auth_entity_assembly_ID_2', 'CSA_auth_seq_ID_2', 'CSA_auth_comp_ID_2', 'CSA_auth_atom_ID_2', 'Entry_ID', 'Cross_correlation_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Cross_correlation_list_ID'],

            }

        },

    'tableNames': ['Cross_correlation_D_CSA_experiment', 'Cross_correlation_D_CSA_software', 'Cross_correlation_D_CSA']

    },

  'order_parameters': {

    'name': 'Order_parameter_list',

    'tags': {

      'Sf_category': ['order_parameters',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Tau_e_val_units': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'Tau_s_val_units': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Tau_e_val_units', 'Tau_s_val_units', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Order_parameter_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Order_parameter_list_ID': [None,returnStarInt,'Order_parameter_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Order_parameter_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Order_parameter_list_ID'],

            },

      'Order_parameter_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Order_parameter_list_ID': [None,returnStarInt,'Order_parameter_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Order_parameter_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Order_parameter_list_ID'],

            },

      'Order_param': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Order_param_val': [None,returnStarFloat,None,True],
          'Order_param_val_fit_err': [None,returnStarFloat,None,False],
          'Tau_e_val': [None,returnStarFloat,None,False],
          'Tau_e_val_fit_err': [None,returnStarFloat,None,False],
          'Rex_val': [None,returnStarFloat,None,False],
          'Rex_val_fit_err': [None,returnStarFloat,None,False],
          'Model_free_sum_squared_errs': [None,returnStarFloat,None,False],
          'Model_fit': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Sf2_val': [None,returnStarFloat,None,False],
          'Sf2_val_fit_err': [None,returnStarFloat,None,False],
          'Ss2_val': [None,returnStarFloat,None,False],
          'Ss2_val_fit_err': [None,returnStarFloat,None,False],
          'Tau_s_val': [None,returnStarFloat,None,False],
          'Tau_s_val_fit_err': [None,returnStarFloat,None,False],
          'SH2_val': [None,returnStarFloat,None,False],
          'SH2_val_fit_err': [None,returnStarFloat,None,False],
          'SN2_val': [None,returnStarFloat,None,False],
          'SN2_val_fit_err': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Order_parameter_list_ID': [None,returnStarInt,'Order_parameter_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Order_param_val', 'Order_param_val_fit_err', 'Tau_e_val', 'Tau_e_val_fit_err', 'Rex_val', 'Rex_val_fit_err', 'Model_free_sum_squared_errs', 'Model_fit', 'Sf2_val', 'Sf2_val_fit_err', 'Ss2_val', 'Ss2_val_fit_err', 'Tau_s_val', 'Tau_s_val_fit_err', 'SH2_val', 'SH2_val_fit_err', 'SN2_val', 'SN2_val_fit_err', 'Resonance_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Order_parameter_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Order_parameter_list_ID'],

            }

        },

    'tableNames': ['Order_parameter_experiment', 'Order_parameter_software', 'Order_param']

    },

  'pH_titration': {

    'name': 'PH_titration_list',

    'tags': {

      'Sf_category': ['pH_titration',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Expt_observed_param': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Expt_observed_param', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'PH_titration_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'PH_titration_list_ID': [None,returnStarInt,'PH_titration_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'PH_titration_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'PH_titration_list_ID'],

            },

      'PH_titration_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'PH_titration_list_ID': [None,returnStarInt,'PH_titration_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'PH_titration_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'PH_titration_list_ID'],

            },

      'PH_titr_result': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Atm_obs_entry_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Atm_obs_entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Atm_obs_entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Atm_obs_comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Atm_obs_seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atm_obs_comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atm_obs_atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atm_obs_atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atm_obs_atom_isotope_number': [None,returnStarInt,None,False],
          'Atm_obs_auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atm_obs_auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atm_obs_auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atm_obs_auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atm_titr_entry_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Atm_titr_entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Atm_titr_entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Atm_titr_comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Atm_titr_seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atm_titr_comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atm_titr_atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atm_titr_atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atm_titr_atom_isotope_number': [None,returnStarInt,None,False],
          'Atm_titr_auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atm_titr_auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atm_titr_auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Atm_titr_auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Hill_coeff_val': [None,returnStarFloat,None,False],
          'Hill_coeff_val_fit_err': [None,returnStarFloat,None,False],
          'High_PH_param_fit_val': [None,returnStarFloat,None,False],
          'High_PH_param_fit_val_err': [None,returnStarFloat,None,False],
          'Low_PH_param_fit_val': [None,returnStarFloat,None,False],
          'Low_PH_param_fit_val_err': [None,returnStarFloat,None,False],
          'PKa_val': [None,returnStarFloat,None,False],
          'PKa_val_fit_err': [None,returnStarFloat,None,False],
          'PHmid_val': [None,returnStarFloat,None,False],
          'PHmid_val_fit_err': [None,returnStarFloat,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'PH_titration_list_ID': [None,returnStarInt,'PH_titration_list.ID',True],

                },

        'tagNames': ['ID', 'Atm_obs_entry_atom_ID', 'Atm_obs_entity_assembly_ID', 'Atm_obs_entity_ID', 'Atm_obs_comp_index_ID', 'Atm_obs_seq_ID', 'Atm_obs_comp_ID', 'Atm_obs_atom_ID', 'Atm_obs_atom_type', 'Atm_obs_atom_isotope_number', 'Atm_obs_auth_entity_assembly_ID', 'Atm_obs_auth_seq_ID', 'Atm_obs_auth_comp_ID', 'Atm_obs_auth_atom_ID', 'Atm_titr_entry_atom_ID', 'Atm_titr_entity_assembly_ID', 'Atm_titr_entity_ID', 'Atm_titr_comp_index_ID', 'Atm_titr_seq_ID', 'Atm_titr_comp_ID', 'Atm_titr_atom_ID', 'Atm_titr_atom_type', 'Atm_titr_atom_isotope_number', 'Atm_titr_auth_entity_assembly_ID', 'Atm_titr_auth_seq_ID', 'Atm_titr_auth_comp_ID', 'Atm_titr_auth_atom_ID', 'Hill_coeff_val', 'Hill_coeff_val_fit_err', 'High_PH_param_fit_val', 'High_PH_param_fit_val_err', 'Low_PH_param_fit_val', 'Low_PH_param_fit_val_err', 'PKa_val', 'PKa_val_fit_err', 'PHmid_val', 'PHmid_val_fit_err', 'Entry_ID', 'PH_titration_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'PH_titration_list_ID'],

            }

        },

    'tableNames': ['PH_titration_experiment', 'PH_titration_software', 'PH_titr_result']

    },

  'pH_param_list': {

    'name': 'PH_param_list',

    'tags': {

      'Sf_category': ['pH_param_list',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'PH_titration_list_ID': [None,returnStarInt,'PH_titration_list.ID',True],
      'PH_titration_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Observed_NMR_param': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'PH_titration_list_ID', 'PH_titration_list_label', 'Observed_NMR_param', 'Data_file_name', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'PH_param': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'PH_titr_result_ID': [None,returnStarInt,'PH_titr_result.ID',True],
          'PH_val': [None,returnStarFloat,None,True],
          'PH_val_err': [None,returnStarFloat,None,False],
          'Observed_NMR_param_val': [None,returnStarFloat,None,True],
          'Observed_NMR_param_val_err': [None,returnStarFloat,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'PH_param_list_ID': [None,returnStarInt,'PH_param_list.ID',True],

                },

        'tagNames': ['ID', 'PH_titr_result_ID', 'PH_val', 'PH_val_err', 'Observed_NMR_param_val', 'Observed_NMR_param_val_err', 'Entry_ID', 'PH_param_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'PH_param_list_ID'],

            }

        },

    'tableNames': ['PH_param']

    },

  'D_H_fractionation_factors': {

    'name': 'D_H_fractionation_factor_list',

    'tags': {

      'Sf_category': ['D_H_fractionation_factors',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'D_H_fract_factor_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'D_H_fractionation_factor_list_ID': [None,returnStarInt,'D_H_fractionation_factor_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'D_H_fractionation_factor_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'D_H_fractionation_factor_list_ID'],

            },

      'D_H_fract_factor_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'D_H_fractionation_factor_list_ID': [None,returnStarInt,'D_H_fractionation_factor_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'D_H_fractionation_factor_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'D_H_fractionation_factor_list_ID'],

            },

      'D_H_fractionation_factor': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Val': [None,returnStarFloat,None,True],
          'Val_err': [None,returnStarFloat,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'D_H_fractionation_factor_list_ID': [None,returnStarInt,'D_H_fractionation_factor_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Val', 'Val_err', 'Resonance_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'D_H_fractionation_factor_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'D_H_fractionation_factor_list_ID'],

            }

        },

    'tableNames': ['D_H_fract_factor_experiment', 'D_H_fract_factor_software', 'D_H_fractionation_factor']

    },

  'deduced_secd_struct_features': {

    'name': 'Deduced_secd_struct_list',

    'tags': {

      'Sf_category': ['deduced_secd_struct_features',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Residue_struct_value_details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Details', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Residue_struct_value_details', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Deduced_secd_struct_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Deduced_secd_struct_list_ID': [None,returnStarInt,'Deduced_secd_struct_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Deduced_secd_struct_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Deduced_secd_struct_list_ID'],

            },

      'Deduced_secd_struct_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Deduced_secd_struct_list_ID': [None,returnStarInt,'Deduced_secd_struct_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Deduced_secd_struct_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Deduced_secd_struct_list_ID'],

            },

      'Deduced_secd_struct_exptl': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_start': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Comp_index_ID_end': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_start': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Seq_ID_end': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
          'Auth_seq_ID_start': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_end': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Code': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
          'Static_field_orientation_angle': [None,returnStarFloat,None,False],
          'Selection_method_ID': [None,returnStarInt,'Method.ID',False],
          'Selection_method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Deduced_secd_struct_list_ID': [None,returnStarInt,'Deduced_secd_struct_list.ID',True],

                },

        'tagNames': ['ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID_start', 'Comp_index_ID_end', 'Seq_ID_start', 'Seq_ID_end', 'Auth_seq_ID_start', 'Auth_seq_ID_end', 'Name', 'Code', 'Static_field_orientation_angle', 'Selection_method_ID', 'Selection_method_label', 'Details', 'Entry_ID', 'Deduced_secd_struct_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Deduced_secd_struct_list_ID'],

            },

      'Deduced_secd_struct_feature': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Chem_comp_struct_val': [None,returnStarFloat,None,True],
          'Chem_comp_struct_element_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Figure_of_merit': [None,returnStarFloat,None,False],
          'Spin_system_ID': [None,returnStarInt,None,False],
          'Auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Deduced_secd_struct_list_ID': [None,returnStarInt,'Deduced_secd_struct_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Chem_comp_struct_val', 'Chem_comp_struct_element_type', 'Figure_of_merit', 'Spin_system_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Deduced_secd_struct_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entity_assembly_ID', 'Comp_index_ID', 'Atom_ID', 'Entry_ID', 'Deduced_secd_struct_list_ID'],

            }

        },

    'tableNames': ['Deduced_secd_struct_software', 'Deduced_secd_struct_experiment', 'Deduced_secd_struct_exptl', 'Deduced_secd_struct_feature']

    },

  'deduced_hydrogen_bonds': {

    'name': 'Deduced_H_bond_list',

    'tags': {

      'Sf_category': ['deduced_hydrogen_bonds',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Details', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Deduced_H_bond_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Deduced_H_bond_list_ID': [None,returnStarInt,'Deduced_H_bond_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Deduced_H_bond_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Deduced_H_bond_list_ID'],

            },

      'Deduced_H_bond_experiment': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Deduced_H_bond_list_ID': [None,returnStarInt,'Deduced_H_bond_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Deduced_H_bond_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Deduced_H_bond_list_ID'],

            },

      'Deduced_H_bond': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Donor_atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Heavy_atom_entry_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Heavy_atom_entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Heavy_atom_entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Heavy_atom_comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Heavy_atom_seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Heavy_atom_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Heavy_atom_atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Heavy_atom_atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Heavy_atom_entry_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Heavy_atom_entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Heavy_atom_entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Heavy_atom_comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Heavy_atom_seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Heavy_atom_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Heavy_atom_atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Heavy_atom_atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Heavy_atom_auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Heavy_atom_auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Heavy_atom_auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Heavy_atom_auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Heavy_atom_auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Heavy_atom_auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Heavy_atom_auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Heavy_atom_auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'H_atom_entry_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'H_atom_entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'H_atom_entity_ID': [None,returnStarInt,'Entity.ID',True],
          'H_atom_comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'H_atom_seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'H_atom_comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'H_atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'H_atom_auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'H_atom_auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'H_atom_auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'H_atom_auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Figure_of_merit': [None,returnStarFloat,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Deduced_hydrogen_bond_list_ID': [None,returnStarInt,'Deduced_H_bond_list.ID',True],

                },

        'tagNames': ['ID', 'Donor_atom_ID', 'Heavy_atom_entry_atom_ID_1', 'Heavy_atom_entity_assembly_ID_1', 'Heavy_atom_entity_ID_1', 'Heavy_atom_comp_index_ID_1', 'Heavy_atom_seq_ID_1', 'Heavy_atom_comp_ID_1', 'Heavy_atom_atom_ID_1', 'Heavy_atom_atom_type_1', 'Heavy_atom_entry_atom_ID_2', 'Heavy_atom_entity_assembly_ID_2', 'Heavy_atom_entity_ID_2', 'Heavy_atom_comp_index_ID_2', 'Heavy_atom_seq_ID_2', 'Heavy_atom_comp_ID_2', 'Heavy_atom_atom_ID_2', 'Heavy_atom_atom_type_2', 'Heavy_atom_auth_entity_assembly_ID_1', 'Heavy_atom_auth_seq_ID_1', 'Heavy_atom_auth_comp_ID_1', 'Heavy_atom_auth_atom_ID_1', 'Heavy_atom_auth_entity_assembly_ID_2', 'Heavy_atom_auth_seq_ID_2', 'Heavy_atom_auth_comp_ID_2', 'Heavy_atom_auth_atom_ID_2', 'H_atom_entry_atom_ID', 'H_atom_entity_assembly_ID', 'H_atom_entity_ID', 'H_atom_comp_index_ID', 'H_atom_seq_ID', 'H_atom_comp_ID', 'H_atom_ID', 'H_atom_auth_entity_assembly_ID', 'H_atom_auth_seq_ID', 'H_atom_auth_comp_ID', 'H_atom_auth_atom_ID', 'Figure_of_merit', 'Entry_ID', 'Deduced_hydrogen_bond_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Deduced_hydrogen_bond_list_ID'],

            }

        },

    'tableNames': ['Deduced_H_bond_software', 'Deduced_H_bond_experiment', 'Deduced_H_bond']

    },

  'conformer_statistics': {

    'name': 'Conformer_stat_list',
    'saveFrameCode': 'conformer_statistics',

    'tags': {

      'Sf_category': ['conformer_statistics',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Conformer_ensemble_only': [None,lambda x = value: returnStarYesNo(x,length = 3),None,True],
      'Both_ensemble_and_rep_conformer': [None,lambda x = value: returnStarYesNo(x,length = 3),None,True],
      'Representative_conformer_only': [None,lambda x = value: returnStarYesNo(x,length = 3),None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],
      'Original_conformer_stats_file_ID': [None,returnStarInt,None,False],
      'Conf_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],
      'Conf_family_coord_set_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Representative_conformer_ID': [None,returnStarInt,'Representative_conformer.ID',False],
      'Representative_conformer_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
      'Conformer_calculated_total_num': [None,returnStarInt,None,True],
      'Conformer_submitted_total_num': [None,returnStarInt,None,True],
      'Conformer_selection_criteria': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Representative_conformer': [None,returnStarInt,None,True],
      'Rep_conformer_selection_criteria': [None,returnStarString,None,False],
      'Statistical_struct_param_details': [None,returnStarString,None,False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Conformer_ensemble_only', 'Both_ensemble_and_rep_conformer', 'Representative_conformer_only', 'Data_file_name', 'Text_data_format', 'Text_data', 'Original_conformer_stats_file_ID', 'Conf_family_coord_set_ID', 'Conf_family_coord_set_label', 'Representative_conformer_ID', 'Representative_conformer_label', 'Conformer_calculated_total_num', 'Conformer_submitted_total_num', 'Conformer_selection_criteria', 'Representative_conformer', 'Rep_conformer_selection_criteria', 'Statistical_struct_param_details', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Conformer_stat_list_ens': {

        'tags': {

          'Stats_not_available': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
          'Ramachan_most_favored_pct': [None,returnStarFloat,None,False],
          'Ramachan_allowed_pct': [None,returnStarFloat,None,False],
          'Ramachan_gen_allowed_pct': [None,returnStarFloat,None,False],
          'Ramachan_disallowed_pct': [None,returnStarFloat,None,False],
          'Total_E_value_': [None,returnStarFloat,None,False],
          'Total_E_value_err': [None,returnStarFloat,None,False],
          'Bond_E_value': [None,returnStarFloat,None,False],
          'Bond_E_value_err': [None,returnStarFloat,None,False],
          'Angle_E_value': [None,returnStarFloat,None,False],
          'Angle_E_value_err': [None,returnStarFloat,None,False],
          'Improper_E_value': [None,returnStarFloat,None,False],
          'Improper_E_value_err': [None,returnStarFloat,None,False],
          'Van_der_Waals_E_value': [None,returnStarFloat,None,False],
          'Van_der_Waals_E_value_err': [None,returnStarFloat,None,False],
          'Torsional_angle_E_value': [None,returnStarFloat,None,False],
          'Torsional_angle_E_value_err': [None,returnStarFloat,None,False],
          'NCS_E_value': [None,returnStarFloat,None,False],
          'NCS_E_value_err': [None,returnStarFloat,None,False],
          'Lennard_Jones_E_value': [None,returnStarFloat,None,False],
          'Lennard_Jones_E_value_err': [None,returnStarFloat,None,False],
          'Covalent_bond_rmsd': [None,returnStarFloat,None,False],
          'Covalent_bond_rmsd_err': [None,returnStarFloat,None,False],
          'Bond_angle_rmsd': [None,returnStarFloat,None,False],
          'Bond_angle_rmsd_err': [None,returnStarFloat,None,False],
          'Dihedral_angle_rmsd': [None,returnStarFloat,None,False],
          'Dihedral_angle_rmsd_err': [None,returnStarFloat,None,False],
          'Improper_torsion_angle_rmsd': [None,returnStarFloat,None,False],
          'Improper_torsion_angle_rmsd_err': [None,returnStarFloat,None,False],
          'Peptide_planarity_rmsd': [None,returnStarFloat,None,False],
          'Peptide_planarity_rmsd_err': [None,returnStarFloat,None,False],
          'Atm_coord_avg_rmsd_calc_method': [None,returnStarFloat,None,False],
          'BB_hvy_atm_coord_avg_rmsd': [None,returnStarFloat,None,False],
          'BB_hvy_atm_coord_avg_rmsd_err': [None,returnStarFloat,None,False],
          'BB_hvy_atm_coord_std_dev': [None,returnStarFloat,None,False],
          'BB_hvy_atm_coord_std_dev_err': [None,returnStarFloat,None,False],
          'BB_hvy_atm_residues_included': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'BB_hvy_atm_exclusions': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'SC_hvy_atm_coord_avg_rmsd': [None,returnStarFloat,None,False],
          'SC_hvy_atm_coord_avg_rmsd_err': [None,returnStarFloat,None,False],
          'SC_hvy_atm_coord_std_dev': [None,returnStarFloat,None,False],
          'SC_hvy_atm_coord_std_dev_err': [None,returnStarFloat,None,False],
          'SC_hvy_atm_residues_included': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'SC_hvy_atm_exclusions': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'All_hvy_atm_coord_avg_rmsd': [None,returnStarFloat,None,False],
          'All_hvy_atm_coord_avg_rmsd_err': [None,returnStarFloat,None,False],
          'All_hvy_atm_coord_std_dev': [None,returnStarFloat,None,False],
          'All_hvy_atm_coord_std_dev_err': [None,returnStarFloat,None,False],
          'All_hvy_atm_residues_included': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'All_hvy_atm_exclusions': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Conformer_stat_list_ID': [None,returnStarInt,'Conformer_stat_list.ID',True],

                },

        'tagNames': ['Stats_not_available', 'Ramachan_most_favored_pct', 'Ramachan_allowed_pct', 'Ramachan_gen_allowed_pct', 'Ramachan_disallowed_pct', 'Total_E_value_', 'Total_E_value_err', 'Bond_E_value', 'Bond_E_value_err', 'Angle_E_value', 'Angle_E_value_err', 'Improper_E_value', 'Improper_E_value_err', 'Van_der_Waals_E_value', 'Van_der_Waals_E_value_err', 'Torsional_angle_E_value', 'Torsional_angle_E_value_err', 'NCS_E_value', 'NCS_E_value_err', 'Lennard_Jones_E_value', 'Lennard_Jones_E_value_err', 'Covalent_bond_rmsd', 'Covalent_bond_rmsd_err', 'Bond_angle_rmsd', 'Bond_angle_rmsd_err', 'Dihedral_angle_rmsd', 'Dihedral_angle_rmsd_err', 'Improper_torsion_angle_rmsd', 'Improper_torsion_angle_rmsd_err', 'Peptide_planarity_rmsd', 'Peptide_planarity_rmsd_err', 'Atm_coord_avg_rmsd_calc_method', 'BB_hvy_atm_coord_avg_rmsd', 'BB_hvy_atm_coord_avg_rmsd_err', 'BB_hvy_atm_coord_std_dev', 'BB_hvy_atm_coord_std_dev_err', 'BB_hvy_atm_residues_included', 'BB_hvy_atm_exclusions', 'SC_hvy_atm_coord_avg_rmsd', 'SC_hvy_atm_coord_avg_rmsd_err', 'SC_hvy_atm_coord_std_dev', 'SC_hvy_atm_coord_std_dev_err', 'SC_hvy_atm_residues_included', 'SC_hvy_atm_exclusions', 'All_hvy_atm_coord_avg_rmsd', 'All_hvy_atm_coord_avg_rmsd_err', 'All_hvy_atm_coord_std_dev', 'All_hvy_atm_coord_std_dev_err', 'All_hvy_atm_residues_included', 'All_hvy_atm_exclusions', 'Entry_ID', 'Conformer_stat_list_ID'],
        'sourcePrimaryKeys': ['Entry_ID', 'Conformer_stat_list_ID'],

            },

      'Conformer_stat_list_rep': {

        'tags': {

          'Stats_not_available': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
          'Ramachan_most_favored_pct': [None,returnStarFloat,None,False],
          'Ramachan_allowed_pct': [None,returnStarFloat,None,False],
          'Ramachan_gen_allowed_pct': [None,returnStarFloat,None,False],
          'Ramachan_disallowed_pct': [None,returnStarFloat,None,False],
          'Total_E_value': [None,returnStarFloat,None,False],
          'Total_E_value_err': [None,returnStarFloat,None,False],
          'Bond_E_value': [None,returnStarFloat,None,False],
          'Bond_E_value_err': [None,returnStarFloat,None,False],
          'Angle_E_value': [None,returnStarFloat,None,False],
          'Angle_E_value_err': [None,returnStarFloat,None,False],
          'Impropers_E_value': [None,returnStarFloat,None,False],
          'Impropers_E_value_err': [None,returnStarFloat,None,False],
          'Van_der_Waals_E_val': [None,returnStarFloat,None,False],
          'Van_der_Waals_E_val_err': [None,returnStarFloat,None,False],
          'NOE_E_value': [None,returnStarFloat,None,False],
          'NOE_E_value_err': [None,returnStarFloat,None,False],
          'Torsional_E_value': [None,returnStarFloat,None,False],
          'Torsional_E_value_err': [None,returnStarFloat,None,False],
          'NCS_E_value': [None,returnStarFloat,None,False],
          'NCS_E_value_err': [None,returnStarFloat,None,False],
          'Lennard_Jones_E_value': [None,returnStarFloat,None,False],
          'Lennard_Jones_E_value_err': [None,returnStarFloat,None,False],
          'Bond_rmsd': [None,returnStarFloat,None,False],
          'Bond_rmsd_err': [None,returnStarFloat,None,False],
          'Angle_rmsd': [None,returnStarFloat,None,False],
          'Angle_rmsd_err': [None,returnStarFloat,None,False],
          'Improper_torsion_angle_rmsd': [None,returnStarFloat,None,False],
          'Improper_torsion_angle_rmsd_err': [None,returnStarFloat,None,False],
          'Peptide_planarity_rmsd': [None,returnStarFloat,None,False],
          'Peptide_planarity_rmsd_err': [None,returnStarFloat,None,False],
          'Struct_figure_of_merit': [None,returnStarFloat,None,False],
          'Struct_figure_of_merit_func_form': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Conformer_stat_list_ID': [None,returnStarInt,'Conformer_stat_list.ID',True],

                },

        'tagNames': ['Stats_not_available', 'Ramachan_most_favored_pct', 'Ramachan_allowed_pct', 'Ramachan_gen_allowed_pct', 'Ramachan_disallowed_pct', 'Total_E_value', 'Total_E_value_err', 'Bond_E_value', 'Bond_E_value_err', 'Angle_E_value', 'Angle_E_value_err', 'Impropers_E_value', 'Impropers_E_value_err', 'Van_der_Waals_E_val', 'Van_der_Waals_E_val_err', 'NOE_E_value', 'NOE_E_value_err', 'Torsional_E_value', 'Torsional_E_value_err', 'NCS_E_value', 'NCS_E_value_err', 'Lennard_Jones_E_value', 'Lennard_Jones_E_value_err', 'Bond_rmsd', 'Bond_rmsd_err', 'Angle_rmsd', 'Angle_rmsd_err', 'Improper_torsion_angle_rmsd', 'Improper_torsion_angle_rmsd_err', 'Peptide_planarity_rmsd', 'Peptide_planarity_rmsd_err', 'Struct_figure_of_merit', 'Struct_figure_of_merit_func_form', 'Entry_ID', 'Conformer_stat_list_ID'],
        'sourcePrimaryKeys': ['Entry_ID', 'Conformer_stat_list_ID'],

            },

      'Conf_stats_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Conformer_stat_list_ID': [None,returnStarInt,'Conformer_stat_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Conformer_stat_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Conformer_stat_list_ID'],

            }

        },

    'tableNames': ['Conformer_stat_list_ens', 'Conformer_stat_list_rep', 'Conf_stats_software']

    },

  'conformer_family_coord_set': {

    'name': 'Conformer_family_coord_set',
    'saveFrameCode': 'ensemble_of_conformers',

    'tags': {

      'Sf_category': ['conformer_family_coord_set',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'File_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Constraints_PDB_file_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
      'PDB_accession_code': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
      'Sample_condition_list_ID': [None,returnStarInt,'Sample_condition_list.ID',True],
      'Sample_condition_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
      'Atom_site_uncertainty_desc': [None,returnStarString,None,False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'File_name', 'Constraints_PDB_file_ID', 'PDB_accession_code', 'Sample_condition_list_ID', 'Sample_condition_list_label', 'Atom_site_uncertainty_desc', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Conformer_family_refinement': {

        'tags': {

          'Refine_method': [None,returnStarString,None,True],
          'Refine_details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

                },

        'tagNames': ['Refine_method', 'Refine_details', 'Entry_ID', 'Conformer_family_coord_set_ID'],
        'sourcePrimaryKeys': ['Refine_method', 'Entry_ID', 'Conformer_family_coord_set_ID'],

            },

      'Conformer_family_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Conformer_family_coord_set_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Conformer_family_coord_set_ID'],

            },

      'Energetic_penalty_function': {

        'tags': {

          'Function': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Description': [None,lambda x = value: returnStarString(x,length = 255),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

                },

        'tagNames': ['Function', 'Description', 'Entry_ID', 'Conformer_family_coord_set_ID'],
        'sourcePrimaryKeys': ['Function', 'Entry_ID', 'Conformer_family_coord_set_ID'],

            },

      'Conformer_family_coord_set_expt': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Conformer_family_coord_set_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Conformer_family_coord_set_ID'],

            },

      'Conf_family_coord_set_constr_list': {

        'tags': {

          'Constraint_list_category': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Constraint_list_ID': [None,returnStarInt,None,True],
          'Constraint_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

                },

        'tagNames': ['Constraint_list_category', 'Constraint_list_ID', 'Constraint_list_label', 'Entry_ID', 'Conformer_family_coord_set_ID'],
        'sourcePrimaryKeys': ['Constraint_list_category', 'Constraint_list_ID', 'Entry_ID', 'Conformer_family_coord_set_ID'],

            },

      'Struct_image': {

        'tags': {

          'File_name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'File_format': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

                },

        'tagNames': ['File_name', 'File_format', 'Details', 'Entry_ID', 'Conformer_family_coord_set_ID'],
        'sourcePrimaryKeys': ['File_name', 'Entry_ID', 'Conformer_family_coord_set_ID'],

            },

      'Local_structure_quality': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 255),None,True],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Asym_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Entity_ID': [None,lambda x = value: returnStarString(x,length = 12),'Entity.ID',False],
          'Entity_label': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Comp_index_ID_begin': [None,returnStarInt,'unknown.ID',True],
          'Comp_index_ID_end': [None,returnStarInt,'unknown.ID',True],
          'Entry_ID': [None,lambda x = value: returnStarString(x,length = 12),'Entry.ID',True],
          'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

                },

        'tagNames': ['ID', 'Type', 'Entity_assembly_ID', 'Asym_ID', 'Entity_ID', 'Entity_label', 'Comp_index_ID_begin', 'Comp_index_ID_end', 'Entry_ID', 'Conformer_family_coord_set_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Conformer_family_coord_set_ID'],

            },

      'Model_type': {

        'tags': {

          'Atom_site_model_ID': [None,returnStarInt,'Atom_site.Model_ID',True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 255),None,True],
          'Entry_ID': [None,lambda x = value: returnStarString(x,length = 12),'Entry.ID',True],
          'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

                },

        'tagNames': ['Atom_site_model_ID', 'Type', 'Entry_ID', 'Conformer_family_coord_set_ID'],
        'sourcePrimaryKeys': ['Atom_site_model_ID', 'Type'],

            },

      'Atom_site': {

        'tags': {

          'Model_ID': [None,returnStarInt,None,True],
          'Model_site_ID': [None,returnStarInt,None,False],
          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Label_entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Label_asym_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Struct_asym.ID',False],
          'Label_entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Label_comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Label_comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Label_seq_ID': [None,returnStarInt,'Entity_poly_seq.Num',True],
          'Label_atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Type_symbol': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Cartn_x': [None,returnStarFloat,None,True],
          'Cartn_y': [None,returnStarFloat,None,True],
          'Cartn_z': [None,returnStarFloat,None,True],
          'Cartn_x_esd': [None,returnStarFloat,None,False],
          'Cartn_y_esd': [None,returnStarFloat,None,False],
          'Cartn_z_esd': [None,returnStarFloat,None,False],
          'Occupancy': [None,returnStarFloat,None,False],
          'Occupancy_esd': [None,returnStarFloat,None,False],
          'Uncertainty': [None,returnStarFloat,None,False],
          'Ordered_flag': [None,returnStarInt,None,False],
          'Footnote_ID': [None,returnStarInt,None,False],
          'PDB_seq_ID': [None,returnStarInt,None,False],
          'PDB_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),'Chem_comp.ID',True],
          'PDB_atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 15),None,False],
          'Auth_asym_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_name': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

                },

        'tagNames': ['Model_ID', 'Model_site_ID', 'ID', 'Assembly_atom_ID', 'Label_entity_assembly_ID', 'Label_asym_ID', 'Label_entity_ID', 'Label_comp_index_ID', 'Label_comp_ID', 'Label_seq_ID', 'Label_atom_ID', 'Type_symbol', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'Cartn_x_esd', 'Cartn_y_esd', 'Cartn_z_esd', 'Occupancy', 'Occupancy_esd', 'Uncertainty', 'Ordered_flag', 'Footnote_ID', 'PDB_seq_ID', 'PDB_comp_ID', 'PDB_atom_ID', 'Auth_asym_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Auth_atom_name', 'Details', 'Entry_ID', 'Conformer_family_coord_set_ID'],
        'sourcePrimaryKeys': ['Model_ID', 'ID', 'Label_asym_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Entry_ID', 'Conformer_family_coord_set_ID'],

            },

      'Atom_sites_footnote': {

        'tags': {

          'Footnote_ID': [None,returnStarInt,None,True],
          'Text': [None,returnStarString,None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

                },

        'tagNames': ['Footnote_ID', 'Text', 'Entry_ID', 'Conformer_family_coord_set_ID'],
        'sourcePrimaryKeys': ['Footnote_ID', 'Entry_ID', 'Conformer_family_coord_set_ID'],

            }

        },

    'tableNames': ['Conformer_family_refinement', 'Conformer_family_software', 'Energetic_penalty_function', 'Conformer_family_coord_set_expt', 'Conf_family_coord_set_constr_list', 'Struct_image', 'Local_structure_quality', 'Model_type', 'Atom_site', 'Atom_sites_footnote']

    },

  'representative_conformer': {

    'name': 'Representative_conformer',
    'saveFrameCode': 'representative_conformer',

    'tags': {

      'Sf_category': ['representative_conformer',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],
      'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Rep_conformer_derivation': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Atom_pos_uncertainty_derivation': [None,returnStarString,None,False],
      'Rep_conformer_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Rep_conformer_original_file': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'PDB_accession_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
      'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],
      'Conformer_family_coord_set_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Details', 'Type', 'Rep_conformer_derivation', 'Atom_pos_uncertainty_derivation', 'Rep_conformer_file_name', 'Rep_conformer_original_file', 'PDB_accession_code', 'Conformer_family_coord_set_ID', 'Conformer_family_coord_set_label'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Rep_conf_refinement': {

        'tags': {

          'Refine_method': [None,returnStarString,None,True],
          'Refine_details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Representative_conformer_ID': [None,returnStarInt,'Representative_conformer.ID',True],

                },

        'tagNames': ['Refine_method', 'Refine_details', 'Entry_ID', 'Representative_conformer_ID'],
        'sourcePrimaryKeys': ['Refine_method', 'Entry_ID', 'Representative_conformer_ID'],

            },

      'Rep_conf_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Representative_conformer_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Representative_conformer_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Representative_conformer_ID'],

            },

      'Terminal_residue': {

        'tags': {

          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Representative_conformer_ID': [None,returnStarInt,'Representative_conformer.ID',True],

                },

        'tagNames': ['Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Entry_ID', 'Representative_conformer_ID'],
        'sourcePrimaryKeys': ['Entity_ID', 'Comp_index_ID', 'Entry_ID', 'Representative_conformer_ID'],

            },

      'Rep_conf': {

        'tags': {

          'Atom_coordinate_ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Atom_site_ID': [None,returnStarInt,'Atom_site.ID',True],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Chem_comp_PDB_ID_code': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_cartn_x': [None,returnStarFloat,None,True],
          'Atom_cartn_y': [None,returnStarFloat,None,True],
          'Atom_cartn_z': [None,returnStarFloat,None,True],
          'Atom_cartn_x_esd': [None,returnStarFloat,None,False],
          'Atom_cartn_y_esd': [None,returnStarFloat,None,False],
          'Atom_cartn_z_esd': [None,returnStarFloat,None,False],
          'Atom_position_uncertainty': [None,returnStarFloat,None,False],
          'Atom_coord_footnote_ID': [None,returnStarInt,'Rep_coordinate_details.Footnote_ID',False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Representative_conformer_ID': [None,returnStarInt,'Representative_conformer.ID',True],

                },

        'tagNames': ['Atom_coordinate_ID', 'Assembly_atom_ID', 'Atom_site_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Comp_ID', 'Chem_comp_PDB_ID_code', 'Seq_ID', 'Atom_ID', 'Atom_type', 'Atom_cartn_x', 'Atom_cartn_y', 'Atom_cartn_z', 'Atom_cartn_x_esd', 'Atom_cartn_y_esd', 'Atom_cartn_z_esd', 'Atom_position_uncertainty', 'Atom_coord_footnote_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Representative_conformer_ID'],
        'sourcePrimaryKeys': ['Atom_coordinate_ID', 'Entry_ID', 'Representative_conformer_ID'],

            },

      'Rep_coordinate_details': {

        'tags': {

          'Footnote_ID': [None,returnStarInt,None,True],
          'Footnote': [None,lambda x = value: returnStarString(x,length = 255),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Representative_conformer_ID': [None,returnStarInt,'Representative_conformer.ID',True],

                },

        'tagNames': ['Footnote_ID', 'Footnote', 'Entry_ID', 'Representative_conformer_ID'],
        'sourcePrimaryKeys': ['Footnote_ID', 'Entry_ID', 'Representative_conformer_ID'],

            }

        },

    'tableNames': ['Rep_conf_refinement', 'Rep_conf_software', 'Terminal_residue', 'Rep_conf', 'Rep_coordinate_details']

    },

  'constraint_statistics': {

    'name': 'Constraint_stat_list',
    'saveFrameCode': 'constraint_statistics',

    'tags': {

      'Sf_category': ['constraint_statistics',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],
      'Stats_not_available': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
      'NOE_interproton_dist_evaluation': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'NOE_pseudoatom_corrections': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'NOE_dist_averaging_method': [None,lambda x = value: returnStarString(x,length = 127),None,False],
      'ROE_interproton_dist_evaluation': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'ROE_pseudoatom_corrections': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'ROE_dist_averaging_method': [None,returnStarFloat,None,False],
      'NOE_tot_num': [None,returnStarInt,None,False],
      'RDC_tot_num': [None,returnStarInt,None,False],
      'Dihedral_angle_tot_num': [None,returnStarInt,None,False],
      'Protein_dihedral_angle_tot_num': [None,returnStarInt,None,False],
      'NA_dihedral_angle_tot_num': [None,returnStarInt,None,False],
      'NOE_intraresidue_tot_num': [None,returnStarInt,None,False],
      'NOE_sequential_tot_num': [None,returnStarInt,None,False],
      'NOE_medium_range_tot_num': [None,returnStarInt,None,False],
      'NOE_long_range_tot_num': [None,returnStarInt,None,False],
      'NOE_unique_tot_num': [None,returnStarInt,None,False],
      'NOE_intraresidue_unique_tot_num': [None,returnStarInt,None,False],
      'NOE_sequential_unique_tot_num': [None,returnStarInt,None,False],
      'NOE_medium_range_unique_tot_num': [None,returnStarInt,None,False],
      'NOE_long_range_unique_tot_num': [None,returnStarInt,None,False],
      'NOE_unamb_intramol_tot_num': [None,returnStarInt,None,False],
      'NOE_unamb_intermol_tot_num': [None,returnStarInt,None,False],
      'NOE_ambig_intramol_tot_num': [None,returnStarInt,None,False],
      'NOE_ambig_intermol_tot_num': [None,returnStarInt,None,False],
      'NOE_interentity_tot_num': [None,returnStarInt,None,False],
      'NOE_other_tot_num': [None,returnStarInt,None,False],
      'ROE_tot_num': [None,returnStarInt,None,False],
      'ROE_intraresidue_tot_num': [None,returnStarInt,None,False],
      'ROE_sequential_tot_num': [None,returnStarInt,None,False],
      'ROE_medium_range_tot_num': [None,returnStarInt,None,False],
      'ROE_long_range_tot_num': [None,returnStarInt,None,False],
      'ROE_unambig_intramol_tot_num': [None,returnStarInt,None,False],
      'ROE_unambig_intermol_tot_num': [None,returnStarInt,None,False],
      'ROE_ambig_intramol_tot_num': [None,returnStarInt,None,False],
      'ROE_ambig_intermol_tot_num': [None,returnStarInt,None,False],
      'ROE_other_tot_num': [None,returnStarInt,None,False],
      'RDC_HH_tot_num': [None,returnStarInt,None,False],
      'RDC_HNC_tot_num': [None,returnStarInt,None,False],
      'RDC_NH_tot_num': [None,returnStarInt,None,False],
      'RDC_CC_tot_num': [None,returnStarInt,None,False],
      'RDC_CN_i_1_tot_num': [None,returnStarInt,None,False],
      'RDC_CAHA_tot_num': [None,returnStarInt,None,False],
      'RDC_HNHA_tot_num': [None,returnStarInt,None,False],
      'RDC_HNHA_i_1_tot_num': [None,returnStarInt,None,False],
      'RDC_CAC_tot_num': [None,returnStarInt,None,False],
      'RDC_CAN_tot_num': [None,returnStarInt,None,False],
      'RDC_intraresidue_tot_num': [None,returnStarInt,None,False],
      'RDC_sequential_tot_num': [None,returnStarInt,None,False],
      'RDC_medium_range_tot_num': [None,returnStarInt,None,False],
      'RDC_long_range_tot_num': [None,returnStarInt,None,False],
      'RDC_other_tot_num': [None,returnStarInt,None,False],
      'RDC_unambig_intramol_tot_num': [None,returnStarInt,None,False],
      'RDC_unambig_intermol_tot_num': [None,returnStarInt,None,False],
      'RDC_ambig_intramol_tot_num': [None,returnStarInt,None,False],
      'RDC_ambig_intermol_tot_num': [None,returnStarInt,None,False],
      'RDC_intermol_tot_num': [None,returnStarInt,None,False],
      'Protein_phi_angle_tot_num': [None,returnStarInt,None,False],
      'Protein_psi_angle_tot_num': [None,returnStarInt,None,False],
      'Protein_chi_one_angle_tot_num': [None,returnStarInt,None,False],
      'Protein_other_angle_tot_num': [None,returnStarInt,None,False],
      'Protein_ambig_dihedral_tot_num': [None,returnStarInt,None,False],
      'Protein_other_tot_num': [None,returnStarInt,None,False],
      'NA_alpha_angle_tot_num': [None,returnStarInt,None,False],
      'NA_beta_angle_tot_num': [None,returnStarInt,None,False],
      'NA_gamma_angle_tot_num': [None,returnStarInt,None,False],
      'NA_delta_angle_tot_num': [None,returnStarInt,None,False],
      'NA_epsilon_angle_tot_num': [None,returnStarInt,None,False],
      'NA_chi_angle_tot_num': [None,returnStarInt,None,False],
      'NA_sugar_pucker_tot_num': [None,returnStarInt,None,False],
      'NA_other_angle_tot_num': [None,returnStarInt,None,False],
      'NA_amb_dihedral_angle_tot_num': [None,returnStarInt,None,False],
      'NA_other_tot_num': [None,returnStarInt,None,False],
      'H_bonds_constrained_tot_num': [None,returnStarInt,None,False],
      'Constr_def_H_bonds_tot_num': [None,returnStarInt,None,False],
      'SS_bonds_constrained_tot_num': [None,returnStarInt,None,False],
      'Constr_def_SS_bonds_tot_num': [None,returnStarInt,None,False],
      'Derived_coupling_const_tot_num': [None,returnStarInt,None,False],
      'Derived_CACB_chem_shift_tot_num': [None,returnStarInt,None,False],
      'Derived_1H_chem_shifts_tot_num': [None,returnStarInt,None,False],
      'Derived_photo_cidnps_tot_num': [None,returnStarInt,None,False],
      'Derived_paramag_relax_tot_num': [None,returnStarInt,None,False],
      'Assumed_distances_tot_num': [None,returnStarInt,None,False],
      'Assumed_angles_tot_num': [None,returnStarInt,None,False],
      'Constraints_per_residue_avg': [None,returnStarFloat,None,False],
      'Constr_violations_per_residue_avg': [None,returnStarFloat,None,False],
      'Dist_constr_violat_stat_calc_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Details', 'Text_data_format', 'Text_data', 'Stats_not_available', 'NOE_interproton_dist_evaluation', 'NOE_pseudoatom_corrections', 'NOE_dist_averaging_method', 'ROE_interproton_dist_evaluation', 'ROE_pseudoatom_corrections', 'ROE_dist_averaging_method', 'NOE_tot_num', 'RDC_tot_num', 'Dihedral_angle_tot_num', 'Protein_dihedral_angle_tot_num', 'NA_dihedral_angle_tot_num', 'NOE_intraresidue_tot_num', 'NOE_sequential_tot_num', 'NOE_medium_range_tot_num', 'NOE_long_range_tot_num', 'NOE_unique_tot_num', 'NOE_intraresidue_unique_tot_num', 'NOE_sequential_unique_tot_num', 'NOE_medium_range_unique_tot_num', 'NOE_long_range_unique_tot_num', 'NOE_unamb_intramol_tot_num', 'NOE_unamb_intermol_tot_num', 'NOE_ambig_intramol_tot_num', 'NOE_ambig_intermol_tot_num', 'NOE_interentity_tot_num', 'NOE_other_tot_num', 'ROE_tot_num', 'ROE_intraresidue_tot_num', 'ROE_sequential_tot_num', 'ROE_medium_range_tot_num', 'ROE_long_range_tot_num', 'ROE_unambig_intramol_tot_num', 'ROE_unambig_intermol_tot_num', 'ROE_ambig_intramol_tot_num', 'ROE_ambig_intermol_tot_num', 'ROE_other_tot_num', 'RDC_HH_tot_num', 'RDC_HNC_tot_num', 'RDC_NH_tot_num', 'RDC_CC_tot_num', 'RDC_CN_i_1_tot_num', 'RDC_CAHA_tot_num', 'RDC_HNHA_tot_num', 'RDC_HNHA_i_1_tot_num', 'RDC_CAC_tot_num', 'RDC_CAN_tot_num', 'RDC_intraresidue_tot_num', 'RDC_sequential_tot_num', 'RDC_medium_range_tot_num', 'RDC_long_range_tot_num', 'RDC_other_tot_num', 'RDC_unambig_intramol_tot_num', 'RDC_unambig_intermol_tot_num', 'RDC_ambig_intramol_tot_num', 'RDC_ambig_intermol_tot_num', 'RDC_intermol_tot_num', 'Protein_phi_angle_tot_num', 'Protein_psi_angle_tot_num', 'Protein_chi_one_angle_tot_num', 'Protein_other_angle_tot_num', 'Protein_ambig_dihedral_tot_num', 'Protein_other_tot_num', 'NA_alpha_angle_tot_num', 'NA_beta_angle_tot_num', 'NA_gamma_angle_tot_num', 'NA_delta_angle_tot_num', 'NA_epsilon_angle_tot_num', 'NA_chi_angle_tot_num', 'NA_sugar_pucker_tot_num', 'NA_other_angle_tot_num', 'NA_amb_dihedral_angle_tot_num', 'NA_other_tot_num', 'H_bonds_constrained_tot_num', 'Constr_def_H_bonds_tot_num', 'SS_bonds_constrained_tot_num', 'Constr_def_SS_bonds_tot_num', 'Derived_coupling_const_tot_num', 'Derived_CACB_chem_shift_tot_num', 'Derived_1H_chem_shifts_tot_num', 'Derived_photo_cidnps_tot_num', 'Derived_paramag_relax_tot_num', 'Assumed_distances_tot_num', 'Assumed_angles_tot_num', 'Constraints_per_residue_avg', 'Constr_violations_per_residue_avg', 'Dist_constr_violat_stat_calc_method'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Constraint_stat_list_ens': {

        'tags': {

          'Constraint_stats_not_available': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
          'Dist_Constraint_violation_max': [None,returnStarFloat,None,False],
          'Upper_dist_constr_violat_max': [None,returnStarFloat,None,False],
          'Lower_dist_constr_violat_max': [None,returnStarFloat,None,False],
          'Dist_Constraint_violation_avg': [None,returnStarFloat,None,False],
          'All_dist_rmsd': [None,returnStarFloat,None,False],
          'All_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Intraresidue_dist_rmsd': [None,returnStarFloat,None,False],
          'Intraresidue_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Sequential_dist_rmsd': [None,returnStarFloat,None,False],
          'Sequential_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Short_range_dist_rmsd': [None,returnStarFloat,None,False],
          'Short_range_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Long_range_dist_rmsd': [None,returnStarFloat,None,False],
          'Long_range_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Unamb_intermol_dist_rmsd': [None,returnStarFloat,None,False],
          'Unamb_intermol_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Amb_intermol_dist_rmsd': [None,returnStarFloat,None,False],
          'Amb_intermol_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Amb_intramol_dist_rmsd': [None,returnStarFloat,None,False],
          'Amb_intramol_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Hydrogen_bond_rmsd': [None,returnStarFloat,None,False],
          'Hydrogen_bond_rmsd_err': [None,returnStarFloat,None,False],
          'Dihedral_const_stat_calc_meth': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Dihedral_const_violat_max': [None,returnStarFloat,None,False],
          'Dihedral_const_violat_avg': [None,returnStarFloat,None,False],
          'Dihedral_angle_rmsd': [None,returnStarFloat,None,False],
          'Dihedral_angle_rmsd_err': [None,returnStarFloat,None,False],
          'Dipolar_1H_1H_rmsd': [None,returnStarFloat,None,False],
          'Dipolar_1H_1H_rmsd_err': [None,returnStarFloat,None,False],
          'Dipolar_1H_15N_rmsd': [None,returnStarFloat,None,False],
          'Dipolar_1H_15N_rmsd_err': [None,returnStarFloat,None,False],
          'Dipolar_1H_13C_rmsd': [None,returnStarFloat,None,False],
          'Dipolar_1H_13C_rmsd_err': [None,returnStarFloat,None,False],
          'Dipolar_13C_13C_rmsd': [None,returnStarFloat,None,False],
          'Dipolar_13C_13C_rmsd_err': [None,returnStarFloat,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Constraint_stat_list_ID': [None,returnStarInt,'Constraint_stat_list.ID',True],

                },

        'tagNames': ['Constraint_stats_not_available', 'Dist_Constraint_violation_max', 'Upper_dist_constr_violat_max', 'Lower_dist_constr_violat_max', 'Dist_Constraint_violation_avg', 'All_dist_rmsd', 'All_dist_rmsd_err', 'Intraresidue_dist_rmsd', 'Intraresidue_dist_rmsd_err', 'Sequential_dist_rmsd', 'Sequential_dist_rmsd_err', 'Short_range_dist_rmsd', 'Short_range_dist_rmsd_err', 'Long_range_dist_rmsd', 'Long_range_dist_rmsd_err', 'Unamb_intermol_dist_rmsd', 'Unamb_intermol_dist_rmsd_err', 'Amb_intermol_dist_rmsd', 'Amb_intermol_dist_rmsd_err', 'Amb_intramol_dist_rmsd', 'Amb_intramol_dist_rmsd_err', 'Hydrogen_bond_rmsd', 'Hydrogen_bond_rmsd_err', 'Dihedral_const_stat_calc_meth', 'Dihedral_const_violat_max', 'Dihedral_const_violat_avg', 'Dihedral_angle_rmsd', 'Dihedral_angle_rmsd_err', 'Dipolar_1H_1H_rmsd', 'Dipolar_1H_1H_rmsd_err', 'Dipolar_1H_15N_rmsd', 'Dipolar_1H_15N_rmsd_err', 'Dipolar_1H_13C_rmsd', 'Dipolar_1H_13C_rmsd_err', 'Dipolar_13C_13C_rmsd', 'Dipolar_13C_13C_rmsd_err', 'Entry_ID', 'Constraint_stat_list_ID'],
        'sourcePrimaryKeys': ['Entry_ID', 'Constraint_stat_list_ID'],

            },

      'Constraint_stat_list_rep': {

        'tags': {

          'Constraint_stats_not_available': [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
          'Dist_constraint_viol_max': [None,returnStarFloat,None,False],
          'Upper_dist_constr_violat_max': [None,returnStarFloat,None,False],
          'Lower_dist_constr_violat_max': [None,returnStarFloat,None,False],
          'Dist_Constraint_violation_avg': [None,returnStarFloat,None,False],
          'Intraresidue_dist_rmsd': [None,returnStarFloat,None,False],
          'Intraresidue_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Sequential_dist_rmsd': [None,returnStarFloat,None,False],
          'Sequential_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Short_range_dist_rmsd': [None,returnStarFloat,None,False],
          'Short_range_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Long_range_dist_rmsd': [None,returnStarFloat,None,False],
          'Long_range_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Unamb_intermol_dist_rmsd': [None,returnStarFloat,None,False],
          'Unamb_intermol_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Amb_intermol_dist_rmsd': [None,returnStarFloat,None,False],
          'Amb_intermol_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Amb_intramol_dist_rmsd': [None,returnStarFloat,None,False],
          'Amb_intramol_dist_rmsd_err': [None,returnStarFloat,None,False],
          'Hydrogen_bond_rmsd': [None,returnStarFloat,None,False],
          'Hydrogen_bond_rmsd_err': [None,returnStarFloat,None,False],
          'Dihedral_angle_rmsd': [None,returnStarFloat,None,False],
          'Dihedral_angle_rmsd_err': [None,returnStarFloat,None,False],
          'Dipolar_1H_1H_rmsd': [None,returnStarFloat,None,False],
          'Dipolar_1H_1H_rmsd_err': [None,returnStarFloat,None,False],
          'Dipolar_1H_13C_rmsd': [None,returnStarFloat,None,False],
          'Dipolar_1H_13C_rmsd_err': [None,returnStarFloat,None,False],
          'Dipolar_1H_15N_rmsd': [None,returnStarFloat,None,False],
          'Dipolar_1H_15N_rmsd_err': [None,returnStarFloat,None,False],
          'Dipolar_13C_13C_rmsd': [None,returnStarFloat,None,False],
          'Dipolar_13C_13C_rmsd_err': [None,returnStarFloat,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Constraint_stat_list_ID': [None,returnStarInt,'Constraint_stat_list.ID',True],

                },

        'tagNames': ['Constraint_stats_not_available', 'Dist_constraint_viol_max', 'Upper_dist_constr_violat_max', 'Lower_dist_constr_violat_max', 'Dist_Constraint_violation_avg', 'Intraresidue_dist_rmsd', 'Intraresidue_dist_rmsd_err', 'Sequential_dist_rmsd', 'Sequential_dist_rmsd_err', 'Short_range_dist_rmsd', 'Short_range_dist_rmsd_err', 'Long_range_dist_rmsd', 'Long_range_dist_rmsd_err', 'Unamb_intermol_dist_rmsd', 'Unamb_intermol_dist_rmsd_err', 'Amb_intermol_dist_rmsd', 'Amb_intermol_dist_rmsd_err', 'Amb_intramol_dist_rmsd', 'Amb_intramol_dist_rmsd_err', 'Hydrogen_bond_rmsd', 'Hydrogen_bond_rmsd_err', 'Dihedral_angle_rmsd', 'Dihedral_angle_rmsd_err', 'Dipolar_1H_1H_rmsd', 'Dipolar_1H_1H_rmsd_err', 'Dipolar_1H_13C_rmsd', 'Dipolar_1H_13C_rmsd_err', 'Dipolar_1H_15N_rmsd', 'Dipolar_1H_15N_rmsd_err', 'Dipolar_13C_13C_rmsd', 'Dipolar_13C_13C_rmsd_err', 'Entry_ID', 'Constraint_stat_list_ID'],
        'sourcePrimaryKeys': ['Entry_ID', 'Constraint_stat_list_ID'],

            },

      'Constraint_stats_constr_list': {

        'tags': {

          'Constraint_list_category': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Constraint_list_ID': [None,returnStarInt,None,True],
          'Constraint_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Constraint_stat_list_ID': [None,returnStarInt,'Constraint_stat_list.ID',True],

                },

        'tagNames': ['Constraint_list_category', 'Constraint_list_ID', 'Constraint_list_label', 'Entry_ID', 'Constraint_stat_list_ID'],
        'sourcePrimaryKeys': ['Constraint_list_category', 'Constraint_list_ID', 'Entry_ID', 'Constraint_stat_list_ID'],

            },

      'Constraint_file': {

        'tags': {

          'ID': [None,returnStarInt,None,False],
          'Constraint_filename': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Software_ID': [None,returnStarInt,'Software.ID',False],
          'Software_label': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Software_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
          'Block_ID': [None,returnStarInt,None,False],
          'Constraint_type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Constraint_subtype': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Constraint_subsubtype': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
          'Constraint_number': [None,returnStarInt,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Constraint_stat_list_ID': [None,returnStarInt,'Constraint_stat_list.ID',True],

                },

        'tagNames': ['ID', 'Constraint_filename', 'Software_ID', 'Software_label', 'Software_name', 'Block_ID', 'Constraint_type', 'Constraint_subtype', 'Constraint_subsubtype', 'Constraint_number', 'Entry_ID', 'Constraint_stat_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Constraint_filename', 'Block_ID', 'Constraint_type', 'Constraint_subtype', 'Entry_ID', 'Constraint_stat_list_ID'],

            }

        },

    'tableNames': ['Constraint_stat_list_ens', 'Constraint_stat_list_rep', 'Constraint_stats_constr_list', 'Constraint_file']

    },

  'force_constants': {

    'name': 'Force_constant_list',

    'tags': {

      'Sf_category': ['force_constants',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Default_software_values_used': [None,lambda x = value: returnStarYesNo(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Default_software_values_used', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Force_constant_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Force_constant_list_ID': [None,returnStarInt,'Force_constant_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Force_constant_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Force_constant_list_ID'],

            },

      'Force_constant': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Term': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Units': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
          'Val': [None,returnStarFloat,None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Force_constant_list_ID': [None,returnStarInt,'Force_constant_list.ID',True],

                },

        'tagNames': ['ID', 'Term', 'Units', 'Val', 'Entry_ID', 'Force_constant_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Force_constant_list_ID'],

            }

        },

    'tableNames': ['Force_constant_software', 'Force_constant']

    },

  'angular_order_parameters': {

    'name': 'Angular_order_parameter_list',

    'tags': {

      'Sf_category': ['angular_order_parameters',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Text_data_format': [None,lambda x = value: returnStarLine(x,length = 31),None,False],
      'Text_data': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Text_data_format', 'Text_data'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Angular_order_param': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Phi_S_angle_val': [None,returnStarFloat,None,False],
          'Psi_S_angle_val': [None,returnStarFloat,None,False],
          'Chi_1_S_angle_val': [None,returnStarFloat,None,False],
          'Chi_2_S_angle_val': [None,returnStarFloat,None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Angular_order_parameter_list_ID': [None,returnStarInt,'Angular_order_parameter_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Phi_S_angle_val', 'Psi_S_angle_val', 'Chi_1_S_angle_val', 'Chi_2_S_angle_val', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Angular_order_parameter_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Angular_order_parameter_list_ID'],

            }

        },

    'tableNames': ['Angular_order_param']

    },

  'tertiary_struct_elements': {

    'name': 'Tertiary_struct_element_list',

    'tags': {

      'Sf_category': ['tertiary_struct_elements',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Tertiary_struct_element_sel': {

        'tags': {

          'Method_ID': [None,returnStarInt,'Method.ID',True],
          'Method': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Tertiary_struct_element_list_ID': [None,returnStarInt,'Tertiary_struct_element_list.ID',True],

                },

        'tagNames': ['Method_ID', 'Method', 'Entry_ID', 'Tertiary_struct_element_list_ID'],
        'sourcePrimaryKeys': ['Method_ID', 'Entry_ID', 'Tertiary_struct_element_list_ID'],

            },

      'Tertiary_struct': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Tertiary_struct_code': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Tertiary_struct_name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Entity_label': [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',False],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Tertiary_struct_element_code': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Details': [None,returnStarString,None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Tertiary_struct_element_list_ID': [None,returnStarInt,'Tertiary_struct_element_list.ID',True],

                },

        'tagNames': ['ID', 'Tertiary_struct_code', 'Tertiary_struct_name', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Entity_label', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Tertiary_struct_element_code', 'Details', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Tertiary_struct_element_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entity_assembly_ID', 'Comp_index_ID', 'Atom_ID', 'Entry_ID', 'Tertiary_struct_element_list_ID'],

            }

        },

    'tableNames': ['Tertiary_struct_element_sel', 'Tertiary_struct']

    },

  'secondary_structs': {

    'name': 'Secondary_struct_list',

    'tags': {

      'Sf_category': ['secondary_structs',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Representative_conformer_ID': [None,returnStarInt,'Representative_conformer.ID',False],
      'Representative_conformer_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Representative_conformer_ID', 'Representative_conformer_label'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Secondary_struct_sel': {

        'tags': {

          'Method_ID': [None,returnStarInt,'Method.ID',True],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Secondary_struct_list_ID': [None,returnStarInt,'Secondary_struct_list.ID',True],

                },

        'tagNames': ['Method_ID', 'Method_label', 'Entry_ID', 'Secondary_struct_list_ID'],
        'sourcePrimaryKeys': ['Method_ID', 'Entry_ID', 'Secondary_struct_list_ID'],

            },

      'Secondary_struct': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Conf_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',False],
          'Conf_family_coord_set_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_start': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Comp_index_ID_end': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_start': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Seq_ID_end': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Auth_seq_ID_start': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_end': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Code': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Details': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Secondary_struct_list_ID': [None,returnStarInt,'Secondary_struct_list.ID',True],

                },

        'tagNames': ['ID', 'Conf_family_coord_set_ID', 'Conf_family_coord_set_label', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID_start', 'Comp_index_ID_end', 'Seq_ID_start', 'Seq_ID_end', 'Auth_seq_ID_start', 'Auth_seq_ID_end', 'Name', 'Code', 'Details', 'Entry_ID', 'Secondary_struct_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Secondary_struct_list_ID'],

            }

        },

    'tableNames': ['Secondary_struct_sel', 'Secondary_struct']

    },

  'bond_annotation': {

    'name': 'Bond_annotation_list',

    'tags': {

      'Sf_category': ['bond_annotation',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Bond_annotation': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Bond_type': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
          'Bond_order': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
          'Selection_method_ID': [None,returnStarInt,'Method.ID',False],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_1': [None,returnStarInt,None,False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_2': [None,returnStarInt,None,False],
          'Details': [None,returnStarString,None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Bond_annotation_list_ID': [None,returnStarInt,'Bond_annotation_list.ID',True],

                },

        'tagNames': ['ID', 'Bond_type', 'Bond_order', 'Selection_method_ID', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Atom_isotope_number_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Atom_isotope_number_2', 'Details', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Entry_ID', 'Bond_annotation_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Bond_annotation_list_ID'],

            },

      'Bond_observed_conformer': {

        'tags': {

          'Bond_annotation_ID': [None,returnStarInt,'Bond_annotation.ID',True],
          'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],
          'Atom_site_model_ID': [None,returnStarInt,'Atom_site.Model_ID',True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Bond_annotation_list_ID': [None,returnStarInt,'Bond_annotation_list.ID',True],

                },

        'tagNames': ['Bond_annotation_ID', 'Conformer_family_coord_set_ID', 'Atom_site_model_ID', 'Entry_ID', 'Bond_annotation_list_ID'],
        'sourcePrimaryKeys': ['Bond_annotation_ID', 'Conformer_family_coord_set_ID', 'Atom_site_model_ID', 'Entry_ID', 'Bond_annotation_list_ID'],

            }

        },

    'tableNames': ['Bond_annotation', 'Bond_observed_conformer']

    },

  'structure_interactions': {

    'name': 'Structure_interaction_list',

    'tags': {

      'Sf_category': ['structure_interactions',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Structure_interaction': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Type': [None,returnStarInt,None,True],
          'Selection_method_ID': [None,returnStarInt,'Method.ID',False],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Details': [None,returnStarString,None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Structure_interaction_list_ID': [None,returnStarInt,'Structure_interaction_list.ID',True],

                },

        'tagNames': ['ID', 'Type', 'Selection_method_ID', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Details', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Entry_ID', 'Structure_interaction_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Structure_interaction_list_ID'],

            },

      'Observed_conformer': {

        'tags': {

          'Structure_interaction_ID': [None,returnStarInt,'Structure_interaction.ID',True],
          'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],
          'Atom_site_model_ID': [None,returnStarInt,'Atom_site.Model_ID',True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Structure_interaction_list_ID': [None,returnStarInt,'Structure_interaction_list.ID',True],

                },

        'tagNames': ['Structure_interaction_ID', 'Conformer_family_coord_set_ID', 'Atom_site_model_ID', 'Entry_ID', 'Structure_interaction_list_ID'],
        'sourcePrimaryKeys': ['Structure_interaction_ID', 'Conformer_family_coord_set_ID', 'Atom_site_model_ID', 'Entry_ID', 'Structure_interaction_list_ID'],

            }

        },

    'tableNames': ['Structure_interaction', 'Observed_conformer']

    },

  'other_struct_features': {

    'name': 'Other_struct_feature_list',

    'tags': {

      'Sf_category': ['other_struct_features',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Struct_feature_name': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Struct_feature_definition': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],
      'Selection_method': [None,lambda x = value: returnStarString(x,length = 127),None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Struct_feature_name', 'Struct_feature_definition', 'Details', 'Selection_method'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Other_struct_feature': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Code': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Other_struct_feature_list_ID': [None,returnStarInt,'Other_struct_feature_list.ID',True],

                },

        'tagNames': ['ID', 'Entity_assembly_ID', 'Entity_ID', 'Seq_ID', 'Auth_seq_ID', 'Code', 'Entry_ID', 'Other_struct_feature_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Other_struct_feature_list_ID'],

            }

        },

    'tableNames': ['Other_struct_feature']

    },

  'distance_constraints': {

    'name': 'Distance_constraint_list',

    'tags': {

      'Sf_category': ['distance_constraints',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Constraint_type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Constraint_file_ID': [None,returnStarInt,'Constraint_file.ID',False],
      'Block_ID': [None,returnStarInt,'Constraint_file.Block_ID',False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Constraint_type', 'Constraint_file_ID', 'Block_ID', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Distance_constraint_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Distance_constraint_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Distance_constraint_list_ID'],

            },

      'Dist_constr_software_setting': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Value': [None,returnStarFloat,None,False],
          'Range': [None,lambda x = value: returnStarCode(x,length = 31),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Type', 'Value', 'Range', 'Entry_ID', 'Distance_constraint_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Type', 'Entry_ID', 'Distance_constraint_list_ID'],

            },

      'Distance_constraint_expt': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Method_ID', 'Method_label', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Distance_constraint_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Distance_constraint_list_ID'],

            },

      'Dist_constraint_tree': {

        'tags': {

          'Constraint_ID': [None,returnStarInt,None,True],
          'Node_ID': [None,returnStarInt,None,True],
          'Down_node_ID': [None,returnStarInt,None,False],
          'Right_node_ID': [None,returnStarInt,None,False],
          'Logic_operation': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

                },

        'tagNames': ['Constraint_ID', 'Node_ID', 'Down_node_ID', 'Right_node_ID', 'Logic_operation', 'Entry_ID', 'Distance_constraint_list_ID'],
        'sourcePrimaryKeys': ['Constraint_ID', 'Node_ID', 'Entry_ID', 'Distance_constraint_list_ID'],

            },

      'Dist_constraint': {

        'tags': {

          'Tree_node_member_constraint_ID': [None,returnStarInt,'Dist_constraint_tree.Constraint_ID',True],
          'Tree_node_member_node_ID': [None,returnStarInt,'Dist_constraint_tree.Node_ID',True],
          'Contribution_fractional_val': [None,returnStarFloat,None,False],
          'Constraint_tree_node_member_ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Auth_asym_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

                },

        'tagNames': ['Tree_node_member_constraint_ID', 'Tree_node_member_node_ID', 'Contribution_fractional_val', 'Constraint_tree_node_member_ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Resonance_ID', 'Auth_asym_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'Distance_constraint_list_ID'],
        'sourcePrimaryKeys': ['Tree_node_member_constraint_ID', 'Tree_node_member_node_ID', 'Entry_ID', 'Distance_constraint_list_ID'],

            },

      'Dist_constraint_value': {

        'tags': {

          'Constraint_ID': [None,returnStarInt,'Dist_constraint_tree.Constraint_ID',True],
          'Tree_node_ID': [None,returnStarInt,'Dist_constraint_tree.Node_ID',True],
          'Source_experiment_ID': [None,returnStarInt,'unknown.ID',False],
          'Spectral_peak_ID': [None,returnStarInt,'Peak.ID',False],
          'Intensity_val': [None,returnStarFloat,None,False],
          'Intensity_lower_val_err': [None,returnStarFloat,None,False],
          'Intensity_upper_val_err': [None,returnStarFloat,None,False],
          'Distance_val': [None,returnStarFloat,None,False],
          'Distance_lower_bound_val': [None,returnStarFloat,None,False],
          'Distance_upper_bound_val': [None,returnStarFloat,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

                },

        'tagNames': ['Constraint_ID', 'Tree_node_ID', 'Source_experiment_ID', 'Spectral_peak_ID', 'Intensity_val', 'Intensity_lower_val_err', 'Intensity_upper_val_err', 'Distance_val', 'Distance_lower_bound_val', 'Distance_upper_bound_val', 'Entry_ID', 'Distance_constraint_list_ID'],
        'sourcePrimaryKeys': ['Constraint_ID', 'Tree_node_ID', 'Entry_ID', 'Distance_constraint_list_ID'],

            },

      'Dist_constraint_comment_org': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Comment_text': [None,returnStarString,None,True],
          'Comment_begin_line': [None,returnStarInt,None,False],
          'Comment_begin_column': [None,returnStarInt,None,False],
          'Comment_end_line': [None,returnStarInt,None,False],
          'Comment_end_column': [None,returnStarInt,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Comment_text', 'Comment_begin_line', 'Comment_begin_column', 'Comment_end_line', 'Comment_end_column', 'Entry_ID', 'Distance_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Distance_constraint_list_ID'],

            },

      'Dist_constraint_parse_err': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Content': [None,returnStarString,None,True],
          'Begin_line': [None,returnStarInt,None,False],
          'Begin_column': [None,returnStarInt,None,False],
          'End_line': [None,returnStarInt,None,False],
          'End_column': [None,returnStarInt,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Content', 'Begin_line', 'Begin_column', 'End_line', 'End_column', 'Entry_ID', 'Distance_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Distance_constraint_list_ID'],

            },

      'Dist_constraint_parse_file': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Name', 'Entry_ID', 'Distance_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Distance_constraint_list_ID'],

            },

      'Dist_constraint_conv_err': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Dist_constraint_parse_file_ID': [None,returnStarInt,'Dist_constraint_parse_file.ID',True],
          'Parse_file_constraint_ID': [None,returnStarInt,'Dist_constraint_tree.Constraint_ID',True],
          'Conv_error_type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Conv_error_note': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Dist_constraint_parse_file_ID', 'Parse_file_constraint_ID', 'Conv_error_type', 'Conv_error_note', 'Entry_ID', 'Distance_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Distance_constraint_list_ID'],

            }

        },

    'tableNames': ['Distance_constraint_software', 'Dist_constr_software_setting', 'Distance_constraint_expt', 'Dist_constraint_tree', 'Dist_constraint', 'Dist_constraint_value', 'Dist_constraint_comment_org', 'Dist_constraint_parse_err', 'Dist_constraint_parse_file', 'Dist_constraint_conv_err']

    },

  'floating_chiral_stereo_assign': {

    'name': 'Floating_chirality_assign',

    'tags': {

      'Sf_category': ['floating_chiral_stereo_assign',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],
      'Stereo_count': [None,returnStarInt,None,True],
      'Stereo_assigned_count': [None,returnStarInt,None,True],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Details', 'Stereo_count', 'Stereo_assigned_count'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Floating_chirality_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Floating_chirality_assign_ID': [None,returnStarInt,'Floating_chirality_assign.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Floating_chirality_assign_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Floating_chirality_assign_ID'],

            },

      'Floating_chirality': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Stereospecific_assignment_code': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
          'Auth_asym_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_asym_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Floating_chirality_assign_ID': [None,returnStarInt,'Floating_chirality_assign.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Resonance_ID_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Resonance_ID_2', 'Stereospecific_assignment_code', 'Auth_asym_ID_1', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_asym_ID_2', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Entry_ID', 'Floating_chirality_assign_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Floating_chirality_assign_ID'],

            }

        },

    'tableNames': ['Floating_chirality_software', 'Floating_chirality']

    },

  'torsion_angle_constraints': {

    'name': 'Torsion_angle_constraint_list',

    'tags': {

      'Sf_category': ['torsion_angle_constraints',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Constraint_file_ID': [None,returnStarInt,'Constraint_file.ID',False],
      'Block_ID': [None,returnStarInt,'Constraint_file.Block_ID',False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Constraint_file_ID', 'Block_ID', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Torsion_angle_constraint_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Torsion_angle_constraint_list_ID': [None,returnStarInt,'Torsion_angle_constraint_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],

            },

      'Torsion_angle_constraints_expt': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Torsion_angle_constraint_list_ID': [None,returnStarInt,'Torsion_angle_constraint_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Method_ID', 'Method_label', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],

            },

      'Karplus_equation': {

        'tags': {

          'Torsion_angle_code': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
          'Citation_ID': [None,returnStarInt,'Citation.ID',False],
          'Citation_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Torsion_angle_constraint_list_ID': [None,returnStarInt,'Torsion_angle_constraint_list.ID',True],

                },

        'tagNames': ['Torsion_angle_code', 'Citation_ID', 'Citation_label', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],
        'sourcePrimaryKeys': ['Torsion_angle_code', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],

            },

      'Torsion_angle_constraint': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Torsion_angle_name': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Assembly_atom_ID_3': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_3': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_3': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_3': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_3': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_3': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_3': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_3': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_3': [None,returnStarInt,'Resonance.ID',False],
          'Assembly_atom_ID_4': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_4': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_4': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_4': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_4': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_4': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_4': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_4': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_4': [None,returnStarInt,'Resonance.ID',False],
          'Angle_lower_bound_val': [None,returnStarFloat,None,True],
          'Angle_upper_bound_val': [None,returnStarFloat,None,True],
          'Source_experiment_ID': [None,returnStarInt,'unknown.ID',False],
          'Auth_asym_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_asym_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_asym_ID_3': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_3': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_3': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_3': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_asym_ID_4': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_4': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_4': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_4': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Torsion_angle_constraint_list_ID': [None,returnStarInt,'Torsion_angle_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Torsion_angle_name', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Resonance_ID_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Resonance_ID_2', 'Assembly_atom_ID_3', 'Entity_assembly_ID_3', 'Entity_ID_3', 'Comp_index_ID_3', 'Seq_ID_3', 'Comp_ID_3', 'Atom_ID_3', 'Atom_type_3', 'Resonance_ID_3', 'Assembly_atom_ID_4', 'Entity_assembly_ID_4', 'Entity_ID_4', 'Comp_index_ID_4', 'Seq_ID_4', 'Comp_ID_4', 'Atom_ID_4', 'Atom_type_4', 'Resonance_ID_4', 'Angle_lower_bound_val', 'Angle_upper_bound_val', 'Source_experiment_ID', 'Auth_asym_ID_1', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_asym_ID_2', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Auth_asym_ID_3', 'Auth_seq_ID_3', 'Auth_comp_ID_3', 'Auth_atom_ID_3', 'Auth_asym_ID_4', 'Auth_seq_ID_4', 'Auth_comp_ID_4', 'Auth_atom_ID_4', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],

            },

      'TA_constraint_comment_org': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Comment_text': [None,returnStarString,None,True],
          'Comment_begin_line': [None,returnStarInt,None,False],
          'Comment_begin_column': [None,returnStarInt,None,False],
          'Comment_end_line': [None,returnStarInt,None,False],
          'Comment_end_column': [None,returnStarInt,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Torsion_angle_constraint_list_ID': [None,returnStarInt,'Torsion_angle_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Comment_text', 'Comment_begin_line', 'Comment_begin_column', 'Comment_end_line', 'Comment_end_column', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],

            },

      'TA_constraint_parse_err': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Content': [None,returnStarString,None,True],
          'Begin_line': [None,returnStarInt,None,False],
          'Begin_column': [None,returnStarInt,None,False],
          'End_line': [None,returnStarInt,None,False],
          'End_column': [None,returnStarInt,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Torsion_angle_constraint_list_ID': [None,returnStarInt,'Torsion_angle_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Content', 'Begin_line', 'Begin_column', 'End_line', 'End_column', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],

            },

      'TA_constraint_parse_file': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Torsion_angle_constraint_list_ID': [None,returnStarInt,'Torsion_angle_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Name', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],

            },

      'TA_constraint_conv_err': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Parse_file_ID': [None,returnStarInt,'TA_constraint_parse_file.ID',True],
          'Parse_file_constraint_ID': [None,returnStarInt,'Torsion_angle_constraint.ID',True],
          'Conv_error_type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Conv_error_note': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Torsion_angle_constraint_list_ID': [None,returnStarInt,'Torsion_angle_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Parse_file_ID', 'Parse_file_constraint_ID', 'Conv_error_type', 'Conv_error_note', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'Torsion_angle_constraint_list_ID'],

            }

        },

    'tableNames': ['Torsion_angle_constraint_software', 'Torsion_angle_constraints_expt', 'Karplus_equation', 'Torsion_angle_constraint', 'TA_constraint_comment_org', 'TA_constraint_parse_err', 'TA_constraint_parse_file', 'TA_constraint_conv_err']

    },

  'RDC_constraints': {

    'name': 'RDC_constraint_list',

    'tags': {

      'Sf_category': ['RDC_constraints',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Constraint_file_ID': [None,returnStarInt,'Constraint_file.ID',False],
      'Block_ID': [None,returnStarInt,'Constraint_file.Block_ID',False],
      'Dipolar_constraint_calib_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Mol_align_tensor_axial_sym_mol': [None,returnStarFloat,None,False],
      'Mol_align_tensor_rhombic_mol': [None,returnStarFloat,None,False],
      'General_order_param_int_motions': [None,returnStarFloat,None,False],
      'Assumed_H_N_bond_length': [None,returnStarFloat,None,False],
      'Assumed_H_C_bond_length': [None,returnStarFloat,None,False],
      'Assumed_C_N_bond_length': [None,returnStarFloat,None,False],
      'Data_file_format': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Constraint_file_ID', 'Block_ID', 'Dipolar_constraint_calib_method', 'Mol_align_tensor_axial_sym_mol', 'Mol_align_tensor_rhombic_mol', 'General_order_param_int_motions', 'Assumed_H_N_bond_length', 'Assumed_H_C_bond_length', 'Assumed_C_N_bond_length', 'Data_file_format', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'RDC_constraint_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'RDC_constraint_list_ID': [None,returnStarInt,'RDC_constraint_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'RDC_constraint_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'RDC_constraint_list_ID'],

            },

      'RDC_constraint_expt': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'RDC_constraint_list_ID': [None,returnStarInt,'RDC_constraint_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Method_ID', 'Method_label', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'RDC_constraint_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'RDC_constraint_list_ID'],

            },

      'RDC_constraint': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_1': [None,returnStarInt,None,False],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number_2': [None,returnStarInt,None,False],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'RDC_val': [None,returnStarFloat,None,True],
          'RDC_lower_bound': [None,returnStarFloat,None,False],
          'RDC_upper_bound': [None,returnStarFloat,None,False],
          'RDC_val_err': [None,returnStarFloat,None,False],
          'RDC_val_scale_factor': [None,returnStarFloat,None,False],
          'Source_experiment_ID': [None,returnStarInt,'unknown.ID',False],
          'Auth_asym_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_asym_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'RDC_constraint_list_ID': [None,returnStarInt,'RDC_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Atom_isotope_number_1', 'Resonance_ID_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Atom_isotope_number_2', 'Resonance_ID_2', 'RDC_val', 'RDC_lower_bound', 'RDC_upper_bound', 'RDC_val_err', 'RDC_val_scale_factor', 'Source_experiment_ID', 'Auth_asym_ID_1', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_asym_ID_2', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Entry_ID', 'RDC_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'RDC_constraint_list_ID'],

            },

      'RDC_constraint_comment_org': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Comment_text': [None,returnStarString,None,True],
          'Comment_begin_line': [None,returnStarInt,None,False],
          'Comment_begin_column': [None,returnStarInt,None,False],
          'Comment_end_line': [None,returnStarInt,None,False],
          'Comment_end_column': [None,returnStarInt,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'RDC_constraint_list_ID': [None,returnStarInt,'RDC_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Comment_text', 'Comment_begin_line', 'Comment_begin_column', 'Comment_end_line', 'Comment_end_column', 'Entry_ID', 'RDC_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'RDC_constraint_list_ID'],

            },

      'RDC_constraint_parse_err': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Content': [None,returnStarString,None,True],
          'Begin_line': [None,returnStarInt,None,False],
          'Begin_column': [None,returnStarInt,None,False],
          'End_line': [None,returnStarInt,None,False],
          'End_column': [None,returnStarInt,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'RDC_constraint_list_ID': [None,returnStarInt,'RDC_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Content', 'Begin_line', 'Begin_column', 'End_line', 'End_column', 'Entry_ID', 'RDC_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'RDC_constraint_list_ID'],

            },

      'RDC_constraint_parse_file': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Name': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'RDC_constraint_list_ID': [None,returnStarInt,'RDC_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Name', 'Entry_ID', 'RDC_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'RDC_constraint_list_ID'],

            },

      'RDC_constraint_conv_err': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'RDC_constr_parse_file_ID': [None,returnStarInt,'RDC_constraint_parse_file.ID',True],
          'Parse_file_constraint_ID': [None,returnStarInt,'RDC_constraint.ID',True],
          'Conv_error_type': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Conv_error_note': [None,returnStarString,None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'RDC_constraint_list_ID': [None,returnStarInt,'RDC_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'RDC_constr_parse_file_ID', 'Parse_file_constraint_ID', 'Conv_error_type', 'Conv_error_note', 'Entry_ID', 'RDC_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'RDC_constraint_list_ID'],

            }

        },

    'tableNames': ['RDC_constraint_software', 'RDC_constraint_expt', 'RDC_constraint', 'RDC_constraint_comment_org', 'RDC_constraint_parse_err', 'RDC_constraint_parse_file', 'RDC_constraint_conv_err']

    },

  'J_three_bond_constraints': {

    'name': 'J_three_bond_constraint_list',

    'tags': {

      'Sf_category': ['J_three_bond_constraints',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Data_file_format': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Constraint_file_ID': [None,returnStarInt,'Constraint_file.ID',False],
      'Block_ID': [None,returnStarInt,'Constraint_file.Block_ID',False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Data_file_format', 'Constraint_file_ID', 'Block_ID', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'J_three_bond_constraint_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'J_three_bond_constraint_list_ID': [None,returnStarInt,'J_three_bond_constraint_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'J_three_bond_constraint_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'J_three_bond_constraint_list_ID'],

            },

      'J_three_bond_constraint_expt': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'J_three_bond_constraint_list_ID': [None,returnStarInt,'J_three_bond_constraint_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Method_ID', 'Method_label', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'J_three_bond_constraint_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'J_three_bond_constraint_list_ID'],

            },

      'J_three_bond_constraint': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Assembly_atom_ID_3': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_3': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_3': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_3': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_3': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_3': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_3': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_3': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_3': [None,returnStarInt,'Resonance.ID',False],
          'Assembly_atom_ID_4': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_4': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_4': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_4': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_4': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_4': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_4': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_4': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_4': [None,returnStarInt,'Resonance.ID',False],
          'Coupling_constant_val': [None,returnStarFloat,None,False],
          'Coupling_constant_lower_bound': [None,returnStarFloat,None,False],
          'Coupling_constant_upper_bound': [None,returnStarFloat,None,False],
          'Coupling_constant_err': [None,returnStarFloat,None,False],
          'Source_experiment_ID': [None,returnStarInt,'unknown.ID',False],
          'Auth_asym_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_asym_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_asym_ID_3': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_3': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_3': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_3': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_asym_ID_4': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_4': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_4': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_4': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'J_three_bond_constraint_list_ID': [None,returnStarInt,'J_three_bond_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Resonance_ID_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Resonance_ID_2', 'Assembly_atom_ID_3', 'Entity_assembly_ID_3', 'Entity_ID_3', 'Comp_index_ID_3', 'Seq_ID_3', 'Comp_ID_3', 'Atom_ID_3', 'Atom_type_3', 'Resonance_ID_3', 'Assembly_atom_ID_4', 'Entity_assembly_ID_4', 'Entity_ID_4', 'Comp_index_ID_4', 'Seq_ID_4', 'Comp_ID_4', 'Atom_ID_4', 'Atom_type_4', 'Resonance_ID_4', 'Coupling_constant_val', 'Coupling_constant_lower_bound', 'Coupling_constant_upper_bound', 'Coupling_constant_err', 'Source_experiment_ID', 'Auth_asym_ID_1', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_asym_ID_2', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Auth_asym_ID_3', 'Auth_seq_ID_3', 'Auth_comp_ID_3', 'Auth_atom_ID_3', 'Auth_asym_ID_4', 'Auth_seq_ID_4', 'Auth_comp_ID_4', 'Auth_atom_ID_4', 'Entry_ID', 'J_three_bond_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'J_three_bond_constraint_list_ID'],

            }

        },

    'tableNames': ['J_three_bond_constraint_software', 'J_three_bond_constraint_expt', 'J_three_bond_constraint']

    },

  'CA_CB_chem_shift_constraints': {

    'name': 'CA_CB_constraint_list',

    'tags': {

      'Sf_category': ['CA_CB_chem_shift_constraints',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Data_file_format': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Constraint_file_ID': [None,returnStarInt,'Constraint_file.ID',False],
      'Block_ID': [None,returnStarInt,'Constraint_file.Block_ID',False],
      'Units': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Data_file_name', 'Data_file_format', 'Constraint_file_ID', 'Block_ID', 'Units', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'CA_CB_constraint_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'CA_CB_constraint_list_ID': [None,returnStarInt,'CA_CB_constraint_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'CA_CB_constraint_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'CA_CB_constraint_list_ID'],

            },

      'CA_CB_constraint_expt': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'CA_CB_constraint_list_ID': [None,returnStarInt,'CA_CB_constraint_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Method_ID', 'Method_label', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'CA_CB_constraint_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'CA_CB_constraint_list_ID'],

            },

      'CA_CB_constraint': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_1': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_1': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_1': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_1': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_1': [None,returnStarInt,'Resonance.ID',False],
          'Assembly_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_2': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_2': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_2': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_2': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_2': [None,returnStarInt,'Resonance.ID',False],
          'Assembly_atom_ID_3': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_3': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_3': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_3': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_3': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_3': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_3': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_3': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_3': [None,returnStarInt,'Resonance.ID',False],
          'Assembly_atom_ID_4': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_4': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_4': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_4': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_4': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_4': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_4': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_4': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_4': [None,returnStarInt,'Resonance.ID',False],
          'Assembly_atom_ID_5': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID_5': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID_5': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID_5': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID_5': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID_5': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID_5': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type_5': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Resonance_ID_5': [None,returnStarInt,'Resonance.ID',False],
          'CA_chem_shift_val': [None,returnStarFloat,None,True],
          'CA_chem_shift_val_err': [None,returnStarFloat,None,True],
          'CB_chem_shift_val': [None,returnStarFloat,None,True],
          'CB_chem_shift_val_err': [None,returnStarFloat,None,True],
          'Source_experiment_ID': [None,returnStarInt,'unknown.ID',False],
          'Auth_asym_ID_1': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_asym_ID_2': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_asym_ID_3': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_3': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_3': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_3': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_asym_ID_4': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_4': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_4': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_4': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_asym_ID_5': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID_5': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID_5': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID_5': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'CA_CB_constraint_list_ID': [None,returnStarInt,'CA_CB_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID_1', 'Entity_assembly_ID_1', 'Entity_ID_1', 'Comp_index_ID_1', 'Seq_ID_1', 'Comp_ID_1', 'Atom_ID_1', 'Atom_type_1', 'Resonance_ID_1', 'Assembly_atom_ID_2', 'Entity_assembly_ID_2', 'Entity_ID_2', 'Comp_index_ID_2', 'Seq_ID_2', 'Comp_ID_2', 'Atom_ID_2', 'Atom_type_2', 'Resonance_ID_2', 'Assembly_atom_ID_3', 'Entity_assembly_ID_3', 'Entity_ID_3', 'Comp_index_ID_3', 'Seq_ID_3', 'Comp_ID_3', 'Atom_ID_3', 'Atom_type_3', 'Resonance_ID_3', 'Assembly_atom_ID_4', 'Entity_assembly_ID_4', 'Entity_ID_4', 'Comp_index_ID_4', 'Seq_ID_4', 'Comp_ID_4', 'Atom_ID_4', 'Atom_type_4', 'Resonance_ID_4', 'Assembly_atom_ID_5', 'Entity_assembly_ID_5', 'Entity_ID_5', 'Comp_index_ID_5', 'Seq_ID_5', 'Comp_ID_5', 'Atom_ID_5', 'Atom_type_5', 'Resonance_ID_5', 'CA_chem_shift_val', 'CA_chem_shift_val_err', 'CB_chem_shift_val', 'CB_chem_shift_val_err', 'Source_experiment_ID', 'Auth_asym_ID_1', 'Auth_seq_ID_1', 'Auth_comp_ID_1', 'Auth_atom_ID_1', 'Auth_asym_ID_2', 'Auth_seq_ID_2', 'Auth_comp_ID_2', 'Auth_atom_ID_2', 'Auth_asym_ID_3', 'Auth_seq_ID_3', 'Auth_comp_ID_3', 'Auth_atom_ID_3', 'Auth_asym_ID_4', 'Auth_seq_ID_4', 'Auth_comp_ID_4', 'Auth_atom_ID_4', 'Auth_asym_ID_5', 'Auth_seq_ID_5', 'Auth_comp_ID_5', 'Auth_atom_ID_5', 'Entry_ID', 'CA_CB_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'CA_CB_constraint_list_ID'],

            }

        },

    'tableNames': ['CA_CB_constraint_software', 'CA_CB_constraint_expt', 'CA_CB_constraint']

    },

  'H_chem_shift_constraints': {

    'name': 'H_chem_shift_constraint_list',

    'tags': {

      'Sf_category': ['H_chem_shift_constraints',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Units': [None,lambda x = value: returnStarCode(x,length = 31),None,True],
      'Data_file_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Data_file_format': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Constraint_file_ID': [None,returnStarInt,'Constraint_file.ID',False],
      'Block_ID': [None,returnStarInt,'Constraint_file.Block_ID',False],
      'Details': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Units', 'Data_file_name', 'Data_file_format', 'Constraint_file_ID', 'Block_ID', 'Details'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'H_chem_shift_constraint_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'H_chem_shift_constraint_list_ID': [None,returnStarInt,'H_chem_shift_constraint_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'H_chem_shift_constraint_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'H_chem_shift_constraint_list_ID'],

            },

      'H_chem_shift_constraint_expt': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'H_chem_shift_constraint_list_ID': [None,returnStarInt,'H_chem_shift_constraint_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Method_ID', 'Method_label', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'H_chem_shift_constraint_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'H_chem_shift_constraint_list_ID'],

            },

      'H_chem_shift_constraint': {

        'tags': {

          'ID': [None,returnStarInt,None,True],
          'Assembly_atom_ID': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
          'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
          'Entity_ID': [None,returnStarInt,'Entity.ID',True],
          'Comp_index_ID': [None,returnStarInt,'Entity_comp_index.ID',True],
          'Seq_ID': [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
          'Comp_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
          'Atom_ID': [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
          'Atom_type': [None,lambda x = value: returnStarCode(x,length = 15),None,True],
          'Atom_isotope_number': [None,returnStarInt,None,False],
          'Resonance_ID': [None,returnStarInt,'Resonance.ID',False],
          'Chem_shift_val': [None,returnStarFloat,None,True],
          'Chem_shift_val_err': [None,returnStarFloat,None,True],
          'Source_experiment_ID': [None,returnStarInt,'unknown.ID',False],
          'Auth_asym_ID': [None,lambda x = value: returnStarCode(x,length = 12),None,False],
          'Auth_seq_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_comp_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Auth_atom_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'H_chem_shift_constraint_list_ID': [None,returnStarInt,'H_chem_shift_constraint_list.ID',True],

                },

        'tagNames': ['ID', 'Assembly_atom_ID', 'Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Atom_isotope_number', 'Resonance_ID', 'Chem_shift_val', 'Chem_shift_val_err', 'Source_experiment_ID', 'Auth_asym_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Entry_ID', 'H_chem_shift_constraint_list_ID'],
        'sourcePrimaryKeys': ['ID', 'Entry_ID', 'H_chem_shift_constraint_list_ID'],

            }

        },

    'tableNames': ['H_chem_shift_constraint_software', 'H_chem_shift_constraint_expt', 'H_chem_shift_constraint']

    },

  'other_constraints': {

    'name': 'Other_constraint_list',

    'tags': {

      'Sf_category': ['other_constraints',lambda x = value: returnStarCode(x,length = 31),None,True],
      'Sf_framecode': [None,lambda x = value: returnStarCode(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Details': [None,returnStarString,None,False],
      'Entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
      'Entity_ID': [None,returnStarInt,'Entity.ID',True],
      'Type': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Subtype': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'File_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'File_format': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
      'Constraint_file_ID': [None,returnStarInt,'Constraint_file.ID',False],
      'Block_ID': [None,returnStarInt,'Constraint_file.Block_ID',False],
      'Text': [None,returnStarString,None,True],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Details', 'Entity_assembly_ID', 'Entity_ID', 'Type', 'Subtype', 'File_name', 'File_format', 'Constraint_file_ID', 'Block_ID', 'Text'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    'tables': {

      'Other_constraint_expt': {

        'tags': {

          'Experiment_ID': [None,returnStarInt,'unknown.ID',True],
          'Experiment_name': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Sample_ID': [None,returnStarInt,'Sample.ID',True],
          'Sample_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Sample_state': [None,lambda x = value: returnStarLine(x,length = 31),None,True],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Other_constraint_list_ID': [None,returnStarInt,'Other_constraint_list.ID',True],

                },

        'tagNames': ['Experiment_ID', 'Experiment_name', 'Method_ID', 'Method_label', 'Sample_ID', 'Sample_label', 'Sample_state', 'Entry_ID', 'Other_constraint_list_ID'],
        'sourcePrimaryKeys': ['Experiment_ID', 'Sample_ID', 'Sample_state', 'Entry_ID', 'Other_constraint_list_ID'],

            },

      'Other_constraint_software': {

        'tags': {

          'Software_ID': [None,returnStarInt,'Software.ID',True],
          'Software_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
          'Method_ID': [None,returnStarInt,'Method.ID',False],
          'Method_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
          'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
          'Other_constraint_list_ID': [None,returnStarInt,'Other_constraint_list.ID',True],

                },

        'tagNames': ['Software_ID', 'Software_label', 'Method_ID', 'Method_label', 'Entry_ID', 'Other_constraint_list_ID'],
        'sourcePrimaryKeys': ['Software_ID', 'Entry_ID', 'Other_constraint_list_ID'],

            }

        },

    'tableNames': ['Other_constraint_expt', 'Other_constraint_software']

    },

  'org_constr_file_comment': {

    'name': 'Org_constr_file_comment',

    'tags': {

      'Sf_category': ['org_constr_file_comment',lambda x = value: returnStarLine(x,length = 127),None,False],
      'Sf_framecode': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
      'Entry_ID': [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
      'ID': [None,returnStarInt,None,True],
      'Constraint_file_ID': [None,returnStarInt,'Constraint_file.ID',False],
      'Block_ID': [None,returnStarInt,'Constraint_file.Block_ID',False],
      'Details': [None,returnStarString,None,False],
      'Comment': [None,returnStarString,None,False],

            },

    'tagNames': ['Sf_category', 'Sf_framecode', 'Entry_ID', 'ID', 'Constraint_file_ID', 'Block_ID', 'Details', 'Comment'],
    'sourcePrimaryKeys': ['Entry_ID', 'ID'],

    }
}

sfList = ['entry_interview', 'deposited_data_files', 'study_list', 'entry_information', 'citations', 'assembly', 'assembly_annotation', 'assembly_subsystems', 'entity', 'natural_source', 'experimental_source', 'chem_comp', 'sample', 'sample_conditions', 'molecule_purity', 'software', 'method', 'Mass_spectrometer', 'Mass_spectrometer_list', 'Mass_spec_ref_compd', 'chromatographic_system', 'chromatographic_system', 'NMR_spectrometer', 'NMR_spectrometer_list', 'NMR_spectrometer_probe', 'experiment_list', 'NMR_spectrometer_expt', 'NMR_spectral_processing', 'computer', 'chem_shift_reference', 'assigned_chemical_shifts', 'coupling_constants', 'spectral_peak_list', 'resonance_linker', 'chem_shift_isotope_effect', 'chem_shift_interaction_diff', 'chem_shift_anisotropy', 'chem_shifts_calc_type', 'shielding_tensors', 'theoretical_chem_shifts', 'RDCs', 'dipolar_couplings', 'spectral_density_values', 'other_data_types', 'H_exch_rates', 'H_exch_protection_factors', 'homonucl_NOEs', 'heteronucl_NOEs', 'heteronucl_T1_relaxation', 'heteronucl_T1rho_relaxation', 'heteronucl_T2_relaxation', 'dipole_dipole_relaxation', 'dipole_dipole_cross_correlations', 'dipole_CSA_cross_correlations', 'order_parameters', 'pH_titration', 'pH_param_list', 'D_H_fractionation_factors', 'deduced_secd_struct_features', 'deduced_hydrogen_bonds', 'conformer_statistics', 'conformer_family_coord_set', 'representative_conformer', 'constraint_statistics', 'force_constants', 'angular_order_parameters', 'tertiary_struct_elements', 'secondary_structs', 'bond_annotation', 'structure_interactions', 'other_struct_features', 'distance_constraints', 'floating_chiral_stereo_assign', 'torsion_angle_constraints', 'RDC_constraints', 'J_three_bond_constraints', 'CA_CB_chem_shift_constraints', 'H_chem_shift_constraints', 'other_constraints', 'org_constr_file_comment']



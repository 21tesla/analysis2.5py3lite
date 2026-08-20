aditMandDict = {

  'assembly': {

    'tags': {

      'Name': ['', '', '', 'Molecular assembly name '],

      },

    'tagNames': ['Name'],

    'tables': {

      'Entity_assembly': {

        'tags': {

          'Entity_assembly_name': ['', '', '', 'Assembly component name'],
          'Experimental_data_reported': ['', '', 'Y', 'Experimental data reported'],
          'Conformational_isomer': ['', '', 'Y', 'NMR observable conformational isomer '],
          'Chemical_exchange_state': ['', '', 'Y', 'Chemical exchange state'],
          'Physical_state': ['', '', '', 'Physical state (biopolymers only) '],
          'Entity_label': ['', '', '', 'Molecular entity name to assign to this assembly component '],

          },

        'tagNames': ['Chemical_exchange_state', 'Conformational_isomer', 'Entity_assembly_name', 'Entity_label', 'Experimental_data_reported', 'Physical_state'],

        },
      },

    'tableNames': ['Entity_assembly'],

    },

  'experiment_list': {

    'tables': {

      'Experiment': {

        'tags': {

          'Sample_label': ['', '', '', 'Sample label '],
          'Sample_condition_list_label': ['', '', '', 'Sample conditions label '],
          'Sample_state': ['', '', '', 'Sample state'],
          'Name': ['', '', '', 'NMR experiment name'],
          'Raw_data_flag': ['', '', 'Y', 'Time-domain data provided'],

          },

        'tagNames': ['Name', 'Raw_data_flag', 'Sample_condition_list_label', 'Sample_label', 'Sample_state'],

        },
      },

    'tableNames': ['Experiment'],

    },

  'NMR_spectrometer_list': {

    'tables': {

      'NMR_spectrometer_view': {

        'tags': {

          'Field_strength': ['', '', '', 'Field strength (MHz)'],
          'Model': ['', '', '', 'Model '],
          'Name': ['', '', '', 'Name for the spectrometer'],
          'Manufacturer': ['', '', '', 'Manufacturer'],

          },

        'tagNames': ['Field_strength', 'Manufacturer', 'Model', 'Name'],

        },
      },

    'tableNames': ['NMR_spectrometer_view'],

    },

  'entity': {

    'tags': {

      'Ambiguous_chem_comp_sites': ['', '', 'Y', 'Ambiguous residue sites exist '],
      'Ambiguous_conformational_states': ['', '', 'Y', 'Ambiguous conformational states are present '],
      'Nstd_chirality': ['', '', 'Y', 'Non-standard chirality exists in the molecule'],
      'Nstd_linkage': ['', '', 'Y', 'A non-standard bond exists in the molecule'],
      'Nstd_monomer': ['', '', 'Y', 'A non-standard monomer is present'],
      'Paramagnetic': ['', '', 'Y', 'Is the molecule paramagnetic'],
      'Polymer_seq_one_letter_code': ['', '', '', 'Enter the complete polymer sequence (one-letter code)'],
      'Polymer_type': ['', '', 'Y', 'Polymer type'],
      'Sf_framecode': ['', '', '', 'Brief molecule name or acronym (label)'],
      'Thiol_state': ['', '', '', 'What is the state of the thiol groups in the molecule'],
      'Type': ['', '', 'Y', 'Molecule type'],

      },

    'tagNames': ['Ambiguous_chem_comp_sites', 'Ambiguous_conformational_states', 'Nstd_chirality', 'Nstd_linkage', 'Nstd_monomer', 'Paramagnetic', 'Polymer_seq_one_letter_code', 'Polymer_type', 'Sf_framecode', 'Thiol_state', 'Type'],

    },

  'sample': {

    'tags': {

      'Sf_framecode': ['', '', '', 'Sample label'],
      'Solvent_system': ['', '', '', 'Solvent system'],
      'Type': ['', '', '', 'Sample type '],

      },

    'tagNames': ['Sf_framecode', 'Solvent_system', 'Type'],

    'tables': {

      'Sample_component': {

        'tags': {

          'Concentration_val_units': ['', '', '', 'Concentration value units'],
          'Isotopic_labeling': ['', '', '', 'Isotopic labeling '],
          'Mol_common_name': ['', '', '', 'Molecule common name'],
          'Concentration_val': ['', '', '', 'Concentration value (or range)'],

          },

        'tagNames': ['Concentration_val', 'Concentration_val_units', 'Isotopic_labeling', 'Mol_common_name'],

        },
      },

    'tableNames': ['Sample_component'],

    },

  'sample_conditions': {

    'tags': {

      'Sf_framecode': ['', '', '', 'Sample conditions label'],

      },

    'tagNames': ['Sf_framecode'],

    'tables': {

      'Sample_condition_variable': {

        'tags': {

          'Val_units': ['', '', '', 'Variable value units'],
          'Type': ['', '', 'Y', 'Variable type'],
          'Val': ['', '', '', 'Variable value'],

          },

        'tagNames': ['Type', 'Val', 'Val_units'],

        },
      },

    'tableNames': ['Sample_condition_variable'],

    },

  'citations': {

    'tags': {

      'Class': ['', '', 'Y', 'Entry or reference citation'],
      'Journal_abbrev': ['Y', "journal' value in _Citation.Type", '', 'Journal abbreviation'],
      'Sf_framecode': ['', '', '', 'Citation label'],
      'Status': ['', '', 'Y', 'Citation status'],
      'Title': ['', '', '', 'Citation title'],
      'Type': ['', '', 'Y', 'Citation type'],

      },

    'tagNames': ['Class', 'Journal_abbrev', 'Sf_framecode', 'Status', 'Title', 'Type'],

    'tables': {

      'Citation_author': {

        'tags': {

          'Family_name': ['', '', '', "Author's last name"],
          'Given_name': ['', '', '', "Author's first name"],

          },

        'tagNames': ['Family_name', 'Given_name'],

        },
      },

    'tableNames': ['Citation_author'],

    },

  'entry_interview': {

    'tags': {

      'Assigned_chem_shifts': ['', '', 'Y', 'Assigned NMR chemical shifts '],
      'Chem_shift_anisotropy': ['', '', 'Y', 'Chemical shift anisotropy '],
      'Coupling_constants': ['', '', 'Y', 'Scalar coupling constants '],
      'DD_cross_correlation': ['', '', 'Y', 'Dipole-dipole cross correlation data'],
      'D_H_fractionation_factors': ['', '', 'Y', 'NMR-derived D/H fractionation factors '],
      'Dipole_CSA_cross_correlation': ['', '', 'Y', 'Dipole-CSA cross correlation data'],
      'Dipole_dipole_relaxation': ['', '', 'Y', 'Dipole-dipole relaxation data '],
      'H_exchange_protection_factors': ['', '', 'Y', 'Hydrogen exchange protection data'],
      'H_exchange_rate': ['', '', 'Y', 'Hydrogen exchange rates'],
      'Heteronucl_NOEs': ['', '', 'Y', 'Heteronuclear NOEs'],
      'Heteronucl_T1_relaxation': ['', '', 'Y', 'T1 (R1) NMR relaxation data'],
      'Heteronucl_T1rho_relaxation': ['', '', 'Y', 'T1rho (R1rho) NMR relaxation data '],
      'Heteronucl_T2_relaxation': ['', '', 'Y', 'T2 (R2) NMR relaxation data'],
      'Homonucl_NOEs': ['', '', 'Y', '"Homonuclear NOEs (not for NOE restraints, restraints are uploaded as part of a PDB deposition) "'],
      'Ligands': ['', '', 'Y', 'The system contains a metal and/or a non-polymeric organic ligand '],
      'Non_standard_residues': ['', '', 'Y', 'At least one biopolymer contains a non-standard residue '],
      'Order_parameters': ['', '', 'Y', 'Order parameters'],
      'Other_kind_of_data': ['', '', 'Y', 'Other kinds of data'],
      'PKa_value_data_set': ['', '', 'Y', "NMR-derived pH transitions (pKa's; pHmid's) "],
      'Residual_dipolar_couplings': ['', '', 'Y', 'Residual dipolar couplings '],
      'Spectral_density_values': ['', '', 'Y', 'Spectral density factors'],
      'Spectral_peak_lists': ['', '', 'Y', 'Spectral peak lists '],
      'Theoretical_chem_shifts': ['', '', 'Y', 'Theoretical (calculated) chemical shift values '],

      },

    'tagNames': ['Assigned_chem_shifts', 'Chem_shift_anisotropy', 'Coupling_constants', 'DD_cross_correlation', 'D_H_fractionation_factors', 'Dipole_CSA_cross_correlation', 'Dipole_dipole_relaxation', 'H_exchange_protection_factors', 'H_exchange_rate', 'Heteronucl_NOEs', 'Heteronucl_T1_relaxation', 'Heteronucl_T1rho_relaxation', 'Heteronucl_T2_relaxation', 'Homonucl_NOEs', 'Ligands', 'Non_standard_residues', 'Order_parameters', 'Other_kind_of_data', 'PKa_value_data_set', 'Residual_dipolar_couplings', 'Spectral_density_values', 'Spectral_peak_lists', 'Theoretical_chem_shifts'],

    },

  'experimental_source': {

    'tables': {

      'Entity_experimental_src': {

        'tags': {

          'Vector_type': ['', '', '', 'Type of vector used (if applicable) '],
          'Entity_label': ['', '', '', 'Label for a molecular entity '],
          'Host_org_scientific_name': ['', '', '', "Host organism's scientific name"],
          'Vector_name': ['', '', '', 'Vector name (if applicable) '],
          'Production_method': ['', '', '', 'Production method '],

          },

        'tagNames': ['Entity_label', 'Host_org_scientific_name', 'Production_method', 'Vector_name', 'Vector_type'],

        },
      },

    'tableNames': ['Entity_experimental_src'],

    },

  'entry_information': {

    'tags': {

      'Dep_release_code_nmr_exptl': ['', '', 'Y', 'Release status for NMR experiment'],
      'Dep_release_code_sequence': ['', '', 'Y', 'Prerelease status for polymer sequences'],
      'Experimental_method_subtype': ['', '', '', 'Type of NMR method'],
      'Title': ['', '', '', 'Entry title'],

      },

    'tagNames': ['Dep_release_code_nmr_exptl', 'Dep_release_code_sequence', 'Experimental_method_subtype', 'Title'],

    'tables': {

      'Contact_person': {

        'tags': {

          'State_province': ['', '', '', 'State or province'],
          'City': ['', '', '', 'City'],
          'Family_name': ['', '', '', "Contact person's last name"],
          'Country': ['', '', '', 'Country'],
          'Phone_number': ['', '', '', 'Telephone number'],
          'Postal_code': ['', '', '', 'Postal code'],
          'Given_name': ['', '', '', "Contact person's first name"],
          'Organization_type': ['', '', 'Y', "Contact person's organization type"],
          'Role': ['', '', 'Y', 'Role of this person in the project'],
          'Email_address': ['', '', '', 'E-mail address'],
          'Address_1': ['', '', '', 'Mailing address (line 1)'],

          },

        'tagNames': ['Address_1', 'City', 'Country', 'Email_address', 'Family_name', 'Given_name', 'Organization_type', 'Phone_number', 'Postal_code', 'Role', 'State_province'],

        },

      'Entry_author': {

        'tags': {

          'Family_name': ['', '', '', "Entry author's last name"],
          'Given_name': ['', '', '', "Entry author's first name"],

          },

        'tagNames': ['Family_name', 'Given_name'],

        },
      },

    'tableNames': ['Contact_person', 'Entry_author'],

    },

  'natural_source': {

    'tables': {

      'Entity_natural_src': {

        'tags': {

          'Entity_label': ['', '', '', 'Label assigned to the molecular entity'],
          'Type': ['', '', '', 'Type of natural source'],
          'Organism_name_scientific': ['', '', '', '"Scientific name for the natural source (i.e., genus and species) "'],

          },

        'tagNames': ['Entity_label', 'Organism_name_scientific', 'Type'],

        },
      },

    'tableNames': ['Entity_natural_src'],

    },

  'chem_comp': {

    'tags': {

      'Aromatic': ['Y', "'yes' value in _Entry_interview.Ligands OR 'yes' value in _Entry_interview.Non_standard_residues", '', 'Is the compound aromatic'],
      'Name': ['Y', "yes' value in _Entry_interview.Ligands OR 'yes' value in _Entry_interview.Non_standard_residues", '', 'Ligand or chemical component name'],
      'Paramagnetic': ['Y', "'yes' value in _Entry_interview.Ligands OR 'yes' value in _Entry_interview.Non_standard_residues", '', 'Is the compound paramagnetic'],

      },

    'tagNames': ['Aromatic', 'Name', 'Paramagnetic'],

    },

  'software': {

    'tags': {

      'Sf_framecode': ['', '', '', 'Software name'],

      },

    'tagNames': ['Sf_framecode'],

    'tables': {

      'Task': {

        'tags': {

          'Task': ['', '', '', 'Task '],

          },

        'tagNames': ['Task'],

        },

      'Vendor': {

        'tags': {

          'Name': ['', '', '', 'Author (vendor) '],

          },

        'tagNames': ['Name'],

        },
      },

    'tableNames': ['Task', 'Vendor'],

    },

  }

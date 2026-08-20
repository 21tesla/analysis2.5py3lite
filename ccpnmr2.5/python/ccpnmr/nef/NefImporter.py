"""
NefImporter - a series of routines for reading a Nef file and examining the contents.

Module Contents
===============

Introduction
Error handling
Examples
Nef File Contents


Introduction
------------

NefImporter consists of two classes: NefImporter - a class for handling the top-level object, and
NefDict for handling individual saveFrames in the dictionary.

NefImporter contains:

  initialise          initialise a new dictionary
  loadFile            read in the contents of a .nef file
  saveFile            save the dictionary to a .nef file

  getCategories       return the current categories defined in the Nef structure
  getSaveFrameNames   return the names of the saveFrames with the file
  hasSaveFrame        return True if the saveFrame exists
  getSaveFrame        return saveFrame of the given name
  addSaveFrame        add a new saveFrame to the dictionary


  get<name>           return the relevant structures of the Nef file
                      defined by the available categories, where <name> can be:

                          NmrMetaData
                          MolecularSystems
                          ChemicalShiftLists
                          DistanceRestraintLists
                          DihedralRestraintLists
                          RdcRestraintLists
                          NmrSpectra
                          PeakRestraintLinks

                      e.g. yourImport.getChemicalShiftLists()

  add<name>           add a new saveFrame to the dictionary
                      defined by the available categories, where <name> can be:

                          ChemicalShiftList
                          DistanceRestraintList
                          DihedralRestraintList
                          RdcRestraintList
                          NmrSpectra
                          PeakLists
                          LinkageTables

                      e.g. addChemicalShiftList

  toString            convert Nef dictionary to a string that can be written to a file
  fromString          convert string to Nef dictionary

  getAttributeNames   get a list of the attributes attached to the dictionary
  getAttribute        return the value of the attribute
  hasAttribute        return True if the attribute Exists

  lastError           error code of the last operation
  lastErrorString     error string of the last operation


NefDict contains handling routines:

  getTableNames   return a list of the tables in the saveFrame
  getTable        return table from the saveFrame, it can be returned as an OrderedDict
                  or as a Pandas DataFrame
  hasTable        return true of the table exists
  setTable        set the table - currently not implemented

  getAttributeNames   get a list of the attributes attached to the saveFrame
  getAttribute        return the value of the attribute
  hasAttribute        return True if the attribute Exists

  lastError           error code of the last operation
  lastErrorString     error string of the last operation


Error Handling
--------------

Errors can be handled in three different modes:

  'silent'            errors are handled internally and can be interrogated with saveFrame.lastError
                      with no logging to the stderr

  'standard'          errors are handled internally, error messages are logged to stderr.

  'strict'            errors message are logged to stderr and errors are raised to be trapped by
                      the calling functions

  error handling mode can be set at the instantiation of the object, e.g.

    newObject = NefImporter(errorLogging='standard')


Examples
--------
Here are a few examples of using the classes:

  # load a Nef file
  test = NefImporter(errorLogging=NEF_STANDARD)
  test.loadFile('/Users/account/Projects/NefFile.nef')

  # get categories
  print (test.getCategories())

  # get saveFrame names
  names = test.getSaveFrameNames(); print(names)
  names = test.getSaveFrameNames(returnType=NEF_RETURNALL); print(names)
  names = test.getSaveFrameNames(returnType=NEF_RETURNNEF); print (names)
  names = test.getSaveFrameNames(returnType=NEF_RETURNOTHER); print (names)

  # get a particular saveFrame
  sf1 = test.getSaveFrame(names[0])

  # convert NefImporter into a string for saving
  ts = test.toString()
  test.fromString(ts)

  # getting tables from a saveFrame
  print (sf1.getTableNames())
  table = sf1.getTable('nmr_atom', asPandas=True)
  print (table)
  print (sf1.hasTable('nmr_residue'))
  print (sf1.getAttributeNames())
  print (sf1.hasAttribute('sf_framecode'))
  print (sf1.hasAttribute('nothing'))
  print (test.getSaveFrame(name='ccpn_assignment').getTable(name='nmr_residue', asPandas=True))

  # saving a file
  print ('SAVE ', test.saveFile('/Users/ejb66/PycharmProjects/Sec5Part3testing.nef'))
  print (test.lastError)

  # test meta creation of category names
  print (test.getMolecularSystems())

There are more examples i the __main__ function at the bottom of the module


Nef File Contents
-----------------

More details of the contents of Nef files can be found in GenericStarParser
The general structure of a Nef file is:

::

    DataExtent
      DataBlock
        Item
        Loop
        SaveFrame
          Item
          Loop

DataExtent, DataBlock and SaveFrame are Python OrderedDict with an additional 'name' attribute
DataBlocks and SaveFrames are entered in their container using their name as the key.

Loop is an object with a 'columns' list, a 'data' list-of-row-OrderedDict, and a name attribute
set equal to the name of the first column. A loop is entered in its container under each
column name, so that e.g. aSaveFrame['_Loopx.loopcol1'] and aSaveFrame['_Loopx.loopcol2'] both
exist and both correspond to the same loop object.

Items are entered as a string key - string value pair.
                                    the string value can be a dictionary
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2019"
__credits__ = ("Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2017-07-07 16:32:41 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.0 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: Ed Brooksbank $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================


import os
import sys
import numpy as np
from collections import OrderedDict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this is a fix to get the import to work when running as a standalone
# when importing into your own code, it can be safely removed

def import_parents(level=1):
    global __package__
    file = Path(__file__).resolve()
    parent, top = file.parent, file.parents[level]

    sys.path.append(str(top))
    try:
        sys.path.remove(str(parent))
    except ValueError:  # already removed
        pass

    __package__ = '.'.join(parent.parts[len(top.parts):])
    importlib.import_module(__package__)  # won't be needed after that


if __name__ == '__main__' and __package__ is None:
    import importlib
    from pathlib import Path

    import_parents()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


from . import StarIo
from . import ErrorLog as el
from . import Validator
from . import Specification


MAJOR_VERSION = '0'
MINOR_VERSION = '8'
PATCH_LEVEL = '4'
__nef_version__ = '.'.join((MAJOR_VERSION, MINOR_VERSION))
# __version__ = '.'.join( (__nef_version__, PATCH_LEVEL) )


NEF_CATEGORIES = [('nef_nmr_meta_data', 'get_nmr_meta_data'),
                  ('nef_molecular_system', 'get_molecular_systems'),
                  ('nef_chemical_shift_list', 'get_chemical_shift_lists'),
                  ('nef_distance_restraint_list', 'get_distance_restraint_lists'),
                  ('nef_dihedral_restraint_list', 'get_dihedral_restraint_lists'),
                  ('nef_rdc_restraint_list', 'get_rdc_restraint_lists'),
                  ('nef_nmr_spectrum', 'get_nmr_spectra'),
                  ('nef_peak_restraint_links', 'get_peak_restraint_links')]

NEF_REQUIRED_SAVEFRAME_BY_FRAMECODE = ['nef_nmr_meta_data',
                                       'nef_molecular_system']
NEF_REQUIRED_SAVEFRAME_BY_CATEGORY = ['nef_chemical_shift_list', ]

NEF_ALL_SAVEFRAME_REQUIRED_FIELDS = ['sf_category',
                                     'sf_framecode', ]

MD_REQUIRED_FIELDS = ['sf_category',
                      'sf_framecode',
                      'format_name',
                      'format_version',
                      'program_name',
                      'program_version',
                      'creation_date',
                      'uuid']
MD_OPTIONAL_FIELDS = ['coordinate_file_name', ]
MD_OPTIONAL_LOOPS = ['nef_related_entries',
                     'nef_program_script',
                     'nef_run_history']
MD_RE_REQUIRED_FIELDS = ['database_name',
                         'database_accession_code']
MD_PS_REQUIRED_FIELDS = ['program_name', ]
MD_RH_REQUIRED_FIELDS = ['run_ordinal',
                         'program_name']
MD_RH_OPTIONAL_FIELDS = ['program_version',
                         'script_name',
                         'script']

MS_REQUIRED_FIELDS = ['sf_category',
                      'sf_framecode']
MS_REQUIRED_LOOPS = ['nef_sequence']
MS_OPTIONAL_LOOPS = ['nef_covalent_links']
MS_NS_REQUIRED_FIELDS = ['chain_code',
                         'sequence_code',
                         'residue_type',
                         'linking',
                         'residue_variant']
MS_CL_REQUIRED_FIELDS = ['chain_code_1',
                         'sequence_code_1',
                         'residue_type_1',
                         'atom_name_1',
                         'chain_code_2',
                         'sequence_code_2',
                         'residue_type_2',
                         'atom_name_2']

CSL_REQUIRED_FIELDS = ['sf_category',
                       'sf_framecode',
                       'atom_chem_shift_units']
CSL_REQUIRED_LOOPS = ['nef_chemical_shift']
CSL_CS_REQUIRED_FIELDS = ['chain_code',
                          'sequence_code',
                          'residue_type',
                          'atom_name',
                          'value']
CSL_CS_OPTIONAL_FIELDS = ['value_uncertainty', ]

DRL_REQUIRED_FIELDS = ['sf_category',
                       'sf_framecode',
                       'potential_type']
DRL_REQUIRED_LOOPS = ['nef_distance_restraint']
DRL_OPTIONAL_FIELDS = ['restraint_origin', ]
DRL_DR_REQUIRED_FIELDS = ['ordinal',
                          'restraint_id',
                          'chain_code_1',
                          'sequence_code_1',
                          'residue_type_1',
                          'atom_name_1',
                          'chain_code_2',
                          'sequence_code_2',
                          'residue_type_2',
                          'atom_name_2',
                          'weight']
DRL_DR_OPTIONAL_FIELDS = ['restraint_combination_id',
                          'target_value',
                          'target_value_uncertainty',
                          'lower_linear_limit',
                          'lower_limit',
                          'upper_limit',
                          'upper_linear_limit']

DIHRL_REQUIRED_FIELDS = ['sf_category',
                         'sf_framecode',
                         'potential_type']
DIHRL_REQUIRED_LOOPS = ['nef_dihedral_restraint']
DIHRL_OPTIONAL_FIELDS = ['restraint_origin', ]
DIHRL_DIHR_REQUIRED_FIELDS = ['ordinal',
                              'restraint_id',
                              'restraint_combination_id',
                              'chain_code_1',
                              'sequence_code_1',
                              'residue_type_1',
                              'atom_name_1',
                              'chain_code_2',
                              'sequence_code_2',
                              'residue_type_2',
                              'atom_name_2',
                              'chain_code_3',
                              'sequence_code_3',
                              'residue_type_3',
                              'atom_name_3',
                              'chain_code_4',
                              'sequence_code_4',
                              'residue_type_4',
                              'atom_name_4',
                              'weight']
DIHRL_DIHR_OPTIONAL_FIELDS = ['target_value',
                              'target_value_uncertainty',
                              'lower_linear_limit',
                              'lower_limit',
                              'upper_limit',
                              'upper_linear_limit',
                              'name']

RRL_REQUIRED_FIELDS = ['sf_category',
                       'sf_framecode',
                       'potential_type']
RRL_REQUIRED_LOOPS = ['nef_rdc_restraint']
RRL_OPTIONAL_FIELDS = ['restraint_origin',
                       'tensor_magnitude',
                       'tensor_rhombicity',
                       'tensor_chain_code',
                       'tensor_sequence_code',
                       'tensor_residue_type', ]
RRL_RR_REQUIRED_FIELDS = ['ordinal',
                          'restraint_id',
                          'chain_code_1',
                          'sequence_code_1',
                          'residue_type_1',
                          'atom_name_1',
                          'chain_code_2',
                          'sequence_code_2',
                          'residue_type_2',
                          'atom_name_2',
                          'weight']
RRL_RR_OPTIONAL_FIELDS = ['restraint_combination_id',
                          'target_value',
                          'target_value_uncertainty',
                          'lower_linear_limit',
                          'lower_limit',
                          'upper_limit',
                          'upper_linear_limit',
                          'scale',
                          'distance_dependent', ]

PL_REQUIRED_FIELDS = ['sf_category',
                      'sf_framecode',
                      'num_dimensions',
                      'chemical_shift_list']
PL_REQUIRED_LOOPS = ['nef_spectrum_dimension',
                     'nef_spectrum_dimension_transfer',
                     'nef_peak']
PL_OPTIONAL_FIELDS = ['experiment_classification',
                      'experiment_type']
PL_SD_REQUIRED_FIELDS = ['dimension_id',
                         'axis_unit',
                         'axis_code']
PL_SD_OPTIONAL_FIELDS = ['spectrometer_frequency',
                         'spectral_width',
                         'value_first_point',
                         'folding',
                         'absolute_peak_positions',
                         'is_acquisition', ]
PL_SDT_REQUIRED_FIELDS = ['dimension_1',
                          'dimension_2',
                          'transfer_type']
PL_SDT_OPTIONAL_FIELDS = ['is_indirect', ]
PL_P_REQUIRED_FIELDS = ['ordinal',
                        'peak_id']
PL_P_REQUIRED_ALTERNATE_FIELDS = [['height', 'volume'], ]
PL_P_REQUIRED_FIELDS_PATTERN = ['position_{}',
                                'chain_code_{}',
                                'sequence_code_{}',
                                'residue_type_{}',
                                'atom_name_{}', ]
PL_P_OPTIONAL_ALTERNATE_FIELDS = {r'(height)'         : ['{}_uncertainty', ],
                                  r'(volume)'         : ['{}_uncertainty', ],
                                  r'position_([0-9]+)': ['position_uncertainty_{}', ],
                                  }
PL_P_OPTIONAL_FIELDS_PATTERN = ['position_uncertainty_{}', ]

PRLS_REQUIRED_FIELDS = ['sf_category',
                        'sf_framecode']
PRLS_REQUIRED_LOOPS = ['nef_peak_restraint_link']
PRLS_PRL_REQUIRED_FIELDS = ['nmr_spectrum_id',
                            'peak_id',
                            'restraint_list_id',
                            'restraint_id']

NEF_CATEGORIES_REMOVEPREFIX = {'nef_distance_restraint'         : 'distance_restraint',
                               'nef_molecular_system'           : 'molecular_system',
                               'nef_covalent_links'             : 'covalent_links',
                               'nef_peak_restraint_links'       : 'peak_restraint_links',
                               'nef_run_history'                : 'run_history',
                               'nef_nmr_meta_data'              : 'nmr_meta_data',
                               'nef_rdc_restraint_list'         : 'rdc_restraint_list',
                               'nef_peak_restraint_link'        : 'peak_restraint_link',
                               'nef_nmr_spectrum'               : 'nmr_spectrum',
                               'nef_spectrum_dimension'         : 'spectrum_dimension',
                               'nef_chemical_shift_list'        : 'chemical_shift_list',
                               'nef_sequence'                   : 'sequence',
                               'nef_program_script'             : 'program_script',
                               'nef_related_entries'            : 'related_entries',
                               'nef_distance_restraint_list'    : 'distance_restraint_list',
                               'nef_rdc_restraint'              : 'rdc_restraint',
                               'nef_chemical_shift'             : 'chemical_shift',
                               'nef_spectrum_dimension_transfer': 'spectrum_dimension_transfer',
                               'nef_dihedral_restraint_list'    : 'dihedral_restraint_list',
                               'nef_peak'                       : 'peak',
                               'nef_dihedral_restraint'         : 'dihedral_restraint'}

NEF_CATEGORIES_INSERTPREFIX = {'distance_restraint'         : 'nef_distance_restraint',
                               'molecular_system'           : 'nef_molecular_system',
                               'covalent_links'             : 'nef_covalent_links',
                               'peak_restraint_links'       : 'nef_peak_restraint_links',
                               'run_history'                : 'nef_run_history',
                               'nmr_meta_data'              : 'nef_nmr_meta_data',
                               'rdc_restraint_list'         : 'nef_rdc_restraint_list',
                               'peak_restraint_link'        : 'nef_peak_restraint_link',
                               'nmr_spectrum'               : 'nef_nmr_spectrum',
                               'spectrum_dimension'         : 'nef_spectrum_dimension',
                               'chemical_shift_list'        : 'nef_chemical_shift_list',
                               'sequence'                   : 'nef_sequence',
                               'program_script'             : 'nef_program_script',
                               'related_entries'            : 'nef_related_entries',
                               'distance_restraint_list'    : 'nef_distance_restraint_list',
                               'rdc_restraint'              : 'nef_rdc_restraint',
                               'chemical_shift'             : 'nef_chemical_shift',
                               'spectrum_dimension_transfer': 'nef_spectrum_dimension_transfer',
                               'dihedral_restraint_list'    : 'nef_dihedral_restraint_list',
                               'peak'                       : 'nef_peak',
                               'dihedral_restraint'         : 'nef_dihedral_restraint'}

NEF_RETURNALL = 'all'
NEF_RETURNNEF = 'nef_'
NEF_RETURNOTHER = 'other'
NEF_PREFIX = 'nef_'


class NefImporter(el.ErrorLog):
    """Top level data block for accessing object tree"""

    # put functions in here to read the contents of the dict.
    # superclassed from DataBlock which is of type StarContainer
    def __init__(self, name=None,
                 programName='Unknown',
                 programVersion='Unknown',
                 initialise=True,
                 errorLogging=el.NEF_STANDARD,
                 hidePrefix=True):

        el.ErrorLog.__init__(self, loggingMode=errorLogging)

        self.name = name
        self._nefDict = StarIo.NmrDataBlock()
        self._validateNefDict = StarIo.NmrDataBlock()
        self.programName = programName
        self.programVersion = programVersion
        self._hidePrefix = hidePrefix
        self._saveFrameNames = {}
        self._valid = Validator.Validator(self._nefDict)

        if initialise:
            # initialise a basic object
            self.initialise()

    @el.ErrorLog(errorCode=el.NEFERROR_ERRORLOADINGFILE)
    def loadValidateDictionary(self, fileName=None, mode='standard'):
        if not os.path.isfile(fileName):
            raise ValueError('Error: file does not exist')

        with open(fileName) as fp:
            data = fp.read()
        converter = Specification.CifDicConverter(data)
        converter.convertToNef()
        self._validateNefDict = converter.result

        return True

    @property
    def isValid(self):
        """
        Check whether the Nef object contains the required information
        :return True or False:
        """
        return self._valid.isValid(self._nefDict, self._validateNefDict)

    @property
    def validErrorLog(self):
        """
        Return the error log from checking validity
        :return dict:
        """
        return self._valid._validation_errors

    def _namedToNefDict(self, frame):
        # change a saveFrame into a normal OrderedDict
        newItem = NefDict(inFrame=frame, errorLogging=self.loggingMode)
        for ky in frame.keys():
            newItem[ky] = frame[ky]
        return newItem

    def _removePrefix(self, name):
        if self._hidePrefix:
            for db in NEF_CATEGORIES_REMOVEPREFIX.keys():
                if name.startswith(db):
                    name = name.replace(db, NEF_CATEGORIES_REMOVEPREFIX[db], 1)
                    break
        return name

    def _insertPrefix(self, name):
        if self._hidePrefix:
            for db in NEF_CATEGORIES_INSERTPREFIX.keys():
                if name.startswith(db):
                    name = name.replace(db, NEF_CATEGORIES_INSERTPREFIX[db], 1)
                    break
        return name

    @el.ErrorLog(errorCode=el.NEFERROR_BADLISTTYPE)
    def _getListType(self, _listType):
        # return a list of '_listType' from the saveFrame,
        # used with nefCategory routines below
        if self._nefDict and isinstance(self._nefDict, OrderedDict):
            sfList = [self._nefDict[db] for db in self._nefDict.keys() if _listType in db]
            sfList = [self._namedToNefDict(sf) for sf in sfList]

            # if there is only one item then return it, otherwise return the list
            if len(sfList) > 1:
                return sfList
            elif sfList:
                return sfList[0]

        return None

    # routines to get the Nef specific data from the dictionary
    def getNmrMetaData(self):
        """
        Return the nef_nmr_meta_data saveFrames
        :return list or single item:
        """
        return self._getListType(NEF_CATEGORIES[0][0])

    def getMolecularSystems(self):
        """
        Return the nef_molecular_system saveFrames
        :return list or single item:
        """
        return self._getListType(NEF_CATEGORIES[1][0])

    def getChemicalShiftLists(self):
        """
        Return the nef_chemical_shift_list saveFrames
        :return list or single item:
        """
        return self._getListType(NEF_CATEGORIES[2][0])

    def getDistanceRestraintLists(self):
        """
        Return the nef_distance_restraint_list saveFrames
        :return list or single item:
        """
        return self._getListType(NEF_CATEGORIES[3][0])

    def getDihedralRestraintLists(self):
        """
        Return the nef_dihedral_restraint_list saveFrames
        :return list or single item:
        """
        return self._getListType(NEF_CATEGORIES[4][0])

    def getRdcRestraintLists(self):
        """
        Return the nef_rdc_restraint_list saveFrames
        :return list or single item:
        """
        return self._getListType(NEF_CATEGORIES[5][0])

    def getNmrSpectra(self):
        """
        Return the nef_nmr_spectrum saveFrames
        :return list or single item:
        """
        return self._getListType(NEF_CATEGORIES[6][0])

    def getPeakRestraintLinks(self):
        """
        Return the nef_peak_restraint_link saveFrames
        :return list or single item:
        """
        return self._getListType(NEF_CATEGORIES[7][0])

    def initialise(self):
        """
        Initialise a new NefImporter object with a starting saveFrame
        """
        nefNmr = 'nef_nmr_meta_data'
        nefMol = 'nef_molecular_system'
        nefChem = 'nef_chemical_shift_list_1'

        self._nefDict[nefNmr] = StarIo.NmrDataBlock()
        self._nefDict[nefNmr].update({k: '' for k in MD_REQUIRED_FIELDS})
        self._nefDict[nefNmr]['sf_category'] = 'nef_nmr_meta_data'
        self._nefDict[nefNmr]['sf_framecode'] = 'nef_nmr_meta_data'
        self._nefDict[nefNmr]['format_name'] = 'Nmr_Exchange_Format'
        self._nefDict[nefNmr]['format_version'] = __nef_version__
        self._nefDict[nefNmr]['program_name'] = self.programName
        self._nefDict[nefNmr]['program_version'] = self.programVersion

        self._nefDict[nefMol] = StarIo.NmrDataBlock()
        self._nefDict[nefMol].update({k: '' for k in MS_REQUIRED_FIELDS})
        self._nefDict[nefMol]['sf_category'] = 'nef_molecular_system'
        self._nefDict[nefMol]['sf_framecode'] = 'nef_molecular_system'
        for l in MS_REQUIRED_LOOPS:
            self._nefDict['nef_molecular_system'][l] = []
            self.addChemicalShiftList(nefChem, 'ppm')

        self._logError(errorCode=el.NEFVALID)

    @el.ErrorLog(errorCode=el.NEFERROR_BADADDSAVEFRAME)
    def addSaveFrame(self, name, category, required_fields=None, required_loops=None):
        """
        Add a new saveFrame to NefImporter
        :param name:
        :param category:
        :param required_fields:
        :param required_loops:
        """
        name = self._insertPrefix(name)

        self._nefDict[name] = StarIo.NmrSaveFrame()
        if required_fields is not None:
            self._nefDict[name].update({k: '' for k in required_fields})
            self._nefDict[name]['sf_category'] = category
            self._nefDict[name]['sf_framecode'] = name
        if required_loops is not None:
            for l in required_loops:
                self._nefDict[name][l] = []
        return self._nefDict[name]

    @el.ErrorLog(errorCode=el.NEFERROR_BADADDSAVEFRAME)
    def addChemicalShiftList(self, name, cs_units='ppm'):
        name = self._insertPrefix(name)

        category = 'nef_chemical_shift_list'
        self.addSaveFrame(name=name, category=category,
                          required_fields=CSL_REQUIRED_FIELDS,
                          required_loops=CSL_REQUIRED_LOOPS)
        self._nefDict[name]['atom_chem_shift_units'] = cs_units
        return self._nefDict[name]

    @el.ErrorLog(errorCode=el.NEFERROR_BADADDSAVEFRAME)
    def addDistanceRestraintList(self, name, potential_type,
                                 restraint_origin=None):
        name = self._insertPrefix(name)

        category = 'nef_distance_restraint_list'
        self.addSaveFrame(name=name, category=category,
                          required_fields=DRL_REQUIRED_FIELDS,
                          required_loops=DRL_REQUIRED_LOOPS)
        self._nefDict[name]['potential_type'] = potential_type
        if restraint_origin is not None:
            self._nefDict[name]['restraint_origin'] = restraint_origin

        return self._nefDict[name]

    @el.ErrorLog(errorCode=el.NEFERROR_BADADDSAVEFRAME)
    def addDihedralRestraintList(self, name, potential_type,
                                 restraint_origin=None):
        name = self._insertPrefix(name)

        category = 'nef_dihedral_restraint_list'
        self.addSaveFrame(name=name, category=category,
                          required_fields=DIHRL_REQUIRED_FIELDS,
                          required_loops=DIHRL_REQUIRED_LOOPS)
        self._nefDict[name]['potential_type'] = potential_type
        if restraint_origin is not None:
            self._nefDict[name]['restraint_origin'] = restraint_origin

        return self._nefDict[name]

    @el.ErrorLog(errorCode=el.NEFERROR_BADADDSAVEFRAME)
    def addRdcRestraintList(self, name, potential_type,
                            restraint_origin=None, tensor_magnitude=None,
                            tensor_rhombicity=None, tensor_chain_code=None,
                            tensor_sequence_code=None, tensor_residue_type=None):
        name = self._insertPrefix(name)

        category = 'nef_rdc_restraint_list'
        self.addSaveFrame(name=name, category=category,
                          required_fields=DIHRL_REQUIRED_FIELDS,
                          required_loops=RRL_REQUIRED_LOOPS)
        self._nefDict[name]['potential_type'] = potential_type
        if restraint_origin is not None:
            self._nefDict[name]['restraint_origin'] = restraint_origin
        if tensor_magnitude is not None:
            self._nefDict[name]['tensor_magnitude'] = tensor_magnitude
        if tensor_rhombicity is not None:
            self._nefDict[name]['tensor_rhombicity'] = tensor_rhombicity
        if tensor_chain_code is not None:
            self._nefDict[name]['tensor_chain_code'] = tensor_chain_code
        if tensor_sequence_code is not None:
            self._nefDict[name]['tensor_sequence_code'] = tensor_sequence_code
        if tensor_residue_type is not None:
            self._nefDict[name]['tensor_residue_type'] = tensor_residue_type

        return self._nefDict[name]

    @el.ErrorLog(errorCode=el.NEFERROR_BADADDSAVEFRAME)
    def addPeakList(self, name, num_dimensions, chemical_shift_list,
                    experiment_classification=None,
                    experiment_type=None):
        name = self._insertPrefix(name)

        category = 'nef_nmr_spectrum'
        if chemical_shift_list in self:
            if self._nefDict[chemical_shift_list]['sf_category'] == 'nef_chemical_shift_list':
                self.addSaveFrame(name=name, category=category,
                                  required_fields=PL_REQUIRED_FIELDS,
                                  required_loops=PL_REQUIRED_LOOPS)
                self._nefDict[name]['num_dimensions'] = num_dimensions
                self._nefDict[name]['chemical_shift_list'] = chemical_shift_list
                if experiment_classification is not None:
                    self._nefDict[name]['experiment_classification'] = experiment_classification
                if experiment_type is not None:
                    self._nefDict[name]['experiment_type'] = experiment_type

                return self._nefDict[name]
            raise Exception('{} is not a nef_chemical_shift_list.'.format(chemical_shift_list))
        raise Exception('{} does not exist.'.format(chemical_shift_list))

    @el.ErrorLog(errorCode=el.NEFERROR_BADADDSAVEFRAME)
    def addLinkageTable(self):
        name = category = 'nef_peak_restraint_links'
        return self.addSaveFrame(name=name, category=category,
                                 required_fields=PRLS_REQUIRED_FIELDS,
                                 required_loops=PL_REQUIRED_LOOPS)

    @el.ErrorLog(errorCode=el.NEFERROR_BADTOSTRING)
    def toString(self):
        return self._nefDict.toString()

    @el.ErrorLog(errorCode=el.NEFERROR_BADFROMSTRING)
    def fromString(self, text, mode='standard'):
        # set the Nef from the contents of the string, opposite of toString
        dataExtent = StarIo.parseNef(text=text, mode=mode)
        if dataExtent:
            dbs = [dataExtent[db] for db in dataExtent.keys()]
            if dbs:
                self._nefDict = dbs[0]
        else:
            self._logError(errorCode=el.NEFERROR_BADFROMSTRING)
            self._nefDict = None

    @el.ErrorLog(errorCode=el.NEFERROR_ERRORLOADINGFILE)
    def loadFile(self, fileName=None, mode='standard'):
        nefDataExtent = StarIo.parseNefFile(fileName=fileName, mode=mode)
        self._nefDict = list(nefDataExtent.values())
        if len(self._nefDict) > 1:
            print('More than one datablock in a NEF file is not allowed.  Using the first and discarding the rest.')
        self._nefDict = self._nefDict[0]

        return True

    @el.ErrorLog(errorCode=el.NEFERROR_ERRORSAVINGFILE)
    def saveFile(self, fileName=None):
        with open(fileName, 'w') as op:
            op.write(self._nefDict.toString())

        return True

    @el.ErrorLog(errorCode=el.NEFERROR_BADCATEGORIES)
    def getCategories(self):
        # return a list of the categories available in a Nef file
        return tuple([self._removePrefix(nm[0]) for nm in NEF_CATEGORIES])

    @el.ErrorLog(errorCode=el.NEFERROR_SAVEFRAMEDOESNOTEXIST)
    def getSaveFrameNames(self, returnType=NEF_RETURNALL):
        # return a list of the saveFrames in the file
        names = [db for db in self._nefDict.keys()
                 if isinstance(self._nefDict[db], StarIo.NmrSaveFrame)]

        if returnType == NEF_RETURNNEF:
            names = [self._removePrefix(nm) for nm in names if nm and nm.startswith(NEF_PREFIX)]
        elif returnType == NEF_RETURNOTHER:
            names = [nm for nm in names if nm and not nm.startswith(NEF_PREFIX)]
        else:
            names = [self._removePrefix(nm) for nm in names]

        return tuple(names)

    @el.ErrorLog(errorCode=el.NEFERROR_SAVEFRAMEDOESNOTEXIST)
    def hasSaveFrame(self, name):
        # return True if the saveFrame exists, else False
        name = self._insertPrefix(name)
        return name in self._nefDict

    @el.ErrorLog(errorCode=el.NEFERROR_SAVEFRAMEDOESNOTEXIST)
    def getSaveFrame(self, name):
        # return the saveFrame 'name'
        name = self._insertPrefix(name)
        return NefDict(self._nefDict[name], errorLogging=self.loggingMode, hidePrefix=self._hidePrefix)

    @el.ErrorLog(errorCode=el.NEFERROR_READATTRIBUTENAMES)
    def getAttributeNames(self):
        names = [self._removePrefix(db) for db in self._nefDict.keys()
                 if not isinstance(self._nefDict[db], StarIo.NmrSaveFrame)]
        return tuple(names)

    @el.ErrorLog(errorCode=el.NEFERROR_READATTRIBUTE)
    def getAttribute(self, name):
        name = self._insertPrefix(name)
        return self._nefDict[name]

    @el.ErrorLog(errorCode=el.NEFERROR_READATTRIBUTE)
    def hasAttribute(self, name):
        name = self._insertPrefix(name)
        return name in self._nefDict

    @property
    def hidePrefix(self):
        # return the current hidePrefix state
        # True - Nef prefixes 'nef_' are hidden
        # False - Nef prefixes 'nef_' can be seen
        # prefixes are still used in the saveFrames bit not seen in general use
        return self._hidePrefix

    @hidePrefix.setter
    def hidePrefix(self, newPrefix):
        # set the current hidePrefix state
        # True - Nef prefixes 'nef_' are hidden
        # False - Nef prefixes 'nef_' can be seen
        # prefixes are still used in the saveFrames bit not seen in general use
        if isinstance(newPrefix, bool):
            self._hidePrefix = newPrefix


class NefDict(StarIo.NmrSaveFrame, el.ErrorLog):
    """
    An orderedDict saveFrame object for extracting information from the NefImporter
    """

    def __init__(self, inFrame, errorLogging=el.NEF_STANDARD, hidePrefix=True):
        """
        Initialise a NefDict orderedDict object from a given saveFrame
        :param inFrame:
        :param errorLogging:
        :param hidePrefix:
        """
        StarIo.NmrSaveFrame.__init__(self, name=inFrame.name)
        el.ErrorLog.__init__(self, loggingMode=errorLogging)

        self._nefFrame = inFrame
        self._hidePrefix = hidePrefix

    def _removePrefix(self, name):
        """
        Remove the prefix 'nef_' from the saveFrame category name
        :param nef_name:
        :return <name>:
        """
        if self._hidePrefix:
            for db in NEF_CATEGORIES_REMOVEPREFIX.keys():
                if name.startswith(db):
                    name = name.replace(db, NEF_CATEGORIES_REMOVEPREFIX[db], 1)
                    break
        return name

    def _insertPrefix(self, name):
        """
        Insert the prefix 'nef_' into the saveFrame category name
        :param name:
        :return nef_<name>:
        """
        if self._hidePrefix:
            for db in NEF_CATEGORIES_INSERTPREFIX.keys():
                if name.startswith(db):
                    name = name.replace(db, NEF_CATEGORIES_INSERTPREFIX[db], 1)
                    break
        return name

    @el.ErrorLog(errorCode=el.NEFERROR_BADKEYS)
    def _namedToOrderedDict(self, frame):
        """
        Convert a named saveFrame to a standard orderedDict
        :param frame:
        :return orderedDict:
        """
        newItem = OrderedDict()
        for ky in frame.keys():
            newItem[ky] = frame[ky]
        return newItem

    @el.ErrorLog(errorCode=el.NEFERROR_BADTABLENAMES)
    def getTableNames(self):
        """
        Return list of attributes in the saveFrame
        :return list or None:
        """
        return tuple([self._removePrefix(db) for db in self._nefFrame.keys()
                      if isinstance(self._nefFrame[db], StarIo.NmrLoop)])

    @el.ErrorLog(errorCode=el.NEFERROR_GENERICGETTABLEERROR)
    def getTable(self, name=None, asPandas=False):
        """
        Return the table 'name' from the saveFrame if it exists
        if asPandas is True then return as a pandas dataFrame,
        otherwise return a list of saveFrames
        :param name:
        :param asPandas:
        :return saveFrames, dataFrames or None:
        """
        # return table 'name' if exists else None
        thisFrame = None
        if name:
            if name in self._nefFrame:
                thisFrame = self._nefFrame[name]
            else:
                # table not found
                name = self._insertPrefix(name)
                if name in self._nefFrame:
                    thisFrame = self._nefFrame[name]
                else:
                    self._logError(errorCode=el.NEFERROR_TABLEDOESNOTEXIST)
                    return None
        else:
            tables = self.getTableNames()
            if tables:
                thisFrame = self._nefFrame[tables[0]]

        if asPandas:
            return self._convertToPandas(thisFrame)
        else:
            return [self._namedToOrderedDict(sf) for sf in thisFrame.data]

    @el.ErrorLog(errorCode=el.NEFERROR_BADMULTICOLUMNVALUES)
    def multiColumnValues(self, column=None):
        """
        Return tuple of orderedDict of values for columns.
        Will work whether columns are in a loop or single values
        If columns match a single loop or nothing, return the loop data.
        Otherwise return a tuple with a single OrderedDict.
        If no column matches return None
        If columns match more than one loop throw an error
        :param column:
        :return orderedDicts:
        """
        return self._nefFrame.multiColumnValues(column=column)

    @el.ErrorLog(errorCode=el.NEFERROR_TABLEDOESNOTEXIST)
    def hasTable(self, name):
        """
        Return True if table 'name' exists in the saveFrame
        :param name:
        :return True or False:
        """
        # return True is the table exists in the saveFrame
        name = self._insertPrefix(name)
        return name in self._nefFrame

    def setTable(self, name):
        # add the table 'name' to the saveFrame, or replace the existing
        # does this need to be here or in the main class?
        print('Not implemented yet.')

    def _convertToPandas(self, sf):
        """
        Convert the saveFrame to a Pandas dataFrame
        :return dataFrame or None on error:
        """
        try:
            import pandas as pd

            df = pd.DataFrame(data=sf.data, columns=sf.columns)
            df.replace({'.': np.NAN, 'true': True, 'false': False}, inplace=True)
            return df
        except:
            return None

    @el.ErrorLog(errorCode=el.NEFERROR_READATTRIBUTENAMES)
    def getAttributeNames(self):
        """
        Return list of attributes in the saveFrame
        :return list or None:
        """
        return tuple([self._removePrefix(db) for db in self._nefFrame.keys()
                      if not isinstance(self._nefFrame[db], StarIo.NmrLoop)])

    @el.ErrorLog(errorCode=el.NEFERROR_READATTRIBUTE)
    def getAttribute(self, name):
        """
        Return attribute 'name' if in the saveFrame
        :param name:
        :return attribute:
        """
        name = self._insertPrefix(name)
        return self._nefFrame[name]

    @el.ErrorLog(errorCode=el.NEFERROR_READATTRIBUTE)
    def hasAttribute(self, name):
        """
        Return True if attribute 'name' is in the saveFrame
        :param name:
        :return True or False:
        """
        name = self._insertPrefix(name)
        return name in self._nefFrame

    @property
    def hidePrefix(self):
        # return the current hidePrefix state
        # True - Nef prefixes 'nef_' are hidden
        # False - Nef prefixes 'nef_' can be seen
        # prefixes are still used in the saveFrames bit not seen in general use
        return self._hidePrefix

    @hidePrefix.setter
    def hidePrefix(self, newPrefix):
        # set the current hidePrefix state
        # True - Nef prefixes 'nef_' are hidden
        # False - Nef prefixes 'nef_' can be seen
        # prefixes are still used in the saveFrames bit not seen in general use
        if isinstance(newPrefix, bool):
            self._hidePrefix = newPrefix


if __name__ == '__main__':

    # define the NefImporter with standard logging
    test = NefImporter(errorLogging=el.NEF_STANDARD)
    testPath = os.path.join(os.getcwd(), 'NEF', 'data_1_1', 'CCPN_Commented_Example.nef')
    testPathOut = os.path.join(os.getcwd(), 'NEF', 'data_1_1', 'CCPN_Commented_Example_Out.nef')

    if not os.path.exists(testPath):
        raise RuntimeError('Error: %s not found' % testPath)

    # load a testFile
    test.loadFile(testPath)

    # print the categories contained in the nef
    print(test.getCategories())
    names = test.getSaveFrameNames()

    validNef = NefImporter(errorLogging=el.NEF_STANDARD)

    # oad the validation dictionary
    infile = 'mmcif_nef.dic'
    filePath = os.path.join(os.getcwd(), 'NEF', 'specification', infile)
    if not os.path.exists(filePath):
        raise RuntimeError('Error: %s not found' % infile)

    test.loadValidateDictionary(filePath)
    validCheck = test.isValid

    # print out any validation errors
    print('~~~~~~~~~~~~~~~~~~~~~~\nValid Nef:', validCheck)
    if not validCheck:
        print('Error Nef:', test.validErrorLog)
        print('Error Nef:')
        for k, v in test.validErrorLog.items():
            print('  >>>', k)
            for err in v:
                print('  >>>        ', err)
    print('~~~~~~~~~~~~~~~~~~~~~~')

    # examples of printing saveFrames information
    print(names)
    names = test.getSaveFrameNames(returnType=NEF_RETURNALL)
    print(names)
    names = test.getSaveFrameNames(returnType=NEF_RETURNNEF)
    print(names)
    names = test.getSaveFrameNames(returnType=NEF_RETURNOTHER)
    print(names)

    # testing the different types of errors trapping - try to find a saveFrame called 'notFound'
    sf1 = None
    if names:
        sf1 = test.getSaveFrame(names[0])
    sf2 = test.getSaveFrame('notFound')

    print(test.hasSaveFrame('notFound'))
    if names:
        print(test.hasSaveFrame(names[0]))

    ts = test.toString()
    test.fromString(ts)

    # testing the different errors from searching for saveFrames
    if sf1 is not None:

        # get the list of tables
        print(sf1.name)
        print(sf1.getTableNames())

        # get a table as a pandas dataframe
        table = sf1.getTable()
        sf1.getTable('nmr_atom', asPandas=True)
        table = sf1.getTable('nmr_atom', asPandas=True)

        print(table)

        print(sf1.hasTable('notFound'))
        print(sf1.hasTable('nmr_residue'))

        print(sf1.getAttributeNames())
        print(sf1.hasAttribute('sf_framecode'))
        print(sf1.hasAttribute('nothing'))
        print(sf1.getAttribute('notHere'))
        print(sf1.getAttribute('sf_category'))

    print('Testing getTable')
    try:
        print(test.getSaveFrame(name='ccpn_assignment').getTable())
    except Exception as es:
        print('Error: %s' % str(es))

    try:
        print(test.getSaveFrame(name='ccpn_assignment').getTable(name='nmr_residue', asPandas=True))
    except Exception as es:
        print('Error: %s' % str(es))

    try:
        print(test.getSaveFrame(name='ccpn_assignment').getTable(name='notFound', asPandas=True))
    except Exception as es:
        print('Error: %s' % str(es))

    print('Testing saveFile')
    print('SAVE ', test.saveFile(testPathOut))
    print(test.lastError)

    # test meta creation of category names
    print(test.getNmrMetaData())
    print(test.getMolecularSystems())
    print(test.getDihedralRestraintLists())
    print(test.getDistanceRestraintLists())
    print(test.getChemicalShiftLists())
    print(test.getNmrSpectra())
    print(test.getPeakRestraintLinks())
    print(test.getRdcRestraintLists())

    if sf1:
        sf1.loggingMode = el.NEF_STRICT
        sf1.loggingMode = el.NEF_STANDARD
        sf1.loggingMode = el.NEF_SILENT
        print(sf1.loggingMode)
        print(sf1.hidePrefix)

    if test:
        test.loggingMode = el.NEF_STANDARD
        test.loggingMode = el.NEF_SILENT
        test.loggingMode = el.NEF_STRICT
        print(test.getAttributeNames())
        print(test.loggingMode)
        print(test.hidePrefix)

    print('Valid Nef:', test.isValid)
    print('Valid Nef:', test.validErrorLog)

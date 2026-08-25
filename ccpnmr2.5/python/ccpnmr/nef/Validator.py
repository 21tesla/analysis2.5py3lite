"""
Module Documentation here
"""

from __future__ import unicode_literals, print_function, absolute_import, division


#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2020"
__credits__ = ("Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2020-07-06 14:28:15 +0100 (Mon, July 06, 2020) $"
__version__ = "$Revision: 3.0.1 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: TJ Ragan $"
__date__ = "$Date: 2017-03-23 16:50:22 +0000 (Thu, March 23, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import re

NMR_EXCHANGE_FORMAT = 'nmr_exchange_format'
FRAME_PREFIX = 'nef_saveframe_'
FRAME_SEARCH = r'{}(\w+)'.format(FRAME_PREFIX)
SF_CATEGORY = 'sf_category'
SF_FRAMECODE = 'sf_framecode'
NAME = 'name'
NEF_ITEM = 'nef_item'
IS_MANDATORY = 'is_mandatory'
IS_KEY = 'is_key'
LOOP_CATEGORY = 'loop_category'
CATEGORY = 'category'
SAVEFRAME = 'saveframe'
NEF_LOOP = 'nef_loop'
METADATA_KEY = 'METADATA'
METADATA_DICTKEY = 'nef_nmr_meta_data'
SPECIFICATION_KEY = 'nef_specification'
FORMAT_NAME = 'format_name'
FORMAT_VERSION = 'format_version'
CREATION_DATE = 'creation_date'
UUID = 'uuid'
VERSION = 'version'
CCPN_PREFIX = 'ccpn_'


class Validator(object):

    def __init__(self, nef=None, validateNefDict=None):
        self.nef = nef
        self.validateNefDict = validateNefDict
        self._validation_errors = None

    def isValid(self, nef=None, validNef=None):
        """Return whether the Nef file is valid
        """
        if nef is None:
            nef = self.nef
        if validNef is None:
            validNef = self.validateNefDict

        self._validateAll(nef, validNef)

        v = list(self._validation_errors.values())
        return not any(v)

    @property
    def validationErrors(self):
        """Return the dict of validation errors
        """
        if not self._validation_errors:
            self.isValid(self.nef, self.validateNefDict)
        return self._validation_errors

    def _validateAll(self, nef=None, validNef=None):
        """Validate a nef file (nef) against a nef dictionary (validNef)
        """
        if nef is None:
            nef = self.nef
        if validNef is None:
            validNef = self.validateNefDict
        if not nef:
            raise RuntimeError('Error: nef not defined')
        if not validNef:
            raise RuntimeError('Error: validateNefDict not defined')

        self._validation_errors = dict()
        self._validation_errors[SAVEFRAME] = []

        # validate meta_data
        e = self._validation_errors[SAVEFRAME]
        e += self._validate_nmr_meta_data(nef, validNef)

        # go through all the saveframes in the Nef object
        for sf_name, saveframe in nef.items():

            if SF_FRAMECODE not in saveframe or saveframe.name != saveframe[SF_FRAMECODE]:
                e = self._validation_errors[SAVEFRAME]
                e += ["Saveframe.name for sf_framecode '{}' is not defined correctly.".format(saveframe[SF_FRAMECODE]), ]
                break

            # check against the validation dictionary
            for vName, validFrame in validNef.items():

                # get the actual name from the end the the name - may need to be more complex later
                checkName = re.findall(FRAME_SEARCH, validFrame.name)

                if checkName and SF_CATEGORY in saveframe and saveframe[SF_CATEGORY] == checkName[0]:

                    ERROR_KEY = checkName[0]
                    if ERROR_KEY not in self._validation_errors:
                        self._validation_errors[ERROR_KEY] = []
                    e = self._validation_errors[ERROR_KEY]

                    # check items against loop_category = None, i.e., this saveframe
                    mandatoryFields = [nm[NAME] for nm in validFrame[NEF_ITEM].data if nm[IS_MANDATORY] is True and nm[LOOP_CATEGORY] is None]
                    optionalFields = [nm[NAME] for nm in validFrame[NEF_ITEM].data if nm[IS_MANDATORY] is False and nm[LOOP_CATEGORY] is None]
                    loopNames = [nm['category'] for nm in validFrame[NEF_LOOP].data]

                    # check for missing words/framecode is not correct/category is mismatched/bad fields (keys)
                    e += self._sf_framecode_name_mismatch(saveframe, sf_name)
                    e += self._dict_missing_keys(saveframe, mandatoryFields, label=sf_name)
                    e += self._sf_category_name_mismatch(saveframe, checkName[0])
                    e += self._dict_nonallowed_keys(saveframe, mandatoryFields + optionalFields + loopNames, label=sf_name)

                    loops = [kk for kk in saveframe.keys() if kk in loopNames]

                    # check that all the mandatory loops have been included
                    mandatoryLoops = [nm[CATEGORY] for nm in validFrame[NEF_LOOP].data if nm[IS_MANDATORY] is True]
                    e += self._dict_missing_keys(loops, mandatoryLoops, label=sf_name, keyType='loop')

                    # iterate through loops
                    for loop in loops:

                        # get the keys that belong to this loop
                        mandatoryLoopFields = [nm[NAME] for nm in validFrame[NEF_ITEM].data if nm[IS_MANDATORY] is True and nm[LOOP_CATEGORY] == loop]
                        optionalLoopFields = [nm[NAME] for nm in validFrame[NEF_ITEM].data if nm[IS_MANDATORY] is False and nm[LOOP_CATEGORY] == loop]

                        if saveframe[loop]:
                            # NOTE:ED - changed to allow empty loops
                            if saveframe[loop].data:
                                # check for missing words/bad fields (keys)/malformed loops
                                e += self._dict_missing_keys(saveframe[loop].data[0], mandatoryLoopFields, label='{}:{}'.format(sf_name, loop))
                                e += self._dict_nonallowed_keys(saveframe[loop].data[0], mandatoryLoopFields + optionalLoopFields, label='{}:{}'.format(sf_name, loop))
                                e += self._loop_entries_inconsistent_keys(saveframe[loop].data, label='{}:{}'.format(sf_name, loop))
                            else:
                                # there should not be any loops without data - could be mandatory loops
                                # e += ["Loop '{}' contains no data.".format(loop), ]
                                pass

                        else:
                            # this error is a catch-all as loadFile should test the integrity of the nef file before validation
                            e += ["Error reading loop '{}'.".format(loop), ]

                    break

                elif SF_CATEGORY in saveframe and saveframe[SF_CATEGORY].startswith(CCPN_PREFIX):

                    # skip ccpn_ specific categories for now.
                    break

            else:
                e = self._validation_errors[SAVEFRAME]
                e += ["No sf_category '{}' found (possibly bad name defined).".format(saveframe[SF_CATEGORY]),]

        return self._validation_errors

    def _validate_nmr_meta_data(self, nef, validNef):
        """Check that the information in the meta_data is correct for this version
        """
        DICT_KEY = METADATA_DICTKEY
        VALID_KEY = SPECIFICATION_KEY

        if DICT_KEY not in nef:
            return ['No {} saveframe.'.format(DICT_KEY)]
        else:

            # TODO:ED format name should be defined in the mmcif_nef.dic
            format_name = NMR_EXCHANGE_FORMAT

            if VALID_KEY not in validNef and VERSION not in validNef[VALID_KEY]:
                return ['Version not found']
            format_version = float(validNef[VALID_KEY][VERSION])

            md = nef[DICT_KEY]
            e = []

            if FORMAT_NAME in md:
                if md[FORMAT_NAME] != format_name:
                    e.append("{} must be '{}'.".format(format_name, NMR_EXCHANGE_FORMAT))

            if FORMAT_VERSION in md:
                mdVersion = float(md[FORMAT_VERSION])
                if mdVersion < format_version:
                    e.append('This reader (version {}) does not support {}.'.format(format_version, mdVersion))

            if CREATION_DATE in md:
                # TODO: ED How to validate the creation date?
                pass

            if UUID in md:
                # TODO: ED How to validate the uuid?
                pass

            return e

    def _dict_missing_keys(self, dct, required_keys, label=None, keyType='label'):
        if label is None:
            return ['Missing {} {}.'.format(key, keyType) for key in required_keys if key not in dct]
        return ['{}: missing {} {}.'.format(label, key, keyType) for key in required_keys if key not in dct]

    def _dict_duplicate_keys(self, dct, label=None, keyType='label'):
        dctSet = set(dct)
        if len(dct) == len(dctSet):
            return []

        if label is None:
            return ['Duplicated {} {}.'.format(key, keyType) for key in dctSet if dct.count(key) > 1]
        return ['{}: Duplicated {} {}.'.format(label, key, keyType) for key in dctSet if dct.count(key) > 1]

    def _dict_missing_value_with_key(self, dct, keys):
        errors = []
        for key in keys:
            found_key = False
            for k, v in dct.items():
                if ('sf_category' in v) and (v['sf_category'] == key):
                    found_key = True
            if not found_key:
                errors.append('No saveframes with sf_category: {}.'.format(key))
        return errors

    def _sf_framecode_name_mismatch(self, dct, sf_framecode):
        if 'sf_framecode' in dct:
            if dct['sf_framecode'] != sf_framecode:
                return ["sf_framecode {} must match key {}.".format(dct['sf_framecode'], sf_framecode)]
        return []

    def _sf_category_name_mismatch(self, dct, sf_category):
        if 'sf_category' in dct:
            if dct['sf_category'] != sf_category:
                return ["sf_category {} must be {}.".format(dct['sf_category'], sf_category)]
        # else:
        #     return ["No sf_category.",]
        return []

    def _loop_entries_inconsistent_keys(self, loop, label):
        errors = []
        if len(loop) > 0:
            fields = list(loop[0].keys())
            fields_count = len(fields)
            finished = False
            while not finished:
                errors = []
                finished = True
                for i, entry in enumerate(loop):
                    for field in entry:
                        if field not in fields:
                            fields.append(field)
                    if len(fields) > fields_count:
                        fields_count = len(fields)
                        finished = False
                        break
                    errors += self._dict_missing_keys(entry, fields, label=label + ' item {}'
                                                      .format(i))
        return errors

    def _dict_nonallowed_keys(self, dct, allowed_keys, label=None):
        if label is None:
            return ["Field '{}' not allowed.".format(key)
                    for key in dct.keys() if key not in allowed_keys and not key.startswith(CCPN_PREFIX)]
        return ["Field '{}' not allowed in {}.".format(key, label)
                for key in dct.keys() if key not in allowed_keys and not key.startswith(CCPN_PREFIX)]

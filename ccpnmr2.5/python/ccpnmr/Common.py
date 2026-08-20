"""Miscellaneous common utilities
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
__modifiedBy__ = "$modifiedBy: CCPN $"
__dateModified__ = "$dateModified: 2017-07-07 16:32:57 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.0 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================
"""Common utilities
"""

import datetime
import os
import random
import sys
import string
import itertools
from collections import Iterable
from . import Constants
from collections import OrderedDict


# Max value used for random integer. Set to be expressible as a signed 32-bit integer.
maxRandomInt = 2000000000

WHITESPACE_AND_NULL = {'\x00', '\t', '\n', '\r', '\x0b', '\x0c'}

# valid characters for file names
# NB string.ascii_letters and string.digits are not compatible
# with Python 2.1 (used in ObjectDomain)
defaultFileNameChar = '_'
separatorFileNameChar = '+'
validFileNamePartChars = ('abcdefghijklmnopqrstuvwxyz'
                          'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
                          + defaultFileNameChar)
validCcpnFileNameChars = validFileNamePartChars + '-.' + separatorFileNameChar


# # Not used - Rasmus 20/2/2017
# Sentinel = collections.namedtuple('Sentinel', ['value'])


def convertStringToFileName(fileNameString, validChars=validCcpnFileNameChars,
                            defaultChar=defaultFileNameChar):
    ll = [x for x in fileNameString]
    for ii, char in enumerate(ll):
        if char not in validChars:
            ll[ii] = defaultChar
    #
    return ''.join(ll)


def getCcpFileString(fileNameString):
    """
    Changes an input string to the one used for a component of file names.
    """

    return convertStringToFileName(fileNameString, validFileNamePartChars,
                                   defaultFileNameChar)


def incrementName(name):
    """Add '_1' to name or change suffix '_n' to '_(n+1) """
    ll = name.rsplit('_', 1)
    if len(ll) == 2:
        try:
            ll[1] = str(int(ll[1]) + 1)
            return '_'.join(ll)

        except ValueError:
            pass

    return name + '_1'


def recursiveImport(dirname, modname=None, ignoreModules=None, force=False):
    """ recursively import all .py files
    (not starting with '__' and not containing internal '.' in their name)
    from directory dirname and all its subdirectories, provided they
    contain '__init__.py'
    Serves to check that files compile without error

    modname is the module name (dot-separated) corresponding to the directory
    dirName.
    If modname is None, dirname must be on the pythonPath

    Note that there are potential problems if the files we want are not
    the ones encountered first on the pythonPath
    """

    # Must be imported here, as entire file must be importable from Python 2 NefIo
    from . import Path

    listdir = os.listdir(dirname)
    try:
        listdir.remove('__init__.py')
    except ValueError:
        if not force:
            return

    files = []

    if ignoreModules is None:
        ignoreModules = []

    if modname is None:
        prefix = ''
    else:
        prefix = modname + '.'

    listdir2 = []
    for name in listdir:
        head, ext = os.path.splitext(name)
        if (prefix + head) in ignoreModules:
            pass
        elif ext == '.py' and head.find('.') == -1:
            files.append(head)
        else:
            listdir2.append(name)

    # import directory and underlying directories
    if modname:
        # Note that files is never empty, so module is lowest level not toplevel
        for ff in files:
            try:
                __import__(modname, {}, {}, [ff])
            except:
                # We want log output, not an Exception in all cases here
                from .Logging import getLogger

                getLogger().warning("Import failed for %s.%s" % (modname, ff))

    for name in listdir2:
        newdirname = Path.joinPath(dirname, name)
        if os.path.isdir(newdirname) and name.find('.') == -1:
            recursiveImport(newdirname, prefix + name, ignoreModules)


def isWindowsOS():
    return sys.platform[:3].lower() == 'win'


def parseSequenceCode(value):
    """split sequence code into (seqCode,seqInsertCode, offset) tuple"""

    # sequenceCodePattern = re.compile('(\d+)?(.*?)(\+\d+|\-\d+)?$')

    tt = Constants.sequenceCodePattern.match(value.strip()).groups()

    if tt[0] is None and not tt[1]:
        # special case: entire string matches offset modifier and is misread
        return None, tt[2], None
    else:
        return (
            tt[0] and int(tt[0]),  # None or an integer
            tt[1],  # Text string, possibly empty
            tt[2] and int(tt[2]),  # None or an integer
            )


def splitIntFromChars(value):
    """convert a string with a leading integer optionally followed by characters
    into an (integer,string) tuple"""

    value = value.strip()

    for ii in reversed(range(1, len(value) + 1)):
        try:
            number = int(value[:ii])
            chars = value[ii:]
            break
        except ValueError:
            continue
    else:
        number = None
        chars = value

    return number, chars


def dictionaryProduct(dict1, dict2):
    """multiply input {a:x}, {b:y} to result {(a,b):x*y} dictionary"""
    result = {}
    for key1, val1 in dict1.items():
        for key2, val2 in dict2.items():
            result[(key1, key2)] = val1 * val2
    #
    return result


def uniquify(sequence):
    """Get list of unique elements in sequence, in order of first appearance"""
    seen = set()
    seen_add = seen.add
    return [x for x in sequence if x not in seen and not seen_add(x)]


def isClose(a, b, relTolerance=1e-05, absTolerance=1e-08):
    """Are a and b identical within reasonable floating point tolerance?
    Uses sum of relative (relTolerance) and absolute (absTolerance) difference

    Inspired by numpy.isclose()"""
    return (abs(a - b) <= (absTolerance + relTolerance * abs(b)))

# NOTE:ED not 2.7 compatible
def isIterable(obj):
# def isIterable(obj) -> bool:
    "Returns True if obj is iterable"
    try:
        iter(obj)
        return True
    except TypeError:
        pass
    return False


def getTimeStamp():
    """Get iso-formtted timestamp"""
    return datetime.datetime.today().isoformat()


def getUuid(programName, timeStamp=None):
    """Get UUid following the NEF convention"""
    if timeStamp is None:
        timeStamp = getTimeStamp()
    return '%s-%s-%s' % (programName, timeStamp, random.randint(0, maxRandomInt))


def name2IsotopeCode(name=None):
    """Get standard isotope code matching atom name or axisCode string

    """
    if not name:
        return None

    result = Constants.DEFAULT_ISOTOPE_DICT.get(name[0])
    if result is None:
        if name[0].isdigit():
            ss = name.title()
            for key in Constants.isotopeRecords:
                if ss.startswith(key):
                    if name[:len(key)].isupper():
                        result = key
                    break
        else:
            result = Constants.DEFAULT_ISOTOPE_DICT.get(name[:2])
    #
    return result


def isotopeCode2Nucleus(isotopeCode=None):
    if not isotopeCode:
        return None

    record = Constants.isotopeRecords.get(isotopeCode)
    if record is None:
        return None
    else:
        return record.symbol.upper()

    # for tag,val in sorted(Constants.DEFAULT_ISOTOPE_DICT.items()):
    #   if val == isotopeCode:
    #     return tag
    # else:
    #   return None


def name2ElementSymbol(name):
    """Get standard element symbol matching name or axisCode

    NB, the first letter takes precedence, so e.g. 'CD' returns 'C' (carbon)
    rather than 'CD' (Cadmium)"""

    # NB, We deliberately do NOT use 'value in Constants.DEFAULT_ISOTOPE_DICT'
    # We want to avoid elements that are in the dict but have value None.
    if not name:
        return None
    elif Constants.DEFAULT_ISOTOPE_DICT.get(name[0]) is not None:
        return name[0]
    elif Constants.DEFAULT_ISOTOPE_DICT.get(name[:2]) is not None:
        return name[:2]
    elif name[0].isdigit():
        ss = name.title()
        for key, record in Constants.isotopeRecords.items():
            if ss.startswith(key):
                if name[:len(key)].isupper():
                    return record.symbol.upper()
                break
    #
    return None


def checkIsotope(text):
    """Convert isotope specifier string to most probable isotope code - defaulting to '1H'

    This function is intended for external format isotope specifications, *not* for
    axisCodes or atom names, hence the difference to name2ElementSymbol.
    """
    defaultIsotope = '1H'
    name = text.strip().upper()

    if not name:
        return defaultIsotope

    if name in Constants.isotopeRecords:
        # Superfluous but should speed things up
        return name

    for isotopeCode in Constants.isotopeRecords:
        # NB checking this first means that e.g. 'H13C' returns '13C' rather than '1H'
        if isotopeCode.upper() in name:
            return isotopeCode

    # NB order of checking means that e.g. 'CA' returns Calcium rather than Carbon
    result = (Constants.DEFAULT_ISOTOPE_DICT.get(name[:2])
              or Constants.DEFAULT_ISOTOPE_DICT.get(name[0]))

    if result is None:
        if name == 'D':
            # special case
            result = '2H'
        else:
            result = defaultIsotope
    #
    return result


def axisCodeMatch(axisCode, refAxisCodes):
    """Get refAxisCode that best matches axisCode """
    for ii, indx in enumerate(_axisCodeMapIndices([axisCode], refAxisCodes)):
        if indx == 0:
            # We have a match
            return refAxisCodes[ii]
    else:
        return None


def axisCodeMapping(axisCodes, refAxisCodes):
    """get {axisCode:refAxisCode} mapping dictionary
    all axisCodes must match, or dictionary will be empty
    NB a series of single-letter axisCodes (e.g. 'N', 'HCN') can be passed in as a string"""
    result = {}
    mapIndices = _axisCodeMapIndices(axisCodes, refAxisCodes)
    if mapIndices:
        for ii, refAxisCode in enumerate(refAxisCodes):
            indx = mapIndices[ii]
            if indx is not None:
                result[axisCodes[indx]] = refAxisCode
    #
    return result


def reorder(values, axisCodes, refAxisCodes):
    """reorder values in axisCode order to refAxisCode order, by matching axisCodes

    NB, the result will be the length of refAxisCodes, with additional Nones inserted
    if this is longer than the values.

    NB if there are multiple matches possible, one is chosen by heuristics"""
    if len(values) != len(axisCodes):
        raise ValueError("Length mismatch between %s and %s" % (values, axisCodes))
    remapping = _axisCodeMapIndices(axisCodes, refAxisCodes)
    result = list(values[x] for x in remapping)
    #
    return result


def _axisCodeMapIndices(axisCodes, refAxisCodes):
    """get mapping tuple so that axisCodes[result[ii]] matches refAxisCodes[ii]
    all axisCodes must match, but result can contain None if refAxisCodes is longer
    if axisCodes contain duplicates, you will get one of possible matches"""

    #CCPNINTERNAL - used in multiple places to map display order and spectrum order

    lenDifference = len(refAxisCodes) - len(axisCodes)
    if lenDifference < 0:
        return None

    # Set up match matrix
    matches = []
    for code in axisCodes:
        matches.append([axisCodesCompare(code, x, mismatch=-999999) for x in refAxisCodes])

    # find best mapping
    maxScore = sum(len(x) for x in axisCodes)
    bestscore = -1
    results = None
    values = list(range(len(axisCodes))) + [None] * lenDifference
    for permutation in itertools.permutations(values):
        score = 0
        for ii, jj in enumerate(permutation):
            if jj is not None:
                score += matches[jj][ii]
        if score > bestscore:
            bestscore = score
            results = [permutation]
        elif score == maxScore:
            # it cannot get any higher
            results.append(permutation)
    #
    if results:
        # Pick the first of the equally good matches as the default answer
        result = results[0]
        if len(results) > 1:
            # Multiple matches - try to select on pairs of bound atoms
            # NB we do not need to be rigorous as this is only used to resolve ambiguity
            # NB this is necessary to do a correct match of e.g. Hc, Ch, H to Hc, Hc1, Ch1
            boundCodeDicts = []
            for tryCodes in axisCodes, refAxisCodes:
                boundCodes = {}
                boundCodeDicts.append(boundCodes)
                for ii, code in enumerate(tryCodes):
                    if len(code) > 1:
                        for jj in range(ii + 1, len(tryCodes)):
                            code2 = tryCodes[jj]
                            if len(code2) > 1:
                                if (code[0].isupper() and code[0].lower() == code2[1] and
                                        code2[0].isupper() and code2[0].lower() == code[1] and
                                        code[2:] == code2[2:]):
                                    # Matches pair of bound atoms - e.g. Hc1, Ch1
                                    boundCodes[tryCodes.index(code)] = tryCodes.index(code2)
                                    boundCodes[tryCodes.index(code2)] = tryCodes.index(code)

            if boundCodeDicts[0] and boundCodeDicts[1]:
                # bound pairs on both sides - check for matching pairs
                bestscore = -1
                for permutation in results:
                    score = 0
                    for idx1, idx2 in boundCodeDicts[1].items():
                        target = permutation[idx1]
                        if target is not None and target == boundCodeDicts[0].get(permutation[idx2]):
                            score += 1
                    if score > bestscore:
                        bestscore = score
                        result = permutation
    else:
        result = None
    #
    return result


def axisCodesCompare(code, code2, mismatch=0):
    """Score code, code2 for matching. Score is length of common prefix, or 'mismatch' if None"""

    if not code or not code2 or code[0] != code2[0]:
        score = mismatch
    elif code == code2:
        score = len(code)
    elif code[0].islower():
        # 'fidX...' 'delay', etc. must match exactly
        score = mismatch
    elif code.startswith('MQ'):
        # 'MQxy...' must match exactly
        score = mismatch
    elif len(code) == 1 or code[1].isdigit() or len(code2) == 1 or code2[1].isdigit():
        # Match against a single upper-case letter on one side. Always OK
        score = 1
    else:
        # Partial match of two strings with at least two significant chars each
        score = len(os.path.commonprefix((code, code2))) or mismatch
        if score == 1:
            # Only first letter matches, second does not
            if ((code.startswith('Hn') and code2.startswith('Hcn')) or
                    (code.startswith('Hcn') and code2.startswith('Hn'))):
                # Hn must matches Hcn
                score = 2
            else:
                # except as above we need at least two char match

                # TODO:ED check matching codes when each contain more than 1 character; dict?
                score = 3  #mismatch
        elif code.startswith('J') and score == 2:
            # 'Jab' matches 'J' or 'Jab...', but NOT 'Ja...'
            score = mismatch
    #
    return score


def doAxisCodesMatch(axisCodes, refAxisCodes):
    """Return True if axisCodes match refAxisCodes else False"""
    if len(axisCodes) != len(refAxisCodes):
        return False

    for ii, code in enumerate(axisCodes):
        if not axisCodesCompare(code, refAxisCodes[ii]):
            return False
    #
    return True


def stringifier(*fields, **options):
    """Get stringifier function, that will format an object x according to

    <str(x): field1=x.field1, field2=x.field2, ...>

    All floating point values encountered will be formatted according to floatFormat"""

    # Unfortunately necessary as this package must be read from v2io
    # and python 2 does not have keyword-only arguments
    # What we should do is the function definition below:
    # def stringifier(*fields, floatFormat=None):
    if 'floatFormat' in options:
        floatFormat = options.pop('floatFormat')
    else:
        floatFormat = None
    if options:
        raise ValueError("Unknown options: %s" % ', '.join(sorted(options.keys())))

    # Proper body of function starts here
    if floatFormat is None:
        # use default formatter, avoiding continuous creation of new ones
        localFormatter = stdLocalFormatter
    else:
        localFormatter = LocalFormatter(overrideFloatFormat=floatFormat)

    fieldFormats = []
    for field in fields:
        # String will be 'field1={_obj.field1}
        fieldFormats.append('{0}={{_obj.{0}!r}}'.format(field))

    formatString = '<{_obj.pid!s}| ' + ', '.join(fieldFormats) + '>'

    def formatter(x):
        return localFormatter.format(format_string=formatString, _obj=x)

    return formatter


def contains_whitespace(s):
    return True in [c in s for c in string.whitespace]


def resetSerial(apiObject, newSerial):
    """ADVANCED Reset serial of object to newSerial, resetting parent link
    and the nextSerial of the parent.

    Raises ValueError for objects that do not have a serial
    (or, more precisely, where the _wrappedData does not have a serial)."""

    # NB, needed both from V2 NefIo and V3, hence putting nit here,
    # even though it uses V2 objects

    if not hasattr(apiObject, 'serial'):
        raise ValueError("Cannot reset serial, %s does not have a 'serial' attribute"
                         % apiObject)
    downlink = apiObject.__class__._metaclass.parentRole.otherRole.name

    parentDict = apiObject.parent.__dict__
    downdict = parentDict[downlink]
    oldSerial = apiObject.serial
    serialDict = parentDict['_serialDict']

    if newSerial == oldSerial:
        return

    elif newSerial in downdict:
        raise ValueError("Cannot reset serial to %s - value already in use" % newSerial)

    else:
        maxSerial = serialDict[downlink]
        apiObject.__dict__['serial'] = newSerial
        downdict[newSerial] = apiObject
        del downdict[oldSerial]
        if newSerial > maxSerial:
            serialDict[downlink] = newSerial
        elif oldSerial == maxSerial:
            serialDict[downlink] = max(downdict)


def makeIterableList(inList=None):
    """
    Take a list of lists and concatenate into a single list.
    Remove any Nones from the list
    :param inList:
    :return single list:
    """
    if isinstance(inList, Iterable) and not isinstance(inList, str):
        return [y for x in inList for y in makeIterableList(x) if inList]
    else:
        if inList is not None:
            return [inList]
        else:
            return []


def _truncateText(text, splitter=' , ', maxWords=4):
    "Splits the text by the given splitter. If more then maxWords, it return the maxWord plus dots, otherwise just the text"
    words = text.split(splitter)
    if len(words) > maxWords:
        return splitter.join(words[:maxWords]) + ' ...'
    else:
        return text


def _traverse(obj, tree_types=(list, tuple)):
    """
    used to flat the state in a long list
    """
    if isinstance(obj, tree_types):
        for value in obj:
            for subvalue in _traverse(value, tree_types):
                yield subvalue
    else:
        yield obj


def _getChildren(obj, path=None):
    """
    Walks in a tree like obj and put all children/parents in list of list eg: [[Parent,child...,],...]
    """
    children = []
    if path is None:
        path = []
    path.append(obj)
    if obj._childClasses:
        for att in obj._childClasses:
            for child in getattr(obj, att._pluralLinkName):
                children.extend(_getChildren(child, path[:]))
    else:
        children.append(path)
    return children


def percentage(percent, whole):
    return (percent * whole) / 100.0


def splitDataFrameWithinRange(dataframe, column1, column2, minX, maxX, minY, maxY):
    """
    :param dataframe: dataframe with index a pid type, columns str, values floats or ints
    :param column1: label1 , eg PC1
    :param column2: label1 , eg PC2
    :param minX:  min value for Y
    :param maxX:  Max value for X
    :param minY: min value for Y
    :param maxY: max value for Y
    :return:  inners  a dataframe like the unput  but containing only the values within the ranges  and
              outers (rest) not included in inners
    """

    bools = dataframe[column1].between(minX, maxX, inclusive=True) & dataframe[column2].between(minY, maxY, inclusive=True)
    inners = dataframe[bools]
    outers = dataframe[-bools]
    filteredInners = inners.filter(items=[column1, column2])
    filteredOuters = outers.filter(items=[column1, column2])

    return filteredInners, filteredOuters


class LocalFormatter(string.Formatter):
    """Overrides the string formatter to change the float formatting"""

    def __init__(self, overrideFloatFormat='.6g'):
        super(LocalFormatter, self).__init__()
        self.overrideFloatFormat = overrideFloatFormat

    def convert_field(self, value, conversion):
        # do any conversion on the resulting object
        # NB, conversion parameter is not used

        if hasattr(value, 'pid'):
            return str(value)
        elif isinstance(value, float):
            return format(value, self.overrideFloatFormat)
        elif type(value) == tuple:
            # Deliberate. We do NOT want to catch tuple subtypes here
            end = ',)' if len(value) == 1 else ')'
            return '(' + ', '.join(self.convert_field(x, 'r') for x in value) + end
        elif type(value) == list:
            # Deliberate. We do NOT want to catch list subtypes here
            return '[' + ', '.join(self.convert_field(x, 'r') for x in value) + ']'
        elif conversion is None:
            return value
        elif conversion == 's':
            return str(value)
        elif conversion == 'r':
            return repr(value)
        elif conversion == 'a':
            try:
                return ascii(value)
            except NameError:
                # Likely we are in Python 2.
                # As ascii behaves like Python 2 repr, this should be the correct workaround
                return repr(value)
        raise ValueError("Unknown conversion specifier {0!s}".format(conversion))


stdLocalFormatter = LocalFormatter()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 20190507:ED new routines to match axis codes and return dict or indices

# NOTE:ED not 2.7 compatible
def _matchSingleAxisCode(code1=None, code2=None, exactMatch=False, allowLowercase=True):
# def _matchSingleAxisCode(code1: str = None, code2: str = None, exactMatch: bool = False, allowLowercase=True) -> int:
    """number of matching characters
    code1, code2 = strings
    e.g. 'Hn1', 'H1'

    Compare single axis codes

    Must always be upper case first letter

    more matching letters = higher code
    difference in length reduces match

    'Jab' matches 'J' or 'Jab...', but NOT 'Ja...'      ie, 1 or 3 or more letter match

    MQ gets no match

    Hn* always matches Hcn*

    :param code1: first axis code to compare
    :param code2: second axis code to compare
    :param exactMatch: only allow exact matches, True/False
    :return: score based on the match
    """
    # undefined codes
    if not code1 or not code2 or (code1[0].islower() and code2[0].islower() and not allowLowercase):
        return 0

    # if exactMatch is True then only test for exact match
    if exactMatch:
        return code1 == code2

    ms = [a for a in zip(code1, code2)]  # zips to the shortest string
    ss = 0

    # add extra tests from v2.4
    if code1.startswith('MQ') or code2.startswith('MQ'):
        return 0
    # char followed by digit already accounted for

    # get count of matching characters - more characters -> higher score
    for a, b in ms:
        if a != b:
            break
        ss += 1

    # another v2.4 test
    if ss:
        if ((code1.startswith('Hn') and code2.startswith('Hcn')) or
                (code1.startswith('Hcn') and code2.startswith('Hn'))):
            # Hn must always match Hcn, give it a high score
            ss += 500

        if code1.startswith('J'):
            if ss == 2:  # must be a 1, or (3 or more) letter match
                return 0

        ss += _matchSingleAxisCodeLength(code1, code2)
    return (1000 + ss) if ss else 0


def _matchSingleAxisCodeLength(code1, code2):
    """return a score based on the mismatch in length
    """
    lenDiff = abs(len(code1) - len(code2))

    return (100 + 800 // (lenDiff + 1))


def _SortByMatch(item):
    """quick sorting key for axisCode match tuples
    """
    return -item[2]  # sort from high to low


# NOTE:ED not 2.7 compatible
def getAxisCodeMatch(axisCodes, refAxisCodes, allMatches=False, exactMatch=False):
# def getAxisCodeMatch(axisCodes, refAxisCodes, allMatches=False, exactMatch=False) -> OrderedDict:
    """Return an OrderedDict containing the mapping from the refAxisCodes to axisCodes

    There may be multiple matches, or None for each axis code.

    Set allMatches to True to return all, or False for only the best match in each case

    e.g. for unique axis codes:

        getAxisCodeMatch(('Hn', 'Nh', 'C'), ('Nh', 'Hn'), allMatches=False)

        ->  { 'Hn':   'Hn'
              'Nh':   'Nh'
              'C' :   None
            }

        getAxisCodeMatch(('Hn', 'Nh', 'C'), ('Nh', 'Hn'), allMatches=True)

        ->  { 'Hn':   ('Hn',)
              'Nh':   ('Nh',)
              'C' :   ()
            }

    for similar repeated axis codes, possibly from matching isotopeCodes:

        getAxisCodeMatch(('Nh', 'H'), ('H', 'H1', 'N'), allMatches=True)

        ->  { 'Nh':   'N'
              'H' :   'H'
            }

        getAxisCodeMatch(('Nh', 'H'), ('H', 'H1', 'N'), allMatches=True)

        ->  { 'Nh':   ('N',)
              'H' :   ('H', 'H1')       <- in this case the first match is always the highest
            }
    """

    found = OrderedDict()
    for ii, code1 in enumerate(axisCodes):
        foundCodes = []
        for jj, code2 in enumerate(refAxisCodes):

            match = _matchSingleAxisCode(code1, code2, exactMatch=exactMatch)
            if match:
                foundCodes.append((code2, jj, match))

        if allMatches:
            found[code1] = tuple(mm[0] for mm in sorted(foundCodes, key=_SortByMatch))
        else:
            found[code1] = sorted(foundCodes, key=_SortByMatch)[0][0] if foundCodes else None

    return found


def getAxisCodeMatchIndices(axisCodes, refAxisCodes, exactMatch=False):
    """Return a tuple containing the indices for each axis code in axisCodes in refAxisCodes

    Only the best match is returned for each code, elements not found in refAxisCodes will be marked as 'None'

    e.g. for unique axis codes:

        indices = getAxisCodeMatchIndices(('Hn', 'Nh', 'C'), ('Nh', 'Hn'))

        ->  (1, 0, None)

                i.e axisCodes[0] = 'Hn' which maps to refAxisCodes[indices[0]] = 'Hn'

    for similar repeated axis codes, possibly from matching isotopeCodes:

        getAxisCodeMatchIndices(('Nh', 'H'), ('H', 'H1', 'N'))

        ->  (2, 0)

    """

    found = []
    for ii, code1 in enumerate(axisCodes):
        foundCodes = []
        for jj, code2 in enumerate(refAxisCodes):

            match = _matchSingleAxisCode(code1, code2, exactMatch=exactMatch)
            if match:
                foundCodes.append((code2, jj, match))

        found.append(sorted(foundCodes, key=_SortByMatch)[0][1] if foundCodes else None)

    return tuple(found)

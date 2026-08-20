"""
compareNef - a series of routines to compare the contents of two Nef files

Command Line Usage:
  compareNef for execution from the command line with a suitable script
  An example can be found in AnalysisV3/bin/compareNef:

      #!/usr/bin/env sh
      export CCPNMR_TOP_DIR="$(dirname $(cd $(dirname "$0"); pwd))"
      export ANACONDA3=${CCPNMR_TOP_DIR}/miniconda
      export PYTHONPATH=${CCPNMR_TOP_DIR}/src/python:${CCPNMR_TOP_DIR}/src/c
      ${ANACONDA3}/bin/python ${CCPNMR_TOP_DIR}/src/python/ccpn/util/nef/compareNef.py $*

  Usage:  compareNef [options]

  optional arguments:
    -h, --help              show this help message
    -H, --Help              Show detailed help

    -i, --ignoreblockname   Ignore the blockname when comparing two Nef files
                            May be required when converting Nef files through
                            different applications.
                            May be used with -f and -b.

    -f inFile1 inFile2, --file inFile1 inFile2
                            Compare two Nef files and print the results to the
                            screen.

    -b inDir1 inDir2 outDir, --block inDir1 inDir2 outDir
                            compare Nef files common to directories
                            inDir1 and inDir2. Write output *.txt for each
                            file into the outDir directory.

    -s, --screen            Output batch processing to screen, default is to .txt files
                            may be used with -b

    -o, --overwrite         Overwrite existing .txt files. If false then files are
                            appended with '(n)' before the extension, where n is
                            the next available number.

    -c, --create            Automatically create directories as required.

    -I, --ignorecase        Ignore case when comparing items.

Details of the contents of Nef files can be found in GenericStarParser
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

Module Contents
===============

compareNef contains the following routines:
  compareNefFiles      compare two Nef files and return a compare object

    Searches through all objects: dataExtents, dataBlocks, saveFrames and Loops within the files.
    Comparisons are made for all data structures that have the same name.
    Differences for items within a column are listed in the form:
      dataExtent:dataBlock:saveFrame:Loop:  <Column>: columnName  <rowIndex>: row  -->  value1 != value2

    dataExtents, dataBlocks, saveFrames, Loops, columns present only in one file are listed, in the form:
      dataExtent:dataBlock:saveFrame:Loop: contains --> parameter1
                                                        parameter2
                                                        ...
                                                        parameterN

  A compare object is a list of nefItems of the form:

  ::

      NefItem
        inWhich         a flag labelling which file the item was found in
                        1 = found in first file, 2 = found in second file, 3 = common to both
        List
          Item          multiple strings containing the comparison tree
          (,List)       the last item of which may be a list of items common to the tree

  e.g., for parameters present in the first file:
        [
          inWhich=1
          list=[dataExtent1, dataBlock1, saveFrame1, Loop1, [parameter1, parameter2, parameter3]]
        ]

  compareDataExtents    compare two DataExtent objects and return a compare list as above
  compareDataBlocks     compare two DataBlock objects and return a compare list as above
  compareSaveFrames     compare two SaveFrame objects and return a compare list as above
  compareLoops          compare two Loop objects and return a compare list as above

  compareNefFiles       compare two Nef files and return a compare list as above
  batchCompareNefFiles  compare two directories of Nef files.
                        Nef Files common to specified directories are compared and the compare
                        lists are written to the third directory as .txt

  printCompareList      print the compare list to the screen
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


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
__dateModified__ = "$dateModified: 2020-01-14 11:49:36 +0000 (Tue, January 14, 2020) $"
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
import copy
import sys


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this is a fix to get the import to work when running as a standalone
# when importing into your own code, with PYTHON_PATH defined it can be safely removed

def import_parents(level=1):
    global __package__

    import sys
    from os import path
    import importlib

    # pathlib does all this a lot nicer, but don't think it's in python2.7
    top = parent = path.dirname(path.abspath(__file__))
    package = []
    for t in range(level):
        package.insert(0, os.path.basename(top))
        top = path.dirname(top)

    sys.path.append(str(top))
    try:
        sys.path.remove(str(parent))
    except ValueError:  # already removed
        pass

    __package__ = str('.'.join(package))
    importlib.import_module(__package__)


if __name__ == '__main__' and __package__ is None:
    import_parents(level=1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from . import GenericStarParser, StarIo
from .SafeOpen import safeOpen
import unittest
from ast import literal_eval
from os import listdir
from os.path import isfile, join
import re


try:
    # Python 3
    from itertools import zip_longest
except:
    # python 2.7
    from itertools import izip_longest as zip_longest

DATAEXTENT = ''
DATABLOCK = ''
SAVEFRAME = ''
LOOP = ''
COLUMN = ''


def defineArguments():
    """Define the arguments of the program

    :return argparse instance
    """
    import argparse

    parser = argparse.ArgumentParser(description='Compare the contents of Nef files', prog='compareNef', usage='%(prog)s [options]',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-H', '--Help', dest='help', action='store_true', default=False, help='Show detailed help')
    parser.add_argument('-i', '--ignoreblockname', dest='ignoreBlockName', action='store_true', default=False,
                        help='Ignore the blockname when comparing two Nef files')
    parser.add_argument('-f', '--file', dest='inFiles', nargs=2, metavar=('inFile1', 'inFile2'), default=None,
                        help='Compare two Nef files and print the results to the screen')
    parser.add_argument('-b', '--block', dest='blockDirs', nargs=3, metavar=('inDir1', 'inDir2', 'outDir'), default=None,
                        help='Batch mode: compare the contents of two directories')
    parser.add_argument('-s', '--screen', dest='screen', action='store_true', default=False, help='Output batch processing to screen')
    parser.add_argument('-o', '--overwrite', dest='overwriteExisting', action='store_true', default=False,
                        help='Overwrite existing .txt files. If false, new files are appended with "(n)"')
    parser.add_argument('-c', '--create', dest='createDirs', action='store_true', default=False,
                        help='Create directories as required')
    parser.add_argument('-I', '--ignorecase', dest='ignoreCase', action='store_true', default=False,
                        help='Ignore case when comparing items.')

    return parser


#=========================================================================================
# nefItem
#=========================================================================================

class nefItem(object):
    """Holds the contents of a single Nef comparison
    inWhich   a flag labelling which file the item was found in
              1 = found in the first file, 2 = found on the second file, 3 = common to both
    list      a list of strings containing the comparison information
    """

    def __init__(self, cItem=None):
        self.inWhich = None
        self.list = []

        if cItem is not None:
            self.list = copy.deepcopy(cItem.list)
            self.inWhich = cItem.inWhich


#=========================================================================================
# _loadGeneralFile
#=========================================================================================

def _loadGeneralFile(path=None):
    """Load a file with the given pathname and return a dict of the contents

    :return entry:dict
    """
    usePath = path if path.startswith('/') else os.path.join(os.getcwd(), path)
    entry = StarIo.parseNefFile(usePath)  # 'lenient')
    print(' %s' % path)
    return entry


#=========================================================================================
# printFile
#=========================================================================================

def printFile(thisFile):
    """Print a file to the screen
    """
    print('~' * 80)
    print(thisFile)
    for i, val in enumerate(thisFile):
        print(i, thisFile[val])
        sub = thisFile[val]
        for j, valj in enumerate(sub):
            print('  ', j, sub[valj])
            if j > 3:
                break

            sub2 = sub[valj]
            for k, valk in enumerate(sub2):
                loopType = sub2[valk]
                if isinstance(loopType, GenericStarParser.Loop):
                    print('    ', k, 'LOOP', loopType)
                else:
                    print('    ', k, loopType)

                if k > 3:
                    break


#=========================================================================================
# sizeNefList
#=========================================================================================

def sizeNefList(nefList, whichType=0):
    """List only those items that are of type whichType

    :param nefList: list to print
    :param whichType: type to print
    """
    count = 0
    if nefList is not None:
        for cCount, cc in enumerate(nefList):
            if cc.inWhich == whichType:
                count += 1
    return count


#=========================================================================================
# printWhichList
#=========================================================================================

def printWhichList(nefList, whichType=0):
    """List only those items that are of type whichType

    :param nefList: list to print
    :param whichType: type to print
    """
    for cCount, cc in enumerate(nefList):
        if cc.inWhich == whichType:

            if isinstance(cc.list[-1], str):
                print('  ' + ':'.join(cc.list[:]))
            else:
                outStr = '  ' + ':'.join(cc.list[:-1]) + ': contains --> '
                lineTab = '\n' + ' ' * len(outStr)
                print(outStr + lineTab.join(cc.list[-1]))


#=========================================================================================
# printCompareList
#=========================================================================================

def printCompareList(nefList, inFile1, inFile2):
    """Print the contents of the nef compare list to the screen

    Output is in three parts:
      - items that are present only in the first file
      - items that are only in the second file
      - differences between objects that are common in both files

    :param nefList: list to print
    :param inFile1: name of the first file
    :param inFile2: name of the second file
    """

    if not isinstance(inFile1, str):
        print('TypeError: inFile1 must be a string.')
        return
    if not isinstance(inFile2, str):
        print('TypeError: inFile2 must be a string.')
        return

    # print the items that are only present in the first nefFile
    if sizeNefList(nefList, whichType=1) > 0:
        print('\nItems that are only present in ' + inFile1 + ':')
        printWhichList(nefList, 1)

    # print the items that are only present in the second nefFile
    if sizeNefList(nefList, whichType=2) > 0:
        print('\nItems that are only present in ' + inFile2 + ':')
        printWhichList(nefList, 2)

    # print the common items
    if sizeNefList(nefList, whichType=3) > 0:
        print('\nItems that are present in both files:')
        printWhichList(nefList, 3)


#=========================================================================================
# _filterName
#=========================================================================================

def _filterName(inName):
    """Remove rogue `n` quotes from the names.
    (This is currently only a test)

    :param inName:
    :return:
    """
    # ejb - need to remove the rogue `n` at the beginning of the name if it exists
    #       as it is passed into the namespace and gets added iteratively every save
    #       next three lines remove all occurrences of `n` from name
    regex = u'\`\d*`+?'
    return re.sub(regex, '', inName)  # substitute with ''


#=========================================================================================
# addToList
#=========================================================================================

def addToList(inList, cItem, nefList):
    """Append cItem to the compare list
    Currently adds one cItem with a list as the last element

    :param inList: a list of items to add to the end of cItem
    :param cItem: object containing the current tree to add to the list
    :param nefList: current list of comparisons
    :return: list of type nefItem
    """
    if len(inList) > 0:
        cItem3 = copy.deepcopy(cItem)
        cItem3.list.append(list(inList))  # nest the list within the cItem
        nefList.append(nefItem(cItem=cItem3))

    return nefList


#=========================================================================================
# _compareDicts
#=========================================================================================

def _compareDicts(dd1, dd2, options):
    """Compare the contents of two dictionaries

    :param dd1: dict containing keys from first nef file
    :param dd2:  dict containing keys from second nef file
    :param options: nameSpace holding the commandLineArguments
    :return:
    """
    if dd1 == dd2:
        return True

    if options.ignoreCase:
        # simple compare - make lowercase dictionaries and check equivalent
        # not perfect as may have mix of upper/lowercase keys that are the same
        try:
            lowerDd1 = literal_eval(str(dd1).lower())
            lowerDd2 = literal_eval(str(dd2).lower())

            if len(dd1.keys()) != len(lowerDd1.keys()):
                return False
            if len(dd2.keys()) != len(lowerDd2.keys()):
                return False

            if lowerDd1 == lowerDd2:
                return True

        except (SyntaxError, ValueError, AssertionError):
            # trap errors generated from a bad literal_eval
            return False

    return False


#=========================================================================================
# compareLoops
#=========================================================================================

def compareLoops(loop1, loop2, options, cItem=None, nefList=None):
    """Compare two Loops

    :param loop1: first Loop object, of type GenericStarParser.Loop
    :param loop2: second Loop object, of type GenericStarParser.Loop
    :param options: nameSpace holding the commandLineArguments
    :param cItem: list of str describing differences between nefItems
    :param nefList: input of nefItems
    :return: list of type nefItem
    """
    if cItem is None:
        cItem = nefItem()
    if nefList is None:
        nefList = []

    lSet = [bl for bl in loop1.columns]
    rSet = [bl for bl in loop2.columns]
    inLeft = set(lSet).difference(rSet)
    dSet = set(lSet).intersection(rSet)
    inRight = set(rSet).difference(lSet)

    cItem1 = copy.deepcopy(cItem)
    cItem1.list.append(LOOP + loop1.name)
    cItem1.inWhich = 1
    addToList(inLeft, cItem=cItem1, nefList=nefList)

    cItem2 = copy.deepcopy(cItem)
    cItem2.list.append(LOOP + loop2.name)
    cItem2.inWhich = 2
    addToList(inRight, cItem=cItem2, nefList=nefList)

    if loop1.data and loop2.data:
        rowRange = min(len(loop1.data), len(loop2.data))

        if len(loop1.data) != len(loop2.data):  # simple compare, same length tables
            cItem3 = copy.deepcopy(cItem)
            cItem3.list.append(LOOP + loop1.name)
            cItem3.list.append([' <rowLength>:  ' + str(len(loop1.data)) + ' != ' + str(len(loop2.data))])
            cItem3.inWhich = 3
            nefList.append(nefItem(cItem=cItem3))

        # carry on and compare the common table

        for compName in dSet:
            for rowIndex in range(rowRange):

                loopValue1 = loop1.data[rowIndex][compName]
                loopValue2 = loop2.data[rowIndex][compName]

                if not ((loopValue1 == loopValue2) or
                        ((str(loopValue1).lower() == str(loopValue2).lower()) and options.ignoreCase)):

                    # The value_strings are different
                    # Check to see if they are dictionaries
                    # and compare contents

                    try:
                        # these 2 lines will crash as a ValueError: malformed string if cannot
                        # be evaluated as:
                        #   - strings, numbers, tuples, lists, dicts, booleans, and None.

                        loopValue1 = literal_eval(loopValue1)
                        loopValue2 = literal_eval(loopValue2)

                        if isinstance(loopValue1, dict) and isinstance(loopValue2, dict):
                            if not _compareDicts(loopValue1, loopValue2, options):
                                _addItem(cItem, compName, loop1, loopValue1, loopValue2, nefList, rowIndex, inWhich=3)

                        else:
                            # not both dicts so compare is applicable
                            _addItem(cItem, compName, loop1, loopValue1, loopValue2, nefList, rowIndex, inWhich=3)

                    except (SyntaxError, ValueError, AssertionError):

                        # loopvalues cannot be converted to proper values
                        # need to check that comments are being loaded correctly
                        _addItem(cItem, compName, loop1, loopValue1, loopValue2, nefList, rowIndex, inWhich=3)

                else:

                    # nothing for the minute as identical already but may want to keep a log
                    pass

        #TODO
        # need to add a further test here, could do a diff on the tables which would pick up
        # insertions to the table - the columns would need to be reordered for this to work
        # what if there are a different number of columns?
        # also check for Mandatory items

    else:
        # can't compare non-existent loopdata
        if loop1.data is None:
            cItem3 = copy.deepcopy(cItem)
            cItem3.list.append(LOOP + loop1.name)
            cItem3.list.append([' <Contains no data>'])
            cItem3.inWhich = 1
            nefList.append(nefItem(cItem=cItem3))
        if loop2.data is None:
            cItem3 = copy.deepcopy(cItem)
            cItem3.list.append(LOOP + loop2.name)
            cItem3.list.append([' <Contains no data>'])
            cItem3.inWhich = 2
            nefList.append(nefItem(cItem=cItem3))

    return nefList


#=========================================================================================
# _addItem to the NefList or append to existing
#=========================================================================================

def _addItem(cItem, compName, loop1, loopValue1, loopValue2, nefList, rowIndex, inWhich):
    """Check the list of already added items and append to the end OR create a new item
    """
    cItem3 = copy.deepcopy(cItem)

    # iterate through existing items
    for itm in nefList:
        for a, b in zip_longest(cItem3.list[:] + [LOOP + loop1.name],
                                itm.list[:-1]):

            # check that the tree of saveFrame names matches
            if a != b:
                break
        else:
            # check that it is the correct frame type (1=inleft only, 2=inRight only, 3=inBoth)
            if itm.inWhich == inWhich:
                itm.list[-1].append(' <Column>: ' + compName + '  <rowIndex>: ' \
                                    + str(rowIndex) + '  -->  ' \
                                    + str(loopValue1) + ' != ' \
                                    + str(loopValue2))
                break

    else:
        # create a new item
        cItem3.list.append(LOOP + loop1.name)
        cItem3.list.append([' <Column>: ' + compName + '  <rowIndex>: ' \
                            + str(rowIndex) + '  -->  ' \
                            + str(loopValue1) + ' != ' \
                            + str(loopValue2)])
        cItem3.inWhich = inWhich
        nefList.append(nefItem(cItem=cItem3))


#=========================================================================================
# compareSaveFrames
#=========================================================================================

def compareSaveFrames(saveFrame1, saveFrame2, options, cItem=None, nefList=None):
    """Compare two saveFrames, if they have the same name then check their contents

    :param saveFrame1: first SaveFrame object, of type GenericStarParser.SaveFrame
    :param saveFrame2: second SaveFrame object, of type GenericStarParser.SaveFrame
    :param options: nameSpace holding the commandLineArguments
    :param cItem: list of str describing differences between nefItems
    :param nefList: input of nefItems
    :return: list of type nefItem
    """
    if cItem is None:
        cItem = nefItem()
    if nefList is None:
        nefList = []

    lSet = [' ' if not isinstance(saveFrame1[bl], GenericStarParser.Loop) else saveFrame1[bl].name for bl in saveFrame1]
    rSet = [' ' if not isinstance(saveFrame2[bl], GenericStarParser.Loop) else saveFrame2[bl].name for bl in saveFrame2]
    inLeft = set(lSet).difference(rSet).difference({' '})
    dSet = set(lSet).intersection(rSet).difference({' '})
    inRight = set(rSet).difference(lSet).difference({' '})

    lVSet = [str(bl) if not isinstance(saveFrame1[bl], GenericStarParser.Loop) else ' ' for bl in saveFrame1]
    rVSet = [str(bl) if not isinstance(saveFrame2[bl], GenericStarParser.Loop) else ' ' for bl in saveFrame2]
    inVLeft = set(lVSet).difference(rVSet).difference({' '})
    dVSet = set(lVSet).intersection(rVSet).difference({' '})
    inVRight = set(rVSet).difference(lVSet).difference({' '})

    # list everything only present in the first saveFrame

    cItem1 = copy.deepcopy(cItem)
    cItem1.list.append(SAVEFRAME + saveFrame1.name)
    cItem1.inWhich = 1
    addToList(inLeft, cItem=cItem1, nefList=nefList)
    addToList(inVLeft, cItem=cItem1, nefList=nefList)

    # list everything only present in the second saveFrame

    cItem2 = copy.deepcopy(cItem)
    cItem2.list.append(SAVEFRAME + saveFrame2.name)
    cItem2.inWhich = 2
    addToList(inRight, cItem=cItem2, nefList=nefList)
    addToList(inVRight, cItem=cItem2, nefList=nefList)

    # compare the common items

    cItem3 = copy.deepcopy(cItem)
    cItem3.list.append(SAVEFRAME + saveFrame1.name)
    cItem3.inWhich = 3
    for compName in dSet:
        # compare the loop items of the matching saveFrames

        compareLoops(saveFrame1[compName], saveFrame2[compName], options, cItem=cItem3, nefList=nefList)

    for compName in dVSet:
        # compare the other items in the saveFrames
        #   mandatory/optional parameters

        if saveFrame1[compName] != saveFrame2[compName]:
            cItem3 = copy.deepcopy(cItem)
            cItem3.list.append(SAVEFRAME + saveFrame2.name)
            cItem3.list.append([' <Value>:  ' + compName + '  -->  ' \
                                + str(saveFrame1[compName]) + ' != ' \
                                + str(saveFrame2[compName])])
            cItem3.inWhich = 3
            nefList.append(nefItem(cItem=cItem3))

    return nefList


#=========================================================================================
# compareDataBlocks
#=========================================================================================

def compareDataBlocks(dataBlock1, dataBlock2, options, cItem=None, nefList=None):
    """Compare two dataBlocks, if they have the same name then check their contents

    :param dataBlock1: first DataBlock object, of type GenericStarParser.DataBlock
    :param dataBlock2: second DataBlock object, of type GenericStarParser.DataBlock
    :param options: nameSpace holding the commandLineArguments
    :param cItem: list of str describing differences between nefItems
    :param nefList: input of nefItems
    :return: list of type nefItem
    """
    if cItem is None:
        cItem = nefItem()
    if nefList is None:
        nefList = []

    lSet = [dataBlock1[bl].name for bl in dataBlock1]
    rSet = [dataBlock2[bl].name for bl in dataBlock2]
    inLeft = set(lSet).difference(rSet)
    dSet = set(lSet).intersection(rSet)
    inRight = set(rSet).difference(lSet)

    # list everything only present in the first DataBlock

    cItem1 = copy.deepcopy(cItem)
    cItem1.list.append(DATABLOCK + dataBlock1.name)
    cItem1.inWhich = 1
    addToList(inLeft, cItem=cItem1, nefList=nefList)

    # list everything only present in the second DataBlock

    cItem2 = copy.deepcopy(cItem)
    cItem2.list.append(DATABLOCK + dataBlock2.name)
    cItem2.inWhich = 2
    addToList(inRight, cItem=cItem2, nefList=nefList)

    # compare the common items - strictly there should only be one DataBlock

    cItem3 = copy.deepcopy(cItem)
    cItem3.list.append(DATABLOCK + dataBlock1.name)
    cItem3.inWhich = 3
    for compName in dSet:
        compareSaveFrames(dataBlock1[compName], dataBlock2[compName], options, cItem=cItem3, nefList=nefList)

    return nefList


#=========================================================================================
# compareDataExtents
#=========================================================================================

def compareDataExtents(dataExt1, dataExt2, options, cItem=None, nefList=None):
    """Compare two dataExtents, if they have the same name then check their contents

    :param dataExt1: first DataExtent object, of type GenericStarParser.DataExtent
    :param dataExt2: second DataExtent object, of type GenericStarParser.DataExtent
    :param options: nameSpace holding the commandLineArguments
    :param cItem: list of str describing differences between nefItems
    :param nefList: input of nefItems
    :return: list of type nefItem
    """
    if cItem is None:
        cItem = nefItem()
    if nefList is None:
        nefList = []

    lSet = [dataExt1[bl].name for bl in dataExt1]
    rSet = [dataExt2[bl].name for bl in dataExt2]
    inLeft = set(lSet).difference(rSet)
    dSet = set(lSet).intersection(rSet)
    inRight = set(rSet).difference(lSet)

    # list everything only present in the first DataExtent

    cItem1 = copy.deepcopy(cItem)
    cItem1.list = [DATAEXTENT + dataExt1.name]
    cItem1.inWhich = 1  # left
    addToList(inLeft, cItem=cItem1, nefList=nefList)

    # list everything only present in the second DataExtent

    cItem2 = copy.deepcopy(cItem)
    cItem2.list = [DATAEXTENT + dataExt2.name]
    cItem2.inWhich = 2  # right
    addToList(inRight, cItem=cItem2, nefList=nefList)

    # compare the common items - strictly there should only be one DataExtent

    cItem3 = copy.deepcopy(cItem)
    cItem3.list = [DATAEXTENT + dataExt1.name]
    cItem3.inWhich = 3  # both
    for compName in dSet:
        compareDataBlocks(dataExt1[compName], dataExt2[compName], options, cItem=cItem3, nefList=nefList)

    return nefList


#=========================================================================================
# compareFiles
#=========================================================================================

def compareNefFiles(inFile1, inFile2, options, cItem=None, nefList=None):
    """Compare two Nef files and return comparison as a nefItem list

    :param inFile1: name of the first file
    :param inFile2: name of the second file
    :param options: nameSpace holding the commandLineArguments
    :param cItem: list of str describing differences between nefItems
    :param nefList: input of nefItems
    :return: list of type nefItem
    """
    if cItem is None:
        cItem = nefItem()
    if nefList is None:
        nefList = []

    if not os.path.isfile(inFile1):
        print('File Error:', inFile1)
    elif not os.path.isfile(inFile2):
        print('File Error:', inFile2)
    else:
        try:
            NefData1 = _loadGeneralFile(path=inFile1)
        except Exception as e:
            print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e), e)
            return None

        try:
            NefData2 = _loadGeneralFile(path=inFile2)
        except Exception as e:
            print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e), e)
            return None

        if options.ignoreBlockName is False:
            compareDataExtents(NefData1, NefData2, options, cItem=cItem, nefList=nefList)
        else:

            # assumes that there is only one block in a file
            # but this may change

            compList1 = [cn for cn in NefData1]
            compList2 = [cn for cn in NefData2]
            compareDataBlocks(NefData1[compList1[0]], NefData2[compList2[0]], options, cItem=cItem, nefList=nefList)

    return nefList


#=========================================================================================
# Test_Compare_Files
#=========================================================================================

def batchCompareNefFiles(inDir1, inDir2, outDir, options):
    """Batch compare the Nef files common to the two directories
    For each file found, write the compare log to the corresponding .txt file

    :param inDir1:
    :param inDir2:
    :param outDir:
    :param options: nameSpace holding the commandLineArguments
    """
    inFileList = [f for f in listdir(inDir1) if isfile(join(inDir1, f)) and f[-4:] == '.nef']
    outFileList = [f for f in listdir(inDir2) if isfile(join(inDir2, f)) and f[-4:] == '.nef']

    if not options.screen:
        if options.createDirs is True and not os.path.exists(outDir):
            os.mkdir(outDir)
        if not (os.path.exists(outDir) and os.path.isdir(outDir)):
            print('Error: No such directory:', str(outDir))
            return

    commonFiles = set(listdir(inDir1)) & set(listdir(inDir2))
    if not commonFiles:
        # if no files found then write message to the screen or log.tx in the out folder
        if options.screen is True:
            # strip the .nef from the end
            print('inDir1: %s' % str(inDir1))
            print('inDir2: %s' % str(inDir2))
            print('No common files found')
        else:
            outFileName = join(outDir, 'log.txt')
            if options.overwriteExisting is False:
                with safeOpen(outFileName, 'w') as (outLog, safeFileName):
                    outLog.write('inDir1: %s\n' % str(inDir1))
                    outLog.write('inDir2: %s\n' % str(inDir2))
                    outLog.write('No common files found')
            else:
                with open(outFileName, 'w') as outLog:
                    sys.stdout = outLog
                    outLog.write('inDir1: %s\n' % str(inDir1))
                    outLog.write('inDir2: %s\n' % str(inDir2))
                    outLog.write('No common files found')
        return

    for fl in inFileList:
        if fl in outFileList:

            if options.screen is True:

                # strip the .nef from the end
                outFileName = join(outDir, fl[:-4] + '.txt')
                print('Batch processing %s > %s' % (fl, outFileName))

                nefList = compareNefFiles(join(inDir1, fl), join(inDir2, fl), options)
                printCompareList(nefList, join(inDir1, fl), join(inDir2, fl))

            else:
                # strip the .nef from the end
                outFileName = join(outDir, fl[:-4] + '.txt')
                # keep the old output stream
                stdOriginal = sys.stdout

                if options.overwriteExisting is False:

                    with safeOpen(outFileName, 'w') as (outLog, safeFileName):
                        sys.stdout = outLog
                        nefList = compareNefFiles(join(inDir1, fl), join(inDir2, fl), options)

                        print('Batch processing %s > %s' % (fl, os.path.basename(safeFileName)))
                        print(join(inDir1, fl))
                        print(join(inDir2, fl))
                        printCompareList(nefList, join(inDir1, fl), join(inDir2, fl))

                else:
                    with open(outFileName, 'w') as outLog:
                        sys.stdout = outLog
                        nefList = compareNefFiles(join(inDir1, fl), join(inDir2, fl), options)

                        print('Batch processing %s > %s' % (fl, outFileName))
                        print(join(inDir1, fl))
                        print(join(inDir2, fl))
                        printCompareList(nefList, join(inDir1, fl), join(inDir2, fl))
                sys.stdout = stdOriginal


#=========================================================================================
# Test_Compare_Files
#=========================================================================================

class Test_Compare_Files(unittest.TestCase):
    """Test the comparison of nef files and print the results
    """

    def test_Compare_Files(self):
        """Load two files and compare
        """
        # define arguments to simulate command line
        parser = defineArguments()
        options = parser.parse_args([])

        # set the two files to compare
        inFile1 = str(os.path.join('.', 'nef', 'testdata', 'Commented_Example.nef'))
        inFile2 = str(os.path.join('.', 'nef', 'testdata', 'Commented_Example_Change.nef'))

        print('\nTEST COMPARISON')
        print('   file1 = ' + inFile1)
        print('   file2 = ' + inFile2)
        print('Loading...')

        # load and output results
        nefList = compareNefFiles(inFile1, inFile2, options)
        printCompareList(nefList, inFile1, inFile2)

    def test_Compare_BatchFiles(self):
        """Compare the Nef files in two directories
        """
        # define arguments to simulate command line
        parser = defineArguments()
        options = parser.parse_args([])
        options.createDirs = True
        options.overwriteExisting = False

        inDir1 = os.path.join('testdata', 'testinfolder1')
        inDir2 = os.path.join('testdata', 'testinfolder2')
        outDir = os.path.join('testdata', 'testoutfolder')

        batchCompareNefFiles(inDir1, inDir2, outDir, options)


#=========================================================================================
# __main__
#=========================================================================================

if __name__ == '__main__':
    helpText = """
  Command Line Usage:
    compareNef for execution from the command line with a suitable script
    An example can be found in AnalysisV3/bin/compareNef:
  
        #!/usr/bin/env sh
        export CCPNMR_TOP_DIR="$(dirname $(cd $(dirname "$0"); pwd))"
        export ANACONDA3=${CCPNMR_TOP_DIR}/miniconda
        export PYTHONPATH=${CCPNMR_TOP_DIR}/src/python:${CCPNMR_TOP_DIR}/src/c
        ${ANACONDA3}/bin/python ${CCPNMR_TOP_DIR}/src/python/ccpn/util/nef/compareNef.py $*
  
    Usage:  compareNef [options]
                                                  
    optional arguments:
      -h, --help              show this help message
      -H, --Help              Show detailed help

      -i, --ignoreblockname   Ignore the blockname when comparing two Nef files
                              May be required when converting Nef files through
                              different applications
                              May be used with -f and -b
                              
      -f inFile1 inFile2, --file inFile1 inFile2
                              Compare two Nef files and print the results to the
                              screen
                              
      -b inDir1 inDir2 outDir, --block inDir1 inDir2 outDir
                              compare Nef files common to directories
                              inDir1 and inDir2. Write output *.txt for each
                              file into the outDir directory.

      -s, --screen            Output batch processing to screen, default is to .txt files
                              may be used with -b

      -o, --overwrite         Overwrite existing .txt files. If false then files are
                              appended with '(n)' before the extension, where n is 
                              the next available number.
                              
      -c, --create            Automatically create directories as required. 
      
      -I, --ignorecase        Ignore case when comparing items.
                          
    Searches through all objects: dataExtents, dataBlocks, saveFrames and Loops within the files.
    Comparisons are made for all data structures that have the same name.
    Differences for items within a column are listed in the form:
      dataExtent:dataBlock:saveFrame:Loop:  <Column>: columnName  <rowIndex>: row  -->  value1 != value2
      
    dataExtents, dataBlocks, saveFrames, Loops, columns present only in one file are listed.
      dataExtent:dataBlock:saveFrame:Loop: contains --> parameter1
                                                        parameter2
                                                        ...
                                                        parameterN
    
    Please see the contents of CompareNef.py for functions available to python.
  """

    parser = defineArguments()
    commandLineArguments = parser.parse_args()

    if commandLineArguments.help:
        print(helpText)
    else:
        if commandLineArguments.inFiles is not None:
            inFile0 = commandLineArguments.inFiles[0]
            inFile1 = commandLineArguments.inFiles[1]

            print()
            print('Loading Nef Files...')
            nefList = compareNefFiles(inFile0, inFile1, commandLineArguments)
            printCompareList(nefList, inFile0, inFile1)

        elif commandLineArguments.blockDirs is not None:
            inDir0 = commandLineArguments.blockDirs[0]
            inDir1 = commandLineArguments.blockDirs[1]
            outDir = commandLineArguments.blockDirs[2]

            batchCompareNefFiles(inDir0, inDir1, outDir, commandLineArguments)

        else:
            print('Incorrect arguments, use compareNef -h')

"""
This documentation has now been moved to file ../CompareNef.py
"""
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
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import os
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
    import_parents(level=2)         # 2 because need to import with 2 dots below
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from .Paths import TEST_FILE_PATH
from ..CompareNef import compareNefFiles, printCompareList, defineArguments


if __name__ == '__main__':
    """
    Load two files and compare
    """

    # define arguments to simulate command line
    parser = defineArguments()

    # set ignoreBlockName flag to True
    commandLineArguments = parser.parse_args()

    inFile1 = os.path.join(TEST_FILE_PATH, 'Commented_Example.nef')
    inFile2 = os.path.join(TEST_FILE_PATH, 'Commented_Example_Change.nef')

    print('\nTEST COMPARISON')
    print('   file1 = ' + inFile1)
    print('   file2 = ' + inFile2)
    print('Loading...')
    nefList = compareNefFiles(inFile1, inFile2, commandLineArguments)
    printCompareList(nefList, inFile1, inFile2)

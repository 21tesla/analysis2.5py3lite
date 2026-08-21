"""
Unit test execute as:
python -u $CINGROOT/python/cing/Scripts/test/test_linkdataNRGforGeerten.py
"""
import unittest
import os
from unittest import TestCase

from cing import cingDirTmp  #@UnusedImport
from cing.Libs.NTutils import *  #@UnusedWildImport
from cing.Scripts.linkdataNRGforGeerten import *  #@UnusedWildImport


class AllChecks(TestCase):

    def _test_linkdataNRGforGeerten(self):
        cingDirTmpTest = os.path.join( cingDirTmp, getCallerName() )
        os.makedirs( cingDirTmpTest , exist_ok=True)
        self.assertFalse(os.chdir(cingDirTmpTest), msg =
            "Failed to change to test directory for files: " + cingDirTmpTest)

        if os.path.exists(outputFnLinkData):
            os.remove(outputFnLinkData)
        self.assertFalse(linkdataNRG())
        self.assertTrue(os.path.exists(outputFnLinkData))

if __name__ == "__main__":
    cing.verbosity = verbosityDebug
    unittest.main()

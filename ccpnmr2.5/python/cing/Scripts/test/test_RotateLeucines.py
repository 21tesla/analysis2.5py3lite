"""
Unit test execute as:
python -u $CINGROOT/python/cing/Scripts/test/test_RotateLeucines.py
"""
import unittest
import os
from unittest import TestCase

from unittest import SkipTest

from cing import cingDirTestsData  #@UnusedImport
from cing.Libs.NTutils import *  #@UnusedWildImport
from cing.PluginCode.required.reqYasara import YASARA_STR

# Import using optional plugins.
try:
    # A bit redundant with above line.
    from cing.Scripts.rotateLeucines import *  #@UnusedWildImport Relies on Yasara as well.
except ImportWarning as extraInfo: # Disable after done debugging; can't use nTdebug yet.
    print("Got ImportWarning %-10s Skipping unit check %s." % ( YASARA_STR, getCallerFileName() ))
    raise SkipTest(YASARA_STR)
# end try

class AllChecks(TestCase):
    def _test_rotateLeucinesInYasara(self):
        '''
        This unit test is by default disabled because we haven't figured out yet how to disable the output from
        Yasara yet. 
        '''
        cingDirTmpTest = os.path.join( cingDirTmp, getCallerName() )
        os.makedirs( cingDirTmpTest , exist_ok=True)
        self.assertFalse(os.chdir(cingDirTmpTest), msg =
            "Failed to change to test directory for files: " + cingDirTmpTest)

        entryId = '1brv'
#        entryId = 'H2_2Ca_64_100'
        inputArchiveDir = os.path.join(cingDirTestsData, "cing")
        self.assertFalse( runRotateLeucines(cingDirTmpTest, inputArchiveDir, entryId, useAll = True))
    # end def
# end class

if __name__ == "__main__":
    cing.verbosity = verbosityDebug
    unittest.main()
# end if

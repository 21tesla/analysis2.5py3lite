"""python $CINGROOT/python/cing/PluginCode/test/test_shiftx.py
"""
import unittest
import os
from unittest import TestCase

from unittest import SkipTest

from cing import cingDirTestsData, cingDirTmp
from cing.core.classes import Project
from cing.core.constants import *  #@UnusedWildImport
from cing.Libs.NTutils import *  #@UnusedWildImport
from cing.PluginCode.required.reqShiftx import SHIFTX_STR

# Import using optional plugins.
try:
    pass
except ImportWarning as extraInfo: # Disable after done debugging; can't use nTdebug yet.
    print("Got ImportWarning %-10s Skipping unit check %s." % ( SHIFTX_STR, getCallerFileName() ))
    raise SkipTest(SHIFTX_STR)
# end try

class AllChecks(TestCase):
    def test_shiftx(self):
#        entryId = "1brv" # Small much studied PDB NMR entry
#        entryId = "2hgh_1model"  RNA-protein complex.
        entryId = "1brv"
#        entryId = "1tgq_1model" # withdrawn entry
        cingDirTmpTest = os.path.join( cingDirTmp, getCallerName() )
        os.makedirs( cingDirTmpTest , exist_ok=True)
        self.assertFalse(os.chdir(cingDirTmpTest), msg =
            "Failed to change to test directory for files: " + cingDirTmpTest)

        project = Project( entryId )
        self.assertFalse( project.removeFromDisk())
        project = Project.open( entryId, status='new' )
        cyanaFile = os.path.join(cingDirTestsData, "cyana", entryId + ".cyana.tgz")
        self.assertTrue(project.initCyana(cyanaFolder = cyanaFile))
        project.runShiftx()

if __name__ == "__main__":
    cing.verbosity = verbosityDebug
    unittest.main()

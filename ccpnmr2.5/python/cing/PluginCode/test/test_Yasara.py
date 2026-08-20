"""
Unit test
python $CINGROOT/python/cing/PluginCode/test/test_Yasara.py

Disabled because shell is hard to test; see $CINGROOT/python/cing/Scripts/test/test_RotateLeucines.py
"""
import unittest
from unittest import TestCase

from nose.plugins.skip import SkipTest

from cing import cingDirTmp
from cing.core.classes import Project
from cing.core.constants import *  #@UnusedWildImport
from cing.Libs.NTutils import *  #@UnusedWildImport
from cing.PluginCode.required.reqYasara import YASARA_STR

# Import using optional plugins.
try:
    import yasara  #@UnresolvedImport
except ImportWarning as extraInfo: # Disable after done debugging; can't use nTdebug yet.
    print("Got ImportWarning %-10s Skipping unit check %s." % ( YASARA_STR, getCallerFileName() ))
    raise SkipTest(YASARA_STR)
# end try

class AllChecks(TestCase):

    def _test_Yasara(self):
        'Test is incomplete; do NOT use.'
        entryId = "testYasara"
        cingDirTmpTest = os.path.join( cingDirTmp, getCallerName() )
        mkdirs( cingDirTmpTest )
        self.assertFalse(os.chdir(cingDirTmpTest), msg =
            "Failed to change to test directory for files: " + cingDirTmpTest)
        project = Project( entryId )
#            project = Project.open( entryId, status='old' )
        project.yasaraShell()
        yasara.Exit() # FAILS to ext really.
    # end def
# end class

if __name__ == "__main__":
    cing.verbosity = verbosityDebug
    unittest.main()

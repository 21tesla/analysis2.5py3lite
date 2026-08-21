"""
Unit test execute as:
python $CINGROOT/python/cing/Database/Scripts/test/test_createStarUserLib.py
"""
import unittest
import os
from unittest import TestCase

import cing
from cing import cingDirTmp
from cing.Database.Scripts.createStarUserLib import createStarUserLib
from cing.Libs.NTutils import *  #@UnusedWildImport


class AllChecks(TestCase):

    cingDirTmpTest = os.path.join( cingDirTmp, 'test_createStarUserLib' )
    os.makedirs( cingDirTmpTest , exist_ok=True)
    os.chdir(cingDirTmpTest)

    def _test_createSimpleFastProject(self):
        createStarUserLib()
    # end def
# end class

if __name__ == "__main__":
    cing.verbosity = cing.verbosityDebug
    unittest.main()

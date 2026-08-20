"""
Unit test execute as:
python $CINGROOT/python/cing/PluginCode/test/test_x3dna.py

Open the ???.r3d files in pymol or so; they're nice.
"""
import shutil
import unittest
from unittest import TestCase

from nose.plugins.skip import SkipTest

from cing import cingDirTestsData, cingDirTmp, osType
from cing.core.classes import Project
from cing.Libs.NTutils import *  #@UnusedWildImport
from cing.PluginCode.required.reqCcpn import CCPN_STR
from cing.PluginCode.required.reqX3dna import X3DNA_STR

# Import using optional plugins.
try:
    pass
except ImportWarning as extraInfo: # Disable after done debugging; can't use nTdebug yet.
    print("Got ImportWarning %-10s Skipping unit check %s." % ( CCPN_STR, getCallerFileName() ))
    raise SkipTest(CCPN_STR)
# end try

class AllChecks(TestCase):

    def test_x3dna(self):
        if osType != OS_TYPE_MAC: # only installed for mac os currently.
            return
        try:
            from cing.PluginCode.x3dna import createHtmlX3dna
        except:
            nTerror("Failed to import x3dna on a mac osx.")
            return
        entryList = "1b4y".split() # triple helix but it only gets analyzed to a double helix.
#        entryList = "1cjg".split()
#        entryList = "2hgh".split()
#        entryList = "1a4d".split()
#        entryList = "1brv".split() # protein only entry
#        entryList = ["SRYBDNA"]

        useNrgArchive = False
        showValues = True
        cingDirTmpTest = os.path.join( cingDirTmp, getCallerName() )
        mkdirs( cingDirTmpTest )
        self.assertFalse(os.chdir(cingDirTmpTest), msg =
            "Failed to change to test directory for files: " + cingDirTmpTest)

        for entryId in entryList:
            project = Project.open(entryId, status = 'new')
            self.assertTrue(project, 'Failed opening project: ' + entryId)

            if useNrgArchive: # default is False
                inputArchiveDir = os.path.join('/Library/WebServer/Documents/NRG-CING/recoordSync', entryId)
            else:
                inputArchiveDir = os.path.join(cingDirTestsData, "ccpn")

            ccpnFile = os.path.join(inputArchiveDir, entryId + ".tgz")
            self.assertTrue(project.initCcpn(ccpnFolder = ccpnFile))
            project.save()
            self.assertTrue(project.runX3dna())
            project.save()
            self.assertFalse(createHtmlX3dna(project))

            if showValues:
                for coplanar in project.coplanars[0]:
#                    nTdebug(coplanar.format())
                    nTdebug('%r' % coplanar[X3DNA_STR])

            # Do not leave the old CCPN directory laying around since it might get added to by another test.
            if os.path.exists(entryId):
                self.assertFalse(shutil.rmtree(entryId))

if __name__ == "__main__":
    cing.verbosity = verbosityDebug
    unittest.main()


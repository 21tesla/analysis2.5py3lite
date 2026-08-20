import unittest
from unittest import TestCase

from cing import cingDirTmp, cingRoot
from cing.Libs.disk import globLast, tail
from cing.Libs.NTutils import *  #@UnusedWildImport


class AllChecks(TestCase):
    # important to switch to temp space before starting to generate files for the project.
    cingDirTmpTest = os.path.join( cingDirTmp, 'test_disk' )
    mkdirs( cingDirTmpTest )
    os.chdir(cingDirTmpTest)

    def testDisk(self):
        doneFileName = "DONE"
        f = open(doneFileName,"w")
        for i in range(10):
            f.write("Line %d\n" % i)
        f.close()
        f2 = open(doneFileName)
        lastLineList = tail(f2,1)
        lastLine = lastLineList[0]
        self.assertEqual( "Line 9", lastLine )
        self.assertEqual( "['Line 9']", repr(lastLineList) ) # not necessary a test.

    def testGlobLast(self):
        globPattern = os.path.join(cingRoot, '*.txt')
        lastFile = globLast(globPattern)
        nTdebug('lastFile: %s' % lastFile)
        d, _basename, extension = nTpath(lastFile)
        self.assertTrue(lastFile)
        self.assertEqual(d, cingRoot)
        self.assertEqual(extension, '.txt')

        globPattern = os.path.join(cingRoot, '*.xyz')
        lastFile = globLast(globPattern)
        nTdebug('lastFile 2: %s' % lastFile)
        self.assertFalse(lastFile)

if __name__ == "__main__":
    cing.verbosity = verbosityDebug
    unittest.main()

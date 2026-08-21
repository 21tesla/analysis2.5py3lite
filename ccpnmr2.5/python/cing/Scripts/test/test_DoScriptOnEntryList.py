"""
Unit test execute as:
python $CINGROOT/python/cing/Scripts/test/test_DoScriptOnEntryList.py
"""
import unittest
import os
from unittest import TestCase

from cing import cingDirScripts, cingDirTmp
from cing.Libs.NTutils import *  #@UnusedWildImport
from cing.Scripts.doScriptOnEntryList import doFunctionOnEntryList, doScriptOnEntryList
from cing.Scripts.validateEntry import ARCHIVE_TYPE_BY_CH23, PROJECT_TYPE_PDB


def sleepy(sleepTime, bogusArgumentList = [] ):
    nTdebug("Will sleep for %s ignoring bogusArgumentList: %s" % (sleepTime, str(bogusArgumentList)))
    time.sleep(float(sleepTime))
# end def

class AllChecks(TestCase):

    def test_DoScriptOnEntryList(self):

        cingDirTmpTest = os.path.join( cingDirTmp, getCallerName() )
        os.makedirs( cingDirTmpTest , exist_ok=True)
        self.assertFalse(os.chdir(cingDirTmpTest), msg =
            "Failed to change to test directory for files: " + cingDirTmpTest)
        entryListFileName = "entry_list_todo.csv"
        entry_list_todo = [ 0,1,2,3,4,5,6,7,8,9 ]
        writeTextToFile(entryListFileName, toCsv(entry_list_todo))

        pythonScriptFileName = os.path.join(cingDirScripts, 'doNothing.py')
        extraArgList = ('.', '.', '.', '.', ARCHIVE_TYPE_BY_CH23, PROJECT_TYPE_PDB)

        self.assertFalse(
            doScriptOnEntryList(pythonScriptFileName,
                            entryListFileName,
                            '.',
                            processes_max = 8,
                            delay_between_submitting_jobs = 5,
                            max_time_to_wait = 20,
                            start_entry_id = 0,
                            max_entries_todo = 1,
                            extraArgList = extraArgList,
                            shuffleBeforeSelecting = True ))
    # end def

    def test_DoFunctionOnEntryList(self):
        cingDirTmpTest = os.path.join( cingDirTmp, getCallerName() )
        os.makedirs( cingDirTmpTest , exist_ok=True)
        self.assertFalse(os.chdir(cingDirTmpTest), msg =
            "Failed to change to test directory for files: " + cingDirTmpTest)
        entryListFileName = 'entryListFileName.csv'
        writeTextToFile(entryListFileName, '\n'.join('0.1 0.2'.split()))
        doFunctionOnEntryList(sleepy, entryListFileName)
    # end def

if __name__ == "__main__":
    cing.verbosity = verbosityDebug
    unittest.main()

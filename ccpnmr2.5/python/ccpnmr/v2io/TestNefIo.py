"""
Module Documentation here
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = ""
__credits__ = ""
__licence__ = ("")
__reference__ = ("")
#=========================================================================================
# Last code modification:
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified$"
__version__ = "$Revision$"
#=========================================================================================
# Created:
#=========================================================================================
__author__ = "$Author: Ed Brooksbank $"
__date__ = "$Date$"
#=========================================================================================
# Start of code
#=========================================================================================

import memops.gui.MessageReporter
import memops.universal.MessageReporter


memops.gui.MessageReporter.showWarning = memops.universal.MessageReporter.showWarning
memops.gui.MessageReporter.showOkCancel = memops.universal.MessageReporter.showOkCancel
memops.gui.MessageReporter.showYesNo = memops.universal.MessageReporter.showYesNo

import re
import os
import unittest
import decorator
from unittest import TestCase, skip
from ccpnmr.v2io.NefIo import loadNefFile
from ccpnmr.nef.NefImporter import NefImporter


CCPNTestPath = '/Users/ejb66/PycharmProjects/Git/NEF/data_1_1/'
CCPNTestFiles = ('CCPN_1nk2_docr.nef',
                 'CCPN_2kko_docr.nef',
                 'CCPN_2mqq_docr.nef',
                 'CCPN_2mtv_docr.nef',
                 'CCPN_Commented_Example.nef',
                 'CCPN_H1GI_clean_extended.nef',
                 'CCPN_H1GI_clean.nef',
                 'CCPN_Sec5Part3.nef',
                 'CCPN_XPLOR_test1.nef',
                 'CSROSETTA_test1.nef',
                 'XPLOR_test1.nef',
                 )

CASDTestFolder = 'CCPN_CASD'
CASDTestPath = os.path.join(CCPNTestPath, CASDTestFolder)
CASDTestFiles = ('2l9r_docr.nef',
                 '2la6_docr.nef',
                 '2lah_docr.nef',
                 '2lci_docr.nef',
                 '2ln3_docr.nef',
                 '2loj_docr.nef',
                 '2ltl_docr.nef',
                 '2ltm_docr.nef',
                 '2m2e_docr.nef',
                 '2m5o_docr.nef',
                 )

PDBStatTestFolder = 'PDBStat'
PDBStatTestPath = os.path.join(CCPNTestPath, PDBStatTestFolder)
PDBStatTestFiles = ('BeR31/BeR31_2k2e.nef',
                    'CsR4/CsR4_2jr2.nef',
                    'CtR107/CtR107_2kcu.nef',
                    'CtR148A/CtR148A_2ko1.nef',
                    'DhR8C/DhR8C_2kyi.nef',
                    'DhR29B/DhR29B_2kpu.nef',
                    'HR41/HR41_2ko7.nef',
                    'L22_T12/L22_rt.nef',
                    'MbR242E/MbR242E_2kko.nef',
                    'MiR12/MiR12_2luz.nef',
                    'SgR145/SgR145_2kw5.nef',
                    'SoR77/SoR77_2juw.nef',
                    'SPR104/SPR104_2l3a.nef',
                    'SR10/SR10_2kzn.nef',
                    'WR73/WR73_2loy.nef',
                    'ZR18/ZR18_1pqx.nef',
                    )
NEFVALIDATIONFILE = '/Users/ejb66/PycharmProjects/Git/NEF/specification/mmcif_nef.dic'


def selectFile():
    """A decorator to wrap a state change event with a verify function
    """

    @decorator.decorator
    def theDecorator(*args, **kwds):
        func = args[0]
        args = args[1:]  # Optional 'self' is now args[0]
        self = args[0]

        # get filename from the method name - test_load<testName><testNum>
        testID = re.findall(r'test_load(\w+?)(\d+)', str(func.__name__))
        testList = globals()[testID[0][0]+'TestFiles']
        testPath = globals()[testID[0][0]+'TestPath']
        testFile = testList[int(testID[0][1])]

        filePath = os.path.join(testPath, testFile)
        print('Running Test:', func.__name__, testFile)
        if os.path.isfile(filePath):

            test = NefImporter(errorLogging='silent')
            test.loadFile(filePath)
            test.loadValidateDictionary(NEFVALIDATIONFILE)
            valid = test.isValid
            if not valid:
                print('~~~~~~~~~~~~~~~~~~~~~~\nValid Nef:', valid)
                print('Error in Nef:', test.validErrorLog)
                print('Error in Nef:')
                for k, v in test.validErrorLog.items():
                    print('  >>>', k)
                    for err in v:
                        print('  >>>        ', err)
                print('~~~~~~~~~~~~~~~~~~~~~~')
            self.assertTrue(valid, 'Error validating nef file.')

            loadNefFile(path=filePath, overwriteExisting=True)

        return

    return theDecorator


class AllCCPNChecks(TestCase):

    def setUp(self):
        pass

    @selectFile()
    def test_loadCCPN0(self):
        """CCPN_1nk2_docr.nef"""
        pass

    @selectFile()
    def test_loadCCPN1(self):
        """CCPN_2kko_docr.nef"""
        pass

    @selectFile()
    def test_loadCCPN2(self):
        """CCPN_2mqq_docr.nef"""
        pass

    @selectFile()
    def test_loadCCPN3(self):
        """CCPN_2mtv_docr.nef"""
        pass

    @selectFile()
    def test_loadCCPN4(self):
        """CCPN_Commented_Example.nef"""
        pass

    @selectFile()
    def test_loadCCPN5(self):
        """CCPN_H1GI_clean_extended.nef"""
        pass

    @selectFile()
    def test_loadCCPN6(self):
        """CCPN_H1GI_clean.nef"""
        pass

    @selectFile()
    def test_loadCCPN7(self):
        """CCPN_Sec5Part3.nef"""
        pass

    @selectFile()
    def test_loadCCPN8(self):
        """CCPN_XPLOR_test1.nef"""
        pass

    @selectFile()
    def test_loadCCPN9(self):
        """CSROSETTA_test1.nef"""
        pass

    @selectFile()
    def test_loadCCPN10(self):
        """XPLOR_test1.nef"""
        pass

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class AllCASDChecks(TestCase):

    def setUp(self):
        pass

    @selectFile()
    def test_loadCASD0(self):
        """2l9r_docr.nef"""
        pass

    @selectFile()
    def test_loadCASD1(self):
        """2la6_docr.nef"""
        pass

    @selectFile()
    def test_loadCASD2(self):
        """2lah_docr.nef"""
        pass

    @selectFile()
    def test_loadCASD3(self):
        """2lci_docr.nef"""
        pass

    @selectFile()
    def test_loadCASD4(self):
        """2ln3_docr.nef"""
        pass

    @selectFile()
    def test_loadCASD5(self):
        """2loj_docr.nef"""
        pass

    @selectFile()
    def test_loadCASD6(self):
        """2ltl_docr.nef"""
        pass

    @selectFile()
    def test_loadCASD7(self):
        """2ltm_docr.nef"""
        pass

    @selectFile()
    def test_loadCASD8(self):
        """2m2e_docr.nef"""
        pass

    @selectFile()
    def test_loadCASD9(self):
        """2m5o_docr.nef"""
        pass

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class AllPDBChecks(TestCase):

    def setUp(self):
        pass

    @selectFile()
    def test_loadPDBStat0(self):
        """BeR31/BeR31_2k2e.nef"""
        pass

    @selectFile()
    def test_loadPDBStat1(self):
        """CsR4/CsR4_2jr2.nef"""
        pass

    @selectFile()
    def test_loadPDBStat2(self):
        """CtR107/CtR107_2kcu.nef"""
        pass

    @selectFile()
    def test_loadPDBStat3(self):
        """CtR148A/CtR148A_2ko1.nef"""
        pass

    @selectFile()
    def test_loadPDBStat4(self):
        """DhR8C/DhR8C_2kyi.nef"""
        pass

    @selectFile()
    def test_loadPDBStat5(self):
        """DhR29B/DhR29B_2kpu.nef"""
        pass

    @selectFile()
    def test_loadPDBStat6(self):
        """HR41/HR41_2ko7.nef"""
        pass

    @selectFile()
    def test_loadPDBStat7(self):
        """L22_T12/L22_rt.nef"""
        pass

    @selectFile()
    def test_loadPDBStat8(self):
        """MbR242E/MbR242E_2kko.nef"""
        pass

    @selectFile()
    def test_loadPDBStat9(self):
        """MiR12/MiR12_2luz.nef"""
        pass

    @selectFile()
    def test_loadPDBStat10(self):
        """SgR145/SgR145_2kw5.nef"""
        pass

    @selectFile()
    def test_loadPDBStat11(self):
        """SoR77/SoR77_2juw.nef"""
        pass

    @selectFile()
    def test_loadPDBStat12(self):
        """SPR104/SPR104_2l3a.nef"""
        pass

    @selectFile()
    def test_loadPDBStat13(self):
        """SR10/SR10_2kzn.nef"""
        pass

    @selectFile()
    def test_loadPDBStat14(self):
        """WR73/WR73_2loy.nef"""
        pass

    @selectFile()
    def test_loadPDBStat15(self):
        """ZR18/ZR18_1pqx.nef"""
        pass


if __name__ == "__main__":
    unittest.main()

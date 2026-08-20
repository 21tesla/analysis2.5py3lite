"""
Unit test execute as:
python -u $CINGROOT/python/cing/Libs/test/test_svd.py
"""

import unittest
from unittest import TestCase

from cing.Libs.NTutils import *  #@UnusedWildImport
from cing.Libs.svd import *  #@UnusedWildImport


class AllChecks(TestCase):
    pass
# end class

if __name__ == "__main__":
    cing.verbosity = cing.verbosityDebug
    unittest.main()

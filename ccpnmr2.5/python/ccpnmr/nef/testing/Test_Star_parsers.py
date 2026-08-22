"""Module Documentation here"""

# =========================================================================================
# Licence, Reference and Credits
# =========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2020"
__credits__ = "Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister"
__licence__ = "CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license"
__reference__ = (
    "Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
    "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
    "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y",
)
# =========================================================================================
# Last code modification
# =========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2020-01-14 11:49:36 +0000 (Tue, January 14, 2020) $"
__version__ = "$Revision: 3.0.0 $"
# =========================================================================================
# Created
# =========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
# =========================================================================================
# Start of code
# =========================================================================================

import os
import sys
import time
import unittest

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this is a fix to get the import to work when running as a standalone
# when importing into your own code, with PYTHON_PATH defined it can be safely removed


def import_parents(level=1):
    global __package__

    import importlib
    from os import path

    # pathlib does all this a lot nicer, but don't think it's in python2.7
    top = parent = path.dirname(path.abspath(__file__))
    package = []
    for t in range(level):
        package.insert(0, os.path.basename(top))
        top = path.dirname(top)

    sys.path.append(str(top))
    try:
        sys.path.remove(str(parent))
    except ValueError:  # already removed
        pass

    __package__ = str(".".join(package))
    importlib.import_module(__package__)


if __name__ == "__main__" and __package__ is None:
    import_parents(level=2)  # 2 because need to import with 2 dots below
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


from .. import GenericStarParser, StarIo
from .Paths import TEST_FILE_PATH


def _requireFile(path):
    """Resolve *path* (relative paths against TEST_FILE_PATH); skip when absent.

    The parser-sample tests below reference real CCPN/BMRB sample files, but the
    distribution only bundles some of them (see nef/testdata/).  A test whose
    required sample is absent is skipped with a clear reason instead of failing,
    so the suite stays green from source, wheel, sdist and CI alike.  Use the
    CCP_TESTDATA env var to point at a full local dataset (see Paths.py).
    """
    usePath = path if path.startswith("/") else os.path.join(TEST_FILE_PATH, path)
    if not os.path.exists(usePath):
        raise unittest.SkipTest(
            "sample data %r not bundled in this distribution (looked in %r); "
            "set CCP_TESTDATA to a directory that contains it" % (path, TEST_FILE_PATH)
        )
    return usePath


def _loadGeneralFile(path):
    usePath = _requireFile(path)
    t0 = time.time()
    entry = GenericStarParser.parseFile(usePath)  # 'lenient')
    print("Parsing time %s for %s" % (time.time() - t0, path))
    return entry


def _loadNmrStarFile(path):
    usePath = _requireFile(path)
    t0 = time.time()
    entry = StarIo.parseNmrStarFile(usePath)  # 'lenient')
    print("Parsing time %s for %s" % (time.time() - t0, path))
    return entry


def _loadNefFile(path):
    usePath = _requireFile(path)
    t0 = time.time()
    entry = StarIo.parseNefFile(usePath)  # 'lenient')
    print("Parsing time %s for %s" % (time.time() - t0, path))
    return entry


def _printContents(object, indent=0):
    indentStep = 4
    indent += indentStep
    print(" " * indent, object, "Tags: %s" % len(object))
    for tag, obj in object.items():
        if isinstance(obj, GenericStarParser.Loop):
            if tag == obj.columns[0]:
                print(" " * (indent + indentStep), obj, "Columns: %s" % len(obj.columns))
        elif not isinstance(obj, str):
            _printContents(obj, indent)


def test_nmrstar_4267():
    print("\n\n", "# nmrstar_4267", "#" * 60, "\n")
    _loadGeneralFile("4267_example.str")
    _loadNmrStarFile("4267_example.str")


def test_nef_commented_example():
    _loadGeneralFile("CCPN_Commented_Example.nef")
    _loadNefFile("CCPN_Commented_Example.nef")


def test_nef_2l9r_Paris_155():
    print("\n\n", "# Paris_155_nef", "#" * 60, "\n")
    _loadGeneralFile("CCPN_2l9r_Paris_155.nef")
    _loadNefFile("CCPN_2l9r_Paris_155.nef")


def test_nef_1lci_Piscataway_179():
    _loadGeneralFile("CCPN_2lci_Piscataway_179.nef")
    _loadNefFile("CCPN_2lci_Piscataway_179.nef")


def test_nef_H1GI():
    _loadGeneralFile("CCPN_H1GI.nef")
    _loadNefFile("CCPN_H1GI.nef")


def test_mmcif_1bgl_1bgm():
    _loadGeneralFile("1bgl_1bgm.cif")


def test_dic_mmcif_nef():
    print("\n\n", "# nef_dic", "#" * 60, "\n")
    _loadGeneralFile("mmcif_nef.dic")


def test_dic_mmcif_nmr_star():
    _loadGeneralFile("mmcif_nmr-star.dic")


def test_dic_mmcif_std():
    print("\n\n", "# mmcif_dic", "#" * 60, "\n")
    _loadGeneralFile("mmcif_std.dic")


def test_dic_mmcif_pdbx_v40():
    _loadGeneralFile("mmcif_pdbx_v40.dic")


def test_parse_bundled_nef_samples():
    """Parse every NEF sample that ships in the distribution (nef/testdata/).

    Exercises the NEF parser against the *shipped* samples, so this passes from
    source, wheel and sdist alike -- no external BMRB dataset required.  Unlike
    the external-sample tests above (skipped when their file is not bundled),
    this always has real data.  Skips only if the bundled sample directory is
    absent entirely.
    """
    if not os.path.isdir(TEST_FILE_PATH):
        raise unittest.SkipTest("bundled NEF sample directory not found: %s" % TEST_FILE_PATH)
    samples = sorted(fn for fn in os.listdir(TEST_FILE_PATH) if fn.lower().endswith(".nef"))
    if not samples:
        raise unittest.SkipTest("no bundled .nef samples found in %s" % TEST_FILE_PATH)
    parsed = []
    for fn in samples:
        entry = GenericStarParser.parseFile(os.path.join(TEST_FILE_PATH, fn))
        assert entry is not None, "parseFile returned None for %s" % fn
        # A NEF/mmcif sample parses into an object that exposes a tag mapping.
        assert getattr(entry, "items", None) is not None, "%s parsed with no tag mapping" % fn
        parsed.append(fn)
    assert parsed, "no .nef samples were parsed"
    print("parsed %d bundled NEF sample(s): %s" % (len(parsed), ", ".join(parsed)))


if __name__ == "__main__":
    # load and run a test cases
    test_nef_commented_example()
    test_parse_bundled_nef_samples()

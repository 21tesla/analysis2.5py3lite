"""
Module Documentation here
"""

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
__author__ = "$Author: Ed Brooksbank $"
__date__ = "$Date: 2020-01-13 17:46:24 +0000 (Mon, January 13, 2020) $"
# =========================================================================================
# Start of code
# =========================================================================================

import os

# ---------------------------------------------------------------------------
# Location of the NEF test-data files.
#
# The previous value anchored on the working directory:
#       os.path.join(os.path.dirname(os.getcwd()), "testdata")
# which only resolved when the test runner's cwd happened to be the *parent*
# of a directory literally named "testdata".  That made these tests data-gated
# for every consumer of the distribution -- and unrelated to where the bundled
# samples actually live.
#
# Anchor on the package instead: the samples ship with the code in
# .../ccpnmr/nef/testdata/, right next to this module (nef/testing/).  An
# explicit CCP_TESTDATA env var takes precedence (use it to point at a full
# local dataset, e.g. a BMRB checkout); the legacy cwd-based location is kept
# only as a last-resort fallback so no existing setup silently loses its data.
# ---------------------------------------------------------------------------
_PKG_TESTDATA = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "testdata"  # .../ccpnmr/nef/testdata
)
_LEGACY_TESTDATA = os.path.join(os.path.dirname(os.getcwd()), "testdata")


def _resolve_test_file_path():
    for cand in (os.environ.get("CCP_TESTDATA"), _PKG_TESTDATA, _LEGACY_TESTDATA):
        if cand and os.path.isdir(cand):
            return cand
    return _PKG_TESTDATA  # stable, well-defined default even when no dir exists


TEST_FILE_PATH = _resolve_test_file_path()

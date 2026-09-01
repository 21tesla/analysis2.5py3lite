"""Unit tests for the NMRdraw .tab import core (Add Peaks feature).

The pure file-parsing and assignment-string helpers of
``ccpnmr.analysis.core.NmrdrawImport`` are covered here headlessly; the
model-mutating import path (peaks, resonances, shift lists) is exercised
by running ``importTabPeaks`` against a project copy (see the session
E2E script), as the full model is heavier than this unit suite's norm.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ccpnmr.analysis.core.PeakListImport import (  # noqa: E402
    _assResidueCode,
    _isotopeElement,
    _rowIntensity,
    parseAssString,
    parseTabFile,
)

TAB_SAMPLE = """VARS   INDEX X_AXIS Y_AXIS X_PPM Y_PPM HEIGHT TYPE ASS CLUSTID MEMCNT
FORMAT %5d %9.3f %9.3f %8.3f %8.3f %e %d %s %4d %4d

NULLVALUE -666
NULLSTRING *

    1    27.858    17.204  10.213  130.039 +1.7e+12 1 W81-HE1    1    1
    2    33.619    22.174  10.152  129.107 -666 1 *    1    1
    3    64.731    68.316   9.824  120.455 +8.0e+11 1 Q15-HN    1    1
"""


def _write_tab(tmp_path, content):
    path = os.path.join(str(tmp_path), "test.tab")
    with open(path, "w") as fp:
        fp.write(content)
    return path


def test_parse_tab_file_data_rows(tmp_path):
    varNames, rows = parseTabFile(_write_tab(tmp_path, TAB_SAMPLE))
    assert varNames[0] == "INDEX"
    assert "X_PPM" in varNames and "ASS" in varNames
    assert len(rows) == 3
    assert rows[0]["INDEX"] == "1"
    assert float(rows[0]["X_PPM"]) == pytest.approx(10.213)
    assert rows[0]["ASS"] == "W81-HE1"
    assert float(rows[0]["HEIGHT"]) == pytest.approx(1.7e12)


def test_parse_tab_file_nulls(tmp_path):
    _, rows = parseTabFile(_write_tab(tmp_path, TAB_SAMPLE))
    # NULLVALUE -666 and NULLSTRING * become None
    assert rows[1]["HEIGHT"] is None
    assert rows[1]["ASS"] is None


def test_parse_tab_file_rejects_garbage(tmp_path):
    path = _write_tab(tmp_path, "not a tab file at all\n1 2 3\n")
    with pytest.raises(ValueError):
        parseTabFile(path)


def test_parse_tab_file_rejects_empty(tmp_path):
    path = _write_tab(tmp_path, "")
    with pytest.raises(ValueError):
        parseTabFile(path)


def test_parse_ass_basic():
    info = parseAssString("W81-HE1")
    assert info == {"chain": None, "res": "W", "num": 81, "atom": "HE1"}


def test_parse_ass_with_chain():
    info = parseAssString("A-W81-HE1")
    assert info is not None and info["chain"] == "A" and info["num"] == 81 and info["atom"] == "HE1"


def test_parse_ass_numbers_only_atom():
    info = parseAssString("Q15-HN")
    assert info["res"] == "Q" and info["atom"] == "HN"


def test_parse_ass_none_markers():
    assert parseAssString(None) is None
    assert parseAssString("") is None
    assert parseAssString("None") is None
    assert parseAssString("  ") is None


def test_parse_ass_unparseable():
    assert parseAssString("W81") is None  # missing atom
    assert parseAssString("??-99") is None  # no residue type
    assert parseAssString("10.213") is None  # a bare number


def test_ass_residue_codes():
    assert _assResidueCode("W") == "Trp"
    assert _assResidueCode("w") == "Trp"
    assert _assResidueCode("Q") == "Gln"
    assert _assResidueCode("TRP") == "Trp"
    assert _assResidueCode("GLN") == "Gln"
    assert _assResidueCode("ZQZ") is None  # unknown code -> no type constraint
    assert _assResidueCode("WL") is None  # two letters is not a code


def test_isotope_element():
    assert _isotopeElement("1H") == "H"
    assert _isotopeElement("15N") == "N"
    assert _isotopeElement("13C") == "C"
    assert _isotopeElement(None) is None
    assert _isotopeElement("unknown") is None


def test_row_intensity():
    row = {"HEIGHT": "123.4", "VOL": None, "WEIRD": "not-a-number", "NEG": "-666"}
    assert _rowIntensity(row, "HEIGHT") == pytest.approx(123.4)
    assert _rowIntensity(row, "VOL") is None
    assert _rowIntensity(row, "WEIRD") is None
    assert _rowIntensity(row, "MISSING") is None
    assert _rowIntensity(row, "NEG") == pytest.approx(-666.0)  # -666 is a legal int here


def _clean_project(nmr_project, spectrum):
    for peak_list in list(spectrum.peakLists):
        peak_list.delete()
    for resonance in list(nmr_project.resonances):
        resonance.delete()
    for resonance_group in list(nmr_project.resonanceGroups):
        resonance_group.delete()


def test_import_tab_peaks_integration():
    from memops.general import Io as memopsIo
    from ccpnmr.analysis.core.PeakListImport import importTabPeaks
    
    # Locate project directory relative to tests directory
    proj_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "addpeaks")
    tab_file = os.path.join(os.path.dirname(__file__), "..", "..", "..", "sswt_assigned.tab")
    
    root = memopsIo.loadProject(proj_dir, projectName="allpeaks2")
    nmr_project = root.currentNmrProject
    spectrum = nmr_project.findFirstExperiment(name="sswt").findFirstDataSource(name="sswt-298K-hsqc-1016")
    
    _clean_project(nmr_project, spectrum)
    report = importTabPeaks(root, tab_file, spectrum)
    assert report["error"] is None
    assert report["peaksAdded"] == 3
    assert report["resonancesCreated"] == 6


def test_import_nef_peaks_integration():
    from memops.general import Io as memopsIo
    from ccpnmr.analysis.core.PeakListImport import importTabPeaks
    
    proj_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "addpeaks")
    nef_file = os.path.join(os.path.dirname(__file__), "..", "..", "..", "sswt_assigned.nef")
    
    root = memopsIo.loadProject(proj_dir, projectName="allpeaks2")
    nmr_project = root.currentNmrProject
    spectrum = nmr_project.findFirstExperiment(name="sswt").findFirstDataSource(name="sswt-298K-hsqc-1016")
    
    _clean_project(nmr_project, spectrum)
    report = importTabPeaks(root, nef_file, spectrum)
    assert report["error"] is None
    assert report["peaksAdded"] == 116
    assert report["resonancesCreated"] == 231


def test_add_peaks_popup_select_file(monkeypatch):
    from ccpnmr.analysis.popups.AddPeaksPopup import AddPeaksPopup
    from unittest.mock import MagicMock, patch
    
    # Mock BasePopup.__init__ so it doesn't try to open a Tkinter window
    monkeypatch.setattr(AddPeaksPopup, "__init__", lambda self, *args, **kwargs: None)
    
    # Create popup instance
    popup = AddPeaksPopup(None)
    
    # Mock spectrum and dataStore with fullPath
    mock_data_store = MagicMock()
    mock_data_store.fullPath = "/path/to/some/spectrum/file.ft2"
    mock_spectrum = MagicMock()
    mock_spectrum.dataStore = mock_data_store
    
    popup.spectrum = mock_spectrum
    popup.fileEntry = MagicMock()
    popup.nameEntry = MagicMock()
    
    # Patch FileSelectPopup to verify the directory argument and simulate file selection
    dir_passed = []
    class MockFileSelectPopup:
        def __init__(self, parent, directory, file_types):
            dir_passed.append(directory)
        def getFile(self):
            return "/path/to/some/spectrum/chosen_peak_list.nef"
        def destroy(self):
            pass
            
    with patch("ccpnmr.analysis.popups.AddPeaksPopup.FileSelectPopup", MockFileSelectPopup):
        popup.selectFile()
        
    assert len(dir_passed) == 1
    assert dir_passed[0] == "/path/to/some/spectrum"
    assert popup.fileEntry.set.called
    assert popup.fileEntry.set.call_args[0][0] == "/path/to/some/spectrum/chosen_peak_list.nef"


def test_import_xeasy_peaks_integration():
    from memops.general import Io as memopsIo
    from ccpnmr.analysis.core.PeakListImport import importTabPeaks
    
    proj_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "addpeaks")
    peaks_file = os.path.join(os.path.dirname(__file__), "..", "..", "..", "sswt_assigned.peaks")
    
    root = memopsIo.loadProject(proj_dir, projectName="allpeaks2")
    nmr_project = root.currentNmrProject
    spectrum = nmr_project.findFirstExperiment(name="sswt").findFirstDataSource(name="sswt-298K-hsqc-1016")
    
    _clean_project(nmr_project, spectrum)
    report = importTabPeaks(root, peaks_file, spectrum)
    assert report["error"] is None
    assert report["peaksAdded"] == 116
    assert report["resonancesCreated"] == 222


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))

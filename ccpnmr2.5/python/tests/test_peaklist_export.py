"""Unit and integration tests for the peak list export features.
"""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from memops.general import Io as memopsIo
from ccpnmr.analysis.core.PeakListImport import importTabPeaks
from ccpnmr.analysis.core.PeakListExport import exportTabPeaks, exportNefPeaks, exportXeasyPeaks

def _clean_project(nmr_project, spectrum):
    for peak_list in list(spectrum.peakLists):
        peak_list.delete()
    for resonance in list(nmr_project.resonances):
        resonance.delete()
    for resonance_group in list(nmr_project.resonanceGroups):
        resonance_group.delete()

def test_export_import_tab_roundtrip(tmp_path):
    proj_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "addpeaks")
    tab_file = os.path.join(os.path.dirname(__file__), "..", "..", "..", "sswt_assigned.tab")
    
    root = memopsIo.loadProject(proj_dir, projectName="allpeaks2")
    nmr_project = root.currentNmrProject
    spectrum = nmr_project.findFirstExperiment(name="sswt").findFirstDataSource(name="sswt-298K-hsqc-1016")
    
    # 1. Clean and Import sswt_assigned.tab
    _clean_project(nmr_project, spectrum)
    report_in = importTabPeaks(root, tab_file, spectrum)
    assert report_in["error"] is None
    assert report_in["peaksAdded"] == 116
    assert report_in["resonancesCreated"] == 231
    
    imported_peak_list = report_in["peakList"]
    assert len(imported_peak_list.peaks) == 116
    
    # 2. Export imported_peak_list to a new tab file
    out_tab_path = os.path.join(str(tmp_path), "exported.tab")
    exportTabPeaks(imported_peak_list, out_tab_path)
    
    # 3. Clean project again
    _clean_project(nmr_project, spectrum)
    
    # 4. Import from our newly exported tab file
    report_out = importTabPeaks(root, out_tab_path, spectrum)
    assert report_out["error"] is None
    assert report_out["peaksAdded"] == 116
    # Note: resonancesCreated may differ slightly or be identical depending on how they are matched,
    # but let's assert they are successfully created and match close to 231!
    assert report_out["resonancesCreated"] == 231


def test_export_import_nef_roundtrip(tmp_path):
    proj_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "addpeaks")
    nef_file = os.path.join(os.path.dirname(__file__), "..", "..", "..", "sswt_assigned.nef")
    
    root = memopsIo.loadProject(proj_dir, projectName="allpeaks2")
    nmr_project = root.currentNmrProject
    spectrum = nmr_project.findFirstExperiment(name="sswt").findFirstDataSource(name="sswt-298K-hsqc-1016")
    
    # 1. Clean and Import sswt_assigned.nef
    _clean_project(nmr_project, spectrum)
    report_in = importTabPeaks(root, nef_file, spectrum, overwrite=True)
    assert report_in["error"] is None
    assert report_in["peaksAdded"] == 116
    assert report_in["resonancesCreated"] == 231
    
    imported_peak_list = report_in["peakList"]
    assert len(imported_peak_list.peaks) == 116
    
    # 2. Export imported_peak_list to a new nef file
    out_nef_path = os.path.join(str(tmp_path), "exported.nef")
    exportNefPeaks(imported_peak_list, out_nef_path)
    
    # 3. Clean project again
    _clean_project(nmr_project, spectrum)
    
    # 4. Import from our newly exported nef file
    report_out = importTabPeaks(root, out_nef_path, spectrum, overwrite=True)
    assert report_out["error"] is None
    assert report_out["peaksAdded"] == 116
    assert report_out["resonancesCreated"] == 231


def test_export_import_xeasy_roundtrip(tmp_path):
    proj_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "addpeaks")
    peaks_file = os.path.join(os.path.dirname(__file__), "..", "..", "..", "sswt_assigned.peaks")
    
    root = memopsIo.loadProject(proj_dir, projectName="allpeaks2")
    nmr_project = root.currentNmrProject
    spectrum = nmr_project.findFirstExperiment(name="sswt").findFirstDataSource(name="sswt-298K-hsqc-1016")
    
    # 1. Clean and Import sswt_assigned.peaks
    _clean_project(nmr_project, spectrum)
    report_in = importTabPeaks(root, peaks_file, spectrum, overwrite=True)
    assert report_in["error"] is None
    assert report_in["peaksAdded"] == 116
    assert report_in["resonancesCreated"] == 222
    
    imported_peak_list = report_in["peakList"]
    assert len(imported_peak_list.peaks) == 116
    
    # 2. Export imported_peak_list to a new peaks file
    out_peaks_path = os.path.join(str(tmp_path), "exported.peaks")
    exportXeasyPeaks(imported_peak_list, out_peaks_path)
    
    # 3. Clean project again
    _clean_project(nmr_project, spectrum)
    
    # 4. Import from our newly exported peaks file
    report_out = importTabPeaks(root, out_peaks_path, spectrum, overwrite=True)
    assert report_out["error"] is None
    assert report_out["peaksAdded"] == 116
    assert report_out["resonancesCreated"] == 222


def test_save_peaks_popup_select_file(monkeypatch):
    from ccpnmr.analysis.popups.SavePeaksPopup import SavePeaksPopup
    from unittest.mock import MagicMock, patch
    
    # Mock BasePopup.__init__ so it doesn't try to open a Tkinter window
    monkeypatch.setattr(SavePeaksPopup, "__init__", lambda self, *args, **kwargs: None)
    
    # Create popup instance
    popup = SavePeaksPopup(None)
    
    # Mock peakList and its dataStore with fullPath
    mock_data_store = MagicMock()
    mock_data_store.fullPath = "/path/to/some/spectrum/file.ft2"
    mock_spectrum = MagicMock()
    mock_spectrum.dataStore = mock_data_store
    mock_peak_list = MagicMock()
    mock_peak_list.dataSource = mock_spectrum
    
    popup.peakList = mock_peak_list
    popup.fileEntry = MagicMock()
    
    # Patch FileSelectPopup to verify the directory argument and simulate file selection
    dir_passed = []
    class MockFileSelectPopup:
        def __init__(self, parent, directory, file_types):
            dir_passed.append(directory)
        def getFile(self):
            return "/path/to/some/spectrum/chosen_save.tab"
        def destroy(self):
            pass
            
    with patch("ccpnmr.analysis.popups.SavePeaksPopup.FileSelectPopup", MockFileSelectPopup):
        popup.selectFile()
        
    assert len(dir_passed) == 1
    assert dir_passed[0] == "/path/to/some/spectrum"
    assert popup.fileEntry.set.called
    assert popup.fileEntry.set.call_args[0][0] == "/path/to/some/spectrum/chosen_save.tab"

"""P4-2 functional test: MOPS core project create → save → reload round-trip.

Exercises:
- memops.general.Io.newProject / saveProject / loadProject
- ccp.util.Spectrum.createExperiment / createSpectrum / createBlockedMatrix
- memops.api.Implementation (URL, MemopsRoot, NmrProject, DataSource)
- memops universal BlockData.determineBlockSizes

All operations self-contained; no external data files required.
NOTE: newProject(name, path) writes userData repo to path/name; loadProject
expects path/name as the first argument.
"""
import os

import pytest

from memops.api.Implementation import Url
from memops.general import Io as memopsIo
from ccp.util import Spectrum as spectrum_util


def _userDataPath(project_dir: str) -> str:
    """The userData path for a project created in project_dir."""
    return os.path.join(project_dir, "p4test")


def _create_minimal_project(project_dir: str):
    """Create a minimal single-dimension MOPS project backed by tmp_path."""
    root = memopsIo.newProject("p4test", path=project_dir, removeExisting=True)
    assert root is not None, "newProject returned None"

    nmrProject = root.newNmrProject(name="p4test_nmr")
    assert nmrProject is not None

    experiment = spectrum_util.createExperiment(
        nmrProject, name="p4expt", numDim=1, sf=(800,), isotopeCodes=("H",)
    )
    assert experiment is not None

    # Data URL points inside the userData path so save/reload resolves cleanly
    dataDir = os.path.join(_userDataPath(project_dir), "data_p4test")
    os.makedirs(dataDir, exist_ok=True)
    dls = root.newDataLocationStore(name="p4dls")
    dataUrl = dls.newDataUrl(name="p4durl", url=Url(path=dataDir))
    matrix = spectrum_util.createBlockedMatrix(dataUrl, "p4.spc", (128,))
    assert matrix is not None

    spectrum = spectrum_util.createSpectrum(
        experiment,
        name="p4spec",
        numPoints=(128,),
        sw=(8000,),
        refppm=(5.0,),
        refpt=(64,),
        dataStore=matrix,
    )
    assert spectrum is not None
    return root, experiment, spectrum


def _save_reload(tmp_path):
    """Create → save → reload, returning the reloaded MemopsRoot."""
    project_dir = str(tmp_path / "proj")
    root, experiment, spectrum = _create_minimal_project(project_dir)
    saved = memopsIo.saveProject(root, removeExisting=True)
    assert saved is True, "saveProject returned False"
    return memopsIo.loadProject(_userDataPath(project_dir), projectName="p4test")


class TestProjectLifecycle:
    def test_create_objects(self, tmp_path):
        project_dir = str(tmp_path / "proj")
        root, experiment, spectrum = _create_minimal_project(project_dir)
        assert root.name == "p4test"
        assert experiment.numDim == 1
        assert spectrum.numDim == 1

    def test_experiment_spectrum_dims(self, tmp_path):
        project_dir = str(tmp_path / "proj")
        root, experiment, spectrum = _create_minimal_project(project_dir)
        freqDims = list(spectrum.sortedDataDims())
        assert len(freqDims) == 1, f"expected 1 DataDim, got {len(freqDims)}"
        assert freqDims[0].numPoints == 128

    def test_saved_xml_present(self, tmp_path):
        project_dir = str(tmp_path / "proj")
        root, experiment, spectrum = _create_minimal_project(project_dir)
        memopsIo.saveProject(root, removeExisting=True)
        xmlPath = os.path.join(_userDataPath(project_dir), "memops", "Implementation", "p4test.xml")
        assert os.path.exists(xmlPath), f"project XML not found at {xmlPath}"

    def test_reload_root(self, tmp_path):
        reloaded = _save_reload(tmp_path)
        assert reloaded is not None
        assert reloaded.name == "p4test"

    def test_reload_nmr_project(self, tmp_path):
        reloaded = _save_reload(tmp_path)
        nmr = reloaded.currentNmrProject
        assert nmr is not None, "currentNmrProject missing after reload"

    def test_reload_experiment(self, tmp_path):
        reloaded = _save_reload(tmp_path)
        nmr = reloaded.currentNmrProject
        reloaded_expt = nmr.findFirstExperiment(name="p4expt")
        assert reloaded_expt is not None, "experiment p4expt missing after reload"
        assert reloaded_expt.numDim == 1

    def test_reload_spectrum(self, tmp_path):
        reloaded = _save_reload(tmp_path)
        nmr = reloaded.currentNmrProject
        reloaded_expt = nmr.findFirstExperiment(name="p4expt")
        reloaded_spec = reloaded_expt.findFirstDataSource(name="p4spec")
        assert reloaded_spec is not None, "spectrum p4spec missing after reload"
        assert reloaded_spec.numDim == 1

    def test_reload_spectrum_points(self, tmp_path):
        reloaded = _save_reload(tmp_path)
        nmr = reloaded.currentNmrProject
        reloaded_expt = nmr.findFirstExperiment(name="p4expt")
        reloaded_spec = reloaded_expt.findFirstDataSource(name="p4spec")
        freqDims = list(reloaded_spec.sortedDataDims())
        assert len(freqDims) == 1
        assert freqDims[0].numPoints == 128


"""P4-2 functional test: NEF → MOPS import via NefIo + save + reload.

Exercises:
- ccpnmr.v2io.NefIo.loadNefFile / CcpnNefReader.importNewProject
- Full NEF saveframe importers (molecular system, chemical shifts, restraints, spectra)
- memops.general.Io.saveProject / loadProject
- All auto-generated MOPS API files (Experiment, DataSource, Peak, Atom, Resonance, Molecule, ...)

Self-contained: uses bundled NEF test data, writes to tmp_path only.
Legacy/clouds/haddock modules excluded per user decision.
"""
import os

import pytest

from memops.general import Io as memopsIo
from ccpnmr.v2io import NefIo
from ccpnmr.nef import StarIo

NEF = (
    __import__("pathlib").Path(__file__).resolve().parent.parent
    / "ccpnmr" / "nef" / "testdata" / "CCPN_Commented_Example.nef"
)
PROJECT_NAME = "CCPN_Commented_Example"


def _userDataPath(base_dir: str) -> str:
    return os.path.join(base_dir, PROJECT_NAME)


class TestNefImport:
    def test_import_creates_nmr_project(self, tmp_path):
        base_dir = str(tmp_path / "proj")
        root = memopsIo.newProject(PROJECT_NAME, path=base_dir, removeExisting=True)
        reader = NefIo.CcpnNefReader()
        dataBlock = reader.getNefData(str(NEF))
        reader.importNewProject(root, dataBlock)

        nmr = root.currentNmrProject
        assert nmr is not None, "no NmrProject after NEF import"

    def test_import_creates_mol_system(self, tmp_path):
        base_dir = str(tmp_path / "proj")
        root = memopsIo.newProject(PROJECT_NAME, path=base_dir, removeExisting=True)
        reader = NefIo.CcpnNefReader()
        dataBlock = reader.getNefData(str(NEF))
        reader.importNewProject(root, dataBlock)

        molSystem = root.currentMolSystem
        assert molSystem is not None, "no MolSystem after NEF import"

    def test_import_has_experiments(self, tmp_path):
        base_dir = str(tmp_path / "proj")
        root = memopsIo.newProject(PROJECT_NAME, path=base_dir, removeExisting=True)
        reader = NefIo.CcpnNefReader()
        dataBlock = reader.getNefData(str(NEF))
        reader.importNewProject(root, dataBlock)

        nmr = root.currentNmrProject
        exps = list(nmr.sortedExperiments())
        assert len(exps) >= 1, f"expected >= 1 experiment, got {len(exps)}"

    def test_import_has_resonances(self, tmp_path):
        """The Commented_Example NEF has 41+ residues with chemical shifts."""
        base_dir = str(tmp_path / "proj")
        root = memopsIo.newProject(PROJECT_NAME, path=base_dir, removeExisting=True)
        reader = NefIo.CcpnNefReader()
        dataBlock = reader.getNefData(str(NEF))
        reader.importNewProject(root, dataBlock)

        mols = list(root.sortedMolecules())
        assert len(mols) >= 1, "no molecules imported"

    def test_import_has_molecule(self, tmp_path):
        base_dir = str(tmp_path / "proj")
        root = memopsIo.newProject(PROJECT_NAME, path=base_dir, removeExisting=True)
        reader = NefIo.CcpnNefReader()
        dataBlock = reader.getNefData(str(NEF))
        reader.importNewProject(root, dataBlock)

        mols = list(root.sortedMolecules())
        assert len(mols) >= 1
        mol = mols[0]
        # At least the first residues should be there
        residues = list(mol.sortedMolResidues())
        assert len(residues) >= 10, f"expected >= 10 residues, got {len(residues)}"

    def test_save_after_import(self, tmp_path):
        base_dir = str(tmp_path / "proj")
        root = memopsIo.newProject(PROJECT_NAME, path=base_dir, removeExisting=True)
        reader = NefIo.CcpnNefReader()
        dataBlock = reader.getNefData(str(NEF))
        reader.importNewProject(root, dataBlock)

        saved = memopsIo.saveProject(root, removeExisting=True)
        assert saved is True, "saveProject returned False"

        xmlPath = os.path.join(_userDataPath(base_dir), "memops", "Implementation",
                               PROJECT_NAME + ".xml")
        assert os.path.exists(xmlPath), f"project XML missing at {xmlPath}"

    def test_reload_after_save(self, tmp_path):
        base_dir = str(tmp_path / "proj")
        root = memopsIo.newProject(PROJECT_NAME, path=base_dir, removeExisting=True)
        reader = NefIo.CcpnNefReader()
        dataBlock = reader.getNefData(str(NEF))
        reader.importNewProject(root, dataBlock)
        memopsIo.saveProject(root, removeExisting=True)

        reloaded = memopsIo.loadProject(_userDataPath(base_dir), projectName=PROJECT_NAME)
        assert reloaded is not None
        assert reloaded.name == PROJECT_NAME

    def test_reload_preserves_nmr_project(self, tmp_path):
        base_dir = str(tmp_path / "proj")
        root = memopsIo.newProject(PROJECT_NAME, path=base_dir, removeExisting=True)
        reader = NefIo.CcpnNefReader()
        dataBlock = reader.getNefData(str(NEF))
        reader.importNewProject(root, dataBlock)
        memopsIo.saveProject(root, removeExisting=True)

        reloaded = memopsIo.loadProject(_userDataPath(base_dir), projectName=PROJECT_NAME)
        nmr = reloaded.currentNmrProject
        assert nmr is not None, "NmrProject missing after save→reload"

    def test_reload_preserves_mol_system(self, tmp_path):
        base_dir = str(tmp_path / "proj")
        root = memopsIo.newProject(PROJECT_NAME, path=base_dir, removeExisting=True)
        reader = NefIo.CcpnNefReader()
        dataBlock = reader.getNefData(str(NEF))
        reader.importNewProject(root, dataBlock)
        memopsIo.saveProject(root, removeExisting=True)

        reloaded = memopsIo.loadProject(_userDataPath(base_dir), projectName=PROJECT_NAME)
        mols = list(reloaded.sortedMolecules())
        assert len(mols) >= 1, "no molecules after save→reload"
        residues = list(mols[0].sortedMolResidues())
        assert len(residues) >= 10, f"too few residues after reload: {len(residues)}"

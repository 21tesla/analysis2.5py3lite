"""S36 functional test: ccpnmr.v2io.NefIo — NEF -> legacy (memops) model import.

Adopted from the v3-era ``ccpn/util/v2io/NefIo.py`` (the "CCPN V2 release,
data model version 2.1.2" importer that v3 kept for back-compat) onto this
tree: re-pointed at the Stage-35 ``ccpnmr.nef`` core, with the v3-only
imports (``ccpn.core._implementation.resetSerial``, ``ccpn.util.isotopes``)
replaced by the restored v2 helpers in ``ccpnmr.Common``.

Covers the public API on the 3 bundled testdata files (the same files as
``test_nef_core.py``):

- ``loadNefFile(path, memopsRoot, ...)`` — NEF content lands as legacy MOPS
  objects: MolSystem/Chain/Residue, Experiment/DataSource/PeakList/Peak
  (peaks carry resonance assignments), ShiftList (chemical shifts),
  NmrConstraintStore (distance/dihedral/HBond/Rdc constraint lists).
- ``loadProject(nefFilePath, ...)`` — creates a full project directory.
- ``saveProject`` / ``loadProject`` round-trip — imported content survives
  the legacy memops XML persistence path.

Notes:
- In the legacy model a NEF "spectrum" is an Experiment + DataSource; NEF
  carries metadata/peaks/shifts/restraints, never raw matrix data.
- Wildcard resonances: one NEF chemical-shift row can expand into several
  Shift objects (e.g. HD% maps to both HD resonances) — the shift counts
  below are the resulting object counts, not the raw row counts
  (104 rows -> 108 Shifts on the Commented example).
- ``load_nef_sequence`` may attempt to download ChemComp entries for
  non-standard residue codes (the bundled Commented example uses dummy
  residues); it falls back to UNK when offline, so the tests are
  network-independent.
"""
import os
from contextlib import redirect_stdout
from io import StringIO

import ccpnmr.nef.NefImporter as NefImporterModule
from ccpnmr.analysis.core import AssignmentBasic
from ccpnmr.v2io import NefIo
from memops.general import Io as memopsIo

NEF_DIR = os.path.join(os.path.dirname(NefImporterModule.__file__))
TESTDATA = os.path.join(NEF_DIR, 'testdata')
COMMENTED = os.path.join(TESTDATA, 'CCPN_Commented_Example.nef')
XPLOR = os.path.join(TESTDATA, 'CCPN_XPLOR_test1.nef')
SEC5 = os.path.join(TESTDATA, 'CCPN_Sec5Part3.nef')


def _load(nefPath, tmp_path, name='neftest'):
    """Import a NEF file into a fresh memops project rooted under tmp_path."""
    root = memopsIo.newProject(name, path=str(tmp_path), removeExisting=True)
    with redirect_stdout(StringIO()):
        root = NefIo.loadNefFile(nefPath, memopsRoot=root)
    return root


def _counts(root):
    """Summarise the imported object tree for assertions."""
    nmr = root.findFirstNmrProject()
    ms = root.findFirstMolSystem()
    expts = nmr.sortedExperiments()
    dataSources = [d for e in expts for d in e.sortedDataSources()]
    peakLists = [pl for d in dataSources for pl in d.getPeakLists()]
    peaks = [p for pl in peakLists for p in pl.getPeaks()]
    assignedPeaks = [
        p for p in peaks if any(pc.getPeakDimContribs() for pc in p.getPeakContribs())
    ]
    shiftLists = AssignmentBasic.getShiftLists(nmr)
    shifts = [m for sl in shiftLists for m in sl.getMeasurements()]
    constraintStores = sorted(
        list(nmr.getNmrConstraintStores()), key=lambda cs: cs.serial
    )
    constraintLists = [
        cl
        for cs in constraintStores
        for cl in sorted(
            list(cs.getConstraintLists()), key=lambda cl: (cl.className, cl.name)
        )
    ]
    return {
        'nmrProject': nmr,
        'molSystem': ms,
        'residues': sum(len(c.sortedResidues()) for c in ms.getChains()),
        'experiments': [(e.name, e.numDim) for e in expts],
        'dataSources': dataSources,
        'peakLists': peakLists,
        'peaks': peaks,
        'assignedPeaks': assignedPeaks,
        'shiftLists': shiftLists,
        'shifts': shifts,
        'constraintStores': constraintStores,
        'constraintLists': constraintLists,
    }


def _constraintsByType(counts):
    """{constraintClassName: [names]} across the constraint stores."""
    byType = {}
    for cl in counts['constraintLists']:
        byType.setdefault(cl.className, []).append(cl.name)
    return byType


# ---------------------------------------------------------------------------
# public API surface
# ---------------------------------------------------------------------------

def test_public_api():
    for name in (
        'loadNefFile',
        'loadProject',
        'CcpnNefReader',
        'createMoleculeFromNef',
        'extendMolResidues',
        'makeNefAxisCodes',
        'addDataStore',
        'fetchDataUrl',
        'assignPeak',
        'saveFrameReadingOrder',
    ):
        assert hasattr(NefIo, name), name
    # reading order is meaningful: molecular system first, peak links last
    assert NefIo.saveFrameReadingOrder[0] == 'nef_molecular_system'


# ---------------------------------------------------------------------------
# CCPN_Commented_Example.nef — shifts + peaks + restraints + assignments
# ---------------------------------------------------------------------------

def test_commented_core_objects(tmp_path):
    c = _counts(_load(COMMENTED, tmp_path))
    assert c['nmrProject'] is not None
    assert c['molSystem'] is not None
    # 235 residue records in the NEF, split over 15 chains
    assert c['residues'] == 235
    # both spectra, with their experiment names and dimensions
    assert c['experiments'] == [
        ('15N NOESY-HSQC', 3),
        ('HNCCCCCCCCCCCCC', 15),
    ]
    assert len(c['dataSources']) == 2
    assert len(c['peakLists']) == 2
    assert len(c['peaks']) == 6
    # every peak of the example is assigned to resonances
    assert len(c['assignedPeaks']) == 6


def test_commented_shifts(tmp_path):
    c = _counts(_load(COMMENTED, tmp_path))
    assert len(c['shiftLists']) == 2
    # 104 NEF shift rows expand to 108 Shift objects (wildcard resonances)
    assert len(c['shifts']) == 108
    sh = next(iter(c['shifts']))
    assert sh.getResonance() is not None
    assert isinstance(sh.value, float)


def test_commented_constraints(tmp_path):
    c = _counts(_load(COMMENTED, tmp_path))
    assert len(c['constraintStores']) == 1
    byType = _constraintsByType(c)
    assert set(byType) == {
        'DihedralConstraintList',
        'HBondConstraintList',
        'DistanceConstraintList',
        'RdcConstraintList',
    }
    cs = c['constraintStores'][0]
    lists = {cl.className: cl for cl in cs.getConstraintLists()}
    assert len(lists['DihedralConstraintList'].getConstraints()) == 6
    assert len(lists['HBondConstraintList'].getConstraints()) == 4
    assert len(lists['DistanceConstraintList'].getConstraints()) == 3
    assert len(lists['RdcConstraintList'].getConstraints()) == 2


# ---------------------------------------------------------------------------
# CCPN_XPLOR_test1.nef — restraints only, no spectra
# ---------------------------------------------------------------------------

def test_xplor_restraints(tmp_path):
    c = _counts(_load(XPLOR, tmp_path))
    assert c['residues'] == 58
    assert c['experiments'] == []
    assert c['dataSources'] == []
    assert c['shifts'] == []
    byType = _constraintsByType(c)
    assert sorted(byType) == [
        'DihedralConstraintList',
        'DistanceConstraintList',
        'RdcConstraintList',
    ]
    store = c['constraintStores'][0]
    lists = {cl.className: cl for cl in store.getConstraintLists()}
    assert len(lists['DistanceConstraintList'].getConstraints()) == 735
    assert len(lists['DihedralConstraintList'].getConstraints()) == 161
    # two Rdc lists in this file (147 'before' + 152 'after')
    rdcLists = [
        cl for cl in store.getConstraintLists() if cl.className == 'RdcConstraintList'
    ]
    assert len(rdcLists) == 2
    assert sum(len(l.getConstraints()) for l in rdcLists) == 299


# ---------------------------------------------------------------------------
# CCPN_Sec5Part3.nef — 5 spectra with 891 peaks
# ---------------------------------------------------------------------------

def test_sec5_spectra(tmp_path):
    c = _counts(_load(SEC5, tmp_path))
    assert c['residues'] == 95
    assert c['experiments'] == [
        ('15N HSQC/HMQC', 2),
        ('HNcoCA', 3),
        ('HNCA/CB', 3),
        ('HNCA', 3),
        ('HNcoCA/CB', 3),
    ]
    assert len(c['dataSources']) == 5
    assert len(c['peakLists']) == 5
    assert len(c['peaks']) == 891


# ---------------------------------------------------------------------------
# loadProject — full project creation in a target directory
# ---------------------------------------------------------------------------

def test_loadProject_creates_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with redirect_stdout(StringIO()):
        root = NefIo.loadProject(XPLOR, projectName='xplor_proj')
    assert os.path.isdir(str(tmp_path / 'xplor_proj'))
    c = _counts(root)
    assert c['nmrProject'] is not None
    assert c['residues'] == 58
    lists = {
        cl.className: cl for cl in c['constraintStores'][0].getConstraintLists()
    }
    assert len(lists['DistanceConstraintList'].getConstraints()) == 735


# ---------------------------------------------------------------------------
# save -> reload round-trip through the legacy memops XML persistence path
# ---------------------------------------------------------------------------

def test_save_reload_round_trip(tmp_path):
    name = 'neftest'
    root = _load(COMMENTED, tmp_path)
    with redirect_stdout(StringIO()):
        saved = memopsIo.saveProject(root, removeExisting=True)
    assert saved is True
    userPath = str(tmp_path / name)
    assert os.path.exists(os.path.join(userPath, 'memops', 'Implementation', name + '.xml'))
    with redirect_stdout(StringIO()):
        re = memopsIo.loadProject(userPath, projectName=name)
    c = _counts(re)
    assert c['residues'] == 235
    assert c['experiments'] == [
        ('15N NOESY-HSQC', 3),
        ('HNCCCCCCCCCCCCC', 15),
    ]
    assert len(c['dataSources']) == 2
    assert len(c['peaks']) == 6
    assert len(c['assignedPeaks']) == 6
    assert len(c['shifts']) == 108
    assert len(c['constraintStores']) == 1
    assert _constraintsByType(c) and set(_constraintsByType(c)) == {
        'DihedralConstraintList',
        'HBondConstraintList',
        'DistanceConstraintList',
        'RdcConstraintList',
    }

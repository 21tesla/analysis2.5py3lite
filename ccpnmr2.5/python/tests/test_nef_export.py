"""S37 functional test: ccpnmr.nefExport — legacy (memops) model -> NEF export.

Write-side companion tests for ``ccpnmr/v2io/NefIo.py`` (see
``test_v2io_nef.py`` for the read side).  Covers the public API
(``makeNefDataBlock`` / ``exportProject``) and the full
import -> export -> reimport round trip on the 3 bundled testdata files:

- ``CCPN_Commented_Example.nef`` — shifts + peaks + restraints + assignments
- ``CCPN_XPLOR_test1.nef`` — restraints only (no spectra, no shift lists)
- ``CCPN_Sec5Part3.nef`` — 5 spectra / 891 peaks + one shift list that
  contains unassigned-chain rows (``element@serial`` atom names)

The round-trip invariants (all verified against the live model before and
after the export): object counts (residues / experiments / data sources /
peaks / shifts / constraint lists), multiset equality of shift values,
peak positions, and per-constraint target/limit values.

Non-trivial export conventions the round trip pins down:

- ``residue_name`` must be the standard 3-letter code (the importer maps it
  through ``v2io.Constants.residueName2chemCompId``); the legacy title-case
  ``ccpCode`` is not a valid key.  For residues without a standard code the
  file falls back to ``UNK`` (the importer's own offline fallback).
- one NEF shift row may back several Shift objects: an ambiguous atom set
  (e.g. ``HG%``) expands to one resonance per atom set on import, so shifts
  of one row carry different name spellings of the same atom and must be
  exported as a single row in the canonical upper-case '%' form;
- dihedral constraints carry their limits on the *items* and the importer
  creates exactly one item per row, so every item gets its own row;
- shift rows for resonances the importer deliberately names ``None``
  (the reserved ``element@serial`` atom form) must recreate that form.

Notes:
- NEF carries metadata / peaks / shifts / restraints, never raw matrices,
  so DataStore content is out of scope here.
- The 4 ``dummy``-linked linker residues of the Commented example are
  imported as ``UNK`` (offline fallback) and therefore re-export as ``UNK``
  as well: that is a property of the import, not of the export.
- The tests are network-independent in the same sense as ``test_v2io_nef``
  (unknown codes fall back to UNK).
"""
import os
from collections import Counter
from contextlib import redirect_stdout
from io import StringIO

import ccpnmr.nef.NefImporter as NefImporterModule
from ccpnmr import nefExport
from ccpnmr.analysis.core import AssignmentBasic
from ccpnmr.nef import StarIo
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


def _export(root, tmp_path, name='out'):
    """Export the project to a NEF file under tmp_path; return the path."""
    path = os.path.join(str(tmp_path), f'{name}.nef')
    with redirect_stdout(StringIO()):
        nefExport.exportProject(root, path)
    return path


def _counts(root):
    """Summarise the object tree, mirroring test_v2io_nef._counts."""
    nmr = root.findFirstNmrProject()
    ms = root.findFirstMolSystem()
    expts = nmr.sortedExperiments()
    dataSources = [d for e in expts for d in e.sortedDataSources()]
    peaks = [p for d in dataSources for pl in d.getPeakLists() for p in pl.getPeaks()]
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
        'residues': sum(len(c.sortedResidues()) for c in ms.getChains()) if ms else 0,
        'experiments': [(e.name, e.numDim) for e in expts],
        'dataSources': dataSources,
        'peaks': peaks,
        'assignedPeaks': assignedPeaks,
        'shiftLists': shiftLists,
        'shifts': shifts,
        'constraintStores': constraintStores,
        'constraintLists': constraintLists,
    }


def _shiftValues(c):
    return sorted(round(s.value, 6) for s in c['shifts'])


def _peakPositions(c):
    out = []
    for p in c['peaks']:
        out.append(tuple(sorted(round(pd.value, 6) for pd in p.sortedPeakDims())))
    return sorted(out)


def _constraintValues(c):
    """Multiset of (list class, list name, target, lower, upper) per constraint."""
    out = []
    for cl in c['constraintLists']:
        for con in cl.getConstraints():
            if cl.className == 'DihedralConstraintList':
                items = sorted(
                    list(con.items), key=lambda x: getattr(x, 'serial', 0)
                )
                src = items[0] if items else None
            else:
                src = con
            if src is None:
                continue
            out.append((
                cl.className, cl.name, src.targetValue,
                src.lowerLimit, src.upperLimit,
            ))
    return Counter(out)


def _saveframeRows(path):
    """Parse an exported NEF file; {framecode: {loopName: numRows}}."""
    imp = NefImporterModule.NefImporter()
    imp.loadFile(path)
    out = {}
    for name, sf in imp.data.items():
        if not isinstance(sf, StarIo.NmrSaveFrame):
            continue
        out[name] = {
            key: len(loop.data)
            for key, loop in sf.items() if isinstance(loop, StarIo.NmrLoop)
        }
    return out


# ---------------------------------------------------------------------------
# public API surface
# ---------------------------------------------------------------------------

def test_public_api():
    for name in ('makeNefDataBlock', 'exportProject'):
        assert hasattr(nefExport, name), name


# ---------------------------------------------------------------------------
# exported file structure (Commented example)
# ---------------------------------------------------------------------------

def test_commented_export_file_structure(tmp_path):
    root = _load(COMMENTED, tmp_path)
    path = _export(root, tmp_path)
    assert os.path.exists(path)
    rows = _saveframeRows(path)

    # meta + molecular system (235 residues in 15 chains)
    assert rows['nef_nmr_meta_data']['nef_program_script'] == 1
    assert rows['nef_molecular_system']['nef_sequence'] == 235

    # 2 shift lists: the 4 ambiguous '%' atom-set rows are collapsed to one
    # row each (108 live shifts -> 104 rows, the original file's row count)
    assert rows['nef_chemical_shift_list_1']['nef_chemical_shift'] == 93
    assert rows['nef_chemical_shift_list_2']['nef_chemical_shift'] == 11

    # restraints: one row per constraint, EXCEPT dihedrals, which carry their
    # limits on the items (1 row per item: 1+4+2+2+1+1 = 11 rows for 6 restraints)
    assert rows['nef_distance_restraint_list_L1']['nef_distance_restraint'] == 3
    assert rows['nef_distance_restraint_list_hbond1']['nef_distance_restraint'] == 4
    assert rows['nef_dihedral_restraint_list_L2']['nef_dihedral_restraint'] == 11
    assert rows['nef_rdc_restraint_list_3']['nef_rdc_restraint'] == 2

    # 2 spectra: 3-dim (6 peaks, one with 3 alternative assignments = 9 rows)
    # and the 15-dim dummy (1 unassigned peak)
    cnoesy = rows['nef_nmr_spectrum_cnoesy1']
    assert cnoesy['nef_spectrum_dimension'] == 3
    assert cnoesy['nef_spectrum_dimension_transfer'] == 2
    assert cnoesy['nef_peak'] == 9
    dummy = rows['nef_nmr_spectrum_dummy15d']
    assert dummy['nef_spectrum_dimension'] == 15
    assert dummy['nef_peak'] == 1

    assert rows['nef_peak_restraint_links']['nef_peak_restraint_link'] == 4


# ---------------------------------------------------------------------------
# full round trips
# ---------------------------------------------------------------------------

def test_commented_round_trip(tmp_path):
    c1 = _counts(_load(COMMENTED, tmp_path, 'rtA'))
    path = _export(_load(COMMENTED, tmp_path, 'rtB'), tmp_path)
    c2 = _counts(_load(path, tmp_path, 'rtC'))

    assert c2['residues'] == 235
    assert c2['experiments'] == [
        ('15N NOESY-HSQC', 3),
        ('HNCCCCCCCCCCCCC', 15),
    ]
    assert len(c2['dataSources']) == 2
    assert len(c2['peaks']) == 6
    assert len(c2['assignedPeaks']) == 6
    assert len(c2['shiftLists']) == 2
    # 104 NEF shift rows expand to 108 Shift objects (wildcard resonances)
    assert len(c2['shifts']) == len(c1['shifts']) == 108
    assert _shiftValues(c1) == _shiftValues(c2)
    assert _peakPositions(c1) == _peakPositions(c2)

    lists = {cl.className: len(cl.getConstraints()) for cl in c2['constraintLists']}
    assert lists == {
        'DihedralConstraintList': 6,
        'HBondConstraintList': 4,
        'DistanceConstraintList': 3,
        'RdcConstraintList': 2,
    }
    assert _constraintValues(c1) == _constraintValues(c2)


def test_xplor_round_trip(tmp_path):
    c1 = _counts(_load(XPLOR, tmp_path, 'rtA'))
    path = _export(_load(XPLOR, tmp_path, 'rtB'), tmp_path)
    c2 = _counts(_load(path, tmp_path, 'rtC'))

    assert c2['residues'] == 58
    assert c2['experiments'] == []
    assert c2['dataSources'] == []
    assert c2['shifts'] == []

    store = c2['constraintStores'][0]
    lists = {cl.className: cl for cl in store.getConstraintLists()}
    assert len(lists['DistanceConstraintList'].getConstraints()) == 735
    assert len(lists['DihedralConstraintList'].getConstraints()) == 161
    rdcLists = [cl for cl in store.getConstraintLists() if cl.className == 'RdcConstraintList']
    assert len(rdcLists) == 2
    assert sum(len(l.getConstraints()) for l in rdcLists) == 299
    assert _constraintValues(c1) == _constraintValues(c2)


def test_sec5_round_trip(tmp_path):
    c1 = _counts(_load(SEC5, tmp_path, 'rtA'))
    path = _export(_load(SEC5, tmp_path, 'rtB'), tmp_path)
    c2 = _counts(_load(path, tmp_path, 'rtC'))

    assert c2['residues'] == 95
    assert c2['experiments'] == [
        ('15N HSQC/HMQC', 2),
        ('HNcoCA', 3),
        ('HNCA/CB', 3),
        ('HNCA', 3),
        ('HNcoCA/CB', 3),
    ]
    assert len(c2['peaks']) == 891
    assert len(c2['shiftLists']) == 1
    # the shift list carries 'element@serial' atom names for unassigned chains
    assert len(c2['shifts']) == len(c1['shifts']) == 542
    assert _shiftValues(c1) == _shiftValues(c2)
    assert _peakPositions(c1) == _peakPositions(c2)

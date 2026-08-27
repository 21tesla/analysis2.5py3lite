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

Stage 39c adds the one exception: a DataSource LINKED to a data file
exports the reference (ccpn_spectrum_file_path + ccpn_file_* + the
data-dimension point counts) so a plain same-machine reimport auto-links
it (test_linked_datasource_export_reimport); unlinked DataSources export
exactly as before (test_unlinked_export_unchanged).

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

import pytest

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


def test_duplicate_datasource_names_export(tmp_path):
    """Two DataSources sharing a name (typical of NMRpipe imports, e.g.
    two datasets both named 'ftt') used to crash the export with
    'duplicate key name nef_nmr_spectrum_ftt' - framecodes must be
    unique per file, so the second one gets a '_2' suffix."""
    root = _load(COMMENTED, tmp_path)
    c1 = _counts(root)
    expts = root.findFirstNmrProject().sortedExperiments()
    ds_a = expts[0].findFirstDataSource()
    expts[1].findFirstDataSource().name = ds_a.name
    path = _export(root, tmp_path)
    rows = _saveframeRows(path)
    base = f'nef_nmr_spectrum_{ds_a.name}'
    spectra = sorted(k for k in rows if k.startswith('nef_nmr_spectrum_'))
    assert spectra == [base, f'{base}_2']

    # the exported file re-imports with both data sources intact
    c2 = _counts(_load(path, tmp_path, 'dup'))
    assert len(c2['dataSources']) == len(c1['dataSources']) == 2
    assert _peakPositions(c1) == _peakPositions(c2)
    assert _shiftValues(c1) == _shiftValues(c2)


def test_export_no_shiftlists_shared_placeholder(tmp_path, monkeypatch):
    """A project without any shift lists but with several spectra used to
    synthesize one empty placeholder shift-list frame per spectrum,
    crashing on the second identical framecode - it must be created once
    and shared by every spectrum saveframe."""
    root = _load(COMMENTED, tmp_path)
    monkeypatch.setattr(nefExport, '_shiftLists', lambda nmrProject: [])
    path = _export(root, tmp_path)

    imp = NefImporterModule.NefImporter()
    imp.loadFile(path)
    frames = {
        sf['sf_framecode']: sf
        for sf in imp.data.values()
        if isinstance(sf, StarIo.NmrSaveFrame)
    }
    placeholders = [k for k in frames if k.startswith('nef_chemical_shift_list_')]
    assert placeholders == ['nef_chemical_shift_list_1']
    assert len(frames[placeholders[0]]['nef_chemical_shift'].data) == 0
    spectra = [sf for sf in frames.values() if sf['sf_category'] == 'nef_nmr_spectrum']
    assert len(spectra) == 2
    for sf in spectra:
        assert sf['chemical_shift_list'] == placeholders[0]


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


def test_native_legacy_project_round_trip(tmp_path):
    """Regression: in a NATIVE legacy project the ResonanceGroups carry
    no name - their identity is the linked MolSystem residue - so the
    export must derive chain/sequence from it (and fall back to the
    reader's '@' unassigned-chain form for groups it cannot place),
    otherwise the reader drops every peak assignment and collapses all
    shifts into the default-chain group on reimport."""

    class _Stub:
        """Minimal attribute bag for _resonanceIdentity branch tests."""
        def __init__(self, **kw):
            self.__dict__.update(kw)

    # branch coverage for _resonanceIdentity
    namedGroup = _Stub(name='A.63', ccpCode='Gly', residue=None, serial=5)
    reso = _Stub(name='H', isotopeCode='1H', serial=1, resonanceGroup=namedGroup)
    assert nefExport._resonanceIdentity(reso) == ('A', '63', 'Gly', 'H')

    chainlessGroup = _Stub(name='63', ccpCode='Gly', residue=None, serial=5)
    reso = _Stub(name='H', isotopeCode='1H', serial=1, resonanceGroup=chainlessGroup)
    assert nefExport._resonanceIdentity(reso) == ('@', '63', 'Gly', 'H')

    unnamedGroup = _Stub(name=None, ccpCode='Val', serial=5, residue=_Stub(
        chain=_Stub(code='A'), seqCode=100, seqInsertCode=' ', chemCompVar=None))
    reso = _Stub(name='H', isotopeCode='1H', serial=1, resonanceGroup=unnamedGroup)
    assert nefExport._resonanceIdentity(reso) == ('A', '100', 'Val', 'H')

    unnamedGroup = _Stub(name=None, ccpCode='Ala', serial=5, residue=_Stub(
        chain=_Stub(code='B'), seqCode=7, seqInsertCode='A',
        chemCompVar=_Stub(chemComp=_Stub(code3Letter='ALA'))))
    reso = _Stub(name='H', isotopeCode='1H', serial=1, resonanceGroup=unnamedGroup)
    assert nefExport._resonanceIdentity(reso) == ('B', '7A', 'ALA', 'H')

    # group the reader cannot place: pinned by the group's own serial
    placelessGroup = _Stub(name=None, ccpCode='Gly', residue=None, serial=7)
    reso = _Stub(name='H', isotopeCode='1H', serial=1, resonanceGroup=placelessGroup)
    assert nefExport._resonanceIdentity(reso) == ('@', '@7', 'Gly', 'H')

    # resonance without any group: pinned by the resonance serial (the atom
    # form carries the element symbol, not the isotope code: 15N -> N@278)
    reso = _Stub(name=None, isotopeCode='15N', serial=278, resonanceGroup=None)
    assert nefExport._resonanceIdentity(reso) == ('@', '@278', None, 'N@278')

    # no group, no serial: nothing to derive
    reso = _Stub(name='H', isotopeCode='1H', serial=None, resonanceGroup=None)
    assert nefExport._resonanceIdentity(reso) == (None, None, None, 'H')

    # integration: strip every group name (the native-project state) and
    # verify a lossless round trip
    def stripRgNames(rt):
        n = 0
        for rg in rt.findFirstNmrProject().resonanceGroups:
            rg.name = None
            n += 1
        return n

    def shiftIdentityMultiset(rt):
        """(chain, seq, residue, atom, value) -> count over the shift list.

        The identity the export derives (the fix under test) is checked
        value-for-value: the writer's 10-significant-figure float format
        cannot perturb the 4th decimal, so quantising to 1e-4 is exact.
        (Peak-dimension identities are NOT compared: the importer long
        pre-dates this fix and deliberately collapses alternatives that
        differ in a single dimension, which the round trip does not
        preserve - a property of the reader, not of the export.)
        """
        from collections import Counter

        def rid(reso):
            rg = reso.resonanceGroup
            res = rg.residue if rg is not None else None
            chain = res.chain.code if res is not None and res.chain is not None else '?'
            seq = str(res.seqCode) if res is not None else '?'
            resname = (res.ccpCode if res is not None
                       else (rg.ccpCode if rg is not None else '?')) or '?'
            atom = reso.name if reso.name is not None else f'{reso.isotopeCode}@{reso.serial}'
            return (chain, seq, resname, atom)

        out = Counter()
        for sl in AssignmentBasic.getShiftLists(rt.findFirstNmrProject()):
            for s in sl.getMeasurements():
                out[(rid(s.resonance), round(s.value, 4))] += 1
        return out

    root = _load(COMMENTED, tmp_path)
    c1 = _counts(root)
    assert stripRgNames(root) > 0

    path = _export(root, tmp_path)

    # every exported assignment row carries an interpretable identity
    imp = NefImporterModule.NefImporter()
    imp.loadFile(path)
    for sf in imp.data.values():
        if not isinstance(sf, StarIo.NmrSaveFrame):
            continue
        if sf['sf_category'] == 'nef_chemical_shift_list':
            for row in sf['nef_chemical_shift'].data:
                assert row['chain_code'] is not None, row
                assert row['sequence_code'] is not None, row
        elif sf['sf_category'] == 'nef_nmr_spectrum':
            loop = sf.get('nef_peak')
            if loop is None:
                continue
            for row in loop.data:
                for col in loop.columns:
                    if col.startswith('atom_name_') and row[col] is not None:
                        dim = col.rsplit('_', 1)[-1]
                        assert row[f'chain_code_{dim}'] is not None, row
                        assert row[f'sequence_code_{dim}'] is not None, row
                        break

    # reimport: no dropped assignments, no warnings, same identities
    reRoot = memopsIo.newProject(
        'nativeRe', path=os.path.join(str(tmp_path), 'reimport'), removeExisting=True)
    log = StringIO()
    with redirect_stdout(log):
        NefIo.loadNefFile(path, memopsRoot=reRoot)
    assert 'Uninterpretable Peak assignment' not in log.getvalue()

    c2 = _counts(reRoot)
    assert len(c2['peaks']) == len(c1['peaks'])
    assert len(c2['assignedPeaks']) == len(c1['assignedPeaks'])
    assert len(c2['shifts']) == len(c1['shifts'])
    assert shiftIdentityMultiset(root) == shiftIdentityMultiset(reRoot)
    assert _peakPositions(c1) == _peakPositions(c2)


# ---------------------------------------------------------------------------
# stage 39c: file reference in the NEF -> same-machine auto-relink
# ---------------------------------------------------------------------------

def _spectrumSaveframes(path):
    """{framecode: NmrSaveFrame} of the exported nef_nmr_spectrum frames."""
    imp = NefImporterModule.NefImporter()
    imp.loadFile(path)
    return {
        name: sf for name, sf in imp.data.items()
        if isinstance(sf, StarIo.NmrSaveFrame)
        and sf['sf_category'] == 'nef_nmr_spectrum'
    }


def test_linked_datasource_export_reimport(tmp_path):
    """39c round trip: relink -> export -> PLAIN reimport (no relink step)
    restores the data file link, the restored data-dim geometry (via
    point_count) and the peak positions on that grid."""
    from ccpnmr import nefRelink
    from tests.test_nef_relink import makePipeFile

    root = _load(COMMENTED, tmp_path)
    ds3 = None
    for expt in root.findFirstNmrProject().sortedExperiments():
        if expt.numDim == 3:
            ds3 = expt.findFirstDataSource()
    assert ds3 is not None and ds3.dataStore is None
    # the ppm values carried by the NEF - they must survive the whole
    # relink -> export -> import cycle
    valuesBefore = {
        (peak.serial, peakDim.dim): peakDim.value
        for peakList in ds3.getPeakLists()
        for peak in peakList.getPeaks()
        for peakDim in peak.sortedPeakDims()
    }

    base = tmp_path / 'yb-demo' / 'noesyn'
    base.mkdir(parents=True)
    dataFile = str(base / 'cnoesy1.ft3')
    makePipeFile(dataFile, (64, 32, 16),
                 sf=(600.0, 150.0, 100.0), sw=(8000.0, 2000.0, 1200.0),
                 nuc=('1H', '15N', '13C'))
    report = nefRelink.relinkSpectra(root, str(tmp_path / 'yb-demo'))
    assert [e['name'] for e in report['linked']] == ['cnoesy1']

    path = _export(root, tmp_path)
    frames = _spectrumSaveframes(path)
    cnoesy = frames['nef_nmr_spectrum_cnoesy1']

    # the linked spectrum carries the file reference + format items and
    # the point counts on the dimension rows
    assert cnoesy.get('ccpn_spectrum_file_path') == dataFile
    assert cnoesy.get('ccpn_file_type') == 'NmrPipe'
    assert cnoesy.get('ccpn_file_header_size') == 2048
    assert cnoesy.get('ccpn_file_byte_number') == 4
    assert cnoesy.get('ccpn_file_number_type') == 'float'
    assert cnoesy.get('ccpn_file_is_big_endian') is False
    # the point counts ride on the ccpn extension loop (the importer
    # reads them from there, not from nef_spectrum_dimension)
    pts = {row['dimension_id']: (row['point_count'], row['total_point_count'])
           for row in cnoesy['ccpn_spectrum_dimension'].data}
    assert pts == {1: (64, 64), 2: (32, 32), 3: (16, 16)}

    # the unlinked 15-dim dummy exports exactly as before 39c
    dummy = frames['nef_nmr_spectrum_dummy15d']
    assert dummy.get('ccpn_spectrum_file_path') is None
    assert dummy.get('ccpn_file_type') is None
    assert dummy.get('ccpn_spectrum_dimension') is None
    for row in dummy['nef_spectrum_dimension'].data:
        assert row.get('point_count') is None

    # plain reimport: the link comes WITH the file, 0 import warnings
    reRoot = memopsIo.newProject(
        'linkedRe', path=os.path.join(str(tmp_path), 'reimport'),
        removeExisting=True)
    log = StringIO()
    with redirect_stdout(log):
        NefIo.loadNefFile(path, memopsRoot=reRoot)
    assert '====>' not in log.getvalue()

    reDs = None
    for expt in reRoot.findFirstNmrProject().sortedExperiments():
        if expt.numDim == 3:
            reDs = expt.findFirstDataSource()
    assert reDs.dataStore is not None
    assert reDs.dataStore.fileType == 'NmrPipe'
    assert reDs.dataStore.fullPath == dataFile
    assert os.path.exists(reDs.dataStore.fullPath)
    dims = reDs.sortedDataDims()
    assert [d.numPoints for d in dims] == [64, 32, 16]
    assert [d.valuePerPoint for d in dims] == pytest.approx([125.0, 62.5, 75.0])
    # restored grid: model constraint satisfied, ppm values preserved
    for peakList in reDs.getPeakLists():
        for peak in peakList.getPeaks():
            for peakDim in peak.sortedPeakDims():
                assert 1.0 <= peakDim.position < \
                    dims[peakDim.dim - 1].numPoints + 1.0
                assert peakDim.value == pytest.approx(
                    valuesBefore[(peak.serial, peakDim.dim)], abs=1e-9)


def test_unlinked_export_unchanged(tmp_path):
    """The 39c guard: with no linked DataSources the spectrum frames
    carry no file items and no point_count columns (pre-39c output)."""
    root = _load(COMMENTED, tmp_path)
    path = _export(root, tmp_path)
    frames = _spectrumSaveframes(path)
    assert set(frames) == {'nef_nmr_spectrum_cnoesy1', 'nef_nmr_spectrum_dummy15d'}
    for sf in frames.values():
        assert sf.get('ccpn_spectrum_file_path') is None
        assert sf.get('ccpn_file_type') is None
        assert sf.get('ccpn_spectrum_dimension') is None
        for row in sf['nef_spectrum_dimension'].data:
            assert row.get('point_count') is None
            assert row.get('total_point_count') is None

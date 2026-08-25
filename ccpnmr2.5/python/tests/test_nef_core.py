"""S35 functional test: ccpnmr.nef — model-free BMRB NEF (v1.1) format core.

Adopted from CCPNMR v3 (ccpn/util/nef).  Covers the model-free stack:
StarTokeniser -> GenericStarParser -> StarIo -> NefImporter, plus the
dictionary-based Validator, using the 3 bundled testdata files.
"""
import os

import pytest

from ccpnmr.nef import StarIo
from ccpnmr.nef import NefImporter as NefImporterModule
from ccpnmr.nef.NefImporter import NefImporter

NEF_DIR = os.path.join(os.path.dirname(NefImporterModule.__file__))
TESTDATA = os.path.join(NEF_DIR, 'testdata')
COMMENTED = os.path.join(TESTDATA, 'CCPN_Commented_Example.nef')
XPLOR = os.path.join(TESTDATA, 'CCPN_XPLOR_test1.nef')
SEC5 = os.path.join(TESTDATA, 'CCPN_Sec5Part3.nef')


def _rows(loop):
    """Row count of an NmrLoop (its .data list)."""
    if loop is None:
        return 0
    if isinstance(loop, StarIo.NmrLoop):
        return len(loop.data)
    if isinstance(loop, (list, tuple)):
        return len(loop)
    return 0


@pytest.fixture(scope='module')
def commented():
    imp = NefImporter()
    imp.loadFile(COMMENTED)
    return imp


@pytest.fixture(scope='module')
def xplor():
    imp = NefImporter()
    imp.loadFile(XPLOR)
    return imp


@pytest.fixture(scope='module')
def sec5():
    imp = NefImporter()
    imp.loadFile(SEC5)
    return imp


# ---------------------------------------------------------------------------
# parsing the bundled files
# ---------------------------------------------------------------------------

def test_all_bundled_files_parse():
    for path in (COMMENTED, XPLOR, SEC5):
        imp = NefImporter()
        imp.loadFile(path)
        names = imp.getSaveFrameNames()
        assert names, f'no saveFrames parsed from {path}'
        assert 'nef_nmr_meta_data' in imp.data
        assert 'nef_molecular_system' in imp.data
        assert os.path.samefile(imp.path, path)


def test_meta_and_molecular_system(commented):
    meta = commented.data['nef_nmr_meta_data']
    assert meta.get('format_name') == 'nmr_exchange_format'
    # StarIo coerces '1.1' to float 1.1
    assert str(meta.get('format_version')) in ('1.1', '1.2')
    seq = commented.data['nef_molecular_system'].get('nef_sequence')
    assert _rows(seq) == 235  # the commented example: 235 residue records


def test_chemical_shift_accessors(commented):
    csls = commented.getChemicalShiftLists()
    assert isinstance(csls, list) and len(csls) == 2
    shifts = sum(_rows(c.get('nef_chemical_shift')) for c in csls)
    assert shifts == 104


def test_nmr_spectrum_accessor(commented):
    ns = commented.getNmrSpectra()
    assert isinstance(ns, list) and len(ns) == 2
    names = {getattr(fr, 'name', '') for fr in ns}
    assert any('cnoesy' in n.lower() for n in names)


def test_xplor_has_no_spectra_but_restraints(xplor):
    assert xplor.getNmrSpectra() is None
    assert xplor.getChemicalShiftLists() is None
    drls = xplor.getDistanceRestraintLists()
    if not isinstance(drls, list):
        drls = [drls]
    assert sum(_rows(d.get('nef_distance_restraint')) for d in drls) == 735


def test_sec5_spectra(sec5):
    ns = sec5.getNmrSpectra()
    assert isinstance(ns, list) and len(ns) == 5
    assert sec5.isValid is True  # validator accepts it (5 informational notes)
    assert len(sec5.validErrorLog) == 5


def test_categories_listed(commented):
    assert commented.getCategories() == (
        'nmr_meta_data', 'molecular_system', 'chemical_shift_list',
        'distance_restraint_list', 'dihedral_restraint_list',
        'rdc_restraint_list', 'nmr_spectrum', 'peak_restraint_links')


# ---------------------------------------------------------------------------
# saveframe management & round-trips
# ---------------------------------------------------------------------------

def test_saveframe_management(commented):
    assert commented.hasSaveFrame('nmr_meta_data')
    assert commented.hasSaveFrame('nef_nmr_meta_data')  # prefixed form accepted
    fr = commented.getSaveFrame('chemical_shift_list_1')
    assert fr is not None
    assert 'chemical_shift' in fr.getTableNames()  # 'nef_' prefix hidden
    tbl = fr.getTable('chemical_shift')
    assert tbl, 'chemical_shift loop should have rows'
    assert 'chain_code' in tbl[0]


def test_tostring_fromstring_roundtrip(commented, tmp_path):
    text = commented.toString()
    reparsed = NefImporter()
    reparsed.fromString(text)
    assert set(reparsed.getSaveFrameNames()) == set(commented.getSaveFrameNames())
    csls = reparsed.getChemicalShiftLists()
    if not isinstance(csls, list):
        csls = [csls] if csls else []
    assert sum(_rows(c.get('nef_chemical_shift')) for c in csls) == 104


def test_savefile_loadfile_roundtrip(commented, tmp_path):
    out = str(tmp_path / 'copied.nef')
    assert commented.saveFile(out)
    reparsed = NefImporter()
    reparsed.loadFile(out)
    assert set(reparsed.getSaveFrameNames()) == set(commented.getSaveFrameNames())


def test_stariot_parse_direct():
    extent = StarIo.parseNefFile(fileName=SEC5)
    dbs = list(extent.values())
    assert len(dbs) == 1
    db = dbs[0]
    assert 'nef_nmr_spectrum_hsqc' in db
    assert 'nef_peak' in db['nef_nmr_spectrum_hsqc']  # NEF v1.1 peak loop


def test_add_and_delete_saveframe():
    imp = NefImporter()  # fresh, empty
    sf = imp.addSaveFrame('my_test_frame', 'nef_nmr_meta_data')
    assert imp.hasSaveFrame('my_test_frame')
    assert imp.deleteSaveFrame('my_test_frame') is True
    assert not imp.hasSaveFrame('my_test_frame')

"""S39a functional test: ccpnmr.nefRelink - automatic relinking of spectrum
data files for NEF-imported projects.

Covers the public API (``scanSpectrumFiles`` / ``matchSpectra`` /
``relinkSpectra``) with synthetic NMRpipe files: ``_makePipeHeader``
crafts the minimal 2048-byte header layout that
``NmrPipeParams.getHeader`` reads (magic + fixed float-word indices), so
no real spectrum data is needed.  The end-to-end test imports the
bundled Commented example, relinks its 3-dim data source against a
synthetic dataset directory, and checks the dataStore link, the restored
data-dimension geometry, and the >4-dim skip / second-pass idempotency.
"""
import os
import struct
from contextlib import redirect_stdout
from io import StringIO

import ccpnmr.nef.NefImporter as NefImporterModule
from ccp.format.spectra.params import NmrPipeParams
from ccpnmr import nefRelink
from ccpnmr.v2io import NefIo
from memops.general import Io as memopsIo

NEF_DIR = os.path.dirname(NefImporterModule.__file__)
TESTDATA = os.path.join(NEF_DIR, 'testdata')
COMMENTED = os.path.join(TESTDATA, 'CCPN_Commented_Example.nef')


def _load(nefPath, tmp_path, name='neftest'):
    """Import a NEF file into a fresh memops project rooted under tmp_path."""
    root = memopsIo.newProject(name, path=str(tmp_path), removeExisting=True)
    with redirect_stdout(StringIO()):
        root = NefIo.loadNefFile(nefPath, memopsRoot=root)
    return root


# ---------------------------------------------------------------------------
# synthetic NMRpipe files
# ---------------------------------------------------------------------------

def _makePipeHeader(npts, sw=None, sf=None, nuc=None):
    """A minimal 2048-byte (512 x float32) little-endian NMRpipe header."""
    P = NmrPipeParams
    ndim = len(npts)
    words = [0.0] * (P.head // 4)
    words[P.ndim_index] = ndim
    for i in range(ndim):
        words[P.order_index[i]] = i + 1  # identity dimension order
        words[P.complex_index[i]] = 1.0  # real data (must be non-zero)
        words[P.npts_index[i]] = npts[i]
        words[P.sw_index[i]] = (sw or [10000.0] * ndim)[i]
        words[P.sf_index[i]] = (sf or [600.0, 150.0, 100.0][:ndim])[i]
        words[P.origin_index[i]] = 0.0
    s = struct.pack(f'<{P.head // 4}f', *words)
    # the byte-order magic occupies raw bytes 8-11; for a little-endian
    # file it is the reversed big-endian marker
    s = s[:8] + bytes([0x7B, 0x14, 0x16, 0x40]) + s[12:]
    # nucleus strings sit at raw byte offsets 4*nuc_index[i] (ASCII,
    # independent of the float endianness)
    for i, n in enumerate((nuc or ['1H'] * ndim)[:ndim]):
        blob = n.encode('ascii').ljust(4)
        off = 4 * P.nuc_index[i]
        s = s[:off] + blob + s[off + 4:]
    return s


def makePipeFile(path, npts, sw=None, sf=None, nuc=None):
    """Write a synthetic NMRpipe data file (header + stub data region)."""
    header = _makePipeHeader(npts, sw=sw, sf=sf, nuc=nuc)
    size = 4
    for n in npts:
        size *= n
    with open(path, 'wb') as f:
        f.write(header)
        f.write(b'\x00' * size)  # stub matrix: only the header is read


# ---------------------------------------------------------------------------
# public API surface + header fixture sanity
# ---------------------------------------------------------------------------

def test_public_api():
    for name in ('scanSpectrumFiles', 'matchSpectra', 'relinkSpectra'):
        assert hasattr(nefRelink, name), name


def test_makePipeFile_readable(tmp_path):
    path = str(tmp_path / 'h.ft3')
    makePipeFile(path, (16, 8, 4),
                 sf=(600.0, 150.0, 100.0), nuc=('1H', '15N', '13C'))
    p = NmrPipeParams.NmrPipeParams(path)
    assert p.ndim == 3
    assert p.npts[:3] == [16, 8, 4]
    assert p.block[:3] == [16, 1, 1]
    assert p.sf[:3] == [600.0, 150.0, 100.0]
    assert p.nuc[:3] == ['1H', '15N', '13C']
    assert p.big_endian is False
    assert p.head == 2048


# ---------------------------------------------------------------------------
# scanSpectrumFiles
# ---------------------------------------------------------------------------

def test_scanSpectrumFiles(tmp_path):
    d = tmp_path / 'yb' / 'sswt-hsqc'
    d.mkdir(parents=True)
    makePipeFile(str(d / 'sswt-298K-hsqc-1016.ft2'), (427, 160))
    (tmp_path / 'yb' / 'sswt-hsqc' / 'notes.txt').write_text('not a spectrum')
    # deterministic non-NMRpipe bytes (first float word != 0)
    (tmp_path / 'yb' / 'sswt-hsqc' / 'garbage.bin').write_bytes(
        b'\x00' + b'\xff' * 2999)

    cands = nefRelink.scanSpectrumFiles(str(tmp_path))
    assert len(cands) == 1
    c = cands[0]
    assert c['stem'] == 'sswt-298K-hsqc-1016'
    assert c['relDir'] == os.path.join('yb', 'sswt-hsqc')
    assert c['ndim'] == 2
    assert c['npts'] == [427, 160]
    assert c['block'] == [427, 1]


# ---------------------------------------------------------------------------
# matchSpectra
# ---------------------------------------------------------------------------

def _ds3(root):
    """The 3-dim data source of the Commented example ('cnoesy1')."""
    nmr = root.findFirstNmrProject()
    for expt in nmr.sortedExperiments():
        if expt.numDim == 3:
            return expt.findFirstDataSource()
    raise AssertionError('no 3-dim experiment found')


def test_matchSpectra_exactStemAndDimGate(tmp_path):
    root = _load(COMMENTED, tmp_path)
    ds = _ds3(root)
    assert ds.name == 'cnoesy1' and ds.numDim == 3

    near = str(tmp_path / 'cnoesy1.ft3')
    makePipeFile(near, (64, 32, 16))
    far = str(tmp_path / 'xb' / 'sswt' / 'cnoesy1.ft3')
    os.makedirs(os.path.dirname(far))
    makePipeFile(far, (64, 32, 16))
    # same stem, wrong dimension count: the ndim gate must reject it
    wrongDim = str(tmp_path / 'cnoesy1.ft2')
    makePipeFile(wrongDim, (64, 32))

    cands = nefRelink.scanSpectrumFiles(str(tmp_path))
    assert len(cands) == 3
    # exact stem beats the path-order tie: both 3-dim files score 100,
    # the shorter path comes first in the path-sorted candidate list
    assert nefRelink.matchSpectra(ds, cands)['path'] == near


def test_matchSpectra_genericFtt_byExperimentDir(tmp_path):
    root = _load(COMMENTED, tmp_path)
    ds = _ds3(root)
    ds.name = 'ftt'  # the generic NMRpipe data-set file convention

    a = tmp_path / 'yb-x-noesyn'
    a.mkdir()
    fa = str(a / 'ftt.ft3')
    makePipeFile(fa, (64, 32, 16))
    b = tmp_path / 'yb-x-cpmg'
    b.mkdir()
    fb = str(b / 'ftt.ft3')
    makePipeFile(fb, (64, 32, 16))

    cands = nefRelink.scanSpectrumFiles(str(tmp_path))
    # both stems are 'ftt' (100 tie); the dataset DIRECTORY carries the
    # identity: the experiment '15N NOESY-HSQC' token 'noesy' occurs in
    # 'yb-x-noesyn' (score 110) and in no other directory (100)
    assert nefRelink.matchSpectra(ds, cands)['path'] == fa


def test_matchSpectra_fusedExperimentLabel(tmp_path):
    root = _load(COMMENTED, tmp_path)
    ds = _ds3(root)
    ds.name = 'ftt_2'
    # fused experiment label ('sse' + 'noesyn' glued): a token-wise
    # comparison would never see the 'sse' / 'noesyn' parts
    ds.experiment.name = 'ssenoesyn'

    sse = tmp_path / 'yb-sse-noesyn'
    sse.mkdir()
    fsse = str(sse / 'ftt.ft3')
    makePipeFile(fsse, (64, 32, 16))
    ssd = tmp_path / 'yb-ssd-noesyn'
    ssd.mkdir()
    fssd = str(ssd / 'ftt.ft3')
    makePipeFile(fssd, (64, 32, 16))

    cands = nefRelink.scanSpectrumFiles(str(tmp_path))
    # stem 'ftt' is a subset of the name 'ftt_2'; in the dir score the
    # BOTH dir tokens 'sse' and 'noesyn' are contained in the fused
    # label 'ssenoesyn' (vs only 'noesyn' for 'yb-ssd-noesyn')
    assert nefRelink.matchSpectra(ds, cands)['path'] == fsse


# ---------------------------------------------------------------------------
# relinkSpectra end to end
# ---------------------------------------------------------------------------

def test_relinkSpectra_e2e(tmp_path):
    root = _load(COMMENTED, tmp_path)
    ds = _ds3(root)
    assert ds.dataStore is None
    # the NEF import left the importer's bogus point-count defaults on
    # the data dimensions (1280/2560)
    before = [dd.numPoints for dd in ds.sortedDataDims()]
    assert before != [64, 32, 16]

    base = tmp_path / 'yb-demo' / 'noesyn'
    base.mkdir(parents=True)
    makePipeFile(str(base / 'cnoesy1.ft3'), (64, 32, 16),
                 sf=(600.0, 150.0, 100.0), sw=(8000.0, 2000.0, 1200.0),
                 nuc=('1H', '15N', '13C'))
    (tmp_path / 'yb-demo' / 'README.txt').write_text('skip me')

    report = nefRelink.relinkSpectra(root, str(tmp_path / 'yb-demo'))
    assert report['candidates'] == 1
    assert [e['name'] for e in report['linked']] == ['cnoesy1']
    # the 15-dim dummy spectrum cannot be backed by an NMRpipe file
    assert [e['name'] for e in report['skipped']] == ['dummy15d']
    assert report['unlinked'] == []

    assert ds.dataStore is not None
    assert ds.dataStore.fileType == 'NmrPipe'
    assert ds.dataStore.headerSize == 2048
    assert ds.dataStore.fullPath == str(base / 'cnoesy1.ft3')
    assert os.path.exists(ds.dataStore.fullPath)
    assert tuple(ds.dataStore.numPoints) == (64, 32, 16)
    assert tuple(ds.dataStore.blockSizes) == (64, 1, 1)

    dims = ds.sortedDataDims()
    assert [d.numPoints for d in dims] == [64, 32, 16]
    assert [d.numPointsOrig for d in dims] == [64, 32, 16]
    assert [d.valuePerPoint for d in dims] == [
        8000.0 / 64, 2000.0 / 32, 1200.0 / 16]

    # second pass: everything already linked, nothing to do
    report2 = nefRelink.relinkSpectra(root, str(tmp_path / 'yb-demo'))
    assert report2['alreadyLinked'] == 1
    assert report2['linked'] == []
    assert report2['unlinked'] == []

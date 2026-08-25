"""S33 functional test: ccpnmr.exportNmrData spectrum + peak-list exporters.

Covers all 4 target formats (bruker / varian / unf / sparky) for both the
1D spectrum and the 1D peak list.  Verification mixes independent parsers
(text columns, a from-scratch UNF v2.1 binary reader) with nmrglue's own
readers (Varian fid round-trip, skipped when nmrglue is absent — it is an
optional dependency).
"""
import os
import struct

import numpy as np
import pytest

from ccpnmr import exportNmrData as x

N = 32
FREQ = np.linspace(499.5e6, 500.5e6, N)        # 500 MHz, 1 MHz wide
DATA = (np.sin(np.linspace(0, 6 * np.pi, N)) * 1e4 + 5e4).astype(np.float64)
PEAKS = [(FREQ[i], 20.0 + 10.0 * i) for i in range(0, N, 7)]


# ---------------------------------------------------------------------------
# independent readers
# ---------------------------------------------------------------------------

def _pairs(file, skip_prefixes=("#", "$", "##$")):
    """Two-column (x, y) parse of a text data file."""
    pts = []
    for line in open(file):
        line = line.strip()
        if not line or line.startswith(skip_prefixes) or line.startswith("##") or line.startswith("$$") or line.startswith("1.0 1.2"):
            continue
        t = line.split()
        if len(t) >= 2 and _isfloat(t[0]) and _isfloat(t[-1]):
            if line.startswith("##XY="):
                continue
            pts.append((float(t[0]), float(t[-1])))
    return np.array(pts) if pts else None


def _isfloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def _parse_bruker_f(path):
    pts = []
    xy = False
    for line in open(path):
        line = line.strip()
        if line.startswith("##XY="):
            xy = True
            continue
        if line.startswith("##") or line.startswith("$$"):
            continue
        if xy:
            t = line.split()
            pts.append((float(t[0]), float(t[1])))
    return np.array(pts)


def _parse_bruker_procs(path):
    d = {}
    for line in open(path):
        line = line.strip()
        if line.startswith("##$") and "=" in line:
            k, v = line[3:].split("=", 1)
            d[k] = v.strip().strip("<>")
    return d


def _parse_unf(path):
    """Minimal UNF v2.1 reader (header + records + float32 data)."""
    blob = open(path, "rb").read()
    assert blob[:4] == b"UNF ", "missing UNF magic"
    version, minor, reserved = struct.unpack_from("<HBB", blob, 4)
    assert version == 2, "expected UNF v2"
    flags, nrec = struct.unpack_from("<HH", blob, 8)
    rec_len, data_len = struct.unpack_from("<ii", blob, 12)
    dtype = blob[20]
    assert flags & 1 == 1, "records bit not set"
    records = {}
    off = 512
    for _ in range(nrec):
        rlen = struct.unpack_from("<i", blob, off)[0]
        p = off + 4
        end = off + rlen
        while p < end:
            nchars = struct.unpack_from("<i", blob, p)[0]
            p += 4
            if nchars == 0:
                break
            key = blob[p:p + nchars].decode("ascii").strip()
            p += ((nchars + 3) // 4) * 4
            nbytes = struct.unpack_from("<i", blob, p)[0]
            p += 4
            val = blob[p:p + nbytes].decode("ascii").strip()
            p += ((nbytes + 3) // 4) * 4
            records[key] = val
        off += rlen
    d = np.frombuffer(blob[off:off + data_len], dtype="<f4" if dtype == 1 else "<f8")
    return records, rec_len, data_len, d


# ---------------------------------------------------------------------------
# spectrum
# ---------------------------------------------------------------------------

class TestSpectrumExport:
    def test_bruker(self, tmp_path):
        d = str(tmp_path / "bruker")
        ret = x.export_spectrum_1d(DATA, FREQ, "bruker", d)
        assert os.path.basename(ret) == "f"
        pts = _parse_bruker_f(os.path.join(d, "f"))
        assert pts.shape == (N, 2)
        np.testing.assert_allclose(pts[:, 0], FREQ, rtol=2e-8, atol=1e-6)  # ~10 sig-fig text
        np.testing.assert_allclose(pts[:, 1], DATA, rtol=1e-6, atol=1e-3)
        procs = _parse_bruker_procs(os.path.join(d, "procs"))
        assert procs["NS"] == str(N) or float(procs["NS"]) == N
        assert abs(float(procs["SF"]) - 500.0) < 1e-9       # obs derived as centre (MHz)
        assert float(procs["NSC"]) == 0                      # real spectrum

    def test_varian_roundtrip(self, tmp_path):
        pytest.importorskip("nmrglue")
        from nmrglue import fileio
        d = str(tmp_path / "varian")
        ret = x.export_spectrum_1d(DATA, FREQ, "varian", d)
        assert os.path.basename(ret) == "fid"
        assert os.path.isfile(os.path.join(d, "procpar"))
        dic, ddata = fileio.varian.read(d)
        assert dic["S_SPEC"] == 1
        dd = np.asarray(ddata).ravel()
        assert dd.size >= N
        np.testing.assert_allclose(dd[:N].real, DATA, rtol=1e-5, atol=1e-2)

    def test_unf(self, tmp_path):
        f = str(tmp_path / "spec.unf")
        x.export_spectrum_1d(DATA, FREQ, "unf", f)
        recs, rec_len, data_len, d = _parse_unf(f)
        assert set(k for k in recs) >= {"1-UNF-LABELS", "1-UNF-CREATOR"}
        assert "HZ" in recs["1-UNF-LABELS"]
        assert data_len == 4 * N
        np.testing.assert_array_equal(d, DATA.astype(np.float32))
        tot = 512 + rec_len + data_len
        assert os.path.getsize(f) == tot

    def test_sparky(self, tmp_path):
        f = str(tmp_path / "spec_1d")
        x.export_spectrum_1d(DATA, FREQ, "sparky", f)
        pts = _pairs(f)
        assert pts.shape == (N, 2)
        np.testing.assert_allclose(pts[:, 0], FREQ, rtol=2e-8, atol=1e-6)  # ~10 sig-fig text
        np.testing.assert_allclose(pts[:, 1], DATA, rtol=1e-6, atol=1e-3)

    def test_bad_format_raises(self, tmp_path):
        with pytest.raises(ValueError):
            x.export_spectrum_1d(DATA, FREQ, "nope", str(tmp_path / "x"))

    def test_length_mismatch_raises(self, tmp_path):
        with pytest.raises(ValueError):
            x.export_spectrum_1d(DATA[:-1], FREQ, "sparky", str(tmp_path / "x"))

    def test_nonuniform_freq_raises(self, tmp_path):
        f = FREQ.copy()
        f[10] += 12345.0
        with pytest.raises(ValueError):
            x.export_spectrum_1d(DATA, f, "sparky", str(tmp_path / "x"))

    def test_explicit_obs_recorded(self, tmp_path):
        d = str(tmp_path / "bruker2")
        x.export_spectrum_1d(DATA, FREQ, "bruker", d, obs=470.0)
        procs = _parse_bruker_procs(os.path.join(d, "procs"))
        assert abs(float(procs["SF"]) - 470.0) < 1e-9


# ---------------------------------------------------------------------------
# peak list
# ---------------------------------------------------------------------------

class TestPeakListExport:
    def test_bruker(self, tmp_path):
        f = str(tmp_path / "list" / "peak")
        x.export_peak_list(PEAKS, "bruker", f)
        lines = [l for l in open(f).read().splitlines() if l.strip()]
        assert lines[0] == f"1.0 1.2 1 1 {len(PEAKS)}"
        for i, line in enumerate(lines[1:]):
            t = line.split()
            assert len(t) == 3
            assert int(t[0]) == i + 1
            assert abs(float(t[1]) - PEAKS[i][0]) < PEAKS[i][0] * 2e-8
            assert abs(float(t[2]) - PEAKS[i][1]) < 1e-9
        assert len(lines) == 1 + len(PEAKS)

    def test_varian(self, tmp_path):
        f = str(tmp_path / "list")
        x.export_peak_list(PEAKS, "varian", f)
        lines = [l for l in open(f).read().splitlines() if l.strip()]
        assert int(lines[0]) == len(PEAKS)
        for i, line in enumerate(lines[1:]):
            t = line.split()
            assert int(t[0]) == i + 1
            assert abs(float(t[1]) - PEAKS[i][0]) < PEAKS[i][0] * 2e-8
            assert abs(float(t[2]) - PEAKS[i][1]) < 1e-9

    def test_sparky(self, tmp_path):
        f = str(tmp_path / "peaks.lst")
        x.export_peak_list(PEAKS, "sparky", f)
        lines = [l for l in open(f).read().splitlines() if l.strip()]
        assert lines[0] == f"{len(PEAKS)} 1"
        for i, line in enumerate(lines[1:]):
            t = line.split()
            assert int(t[0]) == i + 1
            assert abs(float(t[1]) - PEAKS[i][0]) < PEAKS[i][0] * 2e-8
            assert abs(float(t[2]) - PEAKS[i][1]) < 1e-9
        assert len(lines) == 1 + len(PEAKS)

    def test_unf(self, tmp_path):
        f = str(tmp_path / "peaks.unf")
        x.export_peak_list(PEAKS, "unf", f)
        recs, rec_len, data_len, d = _parse_unf(f)
        assert "1-UNF-LABELS" in recs
        assert data_len == 8 * len(PEAKS)
        assert d.shape == (2 * len(PEAKS),)
        d2 = d.reshape(len(PEAKS), 2)
        np.testing.assert_array_equal(d2[:, 0], np.array([p[0] for p in PEAKS]).astype(np.float32))
        np.testing.assert_array_equal(d2[:, 1], np.array([p[1] for p in PEAKS]).astype(np.float32))

    def test_bad_format_raises(self, tmp_path):
        with pytest.raises(ValueError):
            x.export_peak_list(PEAKS, "nope", str(tmp_path / "x"))

    def test_bad_peaks_raise(self, tmp_path):
        with pytest.raises(ValueError):
            x.export_peak_list([(1.0,), (2.0, 3.0)], "sparky", str(tmp_path / "x"))
        with pytest.raises(ValueError):
            x.export_peak_list([], "sparky", str(tmp_path / "x"))
        with pytest.raises(ValueError):
            x.export_peak_list([("a", 1.0)], "sparky", str(tmp_path / "x"))

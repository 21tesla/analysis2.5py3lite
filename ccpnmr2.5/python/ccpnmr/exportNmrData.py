"""Basic NMR spectrum + peak-list export (Stage 33, 2026-08-25).

Self-contained exporters for the two file kinds the removed FormatConverter
app used to produce, written to the four common interop formats:

  * a 1D processed spectrum (frequency scale + real intensities)
  * a 1D peak list (frequency + intensity pairs)

Target formats
--------------
  brukser  (Bruker/TopSpin)
      spectrum : <dir>/f     JCAMP-DX ``##XY=`` file (F Hz, I) + <dir>/procs
      peak list: <file>      classic ``list/peak`` text (``1.0 1.2 1 1 N``)
  varian  (Varian/Agilent)
      spectrum : <dir>/fid   binary spectrum file (S_SPEC, I=0 pairs) +
                             <dir>/procpar   [via nmrglue]
      peak list: <file>      ``N`` / ``serial freq intensity`` text
  unf     (vendor-neutral)
      spectrum : <file>      UNF v2.1 binary, 1-D real (float32)
      peak list: <file>      UNF v2.1 binary, 2-D table (float32):
                             row = peak, col 0 = frequency (HZ),
                             col 1 = intensity (as given)
  sparky  (Sparky)
      spectrum : <file>      2-column text (F Hz, I)
      peak list: <file>      ``N 1`` / ``serial freq intensity`` text

Conventions
-----------
  * Frequencies are ABSOLUTE, in HZ.  Peak-list and spectrum columns are
    written as given (no rescaling).
  * ``obs`` is the observer frequency in MHz, optional; default is the
    centre of the frequency scale.  It is recorded where the format has a
    field for it (procs SF, procpar ob, UNF OBS).
  * ``path`` is a DIRECTORY for the bruker and varian spectrum outputs and
    a FILE for everything else.  Missing directories are created.
  * The 'varian' spectrum needs the optional nmrglue dependency
    (``pip install ccpnmr[export]`` / ``uv pip install 'nmrglue>=0.12'``);
    all other outputs are dependency-free.

This module is intentionally model-free: it takes plain arrays / peak pairs,
so it can be driven from scripts or from the Analysis app alike without
importing the GUI layer.
"""

from __future__ import annotations

import datetime
import os
import struct
import warnings

import numpy as np

__all__ = [
    "SPECTRUM_FORMATS",
    "PEAKLIST_FORMATS",
    "export_spectrum_1d",
    "export_peak_list",
]

SPECTRUM_FORMATS = ("bruker", "varian", "unf", "sparky")
PEAKLIST_FORMATS = ("bruker", "varian", "unf", "sparky")

_NMRGLUE_HINT = (
    "the '%s' spectrum output requires nmrglue; install the optional extra: "
    "uv pip install 'nmrglue>=0.12'  (or: pip install ccpnmr[export])"
)


# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------

def _check_1d(data, name):
    a = np.asarray(data, dtype=np.float64)
    if a.ndim != 1 or a.size == 0:
        raise ValueError(f"expected a non-empty 1D array for {name!r}")
    return a


def _scale_info(freq, obs):
    """Derive obs (MHz) / carrier (Hz) / sw (Hz) from the freq scale."""
    fmin, fmax = float(freq.min()), float(freq.max())
    carrier = (fmin + fmax) / 2.0
    if obs is None:
        obs = carrier / 1.0e6
    sw = max(fmax - fmin, 0.0)
    return float(obs), float(carrier), float(sw)


def _format_hz(v):
    """10 sig-fig decimal, no scientific notation (vendor text files)."""
    return format(float(v), ".10g")


# ---------------------------------------------------------------------------
# bruker (dependency-free)
# ---------------------------------------------------------------------------

def _write_bruker_dir(path, freq, intensities, obs=None):
    """Write <path>/f (JCAMP ##XY) + <path>/procs (classic TopSpin 1D set)."""
    freq = _check_1d(freq, "freq")
    intensities = _check_1d(intensities, "intensities")
    if freq.size != intensities.size:
        raise ValueError("freq and intensities must have the same length")
    obs_mhz, carrier_hz, sw_hz = _scale_info(freq, obs)

    os.makedirs(path, exist_ok=True)

    # -- procs (pure parameters: plain '##$KEY= value' text) --
    d = {
        "PROCNO": "0",
        "PPARMOD": 0,
        "NS": freq.size,
        "SF": obs_mhz,
        "CARF": carrier_hz,
        "SPW1": sw_hz,
        "P1": obs_mhz * 1.0e6,
        "NSC": 0,
        "PASE": 1,
        "ACQ1S": 0,
        "PROC1S": 1,
        "TUBE": "ccpnmr",
        "TITLE": "CCPNmr exporter: 1D spectrum",
    }
    with open(os.path.join(path, "procs"), "w") as f:
        f.write("$$ process parameters\n")
        for k in sorted(d):
            v = d[k]
            if isinstance(v, str):
                f.write(f"##${k}= <{v}>\n")
            else:
                f.write(f"##${k}= {_format_hz(v)}\n")

    # -- f (JCAMP-DX data file: ##Y= data lives after ##END=) --
    df = (float(freq[-1]) - float(freq[0])) / (freq.size - 1) if freq.size > 1 else 0.0
    with open(os.path.join(path, "f"), "w") as f:
        f.write("##TITLE=<f>\n")
        f.write("##TYPE=CALC\n")
        f.write(f"##NPOINTS={freq.size}\n")
        f.write("##FVAR=1\n")
        f.write(f"##SPVAR={freq.size}\n")
        f.write("##XUNITS=<HZ>\n")
        f.write("##YUNITS=<>\n")
        f.write(f"##XFIRST={_format_hz(freq[0])}\n")
        f.write(f"##XDINCR={_format_hz(df)}\n")
        f.write("##END=\n")
        f.write("##XY=\n")
        for fv, iv in zip(freq, intensities):
            f.write(f" {_format_hz(fv)}  {_format_hz(iv)}\n")
    return os.path.join(path, "f")


def _write_bruker_peaklist(path, peaks):
    """Write the classic 'list/peak' text (1D): header + 'serial f int'."""
    peaks = _check_peaks(peaks)
    _ensure_parent(path)
    with open(path, "w") as f:
        f.write(f"1.0 1.2 1 1 {len(peaks)}\n")
        for i, (fv, iv) in enumerate(peaks, 1):
            f.write(f"{i}  {_format_hz(fv)}  {_format_hz(iv)}\n")
    return path


# ---------------------------------------------------------------------------
# varian (nmrglue for the binary spectrum; text peak list is plain)
# ---------------------------------------------------------------------------

def _write_varian_dir(path, freq, intensities, obs=None):
    """Write <path>/fid (S_SPEC binary) + <path>/procpar via nmrglue."""
    freq = _check_1d(freq, "freq")
    intensities = _check_1d(intensities, "intensities")
    if freq.size != intensities.size:
        raise ValueError("freq and intensities must have the same length")
    obs_mhz, carrier_hz, sw_hz = _scale_info(freq, obs)
    n = freq.size

    try:
        from nmrglue import fileio
    except ImportError as e:
        raise ImportError(_NMRGLUE_HINT % "varian") from e
    pdic_param = getattr(fileio.varian, "create_pdic_param", None)
    if pdic_param is None:
        raise ImportError("installed nmrglue lacks fileio.varian; try 'nmrglue>=0.12'")

    os.makedirs(path, exist_ok=True)

    udic = fileio.fileiobase.create_blank_udic(1)
    udic[0].update(
        size=n,
        obs=obs_mhz,
        sw=sw_hz / 1.0e6,
        car=carrier_hz,
        label="1H",
        time=False,
        freq=True,
        complex=False,
    )
    dic = dict(fileio.varian.create_dic(udic))
    dic["S_SPEC"] = 1
    # dic2fileheader copies 'status' verbatim, so set the byte directly:
    # S_DATA(1) | S_SPEC(2) | S_FLOAT(8) | S_ACQPAR(64) = processed spectrum,
    # single-precision floats, parameters in the side file.
    dic["status"] = 1 + 2 + 8 + 64
    pp = dict(dic["procpar"])
    for key, val in (
        ("ob", [str(obs_mhz)]),
        ("sw", [str(sw_hz / 1.0e6)]),
        ("sfo1", [str(carrier_hz)]),
        ("unit", ["hertz"]),
        ("sdata", ["1"]),
        ("sspec", ["1"]),
        ("ntraces", ["1"]),
        ("scomplex", ["1"]),
        ("np", [str(n)]),
        ("np2", [str(n)]),
        ("ns", [str(n)]),
    ):
        pp[key] = pdic_param(key, val)
    dic["procpar"] = pp

    data32 = np.ascontiguousarray(intensities).astype(np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fileio.varian.write_fid(os.path.join(path, "fid"), dic, data32, overwrite=True)
        fileio.varian.write_procpar(os.path.join(path, "procpar"), pp)
    return os.path.join(path, "fid")


def _write_varian_peaklist(path, peaks):
    """Write the Varian 'list' text: line 1 = N, then 'serial f int'."""
    peaks = _check_peaks(peaks)
    _ensure_parent(path)
    with open(path, "w") as f:
        f.write(f"{len(peaks)}\n")
        for i, (fv, iv) in enumerate(peaks, 1):
            f.write(f"{i}  {_format_hz(fv)}  {_format_hz(iv)}\n")
    return path


# ---------------------------------------------------------------------------
# UNF v2.1 (little-endian binary, float32 data) — written by hand:
# no published nmrglue release ships a UNF module, so this minimal writer
# implements the subset of the spec needed for 1-D real arrays and 2-D
# peak tables (header / records / data block, 4-byte alignment).
# ---------------------------------------------------------------------------

_UNF_DATA_FLOAT32 = 1


def _unf_record(key, value):
    """One UNF v2.1 record: [rlen][nchars key nbytes value]* with (0,0) end.

    The length word covers the record as written (terminator + alignment pad
    included); key pads with spaces, value with zero bytes (spec).
    """
    key_b = str(key).encode("ascii")
    val_b = str(value).encode("ascii")
    if len(key_b) > 255 or len(val_b) > 255:
        raise ValueError("UNF key/value must be <=255 bytes")
    payload = struct.pack("<i", len(key_b)) + key_b + b" " * ((4 - len(key_b)) % 4)
    payload += struct.pack("<i", len(val_b)) + val_b + b"\x00" * ((4 - len(val_b)) % 4)
    payload += struct.pack("<ii", 0, 0)  # terminating pair
    rec_pad = (4 - (len(payload) + 4) % 4) % 4
    return struct.pack("<i", 4 + len(payload) + rec_pad) + payload + b"\x00" * rec_pad


def _unf_write(path, data, labels, extra_records):
    data = np.ascontiguousarray(np.asarray(data, dtype=np.float32)).astype("<f4")
    labels_lines = []
    for rec in extra_records:
        labels_lines.append(str(rec))
    records = {
        "1-CCPNMR-EXPORT": "exportNmrData (Stage 33)",
        "1-UNF-CREATOR": "ccpnmr/exportNmrData",
        "1-UNF-CREATION": datetime.date.today().isoformat(),
        "1-UNF-NAME": "1",
        "1-UNF-FILE": os.path.basename(path),
        "1-UNF-KEY": labels_lines[0] if labels_lines else "nmr data",
        "1-UNF-LABELS": labels,
    }
    nrec = len(records)
    blob = b""
    for k in sorted(records):
        blob += _unf_record(k, records[k])
    rec_len = len(blob)

    data_b = data.tobytes()
    # UNF v2.1 header (21 bytes, little-endian, zero-padded to 512):
    #   4s magic 'UNF ' | H version=2 | B minor=1 | B reserved=0
    #   H flags(1)=records | H nrecords | i record_length
    #   i data_length | B data_type (1 = single precision real)
    header = struct.pack("<4sHBBHHiiB", b"UNF ", 2, 1, 0, 1, nrec, rec_len, len(data_b), _UNF_DATA_FLOAT32)
    header += b"\x00" * (512 - len(header))

    _ensure_parent(path)
    with open(path, "wb") as f:
        f.write(header)
        f.write(blob)
        f.write(data_b)
    return path


def _write_unf_spectrum(path, freq, intensities, obs=None):
    """UNF v2.1, 1-D real float32, dim 0 unit HZ (PROC P)."""
    freq = _check_1d(freq, "freq")
    intensities = _check_1d(intensities, "intensities")
    if freq.size != intensities.size:
        raise ValueError("freq and intensities must have the same length")
    obs_mhz, _, _ = _scale_info(freq, obs)
    labels = (
        "[ { 1-UNF-DIMS: 1 }, "
        f"{{ 2-UNF-NAME: \"f1\", 2-UNF-UNIT: \"HZ\", 2-UNF-PROC: \"P\", "
        f"2-UNF-OBS: \"{_format_hz(obs_mhz)}\", 2-UNF-LABEL: \"frequency (Hz)\" }} ]"
    )
    return _unf_write(path, intensities, labels, ["1D spectrum"])


def _write_unf_peaklist(path, peaks):
    """UNF v2.1, 2-D (N,2) float32: col 0 = HZ, col 1 = intensity."""
    peaks = _check_peaks(peaks)
    arr = np.empty((len(peaks), 2), dtype=np.float32)
    arr[:, 0] = [p[0] for p in peaks]
    arr[:, 1] = [p[1] for p in peaks]
    labels = (
        "[ { 1-UNF-DIMS: 2 }, "
        "{ 2-UNF-NAME: \"peak\", 2-UNF-UNIT: \"DIMLESS\" }, "
        "{ 2-UNF-DIM0: 2, 2-UNF-NAME: \"value\", 2-UNF-UNIT: \"DIMLESS\" } ]"
    )
    extra = [
        "peak list",
        "column 0 = frequency (Hz); column 1 = intensity (as given)",
    ]
    return _unf_write(path, arr, labels, extra)


# ---------------------------------------------------------------------------
# sparky (text interop; nmrglue's 1D spectrum writer is unavailable)
# ---------------------------------------------------------------------------

def _write_sparky_spectrum(path, freq, intensities, obs=None):
    """2-column text (F Hz, I) with '#' comment header — Sparky 1D input."""
    freq = _check_1d(freq, "freq")
    intensities = _check_1d(intensities, "intensities")
    if freq.size != intensities.size:
        raise ValueError("freq and intensities must have the same length")
    obs_mhz, _, sw_hz = _scale_info(freq, obs)
    _ensure_parent(path)
    with open(path, "w") as f:
        f.write("# 1D spectrum (CCPNmr exporter, Stage 33)\n")
        f.write(f"# obs = {_format_hz(obs_mhz)} MHz, sw = {_format_hz(sw_hz)} Hz\n")
        f.write("# freq (Hz)  intensity\n")
        for fv, iv in zip(freq, intensities):
            f.write(f"{_format_hz(fv)}  {_format_hz(iv)}\n")
    return path


def _write_sparky_peaklist(path, peaks):
    """Sparky peak-list text: 'N ncols(=1)' then 'serial shift intensity'."""
    peaks = _check_peaks(peaks)
    _ensure_parent(path)
    with open(path, "w") as f:
        f.write(f"{len(peaks)} 1\n")
        for i, (fv, iv) in enumerate(peaks, 1):
            f.write(f"{i}  {_format_hz(fv)}  {_format_hz(iv)}\n")
    return path


# ---------------------------------------------------------------------------
# shared peak validation
# ---------------------------------------------------------------------------

def _check_peaks(peaks):
    try:
        out = []
        for p in peaks:
            fv, iv = float(p[0]), float(p[1])
            if not (np.isfinite(fv) and np.isfinite(iv)):
                raise ValueError(f"non-finite peak value: {p!r}")
            out.append((fv, iv))
    except (TypeError, ValueError, IndexError) as e:
        if isinstance(e, ValueError) and "non-finite" in str(e):
            raise
        raise ValueError(f"peaks must be an iterable of (freq_hz, intensity) pairs: {e}")
    if not out:
        raise ValueError("empty peak list")
    return out


def _ensure_parent(path):
    parent = os.path.dirname(os.path.abspath(path))
    if parent and not os.path.isdir(parent):
        os.makedirs(parent, exist_ok=True)


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------

def export_spectrum_1d(intensities, freq, fmt, path, obs=None):
    """Write a 1D processed spectrum in the requested format.

    Parameters
    ----------
    intensities : array_like
        Real intensities of the spectrum (as given; not rescaled).
    freq : array_like
        Absolute frequency (Hz) of each point, same length as ``intensities``.
        The scale must be evenly spaced (vendor text formats are XFIRST/XDINCR
        based); the writer verifies linearity up to float tolerance.
    fmt : {'bruker', 'varian', 'unf', 'sparky'}
    path : str
        DIRECTORY for 'bruker' (``f`` + ``procs``) and 'varian' (``fid`` +
        ``procpar``); FILE for 'unf' and 'sparky'.  Created if missing.
    obs : float, optional
        Observer frequency in MHz; default = centre of the freq scale.

    Returns
    -------
    str : path of the primary data file written.
    """
    fmt = str(fmt).lower()
    if fmt not in SPECTRUM_FORMATS:
        raise ValueError(f"unsupported spectrum format {fmt!r}; expected one of {SPECTRUM_FORMATS!r}")

    freq = _check_1d(freq, "freq")
    intensities = _check_1d(intensities, "intensities")
    if freq.size > 2:
        d = np.diff(freq)
        if not np.allclose(d, d[0], rtol=1e-6, atol=abs(d[0]) * 1e-6 + 1e-12):
            raise ValueError("freq must be evenly spaced (vendor formats store XFIRST/XDINCR)")

    if fmt == "bruker":
        return _write_bruker_dir(path, freq, intensities, obs)
    if fmt == "varian":
        return _write_varian_dir(path, freq, intensities, obs)
    if fmt == "unf":
        return _write_unf_spectrum(path, freq, intensities, obs)
    return _write_sparky_spectrum(path, freq, intensities, obs)


def export_peak_list(peaks, fmt, path, obs=None):
    """Write a 1D peak list (freq_hz, intensity pairs) in the given format.

    Parameters
    ----------
    peaks : iterable of (float, float)
        (absolute frequency Hz, intensity as given) per peak.
    fmt : {'bruker', 'varian', 'unf', 'sparky'}
    path : str
        Output FILE (e.g. ``out/list/peak``, ``out/list``, ``out/peaks.unf``,
        ``out/peaks.lst``).  Parent directories are created.
    obs : float, optional
        Observer frequency in MHz (currently recorded only where the record
        layout has room; included for API symmetry with the spectrum path).

    Returns
    -------
    str : the file path written.
    """
    fmt = str(fmt).lower()
    if fmt not in PEAKLIST_FORMATS:
        raise ValueError(f"unsupported peak-list format {fmt!r}; expected one of {PEAKLIST_FORMATS!r}")

    if fmt == "bruker":
        return _write_bruker_peaklist(path, peaks)
    if fmt == "varian":
        return _write_varian_peaklist(path, peaks)
    if fmt == "unf":
        return _write_unf_peaklist(path, peaks)
    return _write_sparky_peaklist(path, peaks)

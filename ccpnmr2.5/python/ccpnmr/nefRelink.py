"""Automatic relinking of spectrum data files for NEF-imported projects.

NEF v1.1 carries peaks, shifts, restraints and metadata - never raw
spectrum matrix data.  After "Load NEF..." every DataSource exists but
its ``dataStore`` (the file holding the matrix) is missing, and the GUI
warns "data file None not accessible" when a spectrum is rendered.

This module re-establishes those links against a directory tree of
spectrum data files - typically the directory the NEF file came from,
which is where the original data lives in the customary NMRpipe layout
(``yb-<sample>-<experiment>/`` dataset subdirectories, ``*.ft2``/``*.ft3``
files).  Public API:

- ``scanSpectrumFiles(baseDir)`` - walk ``baseDir`` (recursively) and
  read the 2048-byte NMRpipe header of every file it can parse;
- ``matchSpectra(dataSource, candidates)`` - pick the one file that
  belongs to a single DataSource (name, then dimension count, then
  dataset directory vs experiment name);
- ``relinkSpectra(memopsRoot, baseDir)`` - link every unlinked
  DataSource through ``ccpnmr.v2io.NefIo.addDataStore`` and restore
  each FreqDataDim's ``numPoints``/``numPointsOrig``/``valuePerPoint``
  from the file header (without ``point_count`` in the NEF the importer
  leaves bogus 1280/2560 defaults on the dimensions), returning a
  linked/unlinked report.

NMRpipe is the only header format parsed; the extension point for more
formats is the param-module list in ``ccp.format.spectra.OpenSpectrum``.
"""

import os
import re

from ccp.format.spectra.params import NmrPipeParams
from ccpnmr.v2io import NefIo

_TOKEN_SPLIT = re.compile(r'[^0-9a-zA-Z]+')

# NmrPipe headers express at most 4 dimensions (NmrPipeParams.parseFile);
# a DataSource with more (e.g. imported 15-dim dummy spectra) cannot be
# backed by an NMRpipe file.
_NMPIPE_MAX_DIM = 4

# spectrometer-frequency agreement tolerance (MHz)
_SF_TOLERANCE = 0.5


def _tokens(text):
    """Lower-case alphanumeric tokens of ``text`` (None-safe)."""
    if not text:
        return []
    return [t for t in _TOKEN_SPLIT.split(text.lower()) if t]


def _fileStem(fileName):
    """The file name with ALL trailing extensions stripped.

    'sswt-298K-hsqc-1016.ft2' -> 'sswt-298K-hsqc-1016'; 'ftt.ft3' -> 'ftt'.
    """
    stem = fileName
    while True:
        base, ext = os.path.splitext(stem)
        if not ext:
            return stem
        stem = base


def _sfOf(dataDim):
    """The spectrometer frequency of a data dimension (via its ExpDimRef),
    or None when unknowable."""
    expDim = dataDim.expDim
    if expDim is None:
        return None
    expDimRef = expDim.findFirstExpDimRef()
    if expDimRef is None:
        return None
    return expDimRef.sf


def _sfAgrees(dataSource, cand):
    """True when every header frequency agrees (within tolerance) with the
    DataSource's dimension frequencies."""
    for dataDim, sf in zip(dataSource.sortedDataDims(), cand['sf']):
        dsSf = _sfOf(dataDim)
        if dsSf is not None and abs(dsSf - sf) > _SF_TOLERANCE:
            return False
    return True


def _matchScore(dataSourceName, stem, relDir, experimentName, projectName):
    """Rank how well one candidate file fits one DataSource (0 = no case).

    Name score (a file must at least name-match the DataSource):
    - 100: the file stem is exactly the DataSource name;
    - 60+ : the name's tokens contain / are contained in the stem's;
    - 20+ : the names share tokens.
    Directory score (added): for generic data-set file names (the
    NMRpipe 'ftt' convention, where several files share the stem), the
    dataset DIRECTORY is where the identity lives - each pair of an
    experiment/project name token and a directory token (both len >= 3)
    that match by substring in EITHER direction adds 10; the either-
    direction form handles fused labels ('ssenoesyn' ~ 'sse-noesyn').
    Two exactly-tied names in sibling dataset dirs are thus resolved by
    which dir the experiment name points at.
    """
    name = (dataSourceName or '').lower()
    if not name:
        return 0
    if name == stem:
        nameScore = 100
    else:
        nameTok, stemTok = set(_tokens(name)), set(_tokens(stem))
        common = nameTok & stemTok
        if nameTok <= stemTok or stemTok <= nameTok:
            nameScore = 60 + 10 * len(common)
        elif common:
            nameScore = 20 + 10 * len(common)
        else:
            return 0
    dirScore = 0
    dirToks = {d for d in _tokens(relDir) if len(d) >= 3}
    for label in (experimentName, projectName):
        for tok in _tokens(label):
            if len(tok) < 3:
                continue
            for d in dirToks:
                if d == tok or tok in d or d in tok:
                    dirScore += 10
    return nameScore + dirScore


def scanSpectrumFiles(baseDir):
    """Recursively scan ``baseDir`` for files with a readable NMRpipe header.

    Returns a list (sorted by path) of dicts with keys ``path``, ``stem``,
    ``relDir``, ``ndim``, ``npts``, ``block``, ``sw``, ``sf``, ``nuc``,
    ``bigEndian``, ``head`` and ``nbytes``.  Files whose headers do not
    parse (anything that is not an NMRpipe data file) are skipped
    silently.
    """
    baseDir = os.path.abspath(baseDir)
    if not os.path.isdir(baseDir):
        raise ValueError(f'not a directory: {baseDir}')
    out = []
    for dirpath, _dirnames, filenames in os.walk(baseDir):
        for fileName in sorted(filenames):
            path = os.path.join(dirpath, fileName)
            try:
                params = NmrPipeParams.NmrPipeParams(path)
                ndim = params.ndim
                npts = list(params.npts[:ndim])
                if None in npts or any(p <= 0 for p in npts):
                    continue
            except Exception:
                continue
            out.append({
                'path': path,
                'stem': _fileStem(fileName),
                'relDir': os.path.relpath(dirpath, baseDir),
                'ndim': ndim,
                'npts': npts,
                'block': list(params.block[:ndim]),
                'sw': [params.sw[i] for i in range(ndim)],
                'sf': [params.sf[i] for i in range(ndim)],
                'nuc': [params.nuc[i] for i in range(ndim)],
                'bigEndian': bool(params.big_endian),
                'head': params.head,
                'nbytes': params.nbytes,
            })
    out.sort(key=lambda c: c['path'])
    return out


def matchSpectra(dataSource, candidates):
    """Pick the candidate file belonging to ``dataSource`` (or None).

    Hard gate: ``candidate ndim == dataSource.numDim``.  Candidates are
    ranked by ``_matchScore`` (file name vs data source name, then the
    dataset directory vs the experiment/project names for generic 'ftt'
    names) and, on equal rank, by spectrometer-frequency agreement.  The
    input list order (path-sorted) breaks any final tie, so the result
    is deterministic.
    """
    project = dataSource.root.findFirstNmrProject() if dataSource.root is not None else None
    experiment = dataSource.experiment
    expName = experiment.name if experiment is not None else ''
    projName = project.name if project is not None else ''

    best, bestKey = None, None
    for cand in candidates:
        if cand['ndim'] != dataSource.numDim:
            continue
        score = _matchScore(
            dataSource.name, cand['stem'], cand['relDir'], expName, projName)
        if score <= 0:
            continue
        key = (score, 1 if _sfAgrees(dataSource, cand) else 0)
        if bestKey is None or key > bestKey:
            best, bestKey = cand, key
    return best


def _restoreDataDims(dataSource, cand):
    """Point the DataSource's data dimensions at the header's geometry.

    The NEF importer creates FreqDataDims without ``point_count`` and
    leaves its 1280/2560 defaults in place; the renderers read
    ``numPoints`` from the dimension, so relinking must fix it.
    ``valuePerPoint`` follows the native convention (sw / numPoints, as
    written by ``ccp.util.Spectrum.createSpectrum``).
    """
    for i, dataDim in enumerate(dataSource.sortedDataDims()):
        if i >= len(cand['npts']):
            break
        numPoints = cand['npts'][i]
        dataDim.numPoints = numPoints
        dataDim.numPointsOrig = numPoints
        sw = cand['sw'][i]
        if sw:
            dataDim.valuePerPoint = sw / numPoints


def relinkSpectra(memopsRoot, baseDir):
    """Link every unlinked DataSource of ``memopsRoot`` to a spectrum file
    under ``baseDir``; return a report dict.

    Report keys: ``baseDir``, ``candidates`` (scan count), ``linked``
    (list of {experiment, name, file, numPoints}), ``unlinked`` (list of
    {experiment, name, numDim}), ``skipped`` (DataSources with more than
    4 dimensions, which NMRpipe headers cannot carry) and
    ``alreadyLinked`` (count).  Each candidate file is consumed at most
    once.
    """
    baseDir = os.path.abspath(baseDir)
    candidates = scanSpectrumFiles(baseDir)
    report = {
        'baseDir': baseDir,
        'candidates': len(candidates),
        'linked': [],
        'unlinked': [],
        'skipped': [],
        'alreadyLinked': 0,
    }
    nmrProject = memopsRoot.findFirstNmrProject()
    if nmrProject is None:
        raise ValueError('no NmrProject in the given memops root')
    remaining = list(candidates)
    for experiment in nmrProject.sortedExperiments():
        for dataSource in experiment.sortedDataSources():
            expName = experiment.name
            if dataSource.dataStore is not None:
                report['alreadyLinked'] += 1
                continue
            if dataSource.numDim > _NMPIPE_MAX_DIM:
                report['skipped'].append({
                    'experiment': expName,
                    'name': dataSource.name,
                    'numDim': dataSource.numDim,
                })
                continue
            cand = matchSpectra(dataSource, remaining)
            if cand is None:
                report['unlinked'].append({
                    'experiment': expName,
                    'name': dataSource.name,
                    'numDim': dataSource.numDim,
                })
                continue
            remaining.remove(cand)
            NefIo.addDataStore(
                dataSource,
                cand['path'],
                numPoints=cand['npts'],
                blockSizes=cand['block'],
                isBigEndian=cand['bigEndian'],
                numberType='float',
                headerSize=cand['head'],
                nByte=cand['nbytes'],
                fileType='NmrPipe',
            )
            _restoreDataDims(dataSource, cand)
            report['linked'].append({
                'experiment': expName,
                'name': dataSource.name,
                'file': cand['path'],
                'numPoints': cand['npts'],
            })
    return report

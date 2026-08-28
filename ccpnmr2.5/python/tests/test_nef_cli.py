"""
Stage-38 app-wiring tests: the ``ccpnmr-nef`` command-line entry point
(``ccpnmr/nefCli.py``) and the import wrapper
``ccp.gui.Io.loadNefProject`` (formerly behind the Project menu's
"Load NEF..." item, which has since been removed).

Covers (bundled testdata, headless, network-independent):

- ``import`` creates a new, SAVED CCPN project directory in the current
  working directory (the memops XML is on disk, so the export subcommand
  and the GUI's Project > Open Project can both open it);
- ``export`` writes a NEF v1.1 file for that project and re-importing the
  exported file reproduces the object tree (count-based round trip on the
  Commented example: 235 residues / 2 data sources / 6 peaks / 108 shifts
  / 4 constraint lists);
- import into an already-existing project directory fails without
  ``--force`` and (with the flag) removes it and re-imports;
- ``loadNefProject`` (the GUI wrapper) returns the in-memory project root
  without saving — mirroring the "Open Project" flow, where the user
  saves explicitly afterwards.

The read/write semantics themselves are pinned down by
``test_v2io_nef.py`` (import) and ``test_nef_export.py`` (export +
value-level round trip); these tests only re-check the counts.
"""
import os
import shutil
from contextlib import redirect_stdout
from io import StringIO

import ccpnmr.nef.NefImporter as NefImporterModule
from ccpnmr import nefCli
from ccpnmr.analysis.core import AssignmentBasic
from ccpnmr.v2io import NefIo
from memops.general import Io as memopsIo
from tests.test_nef_relink import makePipeFile

NEF_DIR = os.path.join(os.path.dirname(NefImporterModule.__file__))
TESTDATA = os.path.join(NEF_DIR, 'testdata')
COMMENTED = os.path.join(TESTDATA, 'CCPN_Commented_Example.nef')


def _counts(root):
    """Summarise the object tree (mirrors test_v2io_nef._counts)."""
    nmr = root.findFirstNmrProject()
    ms = root.findFirstMolSystem()
    expts = nmr.sortedExperiments()
    dataSources = [d for e in expts for d in e.sortedDataSources()]
    peaks = [p for d in dataSources for pl in d.getPeakLists() for p in pl.getPeaks()]
    shiftLists = AssignmentBasic.getShiftLists(nmr)
    shifts = [m for sl in shiftLists for m in sl.getMeasurements()]
    constraintLists = [
        cl
        for cs in sorted(list(nmr.getNmrConstraintStores()), key=lambda cs: cs.serial)
        for cl in sorted(list(cs.getConstraintLists()), key=lambda cl: (cl.className, cl.name))
    ]
    return {
        'residues': sum(len(c.sortedResidues()) for c in ms.getChains()),
        'dataSources': dataSources,
        'peaks': peaks,
        'shifts': shifts,
        'constraintLists': constraintLists,
    }


def _cli(*argv):
    """Run the console entry point headless; return its exit code."""
    with redirect_stdout(StringIO()):
        return nefCli.main(list(argv))


# ---------------------------------------------------------------------------
# import
# ---------------------------------------------------------------------------

def test_cli_import_creates_saved_project(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rc = _cli('import', COMMENTED, '--project-name', 'clit')
    assert rc == 0
    projDir = tmp_path / 'clit'
    assert projDir.is_dir()
    # the save happened: the persistence XML is on disk
    assert (projDir / 'memops' / 'Implementation' / 'clit.xml').is_file()
    with redirect_stdout(StringIO()):
        root = memopsIo.loadProject(str(projDir), projectName='clit')
    c = _counts(root)
    assert c['residues'] == 235
    assert len(c['dataSources']) == 2
    assert len(c['peaks']) == 6
    assert len(c['shifts']) == 108
    assert len(c['constraintLists']) == 4


def test_cli_import_existing_dir_requires_force(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    assert _cli('import', COMMENTED, '--project-name', 'clit') == 0
    # same directory name, without --force: refuse
    assert _cli('import', COMMENTED, '--project-name', 'clit') == 1
    assert 'error' in capsys.readouterr().err
    # with --force: the old project is deleted and re-imported
    assert _cli('import', COMMENTED, '--project-name', 'clit', '--force') == 0
    assert (tmp_path / 'clit' / 'memops' / 'Implementation' / 'clit.xml').is_file()


# ---------------------------------------------------------------------------
# export + round trip
# ---------------------------------------------------------------------------

def test_cli_import_export_round_trip(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert _cli('import', COMMENTED, '--project-name', 'clit') == 0
    outNef = tmp_path / 'clit_out.nef'
    assert _cli('export', str(tmp_path / 'clit'), str(outNef)) == 0
    assert outNef.is_file()
    assert outNef.stat().st_size > 0

    # re-import the exported file into a fresh project and compare counts
    root = memopsIo.newProject('reimport', path=str(tmp_path), removeExisting=True)
    with redirect_stdout(StringIO()):
        root = NefIo.loadNefFile(str(outNef), memopsRoot=root)
    c = _counts(root)
    assert c['residues'] == 235
    assert len(c['dataSources']) == 2
    assert len(c['peaks']) == 6
    assert len(c['shifts']) == 108
    assert len(c['constraintLists']) == 4


def test_cli_export_missing_project_fails(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    assert _cli('export', str(tmp_path / 'no_such_project'), 'out.nef') == 1
    assert 'error' in capsys.readouterr().err
    assert not (tmp_path / 'out.nef').exists()


# ---------------------------------------------------------------------------
# Import wrapper (formerly the Project menu "Load NEF..." item)
# ---------------------------------------------------------------------------

def test_gui_loadNefProject_wrapper(tmp_path, monkeypatch):
    from ccp.gui.Io import loadNefProject

    monkeypatch.chdir(tmp_path)
    with redirect_stdout(StringIO()):
        root = loadNefProject(None, COMMENTED, projectName='gui_nef')
    # in-memory project, not yet saved (the GUI saves via Project > Save,
    # like the "Open Project" flow)
    assert root is not None
    assert not (tmp_path / 'gui_nef' / 'memops' / 'Implementation' / 'gui_nef.xml').exists()
    c = _counts(root)
    assert c['residues'] == 235
    assert len(c['peaks']) == 6
    assert len(c['shifts']) == 108


# ---------------------------------------------------------------------------
# import --relink (spectrum relink, Project menu "Import NEF + relink...")
# ---------------------------------------------------------------------------

def _importedDataSources(projDir, projectName):
    """Load the imported (saved) project; {numDim: DataSource}."""
    with redirect_stdout(StringIO()):
        root = memopsIo.loadProject(str(projDir), projectName=projectName)
    out = {}
    for expt in root.findFirstNmrProject().sortedExperiments():
        for ds in expt.sortedDataSources():
            out[ds.numDim] = ds
    return out


def test_cli_import_relink_default_dir(tmp_path, monkeypatch):
    # the NEF and its dataset dirs in ONE directory (the user's ~/NMR
    # layout): --relink with no value scans the NEF file's directory
    nmrLayout = tmp_path / 'nmr'
    nmrLayout.mkdir()
    nefFile = nmrLayout / 'commented.nef'
    shutil.copyfile(COMMENTED, nefFile)
    dset = nmrLayout / 'yb-demo' / 'noesyn'
    dset.mkdir(parents=True)
    makePipeFile(str(dset / 'cnoesy1.ft3'), (64, 32, 16),
                 sf=(600.0, 150.0, 100.0), sw=(8000.0, 2000.0, 1200.0),
                 nuc=('1H', '15N', '13C'))
    monkeypatch.chdir(tmp_path)
    with redirect_stdout(StringIO()) as log:
        rc = nefCli.main(['import', str(nefFile), '--project-name', 'rel',
                          '--relink'])
    assert rc == 0
    # 1 of the 2 unlinked spectra links; the 15-dim dummy is skipped
    assert 'linked 1 of 2' in log.getvalue()
    assert 'not NMRpipe' in log.getvalue()

    projDir = tmp_path / 'rel'
    assert (projDir / 'memops' / 'Implementation' / 'rel.xml').is_file()
    byDim = _importedDataSources(projDir, 'rel')
    assert byDim[3].dataStore is not None
    assert byDim[3].dataStore.fileType == 'NmrPipe'
    assert tuple(byDim[3].dataStore.numPoints) == (64, 32, 16)
    assert os.path.exists(byDim[3].dataStore.fullPath)
    assert byDim[15].dataStore is None  # >4 dims: not NMRpipe, left off


def test_cli_import_relink_explicit_dir(tmp_path, monkeypatch):
    # data living NEXT TO the NEF directory: an explicit --relink DIR
    nefDir = tmp_path / 'nefs'
    nefDir.mkdir()
    nefFile = nefDir / 'commented.nef'
    shutil.copyfile(COMMENTED, nefFile)
    dataRoot = tmp_path / 'elsewhere' / 'yb-demo'
    dataRoot.mkdir(parents=True)
    makePipeFile(str(dataRoot / 'cnoesy1.ft3'), (64, 32, 16))
    monkeypatch.chdir(tmp_path)
    with redirect_stdout(StringIO()) as log:
        rc = nefCli.main(['import', str(nefFile), '--project-name', 'rel2',
                          '--relink', str(tmp_path / 'elsewhere')])
    assert rc == 0
    assert 'linked 1 of 2' in log.getvalue()

    byDim = _importedDataSources(tmp_path / 'rel2', 'rel2')
    assert byDim[3].dataStore is not None
    assert byDim[3].dataStore.fullPath.startswith(str(tmp_path / 'elsewhere'))


def test_cli_import_relink_no_files(tmp_path, monkeypatch):
    # no spectrum files anywhere under the NEF's directory: the import
    # still succeeds, nothing links, the summary says so
    monkeypatch.chdir(tmp_path)
    with redirect_stdout(StringIO()) as log:
        rc = nefCli.main(['import', COMMENTED, '--project-name', 'rel3',
                          '--relink'])
    assert rc == 0
    assert 'no matching spectrum file found' in log.getvalue()

    byDim = _importedDataSources(tmp_path / 'rel3', 'rel3')
    assert byDim[3].dataStore is None
    assert byDim[15].dataStore is None

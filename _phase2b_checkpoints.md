# Phase 2b Checkpoints — fix remaining import-time failures
# (resume-safe map; update at every session boundary)

## Status (2026-08-21, FINAL — ready to commit)
- Baseline (Phase 2a commit `e22a93f`): import_smoke **1561 OK / 165 FAILED, 81 groups**.
- Final (Phase 2b): import_smoke **1634 OK / 92 FAILED** (+73 modules importable).
- Whole-tree compile: **0 syntax errors** (`python -m compileall -q ccpnmr2.5/python/`).
- Pytest (testpaths `tests/` + `nef/testing/`): **15 pass / 14 skip / 10 fail** — all 10
  failures are the known missing `/home/logan/software/testdata/` data (unchanged baseline).
- **124 phase .py files modified** + checkpoint/recipe/helper files (see staging list below).
- `.qwen/settings.json` + `database.txt` + `dbTable-new` + `dbTable.new` + `fooprof` =
  UNRELATED user files, EXCLUDED from staging (per workflow rule).
- Ruff: 152k pre-existing style findings in touched auto-generated API files — pre-existing,
  ruff `--fix` is Phase 3 scope, NOT a 2b gate.
- COMMIT STATE: pending single phase-boundary commit (message in "Decision record").

## Authoritative verification (run before committing)
```bash
cd /home/logan/software/ccpnmr2.5.2-qwen
# Primary metric (whole tree, groups failures by root cause):
MPLBACKEND=Agg .venv/bin/python import_smoke.py
# Syntax:
.venv/bin/python -m compileall -q ccpnmr2.5/python/ 2>&1 | grep -ci "Sorry"   # expect 0
# Per-module probe:
PYTHONPATH=ccpnmr2.5/python MPLBACKEND=Agg .venv/bin/python -c "import <a.b.c>"
```

## What Phase 2b converted (1561 -> 1634)
Same methodology as Phase 1e (explicit relative / correct-source imports,
removed-stdlib mappings, py2-idiom runtime fixes). Touched across
`ccp/api` (38), `cing/Scripts` (23), `utrecht/haddock` (7), `cing/NRG` (7),
`ccpnmr/api` (6), `cing/Database` (5), and ~20 more files in
`cambridge`, `ccpnmr`, `memops`, `molsim`, `ccp`.
Helper scripts (repo root, untracked; keep for resume): `_probe2b.py` (list modules
failing w/ a needle), `_fix2b_cing.py` (idempotent `import cing` inserter, used on
first 20 cing files).
Also applied this session: `cing/Scripts/CASD/casd3.py` import-source fix
(`NTpath` lives in `cing.Libs.NTutils`, not `cing.Libs.disk`).

## Remaining 92 failures — classified (survey complete 2026-08-21)
The 92 are NOT "py2 files left to convert". Breakdown:

### A. External / missing packages (not in tree, not py3-installable-or-appropriate) — ~28
- yasara (7 files) — commercial YASARA plugin.
- Refine (2) `xplor.py`, test_xplor.py — Xplor/Refine tool sub-repo, absent.
- protocol (2) anneal2.py, cingTest.py — XplorNIH `protocol` module, absent.
- UtilsAnalysis (1) Analysis/mouseBuffer.py — CING analysis sub-repo, absent.
- pdbe2 (1), pdbe.analysis (1) — `pdbe` analysis/adatah sub-repos, absent.
- memops.scripts (1) setLicenses.py — `memops/scripts` sub-package, absent.
- ccpncore (1) V2Upgrade.py — internal CCPN v2 core, absent.
- pymc (1) PeakSeparatorPyMC.py — PyMC2 (py2-only).
- sans (1) bmrb.py — py2 SOAP stack.
- pymol (3 + 1 NameError) — open-source pymol 2.x available on PyPI but old cmd API;
  install optional, verify API before adding.
- cherrypy / pycurl / decorator — optional install candidates from recipe (webFc, toposcmd, TestNefIo).

### B. IPython interactive / one-off data-batch scripts (module-level work on
### context vars `p`/`project`/`pTree`/`m` or hardcoded data paths) — ~25
- NameError `p` (2) analyzeCb2.py (`project = p`), contactDifference.py.
- NameError `project` (1) doValidateiCing.py.
- NameError `cgiDir` (2) doAnnotateNrgCing.py, mergeNrgBmrbShifts.py.
- NameError `date2num` (1) mouseBuffer5.py (matplotlib.dates, interactive).
- NameError `pymol` (1) pyMolWorks.py.
- NameError `refineParameters` (1) parametersTest.py.
- SystemExit(1) (3) dbTableUpdate290908.py, CASP/casp.py, interactive/mouseBuffer.py.
- TypeError `positions should be an iterable of numbers` (1) mouseBuffer2.py.
- ValueError `too many values to unpack` (1) mouseBuffer3.py.
- IndexError `list index out of range` (2) runQueenyEntry.py, pdbj_mine.py.
- TypeError `... not NoneType` (7) NRG/Casd* + doAnnotate* (root path = None; data-batch).
- AttributeError `importNameDefs` not found (1) INITIAL_SCRIPTS/mkresidueDefs.py
  (method signature drift; data-batch).
- FileNotFoundError macOS `//Library/WebServer/...` (5), `/Users/jd/...` (2),
  `/Users/wim/...` (1) — hardcoded per-machine data paths.
- FileNotFoundError `/home/logan/software/ccpnmr2.5.2-qwen/ccpnmr2.5/data/PluginCode...` (2)
  CING_paper_queries.py, d1d2plot.py — data files not in repo.
- HTTPError 404 (2) getRCSB_PDB.py, RESTfulExample.py — live network at import.
- noseTestCing.py (`None.run()`) — nose test-runner entry (nose uninstalled).

### C. Environment-config gated (need env vars / external binaries) — ~6
- KeyError `ISD_ROOT` (1) isd_project_template.py.
- KeyError `CYANA` (1) addCYANA2.py.
- KeyError `ARG+` (1) addSHIFTS.py.
- Exception `CASD_HOME not set` (4) nijmegen/CASD/*.py.
- CASD/Meccano C ext (1) grenoble/meccano — C ext not built (out of 2a scope).

### D. Custom ImportWarning plugin wrappers (raise when external tool absent) — 8
Wattos (2), dssp (2), Yasara (1), Vasco (1), nih (1), procheck (1), shiftx (1),
Molgrap (1) — all `cing/PluginCode/*.py`; by-design raise ImportWarning.

### E. Cython/C extension not compiled (build task) — 6 + 1
- cython.superpose (6) test_vector.py, queeny.py, convertChi1Chi22Db.py,
  convertD1D2_2Db2.py, convertPhiPsi2Db.py, core/validate.py — needs `superpose.pyx`
  compiled to `.so` (Cython + C build; like Phase 2a/1b extensions).
- Meccano (1) grenoble — see C.

### F. API not present in this distribution (missing class/member) — 6
- NmrSimRunFrame (1) CingFrame.py — `EditCalculation.py` has NmrCalcRunFrame only.
- Coplanar/CoplanarList (1) x3dna.py — class nowhere in tree.
- ProjectTree.openCompleteTree (2) casd2.py, nmr_redo_compareProjects.py.
- `NoneType.molecule` (2) mouseBuffer6.py, printResonances.py (Project stub has
  no molecule).
- `NoneType.allResidues` (1) compareShifts.py; `NoneType.run` (1) noseTestCing.py.
- SkipTest (1) Libs/test/test_Imagery.py; SystemExit usage (1) cython/compile.py
  (vendored Cython/setuptools tool, `usage()` at import).

## Deferred (NEEDS user decision — involves env changes)
1. Install optional external deps (pymol/cherrypy/pycurl/decorator) — could lift ~5-10
   modules but adds heavy deps; verify API compatibility first.
2. Build Cython `superpose` ext (`uv add cython`, compile `.pyx`) — lifts 6 modules.
3. Stub missing CING API (Project.molecule, ProjectTree.openCompleteTree,
   Coplanar) — would lift ~5 modules but fabricates API surface.

## Phase 2c — full-attribute build (user directive 2026-08-22: work through each
## bucket so the build retains the software's current attributes; checkpoint after
## EACH bucket)

### Bucket 1 — optional third-party deps — ✅ DONE (docs only, no code change)
Verified venv state: **cherrypy 18.10.0, decorator 5.3.1, pycurl 7.47.0, cython 3.2.9
already installed** (installed during the interrupted 2b session; their importers
`webFc.py`/`TestNefIo.py`/`toposcmd.py`/`PoolDownloader.py` all import OK now).
- **pymol: classified EXTERNAL (deliberate, matches original attrs).**
  - Original distribution did NOT bundle pymol (absent from `environment_Linux.yml`;
    scripts say "Requires the pymol python code etc to be properly installed").
  - `pyMolWorks.py` runs `pymol.finish_launching()` at MODULE level (PyMOL C++ engine
    start) — installing pymol 2.x in this headless env would hang/hard-fail the smoke
    run instead of failing fast; 1.x-era API era mismatch.
- **py2-only / internal sub-repos (NOT pip-installable, keep as-is)**: pymc (PyMC2),
  sans, yasara, ccpncore, `Refine`, `protocol`, `UtilsAnalysis`, `pdbe.analysis`,
  `pdbe2`, `memops.scripts`.
- **Smoke delta: 0 (bucket is classification; importers that COULD be lifted are already
  OK from the 2b session work).**

### Bucket 2 — Cython `superpose` ext build — IN PROGRESS
(6 importers: queeny.py, core/validate.py, convertChi1Chi22Db.py, convertD1D2_2Db2.py,
convertPhiPsi2Db.py, test/test_vector.py). cython 3.2.9 available.

## Decision record (2026-08-21)
**User decision: FINALIZE + COMMIT (scope = core only).** Deferred items above are left
for future sessions (explicit follow-ups: optional deps, superpose build, CING stubs).

### Staging list (commit exactly these, never `git add .`)
- 124 modified `.py` files under `ccpnmr2.5/python/` (see `git diff --name-only` at commit time)
- NEW: `_phase2b_checkpoints.md`, `_phase2b_recipe.md` (tracked, matching `_phase2a_*` precedent)
- NEW: `_probe2b.py`, `_fix2b_cing.py` (phase helpers; tracked, matching `fix_has_key.py`
  / `import_smoke.py` / `migrate_syntax.py` precedent)
- EXCLUDE (user files / scratch): `.qwen/settings.json`, `database.txt`, `dbTable-new`,
  `dbTable.new`, `fooprof`

### Commit message
`Fix remaining import-time failures under Python 3.13 (Phase 2b) (smoke 1561->1634)`

### Post-commit
- [x] Update project memory (project/ccpnmr-modernization.md) Phase 2b section.
- Next session: Phase 3 (ruff --fix safe passes, pytest expansion, testdata-gated tests,
  optional: superpose/Cython build + optional deps + missing-CING-API stubs if desired).

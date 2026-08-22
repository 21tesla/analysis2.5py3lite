# Phase 4 Checkpoints — functional verification + distribution readiness
# (resume-safe map; update at every session boundary)

## Mission (user, 2026-08-21)
"Ensure the software will function. Goal: a distribution in py3 that others can use."
- Legacy modules (clouds, haddock, …) **KEEP in the distribution, EXCLUDED from Phase-4
  functional scope/tests** (user decision — not used by the field anymore).
- **GUI launch verified under Xvfb** (user approved installing `xvfb` via apt).
- **C extensions built at pip-install time** (C compiler + numpy headers in build-env;
  Meccano/GSL OPTIONAL with clear fallback). One install → full functionality on Linux.

## Baseline (start of Phase 4, after Phase 3 close-out `b5c81a4`)
- smoke **1643 OK / 0 FAILED / 83 BY-DESIGN**; compile 0 errors;
  pytest 15 pass / 14 skip / 10 fail (10 = data-gated `/home/logan/software/testdata/` missing, pre-existing).
- Phase 1–3 complete: syntax, runtime, C-ext (30 exts), import surface, ruff-safe + F821 audit.
- **After P4-2 (2026-08-21 `4fba0f9`):** smoke **1646 / 0 / 83**; pytest **43 pass** / 14 skip / 10 fail (10 = same data-gated); **6 py3 runtime bugs fixed**.

## Phase 4 scope decisions (recorded 2026-08-21)
1. Legacy modules: keep in dist, skip in functional tests.
2. GUI tests: Xvfb headless launch of each app.
3. C exts: built into the wheel at install time; Meccano optional (GSL).
venv python = Anaconda 3.13.5 (has Tk 8.6 ✅). `setup.py` has full `ext_modules` (incl Meccano FAM entry, GSL via `CCP_GSL_PREFIX`).

## Buckets

### P4-1 — Entry points: 8 working console scripts — status: DONE (commit pending)
ADDED `def main(argv=None):` to the 7 modules that lacked it (wrapping each existing `__main__`
logic verbatim into main + `if __name__: main()` guard): ccpnmr/eci/EntryCompletionGui.py,
cambridge/dangle/DangleGui.py, ccpnmr/format/gui/DataShifter.py (inline block),
pdbe/deposition/dataFileImport/dataFileImportGui.py, extendNmr/ExtendNmrGui.py,
ccpnmr/format/gui/FormatConverter.py (inline block), ccpnmr/update/UpdateAuto.py.
`ccpnmr.analysis.AnalysisGui` already had main().
**Verify:** all 7 compile; **8/8 entry points resolvable** (importlib + signature check); gate below.
Discovered during verify: Gl/TkHandler broken (stale py2 .so) → see P4-4a.
System deps: `libglut.so.3` missing here (freeglut 3.12 has different SONAME) — user-local shim
`~/local-libs/libglut.so.3`; `xvfb` NOT installed (sudo needs password — use conda-forge env instead).
Only `ccpnmr.analysis.AnalysisGui:main` exists (`def main`, line 84). 7 modules have a
`launch*()`/inline `__main__` block but NO `main()`:
| console script | target | current `__main__` logic |
|---|---|---|
| ccpnmr | ccpnmr.analysis.AnalysisGui:main | ✅ `main(projectDir, max_size, glDirect)` exists |
| ccpnmr-eci | ccpnmr.eci.EntryCompletionGui:main | `launchApplication(projectDir)` (argv[1] opt) |
| ccpnmr-dangle | cambridge.dangle.DangleGui:main | `launchDangle(filename)` (argv[1] opt) |
| ccpnmr-data-shifter | ccpnmr.format.gui.DataShifter:main | INLINE block: load projects → Tk → mainloop |
| ccpnmr-deposition | pdbe.deposition.dataFileImport.dataFileImportGui:main | `launchDataFileImport()` |
| ccpnmr-extend-nmr | extendNmr.ExtendNmrGui:main | `launchApplication(projectDir)` (argv[1] opt) |
| ccpnmr-format-converter | ccpnmr.format.gui.FormatConverter:main | INLINE block: root/initProject → mainloop |
| ccpnmr-update | ccpnmr.update.UpdateAuto:main | `updateAuto()` + `os._exit(0)` |
- Pattern: add `def main(argv=None):` wrapping the existing `__main__` logic verbatim;
  `if __name__ == "__main__": main()`. Preserve argv semantics per module.
- Verify: `pip install -e . --no-deps` in .venv → 8 scripts exist; `hasattr`/signature check;
  GUI ones boot-tested in P4-5.

### P4-4a — Migrate + build the 7 remaining stale-py2 C exts — status: DONE (commit pending)
**RESULT 7/7 import OK + functional smoke OK + the 6 GUI modules that used to print
"Will continue without ... C functionality" now import CLEAN (no degraded path):**
ccpnmr.analysis.frames.WindowFrame, ccpnmr.analysis.Analysis, ccp.gui.ViewStructureFrame,
ccp.gui.ViewChemCompVarFrame, ccpnmr.analysis.core.WindowDraw, ccpnmr.analysis.core.DataAnalysisBasic.
Functional: `StructStructure()` instance API (addAtom/addBond/draw/.../zoom) OK;
`PeakCluster(3,0)`/`(3,2)` OK (addPeak/draw/removePeak/...); ContourFile exposes ContourFile+StoredContourFile.

**Migrated (recipe = `_phase2a_recipe.md`, verbatim):** py_gl_handler.c, py_tk_handler.c,
py_structure.c (PLUS 2 stale `Py_FindMethod(methods, obj, "remove")` → `PyObject_CallMethod(obj,"remove","O",obj)`),
py_win_peak_list.c, py_peak_cluster.c, py_contour_file.c (real getattr dispatch → PyUnicode_AsUTF8+strcmp),
py_slice_file.c (ditto). ALSO py_tk_util.c: `PyString_AsString` → `(char*)PyUnicode_AsUTF8`.
setup.py: +7 FAM entries; new infra: `TKINC` (env CCP_TK_PREFIX), `X11LINK=-l:libX11.so.6`,
`DRAWDEPS`/`DRAWLIBS` (py_draw_handler pulls the whole handler chain gl/tk/pdf/ps/store + cores;
analysis exts statically embed them as the stale .so's did). ContourFile also needs `{G}/contourer.c`.
BUILD: `CC=/usr/bin/gcc CXX=/usr/bin/g++ CCP_EXT=<7 names> .venv/bin/python setup.py build_ext --inplace --force`
(--force REQUIRED when the source list changes; system gcc because anaconda gcc lacks GL/tk in its sysroot;
tk.h/tcl.h from anaconda include). Copy flat .so (`ccpnmr2.5/python/<Name>.cpython-313*.so`) onto c-tree targets:
`c/memops/global/{Gl,Tk}Handler.so`, `c/ccp/structure/StructStructure.so`,
`c/ccpnmr/analysis/{WinPeakList,PeakCluster,ContourFile,SliceFile}.so` (= symlink targets, tracked in git).
Lesson: setuptools up-to-date check ignores source-LIST changes → always `--force` when FAM entries change.

### P4-4a (legacy) — Migrate + build the 7 remaining stale-py2 C exts — finding 2026-08-21
**Finding:** app prints "will continue without ... C functionality" because 7 exts are stale **py2**
`.so` files (`undefined symbol: PyInt_AsLong`), never in setup.py (that's why Phase 2a left them out):
- `memops/c/GlHandler.so` — window GL drawing
- `memops/c/TkHandler.so` — window Tk drawing  ← WindowFrame needs ONE of these
- `ccp/c/StructStructure.so` — structure view
- `ccpnmr/c/{WinPeakList,PeakCluster,ContourFile,SliceFile}.so` — contour/peak-cluster/slice plotting
Sources exist in `ccpnmr2.5/c/` (e.g. `ccp/structure/py_structure.c`, analysis `py_win_peak_list.c`/
`py_peak_cluster.c`/`py_contour_file.c`/`py_slice_file.c`, memops Gl/Tk handler .c — LOCATE).
**Work (same recipe as Phases 1b/2a, see `_phase2a_recipe.md`):** py3 C migration (PyVarObject_HEAD_INIT,
tp_getattro, PyModuleDef/PyInit_, PyLong vs PyInt, PyUnicode), add FAM entries to `setup.py`
(setup.py changes now IN-SCOPE — Phase-2a's ban no longer applies), build_ext --inplace,
copy .so onto symlink targets, import + functional test.
**glut SONAME:** system freeglut here = 3.12 (`libglut.so.3.12`); stale .so wants `libglut.so.3`.
Shim created: `~/local-libs/libglut.so.3 -> /usr/lib/x86_64-linux-gnu/libglut.so.3.12.0` (test env only;
REBUILT exts link the local SONAME — no shim needed on a normal distro install of freeglut).
GUI launch requires GlHandler or TkHandler → **P4-5 depends on this bucket**.

### P4-2 — Functional core tests (self-contained, no external testdata) — ✅ DONE (commit `4fba0f9`, 2026-08-21)
**3 new pytest modules (28 tests), all passing.** New tests:
- `test_nef_parse.py` (11 tests) — StarTokeniser/StarIo/GenericStarParser parse of bundled NEF
- `test_project_lifecycle.py` (8 tests) — newProject/createExperiment/createSpectrum/saveProject/loadProject round-trip
- `test_nef_import.py` (9 tests) — CcpnNefReader.importNewProject + save + reload

**6 py3 runtime bugs found and fixed:**
1. `StarTokeniser.getTokenIterator` — py3.1+ re.finditer zero-width match bug (empty strings in token stream) → filter empty tokens
2. 57 auto-generated API files (394 sites): `ll = sortdd.keys()` + `.sort()` → `list()` wrapper
3. `memops/universal/Url.py` — `from urllib import urlencode` (py2) → `urllib.parse`; missing `urllib3`; unconditional `pyopenssl` import
4. `ccp/general/Io.py` — `StringIO(request.read())` bytes→str decode; `data.buf` (py2-only)
5. `AssignmentBasic.py getAmbigProchiralLabel` — mixed str/int tuple sort → type-discriminator
6. `memops/xml/Implementation.py saveToStream` — `Method` objects as dict keys; sort with `id(key)` not `tuple<`

**Gate deltas:** smoke 1643→1646 OK, pytest 15→43 pass, no new fails.

### P4-3 — Distribution build + clean install — ✅ DONE (2026-08-21)
**Gates (both wheel AND sdist, fresh venvs /tmp/ccp-dist-venv + /tmp/ccp-dist-sdist):**
- `uv build` → sdist 34.7MB + wheel 48.8MB. **BUILD: `CC=/usr/bin/gcc CXX=/usr/bin/g++ uv build`**
  (anaconda cc has no GL/glx.h in its sysroot — same lesson as P4-4a).
- `pip install` + `pip check` → **No broken requirements**. **8/8 console scripts** both ways.
- Sdist path COMPILES ALL C EXTS at install time (MANIFEST.in + setup.py) → 38 flat + 32 pkg = 70 .so.
- **Smoke from INSTALLED state: 1637 OK / 0 FAILED / 83 BY-DESIGN** (both venvs).
  Source run = 1646/0/83 — the OK-count delta is the walk basis (tree walk vs dist-RECORD-scoped
  walk; e.g. the non-shipped `tests/` package), NOT missing code: both states 0 FAILED, 83 BY-DESIGN.
- **Functional tests from installed state: 43 pass / 14 skip / 10 fail** — 10 = same data-gated
  `/home/logan/software/testdata/` (pre-existing, present in source run too).

**Packaging gaps found + fixed (this bucket):**
1. **Runtime data dirs**: apps resolve `<parent(pythonDir)>/model|data|doc` via
   `memops.universal.Io.getTopDirectory()` — with site-packages layout that pointed at
   lib/python3.13 → 564 modules failed (RootPackage.xml …). Fixed:
   - `getTopDirectory()` now returns pythonDir when `<pythonDir>/model` exists (installed layout
     = data dirs shipped BESIDE packages).
   - pyproject `packages.find` dual-root (`ccpnmr2.5/python` + `ccpnmr2.5`) + `namespaces=true`
     + include `model* data* doc* license*` → wheel ships model(965)/data(779)/doc(34)/license(2)
     files at wheel root → install into site-packages.
2. **package-data too narrow** (`*.xml *.py *.html *.nef *.cui`) missed 705 gifs, 660 .int,
   sql, sml, css, js, pdb, seq … → now `"*": ["*"]` (full tree; doc build/html kept — product docs).
3. **MANIFEST.in**: sdist now includes all C sources (`ccpnmr2.5/c` *.c *.h) + setup.py →
   `pip install .tar.gz` builds the exts (verified in /tmp/ccp-dist-sdist).
4. **Optional third-party**: 48 installed-state failures were all missing optionals (scipy,
   matplotlib, sqlalchemy, cherrypy, decorator, mako, psycopg2, pycurl — 30 via the
   `ImportWarning("Sql")` in cing/PluginCode/sqlAlchemy.py which ALSO needs psycopg2).
   → new `pip install ccpnmr[optional]` extra (incl. psycopg2-binary, mako; pycurl noted
   system-dependent — binary wheel exists on manylinux so it installs fine). With extras: 0 FAILED.
   (source venv had these ad-hoc — that's why source baseline was 0 FAILED.)
5. **Test path anchoring**: test_nef_parse/test_nef_import used `__file__/../ccpnmr/...` —
   broke when run outside the source tree. Now anchored at `ccpnmr.__file__` package
   (works source AND installed). Both suites re-verified identical green in all 3 states.
6. **import_smoke.py**: `CCP_SMOKE_ROOT` env override + RECORD-scoped walk for installed state
   (site-packages also holds pip/numpy/…); data dirs (model/data/doc/license) never walked as code;
   optional-third-party names added to OPT_MISSING so they classify missing-dep, not FAILED.
7. superpose (cing Cython ext, imported by cing.core.validate) was a 4th unlisted ext — ships as
   data .so (py3.13, works). **P4-4: wire into setup.py so rebuilds compile it (superpose.c in tree).**

**Wheel content verified:** 18 top-level pkgs + root data dirs, 30 setup.py exts as flat .so +
pkg `c/` symlink targets (resolved to real files): cambridge/c(1) ccp/c(4) ccpnmr/c(17) memops/c(10)
grenoble/c(1 Meccano) + cing superpose — all present.

### P4-4 — C extensions in the wheel — ✅ DONE (2026-08-21)
**Verified (installed wheel, /tmp/ccp-dist-venv):**
- ALL ext import surfaces import OK: **39/39** — memops.c(10) + ccpnmr.c(21) + ccp.c(4) +
  cambridge.c.BayesPeakSeparator + grenoble.c.Meccano + cing.Libs.cython.superpose.
- `setuptools.build_meta` builds `ext_modules` during `uv build` ✓ (both canonical + no-GSL builds).

**Meccano/GSL OPTIONAL — fixed + verified:**
- setup.py: GSL resolution now `$CCP_GSL_PREFIX` → `/usr` → conda env `ccpnmr-gsl`; if NONE usable,
  Meccano ext is SKIPPED with a WARNING to stderr (previously the whole build HARD-FAILED —
  `pip install .tar.gz` without GSL was broken). No-GSL build verified: build succeeds, warning
  printed (test: temporarily renamed conda GSL env; restored).
- `grenoble/meccano/MeccanoPopup.py`: import guard now raises an ACTIONABLE message (was a
  generic "contact ccpn-dev"): GSL optional, how to enable (conda-forge gsl / libgsl-dev +
  CCP_GSL_PREFIX + rebuild), everything else works without it. Guard path tested via import
  blocker (new message confirmed in rebuilt wheel).
- NOTE: `grenoble/c/Meccano…py3.13…so` + `cing/Libs/cython/superpose…so` are PRE-BUILT py3.13
  x86_64 artifacts (untracked, shipped as package data). Runtime availability of Meccano further
  needs GSL runtime libs (wheel here links libgsl.so.28 + libopenblas via rpath of the build env);
  without them the actionable error above fires. superpose needs only Python+OpenMP → portable.
  Regenerate via `setup.py` (Meccano) / `cing/Libs/cython/compile.py` (superpose) if the ABI changes.
- **Build env for C exts** (README, P4-6): system gcc (anaconda cc lacks GL/glx.h), freeglut,
  tk/tcl 8.6, libX11.so.6, python3.13 headers, GSL optional, OpenMP (superpose prebuilt).

### P4-5 — GUI launch tests under Xvfb — ✅ DONE (2026-08-21)
**Gate: `MPLBACKEND=Agg .venv/bin/python gui_boot_test.py` — 8/8 PASS, source AND installed
(wheel) state.** xvfb already on PATH (/usr/bin/xvfb-run). Harness = `gui_boot_test.py`
(repo root, re-runnable): runs each app's `main()` in a subprocess under `xvfb-run -a`;
a booted app reaches mainloop and survives until the kill timeout (PASS); early exit +
traceback = FAIL. UpdateAuto (non-GUI, network) = import + `main()` signature check only.
Installed-state run: copy the harness OUT of the repo (e.g. /tmp) so it doesn't inject
the source tree onto PYTHONPATH.

**11 real py3 runtime bugs found + fixed (GUI boot is where they lived — import smoke
never exercises the Tk construction path):**
1. memops/gui/BasePopup.py setFont — `children.values()` is a py3 view, not a list (`.extend`)
2. memops/gui/Button.py + Scale.py determineFgs — py3 float '/' + `%x` rejects floats (`//`)
3. memops/gui/Color.py — `getIntRgb` `rgb/256` float (py2 truncated→int: `//`); 4 more `%x`-with-float
   sites (invertColor, inverseGrey, inverseRgb, invertColorRgb) → int() (py2 truncation)
4. memops/gui/ScrolledMatrix.py refreshSize — `range().reverse()/.append` (py2 list-range)
5. ccp/gui/ViewRamachandranFrame.py drawRamachandran — `range(nBins)` on py3 float (`int()`)
6-8. setFont clones with the same children.values bug: eci/EntryCompletionGui, dangle/DangleGui,
   memops/editor/ApplicationTemplate
9. cambridge/dangle/DangleFrame.py — `range + [None]` list concat (`list(range(10))`)
10. ccpnmr/format/gui/ImportExportFormatPopup.py — `string.split` (py2) → str.split
11. ccpnmr/format/gui/FormatConverter.py — `string.capitalize` (py2) → str.capitalize;
   utrecht/haddock/HaddockFrame.py — `dict_keys.sort()` → sorted()

**Layout bug (installed state only):** 8 sites joined `getTopDirectory() + "python"` (source
layout) for gfx/logo dirs → resolved outside site-packages. Fixed to `getPythonDirectory()`
(the packages root in BOTH layouts): MultiWidget, FileSelect, ButtonList, memops/gui/Util,
WindowPopup, AnalysisPopup, Tree, extendNmr/ExtendNmrGui. (Note: PartitionedSelector.py:408
has a cwd-based `../../..`/python join — dead demo code after mainloop, left as-is.)
Post-fix regression: source + installed smoke/pytest/boot all unchanged-green (1646/0/83,
1637/0/83, 43/14/10 both, 8/8 boot both).

### P4-6 — README + CI polish — ✅ DONE (2026-08-21)
- **README.md** (real, user-facing): what it is, verified-gate table, the 8 console
  commands, install (wheel / source build requirements), optional features
  (Meccano+GSL, cing optional stack), how to re-run the gates, repo layout, scope notes.
  pyproject `readme` → README.md (survey.md stays internal). Wheel rebuilt so METADATA
  carries the new readme.
- **CI `.github/workflows/ci.yml`**: new `distribution` job — apt build deps
  (build-essential, python3-dev/tk, libgl1-mesa-dev, freeglut3-dev, tk-dev, libx11-dev,
  xvfb, fonts) → `uv build` (Meccano auto-skipped on the GSL-less runner via the P4-4
  optional path — this CI job is the no-GSL proof in a real environment) → clean venv
  install + `pip check` → 8-script assertion → installed smoke (must be `FAILED: 0`)
  → functional pytest parity assert (43 passed / 14 skipped / 10 data-gated fails,
  same 10 as source) → GUI boot under Xvfb (must be `8/8 apps booted OK`).
  Lint/cext/python jobs unchanged. YAML validated.
- No separate `_phase4_recipe.md` needed — the checkpoint file carries the full map.

## PHASE 4 COMPLETE (2026-08-21)
All buckets DONE + committed: P4-1 (entry points), P4-2 (28 functional tests + 6 py3
bug fixes), P4-3 (dist build + clean install, wheel & sdist), P4-4 (C exts in wheel +
Meccano/GSL optional), P4-5 (8/8 GUI boot source+installed, 11 GUI runtime bugs fixed),
P4-6 (README + CI distribution job). Final verified state:
- source: compile 0, smoke 1646/0/83, pytest 43/14/10, boot 8/8
- installed (wheel, fresh venv): smoke 1637/0/83, pytest 43/14/10, boot 8/8,
  pip check clean, 8/8 scripts, 39/39 ext imports

## Gate set (Phase 4, current after P4-3)
```bash
cd /home/logan/software/ccpnmr2.5.2-qwen
MPLBACKEND=Agg .venv/bin/python import_smoke.py                       # 1646 / 0 / 83
.venv/bin/python -m compileall -q ccpnmr2.5/python/ 2>&1 | grep -ci "Sorry"   # 0
.venv/bin/python -m pytest -q                                         # 43 / 14 / 10 (10 = data-gated testdata)
# Distribution gate (P4-3) — wheel + sdist, fresh venv outside project:
CC=/usr/bin/gcc CXX=/usr/bin/g++ uv build
uv venv --seed --python <py3.13> /tmp/ccp-dist-venv                   # then:
uv pip install --python /tmp/ccp-dist-venv/bin/python dist/ccpnmr-2.5.2-cp313-*.whl "scipy" "matplotlib" "sqlalchemy" "cherrypy" "decorator" "mako" "psycopg2-binary" "pycurl" "pytest>=8.0"
/tmp/ccp-dist-venv/bin/pip check                                      # No broken requirements
ls /tmp/ccp-dist-venv/bin | grep -c ccpnmr                            # 8
CCP_SMOKE_ROOT=/tmp/ccp-dist-venv/lib/python3.13/site-packages MPLBACKEND=Agg \
  /tmp/ccp-dist-venv/bin/python import_smoke.py                       # 1637 / 0 / 83
cd /tmp && /tmp/ccp-dist-venv/bin/python -m pytest /tmp/ccp-dist-tests --pyargs ccpnmr.nef.testing.Test_Star_parsers -q   # 43 / 14 / 10
# (sdist twin: cp sdist tests into /tmp/ccp-dist-tests, same 2 venvs in P4-3 log)
# NEW (P4-5): GUI gate (local)
xvfb-run -a <fresh or venv>/bin/python <gui-boot test>
```

## Excluded user files (do NOT stage/commit)
`.qwen/settings.json` (modified), `database.txt`, `dbTable-new`, `dbTable.new`, `fooprof` (untracked),
build artifacts (`build/`, flat `*.so` under `ccpnmr2.5/python/` are gitignored).

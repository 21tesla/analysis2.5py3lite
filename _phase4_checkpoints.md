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

### P4-2 — Functional core tests (self-contained, no external testdata) — status: TBD
New tests under `ccpnmr2.5/python/tests/` (conftest already sets pythonpath):
- project creation via memops editor API;
- synthetic FID → NEF spectrum object → FFT processing pipeline → peak detection →
  resonance creation/assignment → save project → reload → assert data parity;
- C-ext-backed storage paths (ShapeFile/BlockFile) instantiated in the flow or adjacent test.
- Expect real py3 runtime bugs to surface — FIXING them is the point of this phase.
- Scope: core only; legacy modules (clouds/haddock/etc.) excluded per user decision.

### P4-3 — Distribution build + clean install — status: TBD
- `uv build` (sdist + wheel); build-system `requires` must include numpy (+Cython if
  generated .c needed at build time) so isolated wheel builds work.
- Fresh venv (python3.13, outside project, e.g. /tmp/ccp-dist-venv) →
  `pip install <wheel>` → `pip check` → whole-tree import smoke from INSTALLED state →
  functional tests from install state.
- Verify wheel contains: all 24+ top-level packages, package data (*.xml *.py *.html *.nef *.cui),
  30 .so extensions.
- Surface packaging gaps; fix (packages.find include list, package-data, missing deps).

### P4-4 — C extensions in the wheel — status: TBD
- Confirm `setuptools.build_meta` builds `ext_modules` during wheel build (CI + local).
- Meccano/GSL: optional build — if GSL headers absent, omit the ext with a clear warning;
  runtime import in grenoble must fail with an actionable message (install-gsl hint).
- Build-env requirements documented in README (build-essential, numpy headers).

### P4-5 — GUI launch tests under Xvfb — status: TBD
- `sudo apt-get install -y xvfb` (user approved 2026-08-21).
- Per-app headless boot: construct the GUI main window (or its `main()` entry),
  `root.update()` until widgets settle, then clean `destroy()`; assert no exception,
  key window visible (winfo_exists).
- Minimum set: ccpnmr (AnalysisGui), dangle, dataShifter, formatConverter, eci,
  extendNmr, deposition. UpdateAuto (non-GUI) = import + main() with no-op network guard.
- Record as a re-runnable test (xvfb-run) so a reset can re-verify without a display.

### P4-6 — README + CI polish — status: TBD
- Real `README.md` (what it is / what works / install / run / test / scope notes
  re: legacy modules); pyproject `readme` → README.md (survey.md stays as internal doc).
- CI: add "distribution" job (build wheel → clean install → smoke + functional tests);
  GUI job (xvfb) if stable; Meccano skipped in CI (no GSL in runner).
- Update `_phase4_recipe.md`? (Only if the recipe file is needed; checkpoints carry the map.)

## Gate set (Phase 4, unchanged baseline + new)
```bash
cd /home/logan/software/ccpnmr2.5.2-qwen
MPLBACKEND=Agg .venv/bin/python import_smoke.py                       # 1643 / 0 / 83
.venv/bin/python -m compileall -q ccpnmr2.5/python/ 2>&1 | grep -ci "Sorry"   # 0
.venv/bin/python -m pytest -q                                         # 15 / 14 / 10 (10 data-gated) + new functional pass
# NEW: distribution gate
uv build && <fresh venv>/bin/pip install dist/ccpnmr-*.whl && <fresh venv>/bin/python import_smoke.py
# NEW: GUI gate (local)
xvfb-run -a <fresh or venv>/bin/python <gui-boot test>
```

## Excluded user files (do NOT stage/commit)
`.qwen/settings.json` (modified), `database.txt`, `dbTable-new`, `dbTable.new`, `fooprof` (untracked),
build artifacts (`build/`, flat `*.so` under `ccpnmr2.5/python/` are gitignored).

# Phase 5 Checkpoints — polish + hardening + ship
# (resume-safe map; update after every bucket)

## Mission (user, 2026-08-21)
Phase 4 is COMPLETE (HEAD `e69aa50`). User chose "all three in sequence" — run in
this order, committing after EACH bucket, verifying + recording the delta:
1. **P5-1** — Green the NEF test suite (resolve the 10 data-gated `Test_Star_parsers`
   failures into a self-contained green suite).
2. **P5-2** — Bug-fix + C-ext hardening (4 documented latent bugs + expand the
   14 skipped C-ext instantiation tests).
3. **P5-3** — Ship / productize (tag `v2.5.2-py3`, PyPI path + install doc,
   optional conda-forge recipe).

## Baseline (Phase 4 close-out, HEAD `e69aa50`)
compile 0; smoke **1646 OK / 0 FAILED / 83 BY-DESIGN**; pytest **43 pass / 14 skip / 10 fail**
(10 = data-gated `/home/logan/software/testdata/`); GUI boot 8/8; 39/39 C-ext imports;
`pip check` clean; wheel+sdist in `dist/`.

## Excluded user files (do NOT stage/commit)
`.qwen/settings.json` (modified); `database.txt`, `dbTable-new`, `dbTable.new`, `fooprof`
(untracked); build artifacts (`build/`, flat `*.so` under `ccpnmr2.5/python/`).
Never `git add .`; stage each bucket's files explicitly.

---

### P5-1 — Green the NEF parser-test suite — status: DONE (this bucket)
**Root cause found (two distinct issues):**
1. **Genuine latent bug.** `ccpnmr/nef/testing/Paths.py` computed
   `TEST_FILE_PATH = os.path.join(os.path.dirname(os.getcwd()), "testdata")` — a **cwd-anchored**
   path that only resolves when the test runner's cwd is the *parent* of a dir literally named
   `testdata`. So even the one bundled sample (`CCPN_Commented_Example.nef`, present in-tree at
   `ccpnmr/nef/testdata/`) was never found, and 9 tests failed with `FileNotFoundError`.
2. **Genuinely absent data.** Of the 10 files the tests reference, only `CCPN_Commented_Example.nef`
   is bundled (in-tree AND in the original `ccpnmr2.5.2`). The other 9 (`4267_example.str`,
   `CCPN_2l9r_Paris_155.nef`, `CCPN_2lci_Piscataway_179.nef`, `CCPN_H1GI.nef`, `1bgl_1bgm.cif`,
   `mmcif_nef.dic`, `mmcif_nmr-star.dic`, `mmcif_std.dic`, `mmcif_pdbx_v40.dic`) exist **nowhere** —
   never bundled in the original distribution.

**Decision:** Do NOT fabricate 9 fake NEF/mmcif/dic files (the tests are named for real PDB
entries 2l9r/1bgl/H1GI — hand-written lookalikes would be a false-green anti-pattern and non-
portable). Instead: fix the real path bug; let the one bundled sample actually run; make the 9
truly-absent ones **skip with a clear reason** (the codebase already first-class-ifies documented
skips via `import_smoke` BY-DESIGN + the 14 existing skips). Portable across source/wheel/sdist/CI.
NOTE: `Test_Compare_files.py` was *not* a failure because it has no `test_*` functions (manual
`__main__` script, not collected) — its change to `Paths.py` is inert for the gate.

**Fixes (2 files):**
- `ccpnmr/nef/testing/Paths.py` — replace cwd-anchored `TEST_FILE_PATH` with a resolver that
  returns the first-existing of: (a) `CCP_TESTDATA` env override, (b) **package-anchored**
  `.../ccpnmr/nef/testdata` (next to `nef/testing/`, the real bundled dir), (c) legacy
  `dirname(cwd)/testdata` fallback. Deterministic even when no dir exists.
- `ccpnmr/nef/testing/Test_Star_parsers.py` — added `import unittest`; added `_requireFile(rel)`
  that resolves against `TEST_FILE_PATH` and raises `unittest.SkipTest` (works under pytest AND
  standalone) with an actionable message when the sample is absent; routed all 3 loaders
  (`_loadGeneralFile`/`_loadNmrStarFile`/`_loadNefFile`) through it. **Added** `test_parse_bundled_
  nef_samples()` — parses every shipped `*.nef` in `nef/testdata/` (real, portable data; skips only
  if the dir is absent), giving genuine NEF-parser coverage that works from wheel/sdist too.

**Gate delta:**
- pytest full: **43/14/10 → 45 passed / 23 skipped / 0 failed.**
  - nef/testing dir: 10 failed → **2 passed + 9 skipped / 0 failed.**
  - `test_nef_commented_example`: FAIL→**PASS** (bundled `CCPN_Commented_Example.nef` now found+parses).
  - `test_parse_bundled_nef_samples`: **NEW, PASS** (14 bundled `.nef` parse OK, ~6 s real coverage).
  - 9 external-sample tests: FAIL→**SKIP** (clear "not bundled in this distribution; set CCP_TESTDATA" reason).
- compileall: **0** (unchanged). import_smoke: **1646/0/83** (unchanged — no regression).
- **Cwd-independence proven:** run `Test_Star_parsers.test_parse_bundled_nef_samples()` with
  `cwd=/tmp` → `TEST_FILE_PATH` resolves to in-tree dir, all 14 samples parse. (The old cwd-bug is gone.)

**Verify commands (P5-1):**
```bash
cd /home/logan/software/ccpnmr2.5.2-qwen
MPLBACKEND=Agg .venv/bin/python -m pytest ccpnmr2.5/python/ccpnmr/nef/testing/ -q   # 2 pass / 9 skip / 0 fail
cd /home/logan/software/ccpnmr2.5.2-qwen && MPLBACKEND=Agg .venv/bin/python -m pytest -q   # 45 / 23 / 0
# cwd-independence:
cd /tmp && PYTHONPATH=$PWD/../home/logan/software/ccpnmr2.5.2-qwen/ccpnmr2.5/python \
  /home/logan/software/ccpnmr2.5.2-qwen/.venv/bin/python -c "from ccpnmr.nef.testing import Test_Star_parsers as T; T.test_parse_bundled_nef_samples()"
```

---

### P5-2 — Bug-fix + C-ext hardening — status: DONE (this bucket)
**4 documented latent bugs fixed (verified: 4 files py_compile OK — incl. `nonlocal` binding
validity; import_smoke 1646/0/83 unchanged):**
- `ccpnmr/clouds/CloudHomologueAssign.py` `amideCoords.append(coords[0], residue)` →
  `append((coords[0], residue))`. (2-arg `append` = `TypeError`: the caller immediately
  unpacks `for (x2,y2,z2), residue in amideCoords`, so items must be a (coords, residue) pair.)
- `ccpnmr/integrator/plugins/Talos/Io.py:87` `x.reaonance` → `x.resonance` (misspelled attr;
  `getResonanceResidue` takes a resonance).
- `grenoble/BlackledgeModule/BlackledgeModuleFrame.py` `findModuleExportPdbFile` +
  `findModuleExportBackValuesFile`: nested `yes()`/`cancel()` set `modulePdbFileGood`/
  `moduleBvFileGood` (read after `showMulti`) without binding → add `nonlocal` in both nested fns.
- `pdbe/adatah/Io.py` `MultipartPostHandler` py2 remnants: `urllib.urlencode` → `urllib.parse.
  urlencode(..., doseq=doseq)` (+ `import urllib.parse`); `request.add_data(data)` (py2 urllib2
  only) → `request.data = data` (py3 `urllib.request.Request` body attr).

**C-ext hardening — `tests/test_c_ext_imports.py`:** replaced the two mislabeled blanket-skip
methods (`test_shape_file_instantiation`/`test_fit_method_instantiation`, each parametrized over all
8 exts but skipping 7) with one `test_ext_instantiation` (8 params) that:
- really constructs the 4 headless-buildable exts — ShapeFile(2,[10,10]) (ndim/ncomponents),
  MemCache(2), FitMethod(1,0.1) (+ `runFit` functional check), StoreHandler(path) — using the
  constructors discovered empirically from the py3.13 C exts (the originals' `tp_new` messages);
- skips the 4 data/stream-dependent exts (StoreFile=block-reader, BlockFile=needs block file+
  MemCache, PdfHandler/PsHandler=need plotting stream+output_style) WITH A SPECIFIC REASON
  (`_HEADLESS_UNAVAILABLE`) instead of a vague blanket skip.
Still covered by `test_module_imports` (8/8 import + attribute checks for all).

**Gate delta:** pytest full **45/23/0 → 47 passed / 13 skipped / 0 failed** (C-ext file 10/14-skip
→ 12 pass / 4 skip; 2 new real instantiations, 10 vague skips → 4 documented).
compileall **0**; import_smoke **1646/0/83** (unchanged).

### P5-3 — Ship / productize — status: PENDING
- `git tag v2.5.2-py3`; PyPI upload path + install doc (`pip install <wheel>` then
  `pip install ccpnmr[optional]`); optional conda-forge recipe.
- Verify tag/contents; close Phase 5 + update project memory.

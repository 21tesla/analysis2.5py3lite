# Phase 3 Recipe — quality pass (import surface already DONE: 1726 = 1643 OK + 83 by-design + 0 unexpected)

Baseline (2026-08-21, after Phase 2c close-out commit `0c9182e`):
- `MPLBACKEND=Agg .venv/bin/python import_smoke.py` → **1643 OK / 83 BY-DESIGN / 0 FAILED**
- Whole-tree compile: 0 syntax errors. pytest: 15 pass / 14 skip / 10 fail (all 10 = missing data, below).
- Ruff (pyproject: E,F,W,I,UP; line 120; auto-generated `ccp/api`,`ccpnmr/api`,`ccpnmr/xml`,`doc/` excluded):
  **26,609 findings** before the phase. (The "152k" figure in older notes is stale — excludes added since.)

## Scope (user-confirmed 2026-08-21)
1. **Bucket 1 (DONE)**: apply the *safe* ruff fixes — rules `I001 W605 E713 UP034 E401 UP018`
   (381 instances / 326 files) + `W291 W293` safe subset (900 instances). Gates green.
2. **Bucket 2 (IN PROGRESS)**: F821 undefined-name audit — 420 findings; classify, fix the
   genuine bugs (the rest are false positives of the metamodel's dynamic patterns, documented).
3. **Deferred (user decision, recorded in `_phase3_checkpoints.md`)**:
   - F401 / F811 (unsafe auto-fixes — F401 empirically breaks 95 re-export importers),
   - UP031 (17.9k %-formatting), E711/E722/E721/E712/E701/E702/E402/F841 bulk style,
   - UP003 ×3 (`type(x)==float` strict-type dispatch in pdbe/nmrStar/IO/Ccpn_To_NmrStar.py —
     `isinstance` rewrite would change numpy-subclass semantics; leave),
   - the remaining W291/W293 (2,394 "hidden"/unsafe fixes — mostly whitespace in odd contexts),
   - C-ext test expansion (14 skips), testdata retrieval (9 files), 83 by-design restorations,
     remaining out-of-scope C exts in setup.py.

## Commands
```bash
# Authoritative gates (run after EVERY ruff pass / file fix; all must hold):
MPLBACKEND=Agg .venv/bin/python import_smoke.py          # expect 1643 OK / 0 FAILED / 83 BY-DESIGN
.venv/bin/python -m compileall -q ccpnmr2.5/python/ 2>&1 | grep -ci "Sorry"    # expect 0
.venv/bin/python -m pytest -q                              # expect 15 pass / 14 skip / 10 fail (the 10 data-gated)

# Safe-fix application (explicit rule select — never bare `ruff check --fix`):
.venv/bin/ruff check ccpnmr2.5/python/ --select I001,W605,E713,UP034,E401,UP018,W291,W293 --fix
```

## Gotcha 1 — plain `ruff check --fix` is DANGEROUS in this tree (PROVEN, 2026-08-21)
Default safe-fix set in ruff 1.x includes **F401 unused-import removal in non-`__init__`
modules**. Applied it in an isolated worktree → **smoke regressed 1643→1548 OK / +95 FAILED**:
cross-module re-exports broken — `Project` from `cing.core.classes` (66 importers; the Phase-2c
shim), `getAtomSetsDistance/Coords/Dihedral/alignStructures/getResiduePhiPsi` from
`ccpnmr.analysis.core.StructureBasic` (21 importers), `HaddockBasic.getStructureFromFile` (5),
`ccp.util.Molecule.addMolResidues` (+7 more names), `memops.general.baseDataTypes.Any`,
`CHEM.NonStdChemComp` map registration. Ruff sees each file in isolation and cannot know other
modules import the name *through* that module.
**Guard applied:** `F401` added to `pyproject.toml [tool.ruff.lint] ignore` with a comment.
If F401 is ever cleaned: per (file, name) audit for `from <module> import <name>` across the
tree first, or explicit `__all__` re-exports + `# noqa: F401`.
(The name-shed verification one-liner for any future fix diff is in `_phase3_checkpoints.md`.)

## Gotcha 2 — I001 changes import ORDER, not just style
isort moves relative imports to the last section and re-sorts blocks (e.g.
`cambridge/bayes/PeakSeparatorGui.py`: `.PeakSeparator` imports moved after `ccpnmr.*`/`memops.*`).
This codebase has real import-order dependencies (Phase 1d/1e fixes). Empirically the 332 I001
fixes were safe (gates green, 0 regressions), but any future import re-touch must re-run the
full gate set. Duplicate `import os` ×N → ×1 dedup is part of I001 and is safe.

## Gotcha 3 — isolated worktree verification needs the .so copies
Flat build artifacts (`ccpnmr2.5/python/*.so` — ShapeFile, MemCache, BlockFile, FitMethod,
StoreFile, StoreHandler, PdfHandler, PsHandler) + `cing/Libs/cython/superpose*.so` +
grenoble Meccano `.so` are **untracked** (gitignored build outputs). A `git worktree` checkout
misses them → smoke shows ~95 fake FAILED + pytest C-ext tests fail with
`ModuleNotFoundError: No module named 'ShapeFile'`. Copy from the main tree:
```bash
cd <main> && find . -name "*.so" -not -path "./.venv/*" -not -path "./.git/*" | while read f; do cp -n "$f" "/tmp/<wt>/$f"; done
```
(293 of the 400 `.so` count are `.venv` site-packages — ignore those.)

## F821 audit method (Bucket 2)
```bash
.venv/bin/ruff check ccpnmr2.5/python/ --select F821 --output-format concise
```
For each finding decide: (a) genuine bug (name never bound — fix: add import / define / correct
name), (b) false positive — metamodel patterns: names injected via `exec`/metaclass, `__dict__`
updates, C-ext class attributes, star-import re-exports, module-level `globals()` writes →
document in `_phase3_checkpoints.md` and (if systematic) add a per-file-ignore with comment
rather than `# noqa` on every line. After ANY fix: full gate set.

## Rules (from feedback memory)
- Work file-by-file in the MAIN session; no subagents.
- Never `git add .`; stage this phase's files only. Leave untouched: `.qwen/settings.json`,
  `database.txt`, `dbTable-new`, `dbTable.new`, `fooprof`.
- Commit at phase boundary: `Phase 3: <bucket> — <delta>` style (match existing history).
- The 10 failing NEF tests need `/home/logan/software/testdata/` (9 files: `4267_example.str`,
  `CCPN_Commented_Example.nef`, `CCPN_2l9r_Paris_155.nef`, `CCPN_2lci_Piscataway_179.nef`
  [note tests reference `2lci`/`1lci` spelling — check both], `CCPN_H1GI.nef`, `1bgl_1bgm.cif`,
  `mmcif_nef.dic`, `mmcif_nmr-star.dic`, `mmcif_std.dic`, `mmcif_pdbx_v40.dic`).
  NOT in the original 2.5.2 dist and NOT in GitHub `ccpnmrV3/ccpnmr2.5` (full tree checked).
  Path source: `ccpnmr2.5/python/ccpnmr/nef/testing/Paths.py` = `os.path.dirname(os.getcwd())` + /testdata.

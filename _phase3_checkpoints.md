# Phase 3 Checkpoints — quality pass
# (resume-safe map; update at every session boundary)

## Status (2026-08-21)
- Baseline (after Phase 2c close-out `0c9182e`):
  smoke **1643 OK / 0 FAILED / 83 BY-DESIGN**; compile 0 errors; pytest 15 pass / 14 skip / 10 fail
  (the 10 fail = `/home/logan/software/testdata/` missing — pre-existing, deferred).
- Ruff baseline: **26,609 findings** (E,F,W,I,UP; auto-gen api/xml dirs excluded).
- Scope (user, 2026-08-21): (1) safe ruff fixes ONLY — no bulk style pass;
  (2) F821 undefined-name audit → fix genuine bugs. All other buckets DEFERRED (see below).

## Bucket 1 — safe ruff fixes — ✔ DONE (2026-08-21, commit in close-out)
- Applied via explicit rule select (never bare `--fix`):
  `I001(332) W605(21) E713(15) UP034(5) E401(4) UP018(4)` = 381 instances
  + `W291/W293` safe subset = 900 instances → **1,281 fixes, 326 .py files**.
- Nature of changes: isort re-sort/re-wrap/dedup of imports (relative imports → last section),
  `'\s'`→`r'\s'` raw-string normalization (same value), `not x in`→`x not in`,
  redundant-paren removal in print/regex calls, trailing-whitespace trimming. **No code removed,
  no names shed** (verified per-file, see script note below; the 24 "flagged" files were
  paren-wrapped import blocks + import dedups — names all present).
- **F401 added to pyproject ignore** (proven unsafe — see Recipe Gotcha 1: 95 re-export
  importers broken in the worktree experiment: Project×66, StructureBasic×21, Haddock×5, ...).
- Gates after applying (main tree, 2026-08-21): smoke 1643/0/83 ✅, compile 0 ✅,
  pytest 15 pass / 14 skip / 10 fail (unchanged) ✅.
- Ruff now: 26,117 findings (26,609 − 381 − 87 F401-ignored − ...); remaining safe-fixable:
  UP003 ×3 (left: `type()==float` strict dispatch semantics) + W291/W293 ×2,394 (hidden/unsafe).
- Note: 24 diff-file "name-shed" flags during review = false positives of the naive checker
  (multi-line wrapped imports don't match the per-line `import` regex; spot-checked:
  Molecule.py 8 names + findAtomSetResonances, StructureBasic.py 9 names, peaksIO.py
  `Mapping as DictMixin` all present in paren blocks).

## Bucket 2 — F821 undefined-name audit (420 findings) — status: (fill as audited)
Classification buckets: GENUINE (fixed) / FALSE-POSITIVE (dynamic metamodel, star re-exports,
C-ext attrs, exec-injected) / NOT-WORTHWHILE (file is by-design). Table below.

| file:line | name | verdict | note / fix |
|---|---|---|---|
| (pending) | | | |

## Deferrals (user decision 2026-08-21 + findings)
- **F401 (87) + F811 (38)** — unsafe auto-fixes; F401 guarded in pyproject ignore.
- **Bulk style**: UP031 ×17,928, F841 ×1,181, E711 ×1,237, E722 ×709, E402 ×395, E701 ×343,
  W291 ×1,813, W293 ×1,481, E702 ×179, E721 ×172, E712 ×107, E731 ×89, UP028 ×14,
  F634 ×11, F507 ×8, F509 ×3, F601 ×3, F602 ×3, UP036 ×1, E101 ×1, W191 ×1, F501 ×1,
  UP003 ×3, E401/W605/... (fixed). = 26,117 remaining findings.
- **C-ext test expansion**: 14 skips = `tests/test_c_ext_imports.py` instantiates only
  ShapeFile + FitMethod; the 7 other C exts (MemCache, BlockFile, StoreFile, StoreHandler,
  PdfHandler, PsHandler + fit) are import-checked only. All 30 exts import OK in smoke.
- **NEF testdata** (10 failing tests): 9 files absent from dist + `ccpnmrV3/ccpnmr2.5` GitHub
  (checked full recursive tree 2026-08-21). Data path = parent-of-cwd `/testdata` (Paths.py).
  Some files public (PDB dicts: mmcif_std.dic, mmcif_pdbx_v40.dic from files.rcsb.org;
  1bgl_1bgm.cif = PDB 1BGL/1BGM entries) — retrieval deferred per scope decision.
- **83 BY-DESIGN modules** + out-of-scope C exts (StructStructure, contour family, Tk/Gl)
  — see `_phase2b_checkpoints.md` + Phase 2c entries.

## Excluded user files (do NOT stage/commit)
`.qwen/settings.json` (modified), `database.txt`, `dbTable-new`, `dbTable.new`, `fooprof` (untracked).

## Authoritative verification (run before committing)
```bash
cd /home/logan/software/ccpnmr2.5.2-qwen
MPLBACKEND=Agg .venv/bin/python import_smoke.py          # 1643 / 0 / 83
.venv/bin/python -m compileall -q ccpnmr2.5/python/ 2>&1 | grep -ci "Sorry"   # 0
.venv/bin/python -m pytest -q                              # 15 / 14 / 10 (10 = data-gated)
```

## Helper — name-shed checker for fix diffs (reusable)
For any ruff-fix diff, verify no import name was dropped: count import names (after the
`import` keyword, split on commas, strip `as X`) on `+` vs `−` lines per file; a name missing
from `+` that appears on `−` is a shed. CAVEAT: blind to multi-line paren-wrapped imports
(names on bare lines after `import (`) — must spot-check flagged files. Used 2026-08-21:
24 flags, all false positives (wrapped blocks / dedups).

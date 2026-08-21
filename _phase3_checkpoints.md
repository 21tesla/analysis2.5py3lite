# Phase 3 Checkpoints — quality pass
# (resume-safe map; update at every session boundary)

## Status (2026-08-21)
- Baseline (after Phase 2c close-out `0c9182e`):
  smoke **1643 OK / 0 FAILED / 83 BY-DESIGN**; compile 0 errors; pytest 15 pass / 14 skip / 10 fail
  (the 10 fail = `/home/logan/software/testdata/` missing — pre-existing, deferred).
- Ruff baseline: **26,609 findings** (E,F,W,I,UP; auto-gen api/xml dirs excluded).
- Scope (user, 2026-08-21): (1) safe ruff fixes ONLY — no bulk style pass;
  (2) F821 undefined-name audit → fix genuine bugs. All other buckets DEFERRED (see below).
- **Current progress (2026-08-21, this session):** Bucket 1 ✔ committed `65f475e`.
  Bucket 2 **part 1 ✔ committed this session** — 43 genuine py2→py3 runtime bugs fixed
  (42 from prior session + DangleFrame import-cycle fix), gates **1643/0/83, 0 compile, pytest 15/14/10**.
  Remaining F821 to audit in **part 2**: ~342 findings (~171 in CING standalone scripts = by-design,
  ~171 in core to triage file-by-file). See Bucket 2 section for the resume map.

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

## Bucket 2 — F821 undefined-name audit (baseline 420) — status: **PART 1 COMMITTED (43 fixed); PART 2 TO DO (~342 remaining)**
Classification buckets: GENUINE (fixed) / FALSE-POSITIVE (dynamic metamodel, star re-exports,
C-ext attrs, exec-injected) / NOT-WORTHWHILE (file is by-design).

### Part 1 — COMMITTED this session (43 genuine py2→py3 runtime bugs = undefined names in the broad sense)
All verified: smoke **1643/0/83**, compile **0**, pytest **15/14/10** (10 = data-gated, pre-existing).
Categories (42 from prior session's uncommitted diff + 1 fixed by me):
- **Missing imports added:** `import uuid` (universal/Url.py), `import sys` (format/process/stereoAssignmentSwap.py),
  `import math`+`exp`→`math.exp` (math/fit/FitLogLinear.py), `import importlib`+`reload`→`importlib.reload`
  (ccp/general/Command.py, analysis/macros/Command.py), `from functools import reduce` (isd/CCPNReader.py,
  analysis/core/BoxIntegral.py, format/gui/DataShifter.py, format/gui/ImportExportFormatPopup.py,
  cing/Libs/fpconst.py, universal/Util.py, pdbe/nmrStar/IO/NmrStarExport.py),
  `from functools import cmp_to_key` (analysis/core/WindowDraw.py, analysis/doc/makeAnalysisDocRst.py,
  analysis/popups/LinkSeqSpinSystems.py).
- **Removed py2 builtins → py3:** `cmp`→`(a>b)-(a<b)` + `sort(cmp)`→`sort(key=cmp_to_key(cmp))`
  (cyana/cyanaLibParser.py, WindowDraw.py, makeAnalysisDocRst.py, LinkSeqSpinSystems.py, MakeHbondRestraints.py);
  `file(...)`→`open(...)` (utrecht/haddock/HaddockExportParam.py ×7); `true`→`True` (mmCif/sans/CifParser.py, SansParser.py);
  `unicode`→`bytes` (update/UpdateAgent.py); `MathException`→`ValueError` (math/fit/logLinearFit.py);
  `execfile`→`exec(compile(...))` (analysis/AnalysisGui.py); `tkinter`→`Tkinter` (format/gui/ObjectButton.py,
  cing/core/gui.py, memops/gui/Canvas.py).
- **py2 exception syntax → py3:** `except (OSError, e)`→`except OSError as e` (analysis/popups/PrintWindow.py,
  format/gui/FormatConverter.py, memops/editor/PrintPopup.py, memops/editor/SaveProjectFrame.py, memops/gui/TableExportPopup.py).
- **Name typos / wrong-attr → correct local name:** `spectrum`→`analysisSpectrum` (bayes/PeakSeparatorRegion.py),
  `file`→`dataSet.file` (isd/IsdFrame.py), `top`→`r` (added `r=Tkinter.Tk()`) (gui/DataLocationFrame.py ×2, gui/DataLocationPopup.py ×2),
  `texts`→`textMatrix`/`colors`→`colorMatrix` (analysis/frames/ResonanceFrame.py),
  `coord1/2`→`coords1/2`, `M`→`N` (clouds/CloudBasic.py), `N` unassigned→`N=len(...)`, `shift`→`shifts[i]`
  (clouds/NoeRelaxation.py), `file`→`self.fileName` (format/converters/SparkyFormat.py),
  `setMessage`→`print` (eci/CompletenessCheck.py).
- **REGRESSION found+fixed this session (IMPORTANT LEARNING):** prior session added a MODULE-level
  `from cambridge.dangle.DangleGui import DangleGui` to dangle/DangleFrame.py to fix an F821 (`DangleGui` used at
  `testMacro`), but DangleGui already imports DangleFrame → **import cycle** → 5 modules failed to import
  (DangleGui, DangleFrame, analysis/AnalysisGui, analysis/AnalysisPopup, extendNmr/ExtendNmrGui; smoke 1643→1638).
  **Fix:** move the import to a function-local inside `testMacro` (lazy) — breaks the cycle, keeps the dependency
  direction (popup→frame). **Rule: when an F821 import fix risks a module-level cycle, use a function-local import.**

### Part 2 — REMAINING ~342 findings (resume map; audit file-by-file in MAIN session, no subagents)
Run: `.venv/bin/ruff check ccpnmr2.5/python/ --select F821 --output-format concise`
- **NOT-WORTHWHILE cluster (≈171, by-design / out-of-import-surface):** CING standalone pipeline/demo scripts —
  `cing/Scripts/CASD/plot3.py` (148: `plt`/`np`/`results`/`NTvalue`/`rmsdToTarget`/`dataPath`),
  `cing/Scripts/interactive/axesRoutinesModified.py` (21: `np`/`mlab`/`datetime`/`mpath`/`mpatches`/`iq`),
  `cing/Scripts/Analysis/mouseBuffer.py` (1: `top`), `cing/PluginCode/test/parametersTest.py` (1: `refineParameters`).
  These import fine in smoke (undefined names are latent in function bodies); not core library modules.
  **Decision: document as by-design, no fix** (consistent with the 83 BY-DESIGN allowlist philosophy).
- **GENUINE candidates (confirmed, to fix in part 2):** `ccpnmr/format/converters/AmberFormat.py`
  (uses py2-only `sys.exc_type`/`sys.exc_value` + unimported `traceback`,`sys` in a bare `except:` — real crash),
  `ccpnmr/clouds/FilterClouds.py` (`math` unimported — likely real).
- **LIKELY FALSE-POSITIVE (document, don't touch):** metamodel dynamic attrs, web-server globals
  (`webServer/webFc.py`: `formatConvert`/`regenerateThisPage`/`createThisPageFirstTime` = FCGI-injected),
  `self` in generator/lambda bodies (`peaksIO.py`, `StructureIo.py`), `except (ValueError, v)`-style unpacks
  (`chemShiftsIO.py`). **Verify each against the actual code before deciding.**
- **Backups / out-of-scope (by-design, likely skip):** `utrecht/haddock/HaddockExportParam_new.py` (7, `_new` = backup file).
- **METHOD per file:** read the flagged line + surrounding fn; is the name bound at runtime (metaclass/exec/`__dict__`/
  C-ext/star-import) → FALSE-POSITIVE; is it a missing import / py2 builtin / typo / py2-only API → GENUINE (fix, then gates);
  is it a standalone script → NOT-WORTHWHILE. After ANY gate-relevant fix: full gate set (smoke+compile+pytest).

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

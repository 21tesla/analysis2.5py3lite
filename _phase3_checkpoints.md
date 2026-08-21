# Phase 3 Checkpoints — quality pass
# (resume-safe map; update at every session boundary)

## Status (2026-08-21)
- Baseline (after Phase 2c close-out `0c9182e`):
  smoke **1643 OK / 0 FAILED / 83 BY-DESIGN**; compile 0 errors; pytest 15 pass / 14 skip / 10 fail
  (the 10 fail = `/home/logan/software/testdata/` missing — pre-existing, deferred).
- Ruff baseline: **26,609 findings** (E,F,W,I,UP; auto-gen api/xml dirs excluded).
- Scope (user, 2026-08-21): (1) safe ruff fixes ONLY — no bulk style pass;
  (2) F821 undefined-name audit → fix genuine bugs. All other buckets DEFERRED (see below).
- **Current progress (2026-08-21):** Bucket 1 ✔ `65f475e`; Bucket 2 part 1 ✔ `bf290e6` (43 genuine fixes);
  Bucket 2 part 2a ✔ `8af5e45` (17 files); Bucket 2 part 2b ✔ this session (16 files:
  missing imports/in-scope-name fixes/two lazy imports for cycle-risk — each name target verified
  against its source, gates **1643/0/83, 0 compile, pytest 15/14/10**).
  F821 residue now **285**: ~170 in the by-design CING standalone-script cluster (plot3 ×148,
  axesRoutines ×21, mouseBuffer ×1), ~115 across ~44 files still to triage file-by-file.
  See Bucket 2 part-2 section for the resume map.

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

## Bucket 2 — F821 undefined-name audit (baseline 420) — status: **PARTS 1+2a+2b COMMITTED (43+17+16=76 fixed); ~285 RESIDUE (≈170 by-design CING, ≈115 to triage)**
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

### Part 2 — RESIDUE MAP (285 findings after part-2b; audit file-by-file in MAIN session, no subagents)
Run: `.venv/bin/ruff check ccpnmr2.5/python/ --select F821 --output-format concise`
**part-2a committed `8af5e45` (2026-08-21):** 17 files — incl. the two previously-confirmed GENUINE
candidates: `ccpnmr/format/converters/AmberFormat.py` (py2 `sys.exc_type`/`sys.exc_value` → `sys.exc_info()`;
added missing `sys`/`traceback` imports) and `ccpnmr/clouds/FilterClouds.py` (`math.sqrt` → `sqrt` —
module already has `from math import sqrt`). Full list in commit message.
**part-2b committed this session (2026-08-21):** 16 files, 31 findings — each name verified against
source before commit (in-scope variable / class existence / callback semantics):
- `cambridge/bayes/kmeans.py` — `vstack` → `numpy.vstack` (`import numpy` present; `__main__` block).
- `cambridge/dangle/DangleFrame.py` — `readGLE(path)` → `readGLE(path, resNum)` + signature `(gleFile, resNum=None)`
  (resNum was read in the missing-file error message → guaranteed NameError).
- `cambridge/isd/CCPNReader.py` (4) — `make_isd_residue(..., index)` + caller `enumerate(...)`;
  `get_volume(constraint, restraint_number)`/`% ccpn_restraint_number` → `constraint.serial` (both match the
  identical sibling call patterns in the same class); `R.restraints = restraints` → `restaints` (loop var).
- `cambridge/wms/ExtendNmrFrame.py` (2) — top-level imports were commented out (cycle-risk) → added
  function-local imports at use sites: `ccpnmr.eci.EntryCompletionFrame.EntryCompletionFrame` (class exists, line 988),
  `utrecht.haddock.HaddockFrame.HaddockFrame` (class exists, line 134).
- `cambridge/wms/RepositoryProperties.py` — undefined `user` (py2-only builtin) → `username` (defined line 121).
- `cambridge/wms/Task.py` — undefined `input`/`output` → `None` placeholders (attribute slots filled later via metamodel).
- `cambridge/wms/WorkflowFrame.py` — `WorkflowFramePopup` → `WorkflowPopup` (actual class, line 29; `__main__` block).
- `ccp/format/ansig/AnsigSpectrum.py` — bare `readParFile(parFile)` → `spectrum.readParFile(parFile)` (method exists, line 46).
- `ccp/format/pronto/peaksIO.py` — `fout.close()` in `write()` where `fout` never existed (whole method is commented-out
  dead code) → commented the orphan line to match.
- `ccp/gui/ViewChemCompVarFrame.py` — `Geometry.vectorsSubtract/Add` → bare `vectorsSubtract/Add`
  (module imports them from `memops.universal.Geometry`, line 62; `Geometry` itself never imported).
- `ccp/gui/ViewRamachandranFrame.py` — undefined `find_mean_sd` → added module-level helper
  (mean + population-SD of (phi, psi) pair columns; `import math` present line 61; call sites lines 776/789).
- `ccpnmr/analysis/core/SpinSystemTyping.py` — `current` → `cc` (defined line 295 `cc = getNewClassifications(cc0, num)`).
- `ccpnmr/clouds/CloudHomologueAssign.py` — `chain.residues` → `chainH.residues` (fn param is `chainH`);
  `getAtomSetCoords(atomSet, ...)` → `atom.atomSet` (loop var). **Adjacent pre-existing bug OBSERVED, NOT in-scope (F821 only):**
  `amideCoords.append(coords[0], residue)` — list.append with 2 args = guaranteed TypeError; needs `append((coords[0], residue))`.
- `ccpnmr/clouds/CloudThreader.py` (2) — `dshiftList` → `shiftList` (defined line 537);
  `XmlIO.loadProject` → function-local `from memops.general.Io import loadProject` (memops I/O API;
  `showWarning=showWarning` is CORRECT — `loadProject` docstring: showWarning is a callback fn (title,message));
  2 F821 remain in this file (separate names).
- `ccpnmr/clouds/PseudoResonances.py` — `mergePseudoSpinSysts(ss1, ss2)` → `self.mergePseudoSpinSysts(ss1, ss2)`
  (only caller, line 74) + signature `name` → `name=None`.
- `ccpnmr/format/converters/DataFormat.py` (5, 1 remains) — `structure.__class__` → `type(self.chemCompVar)` (invalid-object type);
  `peak.sortedPeakDims()` → `self.peakList.sortedPeaks()[0].sortedPeakDims()[0].dataDimRef` (mirrors the else-branch);
  `coordChain`/`coordResidue` ×4 → `chain`/`residue` (loop vars in scope); `% (title, status, year)` → `(title, className, year)`
  (loop var is `tcitation`; `className` is the in-scope param).
- **NOT-WORTHWHILE cluster (≈170, by-design / out-of-import-surface):** CING standalone pipeline/demo scripts —
  `cing/Scripts/CASD/plot3.py` (148: `plt`/`np`/`results`/`NTvalue`/`rmsdToTarget`/`dataPath`),
  `cing/Scripts/interactive/axesRoutinesModified.py` (21: `np`/`mlab`/`datetime`/`mpath`/`mpatches`/`iq`),
  `cing/Scripts/Analysis/mouseBuffer.py` (1: `top`), `cing/PluginCode/test/parametersTest.py` (1: `refineParameters`).
  These import fine in smoke (undefined names are latent in function bodies); not core library modules.
  **Decision: document as by-design, no fix** (consistent with the 83 BY-DESIGN allowlist philosophy).
- **Per-file residue map (285 total; counts after 2b, excluding the 170 CING by-design):**
  `cyana2ccpn/CyanaParser/CyanaParser.py` 12, `ccp/format/marvin/peaksIO.py` 11,
  `utrecht/haddock/HaddockExportParam_new.py` 7 (backup file), `ccp/format/pistachio/chemShiftsIO.py` 7,
  `ccpnmr/format/webServer/webFc.py` 5 (FCGI-injected), `pdbe/nmrStar/IO/NmrStarHandler.py` 4,
  `memops/format/xml/Compatibility.py` 4, `ccpnmr/format/gui/AcqProcParsEditPopup.py` 4,
  `ccpnmr/clouds/CloudThreaderPopup.py` 4, `pdbe/adatah/Pdb.py` 3, `memops/universal/Geometry.py` 3,
  `gottingen/PalesFrame.py` 3, `ccpnmr/integrator/plugins/Rosetta/write.py` 3,
  `ccpnmr/integrator/core/Io.py` 3, `ccpnmr/clouds/NoeMatrix.py` 3, `ccp/lib/StructureIo.py` 3,
  `pdbe/software/Util.py` 2, `pdbe/deposition/dataFileImport/formatConverterWrapper.py` 2,
  `pdbe/adatah/Io.py` 2, `memops/general/Util.py` 2, `grenoble/BlackledgeModule/BlackledgeModuleFrame.py` 2,
  `cing/Libs/helper.py` 2, `ccpnmr/clouds/CloudThreader.py` 2,
  and 22 files × 1 (incl. `ccpnmr/format/converters/DataFormat.py`, `ccpnmr/analysis/core/StructureBasic.py`,
  `memops/metamodel/ModelTraverse_py_2_1.py`, `ccp/lib/DataConvertLib.py`, `ccpnmr/workflow/Fc.py`,
  `cambridge/wms/ProtocolFrame.py`).
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

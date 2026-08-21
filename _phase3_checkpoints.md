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
  part 2a ✔ `8af5e45` (17 files); part 2b ✔ `3486cbb` (16 files); part 2c + 2d/e ✔ this session
  (29 files: CyanaParser, PalesFrame, Io-core, CloudThreaderPopup, ProtocolFrame, MarvinFormat,
  NmrStarFormat, Conversion, DataShifter ×2, DataFormat, CloudThreader, NmrpipeTableFormat,
  Talos Io/Util, Fc, ModelTraverse_py_2_1, Output, CASD/Util, Bmrb, listCcInfo, NoeMatrix,
  software/Util, formatConverterWrapper, memops/general/Util, cing helper, BlackledgeFrame,
  adatah/Io, Compatibility) — gates **1643/0/83, 0 compile, pytest 15/14/10** throughout.
  F821 now down to **229**, 100% classified: 172 CING standalone-scripts (by-design) +
  34 functional-but-idiomatic FPs (documented) + 23 incomplete-legacy (documented).
  **Bucket 2 (F821 audit) is COMPLETE** — see final-residue table in the part-2 section.

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

## Bucket 2 — F821 undefined-name audit (baseline 420) — status: **COMPLETE (191 fixed across parts 1+2a+2b+2c+2d/2e); final residue 229, all classified (172 CING by-design / 34 documented-FP / 23 incomplete-legacy)**
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
**parts 2c/2d/2e this session (29 files, all residue triaged to completion):**
- `cyana2ccpn/CyanaParser/CyanaParser.py` (12) — `path()` was CING's `cing.Libs.disk.Path` class
  (has exists/isdir/split3/iter — all used) → `from cing.Libs.disk import Path`;
  `protFile`/`f` ×2 → `finalProtFile`/`line` (loop var); `nTerror` → `ntu.nTerror`.
- `gottingen/PalesFrame.py` (3) — `widget.setup(names, values, index)` index undefined →
  `try: index = values.index(value) except ValueError: index = 0` (both getters); `FileType` → imported from `memops.gui.FileSelect`.
- `ccpnmr/integrator/core/Io.py` — `os.path.isfile(propfile)` → `propFile` (typo; line above).
- `ccpnmr/clouds/CloudThreaderPopup.py` — dead `if coord1 and coord2:` block after `continue` (unreachable) → commented out.
- `cambridge/wms/ProtocolFrame.py` — `ProtocolFramePopup` → `ProtocolPopup` (actual class; `__main__`).
- `ccpnmr/format/converters/MarvinFormat.py` — `peakName` assignment was commented out (use still live) → restored the assignment.
- `ccpnmr/format/converters/NmrStarFormat.py` — `newMethod(name=meth.methodLabel)` → `name=methName` (established local).
- `ccpnmr/format/general/Conversion.py` — `dataSource.serial` → `peakList.dataSource.serial` (owner of the peakList).
- `ccpnmr/format/gui/DataShifter.py` (2) — `newObjects.append(...)` result unassigned but printed as `newObject` →
  captured into `newObject`; print's `ccpnObject` → `refObject` (loop var in `setPresetLinks`).
- `ccpnmr/format/converters/DataFormat.py` — `setKeywords` loop body `keyword = rawCitation.keywd` was a broken
  no-op stub (rawCitation nowhere defined) → neutralized to `pass` (zero behavior change; used to crash on call).
- `ccpnmr/clouds/CloudThreader.py` — dead code after `return score` in `getSpinSystemScore` (p/Log undefined) → commented out.
- `ccpnmr/integrator/core/NmrpipeTableFormat.py` — `% res` → `% residue` (param).
- `ccpnmr/integrator/plugins/Talos/Io.py` — `atomNames=atomNames` → `atomNames=self.IOkeywords["atomNames"]` (per the guard).
  **Adjacent pre-existing latent bug OBSERVED (NOT in F821 scope):** line 87 `x.reaonance` (typo of `resonance`).
- `ccpnmr/integrator/plugins/Talos/Util.py` — `% res` → `% chemCompVar` (param).
- `ccpnmr/workflow/Fc.py` — `project.findFirstNmrExpPrototype` → `self.ccpnProject.findFirstNmrExpPrototype`
  (established accessor from `fc__init__`).
- `memops/metamodel/ModelTraverse_py_2_1.py` — `partitionRole[1]` → `partitionRoles[1]` (error-message only).
- `memops/universal/Output.py` — `__main__` demo: `from memops.universal.PostScript import PostScript` (class exists there).
- `nijmegen/CASD/Util.py` — `print(ss)` debug leftover (ss nowhere defined) → removed.
- `pdbe/adatah/Bmrb.py` — `chemCompPath = curChemCompRepository` → `chemCompArchiveDataDir`
  (+import; sibling adatah modules Pdb.py/NmrStar.py both use exactly this constant).
- `pdbe/chemComp/check/listCcInfo.py` — `chemAtomSysName.sysName` → `tempCasn.sysName` (loop var).
- `ccpnmr/clouds/NoeMatrix.py` (2) — `threshold` ×2 → `weightingFactor` (fn param; matches `ppm_weightingfactor` local).
  `symmat_w1` (isWatergate branch) left as incomplete-legacy (see final-residue table).
- `pdbe/software/Util.py` (2) — `getApplResNames`/`getNameInfo` live in `ccpnmr.format.general.Util`
  (checked: no import cycle) → function-local import at use site.
- `pdbe/deposition/dataFileImport/formatConverterWrapper.py` (2) — `existingForceChainMappings` first-use before
  assignment in the format loop (across-format carry-over intent) → `existingForceChainMappings = {}` before the `for importFormatName` loop.
- `memops/general/Util.py` (2) — `% (xx.__class__, id(xx))` → `(obj.__class__, id(obj))` (loop var `obj`).
- `cing/Libs/helper.py` (2) — `IPython` referenced without import → `import IPython` inside the try
  (ImportError now caught by the designed bare-except → warning + None, per docstring "Returns None on error").
- `grenoble/BlackledgeModule/BlackledgeModuleFrame.py` (2) — `moduleBvFileGood` read undefined →
  initialized to `False` (mirrors the sibling Pdb-file fn); `os.path.isfile(modulePdbFile)` → `chosenBvFile`
  (typo; matches sibling `modulePdbFile` pattern). **Deeper pre-existing defect OBSERVED (NOT in-scope):**
  nested `yes()`/`cancel()` set their own locals (no `nonlocal`) so the flag never propagates — same in the
  sibling Pdb fn; predates py3.
- `pdbe/adatah/Io.py` (2) — py2 `type(value) == file` → `hasattr(value, "read")` (file-like check; live code via
  Aria workflow `MultipartPostHandler`); py2 2-arg `raise TypeError(msg, traceback)` + undefined `sys` →
  plain `raise TypeError(msg)`. **Latent py2 remnants OBSERVED (NOT in F821 scope, documented):**
  `urllib.urlencode` (py2; py3 = `urllib.parse.urlencode`), `request.add_data`/`has_header` (urllib2 API) —
  the entire `MultipartPostHandler` is a py2 recipe copy; only reached on Aria multipart upload.
- `memops/format/xml/Compatibility.py` (4) — `linkCodes.append(x.linkCode)` (py2 comprehension-scope leak) →
  `linkCodes.append(linkEnds[0].linkCode)` (guarded by len==1 — the unique match, per intent);
  `loadError` prefix constant missing → added `loadError = "Cannot load XML: "` before `getVersion`;
  `for nameMapping in nameMapping:` → `for nameMapping in nameMappings:` (list defined at line 2170).

### FINAL RESIDUE TABLE (229 findings, all classified — bucket 2 audit COMPLETE)
Run: `.venv/bin/ruff check ccpnmr2.5/python/ --select F821 --output-format concise`
**A) CING standalone scripts (by-design, 172; consistent with 83-module BY-DESIGN allowlist):**
  `cing/Scripts/CASD/plot3.py` 148, `cing/Scripts/interactive/axesRoutinesModified.py` 21,
  `cing/Scripts/Analysis/mouseBuffer.py` 1, `cing/PluginCode/test/parametersTest.py` 1,
  `nijmegen/cing/CingFrame.py` 1 (`structure` in `getResiduesString` — structure selection ambiguous; CING GUI plugin).
**B) Functional-but-idiomatic FALSE-POSITIVES (working code ruff can't see, 34; do not touch):**
  - `globals()[var] = value` dynamic dispatch (names resolve as module globals at call time):
    `ccp/format/pistachio/chemShiftsIO.py` 7 (columnInfo loop), `ccpnmr/format/gui/AcqProcParsEditPopup.py` 4
    (phase0/phase1 loop; latent edge: refDims-falsy path).
  - `webFc.py` 5 — FCGI-injected globals (`formatConvert`, `regenerateThisPage`, `createThisPageFirstTime`).
  - loop-carried names (bound in a prior iteration / else-branch; first use guarded):
    `memops/universal/Geometry.py` 3 (`n`/`p` in `matrixMultiply` — **verified works at runtime**),
    `ccp/lib/StructureIo.py` 3 (`tlc`/`s`/`ins` TER branch — first use guarded by `lastChainCode is not None`),
    `ccp/format/cns/generalIO.py` 1 (`prevSearch` — assigned at loop end line 119; only read from 2nd iteration).
  - `ccp/format/marvin/peaksIO.py` 11 — broken generated `Dump_assignment` class (missing `def`, dead utility;
    names latent in never-instantiable class).
**C) INCOMPLETE-LEGACY (never complete in 2.5 distribution; NOT-WORTHWHILE, 23; document, do not touch):**
  - `utrecht/haddock/HaddockExportParam_new.py` 7 (`_new` = backup file).
  - `pdbe/nmrStar/IO/NmrStarHandler.py` 4 (`getSequence` needs `read.readStarDict`/`definitions.*` APIs absent from the tree).
  - `pdbe/adatah/Pdb.py` 3 (`mergeNmrStarWithPdb` fn + `validCodeSearch`/`ignorePdbCodes` helpers absent from tree).
  - `ccpnmr/integrator/plugins/Rosetta/write.py` 3 (`oldwrite` under "Old code from here on"; zero callers; Talos import commented out).
  - `ccpnmr/integrator/core/Io.py` 2 (`runSingleInteractive.else` — `nmrCalcRun`/`convert` never wired up).
  - `ccp/lib/DataConvertLib.py` 1 (`getStdResNameMap` uses `DataMapper.molTypeOrder`; DataMapper absent; fn has zero callers).
  - `pdbe/chemComp/modify/addSubstituent.py` 1 (`setChemAtomSetLinks`; caller at line 792 also passes no args — doubly broken legacy).
  - `pdbe/chemComp/carbo/makeFullSugar.py` 1 (`bondDict` bondType→priority helper absent from tree).
  - `ccpnmr/clouds/NoeMatrix.py` 1 (`symmat_w1` in interactive-watergate debug branch; never defined; has raw `input()`).
- **Per-file residue map BEFORE final close (285 total; counts after 2b, excluding the 170 CING by-design) — superseded by the table above:**
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

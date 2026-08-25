# Codebase Simplification Plan — CCPN 2.5 py3-lite

Status: **IN PROGRESS** (started 2026-08-24)
Target repo: **`github.com/21tesla/analysis2.5py3lite`** (new repo — one commit + push per stage)
Source repo (read-only reference, local remote `original`): `github.com/21tesla/analysis2.5py3`

## Goal

Strip 14 legacy/tool modules out of the CCPN analysis app (already migrated
py2→py3) and everything that only exists to support them, without breaking the
kept core (peak analysis, assignment, resonance, format conversion, ISD,
prodecomp, auremol, memops, tests).

**Remove (menu items + underlying code):**

| Menu | Item | Command (AnalysisPopup.py) | Implementation |
|---|---|---|---|
| Data Analysis | Heteronuclear NOE | `calcHeteroNoe` | `ccpnmr/analysis/popups/CalcHeteroNoe.py` |
| Data Analysis | 3J H-Hα Coupling | `calcHnHaCoupling` | `ccpnmr/analysis/popups/CalcHnHaCoupling.py` |
| Data Analysis | PALES: Alignment and RDCs | `pales` | `gottingen/` (4 files) |
| Data Analysis | MODULE: Alignment and RDCs | `blackledge_module` | `grenoble/BlackledgeModule/` (4 files) |
| Structure | Structure Viewer | `viewStructure` | `ccpnmr/analysis/popups/ViewStructure.py` |
| Structure | Make H Bond Restraints | `makeHbonds` | `ccpnmr/analysis/popups/MakeHbondRestraints.py` |
| Structure | DANGLE: Predict Dihedrals | `startDangle` | `cambridge/dangle/` (+670 data files) |
| Structure | ARIA: Structure calculation | `startAria` (already `pass` stub) | `paris/aria/` |
| Structure | CYANA (3 items + Cyana submenu) | `setupCyanaCalculation` / `importCyanaData` / `runCyana2Ccpn` | `cyana2ccpn/`, `ccpnmr/analysis/macros/MultiStructure.py`, `ccpnmr/integrator/plugins/Cyana/` |
| Structure | HADDOCK: Structure Docking | `startHaddock` | `utrecht/` |
| Structure | MECCANO: Structures from RDCs | `meccano` | `grenoble/meccano/` + C ext `ccpnmr2.5/c/other/meccano/` + `Meccano*.so` |
| Structure | PyRPF: Validate Peaks vs Structure | `startPyRPF` | `rutgers/` (2 files) |
| Structure | CING: Validate Structures | `submitCing` | `cing/` (340 .py) + `nijmegen/cing/` (4 files) |
| Structure | ECI: Database Deposition | `startECI` | `ccpnmr/eci/` (12 files, **except ReadPdb.py — see hazards**) |
| Structure | Secondary Structure Chart *(dup of Chart menu)* | `secStructureGraph` (kept — Chart menu still uses it) | menu entry only |
| Structure | Ramachandran Plot *(dup of Chart menu)* | `plotRamachandran` (kept — Chart menu still uses it) | menu entry only |

Note: "Secondary Structure Chart" and "Ramachandran Plot" entries in the
**Structure** menu are removed, but `secStructureGraph`/`plotRamachandran`
methods and `SecStructureGraph.py` / `ViewRamachandran.py` stay — the **Chart**
menu still uses them.

Also removed (user decision 2026-08-24): the standalone apps
`extendNmr/`, `cambridge/wms/`, `pdbe/deposition/` and their console entry
points (`ccpnmr-extend-nmr`, `ccpnmr-deposition`, `ccpnmr-eci`).

## Locked decisions (2026-08-24, with user)

1. **Workflow:** new GitHub repo `analysis2.5py3lite`; one commit + push per
   stage; no pushes to the original `analysis2.5py3` repo.
2. **Structure Viewer buttons:** `viewStructure` is called from 3 KEPT
   places — `ccpnmr/analysis/frames/PeakTableFrame.py`
   (`showStructConnections`, `showAllStructConnections`),
   `ccpnmr/analysis/frames/WindowFrame.py` (~L6588),
   `ccpnmr/analysis/popups/CalcShiftDifference.py` (~L464). → **remove those
   buttons/call sites** with the module.
3. **Standalone apps:** `extendNmr/`, `cambridge/wms/`,
   `pdbe/deposition/` — **remove entirely** (they hard-import CING, ARIA,
   ECI, HADDOCK).
4. **Not on the removal list — keep by default:** `nijmegen/CASD/`
   (orphaned by CING removal) and the Dangle/Haddock **metamodel** files
   (`cambridge/api/Dangle.py`, `cambridge/xml/Dangle.py`,
   `utrecht/api/Haddock.py`, `utrecht/xml/Haddock.py`,
   `ccpnmr2.5/model/cambridge/xml/Dangle/`,
   `ccpnmr2.5/model/utrecht/xml/Haddock/`) — generated `ccp.api.*` and
   `molsim` code lazily references them; removing them requires scrubbing
   hundreds of generated lines.
5. **Per-stage checkpoint policy (user, 2026-08-24):** every stage MUST end
   with a commit + push checkpoint before work on the next stage starts —
   code change and the `SIMPLIFICATION_PLAN.md` status/log update in ONE
   commit (keeps the "one commit per stage" revert unit intact). Rationale:
   sessions hit token limits; pushed checkpoints make a fresh session's
   resume trivial and no completed stage work is ever lost.
   **Push gotcha (hit 2026-08-24):** local branch is `master`, remote
   default is `main` — push with
   `git push <url> refs/heads/master:refs/heads/main` (one-shot URL-embedded
   token auth), not a bare `git push`, and keep the remote single-branch.

## Hazards (verified 2026-08-24 against source)

1. **`ccpnmr/eci/ReadPdb.py` is a KEPT feature** — powers
   *"Import PDB 3.20"* (`AnalysisPopup.py:129` import, `:2944` use) and
   `ccpnmr/format/converters/PdbFormat.py` chain. Relocate it **before**
   deleting `ccpnmr/eci/` (Stage 10).
2. **`ccpnmr/integrator/core/Io.py:68`** top-level
   `from cyana2ccpn.cyana2ccpn import importFromCyana` — hard break; plus
   CYANA-only helpers ~L1088-1130 (Stage 4).
3. **Shared KEPT code used by removed modules — do NOT delete:**
   `ccp/util/NmrCalc.py` (kept importers: `EditCalculation.py`,
   `cambridge/isd/NmrCalcExchange.py`),
   `ccpnmr/analysis/popups/EditCalculation.py`,
   `ccp/gui/ViewStructureFrame.py` (subclassed by kept `ViewChemCompVarFrame`,
   `ViewIsotopomerFrame`),
   `ccpnmr/format/converters/{CnsFormat,PdbFormat}.py`,
   `ccp/lib/{MoleculeQuery,StructureIo}.py`, `ccpnmr/analysis/core/*`, all of
   `ccp/api`, `memops`.
4. **Name-collision traps:** KEPT `ccp/format/aria/` (format parser, used by
   `ccpnmr/format/converters/AriaXmlFormat.py`) and KEPT
   `ccpnmr/integrator/plugins/Aria/` are **not** the removed `paris/aria`
   feature — leave them.
5. **`cambridge/` keeps most content:** only `dangle/` (and its data) is
   removed; `api/`, `bayes/`, `isd/`, `c/`, `xml/` stay.
6. **`nijmegen/` keeps `CASD/`** (decision 4); only `cing/` subdir removed.
7. **`import_smoke.py` allowlist:** ~45 `cing.*` "by design" entries die with
   the package (Stage 9) — remove them then or they become dangling.
8. **`setup.py`:** Meccano C-ext + GSL resolver logic must go in Stage 7,
   else builds referencing `c/other/meccano/` break. Release scripts
   (`scripts/*release*.sh`) gate on GSL/Meccano.
9. **`gui_boot_test.py` APPS list:** `cci`, `dangle`, `deposition`,
   `extend-nmr` entries must be removed in their stages or the gate fails.
10. **`pyproject.toml`:** `[project.scripts]` (`ccpnmr-eci`,
    `ccpnmr-dangle`, `ccpnmr-deposition`, `ccpnmr-extend-nmr`),
    `[tool.setuptools.packages.find] include`
    (`cing*`, `cyana2ccpn*`, `gottingen*`, `grenoble*`, `paris*`, `rutgers*`,
    `utrecht*`, `extendNmr*`, `nijmegen*`), and
    `[tool.ruff.lint.isort] known-first-party` references removed packages.

## Baseline (recorded 2026-08-24, pre-removal)

| Gate | Result |
|---|---|
| `python import_smoke.py` | exit 0 — **1729 modules**, 33 pre-existing "unexpected" (most are missing-optional-dep in `cing/*`, `cambridge/isd`, `cambridge/bayes`, `nijmegen/CASD`, `ccp/lib/Bmrb`, `ccp/util/V2Upgrade` + webServer `cherrypy`) |
| `python gui_boot_test.py` | **8/8 apps boot** (ccpnmr, eci, dangle, data-shifter, deposition, extend-nmr, format-converter, update) |
| `python -m pytest ccpnmr2.5/python/tests/` | **45 passed, 4 skipped** |

After each stage: `import_smoke` must stay exit 0 with **no NEW** unexpected
failures (counts may drop as we delete); `gui_boot_test` must stay green with
the app list shrinking as expected; pytest must not regress.
Python: anaconda `python` 3.13.5 (no `.venv`); `xvfb-run` available.

## Stages & status

| # | Scope | Status |
|---|---|---|
| 1 | Data Analysis: NOE, 3J, PALES, MODULE | ✅ 2026-08-24 |
| 2 | Standalone apps: `extendNmr/`, `cambridge/wms/`, `pdbe/deposition/` + scripts + boot entries | ✅ 2026-08-24 |
| 3 | ARIA: `paris/` + menu + methods | ✅ 2026-08-24 |
| 4 | CYANA: `cyana2ccpn/` + `macros/MultiStructure.py` + integrator `Cyana/` + `Io.py` import | ✅ 2026-08-24 |
| 5 | DANGLE: `cambridge/dangle/` + `ccpnmr-dangle` | ✅ 2026-08-24 |
| 6 | HADDOCK: `utrecht/haddock/` + menu + method | ✅ 2026-08-24 |
| 7 | MECCANO: `grenoble/` + `c/other/meccano/` + `setup.py` GSL/Meccano blocks | ✅ 2026-08-24 |
| 8 | PyRPF: `rutgers/` | ✅ 2026-08-24 |
| 9 | CING: `cing/` + `nijmegen/cing/` + smoke allowlist | ✅ 2026-08-24 |
| 10 | ECI: `ccpnmr/eci/` (relocated `ReadPdb.py` to `ccpnmr/format/converters/`) + `ccpnmr-eci` | ✅ 2026-08-24 |
| 11 | Structure Viewer + Make H Bond Restraints popups + remove kept-caller buttons (7 files) | ✅ 2026-08-24 |
| 12 | Cross-cutting sweep: `pyproject.toml`, `bin/`, release scripts, docs, extras, final verification | ✅ 2026-08-24 |

### Stage checklist detail

**Stage 1 — Data Analysis (NOE / 3J / PALES / MODULE)**
- `AnalysisPopup.py`: remove top imports L78-79; popupActions entries
  `calc_hnha_coupling`, `calc_hetero_noe`, `pales`, `blackledge_module`;
  `setDataMenu` 4 `add_command` blocks + 4 `menu_items` entries; method
  definitions `calcHnHaCoupling`, `calcHeteroNoe`, `pales`,
  `blackledge_module`.
- Delete: `popups/CalcHeteroNoe.py`, `popups/CalcHnHaCoupling.py`,
  `gottingen/` (4 files), `grenoble/BlackledgeModule/` (4 files; keep
  `grenoble/` itself — meccano stays until Stage 7).
- `pyproject.toml`: drop `gottingen*` include; isort first-party `gottingen`.

**Stage 10 — ECI** (highest risk): relocate `ccpnmr/eci/ReadPdb.py` to
`ccpnmr/format/` (or `popups/`), repoint `AnalysisPopup.py:129` + any other
importers, **then** delete remaining `ccpnmr/eci/*`.

### Rollback

Each stage is exactly one commit on the new repo → `git revert <sha>` restores
it cleanly.

## Stage log

**Stage 1 — Data Analysis (NOE, 3J, PALES, MODULE) — ✅ 2026-08-24**
- `AnalysisPopup.py`: removed top imports (CalcHeteroNoe, CalcHnHaCoupling),
  4 popupActions entries, 4 `menu.add_command` blocks, 4 `menu_items`
  entries, 4 command methods, and the dead commented
  `activatePales`/`activateModule` blocks.
- Deleted: `popups/CalcHeteroNoe.py`, `popups/CalcHnHaCoupling.py`,
  `gottingen/` (4 files), `grenoble/BlackledgeModule/` (4 files).
- `pyproject.toml`: dropped `gottingen*` package include + isort entry.
- Gates: import_smoke exit 0 (1719 modules; pre-existing 33 failures unchanged),
  gui_boot_test 8/8, pytest 45 passed / 4 skipped.
- Pitfall for later stages: labels like the former "3J H-Hα" hold a
  **literal** `\u03b1` escape in source; the edit tool's escaping pipeline
  can't reproduce it — use an assert-protected index-based Python deletion
  for such lines instead.

**Stage 2 — Recon recorded 2026-08-24 (edits NOT yet made; session ended
here, Stage 2 continues in the next session)**
- Delete: `extendNmr/` (14 files: `ExtendNmrGui.py`, `__init__.py`,
  `images/`), `cambridge/wms/` (30 files), `pdbe/deposition/` (5 files).
- **Keep the rest of `pdbe/`** — `pdbe.nmrStar`, `pdbe.adatah`,
  `pdbe.chemComp`, `pdbe.software`, `pdbe.general`, `pdbe.xml` are imported
  by KEPT modules (`ccp/format/nmrStar/generalIO.py`,
  `ccpnmr/integrator/plugins/NmrStar/Io.py`,
  `ccpnmr/format/general/Conversion.py`, `ccpnmr/workflow/{Aria,Cing}.py`,
  `ccpnmr/format/process/sequenceCompare.py`, `ccp/examples/help_doc/`, and
  `ccpnmr/eci/*` until Stage 10). In pyproject: drop only `extendNmr*`
  include + `extendNmr` isort entry + scripts `ccpnmr-extend-nmr` and
  `ccpnmr-deposition`; keep `pdbe*`, `cambridge*`, `nijmegen*`.
- Nothing outside those dirs imports `extendNmr` or `cambridge.wms`
  (grep-verified 2026-08-24) except `gui_boot_test.py` APPS entries
  `deposition` + `extend-nmr` (drop them → expect 6/6 apps to boot).
- `nijmegen/CASD/casdPipeLine.py:21` imports to-deleted
  `pdbe.deposition...FormatConverterWrapper`; that file is already in the
  33 baseline import-failures → no new failures expected.
- Expected gate deltas: import_smoke TOTAL ~1719 → ~1683 (−~36 modules),
  FAILED unchanged (33); gui_boot_test 8/8 → 6/6; pytest unchanged.

**Stage 2 — Standalone apps (extendNmr, cambridge/wms, pdbe/deposition) — ✅ 2026-08-24**
Recon re-verified pre-edit, identical to above: only external importers are
the already-broken `nijmegen/CASD/casdPipeLine.py` and in-dir cross-imports;
`cambridge/__init__.py` + `pdbe/__init__.py` are bare `pass`; no other
references in py/toml/cfg/in/sh/md (outside `dist/`+`build/`, which are
git-ignored stale copies — left untouched).
- Deleted (49 files, 9393 lines): `ccpnmr2.5/python/extendNmr/` (14),
  `ccpnmr2.5/python/cambridge/wms/` (30),
  `ccpnmr2.5/python/pdbe/deposition/` (5). Rest of `pdbe/` and `cambridge/`
  kept as planned.
- `pyproject.toml`: dropped scripts `ccpnmr-deposition` +
  `ccpnmr-extend-nmr`, package include `extendNmr*`, isort first-party
  `extendNmr`.
- `gui_boot_test.py`: dropped APPS entries `deposition` + `extend-nmr`.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1683** (1719−36, exactly as
    predicted), OK 1567, FAILED **33 unchanged** (30× cing `Sql`,
    2× `cherrypy`, 1× `psycopg2` — pre-existing classes), BY-DESIGN 83.
  - `gui_boot_test.py` **6/6** (8/8 − the 2 removed apps).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
- Residual (by design): `nijmegen/CASD/casdPipeLine.py` still imports the
  now-deleted `pdbe.deposition.*` — already among the 33 import failures, so
  no new failures; goes away with the CING removal in Stage 9.

**Stage 3 — ARIA (paris/, menu, methods) — ✅ 2026-08-24**
Recon (verified pre-edit): `paris/` = 5 tracked files (`paris/__init__.py`,
`paris/aria/__init__.py`, `AriaExtendNmrFrame.py`, `AriaRunFrame.py`,
`CcpnToAriaXml.py`, 1841 lines total incl. headers). Zero external importers
of any `paris.aria` symbol — the only external references were inside
`AnalysisPopup.py` (one commented import + two `pass` stubs). KEPT-by-design
and verified untouched: `ccp/format/aria/` (format parser),
`ccpnmr/integrator/plugins/Aria/`, `ccpnmr/workflow/Aria.py` (imports only
`ccpnmr.workflow.Util` — hazard 4 confirmed, no `paris` import). No `aria`/`paris`
entries in `gui_boot_test.py`, no references in `setup.py`, `MANIFEST.in`,
`import_smoke.py`, `bin/`, or `scripts/`.
- `AnalysisPopup.py`: removed commented import
  `#from paris.aria.AriaExtendNmrFrame import AriaPopup`, the
  "ARIA: Structure calculation" `menu.add_command` block (shortcut "A"),
  and both stub methods `activateAriaSetup` + `startAria` (each a
  `pass #self.openPopup("aria_setup", AriaPopup)`). `activateAriaSetup` was
  kept only "so old tutorial script works" — no tutorial script exists in
  the repo (grep-verified), and the plan removes the ARIA feature entirely.
- Deleted: `ccpnmr2.5/python/paris/` (5 tracked files; untracked
  `__pycache__` residue removed with the dir).
- `pyproject.toml`: dropped package include `paris*` and isort first-party
  `paris`.
- Residual (by design): `ccp/util/NmrCalc.py` header docstring still lists
  "paris/aria" among Oct-2012 consumers — stale historical note (it also
  lists Stage-1-removed BlackledgeModule); left untouched per "don't edit
  unrelated comments".
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1678** (1683−5, exactly the 5 removed
    `paris` modules), OK 1562, FAILED **33 unchanged**, BY-DESIGN 83.
  - `gui_boot_test.py` **6/6** (no ARIA app entry — unchanged).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).

**Stage 4 — CYANA (cyana2ccpn, MultiStructure, Cyana plugin, Io.py) — ✅ 2026-08-24**
Recon (verified pre-edit): hazard 2 confirmed — `ccpnmr/integrator/core/Io.py:68`
top-level `from cyana2ccpn.cyana2ccpn import importFromCyana` plus 6 CYANA-only
helpers calling it. `macros/MultiStructure.py` looked mixed (it also wraps
ROSETTA/UNIO/ASDP protocols) but ALL of its functions have zero callers repo-wide
(grep-verified) → whole file removable. Plugin discovery is LAZY per protocol
name (`intUtil.getIntegratorPlugin(protocolName)` — no startup dir scan), so
deleting `plugins/Cyana/` can't break the 9 kept plugins (Aria, Asdp, Cosmos, Isd,
MultiStruc, NmrStar, Rosetta, Talos, Unio). Protocol JSONs are reached only via
the deleted menu/MultiStructure paths (protocol names `CYANA_SS4`/
`CYANA_PEAKLIST`/`CYANA_SS1` appear nowhere else) → 3 JSONs removed. Other
"Cyana" hits in the tree are atom-naming-system names (`CYANA2.1` in ChemComp
data, `ccp/format/cyana/` format parser, `CyanaFormat.py` converter) — KEPT
features, same category as `ccp/format/aria/`.
- Deleted (13 tracked files): `ccpnmr2.5/python/cyana2ccpn/` (5: `__init__`,
  `cyana2ccpn.py`, `classes4.py`, `CyanaParser/{__init__,CyanaParser}`),
  `ccpnmr/analysis/macros/MultiStructure.py`, `ccpnmr/integrator/plugins/Cyana/`
  (4: `__init__`, `read.py`, `Util.py`, `write.py`), protocol JSONs
  `data/ccpnmr/integrator/{Cyana_SS4,CyanaWF,CyanaWF_SS}.json`.
- `ccpnmr/integrator/core/Io.py`: removed the top-level cyana2ccpn import and
  the contiguous 6-function CYANA block `runCyana2Ccpn`,
  `setupCyana2CcpnDialogue`, `runCyana2CcpnDialogue`,
  `setupPreviousCalculation` + `runPreviousCalculation` (the latter calls
  `importFromCyana` — hard break if kept; only caller was the deleted
  MultiStructure with protocolName CYANA_SS4), `importDataFromCyana`.
  Left untouched: `writeExecuteScript` (shared plugin template; its
  `cyanatable.txt` line is inert legacy for non-Cyana protocols).
- `AnalysisPopup.py`: removed Cyana submenu (3 items: Setup/Import/Run CYANA)
  + the `menu.add_cascade(label="Cyana")` line + methods
  `setupCyanaCalculation`, `importCyanaData`, `runCyana2Ccpn`.
- `pyproject.toml`: dropped package include `cyana2ccpn*` (no isort/scripts/
  smoke-test entries existed).
- Residuals (by design): `data/ccp/cyana/{cyana,cyana2.1}.lib` KEPT — read by
  kept `ccp/format/cyana/cyanaLibParser.py` (format converter). `import_smoke.py`
  has 1 CYANA allowlist entry under `cing.Database.Scripts.addCYANA2` — dies
  with cing in Stage 9.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1668** (1678−10, exactly the 10 removed
    `.py` modules), OK 1552, FAILED **33 unchanged**, BY-DESIGN 83.
  - `gui_boot_test.py` **6/6** (unchanged).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).

**Stage 5 — DANGLE (cambridge/dangle, menu, script, boot entry) — ✅ 2026-08-24**
Recon (verified pre-edit): `cambridge/dangle/` = 680 tracked files (10 `.py` +
670 data files — `.tab`, `Plot_*.int` etc., 918k lines total). All
`DangleGui`/`DangleFrame`/`cambridge.dangle.*` references are either in-package
or the 5 external touch points below. Pure-Python predictor (no external
binaries; no `ccpnmr/workflow/Dangle*`, no `cambridge/c/` dangle ext, no
MANIFEST.in entries). Dangle metamodel files (`cambridge/api/Dangle.py`,
`cambridge/xml/Dangle.py`, `ccpnmr2.5/model/cambridge/xml/Dangle/`) KEPT per
decision 4.
- Deleted: `ccpnmr2.5/python/cambridge/dangle/` (680 files — `DangleGui`,
  `DangleFrame`, `neuralNet/{DangleNN,NeuralNetwork}`, `src/{Predictor,Protein,
  Reference,dangle}` + inits + 670 data files).
- `AnalysisPopup.py`: removed top import `from cambridge.dangle.DangleGui
  import DangleGui`, popupActions entry `"dangle": self.startDangle`, the
  "DANGLE: Predict Dihedrals" Structure-menu block (shortcut "D"), and the
  `startDangle` method.
- `pyproject.toml`: dropped script `ccpnmr-dangle`.
- `gui_boot_test.py`: dropped APPS entry `dangle` (+ fixed the usage example in
  the module docstring `--apps ccpnmr,dangle` → `--apps ccpnmr,eci`).
- Residuals (by design, for Stage 12 cross-cutting sweep): `bin/dangle`,
  `bin/dangle2`, `bin/dangle2.5` shell launchers (plan puts `bin/` in Stage 12;
  note Stages 1-2 also left `bin/extendNmr*` + `bin/depositionFileImporter*`,
  and the release-script `need` sets in `scripts/{linux,macos}_release.sh`
  still list `ccpnmr-dangle` + `ccpnmr-deposition`/`ccpnmr-extend-nmr`);
  help-docstring link defs `.. _DANGLE: DangleGui.html` in
  `BrowseConstraints.py:275` + `MakeHbondRestraints.py:201` (inert text, Stage 1
  left the same for removed popups).
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1658** (1668−10, exactly the 10 removed
    `.py` modules), OK 1542, FAILED **33 unchanged**, BY-DESIGN 83.
  - `gui_boot_test.py` **5/5** (as predicted: 6/6 − the removed `dangle` app).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).

**Stage 6 — HADDOCK (utrecht/haddock, menu, method) — ✅ 2026-08-24**
Recon (verified pre-edit): the HADDOCK **feature** is
`ccpnmr2.5/python/utrecht/haddock/` = 17 tracked files (16 popups/API frames +
`__init__.py`, 10056 lines total: `HaddockPopup`, `HaddockFrame`, `HaddockApi`,
`HaddockBasic`, `HaddockDaniPopup`, `HaddockRdcPopup`, `EditSymmetryPopup`,
`SymmetryPopup`, `HaddockDnaRnaRest`, `HaddockExport{Classic,Param,Param_new}`,
`HaddockImportRunCns`, `HaddockLocal`, `HaddockServerUpload`, `APIexample`).
Single external importer: `ccpnmr/analysis/AnalysisPopup.py` (import + one
popup instantiation). No `ccpnmr/analysis/popups/Haddock*` files exist (the
popups live inside `utrecht/haddock/` itself), no `ccpnmr/workflow/Haddock*`
workflow module, no integrator plugin, no data/ protocol JSONs, no
`ccpnmr-haddock` console script, no `bin/` launcher, no `scripts/` or
`gui_boot_test.py` / `import_smoke.py` entries.
KEPT per decision 4 (verified live consumer): `utrecht/__init__.py`,
`utrecht/api/` (`Haddock.py` metamodel + `doc/` API HTML — same category as
Stage 5 keeping `cambridge/api/doc/`), `utrecht/xml/Haddock.py`
(`memops/xml/Implementation.py:282` does `import utrecht.xml.Haddock`),
`ccpnmr2.5/model/utrecht/xml/Haddock/`. Because the `utrecht` package
survives, `pyproject.toml`'s `utrecht*` package-include + isort `utrecht`
first-party entry STAY (unlike Stage 5, no pyproject change here).
- Deleted: `ccpnmr2.5/python/utrecht/haddock/` (17 files, 10056 lines).
- `AnalysisPopup.py`: removed `from utrecht.haddock.HaddockPopup import
  HaddockPopup`, the "HADDOCK: Structure Docking" Structure-menu
  `add_command` block (shortcut "K"), the `"HADDOCK: Structure Docking"`
  `menu_items[StructureMenu]` entry (list feeds `setMenuState`'s entry-count
  loop), and the `startHaddock` method.
- Residual (by design, Stage 12 sweep): `ccpnmr/workflow/Constants.py`
  `programList` still contains the string `"Haddock"` — generic external
  program-name dropdown list (Stage 3 left `"Aria"` in it for the same reason);
  `memops/xml/Implementation.py` keeps its `import utrecht.xml.Haddock`
  (metamodel, decision 4).
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1641** (1658−17, exactly the 17 removed
    `.py` modules), OK 1525, FAILED **33 unchanged**, BY-DESIGN 83.
  - `gui_boot_test.py` **5/5** (unchanged — no HADDOCK app entry ever).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).

**Stage 7 — MECCANO (grenoble/, C sources, setup.py GSL) — ✅ 2026-08-24**
Recon (verified pre-edit): after Stage 1 (BlackledgeModule removed)
`ccpnmr2.5/python/grenoble/` was MECCANO-only: 8 tracked files
(`__init__`×3, `MeccanoPopup.py`, `c/{copySharedObjs,copySharedObjs.bat,
linkSharedObjs}`, `meccano/data/phi_psi_database_loop_glysymm` — a 41,308-line
phi/psi data file). C sources: `ccpnmr2.5/c/other/meccano/` = 28 tracked
files (13 headers, 11 `src/*.c`, `pysrc/py_meccano.c`, 2 `meccano2_stat_ramaDB_fwd/*.c`,
3 Makefiles) — 36 files/48,654 lines total deleted. Single external
importer: `AnalysisPopup.py` `meccano()` method (lazy import inside a
try/except — no top-level import line). KEPT (verified): `ccpnmr/v2io/
Constants.py` `"GSL": ("other", "Gsl")` entries are **file-extension**
constants (unrelated to GNU GSL); `c/environment*.txt` `LINK_LIBRARIES/COPY_
LIBRARIES = ...SharedObjs` vars + `UpdateAgent.py` references are per-package
build helpers shared with KEPT `cambridge/c` (which has its own copies —
untouched); `ccp/util/NmrCalc.py:5` docstring lists
`grenoble/BlackledgeModule` (stale historical note, Stage-1/3 precedent);
`workflow/Constants.py` `programList` has no "Meccano" string; `c/Makefile`
has no meccano target; no import_smoke allowlist entries, no gui_boot APPS
entry, no tests, no `bin/` launcher, no integrator plugin, no data/ protocols.
- Deleted: `ccpnmr2.5/python/grenoble/` (8) + `ccpnmr2.5/c/other/meccano/`
  (28) + untracked build residue (`grenoble/c/Meccano*.so`, `__pycache__`).
  `c/other/` now contains only `cambridge/`.
- `AnalysisPopup.py`: removed "MECCANO: Structures from RDCs" menu block
  (shortcut "M"), the `menu_items[StructureMenu]` entry, and the `meccano()`
  method (incl. its lazy `from grenoble.meccano.MeccanoPopup import ...`).
- `setup.py`: removed `MEC` path const, the GSL resolver (`_gsl_usable` +
  `GSL = next(...)`), the FAM `"Meccano"` ext entry (last FAM member), and
  the "Meccano is the one OPTIONAL ext" skip-block. No GSL/MEC/meccano
  residue; parses clean. (hazard 8 closed — builds no longer reference
  `c/other/meccano/`.)
- `pyproject.toml`: dropped package include `grenoble*` + isort first-party
  `grenoble` (package fully gone — unlike `utrecht`, which kept living).
- `scripts/copy_cext.sh`: dropped `Meccano)` → `grenoble/c` case branch.
- `scripts/{linux,macos}_release.sh`: removed the `die`-gated GSL resolution
  block (loop + `CCP_GSL_PREFIX` export) + GSL prereq/prefix lines +
  "(incl. Meccano)" build echoes + header comment lines. `bash -n` clean.
  (hazard 8's release-gate part closed — no more hard GSL requirement.)
- `scripts/publish.sh`: dropped the 2 commented Meccano/GSL lines.
- Residual (by design, Stage-12 docs sweep): historical changelog line
  `ccpnmr/analysis/doc/Changes.html:145` ("distribute the C code... install
  GSL") — past release note, left untouched.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1637** (1641−4, exactly the 4 removed
    `grenoble` `.py` modules), OK 1521, FAILED **33 unchanged**, BY-DESIGN 83.
  - `gui_boot_test.py` **5/5** (unchanged — no MECCANO app entry ever).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).

**Stage 8 — PyRPF (rutgers/, menu, method) — ✅ 2026-08-24**
Recon (verified pre-edit): `ccpnmr2.5/python/rutgers/` = **3** tracked files
(the plan table said "2"; actual = `__init__.py`×2 + `rpf/PyRPF.py`, 2927
lines total — `PyRpfPopup(BasePopup)` RPF quality popup + `pyRpfMacro` /
`pyRpfApp` helpers). Single external importer: `AnalysisPopup.py` — a
**top-level** `from rutgers.rpf.PyRPF import PyRpfPopup` (NOT lazy like
Stage-7's meccano), the "PyRPF..." menu block, and the `startPyRPF` method.
No `popupActions` entry, no `menu_items[StructureMenu]` entry (simpler than
Stage 6). The `cing/` `PyRPF`/`RPF` files (`cing/PluginCode/RPF.py`,
`cing/Scripts/Analysis/PyRPF.py` + its test) are Stage-9 territory — kept.
KEPT (verified): `ccpnmr2.5/doc/acknowledgements.html:106` "Cathy Lawson
(Rutgers)" is a **person's name**, not the package; `pdbe/adatah/CasdNmr.py:48`
matches legacy data-file name strings (inert text, not an import). No
metamodel / `model/` / `ccp/api` / `ccpnmr/api` references, no `bin/`
launcher, no `scripts/`, no MANIFEST.in entries, no `import_smoke.py`
allowlist entries, no `gui_boot_test.py` APPS entry, no tests.
- Deleted: `ccpnmr2.5/python/rutgers/` (3 files, 2927 lines) +
  untracked `__pycache__` residue.
- `AnalysisPopup.py`: removed the `rutgers.rpf.PyRPF` top import, the
  "PyRPF: Validate Peaks vs Structure" Structure-menu `add_command` block
  (shortcut "F"), and the `startPyRPF` method.
- `pyproject.toml`: dropped package include `rutgers*` + isort first-party
  `rutgers` (package fully gone — Stage-7 rule).
- Residuals (by design, Stage-12 docs sweep): `Structure.rst` menu-doc line
  `PyRPF: Validate Peaks vs Structure <../popups/PyRpfPopup.rst>` (Stages
  3/5/6 left the DANGLE/ARIA/HADDOCK lines there too);
  `ccpnmr/analysis/doc/Changes.html:89` historical changelog line.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1634** (1637−3, exactly the 3 removed
    `rutgers` `.py` modules), OK 1518, FAILED **33 unchanged**, BY-DESIGN 83.
  - `gui_boot_test.py` **5/5** (unchanged — no PyRPF app entry ever).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).

**Stage 9 — CING (cijermen, nijmegen/cijermen, menu, smoke allowlist) — ✅ 2026-08-24**
Recon (verified pre-edit): `ccpnmr2.5/python/cijermen/` = 340 tracked .py
(+213 data files; 553 tracked) and `nijmegen/cijermen/` = 4 tracked (.py:
`CingPopup`, `CingFrame`, `iCingRobot`, `__init__`) — 557 tracked files total.
Single external importer: `AnalysisPopup.py` (top import `from
nijmegen.cijermen.CingPopup import CingPopup` + menu + `submitCing` +
`menu_items` entry). `ccpnmr/workflow/Cijmegen.py` (WebCing/CijmegenWorkFlow
— remote CING web-service client via HTTP) **KEPT**: Stage-3 precedent
(`workflow/Aria.py` survived ARIA removal), zero importers, outside the
stage-table scope. Cython `cijermen/Libs/cython/superpose` C-ext: ONLY
consumer was `cijermen/PluginCode/queeny.py` (`Rm6dist`, dies with cing);
kept code uses the pure-Python `memops/universal/Geometry.
superposeNewVectorsOnOld` (importers `pdbe/chemComp/.../addSubstituent.py`,
`ccp/util/Molecule.py` — verified). The release-script `[2/4]` superpose
build step (`CY=.../cijermen/Libs/cython` + Cython pip-install) now points
at a dead path — RESIDUAL for the Stage-12 release-script sweep (same
bucket as the S2/S5 deferrals; those scripts were already stale: their
FAILED=0 gate predates the 33-failure baseline and their `need` console-set
lists 4 since-deleted apps). `ccpnmr/format/webServer/{Util,webFc}.py` KEPT
(cherrypy consumers, not on the removal list — they ARE the 2 remaining
FAILED entries). KEPT (verified) inert strings: `nijmegen/CASD/*` `CING`
task-name strings (CASD pipeline, decisions 4/6),
`ccpnmr/workflow/Constants.py` programList "Cing" (S3/S6 precedent),
`ccp/util/Validation.py` `storeRogScores(context="CING")` +
`ViewStructure.py` `context == "CING"` ROG-score handling (kept feature),
`pdbe/adatah/Cijmegen.py` (kept pdbe module), `ccp/util/NmrCijmegen.py:4`
stale consumer docstring (S1/S7 precedent). No integrator plugin, no
data/ protocol JSONs, no setup.py/MANIFEST.in/gui_boot_test entries, no
`bin/` launchers, no tests.
- Deleted: `ccpnmr2.5/python/cijermen/` (340 .py + data) +
  `ccpnmr2.5/python/nijmegen/cijermen/` (4) — 557 tracked files +
  `__pycache__` residue. Includes the S8 PyRPF leftovers
  (`PluginCode/RPF.py`, `Scripts/Analysis/PyRPF.py`) and the S4 CYANA
  allowlist entry (`Database.Scripts.addCYANA2`).
- `AnalysisPopup.py`: removed the `from nijmegen.cijermen.CingPopup import
  CingPopup` import, the "CING: Validate Structures" Structure-menu block
  (shortcut "C"), the `submitCing` method, and the
  `menu_items[StructureMenu]` entry (`setMenuState` uses the list only as a
  `len(...)+4` try/except loop bound — inert; S6 dropped the HADDOCK entry
  the same way). KEPT in that list: "DANGLE: Predict Dihedrals" string (S5
  residual, Stage 12).
- `import_smoke.py`: removed **73** `cijermen.*` entries from
  KNOWN_NON_IMPORTABLE (EXTERNAL 19, PLUGIN 11, INTERACTIVE 14, DATA 25,
  ENV 2, VENDOR 3 — the PLUGIN/INTERACTIVE/DATA/VENDOR sections went empty
  and their category-legend comment lines were dropped too). **Hazard 7
  CLOSED.** 10 entries remain: cambridge.bayes, ccp.lib.Bmrb.bmrb,
  ccp.util.V2Upgrade, pdbe.software.vascoReferenceCheck,
  pdbe.chemComp.export.setLicenses, cambridge.isd, nijmegen.CASD×4.
- `pyproject.toml`: dropped package include `cijermen*` + isort first-party
  `cijermen` (package fully gone). KEPT: `nijmegen*`/`nijmegen` (CASD
  survives), the `[project.optional-dependencies] optional` block (Stage-12
  "extras" bucket — cherrypy/mako still consumed by kept webServer;
  sqlalchemy/psycopg2-binary/decorator now orphaned by kept code) and its
  "used by cing / web-server" comment (Stage 12).
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1290** (1634−344, exactly the 344
    removed modules), OK **1278**, FAILED **2** (33−31: the 30× `Sql`
    ImportWarnings + cing's `psycopg2` all gone; 2× `cherrypy` in KEPT
    `ccpnmr/format/webServer/*` remain), BY-DESIGN **10** (83−73).
  - `gui_boot_test.py` **5/5** (unchanged — no CING app entry ever).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).

**Stage 10 — ECI (ccpnmr/eci/, menu, script, boot entry) — ✅ 2026-08-24**
Recon (verified pre-edit): `ccpnmr2.5/python/ccpnmr/eci/` = 12 tracked files
(`AditMandFields`, `CompletenessCheck`, `EciShiftAnalysis`,
`EntryCompletion{Frame,Gui,Popup}`, `IsotopeLabeling`, `ReadPdb`,
`__init__` (bare `pass`), `_licenseInfo`, `makeAditMandDict`,
`nmrStarDictNew`). **Hazard 1 verified:** `ReadPdb.py` (1679 lines) powers the
KEPT "Import PDB 3.20" feature (`AnalysisPopup.py` top-level import at L126 +
`return ReadPdb(...)` at L2713) and is the ONLY eci module used by kept code.
Key finding: **zero `ccpnmr.eci.*` sibling imports** inside `ReadPdb.py` — all
its imports are kept packages (`ccp.api.general`, `ccp.format.pdb`,
`ccp.general.Io`, `ccpnmr.format.converters.PdbFormat`,
`ccpnmr.format.process.matchResonToMolSys`, `memops.api.Implementation`,
`memops.universal.Util`) → relocatable with **zero content change**;
git records the `git mv` as a rename. Verified: `PdbFormat.py` does NOT import
ReadPdb (it's a consumer, not a dependee — the hazard-1 "chain" is the
feature path, not an import edge). External refs to the rest of the package:
`AnalysisPopup.startECI` (lazy import) + the config/doc items below.
No `import_smoke.py` allowlist entries (grep-verified), no pyproject
package-include/isort entries (eci is a subpackage of KEPT `ccpnmr`), no
MANIFEST.in/setup.py references. KEPT (verified): `memops/general/license/
headers.py` `_licenseInfo` mentions are the generic license-strip mechanism
(not eci-specific).
- Relocated: `ccpnmr/eci/ReadPdb.py` → **`ccpnmr/format/converters/ReadPdb.py`**
  (plan's "relocate to `ccpnmr/format/`" — placed alongside `PdbFormat.py`,
  the class it wraps; all converter readers live there). Repointed
  `AnalysisPopup.py:126` to `from ccpnmr.format.converters.ReadPdb import
  ReadPdb` (kept isort order: after `NmrStarFormat`, before `ccpnmr.nexus`).
  The kept "Import PDB 3.20" call site (L2713) is untouched. **Hazard 1
  CLOSED.**
- Deleted: the remaining 11 `ccpnmr/eci/*` files (incl. `__init__.py`,
  `_licenseInfo.py`) + untracked `__pycache__` residue; `ccpnmr/eci/` dir
  gone.
- `AnalysisPopup.py`: removed the `popupActions` entry
  `"entry_completion_interface": self.startECI`, the "ECI: Database
  Deposition" Structure-menu `add_command` block (shortcut "E") + its
  trailing `menu.add_separator()` (the preceding separator kept — menu now
  reads `... H Bond | separator | Secondary Structure Chart, Ramachandran`),
  the `menu_items[StructureMenu]` "ECI: Database Deposition" entry (S6
  precedent; "DANGLE": Predict Dihedrals" string stays — S5 residual), and
  the `startECI` method.
- `pyproject.toml`: dropped the `ccpnmr-eci` console script.
- `gui_boot_test.py`: dropped the APPS entry `eci` and fixed the docstring
  usage example (`--apps ccpnmr,eci` → `--apps ccpnmr,data-shifter`).
- Residuals (by design, Stage-12 sweep): `bin/eci2.5` launcher (S5
  `bin/dangle*` precedent), `scripts/{linux,macos}_release.sh` `need` sets
  still list `ccpnmr-eci` (same stale bucket as the S2/S5/S9 deferrals),
  doc links `ccpnmr/analysis/doc/source/menu/Structure.rst:23` +
  `.../other/ImportProject.rst:47` (S8 `Structure.rst` precedent),
  `EditExperiment.py:292,331` help link def (S5 precedent),
  `Ccpn2NmrStar.py:224,245` EBI-docs URL docstrings (inert text),
  `survey.md:315` snapshot line.
- Env note (pre-existing, NOT repo state): the anaconda env
  (`/home/logan/software/anaconda3`) contained a **stale installed `cci`
  copy** in site-packages (legacy egg 2.5.2, dated May 2026 — still had
  `cci/`, lacked the relocated `ReadPdb.py`; also owned `ccp`, `memops`,
  `cambridge`, `cijermen`, ... + site-packages-root `data/` `doc/`
  `model/` `license/` trees). All three gates resolve to the source tree
  (import_smoke `ROOT` insert, gui_boot `PYTHONPATH=source-tree`,
  `tests/conftest.py` path insert), but a bare `python -c "import
  ccpnmr..."` from a cwd silently imported the stale copy. **RESOLVED
  2026-08-24:** `pip uninstall ccpnmr` removed exactly the egg's registry
  footprint — `import ccpnmr` now fails cleanly, no stray top-level
  packages remain (`pdbeccdutils`, an unrelated package, intact); all
  three gates re-verified green afterwards (TOTAL 1279 / OK 1267 /
  FAILED 2 / BY-DESIGN 10; gui 4/4; pytest 45 passed / 4 skipped).
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1279** (1290−11, exactly the 11
    removed eci `.py` modules; the relocated ReadPdb is still walked as
    `ccpnmr.format.converters.ReadPdb`), OK **1267** (1278−11),
    FAILED **2 unchanged** (2× `cherrypy` in KEPT `ccpnmr/format/webServer/*`),
    BY-DESIGN **10 unchanged**.
  - `gui_boot_test.py` **4/4** (5/5 − the removed `cci` app, exactly as
    predicted).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - Direct check: `from ccpnmr.format.converters.ReadPdb import ReadPdb`
    OK from the source tree; zero `ccpnmr.eci` refs anywhere in
    AnalysisPopup.py; repo-wide `ccpnmr.eci` grep clean outside
    `dist/`/`build/`/Stage-12 buckets.

**Stage 11 — Structure Viewer + Make H Bond Restraints (popups, menu, kept callers) — ✅ 2026-08-24**
Recon (verified pre-edit): popups `ccpnmr/analysis/popups/ViewStructure.py`
(959 lines) + `MakeHbondRestraints.py` (2357 lines) = 3316 lines total.
External touch points: `AnalysisPopup.py` (2 imports + popupActions ×2 +
2 menu blocks + menu_items ×2 + `viewStructure`/`makeHbonds` methods) and
**7 kept caller files**. Locked decision 2 named 3 of them (PeakTableFrame,
WindowFrame, CalcShiftDifference); recon found 4 more (EditAssignment,
EditStructures, LinkNoeResonances, BrowseConstraints) — each opens the
viewer via `guiParent.viewStructure(...)` + `popups["view_structure"]` and
draws connections/highlights; leaving them would crash kept GUIs with
AttributeError/KeyError (gates are import-time only, so they'd ship broken)
→ removed in the same mode as the named 3.
- `PeakTableFrame.py`: "Show On Structure" button (tip/texts/commands) +
  `showStructConnections` + `showAllStructConnections` + the
  `updateButtons` `buttons[12].disable()` gate (index verified = that
  button: bottomButtons1 has 9, so index 12 = 4th of bottomButtons2).
  KEPT: `structPulldown`/`self.structure`/`updateStructures` — they feed
  the peak-table distance column too (`getDistanceDimensions`).
- `WindowFrame.py`: "Structure connections" context-menu entry (shortcut
  "c") + `showStructConnections` method.
- `CalcShiftDifference.py`: "Show On Structure" button entry +
  `showStructure` method (`displayAtomParamsList`).
- `EditAssignment.py`: "Show On Structure" button entry +
  `showStructConnections` + its sole helper `haveCommonPeakContrib`;
  dropped `getThroughSpaceDataDims` from the shared import (kept
  `getOnebondExpDimRefs`, still used); removed `structButton` ×3 refs and
  re-indexed `peaksButton`/`mergeButton`/`infoButton`/`predictButton`
  (button shifted left by one). KEPT: `structurePulldown`/`structure`
  machinery (used by non-viewer resonance handling at L2019/2035/2203).
- `EditStructures.py`: residue-tab "View Residue" + "Display Params" buttons
  + `viewResidue`/`displayStrucParams`/`viewStruct` methods + "Viewer"
  side-tab ButtonList + its `updateButtons` enable block;
  `strucParamPulldown` tip reworded to deletion-only (pulldown +
  `deleteStrucParams` KEPT). `structButtons` row verified viewer-unrelated.
- `LinkNoeResonances.py`: "Structure Display:" pulldown + "Focus
  Structure:" checkbox + `updateAssignments` `showConn` open/draw blocks +
  `coordAtom` compute (incl. the per-assignment atom-count loop — verified
  its only consumer) — all viewer-specific. KEPT: `structurePulldown`/
  `self.structure` (atomic distances feed the table) + unrelated
  `focusSelectC/H` in the spin-systems popup.
- `BrowseConstraints.py`: "Show Selected On Structure" button entry +
  `showStructConnections` method. KEPT: `self.structure` (violations) +
  L274 help-link def (below).
KEPT (hazard 3 re-verified live): `ccp/gui/ViewStructureFrame.py` —
`ViewIsotopomerFrame.py:43` still imports `symbolMultiplier` from it (also
`ViewChemCompVarFrame` subclass). No import_smoke allowlist, pyproject,
gui_boot, tests, `bin/`, or `scripts/` references.
- Deleted: `popups/ViewStructure.py` + `popups/MakeHbondRestraints.py`
  (3316 lines) + `__pycache__` residue.
- `AnalysisPopup.py`: removed 2 imports, popupActions
  `view_structure` + `make_hydrogen_bonds`, "Structure Viewer" (V) +
  "Make H Bond Restraints" (H) menu blocks, menu_items ×2, `viewStructure`
  + `makeHbonds` methods. Structure menu now reads: … Structures | Make
  Distance Restraints | Secondary Structure Chart, Ramachandran.
- Residuals (by design, Stage-12 sweep): doc menu links
  `Structure.rst:14,16` (ViewStructurePopup/MakeHbondRestraintsPopup —
  the rst target files don't exist, links already dangling),
  `Changes.html:1279` changelog line, `BrowseConstraints.py:274` help-link
  def (S5 precedent).
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1277** (1279−2, exactly the 2
    removed popup modules), OK **1265** (1267−2), FAILED **2 unchanged**
    (2× `cherrypy` in KEPT webServer), BY-DESIGN **10 unchanged**.
  - `gui_boot_test.py` **4/4** (unchanged).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - Direct checks: repo-wide `view_structure`/`make_hydrogen_bonds`/
    popup-class-name grep clean (only the 3 Stage-12 residuals above);
    `py_compile` clean on all 9 edited files; diff = 12+/391− across the
    8 edited files + 3316 lines in the 2 deleted popups.

**Stage 12 — Cross-cutting sweep (LAST stage) — ✅ 2026-08-24**

Final verification of 1—11; removed every residual reference to the 14
legacy tools + their orphaned third-party deps. PLAN CLOSED (all 12 stages
done, gates green, pushed).
- `bin/`: deleted 12 dead launchers — `dangle{,2,2.5}`, `eci{,2,2.5}`,
  `extendNmr{,2,2.5}`, `depositionFileImporter{,2,2.5}` (all pointed at
  deleted `python/<pkg>/…Gui.py` targets). KEPT and re-verified live:
  `pipe2azara*` + `xeasy2azara*` (targets `ccp/format/spectra/params/
  {NmrPipeData,XeasyData}.py` still exist) and the 6 live app launchers.
- `pyproject.toml` `optional` extras — recon by whole-tree grep (0 importers
  of a package = orphan, since the smoke gate requires kept importers to be
  importable in a clean install-venv):
  - REMOVED (zero importers in the live tree): `scipy`, `sqlalchemy`,
    `psycopg2-binary`; also dropped the stale `# used by cing / web-server /
    advanced I/O` comment and the now-pointless pycurl comment line (pycurl
    too has ZERO importers).
  - KEPT with importers (each verified): `matplotlib`
    (pdbe/software/violationStatistics.py + cambridge/bayes/kmeans.py, both
    by-design-kept), `cherrypy` + `mako` (kept webServer), `decorator` —
    **deviation from the stage checklist**: `ccpnmr/v2io/TestNefIo.py`
    (kept reference package, by-design) does a TOP-LEVEL `import decorator`;
    dropping it would break the clean-venv release gate (FAILED must be 0)
    and the v2io test, so it stays.
- Release scripts (`linux_release.sh`, `macos_release.sh`,
  `make_standalone_release.sh`): removed the dead superpose C-ext build step
  (`ccpnmr2.5/python/cing/Libs/cython` — deleted in S9), the Cython
  prerequisite install, and the `scipy sqlalchemy psycopg2-binary pycurl`
  optional-stack pip line (now `matplotlib cherrypy decorator mako`);
  `need` entry-point sets 8→4 (`ccpnmr-eci/-dangle/-deposition/-extend-nmr`
  were removed as `[project.scripts]` in S5/S10); step labels [n/4]→[n/3];
  "8 console entry points" → 4. `bash -n` clean on all three; zero
  superpose/Libs/cing/scipy/psycopg residue. `copy_cext.sh`/`publish.sh`
  already clean (S7).
- Docs/help-links (all listed residuals): `Structure.rst` toctree 8 lines
  (Viewer, MakeHbond, DANGLE, ARIA, HADDOCK, PyRPF, CING, ECI — the
  `source/popups/` dir doesn't even exist, so all were dangling);
  `ImportProject.rst` 2× "viewed via the ECI_ option" phrases + the
  `.. _ECI:` link target; `Changes.html` external-programs `<li>` (PALES/
  MODULE2/MECCANO) + the ViewStructurePopup `<LI>` block;
  `EditExperiment.py` "; initially using the str(CcpNmr ECI)_" + its link
  def; `BrowseConstraints.py` MakeHbond/DANGLE/3J-Coupling help text +
  3 link defs (3J popup died S1 — same category); `Ccpn2NmrStar.py` 2×
  "use ECI (http://…eci.html)" error strings → neutral wording;
  `survey.md` `ccpnmr-eci` sample line; `NmrCalc.py` stale "Used in
  paris/aria, nijmegen/cing, grenoble/BlackledgeModule" docstring →
  current importers (memops/api/Implementation, ccp/api, cambridge/isd,
  EditCalculation). `NmrCijmegen.py:4` residual from the checklist: file no
  longer exists (removed during an earlier stage) — nothing to do.
- RECON FINDINGS beyond the checklist (same category, swept):
  - `README.md`: Phase-4-era doc — refreshed the "What works" gate table to
    current numbers (1277/1265/2/10, 45/4, 4/4), console-command table 8→4,
    optional-stack install line (dropped scipy/sqlalchemy/psycopg2/pycurl),
    ALL GSL/Meccano build instructions (gone since S7), `CCP_GSL_PREFIX`,
    "C/Cython" → "C" (superpose gone S9), wrong clone URL
    (`21tesla/analysis2.5py3` → `21tesla/analysis2.5py3lite` per locked
    decision 1), layout list (dropped `cing`, "30 C extensions" count),
    cing scope note → 14-tools-removed note pointing at this plan,
    by-design count 83→10.
  - `docs/PUBLISHING.md`: optional-extras comment, 8→4 console list, GSL
    paragraph, `-> 8`/`8/8 booted` gate comments, GSL/cython recipe line.
  - `INSTALL.md`: 4 launcher bullets for removed apps (ECI/DANGLE/
    deposition/EXTEND-NMR) — data-shifter/format-converter/update kept.
  - `recipe/meta.yaml` + `recipe/README.md`: GSL + cython host deps,
    Meccano-optional prose.
  - Historical records deliberately UNTOUCHED: `Changes.html` other legacy
    entries (changelog), `MACOS_3.13_BUILD_FIXES.md` (build log),
    `_phase*_checkpoints.md`, `.aider.chat.history.md`,
    `SIMPLIFICATION_PLAN.md`'s own earlier log entries.
- Gates (all green, floors held):
  - `import_smoke.py` exit 0 — TOTAL **1277**, OK **1265**, FAILED **2**
    (2× `cherrypy`, unchanged), BY-DESIGN **10** (unchanged).
  - `xvfb-run -a python gui_boot_test.py` **4/4** (unchanged).
  - `python -m pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped**
    (unchanged).
  - Whole-tree grep for the 4 removed entry-point names + `entrycompletion`
    + `meccano`/`blackledge`: clean outside the historical records above.

**PLAN COMPLETE (2026-08-24): 12/12 stages done and pushed to
`21tesla/analysis2.5py3lite:main`.**

---

# Menu Removal Plan — Stages 13-16 (added 2026-08-24)

Status: **COMPLETE (2026-08-24)** — all 4 menu-removal stages (13-16) done;
16/16 across both plans (12/12 Simplification + 4/4 Menu Removal), all pushed
Same repo, same checkpoint policy: ONE commit per stage (code + this log
update in the same commit) + push to `main`. Python: anaconda `python`
3.13.5; `xvfb-run` available.

## Goal

Remove 4 menus/methods from the main (Project) window
(`ccpnmr2.5/python/ccpnmr/analysis/AnalysisPopup.py`) — the CCPN server that
served program updates is dead, macros are out of scope for this py3-lite
build, and Prodecomp/CLOUDS are no longer offered — plus the orphaned code
that becomes unreachable (stages 1-12 pattern).

| # | Menu item | Builder (AnalysisPopup.py) | Command method(s) | Orphan code removed with it |
|---|---|---|---|---|
| 1 | Project ▶ **Updates** | `setProjectMenu` | `updateAnalysis` | `ccpnmr/update/` (6 files) + `ccpnmr-update` script + gui_boot_test NON_GUI entry |
| 2 | **Macro** (top-level) | `setMacroMenu` (+2 notify hooks) | `setMacroMenu`, `reloadMenuMacros`, `runMacro`, `editMacros` | — (see locked decision 1: core engine KEPT) |
| 3 | Other ▶ **Prodecomp** + **CLOUDS** | `setOtherMenu` | `startProdecomp`; `setupClouds`/`setupBacus`/`setupMidge`/`setupHcloudsMd`/`setupFilterClouds`/`setupCloudThreading`/`setupCloudHomologue` | `gothenburg/prodecomp/` (7 files) + 14 `ccpnmr/clouds/` modules |
| 4 | Project ▶ **Help** (Version/About/Help) | `setProjectMenu` | `showVersion`, `showAbout`, `showHelp` | — |

## Locked decisions (2026-08-24, with user)

1. **Macro = menu layer ONLY.** Remove the top-level Macro menu + its four
   methods + the 2 `notify(self.setMacroMenu, "ccpnmr.AnalysisProfile.Macro", …)`
   hooks + the dangling "In Menu" toggle in the EditProfiles Macros tab.
   **KEEP:** the EditProfiles "Macros" tab itself (and its "In Mouse Menu"
   column), the window right-click "Macros" submenu (`WindowFrame.py`),
   the `ccpnmr/analysis/macros/` engine (ArgumentServer, …), the generated
   `Macro` API/model class. `EditProfiles.py` keeps its
   `from ccpnmr.analysis.core.Util import … reloadMacro, runMacro` import
   (still used by `reloadSelectedMacro` + the WindowFrame mouse menu).
2. **Remove orphan code too** (stages 1-12 pattern): `ccpnmr/update/`
   fully (all external references: AnalysisPopup L160 import, L666
   commented import, L3164 `dataFile` in the `updateAnalysis` method,
   `gui_boot_test.py` NON_GUI entry, `pyproject.toml` script);
   `gothenburg/prodecomp/` fully (sole external import: AnalysisPopup L129
   `ProdecompPopup` import + `startProdecomp` method); `ccpnmr/clouds/`
   **partially** (see decision 3 — two of its modules are LIVE).
3. **`ccpnmr/clouds/` — KEEP (verified live importers, grep-verified
   2026-08-24):** `ResonanceIdentification.py` (imported by KEPT
   `popups/CalcDistConstraints.py:907` `makeNoeAdcs`), `FilterClouds.py`
   (imported by KEPT `popups/EditResStructures.py:247,261`),
   `PseudoResonances.py` + `CloudBasic.py` (their dependencies),
   `__init__.py`, `_licenseInfo.py`. **REMOVE (14 files):** the 7 GUI
   popups `BacusPopup`, `CloudsPopup`, `HcloudsMdPopup`, `MidgePopup`,
   `FilterCloudsPopup`, `CloudThreaderPopup`, `CloudHomologueAssignPopup`,
   plus the helpers only they reach: `Clouds`, `CloudHomologueAssign`,
   `CloudThreader`, `FileIO`, `HydrogenDynamics`, `NoeMatrix`,
   `NoeRelaxation`. (Re-verify `FileIO`/`NoeRelaxation` consumers with a
   broad grep at implementation time — both map to deleted popups only.)
4. **Help = main window only.** The FormatConverter window keeps its own
   top-level Help menu (Glossary/Menus/Quick Start/Tutorial/About/Version).
5. **Registration stays — user: out of scope.** Do NOT touch
   `ccpnmr/analysis/core/Register.py`, `popups/Register.py`,
   `registerAnalysis` / `checkRegistration`, the "Register" Project-menu
   item — even though `updateRegister` POSTs to dead
   `http://www.ccpn.ac.uk/cgi-bin/register/register`.

## Menu-index bookkeeping (Tkinter `entryconfig` — separators COUNT)

Project menu construction order (current):
`0 New, 1 Open Project, 2 Open Spectra, 3 Load Nef, 4 Save, 5 Save As,
6 Import, 7 Close, 8 Quit, 9 sep, 10 Summary, 11 Preferences, 12 Register,
13 Validate, 14 Backup, 15 Archive, 16 Updates, 17 sep, 18 Help`.
Current `fixedActiveMenus[(ProjectMenu,…)] = (0,1,2,3,8,16,18)` (items that
stay enabled with no project open).
- After **Stage 13** (Updates removed; its separator stays and now separates
  Archive from Help): items 17/18 shift down by one →
  `fixedActiveMenus = (0,1,2,3,8,17)`.
- After **Stage 16** (Help + that separator removed):
  `fixedActiveMenus = (0,1,2,3,8)`.
`menu_items[ProjectMenu]` (17 entries incl. `updateText` + `"Help"`):
Stage 13 drops `updateText` → 16 entries; Stage 16 drops `"Help"` → 15.
`Other` menu after **Stage 15**: `0 NMR Calculations, 1 Widget Counter,
2 Format Converter` → `fixedActiveMenus[(OtherMenu, 2)] = True` (Format
Converter) stays correct as-is.

## Baseline (from Stage 12 close, 2026-08-24)

| Gate | Result |
|---|---|
| `python import_smoke.py` | exit 0 — TOTAL **1277**, OK 1265, FAILED **2** (2× `cherrypy`), BY-DESIGN **10** |
| `xvfb-run -a python gui_boot_test.py` | **4/4** (ccpnmr, data-shifter, format-converter + NON_GUI `ccpnmr-update`) |
| `python -m pytest ccpnmr2.5/python/tests/` | **45 passed, 4 skipped** |

Per-stage rules: import_smoke exit 0, no NEW unexpected failures, TOTAL
drops by exactly the count of removed `.py` modules when packages go;
gui_boot_test green (−1 entry in Stage 13 when the `ccpnmr-update`
NON_GUI check is deleted); pytest floors held; `ruff check` on edited files
(F401 is project-ignored but dead imports still get deleted).

## Stages & status

| # | Scope | Status |
|---|---|---|
| 13 | Project ▶ Updates + delete `ccpnmr/update/` + `ccpnmr-update` script + boot-test entry | ✅ 2026-08-24 |
| 14 | Macro menu (menu layer only) | ✅ 2026-08-24 |
| 15 | Other ▶ Prodecomp + CLOUDS + delete orphan `gothenburg/prodecomp/` (7) + 14 `ccpnmr/clouds/` modules | ✅ 2026-08-24 |
| 16 | Project ▶ Help (Version/About/Help) | ✅ 2026-08-24 |

## Stage checklist detail

**Stage 13 — Updates menu + `ccpnmr/update/` removal**
- `AnalysisPopup.py`
  - `setProjectMenu`: drop the docstring-commented `UpdateAgent` import +
    the `numUpdates` / `updateText` block (L664-673); remove the
    `menu.add_command(label=updateText, shortcut="U", …,
    command=self.updateAnalysis, …)` block and keep the following
    separator (it now separates Archive from Help); remove the
    `updateText` entry from `menu_items[ProjectMenu]`.
  - `fixedActiveMenus[(ProjectMenu, …)]` → `(0,1,2,3,8,17)`.
  - Delete the `updateAnalysis` method; delete the
    `from ccpnmr.update.UpdatePopup import UpdatePopup` import (L160).
  - `LOCAL_HELP_DOC_DIR` / `getTopDirectory` imports stay until Stage 16
    (still used by `showAbout`/`showHelp`).
- Delete `ccpnmr2.5/python/ccpnmr/update/` (6 files: `__init__`,
  `_licenseInfo`, `UpdateAdministratorPopup`, `UpdateAgent`, `UpdateAuto`,
  `UpdatePopup`) + `__pycache__` residue.
- `pyproject.toml`: drop `ccpnmr-update = "ccpnmr.update.UpdateAuto:main"`.
- `gui_boot_test.py`: drop the NON_GUI `update` entry + the docstring
  paragraph about UpdateAuto.
- Gates expected: TOTAL 1277→**1271** (−6), FAILED 2, BY-DESIGN 10,
  gui_boot_test 3/3 (GUI apps) with no update line, pytest unchanged.
- Commit + push; stage log; mark ✅.

**Stage 14 — Macro menu (menu layer only)**
- `AnalysisPopup.py`: delete the `MacroMenu = "Macro"` constant; the two
  `self.setMacroMenu()` calls (`__init__` L394, `initProject` L1947); the
  two `notify(self.setMacroMenu, "ccpnmr.AnalysisProfile.Macro", …)`
  lines in `curatePopupNotifiers` (L650-651); the `setMacroMenu`,
  `reloadMenuMacros`, `runMacro`, `editMacros` methods.
  If the `Command` import (L49 `from ccp.general.Command import Command`)
  is used nowhere else, drop it (verify first — `runMacro` was its only
  known user).
- `EditProfiles.py`: delete `toggleMacroInMenu` (it calls
  `self.parent.setMacroMenu()` — would AttributeError) + its "In Menu"
  column wiring in the macro table (de-index the remaining columns; the
  "In Mouse Menu" column + `toggleMacroInMouseMenu` STAY).
- Grep-verify zero remaining `setMacroMenu|editMacros|MacroMenu|
  reloadMenuMacros|AnalysisProfile.Macro.*delete/setName` notify hooks.
- Gates: unchanged floors; TOTAL still 1271 (no `.py` deleted);
  gui_boot_test 3/3. Commit + push; stage log; mark ✅.

**Stage 15 — Prodecomp + CLOUDS menus + orphan removal**
- `AnalysisPopup.py`
  - `setOtherMenu`: remove the whole `cloudsMenu` construction (6
    commands), the Prodecomp `menu.add_command` block, the
    `menu.add_cascade(label="CLOUDS", …)` line; `menu_items[OtherMenu]`
    → `["NMR Calculations", "Widget Counter", "Format Converter"]` (the
    `fixedActiveMenus[(OtherMenu, 2)] = True` line stays valid).
  - Delete methods `startProdecomp`, `setupClouds`, `setupBacus`,
    `setupMidge`, `setupHcloudsMd`, `setupFilterClouds`,
    `setupCloudThreading`, `setupCloudHomologue`.
  - Remove the 7 clouds entries from the `self.popups` dict
    (`setup_clouds`, `setup_bacus`, `setup_midge`, `setup_hcloudsmd`,
    `setup_filter_clouds`, `setup_cloud_threader`, `setup_cloud_homologue`).
  - Delete `from gothenburg.prodecomp.ProdecompFrame import ProdecompPopup`
    (L129) + the `self.iconClouds = self.icons["weather-overcast"]` line.
    (Tip: the "NMR Calculations" tipText still names the removed CING/ARIA
    — reword to neutral "dispatch calculation jobs to external programs"
    as part of this menu touch.)
- Delete `ccpnmr2.5/python/gothenburg/prodecomp/` (7 files) and the 14
  `ccpnmr/clouds/` modules per locked decision 3 (KEEP the 6). Re-verify
  with a broad grep that the 7 KEEPED modules have no import of a deleted
  one (`ResonanceIdentification` → CloudBasic/PseudoResonances only ✓).
- Gates expected: TOTAL 1271→**1247** (−7−14 = −24), FAILED 2,
  BY-DESIGN 10, gui_boot_test 3/3, pytest unchanged; grep
  `ccpnmr\.clouds` outside the package → exactly
  `CalcDistConstraints.py` + `EditResStructures.py`; grep
  `gothenburg.prodecomp|ProdecompPopup|startProdecomp|setupBacus|setupMidge|
  setupHcloudsMd|setupFilterClouds|setupCloudThreading|setupCloudHomologue`
  → clean. Commit + push; stage log; mark ✅.

**Stage 16 — Help menu (Project)**
- `AnalysisPopup.py`: delete the `helpMenu` construction (Version/About/
  Help commands), the `menu.add_cascade(label="Help", …)` line + the
  separator immediately before it, the `"Help"` `menu_items[ProjectMenu]`
  entry; `fixedActiveMenus[(ProjectMenu, …)]` → `(0,1,2,3,8)`; delete
  `showVersion`, `showAbout`, `showHelp` methods; delete
  `self.iconHelp = self.icons["help-browser"]` (only other use is the
  deleted cascade); if now unused, drop `LOCAL_HELP_DOC_DIR` from the
  Analysis import (L53) and `getTopDirectory` from the
  `memops.universal.Io` import (L146) — verify no other users in-file.
- Gates: unchanged floors; TOTAL still 1247; gui_boot_test 3/3.
- Commit + push; stage log; mark the Menu Removal Plan **COMPLETE**
  (16/16 across both plans) in the status line + this stage's log entry.

## Rollback

Each stage is exactly one commit → `git revert <sha>` restores it cleanly.

## Stage log (Menu Removal)

**Stage 13 — Project ▶ Updates menu + `ccpnmr/update/` removal — ✅ 2026-08-24**
Recon (verified pre-edit): `ccpnmr2.5/python/ccpnmr/update/` = 7 tracked
files: the 6 `.py` (`__init__` (bare `pass`), `_licenseInfo`,
`UpdateAdministratorPopup`, `UpdateAgent`, `UpdateAuto`, `UpdatePopup`) plus
`uploadFile` — a py-2 CGI script for the dead update server
(`except Exception, e:` syntax; not importable, zero inbound references) that
the plan's "6 files" count missed (+ untracked `__pycache__` residue).
External refs: only `AnalysisPopup.py` (the L160 `UpdatePopup` import, the
`setProjectMenu` docstring-comment import + `numUpdates`/`updateText` block,
the Updates `add_command` block, the `menu_items` entry, and the
`updateAnalysis` method), the `pyproject.toml` `ccpnmr-update` script, and
the `gui_boot_test.py` NON_GUI entry. Kept: `iconRefresh` (still used at
L1775), `LOCAL_HELP_DOC_DIR`/`getTopDirectory` imports (used by
`showAbout`/`showHelp` until Stage 16), `Copyright` import (`__init__`
L322/324).
- `AnalysisPopup.py`: dropped the `UpdatePopup` import (L160); dropped the
  `setProjectMenu` docstring (UpdateAgent comment) + the
  `numUpdates`/`updateText` block; dropped the `Updates`
  `menu.add_command` block (the following separator KEPT — now separates
  Archive from Help); dropped the `updateText` entry from
  `menu_items[ProjectMenu]` (17 → 16); `fixedActiveMenus[(ProjectMenu, …)]`
  → `(0,1,2,3,8,17)` (old 17/18 shift down one — separator counts as
  index, so Help is now index 17); deleted the `updateAnalysis` method.
- Deleted `ccpnmr2.5/python/ccpnmr/update/` entirely (7 tracked files +
  `__pycache__`).
- `pyproject.toml`: dropped `ccpnmr-update = "ccpnmr.update.UpdateAuto:main"`.
- `gui_boot_test.py`: emptied the NON_GUI `update` entry
  (`NON_GUI = []`; generic loop + `total` computation unchanged, now a no-op)
  + removed the UpdateAuto docstring paragraph.
- Grep-verified: repo-wide `ccpnmr.update|ccpnmr/update|ccpnmr-update` clean
  outside the gitignored `dist/` build snapshot, `.aider.chat.history.md`
  chat history, and this plan doc; the only remaining
  `updateAnalysis*` matches are the UNRELATED `updateAnalysisSpectra*` /
  `updateAnalysisPeakList` methods in `EditSpectrum.py` / `EditPeakLists.py`
  (different symbols).
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1271** (1277−6, exactly the 6 removed
    `.py` modules), OK **1259** (1265−6), FAILED **2 unchanged** (2×
    `cherrypy`), BY-DESIGN **10 unchanged**.
  - `gui_boot_test.py` **3/3** (4/4 − the removed `update` entry, exactly as
    planned; no update line in the output).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `ruff check` on both edited `.py` files: zero NEW violations vs HEAD —
    AnalysisPopup UP031 19→17 (the two removed %-format strings live in the
    deleted code), all other rules unchanged (F841 11, E722 10, E731 2,
    W293 1, F811 1, E721 1); `gui_boot_test.py` identical (UP031 8).

**Stage 14 — Macro menu (menu layer only) — ✅ 2026-08-24**
Recon (verified pre-edit): `AnalysisPopup.py` targets at (line numbers post
Stage 13): `MacroMenu = "Macro"` L167; `self.setMacroMenu()` in `__init__`
L392 and `initProject` L1925; two `notify(self.setMacroMenu,
"ccpnmr.AnalysisProfile.Macro", "delete"/"setName")` hooks in
`curatePopupNotifiers` L648-649; contiguous block `setMacroMenu` L1731 /
`reloadMenuMacros` L1774 / `runMacro` L1780 (immediately before
`setOtherMenu`); `editMacros` L2919 (one-liner `popupEditProfiles(tab=1)`).
Checklist's "verify first" on the `Command` import (L49): case-sensitive
` Command(` shows `runMacro` (L1784) is the ONLY user (every other `command`
hit is the `add_command`/`command=` keyword) → import dropped. `iconRefresh`
(L305): sole usage was `setMacroMenu`'s "Reload Menu Macros" item (L1753) —
the Stage 13 log's "still used" note pointed exactly at this macro code, and
grep confirms no other user in the tree → assignment dropped as orphan of
this stage (the other macro-menu icon `iconTable` has 20+ live uses — kept).
`sys` (L44) + `Util` imports both have other users (`sys.exit` L3036;
`Util.setTopObjectAnalysisSaveTime`, `Util.getFormatConverterThreading`) →
kept. Cross-file: ONLY `EditProfiles.py:1173`
(`self.parent.setMacroMenu()` inside `toggleMacroInMenu`) — the AttributeError
hazard the checklist names; `WindowFrame.py` macro paths use
`isInMouseMenu` + `Util.runMacro`/`Util.reloadMacro` (KEPT, untouched);
`OpenMacro.py` clean; no `fixedActiveMenus`/`menu_items` entries for Macro
(`setMenuState` iterates `self.menus` generically — the Macro key simply
stops existing). EditProfiles macro table had 9 columns with col 3
"In main\nmenu?" (get=`toggleMacroInMenu`, set=None); `updateMacros` built 9
tuples, `blankColors = [None] * 9`, colored col 3 (in-menu Yes) + col 4
(in-mouse Yes). KEPT per locked decision 1: the Macros tab itself, the
"In mouse\nmenu?" column + `toggleMacroInMouseMenu`, the
`from ccpnmr.analysis.core.Util import reloadMacro, runMacro` import (live in
`reloadSelectedMacro`/`runSelectedMacro` + WindowFrame), the
`ccpnmr/analysis/macros/` engine, generated `Macro` API/model incl. the
`isInMenu` attribute (now unsurfaced from UI but still valid model state),
and the generic `administerNotifiers` Macro subscription (L681 — serves
`updateMacrosAfter`).
- `AnalysisPopup.py` (−70 lines): removed `from ccp.general.Command import
  Command`; `MacroMenu = "Macro"`; `iconRefresh` assignment; the `__init__`
  `self.setMacroMenu()` call; the two `AnalysisProfile.Macro` notify hooks in
  `curatePopupNotifiers`; the `setMacroMenu` + `reloadMenuMacros` + `runMacro`
  method block; the `initProject` `self.setMacroMenu()` call; the `editMacros`
  method.
- `EditProfiles.py` (−19/+12): removed `toggleMacroInMenu` + the "In main
  \nmenu?" column from `headingList`, its tipText, and the three per-column
  lists (`editWidgets`/`editGetCallbacks`/`editSetCallbacks` — each 9→8);
  `updateMacros` now builds 8-element rows, `blankColors = [None] * 8`,
  dropped the `macro.isInMenu` text cell + `colors[3]` highlight, and the
  mouse-menu highlight moved col 4 → col 3.
- Grep-verified clean: `\b(setMacroMenu|reloadMenuMacros|editMacros|
  MacroMenu|toggleMacroInMenu)\b` → zero matches in `**/python/**/*.py`;
  `isInMenu` no longer appears in EditProfiles.py (remains only in generated
  `ccpnmr/api` / `ccpnmr/xml` / `model` — KEPT); `python -m py_compile` OK on
  both files.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1271** (unchanged — no `.py` removed,
    exactly as planned), FAILED **2** unchanged (2× `cherrypy`), BY-DESIGN
    **10** unchanged.
  - `gui_boot_test.py` **3/3** (the CCPN main window boots through the
    edited `__init__` menu-construction path).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `uvx ruff check` worktree vs HEAD: AnalysisPopup 43 = 43 (E722 10,
    F841 11, UP031 17, E731 2, E721 1, F811 1, W293 1 — identical mix),
    EditProfiles 10 = 10 (UP031 6, F841 3, E722 1 — identical mix) — zero
    NEW violations.

**Stage 15 — Other ▶ Prodecomp + CLOUDS menus + orphan removal — ✅ 2026-08-24**
Recon (verified pre-edit): `AnalysisPopup.py` targets at (line numbers post
Stage 14): `ProdecompPopup` import L129 under a one-off `# NB new` marker
(unique in the file — dropped with the import); 7 `self.popups` entries
L277-283 (NO pre-existing `"prodecomp"` entry — `startProdecomp` creates
it at call time, L2368); `self.iconClouds` assignment L303 (sole other use
is the CLOUDS cascade itself); `setOtherMenu` L1725-1809 — the whole
`cloudsMenu` 6-command build + a commented-out `entryconfig` stub, the
Prodecomp `add_command` block, the `add_cascade(CLOUDS)` line, and
`menu_items[OtherMenu]` 5→3 items; `startProdecomp` L2366-2368; 7 `setup*`
methods L2881-2921 (each: local popup import + `openPopup`).
`gothenburg/prodecomp/` = 7 tracked `.py` (CcpnProdecomp, generateInterval,
PeaksToInterval, ProdecompFrame, Projection, prodecomp, `__init__`) +
`__pycache__` residue; parent `gothenburg/` KEEPS `mdd/` + `__init__.py`,
so the two `gothenburg` package entries in `pyproject.toml` (L86
`gothenburg*` data dir, L149 first-party list) stay. `ccpnmr/clouds/` =
20 tracked `.py`; the 6 KEPT modules' cross-imports verified:
`ResonanceIdentification` → `CloudBasic` + `PseudoResonances` (both KEPT),
`FilterClouds` → `CloudBasic` (KEPT); the only relative imports reaching
deleted modules (`from .HydrogenDynamics`, `from .NoeRelaxation`) sit in
`Clouds.py` itself (removed); repo-wide `openPopup("setup_*")` +
`self.popups["prodecomp"]` writers are ONLY the 9 methods removed; the only
cross-file importers of the 14 removed modules are the 7 `setup*` methods
(`EditResStructures.py` → `FilterClouds` KEPT, `CalcDistConstraints.py` →
`ResonanceIdentification` KEPT). REVIEWED + DEFERRED (consistent with
Stage 14, which left `menu/Macro.rst` in place): `doc/source/menu/Clouds.rst`
+ the `menu/Other.rst` toctree line now dangle (the doc tree's `popups/`
target dir doesn't exist — most menu doc links are already dead);
`EditSpectrum.py:238`, `Nmr.py:31183`, `BrukerParams.py:189/194` mention
"Prodecomp/PRODECOMP" as the third-party program name in docstrings/
comments (not the removed package; no imports); `setup.py:29,156` +
`linkSharedObjs`/`copySharedObjsMac` reference `c/ccpnmr/clouds/` — the C
extension build dir, out of scope.
- `AnalysisPopup.py` (1+/104−): dropped the `gothenburg.prodecomp` import
  + `# NB new` marker (the blank line ruff's single first-party isort group
  left behind was folded back — otherwise one NEW I001); 7 `self.popups`
  entries; the `iconClouds` assignment; in `setOtherMenu`: the whole
  `cloudsMenu` build + commented stub, the Prodecomp `add_command` block,
  the `CLOUDS` cascade; `menu_items[OtherMenu]` → `["NMR Calculations",
  "Widget Counter", "Format Converter"]` (`fixedActiveMenus[(OtherMenu, 2)]`
  stays valid); `startProdecomp`; the 7 `setup*` methods; and the
  "NMR Calculations" tipText — `"sent o external programs like, CING or
  ARIA"` (dead CING/ARIA + typo) → `"Curate and manage calculation jobs
  dispatched to external programs"` per the checklist.
- Deleted `gothenburg/prodecomp/` (7 tracked + `__pycache__`) and the 14
  `ccpnmr/clouds/*.py` (BacusPopup, Clouds, CloudsPopup,
  CloudHomologueAssign, CloudHomologueAssignPopup, CloudThreader,
  CloudThreaderPopup, FileIO, FilterCloudsPopup, HcloudsMdPopup,
  HydrogenDynamics, MidgePopup, NoeMatrix, NoeRelaxation) + clouds
  `__pycache__` + stale `gothenburg/__pycache__`.
- Grep-verified clean: `gothenburg\.prodecomp|ProdecompPopup|
  startProdecomp|iconClouds|cloudsMenu|setup_clouds|setup_bacus|setup_midge|
  setup_hcloudsmd|setup_filter_clouds|setup_cloud_threader|
  setup_cloud_homologue` → zero matches in code (only this plan doc);
  `from ccpnmr.clouds` outside the package → exactly `CalcDistConstraints.py`
  + `EditResStructures.py` (plus one pre-existing COMMENT line for the C
  class `AtomCoordList` — unrelated); repo-wide `NoeMatrix` hits are the
  DIFFERENT kept module `ccpnmr.analysis.frames.NoeMatrix` (and the
  `ViewNoeMatrix` popup) — not the removed `ccpnmr.clouds.NoeMatrix`.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1250** = 1271−**21** (7 prodecomp +
    14 clouds — exactly the 21 removed `.py` files; the checklist's "1247
    (−7−14 = −24)" was an arithmetic slip, 7+14=21), OK **1238**
    (1259−21), FAILED **2 unchanged** (2× `cherrypy`), BY-DESIGN **10
    unchanged** (5 ENV + 5 EXTERNAL).
  - `gui_boot_test.py` **3/3** (the CCPN main window boots through the
    edited `setOtherMenu` path).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `uvx ruff check` worktree vs HEAD: AnalysisPopup 43 = 43 (UP031 17,
    F841 11, E722 10, E731 2, E721 1, F811 1, W293 1 — identical mix; the
    transient I001 from the import removal was fixed by folding the blank
    line) — zero NEW violations. `python -m py_compile` OK.

**Stage 16 — Project ▶ Help menu (Version/About/Help) — ✅ 2026-08-24**
Recon (verified pre-edit, line numbers pre-Stage-16): `AnalysisPopup.py`
targets: L52 `LOCAL_HELP_DOC_DIR` import; L142 `getTopDirectory` (same
import line also carries `getPythonDirectory` + `joinPath`); L290
`self.iconHelp` assignment; L701-705 the `# Help Submenu` + 3-command
`helpMenu` build (followed by a stray lone-`#` leftover line); L824-825
separator + `menu.add_cascade(label="Help", …)`; L848 `"Help"` in
`menu_items[ProjectMenu]` (16 entries); L854 `fixedActiveMenus` loop
`(0, 1, 2, 3, 8, 17)`; L2953-2967 `showVersion`/`showAbout`/`showHelp`.
Orphan checks (grep, in-file + repo-wide): `LOCAL_HELP_DOC_DIR` and
`getTopDirectory` used ONLY by `showAbout`/`showHelp` → both dropped;
`getPythonDirectory` (L173 `GFX_DIR`) + `joinPath` (L2569) still live →
import line keeps them; `showInfo` still used at L2468/2496/2562 → import
KEPT; `self.versionInfo` only READ by `showVersion` here (assigned in
parent `Analysis` — out of scope, untouched); `iconHelp` other uses =
none after the cascade drops; zero repo-wide external callers of
`.showVersion(`/`.showAbout(`/`.showHelp(` (the memops editor's own
`showVersion` and `memops.gui.HelpPopup.showHelp*` are unrelated symbols
in other files — KEPT, as is the FormatConverter window's own Help menu
per locked decision 4).
- `AnalysisPopup.py` (3+/33−, net −30 lines): dropped `LOCAL_HELP_DOC_DIR`
  from the `ccpnmr.analysis.Analysis` import; `getTopDirectory` from the
  `memops.universal.Io` import; the `iconHelp` assignment; the whole
  `helpMenu` build (`# Help Submenu` comment + Version/About/Help commands
  + the stray lone-`#` leftover line between it and `menu = Menu(…)`); the
  `menu.add_separator()` + `menu.add_cascade(label="Help", …)` pair (the
  separator Stage 13 kept to separate Archive from Help — now both
  dropped; the Project menu ends at "Archive"); the `"Help"`
  `menu_items[ProjectMenu]` entry (16 → 15); `fixedActiveMenus` loop
  `(0, 1, 2, 3, 8, 17)` → `(0, 1, 2, 3, 8)`; the `showVersion`,
  `showAbout`, `showHelp` methods. Kept the historical
  `# for ii in (0,1,2,7,15,17):` comment line (pre-existing, unrelated
  text).
- Grep-verified clean in-file: `showVersion|showAbout|showHelp|iconHelp|
  helpMenu|LOCAL_HELP_DOC_DIR|getTopDirectory` → zero matches;
  `python -m py_compile` OK.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1250** (unchanged — no `.py`
    removed, exactly as planned), OK **1238**, FAILED **2 unchanged**
    (2× `cherrypy`), BY-DESIGN **10 unchanged** (5 ENV + 5 EXTERNAL).
  - `gui_boot_test.py` **3/3** (the CCPN main window boots through the
    edited `setProjectMenu` path).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `uvx ruff check --statistics`: AnalysisPopup 43 = 43 (UP031 17, F841
    11, E722 10, E731 2, E721 1, F811 1, W293 1 — identical mix to the
    Stage 13/14/15 baseline; the three dropped methods carried no
    %-format strings, so UP031 is flat) — zero NEW violations.
- **MENU REMOVAL PLAN COMPLETE — 16/16** (Simplification 1-12 + Menu
  Removal 13-16), all stages committed + pushed to
  `21tesla/analysis2.5py3lite:main`.

---

# Assignment Menu Removal Plan — Stages 17-21 (added 2026-08-24)

Status: **COMPLETE (2026-08-24)** — all 5 assignment-menu stages (17-21)
done; 21/21 across all three plans (12/12 Simplification + 4/4 Menu
Removal + 5/5 Assignment Menu), all pushed
Same repo, same checkpoint policy: ONE commit per stage (code + this log
update in the same commit) + push to `main`. Python: anaconda `python`
3.13.5; `xvfb-run` available.

## Goal

Remove 5 items from the Assignment menu (`setAssignMenu` in
`ccpnmr2.5/python/ccpnmr/analysis/AnalysisPopup.py`) plus their popup
modules and the orphan code that becomes unreachable (stages 1-16 pattern).
Per user request (2026-08-24).

| # | Menu item | Command method (AnalysisPopup.py) | Popup/module removed | Orphans removed with it |
|---|---|---|---|---|
| 17 | Assignment ▶ **Initialise Root Resonances** | `initialiseRootSpectra` | `popups/InitRootAssignments.py` (756) | — |
| 18 | Assignment ▶ **Pick & Assign From Roots** | `linkPeakLists` | `popups/LinkPeakLists.py` (1356) | the S17 cross-references inside it |
| 19 | Assignment ▶ **Protein Sequence Assignment** | `linkSeqSpinSystems` | `popups/LinkSeqSpinSystems.py` (3173) | `EditSpinSystem.py:165` dangling help link |
| 20 | Assignment ▶ **Automated Seq. Assignment** (Nexus/MARS/PSIPRED) | `autoBackboneAssign` | `ccpnmr/nexus/` (5 files, 3538 lines) | `wrappers/Mars.py` (507) + `wrappers/Psipred.py` (188) |
| 21 | Assignment ▶ **NOE Contributions** | `linkNoeResonances` | `popups/LinkNoeResonances.py` (1419) | — |

Kept in the Assignment menu: Assignment Panel, Copy Assignments,
Spin System Typing, Assignment Graph, Quality Reports.

## Locked decisions (2026-08-24)

1. **One stage per menu item, menu order 17→21.** `NexusBasic.py` is shared
   by `LinkSeqSpinSystems` (S19) and `AutoBackbonePopup` (S20) — it goes with
   S20, the stage that removes its LAST consumer.
2. **`wrappers/` — remove `Mars.py` + `Psipred.py` only** (S20): `Mars.py`'s
   sole importer is `AutoBackbonePopup` (`runMars`); `Psipred.py`'s sole
   importer is `Mars.py` (lazy `psipredCcpn`). KEPT: `wrappers/CamCoil.py`
   (imported by KEPT `popups/SequenceShiftPredict.py:44`), `wrappers/D2D.py`
   (imported by KEPT `popups/SecStructurePredict.py:44`); `wrappers/Shiftx.py`
   has no external importers but is a separate SHIFTX-prediction wrapper (an
   `argServer` macro entry point), NOT one of the 5 requested items — kept
   (recorded as observation).
3. **`macros/AssignmentMacros.py` KEPT** (Stage-14 precedent: macro surface
   outlives menu removal). Its `initialiseHSQC` /
   `initialiseHNCOorHNCOCA` / `pickAssignSpecFrom*` are standalone
   `argServer` user-macro entry points duplicating (not shared with) the
   removed popups; the file also holds 10+ generic assignment macros. The
   KEPT core `ccpnmr.analysis.core.AssignmentAdvanced` (incl. its
   `pickAssignSpecFromRoot` / `assignSpecNonRootResonances` used by
   `LinkPeakLists`) stays either way.
4. **Generated API attributes `AnalysisProject.linkPeakListsData` /
   `linkSeqSpinSystemsData` KEPT** (`ccpnmr/api` + `ccpnmr/xml` +
   `memops/api` + `ccpnmr2.5/model` — Stage-14 precedent: generated model
   keeps attributes after their UI is removed; plain project data fields).
5. **`HAVE_NUMPY` (AnalysisPopup L145-148) KEPT** — also used by the KEPT
   `peakSeparatorParams` method. KEPT: `showWarning` (many users),
   `confirm_seq_spin_systems` / `type_spin_systems` / `edit_assignment` /
   `quality_reports` popupActions + methods (their menu items stay).
6. **Docs:** each stage drops its own `doc/source/menu/Assignment.rst`
   toctree line (targets are already dangling — the `source/popups/` dir
   doesn't exist; Stage-12 precedent); `doc/Changes.html` changelog history
   UNTOUCHED (Stage-12 precedent); `doc/Readme.txt:66` (`nexus/` layout line)
   dropped in S20.

## Menu-index bookkeeping

Assignment menu construction order (current, separators count):
`0 Assignment Panel, 1 Copy Assignments, 2 Spin System Typing, 3 sep,
4 Initialise Root Resonances, 5 Pick & Assign From Roots,
6 Protein Sequence Assignment, 7 Automated Seq. Assignment,
8 NOE Contributions, 9 sep, 10 Assignment Graph, 11 Quality Reports`.
No `fixedActiveMenus` entry for AssignMenu (grep-verified: only Project /
Experiment / Other have them) → no re-indexing needed. `menu_items[AssignMenu]`
is built via `menuNames.append` and `setMenuState` (L1871+) iterates it
generically (Stage-14 precedent) → only the entries for the 5 items disappear.
After S21: `0 Assignment Panel, 1 Copy Assignments, 2 Spin System Typing,
3 sep, 4 Assignment Graph, 5 Quality Reports` (S21 collapses the pair of
separators left adjacent after the 5 items are gone).

## Baseline (from Stage 16 close, 2026-08-24)

| Gate | Result |
|---|---|
| `python import_smoke.py` | exit 0 — TOTAL **1250**, OK 1238, FAILED **2** (2× `cherrypy`), BY-DESIGN **10** |
| `xvfb-run -a python gui_boot_test.py` | **3/3** (ccpnmr, data-shifter, format-converter) |
| `python -m pytest ccpnmr2.5/python/tests/` | **45 passed, 4 skipped** |
| `ruff check AnalysisPopup.py` | **43** (UP031 17, F841 11, E722 10, E731 2, E721 1, F811 1, W293 1) |

Per-stage rules: import_smoke exit 0 with no NEW unexpected failures and
TOTAL dropping by exactly the removed-module count; gui_boot_test green
(3/3 throughout); pytest floors held; `uvx ruff check` zero NEW.

## Stages & status

| # | Scope | Status |
|---|---|---|
| 17 | Assignment ▶ Initialise Root Resonances + `popups/InitRootAssignments.py` + doc line | ✅ 2026-08-24 |
| 18 | Assignment ▶ Pick & Assign From Roots + `popups/LinkPeakLists.py` + doc line | ✅ 2026-08-24 |
| 19 | Assignment ▶ Protein Sequence Assignment + `popups/LinkSeqSpinSystems.py` + doc line + `EditSpinSystem.py` help link | ✅ 2026-08-24 |
| 20 | Assignment ▶ Automated Seq. Assignment + `ccpnmr/nexus/` (5) + `wrappers/{Mars,Psipred}.py` + doc lines | ✅ 2026-08-24 |
| 21 | Assignment ▶ NOE Contributions + `popups/LinkNoeResonances.py` + doc line + separator collapse | ✅ 2026-08-24 |

## Stage checklist detail

**Stage 17 — Initialise Root Resonances**
- `AnalysisPopup.py`: `InitRootAssignments` import (pre-edit L103);
  `popupActions` entry `"initialise_root_spectra"` (L237); the
  "Initialise Root Resonances" `add_command` block (shortcut "I") +
  its `menuNames.append`; the `initialiseRootSpectra` method.
- Delete `ccpnmr2.5/python/ccpnmr/analysis/popups/InitRootAssignments.py`
  (756 lines, incl. `testInitialiseRootPeakListPopup` argServer macro).
- `doc/source/menu/Assignment.rst`: drop the "Initialise Root Resonances"
  toctree line.
- Gates expected: TOTAL 1250→**1249** (−1), FAILED 2, BY-DESIGN 10,
  gui 3/3, pytest unchanged, ruff zero NEW.
- Commit + push; stage log; mark ✅.

**Stage 18 — Pick & Assign From Roots**
- `AnalysisPopup.py`: `LinkPeakLists` import (pre-edit L106); `popupActions`
  `"link_peaklists"` (L238); menu block (shortcut "P"); `linkPeakLists`
  method.
- Delete `popups/LinkPeakLists.py` (1356 lines — contains the S17
  cross-references at L104/L178).
- `doc/source/menu/Assignment.rst`: drop the "Pick & Assign From Roots" line.
- Gates expected: TOTAL 1249→**1248** (−1).
- Commit + push; stage log; mark ✅.

**Stage 19 — Protein Sequence Assignment**
- `AnalysisPopup.py`: `LinkSeqSpinSystems` import (pre-edit L107);
  `popupActions` `"link_seq_spin_systems"` (L239); menu block (shortcut "S");
  `linkSeqSpinSystems` method.
- Delete `popups/LinkSeqSpinSystems.py` (3173 lines — its L220-221 help-link
  defs to AutoBackbone/LinkPeakLists die with it; incl.
  `LinkSeqSpinSystemsTestMacro`).
- `popups/EditSpinSystem.py`: drop `.. _str(Protein Sequence Assignment):
  LinkSeqSpinSystemsPopup.html` (L165) + the docstring sentence using
  `str(Protein Sequence Assignment)_` (L117 region — verify exact text) if it
  becomes dangling (Stage-12 sweep precedent; keep the rest of the docstring).
- `doc/source/menu/Assignment.rst`: drop the "Protein Sequence Assignment"
  line.
- Gates expected: TOTAL 1248→**1247** (−1).
- Commit + push; stage log; mark ✅.

**Stage 20 — Automated Seq. Assignment (Nexus/MARS/PSIPRED)**
- `AnalysisPopup.py`: `AutoBackbonePopup` import (pre-edit L125);
  `popupActions` "auto_backbone_assign" (L246) + the commented
  `#'auto_backbone_assign': self.activateMars` line (L247); menu block
  (shortcut "u"); `autoBackboneAssign` method + the commented `activateMars`
  block after it.
- Delete `ccpnmr2.5/python/ccpnmr/nexus/` (5 files: `__init__` (bare pass),
  `_licenseInfo`, `NexusBasic` 806, `AutoBackboneNexus` 639,
  `AutoBackbonePopup` 1920 — the last `NexusBasic` consumers,
  `LinkSeqSpinSystems` L70 and `AutoBackbonePopup` L57/L33, are both gone by
  then).
- Delete `ccpnmr2.5/python/ccpnmr/analysis/wrappers/Mars.py` (507; sole
  importer was `AutoBackbonePopup`) + `wrappers/Psipred.py` (188; sole
  importer was `Mars.py`). `wrappers/__init__.py` + `CamCoil` + `D2D` +
  `Shiftx` stay.
- `doc/source/menu/Assignment.rst`: drop the "Automated Seq. Assignment"
  line; `ccpnmr/analysis/doc/Readme.txt:66`: drop the `nexus/` layout line.
- Gates expected: TOTAL 1247→**1240** (−7: 5 nexus + 2 wrappers).
- Commit + push; stage log; mark ✅.

**Stage 21 — NOE Contributions**
- `AnalysisPopup.py`: `LinkNoeResonances` import (pre-edit L105);
  `popupActions` `"link_noe_resonances"` (L240); menu block (shortcut "N");
  `linkNoeResonances` method; collapse the TWO now-adjacent
  `menu.add_separator()` calls left between "Spin System Typing" and
  "Assignment Graph" into one.
- Delete `popups/LinkNoeResonances.py` (1419 lines; its stage-11
  `viewStructure` call site was already removed in Stage 11).
- `doc/source/menu/Assignment.rst`: drop the "NOE Contributions" line.
- Gates expected: TOTAL 1240→**1239** (−1).
- Commit + push; stage log; mark the plan **COMPLETE** (21/21 across all
  three plans).

## Rollback

Each stage is exactly one commit on the new repo → `git revert <sha>` restores
it cleanly.

## Stage log (Assignment Menu Removal)

**Stage 17 — Initialise Root Resonances — ✅ 2026-08-24**
Recon (verified pre-edit): `popups/InitRootAssignments.py` (756 lines) —
`InitRootAssignmentsPopup(BasePopup)` "Initialise Root Resonances" + the
`testInitialiseRootPeakListPopup` argServer macro inside it. Sole external
importer: `AnalysisPopup.py` (import, popupActions entry, menu block
shortcut "I", 3-line `initialiseRootSpectra` method). The in-file test macro
is a standalone `argServer` entry point — no registration table to update
(Stage-14 precedent: macro surface lives in the module). Cross-references in
`LinkPeakLists.py:104,178` (docstring usage sentence + help-link def) die
with that file in Stage 18. `doc/source/menu/Assignment.rst:15` toctree line
(the `source/popups/` target dir doesn't exist — link already dangling).
Zero references in pyproject/MANIFEST/import_smoke/gui_boot_test/bin/,
scripts/, README, INSTALL (grep-verified). `AnalysisProject` model has NO
`initRoot*` data field (unlike LinkPeakLists/LinkSeqSpinSystems — the
`linkPeakListsData` / `linkSeqSpinSystemsData` fields are KEPT, locked
decision 4).
- `AnalysisPopup.py` (−18 lines): dropped the
  `ccpnmr.analysis.popups.InitRootAssignments` import; the
  `"initialise_root_spectra"` popupActions entry; the "Initialise Root
  Resonances" add_command block (label + 6-arg block + `menuNames.append`);
  the `initialiseRootSpectra` method (its `popup =` line is the F841 the
  ruff baseline carried).
- Deleted `ccpnmr2.5/python/ccpnmr/analysis/popups/InitRootAssignments.py`
  (756 lines).
- `doc/source/menu/Assignment.rst`: dropped the toctree line (8 → 7 lines in
  the middle block).
- Grep-verified: repo-wide `InitRootAssignments|initialiseRootSpectra|
  initialise_root_spectra|Initialise Root Resonances` → only the two S18
  cross-references in LinkPeakLists.py (dies next stage) + this plan doc +
  the gitignored `dist/` snapshot.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1249** (1250−1, exactly the one
    removed module), OK **1237**, FAILED **2** unchanged (2× `cherrypy`),
    BY-DESIGN **10** unchanged.
  - `gui_boot_test.py` **3/3** (the CCPN main window boots through the
    edited `setAssignMenu` path).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `uvx ruff check`: AnalysisPopup 43→**42** (F841 11→10 — the removed
    unused `popup =`; UP031 17 / E722 10 / E731 2 / E721 1 / F811 1 / W293 1
    all flat) — zero NEW. `python -m py_compile` OK.
**Stage 18 — Pick & Assign From Roots — ✅ 2026-08-24**
Recon (verified pre-edit): `popups/LinkPeakLists.py` (1356 lines) —
`LinkPeakListsPopup(BasePopup)` + `testPopup` argServer macro; imports
`AssignmentAdvanced.{assignSpecNonRootResonances, pickAssignSpecFromRoot}`
(KEPT core — hazard closed, those functions survive in core along with their
other consumers) and the standard `core/{Assignment,Experiment,Mark,Util,
Window}Basic` helpers (all KEPT, used by many kept popups). Sole external
importer: `AnalysisPopup.py` (import, popupActions `"link_peaklists"`, menu
block shortcut "P", 3-line `linkPeakLists` method). Contains the S17
cross-references (docstring L104 + help-link def L178 to
InitRootAssignmentsPopup — both die with the file). `doc/source/menu/
Assignment.rst` toctree line (already-dangling target). NOTE:
`AnalysisProject.linkPeakListsData` model field is KEPT (locked decision 4 —
generated API/xml/model, plain data attribute; its popup just stops existing).
- `AnalysisPopup.py` (−17 lines): dropped the `ccpnmr.analysis.popups.
  LinkPeakLists` import; the `"link_peaklists"` popupActions entry; the
  "Pick & Assign From Roots" add_command block + `menuNames.append`; the
  `linkPeakLists` method (its unused `popup =` line is the F841 this stage
  removed).
- Deleted `ccpnmr2.5/python/ccpnmr/analysis/popups/LinkPeakLists.py`
  (1356 lines).
- `doc/source/menu/Assignment.rst`: dropped the toctree line (middle block
  now 5 lines).
- Grep-verified: repo-wide `LinkPeakLists|linkPeakLists|link_peaklists|
  Pick & Assign From Roots` (excluding the KEPT `linkPeakListsData` API
  attribute) → only the two S19 cross-references in LinkSeqSpinSystems.py
  (L136 docstring usage sentence + L221 help-link def — die next stage) +
  this plan doc + gitignored `dist/`.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1248** (1249−1, exactly as predicted),
    OK **1236**, FAILED **2** unchanged (2× `cherrypy`), BY-DESIGN **10
    unchanged**.
  - `gui_boot_test.py` **3/3** (the CCPN main window boots through the
    edited `setAssignMenu` path).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `uvx ruff check`: AnalysisPopup 42→**41** (F841 10→9 — the removed
    unused `popup =`; UP031 17 / E722 10 / E731 2 / E721 1 / F811 1 / W293 1
    flat) — zero NEW. `python -m py_compile` OK.
**Stage 19 — Protein Sequence Assignment — ✅ 2026-08-24**
Recon (verified pre-edit, cross-refs confirmed against the HEAD copy of the
removed file): `popups/LinkSeqSpinSystems.py` (3173 lines) —
`LinkSeqSpinSystemsPopup(BasePopup)` (title "Assignment : Protein Sequence
Assignment") + the `LinkSeqSpinSystemsTestMacro` argServer macro (L100);
imports `ccpnmr.nexus.NexusBasic` (L70, `linkSpinSystemInterIntraResonances`
— the SHARED consumer; per locked decision 1 `NexusBasic` goes with S20, the
stage removing its last consumer); reads/writes the KEPT
`AnalysisProject.linkSeqSpinSystemsData` model field (L1231/L1341, locked
decision 4). Sole external importer: `AnalysisPopup.py` (import, popupActions
`"link_seq_spin_systems"`, menu block shortcut "S", 3-line
`linkSeqSpinSystems` method). Carried the two S18 cross-references — L136
docstring usage sentence + L221 help-link def to LinkPeakListsPopup — both
die with the file. Orphans created by the removal: `popups/EditSpinSystem.py`
docstring reference `str(Protein Sequence Assignment)_` (L112) + its
`.. _str(...)` link def (L165); `doc/source/menu/Assignment.rst:15` toctree
line (target dir `source/popups/` doesn't exist — link already dangling,
Stage-12 precedent). Repo-wide grep: beyond the file itself + these orphans,
only the KEPT generated `linkSeqSpinSystemsData` surface (ccpnmr/api,
ccpnmr/xml, memops/api, model XML, generated API docs), a plain-prose phrase
in `BrowseReferenceShifts.py` ("performing protein sequence assignment" —
no role ref, KEPT), this plan doc, and the gitignored `dist/` snapshot.
- `AnalysisPopup.py` (−16 lines): dropped the
  `ccpnmr.analysis.popups.LinkSeqSpinSystems` import; the
  `"link_seq_spin_systems"` popupActions entry; the "Protein Sequence
  Assignment" add_command block (label + 6-arg block + `menuNames.append`) —
  the preceding `menu.add_separator()` KEPT (after S20/S21 the two separators
  around the removed middle block sit adjacent; Stage 21 collapses the pair,
  per the menu-index bookkeeping); the `linkSeqSpinSystems` method (its
  unused `popup =` line carried the F841 this stage removed).
- Deleted `ccpnmr2.5/python/ccpnmr/analysis/popups/LinkSeqSpinSystems.py`
  (3173 lines) via `git rm`.
- `popups/EditSpinSystem.py` (+1/−4): dropped the dangling docstring
  sentence "Such links are independent of a full residue assignment, and
  usually derive from the peak matching performed by tools like the
  str(Protein Sequence Assignment)_ option." (rest of the "Seq. Links"
  bullet kept, per the Stage-12 sweep precedent); dropped the
  `.. _str(Protein Sequence Assignment): LinkSeqSpinSystemsPopup.html`
  link def (the `.. _str(Assignment Panel)` def stays — its target is a
  KEPT popup).
- `doc/source/menu/Assignment.rst` (−1): dropped the toctree line (second
  group 3 → 2 lines; blank-line grouping preserved).
- Grep-verified: repo-wide
  `LinkSeqSpinSystems|linkSeqSpinSystems|link_seq_spin_systems` → only the
  KEPT generated `linkSeqSpinSystemsData` attribute surface + this plan doc
  + gitignored `dist/`.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1247** (1248−1, exactly the one
    removed module), OK **1235**, FAILED **2** unchanged (2× `cherrypy`),
    BY-DESIGN **10** unchanged.
  - `gui_boot_test.py` **3/3** (the CCPN main window boots through the
    edited `setAssignMenu` path).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `uvx ruff check AnalysisPopup.py`: 41→**40** (F841 9→8 — the removed
    unused `popup =`; UP031 17 / E722 10 / E731 2 / E721 1 / F811 1 / W293 1
    flat) — zero NEW. `python -m py_compile` OK (AnalysisPopup +
    EditSpinSystem).
**Stage 20 — Automated Seq. Assignment (Nexus/MARS/PSIPRED) — ✅ 2026-08-24**
Recon (verified pre-edit): `ccpnmr/nexus/` (5 files, 3428 lines —
`__init__` 1 / `_licenseInfo` 62 / `NexusBasic` 806 / `AutoBackboneNexus` 639
/ `AutoBackbonePopup` 1920) — `AutoBackbonePopup(BasePopup)` (title
"Assignment : Automated Seq. Assignment"); imports `ccpnmr.analysis.wrappers.
Mars.runMars` (L56) + `ccpnmr.nexus.NexusBasic` (L57) + lazy
`AutoBackboneNexus.autoBackboneNexus` (L33). With the S19 removal of
`LinkSeqSpinSystems` (its L70 `NexusBasic` import), `AutoBackbonePopup` was
`NexusBasic`'s LAST consumer (locked decision 1) → `nexus/` dies here. No
`argServer` macro inside `nexus/` (grep-verified). `wrappers/Mars.py` (507):
SOLE importer `AutoBackbonePopup` (`runMars`); carries the `tesMars`
`argServer` macro (L81) — a standalone in-file entry point, no registration
table to update (Stage-14 precedent); `__main__` self-test at L88.
`wrappers/Psipred.py` (188): SOLE importer `Mars.py`'s lazy `psipredCcpn`
(L223); carries the `testCcpnPsipred` `argServer` macro (L54) — same
standalone-macro pattern. Repo-wide grep (`AutoBackbone|auto_backbone_
assign|NexusBasic|AutoBackboneNexus|runMars|Psipred|psipred|ccpnmr\.nexus`):
beyond these files, only `AnalysisPopup.py` (the 5 items below), the two
doc lines, untracked/generated `build/`+`dist/`+`ccpnmr.egg-info/` snapshots
(untracked — `git ls-files | grep -c egg-info` = 0 — same as S17-19, no
stage-commit scope), and the unrelated `cambridge/dangle` commented-out path
+ `memops/general/license` `_licenseInfo` machinery (a DIFFERENT
`_licenseInfo` concept — the memops license generator). `wrappers/{CamCoil,
D2D,Shiftx}.py` + `wrappers/__init__.py` KEPT (locked decision 2 — CamCoil
→ `SequenceShiftPredict` and D2D → `SecStructurePredict` are kept popups).
- `AnalysisPopup.py` (−30 lines, 2987→2957): dropped the
  `ccpnmr.nexus.AutoBackbonePopup` import; the `"auto_backbone_assign":
  self.autoBackboneAssign` popupActions entry + the commented
  `#'auto_backbone_assign': self.activateMars` line; the "Automated Seq.
  Assignment" add_command block (label + 6-arg block + `menuNames.append`) —
  the preceding `menu.add_separator()` KEPT (after S21 the two separators
  around the removed middle block sit adjacent; Stage 21 collapses the pair,
  per the menu-index bookkeeping); the `autoBackboneAssign` method + the
  commented `activateMars` block after it (7 commented lines incl. the
  commented `popup = self.popups['auto_backbone_assign']`). `HAVE_NUMPY` KEPT
  (locked decision 5 — `peakSeparatorParams` uses it); `showWarning` KEPT
  (6 remaining uses in-file).
- `git rm` (7 files, 4123 lines): `ccpnmr/nexus/{__init__,_licenseInfo,
  NexusBasic,AutoBackboneNexus,AutoBackbonePopup}.py` (3428) +
  `ccpnmr/analysis/wrappers/{Mars,Psipred}.py` (507+188).
- `doc/source/menu/Assignment.rst` (−1): dropped the toctree line (middle
  group 2 → 1 line — "NOE Contributions" now stands alone; blank-line
  grouping preserved; its line dies in S21).
- `ccpnmr/analysis/doc/Readme.txt` (−2): dropped the `nexus/` layout line +
  its separating blank line. (Observation, out of scope: the `clouds/`
  line in the same block survived S15 — pre-existing, untouched.)
- Grep-verified post-edit: zero `AutoBackbone|auto_backbone|autoBackbone|
  activateMars|nexus` references in `AnalysisPopup.py`.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1240** (1247−7: 5 nexus + 2 wrappers,
    exactly as predicted), OK **1228**, FAILED **2** unchanged (2×
    `cherrypy`), BY-DESIGN **10** unchanged.
  - `gui_boot_test.py` **3/3** (the CCPN main window boots through the
    edited `setAssignMenu` path).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `uvx ruff check AnalysisPopup.py --statistics`: 40→**40** — IDENTICAL
    mix (UP031 17 / E722 10 / F841 8 / E731 2 / E721 1 / F811 1 / W293 1) —
    zero NEW (no F841 delta: the removed method had no unused `popup =`,
    and no import became orphaned). `python -m py_compile` OK.
**Stage 21 — NOE Contributions — ✅ 2026-08-24**
Recon (verified pre-edit): `popups/LinkNoeResonances.py` (1419 lines) —
`LinkNoeResonancesPopup(BasePopup)` (title "Assignment : NOE
Contributions") + the `testNoePopup` `argServer` macro (L74 — standalone
in-file entry point, no registration table to update; Stage-14
precedent). Sole external importer: `AnalysisPopup.py` (import L104,
popupActions `"link_noe_resonances"` L233, menu block L1191-1202 with
shortcut "N", 2-line `linkNoeResonances` method L2156-2158). Its
stage-11 `viewStructure` call site was already removed in Stage 11.
`doc/source/menu/Assignment.rst:15` toctree line (target dir
`source/popups/` doesn't exist — link already dangling; S17-20
precedent). `doc/Changes.html` mentions (L243/L411/L437) are changelog
history — KEPT (locked decision 6, Stage-12 precedent). Repo-wide grep
(`LinkNoeResonances|linkNoeResonances|link_noe_resonances|NOE
Contributions`): beyond the 5 items + doc line above + these, only this
plan doc and the gitignored `dist/` snapshot.
- `AnalysisPopup.py` (−17 lines, 2957→2940): dropped the
  `ccpnmr.analysis.popups.LinkNoeResonances` import; the
  `"link_noe_resonances"` popupActions entry; the "NOE Contributions"
  add_command block (label + 6-arg block + `menuNames.append`); the
  `linkNoeResonances` method (its unused `popup =` line carried the F841
  this stage removed). **Separator collapse (S19/S20 bookkeeping):** the
  TWO now-adjacent `menu.add_separator()` calls between "Spin System
  Typing" and "Assignment Graph" collapsed into ONE — the Assignment menu
  is now `0 Assignment Panel, 1 Copy Assignments, 2 Spin System Typing,
  3 sep, 4 Assignment Graph, 5 Quality Reports`, exactly the plan's
  "After S21" bookkeeping line.
- Deleted `ccpnmr2.5/python/ccpnmr/analysis/popups/LinkNoeResonances.py`
  (1419 lines) via `git rm`.
- `doc/source/menu/Assignment.rst` (−2): dropped the "NOE Contributions"
  toctree line + one separating blank line (it stood alone between two
  blank lines after S20) — the toctree keeps two well-grouped sections:
  (Assignment Panel / Copy Assignments / Spin System Typing), then
  (Assignment Graph / Quality Reports).
- Grep-verified post-edit: zero
  `LinkNoeResonances|linkNoeResonances|link_noe_resonances` references in
  any code or doc outside this plan doc.
- ENV note (not a code regression): the venv had lost the dev/test extras
  (`pytest`/`pytest-cov` from `testing`, `decorator` from `optional`)
  since the S20 baseline — this session's first import_smoke run showed
  FAILED 6 (the known 2× `cherrypy` + 3× `tests/test_*.py` on missing
  `pytest` + `TestNefIo.py` on missing `decorator`; none of the 4 extra
  modules relate to this stage's removal). Restored via
  `uv pip install decorator pytest pytest-cov` (no `uv.lock` touched);
  the gate run below is the post-restore result.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1239** (1240−1, exactly the one
    removed module), OK **1227**, FAILED **2** unchanged (2×
    `cherrypy`), BY-DESIGN **10** unchanged.
  - `gui_boot_test.py` **3/3** (the CCPN main window boots through the
    edited `setAssignMenu` path).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `uvx ruff check AnalysisPopup.py --statistics`: 40→**39** (F841 8→7 —
    the removed unused `popup =`; UP031 17 / E722 10 / E731 2 / E721 1 /
    F811 1 / W293 1 flat) — zero NEW. `python -m py_compile` OK.
- **ASSIGNMENT MENU REMOVAL PLAN COMPLETE — 21/21** (Simplification 1-12
  + Menu Removal 13-16 + Assignment Menu 17-21), all stages committed +
  pushed to `21tesla/analysis2.5py3lite:main`.

## Menu Removal Round 2 (Stages 22–23) — requested 2026-08-24

| # | Scope | Status |
|---|---|---|
| 22 | Molecule menu ▶ "Isotope Labelling" + "Reference Isotope Schemes" + `popups/{EditMolLabelling,IsotopeSchemeEditor}.py` + doc lines | ✅ 2026-08-24 |
| 23 | Assignment menu ▶ "Spin System Typing" + `SpinSystemTypingPopup`/`TypingEnsemblePopup` classes + `core/SpinSystemTyping.py` + doc line (`SpinSystemTypeScoresPopup` KEPT — live callers below) | ✅ 2026-08-24 |

### Stage 22 — Molecule menu: Isotope Labelling + Reference Isotope Schemes

Recon (verified 2026-08-24):
- `AnalysisPopup.py` items: `EditMolLabellingPopup` import (L91);
  `IsotopeSchemeEditor` import (L103); `popupActions` `"isotope_scheme_editor"`
  (L225) + `"isotope_labelling"` (L226); the two `setMoleculeMenu`
  `add_command` blocks ("Isotope Labelling" L1082-1089 with shortcut "L",
  "Reference Isotope Schemes" L1090-1097 with shortcut "I"); the two
  `menu_items[MoleculeMenu]` entries (L1129-1130); the `isotopomerEditor`
  (L2102-2104) + `editIsotopeLabelling` (L2106-2108) methods.
- `popups/EditMolLabelling.py` (1860 lines): `EditMolLabellingPopup` + 4
  module-level helpers (`getMolLabelFromScheme`, `updateResLabelFractions`,
  `setResLabelNaturalAbundance`, `getMolLabelMass`) — each used ONLY
  in-file (grep-verified, `ccpnmr/` + `memops/`, .py, minus pycache).
- `popups/IsotopeSchemeEditor.py` (1304 lines): `IsotopeSchemeEditor` +
  helpers `getSortedIsotopes` (sole importer `EditMolLabelling` L82/L798/
  L1093 — dead pair), `getPrimaryIsotope` + `isSchemeEditable` (in-file
  only); the `testChemCompLabelEditorMacro` `argServer` macro (L79) and the
  `__main__` self-test (L1269) die with the file (standalone-macro
  precedent, no registration table to update).
- `doc/source/menu/Molecule.rst` L12-13: two toctree lines (target dir
  `source/popups/` doesn't exist — links already dangling; S17-21
  precedent).
- Grep-verified: repo-wide
  `isotopomerEditor|editIsotopeLabelling|isotope_labelling|isotope_scheme_
  editor|EditMolLabelling|IsotopeSchemeEditor` → only the items above +
  the two in-file docstring help-link cross-refs (L231/L236, dead pair) +
  the gitignored `build/`/`dist/` snapshots.
- Orphan OBSERVATION (NOT removed — locked hazard-3 "keep `ccpnmr/
  analysis/core/*` + related frames" policy): `frames/ViewIsotopomerFrame.
  py` (subclasses KEPT `ViewChemCompVarFrame`) loses its only importer;
  candidate for a future stage if the policy changes.
- Gates expected: TOTAL 1239→**1237** (−2), FAILED 2, BY-DESIGN 10, gui
  3/3, pytest 45/4, ruff 39 flat (no F841 delta — both removed methods
  use bare `self.openPopup(...)`, no unused `popup =`).
- Commit + push; stage log; mark ✅.

### Stage 23 — Assignment menu: Spin System Typing

Recon (verified 2026-08-24):
- `AnalysisPopup.py` items: import (L112) `from ccpnmr.analysis.popups.
  SpinSystemTyping import SpinSystemTypeScoresPopup, SpinSystemTypingPopup`
  → trim to `SpinSystemTypeScoresPopup` ONLY; `popupActions`
  `"type_spin_systems"` (L233) — `"type_spin_system"` (L234) KEPT; the
  "Spin System Typing" `add_command` block (label L1178, shortcut "T", +
  `menuNames.append`); the `typeSpinSystems` method (L2147-2149 — its
  unused `popup =` carries the F841 this stage removes).
- `popups/SpinSystemTyping.py` (1300 lines): remove
  `class SpinSystemTypingPopup` (L698-1223) + `class TypingEnsemblePopup`
  (L1225-1300, ends the file) + `COLOR_DICT` (L68, only use L825 in the
  removed class) + now-orphan imports: `core.SpinSystemTyping.
  getSpinSystemTypes` (L53), `MoleculeBasic.getResidueCode` (L52 — keep
  `DEFAULT_ISOTOPES`), `memops.gui.{CheckButton,FloatEntry,IntEntry,
  PartitionedSelector,ProgressBar,ScrolledGraph}`. No `argServer` macro or
  `__main__` block in-file.
- `git rm` `core/SpinSystemTyping.py` (541 lines; `getSpinSystemTypes`
  L121 has no importer beyond the removed popup class — grep-verified).
- **KEPT (verified live call sites — do NOT touch):**
  `SpinSystemTypeScoresPopup` (L71-696, the file's first class) — called
  by KEPT `EditSpinSystem.predictType` (L606 →
  `self.guiParent.typeSpinSystem(self.spinSystem)`) and KEPT
  `WindowFrame.predictSpinSystemType` (L6510 →
  `self.topPopup.typeSpinSystem(spinSystem=..., shiftList=...)`);
  `AnalysisPopup.typeSpinSystem` (singular, L2151); `MoleculeBasic.
  DEFAULT_ISOTOPES` (still used popup-file L284/L302); the
  "Spin System Type Scores" toctree lines in `doc/source/menu/{
  Assignment,Resonance}.rst` (popup still reachable); `Changes.html`
  changelog mentions (history — locked decision 6).
- `doc/source/menu/Assignment.rst`: drop the "Spin System Typing" toctree
  line (L13; first group 3→2 lines).
- Assignment-menu bookkeeping after S23: `0 Assignment Panel, 1 Copy
  Assignments, 2 sep, 3 Assignment Graph, 4 Quality Reports` (no adjacent
  separators created).
- Gates expected: TOTAL 1237→**1236** (−1), FAILED 2, BY-DESIGN 10, gui
  3/3, pytest 45/4, ruff 39→**38** (F841 7→6).
- Commit + push; stage log; mark ✅.

## Stage log (Menu Removal Round 2)

**Stage 22 — Molecule menu: Isotope Labelling + Reference Isotope Schemes — ✅ 2026-08-24**
- `AnalysisPopup.py` (−30 lines, 2940→2910): dropped the
  `EditMolLabellingPopup` import (pre L91); the `IsotopeSchemeEditor` import
  (pre L103); the `popupActions` entries `"isotope_scheme_editor"` +
  `"isotope_labelling"` (pre L225-226); the two `setMoleculeMenu`
  `add_command` blocks ("Isotope Labelling" shortcut "L" + "Reference
  Isotope Schemes" shortcut "I") — the menu is now `Molecules | sep | Atom
  Browser, Add Sequence | sep | Residue Information` with NO adjacent
  separators; the two `menu_items[MoleculeMenu]` entries (6→4); the
  `isotopomerEditor` (pre L2102-2104) + `editIsotopeLabelling` (pre
  L2106-2108) methods (both bare `self.openPopup(...)` — zero F841 delta).
  The `# Isotopomer scheme tidy` TBD planning comment in `setMoleculeMenu`
  KEPT (non-functional note, minimal-scope precedent).
- `git rm` (2 files, 3164 lines): `popups/EditMolLabelling.py` (1860 —
  incl. the 4 in-only module helpers + the L236 `.. _str(Reference Isotope
  Schemes)` help-link def + its `getSortedIsotopes` import, all dead with
  the file) + `popups/IsotopeSchemeEditor.py` (1304 — incl. the
  `testChemCompLabelEditorMacro` `argServer` macro (L79) and the `__main__`
  self-test (L1269) — both die with the file (standalone-macro
  precedent) — the in-only `getSortedIsotopes`/`getPrimaryIsotope`/
  `isSchemeEditable`
  helpers (sole importer of `getSortedIsotopes` was `EditMolLabelling` —
  dead pair), and the L231 `.. _str(Isotope Labelling)` help-link def).
- `doc/source/menu/Molecule.rst` (−2): dropped the two toctree lines —
  "Molecules" now stands alone (its target dir `source/popups/` didn't
  exist; link already dangling, S17-21 precedent); "Atom Browser / Add
  Sequence" group intact.
- Grep-verified post-edit: `EditMolLabelling|IsotopeSchemeEditor|
  isotopomerEditor|editIsotopeLabelling|isotope_labelling|isotope_scheme_
  editor` across `ccpnmr2.5/python/**/*.py` + `**/*.rst` → zero hits
  (pycache aside).
- Orphan OBSERVATION (NOT removed — locked core/frame-keep policy):
  `frames/ViewIsotopomerFrame.py` (subclasses KEPT
  `ViewChemCompVarFrame`) is now importerless (its only importer was
  `IsotopeSchemeEditor`) — flagged as a candidate for a future stage if
  that policy changes.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1237** (1239−2, exactly the two
    removed modules), OK **1225**, FAILED **2** unchanged (2× `cherrypy`),
    BY-DESIGN **10** unchanged.
  - `gui_boot_test.py` **3/3** (the CCPN main window boots through the
    edited `setMoleculeMenu` path).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `uvx ruff check AnalysisPopup.py --statistics`: 39→**39** — IDENTICAL
    mix (UP031 17 / E722 10 / F841 7 / E731 2 / E721 1 / F811 1 / W293 1) —
    zero NEW. `python -m py_compile` OK.

**Stage 23 — Assignment menu: Spin System Typing — ✅ 2026-08-24**
- `AnalysisPopup.py` (−15 lines, 2910→2895): trimmed the import (pre
  L112: `SpinSystemTypeScoresPopup, SpinSystemTypingPopup` →
  `SpinSystemTypeScoresPopup` ONLY); dropped the `popupActions` entry
  `"type_spin_systems"` (pre L233; `"type_spin_system"` KEPT); dropped the
  "Spin System Typing" `add_command` block (label + shortcut "T" +
  `menuNames.append`); dropped the `typeSpinSystems` method (pre
  L2147-2149 — its unused `popup =` carried the F841 this stage removed).
  Assignment menu is now `0 Assignment Panel, 1 Copy Assignments, 2 sep,
  3 Assignment Graph, 4 Quality Reports` — no adjacent separators.
- `popups/SpinSystemTyping.py` (1300→687, −613): removed
  `class SpinSystemTypingPopup` (pre L698-1223) + `class
  TypingEnsemblePopup` (pre L1225-1300, ended the file) via an
  assert-guarded line-range truncation (boundary verified: exactly two
  blank lines after the KEPT scores popup's `BasePopup.destroy(self)`);
  dropped `COLOR_DICT` (pre L68 — only use L825, inside the removed
  class); dropped the 7 now-orphan imports: `core.SpinSystemTyping.
  getSpinSystemTypes`; `MoleculeBasic.getResidueCode` (`DEFAULT_ISOTOPES`
  on the same line KEPT — still used at scores-popup L284/L302);
  `memops.gui.{CheckButton,FloatEntry,IntEntry,PartitionedSelector,
  ProgressBar,ScrolledGraph}`. `SpinSystemTypeScoresPopup` (pre
  L71-696) intact — now the file's only class. Its KEPT callers
  `EditSpinSystem.predictType` (L606) + `WindowFrame.
  predictSpinSystemType` (L6510) both route through `AnalysisPopup.
  typeSpinSystem` (singular, kept).
- `git rm` `ccpnmr2.5/python/ccpnmr/analysis/core/SpinSystemTyping.py`
  (541 lines — `getSpinSystemTypes` had no importer beyond the removed
  popup class; `core/__init__.py` is empty — no re-exports to clean).
- `doc/source/menu/Assignment.rst` (−1): dropped the "Spin System
  Typing" toctree line (first group 3→2 lines). "Components available
  indirectly → Spin System Type Scores" section KEPT in both
  `Assignment.rst` and `Resonance.rst` (popup still live); `Changes.html`
  mentions are changelog history (locked decision 6).
- Grep-verified post-edit:
  `SpinSystemTypingPopup|TypingEnsemblePopup|typeSpinSystems|
  core.SpinSystemTyping|getSpinSystemTypes|type_spin_systems` across
  `ccpnmr/**/*.py` + `**/*.rst` → ZERO hits (pycache aside); the
  kept-surface grep shows exactly the 6 expected live refs (class def,
  the two frame/popup callers, the AnalysisPopup import + method).
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1236** (1237−1, exactly the one
    removed module), OK **1224**, FAILED **2** unchanged (2× `cherrypy`),
    BY-DESIGN **10** unchanged.
  - `gui_boot_test.py` **3/3** (the CCPN main window boots through the
    edited `setAssignMenu` path).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `uvx ruff check AnalysisPopup.py --statistics`: 39→**38** (F841 7→6 —
    the removed unused `popup =`; UP031 17 / E722 10 / E731 2 / E721 1 /
    F811 1 / W293 1 flat) — zero NEW. `popups/SpinSystemTyping.py` (newly
    in gate scope, KEPT portion only): 10 PRE-EXISTING (UP031 9 / F841 1)
    — zero NEW, no F401 (import scrub complete). `python -m py_compile`
    OK (both files).
- **MENU REMOVAL ROUND 2 COMPLETE — 2/2** (Stages 22-23), both stages
  committed + pushed to `21tesla/analysis2.5py3lite:main`.

**Stage 24 — Deep scan: module / C-ext / class liveness audit (READ-ONLY) — ✅ 2026-08-24**
- Scope (locked with user 2026-08-24): (1) the 3 generated API files —
  **audit first, then prune**; (2) removals = fully-orphaned Python modules +
  dead C extensions (NO symbol-level sweep inside kept files); (3) same
  checkpoint policy as S1-23 (gates green → one commit code+plan doc → push to
  `21tesla/analysis2.5py3lite:main` via `refs/heads/master:refs/heads/main`).
  **This stage changed no code** — audit only.
- Method: static AST import graph over git-tracked `.py` modules (1,235) +
  string/token scan; live entry roots: `ccpnmr.analysis.AnalysisGui`,
  `ccpnmr.format.gui.DataShifter`, `ccpnmr.format.gui.FormatConverter`,
  `ccpnmr.update.UpdatePopup` + `tests/`; 37 C-ext `.so` bridges audited by
  static importer + bridge-dir mapping (`ccpnmr2.5/c/` sources); generated-API
  classes audited per-class. Artifacts: `/tmp/deepscan/{audit.py,report.json}`
  (re-runnable: `python3 /tmp/deepscan/audit.py`).
- Result: **553 live** (import closure from the live entry roots), **333 strong
  orphans** (zero static, dynamic or data refs), 20 dynamic-ref-only (mostly
  `import_smoke.py` allowlist self-mentions + cross-refs inside dead trees).
- **False-positive classes caught & corrected during audit (do NOT re-flag as
  dead in later stages):**
  - `from X import Y` where `X.Y` IS the submodule (ImportFrom → subpackage edge).
  - **Format-registry dynamic loading**: `__import__("ccpnmr.format.converters.%sFormat"
    % label)` at `ccpnmr/format/general/Conversion.py:274`,
    `general/Util.py:163`, `gui/ImportExportFormatPopup.py:276`,
    `gui/ProcFilePopup.py:302`, `converters/Tool.py:143`, labels sourced from live
    `ccpnmr/format/general/Constants.allFormatsDict` → the whole `ccp/format/*` (93)
    + `ccpnmr/format/converters/*` (38) tree is **FUNCTIONAL LIVE SURFACE** of the
    kept FormatConverter/DataShifter apps — NOT removable without an explicit
    "drop format X" decision.
  - C-exts imported `from <pkg>.c import <Name>` (`.so` ∉ module graph) — check by
    grep, not the graph.
  - `_licenseInfo.py` family is dynamic-live via `memops/general/license/headers.py`
    (`infoFileName="_licenseInfo"`; Credits/TextWriter_py_2_1 live) — KEEP.
- C extensions (37 audited):
  - **DEAD: `ccpnmr2.5/c/ccpnmr/dynamics/` (24 .c/.h sources) + bridge
    `python/ccpnmr/c/{DyAtomCoord,DyAtomCoordList,DyDistConstraint,
    DyDistConstraintList,DyDistForce,DyDynamics,Dynamics}.so`** — zero Python
    importers repo-wide; no `ccpnmr/dynamics/` python package exists.
  - **DEAD (isolated): `ccpnmr.c.CloudUtil` bridge + sources under
    `ccpnmr2.5/c/ccpnmr/clouds/`** — zero importers; C-level `#include` coupling
    from other `c/ccpnmr/clouds/` files to be verified at stage time before removal.
  - **ZERO-IMPORTER BRIDGES NOT IN THE SIGNED-OFF S25 SCOPE** (audit-found,
    verify-by-grep at stage time): `ccpnmr.c.{AtomCoord, DistConstraint,
    DistConstraintList, DistForce, Midge}`. Note `ccpnmr.c.AtomCoordList` HAS one
    importer — `ccpnmr/analysis/popups/EditResStructures.py` (itself a S27
    removal) → verify at S25/S27.
  - KEPT (live importers verified): `ccp/structure` (StructAtom/Bond/Structure ←
    live gui/lib/core; StructUtil ×2), `ccpnmr/analysis` (ContourFile/Levels/Style,
    PeakList/PeakCluster, SliceFile, WinPeakList), `ccpnmr/clouds` minus CloudUtil
    (AtomCoordList, Bacus ← live `clouds/ResonanceIdentification.py:392`),
    `memops/global` (MemCache, FitMethod, GlHandler, PdfHandler, PsHandler,
    BlockFile, ShapeFile, StoreFile, StoreHandler, TkHandler),
    `other/cambridge` (BayesPeakSeparator ← live `cambridge/bayes/
    PeakSeparator.py:14`, consumed by live AnalysisPopup).
- 3 generated API files (audit-first per locked scope) — **VERDICT: HOLD, NO safe
  class prunes**:
  - `ccp/api/nmr/Nmr.py` (150,534 ln / 85 classes), `ccp/api/nmr/
    NmrConstraint.py` (47,375 / 35), `ccpnmr/api/Analysis.py` (42,007 / 28) — all
    AUTOGENERATED (PyFileApiGen 1.57.2.1, Data Model 2.1.2); model XML in
    `ccpnmr2.5/model/`; the generator is NOT in the repo.
  - Per-class audit of all 148 classes: EVERY class is referenced from live code —
    sibling generated files (`ccp/api/general/{Method,Citation,Template}.py`,
    `ccp/api/molecule/{MolSystem,ChemComp,MolStructure}.py`,
    `ccp/api/nmr/{NmrCalc,NmrScreen,NmrExpPrototype,NmrReference}.py`,
    `ccpnmr/api/{AnalysisLayout,AnalysisProfile,AnalysisV3,AnalysisWindow}.py`),
    the live `ccp/xml/*` + `ccpnmr/xml/*` builder layer, and live app core.
  - Pruning any class = a metamodel-level campaign (model XML + api + xml-builders
    + compatibility MapInfo + docs) that breaks old-project loading. Optional
    follow-up campaign on app-unreferenced families (Nmr.py *Derivation*/Prob/
    ChainState groups) ONLY with an explicit data-compat break.
- Orphan clusters (333, grouped for the signed-off removal stages):
  - `ccpnmr/integrator/` core + ALL plugins (~26: Aria, Asdp, Cosmos, Isd,
    MultiStruc, Rosetta, Talos, Unio) — stranded when the S3/S4/S6/S7 menus went;
    supersedes S3 "KEPT plugins/Aria" + S9 `workflow/Cing` (both now orphaned;
    `ccpnmr/workflow/{Aria,Cing}` included) → **S26**.
  - `gothenburg/` (init + `mdd`) — only referrer is live `popups/
    EditExperiment.py:1259 from gothenburg import Usf3Io`, **Usf3Io does not
    exist** → already-broken dead code (remove gothenburg + the broken import)
    → **S26**.
  - dead popups: `popups/{BrowseStrucGen, CalcCouplings, EditNoeClasses,
    EditResStructures, LinkSideChains, PredictKarplus, SecStructurePredict,
    SequenceShiftPredict}` + `wrappers/Shiftx` → **S27** (note: `EditResStructures`
    is the sole static importer of live-kept C-ext `ccpnmr.c.AtomCoordList` — its
    death makes that bridge an orphan too; decide at S27).
  - `ccpnmr/analysis/frames/` **21** modules (AxisLabelList … ViewIsotopomerFrame) —
    the rest of the frames package is live (8+ importers) → **S27**
    (ViewIsotopomerFrame: already flagged importerless-kept in S22).
  - `ccpnmr/analysis/macros/` 12 (T2Macro, RelaxationAnalysis, …) — macro ENGINE
    is live (data-model-driven: project Macro entities name modules at run time,
    NO `.mac` files tracked) → removal safe for fresh projects but **breaks old
    projects referencing those names** (data-compat decision) → **S31 HELD**.
  - `ccp/examples/` 35 + `ccpnmr/format/examples` 6 (workshops/help_doc) → **S28**.
  - `ccpnmr/nef/testing` 2 (nose-era; live pytest covers nef) → **S28**.
  - `nijmegen/CASD` 4 (BY-DESIGN; broken — imports deleted `pdbe.deposition`) +
    `pdbe/adatah/{CasdNmr,Pdb}`; `cambridge/isd` 2 (BY-DESIGN ISD_ROOT) → **S29**
    (+ `import_smoke.py` allowlist shrink 10 → ≤4, TOTAL floor drops).
  - `pdbe` tooling: `chemComp` 9, rest of `adatah`, `nmrStar/IO/Ccpn_To_NmrStar`,
    `software` 1, `xml/Util` 1 — but `pdbe.nmrStar.IO` core (Util/Ccpn2NmrStar
    etc.) has built-in `__import__("%s" % name)` chains → verify per-module →
    **S30**; + `regensburg/auremol/AuremolFrame` (its IO/glue ARE registry-live →
    remove Frame only or nothing) → **S30**.
  - `ccpnmr` root helpers: `Common`, `Constants`, `SafeFilename`, `_serverCheck`
    (zero importers) → **S27**; `ccpnmr/v2io/{TestNefIo,Constants}` (S12-kept
    reference — decide) → **S28**; `memops/{editor(4), gui(9), general(6)}` +
    `ccp.{lib(4), util(3), general(1), gui(1), math(1)}` → **S30**.
  - `memops/format/compatibility` 41 + `memops/format/xml.{Compatibility,XmlGen,
    XmlIO}` — project version-migration surface (old saved projects) —
    **S31 HELD**.
- **STAGE PLAN (user sign-off 2026-08-24):** S25 dead C exts; S26
  integrator/workflow/gothenburg; S27 dead popups + frames 21 + root helpers;
  S28 examples 41 + nef.testing; S29 CASD/ISD (+allowlist shrink); S30 pdbe
  tooling + regensburg + memops helpers. **S31 (macros 12 + compatibility/xml 44)
  HELD — data-compat break (old saved projects) — requires an explicit decision.**
- Gate baseline (post-S23, unchanged — no code touched): `import_smoke.py`
  TOTAL **1236** / OK **1224** / FAILED **2** (cherrypy) / BY-DESIGN **10**;
  `gui_boot_test.py` **3/3**; `pytest ccpnmr2.5/python/tests/` **45 passed,
  4 skipped**; `uvx ruff` **0.16.4** (per-file baselines kept in S22/S23 entries).

**Stage 25 — Dead C extensions: dynamics tree + CloudUtil + 5 extra dead clouds exts — ✅ 2026-08-24**
- Scope: the signed-off S25 "dead C exts (dynamics 24 + CloudUtil)" — plus the
  6 extra zero-importer bridges the S24 audit flagged as verify-by-grep
  (`AtomCoord`, `DistConstraint`, `DistConstraintList`, `DistForce`, `Midge`,
  `Dynamics` (non-Dy)) — all 13 re-verified this stage with a repo-wide
  `grep -E` over every `*.py` (python tree, tests, import_smoke, gui_boot):
  **zero importers** → included under the signed-off "dead C extensions" scope.
- Verification before removal:
  - KEPT-source `#include` check (`py_atom_coord.c`, `py_atom_coord_list.c`,
    `py_bacus.c`, `atom_coord_list.c`): no coupling to any removed header.
  - Build-config sweep: `setup.py` (hit — see below), `MANIFEST.in` (no
    refs), `pyproject.toml` (no refs), `scripts/` (no refs),
    `c/environment*.txt` (no refs). Also found: `memops/c/copySharedObjsMac`
    (missed by the earlier 3-script sweep — trimmed).
- **Removed (59 git-tracked files):**
  - `ccpnmr2.5/c/ccpnmr/dynamics/` — ENTIRE tree, 31 files: 24 .c/.h sources
    (12 pairs), `Makefile`, 6 prebuilt `Dy*.so`.
  - `ccpnmr2.5/c/ccpnmr/clouds/` — 28 files: 7 `.so` (AtomCoord, CloudUtil,
    DistConstraint, DistConstraintList, DistForce, Dynamics, Midge) + 21
    sources (`py_cloud_util.c`; the py_/core `.c/.h` pairs for
    dist_constraint, dist_constraint_list, dist_force, dynamics, midge).
- **Removed (26 UNTRACKED local build artifacts):** the 13 dead-name `.so`
  bridges under `ccpnmr2.5/python/ccpnmr/c/` (7 names × plain + cpython-313
  + 6 Dy* × 2) — plain `rm`, never git-tracked (gitignored by `*.so`);
  remaining `.so` set = exactly the 9 kept exts (AtomCoordList, Bacus,
  ContourFile/Levels/Style, PeakList, PeakCluster, SliceFile, WinPeakList).
- **Edited:**
  - `setup.py` (−120, 316→196): docstring family list "clouds / dynamics /
    analysis…" → "clouds / analysis…"; dropped `DYN =` source dir; FAM:
    dropped all 13 dead entries (+ the `# --- dynamics` comment line). FAM
    clouds section is now exactly `AtomCoordList` + `Bacus`.
  - `python/ccpnmr/c/linkSharedObjs` (−13 lines): 7 analysis + `clouds/
    AtomCoordList` + `clouds/Bacus`.
  - `python/ccpnmr/c/copySharedObjs` (−12), `copySharedObjs.bat` (−6),
    `memops/c/copySharedObjsMac` (−12): same 9-line kept set.
  - `c/ccpnmr/clouds/Makefile` (150→49): `all:` = `AtomCoordList` + `Bacus`
    only; dropped DYNAMICS_OBJS / PY_DIST_* / PY_MIDGE / PY_CLOUD_UTIL groups,
    the `global_diag_objects` / `global_random_objects` phony targets, and the
    dead `py_*` compile rules.
- **KEPT (deliberate, with reason):**
  - `c/ccpnmr/clouds/{atom_coord,py_atom_coord}.{c,h}` — BUILD DEPENDENCY of
    the kept `AtomCoordList` ext (setup.py source list). Its sole importer is
    `popups/EditResStructures.py` — itself a **Stage 27 removal** → when it
    dies, the `AtomCoordList` bridge + these 4 sources also become orphans;
    decide at S27 (flagged in S24 entry).
  - `Bacus` (live `ccpnmr.clouds.ResonanceIdentification`), `AtomCoordList`
    (see above).
  - `c/ccpnmr/clouds/_licenseInfo.py` (stray in the C dir) — `_licenseInfo`
    family is dynamic-live by design (S24 finding); out of scope.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1236** / OK **1224** / FAILED **2**
    unchanged (cherrypy) / BY-DESIGN **10** unchanged (no import-graph `.py`
    touched).
  - `gui_boot_test.py` **3/3** (all apps boot through the live C-ext set).
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `uvx ruff check setup.py --statistics`: 1 PRE-EXISTING I001 (this diff
    touches no import lines) — zero NEW; `python -m py_compile setup.py` OK.
  - Residual grep (excl. `dist/` build snapshot + plan doc): ZERO references
    to any of the 13 dead names anywhere in the tree.

**Stage 26 — Dead clusters: integrator/ (40) + workflow/ (6) + gothenburg/ (2) — ✅ 2026-08-24**
- Scope: signed-off S26 "integrator/workflow/gothenburg" —
  supersedes S3 "KEPT plugins/Aria" and S9-kept `workflow/Cing`, both now
  fully orphaned (verified this stage).
- Verification before removal (repo-wide, outside the target trees):
  - `integrator.*`: ZERO importers — the only external hit is a COMMENT
    (`converters/PalesFormat.py:146` mentions "workflow" in prose).
    Note: `ccp/format/cosmos/*` + `converters/CosmosFormat.py` (registry-LIVE
    "Cosmos" format) do NOT import the integrator — own IO stack — kept.
  - `workflow.*`: sole importer repo-wide =
    `nijmegen/CASD/convertCasdNmrToCcpn.py:13
    from ccpnmr.workflow.Fc import FcWorkFlow` — an S29 removal target whose
    module is already broken (imports deleted `pdbe.deposition`; BY-DESIGN
    in import_smoke).
  - `gothenburg`: sole referrer = a COMMENTED-OUT `readMdd` stub in live
    `popups/EditExperiment.py` (string literal, never executed) whose
    `from gothenburg import Usf3Io` references a symbol that does not exist
    in the package (`gothenburg/` held only `__init__.py` + `mdd/__init__.py`)
    → broken by construction.
  - Build/config sweep: `pyproject.toml` hit (2 lines), `setup.py` /
    `MANIFEST.in` / `import_smoke.py` / `gui_boot_test.py`: no refs.
    Non-py data (xml/html/rst): no refs.
- **Scope extension (flagged, within signed-off "fully-orphaned modules"
  rule):** the audit's orphan list named only `workflow.{Aria,Cing}` —
  `workflow.{Constants,Fc,Util,__init__}` were NOT individually flagged
  (their sole referrer, CASD, is itself an S29 removal). Removed TOGETHER
  with Aria/Cing this stage: their entire reference graph is dead code, and
  leaving them would strand a package whose every importer is scheduled for
  deletion. Kept the removal inside this stage rather than a 4th pass in S29.
- **Removed (48 git-tracked files):**
  - `ccpnmr/integrator/` — ENTIRE tree, 40 files (core 9: Io,
    NmrpipeTableFormat, ParameterEditor, TabularFormat, Util,
    jsonToNmrCalc, projectToJson, __init__×2; plugins 31: Aria 4, Asdp 4,
    Cosmos 3, Isd 2, MultiStruc 3, NmrStar 2, Rosetta 4, Talos 4, Unio 4,
    plugins/__init__).
  - `ccpnmr/workflow/` — 6 files (Aria, Cing, Constants, Fc, Util, __init__).
  - `gothenburg/` — 2 files (`__init__.py`, `mdd/__init__.py`).
- **Edited:**
  - `popups/EditExperiment.py` (3213→3201, −12): removed the commented-out
    `readMdd` stub (its `from gothenburg import Usf3Io` + `Usf3Io.
    readDataSource` call). Adjacent commented-out `editExperimentTypes`
    stub KEPT (references live `guiParent.editExpType`).
  - `pyproject.toml` (173→171, −2): `"gothenburg*"` (packages.find) +
    `"gothenburg"` (isort known-first-party). (Note: stale `nijmegen*`/
    `utrecht*`/`molsim*` include entries from S1-12 removals remain —
    harmless: `packages.find` matches nothing; flagged for final cleanup.)
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1188** (1236−48, exactly the removed
    modules incl. `__init__.py`s), OK **1176**, FAILED **2** unchanged
    (cherrypy), BY-DESIGN **10** unchanged (CASD/ISD die with S29).
  - `gui_boot_test.py` **3/3** — `ccpnmr` app boots through the live popup
    path including the edited `EditExperiment`.
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped** (unchanged).
  - `uvx ruff check popups/EditExperiment.py --statistics`: 56 PRE-EXISTING
    (UP031 48 / F841 6 / W293 2) — the removed block was docstring interior,
    which ruff does not lint → zero NEW. `python -m py_compile` OK.
  - Residual grep: zero refs to any S26 name outside (a) `dist/` build
    snapshot, (b) `ccpnmr/analysis/doc/Readme.txt:86` historic directory
    listing (still lists extendNmr/gottingen/paris — doc history, precedent
    S22-23), (c) the S29-target CASD module (dies next-but-one), (d)
    untracked `ccpnmr.egg-info/` build metadata.

**Stage 27 — Dead popups (8) + wrappers/ (4) + frames (17) + root helpers (2) + AtomCoordList chain — ✅ 2026-08-24**
- Scope: signed-off S27 "dead popups + wrappers/Shiftx + frames + root
  helpers" — executed with 3 in-stage refinements (all verified, all
  flagged):
  1. Audit's "frames 21" resolved to **17 orphan modules** (the other 10 of
     the 28-file package are live: ExperimentList, KeysymList,
     NoeDistParamsFrame, NoeMatrix, PeakListList, PeakTableFrame,
     ResonanceFrame, SpectrumList, ViewResidueFrame, WindowFrame — kept).
  2. `ccpnmr.{Common,Constants}` (S24 memory said "zero importers") are in
     fact **LIVE** — `ccpnmr/v2io/NefIo.py:67-68
     from .. import Common as commonUtil / from .. import Constants as
     genConstants` → KEPT. Only `SafeFilename` + `_serverCheck` actually
     have zero importers (audit-confirmed; grep-confirmed — the only
     "SafeFilename" hit elsewhere is an unrelated same-named function
     `nef/SafeOpen.getSafeFilename`).
  3. **Scope extension (S26 precedent, flagged):** `wrappers/{CamCoil,D2D}`
     — audit did NOT orphan them (their two importers
     `popups/{SequenceShiftPredict,SecStructurePredict}` are S27 removals,
     so they escaped the zero-ref test). Re-verified this stage: their only
     other refs are COMMENTED-OUT stubs in `AnalysisPopup.py` (L241/1232/
     1404/2151) and two DATA strings (`ccp/general/ChemCompOverview.py:28019
     "D2D"` chem-comp code; `ccpnmr/v2io/Constants.py:15643 "D2D"` isotope
     map) — not imports. Dead by same rule as S26 `workflow.Fc` → included.
- Verification sweeps (outside target files): the 8 popups — ZERO refs
  (no AnalysisPopup import, no `popups/__init__` re-export, no rst);
  the 17 frames — zero importers, zero cross-imports from the 10 kept
  frames (`FontList` name hits are the SEPARATE `memops.gui.FontList`);
  AtomCoordList — sole importer was `EditResStructures.py:214` (removed
  this stage).
- **Removed (40 git-tracked + 2 untracked .so):**
  - popups 8: BrowseStrucGen, CalcCouplings, EditNoeClasses,
    EditResStructures, LinkSideChains, PredictKarplus,
    SecStructurePredict, SequenceShiftPredict.
  - frames 17: AxisLabelList, AxisTypeList, AxisUnitList, ChainList,
    ColorList, ColorSchemeList, ExptSpectrumPeakList, ExptSpectrumRows,
    FontList, MeasurementTypeList, MolSystemList, PanelTypeList,
    ReferenceFrame, SetupStructureCalcFrame, SpectrumViewList, SymbolList,
    ViewIsotopomerFrame (S22's importerless-kept flag now retired).
  - `wrappers/` ENTIRE package (Shiftx, CamCoil, D2D, `__init__` — empty).
  - root helpers 2: `ccpnmr/SafeFilename.py`, `ccpnmr/_serverCheck.py`.
  - **AtomCoordList chain (S25 flag "decide at S27", decided — remove):**
    its last importer (EditResStructures) is gone →
    `c/ccpnmr/clouds/{atom_coord,py_atom_coord,atom_coord_list,
    py_atom_coord_list}.{c,h}` (8) + `c/ccpnmr/clouds/AtomCoordList.so` +
    the 2 untracked python bridges. `c/ccpnmr/clouds/` is now exactly the
    Bacus surface (`{bacus,py_bacus}.{c,h}`, `Bacus.so`, `_licenseInfo.py`).
- **Edited:** `setup.py` FAM (drop AtomCoordList — clouds section now
  Bacus-only), `linkSharedObjs` (−1), `copySharedObjs` (−1),
  `copySharedObjs.bat` (−1), `memops/c/copySharedObjsMac` (−1) — all now
  `= 7 analysis lines + clouds/Bacus`; `c/ccpnmr/clouds/Makefile` →
  `all: Bacus` only.
- KEPT (deliberate): the "Shiftx" FORMAT (`ccp/format/shiftx/*` +
  `Constants.py:197 "shiftx": "Shiftx"` registry entry) is independent of
  the removed wrapper — registry-live, untouched.
  `ccpnmr.{Common,Constants}` (NefIo importers — live, see refinement 2).
  `AnalysisPopup.py` commented CamCoil/D2D stubs (comments only — left
  per "don't edit what isn't code" minimalism; revert-unit safe).
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1157** (1188−31, exactly the removed
    `.py` modules), OK **1145**, FAILED **2** unchanged (cherrypy),
    BY-DESIGN **10** unchanged.
  - `gui_boot_test.py` **3/3**.
  - `pytest ccpnmr2.5/python/tests/` **45 passed, 4 skipped**.
  - `uvx ruff check setup.py`: 1 PRE-EXISTING I001, zero NEW; `py_compile`
    OK. No other `.py` modified this stage → lint state otherwise
    structurally unchanged.
  - Residual sweep: zero hits for any removed name except
    `EditNoeClassesFrame` (an unrelated SAME-NAMED class defined inside
    KEPT `NoeDistParamsFrame.py`) and the two "D2D" data strings noted
    above.

**Stage 28 — Examples (47 py) + NEF import feature (menu + nef/ + v2io/ + tests) + root-helper collateral — ✅ 2026-08-24**
- Scope: signed-off S28 "examples + nef.testing" — plus, **decided by the
  user at stage time (2026-08-24)**: REMOVE (not hold) the "Load Nef" menu
  entry and the whole NEF (CCPNmr v2 project file) import feature, since
  FormatConverter already covers data import/export in mainstream formats
  (NMR-STAR, NMRDraw, SPARKY, PDB…). Two in-stage refinements (S26/S27
  collateral-orphan precedent, flagged): (1) the `ccpnmr/nef/` package and
  `ccpnmr/v2io/` package die WHOLE — the signed scope named only
  `nef/testing` + `v2io/{TestNefIo,Constants}`, but `NefIo` (the live
  importer of `v2io.Constants`) and the `nef/` parser stack
  (StarIo/NefImporter/…) have their entire remaining reference graph inside
  this feature (menu method + 2 NEF tests); (2) root helpers
  `ccpnmr/Common.py` + `ccpnmr/Constants.py` — S27 KEPT them solely because
  their only importer (verified this stage: `NefIo.py:67-68
  from .. import Common/Constants`) was this feature — once NefIo dies they
  are exactly the zero-importer state of the S27-removed `SafeFilename`/
  `_serverCheck`, so they follow them out.
- Verification sweeps (repo-wide, outside targets): all three example/test
  trees — ZERO importers. `NefIo` — only importers are the removed
  `AnalysisPopup.loadNefFile` (lazy, L1901) + the 2 removed NEF tests.
  `v2io.Constants` — sole importer `NefIo.py:70 from . import Constants`.
  `StarIo`/`NefImporter` — importers: NefIo, TestNefIo, the 2 NEF tests.
  `ccpnmr.Common`/`ccpnmr.Constants` — comprehensive py+xml+txt sweep:
  ONLY NefIo. `SafeOpen`/`getSafeFilename` (nef/) — zero kept users
  (S27 trap re-confirmed). `decorator` (py package) — sole importer
  `TestNefIo.py:39`. NEF in the format registry
  (`allFormatsDict`)/converters — NONE (the NEF import is a PROJECT import,
  independent of the live FormatConverter surface). `FileType`/
  `FileSelectPopup` (used by loadNefFile) — 6 OTHER live uses in the same
  file (importCoordinates/importPdb/…) → imports kept. KEPT deliberately:
  `ccp/util/V2Upgrade.py` (zero importer, BY-DESIGN `ccpncore` — signed
  S30 `ccp.util` scope, NOT touched); `survey.md` NEF lines + 2 `data/`
  XML comments + `ccp/format/nmrView/projectIO.py:286` comment (all inert
  text/history — S27 precedent); `Tests` elsewhere (test_c_ext_imports,
  test_memops_*, test_project_lifecycle) — zero NEF refs.
- Removed (123 git-tracked files, 69 `.py`):
  - `ccp/examples/` — 68 (41 `.py` + 27 data; workshops/help_doc — audit
    count "35" was strong-orphans-of-41).
  - `ccpnmr/format/examples/` — 14 (6 `.py` + 8 data).
  - `ccpnmr/nef/` — 33 (14 `.py` incl. `testing/` 4 nose-era + parser
    stack 10; 17 `testdata/*.nef` + README + `.gitignore`).
  - `ccpnmr/v2io/` — 4 (NefIo, TestNefIo, Constants 36 312 ln, `__init__`).
  - `ccpnmr/Common.py` (876) + `ccpnmr/Constants.py` (595) — NEF collateral.
  - `tests/test_nef_import.py` (9 tests, P4-2 feature test) +
    `tests/test_nef_parse.py` (11 tests).
- Edited:
  - `AnalysisPopup.py` (2895→2845, −50): removed the "Load Nef"
    `menu.add_command` block (Project menu, between "Open Spectra" and
    "Save"), the `"Load Nef"` `menu_items[ProjectMenu]` entry (15→14), the
    `fixedActiveMenus` loop `(0,1,2,3,8)` → `(0,1,2,7)` (item 8 was Quit —
    `setMenuState` maps these onto real Tk entry indices via
    `entryconfig(n, …)`, so the Quit index shifts 8→7; New/Open Project/
    Open Spectra/Quit stay active without a project), the `loadNefFile`
    method (lazy `NefIo.loadNefFile` + FileSelectPopup *.nef), and the two
    dead COMMENTED NEF stub groups (4-line imports/exports menu stub
    "# NOTE:ED not needed"; 6-line `importNefFile`/`exportNefFile` py2
    comment block).
  - `pyproject.toml` (173→170, −3): `optional` extras — dropped
    `decorator>=5.0` + its "(kept v2io test)" comment clause (sole user
    TestNefIo is gone); `testpaths` — dropped the
    `ccpnmr2.5/python/ccpnmr/nef/testing` line (would be a dangling dir on
    a bare `pytest`); `python_files` — `["test_*.py","Test_*.py"]` →
    `["test_*.py"]` (the only tracked `Test_*.py` files were the two
    removed `nef/testing` ones).
  - `scripts/linux_release.sh` + `scripts/macos_release.sh` (−1 each):
    optional-stack pip line `matplotlib cherrypy decorator mako` →
    `matplotlib cherrypy mako` (S12 had kept `decorator` specifically for
    TestNefIo). `bash -n` clean.
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1088** (1157−69, exactly the removed
    `.py` modules), OK **1076**, FAILED **2** unchanged (2× `cherrypy`,
    pre-existing), BY-DESIGN **10** unchanged (V2Upgrade stays — signed
    S30 scope).
  - `gui_boot_test.py` **3/3** (the CCPN main window boots through the
    edited `setProjectMenu` path; data-shifter + format-converter — the
    user-cited NEF replacement — untouched).
  - `pytest ccpnmr2.5/python/tests/` **25 passed, 4 skipped** (45−20,
    exactly the 9+11 removed NEF tests; the 4 skips are the non-NEF ones).
  - `uvx ruff check AnalysisPopup.py --statistics`: **38 — IDENTICAL mix**
    to post-S23 (UP031 17 / E722 10 / F841 6 / E731 2 / E721 1 / F811 1 /
    W293 1) — zero NEW (the removed method's `except Exception as es:` was
    not a flagged pattern). `python -m py_compile` OK; `bash -n` both
    release scripts OK.
  - Residual sweep (py/toml/sh/in/cfg, excl. dist/build/egg-info/pycache):
    exactly ONE hit — the inert comment
    `ccp/format/nmrView/projectIO.py:286`
    (`#'python/ccpnmr/format/examples/data/…str'`); plus the known
    `survey.md` history lines. Zero functional refs anywhere.

**Stage 29 — Dead CASD/ISD orphan cluster (nijmegen/ whole + adatah pair + cambridge/isd whole) + allowlist shrink 10→5 — ✅ 2026-08-24**
- Scope: signed-off S29 "CASD/ISD (+allowlist shrink)" — two whole dead
  packages + one orphan module pair. Pre-stage verification:
  `nijmegen/` contained ONLY `__init__.py` (empty) + `CASD/` → the package
  dies whole, exactly the signed `nijmegen/CASD` scope. `cambridge/isd`:
  every cross-file reference (IsdFrame↔NmrCalcExchange↔CCPNReader↔
  isd_project_template) is inside the dir; zero external importers repo-wide
  (`ccpnmr/analysis` has NO `isd` references at all — no menu wiring, no
  `popupActions` entry; the `testIsdPopup(argServer)` macro-era entry point
  is referenced by nothing). `cambridge/__init__.py` is just `pass`.
  `pdbe.adatah.{CasdNmr,Pdb}`: sole references are TWO already-dead
  comments inside S29-itself-removed
  `nijmegen/CASD/convertCasdNmrToCcpn.py` (L211 `#from pdbe.adatah.Pdb …`,
  L997 `#from pdbe.adatah.CasdNmr …`) — no live importer anywhere else.
  `casdPipeLine.py`'s own imports prove orphan status:
  `pdbe.deposition.dataFileImport.formatConverterWrapper` (removed in an
  earlier stage) + `ccpnmr.workflow.Fc` (S26) + `nijmegen.CASD.Constants`
  (raises at import if `CASD_HOME` unset — the original BY-DESIGN reason).
  Zero tests reference any target. `memops.scripts` / `scripts/` /
  `MANIFEST.in` / `setup.py`: no refs.
- Removed (15 git-tracked files, 14 `.py` + 1 README):
  - `ccpnmr2.5/python/nijmegen/` — ENTIRE package, 6: empty `__init__.py`
    + `CASD/` 5 (`__init__`, `Constants` 46 — `raise Exception("Environment
    variable CASD_HOME not set")` at import, `Util` 625, `casdPipeLine`
    ~120, `convertCasdNmrToCcpn` 660).
  - `ccpnmr2.5/python/cambridge/isd/` — ENTIRE dir, 7: `__init__`,
    `isd_project_template` (import-time `ISD_ROOT` check — the other
    BY-DESIGN entry), `NmrCalcExchange`, `IsdPopup`, `IsdFrame` (1200+),
    `CCPNReader` (1200+), `README`.
  - `ccpnmr2.5/python/pdbe/adatah/{CasdNmr,Pdb}.py` — CASD's orphan pair
    (rest of `adatah` stays — signed S30 scope).
- Edited:
  - `import_smoke.py`: `KNOWN_NON_IMPORTABLE` 10 → **5** — deleted the
    whole former ENV category (`nijmegen.CASD.{Constants,Util,casdPipeLine,
    convertCasdNmrToCcpn}` + `cambridge.isd.isd_project_template`); the 5
    remaining entries are all EXTERNAL (PyMC2/sans/ccpncore/pdbe-analysis/
    memops.scripts — S30 re-verify still pending). Comment block updated to
    note the ENV category died here.
  - `pyproject.toml` (−2 lines): `nijmegen*` out of `packages.find.include`
    + `nijmegen` out of `isort.known-first-party` — the package no longer
    exists (S24 memory note: "S29's CASD removal can take `nijmegen*`" —
    taken). `cambridge*` stays (bayes/… survive).
- Gates (all green):
  - `import_smoke.py` exit 0 — TOTAL **1074** (1088−14: 5 `nijmegen.*` +
    package + 5 `cambridge.isd.*` + package − dedup = exactly the 14 removed
    modules; measured), OK **1067** (1074−2−5), FAILED **2** unchanged
    (2× `cherrypy`, pre-existing), BY-DESIGN **5** (all EXTERNAL).
  - `gui_boot_test.py` **3/3** (ccpnmr / data-shifter / format-converter).
  - `pytest ccpnmr2.5/python/tests/` **25 passed, 4 skipped** — unchanged
    (no CASD/ISD tests existed).
  - `uvx ruff check AnalysisPopup.py --statistics`: **38 — IDENTICAL mix**
    (UP031 17 / E722 10 / F841 6 / E731 2 / E721 1 / F811 1 / W293 1).
    `python -m py_compile import_smoke.py` OK.
  - Residual sweep (py/toml/sh/in/cfg, excl. dist/build/egg-info/pycache):
    5 inert hits — `pdbe/nmrStar/IO/nmrStarDict.py:231` (`CASD$` NMR-STAR
    data-dict string), `ccp/format/molmol/{sequenceIO:141,
    coordinatesIO:266}` (`markNijmegen` PDB file paths — a PERSON's name,
    same established trap as S27's "Shiftx"), `ccp/util/NmrCalc.py:4`
    (docstring mention), `ccp/general/ChemCompOverview.py:36500` ("Isd"
    chem-comp code); plus known `README.md:163` ("some CASD/education
    scripts … kept in" — now stale text, same precedent as the `survey.md`
    history lines). Zero functional refs anywhere.
- Kept deliberately (signed S30 scope, untouched): the rest of
  `pdbe/adatah` (Io/Util/Bmrb/Generic/… — `Util` has live importer
  `ccpnmr/format/process/sequenceCompare.py:4`), `ccp/util/V2Upgrade`
  (BY-DESIGN EXTERNAL).

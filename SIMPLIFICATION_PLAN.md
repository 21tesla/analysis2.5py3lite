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

Status: **IN PROGRESS** (started 2026-08-24)
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
| 20 | Assignment ▶ Automated Seq. Assignment + `ccpnmr/nexus/` (5) + `wrappers/{Mars,Psipred}.py` + doc lines | ⬜ |
| 21 | Assignment ▶ NOE Contributions + `popups/LinkNoeResonances.py` + doc line + separator collapse | ⬜ |

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

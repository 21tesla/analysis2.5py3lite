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
| 10 | ECI: `ccpnmr/eci/` (relocate `ReadPdb.py` first) + `ccpnmr-eci` | ⬜ pending |
| 11 | Structure Viewer + Make H Bond Restraints popups + remove 3 kept callers | ⬜ pending |
| 12 | Cross-cutting sweep: `pyproject.toml`, `bin/`, release scripts, docs, extras, final verification | ⬜ pending |

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

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
| 3 | ARIA: `paris/` + menu + methods | ⬜ pending |
| 4 | CYANA: `cyana2ccpn/` + `macros/MultiStructure.py` + integrator `Cyana/` + `Io.py` import | ⬜ pending |
| 5 | DANGLE: `cambridge/dangle/` + `ccpnmr-dangle` | ⬜ pending |
| 6 | HADDOCK: `utrecht/` | ⬜ pending |
| 7 | MECCANO: `grenoble/meccano/` + C sources + `setup.py` GSL block | ⬜ pending |
| 8 | PyRPF: `rutgers/` | ⬜ pending |
| 9 | CING: `cing/` + `nijmegen/cing/` + smoke allowlist | ⬜ pending |
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

# NEF Integration Plan — adopt CCPNMR v3 NEF + project code into py3lite

**Status: STAGE 37 DONE** (export round-trip green on all 3 testdata; Stage 38 next = app wiring)

## Context (user directive 2026-08-25)

> "The directory `../ccpnmr` contains version 3 of this software that was written in
> Python 3.11 and PyQt. Inspect the code for the manner in which it reads and writes
> NEF files. Also how it imports a project of same type as this present software and how
> it converts it into the mainstream format that version 3 uses. The goal of this exercise
> is to adopt and integrate the ccpnmr version 3 NEF and project code into this optimized
> legacy code."

Source (v3): `/home/logan/software/ccpnmr` (dev snapshot, NOT a git repo, read-only for us).
Target: this repo (`analysis2.5py3lite`), code root `ccpnmr2.5/python/`.

Background: the legacy repo previously shipped a v2-era NEF stack (`ccpnmr/nef/`,
`ccpnmr/v2io/` — removed in Stage 28 as part of the simplification). This effort adopts
v3's maintained NEF stack back into the legacy, as self-contained modules (house pattern,
S33 `exportNmrData` style): model-free core → model binding → app wiring.

## Architecture map (recon 2026-08-25, verified in source)

### v3 NEF stack — three layers

1. **Model-free format core** — `src/python/ccpn/util/nef/`:
   - `StarTokeniser.py` (151) — STAR tokenizer; deps: stdlib only.
   - `GenericStarParser.py` (1067) — Star → object tree (DataExtent/DataBlock/SaveFrame/Loop/LoopRow, ordered dicts); deps: stdlib only.
   - `StarIo.py` (622) — NEF-restricted layer: `NmrDataBlock`/`NmrSaveFrame`/`NmrLoop`/`NmrLoopRow`, `_StarDataConverter` (numeric heuristics), `parseNef`/`parseNefFile`, `toString`; deps: stdlib only.
   - `ErrorLog.py` (218) — error-code framework (`NEFERROR_*`, `NEF_STANDARD` modes); deps: stdlib only.
   - `Specification.py` (396) — `CifDicConverter`: `mmcif_nef_v1_1.dic` → validation dict; deps: stdlib only.
   - `Validator.py` (283) — `isValid(dataBlock, specDict)`.
   - `NefImporter.py` (1473) — `NefImporter` (loadFile/saveFile/toString/fromString,
     get*/add* per category, saveframe management, optional reader/verifier/content
     attachment) + `NefDict` (saveframe API; optional pandas in `getTable(asPandas=True)`);
     deps: stdlib + **numpy** (one import; make optional).
   - `SafeOpen.py` (115) — retry-open helper. `nef.py` (1535)/`CompareNef.py` (1033) — file
     compare CLI (NOT needed here). `mmcif_nef_v1_1.dic` (141KB) — the format dictionary.
   - `testdata/` — 16 `.nef` samples.
   - NEF v1.1 (BMRB NMR Enhanced Format, `Nmr_Exchange_Format`) categories
     (`NefImporter.NEF_CATEGORIES`): `nef_nmr_meta_data` (mandatory),
     `nef_molecular_system` (mandatory; `nef_sequence` loop), `nef_chemical_shift_list`,
     `nef_distance_restraint_list`, `nef_dihedral_restraint_list`, `nef_rdc_restraint_list`,
     `nef_nmr_spectrum` (spectrum + `ccpn_peak_list`), `nef_peak_restraint_links`.
   - **NEF files carry metadata/peaks/shifts/restraints — no raw matrix data.**

2. **v3 model binding** — `src/python/ccpn/framework/lib/ccpnNef/` (bound to
   `ccpn.core.*` wrapper model): `CcpnNefReader` (`importNewProject`, ~60
   `load_nef_*`/`load_ccpn_*` loaders, `produceNmrChain/Residue/Atom`), `CcpnNefWriter`
   (`exportProject`; `saveNefProject`/`exportNef`/`convertToDataBlock`/`writeDataBlock`;
   `makeNefMetaData`/`chains2Nef`/`chemicalShiftList2Nef`/`restraintTable2Nef`/
   `peakList2Nef`/`spectrumGroup2Nef` …), `CcpnNefCommon` (`nef2CcpnMap`,
   `saveFrameReadingOrder`), `CcpnNefContent`, `CcpnNefImporter` (thin glue used by
   `NefDataLoader` → `Framework._loadNefFile`).

3. **Legacy-model binding (v2, data model 2.1.2)** — `src/python/ccpn/util/v2io/NefIo.py`
   (1996) — v3 kept it "GWV 02 Feb 2022: Likely obsolete code" for compat. **It is already
   written against THE LEGACY MODEL**: `loadNefFile(path, memopsRoot=None)`,
   `loadProject(nefFilePath, pdbFilePaths=None, projectName=None, pdbFileType='pdb')`,
   `CcpnNefReader` (legacy), `createMoleculeFromNef`, `extendMolResidues`,
   `makeNefAxisCodes`, `addDataStore(dataSource, spectrumPath, **params)`,
   `fetchDataUrl(memopsRoot, fullPath)`; + `TestNefIo.py` (340; 45 tests gated on
   external data dirs `/Users/ejb66/...` — not shippable, adapt to bundled testdata).
   v3 also kept the root-helper pair `ccpnmr/Common.py`+`Constants.py` that this repo
   removed in S28 (v3 copy: `ccpn/util/v2io/Constants.py` 36k gen. — not needed here).

### v3 project layer (the "mainstream format")

- v3 mainstream project = **memops-based XML** (same engine family as legacy; v3
  `ccpnmodel/` ships versioned upgraders `v_2_0_4 … v_2_1_2, v_3_0_2, v_3_0_a1, v_3_1_0`)
  under a new wrapper-object runtime (`ccpn.core.*`: `AbstractWrapperObject`, `Pid`,
  `Updater`), with HDF5 as an additional data store (optional; NEF itself has none).
- v2 import path: `CcpNmrV2ProjectDataLoader` (requires `memops/<Implementation>` subdir)
  → `Framework._loadV2Project` → `ccpn.core.Project._loadProject`
  (`ccpn/core/Project.py:3395`) → `XmlLoader(path, readOnly=True)` → `xmlLoader.isV2` →
  copy project to new `<basename>.ccpn` path via `XmlLoader.newFromLoader` →
  `updateProject_fromV2(project)` (`ccpn/core/_implementation/updates/update_v2.py`, 61L)
  → `project.save()`.
- File-type dispatch: `DataLoaderABC._registerFormat()` registry (`NefDataLoader`
  suffixes `['.nef']` → `Framework._loadNefFile` → `CcpnNefReader.importExistingProject`).

### Legacy baseline (verified)

- Model `ccpnmr2.5/python/ccp/api/nmr/Nmr.py` (generated, 150kL): `NmrProject` (:43024),
  `DataSource` (:14469 — the "Spectrum": `dataStore`/`dataType`/`numDim`/`dataDims`/
  `recordNumber`/`peakLists`), `Peak` (:48253), `PeakList` (:61830), `Resonance`
  (:67208), `NmrConstraint`.
- Project IO: `memops.general.Io.loadProject` (+ S34-restored version-compat island for
  pre-2.1.2 files); dataStore wiring `ccp/general/Io.py:717`
  (`getDataSourceFileName`/`setDataStoreFileName`).
- Zero NEF support today (grep-verified; S28 removed the v2-era stack).
- House pattern: self-contained module + tests — `ccpnmr/exportNmrData.py` +
  `tests/test_exportNmrData.py` (S33).
- Gates: `import_smoke.py` (baseline S34: TOTAL 910 / OK 909 / FAILED 0 / BY-DESIGN 1),
  `pytest` (39 passed, 4 skipped), `uvx ruff check` (38 + S34 set, zero-new policy),
  `gui_boot_test.py` (3/3).

## Locked decisions

1. **Adopt v3's modern core** (`util/nef`), not the v2-era parser — same feature,
   maintained, model-free. Restores (with replacement) what S28 removed.
2. **NEF scope = BMRB NEF v1.1 (Nmr_Exchange_Format)** — the 8 categories above.
   Raw spectrum matrix data is out of format (NEF = metadata/peaks/shifts/restraints).
3. **Model binding for import = `ccpn/util/v2io/NefIo.py`** (v3's copy of the v2 NEF
   project-importer) adapted to this tree — it is already legacy-model code
   (`memopsRoot`, `DataSource`, `addDataStore`).
4. **Model binding for export = new module** using the Stage-35 core, using v3's
   `CcpnNefWriter` methods as the semantic reference (they are v3-model code and cannot
   be ported directly).
5. **`updateProject_fromV2` / v3 project-directory opening OUT OF SCOPE** — it is v3
   wrapper-runtime code (`ccpn.core.*`, `Pid`, `Updater`, new `XmlLoader`) and porting it
   would rewrite the legacy runtime. NEF IS v3's project interchange format: importing
   v3 project *content* into the legacy happens via NEF import (St. 36); exporting legacy
   projects as NEF (St. 37) is the reverse. Re-scope if the user wants `.ccpn` dir opens.
6. **No new hard deps.** numpy import in `NefImporter` made optional; pandas stays
   optional (only `NefDict.getTable(asPandas=True)`).
7. **Packaging:** NEF core + dictionary + bundled testdata under `ccpnmr2.5/python/ccpnmr/…`
   (setuptools discovery already covers `ccpnmr*`); dictionary + testdata as data files.
8. **Checkpoint policy (user directive):** every stage ends with commit + push
   (code + plan status/log in ONE commit) + memory update, then the next stage starts.

## Stages

**Stage 35 — NEF format core (model-free) + tests.**
Port `ErrorLog, StarTokeniser, GenericStarParser, StarIo, Specification, Validator,
NefImporter, SafeOpen` → `ccpnmr2.5/python/ccpnmr/nef/` (py3-compat fixes; optional
numpy; `NEF_ROOT_PATH`-anchored dictionary load:
`ccpnmr/nef/mmcif_nef_v1_1.dic`). Bundle 3 small testdata files
(`CCPN_Commented_Example.nef`, `CCPN_XPLOR_test1.nef`, `CCPN_Sec5Part3.nef`).
New `tests/test_nef_core.py`: parse all 3 (validate + `getNmrSpectra()` + categories),
`toString`→`fromString` round-trip, `saveFile`→`loadFile` equality on one file.
Accept: import_smoke FAILED 0; pytest green; ruff zero-new.

**Stage 36 — NEF import → legacy model (the "NEF→MOPS" path).**
Adapt v3 `ccpn/util/v2io/NefIo.py` → `ccpnmr2.5/python/ccpnmr/v2io/NefIo.py`
(`loadNefFile`, `loadProject(nefFilePath, pdbFilePaths, projectName, pdbFileType)`,
legacy `CcpnNefReader`, `addDataStore`/`fetchDataUrl`, `createMoleculeFromNef`,
`extendMolResidues`). Re-point imports at Stage-35 core; drop `decorator`; headless-safe
MessageReporter shims; py3 fixes. Tests `tests/test_v2io_nef.py` on bundled testdata
(`loadNefFile` into a fresh memops root: DataSources/peaks/shifts present;
`loadProject` with a PDB: project dir + MolSystem + DataSource).
Accept: import_smoke FAILED 0; pytest green; ruff zero-new.

**Stage 37 — NEF export ← legacy model + round-trip.**
New `ccpnmr2.5/python/ccpnmr/nefExport.py`: NmrProject → NmrDataBlock
(metadata + molecular_system + chemical_shift_list + restraint lists + nmr_spectra with
peak lists) → `saveFile`. Port the v3 `CcpnNefWriter` semantics (makeNefMetaData,
chains2Nef, chemicalShiftList2Nef, restraintTable2Nef, peakList2Nef) against legacy
classes (`NmrProject`/`Resonance`/`NmrConstraint`/`DataSource`/`PeakList`). Test
round-trip: import(export(project)) preserves spectra/peaks/shifts (tolerant compare).
Accept: import_smoke FAILED 0; pytest green; ruff zero-new.

**Stage 38 — app wiring + docs + final gates.**
GUI: "Load NEF…" (Project menu) via `ccp.gui.Io` → `NefIo.loadProject`; "Export NEF…"
popup (Project menu) → `nefExport`. CLI: console entry (check pyproject entry-point
style) for NEF import/export. Docs: NEF support section in README/INSTALL. Close-out:
plan marked COMPLETE + final gates (import_smoke / pytest / ruff / gui_boot 3/3) +
push. Accept: all gates ≥ baseline, zero new ruff.

## Stage log (append per stage — commit + push + memory each)

**Stage 35 — NEF format core (model-free) + tests — ✅ 2026-08-25**
- Ported 8 v3 modules from `ccpn/util/nef/` → `ccpnmr2.5/python/ccpnmr/nef/`:
  `__init__` (`NEF_ROOT_PATH`), `ErrorLog`, `StarTokeniser`, `GenericStarParser`,
  `StarIo`, `Specification`, `Validator`, `NefImporter`, `SafeOpen`. Package name
  `ccpnmr.nef` restores the v2-era import path the v2io tests use.
- Functional edits (only these): (1) `NefImporter` — top-level `import numpy as np`
  removed → lazy inside `_convertToPandas` try-block (sole use site; numpy+pandas both
  optional now); (2) `GenericStarParser` — 3 invalid `\s` escapes → raw strings
  (SyntaxWarnings on py3.13); (3) verified ZERO `ccpn.*` external refs in all 8.
- Bundled with package: `mmcif_nef_v1_1.dic` (141KB, validator dictionary) + 3 small
  testdata files (`CCPN_Commented_Example`, `CCPN_XPLOR_test1`, `CCPN_Sec5Part3`;
  the 2 H1GI 7-10MB files + docr set deliberately NOT bundled — dist weight).
- NOT ported: `nef.py`/`CompareNef.py` (file-compare CLI, not needed), `testing/`
  (nose-era), `NEF/` (license/charter docs).
- Tests: new `ccpnmr2.5/python/tests/test_nef_core.py` — **12 tests**: parse all 3
  bundled files (mandatory meta+ms frames + samefile path), ground truths (Commented:
  104 shifts / 2 spectra [cnoesy1+dummy15d] / 235 seq residues / 9 validator errors;
  XPLOR: no spectra (getNmrSpectra is None) / 735 dist-restraints; Sec5: 5 spectra /
  `isValid True` with 5 informational notes), 8 categories listed, saveframe mgmt
  (add/get/has/delete, prefixed forms), toString→fromString + saveFile→loadFile
  round-trips, direct StarIo parse (`nef_peak` = the NEF v1.1 peak loop). API learns:
  StarIo coerces `'1.1'`→float 1.1; `getTableNames()` hides `nef_` prefixes; loop rows
  via `NmrLoop.data`.
- Gates: import_smoke TOTAL **920** (=910 + 9 nef modules + 1 test file; the smoke
  walk includes tests/) OK 919 / FAILED **0** / BY-DESIGN 1 (unchanged, PyMC bayes);
  pytest **51 passed, 4 skipped** (39+4 baseline + 12 new, all green); ruff (0.16.3)
  new package: 174 findings, **all pre-existing in the v3 source** (UP031 83 —
  deferred per the user's Phase-3 "no bulk style" decision; UP032 23; I001 17; UP010
  17; UP040 6; E402 5 [v3 standalone-bootstrap block]; E713 5; F541 3; F841 3; UP008
  3; E722 2; ≤2 each of UP004/UP009/UP030/UP034) — ZERO live-file changes, zero new
  on existing files; gui_boot_test untouched (no GUI edits this stage).
- Packaging: no pyproject change needed (`ccpnmr*` include + `"*": ["*"]`
  package-data already cover `.dic`/`.nef`).

**Stage 36 — NEF import → legacy model (the "NEF→MOPS" path) — ✅ 2026-08-25**
- PROVENANCE (verified byte-identical, `diff=0` vs git history `bc4059f1^` —
  the S28-removed set, S34-style restore): `ccpnmr/v2io/NefIo.py` (2045; the
  original 2.5 "CCPN V2, data model 2.1.2" importer — the same file v3 retained
  in `ccpn/util/v2io/` for back-compat), `ccpnmr/v2io/Constants.py` (36312,
  generated chemComp tables), `ccpnmr/v2io/__init__.py` (27),
  `ccpnmr/Common.py` (876: `resetSerial`, `name2IsotopeCode`,
  `isotopeCode2Nucleus`, axis-code utils), `ccpnmr/Constants.py` (595:
  `DEFAULT_ISOTOPE_DICT`), `ccp/lib/NmrExpPrototype.py` (348:
  `refExpDimRefCodeMap`, `getAtomSiteAxisCode`, …). NOT restored: legacy
  `v2io/TestNefIo.py` (gated on `/Users/ejb66/...` external data — would pass
  vacuously; replaced by the new test below).
- NO edits needed: the legacy file is already (a) written against this
  exact model, (b) import-free of anything removed since S28/S30 except the
  three helpers above (all present again), and (c) already imports
  `..nef.StarIo` — which Stage 35 restored as `ccpnmr.nef`, so the legacy
  binding meets the new core at the original v2-era import path. Public API
  as-is: `loadNefFile(path, memopsRoot, overwriteExisting)`,
  `loadProject(nefFilePath, pdbFilePaths, projectName, pdbFileType)`,
  `CcpnNefReader`, `createMoleculeFromNef`, `extendMolResidues`,
  `makeNefAxisCodes`, `addDataStore`, `fetchDataUrl`, `assignPeak`,
  `saveFrameReadingOrder`.
- NEW tests `ccpnmr2.5/python/tests/test_v2io_nef.py` — **8 tests** (pytest,
  bundled testdata, network-independent: ChemComp download for the dummy
  residues falls back to UNK offline), ground truths probed on the live
  legacy model:
  - public API surface + `saveFrameReadingOrder` order
  - Commented: 235 residues / 15 chains; experiments
    `('15N NOESY-HSQC', 3)` + `('HNCCCCCCCCCCCCC', 15)`; 2 DataSources,
    2 PeakLists, 6 peaks all resonance-assigned (peakContribs →
    peakDimContribs); 2 ShiftLists / **108 Shifts** (104 NEF rows expand
    via wildcard resonances, e.g. HD%); 1 NmrConstraintStore with
    Dihedral(6)/Distance(3)/HBond(4)/Rdc(2)
  - XPLOR: restraints only (no spectra): Distance **735**, Dihedral 161,
    Rdc 147+152 (two lists)
  - Sec5: 5 spectra (`15N HSQC/HMQC` 2D + 4 3D), 5 PeakLists,
    **891 peaks**
  - `loadProject(nef, projectName=...)` creates the full project dir
  - **saveProject → loadProject round-trip** preserves 235 residues / 6
    peaks / 108 shifts / 4 constraint lists — exercises the legacy memops
    XML persistence path (incl. the S34 version-compat island)
- Legacy-model API notes (needed by Stage 37-38): MolSystem residues live
  under `MolSystem.getChains()` → `chain.sortedResidues()` (no direct
  residue relation); "spectrum" = `NmrProject.sortedExperiments()` →
  `Experiment.sortedDataSources()` → `DataSource.getPeakLists()` →
  `PeakList.getPeaks()`; shifts = `AssignmentBasic.getShiftLists(nmr)`
  (=`findAllMeasurementLists(className='ShiftList')`) →
  `ShiftList.getMeasurements()` (frozensets); constraints =
  `NmrProject.getNmrConstraintStores()` →
  `NmrConstraintStore.getConstraintLists()` → typed list `getConstraints()`.
- Gates: import_smoke TOTAL **928** (=920 + 6 restored modules + 1 test
  + 1 `nefExport.py` WIP-on-disk) FAILED **0** / BY-DESIGN 1 (unchanged,
  PyMC bayes); pytest **59 passed, 4 skipped** (51+4 baseline + 8 new, all
  green); ruff on the 7 in-scope files: **69 findings** (UP031 54, E402 9,
  E721 2, E722 1, UP028 1, F841 1, F601 1) — ALL in the
  history-restored/upstream-source code (verified: legacy S28 `NefIo.py`
  carries the same 40 UP031+F841; UP031 deferred per the Phase-3
  "no bulk style" decision) — ZERO findings in the newly written test,
  ZERO new findings on existing files; gui_boot_test **1/1 PASS** (the
  suite's APPS list currently has one entry, `ccpnmr`).
- Left for Stage 37 (scope: export): `ccpnmr/nefExport.py` (555-line WIP on
  disk from the aborted session, untracked, imports cleanly — NOT part of
  this commit).

**Stage 37 — NEF export ← legacy model + round-trip — ✅ 2026-08-25**
- `ccpnmr2.5/python/ccpnmr/nefExport.py` (the WIP from the aborted session,
  now complete — 623 lines, tracked from this commit):
  `makeNefDataBlock(memopsRoot)` → `ccpnmr.nef.StarIo.NmrDataBlock`
  (contemporary NEF v1.1) + `exportProject(memopsRoot, fileName)` (+ a
  small CLI under `__main__`). Saveframes in the
  `v2io.NefIo.saveFrameReadingOrder` order: meta, molecular_system
  (`nef_sequence`), one chemical-shift-list per ShiftList, restraint lists
  (distance / hbond / rdc / dihedral + `ccpn_restraint_list` for
  JCoupling / ChemShift / Csa), one spectrum per DataSource (dimensions +
  transfers + `nef_peak` rows, incl. the multi-row alternative-assignment
  form), `nef_peak_restraint_links` (only when constraints reference
  peaks). Column sets mirror a real CCPN-exported NEF v1.1 file.
- Fixed 5 bugs in the WIP (each root-caused by probing the live imported
  model before editing):
  1. `residue_name` — WIP wrote the legacy `residue.ccpCode` (title-case
     'Ala'…; 'T' for nucleic acids), which is NOT a key of
     `v2io.Constants.residueName2chemCompId` (reimport died with
     `ValueError: Unknown residueName T`). Now writes the standard
     3-letter code from
     `residue.chemCompVar.chemComp.code3Letter` — verified to reproduce
     231/235 of the original file's names exactly (the other 4 = the
     dummy-linker residues the IMPORT side already fell back to `UNK`;
     `'UNK'` is a valid map key, so they re-import fine). Resonance
     identity rows use the same code when the group is linked to a
     MolResidue, else keep the importer-assigned `rg.ccpCode`
     (preserves e.g. the 'Glx' wildcard on unassigned chains).
  2. `sequence_code` trailing space — the model default for
     `seqInsertCode` is a bare `' '`; now stripped (insert codes like
     `24B` and negative seq codes `-3` round-trip; both sides parse with
     `ccpnmr.Common.parseSequenceCode`).
  3. `%`-atom shift rows — ONE NEF shift row may back several Shift
     objects: an ambiguous atom set (e.g. `HG%`) expands to one resonance
     per atom set on import, so the row's shifts carry different name
     spellings (`HG%`, `Hg*`) of the same atom; one-row-per-shift re-
     expanded on reimport (+4 shifts on Commented). Rows are now grouped
     by (chain, seq, normalized atom, isotope, value); multi-spelling
     groups are written ONCE in the canonical upper-case '%' form,
     same-spelling duplicates (genuine duplicate file rows, e.g. the W.2
     2H ones) are kept.
  4. `element@serial` atom names — the importer deliberately creates
     `name=None` resonances for the reserved `element@serial` atom form
     (Sec5's unassigned-chains, e.g. `H@349`); the exporter now recreates
     `{element}@{resonance.serial}` or reimport crashed with
     `AtomName must be given`.
  5. One row per DIHEDRAL item — dihedral limits live on the items and
     the importer creates exactly one item per row; the WIP wrote one row
     per constraint, silently dropping items 2..N (Commented L2 has
     6 dihedrals with 11 items: 1+4+2+2+1+1 → 11 exported rows, values
     now round-trip). Pairwise restraints stay one row per constraint
     (the importer derives the item product from the row's resonances).
- Round-trip ground truths (import → export → reimport; counts, then
  value multisets — all GREEN):
  - Commented: 235 residues / 15 chains; experiments `('15N NOESY-HSQC',
    3)` + `('HNCCCCCCCCCCCCC', 15)`; 2 DataSources / 2 PeakLists; 6 peaks
    all assigned (positions equal); 2 ShiftLists / 108 shifts (values
    equal; 104 exported rows = the original file's row count, 93+11);
    constraint lists Dihedral(6) / Distance(3) / HBond(4) / Rdc(2),
    per-constraint target/limits equal.
  - XPLOR: 58 residues; no spectra / no shifts; 735 Dist / 161 Dih /
    147+152 Rdc (two lists); 1195 per-constraint values equal.
  - Sec5: 95 residues; 5 spectra; 891 peaks (positions equal, 813
    assigned); 542 shifts (values equal; exercises the
    `element@serial` rows).
- New tests `ccpnmr2.5/python/tests/test_nef_export.py` — **5 tests**:
  public API surface; exported-file structure (saveframe set + row
  counts, incl. the 93/11 shift-list split and the 11-row L2); and
  import→export→reimport round-trip count + value tests (shift-value /
  peak-position / constraint-value multisets) for all 3 bundled files —
  mirroring `test_v2io_nef.py` style (`_load` via
  `memopsIo.newProject` + `NefIo.loadNefFile` under `redirect_stdout`).
- Gates: import_smoke TOTAL **929** (=928 + 1 test file) FAILED **0** /
  BY-DESIGN 1 (unchanged, PyMC bayes); pytest **64 passed, 4 skipped**
  (59+4 baseline + 5 new, all green); ruff: `nefExport.py` **0**
  (WIP entered with 17 findings: 16×UP031 + 1×UP004 — all converted to
  f-strings / bare class), test file **0**; gui_boot **1/1** (NO GUI
  edits this stage — the menu wiring is Stage 38).
- Left for Stage 38: GUI "Load NEF…" / "Export NEF…" Project-menu items
  (via `ccp.gui.Io` → `NefIo.loadProject` / `nefExport.exportProject`),
  CLI console entry (check pyproject `[project.scripts]` style),
  README/INSTALL NEF section, final gates + plan COMPLETE + push.


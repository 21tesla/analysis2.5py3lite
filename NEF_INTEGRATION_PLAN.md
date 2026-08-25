# NEF Integration Plan — adopt CCPNMR v3 NEF + project code into py3lite

**Status: STAGE 35 IN PROGRESS** (plan + checkpoint 0 committed)

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

<!-- Stage 35 log appended here -->

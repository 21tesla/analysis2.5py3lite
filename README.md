# CCPNMR 2.5.2 (Python 3.13)

CCPNMR Analysis is a workbench for NMR (nuclear magnetic resonance) data analysis:
peak picking/fitting, assignment, restraints, structure calculation and validation —
built on the CCPN data model (MOPS) with a Tkinter GUI. This tree is the 2.5.2
release (circa 2020) modernized to run cleanly on **Python 3.13**, with the C/Cython
extension layers rebuilt and a pip-installable distribution.

## What works (verified gates)

All of the following are green on this tree (see `_phase4_checkpoints.md` for the
audit trail):

| Gate | Result |
|---|---|
| Whole-tree compile (Python 3.13) | 0 syntax errors |
| Import smoke — every module in the tree | 1646 OK / 0 failed / 83 documented-by-design |
| Functional test suite (pytest) | 43 passed / 14 skipped / 10 data-gated (need a `testdata` corpus, pre-existing) |
| C extensions | 30 setup.py exts + cing `superpose` all import and pass functional checks |
| GUI launch under Xvfb | **8/8 console apps boot** (source and installed states) |
| Distribution | sdist + wheel build; clean-venv install passes all of the above |

## The eight console commands

| Command | App |
|---|---|
| `ccpnmr` | CCPNMR Analysis (main workbench) |
| `ccpnmr-eci` | Entry Completion Interface |
| `ccpnmr-dangle` | DANGLE side-chain validation (Cambridge) |
| `ccpnmr-data-shifter` | Project data shifter |
| `ccpnmr-deposition` | PDB/Pdbe deposition data-file importer |
| `ccpnmr-extend-nmr` | EXTEND-NMR (ARIA/CING/HADDOK extenders) |
| `ccpnmr-format-converter` | CCPN project format converter |
| `ccpnmr-update` | Automatic project update (non-GUI) |

Each app's `main()` entry was added/repaired in Phase 4 (P4-1) so the entry points
actually resolve.

## Installation

### Prebuilt distribution (recommended)

```sh
uv venv --python 3.13 .venv            # or: python3.13 -m venv .venv
uv pip install --python .venv/bin/python <ccpnmr-2.5.2-...whl>
# optional stacks (see below):
uv pip install --python .venv/bin/python "scipy" "matplotlib" "sqlalchemy" \
    "cherrypy" "decorator" "mako" "psycopg2-binary" "pycurl"
```

### Build from source

```sh
uv build                                # -> dist/*.whl + dist/*.tar.gz
```

Build-time requirements (Linux):

- C compiler with GL/Tk in the sysroot (on this host: `CC=/usr/bin/gcc CXX=/usr/bin/g++`
  — Anaconda's bundled compiler lacks `GL/glx.h`)
- `freeglut` (glut), Tk/Tcl 8.6 headers + libs, `libX11.so.6`, Python 3.13 headers,
  system GL (`libGL`)
- **optional** GNU Scientific Library for Meccano (see Optional features)

`pip install <sdist>` compiles the C extensions at install time (MANIFEST.in ships the
sources). The wheel ships the compiled extensions (cp313, linux_x86_64).

### Optional features

- **Meccano** (grenoble, structure restraints) — `pip install ccpnmr[optional]` + GSL
  at build time (`conda install -c conda-forge gsl` or `apt install libgsl-dev`,
  `CCP_GSL_PREFIX=<prefix>` when building). Without it the build still succeeds
  (warning) and `grenoble.meccano.MeccanoPopup` raises an actionable import error;
  every other feature works.
- **cing / NRG / web modules** — pull the `optional` extras above (scipy, matplotlib,
  sqlalchemy, cherrypy, mako, psycopg2, pycurl). With them, the import smoke is
  fully green (0 unexpected failures) in both source and installed states.

## Running the gates (how to re-verify)

```sh
# 1. syntax
.venv/bin/python -m compileall -q ccpnmr2.5/python/          # 0 errors
# 2. import smoke (source tree)
MPLBACKEND=Agg .venv/bin/python import_smoke.py             # 1646 / 0 / 83
# 3. functional tests
.venv/bin/python -m pytest -q                                # 43 / 14 / 10 (10 = data-gated)
# 4. GUI launch (needs xvfb: apt install xvfb / conda-forge xvfb)
MPLBACKEND=Agg .venv/bin/python gui_boot_test.py            # 8/8 PASS
# 5. installed-state (fresh venv):
#    CCP_SMOKE_ROOT=<site-packages> import_smoke.py  -> 1637 / 0 / 83
#    pytest (copied test files) + gui_boot_test.py   -> 43/14/10, 8/8
```

`import_smoke.py` classifies the 83 modules that are non-importable *by design*
(commercial tools, macOS-only plugins, live-service scripts) — see its
`KNOWN_NON_IMPORTABLE` table.

## Repository layout

```
ccpnmr2.5/python/   the Python packages (memops, ccp, ccpnmr, cambridge, cing, ...)
ccpnmr2.5/model/    XML data model (loaded at import time)
ccpnmr2.5/data/     runtime data
ccpnmr2.5/doc/      in-program HTML documentation
ccpnmr2.5/c/        C sources for the 30 C extensions
import_smoke.py     whole-tree import smoke harness
gui_boot_test.py    Xvfb GUI launch harness
tests (ccpnmr2.5/python/tests)  functional pytest suite
```

## Scope notes

- Legacy single-purpose modules (clouds, haddock extender, some CASD/education
  scripts, ...) are **kept in the distribution but excluded from the functional
  test scope** — they import cleanly (smoke) but are not exercised end-to-end.
- `cing` remains part of the tree and its core modules import; the full cing
  stack (PyMOL, YASARA, Wattos, …) depends on external software that is not
  bundled (documented per-module in the smoke table).
- `survey.md` is the internal migration survey; the phase-by-phase audit trail
  lives in `_phase1a..4_checkpoints.md` / `_phase*_recipe.md`.

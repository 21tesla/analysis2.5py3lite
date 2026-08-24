# CCPNMR 2.5.2 (Python 3.13)

CCPNMR Analysis is a workbench for NMR (nuclear magnetic resonance) data analysis:
peak picking/fitting, assignment, restraints, structure calculation and validation —
built on the CCPN data model (MOPS) with a Tkinter GUI. This tree is the 2.5.2
release (circa 2020) modernized to run cleanly on **Python 3.13**, with the C
extension layers rebuilt and a pip-installable distribution.

## What works (verified gates)

All of the following are green on this tree (see `_phase4_checkpoints.md` for the
audit trail):

| Gate | Result |
|---|---|
| Whole-tree compile (Python 3.13) | 0 syntax errors |
| Import smoke — every module in the tree | 1277 attempted: 1265 OK / 2 failed (`cherrypy` web module, optional dep) / 10 documented-by-design |
| Functional test suite (pytest) | 45 passed / 4 skipped |
| C extensions | all setup.py C extensions import and pass functional checks |
| GUI launch under Xvfb | **4/4 console apps boot** (source and installed states) |
| Distribution | sdist + wheel build; clean-venv install passes all of the above |

## The four console commands

| Command | App |
|---|---|
| `ccpnmr` | CCPNMR Analysis (main workbench) |
| `ccpnmr-data-shifter` | Project data shifter |
| `ccpnmr-format-converter` | CCPN project format converter |
| `ccpnmr-update` | Automatic project update (non-GUI) |

Each app's `main()` entry was added/repaired in Phase 4 (P4-1) so the entry points
actually resolve.

## Installation

### Prebuilt distribution (recommended)

```sh
uv venv --python 3.13 .venv            # or: python3.13 -m venv .venv
uv pip install --python .venv/bin/python <ccpnmr-2.5.2-...whl>
# optional web/plotting stack (see "Optional features"):
uv pip install --python .venv/bin/python "matplotlib" "cherrypy" \
    "decorator" "mako"
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

`sdist install` gives the Python tree + data model. The C extensions are built
against the **running interpreter** via the root `setup.py` (MANIFEST.in ships
all C sources so any checkout can build):

```sh
python setup.py build_ext --inplace   # builds every ext
./scripts/copy_cext.sh               # places the .so files at their import sites
```

The Linux wheel ships the compiled extensions (cp313, linux_x86_64) and works
as-is; other platforms build from source as above.

### macOS

The Linux wheel does not install on macOS — the C extensions must be compiled
on the Mac (they are the only platform-compiled part). `setup.py`
auto-detects XQuartz, Homebrew (`tcl-tk`, `gsl`) and conda layouts, so a
naive user needs **no environment variables**.

**The easy way — Homebrew only, no conda/uv:**

```sh
xcode-select --install                       # one-time: cc
brew install --cask xquartz                  # one-time: X11 — log out/in after
brew install python-tk@3.13                  # Python 3.13 with Tk
python3.13 -m venv ~/ccpnmr && source ~/ccpnmr/bin/activate
pip install --upgrade pip setuptools
pip install numpy pandas PyOpenGL Pillow olefile requests python-dateutil pytz

git clone https://github.com/21tesla/analysis2.5py3lite.git && cd analysis2.5py3lite
python setup.py build_ext --inplace          # compiles every ext
./scripts/copy_cext.sh                       # places the .so files at import sites
pip install .                                # → the 4 console commands
ccpnmr                                       # main workbench (GUI)
```

**Conda variant** (e.g. if you are already a conda user):

```sh
conda create -n ccpnmr3 -c conda-forge python=3.13 tk tcl -y
conda activate ccpnmr3
git clone https://github.com/21tesla/analysis2.5py3lite.git && cd analysis2.5py3lite
python setup.py build_ext --inplace
./scripts/copy_cext.sh
pip install .            # or just run ./bin/analysis straight from the tree
```

macOS notes:

- macOS has no GLX: the GL window-handler extensions compile with `IGNORE_GL`
  (2D drawing, data layer and fitting keep full function; 3D GL structure
  rendering is off).
- Unusual layouts: override prefixes with `CCP_TK_PREFIX` (dir with
  `include/tk.h` + `lib/libtk8.6*`) and `CCP_X11_PREFIX` (default `/opt/X11`).

### Optional features

- **Web server + plotting extras** — `pip install ccpnmr[optional]`
  (matplotlib, cherrypy, decorator, mako). These cover the kept web-server
  module and the plotting code exercised by the by-design-kept packages.

## Running the gates (how to re-verify)

```sh
# 1. syntax
.venv/bin/python -m compileall -q ccpnmr2.5/python/          # 0 errors
# 2. import smoke (source tree)
MPLBACKEND=Agg .venv/bin/python import_smoke.py             # 1265 ok / 2 failed (cherrypy) / 10 by-design
# 3. functional tests
.venv/bin/python -m pytest -q                                # 45 passed / 4 skipped
# 4. GUI launch (needs xvfb: apt install xvfb / conda-forge xvfb)
MPLBACKEND=Agg .venv/bin/python gui_boot_test.py            # 4/4 PASS
# 5. installed-state (fresh venv + the optional extras): same gates, with
#    CCP_SMOKE_ROOT pointed at site-packages and FAILED required to be 0.
```

`import_smoke.py` classifies the 10 modules that are non-importable *by design*
(live/env-gated services, py2-only external deps) — see its
`KNOWN_NON_IMPORTABLE` table.

## Repository layout

```
ccpnmr2.5/python/   the Python packages (memops, ccp, ccpnmr, cambridge, ...)
ccpnmr2.5/model/    XML data model (loaded at import time)
ccpnmr2.5/data/     runtime data
ccpnmr2.5/doc/      in-program HTML documentation
ccpnmr2.5/c/        C sources for the C extensions
import_smoke.py     whole-tree import smoke harness
gui_boot_test.py    Xvfb GUI launch harness
tests (ccpnmr2.5/python/tests)  functional pytest suite
```

## Publishing (upstream)

- **Tag:** `v2.5.2-py3` (this release).
- **Install / publish docs:** [`docs/PUBLISHING.md`](./docs/PUBLISHING.md) — build env,
  gate set, `scripts/publish.sh` (build → verify → `twine check`/`upload`),
  plus a PyPI name caveat.
- **Conda-forge recipe (starting point):** [`recipe/meta.yaml`](./recipe/meta.yaml).

## Scope notes

- Legacy single-purpose modules (some CASD/education scripts, ...) are **kept in
  the distribution but excluded from the functional test scope** — they import
  cleanly (smoke) but are not exercised end-to-end.
- The 14 legacy third-party tool packages (PALES, MODULE, NOE/3J, DANGLE, ARIA,
  CYANA, HADDOCK, MECCANO, PyRPF, CING, ECI, Structure Viewer, H-bond
  restraints, ...) were removed in the 2026-08 simplification pass — see
  `SIMPLIFICATION_PLAN.md` and the git history for the stage-by-stage record.
- `survey.md` is the internal migration survey; the phase-by-phase audit trail
  lives in `_phase1a..4_checkpoints.md` / `_phase*_recipe.md`.

# CCPNMR 2.5.2 Codebase Survey & Modernization Plan

**Date:** 2026-08-20
**Decision:** Modernize to Python 3.13, keep Tkinter

---

## 1. What This Codebase Is

CCPNMR (CCP NMR) is an NMR (Nuclear Magnetic Resonance) data analysis application
for biochemists. It was started at the University of Cambridge in 2003 and now
has contributions from 10+ university groups.

The code is a layered framework:

```
memops/              — CINT core: MOPS (shared-memory object model),
                        XML metamodel, 74 Tkinter widget classes, editors
ccp/                 — NMR domain: NMR API (auto-generated, Nmr.py = 150K lines),
                        format readers/writers, math
ccpnmr/              — Application: main analysis GUI, ECI, NEXUS, CLOUDS, NEF
university plugins/  — 11 external tool packages (cing, cambridge, utrecht, ...)
```

The app is distributed as raw source files with a bundled conda Python 2.7
interpreter. There is no `setup.py`, no `pyproject.toml`, no wheels.

---

## 2. Current State — Critical Facts

### Environment
| Item | Value |
|---|---|
| Python | 2.7.17 (pinned in `ccpnmr2.5/c/environment_Linux.yml`) |
| Conda env | `miniconda/` (1.1 GB) checked into the repo |
| Distribution | `PYTHONPATH`-based via `bin/paths.sh`, shell-scripts in `bin/` |
| Linting | None |
| CI | None |
| Test framework | `nose 1.3.7` + stdlib `unittest` |
| Test coverage | ~88 files in `cing/` and `nef/`; **core `memops/` and `ccp/` have zero tests** |

### Dependency Pins (Linux, from `environment_Linux.yml`)
| Package | Version |
|---|---|
| python | 2.7.17 |
| numpy / numpy-base | 1.16.6 |
| pandas | 0.24.2 |
| pyopengl | 3.1.5 |
| pyopengl-accelerate | 3.1.5 |
| requests | 2.22.0 |
| six | 1.13.0 (installed but never imported) |
| nose | 1.3.7 |
| tk | 8.6.9 (from rare `smithsp` channel) |
| pillow | 6.2.1 |
| pyopenssl | 19.1.0 |
| openssl | 1.1.1g |
| cryptography | 2.8 |

### GUI Stack
- **Tkinter** — 241 Python files use `import Tkinter` (capital T, Python 2 only)
- **PyOpenGL 3.1.5** + `GlHandler.so` for 2D spectrum and 3D structure rendering
- **74 custom widget wrapper classes** in `memops/gui/` (Button, Frame, Canvas,
  TabbedFrame, PulldownMenu, ItemSelectPopup, ToolTip, DataEntry, ScrolledGraph, etc.)
- **No Qt** — zero PyQt/PySide/Qt imports anywhere in the codebase
- `TkHandler.so` is a native C extension that wraps Tkinter functionality

### Native C Extensions (10)
All compiled against the Python 2.7 C API. Sources in `ccpnmr2.5/c/`.

| .so | Purpose | Migration Note |
|---|---|---|
| `MemCache.so` | **MOPS** — shared-memory object model (data backbone) | Re-compile for py3.13 |
| `BlockFile.so` | Block file I/O | Re-compile |
| `StoreFile.so` | Store persistence | Re-compile |
| `StoreHandler.so` | Store interface | Re-compile |
| `FitMethod.so` | Fitting algorithms | Re-compile |
| `ShapeFile.so` | Shape file I/O (molsim) | Re-compile (good spike candidate) |
| `GlHandler.so` | OpenGL wrapper | Keep, re-compile (or replace with QGLWidget if ever Qt) |
| `PdfHandler.so` | PDF output | Re-compile or replace with PyPDF |
| `PsHandler.so` | PostScript output | Re-compile or replace |
| `TkHandler.so` | Tkinter wrapper | Keep, re-compile (or remove if GUI is refactored) |

### Code Scale
| Measure | Value |
|---|---|
| Python files | 1,720 |
| Total lines (~1.5M) | ~1.35M auto-generated API + ~736K hand-written |
| Largest file | `ccp/api/nmr/Nmr.py` — **150,534 lines** (auto-generated) |
| Second | `ccp/api/nmr/NmrConstraint.py` — 47,375 lines (auto-generated) |
| Third | `ccpnmr/api/Analysis.py` — 42,007 lines (auto-generated) |
| Wildcard imports (`from X import *`) | 262 files |
| Modules that are both importable and `__main__` scripts | 463 |
| Auto-generated `.pyo` files checked into repo | widespread |

---

## 3. Python 2.7 → 3.13 Migration: Full Inventory

### Syntactic (scriptable with `pyupgrade` / regex)
| Artifact | Files | Count | Target |
|---|---|---|---|
| `print 'x'` (statement, no parens) | widespread | **~4,165** | `print('x')` |
| `import Tkinter` (capital T) | **241** | 241 | `import tkinter` (lowercase) |
| `except E, e:` (removed in py3) | **96** | **163** | `except E as e:` |
| `apply(func, args)` (removed in py3) | **87** | **146** | `func(*args, **kwargs)` |
| `cStringIO` (py2 module) | **55** | 55 | `from io import StringIO` |
| `raw_input()` (removed in py3) | ~9 | 12 | `input()` |
| `cPickle` (py2 module) | **11** | 11 | `pickle` |
| `basestring` (removed in py3) | ~10 | ~10 | `str` |
| `unicode` (removed in py3) | ~15 | ~15 | `str` |
| `xrange` (removed in py3) | ~10 | ~10 | `range` |
| `<>` (removed in py3) | 7 | 12 | `!=` |
| `itertools.izip` (removed in py3) | 1 | 1 | `zip` |
| `ConfigParser` (case change in py3) | 2 | 2 | `configparser` |
| `httplib` (case change in py3) | 2 | 2 | `http.client` |
| `urllib2` (py2 module) | 4 | 4 | `urllib.request` |
| `cProfile` (case change) | 1 | 1 | `profile` |
| `execfile()` (removed in py3) | few | few | `exec(open(f).read())` |
| `sys.maxint` | few | few | `sys.maxsize` |
| `dict.has_key(k)` | few | few | `k in dict` |
| `library.zip` path hack | 1 | 1 | Remove (AnalysisGui.py) |

### Things That Are Already Partially OK
- `ccpnmr/Common.py`, `ccpnmr/v2io/`, `ccpnmr/nef/*.py` already have
  `from __future__ import print_function, unicode_literals, absolute_import`
- The codebase uses `threading` correctly (not the old `thread` module)
- `six` is installed in the env (1.13.0) but never actually imported — harmless

### C Extension Migration (for py3.13)
C API changes between Python 2.7 and 3.13:
- `Py_InitModule(name, methods)` → `PyModuleDef` + `PyModule_Create()`
- `PyInt_*` functions removed → `PyLong_*`
- `PyString_*` → `PyUnicode_*` or `PyBytes_*` (unicode vs bytes distinction)
- Implicit string encoding: `PyUnicode_AsUTF8()` for C-string extraction
- `#include <Python.h>` — verify path with py3.13 dev headers
- `PyDict_*` and `PyList_*` — mostly unchanged, verify
- Build system: `distutils` (removed in py3.12+) → `setuptools` (py≥3.12 compatible) or `meson`

### Dependency Updates (for py3.13)
| Current | Target | Notes |
|---|---|---|
| numpy 1.16.6 | ≥2.0 | `np.float_` → `np.float64`, `np.unicode_` → `np.str_` |
| pandas 0.24.2 | ≥2.0 | API changes in date handling, `.append()` removed |
| PyOpenGL 3.1.5 | ≥3.1.9 | Check GL context creation API |
| requests 2.22.0 | ≥2.31 | Mostly compatible |
| six 1.13.0 | Remove | Not used; can drop |
| nose 1.3.7 | pytest | `nose` is unmaintained; migrate 88 test files to pytest |
| smi
sp tk build | system tk | Python 3.13 ships its own tkinter bindings |
| enum34, ipaddress (py2 backports) | drop | Built into py3 stdlib |

---

## 4. Code Smells (Prioritized)

| Rank | Smell | Severity | Notes |
|---|---|---|---|
| 1 | **Python 2.7 syntax throughout** | Critical | 4165 print statements, 96 files `except E,`, `apply()` in 87 files |
| 2 | **Zero test coverage on core framework** | Critical | `memops/` and `ccp/` have no tests; 88 tests all in plugins |
| 3 | **10 C extensions on py2.7 C API** | Critical | Must be recompiled; `TkHandler.so` is Tk-specific |
| 4 | **`miniconda/` in git repo (1.1 GB)** | High | Not reproducible; bloats every clone |
| 5 | **No packaging infrastructure** | High | No `setup.py`, `pyproject.toml`, or wheels |
| 6 | **150K-line auto-generated `Nmr.py`** | High | Not testable, not type-checkable; generated from XML |
| 7 | **262 `from X import *`** | Medium | Breaks type checking; namespace pollution |
| 8 | **`AnalysisPopup.py` (~3000 lines, single file)** | Medium | Hard to maintain; no modular boundaries |
| 9 | **No CI, no linter** | Medium | No safety net for any change |
| 10 | **Auto-generated `.pyo` files in repo** | Low | Outdated py2 bytecode; should be gitignored |
| 11 | **463 dual-purpose files** (module + script) | Low | Works but unconventional; hard to test |
| 12 | **`miniconda/`, `miniconda/` pinned to exact 2019 builds** | Low | Fragile; `smithsp` channel may disappear |

---

## 5. Refactoring Plan

### Phase 0 — Foundation (Weeks 1–2)

**Goal: Set up the development environment so Phase 1 can be done safely.**

- [ ] `git rm -r miniconda/` + add `miniconda/` to `.gitignore`
- [ ] Remove all `.pyo` files: `find . -name "*.pyo" -delete`
- [ ] Create `pyproject.toml`:
  - `[build-system]` with `setuptools` backend
  - `[project]` name, version, python `>=3.13`
  - `[project.entry-points]` for the 37 bin scripts (replaces shell wrappers)
  - `meson` or CMake for C extensions (see Phase 1c)
- [ ] Create `uv` lockfile from `ccpnmr2.5/c/environment_Linux.yml` (mapped to py3.13)
- [ ] Set up minimal CI (GitHub Actions):
  - `uv sync` (create py3.13 env)
  - `ruff check ccpnmr2.5/python/` (lint check only)
  - `ruff format --check ccpnmr2.5/python/` (format check)
  - `pytest ccpnmr2.5/python/ccpnmr/nef/test/ -x` (run the best existing test suite)
- [ ] **C extension spike:** compile `shape_file.c` (smallest) against Python 3.13 C API,
  import it from py3.13, verify it works. *(~1–2 days)*

**Exit criteria:** Clean CI passes on Python 3.13 with `ruff` lint checks;
one C extension compiles and imports on py3.13.

---

### Phase 1 — Python 3.13 Port (Weeks 3–8)

**Goal: Application launches and core workflows work on Python 3.13 with Tkinter.**

#### 1a. Mechanical Syntax Fixes (Weeks 3–4)

Use `pyupgrade --py313-plus` + targeted sed/regex passes. Run on the full tree:

```bash
# Step 1: pyupgrade handles most of it
pyupgrade --py313-plus cc
pnmr2.5/python/**/*.py

# Step 2: Remaining manual fixes
# - except E, e: → except E as e:          (96 files / 163 hits)
# - apply(func, args) → func(*args)        (87 files / 146 hits)
# - from cStringIO import StringIO → from io import StringIO  (55 files)
# - import cPickle → import pickle         (11 files)
# - import Tkinter → import tkinter        (241 files)
# - raw_input → input                      (9 files)
# - xrange → range                         (10 files)
# - basestring → str                       (10 files)
# - unicode → str                          (15 files)
# - <> → !=                                (7 files)
# - from ConfigParser → from configparser  (2 files)
# - import httplib → import http.client    (2 files)
# - import urllib2 → import urllib.request (4 files)
# - execfile(f) → exec(open(f, encoding='utf-8').read())  (few)
```

**Verification:** `python -m compileall ccpnmr2.5/python/` must pass with zero errors.

#### 1b. C Extension Port (Weeks 4–7)

For each of the 10 C extensions:
1. Migrate to modern Python C API (`PyModuleDef`, `PyUnicode_*`, `PyLong_*`)
2. Build with `setuptools` (or `meson`)
3. Add to `pyproject.toml` build
4. Write a trivial import test

Priority order:
1. `ShapeFile.so` — smallest, good for validating approach *(sprint 1)*
2. `MemCache.so` — **data backbone**, must work before app can load projects
3. `BlockFile.so` — file I/O, needed for project open/save
4. `StoreFile.so` / `StoreHandler.so` — store persistence
5. `FitMethod.so` — spectrum fitting
6. `GlHandler.so` — OpenGL (verify PyOpenGL ≥3.1.9 works on py3.13; if so,
   this may be replaceable with pure-Python PyOpenGL calls)
7. `PdfHandler.so` / `PsHandler.so` — can be deferred (replace with PyPDF)
8. `TkHandler.so` — verify what it does; may be replaceable with `tkinter` stdlib
9. `copySharedObjs` — MOPS support, compile alongside `MemCache`

#### 1c. Dependency Updates (Weeks 5–8)
- [ ] `numpy ≥ 2.0` — test `np.float_` → `np.float64` migrations
- [ ] `pandas ≥ 2.0` — test `.append()` → `pd.concat()` in affected files
- [ ] `PyOpenGL ≥ 3.1.9` — verify `pyopengl-accelerate` is still available for py3.13 wheels
- [ ] `requests ≥ 2.31`
- [ ] Remove `six`, `enum34`, `ipaddress` (py2 backports, now stdlib)
- [ ] Replace `nose` test runner with `pytest`

#### 1d. Verification (Weeks 7–8)

Manual smoke tests:
- [ ] `analysis2.5` launches, main window appears, menus work
- [ ] Open a sample NMR project (data tree populates)
- [ ] Peak list displays, click a peak, spectrum canvas renders
- [ ] Run one relaxation analysis macro (uses `FitMethod.so`)
- [ ] ECI: `eci2.5` launches, opens a BMRB entry
- [ ] Run `nef/` test suite (88 tests)
- [ ] Run `cing/` test suite

**Exit criteria:** Full app works on Python 3.13 + Tkinter; test suites pass.

---

### Phase 2 — Quality Foundation (Weeks 9–16, overlaps Phase 3)

**Goal: Build enough safety net that the codebase is safe to refactor further.**

- [ ] Add `pytest` to CI (replaces `nose`)
- [ ] Write tests for `memops` core:
  - `memops/metamodel/MetaModel.py`
  - `memops/general/Io.py` (normalisePath, etc.)
  - `memops/general/Implementation.py` (ApiError, Application)
- [ ] Write tests for `ccp` core:
  - `ccp/api/general/Citation.py`
  - `ccp/api/molecule/ChemComp.py`
- [ ] Write tests for `ccpnmr` core:
  - `ccpnmr/analysis/core/` helpers
  - `ccpnmr/nef/NefIo.py` (already has tests, extend them)
- [ ] Target **>60% coverage** on `memops/` + `ccp/` non-auto-generated code
- [ ] Replace **262 `from X import *`** with explicit imports
- [ ] Remove all remaining `apply(func, args)` if `pyupgrade` missed any
- [ ] Run `ruff format` across the tree (standardises style)
- [ ] Add `py.typed` marker file
- [ ] Type annotations on `ccpnmr/api/Analysis.py` public API surface
- [ ] Run `mypy --check-untyped-defs memops/ ccp/` (expect errors; fix top-50)
- [ ] Decompose `AnalysisPopup.py` into sub-modules:
  - `analysis/menubar.py` — menu construction
  - `analysis/toolbar.py` — toolbar construction
  - `analysis/statusbar.py` — status bar
  - `analysis/callbacks.py` — business logic handlers

---

### Phase 3 — Packaging & Build (Weeks 12–16, parallel with Phase 2)

**Goal: Distribute as a proper Python package, not a source directory + shell scripts.**

- [ ] `pyproject.toml` fully populated:
  - `[project.scripts]` entry points (replaces `bin/*.5` shell scripts):
    ```toml
    [project.scripts]
    ccpnmr-analysis = "ccpnmr.analysis.AnalysisGui:main"
    ccpnmr-eci = "ccpnmr.eci.EntryCompletionGui:main"
    ccpnmr-format-converter = "ccpnmr.format.gui.FormatConverter:main"
    ```
  - `[build-system]` — `meson-python` (replaces `Makefile` for C ext)
- [ ] `CMakeLists.txt` / `meson.build` for all C extensions
- [ ] Build distributable wheels:
  ```bash
  uv build --wheel  # → dist/ccpnmr-2.5.3-py3-none-*.whl  (C ext per-platform)
  ```
- [ ] Dockerfile (multi-stage: build env → runtime):
  ```dockerfile
  FROM python:3.13-slim AS build
  ...
  FROM python:3.13-slim
  ...
  ```
- [ ] `pre-commit` config: `ruff`, `mypy` (core only)
- [ ] `README.md` rewrite: install instructions for py3.13 + `pip`

---

## 6. Effort Summary

| Phase | Duration | What You Get |
|---|---|---|
| **0 — Foundation** | 1–2 wks | Clean CI, py3.13 env, one C ext compiled |
| **1 — Python 3.13 port** | 3–6 wks | **App runs on py3.13 + Tkinter, all tests pass** |
| **2 — Quality foundation** | 4–8 wks | Test coverage, type hints, modularised GUI shell |
| **3 — Packaging/build** | 1–2 wks | `pip install`, wheels, Docker, no shell scripts |

**Total to a working Python 3.13 + Tkinter release: ~8–14 weeks**

(No Qt6 — Tkinter is available and working on Python 3.13; keeping it avoids a
3–6 month GUI rewrite for no functional gain in the short term.)

---

## 7. Key Files for Any Developer Starting Work

```
bin/paths.sh
  → How the environment is set up

ccpnmr2.5/c/environment_Linux.yml
  → Exact dependency pins (pre-migration)

ccpnmr2.5/python/ccpnmr/analysis/AnalysisGui.py
  → App bootstrap: creates Tk root, loads AnalysisPopup

ccpnmr2.5/python/ccpnmr/analysis/AnalysisPopup.py
  → Main GUI shell (~3,000 lines, all menus, toolbar, callbacks)

ccpnmr2.5/python/memops/general/Application.py
  → CINT framework app core

ccpnmr2.5/python/memops/gui/Button.py  (and 73 sibling files in memops/gui/)
  → Custom widget abstraction layer (Tkinter-based)

ccpnmr2.5/c/memops/global/shape_file.c
  → Simplest C extension, good first target for C API migration

ccpnmr2.5/python/ccp/api/nmr/Nmr.py
  → 150K-line auto-generated API class (generated from XML metamodel)
```

---

## 8. Risk Register

| Risk | Likelihood | Mitigation |
|---|---|---|
| C extension API incompatibility on py3.13 | Medium | Spike on `ShapeFile.c` in Phase 0, day 2–3 |
| `MemCache.so` (MOPS) breaks on py3.13 | Medium | Priority 1 in C ext migration; test project loading early |
| `GlHandler.so` + `pyopengl-accelerate` no py3.13 wheel | Low | Fall back to pure-Python PyOpenGL (slower but functional) |
| Wildcard imports mask circular dependencies at import time | Medium | Run `ruff` `F403` check early; fix one package at a time |
| Auto-generated API classes (150K lines) have py2 syntax missed by tools | Low | Regenerate from XML metamodel after port if tooling misses them |
| University plugin packages unmaintained — may not adapt | Medium | Scope Phase 1 to core `memops`/`ccp`/`ccpnmr`;
  plugins are best-effort after core works |
| `smithsp` channel disappears (provides the only working tk build) | Low | Phase 0: switch to system tkinter (py3.13 ships it) |

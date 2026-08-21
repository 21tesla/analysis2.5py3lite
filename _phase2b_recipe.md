# Phase 2b Recipe — fix remaining import-time failures (165 / 1726 failed at baseline)

Baseline (2026-08-21, after Phase 2a commit `e22a93f`):
- `MPLBACKEND=Agg .venv/bin/python import_smoke.py` → **1561 OK / 165 FAILED, 81 groups**
- Whole-tree compile: 0 syntax errors. Core tests 15 pass.

## Scope decisions (apply top-down)
1. **First, check whether the missing name is a FIRST-PARTY module that EXISTS in the tree**
   (py2 implicit-relative-import style: `import Foo` / `from Foo import X` where
   `Foo.py` sits next to the importing file). If it exists → make the import
   **relative** (`from .Foo import X` / `from ccpnmr.clouds.HydrogenDynamics import ...`
   matching the package root `ccpnmr2.5/python`). This was the Phase-1e pattern.
   Verify with: `find ccpnmr2.5/python -name "Foo.py"`.
2. **Pure py2→py3 code fixes**: removed stdlib mappings —
   - `cookielib` → `http.cookiejar`
   - `urlparse` → `urllib.parse`
   - `sgmllib` → `html.parser` (HTMLParser subclass) — rewrite, not a shim
   - `anydbm` → `dbm` (try dbm.gnu/dbm.ndbm fallback)
   - `module 'string' has no attribute 'replace'` → `str.translate` or plain replace
   - `'dict_items' object is not subscriptable` → index into `list(d.items())[i]` or use next(iter(...))
   - `dictionary changed size during iteration` → iterate `list(d.items())` snapshot
   - `'range' + 'list'` → `list(range(...)) + list`
   - `AwkLikeS.next`/`AwkLike.next` AttributeError → py2 `.next()` method removed; use `next(obj)` / `__next__`
   - `numpy.testing.utils` → `numpy.testing` (assert functions live on `np.testing` top level in numpy≥1.22)
   - `numpy.lib.twodim_base` → re-export shim removed; `np.triu`/`np.tril`/`np.diag` top-level
   - `numpy.lib.index_tricks` → `numpy.lib.index_tricks` may exist as `np.lib.index_tricks`; fallback `np.tri_indices` etc. Check runtime.
   - `NameError: cing` → add explicit `import cing` (Phase-1e pattern, `cing/core/parameters.py` precedent)
   - `NameError: Locator` (37) → investigate: likely `from SharedBeanService_services import Locator`
     exists (stub created Phase 1e) but import is implicit-relative or missing. Fix the import statement.
   - `ImportError: cannot_import X :: name` (3 groups) → module exists but `name` missing:
     either the name was renamed / removed in our port, or it lives in a different submodule —
     find the real location and fix the import (or add a re-export shim consistent with Phase-1e).
3. **Installable py3 external packages** (uv add, only if API is actually used compatibly):
   - `cherrypy` 18.x (webFc.py) — verify API use before installing (old cherrypy 3.x API differs)
   - `decorator` (pure python, stable API — TestNefIo.py)
   - `pycurl` (wheel available, API stable — toposcmd.py)
   - `pymol` (pypi `pymol` open-source wheel — 3 files)
   Skip install if the module uses a removed API that would break at runtime anyway;
   then classify as external.
4. **NOT installable on py3.13 (leave as-is, record in checkpoints)**:
   - `pymc` (PyMC2, py2-only), `sans` (py2 SOAP), `ccpncore` (internal), `yasara` (commercial, 7 files)
   - `CASD_HOME` env (4) — runtime config, not code (Phase-1e precedent)
   - `KeyError: ISD_ROOT / CYANA / ARG+`, `FileNotFoundError /Library/WebServer/...` (2),
     `SystemExit` (1), `SkipTest` (2) — likewise environment/data-gated.

## Rules (from feedback memory)
- Work file-by-file in the MAIN session; no subagents.
- After EVERY file edit: `PYTHONPATH=ccpnmr2.5/python .venv/bin/python -c "import <module>"` to confirm.
- At END of phase: full `import_smoke.py` run, record old→new counts in `_phase2b_checkpoints.md`.
- Stage ONLY files belonging to this phase. Never `git add .`. Leave `.qwen/settings.json`,
  `database.txt`, `dbTable-new`, `dbTable.new` untouched.
- Commit at phase boundary with message: `Fix remaining import-time failures under Python 3.13 (Phase 2b)`
  listing the biggest gains (smoke N→M).

## Commands
```bash
# Authoritative smoke (primary metric)
MPLBACKEND=Agg .venv/bin/python import_smoke.py
# Per-module probe
PYTHONPATH=ccpnmr2.5/python MPLBACKEND=Agg .venv/bin/python -c "import <top.package.name>"
# Find a first-party module
find ccpnmr2.5/python -name "ModuleName.py"
# Baseline snapshot (saved at phase start)
/tmp/smoke_2b_baseline.txt
```

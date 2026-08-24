#!/usr/bin/env python
"""Import-smoke test: import every module under ccpnmr2.5/python and classify
failures. Distinguishes real first-party migration bugs from missing optional
third-party deps (matplotlib/scipy). Run: .venv/bin/python import_smoke.py"""
import importlib
import os
import re
import sys
import traceback

# Default: the in-repo source tree next to this script. Override with
# CCP_SMOKE_ROOT to smoke-test an INSTALLED tree (site-packages) — used by the
# Phase-4 distribution gate (fresh venv -> pip install wheel -> smoke).
ROOT = os.environ.get(
    'CCP_SMOKE_ROOT', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ccpnmr2.5', 'python'))
sys.path.insert(0, ROOT)

# INSTALLED state: site-packages also holds third-party code (pip, numpy, ...).
# Restrict the walk to top-level entries of the installed ccpnmr distribution
# itself, read from its dist-info RECORD.
ALLOWED_TOP = None
if os.environ.get('CCP_SMOKE_ROOT'):
    record = os.path.join(ROOT, 'ccpnmr-2.5.2.dist-info', 'RECORD')
    ALLOWED_TOP = set()
    with open(record) as f:
        for line in f:
            path = line.split(',', 1)[0].strip()
            if path:
                top = path.split('/')[0]
                if not top.endswith('.dist-info'):
                    ALLOWED_TOP.add(top)

# Third-party modules that are NOT installed in this venv -> "missing dep", not a code bug.
MISSING_DEPS = {'scipy', 'tkinter', 'Tkinter'}  # Tkinter kept: we WANT to flag capital-Tkinter as a code smell? No -> see below.
# matplotlib is genuinely optional/heavy; missing = dep, not code bug.
OPT_MISSING = {'matplotlib', 'scipy', 'PIL', 'reportlab', 'pyproj', 'olefile', 'nose',
               # Optional third-party (not in the core distribution deps):
               'sqlalchemy', 'cherrypy', 'decorator', 'mako', 'psycopg2', 'pycurl'}

# ============================================================================
# Modules that are NON-IMPORTABLE BY DESIGN (NOT code regressions).
# Phase 2c, 2026-08-22: each of the remaining import failures was surveyed and
# given a concrete reason. Categories:
#   ENV    — needs an environment variable (CASD_HOME, ISD_ROOT) or live
#            network access at import time.
#   EXTERNAL — optional third-party software NOT bundled in the 2.5.2
#            distribution: py2-only packages (PyMC2, sans), commercial tools
#            (YASARA, PyMOL C++ engine), or internal sub-repos absent here
#            (Refine/protocol/UtilsAnalysis/pdbe-analysis/memops.scripts/...).
# If any of these ever starts importing cleanly, smoke prints a NOTE and the
# entry should be DELETED (a successful import then counts as OK).
KNOWN_NON_IMPORTABLE = {
    # --- EXTERNAL: optional software not in the original distribution ---
    "cambridge.bayes.PeakSeparatorPyMC": "EXTERNAL: PyMC2 — py2-only package, no py3 port",
    "ccp.lib.Bmrb.bmrb": "EXTERNAL: `sans` — py2-only SOAP stack",
    "ccp.util.V2Upgrade": "EXTERNAL: `ccpncore` — CCPN v2 internal core, not in this distribution",
    "pdbe.software.vascoReferenceCheck": "EXTERNAL: `pdbe.analysis` sub-repo not in this distribution",
    "pdbe.chemComp.export.setLicenses": "EXTERNAL: `memops.scripts` sub-package not in this distribution",
    # --- ENV: env vars / live network ---
    "cambridge.isd.isd_project_template": "ENV: ISD_ROOT environment variable",
    "nijmegen.CASD.Constants": "ENV: CASD_HOME environment variable",
    "nijmegen.CASD.Util": "ENV: CASD_HOME environment variable",
    "nijmegen.CASD.casdPipeLine": "ENV: CASD_HOME environment variable",
    "nijmegen.CASD.convertCasdNmrToCcpn": "ENV: CASD_HOME environment variable",
}

# Modules we can't expect to import cleanly in a headless/no-GUI env (GUI entry points, etc.)
# We still try them, but treat ImportError on tkinter as an environment note.

def classify(exc):
    """Return a (category, key) that groups failures by their concrete cause."""
    name = type(exc).__name__
    msg = str(exc)
    if name == 'ModuleNotFoundError':
        mod = getattr(exc, 'name', None) or 'unknown'
        top = mod.split('.')[0]
        if top in OPT_MISSING:
            return ('missing-dep', top)
        if top == 'tkinter':
            return ('env-gui', 'tkinter')
        return ('first-party-missing', top)
    if name == 'NameError':
        m = re.search(r"name '(\w+)' is not defined", msg)
        key = 'Name: ' + (m.group(1) if m else msg)
        return ('NameError', key)
    if name == 'ImportError' and 'cannot import name' in msg:
        m = re.search(r"cannot import name '(\w+)' from '([^']+)'", msg)
        if m:
            tgt = m.group(2).split('/')[-1] + ' :: ' + m.group(1)
            return ('ImportError', 'cannot_import ' + tgt)
        return ('ImportError', msg[:60])
    if name == 'ImportError' and 'undefined symbol' in msg:
        return ('stale-So', msg.split('undefined symbol:')[-1][:40])
    if name == 'AttributeError':
        m = re.search(r"'(\w+)' object has no attribute '(\w+)'", msg)
        if m:
            return ('AttributeError', f'{m.group(1)}.{m.group(2)}')
        return ('AttributeError', msg[:60])
    return (name, msg[:60])


files = []
for dirpath, dirnames, filenames in os.walk(ROOT):
    dirnames[:] = [d for d in dirnames if d not in ('__pycache__', 'build', 'node_modules')]
    for fn in filenames:
        if fn.endswith('.py'):
            full = os.path.join(dirpath, fn)
            rel = os.path.relpath(full, ROOT)
            modname = rel[:-3].replace(os.sep, '.')
            if modname.endswith('.__init__'):
                modname = modname[:-len('.__init__')]
            files.append((rel, modname))

files.sort(key=lambda x: x[1])

ok = 0
failures = []  # (rel, modname, category, top, exc_str)
allowed = []  # (rel, reason) — failed, but documented-by-design (see KNOWN_NON_IMPORTABLE)
allowed_now_ok = []  # allowlisted modules that now import cleanly (remove their entry!)
seen_mod = set()
for rel, modname in files:
    if modname in seen_mod:
        continue
    seen_mod.add(modname)
    if ALLOWED_TOP is not None:
        top = modname.split('.')[0]
        if top in ('model', 'data', 'doc', 'license'):
            continue  # pure runtime data dirs shipped at the wheel root — not code
        if top not in ALLOWED_TOP:
            continue
    try:
        importlib.import_module(modname)
        if modname in KNOWN_NON_IMPORTABLE:
            allowed_now_ok.append((rel, modname))
        ok += 1
    except BaseException as e:  # noqa: BLE001
        if modname in KNOWN_NON_IMPORTABLE:
            allowed.append((rel, KNOWN_NON_IMPORTABLE[modname]))
            continue
        cat, top = classify(e)
        exc_str = f'{type(e).__name__}: {e}'
        exc_str = exc_str.replace(os.path.abspath(ROOT), ROOT)
        exc_str = exc_str.replace(ROOT + '/', '')
        failures.append((rel, cat, top, exc_str))

# --- Report ---
print(f'TOTAL modules attempted: {len(seen_mod)}')
print(f'  OK:            {ok}')
print(f'  FAILED:        {len(failures)}   <- unexpected, investigate these')
print(f'  BY-DESIGN:     {len(allowed)}   (documented in KNOWN_NON_IMPORTABLE: not regressions)')
if allowed_now_ok:
    print()
    for rel, modname in allowed_now_ok:
        print(f'  NOTE: allowlisted module now imports cleanly — DELETE its entry: {rel} ({modname})')
print()

from collections import Counter, defaultdict
by_top = defaultdict(list)
for rel, cat, top, exc_str in failures:
    by_top[(cat, top)].append((rel, exc_str))

print('=== Failures grouped by (category, missing_module) ===')
for (cat, top), items in sorted(by_top.items(), key=lambda kv: (-len(kv[1]), kv[0])):
    print(f'\n### cat={cat}  missing={top!r}   ({len(items)} files)')
    for rel, exc_str in items[:8]:
        print(f'    {rel}  ::  {exc_str[:120]}')
    if len(items) > 8:
        print(f'    ... and {len(items)-8} more')

# --- By-design (allowlisted) modules, grouped by category (audit trail) ---
allowed_by_cat = defaultdict(list)
for rel, reason in allowed:
    allowed_by_cat[reason.split(':', 1)[0].strip()].append(rel)
if allowed_by_cat:
    print('\n=== BY-DESIGN non-importable modules (documented, NOT regressions) ===')
    for cat in sorted(allowed_by_cat):
        print(f'\n### {cat}   ({len(allowed_by_cat[cat])})')
        for rel in allowed_by_cat[cat]:
            print(f'    {rel}')

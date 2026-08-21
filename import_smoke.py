#!/usr/bin/env python
"""Import-smoke test: import every module under ccpnmr2.5/python and classify
failures. Distinguishes real first-party migration bugs from missing optional
third-party deps (matplotlib/scipy). Run: .venv/bin/python import_smoke.py"""
import importlib
import os
import re
import sys
import traceback

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ccpnmr2.5', 'python')
sys.path.insert(0, ROOT)

# Third-party modules that are NOT installed in this venv -> "missing dep", not a code bug.
MISSING_DEPS = {'scipy', 'tkinter', 'Tkinter'}  # Tkinter kept: we WANT to flag capital-Tkinter as a code smell? No -> see below.
# matplotlib is genuinely optional/heavy; missing = dep, not code bug.
OPT_MISSING = {'matplotlib', 'scipy', 'PIL', 'reportlab', 'pyproj', 'olefile'}

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
seen_mod = set()
for rel, modname in files:
    if modname in seen_mod:
        continue
    seen_mod.add(modname)
    try:
        importlib.import_module(modname)
        ok += 1
    except BaseException as e:  # noqa: BLE001
        cat, top = classify(e)
        exc_str = f'{type(e).__name__}: {e}'
        exc_str = exc_str.replace(os.path.abspath(ROOT), ROOT)
        exc_str = exc_str.replace(ROOT + '/', '')
        failures.append((rel, cat, top, exc_str))

# --- Report ---
print(f'TOTAL modules attempted: {len(seen_mod)}')
print(f'  OK:            {ok}')
print(f'  FAILED:        {len(failures)}')
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

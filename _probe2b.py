"""Phase 2b probe: single pass over the tree (same methodology as import_smoke),
list every module whose import fails with an error matching <needle>.
Usage: _probe2b.py <needle>
"""
import importlib
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ccpnmr2.5", "python")
sys.path.insert(0, ROOT)
os.environ.setdefault("MPLBACKEND", "Agg")
needle = sys.argv[1]

files = []
for dirpath, dirnames, filenames in os.walk(ROOT):
    dirnames[:] = [d for d in dirnames if d not in ("__pycache__", "build", "node_modules")]
    for fn in filenames:
        if fn.endswith(".py"):
            full = os.path.join(dirpath, fn)
            rel = os.path.relpath(full, ROOT)
            modname = rel[:-3].replace(os.sep, ".")
            if modname.endswith(".__init__"):
                modname = modname[: -len(".__init__")]
            files.append((rel, modname))
files.sort(key=lambda x: x[1])

seen = set()
hits = []
for rel, modname in files:
    if modname in seen:
        continue
    seen.add(modname)
    try:
        importlib.import_module(modname)
    except BaseException as e:  # noqa: BLE001
        s = f"{type(e).__name__}: {e}"
        if needle in s:
            print(f"{rel}  ::  {s[:130]}")
            hits.append(modname)
print(f"TOTAL {len(hits)}")

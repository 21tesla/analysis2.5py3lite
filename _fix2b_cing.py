"""Phase 2b: add missing `import cing` to files that use `cing.` at module level
but only do `from cing.X import Y` (py2 star-import leak no longer reliable).
Deterministic + idempotent. Verifies each file compiles afterwards.
"""
import ast
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ccpnmr2.5", "python")
sys.path.insert(0, ROOT)
os.environ.setdefault("MPLBACKEND", "Agg")

files = [
    "cing/Database/Scripts/createStarUserLib.py",
    "cing/Database/Scripts/mkINTERNAL_1.py",
    "cing/Database/Scripts/test/test_createStarUserLib.py",
    "cing/NRG/doAnnotateCaspNmrLoop.py",
    "cing/NRG/doAnnotateNrgCing.py",
    "cing/NRG/getRCSB_PDB.py",
    "cing/NRG/runQueenyAll.py",
    "cing/Scripts/CING_paper_queries.py",
    "cing/Scripts/FC/convertCyana2Ccpn.py",
    "cing/Scripts/FC/convertPdb2Ccpn.py",
    "cing/Scripts/FC/convertStar2Ccpn.py",
    "cing/Scripts/FC/convertXplor2Ccpn.py",
    "cing/Scripts/FC/mergeNrgBmrbShifts.py",
    "cing/Scripts/interactive/mouseBuffer.py",
    "cing/Scripts/interactive/mouseBuffer2.py",
    "cing/Scripts/interactive/mouseBuffer3.py",
    "cing/Scripts/publishVC.py",
    "cing/Scripts/updateXplorConv.py",
    "cing/Scripts/validateForNmrCmbi.py",
    "cing/Scripts/validateForProteinsDotDynDnsDotOrg.py",
]

import re
from cing import (  # noqa: F401  (sanity: package exposes the attrs used)
    verbosityDebug,
    verbosity,
)

def find_insert_index(lines):
    """Index of line before which `import cing` should be inserted.
    Prefer: immediately before first col-0 `from cing`/`import cing.` line.
    Fallback: after last col-0 import statement."""
    first_cing_imp = None
    last_imp = None
    for i, ln in enumerate(lines):
        s = ln.rstrip("\n")
        if first_cing_imp is None and re.match(r"^(from cing[. ]|import cing[. ])", s):
            first_cing_imp = i
        if re.match(r"^(import \w|from \w)", s):
            last_imp = i
    if first_cing_imp is not None:
        return first_cing_imp
    if last_imp is not None:
        return last_imp + 1
    return None

fixed, skipped = [], []
for rel in files:
    path = os.path.join(ROOT, rel)
    with open(path, encoding="utf-8") as fh:
        src = fh.read()
    if re.search(r"^import cing\s*$", src, flags=re.M):
        skipped.append(rel)
        continue
    lines = src.splitlines(keepends=True)
    idx = find_insert_index(lines)
    if idx is None:
        print("NO-INSERT-POINT", rel)
        continue
    # ensure a blank line above if previous line is non-blank
    if idx > 0 and lines[idx - 1].strip() != "":
        lines.insert(idx, "\n")
        idx += 1
    lines.insert(idx, "import cing\n")
    new = "".join(lines)
    ast.parse(new)  # must compile
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(new)
    fixed.append(rel)

print(f"fixed {len(fixed)}, already-had-import {len(skipped)}")

#!/usr/bin/env bash
# Place the freshly built C extensions where the importers find them.
#
# `python setup.py build_ext --inplace` (root setup.py) emits flat,
# top-level-named shared objects (e.g. ShapeFile, PeakList) in the project
# root, because the extension names are top-level.  Run THIS script after
# the build to install them at their two import sites:
#
#   1. onto the per-package symlink targets
#      ccpnmr2.5/python/<pkg>/c/<Name>.so  ->  c/.../<Name>.so
#      (package imports: memops.c.StoreFile, ccp.c.StructAtom, ccpnmr.c.PeakList, ...)
#   2. into ccpnmr2.5/python/               (flat imports: import ShapeFile)
#
set -euo pipefail
cd "$(dirname "$0")/.."

n=0
for so in *.so; do
  [ -f "$so" ] || continue        # no match: glob stayed literal
  name="${so%%.*}"                # module name = first dot-separated component
  for link in $(find ccpnmr2.5/python -path "*/c/${name}.so" -type l); do
    cp -f "$so" "$link"           # writes through the symlink to its .so target
    n=$((n+1))
  done
  cp -f "$so" "ccpnmr2.5/python/$so"
  n=$((n+1))
done
echo "placed $n built extension(s) at their import sites"

#!/usr/bin/env bash
# Place the freshly built C extensions where the importers find them.
#
# `python setup.py build_ext --inplace` (root setup.py) emits flat,
# top-level-named shared objects (e.g. ShapeFile, PeakList) in the project
# root, because the extension names are top-level.  Run THIS script after
# the build to install them at their import sites:
#
#   1. into their per-package c/ directory:
#      ccpnmr2.5/python/<pkg>/c/<Name>.cpython-313-darwin.so
#      (package imports: memops.c.StoreFile, ccp.c.StructAtom, ccpnmr.c.PeakList, ...)
#   2. into ccpnmr2.5/python/               (flat imports: import ShapeFile)
#
set -euo pipefail
cd "$(dirname "$0")/.."

get_target_dir() {
  case "$1" in
    StructAtom|StructBond|StructUtil|StructStructure)
      echo "ccpnmr2.5/python/ccp/c" ;;
    BayesPeakSeparator)
      echo "ccpnmr2.5/python/cambridge/c" ;;
    ShapeFile|MemCache|BlockFile|FitMethod|StoreFile|StoreHandler|PdfHandler|PsHandler|GlHandler|TkHandler)
      echo "ccpnmr2.5/python/memops/c" ;;
    *)
      echo "ccpnmr2.5/python/ccpnmr/c" ;;
  esac
}

n=0
for so in *.so; do
  [ -f "$so" ] || continue        # no match: glob stayed literal
  name="${so%%.*}"                # module name = first dot-separated component
  
  # Copy to the flat python dir
  cp -f "$so" "ccpnmr2.5/python/$so"
  n=$((n+1))
  
  # Copy to the package c/ dir
  target_dir="$(get_target_dir "$name")"
  mkdir -p "$target_dir"
  cp -f "$so" "$target_dir/$so"
  # Also copy without the ABI tag (just as .so) to support any legacy paths looking for it
  cp -f "$so" "$target_dir/${name}.so"
  n=$((n+2))
done
echo "placed $n built extension(s) at their import sites"

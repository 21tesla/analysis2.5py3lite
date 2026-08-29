#!/usr/bin/env bash
# Build the Linux x86_64 "standalone" (extract-and-run) distribution:
#
#   dist/ccpnmr-<ver>-linux-x86_64-standalone.tar.gz
#
# A fresh user unpacks the tarball and runs ./bin/analysis - no system
# python, no pip, no venv: the tree embeds a private CPython 3.14 runtime
# (from-source, Tcl/Tk 9-linked) with the ccpnmr wheel (and all deps)
# installed into its own site-packages, plus the Tcl/Tk 9 shared libraries
# bundled under lib/tcl9.  Same model as the macOS portable distribution
# described in INSTALL.md.
#
# Inputs (overridable):
#   CCP_PYTHON_PREFIX  CPython 3.14 prefix built with Tcl/Tk 9
#                      (default /usr/local; must contain bin/python3.14
#                      whose _tkinter links the tcl9 libraries)
#   CCP_TCL9_LIBS_DIR  dir holding libtcl9.0.so + libtcl9tk9.0.so
#                      (default /usr/local/lib)
#
# Usage:  ./make-standalone-linux.sh   (needs: uv, gcc not required - the
#         C extensions are pre-built and ride in the wheel)
set -euo pipefail

# the script lives in the repo root; dirname handles both `./script` and
# absolute-path invocations
TOP="$(cd "$(dirname "$0")" && pwd)"
cd "$TOP"

TCL9_PREFIX="${CCP_PYTHON_PREFIX:-/usr/local}"
TCL9_LIBS_DIR="${CCP_TCL9_LIBS_DIR:-/usr/local/lib}"

# --- version + architecture ---------------------------------------------------
VER="$(grep -m1 '^version' pyproject.toml | cut -d'"' -f2)"
ARCH="x86_64"
NAME="ccpnmr-${VER}-linux-${ARCH}-standalone"
STAGE="dist/${NAME}"
PYTAG="3.14"                       # must match the wheel's cp314 tags
CPYPYBIN="${TCL9_PREFIX}/bin/python${PYTAG}"

echo "==> standalone dist: ${NAME}"
[ -x "$CPYPYBIN" ] || { echo "missing interpreter: $CPYPYBIN (set CCP_PYTHON_PREFIX)"; exit 1; }
[ -f "$TCL9_LIBS_DIR/libtcl9.0.so" ] && [ -f "$TCL9_LIBS_DIR/libtcl9tk9.0.so" ] || {
  echo "missing Tcl/Tk 9 libs in $TCL9_LIBS_DIR (set CCP_TCL9_LIBS_DIR)"; exit 1; }

# --- 1. fresh wheel from the working tree -------------------------------------
echo "==> building wheel (uv build --wheel --python $CPYPYBIN)"
uv build --wheel --python "$CPYPYBIN"
WHEEL="dist/ccpnmr-${VER}-cp314-cp314-linux_${ARCH}.whl"
[ -f "$WHEEL" ] || { echo "wheel not found: $WHEEL"; exit 1; }

# --- 2. private relocatable CPython 3.14 runtime ------------------------------
# A from-source prefix install is relocatable: bin/ + lib/ keep their relative
# layout, so the staged tree stands up anywhere.  Stage only what the runtime
# needs (static libpython is baked into the binary; include/ is build-time).
echo "==> runtime: $TCL9_PREFIX (Python 3.14, Tcl/Tk 9)"
"$CPYPYBIN" -c "import _tkinter" || { echo "interpreter lacks _tkinter - abort"; exit 1; }

rm -rf "$STAGE"
mkdir -p "$STAGE/bin" "$STAGE/runtime/bin" "$STAGE/lib/tcl9"
cp -a "$TCL9_PREFIX/bin/python${PYTAG}" "$STAGE/runtime/bin/python${PYTAG}"
ln -s "python${PYTAG}" "$STAGE/runtime/bin/python3"
ln -s "python3"        "$STAGE/runtime/bin/python"
mkdir -p "$STAGE/runtime/lib"
cp -a "$TCL9_PREFIX/lib/python${PYTAG}" "$STAGE/runtime/lib/python${PYTAG}"
chmod -R u+w "$STAGE/runtime"            # the source prefix may be root-owned
# Tcl/Tk 9 shared libraries (the _tkinter + C-extension NEEDED names):
cp -a "$TCL9_LIBS_DIR/libtcl9.0.so" "$TCL9_LIBS_DIR/libtcl9tk9.0.so" "$STAGE/lib/tcl9/"

# --- 3. install the wheel (+ deps) into the private runtime --------------------
echo "==> installing wheel into private runtime"
uv pip install --python "$STAGE/runtime/bin/python" "$WHEEL" 2>/dev/null \
  || uv pip install --python "$STAGE/runtime/bin/python" --break-system-packages "$WHEEL"
# sanity: the package must import from INSIDE the staged runtime, full stop
# (cd / proves independence from the source tree - hence absolute paths);
# the bundled tcl9 libs are on LD_LIBRARY_PATH exactly as the launchers set it
(cd / && LD_LIBRARY_PATH="$TOP/$STAGE/lib/tcl9" "$TOP/$STAGE/runtime/bin/python" -W ignore -c "
import tkinter, ccpnmr.analysis.AnalysisGui as G, ccpnmr.nefCli as N
assert G.__file__.startswith('$TOP/$STAGE/runtime')
print('runtime import OK (tkinter tcl %s):' % (tkinter.TclVersion,), G.__file__)")

# --- 4. portable launcher -------------------------------------------------------
# The wheel ships the packages AND the model/doc data at the package root
# (site-packages), so no source tree / PYTHONPATH is needed - getTopDirectory()
# understands the installed layout.  -m keeps the entry layout-independent.
# lib/tcl9 goes first on LD_LIBRARY_PATH so the bundled Tcl/Tk 9 wins even on
# hosts that also carry it in /usr/local (and satisfies it on hosts that do not).
cat > "$STAGE/bin/paths.sh" <<'EOF'
#!/usr/bin/env bash
# Standalone (portable) layout: everything lives under the unpacked tree.
CCPNMR_TOP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit; pwd)"
export CCPNMR_TOP_DIR
export CONDA="${CCPNMR_TOP_DIR}"/runtime
export LD_LIBRARY_PATH="${CCPNMR_TOP_DIR}/lib/tcl9${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
EOF
cat > "$STAGE/bin/analysis" <<'EOF'
#!/usr/bin/env bash
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" || exit && pwd)"
source "${SCRIPTDIR}/paths.sh"
exec "${CONDA}"/bin/python -O -W ignore -m ccpnmr.analysis.AnalysisGui "$@"
EOF
chmod +x "$STAGE/bin/analysis" "$STAGE/bin/paths.sh"

# pip bakes the STAGING path into console-script shebangs, which dangle once
# the user unpacks somewhere else.  Replace the two ccpnmr-generated scripts
# with self-locating wrappers (both modules have __main__ guards).
for spec in ccpnmr:ccpnmr.analysis.AnalysisGui ccpnmr-nef:ccpnmr.nefCli; do
  cmd="${spec%%:*}"; mod="${spec#*:}"
  cat > "$STAGE/runtime/bin/$cmd" <<EOF
#!/bin/sh
DIR=\$(CDPATH= cd -- "\$(dirname -- "\$0")" && pwd)
TOP=\$(CDPATH= cd -- "\$DIR/../.." && pwd)
export LD_LIBRARY_PATH="\$TOP/lib/tcl9\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
exec "\$DIR/python" -m ${mod} "\$@"
EOF
  chmod +x "$STAGE/runtime/bin/$cmd"
done

cat > "$STAGE/README-PORTABLE.txt" <<EOF
CCPNMR Analysis ${VER} - Linux x86_64 standalone
================================================

Self-contained: no system Python, no pip, no installation step.
The private CPython ${PYTAG} runtime (built against Tcl/Tk 9), all
dependencies, and the Tcl/Tk 9 shared libraries (lib/tcl9) are embedded
in this tree.

Run:
  cd $(basename "$STAGE")
  ./bin/analysis                 (optionally: ./bin/analysis /path/to/project)

Requirements on the host (a normal Linux desktop):
  - graphics: X11 libraries (libX11/libGL/libglut/libfontconfig/libXft -
    present on any desktop install).  Nothing else is read from or
    written to the system (except the project you open/save and your
    default browser for Project > Summary)

CLI utilities (non-GUI):
  ./runtime/bin/ccpnmr-nef import file.nef [--project-name NAME] [--force]
  ./runtime/bin/ccpnmr-nef export <project-directory> <output.nef>
EOF

# --- 5. pack --------------------------------------------------------------------
echo "==> packing"
tar czf "dist/${NAME}.tar.gz" -C dist "$NAME"
rm -rf "$STAGE"
ls -lh "dist/${NAME}.tar.gz"
echo "==> done: dist/${NAME}.tar.gz"

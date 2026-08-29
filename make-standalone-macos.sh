#!/usr/bin/env bash
# Build the macOS (arm64 / x86_64) "standalone" (extract-and-run) distribution:
#
#   dist/ccpnmr-<ver>-macos-<arch>-standalone.tar.gz
#
# A fresh user unpacks the tarball and runs ./bin/analysis - no system python,
# no pip, no venv: the tree embeds a private CPython 3.13 runtime (uv-managed
# python-build-standalone, relocatable) with the ccpnmr wheel (dependencies
# included) installed into its own site-packages.  Same model as the Linux
# standalone (make-standalone-linux.sh) and the portable macOS distribution
# described in INSTALL.md.
#
# RUNS ON THE TARGET MAC.  Unlike the Linux build, the 30 C extensions CANNOT
# ride prebuilt from the tree (the in-tree .so files are Linux binaries) -
# `uv build` compiles them here with the host compiler.  The archive therefore
# only runs on the architecture it was built on.
#
# Prerequisites on the build Mac:
#   - uv            (curl -LsSf https://astral.sh/uv/install.sh | sh)
#   - C compiler    (Xcode Command Line Tools: xcode-select --install)
#   - Tcl/Tk 8 headers for the C extensions (section 1b): a conda env at
#     <conda>/envs/analysis (`conda create -n analysis python=3.13 tcl-tk=8`),
#     auto-detected - or set CCP_TK_PREFIX to any prefix with
#     include/tcl-tk/tk.h (Tcl/Tk 8).  The C extensions LINK the private
#     runtime's OWN bundled Tcl/Tk 8.6 dylibs (not the conda env's), so a
#     fresh user's mac needs neither conda nor tcl-tk - XQuartz only.
#
# Usage:  ./make-standalone-macos.sh
set -euo pipefail

case "$(uname -s)" in
  Darwin*) ;;
  *) echo "error: this script must be run ON macOS (it compiles the C extensions for the host)"; exit 1 ;;
esac

TOP="$(cd "$(dirname "$0")" && pwd)"
cd "$TOP"

command -v uv >/dev/null 2>&1 || { echo "error: uv is required (curl -LsSf https://astral.sh/uv/install.sh | sh)"; exit 1; }
command -v cc >/dev/null 2>&1 || { echo "error: no C compiler found - install the Xcode Command Line Tools: xcode-select --install"; exit 1; }

# --- version + architecture ---------------------------------------------------
# awk instead of grep|cut: BSD-grep-safe and pipefail-clean
VER="$(awk -F'"' '/^version/{print $2; exit}' pyproject.toml)"
ARCH="$(uname -m)"                    # arm64 (Apple Silicon) or x86_64 (Intel)
NAME="ccpnmr-${VER}-macos-${ARCH}-standalone"
STAGE="dist/${NAME}"
PYTAG="3.13"                          # must match the wheel's cp313 tags

echo "==> standalone dist: ${NAME} (built on $(sw_vers -productVersion 2>/dev/null || echo 'macOS'))"

# --- 1. private relocatable interpreter (also pins the BUILD interpreter) ------
# The venv/homebrew/pythons here are NOT relocatable; the runtime must be a
# uv-managed python-build-standalone (self-contained, moves anywhere).
UVPY="$(uv python find "$PYTAG" --system --python-preference only-managed 2>/dev/null || true)"
if [ -z "$UVPY" ]; then
  echo "==> installing managed CPython ${PYTAG} (one-time network download)"
  uv python install "$PYTAG"
  UVPY="$(uv python find "$PYTAG" --system --python-preference only-managed)"
fi
UVHOME="$(cd "$(dirname "$UVPY")/.." && pwd)"
echo "==> runtime: $UVHOME"
"$UVPY" -c "import _tkinter" || { echo "error: managed python lacks _tkinter - abort"; exit 1; }

# --- 1b. Tcl/Tk 8 headers for the C extensions --------------------------------
# The C draw extensions must link the SAME Tcl/Tk the runtime _tkinter loads -
# that is this uv-managed CPython's bundled Tcl/Tk 8.6 dylibs ($UVHOME/lib,
# which becomes the staged runtime's lib; setup.py puts sys.base_prefix/lib
# first in the linker search, so -ltk8.6/-ltcl8.6 pick those exact dylibs and
# the install names recorded in the .so match _tkinter's: ONE Tk copy per
# process).  A second Tk (Homebrew's tcl9, or the conda env's own dylibs)
# makes the canvas-XOR crosshair ghost over the spectrum.  python-build-
# standalone ships NO tcl/tk headers, so take them from a Tcl/Tk 8 prefix:
# the conda 'analysis' env (auto-detect) or $CCP_TK_PREFIX.  The conda python
# itself is not used for the build.
if [ -z "${CCP_TK_PREFIX:-}" ]; then
  for c in "$HOME/miniconda3/envs/analysis" "$HOME/anaconda3/envs/analysis" \
           "$HOME/miniforge3/envs/analysis" /opt/miniconda3/envs/analysis \
           /opt/anaconda3/envs/analysis; do
    if [ -f "$c/include/tcl-tk/tk.h" ]; then
      CCP_TK_PREFIX="$c"
      break
    fi
  done
fi
if [ -z "${CCP_TK_PREFIX:-}" ]; then
  echo "error: no Tcl/Tk 8 header prefix found (needed for the C extensions)."
  echo "  create the conda env (used for its headers only):"
  echo "    conda create -n analysis python=3.13 tcl-tk=8"
  echo "  or point CCP_TK_PREFIX at any prefix with include/tcl-tk/tk.h (Tcl/Tk 8)."
  exit 1
fi
[ -f "$CCP_TK_PREFIX/include/tcl-tk/tk.h" ] || { echo "error: $CCP_TK_PREFIX has no include/tcl-tk/tk.h"; exit 1; }
echo "==> Tcl/Tk headers: $CCP_TK_PREFIX (dylibs: the runtime python's own)"
export CCP_TK_PREFIX

# --- 2. fresh wheel, compiled for THIS mac --------------------------------------
# --python pins the build interpreter to the managed 3.13: a newer Homebrew
# python on the mac would otherwise yield cp3xx extensions the private runtime
# (3.13) can't load.  This is the slow step (30 C extensions).
echo "==> building wheel (compiles the C extensions for macos/${ARCH})"
uv build --wheel --python "$UVPY"
# setuptools platform tags vary (macosx_11_0_arm64 vs macosx_10_9_x86_64)
WHEEL="$(ls dist/ccpnmr-${VER}-cp313-*.whl 2>/dev/null | grep 'macosx' | head -n1 || true)"
[ -n "$WHEEL" ] && [ -f "$WHEEL" ] || { echo "error: no macos wheel found in dist/ (expected ccpnmr-${VER}-cp313-*-macosx*.whl)"; exit 1; }
echo "==> wheel: $WHEEL"

# --- 3. stage the tree -----------------------------------------------------------
rm -rf "$STAGE"
mkdir -p "$STAGE/bin"
cp -a "$UVHOME" "$STAGE/runtime"
chmod -R u+w "$STAGE/runtime"      # uv's pythons are r-x; pip needs to write

# --- 4. install the wheel (+ deps) into the private runtime ----------------------
echo "==> installing wheel into private runtime"
PYTHONHOME="$TOP/$STAGE/runtime" "$TOP/$STAGE/runtime/bin/python" -m pip install --force-reinstall "$WHEEL" --break-system-packages
# sanity: the package must import from INSIDE the staged runtime, full stop
# (cd / proves independence from the source tree - hence absolute paths)
(cd / && PYTHONHOME="$TOP/$STAGE/runtime" "$TOP/$STAGE/runtime/bin/python" -c "
import ccpnmr.analysis.AnalysisGui as G, ccpnmr.nefCli as N, ccp.gui.Io as Io
assert G.__file__.startswith('$TOP/$STAGE/runtime')
assert hasattr(Io, 'loadNefProject')
print('runtime import OK:', G.__file__)")

# --- 4b. one-Tk-per-process audit ----------------------------------------------
# Proves 1b held: every tcl/tk reference in the installed C extensions must
# resolve to the runtime's own Tcl/Tk 8.6 copy.  @rpath-family names resolve
# through the python binary's rpath to runtime/lib.  Any OTHER tcl/tk
# (Homebrew tcl9, the conda env's dylibs, ...) = a second Tk in the process =
# the crosshair ghosting.
echo "==> auditing tcl/tk linkage of the installed C extensions"
AUDIT_TMP="$(mktemp)"
find "$STAGE/runtime/lib/python3.13/site-packages" \( -name '*.so' -o -name '*.dylib' \) -type f > "$AUDIT_TMP" 2>/dev/null || true
tk_link_fail=0
while IFS= read -r so; do
  for ref in $(otool -L "$so" 2>/dev/null | awk 'NR>1{print $1}' | grep -E 'libtcl8\.6|libtk8\.6|tcl9|tk9' || true); do
    case "$ref" in
      @rpath/*|@loader_path/*|@executable_path/*) : ;;
      "$STAGE/runtime/"*) : ;;
      *) echo "error: $(basename "$so") links tcl/tk outside the private runtime: $ref"
         tk_link_fail=1 ;;
    esac
  done
done < "$AUDIT_TMP"
rm -f "$AUDIT_TMP"
[ "$tk_link_fail" -eq 0 ] || { echo "two Tcl/Tk copies would load in one process (ghosting) - rebuild with the section-1b Tcl/Tk 8 header prefix"; exit 1; }
echo "==> tcl/tk linkage OK (single Tk copy per process)"

# --- 5. portable launcher ---------------------------------------------------------
# The wheel ships the packages AND the model/doc data at the package root
# (site-packages), so no source tree / PYTHONPATH is needed - getTopDirectory()
# understands the installed layout.  -m keeps the entry layout-independent.
cat > "$STAGE/bin/paths.sh" <<'EOF'
#!/usr/bin/env bash
# Standalone (portable) layout: everything lives under the unpacked tree.
CCPNMR_TOP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit; pwd)"
export CCPNMR_TOP_DIR
export CONDA="${CCPNMR_TOP_DIR}"/runtime
export PYTHONHOME="${CONDA}"
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
exec "\$DIR/python" -m ${mod} "\$@"
EOF
  chmod +x "$STAGE/runtime/bin/$cmd"
done

# --- 6. README --------------------------------------------------------------------
cat > "$STAGE/README-PORTABLE.txt" <<EOF
CCPNMR Analysis ${VER} - macOS ${ARCH} standalone
=================================================

Self-contained: no system Python, no pip, no installation step.
The private CPython ${PYTAG} runtime and all dependencies (including the
compiled C extensions) are embedded in this tree.  Built for ${ARCH} -
run it on the same architecture.

Run:
  cd $(basename "$STAGE")
  ./bin/analysis                 (optionally: ./bin/analysis /path/to/project)

Host requirements:
  - macOS ${ARCH}
  - XQuartz (Tk/X11 for the GUI: brew install --cask xquartz, then
    re-login or restart so /usr/X11 is initialised)
  - nothing else is read from or written to the system (except the project
    you open/save and your default browser for Project > Summary).  Tcl/Tk
    8.6 is embedded in the runtime (runtime/lib) - a matching build-time
    tcl-tk (conda 'analysis' env / CCP_TK_PREFIX) is a BUILD-machine
    requirement only.

CLI utilities (non-GUI):
  ./runtime/bin/ccpnmr-nef import file.nef [--project-name NAME] [--force]
  ./runtime/bin/ccpnmr-nef export <project-directory> <output.nef>
EOF

# --- 7. pack ------------------------------------------------------------------------
echo "==> packing"
tar czf "dist/${NAME}.tar.gz" -C dist "$NAME"
rm -rf "$STAGE"
ls -lh "dist/${NAME}.tar.gz"
echo "==> done: dist/${NAME}.tar.gz"

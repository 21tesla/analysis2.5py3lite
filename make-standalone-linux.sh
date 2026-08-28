#!/usr/bin/env bash
# Build the Linux x86_64 "standalone" (extract-and-run) distribution:
#
#   dist/ccpnmr-<ver>-linux-x86_64-standalone.tar.gz
#
# A fresh user unpacks the tarball and runs ./bin/analysis - no system
# python, no pip, no venv: the tree embeds a private CPython 3.13 runtime
# (uv-managed python-build-standalone, relocatable) with the ccpnmr wheel
# (and all dependencies) installed into its own site-packages.  Same model
# as the macOS portable distribution described in INSTALL.md.
#
# Usage:  ./make-standalone-linux.sh   (needs: uv, gcc not required - the
#         30 C extensions are pre-built in-tree and ride in the wheel)
set -euo pipefail

# the script lives in the repo root; dirname handles both `./script` and
# absolute-path invocations
TOP="$(cd "$(dirname "$0")" && pwd)"
cd "$TOP"

# --- version + architecture ---------------------------------------------------
VER="$(grep -m1 '^version' pyproject.toml | cut -d'"' -f2)"
ARCH="x86_64"
NAME="ccpnmr-${VER}-linux-${ARCH}-standalone"
STAGE="dist/${NAME}"
PYTAG="3.13"                       # must match the wheel's cp313 tags

echo "==> standalone dist: ${NAME}"

# --- 1. fresh wheel from the working tree -------------------------------------
echo "==> building wheel (uv build --wheel)"
uv build --wheel
WHEEL="dist/ccpnmr-${VER}-cp313-cp313-linux_${ARCH}.whl"
[ -f "$WHEEL" ] || { echo "wheel not found: $WHEEL"; exit 1; }

# --- 2. private relocatable interpreter ---------------------------------------
# The venv/conda pythons on this box are NOT relocatable; the runtime must be
# a uv-managed python-build-standalone (self-contained, relative-internal).
UVPY="$(uv python find "$PYTAG" --python-preference only-managed 2>/dev/null || true)"
if [ -z "$UVPY" ]; then
  echo "==> installing managed CPython ${PYTAG} (one-time network download)"
  uv python install "$PYTAG"
  UVPY="$(uv python find "$PYTAG" --python-preference only-managed)"
fi
UVHOME="$(cd "$(dirname "$UVPY")/.." && pwd)"
echo "==> runtime: $UVHOME"
"$UVPY" -c "import _tkinter" || { echo "managed python lacks _tkinter - abort"; exit 1; }

rm -rf "$STAGE"
mkdir -p "$STAGE/bin"
cp -a "$UVHOME" "$STAGE/runtime"
chmod -R u+w "$STAGE/runtime"      # uv's pythons are r-x; pip needs to write

# --- 3. install the wheel (+ deps) into the private runtime --------------------
echo "==> installing wheel into private runtime"
uv pip install --python "$STAGE/runtime/bin/python" "$WHEEL" 2>/dev/null \
  || uv pip install --python "$STAGE/runtime/bin/python" --break-system-packages "$WHEEL"
# sanity: the package must import from INSIDE the staged runtime, full stop
# (cd / proves independence from the source tree - hence absolute paths)
(cd / && "$TOP/$STAGE/runtime/bin/python" -c "
import ccpnmr.analysis.AnalysisGui as G, ccpnmr.nefCli as N, ccp.gui.Io as Io
assert G.__file__.startswith('$TOP/$STAGE/runtime')
assert hasattr(Io, 'loadNefProject')
print('runtime import OK:', G.__file__)")

# --- 4. portable launcher -------------------------------------------------------
# The wheel ships the packages AND the model/doc data at the package root
# (site-packages), so no source tree / PYTHONPATH is needed - getTopDirectory()
# understands the installed layout.  -m keeps the entry layout-independent.
cat > "$STAGE/bin/paths.sh" <<'EOF'
#!/usr/bin/env bash
# Standalone (portable) layout: everything lives under the unpacked tree.
CCPNMR_TOP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit; pwd)"
export CCPNMR_TOP_DIR
export CONDA="${CCPNMR_TOP_DIR}"/runtime
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

cat > "$STAGE/README-PORTABLE.txt" <<EOF
CCPNMR Analysis ${VER} - Linux x86_64 standalone
================================================

Self-contained: no system Python, no pip, no installation step.
The private CPython ${PYTAG} runtime and all dependencies (including the
compiled C extensions) are embedded in this tree.

Run:
  cd $(basename "$STAGE")
  ./bin/analysis                 (optionally: ./bin/analysis /path/to/project)

Requirements on the host (a normal Linux desktop):
  - graphics: X11 libraries (libX11/libGL - present on any desktop install)
  - nothing else is read from or written to the system (except the project
    you open/save and your default browser for Project > Summary)

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

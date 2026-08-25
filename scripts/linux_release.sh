#!/usr/bin/env bash
# ============================================================================
# Linux release builder for CCPNMR Analysis 2.5.2 (Python 3).
#
# Builds the full C-extension set (memops backbone + per-package FAM exts),
# packages the wheel, and verifies it the same way the Linux wheel was
# shipped: clean venv -> pip install -> pip check -> 4 console entry points
# -> whole-tree import smoke (FAILED must be 0).
#
# Usage:
#   ./scripts/linux_release.sh                  # build + verify the wheel
#   ./scripts/linux_release.sh --release        # ...and create the GitHub Release (needs gh)
#   ./scripts/linux_release.sh --tag NAME       # custom release tag
#   PYTHON=python3.13 ./scripts/linux_release.sh    # pick the interpreter
#
# Prereqs (auto-detected, no env vars needed normally):
#   Linux, gcc (cc), Python >=3.13 with Tk/Tcl development packages,
#   X11 libraries, OpenGL/Mesa + GLUT development packages.
#   Prefix overrides: CCP_TK_PREFIX.
# ============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

RELEASE=0
TAG=""
PYTHON="${PYTHON:-python3}"
while [ $# -gt 0 ]; do
  case "$1" in
    --release) RELEASE=1 ;;
    --tag)     shift; TAG="${1:-}" ;;
    --python)  shift; PYTHON="${1:-}" ;;
    -h|--help) sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown argument: $1 (see --help)" >&2; exit 2 ;;
  esac
  shift
done

die() { echo; echo "ERROR: $*" >&2; echo; exit 1; }

# Resolve the interpreter to an absolute path — the script cd's around, and
# relative PYTHON values (e.g. .venv/bin/python) would break in subshells.
if [ -x "$PYTHON" ]; then
  case "$PYTHON" in /*) ;; *) PYTHON="$PWD/$PYTHON" ;; esac
elif command -v "$PYTHON" >/dev/null 2>&1; then
  PYTHON="$(command -v "$PYTHON")"
fi

echo "=== ccpnmr 2.5.2 (py3) Linux release builder ==="
echo "repo:     $ROOT @ $(git rev-parse --short HEAD 2>/dev/null || echo '?')"
echo "platform: $(uname -s) $(uname -m)"
echo "python:   $PYTHON"
echo

# ---------------------------------------------------------------------------
# Prerequisites
# ---------------------------------------------------------------------------
[ "$(uname -s)" = "Linux" ] && echo "ok: Linux" || echo "note: not Linux (test/CI mode; the wheel will be tagged for this platform)"
command -v cc >/dev/null || die "cc (compiler) not found — run: apt-get install build-essential (or equivalent)"

"$PYTHON" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3,13) else 1)' \
  || die "$PYTHON is not Python 3.13+ — try: PYTHON=python3.13 ./scripts/linux_release.sh"

# pip (uv-created venvs may lack it) + setuptools (build)
"$PYTHON" -m pip --version >/dev/null 2>&1 || "$PYTHON" -m ensurepip --upgrade >/dev/null 2>&1 \
  || die "pip (and ensurepip) unavailable in $PYTHON — use a standard python"
"$PYTHON" -c 'import setuptools' 2>/dev/null || "$PYTHON" -m pip install --quiet setuptools

# Tk headers — same candidate order as setup.py::_tkinc()
"$PYTHON" -c "
import os, sys, sysconfig
cands = []
if os.environ.get(\"CCP_TK_PREFIX\"):
    cands.append(os.path.join(os.environ[\"CCP_TK_PREFIX\"], \"include\", \"tcl-tk\"))
    cands.append(os.path.join(os.environ[\"CCP_TK_PREFIX\"], \"include\"))
inc = sysconfig.get_paths()[\"include\"]
cands.append(os.path.dirname(inc)) # direct parent of include dir, e.g. /conda/include
cands.append(os.path.join(os.path.dirname(inc), \"include\"))
cfg = os.path.join(os.path.dirname(os.path.dirname(sys.executable)), \"pyvenv.cfg\")
if os.path.exists(cfg):
    for line in open(cfg):
        if line.startswith(\"home = \"):
            base_dir = os.path.dirname(line.split(\"=\", 1)[1].strip())
            cands.append(base_dir)
            cands.append(os.path.join(base_dir, \"include\"))
            break
if os.environ.get(\"CONDA_PREFIX\"):
    cands.append(os.path.join(os.environ[\"CONDA_PREFIX\"], \"include\"))
cands += [
    \"/usr/include\",
    \"/usr/include/tcl-tk\",
    inc
]
hit = next((c for c in cands if os.path.exists(os.path.join(c, \"tk.h\"))), None)
if not hit:
    sys.exit(\"tk.h not found. Install Tk for Python or set CCP_TK_PREFIX=<prefix>\")
print(\"ok: tk.h in \" + hit)
" || die "tk.h not found — install Tk (see above)"

# OpenGL + X11 headers and libraries check on Linux
if [ "$(uname -s)" = "Linux" ]; then
  # Check for GL/glx.h, GL/glut.h
  [ -f "/usr/include/GL/glx.h" ] || die "GL/glx.h not found — install OpenGL headers (e.g., apt install libgl1-mesa-dev or freeglut3-dev)"
  [ -f "/usr/include/GL/glut.h" ] || die "GL/glut.h not found — install GLUT headers (e.g., apt install freeglut3-dev)"
  
  # Check for X11 libs
  X11_LIB=""
  for lib in "/usr/lib/x86_64-linux-gnu" "/usr/lib" "/usr/lib64"; do
    if [ -d "$lib" ] && ls "$lib"/libX11.so* >/dev/null 2>&1; then
      X11_LIB="$lib"; break
    fi
  done
  [ -n "$X11_LIB" ] || die "libX11.so not found — install X11 (e.g., apt install libx11-dev)"
  echo "ok: X11 libs found in $X11_LIB"
  
  # Check for GL/glut libs
  GL_LIB=""
  for lib in "/usr/lib/x86_64-linux-gnu" "/usr/lib" "/usr/lib64"; do
    if [ -d "$lib" ] && ls "$lib"/libGL.so* >/dev/null 2>&1 && ls "$lib"/libglut.so* >/dev/null 2>&1; then
      GL_LIB="$lib"; break
    fi
  done
  [ -n "$GL_LIB" ] || die "libGL.so or libglut.so not found — install OpenGL/glut (e.g., apt install libgl1-mesa-dev freeglut3-dev)"
  echo "ok: OpenGL/glut libs found in $GL_LIB"
fi

echo

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
echo "=== [1/3] building C extensions ==="
"$PYTHON" setup.py build_ext --inplace

echo "=== [2/3] wheel ==="
./scripts/copy_cext.sh
"$PYTHON" -m pip wheel . --no-deps -w dist --quiet
WHEEL="$(ls -t dist/ccpnmr-2.5.2-*.whl 2>/dev/null | head -1 || true)"
[ -n "$WHEEL" ] || die "no wheel found in dist/"
echo "wheel:    $WHEEL ($(du -h "$WHEEL" | cut -f1))"
echo

# ---------------------------------------------------------------------------
# Verify in a clean venv (the Linux-release gate)
# ---------------------------------------------------------------------------
VERIFY="$ROOT/build/release-verify-venv"
mkdir -p "$ROOT/build"
rm -rf "$VERIFY"
echo "=== [3/3] verification (clean venv: install + gates) ==="
"$PYTHON" -m venv "$VERIFY"
VP="$VERIFY/bin/python"
"$VP" -m pip install --quiet --upgrade pip
# Core distribution deps (what a naive user installs), then the optional
# web/plotting extras — the release gate (import-smoke FAILED=0) matches the
# Phase-4 installed-state gate, which ran green with the full stack.
"$VP" -m pip install --quiet "$WHEEL" numpy pandas PyOpenGL Pillow olefile requests python-dateutil pytz
"$VP" -m pip install --quiet matplotlib cherrypy mako
"$VP" -m pip check && echo "ok: pip check clean"

"$VP" - <<'EOF'
from importlib.metadata import entry_points
need = {"ccpnmr", "ccpnmr-data-shifter", "ccpnmr-format-converter", "ccpnmr-update"}
have = {e.name for e in entry_points(group="console_scripts")}
missing = need - have
assert not missing, f"missing console entry points: {sorted(missing)}"
print("ok: all 4 console entry points present")
EOF

SP="$("$VP" -c 'import site; print(site.getsitepackages()[0])')"
MPLBACKEND=Agg CCP_SMOKE_ROOT="$SP" "$VP" "$ROOT/import_smoke.py" > build/release-smoke.log 2>&1 || true
read -r OK FAILED DESIGNED < <(awk '/^  OK:/{ok=$2} /^  FAILED:/{f=$2} /^  BY-DESIGN:/{d=$2} END{print ok, f, d}' build/release-smoke.log)
echo "import-smoke (installed): OK=${OK:-?} FAILED=${FAILED:-?} BY-DESIGN=${DESIGNED:-?}   (log: build/release-smoke.log)"
[ "${FAILED:-1}" = "0" ] || { echo; tail -80 build/release-smoke.log; die "import-smoke FAILED != 0 — wheel is NOT release-ready"; }

if "$VP" -m pip install --quiet twine 2>/dev/null; then
  "$VP" -m twine check "$WHEEL" && echo "ok: twine check"
else
  echo "warn: twine check skipped (unavailable)"
fi

echo
echo "=== RELEASE READY ==="
WHEEL_BASE="$(basename "$WHEEL")"
if [ "$RELEASE" = 1 ]; then
  ARCH="$(uname -m)"
  TAG="${TAG:-v2.5.2-py3-linux-$ARCH}"
  command -v gh >/dev/null || die "gh CLI not found — install gh — or drop --release and upload $WHEEL to a GitHub Release manually"
  gh auth status >/dev/null 2>&1 || die "not authenticated — run: gh auth login"
  REPO="$(gh repo view --json nameWithOwner -q .nameWithOwner)"
  gh release create "$TAG" "$WHEEL" \
    --title "$TAG — CCPNMR Analysis 2.5.2, Python 3, Linux $ARCH" \
    --notes "CCPNMR Analysis 2.5.2 modernized to Python 3 — Linux $ARCH wheel (tag $(git rev-parse --short HEAD)).

Install:
  python3.13 -m venv ~/ccpnmr && source ~/ccpnmr/bin/activate
  pip install $WHEEL_BASE numpy pandas PyOpenGL Pillow olefile requests python-dateutil pytz
  ccpnmr"
  URL="https://github.com/$REPO/releases/download/$TAG/$WHEEL_BASE"
  echo
  echo "published: $TAG on $REPO"
  echo "naive install URL:"
  echo "  pip install \"$URL\""
else
  echo "wheel ready at: $WHEEL"
  echo "to publish:     ./scripts/linux_release.sh --release   (creates the GitHub Release + upload)"
fi

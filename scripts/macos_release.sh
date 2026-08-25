#!/usr/bin/env bash
# ============================================================================
# macOS release builder for CCPNMR Analysis 2.5.2 (Python 3).
#
# Builds the full C-extension set (memops backbone + per-package FAM exts),
# packages the wheel, and verifies it the same way the wheel was shipped:
# clean venv -> pip install -> pip check -> 4 console entry points ->
# whole-tree import smoke (FAILED must be 0).
#
# Usage:
#   ./scripts/macos_release.sh                  # build + verify the wheel
#   ./scripts/macos_release.sh --release        # ...and create the GitHub Release (needs gh)
#   ./scripts/macos_release.sh --tag NAME       # custom release tag
#   PYTHON=python3.13 ./scripts/macos_release.sh    # pick the interpreter
#
# Prereqs (auto-detected, no env vars needed normally):
#   macOS, Xcode CLT (cc), Python >=3.13 with Tk (e.g. brew install python-tk@3.13),
#   XQuartz (brew install --cask xquartz).
#   Prefix overrides: CCP_TK_PREFIX, CCP_X11_PREFIX (default /opt/X11).
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

echo "=== ccpnmr 2.5.2 (py3) release builder ==="
echo "repo:     $ROOT @ $(git rev-parse --short HEAD 2>/dev/null || echo '?')"
echo "platform: $(uname -s) $(uname -m)"
echo "python:   $PYTHON"
echo

# ---------------------------------------------------------------------------
# Prerequisites
# ---------------------------------------------------------------------------
[ "$(uname -s)" = "Darwin" ] && echo "ok: Darwin" || echo "note: not Darwin (test/CI mode; the wheel will be tagged for this platform)"
command -v cc >/dev/null || die "cc not found — run: xcode-select --install"

"$PYTHON" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3,13) else 1)' \
  || die "$PYTHON is not Python 3.13+ — try: PYTHON=python3.13 ./scripts/macos_release.sh"

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
cands.append(os.path.join(os.path.dirname(inc), \"include\"))
cfg = os.path.join(os.path.dirname(os.path.dirname(sys.executable)), \"pyvenv.cfg\")
if os.path.exists(cfg):
    for line in open(cfg):
        if line.startswith(\"home = \"):
            cands.append(os.path.join(os.path.dirname(line.split(\"=\", 1)[1].strip()), \"include\"))
            break
cands += [
    \"/opt/homebrew/opt/tcl-tk/include/tcl-tk\",
    \"/opt/homebrew/opt/tcl-tk/include\",
    \"/usr/local/opt/tcl-tk/include/tcl-tk\",
    \"/usr/local/opt/tcl-tk/include\",
    inc
]
hit = next((c for c in cands if os.path.exists(os.path.join(c, \"tk.h\"))), None)
if not hit:
    sys.exit(\"tk.h not found. Install Tk for Python (brew install python-tk@3.13, or conda: python=3.13 tk) or set CCP_TK_PREFIX=<prefix>\")
print(\"ok: tk.h in \" + hit)
" || die "tk.h not found — install Tk (see above)"

# X11 (macOS: XQuartz supplies libX11 for the Tk window handler)
if [ "$(uname -s)" = "Darwin" ]; then
  X11P="${CCP_X11_PREFIX:-/opt/X11}"
  [ -d "$X11P/lib" ] || die "X11 libs not found at $X11P — install XQuartz: https://www.xquartz.org (or brew install --cask xquartz), then log out/in"
  echo "ok: X11 at $X11P"
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
  TAG="${TAG:-v2.5.2-py3-macos-$ARCH}"
  command -v gh >/dev/null || die "gh CLI not found — brew install gh — or drop --release and upload $WHEEL to a GitHub Release manually"
  gh auth status >/dev/null 2>&1 || die "not authenticated — run: gh auth login"
  REPO="$(gh repo view --json nameWithOwner -q .nameWithOwner)"
  gh release create "$TAG" "$WHEEL" \
    --title "$TAG — CCPNMR Analysis 2.5.2, Python 3, macOS $ARCH" \
    --notes "CCPNMR Analysis 2.5.2 modernized to Python 3 — macOS $ARCH wheel (tag $(git rev-parse --short HEAD)).

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
  echo "to publish:     ./scripts/macos_release.sh --release   (creates the GitHub Release + upload)"
fi

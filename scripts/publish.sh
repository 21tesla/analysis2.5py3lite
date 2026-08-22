#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Build + verify (and optionally upload) the ccpnmr 2.5.2 (py3.13) distribution.
#
#   bash scripts/publish.sh            build + verify + `twine check` (NO upload)
#   bash scripts/publish.sh --upload   ... then `twine upload`
#
# See docs/PUBLISHING.md for the full rationale, build-dep list, and the note
# about the possibly-upstream-owned `ccpnmr` PyPI name.
# -----------------------------------------------------------------------------
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

# Use a real C compiler (the Anaconda-provided cc lacks GL/glx.h — see P4-4a).
export CC="${CC:-/usr/bin/gcc}"
export CXX="${CXX:-/usr/bin/g++}"
# Optional: Meccano/GSL. Uncomment to build the grenoble Meccano C ext.
# export CCP_GSL_PREFIX="${CCP_GSL_PREFIX:-/opt/conda/envs/ccpnmr-gsl}"

echo "==> Building sdist + wheel"
python -m build   # -> dist/ccpnmr-2.5.2-*.whl + dist/ccpnmr-2.5.2.tar.gz

echo "==> Verifying a clean install"
PUBVENV="${TMPDIR:-/tmp}/ccpnmr-pub-check"
python -m venv --clear "$PUBVENV"
"$PUBVENV/bin/pip" install --quiet --upgrade pip
"$PUBVENV/bin/pip" install --quiet dist/ccpnmr-2.5.2-*-linux_*.whl || \
  "$PUBVENV/bin/pip" install --quiet dist/ccpnmr-2.5.2-*.whl
"$PUBVENV/bin/pip" check
n_scripts="$(ls "$PUBVENV/bin" | grep -c '^ccpnmr' || true)"
echo "console scripts installed: $n_scripts (expect 8)"
[ "$n_scripts" -eq 8 ]

echo "==> twine check (metadata)"
command -v twine >/dev/null 2>&1 || python -m pip install --quiet --user twine
python -m twine check dist/ccpnmr-*

if [[ "${1:-}" == "--upload" ]]; then
  echo "==> twine upload"
  python -m twine upload ${REPO:+--repository "$REPO"} dist/ccpnmr-*
else
  echo "Build + checks complete. Re-run with --upload (and PYPI_API_TOKEN set) to publish."
fi

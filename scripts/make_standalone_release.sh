#!/usr/bin/env bash
# ============================================================================
# Creates a standalone pre-built tar.gz release of the repository.
# This builds the C and Cython extensions, then packages the files needed
# to run the application natively without needing `pip install`.
#
# Usage:
#   ./scripts/make_standalone_release.sh
#   PYTHON=python3.13 ./scripts/make_standalone_release.sh
# ============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PYTHON="${PYTHON:-python3}"
VERSION="2.5.2"
ARCH="$(uname -s | tr '[:upper:]' '[:lower:]')-$(uname -m)"
REL_NAME="ccpnmr-${VERSION}-standalone-${ARCH}"

# Resolve the interpreter to an absolute path
if [ -x "$PYTHON" ]; then
  case "$PYTHON" in /*) ;; *) PYTHON="$PWD/$PYTHON" ;; esac
elif command -v "$PYTHON" >/dev/null 2>&1; then
  PYTHON="$(command -v "$PYTHON")"
fi

echo "=== Building C extensions ==="
"$PYTHON" setup.py build_ext --inplace

CY="$ROOT/ccpnmr2.5/python/cing/Libs/cython"
echo "=== Building Cython superpose (cing) ==="
if [ -d "$CY" ] && [ -f "$CY/compile.py" ]; then
    ( cd "$CY" && "$PYTHON" compile.py build_ext --inplace )
else
    echo "Warning: Cython superpose path not found ($CY)"
fi

echo "=== Copying C extensions to target directories ==="
if [ -f "./scripts/copy_cext.sh" ]; then
    ./scripts/copy_cext.sh
else
    echo "Warning: ./scripts/copy_cext.sh not found."
fi

echo "=== Packaging standalone tar.gz ==="
TMP_DIR="$(mktemp -d)"
TARGET="${TMP_DIR}/${REL_NAME}"

mkdir -p "${TARGET}"

# Copy necessary directories and files
echo "Copying runtime files..."
cp -a bin "${TARGET}/"
cp -a ccpnmr2.5 "${TARGET}/"
[ -d doc ] && cp -a doc "${TARGET}/" || true
[ -d testproject ] && cp -a testproject "${TARGET}/" || true
[ -f database.txt ] && cp -a database.txt "${TARGET}/" || true
[ -f dbTable-new ] && cp -a dbTable-new "${TARGET}/" || true
[ -f dbTable.new ] && cp -a dbTable.new "${TARGET}/" || true
[ -f README.md ] && cp -a README.md "${TARGET}/" || true
[ -f INSTALL.md ] && cp -a INSTALL.md "${TARGET}/" || true
[ -f LICENSE ] && cp -a LICENSE "${TARGET}/" || true

# Clean up pycache and unneeded build artifacts from the package
echo "Cleaning up __pycache__ and git metadata..."
find "${TARGET}" -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
find "${TARGET}" -name "*.pyc" -delete

TARBALL="${ROOT}/${REL_NAME}.tar.gz"
echo "Creating tarball: ${REL_NAME}.tar.gz"
tar -czf "${TARBALL}" -C "${TMP_DIR}" "${REL_NAME}"

rm -rf "${TMP_DIR}"
echo "=== Done: ${TARBALL} ==="

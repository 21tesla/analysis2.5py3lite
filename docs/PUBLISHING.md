# Publishing & installing ccpnmr 2.5.2 (Python 3.13)

This documents the two things "shipped" by Phase 5:

- **How downstream users install it** (from PyPI, a local wheel, or an sdist).
- **How you publish the PyPI release** (build → verify → upload), and the
  optional conda-forge recipe starting point in [`recipe/`](./recipe/meta.yaml).

Release artifact: **`ccpnmr 2.5.2`** for **Python ≥ 3.13**, tag **`v2.5.2-py3`**.

---

## 1. Installing (for downstream users)

The verified install options, best-first:

### a) From a source-built wheel (most robust)
```bash
pip install <path>/ccpnmr-2.5.2-cp313-*-linux_x86_64.whl
# for full cing / web-server / advanced-I/O coverage (optional third-party):
pip install "ccpnmr[optional]"
```

### b) From the sdist (compiles the C extensions at install time)
```bash
pip install <path>/ccpnmr-2.5.2.tar.gz
```
This builds the C extensions locally — it needs the build dependencies below.

### c) From PyPI (once the release is published)
```bash
pip install "ccpnmr>=2.5.2"
pip install "ccpnmr[optional]"   # optional third-party (scipy, matplotlib, sqlalchemy, ...)
```

The 8 console commands become available on `PATH`:
`ccpnmr`, `ccpnmr-eci`, `ccpnmr-dangle`, `ccpnmr-data-shifter`, `ccpnmr-deposition`,
`ccpnmr-extend-nmr`, `ccpnmr-format-converter`, `ccpnmr-update`.

---

## 2. Publishing (build → verify → upload)

### 2.1 Build-Environment requirements
The C extensions need a real compiler toolchain. **Use the system gcc** (the
Anaconda-provided `cc` lacks GL / `glx.h` in its sysroot — Phase-4 lesson):

```bash
export CC=/usr/bin/gcc CXX=/usr/bin/g++
```
System packages (Linux): `build-essential`, `python3-dev`, `python3-tk`, `libgl1-mesa-dev`,
`freeglut3-dev`, `tk-dev`, `libx11-dev`. Python 3.13 headers present in the venv.

**GSL is optional** (enables the grenoble *Meccano* C extension). Without it the build
succeeds and a warning is printed; with it set `CCP_GSL_PREFIX` (conda `gsl` from the
`ccpnmr-gsl` env, or `apt install libgsl-dev`).

### 2.2 Build
```bash
CC=/usr/bin/gcc CXX=/usr/bin/g++ uv build     # -> dist/ccpnmr-2.5.2-*.whl + dist/ccpnmr-2.5.2.tar.gz
```

### 2.3 Verify (the gate set)
```bash
# clean install into a throwaway env, then confirm the distribution is self-contained
uv venv --seed --python 3.13 /tmp/ccpnmr-pub
uv pip install --python /tmp/ccpnmr-pub/bin/python dist/ccpnmr-2.5.2-cp313-*.whl
/tmp/ccpnmr-pub/bin/pip check                     # -> No broken requirements
ls /tmp/ccpnmr-pub/bin | grep -c '^ccpnmr'        # -> 8 (the console scripts)
CCP_SMOKE_ROOT=/tmp/ccpnmr-pub/lib/python3.13/site-packages MPLBACKEND=Agg \
  /tmp/ccpnmr-pub/bin/python import_smoke.py      # -> FAILED: 0
cd /tmp && /tmp/ccpnmr-pub/bin/python -m pytest --pyargs ccpnmr -q   # -> 0 failures
xvfb-run -a /tmp/ccpnmr-pub/bin/python gui_boot_test.py              # -> 8/8 booted
```

### 2.4 Upload
```bash
pip install --user twine
python -m twine check dist/ccpnmr-*                       # validate metadata
python -m twine upload --repository testpypi dist/ccpnmr-*  # dry-run against TestPyPI first
python -m twine upload dist/ccpnmr-*                       # -> pypi.org (uses $PYPI_API_TOKEN)
```
`scripts/publish.sh` bundles build → verify → (opt-in) upload:
```bash
bash scripts/publish.sh                 # build + verify + twine check (no upload)
bash scripts/publish.sh --upload        # ... then twine upload
```

> **Name caveat:** `ccpnmr` is the canonical CCPN project name and may be owned by the
> upstream project on PyPI. If the name is already claimed, publish your 3.13 port under a
> scoped index / a fork name (e.g. `ccpnmr-py3`) via `twine upload --repository-url ...`,
> or get a name transfer / re-release from the upstream maintainers. Verify with
> `pip index versions ccpnmr` before publishing.

---

## 3. Conda-forge (optional)

`recipe/` holds a **starting-point** conda-forge recipe (`meta.yaml`) that builds from the
sdist with the C compiler + numpy/cython host deps. To publish to conda-forge you must open a
`conda-forge` PR that includes `recipe/`, fill in the `sha256`, and ensure the system
GL/Tk/GSL deps are available in the conda-forge build env (see `recipe/README.md`).
GSL / Meccano remains optional there.

"""Verify all 8 C extensions import and expose their expected APIs."""
import importlib
import sys
from pathlib import Path

import pytest

TEST_FILES = [
    ("ShapeFile", "ShapeFile", ["ShapeFile"]),
    ("MemCache", "MemCache", ["MemCache"]),
    ("BlockFile", "BlockFile", ["BlockFile", "ShapeBlockFile", "compareSlices"]),
    ("FitMethod", "FitMethod", ["FitMethod", "runFit", "fitData", "bootstrapData"]),
    ("StoreFile", "StoreFile", ["StoreFile"]),
    ("StoreHandler", "StoreHandler", ["StoreHandler"]),
    ("PdfHandler", "PdfHandler", ["PdfHandler", "error"]),
    ("PsHandler", "PsHandler", ["PsHandler", "error"]),
]


def _import_cc_ext(name):
    """Import a C extension by module name, falling back to the pre-built path."""
    try:
        return importlib.import_module(name)
    except ImportError:
        # Look in ccpnmr2.5/python/ for the .so file
        lib_dir = Path(__file__).resolve().parent.parent
        for so_file in lib_dir.glob(f"{name}.cpython-*.so"):
            spec = importlib.util.spec_from_file_location(name, so_file)
            mod = importlib.util.module_from_spec(spec)
            sys.modules[name] = mod
            spec.loader.exec_module(mod)
            return mod
        raise


class TestCEstensions:
    @pytest.mark.parametrize(
        "mod_name, py_name, attrs",
        [(m, p, a) for m, p, a in TEST_FILES],
        ids=[m for m, p, a in TEST_FILES],
    )
    def test_module_imports(self, mod_name, py_name, attrs):
        mod = _import_cc_ext(mod_name)
        assert hasattr(mod, "error"), f"{mod_name} missing 'error' exception"
        for attr in attrs:
            assert hasattr(mod, attr), f"{mod_name} missing attribute '{attr}'"

    # Extensions that cannot be built as a bare headless test object: they are
    # store/block readers (need a pre-existing data file) or plotting-stream
    # handlers (need an output stream + style). They are still import/attribute
    # checked by test_module_imports above and exercised in the GUI plot path.
    _HEADLESS_UNAVAILABLE = {
        "StoreFile": "block-file reader; needs a pre-existing valid .stc data file",
        "BlockFile": "needs a pre-existing block file + a MemCache object",
        "PdfHandler": "plotting handler; needs an output stream + output_style (GUI plot path)",
        "PsHandler": "plotting handler; needs an output stream + output_style (GUI plot path)",
    }

    @pytest.mark.parametrize(
        "mod_name, py_name, attrs",
        [(m, p, a) for m, p, a in TEST_FILES],
        ids=[f"inst_{m}" for m, p, a in TEST_FILES],
    )
    def test_ext_instantiation(self, mod_name, py_name, attrs, tmp_path):
        """Instantiate every C extension that builds without external data, and
        add a cheap functional check where one exists (ShapeFile dims, FitMethod
        runFit).  Exts that need runtime data/streams skip with a specific reason
        (see _HEADLESS_UNAVAILABLE) instead of the previous blanket per-param skip,
        so all 8 parameters are either exercised (ShapeFile, MemCache, FitMethod,
        StoreHandler) or documented (StoreFile, BlockFile, PdfHandler, PsHandler).
        """
        mod = _import_cc_ext(mod_name)
        if mod_name == "ShapeFile":
            s = mod.ShapeFile(2, [10, 10])
            assert s.ndim == 2
            assert s.ncomponents == 2
        elif mod_name == "MemCache":
            assert mod.MemCache(2) is not None
        elif mod_name == "FitMethod":
            assert mod.FitMethod(1, 0.1) is not None  # LINEAR_FIT = 1
            params, params_dev, y_fit, chisq = mod.runFit(1, 0.1, [0, 1, 2, 3], [1, 2, 4, 7])
            assert len(params) == 2  # slope + intercept
            assert chisq >= 0
        elif mod_name == "StoreHandler":
            assert mod.StoreHandler(str(tmp_path / "store.stc")) is not None
        else:
            pytest.skip(self._HEADLESS_UNAVAILABLE.get(mod_name, "requires external data/streams"))

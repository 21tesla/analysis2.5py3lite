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

    @pytest.mark.parametrize(
        "mod_name, py_name, attrs",
        [(m, p, a) for m, p, a in TEST_FILES],
        ids=[f"inst_{m}" for m, p, a in TEST_FILES],
    )
    def test_shape_file_instantiation(self, mod_name, py_name, attrs):
        if mod_name != "ShapeFile":
            pytest.skip("only ShapeFile gets full instantiation check here")
        mod = _import_cc_ext(mod_name)
        s = mod.ShapeFile(2, [10, 10])
        assert s.ndim == 2
        assert s.ncomponents == 2

    @pytest.mark.parametrize(
        "mod_name, py_name, attrs",
        [(m, p, a) for m, p, a in TEST_FILES],
        ids=[f"fit_{m}" for m, p, a in TEST_FILES],
    )
    def test_fit_method_instantiation(self, mod_name, py_name, attrs):
        if mod_name != "FitMethod":
            pytest.skip("only FitMethod gets instantiation check here")
        mod = _import_cc_ext(mod_name)
        # LINEAR_FIT = 1
        f = mod.FitMethod(1, 0.1)
        assert f is not None
        # runFit should return a 4-tuple
        params, params_dev, y_fit, chisq = mod.runFit(1, 0.1, [0, 1, 2, 3], [1, 2, 4, 7])
        assert len(params) == 2  # slope + intercept
        assert chisq >= 0

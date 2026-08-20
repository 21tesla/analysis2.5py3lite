"""Smoke tests for memops.general.Util and helper modules."""
from memops.general.Constants import currentModelVersion


class TestConstants:
    def test_ccpn_version_exists(self):
        assert currentModelVersion is not None
        assert len(str(currentModelVersion)) > 0
        # Version should be a dotted numeric string (e.g. "2.1.2")
        parts = str(currentModelVersion).split(".")
        assert len(parts) >= 2, f"Unexpected version format: {currentModelVersion}"
        assert all(p.isdigit() for p in parts[:2]), f"Version has non-numeric component: {currentModelVersion}"


class TestUtilFunctions:
    def test_util_module_imports(self):
        import memops.general.Util as Util

        # These are the core helper functions that must be importable
        assert hasattr(Util, "copySubTree") or hasattr(Util, "setLinks") or hasattr(Util, "getLinkData")

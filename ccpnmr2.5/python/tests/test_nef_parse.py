"""P4-2 functional test: NEF text parsing (pure Python, no C exts required).

Verifies that the bundled NEF test file can be parsed and that the structural
objects (DataBlock, SaveFrame, Loop) expose the expected content.
"""
from pathlib import Path

import ccpnmr
from ccpnmr.nef import StarIo

# Anchored at the ccpnmr package: identical for source tree and installed dist.
NEF = Path(ccpnmr.__file__).resolve().parent / "nef" / "testdata" / "CCPN_Commented_Example.nef"


def _parse():
    return StarIo.parseNefFile(str(NEF))


def _block(extent):
    return list(extent.values())[0]


class TestNefParse:
    def test_single_data_block(self):
        nef = _parse()
        blocks = list(nef.values())
        assert len(blocks) == 1, f"expected 1 data block, got {len(blocks)}"

    def test_data_block_name(self):
        block = _block(_parse())
        assert block.name.startswith("nef_my_nmr_project"), f"got {block.name!r}"

    def test_meta_saveframe(self):
        block = _block(_parse())
        assert "nef_nmr_meta_data" in block, f"missing meta-saveframe; keys={list(block.keys())}"
        meta = block["nef_nmr_meta_data"]
        assert meta.get("sf_category") == "nef_nmr_meta_data"
        fmt = meta.get("format_name")
        assert fmt in ("nmr_exchange_format", "nif"), f"format_name={fmt!r}"

    def test_molecular_system_chains(self):
        block = _block(_parse())
        assert "nef_molecular_system" in block, f"keys={list(block.keys())}"
        mol_sf = block["nef_molecular_system"]
        assert "nef_sequence" in mol_sf, f"keys={list(mol_sf.keys())}"
        seq_loop = mol_sf["nef_sequence"]
        assert seq_loop.data, "nef_sequence loop is empty"
        chain_codes = set(row["chain_code"] for row in seq_loop.data)
        assert "A" in chain_codes, f"chain A missing; codes={chain_codes}"

    def test_mol_system_residue_count(self):
        block = _block(_parse())
        seq_loop = block["nef_molecular_system"]["nef_sequence"]
        assert len(seq_loop.data) >= 40, f"expected >= 40 residues, got {len(seq_loop.data)}"

    def test_chemical_shift_list_present(self):
        block = _block(_parse())
        cs_frames = [
            k
            for k, v in block.items()
            if hasattr(v, "category") and v.category and "chemical_shift" in v.category
        ]
        assert cs_frames, f"no chemical-shift saveframe; keys={list(block.keys())}"

    def test_chemical_shift_rows(self):
        block = _block(_parse())
        cs_saveframe = None
        for k, v in block.items():
            if hasattr(v, "category") and v.category and "chemical_shift" in v.category:
                cs_saveframe = v
                break
        if cs_saveframe is None:
            return  # already covered by test_chemical_shift_list_present
        loop_names = [k for k in cs_saveframe if hasattr(cs_saveframe[k], "data")]
        assert loop_names, f"no loops in chemical-shift frame; keys={list(cs_saveframe.keys())}"
        loop = cs_saveframe[loop_names[0]]
        assert len(loop.data) >= 5, f"expected >= 5 CS rows, got {len(loop.data)}"

    def test_saveframe_categories_all_present(self):
        """All required NEF saveframe categories appear in the parsed structure."""
        block = _block(_parse())
        cats = set(
            v.category for v in block.values() if hasattr(v, "category") and v.category
        )
        expected = {
            "nef_nmr_meta_data",
            "nef_molecular_system",
            "nef_chemical_shift_list",
        }
        # NOTE: category prefixes are checked with `in` because
        # categories embed a per-saveframe suffix (e.g. nef_chemical_shift_list_1)
        matched = {
            c for c in cats
            if any(e in c for e in expected)
        }
        for e in expected:
            assert matched and any(e in m for m in matched), (
                f"category {e!r} not found in parsed categories: {cats}"
            )
    def test_distance_restraint_saveframe(self):
        block = _block(_parse())
        dr_frames = [
            k
            for k, v in block.items()
            if hasattr(v, "category") and v.category and "distance_restraint" in v.category
        ]
        assert dr_frames, f"no distance-restraint frame; keys={list(block.keys())}"

    def test_dihedral_restraint_saveframe(self):
        block = _block(_parse())
        dih_frames = [
            k
            for k, v in block.items()
            if hasattr(v, "category") and v.category and "dihedral_restraint" in v.category
        ]
        assert dih_frames, f"no dihedral-restraint frame; keys={list(block.keys())}"

    def test_meta_format_version(self):
        block = _block(_parse())
        meta = block["nef_nmr_meta_data"]
        ver = meta.get("format_version")
        assert ver is not None, "format_version not in meta saveframe"
        # Converter coerces numeric-looking values to int/float; 1.1 → 1.1 (float)
        assert str(ver).startswith("1."), f"format_version={ver!r}"

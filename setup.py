"""Build script for CCPNMR C extensions (Python 3.13).

Two families of C extension modules live in this tree:

  1. The memops/MOPS "data backbone" (8 exts) — already migrated in Phase 1b.
     Imported flat/top-level (e.g. ``import ShapeFile``) and built into
     ``ccpnmr2.5/python/``.

  2. The per-package analysis exts (clouds / dynamics / analysis / ccp-structure /
     cambridge) — imported as *package submodules* (``ccpnmr.c.PeakList``,
     ``ccp.c.StructAtom``, ``cambridge.c.BayesPeakSeparator``).  Each one is
     exposed through a symlink ``<pkg>/c/<Name>.so -> c/<...>/<Name>.so``.
     We therefore build a flat ``<Name>`` extension and copy the resulting
     shared object onto that symlink target (overwriting the stale Py2 build).

Build filter:
    CCP_EXT=ShapeFile,PeekList python setup.py build_ext --inplace
    (comma-separated module names; defaults to building everything defined
    below).  This lets you build a single extension for incremental work.
"""
import os
import sys
import sysconfig

from setuptools import setup, Extension

# Source directories -----------------------------------------------------------------
G     = "ccpnmr2.5/c/memops/global"          # shared helpers + the backbone exts
CLOUD = "ccpnmr2.5/c/ccpnmr/clouds"
DYN   = "ccpnmr2.5/c/ccpnmr/dynamics"
ANA   = "ccpnmr2.5/c/ccpnmr/analysis"
STR   = "ccpnmr2.5/c/ccp/structure"
BAYES = "ccpnmr2.5/c/other/cambridge/bayes"
MEC   = "ccpnmr2.5/c/other/meccano"

# GSL (GNU Scientific Library) prefix for the Meccano ext.  Meccano is
# OPTIONAL: if no usable GSL is found the whole distribution still builds
# and Meccano is omitted with a warning (grenoble.c.Meccano import then fails
# with an actionable message).  Resolution order: $CCP_GSL_PREFIX, /usr,
# isolated conda env (conda create -n ccpnmr-gsl -c conda-forge gsl).
def _gsl_usable(prefix):
    if not prefix or not os.path.isdir(os.path.join(prefix, "include", "gsl")):
        return False
    libdir = os.path.join(prefix, "lib")
    return os.path.isdir(libdir) and any(n.startswith("libgsl") for n in os.listdir(libdir))

GSL = next(
    (p for p in (os.environ.get("CCP_GSL_PREFIX"), "/usr",
                 "/home/logan/software/anaconda3/envs/ccpnmr-gsl",
                 os.environ.get("CONDA_PREFIX"),
                 "/opt/homebrew/opt/gsl",   # Homebrew (Apple silicon)
                 "/usr/local/opt/gsl")      # Homebrew (Intel)
     if _gsl_usable(p)),
    None,
)

DARWIN = sys.platform == "darwin"


def _venv_base_prefix():
    """For venvs: the base interpreter prefix (pyvenv.cfg `home =`), or None."""
    cfg = os.path.join(os.path.dirname(os.path.dirname(sys.executable)), "pyvenv.cfg")
    if os.path.exists(cfg):
        with open(cfg) as f:
            for line in f:
                if line.startswith("home = "):
                    return os.path.dirname(line.split("=", 1)[1].strip())
    return None


def _tkinc():
    """Find a dir containing tk.h, in preference order: $CCP_TK_PREFIX/include,
    the conda-style <prefix>/include of the running (or venv BASE) interpreter
    (conda envs ship tk.h there), Homebrew tcl-tk, then the raw interpreter
    include dir as last resort."""
    cands = []
    if os.environ.get("CCP_TK_PREFIX"):
        cands.append(os.path.join(os.environ["CCP_TK_PREFIX"], "include", "tcl-tk"))
        cands.append(os.path.join(os.environ["CCP_TK_PREFIX"], "include"))
    inc = sysconfig.get_paths()["include"]
    cands.append(os.path.join(os.path.dirname(inc), "include"))   # conda layout
    base = _venv_base_prefix()
    if base:
        cands.append(os.path.join(base, "include"))               # venv -> base python
    cands += [
        "/opt/homebrew/opt/tcl-tk/include/tcl-tk",
        "/opt/homebrew/opt/tcl-tk/include",                       # Homebrew (Apple silicon)
        "/usr/local/opt/tcl-tk/include/tcl-tk",
        "/usr/local/opt/tcl-tk/include",                          # Homebrew (Intel)
        inc
    ]
    for c in cands:
        if os.path.exists(os.path.join(c, "tk.h")):
            return c
    return cands[-1]  # best guess; the compile error will name the missing header


TKINC = _tkinc()

# GL-context support differs per platform:
#  * Linux: OpenGL via GLX from an X11 GL stack (freeglut/mesa provide glx.h);
#    X11 symbols resolved via a versioned direct link (no -dev symlink needed).
#  * macOS: there is no GLX — gl_handler.c compiles its GLX context code out
#    under IGNORE_GL (the guard wraps essentially the whole file), so the
#    GL-dependent window handlers degrade to the Tk path while the data layer,
#    fitting and 2D drawing keep full function.  tk_handler.c still calls
#    Xlib directly (XDrawLine, XFillArc, ...) — XQuartz (/opt/X11) supplies
#    those headers + libX11.  GLUT is then not needed at all.
if DARWIN:
    GLX_DEFINE = (("IGNORE_GL", "1"),)
    _x11p = os.environ.get("CCP_X11_PREFIX", "/opt/X11")
    GLX_INC = [os.path.join(_x11p, "include")]
    # Library search paths (must be -L library_dirs, emitted BEFORE the -ltk8.6
    # etc. references — macOS ld resolves -l left-to-right).  Tk/Tcl live in
    # the prefix that supplied tk.h (conda env, Homebrew tcl-tk, or
    # $CCP_TK_PREFIX); X11 in XQuartz.
    if TKINC.endswith("/include/tcl-tk"):
        _tk_base = TKINC[:-15]
    elif TKINC.endswith("/include"):
        _tk_base = TKINC[:-8]
    else:
        _tk_base = os.path.dirname(TKINC)
    GLX_LIBDIRS = [os.path.join(_tk_base, "lib"),
                   os.path.join(_x11p, "lib")]
    GLX_LINK = ["-lX11"]
else:
    GLX_DEFINE = ()
    GLX_INC = []
    GLX_LIBDIRS = []
    GLX_LINK = ["-l:libX11.so.6"]

CFLAGS = ["-Wall", "-Wno-unused-function", "-Wno-unused-variable", "-Wno-error=incompatible-function-pointer-types"]


def mk(name, sources, include, libs=(), libdirs=(), define=(), link=()):
    """Build an Extension with explicit source paths + include/link settings."""
    return Extension(
        name,
        sources=sources,
        include_dirs=include,
        define_macros=list(define),
        libraries=list(libs),
        library_dirs=list(libdirs),
        extra_compile_args=CFLAGS,
        extra_link_args=list(link),
    )


# Reusable source groups (memops/global helpers that are already Py3-ported) -----
GU   = [f"{G}/utility.c", f"{G}/python_util.c"]                 # util + py_util
GMEM = [f"{G}/hash_list.c", f"{G}/hash_table.c", f"{G}/mem_cache.c",
        f"{G}/mutex.c", f"{G}/py_mem_cache.c"]                   # +GU below
GBLK = GMEM + [f"{G}/block_file.c", f"{G}/shape_file.c", f"{G}/int_array.c",
               f"{G}/py_block_file.c", f"{G}/py_shape_file.c"]   # block-file I/O

# Shared drawing/IO dep set.  py_draw_handler wraps every handler backend, so
# pull in the full handler chain (gl/tk/pdf/ps/store) with cores:
DRAWDEPS = [f"{G}/py_draw_handler.c",
            f"{G}/py_store_handler.c", f"{G}/store_handler.c",
            f"{G}/py_pdf_handler.c", f"{G}/pdf_handler.c",
            f"{G}/py_ps_handler.c", f"{G}/ps_handler.c",
            f"{G}/py_gl_handler.c", f"{G}/gl_handler.c",
            f"{G}/py_tk_handler.c", f"{G}/tk_handler.c", f"{G}/py_tk_util.c",
            f"{G}/clipping.c"]
if DARWIN:
    # IGNORE_GL removes every GL/glut reference; macOS has no GLUT by default.
    DRAWLIBS = ["tk8.6", "tcl8.6", "m"]
else:
    DRAWLIBS = ["GL", "glut", "tk8.6", "tcl8.6", "m"]
DRAWINC_EXTRA = [TKINC] + GLX_INC

# ------------------------------------------------------------------ family defs
# name -> (sources, include_dirs, libs).  Tier-1 = no GL/Tk/X11.
FAM = {
    # --- clouds (import:  ccpnmr.c.<Name>) ----------------------------------
    "CloudUtil":       ([f"{CLOUD}/py_cloud_util.c"], [CLOUD, G], ["m"]),
    "AtomCoord":       ([f"{CLOUD}/py_atom_coord.c", f"{CLOUD}/atom_coord.c"] + GU,
                        [CLOUD, G], []),
    "AtomCoordList":   ([f"{CLOUD}/py_atom_coord_list.c", f"{CLOUD}/atom_coord_list.c",
                         f"{CLOUD}/py_atom_coord.c", f"{CLOUD}/atom_coord.c"] + GU,
                        [CLOUD, G], []),
    "DistConstraint":  ([f"{CLOUD}/py_dist_constraint.c", f"{CLOUD}/dist_constraint.c"] + GU,
                        [CLOUD, G], []),
    "DistConstraintList": ([f"{CLOUD}/py_dist_constraint_list.c",
                            f"{CLOUD}/dist_constraint_list.c",
                            f"{CLOUD}/py_dist_constraint.c", f"{CLOUD}/dist_constraint.c"] + GU,
                        [CLOUD, G], []),
    "DistForce":       ([f"{CLOUD}/py_dist_force.c", f"{CLOUD}/dist_force.c"] + GU,
                        [CLOUD, G], []),
    "Dynamics":        ([f"{CLOUD}/py_dynamics.c", f"{CLOUD}/dynamics.c",
                         f"{CLOUD}/py_atom_coord_list.c", f"{CLOUD}/atom_coord_list.c",
                         f"{CLOUD}/py_atom_coord.c", f"{CLOUD}/atom_coord.c",
                         f"{CLOUD}/py_dist_constraint_list.c", f"{CLOUD}/dist_constraint_list.c",
                         f"{CLOUD}/py_dist_constraint.c", f"{CLOUD}/dist_constraint.c",
                         f"{CLOUD}/py_dist_force.c", f"{CLOUD}/dist_force.c"]
                        + GU + [f"{G}/random.c"], [CLOUD, G], ["m"]),
    "Midge":           ([f"{CLOUD}/py_midge.c", f"{CLOUD}/midge.c"]
                        + GU + [f"{G}/diag_dbl.c"], [CLOUD, G], ["m"]),
    "Bacus":           ([f"{CLOUD}/py_bacus.c", f"{CLOUD}/bacus.c"] + GU,
                        [CLOUD, G], ["m"]),

    # --- dynamics (import:  ccpnmr.c.<Name>) --------------------------------
    "DyAtomCoord":     ([f"{DYN}/py_atom_coord.c", f"{DYN}/atom_coord.c"] + GU,
                        [DYN, G], []),
    "DyAtomCoordList": ([f"{DYN}/py_atom_coord_list.c", f"{DYN}/atom_coord_list.c",
                         f"{DYN}/py_atom_coord.c", f"{DYN}/atom_coord.c"] + GU,
                        [DYN, G], []),
    "DyDistConstraint": ([f"{DYN}/py_dist_constraint.c", f"{DYN}/dist_constraint.c"] + GU,
                        [DYN, G], []),
    "DyDistConstraintList": ([f"{DYN}/py_dist_constraint_list.c",
                            f"{DYN}/dist_constraint_list.c",
                            f"{DYN}/py_dist_constraint.c", f"{DYN}/dist_constraint.c"] + GU,
                        [DYN, G], []),
    "DyDistForce":     ([f"{DYN}/py_dist_force.c", f"{DYN}/dist_force.c"] + GU,
                        [DYN, G], []),
    "DyDynamics":      ([f"{DYN}/py_dynamics.c", f"{DYN}/dynamics.c",
                         f"{DYN}/py_atom_coord_list.c", f"{DYN}/atom_coord_list.c",
                         f"{DYN}/py_atom_coord.c", f"{DYN}/atom_coord.c",
                         f"{DYN}/py_dist_constraint_list.c", f"{DYN}/dist_constraint_list.c",
                         f"{DYN}/py_dist_constraint.c", f"{DYN}/dist_constraint.c",
                         f"{DYN}/py_dist_force.c", f"{DYN}/dist_force.c"]
                        + GU + [f"{G}/random.c"], [DYN, G], ["m"]),

    # --- analysis, Tier-1 (import:  ccpnmr.c.<Name>) ------------------------
    "ContourLevels":   ([f"{ANA}/py_contour_levels.c", f"{ANA}/contour_levels.c"] + GU,
                        [ANA, G], []),
    "ContourStyle":    ([f"{ANA}/py_contour_style.c", f"{ANA}/contour_style.c"] + GU,
                        [ANA, G], []),
    "PeakList":        ([f"{ANA}/method.c", f"{ANA}/peak.c", f"{ANA}/peak_list.c",
                         f"{ANA}/symbol.c", f"{ANA}/py_peak.c", f"{ANA}/py_peak_list.c"]
                        + GU + GBLK + [f"{G}/nonlinear_model.c", f"{G}/gauss_jordan.c"],
                        [ANA, G], ["m"]),

    # --- analysis, Py3-migrated Phase 4 (import:  ccpnmr.c.<Name>) ----------
    "WinPeakList":     ([f"{ANA}/py_win_peak_list.c", f"{ANA}/win_peak_list.c",
                         f"{ANA}/method.c", f"{ANA}/peak.c", f"{ANA}/peak_list.c",
                         f"{ANA}/symbol.c", f"{ANA}/py_peak.c", f"{ANA}/py_peak_list.c"]
                        + DRAWDEPS + [f"{G}/nonlinear_model.c", f"{G}/gauss_jordan.c"] + GU + GBLK,
                        [ANA, G] + DRAWINC_EXTRA, DRAWLIBS, GLX_LIBDIRS, GLX_DEFINE, GLX_LINK),
    "PeakCluster":     ([f"{ANA}/py_peak_cluster.c", f"{ANA}/peak_cluster.c",
                         f"{ANA}/method.c", f"{ANA}/peak.c", f"{ANA}/symbol.c",
                         f"{ANA}/py_peak.c", f"{G}/nonlinear_model.c",
                         f"{G}/gauss_jordan.c"]
                        + DRAWDEPS + GU + GBLK, [ANA, G] + DRAWINC_EXTRA, DRAWLIBS,
                        GLX_LIBDIRS, GLX_DEFINE, GLX_LINK),
    "ContourFile":     ([f"{ANA}/py_contour_file.c", f"{ANA}/contour_file.c",
                         f"{ANA}/contour_data.c", f"{ANA}/contour_levels.c",
                         f"{ANA}/contour_style.c", f"{ANA}/py_contour_levels.c",
                         f"{ANA}/py_contour_style.c"]
                        + DRAWDEPS + GU + GBLK + [f"{G}/store_file.c", f"{G}/py_store_file.c",
                                                   f"{G}/contourer.c"],
                        [ANA, G] + DRAWINC_EXTRA, DRAWLIBS, GLX_LIBDIRS, GLX_DEFINE, GLX_LINK),
    "SliceFile":       ([f"{ANA}/py_slice_file.c", f"{ANA}/slice_file.c"]
                        + DRAWDEPS + GU + GBLK,
                        [ANA, G] + DRAWINC_EXTRA, DRAWLIBS, GLX_LIBDIRS, GLX_DEFINE, GLX_LINK),

    # --- ccp structure, Tier-1 (import:  ccp.c.<Name>) ----------------------
    "StructAtom":      ([f"{STR}/py_atom.c", f"{STR}/atom.c", f"{STR}/bond.c"]
                        + [f"{G}/color.c"] + GU, [STR, G], []),
    "StructBond":      ([f"{STR}/py_bond.c", f"{STR}/py_atom.c", f"{STR}/atom.c",
                         f"{STR}/bond.c"] + [f"{G}/color.c"] + GU, [STR, G], []),
    "StructUtil":      ([f"{STR}/py_struct_util.c", f"{STR}/struct_util.c"]
                        + [f"{G}/geometry.c", f"{G}/eigenvalue.c", f"{G}/linalg.c"] + GU,
                        [STR, G], ["m"]),
    "StructStructure": ([f"{STR}/py_structure.c", f"{STR}/structure.c",
                         f"{STR}/atom.c", f"{STR}/bond.c", f"{STR}/struct_util.c",
                         f"{STR}/py_atom.c", f"{STR}/py_bond.c"]
                        + DRAWDEPS + [f"{G}/color.c", f"{G}/geometry.c", f"{G}/eigenvalue.c",
                                      f"{G}/linalg.c", f"{G}/sorts.c"] + GU,
                        [STR, G] + DRAWINC_EXTRA, DRAWLIBS, GLX_LIBDIRS, GLX_DEFINE, GLX_LINK),

    # --- cambridge bayes (import:  cambridge.c.BayesPeakSeparator) ----------
    "BayesPeakSeparator": ([f"{BAYES}/py_bayes.c", f"{BAYES}/bayes_nmr.c",
                        f"{BAYES}/app.c", f"{BAYES}/distribution.c", f"{BAYES}/random.c",
                        f"{BAYES}/hilbert.c", f"{BAYES}/bayesys3.c"]
                        + GU + GBLK, [BAYES, G], ["m"]),

    # --- memops window handlers (Tier-2: needs Tk/TK headers + GL/glut/X11) --
    # import:  memops.c.GlHandler / memops.c.TkHandler
    "GlHandler":       ([f"{G}/py_gl_handler.c", f"{G}/gl_handler.c", f"{G}/py_tk_util.c",
                         f"{G}/clipping.c"]
                        + GU, [G, TKINC] + GLX_INC, DRAWLIBS,
                        GLX_LIBDIRS, GLX_DEFINE, GLX_LINK),
    "TkHandler":       ([f"{G}/py_tk_handler.c", f"{G}/tk_handler.c", f"{G}/py_tk_util.c",
                         f"{G}/clipping.c"]
                        + GU, [G, TKINC] + GLX_INC, DRAWLIBS,
                        GLX_LIBDIRS, GLX_DEFINE, GLX_LINK),

    # --- grenoble Meccano (import:  grenoble.c.Meccano; needs GSL) ----------
    "Meccano":         ([f"{MEC}/pysrc/py_meccano.c",
                        f"{MEC}/meccano2_stat_ramaDB_fwd/meccano2_stat_ramaDB_fwd.c",
                        f"{MEC}/src/myRDC.c", f"{MEC}/src/myGEOMETRY.c",
                        f"{MEC}/src/mySTRUCT.c", f"{MEC}/src/myDAT.c",
                        f"{MEC}/src/myRAMACHANDRAN.c", f"{MEC}/src/myHBSC.c",
                        f"{MEC}/src/myJUNC.c", f"{MEC}/src/myMINIMISATION.c",
                        f"{MEC}/src/myPPGEO.c", f"{MEC}/src/gradient.c"]
                        + GU, [f"{MEC}/inc", G, f"{GSL}/include"],
                        ["gsl", "gslcblas", "m"], [f"{GSL}/lib"], (),
                        [f"-Wl,-rpath,{GSL}/lib"]),
}

# ---------------------------------------------------------------------------
# The 8 backbone extensions (Phase 1b) — kept as-is, imported top-level.
def ext(name, sources, extra_sources=None):
    srcs = [os.path.join(G, s) for s in sources]
    if extra_sources:
        srcs += [os.path.join(G, s) for s in extra_sources]
    return Extension(
        name, sources=srcs, include_dirs=[G], extra_compile_args=CFLAGS)


BACKBONE = [
    ext("ShapeFile", ["py_shape_file.c", "shape_file.c", "python_util.c", "utility.c"]),
    ext("MemCache", ["py_mem_cache.c", "mem_cache.c", "hash_list.c", "hash_table.c",
                     "int_array.c", "list.c", "mutex.c", "python_util.c", "utility.c"]),
    ext("BlockFile", ["py_block_file.c", "block_file.c", "py_mem_cache.c", "py_shape_file.c",
                      "hash_list.c", "hash_table.c", "int_array.c", "list.c", "mutex.c",
                      "mem_cache.c", "shape_file.c", "python_util.c", "utility.c"]),
    ext("FitMethod", ["py_fit.c", "fit.c", "fit1d.c", "nonlinear_model.c", "cpmg.c",
                      "line_fit.c", "random.c", "gauss_jordan.c", "gamma.c",
                      "python_util.c", "utility.c"]),
    ext("StoreFile", ["py_store_file.c", "store_file.c", "python_util.c", "utility.c"]),
    ext("StoreHandler", ["py_store_handler.c", "store_handler.c", "python_util.c", "utility.c"]),
    ext("PdfHandler", ["py_pdf_handler.c", "pdf_handler.c", "clipping.c",
                       "python_util.c", "utility.c"]),
    ext("PsHandler", ["py_ps_handler.c", "ps_handler.c", "clipping.c",
                      "python_util.c", "utility.c"]),
]

# ---------------------------------------------------------------------------
all_exts = list(BACKBONE)
# Meccano is the one OPTIONAL ext: needs GSL (see resolution above).
if "Meccano" in FAM and GSL is None:
    print("WARNING: no usable GSL found ($CCP_GSL_PREFIX, /usr, conda envs) — "
          "skipping OPTIONAL ext 'Meccano'. The distribution builds without it; "
          "grenoble.c.Meccano imports fail with an actionable hint. Install GSL "
          "(e.g. conda install -c conda-forge gsl / apt install libgsl-dev) and "
          "rebuild to enable Meccano.", file=sys.stderr)
    FAM = {n: s for n, s in FAM.items() if n != "Meccano"}
# FAM specs are (srcs, inc, libs) or extended (srcs, inc, libs, libdirs, define, link)
all_exts += [mk(name, *spec) for name, spec in FAM.items()]

# Build filter (CCP_EXT=Name1,Name2) ----------------------------------------
_filter = os.environ.get("CCP_EXT", "").strip()
if _filter:
    keep = {n.strip() for n in _filter.split(",") if n.strip()}
    all_exts = [e for e in all_exts if e.name in keep]

setup(
    name="ccpnmr-ext",
    version="2.5.2",
    ext_modules=all_exts,
    zip_safe=False,
)

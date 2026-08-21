# Phase 2a — CHECKPOINT / RESUME MAP (durable, keeps working across session resets)

> **STATUS: ✅ COMPLETE (2026-08-21).** All 17 in-scope binding files DONE (rows 1–17).
> All 22 FAM exts + 8 backbone built (`setup.py build_ext --inplace`) and copied onto
> the package `c/` symlink targets. 22/22 import OK, functional smoke pass,
> `import_smoke` 1556→1561 OK (+5). No remaining PENDING rows — do NOT redo on reset.

Goal: migrate the remaining Python-2 C bindings so all 22 `setup.py` FAM exts build
under Python 3.13, then build + import-test. Work is done **file by file, by the
main session (NO sub-agents — reset-safe)**. This file is updated after every
completed file so a reset always lands on a clean file boundary.

Repo root: `/home/logan/software/ccpnmr2.5.2-qwen`
Recipe: `_phase2a_recipe.md` (canonical py2→py3 binding recipe — VERIFIED)
References (already py3): `ccpnmr2.5/c/ccpnmr/clouds/py_midge.c` (type-based),
`py_atom_coord.c` (type w/ real get/set), `py_cloud_util.c` (function-only).

## 0. HOW TO RESUME after a reset (do this first)
1. Read THIS file. Then re-establish which files are DONE vs PENDING by GROUND
   TRUTH (does not trust the table below if in doubt):
   ```
   cd /home/logan/software/ccpnmr2.5.2-qwen
   # A file is DONE iff it has a py3 init AND compiles. PENDING iff it still has py2 init.
   for f in ccpnmr2.5/c/ccpnmr/clouds/*.c ccpnmr2.5/c/ccpnmr/dynamics/*.c \
            ccpnmr2.5/c/ccpnmr/analysis/*.c ccpnmr2.5/c/ccp/structure/*.c \
            ccpnmr2.5/c/other/cambridge/bayes/*.c; do
     p=$(grep -c "PyInit_" "$f" 2>/dev/null); q=$(grep -cE "Py_InitModule|PY_MOD_INIT_FUNC" "$f" 2>/dev/null)
     printf "%-70s PyInit=%s py2init=%s\n" "$f" "$p" "$q"
   done
   ```
   (Ignore non-binding files like `midge.c`, `atom_coord.c` etc. — they have no
   `Py_InitModule` or `PyInit_` and are plain C impls; only `py_*.c` bindings matter.)
2. Pick the first PENDING `py_*.c` binding file, follow **§3** to migrate it,
   self-check per **§4**, then flip its row to DONE in **§2**.
3. When all 17 §2 rows are DONE → run **§5 build** + **§6 import test**, then §7.

## 1. CONSTRAINTS (never break these)
- Do NOT modify `setup.py`; do NOT run `setup.py build_ext` until the final §5 build (parent only).
- Do NOT build or copy any `.so` during per-file work.
- Do NOT edit shared helpers `memops/global/{python_util.c,utility.c,python_util.h,utility.h}` (already py3).
- Keep license headers, comments, and unrelated C logic intact — change only the py2→py3 binding parts.
- SAFETY NET: if a file hits a py2 token not in the recipe (e.g. `Py_FindMethod` used for
  method-lookup INSIDE a body, `PyList2DoubleSeq`, `x_item`/`x_slice`, 13-field PyTypeObject),
  STOP that file, mark it FAILED in §2 with the file:line+token — do NOT guess.
- PY2→PY3 API RENAMES (do these, NOT a stop): `PyInt_AsLong`→`PyLong_AsLong` (the exact
  "undefined symbol" that broke the stale .so — appears in setattr/getattr bodies). Real
  `print_py_X(FILE*,...fprintf...)` → `repr_py_X(PyObject*)`: `char buf[N]; snprintf(buf,...);
  return PyUnicode_FromString(buf);` and wire `(reprfunc) repr_py_X` into tp_repr. If the
  `print_py_X` is already commented-out, set tp_repr=0 instead. RELIABLE 4-step when print is REAL:
  (a) rename signature line; (b) `fprintf(fp,`→`snprintf(buf, sizeof buf,` (unique token, NO context);
  (c) prepend `char buf[N];` by anchoring on `snprintf(buf, sizeof buf, "<first-chars>`; (d) the fn's
  trailing `return 0;`→`return PyUnicode_FromString(buf);`. ⚠ These C files put a TRAILING SPACE on the
  blank lines, so multi-line anchors like `decl\n\n  fprintf` FAIL — use the single-token anchors above.

## 2. FILE STATUS TABLE (update Status column as you go)
Legend kind: T = type-based (steps 1–7), F = function-only (step 7 only).
Module name + methods-table id MUST be read from each file's own `Py_InitModule("MODULE", TABLE)`.

| # | File (relative to repo root) | MODULE | TABLE | kind | Status |
|---|---|---|---|---|---|
| 1 | ccpnmr2.5/c/ccpnmr/clouds/py_atom_coord_list.c | AtomCoordList | tp=py_handler_methods / mod=Atom_coord_list_type_methods | T | ✅ DONE compiles clean |
| 2 | ccpnmr2.5/c/ccpnmr/clouds/py_bacus.c | Bacus | bacus_type_methods | F | ✅ DONE compiles clean |
| 3 | ccpnmr2.5/c/ccpnmr/clouds/py_dist_constraint.c | DistConstraint | tp=py_handler_methods / mod=Dist_constraint_type_methods | T | ✅ DONE compiles clean |
| 4 | ccpnmr2.5/c/ccpnmr/clouds/py_dist_constraint_list.c | DistConstraintList | tp=py_handler_methods / mod=Dist_constraint_list_type_methods | T | ✅ DONE compiles clean |
| 5 | ccpnmr2.5/c/ccpnmr/clouds/py_dist_force.c | DistForce | tp=py_handler_methods / mod=Dist_force_type_methods | T | ✅ DONE compiles clean |
| 6 | ccpnmr2.5/c/ccpnmr/clouds/py_dynamics.c | Dynamics | tp=py_handler_methods / mod=Dynamics_type_methods | T | ✅ DONE compiles clean |
| 7 | ccpnmr2.5/c/ccpnmr/dynamics/py_atom_coord.c | DyAtomCoord | tp=py_handler_methods / mod=Atom_coord_type_methods (tp_name AtomCoord) | T | ✅ DONE compiles clean |
| 8 | ccpnmr2.5/c/ccpnmr/dynamics/py_atom_coord_list.c | DyAtomCoordList | tp=py_handler_methods / mod=Atom_coord_list_type_methods (tp_name AtomCoordList) | T | ✅ DONE compiles clean |
| 9 | ccpnmr2.5/c/ccpnmr/dynamics/py_dist_constraint.c | DyDistConstraint | tp=py_handler_methods / mod=Dist_constraint_type_methods (tp_name DistConstraint, tp_repr=0: dynamic multi-fprintf print dropped) | T | ✅ DONE compiles clean |
| 10 | ccpnmr2.5/c/ccpnmr/dynamics/py_dist_constraint_list.c | DyDistConstraintList | tp=py_handler_methods / mod=Dist_constraint_list_type_methods (tp_name DistConstraintList) | T | ✅ DONE compiles clean |
| 11 | ccpnmr2.5/c/ccpnmr/dynamics/py_dist_force.c | DyDistForce | tp=py_handler_methods / mod=Dist_force_type_methods (tp_name DistForce, +dist_power attr) | T | ✅ DONE compiles clean |
| 12 | ccpnmr2.5/c/ccpnmr/dynamics/py_dynamics.c | DyDynamics | tp=py_handler_methods / mod=Dynamics_type_methods (tp_name Dynamics) | T | ✅ DONE compiles clean |
| 13 | ccpnmr2.5/c/ccpnmr/analysis/py_contour_style.c | ContourStyle | Contour_style_type_methods | T | ✅ DONE — was 90% migrated w/ broken dangling getattr; restored `getattr_py_contour_style(...)` sig+brace (build-surfaced) |
| 14 | ccpnmr2.5/c/ccp/structure/py_atom.c | StructAtom | Atom_type_methods / tp=py_handler_methods | T | ✅ DONE this session (PyInit_StructAtom) |
| 15 | ccpnmr2.5/c/ccp/structure/py_bond.c | StructBond | Bond_type_methods / tp=py_handler_methods | T | ✅ DONE this session (PyInit_StructBond) |
| 16 | ccpnmr2.5/c/ccp/structure/py_struct_util.c | StructUtil | Structure_type_methods | F | ✅ DONE this session (PyInit_StructUtil, function-only) |
| 17 | ccpnmr2.5/c/other/cambridge/bayes/py_bayes.c | BayesPeakSeparator | BayesPeakSeparator_type_methods | F | ✅ DONE this session (PyInit_BayesPeakSeparator, function-only; vestigial type validated) |

### Already DONE (verified done — no action; re-check compiles at §5 build)
| File | init |
|---|---|
| ccpnmr2.5/c/ccpnmr/clouds/py_midge.c | PyInit_Midge |
| ccpnmr2.5/c/ccpnmr/clouds/py_atom_coord.c | PyInit_AtomCoord |
| ccpnmr2.5/c/ccpnmr/clouds/py_cloud_util.c | PyInit_CloudUtil |
| ccpnmr2.5/c/ccpnmr/analysis/py_contour_levels.c | PyInit_ContourLevels |
| ccpnmr2.5/c/ccpnmr/analysis/py_peak_list.c | PyInit_PeakList |
| ccpnmr2.5/c/ccpnmr/analysis/py_peak.c | (Peak type, registered inside PeakList) |
| ccpnmr2.5/c/memops/global/* (8 backbone exts) | Phase 1b done |

### OUT OF BUILD SCOPE (NOT in setup.py FAM → do NOT migrate unless user asks; flag in report)
- analysis/py_contour_file.c (ContourFile), py_peak_cluster.c (PeakCluster),
  py_slice_file.c (SliceFile), py_win_peak_list.c (WinPeakList) — imported by
  Analysis.py/ContourStore.py/WindowDraw.py/WindowFrame.py but have NO setup.py
  ext, so they'd need a setup.py change (forbidden here) to ever import.
- ccp/structure/py_structure.c (StructStructure) — not in FAM; unusual Py_FindMethod-in-body.
- other/meccano/pysrc/py_meccano.c; memops/global/py_{tk_util,draw_handler,gl_handler,tk_handler}.c — not in build.

## 3. PER-FILE MIGRATION (type-based "T")
Read the whole file, apply every applicable step, then self-check (§4).
1. `PY_MALLOC(<v>, struct Py_X, &X_type);` → `<v> = (Py_X) PyObject_New(struct Py_X, &X_type);`
2. destructor: `PY_FREE(self);` → `Py_TYPE(self)->tp_free(self);`
3. `print_py_X(FILE *fp...)` fprintf fn → `static PyObject *repr_py_X(PyObject *self)` returning
   `PyUnicode_FromString(buf)` (snprintf local `char buf[...];`); keep same format; pointer → `(void *)`+`%p`.
4. getattr: `(self,char *name)` → `(self,PyObject *attr_name)`. Body only forwards → whole body
   = `return PyObject_GenericGetAttr(self, attr_name);` (drop unused cast local). Else wrap
   strcmps in `if (PyUnicode_Check(attr_name)) { const char *name = PyUnicode_AsUTF8(attr_name); if (name) {...} }`,
   `equal_strings(name,"a")`→`strcmp(name,"a")==0`, last line `return PyObject_GenericGetAttr(self, attr_name);`
5. setattr: `(self,char *name,PyObject *v)` → `(self,PyObject *attr_name,PyObject *v)`. Body `return 0;` →
   just fix sig. Else add: `if (!PyUnicode_Check(attr_name)){ PyErr_SetString(ErrorObject,"attribute name must be a string"); return -1; }`
   `const char *name = PyUnicode_AsUTF8(attr_name); if (!name) return -1;`; `equal_strings`→`strcmp ... ==0`.
6. Replace the py2 `PyTypeObject X_type` literal with the EXACT 40-field layout below (§3-T).
   ⚠ TWO DIFFERENT PyMethodDef tables (file-1 trap): `tp_methods` = `py_handler_methods`
   (instance methods add/append/item, defined BEFORE the type); module-def = `<X>_type_methods`
   (the class constructor, defined AFTER the type). Pointing tp_methods at `<X>_type_methods`
   → "undeclared" ERROR. Keep leftover `PySequenceMethods X_sequence_methods` struct (matches
   done py_peak_list.c) but set tp_as_sequence=0. No `print_py_X`? set tp_repr=0. Dead `sq_`
   init emits harmless -Wincompatible-pointer-types warnings (done files have them) — OK.
7. Replace py2 `PY_MOD_INIT_FUNC initX(...)`+`Py_InitModule` with the §3-T module-init block.
For function-only "F" (bacus): only step 7 (use §3-F init; NO PyType_Ready; keep all PyMethodDef/bare fns as-is).
   Note check whether bacus declares `static PyObject *ErrorObject;` — if not, mirror `py_cloud_util.c`
   for the exact init/error handling it uses.

### §3-T CANONICAL PY3 TEMPLATE (copy VERBATIM, substitute X / struct Py_X / methods name)
```c
static PyTypeObject X_type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    "X",                                  /* tp_name */
    sizeof(struct Py_X),                  /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor) delete_py_X,             /* tp_dealloc */
    0,                                    /* tp_vectorcall */
    0,                                    /* tp_getattr */
    0,                                    /* tp_setattr */
    0,                                    /* tp_as_async */
    (reprfunc) repr_py_X,                 /* tp_repr */
    0,                                    /* tp_as_number */
    0,                                    /* tp_as_sequence */
    0,                                    /* tp_as_mapping */
    0,                                    /* tp_hash */
    0,                                    /* tp_call */
    0,                                    /* tp_str */
    (getattrofunc) getattr_py_X,          /* tp_getattro */
    (setattrofunc) setattr_py_X,          /* tp_setattro */
    0,                                    /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                   /* tp_flags */
    "X -- <short tp_doc>",                /* tp_doc */
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                    /* tp_weaklistoffset */
    0,                                    /* tp_iter */
    0,                                    /* tp_iternext */
    py_handler_methods,                   /* tp_methods = FILE's INSTANCE table (usually `py_handler_methods`) */
    0,                                    /* tp_members */
    0,                                    /* tp_getset */
    0,                                    /* tp_base */
    0,                                    /* tp_dict */
    0,                                    /* tp_descr_get */
    0,                                    /* tp_descr_set */
    0,                                    /* tp_dictoffset */
    0,                                    /* tp_init */
    0,                                    /* tp_alloc */
    0,                                    /* tp_new */
};

static struct PyModuleDef X_module_def =
{
    PyModuleDef_HEAD_INIT,
    "X",
    "CCPNMR X module (Python 3 compatible)",
    -1,
    X_type_methods   /* FILE's MODULE-constructor table (usually `<X>_type_methods`) */
};

PyMODINIT_FUNC PyInit_X(void)
{
    if (PyType_Ready(&X_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&X_module_def);
    if (!m)
        return NULL;

    ErrorObject = PyErr_NewException("X.error", NULL, NULL);
    if (!ErrorObject)
    {
        Py_DECREF(m);
        return NULL;
    }
    Py_INCREF(ErrorObject);
    if (PyDict_SetItemString(PyModule_GetDict(m), "error", ErrorObject) < 0)
    {
        Py_DECREF(ErrorObject);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
```
### §3-F CANONICAL PY3 TEMPLATE (function-only, e.g. Bacus)
```c
static struct PyModuleDef bacus_module_def =
{
    PyModuleDef_HEAD_INIT,
    "Bacus",
    "CCPNMR Bacus module (Python 3 compatible)",
    -1,
    bacus_type_methods
};

PyMODINIT_FUNC PyInit_Bacus(void)
{
    PyObject *m = PyModule_Create(&bacus_module_def);
    if (!m)
        return NULL;

    ErrorObject = PyErr_NewException("Bacus.error", NULL, NULL);
    if (!ErrorObject) { Py_DECREF(m); return NULL; }
    Py_INCREF(ErrorObject);
    if (PyDict_SetItemString(PyModule_GetDict(m), "error", ErrorObject) < 0)
    { Py_DECREF(ErrorObject); Py_DECREF(m); return NULL; }

    return m;
}
```

## 4. SELF-CHECK per file (from repo root) — fixes COMPILE errors; warnings OK
```
PYI=$( .venv/bin/python -c "import sysconfig;print(sysconfig.get_paths()['include'])" )
PYP=$( .venv/bin/python -c "import sysconfig;print(sysconfig.get_path if 0 else sysconfig.get_paths()['platinclude'])" )   # (platinclude)
cd /home/logan/software/ccpnmr2.5.2-qwen
# INCLUDE = one of:  ccpnmr2.5/c/ccpnmr/clouds | ccpnmr2.5/c/ccpnmr/dynamics |
#                   ccpnmr2.5/c/ccpnmr/analysis | ccpnmr2.5/c/ccp/structure |
#                   ccpnmr2.5/c/other/cambridge/bayes
gcc -fsyntax-only -I<INCLUDE> -Iccpnmr2.5/c/memops/global -I"$PYI" -I"$PYP" <FAMILYFILE>
```
Correct PYI/PYP one-liner (verified): `PYP=$( .venv/bin/python -c "import sysconfig;print(sysconfig.get_paths()['platinclude'])" )`

## 5. AUTHORITATIVE BUILD (only after all 17 DONE) — parent only
```
cd /home/logan/software/ccpnmr2.5.2-qwen
.venv/bin/python setup.py build_ext --inplace 2>&1 | tee /tmp/p2a_build.log
```
Builds flat `<Name>.so`. Then copy each built `.so` onto its symlink target, e.g.:
`clouds/AtomCoord.so`, `atom_coord_list.so`... — verify symlink layout:
`ls -l ccpnmr2.5/python/ccpnmr/c/ ccp/c/ cambridge/c/ 2>/dev/null` and
`readlink -f` the `.so` names to confirm where the build output lands vs. where the
package symlink expects it. Copy/overwrite the stale py2 `.so` at each target.

## 6. IMPORT TEST
```
cd /home/logan/software/ccpnmr2.5.2-qwen
.venv/bin/python -c "import ccpnmr.c.AtomCoordList, ccpnmr.c.Midge, ccpnmr.c.DyDynamics, ccpnmr.c.StructAtom, ccpnmr.c.BayesPeakSeparator; print('OK: ext imports')"
# Then re-run the tree-wide smoke and diff the OK count upward (was 1556/1726 before phase2a):
MPLBACKEND=Agg .venv/bin/python import_smoke.py 2>&1 | tail -40
```

## 7. FINALIZE
- Update memory `project/ccpnmr-modernization.md` with phase2a outcomes (which exts, import OK delta).
- Report: DONE/FAILED per file, build result, import-smoke delta, and the OUT-OF-SCOPE list (§2)
  for the user to decide (setup.py change needed to build them).

## LOG (append one line per completed file, newest last)
- 2026-08-21 🏁 **PHASE 2a COMPLETE**: 4 files migrated this session (StructAtom/StructBond/StructUtil/BayesPeakSeparator). Build surfaced + I fixed 2 pre-existing bugs: `py_contour_style.c` dangling `getattr` (deleted sig) and `py_peak_list.c` `PyInt_Check/PyInt_AsLong`→`PyLong_*` (was the "undefined symbol: PyInt_AsLong" import error). `setup.py build_ext --inplace` → all 22 FAM + 8 backbone built; copied 30 `.so` onto package `c/` symlink targets (7 out-of-scope skipped). 22/22 import OK, functional smoke (construct/setters/dealloc + BayesPeakSeparator.test_bayes + ShapeFile) PASS, import_smoke **1556→1561 OK (+5)**. `.gitignore` keeps build/+*.so ignored, `!ccpnmr2.5/c/**/*.so` keeps the tracked family-dir .so.
- 2026-08-21 ⚛ **dynamics 6/6 DONE**: #7 py_atom_coord(PyInit_DyAtomCoord) #8 py_atom_coord_list(PyInit_DyAtomCoordList) #9 py_dist_constraint(PyInit_DyDistConstraint, repr dropped) #10 py_dist_constraint_list(PyInit_DyDistConstraintList) #11 py_dist_force(PyInit_DyDistForce) #12 py_dynamics(PyInit_DyDynamics, PyLong_AsLong). All gcc clean. tp_name stays unprefixed (AtomCoord/Dynamics/...); only the MODULE name is Dy* -prefixed.
- 2026-08-21 ☁ **clouds 6/6 DONE**: #1 atom_coord_list(PyInit_AtomCoordList) #2 bacus(PyInit_Bacus,F) #3 dist_constraint(PyInit_DistConstraint) #4 dist_constraint_list(PyInit_DistConstraintList) #5 dist_force(PyInit_DistForce) #6 dynamics(PyInit_Dynamics). All gcc -fsyntax-only clean. Traps hit: tp_methods=py_handler_methods(not ctor table); PyInt_AsLong→PyLong_AsLong; real print→repr(snprintf+PyUnicode_FromString, wire tp_repr); commented print→tp_repr=0.
- (init) checkpoint created 2026-08-21; 17 files PENDING.

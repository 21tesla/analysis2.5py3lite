# Phase 2a — C extension py2→py3 migration recipe (CANONICAL)

You are migrating CCPNMR C extension Python-2 → Python-3 *binding* files. The recipe
below is VERIFIED — it produced working `py_midge.c` and `py_atom_coord.c`. Do NOT
guess; mirror the references.

## References (already migrated & proven — READ them first)
- Type-based:   `ccpnmr2.5/c/ccpnmr/clouds/py_midge.c`  (simple type, one instance method)
- Type-based:   `ccpnmr2.5/c/ccpnmr/clouds/py_atom_coord.c` (type with real get/set attrs)
- Function-only:`ccpnmr2.5/c/ccpnmr/clouds/py_cloud_util.c` (PyMethodDef table, no type)

Repo root: `/home/logan/software/ccpnmr2.5.2-qwen`

## THE RECIPE — type-based file (each step is mechanical)

1. **new_py_X** — replace
   `PY_MALLOC(<v>, struct Py_X, &X_type);`
   → `<v> = (Py_X) PyObject_New(struct Py_X, &X_type);`

2. **delete_py_X (destructor)** — replace `PY_FREE(self);` → `Py_TYPE(self)->tp_free(self);`

3. **print → repr** — replace
   `static int print_py_X(PyObject *self, FILE *fp, int flags){ ... fprintf(fp, "<...>", args); return 0; }`
   →
   `static PyObject *repr_py_X(PyObject *self){ ... char buf[64]; snprintf(buf, sizeof buf, "<...>", args); return PyUnicode_FromString(buf); }`
   Keep the same format string. If the format printed a POINTER via an `(int)` cast, use `(void *)` + `%p`.

4. **getattr** — signature `(PyObject *self, char *name)` → `(PyObject *self, PyObject *attr_name)`.
   - If the body ONLY forwards (`return Py_FindMethod(py_handler_methods, self, name);`),
     replace the whole body with: `return PyObject_GenericGetAttr(self, attr_name);`
     and REMOVE any unused `Py_X py_... = (Py_X) self;` local.
   - Otherwise wrap comparisons in
     `if (PyUnicode_Check(attr_name)) { const char *name = PyUnicode_AsUTF8(attr_name); if (name) { ... } }`
     with the `...` being your comparisons using `strcmp(name, "<a>") == 0` (was `equal_strings(name, "<a>")`),
     and the final line `return PyObject_GenericGetAttr(self, attr_name);` (was `return Py_FindMethod(py_handler_methods, self, name);`).

5. **setattr** — signature `(PyObject *self, char *name, PyObject *value)` → `(PyObject *self, PyObject *attr_name, PyObject *value)`.
   - If the body is `return 0;`, leave the body (just fix the signature).
   - Otherwise add right after the local declarations:
     `if (!PyUnicode_Check(attr_name)) { PyErr_SetString(ErrorObject, "attribute name must be a string"); return -1; }`
     `const char *name = PyUnicode_AsUTF8(attr_name); if (!name) return -1;`
     and replace `equal_strings(name, "<a>")` → `strcmp(name, "<a>") == 0`.

6. **PyTypeObject literal (THE BIG ONE)** — replace the py2 block:
   ```c
   static PyTypeObject X_type =
   {
   #ifdef WIN64
       1, NULL,
   #else
       PyObject_HEAD_INIT(&PyType_Type)
   #endif
       0,
       "X", /* name */
       sizeof(struct Py_X), /* basicsize */
       0, /* itemsize */
       delete_py_X, /* destructor */
       print_py_X, /* printfunc */
       getattr_py_X, /* getattr */
       setattr_py_X, /* setattr */
       0, /* cmpfunc */
       0, /* reprfunc */
       0, /* PyNumberMethods */
       /* ... */ /* PySequenceMethods */
   };
   ```
   with the EXACT py3 40-field layout used in `py_midge.c` (copy it VERBATIM), substituting
   `X_type`, `sizeof(struct Py_X)`, `delete_py_X`, `(reprfunc) repr_py_X`,
   `(getattrofunc) getattr_py_X`, `(setattrofunc) setattr_py_X`, `py_handler_methods` at `tp_methods`,
   and a short `tp_doc`. Field ORDER/POSITIONS ARE CRITICAL — copy the reference, do not invent fields.

7. **module init** — replace the py2
   ```c
   PY_MOD_INIT_FUNC initX(void){ ... m = Py_InitModule("X", <methods>); ErrorObject = PyErr_NewException("X.error", NULL, NULL); Py_INCREF(ErrorObject); PyModule_AddObject(m, "error", ErrorObject); if (PyErr_Occurred()) Py_FatalError(...); }
   ```
   with:
   ```c
   static struct PyModuleDef X_module_def =
   {
       PyModuleDef_HEAD_INIT,
       "X",
       "CCPNMR X module (Python 3 compatible)",
       -1,
       <methods>
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
   Use the CORRECT module name `"X"` — read it from the file's py2 `Py_InitModule("X", ...)` —
   and the correct methods-table identifier (the `PyMethodDef` array name in the file).

## FUNCTION-ONLY file (NO PyTypeObject, e.g. py_bacus.c)
- Keep every `PyMethodDef` table and every `static PyObject *xxx(PyObject *self, PyObject *args)`
  function EXACTLY as-is (they already call py3-compatible `PyArg_ParseTuple`/`PyBuildValue`).
- Apply ONLY step 7 → `PyInit_<Name>` (name from the file's py2 `Py_InitModule("...", ...)`).
- Do NOT add `PyType_Ready` (no type). Follow `py_cloud_util.c` exactly.

## CONSTRAINTS
- Do NOT modify `setup.py`, do NOT run `setup.py`, do NOT build, do NOT copy any `.so`.
  The parent does the authoritative build + import test once all files are edited.
- Keep all license headers, comments, and unrelated C logic intact. Change only the binding parts above.
- The shared helpers `python_util.c`, `utility.c`, `python_util.h`, `utility.h` are ALREADY py3 — DO NOT edit them.

## SELF-CHECK after each file (from repo root)
```
PYI=$( /home/logan/software/ccpnmr2.5.2-qwen/.venv/bin/python -c "import sysconfig;print(sysconfig.get_paths()['include'])" )
PYP=$( /home/logan/software/ccpnmr2.5.2-qwen/.venv/bin/python -c "import sysconfig;print(sysconfig.get_paths()['platinclude'])" )
cd /home/logan/software/ccpnmr2.5.2-qwen
gcc -fsyntax-only -Iccpnmr2.5/c/ccpnmr/<FAMILY> -Iccpnmr2.5/c/memops/global -I"$PYI" -I"$PYP" \
    ccpnmr2.5/c/ccpnmr/<FAMILY>/<FILENAME>
```
`<FAMILY>` is `clouds`, `dynamics`, or `analysis`. Fix COMPILE errors (they mean a bad edit).
Warnings are acceptable.

## SAFETY NET
If you hit ANY other py2-only C API token not covered above (e.g. `Py_FindMethod`, `equal_strings`,
`PyObject_HEAD_INIT`, `Py_InitModule`, `PY_MOD_INIT`, a 13-field `PyTypeObject`, py2 sequence/number
methods `x_item`/`x_slice`, `PyArg_VaryArguments`, `PyList2DoubleSeq`, or any `PY_XXX` macro),
STOP, do not guess — record it in your report.

## REPORT
One line per file: `DONE <file>  (PyInit_<Name>)`  or  `FAILED <file>: <reason>`.
Then a one-paragraph summary of anything unusual.

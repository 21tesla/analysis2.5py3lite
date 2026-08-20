/*
======================COPYRIGHT/LICENSE START==========================

py_mem_cache.c: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================
*/
#include "py_mem_cache.h"

#include "python_util.h"

/* Locally-raised exception type */
static PyObject *ErrorObject;

static PyTypeObject Mem_cache_type;

Bool is_py_mem_cache(PyObject *obj)
{
    return valid_py_object(obj, &Mem_cache_type);
}

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *resize(PyObject *self, PyObject *args)
{
    Py_Mem_cache obj = (Py_Mem_cache) self;
    Mem_cache mem_cache = obj->mem_cache;
    int max_size;

    if (!PyArg_ParseTuple(args, "i", &max_size))
        RETURN_OBJ_ERROR("need one argument: max_size");

    if (resize_mem_cache(mem_cache, max_size) == CCPN_ERROR)
        RETURN_OBJ_ERROR("resizing Mem_cache");

    Py_RETURN_NONE;
}

static PyObject *clear(PyObject *self, PyObject *args)
{
    Py_Mem_cache obj = (Py_Mem_cache) self;
    Mem_cache mem_cache = obj->mem_cache;

    clear_mem_cache(mem_cache);

    Py_RETURN_NONE;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "resize",		resize,		METH_VARARGS },
    { "clear",		clear,		METH_VARARGS },
    { NULL,		NULL,		0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_mem_cache(int max_size)
{
    Py_Mem_cache obj;
    Mem_cache mem_cache;

    mem_cache = new_mem_cache(max_size, NULL, NULL);

    if (!mem_cache)
        RETURN_OBJ_ERROR("allocating Mem_cache object");

    obj = (Py_Mem_cache) PyObject_New(struct Py_Mem_cache, &Mem_cache_type);

    if (!obj)
    {
        delete_mem_cache(mem_cache);
        RETURN_OBJ_ERROR("allocating Py_Mem_cache object");
    }

    obj->mem_cache = mem_cache;
    return (PyObject *) obj;
}

static void delete_py_mem_cache(PyObject *self)
{
    Py_Mem_cache obj = (Py_Mem_cache) self;
    Mem_cache mem_cache = obj->mem_cache;

    delete_mem_cache(mem_cache);
    Py_TYPE(self)->tp_free(self);
}

/* Python 3 tp_getattro — receives a unicode object */
static PyObject *getattr_py_mem_cache(PyObject *self, PyObject *attr_name)
{
    /* Fall back to PyObject_GenericGetAttr which handles tp_methods */
    return PyObject_GenericGetAttr(self, attr_name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

static PyObject *mem_cache_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = { "max_size", NULL };
    int max_size;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i", kwlist, &max_size))
        RETURN_OBJ_ERROR("must have one argument: max_size");

    return new_py_mem_cache(max_size);
}

static PyTypeObject Mem_cache_type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    "MemCache",                              /* tp_name */
    sizeof(struct Py_Mem_cache),             /* tp_basicsize */
    0,                                       /* tp_itemsize */
    (destructor) delete_py_mem_cache,        /* tp_dealloc */
    0,                                       /* tp_vectorcall */
    0,                                       /* tp_getattr */
    0,                                       /* tp_setattr */
    0,                                       /* tp_as_async */
    0,                                       /* tp_repr */
    0,                                       /* tp_as_number */
    0,                                       /* tp_as_sequence */
    0,                                       /* tp_as_mapping */
    0,                                       /* tp_hash */
    0,                                       /* tp_call */
    0,                                       /* tp_str */
    (getattrofunc) getattr_py_mem_cache,     /* tp_getattro */
    0,                                       /* tp_setattro */
    0,                                       /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                      /* tp_flags */
    "MemCache — MOPS shared-memory model",   /* tp_doc */
    0,                                       /* tp_traverse */
    0,                                       /* tp_clear */
    0,                                       /* tp_richcompare */
    0,                                       /* tp_weaklistoffset */
    0,                                       /* tp_iter */
    0,                                       /* tp_iternext */
    py_handler_methods,                      /* tp_methods */
    0,                                       /* tp_members */
    0,                                       /* tp_getset */
    0,                                       /* tp_base */
    0,                                       /* tp_dict */
    0,                                       /* tp_descr_get */
    0,                                       /* tp_descr_set */
    0,                                       /* tp_dictoffset */
    0,                                       /* tp_init */
    0,                                       /* tp_alloc */
    mem_cache_new,                           /* tp_new */
};

/*****************************************************************************
 * MODULE LOGIC (Python 3 API)
 *****************************************************************************/

static struct PyMethodDef Mem_cache_type_methods[] =
{
    { "MemCache",	(PyCFunction) mem_cache_new,	METH_VARARGS | METH_KEYWORDS },
    { NULL,		NULL,			0 }
};

static struct PyModuleDef mem_cache_module_def =
{
    PyModuleDef_HEAD_INIT,
    "MemCache",
    "CCPNMR MOPS shared-memory model (Python 3 compatible)",
    -1,
    Mem_cache_type_methods
};

PyMODINIT_FUNC PyInit_MemCache(void)
{
    if (PyType_Ready(&Mem_cache_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&mem_cache_module_def);
    if (!m)
        return NULL;

    PyObject *d = PyModule_GetDict(m);
    if (PyDict_SetItemString(d, "MemCache", (PyObject *) &Mem_cache_type) < 0)
    {
        Py_DECREF(m);
        return NULL;
    }

    /* Create exception object and add to module */
    ErrorObject = PyErr_NewException("MemCache.error", NULL, NULL);
    if (ErrorObject != NULL)
    {
        Py_INCREF(ErrorObject);
        if (PyDict_SetItemString(d, "error", ErrorObject) < 0)
        {
            Py_DECREF(m);
            return NULL;
        }
    }

    return m;
}

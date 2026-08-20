
/*
======================COPYRIGHT/LICENSE START==========================

py_store_file.c: Part of the CcpNmr Analysis program

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
#include "py_store_file.h"

#include "python_util.h"
#include "utility.h"

#define  NCOLORS  3

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Store_file_type;

Bool is_py_store_file(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Store_file_type
    return (obj->ob_type == &Store_file_type);
*/
    return valid_py_object(obj, &Store_file_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static struct PyMethodDef py_file_methods[] =
{
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_store_file(CcpnString file_name, int ndim,
			int xdim, int ydim, PyObject *block_size_obj)
{
    int n;
    Py_Store_file obj;
    Store_file store_file;
    int block_size[MAX_NDIM];
    Line error_msg;

    if ((get_python_int_array(block_size_obj, MAX_NDIM, &n, block_size,
                                                error_msg) == CCPN_ERROR)
        || (n != ndim))
    {
        sprintf(error_msg, "block_size must be int list or tuple of size %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    store_file = new_store_file(file_name, ndim, xdim, ydim, block_size, error_msg);

    if (!store_file)
	RETURN_OBJ_ERROR(error_msg);

    obj = (Py_Store_file) PyObject_New(struct Py_Store_file, &Store_file_type);

    if (!obj)
    {
	delete_store_file(store_file);

	RETURN_OBJ_ERROR("allocating Py_Store_file object");
    }

    obj->store_file = store_file;

    return (PyObject *) obj;
}

static void delete_py_store_file(PyObject *self)
{
    Py_Store_file obj = (Py_Store_file) self;
    Store_file store_file = obj->store_file;

    delete_store_file(store_file);

    Py_TYPE(self)->tp_free(self);
}

/*
static int print_py_store_file(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_file\n");

    return 0;
}
*/

static PyObject *getattr_py_store_file(PyObject *self, PyObject *attr_name)
{
    Py_Store_file obj = (Py_Store_file) self;
    Store_file store_file = obj->store_file;

    const char *name = PyUnicode_AsUTF8(attr_name);
    if (!name)
        return NULL;

    if (strcmp(name, "have_pos") == 0)
        return PyLong_FromLong((long) store_file->have_pos);
    else if (strcmp(name, "have_neg") == 0)
        return PyLong_FromLong((long) store_file->have_neg);
    else if (strcmp(name, "dir_size") == 0)
        return PyLong_FromLong((long) store_file->dir_size);

    return PyObject_GenericGetAttr(self, attr_name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Store_file_sequence_methods =
{
    Store_file_length,
    Store_file_concat,
    Store_file_repeat,
    Store_file_item,
    Store_file_slice,
    Store_file_ass_item,
    Store_file_ass_slice
};

static PySequenceMethods Store_file_sequence_methods =
{
    Store_file_length,
    0,
    0,
    Store_file_item,
    0,
    Store_file_ass_item,
    0
};
*/

static PyTypeObject Store_file_type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    "StoreFile",                              /* tp_name */
    sizeof(struct Py_Store_file),             /* tp_basicsize */
    0,                                        /* tp_itemsize */
    (destructor) delete_py_store_file,        /* tp_dealloc */
    0,                                        /* tp_vectorcall */
    0,                                        /* tp_getattr */
    0,                                        /* tp_setattr */
    0,                                        /* tp_as_async */
    0,                                        /* tp_repr */
    0,                                        /* tp_as_number */
    0,                                        /* tp_as_sequence */
    0,                                        /* tp_as_mapping */
    0,                                        /* tp_hash */
    0,                                        /* tp_call */
    0,                                        /* tp_str */
    (getattrofunc) getattr_py_store_file,     /* tp_getattro */
    0,                                        /* tp_setattro */
    0,                                        /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                       /* tp_flags */
    "StoreFile -- NMR store file model",      /* tp_doc */
    0,                                        /* tp_traverse */
    0,                                        /* tp_clear */
    0,                                        /* tp_richcompare */
    0,                                        /* tp_weaklistoffset */
    0,                                        /* tp_iter */
    0,                                        /* tp_iternext */
    py_file_methods,                          /* tp_methods */
    0,                                        /* tp_members */
    0,                                        /* tp_getset */
    0,                                        /* tp_base */
    0,                                        /* tp_dict */
    0,                                        /* tp_descr_get */
    0,                                        /* tp_descr_set */
    0,                                        /* tp_dictoffset */
    0,                                        /* tp_init */
    0,                                        /* tp_alloc */
    0,                                        /* tp_new */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Store_file(PyObject *self, PyObject *args)
{
    int ndim, xdim, ydim;
    CcpnString file_name;
    PyObject *block_size_obj;

    if (!PyArg_ParseTuple(args, "siiiO", &file_name, &ndim, &xdim, &ydim, &block_size_obj))
        RETURN_OBJ_ERROR("must have arguments: file_name, ndim, xdim, ydim, block_size");

    return new_py_store_file(file_name, ndim, xdim, ydim, block_size_obj);
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Store_file_type_methods[] =
{
    { "StoreFile",	(PyCFunction) init_Py_Store_file,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import StoreFile" in 
* a Python program. The function is usually called "initStore_file": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

static struct PyModuleDef store_file_module_def =
{
    PyModuleDef_HEAD_INIT,
    "StoreFile",
    "CCPNMR Store File module (Python 3 compatible)",
    -1,
    Store_file_type_methods
};

PyMODINIT_FUNC PyInit_StoreFile(void)
{
    if (PyType_Ready(&Store_file_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&store_file_module_def);
    if (!m)
        return NULL;

    ErrorObject = PyErr_NewException("StoreFile.error", NULL, NULL);
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

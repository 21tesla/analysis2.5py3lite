
/*
======================COPYRIGHT/LICENSE START==========================

py_store_handler.c: Part of the CcpNmr Analysis program

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
#include "py_store_handler.h"

#include "python_util.h"

#define  NCOLORS  3

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Store_handler_type;

Bool is_py_store_handler(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Store_handler_type
    return (obj->ob_type == &Store_handler_type);
*/
    return valid_py_object(obj, &Store_handler_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static struct PyMethodDef py_handler_methods[] =
{
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_store_handler(CcpnString file_name, Bool swap)
{
    Py_Store_handler obj;
    Store_handler store_handler;
    Line error_msg;

    store_handler = new_store_handler(file_name, swap, error_msg);

    if (!store_handler)
	RETURN_OBJ_ERROR(error_msg);

    PY_MALLOC(obj, struct Py_Store_handler, &Store_handler_type);

    if (!obj)
    {
	delete_store_handler(store_handler);

	RETURN_OBJ_ERROR("allocating Py_Store_handler object");
    }

    obj->store_handler = store_handler;

    return (PyObject *) obj;
}

static void delete_py_store_handler(PyObject *self)
{
    Py_Store_handler obj = (Py_Store_handler) self;
    Store_handler store_handler = obj->store_handler;

    delete_store_handler(store_handler);

    PY_FREE(self);
}

/*
static int print_py_store_handler(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_store_handler(PyObject *self, char *name)
{
/*
    Store_handler *obj = (Store_handler *) self;
    Random_access *a = obj->py_handler;

    if (equal_strings(name, "par_file"))
	return Py_BuildValue("s", a->par_file);
    else if (equal_strings(name, "access_method"))
	return Py_BuildValue("s", access_method_name(a->access_method));
    else if (equal_strings(name, "data_format"))
	return get_Store_handler_format(a);
    else
*/
	return Py_FindMethod(py_handler_methods, self, name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Store_handler_sequence_methods =
{
    Store_handler_length,
    Store_handler_concat,
    Store_handler_repeat,
    Store_handler_item,
    Store_handler_slice,
    Store_handler_ass_item,
    Store_handler_ass_slice
};

static PySequenceMethods Store_handler_sequence_methods =
{
    Store_handler_length,
    0,
    0,
    Store_handler_item,
    0,
    Store_handler_ass_item,
    0
};
*/

static PyTypeObject Store_handler_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "StoreHandler", /* name */
    sizeof(struct Py_Store_handler), /* basicsize */
    0, /* itemsize */
    delete_py_store_handler, /* destructor */
    0, /* printfunc */
    getattr_py_store_handler, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Store_handler_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Store_handler(PyObject *self, PyObject *args)
{
    int swap = 0;
    CcpnString file_name;

    if (!PyArg_ParseTuple(args, "s|i", &file_name, &swap))
        RETURN_OBJ_ERROR("must have arguments: file_name [ swap ]");

    return new_py_store_handler(file_name, (Bool) swap);
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Store_handler_type_methods[] =
{
    { "StoreHandler",	(PyCFunction) init_Py_Store_handler,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import StoreHandler" in 
* a Python program. The function is usually called "initStore_handler": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initStoreHandler(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Store_handler_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("StoreHandler", Store_handler_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("StoreHandler.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module StoreHandler");
}

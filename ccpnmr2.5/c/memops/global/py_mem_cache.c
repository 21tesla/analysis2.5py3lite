
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

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Mem_cache_type;

Bool is_py_mem_cache(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Mem_cache_type
    return (obj->ob_type == &Mem_cache_type);
*/
    return valid_py_object(obj, &Mem_cache_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

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

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *clear(PyObject *self, PyObject *args)
{
    Py_Mem_cache obj = (Py_Mem_cache) self;
    Mem_cache mem_cache = obj->mem_cache;

    clear_mem_cache(mem_cache);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "resize",		resize,			METH_VARARGS },
    { "clear",		clear,			METH_VARARGS },
    { NULL,		NULL,			0 }
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

    PY_MALLOC(obj, struct Py_Mem_cache, &Mem_cache_type);

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

/*
    printf("in delete_py_mem_cache\n");
*/

    delete_mem_cache(mem_cache);

    PY_FREE(self);
}

/*
static int print_py_mem_cache(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_mem_cache(PyObject *self, char *name)
{
/*
    Mem_cache *obj = (Mem_cache *) self;
    Random_access *a = obj->py_handler;

    if (equal_strings(name, "par_file"))
	return Py_BuildValue("s", a->par_file);
    else if (equal_strings(name, "access_method"))
	return Py_BuildValue("s", access_method_name(a->access_method));
    else if (equal_strings(name, "data_format"))
	return get_Mem_cache_format(a);
    else
*/
	return Py_FindMethod(py_handler_methods, self, name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Mem_cache_sequence_methods =
{
    Mem_cache_length,
    Mem_cache_concat,
    Mem_cache_repeat,
    Mem_cache_item,
    Mem_cache_slice,
    Mem_cache_ass_item,
    Mem_cache_ass_slice
};

static PySequenceMethods Mem_cache_sequence_methods =
{
    Mem_cache_length,
    0,
    0,
    Mem_cache_item,
    0,
    Mem_cache_ass_item,
    0
};
*/

static PyTypeObject Mem_cache_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "MemCache", /* name */
    sizeof(struct Py_Mem_cache), /* basicsize */
    0, /* itemsize */
    delete_py_mem_cache, /* destructor */
    0, /* printfunc */
    getattr_py_mem_cache, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Mem_cache_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Mem_cache(PyObject *self, PyObject *args)
{
    int max_size;

    if (!PyArg_ParseTuple(args, "i", &max_size))
        RETURN_OBJ_ERROR("must have one argument: max_size");

    return new_py_mem_cache(max_size);
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Mem_cache_type_methods[] =
{
    { "MemCache",	(PyCFunction) init_Py_Mem_cache,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import MemCache" in 
* a Python program. The function is usually called "initMem_cache": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initMemCache(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Mem_cache_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("MemCache", Mem_cache_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("MemCache.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module MemCache");
}

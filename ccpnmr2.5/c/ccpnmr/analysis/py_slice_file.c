/*
======================COPYRIGHT/LICENSE START==========================

py_slice_file.c: Part of the CcpNmr Analysis program

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
#include "py_slice_file.h"

#include "py_block_file.h"
#include "py_draw_handler.h"
#include "py_mem_cache.h"
#include "python_util.h"
#include "utility.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Slice_file_type;

Bool is_py_slice_file(PyObject *obj)
{
    /*  below does not work because different *.so files end up
        with different addresses for Slice_file_type
        return (obj->ob_type == &Slice_file_type);
    */
    return valid_py_object(obj, &Slice_file_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *draw(PyObject *self, PyObject *args)
{
    int n, first, last;
    float position[MAX_NDIM];
    PyObject *position_obj, *handler_obj;
    Py_Slice_file obj = (Py_Slice_file) self;
    Slice_file slice_file = obj->slice_file;
    int ndim = slice_file->block_file->ndim;
    Py_draw_handler py_draw_handler;
    Drawing_funcs *drawing_funcs;
    Generic_ptr handler;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "OiiO", &handler_obj, &first, &last, &position_obj))
        RETURN_OBJ_ERROR("need four arguments: handler, first (int), last (int), position (float tuple)");

    py_draw_handler = new_py_draw_handler(handler_obj);
    if (!py_draw_handler)
        RETURN_OBJ_ERROR("first argument must be handler object");

    handler = py_draw_handler->handler;
    drawing_funcs = py_draw_handler->drawing_funcs;
    delete_py_draw_handler(py_draw_handler);

    if ((get_python_float_array(position_obj, MAX_NDIM, &n, position, error_msg) == CCPN_ERROR)
            || (n != ndim))
    {
        sprintf(error_msg,
                "fourth argument, position, must be list or tuple of size %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (draw_slice_file(slice_file, first, last, position,
                        drawing_funcs, handler, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawAll(PyObject *self, PyObject *args)
{
    int n, first[MAX_NDIM], last[MAX_NDIM], ncomponents = 0, *components = NULL;
    PyObject *first_obj, *last_obj, *handler_obj, *components_obj = NULL;
    Py_Slice_file obj = (Py_Slice_file) self;
    Slice_file slice_file = obj->slice_file;
    int ndim = slice_file->block_file->ndim;
    Py_draw_handler py_draw_handler;
    Drawing_funcs *drawing_funcs;
    Generic_ptr handler;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "OOO|O", &handler_obj, &first_obj, &last_obj, &components_obj))
        RETURN_OBJ_ERROR("need arguments: handler, first, last [, components]");

    py_draw_handler = new_py_draw_handler(handler_obj);
    if (!py_draw_handler)
        RETURN_OBJ_ERROR("first argument must be handler object");

    handler = py_draw_handler->handler;
    drawing_funcs = py_draw_handler->drawing_funcs;
    delete_py_draw_handler(py_draw_handler);

    if ((get_python_int_array(first_obj, MAX_NDIM, &n, first, error_msg) == CCPN_ERROR)
            || (n != ndim))
    {
        sprintf(error_msg,
                "second argument, first, must be list or tuple of size %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_int_array(last_obj, MAX_NDIM, &n, last, error_msg) == CCPN_ERROR)
            || (n != ndim))
    {
        sprintf(error_msg,
                "second argument, last, must be list or tuple of size %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (components_obj && (components_obj != Py_None))
    {
        if (get_python_int_alloc_array(components_obj, &ncomponents, &components, error_msg) == CCPN_ERROR)
            RETURN_OBJ_ERROR(error_msg);
    }

    if (draw_all_slice_file(slice_file, first, last, ncomponents, components,
                            drawing_funcs, handler, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "draw",		draw,			METH_VARARGS },
    { "drawAll",	drawAll,		METH_VARARGS },
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_slice_file(int orient, int dim,
                                   PyObject *block_file_obj, PyObject *mem_cache_obj)
{
    Slice_file slice_file;
    Py_Slice_file obj;
    Block_file block_file;

    if (!is_py_some_block_file(block_file_obj))
        RETURN_OBJ_ERROR("must pass block_file object");

    if (!is_py_mem_cache(mem_cache_obj))
        RETURN_OBJ_ERROR("must pass mem_cache object");

    if (is_py_block_file(block_file_obj))
        block_file = ((Py_Block_file) block_file_obj)->block_file;
    else
        block_file = (Block_file) (((Py_Shape_block_file) block_file_obj)->shape_block_file);

    slice_file = new_slice_file(orient, dim, block_file,
                                ((Py_Mem_cache) mem_cache_obj)->mem_cache);

    if (!slice_file)
        RETURN_OBJ_ERROR("allocating Slice_file object");

    PY_MALLOC(obj, struct Py_Slice_file, &Slice_file_type);

    if (!obj)
    {
        delete_slice_file(slice_file);

        RETURN_OBJ_ERROR("allocating Py_Slice_file object");
    }

    obj->slice_file = slice_file;

    return (PyObject *) obj;
}

static void delete_py_slice_file(PyObject *self)
{
    Py_Slice_file obj = (Py_Slice_file) self;
    Slice_file slice_file = obj->slice_file;

    /*
        printf("in delete_py_slice_file\n");
    */

    delete_slice_file(slice_file);

    PY_FREE(self);
}

/*
static int print_py_slice_file(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_slice_file(PyObject *self, char *name)
{
    Py_Slice_file obj = (Py_Slice_file) self;
    Slice_file slice_file = obj->slice_file;

    if (equal_strings(name, "dim"))
        return Py_BuildValue("i", slice_file->dim);
    else
        return Py_FindMethod(py_handler_methods, self, name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Slice_file_sequence_methods =
{
    Slice_file_length,
    Slice_file_concat,
    Slice_file_repeat,
    Slice_file_item,
    Slice_file_slice,
    Slice_file_ass_item,
    Slice_file_ass_slice
};

static PySequenceMethods Slice_file_sequence_methods =
{
    Slice_file_length,
    0,
    0,
    Slice_file_item,
    0,
    Slice_file_ass_item,
    0
};
*/

static PyTypeObject Slice_file_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "SliceFile", /* name */
    sizeof(struct Py_Slice_file), /* basicsize */
    0, /* itemsize */
    delete_py_slice_file, /* destructor */
    0, /* printfunc */
    getattr_py_slice_file, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Slice_file_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Slice_file(PyObject *self, PyObject *args)
{
    int dim, orient;
    PyObject *block_file_obj, *mem_cache_obj, *obj;

    if (!PyArg_ParseTuple(args, "iiOO", &orient, &dim,
                          &block_file_obj, &mem_cache_obj))
        RETURN_OBJ_ERROR("must have three arguments: dim, block_file, mem_cache");

    obj = new_py_slice_file(orient, dim, block_file_obj, mem_cache_obj);

    return obj;
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Slice_file_type_methods[] =
{
    { "SliceFile",	(PyCFunction) init_Py_Slice_file,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import SliceFile" in
* a Python program. The function is usually called "initSlice_file": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initSliceFile(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Slice_file_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("SliceFile", Slice_file_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("SliceFile.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);

    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module SliceFile");
}

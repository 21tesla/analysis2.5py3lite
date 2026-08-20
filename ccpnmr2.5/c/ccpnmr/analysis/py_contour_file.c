/*
======================COPYRIGHT/LICENSE START==========================

py_contour_file.c: Part of the CcpNmr Analysis program

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
#include "py_contour_file.h"

#include "py_block_file.h"
#include "py_contour_levels.h"
#include "py_contour_style.h"
#include "py_draw_handler.h"
#include "py_mem_cache.h"
#include "py_store_file.h"
#include "python_util.h"
#include "utility.h"

/* TBD: temporary until remove static variables */
#include "contour_data.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Contour_file_type;

Bool is_py_contour_file(PyObject *obj)
{
    /*  below does not work because different *.so files end up
        with different addresses for Contour_file_type
        return (obj->ob_type == &Contour_file_type);
    */
    return valid_py_object(obj, &Contour_file_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *draw(PyObject *self, PyObject *args)
{
    int n, ndim, first[MAX_NDIM], last[MAX_NDIM], ncomponents = 0, *components = NULL;
    PyObject *first_obj, *last_obj, *contour_levels_obj, *contour_style_obj;
    PyObject *handler_obj, *components_obj = NULL;
    Py_Contour_file obj = (Py_Contour_file) self;
    Contour_file contour_file = obj->contour_file;
    Py_draw_handler py_draw_handler;
    Drawing_funcs *drawing_funcs;
    Generic_ptr handler;
    CcpnStatus status;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "OOOOO|O", &handler_obj, &first_obj, &last_obj,
                          &contour_levels_obj, &contour_style_obj, &components_obj))
        RETURN_OBJ_ERROR("need arguments: handler, first, last, contour_levels, contour_style [, components ]");

    py_draw_handler = new_py_draw_handler(handler_obj);
    if (!py_draw_handler)
        RETURN_OBJ_ERROR("first argument must be handler object");

    handler = py_draw_handler->handler;
    drawing_funcs = py_draw_handler->drawing_funcs;
    delete_py_draw_handler(py_draw_handler);

    if (contour_file->block_file)
        ndim = contour_file->block_file->ndim;
    else
        ndim = contour_file->store_file->ndim;

    if ((get_python_int_array(first_obj, MAX_NDIM, &n, first, error_msg) == CCPN_ERROR)
            || (n != ndim))
    {
        /*
        printf("contour_file->block_file = %d\n", contour_file->block_file);
        printf("contour_file->store_file = %d\n", contour_file->store_file);
        printf("n = %d\n", n);
        printf("ndim = %d\n", ndim);
        printf("xdim = %d\n", contour_file->store_file->xdim);
        printf("ydim = %d\n", contour_file->store_file->ydim);
        */
        sprintf(error_msg,
                "second argument, first, must be list or tuple of size %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_int_array(last_obj, MAX_NDIM, &n, last, error_msg) == CCPN_ERROR)
            || (n != ndim))
    {
        sprintf(error_msg,
                "third argument, last, must be list or tuple of size %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (!is_py_contour_levels(contour_levels_obj))
        RETURN_OBJ_ERROR("fourth argument must be contour_levels object");

    if (!is_py_contour_style(contour_style_obj))
        RETURN_OBJ_ERROR("fifth argument must be contour_style object");

    if (components_obj && (components_obj != Py_None))
    {
        if (get_python_int_alloc_array(components_obj, &ncomponents, &components, error_msg) == CCPN_ERROR)
            RETURN_OBJ_ERROR(error_msg);
    }

    /*  TODO: abort_funcs  */
    status = region_contour_file(contour_file, first, last, no_abort_func,
                                 ((Py_Contour_levels) contour_levels_obj)->contour_levels,
                                 ((Py_Contour_style) contour_style_obj)->contour_style,
                                 ncomponents, components, drawing_funcs, handler, error_msg);

    FREE(components, int);

    if (status == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "draw",		draw,			METH_VARARGS },
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_contour_file(int xdim, int ydim,
                                     Block_file block_file, Store_file store_file,
                                     PyObject *mem_cache_obj, Bool transposed)
{
    Contour_file contour_file;
    Py_Contour_file obj;
    Line error_msg;

    if (!is_py_mem_cache(mem_cache_obj))
        RETURN_OBJ_ERROR("must pass mem_cache object");

    contour_file = new_contour_file(xdim, ydim, block_file, store_file,
                                    ((Py_Mem_cache) mem_cache_obj)->mem_cache, transposed, error_msg);

    if (!contour_file)
        RETURN_OBJ_ERROR(error_msg);

    PY_MALLOC(obj, struct Py_Contour_file, &Contour_file_type);

    if (!obj)
    {
        delete_contour_file(contour_file);

        RETURN_OBJ_ERROR("allocating Py_Contour_file object");
    }

    obj->contour_file = contour_file;

    return (PyObject *) obj;
}

static void delete_py_contour_file(PyObject *self)
{
    Py_Contour_file obj = (Py_Contour_file) self;
    Contour_file contour_file = obj->contour_file;

    /*
        printf("in delete_py_contour_file\n");
    */

    delete_contour_file(contour_file);

    PY_FREE(self);
}

/*
static int print_py_contour_file(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_contour_file(PyObject *self, char *name)
{
    Py_Contour_file obj = (Py_Contour_file) self;
    Contour_file contour_file = obj->contour_file;

    if (equal_strings(name, "xdim"))
    {
        return Py_BuildValue("i", contour_file->xdim);
    }
    else if (equal_strings(name, "ydim"))
    {
        return Py_BuildValue("i", contour_file->ydim);
    }
    else if (equal_strings(name, "have_neg"))
    {
        if (contour_file->store_file)
            return Py_BuildValue("i", contour_file->store_file->have_neg);
        else if (contour_file->contour_levels)
            return Py_BuildValue("i", have_neg_contour_levels(contour_file->contour_levels));
        else
            RETURN_OBJ_ERROR("do not have valid attribute have_neg");
    }
    else if (equal_strings(name, "have_pos"))
    {
        if (contour_file->store_file)
            return Py_BuildValue("i", contour_file->store_file->have_pos);
        else if (contour_file->contour_levels)
            return Py_BuildValue("i", have_pos_contour_levels(contour_file->contour_levels));
        else
            RETURN_OBJ_ERROR("do not have valid attribute have_pos");
    }
    else
    {
        return Py_FindMethod(py_handler_methods, self, name);
    }
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Contour_file_sequence_methods =
{
    Contour_file_length,
    Contour_file_concat,
    Contour_file_repeat,
    Contour_file_item,
    Contour_file_slice,
    Contour_file_ass_item,
    Contour_file_ass_slice
};

static PySequenceMethods Contour_file_sequence_methods =
{
    Contour_file_length,
    0,
    0,
    Contour_file_item,
    0,
    Contour_file_ass_item,
    0
};
*/

static PyTypeObject Contour_file_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "ContourFile", /* name */
    sizeof(struct Py_Contour_file), /* basicsize */
    0, /* itemsize */
    delete_py_contour_file, /* destructor */
    0, /* printfunc */
    getattr_py_contour_file, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Contour_file_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Contour_file(PyObject *self, PyObject *args)
{
    int xdim, ydim;
    PyObject *block_file_obj, *mem_cache_obj, *obj;
    Block_file block_file;

    if (!PyArg_ParseTuple(args, "iiOO", &xdim, &ydim,
                          &block_file_obj, &mem_cache_obj))
        RETURN_OBJ_ERROR("must have arguments: xdim, ydim, block_file, mem_cache");

    if (!is_py_some_block_file(block_file_obj))
        RETURN_OBJ_ERROR("must pass block_file object");

    if (is_py_block_file(block_file_obj))
        block_file = ((Py_Block_file) block_file_obj)->block_file;
    else
        block_file = (Block_file) (((Py_Shape_block_file) block_file_obj)->shape_block_file);

    obj = new_py_contour_file(xdim, ydim, block_file, NULL,
                              mem_cache_obj, CCPN_FALSE);

    return obj;
}

static PyObject *init_Py_Stored_Contour_file(PyObject *self, PyObject *args)
{
    int xdim, ydim, transposed = 0;
    Store_file store_file;
    PyObject *store_file_obj, *mem_cache_obj, *obj;

    if (!PyArg_ParseTuple(args, "OO|i",
                          &store_file_obj, &mem_cache_obj, &transposed))
        RETURN_OBJ_ERROR("must have arguments: store_file, mem_cache [ transposed ]");

    if (!is_py_store_file(store_file_obj))
        RETURN_OBJ_ERROR("must pass store_file object");

    store_file = ((Py_Store_file) store_file_obj)->store_file;

    obj = new_py_contour_file(store_file->xdim, store_file->ydim,
                              NULL, store_file, mem_cache_obj, (Bool) transposed);

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


static struct PyMethodDef Contour_file_type_methods[] =
{
    { "ContourFile",	(PyCFunction) init_Py_Contour_file,	METH_VARARGS },
    { "StoredContourFile",	(PyCFunction) init_Py_Stored_Contour_file,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import ContourFile" in
* a Python program. The function is usually called "initContour_file": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initContourFile(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Contour_file_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("ContourFile", Contour_file_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("ContourFile.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);

    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module ContourFile");
}

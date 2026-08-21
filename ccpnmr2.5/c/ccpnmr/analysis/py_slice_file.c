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

    obj = (Py_Slice_file) PyObject_New(struct Py_Slice_file, &Slice_file_type);

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

    Py_TYPE(self)->tp_free(self);
}

/*
static int print_py_slice_file(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_slice_file(PyObject *self, PyObject *attr_name)
{
    Py_Slice_file obj = (Py_Slice_file) self;
    Slice_file slice_file = obj->slice_file;

    if (PyUnicode_Check(attr_name))
    {
        const char *name = PyUnicode_AsUTF8(attr_name);
        if (name && strcmp(name, "dim") == 0)
            return Py_BuildValue("i", slice_file->dim);
    }

    return PyObject_GenericGetAttr(self, attr_name);
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
    PyVarObject_HEAD_INIT(NULL, 0)
    "SliceFile",                                 /* tp_name */
    sizeof(struct Py_Slice_file),                 /* tp_basicsize */
    0,                                            /* tp_itemsize */
    (destructor) delete_py_slice_file,            /* tp_dealloc */
    0,                                            /* tp_vectorcall */
    0,                                            /* tp_getattr */
    0,                                            /* tp_setattr */
    0,                                            /* tp_as_async */
    0,                                            /* tp_repr */
    0,                                            /* tp_as_number */
    0,                                            /* tp_as_sequence */
    0,                                            /* tp_as_mapping */
    0,                                            /* tp_hash */
    0,                                            /* tp_call */
    0,                                            /* tp_str */
    (getattrofunc) getattr_py_slice_file,         /* tp_getattro */
    0,                                            /* tp_setattro */
    0,                                            /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                           /* tp_flags */
    "SliceFile -- contour slice file",            /* tp_doc */
    0,                                            /* tp_traverse */
    0,                                            /* tp_clear */
    0,                                            /* tp_richcompare */
    0,                                            /* tp_weaklistoffset */
    0,                                            /* tp_iter */
    0,                                            /* tp_iternext */
    py_handler_methods,                           /* tp_methods */
    0,                                            /* tp_members */
    0,                                            /* tp_getset */
    0,                                            /* tp_base */
    0,                                            /* tp_dict */
    0,                                            /* tp_descr_get */
    0,                                            /* tp_descr_set */
    0,                                            /* tp_dictoffset */
    0,                                            /* tp_init */
    0,                                            /* tp_alloc */
    0,                                            /* tp_new */
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

static struct PyModuleDef slice_file_module_def =
{
    PyModuleDef_HEAD_INIT,
    "SliceFile",
    "CCPNMR SliceFile module (Python 3 compatible)",
    -1,
    Slice_file_type_methods
};

PyMODINIT_FUNC PyInit_SliceFile(void)
{
    if (PyType_Ready(&Slice_file_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&slice_file_module_def);
    if (!m)
        return NULL;

    ErrorObject = PyErr_NewException("SliceFile.error", NULL, NULL);
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

/*
======================COPYRIGHT/LICENSE START==========================

py_contour_style.c: Part of the CcpNmr Analysis program

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
#include "py_contour_style.h"

#include "python_util.h"

#define  NCOLORS  3

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Contour_style_type;

Bool is_py_contour_style(PyObject *obj)
{
    /*  below does not work because different *.so files end up
        with different addresses for Contour_style_type
        return (obj->ob_type == &Contour_style_type);
    */
    return valid_py_object(obj, &Contour_style_type);
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

static PyObject *new_py_contour_style(PyObject *pos_color_obj,
                                      PyObject *neg_color_obj, int pos_line_style, int neg_line_style)
{
    int n, npos_colors, nneg_colors;
    float **pos_colors, **neg_colors;
    Contour_style contour_style;
    Py_Contour_style obj;
    Line error_msg;

    if ((get_python_float_alloc_matrix(pos_color_obj, &npos_colors, &n,
                                       &pos_colors, error_msg) == CCPN_ERROR)
            || (n != NCOLORS))
    {
        sprintf(error_msg, "pos_colors must be list or tuple with elements of size %d", NCOLORS);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_float_alloc_matrix(neg_color_obj, &nneg_colors, &n,
                                       &neg_colors, error_msg) == CCPN_ERROR)
            || (n != NCOLORS))
    {
        sprintf(error_msg, "neg_colors must be list or tuple with elements of size %d", NCOLORS);
        RETURN_OBJ_ERROR(error_msg);
    }

    contour_style = new_contour_style(npos_colors, pos_colors, nneg_colors, neg_colors,
                                      pos_line_style, neg_line_style);

    if (!contour_style)
        RETURN_OBJ_ERROR("allocating Contour_style object");

    obj = (Py_Contour_style) PyObject_New(struct Py_Contour_style, &Contour_style_type);

    if (!obj)
    {
        delete_contour_style(contour_style);

        RETURN_OBJ_ERROR("allocating Py_Contour_style object");
    }

    obj->contour_style = contour_style;

    return (PyObject *) obj;
}

static void delete_py_contour_style(PyObject *self)
{
    Py_Contour_style obj = (Py_Contour_style) self;
    Contour_style contour_style = obj->contour_style;

    /*
        printf("in delete_py_contour_style\n");
    */

    delete_contour_style(contour_style);

    Py_TYPE(self)->tp_free(self);
}

/*
static int print_py_contour_style(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_contour_style(PyObject *self, PyObject *attr_name)
{
    return PyObject_GenericGetAttr(self, attr_name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Contour_style_sequence_methods =
{
    Contour_style_length,
    Contour_style_concat,
    Contour_style_repeat,
    Contour_style_item,
    Contour_style_slice,
    Contour_style_ass_item,
    Contour_style_ass_slice
};

static PySequenceMethods Contour_style_sequence_methods =
{
    Contour_style_length,
    0,
    0,
    Contour_style_item,
    0,
    Contour_style_ass_item,
    0
};
*/

static PyTypeObject Contour_style_type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    "ContourStyle",                         /* tp_name */
    sizeof(struct Py_Contour_style),        /* tp_basicsize */
    0,                                      /* tp_itemsize */
    (destructor) delete_py_contour_style,   /* tp_dealloc */
    0,                                      /* tp_vectorcall */
    0,                                      /* tp_getattr */
    0,                                      /* tp_setattr */
    0,                                      /* tp_as_async */
    0,                                      /* tp_repr */
    0,                                      /* tp_as_number */
    0,                                      /* tp_as_sequence */
    0,                                      /* tp_as_mapping */
    0,                                      /* tp_hash */
    0,                                      /* tp_call */
    0,                                      /* tp_str */
    (getattrofunc) getattr_py_contour_style, /* tp_getattro */
    0,                                      /* tp_setattro */
    0,                                      /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                     /* tp_flags */
    "ContourStyle -- contour style",        /* tp_doc */
    0,                                      /* tp_traverse */
    0,                                      /* tp_clear */
    0,                                      /* tp_richcompare */
    0,                                      /* tp_weaklistoffset */
    0,                                      /* tp_iter */
    0,                                      /* tp_iternext */
    py_handler_methods,                     /* tp_methods */
    0,                                      /* tp_members */
    0,                                      /* tp_getset */
    0,                                      /* tp_base */
    0,                                      /* tp_dict */
    0,                                      /* tp_descr_get */
    0,                                      /* tp_descr_set */
    0,                                      /* tp_dictoffset */
    0,                                      /* tp_init */
    0,                                      /* tp_alloc */
    0,                                      /* tp_new */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Contour_style(PyObject *self, PyObject *args)
{
    int pos_line_style, neg_line_style;
    PyObject *pos_color_obj, *neg_color_obj, *obj;

    if (!PyArg_ParseTuple(args, "OOii", &pos_color_obj, &neg_color_obj,
                          &pos_line_style, &neg_line_style))
        RETURN_OBJ_ERROR("must have four arguments: pos color, neg color, pos line style, neg line style");

    obj = new_py_contour_style(pos_color_obj, neg_color_obj, pos_line_style, neg_line_style);

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


static struct PyMethodDef Contour_style_type_methods[] =
{
    { "ContourStyle",	(PyCFunction) init_Py_Contour_style,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import ContourStyle" in
* a Python program. The function is usually called "initContour_style": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

static struct PyModuleDef ContourStyle_module_def =
{
    PyModuleDef_HEAD_INIT,
    "ContourStyle",
    "CCPNMR ContourStyle module (Python 3 compatible)",
    -1,
    Contour_style_type_methods
};

PyMODINIT_FUNC PyInit_ContourStyle(void)
{
    if (PyType_Ready(&Contour_style_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&ContourStyle_module_def);
    if (!m)
        return NULL;

    ErrorObject = PyErr_NewException("ContourStyle.error", NULL, NULL);
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

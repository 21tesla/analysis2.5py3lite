
/*
======================COPYRIGHT/LICENSE START==========================

py_ps_handler.c: Part of the CcpNmr Analysis program

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
#include "py_ps_handler.h"

#include "python_util.h"

#define  NCOLORS  3

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Ps_handler_type;

Bool is_py_ps_handler(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Ps_handler_type
    return (obj->ob_type == &Ps_handler_type);
*/
    return valid_py_object(obj, &Ps_handler_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *newRange(PyObject *self, PyObject *args)
{
    Py_Ps_handler obj = (Py_Ps_handler) self;
    Ps_handler ps_handler = obj->ps_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    new_range_ps_handler(ps_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *clipRange(PyObject *self, PyObject *args)
{
    Py_Ps_handler obj = (Py_Ps_handler) self;
    Ps_handler ps_handler = obj->ps_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    clip_range_ps_handler(ps_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawText(PyObject *self, PyObject *args)
{
    Py_Ps_handler obj = (Py_Ps_handler) self;
    Ps_handler ps_handler = obj->ps_handler;
    float x, y, a, b;
    CcpnString text;

    if (!PyArg_ParseTuple(args, "sffff", &text, &x, &y, &a, &b))
        RETURN_OBJ_ERROR("need five arguments: text, x, y, a, b");

    draw_text_ps_handler(ps_handler, text, x, y, a, b);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawLine(PyObject *self, PyObject *args)
{
    Py_Ps_handler obj = (Py_Ps_handler) self;
    Ps_handler ps_handler = obj->ps_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    draw_line_ps_handler(ps_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawClippedLine(PyObject *self, PyObject *args)
{
    Py_Ps_handler obj = (Py_Ps_handler) self;
    Ps_handler ps_handler = obj->ps_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    draw_clipped_line_ps_handler(ps_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawDashLine(PyObject *self, PyObject *args)
{
    Py_Ps_handler obj = (Py_Ps_handler) self;
    Ps_handler ps_handler = obj->ps_handler;
    float x0, y0, x1, y1;
    int dash_length, gap_length;

    if (!PyArg_ParseTuple(args, "ffffii", &x0, &y0, &x1, &y1, &dash_length, &gap_length))
        RETURN_OBJ_ERROR("need six arguments: x0, y0, x1, y1, dash_length, gap_length");

    if (dash_length < 0)
        RETURN_OBJ_ERROR("dash_length < 0");

    if (gap_length < 0)
        RETURN_OBJ_ERROR("gap_length < 0");

    draw_dash_line_ps_handler(ps_handler, x0, y0, x1, y1, dash_length, gap_length);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawDashBox(PyObject *self, PyObject *args)
{
    Py_Ps_handler obj = (Py_Ps_handler) self;
    Ps_handler ps_handler = obj->ps_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    draw_dash_box_ps_handler(ps_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setColor(PyObject *self, PyObject *args)
{
    int n;
    float color[NCOLORS];
    PyObject *color_obj;
    Py_Ps_handler obj = (Py_Ps_handler) self;
    Ps_handler ps_handler = obj->ps_handler;
    Line error_msg;
 
    if (!PyArg_ParseTuple(args, "O", &color_obj))
	RETURN_OBJ_ERROR("need one argument: color");
 
    if ((get_python_float_array(color_obj, NCOLORS, &n, color,
						error_msg) == CCPN_ERROR)
	|| (n != NCOLORS))
    {
	sprintf(error_msg, "color must be list or tuple of size %d", NCOLORS);
	RETURN_OBJ_ERROR(error_msg);
    }

    set_color_ps_handler(ps_handler, color);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setBlack(PyObject *self, PyObject *args)
{
    Py_Ps_handler obj = (Py_Ps_handler) self;
    Ps_handler ps_handler = obj->ps_handler;

    set_black_ps_handler(ps_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setLineWidth(PyObject *self, PyObject *args)
{
    Py_Ps_handler obj = (Py_Ps_handler) self;
    Ps_handler ps_handler = obj->ps_handler;
    float line_width;

    if (!PyArg_ParseTuple(args, "f", &line_width))
	RETURN_OBJ_ERROR("need one argument: linewidth");
 
    set_line_width_ps_handler(ps_handler, line_width);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setFont(PyObject *self, PyObject *args)
{
    Py_Ps_handler obj = (Py_Ps_handler) self;
    Ps_handler ps_handler = obj->ps_handler;
    int size;
    CcpnString name;

    if (!PyArg_ParseTuple(args, "si", &name, &size))
        RETURN_OBJ_ERROR("need two arguments: name, size");

    if (size < 1)
        RETURN_OBJ_ERROR("font size needs to be at least 1");

    set_font_ps_handler(ps_handler, name, size);

/*
    Py_DECREF(name);
*/

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "newRange",	newRange,		METH_VARARGS },
    { "clipRange",	clipRange,		METH_VARARGS },
    { "drawText",	drawText,		METH_VARARGS },
    { "drawLine",	drawLine,		METH_VARARGS },
    { "drawClippedLine",drawClippedLine,	METH_VARARGS },
    { "drawDashLine",	drawDashLine,		METH_VARARGS },
    { "drawDashBox",	drawDashBox,		METH_VARARGS },
    { "setColor",	setColor,		METH_VARARGS },
    { "setBlack",	setBlack,		METH_VARARGS },
    { "setLineWidth",	setLineWidth,		METH_VARARGS },
    { "setFont",	setFont,		METH_VARARGS },
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_ps_handler(CcpnString *file_name, float width, float height,
                                                        CcpnString output_style)
{
    Py_Ps_handler obj;
    Ps_handler ps_handler;
//    FILE *fp;

//    if (!PyFile_Check(fp_obj))
//	RETURN_OBJ_ERROR("argument not a Python file object");
 
//    fp = PyFile_AsFile(fp_obj);

    ps_handler = new_ps_handler(file_name, width, height, output_style);

    if (!ps_handler)
	 RETURN_OBJ_ERROR("allocating Ps_handler object");

    obj = (Py_Ps_handler) PyObject_New(struct Py_Ps_handler, &Ps_handler_type);

    if (!obj)
    {
	delete_ps_handler(ps_handler);

	RETURN_OBJ_ERROR("allocating Py_Ps_handler object");
    }

    obj->ps_handler = ps_handler;

    return (PyObject *) obj;
}

static void delete_py_ps_handler(PyObject *self)
{
    Py_Ps_handler obj = (Py_Ps_handler) self;
    Ps_handler ps_handler = obj->ps_handler;

    delete_ps_handler(ps_handler);

    Py_TYPE(self)->tp_free(self);
}

/*
static int print_py_ps_handler(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_ps_handler(PyObject *self, PyObject *attr_name)
{
    return PyObject_GenericGetAttr(self, attr_name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Ps_handler_sequence_methods =
{
    Ps_handler_length,
    Ps_handler_concat,
    Ps_handler_repeat,
    Ps_handler_item,
    Ps_handler_slice,
    Ps_handler_ass_item,
    Ps_handler_ass_slice
};

static PySequenceMethods Ps_handler_sequence_methods =
{
    Ps_handler_length,
    0,
    0,
    Ps_handler_item,
    0,
    Ps_handler_ass_item,
    0
};
*/

static PyTypeObject Ps_handler_type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    "PsHandler",                              /* tp_name */
    sizeof(struct Py_Ps_handler),             /* tp_basicsize */
    0,                                        /* tp_itemsize */
    (destructor) delete_py_ps_handler,        /* tp_dealloc */
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
    (getattrofunc) getattr_py_ps_handler,     /* tp_getattro */
    0,                                        /* tp_setattro */
    0,                                        /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                       /* tp_flags */
    "PsHandler -- NMR PostScript drawing",    /* tp_doc */
    0,                                        /* tp_traverse */
    0,                                        /* tp_clear */
    0,                                        /* tp_richcompare */
    0,                                        /* tp_weaklistoffset */
    0,                                        /* tp_iter */
    0,                                        /* tp_iternext */
    py_handler_methods,                       /* tp_methods */
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

static PyObject *init_Py_Ps_handler(PyObject *self, PyObject *args)
{
    float width, height;
    CcpnString output_style;
//    PyObject *fp_obj;
    CcpnString file_name;

//    if (!PyArg_ParseTuple(args, "Offs", &fp_obj, &width, &height, &output_style))
//        RETURN_OBJ_ERROR("must have four arguments: stream width height output_style");
    if (!PyArg_ParseTuple(args, "sffs", &file_name, &width, &height, &output_style))
        RETURN_OBJ_ERROR("must have four arguments: stream width height output_style");

    return new_py_ps_handler(file_name, width, height, output_style);
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Ps_handler_type_methods[] =
{
    { "PsHandler",	(PyCFunction) init_Py_Ps_handler,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import PsHandler" in 
* a Python program. The function is usually called "initPs_handler": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

static struct PyModuleDef ps_handler_module_def =
{
    PyModuleDef_HEAD_INIT,
    "PsHandler",
    "CCPNMR Ps Handler module (Python 3 compatible)",
    -1,
    Ps_handler_type_methods
};

PyMODINIT_FUNC PyInit_PsHandler(void)
{
    if (PyType_Ready(&Ps_handler_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&ps_handler_module_def);
    if (!m)
        return NULL;

    ErrorObject = PyErr_NewException("PsHandler.error", NULL, NULL);
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

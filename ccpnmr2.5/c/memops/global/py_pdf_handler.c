
/*
======================COPYRIGHT/LICENSE START==========================

py_pdf_handler.c: Part of the CcpNmr Analysis program

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
#include "py_pdf_handler.h"

#include "python_util.h"

#define  NCOLORS  3

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Pdf_handler_type;

Bool is_py_pdf_handler(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Pdf_handler_type
    return (obj->ob_type == &Pdf_handler_type);
*/
    return valid_py_object(obj, &Pdf_handler_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *newRange(PyObject *self, PyObject *args)
{
    Py_Pdf_handler obj = (Py_Pdf_handler) self;
    Pdf_handler pdf_handler = obj->pdf_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    new_range_pdf_handler(pdf_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *clipRange(PyObject *self, PyObject *args)
{
    Py_Pdf_handler obj = (Py_Pdf_handler) self;
    Pdf_handler pdf_handler = obj->pdf_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    clip_range_pdf_handler(pdf_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawText(PyObject *self, PyObject *args)
{
    Py_Pdf_handler obj = (Py_Pdf_handler) self;
    Pdf_handler pdf_handler = obj->pdf_handler;
    float x, y, a, b;
    CcpnString text;

    if (!PyArg_ParseTuple(args, "sffff", &text, &x, &y, &a, &b))
        RETURN_OBJ_ERROR("need five arguments: text, x, y, a, b");

    draw_text_pdf_handler(pdf_handler, text, x, y, a, b);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawLine(PyObject *self, PyObject *args)
{
    Py_Pdf_handler obj = (Py_Pdf_handler) self;
    Pdf_handler pdf_handler = obj->pdf_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    draw_line_pdf_handler(pdf_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawDashLine(PyObject *self, PyObject *args)
{
    Py_Pdf_handler obj = (Py_Pdf_handler) self;
    Pdf_handler pdf_handler = obj->pdf_handler;
    float x0, y0, x1, y1;
    int dash_length, gap_length;

    if (!PyArg_ParseTuple(args, "ffffii", &x0, &y0, &x1, &y1, &dash_length, &gap_length))
        RETURN_OBJ_ERROR("need six arguments: x0, y0, x1, y1, dash_length, gap_length");

    if (dash_length < 0)
        RETURN_OBJ_ERROR("dash_length < 0");

    if (gap_length < 0)
        RETURN_OBJ_ERROR("gap_length < 0");

    draw_dash_line_pdf_handler(pdf_handler, x0, y0, x1, y1, dash_length, gap_length);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawDashBox(PyObject *self, PyObject *args)
{
    Py_Pdf_handler obj = (Py_Pdf_handler) self;
    Pdf_handler pdf_handler = obj->pdf_handler;
    float x0, y0, x1, y1;
 
    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");
 
    draw_dash_box_pdf_handler(pdf_handler, x0, y0, x1, y1);
 
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setColor(PyObject *self, PyObject *args)
{
    int n;
    float color[NCOLORS];
    PyObject *color_obj;
    Py_Pdf_handler obj = (Py_Pdf_handler) self;
    Pdf_handler pdf_handler = obj->pdf_handler;
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

    set_color_pdf_handler(pdf_handler, color);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setBlack(PyObject *self, PyObject *args)
{
    Py_Pdf_handler obj = (Py_Pdf_handler) self;
    Pdf_handler pdf_handler = obj->pdf_handler;

    set_black_pdf_handler(pdf_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setLineWidth(PyObject *self, PyObject *args)
{
    Py_Pdf_handler obj = (Py_Pdf_handler) self;
    Pdf_handler pdf_handler = obj->pdf_handler;
    float line_width;

    if (!PyArg_ParseTuple(args, "f", &line_width))
	RETURN_OBJ_ERROR("need one argument: linewidth");
 
    set_line_width_pdf_handler(pdf_handler, line_width);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setFont(PyObject *self, PyObject *args)
{
    Py_Pdf_handler obj = (Py_Pdf_handler) self;
    Pdf_handler pdf_handler = obj->pdf_handler;
    int size;
    CcpnString name;

    if (!PyArg_ParseTuple(args, "si", &name, &size))
        RETURN_OBJ_ERROR("need two arguments: name, size");

    if (size < 1)
        RETURN_OBJ_ERROR("font size needs to be at least 1");

    set_font_pdf_handler(pdf_handler, name, size);

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
    { "drawDashLine",	drawDashLine,		METH_VARARGS },
    { "drawDashBox",	drawDashBox,		METH_VARARGS },
    { "setColor",	setColor,		METH_VARARGS },
    { "setBlack",	setBlack,		METH_VARARGS },
    { "setLineWidth",	setLineWidth,		METH_VARARGS },
    { "setFont",        setFont,                METH_VARARGS },
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_pdf_handler(CcpnString *file_name, float width, float height,
                                                        CcpnString output_style)
{
    Py_Pdf_handler obj;
    Pdf_handler pdf_handler;
//    FILE *fp;

//    if (!PyFile_Check(fp_obj))
//	RETURN_OBJ_ERROR("argument not a Python file object");
 
//    fp = PyFile_AsFile(fp_obj);

    pdf_handler = new_pdf_handler(file_name, width, height, output_style);

    if (!pdf_handler)
	 RETURN_OBJ_ERROR("allocating Pdf_handler object");

    PY_MALLOC(obj, struct Py_Pdf_handler, &Pdf_handler_type);

    if (!obj)
    {
	delete_pdf_handler(pdf_handler);

	RETURN_OBJ_ERROR("allocating Py_Pdf_handler object");
    }

    obj->pdf_handler = pdf_handler;

    return (PyObject *) obj;
}

static void delete_py_pdf_handler(PyObject *self)
{
    Py_Pdf_handler obj = (Py_Pdf_handler) self;
    Pdf_handler pdf_handler = obj->pdf_handler;

    delete_pdf_handler(pdf_handler);

    PY_FREE(self);
}

/*
static int print_py_pdf_handler(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_pdf_handler(PyObject *self, char *name)
{
/*
    Pdf_handler *obj = (Pdf_handler *) self;
    Random_access *a = obj->py_handler;

    if (equal_strings(name, "par_file"))
	return Py_BuildValue("s", a->par_file);
    else if (equal_strings(name, "access_method"))
	return Py_BuildValue("s", access_method_name(a->access_method));
    else if (equal_strings(name, "data_format"))
	return get_Pdf_handler_format(a);
    else
*/
	return Py_FindMethod(py_handler_methods, self, name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Pdf_handler_sequence_methods =
{
    Pdf_handler_length,
    Pdf_handler_concat,
    Pdf_handler_repeat,
    Pdf_handler_item,
    Pdf_handler_slice,
    Pdf_handler_ass_item,
    Pdf_handler_ass_slice
};

static PySequenceMethods Pdf_handler_sequence_methods =
{
    Pdf_handler_length,
    0,
    0,
    Pdf_handler_item,
    0,
    Pdf_handler_ass_item,
    0
};
*/

static PyTypeObject Pdf_handler_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "PdfHandler", /* name */
    sizeof(struct Py_Pdf_handler), /* basicsize */
    0, /* itemsize */
    delete_py_pdf_handler, /* destructor */
    0, /* printfunc */
    getattr_py_pdf_handler, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Pdf_handler_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Pdf_handler(PyObject *self, PyObject *args)
{
    float width, height;
    CcpnString output_style;
//    PyObject *fp_obj;
    CcpnString file_name;

//    if (!PyArg_ParseTuple(args, "Offs", &fp_obj, &width, &height, &output_style))
//        RETURN_OBJ_ERROR("must have four arguments: stream width height output_style");
    if (!PyArg_ParseTuple(args, "sffs", &file_name, &width, &height, &output_style))
        RETURN_OBJ_ERROR("must have four arguments: stream width height output_style");

    return new_py_pdf_handler(file_name, width, height, output_style);
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Pdf_handler_type_methods[] =
{
    { "PdfHandler",	(PyCFunction) init_Py_Pdf_handler,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import PdfHandler" in 
* a Python program. The function is usually called "initPdf_handler": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initPdfHandler(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Pdf_handler_type.ob_type = &PyType_Type;
#endif

    /* create the module and add the functions */
    m = Py_InitModule("PdfHandler", Pdf_handler_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("PdfHandler.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module PdfHandler");
}

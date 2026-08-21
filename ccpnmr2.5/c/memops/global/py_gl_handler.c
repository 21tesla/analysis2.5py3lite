
/*
======================COPYRIGHT/LICENSE START==========================

py_gl_handler.c: Part of the CcpNmr Analysis program

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
#include "py_gl_handler.h"

#ifdef WIN64
#include <windows.h>
#endif /* end WIN32 */

#ifndef IGNORE_GL
#include <GL/gl.h>

#include "python_util.h"
#include "py_tk_util.h"

#define  NCOLORS  3

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Gl_handler_type;

Bool is_py_gl_handler(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Gl_handler_type
    return (obj->ob_type == &Gl_handler_type);
*/
    return valid_py_object(obj, &Gl_handler_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *resize(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    int w, h;

    if (!PyArg_ParseTuple(args, "ii", &w, &h))
        RETURN_OBJ_ERROR("need two arguments: w, h");

    resize_gl_handler(gl_handler, w, h);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *expose(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    int x, y, w, h;

    if (!PyArg_ParseTuple(args, "iiii", &x, &y, &w, &h))
        RETURN_OBJ_ERROR("need four arguments: x, y, w, h");

    expose_gl_handler(gl_handler, x, y, w, h);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *flush(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    flush_gl_handler(gl_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *mapRanges(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    float x0, y0, x1, y1, a0, b0, a1, b1;

    if (!PyArg_ParseTuple(args, "ffffffff", &x0, &y0, &x1, &y1, &a0, &b0, &a1, &b1))
        RETURN_OBJ_ERROR("need eight arguments: x0, y0, x1, y1, a0, b0, a1, b1");

    map_ranges_gl_handler(gl_handler, x0, y0, x1, y1, a0, b0, a1, b1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *makeCurrent(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    return Py_BuildValue("i", (int) make_current_gl_handler(gl_handler));
}

static PyObject *swapBuffers(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    swap_buffers_gl_handler(gl_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawBox(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    draw_box_gl_handler(gl_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawDashBox(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    draw_dash_box_gl_handler(gl_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawXorBox(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    draw_xor_box_gl_handler(gl_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *startXor(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    start_xor_gl_handler(gl_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *finishXor(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    finish_xor_gl_handler(gl_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *resetXor(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    reset_xor_gl_handler(gl_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *clearXor(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    clear_xor_gl_handler(gl_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *startFront(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    start_front_gl_handler(gl_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *startBack(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    start_back_gl_handler(gl_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *endBack(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    end_back_gl_handler(gl_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setBackground(PyObject *self, PyObject *args)
{
    int n;
    float background[NCOLORS];
    PyObject *background_obj;
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    Line error_msg;
 
    if (!PyArg_ParseTuple(args, "O", &background_obj))
	RETURN_OBJ_ERROR("need one argument: background");
 
    if ((get_python_float_array(background_obj, NCOLORS, &n, background,
						error_msg) == CCPN_ERROR)
	|| (n != NCOLORS))
    {
	sprintf(error_msg, "background must be list or tuple of size %d", NCOLORS);
	RETURN_OBJ_ERROR(error_msg);
    }

    set_background_gl_handler(gl_handler, background);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawText(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    float x, y, a, b;
    CcpnString text;

    if (!PyArg_ParseTuple(args, "sffff", &text, &x, &y, &a, &b))
        RETURN_OBJ_ERROR("need five arguments: text, x, y, a, b");

    draw_text_gl_handler(gl_handler, text, x, y, a, b);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawLine(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    draw_line_gl_handler(gl_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawClippedLine(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    float x0, y0, x1, y1;

    if (!PyArg_ParseTuple(args, "ffff", &x0, &y0, &x1, &y1))
        RETURN_OBJ_ERROR("need four arguments: x0, y0, x1, y1");

    draw_clipped_line_gl_handler(gl_handler, x0, y0, x1, y1);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawDashLine(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    float x0, y0, x1, y1;
    int dash_length, gap_length;

    if (!PyArg_ParseTuple(args, "ffffii", &x0, &y0, &x1, &y1, &dash_length, &gap_length))
        RETURN_OBJ_ERROR("need six arguments: x0, y0, x1, y1, dash_length, gap_length");

    if (dash_length < 0)
        RETURN_OBJ_ERROR("dash_length < 0");

    if (gap_length < 0)
        RETURN_OBJ_ERROR("gap_length < 0");

    draw_dash_line_gl_handler(gl_handler, x0, y0, x1, y1, dash_length, gap_length);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setColor(PyObject *self, PyObject *args)
{
    int n;
    float color[NCOLORS];
    PyObject *color_obj;
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
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

    set_color_gl_handler(gl_handler, color);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setBlack(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    if (!PyArg_ParseTuple(args, ""))
	RETURN_OBJ_ERROR("should have no arguments");

    set_black_gl_handler(gl_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setWhite(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    if (!PyArg_ParseTuple(args, ""))
	RETURN_OBJ_ERROR("should have no arguments");

    set_white_gl_handler(gl_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setLineWidth(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    float line_width;

    if (!PyArg_ParseTuple(args, "f", &line_width))
	RETURN_OBJ_ERROR("need one argument: linewidth");
 
    set_line_width_gl_handler(gl_handler, line_width);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *resetLineWidth(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    reset_line_width_gl_handler(gl_handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *hasFont(PyObject *self)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    int size;
    Line name;

    if (has_font_gl_handler(gl_handler, &name, &size) == CCPN_TRUE) {
//        Py_INCREF(Py_True);
//        return Py_True;
        return Py_BuildValue("(si)", name, size);
    }
    else {
        Py_INCREF(Py_None);
        return Py_None;
    }
}

static PyObject *setFont(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    int size;
    CcpnString name;

    if (!PyArg_ParseTuple(args, "si", &name, &size))
        RETURN_OBJ_ERROR("need two arguments: name, size");

    if (size < 1)
        RETURN_OBJ_ERROR("font size needs to be at least 1");

    if (set_font_gl_handler(gl_handler, name, size) == CCPN_ERROR)
        RETURN_OBJ_ERROR("setting font");

/*
    Py_DECREF(name);
*/

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *getRegion(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    float x0, x1, y0, y1;

    if (!PyArg_ParseTuple(args, ""))
        RETURN_OBJ_ERROR("should have no arguments");

    get_region_gl_handler(gl_handler, &x0, &y0, &x1, &y1);

    return Py_BuildValue("(ffff)", x0, y0, x1, y1);
}

static PyObject *setIsDoubleBuffer(PyObject *self, PyObject *args)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;
    Bool is_double_buffer;

    if (!PyArg_ParseTuple(args, "i", &is_double_buffer))
        RETURN_OBJ_ERROR("need one argument: is_double_buffer");

    set_is_double_buffer_gl_handler(gl_handler, is_double_buffer);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "resize",		resize,			METH_VARARGS },
    { "expose",		expose,			METH_VARARGS },
    { "flush",		flush,			METH_VARARGS },
    { "mapRanges",	mapRanges,		METH_VARARGS },
    { "makeCurrent",	makeCurrent,		METH_VARARGS },
    { "swapBuffers",	swapBuffers,		METH_VARARGS },
    { "drawBox",	drawBox,		METH_VARARGS },
    { "drawDashBox",	drawDashBox,		METH_VARARGS },
    { "drawXorBox",	drawXorBox,		METH_VARARGS },
    { "startXor",	startXor,		METH_VARARGS },
    { "finishXor",	finishXor,		METH_VARARGS },
    { "resetXor",	resetXor,		METH_VARARGS },
    { "clearXor",	clearXor,		METH_VARARGS },
    { "startFront",	startFront,		METH_VARARGS },
    { "startBack",	startBack,		METH_VARARGS },
    { "endBack",	endBack,		METH_VARARGS },
    { "setBackground",	setBackground,		METH_VARARGS },
    { "drawText",	drawText,		METH_VARARGS },
    { "drawLine",	drawLine,		METH_VARARGS },
    { "drawClippedLine",drawClippedLine,	METH_VARARGS },
    { "drawDashLine",	drawDashLine,		METH_VARARGS },
    { "setColor",	setColor,		METH_VARARGS },
    { "setBlack",	setBlack,		METH_VARARGS },
    { "setWhite",	setWhite,		METH_VARARGS },
    { "setLineWidth",	setLineWidth,		METH_VARARGS },
    { "resetLineWidth",	resetLineWidth,		METH_VARARGS },
    { "hasFont",        hasFont,                METH_VARARGS },
    { "setFont",        setFont,                METH_VARARGS },
    { "getRegion",      getRegion,              METH_VARARGS },
    { "setIsDoubleBuffer",      setIsDoubleBuffer,              METH_VARARGS },
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_gl_handler(PyObject *widget, Bool direct, Bool dbl_buff)
{
    Py_Gl_handler obj;
    Gl_handler gl_handler;
    Tcl_Interp *tcl_interp;
    Tk_Window tk_display_win;
    Line error_msg;

/*  cannot do this when using Anaconda Python, it doesn't work
    if (!PyInstance_Check(widget))
	    RETURN_OBJ_ERROR("argument not a Python object");
*/
    tcl_interp = get_tcl_interp(widget, error_msg);
    if (!tcl_interp)
        RETURN_OBJ_ERROR(error_msg);

    tk_display_win = get_tk_window(widget, tcl_interp, error_msg);
    if (!tk_display_win)
        RETURN_OBJ_ERROR(error_msg);

    // NOTE:ED - this is where it crashes glX functions have bad display, segfault results

    gl_handler = new_gl_handler(tcl_interp, tk_display_win, direct, dbl_buff);
    if (!gl_handler)
	    RETURN_OBJ_ERROR("allocating Gl_handler object");

    obj = (Py_Gl_handler) PyObject_New(struct Py_Gl_handler, &Gl_handler_type);

    if (!obj)
    {
	delete_gl_handler(gl_handler);

	RETURN_OBJ_ERROR("allocating Py_Gl_handler object");
    }

    obj->gl_handler = gl_handler;

    return (PyObject *) obj;
}

static void delete_py_gl_handler(PyObject *self)
{
    Py_Gl_handler obj = (Py_Gl_handler) self;
    Gl_handler gl_handler = obj->gl_handler;

    delete_gl_handler(gl_handler);

    Py_TYPE(self)->tp_free(self);
}

/*
static int print_py_gl_handler(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_gl_handler(PyObject *self, PyObject *attr_name)
{
/*
    Gl_handler *obj = (Gl_handler *) self;
    Random_access *a = obj->py_handler;

    if (equal_strings(name, "par_file"))
	return Py_BuildValue("s", a->par_file);
    else if (equal_strings(name, "access_method"))
	return Py_BuildValue("s", access_method_name(a->access_method));
    else if (equal_strings(name, "data_format"))
	return get_Gl_handler_format(a);
    else
*/
	return PyObject_GenericGetAttr(self, attr_name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Gl_handler_sequence_methods =
{
    Gl_handler_length,
    Gl_handler_concat,
    Gl_handler_repeat,
    Gl_handler_item,
    Gl_handler_slice,
    Gl_handler_ass_item,
    Gl_handler_ass_slice
};

static PySequenceMethods Gl_handler_sequence_methods =
{
    Gl_handler_length,
    0,
    0,
    Gl_handler_item,
    0,
    Gl_handler_ass_item,
    0
};
*/

static PyTypeObject Gl_handler_type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    "GlHandler",                           /* tp_name */
    sizeof(struct Py_Gl_handler),           /* tp_basicsize */
    0,                                      /* tp_itemsize */
    (destructor) delete_py_gl_handler,      /* tp_dealloc */
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
    (getattrofunc) getattr_py_gl_handler,   /* tp_getattro */
    0,                                      /* tp_setattro */
    0,                                      /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                     /* tp_flags */
    "GlHandler -- GL window drawing handler", /* tp_doc */
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

static PyObject *init_Py_Gl_handler(PyObject *self, PyObject *args)
{
    PyObject *widget;
    int dbl_buff = 1;         /* DOUBLE_BUFFER always on unless specified off */
#ifdef USE_GL_TRUE
    int direct = GL_TRUE;
#else
    int direct = GL_FALSE;
#endif

    if (!PyArg_ParseTuple(args, "O|ii", &widget, &direct, &dbl_buff))
        RETURN_OBJ_ERROR("must have one argument: widget");

    return new_py_gl_handler(widget, direct ? CCPN_TRUE : CCPN_FALSE, 
                                   dbl_buff ? CCPN_TRUE : CCPN_FALSE );
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Gl_handler_type_methods[] =
{
    { "GlHandler",	(PyCFunction) init_Py_Gl_handler,	METH_VARARGS },
    { NULL,		NULL,			0 }
};

/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import GlHandler" in 
* a Python program. The function is usually called "initGl_handler": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

static struct PyModuleDef gl_handler_module_def =
{
    PyModuleDef_HEAD_INIT,
    "GlHandler",
    "CCPNMR GlHandler module (Python 3 compatible)",
    -1,
    Gl_handler_type_methods
};

PyMODINIT_FUNC PyInit_GlHandler(void)
{
    if (PyType_Ready(&Gl_handler_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&gl_handler_module_def);
    if (!m)
        return NULL;

    ErrorObject = PyErr_NewException("GlHandler.error", NULL, NULL);
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

#endif /* IGNORE_GL */

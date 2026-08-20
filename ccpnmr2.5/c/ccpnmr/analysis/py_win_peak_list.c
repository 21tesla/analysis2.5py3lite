/*
======================COPYRIGHT/LICENSE START==========================

py_win_peak_list.c: Part of the CcpNmr Analysis program

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
#include "py_win_peak_list.h"

#include "py_draw_handler.h"
#include "py_peak_list.h"
#include "python_util.h"

#define  NCOLORS  3

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Win_peak_list_type;

Bool is_py_win_peak_list(PyObject *obj)
{
    /*  below does not work because different *.so files end up
        with different addresses for Win_peak_list_type
        return (obj->ob_type == &Win_peak_list_type);
    */
    return valid_py_object(obj, &Win_peak_list_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *setIsSymbolDrawn(PyObject *self, PyObject *args)
{
    Bool isSymbolDrawn;
    Py_Win_peak_list obj = (Py_Win_peak_list) self;
    Win_peak_list win_peak_list = obj->win_peak_list;

    if (!PyArg_ParseTuple(args, "i", &isSymbolDrawn))
        RETURN_OBJ_ERROR("need one argument: isSymbolDrawn");

    is_symbol_drawn_win_peak_list(win_peak_list, isSymbolDrawn);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setIsTextDrawn(PyObject *self, PyObject *args)
{
    Bool isTextDrawn;
    Py_Win_peak_list obj = (Py_Win_peak_list) self;
    Win_peak_list win_peak_list = obj->win_peak_list;

    if (!PyArg_ParseTuple(args, "i", &isTextDrawn))
        RETURN_OBJ_ERROR("need one argument: isTextDrawn");

    is_text_drawn_win_peak_list(win_peak_list, isTextDrawn);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setIsTextPointerDrawn(PyObject *self, PyObject *args)
{
    Bool isTextPointerDrawn;
    Py_Win_peak_list obj = (Py_Win_peak_list) self;
    Win_peak_list win_peak_list = obj->win_peak_list;

    if (!PyArg_ParseTuple(args, "i", &isTextPointerDrawn))
        RETURN_OBJ_ERROR("need one argument: isTextPointerDrawn");

    is_text_pointer_drawn_win_peak_list(win_peak_list, isTextPointerDrawn);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *drawPeaks(PyObject *self, PyObject *args)
{
    int n, xdim, ydim, rdim, draw_method;
    int tile[MAX_NDIM];
    float xscale, yscale, first[MAX_NDIM+1], last[MAX_NDIM+1];
    float intensity_max, volume_max, xpix, ypix;
    float bg_color[NCOLORS], center[MAX_NDIM], thickness[MAX_NDIM];
    PyObject *first_obj, *last_obj, *handler_obj;
    PyObject *bg_color_obj, *tile_obj, *center_obj, *thickness_obj;
    Py_Win_peak_list obj = (Py_Win_peak_list) self;
    Win_peak_list win_peak_list = obj->win_peak_list;
    int ndim = win_peak_list->peak_list->ndim;
    Py_draw_handler py_draw_handler;
    Generic_ptr handler;
    Drawing_funcs *drawing_funcs;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "OiiffffOOiffOOOO", &handler_obj, &xdim, &ydim,
                          &xpix, &ypix, &xscale, &yscale,
                          &first_obj, &last_obj,
                          &draw_method, &intensity_max, &volume_max,
                          &bg_color_obj, &center_obj, &thickness_obj,
                          &tile_obj))
        RETURN_OBJ_ERROR("need arguments: handler, xdim, ydim, xpix, ypix, xscale, yscale, first, last, draw_method, intensity_max, volume_max, bg_color, center, thickness, tile");

    py_draw_handler = new_py_draw_handler(handler_obj);
    if (!py_draw_handler)
        RETURN_OBJ_ERROR("first argument must be handler object");

    handler = py_draw_handler->handler;
    drawing_funcs = py_draw_handler->drawing_funcs;
    delete_py_draw_handler(py_draw_handler);

    if ((xdim < 0) || (xdim >= ndim))
    {
        sprintf(error_msg, "xdim must be >= 0 and < %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (!win_peak_list->hasValueAxis)
    {
        if ((ydim < 0) || (ydim >= ndim))
        {
            sprintf(error_msg, "ydim must be >= 0 and < %d", ndim);
            RETURN_OBJ_ERROR(error_msg);
        }
    }

    rdim = ndim;
    if (win_peak_list->hasValueAxis)
        rdim++;

    if ((get_python_float_array(first_obj, MAX_NDIM+1, &n, first, error_msg) == CCPN_ERROR)
            || (n != rdim))
    {
        sprintf(error_msg, "first must be float list or tuple of size %d (is %d)", rdim, n);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_float_array(last_obj, MAX_NDIM+1, &n, last, error_msg) == CCPN_ERROR)
            || (n != rdim))
    {
        sprintf(error_msg, "last must be float list or tuple of size %d (is %d)", rdim, n);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_float_array(bg_color_obj, NCOLORS, &n, bg_color, error_msg) == CCPN_ERROR)
            || (n != NCOLORS))
    {
        sprintf(error_msg, "bg_color must be float list or tuple of size %d (is %d)", NCOLORS, n);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_float_array(center_obj, MAX_NDIM, &n, center, error_msg) == CCPN_ERROR)
            || (n != ndim))
    {
        sprintf(error_msg, "center must be float list or tuple of size %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_float_array(thickness_obj, MAX_NDIM, &n, thickness, error_msg) == CCPN_ERROR)
            || (n != ndim))
    {
        sprintf(error_msg, "thickness must be float list or tuple of size %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_int_array(tile_obj, MAX_NDIM, &n, tile, error_msg) == CCPN_ERROR)
            || (n != ndim))
    {
        sprintf(error_msg, "tile must be int list or tuple of size %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (draw_win_peak_list(win_peak_list, xdim, ydim, xpix, ypix,
                           xscale, yscale, first, last,
                           draw_method, intensity_max, volume_max,
                           bg_color, center, thickness, tile,
                           drawing_funcs, handler, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "setIsSymbolDrawn",	setIsSymbolDrawn,	METH_VARARGS },
    { "setIsTextDrawn",		setIsTextDrawn,		METH_VARARGS },
    { "setIsTextPointerDrawn",	setIsTextPointerDrawn,	METH_VARARGS },
    { "drawPeaks",		drawPeaks,		METH_VARARGS },
    { NULL,			NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_win_peak_list(Peak_list peak_list, Bool hasValueAxis)
{
    Win_peak_list win_peak_list;
    Py_Win_peak_list obj;

    win_peak_list = new_win_peak_list(peak_list, hasValueAxis);

    if (!win_peak_list)
        RETURN_OBJ_ERROR("allocating Win_peak_list object");

    PY_MALLOC(obj, struct Py_Win_peak_list, &Win_peak_list_type);

    if (!obj)
    {
        delete_win_peak_list(win_peak_list);

        RETURN_OBJ_ERROR("allocating Py_Win_peak_list object");
    }

    obj->win_peak_list = win_peak_list;

    return (PyObject *) obj;
}

static void delete_py_win_peak_list(PyObject *self)
{
    Py_Win_peak_list obj = (Py_Win_peak_list) self;
    Win_peak_list win_peak_list = obj->win_peak_list;

    /*
        printf("in delete_py_win_peak_list\n");
    */

    delete_win_peak_list(win_peak_list);

    PY_FREE(self);
}

/*
static int print_py_win_peak_list(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_win_peak_list(PyObject *self, char *name)
{
    /*
        Win_peak_list *obj = (Win_peak_list *) self;
        Random_access *a = obj->py_handler;

        if (equal_strings(name, "par_file"))
    	return Py_BuildValue("s", a->par_file);
        else if (equal_strings(name, "access_method"))
    	return Py_BuildValue("s", access_method_name(a->access_method));
        else if (equal_strings(name, "data_format"))
    	return get_Win_peak_list_format(a);
        else
    */
    return Py_FindMethod(py_handler_methods, self, name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Win_peak_list_sequence_methods =
{
    Win_peak_list_length,
    Win_peak_list_concat,
    Win_peak_list_repeat,
    Win_peak_list_item,
    Win_peak_list_slice,
    Win_peak_list_ass_item,
    Win_peak_list_ass_slice
};

static PySequenceMethods Win_peak_list_sequence_methods =
{
    Win_peak_list_length,
    0,
    0,
    Win_peak_list_item,
    0,
    Win_peak_list_ass_item,
    0
};
*/

static PyTypeObject Win_peak_list_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "WinPeakList", /* name */
    sizeof(struct Py_Win_peak_list), /* basicsize */
    0, /* itemsize */
    delete_py_win_peak_list, /* destructor */
    0, /* printfunc */
    getattr_py_win_peak_list, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Win_peak_list_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Win_peak_list(PyObject *self, PyObject *args)
{
    Bool hasValueAxis;
    PyObject *obj, *peak_list_obj;

    if (!PyArg_ParseTuple(args, "Oi", &peak_list_obj, &hasValueAxis))
        RETURN_OBJ_ERROR("must have two arguments: peak list and hasValueAxis");

    if (!is_py_peak_list(peak_list_obj))
        RETURN_OBJ_ERROR("first argument must be peak list");

    obj = new_py_win_peak_list(((Py_Peak_list) peak_list_obj)->peak_list,
                               hasValueAxis);

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


static struct PyMethodDef Win_peak_list_type_methods[] =
{
    { "WinPeakList",	(PyCFunction) init_Py_Win_peak_list,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import WinPeakList" in
* a Python program. The function is usually called "initWin_peak_list": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initWinPeakList(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Win_peak_list_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("WinPeakList", Win_peak_list_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("WinPeakList.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);

    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module WinPeakList");
}

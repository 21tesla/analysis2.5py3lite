/*
======================COPYRIGHT/LICENSE START==========================

py_peak.c: Part of the CcpNmr Analysis program

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
#include "py_peak.h"

#include "py_block_file.h"
#include "python_util.h"

#define RETURN_CHECK_OBJ_ERROR(message) \
	{  if (!ErrorObject) \
           {  ErrorObject = PyErr_NewException("Peak.error", NULL, NULL); \
              Py_INCREF(ErrorObject); } \
           RETURN_OBJ_ERROR(message); }

#define  NCOLORS  3

static PyObject *ErrorObject = NULL;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Peak_type;

Bool is_py_peak(PyObject *obj)
{
    /*  below does not work because different *.so files end up
        with different addresses for Peak_type
        return (obj->ob_type == &Peak_type);
    */
    return valid_py_object(obj, &Peak_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *getIsSelected(PyObject *self, PyObject *args)
{
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;

    if (!PyArg_ParseTuple(args, ""))
        RETURN_CHECK_OBJ_ERROR("should have no arguments");

    return Py_BuildValue("i", get_is_selected_peak(peak));
}

static PyObject *setIsSelected(PyObject *self, PyObject *args)
{
    Bool isSelected;
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;

    if (!PyArg_ParseTuple(args, "i", &isSelected))
        RETURN_CHECK_OBJ_ERROR("need one argument: isSelected");

    set_is_selected_peak(peak, isSelected);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setText(PyObject *self, PyObject *args)
{
    char *text;
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "s", &text))
        RETURN_CHECK_OBJ_ERROR("need one argument: text");

    if (set_text_peak(peak, text, error_msg) == CCPN_ERROR)
        RETURN_CHECK_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setTextOffset(PyObject *self, PyObject *args)
{
    int dim;
    float offset;
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "if", &dim, &offset))
        RETURN_CHECK_OBJ_ERROR("need two arguments: dim, offset");

    set_text_offset_peak(peak, dim, offset);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *resetTextOffset(PyObject *self, PyObject *args)
{
    int dim;
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "i", &dim))
        RETURN_CHECK_OBJ_ERROR("need one argument: dim");

    reset_text_offset_peak(peak, dim);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setPosition(PyObject *self, PyObject *args)
{
    int n;
    float position[MAX_NDIM];
    PyObject *position_obj;
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &position_obj))
        RETURN_CHECK_OBJ_ERROR("need one argument: position");

    if ((get_python_float_array(position_obj, MAX_NDIM, &n, position,
                                error_msg) == CCPN_ERROR)
            || (n != peak->ndim))
    {
        sprintf(error_msg, "position must be list or tuple of size %d",
                peak->ndim);
        RETURN_CHECK_OBJ_ERROR(error_msg);
    }

    set_position_peak(peak, position);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setNumAliasing(PyObject *self, PyObject *args)
{
    int n, num_aliasing[MAX_NDIM];
    PyObject *num_aliasing_obj;
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &num_aliasing_obj))
        RETURN_CHECK_OBJ_ERROR("need one argument: num_aliasing");

    if ((get_python_int_array(num_aliasing_obj, MAX_NDIM, &n, num_aliasing,
                              error_msg) == CCPN_ERROR)
            || (n != peak->ndim))
    {
        sprintf(error_msg, "num_aliasing must be int list or tuple of size %d",
                peak->ndim);
        RETURN_CHECK_OBJ_ERROR(error_msg);
    }

    set_num_aliasing_peak(peak, num_aliasing);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *getIntensity(PyObject *self, PyObject *args)
{
    float intensity;
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;
    PyObject *block_file_obj;
    Block_file block_file;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &block_file_obj))
        RETURN_CHECK_OBJ_ERROR("need argument: block_file");

    if (!is_py_some_block_file(block_file_obj))
        RETURN_CHECK_OBJ_ERROR("must have BlockFile object");

    if (is_py_block_file(block_file_obj))
        block_file = ((Py_Block_file) block_file_obj)->block_file;
    else
        block_file = (Block_file) (((Py_Shape_block_file) block_file_obj)->shape_block_file);

    if (get_intensity_peak(peak, &intensity, block_file, error_msg) == CCPN_ERROR)
        RETURN_CHECK_OBJ_ERROR(error_msg);

    return Py_BuildValue("f", intensity);
}

static PyObject *setIntensity(PyObject *self, PyObject *args)
{
    float intensity;
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;

    if (!PyArg_ParseTuple(args, "f", &intensity))
        RETURN_CHECK_OBJ_ERROR("need argument: intensity");

    set_intensity_peak(peak, intensity);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setLineWidth(PyObject *self, PyObject *args)
{
    int dim;
    float line_width;
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "if", &dim, &line_width))
        RETURN_CHECK_OBJ_ERROR("need two arguments: dim, line_width");

    set_line_width_peak(peak, dim, line_width);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *fitVolume(PyObject *self, PyObject *args)
{
    int i, n, method;
    Bool valid;
    float volume;
    int dd[MAX_NDIM];
    Bool dim_done[MAX_NDIM];
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;
    PyObject *block_file_obj, *dim_done_obj = NULL;
    Block_file block_file;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "iO|O", &method, &block_file_obj, &dim_done_obj))
        RETURN_CHECK_OBJ_ERROR("need arguments: method block_file [ dim_done ]");

    if (!is_py_some_block_file(block_file_obj))
        RETURN_CHECK_OBJ_ERROR("must have BlockFile object");

    if (is_py_block_file(block_file_obj))
        block_file = ((Py_Block_file) block_file_obj)->block_file;
    else
        block_file = (Block_file) (((Py_Shape_block_file) block_file_obj)->shape_block_file);

    if (dim_done_obj)
    {
        if ((get_python_int_array(dim_done_obj, MAX_NDIM, &n,
                                  dd, error_msg) == CCPN_ERROR) || (n != peak->ndim))
        {
            sprintf(error_msg, "dim_done must be int list or tuple of size %d", peak->ndim);
            RETURN_OBJ_ERROR(error_msg);
        }

        for (i = 0; i < peak->ndim; i++)
            dim_done[i] = dd[i];
    }
    else
    {
        for (i = 0; i < peak->ndim; i++)
            dim_done[i] = CCPN_TRUE;
    }

    if (fit_volume_peak(peak, &volume, &valid, method, block_file, dim_done, error_msg) == CCPN_ERROR)
        RETURN_CHECK_OBJ_ERROR(error_msg);

    if (!valid)
    {
        Py_INCREF(Py_None);
        return Py_None;
    }

    return Py_BuildValue("f", volume);
}

static PyObject *setVolume(PyObject *self, PyObject *args)
{
    float volume;
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;

    if (!PyArg_ParseTuple(args, "f", &volume))
        RETURN_CHECK_OBJ_ERROR("need argument: volume");

    set_volume_peak(peak, volume);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *fitCenter(PyObject *self, PyObject *args)
{
    int i, n, method;
    Bool valid;
    int dd[MAX_NDIM];
    Bool dim_done[MAX_NDIM];
    float center[MAX_NDIM];
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;
    PyObject *block_file_obj, *dim_done_obj = NULL;
    Block_file block_file;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "iO|O", &method, &block_file_obj, &dim_done_obj))
        RETURN_CHECK_OBJ_ERROR("need arguments: method block_file [ dim_done ]");

    if (!is_py_some_block_file(block_file_obj))
        RETURN_CHECK_OBJ_ERROR("must have BlockFile object");

    if (is_py_block_file(block_file_obj))
        block_file = ((Py_Block_file) block_file_obj)->block_file;
    else
        block_file = (Block_file) (((Py_Shape_block_file) block_file_obj)->shape_block_file);

    if (dim_done_obj)
    {
        if ((get_python_int_array(dim_done_obj, MAX_NDIM, &n,
                                  dd, error_msg) == CCPN_ERROR) || (n != peak->ndim))
        {
            sprintf(error_msg, "dim_done must be int list or tuple of size %d", peak->ndim);
            RETURN_OBJ_ERROR(error_msg);
        }

        for (i = 0; i < peak->ndim; i++)
            dim_done[i] = dd[i];
    }
    else
    {
        for (i = 0; i < peak->ndim; i++)
            dim_done[i] = CCPN_TRUE;
    }

    if (fit_center_peak(peak, center, &valid, method, block_file, dim_done, error_msg) == CCPN_ERROR)
        RETURN_CHECK_OBJ_ERROR(error_msg);

    if (!valid)
    {
        Py_INCREF(Py_None);
        return Py_None;
    }

    return get_python_float_list(peak->ndim, center);
}

/* this fits linewidth in points */
static PyObject *fitLinewidth(PyObject *self, PyObject *args)
{
    int i, n;
    int dd[MAX_NDIM];
    Bool dim_done[MAX_NDIM];
    float linewidth[MAX_NDIM];
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;
    PyObject *block_file_obj, *dim_done_obj = NULL;
    Block_file block_file;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O|O", &block_file_obj, &dim_done_obj))
        RETURN_CHECK_OBJ_ERROR("need arguments: block_file [ dim_done ]");

    if (!is_py_some_block_file(block_file_obj))
        RETURN_CHECK_OBJ_ERROR("must have BlockFile object");

    if (is_py_block_file(block_file_obj))
        block_file = ((Py_Block_file) block_file_obj)->block_file;
    else
        block_file = (Block_file) (((Py_Shape_block_file) block_file_obj)->shape_block_file);

    if (dim_done_obj)
    {
        if ((get_python_int_array(dim_done_obj, MAX_NDIM, &n,
                                  dd, error_msg) == CCPN_ERROR) || (n != peak->ndim))
        {
            sprintf(error_msg, "dim_done must be int list or tuple of size %d", peak->ndim);
            RETURN_OBJ_ERROR(error_msg);
        }

        for (i = 0; i < peak->ndim; i++)
            dim_done[i] = dd[i];
    }
    else
    {
        for (i = 0; i < peak->ndim; i++)
            dim_done[i] = CCPN_TRUE;
    }

    if (fit_linewidth_peak(peak, linewidth, block_file, dim_done, error_msg) == CCPN_ERROR)
        RETURN_CHECK_OBJ_ERROR(error_msg);

    return get_python_float_list(peak->ndim, linewidth);
}

/*
static PyObject *setQuality(PyObject *self, PyObject *args)
{
    float quality;
    Py_Peak obj = (Py_Peak) self;
    Peak peak = obj->peak;

    if (!PyArg_ParseTuple(args, "f", &quality))
        RETURN_CHECK_OBJ_ERROR("need argument: quality");

    set_quality_peak(peak, quality);

    Py_INCREF(Py_None);
    return Py_None;
}
*/

static struct PyMethodDef py_handler_methods[] =
{
    { "getIsSelected",		getIsSelected,		METH_VARARGS },
    { "setIsSelected",		setIsSelected,		METH_VARARGS },
    { "setText",		setText,		METH_VARARGS },
    { "setTextOffset",		setTextOffset,		METH_VARARGS },
    { "resetTextOffset",	resetTextOffset,	METH_VARARGS },
    { "setPosition",		setPosition,		METH_VARARGS },
    { "setNumAliasing",		setNumAliasing,		METH_VARARGS },
    { "getIntensity",		getIntensity,		METH_VARARGS },
    { "setIntensity",		setIntensity,		METH_VARARGS },
    { "setLineWidth",		setLineWidth,		METH_VARARGS },
    { "fitVolume",		fitVolume,		METH_VARARGS },
    { "setVolume",		setVolume,		METH_VARARGS },
    { "fitCenter",		fitCenter,		METH_VARARGS },
    { "fitLinewidth",		fitLinewidth,		METH_VARARGS },
    /*
        { "setQuality",		setQuality,		METH_VARARGS },
    */
    { NULL,			NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

PyObject *new_py_peak(Peak peak)
{
    Py_Peak obj;

    if (!peak)
        RETURN_CHECK_OBJ_ERROR("allocating Peak object");

#ifdef WIN64
    Peak_type.ob_type = &PyType_Type;
#endif

    PY_MALLOC(obj, struct Py_Peak, &Peak_type);

    if (!obj)
    {
        delete_peak(peak);

        RETURN_CHECK_OBJ_ERROR("allocating Py_Peak object");
    }

    obj->peak = peak;

    return (PyObject *) obj;
}

void delete_py_peak(PyObject *self)
{
    /*
        Py_Peak obj = (Py_Peak) self;
        Peak peak = obj->peak;
    */
    /*
        printf("in delete_py_peak\n");
    */
    /*  delete_peak done in remove_peak_peak_list() in peak_list.c, so do not do here  */
    /*
        delete_peak(peak);
    */

    PY_FREE(self);
}

/*
static int print_py_peak(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_peak(PyObject *self, char *name)
{
    /*
        Peak *obj = (Peak *) self;
        Random_access *a = obj->py_handler;

        if (equal_strings(name, "par_file"))
    	return Py_BuildValue("s", a->par_file);
        else if (equal_strings(name, "access_method"))
    	return Py_BuildValue("s", access_method_name(a->access_method));
        else if (equal_strings(name, "data_format"))
    	return get_Peak_format(a);
        else
    */
    return Py_FindMethod(py_handler_methods, self, name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Peak_sequence_methods =
{
    Peak_length,
    Peak_concat,
    Peak_repeat,
    Peak_item,
    Peak_slice,
    Peak_ass_item,
    Peak_ass_slice
};

static PySequenceMethods Peak_sequence_methods =
{
    Peak_length,
    0,
    0,
    Peak_item,
    0,
    Peak_ass_item,
    0
};
*/

static PyTypeObject Peak_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "Peak", /* name */
    sizeof(struct Py_Peak), /* basicsize */
    0, /* itemsize */
    delete_py_peak, /* destructor */
    0, /* printfunc */
    getattr_py_peak, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Peak_sequence_methods*/ /* PySequenceMethods */
};


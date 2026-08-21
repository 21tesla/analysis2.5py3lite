/*
======================COPYRIGHT/LICENSE START==========================

py_peak_cluster.c: Part of the CcpNmr Analysis program

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
#include "py_peak_cluster.h"

#include "py_draw_handler.h"
#include "py_peak.h"
#include "python_util.h"

static PyObject *ErrorObject = NULL;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Peak_cluster_type;

Bool is_py_peak_cluster(PyObject *obj)
{
    /*  below does not work because different *.so files end up
        with different addresses for Peak_cluster_type
        return (obj->ob_type == &Peak_cluster_type);
    */
    return valid_py_object(obj, &Peak_cluster_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *addPeak(PyObject *self, PyObject *args)
{
    Py_Peak_cluster obj = (Py_Peak_cluster) self;
    Peak_cluster peak_cluster = obj->peak_cluster;
    PyObject *peak_obj;
    Peak peak;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &peak_obj))
        RETURN_OBJ_ERROR("need argument: peak");

    if (!is_py_peak(peak_obj))
        RETURN_OBJ_ERROR("argument must be a peak");

    peak = ((Py_Peak) peak_obj)->peak;
    if (add_peak_peak_cluster(peak_cluster, peak, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *removePeak(PyObject *self, PyObject *args)
{
    Py_Peak_cluster obj = (Py_Peak_cluster) self;
    Peak_cluster peak_cluster = obj->peak_cluster;
    PyObject *peak_obj;
    Peak peak;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &peak_obj))
        RETURN_OBJ_ERROR("need argument: peak");

    if (!is_py_peak(peak_obj))
        RETURN_OBJ_ERROR("argument must be a peak");

    peak = ((Py_Peak) peak_obj)->peak;
    remove_peak_peak_cluster(peak_cluster, peak);

    Py_INCREF(Py_None);
    return Py_None;
}

static Peak *alloc_peaks(int npeaks)
{
    Peak *peaks;

    MALLOC_NEW(peaks, Peak, npeaks);

    return peaks;
}

static PyObject *setPeaks(PyObject *self, PyObject *args)
{
    Py_Peak_cluster obj = (Py_Peak_cluster) self;
    Peak_cluster peak_cluster = obj->peak_cluster;
    int i, npeaks;
    PyObject *peak_list_obj, *peak_obj;
    Peak *peaks;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &peak_list_obj))
        RETURN_OBJ_ERROR("need argument: list_of_peaks");

    npeaks = get_python_list_size(peak_list_obj);
    if (npeaks < 0)
        RETURN_OBJ_ERROR("argument should be a list or tuple of peaks");

    peaks = alloc_peaks(npeaks);
    if (!peaks)
        RETURN_OBJ_ERROR("allocating peaks memory");

    for (i = 0; i < npeaks; i++)
    {
        peak_obj = get_python_object_by_index(peak_list_obj, i);
        if (!is_py_peak(peak_obj))
        {
            FREE(peaks, Peak);
            sprintf(error_msg, "element %d of list/tuple is not a peak", i);
            RETURN_OBJ_ERROR(error_msg);
        }

        peaks[i] = ((Py_Peak) peak_obj)->peak;
    }

    set_peaks_peak_cluster(peak_cluster, npeaks, npeaks, peaks);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *clearPeaks(PyObject *self, PyObject *args)
{
    Py_Peak_cluster obj = (Py_Peak_cluster) self;
    Peak_cluster peak_cluster = obj->peak_cluster;
    PyObject *peak_obj;
    Peak peak;
    Line error_msg;

    if (!PyArg_ParseTuple(args, ""))
        RETURN_OBJ_ERROR("no arguments should be given");

    clear_peaks_peak_cluster(peak_cluster);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setText(PyObject *self, PyObject *args)
{
    char *text;
    Py_Peak_cluster obj = (Py_Peak_cluster) self;
    Peak_cluster peak_cluster = obj->peak_cluster;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "s", &text))
        RETURN_OBJ_ERROR("need argument: text");

    if (set_text_peak_cluster(peak_cluster, text, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setDimText(PyObject *self, PyObject *args)
{
    int dim;
    char *text;
    Py_Peak_cluster obj = (Py_Peak_cluster) self;
    Peak_cluster peak_cluster = obj->peak_cluster;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "is", &dim, &text))
        RETURN_OBJ_ERROR("need arguments: dim, text");

    if (set_dim_text_peak_cluster(peak_cluster, dim, text, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *draw(PyObject *self, PyObject *args)
{
    Py_Peak_cluster obj = (Py_Peak_cluster) self;
    Peak_cluster peak_cluster = obj->peak_cluster;
    int xdim, ydim, ndim = peak_cluster->ndim;
    PyObject *handler_obj;
    Py_draw_handler py_draw_handler;
    Generic_ptr handler;
    Drawing_funcs *drawing_funcs;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "Oii", &handler_obj, &xdim, &ydim))
        RETURN_OBJ_ERROR("need arguments: handler, xdim, ydim");

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

    if ((ydim < 0) || (ydim >= ndim))
    {
        sprintf(error_msg, "ydim must be >= 0 and < %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    draw_peak_cluster(peak_cluster, xdim, ydim, drawing_funcs, handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "addPeak",		addPeak,		METH_VARARGS },
    { "removePeak",		removePeak,		METH_VARARGS },
    { "clearPeaks",		clearPeaks,		METH_VARARGS },
    { "setPeaks",		setPeaks,		METH_VARARGS },
    { "setText",		setText,		METH_VARARGS },
    { "setDimText",		setDimText,		METH_VARARGS },
    { "draw",			draw,			METH_VARARGS },
    { NULL,			NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

PyObject *new_py_peak_cluster(int ndim, int type)
{
    Peak_cluster peak_cluster;
    Py_Peak_cluster obj;
    PeakClusterType cluster_type;

    if ((type < 0) || (type > PEAK_CLUSTER_SYMMETRY))
        RETURN_OBJ_ERROR("cluster type not in correct range 0 to 3");

    cluster_type = type;
    peak_cluster = new_peak_cluster(ndim, cluster_type);

    if (!peak_cluster)
        RETURN_OBJ_ERROR("allocating Peak_cluster object");

#ifdef WIN64
    Peak_cluster_type.ob_type = &PyType_Type;
#endif

    obj = (Py_Peak_cluster) PyObject_New(struct Py_Peak_cluster, &Peak_cluster_type);

    if (!obj)
    {
        delete_peak_cluster(peak_cluster);

        RETURN_OBJ_ERROR("allocating Py_Peak_cluster object");
    }

    obj->peak_cluster = peak_cluster;

    return (PyObject *) obj;
}

void delete_py_peak_cluster(PyObject *self)
{
    Py_Peak_cluster obj = (Py_Peak_cluster) self;
    Peak_cluster peak_cluster = obj->peak_cluster;
    /*
        printf("in delete_py_peak_cluster\n");
    */
    delete_peak_cluster(peak_cluster);

    Py_TYPE(self)->tp_free(self);
}

/*
static int print_py_peak_cluster(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_peak_cluster(PyObject *self, PyObject *attr_name)
{
    /*
        Peak_cluster *obj = (Peak_cluster *) self;
        Random_access *a = obj->py_handler;

        if (equal_strings(name, "par_file"))
    	return Py_BuildValue("s", a->par_file);
        else if (equal_strings(name, "access_method"))
    	return Py_BuildValue("s", access_method_name(a->access_method));
        else if (equal_strings(name, "data_format"))
    	return get_Peak_format(a);
        else
    */
    return PyObject_GenericGetAttr(self, attr_name);
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

static PyTypeObject Peak_cluster_type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    "Peak",                                      /* tp_name */
    sizeof(struct Py_Peak_cluster),              /* tp_basicsize */
    0,                                           /* tp_itemsize */
    (destructor) delete_py_peak_cluster,         /* tp_dealloc */
    0,                                           /* tp_vectorcall */
    0,                                           /* tp_getattr */
    0,                                           /* tp_setattr */
    0,                                           /* tp_as_async */
    0,                                           /* tp_repr */
    0,                                           /* tp_as_number */
    0,                                           /* tp_as_sequence */
    0,                                           /* tp_as_mapping */
    0,                                           /* tp_hash */
    0,                                           /* tp_call */
    0,                                           /* tp_str */
    (getattrofunc) getattr_py_peak_cluster,      /* tp_getattro */
    0,                                           /* tp_setattro */
    0,                                           /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                          /* tp_flags */
    "Peak -- peak cluster",                      /* tp_doc */
    0,                                           /* tp_traverse */
    0,                                           /* tp_clear */
    0,                                           /* tp_richcompare */
    0,                                           /* tp_weaklistoffset */
    0,                                           /* tp_iter */
    0,                                           /* tp_iternext */
    py_handler_methods,                          /* tp_methods */
    0,                                           /* tp_members */
    0,                                           /* tp_getset */
    0,                                           /* tp_base */
    0,                                           /* tp_dict */
    0,                                           /* tp_descr_get */
    0,                                           /* tp_descr_set */
    0,                                           /* tp_dictoffset */
    0,                                           /* tp_init */
    0,                                           /* tp_alloc */
    0,                                           /* tp_new */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Peak_cluster(PyObject *self, PyObject *args)
{
    int ndim, type;
    PeakClusterType cluster_type;
    PyObject *obj;

    if (!PyArg_ParseTuple(args, "ii", &ndim, &type))
        RETURN_OBJ_ERROR("must have arguments: ndim, cluster_type");

    obj = new_py_peak_cluster(ndim, type);

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

static struct PyMethodDef Peak_cluster_type_methods[] =
{
    { "PeakCluster",	(PyCFunction) init_Py_Peak_cluster,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import PeakCluster" in
* a Python program. The function is usually called "initWin_peak_list": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

static struct PyModuleDef peak_cluster_module_def =
{
    PyModuleDef_HEAD_INIT,
    "PeakCluster",
    "CCPNMR PeakCluster module (Python 3 compatible)",
    -1,
    Peak_cluster_type_methods
};

PyMODINIT_FUNC PyInit_PeakCluster(void)
{
    if (PyType_Ready(&Peak_cluster_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&peak_cluster_module_def);
    if (!m)
        return NULL;

    ErrorObject = PyErr_NewException("PeakCluster.error", NULL, NULL);
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


/*
======================COPYRIGHT/LICENSE START==========================

py_cloud_util.c: Part of the CcpNmr Analysis program

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
#include "python_util.h"

static PyObject *ErrorObject;   /* locally-raised exception */

static double py_distance2(PyObject *obj_coords1, PyObject *obj_coords2)
{
    int i, size;
    double c1, c2, dc, d2 = 0;

    size = PyList_GET_SIZE(obj_coords1);
    if (size != PyList_GET_SIZE(obj_coords2))
	return -1.0;

    for (i = 0; i < size; i++)
    {
	c1 = PyFloat_AS_DOUBLE(PyList_GET_ITEM(obj_coords1, i));
	c2 = PyFloat_AS_DOUBLE(PyList_GET_ITEM(obj_coords2, i));
	dc = c2 - c1;
	d2 += dc * dc;
    }

    return d2;
}

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *py_distance_list(PyObject *self, PyObject *args)
{
    int i, size;
    double d2;
    Line error_msg;
    PyObject *obj_coords_list1, *obj_coords_list2, *dist_list;

    if (!PyArg_ParseTuple(args, "OO", &obj_coords_list1, &obj_coords_list2))
        RETURN_OBJ_ERROR("must have two arguments: objCoords1, objCoords2");

    if (!PyList_Check(obj_coords_list1))
        RETURN_OBJ_ERROR("first argument must be list");

    if (!PyList_Check(obj_coords_list2))
        RETURN_OBJ_ERROR("second argument must be list");

    size = PyList_GET_SIZE(obj_coords_list1);
    if (size != PyList_GET_SIZE(obj_coords_list2))
        RETURN_OBJ_ERROR("two arguments must be lists of equal length");

    dist_list = PyList_New(size);
    for (i = 0; i < size; i++)
    {
	d2 = py_distance2(PyList_GET_ITEM(obj_coords_list1, i),
				PyList_GET_ITEM(obj_coords_list2, i));
	if (d2 < 0)
	{
	    sprintf(error_msg, "size mismatch for entry %d of two lists", i);
	    RETURN_OBJ_ERROR(error_msg);
	}

	PyList_SET_ITEM(dist_list, i, Py_BuildValue("d", sqrt(d2)));
    }

    return dist_list;
}

static PyObject *py_distance_stats(PyObject *self, PyObject *args)
{
    int i, size;
    double d, d2, avg, rms;
    Line error_msg;
    PyObject *obj_coords_list1, *obj_coords_list2, *dist_list, *stats_list;

    if (!PyArg_ParseTuple(args, "OO", &obj_coords_list1, &obj_coords_list2))
        RETURN_OBJ_ERROR("must have two arguments: objCoords1, objCoords2");

    if (!PyList_Check(obj_coords_list1))
        RETURN_OBJ_ERROR("first argument must be list");

    if (!PyList_Check(obj_coords_list2))
        RETURN_OBJ_ERROR("second argument must be list");

    size = PyList_GET_SIZE(obj_coords_list1);
    if (size != PyList_GET_SIZE(obj_coords_list2))
        RETURN_OBJ_ERROR("two arguments must be lists of equal length");

    stats_list = PyTuple_New(3);
    dist_list = PyList_New(size);
    avg = rms = 0;
    for (i = 0; i < size; i++)
    {
	d2 = py_distance2(PyList_GET_ITEM(obj_coords_list1, i),
				PyList_GET_ITEM(obj_coords_list2, i));
	if (d2 < 0)
	{
	    sprintf(error_msg, "size mismatch for entry %d of two lists", i);
	    RETURN_OBJ_ERROR(error_msg);
	}

	rms += d2;
	d = sqrt(d2);
	avg += d;
	PyList_SET_ITEM(dist_list, i, Py_BuildValue("d", d));
    }

    avg /= size;
    rms /= size;
    rms = sqrt(MAX(0, rms-avg*avg));
    PyTuple_SET_ITEM(stats_list, 0, dist_list);
    PyTuple_SET_ITEM(stats_list, 1, Py_BuildValue("d", avg));
    PyTuple_SET_ITEM(stats_list, 2, Py_BuildValue("d", rms));

    return stats_list;
}

static PyObject *py_distance_prob(PyObject *self, PyObject *args)
{
    int i, n, dist_size, prob_size;
    double d, s, grid;
    PyObject *dist_list, *prob_list;

    if (!PyArg_ParseTuple(args, "OOd", &dist_list, &prob_list, &grid))
        RETURN_OBJ_ERROR("must have three arguments: distanceList, probDistrib, grid");

    if (!PyList_Check(dist_list))
        RETURN_OBJ_ERROR("first argument must be list");

    if (!PyList_Check(prob_list))
        RETURN_OBJ_ERROR("second argument must be list");

    dist_size = PyList_GET_SIZE(dist_list);
    prob_size = PyList_GET_SIZE(prob_list);

    s = 0;
    for (i = 0; i < dist_size; i++)
    {
	d = PyFloat_AsDouble(PyList_GET_ITEM(dist_list, i));
	d /= grid;
	n = FLOOR(d);
	n = MAX(n, 0);
	n = MIN(n, prob_size-1);
	s += PyFloat_AsDouble(PyList_GET_ITEM(prob_list, n));
    }

    s /= dist_size;

    return Py_BuildValue("d", s);
}

static PyObject *py_reqd_nearness(PyObject *self, PyObject *args)
{
    int i, size, reqd, result = 0;
    double d2, dmin, dmax, d2min, d2max;
    Line error_msg;
    PyObject *obj_coords_list1, *obj_coords_list2;

    if (!PyArg_ParseTuple(args, "OOddi", &obj_coords_list1, &obj_coords_list2,
		&dmin, &dmax, &reqd))
        RETURN_OBJ_ERROR("must have five arguments: objCoords1, objCoords2, dmin, dmax, reqd");

    if (!PyList_Check(obj_coords_list1))
        RETURN_OBJ_ERROR("first argument must be list");

    if (!PyList_Check(obj_coords_list2))
        RETURN_OBJ_ERROR("second argument must be list");

    size = PyList_GET_SIZE(obj_coords_list1);
    if (size != PyList_GET_SIZE(obj_coords_list2))
        RETURN_OBJ_ERROR("two arguments must be lists of equal length");

    if (reqd <= 0)
    {
	result = 1;
    }
    else
    {
	d2min = dmin * dmin;
	d2max = dmax * dmax;
	for (i = 0; i < size; i++)
	{
	    if (reqd > (size-i))
	    break;

	    d2 = py_distance2(PyList_GET_ITEM(obj_coords_list1, i),
				PyList_GET_ITEM(obj_coords_list2, i));
	    if (d2 < 0)
	    {
		sprintf(error_msg, "size mismatch for entry %d of two lists", i);
		RETURN_OBJ_ERROR(error_msg);
	    }

	    if ((d2 > d2min) && (d2 < d2max))
	    {
		reqd--;
		if (reqd == 0)
		{
		    result = 1;
		    break;
		}
	    }
	}
    }

    return Py_BuildValue("i", result);
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Cloud_util_type_methods[] =
{
    { "distanceList",	(PyCFunction) py_distance_list,	METH_VARARGS },
    { "distanceStats",	(PyCFunction) py_distance_stats,METH_VARARGS },
    { "distanceProb",	(PyCFunction) py_distance_prob,	METH_VARARGS },
    { "reqdNearness",	(PyCFunction) py_reqd_nearness,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import CloudUtil" in 
* a Python program. The function is usually called "initCloud_util": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initCloudUtil(void)
{
    PyObject *m, *d;

    /* create the module and add the functions */
    m = Py_InitModule("CloudUtil", Cloud_util_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("CloudUtil.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module CloudUtil");
}

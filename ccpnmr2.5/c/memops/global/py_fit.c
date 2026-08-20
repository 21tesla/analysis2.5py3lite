
/*
======================COPYRIGHT/LICENSE START==========================

py_fit.c: Part of the CcpNmr Analysis program

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
#include "py_fit.h"

#include "fit.h"

#include "python_util.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Fit_method_type;

Bool is_py_fit(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Fit_method_type
    return (obj->ob_type == &Fit_method_type);
*/
    return valid_py_object(obj, &Fit_method_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static CcpnStatus alloc_fit_memory(int nparams, int n, float **params_fit,
			float **params_dev, float **y_fit)
{
    MALLOC(*params_fit, float, nparams);
    MALLOC(*params_dev, float, nparams);
    MALLOC(*y_fit, float, n);

    return CCPN_OK;
}

static PyObject *run(PyObject *self, PyObject *args)
{
    Py_Fit_method py_fit = (Py_Fit_method) self;
    Fit_method fit = py_fit->fit;
    int method = fit->method;
    int nparams = get_method_nparams(method);
    int n1, n2;
    float *x, *y, *params_fit, *params_dev, *y_fit, chisq;
    PyObject *x_obj, *y_obj, *params_obj, *params_dev_obj, *y_fit_obj, *result;
    Line error_msg;
    CcpnStatus status;
 
    if (!PyArg_ParseTuple(args, "OO", &x_obj, &y_obj))
        RETURN_OBJ_ERROR("need two arguments: x, y,");

    if (get_python_float_alloc_array(x_obj, &n1, &x, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    if (get_python_float_alloc_array(y_obj, &n2, &y, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    if (n1 != n2)
    {
        sprintf(error_msg, "x any y must both be of size %d", n1);
	FREE(x, float);
	FREE(y, float);
        RETURN_OBJ_ERROR(error_msg);
    }
 
    if (alloc_fit_memory(nparams, n1, &params_fit, &params_dev, &y_fit) == CCPN_ERROR)
    {
	FREE(x, float);
	FREE(y, float);
        RETURN_OBJ_ERROR("allocating fit memory");
    }

    status = run_fit(fit, n1, x, y, params_fit, params_dev, y_fit, &chisq, error_msg);

    if (status == CCPN_OK)
    {
	params_obj = get_python_float_list(nparams, params_fit);
	params_dev_obj = get_python_float_list(nparams, params_dev);
        y_fit_obj  = get_python_float_list(n1, y_fit);
    }

    FREE(x, float);
    FREE(y, float);
    FREE(params_fit, float);
    FREE(params_dev, float);
    FREE(y_fit, float);

    if (status == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    result = Py_BuildValue("(OOOf)", params_obj, params_dev_obj, y_fit_obj, chisq);
    Py_DECREF(params_obj);
    Py_DECREF(params_dev_obj);
    Py_DECREF(y_fit_obj);

    return result;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "run",		run,		METH_VARARGS },
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static Py_Fit_method new_py_fit(int method, float noise)
{
    Py_Fit_method py_fit;
    Fit_method fit;
    Line error_msg;

    if ((method < 0) || (method >= NFIT_METHODS))
    {
	sprintf(error_msg, "method must be between 0 and %d", NFIT_METHODS-1);
	RETURN_OBJ_ERROR(error_msg);
    }

    fit = new_fit(method, noise);

    if (!fit)
	 RETURN_OBJ_ERROR("allocating Fit_method object");

    PY_MALLOC(py_fit, struct Py_Fit_method, &Fit_method_type);

    if (!py_fit)
    {
	delete_fit(fit);

	RETURN_OBJ_ERROR("allocating Py_Fit_method object");
    }

    py_fit->fit = fit;

    return py_fit;
}

static void delete_py_fit(PyObject *self)
{
    Py_Fit_method py_fit = (Py_Fit_method) self;
    Fit_method fit = py_fit->fit;

/*
    printf("in delete_py_fit\n");
*/

    delete_fit(fit);

    PY_FREE(self);
}

static int print_py_fit(PyObject *self, FILE *fp, int flags)
{
    Py_Fit_method py_fit = (Py_Fit_method) self;
    Fit_method fit = py_fit->fit;

    fprintf(fp, "<Fit_method object %p>", fit);

    return 0;
}

static PyObject *getattr_py_fit(PyObject *self, char *name)
{
/*
    Py_Fit_method py_fit = (Py_Fit_method) self;
    Fit_method fit = py_fit->fit;
*/

/*
    if (equal_strings(name, "rp_force_const"))
	return Py_BuildValue("f", fit->rp_force_const);
    else
*/
	return Py_FindMethod(py_handler_methods, self, name);
}

static int setattr_py_fit(PyObject *self, char *name, PyObject *value)
{
    return 0;
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Fit_method_sequence_methods =
{
    Fit_method_length,
    Fit_method_concat,
    Fit_method_repeat,
    Fit_method_item,
    Fit_method_slice,
    Fit_method_ass_item,
    Fit_method_ass_slice
};

static PySequenceMethods Fit_method_sequence_methods =
{
    Fit_method_length,
    0,
    0,
    Fit_method_item,
    0,
    Fit_method_ass_item,
    0
};
*/

static PyTypeObject Fit_method_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "FitMethod", /* name */
    sizeof(struct Py_Fit_method), /* basicsize */
    0, /* itemsize */
    delete_py_fit, /* destructor */
    print_py_fit, /* printfunc */
    getattr_py_fit, /* getattr */
    setattr_py_fit, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Fit_method_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Fit_method(PyObject *self, PyObject *args)
{
    int method;
    float noise;

    if (!PyArg_ParseTuple(args, "if", &method, &noise))
        RETURN_OBJ_ERROR("need two arguments: method, noise");

    return (PyObject *) new_py_fit(method, noise);
}

static CcpnStatus alloc_fit_data_memory(int nparams, int n, float **xw, float **yw,
	float **params_fit, float **params_avg, float **params_dev, float **y_fit)
{
    if (xw != NULL)
	MALLOC(*xw, float, n);

    if (yw != NULL)
	MALLOC(*yw, float, n);

    if (params_avg != NULL)
	MALLOC(*params_avg, float, nparams);

    if (params_dev != NULL)
	MALLOC(*params_dev, float, nparams);

    if (params_fit != NULL)
        MALLOC(*params_fit, float, nparams);

    if (y_fit != NULL)
        MALLOC(*y_fit, float, n);

    return CCPN_OK;
}

static void free_fit_data_memory(float *x, float *y, float *x_dev, float *y_dev,
		float *xw, float *yw, float *params_fit,
		float *params_avg, float *params_dev, float *y_fit)
{
    FREE(x, float);
    FREE(y, float);
    FREE(x_dev, float);
    FREE(y_dev, float);
    FREE(xw, float);
    FREE(yw, float);
    FREE(params_fit, float);
    FREE(params_avg, float);
    FREE(params_dev, float);
    FREE(y_fit, float);
}

static PyObject *runFit(PyObject *self, PyObject *args)
{
    int n1, n2, method, niter, nparams;
    float noise;
    float *x, *y;
    float *params_fit, *params_dev, *y_fit, chisq;
    PyObject *x_obj, *y_obj;
    PyObject *params_obj, *params_dev_obj, *y_fit_obj, *result;
    Fit_method fit;
    Line error_msg;
    CcpnStatus status;

    x = y = params_fit = y_fit = NULL;

    if (!PyArg_ParseTuple(args, "ifOO", &method, &noise, &x_obj, &y_obj))
        RETURN_OBJ_ERROR("need four arguments: method, noise, x, y");

    if ((method < 0) || (method >= NFIT_METHODS))
    {
	sprintf(error_msg, "method must be between 0 and %d", NFIT_METHODS-1);
	RETURN_OBJ_ERROR(error_msg);
    }

    if (get_python_float_alloc_array(x_obj, &n1, &x, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    if (get_python_float_alloc_array(y_obj, &n2, &y, error_msg) == CCPN_ERROR)
    {
        free_fit_data_memory(x, y, NULL, NULL, NULL, NULL,
				params_fit, NULL, NULL, y_fit);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (n1 != n2)
    {
        sprintf(error_msg, "x any y must both be of size %d", n1);
        free_fit_data_memory(x, y, NULL, NULL, NULL, NULL,
				params_fit, NULL, NULL, y_fit);
        RETURN_OBJ_ERROR(error_msg);
    }
 
    nparams = get_method_nparams(method);
    if (alloc_fit_data_memory(nparams, n1, NULL, NULL, &params_fit, NULL,
					&params_dev, &y_fit) == CCPN_ERROR)
    {
        free_fit_data_memory(x, y, NULL, NULL, NULL, NULL,
				params_fit, NULL, params_dev, y_fit);
        RETURN_OBJ_ERROR("allocating fit memory");
    }

    fit = new_fit(method, noise);

    if (!fit)
    {
        free_fit_data_memory(x, y, NULL, NULL, NULL, NULL,
				params_fit, NULL, params_dev, y_fit);
	RETURN_OBJ_ERROR("allocating Fit_method object");
    }

    status = run_fit(fit, n1, x, y, params_fit, params_dev, y_fit, &chisq, error_msg);

    if (status == CCPN_OK)
    {
	params_obj = get_python_float_list(nparams, params_fit);
	params_dev_obj = get_python_float_list(nparams, params_dev);
        y_fit_obj  = get_python_float_list(n1, y_fit);
    }

    free_fit_data_memory(x, y, NULL, NULL, NULL, NULL,
				params_fit, NULL, params_dev, y_fit);
    delete_fit(fit);

    if (status == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    result = Py_BuildValue("(OOOf)", params_obj, params_dev_obj, y_fit_obj, chisq);
    Py_DECREF(params_obj);
    Py_DECREF(params_dev_obj);
    Py_DECREF(y_fit_obj);

    return result;
}

static PyObject *fitData(PyObject *self, PyObject *args)
{
    int n1, n2, method, niter, nparams;
    float noise;
    float *x, *y, *x_dev, *y_dev, *xw, *yw;
    float *params_fit, *params_avg, *params_dev, *y_fit, chisq;
    PyObject *x_obj, *y_obj, *x_dev_obj, *y_dev_obj;
    PyObject *params_obj, *params_dev_obj, *y_fit_obj, *result;
    Fit_method fit;
    Line error_msg;
    CcpnStatus status;

    x = y = x_dev = y_dev = xw = yw = params_fit = params_avg = params_dev = y_fit = NULL;

    if (!PyArg_ParseTuple(args, "iifOOOO", &method, &niter, &noise, &x_obj,
						&y_obj, &x_dev_obj, &y_dev_obj))
        RETURN_OBJ_ERROR("need seven arguments: method, niter, noise, x, y, x_dev, y_dev");

    if ((method < 0) || (method >= NFIT_METHODS))
    {
	sprintf(error_msg, "method must be between 0 and %d", NFIT_METHODS-1);
	RETURN_OBJ_ERROR(error_msg);
    }

    if (get_python_float_alloc_array(x_obj, &n1, &x, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    if (get_python_float_alloc_array(y_obj, &n2, &y, error_msg) == CCPN_ERROR)
    {
        free_fit_data_memory(x, y, x_dev, y_dev, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (n1 != n2)
    {
        sprintf(error_msg, "x any y must both be of size %d", n1);
        free_fit_data_memory(x, y, x_dev, y_dev, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
        RETURN_OBJ_ERROR(error_msg);
    }
 
    if (get_python_float_alloc_array(x_dev_obj, &n2, &x_dev, error_msg) == CCPN_ERROR)
    {
        free_fit_data_memory(x, y, x_dev, y_dev, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (n1 != n2)
    {
        sprintf(error_msg, "x any x_dev must both be of size %d", n1);
        free_fit_data_memory(x, y, x_dev, y_dev, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
        RETURN_OBJ_ERROR(error_msg);
    }
 
    if (get_python_float_alloc_array(y_dev_obj, &n2, &y_dev, error_msg) == CCPN_ERROR)
    {
        free_fit_data_memory(x, y, x_dev, y_dev, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (n1 != n2)
    {
        sprintf(error_msg, "y any y_dev must both be of size %d", n1);
        free_fit_data_memory(x, y, x_dev, y_dev, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
        RETURN_OBJ_ERROR(error_msg);
    }
 
    nparams = get_method_nparams(method);
    if (alloc_fit_data_memory(nparams, n1, &xw, &yw, &params_fit, &params_avg,
					&params_dev, &y_fit) == CCPN_ERROR)
    {
        free_fit_data_memory(x, y, x_dev, y_dev, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
        RETURN_OBJ_ERROR("allocating fit memory");
    }

    fit = new_fit(method, noise);

    if (!fit)
    {
        free_fit_data_memory(x, y, x_dev, y_dev, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
	RETURN_OBJ_ERROR("allocating Fit_method object");
    }

    status = run_fit_data(fit, niter, n1, x, y, x_dev, y_dev, xw, yw,
		params_fit, params_avg, params_dev, y_fit, &chisq, error_msg);

    if (status == CCPN_OK)
    {
	params_obj = get_python_float_list(nparams, params_fit);
	params_dev_obj = get_python_float_list(nparams, params_dev);
        y_fit_obj  = get_python_float_list(n1, y_fit);
    }

    free_fit_data_memory(x, y, x_dev, y_dev, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
    delete_fit(fit);

    if (status == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    result = Py_BuildValue("(OOOf)", params_obj, params_dev_obj, y_fit_obj, chisq);
    Py_DECREF(params_obj);
    Py_DECREF(params_dev_obj);
    Py_DECREF(y_fit_obj);

    return result;
}

static PyObject *bootstrapData(PyObject *self, PyObject *args)
{
    int n1, n2, method, niter, nparams;
    float noise;
    float *x, *y, *xw, *yw;
    float *params_fit, *params_avg, *params_dev, *y_fit, chisq;
    PyObject *x_obj, *y_obj;
    PyObject *params_obj, *params_dev_obj, *y_fit_obj, *result;
    Fit_method fit;
    Line error_msg;
    CcpnStatus status;

    x = y = xw = yw = params_fit = params_avg = params_dev = y_fit = NULL;

    if (!PyArg_ParseTuple(args, "iifOO", &method, &niter, &noise, &x_obj, &y_obj))
        RETURN_OBJ_ERROR("need five arguments: method, niter, noise, x, y");

    if ((method < 0) || (method >= NFIT_METHODS))
    {
	sprintf(error_msg, "method must be between 0 and %d", NFIT_METHODS-1);
	RETURN_OBJ_ERROR(error_msg);
    }

    if (get_python_float_alloc_array(x_obj, &n1, &x, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    if (get_python_float_alloc_array(y_obj, &n2, &y, error_msg) == CCPN_ERROR)
    {
        free_fit_data_memory(x, y, NULL, NULL, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (n1 != n2)
    {
        sprintf(error_msg, "x any y must both be of size %d", n1);
        free_fit_data_memory(x, y, NULL, NULL, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
        RETURN_OBJ_ERROR(error_msg);
    }
 
    nparams = get_method_nparams(method);
    if (alloc_fit_data_memory(nparams, n1, &xw, &yw, &params_fit, &params_avg,
					&params_dev, &y_fit) == CCPN_ERROR)
    {
        free_fit_data_memory(x, y, NULL, NULL, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
        RETURN_OBJ_ERROR("allocating fit memory");
    }

    fit = new_fit(method, noise);

    if (!fit)
    {
        free_fit_data_memory(x, y, NULL, NULL, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
	RETURN_OBJ_ERROR("allocating Fit_method object");
    }

    status = bootstrap_fit_data(fit, niter, n1, x, y, xw, yw,
		params_fit, params_avg, params_dev, y_fit, &chisq, error_msg);

    if (status == CCPN_OK)
    {
	params_obj = get_python_float_list(nparams, params_fit);
	params_dev_obj = get_python_float_list(nparams, params_dev);
        y_fit_obj  = get_python_float_list(n1, y_fit);
    }

    free_fit_data_memory(x, y, NULL, NULL, xw, yw,
				params_fit, params_avg, params_dev, y_fit);
    delete_fit(fit);

    if (status == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    result = Py_BuildValue("(OOOf)", params_obj, params_dev_obj, y_fit_obj, chisq);
    Py_DECREF(params_obj);
    Py_DECREF(params_dev_obj);
    Py_DECREF(y_fit_obj);

    return result;
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Fit_method_type_methods[] =
{
    { "FitMethod",	(PyCFunction) init_Py_Fit_method,	METH_VARARGS },
    { "runFit",		(PyCFunction) runFit,			METH_VARARGS },
    { "fitData",	(PyCFunction) fitData,			METH_VARARGS },
    { "bootstrapData",	(PyCFunction) bootstrapData,			METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import Fit_method" in 
* a Python program. The function is usually called "initFit_method": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initFitMethod(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Fit_method_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("FitMethod", Fit_method_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("FitMethod.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module FitMethod");
}

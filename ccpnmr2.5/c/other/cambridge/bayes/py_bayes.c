/* ====================================================== 
======================COPYRIGHT/LICENSE START==========================

py_bayes.c: Part of the Bayes program

Copyright (C) 2003-2010 Daniel O'Donovan (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the author: djo35@bioc.cam.ac.uk
=======================================================================
===========================REFERENCE START=============================
===========================REFERENCE END===============================

*/
/* = Once again, shamelessly ripped out of block_file.c = */
/* Created by Daniel O'Donovan on 2008-09-22.
 * Copyright (c) 2008 University of Cambridge. 
 * All rights reserved.                                   */
/* ====================================================== */

#include "py_bayes.h"

#include "python_util.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject BayesPeakSeparator_type;

Bool is_py_bayes(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Block_file_type
    return (obj->ob_type == &Block_file_type);
*/
    return valid_py_object(obj, &BayesPeakSeparator_type);
}

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

static PyTypeObject BayesPeakSeparator_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "BayesPeakSeparator", /* name */
    sizeof(struct Py_BayesPeakSeparator), /* basicsize */
    0, /* itemsize */
    0, /* destructor */
    0, /* printfunc */
    0, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *test_bayes(PyObject *self, PyObject *args)
{
    PyObject * obj_in, * obj_out;
    int i, n;
    double * array;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &obj_in))
        RETURN_OBJ_ERROR("oops\n");

    if (get_python_double_alloc_array(obj_in, &n, &array,
                                        error_msg) == CCPN_ERROR)
    {
        sprintf(error_msg, "max_sigma must be list or tuple");
        RETURN_OBJ_ERROR(error_msg);
    }

    printf("parsed and got arrays\n");

    for( i = 0; i < n; i++ )
    {
        array[i] *= 2;
    }

    obj_out = get_python_double_list(n, array);

    return obj_out;
}


static PyObject *run_bayes(PyObject *self, PyObject *args)
{
    int n;

    /* input from Python Args */
    char * spc_file_obj;
    int ndim_obj, endianess_obj;
    PyObject *npoints_obj, * block_size_obj, *sample_start_obj, *sample_end_obj;
    float max_height_obj, min_height_obj;
    PyObject *max_sigma_obj, *min_sigma_obj;
    double max_Q_obj, min_Q_obj;
    PyObject *dim_wrapped_int_obj;
    int shape_obj, pos_peaks_obj, min_atoms_obj, max_atoms_obj;
    double rate_obj;

    /* to be set locally */
    char * spc_file;
    int ndim, endianess, *npoints, * block_size, *sample_start, *sample_end;
    float max_height, min_height;
    double *max_sigma, *min_sigma;
    double max_Q, min_Q;
    int *dim_wrapped_int;
    int shape, pos_peaks, min_atoms, max_atoms;
    double rate;


    PyObject *py_list;
    
    Line error_msg;
    int return_code;

    /* parse the python args and test */
    if (!PyArg_ParseTuple(args, "siiOOOOffOOddOiiiid", 
                        &spc_file_obj,                                              /* char * */
                        &ndim_obj,                                                  /* int */
                        &endianess_obj,                                             /* int */
                        &npoints_obj, &block_size_obj,                              /* PyObject * (int) */
                        &sample_start_obj, &sample_end_obj,
                        &max_height_obj, &min_height_obj,                           /* float */
                        &max_sigma_obj, &min_sigma_obj,                             /* PyObject * (double) */                
                        &max_Q_obj, &min_Q_obj,                                     /* double */
                        &dim_wrapped_int_obj,                                       /* PyObject * (int) */
                        &shape_obj, &pos_peaks_obj, &min_atoms_obj, &max_atoms_obj, /* int */
                        &rate_obj))                                              
        RETURN_OBJ_ERROR("run_bayes needs 19 arguments!");

    /* set everything not an array (easy and can be done in ParseTuple)*/
    spc_file = spc_file_obj;
    ndim = ndim_obj;
    endianess = endianess_obj;
    max_height = max_height_obj;
    min_height = min_height_obj;
    max_Q = max_Q_obj;
    min_Q = min_Q_obj;
    shape = shape_obj;
    pos_peaks = pos_peaks_obj;
    min_atoms = min_atoms_obj;
    max_atoms = max_atoms_obj;
    rate = rate_obj;

    /* create the int arrays */
    if (get_python_int_alloc_array(npoints_obj, &n, &npoints,
                                    error_msg) == CCPN_ERROR)
    {
        sprintf(error_msg, "points must be list or tuple");
        RETURN_OBJ_ERROR(error_msg);
    }

    if (get_python_int_alloc_array(block_size_obj, &n, &block_size,
                                    error_msg) == CCPN_ERROR)
    {
        sprintf(error_msg, "block_size must be list or tuple");
        RETURN_OBJ_ERROR(error_msg);
    }

    if (get_python_int_alloc_array(sample_start_obj, &n, &sample_start,
                                    error_msg) == CCPN_ERROR)
    {
        sprintf(error_msg, "sample_start must be list or tuple");
        RETURN_OBJ_ERROR(error_msg);
    }

    if (get_python_int_alloc_array(sample_end_obj, &n, &sample_end,
                                    error_msg) == CCPN_ERROR)
    {
        sprintf(error_msg, "sample_end must be list or tuple");
        RETURN_OBJ_ERROR(error_msg);
    }

    /* these are double arrays */
    if (get_python_double_alloc_array(max_sigma_obj, &n, &max_sigma,
                                    error_msg) == CCPN_ERROR)
    {
        sprintf(error_msg, "max_sigma must be list or tuple");
        RETURN_OBJ_ERROR(error_msg);
    }

    if (get_python_double_alloc_array(min_sigma_obj, &n, &min_sigma,
                                    error_msg) == CCPN_ERROR)
    {
        sprintf(error_msg, "min_sigma must be list or tuple");
        RETURN_OBJ_ERROR(error_msg);
    }

    if (get_python_int_alloc_array(dim_wrapped_int_obj, &n, &dim_wrapped_int,
                                    error_msg) == CCPN_ERROR)
    {
        sprintf(error_msg, "dim_wrapped must be list or tuple");
        RETURN_OBJ_ERROR(error_msg);
    }

    py_list = (PyObject* ) PyList_New( 0 );

    return_code = bayesNMR(spc_file, ndim, endianess, npoints, block_size, 
                sample_start, sample_end, max_height, min_height, max_sigma, min_sigma, max_Q, min_Q,
                dim_wrapped_int, shape, pos_peaks, min_atoms, max_atoms, rate, py_list);

    /* 8 free calls */
    FREE_TYPE(npoints, int);
    FREE_TYPE(block_size, int);
    FREE_TYPE(sample_start, int);
    FREE_TYPE(sample_end, int);

    FREE_TYPE(max_sigma, double);
    FREE_TYPE(min_sigma, double);

    FREE_TYPE(dim_wrapped_int, int);

    return py_list;
} 


/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef BayesPeakSeparator_type_methods[] =
{
    { "test_bayes",     test_bayes, METH_VARARGS    },
    { "run_bayes",      run_bayes,  METH_VARARGS    },
    { NULL,             NULL,       0               }         /* sentinel */
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import BlockFile" in 
* a Python program. The function is usually called "initBlock_file": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initBayesPeakSeparator( void )
{
    PyObject *m, *d;

#ifdef WIN64
    BayesPeakSeparator_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("BayesPeakSeparator", BayesPeakSeparator_type_methods);

    /* add symbolic constants to the module */
    d = PyModule_GetDict(m);
    ErrorObject = Py_BuildValue("s", "BayesPeakSeparator.error");
    PyDict_SetItemString(d, "error", ErrorObject);

    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module BayesPeakSeparator");
}

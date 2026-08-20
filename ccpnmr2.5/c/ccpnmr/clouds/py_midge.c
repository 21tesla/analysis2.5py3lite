
/*
======================COPYRIGHT/LICENSE START==========================

py_midge.c: Part of the CcpNmr Analysis program

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
#include "py_midge.h"

#include "midge.h"

#include "python_util.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Midge_type;

Bool is_py_midge(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Midge_type
    return (obj->ob_type == &Midge_type);
*/
    return valid_py_object(obj, &Midge_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *run(PyObject *self, PyObject *args)
{
    Py_Midge py_midge = (Py_Midge) self;
    Midge midge = py_midge->midge;
    int n = midge->n, n1, n2, n15_lab, c13_lab, max_iter;
    PyObject *amat_obj, *rmat_obj;
    double **amat, **rmat;
    float sf, tmix, tcor, rleak, err;
    Bool n15_labelled, c13_labelled;
    Line error_msg;
    CcpnStatus status;
 
    if (!PyArg_ParseTuple(args, "OOiffffii", &amat_obj, &rmat_obj,
			&max_iter, &sf, &tmix, &tcor, &rleak, &n15_lab, &c13_lab))
        RETURN_OBJ_ERROR("need nine arguments: amat, rmat, max_iter, sf, tmix, tcor, rleak, n15_labelled, c13_labelled");

    if (get_python_double_alloc_matrix(amat_obj, &n1, &n2, &amat, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    if ((n1 != n) || (n2 != n))
    {
        sprintf(error_msg, "amat must be square matrix of size %d", n);
	FREE2(amat, double, n1);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (get_python_double_alloc_matrix(rmat_obj, &n1, &n2, &rmat, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    if ((n1 != n) || (n2 != n))
    {
        sprintf(error_msg, "rmat must be square matrix of size %d", n);
	FREE2(amat, double, n);
	FREE2(rmat, double, n1);
        RETURN_OBJ_ERROR(error_msg);
    }

    n15_labelled = n15_lab;
    c13_labelled = c13_lab;
    status = run_midge(midge, amat, rmat, max_iter, sf, tmix, tcor, rleak,
		n15_labelled, c13_labelled, &err, error_msg);

    if (status == CCPN_OK)
    {
	status = set_python_double_matrix(amat_obj, n, n, amat, error_msg);
        if (status == CCPN_OK)
	    status = set_python_double_matrix(rmat_obj, n, n, rmat, error_msg);
    }

    FREE2(amat, double, n);
    FREE2(rmat, double, n);

    if (status == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    return Py_BuildValue("f", err);
}

static struct PyMethodDef py_handler_methods[] =
{
    { "run",		run,		METH_VARARGS },
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static Py_Midge new_py_midge(int n, int *nhs, int *types)
{
    Py_Midge py_midge;
    Midge midge;

    midge = new_midge(n, nhs, types);

    if (!midge)
	 RETURN_OBJ_ERROR("allocating Midge object");

    PY_MALLOC(py_midge, struct Py_Midge, &Midge_type);

    if (!py_midge)
    {
	delete_midge(midge);

	RETURN_OBJ_ERROR("allocating Py_Midge object");
    }

    py_midge->midge = midge;

    return py_midge;
}

static void delete_py_midge(PyObject *self)
{
    Py_Midge py_midge = (Py_Midge) self;
    Midge midge = py_midge->midge;

/*
    printf("in delete_py_midge\n");
*/

    delete_midge(midge);

    PY_FREE(self);
}

static int print_py_midge(PyObject *self, FILE *fp, int flags)
{
    Py_Midge py_midge = (Py_Midge) self;
    Midge midge = py_midge->midge;

    fprintf(fp, "<Midge object %d>", (int) midge);

    return 0;
}

static PyObject *getattr_py_midge(PyObject *self, char *name)
{
    Py_Midge py_midge = (Py_Midge) self;
/*
    Midge midge = py_midge->midge;

    if (equal_strings(name, "rp_force_const"))
	return Py_BuildValue("f", midge->rp_force_const);
    else
*/
	return Py_FindMethod(py_handler_methods, self, name);
}

static int setattr_py_midge(PyObject *self, char *name, PyObject *value)
{
    return 0;
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Midge_sequence_methods =
{
    Midge_length,
    Midge_concat,
    Midge_repeat,
    Midge_item,
    Midge_slice,
    Midge_ass_item,
    Midge_ass_slice
};

static PySequenceMethods Midge_sequence_methods =
{
    Midge_length,
    0,
    0,
    Midge_item,
    0,
    Midge_ass_item,
    0
};
*/

static PyTypeObject Midge_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "Midge", /* name */
    sizeof(struct Py_Midge), /* basicsize */
    0, /* itemsize */
    delete_py_midge, /* destructor */
    print_py_midge, /* printfunc */
    getattr_py_midge, /* getattr */
    setattr_py_midge, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Midge_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Midge(PyObject *self, PyObject *args)
{
    int m, n, *nhs, *types;
    PyObject *nhs_obj, *types_obj;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "OO", &nhs_obj, &types_obj))
        RETURN_OBJ_ERROR("need two arguments: nhs, types");

    if (get_python_int_alloc_array(nhs_obj, &n, &nhs, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);
	
    if (get_python_int_alloc_array(types_obj, &m, &types, error_msg) == CCPN_ERROR)
    {
	FREE(nhs, int);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (m != n)
    {
	FREE(nhs, int);
	FREE(types, int);
	RETURN_OBJ_ERROR("nhs and types must be lists of same length");
    }

    return (PyObject *) new_py_midge(n, nhs, types);
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Midge_type_methods[] =
{
    { "Midge",	(PyCFunction) init_Py_Midge,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import Midge" in 
* a Python program. The function is usually called "initMidge": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initMidge(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Midge_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("Midge", Midge_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("Midge.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module Midge");
}

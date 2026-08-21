
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

    py_midge = (Py_Midge) PyObject_New(struct Py_Midge, &Midge_type);

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

    Py_TYPE(self)->tp_free(self);
}

static PyObject *repr_py_midge(PyObject *self)
{
    Py_Midge py_midge = (Py_Midge) self;
    Midge midge = py_midge->midge;

    char buf[48];
    snprintf(buf, sizeof buf, "<Midge object %p>", (void *) midge);
    return PyUnicode_FromString(buf);
}

static PyObject *getattr_py_midge(PyObject *self, PyObject *attr_name)
{
    return PyObject_GenericGetAttr(self, attr_name);
}

static int setattr_py_midge(PyObject *self, PyObject *attr_name, PyObject *value)
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
    PyVarObject_HEAD_INIT(NULL, 0)
    "Midge",                                  /* tp_name */
    sizeof(struct Py_Midge),                  /* tp_basicsize */
    0,                                        /* tp_itemsize */
    (destructor) delete_py_midge,             /* tp_dealloc */
    0,                                        /* tp_vectorcall */
    0,                                        /* tp_getattr */
    0,                                        /* tp_setattr */
    0,                                        /* tp_as_async */
    (reprfunc) repr_py_midge,                 /* tp_repr */
    0,                                        /* tp_as_number */
    0,                                        /* tp_as_sequence */
    0,                                        /* tp_as_mapping */
    0,                                        /* tp_hash */
    0,                                        /* tp_call */
    0,                                        /* tp_str */
    (getattrofunc) getattr_py_midge,          /* tp_getattro */
    (setattrofunc) setattr_py_midge,          /* tp_setattro */
    0,                                        /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                       /* tp_flags */
    "Midge -- NOE relaxation simulation",     /* tp_doc */
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

static struct PyModuleDef midge_module_def =
{
    PyModuleDef_HEAD_INIT,
    "Midge",
    "CCPNMR Midge module (Python 3 compatible)",
    -1,
    Midge_type_methods
};

PyMODINIT_FUNC PyInit_Midge(void)
{
    if (PyType_Ready(&Midge_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&midge_module_def);
    if (!m)
        return NULL;

    ErrorObject = PyErr_NewException("Midge.error", NULL, NULL);
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

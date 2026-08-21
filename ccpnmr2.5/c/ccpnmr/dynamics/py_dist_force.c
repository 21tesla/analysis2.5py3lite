
/*
======================COPYRIGHT/LICENSE START==========================

py_dist_force.c: Part of the CcpNmr Analysis program

Copyright (C) 2005 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in ../../../license/CCPN.license.

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
#include "py_dist_force.h"

#include "python_util.h"

#include "utility.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Dist_force_type;

Bool is_py_dist_force(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Dist_force_type
    return (obj->ob_type == &Dist_force_type);
*/
    return valid_py_object(obj, &Dist_force_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static struct PyMethodDef py_handler_methods[] =
{
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static Py_Dist_force new_py_dist_force(float force_const, float exponent,
		float soft_exponent, float r_switch, float asymptote, float dist_power)
{
    Py_Dist_force obj;
    Dist_force dist_force;

    if (exponent != 2)
	 RETURN_OBJ_ERROR("exponent must be 2 (for now)");

    if (soft_exponent != 1)
	 RETURN_OBJ_ERROR("soft_exponent must be 1 (for now)");

    if (dist_power != 4)
	 RETURN_OBJ_ERROR("dist_power must be 4 (for now)");

    dist_force = new_dist_force(force_const, exponent, soft_exponent,
					r_switch, asymptote, dist_power);

    if (!dist_force)
	 RETURN_OBJ_ERROR("allocating Dist_force object");

    obj = (Py_Dist_force) PyObject_New(struct Py_Dist_force, &Dist_force_type);

    if (!obj)
    {
	delete_dist_force(dist_force);

	RETURN_OBJ_ERROR("allocating Py_Dist_force object");
    }

    obj->dist_force = dist_force;

    return obj;
}

static void delete_py_dist_force(PyObject *self)
{
    Py_Dist_force obj = (Py_Dist_force) self;
    Dist_force dist_force = obj->dist_force;

/*
    printf("in delete_py_dist_force\n");
*/

    delete_dist_force(dist_force);

    Py_TYPE(self)->tp_free(self);
}

/*
static int print_py_dist_force(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_dist_force(PyObject *self, PyObject *attr_name)
{
    Py_Dist_force obj = (Py_Dist_force) self;
    Dist_force dist_force = obj->dist_force;

    if (PyUnicode_Check(attr_name))
    {
        const char *name = PyUnicode_AsUTF8(attr_name);
        if (name)
        {
            if (strcmp(name, "force_const") == 0)
                return Py_BuildValue("f", dist_force->force_const);
            else if (strcmp(name, "exponent") == 0)
                return Py_BuildValue("f", dist_force->exponent);
            else if (strcmp(name, "soft_exponent") == 0)
                return Py_BuildValue("f", dist_force->soft_exponent);
            else if (strcmp(name, "r_switch") == 0)
                return Py_BuildValue("f", dist_force->r_switch);
            else if (strcmp(name, "asymptote") == 0)
                return Py_BuildValue("f", dist_force->asymptote);
            else if (strcmp(name, "dist_power") == 0)
                return Py_BuildValue("f", dist_force->dist_power);
        }
    }
    return PyObject_GenericGetAttr(self, attr_name);
}

static int setattr_py_dist_force(PyObject *self, PyObject *attr_name, PyObject *value)
{
    Py_Dist_force obj = (Py_Dist_force) self;
    Dist_force dist_force = obj->dist_force;
    float v = (float) PyFloat_AsDouble(value);

    if (PyErr_Occurred())
	RETURN_INT_ERROR("must have float value");

    if (!PyUnicode_Check(attr_name))
    {
        PyErr_SetString(ErrorObject, "attribute name must be a string");
        return -1;
    }
    const char *name = PyUnicode_AsUTF8(attr_name);
    if (!name)
        return -1;

    if (strcmp(name, "force_const") == 0)
    {
	dist_force->force_const = v;
    }
    else if (strcmp(name, "exponent") == 0)
    {
        if (v != 2)
	    RETURN_INT_ERROR("exponent must be 2 (for now)");
	dist_force->exponent = v;
    }
    else if (strcmp(name, "soft_exponent") == 0)
    {
        if (v != 1)
	    RETURN_INT_ERROR("soft_exponent must be 1 (for now)");
	dist_force->soft_exponent = v;
    }
    else if (strcmp(name, "r_switch") == 0)
    {
	dist_force->r_switch = v;
    }
    else if (strcmp(name, "asymptote") == 0)
    {
	dist_force->asymptote = v;
    }
    else if (strcmp(name, "dist_power") == 0)
    {
        if (v != 4)
	    RETURN_INT_ERROR("dist_power must be 4 (for now)");
	dist_force->dist_power = v;
    }
    else
    {
	RETURN_INT_ERROR("unknown attribute name");
    }

    return 0;
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Dist_force_sequence_methods =
{
    Dist_force_length,
    Dist_force_concat,
    Dist_force_repeat,
    Dist_force_item,
    Dist_force_slice,
    Dist_force_ass_item,
    Dist_force_ass_slice
};

static PySequenceMethods Dist_force_sequence_methods =
{
    Dist_force_length,
    0,
    0,
    Dist_force_item,
    0,
    Dist_force_ass_item,
    0
};
*/

static PyTypeObject Dist_force_type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    "DistForce",                        /* tp_name */
    sizeof(struct Py_Dist_force),       /* tp_basicsize */
    0,                                  /* tp_itemsize */
    (destructor) delete_py_dist_force,  /* tp_dealloc */
    0,                                  /* tp_vectorcall */
    0,                                  /* tp_getattr */
    0,                                  /* tp_setattr */
    0,                                  /* tp_as_async */
    0,                                  /* tp_repr */
    0,                                  /* tp_as_number */
    0,                                  /* tp_as_sequence */
    0,                                  /* tp_as_mapping */
    0,                                  /* tp_hash */
    0,                                  /* tp_call */
    0,                                  /* tp_str */
    (getattrofunc) getattr_py_dist_force,  /* tp_getattro */
    (setattrofunc) setattr_py_dist_force,  /* tp_setattro */
    0,                                  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                 /* tp_flags */
    "DistForce -- force curve parameters", /* tp_doc */
    0,                                  /* tp_traverse */
    0,                                  /* tp_clear */
    0,                                  /* tp_richcompare */
    0,                                  /* tp_weaklistoffset */
    0,                                  /* tp_iter */
    0,                                  /* tp_iternext */
    py_handler_methods,                 /* tp_methods */
    0,                                  /* tp_members */
    0,                                  /* tp_getset */
    0,                                  /* tp_base */
    0,                                  /* tp_dict */
    0,                                  /* tp_descr_get */
    0,                                  /* tp_descr_set */
    0,                                  /* tp_dictoffset */
    0,                                  /* tp_init */
    0,                                  /* tp_alloc */
    0,                                  /* tp_new */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Dist_force(PyObject *self, PyObject *args)
{
    float force_const, exponent, soft_exponent, r_switch, asymptote, dist_power;

    if (!PyArg_ParseTuple(args, "ffffff", &force_const, &exponent,
		&soft_exponent, &r_switch, &asymptote, &dist_power))
        RETURN_OBJ_ERROR("must have arguments: force_const, exponent, soft_exponent, r_switch, asymptote, dist_power");

    return (PyObject *) new_py_dist_force(force_const,
		exponent, soft_exponent, r_switch, asymptote, dist_power);
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Dist_force_type_methods[] =
{
    { "DistForce",	(PyCFunction) init_Py_Dist_force,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import DistForce" in 
* a Python program. The function is usually called "initDist_force": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

static struct PyModuleDef DyDistForce_module_def =
{
    PyModuleDef_HEAD_INIT,
    "DyDistForce",
    "CCPNMR DyDistForce module (Python 3 compatible)",
    -1,
    Dist_force_type_methods
};

PyMODINIT_FUNC PyInit_DyDistForce(void)
{
    if (PyType_Ready(&Dist_force_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&DyDistForce_module_def);
    if (!m)
        return NULL;

    ErrorObject = PyErr_NewException("DyDistForce.error", NULL, NULL);
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

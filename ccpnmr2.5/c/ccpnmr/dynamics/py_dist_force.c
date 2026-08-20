
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

    PY_MALLOC(obj, struct Py_Dist_force, &Dist_force_type);

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

    PY_FREE(self);
}

/*
static int print_py_dist_force(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_dist_force(PyObject *self, char *name)
{
    Py_Dist_force obj = (Py_Dist_force) self;
    Dist_force dist_force = obj->dist_force;

    if (equal_strings(name, "force_const"))
	return Py_BuildValue("f", dist_force->force_const);
    else if (equal_strings(name, "exponent"))
	return Py_BuildValue("f", dist_force->exponent);
    else if (equal_strings(name, "soft_exponent"))
	return Py_BuildValue("f", dist_force->soft_exponent);
    else if (equal_strings(name, "r_switch"))
	return Py_BuildValue("f", dist_force->r_switch);
    else if (equal_strings(name, "asymptote"))
	return Py_BuildValue("f", dist_force->asymptote);
    else if (equal_strings(name, "dist_power"))
	return Py_BuildValue("f", dist_force->dist_power);
    else
	return Py_FindMethod(py_handler_methods, self, name);
}

static int setattr_py_dist_force(PyObject *self, char *name, PyObject *value)
{
    Py_Dist_force obj = (Py_Dist_force) self;
    Dist_force dist_force = obj->dist_force;
    float v = (float) PyFloat_AsDouble(value);

    if (PyErr_Occurred())
	RETURN_INT_ERROR("must have float value");

    if (equal_strings(name, "force_const"))
    {
	dist_force->force_const = v;
    }
    else if (equal_strings(name, "exponent"))
    {
        if (v != 2)
	    RETURN_INT_ERROR("exponent must be 2 (for now)");
	dist_force->exponent = v;
    }
    else if (equal_strings(name, "soft_exponent"))
    {
        if (v != 1)
	    RETURN_INT_ERROR("soft_exponent must be 1 (for now)");
	dist_force->soft_exponent = v;
    }
    else if (equal_strings(name, "r_switch"))
    {
	dist_force->r_switch = v;
    }
    else if (equal_strings(name, "asymptote"))
    {
	dist_force->asymptote = v;
    }
    else if (equal_strings(name, "dist_power"))
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
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "DistForce", /* name */
    sizeof(struct Py_Dist_force), /* basicsize */
    0, /* itemsize */
    delete_py_dist_force, /* destructor */
    0, /* printfunc */
    getattr_py_dist_force, /* getattr */
    setattr_py_dist_force, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Dist_force_sequence_methods*/ /* PySequenceMethods */
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

PY_MOD_INIT_FUNC initDyDistForce(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Dist_force_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("DyDistForce", Dist_force_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("DyDistForce.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module DyDistForce");
}

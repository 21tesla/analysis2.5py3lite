
/*
======================COPYRIGHT/LICENSE START==========================

py_dist_constraint.c: Part of the CcpNmr Analysis program

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
#include "py_dist_constraint.h"

#include "python_util.h"

#include "utility.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Dist_constraint_type;

Bool is_py_dist_constraint(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Dist_constraint_type
    return (obj->ob_type == &Dist_constraint_type);
*/
    return valid_py_object(obj, &Dist_constraint_type);
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

Py_Dist_constraint new_py_dist_constraint(int atom0, int atom1,
				float dist_lower, float dist_upper)
{
    Py_Dist_constraint py_dist_constraint;
    Dist_constraint dist_constraint;

    dist_constraint = new_dist_constraint(atom0, atom1, dist_lower, dist_upper);

    if (!dist_constraint)
	 RETURN_OBJ_ERROR("allocating Dist_constraint object");

    PY_MALLOC(py_dist_constraint, struct Py_Dist_constraint, &Dist_constraint_type);

    if (!py_dist_constraint)
    {
	delete_dist_constraint(dist_constraint);

	RETURN_OBJ_ERROR("allocating Py_Dist_constraint object");
    }

    py_dist_constraint->dist_constraint = dist_constraint;

    return py_dist_constraint;
}

static void delete_py_dist_constraint(PyObject *self)
{
    Py_Dist_constraint py_dist_constraint = (Py_Dist_constraint) self;
    Dist_constraint dist_constraint = py_dist_constraint->dist_constraint;

/*
    printf("in delete_py_dist_constraint\n");
*/

    delete_dist_constraint(dist_constraint);

    PY_FREE(self);
}

static int print_py_dist_constraint(PyObject *self, FILE *fp, int flags)
{
    Py_Dist_constraint py_dist_constraint = (Py_Dist_constraint) self;
    Dist_constraint dist_constraint = py_dist_constraint->dist_constraint;

    fprintf(fp, "<atom0=%d, atom1=%d, dist_lower=%3.2e, dist_upper=%3.2e>",
                dist_constraint->atom0, dist_constraint->atom1,
                dist_constraint->dist_lower, dist_constraint->dist_upper);

    return 0;
}

static PyObject *getattr_py_dist_constraint(PyObject *self, char *name)
{
    Py_Dist_constraint py_dist_constraint = (Py_Dist_constraint) self;
    Dist_constraint dist_constraint = py_dist_constraint->dist_constraint;

    if (equal_strings(name, "atom0"))
	return Py_BuildValue("i", dist_constraint->atom0);
    else if (equal_strings(name, "atom1"))
	return Py_BuildValue("i", dist_constraint->atom1);
    else if (equal_strings(name, "dist_lower"))
	return Py_BuildValue("f", dist_constraint->dist_lower);
    else if (equal_strings(name, "dist_upper"))
	return Py_BuildValue("f", dist_constraint->dist_upper);
    else
	return Py_FindMethod(py_handler_methods, self, name);
}

static int setattr_py_dist_constraint(PyObject *self, char *name, PyObject *value)
{
    Py_Dist_constraint py_dist_constraint = (Py_Dist_constraint) self;
    Dist_constraint dist_constraint = py_dist_constraint->dist_constraint;
    int vi;
    float vf;

    if (equal_strings(name, "atom0"))
    {
	vi = (int) PyInt_AsLong(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have int value");

	dist_constraint->atom0 = vi;
    }
    else if (equal_strings(name, "atom1"))
    {
	vi = (int) PyInt_AsLong(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have int value");

	dist_constraint->atom1 = vi;
    }
    else if (equal_strings(name, "dist_lower"))
    {
	vf = (float) PyFloat_AsDouble(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have float value");

	dist_constraint->dist_lower = vf;
    }
    else if (equal_strings(name, "dist_upper"))
    {
	vf = (float) PyFloat_AsDouble(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have float value");

	dist_constraint->dist_upper = vf;
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
static PySequenceMethods Dist_constraint_sequence_methods =
{
    Dist_constraint_length,
    Dist_constraint_concat,
    Dist_constraint_repeat,
    Dist_constraint_item,
    Dist_constraint_slice,
    Dist_constraint_ass_item,
    Dist_constraint_ass_slice
};

static PySequenceMethods Dist_constraint_sequence_methods =
{
    Dist_constraint_length,
    0,
    0,
    Dist_constraint_item,
    0,
    Dist_constraint_ass_item,
    0
};
*/

static PyTypeObject Dist_constraint_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "DistConstraint", /* name */
    sizeof(struct Py_Dist_constraint), /* basicsize */
    0, /* itemsize */
    delete_py_dist_constraint, /* destructor */
    0, /* printfunc */
    getattr_py_dist_constraint, /* getattr */
    setattr_py_dist_constraint, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Dist_constraint_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Dist_constraint(PyObject *self, PyObject *args)
{
    int atom0, atom1;
    float dist_lower, dist_upper;

    if (!PyArg_ParseTuple(args, "iiff", &atom0, &atom1,
					&dist_lower, &dist_upper))
        RETURN_OBJ_ERROR("must have four arguments:  atom0, atom1, dist_lower, dist_upper");

    return (PyObject *) new_py_dist_constraint(atom0, atom1, dist_lower, dist_upper);
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Dist_constraint_type_methods[] =
{
    { "DistConstraint",	(PyCFunction) init_Py_Dist_constraint,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import DistConstraint" in 
* a Python program. The function is usually called "initDist_constraint": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initDistConstraint(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Dist_constraint_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("DistConstraint", Dist_constraint_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("DistConstraint.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module DistConstraint");
}

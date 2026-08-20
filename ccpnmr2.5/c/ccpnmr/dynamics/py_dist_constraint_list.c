
/*
======================COPYRIGHT/LICENSE START==========================

py_dist_constraint_list.c: Part of the CcpNmr Analysis program

Copyright (C) 2009 Wayne Boucher and Tim Stevens (University of Cambridge)

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
#include "py_dist_constraint_list.h"

#include "python_util.h"

#define  NALLOC  100

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Dist_constraint_list_type;

Bool is_py_dist_constraint_list(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Dist_constraint_list_type
    return (obj->ob_type == &Dist_constraint_list_type);
*/
    return valid_py_object(obj, &Dist_constraint_list_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

Dist_constraint_list get_dist_constraint_list(Py_Dist_constraint_list py_dist_constraint_list)
{
    int i;
    Dist_constraint_list dist_constraint_list = new_dist_constraint_list(0, NULL, NULL, NULL, NULL, NULL);
    Dist_constraint dist_constraint;
    Line error_msg;
 
    if (!dist_constraint_list)
        return NULL;
 
    for (i = 0; i < py_dist_constraint_list->ndist_constraints; i++)
    {
        dist_constraint = py_dist_constraint_list->py_dist_constraints[i]->dist_constraint;
        if (append_dist_constraint_list(dist_constraint_list,
				dist_constraint, error_msg) == CCPN_ERROR)
        {
            clear_dist_constraint_list(dist_constraint_list);
            return NULL;
        }
    }
 
    return dist_constraint_list;
}

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static CcpnStatus add_py_dist_constraint(Py_Dist_constraint_list py_dist_constraint_list,
			Py_Dist_constraint py_dist_constraint, CcpnString error_msg)
{
    int n;
    int natoms = py_dist_constraint_list->ndist_constraints;
    int nalloc = py_dist_constraint_list->ndist_constraints_alloc;
 
    if (natoms >= nalloc)
    {
        sprintf(error_msg, "allocating py atom coords memory");
        if (natoms == 0)
        {
            n = NALLOC;
            MALLOC(py_dist_constraint_list->py_dist_constraints, Py_Dist_constraint, n);
        }
        else
        {
            n = nalloc + NALLOC;
            REALLOC(py_dist_constraint_list->py_dist_constraints, Py_Dist_constraint, n);
        }
 
        py_dist_constraint_list->ndist_constraints_alloc = n;
    }

    py_dist_constraint_list->py_dist_constraints[natoms] = py_dist_constraint;
    py_dist_constraint_list->ndist_constraints++;

    return CCPN_OK;
}

static PyObject *add(PyObject *self, PyObject *args)
{
    Py_Dist_constraint_list py_dist_constraint_list = (Py_Dist_constraint_list) self;
    Py_Dist_constraint py_dist_constraint;
    PyObject *atoms_obj0, *atoms_obj1;
    int natom_pairs, n, *atoms0, *atoms1;
    float dist_lower, dist_upper;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "OOff", &atoms_obj0, &atoms_obj1,
                                        &dist_lower, &dist_upper))
        RETURN_OBJ_ERROR("must have four arguments:  atom0, atom1, dist_lower, dist_upper"); 

    if (get_python_int_alloc_array(atoms_obj0, &natom_pairs, &atoms0, error_msg) == CCPN_ERROR) 
    {
        RETURN_OBJ_ERROR("first argument must be list of atoms0");
    }

    if (get_python_int_alloc_array(atoms_obj1, &n, &atoms1, error_msg) == CCPN_ERROR)
    {
        FREE(atoms0, int);
        RETURN_OBJ_ERROR("second argument must be list of atoms1");
    }

    if (n != natom_pairs)
    {
        FREE(atoms0, int);
        FREE(atoms1, int);
        sprintf(error_msg, "len(atoms0) = %d != %d = len(atoms1)", natom_pairs, n);
        RETURN_OBJ_ERROR(error_msg);
    }

    py_dist_constraint = new_py_dist_constraint(natom_pairs, atoms0, atoms1,
						dist_lower, dist_upper);

    if (!py_dist_constraint)
	RETURN_OBJ_ERROR("allocating Py_Dist_constraint object");

    if (add_py_dist_constraint(py_dist_constraint_list, py_dist_constraint,
							error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);
 
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *append(PyObject *self, PyObject *args)
{
    Py_Dist_constraint_list obj = (Py_Dist_constraint_list) self;
    PyObject *dist_constraint_obj;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &dist_constraint_obj))
        RETURN_OBJ_ERROR("must have one argument: distConstraint");

    if (!is_py_dist_constraint(dist_constraint_obj))
        RETURN_OBJ_ERROR("argument must be DistConstraint object");

    if (add_py_dist_constraint(obj, (Py_Dist_constraint) dist_constraint_obj, error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);
 
    Py_INCREF(dist_constraint_obj);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "add",		add,		METH_VARARGS },
    { "append",		append,		METH_VARARGS },
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static Py_Dist_constraint_list new_py_dist_constraint_list(void)
{
    Py_Dist_constraint_list py_dist_constraint_list;

    PY_MALLOC(py_dist_constraint_list, struct Py_Dist_constraint_list, &Dist_constraint_list_type);

    if (!py_dist_constraint_list)
	RETURN_OBJ_ERROR("allocating Py_Dist_constraint_list object");

    py_dist_constraint_list->py_dist_constraints = NULL;
    py_dist_constraint_list->ndist_constraints = 0;
    py_dist_constraint_list->ndist_constraints_alloc = 0;

    return py_dist_constraint_list;
}

static void delete_py_dist_constraint_list(PyObject *self)
{
    int i;
    Py_Dist_constraint_list py_dist_constraint_list = (Py_Dist_constraint_list) self;

/*
    printf("in delete_py_dist_constraint_list\n");
*/

    for (i = 0; i < py_dist_constraint_list->ndist_constraints; i++)
	Py_DECREF((PyObject*) py_dist_constraint_list->py_dist_constraints[i]);

    FREE(py_dist_constraint_list->py_dist_constraints, Py_Dist_constraint);

    PY_FREE(self);
}

/*
static int print_py_dist_constraint_list(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static int Dist_constraint_list_length(PyObject *self)
{
    Py_Dist_constraint_list py_dist_constraint_list = (Py_Dist_constraint_list) self;

    return py_dist_constraint_list->ndist_constraints;
}

static PyObject *Dist_constraint_list_item(PyObject *self, int i)
{
    Py_Dist_constraint_list py_dist_constraint_list = (Py_Dist_constraint_list) self;
    Py_Dist_constraint py_dist_constraint;

    if (i < 0 || i >= py_dist_constraint_list->ndist_constraints)
    {
	PyErr_SetString(PyExc_IndexError, "array index out of range");
	return NULL;
    }

    py_dist_constraint = py_dist_constraint_list->py_dist_constraints[i];
    Py_INCREF(py_dist_constraint);
 
    return (PyObject *) py_dist_constraint;
}

static PyObject *getattr_py_dist_constraint_list(PyObject *self, char *name)
{
/*
    Py_Dist_constraint_list obj = (Py_Dist_constraint_list) self;
    Dist_constraint_list dist_constraint_list = obj->dist_constraint_list;

    if (equal_strings(name, "mass"))
	return Py_BuildValue("f", dist_constraint_list->mass);
    else if (equal_strings(name, "x"))
	return Py_BuildValue("f", dist_constraint_list->x);
    else if (equal_strings(name, "y"))
	return Py_BuildValue("f", dist_constraint_list->y);
    else if (equal_strings(name, "z"))
	return Py_BuildValue("f", dist_constraint_list->z);
    else
*/
	return Py_FindMethod(py_handler_methods, self, name);
}

static int setattr_py_dist_constraint_list(PyObject *self, char *name, PyObject *value)
{
/*
    Py_Dist_constraint_list obj = (Py_Dist_constraint_list) self;
    Dist_constraint_list dist_constraint_list = obj->dist_constraint_list;
    float v = (float) PyFloat_AsDouble(value);

    if (PyErr_Occurred())
	RETURN_INT_ERROR("must have float value");

    if (equal_strings(name, "mass"))
	dist_constraint_list->mass = v;
    else if (equal_strings(name, "x"))
	dist_constraint_list->x = v;
    else if (equal_strings(name, "y"))
	dist_constraint_list->y = v;
    else if (equal_strings(name, "z"))
	dist_constraint_list->z = v;
    else
*/
	RETURN_INT_ERROR("unknown attribute name");

/*
    return 0;
*/
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Dist_constraint_list_sequence_methods =
{
    Dist_constraint_list_length,
    Dist_constraint_list_concat,
    Dist_constraint_list_repeat,
    Dist_constraint_list_item,
    Dist_constraint_list_slice,
    Dist_constraint_list_ass_item,
    Dist_constraint_list_ass_slice
};
*/

static PySequenceMethods Dist_constraint_list_sequence_methods =
{
    Dist_constraint_list_length,
    0,
    0,
    Dist_constraint_list_item,
    0,
    0,
    0
};

static PyTypeObject Dist_constraint_list_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "DistConstraintList", /* name */
    sizeof(struct Py_Dist_constraint_list), /* basicsize */
    0, /* itemsize */
    delete_py_dist_constraint_list, /* destructor */
    0, /* printfunc */
    getattr_py_dist_constraint_list, /* getattr */
    setattr_py_dist_constraint_list, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    &Dist_constraint_list_sequence_methods /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Dist_constraint_list(PyObject *self, PyObject *args)
{
    PyObject *obj;

    if (!PyArg_ParseTuple(args, ""))
        RETURN_OBJ_ERROR("must have no arguments");

    obj = (PyObject *) new_py_dist_constraint_list();

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


static struct PyMethodDef Dist_constraint_list_type_methods[] =
{
    { "DistConstraintList",	(PyCFunction) init_Py_Dist_constraint_list,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import DistConstraintList" in 
* a Python program. The function is usually called "initDist_constraint_list": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initDyDistConstraintList(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Dist_constraint_list_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("DyDistConstraintList", Dist_constraint_list_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("DyDistConstraintList.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module DyDistConstraintList");
}

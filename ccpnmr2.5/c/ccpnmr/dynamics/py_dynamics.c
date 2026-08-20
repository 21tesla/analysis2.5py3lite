
/*
======================COPYRIGHT/LICENSE START==========================

py_dynamics.c: Part of the CcpNmr Analysis program

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
#include "py_dynamics.h"

#include "py_atom_coord_list.h"
#include "py_dist_constraint_list.h"
#include "py_dist_force.h"

#include "python_util.h"

#include "utility.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Dynamics_type;

Bool is_py_dynamics(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Dynamics_type
    return (obj->ob_type == &Dynamics_type);
*/
    return valid_py_object(obj, &Dynamics_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *run(PyObject *self, PyObject *args)
{
    Py_Dynamics py_dynamics = (Py_Dynamics) self;
    PyObject *atom_coord_list_obj, *noe_list_obj, *noe_force_obj;
    Py_Atom_coord_list py_atom_coord_list;
    Py_Dist_constraint_list py_noe_list;
    Py_Dist_force py_noe_force;
    Atom_coord_list atom_coord_list;
    Dist_constraint_list noe_list;
    Line error_msg;
    CcpnStatus status;
 
    if (!PyArg_ParseTuple(args, "OOO", &atom_coord_list_obj,
				&noe_list_obj, &noe_force_obj))
        RETURN_OBJ_ERROR("must have three arguments: atomCoordList, noeList, noeForce");
 
    if (!is_py_atom_coord_list(atom_coord_list_obj))
        RETURN_OBJ_ERROR("first argument must be AtomCoordList object");
 
    if (!is_py_dist_constraint_list(noe_list_obj))
        RETURN_OBJ_ERROR("second argument must be DistConstraintList object");
 
    if (!is_py_dist_force(noe_force_obj))
        RETURN_OBJ_ERROR("third argument must be DistForce object");
 
    py_atom_coord_list = (Py_Atom_coord_list) atom_coord_list_obj;
    py_noe_list = (Py_Dist_constraint_list) noe_list_obj;
    py_noe_force = (Py_Dist_force) noe_force_obj;

    atom_coord_list = get_atom_coord_list(py_atom_coord_list);
    if (!atom_coord_list)
	RETURN_OBJ_ERROR("allocating atom coord list");

    noe_list = get_dist_constraint_list(py_noe_list);
    if (!noe_list)
    {
	clear_atom_coord_list(atom_coord_list);
	RETURN_OBJ_ERROR("allocating noe list");
    }

    status = run_dynamics(py_dynamics->dynamics, atom_coord_list, noe_list,
				py_noe_force->dist_force, error_msg);

    clear_atom_coord_list(atom_coord_list);
    clear_dist_constraint_list(noe_list);

    if (status == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "run",		run,		METH_VARARGS },
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static Py_Dynamics new_py_dynamics(float rp_force_const, float beta,
		float rmin, float drzap, float tref, float tau,
		float elapsed_time, int nsteps, int nprint)
{
    Py_Dynamics py_dynamics;
    Dynamics dynamics;

    dynamics = new_dynamics(rp_force_const, beta, rmin, drzap, tref,
				tau, elapsed_time, nsteps, nprint);

    if (!dynamics)
	 RETURN_OBJ_ERROR("allocating Dynamics object");

    PY_MALLOC(py_dynamics, struct Py_Dynamics, &Dynamics_type);

    if (!py_dynamics)
    {
	delete_dynamics(dynamics);

	RETURN_OBJ_ERROR("allocating Py_Dynamics object");
    }

    py_dynamics->dynamics = dynamics;

    return py_dynamics;
}

static void delete_py_dynamics(PyObject *self)
{
    Py_Dynamics py_dynamics = (Py_Dynamics) self;
    Dynamics dynamics = py_dynamics->dynamics;

/*
    printf("in delete_py_dynamics\n");
*/

    delete_dynamics(dynamics);

    PY_FREE(self);
}

static int print_py_dynamics(PyObject *self, FILE *fp, int flags)
{
    Py_Dynamics py_dynamics = (Py_Dynamics) self;
    Dynamics dynamics = py_dynamics->dynamics;

    fprintf(fp, "<rp_force_const=%3.2e, beta=%3.2e, rmin=%3.2e, drzap=%3.2e, tref=%3.2e, tau=%3.2e, elapsed_time=%3.2e, nsteps=%d, nprint=%d>",
		dynamics->rp_force_const, dynamics->beta, dynamics->rmin,
		dynamics->drzap, dynamics->tref, dynamics->tau,
		dynamics->elapsed_time, dynamics->nsteps, dynamics->nprint);

    return 0;
}

static PyObject *getattr_py_dynamics(PyObject *self, char *name)
{
    Py_Dynamics py_dynamics = (Py_Dynamics) self;
    Dynamics dynamics = py_dynamics->dynamics;

    if (equal_strings(name, "rp_force_const"))
	return Py_BuildValue("f", dynamics->rp_force_const);
    else if (equal_strings(name, "beta"))
	return Py_BuildValue("f", dynamics->beta);
    else if (equal_strings(name, "rmin"))
	return Py_BuildValue("f", dynamics->rmin);
    else if (equal_strings(name, "drzap"))
	return Py_BuildValue("f", dynamics->drzap);
    else if (equal_strings(name, "tref"))
	return Py_BuildValue("f", dynamics->tref);
    else if (equal_strings(name, "tau"))
	return Py_BuildValue("f", dynamics->tau);
    else if (equal_strings(name, "elapsed_time"))
	return Py_BuildValue("f", dynamics->elapsed_time);
    else if (equal_strings(name, "nsteps"))
	return Py_BuildValue("i", dynamics->nsteps);
    else if (equal_strings(name, "nprint"))
	return Py_BuildValue("i", dynamics->nprint);
    else
	return Py_FindMethod(py_handler_methods, self, name);
}

static int setattr_py_dynamics(PyObject *self, char *name, PyObject *value)
{
    Py_Dynamics py_dynamics = (Py_Dynamics) self;
    Dynamics dynamics = py_dynamics->dynamics;
    int vi;
    float vf;

    if (equal_strings(name, "rp_force_const"))
    {
	vf = (float) PyFloat_AsDouble(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have float value");

	dynamics->rp_force_const = vf;
    }
    else if (equal_strings(name, "beta"))
    {
	vf = (float) PyFloat_AsDouble(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have float value");

	dynamics->beta = vf;
    }
    else if (equal_strings(name, "rmin"))
    {
	vf = (float) PyFloat_AsDouble(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have float value");

	dynamics->rmin = vf;
    }
    else if (equal_strings(name, "drzap"))
    {
	vf = (float) PyFloat_AsDouble(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have float value");

	dynamics->drzap = vf;
    }
    else if (equal_strings(name, "tref"))
    {
	vf = (float) PyFloat_AsDouble(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have float value");

	dynamics->tref = vf;
    }
    else if (equal_strings(name, "tau"))
    {
	vf = (float) PyFloat_AsDouble(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have float value");

	dynamics->tau = vf;
    }
    else if (equal_strings(name, "elapsed_time"))
    {
	vf = (float) PyFloat_AsDouble(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have float value");

	dynamics->elapsed_time = vf;
    }
    else if (equal_strings(name, "nsteps"))
    {
	vi = (int) PyInt_AsLong(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have int value");

	dynamics->nsteps = vi;
    }
    else if (equal_strings(name, "nprint"))
    {
	vi = (int) PyInt_AsLong(value);
	if (PyErr_Occurred())
	    RETURN_INT_ERROR("must have int value");

	dynamics->nprint = vi;
    }
    else
	RETURN_INT_ERROR("unknown attribute name");

    return 0;
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Dynamics_sequence_methods =
{
    Dynamics_length,
    Dynamics_concat,
    Dynamics_repeat,
    Dynamics_item,
    Dynamics_slice,
    Dynamics_ass_item,
    Dynamics_ass_slice
};

static PySequenceMethods Dynamics_sequence_methods =
{
    Dynamics_length,
    0,
    0,
    Dynamics_item,
    0,
    Dynamics_ass_item,
    0
};
*/

static PyTypeObject Dynamics_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "Dynamics", /* name */
    sizeof(struct Py_Dynamics), /* basicsize */
    0, /* itemsize */
    delete_py_dynamics, /* destructor */
    print_py_dynamics, /* printfunc */
    getattr_py_dynamics, /* getattr */
    setattr_py_dynamics, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Dynamics_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Dynamics(PyObject *self, PyObject *args, PyObject *keywds)
{
    float rp_force_const = 1, beta = 10, rmin = 2.25, drzap = 2;
    float tref = 1000, tau = 0.001, elapsed_time = 0;
    int nsteps = 1000, nprint = 3000;
    static char *kwlist[] = { "rp_force_const", "beta", "rmin",
	"drzap", "tref", "tau", "elapsed_time", "nsteps", "nprint", NULL };

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "|fffffffii", kwlist,
		&rp_force_const, &beta, &rmin, &drzap, &tref, &tau,
		&elapsed_time, &nsteps, &nprint))
        RETURN_OBJ_ERROR("need arguments: [ rp_force_const, beta, rmin, drzap, tref, tau, elapsed_time, nsteps, nprint ]");

    return (PyObject *) new_py_dynamics(rp_force_const, beta, rmin,
			drzap, tref, tau, elapsed_time, nsteps, nprint);
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Dynamics_type_methods[] =
{
    { "Dynamics",	(PyCFunction) init_Py_Dynamics,	METH_VARARGS|METH_KEYWORDS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import Dynamics" in 
* a Python program. The function is usually called "initDynamics": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initDyDynamics(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Dynamics_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("DyDynamics", Dynamics_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("DyDynamics.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module DyDynamics");
}

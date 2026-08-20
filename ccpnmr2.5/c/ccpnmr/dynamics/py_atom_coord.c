
/*
======================COPYRIGHT/LICENSE START==========================

py_atom_coord.c: Part of the CcpNmr Analysis program

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
#include "py_atom_coord.h"

#include "python_util.h"

#include "utility.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Atom_coord_type;

Bool is_py_atom_coord(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Atom_coord_type
    return (obj->ob_type == &Atom_coord_type);
*/
    return valid_py_object(obj, &Atom_coord_type);
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

Py_Atom_coord new_py_atom_coord(float mass, float x, float y, float z, int isFixed)
{
    Py_Atom_coord py_atom_coord;
    Atom_coord atom_coord;

    atom_coord = new_atom_coord(mass, x, y, z, isFixed);

    if (!atom_coord)
	 RETURN_OBJ_ERROR("allocating Atom_coord object");

    PY_MALLOC(py_atom_coord, struct Py_Atom_coord, &Atom_coord_type);

    if (!py_atom_coord)
    {
	delete_atom_coord(atom_coord);

	RETURN_OBJ_ERROR("allocating Py_Atom_coord object");
    }

    py_atom_coord->atom_coord = atom_coord;

    return py_atom_coord;
}

static void delete_py_atom_coord(PyObject *self)
{
    Py_Atom_coord py_atom_coord = (Py_Atom_coord) self;
    Atom_coord atom_coord = py_atom_coord->atom_coord;

/*
    printf("in delete_py_atom_coord\n");
*/

    delete_atom_coord(atom_coord);

    PY_FREE(self);
}

static int print_py_atom_coord(PyObject *self, FILE *fp, int flags)
{
    Py_Atom_coord py_atom_coord = (Py_Atom_coord) self;
    Atom_coord atom_coord = py_atom_coord->atom_coord;
 
    fprintf(fp, "<mass=%3.2e, x=%3.2e, y=%3.2e, z=%3.2e fixed=%d>", atom_coord->mass,
			atom_coord->x, atom_coord->y, atom_coord->z, atom_coord->isFixed);

    return 0;
}

static PyObject *getattr_py_atom_coord(PyObject *self, char *name)
{
    Py_Atom_coord py_atom_coord = (Py_Atom_coord) self;
    Atom_coord atom_coord = py_atom_coord->atom_coord;

    if (equal_strings(name, "mass"))
	return Py_BuildValue("f", atom_coord->mass);
    else if (equal_strings(name, "x"))
	return Py_BuildValue("f", atom_coord->x);
    else if (equal_strings(name, "y"))
	return Py_BuildValue("f", atom_coord->y);
    else if (equal_strings(name, "z"))
	return Py_BuildValue("f", atom_coord->z);
    else if (equal_strings(name, "isFixed"))
	return Py_BuildValue("i", atom_coord->isFixed);
    else
	return Py_FindMethod(py_handler_methods, self, name);
}

static int setattr_py_atom_coord(PyObject *self, char *name, PyObject *value)
{
    Py_Atom_coord py_atom_coord = (Py_Atom_coord) self;
    Atom_coord atom_coord = py_atom_coord->atom_coord;
    float v = (float) PyFloat_AsDouble(value);

    if (PyErr_Occurred())
	RETURN_INT_ERROR("must have float value");

    if (equal_strings(name, "mass"))
	atom_coord->mass = v;
    else if (equal_strings(name, "x"))
	atom_coord->x = v;
    else if (equal_strings(name, "y"))
	atom_coord->y = v;
    else if (equal_strings(name, "z"))
	atom_coord->z = v;
    else if (equal_strings(name, "isFixed"))
	atom_coord->isFixed = v;
    else
	RETURN_INT_ERROR("unknown attribute name");

    return 0;
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Atom_coord_sequence_methods =
{
    Atom_coord_length,
    Atom_coord_concat,
    Atom_coord_repeat,
    Atom_coord_item,
    Atom_coord_slice,
    Atom_coord_ass_item,
    Atom_coord_ass_slice
};

static PySequenceMethods Atom_coord_sequence_methods =
{
    Atom_coord_length,
    0,
    0,
    Atom_coord_item,
    0,
    Atom_coord_ass_item,
    0
};
*/

static PyTypeObject Atom_coord_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "AtomCoord", /* name */
    sizeof(struct Py_Atom_coord), /* basicsize */
    0, /* itemsize */
    delete_py_atom_coord, /* destructor */
    print_py_atom_coord, /* printfunc */
    getattr_py_atom_coord, /* getattr */
    setattr_py_atom_coord, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Atom_coord_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Atom_coord(PyObject *self, PyObject *args)
{
    float mass, x, y, z;
    int isFixed;

    if (!PyArg_ParseTuple(args, "ffffi", &mass, &x, &y, &z, &isFixed))
        RETURN_OBJ_ERROR("must have five arguments: mass, x, y, z, isFixed");

    return (PyObject *) new_py_atom_coord(mass, x, y, z, isFixed);
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Atom_coord_type_methods[] =
{
    { "AtomCoord",	(PyCFunction) init_Py_Atom_coord,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import AtomCoord" in 
* a Python program. The function is usually called "initAtom_coord": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initDyAtomCoord(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Atom_coord_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("DyAtomCoord", Atom_coord_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("DyAtomCoord.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module DyAtomCoord");
}

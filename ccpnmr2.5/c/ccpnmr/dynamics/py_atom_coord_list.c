
/*
======================COPYRIGHT/LICENSE START==========================

py_atom_coord_list.c: Part of the CcpNmr Analysis program

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
#include "py_atom_coord_list.h"

#include "python_util.h"

#define  NALLOC  100

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Atom_coord_list_type;

Bool is_py_atom_coord_list(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Atom_coord_list_type
    return (obj->ob_type == &Atom_coord_list_type);
*/
    return valid_py_object(obj, &Atom_coord_list_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

Atom_coord_list get_atom_coord_list(Py_Atom_coord_list py_atom_coord_list)
{
    int i;
    Atom_coord_list atom_coord_list = new_atom_coord_list(0, NULL, NULL, NULL, NULL, NULL);
    Atom_coord atom_coord;
    Line error_msg;

    if (!atom_coord_list)
	return NULL;

    for (i = 0; i < py_atom_coord_list->natom_coords; i++)
    {
	atom_coord = py_atom_coord_list->py_atom_coords[i]->atom_coord;
	if (append_atom_coord_list(atom_coord_list, atom_coord, error_msg) == CCPN_ERROR)
	{
	    clear_atom_coord_list(atom_coord_list);
	    return NULL;
	}
    }

    return atom_coord_list;
}

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static CcpnStatus add_py_atom_coord(Py_Atom_coord_list py_atom_coord_list,
			Py_Atom_coord py_atom_coord, CcpnString error_msg)
{
    int n;
    int natoms = py_atom_coord_list->natom_coords;
    int nalloc = py_atom_coord_list->natom_coords_alloc;
 
    if (natoms >= nalloc)
    {
        sprintf(error_msg, "allocating py atom coords memory");
        if (natoms == 0)
        {
            n = NALLOC;
            MALLOC(py_atom_coord_list->py_atom_coords, Py_Atom_coord, n);
        }
        else
        {
            n = nalloc + NALLOC;
            REALLOC(py_atom_coord_list->py_atom_coords, Py_Atom_coord, n);
        }
 
        py_atom_coord_list->natom_coords_alloc = n;
    }

    py_atom_coord_list->py_atom_coords[natoms] = py_atom_coord;
    py_atom_coord_list->natom_coords++;

    return CCPN_OK;
}

static PyObject *add(PyObject *self, PyObject *args)
{
    Py_Atom_coord_list py_atom_coord_list = (Py_Atom_coord_list) self;
    Py_Atom_coord py_atom_coord;
    float mass, x, y, z;
    int isFixed;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "ffffi", &mass, &x, &y, &z, &isFixed))
        RETURN_OBJ_ERROR("must have five arguments: mass, x, y, z, isFixed");

    py_atom_coord = new_py_atom_coord(mass, x, y, z, isFixed);
    if (!py_atom_coord)
	RETURN_OBJ_ERROR("allocating Py_Atom_coord object");

    if (add_py_atom_coord(py_atom_coord_list, py_atom_coord, error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);
 
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *append(PyObject *self, PyObject *args)
{
    Py_Atom_coord_list obj = (Py_Atom_coord_list) self;
    PyObject *atom_coord_obj;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &atom_coord_obj))
        RETURN_OBJ_ERROR("must have one argument: atomCoord");

    if (!is_py_atom_coord(atom_coord_obj))
        RETURN_OBJ_ERROR("argument must be AtomCoord object");

    if (add_py_atom_coord(obj, (Py_Atom_coord) atom_coord_obj, error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);
 
    Py_INCREF(atom_coord_obj);

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

static Py_Atom_coord_list new_py_atom_coord_list(void)
{
    Py_Atom_coord_list py_atom_coord_list;

    PY_MALLOC(py_atom_coord_list, struct Py_Atom_coord_list, &Atom_coord_list_type);

    if (!py_atom_coord_list)
	RETURN_OBJ_ERROR("allocating Py_Atom_coord_list object");

    py_atom_coord_list->py_atom_coords = NULL;
    py_atom_coord_list->natom_coords = 0;
    py_atom_coord_list->natom_coords_alloc = 0;

    return py_atom_coord_list;
}

static void delete_py_atom_coord_list(PyObject *self)
{
    int i;
    Py_Atom_coord_list py_atom_coord_list = (Py_Atom_coord_list) self;

/*
    printf("in delete_py_atom_coord_list\n");
*/

    for (i = 0; i < py_atom_coord_list->natom_coords; i++)
	Py_DECREF((PyObject*) py_atom_coord_list->py_atom_coords[i]);

    FREE(py_atom_coord_list->py_atom_coords, Py_Atom_coord);

    PY_FREE(self);
}

/*
static int print_py_atom_coord_list(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static int Atom_coord_list_length(PyObject *self)
{
    Py_Atom_coord_list py_atom_coord_list = (Py_Atom_coord_list) self;

    return py_atom_coord_list->natom_coords;
}

static PyObject *Atom_coord_list_item(PyObject *self, int i)
{
    Py_Atom_coord_list py_atom_coord_list = (Py_Atom_coord_list) self;
    Py_Atom_coord py_atom_coord;

    if (i < 0 || i >= py_atom_coord_list->natom_coords)
    {
	PyErr_SetString(PyExc_IndexError, "array index out of range");
	return NULL;
    }

    py_atom_coord = py_atom_coord_list->py_atom_coords[i];
    Py_INCREF(py_atom_coord);
 
    return (PyObject *) py_atom_coord;
}

static PyObject *getattr_py_atom_coord_list(PyObject *self, char *name)
{
/*
    Py_Atom_coord_list obj = (Py_Atom_coord_list) self;
    Atom_coord_list atom_coord_list = obj->atom_coord_list;

    if (equal_strings(name, "mass"))
	return Py_BuildValue("f", atom_coord_list->mass);
    else if (equal_strings(name, "x"))
	return Py_BuildValue("f", atom_coord_list->x);
    else if (equal_strings(name, "y"))
	return Py_BuildValue("f", atom_coord_list->y);
    else if (equal_strings(name, "z"))
	return Py_BuildValue("f", atom_coord_list->z);
    else
*/
	return Py_FindMethod(py_handler_methods, self, name);
}

static int setattr_py_atom_coord_list(PyObject *self, char *name, PyObject *value)
{
/*
    Py_Atom_coord_list obj = (Py_Atom_coord_list) self;
    Atom_coord_list atom_coord_list = obj->atom_coord_list;
    float v = (float) PyFloat_AsDouble(value);

    if (PyErr_Occurred())
	RETURN_INT_ERROR("must have float value");

    if (equal_strings(name, "mass"))
	atom_coord_list->mass = v;
    else if (equal_strings(name, "x"))
	atom_coord_list->x = v;
    else if (equal_strings(name, "y"))
	atom_coord_list->y = v;
    else if (equal_strings(name, "z"))
	atom_coord_list->z = v;
    else if (equal_strings(name, "isFixed"))
	atom_coord_list->isFixed = v;
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
static PySequenceMethods Atom_coord_list_sequence_methods =
{
    Atom_coord_list_length,
    Atom_coord_list_concat,
    Atom_coord_list_repeat,
    Atom_coord_list_item,
    Atom_coord_list_slice,
    Atom_coord_list_ass_item,
    Atom_coord_list_ass_slice
};
*/

static PySequenceMethods Atom_coord_list_sequence_methods =
{
    Atom_coord_list_length,
    0,
    0,
    Atom_coord_list_item,
    0,
    0,
    0
};

static PyTypeObject Atom_coord_list_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "AtomCoordList", /* name */
    sizeof(struct Py_Atom_coord_list), /* basicsize */
    0, /* itemsize */
    delete_py_atom_coord_list, /* destructor */
    0, /* printfunc */
    getattr_py_atom_coord_list, /* getattr */
    setattr_py_atom_coord_list, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    &Atom_coord_list_sequence_methods /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Atom_coord_list(PyObject *self, PyObject *args)
{
    PyObject *obj;

    if (!PyArg_ParseTuple(args, ""))
        RETURN_OBJ_ERROR("must have no arguments");

    obj = (PyObject *) new_py_atom_coord_list();

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


static struct PyMethodDef Atom_coord_list_type_methods[] =
{
    { "AtomCoordList",	(PyCFunction) init_Py_Atom_coord_list,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import AtomCoordList" in 
* a Python program. The function is usually called "initAtom_coord_list": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initDyAtomCoordList(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Atom_coord_list_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("DyAtomCoordList", Atom_coord_list_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("DyAtomCoordList.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module DyAtomCoordList");
}

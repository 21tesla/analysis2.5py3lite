
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

    py_atom_coord = (Py_Atom_coord) PyObject_New(struct Py_Atom_coord, &Atom_coord_type);

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

    Py_TYPE(self)->tp_free(self);
}

static PyObject *repr_py_atom_coord(PyObject *self)
{
    Py_Atom_coord py_atom_coord = (Py_Atom_coord) self;
    Atom_coord atom_coord = py_atom_coord->atom_coord;
 
    char buf[128];

    snprintf(buf, sizeof buf, "<mass=%3.2e, x=%3.2e, y=%3.2e, z=%3.2e fixed=%d>", atom_coord->mass,
			atom_coord->x, atom_coord->y, atom_coord->z, atom_coord->isFixed);

    return PyUnicode_FromString(buf);
}

static PyObject *getattr_py_atom_coord(PyObject *self, PyObject *attr_name)
{
    Py_Atom_coord py_atom_coord = (Py_Atom_coord) self;
    Atom_coord atom_coord = py_atom_coord->atom_coord;

    if (PyUnicode_Check(attr_name))
    {
        const char *name = PyUnicode_AsUTF8(attr_name);
        if (name)
        {
            if (strcmp(name, "mass") == 0)
                return Py_BuildValue("f", atom_coord->mass);
            else if (strcmp(name, "x") == 0)
                return Py_BuildValue("f", atom_coord->x);
            else if (strcmp(name, "y") == 0)
                return Py_BuildValue("f", atom_coord->y);
            else if (strcmp(name, "z") == 0)
                return Py_BuildValue("f", atom_coord->z);
            else if (strcmp(name, "isFixed") == 0)
                return Py_BuildValue("i", atom_coord->isFixed);
        }
    }
    return PyObject_GenericGetAttr(self, attr_name);
}

static int setattr_py_atom_coord(PyObject *self, PyObject *attr_name, PyObject *value)
{
    Py_Atom_coord py_atom_coord = (Py_Atom_coord) self;
    Atom_coord atom_coord = py_atom_coord->atom_coord;
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

    if (strcmp(name, "mass") == 0)
        atom_coord->mass = v;
    else if (strcmp(name, "x") == 0)
        atom_coord->x = v;
    else if (strcmp(name, "y") == 0)
        atom_coord->y = v;
    else if (strcmp(name, "z") == 0)
        atom_coord->z = v;
    else if (strcmp(name, "isFixed") == 0)
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
    PyVarObject_HEAD_INIT(NULL, 0)
    "AtomCoord",                         /* tp_name */
    sizeof(struct Py_Atom_coord),        /* tp_basicsize */
    0,                                   /* tp_itemsize */
    (destructor) delete_py_atom_coord,   /* tp_dealloc */
    0,                                   /* tp_vectorcall */
    0,                                   /* tp_getattr */
    0,                                   /* tp_setattr */
    0,                                   /* tp_as_async */
    (reprfunc) repr_py_atom_coord,       /* tp_repr */
    0,                                   /* tp_as_number */
    0,                                   /* tp_as_sequence */
    0,                                   /* tp_as_mapping */
    0,                                   /* tp_hash */
    0,                                   /* tp_call */
    0,                                   /* tp_str */
    (getattrofunc) getattr_py_atom_coord,    /* tp_getattro */
    (setattrofunc) setattr_py_atom_coord,    /* tp_setattro */
    0,                                   /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                  /* tp_flags */
    "AtomCoord -- cartesian atom coordinate", /* tp_doc */
    0,                                   /* tp_traverse */
    0,                                   /* tp_clear */
    0,                                   /* tp_richcompare */
    0,                                   /* tp_weaklistoffset */
    0,                                   /* tp_iter */
    0,                                   /* tp_iternext */
    py_handler_methods,                  /* tp_methods */
    0,                                   /* tp_members */
    0,                                   /* tp_getset */
    0,                                   /* tp_base */
    0,                                   /* tp_dict */
    0,                                   /* tp_descr_get */
    0,                                   /* tp_descr_set */
    0,                                   /* tp_dictoffset */
    0,                                   /* tp_init */
    0,                                   /* tp_alloc */
    0,                                   /* tp_new */
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

static struct PyModuleDef DyAtomCoord_module_def =
{
    PyModuleDef_HEAD_INIT,
    "DyAtomCoord",
    "CCPNMR DyAtomCoord module (Python 3 compatible)",
    -1,
    Atom_coord_type_methods
};

PyMODINIT_FUNC PyInit_DyAtomCoord(void)
{
    if (PyType_Ready(&Atom_coord_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&DyAtomCoord_module_def);
    if (!m)
        return NULL;

    ErrorObject = PyErr_NewException("DyAtomCoord.error", NULL, NULL);
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

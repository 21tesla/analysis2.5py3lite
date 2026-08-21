
/*
======================COPYRIGHT/LICENSE START==========================

py_bond.c: Part of the CcpNmr Analysis program

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
#include "py_bond.h"

#include "py_atom.h"
#include "python_util.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Bond_type;

Bool is_py_bond(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Bond_type
    return (obj->ob_type == &Bond_type);
*/
    return valid_py_object(obj, &Bond_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *setColor(PyObject *self, PyObject *args)
{
    int n;
    Py_Bond obj = (Py_Bond) self;
    Bond bond = obj->bond;
    PyObject *color_obj = NULL;
    float color[ATOM_NCOLORS];
    Line error_msg;

    if (!PyArg_ParseTuple(args, "|O", &color_obj))
        RETURN_OBJ_ERROR("must have zero arguments or one argument (color)");

    if (color_obj)
    {
        if ((get_python_float_array(color_obj, BOND_NCOLORS, &n, color,
                                                error_msg) == CCPN_ERROR)
            || (n != ATOM_NCOLORS))
        {
            sprintf(error_msg, "color must be list or tuple of size %d", BOND_NCOLORS);
            RETURN_OBJ_ERROR(error_msg);
        }

        set_color_bond(bond, color);
    }
    else
    {
        set_color_bond(bond, NULL);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setLineWidth(PyObject *self, PyObject *args)
{
    float line_width;
    Py_Bond obj = (Py_Bond) self;
    Bond bond = obj->bond;

    if (!PyArg_ParseTuple(args, "f", &line_width))
        RETURN_OBJ_ERROR("must have one argument: line_width");

    set_line_width_bond(bond, line_width);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setLineStyle(PyObject *self, PyObject *args)
{
    int line_style;
    Py_Bond obj = (Py_Bond) self;
    Bond bond = obj->bond;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "i", &line_style))
        RETURN_OBJ_ERROR("must have one argument: line_style");

    if ((line_style < 0) || (line_style >= NLINE_STYLES))
    {
        sprintf(error_msg, "line_style must be >= 0 and < %d", NLINE_STYLES);
        RETURN_OBJ_ERROR(error_msg);
    }

    set_line_style_bond(bond, line_style);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setAnnotation(PyObject *self, PyObject *args)
{
    Py_Bond obj = (Py_Bond) self;
    Bond bond = obj->bond;
    char *annotation;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "s", &annotation))
        RETURN_OBJ_ERROR("need argument: annotation");

    if (set_annotation_bond(bond, annotation) == CCPN_ERROR)
    {
        sprintf(error_msg, "setting annotation");
        RETURN_OBJ_ERROR(error_msg);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "setColor",               setColor,               METH_VARARGS },
    { "setLineWidth",           setLineWidth,           METH_VARARGS },
    { "setLineStyle",           setLineStyle,           METH_VARARGS },
    { "setAnnotation",          setAnnotation,          METH_VARARGS },
    { NULL,		        NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_bond(PyObject *atom1_obj, PyObject *atom2_obj,
                                                PyObject *color_obj)
{
    Bond bond;
    Atom atom1 = ((Py_Atom) atom1_obj)->atom;
    Atom atom2 = ((Py_Atom) atom2_obj)->atom;
    int n;
    float color[BOND_NCOLORS];
    Py_Bond obj;
    Line error_msg;
 
    if (color_obj)
    {
        if ((get_python_float_array(color_obj, BOND_NCOLORS, &n, color,
                                                        error_msg) == CCPN_ERROR)
            || (n != BOND_NCOLORS))
        {
            sprintf(error_msg, "color must be list or tuple of size %d", BOND_NCOLORS);
            RETURN_OBJ_ERROR(error_msg);
        }

        bond = new_bond(atom1, atom2, color);
    }
    else
    {
        bond = new_bond(atom1, atom2, NULL);
    }

    if (!bond)
	RETURN_OBJ_ERROR("allocating Bond object");

    obj = (Py_Bond) PyObject_New(struct Py_Bond, &Bond_type);

    if (!obj)
    {
	delete_bond(bond);

	RETURN_OBJ_ERROR("allocating Py_Bond object");
    }

    obj->bond = bond;
    obj->atom1_obj = atom1_obj;
    obj->atom2_obj = atom2_obj;

    Py_INCREF(atom1_obj);
    Py_INCREF(atom2_obj);

    return (PyObject *) obj;
}

static void delete_py_bond(PyObject *self)
{
    Py_Bond obj = (Py_Bond) self;
    Bond bond = obj->bond;

/*
    printf("in delete_py_bond\n");
*/

    delete_bond(bond);
    Py_DECREF(obj->atom1_obj);
    Py_DECREF(obj->atom2_obj);

    Py_TYPE(self)->tp_free(self);
}

/*
static int print_py_bond(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_bond(PyObject *self, PyObject *attr_name)
{
/*
    Bond *obj = (Bond *) self;
    Random_access *a = obj->py_handler;

    if (equal_strings(name, "par_file"))
	return Py_BuildValue("s", a->par_file);
    else if (equal_strings(name, "access_method"))
	return Py_BuildValue("s", access_method_name(a->access_method));
    else if (equal_strings(name, "data_format"))
	return get_Bond_format(a);
    else
*/
	return PyObject_GenericGetAttr(self, attr_name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Bond_sequence_methods =
{
    Bond_length,
    Bond_concat,
    Bond_repeat,
    Bond_item,
    Bond_slice,
    Bond_ass_item,
    Bond_ass_slice
};

static PySequenceMethods Bond_sequence_methods =
{
    Bond_length,
    0,
    0,
    Bond_item,
    0,
    Bond_ass_item,
    0
};
*/

static PyTypeObject Bond_type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    "StructBond",                        /* tp_name */
    sizeof(struct Py_Bond),              /* tp_basicsize */
    0,                                   /* tp_itemsize */
    (destructor) delete_py_bond,         /* tp_dealloc */
    0,                                   /* tp_vectorcall */
    0,                                   /* tp_getattr */
    0,                                   /* tp_setattr */
    0,                                   /* tp_as_async */
    0,                                   /* tp_repr (print commented out) */
    0,                                   /* tp_as_number */
    0,                                   /* tp_as_sequence */
    0,                                   /* tp_as_mapping */
    0,                                   /* tp_hash */
    0,                                   /* tp_call */
    0,                                   /* tp_str */
    (getattrofunc) getattr_py_bond,      /* tp_getattro */
    0,                                   /* tp_setattro */
    0,                                   /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                  /* tp_flags */
    "StructBond -- CCPNMR structure Bond", /* tp_doc */
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

static PyObject *init_Py_Bond(PyObject *self, PyObject *args)
{
    PyObject *atom1_obj, *atom2_obj, *color_obj = NULL, *obj;

    if (!PyArg_ParseTuple(args, "OO|O", &atom1_obj, &atom2_obj, &color_obj))
    {
        RETURN_OBJ_ERROR("must have arguments: atom1 atom2 [color]");
    }

    if (!is_py_atom(atom1_obj) || !is_py_atom(atom2_obj))
        RETURN_OBJ_ERROR("first two arguments must be atoms");

    obj = new_py_bond(atom1_obj, atom2_obj, color_obj);

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


static struct PyMethodDef Bond_type_methods[] =
{
    { "StructBond",	(PyCFunction) init_Py_Bond,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import StructBond" in 
* a Python program. The function is usually called "initStructBond": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

static struct PyModuleDef Bond_module_def =
{
    PyModuleDef_HEAD_INIT,
    "StructBond",
    "CCPNMR StructBond module (Python 3 compatible)",
    -1,
    Bond_type_methods
};

PyMODINIT_FUNC PyInit_StructBond(void)
{
    if (PyType_Ready(&Bond_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&Bond_module_def);
    if (!m)
        return NULL;

    ErrorObject = PyErr_NewException("StructBond.error", NULL, NULL);
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

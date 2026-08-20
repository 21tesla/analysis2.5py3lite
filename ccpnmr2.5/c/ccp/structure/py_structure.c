
/*
======================COPYRIGHT/LICENSE START==========================

py_structure.c: Part of the CcpNmr Analysis program

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
#include "py_structure.h"

#include "py_atom.h"
#include "py_bond.h"
#include "py_draw_handler.h"
#include "python_util.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Structure_type;

Bool is_py_structure(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Structure_type
    return (obj->ob_type == &Structure_type);
*/
    return valid_py_object(obj, &Structure_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *addAtom(PyObject *self, PyObject *args)
{
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;
    PyObject *atom_obj;
    Atom atom;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &atom_obj))
	RETURN_OBJ_ERROR("need argument: atom");

    if (!is_py_atom(atom_obj))
	RETURN_OBJ_ERROR("need argument: atom");

    atom = ((Py_Atom) atom_obj)->atom;

    if (add_atom_structure(structure, atom) == CCPN_ERROR)
    {
	sprintf(error_msg, "adding atom");
	RETURN_OBJ_ERROR(error_msg);
    }

    if (PyList_Append(obj->atom_list, atom_obj) != 0)
    {
	sprintf(error_msg, "adding atom_obj");
	RETURN_OBJ_ERROR(error_msg);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *removeAtom(PyObject *self, PyObject *args)
{
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;
    PyObject *atom_obj, *atom_list = obj->atom_list;
    PyMethodDef *methods;
    PyCFunction func;
    Atom atom;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &atom_obj))
	RETURN_OBJ_ERROR("need argument: atom");

    if (!is_py_atom(atom_obj))
	RETURN_OBJ_ERROR("need argument: atom");

    atom = ((Py_Atom) atom_obj)->atom;

    if (remove_atom_structure(structure, atom) == CCPN_ERROR)
    {
	sprintf(error_msg, "removing atom");
	RETURN_OBJ_ERROR(error_msg);
    }

    methods = atom_list->ob_type->tp_methods;
    func = PyCFunction_GetFunction(Py_FindMethod(methods, atom_list, "remove"));
    if ((func)(atom_list, atom_obj) == NULL)
    {
	sprintf(error_msg, "removing atom_obj");
	RETURN_OBJ_ERROR(error_msg);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *nearestAtom(PyObject *self, PyObject *args)
{
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;
    PyObject *aobj, *atom_obj, *atom_list = obj->atom_list;
    PySequenceMethods *seq;
    Atom atom;
    int i;
    float x, y, tol;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "fff", &x, &y, &tol))
	RETURN_OBJ_ERROR("need arguments: x, y, tolerance");

    atom = nearest_atom_structure(structure, x, y, tol);

    atom_obj = Py_None;
    if (atom != NULL)
    {
	seq = atom_list->ob_type->tp_as_sequence;
	for (i = 0; i < seq->sq_length(atom_list); i++)
	{
	    /* note that sq_item increments aobj reference count */
	    aobj = seq->sq_item(atom_list, i);
	    if (((Py_Atom) aobj)->atom == atom)
	    {
		atom_obj = aobj;
		break;
	    }

	    /* note that atom_obj is not decremented because of break */
	    Py_DECREF(aobj);
	}
    }

    return atom_obj;
}

static PyObject *addBond(PyObject *self, PyObject *args)
{
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;
    PyObject *bond_obj;
    Bond bond;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &bond_obj))
	RETURN_OBJ_ERROR("need argument: bond");

    if (!is_py_bond(bond_obj))
	RETURN_OBJ_ERROR("need argument: bond");

    bond = ((Py_Bond) bond_obj)->bond;

    if (add_bond_structure(structure, bond) == CCPN_ERROR)
    {
	sprintf(error_msg, "adding bond");
	RETURN_OBJ_ERROR(error_msg);
    }

    if (PyList_Append(obj->bond_list, bond_obj) != 0)
    {
	sprintf(error_msg, "adding bond_obj");
	RETURN_OBJ_ERROR(error_msg);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *removeBond(PyObject *self, PyObject *args)
{
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;
    PyObject *bond_obj, *bond_list = obj->bond_list;
    PyMethodDef *methods;
    PyCFunction func;
    Bond bond;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &bond_obj))
	RETURN_OBJ_ERROR("need argument: bond");

    if (!is_py_bond(bond_obj))
	RETURN_OBJ_ERROR("need argument: bond");

    bond = ((Py_Bond) bond_obj)->bond;

    if (remove_bond_structure(structure, bond) == CCPN_ERROR)
    {
	sprintf(error_msg, "removing bond");
	RETURN_OBJ_ERROR(error_msg);
    }

    methods = bond_list->ob_type->tp_methods;
    func = PyCFunction_GetFunction(Py_FindMethod(methods, bond_list, "remove"));
    if ((func)(bond_list, bond_obj) == NULL)
    {
	sprintf(error_msg, "removing bond_obj");
	RETURN_OBJ_ERROR(error_msg);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *nearestBond(PyObject *self, PyObject *args)
{
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;
    PyObject *bobj, *bond_obj, *bond_list = obj->bond_list;
    PySequenceMethods *seq;
    Bond bond;
    int i;
    float x, y, tol;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "fff", &x, &y, &tol))
	RETURN_OBJ_ERROR("need arguments: x, y, tolerance");

    bond = nearest_bond_structure(structure, x, y, tol);

    bond_obj = Py_None;
    if (bond != NULL)
    {
	seq = bond_list->ob_type->tp_as_sequence;
	for (i = 0; i < seq->sq_length(bond_list); i++)
	{
	    /* note that sq_item increments bobj reference count */
	    bobj = seq->sq_item(bond_list, i);
	    if (((Py_Bond) bobj)->bond == bond)
	    {
		bond_obj = bobj;
		break;
	    }

	    /* note that bond_obj is not decremented because of break */
	    Py_DECREF(bobj);
	}
    }

    return bond_obj;
}

static PyObject *draw(PyObject *self, PyObject *args)
{
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;
    PyObject *handler_obj;
    Py_draw_handler py_draw_handler;
    Drawing_funcs *drawing_funcs;
    Generic_ptr handler;

    if (!PyArg_ParseTuple(args, "O", &handler_obj))
	RETURN_OBJ_ERROR("need argument: handler");

    py_draw_handler = new_py_draw_handler(handler_obj);
    if (!py_draw_handler)
        RETURN_OBJ_ERROR("argument must be handler object");

    handler = py_draw_handler->handler;
    drawing_funcs = py_draw_handler->drawing_funcs;
    delete_py_draw_handler(py_draw_handler);

    draw_structure(structure, drawing_funcs, (Generic_ptr) handler);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *translate(PyObject *self, PyObject *args)
{
    int n;
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;
    PyObject *delta_obj;
    float delta[ATOM_NDIMS];
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &delta_obj))
        RETURN_OBJ_ERROR("must have argument: delta");

    if ((get_python_float_array(delta_obj, ATOM_NDIMS, &n, delta,
						error_msg) == CCPN_ERROR)
	|| (n != ATOM_NDIMS))
    {
	sprintf(error_msg, "delta must be list or tuple of size %d", ATOM_NDIMS);
	RETURN_OBJ_ERROR(error_msg);
    }

    translate_structure(structure, delta);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *rotate(PyObject *self, PyObject *args)
{
    int i, n, n1, n2;
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;
    PyObject *matrix_obj, *origin_obj = NULL;
    float **matrix;
    float origin[ATOM_NDIMS];
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O|O", &matrix_obj, &origin_obj))
        RETURN_OBJ_ERROR("must have one or two arguments: matrix [ origin ]");

    if ((get_python_float_alloc_matrix(matrix_obj, &n1, &n2, &matrix,
						error_msg) == CCPN_ERROR)
	|| (n1 != ATOM_NDIMS) || (n2 != ATOM_NDIMS))
    {
	sprintf(error_msg, "matrix must be %d x %d", ATOM_NDIMS, ATOM_NDIMS);
	RETURN_OBJ_ERROR(error_msg);
    }

    if (origin_obj)
    {
        if ((get_python_float_array(origin_obj, ATOM_NDIMS, &n, origin,
						error_msg) == CCPN_ERROR)
	    || (n != ATOM_NDIMS))
        {
	    sprintf(error_msg, "origin must be list or tuple of size %d", ATOM_NDIMS);
	    RETURN_OBJ_ERROR(error_msg);
        }
    }
    else
    {
        ZERO_VECTOR(origin, ATOM_NDIMS);
    }

    rotate_structure(structure, matrix, origin);

    for (i = 0; i < ATOM_NDIMS; i++)
        FREE(matrix[i], float);
    FREE(matrix, float *);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *zoom(PyObject *self, PyObject *args)
{
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;
    float delta;
    Line error_msg;
    sprintf(error_msg, "Camera IN %f", delta);

    if (!PyArg_ParseTuple(args, "f", &delta))
        RETURN_OBJ_ERROR("must have argument: delta");

    zoom_structure(structure, delta);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *moveToCenter(PyObject *self, PyObject *args)
{
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;

    if (!PyArg_ParseTuple(args, ""))
	RETURN_OBJ_ERROR("should have no arguments");


    move_to_center_structure(structure);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *focusOnAtom(PyObject *self, PyObject *args)
{
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;
    PyObject *atom_obj;
    Atom atom;

    if (!PyArg_ParseTuple(args, "O", &atom_obj))
	RETURN_OBJ_ERROR("need argument: atom");

    if (!is_py_atom(atom_obj))
	RETURN_OBJ_ERROR("need argument: atom");

    atom = ((Py_Atom) atom_obj)->atom;

    focus_on_atom_structure(structure, atom);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "addAtom",		addAtom,		METH_VARARGS },
    { "removeAtom",		removeAtom,		METH_VARARGS },
    { "nearestAtom",		nearestAtom,		METH_VARARGS },
    { "addBond",		addBond,		METH_VARARGS },
    { "removeBond",		removeBond,		METH_VARARGS },
    { "nearestBond",		nearestBond,		METH_VARARGS },
    { "draw",	                draw,   		METH_VARARGS },
    { "translate",		translate,		METH_VARARGS },
    { "rotate",		        rotate, 		METH_VARARGS },
    { "zoom",	        	zoom,		        METH_VARARGS },
    { "moveToCenter",		moveToCenter,		METH_VARARGS },
    { "focusOnAtom",		focusOnAtom,		METH_VARARGS },
    { NULL,		        NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_structure(void)
{
    Structure structure;
    Py_Structure obj;
 
    structure = new_structure();

    if (!structure)
	RETURN_OBJ_ERROR("allocating Structure object");

    PY_MALLOC(obj, struct Py_Structure, &Structure_type);

    obj->atom_list = PyList_New(0);
    obj->bond_list = PyList_New(0);

    if (!obj || !obj->atom_list || !obj->bond_list)
    {
	delete_structure(structure);

	RETURN_OBJ_ERROR("allocating Py_Structure object");
    }

    obj->structure = structure;

    return (PyObject *) obj;
}

static void delete_py_structure(PyObject *self)
{
    Py_Structure obj = (Py_Structure) self;
    Structure structure = obj->structure;

/*
    printf("in delete_py_structure\n");
*/

    delete_structure(structure);

    Py_DECREF(obj->atom_list);
    Py_DECREF(obj->bond_list);
    PY_FREE(self);
}

/*
static int print_py_structure(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_structure(PyObject *self, char *name)
{
/*
    Structure *obj = (Structure *) self;
    Random_access *a = obj->py_handler;

    if (equal_strings(name, "par_file"))
	return Py_BuildValue("s", a->par_file);
    else if (equal_strings(name, "access_method"))
	return Py_BuildValue("s", access_method_name(a->access_method));
    else if (equal_strings(name, "data_format"))
	return get_Structure_format(a);
    else
*/
	return Py_FindMethod(py_handler_methods, self, name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Structure_sequence_methods =
{
    Structure_length,
    Structure_concat,
    Structure_repeat,
    Structure_item,
    Structure_slice,
    Structure_ass_item,
    Structure_ass_slice
};

static PySequenceMethods Structure_sequence_methods =
{
    Structure_length,
    0,
    0,
    Structure_item,
    0,
    Structure_ass_item,
    0
};
*/

static PyTypeObject Structure_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "StructStructure", /* name */
    sizeof(struct Py_Structure), /* basicsize */
    0, /* itemsize */
    delete_py_structure, /* destructor */
    0, /* printfunc */
    getattr_py_structure, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Structure_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Structure(PyObject *self, PyObject *args)
{
    PyObject *obj;

    if (!PyArg_ParseTuple(args, ""))
        RETURN_OBJ_ERROR("must have no arguments");

    obj = new_py_structure();

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


static struct PyMethodDef Structure_type_methods[] =
{
    { "StructStructure",	(PyCFunction) init_Py_Structure,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import StructStructure" in 
* a Python program. The function is usually called "initStructStructure": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initStructStructure(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Structure_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("StructStructure", Structure_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("StructStructure.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module StructStructure");
}

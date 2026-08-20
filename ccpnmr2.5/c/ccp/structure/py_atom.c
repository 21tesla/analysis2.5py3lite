
/*
======================COPYRIGHT/LICENSE START==========================

py_atom.c: Part of the CcpNmr Analysis program

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
#include "py_atom.h"

#include "python_util.h"

#include "bond.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Atom_type;

Bool is_py_atom(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Atom_type
    return (obj->ob_type == &Atom_type);
*/
    return valid_py_object(obj, &Atom_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *setSize(PyObject *self, PyObject *args)
{
    float size;
    Py_Atom obj = (Py_Atom) self;
    Atom atom = obj->atom;

    if (!PyArg_ParseTuple(args, "f", &size))
	RETURN_OBJ_ERROR("need argument: size");

    set_size_atom(atom, size);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setSymbol(PyObject *self, PyObject *args)
{
    Py_Atom obj = (Py_Atom) self;
    Atom atom = obj->atom;
    char *symbol;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "s", &symbol))
	RETURN_OBJ_ERROR("need argument: symbol");

    if (set_symbol_atom(atom, symbol) == CCPN_ERROR)
    {
	sprintf(error_msg, "setting symbol");
	RETURN_OBJ_ERROR(error_msg);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setAnnotation(PyObject *self, PyObject *args)
{
    Py_Atom obj = (Py_Atom) self;
    Atom atom = obj->atom;
    char *annotation;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "s", &annotation))
	RETURN_OBJ_ERROR("need argument: annotation");

    if (set_annotation_atom(atom, annotation) == CCPN_ERROR)
    {
	sprintf(error_msg, "setting annotation");
	RETURN_OBJ_ERROR(error_msg);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setColor(PyObject *self, PyObject *args)
{
    int n;
    Py_Atom obj = (Py_Atom) self;
    Atom atom = obj->atom;
    PyObject *color_obj;
    float color[ATOM_NCOLORS];
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &color_obj))
        RETURN_OBJ_ERROR("must have argument: color");

    if ((get_python_float_array(color_obj, ATOM_NCOLORS, &n, color,
						error_msg) == CCPN_ERROR)
	|| (n != ATOM_NCOLORS))
    {
	sprintf(error_msg, "color must be list or tuple of size %d", ATOM_NCOLORS);
	RETURN_OBJ_ERROR(error_msg);
    }

    set_color_atom(atom, color);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *turnOn(PyObject *self, PyObject *args)
{
    Py_Atom obj = (Py_Atom) self;
    Atom atom = obj->atom;

    if (!PyArg_ParseTuple(args, ""))
	RETURN_OBJ_ERROR("should have no arguments");

    turn_on_atom(atom);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *turnOff(PyObject *self, PyObject *args)
{
    Py_Atom obj = (Py_Atom) self;
    Atom atom = obj->atom;

    if (!PyArg_ParseTuple(args, ""))
	RETURN_OBJ_ERROR("should have no arguments");

    turn_off_atom(atom);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setAnnotationColor(PyObject *self, PyObject *args)
{
    int n;
    Py_Atom obj = (Py_Atom) self;
    Atom atom = obj->atom;
    PyObject *color_obj = NULL;
    float color[ATOM_NCOLORS];
    Line error_msg;

    if (!PyArg_ParseTuple(args, "|O", &color_obj))
        RETURN_OBJ_ERROR("must have zero arguments or one argument (color)");

    if (color_obj)
    {
        if ((get_python_float_array(color_obj, ATOM_NCOLORS, &n, color,
						        error_msg) == CCPN_ERROR)
	        || (n != ATOM_NCOLORS))
        {
	    sprintf(error_msg, "color must be list or tuple of size %d", ATOM_NCOLORS);
	    RETURN_OBJ_ERROR(error_msg);
        }

        set_annotation_color_atom(atom, color);
    }
    else
    {
        set_annotation_color_atom(atom, NULL);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *getCoords(PyObject *self, PyObject *args)
{
    int n;
    Py_Atom obj = (Py_Atom) self;
    Atom atom = obj->atom;
    PyObject *coords_obj;
    Line error_msg;

    if (!PyArg_ParseTuple(args, ""))
        RETURN_OBJ_ERROR("must have zero arguments");

    coords_obj = get_python_float_list(ATOM_NDIMS, atom->x);

    return coords_obj;
}

static PyObject *setCoords(PyObject *self, PyObject *args)
{
    int n;
    Py_Atom obj = (Py_Atom) self;
    Atom atom = obj->atom;
    PyObject *coords_obj = NULL;
    float coords[ATOM_NDIMS];
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &coords_obj))
        RETURN_OBJ_ERROR("must have zero arguments or one argument (coords)");

    if ((get_python_float_array(coords_obj, ATOM_NCOLORS, &n, coords,
						        error_msg) == CCPN_ERROR)
	        || (n != ATOM_NDIMS))
    {
	sprintf(error_msg, "coords must be list or tuple of size %d", ATOM_NDIMS);
	RETURN_OBJ_ERROR(error_msg);
    }

    set_coords_atom(atom, coords);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "setSize",		setSize,		METH_VARARGS },
    { "setSymbol",		setSymbol,		METH_VARARGS },
    { "setAnnotation",		setAnnotation,		METH_VARARGS },
    { "setColor",	        setColor,		METH_VARARGS },
    { "turnOn",                 turnOn,                 METH_VARARGS },
    { "turnOff",                turnOff,                METH_VARARGS },
    { "setAnnotationColor",	setAnnotationColor,	METH_VARARGS },
    { "getCoords",		getCoords,		METH_VARARGS },
    { "setCoords",		setCoords,		METH_VARARGS },
    { NULL,		        NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_atom(float size, char *symbol, char *annotation,
		PyObject *x_obj, PyObject *color_obj)
{
    int n;
    float x[ATOM_NDIMS], color[ATOM_NCOLORS];
    Atom atom;
    Py_Atom obj;
    Line error_msg;
 
    if ((get_python_float_array(x_obj, ATOM_NDIMS, &n, x,
						error_msg) == CCPN_ERROR)
	|| (n != ATOM_NDIMS))
    {
	sprintf(error_msg, "x must be list or tuple of size %d", ATOM_NDIMS);
	RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_float_array(color_obj, ATOM_NDIMS, &n, color,
						error_msg) == CCPN_ERROR)
	|| (n != ATOM_NCOLORS))
    {
	sprintf(error_msg, "color must be list or tuple of size %d", ATOM_NCOLORS);
	RETURN_OBJ_ERROR(error_msg);
    }

    atom = new_atom(size, symbol, annotation, x, color);

    if (!atom)
	RETURN_OBJ_ERROR("allocating Atom object");

    PY_MALLOC(obj, struct Py_Atom, &Atom_type);

    if (!obj)
    {
	delete_atom(atom);

	RETURN_OBJ_ERROR("allocating Py_Atom object");
    }

    obj->atom = atom;

    return (PyObject *) obj;
}

static void delete_py_atom(PyObject *self)
{
    Py_Atom obj = (Py_Atom) self;
    Atom atom = obj->atom;

/*
    printf("in delete_py_atom\n");
*/

    delete_atom(atom);

    PY_FREE(self);
}

/*
static int print_py_atom(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_atom(PyObject *self, char *name)
{
/*
    Atom *obj = (Atom *) self;
    Random_access *a = obj->py_handler;

    if (equal_strings(name, "par_file"))
	return Py_BuildValue("s", a->par_file);
    else if (equal_strings(name, "access_method"))
	return Py_BuildValue("s", access_method_name(a->access_method));
    else if (equal_strings(name, "data_format"))
	return get_Atom_format(a);
    else
*/
	return Py_FindMethod(py_handler_methods, self, name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Atom_sequence_methods =
{
    Atom_length,
    Atom_concat,
    Atom_repeat,
    Atom_item,
    Atom_slice,
    Atom_ass_item,
    Atom_ass_slice
};

static PySequenceMethods Atom_sequence_methods =
{
    Atom_length,
    0,
    0,
    Atom_item,
    0,
    Atom_ass_item,
    0
};
*/

static PyTypeObject Atom_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "StructAtom", /* name */
    sizeof(struct Py_Atom), /* basicsize */
    0, /* itemsize */
    delete_py_atom, /* destructor */
    0, /* printfunc */
    getattr_py_atom, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Atom_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Atom(PyObject *self, PyObject *args)
{
    float size;
    char *symbol, *annotation;
    PyObject *x_obj, *color_obj, *obj;

    if (!PyArg_ParseTuple(args, "fssOO", &size, &symbol, &annotation,
		&x_obj, &color_obj))
        RETURN_OBJ_ERROR("must have arguments: size, symbol, annotation, x, color");

    obj = new_py_atom(size, symbol, annotation, x_obj, color_obj);

/*  this happens automatically?? (I think yes)
    FREE(symbol, char);
    FREE(annotation, char);
*/

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


static struct PyMethodDef Atom_type_methods[] =
{
    { "StructAtom",	(PyCFunction) init_Py_Atom,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import StructAtom" in 
* a Python program. The function is usually called "initStructAtom": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initStructAtom(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Atom_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("StructAtom", Atom_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("StructAtom.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module StructAtom");
}

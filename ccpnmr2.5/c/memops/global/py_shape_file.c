
/*
======================COPYRIGHT/LICENSE START==========================

py_shape_file.c: Part of the CcpNmr Analysis program

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
#include "py_shape_file.h"

#include "py_mem_cache.h"
#include "python_util.h"
#include "utility.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Shape_file_type;

Bool is_py_shape_file(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Shape_file_type
    return (obj->ob_type == &Shape_file_type);
*/
    return valid_py_object(obj, &Shape_file_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *setComponentAmplitude(PyObject *self, PyObject *args)
{
    Py_Shape_file obj = (Py_Shape_file) self;
    Shape_file shape_file = obj->shape_file;
    int comp;
    float amplitude;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "if", &comp, &amplitude))
	RETURN_OBJ_ERROR("need arguments: component, amplitude");

    if (set_amplitude_shape_component(shape_file, comp, amplitude, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setShapeData(PyObject *self, PyObject *args)
{
    Py_Shape_file obj = (Py_Shape_file) self;
    Shape_file shape_file = obj->shape_file;
    PyObject *values_obj;
    int comp, shape, size, offset;
    float *values;
    CcpnStatus status;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "iiiO", &comp, &shape, &offset, &values_obj))
	RETURN_OBJ_ERROR("need arguments: component, shape, offset, values");

    sprintf(error_msg, "values: ");
    if (get_python_float_alloc_array(values_obj, &size, &values,
				error_msg+strlen(error_msg)) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    status = set_data_shape_shape(shape_file, comp, shape, size, offset, values, error_msg);

    FREE(values, float);

    if (status == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "setComponentAmplitude",	setComponentAmplitude,	METH_VARARGS },
    { "setShapeData",		setShapeData,		METH_VARARGS },
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_shape_file(int ncomponents, PyObject *points_obj)
{
    int ndim, points[MAX_NDIM];
    Shape_file shape_file;
    Py_Shape_file obj;
    Line error_msg;
 
    sprintf(error_msg, "points: ");
    if (get_python_int_array(points_obj, MAX_NDIM, &ndim, points,
				error_msg+strlen(error_msg)) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    shape_file = new_shape_file(ndim, ncomponents, points);

    if (!shape_file)
	RETURN_OBJ_ERROR("allocating Shape_file object");

    PY_MALLOC(obj, struct Py_Shape_file, &Shape_file_type);

    if (!obj)
    {
	delete_shape_file(shape_file);

	RETURN_OBJ_ERROR("allocating Py_Shape_file object");
    }

    obj->shape_file = shape_file;

    return (PyObject *) obj;
}

static void delete_py_shape_file(PyObject *self)
{
    Py_Shape_file obj = (Py_Shape_file) self;
    Shape_file shape_file = obj->shape_file;

/*
    printf("in delete_py_shape_file\n");
*/

    delete_shape_file(shape_file);

    PY_FREE(self);
}

/*
static int print_py_shape_file(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_shape_file(PyObject *self, char *name)
{
    Py_Shape_file obj = (Py_Shape_file) self;
    Shape_file shape_file = obj->shape_file;

    if (equal_strings(name, "ndim") || equal_strings(name, "nshapes"))
	return Py_BuildValue("i", shape_file->ndim);
    else if (equal_strings(name, "ncomponents"))
	return Py_BuildValue("i", shape_file->ncomponents);
    else if (equal_strings(name, "points"))
	return get_python_int_list(shape_file->ndim, shape_file->points);
    else
	return Py_FindMethod(py_handler_methods, self, name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Shape_file_sequence_methods =
{
    Shape_file_length,
    Shape_file_concat,
    Shape_file_repeat,
    Shape_file_item,
    Shape_file_slice,
    Shape_file_ass_item,
    Shape_file_ass_slice
};

static PySequenceMethods Shape_file_sequence_methods =
{
    Shape_file_length,
    0,
    0,
    Shape_file_item,
    0,
    Shape_file_ass_item,
    0
};
*/

static PyTypeObject Shape_file_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "ShapeFile", /* name */
    sizeof(struct Py_Shape_file), /* basicsize */
    0, /* itemsize */
    delete_py_shape_file, /* destructor */
    0, /* printfunc */
    getattr_py_shape_file, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Shape_file_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Shape_file(PyObject *self, PyObject *args)
{
    int ncomponents;
    PyObject *points_obj, *obj;

    if (!PyArg_ParseTuple(args, "iO", &ncomponents, &points_obj))
        RETURN_OBJ_ERROR("must have arguments: ncomponents, points");

    obj = new_py_shape_file(ncomponents, points_obj);

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


static struct PyMethodDef Shape_file_type_methods[] =
{
    { "ShapeFile",	(PyCFunction) init_Py_Shape_file,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import ShapeFile" in 
* a Python program. The function is usually called "initShape_file": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initShapeFile(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Shape_file_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("ShapeFile", Shape_file_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("ShapeFile.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module ShapeFile");
}

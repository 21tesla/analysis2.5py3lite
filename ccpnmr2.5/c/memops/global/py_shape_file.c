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

#include "python_util.h"
#include "utility.h"

/* Locally-raised exception type */
static PyObject *ErrorObject;

static PyTypeObject Shape_file_type;

Bool is_py_shape_file(PyObject *obj)
{
    return valid_py_object(obj, &Shape_file_type);
}

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

    Py_RETURN_NONE;
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

    Py_RETURN_NONE;
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

    obj = (Py_Shape_file) PyObject_New(struct Py_Shape_file, &Shape_file_type);

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

    delete_shape_file(shape_file);

    Py_TYPE(self)->tp_free(self);
}

/* Python 3 tp_getattro — receives a unicode object, not a char* */
static PyObject *getattr_py_shape_file(PyObject *self, PyObject *attr_name)
{
    Py_Shape_file obj = (Py_Shape_file) self;
    Shape_file shape_file = obj->shape_file;

    if (!PyUnicode_Check(attr_name))
    {
        PyErr_Format(PyExc_TypeError, "attribute name must be string, not '%.200s'",
                     Py_TYPE(attr_name)->tp_name);
        return NULL;
    }

    const char *name = PyUnicode_AsUTF8(attr_name);
    if (!name)
        return NULL;

    if (strcmp(name, "ndim") == 0 || strcmp(name, "nshapes") == 0)
        return PyLong_FromLong((long) shape_file->ndim);
    else if (strcmp(name, "ncomponents") == 0)
        return PyLong_FromLong((long) shape_file->ncomponents);
    else if (strcmp(name, "points") == 0)
        return get_python_int_list(shape_file->ndim, shape_file->points);

    /* Fall back to PyObject_GenericGetAttr which checks tp_methods */
    return PyObject_GenericGetAttr(self, attr_name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

static PyObject *shape_file_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = { "ncomponents", "points", NULL };
    int ncomponents;
    PyObject *points_obj;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iO", kwlist,
                                     &ncomponents, &points_obj))
    {
        RETURN_OBJ_ERROR("must have arguments: ncomponents, points");
    }

    return new_py_shape_file(ncomponents, points_obj);
}

static PyTypeObject Shape_file_type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    "ShapeFile",                          /* tp_name */
    sizeof(struct Py_Shape_file),         /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor) delete_py_shape_file,    /* tp_dealloc */
    0,                                    /* tp_vectorcall */
    0,                                    /* tp_getattr */
    0,                                    /* tp_setattr */
    0,                                    /* tp_as_async */
    0,                                    /* tp_repr */
    0,                                    /* tp_as_number */
    0,                                    /* tp_as_sequence */
    0,                                    /* tp_as_mapping */
    0,                                    /* tp_hash */
    0,                                    /* tp_call */
    0,                                    /* tp_str */
    (getattrofunc) getattr_py_shape_file, /* tp_getattro */
    0,                                    /* tp_setattro */
    0,                                    /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                   /* tp_flags */
    "ShapeFile — NMR line shape model",   /* tp_doc */
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                    /* tp_weaklistoffset */
    0,                                    /* tp_iter */
    0,                                    /* tp_iternext */
    py_handler_methods,                   /* tp_methods */
    0,                                    /* tp_members */
    0,                                    /* tp_getset */
    0,                                    /* tp_base */
    0,                                    /* tp_dict */
    0,                                    /* tp_descr_get */
    0,                                    /* tp_descr_set */
    0,                                    /* tp_dictoffset */
    0,                                    /* tp_init */
    0,                                    /* tp_alloc */
    shape_file_new,                       /* tp_new */
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

static struct PyMethodDef Shape_file_type_methods[] =
{
    { "ShapeFile",	(PyCFunction) init_Py_Shape_file,	METH_VARARGS },
    { NULL,		NULL,			0 }
};

/*
 * Python 3 module definition (replaces Py2 Py_InitModule)
 */
static struct PyModuleDef shape_file_module_def =
{
    PyModuleDef_HEAD_INIT,
    "ShapeFile",
    "CCPNMR Line Shape File module (Python 3 compatible)",
    -1,
    Shape_file_type_methods
};

PyMODINIT_FUNC PyInit_ShapeFile(void)
{
    if (PyType_Ready(&Shape_file_type) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&shape_file_module_def);
    if (!m)
        return NULL;

    /* Add the custom type to the module dict */
    PyObject *d = PyModule_GetDict(m);
    if (PyDict_SetItemString(d, "ShapeFile", (PyObject *) &Shape_file_type) < 0)
    {
        Py_DECREF(m);
        return NULL;
    }

    /* Create exception object and add to module */
    ErrorObject = PyErr_NewException("ShapeFile.error", NULL, NULL);
    if (ErrorObject != NULL)
    {
        Py_INCREF(ErrorObject);
        if (PyDict_SetItemString(d, "error", ErrorObject) < 0)
        {
            Py_DECREF(m);
            return NULL;
        }
    }

    return m;
}

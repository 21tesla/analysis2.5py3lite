/*
======================COPYRIGHT/LICENSE START==========================

py_contour_levels.c: Part of the CcpNmr Analysis program

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
#include "py_contour_levels.h"

#include "python_util.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Contour_levels_type;

Bool is_py_contour_levels(PyObject *obj)
{
    /*  below does not work because different *.so files end up
        with different addresses for Contour_levels_type
        return (obj->ob_type == &Contour_levels_type);
    */
    return valid_py_object(obj, &Contour_levels_type);
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

static PyObject *new_py_contour_levels(PyObject *levels_obj)
{
    int n;
    float *levels;
    Contour_levels contour_levels;
    Py_Contour_levels obj;
    Line error_msg;

    if (get_python_float_alloc_array(levels_obj, &n, &levels,
                                     error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    contour_levels = new_contour_levels(n, levels);

    /*  new_contour_levels creates own copy  */
    FREE(levels, float);

    if (!contour_levels)
        RETURN_OBJ_ERROR("allocating Contour_levels object");

    PY_MALLOC(obj, struct Py_Contour_levels, &Contour_levels_type);

    if (!obj)
    {
        delete_contour_levels(contour_levels);

        RETURN_OBJ_ERROR("allocating Py_Contour_levels object");
    }

    obj->contour_levels = contour_levels;

    return (PyObject *) obj;
}

static void delete_py_contour_levels(PyObject *self)
{
    Py_Contour_levels obj = (Py_Contour_levels) self;
    Contour_levels contour_levels = obj->contour_levels;

    /*
        printf("in delete_py_contour_levels\n");
    */

    delete_contour_levels(contour_levels);

    PY_FREE(self);
}

/*
static int print_py_contour_levels(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_contour_levels(PyObject *self, char *name)
{
    /*
        Contour_levels *obj = (Contour_levels *) self;
        Random_access *a = obj->py_handler;

        if (equal_strings(name, "par_file"))
    	return Py_BuildValue("s", a->par_file);
        else if (equal_strings(name, "access_method"))
    	return Py_BuildValue("s", access_method_name(a->access_method));
        else if (equal_strings(name, "data_format"))
    	return get_Contour_levels_format(a);
        else
    */
    return Py_FindMethod(py_handler_methods, self, name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Contour_levels_sequence_methods =
{
    Contour_levels_length,
    Contour_levels_concat,
    Contour_levels_repeat,
    Contour_levels_item,
    Contour_levels_slice,
    Contour_levels_ass_item,
    Contour_levels_ass_slice
};

static PySequenceMethods Contour_levels_sequence_methods =
{
    Contour_levels_length,
    0,
    0,
    Contour_levels_item,
    0,
    Contour_levels_ass_item,
    0
};
*/

static PyTypeObject Contour_levels_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "ContourLevels", /* name */
    sizeof(struct Py_Contour_levels), /* basicsize */
    0, /* itemsize */
    delete_py_contour_levels, /* destructor */
    0, /* printfunc */
    getattr_py_contour_levels, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Contour_levels_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Contour_levels(PyObject *self, PyObject *args)
{
    PyObject *levels_obj, *obj;

    if (!PyArg_ParseTuple(args, "O", &levels_obj))
        RETURN_OBJ_ERROR("must have one argument: levels");

    obj = new_py_contour_levels(levels_obj);

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


static struct PyMethodDef Contour_levels_type_methods[] =
{
    { "ContourLevels",	(PyCFunction) init_Py_Contour_levels,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import ContourLevels" in
* a Python program. The function is usually called "initContour_levels": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initContourLevels(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Contour_levels_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("ContourLevels", Contour_levels_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("ContourLevels.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);

    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module ContourLevels");
}


/*
======================COPYRIGHT/LICENSE START==========================

py_block_file.c: Part of the CcpNmr Analysis program

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
#include "py_block_file.h"

#include "py_mem_cache.h"
#include "python_util.h"
#include "utility.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Block_file_type;
static PyTypeObject Shape_block_file_type;

Bool is_py_block_file(PyObject *obj)
{
/*  below does not work because different *.so files end up
    with different addresses for Block_file_type
    return (obj->ob_type == &Block_file_type);
*/
    return valid_py_object(obj, &Block_file_type);
}

Bool is_py_shape_block_file(PyObject *obj)
{
    return valid_py_object(obj, &Shape_block_file_type);
}

Bool is_py_some_block_file(PyObject *obj)
{
    return is_py_block_file(obj) || is_py_shape_block_file(obj);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *openFile(PyObject *self, PyObject *args)
{
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    Line error_msg;

    if (open_block_file(block_file) == CCPN_ERROR)
    {
	sprintf(error_msg, "opening block file '%s'", block_file->file);
	RETURN_OBJ_ERROR(error_msg);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *closeFile(PyObject *self, PyObject *args)
{
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;

    close_block_file(block_file);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *havePeak(PyObject *self, PyObject *args, PyObject *keywds)
{
    int n, first[MAX_NDIM], last[MAX_NDIM];
    Bool have_peak, have_low = CCPN_FALSE, have_high = CCPN_FALSE;
    float low = 0, high = 0;
    PyObject *first_obj, *last_obj;
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    static char *kwlist[] = { "first", "last", "have_low", "have_high",
					"low", "high", NULL };
    Line error_msg;
 
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO|iiff", kwlist,
		&first_obj, &last_obj,
		&have_low, &have_high, &low, &high))
	RETURN_OBJ_ERROR("need arguments: first, last, [ have_low, have_high, low, high ]");

    if ((get_python_int_array(first_obj, MAX_NDIM, &n, first,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "first must be int list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if ((get_python_int_array(last_obj, MAX_NDIM, &n, last,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "last must be int list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }

    if (have_peak_block_file(block_file, first, last, have_low, have_high,
			low, high, &have_peak, error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    return Py_BuildValue("i", (int) have_peak);
}

static PyObject *maxValue(PyObject *self, PyObject *args)
{
    int n, first[MAX_NDIM], last[MAX_NDIM];
    float max_value;
    PyObject *first_obj, *last_obj;
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "OO", &first_obj, &last_obj))
	RETURN_OBJ_ERROR("need arguments: first, last");

    if ((get_python_int_array(first_obj, MAX_NDIM, &n, first,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "first must be int list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if ((get_python_int_array(last_obj, MAX_NDIM, &n, last,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "last must be int list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }

    if (max_value_block_file(block_file, first, last, &max_value,
						error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    return Py_BuildValue("f", max_value);
}

static PyObject *minValue(PyObject *self, PyObject *args)
{
    int n, first[MAX_NDIM], last[MAX_NDIM];
    float min_value;
    PyObject *first_obj, *last_obj;
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "OO", &first_obj, &last_obj))
	RETURN_OBJ_ERROR("need arguments: first, last");

    if ((get_python_int_array(first_obj, MAX_NDIM, &n, first,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "first must be int list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if ((get_python_int_array(last_obj, MAX_NDIM, &n, last,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "last must be int list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }

    if (min_value_block_file(block_file, first, last, &min_value,
						error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    return Py_BuildValue("f", min_value);
}

static PyObject *getValue(PyObject *self, PyObject *args)
{
    int n;
    float value, position[MAX_NDIM];
    PyObject *posn_obj;
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &posn_obj))
	RETURN_OBJ_ERROR("need argument: position");

    if ((get_python_float_array(posn_obj, MAX_NDIM, &n, position,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "position must be float list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if (get_position_block_file(block_file, &value, position,
						error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    return Py_BuildValue("f", value);
}

static PyObject *getPointValue(PyObject *self, PyObject *args)
{
    int n, position[MAX_NDIM];
    float value;
    PyObject *posn_obj;
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &posn_obj))
	RETURN_OBJ_ERROR("need argument: position");

    if ((get_python_int_array(posn_obj, MAX_NDIM, &n, position,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "position must be int list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if (get_point_block_file(block_file, &value, position,
						error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    return Py_BuildValue("f", value);
}

static PyObject *getValues(PyObject *self, PyObject *args)
{
    int n, npoints, box_min[MAX_NDIM], box_max[MAX_NDIM];
    float *values;
    PyObject *box_min_obj, *box_max_obj, *values_obj;
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "OO", &box_min_obj, &box_max_obj))
	RETURN_OBJ_ERROR("need two arguments: box_min, box_max");

    if ((get_python_int_array(box_min_obj, MAX_NDIM, &n, box_min,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "box_min must be int list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if ((get_python_int_array(box_max_obj, MAX_NDIM, &n, box_max,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "box_max must be int list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if (get_box_block_file(block_file, &npoints, &values,
                                box_min, box_max, error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    values_obj = get_python_float_list(npoints, values);
    FREE(values, float);

    if (values_obj == NULL)
        RETURN_OBJ_ERROR("allocating values memory");

    return values_obj;
}

static PyObject *setValue(PyObject *self, PyObject *args)
{
    int n, position[MAX_NDIM];
    float value;
    PyObject *posn_obj;
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "Of", &posn_obj, &value))
	RETURN_OBJ_ERROR("need arguments: position value");

    if ((get_python_int_array(posn_obj, MAX_NDIM, &n, position,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "position must be int list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if (set_point_block_file(block_file, value, position, error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setValues(PyObject *self, PyObject *args)
{
    int n, npoints, box_min[MAX_NDIM], box_max[MAX_NDIM];
    float *values;
    PyObject *box_min_obj, *box_max_obj, *values_obj;
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "OOO", &box_min_obj, &box_max_obj, &values_obj))
	RETURN_OBJ_ERROR("need arguments: box_min, box_max, values");

    if ((get_python_int_array(box_min_obj, MAX_NDIM, &n, box_min,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "box_min must be int list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if ((get_python_int_array(box_max_obj, MAX_NDIM, &n, box_max,
						error_msg) == CCPN_ERROR)
	|| (n != block_file->ndim))
    {
	sprintf(error_msg, "box_max must be int list or tuple of size %d",
							block_file->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if (get_python_float_alloc_array(values_obj, &npoints, &values,
						error_msg) == CCPN_ERROR)
    {
	sprintf(error_msg, "values must be float list or tuple");
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if (set_box_block_file(block_file, npoints, values,
                                box_min, box_max, error_msg) == CCPN_ERROR)
    {
        FREE(values, float);
	RETURN_OBJ_ERROR(error_msg);
    }

    FREE(values, float);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *save(PyObject *self, PyObject *args)
{
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    Line error_msg;

    if (!PyArg_ParseTuple(args, ""))
	RETURN_OBJ_ERROR("should have no arguments");

    if (save_block_file(block_file, error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *check(PyObject *self, PyObject *args)
{
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    Line error_msg;

    if (check_block_file(block_file, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *planeProduct(PyObject *self, PyObject *args)
{
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    int i, n, nplanes;
    float scale = 1;
    Block_plane *block_planes, *bp;
    PyObject *planesObj, *tobj, *blockObj, *planeObj;
    CcpnStatus status;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "Of", &planesObj, &scale))
	RETURN_OBJ_ERROR("need two arguments: list/tuple of (block_file, xdim, xmin, xmax, ydim, ymin, ymax, plane) and scale");

    nplanes = get_python_list_size(planesObj);
    if (nplanes < 0)
	RETURN_OBJ_ERROR("first argument needs to be list/tuple of (block_file, xdim, xmin, xmax, ydim, ymin, ymax, plane)");

    block_planes = alloc_block_planes(nplanes);
    if (!block_planes)
	RETURN_OBJ_ERROR("allocating block_planes");

    for (i = 0; i < nplanes; i++)
    {
	bp = block_planes + i;

        tobj = get_python_object_by_index(planesObj, i);
        if (!PyArg_ParseTuple(tobj, "OiffiffO", &blockObj, &bp->xdim, &bp->xmin, &bp->xmax, &bp->ydim, &bp->ymin, &bp->ymax, &planeObj))
        {
	    free_block_planes(block_planes);
            sprintf(error_msg, "Element %d: must be tuple (blockFile, xdim, xmin, xmax, ydim, ymin, ymax, planes)", i);
            RETURN_OBJ_ERROR(error_msg);
        }

	if (!is_py_block_file(blockObj))
        {
	    free_block_planes(block_planes);
            sprintf(error_msg, "Element %d: first element of this must be blockFile", i);
            RETURN_OBJ_ERROR(error_msg);
        }

	bp->block_file = ((Py_Block_file) blockObj)->block_file;

	if (get_python_float_array(planeObj, MAX_NDIM, &n, bp->plane, error_msg) == CCPN_ERROR)
	{
	    free_block_planes(block_planes);
            RETURN_OBJ_ERROR(error_msg);
	}

	if (n != bp->block_file->ndim-2)
	{
	    free_block_planes(block_planes);
	    sprintf(error_msg, "Element %d: plane has size %d, should be %d (ndim-2)", i, n, bp->block_file->ndim-2);
            RETURN_OBJ_ERROR(error_msg);
	}
    }

    status = plane_product_block_files(block_file, nplanes, block_planes, scale, error_msg);

    free_block_planes(block_planes);

    if (status == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef py_handler_methods[] =
{
    { "open",		openFile,		METH_VARARGS },
    { "close",		closeFile,		METH_VARARGS },
    { "havePeak",    (PyCFunction) havePeak,	METH_VARARGS|METH_KEYWORDS },
    { "maxValue",	maxValue,		METH_VARARGS },
    { "minValue",	minValue,		METH_VARARGS },
    { "getValue",	getValue,		METH_VARARGS },
    { "getPointValue",	getPointValue,		METH_VARARGS },
    { "getValues",	getValues,		METH_VARARGS },
    { "setValue",	setValue,		METH_VARARGS },
    { "setValues",	setValues,		METH_VARARGS },
    { "save",		save,			METH_VARARGS },
    { "check",		check,			METH_VARARGS },
    { "planeProduct",	planeProduct,		METH_VARARGS },
    { NULL,		NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_block_file(char *file, int ndim,
		PyObject *points_obj, PyObject *block_size_obj,
		PyObject *dim_wrapped_obj,
		PyObject *mem_cache_obj, int bytes_per_point,
		Bool big_endian, Bool padded, int header,
		Bool integer, Bool writeable, int blockHeader)
{
    int i, n, points[MAX_NDIM], block_size[MAX_NDIM];
    int dim_wrapped_int[MAX_NDIM];
    Bool dim_wrapped[MAX_NDIM];
    Block_file block_file;
    Py_Block_file obj;
    Line error_msg;
 
    if (!is_py_mem_cache(mem_cache_obj))
	RETURN_OBJ_ERROR("must pass mem_cache object");

    if ((get_python_int_array(points_obj, MAX_NDIM, &n, points,
						error_msg) == CCPN_ERROR)
	|| (n != ndim))
    {
	sprintf(error_msg, "points must be list or tuple of size %d", ndim);
	RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_int_array(block_size_obj, MAX_NDIM, &n, block_size,
						error_msg) == CCPN_ERROR)
	|| (n != ndim))
    {
	sprintf(error_msg, "block_size must be list or tuple of size %d", ndim);
	RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_int_array(dim_wrapped_obj, MAX_NDIM, &n, dim_wrapped_int,
						error_msg) == CCPN_ERROR)
	|| (n != ndim))
    {
	sprintf(error_msg, "dim_wrapped must be list or tuple of size %d", ndim);
	RETURN_OBJ_ERROR(error_msg);
    }

    for (i = 0; i < ndim; i++)
	dim_wrapped[i] = dim_wrapped_int[i] ? CCPN_TRUE : CCPN_FALSE;

    block_file = new_block_file(file, ndim, points, block_size, dim_wrapped,
			((Py_Mem_cache) mem_cache_obj)->mem_cache,
			bytes_per_point, big_endian, padded,
			header, integer, writeable, blockHeader);

    if (!block_file)
	RETURN_OBJ_ERROR("allocating Block_file object");

    PY_MALLOC(obj, struct Py_Block_file, &Block_file_type);

    if (!obj)
    {
	delete_block_file(block_file);

	RETURN_OBJ_ERROR("allocating Py_Block_file object");
    }

    obj->block_file = block_file;

    return (PyObject *) obj;
}

static PyObject *new_py_shape_block_file(char *file, int ndim,
		PyObject *points_obj, PyObject *block_size_obj,
		PyObject *dim_wrapped_obj, PyObject *shape_file_obj)
{
    int i, n, points[MAX_NDIM], block_size[MAX_NDIM];
    int dim_wrapped_int[MAX_NDIM];
    Bool dim_wrapped[MAX_NDIM];
    Block_file block_file;
    Shape_file shape_file;
    Py_Shape_block_file obj;
    Line error_msg;
 
    if (!is_py_shape_file(shape_file_obj))
	RETURN_OBJ_ERROR("must pass shape_file object");

    if ((get_python_int_array(points_obj, MAX_NDIM, &n, points,
						error_msg) == CCPN_ERROR)
	|| (n != ndim))
    {
	sprintf(error_msg, "points must be list or tuple of size %d", ndim);
	RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_int_array(block_size_obj, MAX_NDIM, &n, block_size,
						error_msg) == CCPN_ERROR)
	|| (n != ndim))
    {
	sprintf(error_msg, "block_size must be list or tuple of size %d", ndim);
	RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_int_array(dim_wrapped_obj, MAX_NDIM, &n, dim_wrapped_int,
						error_msg) == CCPN_ERROR)
	|| (n != ndim))
    {
	sprintf(error_msg, "dim_wrapped must be list or tuple of size %d", ndim);
	RETURN_OBJ_ERROR(error_msg);
    }

    for (i = 0; i < ndim; i++)
	dim_wrapped[i] = dim_wrapped_int[i] ? CCPN_TRUE : CCPN_FALSE;

    block_file = new_shape_block_file(file, ndim, points, block_size, dim_wrapped,
			((Py_Shape_file) shape_file_obj)->shape_file);

    if (!block_file)
	RETURN_OBJ_ERROR("allocating Shape_block_file object");

    PY_MALLOC(obj, struct Py_Shape_block_file, &Shape_block_file_type);

    if (!obj)
    {
	delete_block_file(block_file);

	RETURN_OBJ_ERROR("allocating Py_Block_file object");
    }

    obj->shape_block_file = (Shape_block_file) block_file;
    obj->shape_file_obj = (Py_Shape_file) shape_file_obj;
    Py_INCREF(shape_file_obj);

    return (PyObject *) obj;
}

static void delete_py_block_file(PyObject *self)
{
    Py_Block_file obj = (Py_Block_file) self;
    Block_file block_file = obj->block_file;
    Py_Shape_block_file shape_block_file_obj;
    Shape_block_file shape_block_file;

/*
    printf("in delete_py_block_file\n");
*/

    if (block_file->block_file_kind == SHAPE_FILE)
    {
	/* just decref Python object, do not clean up shape_file memory */
	/* so make shape_file NULL so doesn't get cleaned up in delete here */
	shape_block_file_obj = (Py_Shape_block_file) self;
	shape_block_file = (Shape_block_file) block_file;
        Py_DECREF(shape_block_file_obj->shape_file_obj);
	shape_block_file->shape_file = NULL;
    }

    delete_block_file(block_file);

    PY_FREE(self);
}

/*
static int print_py_block_file(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_block_file(PyObject *self, char *name)
{
    Py_Block_file obj = (Py_Block_file) self;
    Py_Shape_block_file shape_block_file_obj;
    Py_Shape_file shape_file_obj;
    Block_file block_file = obj->block_file;

    if (equal_strings(name, "ndim"))
	return Py_BuildValue("i", block_file->ndim);

    else if ((block_file->block_file_kind == SHAPE_FILE) && equal_strings(name, "shapeFile"))
    {
	shape_block_file_obj = (Py_Shape_block_file) obj;
	shape_file_obj = shape_block_file_obj->shape_file_obj;
	Py_INCREF(shape_file_obj);

	return (PyObject *) shape_file_obj;
    }

    else
	return Py_FindMethod(py_handler_methods, self, name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Block_file_sequence_methods =
{
    Block_file_length,
    Block_file_concat,
    Block_file_repeat,
    Block_file_item,
    Block_file_slice,
    Block_file_ass_item,
    Block_file_ass_slice
};

static PySequenceMethods Block_file_sequence_methods =
{
    Block_file_length,
    0,
    0,
    Block_file_item,
    0,
    Block_file_ass_item,
    0
};
*/

static PyTypeObject Block_file_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "BlockFile", /* name */
    sizeof(struct Py_Block_file), /* basicsize */
    0, /* itemsize */
    delete_py_block_file, /* destructor */
    0, /* printfunc */
    getattr_py_block_file, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Block_file_sequence_methods*/ /* PySequenceMethods */
};

static PyTypeObject Shape_block_file_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "ShapeBlockFile", /* name */
    sizeof(struct Py_Shape_block_file), /* basicsize */
    0, /* itemsize */
    delete_py_block_file, /* destructor */
    0, /* printfunc */
    getattr_py_block_file, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Block_file_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Block_file(PyObject *self, PyObject *args)
{
    int ndim, bytes_per_point, big_endian, padded, header, integer, writeable = 0, block_header = 0;
    PyObject *points_obj, *block_size_obj, *dim_wrapped_obj, *mem_cache_obj, *obj;
    char *file;

    if (!PyArg_ParseTuple(args, "siOOOOiiiii|ii", &file, &ndim,
		&points_obj, &block_size_obj, &dim_wrapped_obj,
		&mem_cache_obj, &bytes_per_point, &big_endian,
		&padded, &header, &integer, &writeable, &block_header))
        RETURN_OBJ_ERROR("must have arguments: file, ndim, points, block_size, dim_wrapped, mem_cache, bytes_per_point, big_endian, padded, header, integer [, writeable]");

    obj = new_py_block_file(file, ndim, points_obj, block_size_obj, dim_wrapped_obj,
			mem_cache_obj, bytes_per_point, big_endian?CCPN_TRUE:CCPN_FALSE,
			padded?CCPN_TRUE:CCPN_FALSE, header, integer?CCPN_TRUE:CCPN_FALSE,
                        writeable?CCPN_TRUE:CCPN_FALSE, block_header);

/*  this happens automatically??
    FREE(file, char);
*/

    return obj;
}

static PyObject *init_Py_Shape_block_file(PyObject *self, PyObject *args)
{
    int ndim;
    PyObject *points_obj, *block_size_obj, *dim_wrapped_obj, *shape_file_obj, *obj;
    char *file;

    if (!PyArg_ParseTuple(args, "siOOOO", &file, &ndim,
		&points_obj, &block_size_obj, &dim_wrapped_obj, &shape_file_obj))
        RETURN_OBJ_ERROR("must have arguments: file, ndim, points, block_size, dim_wrapped, shape_file");

    obj = new_py_shape_block_file(file, ndim, points_obj, block_size_obj, dim_wrapped_obj, shape_file_obj);

/*  this happens automatically??
    FREE(file, char);
*/

    return obj;
}

static PyObject *compareSlices(PyObject *self, PyObject *args)
{
    PyObject *block_file1_obj, *block_file2_obj, *slice1_obj, *slice2_obj;
    Block_file block_file1, block_file2;
    int n, dim1, dim2, slice1[MAX_NDIM], slice2[MAX_NDIM];
    float v, scale, offset;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "OiOOiOff", &block_file1_obj, &dim1, &slice1_obj,
			&block_file2_obj, &dim2, &slice2_obj, &scale, &offset))
	RETURN_OBJ_ERROR("need arguments: block_file1, dim1, slice1, block_file2,dim2, slice2, scale, offset");

    if (!is_py_block_file(block_file1_obj)) 
        RETURN_OBJ_ERROR("block_file1 (first argument) must be block_file object");

    block_file1 = ((Py_Block_file) block_file1_obj)->block_file;

    if (!is_py_block_file(block_file2_obj)) 
        RETURN_OBJ_ERROR("block_file2 (fourth argument) must be block_file object");

    block_file2 = ((Py_Block_file) block_file2_obj)->block_file;

    if ((get_python_int_array(slice1_obj, MAX_NDIM, &n, slice1,
						error_msg) == CCPN_ERROR)
	|| (n != block_file1->ndim))
    {
	sprintf(error_msg, "slice1 (third argument) must be int list or tuple of size %d",
							block_file1->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if ((get_python_int_array(slice2_obj, MAX_NDIM, &n, slice2,
						error_msg) == CCPN_ERROR)
	|| (n != block_file2->ndim))
    {
	sprintf(error_msg, "slice2 (sixth argument) must be int list or tuple of size %d",
							block_file2->ndim);
	RETURN_OBJ_ERROR(error_msg);
    }
 
    if (compare_slices_block_files(block_file1, dim1, slice1,
		block_file2, dim2, slice2, &v, scale, offset, error_msg) == CCPN_ERROR)
	RETURN_OBJ_ERROR(error_msg);

    return Py_BuildValue("f", v);
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Block_file_type_methods[] =
{
    { "BlockFile",	(PyCFunction) init_Py_Block_file,	METH_VARARGS },
    { "ShapeBlockFile",	(PyCFunction) init_Py_Shape_block_file,	METH_VARARGS },
    { "compareSlices",	(PyCFunction) compareSlices,		METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import BlockFile" in 
* a Python program. The function is usually called "initBlock_file": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initBlockFile(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Block_file_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("BlockFile", Block_file_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("BlockFile.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module BlockFile");
}

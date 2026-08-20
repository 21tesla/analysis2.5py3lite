/*
======================COPYRIGHT/LICENSE START==========================

py_peak_list.c: Part of the CcpNmr Analysis program

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
#include "py_peak_list.h"

#include "py_block_file.h"
#include "py_peak.h"
#include "py_pdf_handler.h"
#include "py_ps_handler.h"
#include "python_util.h"

#define  NCOLORS  3

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

static PyTypeObject Peak_list_type;

Bool is_py_peak_list(PyObject *obj)
{
    /*  below does not work because different *.so files end up
        with different addresses for Peak_list_type
        return (obj->ob_type == &Peak_list_type);
    */
    return valid_py_object(obj, &Peak_list_type);
}

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

static CcpnStatus alloc_diagonal_memory(Diagonal_exclusion *diagonal_exclusion, int nexclusions)
{
    Diagonal_exclusion de;

    if (nexclusions == 0)
        return CCPN_OK;

    MALLOC(de, struct Diagonal_exclusion, nexclusions);

    *diagonal_exclusion = de;

    return CCPN_OK;
}

static void free_diagonal_memory(Diagonal_exclusion diagonal_exclusion, int nexclusions)
{
    if (nexclusions == 0)
        return;

    FREE(diagonal_exclusion, struct Diagonal_exclusion);
}

static CcpnStatus alloc_exclude_memory(float ****regions, int nregions, int ndim)
{
    int i, j;
    float ***r;

    if (nregions == 0)
        return CCPN_OK;

    MALLOC(r, float **, nregions);
    for (i = 0; i < nregions; i++)
    {
        MALLOC(r[i], float *, ndim);
        for (j = 0; j < ndim; j++)
            MALLOC(r[i][j], float, 2);
    }

    *regions = r;

    return CCPN_OK;
}

static void free_exclude_memory(float ***regions, int nregions, int ndim)
{
    int i, j;

    if (nregions == 0)
        return;

    for (i = 0; i < nregions; i++)
    {
        for (j = 0; j < ndim; j++)
            FREE(regions[i][j], float);

        FREE(regions[i], float *);
    }

    FREE(regions, float **);
}

static CcpnStatus get_int(PyObject *w, Bool have_list, int n, int *value, CcpnString error_msg)
{
    PyObject *z;

    if (have_list)
        z = PyList_GetItem(w, n);
    else
        z = PyTuple_GetItem(w, n);

    if (!PyInt_Check(z))
    {
        sprintf(error_msg, "item %d must be an int", n);
        return CCPN_ERROR;
    }

    *value = (int) PyInt_AsLong(z);

    return CCPN_OK;
}

static CcpnStatus get_float(PyObject *w, Bool have_list, int n, float *value, CcpnString error_msg)
{
    PyObject *z;

    if (have_list)
        z = PyList_GetItem(w, n);
    else
        z = PyTuple_GetItem(w, n);

    if (!PyFloat_Check(z) && !PyInt_Check(z))
    {
        sprintf(error_msg, "item %d must be an int or float", n);
        return CCPN_ERROR;
    }

    if (PyFloat_Check(z))
        *value = (float) PyFloat_AsDouble(z);
    else
        *value = (float) PyInt_AsLong(z);

    return CCPN_OK;
}

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

static PyObject *setIsSelected(PyObject *self, PyObject *args)
{
    Bool isSelected;
    Py_Peak_list obj = (Py_Peak_list) self;
    Peak_list peak_list = obj->peak_list;

    if (!PyArg_ParseTuple(args, "i", &isSelected))
        RETURN_OBJ_ERROR("need one argument: isSelected");

    is_selected_peak_list(peak_list, isSelected);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setColor(PyObject *self, PyObject *args)
{
    int n;
    float color[NCOLORS];
    PyObject *color_obj;
    Py_Peak_list obj = (Py_Peak_list) self;
    Peak_list peak_list = obj->peak_list;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &color_obj))
        RETURN_OBJ_ERROR("need one argument: color");

    if ((get_python_float_array(color_obj, NCOLORS, &n, color,
                                error_msg) == CCPN_ERROR)
            || (n != NCOLORS))
    {
        sprintf(error_msg, "color must be list or tuple of size %d", NCOLORS);
        RETURN_OBJ_ERROR(error_msg);
    }

    set_color_peak_list(peak_list, color);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *setSymbol(PyObject *self, PyObject *args)
{
    int symbol;
    Py_Peak_list obj = (Py_Peak_list) self;
    Peak_list peak_list = obj->peak_list;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "i", &symbol))
        RETURN_OBJ_ERROR("need one argument: symbol");

    CHECK_OBJ_STATUS(set_symbol_peak_list(peak_list, symbol, error_msg));

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *addPeak(PyObject *self, PyObject *args)
{
    Py_Peak_list obj = (Py_Peak_list) self;
    Peak_list peak_list = obj->peak_list;
    Peak peak;
    Line error_msg;

    CHECK_OBJ_STATUS(add_peak_peak_list(peak_list, &peak, error_msg));

    return new_py_peak(peak);
}

static PyObject *removePeak(PyObject *self, PyObject *args)
{
    PyObject *peak_obj;
    Py_Peak_list obj = (Py_Peak_list) self;
    Peak_list peak_list = obj->peak_list;
    Peak peak;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &peak_obj))
        RETURN_OBJ_ERROR("need one argument: peak");

    if (!is_py_peak(peak_obj))
        RETURN_OBJ_ERROR("need one argument: peak");

    peak = ((Py_Peak) peak_obj)->peak;

    CHECK_OBJ_STATUS(remove_peak_peak_list(peak_list, peak, error_msg));

    /*  this done automatically by garbage collection
        delete_py_peak(peak_obj);
    */

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *removeSelectedPeaks(PyObject *self, PyObject *args)
{
    Py_Peak_list obj = (Py_Peak_list) self;
    Peak_list peak_list = obj->peak_list;

    remove_selected_peak_list(peak_list);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *unselectSelectedPeaks(PyObject *self, PyObject *args)
{
    Py_Peak_list obj = (Py_Peak_list) self;
    Peak_list peak_list = obj->peak_list;

    unselect_selected_peak_list(peak_list);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *findPeaks(PyObject *self, PyObject *args, PyObject *keywds)
{
    int i, j, n, first[MAX_NDIM], last[MAX_NDIM], buffer[MAX_NDIM], npeaks;
    int ndim, nexcluded_regions = 0, ndiagonal_exclusions = 0;
    int ignore_peak = -1;
    int dc[MAX_NDIM];
    Bool have_list, have_list2, dim_checked[MAX_NDIM];
    float min_linewidth[MAX_NDIM], ***excluded_regions;
    Bool nonadjacent = CCPN_FALSE, have_low = CCPN_FALSE, have_high = CCPN_FALSE;
    float low = 0, high = 0, drop_factor = 0;
    PyObject *first_obj, *last_obj, *block_file_obj;
    PyObject *min_linewidth_obj = NULL;
    PyObject *diagonal_exclusion_obj = NULL, *excluded_region_obj = NULL;
    PyObject *dim_checked_obj = NULL;
    Diagonal_exclusion de, diagonal_exclusions;
    PyObject *buffer_obj = NULL, *list, *p, *z, *z2;
    Py_Peak_list obj = (Py_Peak_list) self;
    Peak_list peak_list = obj->peak_list;
    static char *kwlist[] = { "first", "last", "block_file",
                              "have_low", "have_high", "low", "high", "buffer", "nonadjacent",
                              "drop_factor", "min_linewidth", "diagonal_exclusions",
                              "excluded_regions", "dim_checked", "ignore_peak", NULL
                            };
    Block_file block_file;
    Line error_msg;

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OOO|iiffOifOOOOi", kwlist,
                                     &first_obj, &last_obj, &block_file_obj,
                                     &have_low, &have_high, &low, &high, &buffer_obj, &nonadjacent,
                                     &drop_factor, &min_linewidth_obj, &diagonal_exclusion_obj,
                                     &excluded_region_obj, &dim_checked_obj, &ignore_peak))
        RETURN_OBJ_ERROR("need arguments: first, last, block_file, [ have_low, have_high, low, high, buffer, nonadjacent, drop_factor, min_linewidth, diagonal_exclusions, excluded_regions, dim_checked, ignore_peak ]");

    ndim = peak_list->ndim;

    if ((get_python_int_array(first_obj, MAX_NDIM, &n, first,
                              error_msg) == CCPN_ERROR)
            || (n != ndim))
    {
        sprintf(error_msg, "first must be int list or tuple of size %d",
                ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_int_array(last_obj, MAX_NDIM, &n, last,
                              error_msg) == CCPN_ERROR)
            || (n != ndim))
    {
        sprintf(error_msg, "last must be int list or tuple of size %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (!is_py_some_block_file(block_file_obj))
        RETURN_OBJ_ERROR("must have BlockFile object");

    if (is_py_block_file(block_file_obj))
        block_file = ((Py_Block_file) block_file_obj)->block_file;
    else
        block_file = (Block_file) (((Py_Shape_block_file) block_file_obj)->shape_block_file);

    if (buffer_obj)
    {
        if ((get_python_int_array(buffer_obj, MAX_NDIM, &n, buffer,
                                  error_msg) == CCPN_ERROR)
                || (n != ndim))
        {
            sprintf(error_msg, "buffer must be int list or tuple of size %d", ndim);
            RETURN_OBJ_ERROR(error_msg);
        }
    }
    else
    {
        ZERO_VECTOR(buffer, ndim);
    }

    if (min_linewidth_obj)
    {
        if ((get_python_float_array(min_linewidth_obj, MAX_NDIM, &n,
                                    min_linewidth, error_msg) == CCPN_ERROR)
                || (n != ndim))
        {
            sprintf(error_msg, "min_linewidth must be float list or tuple of size %d", ndim);
            RETURN_OBJ_ERROR(error_msg);
        }
    }
    else
    {
        ZERO_VECTOR(min_linewidth, ndim);
    }

    if (diagonal_exclusion_obj)
    {
        if (PyList_Check(diagonal_exclusion_obj))
            have_list = CCPN_TRUE;
        else if (PyTuple_Check(diagonal_exclusion_obj))
            have_list = CCPN_FALSE;
        else
            RETURN_OBJ_ERROR("diagonal_exclusions must be list or tuple");

        if (have_list)
            ndiagonal_exclusions = PyList_Size(diagonal_exclusion_obj);
        else
            ndiagonal_exclusions = PyTuple_Size(diagonal_exclusion_obj);

        if (alloc_diagonal_memory(&diagonal_exclusions, ndiagonal_exclusions) == CCPN_ERROR)
            RETURN_OBJ_ERROR("allocating exclude memory");

        for (i = 0; i < ndiagonal_exclusions; i++)
        {
            if (have_list)
                z = PyList_GetItem(diagonal_exclusion_obj, i);
            else
                z = PyTuple_GetItem(diagonal_exclusion_obj, i);

            if (PyList_Check(z))
                have_list2 = CCPN_TRUE;
            else if (PyTuple_Check(z))
                have_list2 = CCPN_FALSE;
            else
            {
                free_diagonal_memory(diagonal_exclusions, ndiagonal_exclusions);
                RETURN_OBJ_ERROR("diagonal_exclusions must be list or tuple");
            }

            if (have_list2)
                n = PyList_Size(z);
            else
                n = PyTuple_Size(z);

            if (n != 6)
            {
                free_diagonal_memory(diagonal_exclusions, ndiagonal_exclusions);
                sprintf(error_msg, "diagonal_exclusions element %d must be list or tuple of size 6, found %d",
                        i, n);
                RETURN_OBJ_ERROR(error_msg);
            }

            de = diagonal_exclusions + i;
            CHECK_OBJ_STATUS(get_int(z, have_list2, 0, &de->dim1, error_msg));
            CHECK_OBJ_STATUS(get_int(z, have_list2, 1, &de->dim2, error_msg));
            CHECK_OBJ_STATUS(get_float(z, have_list2, 2, &de->a1, error_msg));
            CHECK_OBJ_STATUS(get_float(z, have_list2, 3, &de->a2, error_msg));
            CHECK_OBJ_STATUS(get_float(z, have_list2, 4, &de->b12, error_msg));
            CHECK_OBJ_STATUS(get_float(z, have_list2, 5, &de->d, error_msg));
        }
    }

    if (excluded_region_obj)
    {
        if (PyList_Check(excluded_region_obj))
            have_list = CCPN_TRUE;
        else if (PyTuple_Check(excluded_region_obj))
            have_list = CCPN_FALSE;
        else
            RETURN_OBJ_ERROR("excluded_regions must be list or tuple");

        if (have_list)
            nexcluded_regions = PyList_Size(excluded_region_obj);
        else
            nexcluded_regions = PyTuple_Size(excluded_region_obj);

        if (alloc_exclude_memory(&excluded_regions, nexcluded_regions, ndim) == CCPN_ERROR)
            RETURN_OBJ_ERROR("allocating exclude memory");

        for (i = 0; i < nexcluded_regions; i++)
        {
            if (have_list)
                z = PyList_GetItem(excluded_region_obj, i);
            else
                z = PyTuple_GetItem(excluded_region_obj, i);

            if (PyList_Check(z))
                have_list2 = CCPN_TRUE;
            else if (PyTuple_Check(z))
                have_list2 = CCPN_FALSE;
            else
            {
                free_exclude_memory(excluded_regions, nexcluded_regions, ndim);
                RETURN_OBJ_ERROR("excluded_regions must be list or tuple");
            }

            if (have_list2)
                n = PyList_Size(z);
            else
                n = PyTuple_Size(z);

            if (n != ndim)
            {
                free_exclude_memory(excluded_regions, nexcluded_regions, ndim);
                sprintf(error_msg, "excluded_regions element %d must be list or tuple of size %d, found %d",
                        i, ndim, n);
                RETURN_OBJ_ERROR(error_msg);
            }

            for (j = 0; j < ndim; j++)
            {
                if (have_list2)
                    z2 = PyList_GetItem(z, j);
                else
                    z2 = PyTuple_GetItem(z, j);

                if ((get_python_float_array(z2, 2, &n,
                                            excluded_regions[i][j], error_msg) == CCPN_ERROR)
                        || (n != 2))
                {
                    free_exclude_memory(excluded_regions, nexcluded_regions, ndim);
                    sprintf(error_msg, "excluded_regions element %d %d must be list or tuple of size 2",
                            i, j);
                    RETURN_OBJ_ERROR(error_msg);
                }
            }
        }
    }

    if (dim_checked_obj)
    {
        if ((get_python_int_array(dim_checked_obj, MAX_NDIM, &n,
                                  dc, error_msg) == CCPN_ERROR) || (n != ndim))
        {
            sprintf(error_msg, "dim_checked must be int list or tuple of size %d", ndim);
            RETURN_OBJ_ERROR(error_msg);
        }

        for (i = 0; i < ndim; i++)
            dim_checked[i] = dc[i];
    }
    else
    {
        for (i = 0; i < ndim; i++)
            dim_checked[i] = CCPN_TRUE;
    }

    /* can only pass back positions, peaks must be created in Python world */
    /* so add peaks to peak_list and then delete afterwards */
    npeaks = peak_list->npeaks;

    if (find_peak_list(peak_list, first, last, block_file,
                       have_low, have_high, low, high, buffer, nonadjacent,
                       drop_factor, min_linewidth, ndiagonal_exclusions,
                       diagonal_exclusions, nexcluded_regions,
                       excluded_regions, dim_checked, ignore_peak,
                       error_msg) == CCPN_ERROR)
    {
        free_diagonal_memory(diagonal_exclusions, ndiagonal_exclusions);
        free_exclude_memory(excluded_regions, nexcluded_regions, ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    free_diagonal_memory(diagonal_exclusions, ndiagonal_exclusions);
    free_exclude_memory(excluded_regions, nexcluded_regions, ndim);

    n = peak_list->npeaks - npeaks;

    list = PyList_New(n);
    if (!list)
        RETURN_OBJ_ERROR("allocating peak list memory");

    for (i = 0; i < n; i++)
    {
        p = get_python_float_list(peak_list->ndim,
                                  peak_list->peaks[npeaks+i]->position);
        if (!p)
            RETURN_OBJ_ERROR("allocating peak position memory");
        PyList_SetItem(list, i, p);
        delete_peak(peak_list->peaks[npeaks+i]);
    }

    peak_list->npeaks = npeaks;

    return list;
}

static PyObject *searchPeaks(PyObject *self, PyObject *args, PyObject *keywds)
{
    int i, n, npeaks, *index_list;
    float first[MAX_NDIM], last[MAX_NDIM];
    int aa[MAX_NDIM];
    Bool allow_aliasing[MAX_NDIM];
    PyObject *first_obj, *last_obj, *allow_aliasing_obj = NULL, *list;
    Py_Peak_list obj = (Py_Peak_list) self;
    Peak_list peak_list = obj->peak_list;
    static char *kwlist[] = { "first", "last", "allow_aliasing", NULL };
    Line error_msg;

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO|O", kwlist,
                                     &first_obj, &last_obj, &allow_aliasing_obj))
        RETURN_OBJ_ERROR("need arguments: first, last, [ allow_aliasing ]");

    if ((get_python_float_array(first_obj, MAX_NDIM, &n, first,
                                error_msg) == CCPN_ERROR)
            || (n != peak_list->ndim))
    {
        sprintf(error_msg, "first must be float list or tuple of size %d",
                peak_list->ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_float_array(last_obj, MAX_NDIM, &n, last,
                                error_msg) == CCPN_ERROR)
            || (n != peak_list->ndim))
    {
        sprintf(error_msg, "last must be float list or tuple of size %d",
                peak_list->ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (allow_aliasing_obj)
    {
        if ((get_python_int_array(allow_aliasing_obj, MAX_NDIM, &n,
                                  aa, error_msg) == CCPN_ERROR) || (n != peak_list->ndim))
        {
            sprintf(error_msg, "allow_aliasing must be int list or tuple of size %d", peak_list->ndim);
            RETURN_OBJ_ERROR(error_msg);
        }

        for (i = 0; i < peak_list->ndim; i++)
            allow_aliasing[i] = aa[i];
    }
    else
    {
        for (i = 0; i < peak_list->ndim; i++)
            allow_aliasing[i] = CCPN_TRUE;
    }

    if (search_region_peak_list(peak_list, first, last, allow_aliasing,
                                &npeaks, &index_list, error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    list = get_python_int_list(npeaks, index_list);
    FREE(index_list, int);

    if (!list)
        RETURN_OBJ_ERROR("allocating search list memory");

    return list;
}

static PyObject *nearestPeak(PyObject *self, PyObject *args, PyObject *keywds)
{
    int i, n, index_peak, xdim, ydim;
    float xscale, yscale, first[MAX_NDIM], last[MAX_NDIM], d2_min;
    double d2_mind;
    int aa[MAX_NDIM];
    Bool allow_aliasing[MAX_NDIM];
    PyObject *first_obj, *last_obj, *allow_aliasing_obj = NULL;
    Py_Peak_list obj = (Py_Peak_list) self;
    Peak_list peak_list = obj->peak_list;
    static char *kwlist[] = { "xdim", "ydim", "xscale", "yscale",
                              "first", "last", "allow_aliasing", NULL
                            };
    Line error_msg;

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "iiffOO|O", kwlist,
                                     &xdim, &ydim, &xscale, &yscale, &first_obj, &last_obj, &allow_aliasing_obj))
        RETURN_OBJ_ERROR("need arguments: xdim, ydim, xscale, yscale, first, last, [ allow_aliasing ]");

    if ((get_python_float_array(first_obj, MAX_NDIM, &n, first,
                                error_msg) == CCPN_ERROR)
            || (n != peak_list->ndim))
    {
        sprintf(error_msg, "first must be float list or tuple of size %d",
                peak_list->ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_float_array(last_obj, MAX_NDIM, &n, last,
                                error_msg) == CCPN_ERROR)
            || (n != peak_list->ndim))
    {
        sprintf(error_msg, "last must be float list or tuple of size %d",
                peak_list->ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (allow_aliasing_obj)
    {
        if ((get_python_int_array(allow_aliasing_obj, MAX_NDIM, &n,
                                  aa, error_msg) == CCPN_ERROR) || (n != peak_list->ndim))
        {
            sprintf(error_msg, "allow_aliasing must be int list or tuple of size %d", peak_list->ndim);
            RETURN_OBJ_ERROR(error_msg);
        }

        for (i = 0; i < peak_list->ndim; i++)
            allow_aliasing[i] = aa[i];
    }
    else
    {
        for (i = 0; i < peak_list->ndim; i++)
            allow_aliasing[i] = CCPN_TRUE;
    }

    /*
    printf("nearestPeak: xdim = %d, ydim = %d\n", xdim, ydim);
    printf("nearestPeak: xscale = %2.1f, yscale = %2.1f\n", xscale, yscale);
    printf("nearestPeak: first = %2.1f %2.1f %2.1f\n", first[0], first[1], first[2]);
    printf("nearestPeak: last = %2.1f %2.1f %2.1f\n", last[0], last[1], last[2]);
    */

    search_nearest_peak_list(peak_list, xdim, ydim, xscale, yscale,
                             first, last, allow_aliasing, &index_peak, &d2_min);

    d2_mind = d2_min;

    return Py_BuildValue("if", index_peak, d2_mind);
}

static struct PyMethodDef py_handler_methods[] =
{
    { "setIsSelected",		setIsSelected,		METH_VARARGS },
    { "setColor",		setColor,		METH_VARARGS },
    { "setSymbol",		setSymbol,		METH_VARARGS },
    { "addPeak",		addPeak,		METH_VARARGS },
    { "removePeak",		removePeak,		METH_VARARGS },
    { "removeSelectedPeaks",	removeSelectedPeaks,	METH_VARARGS },
    { "unselectSelectedPeaks",	unselectSelectedPeaks,	METH_VARARGS },
    { "findPeaks",	(PyCFunction) findPeaks,	METH_VARARGS|METH_KEYWORDS },
    { "searchPeaks",	(PyCFunction) searchPeaks,	METH_VARARGS|METH_KEYWORDS },
    { "nearestPeak",	(PyCFunction) nearestPeak,	METH_VARARGS|METH_KEYWORDS },
    { NULL,			NULL,			0 }
};

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

static PyObject *new_py_peak_list(int ndim, int *npoints)
{
    Peak_list peak_list;
    Py_Peak_list obj;

    peak_list = new_peak_list(ndim, npoints);

    if (!peak_list)
        RETURN_OBJ_ERROR("allocating Peak_list object");

    PY_MALLOC(obj, struct Py_Peak_list, &Peak_list_type);

    if (!obj)
    {
        delete_peak_list(peak_list);

        RETURN_OBJ_ERROR("allocating Py_Peak_list object");
    }

    obj->peak_list = peak_list;

    return (PyObject *) obj;
}

static void delete_py_peak_list(PyObject *self)
{
    Py_Peak_list obj = (Py_Peak_list) self;
    Peak_list peak_list = obj->peak_list;

    /*
        printf("in delete_py_peak_list\n");
    */

    delete_peak_list(peak_list);

    PY_FREE(self);
}

/*
static int print_py_peak_list(PyObject *self, FILE *fp, int flags)
{
    printf("in print_py_handler\n");

    return 0;
}
*/

static PyObject *getattr_py_peak_list(PyObject *self, char *name)
{
    /*
        Peak_list *obj = (Peak_list *) self;
        Random_access *a = obj->py_handler;

        if (equal_strings(name, "par_file"))
    	return Py_BuildValue("s", a->par_file);
        else if (equal_strings(name, "access_method"))
    	return Py_BuildValue("s", access_method_name(a->access_method));
        else if (equal_strings(name, "data_format"))
    	return get_Peak_list_format(a);
        else
    */
    return Py_FindMethod(py_handler_methods, self, name);
}

/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

/*  if implementing more...
static PySequenceMethods Peak_list_sequence_methods =
{
    Peak_list_length,
    Peak_list_concat,
    Peak_list_repeat,
    Peak_list_item,
    Peak_list_slice,
    Peak_list_ass_item,
    Peak_list_ass_slice
};

static PySequenceMethods Peak_list_sequence_methods =
{
    Peak_list_length,
    0,
    0,
    Peak_list_item,
    0,
    Peak_list_ass_item,
    0
};
*/

static PyTypeObject Peak_list_type =
{
#ifdef WIN64
    1, NULL,
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,
    "PeakList", /* name */
    sizeof(struct Py_Peak_list), /* basicsize */
    0, /* itemsize */
    delete_py_peak_list, /* destructor */
    0, /* printfunc */
    getattr_py_peak_list, /* getattr */
    0, /* setattr */
    0, /* cmpfunc */
    0, /* reprfunc */
    0, /* PyNumberMethods */
    /*&Peak_list_sequence_methods*/ /* PySequenceMethods */
};

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *init_Py_Peak_list(PyObject *self, PyObject *args)
{
    int ndim, npoints[MAX_NDIM];
    PyObject *obj, *npoints_obj;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "O", &npoints_obj))
        RETURN_OBJ_ERROR("must have one argument: npoints");

    if (get_python_int_array(npoints_obj, MAX_NDIM, &ndim, npoints,
                             error_msg) == CCPN_ERROR)
        RETURN_OBJ_ERROR("argument must be int list or tuple");

    obj = new_py_peak_list(ndim, npoints);

    return obj;
}

static PyObject *fitPeaksInRegion(PyObject *self, PyObject *args)
{
    int i, j, k, n, method, ndim, nparams, npeaks, first[MAX_NDIM], last[MAX_NDIM], dd[MAX_NDIM];
    float height, *params, *params_dev;
    Bool dim_done[MAX_NDIM];
    PyObject *peaks_obj, *block_file_obj, *first_obj, *last_obj, *peak_obj, *dim_done_obj, *fit_list, *fit_obj, *posn_obj, *lw_obj;
    Block_file block_file;
    Peak peak, *peaks;
    Line error_msg;
    CcpnStatus status;

    if (!PyArg_ParseTuple(args, "iOOOOO", &method, &peaks_obj, &block_file_obj, &first_obj, &last_obj, &dim_done_obj))
        RETURN_OBJ_ERROR("must have arguments: method, peaks, block_file, first, last, dim_done");

    if ((method != GAUSSIAN_METHOD) && (method != LORENTZIAN_METHOD))
        RETURN_ERROR_MSG("method must be GAUSSIAN_METHOD or LORENTZIAN_METHOD");

    npeaks = get_python_list_size(peaks_obj);
    if (npeaks < 0)
        RETURN_OBJ_ERROR("must have peaks list or tuple");

    if (npeaks == 0)
        RETURN_OBJ_ERROR("must have at least one peak");

    peaks = (Peak *) malloc(npeaks*sizeof(Peak));
    if (!peaks)
        RETURN_OBJ_ERROR("allocating peaks memory");

    for (i = 0; i < npeaks; i++)
    {
        peak_obj = get_python_object_by_index(peaks_obj, i);
        if (!is_py_peak(peak_obj))
        {
            FREE(peaks, Peak);
            sprintf(error_msg, "Object %d in peaks list is not a peak", i+1);
            RETURN_OBJ_ERROR(error_msg);
        }

        peak = ((Py_Peak)peak_obj)->peak;
        if (i == 0)
        {
            ndim = peak->ndim;
        }
        else if (ndim != peak->ndim)
        {
            FREE(peaks, Peak);
            sprintf(error_msg, "Peak 1 has ndim = %d, Peak %d has ndim = %d", ndim, i+1, peak->ndim);
            RETURN_OBJ_ERROR(error_msg);
        }

        peaks[i] = peak;
    }

    if (!is_py_some_block_file(block_file_obj))
    {
        FREE(peaks, Peak);
        RETURN_OBJ_ERROR("must have BlockFile object");
    }

    if (is_py_block_file(block_file_obj))
        block_file = ((Py_Block_file) block_file_obj)->block_file;
    else
        block_file = (Block_file) (((Py_Shape_block_file) block_file_obj)->shape_block_file);

    if ((get_python_int_array(first_obj, MAX_NDIM, &n, first,
                              error_msg) == CCPN_ERROR)
            || (n != ndim))
    {
        FREE(peaks, Peak);
        sprintf(error_msg, "first must be float list or tuple of size %d",
                ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_int_array(last_obj, MAX_NDIM, &n, last,
                              error_msg) == CCPN_ERROR)
            || (n != ndim))
    {
        FREE(peaks, Peak);
        sprintf(error_msg, "last must be float list or tuple of size %d",
                ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    if ((get_python_int_array(dim_done_obj, MAX_NDIM, &n,
                              dd, error_msg) == CCPN_ERROR) || (n != ndim))
    {
        FREE(peaks, Peak);
        sprintf(error_msg, "dim_done must be int list or tuple of size %d", ndim);
        RETURN_OBJ_ERROR(error_msg);
    }

    nparams = 1; /* height */
    for (i = 0; i < ndim; i++)
    {
        dim_done[i] = dd[i];
        nparams += 2; /* position, linewidth */
    }

    nparams *= npeaks;
    params = (float *) malloc(nparams*sizeof(float));
    if (!params)
    {
        FREE(peaks, Peak);
        RETURN_OBJ_ERROR("allocating params memory");
    }

    params_dev = (float *) malloc(nparams*sizeof(float));
    if (!params_dev)
    {
        FREE(params, float);
        FREE(peaks, Peak);
        RETURN_OBJ_ERROR("allocating params_dev memory");
    }

    status = fit_peaks_in_region(npeaks, peaks, method, first, last, block_file, dim_done, params, params_dev, error_msg);
    if (status == CCPN_ERROR)
    {
        FREE(params, float);
        FREE(params_dev, float);
        FREE(peaks, Peak);
        RETURN_OBJ_ERROR(error_msg);
    }

    fit_list = PyList_New(0);
    if (!fit_list)
    {
        FREE(params, float);
        FREE(params_dev, float);
        FREE(peaks, Peak);
        RETURN_OBJ_ERROR("allocating memory for fit list");
    }

    k = 0;
    for (j = 0; j < npeaks; j++)
    {
        fit_obj = PyTuple_New(3); // height, position, linewidth
        if (!fit_obj)
        {
            FREE(params, float);
            FREE(params_dev, float);
            FREE(peaks, Peak);
            RETURN_OBJ_ERROR("allocating fit data");
        }

        posn_obj = PyTuple_New(ndim);
        if (!posn_obj)
        {
            FREE(params, float);
            FREE(params_dev, float);
            FREE(peaks, Peak);
            RETURN_OBJ_ERROR("allocating position");
        }

        lw_obj = PyTuple_New(ndim);
        if (!lw_obj)
        {
            FREE(params, float);
            FREE(params_dev, float);
            FREE(peaks, Peak);
            RETURN_OBJ_ERROR("allocating linewidth");
        }

        height = params[k++];
        PyTuple_SetItem(fit_obj, 0, PyFloat_FromDouble((double) height));
        PyTuple_SetItem(fit_obj, 1, posn_obj);
        PyTuple_SetItem(fit_obj, 2, lw_obj);

        for (i = 0; i < ndim; i++)
        {
            if (dim_done[i])
                PyTuple_SetItem(posn_obj, i, PyFloat_FromDouble((double) params[k++]));
        }

        for (i = 0; i < ndim; i++)
        {
            if (dim_done[i])
                PyTuple_SetItem(lw_obj, i, PyFloat_FromDouble((double) params[k++]));
        }

        if (PyList_Append(fit_list, fit_obj))
        {
            FREE(params, float);
            FREE(params_dev, float);
            FREE(peaks, Peak);
            RETURN_OBJ_ERROR("appending fit data to fit list");
        }

        Py_DECREF(fit_obj);
    }

    FREE(params, float);
    FREE(params_dev, float);
    FREE(peaks, Peak);

    return fit_list;
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Peak_list_type_methods[] =
{
    { "PeakList",	(PyCFunction) init_Py_Peak_list,	METH_VARARGS },
    { "fitPeaksInRegion",	fitPeaksInRegion,		METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import PeakList" in
* a Python program. The function is usually called "initPeak_list": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initPeakList(void)
{
    PyObject *m, *d;

#ifdef WIN64
    Peak_list_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("PeakList", Peak_list_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("PeakList.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);

    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module PeakList");
}

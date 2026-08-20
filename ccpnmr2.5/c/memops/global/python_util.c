
/*
======================COPYRIGHT/LICENSE START==========================

python_util.c: Part of the CcpNmr Analysis program

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
#include "python_util.h"

#include "utility.h"

Bool valid_py_object(PyObject *obj, PyTypeObject *type_obj)
{
    return equal_strings(obj->ob_type->tp_name, type_obj->tp_name);
}

int get_python_list_size(PyObject *list)
{
    int size;

    if (PyList_Check(list))
	size = PyList_Size(list);
    else if (PyTuple_Check(list))
	size = PyTuple_Size(list);
    else
	size = -1;

    return size;
}

PyObject *get_python_object_by_index(PyObject *list, int ind)
{
    PyObject *obj;

    if (PyList_Check(list))
	obj = PyList_GetItem(list, ind);
    else if (PyTuple_Check(list))
	obj = PyTuple_GetItem(list, ind);
    else if (PyTuple_Check(list))
	obj = NULL;

    return obj;
}

CcpnStatus get_python_int_array(PyObject *list, int nalloc, int *n, int *value,
							CcpnString error_msg)
{
    int i, m;
    PyObject *z;
    Bool have_list = CCPN_FALSE;

    if (!list)
    {
	*n = -1;

	return CCPN_OK;
    }

    if (PyList_Check(list))
	have_list = CCPN_TRUE;
    else if (!PyTuple_Check(list))
	RETURN_ERROR_MSG("require list or tuple for int array");

    if (have_list)
	m = PyList_Size(list);
    else
	m = PyTuple_Size(list);

    if (m == 0)
    {
	*n = 0;

	return CCPN_OK;
    }

    if (nalloc < m)
	RETURN_ERROR_MSG("int array not large enough");

    for (i = 0; i < m; i++)
    {
	if (have_list)
	    z = PyList_GetItem(list, i);
	else
	    z = PyTuple_GetItem(list, i);

	if (!PyInt_Check(z))
	    RETURN_ERROR_MSG("data must be a list or tuple of integers");

	value[i] = (int) PyInt_AsLong(z);
    }

    *n = m;

    return CCPN_OK;
}

CcpnStatus get_python_int_alloc_array(PyObject *list, int *n, int **value,
							CcpnString error_msg)
{
    int i, m;
    PyObject *z;
    Bool have_list = CCPN_FALSE;

    if (!list)
    {
    	*n = -1;
    	*value = NULL;

    	return CCPN_OK;
    }

    if (PyList_Check(list))
	    have_list = CCPN_TRUE;
    else if (!PyTuple_Check(list))
	    RETURN_ERROR_MSG("require list or tuple for int array");

    if (have_list)
	    m = PyList_Size(list);
    else
	    m = PyTuple_Size(list);

    if (m == 0)
    {
    	*n = 0;
    	*value = NULL;

    	return CCPN_OK;
    }

    sprintf(error_msg, "allocating int array memory");
    MALLOC(*value, int, m);

    for (i = 0; i < m; i++)
    {
    	if (have_list)
    	    z = PyList_GetItem(list, i);
    	else
    	    z = PyTuple_GetItem(list, i);

    	if (!PyInt_Check(z))
    	{
    	    FREE(*value, int);

    	    RETURN_ERROR_MSG("data must be a list or tuple of integers");
    	}

    	(*value)[i] = (int) PyInt_AsLong(z);
    }

    *n = m;

    return CCPN_OK;
}

CcpnStatus get_python_float_array(PyObject *list, int nalloc, int *n, float *value,
							CcpnString error_msg)
{
    int i, m;
    PyObject *z;
    Bool have_list = CCPN_FALSE;

    if (!list)
    {
	*n = -1;

	return CCPN_OK;
    }

    if (PyList_Check(list))
	have_list = CCPN_TRUE;
    else if (!PyTuple_Check(list))
	RETURN_ERROR_MSG("require list or tuple for float array");

    if (have_list)
	m = PyList_Size(list);
    else
	m = PyTuple_Size(list);

    if (m == 0)
    {
	*n = 0;

	return CCPN_OK;
    }

    if (nalloc < m)
	RETURN_ERROR_MSG("float array not large enough");

    for (i = 0; i < m; i++)
    {
	if (have_list)
	    z = PyList_GetItem(list, i);
	else
	    z = PyTuple_GetItem(list, i);

	if (!PyFloat_Check(z) && !PyInt_Check(z))
	    RETURN_ERROR_MSG("data must be a list or tuple of floats");

	if (PyFloat_Check(z))
	    value[i] = (float) PyFloat_AsDouble(z);
	else
	    value[i] = (float) PyInt_AsLong(z);
    }

    *n = m;

    return CCPN_OK;
}

CcpnStatus get_python_float_alloc_array(PyObject *list, int *n, float **value,
							CcpnString error_msg)
{
    int i, m;
    PyObject *z;
    Bool have_list = CCPN_FALSE;

    if (!list)
    {
	*n = -1;
	*value = NULL;

	return CCPN_OK;
    }

    if (PyList_Check(list))
	have_list = CCPN_TRUE;
    else if (!PyTuple_Check(list))
	RETURN_ERROR_MSG("require list or tuple for float array");

    if (have_list)
	m = PyList_Size(list);
    else
	m = PyTuple_Size(list);

    if (m == 0)
    {
	*n = 0;
	*value = NULL;

	return CCPN_OK;
    }

    sprintf(error_msg, "allocating float array memory");
    MALLOC(*value, float, m);

    for (i = 0; i < m; i++)
    {
	if (have_list)
	    z = PyList_GetItem(list, i);
	else
	    z = PyTuple_GetItem(list, i);

	if (!PyFloat_Check(z) && !PyInt_Check(z))
	{
	    FREE(*value, float);

	    RETURN_ERROR_MSG("data must be a list or tuple of floats");
	}

	if (PyFloat_Check(z))
	    (*value)[i] = (float) PyFloat_AsDouble(z);
	else
	    (*value)[i] = (float) PyInt_AsLong(z);
    }

    *n = m;

    return CCPN_OK;
}

CcpnStatus get_python_float_alloc_matrix(PyObject *list, int *n1, int *n2,
					float ***value, CcpnString error_msg)
{
    int i, m, n;
    PyObject *z;
    Bool have_list = CCPN_FALSE;

    if (!list)
    {
	*n1 = *n2 = -1;
	*value = NULL;

	return CCPN_OK;
    }

    if (PyList_Check(list))
	have_list = CCPN_TRUE;
    else if (!PyTuple_Check(list))
	RETURN_ERROR_MSG("require list or tuple for float matrix");

    if (have_list)
	m = PyList_Size(list);
    else
	m = PyTuple_Size(list);

    if (m == 0)
    {
	*n1 = *n2 = 0;
	*value = NULL;

	return CCPN_OK;
    }

    sprintf(error_msg, "allocating float matrix memory");
    MALLOC(*value, float *, m);

    for (i = 0; i < m; i++)
    {
	if (have_list)
	    z = PyList_GetItem(list, i);
	else
	    z = PyTuple_GetItem(list, i);

	CHECK_STATUS(get_python_float_alloc_array(z, &n, (*value)+i, error_msg));

	if (i == 0)
	{
	    *n2 = n;
	}
	else if (n != *n2)
	{
	    sprintf(error_msg, "inconsistent second matrix dimension %d and %d",
								*n2, n);
	    return CCPN_ERROR;
	}
    }

    *n1 = m;

    return CCPN_OK;
}

CcpnStatus get_python_float_alloc_cube(PyObject *list, int *n1, int *n2, int *n3,
					float ****value, CcpnString error_msg)
{
    int i, m, nn2, nn3;
    PyObject *z;
    Bool have_list = CCPN_FALSE;

    if (!list)
    {
	*n1 = *n2 = *n3 = -1;
	*value = NULL;

	return CCPN_OK;
    }

    if (PyList_Check(list))
	have_list = CCPN_TRUE;
    else if (!PyTuple_Check(list))
	RETURN_ERROR_MSG("require list or tuple for float cube");

    if (have_list)
	m = PyList_Size(list);
    else
	m = PyTuple_Size(list);

    if (m == 0)
    {
	*n1 = *n2 = *n3 = 0;
	*value = NULL;

	return CCPN_OK;
    }

    sprintf(error_msg, "allocating float cube memory");
    MALLOC(*value, float **, m);

    for (i = 0; i < m; i++)
    {
	if (have_list)
	    z = PyList_GetItem(list, i);
	else
	    z = PyTuple_GetItem(list, i);

	CHECK_STATUS(get_python_float_alloc_matrix(z, &nn2, &nn3, (*value)+i, error_msg));

	if (i == 0)
	{
	    *n2 = nn2;
	    *n3 = nn3;
	}
	else if ((nn2 != *n2) || (nn3 != *n3))
	{
	    sprintf(error_msg, "inconsistent second and third cube dimension (%d, %d) versus (%d, %d)",
								*n2, *n3, nn2, nn3);
	    return CCPN_ERROR;
	}
    }

    *n1 = m;

    return CCPN_OK;
}

CcpnStatus get_python_double_array(PyObject *list, int nalloc, int *n, double *value,
							CcpnString error_msg)
{
    int i, m;
    PyObject *z;
    Bool have_list = CCPN_FALSE;

    if (!list)
    {
	*n = -1;

	return CCPN_OK;
    }

    if (PyList_Check(list))
	have_list = CCPN_TRUE;
    else if (!PyTuple_Check(list))
	RETURN_ERROR_MSG("require list or tuple for int array");

    if (have_list)
	m = PyList_Size(list);
    else
	m = PyTuple_Size(list);

    if (m == 0)
    {
	*n = 0;

	return CCPN_OK;
    }

    if (nalloc < m)
	RETURN_ERROR_MSG("double array not large enough");

   for (i = 0; i < m; i++)
    {
	if (have_list)
	    z = PyList_GetItem(list, i);
	else
	    z = PyTuple_GetItem(list, i);

	if (!PyFloat_Check(z) && !PyInt_Check(z))
	    RETURN_ERROR_MSG("data must be a list or tuple of floats");

	if (PyFloat_Check(z))
	    value[i] = (double) PyFloat_AsDouble(z);
	else
	    value[i] = (double) PyInt_AsLong(z);
    }

    *n = m;

    return CCPN_OK;
}

CcpnStatus get_python_double_alloc_array(PyObject *list, int *n, double **value,
							CcpnString error_msg)
{
    int i, m;
    PyObject *z;
    Bool have_list = CCPN_FALSE;

    if (!list)
    {
    	*n = -1;
    	*value = NULL;

    	return CCPN_OK;
    }

    if (PyList_Check(list))
	    have_list = CCPN_TRUE;
    else if (!PyTuple_Check(list))
	    RETURN_ERROR_MSG("require list or tuple for double array");

    if (have_list)
	    m = PyList_Size(list);
    else
	    m = PyTuple_Size(list);

    if (m == 0)
    {
    	*n = 0;
    	*value = NULL;

    	return CCPN_OK;
    }

    sprintf(error_msg, "allocating double array memory");
    MALLOC(*value, double, m);

    for (i = 0; i < m; i++)
    {
    	if (have_list)
    	    z = PyList_GetItem(list, i);
    	else
    	    z = PyTuple_GetItem(list, i);

    	if (!PyFloat_Check(z) && !PyInt_Check(z))
    	{
    	    FREE(*value, double);

    	    RETURN_ERROR_MSG("data must be a list or tuple of floats");
    	}

    	if (PyFloat_Check(z))
    	    (*value)[i] = (double) PyFloat_AsDouble(z);
    	else
    	    (*value)[i] = (double) PyInt_AsLong(z);
    }

    *n = m;

    return CCPN_OK;
}

CcpnStatus get_python_double_alloc_matrix(PyObject *list, int *n1, int *n2,
					double ***value, CcpnString error_msg)
{
    int i, m, n;
    PyObject *z;
    Bool have_list = CCPN_FALSE;

    if (!list)
    {
	*n1 = -1;
	*value = NULL;

	return CCPN_OK;
    }

    if (PyList_Check(list))
	have_list = CCPN_TRUE;
    else if (!PyTuple_Check(list))
	RETURN_ERROR_MSG("require list or tuple for double matrix");

    if (have_list)
	m = PyList_Size(list);
    else
	m = PyTuple_Size(list);

    if (m == 0)
    {
	*n1 = 0;
	*value = NULL;

	return CCPN_OK;
    }

    sprintf(error_msg, "allocating double matrix memory");
    MALLOC(*value, double *, m);

    for (i = 0; i < m; i++)
    {
	if (have_list)
	    z = PyList_GetItem(list, i);
	else
	    z = PyTuple_GetItem(list, i);

	CHECK_STATUS(get_python_double_alloc_array(z, &n, (*value)+i, error_msg));

	if (i == 0)
	{
	    *n2 = n;
	}
	else if (n != *n2)
	{
	    sprintf(error_msg, "inconsistent second matrix dimension %d and %d",
								*n2, n);
	    return CCPN_ERROR;
	}
    }

    *n1 = m;

    return CCPN_OK;
}

PyObject *get_python_int_list(int n, int *array)
{
    int i;
    PyObject *list, *z;

    list = PyList_New(n);
    if (!list)
      return NULL;

    for (i = 0; i < n; i++)
    {
	z = PyInt_FromLong((long) array[i]);

	PyList_SetItem(list, i, z);
    }

    return list;
}

PyObject *get_python_float_list(int n, float *array)
{
    int i;
    PyObject *list, *z;

    list = PyList_New(n);
    if (!list)
      return NULL;

    for (i = 0; i < n; i++)
    {
	z = PyFloat_FromDouble((double) array[i]);

	PyList_SetItem(list, i, z);
    }

    return list;
}

PyObject *get_python_double_list(int n, double *array)
{
    int i;
    PyObject *list, *z;

    list = PyList_New(n);
    if (!list)
      return NULL;

    for (i = 0; i < n; i++)
    {
	z = PyFloat_FromDouble(array[i]);

	PyList_SetItem(list, i, z);
    }

    return list;
}


PyObject *get_python_float_matrix(int n1, int n2, float **matrix)
{
    int i, j;
    PyObject *mat, *list, *z;

    mat = PyList_New(n1);
    if (!mat)
      return NULL;

    for (i = 0; i < n1; i++)
    {
      list = PyList_New(n2);
      if (!list)
        return NULL;

      PyList_SetItem(mat, i, list);

      for (j = 0; j < n2; j++)
      {
	z = PyFloat_FromDouble((double) matrix[i][j]);

	PyList_SetItem(list, j, z);
      }
    }

    return mat;
}

CcpnStatus set_python_float_array(PyObject *obj, int n, float *array,
							CcpnString error_msg)
{
    int i, m;
    PyObject *z;

    if (!PyList_Check(obj))
	RETURN_ERROR_MSG("require list for float array");

    m = PyList_Size(obj);

    if (m != n)
    {
	sprintf(error_msg,
	    "inconsistent size for float array, found %d, should be %d", m, n);
	return CCPN_ERROR;
    }

    for (i = 0; i < n; i++)
    {
	z = PyFloat_FromDouble((double) array[i]);

	PyList_SetItem(obj, i, z);
    }

    return CCPN_OK;
}

CcpnStatus create_python_float_array(PyObject **obj, int n, float *array,
							CcpnString error_msg)
{
    *obj = PyList_New(n);

    if (*obj == NULL)
    {
	sprintf(error_msg, "creating float array object");
	return CCPN_ERROR;
    }

    return set_python_float_array(*obj, n, array, error_msg);
}

CcpnStatus set_python_double_array(PyObject *obj, int n, double *array,
							CcpnString error_msg)
{
    int i, m;
    PyObject *z;

    if (!PyList_Check(obj))
	RETURN_ERROR_MSG("require list for double array");

    m = PyList_Size(obj);

    if (m != n)
    {
	sprintf(error_msg,
	    "inconsistent size for double array, found %d, should be %d", m, n);
	return CCPN_ERROR;
    }

    for (i = 0; i < n; i++)
    {
	z = PyFloat_FromDouble((double) array[i]);

	PyList_SetItem(obj, i, z);
    }

    return CCPN_OK;
}

CcpnStatus set_python_float_matrix(PyObject *obj, int n1, int n2,
				float **mat, CcpnString error_msg)
{
    int i, m;
    Bool have_list = CCPN_FALSE;
    PyObject *z;

    if (PyList_Check(obj))
	have_list = CCPN_TRUE;
    else if (!PyTuple_Check(obj))
	RETURN_ERROR_MSG("require list or tuple for float matrix");

    if (have_list)
	m = PyList_Size(obj);
    else
	m = PyTuple_Size(obj);

    if (m != n1)
    {
	sprintf(error_msg,
	    "inconsistent size for float matrix, found %d, should be %d", m, n1);
	return CCPN_ERROR;
    }

    for (i = 0; i < n1; i++)
    {
	if (have_list)
	    z = PyList_GetItem(obj, i);
	else
	    z = PyTuple_GetItem(obj, i);

	CHECK_STATUS(set_python_float_array(z, n2, mat[i], error_msg));
    }

    return CCPN_OK;
}

CcpnStatus set_python_double_matrix(PyObject *obj, int n1, int n2,
				double **mat, CcpnString error_msg)
{
    int i, m;
    Bool have_list = CCPN_FALSE;
    PyObject *z;

    if (PyList_Check(obj))
	have_list = CCPN_TRUE;
    else if (!PyTuple_Check(obj))
	RETURN_ERROR_MSG("require list or tuple for double matrix");

    if (have_list)
	m = PyList_Size(obj);
    else
	m = PyTuple_Size(obj);

    if (m != n1)
    {
	sprintf(error_msg,
	    "inconsistent size for double matrix, found %d, should be %d", m, n1);
	return CCPN_ERROR;
    }

    for (i = 0; i < n1; i++)
    {
	if (have_list)
	    z = PyList_GetItem(obj, i);
	else
	    z = PyTuple_GetItem(obj, i);

	CHECK_STATUS(set_python_double_array(z, n2, mat[i], error_msg));
    }

    return CCPN_OK;
}

CcpnStatus set_python_float_cube(PyObject *obj, int n1, int n2, int n3,
				float ***cube, CcpnString error_msg)
{
    int i, m;
    Bool have_list = CCPN_FALSE;
    PyObject *z;

    if (PyList_Check(obj))
	have_list = CCPN_TRUE;
    else if (!PyTuple_Check(obj))
	RETURN_ERROR_MSG("require list or tuple for float cube");

    if (have_list)
	m = PyList_Size(obj);
    else
	m = PyTuple_Size(obj);

    if (m != n1)
    {
	sprintf(error_msg,
	    "inconsistent size for float cube, found %d, should be %d", m, n1);
	return CCPN_ERROR;
    }

    for (i = 0; i < n1; i++)
    {
	if (have_list)
	    z = PyList_GetItem(obj, i);
	else
	    z = PyTuple_GetItem(obj, i);

	CHECK_STATUS(set_python_float_matrix(z, n2, n3, cube[i], error_msg));
    }

    return CCPN_OK;
}


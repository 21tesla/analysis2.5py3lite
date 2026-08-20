
/*
======================COPYRIGHT/LICENSE START==========================

python_util.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_python_util
#define _incl_python_util

#include "Python.h"

#include "macros.h"
#include "types.h"

#define RETURN_OBJ_ERROR(message) \
	{ PyErr_SetString(ErrorObject, message); return NULL; }

#define RETURN_INT_ERROR(message) \
	{ PyErr_SetString(ErrorObject, message); return -1; }

#define CHECK_OBJ_STATUS(status) \
	{   if ((status) == CCPN_ERROR)  RETURN_OBJ_ERROR(error_msg);   }

#define PY_MALLOC(obj, type, typeobj) \
	{   obj = PyObject_NEW(type, typeobj); \
            DEBUG_CODE_HEX("py_malloc", obj);   }

#define PY_FREE(obj) \
	{   DEBUG_CODE_HEX("py_free", obj); \
            PyObject_DEL(obj);   }

#ifdef WIN64
#define PY_MOD_INIT_FUNC PyMODINIT_FUNC
#else
#define PY_MOD_INIT_FUNC PyMODINIT_FUNC
#endif
            
extern Bool valid_py_object
	(PyObject *obj, PyTypeObject *type_obj);

extern int get_python_list_size
	(PyObject *list);
 
/* get ind element of list */
extern PyObject *get_python_object_by_index
	(PyObject *list, int ind);

extern CcpnStatus get_python_int_array
	(PyObject *list, int nalloc, int *n, int *value, CcpnString error_msg);

extern CcpnStatus get_python_int_alloc_array
	(PyObject *list, int *n, int **value, CcpnString error_msg);

extern CcpnStatus get_python_float_array
	(PyObject *list, int nalloc, int *n, float *value, CcpnString error_msg);

extern CcpnStatus get_python_float_alloc_array
	(PyObject *list, int *n, float **value, CcpnString error_msg);

extern CcpnStatus get_python_float_alloc_matrix
	(PyObject *list, int *n1, int *n2, float ***value, CcpnString error_msg);

extern CcpnStatus get_python_float_alloc_cube
	(PyObject *list, int *n1, int *n2, int *n3, float ****value, CcpnString error_msg);

extern CcpnStatus get_python_double_array
	(PyObject *list, int nalloc, int *n, double *value, CcpnString error_msg);

extern CcpnStatus get_python_double_alloc_array
	(PyObject *list, int *n, double **value, CcpnString error_msg);

extern CcpnStatus get_python_double_alloc_matrix
	(PyObject *list, int *n1, int *n2, double ***value, CcpnString error_msg);

extern PyObject *get_python_int_list
	(int n, int *array);

extern PyObject *get_python_float_list
	(int n, float *array);

extern PyObject *get_python_double_list
	(int n, double *array);

extern PyObject *get_python_float_matrix
	(int n1, int n2, float **matrix);

extern CcpnStatus set_python_float_array
	(PyObject *obj, int n, float *array, CcpnString error_msg);

extern CcpnStatus create_python_float_array
	(PyObject **obj, int n, float *array, CcpnString error_msg);

extern CcpnStatus set_python_double_array
	(PyObject *obj, int n, double *array, CcpnString error_msg);

extern CcpnStatus set_python_float_matrix
	(PyObject *obj, int n1, int n2, float **mat, CcpnString error_msg);

extern CcpnStatus set_python_double_matrix
	(PyObject *obj, int n1, int n2, double **mat, CcpnString error_msg);

extern CcpnStatus set_python_float_cube
	(PyObject *list, int n1, int n2, int n3, float ***cube, CcpnString error_msg);

#endif /* _incl_python_util */

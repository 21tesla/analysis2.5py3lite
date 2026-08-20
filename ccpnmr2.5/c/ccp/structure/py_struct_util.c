
/*
======================COPYRIGHT/LICENSE START==========================

py_struct_util.c: Part of the CcpNmr Analysis program

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
#include "py_struct_util.h"

#include "struct_util.h"
#include "python_util.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * TYPE INFORMATION
 *****************************************************************************/

/*****************************************************************************
 * MISCELLANEOUS METHODS
 *****************************************************************************/

/*****************************************************************************
 * INSTANCE METHODS
 *****************************************************************************/

/*****************************************************************************
 * BASIC TYPE-OPERATIONS
 *****************************************************************************/

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static void free_coords_memory(int natoms, float **coords1, float **coords2, float *weight)
{
    int i;

    for (i = 0; i < natoms; i++)
    {
	if (coords1)
	    FREE(coords1[i], float);

	if (coords2)
	    FREE(coords2[i], float);
    }

    FREE(coords1, float *);
    FREE(coords2, float *);

    FREE(weight, float);
}

/* Data comes in as follows:
     tuple/list of coords1_obj, coords2_obj of length natoms
     for each coords obj have list [x, y, z]
   Note that the coords1 coordinates must be in a list (they are modified in place)
*/
static PyObject *py_align_coordinates(PyObject *self, PyObject *args)
{
    int i, natoms, n1, n2;
    float **coords1 = NULL, **coords2 = NULL, *weight = NULL, err;
    PyObject *obj, *rot_obj, *err_obj, *coords1_obj, *coords2_obj, *weight_obj = NULL;
    static float r0[COORD_NDIMS], r1[COORD_NDIMS], r2[COORD_NDIMS];
    float *rotation[] = {r0, r1, r2};
    CcpnStatus status;
    Line error_msg;

    if (!PyArg_ParseTuple(args, "OO|O", &coords1_obj, &coords2_obj, &weight_obj))
        RETURN_OBJ_ERROR("must have arguments: coords1 coords2 [ weight ]");

    if ((get_python_float_alloc_matrix(coords1_obj, &n1, &n2, &coords1,
                                                error_msg) == CCPN_ERROR)
            || (n2 != COORD_NDIMS))
    {
        sprintf(error_msg,
            "coords1 must be matrix of size n x %d for some n", COORD_NDIMS);
	free_coords_memory(natoms, coords1, coords2, weight);
        RETURN_OBJ_ERROR(error_msg);
    }

    natoms = n1;

    if ((get_python_float_alloc_matrix(coords2_obj, &n1, &n2, &coords2,
                                                error_msg) == CCPN_ERROR)
            || (n1 != natoms) || (n2 != COORD_NDIMS))
    {
        sprintf(error_msg,
            "coords2 must be matrix of size %d x %d (first dim size from coords1)",
            natoms, COORD_NDIMS);
	free_coords_memory(natoms, coords1, coords2, weight);
        RETURN_OBJ_ERROR(error_msg);
    }

    if (weight_obj)
    {
        if ((get_python_float_alloc_array(weight_obj, &n1, &weight,
                                                error_msg) == CCPN_ERROR)
            || (n1 != natoms))
        {
            sprintf(error_msg,
                "weight must be array of size %d (from coords1, coords2)",
                natoms);
	    free_coords_memory(natoms, coords1, coords2, weight);
            RETURN_OBJ_ERROR(error_msg);
        }
    }

    status = align_coordinates(natoms, coords1, coords2, weight, rotation,
                                &err, error_msg);

    if (status == CCPN_OK)
	status = set_python_float_matrix(coords1_obj, natoms, COORD_NDIMS, coords1, error_msg);

    free_coords_memory(natoms, coords1, coords2, weight);

    if (status == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    rot_obj = get_python_float_matrix(COORD_NDIMS, COORD_NDIMS, rotation);
    err_obj = Py_BuildValue("f", err);
    obj = PyTuple_New(2);
    PyTuple_SetItem(obj, 0, rot_obj);
    PyTuple_SetItem(obj, 1, err_obj);

    return obj;
}

static CcpnStatus alloc_ensemble_memory(int nensembles, int natoms, float **weight,
		float **atom_rmsd, float **ensemble_rmsd, CcpnString error_msg)
{
    sprintf(error_msg, "allocating ensemble memory");

    MALLOC(*weight, float, natoms);
    MALLOC(*atom_rmsd, float, natoms);
    MALLOC(*ensemble_rmsd, float, nensembles);

    return CCPN_OK;
}

static void free_ensemble_memory(int nensembles, int natoms, float ***coords,
		float *weight, float *atom_rmsd, float *ensemble_rmsd)
{
    int i, j;

    if (coords)
    {
	for (i = 0; i < nensembles; i++)
	{
	    if (coords[i])
	    {
		for (j = 0; j < natoms; j++)
		    FREE(coords[i][j], float);

		FREE(coords[i], float*);
	    }
	}

	FREE(coords, float**);
    }

    FREE(weight, float);
    FREE(atom_rmsd, float);
    FREE(ensemble_rmsd, float);
}

static PyObject *py_align_ensemble(PyObject *self, PyObject *args)
{
    int i, m, n, natoms, nensembles;
    float ***coords, *weight, *atom_rmsd, *ensemble_rmsd, toterr;
    PyObject *obj, *weight_obj = NULL, *ensemble_obj, *coords_obj;
    PyObject *atom_rmsd_obj, *ensemble_rmsd_obj, *err_obj;
    CcpnStatus status;
    Line error_msg;

/* Data comes in as follows:
     tuple/list of ensemble objs of length nensembles
     for each ensemble obj have tuple/list ([x1, y1, z1], ..., [xn, yn, zn])
     where n = natoms
   Note that the coordinates must be in a list (they are modified in place)
   There is also an optional weight passed in, of length natoms
*/

    if (!PyArg_ParseTuple(args, "O|O", &ensemble_obj, &weight_obj))
        RETURN_OBJ_ERROR("must have arguments: ensemble [ weight ]");

    if (get_python_float_alloc_cube(ensemble_obj, &nensembles, &natoms, &n, &coords, error_msg) == CCPN_ERROR)
       	RETURN_OBJ_ERROR(error_msg);

    if (n != COORD_NDIMS)
    {
	sprintf(error_msg, "coordinates each of size %d, must be of size %d", n, COORD_NDIMS);
       	RETURN_OBJ_ERROR(error_msg);
    }

    if (alloc_ensemble_memory(nensembles, natoms, &weight,
			&atom_rmsd, &ensemble_rmsd, error_msg) == CCPN_ERROR)
    {
	free_ensemble_memory(nensembles, natoms, coords, weight, atom_rmsd, ensemble_rmsd);
       	RETURN_OBJ_ERROR(error_msg);
    }

    if (weight_obj)
    {
	if (get_python_float_array(weight_obj, natoms, &n, weight, error_msg) == CCPN_ERROR)
       	RETURN_OBJ_ERROR(error_msg);

	if (n != natoms)
	{
	    sprintf(error_msg, "weight is of size %d, must be of size %d", n, natoms);
	    free_ensemble_memory(nensembles, natoms, coords, weight, atom_rmsd, ensemble_rmsd);
       	    RETURN_OBJ_ERROR(error_msg);
	}
    }
    else
    {
	for (i = 0; i < natoms; i++)
	    weight[i] = 1;
    }

    status = align_ensemble(nensembles, natoms, coords, weight, atom_rmsd, ensemble_rmsd,
                                &toterr, error_msg);

    if (status == CCPN_OK)
	status = set_python_float_cube(ensemble_obj, nensembles, natoms, COORD_NDIMS, coords, error_msg);

    if (status == CCPN_OK)
	status = create_python_float_array(&atom_rmsd_obj, natoms, atom_rmsd, error_msg);

    if (status == CCPN_OK)
	status = create_python_float_array(&ensemble_rmsd_obj, nensembles, ensemble_rmsd, error_msg);

    free_ensemble_memory(nensembles, natoms, coords, weight, atom_rmsd, ensemble_rmsd);

    if (status == CCPN_ERROR)
        RETURN_OBJ_ERROR(error_msg);

    err_obj = Py_BuildValue("f", toterr);
    obj = PyTuple_New(3);
    PyTuple_SetItem(obj, 0, err_obj);
    PyTuple_SetItem(obj, 1, atom_rmsd_obj);
    PyTuple_SetItem(obj, 2, ensemble_rmsd_obj);

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
    { "alignCoordinates",	(PyCFunction) py_align_coordinates,	METH_VARARGS },
    { "alignEnsemble",		(PyCFunction) py_align_ensemble,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import StructUtil" in 
* a Python program. The function is usually called "initStructUtil": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initStructUtil(void)
{
    PyObject *m, *d;

    /* create the module and add the functions */
    m = Py_InitModule("StructUtil", Structure_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("StructUtil.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module StructUtil");
}

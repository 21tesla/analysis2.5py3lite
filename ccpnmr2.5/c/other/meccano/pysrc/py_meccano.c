
/*
======================COPYRIGHT/LICENSE START==========================

py_meccano.c: Python/C interface to Meccano code

Copyright (C) 2008 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in ../../../license/CCPN.license.

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
#include "meccano.h"

#include "python_util.h"

static PyObject *ErrorObject;   /* locally-raised exception */

/*****************************************************************************
 * INTERNAL CODE
 *****************************************************************************/

static PyTypeObject Meccano_type;

static PyObject *getCoordTuple(char *atomType, int resNb, gsl_vector *v)
{
    return Py_BuildValue("(siddd)", atomType, resNb, XYZ(v));
}

static PyObject *createCoordinates(polyPeptide_t * polyP)
{
    int             i = 1, l, resNb;
    int             windowSize = polyP->windowSize;
    int             windowOffset = polyP->windowOffset;
    acAm_t          *acAm;
    PyObject        *result;

    result = PyList_New(0);

    for (l = windowOffset; l <= windowOffset + windowSize; l++) {
	resNb = polyP->acAm[l].resNb;
	acAm = polyP->acAm + l;
        if (acAm->n_Flag)
            PyList_Append(result, getCoordTuple("N", resNb, acAm->n_));

        if (polyP->acAm[l].h_Flag)
            PyList_Append(result, getCoordTuple("H", resNb, acAm->h_));

        if (polyP->acAm[l].caFlag)
            PyList_Append(result, getCoordTuple("CA", resNb, acAm->ca));

        if (polyP->acAm[l].haFlag)
            PyList_Append(result, getCoordTuple("HA", resNb, acAm->ha));

        if (polyP->acAm[l].cbFlag)
            PyList_Append(result, getCoordTuple("CB", resNb, acAm->cb));

        if (polyP->acAm[l].c_Flag) 
            PyList_Append(result, getCoordTuple("C", resNb, acAm->c_));

        if (polyP->acAm[l].o_Flag)
            PyList_Append(result, getCoordTuple("O", resNb, acAm->o_));
    }

    return result;
}

static PyObject *runFwdParams(param_t *params, int firstPPlaneFrag, int lastPPlaneFrag,
	int ppNbMin, int minValueBest, int maxValueBest, char *ramaDBFileName,
	PyObject *seqData, PyObject *mediumData, PyObject *phipsiData,
        PyObject *hbscData, char *error_msg)
{
    int i, m, n, seq_num, mediaNb, num, num1, num2;
    int fPepPla = 1000000, lPepPla = 0, totalPepPlaNb, firstPepPla;
    char *seq_name, *atom1, *atom2, *tensorName;
    double Aa, Ar, alpha, beta, gamma, r2, r3, r4, r5, y, sigma;
    PyObject *tdObj, *tsObj, *rObj, *obj, *result;

    params->ramaDB = ReadRamaDBFile(ramaDBFileName);
    if (!params->ramaDB)
    {
	sprintf(error_msg, "opening ramaDBFileName '%s' for reading", ramaDBFileName);
	RETURN_OBJ_ERROR(error_msg);
    }

    mediaNb = get_python_list_size(mediumData);
    if (mediaNb < 1)
	RETURN_OBJ_ERROR("mediumData must be list or tuple of size at least 1");

    for (m = 0; m < mediaNb; m++) {
	obj = get_python_object_by_index(mediumData, m);
        if (!PyArg_ParseTuple(obj, "OOO", &tdObj, &tsObj, &rObj))
	{
	    sprintf(error_msg, "Element %d of mediumData must be tuple (tensorDyn, tensorStat, rdc)", m);
            RETURN_OBJ_ERROR(error_msg);
	}

        n = get_python_list_size(rObj);
        if (n < 0)
	{
	    sprintf(error_msg, "Third part of element %d of mediumData must be rdcData list or tuple of size at least 1", m);
            RETURN_OBJ_ERROR(error_msg);
	}

	for (i = 0; i < n; i++)
	{
	    obj = get_python_object_by_index(rObj, i);
            if (!PyArg_ParseTuple(obj, "isisdd", &num1, &atom1, &num2, &atom2, &y, &sigma))
	    {
	        sprintf(error_msg, "Element %d of third part of element %d of mediumData must be rdcData tuple (num1, atom1, num2, atom2, y, sigma)", i, m);
                RETURN_OBJ_ERROR(error_msg);
	    }

	    UpdatePepPlaNbDet(num1, atom1, num2, atom2, &fPepPla, &lPepPla);
	}
    }

    totalPepPlaNb = lPepPla - fPepPla + 1;
    firstPepPla = fPepPla;

    params->pPlaneRef = InitPepPlaFwd();
    if (!params->pPlaneRef)
	RETURN_OBJ_ERROR("allocating Meccano pPlaneRef");

    params->polyP = InitPolyPeptide(totalPepPlaNb, firstPepPla);
    if (!params->polyP)
	RETURN_OBJ_ERROR("allocating Meccano polyP");

    params->data = InitData(mediaNb, totalPepPlaNb, firstPepPla);
    if (!params->data)
	RETURN_OBJ_ERROR("allocating Meccano data");

    params->tensorStat = InitTensor(mediaNb);
    if (!params->tensorStat)
	RETURN_OBJ_ERROR("allocating Meccano tensorStat");

    params->tensorDyn = InitTensor(mediaNb);
    if (!params->tensorDyn)
	RETURN_OBJ_ERROR("allocating Meccano tensorDyn");

    n = get_python_list_size(seqData);
    if (n < 0)
	RETURN_OBJ_ERROR("seqData must be list or tuple");

    for (i = 0; i < n; i++)
    {
	obj = get_python_object_by_index(seqData, i);
        if (!PyArg_ParseTuple(obj, "is", &seq_num, &seq_name))
	{
	    sprintf(error_msg, "Element %d of seqData must be tuple (seq_num, seq_name)", i);
            RETURN_OBJ_ERROR(error_msg);
	}

	SetSequence(params->polyP, seq_num, seq_name);
    }

    for (m = 0; m < mediaNb; m++) {
	obj = get_python_object_by_index(mediumData, m);
        PyArg_ParseTuple(obj, "OOO", &tdObj, &tsObj, &rObj);

        if (!PyArg_ParseTuple(tdObj, "sddddd", &tensorName, &Aa, &Ar, &alpha, &beta, &gamma))
	{
	    sprintf(error_msg, "First part of element %d of mediumData must be tuple tensorDyn (tensorName, Aa, Ar, alpha, beta, gamma)", i);
            RETURN_OBJ_ERROR(error_msg);
	}

	SetTensor(tensorName, Aa, Ar, alpha, beta, gamma, params->tensorDyn+m);

        if (!PyArg_ParseTuple(tsObj, "sddddd", &tensorName, &Aa, &Ar, &alpha, &beta, &gamma))
	{
	    sprintf(error_msg, "Second part of element %d of mediumData must be tuple tensorStat (tensorName, Aa, Ar, alpha, beta, gamma)", m);
            RETURN_OBJ_ERROR(error_msg);
	}

	SetTensor(tensorName, Aa, Ar, alpha, beta, gamma, params->tensorStat+m);

        n = get_python_list_size(rObj);
	for (i = 0; i < n; i++)
	{
	    obj = get_python_object_by_index(rObj, i);
            PyArg_ParseTuple(obj, "isisdd", &num1, &atom1, &num2, &atom2, &y, &sigma);
	    SetRdcData(m, num1, atom1, num2, atom2, y, sigma, params->data);
	}
    }

    n = get_python_list_size(phipsiData);
    if (n < 0)
	RETURN_OBJ_ERROR("phipsiData must be list or tuple");

    for (i = 0; i < n; i++)
    {
	obj = get_python_object_by_index(phipsiData, i);
        if (!PyArg_ParseTuple(obj, "idddd", &num, &r2, &r3, &r4, &r5))
	{
	    sprintf(error_msg, "Element %d of phipsiData must be tuple (num, r2, r3, r4, r5)", i);
            RETURN_OBJ_ERROR(error_msg);
	}

	SetPhiPsi(num, r2, r3, r4, r5, params->data);
    }

    if (hbscData)
    {
        n = get_python_list_size(hbscData);
        if (n < 0)
	{
	    sprintf(error_msg, "hbscData must be list or tuple");
            RETURN_OBJ_ERROR(error_msg);
	}

	for (i = 0; i < n; i++)
	{
	    obj = get_python_object_by_index(hbscData, i);
            PyArg_ParseTuple(obj, "isisdd", &num1, &atom1, &num2, &atom2, &y, &sigma);
	    SetHbscData(num1, atom1, num2, atom2, y, sigma, params->data);
	}
    }

    RunMeccanoFromParams(firstPPlaneFrag, lastPPlaneFrag, ppNbMin,
        minValueBest, maxValueBest, params);

    result = createCoordinates(params->polyP);

    return result;
}

/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static PyObject *runFwd(PyObject *self, PyObject *args)
{
    int firstPPlaneFrag, lastPPlaneFrag, ppNbMin, minValueBest, maxValueBest;
    char *ramaDBFileName;
    PyObject *seqData, *mediumData, *phipsiData, *hbscData=NULL, *result;
    param_t *params;
    Long_line error_msg;

    if (!PyArg_ParseTuple(args, "iiiiisOOO|O", &firstPPlaneFrag, &lastPPlaneFrag, &ppNbMin, &minValueBest, &maxValueBest, &ramaDBFileName, &seqData, &mediumData, &phipsiData, &hbscData))
        RETURN_OBJ_ERROR("need arguments: firstPPlaneFrag, lastPPlaneFrag, ppNbMin, minValueBest, maxValueBest, ramaDBFileName, seqData, mediumData, phipsiData, [hbscData]");

    if (lastPPlaneFrag - firstPPlaneFrag < 0)
        RETURN_OBJ_ERROR("lastPPlaneFrag < firstPPlaneFrag");

    params = NewParams();
    if (!params)
        RETURN_OBJ_ERROR("allocating Meccano params");

    result = runFwdParams(params, firstPPlaneFrag, lastPPlaneFrag, ppNbMin,
	minValueBest, maxValueBest, ramaDBFileName, seqData, mediumData, phipsiData, hbscData, error_msg);

    DelParams(params);

    return result;
}

/******************************************************************************
* METHOD REGISTRATION TABLE: NAME-STRING -> FUNCTION-POINTER
*
* List of functions defined in the module. A name->address method map, used
* to build-up the module's dictionary in "Py_InitModule". Once imported, this
* module acts just like it's coded in Python. The method functions handle
* converting data from/to python objects, and linkage to other C functions.
******************************************************************************/


static struct PyMethodDef Meccano_type_methods[] =
{
    { "runFwd",		(PyCFunction) runFwd,	METH_VARARGS },
    { NULL,		NULL,			0 }
};


/******************************************************************************
* INITIALIZATION FUNCTION (IMPORT-TIME)
*
* Initialization function for the module. Called on first "import Meccano" in 
* a Python program. The function is usually called "initMeccano": this name's
* added to the built-in module table in config.c statically (if added to file
* Module/Setup), or called when the module's loaded dynamically as a shareable 
* object-file found on PYTHONPATH. File and function names matter if dynamic.
******************************************************************************/

PY_MOD_INIT_FUNC initMeccano()
{
    PyObject *m, *d;

#ifdef WIN64
    Meccano_type.ob_type = &PyType_Type;
#endif
    /* create the module and add the functions */
    m = Py_InitModule("Meccano", Meccano_type_methods);

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("Meccano.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
    
    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module Meccano");
}


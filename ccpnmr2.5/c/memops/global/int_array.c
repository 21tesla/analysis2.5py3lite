
/*
======================COPYRIGHT/LICENSE START==========================

int_array.c: Part of the CcpNmr Analysis program

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
#include "int_array.h"

#include "hash_table.h"

Int_array new_int_array(int ndim, int *values)
{
    Int_array int_array = (Int_array) NULL;

    if (ndim > MAX_NDIM)
	return NULL;

    MALLOC_NEW(int_array, struct Int_array, 1);

    int_array->ndim = ndim;
    COPY_VECTOR(int_array->values, values, ndim);

    return int_array;
}

Int_array copy_int_array(Int_array int_array)
{
    Int_array int_array_new = (Int_array) NULL;

    if (int_array->ndim > MAX_NDIM)
	return NULL;

    MALLOC_NEW(int_array_new, struct Int_array, 1);

    int_array_new->ndim = int_array->ndim;
    COPY_VECTOR(int_array_new->values, int_array->values, int_array->ndim);

    return int_array_new;
}

void delete_int_array(Int_array int_array)
{
    FREE(int_array, struct Int_array);
}

int get_int_array(Int_array int_array, int dim)
{
    if ((dim >= 0) && (dim < int_array->ndim))
	return int_array->values[dim];

    return 0; /* arbitrary */
}

void set_int_array(Int_array int_array, int dim, int value)
{
    if ((dim >= 0) && (dim < int_array->ndim))
	int_array->values[dim] = value;
}

Bool equal_int_array(Int_array int_array1, Int_array int_array2)
{
    int i;
 
    if (int_array1->ndim != int_array2->ndim)
        return CCPN_FALSE;
 
    for (i = 0; i < int_array1->ndim; i++)
    {
        if (int_array1->values[i] != int_array2->values[i])
            return CCPN_FALSE;
    }
 
    return CCPN_TRUE;
}

Hash_value hash_int_array(Int_array int_array)
{
    return array_hash(int_array->ndim, int_array->values);
}


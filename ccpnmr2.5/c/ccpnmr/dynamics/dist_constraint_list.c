
/*
======================COPYRIGHT/LICENSE START==========================

dist_constraint_list.c: Part of the CcpNmr Analysis program

Copyright (C) 2005 Wayne Boucher and Tim Stevens (University of Cambridge)

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
#include "dist_constraint_list.h"

#define  NALLOC  100

Dist_constraint_list new_dist_constraint_list(int nconstraints,
	int *natom_pairs, int **atoms0, int **atoms1,
	float *dist_lower, float *dist_upper)
{
    int i;
    Dist_constraint dist_constraint;
    Dist_constraint_list dist_constraint_list;

    MALLOC_NEW(dist_constraint_list, struct Dist_constraint_list, 1);

    if (nconstraints > 0)
    {
	MALLOC_NEW(dist_constraint_list->dist_constraints, Dist_constraint, nconstraints);
	dist_constraint_list->ndist_constraints_alloc = nconstraints;
    }
    else
    {
	dist_constraint_list->dist_constraints = NULL;
	dist_constraint_list->ndist_constraints_alloc = 0;
    }

    dist_constraint_list->ndist_constraints = nconstraints;

    for (i = 0; i < dist_constraint_list->ndist_constraints_alloc; i++)
	dist_constraint_list->dist_constraints[i] = NULL;

    for (i = 0; i < nconstraints; i++)
    {
	dist_constraint = new_dist_constraint(natom_pairs[i], atoms0[i], atoms1[i],
					dist_lower[i], dist_upper[i]);
	if (!dist_constraint)
	    return NULL;

	dist_constraint_list->dist_constraints[i] = dist_constraint;
    }

    return dist_constraint_list;
}
 
void delete_dist_constraint_list(Dist_constraint_list dist_constraint_list)
{
    int i;

    if (dist_constraint_list)
    {
	for (i = 0; i < dist_constraint_list->ndist_constraints; i++)
	    delete_dist_constraint(dist_constraint_list->dist_constraints[i]);

	FREE(dist_constraint_list->dist_constraints, Dist_constraint);

	FREE(dist_constraint_list, struct Dist_constraint_list);
    }
}

void clear_dist_constraint_list(Dist_constraint_list dist_constraint_list)
{
    if (dist_constraint_list)
    {
	FREE(dist_constraint_list->dist_constraints, Dist_constraint);

	FREE(dist_constraint_list, struct Dist_constraint_list);
    }
}

CcpnStatus add_dist_constraint_list(Dist_constraint_list dist_constraint_list,
	int natom_pairs, int *atoms0, int *atoms1,
	float dist_lower, float dist_upper, CcpnString error_msg)
{
    int n;
    int nconstraints = dist_constraint_list->ndist_constraints;
    int nalloc = dist_constraint_list->ndist_constraints_alloc;
    Dist_constraint dist_constraint;

    if (nconstraints >= nalloc)
    {
	sprintf(error_msg, "allocating atom coords memory");
	if (nconstraints == 0)
	{
	    n = NALLOC;
	    MALLOC(dist_constraint_list->dist_constraints, Dist_constraint, n);
	}
	else
	{
	    n = nalloc + NALLOC;
	    REALLOC(dist_constraint_list->dist_constraints, Dist_constraint, n);
	}

	dist_constraint_list->ndist_constraints_alloc = n;
    }

    dist_constraint = new_dist_constraint(natom_pairs, atoms0, atoms1, dist_lower, dist_upper);

    if (!dist_constraint)
	RETURN_ERROR_MSG("allocating dist constraint memory");

    dist_constraint_list->dist_constraints[nconstraints] = dist_constraint;
    dist_constraint_list->ndist_constraints++;

    return CCPN_OK;
}

CcpnStatus append_dist_constraint_list(Dist_constraint_list dist_constraint_list,
			Dist_constraint dist_constraint, CcpnString error_msg)
{
    int n;
    int nconstraints = dist_constraint_list->ndist_constraints;
    int nalloc = dist_constraint_list->ndist_constraints_alloc;

    if (nconstraints >= nalloc)
    {
	sprintf(error_msg, "allocating atom coords memory");
	if (nconstraints == 0)
	{
	    n = NALLOC;
	    MALLOC(dist_constraint_list->dist_constraints, Dist_constraint, n);
	}
	else
	{
	    n = nalloc + NALLOC;
	    REALLOC(dist_constraint_list->dist_constraints, Dist_constraint, n);
	}

	dist_constraint_list->ndist_constraints_alloc = n;
    }

    dist_constraint_list->dist_constraints[nconstraints] = dist_constraint;
    dist_constraint_list->ndist_constraints++;

    return CCPN_OK;
}


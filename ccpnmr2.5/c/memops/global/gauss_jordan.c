/*
======================COPYRIGHT/LICENSE START==========================

gauss_jordan.c: Part of the CcpNmr Analysis program

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
#include "gauss_jordan.h"

/* Idea for code comes from Numerical Recipes, pp 24 - 29. */

/* gauss_jordan_vector solves (one) linear equation, a x = b. */

/* On input have a[n][n], b[n]. */
/* On output these replaced by a_inverse[n][n], x[n]. */

/* piv[n], row[n], col[n] (all ints) are used for storage. */

static CcpnStatus find_pivot(float **a, int *piv, int n, int *max_row, int *max_col)
{
    int i, j;
    float max;

    max = 0;

    for (i = 0; i < n; i++)
    {
	if (piv[i] != 1)
	{
	    for (j = 0; j < n; j++)
	    {
		if (piv[j] == 0)
		{
		    if (ABS(a[i][j]) >= max)
		    {
			max = ABS(a[i][j]);
			*max_row = i;
			*max_col = j;
		    }
		}
		else if (piv[j] > 1)
		{
		    return  CCPN_ERROR;
		}
	    }
	}
    }

    piv[*max_col]++;

    return  CCPN_OK;
}

static void interchange_rows_vector(float **a, float *b, int n,
						int max_row, int max_col)
{
    int j;

    for (j = 0; j < n; j++)
	SWAP(a[max_row][j], a[max_col][j], float);

    SWAP(b[max_row], b[max_col], float);
}

static CcpnStatus pivot_vector(float **a, float *b, int n, int max_col)
{
    int i, j;
    float x, piv_inv = a[max_col][max_col];

    if (piv_inv == 0)
	return  CCPN_ERROR;

    piv_inv = 1 / piv_inv;

    a[max_col][max_col] = 1;

    SCALE_VECTOR(a[max_col], a[max_col], piv_inv, n);
    b[max_col] *= piv_inv;

    for (i = 0; i < n; i++)
    {
	if (i != max_col)
	{
	    x = a[i][max_col];
	    a[i][max_col] = 0;

	    for (j = 0; j < n; j++)
		a[i][j] -= x * a[max_col][j];

	    b[i] -= x * b[max_col];
	}
    }

    return  CCPN_OK;
}

static void unscramble_vector(float **a, int n, int *row, int *col)
{
    int i, j;

    for (j = n-1; j >= 0; j--)
    {
	if (row[j] != col[j])
	{
	    for (i = 0; i < n; i++)
		SWAP(a[i][row[j]], a[i][col[j]], float);
	}
    }
}

void gauss_jordan_vector(float **a, float *b, int n, int *piv,
				int *row, int *col, Bool *singular)
{
    int i, max_row, max_col;

    ZERO_VECTOR(piv, n);

    for (i = 0; i < n; i++)
    {
	if (find_pivot(a, piv, n, &max_row, &max_col) == CCPN_ERROR)
	{
	    *singular = CCPN_TRUE;
	    return;
	}

	if (max_row != max_col)
	    interchange_rows_vector(a, b, n, max_row, max_col);

	row[i] = max_row;
	col[i] = max_col;

	if (pivot_vector(a, b, n, max_col) == CCPN_ERROR)
	{
	    *singular = CCPN_TRUE;
	    return;
	}

    }

    unscramble_vector(a, n, row, col);
    *singular = CCPN_FALSE;
}

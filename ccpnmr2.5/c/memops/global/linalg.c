
/*
======================COPYRIGHT/LICENSE START==========================

linalg.c: Part of the CcpNmr Analysis program

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

#include "linalg.h"

void matrix_matrix_multiply(float **matrix1, float **matrix2,
                        float **matrix3, int n1, int n2, int n3)
{
    int i, j, k;
    float t;

    for (i = 0; i < n1; i++)
    {
        for (j = 0; j < n2; j++)
        {
            t = 0;
            for (k = 0; k < n3; k++)
                t += matrix2[i][k] * matrix3[k][j];
            matrix1[i][j] = t;
        }
    }
}

void matrix_vector_multiply(float *vector1, float **matrix, float *vector2,
                                                                int n1, int n2)
{
    int i, j;
    float t;

    for (i = 0; i < n1; i++)
    {
        t = 0;
        for (j = 0; j < n2; j++)
            t += matrix[i][j] * vector2[j];
        vector1[i] = t;
    }
}

void matrix_transpose(float **matrix1, float **matrix2, int n1, int n2)
{
    int i, j;

    for (i = 0; i < n1; i++)
    {
        for (j = 0; j < n2; j++)
            matrix1[i][j] = matrix2[j][i];
    }
}

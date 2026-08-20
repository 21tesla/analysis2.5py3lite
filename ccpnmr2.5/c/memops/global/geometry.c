
/*
======================COPYRIGHT/LICENSE START==========================

geometry.c: Part of the CcpNmr Analysis program

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
#include "geometry.h"

/* return |v| */
float vector_length(float *v, int n)
{
    float d = inner_product(v, v, n);

    return (float) sqrt((double) d);
}

/* return v1 . v2 */
float inner_product(float *v1, float *v2, int n)
{
    int i;
    float d = 0;

    for (i = 0; i < n; i++)
        d += v1[i] * v2[i];

    return d;
}

/* set v1 = v2 ^ v3 */
void cross_product(float *v1, float *v2, float *v3)
{
    v1[0] = v2[1]*v3[2] - v2[2]*v3[1];
    v1[1] = v2[2]*v3[0] - v2[0]*v3[2];
    v1[2] = v2[0]*v3[1] - v2[1]*v3[0];
}

/* normalise v to have unit length */
void normalise_vector(float *v, int n)
{
    float d = vector_length(v, n);

    if (d > 0)
    {
        d = 1.0 / d;
        SCALE_VECTOR(v, v, d, n);
    }
}

/* return angle (radians) between v1 and v2 */
float vectors_angle(float *v1, float *v2, int n)
{
    float d1 = vector_length(v1, n);
    float d2 = vector_length(v2, n);
    float d, phi;

    if ((d1 > 0) && (d2 > 0))
    {
        d = inner_product(v1, v2, n) / (d1 * d2);
        d = MAX(-1.0, MIN(1.0, d));
        phi = (float) acos((double) d);
    }
    else
    {
        phi = 0; /* arbitrary */
    }

    return phi;
}

/* set matrix = rotation matrix for angle about axis */
void rotation_matrix(float **matrix, float *axis, float angle)
{
    int i, j, n = 3;
    float d;
    float c = (float) cos((double) angle);
    float s = (float) sin((double) angle);

    normalise_vector(axis, n);

    for (i = 0; i < n; i++)
    {
        ZERO_VECTOR(matrix[i], n);
        matrix[i][i] = c;

        for (j = 0; j < n; j++)
            matrix[i][j] += (1-c) * axis[i] * axis[j];
    }

    SCALE_VECTOR(axis, axis, s, n);

    matrix[0][1] -= axis[2];
    matrix[1][0] += axis[2];
    matrix[1][2] -= axis[0];
    matrix[2][1] += axis[0];
    matrix[2][0] -= axis[1];
    matrix[0][2] += axis[1];
}

#define  SMALL_ANGLE  1.0e-6

void rotation_matrix_vector_to_vector(float **matrix, float *v1, float *v2)
{
    int n = 3;
    float d1 = vector_length(v1, n);
    float d2 = vector_length(v2, n);
    float angle = vectors_angle(v1, v2, n);
    float axis[3];

    if (d1 > 0 && d2 > 0 && angle > 0)
    {
        if (angle < PI - SMALL_ANGLE)
        {
            cross_product(axis, v1, v2);
            normalise_vector(axis, n);
        }
        else /* anti parallel */
        {
            if (v1[0] != 0 || v1[1] != 0)
            {
                axis[0] = v1[1];
                axis[1] = -v1[0];
                axis[2] = 0;
            }
            else
            {
                axis[0] = 1;
                axis[1] = axis[2] = 0;
            }
        }
    }
    else /* parallel, do nothing */
    {
        angle = 0;
        axis[0] = 1;
        axis[1] = axis[2] = 0;
    }

    rotation_matrix(matrix, axis, angle);
}

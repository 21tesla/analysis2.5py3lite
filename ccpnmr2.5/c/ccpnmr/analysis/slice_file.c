/*
======================COPYRIGHT/LICENSE START==========================

slice_file.c: Part of the CcpNmr Analysis program

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
#include "slice_file.h"

Slice_file new_slice_file(int orient, int dim, Block_file block_file,
                          Mem_cache mem_cache)
{
    int ndim = block_file->ndim;
    Slice_file slice_file = (Slice_file) NULL;

    if ((dim < 0) || (dim >= ndim))
        return (Slice_file) NULL;

    MALLOC_NEW(slice_file, struct Slice_file, 1);

    slice_file->orient = orient;
    slice_file->dim = dim;
    slice_file->block_file = block_file;
    slice_file->mem_cache = mem_cache;

    return slice_file;
}

void delete_slice_file(Slice_file slice_file)
{
    FREE(slice_file, struct Slice_file);
}

CcpnStatus draw_slice_file(Slice_file slice_file,
                           int first, int last, float *position,
                           Drawing_funcs *drawing_funcs, Generic_ptr data,
                           CcpnString error_msg)
{
    Block_file block_file = slice_file->block_file;
    int ndim = block_file->ndim;
    int dim = slice_file->dim;
    int array[MAX_NDIM], cumul[MAX_NDIM], point[MAX_NDIM];
    int total, i, m, p, t;
    float v, value, weight, a0, b0, a1, b1;

    if (first < 0)
    {
        sprintf(error_msg, "first = %d < 0\n", first);
        return CCPN_ERROR;
    }

    if (last > block_file->points[dim])
    {
        sprintf(error_msg, "last = %d > points = %d\n",
                last, block_file->points[dim]);
        return CCPN_ERROR;
    }

    if (first >= last)
    {
        sprintf(error_msg, "first = %d >= last = %d\n", first, last);
        return CCPN_ERROR;
    }

    for (i = 0; i < ndim; i++)
    {
        if (i == dim)
            continue;

        if (position[i] < 0)
        {
            sprintf(error_msg, "dim %d: position = %4.3f < 0\n",
                    i+1, position[i]);
            return CCPN_ERROR;
        }

        if (position[i] >= block_file->points[i])
        {
            sprintf(error_msg, "dim %d: position = %4.3f >= points = %d\n",
                    i+1, position[i], block_file->points[i]);
            return CCPN_ERROR;
        }
    }

    total = 1;
    for (i = 0; i < ndim; i++)
    {
        if ((i == dim) || (position[i] >= block_file->points[i]-1))
            m = 1;
        else
            m = 2;

        cumul[i] = total;
        total *= m;
    }

    for (p = first; p < last; p++)
    {
        point[dim] = p;

        value = 0.0;
        for (t = 0; t < total; t++)
        {
            weight = 1.0;
            ARRAY_OF_INDEX(array, t, cumul, ndim);

            for (i = 0; i < ndim; i++)
            {
                if (i == dim)
                    continue;

                point[i] = floor(position[i]);
                if (array[i])
                {
                    point[i]++;
                    weight *= 1.0 - point[i] + position[i];
                }
                else
                {
                    weight *= 1.0 + point[i] - position[i];
                }
            }

            CHECK_STATUS(get_point_block_file(block_file, &v, point, error_msg));

            value += weight * v;
        }

        a1 = p;
        b1 = value;

        if (p > first)
        {
            if (slice_file->orient) /* horizontal */
                (drawing_funcs->draw_line)(data, a0, b0, a1, b1);
            else /* vertical */
                (drawing_funcs->draw_line)(data, b0, a0, b1, a1);
            /*
            printf("draw_slice_file: a0=%3.2f, b0=%3.2f, a1=%3.2f, b1=%3.2f\n", a0, b0, a1, b1);
            */
        }

        a0 = a1;
        b0 = b1;
    }

    return CCPN_OK;
}

CcpnStatus draw_all_slice_file(Slice_file slice_file,
                               int *first, int *last,
                               int ncomponents, int *components,
                               Drawing_funcs *drawing_funcs, Generic_ptr data,
                               CcpnString error_msg)
{
    Block_file block_file = slice_file->block_file;
    int ndim = block_file->ndim;
    int dim = slice_file->dim;
    int cumul[MAX_NDIM], point[MAX_NDIM];
    int total, i, m, p, t;
    float v, a0, b0, a1, b1;

    /*
    printf("draw_all_slice_file: first=%d, %d, %d, last=%d, %d, %d\n",
    first[0], first[1], first[2], last[0], last[1], last[2]);
    */
    for (i = 0; i < ndim; i++)
    {
        if (first[i] < 0)
        {
            sprintf(error_msg, "dim %d: first = %d < 0\n",
                    i+1, first[i]);
            return CCPN_ERROR;
        }

        if (last[i] > block_file->points[i])
        {
            sprintf(error_msg, "dim %d: last = %d > points = %d\n",
                    i+1, last[i], block_file->points[i]);
            return CCPN_ERROR;
        }

        if (first[i] >= last[i])
        {
            sprintf(error_msg, "dim %d: first = %d >= last = %d\n",
                    i+1, first[i], last[i]);
            return CCPN_ERROR;
        }
    }

    total = 1;
    for (i = 0; i < ndim; i++)
    {
        if (i == dim)
            m = 1;
        else
            m = last[i] - first[i];

        cumul[i] = total;
        total *= m;
    }

    for (t = 0; t < total; t++)
    {
        ARRAY_OF_INDEX(point, t, cumul, ndim);
        ADD_VECTORS(point, point, first, ndim);

        for (p = first[dim]; p < last[dim]; p++)
        {
            point[dim] = p;

            CHECK_STATUS(get_components_point_block_file(block_file, &v, point, ncomponents, components, error_msg));

            a1 = p;
            b1 = v;

            if (p > first[dim])
            {
                if (slice_file->orient) /* horizontal */
                    (drawing_funcs->draw_line)(data, a0, b0, a1, b1);
                else /* vertical */
                    (drawing_funcs->draw_line)(data, b0, a0, b1, a1);
            }
            /*
            printf("draw_all_slice_file: a0=%3.2f, b0=%3.2f, a1=%3.2f, b1=%3.2f\n", a0, b0, a1, b1);
            */
            a0 = a1;
            b0 = b1;
        }
    }

    return CCPN_OK;
}

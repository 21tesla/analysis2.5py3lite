/*
======================COPYRIGHT/LICENSE START==========================

contour_data.c: Part of the CcpNmr Analysis program

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
#include "contour_data.h"

#include "contourer.h"
#include "utility.h"

typedef struct
{
    Integer contour_xdelta;
    Integer contour_ydelta;
    Integer contour_index;
    Contour_data contour_data;
    void *contourer_info;
    Contour_vertices contour_vertices;
    float *contour_block_data;
    int *npoints;
    float *contour_data1;
    float *contour_data2;
}   Contourer_data;

static void copy_data0(Block_data block_data, Block_data old_block_data)
{
    int ndim = block_data->ndim, array[MAX_NDIM];
    Integer i, j;

    for (i = 0; i < old_block_data->size; i++)
    {
        ARRAY_OF_INDEX(array, i, old_block_data->cum_block_size, ndim);
        INDEX_OF_ARRAY(j, array, block_data->cum_block_size, ndim);
        block_data->data[j] = old_block_data->data[i];
        /*
        	printf("array0: %d %d %d\n", array[0], array[1], array[2]);
        	printf("copy_data0: i = %3d, j = %3d, data = %9.1f\n",
        					i, j,  old_block_data->data[i]);
        */
    }
}

static void copy_data1(int dim, Block_data block_data, Block_data old_block_data)
{
    int ndim = block_data->ndim, array[MAX_NDIM];
    Integer i, j, k, n, cum_block_size[MAX_NDIM];

    n = 1;
    for (i = 0; i < ndim; i++)
    {
        cum_block_size[i] = n;

        if (i != dim)
            n *= old_block_data->block_size[i];
    }

    for (i = 0; i < n; i++)
    {
        ARRAY_OF_INDEX(array, i, cum_block_size, ndim);
        INDEX_OF_ARRAY(k, array, old_block_data->cum_block_size, ndim);
        array[dim] = old_block_data->block_size[dim];
        INDEX_OF_ARRAY(j, array, block_data->cum_block_size, ndim);
        block_data->data[j] = old_block_data->data[k];
        /*
        	printf("array1: %d %d %d\n", array[0], array[1], array[2]);
        	printf("copy_data1: i = %3d, j = %3d, k = %3d, data = %9.1f\n",
        					i, j,  k, old_block_data->data[k]);
        */
    }
}

static void copy_data2(int dim1, int dim2, Block_data block_data,
                       Block_data old_block_data)
{
    int ndim = block_data->ndim, array[MAX_NDIM];
    Integer i, j, k, n, cum_block_size[MAX_NDIM];

    n = 1;
    for (i = 0; i < ndim; i++)
    {
        cum_block_size[i] = n;

        if ((i != dim1) && (i != dim2))
            n *= old_block_data->block_size[i];
    }

    for (i = 0; i < n; i++)
    {
        ARRAY_OF_INDEX(array, i, cum_block_size, ndim);
        INDEX_OF_ARRAY(k, array, old_block_data->cum_block_size, ndim);
        array[dim1] = old_block_data->block_size[dim1];
        array[dim2] = old_block_data->block_size[dim2];
        INDEX_OF_ARRAY(j, array, block_data->cum_block_size, ndim);
        block_data->data[j] = old_block_data->data[k];
        /*
        	printf("array2: %d %d %d\n", array[0], array[1], array[2]);
        	printf("copy_data2: i = %3d, j = %3d, k = %3d, data = %9.1f\n",
        					i, j,  k, old_block_data->data[k]);
        */
    }
}

static CcpnStatus alloc_block_data(Block_data block_data)
{
    MALLOC(block_data->data, Float, block_data->size);

    return CCPN_OK;
}

static Block_data get_extended_block(int xdim, int ydim,
                                     Block_file block_file, Int_array block,
                                     int ncomponents, int *components, CcpnString error_msg)
{
    int ndim, xb, yb;
    Block_data block_data = (Block_data) NULL, old_block_data;
    Bool xextended = CCPN_FALSE, yextended = CCPN_FALSE;
    Bool doUnlock = (block_file->block_file_kind == BLOCK_FILE);

    ndim = block_file->ndim;

    sprintf(error_msg, "allocating block data");
    MALLOC_NEW(block_data, struct Block_data, 1);
    block_data->ndim = ndim;
    block_data->data = (Float *) NULL;
    block_data->block = copy_int_array(block);

    if (!block_data->block)
    {
        delete_block_data(block_data);
        sprintf(error_msg, "allocating block data");
        return (Block_data) NULL;
    }

    old_block_data = get_components_block_data(block_file, block->values, ncomponents, components);
    if (!old_block_data)
    {
        sprintf(error_msg, "getting normal block data");
        delete_block_data(block_data);
        return (Block_data) NULL;
    }

    COPY_VECTOR(block_data->block_size, old_block_data->block_size, ndim);

    if (block_file->dim_wrapped[xdim]
            || (block_data->block->values[xdim] != (block_file->blocks[xdim]-1)))
    {
        xextended = CCPN_TRUE;
        block_data->block_size[xdim]++;
    }

    if (block_file->dim_wrapped[ydim]
            || (block_data->block->values[ydim] != (block_file->blocks[ydim]-1)))
    {
        yextended = CCPN_TRUE;
        block_data->block_size[ydim]++;
    }

    CUMULATIVE(block_data->cum_block_size, block_data->block_size, block_data->size, ndim);
    /*
        printf("block_data->cum_block_size = %d %d %d\n",
    			block_data->cum_block_size[0],
    			block_data->cum_block_size[1],
    			block_data->cum_block_size[2]);
    */

    if (alloc_block_data(block_data) == CCPN_ERROR)
    {
        sprintf(error_msg, "allocating block data");
        if (doUnlock)
            (void) unlock_block_data(block_file, old_block_data);
        delete_block_data(block_data);
        return (Block_data) NULL;
    }

    copy_data0(block_data, old_block_data);
    if (doUnlock && unlock_block_data(block_file, old_block_data) == CCPN_ERROR)
    {
        sprintf(error_msg, "unlocking block data");
        delete_block_data(block_data);
        return (Block_data) NULL;
    }

    if (xextended)
    {
        xb = block->values[xdim];
        block->values[xdim] = (xb+1) % block_file->blocks[xdim];
        old_block_data = get_components_block_data(block_file, block->values, ncomponents, components);
        block->values[xdim] = xb;

        if (!old_block_data)
        {
            sprintf(error_msg, "getting xextended block data");
            delete_block_data(block_data);
            return (Block_data) NULL;
        }

        copy_data1(xdim, block_data, old_block_data);
        if (doUnlock && unlock_block_data(block_file, old_block_data) == CCPN_ERROR)
        {
            sprintf(error_msg, "unlocking block data");
            delete_block_data(block_data);
            return (Block_data) NULL;
        }
    }

    if (yextended)
    {
        yb = block->values[ydim];
        block->values[ydim] = (yb+1) % block_file->blocks[ydim];
        old_block_data = get_components_block_data(block_file, block->values, ncomponents, components);
        block->values[ydim] = yb;

        if (!old_block_data)
        {
            sprintf(error_msg, "getting yextended block data");
            delete_block_data(block_data);
            return (Block_data) NULL;
        }

        copy_data1(ydim, block_data, old_block_data);
        if (doUnlock && unlock_block_data(block_file, old_block_data) == CCPN_ERROR)
        {
            sprintf(error_msg, "unlocking block data");
            delete_block_data(block_data);
            return (Block_data) NULL;
        }
    }

    if (xextended && yextended)
    {
        xb = block->values[xdim];
        block->values[xdim] = (xb+1) % block_file->blocks[xdim];
        yb = block->values[ydim];
        block->values[ydim] = (yb+1) % block_file->blocks[ydim];
        old_block_data = get_components_block_data(block_file, block->values, ncomponents, components);
        block->values[xdim] = xb;
        block->values[ydim] = yb;

        if (!old_block_data)
        {
            sprintf(error_msg, "getting xextended and yextended block data");
            delete_block_data(block_data);
            return (Block_data) NULL;
        }

        copy_data2(xdim, ydim, block_data, old_block_data);
        if (doUnlock && unlock_block_data(block_file, old_block_data) == CCPN_ERROR)
        {
            sprintf(error_msg, "unlocking block data");
            delete_block_data(block_data);
            return (Block_data) NULL;
        }
    }

    return block_data;
}

static CcpnStatus count_polys(void *user_data, int nvertices,
                              Contour_vertex first_vertex, CcpnString error_msg)
{
    int *npolys = (int *) user_data;

    (*npolys)++;

    return CCPN_OK;
}

static CcpnStatus set_poly(void *user_data, int nvertices,
                           Contour_vertex first_vertex, CcpnString error_msg)
{
    int i;
    Poly_line **polylines = (Poly_line **) user_data;
    Poly_line polyline = **polylines;
    Contour_vertex v = first_vertex;
    Bool closed;

    if (v->v1)
        closed = CCPN_TRUE;
    else
        closed = CCPN_FALSE;

    MALLOC(polyline->vertices, Point2f, nvertices);
    polyline->closed = closed;
    polyline->nvertices = nvertices;

    for (i = 0; i < nvertices; i++)
    {
        polyline->vertices[i].x = v->x[0];
        polyline->vertices[i].y = v->x[1];
        v = v->v2;
    }

    (*polylines)++;

    return CCPN_OK;
}

static CcpnStatus construct_polys(Contourer_data *contourer_data, int plane, int level,
                                  CcpnString error_msg)
{
    int i, npolys;
    Contour_data contour_data = contourer_data->contour_data;
    Contour_vertices contour_vertices = contourer_data->contour_vertices;
    Poly_line polyline, *polylines;

    npolys = 0;
    CHECK_STATUS(process_chains(contour_vertices, &npolys, count_polys, error_msg));

    if (npolys == 0)
    {
        /* not sure really need this */
        contour_data->npolylines[plane][level] = npolys;
        contour_data->polylines[plane][level] = (Poly_line *) NULL;
        return CCPN_OK;
    }

    sprintf(error_msg, "allocating memory for contour polys");
    MALLOC_ZERO(contour_data->polylines[plane][level], Poly_line, npolys);

    for (i = 0; i < npolys; i++)
    {
        MALLOC(contour_data->polylines[plane][level][i], struct Poly_line, 1);
        contour_data->polylines[plane][level][i]->vertices = (Point2f *) NULL;
    }

    contour_data->npolylines[plane][level] = npolys;

    polylines = contour_data->polylines[plane][level];
    sprintf(error_msg, "allocating polyline memory");
    CHECK_STATUS(process_chains(contour_vertices, &polylines, set_poly, error_msg));

    for (i = 0; i < npolys; i++)
        contour_data->nvertices += contour_data->polylines[plane][level][i]->nvertices;

    return CCPN_OK;
}

static CcpnStatus contour_plane_data(Contourer_data *contourer_data,
                                     int plane, CcpnString error_msg)
{
    int i;
    Contours contours;
    CcpnStatus status = CCPN_OK;

    if (!(contours = calculate_contours(contourer_data->contourer_info, error_msg)))
        return CCPN_ERROR;

    for (i = 0; i < contourer_data->contour_data->nlevels; i++)
    {
        contourer_data->contour_vertices = contours->vertices[i];
        status = construct_polys(contourer_data, plane, i, error_msg);
        if (status == CCPN_ERROR)
            break;
    }

    delete_contours(contours);

    return status;
}

static CcpnStatus get_data_row(void *user_data, float **p_data, CcpnString error_msg)
{
    Contourer_data *c = (Contourer_data *) user_data;
    int i, j = c->contour_index;
    float *data;

    *p_data = data = c->contour_data1;

    for (i = 0; i < c->npoints[0]; i++)
    {
        data[i] = c->contour_block_data[j];
        j += c->contour_xdelta;
    }

    c->contour_index += c->contour_ydelta;

    SWAP(c->contour_data1, c->contour_data2, float *);

    return CCPN_OK;
}

static CcpnStatus alloc_contourer_data(Contourer_data *contourer_data)
{
    int n;

    n = contourer_data->npoints[0];

    MALLOC(contourer_data->contour_data1, float, n);
    MALLOC(contourer_data->contour_data2, float, n);

    return CCPN_OK;
}

static void free_contourer_data(Contourer_data *contourer_data)
{
    if (!contourer_data)
        return;

    FREE(contourer_data->contour_data1, float);
    FREE(contourer_data->contour_data2, float);
}

static Contour_data alloc_contour_data(Int_array block,
                                       int xdim, int ydim,
                                       int *npoints, int *nblocks, int *block_size,
                                       int nlevels, CcpnString error_msg)
{
    int i, j, m, nplanes, ndim = block->ndim;
    Contour_data contour_data;

    sprintf(error_msg, "allocating contour memory");
    MALLOC_NEW(contour_data, struct Contour_data, 1);

    nplanes = 1;
    for (i = 0; i < ndim; i++)
    {
        contour_data->cum_planes[i] = nplanes;

        if ((i != xdim) && (i != ydim))
        {
            m = npoints[i] % block_size[i];

            if ((block->values[i] < (nblocks[i]-1)) || !m)
                m = block_size[i];

            nplanes *= m;
        }
    }

    contour_data->block = block;
    contour_data->nplanes = 0;
    contour_data->nlevels = 0;
    MALLOC_NEW(contour_data->npolylines, int *, nplanes);
    MALLOC_NEW(contour_data->polylines, Poly_line **, nplanes);
    contour_data->nplanes = nplanes;

    for (i = 0; i < nplanes; i++)
    {
        MALLOC_NEW(contour_data->npolylines[i], int, nlevels);
        MALLOC_NEW(contour_data->polylines[i], Poly_line *, nlevels);
        ZERO_VECTOR(contour_data->npolylines[i], nlevels);

        for (j = 0; j < nlevels; j++)
        {
            contour_data->npolylines[i][j] = 0;
            contour_data->polylines[i][j] = (Poly_line *) NULL;
        }
    }

    contour_data->nlevels = nlevels;

    return contour_data;
}

static Contour_data calculate_contour_data(Block_data block_data,
        int xdim, int ydim, Block_file block_file,
        Int_array block, Contour_levels contour_levels,
        CcpnString error_msg)
{
    int i, j, m, ndim = block_data->ndim;
    int npoints[2], array[MAX_NDIM];
    Contour_data contour_data;
    void *contourer_info;
    static float offset[2];
    static float scale[] = { 1, 1 };
    Contourer_data contourer_data;
    Integer ind;

    contour_data = alloc_contour_data(block, xdim, ydim, block_file->points,
                                      block_file->blocks, block_file->block_size,
                                      contour_levels->nlevels, error_msg);

    if (!contour_data)
        return NULL;

    m = block_file->points[xdim] % block_file->block_size[xdim];

    if ((block_data->block->values[xdim] < (block_file->blocks[xdim]-1)) || !m)
        m = block_data->block_size[xdim];

    npoints[0] = m;

    m = block_file->points[ydim] % block_file->block_size[ydim];

    if ((block_data->block->values[ydim] < (block_file->blocks[ydim]-1)) || !m)
        m = block_data->block_size[ydim];

    npoints[1] = m;

    offset[0] = block_data->block->values[xdim] * block_file->block_size[xdim];
    offset[1] = block_data->block->values[ydim] * block_file->block_size[ydim];

    if (!(contourer_info = new_contourer_info(&contourer_data,
                           contour_levels->nlevels, contour_levels->levels,
                           npoints, offset, scale, get_data_row, error_msg)))
    {
        delete_contour_data(contour_data);
        return (Contour_data) NULL;
    }

    contourer_data.contourer_info = contourer_info;
    contourer_data.contour_data = contour_data;
    contourer_data.npoints = npoints;
    contourer_data.contour_block_data = block_data->data;
    contourer_data.contour_xdelta = block_data->cum_block_size[xdim];
    contourer_data.contour_ydelta = block_data->cum_block_size[ydim];

    if (alloc_contourer_data(&contourer_data) == CCPN_ERROR)
    {
        sprintf(error_msg, "allocating contourer data");
        delete_contour_data(contour_data);
        delete_contourer_info(contourer_info);
        return (Contour_data) NULL;
    }

    contour_data->nvertices = 0;
    contour_data->npolys = 0;

    for (i = 0; i < contour_data->nplanes; i++)
    {
        ARRAY_OF_INDEX(array, i, contour_data->cum_planes, ndim);
        INDEX_OF_ARRAY(ind, array, block_data->cum_block_size, ndim);
        contourer_data.contour_index = ind;

        if (contour_plane_data(&contourer_data, i, error_msg) == CCPN_ERROR)
        {
            delete_contour_data(contour_data);
            delete_contourer_info(contourer_info);
            free_contourer_data(&contourer_data);
            return (Contour_data) NULL;
        }
    }

    delete_contourer_info(contourer_info);
    free_contourer_data(&contourer_data);

    return contour_data;
}

static CcpnStatus process_polys(void *user_data,
                                int *block, int *plane, int level, int npolylines,
                                Poly_line *polylines, CcpnString error_msg)
{
    Contour_data contour_data = (Contour_data) user_data;
    int ind, ndim = contour_data->block->ndim;

    INDEX_OF_ARRAY(ind, plane, contour_data->cum_planes, ndim);
    contour_data->npolylines[ind][level] = npolylines;
    contour_data->polylines[ind][level] = polylines;

    return CCPN_OK;
}

static Contour_data get_stored_contour_data(int xdim, int ydim,
        Store_file store_file, Int_array block, Bool transposed,
        CcpnString error_msg)
{
    Contour_data contour_data;
    int i, j, b, p, nlevels, ndim = block->ndim;
    int plane[MAX_NDIM];
    Bool do_plane;

    sprintf(error_msg, "inconsistent xdim");
    if (xdim != store_file->xdim)
        return NULL;

    sprintf(error_msg, "inconsistent ydim");
    if (ydim != store_file->ydim)
        return NULL;

    nlevels = 0;
    if (store_file->have_neg)
        nlevels++;

    if (store_file->have_pos)
        nlevels++;

    contour_data = alloc_contour_data(block, xdim, ydim, store_file->npoints,
                                      store_file->nblocks, store_file->block_size,
                                      nlevels, error_msg);

    if (!contour_data)
        return NULL;

    for (i = 0; i < contour_data->nplanes; i++)
    {
        do_plane = CCPN_TRUE;
        ARRAY_OF_INDEX(plane, i, contour_data->cum_planes, ndim);
        for (j = 0; j < ndim; j++)
        {
            if ((j != xdim) && (j != ydim))
            {
                b = block->values[j]*store_file->block_size[j];
                p = b + plane[j];
                if ((p < store_file->first[j]) || (p >= store_file->last[j]))
                {
                    do_plane = CCPN_FALSE;
                    break;
                }
            }
        }

        if (do_plane)
        {
            if (process_contours_store_file(store_file, block->values,
                                            plane, process_polys, contour_data, transposed,
                                            error_msg) == CCPN_ERROR)
            {
                delete_contour_data(contour_data);
                return NULL;
            }
        }
    }

    return contour_data;
}

Contour_data new_contour_data(int xdim, int ydim,
                              Block_file block_file, Store_file store_file, Int_array block,
                              Contour_levels contour_levels, int ncomponents, int *components,
                              Bool transposed, CcpnString error_msg)
{
    int ndim;
    Block_data block_data;
    Contour_data contour_data = NULL;

    if (block_file)
        ndim = block_file->ndim;
    else
        ndim = store_file->ndim;

    sprintf(error_msg, "illegal x dim = %d or y dim = %d", xdim, ydim);

    if ((xdim < 0) || (xdim >= ndim))
        return (Contour_data) NULL;

    if ((ydim < 0) || (ydim >= ndim))
        return (Contour_data) NULL;

    if (xdim == ydim)
        return (Contour_data) NULL;

    if (block_file)
    {
        block_data = get_extended_block(xdim, ydim, block_file, block, ncomponents, components, error_msg);

        if (!block_data)
            return (Contour_data) NULL;

        contour_data = calculate_contour_data(block_data, xdim, ydim,
                                              block_file, block, contour_levels, error_msg);

        delete_block_data(block_data);
    }
    else if (store_file)
    {
        contour_data = get_stored_contour_data(xdim, ydim, store_file, block,
                                               transposed, error_msg);
    }

    return contour_data;
}

void delete_contour_data(Contour_data contour_data)
{
    int i, j, k;

    for (i = 0; i < contour_data->nplanes; i++)
    {
        for (j = 0; j < contour_data->nlevels; j++)
        {
            for (k = 0; k < contour_data->npolylines[i][j]; k++)
            {
                FREE(contour_data->polylines[i][j][k]->vertices, Point2f);
                FREE(contour_data->polylines[i][j][k], struct Poly_line);
            }

            FREE(contour_data->polylines[i][j], Poly_line);
        }

        FREE(contour_data->polylines[i], Poly_line *);
        FREE(contour_data->npolylines[i], int);
    }

    FREE(contour_data->polylines, Poly_line **);
    FREE(contour_data->npolylines, int *);

    delete_int_array(contour_data->block);

    FREE(contour_data, struct Contour_data);
}

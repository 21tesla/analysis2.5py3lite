/*
======================COPYRIGHT/LICENSE START==========================

contour_file.c: Part of the CcpNmr Analysis program

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
#include "contour_file.h"

#include "int_array.h"
#include "store_handler.h"

static Bool equal_contour(Hash_key key1, Hash_key key2)
{
    Int_array block1 = (Int_array) key1;
    Int_array block2 = (Int_array) key2;

    return equal_int_array(block1, block2);
}

static Hash_value hash_contour(Hash_key key)
{
    Int_array block = (Int_array) key;

    return hash_int_array(block);
}

static CcpnStatus lock_contour_data(Contour_file contour_file, Contour_data contour_data)
{
    return lock_mem_cache(contour_file->mem_cache, (Generic_ptr) contour_data);
}

static CcpnStatus unlock_contour_data(Contour_file contour_file,
                                      Contour_data contour_data)
{
    return unlock_mem_cache(contour_file->mem_cache, (Generic_ptr) contour_data);
}

static void delete_data_contour(Generic_ptr object, Generic_ptr data)
{
    Contour_data contour_data = (Contour_data) object;
    Contour_file contour_file = (Contour_file) data;

    /*  if clearing then do not call remove_hash_table because
        this is already being done and would cause problems  */
    if (!contour_file->clearing)
        (void) remove_hash_table(contour_file->contour_table,
                                 (Hash_key) contour_data->block);
    delete_contour_data(contour_data);
}

static void delete_contour_cache(Hash_key key, Hash_data data,
                                 Generic_ptr user_data)
{
    Contour_data contour_data = (Contour_data) data;
    Contour_file contour_file = (Contour_file) user_data;

    /*  below automatically calls delete_data_contour  */
    (void) remove_mem_cache(contour_file->mem_cache,
                            (Generic_ptr) contour_data);
}

static Contour_data get_contour_data(Contour_file contour_file,
                                     int *block_array, CcpnString error_msg)
{
    Hash_data data;
    Contour_data contour_data;
    Int_array block;
    int ndim;

    if (contour_file->block_file)
        ndim = contour_file->block_file->ndim;
    else
        ndim = contour_file->store_file->ndim;

    block = new_int_array(ndim, block_array);

    if (!block)
    {
        sprintf(error_msg, "allocating contour data");
        return (Contour_data) NULL;
    }

    if (find_hash_table(contour_file->contour_table, (Hash_key) block, &data))
    {
        contour_data = (Contour_data) data;
        if (lock_contour_data(contour_file, contour_data) == CCPN_ERROR)
        {
            sprintf(error_msg, "locking contour data");
            return (Contour_data) NULL;
        }
    }

    /* do it again in case deleted before lock obtained */
    if (data &&
            find_hash_table(contour_file->contour_table, (Hash_key) block, &data))
    {
        delete_int_array(block);
    }
    else
    {
        contour_data = new_contour_data(contour_file->xdim, contour_file->ydim,
                                        contour_file->block_file, contour_file->store_file, block,
                                        contour_file->contour_levels, contour_file->ncomponents,
                                        contour_file->components, contour_file->transposed, error_msg);
        if (!contour_data)
            return (Contour_data) NULL;

        if (insert_hash_table(contour_file->contour_table,
                              (Hash_key) contour_data->block, (Hash_data) contour_data) == CCPN_ERROR)
        {
            delete_contour_data(contour_data);
            return (Contour_data) NULL;
        }

        if (add_mem_cache(contour_file->mem_cache, (Generic_ptr) contour_data,
                          sizeof(Point2f)*contour_data->nvertices,
                          delete_data_contour, contour_file) == CCPN_ERROR)
        {
            delete_data_contour((Generic_ptr) contour_data,
                                (Generic_ptr) contour_file);
            return (Contour_data) NULL;
        }
    }

    return contour_data;
}

static delete_components(Contour_file contour_file)
{
    FREE(contour_file->components, int);
    contour_file->ncomponents = 0;
}

static Bool equal_components(Contour_file contour_file, int ncomponents, int *components)
{
    int i;

    if (!contour_file->components && !components)
        return CCPN_TRUE;

    if (contour_file->components && !components)
    {
        if (contour_file->ncomponents == 0)
            return CCPN_TRUE;
        else
            return CCPN_FALSE;
    }

    if (!contour_file->components && components)
    {
        if (ncomponents == 0)
            return CCPN_TRUE;
        else
            return CCPN_FALSE;
    }

    if (contour_file->ncomponents != ncomponents)
        return CCPN_FALSE;

    for (i = 0; i < ncomponents; i++)
    {
        if (contour_file->components[i] != components[i])
            return CCPN_FALSE;
    }

    return CCPN_TRUE;
}

static CcpnStatus copy_components(Contour_file contour_file, int ncomponents,
                                  int *components, CcpnString error_msg)
{
    delete_components(contour_file);

    if (!components || (ncomponents <= 0))
        return CCPN_OK;

    sprintf(error_msg, "allocating components");
    MALLOC(contour_file->components, int, ncomponents);
    COPY_VECTOR(contour_file->components, components, ncomponents);
    contour_file->ncomponents = ncomponents;

    return CCPN_OK;
}

Contour_file new_contour_file(int xdim, int ydim,
                              Block_file block_file, Store_file store_file,
                              Mem_cache mem_cache, Bool transposed, CcpnString error_msg)
{
    int ndim;
    Contour_file contour_file = (Contour_file) NULL;

    if (block_file && store_file)
    {
        sprintf(error_msg, "both block_file and store_file are set");
        return (Contour_file) NULL;
    }

    if (!block_file && !store_file)
    {
        sprintf(error_msg, "neither block_file nor store_file is set");
        return (Contour_file) NULL;
    }

    if (block_file)
        ndim = block_file->ndim;
    else
        ndim = store_file->ndim;

    if ((xdim < 0) || (xdim >= ndim))
    {
        sprintf(error_msg, "xdim = %d, should be >= 0 and < %d", xdim, ndim);
        return (Contour_file) NULL;
    }

    if ((ydim < 0) || (ydim >= ndim))
    {
        sprintf(error_msg, "ydim = %d, should be >= 0 and < %d", ydim, ndim);
        return (Contour_file) NULL;
    }

    if (xdim == ydim)
    {
        sprintf(error_msg, "xdim = ydim = %d", xdim);
        return (Contour_file) NULL;
    }

    sprintf(error_msg, "allocating Contour_file");
    MALLOC_NEW(contour_file, struct Contour_file, 1);

    contour_file->xdim = xdim;
    contour_file->ydim = ydim;
    contour_file->block_file = block_file;
    contour_file->store_file = store_file;
    contour_file->mem_cache = mem_cache;
    contour_file->transposed = transposed;
    contour_file->contour_levels = (Contour_levels) NULL;
    contour_file->ncomponents = 0;
    contour_file->components = (int *) NULL;
    contour_file->clearing = CCPN_FALSE;

    contour_file->contour_table = new_hash_table(equal_contour, hash_contour);

    if (!contour_file->contour_table)
    {
        delete_contour_file(contour_file);
        return NULL;
    }

    return contour_file;
}

void delete_contour_file(Contour_file contour_file)
{
    contour_file->clearing = CCPN_TRUE;
    clear_hash_table(contour_file->contour_table,
                     delete_contour_cache, (Generic_ptr) contour_file);
    delete_hash_table(contour_file->contour_table);
    delete_contour_levels(contour_file->contour_levels);
    delete_components(contour_file);
    FREE(contour_file, struct Contour_file);
}

static void draw_plane(int plane, Contour_data contour_data,
                       Contour_levels contour_levels, Contour_style contour_style,
                       Drawing_funcs *drawing_funcs, Generic_ptr data, Bool *do_level)
{
    int i, j, m, n;
    int npos_colors = contour_style->npos_colors, pos_ind = 0;
    int nneg_colors = contour_style->nneg_colors, neg_ind = 0;

    for (i = 0; i < contour_levels->nlevels; i++)
    {
        if (do_level && !do_level[i])
            continue;

        if (drawing_funcs->draw_medium == STORE_DRAWING)
            init_store_level((Store_handler) data, contour_levels->levels[i]);

        n = contour_data->npolylines[plane][i];

        if (n > 0)
        {
            if (contour_levels->levels[i] > 0)
            {
                m = (pos_ind++) % npos_colors;
                (drawing_funcs->set_draw_color)(data, contour_style->pos_colors[m]);
                (drawing_funcs->set_line_style)(data, contour_style->pos_line_style);
            }
            else
            {
                m = (neg_ind++) % nneg_colors;
                (drawing_funcs->set_draw_color)(data, contour_style->neg_colors[m]);
                (drawing_funcs->set_line_style)(data, contour_style->neg_line_style);
            }
        }

        for (j = 0; j < n; j++)
        {
            (drawing_funcs->draw_clipped_polyline)(data,
                                                   contour_data->polylines[plane][i][j]);
        }
    }
}

CcpnStatus region_contour_file(Contour_file contour_file,
                               int *first, int *last, Abort_func abort_func,
                               Contour_levels contour_levels, Contour_style contour_style,
                               int ncomponents, int *components, Drawing_funcs *drawing_funcs,
                               Generic_ptr data, CcpnString error_msg)
{
    Block_file block_file = contour_file->block_file;
    Store_file store_file = contour_file->store_file;
    int i, s;
    int ndim, *block_size;
    int xdim = contour_file->xdim, ydim = contour_file->ydim;
    int cum_blocks[MAX_NDIM], block[MAX_NDIM];
    int block_min[MAX_NDIM], nblocks[MAX_NDIM], block_max[MAX_NDIM];
    int plane_min[MAX_NDIM], cum_nplanes[MAX_NDIM], plane_max;
    int point_min, point_max;
    int array[MAX_NDIM];
    int b, total_nblocks, p, q, n, total_nplanes;
    Contour_data contour_data;
    static float stored_neg_levels[] = {-1};
    static float stored_pos_levels[] = {1};
    static float stored_neg_pos_levels[] = {-1,1};
    int nlevels;
    float *levels;
    Bool stored_do_level[2], *do_level;

    if (block_file)
    {
        ndim = block_file->ndim;
        block_size = block_file->block_size;
    }
    else if (store_file)
    {
        ndim = store_file->ndim;
        block_size = store_file->block_size;
    }

    for (i = 0; i < ndim; i++)
    {
        if (first[i] < 0)
        {
            sprintf(error_msg, "dim %d: region min = %d < 0\n",
                    i+1, first[i]);
            return CCPN_ERROR;
        }

        if (block_file)
            n = block_file->points[i];
        else
            n = store_file->npoints[i];

        if (last[i] > n)
        {
            sprintf(error_msg, "dim %d: region max = %d > points = %d\n",
                    i+1, last[i], n);
            return CCPN_ERROR;
        }

        if (first[i] >= last[i])
        {
            sprintf(error_msg, "dim %d: region min = %d >= max = %d\n",
                    i+1, first[i], last[i]);
            return CCPN_ERROR;
        }
    }

    if (store_file)
    {
        for (i = 0; i < ndim; i++)
        {
            if (first[i] >= store_file->last[i])
                return CCPN_OK; /* nothing to draw */

            if (last[i] <= store_file->first[i])
                return CCPN_OK; /* nothing to draw */

            first[i] = MAX(first[i], store_file->first[i]);
            last[i] = MIN(last[i], store_file->last[i]);
        }
    }

    if (block_file)
    {
        if (!contour_file->contour_levels ||
                !equal_contour_levels(contour_file->contour_levels, contour_levels))
        {
            contour_file->clearing = CCPN_TRUE;
            clear_hash_table(contour_file->contour_table,
                             delete_contour_cache, (Generic_ptr) contour_file);
            contour_file->clearing = CCPN_FALSE;

            delete_contour_levels(contour_file->contour_levels);
            /* note that take copy instead of just using pointer */
            contour_file->contour_levels = copy_contour_levels(contour_levels);

            if (!contour_file->contour_levels)
            {
                sprintf(error_msg, "allocating contour levels");
                return CCPN_ERROR;
            }
        }

        if (!equal_components(contour_file, ncomponents, components))
        {
            contour_file->clearing = CCPN_TRUE;
            clear_hash_table(contour_file->contour_table,
                             delete_contour_cache, (Generic_ptr) contour_file);
            contour_file->clearing = CCPN_FALSE;

            /* note that take copy instead of just using pointer */
            CHECK_STATUS(copy_components(contour_file, ncomponents, components, error_msg));
        }

        do_level = NULL;
    }
    else /* store_file */
    {
        if (!contour_file->contour_levels)
        {
            if (store_file->have_neg && store_file->have_pos)
            {
                nlevels = 2;
                levels = stored_neg_pos_levels;
            }
            else if (store_file->have_neg)
            {
                nlevels = 1;
                levels = stored_neg_levels;
            }
            else /* store_file->have_pos */
            {
                nlevels = 1;
                levels = stored_pos_levels;
            }

            contour_file->contour_levels = new_contour_levels(nlevels, levels);

            if (!contour_file->contour_levels)
            {
                sprintf(error_msg, "allocating contour levels");
                return CCPN_ERROR;
            }
        }

        do_level = stored_do_level;
        if (store_file->have_neg && store_file->have_pos)
        {
            do_level[0] = have_neg_contour_levels(contour_levels);
            do_level[1] = have_pos_contour_levels(contour_levels);
        }
        else if (store_file->have_neg)
        {
            if (have_neg_contour_levels(contour_levels))
                do_level[0] = CCPN_TRUE;
            else
                return CCPN_OK;
        }
        else /* store_file->have_pos */
        {
            if (have_pos_contour_levels(contour_levels))
                do_level[0] = CCPN_TRUE;
            else
                return CCPN_OK;
        }
    }

    total_nblocks = 1;
    for (i = 0; i < ndim; i++)
    {
        block_min[i] = first[i] / block_size[i];
        block_max[i] = (last[i]-1) / block_size[i];
        nblocks[i] =  block_max[i] - block_min[i] + 1;
        cum_blocks[i] = total_nblocks;
        total_nblocks *= nblocks[i];
    }

    if (drawing_funcs->draw_medium == STORE_DRAWING)
    {
        CHECK_STATUS(init_store_save((Store_handler) data, ndim, xdim, ydim, block_file->points, first, last, block_size, nblocks, contour_levels->nlevels, contour_levels->levels, error_msg));
        for (i = 0; i < ndim; i++)
            ((Store_handler) data)->first_plane[i] = first[i] % block_size[i];
    }

    for (b = 0; b < total_nblocks; b++)
    {
        ARRAY_OF_INDEX(block, b, cum_blocks, ndim);

        if (drawing_funcs->draw_medium == STORE_DRAWING)
            init_store_block((Store_handler) data, block);

        ADD_VECTORS(block, block, block_min, ndim);

        contour_data = get_contour_data(contour_file, block, error_msg);

        if (!contour_data)
            return CCPN_ERROR;

        total_nplanes = 1;
        for (i = 0; i < ndim; i++)
        {
            cum_nplanes[i] = total_nplanes;

            if ((i != xdim) && (i != ydim))
            {
                point_min = block[i] * block_size[i];

                if (block[i] == block_min[i])
                    point_min += first[i] % block_size[i];

                point_max = (block[i] + 1) * block_size[i];

                if (block[i] == block_max[i])
                {
                    s =  last[i] % block_size[i];
                    s = block_size[i] - s;
                    point_max -= block_size[i] - last[i] % block_size[i];
                }

                plane_min[i] = point_min % block_size[i];

                plane_max = point_max % block_size[i];
                if (plane_max == 0)
                    plane_max = block_size[i];

                n = plane_max - plane_min[i];
                total_nplanes *= n;
            }
            else
            {
                plane_min[i] = 0;
            }
        }

        for (p = 0; p < total_nplanes; p++)
        {
            if (abort_func())
                break;

            ARRAY_OF_INDEX(array, p, cum_nplanes, ndim);

            if (drawing_funcs->draw_medium == STORE_DRAWING)
                init_store_plane((Store_handler) data, array);

            ADD_VECTORS(array, array, plane_min, ndim);
            INDEX_OF_ARRAY(q, array, contour_data->cum_planes, ndim);

            draw_plane(q, contour_data, contour_file->contour_levels,
                       contour_style, drawing_funcs, data, do_level);

            if (drawing_funcs->draw_medium == STORE_DRAWING)
                end_store_plane((Store_handler) data);
        }

        if (unlock_contour_data(contour_file, contour_data) == CCPN_ERROR)
        {
            sprintf(error_msg, "unlocking contour data");
            return CCPN_ERROR;
        }
    }

    if (drawing_funcs->draw_medium == STORE_DRAWING)
        end_store_save((Store_handler) data);

    return CCPN_OK;
}

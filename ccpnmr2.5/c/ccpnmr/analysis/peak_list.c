/*
======================COPYRIGHT/LICENSE START==========================

peak_list.c: Part of the CcpNmr Analysis program

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
#include "peak_list.h"

#include "symbol.h"

#define  NCOLORS  3

#define  NALLOC  50

Peak_list new_peak_list(int ndim, int *npoints)
{
    Peak_list peak_list;

    MALLOC_NEW(peak_list, struct Peak_list, 1);
    MALLOC_NEW(peak_list->npoints, int, ndim);
    MALLOC_NEW(peak_list->color, float, NCOLORS);

    COPY_VECTOR(peak_list->npoints, npoints, ndim);
    ZERO_VECTOR(peak_list->color, NCOLORS); /* make black default */
    peak_list->ndim = ndim;
    peak_list->symbol = DEFAULT_SYMBOL;
    peak_list->npeaks = 0;
    peak_list->npeaks_alloc = 0;
    peak_list->peaks = NULL;

    return peak_list;
}

void delete_peak_list(Peak_list peak_list)
{
    int i;

    /*
        printf("entering delete_peak_list\n");
    */

    if (peak_list)
    {
        for (i = 0; i < peak_list->npeaks; i++)
            delete_peak(peak_list->peaks[i]);

        FREE(peak_list->npoints, int);
        FREE(peak_list->color, float);
        FREE(peak_list->peaks, Peak);
    }

    FREE(peak_list, struct Peak_list);
}

void is_selected_peak_list(Peak_list peak_list, Bool isSelected)
{
    int i;

    for (i = 0; i < peak_list->npeaks; i++)
        set_is_selected_peak(peak_list->peaks[i], isSelected);
}

void set_color_peak_list(Peak_list peak_list, float *color)
{
    COPY_VECTOR(peak_list->color, color, NCOLORS);
}

CcpnStatus set_symbol_peak_list(Peak_list peak_list, int symbol,
                                CcpnString error_msg)
{
    /* TBD: check if symbol in range */
    peak_list->symbol = symbol;

    return CCPN_OK;
}

CcpnStatus add_peak_peak_list(Peak_list peak_list, Peak *p_peak, CcpnString error_msg)
{
    int n;
    Peak peak;

    if (peak_list->npeaks == peak_list->npeaks_alloc)
    {
        sprintf(error_msg, "allocating peak memory");
        if (peak_list->npeaks_alloc == 0)
        {
            n = NALLOC;
            MALLOC(peak_list->peaks, Peak, n);
        }
        else
        {
            n = peak_list->npeaks_alloc + NALLOC;
            REALLOC(peak_list->peaks, Peak, n);
        }

        peak_list->npeaks_alloc = n;
    }

    peak = new_peak(peak_list->ndim);

    if (!peak)
        RETURN_ERROR_MSG("allocating peak memory");

    peak_list->peaks[peak_list->npeaks] = peak;
    peak_list->npeaks++;

    *p_peak = peak;

    return CCPN_OK;
}

CcpnStatus remove_peak_peak_list(Peak_list peak_list, Peak peak,
                                 CcpnString error_msg)
{
    int i;

    for (i = 0; i < peak_list->npeaks; i++)
    {
        if (peak == peak_list->peaks[i])
            break;
    }

    if (i == peak_list->npeaks)
        RETURN_ERROR_MSG("unknown peak");

    delete_peak(peak);

    for (; i < peak_list->npeaks-1; i++)
        peak_list->peaks[i] = peak_list->peaks[i+1];

    peak_list->npeaks--;

    return CCPN_OK;
}

void remove_selected_peak_list(Peak_list peak_list)
{
    int i, j;
    Peak peak;

    for (i = j = 0; i < peak_list->npeaks; i++)
    {
        peak = peak_list->peaks[i];
        if (peak->isSelected)
            delete_peak(peak_list->peaks[i]);
        else
            peak_list->peaks[j++] = peak;
    }

    peak_list->npeaks = j;
}

void unselect_selected_peak_list(Peak_list peak_list)
{
    int i;

    for (i = 0; i < peak_list->npeaks; i++)
        set_is_selected_peak(peak_list->peaks[i], CCPN_FALSE);
}

static void find_point(int ndim, int n, int *point, int *cum_points, int *offset)
{
    ARRAY_OF_INDEX(point, n, cum_points, ndim);
    ADD_VECTORS(point, point, offset, ndim);
}

static Bool peak_within_buffer(int ndim, int *point, Peak peak, int *buffer)
{
    int i, d;

    if (!peak->position)
        return CCPN_FALSE;

    for (i = 0; i < ndim; i++)
    {
        d = ABS(point[i] - peak->position[i] + 1);
        /* +1 because positions start at 1 */
        if (d > buffer[i])
            return CCPN_FALSE;
    }

    return CCPN_TRUE;
}

static CcpnStatus check_new_peak(Peak_list peak_list, int *buffer,
                                 int *point, int ignore_peak, CcpnString error_msg)
{
    int i, ndim = peak_list->ndim;
    float position[MAX_NDIM];
    Peak peak;

    /* check that not within buffer */
    /* could do this first if implement efficient way of doing so */

    for (i = 0; i < peak_list->npeaks; i++)
    {
        if (ignore_peak == i)
            continue;

        if (peak_within_buffer(ndim, point, peak_list->peaks[i], buffer))
            return CCPN_OK;
    }

    /* have new peak, add to peak_list */

    CHECK_STATUS(add_peak_peak_list(peak_list, &peak, error_msg));

    for (i = 0; i < ndim; i++)
        position[i] = point[i] + 1;

    set_position_peak(peak, position);

    return CCPN_OK;
}

static CcpnStatus check_nonadjacent_points(Block_file block_file,
        Bool find_maximum, int *buffer, float v, int *point,
        int npoints, int *cumulative, Bool *ok_extreme,
        Bool *dim_checked, CcpnString error_msg)
{
    int i, n, zero_index, p[MAX_NDIM], ndim = block_file->ndim;
    Bool do_point;
    float v2;

    *ok_extreme = CCPN_FALSE;

    zero_index = 1;
    for (i = 0; i < ndim; i++)
        zero_index *= 3;
    zero_index = (zero_index - 1) / 2;

    /* check that local extremum */

    for (n = 0; n < npoints; n++)
    {
        if (n == zero_index) /* this is the central point */
            continue;

        ARRAY_OF_INDEX(p, n, cumulative, ndim);

        do_point = CCPN_TRUE;
        for (i = 0; i < ndim; i++)
        {
            /* would have been more efficient to do this before this function called... */
            if (!dim_checked[i] && (p[i] != 0))
            {
                do_point = CCPN_FALSE;
                break;
            }

            p[i] += point[i] - 1;
            /* p initially goes 0, 1, 2 so extra -1 makes it -1, 0, 1 */

            if (block_file->dim_wrapped[i])
            {
                if (p[i] < 0)
                    p[i] = block_file->points[i] - 1;
                else if (p[i] >= block_file->points[i])
                    p[i] = 0;
            }
            else if ((p[i] < 0) || (p[i] >= block_file->points[i]))
            {
                do_point = CCPN_FALSE;
                break;
            }
        }

        if (!do_point)
            continue;

        CHECK_STATUS(get_point_block_file(block_file, &v2, p, error_msg));

        if (find_maximum)
        {
            if (v2 > v)
                return CCPN_OK;
        }
        else
        {
            if (v2 < v)
                return CCPN_OK;
        }
    }

    *ok_extreme = CCPN_TRUE;

    return CCPN_OK;
}

/* TBD: ignores aliasing so does not work correctly on boundaries */
static CcpnStatus check_adjacent_points(Block_file block_file, Bool find_maximum,
                                        int *buffer, float v, int *point, Bool *ok_extreme,
                                        Bool *dim_checked, CcpnString error_msg)
{
    int i, p, ndim = block_file->ndim;
    float v2;
    CcpnStatus status;

    *ok_extreme = CCPN_FALSE;

    /* check that local extremum */

    for (i = 0; i < ndim; i++)
    {
        if (!dim_checked[i])
            continue;

        if (block_file->dim_wrapped[i] || point[i] > 0)
        {
            p = point[i];
            if (p == 0)
                point[i] = block_file->points[i] - 1;
            else
                point[i] = p - 1;
            status = get_point_block_file(block_file, &v2, point, error_msg);
            point[i] = p;

            CHECK_STATUS(status);

            if (find_maximum)
            {
                if (v2 > v)
                    return CCPN_OK;
            }
            else
            {
                if (v2 < v)
                    return CCPN_OK;
            }
        }

        if (block_file->dim_wrapped[i] || (point[i] < (block_file->points[i] - 1)))
        {
            p = point[i];
            point[i] = (p + 1) % block_file->points[i];
            status = get_point_block_file(block_file, &v2, point, error_msg);
            point[i] = p;

            CHECK_STATUS(status);

            if (find_maximum)
            {
                if (v2 > v)
                    return CCPN_OK;
            }
            else
            {
                if (v2 < v)
                    return CCPN_OK;
            }
        }
    }

    *ok_extreme = CCPN_TRUE;

    return CCPN_OK;
}

static CcpnStatus drops_in_direction(Block_file block_file, Bool find_maximum,
                                     float drop, float v, int *point, int dim, int dirn,
                                     Bool *ok_drop, CcpnString error_msg)
{
    int i, i_start, i_end, i_step;
    float v_prev = v, v_this;
    int q[MAX_NDIM];

    /*
    printf("drops_in_direction0: dim = %d, dirn = %d, v = %3.2f, point = %d %d %d\n",
    dim, dirn, v, point[0], point[1], point[2]);
    */
    if (dirn == 1)
    {
        i_start = point[dim] + 1;
        if (block_file->dim_wrapped[dim])
            i_end = i_start + block_file->points[dim] - 1;
        else
            i_end = block_file->points[dim];
        i_step = 1;
    }
    else
    {
        i_start = point[dim] - 1;
        if (block_file->dim_wrapped[dim])
            i_end = i_start - block_file->points[dim] + 1;
        else
            i_end = -1;
        i_step = -1;
    }

    /*
    printf("drops_in_direction1: i_start = %d, i_end = %d, i_step = %d\n",
    i_start, i_end, i_step);
    */
    COPY_VECTOR(q, point, block_file->ndim);

    for (i = i_start; i != i_end; i += i_step)
    {
        q[dim] = (i + block_file->points[dim]) % block_file->points[dim];
        /*
        printf("drops_in_direction2: i = %d, q = %d %d %d\n", i, q[0], q[1], q[2]);
        */
        CHECK_STATUS(get_point_block_file(block_file, &v_this, q, error_msg));
        /*
        printf("drops_in_direction3: i = %d, v_this = %3.2f, v_prev = %3.2f\n",
        i, v_this, v_prev);
        */

        if (find_maximum)
        {
            if (v_this > v_prev)
            {
                *ok_drop = CCPN_FALSE;
                break;
            }
            else if ((v-v_this) >= drop)
            {
                break;
            }
        }
        else
        {
            if (v_this < v_prev)
            {
                *ok_drop = CCPN_FALSE;
                break;
            }
            else if ((v_this-v) >= drop)
            {
                break;
            }
        }

        v_prev = v_this;
    }

    return CCPN_OK;
}

static CcpnStatus check_drop(Block_file block_file, Bool find_maximum,
                             float drop_factor, float v, int *point, Bool *ok_drop,
                             CcpnString error_msg)
{
    int i, ndim = block_file->ndim;
    float drop = drop_factor * ABS(v);

    *ok_drop = CCPN_TRUE;

    if (drop_factor <= 0)
        return CCPN_OK;

    /*
    printf("check_drop1: v = %3.2f, drop = %3.2f\n", v, drop);
    */
    for (i = 0; i < ndim; i++)
    {
        /*
        printf("check_drop2: i = %d\n", i);
        */
        if (drops_in_direction(block_file, find_maximum, drop, v, point,
                               i, 1, ok_drop, error_msg) == CCPN_ERROR)
            return CCPN_ERROR;

        /*
        printf("check_drop3: i = %d\n", i);
        */
        if (*ok_drop == CCPN_FALSE)
            return CCPN_OK;

        if (drops_in_direction(block_file, find_maximum, drop, v, point,
                               i, -1, ok_drop, error_msg) == CCPN_ERROR)
            return CCPN_ERROR;

        /*
        printf("check_drop4: i = %d\n", i);
        */
        if (*ok_drop == CCPN_FALSE)
            return CCPN_OK;
    }

    return CCPN_OK;
}

static CcpnStatus half_max_position(Block_file block_file, Bool find_maximum,
                                    float v, int *point, int dim, int dirn,
                                    float *half_max, CcpnString error_msg)
{
    int i, i_start, i_end, i_step;
    float v_half = 0.5 * v, v_prev = v, v_this;
    int q[MAX_NDIM];

    /*
    printf("half_max_position1: v = %3.2f, dim = %d, dirn = %d, point = %d %d %d\n",
    v, dim, dirn, point[0], point[1], point[2]);
    */
    if (dirn == 1)
    {
        i_start = point[dim] + 1;
        if (block_file->dim_wrapped[dim])
            i_end = i_start + block_file->points[dim] - 1;
        else
            i_end = block_file->points[dim];
        i_step = 1;
    }
    else
    {
        i_start = point[dim] - 1;
        if (block_file->dim_wrapped[dim])
            i_end = i_start - block_file->points[dim] + 1;
        else
            i_end = -1;
        i_step = -1;
    }
    /*
    printf("half_max_position2: i_start = %d, i_end = %d, i_step = %d\n",
    i_start, i_end, i_step);
    */

    COPY_VECTOR(q, point, block_file->ndim);

    for (i = i_start; i != i_end; i += i_step)
    {
        q[dim] = (i + block_file->points[dim]) % block_file->points[dim];
        /*
        printf("half_max_position3: i = %d, q = %d %d %d\n", i, q[0], q[1], q[2]);
        */
        CHECK_STATUS(get_point_block_file(block_file, &v_this, q, error_msg));
        /*
        printf("half_max_position4: i = %d, v_this = %3.2f, v_prev = %3.2f\n",
        i, v_this, v_prev);
        */

        if (find_maximum)
        {
            if (v_this < v_half)
            {
                *half_max = i - i_step*(v_half-v_this)/(v_prev-v_this);
                return CCPN_OK;
            }
        }
        else
        {
            if (v_this > v_half)
            {
                *half_max = i - i_step*(v_half-v_this)/(v_prev-v_this);
                return CCPN_OK;
            }
        }

        v_prev = v_this;
    }

    if (dirn == 1)
        *half_max = block_file->points[i] - 1.0;
    else
        *half_max = 1.0;

    return CCPN_OK;
}

static CcpnStatus check_dim_linewidth(Block_file block_file, Bool find_maximum,
                                      float min_linewidth, float v, int *point, int dim,
                                      Bool *ok_linewidth, CcpnString error_msg)
{
    float linewidth;

    /*
    printf("check_dim_linewidth1: v = %3.2f, dim = %d, min_linewidth = %3.2f, point = %d %d %d\n",
    v, dim, min_linewidth, point[0], point[1], point[2]);
    */
    if (linewidth_block_file(block_file, find_maximum, v, point, dim,
                             &linewidth, error_msg) == CCPN_ERROR)
        return CCPN_ERROR;

    /*
    printf("check_dim_linewidth2: linewidth = %3.2f\n", linewidth);
    */
    if (linewidth < min_linewidth)
        *ok_linewidth = CCPN_FALSE;
    /* do not need the below
        else
    	*ok_linewidth = CCPN_TRUE;
    */

    return CCPN_OK;
}

static CcpnStatus check_linewidth(Block_file block_file, Bool find_maximum,
                                  float *min_linewidth, float v, int *point, Bool *ok_linewidth,
                                  CcpnString error_msg)
{
    int i, ndim = block_file->ndim;

    *ok_linewidth = CCPN_TRUE;

    for (i = 0; i < ndim; i++)
    {
        if (min_linewidth[i] <= 0)
            continue;

        if (check_dim_linewidth(block_file, find_maximum, min_linewidth[i],
                                v, point, i, ok_linewidth, error_msg) == CCPN_ERROR)
            return CCPN_ERROR;

        if (*ok_linewidth == CCPN_FALSE)
            return CCPN_OK;
    }

    return CCPN_OK;
}

CcpnStatus find_peak_list(Peak_list peak_list, int *first, int *last,
                          Block_file block_file, Bool have_low, Bool have_high,
                          float low, float high, int *buffer, Bool nonadjacent,
                          float drop_factor, float *min_linewidth,
                          int ndiagonal_exclusions, Diagonal_exclusion diagonal_exclusions,
                          int nexcluded_regions, float ***excluded_regions,
                          Bool *dim_checked, int ignore_peak, CcpnString error_msg)
{
    int i, j, k, npoints, nadj_points, ndim = peak_list->ndim, dim1, dim2;
    int cum_points[MAX_NDIM], point[MAX_NDIM], cumulative[MAX_NDIM];
    float v, a1, a2, b12, d, delta;
    Bool find_maximum, ok_extreme, ok_drop, ok_linewidth;
    Diagonal_exclusion de;

    /*
        printf("find_peak_list: first = %d %d %d\n", first[0], first[1], first[2]);
        printf("find_peak_list: last = %d %d %d\n", last[0], last[1], last[2]);
        printf("find_peak_list: drop_factor = %3.2f, min_linewidth = %3.2f %3.2f %3.2f\n",
    	drop_factor, min_linewidth[0], min_linewidth[1], min_linewidth[2]);
    */
    if (!have_low && !have_high)
        return CCPN_OK;

    if (nonadjacent)
    {
        nadj_points = 1;
        for (i = 0; i < ndim; i++)
        {
            cumulative[i] = nadj_points;
            nadj_points *= 3;
        }
    }

    npoints = 1;
    for (i = 0; i < ndim; i++)
    {
        cum_points[i] = npoints;
        npoints *= last[i] - first[i];
    }

    for (i = 0; i < npoints; i++)
    {
        find_point(ndim, i, point, cum_points, first);

        if (ndiagonal_exclusions > 0)
        {
            for (j = 0; j < ndiagonal_exclusions; j++)
            {
                de = diagonal_exclusions + j;
                dim1 = de->dim1;
                dim2 = de->dim2;
                a1 = de->a1;
                a2 = de->a2;
                b12 = de->b12;
                d = de->d;

                delta = a1*point[dim1] - a2*point[dim2] + b12;
                if (ABS(delta) < d)
                    break;
            }

            if (j < ndiagonal_exclusions)
                continue; /* point is in some excluded diagonal */
        }

        if (nexcluded_regions > 0)
        {
            for (j = 0; j < nexcluded_regions; j++)
            {
                for (k = 0; k < ndim; k++)
                {
                    if ((point[k] < excluded_regions[j][k][0]) ||
                            (point[k] > excluded_regions[j][k][1]))
                        break; /* point is not in excluded region in this dim */
                }

                if (k == ndim)
                    break; /* point is in this excluded region */
            }

            if (j < nexcluded_regions)
                continue; /* point is in some excluded region */
        }

        CHECK_STATUS(get_point_block_file(block_file, &v, point, error_msg));

        if (have_high && (v >= high))
            find_maximum = CCPN_TRUE;
        else if (have_low && (v <= low))
            find_maximum = CCPN_FALSE;
        else
            continue;

        if (nonadjacent)
        {
            if (check_nonadjacent_points(block_file, find_maximum,
                                         buffer, v, point, nadj_points, cumulative,
                                         &ok_extreme, dim_checked, error_msg) == CCPN_ERROR)
                return CCPN_ERROR;
        }
        else
        {
            if (check_adjacent_points(block_file, find_maximum,
                                      buffer, v, point, &ok_extreme, dim_checked,
                                      error_msg) == CCPN_ERROR)
                return CCPN_ERROR;
        }

        if (!ok_extreme)
            continue;

        if (check_drop(block_file, find_maximum,
                       drop_factor, v, point, &ok_drop, error_msg) == CCPN_ERROR)
            return CCPN_ERROR;

        if (!ok_drop)
            continue;

        if (check_linewidth(block_file, find_maximum,
                            min_linewidth, v, point, &ok_linewidth, error_msg) == CCPN_ERROR)
            return CCPN_ERROR;

        if (!ok_linewidth)
            continue;

        CHECK_STATUS(check_new_peak(peak_list, buffer, point, ignore_peak, error_msg));
    }

    return CCPN_OK;
}

static CcpnStatus alloc_index_list_memory(int npeaks, int **index_list)
{
    MALLOC(*index_list, int, npeaks);

    return CCPN_OK;
}

CcpnStatus search_region_peak_list(Peak_list peak_list, float *first, float *last,
                                   Bool *allow_aliasing, int *npeaks, int **index_list, CcpnString error_msg)
{
    int i, n, npeaks_max = peak_list->npeaks;
    int *npoints = peak_list->npoints;
    Peak peak;

    if (npeaks_max == 0)
    {
        *npeaks = 0;
        *index_list = NULL;
        return CCPN_OK;
    }

    sprintf(error_msg, "allocating peak index list memory");
    CHECK_STATUS(alloc_index_list_memory(npeaks_max, index_list));

    n = 0;
    for (i = 0; i < npeaks_max; i++)
    {
        peak = peak_list->peaks[i];
        if (is_in_region_peak(peak, first, last, npoints, allow_aliasing, NULL))
            (*index_list)[n++] = i;
    }

    *npeaks = n;

    return CCPN_OK;
}

void search_nearest_peak_list(Peak_list peak_list, int xdim, int ydim,
                              float xscale, float yscale, float *first, float *last,
                              Bool *allow_aliasing, int *index_peak, float *d2_min)
{
    int i, npeaks = peak_list->npeaks, nearest = -1;
    int *npoints = peak_list->npoints;
    float d2, d2_array[MAX_NDIM], intensity_max, volume_max, scale;
    float xfirst, xlast, yfirst, ylast;
    Peak peak;

    *d2_min = 0; /* arbitrary */

    xfirst = first[xdim];
    xlast = last[xdim];
    yfirst = first[ydim];
    ylast = last[ydim];

    /*
    printf("search_nearest_peak_list: xdim = %d, ydim = %d\n", xdim, ydim);
    printf("search_nearest_peak_list: xfirst = %2.1f, xlast = %2.1f\n", xfirst, xlast);
    printf("search_nearest_peak_list: yfirst = %2.1f, ylast = %2.1f\n", yfirst, ylast);
    printf("search_nearest_peak_list: xscale = %2.1f, yscale = %2.1f\n", xscale, yscale);
    */

    determine_peak_list_max(peak_list, &intensity_max, &volume_max);

    for (i = 0; i < npeaks; i++)
    {
        first[xdim] = xfirst;
        last[xdim] = xlast;
        first[ydim] = yfirst;
        last[ydim] = ylast;

        peak = peak_list->peaks[i];
        /*
        	scale = determine_peak_scale(peak, intensity_max, volume_max);
        */
        scale = 1.0;
        find_scaled_region_peak(peak, xdim, ydim, scale*xscale, scale*yscale,
                                first, last);
        /*
        printf("search_nearest_peak_list: i = %d, scale = %2.1f\n", i, scale);
        printf("search_nearest_peak_list: first = %2.1f %2.1f %2.1f\n", first[0], first[1], first[2]);
        printf("search_nearest_peak_list: last = %2.1f %2.1f %2.1f\n", last[0], last[1], last[2]);
        printf("search_nearest_peak_list: position = %2.1f %2.1f %2.1f\n",
        peak->position[0]-1, peak->position[1]-1, peak->position[2]-1);
        */
        if (is_in_region_peak(peak, first, last, npoints, allow_aliasing, d2_array))
        {
            /*
            printf("search_nearest_peak_list: d2_array = %2.1f %2.1f %2.1f\n", d2_array[0], d2_array[1], d2_array[2]);
            */
            d2 = d2_array[xdim] + d2_array[ydim];
            if ((nearest < 0) || (d2 < *d2_min))
            {
                *d2_min = d2;
                nearest = i;
            }
        }
    }

    /*
    printf("search_nearest_peak_list: nearest = %d, d2_min = %2.1f\n", nearest, d2_min);
    */

    *index_peak = nearest;
}

float determine_peak_scale(Peak peak, float intensity_max, float volume_max)
{
    float z, zmax, scale;

    z = ABS(peak->volume);
    if (z == 0)
    {
        zmax = intensity_max;
        z = ABS(peak->intensity);
    }
    else
    {
        zmax = volume_max;
    }

    if (z == 0)
        scale = 0.5; /* arbitrary */
    else if (z >= zmax)
        scale = 1.0;
    else
        scale = 1.0 / (1.0 - log((double) (z/zmax))/log(10.0));

    return scale;
}

void determine_peak_list_max(Peak_list peak_list, float *intensity_max,
                             float *volume_max)
{
    int i;
    Peak peak;

    *intensity_max = *volume_max = 0;
    for (i = 0; i < peak_list->npeaks; i++)
    {
        peak = peak_list->peaks[i];
        *intensity_max = MAX(*intensity_max, ABS(peak->intensity));
        *volume_max = MAX(*volume_max, ABS(peak->volume));
    }
}

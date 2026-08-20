/*
======================COPYRIGHT/LICENSE START==========================

peak.c: Part of the CcpNmr Analysis program

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
#include "peak.h"

#include "method.h"
#include "nonlinear_model.h"
#include "symbol.h"

#define LARGE_NUMBER 1.0e20

/* don't use #define because of rounding */
static float DEFAULT_TEXT_OFFSET = -999.9;

typedef struct _FitPeak
{
    int ndim;
    int npeaks;
    int *region_offset;
    int *region_end;
    int *cumul_region;
    int method;
    Bool *dim_done;
}   FitPeak;

Peak new_peak(int ndim)
{
    int i;
    Peak peak;

    MALLOC_NEW(peak, struct Peak, 1);
    MALLOC_NEW(peak->position, float, ndim);
    ZERO_VECTOR(peak->position, ndim);
    MALLOC_NEW(peak->text_offset, float, ndim);
    for (i = 0; i < ndim; i++)
        peak->text_offset[i] = DEFAULT_TEXT_OFFSET;
    peak->value_offset = DEFAULT_TEXT_OFFSET;
    MALLOC_NEW(peak->num_aliasing, int, ndim);
    ZERO_VECTOR(peak->num_aliasing, ndim);
    MALLOC_NEW(peak->line_width, float, ndim);
    ZERO_VECTOR(peak->line_width, ndim);
    /*
        MALLOC_NEW(peak->box_width, float, ndim);
        ZERO_VECTOR(peak->box_width, ndim);
    */

    peak->ndim = ndim;
    MALLOC_NEW(peak->text, char, 1);
    *(peak->text) = 0;
    peak->isSelected = CCPN_FALSE;
    peak->intensity = 0.0;
    peak->volume = 0.0;
    /*
        peak->quality = 1.0;
    */

    return peak;
}

void delete_peak(Peak peak)
{
    if (peak)
    {
        FREE(peak->position, float);
        FREE(peak->text_offset, float);
        FREE(peak->num_aliasing, int);
        FREE(peak->text, char);
        FREE(peak->line_width, float);
        /*
        	FREE(peak->box_width, float);
        */
    }

    FREE(peak, struct Peak);
}

Bool get_is_selected_peak(Peak peak)
{
    return peak->isSelected;
}

void set_is_selected_peak(Peak peak, Bool isSelected)
{
    peak->isSelected = isSelected;
}

CcpnStatus set_text_peak(Peak peak, CcpnString text, CcpnString error_msg)
{
    FREE(peak->text, char);

    sprintf(error_msg, "allocating peak text memory");
    STRING_MALLOC_COPY(peak->text, text);

    return CCPN_OK;
}

void set_position_peak(Peak peak, float *position)
{
    COPY_VECTOR(peak->position, position, peak->ndim);
}

void set_text_offset_peak(Peak peak, int dim, float text_offset)
{
    if ((dim >= 0) && dim < (peak->ndim))
        peak->text_offset[dim] = text_offset;
    else if (dim == -1)
        peak->value_offset = text_offset;
}

void reset_text_offset_peak(Peak peak, int dim)
{
    if ((dim >= 0) && dim < (peak->ndim))
        peak->text_offset[dim] = DEFAULT_TEXT_OFFSET;
    else if (dim == -1)
        peak->value_offset = DEFAULT_TEXT_OFFSET;
}

void set_num_aliasing_peak(Peak peak, int *num_aliasing)
{
    COPY_VECTOR(peak->num_aliasing, num_aliasing, peak->ndim);
}

/* xpix, ypix are points per pixel */
/* x, y, xoff, yoff are in points */
static void draw_peak_text(Peak peak, float x, float y, float xoff, float yoff,
                           float xpix, float ypix, float xscale, float yscale,
                           Bool isPointerDrawn, Drawing_funcs *drawing_funcs, Generic_ptr data)
{
    static float rx = 1.5, ry = 1.5, tsize = 2.0;
    float a = 0, b = 0, d;
    float x0, y0, x1, y1, x2, y2;
    CcpnString text = peak->text;

    if (!text || !*text)
        return;

    if (xoff == DEFAULT_TEXT_OFFSET)
        xoff = rx * xscale;

    if (yoff == DEFAULT_TEXT_OFFSET)
        yoff = ry * yscale;

    a = 0.5 * (1 - xoff/5.0);
    a = MAX(0, MIN(1, a));

    b = 0.5 * (1 - yoff/5.0);
    b = MAX(0, MIN(1, b));

    (drawing_funcs->draw_text)(data, text, x+xoff, y+yoff, a, b);

    if (isPointerDrawn)
    {
        d = xoff*xoff/(xpix*xpix) + yoff*yoff/(ypix*ypix);

        /*
            	printf("HERE0: %3.2f %3.2f\n", xpix, ypix);
            	printf("HERE1: %3.2f %3.2f %3.2f\n", sqrt(d), a, b);
            	printf("HERE2: %3.2f %3.2f\n", xscale, yscale);
            	printf("HERE3: %3.2f %3.2f\n", xoff, yoff);
        */

        if (d > 500)
        {
            x0 = x;
            y0 = y;
            x += xoff;
            y += yoff;
            if (a < 0.05 || a > 0.95)
            {
                d = tsize * ypix;
                x1 = x;
                y1 = y+d;
                x2 = x;
                y2 = y-d;
            }
            else
            {
                d = tsize * xpix;
                x1 = x+d;
                y1 = y;
                x2 = x-d;
                y2 = y;
            }
            /*
            	    printf("HERE4: %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f\n", a, d, x0, y0, x1, y1, x2, y2);
            */
            (drawing_funcs->fill_triangle)(data, x0, y0, x1, y1, x2, y2);
        }
    }
}

void draw_peak(Peak peak, int xdim, int ydim, float *first, float *last,
               float xpix, float ypix, float xscale, float yscale,
               int symbol, Bool isSymbolDrawn, Bool isTextDrawn,
               Bool isTextPointerDrawn, Bool hasValueAxis,
               int *tile, Drawing_funcs *drawing_funcs, Generic_ptr data)
{
    int i;
    float x, y, xoff, yoff, dx, dy, r;
    Bool isAliased;

    if (!peak->position)
        return;

    /*
    printf("draw_peak: x = %3.2f, y = %3.2f\n", peak->position[xdim], peak->position[ydim]);
    printf("draw_peak: xoff = %3.2f, yoff = %3.2f\n", peak->text_offset[xdim], peak->text_offset[ydim]);
    printf("draw_peak: xscale = %3.2f, yscale = %3.2f\n", xscale, yscale);
    */

    x = peak->position[xdim] - 1;
    xoff = peak->text_offset[xdim];
    if (hasValueAxis)
    {
        /*
        printf("drawing peak line: x = %3.2f, intensity = %3.2f\n", x, peak->intensity);
        */
        yoff = peak->value_offset;
        if (xoff == DEFAULT_TEXT_OFFSET)
            dx = 0.03 * (last[0] - first[0]);
        else
            dx = xoff;
        if (yoff == DEFAULT_TEXT_OFFSET)
            dy = 0.06 * (last[1] - first[1]);
        else
            dy = yoff;
        if (isSymbolDrawn)
        {
            y = peak->intensity;
            /*
            	    (drawing_funcs->set_line_style)(data, DASHED_LINE_STYLE);
            */
            if (peak->isSelected)
                (drawing_funcs->set_line_width)(data, 2);
            (drawing_funcs->draw_line)(data, x, y, x+dx, y+dy);
            if (peak->isSelected)
                (drawing_funcs->set_line_width)(data, 1);
            /*
            	    (drawing_funcs->set_line_style)(data, NORMAL_LINE_STYLE);
            */
            /*  below draws line from 0 line to top of peak, does not look very good
            	    (drawing_funcs->draw_line)(data, x, 0.0, x, peak->intensity);
            */
        }

        if (isTextDrawn)
        {
            (drawing_funcs->draw_text)(data, peak->text, x+dx, y+dy, -0.05, 0.5);
        }
    }
    else
    {
        y = peak->position[ydim] - 1;
        yoff = peak->text_offset[ydim];

        /*
        printf("draw_peak: x = %2.1f, y = %2.1f\n", x, y);
        */
        if (isSymbolDrawn)
        {
            isAliased = CCPN_FALSE;
            for (i = 0; i < peak->ndim; i++)
            {
                if (tile[i] != peak->num_aliasing[i])
                {
                    isAliased = CCPN_TRUE;
                    break;
                }
            }

            draw_symbol(symbol, x, y, xscale, yscale, isAliased,
                        drawing_funcs, data);
            if (peak->isSelected)
                draw_symbol(BOX_SYMBOL, x, y, xscale, yscale, CCPN_FALSE,
                            drawing_funcs, data);
        }

        if (isTextDrawn)
            draw_peak_text(peak, x, y, xoff, yoff, xpix, ypix, xscale, yscale, isTextPointerDrawn, drawing_funcs, data);
    }
}

static CcpnStatus get_peak_point(Peak peak, int *point, int *points, CcpnString error_msg)
{
    int i;

    if (!peak->position)
        RETURN_ERROR_MSG("peak position not set");

    for (i = 0; i < peak->ndim; i++)
    {
        point[i] = peak->position[i] - 0.5;
        /* round to nearest integer but also -1 because position starts at 1 */

        /* TBD: note that below assumes that npoints = npoints_orig, which might be false */
        point[i] %= points[i];
        if (point[i] < 0)
            point[i] += points[i];
    }

    return CCPN_OK;
}

static CcpnStatus get_nearby_values(Peak peak, float *y, float *ym, float *yp,
                                    Block_file block_file, Bool *dim_done, CcpnString error_msg)
{
    int i, point[MAX_NDIM];

    CHECK_STATUS(get_peak_point(peak, point, block_file->points, error_msg));
    CHECK_STATUS(get_point_block_file(block_file, y, point, error_msg));

    for (i = 0; i < peak->ndim; i++)
    {
        if (!dim_done[i])
            continue;

        ym[i] = *y; /* safety in case npoints = 1 in this dim */

        if (point[i] > 0)
        {
            point[i]--;
            CHECK_STATUS(get_point_block_file(block_file, ym+i, point, error_msg));
            point[i]++;
        }

        if (point[i] < (block_file->points[i] - 1))
        {
            point[i]++;
            CHECK_STATUS(get_point_block_file(block_file, yp+i, point, error_msg));
            point[i]--;
        }
        else
        {
            yp[i] = ym[i]; /* not ideal */
        }

        if (point[i] == 0)
            ym[i] = yp[i]; /* not ideal */
    }

    return CCPN_OK;
}

CcpnStatus get_intensity_peak(Peak peak, float *intensity,
                              Block_file block_file, CcpnString error_msg)
{
    int point[MAX_NDIM];

    CHECK_STATUS(get_peak_point(peak, point, block_file->points, error_msg));
    CHECK_STATUS(get_point_block_file(block_file, intensity, point, error_msg));

    set_intensity_peak(peak, *intensity);

    return CCPN_OK;
}

void set_intensity_peak(Peak peak, float intensity)
{
    peak->intensity = intensity;
}

static Bool valid_parabolic_values(int ndim, float y, float *ym, float *yp, Bool *dim_done)
{
    int i;

    if (y == 0)
        return CCPN_FALSE;

    if (y > 0)
    {
        for (i = 0; i < ndim; i++)
        {
            if (dim_done[i] && ((ym[i] > y) || (yp[i] > y)))
                return CCPN_FALSE;
        }
    }
    else
    {
        for (i = 0; i < ndim; i++)
        {
            if (dim_done[i] && ((ym[i] < y) || (yp[i] < y)))
                return CCPN_FALSE;
        }
    }

    return CCPN_TRUE;
}

static CcpnStatus fit_volume_gaussian3(Peak peak, float *volume, Bool *valid,
                                       Block_file block_file, Bool *dim_done, CcpnString error_msg)
{
    float y, ym[MAX_NDIM], yp[MAX_NDIM];

    CHECK_STATUS(get_nearby_values(peak, &y, ym, yp, block_file, dim_done, error_msg));
    /*  13 Aug 2004: remove check because fiddled in fit_volume_gaussian3_method now
        if (!valid_parabolic_values(peak->ndim, y, ym, yp, dim_done))
        {
    	*valid = CCPN_FALSE;
    	return CCPN_OK;
        }
    */

    *volume = fit_volume_gaussian3_method(peak->ndim, y, ym, yp, dim_done);

    return CCPN_OK;
}

static CcpnStatus fit_center_parabolic(Peak peak, float *center, Bool *valid,
                                       Block_file block_file, Bool *dim_done, CcpnString error_msg)
{
    float y, ym[MAX_NDIM], yp[MAX_NDIM];

    CHECK_STATUS(get_nearby_values(peak, &y, ym, yp, block_file, dim_done, error_msg));
    if (!valid_parabolic_values(peak->ndim, y, ym, yp, dim_done))
    {
        *valid = CCPN_FALSE;
        return CCPN_OK;
    }

    fit_center_parabolic_method(peak->ndim, center, y, ym, yp, dim_done);

    return CCPN_OK;
}

static CcpnStatus fit_center_gaussian3(Peak peak, float *center, Bool *valid,
                                       Block_file block_file, Bool *dim_done, CcpnString error_msg)
{
    float y, ym[MAX_NDIM], yp[MAX_NDIM];

    CHECK_STATUS(get_nearby_values(peak, &y, ym, yp, block_file, dim_done, error_msg));
    if (!valid_parabolic_values(peak->ndim, y, ym, yp, dim_done))
    {
        *valid = CCPN_FALSE;
        return CCPN_OK;
    }

    fit_center_gaussian3_method(peak->ndim, center, y, ym, yp, dim_done);

    return CCPN_OK;
}

CcpnStatus fit_volume_peak(Peak peak, float *volume, Bool *valid, int method,
                           Block_file block_file, Bool *dim_done, CcpnString error_msg)
{
    *valid = CCPN_TRUE;

    if (method == FIT_VOLUME_GAUSSIAN3_METHOD)
    {
        CHECK_STATUS(fit_volume_gaussian3(peak, volume, valid, block_file, dim_done, error_msg));
    }
    else
    {
        RETURN_ERROR_MSG("illegal fit method");
    }

    if (*valid)
        set_volume_peak(peak, *volume);

    return CCPN_OK;
}

void set_volume_peak(Peak peak, float volume)
{
    peak->volume = volume;
}

CcpnStatus fit_center_peak(Peak peak, float *center, Bool *valid, int method,
                           Block_file block_file, Bool *dim_done, CcpnString error_msg)
{
    *valid = CCPN_TRUE;

    if (method == FIT_CENTER_PARABOLIC_METHOD)
    {
        CHECK_STATUS(fit_center_parabolic(peak, center, valid, block_file, dim_done, error_msg));
    }
    else if (method == FIT_CENTER_GAUSSIAN3_METHOD)
    {
        CHECK_STATUS(fit_center_gaussian3(peak, center, valid, block_file, dim_done, error_msg));
    }
    else
    {
        RETURN_ERROR_MSG("illegal fit method");
    }

    return CCPN_OK;
}

Bool is_in_region_peak(Peak peak, float *first, float *last,
                       int *npoints, Bool *allow_aliasing, float *d2_array)
{
    int i, m;
    float p, d;

    if (!peak->position)
        return CCPN_FALSE;

    for (i = 0; i < peak->ndim; i++)
    {
        p = peak->position[i] - 1;

        if (allow_aliasing[i])
        {
            m = (int) floor((p - first[i]) / npoints[i]);
            p -= m * npoints[i];
            if (p > last[i])
                return CCPN_FALSE;
        }
        else
        {
            p += peak->num_aliasing[i] * npoints[i];
            if ((p < first[i]) || (p > last[i]))
                return CCPN_FALSE;
        }

        if (d2_array)
        {
            d = p - 0.5 * (first[i] + last[i]);
            d2_array[i] = d * d;
        }
    }

    return CCPN_TRUE;
}

void find_scaled_region_peak(Peak peak, int xdim, int ydim,
                             float xscale, float yscale, float *first, float *last)
{
    first[xdim] -= xscale;
    first[ydim] -= yscale;
    last[xdim] += xscale;
    last[ydim] += yscale;
}

CcpnStatus fit_linewidth_peak(Peak peak, float *linewidth, Block_file block_file,
                              Bool *dim_done, CcpnString error_msg)
{
    int i, point[MAX_NDIM];
    float v;
    Bool find_maximum;

    CHECK_STATUS(get_peak_point(peak, point, block_file->points, error_msg));

    /* could use peak->intensity instead */
    CHECK_STATUS(get_point_block_file(block_file, &v, point, error_msg));

    /* slightly dangerous, this assumption */
    find_maximum = (v > 0);

    for (i = 0; i < peak->ndim; i++)
    {
        if (!dim_done[i])
        {
            linewidth[i] = 0;
            continue;
        }

        if (linewidth_block_file(block_file, find_maximum, v, point, i,
                                 linewidth+i, error_msg) == CCPN_ERROR)
            return CCPN_ERROR;
    }

    return CCPN_OK;
}

/*
void set_quality_peak(Peak peak, float quality)
{
    peak->quality = quality;
}
*/

static float gaussian(int ndim, int ndim_done, Bool *dim_done, int *x, float *a, float *dy_da)
{
    float h = a[0], *position = a+1, *linewidth = a+1+ndim_done;
    float *dy_dh, *dy_dp, *dy_dl;
    float lw, dx, y = h;
    int i, j;

    if (dy_da)
    {
        dy_dh = dy_da;
        dy_dp = dy_da+1;
        dy_dl = dy_da+1+ndim_done;
    }

    j = 0;
    for (i = 0; i < ndim; i++)
    {
        if (dim_done[i])
        {
            dx = x[i] - position[j];
            lw = linewidth[j];
            y *= exp(-4*log(2)*dx*dx/(lw*lw));
            if (dy_da)
            {
                dy_dp[j] = 8*log(2)*dx/(lw*lw);
                dy_dl[j] = 8*log(2)*dx*dx/(lw*lw*lw);
            }

            j++;
        }
    }

    if (dy_da)
    {
        *dy_dh = y / h;
        SCALE_VECTOR(dy_dp, dy_dp, y, ndim_done);
        SCALE_VECTOR(dy_dl, dy_dl, y, ndim_done);
    }

    return y;
}

static float lorentzian(int ndim, int ndim_done, Bool *dim_done, int *x, float *a, float *dy_da)
{
    float h = a[0], *position = a+1, *linewidth = a+1+ndim_done;
    float *dy_dh, *dy_dp, *dy_dl;
    float lw, dx, d, y = h;
    int i, j;

    if (dy_da)
    {
        dy_dh = dy_da;
        dy_dp = dy_da+1;
        dy_dl = dy_da+1+ndim_done;
    }

    j = 0;
    for (i = 0; i < ndim; i++)
    {
        if (dim_done[i])
        {
            dx = x[i] - position[j];
            lw = linewidth[j];
            d = lw*lw+4*dx*dx;
            y *= lw*lw/d;
            dy_dp[j] = 8*dx/d;
            dy_dl[j] = 8*dx*dx/(lw*d);

            j++;
        }
    }

    if (dy_da)
    {
        *dy_dh = y / h;
        SCALE_VECTOR(dy_dp, dy_dp, y, ndim_done);
        SCALE_VECTOR(dy_dl, dy_dl, y, ndim_done);
    }

    return y;
}

static void _fitting_func(float xind, float *a, float *y_fit, float *dy_da, void *user_data)
{
    FitPeak *fitPeak = (FitPeak *) user_data;
    int ndim = fitPeak->ndim;
    int npeaks = fitPeak->npeaks;
    int *region_offset = fitPeak->region_offset;
    int *region_end = fitPeak->region_end;
    int *cumul_region = fitPeak->cumul_region;
    int method = fitPeak->method;
    Bool *dim_done = fitPeak->dim_done;
    int ind = NEAREST_INTEGER(xind);
    int nparams_per_peak, ndim_done;
    int i, j, x[MAX_NDIM];
    float *position = a+1;

    ARRAY_OF_INDEX(x, ind, cumul_region, ndim);
    ADD_VECTORS(x, x, region_offset, ndim);

    // check whether position is outside the intended fitting region
    for (i = 0; i < ndim; i++)
    {
        if (dim_done[i])
        {
            if ((position[i] < region_offset[i]) || (position[i] >= region_end[i]))
            {
                *y_fit = LARGE_NUMBER;
                ZERO_VECTOR(dy_da, npeaks*nparams_per_peak); // arbitrary, hopefully ok
                return;
            }
        }
    }

    nparams_per_peak = 1;
    ndim_done = 0;
    for (i = 0; i < ndim; i++)
    {
        if (dim_done[i])
        {
            ndim_done++;
            nparams_per_peak += 2;
        }
    }

    *y_fit = 0;
    for (j = 0; j < npeaks; j++)
    {
        if (method == GAUSSIAN_METHOD)
            *y_fit += gaussian(ndim, ndim_done, dim_done, x, a, dy_da);
        else
            *y_fit += lorentzian(ndim, ndim_done, dim_done, x, a, dy_da);

        a += nparams_per_peak;
        if (dy_da)
            dy_da += nparams_per_peak;
    }
}

CcpnStatus fit_peaks_in_region(int npeaks, Peak *peaks, int method,
                               int *first, int *last, Block_file block_file, Bool *dim_done,
                               float *params, float *params_dev, CcpnString error_msg)
{
    int i, j, k, n, ndim, total_region_size, nparams, np;
    int region_size[MAX_NDIM], cumul_region[MAX_NDIM], array[MAX_NDIM];
    float *x, *y, *w = NULL, *y_fit = NULL;
    float chisq, scale = 0, max_iter = 0, noise = 0;
    float linewidth[MAX_NDIM];
    Peak peak;
    FitPeak fitPeak;
    CcpnStatus status;

    if (npeaks < 1)
        return CCPN_OK;

    if ((method != GAUSSIAN_METHOD) && (method != LORENTZIAN_METHOD))
        RETURN_ERROR_MSG("method must be GAUSSIAN_METHOD or LORENTZIAN_METHOD");

    /* it is assumed that all peaks come from the same spectrum and so have the
       same ndim and also the same block_file */

    ndim = peaks[0]->ndim;

    for (i = 0; i < ndim; i++)
    {
        /*
        printf("first[%d] = %d\n", i, first[i]);
        */
        if (first[i] < 0)
        {
            sprintf(error_msg, "first[%d] = %d < 0", i, first[i]);
            return CCPN_ERROR;
        }

        if (first[i] >= last[i])
        {
            sprintf(error_msg, "first[%d] = %d < %d = last[%d]", i, first[i], last[i], i);
            return CCPN_ERROR;
        }

        if (last[i] > block_file->points[i])
        {
            sprintf(error_msg, "last[%d] = %d > %d = npts", i, last[i], block_file->points[i]);
            return CCPN_ERROR;
        }
    }

    for (i = 0; i < ndim; i++)
        region_size[i] = last[i] - first[i];

    CUMULATIVE(cumul_region, region_size, total_region_size, ndim);

    sprintf(error_msg, "allocating memory for x, y");
    MALLOC(x, float, total_region_size);
    MALLOC(y, float, total_region_size);

    for (j = 0; j < total_region_size; j++)
    {
        x[j] = j;  // the real x is multidimensional so have to use index into it
        ARRAY_OF_INDEX(array, j, cumul_region, ndim);
        ADD_VECTORS(array, array, first, ndim);
        CHECK_STATUS(get_point_block_file(block_file, y+j, array, error_msg));
        scale = MAX(scale, ABS(y[j]));
        /*
        printf("Y[%d] = %f\n", j, y[j]);
        */
    }

    if (scale <= 0)
        scale = 1.0;
    /* scale the intensities down */
    SCALE_VECTOR(y, y, 1.0/scale, total_region_size);

    nparams = 0;
    for (j = 0; j < npeaks; j++)
    {
        peak = peaks[j];
        params[nparams++] = peak->intensity / scale;
        for (i = 0; i < ndim; i++)
        {
            if (dim_done[i])
                params[nparams++] = peak->position[i] - 1;
        }
        CHECK_STATUS(fit_linewidth_peak(peak, linewidth, block_file, dim_done, error_msg));
        for (i = 0; i < ndim; i++)
        {
            if (dim_done[i])
                params[nparams++] = linewidth[i];
        }
    }
    /*
    for (i = 0; i < nparams; i++)
    printf("PARAM %d before: %f\n", i, params[i]);
    */

    fitPeak.ndim = ndim;
    fitPeak.npeaks = npeaks;
    fitPeak.region_offset = first;
    fitPeak.region_end = last;
    fitPeak.cumul_region = cumul_region;
    fitPeak.method = method;
    fitPeak.dim_done = dim_done;

    status = nonlinear_fit(total_region_size, x, y, w, y_fit,
                           nparams, params, params_dev,
                           max_iter, noise, &chisq,
                           _fitting_func, (void *) &fitPeak, error_msg);

    np = 0;
    for (j = 0; j < npeaks; j++)
    {
        params[np++] *= scale;
        for (i = 0; i < ndim; i++)
        {
            if (dim_done[i])
                params[np++] += 1;
        }
        for (i = 0; i < ndim; i++)
        {
            if (dim_done[i])
                np++;
        }
    }
    /*
    for (i = 0; i < nparams; i++)
    printf("PARAM %d after: %f\n", i, params[i]);
    */

    FREE(x, float);
    FREE(y, float);

    return status;
}

void set_line_width_peak(Peak peak, int dim, float line_width)
{
    if ((dim >= 0) && dim < (peak->ndim))
        peak->line_width[dim] = line_width;
}

/*
void set_box_width_peak(Peak peak, int dim, float box_width)
{
    if ((dim >= 0) && dim < (peak->ndim))
        peak->box_width[dim] = box_width;
}
*/


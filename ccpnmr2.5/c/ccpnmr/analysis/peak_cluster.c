/*
======================COPYRIGHT/LICENSE START==========================

peak_cluster.c: Part of the CcpNmr Analysis program

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
#include "peak_cluster.h"

#define NPEAKS_ALLOC  10

static CcpnStatus check_alloc(Peak_cluster peak_cluster, CcpnString error_msg)
{
    int npeaks = peak_cluster->npeaks;
    int nalloc = peak_cluster->nalloc;

    sprintf(error_msg, "allocating peak_cluster peak memory");
    if (nalloc == 0)
    {
        MALLOC(peak_cluster->peaks, Peak, NPEAKS_ALLOC);
        peak_cluster->nalloc = NPEAKS_ALLOC;
    }
    else if (npeaks >= nalloc)
    {
        nalloc += NPEAKS_ALLOC;
        REALLOC(peak_cluster->peaks, Peak, nalloc);
        peak_cluster->nalloc = nalloc;
    }

    return CCPN_OK;
}

Peak_cluster new_peak_cluster(int ndim, int cluster_type)
{
    int i;
    Peak_cluster peak_cluster;

    MALLOC_NEW(peak_cluster, struct Peak_cluster, 1);

    peak_cluster->ndim = ndim;
    peak_cluster->cluster_type = cluster_type;
    peak_cluster->npeaks = 0;
    peak_cluster->nalloc = 0;
    peak_cluster->peaks = NULL;

    MALLOC_NEW(peak_cluster->text, char, ndim);
    *(peak_cluster->text) = 0;

    MALLOC_NEW(peak_cluster->dim_text, CcpnString, ndim);
    for (i = 0; i < ndim; i++)
    {
        MALLOC_NEW(peak_cluster->dim_text[i], char, 1);
        *(peak_cluster->dim_text[i]) = 0;
    }

    return peak_cluster;
}

void delete_peak_cluster(Peak_cluster peak_cluster)
{
    int i;

    if (peak_cluster)
    {
        for (i = 0; i < peak_cluster->ndim; i++)
            FREE(peak_cluster->dim_text[i], char);
        FREE(peak_cluster->dim_text, CcpnString);

        FREE(peak_cluster->peaks, Peak);
    }

    FREE(peak_cluster, struct Peak_cluster);
}

CcpnStatus add_peak_peak_cluster(Peak_cluster peak_cluster, Peak peak, CcpnString error_msg)
{
    if (peak->ndim != peak_cluster->ndim)
    {
        sprintf(error_msg, "peak_cluster ndim = %d != %d = peak ndim", peak_cluster->ndim, peak->ndim);
        return CCPN_ERROR;
    }

    CHECK_STATUS(check_alloc(peak_cluster, error_msg));

    peak_cluster->peaks[peak_cluster->npeaks++] = peak;

    return CCPN_OK;
}

void remove_peak_peak_cluster(Peak_cluster peak_cluster, Peak peak)
{
    int i, j, npeaks = peak_cluster->npeaks;
    Peak *peaks = peak_cluster->peaks;

    for (i = 0; i < npeaks; i++)
    {
        if (peaks[i] == peak)
        {
            for (j = i; j < npeaks-1; j++)
                peaks[j] = peaks[j+1];
            peak_cluster->npeaks--;
            break;
        }
    }
}

void set_peaks_peak_cluster(Peak_cluster peak_cluster, int npeaks, int nalloc, Peak *peaks)
{
    FREE(peak_cluster->peaks, Peak);

    peak_cluster->npeaks = npeaks;
    peak_cluster->nalloc = nalloc;
    peak_cluster->peaks = peaks;
}

void clear_peaks_peak_cluster(Peak_cluster peak_cluster)
{
    /*  don't free memory, just reset counter  */
    peak_cluster->npeaks = 0;
}

CcpnStatus set_text_peak_cluster(Peak_cluster peak_cluster, CcpnString text, CcpnString error_msg)
{
    FREE(peak_cluster->text, char);

    sprintf(error_msg, "allocating peak_cluster text memory");
    STRING_MALLOC_COPY(peak_cluster->text, text);

    return CCPN_OK;
}

CcpnStatus set_dim_text_peak_cluster(Peak_cluster peak_cluster, int dim, CcpnString text, CcpnString error_msg)
{
    if ((dim < 0) || (dim >= peak_cluster->ndim))
    {
        sprintf(error_msg, "dim = %d, must be between 0 and %d", dim, peak_cluster->ndim-1);
        return CCPN_ERROR;
    }

    FREE(peak_cluster->dim_text[dim], char);

    sprintf(error_msg, "allocating peak_cluster text memory");
    STRING_MALLOC_COPY(peak_cluster->dim_text[dim], text);

    return CCPN_OK;
}

void draw_peak_cluster(Peak_cluster peak_cluster, int xdim, int ydim,
                       Drawing_funcs *drawing_funcs, Generic_ptr data)
{
    int i, ndim = peak_cluster->ndim, npeaks = peak_cluster->npeaks;
    float xmin, xmax, ymin, ymax, p, x, y, dx, dy, offset;
    CcpnString text;
    Peak *peaks = peak_cluster->peaks;

    if (npeaks == 0)
        return;

    if ((xdim < 0) || (xdim >= ndim))
        return;

    if ((ydim < 0) || (ydim >= ndim))
        return;

    xmin = xmax = peaks[0]->position[xdim];
    ymin = ymax = peaks[0]->position[ydim];
    for (i = 1; i < npeaks; i++)
    {
        p = peaks[i]->position[xdim];
        xmin = MIN(xmin, p);
        xmax = MAX(xmax, p);
        p = peaks[i]->position[ydim];
        ymin = MIN(ymin, p);
        ymax = MAX(ymax, p);
    }

    (drawing_funcs->draw_line)(data, xmin, ymin, xmax, ymin);
    (drawing_funcs->draw_line)(data, xmax, ymin, xmax, ymax);
    (drawing_funcs->draw_line)(data, xmax, ymax, xmin, ymax);
    (drawing_funcs->draw_line)(data, xmin, ymax, xmin, ymin);

    x = HALF * (xmin + xmax);
    y = HALF * (ymin + ymax);
    dx = xmax - xmin;
    dy = ymax - ymin;

    text = peak_cluster->dim_text[xdim];
    if (*text)
    {
        /* TBD: offset?? */
        (drawing_funcs->draw_text)(data, text, x, ymax, 0.0, 0.0);
    }

    text = peak_cluster->dim_text[ydim];
    if (*text)
    {
        /* TBD: offset?? */
        (drawing_funcs->draw_text)(data, text, xmax, y, 0.0, 0.0);
    }

    text = peak_cluster->text;
    if (*text)
    {
        /* TBD: offset?? */
        (drawing_funcs->draw_text)(data, text, x, y, 0.0, 0.0);
    }
}


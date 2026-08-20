/*
======================COPYRIGHT/LICENSE START==========================

win_peak_list.c: Part of the CcpNmr Analysis program

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
#include "win_peak_list.h"

#include "symbol.h"

#define  NCOLORS  3

Win_peak_list new_win_peak_list(Peak_list peak_list, Bool hasValueAxis)
{
    Win_peak_list win_peak_list;

    MALLOC_NEW(win_peak_list, struct Win_peak_list, 1);

    /*  win_peak_list does not own peak_list  */
    win_peak_list->peak_list = peak_list;
    win_peak_list->hasValueAxis = hasValueAxis;
    win_peak_list->isSymbolDrawn = CCPN_TRUE;
    win_peak_list->isTextDrawn = CCPN_TRUE;
    win_peak_list->isTextPointerDrawn = CCPN_TRUE;

    return win_peak_list;
}

void delete_win_peak_list(Win_peak_list win_peak_list)
{
    FREE(win_peak_list, struct Win_peak_list);
}

void is_symbol_drawn_win_peak_list(Win_peak_list win_peak_list, Bool isSymbolDrawn)
{
    win_peak_list->isSymbolDrawn = isSymbolDrawn;
}

void is_text_drawn_win_peak_list(Win_peak_list win_peak_list, Bool isTextDrawn)
{
    win_peak_list->isTextDrawn = isTextDrawn;
}

void is_text_pointer_drawn_win_peak_list(Win_peak_list win_peak_list, Bool isTextPointerDrawn)
{
    win_peak_list->isTextPointerDrawn = isTextPointerDrawn;
}

static Bool is_peak_drawn(int ndim, Peak peak, float *first, float *last)
{
    int i;
    float p;

    if (!peak->position)
        return CCPN_FALSE;

    /*
    printf("is_peak_drawn: ndim=%d\n", ndim);
    printf("is_peak_drawn: position=%2.1f %2.1f %2.1f\n",
    peak->position[0], peak->position[1], peak->position[2]);
    printf("is_peak_drawn: first=%2.1f %2.1f %2.1f\n", first[0], first[1], first[2]);
    printf("is_peak_drawn: last=%2.1f %2.1f %2.1f\n", last[0], last[1], last[2]);
    */

    /* TBD: aliasing of peak position */

    for (i = 0; i < ndim; i++)
    {
        p = peak->position[i] - 1;
        /*
        printf("is_peak_drawn: i=%d, p=%2.1f, first=%d, last=%d\n", i, p, first[i], last[i]);
        */
        if ((p < first[i]) || (p >= last[i]))
            return CCPN_FALSE;
    }

    return CCPN_TRUE;
}

static void invert_color(float *inverted_color, float *color)
{
    int i;

    for (i = 0; i < NCOLORS; i++)
        inverted_color[i] = 1 - color[i];
}

static void depth_cue_peak(Peak peak, int xdim, int ydim,
                           int *npoints, float *center, float *thickness,
                           int *tile, float *fg_color, float *bg_color,
                           Drawing_funcs *drawing_funcs, Generic_ptr data)
{
    int i;
    float p, s = 1;
    float color[NCOLORS];

    for (i = 0; i < peak->ndim; i++)
    {
        /* thickness = 0 for pseudo-3D */
        if ((i != xdim) && (i != ydim) && (thickness[i] > 0))
        {
            p = peak->position[i] + tile[i] * npoints[i];
            s *= 1.0 - MIN(1.0, ABS(p - center[i]) / thickness[i]);
            /*
            printf("depth_cue_peak1: i = %d, s = %4.3f, pos = %4.3f, tile = %d, p = %4.3f, center = %4.3f\n", i, s, peak->position[i], tile[i], p, center[i]);
            */
        }
    }

    /* before: s = 1 in center, 0 at edge */
    /* after: s = 0 in center, 1 at edge */
    s = 1 - s;

    s *= 0.5;
    /*
    printf("depth_cue_peak2: s = %4.3f\n", s);
    */

    for (i = 0; i < NCOLORS; i++)
        color[i] = (1-s)*fg_color[i] + s*bg_color[i];

    (drawing_funcs->set_draw_color)(data, color);
}

CcpnStatus draw_win_peak_list(Win_peak_list win_peak_list, int xdim, int ydim,
                              float xpix, float ypix,  float xscale, float yscale,
                              float *first, float *last,
                              int draw_method, float intensity_max, float volume_max,
                              float *bg_color, float *center, float *thickness, int *tile,
                              Drawing_funcs *drawing_funcs, Generic_ptr data,
                              CcpnString error_msg)
{
    int i, font_size = 10;
    Peak peak;
    float scale = 1.0;
    float *fg_color, inverted_bg_color[NCOLORS];
    float xsc, ysc;
    Peak_list peak_list = win_peak_list->peak_list;

    invert_color(inverted_bg_color, bg_color);

    if (!win_peak_list->isSymbolDrawn && !win_peak_list->isTextDrawn)
        return CCPN_OK;

    if (peak_list->npeaks == 0)
        return CCPN_OK;

    /*
    printf("draw_win_peak_list: npeaks = %d, xdim=%d, ydim=%d, xscale=%f, yscale=%f\n",
    peak_list->npeaks, xdim, ydim, xscale, yscale);
    printf("draw_win_peak_list: first=%2.1f %2.1f %2.1f\n", first[0], first[1], first[2]);
    printf("draw_win_peak_list: last=%2.1f %2.1f %2.1f\n", last[0], last[1], last[2]);
    */

    if (draw_method == DRAW_PEAKLIST_METHOD)
        determine_peak_list_max(peak_list, &intensity_max, &volume_max);

    /*
    printf("draw_win_peak_list: intensity_max=%f, volume_max=%f\n",
    intensity_max, volume_max);
    */

    (drawing_funcs->set_draw_color)(data, peak_list->color);
    /*
        (drawing_funcs->set_draw_font)(data, "Times-Roman", font_size);
        if (draw_method == DRAW_UNIFORM_METHOD)
            font_size = 10;
        else
            font_size = CEILING(10 * MAX(xscale, yscale));
    */
    /*
        (drawing_funcs->set_draw_font)(data, "Helvetica", 10);
    */

    xsc = xscale;
    ysc = yscale;
    for (i = 0; i < peak_list->npeaks; i++)
    {
        peak = peak_list->peaks[i];
        /*
        printf("draw_win_peak_list: peak %d position =%2.1f %2.1f %2.1f\n", i,
        peak->position[0], peak->position[1], peak->position[2]);
        */
        if (is_peak_drawn(peak_list->ndim, peak, first, last))
        {
            if (win_peak_list->hasValueAxis && peak->isSelected)
                fg_color = inverted_bg_color;
            else
                fg_color = peak_list->color;

            if (draw_method == DRAW_LINE_WIDTH_METHOD)
            {
                xsc = xscale * peak->line_width[xdim];
                ysc = yscale * peak->line_width[ydim];
                /*
                printf("draw_win_peak_list: peak = %d, xsc=%f, ysc=%f\n", i, xsc, ysc);
                printf("draw_win_peak_list: xlinewidth=%f, ylinewidth=%f\n", peak->line_width[xdim], peak->line_width[ydim]);
                */
            }
            else if (draw_method != DRAW_UNIFORM_METHOD)
            {
                scale = determine_peak_scale(peak, intensity_max, volume_max);
                xsc = scale * xscale;
                ysc = scale * yscale;
            }
            /*
            printf("draw_win_peak_list: peak = %d, scale=%f\n", i, scale);
            */
            if (win_peak_list->hasValueAxis)
                (drawing_funcs->set_draw_color)(data, fg_color);
            else
                depth_cue_peak(peak, xdim, ydim,
                               peak_list->npoints, center, thickness, tile,
                               fg_color, bg_color, drawing_funcs, data);

            draw_peak(peak, xdim, ydim, first, last,
                      xpix, ypix, xsc, ysc,
                      peak_list->symbol, win_peak_list->isSymbolDrawn,
                      win_peak_list->isTextDrawn, win_peak_list->isTextPointerDrawn,
                      win_peak_list->hasValueAxis, tile, drawing_funcs, data);
        }
    }

    return CCPN_OK;
}

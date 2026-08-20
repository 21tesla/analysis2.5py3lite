/*
======================COPYRIGHT/LICENSE START==========================

symbol.c: Part of the CcpNmr Analysis program

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
#include "symbol.h"

static void draw_cross(float x, float y, float xscale, float yscale,
                       Drawing_funcs *drawing_funcs, Generic_ptr data)
{
    static float r = 1;
    float dx = r*xscale, dy = r*yscale;

    (drawing_funcs->draw_line)(data, x-dx, y-dy, x+dx, y+dy);
    (drawing_funcs->draw_line)(data, x-dx, y+dy, x+dx, y-dy);
}

static void draw_plus(float x, float y, float xscale, float yscale,
                      Drawing_funcs *drawing_funcs, Generic_ptr data)
{
    static float r = 1;
    float dx = r*xscale, dy = r*yscale;

    (drawing_funcs->draw_line)(data, x, y-dy, x, y+dy);
    (drawing_funcs->draw_line)(data, x-dx, y, x+dx, y);
}

static void draw_circle(float x, float y, float xscale, float yscale,
                        Drawing_funcs *drawing_funcs, Generic_ptr data)
{
    static float r = 1;
    float dx = r*xscale, dy = r*yscale;
    /*
        float dr = MIN(dx, dy);

        (drawing_funcs->draw_circle)(data, x, y, dr);
    */
    (drawing_funcs->draw_ellipse)(data, x, y, dx, dy);
}

static void draw_disk(float x, float y, float xscale, float yscale,
                      Drawing_funcs *drawing_funcs, Generic_ptr data)
{
    static float r = 1;
    float dx = r*xscale, dy = r*yscale;
    /*
        float dr = MIN(dx, dy);

        (drawing_funcs->fill_circle)(data, x, y, dr);
    */
    (drawing_funcs->fill_ellipse)(data, x, y, dx, dy);
}

static void draw_box(float x, float y, float xscale, float yscale,
                     Drawing_funcs *drawing_funcs, Generic_ptr data)
{
    static float r = 1;
    float dx = r*xscale, dy = r*yscale;

    (drawing_funcs->draw_line)(data, x-dx, y-dy, x+dx, y-dy);
    (drawing_funcs->draw_line)(data, x+dx, y-dy, x+dx, y+dy);
    (drawing_funcs->draw_line)(data, x+dx, y+dy, x-dx, y+dy);
    (drawing_funcs->draw_line)(data, x-dx, y+dy, x-dx, y-dy);
}

void draw_symbol(int symbol, float x, float y,
                 float xscale, float yscale, Bool isAliased,
                 Drawing_funcs *drawing_funcs, Generic_ptr data)
{
    if (isAliased)
        (drawing_funcs->set_line_style)(data, DASHED_LINE_STYLE);

    if (symbol == CROSS_SYMBOL)
        draw_cross(x, y, xscale, yscale, drawing_funcs, data);
    else if (symbol == PLUS_SYMBOL)
        draw_plus(x, y, xscale, yscale, drawing_funcs, data);
    else if (symbol == CIRCLE_SYMBOL)
        draw_circle(x, y, xscale, yscale, drawing_funcs, data);
    else if (symbol == DISK_SYMBOL)
        draw_disk(x, y, xscale, yscale, drawing_funcs, data);
    else if (symbol == BOX_SYMBOL)
        draw_box(x, y, xscale, yscale, drawing_funcs, data);

    if (isAliased)
        (drawing_funcs->set_line_style)(data, NORMAL_LINE_STYLE);
}

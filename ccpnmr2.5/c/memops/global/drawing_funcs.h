
/*
======================COPYRIGHT/LICENSE START==========================

drawing_funcs.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_drawing_funcs
#define _incl_drawing_funcs

#include "macros.h"
#include "types.h"

#include "polyline.h"

#define  SCREEN_DRAWING		0
#define  PAPER_DRAWING		1
#define  STORE_DRAWING		2

#define  NORMAL_LINE_STYLE	0
#define  DASHED_LINE_STYLE	1
#define  NLINE_STYLES           2

/* it's possible do not need start_draw, end_draw, new_draw_range */
/* relevant information can now be stored in data */
typedef struct
{
    int draw_medium;
    void (*start_draw)(Generic_ptr data);
    void (*end_draw)(Generic_ptr data);
    void (*new_draw_range)(Generic_ptr data, float x0, float y0, float x1, float y1, Bool clip);
    void (*draw_line)(Generic_ptr data, float x0, float y0, float x1, float y1);
    void (*draw_clipped_line)(Generic_ptr data, float x0, float y0, float x1, float y1);
    void (*draw_polyline)(Generic_ptr data, Poly_line polyline);
    void (*draw_clipped_polyline)(Generic_ptr data, Poly_line polyline);
    void (*draw_text)(Generic_ptr data, CcpnString text, float x, float y, float a, float b);
    void (*set_draw_color)(Generic_ptr data, float *color);
    void (*set_draw_font)(Generic_ptr data, CcpnString name, int size);
    void (*set_line_style)(Generic_ptr data, int line_style);
    void (*set_line_width)(Generic_ptr data, float line_width);
    void (*fill_circle)(Generic_ptr data, float x, float y, float r);
    void (*fill_ellipse)(Generic_ptr data, float x, float y, float rx, float ry);
    void (*draw_circle)(Generic_ptr data, float x, float y, float r);
    void (*draw_ellipse)(Generic_ptr data, float x, float y, float rx, float ry);
    void (*get_background)(Generic_ptr data, float *color);
    void (*get_region)(Generic_ptr data, float *x0, float *y0, float *x1, float *y1);
    void (*get_text_size)(Generic_ptr data, CcpnString text, float *w, float *h);
    void (*fill_triangle)(Generic_ptr data, float x0, float y0,
		    		float x1, float y1, float x2, float y2);
}   Drawing_funcs;

#endif /* _incl_drawing_funcs */

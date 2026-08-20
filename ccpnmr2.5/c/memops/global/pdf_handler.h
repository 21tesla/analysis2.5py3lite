
/*
======================COPYRIGHT/LICENSE START==========================

pdf_handler.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_pdf_handler
#define _incl_pdf_handler

#include "macros.h"
#include "types.h"

#include "drawing_funcs.h"

typedef struct Pdf_handler
{
//    FILE *fp;
    Long_line fp;
    float width, height;
    float x0, y0, x1, y1;
    float ax, ay, bx, by;
    int fontsize;
    Line font;
    Bool isColor;
}   *Pdf_handler;

extern Drawing_funcs *pdf_drawing_funcs(void);

//extern Pdf_handler new_pdf_handler(FILE *fp, float width, float height,
//                                                CcpnString output_style);
extern Pdf_handler new_pdf_handler(CcpnString *file_name, float width, float height,
                                                CcpnString output_style);

extern void delete_pdf_handler(Pdf_handler pdf_handler);

extern void new_range_pdf_handler(Pdf_handler pdf_handler,
		float x0, float y0, float x1, float y1);

extern void clip_range_pdf_handler(Pdf_handler pdf_handler,
		float x0, float y0, float x1, float y1);

extern void draw_line_pdf_handler(Pdf_handler pdf_handler,
				float x0, float y0, float x1, float y1);

extern void draw_clipped_line_pdf_handler(Pdf_handler pdf_handler,
				float x0, float y0, float x1, float y1);

extern void fill_circle_pdf_handler(Pdf_handler pdf_handler,
                                    float x, float y, float r);

extern void fill_ellipse_pdf_handler(Pdf_handler pdf_handler,
                            float x, float y, float rx, float ry);

extern void draw_circle_pdf_handler(Pdf_handler pdf_handler,
                                    float x, float y, float r);

extern void draw_ellipse_pdf_handler(Pdf_handler pdf_handler,
                            float x, float y, float rx, float ry);

extern void draw_polyline_pdf_handler(Pdf_handler pdf_handler,
						Poly_line polyline);

extern void draw_clipped_polyline_pdf_handler(Pdf_handler pdf_handler,
						Poly_line polyline);

extern void draw_dash_line_pdf_handler(Pdf_handler pdf_handler,
				float x0, float y0, float x1, float y1,
				int dash_length, int gap_length);

extern void draw_clipped_dash_line_pdf_handler(Pdf_handler pdf_handler,
				float x0, float y0, float x1, float y1,
				int dash_length, int gap_length);

extern void draw_text_pdf_handler(Pdf_handler pdf_handler, CcpnString text,
				float x, float y, float a, float b);

extern void draw_dash_box_pdf_handler(Pdf_handler pdf_handler,
				float x0, float y0, float x1, float y1);

extern void set_color_pdf_handler(Pdf_handler pdf_handler, float *color);

extern void set_black_pdf_handler(Pdf_handler pdf_handler);

extern void set_font_pdf_handler(Pdf_handler pdf_handler, CcpnString name, int size);

extern void set_line_style_pdf_handler(Pdf_handler pdf_handler, int line_style);

extern void set_line_width_pdf_handler(Pdf_handler pdf_handler, float line_width);

extern void get_region_pdf_handler(Pdf_handler pdf_handler, float *x0, float *y0,
                                                        float *x1, float *y1);

extern void get_text_size_pdf_handler(Pdf_handler pdf_handler, CcpnString text,
                                                        float *w, float *h);

extern void fill_triangle_pdf_handler(Pdf_handler pdf_handler, float x0, float y0,
                                float x1, float y1, float x2, float y2);

#endif /* _incl_pdf_handler */

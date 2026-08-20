
/*
======================COPYRIGHT/LICENSE START==========================

tk_handler.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_tk_handler
#define _incl_tk_handler

#include <tk.h>

#include "macros.h"
#include "types.h"

#include "drawing_funcs.h"

typedef void* Tk_handler;

extern Drawing_funcs *tk_drawing_funcs(void);

extern Tk_handler new_tk_handler(Tcl_Interp *interp, Tk_Window tk_win);

extern void delete_tk_handler(Tk_handler tk_handler);

extern void resize_tk_handler(Tk_handler tk_handler, int width, int height);

extern void expose_tk_handler(Tk_handler tk_handler, int x, int y, int w, int h);

extern void flush_tk_handler(Tk_handler tk_handler);

/* map_range_tk_handler sets up GL so that spectrum point (a0, b0)
   gets mapped to window point (x0, y0), where latter in range (0, 1),
   and similarly for (a1, b1) and (x1, y1) */
extern void map_ranges_tk_handler(Tk_handler tk_handler,
	float x0, float y0, float x1, float y1,
	float a0, float b0, float a1, float b1);

extern Bool make_current_tk_handler(Tk_handler tk_handler);

extern void swap_buffers_tk_handler(Tk_handler tk_handler);

extern void draw_box_tk_handler(Tk_handler tk_handler,
			float x0, float y0, float x1, float y1);

extern void draw_dash_box_tk_handler(Tk_handler tk_handler,
			float x0, float y0, float x1, float y1);

extern void draw_xor_box_tk_handler(Tk_handler tk_handler,
			float x0, float y0, float x1, float y1);

extern void start_xor_tk_handler(Tk_handler tk_handler);

extern void finish_xor_tk_handler(Tk_handler tk_handler);

extern void reset_xor_tk_handler(Tk_handler tk_handler);

extern void clear_xor_tk_handler(Tk_handler tk_handler);

extern void start_front_tk_handler(Tk_handler tk_handler);

extern void start_back_tk_handler(Tk_handler tk_handler);

/* tk_handler takes copy of background */
extern void set_background_tk_handler(Tk_handler tk_handler, float *background);

extern void draw_text_tk_handler(Tk_handler tk_handler, CcpnString text,
				float x, float y, float a, float b);

extern void draw_line_tk_handler(Tk_handler tk_handler,
				float x0, float y0, float x1, float y1);

extern void draw_clipped_line_tk_handler(Tk_handler tk_handler,
				float x0, float y0, float x1, float y1);

extern void fill_circle_tk_handler(Tk_handler tk_handler,
                                                float x, float y, float r);

extern void fill_ellipse_tk_handler(Tk_handler tk_handler,
                                        float x, float y, float rx, float ry);

extern void draw_circle_tk_handler(Tk_handler tk_handler,
                                                float x, float y, float r);

extern void draw_ellipse_tk_handler(Tk_handler tk_handler,
                                        float x, float y, float rx, float ry);

extern void draw_polyline_tk_handler(Tk_handler tk_handler,
						Poly_line polyline);

extern void draw_clipped_polyline_tk_handler(Tk_handler tk_handler,
						Poly_line polyline);

extern void draw_dash_line_tk_handler(Tk_handler tk_handler,
				float x0, float y0, float x1, float y1,
				int dash_length, int gap_length);

extern void draw_clipped_dash_line_tk_handler(Tk_handler tk_handler,
				float x0, float y0, float x1, float y1,
				int dash_length, int gap_length);

extern void set_color_tk_handler(Tk_handler tk_handler, float *color);

extern void set_black_tk_handler(Tk_handler tk_handler);

extern void set_white_tk_handler(Tk_handler tk_handler);

extern void set_line_width_tk_handler(Tk_handler tk_handler, float line_width);

extern void reset_line_width_tk_handler(Tk_handler tk_handler);

extern CcpnStatus set_font_tk_handler(Tk_handler tk_handler, CcpnString name, int size);

extern void get_region_tk_handler(Tk_handler tk_handler, float *x0, float *y0,
                                                        float *x1, float *y1);

extern void get_text_size_tk_handler(Tk_handler tk_handler, CcpnString text,
							float *w, float *h);

extern void fill_triangle_tk_handler(Tk_handler tk_handler, float x0, float y0,
                                float x1, float y1, float x2, float y2);

extern void set_is_double_buffer_tk_handler(Tk_handler tk_handler, Bool is_double_buffer);

#endif /* _incl_tk_handler */


/*
======================COPYRIGHT/LICENSE START==========================

gl_handler.h: Part of the CcpNmr Analysis program

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
#ifdef WIN64
#include <Windows.h>
#endif

#ifndef _incl_gl_handler
#define _incl_gl_handler

#include <tk.h>

#include "macros.h"
#include "types.h"

#include "drawing_funcs.h"

typedef void* Gl_handler;

#ifndef IGNORE_GL

extern Drawing_funcs *gl_drawing_funcs(void);

extern Gl_handler new_gl_handler(Tcl_Interp *interp, Tk_Window tk_win, Bool direct, Bool dbl_buff);

extern void delete_gl_handler(Gl_handler gl_handler);

extern void resize_gl_handler(Gl_handler gl_handler, int width, int height);

extern void expose_gl_handler(Gl_handler gl_handler, int x, int y, int w, int h);

extern void flush_gl_handler(Gl_handler gl_handler);

/* map_range_gl_handler sets up GL so that spectrum point (a0, b0)
   gets mapped to window point (x0, y0), where latter in range (0, 1),
   and similarly for (a1, b1) and (x1, y1) */
extern void map_ranges_gl_handler(Gl_handler gl_handler,
	float x0, float y0, float x1, float y1,
	float a0, float b0, float a1, float b1);

extern Bool make_current_gl_handler(Gl_handler gl_handler);

extern void swap_buffers_gl_handler(Gl_handler gl_handler);

extern void test_gl_handler(Gl_handler gl_handler);

extern void draw_box_gl_handler(Gl_handler gl_handler,
			float x0, float y0, float x1, float y1);

extern void draw_dash_box_gl_handler(Gl_handler gl_handler,
			float x0, float y0, float x1, float y1);

extern void draw_xor_box_gl_handler(Gl_handler gl_handler,
			float x0, float y0, float x1, float y1);

extern void start_xor_gl_handler(Gl_handler gl_handler);

extern void finish_xor_gl_handler(Gl_handler gl_handler);

extern void reset_xor_gl_handler(Gl_handler gl_handler);

extern void clear_xor_gl_handler(Gl_handler gl_handler);

extern void start_front_gl_handler(Gl_handler gl_handler);

extern void start_back_gl_handler(Gl_handler gl_handler);

extern void end_back_gl_handler(Gl_handler gl_handler);

/* gl_handler takes copy of background */
extern void set_background_gl_handler(Gl_handler gl_handler, float *background);

extern void draw_text_gl_handler(Gl_handler gl_handler, CcpnString text,
				float x, float y, float a, float b);

extern void draw_line_gl_handler(Gl_handler gl_handler,
				float x0, float y0, float x1, float y1);

extern void draw_clipped_line_gl_handler(Gl_handler gl_handler,
				float x0, float y0, float x1, float y1);

extern void fill_circle_gl_handler(Gl_handler gl_handler,
                                float x, float y, float r);

extern void fill_ellipse_gl_handler(Gl_handler gl_handler,
                                float x, float y, float rx, float ry);

extern void draw_circle_gl_handler(Gl_handler gl_handler,
                                float x, float y, float r);

extern void draw_ellipse_gl_handler(Gl_handler gl_handler,
                                float x, float y, float rx, float ry);

extern void draw_polyline_gl_handler(Gl_handler gl_handler,
						Poly_line polyline);

extern void draw_clipped_polyline_gl_handler(Gl_handler gl_handler,
						Poly_line polyline);

extern void draw_dash_line_gl_handler(Gl_handler gl_handler,
				float x0, float y0, float x1, float y1,
				int dash_length, int gap_length);

extern void draw_clipped_dash_line_gl_handler(Gl_handler gl_handler,
				float x0, float y0, float x1, float y1,
				int dash_length, int gap_length);

extern void set_color_gl_handler(Gl_handler gl_handler, float *color);

extern void set_black_gl_handler(Gl_handler gl_handler);

extern void set_white_gl_handler(Gl_handler gl_handler);

extern void set_line_width_gl_handler(Gl_handler gl_handler, float line_width);

extern void reset_line_width_gl_handler(Gl_handler gl_handler);

extern CcpnStatus has_font_gl_handler(Gl_handler gl_handler, CcpnString name, int *size);

extern CcpnStatus set_font_gl_handler(Gl_handler gl_handler, CcpnString name, int size);

extern void get_region_gl_handler(Gl_handler gl_handler, float *x0, float *y0,
                                                        float *x1, float *y1);
extern void get_text_size_gl_handler(Gl_handler gl_handler, CcpnString text,
							float *w, float *h);

extern void fill_triangle_gl_handler(Gl_handler gl_handler, float x0, float y0,
				float x1, float y1, float x2, float y2);

extern void set_is_double_buffer_gl_handler(Gl_handler gl_handler, Bool is_double_buffer);

#endif /* IGNORE_GL */

#endif /* _incl_gl_handler */

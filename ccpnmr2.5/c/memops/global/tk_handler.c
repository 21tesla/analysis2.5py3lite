
/*
======================COPYRIGHT/LICENSE START==========================

tk_handler.c: Part of the CcpNmr Analysis program

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
#include "tk_handler.h"

#include "clipping.h"
#include "utility.h"

#define  NCOLORS  3

typedef struct Tk_handler_p
{
    Tk_Window tk_win;
    Tcl_Interp *interp;
    Display *display;
    GC gc;
    Pixmap pixmap;
    Drawable drawable;
    int width;
    int height;
    float background[NCOLORS];
    float x0, y0, x1, y1;
    float sx, sy, tx, ty;
    int line_style;
    CcpnString font_name;
    int font_size;
    Tk_Font font;
    float color[NCOLORS];
    XColor *xcolor;
    Bool is_double_buffer;
}   *Tk_handler_p;

#define  SCALE_X(x)  ((int) (tk_handler_p->width * (tk_handler_p->sx*(x))))
#define  SCALE_Y(y)  ((int) (tk_handler_p->height * (tk_handler_p->sy*(y))))
#define  CONVERT_X(x)  ((int) (tk_handler_p->width * (tk_handler_p->sx*(x) + tk_handler_p->tx)))
#define  CONVERT_Y(y)  (tk_handler_p->height - 1 - (int) (tk_handler_p->height * (tk_handler_p->sy*(y) + tk_handler_p->ty)))

#define  DEFAULT_LINEWIDTH  0.0

static char default_font_name[] = "Helvetica";

static int default_font_size = 10;

static void tk_start_draw(Generic_ptr data)
{
/*
    Tk_handler tk_handler = (Tk_handler) data;

    make_current_tk_handler(tk_handler);
*/
}

static void tk_end_draw(Generic_ptr data)
{
/* TODO */
}

static void tk_new_draw_range(Generic_ptr data, float x0, float y0,
					float x1, float y1, Bool clip)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) data;

    tk_handler_p->x0 = MIN(x0, x1);
    tk_handler_p->y0 = MIN(y0, y1);
    tk_handler_p->x1 = MAX(x0, x1);
    tk_handler_p->y1 = MAX(y0, y1);
}

static void tk_draw_line(Generic_ptr data, float x0, float y0, float x1, float y1)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) data;
    static int dash_length = 2, gap_length = 2;

    if ((tk_handler_p != NULL) && (tk_handler_p->line_style == DASHED_LINE_STYLE))
	draw_dash_line_tk_handler((Tk_handler) data, x0, y0, x1, y1,
					dash_length, gap_length);
    else
	draw_line_tk_handler((Tk_handler) data, x0, y0, x1, y1);
}

static void tk_draw_clipped_line(Generic_ptr data, float x0, float y0, float x1, float y1)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) data;
    static int dash_length = 2, gap_length = 2;

    if ((tk_handler_p != NULL) && (tk_handler_p->line_style == DASHED_LINE_STYLE))
	draw_clipped_dash_line_tk_handler((Tk_handler) data, x0, y0, x1, y1,
					dash_length, gap_length);
    else
	draw_clipped_line_tk_handler((Tk_handler) data, x0, y0, x1, y1);
}

static void tk_draw_polyline(Generic_ptr data, Poly_line polyline)
{
    draw_polyline_tk_handler((Tk_handler) data, polyline);
}

static void tk_draw_clipped_polyline(Generic_ptr data, Poly_line polyline)
{
    draw_clipped_polyline_tk_handler((Tk_handler) data, polyline);
}

static void tk_draw_text(Generic_ptr data, CcpnString text, float x, float y,
							float a, float b)
{
    draw_text_tk_handler((Tk_handler) data, text, x, y, a, b);
}

static void tk_set_draw_color(Generic_ptr data, float *color)
{
    set_color_tk_handler((Tk_handler) data, color);
}

static void tk_set_draw_font(Generic_ptr data, CcpnString name, int size)
{
    set_font_tk_handler((Tk_handler) data, name, size);
}

static void tk_set_line_style(Generic_ptr data, int line_style)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) data;

    tk_handler_p->line_style = line_style;
}

static void tk_set_line_width(Generic_ptr data, float line_width)
{
    set_line_width_tk_handler((Tk_handler) data, line_width);
}

static void tk_fill_circle(Generic_ptr data, float x, float y, float r)
{
    fill_circle_tk_handler((Tk_handler) data, x, y, r);
}

static void tk_fill_ellipse(Generic_ptr data, float x, float y, float rx, float ry)
{
    fill_ellipse_tk_handler((Tk_handler) data, x, y, rx, ry);
}

static void tk_draw_circle(Generic_ptr data, float x, float y, float r)
{
    draw_circle_tk_handler((Tk_handler) data, x, y, r);
}

static void tk_draw_ellipse(Generic_ptr data, float x, float y, float rx, float ry)
{
    draw_ellipse_tk_handler((Tk_handler) data, x, y, rx, ry);
}

static void tk_get_background(Generic_ptr data, float *color)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) data;

    COPY_VECTOR(color, tk_handler_p->background, NCOLORS);
}

static void tk_get_region(Generic_ptr data, float *x0, float *y0,
                                                    float *x1, float *y1)
{
    get_region_tk_handler((Tk_handler) data, x0, y0, x1, y1);
}

static void tk_get_text_size(Generic_ptr data, CcpnString text, float *w, float *h)
{
    get_text_size_tk_handler((Tk_handler) data, text, w, h);
}

static void tk_fill_triangle(Generic_ptr data, float x0, float y0,
				float x1, float y1, float x2, float y2)
{
    fill_triangle_tk_handler((Tk_handler) data, x0, y0, x1, y1, x2, y2);
}

static Drawing_funcs drawing_funcs = { SCREEN_DRAWING,
	tk_start_draw, tk_end_draw, tk_new_draw_range,
	tk_draw_line, tk_draw_clipped_line,
	tk_draw_polyline, tk_draw_clipped_polyline, tk_draw_text,
	tk_set_draw_color, tk_set_draw_font, tk_set_line_style,
        tk_set_line_width, tk_fill_circle, tk_fill_ellipse,
	tk_draw_circle, tk_draw_ellipse, tk_get_background,
	tk_get_region, tk_get_text_size, tk_fill_triangle };

Drawing_funcs *tk_drawing_funcs()
{
    return &drawing_funcs;
}

static Bool init_tk_handler(Tk_handler tk_handler)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    int i;

/*
    if (!make_current_tk_handler(tk_handler))
      return CCPN_FALSE;
*/

    for (i = 0; i < NCOLORS; i++)
	tk_handler_p->background[i] = 1;  /* default = white */

    /* arbitrary */
    tk_handler_p->width = 1;
    tk_handler_p->height = 1;
    tk_handler_p->sx = 1.0;
    tk_handler_p->sy = 1.0;
    tk_handler_p->tx = 0.0;
    tk_handler_p->ty = 0.0;
    tk_handler_p->line_style = NORMAL_LINE_STYLE;

    return CCPN_TRUE;
}

static void end_tk_handler(Tk_handler tk_handler)
{
/*
    printf("end_tk_handler\n");
*/
}

Tk_handler new_tk_handler(Tcl_Interp *interp, Tk_Window tk_win)
{
    int i;
    Tk_handler tk_handler;
    Tk_handler_p tk_handler_p;
    Display *display;
    GC gc;
    XGCValues gcv;
    static char stipple_data[] = {
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22,
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22
    };
    
    Tk_MakeWindowExist(tk_win);
    display = Tk_Display(tk_win);

    Tk_DefineBitmap(interp, Tk_GetUid("stipple_data"), (char *) stipple_data, 32, 32);
    gcv.stipple = Tk_GetBitmap(interp, tk_win, Tk_GetUid("stipple_data"));

    gc = XCreateGC(display, Tk_WindowId(tk_win), GCForeground | GCBackground | GCStipple, &gcv);

    if (!gc)
        return NULL;

    MALLOC_NEW(tk_handler_p, struct Tk_handler_p, 1);

    tk_handler_p->tk_win = tk_win;
    tk_handler_p->interp = interp;
    tk_handler_p->display = display;
    tk_handler_p->pixmap = (Pixmap) 0;
    tk_handler_p->drawable = Tk_WindowId(tk_win);

    tk_handler_p->gc = gc;
    tk_handler_p->font_name = NULL;
    tk_handler_p->font_size = 10;
    tk_handler_p->font = NULL;
    for (i = 0; i < NCOLORS; i++)
        tk_handler_p->color[i] = -1;
    tk_handler_p->xcolor = NULL;
    tk_handler_p->is_double_buffer = CCPN_TRUE;

    tk_handler = (Tk_handler) tk_handler_p;

    if (!init_tk_handler(tk_handler))
    {
        delete_tk_handler(tk_handler);

        return NULL;
    }

/*
    printf("new_tk_handler %x\n", tk_handler);
*/

    return tk_handler;
}

void delete_tk_handler(Tk_handler tk_handler)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

/*
    printf("delete_tk_handler %x\n", tk_handler);
*/

    if (tk_handler)
    {
	end_tk_handler(tk_handler);

	XFreeGC(tk_handler_p->display, tk_handler_p->gc);

        if (tk_handler_p->pixmap)
        {
            Tk_FreePixmap(tk_handler_p->display, tk_handler_p->pixmap);
	}

        if (tk_handler_p->xcolor)
        {
            Tk_FreeColor(tk_handler_p->xcolor);
	}

        if (tk_handler_p->font)
        {
            Tk_FreeFont(tk_handler_p->font);
            tk_handler_p->font = NULL;
        }

	FREE(tk_handler_p->font_name, char);

	FREE(tk_handler_p, struct Tk_handler_p);
    }
}

void resize_tk_handler(Tk_handler tk_handler, int width, int height)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    Window win = Tk_WindowId(tk_handler_p->tk_win);

/*
    if (!make_current_tk_handler(tk_handler))
      return;

    XResizeWindow(tk_handler->display, Tk_WindowId(tk_handler->tk_win),
							width, height);
    printf("resize_tk_handler1: width = %d, height = %d\n", width, height);
    printf("resize_tk_handler2: width = %d, height = %d\n", Tk_Width(tk_handler_p->tk_win), Tk_Height(tk_handler_p->tk_win));
*/

    if (tk_handler_p->pixmap && (tk_handler_p->width == width)
                                && (tk_handler_p->height == height))
        return;

    tk_handler_p->width = width;
    tk_handler_p->height = height;

    if (tk_handler_p->pixmap)
        Tk_FreePixmap(tk_handler_p->display, tk_handler_p->pixmap);

    if (width && height)
        tk_handler_p->pixmap = Tk_GetPixmap(tk_handler_p->display, win, width, height,
		                            Tk_Depth(tk_handler_p->tk_win));
/* below crashes on OSX native (i.e. Aqua)
                            DefaultDepthOfScreen(Tk_Screen(tk_handler_p->tk_win)));
*/

/*
    printf("resize_tk_handler: pixmap = %d\n", tk_handler_p->pixmap);
*/
}

void expose_tk_handler(Tk_handler tk_handler, int x, int y, int w, int h)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

/*
    printf("expose_tk_handler: x = %d, y = %d, w = %d, h = %d\n", x, y, w, h);
*/
    if (!make_current_tk_handler(tk_handler))
      return;

    set_color_tk_handler(tk_handler, tk_handler_p->background);
    XFillRectangle(tk_handler_p->display, tk_handler_p->drawable,
                                tk_handler_p->gc, x, y, w, h);
}

void flush_tk_handler(Tk_handler tk_handler)
{
}

void map_ranges_tk_handler(Tk_handler tk_handler,
		float x0, float y0, float x1, float y1,
		float a0, float b0, float a1, float b1)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    float sx, sy;

    if ((x0 == x1) || (y0 == y1) || (a0 == a1) || (b0 == b1))
	return;

    if (a0 > a1)
    {
	SWAP(a0, a1, float);
	SWAP(x0, x1, float);
    }

    if (b0 > b1)
    {
	SWAP(b0, b1, float);
	SWAP(y0, y1, float);
    }
/*
    if (!make_current_tk_handler(tk_handler))
      return;
*/

    sx = tk_handler_p->sx = (x1 - x0) / (a1 - a0);
    sy = tk_handler_p->sy = (y1 - y0) / (b1 - b0);
    tk_handler_p->tx = (x0*a1-x1*a0) / (a1 - a0);
    tk_handler_p->ty = (y0*b1-y1*b0) / (b1 - b0);

    a0 -= x0 / sx;
    a1 += (1-x1) / sx;
    b0 -= y0 / sy;
    b1 += (1-y1) / sy;
    tk_new_draw_range(tk_handler, a0, b0, a1, b1, CCPN_TRUE);
/*
    printf("map_ranges_tk_handler1: %2.1f, %2.1f, %2.1f, %2.1f\n", x0, y0, x1, y1);
    printf("map_ranges_tk_handler2: %2.1f, %2.1f, %2.1f, %2.1f\n", a0, b0, a1, b1);
    printf("map_ranges_tk_handler3: %f, %f, %f, %f\n", sx, sy, tk_handler_p->tx, tk_handler_p->ty);
*/
}

Bool make_current_tk_handler(Tk_handler tk_handler)
{
/*
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
*/

    return CCPN_TRUE;
}

void swap_buffers_tk_handler(Tk_handler tk_handler)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    Window win = Tk_WindowId(tk_handler_p->tk_win);

/*
    printf("swap_buffers_tk_handler: %d %d %d\n", win, tk_handler_p->drawable, tk_handler_p->pixmap);
    if (!make_current_tk_handler(tk_handler))
      return;
*/
    if (tk_handler_p->drawable == tk_handler_p->pixmap)
        XCopyArea(tk_handler_p->display, tk_handler_p->pixmap, win, tk_handler_p->gc,
                        0, 0, tk_handler_p->width, tk_handler_p->height, 0, 0);
}

void draw_box_tk_handler(Tk_handler tk_handler,
			float x0, float y0, float x1, float y1)
{
/*
    printf("draw_box_tk_handler: %4.3f, %4.3f, %4.3f, %4.3f\n", x0, y0, x1, y1);
*/
    draw_clipped_line_tk_handler(tk_handler, x0, y0, x1, y0);
    draw_clipped_line_tk_handler(tk_handler, x1, y0, x1, y1);
    draw_clipped_line_tk_handler(tk_handler, x1, y1, x0, y1);
    draw_clipped_line_tk_handler(tk_handler, x0, y1, x0, y0);
}

void draw_dash_box_tk_handler(Tk_handler tk_handler,
			float x0, float y0, float x1, float y1)
{
/*
    printf("draw_dash_box_tk_handler: %4.3f, %4.3f, %4.3f, %4.3f\n", x0, y0, x1, y1);
*/

    draw_clipped_dash_line_tk_handler(tk_handler, x0, y0, x1, y0, 2, 6);
    draw_clipped_dash_line_tk_handler(tk_handler, x1, y0, x1, y1, 2, 6);
    draw_clipped_dash_line_tk_handler(tk_handler, x1, y1, x0, y1, 2, 6);
    draw_clipped_dash_line_tk_handler(tk_handler, x0, y1, x0, y0, 2, 6);
}

void draw_xor_box_tk_handler(Tk_handler tk_handler,
			float x0, float y0, float x1, float y1)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    int xx0 = CONVERT_X(x0);
    int yy0 = CONVERT_Y(y0);
    int xx1 = CONVERT_X(x1);
    int yy1 = CONVERT_Y(y1);
    int w, h;
    XGCValues gcv;

    if (xx0 > xx1)
	SWAP(xx0, xx1, int);

    if (yy0 > yy1)
	SWAP(yy0, yy1, int);

    w = xx1 - xx0;
    h = yy1 - yy0;

    if (!(w && h))
        return;
/*
    printf("draw_xor_box_tk_handler: %4.3f, %4.3f, %4.3f, %4.3f\n", x0, y0, x1, y1);
*/
    start_xor_tk_handler(tk_handler);
/*
    draw_box_tk_handler(tk_handler, x0, y0, x1, y1);
*/

    // NOTE:ED - change to draw a series of dotted lines to simulate stipple - doesn't work in MacOS
    int dash_length = 1;
    int gap_length = 3;
    int dash_offset = (gap_length+dash_length)/2;
    int ndashes = 2;
    char dash_list[2];

    dash_list[0] = (char) dash_length;
    dash_list[1] = (char) gap_length;

    gcv.line_style = LineOnOffDash;
    XChangeGC(tk_handler_p->display, tk_handler_p->gc, GCLineStyle, &gcv);

    // draw alternating lines to give offset to stippling
    XSetDashes(tk_handler_p->display, tk_handler_p->gc, 0, dash_list, sizeof(dash_list));
    int x;
    for (x=xx0+1; x<xx1; x+=(gap_length+dash_length))
        XDrawLine(tk_handler_p->display, tk_handler_p->drawable,
                                tk_handler_p->gc, x, yy0+1, x, yy1-1);

    // change the start offset for the stipple
    XSetDashes(tk_handler_p->display, tk_handler_p->gc, dash_offset, dash_list, sizeof(dash_list));
    for (x=xx0+dash_offset+1; x<xx1; x+=(gap_length+dash_length))
        XDrawLine(tk_handler_p->display, tk_handler_p->drawable,
                                tk_handler_p->gc, x, yy0+1, x, yy1-1);

    // draw a solid box around the selection
    gcv.line_style = LineSolid;
    XChangeGC(tk_handler_p->display, tk_handler_p->gc, GCLineStyle, &gcv);
    XDrawRectangle(tk_handler_p->display, tk_handler_p->drawable, tk_handler_p->gc, xx0, yy0, w, h);

    finish_xor_tk_handler(tk_handler);
}

void start_xor_tk_handler(Tk_handler tk_handler)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
/*
    printf("start_xor_tk_handler\n");
    clear_xor_tk_handler(tk_handler);
*/
    swap_buffers_tk_handler(tk_handler);
/* below is because Solaris seems to require equiv to do xor (not sure why) */
#ifdef XOR_IS_EQUIV
    XSetFunction(tk_handler_p->display, tk_handler_p->gc, GXequiv);
#else
    XSetFunction(tk_handler_p->display, tk_handler_p->gc, GXxor);
#endif
    start_front_tk_handler(tk_handler);
}

void finish_xor_tk_handler(Tk_handler tk_handler)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

    XSetFunction(tk_handler_p->display, tk_handler_p->gc, GXcopy);
/*
    swap_buffers_tk_handler(tk_handler);
*/
/* TBD: should only set to start_back if was this to start with */
/* but whole xor'ing procedure implemented here doesn't work unless using back to start with */
    if (tk_handler_p->is_double_buffer)
	start_back_tk_handler(tk_handler);
}

void reset_xor_tk_handler(Tk_handler tk_handler)
{
/*
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

    printf("reset_xor_tk_handler\n");
*/
}

void clear_xor_tk_handler(Tk_handler tk_handler)
{
/*
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

    printf("clear_xor_tk_handler\n");
*/
    start_xor_tk_handler(tk_handler);
    finish_xor_tk_handler(tk_handler);
}

void start_front_tk_handler(Tk_handler tk_handler)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    Window win = Tk_WindowId(tk_handler_p->tk_win);

/*
    printf("start_front_tk_handler\n");
*/
    tk_handler_p->drawable = win;
}

void start_back_tk_handler(Tk_handler tk_handler)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

/*
    printf("start_back_tk_handler\n");
*/
    if (tk_handler_p->pixmap)
        tk_handler_p->drawable = (Drawable) tk_handler_p->pixmap;
    else
        start_front_tk_handler(tk_handler);
}

void end_back_tk_handler(Tk_handler tk_handler)
{
}

void set_background_tk_handler(Tk_handler tk_handler, float *background)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

    COPY_VECTOR(tk_handler_p->background, background, NCOLORS);
}

void draw_text_tk_handler(Tk_handler tk_handler, CcpnString text, float x, float y,
							float a, float b)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    int s = CONVERT_X(x);
    int t = CONVERT_Y(y);
    float w, h;
    Tk_Font font = tk_handler_p->font;

/*
    printf("draw_text_tk_handler: '%s' %s %d %x\n", text, tk_handler_p->font_name, tk_handler_p->font_size, font);
*/

    if (!font)
    {
        if (set_font_tk_handler(tk_handler, default_font_name, default_font_size) == CCPN_ERROR)
            return; /* not much else can do */
    }

    font = tk_handler_p->font;
    w = Tk_TextWidth(font, text, strlen(text));

    get_text_size_tk_handler(tk_handler, text, &w, &h);
    s -= a * w;
    t += b * h;
    Tk_DrawChars(tk_handler_p->display, tk_handler_p->drawable,
                        tk_handler_p->gc, font, text, strlen(text), s, t);

/*
    printf("draw_text_tk_handler2\n");
*/

}

void draw_line_tk_handler(Tk_handler tk_handler,
			float x0, float y0, float x1, float y1)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
/*
    Window win = Tk_WindowId(tk_handler_p->tk_win);
*/
    int s0 = CONVERT_X(x0);
    int t0 = CONVERT_Y(y0);
    int s1 = CONVERT_X(x1);
    int t1 = CONVERT_Y(y1);

/*
    printf("draw_line_tk_handler1: %2.1f %2.1f %2.1f %2.1f\n", x0, y0, x1, y1);
    printf("draw_line_tk_handler2: %d %d %d %d\n", s0, t0, s1, t1);
    printf("draw_line_tk_handler3: %d\n", tk_handler_p->drawable);
*/

    XDrawLine(tk_handler_p->display, tk_handler_p->drawable,
                                tk_handler_p->gc, s0, t0, s1, t1);
/**
    XDrawLine(tk_handler_p->display, win,
                                tk_handler_p->gc, s0, t0, s1, t1);
**/
}

void draw_clipped_line_tk_handler(Tk_handler tk_handler,
                        float x0, float y0, float x1, float y1)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    float xmin = tk_handler_p->x0;
    float ymin = tk_handler_p->y0;
    float xmax = tk_handler_p->x1;
    float ymax = tk_handler_p->y1;

    draw_clipped_line(x0, y0, x1, y1, &drawing_funcs,
                (Generic_ptr) tk_handler, xmin, ymin, xmax, ymax);
}

void fill_circle_tk_handler(Tk_handler tk_handler, float x, float y, float r)
{
    fill_ellipse_tk_handler(tk_handler, x, y, r, r);
}

void fill_ellipse_tk_handler(Tk_handler tk_handler, float x, float y, float rx, float ry)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    int xx = CONVERT_X(x);
    int yy = CONVERT_Y(y);
    int rxx = SCALE_X(rx);
    int ryy = SCALE_Y(ry);

/*
    printf("draw_circle_tk_handler1: %2.1f %2.1f %2.1f\n", x, y, r);
    printf("draw_circle_tk_handler2: %d %d %d %d\n", xx, yy, rx, ry);
*/

    XFillArc(tk_handler_p->display, tk_handler_p->drawable,
                tk_handler_p->gc, xx-rxx, yy-ryy, 2*rxx, 2*ryy, 64*0, 64*360);
}

void draw_circle_tk_handler(Tk_handler tk_handler, float x, float y, float r)
{
    draw_ellipse_tk_handler(tk_handler, x, y, r, r);
}

void draw_ellipse_tk_handler(Tk_handler tk_handler, float x, float y, float rx, float ry)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    int xx = CONVERT_X(x);
    int yy = CONVERT_Y(y);
    int rxx = SCALE_X(rx);
    int ryy = SCALE_Y(ry);

/*
    printf("draw_circle_tk_handler1: %2.1f %2.1f %2.1f\n", x, y, r);
    printf("draw_circle_tk_handler2: %d %d %d %d\n", xx, yy, rx, ry);
*/

    XDrawArc(tk_handler_p->display, tk_handler_p->drawable,
                tk_handler_p->gc, xx-rxx, yy-ryy, 2*rxx, 2*ryy, 64*0, 64*360);
}

void draw_polyline_tk_handler(Tk_handler tk_handler, Poly_line polyline)
{
    int i, n = polyline->nvertices;
    Point2f *v = polyline->vertices;
    float x0, y0, x1, y1;

/*
printf("tk_draw_polyline: n = %2d\n", n);
*/

    x0 = v[0].x;
    y0 = v[0].y;
    for (i = 1; i < n; i++)
    {
	x1 = v[i].x;
	y1 = v[i].y;
        draw_line_tk_handler(tk_handler, x0, y0, x1, y1);
        x0 = x1;
        y0 = y1;
    }

    if (polyline->closed)
        draw_line_tk_handler(tk_handler, x0, y0, v[0].x, v[0].y);

/*
printf("tk_draw_polyline2\n");
*/
}

void draw_clipped_polyline_tk_handler(Tk_handler tk_handler, Poly_line polyline)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    float xmin = tk_handler_p->x0;
    float ymin = tk_handler_p->y0;
    float xmax = tk_handler_p->x1;
    float ymax = tk_handler_p->y1;

    draw_clipped_polyline(polyline, &drawing_funcs,
                (Generic_ptr) tk_handler, xmin, ymin, xmax, ymax);
}

void draw_dash_line_tk_handler(Tk_handler tk_handler,
			float x0, float y0, float x1, float y1,
			int dash_length, int gap_length)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    int dash_offset = 0, ndashes = 2;
    char dash_list[2];
    XGCValues gcv;

    dash_list[0] = (char) dash_length;
    dash_list[1] = (char) gap_length;

    gcv.line_style = LineOnOffDash;
    XChangeGC(tk_handler_p->display, tk_handler_p->gc, GCLineStyle, &gcv);
    XSetDashes(tk_handler_p->display, tk_handler_p->gc, dash_offset, dash_list, ndashes);
    draw_line_tk_handler(tk_handler, x0, y0, x1, y1);
    gcv.line_style = LineSolid;
    XChangeGC(tk_handler_p->display, tk_handler_p->gc, GCLineStyle, &gcv);
}

void draw_clipped_dash_line_tk_handler(Tk_handler tk_handler,
			float x0, float y0, float x1, float y1,
			int dash_length, int gap_length)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    int dash_offset = 0, ndashes = 2;
    char dash_list[2];
    XGCValues gcv;

    dash_list[0] = (char) dash_length;
    dash_list[1] = (char) gap_length;

    gcv.line_style = LineOnOffDash;
    XChangeGC(tk_handler_p->display, tk_handler_p->gc, GCLineStyle, &gcv);
    XSetDashes(tk_handler_p->display, tk_handler_p->gc, dash_offset, dash_list, ndashes);
    draw_clipped_line_tk_handler(tk_handler, x0, y0, x1, y1);
    gcv.line_style = LineSolid;
    XChangeGC(tk_handler_p->display, tk_handler_p->gc, GCLineStyle, &gcv);
}

#define  CONVERT_COLOR(t) ((unsigned short) MIN(65535, 65536*(t)))

void set_color_tk_handler(Tk_handler tk_handler, float *color)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    int i;
    XColor c, *cc;

    if (tk_handler_p->xcolor)
    {
        for (i = 0; i < NCOLORS; i++)
        {
            if (color[i] != tk_handler_p->color[i])
                break;
        }

        if (i == NCOLORS)
            return;

/*      Do not do this, it takes ages and thrashes everything
        Just live with the memory leak (minor, it's only colors)
        Tk_FreeColor(tk_handler_p->xcolor);
*/
    }

/*
    printf("set_color_tk_handler0: %2.1f %2.1f, %2.1f\n", color[0], color[1], color[2]);
*/
    COPY_VECTOR(tk_handler_p->color, color, NCOLORS);

    c.red = CONVERT_COLOR(color[0]);
    c.green = CONVERT_COLOR(color[1]);
    c.blue = CONVERT_COLOR(color[2]);

/*
    printf("set_color_tk_handler1: %d %d, %d\n", c.red, c.green, c.blue);
    printf("set_color_tk_handler2: %d 0x%x\n", c.pixel, &c);
*/
    tk_handler_p->xcolor = cc = Tk_GetColorByValue(tk_handler_p->tk_win, &c);
/*
    printf("set_color_tk_handler3: %d 0x%x\n", cc->pixel, cc);
*/

    XSetForeground(tk_handler_p->display, tk_handler_p->gc, cc->pixel);
}

void set_black_tk_handler(Tk_handler tk_handler)
{
    float color[NCOLORS];

    ZERO_VECTOR(color, NCOLORS);

    set_color_tk_handler(tk_handler, color);
}

void set_white_tk_handler(Tk_handler tk_handler)
{
    int i;
    float color[NCOLORS];

    for (i = 0; i < NCOLORS; i++)
        color[i] = 1;

    set_color_tk_handler(tk_handler, color);
}

void set_line_width_tk_handler(Tk_handler tk_handler, float line_width)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    XGCValues gcv;

    gcv.line_width = (unsigned int) line_width;
    XChangeGC(tk_handler_p->display, tk_handler_p->gc, GCLineWidth, &gcv);
}

void reset_line_width_tk_handler(Tk_handler tk_handler)
{
    set_line_width_tk_handler(tk_handler, DEFAULT_LINEWIDTH);
}

CcpnStatus set_font_tk_handler(Tk_handler tk_handler, CcpnString name, int size)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    Tk_Font font;
    Long_line fontstring;

    if ((size == tk_handler_p->font_size) &&
			equal_strings(name, tk_handler_p->font_name))
        return CCPN_OK;

    sprintf(fontstring, "%s -%d", name, size);

    font = Tk_GetFont(tk_handler_p->interp, tk_handler_p->tk_win, fontstring);
    if (!font)
    {
        name = default_font_name;
		size = default_font_size;
    	if ((size == tk_handler_p->font_size) &&
			equal_strings(name, tk_handler_p->font_name))
            return CCPN_OK;

    	sprintf(fontstring, "%s -%d", name, size);
        font = Tk_GetFont(tk_handler_p->interp, tk_handler_p->tk_win, fontstring);
        if (!font)
            return CCPN_ERROR;
    }

    if (tk_handler_p->font)
        Tk_FreeFont(tk_handler_p->font);

    XSetFont(tk_handler_p->display, tk_handler_p->gc, Tk_FontId(font));

    tk_handler_p->font = font;
    tk_handler_p->font_size = size;

    if (!equal_strings(name, tk_handler_p->font_name))
    {
    	FREE(tk_handler_p->font_name, char);
    	STRING_MALLOC_COPY(tk_handler_p->font_name, name);
    }

    return CCPN_OK;
}

void get_region_tk_handler(Tk_handler tk_handler, float *x0, float *y0,
                                                    float *x1, float *y1)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

    *x0 = tk_handler_p->x0;
    *y0 = tk_handler_p->y0;
    *x1 = tk_handler_p->x1;
    *y1 = tk_handler_p->y1;
}

void get_text_size_tk_handler(Tk_handler tk_handler, CcpnString text,
							float *w, float *h)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    Tk_Font font = tk_handler_p->font;
    Tk_FontMetrics metrics;

/*
    printf("get_text_size_tk_handler: '%s'\n", text);
*/
    if (!font)
    {
        if (set_font_tk_handler(tk_handler, default_font_name, default_font_size) == CCPN_ERROR)
	{
	    *w = *h = 10; /* arbitrary*/
            return; /* not much else can do */
	}
    }

    font = tk_handler_p->font;
    Tk_GetFontMetrics(font, &metrics);
    *w = Tk_TextWidth(font, text, strlen(text));
    *h = metrics.linespace;
}

void fill_triangle_tk_handler(Tk_handler tk_handler, float x0, float y0,
                                float x1, float y1, float x2, float y2)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    int npoints = 3;
    XPoint points[3];

    points[0].x = CONVERT_X(x0);
    points[0].y = CONVERT_Y(y0);
    points[1].x = CONVERT_X(x1);
    points[1].y = CONVERT_Y(y1);
    points[2].x = CONVERT_X(x2);
    points[2].y = CONVERT_Y(y2);

    XFillPolygon(tk_handler_p->display, tk_handler_p->drawable,
	tk_handler_p->gc, points, npoints, Convex, CoordModeOrigin);
}

void set_is_double_buffer_tk_handler(Tk_handler tk_handler, Bool is_double_buffer)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

    tk_handler_p->is_double_buffer = is_double_buffer;
}

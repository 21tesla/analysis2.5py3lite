
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
#include <stdio.h>
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
#ifdef __APPLE__
    /* Native (Aqua) Tk: the widget is NOT an X11 window, so the Xlib path
       below cannot reach it.  All drawing is issued to a canvas child
       widget that fills tk_win:
         back buffer  ->  canvas items tagged "ccpBack"
         front / XOR  ->  canvas items tagged "ccpFront"
         swapBuffers  ->  nothing (items only become visible when the Tk
                          event loop runs, i.e. after the Python draw pass
                          finishes - the same visible timing as the pixmap
                          swap, and no ghosting is possible) */
    CcpnString canvas_path;    /* path of the drawing canvas */
    CcpnString font_spec;      /* canvas -font spec, kept in sync with font */
    char canvas_color[8];      /* current draw color "#rrggbb" */
    float cur_line_width;
    int front_layer;           /* 1: items go on the front layer */
    char dash[16];             /* active dash pattern ("d g") */
    Bool dash_set;
#endif
}   *Tk_handler_p;

#define  SCALE_X(x)  ((int) (tk_handler_p->width * (tk_handler_p->sx*(x))))
#define  SCALE_Y(y)  ((int) (tk_handler_p->height * (tk_handler_p->sy*(y))))
#define  CONVERT_X(x)  ((int) (tk_handler_p->width * (tk_handler_p->sx*(x) + tk_handler_p->tx)))
#define  CONVERT_Y(y)  (tk_handler_p->height - 1 - (int) (tk_handler_p->height * (tk_handler_p->sy*(y) + tk_handler_p->ty)))

#define  DEFAULT_LINEWIDTH  0.0

static char default_font_name[] = "Helvetica";

static int default_font_size = 10;

/*
================ macOS (Aqua) canvas drawing backend ================

On macOS the Tk window is a native (Aqua) widget - there is no X11 window
or Display behind it, so the Xlib primitives further down cannot target
it.  The same C drawing API is served here from a single borderless
canvas child widget that fills tk_win.  Items are tagged by layer
("ccpBack" / "ccpFront"); clearing a tag replaces the X pixmap/XOR
erase semantics.
*/
#ifdef __APPLE__

static const char *ccp_tag_back  = "ccpBack";
static const char *ccp_tag_front = "ccpFront";

static void tkc_fill_color(char *out, const float *c)
{
    int i, v[3];

    for (i = 0; i < 3; i++)
    v[i] = (int) (255.0 * (c[i] < 0.0 ? 0.0 : (c[i] > 1.0 ? 1.0 : c[i])) + 0.5);

    sprintf(out, "#%02x%02x%02x", v[0], v[1], v[2]);
}

/*
    Tcl 9 refcount contract: Tcl_New*Obj objects start at refCount 0 - the
    holder's reference is NOT counted (in Tcl 8 it was counted: refCount 1).
    A command that retains an argument INCRs it (0 → 1), and if the caller
    then DECREFs afterwards "once per object" (the Tcl 8 idiom used by every
    tkc_* caller below), the retained object is freed while the command
    still points at it → use-after-free at canvas render time.  INCRef
    each argument here so the caller's hold is counted (0 → 1): retained
    objects survive the caller's DECRef (1 → 0 happens after the command
    has INCRed itself, i.e. 2 → 1), and non-retained ones are freed exactly
    once by the caller.
*/
static void tkc_eval(Tk_handler_p tk_handler_p, int argc, Tcl_Obj **argv)
{
    Tcl_Interp *interp = tk_handler_p->interp;
    int i;

    for (i = 0; i < argc; i++)
        if (argv[i])
            Tcl_IncrRefCount(argv[i]);

    if (getenv("CCP_TK_DEBUG"))
        for (i = 0; i < argc; i++)
            if (argv[i])
                fprintf(stderr, "  arg[%d]=%p '%s'\n", i, (void *) argv[i],
                        Tcl_GetString(argv[i]));
            else
                fprintf(stderr, "  arg[%d]=<NULL>\n", i);

    if (Tcl_EvalObjv(interp, argc, argv, TCL_EVAL_GLOBAL) != TCL_OK)
    {
        char *msg = Tcl_GetStringResult(interp);

        fprintf(stderr, "TkHandler canvas error: %s\n", msg ? msg : "(no message)");
        /* Tcl 9: the failed eval's result is owned by the interpreter and
           freed on the next eval - nothing to reset */
    }
}

static void tkc_item_options(Tk_handler_p tk_handler_p, Tcl_Obj **argv, int *n,
                             int filled, const char *color)
{
    int lw;

    lw = (int) (tk_handler_p->cur_line_width + 0.5);
    if (lw < 1)
        lw = 1;

    if (filled)
    {
        /* outline = fill: the X11 fills (XFillArc/XFillRectangle/
           XFillPolygon) draw NO border, and an empty-string color is
           not a valid canvas color ("unknown color name") */
        argv[(*n)++] = Tcl_NewStringObj("-fill", -1);
        argv[(*n)++] = Tcl_NewStringObj(color, -1);
        argv[(*n)++] = Tcl_NewStringObj("-outline", -1);
        argv[(*n)++] = Tcl_NewStringObj(color, -1);
    }
    else
    {
        argv[(*n)++] = Tcl_NewStringObj("-outline", -1);
        argv[(*n)++] = Tcl_NewStringObj(color, -1);
        argv[(*n)++] = Tcl_NewStringObj("-width", -1);
        argv[(*n)++] = Tcl_NewIntObj(lw);
    }
    argv[(*n)++] = Tcl_NewStringObj("-tags", -1);
    argv[(*n)++] = Tcl_NewStringObj(tk_handler_p->front_layer ? ccp_tag_front
                                                              : ccp_tag_back, -1);
}

static void tkc_delete_layer(Tk_handler_p tk_handler_p, const char *tag)
{
    Tcl_Obj *argv[3];

    argv[0] = Tcl_NewStringObj(tk_handler_p->canvas_path, -1);
    argv[1] = Tcl_NewStringObj("delete", -1);
    argv[2] = Tcl_NewStringObj(tag, -1);
    tkc_eval(tk_handler_p, 3, argv);
    Tcl_DecrRefCount(argv[0]);
    Tcl_DecrRefCount(argv[1]);
    Tcl_DecrRefCount(argv[2]);
}

static void tkc_line(Tk_handler_p tk_handler_p, int x0, int y0, int x1, int y1)
{
    Tcl_Obj *argv[18] = {0};
    int n = 0;
    int lw;

    lw = (int) (tk_handler_p->cur_line_width + 0.5);
    if (lw < 1)
        lw = 1;

    argv[n++] = Tcl_NewStringObj(tk_handler_p->canvas_path, -1);
    argv[n++] = Tcl_NewStringObj("create", -1);
    argv[n++] = Tcl_NewStringObj("line", -1);
    argv[n++] = Tcl_NewIntObj(x0);
    argv[n++] = Tcl_NewIntObj(y0);
    argv[n++] = Tcl_NewIntObj(x1);
    argv[n++] = Tcl_NewIntObj(y1);
    if (tk_handler_p->dash_set)
    {
        argv[n++] = Tcl_NewStringObj("-dash", -1);
        argv[n++] = Tcl_NewStringObj(tk_handler_p->dash, -1);
    }
    argv[n++] = Tcl_NewStringObj("-fill", -1);
    argv[n++] = Tcl_NewStringObj(tk_handler_p->canvas_color, -1);
    argv[n++] = Tcl_NewStringObj("-width", -1);
    argv[n++] = Tcl_NewIntObj(lw);
    argv[n++] = Tcl_NewStringObj("-tags", -1);
    argv[n++] = Tcl_NewStringObj(tk_handler_p->front_layer ? ccp_tag_front
                                                           : ccp_tag_back, -1);
    tkc_eval(tk_handler_p, n, argv);
    for (n = 0; n < 18; n++)
        if (argv[n])
            Tcl_DecrRefCount(argv[n]);
}

static void tkc_oval(Tk_handler_p tk_handler_p, int x0, int y0, int x1, int y1,
                     int filled)
{
    Tcl_Obj *argv[24] = {0};
    int n = 0;

    argv[n++] = Tcl_NewStringObj(tk_handler_p->canvas_path, -1);
    argv[n++] = Tcl_NewStringObj("create", -1);
    argv[n++] = Tcl_NewStringObj("oval", -1);
    argv[n++] = Tcl_NewIntObj(x0);
    argv[n++] = Tcl_NewIntObj(y0);
    argv[n++] = Tcl_NewIntObj(x1);
    argv[n++] = Tcl_NewIntObj(y1);
    tkc_item_options(tk_handler_p, argv, &n, filled,
                     tk_handler_p->canvas_color);
    tkc_eval(tk_handler_p, n, argv);
    for (n = 0; n < 24; n++)
        if (argv[n])
            Tcl_DecrRefCount(argv[n]);
}

static void tkc_rect(Tk_handler_p tk_handler_p, int x0, int y0, int x1, int y1,
                     int filled, const char *color)
{
    Tcl_Obj *argv[24] = {0};
    int n = 0;

    argv[n++] = Tcl_NewStringObj(tk_handler_p->canvas_path, -1);
    argv[n++] = Tcl_NewStringObj("create", -1);
    argv[n++] = Tcl_NewStringObj("rectangle", -1);
    argv[n++] = Tcl_NewIntObj(x0);
    argv[n++] = Tcl_NewIntObj(y0);
    argv[n++] = Tcl_NewIntObj(x1);
    argv[n++] = Tcl_NewIntObj(y1);
    tkc_item_options(tk_handler_p, argv, &n, filled, color);
    tkc_eval(tk_handler_p, n, argv);
    for (n = 0; n < 24; n++)
        if (argv[n])
            Tcl_DecrRefCount(argv[n]);
}

static void tkc_text(Tk_handler_p tk_handler_p, CcpnString text, int x, int y)
{
/*
    (x, y) is the TOP-LEFT corner of the text box, exactly as Tk_DrawChars
    uses it on X11 - the caller has already applied its anchor factors to
    the coordinates, so the canvas item is always placed with the "nw"
    (box top-left) anchor.
*/
    Tcl_Obj *argv[20] = {0};
    int n = 0;

    argv[n++] = Tcl_NewStringObj(tk_handler_p->canvas_path, -1);
    argv[n++] = Tcl_NewStringObj("create", -1);
    argv[n++] = Tcl_NewStringObj("text", -1);
    argv[n++] = Tcl_NewIntObj(x);
    argv[n++] = Tcl_NewIntObj(y);
    argv[n++] = Tcl_NewStringObj("-text", -1);
    argv[n++] = Tcl_NewStringObj(text, -1);
    argv[n++] = Tcl_NewStringObj("-font", -1);
    argv[n++] = Tcl_NewStringObj(tk_handler_p->font_spec ? tk_handler_p->font_spec
                                                         : "", -1);
    argv[n++] = Tcl_NewStringObj("-fill", -1);
    argv[n++] = Tcl_NewStringObj(tk_handler_p->canvas_color, -1);
    argv[n++] = Tcl_NewStringObj("-anchor", -1);
    argv[n++] = Tcl_NewStringObj("nw", -1);
    argv[n++] = Tcl_NewStringObj("-tags", -1);
    argv[n++] = Tcl_NewStringObj(tk_handler_p->front_layer ? ccp_tag_front
                                                           : ccp_tag_back, -1);
    tkc_eval(tk_handler_p, n, argv);
    for (n = 0; n < 20; n++)
        if (argv[n])
            Tcl_DecrRefCount(argv[n]);
}

#endif /* __APPLE__ */

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

Tk_handler new_tk_handler(Tcl_Interp *interp, Tk_Window tk_win, CcpnString win_path)
{
    int i;
    Tk_handler tk_handler;
    Tk_handler_p tk_handler_p;
#ifdef __APPLE__
    char path[300];
    Tcl_Obj *argv[24] = {0};
    int n;
#else
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
        0x88, 0x88, 0x88, 0x88, 0x22, 0x22, 0x22, 0x22
    };
#endif

    Tk_MakeWindowExist(tk_win);

    MALLOC_NEW(tk_handler_p, struct Tk_handler_p, 1);

    tk_handler_p->tk_win = tk_win;
    tk_handler_p->interp = interp;
    tk_handler_p->display = NULL;
    tk_handler_p->gc = (GC) 0;
    tk_handler_p->pixmap = (Pixmap) 0;
    tk_handler_p->drawable = 0;
    tk_handler_p->font_name = NULL;
    tk_handler_p->font_size = 10;
    tk_handler_p->font = NULL;
    for (i = 0; i < NCOLORS; i++)
        tk_handler_p->color[i] = -1;
    tk_handler_p->xcolor = NULL;
    tk_handler_p->is_double_buffer = CCPN_TRUE;

#ifdef __APPLE__
    {
        static int canvas_seq = 0;
        char child_name[40];

        snprintf(child_name, sizeof child_name, "ccpCanvas%d", ++canvas_seq);
        sprintf(path, "%s.%s", win_path, child_name);

        /* MALLOC_NEW: returns NULL from this function on OOM */
        MALLOC_NEW(tk_handler_p->canvas_path, char, strlen(path) + 1);
        strcpy(tk_handler_p->canvas_path, path);

        tk_handler_p->font_spec = NULL;
        tk_handler_p->front_layer = 0;
        tk_handler_p->cur_line_width = 1.0;
        tk_handler_p->dash_set = CCPN_FALSE;
        sprintf(tk_handler_p->canvas_color, "#ffffff");

        /* borderless canvas filling the widget (place: agnostic to the
           layout manager the parent uses).  Tk 9 class commands take the
           widget's FULL path as the single path argument ("canvas
           <path> ?options?"); the Tk 8 "canvas <parent> <name>" split
           is rejected (the first arg is the path to create). */
        n = 0;
        argv[n++] = Tcl_NewStringObj("canvas", -1);
        argv[n++] = Tcl_NewStringObj(path, -1);
        argv[n++] = Tcl_NewStringObj("-bd", -1);
        argv[n++] = Tcl_NewStringObj("0", -1);
        argv[n++] = Tcl_NewStringObj("-highlightthickness", -1);
        argv[n++] = Tcl_NewStringObj("0", -1);
        argv[n++] = Tcl_NewStringObj("-bg", -1);
        argv[n++] = Tcl_NewStringObj("#ffffff", -1);
        /* The canvas sits over the host and covers it 1:1. It must never take
           keyboard focus, or the app's key commands (routed to the focused
           WindowCanvas frame) would be lost after the first click on it. */
        argv[n++] = Tcl_NewStringObj("-takefocus", -1);
        argv[n++] = Tcl_NewStringObj("0", -1);
        tkc_eval(tk_handler_p, n, argv);
        for (i = 0; i < n; i++)
            Tcl_DecrRefCount(argv[i]);

        /* geometry managers are standalone commands: "place <path> ..."
           (never a method of the widget, on any Tk version) */
        n = 0;
        argv[n++] = Tcl_NewStringObj("place", -1);
        argv[n++] = Tcl_NewStringObj(path, -1);
        argv[n++] = Tcl_NewStringObj("-x", -1);
        argv[n++] = Tcl_NewStringObj("0", -1);
        argv[n++] = Tcl_NewStringObj("-y", -1);
        argv[n++] = Tcl_NewStringObj("0", -1);
        argv[n++] = Tcl_NewStringObj("-relwidth", -1);
        argv[n++] = Tcl_NewStringObj("1.0", -1);
        argv[n++] = Tcl_NewStringObj("-relheight", -1);
        argv[n++] = Tcl_NewStringObj("1.0", -1);
        tkc_eval(tk_handler_p, n, argv);
        for (i = 0; i < n; i++)
            Tcl_DecrRefCount(argv[i]);

        /* The canvas covers the host widget 1:1, so the pointer is always
           over the CANVAS - Tk dispatches pointer events through the
           canvas bindtags only, and the bindings the app registers on the
           host widget (crosshair <Motion>, peak-select
           <ButtonPress-1> + <B1-Motion> + <ButtonRelease-1> and their
           Shift/Control variants, drag-<Button-2/3>, right-click menu,
           <Enter>/<Leave> ...) never fire.  Replay each pointer event on
           the host widget with "event generate": the host's bindings run
           with event.widget = host (the app reads its attributes off the
           host there) and identical coordinates (the canvas sits at 0,0).
           The BASE type to generate is passed per binding - "event
           generate" does not accept modifier-qualified sequences.
           Configure/Expose/keyboard are NOT replayed: the host (not the
           focus-stealing canvas, see -takefocus 0) receives its own real
           ones, so keyboard and geometry events are unaffected. */
        /* Defined with Tcl_Eval on a full source string, NOT tkc_eval /
           Tcl_EvalObjv.  Passing the multi-command BODY to Tcl_EvalObjv as
           one argv word pre-compiles it into a SINGLE literal token, and
           `proc` then stores/executes the whole body as one command
           ("invalid command name <body>") - the crosshair dies silently.
           Tcl_Eval parses the body as an ordinary script.  (The single-level
           canvas/place/bind calls above are fine with Tcl_EvalObjv because
           every argv word there is just a value.)
           The body never touches the $event array: the coordinates,
           button and modifier state are substituted into the binding
           SCRIPT itself (%X %Y %s - the core Tk binding mechanism,
           evaluated where the event is actually dispatched) and passed
           in as ONE braced list "{<win> <seq> x y btn state}".  The
           earlier proc-scope `global event` read of event(x)/event(y)
           resolved to an EMPTY set for real pointer events, so every
           replay arrived at the host with the default x=0 y=0 (the
           stuck-at-origin crosshair).  Only the fields the app's
           handlers read are replayed (-x -y -state -button); -detail /
           -rootx / -rooty are omitted (unused, and event(detail) is an
           integer Notify code that `event generate` rejects - it
           requires the name).  With the CCP_FWD_DIAG env var set, the
           first forwarded events are recorded to
           $TMPDIR/ccp-fwd-diag.log for on-site diagnosis.  `{*}$opts`
           handles the empty-options case that a trailing empty argument
           would not. */
        static const char ccp_forward_proc_src[] =
"proc ccp_canvas_forward {pt} {\n"
"    lassign $pt w type x y btn state\n"
"    if {$state eq {}} {set state 0}\n"
"    if {[string match <B*-Motion> $type]} {\n"
"        # Tk 9.0.4 mac: a generated B-Motion whose -state lacks the\n"
"        # button's mask is demoted to a plain Motion (the host's\n"
"        # <B1-Motion> etc. never fires).  OR the mask in (idempotent if\n"
"        # the real %s already carries it) to keep the drag event type.\n"
"        set state [expr {$state | (1 << (7 + $btn))}]\n"
"    }\n"
"    set opts [list -x $x -y $y]\n"
"    if {$state ne {}} {lappend opts -state $state}\n"
"    if {$btn > 0 && [string match <Button* $type]} {lappend opts -button $btn}\n"
"    if {[info exists ::env(CCP_FWD_DIAG)]} {\n"
"        if {![info exists ::env(TMPDIR)]} {set ::env(TMPDIR) /tmp}\n"
"        if {![info exists ::ccp_fwd_diag_n]} {set ::ccp_fwd_diag_n 0}\n"
"        if {$::ccp_fwd_diag_n < 50} {\n"
"            incr ::ccp_fwd_diag_n\n"
"            catch {\n"
"                set f [open [file join $::env(TMPDIR) ccp-fwd-diag.log] a]\n"
"                puts $f [list $type \"->\" $w \"x=$x\" \"y=$y\" \"button=$btn\" \"state=$state\" \"opts=[$opts]\"]\n"
"                close $f\n"
"            }\n"
"        }\n"
"    }\n"
"    catch {event generate $w $type {*}$opts}\n"
"}\n";
        if (Tcl_Eval(tk_handler_p->interp, ccp_forward_proc_src) != TCL_OK)
        {
            const char *msg = Tcl_GetStringResult(tk_handler_p->interp);

            fprintf(stderr, "TkHandler forwarder-proc error: %s\n",
                    msg ? msg : "(no message)");
        }

        {
            /* Pointer sequences forwarded to the host.  The three button
               families below (press / B-Motion / release for buttons 1-3,
               with none/Shift/Control/Control-Shift modifiers) cover
               click-select, drag/translate and the right-click menu.
               <Double-1> is NOT forwarded (it cannot be `event generate`d;
               the host still receives its Double-click from the two forwarded
               ButtonPress-1 events).  Scroll-wheel is NOT forwarded here: its
               handler lives on an ancestor window (not the host) and the
               direction would be lost in the replay anyway. */
            static const char *mods[] = { "", "Control-", "Shift-", "Control-Shift-" };
            static const char *plain_seq[] = { "<Motion>", "<Enter>", "<Leave>" };
            char seq[80];
            char gseq[40];
            char script[400];
            int mi, pi, b, e3;

            for (pi = 0; pi < (int) (sizeof plain_seq / sizeof plain_seq[0]); pi++)
            {
                n = 0;
                argv[n++] = Tcl_NewStringObj("bind", -1);
                argv[n++] = Tcl_NewStringObj(path, -1);
                argv[n++] = Tcl_NewStringObj(plain_seq[pi], -1);
                /* %X %Y %s are substituted by the binding machinery at
                   dispatch time (no $event in the proc at all) */
                sprintf(script, "ccp_canvas_forward {%s %s %%X %%Y 0 %%s}", win_path, plain_seq[pi]);
                argv[n++] = Tcl_NewStringObj(script, -1);
                tkc_eval(tk_handler_p, n, argv);
                for (i = 0; i < n; i++)
                    Tcl_DecrRefCount(argv[i]);
            }
            for (mi = 0; mi < 4; mi++)
                for (b = 1; b <= 3; b++)
                    for (e3 = 0; e3 < 3; e3++)
                    {
                        if (e3 == 0)
                            snprintf(seq, sizeof seq, "<%sButtonPress-%d>", mods[mi], b);
                        else if (e3 == 1)
                            snprintf(seq, sizeof seq, "<%sB%d-Motion>", mods[mi], b);
                        else
                            snprintf(seq, sizeof seq, "<%sButtonRelease-%d>", mods[mi], b);
                        n = 0;
                        argv[n++] = Tcl_NewStringObj("bind", -1);
                        argv[n++] = Tcl_NewStringObj(path, -1);
                        argv[n++] = Tcl_NewStringObj(seq, -1);
                        if (e3 == 1)
                            snprintf(gseq, sizeof gseq, "<B%d-Motion>", b);
                        else if (e3 == 0)
                            snprintf(gseq, sizeof gseq, "<ButtonPress-%d>", b);
                        else
                            snprintf(gseq, sizeof gseq, "<ButtonRelease-%d>", b);
                        sprintf(script, "ccp_canvas_forward {%s %s %%X %%Y %d %%s}", win_path, gseq, b);
                        argv[n++] = Tcl_NewStringObj(script, -1);
                        tkc_eval(tk_handler_p, n, argv);
                        for (i = 0; i < n; i++)
                            Tcl_DecrRefCount(argv[i]);
                    }
        }

    }
#else
    display = Tk_Display(tk_win);

    Tk_DefineBitmap(interp, Tk_GetUid("stipple_data"), (char *) stipple_data, 32, 32);
    gcv.stipple = Tk_GetBitmap(interp, tk_win, Tk_GetUid("stipple_data"));

    gc = XCreateGC(display, Tk_WindowId(tk_win), GCForeground | GCBackground | GCStipple, &gcv);

    if (!gc)
    {
        FREE(tk_handler_p, struct Tk_handler_p);
        return NULL;
    }

    tk_handler_p->display = display;
    tk_handler_p->drawable = Tk_WindowId(tk_win);
    tk_handler_p->gc = gc;
#endif

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

#ifdef __APPLE__
        /* NOTE: the canvas WIDGET destroy is not done here - it is a Tcl
           call, and this function also runs from the Python destructor
           during Py_Finalize, by which time the interpreter is already
           gone (Tcl calls crash on a finalized interp).  The Python layer
           calls destroy_tk_canvas() for that, guarded by
           Py_IsFinalizing().  Everything here is plain C. */
        if (tk_handler_p->canvas_path)
        {
            FREE(tk_handler_p->canvas_path, char);
            tk_handler_p->canvas_path = NULL;
        }

        if (tk_handler_p->font_spec)
        {
            FREE(tk_handler_p->font_spec, char);
            tk_handler_p->font_spec = NULL;
        }
#else
	XFreeGC(tk_handler_p->display, tk_handler_p->gc);

        if (tk_handler_p->pixmap)
        {
            Tk_FreePixmap(tk_handler_p->display, tk_handler_p->pixmap);
	}

        if (tk_handler_p->xcolor)
        {
            Tk_FreeColor(tk_handler_p->xcolor);
	}
#endif

        if (tk_handler_p->font)
        {
            Tk_FreeFont(tk_handler_p->font);
            tk_handler_p->font = NULL;
        }

	FREE(tk_handler_p->font_name, char);

	FREE(tk_handler_p, struct Tk_handler_p);
    }
}

#ifdef __APPLE__
/*
    Destroy the canvas child widget (only while the interpreter is alive -
    the Python layer checks Py_IsFinalizing() before calling this; Tcl
    calls crash on a finalized interp).  If the parent widget has already
    gone with the window, the error is harmlessly swallowed by tkc_eval.
*/
void destroy_tk_canvas(Tk_handler tk_handler)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

    if (tk_handler && tk_handler_p->canvas_path && tk_handler_p->interp)
    {
        Tcl_Obj *argv[2] = { 0 };

        /* "destroy" is a standalone command, not a widget method:
           "destroy <path>" (the reverse order is "bad option destroy") */
        argv[0] = Tcl_NewStringObj("destroy", -1);
        argv[1] = Tcl_NewStringObj(tk_handler_p->canvas_path, -1);
        tkc_eval(tk_handler_p, 2, argv);
        Tcl_DecrRefCount(argv[0]);
        Tcl_DecrRefCount(argv[1]);
    }
}
#endif

void resize_tk_handler(Tk_handler tk_handler, int width, int height)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

#ifdef __APPLE__
    /* the canvas track-fits tk_win; only the coordinate scaling changes */
    if ((tk_handler_p->width == width) && (tk_handler_p->height == height))
        return;

    tk_handler_p->width = width;
    tk_handler_p->height = height;
#else
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

/*
    printf("resize_tk_handler: pixmap = %d\n", tk_handler_p->pixmap);
*/
#endif
}

void expose_tk_handler(Tk_handler tk_handler, int x, int y, int w, int h)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

/*
    printf("expose_tk_handler: x = %d, y = %d, w = %d, h = %d\n", x, y, w, h);
*/
    if (!make_current_tk_handler(tk_handler))
      return;

#ifdef __APPLE__
    {
        char bg[8];

        tkc_fill_color(bg, tk_handler_p->background);
        tkc_rect(tk_handler_p, x, y, x + w, y + h, CCPN_TRUE, bg);
    }
#else
    set_color_tk_handler(tk_handler, tk_handler_p->background);
    XFillRectangle(tk_handler_p->display, tk_handler_p->drawable,
                                tk_handler_p->gc, x, y, w, h);
#endif
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

#ifdef __APPLE__
    /* canvas items become visible on the next Tk event-loop turn, i.e.
       once this Python draw pass has finished - nothing to copy */
    (void) tk_handler_p;
#else
    Window win = Tk_WindowId(tk_handler_p->tk_win);

/*
    printf("swap_buffers_tk_handler: %d %d %d\n", win, tk_handler_p->drawable, tk_handler_p->pixmap);
    if (!make_current_tk_handler(tk_handler))
      return;
*/
    if (tk_handler_p->drawable == tk_handler_p->pixmap)
        XCopyArea(tk_handler_p->display, tk_handler_p->pixmap, win, tk_handler_p->gc,
                        0, 0, tk_handler_p->width, tk_handler_p->height, 0, 0);
#endif
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
#ifdef __APPLE__
    int x;
#else
    XGCValues gcv;
#endif

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

#ifdef __APPLE__
    /* same look as the X11 version: two offset passes of dashed verticals
       (simulating the stipple) plus a solid border, all on the front layer */
    for (x = xx0 + 1; x < xx1; x += 4)
        tkc_line(tk_handler_p, x, yy0 + 1, x, yy1 - 1);

    for (x = xx0 + 3; x < xx1; x += 4)
        tkc_line(tk_handler_p, x, yy0 + 1, x, yy1 - 1);

    tkc_rect(tk_handler_p, xx0, yy0, xx0 + w, yy0 + h, CCPN_FALSE,
             tk_handler_p->canvas_color);
#else
    // NOTE:ED - change to draw a series of dotted lines to simulate stipple - doesn't work in MacOS
    int dash_length = 1;
    int gap_length = 3;
    int dash_offset = (gap_length+dash_length)/2;
    int ndashes = 2;
    char dash_list[2];
    int x;

    dash_list[0] = (char) dash_length;
    dash_list[1] = (char) gap_length;

    gcv.line_style = LineOnOffDash;
    XChangeGC(tk_handler_p->display, tk_handler_p->gc, GCLineStyle, &gcv);

    // draw alternating lines to give offset to stippling
    XSetDashes(tk_handler_p->display, tk_handler_p->gc, 0, dash_list, sizeof(dash_list));
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
#endif

    finish_xor_tk_handler(tk_handler);
}

void start_xor_tk_handler(Tk_handler tk_handler)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
/*
    printf("start_xor_tk_handler\n");
    clear_xor_tk_handler(tk_handler);
*/
#ifdef __APPLE__
    /* canvas: "XOR" marks are managed front-layer items, so a fresh pass
       just clears the front layer and draws again */
    start_front_tk_handler(tk_handler);
#else
    swap_buffers_tk_handler(tk_handler);
/* below is because Solaris seems to require equiv to do xor (not sure why) */
#ifdef XOR_IS_EQUIV
    XSetFunction(tk_handler_p->display, tk_handler_p->gc, GXequiv);
#else
    XSetFunction(tk_handler_p->display, tk_handler_p->gc, GXxor);
#endif
    start_front_tk_handler(tk_handler);
#endif
}

void finish_xor_tk_handler(Tk_handler tk_handler)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

#ifdef __APPLE__
    /* front items stay on top until the next frame/crosshair pass clears
       the front layer - mirroring how the X11 XOR marks survive until the
       next swap; back-buffer drawing resumes */
    tk_handler_p->front_layer = 0;
#else
    XSetFunction(tk_handler_p->display, tk_handler_p->gc, GXcopy);
/*
    swap_buffers_tk_handler(tk_handler);
*/
/* TBD: should only set to start_back if was this to start with */
/* but whole xor'ing procedure implemented here doesn't work unless using back to start with */
    if (tk_handler_p->is_double_buffer)
	start_back_tk_handler(tk_handler);
#endif
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

#ifdef __APPLE__
/*
    printf("start_front_tk_handler\n");
*/
    tkc_delete_layer(tk_handler_p, ccp_tag_front);
    tk_handler_p->front_layer = 1;
#else
    Window win = Tk_WindowId(tk_handler_p->tk_win);

/*
    printf("start_front_tk_handler\n");
*/
    tk_handler_p->drawable = win;
#endif
}

void start_back_tk_handler(Tk_handler tk_handler)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

#ifdef __APPLE__
/*
    printf("start_back_tk_handler\n");
*/
/*
    A new frame invalidates the back layer; clear the front layer as well
    or stale crosshair items would sit above the new frame (the GUI
    redraws crosshairs after every frame, so none is ever missed).
*/
    tkc_delete_layer(tk_handler_p, ccp_tag_back);
    tkc_delete_layer(tk_handler_p, ccp_tag_front);
    tk_handler_p->front_layer = 0;
#else
/*
    printf("start_back_tk_handler\n");
*/
    if (tk_handler_p->pixmap)
        tk_handler_p->drawable = (Drawable) tk_handler_p->pixmap;
    else
        start_front_tk_handler(tk_handler);
#endif
}

void end_back_tk_handler(Tk_handler tk_handler)
{
}

void set_background_tk_handler(Tk_handler tk_handler, float *background)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

    COPY_VECTOR(tk_handler_p->background, background, NCOLORS);

#ifdef __APPLE__
    {
        Tcl_Obj *argv[4] = {0};
        char bg[8];
        int n = 0;

        tkc_fill_color(bg, background);
        argv[n++] = Tcl_NewStringObj(tk_handler_p->canvas_path, -1);
        argv[n++] = Tcl_NewStringObj("configure", -1);
        argv[n++] = Tcl_NewStringObj("-bg", -1);
        argv[n++] = Tcl_NewStringObj(bg, -1);
        tkc_eval(tk_handler_p, n, argv);
        for (n = 0; n < 4; n++)
            if (argv[n])
                Tcl_DecrRefCount(argv[n]);
    }
#endif
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

#ifdef __APPLE__
    tkc_text(tk_handler_p, text, s, t);
#else
    Tk_DrawChars(tk_handler_p->display, tk_handler_p->drawable,
                        tk_handler_p->gc, font, text, strlen(text), s, t);
#endif

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

#ifdef __APPLE__
    tkc_line(tk_handler_p, s0, t0, s1, t1);
#else
    XDrawLine(tk_handler_p->display, tk_handler_p->drawable,
                                tk_handler_p->gc, s0, t0, s1, t1);
/**
    XDrawLine(tk_handler_p->display, win,
                                tk_handler_p->gc, s0, t0, s1, t1);
**/
#endif
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

#ifdef __APPLE__
    if (rxx > 0 && ryy > 0)
        tkc_oval(tk_handler_p, xx - rxx, yy - ryy, xx + rxx, yy + ryy, CCPN_TRUE);
#else
    XFillArc(tk_handler_p->display, tk_handler_p->drawable,
                tk_handler_p->gc, xx-rxx, yy-ryy, 2*rxx, 2*ryy, 64*0, 64*360);
#endif
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

#ifdef __APPLE__
    if (rxx > 0 && ryy > 0)
        tkc_oval(tk_handler_p, xx - rxx, yy - ryy, xx + rxx, yy + ryy, CCPN_FALSE);
#else
    XDrawArc(tk_handler_p->display, tk_handler_p->drawable,
                tk_handler_p->gc, xx-rxx, yy-ryy, 2*rxx, 2*ryy, 64*0, 64*360);
#endif
}

void draw_polyline_tk_handler(Tk_handler tk_handler, Poly_line polyline)
{
    int i, n = polyline->nvertices;
    Point2f *v = polyline->vertices;
#ifdef __APPLE__
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    int lw;

    /* one canvas item for the WHOLE polyline - far fewer Tcl calls than
       one XDrawLine (or one create_line) per segment */
    {
        int npoints = polyline->closed ? (2 * n + 2) : 2 * n;
        Tcl_Obj **argv;

        lw = (int) (tk_handler_p->cur_line_width + 0.5);
        if (lw < 1)
            lw = 1;

        /* plain malloc: this function returns void, the MALLOC_* macros
           would `return` on OOM */
        argv = (Tcl_Obj **) malloc((3 + npoints + 10) * sizeof (Tcl_Obj *));
        if (!argv)
            return;
        argv[0] = Tcl_NewStringObj(tk_handler_p->canvas_path, -1);
        argv[1] = Tcl_NewStringObj("create", -1);
        argv[2] = Tcl_NewStringObj("line", -1);
        for (i = 0; i < n; i++)
        {
            argv[3 + 2 * i] = Tcl_NewIntObj(CONVERT_X(v[i].x));
            argv[3 + 2 * i + 1] = Tcl_NewIntObj(CONVERT_Y(v[i].y));
        }
        {
            int k = 3 + npoints;

            if (polyline->closed)
            {
                argv[k++] = Tcl_NewIntObj(CONVERT_X(v[0].x));
                argv[k++] = Tcl_NewIntObj(CONVERT_Y(v[0].y));
            }
            argv[k++] = Tcl_NewStringObj("-fill", -1);
            argv[k++] = Tcl_NewStringObj(tk_handler_p->canvas_color, -1);
            argv[k++] = Tcl_NewStringObj("-width", -1);
            argv[k++] = Tcl_NewIntObj(lw);
            argv[k++] = Tcl_NewStringObj("-tags", -1);
            argv[k++] = Tcl_NewStringObj(tk_handler_p->front_layer ? ccp_tag_front
                                                                   : ccp_tag_back, -1);
            tkc_eval(tk_handler_p, k, argv);
            for (i = 0; i < k; i++)
                Tcl_DecrRefCount(argv[i]);
        }
        FREE(argv, Tcl_Obj *);
    }
#else
    float x0, y0, x1, y1;

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
#endif

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

#ifdef __APPLE__
    /* canvas: the dash pattern is per-item, carried as transient context
       (clipped re-entry goes through tk_draw_line -> draw_line) */
    snprintf(tk_handler_p->dash, sizeof tk_handler_p->dash, "%d %d",
             dash_length, gap_length);
    tk_handler_p->dash_set = CCPN_TRUE;
    draw_line_tk_handler(tk_handler, x0, y0, x1, y1);
    tk_handler_p->dash_set = CCPN_FALSE;
#else
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
#endif
}

void draw_clipped_dash_line_tk_handler(Tk_handler tk_handler,
			float x0, float y0, float x1, float y1,
			int dash_length, int gap_length)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

#ifdef __APPLE__
    snprintf(tk_handler_p->dash, sizeof tk_handler_p->dash, "%d %d",
             dash_length, gap_length);
    tk_handler_p->dash_set = CCPN_TRUE;
    draw_clipped_line_tk_handler(tk_handler, x0, y0, x1, y1);
    tk_handler_p->dash_set = CCPN_FALSE;
#else
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
#endif
}

#define  CONVERT_COLOR(t) ((unsigned short) MIN(65535, 65536*(t)))

void set_color_tk_handler(Tk_handler tk_handler, float *color)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;
    int i;
#ifndef __APPLE__
    XColor c, *cc;
#endif

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

#ifdef __APPLE__
    tkc_fill_color(tk_handler_p->canvas_color, color);
#else
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
#endif
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

#ifdef __APPLE__
    tk_handler_p->cur_line_width = line_width;
#else
    XGCValues gcv;

    gcv.line_width = (unsigned int) line_width;
    XChangeGC(tk_handler_p->display, tk_handler_p->gc, GCLineWidth, &gcv);
#endif
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

#ifdef __APPLE__
    /* canvas -font takes the same spec Tk_GetFont used */
    if (tk_handler_p->font_spec)
        FREE(tk_handler_p->font_spec, char);
    STRING_MALLOC_COPY(tk_handler_p->font_spec, fontstring);
#else
    XSetFont(tk_handler_p->display, tk_handler_p->gc, Tk_FontId(font));
#endif

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
#ifdef __APPLE__
    Tcl_Obj *argv[18] = {0};
    int n = 0;

    argv[n++] = Tcl_NewStringObj(tk_handler_p->canvas_path, -1);
    argv[n++] = Tcl_NewStringObj("create", -1);
    argv[n++] = Tcl_NewStringObj("polygon", -1);
    argv[n++] = Tcl_NewIntObj(CONVERT_X(x0));
    argv[n++] = Tcl_NewIntObj(CONVERT_Y(y0));
    argv[n++] = Tcl_NewIntObj(CONVERT_X(x1));
    argv[n++] = Tcl_NewIntObj(CONVERT_Y(y1));
    argv[n++] = Tcl_NewIntObj(CONVERT_X(x2));
    argv[n++] = Tcl_NewIntObj(CONVERT_Y(y2));
    argv[n++] = Tcl_NewStringObj("-fill", -1);
    argv[n++] = Tcl_NewStringObj(tk_handler_p->canvas_color, -1);
    argv[n++] = Tcl_NewStringObj("-outline", -1);
    argv[n++] = Tcl_NewStringObj(tk_handler_p->canvas_color, -1);
    argv[n++] = Tcl_NewStringObj("-tags", -1);
    argv[n++] = Tcl_NewStringObj(tk_handler_p->front_layer ? ccp_tag_front
                                                           : ccp_tag_back, -1);
    tkc_eval(tk_handler_p, n, argv);
    for (n = 0; n < 18; n++)
        if (argv[n])
            Tcl_DecrRefCount(argv[n]);
#else
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
#endif
}

void set_is_double_buffer_tk_handler(Tk_handler tk_handler, Bool is_double_buffer)
{
    Tk_handler_p tk_handler_p = (Tk_handler_p) tk_handler;

    tk_handler_p->is_double_buffer = is_double_buffer;
}

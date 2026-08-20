
/*
======================COPYRIGHT/LICENSE START==========================

ps_handler.c: Part of the CcpNmr Analysis program

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
#include "ps_handler.h"

#include "clipping.h"
#include "utility.h"

#define FMT  "%.2f "
#define FMT2 "%.2f "

#define NCOLORS 3

static void ps_start_draw(Generic_ptr data) { /* TODO */
}

static void ps_end_draw(Generic_ptr data) { /* TODO */
}

static void ps_new_draw_range(Generic_ptr data, float x0, float y0, float x1, float y1, Bool clip) {
    new_range_ps_handler((Ps_handler)data, x0, y0, x1, y1);

    if (clip) clip_range_ps_handler(data, x0, y0, x1, y1);
}

static void ps_draw_line(Generic_ptr data, float x0, float y0, float x1, float y1) {
    draw_line_ps_handler((Ps_handler)data, x0, y0, x1, y1);
}

static void ps_draw_clipped_line(Generic_ptr data, float x0, float y0, float x1, float y1) {
    draw_clipped_line_ps_handler((Ps_handler)data, x0, y0, x1, y1);
}

static void ps_draw_polyline(Generic_ptr data, Poly_line polyline) { draw_polyline_ps_handler((Ps_handler)data, polyline); }

static void ps_draw_clipped_polyline(Generic_ptr data, Poly_line polyline) {
    draw_clipped_polyline_ps_handler((Ps_handler)data, polyline);
}

static void ps_draw_text(Generic_ptr data, CcpnString text, float x, float y, float a, float b) {
    draw_text_ps_handler((Ps_handler)data, text, x, y, a, b);
}

static void ps_set_draw_color(Generic_ptr data, float *color) { set_color_ps_handler((Ps_handler)data, color); }

static void ps_set_draw_font(Generic_ptr data, CcpnString name, int size) { set_font_ps_handler((Ps_handler)data, name, size); }

static void ps_set_line_style(Generic_ptr data, int line_style) { set_line_style_ps_handler((Ps_handler)data, line_style); }

static void ps_set_line_width(Generic_ptr data, float line_width) { set_line_width_ps_handler((Ps_handler)data, line_width); }

static void ps_fill_circle(Generic_ptr data, float x, float y, float r) { fill_circle_ps_handler((Ps_handler)data, x, y, r); }

static void ps_fill_ellipse(Generic_ptr data, float x, float y, float rx, float ry) {
    fill_ellipse_ps_handler((Ps_handler)data, x, y, rx, ry);
}

static void ps_draw_circle(Generic_ptr data, float x, float y, float r) { draw_circle_ps_handler((Ps_handler)data, x, y, r); }

static void ps_draw_ellipse(Generic_ptr data, float x, float y, float rx, float ry) {
    draw_ellipse_ps_handler((Ps_handler)data, x, y, rx, ry);
}

static void ps_get_background(Generic_ptr data, float *color) {
    int i;

    for (i = 0; i < NCOLORS; i++) color[i] = 1.0; /* white */
}

static void ps_get_region(Generic_ptr data, float *x0, float *y0, float *x1, float *y1) {
    get_region_ps_handler((Ps_handler)data, x0, y0, x1, y1);
}

static void ps_get_text_size(Generic_ptr data, CcpnString text, float *w, float *h) {
    get_text_size_ps_handler((Ps_handler)data, text, w, h);
}

static void ps_fill_triangle(Generic_ptr data, float x0, float y0, float x1, float y1, float x2, float y2) {
    fill_triangle_ps_handler((Ps_handler)data, x0, y0, x1, y1, x2, y2);
}

static Drawing_funcs drawing_funcs = {PAPER_DRAWING,     ps_start_draw,        ps_end_draw,      ps_new_draw_range,
                                      ps_draw_line,      ps_draw_clipped_line, ps_draw_polyline, ps_draw_clipped_polyline,
                                      ps_draw_text,      ps_set_draw_color,    ps_set_draw_font, ps_set_line_style,
                                      ps_set_line_width, ps_fill_circle,       ps_fill_ellipse,  ps_draw_circle,
                                      ps_draw_ellipse,   ps_get_background,    ps_get_region,    ps_get_text_size,
                                      ps_fill_triangle};

Drawing_funcs *ps_drawing_funcs(void) { return &drawing_funcs; }

Ps_handler new_ps_handler(CcpnString *file_name, float width, float height, CcpnString output_style) {
    Ps_handler ps_handler;

    MALLOC_NEW(ps_handler, struct Ps_handler, 1);

    strcpy(ps_handler->fp, file_name);

    //    ps_handler->fp = file_name;
    ps_handler->width = width;
    ps_handler->height = height;

    ps_handler->ax = 1.0;
    ps_handler->ay = 1.0;
    ps_handler->bx = 0.0;
    ps_handler->by = 0.0;

    if (equal_strings(output_style, "Color"))
        ps_handler->isColor = CCPN_TRUE;
    else
        ps_handler->isColor = CCPN_FALSE;

    strcpy(ps_handler->font, "Times-Roman"); /* arbitrary */
    ps_handler->fontsize = 12;               /* arbitrary */

    return ps_handler;
}

void delete_ps_handler(Ps_handler ps_handler) { FREE(ps_handler, struct Ps_handler); }

void new_range_ps_handler(Ps_handler ps_handler, float x0, float y0, float x1, float y1) {
    //    FILE *fp = ps_handler->fp;
    float width = ps_handler->width;
    float height = ps_handler->height;

    FILE *fp;
    OPEN_FOR_APPENDING(fp, ps_handler->fp);
    fprintf(fp, "PS_graphics_restore\n");
    fprintf(fp, "PS_graphics_save\n");

    /*
        printf("new_range_ps_handler1: %3.2f %3.2f %3.2f %3.2f\n", x0, y0, x1, y1);
    */
    ps_handler->x0 = MIN(x0, x1);
    ps_handler->y0 = MIN(y0, y1);
    ps_handler->x1 = MAX(x0, x1);
    ps_handler->y1 = MAX(y0, y1);

    x0 = ps_handler->x0;
    y0 = ps_handler->y0;
    x1 = ps_handler->x1;
    y1 = ps_handler->y1;

    ps_handler->ax = width / (x1 - x0);
    ps_handler->ay = height / (y1 - y0);
    ps_handler->bx = -ps_handler->ax * x0;
    ps_handler->by = -ps_handler->ay * y0;
    /*
        printf("new_range_ps_handler2: %3.2f %3.2f %3.2f %3.2f\n",
            ps_handler->ax, ps_handler->ay, ps_handler->bx, ps_handler->by);
    */

    fprintf(fp, FMT FMT "PS_translate\n", ps_handler->bx, ps_handler->by);

    /*  scaling causes problems with fonts, think again about this later
        fprintf(fp, FMT FMT "PS_scale\n", ps_handler->ax, ps_handler->ay);
    */
    FCLOSE(fp);
}

void clip_range_ps_handler(Ps_handler ps_handler, float x0, float y0, float x1, float y1) {
    FILE *fp;
    OPEN_FOR_APPENDING(fp, ps_handler->fp);
    //    fprintf(ps_handler->fp, FMT FMT FMT FMT "PS_rectangle_clip\n",
    //				ps_handler->ax*x0, ps_handler->ay*y0,
    //				ps_handler->ax*(x1-x0), ps_handler->ay*(y1-y0));
    fprintf(fp, FMT FMT FMT FMT "PS_rectangle_clip\n", ps_handler->ax * x0, ps_handler->ay * y0, ps_handler->ax * (x1 - x0),
            ps_handler->ay * (y1 - y0));
    FCLOSE(fp);
}

void draw_line_ps_handler(Ps_handler ps_handler, float x0, float y0, float x1, float y1) {
    float ax = ps_handler->ax;
    float ay = ps_handler->ay;
    //    FILE *fp = ps_handler->fp;

    FILE *fp;
    OPEN_FOR_APPENDING(fp, ps_handler->fp);
    fprintf(fp, "NP " FMT FMT "M ", ax * x0, ay * y0);
    fprintf(fp, FMT FMT "L ", ax * x1, ay * y1);
    fprintf(fp, "S\n");
    /*
        fprintf(ps_handler->fp, FMT FMT FMT FMT "PS_draw_line\n",
                                    ps_handler->ax*x0, ps_handler->ay*y0,
                                    ps_handler->ax*x1, ps_handler->ay*y1);
    */
    FCLOSE(fp);
}

void draw_clipped_line_ps_handler(Ps_handler ps_handler, float x0, float y0, float x1, float y1) {
    float xmin = ps_handler->x0;
    float ymin = ps_handler->y0;
    float xmax = ps_handler->x1;
    float ymax = ps_handler->y1;

    draw_clipped_line(x0, y0, x1, y1, &drawing_funcs, (Generic_ptr)ps_handler, xmin, ymin, xmax, ymax);
}

void fill_circle_ps_handler(Ps_handler ps_handler, float x, float y, float r) {
    float s = (float)sqrt((double)(ps_handler->ax * ps_handler->ay));

    x *= ps_handler->ax;
    y *= ps_handler->ay;
    r *= s;

    FILE *fp;
    OPEN_FOR_APPENDING(fp, ps_handler->fp);
    //    fprintf(ps_handler->fp, FMT FMT FMT "PS_fill_circle\n", x, y, r);
    fprintf(fp, FMT FMT FMT "PS_fill_circle\n", x, y, r);
    FCLOSE(fp)
}

void fill_ellipse_ps_handler(Ps_handler ps_handler, float x, float y, float rx, float ry) {
    float r = (float)sqrt((double)(rx * ry));

    fill_circle_ps_handler(ps_handler, x, y, r);
}

void draw_circle_ps_handler(Ps_handler ps_handler, float x, float y, float r) {
    float s = (float)sqrt((double)(ps_handler->ax * ps_handler->ay));

    x *= ps_handler->ax;
    y *= ps_handler->ay;
    r *= s;

    FILE *fp;
    OPEN_FOR_APPENDING(fp, ps_handler->fp);
    //    fprintf(ps_handler->fp, FMT FMT FMT "PS_draw_circle\n", x, y, r);
    fprintf(fp, FMT FMT FMT "PS_draw_circle\n", x, y, r);
    FCLOSE(fp)
}

void draw_ellipse_ps_handler(Ps_handler ps_handler, float x, float y, float rx, float ry) {
    float r = (float)sqrt((double)(rx * ry));

    draw_circle_ps_handler(ps_handler, x, y, r);
}

void draw_polyline_ps_handler(Ps_handler ps_handler, Poly_line polyline) {
    int i, n = polyline->nvertices;
    Point2f *v = polyline->vertices;
    float ax = ps_handler->ax;
    float ay = ps_handler->ay;
    //    FILE *fp = ps_handler->fp;

    if (n < 2) return;

    FILE *fp;
    OPEN_FOR_APPENDING(fp, ps_handler->fp);
    fprintf(fp, "NP " FMT FMT "M ", ax * v[0].x, ay * v[0].y);

    for (i = 1; i < n; i++) fprintf(fp, FMT FMT "L ", ax * v[i].x, ay * v[i].y);

    if (polyline->closed) fprintf(fp, "CP ");

    fprintf(fp, "S\n");
    /*
        for (i = 0; i < n-1; i++)
            draw_line_ps_handler(ps_handler, v[i].x, v[i].y, v[i+1].x, v[i+1].y);

        if (polyline->closed)
            draw_line_ps_handler(ps_handler, v[n-1].x, v[n-1].y, v[0].x, v[0].y);
    */
    FCLOSE(fp)
}

void draw_clipped_polyline_ps_handler(Ps_handler ps_handler, Poly_line polyline) {
    float xmin = ps_handler->x0;
    float ymin = ps_handler->y0;
    float xmax = ps_handler->x1;
    float ymax = ps_handler->y1;

    draw_clipped_polyline(polyline, &drawing_funcs, (Generic_ptr)ps_handler, xmin, ymin, xmax, ymax);
}

void draw_dash_line_ps_handler(Ps_handler ps_handler, float x0, float y0, float x1, float y1, int dash_length, int gap_length) {
    /* TODO */
    /* below is not right but it's better than nothing */
    set_line_style_ps_handler(ps_handler, DASHED_LINE_STYLE);
    draw_line_ps_handler(ps_handler, x0, y0, x1, y1);
    set_line_style_ps_handler(ps_handler, NORMAL_LINE_STYLE);
}

void draw_clipped_dash_line_ps_handler(Ps_handler ps_handler, float x0, float y0, float x1, float y1, int dash_length,
                                       int gap_length) {
    /* TODO */
    /* below is not right but it's better than nothing */
    set_line_style_ps_handler(ps_handler, DASHED_LINE_STYLE);
    draw_clipped_line_ps_handler(ps_handler, x0, y0, x1, y1);
    set_line_style_ps_handler(ps_handler, NORMAL_LINE_STYLE);
}

#define PS_TEXT_ESCAPE '\\'
#define PS_TEXT_START  '('
#define PS_TEXT_END    ')'

static void string_c_to_ps(CcpnString dest, CcpnString src) {
    while (*src) {
        if ((*src == PS_TEXT_ESCAPE) || (*src == PS_TEXT_START) || (*src == PS_TEXT_END)) *dest++ = PS_TEXT_ESCAPE;

        *dest++ = *src++;
    }

    *dest = 0;
}

void draw_text_ps_handler(Ps_handler ps_handler, CcpnString text, float x, float y, float a, float b) {
    Line line;

    FILE *fp;
    OPEN_FOR_APPENDING(fp, ps_handler->fp);
    string_c_to_ps(line, text);
    //    fprintf(ps_handler->fp, FMT FMT FMT2 FMT2 "(%s) PS_text\n",
    //		ps_handler->ax*x, ps_handler->ay*y, a, b, line);
    fprintf(fp, FMT FMT FMT2 FMT2 "(%s) PS_text\n", ps_handler->ax * x, ps_handler->ay * y, a, b, line);
    FCLOSE(fp)
}

void draw_dash_box_ps_handler(Ps_handler ps_handler, float x0, float y0, float x1, float y1) {
    int dash_length = 2, gap_length = 6;

    draw_clipped_dash_line_ps_handler(ps_handler, x0, y0, x0, y1, dash_length, gap_length);
    draw_clipped_dash_line_ps_handler(ps_handler, x0, y1, x1, y1, dash_length, gap_length);
    draw_clipped_dash_line_ps_handler(ps_handler, x1, y1, x1, y0, dash_length, gap_length);
    draw_clipped_dash_line_ps_handler(ps_handler, x1, y0, x0, y0, dash_length, gap_length);
}

void set_color_ps_handler(Ps_handler ps_handler, float *color) {
    float r = color[0], g = color[1], b = color[2];

    if (!ps_handler->isColor) {
        if (r + g + b < 2.99)
            r = g = b = 0.0;
        else
            r = g = b = 1.0;
    }

    FILE *fp;
    OPEN_FOR_APPENDING(fp, ps_handler->fp);
    //    fprintf(ps_handler->fp, FMT2 FMT2 FMT2 "PS_set_color\n", r, g, b);
    fprintf(fp, FMT2 FMT2 FMT2 "PS_set_color\n", r, g, b);
    FCLOSE(fp)
}

void set_black_ps_handler(Ps_handler ps_handler) {
    float color[NCOLORS];

    ZERO_VECTOR(color, NCOLORS);

    set_color_ps_handler(ps_handler, color);
}

void set_font_ps_handler(Ps_handler ps_handler, CcpnString name, int size) {
    if (equal_strings(name, "Times"))
        strcpy(ps_handler->font, "Times-Roman");
    else
        strncpy(ps_handler->font, name, LINE_SIZE - 1);
    ps_handler->font[LINE_SIZE - 1] = 0;
    ps_handler->fontsize = size;

    FILE *fp;
    OPEN_FOR_APPENDING(fp, ps_handler->fp);
    //    fprintf(ps_handler->fp, "/%s %d PS_select_font\n", ps_handler->font, size);
    fprintf(fp, "/%s %d PS_select_font\n", ps_handler->font, size);
    FCLOSE(fp)
}

void set_line_style_ps_handler(Ps_handler ps_handler, int line_style) {
    FILE *fp;
    OPEN_FOR_APPENDING(fp, ps_handler->fp);
    //    if (line_style == NORMAL_LINE_STYLE)
    //	fprintf(ps_handler->fp, "PS_nondashed\n");
    //    else /* line_style == DASHED_LINE_STYLE */
    //	fprintf(ps_handler->fp, "PS_dashed\n");
    if (line_style == NORMAL_LINE_STYLE)
        fprintf(fp, "PS_nondashed\n");
    else /* line_style == DASHED_LINE_STYLE */
        fprintf(fp, "PS_dashed\n");
    FCLOSE(fp)
}

void set_line_width_ps_handler(Ps_handler ps_handler, float line_width) {
    FILE *fp;
    OPEN_FOR_APPENDING(fp, ps_handler->fp);
    //    fprintf(ps_handler->fp, FMT "PS_set_linewidth\n", line_width);
    fprintf(fp, FMT "PS_set_linewidth\n", line_width);
    FCLOSE(fp)
}

void get_region_ps_handler(Ps_handler ps_handler, float *x0, float *y0, float *x1, float *y1) {
    *x0 = ps_handler->x0;
    *y0 = ps_handler->y0;
    *x1 = ps_handler->x1;
    *y1 = ps_handler->y1;
}

void get_text_size_ps_handler(Ps_handler ps_handler, CcpnString text, float *w, float *h) {
    /* cannot really do this properly */
    int size = ps_handler->fontsize;

    *w = size * strlen(text);
    *h = size;
}

void fill_triangle_ps_handler(Ps_handler ps_handler, float x0, float y0, float x1, float y1, float x2, float y2) {
    float _x0 = x0 * ps_handler->ax;
    float _y0 = y0 * ps_handler->ay;
    float _x1 = x1 * ps_handler->ax;
    float _y1 = y1 * ps_handler->ay;
    float _x2 = x2 * ps_handler->ax;
    float _y2 = y2 * ps_handler->ay;

    //    printf("HEREafterTri: %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f\n", _x0, _y0, _x1, _y1, _x2, _y2);

    FILE *fp;
    OPEN_FOR_APPENDING(fp, ps_handler->fp);
    //    fprintf(ps_handler->fp, FMT FMT FMT FMT FMT FMT "PS_fill_triangle\n", x0, y0, x1, y1, x2, y2);
    fprintf(fp, FMT FMT FMT FMT FMT FMT "PS_fill_triangle\n", _x0, _y0, _x1, _y1, _x2, _y2);
    FCLOSE(fp)
}


/*
======================COPYRIGHT/LICENSE START==========================

pdf_handler.c: Part of the CcpNmr Analysis program

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
#include "pdf_handler.h"

#include "clipping.h"
#include "utility.h"

#define  FMT	"%.2f "
#define  FMT2	"%.2f "

#define  NCOLORS  3

static void pdf_start_draw(Generic_ptr data)
{
/* TODO */
}

static void pdf_end_draw(Generic_ptr data)
{
/* TODO */
}

static void pdf_new_draw_range(Generic_ptr data, float x0, float y0,
				float x1, float y1, Bool clip)
{
    new_range_pdf_handler((Pdf_handler) data, x0, y0, x1, y1);

    if (clip)
        clip_range_pdf_handler(data, x0, y0, x1, y1);
}

static void pdf_draw_line(Generic_ptr data, float x0, float y0, float x1, float y1)
{
    draw_line_pdf_handler((Pdf_handler) data, x0, y0, x1, y1);
}

static void pdf_draw_clipped_line(Generic_ptr data, float x0, float y0, float x1, float y1)
{
    draw_clipped_line_pdf_handler((Pdf_handler) data, x0, y0, x1, y1);
}

static void pdf_draw_polyline(Generic_ptr data, Poly_line polyline)
{
    draw_polyline_pdf_handler((Pdf_handler) data, polyline);
}

static void pdf_draw_clipped_polyline(Generic_ptr data, Poly_line polyline)
{
    draw_clipped_polyline_pdf_handler((Pdf_handler) data, polyline);
}

static void pdf_draw_text(Generic_ptr data, CcpnString text, float x, float y,
                                                        float a, float b)
{
    draw_text_pdf_handler((Pdf_handler) data, text, x, y, a, b);
}
 
static void pdf_set_draw_color(Generic_ptr data, float *color)
{
    set_color_pdf_handler((Pdf_handler) data, color);
}
 
static void pdf_set_draw_font(Generic_ptr data, CcpnString name, int size)
{
    set_font_pdf_handler((Pdf_handler) data, name, size);
}
 
static void pdf_set_line_style(Generic_ptr data, int line_style)
{
    set_line_style_pdf_handler((Pdf_handler) data, line_style);
}

static void pdf_set_line_width(Generic_ptr data, float line_width)
{
    set_line_width_pdf_handler((Pdf_handler) data, line_width);
}

static void pdf_fill_circle(Generic_ptr data, float x, float y, float r)
{
    fill_circle_pdf_handler((Pdf_handler) data, x, y, r);
}
 
static void pdf_fill_ellipse(Generic_ptr data, float x, float y, float rx, float ry)
{
    fill_ellipse_pdf_handler((Pdf_handler) data, x, y, rx, ry);
}
 
static void pdf_draw_circle(Generic_ptr data, float x, float y, float r)
{
    draw_circle_pdf_handler((Pdf_handler) data, x, y, r);
}
 
static void pdf_draw_ellipse(Generic_ptr data, float x, float y, float rx, float ry)
{
    draw_ellipse_pdf_handler((Pdf_handler) data, x, y, rx, ry);
}
 
static void pdf_get_background(Generic_ptr data, float *color)
{
    int i;

    for (i = 0; i < NCOLORS; i++)
        color[i] = 1.0; /* white */
}

static void pdf_get_region(Generic_ptr data, float *x0, float *y0,
                                                    float *x1, float *y1)
{
    get_region_pdf_handler((Pdf_handler) data, x0, y0, x1, y1);
}

static void pdf_get_text_size(Generic_ptr data, CcpnString text, float *w, float *h)
{
    get_text_size_pdf_handler((Pdf_handler) data, text, w, h);
}

static void pdf_fill_triangle(Generic_ptr data, float x0, float y0,
                                float x1, float y1, float x2, float y2)
{
    fill_triangle_pdf_handler((Pdf_handler) data, x0, y0, x1, y1, x2, y2);
}

static Drawing_funcs drawing_funcs = { PAPER_DRAWING,
	pdf_start_draw, pdf_end_draw, pdf_new_draw_range,
	pdf_draw_line, pdf_draw_clipped_line,
	pdf_draw_polyline, pdf_draw_clipped_polyline, pdf_draw_text,
	pdf_set_draw_color, pdf_set_draw_font, pdf_set_line_style,
        pdf_set_line_width, pdf_fill_circle, pdf_fill_ellipse,
	pdf_draw_circle, pdf_draw_ellipse, pdf_get_background,
	pdf_get_region, pdf_get_text_size, pdf_fill_triangle };
 
Drawing_funcs *pdf_drawing_funcs(void)
{
    return &drawing_funcs;
}

Pdf_handler new_pdf_handler(CcpnString *file_name, float width, float height,
                                                CcpnString output_style)
{
    Pdf_handler pdf_handler;

    MALLOC_NEW(pdf_handler, struct Pdf_handler, 1);

    strcpy(pdf_handler->fp, file_name);
//    strcat(pdf_handler->fp, "0");

//    pdf_handler->fp = fp;
    pdf_handler->width = width;
    pdf_handler->height = height;

    pdf_handler->ax = 1.0;
    pdf_handler->ay = 1.0;
    pdf_handler->bx = 0.0;
    pdf_handler->by = 0.0;

    if (equal_strings(output_style, "Color"))
        pdf_handler->isColor = CCPN_TRUE;
    else
        pdf_handler->isColor = CCPN_FALSE;

    set_font_pdf_handler(pdf_handler, "Times-Roman", 10);

    return pdf_handler;
}

void delete_pdf_handler(Pdf_handler pdf_handler)
{
    FREE(pdf_handler, struct Pdf_handler);
}

void new_range_pdf_handler(Pdf_handler pdf_handler, float x0, float y0,
					float x1, float y1)
{
//    FILE *fp = pdf_handler->fp;
    float width = pdf_handler->width;
    float height = pdf_handler->height;

    FILE *fp;
    OPEN_FOR_APPENDING(fp, pdf_handler->fp);
    fprintf(fp, "Q\n"); /* restore */
    fprintf(fp, "q\n"); /* save */

/*
    {
	fprintf(fp, "0 0 m\n");
	fprintf(fp, FMT "0 l\n", width);
	fprintf(fp, FMT FMT "l\n", width, height);
	fprintf(fp, "0 " FMT "l\n", height);
	fprintf(fp, "0 0 l\n");
	fprintf(fp, "W n\n");
    }
*/

/*
    printf("new_range_pdf_handler1: %3.2f %3.2f %3.2f %3.2f\n", x0, y0, x1, y1);
*/
    pdf_handler->x0 = MIN(x0, x1);
    pdf_handler->y0 = MIN(y0, y1);
    pdf_handler->x1 = MAX(x0, x1);
    pdf_handler->y1 = MAX(y0, y1);

    x0 = pdf_handler->x0;
    y0 = pdf_handler->y0;
    x1 = pdf_handler->x1;
    y1 = pdf_handler->y1;

    pdf_handler->ax = width / (x1 - x0);
    pdf_handler->ay = height / (y1 - y0);
    pdf_handler->bx = - pdf_handler->ax * x0;
    pdf_handler->by = - pdf_handler->ay * y0;
/*
    printf("new_range_pdf_handler2: %3.2f %3.2f %3.2f %3.2f\n", 
	pdf_handler->ax, pdf_handler->ay, pdf_handler->bx, pdf_handler->by);
*/

    fprintf(fp, "1 0 0 1 %3.2f %3.2f cm\n", pdf_handler->bx, pdf_handler->by);
    FCLOSE(fp)
}

void clip_range_pdf_handler(Pdf_handler pdf_handler, float x0, float y0,
					float x1, float y1)
{
//    FILE *fp = pdf_handler->fp;

    x0 *= pdf_handler->ax;
    x1 *= pdf_handler->ax;
    y0 *= pdf_handler->ay;
    y1 *= pdf_handler->ay;

    FILE *fp;
    OPEN_FOR_APPENDING(fp, pdf_handler->fp);
    fprintf(fp, FMT FMT "m\n", x0, y0);
    fprintf(fp, FMT FMT "l\n", x1, y0);
    fprintf(fp, FMT FMT "l\n", x1, y1);
    fprintf(fp, FMT FMT "l\n", x0, y1);
    fprintf(fp, FMT FMT "l\n", x0, y0);
    fprintf(fp, "W n\n");
    FCLOSE(fp)
}

void draw_line_pdf_handler(Pdf_handler pdf_handler, float x0, float y0,
							float x1, float y1)
{
    float ax = pdf_handler->ax;
    float ay = pdf_handler->ay;
//    FILE *fp = pdf_handler->fp;
/*
    printf("draw_line_pdf_handler0: %3.2f %3.2f\n", pdf_handler->ax, pdf_handler->ay);
    printf("draw_line_pdf_handler1: %3.2f %3.2f %3.2f %3.2f\n", x0, y0, x1, y1);
    printf("draw_line_pdf_handler2: %3.2f %3.2f %3.2f %3.2f\n",
	pdf_handler->ax*x0, pdf_handler->ay*y0, pdf_handler->ax*x1, pdf_handler->ay*y1);
*/

    FILE *fp;
    OPEN_FOR_APPENDING(fp, pdf_handler->fp);
    fprintf(fp, FMT FMT "m ", ax*x0, ay*y0);
    fprintf(fp, FMT FMT "l ", ax*x1, ay*y1);
    fprintf(fp, "S\n");
    FCLOSE(fp)
}

void draw_clipped_line_pdf_handler(Pdf_handler pdf_handler,
                        float x0, float y0, float x1, float y1)
{
    float xmin = pdf_handler->x0;
    float ymin = pdf_handler->y0;
    float xmax = pdf_handler->x1;
    float ymax = pdf_handler->y1;

    draw_clipped_line(x0, y0, x1, y1, &drawing_funcs,
                (Generic_ptr) pdf_handler, xmin, ymin, xmax, ymax);
}

#define  NCIRCLE_PTS  20
#define  DEGREES_PER_PT  (360 / NCIRCLE_PTS)
#define  RADIANS_PER_PT  (DEGREES_PER_PT * RADIAN_SCALE)

void fill_circle_pdf_handler(Pdf_handler pdf_handler, float x, float y, float r)
{
    fill_ellipse_pdf_handler(pdf_handler, x, y, r, r);
}

void fill_ellipse_pdf_handler(Pdf_handler pdf_handler, float x, float y, float rx, float ry)
{
    int i;
    float xx, yy, angle;
    FILE *fp;
    OPEN_FOR_APPENDING(fp, pdf_handler->fp);

    xx = x + rx;
    yy = y;
//    fprintf(pdf_handler->fp, FMT FMT " m\n",
//			pdf_handler->ax*xx, pdf_handler->ay*yy);
    fprintf(fp, FMT FMT " m\n",
			pdf_handler->ax*xx, pdf_handler->ay*yy);
    for (i = 1; i <= NCIRCLE_PTS; i++)
    {
        angle = i * RADIANS_PER_PT;
        xx = x + rx * (float) cos((double) angle);
        yy = y + ry * (float) sin((double) angle);
//        fprintf(pdf_handler->fp, FMT FMT " l\n",
//			pdf_handler->ax*xx, pdf_handler->ay*yy);
        fprintf(fp, FMT FMT " l\n",
			pdf_handler->ax*xx, pdf_handler->ay*yy);
    }

//    fprintf(pdf_handler->fp, "f\n");
    fprintf(fp, "f\n");
    FCLOSE(fp)
}

void draw_circle_pdf_handler(Pdf_handler pdf_handler, float x, float y, float r)
{
    draw_ellipse_pdf_handler(pdf_handler, x, y, r, r);
}

void draw_ellipse_pdf_handler(Pdf_handler pdf_handler, float x, float y, float rx, float ry)
{
    int i;
    float xx, yy, angle;
    FILE *fp;
    OPEN_FOR_APPENDING(fp, pdf_handler->fp);

    xx = x + rx;
    yy = y;
//    fprintf(pdf_handler->fp, FMT FMT " m\n",
//			pdf_handler->ax*xx, pdf_handler->ay*yy);
    fprintf(fp, FMT FMT " m\n",
			pdf_handler->ax*xx, pdf_handler->ay*yy);
    for (i = 1; i <= NCIRCLE_PTS; i++)
    {
        angle = i * RADIANS_PER_PT;
        xx = x + rx * (float) cos((double) angle);
        yy = y + ry * (float) sin((double) angle);
//        fprintf(pdf_handler->fp, FMT FMT " l\n",
//			pdf_handler->ax*xx, pdf_handler->ay*yy);
        fprintf(fp, FMT FMT " l\n",
			pdf_handler->ax*xx, pdf_handler->ay*yy);
    }

//    fprintf(pdf_handler->fp, "S\n");
    fprintf(fp, "S\n");
    FCLOSE(fp)
}

void draw_polyline_pdf_handler(Pdf_handler pdf_handler, Poly_line polyline)
{
    int i, n = polyline->nvertices;
    Point2f *v = polyline->vertices;
    float ax = pdf_handler->ax;
    float ay = pdf_handler->ay;
//    FILE *fp = pdf_handler->fp;
    FILE *fp;
    OPEN_FOR_APPENDING(fp, pdf_handler->fp);

/*
    for (i = 0; i < n-1; i++)
	draw_line_pdf_handler(pdf_handler, v[i].x, v[i].y, v[i+1].x, v[i+1].y);

    if (polyline->closed)
	draw_line_pdf_handler(pdf_handler, v[n-1].x, v[n-1].y, v[0].x, v[0].y);
*/
    if (n < 2)
        return;

    fprintf(fp, FMT FMT "m ", ax*v[0].x, ay*v[0].y);

    for (i = 1; i < n; i++)
        fprintf(fp, FMT FMT "l ", ax*v[i].x, ay*v[i].y);

    if (polyline->closed)
        fprintf(fp, "h ");
        //fprintf(fp, FMT FMT "l ", ax*v[0].x, ay*v[0].y);

    fprintf(fp, "S\n");
    FCLOSE(fp)
}

void draw_clipped_polyline_pdf_handler(Pdf_handler pdf_handler, Poly_line polyline)
{   
    float xmin = pdf_handler->x0;
    float ymin = pdf_handler->y0;
    float xmax = pdf_handler->x1;
    float ymax = pdf_handler->y1;

    draw_clipped_polyline(polyline, &drawing_funcs,
                (Generic_ptr) pdf_handler, xmin, ymin, xmax, ymax);
}

void draw_dash_line_pdf_handler(Pdf_handler pdf_handler,
			float x0, float y0, float x1, float y1,
			int dash_length, int gap_length)
{
/* TODO */
/* below is not right but it's better than nothing */
   set_line_style_pdf_handler(pdf_handler, DASHED_LINE_STYLE);
   draw_line_pdf_handler(pdf_handler, x0, y0, x1, y1);
   set_line_style_pdf_handler(pdf_handler, NORMAL_LINE_STYLE);
}

void draw_clipped_dash_line_pdf_handler(Pdf_handler pdf_handler,
			float x0, float y0, float x1, float y1,
			int dash_length, int gap_length)
{
/* TODO */
/* below is not right but it's better than nothing */
   set_line_style_pdf_handler(pdf_handler, DASHED_LINE_STYLE);
   draw_clipped_line_pdf_handler(pdf_handler, x0, y0, x1, y1);
   set_line_style_pdf_handler(pdf_handler, NORMAL_LINE_STYLE);
}

#define  PDF_TEXT_ESCAPE	'\\'
#define  PDF_TEXT_START		'('
#define  PDF_TEXT_END		')'

static void string_c_to_pdf(CcpnString dest, CcpnString src)
{
    while (*src)
    {
	if ((*src == PDF_TEXT_ESCAPE) || (*src == PDF_TEXT_START)
						|| (*src == PDF_TEXT_END))
	    *dest++ = PDF_TEXT_ESCAPE;

	*dest++ = *src++;
    }

    *dest = 0;
}

/* fudge factors needed to estimate bounding box of text
   so that can try to correctly place text for a, b non-zero */
#define XFRACTION  0.45
#define YFRACTION  0.9

void draw_text_pdf_handler(Pdf_handler pdf_handler, CcpnString text,
				float x, float y, float a, float b)
{
    Line line;

    x *= pdf_handler->ax;
    y *= pdf_handler->ay;

    if (a != 0)
        x -= a * XFRACTION * pdf_handler->fontsize * strlen(text);

    if (b != 0)
        y -= b * YFRACTION * pdf_handler->fontsize;

    string_c_to_pdf(line, text);

    FILE *fp;
    OPEN_FOR_APPENDING(fp, pdf_handler->fp);

//    fprintf(pdf_handler->fp, "BT\n");
//    fprintf(pdf_handler->fp, "/%s %d Tf\n", pdf_handler->font,
//						pdf_handler->fontsize);
//    fprintf(pdf_handler->fp, FMT FMT "Td\n", x, y);
//    fprintf(pdf_handler->fp, "(%s) Tj\n", line);
//    fprintf(pdf_handler->fp, "ET\n");
    fprintf(fp, "BT\n");
    fprintf(fp, "/%s %d Tf\n", pdf_handler->font,
						pdf_handler->fontsize);
    fprintf(fp, FMT FMT "Td\n", x, y);
    fprintf(fp, "(%s) Tj\n", line);
    fprintf(fp, "ET\n");
    FCLOSE(fp)
}

void draw_dash_box_pdf_handler(Pdf_handler pdf_handler,
			float x0, float y0, float x1, float y1)
{
   int dash_length = 2, gap_length = 6;

   draw_dash_line_pdf_handler(pdf_handler, x0, y0, x0, y1, dash_length, gap_length);
   draw_dash_line_pdf_handler(pdf_handler, x0, y1, x1, y1, dash_length, gap_length);
   draw_dash_line_pdf_handler(pdf_handler, x1, y1, x1, y0, dash_length, gap_length);
   draw_dash_line_pdf_handler(pdf_handler, x1, y0, x0, y0, dash_length, gap_length);
}
void set_color_pdf_handler(Pdf_handler pdf_handler, float *color)
{
    float r = color[0], g = color[1], b = color[2];

    if (!pdf_handler->isColor)
    {
        if (r + g + b < 2.99)
            r = g = b = 0.0;
        else
            r = g = b = 1.0;
    }

    FILE *fp;
    OPEN_FOR_APPENDING(fp, pdf_handler->fp);
//    fprintf(pdf_handler->fp, FMT2 FMT2 FMT2 "RG\n", r, g, b);
//    fprintf(pdf_handler->fp, FMT2 FMT2 FMT2 "rg\n", r, g, b);
    fprintf(fp, FMT2 FMT2 FMT2 "RG\n", r, g, b);
    fprintf(fp, FMT2 FMT2 FMT2 "rg\n", r, g, b);
    FCLOSE(fp)
}

void set_black_pdf_handler(Pdf_handler pdf_handler)
{
    float color[NCOLORS];
 
    ZERO_VECTOR(color, NCOLORS);
 
    set_color_pdf_handler(pdf_handler, color);
}

void set_font_pdf_handler(Pdf_handler pdf_handler, CcpnString name, int size)
{
    if (equal_strings(name, "Times"))
	strcpy(pdf_handler->font, "Times-Roman");
    else
    	strncpy(pdf_handler->font, name, LINE_SIZE-1);
    pdf_handler->font[LINE_SIZE-1] = 0;
    pdf_handler->fontsize = size;
}

void set_line_style_pdf_handler(Pdf_handler pdf_handler, int line_style)
{
/*  TODO
    if (line_style == NORMAL_LINE_STYLE)
	fprintf(pdf_handler->fp, "PS_nondashed\n");
    else *//* line_style == DASHED_LINE_STYLE *//*
	fprintf(pdf_handler->fp, "PS_dashed\n");

    FILE *fp;
    OPEN_FOR_APPENDING(fp, pdf_handler->fp);
    if (line_style == NORMAL_LINE_STYLE)
	fprintf(fp, "PS_nondashed\n");
    else *//* line_style == DASHED_LINE_STYLE *//*
	fprintf(fp, "PS_dashed\n");
	FCLOSE(fp)
*/
}

void set_line_width_pdf_handler(Pdf_handler pdf_handler, float line_width)
{
/* TODO */
}


void get_region_pdf_handler(Pdf_handler pdf_handler, float *x0, float *y0,
                                                        float *x1, float *y1)
{
    *x0 = pdf_handler->x0;
    *y0 = pdf_handler->y0;
    *x1 = pdf_handler->x1;
    *y1 = pdf_handler->y1;
}

void get_text_size_pdf_handler(Pdf_handler pdf_handler, CcpnString text,
                                                        float *w, float *h)
{
/* cannot really do this properly */
    int size = pdf_handler->fontsize;

    *w = size * strlen(text);
    *h = size;
}

void fill_triangle_pdf_handler(Pdf_handler pdf_handler, float x0, float y0,
                                float x1, float y1, float x2, float y2)
{
    x0 *= pdf_handler->ax;
    y0 *= pdf_handler->ay;
    x1 *= pdf_handler->ax;
    y1 *= pdf_handler->ay;
    x2 *= pdf_handler->ax;
    y2 *= pdf_handler->ay;

    FILE *fp;
    OPEN_FOR_APPENDING(fp, pdf_handler->fp);
//    fprintf(pdf_handler->fp, FMT FMT " m\n", x0, y0);
//    fprintf(pdf_handler->fp, FMT FMT " l\n", x1, y1);
//    fprintf(pdf_handler->fp, FMT FMT " l\n", x2, y2);
//    fprintf(pdf_handler->fp, "h\n");
//    fprintf(pdf_handler->fp, "f\n");
    fprintf(fp, FMT FMT " m\n", x0, y0);
    fprintf(fp, FMT FMT " l\n", x1, y1);
    fprintf(fp, FMT FMT " l\n", x2, y2);
    fprintf(fp, "h\n");
    fprintf(fp, "f\n");
    FCLOSE(fp)
}

